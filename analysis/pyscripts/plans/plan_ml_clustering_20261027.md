# Plan: ML in clustering and cluster-level energy regression (HyCal PbWO4, elastic e-p)

**User decisions (2026-10-27)**

- Regression: augment-first (correct on top of E_R3); absolute energy excluded as a feature; shuffled-truth control mandatory — confirmed.
- Order: e-p run types first (0.73/2.24/3.49 H2, 2.24 C); X17/Ta per run afterwards.
- X17 needs three-cluster detection (e+e− pair plus scattered electron). The user recalls an algorithm named roughly 'millepede' for this; note: Millepede-II is a track-alignment least-squares package, so the intended algorithm is still to be identified — record the name as given and resolve before planning the X17 unit. Candidates to check: PRad island clustering with shower-shape splitting, CMS ECAL PF/mustache clustering, PandoraPFA-style topological clustering.
- ML-assisted clustering for X17 stays outside the σ_E/E paper's e-p scope but is retained in this plan as the X17 follow-on.

**Date:** 2026-10-27. **Status:** draft for user approval; nothing implemented.
**Reads first:** `notes/PbWO4_methods_summary_20261027.md` §4–§6, §8; `notes/PbWO4_energy_resolution_summary.md` §3.4, §3.5, §3.6, §4; `prad2det/src/HyCalCluster.cpp` (read only); `position_corr.py`; `percrystal_calib.py`; the R3-bank per-cluster dumps `output/scan_20261024/spl_percrystal/scan_<run>/split_<n>/result_dump.npz` (118) plus `output/percrystal_spl/dumps_extra/024386_00060_dump.npz` (the best-50 dumps under `output/percrystal_best50_20261026/dumps/` carry the other bank and are not used).
Paths are relative to `analysis/pyscripts/`. Standing constraints: Python only, no C++ change without approval; no commits; `database/` never written; Ta runs 25800/26138/26230 per run, never pooled; outputs under `output/ml_cluster_20261027/`; no package installs without approval (§8).

## 1. Current state (from code)

- Clustering (`det.HyCalCluster`, island algorithm; methods note §4): hits sorted by E; seed = largest unconsumed hit ≥ 10 MeV; BFS growth over edge neighbours (corner_conn off), one pulse per module, min_module_energy 1 MeV; local maxima = hits ≥ 10 MeV with no higher neighbour; profile-based split (SimpleProfile in the benchmark, Geant4 profile in production), 6 iterations, drop fraction < 0.01, no split if ≥ 100 hits or ≥ 10 maxima; energy = sum of split-weighted hit energies; position = log-weighted centroid over the seed 3×3 (w = max(0, 3.6 + ln(E_i/E_tot))); leakage and depth corrections exist in C++ but are off / not called; non-linearity is an identity for PbWO4.
- Cluster energy chain (R3): E_u = A·I_int·f(A) per pulse → LMS (null) → cluster sum with per-crystal gains g_i → position correction E/R(u, v) (hybrid per-module poly8 at 3.49 H2, 2.24 H2, C; bilinear at 0.73 and Ta) → seed-indexed elastic-window gain g_M → spectrum fit (CB or ±2σ Gaussian) per split.
- Reference numbers (R3, legacy, 118 splits): σ_E/E 1.492 / 1.854 / 1.797 / 3.718 % at 3.49 H2 / 2.24 H2 / 2.24 C / 0.73 H2; Ta 1.986 ± 0.028 (25800), 1.950 ± 0.114 (26138), 1.983 ± 0.066 (26230). 2-term fit (CB for all, H2): a = 2.49 ± 0.02 %·√GeV, b = 0.68 ± 0.03 %; module medians 1.52 / 1.84 / 1.81 / 3.50 %.
- What is left after R3: the residual position map is flat (after-map RMS 0.00–0.05 % at 2.24, 3.49, C), but the edge-vs-centre width excess remains (1.85 vs 1.55 % at 3.49 GeV, 0.2–0.3 pp; ≈ 1 pp at 0.73 GeV), ring-1 modules are worse by 0.11 pp at 3.49 GeV, one map per run type corrects the average only (ring/side differences 0.1–0.3 %), and the hot spot W461/W427/W462 is unexplained (summary §3.5). These are the targets of a cluster-level regression.

## 2. Target (a): ML-assisted clustering. Does PRad need it?

- Physics case: elastic e-p with one cluster per event; the resolution pool requires exactly one cluster > 500·kE MeV. Møller pairs are back-to-back at the ±20° level and do not overlap in the inner region; the seed-time gate (overlap of out-of-time showers) is a null on e-p and C and gains −0.048 / −0.030 / −0.155 pp on Ta runs 25800 / 26138 / 26230 only (20261027 §7). The split algorithm only acts when two local maxima sit in one island.
- Conclusion to verify, not assume: for the σ_E/E figure of merit, shower splitting/merging and seed selection are not a resolution lever. The one bounded check is V4 (§6): the fraction of pool clusters flagged kSplit or with ≥ 2 local maxima, per run type and per Ta run. If < 1 % on H2/C, target (a) is closed for this paper with one sentence; Ta and the X17 two-body final states are a separate analysis (Open Question 4).
- If not negligible: a learned split-weight (replacing SimpleProfile's fraction for shared hits) is the only candidate; separate plan.
- Revisited for X17 (three clusters per event) after the e-p units; see User decisions.

## 3. Target (b): cluster-level energy regression

- Goal: a per-run-type regressor f(features) → E′(θ_seed)/E_cluster that absorbs what the hybrid map and the seed-indexed gain cannot (edge width, ring/side dependence, shape-dependent leakage), evaluated on the R3 pool with the R3 estimators.
- Two arms: **augment** (base E = E_R3 after map and g_M; the regressor learns the residual; cheap, directly comparable) and **replace** (base E = Σ g_i E_i before the map; the regressor must rediscover R(u, v)). Augment first; replace only as a cross-check that the regressor is not just refitting the map.
- Models: ridge regression on engineered features (numpy only; V0); gradient-boosted trees (`sklearn.ensemble.HistGradientBoostingRegressor`, ≈ 200 trees, depth ≤ 6; V1); a 3-layer MLP (64–64–1, torch CPU) on the same features plus the 5×5 energy-fraction map (V1b). Monotone-in-nothing; the target is a ratio near 1, Huber loss.
- Features (per cluster, all from the dump schema: `E_cluster, E_seed, x_cl, y_cl, x_seed, y_seed, size_x, size_y, row, col, theta_deg, f_lms, ptr→midx/E`): u, v = (centroid − seed centre)/size; seed module id (categorical, or ring index and side); E_seed/E_cluster; 3×3 and 5×5 fractions; n_members; second moments of the member energies about the centroid; max neighbour fraction and its direction; run type. **Not** the absolute E_cluster, not E_cluster/E′ (§4). Ta: run index as a feature (report-only, §8).
- Truth: `eprime_MeV` = E′(θ) at the seed-module centre minus energy loss, as the per-crystal fit uses (methods §3); E′ varies by 3·10⁻⁴–10⁻³ across a crystal (largest at 3.49 GeV in the outer ring; estimate, unsourced), so V0 also tests E′ at the cluster centroid as the target. E_beam itself is not used.

## 4. Truth-leakage protocol (the elastic selection uses the energy)

- Training pool: the R3 pool (one cluster > 500·kE, seed r < 200 mm, seed E ≥ 200·kE, rectangle out) restricted to r = E_R3/E′ ∈ [0.85, 1.12] (the window `position_corr.py` already uses), so the radiative tail and Møller leak-through do not set the fit. The evaluation pool is the untouched R3 pool with the R3 estimator; no window.
- Sculpting controls, all mandatory: (i) no absolute-energy feature (§3); (ii) shuffled-truth control: refit with E′ permuted within seed module → Δσ/E must be 0 within SE; (iii) ε_peak, ε_yield and the E/E′ < 0.90 tail fraction must not move by more than 2 SE from R3; (iv) training on even split files, evaluation on odd, and a leave-one-run-out pass at 2.24 GeV H2 (5 runs) and 3.49 GeV H2 (3 runs); (v) the per-module medians must improve at the edge bins (|u|, |v| > 0.4) more than at the centre, or the gain is suspect.
- The seed-indexed gain g_M is refitted after the regression in the augment arm only if the regressor shifts module peaks by > 0.2 % (it should not, since module id is a feature).

## 5. Metrics

- Per-split σ_E/E on all 118 splits (119 with 24386.00060 where the dump exists), paired Δ vs R3 ± SE, legacy path (gated differs by ≤ 0.01 pp on H2/C), preferred estimator and CB-for-all; Ta per run.
- a/b decomposition (CB for all, H2 only) with the fit's χ²; per-module medians (N > 1000) and the edge/centre widths; the five controls of §4.
- Adoption rule as before: no H2/C run type worse beyond its SE, and a gain of ≥ 0.02 pp at > 2 SE in at least two run types.

## 6. Steps (bounded units, each with an acceptance criterion)

- V0, cheapest unit that proves or kills the idea (numpy only, < 1 h): `cluster_regression.py` reads the dumps, gain files and R3 maps named in §9, rebuilds E_R3 per cluster, and (gate) reproduces the R3 per-split σ_E/E within 0.005 pp on 4 splits per run type; then ridge regression on the §3 features, even→odd, augment arm, 3.49 and 2.24 GeV H2 first, with the shuffled-truth control. Accept: Δ with SE and the control at 0. Kill: if ridge and the oracle bound below both give |Δ| < 0.02 pp at 3.49 GeV.
- V0b oracle bound (same script): in-sample per-(seed module, u-bin, v-bin, E_seed/E_cl-bin) lookup correction; its in-sample Δ is the ceiling any regressor can reach on these features (overfit by construction). Accept: a ceiling per run type.
- V1: HistGradientBoosting (needs sklearn, §8) and V1b MLP (needs torch) on the same split; model cards. Accept: validation Δ within 2 SE of training Δ; §4 controls pass.
- V2: full 118-split evaluation of the best model, leave-one-run-out, a/b fit, module medians, edge/centre table; Ta per run. Accept: tables in `output/ml_cluster_20261027/compare/`.
- V3: replace arm as a cross-check (one run type). Accept: within 0.01 pp of augment or explained.
- V4: target (a) check: one split per run type (Ta per run) through the benchmark, reading `cluster.flag` (kSplit bit) and `cluster.hits` (local maxima recounted in Python) for the pool clusters; add the two columns to the cluster TSV if absent (Python only). Accept: one table of fractions and a decision sentence.
- V5: paper material: one table (R3 vs regression per run type, a/b), one figure (per-module σ_M/μ_M map before/after, or edge/centre widths), one paragraph on the protocol and the controls.

Compute (18 cores, 48 GB, no GPU; estimate, unsourced): dumps hold ≈ 10⁵ clusters per split, ≈ 10⁷ in total; ridge seconds; HGB on 10⁶ rows × 20 features ≈ 1–5 min per run type; MLP ≈ 10–30 min per run type; the 118-split evaluation is dump-based (no evio rescan).

## 7. Improvement estimate from other experiments

Survey (researcher-worker, 2026-10-27, web). Few papers isolate the gain of a regression over a per-crystal calibration plus analytic corrections, none on PbWO4 at GeV fixed-target energies, and none reports the change of the constant term b. What was found:

| Source | Detector, range | Baseline → regression | Reported gain | Truth |
|---|---|---|---|---|
| ATLAS Run 1 e/γ calibration, EPJC 74 (2014) 3071, arXiv:1407.5063 §5.3 | LAr EM, E_T 5 GeV–TeV | initial calibration → MVA (BDT) on E_calo, presampler fraction, depth, η, barycentre, R_conv, binned in η, E_T | resolution better by 3–10 % (barrel) and 10–15 % (endcap) for unconverted γ, ≈ 20 % converted γ, "a few percent" for electrons (10–30 % in the crack); non-linearity better by > 2×; H→γγ mass 10 % (max 15 %); the estimator reproduces the expected sampling term for unconverted γ, so the gain sits in tails and material effects (inference) | simulation E_true |
| ATLAS Run 2, JINST 14 (2019) P03017, arXiv:1812.03848 §4 | same | BDT without → with E4 gap-scintillator inputs | ≈ 20 % (e), 5 % (unconverted γ) at 1.4 < \|η\| < 1.6, 50–100 GeV; no overall MVA-vs-old number | simulation |
| CMS ECAL, JINST 8 (2013) P09009, arXiv:1306.2016 §4.4; photons arXiv:1502.02702 §4.3; Run 2 e/γ arXiv:2012.06888 §6.1 | PbWO4, 5 GeV–TeV | factorised parametric corrections → semi-parametric BDT (seed position in crystal, E3×3/E5×5, R9, widths, ρ, N_vtx) | "better energy resolution" and "significantly improved"; no number isolates the regression; H→γγ 8 TeV σ_eff −20 to −30 % with the intercalibration included (arXiv:1407.0558) | simulation E_true |
| Simkina, Couderc, Malclès, Sahin, EPJC 84 (2024) 639, arXiv:2311.17914, doi:10.1140/epjc/s10052-024-12978-1 | toy PbWO4 51×51 crystals (CMS-barrel-like: a 3 %, noise 0.167 GeV, c 0.35 %), 1–100 GeV | raw PFCluster → BDT (E_max, E_2nd, asymmetries, covariances, E3×3, position): −26 to −43 % by energy (App. A); BDT → CNN/GNN on the crystal grid: ≈ −10 % single γ, ≈ −60 % for two overlapping γ (§6.1) | simulation E_gen |
| CMS ECAL GNN calibration, Rothman, NeurIPS ML4PS 2021 (neurips.cc/virtual/2021/37109) | PbWO4 | legacy CMS → dynamic GNN | "up to 10 %" for calibration electrons; snippet-level | not stated |
| Belayneh et al., EPJC 80 (2020) 688, arXiv:1912.06794 | cell-level, resampled ATLAS/CMS-like | XGBoost → DNN/CNN/GN | surpass on fine granularity, "match but do not improve" on coarse CMS-like granularity; no number | simulation |

No source states that the gain is mainly in the constant term; the closest is Bendavid, CHEP 2016 (indico.cern.ch/event/472938/contributions/1150753, slides): position- and angle-dependent containment losses do not factorise, which is what a regression recovers. Belle II ECL uses a position-binned leakage lookup (no MVA number found); LHCb, PANDA, PRad: nothing quantitative.

Translation to our setting (R3 CB-for-all fit: a = 2.49 %·√GeV, b = 0.68 %, E = 0.727 / 2.236 / 3.484 GeV; σ = a/√E ⊕ b). b² is 5 / 14 / 21 % of σ² at the three energies, so a regression that acts only on b has limited headroom:

| b (%) | σ_E/E at 0.73 / 2.24 / 3.49 GeV (%) | Δ vs R3 (pp) | relative |
|---|---|---|---|
| 0.68 (R3) | 2.998 / 1.799 / 1.497 | — | — |
| 0.55 | 2.972 / 1.754 / 1.443 | −0.03 / −0.05 / −0.05 | −1 / −2.5 / −3.6 % |
| 0.45 | 2.955 / 1.725 / 1.408 | −0.04 / −0.07 / −0.09 | −1.5 / −4 / −6 % |
| 0.35 (the CMS-like toy's c) | 2.941 / 1.702 / 1.379 | −0.06 / −0.10 / −0.12 | −2 / −5 / −8 % |
| 0 (floor) | 2.920 / 1.665 / 1.334 | −0.08 / −0.13 / −0.16 | −3 / −7 / −11 % |

- Expected: our baseline already holds per-crystal gains, a per-module position map and a seed-indexed gain, i.e. the "analytic corrections" stage of ATLAS/CMS. The analogues are ATLAS's "a few percent" for electrons and Simkina's ≈ 10 % for a network over a BDT (on a full 51×51 grid input with overlapping showers; our events are single, isolated showers, so the lower end applies). A realistic target is b 0.68 → 0.45–0.55 %, i.e. σ_E/E 1.49 → 1.41–1.44 % at 3.49 GeV, 1.81 → 1.73–1.75 % at 2.24 GeV, 2.97 → 2.96–2.95 % at 0.73 GeV (estimate, unsourced; the own-data evidence is the edge/centre width excess 0.2–0.3 pp at 3.49 GeV, the ring-1 excess 0.11 pp and the unexplained hot spot, summary §3.5, which are exactly the features the regression sees).
- Honest floor: b carries pieces no cluster feature sees: the in-sample-vs-held-out gain difference 0.01–0.02 pp, module gain stability (run-to-run SD 0.16 %), the ±0.4 % sub-sample-phase MF bias (summary §3.1), and the preferred-estimator dependence (±2σ Gaussian vs CB: b 0.68 vs −0.98, unphysical). The sampling term a is untouched by any cluster regression (ATLAS reproduces its sampling term; inference). The 0.73 GeV point gains ≤ 0.06 pp whatever happens to b. Any Δ larger than the b = 0 row at 3.49 GeV (−0.16 pp) is sculpting, not physics, and fails §4.

## 8. Risks and compatibility

- Installs: sklearn and torch are not present (numpy 2.4.6, scipy 1.18.0 are). V0/V0b need nothing. V1 needs `scikit-learn`; V1b needs CPU torch; both go into a separate environment after approval (`python-env` skill).
- In-sample gains and maps: R3's own gains and hybrid maps are in-sample (methods §9); the regression's held-out Δ is therefore compared against an in-sample reference, which favours R3 by 0.01–0.02 pp (summary §5.1 held-out gains). Report both.
- Dump coverage: the dumps were produced for the per-crystal fit; if the pool in the dump differs from the scan's pool, the V0 gate fails and the dump stage must be rerun (Python, ≈ 45 min at 12 parallel; no evio needed beyond the local splits).
- Ta: ≈ 16 calibrated modules per split; a per-run regressor will not train. Ta is report-only with the 2.24 GeV H2 model applied, per run.
- No C++ change: the regression lives in Python after the benchmark, like the position map; if adopted, a `--cluster-corr FILE` hook in the benchmark (Python) is the integration point.

## 9. Acceptance criteria and validation commands

- `python3 cluster_regression.py rebuild --dumps 'output/scan_20261024/spl_percrystal/scan_<run>/split_*/result_dump.npz' output/percrystal_spl/dumps_extra/024386_00060_dump.npz --gains output/percrystal_spl/<run>/gains_<run>_all.tsv --poscorr output/poscorr_20261025/hybrid/<rt>/poscorr_hybrid_<rt>_all.json (3.49 H2, 2.24 H2, C) | output/poscorr_20261025/<rt>/poscorr_<rt>_all.json (0.73 H2, Ta) --out output/ml_cluster_20261027/v0/` reproduces R3 σ/E on 4 splits per type within 0.005 pp (gate).
- `python3 cluster_regression.py fit --model ridge --arm augment --train even --eval odd --control shuffle ...` writes `v0/ridge_<slug>.md` with Δ ± SE and the control.
- `python3 compare_scan_dirs.py`-style paired tables from the per-split outputs; `plot_module_resolution.py --group-runs` on the corrected per-cluster energies.

## Open Questions

1. May `scikit-learn` (and a CPU `torch`) be installed in a separate environment for V1/V1b?
2. Augment-first (recommended) or replace-first?
3. Is the no-absolute-energy feature rule acceptable, or should E_cluster enter with the sculpting controls as the only guard? (Recommendation: exclude.)
4. Target (a) for the X17/Ta two-body final states is outside this paper's σ_E/E scope; confirm it stays out.

## Todo

- [ ] V0 rebuild gate (R3 reproduced from dumps) and ridge regression with shuffled-truth control, 3.49 and 2.24 GeV H2
- [ ] V0b oracle ceiling per run type
- [ ] V1 HGB and V1b MLP (after Open Question 1)
- [ ] V2 118-split evaluation, leave-one-run-out, a/b, module medians, edge/centre, Ta per run
- [ ] V3 replace-arm cross-check
- [ ] V4 clustering-need check (split/local-maxima fractions)
- [ ] V5 paper table, figure, protocol paragraph
