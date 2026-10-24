# PbWO4 energy resolution: per-crystal calibration, seed-time neighbours, free-τ pulse fit, pile-up

**Date:** 2026-10-23
**Branch:** `waveform_jz` (HEAD `5359d63`, pushed; everything after it is uncommitted, §1)
**Plan:** `analysis/pyscripts/plans/plan_20261023.md` (Phases 0, 1, 3; Round 2; "Study: free-tau pulse fit"; decisions D1–D5)
**Previous notes:** `analysis/pyscripts/notes/PbWO4_energy_resolution_20261021.md` (LMS, elastic-window refinement, Legacy vs New), `analysis/pyscripts/notes/PbWO4_energy_resolution_20261022.md` (mf template method, gate scans, working point)

Paths are relative to `analysis/pyscripts/` in `/Users/L00338853/work/PRad/prad2evviewer` unless they start with `prad2det/`, `prad2dec/`, `database/` or `docs/` (repo root), or are absolute. "mf", "rawsum", Legacy and gated (New) are as in the 20261022 note. σ/E is post-refinement σ_E/E, inner ring r < 200 mm, preferred fit (D5: CB for 3.49 GeV H2, C, Ta; Gaussian for 2.24 and 0.73 GeV H2), legacy path unless noted. Δ in pp; "paired" means mean ± SE of per-split differences. "mf baseline" is the mf working point `mf_traw_w16_p100-200` (20261022 note §7).

## 1. Status

| Item | State |
|---|---|
| Benchmark | `benchmark_hycal_timing_test_20261021.py`, unchanged role. |
| New benchmark options (in `5359d63`) | `--cal-method percrystal` (member-level χ² gain fit on the job's own events, writes `<out>_gains.tsv`); `--gain-file PATH` (applies a pooled per-crystal table after clustering, before the refinement); `--cluster-dump {PATH,auto}` (calibration-pool NPZ for `percrystal_calib.py fit`); `--mf-nb-amp {off,seedtime}`, `--mf-nb-window {3,5}`, `--mf-nb-time {fixed,search}` (§5). |
| New benchmark options (uncommitted) | `--mf-shape {fixed,bounded,prior}`, `--mf-shape-nsigma` (3.0), `--mf-shape-weight {uniform,seedmodel}` (§6). Defaults byte-identical. `--exclude-rowcol` now read 1-based inclusive (`15 20 15 20` = centred 6 × 6; §8a.1). |
| Driver (uncommitted) | `benchmark_hycal_timing_scan.py --gain-file-map {off,auto}`, default off. `auto` passes `--gain-file output/percrystal/<run>/gains_<run>_all.tsv`; a missing table skips the split. Output stem `benchmark_hycal_timing_scan_mf_pc_*`. rawsum and mf commands unchanged against HEAD. `--gain-dir DIR` (default `output/percrystal`) sets the table directory for `auto`; `--write-cluster-tsv` passed through to the benchmark. |
| Harness (in `5359d63`) | `gate_scan_20261022.py`: `run --runs R1,R2 --splits {even,odd,all} --tags T1,T2 --configs-file F`; `{run}`, `{runtype}`, `{job_dir}` placeholders in `bench_args`; `restable --configs-file F --baseline-tag T` writes `variant_table_<tag>[_x17].{md,tsv}` with `eps_yield`, `eps_peak` (floor 0.99) and `n_cal`. |
| New scripts | `percrystal_calib.py` (committed; uncommitted change: `fit` provenance records `fixed_windows`, `max_iter`, `argv`); `mf_fit_pulse.py`, `mf_fit_study.py`, `plot_pileup_examples.py`, `plot_module_resolution.py` (untracked). |
| Git | Modified since `5359d63`: the benchmark, the driver, `percrystal_calib.py`, `../plans/plan_20261023.md`. Untracked: the four scripts above and this note. `output/` is gitignored. |

## 2. Resolution decomposition (Phase 0)

Source: `output/resolution_fit_20261023/phase0.py` and `out.txt` (rerun identical). Inputs: per-split `result.tsv` and `result_window.tsv` of `output/scan_20261022_mf/` and `output/scan_20261021_full/`; E = 3.484, 2.236, 0.727 GeV; C and Ta not in the fits.

σ² = a²/E + b² (2-term, 1 dof) and σ² = a²/E + b² + c²/E² (3-term, exactly determined). a in %·√GeV, b in %, c in %·GeV. Errors are SEM-only and the 2-term χ² is 180–820 for 1 dof, so the real uncertainty on the mf legacy preferred fit is about ±0.4 on a and ±0.5 on b.

| Scan / path / estimator | σ/E % at 0.73 / 2.24 / 3.49 | 2-term a, b (χ²) | 3-term a, b, c |
|---|---|---|---|
| mf legacy, preferred | 4.350 / 2.223 / 1.873 | 3.13, 0.77 (178) | 2.19, 1.29, 2.37 |
| mf legacy, CB for all | 3.039 / 1.907 / 1.873 | 1.85, 1.49 (448) | a² < 0; b 1.97, c 2.19 |
| mf legacy, Gaussian for all | 4.350 / 2.223 / 1.943 | 2.83, 1.18 | – |
| mf gated, preferred | 4.352 / 2.220 / 1.873 | 3.12, 0.79 (184) | 2.16, 1.30, 2.39 |
| mf gated, CB for all | 3.042 / 1.905 / 1.873 | 1.89, 1.47 (498) | a² < 0; b 1.98, c 2.21 |
| rawsum legacy, preferred | 4.655 / 2.263 / 1.936 | 3.04, 1.02 (403) | 1.52, 1.54, 2.92 |
| rawsum legacy, CB for all | 3.391 / 2.028 / 1.936 | 2.34, 1.40 (789) | a² < 0; b 1.94, c 2.28 |
| rawsum gated, preferred | 4.680 / 2.293 / 1.971 | 3.02, 1.10 (366) | 1.50, 1.59, 2.93 |
| rawsum gated, CB for all | 3.406 / 2.040 / 1.971 | 2.20, 1.52 (821) | a² < 0; b 2.03, c 2.40 |
| mf legacy, **pre**-refinement, preferred | 4.423 / 3.234 / 3.283 | 3.12, 2.48 (2.7) | – |

| Reading | Detail |
|---|---|
| b depends on the estimator | 0.77 (preferred) to 1.49 % (CB for all). The preferred set mixes a Gaussian at 0.73 and 2.24 GeV (radiative tail) with CB at 3.49 GeV. With CB throughout the 2.24 → 3.49 GeV step is 0.03 pp, which forces b ≥ 1.5 %. CB pair (2.24, 3.49 only): a 0.90, b 1.81. |
| Working figure | b ≈ 1.3 % (3-term; also √(1.87² − 1.39²) = 1.25 % after the 2.6 %/√E beam-test stochastic term). Constant share of σ² at 3.49 GeV: 47 % (17 % if b = 0.77, 63 % if b = 1.49). |
| Rough model of the 20261022 note §8 | a ≈ 2.5, b ≈ 1.2 is not the 2-term result (3.13, 0.77); it is close to the 3-term one, which needs c ≈ 2.4 %·GeV (24 MeV at 0.73 GeV). |
| 3-term shares (mf legacy preferred) | a/√E 2.57 / 1.46 / 1.17, b 1.29, c/E 3.27 / 1.06 / 0.68 % at 0.73 / 2.24 / 3.49 GeV. |
| Refinement | Pre → post b 2.48 → 0.77–1.3 %: the seed-indexed step already removes most of a 2.5 % constant term. |

Per-seed-module gains g_M (`result_window.tsv`, `used = 1`, mf full scan):

| Run type | Modules / split | RMS of g_M | Statistical part | True spread | MAD-based |
|---|---|---|---|---|---|
| 2.24 GeV H2 (25228–25320) | 126–130 | 0.0226 | 0.0028 | **0.0224** | 0.018 |
| 3.49 GeV H2 (24386, 24503, 24512) | 142–150 (24386: 42–123) | 0.0195 (24386: 0.0126) | 0.0034 | **0.0192** | 0.016 |
| 0.73 GeV H2 (24590, 24603) | 96–113 | 0.0121 | 0.0051 | **0.0110** | 0.010 |
| 2.24 GeV C (25308) | 174 | 0.0220 | – | – | 0.018 |

- Run-to-run SD of the per-module run means: 0.0018 (2.24 GeV, 5 runs), 0.0077 (3.49 GeV, driven by 24386), 0.0024 (0.73 GeV); split-to-split SD 0.0032 per module at 2.24 GeV. Legacy and gated g_M agree to 0.0002–0.0005 RMS in mf.
- **Per-crystal lever bound.** The seed carries 0.644 ± 0.124 of the cluster energy (25320 legacy, 2100–2400 MeV). With a 2.2 % module spread, the event-by-event neighbour share leaves 0.25–0.35 % in quadrature after the seed-indexed step: b 1.29 → 1.25 %, i.e. −0.02 to −0.04 pp at 2.24 and 3.49 GeV, if neighbour modules are no worse calibrated than seeds. §3 measures more than this.
- Seed saturation (mf headers): n_sat/n_peaks 2 × 10⁻⁶ at 3.49 GeV; not a lever.
- The exclusion rectangle `15 20 15 20` sits one ring off-centre in the benchmark's 0-based convention (centred: `14 19 14 19`); D1 keeps it for comparability.

## 3. Per-crystal calibration (Phase 1)

### 3.1 Method

| Item | Behaviour |
|---|---|
| Fit | χ²(g) = Σ_k w_k (Σ_i g_i E_ki − E′_k)² + Σ_i λ_i (g_i − g0)², over pool events k inside their seed module's elastic window and member modules i. Sparse normal equations, solved directly. E′ from the seed-centre θ as in the seed-indexed refinement (20261021 note §2.1). |
| Prior | λ_i = the diagonal a module would have from 30 seed events at the mean seed fraction. Prior centre g0 = E′_mean/μ_global (the calibration files carry +8–10 %), not 1. Modules with prior fraction > 0.5 sit near g0; never-seen modules, glass and Veto fixed. |
| Windows | Per seed module, the pre-fit/window/Møller machinery of the seed-indexed step, on the current E_cal. Scan algorithm: windows refit in passes 1–3, then frozen. `fit --fixed-windows` (opt-in): events selected once with the g = 1 windows in pre-refinement energy. |
| Convergence | Every free module \|Δg\| < max(1 × 10⁻³, 0.5 σ_g); at most 6 passes. |
| Benchmark use | `--cal-method percrystal`: in-split fit, E_cal = Σ g_i E_ki per event after clustering (D3). `--gain-file`: pooled table applied the same way, then the default seed-indexed refinement on top. |
| Identity | With g = 1 the member sum equals the cluster energy: 0 violations, max relative difference 5 × 10⁻⁷. |
| Self-test (P1-2) | PASS: 0.25 % RMS against the λ-weighted expectation (0.35 % raw, including prior shrinkage). |

### 3.2 Results (S_dev = 25320, 24503, 24590, 25308; paired vs the mf baseline)

| Run type | In-split `--cal-method percrystal` | **Held-out (P1-4)** | Residual g_M median / RMS (P1-5) | RMS(g_i − g_M) (P1-7) |
|---|---|---|---|---|
| 3.49 GeV H2 | −0.035 ± 0.002 | **−0.066 ± 0.003** | 1.0000 / 0.37 % | 0.61 % |
| 2.24 GeV H2 | −0.062 ± 0.003 | **−0.095 ± 0.003** | 0.9999 / 0.42 % | 0.70 % |
| 0.73 GeV H2 | −0.106 ± 0.036 | **−0.066 ± 0.022** | 1.0014 / 0.71 % | 0.49 % |
| 2.24 GeV C | −0.064 ± 0.003 | **−0.091 ± 0.003** | 0.9998 / 0.40 % | 0.70 % |

- Held-out: pooled gains from the even splits of a run applied to its odd splits and vice versa (`--gain-file`, seed-indexed step on top). Gated agrees within 0.003. The success rule (Δ ≤ −0.02 pp at 2.24 and 3.49 GeV, |Δ| > 2 SE) is met. At 0.73 GeV with CB for all: −0.004 ± 0.014.
- The held-out gain is 2–3× the Phase 0 bound (§2).
- P1-5 target RMS ≤ 0.5 % is met except at 0.73 GeV (0.71 %). P1-7: median g_i − g_M +0.05 to +0.18 %, correlation 0.92–0.97.
- P1-3: 67 of 72 in-split path-fits converged in 4–6 passes; 5 stopped at the 6-pass cap.
- In-split S_full (P1-8, 118 splits): 3.49 GeV −0.017 ± 0.005, 2.24 GeV H2 −0.059 ± 0.002, 0.73 GeV −0.155 ± 0.035, C −0.064 ± 0.003, Ta +0.049 ± 0.031 (reported only). At 0.73 GeV the calibrated block loses 3–4 % of the events (n_cal 0.961), so part of that Δ may be selection. Tables: `output/scan_20261023/variant_table_mf_percrystal[_x17].{md,tsv}`.
- **W428 and W494** differ from g_M by > 2 % in most 2.24 GeV splits in-split (P1-7). With pooled gains (§4) they do not stand out (≤ 1.3 RMS; W428 −0.3 to −0.5 % at 2.24 GeV).

### 3.3 Stability and initialisation (P1-6)

| Check | Per-crystal | Seed-indexed |
|---|---|---|
| Run-to-run SD of per-module run means, 2.24 GeV H2 (5 runs, 99 modules), mean / median | 0.0019 / 0.0013 | 0.0017 / 0.0014 |
| 3.49 GeV, 24503 vs 24512 (125 modules) | 0.0019 / 0.0017 | 0.0018 / 0.0017 |
| 0.73 GeV, 24590 vs 24603 (72 modules) | 0.0019 / 0.0016 | 0.0018 / 0.0015 |
| Statistical expectation on a run mean (2.24 / 3.49 / 0.73) | 0.0009 / 0.0009 / 0.0020 | |
| Split-to-split SD | 0.0027 / 0.0026 / 0.0054 | |

A real run-to-run component of ≈ 0.0016 at 2.24 and 3.49 GeV, the same for both methods.

Initialisation independence (025320.00050, g = 1 vs seed-indexed start), rule max |Δg| < 1 × 10⁻³:

| Algorithm | max \|Δg\| | Note |
|---|---|---|
| Scan algorithm (windows refit in passes 1–3) | 4.0 × 10⁻³ (130 free modules; 4 above 10⁻³; max 2.35 σ_g; median 3.7 × 10⁻⁵) | Fails. Cause: the window refits depend on the start. tol_sigma = 0 gives the same. |
| Windows refit every pass | 1.8 × 10⁻³ | Never converges. |
| `fit --fixed-windows` | 0, 2 passes | Passes. Against the scan algorithm: median \|Δg\| 8.5 × 10⁻⁴ (single split), 6.7 × 10⁻⁴ (pooled 025320 even); RMS 2.6 / 1.9 × 10⁻³; σ/E +0.004 pp (Gauss), +0.007 pp (CB) on 025320.00050. |

`--fixed-windows` stays off by default so the delivered S_dev tables reproduce. The even/odd tables of P1-4 (`gains_<run>_{even,odd}.tsv`, 24503, 24590, 25308, 25320 only) predate it; 24503 odd stopped at the 6-pass cap. The `_all` tables of §4 use it.

## 4. Full scan with pooled per-crystal gains (Round 2)

| Item | Detail |
|---|---|
| Gain tables | `output/percrystal/<run>/gains_<run>_all.tsv` (+ `.json`: `fixed_windows` true, `max_iter`, `argv`), all 14 runs, `percrystal_calib.py fit --splits all --fixed-windows` on the `mf_percrystal` dumps (`output/scan_20261023/mf_percrystal/`; 24386.00060 from `output/percrystal/dumps_extra/`). Every fit converged in 2 passes. 26138 refit reproduces byte for byte. |
| Free modules (prior fraction < 0.5) | 24386 265, 24503 283, 24512 283, 24590 266, 24603 264, 25228–25320 284, 25308 284; X17: 25800 144, 26138 54 (11k pool events), 26230 110. |
| Scan | Driver, `--pulse-method mf --window 16 --pre-window 100 200 --gain-file-map auto`: 119/119 splits, 0 failures, 58 min. `output/scan_20261023_pc/` (`benchmark_hycal_timing_scan_mf_pc_{summary.tsv,findings.md,sigma_E_vs_file_pre.png,sigma_E_vs_file_post.png,all_calibrated_plots.pdf}`). |

σ/E (%), rawsum full scan (`output/scan_20261021_full/`) → mf full scan (`output/scan_20261022_mf/`) → mf + per-crystal. Δ = mf + per-crystal − mf, paired over splits.

| Run type (splits), ★ fit | rawsum | mf | mf + pc **pre** | mf + pc ★ post | Δ legacy | Δ gated |
|---|---|---|---|---|---|---|
| 3.49 GeV H2 (28), CB | 1.936 | 1.873 | 1.846 | **1.826** | −0.048 ± 0.004 | −0.048 ± 0.004 |
| 2.24 GeV H2 (42), Gauss | 2.263 | 2.223 | 2.146 | **2.132** | −0.091 ± 0.001 | −0.090 ± 0.001 |
| 0.73 GeV H2 (16), Gauss | 4.655 | 4.350 | 4.227 | **4.285** | −0.065 ± 0.020 | −0.066 ± 0.020 |
| 2.24 GeV C (8), CB | 2.166 | 2.093 | 2.019 | **2.003** | −0.090 ± 0.004 | −0.090 ± 0.004 |
| Ta X17 (25), CB, reported only | 2.428 | 2.359 | 2.459 | 2.395 | +0.036 ± 0.039 | +0.014 ± 0.031 |

- mf full-scan pre-refinement for comparison: 3.283 / 3.234 / 4.423 / 3.121 / 2.796 (same order). With the gain file the pre spectrum is already close to the post one (25320 CB pre 1.856, post 1.838; mf full CB pre 2.893). At 0.73 GeV the seed-indexed step on top raises σ/E (4.227 → 4.285).
- Per run at 3.49 GeV: 24386 −0.027 ± 0.005, 24503 −0.060, 24512 −0.062. 3.49 GeV has 28 splits here (24386.00060 included) against 27 in the harness tables of §3.
- Residual seed-indexed g_M on top: median 1.000, RMS 0.32–0.39 % (0.60–0.65 % at 0.73 GeV).
- **In-sample caveat.** The `_all` gains are fitted on the same splits they are applied to. The S_dev held-out test (§3.2) gave the same size at 2.24 and 3.49 GeV and C.
- Ta: ≤ 34 modules per run, 54–144 free modules; reported, not interpreted.

## 5. Seed-time neighbour amplitudes (Phase 3)

Definition (D2): in mf, per path, for a single-cluster event, the 3 × 3 PbWO4 neighbours of the seed (seed and dead modules excluded) are re-measured with the template at the seed time plus the per-module offset (`--time-calib auto`), linear amplitude, no threshold; the event is re-clustered. Control `mf_pc_tcal`: the same flags without seedtime.

| Check | Result |
|---|---|
| P3-1 self-test | `amp_at` equals the K = 1 closed form; accuracy 2 × 10⁻⁷ at Δt = 0. Amplitude loss for a wrong onset: 0.2 / 0.75 / 1.6 / 2.7 / 5.2 % at 0.5 / 1 / 1.5 / 2 / 3 ns (planning estimate 0.2 / 1.0 / 2.6 / 5.0 / 8). |
| P3-2 linearity (20k events) | Median A_constrained/A_free 0.9992–0.9993 above 50 MeV, 0.987–0.989 at 3–10 MeV. Seed–neighbour Δt σ_MAD 0.85 ns (25320), 0.83 ns (24590); `search` fallback not needed. Constrained neighbour energy ≈ 0.4 % below the free one (median onset offset +0.13 to +0.23 ns). |
| P3-3 empty channels | Read-out neighbours without a detected pulse: 0.1–0.2 % of the replaced channels, mean +4 to +10 MeV. |
| **P3-4 σ/E (S_dev, vs control)** | 3.49 GeV +0.001 ± 0.002; 2.24 GeV H2 +0.015 ± 0.002; 0.73 GeV −0.029 ± 0.040 (gated −0.065 ± 0.046); C +0.016 ± 0.002. **Fails** (2.24 GeV and C worse by ~7 SE; 0.73 GeV not significant). `--time-calib auto` alone changes nothing (0.000; 0.73 GeV +0.014 ± 0.029). |
| S_full (118 splits), Δ seedtime on top of in-split per-crystal | 3.49 GeV +0.002 ± 0.002; 2.24 GeV H2 +0.016 ± 0.001; 0.73 GeV −0.015 ± 0.034; C +0.016 ± 0.002; Ta −0.006 ± 0.018. `output/scan_20261023/variant_table_mf_pc_st[_x17].{md,tsv}`. |

**Zero suppression.** The FADC readout is sparsified. PbWO4 channels read out per event: 7.1 (0.73 GeV), 15–16 (2.24 GeV H2 and C), 19.9 (3.49 GeV); 99.9 % of them are mf candidates. Of the 3 × 3 neighbours of a single-cluster seed, 48 % (24590.00010), 13 % (25320.00050), 12 % (25308.00050) and 4 % (24503.00010) have no waveform. The option can only re-measure channels that were read out, and nearly every read-out channel already yields an mf pulse, so the premise (energy lost in sub-threshold corner neighbours) does not hold in this data.

**Conclusion.** Seed-time amplitudes stay off (`--mf-nb-amp off`). The 5 × 5 and `search` variants are implemented, not run and not validated.

## 6. Free-τ pulse-fit study

Setup: mf finds the pulse unchanged; single-peak unclipped mf pulses (≈ 91 %) are refit for A, t0, τ_r, τ_f by batched Levenberg–Marquardt (`mf_fit_pulse.py`), seeded at the mf result and the h50 τ. Arms: `bounded` (τ in template ± 3σ, σ = 1.4826·MAD of the h50 table; the user's proposal) and `prior` (Gaussian prior of the same σ). Hit times stay mf. E = A_fit · I_int(fitted shape) · κ, then LMS, pooled per-crystal gains and the seed-indexed step, as in §4. Data: 025320.00050/00150/00250 and 024503.00010/00040/00070. `output/mf_fit_study_20261023/`: `uniform/` = headline arms (one σ per pulse, as the h50 extraction); `res/fixed` reproduces `output/scan_20261023_pc` byte for byte on these 6 splits; the root and `dump/` hold the seed-model-weight variant.

| Item | Result |
|---|---|
| Weighting | Per-sample weights frozen at the seed model trade the peak for pre-onset samples: E/E_raw at A > 1000 ADC +9 to +12 % (bounded), −8 to −11 % (prior); σ/E +1.20 / +1.77 pp (25320) and +1.27 / +1.63 pp (24503). Unusable. With uniform weights E_arm/E_fixed = 1.000 ± 0.007 (median) for A ≳ 50 ADC; below 50 ADC bounded gives 0.966–0.993. |
| Shape variation | Real and far above the fit noise: at A > 300 the within-module spread is 0.6–0.9 ns (τ_r) and 0.4–0.8 ns (τ_f), against per-fit σ 0.02–0.05 ns. |
| χ²/dof (median) | Freeing τ lowers it by only 15–20 %; the two-τ form does not describe large pulses (table below). |
| Against the WaveAnalyzer integral | rms(E/E_raw), 25320: 0.024 (fixed) → 0.012–0.013 at 100–300 ADC; 0.010–0.012 → 0.007–0.010 at 1000–4000 ADC (both runs similar). |
| Bound hits (bounded) | Any bound: 45 % of fits (25320), 49 % (24503). 3.2–5.5 % at A ≥ 500 ADC, 5.3–6.8 % at 200–500; 24–84 % (25320) and 30–86 % (24503) below 100 ADC. |
| Failures (fallback to the fixed energy) | Bounded 0.13 % / 0.06 %; prior 0.16 % / 0.32 % (25320 / 24503); worst bin prior 1.44 % at 5–10 ADC (24503). |
| Yields | n_kept identical; calibrated-block events arm/fixed 0.999–1.001; n_in_fit 0.999–1.023 per run (0.978–1.040 per split). |
| Cost | Fit 35.4 (bounded) / 30.1 (prior) µs per pulse vs mf 7.30 and WaveAnalyzer 7.71 µs per channel (`uniform/cost.txt`). Job wall time over mf: bounded +31 % / +34 %, prior +22 % / +28 % (25320 / 24503). |

χ²/dof, median, status-0 fits (`uniform/diag_summary.txt`):

| A (ADC) | 25320 fixed / bounded / prior | 24503 fixed / bounded / prior |
|---|---|---|
| 100–200 | 9.2 / 7.6 / 7.5 | 10.6 / 8.7 / 8.4 |
| 200–500 | 33 / 27 / 26 | 36 / 31 / 29 |
| 500–1000 | 128 / 100 / 98 | 134 / 104 / 101 |
| 1000–2000 | 331 / 279 / 265 | 414 / 333 / 318 |
| 2000–4000 | 448 / 411 / 392 | 731 / 611 / 576 |

σ/E (%), full splits, Δ vs fixed paired (`uniform/resolution_table.md`):

| Run (fit) | Path | fixed | bounded | Δ bounded | prior | Δ prior |
|---|---|---|---|---|---|---|
| 25320 (Gauss) | legacy | 2.136 | 2.162 | +0.026 ± 0.007 | 2.098 | −0.038 ± 0.006 |
| | gated | 2.133 | 2.159 | +0.026 ± 0.007 | 2.096 | −0.037 ± 0.005 |
| 24503 (CB) | legacy | 1.825 | 1.906 | +0.081 ± 0.003 | 1.824 | −0.002 ± 0.004 |
| | gated | 1.825 | 1.905 | +0.080 ± 0.003 | 1.823 | −0.002 ± 0.003 |

**Recommendation.** Drop the bounded arm: worse at both energies, and about half of its fits end on a bound, mostly sub-100 ADC neighbour pulses where the clipped shape biases the integral. Keep `--mf-shape prior` as an option for one full-scan comparison (−0.04 pp at 2.24 GeV, null at 3.49 GeV, +22–28 % run time). χ²/dof stays at 100–730 at high amplitude with τ free, so the limit is the two-τ functional form; a sampled per-module template is the next step (§9).

### 6.1 Template integral vs raw sum

Q_template = A·I_int (κ removed) against the WaveAnalyzer integral Q_raw, both ADC·ns, same single-peak, unclipped, in-window pulses; 20k events per split (025320 and 024503 × 3 splits, 024590 × 1), uniform weights. Median ratio [16–84 %]. `output/mf_fit_study_20261023/int_vs_rawsum/int_vs_rawsum.md` and the `ratio_vs_amp_`, `ratio_hist_`, `module_ratio_`, `q2d_<run>.png` plots.

| Run | A (ADC) | n | fixed | prior | bounded (median) |
|---|---|---|---|---|---|
| 25320 | 5–20 | 280101 | 1.113 [1.007, 1.260] | 1.105 [1.011, 1.233] | 1.088 |
| | 20–100 | 335475 | 1.009 [0.961, 1.097] | 1.008 [0.977, 1.084] | 1.003 |
| | 100–300 | 69679 | 0.978 [0.954, 1.000] | 0.986 [0.974, 0.999] | 0.984 |
| | 300–1000 | 47377 | 0.980 [0.961, 0.997] | 0.987 [0.978, 0.998] | 0.984 |
| | 1000–4000 | 32009 | 0.993 [0.978, 1.001] | 0.995 [0.983, 1.000] | 0.990 |
| 24503 | 5–20 | 359508 | 1.133 [1.017, 1.291] | 1.119 [1.019, 1.244] | 1.103 |
| | 20–100 | 416998 | 1.019 [0.969, 1.121] | 1.012 [0.980, 1.097] | 1.006 |
| | 100–300 | 81796 | 0.984 [0.960, 1.008] | 0.989 [0.976, 1.001] | 0.985 |
| | 300–1000 | 55893 | 0.983 [0.966, 1.000] | 0.988 [0.979, 0.999] | 0.984 |
| | 1000–4000 | 49647 | 0.995 [0.976, 1.004] | 0.995 [0.982, 1.000] | 0.986 |
| 24590 | 5–20 | 39082 | 1.080 [0.992, 1.185] | 1.074 [0.998, 1.172] | 1.058 |
| | 20–100 | 51574 | 1.004 [0.959, 1.071] | 1.003 [0.975, 1.061] | 0.999 |
| | 100–300 | 12950 | 0.982 [0.956, 1.004] | 0.988 [0.974, 1.000] | 0.985 |
| | 300–1000 | 11118 | 0.994 [0.976, 1.004] | 0.994 [0.982, 1.000] | 0.989 |

| Reading | Detail |
|---|---|
| Raw-sum low-amplitude loss | Ratio > 1 % above its 300–1000 ADC value below 50 / 100 / 200 ADC (25320 / 24503 / 24590); 8–13 % at 5–20 ADC. |
| Template deficit | At 100–1000 ADC the fixed template integral is 1–2 % below the raw sum (0.978–0.994), amplitude-dependent, returning towards 0.99–1.00 at 1000–4000 ADC. κ = 1.010 (20261022 note §3) corrects the pooled average only. |
| Per-pulse spread | 16–84 % half-width at 100–300 ADC: ±2.3 % (fixed) vs ±1.2–1.3 % (free arms). |
| Per-module spread | Module medians at A > 300 (≥ 20 pulses): RMS 0.60 / 0.50 / 0.47 % (25320, fixed / prior / bounded), 0.65 / 0.51 / 0.48 % (24503). |

## 7. Pile-up examples and rates

`plot_pileup_examples.py`, outputs in `output/pileup_examples_20261023/` (`pileup_rate_<run>.txt`, `pileup_candidates_<run>.tsv`, `pileup_examples_<run>.pdf` with a cover page and 10 candidates, `logs/`).

| Item | Rule |
|---|---|
| Events | Physics trigger (SSP0, bits 0x100); no cluster requirement. |
| Channels | PbWO4 only, unclipped. mf with nsigma 5, min_amp 5, kmax 4, `MIN_SEP_NS` 12, r_model = `--secondary-frac` (0.10). Templates: h50 column 25320; X17 runs the combined file (column 27000). E = A·I_int·κ·factor, no LMS. |
| Denominator | Channels whose largest pulse has E_dom > 0.5 E_beam (`--dominant-frac`); channels with an onset before ped_n/2 or a dominant at the window end are excluded and counted on the PDF cover. |
| Secondary | Any other pulse with E/E_dom > 0.10, onset not at the window end. Candidate = denominator channel with ≥ 1 secondary. |
| Scan | First split complete, then further splits until `--n-examples` (10) candidates. |
| Plots | Data with errors, each pulse component, their sum, the pre-window band, residuals; table of t, A, E, E/E_d, A/A_d, Δχ² per secondary (χ² without that pulse minus full fit). |

| Run | Type | Events | Candidates / denominator | Rate |
|---|---|---|---|---|
| 24386 | 3.49 GeV H2 | 118418 | 10 / 81910 | 0.012 % |
| 24503 | 3.49 GeV H2 | 159005 | 13 / 100739 | 0.013 % |
| 24512 | 3.49 GeV H2 | 307949 | 10 / 187583 | 0.005 % |
| 24590 | 0.73 GeV H2 | 197569 | 95 / 122853 | 0.077 % |
| 24603 | 0.73 GeV H2 | 194246 | 99 / 101886 | 0.097 % |
| 25228 | 2.24 GeV H2 | 154996 | 202 / 91793 | 0.220 % |
| 25257 | 2.24 GeV H2 | 160756 | 230 / 95329 | 0.241 % |
| 25260 | 2.24 GeV H2 | 162677 | 213 / 96197 | 0.221 % |
| 25299 | 2.24 GeV H2 | 161217 | 173 / 95329 | 0.182 % |
| 25308 | 2.24 GeV C | 163358 | 196 / 106587 | 0.184 % |
| 25320 | 2.24 GeV H2 | 159596 | 172 / 94411 | 0.182 % |
| 25800 | 2.24 GeV Ta (X17) | 12066 | 111 / 10726 | 1.035 % |
| 26138 | 2.24 GeV Ta (X17) | 5357 | 57 / 5007 | 1.138 % |
| 26230 | 2.24 GeV Ta (X17) | 4302 | 48 / 3839 | 1.250 % |

Secondaries pooled from the candidate TSVs (Δt = t_sec − t_dom):

| Sample | Candidates / secondaries | Δt 16 / 50 / 84 % (ns) | Δt < 0 | A ratio 16 / 50 / 84 % | ≥ 2 / ≥ 3 secondaries |
|---|---|---|---|---|---|
| All | 1629 / 2344 | −36.8 / 38.3 / 88.3 | 29 % | 0.13 / 0.30 / 0.69 | 514 / 201 |
| 2.24 GeV H2 and C | 1186 / 1809 | −32.1 / 36.8 / 84.5 | 27 % | 0.13 / 0.33 / 0.71 | 451 / 172 |
| 3.49 and 0.73 GeV | 227 / 315 | −38.4 / 42.6 / 90.6 | 31 % | 0.12 / 0.22 / 0.58 | 59 / 29 |
| Ta (X17) | 216 / 220 | −76.3 / 40.3 / 107.7 | 40 % | 0.12 / 0.22 / 0.75 | 4 / 0 |

- 41 % of secondaries lie inside the 100–200 ns pre-window; 24 % within |Δt| < 30 ns. 26230 has median Δt −19.9 ns and 54 % with Δt < 0.
- **Several similar pulses.** 503 of 1629 candidates reach the kmax = 4 cap and 201 have ≥ 3 secondaries, mostly at 2.24 GeV, some with comparable ratios (25257 event 16423858 W212: 0.99, 0.57, 0.88; 25257 16340095 W172: 0.60, 0.77, 0.997). Modules W950, W1014, W480, W1025 and W257 recur; W257 shows secondaries near +16, +33 and +80 ns at ratios ≈ 0.25, 0.16, 0.13. Almost absent on Ta (4 candidates with ≥ 2).
- **χ²/dof of large pulses.** The mf joint-fit χ²/dof of candidates grows roughly as A_dom²: median 78 (A_dom < 500 ADC), 107 (500–1000), 720 (1000–2000), 2984 (2000–4000), 66534 (> 4000, 12 candidates). This matches the template misfit of §6.

## 8. Clustering as configured

| Constant | Value | Source |
|---|---|---|
| min_module_energy | 1.0 MeV (C++ struct default 0) | `database/reconstruction_config.json:15`; `prad2det/include/HyCalCluster.h:36` |
| min_center_energy | 10 MeV | JSON:16 |
| min_cluster_energy | 50 MeV | JSON:17 |
| min_cluster_size | 1 | JSON:18 |
| corner_conn | false (4-connected island) | JSON:19 |
| split_iter / least_split | 6 / 0.01 | JSON:27–28 |
| log_weight_thres T | 3.6 (seed 3 × 3 centroid) | JSON:29 |
| seed_time_window | JSON 4.0 ns; C++ default −1 (off); benchmark legacy −1, gated `--window` (default 8, working point 16) | JSON:34; `HyCalCluster.h:76`; benchmark :1473, :1963–1964 |
| non_linear_corr | true (struct default in the benchmark; the `_fadc` files carry no `nl1`/`nl2`, so the factor is 1) | JSON:20; `HyCalCluster.h:50` |
| leakage_correction | false | JSON:21 |
| Profile | Benchmark: SimpleProfile (`set_profile` never called). Production: Geant4 `cluster_profiles/prof_pwo.dat`, `prof_lg.dat` | `prad2det/src/HyCalCluster.cpp:45`; `prad2det/src/PipelineBuilder.cpp:320–328` |
| WaveAnalyzer peak | max(5σ_ped, 10 ADC); min_peak_ratio 0.3; int_tail_ratio 0.1 | `database/daq_config.json:103–106`; `prad2dec/include/WaveAnalyzer.h:55–62` |

| Item | Detail |
|---|---|
| Legacy vs gated | The benchmark builds both clusterers from 8 copied fields (energies, size, corner_conn, split, T); non-linearity, leakage and profile take struct defaults. Legacy: one hit per module (largest pulse in the pre-window), no time cut. Gated: every pre-window pulse; each neighbour contributes its largest pulse within ±W of the seed. |
| Scanned | Seed time window W and pre-window LO/HI (20261022 note §2, §7); calibration method (§3–4); neighbour amplitude definition (§5). |
| Not scanned | min_module_energy (Phase 4, `--cluster-min-module-energy`, not implemented yet), min_center_energy, min_cluster_energy, split_iter, least_split, T, corner_conn, profile (SimpleProfile vs Geant4), cluster size (3 × 3 / 5 × 5 vs island). |

Documentation inconsistencies in `docs/technical_notes/hycal_clustering/hycal_clustering.md`:

| # | Note | Files | Status |
|---|---|---|---|
| 1 | Lines 372–373: "W = 16 ns the production default"; line 528: default −1 | JSON:34 has 4.0; C++ default −1; benchmark default 8 | Real; also 4 ns is below the note's own finding that ±8 ns clips the signal core (line 367). |
| 2 | Line 58: min_module_energy "default 0 MeV"; line 520: 1.0 MeV | Struct 0, JSON 1.0; lines 514–515 say the table matches the JSON | Real; line 58 gives the struct default, not the configured value. |
| 3 | Step 3 (line 111): kSplit set when a hit is shared by two clusters | `HyCalCluster.cpp:324` sets kSplit on every cluster emitted by the splitter that keeps a hit | Real; line 561 ("came out of the splitter") matches the code. |

## 8a. Inner-acceptance fix and module-by-module resolution (2026-10-23/24)

### 8a.1 The bug and the fix

| Item | Detail |
|---|---|
| Bug | `_select_pool` compared `--exclude-rowcol 15 20 15 20` with the C++ 0-based Module.row/column. The beam hole (W561/W562/W595/W596) is rows/cols 16–17 0-based, 17–18 1-based. |
| Old rectangle | 0-based 15–20: one ring on the top/left side of the hole (row/col 15), three on the bottom/right (18–20); x [−41.5, 83.1] mm, y [−83.0, 41.5] mm. Equals `16 21 16 21` in the new reading. |
| Fix | 1-based inclusive (compares row+1/col+1). `15 20 15 20` = centred 6 × 6, 32 modules, two rings on every side, x and y in [−62.3, 62.3] mm. Edges: row 15 W491–W496 (top), row 20 W661–W666 (bottom), col 15 W491…W661, col 20 W496…W666. The driver's hard-coded `15 20 15 20` now gives the centred rectangle; `plot_module_resolution.py` uses the same convention (pool and dashed box). |
| Pool swap (025320.00050) | 11 out (W491–W496, W525, W559, W593, W627, W661), 11 in (W531, W565, W599, W633, W667, W696–W701). |
| Byte identity, no rectangle | Pre-edit vs edited, 30k events of 025320.00050, mf + pooled gains + cluster TSV: window, cluster and LMS TSVs and all PNGs identical; `result.tsv` differs only in `elapsed_s` and the template path (same sha1). |
| Provenance | `percrystal_calib` tables and the benchmark percrystal table carry `rowcol_origin=1`. |
| Still 0-based | `benchmark_hycal_timing.py` (reference, protected), `benchmark_hycal_timing_deconv.py`. |
| Earlier outputs | Scans up to `scan_20261023_pc*` and the tables in `output/percrystal/` used the old reading. D1 (§2) superseded (`plans/plan_20261023.md`, "D1 superseded", "Centred exclusion rectangle"). |

### 8a.2 Gains refit

| Item | Detail |
|---|---|
| Tables | `output/percrystal_rect/<run>/gains_<run>_all.tsv` (+ `.json`; log `output/percrystal_rect/fit_<run>_all_fixedwin.log`); same dumps and flags as `output/percrystal/` (`--splits all --fixed-windows`). All 14 runs converged in 2 iterations, ~1 min total. |
| Pool | 25320 556k → 242k, 24503 644k → 333k, 25800 50k → 30k events. The high-rate ring 2 around the hole at top/left (r ≈ 52 mm) leaves; ring 3 around the hole at bottom/right (r ≈ 73 mm; ring 1 of §8a.4) enters. |
| Δg, modules free in both fits | Median \|Δg\| 0.03–0.10 % (26230 0.32 %); RMS 0.86–1.05 % at 2.24 and 3.49 GeV, 0.13–0.41 % at 0.73 GeV and Ta. |
| 25320 by ring | Rectangle edge ring median −1.1 %, RMS 2.6 %; first ring outside RMS 2.1 %; second ring outward 0.07–0.25 %. |

- W495 (newly excluded): −4.0 / −7.8 / −4.7 % (25320 / 24503 / 25308); partly prior-bound, prior_frac 0.07–0.20 at H2/C.
- W697 (newly included): +4.8 to +8.3 % across the 2.24/3.49 GeV runs; becomes free (prior_frac 0.11–0.46 → ≤ 0.02).
- W598 (inside both rectangles): prior-bound → free (prior_frac 0.97 → 0.07–0.11), −9.2 / −9.8 / −9.7 %; largest \|Δg\| in 9 of 14 runs.
- 0.73 GeV changes ≤ 1.7 % (W700). Hot-spot neighbours: W461 +0.8/+1.4 %, W462 +0.5/+1.0 %, W427 −0.2 %.

### 8a.3 Old → centred rectangle

Scan `output/scan_20261024_pc`: driver `--pulse-method mf --window 16 --pre-window 100 200 --gain-file-map auto --gain-dir output/percrystal_rect --write-cluster-tsv --parallel 6 --timeout-seconds 1800`; 119/119 splits, 0 failures, 1 h 45 min. Old = `output/scan_20261023_pc`. Post σ_E/E (%), legacy, preferred fit, mean over splits; Δ (pp) paired per split, SE and SD over splits; `compare/rect_old_vs_centred.tsv`.

| Run type (fit) | Splits | Old | Centred | Δ ± SE | SD of Δ | Modules / split | Events / split |
|---|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 28 | 1.826 | 1.823 | −0.003 ± 0.010 | 0.052 | 108.9 → 108.9 | 42.5k → 20.7k |
| 2.24 GeV H2 (Gauss) | 42 | 2.132 | 2.071 | −0.061 ± 0.002 | 0.013 | 124.2 → 125.7 | 66.6k → 27.6k |
| 0.73 GeV H2 (Gauss) | 16 | 4.285 | 4.440 | +0.154 ± 0.019 | 0.076 | 105.9 → 106.2 | 26.6k → 28.4k |
| 2.24 GeV C (CB) | 8 | 2.003 | 1.966 | −0.037 ± 0.004 | 0.010 | 167.4 → 166.9 | 81.3k → 41.6k |
| Ta X17 (CB) | 25 | 2.395 | 2.082 | −0.313 ± 0.055 | 0.277 | 20.8 → 15.8 | 2.5k → 1.0k |

| Reading | Detail |
|---|---|
| Gated | Δ within 0.002 pp of legacy except Ta (−0.258 ± 0.054). |
| 2.24, 3.49 GeV and C | Sample change. The centred rectangle removes the high-rate top/left ring 2 (the busiest and worst modules; the old worst lists had W495, W496, W493, W559 there) and halves the calibrated-block events. Per-module medians move ≤ 0.02 pp (§8a.4). |
| 0.73 GeV | Worse by +0.15 pp. Ring 3 at bottom/right enters (W667, W700, W701 become free). Hypothesis, not tested: Møller contamination in the newly included inner modules. W667 is the 2nd worst module at 0.73 GeV (7.65 %). |
| Ta | ~16 modules and 1k events per split; too few to interpret. |

### 8a.4 Module-by-module resolution

| Item | Detail |
|---|---|
| Tool | `plot_module_resolution.py`. Input: a driver scan run with `--write-cluster-tsv` (`result.tsv`, `result_window.tsv`, `result_clusters.tsv` per split). |
| Pool | The benchmark calibration pool: path, seed r < `--inner-radius` 200, seed E ≥ min_seed_energy, outside `--exclude-rowcol` (1-based, default `15 20 15 20`), finite E_cal. Pre = cluster_E_raw, post = cluster_E_cal. |
| Fits | Gaussian and CB per seed module; `--fit auto` = CB for 3.49 GeV, C and Ta, Gaussian for 2.24 and 0.73 GeV H2. `--min-n 1000` (default): fitted only if N > min_n, else low_stat (hatched). |
| `--group-runs` | Pools the events of all runs of a type before the fits; per-run maps still written. |
| Outputs | Per run / per type: sigma_post/pre, resid_post/pre, N, sigma_hist PNGs, `spectra.pdf`, `modules.tsv`. `summary/`: `runs_summary.tsv`, `types_summary.tsv`, `worst_modules.tsv`, `corr_2p24_H2_r.tsv` + scatter. `--path gated` writes to `<out>_gated` (gated maps not run). |
| Runs | After: `output/module_resolution_20261024/` (`scan_20261024_pc`, `--group-runs`, centred). Before, like for like: `output/module_resolution_20261023_n1000/` (`scan_20261023_pc_clusters`, `--group-runs --exclude-rowcol 16 21 16 21`). First per-run maps: `output/module_resolution_20261023/` (`scan_20261023_pc_clusters`, per run, old rectangle, fitted modules down to N = 112). |

| Run type (before → after, `types_summary.tsv`) | Modules N > 1000 / with data | σ_M/μ_M median [16–84 %] (%) | Pooled σ_E/E (%) |
|---|---|---|---|
| 3.49 GeV H2 | 112/171 → 112/171 | 1.72 [1.59, 1.93] → 1.74 [1.59, 1.97] | 1.83 → 1.83 |
| 2.24 GeV H2 | 121/164 → 125/165 | 2.04 [1.94, 2.17] → 2.05 [1.94, 2.15] | 2.13 → 2.07 |
| 0.73 GeV H2 | 96/139 → 97/138 | 3.69 [3.36, 5.55] → 3.70 [3.36, 5.70] | 4.26 → 4.41 |
| 2.24 GeV C | 78/195 → 78/194 | 1.96 [1.85, 2.08] → 1.98 [1.84, 2.06] | 2.00 → 1.97 |
| Ta X17 | 15/42 → 8/41 | 2.33 [2.08, 2.43] → 1.93 [1.90, 2.26] | 2.30 → 2.08 |

| Run type | Worst 10 after, σ_M/μ_M (%), `worst_modules.tsv` |
|---|---|
| 3.49 GeV H2 | W461 2.49, W700 2.34, W692 2.31, W427 2.30, W633 2.28, W699 2.28, W489 2.27, W696 2.23, W462 2.21, W490 2.16 |
| 2.24 GeV H2 | W697 2.57, W461 2.36, W427 2.32, W533 2.28, W396 2.25, W727 2.22, W531 2.22, W590 2.21, W566 2.21, W360 2.20 |
| 0.73 GeV H2 | W490 8.03, W667 7.65, W497 7.41, W426 7.00, W427 6.63, W566 6.36, W425 6.35, W462 6.33, W700 6.25, W600 6.21 |
| 2.24 GeV C | W697 2.36, W461 2.25, W427 2.25, W533 2.23, W625 2.22, W422 2.17, W728 2.17, W532 2.15, W767 2.14, W667 2.13 |
| Ta X17 | only 8 fitted (W697 2.41 worst) |

- Per run, after, N > 1000 / with data: 24386 30/127, 24503 77/171, 24512 52/159, 24590 84/138, 24603 60/118, 25228 70/145, 25257 49/149, 25260 52/149, 25299 52/153, 25308 78/194, 25320 52/147; 25800, 26138, 26230 none. Ta usable only pooled (8 of 41 modules).
- Hot spot W461, W427, W462 remains (outside both rectangles). W461 and W462 less extreme: W461 2.92 → 2.49 (3.49 GeV), 2.57 → 2.36 (2.24 GeV), 2.46 → 2.25 (C); W462 2.49 → 2.21 (3.49 GeV). W427 +0.01 to +0.06 pp. W495/W496 (old worst at 3.49 GeV 2.43/2.54, 2.24 GeV 2.34/2.34) now excluded. Newly included modules in the worst 10: W697, W700, W633, W667, W531 (2.24 GeV H2); also W699, W696 at 3.49 GeV.

Ring structure from the after `modules.tsv`. Ring k = k-th ring of modules outside the centred rectangle; ring 1 = first included ring, r 73–103 mm. Median σ_M/μ_M (%) of fitted modules (count):

| Run type | Ring 1 | Ring 2 | Ring 3 | Ring 4 |
|---|---|---|---|---|
| 3.49 GeV H2 | 1.86 (28) | 1.75 (36) | 1.72 (32) | 1.64 (16) |
| 2.24 GeV H2 | 2.05 (28) | 2.04 (36) | 2.04 (37) | 2.04 (24) |
| 2.24 GeV C | 2.00 (28) | 1.97 (32) | 1.96 (18) | – |
| 0.73 GeV H2 | 6.07 (12) | 3.88 (36) | 3.59 (32) | 3.37 (17) |

Ring 1 by side (6 non-corner modules per side):

| Run type | Overall median | Top | Bottom | Right | Left |
|---|---|---|---|---|---|
| 3.49 GeV H2 | 1.74 | 1.81 | 2.12 | 1.88 | 1.62 |
| 2.24 GeV H2 | 2.05 | 2.11 | 2.08 | 2.11 | 1.99 |
| 2.24 GeV C | 1.98 | 2.00 | 2.05 | 2.03 | 1.84 |

| Reading | Detail |
|---|---|
| Ring 1 | Above the overall median on top, bottom and right; the left side at or below it at all three. Excess over ring 2: +0.11 pp (3.49 GeV), +0.01 (2.24 GeV H2), +0.03 (C), +2.2 pp (0.73 GeV, 12 modules). Worst module at 0.73 GeV: W490 (1-based row 15, col 14, left side). |
| Radial, outside ring 1 | 3.49 GeV ring medians fall 0.11 pp from ring 2 to ring 4 (1.75 → 1.64; ring 4 = 16 modules). 2.24 GeV H2 and C flat within 0.03 pp. 0.73 GeV falls 3.88 → 3.37. |
| Banding, outside ring 1 (3.49, 2.24 H2, C) | SD of row medians and of column medians 0.04–0.11 pp, below the module SD (0.10–0.19 pp); no row/column banding. |
| W428, W494 | W428 normal: 1.71 / 1.99 / 1.85 % (3.49 / 2.24 H2 / C; rank 65 of 112, 90 of 125, 65 of 78). W494 now inside the rectangle (1-based row 15), excluded; before maps 1.73 % (3.49) and 2.19 % (2.24 H2). |
| Run-to-run, 2.24 GeV H2 | First per-run maps (`corr_2p24_H2_r.tsv`): Pearson r of per-module σ_M/μ_M between runs 0.46–0.60. Module-to-module SD per run 0.14–0.17 pp (mean 0.15); the between-run covariance gives a reproducible part of 0.11 pp and a statistical part of 0.10 pp. N > 1000 grouped maps: 0.78–0.93 (before), 0.65–0.92 (after). |

### 8a.5 Decision and paths

- The centred two-ring rectangle (`--exclude-rowcol 15 20 15 20`, 1-based) is the reference acceptance from here on.
- `output/scan_20261024_pc` is the baseline for the spline-template comparison (§9).
- The inner-ring degradation stays in the sample and is reported as a σ-vs-radius statement (ring table, §8a.4).
- Paths: `output/percrystal_rect/`, `output/scan_20261024_pc/` (+ `compare/rect_old_vs_centred{,_splits}.tsv`), `output/module_resolution_20261024/{3p49_H2,2p24_H2,0p73_H2,2p24_C,Ta_X17,summary}/` (+ per-run dirs), `output/module_resolution_20261023_n1000/`, `output/module_resolution_20261023/` (+ `run.log`), `plans/plan_20261023.md` ("D1" superseded entry; "Centred exclusion rectangle").

## 9. Interpretation and next step (in preparation; results to be appended)

| Finding | Evidence |
|---|---|
| The two-τ template is 1–2 % short of the raw sum for mid and large pulses, amplitude-dependent | §6.1; χ²/dof 100–730 at A > 500 ADC with fixed or free τ (§6); χ²/dof ∝ A² in pile-up candidates (§7). |
| Freeing τ narrows the per-pulse spread but does not fix the shape | §6: −0.04 pp at 2.24 GeV (prior), null at 3.49 GeV. |
| Per-crystal gains are the one constant-term lever that worked | §3–4: −0.05 to −0.09 pp at 2.24 and 3.49 GeV and C. |
| Timing is not a lever on e-p | 20261022 note §7 (Δgate ≤ 0.004 pp on H2 and C); seed-time amplitudes fail (§5). |
| Pile-up rates explain why only Ta gains from the gate | Rates 1.0–1.25 % on Ta vs 0.005–0.24 % on H2 and C (§7); Ta Δgate −0.099 ± 0.024 pp (20261022 note §7). |

**Decision (user).** Build sampled per-module templates by spline fit, one set per run type, following Wassim's ROOT macros `best_ref_shape_fit.C` and `hycal_mf_fit_premodule_visual.C` (`/Users/jingyi.zhou/research/PRad-II/Waveform_Wassim/`). Pipeline research: `docs/research-spline-template-pipeline.md` (in progress).

| Planned deliverable | Status |
|---|---|
| ROOT template file and per-module fit PDF, per run type | done: `output/templates_spline/<slug>/templates_<slug>.{npz,root,pdf}` (§9.1) |
| Pile-up rate heat maps (per module) | done: `output/pileup_20261024_fix/` (§9.3) |
| Combined pile-up example PDF | done: `output/pileup_20261024_fix/pileup_examples_all_runs.pdf` |
| Detection efficiency table in the style of Table 4 of `docs/technical_notes/waveform_analysis/fit_pulse_template.md`, from the pile-up generator | done: `output/pileup_eff_20261024_fix/<slug>/eff_<slug>.md` (§9.4) |
| Calibration harness rerun with the sampled templates (mf baseline → + per-crystal) | done: `output/scan_20261024_spline/` (§9.5) |

### 9.1 Results (2026-10-24): sampled per-module templates

**Templates.** `spline_template_extract.py` (select / build / ampcorr / check; plan `../plans/plan_20261024.md` §6 and its progress log). Clean single pulses from all 119 local splits (exact trigger 0x100; Ta any physics trigger 0xf00), our vetoes only (clip, one mf pulse above 10 % of the main one, onset ≥ 60 ns, first threshold ≥ sample 20, contained, pedestal RMS/slope, tail), PASS 1's shape cut dropped (user, Checkpoint 1). Per module: up to 200 phase-aligned pulses, 3 alignment iterations (worst 20 % by χ² dropped), 0.5 ns bins, penalised cubic B-spline (λ = 10⁴, chosen by a held-out phase-quartile rule on the earlier PASS-1-shape-cut pools, whose sampling phase was concentrated, Rayleigh R ≈ 0.75; the final pools without that cut give R ≈ 0.13 on the aligned t0, and G1-5 still passes), unit peak, onset at 0.5 % of the peak, slow tail extended exponentially to 400 ns (the real tail is ≈ 2.5 % of the peak at 160 ns, τ_slow ≈ 125 ns). Energy integral: the phase average of the WaveAnalyzer walk-out on the template (no sub-sample step). Amplitude correction E = A·I_int·f(A): f from noise-free band-template waveforms (100–400, 400–1000, 1000–2500 ADC) fitted with the module template, chained within-module steps, PCHIP in ln A, constant outside the knots, per run type (the five run types are not consistent within 2 SE).

| Run type | fallback 0 / 1 / 2 | f(A) knots (f @ A ADC) | G1-7 Q_template/Q_raw, 100–300 / 300–1000 / 1000–2500 ADC, with f (without) |
|---|---|---|---|
| 3.49 GeV H2 | 858 / 268 / 26 | 1.0052@201, 0.9949@537, 0.9973@1566 | 0.996 / 0.994 / 0.999 (0.993 / 0.998 / 1.003) |
| 2.24 GeV H2 | 1062 / 89 / 1 | 1.0068@196, 0.9976@519, 0.9936@1550 | 0.999 / 0.996 / 0.998 (0.993 / 0.998 / 1.004) |
| 0.73 GeV H2 | 577 / 543 / 32 | 1.0034@191, 0.9962@507 | 0.995 / 0.997 / – (0.993 / 0.999 / –) |
| 2.24 GeV C | 742 / 310 / 100 | 1.0043@193, 0.9984@592, 0.9911@1592 | 0.997 / 0.996 / 0.996 (0.993 / 0.996 / 1.005) |
| Ta X17 | 812 / 297 / 43 | 1.0032@165, 0.9988@583, 0.9936@1365 | 0.998 / 0.998 / 0.999 (0.995 / 0.998 / 1.005) |

Gates (`output/templates_spline/<slug>/check/gates.tsv`): G0 (two-τ default byte-identical) pass; PASS 1 port vs the patched macro on 24386 (0 mismatches over 227 908 waveforms; best reference agrees for 798/798 modules; the macro's best waveform is among our candidates for every module that passes our vetoes); held-out χ²/dof 0.16–0.20 / 0.06–0.07 / 0.03 of the two-τ value at 100–300 / 300–1000 / 1000–2500 ADC; sub-sample rise residual ≤ 0.3 % in every phase quartile; phase-averaged C++ walk-out integral within 0.22 % of the bank's I_int (3–5 modules per run type have no single C++ peak on the noise-free template); G1-7 bin-to-bin spread with f(A) 0.10–0.53 % (without 0.67–1.14 %), not flat by the 2-SE rule (SE ≈ 0.02–0.03 %); amplitude bands differ beyond noise (the reason for f(A)); all files load in a PASS 2-style ROOT reader. Gate list: G0 PASS, G1-1 PASS, G1-2 PASS, G1-3 template-vs-PASS-1 RMS 0.65 % PASS / χ²/ndf (σ = rms) 75 FAIL, G1-4 PASS, G1-5 PASS, G1-6 FAIL (only the no-single-peak modules), G1-7 FAIL (2-SE flatness), G1-8 PASS (2.24 H2, 3.49 H2, C), G1-9 FAIL, G1-10 PASS, G2-1 PASS (172/94411 reproduced), G2-2 FAIL on the new pseudo-data (28/45).

### 9.2 Integral and χ² on the study pulses (§6.1 and §6 tables, sampled arm)

`mf_fit_study.py intraw/diag --arm sampled` (`../output/mf_fit_study_20261023/int_vs_rawsum/int_vs_rawsum_sampled.md`, `uniform/diag_summary_sampled.txt`): median Q_template/Q_raw at 100–300 / 300–1000 / 1000–4000 ADC, two-τ fixed → sampled: 25320 0.978 / 0.980 / 0.993 → 0.997 / 0.991 / 0.998; 24503 0.984 / 0.983 / 0.995 → 0.999 / 0.990 / 1.001; 24590 0.982 / 0.994 / – → 0.991 / 0.995 / –. Spread of module medians at A > 300 (25320): 0.60 % → 0.41 %; the per-pulse 16–84 % half-width stays ≈ ±2.4 % (100–300 ADC; fixed ±2.3 %, free-τ arms ±1.2–1.3 %). The sampled Q_template includes f(A) (§6.1 defines Q_template = A·I_int). mf χ²/dof (σ² + (0.03·model)²) at A > 500 ADC on 25320: 233 (fixed) → 1.08 (sampled).

### 9.3 Pile-up rates and the +12–40 ns excess

With the 50 MeV denominator the rates rose to 0.23–0.60 % (e-p, C) and 10.6–11.2 % (Ta); the level change is the denominator (the two-τ bank gives 0.53 % on the same events). Study `../output/pileup_20261024/artefact/` (`pileup_artefact_study.py`): 37–70 % of the e-p candidates sat at dt = +12–40 ns. Not a template-tail error: the clean single-pulse residual averaged over [12, 40) ns is within ≈ 0.4 % of the peak; per ns it reaches +2.5 / +4 % at τ ≈ 20 ns for A < 150 ADC (2.24 / 3.49 GeV, with negative lobes around it) and ≤ 0.9 % at 400–1000 ADC. Two populations: burst events with ≥ 10 candidate channels (34–1149 channels, dominant at 186–275 ns, χ²/dof ~10³; 25 % of 2.24 GeV candidates) and isolated real second bumps at dt ≈ 23 ns with E_ratio ≈ 0.18, amplitude-proportional, in ≈ 0.1–0.3 % of pulses (≈ 0.1 % at 3.49 GeV), mostly single-channel (origin not identified). Fix: burst veto (≥ 10 candidate channels, removed from numerator and denominator) and a window rule (a secondary with 0 ≤ dt < 40 ns and E_ratio < 0.16 is not counted; 0.17 already loses real ratio-0.2 pile-up). Rates after the fix (`../output/pileup_20261024_fix/pileup_rate_summary.tsv`, heat maps, examples PDF, `pileup_hot_channels.tsv`): 3.49 GeV 0.15–0.23 %, 2.24 GeV H2 0.40–0.41 %, 0.73 GeV 0.26–0.29 %, C 0.36 %, Ta 10.3–11.0 %; hot channels (W335, W1, W1019, W1020, W934; W934 a sporadic +25 ns, ratio 0.29 feature at 0.73 GeV) change the rates by ≤ 0.011 pp. The fixed rate still contains second bumps with E_ratio ≥ 0.16 and is an upper bound on accidental pile-up.

### 9.4 Detection efficiency (Table 4 style)

`benchmark_pileup_mf.py` on sampled-template pseudo-data (`../output/pileup_eff_20261024_fix/<slug>/eff_<slug>.md`; 50 bases per cell, generator dt = injected peak − base 50 % CFD time, true onset separation ≈ dt − 6.4 ns; ±4 ns onset match). At r_model 0.3 (calibration working point): ratio ≥ 0.5 found 100 % from dt = 32 ns (Ta 96 %), ratio 2 from dt = 20 ns; ratio ≤ 0.3 only at dt = 128 ns (the acceptance rule applies inside the 96 ns kernel span); fakes 0 except 0.1–0.8 per event at dt 16–24 ns, ratio 1 (Ta: a 0.06 per-event baseline). At r_model 0.10: ratio ≥ 0.2 found 96–100 % from dt 24–32 ns. The window rule changes eff_match by 0 points in every cell with ratio ≥ 0.2, dt ≥ 24 ns. The WaveAnalyzer "≥ 2 peaks" column reproduces Table 4 in 45/45 cells on the original Table 4 dataset and in 28/45 on the new pseudo-data (different injected shape and bases). Pseudo-data fake rates are optimistic: the injected shape is the reconstruction template.

### 9.5 Calibration with the sampled templates

Harness dump pass (`../output/scan_20261024/spl_percrystal/`, κ = 1.0, f(A) in every energy path), pooled per-crystal gains refitted (`../output/percrystal_spl/`; 14/14 converged; vs `percrystal_rect`: median ratio 0.998–1.004, per-module RMS 0.3–0.5 % for modules used in both fits, 0.16–0.40 % over all modules incl. prior-bound ones), full driver scan `../output/scan_20261024_spline/` (119/119 splits, cluster TSVs), paired against `../output/scan_20261024_pc/` (`compare/arms.md`). σ/E (%), preferred fit, legacy path, mean over splits; A = two-τ + pooled gains + seed-indexed, A′ = A without the seed-indexed step, B = sampled + refitted gains + seed-indexed, C = B without it:

| Run type | A | A′ | B | C | B − A (pp) | chosen |
|---|---|---|---|---|---|---|
| 3.49 GeV H2 (CB) | 1.823 | 1.859 | 1.709 | 1.733 | −0.114 ± 0.006 | B |
| 2.24 GeV H2 (Gauss) | 2.071 | 2.093 | 1.952 | 1.975 | −0.119 ± 0.002 | B |
| 0.73 GeV H2 (Gauss) | 4.440 | 4.354 | 4.366 | 4.289 | −0.074 ± 0.021 | B (B − C +0.077 ± 0.041, < 2 SE) |
| 2.24 GeV C (CB) | 1.966 | 1.987 | 1.858 | 1.883 | −0.108 ± 0.005 | B |
| Ta X17 (CB) | 2.082 | 2.194 | 2.000 | 2.099 | −0.082 ± 0.029 | B (report only) |

Gated path within 0.03 pp of these (Ta −0.107). Yields: ε_yield 0.994–1.001, ε_peak 1.000–1.013. Template success at 2.24 and 3.49 GeV H2 (plan §8.3). Per-module maps (`../output/module_resolution_20261024_spline/`, `--group-runs`, N > 1000) vs `module_resolution_20261024`: per-module median σ/E post 3.49 GeV 1.74 → 1.61 %, 2.24 GeV H2 2.05 → 1.92 %, 0.73 GeV 3.70 → 3.60 %, C 1.98 → 1.83 %, Ta 1.93 → 1.99 % (8 modules).

Caveats: in-sample gains as before; f(A) leaves an amplitude-dependent residual of 0.1–0.5 % on the template pulses (G1-7) and 0.7–1.1 % between 300–1000 and 1000–4000 ADC on the §6.1 study pulses (§9.2); at 0.73 GeV "pooled gains only" is lower than with the seed-indexed step (4.289 vs 4.366) but within 2 SE; Ta has ≈ 16 calibrated modules per split.

## 10. Caveats and open items

| Item | Detail |
|---|---|
| In-sample gains | Held-out tested on S_dev only (§3.2, §4). |
| Even/odd tables | Fitted without `--fixed-windows`; 24503 odd not converged; 25320 even/odd have 276 free modules (284 in `_all`). |
| 0.73 GeV | Per-crystal gain not robust there (§3.2, §4). |
| Ta / X17 | Reported only (§4). New X17 files expected next week. |
| Free-τ figures | Single-run in places: wall time and χ²/dof differ between 25320 and 24503 (§6 tables). |
| Pile-up rates | Scan length differs by run (one full split for 0.73 and 2.24 GeV and 24503; a full plus a partial split for 24386 and 24512; 4–12k events for Ta). Denominators of 25257 and 25299 are both 95329 (event and candidate counts differ); unchecked coincidence. Excluded edge channels are not in the rate. |
| Phases 2 and 4 | Position correction and quick checks (threshold scan, window sums, rectangle, saturation percentile) not run. |
| Planning record | `../plans/plan_20261023.md` §2.1 lists 1.96 for "mf legacy, Gaussian for all" at 3.49 GeV; the file gives 1.943 (used here). Its Results section lists 268 free modules for 24590; the current `_all` table has 266. |

## 11. How to run

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh
E=/Users/L00338853/work/PRad/data/evio; DB=$PRAD2_DATABASE_DIR/calibration
# --exclude-rowcol is 1-based inclusive: 15 20 15 20 = centred 6 x 6 rectangle (§8a.1)
WP="--calibrate --exclude-rowcol 15 20 15 20 --pre-window 100 200 --inner-radius 200 --min-per-module 30 --window 16 --pulse-method mf"

# (a) per-crystal: dump the calibration pool of one split (in-split fit as a side product), then the pooled fit
python3 benchmark_hycal_timing_test_20261021.py $E/prad_025320.evio.00050 output/<tag>/025320_00050 $WP \
  --lms-corr auto --hc-calib-file $DB/calibration_factor_2p2_fadc.json --cal-method percrystal --cluster-dump auto
python3 percrystal_calib.py fit --dumps "output/scan_20261023/mf_percrystal/scan_025320/split_*/result_dump.npz" \
  --splits all --fixed-windows --out output/percrystal/025320/           # --splits even|odd for held-out tables
python3 percrystal_calib.py stability --gains "output/scan_20261023/mf_percrystal/scan_0252*/split_*/result_gains.tsv" \
  "output/scan_20261023/mf_percrystal/scan_025320/split_*/result_gains.tsv"
python3 percrystal_calib.py --self-test

# (b) pooled-gain full scan through the driver (as run for output/scan_20261023_pc)
python3 benchmark_hycal_timing_scan.py --pulse-method mf --window 16 --pre-window 100 200 --gain-file-map auto \
  --parallel 12 --timeout-seconds 1800 --output-dir output/scan_<tag>

# (c) free-τ study: one split per arm (024503: calibration_factor_fadc_3p5.json, --lms-corr off), then the analysis
R=025320; S=00050; O=output/mf_fit_study_<tag>/uniform
for arm in fixed bounded prior; do X=""; [ $arm != fixed ] && X="--mf-shape $arm"; d=$O/res/$arm/${R}_$S; mkdir -p $d
  python3 benchmark_hycal_timing_test_20261021.py $E/prad_$R.evio.$S $d/result $WP --lms-corr auto \
    --hc-calib-file $DB/calibration_factor_2p2_fadc.json --gain-file output/percrystal/$R/gains_${R}_all.tsv $X; done
#   per-pulse dumps for diag/intraw (bounded, prior): add --max-events 20000 --pulse-dump $O/dump/$arm/${R}_$S/pulses.npz
python3 mf_fit_pulse.py --self-test
for c in diag res cost intraw; do python3 mf_fit_study.py $c --root $O; done   # intraw writes to <root>/../int_vs_rawsum

# (d) pile-up examples, one run (defaults: --n-examples 10 --dominant-frac 0.5 --secondary-frac 0.10 --pre-window 100 200)
MPLBACKEND=Agg python3 plot_pileup_examples.py --run 25320
MPLBACKEND=Agg python3 plot_pileup_examples.py --run 24386 --splits 10 20 --max-events 20000

# (e) centred rectangle: gains refit (per run), scan with cluster TSVs, module maps
python3 percrystal_calib.py fit --dumps "output/scan_20261023/mf_percrystal/scan_025320/split_*/result_dump.npz" \
  --splits all --fixed-windows --out output/percrystal_rect/025320/
python3 benchmark_hycal_timing_scan.py --pulse-method mf --window 16 --pre-window 100 200 --gain-file-map auto \
  --gain-dir output/percrystal_rect --write-cluster-tsv --parallel 6 --timeout-seconds 1800 --output-dir output/scan_<tag>
MPLBACKEND=Agg python3 plot_module_resolution.py --scan-dir output/scan_<tag> --out output/module_resolution_<tag> --group-runs
#   --min-n 1000 default; --path gated; --exclude-rowcol 16 21 16 21 for scans before 20261024
```

The harness commands for the S_dev and S_full variant scans of §3 and §5 (`gate_scan_20261022.py run --stage R|RX --configs-file … --tags … --runs … --splits …`, then `restable --configs-file … --baseline-tag mf_traw_w16_p100-200`) are in `../plans/plan_20261023.md` §13; the configs file used is `/private/tmp/claude-505/-Users-L00338853-work-PRad-prad2evviewer-analysis-pyscripts/b1de5a72-6ab3-4e2b-8b9b-77e955030d7b/scratchpad/p23/cfgs.json` (scratch; copy it before the session area is cleared).
