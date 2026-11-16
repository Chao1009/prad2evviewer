# Plan: further HyCal energy-resolution improvements (measure first, then implement) — 2026-11-12

Working dir `analysis/pyscripts` (branch `waveform_jz`); paths are relative to it. Evidence base (not repeated here):
`notes/energy_resolution_improvement_options_20261112.md` ("RA"; its option numbers 1-11 are kept). BM =
`waveform_analysis/energy_resolution_benchmark.py`; OLD = `benchmark_hycal_timing_test_20261021.py`; ER =
`waveform_analysis/energy_resolution`. Environment: `source setup_env.sh` (puts `../../build/python`, which holds
`prad2py`, on `PYTHONPATH`), `MPLBACKEND=Agg`.

## Open Questions

None. The user answered RA questions Q1-Q7 on 11-12 (decisions U1-U7). Two lead choices can be changed at approval
without changing the plan: when the 2.24 GeV scans are re-baselined (D5) and the elastic cut behind the efficiency and
leak numbers (D1).

## Goal

Improve the power to separate ep elastic events from inelastic and Møller events (U7), by finding out which resolution
and tail terms are real and then fixing those. One metric protocol (D1-D3) judges everything. Phase 1 (P0-P5, P12,
P13) only measures and attributes; phase 2 (P6-P11, P14) implements, and a package starts only after its gate in D7
passes at the P5 review.

- Critical case: 3.49 GeV, where the pion-production threshold lies only 2.6 σ_E below the elastic peak (D1).
- Pool (r < 200 mm): close or explain the 2.24 GeV gap to the old chain (1.876 vs 1.794 %), and test every remaining
  lever (each ≤ 0.08 pp in the core, RA §0) for its effect on core and tails.
- Coverage: first numbers for the PbWO4 rings outside the pool, the hole ring and the outermost PbWO4 layers, where
  about 870 modules keep a prior gain. Lead glass is outside the acceptance (U1).

Non-goals: redoing the research; a lead-glass resolution or non-linearity (U1); temperature or HV corrections (U5);
production C++ or `../../database/` changes (except P10, a port to `mf_pulse.py` by definition); template-bank work
(in progress outside this plan); anything in Scope Notes "dropped".

Expected change size. Phase 1: about 300 lines added to BM (`--set`, 7 switches, 6 constants, two new output families;
existing outputs unchanged), one compare script (≤ 80 lines), at most four tool scripts of ≤ 80 lines each, about 45
report runs (5-15 min each), 4 single-run rescans (about 51 min each), one five-run rescan (about 4.5 CPU-h), five
result tables, one run inventory and one budget note. Phase 2: sized by the child plans.

## Approach

- One instrumented benchmark, many variants. BM has no flags (every knob is a module constant), so P0 adds a `--set
  NAME=VALUE` override plus the few switches the tests need. A variant is then one command and one output tree, and
  nothing is edited between runs. Rejected: editing constants per variant (not reproducible, files collide across
  packages); a config file (more code than the override).
- Report-level first. The scan stores hits only and all clustering, selection, gains and fits happen in `report`
  (BM:1847), so most tests need no rescan. Scan-level tests (pulse acceptance, phase correction) use paired rescans of
  run 25320 only, extended to more runs only when the 8-split result warrants it (D5).
- Separation metrics come from the same per-split spectra as the widths, so they cost no extra runs.
- The old chain is not rebuilt; its reference is re-established inside BM (D6). Rejected: regenerating R5nr-ri at 2.24
  GeV (bank rebuild plus about 14 CPU-h, and still not the same bank).
- Packages, not a flat plan: the phase-1 lines are independent once P0 lands and run in parallel, and each phase-2
  line needs its own plan stage (design forks in P7 and P9). P0 owns all BM edits so package file sets stay disjoint.

## Decisions

User answers (2026-11-12, binding):

- U1 (Q1) "Lead glass not in acceptance". A lead-glass resolution or non-linearity is not a goal (RA option 6
  dropped). Kept: a diagnostic of whether clusters seeded in the outermost PbWO4 layers share energy with lead-glass
  modules (factor 0 today), in P2; if they do, P11 finds the minimum lead-glass gain for those clusters or a fiducial
  cut. Regions are PbWO4 rings only (D2).
- U2 (Q2) "Gain refinement can be back for comparison, but I worry about over fitting". Refinement returns as a
  comparison arm only; refinement off (10-29 contract) stays the default. Every gain-solve package (P6, P7, P11) is
  judged out of sample, reports the gap (out-of-sample minus in-sample `res_cb`) as its overfitting measure, states
  the number of fitted parameters against events per module, and is killed when there is no out-of-sample gain or when
  the gap grows by more than max(2 SE, 0.02 pp) over the control's gap (today 0.031 / 0.013 / 0.029 pp at 0.73 / 2.24
  / 3.49 GeV).
  - 2026-11-13: P6 shows the old chain's seed-indexed refinement is a training artefact (0 held out); refinement stays comparison-only, default OFF. P1 residual: 0.029 pp attributed, ~0.05 pp unattributed (old bank/pulse chain); D6 rebuild not run.
  - 2026-11-14: gap clause waived for ring-scoped (region-only) calibrations whose pool Δ is 0 by construction; out-of-sample and leave-one-run-out gains are the test.
- U3 (Q3) "We can do a comparison whether they help". Møller rejection and non-linearity each become an on/off
  comparison in P1 (same data, same splits, one switch at a time). Both stay off (BM:153,155) until the user has seen
  that table.
  - 2026-11-13 (user, via P5 table): Møller rejection and non-linearity stay OFF; both null or worse on H2.
- U4 (Q4) "Snack-scan and Cosmic should exist but I don't know which run, how will they help". P13 locates the runs
  and states what each would give. Snake scan: an absolute gain per module at a known beam energy, independent of ep
  kinematics; coverage of the outer PbWO4 modules that ep and Møller populate poorly; a per-module non-linearity if
  several energies exist; the measured position dependence across a module; and one more energy point for the
  three-term fit, which has 0 degrees of freedom today. Cosmics: a beam-independent relative gain equalisation of all
  modules at the MIP scale and a stability check, but weak at the GeV scale because of the non-linearity
  extrapolation. Nothing in P0-P12 depends on them; their use is P14, blocked until the runs are found.
- U5 (Q5) "They are not precise at all" (temperature and HV logs). No temperature or HV correction. Time dependence is
  tracked from the data: the elastic peak per run and per split (P2), per-run gains (already in BM), LMS.
- U6 (Q6) "Good question I don't know. What is the correlation?" Why the readout threshold matters: energy in cluster
  members below threshold is lost; the mean loss makes the response non-linear and position dependent; its
  event-to-event fluctuation adds width, with the largest relative effect at 0.73 GeV; and a threshold in ADC counts
  with unequal gains is a different MeV threshold in every module. P12 measures the threshold from the data and the
  DAQ configuration and sizes the loss. A KTeV-style below-threshold correction follows only if the loss fluctuation
  is a visible budget term.
- U7 (Q7) "Eventually we want to distinguish e-p elastic from inelastic and Moller, think about it". The figure of
  merit is separation power, not the core width alone: D1 adds the separation, the two tails and an efficiency / leak
  pair; D3 decides on them; D13 re-ranks the packages.

Lead decisions:

- D1 (estimators). Per run type, path and region, in sample and out of sample (`*_oos`: even/odd cross fit per run,
  BM:1867-1877). Paths: `legacy` (headline) and the adopted gate (`g12` for H2, `g8` for Ta), so every report runs
  with gates `8 12`.
  - Core (exist, `fit_path`, BM:1214): `res_cb` (two-peak Crystal Ball core σ/E'), `res_g2` (Gaussian iterated within
    ±2σ), `mu_ep`, n.
  - Separation (new): σ_E = `res_cb`·E' in MeV and `sep_pi` = Δ_π/σ_E, with Δ_W = (W² − M²)/(2M)·E'/E the distance
    from the elastic peak to invariant mass W at the same angle. Δ_π = 145 MeV (π⁰p threshold; 151 MeV for π⁺n),
    within 1 % of this value inside the pool. The Δ(1232) lies 340 MeV below the peak, so its separation is 2.35
    `sep_pi` and is quoted here only.

    | Beam energy | σ_E today (legacy, CB) | `sep_pi` | Δ(1232) separation |
    |---|---|---|---|
    | 0.73 GeV | 24.7 MeV | 5.9 | 13.8 |
    | 2.24 GeV | 42.0 MeV | 3.4 | 8.1 |
    | 3.49 GeV | 55.2 MeV | 2.6 | 6.2 |

    For Ta the same offset is reported for comparison only (the inelastic physics differs).
  - Tails (new; today's `tail`, BM:1236, is two-sided): `tail_lo3` and `tail_hi3`, the % of in-range clusters below μ
    − 3σ_CB and above μ + 3σ_CB; `tail_lo_f` and `tail_hi_f`, the % beyond the fixed offset Δ_π from μ. The high side
    is detector only on H2 (pile-up, secondary pulses, miscalibration) and is the proxy for inelastic events smeared
    up into the elastic window. The low side mixes the radiative tail (and Møller events at 0.73 and 2.24 GeV) with
    detector losses and sets the elastic-cut efficiency; only its change between variant and control is interpreted.
  - End to end (new), at a fixed elastic cut E > μ − c with c = 2 σ_E of the 11-12 baseline, frozen in MeV: 49 / 84 /
    110 MeV at 0.73 / 2.24 / 3.49 GeV (constant `SEP_CUT_MEV`; another cut changes the values, not the plan).
    `eff_cut`: the % of in-range clusters that pass. `leak_cut`: the % with E > μ + (Δ_π − c), i.e. the chance that an
    event at the pion threshold passes the cut if the response has the measured high-side shape; Δ_π − c is 95 / 61 /
    34 MeV. `leak_cut` ranks changes; it is not a contamination, which would need an inelastic and radiative generator
    (outside this plan).
  - Møller (new). `kind = moller2` rows: μ and width of (E₁ + E₂)/E_beam for two-cluster events; this is the
    resolution at Møller energies and it does not depend on the beam offset. Both clusters must pass `SIGNAL_MIN`·kE
    (BM:110), which puts the pair members at r = 126-437 / 72-250 / 58-200 mm at 0.73 / 2.24 / 3.49 GeV. `kind =
    moller` rows: μ of E/E_ee(θ) for kin-tagged single clusters per region, against the ep peak of the same region
    (scale linearity). They exist only inside the same outer radii and only for events whose partner was not
    reconstructed above the threshold (BM:1012). Their width is dominated by the angle (dE/E = 2 dθ/θ); it is quoted
    as an upper bound only if `moller2` turns out not to be usable.
- D2 (regions, PbWO4 only, by seed-module centre radius r and position; D8 gives the implementation). H: the rows
  15-20 × columns 15-20 block (excluded from the pool; holds the r < 85 mm hole ring); I1 r < 120, I2 120-160, I3
  160-200 mm (inner = pool = I1+I2+I3, the continuity number); M 200-300 (mid); T: seeds in the outermost two PbWO4
  layers, split into clusters with and without a lead-glass coincidence (U1): a stored lead-glass hit in the same
  event whose module centre lies within 60 mm of the seed centre (the nearest lead-glass centre is 30-35 mm from the
  outermost layer and 50-54 mm from the second; stored hits lie in 100-200 ns, BM:76, and no further time cut is
  applied, because the lead-glass time sits about 14 ns after the seed time); O: r ≥ 300 mm and not T (outer). Pool
  rows are also split edge/centre by impact position (|u| or |v| > 0.35, `uv`, BM:1345). Outside the pool the spectrum
  is E·⟨E'⟩/E'(θ) pooled over the region's seeds; every region row is the mean ± SE over per-split fits, as for the
  pool, and a split with fewer than 50 clusters is NaN (`cb_single`, BM:1117). A region with fewer than 2000 in-range
  clusters, pooled over the run type, or with no finite split is reported as "ep-starved". T itself is a fitted region
  row; its two sub-rows are not fitted (a hit-level count gives 41-81 coincidence clusters per split): they carry n,
  the share and the median E/E', pooled over the run type.
- D3 (comparison and decision rule). Δ = variant − control, paired per (run, split); SE = std(ddof 1)/√n (`mean_se`,
  BM:1684). Anything fitted on data (gains, maps, refinement, NL, tuned cuts) is judged on the out-of-sample Δ.
  Figures of merit: (1) `res_cb` (it sets `sep_pi`), (2) `tail_hi3` and `leak_cut`. Guards: `eff_cut`, `tail_lo3`, and
  the yield loss 1 − Σ `yield_n`(variant) / Σ `yield_n`(control) from `splits.tsv` (`yield_pct` is normalised to the
  same report's legacy path and cannot show a loss there). Real: |Δ| ≥ 2 SE in a figure of merit. Adopt: Δ ≤ −2 SE in
  (1) or (2) at 3.49 GeV with the same sign in `3p49GeV_H2_old`, or at ≥ 2 of the 3 H2 energies (or in the package's
  target region); the other figure of merit not worse by > 2 SE at any energy; `eff_cut` not down and `tail_lo3` not
  up by > 2 SE; yield loss ≤ 1 %. Null (kill): no figure of merit moves by 2 SE (or by 0.01 pp for `res_cb`) at any H2
  energy. With about 35 comparisons a few false "real" results are expected at 2 SE; adoption therefore needs the
  replication above. Gain-type changes also obey U2.
  - 2026-11-13 amendment (user): for variants that narrow σ, the low-tail guard uses the fixed-offset tail `tail_lo_f` instead of `tail_lo3` (tail_lo3 moves with the 3σ line). Log weight R7ri adopted at 2.24 and 3.49 GeV on this basis.
  - 2026-11-13 (supervisor, same rationale as the user's low-side amendment, flagged for the user): for σ-narrowing variants the high-tail guard uses the fixed-offset `tail_hi_f`, not `tail_hi3`. Applied to P9 shape_all2.
  - 2026-11-14: 2.24 GeV tail_lo_f guard waived for P16 (user).
- D4 (datasets). Scans `output/energy_resolution/scan/<RT>/run_<run>.npz`; one split = one EVIO file.

  | Run type | Runs | Splits | legacy `res_cb` ± SE (11-12) |
  |---|---|---|---|
  | `0p73GeV_H2` | 24590, 24603 | 16 | 3.389 ± 0.014 |
  | `2p24GeV_H2` | 25228, 25257, 25260, 25299, 25320 | 42 | 1.876 ± 0.005 |
  | `3p49GeV_H2` (headline) | 24698, 24911, 24937 | 20, of which 16 finite (run 24937: `n_used` 0) | 1.585 ± 0.006 |
  | `3p49GeV_H2_old` (replication) | 24386, 24503, 24512 | 28 | 1.559 ± 0.006 |
  | `2p24GeV_Ta_X17` (also per run in `runs.tsv`) | 25800, 26138, 26230 | 25 | 2.308 ± 0.035 |

  Raw EVIO: `/Users/L00338853/work/PRad/data/evio/prad_0<run>.evio.<NNNNN>` (all present). A Ta variant needs the
  `2p24GeV_H2` variant in the same `--out` root (borrowed position tables, BM:1604-1611).
  - 2026-11-13 amendment (user): run 24937 dropped from the 3.49 GeV H2 set (~24 k events, 389 pool clusters, gains fit not converged); 3.49 GeV = 24698 + 24911.
- D5 (baseline, scan-level tests, re-baseline). BM was edited after the on-disk reports ran, the `2p24GeV_H2` template
  bank was rebuilt on 11-12 (sha1 7f125b75; the scans used 15dfd867, no longer on disk), and `waveform_analysis/` is
  untracked. So P0 first saves a copy and sha1 of BM, reruns the default report for all five run types into
  `output/energy_resolution/baseline_20261112/pre/` and diffs it against `ER/<RT>/summary.tsv`; that tree is the
  control of every report-level variant.
  - Scan-level variants are compared only with control rescan C0 (run 25320, current code and bank, made by P0), never
    with the existing npz, and control and variant must carry the same `bank_sha1`. Both sides are single-run reports,
    because the position correction is fitted per run type; the C0 report is
    `output/energy_resolution/baseline_20261112/c0/2p24GeV_H2`. A variant with |Δ| ≥ 2 paired SE in a figure of merit
    is extended, control and variant, to the other 2.24 GeV runs and the 3.49 GeV runs (about 9.5 CPU-h for both arms)
    before any verdict.
  - Re-baseline: later, at the P5 gate, by rescanning the five `2p24GeV_H2` runs (about 4.5 CPU-h) into
    `output/energy_resolution/scan_variants/rebaseline/2p24GeV_H2/` and rerunning the H2 and Ta reports.
    `output/energy_resolution/scan/` is never overwritten (its bank is gone). Reason: report-level Δ are paired within
    one scan, the 3.49 GeV scans are not affected, and the bank work is still running. Rejected: now, in P0 (same
    cost; the headline 2.24 GeV number would move before P1 has reconciled it with the old chain).
- D6 (P1 reference). The target 1.794 ± 0.004 % (n 42) is note-only (`notes/PbWO4_energy_resolution_20261105.md:63`);
  `output/r5ri_20261104`, `output/r5_20261028` and the best-50 banks are gone. The reference is re-established as arm
  `ref` in BM on the same 42 splits: Møller rejection with Møller-rejected gains and a single CB when rejection is on
  (OLD:3541). The old chain's fixed-window gain solve is not rerun: the old chain measured it against its default at
  +0.004 pp (`percrystal_calib.py:35-39`), and BM has no equivalent setting. P1 closes when `ref` in-sample `res_cb` ≤
  1.81 % (1.794 + 2 SE of both numbers) with every step attributed, or when it reports the unexplained residual.
  Fallback, not executed without the user: rebuild the old chain for run 25320 only (about 3 CPU-h plus a best-50
  bank; anchor `waveform_analysis/opens.md:38`).
- D7 (phase gate). P5 lists, per phase-2 package, the gate metric, the measured value and go / no-go; the user
  approves that table once. Gates:

  | Package | Starts when | Killed when |
  |---|---|---|
  | P6 gain solve (comparison arm, U2) | P1 closed (solver steps attributed) | out-of-sample Δ > −2 SE at every H2 energy and on Ta, or the U2 gap rule fails |
  | P7 beyond the pool | P2: a region outside the pool has ≥ 2000 ep clusters and σ_CB − pool σ_CB ≥ 0.15 pp, \|μ_region − μ_pool\| ≥ 0.5 %, or `tail_hi3` above the pool's by 2 SE | region figures of merit not improved by 2 SE out of sample, or the U2 gap rule fails |
  | P8 secondary acceptance | P3: a (dt, A2/A1) cell carries a width or tail excess ≥ 2 SE and a rescan variant confirms it | H2 gain < 2 SE in both figures of merit, or Ta gated `res_cb` worse by > 2 SE |
  | P9 edge correction | P2: edge − centre ≥ 0.15 pp at 2.24 and 3.49 GeV (old basis: 0.27 / 0.30), or an edge tail excess ≥ 2 SE | best candidate out-of-sample Δ > −2 SE |
  | P10 I(φ) port | P4a: \|Δ\| ≥ 2 SE and ≥ 0.02 pp | below that |
  | P11 boundary and lead glass | P2: ≥ 1 % of T clusters have a lead-glass coincidence and their median E/E' lies ≥ 0.5 % below that of the T clusters without one | below that (no action needed) |
  | P14 snake-scan / cosmic data | P13 finds runs that can be read | no runs |

- D8 (region evaluation leaves the training untouched). With `EVAL_REGIONS` on, clusters with any PbWO4 seed are kept
  and tagged; module set, gains and position correction are still trained on the pool only (`gains.tsv` identical,
  test in Slice 1). Outside the pool the corrections act as the code applies them today: a module the solve never
  selected gets the global scale g0 (BM:1329) and `pc_eval` (BM:1474) applies the pool's map. This is RA option 2's
  first test.
- D9 (option coverage). 1 → P1. 2 → P2 (measure), P7 (implement). 3 → P3, P8. 4 → P6 (comparison arm, U2), with P1c,d
  as its phase-1 evidence. 5 → P2 (edge rows), P4b (log weight), P9. 6 lead glass: dropped as a goal (U1); only the
  boundary diagnostic (P2) and P11 remain. 7 Ta: merged; per-run gains exist, the gate is adopted, refinement is P6's
  Ta column, the pile-up guard is in P3/P8. 8 → P4a, P10. 9 → P1g (σ with NL on/off, U3), P2 (Møller against ep), P7
  (NL fixed before Møller enters a gain solve). 10 noise: deferred; ≤ 0.01 pp at ≥ 2.24 GeV, and the test needs a new
  decode because no scan stores baseline samples; revisit if P5 leaves > 0.3 % unattributed at 0.73 GeV. 11 templates:
  deferred; ≤ 0.01-0.03 pp, no first test, and the banks are being rebuilt outside this plan.
- D10 (literature techniques). NPS / PRad-I global solve with split-sample validation → P6 (pool), P7 (outward, ep +
  Møller). CMS combined calibration → P7: ep, Møller and LMS constants compared pairwise, the spreads give each
  method's precision. CMS multifit / NPS multi-pulse fit → P3, P8. KTeV below-threshold correction → P12 (sizing),
  then an amendment if visible. NA48 constant-term budget → P5. Snake-scan and cosmic calibration (PRad-I, NPS) → P13,
  P14.
- D11 (compute). 18 CPUs; a report uses 8 workers (`NPROC`, BM:132): at most two reports at once, or one report plus
  rescans.
- D12 (results and provenance). Each phase-1 package writes one table `ER/improvement/P<n>_<name>.md` (+ `.tsv`), one
  row per variant from its `delta.tsv` (Outline), with the `bank_sha1` of the scans (npz key, BM:905) and the override
  string. Variant trees: `output/energy_resolution/variants/P<n>/<variant>/<RT>/`; rescans:
  `output/energy_resolution/scan_variants/<variant>/<RT>/`.
- D13 (ranking by separation, U7). Tail and coverage lines move ahead of pure core-width lines. Phase 1 runs in
  parallel; if compute is short the order is P3, P2, P1, P12, P4. Phase-2 child plans are written in the order of the
  measured effect on `leak_cut`, `tail_hi3` and `sep_pi` at 3.49 GeV in the P5 table; the expected order is P8
  (secondary pulses), P7 (modules on prior gains), P9 (edges), then P6 and P10 (core only), then P11.
  - 2026-11-13 (user): phase-2 control = best calibrated gate per run type (TIME_OFFSETS on; 20/16/16/16/8 ns per P5, 12/16/4/12/12 by the 1-SE rule — equivalent within errors) + log weight R7ri at 2.24/3.49; phase-2 order: control change → P8 → P9 ‖ P17 ‖ P6 → P7 → P16 → P10 → P11.
  - 2026-11-13 (supervisor, under the user's 'everything' instruction): phase-2 control v2 = v1 + Ta gate 4 ns (P17) + PC_EXTRA=shape_all2 at 2.24/3.49 GeV (P9; off at 0.73 pending gain re-anchoring). P16 must refit the shape factor after changing gains.
  - 2026-11-14: control v3 defined (v2 + P16 + M1 + P7 rings, benchmark defaults); legacy reproduction string = the four v2 terms plus `--set "BELOW_THR_BY_RT={}" --set "PHASE_RESID_BY_RT={}" --set "RING_CALIB_BY_RT={}"`.

## Outline (P0 only; phase-2 code is outlined in the child plans)

`waveform_analysis/energy_resolution_benchmark.py`:
- `main()` (BM:2091): `--set NAME=VALUE` options (one assignment each, repeatable) are taken out of the argument list
  first, for `report` and `scan`; the rest of the command line is unchanged (`report [--out DIR] <scan dir> <gates>`,
  `scan <out.npz> <evio ...>`). Values are Python literals; NAME must be in the allow-list `OVERRIDABLE` = the seven
  switches below plus the existing `MOLLER_REJECT`, `NONLINEAR`, `MIN_SEP_NS`, `R_MODEL`. Overrides are applied before
  `report()` / `scan()` (the workers are forked at BM:1863 and inherit them) and recorded in a `# overrides:` first
  line of `regions_splits.tsv` and in the scan npz key `overrides`. Invariant: with no `--set`, every existing output
  file is byte-identical to the pre-edit code; new columns go only into new files.
- Switches, defaults = current behaviour: `EVAL_REGIONS = False` (D2, D8); `TWO_PEAK_WITH_REJECT = True` (False =
  single CB when Møller rejection is on, pattern OLD:3541); `CLUSTER_SET = {}` (clusterer fields set after
  `CLUSTER_FIELDS` in `make`, BM:931: `log_weight_thres`, `min_module_energy`); `MOLLER_BEAM_OFFSET = (0.0, 0.0)` mm,
  one value per report; `MOLLER_PAIRS = False`; `SECONDARY_DUMP = False`; scan: `PHASE_CORR = True` (False = each
  module's `ptab` row replaced by its phase average, at the lookup BM:856-858).
- Constants: `REGION_EDGES = (120, 160, 200, 300)`, `TRANSITION_LAYERS = 2`, `REGION_MIN_N = 2000`, `EDGE_UV = 0.35`
  (D2); `SEP_GAP_MEV = 144.7` and `SEP_CUT_MEV` per run type (49 / 84 / 110; Ta as 2.24 GeV) (D1).
- Regions and separation, one output family per run type, always written: `regions.tsv` (columns: path, kind = `ep` |
  `moller` | `moller2`, region, n, n_modules, prior_frac_mean (a module without a `gains.tsv` row counts as 1),
  res_cb, res_cb_se, res_g2, mu, mu_se, med (median E/E'; with n, the only filled columns of the two T sub-rows), σ_E,
  `sep_pi`, `tail_lo3`, `tail_hi3`, `tail_lo_f`, `tail_hi_f`, `eff_cut`, `leak_cut`, each fitted column also `_oos`)
  and `regions_splits.tsv` (the same per run and split). By default they hold the `pool` row only, computed from the
  arrays `fit_path` already holds, and it must equal `summary.tsv` / `splits.tsv` (`tail_lo3` + `tail_hi3` = `tail`).
  Offsets scale with E'/E. With `EVAL_REGIONS`, `cluster_events` (BM:989) keeps and tags `pool["region"]`, flags the T
  clusters with a lead-glass coincidence (D2; the scans keep lead-glass hits with their time and energy 0), and the
  other region rows are added; `pool_modules` (BM:982) still defines the training pool. `kind = moller` rows fit
  E/E_ee(θ) of kin-tagged clusters (`moller_tag`, BM:1531) with θ shifted by `MOLLER_BEAM_OFFSET`.
- `MOLLER_PAIRS`: `cluster_events` keeps only events with exactly one cluster above `SIGNAL_MIN`·kE (BM:1012); with
  the switch on, events with exactly two PbWO4 clusters above it that are back to back in azimuth within 10° are kept
  in a separate list and give the `moller2` row (region `all`; both clusters get the default-arm gains and position
  map).
- `SECONDARY_DUMP`: `secondaries.npz` per run type, one row per legacy-pool cluster in the used modules: `run`,
  `split`, `e_over_ep`, `n_extra`, and for its largest in-window extra pulse `dt` (ns, extra − main in that module),
  `ar` (A_extra / A_main), `member_rank` (0 = seed), `e_extra` (MeV). Source: `h_t`, `h_a`, `h_extra`, `h_k` of the
  scan npz (`h_a` is NaN only for WaveAnalyzer hits, `h_k = 0`).

`waveform_analysis/compare_reports.py` (new, ≤ 80 lines): `compare_reports.py <control dir> <variant dir> [--paths
legacy g12]` joins the two `regions_splits.tsv` on (run, split, path, kind, region), takes `yield_n` from
`splits.tsv`, and writes `delta.tsv` in the variant dir: path, kind, region, n, Δ ± SE for `res_cb`, `res_cb_oos`,
`res_g2`, `tail_hi3`, `tail_lo3`, `leak_cut`, `eff_cut`, `mu`, the yield loss, the in/oos gap (U2), and columns
`verdict_<figure of merit>` = `real` / `null` for that run type (D3; adoption is decided in the package table). No
work at import (`poscorr_comparison.py` runs at import; do not copy that pattern).

Gotchas: the report's default arm is `all` (BM:160) and equals whichever of `MOLLER_REJECT` and `NONLINEAR` are on, so
every switch is judged as variant report against `pre/`, and the per-arm columns of `splits.tsv` (three quantities
only) are not used, except `mu_ep_none` in P12; lead-glass hits in the existing scans carry `h_e = 0` and `h_a = NaN`
(BM:836-840 with factor 0), so only their presence and time can be read; the pulse acceptance (`NSIGMA, MIN_AMP = 5.0,
5.0`, BM:86) is a second threshold on top of the readout threshold (P12).

## Slices

1. P0: tooling and baseline.
   Units: (a) baseline runs and C0 rescan (before and independent of the edits), (b) BM edits, (c)
   `compare_reports.py`.
   Tests: `baseline_matches_disk: default report of each run type -> res_cb, res_cb_oos within 1e-3 of
   ER/<RT>/summary.tsv, else the difference is reported and pre/ becomes the baseline`; `deterministic: 3p49GeV_H2
   default report twice -> splits.tsv byte-identical`; `default_regression: 3p49GeV_H2 report after the edits, no
   --set -> splits.tsv, summary.tsv, gains.tsv byte-identical to pre/`; `override_roundtrip: --set MOLLER_REJECT=True
   on 2p24GeV_H2 -> res_cb differs from pre/, regions_splits.tsv starts with the override line; --set FOO=1 -> exit
   with an error`; `regions_identity: default report -> pool row of regions_splits.tsv equals splits.tsv in res_cb and
   res_cb_oos, tail_lo3 + tail_hi3 = tail, sep_pi = 144.7 (E'/E) / σ_E; --set EVAL_REGIONS=True -> the same pool row,
   n(I1)+n(I2)+n(I3) = pool n_used, gains.tsv identical to pre/`; `pairs: --set MOLLER_PAIRS=True on 2p24GeV_H2 -> a
   moller2 row with n > 0 and a finite μ; μ is reported as a measurement`; `secondary_fraction: --set SECONDARY_DUMP=True on 2p24GeV_H2 ->
   share of rows with n_extra > 0 is within 0.5 pp of the pu column (24.3 %)`; `compare_self: pre/ vs pre/ -> every Δ
   = 0`; `c0: rescan of run 25320 -> report into baseline_20261112/c0/2p24GeV_H2 with res_cb finite on its 8 splits`.
   Checkpoint: `python3 waveform_analysis/energy_resolution_benchmark.py report --out
   output/energy_resolution/baseline_20261112/post --set EVAL_REGIONS=True output/energy_resolution/scan/3p49GeV_H2 8
   12 && python3 waveform_analysis/compare_reports.py output/energy_resolution/baseline_20261112/pre/3p49GeV_H2
   output/energy_resolution/baseline_20261112/post/3p49GeV_H2` prints Δ = 0 on every pool row.
2. P1-P4, P12 and P13 in parallel (run-only packages; variants listed in the package blocks).
   Tests: `P<n>_complete: every variant of the block -> res_cb and res_cb_oos finite for every split that is finite in
   pre/, one delta.tsv per variant`; `P<n>_table: ER/improvement/P<n>_<name>.md -> one row per variant with its
   delta.tsv values, provenance (D12) and a verdict`.
   Checkpoint: `ls waveform_analysis/energy_resolution/improvement/P{1,2,3,4,12}_*.md
   notes/calibration_runs_snake_cosmic_20261112.md` lists six files, and `grep -c verdict` on each `delta.tsv` header
   is 1.
3. P5: budget note, separation table and gate table; the user approves the gates and decides the U3 defaults.
   Tests: `budget_rows: every [I] row of RA §2 -> replaced by a measured value with its source table, or kept and
   marked unmeasured`; `gate_table: one row per phase-2 package -> metric, value, go / no-go per D7, ordered per D13`.
   Checkpoint: `grep -c '^| P' notes/energy_resolution_budget_20261112.md` ≥ 7.
4. Phase 2: one child plan per approved package (P6-P11, P14), then execute. Each child plan carries its own slices
   and uses D1-D3 and U2 unchanged.
   Tests: `P<n>_oos: the package's adopted setting vs the P5 baseline -> meets the D3 adoption rule, or the package is
   closed as null`.
   Checkpoint: the child plan's verification command.

## Verification

`python3 waveform_analysis/compare_reports.py <P5 baseline> <final tree>` for the five run types, the final tree run
with `EVAL_REGIONS=True`: every adopted change meets D3 out of sample; the `default_regression` test of Slice 1 still
passes with all new switches at their defaults; `notes/energy_resolution_budget_20261112.md` lists each term as
measured, null or deferred, and `sep_pi`, `tail_hi3`, `eff_cut` and `leak_cut` per energy before and after.

## Risks

- P1 may end with an unexplained residual, since the reference cannot be rerun (D6). Phase 2 then proceeds on the new
  chain's own basis.
- The split-level SE understates the run-to-run part for fitted quantities (gains and maps are fitted per run or run
  type); the replication in D3 is the guard.
- Ta is statistics-starved (SE 0.035 pp): Ta verdicts use the gated path and the per-run rows, and never decide alone.
- 0.73 GeV CB and Gaussian widths disagree (3.39 vs 3.88 %): both are reported; 0.73 GeV alone never triggers
  adoption.
- Møller: μ of E/E_ee moves 0.3-0.9 % per mm of beam offset, and `output/nonlin_20261104/step3/beam_offsets.tsv`
  covers only the `0p73GeV_H2`, `2p24GeV_H2` and `3p49GeV_H2_old` runs (spread within a run type ≤ 0.13 mm); the other
  two run types get `moller` rows marked "offset unknown". The old-chain study found few pair-tagged events at r ≥ 85
  mm (`output/nonlin_20261104/step3/findings_step3.md`), so `moller2` is a first test and may end as "not usable".

## Scope Notes

In: P0-P5, P12, P13 as written; child plans and execution of the phase-2 packages that pass D7. Phase-1 packages do
not edit BM; a missing switch goes back to P0 as an `Amend parent:`. Out by user decision: lead-glass resolution and
non-linearity (U1); temperature and HV corrections (U5); gain refinement as a default (U2;
`notes/PbWO4_energy_resolution_20261029.md:24`, `plans/plan_r5_20261028.md:50`); a change of the Møller-rejection or
non-linearity default before the P1 table (U3); the runinfo calibration basis stays (`plans/plan_ri_20261104.md`); no
commits; concise code. Out, dropped by the research (RA §4.3), not to be reopened: crosstalk correction (< 0.03 pp;
rms 2-3 MeV per cluster), saturation (0.0003-0.014 % of clusters), within-run LMS correction (measured +0.000 pp),
tighter time gates and χ² cuts (lose yield or energy) (reopened 2026-11-13 by user decision as P15: with per-module time offsets applied; the raw-time scan is superseded), leakage, Geant4 profile, `min_module_energy` as a σ lever
(inert; it stays in P12 as a threshold probe), integral vs template fit (decided for the template fit). Out of scope
at this stage: running anything, code changes, commits.

## Assumptions

- `report` is deterministic: checked by `deterministic` in Slice 1.
- The phase-averaged table reproduces production `mf_pulse` (`notes/ep_resolution_walkthrough_20261102.md:63`):
  checked in P4a by a mean-energy change < 0.05 %.
- The beam offsets of the old-chain study apply to the new clusters: checked in P2 by a flat E/E_ee against azimuth.
- The high side of the H2 elastic peak holds no physics: checked in P3 by its dependence on the pile-up flag and the
  time gate.
- The raw data hold every channel of an event, so a readout threshold can be read off: checked first in P12.

## Work Packages

Phase 1: measure and attribute.

### P0: protocol tooling and baseline
- Goal: BM takes overrides and writes region, separation, pair and secondary outputs with existing outputs unchanged,
  and a compare script gives the paired Δ; done when Slice 1 passes
- Owns: `waveform_analysis/energy_resolution_benchmark.py`, `waveform_analysis/compare_reports.py`,
  `output/energy_resolution/baseline_20261112/`, `output/energy_resolution/scan_variants/c0/`
- Deps: none
- Stages: execute
- Status: done
- Artifacts: -
- Cost: about 1.5 days of code; 45 min of reports twice; one rescan (51 min). Expected gain: none (enables D1-D3).

### P1: reconcile the 2.24 GeV benchmark with the old chain; Møller rejection and NL comparisons (option 1, U3)
- Goal: attribute 1.876 vs 1.794 % step by step, and give the on/off table for Møller rejection and non-linearity;
  done when `ref` ≤ 1.81 % in sample or the residual is named (D6), row a exists for every run type and row g wherever
  the calibration has nl ≠ 0
- Owns: `output/energy_resolution/variants/P1/`, `ER/improvement/P1_reconcile.md`
- Deps: P0
- Stages: execute
- Status: done
- Artifacts: -
- First test (one switch at a time, same splits, each against `pre/`): (a) `MOLLER_REJECT=True` on all five run types;
  (b) a + `TWO_PEAK_WITH_REJECT=False` on `2p24GeV_H2` (42 splits); (c) read, no run: `ER/poscorr_comparison.tsv`
  gives gains only at 1.9241 % in sample / 1.9346 % out of sample, so the gain step is worth 0.080 pp and the position
  step 0.048 pp here against 0.18 and 0.05 in the old chain (`waveform_analysis/opens.md:38`), and it gives the
  position-form rows; (d) read, no run: fixed gain windows, +0.004 pp in the old chain (D6); (e) `ref` = row b, no
  extra run; (f) pulse model `R_MODEL` 0.3, cited from P3 when available and not needed to close P1; (g)
  `NONLINEAR=True` on `0p73GeV_H2`, `3p49GeV_H2`, `3p49GeV_H2_old` (the only calibration files with nl ≠ 0; at 2.24
  GeV the switch is an identity, which the table states).
- Success / kill: D6; rows a and g get a D3 verdict each. Expected gain: up to −0.08 pp at 2.24 GeV. Cost: 9 reports,
  about 2 h.

### P2: region-resolved resolution and separation (options 2, 5, 9 first tests; U1 boundary diagnostic)
- Goal: core, tails and prior fraction per PbWO4 region (D2) for ep at all energies, the Møller rows of D1, and the
  T-layer split by lead-glass coincidence; done when the table has every region finite or marked ep-starved
- Owns: `output/energy_resolution/variants/P2/`, `ER/improvement/P2_regions.md`, `ER/improvement/tools/p2_*.py`
- Deps: P0
- Stages: execute
- Status: done
- Artifacts: -
- First test: `--set EVAL_REGIONS=True --set MOLLER_PAIRS=True` on the five run types, with `MOLLER_BEAM_OFFSET` = the
  mean over the run type's runs in `output/nonlin_20261104/step3/beam_offsets.tsv` where it exists. Read from the
  outputs: edge vs centre rows of the pool; `moller2` μ and width per energy (resolution at Møller energies);
  μ(Møller) − μ(ep) per region against E_ee (scale linearity); T clusters with vs without a lead-glass coincidence:
  share and median E/E' deficit; `mu` per split against split order per run (drift, U5); prior fraction per region
  from `gains.tsv`.
- Success / kill: a measurement; it sets the P7, P9 and P11 gates (D7). Expected gain: unknown, first numbers. Cost: 5
  reports, about 1.5 h.

### P3: secondary pulses in cluster members (option 3; multi-pulse fits of CMS and NPS)
- Goal: decide whether the 24-34 % flagged secondaries are spurious and which tail they feed; done when the (dt,
  A2/A1) table and the paired rescans give a D3 verdict
- Owns: `output/energy_resolution/variants/P3/`, `output/energy_resolution/scan_variants/{minsep12,rmodel0p3}/`,
  `ER/improvement/P3_secondaries.md`, `ER/improvement/tools/p3_*.py`
- Deps: P0
- Stages: execute
- Status: done
- Artifacts: -
- First test: (A) `SECONDARY_DUMP` on `2p24GeV_H2`, `3p49GeV_H2`, `2p24GeV_Ta_X17`: share of clusters, width excess,
  and low-side and high-side tail excess per (dt, A2/A1) cell against clusters without a secondary. (B) rescans of run
  25320 against C0: `MIN_SEP_NS=12` (the open to-do, `notes/PbWO4_energy_resolution_20261108.md:16`) and `R_MODEL=0.3`
  (the old chain's value; BM uses 0.02, BM:89): flag share, d_pu, `res_cb`, `tail_hi3`, `leak_cut`. Extension per D5.
  Ta guard only for a passing variant: rescan 26138 (4 splits), gated path.
- Success / kill: closes as "real small pulses" if no cell carries an excess and both rescans are null; otherwise gate
  for P8. Expected gain: ≤ 0.04-0.05 pp H2 legacy, ≤ 0.03-0.04 gated, Ta ≤ 0.11 gated; tail effect unknown, first
  numbers. Cost: 5 reports, 2 rescans (1 h wall).

### P4: pool attribution switches (options 5, 8)
- Goal: one D3 verdict per switch; done when rows a and b are in the table
- Owns: `output/energy_resolution/variants/P4/`, `output/energy_resolution/scan_variants/phase_off/`,
  `ER/improvement/P4_switches.md`
- Deps: P0
- Stages: execute
- Status: done
- Artifacts: -
- First test: (a) `PHASE_CORR=False` rescan of run 25320 against C0; (b) `CLUSTER_SET` `log_weight_thres` = the R7ri
  value per energy (3.41 / 3.66 / 3.71, `output/r7_20261104/configs/hycal_<RT>.json`) against 5.2, on the three H2
  energies and `3p49GeV_H2_old`.
- Success / kill: D3 per row. Expected gain: a about −0.04 pp (production only), b ≤ 0.03 pp at 3.49 GeV. Cost: 5
  reports, 1 rescan; about 1.5 h.

### P12: readout threshold and below-threshold loss (U6; KTeV)
- Goal: the effective threshold per module in ADC and MeV, and the energy lost below it per cluster (mean and
  event-to-event fluctuation) against beam energy and impact position; done when the table states whether the loss
  fluctuation is a visible budget term (≥ 10 % of σ_E in quadrature at any energy)
- Owns: `ER/improvement/P12_threshold.md`, `ER/improvement/tools/p12_*.py`, `output/energy_resolution/variants/P12/`
- Deps: P0
- Stages: execute
- Status: done
- Artifacts: -
- First test: (1) from one raw split per energy: channels present per event and per channel (are channels without a
  pulse read out), and the lower edge of each channel's pulse-amplitude distribution; (2) from the configuration:
  `../../database/daq_config.json:90` (TET 10 ADC above pedestal, NSAT 4), `../../database/daq_cnf/`, and the fit's
  own acceptance (BM:86); (3) MeV per module through the calibration factors; (4) loss: `CLUSTER_SET`
  `min_module_energy` = 8, 12, 16 MeV at the three H2 energies (the stored hits start at about 6 MeV, so lower values
  change nothing): the mean loss from `mu_ep_none` of `splits.tsv` (the peak before gains; the default arm refits the
  gains and hides it), `res_cb` and `tail_lo3` from the default arm, each against threshold, extrapolated to zero;
  `res_cb` and `tail_lo3` also edge vs centre.
- Success / kill: a measurement. A KTeV-style correction becomes an `Amend parent:` package only if the term is
  visible. Expected gain: unknown; largest at 0.73 GeV. Cost: 9 reports, three short decodes; 2.5 h.

### P13: locate snake-scan and cosmic runs (U4)
- Goal: a run inventory with, per run, type, beam energy, date, where the data are, and what it would give (U4); done
  when the note lists the runs or states "none found" with every place searched
- Owns: `notes/calibration_runs_snake_cosmic_20261112.md`
- Deps: none
- Stages: research
- Status: done
- Artifacts: `notes/calibration_runs_snake_cosmic_20261112.md`
- First test: search `../../calibration/README.md:81` (snake-scan operator manual), `../../calibration/cosmic_gain/`,
  `../../database/gain_equalization_results.json`, `../../database/gain_factor/`, `../../database/runinfo/`,
  `/Users/L00338853/work/PRad/data/` (`evio`, `raw_root`, `physics_calib_out`), and RCDB
  (`https://clasweb.jlab.org/rcdb/runs/info/<run>`) or logbook exports if the user provides access. Expected gain:
  none by itself (unblocks P14). Cost: under 1 h.

### P5: resolution budget, separation table and phase gate (NA48-style)
- Goal: a current budget and separation table per energy from the P1-P4 and P12 tables, and the D7 gate table ordered
  per D13; done when the user has approved the gates and decided the U3 defaults
- Owns: `notes/energy_resolution_budget_20261112.md`, `output/energy_resolution/scan_variants/rebaseline/`
- Deps: P1, P2, P3, P4, P12
- Stages: execute, write
- Status: done
- Artifacts: -
- First test: one budget row per term of RA §2 from the phase-1 tables; the re-baseline of D5 (rescan of the five
  `2p24GeV_H2` runs, H2 and Ta reports), which gives the phase-2 control. Expected gain: none (decides phase 2). Cost:
  about 4.5 CPU-h of rescans, two reports, half a day of writing.

Phase 2: implement (each starts at its D7 gate, in the D13 order; child plans at
`notes/plan_energy_resolution_improvement_20261112_P<n>.md`).

### P8: secondary-pulse acceptance (option 3)
- Goal: an acceptance setting that removes the spurious secondaries without losing Ta pile-up efficiency
- Owns: pulse-acceptance constants of BM (BM:85-97) and the full rescan
- Deps: P3, P5
- Stages: plan, execute
- Status: done
- Artifacts: `notes/plan_energy_resolution_improvement_20261112_P8.md`
- First test: the P3 winner on all runs; guard: the synthetic two-pulse efficiency stays ≥ 97.5 % at dt ≥ 12 ns, r ≥
  0.1 (`waveform_analysis/pileup_analysis/efficiency.tsv`). Expected gain: as P3. Cost: one full rescan (11.8 CPU-h).

### P7: calibration beyond the pool (option 2; PRad-I ep + Møller, CMS combined constants, option 9 NL)
- Goal: measured gains for the modules now on the prior, validated out of sample per region (U2)
- Owns: a new training-region mode of BM and `ER/improvement/P7_*`
- Deps: P2, P5, P6
- Stages: plan, execute
- Status: done
- Artifacts: `notes/plan_energy_resolution_improvement_20261112_P7.md`
- First test (in the child plan): extend the ep training region ring by ring where P2 shows ≥ 2000 clusters; where ep
  is starved, fix the NL from the Møller-vs-ep peaks first, then add Møller pairs; compare ep, Møller and LMS
  constants pairwise. Design fork for the planner: ep-only extension vs ep + Møller.
- Expected gain: unknown (P2 gives the first number). Cost: medium.

### P9: edge-aware position correction (option 5; CMS regression, HPS edge parametrisation)
- Goal: bring the edge width and tails toward the centre values
- Owns: position-correction part of BM (`fit_pc`, BM:1455; `pc_eval`, BM:1474)
- Deps: P2, P4, P5
- Stages: plan, execute
- Status: done
- Artifacts: `notes/plan_energy_resolution_improvement_20261112_P9.md`
- First test (child plan; one planner-worker per option): ring-dependent maps, an edge-distance parametrisation, a
  shower-shape regression (`cluster_regression.py` exists), W0 per energy from P4b. Expected gain: upper bound 0.1 pp
  at 3.49 GeV [I]. Cost: medium.

### P6: gain solve as a comparison arm (option 4, U2; NPS global solve, split-sample validation; Ta part of option 7)
- Goal: measure whether refinement or a changed ridge solve lowers the per-module scale residual out of sample; done
  when each setting has a D3 verdict and its U2 numbers, with refinement off still the default
- Owns: gain-solve part of BM (`solve_gains`, BM:1286; `corrections`, BM:1573) after P0-P5 are done
- Deps: P1, P5
- Stages: plan, execute
- Status: done
- Artifacts: `notes/plan_energy_resolution_improvement_20261112_P6.md`
- First test: slice 1 ports the seed-indexed refinement (`fit_elastic_window_gains`, OLD:823) as a switched arm on all
  five run types; slice 2 scans the prior strength, the free-module set and the gain-solve sample (seed energy
  fraction, a narrower low-side window) of the ridge solve. Reported per setting: out-of-sample Δ, the gap against the
  control's gap, the number of fitted factors (one per free seed module, 267-284 today) and the events per module
  (minimum 30, `MIN_PER_MODULE`, BM:116; median and minimum quoted).
- Success / kill: U2. Expected gain: −0.02 to −0.04 pp ep, −0.07 to −0.09 pp Ta (old basis, in sample). Cost:
  low-medium.

### P10: onset-phase correction I(φ) in production (option 8)
- Goal: `mf_pulse.py` applies the per-module phase table instead of the phase average
- Owns: `mf_pulse.py`
- Deps: P4, P5
- Stages: plan, execute
- Status: done
- Artifacts: `notes/plan_energy_resolution_improvement_20261112_P10.md`
- First test: P4a. Expected gain: about −0.04 pp in production [I]; the benchmark number does not change. Cost: low.
  Note: `mf_pulse.py` carries an uncommitted `MIN_SEP_NS` 12 → 8 edit that the plan must keep.

### P11: boundary clusters and lead glass (U1)
- Goal: outermost-PbWO4 clusters that share energy with lead glass are either recovered with a minimum lead-glass gain
  or removed by a fiducial cut, whichever keeps more acceptance at equal separation
- Owns: lead-glass branch of the scan (`wa_hit`, BM:836) and `ER/improvement/P11_*`
- Deps: P2, P5
- Stages: plan, execute
- Status: done
- Artifacts: `notes/plan_energy_resolution_improvement_20261112_P11.md`
- First test: rescan one run per energy storing the raw lead-glass integral (the existing scans lost it); one gain per
  lead-glass module group from the T-cluster peak deficit against the integral, judged per U2; against it, a fiducial
  cut on the seed layer. Expected gain: unknown (P2 gives the first number). Cost: medium.

### P14: use of snake-scan and cosmic data (U4)
- Goal: use the located runs for what U4 lists (absolute gains outside the pool, a fourth energy point, position
  dependence, gain equalisation)
- Owns: to be set by its plan
- Deps: P13
- Stages: research, plan, execute
- Status: execute
- Child plan: notes/plan_p14_snake_scan_20261114.md
- First test: set by its research stage (read one located run with the existing decoder). Expected gain and cost:
  O/T ring constants and a per-module single-energy resolution; cost set by its plan once the sample is local.

### P15 calibrated-time gate scan
- Scope: apply per-module time offsets (`../../database/hycal_time_offsets/time_calib.json`, calibrated = raw − offset_ns) in the benchmark's report-side `add_hit` as a default-off `--set TIME_OFFSETS` switch, rerun the gate scan 1–16 ns on all five run types, report the best gate per run type vs legacy and vs the raw best gate
- Owns: `waveform_analysis/energy_resolution_benchmark.py` (default-off `--set TIME_OFFSETS` switch only), `waveform_analysis/energy_resolution/improvement/P15_*`, `output/energy_resolution/variants/P15/`
- Deps: P0
- Status: done
- Artifacts: `waveform_analysis/energy_resolution/improvement/P15_tgate_calibrated.md` / `.tsv`, `output/energy_resolution/variants/P15/`
- Source: `notes/timing_consistency_gate_vs_crosstalk_20261113.md`

### P16 below-threshold (KTeV-style) correction
- Scope: slice 1 — measurement on existing data, modules with threshold < 8 MeV vs > 11 MeV (res_cb, mu_none, edge); slice 2 — the correction (add expected sub-threshold energy per cluster), only if slice 1 is real
- Deps: P12, P5
- Status: done
- Acceptance: D3 out of sample at ≥ 2 H2 energies
- Artifacts: `waveform_analysis/energy_resolution/improvement/P16_threshold_corr.md` / `.tsv`, `output/energy_resolution/variants/P16/`
- Source: `ER/improvement/P12_threshold.md`, `ER/improvement/P5_budget.md`

### P17 asymmetric time window
- Scope: tight early-side / loose late-side member window (e.g. +2–3 ns early, ~12 ns late, Ta-led scan) against the phase-2 control
- Deps: P15, P5
- Status: done
- Acceptance: D3 on Ta X17 and ≥ 1 H2 energy, pile-up flag and tail_hi3 must not worsen
- Artifacts: `waveform_analysis/energy_resolution/improvement/P17_asym_window.md` / `.tsv`, `output/energy_resolution/variants/P17/`
- Source: `waveform_analysis/energy_resolution/improvement/P15_tgate_calibrated.md` (figure section), `ER/improvement/P5_budget.md`

## Log

- 2026-11-12 P0: done; `--set`, seven switches, regions/separation/pair/secondary outputs and compare_reports.py in place; baseline_20261112 frozen (pre/ reproduces ER summaries exactly, default outputs byte-identical after the edits), C0 rescan (bank 7f125b75) and report made; `pairs` criterion amended (moller2 μ −2.54 % is an ep-vs-Møller scale effect, not a bug)
- 2026-11-12 P13: done, notes/calibration_runs_snake_cosmic_20261112.md: snake scan 23699-23723 (2108.85 MeV bremsstrahlung, 1.16e9 events) and cosmics 23741-45, 24154-24185, 24387-24422, 24617 found in public RCDB, all at JLab only (none local; runs ≤23849 need history DAQ map); P14 unblocked for research; deferred: data copy, scan log, HV history and tagger readout (user to supply)
- 2026-11-12 decision: user approved the plan ("Approve the plan"); project store is tracked ("I want it tracked")
- 2026-11-12 decision: P0 `pairs` criterion amended on the P0 lead's finding (supervisor accepted; reported to the user)
- 2026-11-12 decision: user schedules P14 (snake-scan and cosmic data) for 2026-11-13 and asked to be reminded
- 2026-11-13 P2: done; region/Møller/T-lead-glass tables for 5 run types in ER/improvement/P2_regions.md (pool continuity exact; D7 input: T lead-glass coincidence 8-30 %, median deficit 7.7-19.6 % vs 0-2 % without, so P11 gate input positive); deferred: E/E_ee flat-vs-azimuth assumption check (no azimuth output in BM, not in P2 block)
- 2026-11-13 P1: done; ref (row b) res_cb 1.882 ± 0.005 % vs target ≤1.81, residual 0.088 pp (gain-solve step −0.080 vs old −0.18, Møller +0.006 vs old −0.010; unattributable without the D6 fallback); Møller rejection worse (+0.167 pp at 0.73, +0.009 at 2.24, null at 3.49, Ta ~−0.018 ± 0.009 with worse tail_hi3), NL worse at 3.49 (+0.095/+0.103 pp) and a non-replicating −0.033 pp at 0.73, both stay OFF; deferred: D6 fallback (old-chain rebuild, ~3 CPU-h), row f (R_MODEL 0.3, P3 supplies)
- 2026-11-13 P12: threshold = TET 10 ADC = 9.5-9.7 MeV median (5-13 MeV across modules); loss term 0.18-0.64 of res_cb on the amplitude ratio (variance share 0.03-0.41), so visible at 0.73 GeV and under the extrapolation at 2.24/3.49 GeV; KTeV-style correction is an Amend-parent candidate; deferred: module-split edge test (outside block), no production change.
- 2026-11-13 P4: done; a (PHASE_CORR off) res_cb +0.063±0.007 pp, leak_cut +0.44 %, P10 gate met on run 25320 only; b (R7ri log_weight_thres) res_cb −0.026/−0.033 pp and leak_cut −0.37/−0.50 % at 3.49/3.49_old (real, same sign), −0.011 pp at 2.24, null at 0.73, but tail_lo3 +0.07/+0.13 (fails the guard), so no adoption; deferred: multi-run extension of row a (D5), tail_lo3 trade-off decision to the P5 gate
- 2026-11-13 P3: (A) cells with width and tail excess exist (secondaries carry ~1 % dropped energy; hi tail at A2/A1>=0.25 or dt<0); (B) 25320 rescans: MIN_SEP_NS=12 null, R_MODEL=0.3 removes 97 % of secondaries and gives res_cb -0.036±0.008 pp (g12 -0.020), leak_cut -0.16±0.03; P8 gate condition met on one run, not yet a D5 verdict; deferred: D5 extension of R_MODEL=0.3 to the other 2.24 GeV and 3.49 GeV runs (~9.5 CPU-h) and Ta guard (26138), both outside this block
- 2026-11-13 decision: user: P14 deferred to the next few days; git/commits done by the user later.
- 2026-11-13 decision: user: "Go with everything we can do to improve the energy resolution" — proceed through P5 into phase 2 without a gate pause; P11 unblocked on P2's positive D7 input; D5/D6 extensions authorised at P5's judgement.
- 2026-11-13 decision: user: per-module time-offset correction must be applied before any time-gate scan; report the best gate per run type, never a fixed 12 ns. New package P15. Research: notes/timing_consistency_gate_vs_crosstalk_20261113.md (gate used raw MF times; crosstalk study used time_calib.json offsets; Δt robust σ 2.3-3.1 ns raw vs 0.7-1.4 ns calibrated).
- 2026-11-13 P15: calibrated-time gate scan (TIME_OFFSETS switch, default off) done for 5 run types at 10 gates; legacy reproduces baseline; calibrated gate beats raw only below 8 ns, H2 best gate equals raw best within 1 SE (best minus legacy -0.013 to -0.021 pp), Ta X17 -0.30 pp vs legacy and -0.03 pp vs raw best; deferred: seed-member Δt recount, time walk. Addendum: H2 extended to 20 and 24 ns, scan minima now inside range (16 ns at 2.24 and 3.49 GeV, 20 ns at 0.73 GeV), rule-best 12/16/12/12 ns unchanged.
- 2026-11-13 executor: P15 seed-member dt figure and script made (improvement/figures/P15_seed_member_dt.png; raw shoulder/cliff at +8–9 ns is the offset spread; late tail 0.06/1.75/3.06 % H2, 4.2 % Ta beyond 8 ns, untouched by offsets), Figure section appended to P15_tgate_calibrated.md; deferred: none
- 2026-11-13 P5: done; latest table (best calibrated gate 20/16/16/16/8 ns: 3.365/1.856/1.565/1.546/2.004 %; + log weight 1.847/1.538/1.509), budget, separation table and gate table in ER/improvement/P5_budget.md; adopt P15 best gate + P4b log weight (2.24, 3.49); P8 go but R_MODEL 0.3 killed by the Ta guard; P9, P6, P7, P10, P11 go; new P16 (threshold) and P17 (asymmetric window); deferred: PHASE_CORR extension (to P10), D6 rebuild (P6 fallback), D5 re-baseline (shown unnecessary)
- 2026-11-13 P5 compute: R_MODEL 0.3 D5 extension, 7 rescans (3.8 CPU-h incl. the 25228 identity rescan) + 10 reports, existing scans as control
- 2026-11-13 decision: user: D3 low-tail guard amended to tail_lo_f for σ-narrowing variants; phase-2 order approved as P5 recommends; run 24937 dropped from D4; U3 defaults stay OFF. P1 closed with the 0.088 pp residual named.
- 2026-11-13 executor P0-phase2: control set (TIME_OFFSETS on, calibrated gate 20/16/16/16/8 ns, log weight R7ri 3.6576 at 2.24 and 3.7073 at 3.49 GeV, 24937 dropped) reproduces P5 in-sample to 0.001 pp (control 3.365 / 1.846 / 1.538 / 1.509 / 2.004 %) and legacy `--set` is bit-identical to baseline `pre/`; 24937 effect null (≤ 0.0005 pp); Ta oos +0.023 vs P5 (borrowed log-weighted 2.24 position correction); deferred: Ta gate 5 vs 8 ns, `poscorr_comparison.py` defaults check
- 2026-11-13 P17: done; asymmetric window (ASYM_WINDOW switch, late to C++, early python pre-filter) buys at most 0.013 pp (Ta, 1.8 SE) and ~6 % of tail_hi3 on H2, null on res_cb; not adopted; Ta symmetric gate 8 -> 4 ns recommended (tail_hi3 -0.25, leak -0.42, pile-up -0.67, res null); deferred: Ta CONTROL_GATE change and control rebuild (supervisor/user decision)
- 2026-11-13 P6: done; seed-indexed refinement (the old chain's second gain layer) is −0.03 pp H2 and −0.08 pp Ta in sample per split but 0 held out (2.24 +0.003 ± 0.001, 3.49 −0.002 ± 0.003, Ta −0.019 ± 0.010, 0.73 −0.012 ± 0.004) and +0.02 worse cross-fitted, so a training artefact; member-solve settings (prior, passes, windows, seed fraction) all null; 0.029 of P1's 0.088 pp residual attributed (in-sample seed refit), about 0.04 to 0.06 pp unattributed; recommend refinement stays comparison-only, default OFF (U2 = user); deferred: D6 old-chain rebuild (not needed for the verdict), 0.73 GeV run-level −0.012 pp replication, 3.49 GeV old runs.
- 2026-11-13 P8: done; no secondary-pulse acceptance passes D3 at the phase-2 control gate (merges C1/C2/S1: res_cb null, tail_hi3 +0.008-0.030 pp replicated, eff_cut +0.3-0.45, Ta 26138 g8 res_cb +0.055/+0.070 kill; drop C3 worse +0.008-0.014 pp; pre-pulse drop C4 null); mechanism: gated clusterer already absorbs in-gate extra pulses; no production change; `SECONDARY_RULE` switch (default {}) added to BM; deferred: seed-time-aware merge (C++ grow_island, outside P8), 0.73 GeV (no candidate)
- 2026-11-13 P9: shape factor shape_all2 found (oos res_cb −0.057/−0.035 pp at 2.24/3.49, leak_cut −0.4 pp, eff_cut +0.3-0.5 pp); edge−centre gap not closed, edge eff_cut deficit closed <10 %, position-only edge forms null; deferred: tail_hi_f ruling for D3, refit after P6/P16, 0.73 GeV yield guard 2.0 % (peak shift)
- 2026-11-13 decision: supervisor under the user's "everything" instruction: Ta gate 4 ns adopted; shape_all2 adopted into control v2 at 2.24/3.49; tail_hi_f guard applied symmetrically to D3; P6 closed without the D6 rebuild; asymmetric window and secondary-pulse rules not adopted. Flagged for the user.
- 2026-11-13 P11: lead-glass drop (b) and held-out presence-proxy corrections (c0/c1/d1) do not improve T sep_pi or tails; the lead-glass coincidence is outermost-layer (17-62 % of layer 0, under 5 % of layer 1) and a seed-layer cut gives T res_cb −0.20 to −0.43 pp and sep_pi +0.4 at 0.73/2.24/3.49 GeV for 43-47 % of T; pool null 0 differences; keep (a); deferred: (c)/(d) with the real lead-glass integral (the scans store none; rescan ~38 min per run), and the layer-0 fiducial decision to the user (U7).
- 2026-11-13 P7: ring-scoped per-ring calibration trained on pool + ring (all runs combined, held-out by even/odd and by run; pool Δ 0.000 by construction) gives H −0.30 / −0.69 / −0.85 pp and M −0.06 / −0.12 / −0.21 pp (2.24 / 3.49 / 3.49 old; M −0.46 at 0.73 GeV, H +0.12 at 0.73 GeV); O and T null; position map and shape factor on regions harmful; U2 gap clause fires for H and M (user's call); deferred: Møller-pair / snake-scan / cosmic constants for O and T (NL fix first, U3; P14), production consumer for ring-scoped constants.
- 2026-11-13 P16: done; slice 1 loss real at module level (absent modules +0.39 to +0.52/MeV, ring energy −4 to −6 MeV/MeV at 2.24/3.49/3.49old) but class split null in res_cb and wrong-signed in mu (slice 2 gated on module-level criterion, flagged); KTeV-style BELOW_THR correction (default off) cross-fitted oos res_cb −0.215/−0.026/−0.019/−0.018 pp, leak_cut −0.17 to −0.26 (2.24/3.49), all guards pass at 0.73/3.49/3.49old, 2.24 tail_lo_f +0.053±0.006 fails formally; recommend adoption at 0.73/3.49/3.49old, 2.24 pending guard ruling; deferred: Ta table (Ta borrows 2.24 position correction), 0.73 GeV position-correction arm, slice-1 confound identification
- 2026-11-13 P10: done; D5 PHASE_CORR=False worsens res_cb +0.055±0.003 pp (2.24, 42 pairs, all 5 runs) and +0.050±0.005 pp (3.49, 2 runs), leak_cut +0.45/+0.78 %, so keep PHASE_CORR=True; residual first harmonic (−4.0, +1.9) permille per seed fraction survives the correction (about half left, module-specific); refined form M1 (global k=1,2 + per-module k=1 table, cross-fitted) res_cb −0.028/−0.040/−0.034 pp at 2.24/3.49/3.49 old, leak_cut −0.2/−0.4 %, guards ok, replicated on other runs, 0.73 −0.019±0.007; switches PHASE_STUDY/PHASE_RESID (default off) added to BM; compute 6 rescans + ~20 reports (≈5 CPU-h); deferred: production table fitter and storage, M1 default decision (user), cross-energy and Ta tests, 0.73 yield +0.6 %, root cause of the 0.75× rotated I(φ)
- 2026-11-14 decision: P16 below-threshold correction adopted at all four H2 run types (0.73/2.24/3.49/3.49 old); the 2.24 GeV tail_lo_f guard (+0.053±0.006 on a 26.9 % tail) is waived by the user; module-level slice-2 criterion accepted; Ta gets its own table in control v3.
- 2026-11-14 decision: P10 refined phase correction M1 adopted with per-energy tables (no joint table); PHASE_CORR stays True.
- 2026-11-14 decision: P7 U2 gap clause waived for ring-scoped calibrations; adopt M ring at 2.24/3.49 old/0.73 and H ring at ≥2.24 GeV; O and T stay on the prior; pool headline unchanged by construction.
- 2026-11-14 decision: P11 outermost-PbWO4-layer fiducial cut not applied; current T-ring treatment kept; documented as an option.
- 2026-11-14 executor: Phase-2 control v3 = v2 + P16 BELOW_THR (4 H2; Ta own table fails tail_lo_f +0.86/eff_cut −0.64, not adopted) + P10 M1 refit on the chain (all 5, Ta incl.) + P7 ring scope (M 0.73/2.24/3.49old, H 2.24/3.49/3.49old; pool byte-identical), cross-fitted, now benchmark defaults (`*_BY_RT`, tables in ER/tables/); oos res_cb 3.160/1.746/1.470/1.443/1.987 % (v2 3.397/1.802/1.526/1.500/2.028), Δ −0.237/−0.056/−0.057/−0.057/−0.041, additive, guards pass (2.24 tail_lo_f +0.056 waived), control gates confirmed; v3 fit 0.89 ⊕ 2.06/√E ⊕ 1.34/E %; v2 and legacy byte-reproduced (legacy string now needs the three `*_BY_RT={}` terms); deferred: 2.24 tail_lo_f mechanism, Ta P16 cause, 0.73 position arm, M-ring sub-additivity cause, O/T constants (P14), production ring consumer
- 2026-11-14 decision: P14 — HyCal HV changed between the Apr 2026 snake scan/cosmics and production, so scan and cosmic constants transfer only as LMS-normalised ratios (gain_scan × LMS_prod/LMS_scan), never as absolute ADC/MeV; the absolute scale stays with the elastic pool. No tagger in the scan: beam energy from the accelerator (RCDB 2108.85 MeV). First data sample: run 23701 files .00000–.00080 every 10th, run 23744 files .00000/.00030/.00060, to be copied to /Users/L00338853/work/PRad/data/evio/.
- 2026-11-14 decision: P14 — no access to the scan log or logbook; the beam path (module ↔ event-number map, dwell boundaries) is reconstructed from the data (peak-occupancy module per event block, transitions excluded); HV/LMS segments from the LMS and pedestal time series in the same files. RCDB for 23701: 54 min, 54.3 M events, 17 kHz, 81 files, valid end; ~1 min per module across the 22-run scan.
- 2026-11-14 P14: plan written, notes/plan_p14_snake_scan_20261114.md (P14.1-P14.7; run 23701 readable with a derived hycal map via set_hycal_map, EPICS ptrans gives the beam path, 9 local files = 9 whole dwells on rows 17-18 covering T/O/M/H and 7 pool modules, endpoint edge at 18.5 k ADC, only ~25 LMS events locally; transfer as LMS ratio with pool-fixed normalisation, 2 % spread gate, comparison-only GAIN_PRIOR_FILE hook deferred behind P1b); deferred: Q1 farm pass of the full scan, Q2 extra 23701 files (15 GiB free), spot size if known
- 2026-11-14 decision: P14 child plan approved by the user; P14.1–P14.5 start on the local files (23701 ×9, 23744 ×3); P14.7 farm pass approved to run after P14.4 passes the 2 % gate (user submits; agents prepare scripts); no further 23701 files for the first pass; beam spot size unknown.

## Todo

- [x] User: approve this plan
- [x] P0a: copy BM and its sha1 to `output/energy_resolution/baseline_20261112/`; default reports of the five run
      types into `pre/`; diff against `ER/<RT>/summary.tsv`; rescan run 25320 into `scan_variants/c0/2p24GeV_H2/`
- [x] P0b: `--set`, the switches and the region / separation / pair / secondary outputs in BM; P0c:
      `compare_reports.py`; Slice 1 tests; rerun the default reports into `pre/` after `default_regression` passes
      (adds `regions_splits.tsv`)
- [x] P1: row a on five run types, b on `2p24GeV_H2`, g on three; `ER/improvement/P1_reconcile.md`
- [x] P2: `EVAL_REGIONS` and `MOLLER_PAIRS` reports on five run types; `ER/improvement/P2_regions.md`
- [x] P3: step A table; rescans `minsep12`, `rmodel0p3`; `ER/improvement/P3_secondaries.md`; hand row f to P1
- [x] P4: rows a, b; `ER/improvement/P4_switches.md`
- [x] P12: threshold from data and configuration, MeV per module, loss table; `ER/improvement/P12_threshold.md`
- [x] P13: run inventory; `notes/calibration_runs_snake_cosmic_20261112.md`
- [x] P5: budget note, separation table, re-baseline (D5), gate table; user approval of the gates and the U3 defaults
- [x] P15: `--set TIME_OFFSETS` (calibrated module times) in the report-side `add_hit`; gate scan 1-16 ns on five run types; best gate per run type vs legacy and vs raw; `ER/improvement/P15_tgate_calibrated.md`
- [x] P16: below-threshold (KTeV-style) correction, slice 1 measurement, slice 2 only if slice 1 is real (`ER/improvement/P16_threshold_corr.md`)
- [x] P17: asymmetric time window, tight early / loose late member window (`ER/improvement/P17_asym_window.md`)
- [x] P6: gain solve as a comparison arm (U2)
- [x] P8: secondary-pulse acceptance
- [x] P9: edge-aware position correction
- [x] P7: calibration beyond the pool
- [x] P10: onset-phase correction verdict
- [x] P11: boundary clusters and lead glass
- [x] Phase 2: control v3 built; final table in improvement/control_phase2_v3.md
- [ ] Phase 2: child plan, then execute, for each approved package in the D13 order; P11 and P14 when unblocked
