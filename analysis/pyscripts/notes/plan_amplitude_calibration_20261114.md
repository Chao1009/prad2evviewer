# Plan: amplitude-based energy calibration track (comparison against control v3) — 2026-11-14

**Status: parked 2026-11-14 after P1b (user). Resume at P2 with the cfd7/sum5 columns; see Log.**

Working dir `analysis/pyscripts` (branch `waveform_jz`); paths are relative to it. BM =
`waveform_analysis/energy_resolution_benchmark.py`; ER = `waveform_analysis/energy_resolution`; XS =
`crosstalk_study/crosstalk_study.py` (its 11-07 version with the offset solver is
`crosstalk_study/previous_20261107/crosstalk_study.py`); V3 = `output/energy_resolution/control_phase2_v3_20261114/`.
Control = v3 (`ER/improvement/control_phase2_v3.md`; the chain step by step in
`waveform_analysis/writeup/waveform_analysis_prad2_x17.tex` L1146-1352). Parent protocol: D1-D4 and D12 of
`notes/plan_energy_resolution_improvement_20261112.md` (the "parent plan"). Environment: `source setup_env.sh`,
`MPLBACKEND=Agg`.

## Time estimate

| | estimate |
|---|---|
| Agent time | about 3.5 agent-days (P0 0.5, P1 0.5, P2 0.1, P3 0.5, P4 0.5, P5 0.75, P6 0.5), plus review rounds |
| CPU | scans 10.3 CPU-h (run 25320 study 0.85 h + full rescan of the 12 H2 runs, 102 EVIO files, 9.5 h; Ta later, +1.9 h); reports about 45 report calls, 1-10 min each at 8 workers, about 2.5 h wall sequential |
| Wall | 4-5 working days; P1's raw-versus-parabolic result is reported to the user before P2 starts |
| Critical path | P0 code (0.5 d) → P1 estimator study on run 25320 (about 1 h compute + 0.5 d) → P2 full rescan (2.5 h wall at 4 lanes, 4.8 h at 2) → P3 CFD offsets → P4 gate scan and scale → P5 tables and final stack → P6 comparison |

Cost basis: `output/energy_resolution/logs/scan_<run>.log` (270-390 s per EVIO file, one core per scan), V3
`logs/walltimes.txt` (reports 46-547 s), parent-plan P10 log (6 rescans + 20 reports ≈ 5 CPU-h). A scan is
single-process (BM scan loop, no `NPROC`), so lanes are shell processes; 18 CPUs allow 4 lanes plus one report.

## Open Questions

Material decisions for the user; each has a recommended answer. The plan proceeds on the recommendations if the
user approves it as is.

- Q1 and Q2 are closed by the user (U1, U2 in Decisions). P1 settles the raw-versus-parabolic sub-choice by its
  phase-sensitivity study and reports it; no user action is needed unless the two disagree by ≥ 2 SE on
  `res_cb_oos` and the phase-stable one is the worse.
- Q3 (templates). No refit: the estimators do not use the bank, and the bank's unit-peak normalisation is already
  "amplitude normalisation". The bank stays what it is for pulse identification (MF), the CFD search window and
  the integral control. Recommended: no.
- Q4 (gates). CFD times differ from MF times (heavier late tail, 4.75 vs 1.00 % of |dt| > 10 ns raw at 2.24 GeV,
  `crosstalk_study/notes` findings), so the v3 control gates are not transferable. Recommended: re-scan the
  gates (one report call per run type carries the gate list) and report both the control-gate row and the best
  gate by the parent's 1-SE rule.
- Q5 (Ta). Recommended: not in the first pass. The Ta runs (3 runs, 25 files, 1.9 CPU-h) are rescanned and
  calibrated (own M1 table, Ta gate scan, borrowed 2.24 GeV map) only after the H2 verdict in P6. First-pass
  comparison = the four H2 run types.
- Q6 (CFD definition). Recommended: the raw-sample 50 % linear CFD of the crosstalk study (XS:645-663,
  `CFD_FRAC` 0.5, search ± `CFD_SEARCH_NS` 12 ns around the main pulse's template peak, linear interpolation, no
  smoothing), computed in BM's scan `flush()` from `Y`; template-free like the estimator and the CFD the study's
  figures use. The alternative is the WaveAnalyzer log-normal CFD (`prad2dec/src/WaveAnalyzer.cpp:515-542`),
  which needs the C++ path. WaveAnalyzer fallback hits (clipped PbWO4 and lead glass, 1.9e-4 of hits, `h_a` NaN)
  keep their `pk.time`, a native CFD, either way.

## Goal

Build and calibrate from scratch a second energy chain whose hit energy is a template-free pulse amplitude and
whose hit time is a CFD time with its own module offsets, and compare it with control v3 on the same runs, splits
and metrics, out of sample. The product is a verdict table, the diagnostics listed in D6, and a note. Non-goals:
changing any v3 default; production code or `../../database/` changes; template-bank work; O and T ring constants
(P14 of the parent); a new pile-up fitter (`pileup_analysis/` and `exotic_study/` are not touched).

Expected change size: BM about 150 lines (scan: three peak columns and the CFD time computed in `flush()`; report: `ENERGY_DEF`,
`TIME_DEF`, `TIME_OFFSETS_FILE`, `AMP_SCALE` with their `OVERRIDABLE` entries; no existing output changes); one
offsets tool ≤ 120 lines (port of `selfref_offsets`); one estimator-study tool ≤ 100 lines; one ratio-diagnostic
tool ≤ 80 lines (may extend `amp_vs_integral.py`); 5-6 shell drivers in the V3 pattern; 2 table families under
`ER/tables/`; one offsets json; one results note.

## Approach

- Same benchmark, two switches. `ENERGY_DEF` ("integral" default = v3; "rawpk", "parpk") selects
  which stored per-hit column is the hit energy and `TIME_DEF` ("mf" default; "cfd") which time column; everything
  downstream (gate, clusterer, gains, maps, tables, fits) is re-derived by the existing code. Rejected: a second
  script (the comparison needs identical selection and fit code); editing constants per variant (not
  reproducible).
- One rescan that stores the peak columns plus the CFD time next to the v3 hit arrays, so every later choice is
  report-level: the `integral + cfd` arm isolates the timing effect, and the integral control re-run on the same
  scan (`v3r`, 2.24 GeV only, where the bank changed) gives the bank check.
- Estimator study on one run before the full rescan (P1), because the raw-versus-parabolic choice and the phase
  sensitivity are the facts that decide what the rescan must carry.
- Calibration order copied from the V3 drivers (`stage1_p16.sh` → `stage2_m1.sh` → `stage3_v3.sh`): first the
  scale and gate, then the below-threshold table, then the M1 table fitted on the chain with that table, then the
  final stack with rings; gains, position map, shape factor and ring gains are solved inside each report.
- Verdict by the parent's D3, paired per (run, split), out of sample, against v3 on disk (headline) and `v3r`
  (bank check).

## Decisions

User decisions (2026-11-14, binding, relayed by the coordinator):

- U1 (estimator). The hit energy of the track is the raw peak sample, or the local parabolic peak through that
  sample and its two neighbours, of the baseline-subtracted waveform of the identified pulse, independent of the
  template; timing is CFD. Why the template-fit amplitude `h_a` is not used: `h_e = gain · iint_m · ph(φ) · h_a`
  (BM:937), so `h_a / h_e` is constant per module up to the 1.8 % phase wiggle
  (`ER/improvement/tools/amp_vs_integral.py`: median 0.140 per sample, 4.5 % width), and an `h_a` chain is the
  `PHASE_CORR=False` arm already measured (+0.055 ± 0.003 pp at 2.24 GeV, +0.050 at 3.49 GeV, parent P10) with
  the gains re-absorbing `iint_m`; it reproduces v3 otherwise. Raw versus parabolic: P1 settles it by the phase
  study (D2); the parabola is expected to remove most of the 4 ns sampling phase dependence the raw sample
  carries. Rejected: a fixed-shape single-sample fit anchored at the CFD time (a one-shape template fit, half-way
  back to v3); a fixed-window sample sum (a third chain, not the user's question).
- U2 (overlapping pulses). For a hit whose waveform holds a close pile-up (the joint fit found ≥ 2 pulses,
  `h_k ≥ 2`, BM:937), the amplitude is the raw waveform value at the primary's peak sample minus the contribution
  of the *other* fitted pulses at that sample: Σ_{k ≠ main} a_k · s_m(t_i − t0_k), with (a_k, t0_k) the
  multi-pulse fitter's output for the channel (`run_ref` / `MF`, BM:311 ff., 708-716) and s_m the module template
  the fitter already uses; the primary's own template is not used. The subtraction acts on the waveform value
  (and on the two neighbours for the parabola), so the single-pulse case reduces to the plain raw or parabolic
  peak. Diagnostic to carry: the amplitude with and without the subtraction for the multi-pulse subset (`pu`
  flag at cluster level, `h_k ≥ 2` at hit level) against Δt to the secondary (`h_t` of the `h_extra` hit of the
  same module and event).

Lead decisions:

- D1 (where the amplitude comes from today). The scan fits all pulses of a waveform jointly (Gauss-Newton on
  (a_k, t0_k), merge at `MERGE_NS`, shoulders with BIC; `MF` class BM:311 ff.) with the module's unit-peak spline
  template; `h_a` is that amplitude in baseline-subtracted ADC counts, `h_e = gain·iint·ph·h_a` (BM:937; `gain` =
  `Module.energize(1.0)` MeV per integral unit from the runinfo calibration, `iint` = the template's walk-out
  integral, 5.0-9.1 across modules at 2.24 GeV, `ph` = `ptab` BM:853, 1.8 % peak to peak). Both are in every
  scan npz (`h_a` float32, NaN only for the 1.9e-4 WaveAnalyzer fallback hits, BM:915-916, 940-945). The
  "identified pulse" of a module in an event is the largest-`h_e` hit in the 100-200 ns window (`lexsort`,
  BM:974); every later pulse of the same module carries `h_extra` (BM:976), is skipped on the legacy path and kept
  as a hit on gated paths (BM:1214). Overlap is not double counted in `h_a` (joint fit); it is in the raw
  estimators before the U2 subtraction. The per-module amplitude/integral ratio is `1/(gain·iint·ph)` by construction; the diagnostic
  the track reports is therefore `E_est / h_e` per module and versus energy, where `E_est` is the template-free
  estimator (D6).
- D2 (the estimators and their storage). The scan's `flush()` (BM:918-936) holds `Y` (raw uint16 samples, 4
  ns), the pedestal and the fitted `t0` per row; the columns are computed there, next to `h_a`, and written by
  the npz writer (BM:979-990) as `h_apk_raw` and `h_apk_par` (U1, with the U2 subtraction), `h_apk_par0` (no
  subtraction, the U2 diagnostic; equal to `h_apk_par` for single-pulse hits) (float32, ADC; NaN for fallback
  hits and clipped waveforms) and `h_tcfd` (ns from sample 0, the MF time base, NaN when no crossing or clipped,
  XS:662-663). The other pulses' model comes from the fitter's (a_k, t0_k) and the module template, as the
  `model` array of XS:636-637 is built but without the main pulse; it is subtracted from `y` first, and the peak
  sample is then located on that subtracted waveform in the XS:654 window (± 12 ns around the main pulse's
  template peak), so a secondary's sample cannot be picked; `h_apk_par0` uses the same sample on the unsubtracted
  `y`. The subtraction term at the peak sample is the same for the raw and the parabolic peak, so `h_apk_par −
  h_apk_par0` is the U2 diagnostic for either choice. P1 and the full rescan store all three columns (float32,
  about 80 MB per column per run at 2.24 GeV, 20 M in-window hits; four new columns, about 4 GB over the 12 runs
  before compression). `h_t`, `h_e` and `h_a` are unchanged, so the same npz serves the integral control.
  `ENERGY_DEF` takes "integral" (default), "rawpk" or "parpk".
  - 2026-11-14 (P0): CFD computed on the waveform with the other fitted pulses subtracted (extras within ~24 ns otherwise inherit the primary's time); single pulses unchanged.
- D3 (hit energy of the track and the pre-gain scale). The report builds `e = AMP_SCALE[m] · column` at hit load
  (BM:2315, where `h_e` is read today). Every threshold before the gain solve is in MeV on the raw hit energy
  (`MIN_CENTER_E` 30, `SIGNAL_MIN` 500·kE, `SEED_MIN` 200·kE, `min_module_energy` 1.0 from
  `database/reconstruction_config.json:15`; BM:110-117, 1197), so the track needs a pre-gain scale. Two passes:
  pass 1 uses one global factor s = the median of `h_e / column` over the pool hits of the run type (about 1.16
  MeV/ADC for `h_a` at Ta; the P1 tool reports it per column), solves the gains, and writes
  `ER/tables/amp_scale/<RT>_<column>.tsv` (module, s·G_m) for the solved modules; pass 2 reads that table as
  `AMP_SCALE` (modules without a row keep s) and is the chain all later stages run on. Reason: only one global
  number taken from `h_e` enters pass 1, and pass 2 replaces it per module from the track's own gain solve, so
  the pool calibration is the track's own and every MeV constant keeps its meaning. Rejected: the
  per-module pseudo-scale `gain·iint` (reuses the integral calibration); a single global factor throughout
  (pre-gain thresholds then vary by the 1.8× template spread across modules).
- D4 (CFD time and offsets). `TIME_DEF="cfd"` reads `h_tcfd` instead of `h_t` at BM:1202-1203; hits with NaN
  `h_tcfd` fall back to `h_t` (the fallback count is a diagnostic). Offsets: `time_calib.json` is on the
  WaveAnalyzer-CFD basis and global; the track re-derives offsets on its own CFD from the P2 scan with the
  pair-network method of `selfref_offsets` (XS 11-07 version :334, writer :665-678: per (module, seed-module)
  pair the median dt in shrinking windows 20 → 5 → 3 → 3 ns, pairs ≥ 20 hits, √n-weighted least squares), ported
  to BM hit arrays (seed = largest hit of the event, members = the other in-window hits of the same event, legacy
  clusters), gauge re-zeroed to W495 (module id 1495) instead of the median. One table per H2 run type for the
  comparison, then one global table from the three H2 run types; the plan uses the global table if the per-run-
  type tables agree within 0.3 ns rms (else per run type, flagged). Modules with n < 50 get offset 0 and
  `resolved = false`. Format = `time_calib.json` (keys `reference_module`, `modules[]` with `name`, `module_id`,
  `offset_ns`, `offset_error_ns`, `resolved`), written to `ER/tables/time_offsets/cfd_<scope>.json`, read through
  the new `TIME_OFFSETS_FILE` setting at the loader BM:1043-1050. Self-check: the same tool on `h_t` must
  reproduce `time_calib.json` to the 1.7-2.4 ns module scatter the 11-07 study found.
- D5 (what is recalibrated, in order; stages mirror V3 `stage1-3.sh`). Items marked "auto" are solved inside every
  report with no table.

  | step | v3 source | track | re-derived by |
  |---|---|---|---|
  | pulse acceptance `NSIGMA, MIN_AMP` (BM:91) | ADC amplitude | unchanged (acts on the fit) | — |
  | onset-phase table `ptab` | integral only | not applied | — |
  | time offsets | `time_calib.json` (MF basis) | CFD table | P3 |
  | gate | 20/16/16/16/4 ns | own scan 4-24 ns | P4 |
  | pre-gain scale, thresholds | `gain` MeV/integral | `AMP_SCALE` table (D3) | P4 |
  | log weight W0 (`LOG_WEIGHT_BY_RT`, BM:220) | 3.658 / 3.707, dimensionless | kept unless the pass-2 pool fit moves its optimum: a 5-point check at 2.24 and 3.49 GeV, full scan only if the minimum moves | P4 |
  | below-threshold table (`below_thr`, BM:1148) | `ER/tables/below_thr/*.npz`, T_m MeV | refit on the track's `THR_DUMP` with `p16_table.py` | P5 |
  | M1 phase residual (`phase_resid`, BM:2279-2299) | `ER/tables/phase_resid/M1_<RT>.tsv` | refit on the track's `PHASE_STUDY` dump with `p10_phase_fit.py` (mandatory: larger phase dependence expected) | P5 |
  | gains (`solve_gains`, BM:1547), rings (`RING_CALIB_BY_RT`), position map + shape factor (`fit_pc` BM:1759, `fit_extra` BM:1846) | auto | auto | each report |
  | CB fit, tails, `SEP_GAP_MEV`, `SEP_CUT_MEV` | E' scale | unchanged (scale-free once gains map to E') | — |

  The tables live under `ER/tables/amp/<column>_<time>/` (below_thr, phase_resid) so the v3 tables are untouched;
  the report selects them through the existing `*_BY_RT` dicts given on the command line.
- D6 (comparison contract and diagnostics). Runs and splits: the parent D4 H2 set (12 runs; 24937 stays out; Ta per Q5), the
  rescan reusing the split lists recorded in the on-disk npz (`output/energy_resolution/variants/P5/scan_lane.sh`
  pattern). The `3p49GeV_H2_old` runs are on disk as lane files (`run_<run>_{a,b,c}.npz`, 8 files for 3 runs);
  the rescan writes one npz per run with the union of the lane split lists (the driver globs `run_<run>*.npz`),
  so the rescan holds 12 npz. Metrics: parent D1 (`res_cb`, `sep_pi`, `tail_hi3/lo3`, `tail_hi_f/lo_f`, `eff_cut`, `leak_cut`,
  yield), out of sample (even/odd cross fit), paired per (run, split) with `compare_reports.py`. Control: v3 on
  disk (`V3/v3/<RT>`), with `v3r` (integral chain re-run on the 2.24 GeV P2 scan, D7) as the bank check. Rows per
  run type: track at the v3 control gate, track at its own best gate, and the timing arm `integral + cfd` (the
  v3 chain on CFD times and offsets, which isolates the timing effect from the amplitude effect). Verdict (parent D3): the track **wins** when Δ`res_cb_oos` ≤ −2 SE at 3.49 GeV with the same sign
  at `3p49GeV_H2_old`, or at ≥ 2 of the 3 H2 energies, with `tail_hi_f` not worse by > 2 SE, `eff_cut` not down
  and `tail_lo_f` not up by > 2 SE, yield loss ≤ 1 %; **loses** when Δ ≥ +2 SE under the same guards; otherwise
  **null**. Diagnostics reported regardless: `E_est/h_e` per module (median, IQR) and in 8 energy bins (tool on
  the npz); the first-harmonic phase amplitude of E/E' before and after M1 (from the `PHASE_STUDY` dump, as
  `P10_phase_corr.md` reports it); `res_cb` and `tail_hi_f` on the `all` and `nopu` arms and the `pu` share
  (`splits.tsv`); the U2 subtraction size (`h_apk_par` − `h_apk_par0`) for `h_k ≥ 2` hits against Δt to the
  secondary; the tail columns and the CB tail parameters of `summary.tsv`; the CFD fallback fraction and the
  offset table's resolved count and rms against `time_calib.json`.
- D7 (`v3r`). For `2p24GeV_H2` the on-disk scans used bank 15dfd867, which is gone (parent D5); the rescan uses
  7f125b75, so `h_e` changes there and the M1 table must be refit for it (write-up open item 8). `v3r` = the
  integral chain on the 2.24 GeV P2 scan with a refit M1 table (`ENERGY_DEF="integral"`, `TIME_DEF="mf"`; one
  report plus one table fit). Its paired Δ to v3 on disk is the bank effect and is quoted next to the track's Δ.
  The other run types' banks are unchanged, and `scan_additive` covers their reproducibility.
- D8 (provenance). Variant trees `output/energy_resolution/variants/AMP/<stage>/<RT>/`, rescans
  `output/energy_resolution/scan_variants/amp_cfd/<RT>/run_<run>.npz` with the override string in the npz
  (`overrides` key) and `bank_sha1`; each stage driver records its code sha1 as in V3; every `delta.tsv` names the
  control dir. The P1 study tree is `variants/AMP/p1_25320/`.
- D9 (byte-reproducibility of v3). Every new constant defaults to the v3 behaviour (`ENERGY_DEF="integral"`,
  `TIME_DEF="mf"`, `TIME_OFFSETS_FILE=None` → `time_calib.json`, `AMP_SCALE=None`); the P0 checkpoint runs
  `check_identical.sh` between a default report on the on-disk scans and `V3/v3_defaults`.

## Outline

`waveform_analysis/energy_resolution_benchmark.py`:
- constants `ENERGY_DEF = "integral"`, `TIME_DEF = "mf"`, `TIME_OFFSETS_FILE = None`, `AMP_SCALE = None` (path of
  the per-module table or None), added to `OVERRIDABLE`
  (BM:232-236; pattern of `BELOW_THR_BY_RT`).
- scan `flush()` (BM:918-936): extend `run_ref` (BM:708-716) to also return the baseline it subtracted
  (`np.where(far, b_ref, ped)`); after the fit, compute the D2 columns from `y = Y - baseline` with the XS:650-663
  code (peak search window, half crossing, parabola through `p-1, p, p+1`); append to the hit arrays next to `a`
  (BM:933-937); writer BM:979-990 adds the keys. Invariant: `h_t`, `h_e`, `h_a`, `h_k`, `h_extra` byte-identical
  to a scan without the switch (the columns are pure additions).
- report hit load (BM:2315): `e = h_e` for "integral", else `AMP_SCALE[h_mod] · npz[column]`; a missing column
  is an error naming the npz. `TIME_DEF`: `t = h_tcfd` with NaN → `h_t` at BM:1202-1203; the loader at
  BM:1043-1050 takes `TIME_OFFSETS_FILE`. `AMP_SCALE` table writer after `solve_gains` (module, scale) when
  `ENERGY_DEF != "integral"`, in `gains.tsv`'s directory as `amp_scale.tsv`. Gotcha: the `SECONDARY_RULE` ratio
  uses `h_a` (BM:1118) and stays on `h_a`.

`waveform_analysis/energy_resolution/improvement/tools/amp_offsets.py` (new, ≤ 120 lines): `amp_offsets.py
<scan dir> [--time h_tcfd|h_t] [--out json] [--ref W495]`, the port of `selfref_offsets` (D4); output columns
`module n offset_ns se` as a tsv plus the json.

`waveform_analysis/energy_resolution/improvement/tools/amp_estimators.py` (new, ≤ 100 lines): per npz column,
the ratio to `h_e` per module and per energy bin, its phase dependence (first harmonic versus `h_t % 4 ns`), and
the `E_est/h_e` table of D6; may extend `amp_vs_integral.py` rather than copy it.

Drivers `output/energy_resolution/variants/AMP/stage{0_code,1_scan,2_offsets,3_gate_scale,4_tables,5_final}.sh`
and `lib.sh` in the V3 pattern (sha1 per stage, `walltimes.txt`).

## Slices

1. P0 code + P1 study (run 25320 only).
   Tests: `v3_identical: default report of 2p24GeV_H2 on the on-disk scan -> check_identical.sh against
   V3/v3_defaults passes`; `scan_additive: rescan of 25320 split 0 with the new columns -> h_t/h_e/h_a identical to
   `output/energy_resolution/scan_variants/c0/2p24GeV_H2/run_25320.npz` (bank 7f125b75) for that split, new columns finite on ≥ 99.9 % of PbWO4 hits`; `estimator_table:
   amp_estimators.py on the P1 npz -> one row per column with ratio median, IQR, first-harmonic amplitude`;
   `phase_study: first-harmonic amplitude of E_est/h_e versus onset phase -> rawpk vs parpk reported, the
   phase-stable one named`; `overlap_diag: h_apk_par − h_apk_par0 on h_k ≥ 2 hits -> table against Δt to the
   secondary, zero for h_k = 1`; `p1_reports: ENERGY_DEF in (rawpk, parpk) with TIME_DEF=cfd, plus integral with
   TIME_DEF in (mf, cfd), on 25320 -> res_cb_oos finite on 8 splits, delta.tsv against integral_mf`.
   Checkpoint: `python3 waveform_analysis/compare_reports.py output/energy_resolution/variants/AMP/p1_25320/integral_mf
   output/energy_resolution/variants/AMP/p1_25320/parpk_cfd` prints paired Δ on 8 splits.
2. P2 rescan + P3 offsets.
   Tests: `rescan_complete: 12 npz, bank_sha1 per run type equals the on-disk bank, splits equal the baseline
   lists`; `offsets_selfcheck: amp_offsets.py --time h_t -> rms against time_calib.json ≤ 2.5 ns, W495 = 0`;
   `offsets_cfd: resolved count ≥ 1100, per-run-type tables vs global rms reported`.
   Checkpoint: `ls output/energy_resolution/scan_variants/amp_cfd/*/run_*.npz | wc -l` = 12 and
   `ER/tables/time_offsets/cfd_global.json` exists.
3. P4 gate and scale, P5 tables and final stack, P6 comparison.
   Tests: `gate_scan: one report per H2 run type with gates 4 6 8 12 16 20 24 -> best gate by the 1-SE rule`;
   `scale_pass2: amp_scale.tsv rows ≥ pool module count, pass-2 pool μ within 1 % of E'`; `w0_check: 5 values at 2.24 and 3.49 GeV -> v3's W0 kept
   or the new minimum reported`; `p16_amp: table built, THR_DUMP share of empty windows
   reported`; `m1_amp: first harmonic after M1 ≤ half of before`; `final: v3-equivalent stack at the control gate
   and the best gate, the integral + cfd arm, v3r at 2.24 GeV -> delta.tsv each; verdict per D6`.
   Checkpoint: `ER/improvement/AMP_amplitude_track.md` has one verdict row per H2 run type with Δ ± SE and the
   D6 diagnostics.

## Verification

`compare_reports.py V3/v3/<RT> output/energy_resolution/variants/AMP/final/<RT>` for the four H2 run types, plus
the `v3r` pair, each with `verdict_res_cb_oos` filled; `check_identical.sh` of a default report against
`V3/v3_defaults` still passes after all edits; the note lists the win/lose/null verdict, the timing arm,
and the D6 diagnostics.

## Risks

- The raw estimators carry the full sampling-phase dependence (4 ns samples); if M1 cannot reduce it the track
  loses on phase alone. P1 measures this before the rescan; a per-module phase table in amplitude units would be
  an `Amend parent:` package, not part of this plan.
- CFD has a heavier late tail than MF; the best CFD gate may be wider and the Ta 4 ns gate may not transfer. The
  `integral + cfd` arm separates this from the amplitude effect.
- Pile-up: the U2 subtraction leans on the fitter's secondary pulses (missed or merged secondaries, `MERGE_NS`
  6.5 ns, stay in the peak); the `pu`/`nopu` split and the U2 diagnostic bound the effect; Ta (34 % `pu`)
  is deferred (Q5).
- `v3r` differs from v3 at 2.24 GeV by the bank; the Δ is quoted, not corrected.
- Agent time assumes the V3 drivers transfer with path edits; a `p16_table.py` or `p10_phase_fit.py` integral
  assumption would add about half a day.

## Scope Notes

In: P0-P6 as written, H2 first pass; Ta after the verdict (Q5). Out (user constraints): nothing becomes a default; no git
(the user commits); `pileup_analysis/pileup_rates.py` (uncommitted diff) and `exotic_study/` untouched; no
production or `../../database/` changes; no template refit (Q3). Out, parent decisions not reopened: U1-U7 and
D1-D3 of the parent plan; Møller rejection and non-linearity stay off; refinement stays off (U2).

## Assumptions

- Scan cost scales with EVIO files at 270-390 s per file, one core (checked: `logs/scan_<run>.log`).
- `Y` and `t0` are available in `flush()` for every fitted row (checked at BM:918-936); the baseline is returned
  by `run_ref` after the P0 edit (Outline). The executor confirms the clipped-row path.
- `p16_table.py` and `p10_phase_fit.py` read the dumps without an integral-unit assumption (checked in P0 by
  grep; else Risks).
- The `report` gate list runs all gates in one call (checked: V3 `stage3_v3.sh`).

## Work Packages

### P0: switches and columns
- Goal: BM carries `ENERGY_DEF`, `TIME_DEF`, `TIME_OFFSETS_FILE`, `AMP_SCALE` and the scan columns of D2, with v3 byte-identical at the defaults; done when `v3_identical` and `scan_additive` pass
- Owns: BM, `variants/AMP/stage0_code.sh`, `lib.sh`
- Deps: none
- Stages: execute
- Status: done
- Artifacts: ->
- Cost: 0.5 agent-day; one single-split scan (6 min); one report (2 min).

### P1: estimator study on run 25320
- Goal: measured phase dependence, U2 overlap diagnostic, ratio to `h_e`, and paired `res_cb_oos` on 8 splits for rawpk and parpk at CFD timing plus the integral at both timings; done when the raw-versus-parabolic choice is made and reported
- Owns: `tools/amp_estimators.py`, `variants/AMP/p1_25320/`, `ER/improvement/AMP_P1_estimators.md`
- Deps: P0
- Stages: execute, write
- Status: done
- Artifacts: `ER/improvement/AMP_P1_estimators.md`
- Cost: 0.5 agent-day; rescan of 25320 (8 files, 0.85 CPU-h, 7 min wall at 8 lanes by split subset); 4 single-run reports (about 1 min each).

### P2: full rescan
- Goal: the 12 H2 runs rescanned with the peak columns and `h_tcfd`, same splits and banks as the on-disk scans; done when `rescan_complete` passes
- Owns: `scan_variants/amp_cfd/`, `variants/AMP/stage1_scan.sh`
- Deps: P1
- Stages: execute
- Status: parked: user decision 2026-11-14 (revisit later)
- Artifacts: `output/energy_resolution/scan_variants/amp_cfd/<RT>/run_<run>.npz`
- Cost: 9.5 CPU-h; 2.5 h wall at 4 lanes (the longest runs 25228, 25260 and 24590 start first); 0.1 agent-day. Ta later: +1.9 CPU-h.

### P3: CFD module offsets
- Goal: per-run-type and global CFD offset tables, W495 gauge, self-checked against `time_calib.json` on `h_t`; done when `offsets_selfcheck` and `offsets_cfd` pass
- Owns: `tools/amp_offsets.py`, `ER/tables/time_offsets/`, `variants/AMP/stage2_offsets.sh`
- Deps: P2
- Stages: execute
- Status: parked: user decision 2026-11-14 (revisit later)
- Artifacts: `ER/tables/time_offsets/cfd_{0p73GeV_H2,2p24GeV_H2,3p49GeV_H2,global}.json`, `cfd_offsets.md`
- Cost: 0.5 agent-day; minutes of compute.

### P4: gate scan, pre-gain scale, log weight
- Goal: best CFD gate per H2 run type, the pass-2 `AMP_SCALE` tables, and the W0 check; done when `gate_scan`, `scale_pass2` and `w0_check` pass
- Owns: `ER/tables/amp_scale/`, `variants/AMP/stage3_gate_scale.sh`, `variants/AMP/s3_*/`
- Deps: P3
- Stages: execute
- Status: parked: user decision 2026-11-14 (revisit later)
- Artifacts: `ER/improvement/AMP_P4_gate_scale.md`
- Cost: 0.5 agent-day; 4 gate-scan reports (7 gates, about 10 min each), 4 pass-1 and 4 pass-2 reports, up to 10 W0 reports (2.24 and 3.49 GeV, D5): about 1-1.5 h wall.

### P5: below-threshold and M1 tables, final stack
- Goal: track tables for P16 and M1, then the v3-equivalent stack (rings, map, shape factor solved in the report) at the control gate and the best gate, plus the `integral + cfd` arm and `v3r`; done when `final` passes
- Owns: `ER/tables/amp/`, `variants/AMP/stage4_tables.sh`, `stage5_final.sh`, `variants/AMP/{s4_*,final,v3r,arms}/`
- Deps: P4
- Stages: execute
- Status: parked: user decision 2026-11-14 (revisit later)
- Artifacts: `ER/tables/amp/<column>_cfd/{below_thr,phase_resid}/`, `variants/AMP/final/<RT>/delta.tsv`
- Cost: 0.75 agent-day; about 22 reports (THR_DUMP 4, PHASE_STUDY 4, final 4 × 2 gates, timing arm 4, v3r 1 + 1 table refit): about 1 h wall at two at a time.

### P6: comparison note
- Goal: the D6 verdict table and diagnostics for the four H2 run types, against v3 and `v3r`, with the Ta plan; done when the note is reviewed
- Owns: `ER/improvement/AMP_amplitude_track.md` (+ `.tsv`), figures under `ER/improvement/figures/AMP_*`
- Deps: P5
- Stages: write
- Status: parked: user decision 2026-11-14 (revisit later)
- Artifacts: `ER/improvement/AMP_amplitude_track.md`
- Cost: 0.5 agent-day.

## Log

- 2026-11-14 P0: switches `ENERGY_DEF`/`TIME_DEF`/`TIME_OFFSETS_FILE`/`AMP_SCALE`/`PEAK_COLUMNS` and columns `h_apk_raw`/`h_apk_par`/`h_apk_par0`/`h_tcfd` added, v3 byte-identical at defaults (`check_identical.sh` against `v3_defaults`: all 7 tables identical, final code), scan_additive passes on all 8 splits of 25320; CFD is computed on the other-pulse-subtracted waveform (deviation from D2); deferred: AMP_SCALE comment trim, ring-scoped `amp_scale.tsv` guidance (P4 note)
- 2026-11-14 P1: raw peak beats 3-point parabola on phase dependence (E/E' first harmonic 11.9 vs 14.2 permille) and on `res_cb_oos` (2.284 vs 2.473 %, Δ +0.189 ± 0.021 pp for the parabola, paired out of sample); both are far behind the integral without tables (1.827 %, raw +0.456 ± 0.011 pp); the U2 subtraction brings overlapped peaks to within about 5 % of the single-pulse value; CFD and MF times agree per hit (robust sigma 0.13-0.23 ns); verdict goes to the user before P2; deferred: cause of the non-phase excess and a cluster-level U2 test (P5/P6)
- 2026-11-14 P1b: the raw-peak excess of 2.00 ± 0.11 percent² (1.37 % in quadrature) is 41 % calibration-recoverable (phase table 0.83 ± 0.09), the rest a 1/√A single-sample fluctuation (0.44), a shape term (0.18), pile-up members (0.2) and about 0.3 unexplained, with noise only about 1 %; best template-free variant `cfd7` with tables gives 1.935 % against 1.866 % for the integral, about 0.14 pp behind v3; recommend stop, no P2 as planned; deferred: out-of-sample check of the sum variants (needs a second file or P2 data), isolating the 0.44 gain-solve cost.
- 2026-11-14 decision: user stopped the amplitude track after P1b ("stop here for me, we can revisit later"); P0/P1/P1b stand as the documented comparison; switches ENERGY_DEF/TIME_DEF/PEAK_COLUMNS and the P1/P1b scans (scan_variants/amp_p0, amp_p1, amp_p1b) stay for a later resume; resume point = P2 with h_c7/h_s5/h_s7 columns (cfd7/sum5 estimators) if revisited.
- 2026-11-14 lead: figures waveform_analysis/figures/amp_vs_integral.{png,pdf} (fitted amplitude vs integral) and amp_raw_vs_integral.{png,pdf} (raw peak, cfd7, sum5 vs integral; scripts in energy_resolution/improvement/tools/)

## Todo

- [ ] User: approve this plan and answer Q3-Q6
- [x] P0: switches, scan columns, `v3_identical`, `scan_additive`
- [x] P1: 25320 rescan with the peak columns, `amp_estimators.py`, 4 reports, `AMP_P1_estimators.md`; raw vs parabolic settled
- [x] P1b: excess diagnostic (AMP_P1b_excess.md)
- [ ] P2: rescan of the 12 H2 runs into `scan_variants/amp_cfd/` (4 lanes)
- [ ] P3: `amp_offsets.py`, self-check on `h_t`, CFD tables per run type and global, `cfd_offsets.md`
- [ ] P4: gate scan, pass-1/pass-2 scale, W0 check, `AMP_P4_gate_scale.md`
- [ ] P5: `THR_DUMP` + `p16_table.py`, `PHASE_STUDY` + `p10_phase_fit.py`, final stack at two gates, `integral + cfd` arm, `v3r` (2.24 GeV)
- [ ] P6: `AMP_amplitude_track.md` with the D6 verdict and diagnostics; Ta rescan and stages scheduled or dropped (Q5)
