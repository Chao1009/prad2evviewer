# Plan: pile-up fitter validation, T0-T9 (shape error, flat tops, overlay truth, rates, benchmark) — 2026-11-16

Working dir `analysis/pyscripts` (branch `waveform_jz`); paths are relative to it unless they start with `/`. WA =
`waveform_analysis`; PA = `WA/pileup_analysis`; PR = `PA/pileup_rates.py`; BM = `WA/energy_resolution_benchmark.py`;
MS = `output/pileup_rates/monoedge_sample/` (`old/` rule off, `new/` rule on, `splits.txt`; the 2-file sample of
`PA/run_monoedge_sample.sh`); SCAN = `output/pileup_rates/scan/` (13 run npz of `PA/run_pileup_scan.sh`); EVIO =
`/Users/L00338853/work/PRad/data/evio/`; TEX = `WA/writeup/waveform_analysis_prad2_x17.tex`; V3 =
`output/energy_resolution/control_phase2_v3_20261114/`. Environment: `source setup_env.sh` (`PRAD2_DATABASE_DIR`,
`PYTHONPATH=build/python`), `MPLBACKEND=Agg`. Previous evidence: `PA/monoedge_fit_comparison_20261116.pdf/.tsv`
(script `PA/monoedge_fit_comparison.py`), `PA/checks_monoedge_removed_W493_timing_20261116.pdf`,
`PA/pileup_multiplicity_20261115_monoedge.pdf`. User request: "Give me a full plan of testing everything", then "start
and for each stage, generate examples to help demonstration and for me to inspect". The approach below is decided and
approved; this plan pins paths, outputs, tests and gates.

## Time estimate

| | estimate |
|---|---|
| Agent time | T0 10 min, T1 30 min, T2 20 min (shares T1's pass), T3 1-2 h, T4 15 min, T5 1 h, T6 1 h, T7 30 min, T8 2 h run, T9 overnight; about 1.5 agent-days plus review |
| CPU, local | T0-T5: refits of at most a few thousand waveforms (minutes); T6: 2 files x 13 runs x 1 extra variant of the sample scan (about 20 min wall at 16 lanes, as the old/new sample); T8: 2.24 GeV reports 2 x 45 min, 3.49 GeV 2 x 20 min, 0.73 GeV 2 x 20 min; T9: the full scan (13 runs, all local splits, 8 in parallel: several hours) plus the report |
| Critical path | T0/T1/T2 (in flight, parallel) -> gate 1 (user) -> T3 -> T4 and T5 (parallel) -> T6 (T7 in parallel) -> gate 2 (user) -> T8 -> gate 3 (user) -> T9 |

## Goal

Make the fitter's acceptance of secondary pulses rest on a per-point error that covers the coherent template
mismatch, so that fake shoulders (r < 0.1, dt < 40 ns, chi2/dof 0.40 with the shoulder vs 1.49 without) are no
longer accepted while real pile-up (Ta, r >= 0.2) is kept; measure the flat-top population (sampling phase, close
pile-up or analog saturation below 4095) and handle it with a per-module saturation amplitude; verify every step
on inspectable examples; propagate the adopted model to the rates (R_N), the hot-module examples, the benchmark
(only if the oos CB resolution is unchanged or better) and the full scans and write-up.

Non-goals: changing any benchmark default before T8 reports; the ifarm farm pass (separate track); the exotic
study (`PA/exotic_study/` untouched); the noisy-waveform flag's own rule (status recorded in T9 only); deletions.

Expected change size: PR about 150 lines (two switches, one table loader, the A_sat weight), BM about 60 lines
(the same switches in its MF copy, T8), about 8 new scripts of 150-400 lines under `PA/` (one per stage), two
constant tables (`PA/shape_err_20261116.tsv`, `PA/asat_20261116.tsv`), `PA/run_monoedge_sample.sh` +10 lines (third
variant), `PA/exotic_events/hot_modules_chi2_examples.py` +20 lines (model switch), TEX section 5 updates, one note.

## Approach

Measure the actual residual band of clean single pulses as a function of time since onset, amplitude and module
class (T1), replace the assumed shape term by the measured one behind a switch (T3a), add a template correction only
if T1 shows a coherent bias (T3b), weight out samples above a per-module saturation amplitude found in T2 (T3c), and
verify on three sample tiers before the full scans: the 50 comparison waveforms + 50 random two-pulse fits (T4), a
truth test built by overlaying real pulses (T5), and the 2-file rate sample (T6). The mono-edge rule stays as an
option, default off; the benchmark adopts the new model only on a non-regression (T8).

Rejected: keeping the mono-edge rule as the fix (it removes real Ta pile-up at r >= 0.2 and has the rise-stretch
bug when a late pulse is the largest sample); tuning `BIC_PENALTY` (a global penalty cannot separate a coherent
mismatch at the peak from a real shoulder; PR:157 records that 1-3 give the same efficiencies); simulated pulses for
the truth test (real overlaid pulses carry the real template mismatch; the pseudodata of
`pileup_deconvolution_study.py` stays as is).

## Decisions

- D1 (samples). S1 = clean single pulses from MS `old/run_<run>.npz` (channels with `ch_k == 1`, `ch_nclip == 0`,
  `ch_chi2` <= 3, 50 <= A <= 3000 ADC, decoded again from the EVIO files of `MS/splits.txt`), per run type (4) and
  module class (PbWO4 / lead glass from `mod_pwo4`), capped at 2000 waveforms per (run type, class). S2 = the 50 rows
  of `PA/monoedge_fit_comparison_20261116.tsv` (categories A-D of `PA/monoedge_fit_comparison.py:38-44`) plus 50
  random `ch_k == 2` channels of MS `old/` (seeded, 25 H2 / 25 Ta, not saturated, chi2/dof <= 3). S3 = the MS
  2-file sample. S4 = the full scans (SCAN run list, `PA/run_pileup_scan.sh`).
- D2 (error model). T3a: `V_j = sigma_ped^2 + (A s_c(tau))^2` with `tau = t_j - t0` and `s_c` a per-module-class
  table measured in T1 (rms of residual/A per tau bin, pooled over amplitude bins after checking they agree within
  20 %), replacing the `shape_sigma * sqrt(a_band/a)` term at both of its sites, the variance PR:254-261 and the
  later-candidate test of `find_candidates` PR:399-400, when `PILEUP_SHAPE_ERR=1`. The table is a
  file (`PA/shape_err_20261116.tsv`, columns `class tau_ns s`), path overridable by `PILEUP_SHAPE_ERR_FILE`. The
  diagnostic `PILEUP_SHAPE_SCALE` stays (T0 uses it).
- D3 (template bias). T3b only if T1 shows a coherent mean residual: |mean| > 0.5 rms in any tau bin at any
  amplitude bin, consistent in sign across modules (per-module spread of the mean smaller than the mean). Then the
  correction is a per-class additive shape `d_c(tau)` applied to the template when `PILEUP_TEMPLATE_CORR=1`, measured
  from the same T1 pass. Otherwise T3b is skipped and recorded.
- D4 (saturation). T3c: `PILEUP_ASAT=1` reads `PA/asat_20261116.tsv` (columns `module A_sat`), sets weight 0 for
  samples whose raw ADC is >= A_sat of that module (as the 4095 overflow does, PR:652-654 and PR:715-719) and sets a `ch_sat_analog`
  flag in the scan output; modules absent from the table are unchanged. A_sat per module comes from T2: the smallest
  raw-ADC peak value above which the flat-top width exceeds 2 samples in more than 50 % of the pulses; empty table if no
  module qualifies.
- D5 (mono-edge). `PILEUP_MONO_EDGE` stays an option, default off, untouched; T6 reports it as a third column
  for comparison only; T9 records its disposition.
- D6 (acceptance thresholds). T0: A/C -> one pulse, B/D keep two at x2 or x3. T4: A/C to one pulse >= 90 %, B/D
  keep the second >= 90 %, chi2/dof median of the 50 random fits in 0.8-1.2, the random 50 consistent (k change in
  <= 10 % of them, no new pulse at r < 0.05). T5: fake rate < 0.5 % on unoverlaid S1, efficiency > 90 % in the
  resolvable region (dt >= 16 ns, r >= 0.1; the map says where the border is). T6: dropped pulses concentrated at
  r < 0.1, dt < 40 ns; the Ta r >= 0.2 population kept within 10 %; H2 N >= 4 tail falls. T8: adopt only if the
  oos CB of control v3 (3.164 / 1.747 / 1.473 % at 0.73 / 2.24 / 3.49 GeV) is unchanged within its SE or better.
- D7 (every stage ships examples). Each stage writes `PA/<stage>_<topic>_20261116.pdf` and `.tsv` as listed in
  the packages; a stage without its PDF is not done. Example panels follow `PA/monoedge_fit_comparison.py:181`
  (waveform, model, residual strip, fit text) so the user compares across stages.
- D8 (code ownership). T3 edits only PR; BM's MF copy (BM:320, variance BM:347) gets the same switches in T8,
  behind the same env names, default off; `PA/pileup_deconvolution_study.py` is not edited (its copy stays the
  reference of the pseudodata study).
- D9 (gates). User decisions after T1/T2 (T3a vs T3b, A_sat table), after T5/T6 (run T8), after T8 (adopt for
  T9). The plan pauses at each gate; the executor reports and does not continue.
- D10 (in flight). T0, T1 and T2 were started by the supervisor on 2026-11-16 in parallel with this plan (Status
  `execute`); their outputs are
  the ones named in the packages; the executor of those stages owns their scripts and the plan takes them as given.

## Outline

- PR:120-122 switches: add `SHAPE_ERR`, `SHAPE_ERR_FILE`, `TEMPLATE_CORR` (if D3), `ASAT`, `ASAT_FILE`, same
  `os.environ.get(..., "0") == "1"` form; a loader `load_table(path) -> dict` near `Bank.__init__` (PR:169).
- PR:254-261 `MF.variance` and PR:399-400 (`find_candidates`, later-candidate test): branch on `SHAPE_ERR`:
  per-pulse term `a_k * s_class(tau)` with `s` interpolated on the tau grid of the table (linear, constant beyond
  the last bin); `SHAPE_SCALE` multiplies either term at both sites.
- Zero-weight mask `W = (Yf < overflow)` at PR:652-654 (`MF.run`) and PR:715-719 (`run_ref`, the scan's path,
  PR:770): `ASAT` sets `W = 0` for samples with raw ADC `Yf >= A_sat` of the module (same comparison as `overflow`)
  and marks the row `clipped`, so it takes the masked branch (rows without overflow run with `w = None` and would
  never see a sample mask); a new `ch_sat_analog` output key (PR:1105 `ch_*` keys). Gotcha: `bad_clipped`
  PR:629-637 (`< 6 unclipped` window samples) then counts the combined mask.
- `TEMPLATE_CORR` (only if D3 fires): additive `d_c(tau)` per class applied in `Bank.shape` PR:214-215 (not in
  the sigma branch); table `PA/template_corr_20261116.tsv`, columns `class tau_ns d`.
- PR:585-612 `add_shoulders`: unchanged; the BIC at PR:612 sees the new variance through `variance()`.
- BM:347 `variance` and BM's clip handling: port of the same branches (T8), no default change.
- `PA/run_monoedge_sample.sh`: third variant `shape` (`PILEUP_SHAPE_ERR=1 PILEUP_ASAT=1`, `PILEUP_MONO_EDGE=0`) into
  `MS/shape/`, same splits.
- `PA/exotic_events/hot_modules_chi2_examples.py:38` `_fit`: honour the env switches (they are read at PR import;
  the script sets `PILEUP_NOISE_FLAG` at line 15 the same way); output names through the existing `--stem`
  option (line 171), no new option.
- New scripts under `PA/`: `t0_shape_err_sensitivity.py` (in flight), `t1_residual_profile.py` + `t2_flat_top.py`
  (in flight, one decode pass), `T3_fixes_examples.py`, `T4_verification.py`, `T5_overlay_truth.py`,
  `T6_rates.py`, `T8_benchmark_impact.py`, `T9_scan_report.sh`. Decoding and refitting follow
  `PA/monoedge_fit_comparison.py:100-150` (`fit_both`, `make`) and `PA/pileup_monoedge.py:254-313`
  (`event_files`, `extras`); sample loading `PA/pileup_monoedge.py:70` (`load_sample`).
- Data shapes: residual profile table `class amp_bin phase_bin tau_ns mean rms sigma_assumed n spread_mean`;
  flat-top table `run_type module A width_samples ratio_amp_int integral`; overlay table
  `dt_ns r n_true n_found dt_rec_mean dt_rec_rms r_rec_mean r_rec_rms`; rates table one row per (run type, N,
  variant) with R_N and binomial error (format of `PA/pileup_noisybase.py:76` `rn_rows`).
- Invariants: all switches off reproduces the pre-T3 baseline scan (every `ch_*` key `np.array_equal`); nothing written under
  `PA/exotic_study/`; no file deleted; the user commits.

## Slices

1. T0 sensitivity. Delivers `PA/T0_shape_err_sensitivity_20261116.pdf/.tsv`: S2's 50 comparison waveforms refitted
   with `PILEUP_SHAPE_SCALE` 1, 2, 3 and the rule off; per-waveform panels x1/x2/x3. Tests:
   `t0_small: categories A and C at x2 or x3 -> one pulse in >= 80 % of the rows`;
   `t0_large: categories B and D at x3 -> two pulses kept in >= 90 % of the rows`;
   `t0_repro: x1 -> k and chi2/dof equal to k_off, chi2_off of the TSV`. Checkpoint: the TSV has 50 rows x 3
   scales and the PDF opens with a title page stating the pass/fail per category.
2. T1 residual profile and T2 flat tops (one decode pass over S1). Delivers `PA/T1_residual_profile_20261116.pdf/.tsv`
   (residual/A vs tau in 4 amplitude bins 50-3000 ADC, per class, per phase bin: mean, rms, assumed sigma,
   per-module spread of the mean; 20 example waveforms with residual and the measured band) and
   `PA/T2_flat_top_20261116.pdf/.tsv` (flat-top width vs amplitude per module, amp/integral ratio vs integral,
   A_sat table, 20 widest flat tops at high amplitude and 10 at low). Tests:
   `t1_band: tau bins at the peak and fall -> rms of residual/A larger than the assumed sigma/A (the diagnosis)`;
   `t1_decision: the TSV -> the D3 criterion evaluated and printed (bias / jitter)`;
   `t2_width: width vs A per module -> monotone rise only above a threshold, else flat (sampling-phase floor)`;
   `t2_asat: the A_sat table -> empty or every listed module has >= 20 pulses above A_sat`. Checkpoint: both PDFs
   and TSVs exist; the decision line is in the T1 title page. Gate 1: user reads T1/T2 and confirms T3a/T3b/T3c.
3. T3 fixes. Delivers the PR switches of D2-D4, the two tables, and `PA/T3_fixes_examples_20261116.pdf/.tsv` (the
   20 T1 example waveforms with the old and new error band, the 20 T2 flat tops with the A_sat weight-0 region
   marked, fit text per panel). Units: PR (one worker), `T3_fixes_examples.py` + tables (one worker). Tests:
   `off_identity: full rescan of one run of MS/splits.txt (both splits) with all switches off, once before the T3
   edits (baseline, `MS/base/run_<run>.npz`; MS/old predates the PR edits of 2026-11-16) and once after -> every ch_*
   key np.array_equal to the baseline`;
   `shape_err_band: PILEUP_SHAPE_ERR=1 on a T1 example -> V_j at the peak equal to sigma_ped^2 + (A s)^2 within 1e-6`;
   `asat_weight: PILEUP_ASAT=1 on a listed module -> samples above A_sat have w == 0 and ch_sat_analog set`;
   `asat_absent: a module not in the table -> fit identical to switches off`;
   `table_missing: PILEUP_SHAPE_ERR=1 with a wrong path -> clear error, no fallback`;
   `template_corr (only if D3 fires): PILEUP_TEMPLATE_CORR=1 on a T1 example -> template minus the off template equals
   d_c(tau) within 1e-6`. Checkpoint: the tests
   as `python3 -c` one-liners or a `PA/test_T3.py`, printed in the T3 title page.
4. T4 verification on S2. Delivers `PA/T4_verification_20261116.pdf/.tsv`: old vs new in the 50-comparison layout,
   plus the 50 random two-pulse fits. Tests: `t4_small: A and C -> one pulse >= 90 %`;
   `t4_large: B and D -> second pulse kept >= 90 %`; `t4_chi2: random 50 -> chi2/dof median in 0.8-1.2`;
   `t4_random: random 50 -> k unchanged in >= 90 %, no new pulse at r < 0.05`. Checkpoint: pass/fail per criterion
   in the TSV and the title page.
5. T5 overlay truth test. Delivers `PA/T5_overlay_truth_20261116.pdf/.tsv`: S1 pulses with a second real S1 pulse
   of the same module class overlaid at dt in {4, 8, 16, 24, 40, 80, 160, 300} ns and r in {0.02, 0.05, 0.1, 0.2,
   0.5, 1, 2} (the second pulse scaled and shifted by whole samples plus the phase offset of the template), 200
   overlays per cell, plus the unoverlaid S1 for the fake rate; efficiency map, recovered dt/r vs truth, example
   panels at every (dt, r) grid point. Tests: `t5_fake: unoverlaid S1 -> extra pulses < 0.5 %`;
   `t5_eff: resolvable region -> efficiency > 90 %`; `t5_bias: recovered dt and r in the resolvable region -> mean
   bias < 2 ns and < 10 % of r`; `t5_old: the same map with switches off -> fake rate larger (the diagnosis
   quantified)`. Checkpoint: the TSV map and the PDF.
6. T6 rates on S3. Delivers `MS/shape/` scans, `PA/T6_rates_20261116.pdf/.tsv`: R_N per run type for rule off /
   mono-edge / new, dt-r distribution of the pulses dropped by new vs off, example panels of 10 dropped and 10 kept
   pulses per run type. Tests: `t6_dropped: dropped pulses -> >= 80 % at r < 0.1 and dt < 40 ns`;
   `t6_ta: Ta r >= 0.2 population -> kept within 10 % of rule off`; `t6_tail: H2 N >= 4 -> R_4 lower than rule off`.
   Checkpoint: `PA/run_monoedge_sample.sh` third variant finished (`DONE` line in its log) and the TSV. Gate 2:
   user decides on T8.
7. T7 hot modules and in-time shapes. Delivers `PA/T7_hot_modules_20261116.pdf/.tsv`: the W121/W493 chi2 examples
   of `PA/exotic_events/hot_modules_chi2_examples.py` refitted with the new model, with an in-time residual-shape
   check (mean residual/A vs tau for the two modules against the T1 class band). Tests:
   `t7_chi2: the 40 examples -> chi2/dof median closer to 1 than with the old model`;
   `t7_shape (report only): W121/W493 in-time residual -> inside or outside the T1 band stated per module`. Checkpoint: the PDF.
8. T8 benchmark impact. Delivers the BM switches (D8), `output/energy_resolution/variants/T8/<rt>/` for the three
   run types, `PA/T8_benchmark_impact_20261116.pdf/.tsv` (oos CB, tails, eff_cut per run type, new vs V3, with SE;
   10 example waveforms per run type where the pile-up decision changed). Tests:
   `bm_off_identity: BM with switches off -> regions.tsv identical (cmp) to V3/v3_defaults/<rt>/regions.tsv and
   V3/check_identical.sh passes on its file set`;
   `bm_delta: delta rows -> oos CB new - v3 <= +1 SE at every energy (adopt), else record and stop`.
   Checkpoint: `compare_reports.py` as in `notes/plan_p14_snake_scan_20261114.md` D7. Gate 3: user decides on the
   adoption for T9.
9. T9 full scans and write-up. Delivers `output/pileup_rates/scan_shape_20261116/` (S4 with the adopted model,
   `PA/run_pileup_scan.sh` with the switches exported, no change to its default), the report of PR:1604 on it,
   `PA/T9_full_scan_20261116.pdf/.tsv` (final R_N table and 10 examples per run type), TEX section 5.2 (line 611)
   and the constants table (line 1355) updated, and `notes/pileup_fit_validation_20261116.md` recording the
   mono-edge disposition and the noise-flag status. Tests: `t9_rn: final R_N table -> same ordering of run types
   as T6 and the H2 N >= 4 tail lower than the 2026-11-15 table`; `t9_tex: TEX -> compiles once TeX is available
   (not on this machine; the user compiles), the new switch rows present in the constants table`. Checkpoint:
   the report PDF and the note.

## Verification

After slice 9: with every switch off the scan reproduces the pre-T3 baseline (`off_identity`); with the adopted switches the
S2 criteria of D6 hold, the T5 fake rate and efficiency hold, the S3 R_N table shows the fake-shoulder drop with the
Ta population kept, the benchmark is not worse than V3, and the final R_N table, TEX and the note agree.

## Risks

- T1's band may depend on amplitude in a way that a single `s_c(tau)` does not capture (the current model scales
  with `sqrt(a_band/a)`); the amplitude-bin check of D2 says so, and then the table gets an amplitude axis.
- Lead glass uses the bank's average PbWO4 spline (PR:15); its band will be wider and its S1 smaller; the class
  split keeps that from widening the PbWO4 error.
- A larger error near the peak lowers the sensitivity to real close pile-up (dt < 16 ns); T5 maps the loss, and
  T6 checks the Ta population.
- Overlaying two real pulses doubles the pedestal noise and the template mismatch; the T5 fake rate on unoverlaid
  pulses and the T5 efficiency at large dt bound this.
- The benchmark's MF copy differs from PR beyond the variance (its own clip handling); the T8 port is verified by
  the off-identity test only, so a divergence elsewhere stays.
- The scan output grows by one key (`ch_sat_analog`); `PA/pileup_monoedge.py` and `PA/pileup_noisybase.py` read
  known keys and are unaffected.

## Scope Notes

In: T0-T9 as written, with examples at every stage. Out: the ifarm farm pass (separate track); any benchmark
default change before T8 reports and the user decides; deletions anywhere; anything under `PA/exotic_study/`;
the noisy-flag rule itself; `PA/pileup_deconvolution_study.py`; git (the user commits). User decisions the executor
must not reopen: the diagnosis (per-point error underestimates the coherent template mismatch), the fix family
(measured shape error behind a switch, mono-edge stays an option off), the sample tiers S1-S4, the three gates,
the T8 adoption rule. Memory holds no entries on the pile-up fitter (grep of `.supervisor/memory/` 2026-11-16).

## Assumptions

- The MS old/new scans and `MS/splits.txt` are complete for the 13 runs; checked by T0's loader (26 npz).
- The EVIO splits of `MS/splits.txt` are still on disk; checked by T1's decode step before fitting.
- `mod_pwo4` of the scan npz is the module class; checked against the names (W = PbWO4, G = lead glass) in T1.
- `PILEUP_SHAPE_SCALE` scales exactly the shape term at its two sites (PR:261 and PR:399); checked by `t0_repro`.
- `V3/v3_defaults/<rt>/` is the reference for `bm_off_identity`: `regions.tsv` by `cmp`, plus the file set of
  `V3/check_identical.sh` (`summary.tsv splits.tsv regions_splits.tsv gains.tsv poscorr_*.tsv runs.tsv`).

## Work Packages

### T0: shape-error sensitivity
- Goal: S2 comparison waveforms refitted at shape scale x1/x2/x3; done when slice 1 tests are evaluated and the
  per-waveform panels exist
- Owns: `PA/t0_shape_err_sensitivity.py`, `PA/T0_shape_err_sensitivity_20261116.pdf/.tsv`
- Deps: none
- Stages: execute
- Status: execute
- Artifacts: `PA/T0_shape_err_sensitivity_20261116.pdf`
- Cost: 10 min

### T1: residual profile of clean pulses
- Goal: measured residual/A band vs tau per class, amplitude and phase bin, with the bias/jitter decision; done
  when slice 2's t1 tests are evaluated and the 20 examples exist
- Owns: `PA/t1_residual_profile.py`, `PA/T1_residual_profile_20261116.pdf/.tsv`
- Deps: none
- Stages: execute
- Status: execute
- Artifacts: `PA/T1_residual_profile_20261116.pdf`
- Cost: 30 min

### T2: flat tops
- Goal: flat-top width vs amplitude per module, amp/integral cross-check, A_sat table; done when slice 2's t2
  tests are evaluated and the 30 examples exist
- Owns: `PA/t2_flat_top.py`, `PA/T2_flat_top_20261116.pdf/.tsv`
- Deps: none (shares T1's decode pass)
- Stages: execute
- Status: execute
- Artifacts: `PA/T2_flat_top_20261116.pdf`
- Cost: 20 min

### T3: fixes behind switches
- Goal: `PILEUP_SHAPE_ERR`, `PILEUP_ASAT` (and `PILEUP_TEMPLATE_CORR` only if D3 fires) in PR, default off, with the
  two tables; done when slice 3 tests pass and the examples PDF exists
- Owns: PR, `PA/shape_err_20261116.tsv`, `PA/asat_20261116.tsv`, `PA/T3_fixes_examples.py`,
  `PA/T3_fixes_examples_20261116.pdf/.tsv`, `PA/test_T3.py`
- Deps: T0, T1, T2; gate 1 (user)
- Stages: execute
- Status: blocked: execute: gate 1 (T1/T2 results and the user's T3a/T3b/A_sat decision)
- Artifacts: `PA/T3_fixes_examples_20261116.pdf`
- Cost: 1-2 h

### T4: verification on S2
- Goal: old vs new on the 50 comparison waveforms and 50 random two-pulse fits; done when the D6 T4 criteria are
  evaluated in the TSV and the PDF exists
- Owns: `PA/T4_verification.py`, `PA/T4_verification_20261116.pdf/.tsv`
- Deps: T3
- Stages: execute
- Status: blocked: execute: T3
- Artifacts: `PA/T4_verification_20261116.pdf`
- Cost: 15 min

### T5: overlay truth test
- Goal: efficiency and fake-rate map on (dt, r) from real overlaid pulses; done when the D6 T5 criteria are
  evaluated and the grid examples exist
- Owns: `PA/T5_overlay_truth.py`, `PA/T5_overlay_truth_20261116.pdf/.tsv`
- Deps: T3
- Stages: execute
- Status: blocked: execute: T3
- Artifacts: `PA/T5_overlay_truth_20261116.pdf`
- Cost: 1 h

### T6: rates on S3
- Goal: R_N table rule off / mono-edge / new and the dropped-pulse dt-r distribution; done when the D6 T6
  criteria are evaluated and the examples exist
- Owns: `PA/run_monoedge_sample.sh` (third variant), `MS/shape/`, `PA/T6_rates.py`, `PA/T6_rates_20261116.pdf/.tsv`
- Deps: T4, T5
- Stages: execute
- Status: blocked: execute: T4, T5
- Artifacts: `PA/T6_rates_20261116.pdf`
- Cost: 1 h (about 20 min wall of scans at 16 lanes)

### T7: hot modules and in-time shapes
- Goal: W121/W493 chi2 examples with the new model and the in-time residual-shape check; done when slice 7
  tests are evaluated
- Owns: `PA/exotic_events/hot_modules_chi2_examples.py` (switch only; `--stem` names the output), `PA/T7_hot_modules_20261116.pdf/.tsv`
- Deps: T4
- Stages: execute
- Status: blocked: execute: T4
- Artifacts: `PA/T7_hot_modules_20261116.pdf`
- Cost: 30 min

### T8: benchmark impact
- Goal: BM switches and the new-vs-V3 comparison at the three energies with the adoption verdict; done when
  `bm_off_identity` passes and the delta TSV and PDF exist
- Owns: BM (switch port only), `output/energy_resolution/variants/T8/`, `PA/T8_benchmark_impact.py`,
  `PA/T8_benchmark_impact_20261116.pdf/.tsv`
- Deps: T6; gate 2 (user)
- Stages: execute
- Status: blocked: execute: gate 2 (T5/T6 results)
- Artifacts: `PA/T8_benchmark_impact_20261116.pdf`
- Cost: 2 h of reports

### T9: full scans and write-up
- Goal: S4 with the adopted model, final R_N table, TEX update and the note; done when the scan report, the PDF,
  the TEX rows and the note exist
- Owns: `PA/T9_scan_report.sh`, `output/pileup_rates/scan_shape_20261116/`, `PA/T9_full_scan_20261116.pdf/.tsv`,
  TEX section 5.2 and the constants table, `notes/pileup_fit_validation_20261116.md`
- Deps: T8; gate 3 (user)
- Stages: execute, write
- Status: blocked: execute: gate 3 (T8 verdict)
- Artifacts: `PA/T9_full_scan_20261116.pdf`, `notes/pileup_fit_validation_20261116.md`
- Cost: overnight scan

## Log

- 2026-11-16 decision: plan approved (user); T0, T1, T2 started in parallel with the plan; source: user
- 2026-11-16 T1x (extension of T1): chi2 fraction vs amplitude per run type from the full scans, and high-amplitude residual profiles up to 4095 ADC (pileup_analysis/T1x_chi2_vs_amplitude_20261116.pdf, T1x_highamp_residual_20261116.pdf). Curves overlap below 2000 ADC (2p24 H2 / 3p49 H2 ratio 0.86–1.5; Ta 10–50× higher from pile-up). Above 2500 ADC 2p24 H2 rises sharply (27 % chi2>20 at 2500–3000, ~97 % above 3000, from 960 and 240 rows): a rare population above the elastic edge with a second hump (+4 % at peak, −7.5 % at fall), not a bank effect. For normal pulses the coherent mean does not grow with A (3p49 at 3000–4095: −0.03 % peak, −0.23 % fall); deferred: none.
- 2026-11-17 side check (user request, not a package): full-energy single-pulse channels with chi2/dof > 20 (pileup_analysis/single_peak_beam_energy_chi2gt20_20261116.pdf ±10 %, _3sigma_ ±3σ, _inwindow_by_k_ table). Full-energy single-module channels at 2.24/3.49 GeV are a tiny set (~200/100 of 20–45 M) dominated by multi-pulse and saturated waveforms at or above the elastic edge; the chi2 fraction there measures the selection, not the fitter. Unexplained: 21k such channels at 0.73 GeV vs 200 at 2.24 (energy-scale/window check pending).
- 2026-11-17 side check (user request, item 16): fitted integral A·I_m vs raw walk-out sum for monotonic single peaks (pileup_analysis/integral_vs_rawsum_20261117.pdf; notes/note_integral_vs_rawsum_20261117.md). Median ratio 0.994–0.996 (rms 0.011–0.013) at chi2/dof < 5 and 0.998–1.003 at chi2/dof > 20 with complete windows; amplitude drift 0.5 % over 50–3000 ADC. Truncated late pulses (ratio 1.25–1.28) are the fit extrapolating past the record end. Closed by the user 2026-11-17.
- 2026-11-17 decision: gate-1 inputs: T3b (template correction) and T3c (analog saturation) not indicated by T1/T1x/T2; T3a (per-module shape error, PILEUP_SHAPE_SCALE=2 fallback) remains the candidate; user decision on starting T3a pending.
- 2026-11-17 decision: late pulses truncated by the record end get no special treatment (time gate removes them from the energy measurement; no overlap with the in-time template window). At T6 the rate table splits pulses into an overlapping column (onset within about [t_trig − 40 ns, t_trig + 130 ns]) and an occupancy column (all others).
- 2026-11-17 decision: user feedback on the 11/16 deliverables: noisy-baseline flag (removed from pileup_rates.py, patch kept), monotonic-edge rule and its checks (PILEUP_MONO_EDGE stays an off-by-default option), and the flat-top study (T2) "doesn't really work as I expected. I will think more." Parked for the user's rethink; no code change. Left over: ifarm access (ssh -R login to ifarm2402 fails; exact error pending). Item 16 note: notes/note_integral_vs_rawsum_20261117.md.
- 2026-11-17 side track (user request): PRADLOG summary written (notes/elog_pradlog_summary_20261117.md) via logbook API with the user's ~/.netrc credentials; not obtained: attachments (HV/gain spreadsheets), run tables cut by the 4000-char limit (e.g. runs 25995–26048), two special-test reports outside PRADLOG (4515733, 4515793), file availability. Facts bearing on this plan: trigger latency/lookback change 07-03 [4519495] (check the 140–156 ns window per X17 run); inner-module gain loss 50–70 % by 08-08 [4539799] and inner-layer HV raised 08-10 / partly reverted 08-14 [4540517, 4543254], run 26230 inside that window; local run 24578 is a random-trigger run at 0.7 GeV [4498531], not H2 production.
- 2026-11-17 decision (user): "(1) Can you fetch RCDB? (2) go (3) Yes (4) Yes" → fetch RCDB run types; non-production runs: 24578 first (local), then one beam-off and 26266 cosmics; re-examine 26230 vs the HV change; pull the two special-test reports. X17 Møller third-cluster rule: "do it and keep 300 MeV as the comparison" (isolation cut becomes default).
- 2026-11-17 side track (user request): integral vs raw sum over all data: PILEUP_RAWSUM=1 added to pileup_rates.py (off-identity verified), 2-file sample (pileup_analysis/integral_vs_rawsum_full_20261117.pdf) and full 13-run scan (integral_vs_rawsum_full_all_20261117.pdf, output/pileup_rates/rawsum_full). PbWO4 monotonic single pulses, A ≥ 50: median A·I_m/raw 0.9945–0.9964, 16–84 % 0.984–1.009, phase saw-tooth 1.2–1.3 % p-p, per-module sd 0.003–0.005; fraction within ±2 % falls with chi2/dof (93 % at <1, 60 % at 5–20; above 20: Ta 85–93 %, H2 25–55 %). All-cases k≥2 share 4.8/7.1/8.1/32.6 %.
- 2026-11-17 side track (user request): run 24578 random-trigger study (pileup_analysis/random_trigger_24578_20261117.pdf). Readout is zero-suppressed (97.6 % of random events have no HyCal channel), so no pulse-free records; noise/pedestal consistent with production (rms/σ_ref 0.95, b_ref diff 0.00 ADC, no drift); random coincidence predicts ≤ ~10 % of the measured k≥2 at A≥20 in 0.73 production (0.024 % vs 0.262 %); second pulses clustered (17× Poisson). Verdict: not a pile-up-free reference; usable for rates. data/evio/README.md should relabel 24578 as random-trigger.
- 2026-11-17 side track (user request): X17 energy resolution with Møller two-cluster calibration (moller_calib.py, energy_resolution/2p24GeV_Ta_X17_moller/): pair-sum σ/E 2.09 ± 0.04 % at 2.24 GeV (isolation cut 60 mm, default; 300 MeV rule 2.11 ± 0.05 %), symmetric pairs 2.85 ± 0.11 % at 1.12 GeV; ep value 1.987 %; all 15–25 % above the H2 curve. HV check (energy_resolution/2p24GeV_Ta_X17_hvcheck/): local Ta runs 25800, 26138 (before the 08-10 HV raise) and 26230 (inside); HV step visible in calibration files and raw ADC/MeV (2.2×), not in resolution (1.984/1.942/2.003 %); pre-window pooled 1.98 ± 0.03 %. User decision pending on which to quote.
- 2026-11-17 side track (user request): exotic events plan P1–P3 (pileup_analysis/exotic_events/). P1: quantities Q1–Q14 on 12 user examples (types 1–4) + random samples; cuts: f1 Δb>8σ & tail<−5σ, f2 w50 ratio>3, f4 nclip≥1 & (dclip≤8 ns or E_sec>0.5 E_beam). P2: flags behind PILEUP_EXOTIC=1, options PILEUP_BASE_SLOPE=1 (gated), PILEUP_CLIP_SINGLE=1 (8 ns), identity 145/145 keys; stretched-template flag f3 catches types 1 and 2 (one smooth-wide population), not type 3 (indistinguishable from real close pile-up by shape).
- 2026-11-17 decision (user): "Keep flag 3 as the smooth-wide flag for types 1–2 and give me 20 examples each for the identified events; 2 and 3 yes for now but I'll need to think more after seeing more examples flagged" → BASE_SLOPE and CLIP_SINGLE (8 ns) adopted provisionally; P3 (2-file sample with flags, rates, 20 examples per flag per run type) dispatched.
- 2026-11-17 lead: P3 done (sample scan 2.5 h wall; exotic_events/exotic_P3_sample_20261117.{pdf,tsv,_modules.tsv,_examples.tsv}); any-flag rate 0.08–0.15 % H2, 1.08 % Ta X17; f4 leakage into ±3σ elastic window 14–15 % of f4 channels (H2); deferred: CLIP_SINGLE k=0 rule, per-type Ta X17 thresholds, user's decision on flags 2 and 3 after seeing examples
- 2026-11-17 decision (user): P3 feedback — per-run-type flag thresholds: yes; marginal f1 cases (db 8–13σ) stay flagged; exclude window-truncated pulses from the w50/f2 flag; f4 leakage stays per module (no per-cluster); BASE_SLOPE and CLIP_SINGLE stay off, flags are event markers only; items 16–22 inspected
- 2026-11-17 supervisor: P4 dispatched — full scan all run types with PILEUP_EXOTIC=1 (output/pileup_rates/exotic_full/), per-module rate per flag, four module-map pages in exotic_events/exotic_P4_module_rates_20261117.pdf
- 2026-11-17 decision (user): P4 blocked by disk (926 MiB free, scan needs ~11 GB); user removed output/pileup_rates/rawsum_full npz/logs themselves (PDF/TSV outputs kept); P4 restarted with per-file split, Ta files first, 16 workers (est. 37 CPU-h, ~3 h wall)
- 2026-11-18 lead: P4 done. Full exotic scan (13 runs, 103 files, PILEUP_EXOTIC=1 only) took 6.0 h wall at 16 procs (339,000 CPU-s). Per-module rates per flag: any-flag 0.069 / 0.089 / 0.107 % (0.73 / 2.24 / 3.49 H2) and 0.84 % Ta X17 (1.11 % uniform). Per-type thresholds (same-sigma on the k=1 core, H2-pooled reference) move little (only Ta db 8 → 9.76); the Ta excess persists. f2 truncation exclusion (t0 > 284 ns) removed 7–23 % of H2 f2 flags, 0.6 % of Ta. Hot modules: W335 at 0.73 GeV, W613 for f4, W1/W950 for f2/f3. Ta X17 f1/f2/f3 concentrate in the inner ring around the beam hole; H2 f3/f4 rise toward the outer edge. exotic_events/exotic_P4_module_rates_20261117.{pdf,tsv,_thresholds.tsv,_summary.tsv}; deferred: scan npz kept only for 19 Ta X17 files, the other 85 reduced to candidate-row pickles (disk), so any quantity beyond the stored candidate rows needs a rescan
- 2026-11-18 decision (user): P4 feedback — no module masking (per-channel flags only); npz reduction to candidate-row pickles accepted, no rescan; uniform (P2) thresholds remain the reference, per-type set kept in the appendix (user asked what they are; no preference stated)
- 2026-11-18 supervisor: P5 dispatched — maximum amplitude per exotic type (raw peak and largest fitted A) per run type, median and 16/84 % band, exotic_events/exotic_P5_max_amplitude_20261118.{pdf,tsv}
- 2026-11-18 lead: P5 done (exotic_events/exotic_P5_max_amplitude_20261118.{pdf,tsv}, exotic_p5_max_amplitude.py). Max amplitude of flagged channels per run type, raw peak (max sample − pedestal, EVIO re-read, decode only) and fitted A_max (= A_main), median [16 %, 84 %]: f4 raw ≈ 3930 ADC, all saturated, fitted 3600–10000 ADC (≈ 1–2 E_beam, 14 at 0.73 GeV with N=2); f1 raw 59–81 ADC (fitted 71–91; 2–4 % E_beam, 12 % at 0.73 GeV); f2 ≈ 10 raw / 33 fitted ADC (≈ 1 % E_beam) for Ta X17 and 3.49 GeV H2, bimodal with a saturated spike (38.6 %) at 2.24 GeV H2; f3 bimodal in H2 (≈ 10 ADC and ≈ 400 ADC), 16 raw / 32 fitted ADC in Ta X17; all kept ≈ 28 ADC. Sources: H2 from the 2-file sample npz (the H2 pickles carry no amplitudes), Ta X17 from the 19 full npz; 1 % random reference sample; deferred: H2 f2/f4 small-N (0.73 GeV f4 N=2, f2 N=56), larger H2 N needs an amplitude-carrying rescan; critic review skipped by the lead (self-validated)
- 2026-11-18 decision (user): add "clipped ⇒ not f2"; evaluate f4 first, then f2, so saturated records are never double counted; no module masking; npz reduction accepted
- 2026-11-18 lead: f2 now excludes clipped records (nclip ≥ 1) and window-truncated ones and is evaluated after f4, through the shared `f2_eligible` in exotic_p4_module_rates.py used by P4 and P5; clip rule removes f2 flags 0/13/9/0 (uniform; 0.73/2.24/3.49 H2/Ta X17), new f2 rates 0.0012/0.0006/0.0031/0.1108 %, f2∩f4 = 0, f1∩f4 = 0, f3∩f4 = 6/192/42/4 (unchanged); P4 and P5 outputs regenerated in place; P5 f2 at 2.24 GeV H2 still bimodal (q84 ≈ 900 ADC, unclipped rows, 86 saturated rows removed); deferred: f3∩f4 priority, clipped channels in the per-type threshold core; critic skipped by the lead
- 2026-11-18 decision (user): new flag for saturated pulses (f5 = nclip ≥ 1 and not f4, order f4 → f5 → f2); amplitude vs w50 for f2; saturated ratio per run type and per module from the full scan, quoted only where statistics allow, with intervals
- 2026-11-18 lead: f5 added to P4/P5 with saturated fractions per run type and module (Wilson 68 %, quoted at n_kept ≥ 100 and n_sat ≥ 5), radial and edge-ADC tables; new exotic_events/exotic_P6_amp_vs_w50_20261118.{pdf,tsv} (amp-vs-w50 maps, f2 high-A table, 40 examples); saturation 0.00015/0.0038/0.0015/0.00003 % of kept channels (0.73/2.24/3.49 H2/Ta), 0.10 % of E>0.5 E_beam pulses at 2.24 H2, clustered in W613/W307/W335 and rising with radius, not at the inner ring; the 2.24 H2 high-A f2 mode is merged k≥4 pile-up over 18 modules, the 0.73 H2 one mostly W335; w50 = 3 sits in a density gap for H2 but is an arbitrary cut in a continuum for Ta; critic round PASS (4 fixed, 1 rejected); deferred: raw peak for H2 full scan (pickles lack event ids), count criterion on f1–f5 colour maps, run 24937 (3.49 sample) looks Ta-like — check run-type assignment
- 2026-11-18 decision (user): saturated pulses concluded at "up to 0.1 % (of full-energy pulses) for now"; statistics too small for a solid conclusion; keep flag f5 for the full-data scan on ifarm
- 2026-11-18 decision (user): flags f1–f5 become mutually exclusive categories (no double counting; priority f4 → f5 → f1 → f2 → f3, supervisor's order, pending user's objection), so no separate f4-before-f3 rule; w50 > 3 (f2) examples to be added to exotic_P6; run 24937 is "not a good run" — dropped from the exotic study inputs; template FWHM 23.5 ns (PbWO4, all run types), so w50 > 3 ≈ data FWHM > 70 ns; Ta w50 cut kept at 3 (tail fraction, not a distinct population)
- 2026-11-18 decision (user): after P6, one consolidated PDF with the latest numbers and all plots from the previous exotic steps; then proceed with P7–P10 (cluster-level impact packages)
- 2026-11-18 lead: exclusive categories f4 > f5 > f1 > f2 > f3 implemented in the shared flags_of (raw plus exclusive rates and a moved matrix in P4/P5/P6); run 24937 is 3.49 GeV H2 (not 2.24 as the P4 enumeration said) and is excluded (full scan 103 → 99 files; 3.49 f2 1446 → 269 uniform, 2.24 unchanged 562); exclusive any-flag 0.069/0.090/0.073/1.114 % (0.73/2.24/3.49 H2/Ta); Ta f2 falls 0.111 → 0.016 % and f3 0.392 → 0.198 % once f1 takes precedence; 20 uniform w50 > 3 f2 examples per run type added to P6 (Ta f2 looks like baseline-wobble over-fit, H2 f2 a mix of real k 4–6 trains and 10–50 ADC noisy trains); critic PASS; deferred: P1–P3 PDFs still include 24937, 3.49 GeV sample now one run (24911) — thin
- 2026-11-18 lead: P6c done — exotic_events/exotic_summary_20261118.pdf (119 pp, exotic_summary.py): 8 front-matter pages read from the P4/P5/P6 TSVs, then P1–P6 with section pages; P1–P3 unchanged, marked pre-exclusive incl. 24937; critic PASS (6 SHOULD, 3 NIT applied); earlier Log numbers superseded by the TSVs: 3.49 f2 rate 0.0031 → 0.00058 % (24937 dropped), f3∩f4 3.49 42 → 41, truncation removal 23–25 % of H2 f2 (not 7–23); deferred: visual page check, amplitude-carrying H2 rescan, 3.49 sample is one run
- 2026-11-18 lead: P7 done (exotic_events/exotic_P7_cluster_level_20261118.{pdf,tsv,_examples.tsv}, exotic_p7_cluster_level.py). Join: exotic rows → the benchmark's own gated clusters (control v3, same EVIO splits, 100 % event match, Møller pairs reproduced exactly). Clusters with ≥1 flagged module: 0.13/0.30/0.31 % (0.73/2.24/3.49 H2, sample) and 4.3 % (Ta full); elastic 0.10/0.15/0.20/6.0 %. Flagged seeds 0.6/2.1/1.8/22 per 1000 events. Saturated seeds in the elastic window 0/0.0044/0.0086/0 %, single-module fakes (several on W457). No H2 calibration cell reaches 0.1 % mean flagged-energy fraction; Ta elastic f1-neighbour 0.16 % and f3-seed 0.13 % do. 0.73 GeV elastic window restricted to angles where ep and Møller separate (R ≳ 107 mm). Sample run tags: 3.49 GeV H2 = 24698, 24911 (so the 3.49 sample is two runs, not one); 2.24 GeV H2 = 25228–25320. Critic PASS round 3. Cache output/pileup_rates/exotic_p7_cache.npz (312 MB, needed by P8). Deferred: per-module ΔE (P8), verdict on single-module saturated fakes (P10), 191k of 558k Ta f1 rows are outside any cluster
- 2026-11-18 supervisor: P8 dispatched — energy error per category, Ta elastic f1-neighbour and f3-seed cells first, per module; 0.73 GeV elastic sample restricted to R ≳ 107 mm (supervisor's choice, pending user)
- 2026-11-18 decision (user): P7 feedback — single-module saturated clusters at ≈E_beam are treated as fakes (high-gain module + non-shower single-channel pulse; to veto in P10); 0.73 GeV calibration sample restricted to R ≳ 107 mm for P8/P9; Ta f1 rows outside any cluster left to P10; the one-run note on the 3.49 sample in earlier lines is superseded (24698 and 24911 are both 3.49 GeV H2)
- 2026-11-18 decision (user): "By the end of today, we would like to summarize the exotic event study based on the findings and dicussions as one chapter in the write up tex file. Make a new folder under writeup/figs to include the new figures for the writeup and update the tex. This will be done after P10. Make sure the report is concise and logical."; gain check: W457/W564 gains are normal (ratio to median 1.00–1.03), W613/W307/W335 outside the r < 200 mm calibration pool; saturated seeds inside the elastic window are therefore narrow non-template spikes (fitted A ≪ raw peak), not high-gain modules
- 2026-11-18 decision (user): saturated seeds in the elastic window are two populations — (a) fakes: 1–2-member clusters, k 2–6, chi2/dof 10–1000 (W457, W305, W556 at 2.24; W299, W129, W44, W387 at 3.49), fitted A ≈ 0.7–1.0 of the clip level; (b) real elastic showers with a clipped high-gain seed: 15–19 members, seed fraction 0.77–0.81, k=1, chi2/dof 1–2 (W398 at 3.49 GeV, ×3) — not to be vetoed, an energy-estimate issue for P8/P9. P10 veto V3 = saturated seed AND ≤ 2 members (applies to (a) only); (b) reported per module as high-gain candidates; consolidated summary PDF to be updated with P7 and this split (executor dispatched)
- 2026-11-18 lead: summary PDF updated (exotic_summary_20261118.pdf, 164 pp = 11 front + P1–P7): page 9 P7 cluster level, pages 10–11 saturated seeds in the elastic window (41 clusters; new exotic_p7_saturated_seeds.py + exotic_P7_saturated_seeds_20261118.tsv); (a) ≤ 2 members: 23 at 2.24 (W457 16), 9 at 3.49; (b) ≥ 3 members: 9 at 3.49, of which only W398 ×4 (not ×3) are real showers (15–19 members, seed fraction 0.77–0.81, k=1, chi2/dof 1–2) — the other 5 (W613 ×2, W349, W636, W762; seed fraction ≥ 0.985, chi2/dof 43–872) are fakes with tiny neighbours; critic PASS. Supervisor: V3 veto redefined for P10 as saturated seed AND seed energy fraction > 0.95 (catches all 37 fakes, spares the 4 W398 showers); W398 = high-gain candidate
- 2026-11-18 lead: P10 done (exotic_events/exotic_P10_yield_impact_20261118.{pdf,tsv}, exotic_p10_yield_impact.py; cache output/pileup_rates/exotic_p10_cache.pkl 52 MB). Veto loss in the ep elastic window, V3 = saturated seed with seed fraction > 0.95: 0 / 0.0044 / 0.0067 / 0 % (0.73 / 2.24 / 3.49 H2 / Ta), removing 23 and 14 fakes at 2.24 and 3.49 GeV; the 4 W398 clipped-seed showers kept. V4 (+ f1/f3 seeds) costs 0.036 / 0.008 / 0.007 / 0.286 % (Møller pairs 0.122 %, Ta inner ring 0.305 %); V2 (any member) 0.1–0.2 % on H2 and 6.0 % on Ta (1.6 % of Møller pairs), not viable. Ta inner-ring f1 rate follows run-level occupancy (0.034 ± 0.002 %/kept channel, elasticity 1.06, 3 runs), ×1.21 for current ×1.17 (25800 3.0 nA from the elog, RCDB blank); current, occupancy, trigger version and inner-layer HV confounded. Out-of-cluster f1 (Ta): 366k of 557.5k are outside a cluster (the P7 line had the counts reversed), 97.5 % with no own-pulse hit, mostly early pulses ~100 ns before the window; only 120 reach 500 MeV; no action for yields. In-cluster f1: 114.6k of 191.2k are flagged on a pulse other than the clustered hit (input to P8). Critic PASS round 2 (BLOCKING own-pulse match fixed). Recommendation: V3 on 2.24/3.49 H2 (cost < 0.01 %, removes fakes only); V4 only if P9 shows a gain; V2 rejected
- 2026-11-18 decision (user): P10 feedback — V3 (saturated seed with seed fraction > 0.95) adopted as policy on 2.24/3.49 GeV H2; V4 decided after P9; run 25800 beam current taken as 3.0 nA (elog; "looks like yes"); W398 listed as the high-gain module for P8/P9; P11 write-up chapter drafting in parallel (writer dispatched) with the P8/P9 subsection as a placeholder
- 2026-11-18 lead: P8 done (exotic_events/exotic_P8_energy_error_20261118.{pdf,tsv,_examples.tsv}, exotic_p8_energy_error.py; cache output/pileup_rates/exotic_p8_cache.npz, needed by P9). 391,292 flagged cluster-member records refitted (default reproduced 100 %, P7 join 100 %); ΔE on the clustered hit's own pulse. Energy-scale shift: Ta elastic f1 sloped-baseline −0.0144 ± 0.0003 % (inner ring −0.0153 %; hit = main −0.0076, hit ≠ main −0.0068), f3 stretched −0.0028 %; Ta Møller f1 −0.0051 % (inner ring −0.0114 %); Ta inner-ring modules −0.018 to −0.027 % on their seeded elastic clusters; every credible H2 cell < 0.01 %, none ≥ 0.1 %. Better estimates: f1 sloped baseline (credible), f3 stretched in Ta only, default for f4, f2 and H2 f3, f5 k=1 extrapolation for mild clipping (W398: −0.7 to −1.0 % of E_cluster by the integral), deep saturation not correctable → veto. Critic PASS round 2 (double credit across split clusters fixed). Deferred: second-order propagation (reclustering, window membership) to P9; f5 deep-saturation veto to P9/P10; H2 small-N cells
- 2026-11-18 supervisor: P9 dispatched — benchmark three ways (default / V3 veto, V4 variant / corrected: f1 sloped baseline everywhere, f3 stretched in Ta only, f5 k=1 for mild clipping), per-module Δgain, elastic-peak shift, σ/E, tails
- 2026-11-18 writer: P11 draft — write-up chapter "Exotic waveforms" (sec:exotic) added after Section 5 in writeup/waveform_analysis_prad2_x17.tex with 4 figures in writeup/figs/exotic_events/ from writeup/make_figs/exotic_events_figs.py (types 2×2, module maps, amp vs w50, saturated seeds); 6 subsections (types/flags, rates, nature, cluster level, yields with V3 adopted on 2.24/3.49 H2, P8/P9 placeholder + status); 5.2 and Summary pointers updated; static TeX checks pass (no engine on this machine), ≈5.9 pp; critic round 2 PASS; P8 numbers being slotted in; deferred: P9 paragraph, compile check and page count by the user, Summary-table row for the chapter
- 2026-11-18 writer: P11 P8 filled — 6.6 now carries the credible alternatives per category and the energy-scale shifts (Ta elastic f1 −0.0144 ± 0.0003 %, f3 −0.0028 %; Møller f1 −0.0051 %; inner-ring modules −0.018 to −0.027 %; credible H2 < 0.01 %) with a fifth figure exotic_energy_scale.pdf; P9 paragraph remains a TODO placeholder; critic round 3 fixes applied (unsourced "well inside the calibration uncertainty" removed pending a sourced value); est. 6–6.5 pp, not compiled
- 2026-11-18 lead: P9 done (exotic_events/exotic_P9_calibration_impact_20261118.{pdf,tsv}, exotic_p9_calibration_impact.py; opt-in EXOTIC_HOOK added to energy_resolution_benchmark.py and moller_calib.py, default path reproduces v3 exactly; intermediates output/pileup_rates/exotic_p9, 368 MB, ledgered needed-by P11). Arms A / V3 / V4 / C (first-order P8 credible ΔE; ambiguous keys kept default) refit with the benchmark chain on the P7/P8 splits (0.73 GeV at R ≥ 107 mm): no arm moves oos σ/E, the peak or tails beyond 2σ of the paired bootstrap error at any run type or the Ta Møller sum (largest |z| 1.8); no free module gain moves by 2σ_g (largest W659, Ta Møller V4, −0.87 %, −1.98σ); the C scale shift is absorbed by the gains (Ta peak −0.0004 ± 0.0047 %); the Ta inner ring (block H) is prior-held in both Ta calibrations, so its Δgain is not measurable. Critic PASS round 2 (BLOCKING key collision fixed). Deferred: P8 ρ = −1 credibility (81k Ta members), f5 mild threshold (nclip ≤ 3) to confirm, a calibration that frees the H block
- 2026-11-18 gate (P9/P10): 2.24/3.49 GeV H2 keep V3 (0.0044/0.0067 % of the elastic window, fakes only; user policy), no V4, no in-place correction; 0.73 GeV and Ta X17: no veto, no correction; flags remain event markers; f5 kept for the farm scan; writer filling the P9 paragraph (P11)
- 2026-11-18 writer: P11 done — exotic chapter complete in writeup/waveform_analysis_prad2_x17.tex (sec:exotic, 6 subsections, 4 tables, 5 figures in writeup/figs/exotic_events via make_figs/exotic_events_figs.py); P8 and P9 filled from the exotic_P8/P9 TSVs; gate decision recorded; 5.2 pointer, script map, roadmap, tab:scripts and Summary open point updated; critic 4 rounds, final PASS; statics clean, not compiled (no TeX engine), est. 6.3–6.5 pages; deferred: compile check of the page count by the user; next cuts if needed: the f5 method clause in P8 and the second "f2: wide pulses" sentence block
- 2026-11-18 lead: summary PDF extended (exotic_events/exotic_summary_20261118.pdf, 221 pp = 14 front + P1–P7 + P8, P10, P9; exotic_summary.py): new pages 12 P8 energy error, 13 P10 yield impact, 14 P9 calibration impact and gate (numbers read from the three TSVs); page 1 map and page 8 decisions/gate/deferred updated; critic PASS round 1, 6 findings fixed. TSV vs Log corrections (applied to the TeX by fixer): the credible H2 cell maximum is +0.012 ± 0.014 % (f5 integral, 3.49 GeV, 3 flagged), cells with ≥ 10 flagged ≤ 0.001 %; Ta block H is prior-held in the Møller calibration and largely in ep (6/36 free, prior 0.47–0.57); max |z| 1.77 oos (1.95 with rings); the P10 totals 557,506 / 191,223 / 114,637 / 120 and the W398 integral figure come from the PDFs, not the TSVs; deferred: P8 ρ = −1 for the 81k Ta members, f5 mild threshold nclip ≤ 3, a calibration that frees block H, farm scan, H2 small-N cells
- 2026-11-18 decision (user): exotic study conclusion confirmed — no measurable impact on energy resolution (only fake-cluster counting, handled by V3 on 2.24/3.49 H2); writing convention "Purpose, inputs, outputs → Procedure → Results" approved and stored (curator), the two-semicolon convention not approved; f5 mild-clipping threshold nclip ≤ 3 confirmed; the 81k Ta members whose energy the credible f1/f3 alternative removes stay deferred (no correction is applied); user asked to commit everything except large files (git-ops)
- 2026-11-18 decision (user): "Rates at different energies not just per beam energy but also in impact on different energy slices: 10-50 MeV, 50-100MeV, 100-300, 300-500, 500-1000, 1000-1500, 1500-2500, 2500-3500, 3500-4000 etc"; P12 dispatched — per-channel category rates and cluster-level flagged fractions / ⟨fE⟩ per energy slice, exotic_events/exotic_P12_rates_vs_energy_20261118.{pdf,tsv}, figure for the chapter

## Exotic events: impact packages (approved 2026-11-18)

Inputs: exclusive categories f4 → f5 → f1 → f2 → f3 from `exotic_events/exotic_p4_module_rates.py` (`flags_of`), run 24937 excluded, PbWO4 only, flags are markers (BASE_SLOPE/CLIP_SINGLE off).

- [x] P6c consolidation: `exotic_events/exotic_summary_20261118.pdf` — latest number tables (definitions and priority, exclusive rates per run type, saturation, per-module tops, amplitudes, f2 diagnosis, decisions) followed by every P1–P6 plot in order with section pages; P1–P3 marked as pre-exclusive definitions.
- [x] P7 flags at cluster level: join categories to reconstructed clusters; fraction of clusters with ≥1 flagged module by role (seed/neighbour), category, run type, radius ring, energy bin; energy fraction in flagged modules; same for ep elastic and Møller samples; fake clusters (flagged seed).
- [x] P8 energy error per category: alternative energy per flagged channel (f1 sloped baseline, f3 stretched single, f4 clipped main kept, f2 single wide, f5 saturated single), ΔE per category, propagated to ΔE_cluster/E by run type and ring.
- [x] P9 impact on calibration and resolution: ep and Møller benchmark default / flagged clusters vetoed / flagged channels corrected; per-module Δgain (Ta inner ring), elastic-peak shift, σ/E, tails.
- [x] P10 impact on yields: ep and Møller count loss under veto per run type and ring; Ta inner-ring f1 rate vs beam current/occupancy per run.
- [x] Gate after P9/P10: veto + efficiency vs correct-in-place vs ignore (H2) / Ta separately.
- [x] P11 write-up chapter (after P10): one concise, logical chapter on the exotic-event study in the write-up TeX; new figures in a new folder under `writeup/figs/`; TeX updated; `writer` lead with the project's writing rules.
- [ ] P12 rates vs energy slice (10–50, 50–100, 100–300, 300–500, 500–1000, 1000–1500, 1500–2500, 2500–3500, 3500–4000, >4000 MeV): per-channel category rates and cluster-level flagged fraction / ⟨fE⟩ per slice; figure + paragraph in Section 6.2.

## Todo

- [ ] T0 refit S2 at x1/x2/x3; per-waveform panels; t0 tests
- [ ] T1 residual profile per class/amplitude/phase; 20 examples with the band; D3 decision line
- [ ] T2 flat-top width vs A; amp/integral check; A_sat table; 30 examples
- [ ] gate 1: user confirms T3a, T3b (only if D3 fires), T3c
- [ ] T3 PR switches, tables, `test_T3.py`, examples PDF
- [ ] T4 old vs new on S2 + 50 random; D6 criteria
- [ ] T5 overlay grid, map, bias, fake rate; grid examples
- [ ] T6 third sample variant, R_N table, dropped dt-r, examples
- [ ] T7 W121/W493 with the new model; in-time shape check
- [ ] gate 2: user decides on T8
- [ ] T8 BM port, off-identity, variants, delta vs V3, examples
- [ ] gate 3: user decides on adoption
- [ ] T9 full scans, report, TEX, note

## Amendments
