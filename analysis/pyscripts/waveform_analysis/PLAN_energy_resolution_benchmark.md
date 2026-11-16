# Plan: energy_resolution_benchmark.py (standalone legacy vs timing-gate benchmark)

Paths are relative to `analysis/pyscripts`. B = `benchmark_hycal_timing_test_20261021.py`, PR = `waveform_analysis/pileup_analysis/pileup_rates.py`, bank = `waveform_analysis/pulse_templates/<RT>/templates_<RT>.npz`.

## Open Questions (answer, then dispatch `executor`; the recommended answer is the plan's working assumption)
1. Møller rejection (`--moller-reject kin`, part of the reference sample): drop it; the two-peak CB models the Møller peak, the cut widens the CB core by +0.168 pp and 51-65 % of its Gaussian gain is a cut artefact (notes/PbWO4_energy_resolution_20261030.md:50-52); the 0.73 GeV estimator question stays recorded as undecided (notes/PbWO4_energy_resolution_summary.md:52). **Recommended: drop** (revisit in P2). It changes the sample definition versus R5nr-ri, so it is the user's call.

## Goal
One standalone script, `waveform_analysis/energy_resolution_benchmark.py` (scan + report), a shell driver, outputs under `waveform_analysis/energy_resolution/`, reproducing B's legacy-vs-gated PbWO4 energy-resolution comparison with the new MF and bank. Non-goals: per-crystal gains, position maps, Møller rejection, LMS, time offsets, pulse dumps (P2 or dropped). Expected size: script 700-900 lines (about 300 copied MF/Bank, 250 copied selection and fits), driver about 30 lines, no change to existing files.

## Approach
Per run: decode -> MF per PbWO4 channel (WaveAnalyzer for glass and clipped channels) -> hits -> three island clusterers (legacy, gate 8, gate 16) -> per-cluster records in one npz. Per run type: the report step pools the npz files, selects, fits per split and writes tables and a PDF. Caching cluster records separates the expensive MF pass from every later fit variant (the open item of notes/PbWO4_energy_resolution_20261105.md:116) and gives P2 its input. Rejected: a one-shot script fitting in the scan step (every fit change reruns the MF); vendoring mf_pulse (1700 lines, forbidden by the standalone rule); re-implementing the gains and maps now (the plan would exceed the change).

## Decisions
- D0 Settled questions. (a) Per-crystal gains and position correction stay out of P1 and become a report-step package on the cached records (P2): the reference trees (`output/r5ri_20261104`, `output/r5_20261028`) are deleted and both solvers are large; P1 reproduces the frozen benchmark at matching flags, not the R5nr-ri numbers. (b) Gate windows `GATE_NS = (8, 16)`: the brief asks for 8 ns, the reference used `--window 16` (Commands.md:156) and tight gates lose peak yield (notes/PbWO4_energy_resolution_summary.md:136-139); both run in one pass at negligible cost; if the user later picks one, drop the other clusterer and its npz keys. (c) Runs = the bank's run lists (`waveform_analysis/run_template_select.sh:12-16`), Ta per run; 2.24 GeV C (25308) has no bank and is not run; validation on run 25320 (in the bank list, the 11-08 note's to-do run). (d) Resolution: headline = two-peak CB core σ/μ, secondary = iterated ±2σ Gaussian, mean ± SE over splits (notes/poscorr_form_choice_20261106.md:57-71), plus the IQR/1.349 width as a model-free check.
- D1 Layout: `energy_resolution_benchmark.py scan <out.npz> <evio>...` (one run) and `... report <dir>`; positional only, every setting a named constant (pattern PR:60-115, `template_extract.py`). Subcommands match the `select`/`build` pattern of the bank scripts.
- D2 Per-run config: `run_config(run)` copied from PR:457-469 (beam, calibration file, dead list from `database/runinfo/general.json`); `b.set_hycal_calib(DB_DIR/cfg["calib"])`, `b.set_run_number(run)`; abort if `core.run_cfg.Ebeam` differs from `cfg["beam"]` by > 0.1 %. Dead modules get no channel (PR:605). Run type from `RUN_TYPE_OF` (PR:70-78, same lists as the driver); target mass M = M_p (H2), 180.948 u (Ta).
- D3 Event gate: `EventType.Physics` and `trigger_bits == 0x100` (B's `_common.PHYSICS_TRIGGER_BITS`), not PR's mask; LMS events are excluded by the exact test.
- D4 Pulses: `Bank(path, names, pwo4)` and `MF` copied verbatim from PR:99-455 with its constants (KMAX 6, NSIGMA/MIN_AMP 5/5, MIN_SEP_NS 8, MERGE_NS 6.5, R_MODEL 0.02; a change from B's 5/5/4/0.3, not a drop); blocks of 256 events; `mf.run(Y, bank.row[mi], overflow)`, without PR's PRECUT_ADC/A_MIN cuts. Hit time t = t0 + bank `t50[module]` (B's `p_time`, mf_pulse.py:1298); hit energy E = `mod.energize(A * iint_wa[module])` (PR:680 `module_gain`; the reference used `--mf-scale 1.0`, f(A) was neutral, fAri). Glass channels and clipped PbWO4 channels without a valid pulse: `dec.WaveAnalyzer(core.daq_cfg.wave_cfg)` peaks, E = `energize(integral)` (B:1594-1609 fallback).
- D5 Hit rules (restated for the new MF, which keeps small secondaries): a pulse is in-window if `PRE_LO < t < PRE_HI` = (100, 200) ns (strict, no offsets). Legacy: one hit per module, the in-window pulse with the largest A (B:1602-1627). Gated: every in-window pulse is a hit; no amplitude cut beyond the MF's own `max(5 sigma_A, 5 ADC)` and the clusterer's `min_module_energy`; the 8/16 ns coincidence is `seed_time_window` of the clusterer (B:2708). Per module record `npulse` = number of in-window valid pulses.
- D6 Clusterers: copy the builder config fields as B:2700-2725, `min_center_energy = 30` (notes/PbWO4_energy_resolution_summary.md:54; database default is 10), `log_weight_thres` stays 5.2, `seed_time_window` = -1, 8, 16. Keep the builder's Geant4 profile by calling `set_config(base_cfg)` before `set_config(cfg)` as B does.
- D7 Scan npz (per path p in legacy/g8/g16): `cl_<p>_{ev, split, seed_idx, row, col, seed_x, seed_y, x_cl, y_cl, e_cl, n_mem, mem_ptr}` and CSR members `mem_<p>_{idx, e, t, npulse}`; per event `ev_number, n_big_<p>` (clusters with E > SIGNAL_MIN); header keys run, run_type, splits, beam, calib, dead, bank sha1, n_seen/n_pass, wall. Storing the split index per cluster keeps the per-split statistic of the reference.
- D8 Selection (report): events with exactly one cluster above `SIGNAL_MIN = 500 kE`, kE = E_beam/2239.51; pool = seed r < `INNER_R = 200` mm, seed E >= `SEED_MIN = 200 kE`, seed outside `EXCLUDE_ROWCOL = (15, 20, 15, 20)` (1-based inclusive against the 0-based C++ row/column, copy B:615-634 `_select_pool` exactly; notes/PbWO4_energy_resolution_20261023.md:295-299 records the pitfall). E' per seed module from B:257-287 (elastic kinematics minus energy loss, theta from the module centre).
- D9 Module set: copy B:823-955 `fit_elastic_window_gains` (pooled Gaussian then CB -> per-module CB prefit on the Møller-aware range PREFIT_RANGE B:248, median fallback below MIN_FIT_EVENTS 100, gscale gap test >= 3 sigma_M, n(mu_M ± 2 sigma_M) >= `MIN_PER_MODULE = 30`) and discard the gain it returns (refine off, post = pre); the used-module set is then B's by construction.
- D10 Fits per split and path on the used-module pool spectrum, range E' ± 800 kE, bin 20 kS MeV (kS = sqrt(kE)): copy B:957-1013 (`fit_peak`, `_gauss_fixed_hist`, the fixed ±150 kS seed), B:296-510 (`fit_peak_cb`, `fit_cal_cb`: two-peak CB elastic + Møller, fallback single CB, seeded by the fixed Gaussian) and B:1016 (iterated ±2σ Gaussian); add IQR/1.349 and the fraction outside ±3 sigma_CB. Quote mean ± SE over splits, paired Δ(gated − legacy) per split per gate, and the peak yield as B:3623-3645 (pool events within ±2.5 s_ref of the module target, s_ref = the legacy fixed-window σ/μ, modules used in both paths).
- D11 Diagnostics in the report: pile-up systematic = the same fits after dropping clusters with any member `npulse >= 2` (expected 0.3-3 % of clusters, `waveform_analysis/pileup_analysis/pileup_rates_out/summary.txt`); mu_CB/E' − 1 per run (calibration check). Outputs: `energy_resolution/<RT>/run_<run>.npz`, `energy_resolution/report/{resolution.tsv, resolution.md, report.pdf}` (one page per run type: spectra of the three paths with fits, table).
- D12 Driver `waveform_analysis/run_energy_resolution.sh` (pattern `run_template_select.sh`): the same run lists, `xargs -P NPROC`, one scan per run, then `report`.
- D13 Dropped from B (one reason each): rawsum pulse method (MF is the chain); LMS correction and uproot (reference `--lms-corr off`); time offsets (`--time-calib off`); `--mf-scale`/kappa and f(A) (1.0, neutral); gain refinement and pre/post duality (off, post = pre); `--cal-method cb-peak/median/percrystal`, `--gain-file`, `--pos-corr`, `--nl-members` (P2); `--cluster-dump`, `--pulse-dump`, `--write-cluster-tsv` (npz replaces them); `--mf-shape*`, `--mf-nb-*`, `--single-pulse`, `--cb-sigma-tie`, `--target-peak`, `--fit-select` (CB is the headline, the label only); Møller rejection and pair vetoes (Q1); `--veto-*`, `--cluster-*` overrides, GEM/DAQ/map overrides, `--max-events`, CSV output (diagnostic or obsolete).

## Outline
- `waveform_analysis/energy_resolution_benchmark.py`: constants block (D3-D10); `run_config`, `Bank`, `MF` (PR:457, 119, 171); `open_run(evios) -> core, hycal, names, pwo4, bank, dead`; `scan(out, evios)`: event loop as PR:580-620 with D3 and without PR's PRECUT_ADC/A_MIN, block fit, `hits_of(block) -> per module (t, E, A) lists`, `cluster(cl, hits) -> records`, `np.savez_compressed` (D7); `eprime(row, col)` kinematics (B:257-287); `select_pool`, `module_set` (D8-D9); `fit_peak`, `fit_cb2`, `fit_gauss2s`, `robust_width` (D10); `report(dir)`: per run type and split -> rows -> TSV/MD/PDF (matplotlib PdfPages as PR's report). Invariants: MF row index is `bank.row[mi]`, never the module index (PR:546); cluster energies come from C++ with nl (B:2109); hit times are on the t0 + t50 axis the pre-window was defined on.
- `waveform_analysis/run_energy_resolution.sh`: as `run_template_select.sh` with `energy_resolution_benchmark.py scan`/`report`.

## Slices
1. Scan step end to end on one split. Units: the script's scan half. Tests: `scan_smoke: 2 000 events of prad_025320.evio.00050 -> npz with the D7 keys, n_pass > 0, cl_legacy count > 0 and cl_g8 count > 0`; `hits_legacy_one_per_module: a block -> at most one legacy hit per module index per event`; `dead_excluded: run 25800 -> no member with a dead module index`. Checkpoint: `python3 waveform_analysis/energy_resolution_benchmark.py scan /tmp/er_test.npz <split>` (run only when the machine is free) and `python3 -c "import numpy as np; z=np.load('/tmp/er_test.npz'); print(sorted(z.files))"`.
2. Report step on the slice-1 npz. Units: the report half. Tests: `pool_matches_B: the same split through B's _select_pool -> identical seed_id set` (executor writes a 20-line check script under /tmp, not the repo); `cb_fit_recovers: synthetic Gaussian + Møller peak -> sigma within 3 %`; `paired_delta: legacy == gated records -> Δ = 0 exactly`. Checkpoint: `python3 waveform_analysis/energy_resolution_benchmark.py report waveform_analysis/energy_resolution` writes the three report files.
3. Driver and validation. Units: the shell driver. Tests: `driver_dry: NPROC=1 with one run -> one npz, report runs`. Checkpoint: the validation gate below on run 25320.

## Verification
Reference regeneration (machine free): run B (sha1 3293ab96…, Commands.md:9) on the splits of run 25320 with the final-line flags (Commands.md:156) minus `--gain-file`, `--pos-corr`, `--moller-reject*`, `--nl-members`, with a cluster config rebuilt under /tmp from `database/reconstruction_config.json` with `min_center_energy` 30 (notes/ep_resolution_walkthrough_20261102.md:220), the runinfo 2.24 GeV calibration file, the best-50 bank `output/templates_spline_best50_20261026/2p24GeV_H2/templates_2p24GeV_H2.npz`, `--exclude-rowcol 15 20 15 20`, `--window 16`, the default `--cal-method elastic-window`. Gate: per split, |σ/E_CB(new) − σ/E_CB(B)| <= 0.05 pp on legacy and on gate 16 (the bank swap alone moved 3.49 GeV by 0.007 pp, w1ri; per-split SD is 0.03-0.04 pp), |mu_CB ratio − 1| <= 1 %, used-module set differs by <= 3 modules, peak yield ratio within 1 %.

## Risks
- Clipped seeds at 3.49 GeV fall back to WaveAnalyzer integrals (different energy estimator); count them (`n_clip` per path) and report. Lead-glass hits use WaveAnalyzer as before.
- The new MF's extra small pulses can shift the legacy hit choice only when a secondary exceeds the main pulse (rare); the gated paths see more small hits, which `min_module_energy` (1 MeV) does not remove; the pile-up flag quantifies it. The MF change itself is unquantified on real data (notes/PbWO4_energy_resolution_20261108.md: pseudodata only).
- `t50` of the new bank differs from the best-50 bank by about 0.04 ns median (planner-worker npz comparison, 2.24 GeV H2 only; other types unchecked).
- Memory: per-run cluster records of a 13-split Ta run; CSR float32 keeps it under about 200 MB.

## Scope Notes
In: the script, the driver, the report, the validation on run 25320. Out (user-fixed): mf_pulse or any top-level module import; option flags; history comments; per-crystal gains, poscorr, Møller rejection (P2); changes to B or PR; new scans while the pile-up scan runs; the full four-type driver run (the user starts it).

## Assumptions
- `core.run_cfg.Ebeam` equals runinfo `beam_energy` for every run (checked by the D2 abort).
- The builder's Geant4 profile survives the second `set_config` (checked by comparing `cl.get_config()` fields after setup, as B:2718 comments).
- `HyCalCluster.reconstruct_matched()` records expose `hit.x/y/energy/center_id` and `cluster.hits[].index/energy/time` (B:2897-2904; checked at slice 1).

## Work Packages
### P1: script, driver, report, validation
- Goal: the standalone benchmark reproduces B at matching flags on run 25320; done when the Verification gate passes
- Owns: `waveform_analysis/energy_resolution_benchmark.py`, `waveform_analysis/run_energy_resolution.sh`, `waveform_analysis/energy_resolution/`
- Deps: none
- Stages: execute
- Status: execute
- Artifacts: waveform_analysis/PLAN_energy_resolution_benchmark.md
### P2: per-crystal gains and position correction on cached records
- Goal: add member gains (elastic-window or least-squares) and the chosen poscorr forms (notes/poscorr_form_choice_20261106.md §1) as report-step passes over the P1 npz; done when the R5nr-ri numbers are reproduced within 0.05 pp
- Owns: the report half of the script and `energy_resolution/report/`
- Deps: P1
- Stages: plan, execute
- Status: blocked: plan: P1 npz schema and the Q1 answer needed
- Artifacts: plan-er-gains-poscorr.md at the default path

## Improvements (ranked; 1-5 are in P1 by D7, D11, D0b, D10, D11)
1. Cached cluster records: fit variants without rerunning the MF (note 20261105 open item).
2. Pile-up flag systematic from the MF (D11; in-window extra-pulse rates 0.3-3 %).
3. Both gate windows in one pass with the peak-yield check (D0b).
4. Robust width and tail fraction next to CB/Gaussian (D10).
5. mu_CB/E' − 1 per run as a calibration check.
6. Deferred: module-level resolution maps (D7 records allow a later plotter); gate window from the measured sigma_t (0.2-1 ns per module, summary:408; 8 ns is already about 8 sigma); sigma(t)-weighted amplitudes (changes the MF, breaks comparability; dropped).

## Log

## Todo
- [x] Slice 1: scan step (constants, run_config, Bank/MF copy, decode, hits, clusterers, npz)
- [x] Slice 2: report step (selection, module set, fits, tables, PDF)
- [ ] Slice 3: driver; reference regeneration with B on run 25320; validation gate
