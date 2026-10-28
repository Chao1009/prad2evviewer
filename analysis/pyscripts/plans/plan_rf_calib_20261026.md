# Plan: per-run RF time calibration of the PbWO4 modules (`rf_time_calib.py`)

**Status:** executed 2026-10-26 (S0–S7); results in the Progress log, `output/rf_calib_20261026/` and note 20261025 §9.7.
**Branch:** `waveform_jz`. **Inputs:** `notes/PbWO4_energy_resolution_20261025.md` §9.1, §9.5, §10, §11(h); `plot_module_timing_rf.py`, `plot_module_timing.py`, `rf_time_dump.py`, `module_timing_dump.py`; `database/hycal_rf_offsets/24840.json`, `prad2det/include/HyCalRfOffsets.h`, `prad2det/include/RfTime.h`, `prad2det/include/EventData.h:289-304`, `analysis/src/Replay.cpp:1093`; `database/hycal_time_offsets/time_calib.json`; `database/hycal_map.json`; `output/module_timing_20261025/{rf_summary.tsv, module_timing_rf.tsv, module_timing_evtref_tcal.tsv, run_rf.log}`; dump and RF npz headers; a planner diagnostic on runs 25800, 26138, 26230 (read and run 2026-10-26, §2).

Paths are relative to `/Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts/` unless they start with `database/`, `docs/`, `prad2det/` or `analysis/` (repo root). `RF:` = `plot_module_timing_rf.py`, `PMT:` = `plot_module_timing.py`. P = RF bunch period (4.00771 ns measured per split, `RF:` pass 1), fold(x) = x − P·round(x/P) ∈ (−P/2, P/2].

---

## 1. Goal and approach

**Goal.** One script, `rf_time_calib.py`, that derives from the existing dumps (`output/module_timing_20261025/dump/`, 111 splits, 13 runs) joined to the RF dumps (`.../rf/`) a per-run, per-module RF timing offset `o_{m,r}` such that fold(t_module − t_RF − o) is centred at zero module by module; writes one JSON per run in the `database/hycal_rf_offsets/24840.json` schema to `output/rf_calib_20261026/<run>.json`, a per-(run, module) TSV and a per-run TSV, one before/after PDF per run; resolves the ±P/2 branch of the continuous residual to the `time_calib.json` prior with the leave-one-out event median; and validates the result, including a cross-check against `RF:`.

**Approach.** Reuse the join, fold and nearest-hit code of `RF:prep` and the peak estimator `RF:fit_unwrapped` (80-bin circular histogram unwrapped around the smoothed argmax, `PMT:fit_core` ±2σ Gaussian core) exactly as they stand, per (run, module) instead of per (type, module). Two energy tiers (E ≥ 100 MeV plateau, all energies). The JSON offset is a phase and has no branch (the consumer re-folds, `HyCalRfOffsets::apply`); what needs a branch is the *continuous* deviation D from the prior, obtained as fold(o − τ − G_r) + kP with a run-global reference phase G_r and k chosen by the event median e (`PMT:event_refs`, in-window mates, corrected once with the first-pass D). The planner diagnostic (§2) shows that fold(o − τ − e) is the same number for every module of a run to 0.04–0.09 ns MAD, shifted modules included, so G_r needs no "good module" selection and no per-board structure. Run 24386 is detected as bimodal and written as unusable; within-run drift (24503) is flagged, not corrected.

**Rejected alternatives (one line each).**
- Unbinned iterative truncated median as the estimator: equivalent in precision at N ≥ 300, but every σ in the note comes from `fit_core`, and the cross-check against `RF:` is exact only with the same estimator.
- Per-split-range offsets for 24386: the two phases coexist 50/50 within splits (§9.5, finding 3), so no split range is single-phase.
- Per-board reference phase G_b for the branch: on board (3,15) six of eight modules are shifted at 26230, so the board median absorbs the shift (critic C1); the global G_r has MAD 0.04–0.09 ns and G_b − G_r is ≤ 0.1 ns on every uncontaminated board (§2).
- Board-level offsets (one phase per board per run): the 25800 → 26230 jump is global (+2.00 ns, within ±0.07 ns on every board except the two that hold the central block, §2), and the 26-module shift is per module and per run.
- Branch from the prior alone (nearest branch): wrong for W494, W528, W629 at 26230 (RF residual +1.7/+1.8/+1.7 ≡ −2.3/−2.2/−2.3; the event median is −1.97/−1.77/−1.76).
- t_seed as the continuous reference: one module, so a stale seed offset biases every mate; the leave-one-out median is robust and already vectorised.
- Per-board fallback tier for modules below N_MIN: the run-to-run phase is global (§2), so the prediction fold(τ_m + Φ_r) from the prior and the run phase is the simpler fallback; it fires for ≈ 76 outer-ring modules in the 4-split run 26138 and a handful elsewhere.
- A new evio pass: everything is in the dumps.

## 2. Facts the plan relies on (verified 2026-10-26)

| Fact | Where |
|---|---|
| Dump npz per split: `ev_num`, `trig`, `m_ev`, `m_idx`, `m_t` (p_time), `m_E` (MeV), `m_npk`, `m_inwin`, `m_fb`, `module_names` (1735, p.hycal order), `module_is_pwo4`, `templates`, `calib_file`, `run`, `split`. No fallback rows (`m_fb` all False in the four splits checked). `m_inwin` fraction 0.99 / 0.98 / 0.99 / 0.71 (24590 / 25320 / 24503 / 26230) | `module_timing_dump.py` header; npz read |
| RF npz per split: `ev_num`, `trig`, `ts`, `rf_off`, `rf_t`, `rf_ch` (0 = A, 8 = B). Join by `ev_num` matches 100 % | `rf_time_dump.py`; `rf_summary.tsv` |
| Splits per run: 24386 ×11, 24503 ×10, 24512 ×7, 24590 ×10, 24603 ×6, 25228 ×11, 25257 ×7, 25260 ×8, 25299 ×8, 25320 ×8, 25800 ×8, 26138 ×4, 26230 ×13 | `ls dump/` |
| `RF:prep(dp, rp, cfg)` returns `run, split, n_dump, n_matched, n_sel, a_hits, sp_sum, sp_n, P, trig_equal, ri, pm, t, t_rf, phi, d, m_idx`; cfg keys `names_hash, mask, min_e, pmap, npw, off_p`. It does not return the event mask `sel`. `RF:pass1` deletes keys `d, ri, pm, t, t_rf, phi, m_idx` before returning; `prep` is called only by `pass1` and `pass2` | `RF:` lines 110-186 |
| `RF:fit_unwrapped((h, P))` → dict with `ok, mu, sigma, mu_err, sigma_err, lo, hi, amp, a, cen, mu_u`; `RF:unwrap`, `RF:fold`, `RF:circ_median`, `RF:circ_mad`, `RF:_bin`, `RF:NB = 80`, `RF:N_SUB = 32`, `RF:MIN_FIT = 200` | `RF:` |
| `PMT:load_geometry(run)` → `name, id, x, y, sx, sy, pwo4` (needs the C++ pipeline, `source setup_env.sh`); `PMT:load_offsets(path, geo)` → `off, nres, missing` (resolved entries only; W1123, W1124 unresolved); `PMT:fit_core(x, c, bg, smooth_ns, max_it, tol)`; `PMT:event_refs(m_ev, t, E, sel, n_ev, ref_min_e, ref_min_n)` → `t_evt, t_loo` (eligibility `E > ref_min_e`); `PMT:centres(lo, w, nb)`; `PMT:EVT_LO = −60, EVT_NB = 600, W_BIN = 0.2`; `PMT:TYPE_OF_RUN, TAGS, LABEL, RUN_NAME, TRIG_EXCL, NMAX_JOBS, TIME_CALIB, _mpl, text_page` | `PMT:` lines 86-165, 372-437, 448, 601-700 |
| `database/hycal_map.json`: list of `{n, t, geo{x, y, sx, sy, sec, row, col}, daq{crate, slot, channel}}`; 1152 PbWO4 modules on 88 (crate, slot) boards, 1–16 modules per board | json read |
| `database/hycal_rf_offsets/24840.json`: `{"default": 0.0, "modules": [{"name", "offset_ns"}]}` (84 modules, offsets in (−2, 2)). Loader `HyCalRfOffsets.h`: reads `default` if numeric and `modules[].name/offset_ns`; ignores every other key; unknown names skipped; `apply(idx, dt) = FoldRfDelta(dt − off)` with `RF_PERIOD_NS = 4.008`. Consumer `Replay.cpp:1093,1677` (`cl_dt_rf` of the cluster centre module time) | headers |
| `time_calib.json`: 1152 modules, `offset_ns` (−8.9 … +7.0), `module_id`, `resolved`; reference W495; derived on the 2.24 GeV geometry run 25320 | json read; `RF:` docstring |
| Per run, `RF:prep` rows (selected events, PbWO4, finite time), modules with N ≥ 300: 26138 (4 splits, the smallest run) 235 at E ≥ 100 MeV / 1076 at all energies (2 modules with no rows at all); 25800 327 / 1142 | count with `RF.prep` over the 26138 and 25800 dumps (critic, round 2) |
| 26 X17 central modules: W492–W495, W525–W529, W559, W560, W563, W564, W593, W594, W597, W598, W627–W631, W663–W665, W697, on boards (3,18) ×6, (3,8) ×6, (3,15) ×6, (3,6) ×4, (3,14) ×2, (3,5) ×2. Seven halo modules around them: W491 (3,18), W496 and W530 (3,8), W632 (3,6), W666 (3,5), W661 and W662 (3,15). The remaining modules of those boards have type-level `offset_tcal` −0.28 … +0.22 at X17 | `module_timing_rf.tsv`; `hycal_map.json` |
| **Planner diagnostic (10-26)**, per module circular median of φ at E ≥ 100 MeV, N ≥ 300, and e = median of t − τ − t_loo with in-window rows and mates (E > 10 MeV, ≥ 3 mates), per run: (a) 25800 → 26230 RF jump: global +2.00 ns; relative to it every board with ≥ 4 modules is within ±0.07 ns (within-board MAD 0.003–0.08) except (3,15) +1.55 and (3,8) +0.32, the two boards holding most of the central block; (b) fold(ρ − e) over all measured modules: MAD 0.042 / 0.039 / 0.090 ns (25800 / 26138 / 26230); per-board medians minus the global one: max 0.095 / 0.079 / 0.316 (the 0.316 is (3,15)); (c) with the global G_r and pass-A e, \|fold(ρ − G_r) − e\| < 0.5 ns for 100 / 100 / 99.1 % of modules, the 0.9 % being W494, W528 and W629 whose fold sits at +1.78/+1.88/+1.80 against e −1.97/−1.77/−1.76 (branch k = −1); (d) the central block is shifted only in 26230: there ρ − Φ_r = −0.7 … −2.3 (W494, W528, W629 ≈ −2.2) with e = −0.67 … −1.97 of the same sign; at 25800 and 26138 the 26 modules have \|e\| ≤ 0.2 except W525, W528, W597, W697 at e = −0.78 … −1.30 (−1.00 … −1.77 at 26230, where they follow the block); the seven halo modules are at e = −0.25 … −0.69 in 26230 only; 26138 vs 26230 global RF jump −0.02 ns | scratch `board_diag.py`, outputs `board_diag_*.txt` in the planner scratchpad; the script's items (iii)–(v) reproduce them |
| Event-median reference with the prior, type level: `mean_x17` −0.22 … −1.17 ns for the 26 (σ 0.7–1.4), `mean_2p24` −0.27 … +0.05 | `module_timing_evtref_tcal.tsv` |
| `rf_summary.tsv`: `run_offset_tcal` (fit of the pooled tcal phase histogram `ht` per run), `split_offset_span_tcal` (0.003–0.004 for 9 runs; 24512 0.058, 26230 0.042, 24503 0.235, 24386 0.878), `n_mod_jump` (25800 vs reference 26230: 1139/1150; 26138: 31/1132) | file |
| `RF:` full run: 111 files, pass 1 17 s, pass 2 18 s, 14 425 per-(module, run) fits 8 s (31 argmax fallbacks), PDFs 262 s, 12 workers | `run_rf.log` |

## 3. Definitions used by the script

| Symbol | Definition |
|---|---|
| φ | fold(t_m − t_RF) per row, P of the row's split, t_RF = nearest channel-A hit (`RF:prep`, unchanged). Rows: dominant PbWO4 pulses of events with (trig & 0xf00) ≠ 0, (trig & TRIG_EXCL) = 0, matched in the RF file, ≥ 1 A hit; `m_fb` rows excluded (none exist) |
| tier `e100` | rows with m_E ≥ E_MIN = 100 MeV (plateau of the walk, §9.4/§9.5) |
| tier `allE` | all rows |
| halves | split index parity within the run's sorted split list (even / odd); every histogram is kept per half and summed for the full estimate |
| H[r, m, tier, half] | 80-bin circular histogram of φ (bin P/80 = 0.050 ns, `RF:_bin`) |
| c, σ, c_err | `RF:fit_unwrapped((H, P_r))` → c = fold(`mu_u`) (fit `mu`, argmax centre on failure, flag `argmax`), `sigma`, `mu_err`; P_r = entries-weighted run period as in `RF:` |
| o_{m,r} | the RF offset of the JSON: c of tier `e100` if N_e100 ≥ N_MIN = 300; else c of tier `allE` if N_allE ≥ N_MIN (flag `allE`); else, if the module has a prior, the prediction fold(τ_m + Φ_r) (tier and flag `pred`; no σ, no dhalf, excluded from every statistic over "measured" modules); else unresolved (omitted from `modules`, flag `unresolved` in the TSV) |
| τ_m | `time_calib.json` offset (`PMT:load_offsets`; 0 and flag `noprior` for W1123, W1124) |
| ρ_m | fold(o_m − τ_m): residual of the prior in the RF frame (measured tiers only) |
| Φ_r | RF run phase: `circ_median` of ρ over the measured modules of run r. Not the entry-weighted pooled fit `rf_summary.tsv:run_offset_tcal`, which the high-N central block pulls by ≈ 0.8 ns at 26230 |
| ρ^r_m | fold(ρ_m − Φ_r): the module's deviation from the run phase ("before" quantity of item (iii)) |
| Δ_{b,r}, mad_b | board phase: `circ_median` / `circ_mad` of ρ over the measured modules of board b = (crate, slot) in run r, ≥ 4 modules; diagnostic of item (iv) only |
| e_m | continuous reference: `fit_core` mean (bg=False, smooth_ns 1.0) of the 600-bin histogram over [−60, 60) ns (`PMT:EVT_LO/EVT_NB/W_BIN`) of t_m − τ_m − t_loo, where t_loo = `PMT:event_refs` leave-one-out median over the event's other PbWO4 pulses with m_E > 10 MeV and `m_inwin`, ≥ 3 mates, and the module's own row `m_inwin`; same tier rows as o_m; "measured" if N_e ≥ N_MIN, fit ok and σ_e < 2.0 ns. Pass A: mates' times t_j − τ_j. Pass B: mates' times t_j − τ_j − D¹_j for mates with a firm pass-A branch (else − τ_j only) |
| G_r | branch reference phase of the run: `circ_median` over all measured modules with measured e of fold(o_m − τ_m − e_m); its `circ_mad` is the RF-vs-event-median consistency number of item (vi) |
| D_m (`dt_tcal_ns`) | continuous correction to τ_m for run r: fold(o_m − τ_m − G_r) + k·P with k ∈ {−1, 0, +1} = argmin_k \|fold(o_m − τ_m − G_r) + kP − e_m\|; margin = second-best − best distance; `firm` if margin ≥ 1.0 ns; k = 0 and flag `nobranch` when e is not measured. D¹ from pass A, D² (final) from pass B; flag `branch_changed` if \|D² − D¹\| > P/2. D is defined up to the run constant G_r (an overall time shift of the run is unobservable in the event frame); pass B moves G_r by a run constant (critic check: +0.09 ns at 26138, +0.11 ns at 26230), reported as `G1_ns`, `G2_ns` in `rf_calib_runs.tsv` |
| bimodal (module) | on the `allE` histogram: residual = 5-bin smoothed circular histogram − fitted Gaussian core (amp, c, σ); flagged when max residual over \|fold(x − c)\| > 0.75 ns ≥ 0.35 × amp. Calibrated on the data: 24386 (second phase 1.33 ns away, ≈ equal height) must fire, no other run above 10 % of its measured modules; if one does, raise the ratio to the smallest value that separates and record it in `rf_calib_runs.tsv` |
| bimodal (run) | > 30 % of the measured modules flagged bimodal → `usable: false`, `quality: "bimodal"` |
| drift | pass B per split: pooled `allE` residual fold(φ − o_m) over all measured modules, fitted; span of the per-split centres (`split_span_ns`) > 0.10 ns → `quality: "drift"` (usable) |
| held-out | `dhalf_ns` = fold(c_even − c_odd) of the tier used, when both halves have ≥ 150 entries; its error = √(err_even² + err_odd²) |

## 4. Decisions (the five the brief asked for)

1. **Output format.** One JSON per run, `output/rf_calib_20261026/<run>.json` (no zero padding, like `24840.json`), with `default` and `modules[{name, offset_ns, …}]` as the first keys so `HyCalRfOffsets.h` reads it unchanged; extra per-module keys (`module_id`, `offset_error_ns`, `sigma_ns`, `n`, `tier`, `flags`) and top-level `run`, `run_type`, `period_ns`, `usable`, `quality`, `split_span_ns`, `time_definition`, `created` are ignored by the loader. No NaN in any JSON (unresolved modules are omitted from `modules`). Plus `rf_calib_modules.tsv` (one row per run × PbWO4 module), `rf_calib_runs.tsv`, `rf_calib_validation.txt`, `rf_calib_<run>.pdf` ×13, `run.log`. Never under `database/`.
2. **Estimator.** `RF:fit_unwrapped` on the 0.050 ns circular histogram (unwrap around the smoothed argmax, `fit_core` ±2σ Gaussian core, no background) on E ≥ 100 MeV pulses when N ≥ 300 (statistical error ≤ 0.03 ns at σ ≤ 0.5), otherwise on all energies with flag `allE` (walk bias late by the per-type median of c_allE − c_e100, reported), otherwise the prediction fold(τ_m + Φ_r) from the prior and the run phase with tier `pred` (≈ 76 outer-ring modules in the 4-split run 26138, a handful elsewhere; the prior holds to 0.3 ns for ≳ 95 % of modules outside the 26230 central block, §2), otherwise unresolved (W1123, W1124). The circular unwrap makes the estimate exact at ±P/2; the brief's truncated median is rejected in §1.
3. **Fold ambiguity.** The JSON offset is a phase and has no branch (the consumer re-folds). The branch is resolved for the continuous residual D to the prior: (c) a run-global reference phase G_r = median of fold(o − τ − e), which the diagnostic shows is the same for every module of a run to 0.04–0.09 ns MAD (boards are not needed: the run-to-run jump is global and the board-level variant is contaminated where a board is mostly shifted); (b) the leave-one-out event median with τ-corrected, in-window mates picks k, and one pass with the mates corrected by D¹ confirms it; (a) k = 0 only when no continuous measurement exists (flag `nobranch`). Acceptance (gate): in the three X17 runs the 26 central and 7 halo modules all have a firm branch with \|D² − e²\| < 0.5 ns (critic check: max 0.09 ns at 26138, 0.08 at 26230); no module firm in pass A has \|D² − D¹\| > P/2. Expected values (report, from the diagnostic): at 26230 D² ∈ [−2.6, −0.4] for the 33 modules; at 25800 and 26138 \|D²\| ≤ 0.4 for 29 of them and D² ∈ [−1.5, −0.6] for W525, W528, W597, W697; at the 2.24 GeV runs \|D²\| ≤ 0.4 for all 33.
4. **24386 and 24503.** 24386: the bimodal detector must flag it (two phases 1.33 ns apart, 50/50 in most splits); its JSON is written with `usable: false`, `quality: "bimodal"`, empty `modules` (a loader fall-back to 0 ns), and the per-module numbers go to the TSV with flag `bimodal`/`run_unusable`; excluded from every validation statistic. 24503: a single run offset, `quality: "drift"`, `split_span_ns` in the JSON and the per-split centres on a drift page of its PDF; reported separately in the held-out summary (expected span ≈ 0.2 ns, dhalf ≤ 0.1 ns since the halves interleave). Only 24503 is expected to get `drift`.
5. **Validation.** Written by the script to `rf_calib_validation.txt`, gates marked G (the run fails if a gate fails), the rest are reports with an expectation: (i, report) per-run σ per module (tier e100 median per type, expected 0.18–0.30 ns) and the fraction with σ < 0.3 ns; (ii, G) held-out dhalf over the 12 usable runs: median \|dhalf\| < 0.05 ns and ≥ 99 % of measured (run, module) with \|dhalf\| < 0.30 ns; (iii, report) prior consistency: fraction of modules with \|ρ^r\| < 0.3 ns per run (expected ≳ 95 % at the H2 runs, 25800 and 26138; the 33 modules fail at 26230, four of them in every X17 run); (iv, report) board pattern 25800 vs 26230: per board Δ_b(25800), Δ_b(26230), the fold of their difference relative to the global jump, mad_b (expected: global jump ≈ +2.0 ns, every board within ±0.1 ns except (3,15) and (3,8)); (v, G) the 33-module table (ρ^r, e², D², margin, dhalf per X17 run and at 25320) with the §4.3 gate; (vi, G) consistency and cross-check: `circ_mad` of fold(o − τ − e) per usable run < 0.15 ns, and the per-run pooled tcal histogram fit reproduces `rf_summary.tsv:run_offset_tcal` to < 0.002 ns in all 13 runs (same data, same code; rows skipped under `--runs` when the run is absent).

## 5. Steps

### S0. Expose the event mask from `RF:prep` (10 min)

Add `sel=sel` to the `out` dict in `RF:prep` and `"sel"` to the `del` list in `RF:pass1`. No other change. Outputs of `RF:` are unaffected (`pass2` builds a fresh dict; the key never reaches a file).

Acceptance: before the edit, `MPLBACKEND=Agg python3 plot_module_timing_rf.py --runs 26138 --jobs 4 --out-dir <scratch>/rf_before`; after it, the same into `<scratch>/rf_after`; `diff` of `module_timing_rf_test.tsv`, `rf_summary_test.tsv` and `module_timing_rf_vs_integral_test.tsv` is empty.

### S1. `rf_time_calib.py`: loading and pass A (0.5 day)

- CLI: `--runs R …` (subset, same file names), `--jobs 12`, `--out-dir output/rf_calib_20261026`, `--no-pdf`. Module constants: `DUMP_DIR`, `RF_DIR` (as `RF:` defaults), `E_MIN = 100.0`, `N_MIN = 300`, `GEOMETRY_RUN = 25320`, `TIME_CALIB = PMT.TIME_CALIB`, `HYCAL_MAP = database/hycal_map.json`, `XCHECK = output/module_timing_20261025/rf_summary.tsv`.
- Parent: `PMT:load_geometry`, names check against the first dump (as `RF:`), `pmap`/`npw`, `PMT:load_offsets` → τ, board index per PbWO4 module from `hycal_map.json` (`daq.crate`, `daq.slot`; names keyed by `n`; fail if any PbWO4 name is missing), `cfg` as `RF:` with `min_e = 0`.
- Worker `passA(job)`: `o = RF.prep(dp, rp, cfg)`; rows `ri`, `pm`, `phi`, `E = d["m_E"][ri]`, `inwin = d["m_inwin"][ri]`, `fb = d["m_fb"][ri]` (drop fb rows, count them); half = parity of the split's index in the run's list (passed in the job). Fill `H[tier, half]` (npw × 80, int32) with `RF:_bin`. Event reference: `t_evt, t_loo = PMT.event_refs(m_ev, m_t − τ[m_idx], E', sel, n_ev, 10.0, 3)` with `sel = o["sel"]` and `E'` = `m_E` with −1 for rows whose `m_inwin` is False (so they are not eligible); fill `HE[tier]` (npw × 600, int32) with t_m − τ_m − t_loo for selected rows with `m_inwin` and finite t_loo. Also the pooled tcal φ histogram `ht` exactly as `RF:pass1` (for item (vi)). Return run, split, half, P, n counts, H, HE, ht.
- Parent aggregates per run: H (sum over splits, kept per half), HE, ht, P_r (entries-weighted as `RF:`), split list.

Acceptance: `--runs 26138 --no-pdf` finishes in < 1 min; the printed per-run P equals `rf_summary.tsv:P_ns` (4.00771); the sum over halves of H[`allE`] per module equals the row counts; no fb rows.

### S2. Per-(run, module) estimates, tiers, run and board phases, bimodality (0.5 day)

- Fits in a pool: `RF:fit_unwrapped` for every (run, module, tier, half ∈ {even, odd, both}) with N ≥ 150 (half) / N ≥ N_MIN (both). Tier choice and o_m as in §3; `offset_error_ns = mu_err` (argmax fallback: error = P/80).
- Φ_r, ρ^r; `pred` fill; Δ_{b,r}, mad_b per (run, board); unresolved list.
- Bimodal flag per module (`allE` histogram) and per run (> 30 %); `quality`, `usable`.
- Walk bias per run type: median over modules with both tiers measured of fold(c_allE − c_e100); printed and stored in `rf_calib_runs.tsv`.

Acceptance: on the full set, 24386 is the only run with `usable: false`; every other run has ≥ 200 modules in tier `e100` (26138: 235, 25800: 327, §2) and ≥ 1000 measured (`e100` + `allE`; 26138: 1076); `pred` counts printed per run (26138 ≈ 76); the argmax-fallback count is printed and ≤ 1 % of fits.

### S3. Branch: e¹, D¹, pass B, e², D² (0.5 day)

- e¹ from HE of pass A (`fit_core`, bg=False, smooth_ns 1.0, x = `PMT:centres(EVT_LO, W_BIN, EVT_NB)`), G_r, D¹, margin, `firm`.
- Worker `passB(job)` = pass A's event-reference part only, with mates' times t_j − τ_j − D¹_j for firm modules, plus the per-split pooled residual histogram fold(φ − o_m) over measured modules (80 bins, allE tier) for the drift span, and the pooled held-out residual histogram (even-half rows with the odd-half offsets c_odd, odd-half rows with c_even; modules with both halves measured) for the "after" panel of the run PDF. Return HE², per-split pooled histogram, held-out pooled histogram.
- e², G_r², D², margin², flags `branch_weak` (margin < 1.0), `branch_changed` (\|D² − D¹\| > P/2), `nobranch`.
- `split_span_ns` per run from the per-split pooled fits; `drift` flag.

Acceptance: the §4.3 gate on the 33 modules in 25800, 26138, 26230; `branch_changed` empty among modules firm in pass A; G_r² − G_r¹ printed per run (expected ≈ 0.1 ns at X17); only 24503 gets `drift`.

### S4. Outputs: JSON, TSVs, validation text (0.25 day)

- JSON per run (schema §7), `json.dump(..., indent=1, allow_nan=False)`; `modules` sorted by module id; `offset_ns` rounded to 4 decimals (0.1 ps).
- `rf_calib_modules.tsv` (13 × 1152 rows), `rf_calib_runs.tsv` (13 rows), `rf_calib_validation.txt` (items (i)–(vi) of §4.5 with PASS/FAIL per gate, the (iv) board table and the (v) 33-module table), `run.log` via `tee`.

Acceptance: `python3 -c "import json; [json.load(open(f'output/rf_calib_20261026/{r}.json')) for r in (24386, 24503, 26230)]"` passes; `24386.json` has `"usable": false` and `"modules": []`; every other JSON has `default == 0.0`, ≥ 1100 entries in `modules` (measured + `pred`; 26138 ≈ 1150), all `offset_ns` in (−P/2, P/2], no `NaN` string in the file (`grep -c NaN` = 0); a loader-style check in Python (`{m["name"]: m["offset_ns"] for m in j["modules"]}`) resolves every listed name in `module_names`.

### S5. PDFs (0.5 day)

`rf_calib_<run>.pdf` (13 files, `MPLBACKEND=Agg`, rendered in a 3-process pool as `RF:`):
1. Text page: run, type, P, splits, tier counts, quality, split span, walk bias, the run's validation numbers.
2. Pooled phase, three panels: raw φ (no offsets), fold(φ − τ) (prior only), the held-out pooled residual of pass B (after); allE tier, 80 bins, Gaussian core overlaid.
3. Map page: the HyCal x–y map (`gp["x"], gp["y"], sx, sy` as `Rectangle` patches) coloured by ρ^r (before, range ±2.0 ns) and by dhalf (after, range ±0.3 ns).
4. Drift page (runs with `quality: "drift"` only): per-split centres of the pooled residual.
5. Module panels (4 per page, the full-statistics histogram unwrapped around the fit with the fit overlaid, as `RF:timing_panel` but one panel per run, with dhalf, e², D², tier and flags printed): the 33 modules of §2 in the three X17 runs, and every flagged module (`allE` is not a flag for this purpose) in every run.

Acceptance: 13 PDFs exist and open (page count printed by the renderer; `qpdf --check` if available); total render time < 5 min.

### S6. Full run and cross-check (0.25 day)

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh
MPLBACKEND=Agg python3 rf_time_calib.py --jobs 12 2>&1 | tee output/rf_calib_20261026/run.log
```

Acceptance: §4.5 gates (ii), (v), (vi) PASS in `rf_calib_validation.txt`, reports (i), (iii), (iv) present with their numbers; wall time < 8 min; `git status` shows no new entries beyond `rf_time_calib.py`, the S0 edit and the S7 documents; nothing under `database/`; no C++ file.

### S7. Documents (0.25 day)

- `notes/PbWO4_energy_resolution_20261025.md`: new §9.7 "RF time calibration per run (10-26)" after §9.6: method (condensed §3), the per-type σ and tier counts, the global 25800 → 26230 jump and the board table (iv), the 33-module table (v) with the finding that the central block is shifted in 26230 only and W525/W528/W597/W697 in all three X17 runs, 24386/24503 treatment, the (vi) numbers, output list; §10 rows "RF time calibration" → done with the caveats (time definition, template bank, allE walk) and the correction that the run-to-run phase change is global, not per board; §11(h) gets the S0 acceptance command and the S6 command.
- `notes/PbWO4_energy_resolution_summary.md`: §3.9 one finding row ("Per-run RF offsets derived …", with the 26230-only finding), §5.1 row "RF time calibration" → Done with the key numbers and `output/rf_calib_20261026/`; the §3.9 row "The FADC–RF phase changes per run" loses "partly per board".
- `plans/plan_20261025.md`: Todo 34 checked with a pointer to this plan and the output dir; one progress-log entry with the timestamp.
- This plan: Progress log and Todo synced.

Acceptance: every number in the documents is traceable to `rf_calib_validation.txt`, `rf_calib_runs.tsv` or `rf_calib_modules.tsv`; no number typed from memory.

## 6. Files to create or modify

| File | Change |
|---|---|
| `rf_time_calib.py` (new, ≈ 500 lines) | S1–S5 |
| `plot_module_timing_rf.py` | S0: `prep` returns `sel`; `pass1` deletes it (2 lines) |
| `output/rf_calib_20261026/` (new, outside git) | `<run>.json` ×13, `rf_calib_modules.tsv`, `rf_calib_runs.tsv`, `rf_calib_validation.txt`, `rf_calib_<run>.pdf` ×13, `run.log` |
| `notes/PbWO4_energy_resolution_20261025.md` | §9.7 new; §10 and §11(h) rows |
| `notes/PbWO4_energy_resolution_summary.md` | §3.9 rows, §5.1 row |
| `plans/plan_20261025.md` | Todo 34, progress log |
| `plans/plan_rf_calib_20261026.md` | this plan: log and Todo |

Not touched: `database/`, any C++ file, `module_timing_dump.py`, `rf_time_dump.py`, `plot_module_timing.py`, the `module_timing_20261025` outputs.

## 7. Interfaces and data

**JSON (`<run>.json`)**, key order as written:
```json
{
  "default": 0.0,
  "modules": [
    {"name": "W492", "offset_ns": -0.3831, "module_id": 1492, "offset_error_ns": 0.0031,
     "sigma_ns": 0.312, "n": 412345, "tier": "e100", "flags": ""}
  ],
  "run": 26230, "run_type": "x17", "period_ns": 4.00771, "usable": true, "quality": "ok",
  "split_span_ns": 0.042,
  "time_definition": "fold(t_module - t_rf - offset_ns) centred at 0 for pulses >= 100 MeV (tier allE: all pulses); t_module = matched-filter p_time of module_timing_dump 1.0 (20261025), templates <templates path from the dump>, no per-module time offsets applied; t_rf = nearest channel-A RF hit; fold(x) = x - P*round(x/P), P = period_ns",
  "created": "2026-10-26T..."
}
```
`quality` ∈ {`ok`, `drift`, `bimodal`}; `tier` ∈ {`e100`, `allE`, `pred`}; `flags` is a comma-joined subset of {`allE`, `pred`, `noprior`, `bimodal`, `nobranch`, `branch_weak`, `branch_changed`, `argmax`, `run_unusable`}. `pred` entries carry `offset_error_ns`, `sigma_ns` and `n` as `null`. For 24386: `"usable": false`, `"quality": "bimodal"`, `"modules": []`.

**`rf_calib_modules.tsv` columns:** `run type module id crate slot n_e100 n_allE tier offset_ns offset_err_ns sigma_ns tau_ns rho_ns rho_run_ns board_phase_ns board_mad_ns evt_ns dt_tcal_ns margin_ns dhalf_ns dhalf_err_ns flags` (`evt_ns` = e², `dt_tcal_ns` = D²; unresolved modules have `tier = unresolved` and `nan` values). NaN written as `nan`.

**`rf_calib_runs.tsv`:** `run type n_splits P_ns usable quality n_e100 n_allE n_pred n_unresolved walk_bias_ns split_span_ns frac_bimodal bimodal_ratio run_phase_ns G1_ns G2_ns G_mad_ns pooled_tcal_offset_ns rf_summary_run_offset_tcal_ns frac_dhalf_lt_0p3 median_abs_dhalf frac_rhorun_lt_0p3 median_sigma_e100 frac_sigma_lt_0p3`.

**Consumers.** Nothing in the repo reads the new files. The C++ loader would accept the JSON; the Python benchmark does not use RF offsets (out of scope).

## 8. Risks and compatibility

| Risk | Handling |
|---|---|
| Time definition: the offsets centre the matched-filter `p_time` (Python benchmark pipeline); `Replay.cpp` applies `hycal_rf_offsets` to `ClusterHit.time` (WaveAnalyzer CFD). Format-compatible, values not transferable to the C++ branch without a constant per module type | `time_definition` field in the JSON; stated in §9.7 and §10 of the note; adoption into `database/` is an open question |
| The 26230-only central-block shift may be a property of that run's data (rate, DAQ configuration) or of the `2p24GeV_Ta_X17` template bank interacting with it; the per-run offsets absorb it either way, but a bank change invalidates the X17 files | `templates` path recorded in `time_definition`; note caveat; Open Question 4 |
| Fold period: offsets are folded with the measured 4.00771 ns, the C++ with 4.008 ns; for \|t_module − t_RF\| ≤ 64 ns the difference is ≤ 5 ps | documented, no action |
| `allE` tier offsets are late by the walk (expected 0.1–0.3 ns, §9.4) relative to the plateau definition | walk bias per type printed and stored; flag `allE` on every such module |
| The continuous reference uses in-window rows only (X17: 71 % of dominant pulses in window) and is diluted by 0.2–0.3 ns when mates share the module's shift (diagnostic (c)); pass B corrects the mates once and shifts G_r by a run constant ≈ 0.1 ns, which is not an error of D (defined up to that constant) | the \|D² − e²\| < 0.5 ns gate and `branch_changed`; G1/G2 reported |
| `pred` tier offsets inherit the prior's error for that run (≈ 0.3 ns for most modules; wrong for a block like the 26230 central modules, which however always have N ≫ N_MIN) | `pred` flag; count per run reported |
| G_r assumes the FADC–RF phase is the same for all modules of a run up to the module's own continuous shift; the diagnostic confirms it on 25800, 26138, 26230 (MAD ≤ 0.09 ns); on the H2 runs it is checked by gate (vi) | gate (vi) |
| 24503 drift 0.23 ns within the run is not corrected by a per-run constant | `drift` quality, span in the JSON |
| `load_geometry` needs the C++ pipeline (`source setup_env.sh`), like every script of §9 | same environment as `RF:` |
| Memory: per worker one dump (≤ 3.4 M rows) plus H (2 × 2 × 1152 × 80 int32) and HE (2 × 1152 × 600 int32 ≈ 5.5 MB) per split; parent sums 111 returns | fine with 12 workers |
| Git: `plot_module_timing_rf.py` is untracked (new on 10-25); S0 modifies an untracked file | acceptance of S0 compares outputs, not git diffs |

## 9. Acceptance criteria and validation commands (consolidated)

```bash
cd /Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts && source setup_env.sh
SC=<session scratchpad>; O=output/rf_calib_20261026

# S0: RF: outputs unchanged by the prep/pass1 edit
MPLBACKEND=Agg python3 plot_module_timing_rf.py --runs 26138 --jobs 4 --out-dir $SC/rf_before   # before the edit
MPLBACKEND=Agg python3 plot_module_timing_rf.py --runs 26138 --jobs 4 --out-dir $SC/rf_after    # after
for f in module_timing_rf_test.tsv rf_summary_test.tsv module_timing_rf_vs_integral_test.tsv; do diff -q $SC/rf_before/$f $SC/rf_after/$f; done

# S1-S3 smoke (one X17 run, one 2.24 GeV run)
MPLBACKEND=Agg python3 rf_time_calib.py --runs 26138 25320 --jobs 8 --out-dir $SC/rf_smoke --no-pdf

# S6 full run
MPLBACKEND=Agg python3 rf_time_calib.py --jobs 12 2>&1 | tee $O/run.log
python3 - <<'EOF'
import json, glob
for f in sorted(glob.glob("output/rf_calib_20261026/*.json")):
    j = json.load(open(f)); P = j["period_ns"]
    assert j["default"] == 0.0 and all(-P/2 < m["offset_ns"] <= P/2 for m in j["modules"])
    print(f, j["usable"], j["quality"], len(j["modules"]))
EOF
grep -c NaN $O/*.json            # 0 for every file
cat $O/rf_calib_validation.txt   # items (i)-(vi), gates with PASS/FAIL
git status --short               # no database/ change, no C++ change
```

Per-step criteria are in §5. Global: gates (ii), (v), (vi) of §4.5 PASS; 24386 unusable, 12 runs usable; the 33-module table in `rf_calib_validation.txt` shows firm branches with \|D² − e²\| < 0.5 ns in the three X17 runs, and the values are compared with the §2 diagnostic expectations (deviations > 0.3 ns are reported, not gated).

## 10. Out of scope

Applying the offsets in the clusterer or benchmark; any C++ change; writing under `database/` or `database/runinfo`; a `time_calib.json` candidate for the X17 period (the TSV column `dt_tcal_ns` holds the numbers if wanted later); new evio processing; new X17 data; the 24386 DAQ question (Todo 35 of `plan_20261025`); the 131 → 128.25 ns doc fix (Todo 37).

## Open Questions

1. Adoption: once the files exist, which runs (if any) should be copied to `database/hycal_rf_offsets/` and wired in `database/runinfo/general.json` (`time_cuts.hycal_rf_offsets`)? Given the time-definition mismatch with `Replay.cpp` (§8, row 1), the user may prefer to keep them as Python-pipeline inputs only.
2. Should the continuous corrections D² for the X17 runs be turned into a per-period `time_calib.json` candidate (same schema as the current file)? Not planned here.
3. Is E ≥ 100 MeV the right plateau definition for the offsets, or should the threshold follow the cluster-seed energy of each run type?
4. The central-block shift of ≈ −1.6 ns exists in run 26230 only (25800 and 26138 agree with the prior), and W525, W528, W597, W697 are ≈ −1 ns early in all three X17 runs. Does the run log record a configuration change at 26230 (central-region FADC, trigger or threshold settings)? This is a question for the user / DAQ, not for the executor.

## Progress log

- 2026-10-26 — plan drafted; review round 1 (critic + scope-critic) applied: run-global reference phase replaces the per-board one (planner diagnostic on 25800, 26138, 26230), `branch_changed` defined on D, board tier dropped, outputs trimmed to JSON ×13, two TSVs, validation text, 13 PDFs. Round 2 (delta): per-run module counts corrected (26138: 235 / 1076), S2/S4 thresholds lowered, `pred` fallback tier added for modules below N_MIN, the \|D² − D¹\| < 0.10 ns clause dropped (pass B shifts G_r by a run constant), Φ_r decoupled from the pooled fit. Nothing executed.
- 2026-10-26 10:30 — **Executed S0–S7.** `rf_time_calib.py` full run 80 s (12 jobs); gates (ii), (v), (vi) PASS; 12 usable runs, 24386 `usable: false` (bimodal, 78.2 % of measured modules, ratio 0.35; empty `modules` in its JSON), 24503 `drift` (span 0.133 ns); no NaN in any JSON; 13 PDFs open. Outputs: 13 `<run>.json`, `rf_calib_modules.tsv`, `rf_calib_runs.tsv`, `rf_calib_validation.txt`, 13 `rf_calib_<run>.pdf`, `run.log`. Documents: note 20261025 §9.7 new, §10 row, §11(h); summary §3.9 row and §5.1 row; `plan_20261025` Todo 34 and log.
  - **Results.** (ii) median \|dhalf\| 0.0167 ns, 0.9985 below 0.30 ns over 13 011 of 13 254 measured (run, module) pairs (dhalf defined when both halves ≥ 150 rows; counting the undefined as failing: 0.9802); 24503 alone 0.0138 / 0.9981. Median e100 σ 0.183–0.203 ns by type. (iv) global jump +1.990 ns. (v) all 33 modules firm in 25320, 25800, 26138, 26230, max \|D² − e²\| 0.065–0.089 ns; the central block is shifted in 26230 only (D² −2.10 … −0.34), W525, W528, W597, W697 are early in all three X17 runs (§9.7 table). (vi) G_mad 0.020–0.040 (pass B), pooled tcal fit identical to `rf_summary.tsv` in all 13 runs.
  - **Deviations from the expectations.** (1) S2 acceptance: run 24512 has 999 measured modules (< 1000); `pred` counts 94 / 151 / 128 at 24503 / 24512 / 24603 (the plan expected "a handful"). (2) (iii): \|ρ^r\| < 0.3 for 0.876–0.945 per usable run (expected ≳ 0.95); e100 tier alone 0.969–0.984 (26230: 0.894). (3) (iv): boards outside ±0.1 ns of the jump: (3,15) +1.578, (3,8) +0.322 as expected, and (3,18) +0.143 (not expected). (4) (v): at the 2.24 GeV runs \|D²\| ≤ 0.4 for 24 of 33 (expected 33), none beyond 0.7. (5) (ii) gate denominator is the set with dhalf defined, 13 011 of 13 254; the plan's 'of measured' wording gives 0.9802 < 0.99, open for the user. (6) `m_fb` rows exist: 660 in the full set (§2 said none), dropped from the histograms. (7) 26230 D² upper end −0.34 against the expected −0.4 (0.06 ns outside).

## Todo

- [x] 1. S0: `RF:prep` returns `sel`; `RF:pass1` deletes it; before/after `--runs 26138` outputs identical.
- [x] 2. S1: `rf_time_calib.py` skeleton, CLI, geometry, prior, board map, `passA`, per-run aggregation; smoke on 26138.
- [x] 3. S2: fits per (run, module, tier, half), tier choice incl. `pred`, Φ_r, Δ_b, bimodal detector, walk bias; 24386 flagged, others not.
- [x] 4. S3: e¹, G_r, D¹; `passB`; e², D², flags; drift span; 33-module gate.
- [x] 5. S4: JSON ×13, two TSVs, validation text; JSON loads, no NaN, offsets in range.
- [x] 6. S5: 13 per-run PDFs.
- [x] 7. S6: full run with `tee`, gates (ii), (v), (vi) PASS, cross-check vs `rf_summary.tsv`.
- [x] 8. S7: note §9.7, §10, §11(h); summary §3.9, §5.1; `plan_20261025` Todo 34 and log; this plan's log and Todo.
