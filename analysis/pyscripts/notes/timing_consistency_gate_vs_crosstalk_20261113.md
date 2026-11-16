# Timing consistency: time-gated clustering vs the crosstalk study

**Date:** 2026-11-13. **Question (user):** the seed vs non-seed Δt in the crosstalk study is very narrow; are the times computed and corrected the same way as in the time-gated clustering? **Scope:** code tracing plus one read-only recount of the cached crosstalk samples; no code changed. Paths are relative to the repo root (`prad2evviewer/`); `ERB` = `analysis/pyscripts/waveform_analysis/energy_resolution_benchmark.py`, `FB` = `analysis/pyscripts/benchmark_hycal_timing_test_20261021.py` (frozen benchmark).

## Answer

The time is the same quantity, and the difference is the correction. Both use the matched-filter (MF) 50 % time t = onset + t50 of the run type's template, and neither corrects for time walk. The crosstalk study subtracts the per-module offsets of `database/hycal_time_offsets/time_calib.json` from the seed and member times. The gate in `ERB` uses raw times. The seed-member offset difference has an RMS of 2.5 ns, which widens the raw Δt core from about 1 ns to about 2.8 ns (robust σ, 2.24 GeV). Inferred: that accounts for the steep σ/E gain up to about 8 ns. The approach from 8 to about 11 ns is a late tail that the offsets do not touch (89-94 % of the |Δt| > 8 ns members are late at 2.24 and 3.49 GeV; at 0.73 GeV only 38 % are, so there the tail beyond 8 ns is mostly offsets).

## 1. Time in the gated clustering (A)

- Per pulse: `t = t0 + t50[mi]` (`ERB`:879). t0 is the MF onset (parabolic vertex on the 1 ns grid, `ERB`:318-326). t50 is the template's 50 % rise time per module (`template_columns`, `ERB`:705-715). Units are ns from the start of the waveform. Lead glass and clipped PbWO4 channels without an MF pulse use the WaveAnalyzer CFD pick-off `pk.time` instead (`ERB`:860-864; `prad2dec/src/WaveAnalyzer.cpp`:572). The crosstalk note puts that basis about 0.6 ns earlier than MF.
- No offsets: `PRE_LO, PRE_HI = 100.0, 200.0 # ... (no time offsets)` (`ERB`:79). Neither `time_calib`, `t_off` nor `time_offset` appears anywhere in `ERB`.
- Pulse selection: the scan stores every in-window pulse, and `extra` marks every pulse but the largest on a module (`ERB`:924). Legacy feeds only non-extra pulses. Gated paths feed all of them (`ERB`:1074-1075, `cl.add_hit(mod[i], e[i], t[i])`). C++ then takes, per neighbour module, the largest-energy unconsumed pulse with |t − t_seed| ≤ gate (`prad2det/src/HyCalCluster.cpp`:198-209). Every member is compared with the seed's own raw time (`:186`), not with its BFS parent. The same gate is applied again in `calculate_energy_square` (`:336-337`).
- No time-walk correction and no seed-time correction anywhere. `ModuleHit::time` is documented only as "ADC peaking time (ns)" (`prad2det/include/HyCalCluster.h`:85).
- Gate set from python: `cfg.seed_time_window = window` in `clusterers()` (`ERB`:963). `ERB` uses its own python MF copy. The C++ side only clusters the times it is given.
- The C++ builder loads `time_calib.json` into `Module::time_offset` (`prad2det/src/PipelineBuilder.cpp`:310-341; `database/reconstruction_config.json`:13-14, `time_calib: true`). `HyCalCluster` never reads that field: outside the loader, `time_offset` appears only in `analysis/tools/*`, and `python/` has no binding for it. The C++ convention is that the caller subtracts the offset before `AddHit` (`analysis/tools/physics_calib.cpp`:603, `hycal_shower_profile.cpp`:187-192).

## 2. Time in the crosstalk and time studies (B)

- Previous crosstalk study (2026-11-06/07, the source of the narrow Δt): the default time is the WaveAnalyzer log-normal CFD (`--time-source` default `cfd`, :541); the MF variant is the same frozen-benchmark onset + t50 as (A) (`FB`:1693). Offsets come from `time_calib.json` with `t = t − offset` on seed and member (`crosstalk_study/previous_20261107/crosstalk_study.py`:21, 320-331; default `--tcal auto`, :550). The frozen-benchmark runs themselves were raw (`time_study/README.md`:3: "time calib : off").
- Members come from the legacy path (`--path legacy` default, :538), so the narrow Δt is not a selection effect of a gate.
- Reported widths (CFD, the default; the MF equivalent is 1.05 ns robust, §3): calibration narrows the core Gauss σ from 2.37 to 0.68 ns (2.24 GeV e-p; `crosstalk_study/notes/crosstalk_findings_20261107.md`:25). Self-referenced offsets reproduce 0.676 ns (:78). An earlier note gives seed-neighbour σ_MAD 3.10 → 1.05 ns with the offsets (`notes/PbWO4_energy_resolution_20261022.md`:23).
- No walk correction. In the time study, an in-sample walk correction gains ≤ 0.01 ns for seeds ≥ 500 ADC (`time_study/README.md`:28).
- The current `crosstalk_study.py` computes no Δt to the seed and uses no time calibration (`crosstalk_study/README.md`, "Dropped, and why").
- `pileup_rates.py` uses the raw MF onset (not onset + t50) in a 130-200 ns window, with offsets off and no time coincidence for members (`waveform_analysis/pileup_analysis/pileup_rates.py`:38-41, 125).

## 3. Recount on the cached crosstalk samples (verified, read-only)

Source: `crosstalk_study/cache/<rt>/*_legacy.npz`, first 3 splits, all legacy non-seed members. MF times are raw, and the offsets come from `time_calib.json`. Script: scratch `sup-dt.py` (not kept). Robust σ is half the 16-84 % range.

| run type | N members | raw Δt robust σ | calibrated | offset diff robust σ / RMS | members with abs(Δt) > 8 ns, raw / cal | non-seed energy with abs(Δt) > 4, 6, 8, 10 ns, raw | same, calibrated |
|---|---|---|---|---|---|---|---|
| 0.73 GeV H2 | 2.53 M | 2.30 ns | 0.69 ns | 2.16 / 2.45 ns | 0.50 / 0.09 % | n/a | n/a |
| 2.24 GeV H2 | 4.12 M | 2.79 ns | 1.05 ns | 2.20 / 2.50 ns | 2.48 / 1.87 % | 13.8, 4.3, 0.55, 0.17 % | 1.3, 0.66, 0.30, 0.15 % |
| 3.49 GeV H2 | 2.70 M | 3.11 ns | 1.43 ns | 2.26 / 2.53 ns | 4.07 / 3.14 % | 14.5, 5.0, 0.64, 0.18 % | 1.8, 0.96, 0.41, 0.13 % |

The seed-member offset difference exceeds 8 ns for only 0.04-0.09 % of members. Over the whole table, the offsets have σ = 2.69 ns (range −8.9 to +7.0 ns, 1150 resolved). Adjacent 3 × 3 pairs differ by a median of 1.24 ns (p90 3.9 ns, p99 7.4 ns).

## 4. Differences between (A) and (B), and which one sets the 11 ns plateau

| aspect | (A) gate in `ERB` | (B) previous crosstalk study |
|---|---|---|
| time quantity | MF onset + t50 (WaveAnalyzer CFD for LG and clipped fallback) | WaveAnalyzer log-normal CFD (default), or the same MF onset + t50 |
| per-module offset | **none** | `time_calib.json` subtracted |
| time walk | none | none |
| pulse per module | largest within ±gate of seed (gated); largest in window (legacy) | largest in window (legacy members) |
| reference | seed's raw time, every member | seed's calibrated time |

- **Verified:** the offsets are the only correction that differs. They widen the core by a factor of 2.2-3.3, and the raw curve loses 4.3 % of non-seed energy at 6 ns against 0.66 % calibrated (2.24 GeV).
- **Inference:** the steep part of the σ/E scan (2-8 ns) is the uncalibrated offset spread. The 8 → 11 ns approach is the late tail, which offsets do not move: above 10 ns the raw and calibrated losses agree (0.17 vs 0.15 %). This assumes σ/E tracks the lost non-seed energy.
- **Expected with offsets (inference, e) below):** the energy loss at a 5-6 ns calibrated gate matches the raw loss at about 8 ns, so the steep rise should end near 5-6 ns. The final plateau should stay near 10-11 ns because of the late tail. A symmetric gate is wider than the core requires; the tail is late-sided.
- **Prior measurement (verified, WaveAnalyzer `pk.time` basis, not MF):** a raw vs calibrated gate scan found that calibrated times help only at gates ≤ 8 ns, that raw and calibrated agree within 0.002 pp at 16 ns, and that tight calibrated gates are still slightly worse than no gate (`tcal_w2` +0.14, `tcal_w4` +0.03 pp; `notes/PbWO4_energy_resolution_20261022.md`:24, 30). This agrees with the expectation above.

## 5. Applying the offset before the gate scan (binding user decision)

**(a) Where.** Use the python driver, in the report step only, with no rescan: `h_t` is stored raw (`ERB`:924-936). At `ERB`:1075 pass `t[i] - toff[mod[i]]`, with `toff` indexed by `Module.index` and built once per run in `clusterers()` (`ERB`:949-972) or `split_job`. Port `load_time_offsets` from `FB`:1125-1139 (uses `hycal.module_by_id`, which is bound, `python/bind_det.cpp`:673). This follows the C++ convention (caller subtracts, §1). Do not change `HyCalCluster`: its other callers already pass calibrated times (`physics_calib.cpp`:603), and reading `Module::time_offset` there would subtract twice. Same-module differences in `secondary()` (`ERB`:1033-1044) cancel the offset. The 100-200 ns pre-window (`ERB`:79) is far wider than the ±9 ns offsets and can stay raw. Record the table path in the outputs as `FB`:3800 does. The PbWO4 offsets are not centred (mean −1.42 ns, median −1.62 ns; reference W495) and lead glass stays at 0, so calibration shifts PbWO4 by about +1.4 ns relative to lead glass. Transition-region clusters and the LG/T region rows are therefore not comparable between the raw and calibrated scans.

**(b) Table.**
- Production table: `database/hycal_time_offsets/time_calib.json`, reachable through `$PRAD2_DATABASE_DIR`.
- Header: `reference_module` 1495 (`W495`, offset 0), `offset_units` ns, `default` 0.0.
- Rows: `name`, `module_id`, `offset_ns`, `offset_error_ns`, `connected`, `resolved`. 1152 PbWO4 rows, 1150 resolved; there are no lead-glass rows, so lead glass stays at 0.
- Sign convention: `calib_time = raw_time − offset_ns` (`prad2det/include/HyCalTimeCalib.h`:5-7). Use only `resolved` rows, as `FB` does.
- The `crosstalk_study/tcal/*_selfref_offsets.tsv` files are a cross-check, not production. Their columns are `module n offset_ns se`, they hold 0 for unconstrained modules, and the findings note says "The official table is kept" (:78).

**(c) Time walk.** None in the table: it holds a constant `offset_ns` per module. Walk is not needed for the core (≤ 0.01 ns at seed ≥ 500 ADC, `time_study/README.md`:28). For small members (late-tail medians 7-12 ADC) it has not been measured; that is a gap, not a reason to add a correction now.

**(d) Coverage.**
- One run-independent table. It was committed on 2026-08-01 (commit `aaeaaa5`) and derived with `analysis/tools/hycal_module_time_calib.cpp`, a neighbour-pair Δt network on WaveAnalyzer peak times. Its source runs are not recorded in the file. The MF time is a cross-basis use of it, though MF self-referenced offsets agree on pairs within RMS 0.41-0.46 ns (findings :78).
- **Verified for H2 0.73, 2.24 and 3.49 GeV:** it narrows the core 2.3-3.4× in all three (§3).
- **Ta X17 (25800, 26138, 26230), inference:** the calibrated core |Δt| < 2 ns holds 69-83 % of hits there against 78-95 % for H2 (`crosstalk_study/tables/crosstalk_ratio.md`). No Ta raw-vs-calibrated width was recounted here. Raw and calibrated Ta jobs exist (`analysis/pyscripts/Commands.md`:198).
- No per-period re-derivation appears necessary for H2. Run-to-run drift is ≤ 0.1 ns for the veto (`notes/PbWO4_energy_resolution_20261025.md`:682).

**(e) Expected widths and rescan range.** Calibrated seed-member Δt robust σ: 0.7 ns at 0.73 GeV, 1.05 ns at 2.24 GeV, 1.4 ns at 3.49 GeV (Gauss core 0.68 ns at 2.24 GeV). Suggested rescan: 1, 2, 3, 4, 5, 6, 8, 10, 12, 16 ns. Keep the raw scan as reference, and dense below 6 ns where the curve should now bend.

## Risks and gaps

- §3 uses the frozen-benchmark MF cache (older MF version, 3 splits, no e-p cut). The `ERB` MF (`MIN_SEP_NS` 8) may differ slightly.
- The link from lost energy to the σ/E plateau is inferred, not measured. The late tail's origin (crosstalk shoulder, walk of small pulses, or physics) is open.
- Clipped seeds fall back to the WaveAnalyzer time, a basis about 0.6 ns off MF. Their pick-off on clipped samples has not been validated.
