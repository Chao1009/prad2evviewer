# Research: LMS gain tracking for the HyCal waveform timing benchmark

**Date:** 2026-10-21  
**Scope:** What already exists in `prad2evviewer` (C++, Python, config, docs; build dirs and `/Users/L00338853/work/PRad/old` excluded) for measuring HyCal gain drift from LMS events in time bins and applying it as a correction. Also: how LMS events, reference PMTs and timestamps are exposed, how `analysis/pyscripts/benchmark_hycal_timing.py` treats LMS events, which e-p kinematics, geometry and resolution helpers exist for a 2-sigma elastic cut, and what LMS state the FADC calibration corresponds to.  
**Out of scope:** designing or writing `benchmark_hycal_timing_test_20261021.py`.

All paths are relative to `/Users/L00338853/work/PRad/prad2evviewer` unless absolute. Line numbers are as of commit `4379db3` plus the working tree.

---

## 0. Bottom line

- **Q1 verdict: partial, and usable for the C++ replay path. Nothing exists in Python.**
  - The C++ chain `replay_gainCorr` → `GainCorrCompute` already does what the user describes, with one difference: it bins by a fixed number of LMS events, not by fixed time. Each bin gets a per-module, reference-PMT-normalised LMS gain correction, written as a ROOT time series that the C++ replay applies event by event.
  - A file for the user's run already exists: `database/gain_factor/gain_correction/prad_025320_gain_corr.root`. Its bins are about 7.5 minutes long (4000 LMS events each).
  - Re-running with `-b ~1600` gives about 3-minute bins (inference from the measured LMS rate).
  - No Python script computes or applies LMS corrections. `scripts/online_gain_monitor.py` reads the output `gain_corr.root` with uproot, which shows the format can be read from Python.
- The "reference PMT normalisation" in the existing code is **not** a plain division by an LMS reference PMT. Each reference PMT (LMS1/2/3) is itself normalised by the Am-241 alpha source seen on the same PMT. That removes drift of the LMS lamp and of the reference PMT gain. The correction is then taken relative to a **reference run** (`ref_run`), not to the first bin of the current run.
- `benchmark_hycal_timing.py` drops LMS events. It requires `trigger_bits == 0x100` exactly, skips non-HyCal channels (so it never sees LMS1-3), and never reads timestamps.

---

## 1. Q1: Existing LMS handling

### 1.1 Verified facts

| Component | Path | What it does | Input → output | Waveform? | Time binning | Status |
|---|---|---|---|---|---|---|
| `replay_gainCorr` (phase 1: `Replay::Process_LMSgainFactor`) | `analysis/tools/replay_gainCorr.cpp`; `analysis/src/Replay.cpp:1810`, selector at `Replay.cpp:1981-1982` | Runs the soft `WaveAnalyzer` on every FADC channel of LMS and alpha events. Writes per-channel peak integrals and the TI timestamp. | EVIO → `*_lms.root`, tree `lms_gain` (`event_type` 0=LMS 1=alpha, `module_id`, `npeaks`, `peak_integral`, `timestamp`), plus `epics` and `scalers` trees | Yes, FADC waveform integrals | per event | Built (`analysis/CMakeLists.txt:21,54`) |
| `ComputeGainCorrections` (phase 2) | `analysis/src/GainCorrCompute.cpp` (cuts `:302-304`, batch close `:330`, unix time `:75-110`, batch event range `:332-333`) | Groups LMS events into batches of `batch_size`. Per batch it fits each PbWO4 module's LMS peak, each reference's LMS peak, and each reference's alpha peak. It forms ratios and compares them with the reference run. | `*_lms.root` + ref `.dat` → `prad_NNNNNN_gain_corr.root`, tree `gain_corr` | Yes | **N LMS events per batch**; `-b` default 4000 (`replay_gainCorr.cpp:202,222`) | Built. Usable. |
| `GainCorrTimeSeries` / `LoadGainCorrTimeSeries` | `analysis/include/gain_factor.h:206-226` (`GetCorr`), `:279-317` (loader) | Step-function lookup by event number: the last batch with `event_num_start <= event_num`, with no interpolation. Events before batch 0 get batch 0. The trailing partial batch is dropped (no flush after the loop, `GainCorrCompute.cpp:330-357`), so the tail of a run uses the last full batch. | `gain_corr.root` → per-module factor | n/a | per batch | Used by replay |
| Application in replay | `analysis/src/Replay.cpp:474,492,527-529`; same at `:885,921,1003-1004` and `:1402,1459,1562-1563` | `gain_factor = (corr[LMS2] + corr[LMS3]) / 2` for PbWO4; PbGlass gets `g[id].avg` (never filled, so 1.0). `adc = peak_integral * gain_factor`, then `energize()` applies the calibration. | `hycal.gain_factor` branch (`docs/REPLAYED_DATA.md:78`) | Yes | per batch | Usable. PbWO4 only. |
| `refGain_produce` | `analysis/tools/refGain_produce.cpp` (selector `:239-247`, `g[j]` `:369-373`, writer `:337-385`) | One global fit per run, producing the reference table. | EVIO → `database/gain_factor/ref_gain/prad_NNNNNN_LMS.dat` (+ `_hists.root`) | Yes | whole run | Built. 11 ref `.dat` files exist (23915, 24512, 24827, 25293, 25546, 25547, 25560, 25655, 25825, 26094, 26214). |
| `online_gain_monitor_base` + `scripts/online_gain_monitor.py` | `analysis/tools/online_gain_monitor_base.cpp`; `scripts/online_gain_monitor.py:1-10,117,213-256` | Incremental wrapper around the same `ComputeGainCorrections`. The PyQt GUI reads `gain_corr` with **uproot** (`gain_corr_W`, `gain_W`, `unix_time`). | EVIO splits → `gain_corr.root` → plots | Yes | batch; GUI default 1000 events, "about 100 s at 10 Hz" (`docs/analysis_notes/online_gain_monitor.md:238`) | Display tool. Its uproot reader is a working Python template. |
| Server LMS monitor + `resources/lms.js` | `src/app_state.cpp:526-562,648-665`; `src/app_state_api.cpp:297`; `database/monitor_config.json:65-110` | Live drift display. Drift = current / (baseline × lamp_scale), with lamp_scale taken from LMS2. LMS time cut is 160-220 ns. | live events → web UI | Yes | rolling | Display only. Writes no output file. |
| `scripts/hycal_gain_monitor.py` | `scripts/hycal_gain_monitor.py:1-13,171-216` | GUI showing the run-by-run `.dat` gain factors. | `.dat` files → plots | n/a | per run | Display only. Column-order bug, see 1.3. |
| Legacy | `analysis/tools/gain_monitor.cpp` (`.dat` writer commented out at `:310`), `analysis/tools/gain_fitter.cpp`, `analysis/scripts/lms_alpha_normalize.C` (per-event alpha-normalised graphs, hardcoded `/data/prad/run1234` example; `analysis/README.md:506-531`) | Older per-run or per-event studies | — | Yes | none or per run | Superseded |

**Python:** the only LMS hits in `analysis/pyscripts/*.py` and `calibration/*.py` are module-type labels and colours, plus `fit_pulse_template.py --trigger-event-type LMS`. Nothing computes LMS amplitudes, ratios or drift (`git ls-files '*.py' | xargs grep -il lms`).

### 1.2 Definition of the existing correction (verified, `GainCorrCompute.cpp`, branch names checked with ROOT)

The quantities below are computed for each batch b, each PbWO4 module i, and each reference PMT j = LMS1..3.

```
refPMT_ratio[j]    = fit_mean_ref_lms[j] / fit_mean_ref_alpha[j]        # lamp intensity seen by ref j
gain_W[i][j]       = fit_mean_W_lms[i] / refPMT_ratio[j]
gain_W_ref[i][j]   = (same quantity from ref_run .dat: W_peak * alpha_ref[j] / lms_ref[j])
gain_corr_W[i][j]  = gain_W_ref[i][j] / gain_W[i][j]                    # 1.0 if fit failed
applied factor     = (gain_corr_W[i][LMS2] + gain_corr_W[i][LMS3]) / 2  # LMS1 unused (physics_calib.cpp:792)
```

Fits use `gain_hist_fitter` (Gaussian over the region above 10% of the maximum; 600 bins, 0-15000; `analysis/include/gain_factor.h`). Channel cuts are `npeaks == 1` and `peak_integral > 0`.

### 1.3 Content of the existing file for run 025320 (verified, ROOT scan)

`database/gain_factor/gain_correction/prad_025320_gain_corr.root`, tree `gain_corr`, 9 entries:

| batch | event_num range | n_lms | n_alpha | unix_time |
|---|---|---|---|---|
| 0 | 1222 – 9,381,003 | 4000 | 6968 | 1781983239 |
| 1 | 9,387,322 – 17,750,574 | 4000 | 6905 | 1781983691 |
| ... | ... | 4000 | ~7000 | steps of ~450 s |
| 8 | 62,040,879 – 70,146,636 | 4000 | 6988 | 1781986817 |

- `ref_run = 24827`. That is the 3.5 GeV reference, inherited by run 25320 through `database/runinfo/general.json`; see Q7.
- Batch 0 `gain_corr_W` over all 1155 non-zero PbWO4 modules ranges from 0.876 to 1.157 (mean 0.980). 1023 of the 1155 lie in [0.94, 1.00]. The gains therefore differ from run 24827 by a few percent.
- Drift within the run: the applied factor (LMS2+LMS3)/2 changes by at most 1.2% relative to batch 0 across batches 1-8 (about 60 min).
- W1141 is 0.000 in all three columns in every batch. It is absent from `prad_024827_LMS.dat` and listed as dead from run 24850 (`general.json:169`). The C++ loader treats 0 as identity, 1.0 (`gain_factor.h`, `LoadGainCorrTimeSeries`). A raw uproot reader must do the same.
- Implied LMS rate: 4000 / ~450 s ≈ **8.9 Hz**, or about 1600 LMS events per 3 minutes.
- 21 runs have `gain_corr.root` files: 24586, 24827, 24917-24963 (13 files), 25308, 25318, 25319, 25320, 25324, 25326. No X17-period run has one.

### 1.4 Defects found (verified unless marked)

- **Column-order mismatch in display readers.**
  - The writer (`refGain_produce.cpp:343-361`) emits `Name lms_peak lms_sigma lms_chi2 alpha_peak alpha_sigma alpha_chi2` for LMS rows. All 11 current `.dat` files are in this LMS-first order.
  - Both `scripts/hycal_gain_monitor.py:195-201` and `src/app_state_init.cpp:19-52` (`parse_lms_dat`) read the LMS rows with the legacy alpha-first order of `gain_fitter.cpp:80-90`, so those two displays show LMS and alpha swapped.
  - The replay loader `LoadRefGainFile` (`gain_factor.h:89-139`) reads only g1-g3 from the W/G rows, so the correction path is unaffected.
- The `FindGainCorrRootFile` header comment says "largest run <= N". The code matches the exact run only.

### 1.5 Inferences

- The phase 1 integral is a soft-`WaveAnalyzer` `peak_integral` on the same `daq_config.json` `wave_cfg` that the Python pipeline uses (`_common.py:~232`). The factors are therefore probably on the same integral definition as `benchmark_hycal_timing.py`'s `pk.integral`. I did not confirm that the integration settings are identical in both paths.
- Batches are keyed by `event_num`, which comes from `info.event_number` (`Replay.cpp:474`). The Python pipeline has the same field, so the existing `GetCorr` lookup can be reproduced in Python without timestamps.

---

## 2. Q2: How LMS events are identified

### Verified

- `database/trigger_bits.json`:
  - FP bit 24 = LMS, bit 25 = Alpha, bit 15 = Pulser.
  - `trigger_type` 0x39 = LMS, 0x3A = Alpha.
  - `prad2det/include/EventData.h:53-54`: `TBIT_lms = 1u<<24`, `TBIT_alpha = 1u<<25`.
- Python access is through `fadc_evt.info.trigger_bits`, a uint32 FP snapshot, and `fadc_evt.info.trigger_type`, a uint8 TI event type. See `prad2dec/include/Fadc250Data.h:82-87` and `python/bind_dec.cpp:143-150`.
- The **existing C++ selector combines the trigger bit with multiplicity**. Quoted below from `refGain_produce.cpp:239-247`. `Replay.cpp:1981-1982` has the same LMS cut, with the alpha cut written as `nch < 50` (no `>= 1`).
  ```
  // Classify by channel count (trigger_bits unreliable):
  //   LMS flash illuminates all ~1735 channels -> nch > 1000
  //   Alpha source fires only a few crystals  -> 1 <= nch < 50
  const bool is_lms   = (nch > 1000) && trig_lms;
  const bool is_alpha = (nch >= 1 && nch < 50) && trig_alpha;
  ```
  Here `nch` counts readout channels with `nsamples > 0` (`refGain_produce.cpp:227-236`). It is not a software threshold.
  - `gain_monitor.cpp:227-291` uses the bit plus `nIntegrals > 500`.
  - The online monitor uses the bit only (`monitor_config.json` `lms_monitor.accept_trigger_bits: ["LMS"]`).
- Module counts in `database/hycal_map.json` (1735 entries): PbWO4 1152, PbGlass 576, Veto 4, LMS 3. HyCal proper is 1728.
- `analysis/README.md:513-516`: LMS lights "All HyCal modules + LMS 1/2/3". Alpha (Am-241) lights "LMS 1/2/3 only".
- There is an alpha trigger: bit 25 / 0x3A, prescaled at the TS (factor 16385 per `prad_v0.trg` / `banktags.md`). Run 25320 has about 7000 alpha events per 4000-LMS batch (section 1.3).
- **Surprise.** `docs/technical_notes/waveform_analysis/fit_pulse_template.md:38,42` reports that `--trigger-event-type LMS` (a `trigger_type == 0x39` filter) found **no** LMS events in runs 025308 and 025320. Yet `prad_025320_gain_corr.root` holds 36,000 LMS events selected by bit 24 plus `nch > 1000`.

### Inferences

- ">1000 modules fired" is consistent with the code base as a confirmation cut. It is weaker as the only selector.
  - Large SSP raw-sum events can also carry many channels: event sizes up to 71k words versus 25-55k for LMS (`docs/rols/banktags.md:101`).
  - The C++ code always ANDs it with bit 24.
- The template-scan result suggests that in the 2.2 GeV period LMS events do not carry TI `trigger_type` 0x39, even though FP bit 24 is set. That would match the "trigger_bits unreliable" comment.
- Safest selector: `(trigger_bits & (1<<24)) and nch > 1000`, the same as the existing C++. This has not been checked on raw EVIO from 25320.

---

## 3. Q3: Reference PMTs

### Verified

- There are three reference PMTs, LMS1, LMS2 and LMS3. In `database/hycal_map.json:1730-1732` they are crate 3, slot 13, channels 10, 11 and 12, `"t":"LMS"`. In the replay they are module ids 3101-3103 (`refGain_produce.cpp:58`).
- An "LMSPin" channel (id 3100) is referenced (`analysis/include/Replay.h:59`, `Replay.cpp:124`, `EventData.h:72,88`) but is **not mapped** in the current `hycal_map.json`.
- The older maps in `database/history/daq_map_*.json`, used up to run 23849, place LMSP and LMS1-3 in different slots.
- The replay uses LMS2 and LMS3 only. LMS1 is computed but not applied (`physics_calib.cpp:792`). The online monitor uses LMS2 as `drift_ref_channel`.

### Inference

The current map applies to run 25320, which is after 23849.

---

## 4. Q4: Timestamps and LMS rate

### Verified

- `fadc_evt.info.timestamp` is a uint64 holding the 48-bit TI timestamp in 250 MHz ticks, 4 ns each (`Fadc250Data.h:87`; assembled at `prad2dec/src/EvChannel.cpp:489`). It is exposed in Python (`python/bind_dec.cpp:148`).
- Also exposed: `info.event_number`, `info.unix_time` (from the run-info bank, 0 if absent), and the SYNC `unix_time`.
- The C++ tick constant is `TI_TICK_SEC = 4e-9` (`src/viewer_utils.h:17-30`, `analysis/tools/replay_filter.cpp:938`).
- No Python script reads `info.timestamp` today. `benchmark_hycal_timing.py` reads neither the timestamp nor `event_number`.
- Measured LMS rate in run 25320: about 8.9 Hz, from 4000 LMS events per about 450 s (section 1.3). The nominal rate in the docs is 10 Hz (`docs/analysis_notes/online_gain_monitor.md:238,322`). LMS is not prescaled at the TS (`prad_v0.trg:125`, `TS_FP_PRESCALE 24 0`).

### Inference

- A 3-minute bin holds about 1600 LMS events and about 2800 alpha events.
- Elapsed time in seconds is `(timestamp - t_first) * 4e-9`. The 48-bit counter wraps only after about 13 days.

---

## 5. Q5: LMS events in `benchmark_hycal_timing.py`

### Verified

- **Selection.** `benchmark_hycal_timing.py:369-376` skips non-Physics records, then `if fadc_evt.info.trigger_bits != C.PHYSICS_TRIGGER_BITS: continue` with `PHYSICS_TRIGGER_BITS = 0x100` (`_common.py:414`). This exact-equality test drops LMS, alpha, pulser and SSP1-3 events.
  - `_common.py:402-408` provides the recommended `passes_physics_trigger()` (`bits & 0xf00`). Its docstring says to use it "instead of `== PHYSICS_TRIGGER_BITS`".
  - At about line 400, `if mod is None or not mod.is_hycal(): continue` also skips the LMS1-3 channels.
- **Gates that the elastic cut would replace.** Refinement events must have exactly one cluster with `r.hit.energy > --signal-min` (default 3000 MeV; `:249,424-425`). Gains are scaled to a fixed `--target-peak` (default 3488.5; `:262,196`), not to E'(θ). `Commands.md` runs 25320 with `--signal-min 500 --target-peak 2239`.
- **Integration window.** The peak-time gate is `--pre-window` (130, 200) ns (`:245`), applied at `:404` as `if pk.time <= pre_lo or pk.time >= pre_hi`. The integral itself is the `WaveAnalyzer` peak integral with the `int_tail_ratio` rule (`database/daq_config.json` `wave_cfg`, 4 ns samples). It is not a fixed box.
- **LMS pulse time.** The repo's LMS gate is 160-220 ns (`monitor_config.json:65-71` `lms_monitor.time_cut`, applied at `src/app_state.cpp:558-562`). Inside the benchmark gate, only 160-200 ns overlaps it.
  - `database/waveform/pulse_templates_024177.json` has an LMS block with LMS1-3 (585 pulses). Its `t0_ns` values are relative to the fit window, not absolute.
  - No measured per-channel LMS peak time in absolute ns exists in the repo.
- **Seed-module refinement.**
  - `fit_per_module_gains()` is at `:152-205` and is called at about `:503-516`.
  - It uses the per-event tuples `(seed_id, x, y, row, col, seed_E, cluster_E)` collected at about `:453-472` when `--calibrate` is set. Inner radius is 200 mm, seed E is at least 200 MeV, there are 2 iterations of `gain *= target/ref`, and at least 30 events per module are required.
  - The calibration comes in through `setup_pipeline(..., hc_calib_file=...)` (`_common.py:204-250`). Per-hit energy is `mod.energize(pk.integral)` (about `:410-419`).
- **Possible insertion points.**
  1. Scale `pk.integral` per module and per time bin before `energize`. This is the same place the C++ replay applies `gain_factor`.
  2. Store `event_number` or the timestamp in the tuples and apply the factor before `fit_per_module_gains`.

### Inference

If the user's script integrates LMS pulses in the benchmark's 130-200 ns gate, LMS peaks above 200 ns will be lost. The existing C++ LMS path uses 160-220 ns.

---

## 6. Q6: Elastic e-p kinematics, geometry, resolution

### Verified

- `PhysicsTools::ExpectedEnergy(float theta_deg, float Ebeam, const std::string& type)` is in `analysis/include/PhysicsTools.h`. The `"ep"` branch at `analysis/src/PhysicsTools.cpp:291-295` computes:
  ```
  expectE = Ebeam*M_PROTON/(M_PROTON + Ebeam*(1-cos_t));
  return expectE - EnergyLoss(theta_deg, expectE);
  ```
  - `M_PROTON = 938.272` MeV (`:23`). Units are MeV and degrees.
  - `EnergyLoss` (`:316-336`) is about 2 MeV at θ=0, from the Al window, GEM foils and kapton.
  - **Not bound to Python.** `grep PhysicsTools python/` returns nothing.
- An existing use of the same cut idea is `analysis/tools/physics_calib.cpp:571-576`. It computes `theta = atan(r/hycal_z)` on module x,y, then the expected e-p peak, then a Gaussian fit.
- Beam energy and HyCal z come from `database/runinfo/general.json`, chained by `from_run`.
  - `from_run 25200`: `beam_energy 2239.51` MeV.
  - `defaults.hycal.position = [0,0,6267.0]` mm. The X17 period (`from_run 25540`) changes this to 7530.681 mm.
  - The target is at z = 0.
  - Python reads these from `p.geo.Ebeam` and `p.geo.hycal_z` (`python/bind_det.cpp:974-981`, `_common.py:236`).
- Theta in the full replay (`Replay.cpp:1074-1077`) builds the hit from `shower_depth(center_id, E)` (`det.shower_depth`, bound at `bind_det.cpp:476`; about 56 mm for PbWO4 at 2.2 GeV). It then applies `ApplyToLab` and `GetProjection` (`MatchingTools.cpp:41`).
- `benchmark_hycal_timing.py` currently stores only the seed-module x,y. The cluster `r.hit.x/y` is available (`bind_det.cpp:750-755`) but not stored.
- Energy resolution:
  - `database/reconstruction_config.json:26` sets `energy_resolution: [3.3, 0, 0]`, i.e. 3.3%/√E_GeV.
  - It is bound as `HyCalSystem.energy_resolution(E_MeV)` (`bind_det.cpp:646-653`; formula at `HyCalSystem.h:311-318`).
  - At E' ≈ 2.2 GeV, σ ≈ 49 MeV, so 2σ ≈ ±98 MeV.

### Inferences

- Run 25320 uses z = 6267 mm. That follows from the chained config. Print `p.geo.hycal_z` to confirm.
- The real width of E/E' after calibration may be wider than 3.3%/√E. `physics_calib` assumes 3%/√E in its fit-quality cut.

---

## 7. Q7: LMS state behind the absolute calibration

### Verified

- `database/calibration/calibration_factor_2p2_fadc.json` holds 1728 entries `{name, factor, base_energy, non_linear}`. It carries **no** LMS, ref-run or timestamp metadata.
- In `database/runinfo/general.json`, calibrations are paired with an LMS reference run only by appearing in the same configuration entry:
  - `calibration_factor_3p5_June7_nonLinear.json` pairs with `ref_run 24827` (from_run 0, 24185, 24670).
  - `x17_25560/25655/25825/26094/26214.json` each pair with the ref run of the same number.
  - `from_run 25200` (2.2 GeV) sets neither `calibration` nor `gain_factor`, so it **inherits ref_run 24827** and the 3.5 GeV calibration.
- `calibration_factor_2p2_fadc.json` does not appear in `general.json`. It is referenced only in `analysis/pyscripts/Commands.md:145,163`, `docs/analysis_notes/plan_waveform.md:16`, and the headers of `database/daq_cnf/adchycal_gain_2p2.cnf` and `adchycal_gain_v4_new_x17.cnf`.
- `prad_025320_gain_corr.root` is relative to `ref_run 24827`.

### Inference

- There is no reference (t0) LMS ratio tied to `calibration_factor_2p2_fadc.json`. The repo does not record which run it was derived from, or whether LMS correction was applied when it was derived.
- Applying the existing `gain_corr` factors (relative to 24827) on top of it moves the absolute scale by about 2% on average, and by up to about 12% for single modules (batch 0 range 0.876-1.157).
- `ref_gain/prad_025293_LMS.dat` is a reference taken within the 2.2 GeV period (25200 ≤ 25293 < 25540) but is not wired into `general.json`. It is a candidate t0 closer to the 2p2 calibration. Whether it is the calibration's actual source is unverified.
- The seed-module refinement would absorb a constant per-module offset. It would not absorb drift within the run.

---

## 8. Risks, open questions, gaps

### Open questions for the user

1. **Reference point.** Normalise to one of:
   - (a) the first 3-minute bin of the run. This gives pure relative drift and leaves the 2p2 calibration's absolute scale unchanged.
   - (b) the existing `ref_run` convention (24827 for 2.2 GeV). This is consistent with the C++ replay but not with how `calibration_factor_2p2_fadc.json` was made.
   - (c) the 2.2 GeV-period reference `prad_025293_LMS.dat`.
2. **Normalisation.** Use (a) the alpha-normalised reference ratio as the existing code does (`LMS_ref/alpha_ref`), or (b) a plain ratio to the reference PMT, as stated in the request. Option (b) assumes the reference PMT gain is stable. And which reference PMTs: LMS2+LMS3 as in the replay, or all three?
3. **Reuse or reimplement.** Reuse `replay_gainCorr -b ~1600` and read `gain_corr.root` via uproot, which covers PbWO4 only and has a step lookup by event number. Or recompute in Python from waveforms. If reimplementing, which integration window applies to LMS pulses: 160-220 ns as in the C++, or the benchmark's 130-200 ns?
4. **PbGlass.** The existing correction covers PbWO4 only, and PbGlass stays at 1.0. Is the benchmark restricted to PbWO4 (inner radius 200 mm)?

### Risks

- LMS selection by `nch > 1000` alone can admit large physics events. Selecting by `trigger_type` fails in 2.2 GeV runs (see §2).
- Failed fits or missing references appear as 0 in the ROOT file (W1141) and are mapped to 1.0 by the C++ loader. Downstream code should not mistake these for "no drift".

### Coverage gaps (not verified)

- How `trigger_type` and `trigger_bits` are distributed for LMS events in raw 25320 EVIO. No raw data was opened.
- The absolute LMS peak-time distribution per channel.
- Whether the phase 1 `WaveAnalyzer` settings match the Python pipeline exactly.
- Which run produced `calibration_factor_2p2_fadc.json`.
