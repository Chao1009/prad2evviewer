# Waveform Analysis Roadmap

**Author:** Jingyi Zhou (jyzhou)
**Last updated:** 2026-10-20
**Status:** Active development on branch `waveform_jz`

---

## 1. Completed Work

### Template extraction

- Per-material (PbWO4) template extraction pipeline (`fit_pulse_template.py`). Documented in `docs/technical_notes/waveform_analysis/fit_pulse_template.md`.
- χ² amplitude-scaling bias identified and corrected (`--model-err-floor 0.03`).
- Two shape populations (physics vs. background) identified and separated by fitted t₀ (`--t0-min 22.0`).
- Missing Veto enum in pybind11 fixed (`python/bind_det.cpp`, commit `e0ec70d`).
- Trigger-event-type filtering added (`--trigger-event-type`).
- Template stability verified across 3 runs (025308 Carbon, 025320 ep elastic, 026138 X17) and 4 trigger types (SSP_RawSum, SSP_Cluster, Pulser, LMS).
- Pile-up multiplicity distribution measured: 85% single-pulse, 0.7–5.1% two-peak pile-up (run-dependent).
- Six diagnostic plot scripts: `plot_chi2_vs_amp.py`, `plot_template_by_crate.py`, `plot_template_2d_map.py`, `plot_tau_vs_amp.py`, `plot_raw_pulses_by_amp.py`, `plot_chi2_vs_amp.py`.
- Trigger-scan batch runner: `run_trigger_scan.sh` (parallel execution).

### Pile-up benchmarking

- Synthetic pile-up generator (`pileup_generator.py`): 2-pulse injection at controlled (ΔT, ratio) grid or random ΔT. Plan documented in `docs/plans/plan-pileup-generator.md`.
- Deconvolver benchmark (`benchmark_deconv.py`): LM-based `WaveAnalyzer::Deconvolve` tested on 3000 synthetic events. WA separation floor: 24 ns (ratio ≥ 1), 48 ns (ratio = 0.5), >128 ns (ratio ≤ 0.3). Amplitude recovery ±1–5% at ΔT ≥ 48 ns.
- Matched-filter peak finder prototype (`benchmark_matched_filter.py`): cross-correlation with analytic two-τ kernel. MF dramatically outperforms WA at low amplitude ratios (86–96% detection at ratio ≤ 0.3, ΔT ≥ 96 ns vs. 6–18% for WA). WA outperforms MF at equal-amplitude close spacing (100% vs. 68–82% at ΔT = 24–32 ns, ratio = 1).

---

## 2. Near-Term Priorities

### 2a. Cluster-level timing gate (roadmap item 4)

**Goal:** Filter per-module hits by ΔT-to-seed-time at the cluster level. Modules whose pulse time disagrees with the cluster seed time by more than ΔT_max(E) are excluded from the cluster energy sum. This rejects accidental pile-up at the cluster level, complementing the per-channel deconvolution.

**What it needs:**
- Cluster reconstruction running via `prad2py` (exists in `prad2det`, accessible through bindings).
- Per-module pulse times from `WaveAnalyzer` (already available via `peak.time`).
- A configurable energy-dependent gate width ΔT_max(E).
- A script that runs cluster reconstruction with and without the gate on the same events and compares cluster energies.

**Deliverables:**
- `timing_gate_study.py` — apply the timing gate to reconstructed clusters and measure its effect on cluster energy resolution.
- Sweep ΔT_max to find the working point.
- ROC curve: signal efficiency vs. accidental rejection as a function of ΔT_max.

**Data:** Use the same Carbon/ep/X17 runs already analyzed (025308, 025320, 026138). The Population A vs. B t₀ difference (Section 2b of `fit_pulse_template.md`) directly motivates the timing gate.

---

### 2b. Calibration harness (roadmap item 5)

**Goal:** Run three energy estimators (firmware Mode-2 integral, WaveAnalyzer + timing gate, deconv + timing gate) on the same events through the same iteration/cuts/framework. Only the per-module energy estimator changes. Emit side-by-side σ(E)/E plots.

**What it needs:**
- Timing gate working (item 2a above).
- Calibration constants (gain per module) — at minimum, relative gain equalization.
- Energy resolution measurement on known-energy events (ep elastic or Möller peaks).
- A framework that configures and runs the three reconstruction paths.

**Deliverables:**
- `calibration_harness.py` — script implementing the three-method comparison.
- σ(E)/E vs. E for each method.
- Per-method reconstruction time per event.

---

### 2c. Luminosity / pile-up scan (roadmap item 6)

**Goal:** Measure timing-gate and deconvolver performance as a function of pile-up rate.

**Options:**
- Real data at varying beam current (if available).
- Synthetic pile-up overlaid on low-rate data using `pileup_generator.py` (recommended — controlled, reproducible).

**Deliverables:**
- Performance metrics (σ(E)/E, detection efficiency, fake rate) vs. pile-up rate.
- Calibration-constant stability vs. run number for each method.

---

## 3. Medium-Term Improvements

### 3a. Hybrid WA + MF peak finder

**Priority: HIGH.** WA and MF excel in complementary regimes. A hybrid approach would:

1. Run `WaveAnalyzer.findPeaks` as the primary pass (good at resolving close equal-amplitude peaks).
2. Run the matched filter as a second pass, scanning the residual (waveform minus fitted peaks from WA) for small pulses that WA missed.
3. Merge the two peak lists.

**Expected improvement:** Combined separation floor of ~20 ns at ratio ≥ 1 (WA) plus ~48 ns at ratio ≤ 0.3 (MF), vs. >128 ns for WA alone at low ratios.

**Implementation:** Python prototype first (using existing MF code), then C++ if results warrant (new method in `WaveAnalyzer`).

---

### 3b. MF-seeded deconvolution

**Priority: HIGH.** Currently the deconvolver is seeded by WA's peak positions. If WA misses a peak, the deconvolver cannot find it. Use MF-found peaks to seed the deconvolver instead:

1. Run MF to find all peaks (including small ones WA misses).
2. Use MF peak positions and amplitudes as initial guesses for the LM deconvolver.
3. LM refines the amplitudes and times.

This combines MF's detection sensitivity with the deconvolver's amplitude recovery accuracy. The deconvolver benchmark showed ±1–5% accuracy when both peaks are found — the bottleneck is finding them, which MF addresses.

**Implementation:** Modify `benchmark_deconv.py` to accept MF-found peaks as seeds instead of WA peaks. Test on the existing synthetic data. If improvement is confirmed, implement in C++ by adding an MF mode to `WaveAnalyzer`.

---

### 3c. Non-parametric templates

**Priority: MEDIUM.** The two-τ model has a ~3% leading-edge misfit (absorbed by `--model-err-floor 0.03`). Non-parametric templates (storing the full waveform shape as ~100 sample values per material, like Wassim's approach in `best_ref_shape_fit.C`) would:

- Eliminate the leading-edge misfit entirely.
- Give χ²/dof ≈ 1 without needing a model-error floor.
- Provide the natural input for an NNLS deconvolver (item 3d below).

**Implementation:** Extend `fit_pulse_template.py` with a `--template-mode discrete` option that averages sub-sample-aligned clean pulses instead of fitting a parametric model. The C++ `PulseTemplateStore` would need a code path to consume the discrete template.

---

### 3d. NNLS deconvolver (roadmap item 3 from jist.md)

**Priority: MEDIUM.** Replace the current per-event LM fit with non-negative least-squares on a fixed time grid with measured pedestal-autocorrelation noise covariance (CMS-MAHI style). Requires:

- Non-parametric templates (item 3c above).
- Pedestal autocorrelation measurement from real pedestal-only data.
- New C++ class `DeconvAnalyzer` in `prad2dec/`.
- pybind11 binding.
- Benchmark against synthetic pile-up data from `pileup_generator.py`.

---

### 3e. Multi-pulse pile-up generator

**Priority: LOW (for now).** Extend `pileup_generator.py` to inject 2–6 pulses per event, drawing pulse count from the measured npeaks distribution (Section 2c of `fit_pulse_template.md`) and ΔT from a Poisson process at the PRad-II accidental rate. Required for the "Poisson-random pile-up ROC vs. luminosity" plot from jist.md roadmap item 2, and for benchmarking at realistic pile-up multiplicities.

---

## 4. X17-Specific Work

- X17 runs (026138) show ~5% two-peak pile-up rate (8× higher than Carbon/ep elastic runs) and ~21% background contamination (vs. ~5% for other runs).
- X17 three-cluster topology means: three real clusters coincident in time, with accidental pile-up threats on any of them.
- The timing gate (item 2a) is particularly important for X17: it rejects both inter-cluster accidentals and the Population B background identified in `fit_pulse_template.md`.
- X17-specific deliverable: fake-three-cluster rate vs. real-three-cluster efficiency as a function of ΔT_max (ROC curve).

---

## 5. Reference Materials

| Document | Path |
|---|---|
| WaveAnalyzer algorithm documentation | `docs/technical_notes/waveform_analysis/wave_analysis.md` |
| Template extraction findings and benchmark results | `docs/technical_notes/waveform_analysis/fit_pulse_template.md` |
| Pile-up generator design plan (iteration 2, approved) | `docs/plans/plan-pileup-generator.md` |
| Wassim Hamdi's matched-filter code | `/Users/jingyi.zhou/research/PRad-II/Waveform_Wassim/` |

Wassim's code uses non-parametric templates + TSpectrum + TMinuit (NPS-heritage per-module reference extraction + multi-pulse fitting). Relevant for items 3a and 3c.

---

---

# Waveform Analysis — Plan & Progress Log

_Last updated: 2026-10-20_

## Overview

This document tracks the ongoing investigation of PRad-II HyCal waveform templates, per-module calibration quality, and energy spectrum artifacts across multiple runs. The work sits at the intersection of the physics_calib pipeline, waveform-template fitting, and benchmark_hycal_timing diagnostics.

---

## Completed Work

### 1. Pulse Template Comparison (25308 vs 25320)
**Status:** DONE  
**Script:** `analysis/pyscripts/dump_pbwo4_table.py`  
**Output:** `analysis/pyscripts/output/pulse_templates_h50.dat`

Per-module tau_r / tau_f comparison across 1151 PbWO4 modules between run 25308 (Carbon target) and run 25320 (H2 target). Produced a table with 11 columns including relative differences and metadata.

Key findings:
- τ_r relative difference: +1.18% (mean), ±4.49% stdev
- τ_f relative difference: -0.24% (mean), ±0.76% stdev
- Pulse shapes highly consistent between runs at the module level

### 2. Physics_calib Investigation on Run 25320
**Status:** DONE (documented as "insufficient-statistics divergence")

Attempted physics_calib iterative gain refinement on run 25320 with increasing statistics:
- 11 EVIO files (initial): σ_CB/E degraded from 3.75% → 5.89% after 7 iterations
- 88 EVIO files (full run): converged mathematically (ratio_stdev → 0.0018) but still gave σ_CB/E = 5.89%
- Root cause not just low statistics — see Section 4 below for the real cause

Validation: run 24386 (11 files, 3.49 GeV, H2) iter1 gave σ_CB/E = 2.0% cleanly, confirming the physics_calib tool works correctly when input conditions are right.

### 3. Timing-Gate Benchmark Study on Run 24386
**Status:** DONE  
**Script:** `analysis/pyscripts/sweep_pre_window_24386.py`  
**Outputs:**
- `analysis/pyscripts/output/timing_gate_benchmark_24386_summary.tsv`
- `analysis/pyscripts/output/timing_gate_benchmark_24386_curves_*.png` (3 plots)

Swept pre-window LO ∈ {50, 75, 100, 125, 130, 150, 175} ns (fixed HI=200 ns) on prad_024386.evio.00050 with iter1 calibration.

Key findings:
- σ_CB/E ≈ 2.0% across all tested LO values
- Pre-window width has minimal impact on energy resolution
- Runinfo baseline (LO=130, width=70 ns) is reasonable production default

### 4. Multi-Run EVIO Data Quality Scan
**Status:** DONE (older version) / IN PROGRESS (current scan with pre+post CB)  
**Script:** `analysis/pyscripts/benchmark_hycal_timing_scan.py`  
**Outputs:**
- `analysis/pyscripts/output/benchmark_hycal_timing_scan_summary.tsv`
- `analysis/pyscripts/output/benchmark_hycal_timing_scan_sigma_E_vs_file_pre.png`
- `analysis/pyscripts/output/benchmark_hycal_timing_scan_sigma_E_vs_file_post.png`
- `analysis/pyscripts/output/benchmark_hycal_timing_scan_all_calibrated_plots.pdf`
- `analysis/pyscripts/output/benchmark_hycal_timing_scan_findings.md`

Scans all EVIO files in `/Users/L00338853/work/PRad/data/evio/`, runs benchmark_hycal_timing.py (with --calibrate --cal-method cb-peak) on each, parses output TSVs, generates aggregate trend plots and combined PDF.

Current scope: 134 files across 11 runs (prad_024386.evio.00040 quarantined due to subprocess hang bug):
- 24386, 24503, 24512 (3.49 GeV, H2)
- 24590 (0.73 GeV, H2)
- 25260, 25299, 25308 (Carbon), 25320 (H2) — 2.24 GeV runs
- 25800, 26138, 26230 — X17 Ta runs at 2.24 GeV

Processing: full EVIO files (no event cap), parallel=8, pre-window per-run from runinfo.

### 5. Bimodality Analysis on Run 25320
**Status:** DONE  
**Script:** `analysis/pyscripts/scan_module_bimodality.py`  
**Outputs:**
- `analysis/pyscripts/output/module_bimodality_025320_split00010.tsv`
- `analysis/pyscripts/output/module_bimodality_025320_map.png`
- `analysis/pyscripts/output/module_bimodality_025320_spatial.png`
- `analysis/pyscripts/output/module_bimodality_025320_top20.pdf`
- `analysis/pyscripts/output/module_bimodality_025320_findings.md`

Classifies every PbWO4 module as unimodal/bimodal using scipy.signal.find_peaks with peak_2/peak_1 amplitude ratio ≥ 0.30.

Key finding: 33 of 97 modules (34%) in run 25320 show bimodal pre-calibration spectra. Spatial pattern: strongest near beam hole, decreasing with distance. Not seen in run 24386 (which is clean).

### 6. Per-Module Pre/Post Spectrum Diagnostic
**Status:** DONE  
**Script:** `analysis/pyscripts/plot_per_module_pre_post.py`  
**Output:** `analysis/pyscripts/output/per_module_025320_split00010_pre_post_spectra.pdf`

Per-module two-panel plots (pre and post calibration) with Crystal Ball fits. Used to visually confirm the bimodality mechanism — some modules have one peak at the elastic position and a second peak at a lower energy, both static over time.

### 7. Time-Split Diagnostic on W559/W593
**Status:** DONE  
**Script:** `analysis/pyscripts/plot_module_time_split.py`  
**Output:** `analysis/pyscripts/output/per_module_025320_split00010_time_split.pdf`

Split events into 3 equal-count time bins by event_num and overlaid pre-cal spectra per bin to test gain-drift hypothesis.

Result: All three time slices show identical double-peak structure — **not gain drift**. The bimodality is a static feature of the module's response within a given run.

---

## Current Open Investigations

### A. Why does run 25320 show bimodality but 24386 doesn't?
**Status:** OPEN — mechanism unclear

Observations so far:
- Run 25320 (2.24 GeV, H2): 34% of PbWO4 modules bimodal, pattern concentrated near beam hole
- Run 24386 (3.49 GeV, H2): clean, no bimodality
- Beam-hole edge leakage hypothesis: initially plausible but refuted by 24386 being clean (same geometry)
- Over-correction of pre-cal under-calibrated modules hypothesis: refuted by `--cal-method cb-peak` giving same result as median
- Pre-cal spectrum is clean at ~2600 MeV; post-cal shows a secondary bump → the bump is algorithm-induced, not from pile-up

Remaining candidate mechanisms:
1. Different physics process population at 2.24 GeV vs 3.49 GeV (needs investigating)
2. Different trigger composition between runs
3. Something unique to 25320 (not reproducible in other 2.24 GeV runs like 25260/25299)

The current full-detector scan (Section 4 above) with 25260, 25299, 25800, 26230 new H2/X17 2.24 GeV runs will tell us whether the issue is 25320-specific or a 2.24 GeV phenomenon.

### B. Signal-min cut adjustment for 24590
**Status:** DONE — code change applied, scan in progress

The 24590 (729 MeV beam) run was being processed with --signal-min 500 MeV, which truncated the elastic peak's low-energy radiative tail. The Gaussian σ was artificially inflated to ~8%.

Fix: `_SIGNAL_MIN_MAP` added to `benchmark_hycal_timing_scan.py`. 24590 now uses 200 MeV cut; all other runs keep 500 MeV. Scan re-run in progress to apply this to 24590.

### C. Pre-calibration Crystal Ball fit addition
**Status:** DONE — benchmark_hycal_timing.py now fits CB on both pre and post spectra, new TSV columns `mu_cb_pre_MeV`, `sigma_cb_pre_MeV`, `alpha_cb_pre`, `n_cb_pre`, `sigma_over_E_cb_pre` added. Scan script updated to parse new columns and plot both Gaussian + CB on the same chart (solid vs dashed lines). Full re-run needed to populate new columns (in progress with current scan).

### D. Known subprocess hang bug
**Status:** OPEN — workaround applied

File `prad_024386.evio.00040` reliably hangs benchmark_hycal_timing.py when processed. Observed twice independently. CPU time stays near zero while wall time grows indefinitely. Hangs ~2 minutes into processing (consistent).

Workaround: file moved to `/Users/L00338853/work/PRad/data/evio_quarantine/` so it's excluded from the scan.

Investigation needed: is the EVIO file corrupted, or does it trigger a C++ decoder bug? Not blocking current analysis.

### E. ThreadPoolExecutor orchestrator hangs in aggregation phase
**Status:** OPEN — workaround known

The scan orchestrator sometimes hangs after all worker subprocesses complete, before the aggregation phase. Pattern: all per-file TSVs present on disk, orchestrator process alive but 99.9% idle. Workaround: kill orchestrator, re-run with `--parallel 1 --skip-benchmark` to reuse the completed TSVs and aggregate single-threaded.

---

## Key Scripts and Their Purposes

| Script | Purpose |
|--------|---------|
| `benchmark_hycal_timing.py` | Core single-run benchmark — reads EVIO, clusters, fits elastic peak with Gaussian + CB, writes TSV with pre/post calibration metrics and histograms |
| `benchmark_hycal_timing_scan.py` | Driver that runs benchmark_hycal_timing.py on every EVIO file in a directory, aggregates results, detects secondary peaks, generates summary TSV + trend plots + combined PDF |
| `sweep_pre_window_24386.py` | Timing-gate pre-window sweep tool (1-D LO scan with fixed HI=200 ns) |
| `dump_pbwo4_table.py` | Pulse-template per-module comparison across runs |
| `scan_module_bimodality.py` | Classify PbWO4 modules as unimodal/bimodal, map spatially, flag worst offenders |
| `plot_per_module_pre_post.py` | Per-module pre/post spectrum diagnostic (one page per module) with CB fits |
| `plot_module_time_split.py` | Time-bin diagnostic for suspected gain-drift modules |

---

## Open Questions for Next Session

1. **25260/25299 bimodality:** Do the new 2.24 GeV H2 runs (25260, 25299) show the same bimodality pattern as 25320? If yes → 2.24 GeV configuration issue. If no → 25320-specific.
2. **X17 run comparison:** How do 25800, 26138, 26230 (Ta target, 2.24 GeV) compare to the H2 runs?
3. **Scan hang investigation:** Can we identify the mechanism of the orchestrator hang and fix it, or is a known-good `--parallel 1 --skip-benchmark` recovery pattern sufficient?
4. **24386/split_00040 hang:** Is it the file (corruption), the content (some pathological event), or something else? Low priority but needs eventual resolution.

---

## File Paths Reference

- EVIO data: `/Users/L00338853/work/PRad/data/evio/` (quarantine: `/Users/L00338853/work/PRad/data/evio_quarantine/`)
- Raw ROOT (older): `/Users/L00338853/work/PRad/data/raw_root/`
- Physics_calib outputs: `/Users/L00338853/work/PRad/data/physics_calib_out/Physics_calib/`
- Analysis scripts: `/Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts/`
- Analysis outputs: `/Users/L00338853/work/PRad/prad2evviewer/analysis/pyscripts/output/`
- Calibration JSONs: `/Users/L00338853/work/PRad/prad2evviewer/database/calibration/`
- Runinfo: `/Users/L00338853/work/PRad/prad2evviewer/database/runinfo/general.json`

---

## Fit-model note (2026-10-21, user observation from the full scan `output/scan_20261021_full/`)
- **The Crystal Ball (CB) fit describes the elastic peak well** for 3.49 GeV e-p (H2: 24386, 24503, 24512), 2.24 GeV carbon (25308) and the 2.24 GeV X17 Ta runs (25800, 26138, 26230).
- **The Gaussian fit works better** for the other runs: 2.24 GeV e-p H2 (25228, 25257, 25260, 25299, 25320) and 0.73 GeV e-p H2 (24590, 24603).
- Quote σ_E/E from the preferred fit model for each run type.
- **Run 24578 is excluded** from the scan (2026-10-21): its trigger was not the standard one, leaving only ~2000 physics-trigger events per split and no per-module refinement. It is excluded via `_EXCLUDED_RUNS` in `benchmark_hycal_timing_scan.py`.
