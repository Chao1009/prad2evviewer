# Waveform Analysis Roadmap

**Author:** Jingyi Zhou (jyzhou)
**Last updated:** 2026-09-21
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
