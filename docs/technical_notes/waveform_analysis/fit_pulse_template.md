# PbWO4 Pulse Template Extraction: `fit_pulse_template.py`

**Author:** Jingyi Zhou (jyzhou), documented 2026-09-16.

## 1. Purpose

This note documents the clean per-material template extraction for PbWO4, ready for use by the
pile-up deconvolver. For the underlying `WaveAnalyzer::FitPulseShape` algorithm and the parametric
two-tau template model, see [`wave_analysis.md`](wave_analysis.md).

## 2. Results

### 2a. Run Summary

| Run    | Physics                     | Beam Energy  | Target            |
|--------|-----------------------------|--------------|-------------------|
| 025308 | Elastic e-Carbon scattering | 2239.51 MeV  | Carbon            |
| 025320 | ep elastic scattering       | 2239.51 MeV  | Hydrogen          |
| 026138 | X17 search                  | 2239.51 MeV  | X17 physics run   |

All three runs share the same beam energy of 2239.51 MeV.

### 2b. PbWO4 Template Stability Across Runs

| Metric               | Run 025308 | Run 025320 | Run 026138 |
|----------------------|------------|------------|------------|
| τ_r (ns)             | 2.24 ± 0.21 | 2.22 ± 0.20 | 2.36 ± 0.21 |
| τ_f (ns)             | 24.01 ± 1.00 | 23.95 ± 1.02 | 23.26 ± 1.02 |
| χ²/dof (median)      | 0.86       | 0.85       | 0.81       |
| Good channels        | 424/425    | 458/461    | 238/238    |
| n_pulses used        | 133,715    | 169,041    | 57,136     |
| Selection efficiency | ~99.8%     | ~99.3%     | 100%       |

Selection efficiency is `n_pulses_used / n_pulses_attempted` — the fraction of pulses passing the
clean-pulse gate (single peak, no pile-up flag, above `--height-min`, in-window, not overflowed)
whose LM fit also converged and landed inside the `--t0-min` window. Actual per-channel values may
be lower for individual channels — this is the aggregate. Efficiency of the pre-fit clean-pulse gate
itself relative to all `WaveAnalyzer`-found peaks is not computed in the current pipeline.

τ_r between runs 025308 and 025320 agrees to below 1% (2.24 vs. 2.22 ns), but rises to 2.36 ns on
run 026138 — a ~5% shift. This drift is smaller than the per-pulse MAD (~0.21 ns, ~9% of the
median), so it is not statistically significant per-pulse, but the shift is coherent across all 238
contributing modules, pointing to a genuine detector-state change rather than statistical
fluctuation. τ_f is stable to within ~3% across all three runs (24.01, 23.95, 23.26 ns). χ²/dof
stays in the 0.81–0.86 range on all runs — all fits describe the data to sub-noise-floor precision.
The good-channel fraction is essentially 100% on every run. Run 026138 contributes ~40% of the
pulse statistics of the other runs but still yields a well-defined per-material aggregate.

**Per-run summary plots.**

Each per-run summary combines the per-channel-median distributions of the fit parameters (τ_r, τ_f,
t₀, peak amplitude, χ²/dof) split by module type, along with a τ_r-vs-τ_f scatter. All three runs
used the same `--height-min 500 --model-err-floor 0.03 --t0-min 25.0` cuts.

**Run 025308 (Carbon target).**

![Template extraction summary for run 025308](plots/template_summary_025308.png)

**Run 025320 (ep elastic).**

![Template extraction summary for run 025320](plots/template_summary_025320.png)

**Run 026138 (X17 physics).**

![Template extraction summary for run 026138](plots/template_summary_026138.png)

### 2c. JSON Ready for Deconvolution

`output/pulse_templates_025308_h500_t0cut.json` (or the equivalent per-run file) is the template
ready for use by the C++ pile-up deconvolver. To enable it in production, set
`database/daq_config.json`'s `fadc250_waveform.analyzer.nnls_deconv.template_file` to the deployed
template path and set `enabled` to `true`.

## 3. Reproducing the Extraction

All commands run from `analysis/pyscripts/`.

### Environment Setup

```bash
# prad2py binding (adjust path to your build directory)
export PYTHONPATH=$HOME/work/PRad/prad2evviewer/build/python

# Reused config files — explicit absolute paths remove CWD ambiguity
export DAQ_CONFIG=$HOME/work/PRad/prad2evviewer/database/daq_config.json
export HC_MAP=$HOME/work/PRad/prad2evviewer/database/hycal_map.json
```

### Extracting a Template for One Run

Replace `<RUN>` with the run number (e.g. `025308`):

```bash
python3 fit_pulse_template.py \
    ~/work/PRad/data/evio/prad_<RUN>.evio.* \
    -o output/pulse_templates_<RUN>_h500_t0cut.json \
    --max-events 0 \
    --height-min 500 \
    --model-err-floor 0.03 \
    --t0-min 25.0 \
    --plot-dir output/template_plots_<RUN>_h500_t0cut \
    --daq-config $DAQ_CONFIG \
    --hc-map-file $HC_MAP
```

| Flag | Value | Purpose |
|---|---|---|
| `--max-events 0` | all events | Full run statistics |
| `--height-min 500` | 500 ADC | Reject low-amplitude background pulses (§4b) |
| `--model-err-floor 0.03` | 3% | Correct χ² amplitude bias (§4a) |
| `--t0-min 25.0` | 25 ns | Select Population A physics pulses by arrival time (§4b) |
| `--plot-dir ...` | per-run | Enables the per-pulse `.npz` dump and diagnostic PNGs |
| `--daq-config ...` | absolute path | DAQ channel map and analyzer config |
| `--hc-map-file ...` | absolute path | HyCal module geometry for downstream 2D maps |

### Generating the 2D HyCal Map

```bash
python3 plot_template_2d_map.py \
    output/pulse_templates_<RUN>_h500_t0cut.json \
    --out-dir output/template_plots_<RUN>_h500_t0cut
```

Produces six PNG files (τ_r, τ_f, t₀ per material) showing the spatial distribution of template
parameters on the HyCal face.

### Three-Run Comparison One-Liner

```bash
python3 -c "
import json

runs = {
    '025308': 'output/pulse_templates_025308_h500_t0cut.json',
    '025320': 'output/pulse_templates_025320_h500_t0cut.json',
    '026138': 'output/pulse_templates_026138_h500_t0cut.json',
}
data = {r: json.load(open(p))['_by_type']['PbWO4'] for r, p in runs.items()}

print(f'{\"metric\":15s}' + ''.join(f'  {r:>12s}' for r in runs))
print('-' * (15 + 14 * len(runs)))

def row(label, get):
    print(f'{label:15s}' + ''.join(f'  {get(d):>12s}' for d in data.values()))

row('τ_r (ns)',    lambda d: f'{d[\"tau_r_ns\"][\"median\"]:.2f}±{d[\"tau_r_ns\"][\"mad\"]:.2f}')
row('τ_f (ns)',    lambda d: f'{d[\"tau_f_ns\"][\"median\"]:.2f}±{d[\"tau_f_ns\"][\"mad\"]:.2f}')
row('χ²/dof',      lambda d: f'{d[\"chi2_per_dof\"][\"median\"]:.2f}')
row('n_good',      lambda d: f'{d[\"n_channels_good\"]}/{d[\"n_channels\"]}')
row('n_pulses',    lambda d: f'{d[\"n_pulses_total\"]}')
"
```

## 4. Findings and Discussion

### 4a. χ² Gate Had a Hidden Amplitude-Scaling Bias

The per-sample σ in `WaveAnalyzer::FitPulseShape` (`WaveAnalyzer.cpp:1265`) is:

```
sigma = max(ped_rms / peak_amp, model_err_floor)
```

At the default `model_err_floor = 0.01` (1%), σ is amplitude-dependent: dim pulses get a large
relative σ while bright pulses are floored at 1%, which is tighter than the actual 2–3% model
misfit on the leading edge. This drives χ²/dof artificially high for tall pulses.

Measured on run 025308, the ratio of median χ²/dof in the highest-amplitude decile to the
lowest-amplitude decile was:

| Material | Ratio (floor 0.01) | Ratio (floor 0.03) |
|---|---:|---:|
| PbWO4    | 6–9× | 0.87–1.00 |
| PbGlass  | 6–9× | 0.87–1.00 |
| Veto     | ~9×  | 1.00 |
| LMS      | ~1× (χ²~22 uniformly) | ~1× (χ²~2.5 uniformly) |

Raising `--model-err-floor` to 0.03 collapsed the ratio to 0.87–1.00 for PbWO4, PbGlass, and Veto.
LMS had uniformly high χ²/dof (~22) at the default floor — a genuine model-family mismatch for
laser pulses, not an amplitude bias; raising the floor to 0.03 corrected it to ~2.5 uniformly.
Fix: use `--model-err-floor 0.03`.

### 4b. Two Shape Populations Distinguished by t₀

At `--height-min 500`, PbWO4 τ_r and τ_f still showed multi-modal distributions. The 2D HyCal map
of τ_r from `plot_template_2d_map.py` revealed the cause: the populations are spatially organized,
not randomly scattered.

- **Population A (physics)**: central circular region near the beam axis. Small τ_r, small τ_f,
  large t₀ — pulses arrive later in the readout window, consistent with trigger-anchored physics
  timing for electromagnetic showers.
- **Population B (background)**: distributed across the HyCal face outside the physics circle.
  Large τ_r, large τ_f, small t₀ — earlier arrival time, distinct pulse shape.

Because Population B reaches physics-scale amplitudes, amplitude cuts alone cannot remove it. The
distinct t₀ distribution is the discriminating observable. A per-pulse `--t0-min 25` cut cleanly
selects Population A. The different t₀ between populations also directly motivates the
cluster-level timing gate: Population B's systematically earlier arrival is exactly what a ΔT cut
between cluster and seed time exploits.

### 4c. Veto Enum Missing from Python Bindings

`python/bind_det.cpp` was missing the `Veto` entry in the pybind11 `ModuleType` enum. Python calls
to `Module.type.name` for Veto channels returned `"???"`, silently routing Veto channels into the
`"Unknown"` bucket in per-material aggregation. Fixed in commit `e0ec70d` (single-line addition of
`.value("Veto", fdec::ModuleType::Veto)` in `bind_det.cpp`). After the fix, 4 Veto channels
appeared correctly with an amplitude-bias ratio of ~9× at the default floor, collapsing to 1.00×
at the corrected floor, consistent with PbWO4 and PbGlass.

## 5. Diagnostic Scripts

All under `analysis/pyscripts/`. See individual `--help` output for usage.

| Script | Purpose |
|---|---|
| `fit_pulse_template.py` | Template extraction with per-pulse `.npz` and waveform-by-amp-bin dumps |
| `plot_chi2_vs_amp.py` | χ²/dof vs. peak amplitude, per-material scatter + running median. Diagnostic for finding §4a (committed with `bbdc8b1`) |
| `plot_template_by_crate.py` | Per-material params grouped by FADC crate. Diagnostic for testing electronics origin |
| `plot_template_2d_map.py` | 2D spatial maps of τ_r, τ_f, t₀ on HyCal face. Diagnostic for finding §4b |
| `plot_tau_vs_amp.py` | Per-pulse τ_r and τ_f vs. peak amplitude, per-material and per-channel |
| `plot_raw_pulses_by_amp.py` | Per-channel raw waveform overlays split by amplitude bin |

## 6. What Is Next

**Pile-up deconvolution.** Feed the clean per-material template into the existing C++
`WaveAnalyzer::Deconvolve`. The `output/pulse_templates_<RUN>_h500_t0cut.json` files are ready;
the only step required is pointing `daq_config.json` at the deployed template and enabling the
deconvolver (see §2c).

**Synthetic pile-up generation and ROC characterization.** New scripts are needed to inject
synthetic pile-up at known separations, run the deconvolver, and measure detection efficiency vs.
false-positive rate. This characterizes the deconvolver's operating regime before unblinding.

**Cluster-level timing gate.** The different t₀ distributions of Populations A and B (§4b) provide
both the physics motivation and a concrete benchmark dataset for gate efficiency vs. background
rejection. The next component to build is cluster reconstruction with a per-module time comparison
against the seed-cluster time.

**Multi-run stability curve.** Extend the three-run comparison to additional runs as data become
available. The coherent τ_r drift seen on run 026138 should be tracked as a calibration-constant
stability systematic.

## 7. See Also

- [`docs/technical_notes/waveform_analysis/wave_analysis.md`](wave_analysis.md) — parent note
  documenting `WaveAnalyzer` and `Fadc250FwAnalyzer`, including the parametric two-tau template
  model used by `FitPulseShape`.
- [`analysis/pyscripts/fit_pulse_template.py`](../../../analysis/pyscripts/fit_pulse_template.py)
  — the template extraction script.
- [`python/bind_det.cpp`](../../../python/bind_det.cpp) — pybind11 bindings for detector types;
  Veto fix landed in commit `e0ec70d`.
- [`prad2dec/src/WaveAnalyzer.cpp:1241`](../../../prad2dec/src/WaveAnalyzer.cpp) — C++
  implementation of `FitPulseShape` and `FitPulseShapeTwoTauP`.
- [`database/hycal_map.json`](../../../database/hycal_map.json) — module geometry and DAQ mapping
  used by `plot_template_2d_map.py`.
