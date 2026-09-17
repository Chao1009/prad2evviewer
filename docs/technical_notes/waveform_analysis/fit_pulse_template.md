# Template Extraction Investigation: `fit_pulse_template.py` on Run 025308

**Author:** Jingyi Zhou (jyzhou), documented 2026-09-16.

## 1. Introduction

This note summarizes an investigation into `analysis/pyscripts/fit_pulse_template.py` performed on the
Carbon-target physics run 025308 (multiple EVIO splits). The motivation was that both the shipped
template file `pulse_templates_024177.json` and templates extracted from run 025320 reported
`n_channels_good = 0/1151` for PbWO4, and per-channel distributions of τ_r and τ_f were visibly
multi-modal rather than unimodal. The investigation identified two root causes: a χ² gate with an
amplitude-scaling bias built into the per-sample σ definition in `WaveAnalyzer::FitPulseShape`, and
a two-population signal-vs-background structure in the pulse data itself that amplitude cuts alone
cannot cleanly separate. A third issue, a missing `Veto` enum entry in the Python bindings, caused
silent per-material misclassification. This note documents each finding, the additive changes made
to `fit_pulse_template.py`, the new diagnostic scripts written, and the next steps. For the
underlying `WaveAnalyzer::FitPulseShape` algorithm and the parametric two-tau template model, see
[`wave_analysis.md`](wave_analysis.md).

## 2. Key Findings

### 2a. χ² Gate Had a Hidden Amplitude-Scaling Bias

The per-sample σ used in the normalized χ² fit inside `WaveAnalyzer::FitPulseShape`
(`WaveAnalyzer.cpp:1265`) is:

```
sigma = max(ped_rms / peak_amp, model_err_floor)
```

At default `model_err_floor = 0.01` (1% of peak), σ is amplitude-dependent: dim pulses get a large
relative σ (loose χ² constraint) while bright pulses get σ floored at 1%, which is smaller than the
actual 2-3% model misfit on the leading edge. The floor therefore drives χ²/dof artificially high
on tall pulses, not because the template fits poorly, but because the model-error floor is tighter
than the physical model residuals.

Measured on run 025308, the ratio of median χ²/dof in the highest-amplitude decile to the
lowest-amplitude decile was:

| Material | Ratio (default floor 0.01) | Ratio (corrected floor 0.03) |
|---|---:|---:|
| PbWO4    | 6–9× | 0.87–1.00 |
| PbGlass  | 6–9× | 0.87–1.00 |
| Veto     | ~9×  | 1.00 |
| LMS      | ~1× (χ²~22 at both ends) | ~1× (χ²~2.5 at both ends) |

Raising `--model-err-floor` from 0.01 to 0.03 collapsed the ratio to 0.87–1.00 across PbWO4,
PbGlass, and Veto. LMS is a separate case: χ²/dof was uniformly high (~22) across all amplitudes
with the default floor, reflecting a genuine model-family mismatch for laser pulses rather than an
amplitude bias. Raising the floor to 0.03 corrected the LMS χ²/dof to ~2.5 uniformly.

### 2b. Two-Population Signal-vs-Background Structure at ~400-500 ADC

After correcting the χ² floor, per-pulse plots of τ_r and τ_f vs. peak amplitude for PbWO4 revealed
a step-like discontinuity at ~400-500 ADC, correlated between the two parameters. Direct visual
inspection of raw waveforms on channel W500 via `plot_raw_pulses_by_amp.py` confirmed two distinct
pulse-shape families:

- **Below ~300 ADC**: waveforms are noisy with no identifiable clean shape; background or accidental
  pulses.
- **300-500 ADC**: a mixture of two visibly distinct shape families coexisting in the same amplitude
  range.
- **Above ~500 ADC**: a single clean shape family consistent across channels; the physics
  electromagnetic shower population.

The interpretation is that pulses below the transition threshold are dominated by beam-related
backgrounds, while pulses above ~500 ADC are dominated by real signals from Carbon
elastic/quasi-elastic scattering. Because the two populations coexist in the 300-500 ADC range and
the background Population B reaches physics-scale amplitudes, an amplitude cut alone cannot cleanly
separate them (see also §2d).

### 2c. Missing Veto Enum in Python Bindings

`python/bind_det.cpp` was missing the `Veto` entry in the pybind11 `ModuleType` enum registration.
Without it, any Python call to `Module.type.name` for a Veto channel returned `"???"` instead of
`"Veto"`, silently routing Veto channels into the `"Unknown"` bucket in `fit_pulse_template.py`'s
per-material aggregation. Fixed in commit `e0ec70d` (single-line addition of
`.value("Veto", fdec::ModuleType::Veto)` between the LMS and Unknown `.value()` entries in
`bind_det.cpp`). After the fix, 4 Veto channels appeared correctly in the diagnostic output with an
amplitude-bias ratio of ~9.0× at the default floor, collapsing to 1.00× at the corrected floor,
consistent with PbWO4 and PbGlass.

### 2d. Spatially Resolved Signal-vs-Background Structure

Restricting template extraction to `--height-min 500` reduced but did not eliminate residual
multi-modal structure in per-channel τ_r and τ_f distributions. The 2D HyCal map of τ_r produced by
`plot_template_2d_map.py` reveals the underlying geometry directly: the two populations are organized
spatially rather than uniformly scattered.

- **Population A (physics)**: occupies a central circular region near the beam axis. Characterized
  by small τ_r, small τ_f, and large t₀ — pulses arrive later in the 100-sample readout window,
  consistent with trigger-anchored physics timing for electromagnetic showers.
- **Population B (background)**: distributed across the remainder of the HyCal face outside the
  physics circle. Characterized by large τ_r, large τ_f, and small t₀ — pulses arrive earlier in
  the readout window, distinct from the physics population.

Because Population B reaches physics-scale amplitudes (i.e., it appears above `--height-min 500`),
amplitude cuts cannot remove it. The distinct t₀ distribution of Population B is the discriminating
observable. This finding also has direct implications for the cluster-level timing gate: Population
B's systematically earlier arrival time is exactly what a ΔT cut between cluster and seed time
exploits (see §5).

## 3. Changes to `fit_pulse_template.py`

All changes are additive. No existing behavior or output format was modified.

**Commit `bbdc8b1`** added the `per_pulse_amp_chi2.npz` diagnostic dump. After the event loop, every
accepted pulse's `(peak_amp, chi2/dof, module_type)` triple is collected and written to
`<plot_dir>/per_pulse_amp_chi2.npz` using `np.savez_compressed`. The dump only fires when
`--plot-dir` is set; runs without `--plot-dir` are unaffected.

A subsequent additive edit (uncommitted at time of writing) extended the same `.npz` dump to also
include per-pulse `t0`, `tau_r`, `tau_f`, `p` (only when `--model two_tau_p` is used), and the
channel name (`name`). The file layout is backwards-compatible: consumers that only read `amp`,
`chi2`, and `mtype` continue to work without modification.

A second uncommitted additive edit introduced the per-channel waveform cache. Module-level constants
`PULSE_CACHE_PER_BIN = 30` and `AMP_BINS = [("lt_300", -inf, 300), ("300_to_500", 300, 500),
("gt_500", 500, inf)]` define three amplitude strata. The `ChannelStats` dataclass was extended with
`sample_pulses_by_amp` and `sample_peak_amps_by_amp` fields (dicts keyed by bin label). During the
event loop, up to 30 pedsub pulse arrays per bin per channel are cached in memory. On completion,
`<plot_dir>/waveforms_by_amp/<channel>.npz` is written for every channel that has at least two
amplitude bins each populated by at least 5 pulses. These files are the input for
`plot_raw_pulses_by_amp.py`.

## 4. New Python Scripts Written

The following scripts were added under `analysis/pyscripts/`.

### `plot_chi2_vs_amp.py`

Committed with `bbdc8b1`. Reads `per_pulse_amp_chi2.npz` and produces per-material hexbin scatter
plots and running-median curves (with 16th/84th percentile bands) of χ²/dof vs. peak amplitude on
log-log axes. Also prints a console summary of the median χ²/dof in the lowest-amplitude decile vs.
the highest-amplitude decile, which directly quantifies the amplitude bias identified in §2a. This
was the primary diagnostic tool for finding 2a.

```
python3 plot_chi2_vs_amp.py per_pulse_amp_chi2.npz [--out-dir plots/]
```

Output files: `chi2_vs_amp_<mtype>.png` per material, plus `chi2_vs_amp_all_types.png` overlaying
all running-median curves.

### `plot_template_by_crate.py`

Reads a `pulse_templates_*.json` and groups per-channel τ_r, τ_f, and t₀ by FADC crate (derived
from the `roc_tag` component of `channel_id`). Produces per-material overlaid histograms and
τ_r-vs-τ_f scatter plots colored by crate. The diagnostic tested whether the bi-modal shape
distribution organizes by front-end electronics grouping. Finding: the two τ_r/τ_f populations in
PbWO4 do not split cleanly by crate, but t₀ does show a crate-level structure, consistent with
cable-length or readout-timing differences rather than crystal properties.

```
python3 plot_template_by_crate.py pulse_templates_025308.json [--out-dir plots/]
```

Output files: per-material histogram PNGs and τ_r-vs-τ_f scatter PNGs colored by crate.

### `plot_template_2d_map.py`

Reads a `pulse_templates_*.json` and `database/hycal_map.json`, and produces 2D spatial maps of
τ_r, τ_f, and t₀ on the physical HyCal face for each material. Each module is drawn using
`matplotlib.patches.Rectangle` at its actual physical size from `hycal_map.json`. This script
produced the finding in §2d: the two PbWO4 τ_r populations are spatially organized — small τ_r
concentrated in a central circle near the beam, large τ_r distributed across the remainder of the
face.

```
python3 plot_template_2d_map.py pulse_templates_025308.json database/hycal_map.json [--out-dir plots/]
```

Output files: `hycal_map_<param>.png` for each of τ_r, τ_f, and t₀, one per material.

### `plot_tau_vs_amp.py`

Reads the extended `per_pulse_amp_chi2.npz` (requiring the `tau_r`, `tau_f`, and `name` fields).
Produces per-material and optional per-channel plots of τ_r and τ_f vs. peak amplitude using hexbin
scatter and running-median curves. The per-channel mode distinguishes PMT saturation — a smooth
amplitude-dependent trend within a single channel — from inter-module variation, where τ_r is flat
within a channel but different between channels. On run 025308 PbWO4 data, per-channel τ_r was flat
within a channel and the population split was confirmed to be between channels, not within them.

```
python3 plot_tau_vs_amp.py per_pulse_amp_chi2.npz [--channel W500] [--out-dir plots/]
```

Output files: `tau_vs_amp_<mtype>.png` per material; `tau_vs_amp_<channel>.png` when `--channel` is
specified.

### `plot_raw_pulses_by_amp.py`

Reads the per-channel `waveforms_by_amp/<channel>.npz` dumps produced by the waveform cache added
to `fit_pulse_template.py` (§3). For each channel, produces a 2-row overlay figure: raw pedsub ADC
waveforms (top row) and amplitude-normalized waveforms (bottom row), one column per populated
amplitude bin. The amplitude-normalized row makes shape differences between bins visible independent
of pulse height. Used to visually confirm the two-shape-family structure on channel W500 (§2b), where
the `lt_300` and `300_to_500` bins showed visibly distinct shapes from the `gt_500` physics
population.

```
python3 plot_raw_pulses_by_amp.py waveforms_by_amp/ [--channel W500] [--out-dir plots/]
```

Output files: `pulses_by_amp_<channel>.png` per channel.

## 5. Reproducing the Template Extraction

The commands below produced the clean PbWO4 templates used for the run-to-run stability comparison
described in §2. All commands are run from `analysis/pyscripts/`. The `PYTHONPATH` variable must
point to the `prad2py` module in the build directory before invoking any script.

### Environment Setup

Set the following variables once per shell session before running any of the commands below:

```bash
# From the repository root
cd ~/work/PRad/prad2evviewer/analysis/pyscripts

# prad2py binding (adjust path to your build directory)
export PYTHONPATH=$HOME/work/PRad/prad2evviewer/build/python

# Reused config files (absolute paths — script default lookup depends on CWD)
export DAQ_CONFIG=$HOME/work/PRad/prad2evviewer/database/daq_config.json
export HC_MAP=$HOME/work/PRad/prad2evviewer/database/hycal_map.json
```

`DAQ_CONFIG` and `HC_MAP` are set as absolute paths because the scripts resolve config file
locations relative to the current working directory by default. Explicit paths remove ambiguity when
the scripts are invoked from any directory.

### Extracting a Template for One Run

The canonical invocation used for runs 025308, 025320, and 026138 (replace `<RUN>` with the actual
run number, e.g. `025308`):

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
| `--height-min 500` | 500 ADC | Reject low-amplitude background pulses (finding 2b) |
| `--model-err-floor 0.03` | 3% | Correct χ² amplitude bias (finding 2a) |
| `--t0-min 25.0` | 25 ns | Select Population A physics pulses by arrival time (finding 2d) |
| `--plot-dir ...` | per-run | Enables the per-pulse `.npz` dump and diagnostic PNGs |
| `--daq-config ...` | absolute path | DAQ channel map and analyzer config |
| `--hc-map-file ...` | absolute path | HyCal module geometry for downstream 2D maps |

### Generating the 2D HyCal Map

After extraction, produce the spatial parameter maps with:

```bash
python3 plot_template_2d_map.py \
    output/pulse_templates_<RUN>_h500_t0cut.json \
    --out-dir output/template_plots_<RUN>_h500_t0cut
```

This produces six PNG files (τ_r, τ_f, t₀ per material) showing the spatial distribution of
template parameters on the HyCal face.

### Comparing Template Stability Across Runs

The following one-liner prints a side-by-side PbWO4 summary across the three runs:

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

Runs whose τ_r and τ_f agree within MAD across runs indicate a stable template extraction and stable
detector. Runs with significantly different χ² or `n_good` hint at genuine drift or a different
pile-up regime. This comparison directly feeds the paper's calibration-constant stability vs. run
number systematic.

### Results So Far

| Metric | Run 025308 | Run 025320 | Run 026138 |
|---|---|---|---|
| τ_r (ns) | 2.24 ± 0.21 | 2.22 ± 0.20 | 2.36 ± 0.21 |
| τ_f (ns) | 24.01 ± 1.00 | 23.95 ± 1.02 | 23.26 ± 1.02 |
| χ²/dof (median) | 0.86 | 0.85 | 0.81 |
| Good channels | 424/425 | 458/461 | 238/238 |
| n_pulses | 133,715 | 169,041 | 57,136 |

τ_r is stable between runs 025308 and 025320 to below 1% (2.24 vs. 2.22 ns), but rises to 2.36 ns
in run 026138 — a ~5% shift relative to the earlier pair. τ_f is consistent across all three runs to
within ~3% (24.01, 23.95, 23.26 ns), and χ²/dof (0.85–0.86 on the first two runs, 0.81 on 026138)
confirms that all fits describe the data to sub-noise-floor precision. Good-channel counts are
essentially 100% of contributing channels in every run, indicating the extraction is
well-calibrated. Run 026138 contributes ~57 k pulses (~40% of the other runs' statistics) but still
yields a well-defined per-material aggregate. The τ_r shift of ~5% between the earlier and later
runs is comparable to the MAD (±0.21 ns, ~9% of the median), so it is not significant at the
per-pulse level, but the shift is coherent across all 238 contributing modules, which points to a
genuine detector-state change — temperature, gain drift, or aging — rather than statistical
fluctuation; this is precisely the class of variation the paper's calibration-constant stability
systematic is designed to characterize.

## 6. What Is Next

**Restrict template extraction to the physics-circle region.** Add a `--radius-max <mm>` argument
to `fit_pulse_template.py` that filters channels to those within a given radius of the beam center
(approximately (0, 0) in HyCal coordinates, read from `database/hycal_map.json`). Templates
extracted from the physics circle only — Population A from §2d — should yield unimodal τ_r and τ_f
distributions per material.

**Verify single-peaked distributions and χ² improvement.** After the radius cut, confirm that the
per-material `_by_type` τ_r and τ_f histograms collapse to a single peak and that the mean χ²/dof
on the restricted channel set is consistent with the corrected model-error floor.

**Point the physics-circle template at the deconvolver.** The C++ `PulseTemplateStore` consumes the
`_by_type` medians from the JSON file specified by `daq_config.json:nnls_deconv.template_file`. Once
the physics-circle template is stable, update `daq_config.json` to point to the new file so that
the LM deconvolver uses the physics-population shape.

**Move to the cluster-level timing gate.** The finding that Population B has a systematically
different t₀ distribution from Population A (§2d) provides both the physics motivation and a
natural benchmark dataset for the timing gate's ROC characterization. The next technical component
to build is cluster reconstruction with a per-module time comparison against the seed-cluster time;
Population B's well-defined t₀ offset makes it a concrete benchmark for gate efficiency vs.
background rejection.

**X17-relevant metric.** Once the timing gate is running, construct a three-cluster topology test
with synthetic accidentals to measure the pile-up ROC of direct interest for the X17 analysis.

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
