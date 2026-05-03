# Software Waveform Analysis in `prad2dec`

**Author:** Chao Peng (Argonne National Laboratory)

`prad2dec` provides two offline analyzers for raw FADC250 samples
(`uint16_t[nsamples]`, 4 ns/sample at 250 MHz):

| Analyzer | Class | Purpose |
|---|---|---|
| Waveform | `fdec::WaveAnalyzer` | Robust local-maxima peak finding for HyCal energy and time observables. Tolerates noisy pedestals; supports multiple peaks per channel. |
| Firmware emulator | `fdec::Fadc250FwAnalyzer` | Bit-faithful emulation of the FADC250 firmware Mode 1/2/3 (Hall-D V3 with NSAT/NPED/MAXPED extensions) for offline–firmware cross-checks. |

Both analyzers are stack-allocated and execute side-by-side under
`prad2ana_replay_rawdata -p`. The remainder of this note describes the
algorithms on a representative waveform with parameter values from
[`database/daq_config.json`](../../../database/daq_config.json).

## Example waveform

100 samples (400 ns); a bright pulse on a quiet 146 ADC baseline,
followed by a slow scintillation tail.

![overview](plots/fig1_overview.png)

| Feature | Value |
|---|---|
| Length | 100 samples (400 ns) |
| Baseline | 146 ADC (samples 0–29) |
| Pulse onset | sample 30 (t = 120 ns) |
| Peak | sample 32, 1393 ADC |
| Rise time | 8 ns (cross → peak) |
| Tail | exponential decay; ≈ 20 ADC above baseline at sample 99 |

This is representative of a HyCal PbWO₄ signal: a fast leading edge
(≈ 10 ns) followed by a long PMT/scintillator tail.

## Waveform analyzer — `WaveAnalyzer`

### Algorithm summary

The analyzer is tuned for production data with realistic baseline
contamination. Each design choice below addresses a specific failure
mode of the naive "mean ± σ pedestal + first local maximum" recipe:

- **Median + MAD pedestal seed.** A previous-event tail or early
  ringing biases a simple-mean seed and inflates its σ-clip band.
  Median + MAD·1.4826 remains correct under up to 50 % contamination
  ([§ Robustness](#robustness-on-contaminated-pedestals)).
- **Adaptive pedestal window.** When the leading window flags as
  contaminated, the analyzer estimates the pedestal on the trailing
  samples and selects whichever has the lower RMS. The choice is
  recorded in `Q_PED_TRAILING_WINDOW`.
- **Per-channel pedestal quality bitmask.** `Q_PED_*` flags expose the
  pedestal trustworthiness without re-running the analyzer
  ([§ Pedestal quality](#pedestal-quality)).
- **Pedestal slope.** Least-squares slope on the surviving samples
  detects baseline drift that σ-clipping alone misses.
- **Noise-scaled local-maximum tolerance.** `trend()` uses
  `max(0.1, 0.5·rms)` to suppress spurious mini-peaks on quiet channels.
- **N-consecutive tail termination.** Integration ends only after
  `tail_break_n` consecutive sub-threshold samples, preventing
  premature truncation by single noise dips.
- **Inclusive integration bounds.** `peak.left`/`peak.right` are the
  outermost samples summed into `peak.integral`; the JS viewer
  shading uses the same convention.
- **Sub-sample peak time.** A 3-point quadratic vertex fit refines the
  peak time below the 4 ns sample grid.
- **Pile-up flagging.** Peaks whose integration windows lie within
  `peak_pileup_gap` of each other are flagged with `Q_PEAK_PILED`
  ([§ Peak quality](#peak-quality)).

The pipeline below walks through each step on the example trace;
figures and numerical values are reproduced by
[`scripts/plot_wave_analysis.py`](scripts/plot_wave_analysis.py).

### Parameters

Parameters are defined in `fdec::WaveConfig`
(`prad2dec/include/WaveAnalyzer.h`); defaults are shown.

| Field | Default | Unit | Role |
|---|---:|---|---|
| `smooth_order`     |   2 | order | Triangular kernel order; `1` disables smoothing, `N` gives a `2N − 1` tap kernel. |
| `peak_nsigma`      | 5.0 | × pedestal RMS | Sigma-scaled peak detection threshold above the local baseline. |
| `min_peak_height`  | 10.0 | ADC counts | Absolute floor on the detection threshold. Effective threshold is `max(peak_nsigma·rms, min_peak_height)`. |
| `min_peak_ratio`   | 0.3 | fraction | When two peaks overlap, a secondary survives only if its smoothed height is ≥ this fraction of the primary. |
| `int_tail_ratio`   | 0.1 | fraction | Integration stops when the pedsub waveform falls below `r · peak_height`. |
| `tail_break_n`     |   2 | samples | Required consecutive sub-threshold samples to terminate integration. |
| `peak_pileup_gap`  |   2 | samples | Two peaks within this many samples are both flagged `Q_PEAK_PILED`. |
| `ped_nsamples`     |  30 | samples | Window length for the pedestal estimate. |
| `ped_flatness`     | 1.0 | ADC counts | Floor on the σ-clip band: samples kept iff `|s − μ| < max(rms, ped_flatness)`. |
| `ped_max_iter`     |   3 | iterations | Maximum σ-clip passes; terminates early on mask convergence or fewer than 5 survivors. |
| `overflow`         | 4095 | ADC counts | 12-bit overflow value; peaks at this height are tagged. |
| `clk_mhz`          | 250.0 | MHz | Sample rate for time conversion `t = pos · 1000/clk_mhz` (ns). |

#### Configuration in `daq_config.json`

The `fadc250_waveform.analyzer` block in
[`database/daq_config.json`](../../../database/daq_config.json) overrides
any subset of the defaults. Replay, viewer servers, and filter paths
all consume the analyzer through `evc::DaqConfig::wave_cfg`, so a single
edit propagates to every consumer.

```jsonc
"fadc250_waveform": {
    "analyzer": {
        "smooth_order":     2,
        "peak_nsigma":      5.0,
        "min_peak_height":  10.0,
        "min_peak_ratio":   0.3,
        "int_tail_ratio":   0.1,
        "tail_break_n":     2,
        "peak_pileup_gap":  2,
        "ped_nsamples":     30,
        "ped_flatness":     1.0,
        "ped_max_iter":     3,
        "overflow":         4095,
        "clk_mhz":          250.0
    },
    "firmware": { /* see Firmware emulator */ }
}
```

### Pipeline

**1. Triangular smoothing.** With `smooth_order = N`, the kernel is
`buf[i] = (raw[i−1]·w + raw[i] + raw[i+1]·w) / (1 + 2w)` with
`w = 1 − 1/(N + 1)`. `N = 1` disables smoothing; the kernel half-width
is `N − 1`.

**2. Iterative pedestal.** On the first `ped_nsamples` samples of the
smoothed buffer:

- Seed `mean` and `rms` with the median and MAD·1.4826 of the window.
- Iterate up to `ped_max_iter` passes: drop samples deviating by more
  than `max(rms, ped_flatness)` from the mean; recompute statistics on
  the survivors.
- Record `nused` (surviving sample count), a `Q_PED_*` quality
  bitmask, and a least-squares slope (ADC/sample) on the survivors.

For the example trace: `mean = 145.61`, `rms = 0.45`, `nused = 28`,
`slope ≈ 0`, `quality = Q_PED_GOOD`.

**2b. Adaptive window.** If the leading window flags as suspicious
(`Q_PED_NOT_CONVERGED`, `Q_PED_TOO_FEW_SAMPLES`, `Q_PED_OVERFLOW`, or
`nused < ped_nsamples / 2`) and the buffer is long enough that the
trailing window does not overlap the leading one, the analyzer
estimates the pedestal on the trailing samples and adopts whichever
estimate has the lower RMS (with `nused` as tiebreaker). The trailing
choice is flagged with `Q_PED_TRAILING_WINDOW`.

**3. Threshold.** `thr = max(peak_nsigma · rms, min_peak_height)`.

**4. Local-maxima search.** The smoothed buffer is scanned for local
maxima. A peak is accepted if all three conditions hold:

- it is a local maximum (the `trend()` flat-tolerance scales with `rms`
  so noise wiggles do not fragment plateaus);
- its height above the local baseline (linear interpolation between
  surrounding minima) exceeds `thr`;
- its height above the pedestal mean exceeds `thr`.

**5. Integration.** Walking outward from the peak, the analyzer sums
pedsub samples until either the value drops below `tail_cut =
int_tail_ratio · peak_height` or below `rms`. Termination requires
`tail_break_n` consecutive sub-threshold samples; isolated dips are
held in a `pending` queue and either committed on recovery or
discarded at termination. `peak.left` and `peak.right` are the
inclusive integration bounds.

**6. Raw-position correction.** `pos` is recorded as the raw-sample
maximum near the smoothed peak so that `peak.height` reflects the
actual ADC value rather than a smoothed under-estimate.

**7. Sub-sample peak time.** A quadratic vertex through
`raw[pos − 1, pos, pos + 1]` yields
`δ = (h_{−1} − h_{+1}) / (2(h_{−1} − 2h_0 + h_{+1}))`, clamped to
`δ ∈ [−1, 1]` and applied only when the parabola is concave-down.
The reported time is `peak.time = (pos + δ) · 1000 / clk_mhz` (ns).
Time resolution improves from the 4 ns sample grid to ≪ 1 ns for
clean peaks.

**8. Pile-up flagging.** Each newly accepted peak is compared against
all previously accepted peaks; if either window comes within
`peak_pileup_gap` samples of the other, both peaks receive
`Q_PEAK_PILED`.

For the example trace:

| Field | Value |
|---|---|
| `peak.pos`   | 32 |
| `peak.time`  | 126.6 ns (`pos + δ = −0.36`) |
| `peak.height` | 1247 ADC |
| `[peak.left, peak.right]` | `[29, 48]` (20 samples, inclusive) |
| `peak.integral` | 9600 (ADC·sample) |
| `peak.quality` | `Q_PEAK_GOOD` |

![waveform](plots/fig3_soft_analysis.png)

### Peak quality

Each peak carries an 8-bit `quality` mask (defined in
`prad2dec/include/Fadc250Data.h`):

| Bit | Flag | Set when |
|---|---|---|
| `0`     | `Q_PEAK_GOOD`        | clean peak, no flags |
| `1<<0`  | `Q_PEAK_PILED`       | another peak's integration window lies within `peak_pileup_gap` samples; both peaks of the pair are flagged |
| `1<<1`  | `Q_PEAK_DECONVOLVED` | `height` and `integral` have been overwritten by the pile-up deconvolution path ([§ Pile-up deconvolution](#pile-up-deconvolution)) |

Downstream filters can require `peak.quality == Q_PEAK_GOOD` for a
clean-pulse subset, or accept `Q_PEAK_DECONVOLVED` peaks (where the
deconvolution provides corrected amplitudes) while excluding raw
piled peaks.

### Pedestal quality

The analyzer reports four scalars per channel describing the pedestal
estimate, written to the `events` tree as
`hycal.ped_{mean,rms,nused,quality,slope}`
(see [`docs/REPLAYED_DATA.md`](../../REPLAYED_DATA.md)):

| Field | Type | Use |
|---|---|---|
| `ped_mean`    | `float` | Pedestal mean after rejection |
| `ped_rms`     | `float` | RMS after rejection |
| `ped_nused`   | `uint8` | Surviving sample count (vs. `ped_nsamples`) |
| `ped_slope`   | `float` | LSQ drift (ADC/sample) on survivors |
| `ped_quality` | `uint8` | `Q_PED_*` bitmask |

| Bit | Flag | Set when |
|---|---|---|
| `0`     | `Q_PED_GOOD`             | clean estimate |
| `1<<0`  | `Q_PED_NOT_CONVERGED`    | `ped_max_iter` exhausted with kept-mask still moving |
| `1<<1`  | `Q_PED_FLOOR_ACTIVE`     | `rms < ped_flatness` (informational) |
| `1<<2`  | `Q_PED_TOO_FEW_SAMPLES`  | < 5 samples survived rejection |
| `1<<3`  | `Q_PED_PULSE_IN_WINDOW`  | `findPeaks` returned a peak with `pos` inside the pedestal window |
| `1<<4`  | `Q_PED_OVERFLOW`         | a raw window sample reached `cfg.overflow` |
| `1<<5`  | `Q_PED_TRAILING_WINDOW`  | adaptive window logic chose the trailing samples (informational) |

A clean-event filter is `ped_quality == 0`. To exclude only events
where the iterative cut failed,
`ped_quality & (Q_PED_NOT_CONVERGED | Q_PED_TOO_FEW_SAMPLES | Q_PED_PULSE_IN_WINDOW)` suffices.

### Robustness on contaminated pedestals

The example trace has a clean baseline; algorithm robustness becomes
relevant on contaminated baselines. The figure below uses a synthetic
30-sample window where the first 14 samples carry a 2.5 ADC bias
(simulating a previous-event tail) on top of the same ±0.4 ADC noise:

![robustness](plots/fig6_robustness.png)

- The simple-mean seed lies between the two clusters and pulls the
  σ-clip band along with it; iteration either locks onto the
  contaminated subset or rejects most samples and reports the biased
  seed.
- The median + MAD·1.4826 seed sits on the true baseline; a single
  σ-clip pass excludes the contaminated samples.

The contamination level is in the regime where σ-clip is most fragile
(close to the true baseline, more than half the band). Stronger
contamination is rejected by both seeds; the median's advantage lies
on the marginal cases that bias the energy resolution silently.

### Crowded windows and pile-up flagging

The 100-sample readout (400 ns) typically contains more than one
PbWO₄ pulse from accidentals, after-pulses, or beam-related
multi-hits. The figure below shows three synthetic pulses at samples
20, 35, 50 with heights 800/600/350 ADC:

![crowded](plots/fig7_crowded.png)

- Each pulse rises above its surrounding minima by more than `thr`
  and is found independently by the local-maxima search.
- The integration windows of adjacent pulses touch (the slow tail of
  pulse N runs into the rising edge of pulse N+1). Step 8 sets
  `Q_PEAK_PILED` on both peaks of every such pair.
- Downstream code may keep all peaks (the integrals remain correct
  under the tail-cutoff stopping rule) or filter on `Q_PEAK_GOOD`.

`Q_PEAK_PILED` is independent of the `min_peak_ratio` rejection
([§ Pipeline](#pipeline) step 4): the latter decides whether a peak
is recorded at all; the former describes whether a recorded peak has
a close neighbour.

### Pile-up deconvolution

The local-maxima + tail-cutoff integral biases both height and
integral on `Q_PEAK_PILED` peaks: each pulse's tail lifts the
apparent baseline of the next, which in turn distorts its peak
height. When a per-type pulse template is available for the channel's
category (PbGlass / PbWO4 / LMS / Veto), the analyzer recovers the
underlying per-pulse amplitudes by per-event Levenberg–Marquardt fits
of a parametric two-tau template.

**Model.** A piled-up event with `K` peaks at times `t_k` is modeled
as a non-negative linear combination of `K` template instances:

```
s(t_i) = Σ_k a_k · T(t_i − τ_k; τ_r_k, τ_f_k)
```

`T(·; τ_r, τ_f)` is the two-tau pulse fitted by
[`fit_pulse_template.py`](scripts/fit_pulse_template.py). For each
peak the LM frees `(a_k, t0_k, τ_r_k, τ_f_k)` with the per-type
template providing the initial guess and tight bounds; the rise/fall
shapes are constrained to `[τ/factor, τ·factor]` of the template
medians (`shape_window_factor`), and `t0_k` is constrained to lie
within `t0_window_ns` of the WaveAnalyzer-reported peak time.
Non-negativity is enforced through bounded amplitude constraints
(`a_k ∈ [0, amp_max_factor·peak.height]`). `K ≤ MAX_PEAKS = 8` in
production, so the LM solve is negligible per event.

**Per-type template store
([`PulseTemplateStore`](../../../prad2dec/include/PulseTemplateStore.h)).**
Loaded once at startup from
[`fit_pulse_template.py`](scripts/fit_pulse_template.py) output. Two
pieces are consumed:

1. The `_by_type` block — one (τ_r, τ_f) median per category
   (PbGlass / PbWO4 / LMS / Veto). These are the templates the
   deconvolver uses, after re-validation against the analyzer's
   `tau_*_range_ns` gates.
2. Each per-channel record's `module_type` field — used only to
   classify the channel. Per-channel τ_r/τ_f are deliberately
   ignored; pulse shapes group cleanly by crystal type, and a
   well-fit category aggregate is more reliable than thousands of
   single-channel fits with varying statistics.

`Lookup(roc_tag, slot, channel)` resolves the channel's category and
returns the matching template. Channels with `module_type = "Unknown"`
or absent from the JSON receive a `nullptr` and are skipped. File or
parse errors leave the store invalid; analyzers fall back silently to
the local-maxima heights.

**In-place output.** When deconvolution converges, the affected `Peak`
objects are overwritten:

- `peak.height`   ← `a_k · T_max`
- `peak.integral` ← `a_k · Σᵢ T(t_i − τ_k)` over the per-peak window
- `peak.quality |= Q_PEAK_DECONVOLVED`

Failure paths (template missing, template out of validity range, LM
non-convergence) leave the WaveAnalyzer values untouched.

**Configuration.** All knobs live under
`fadc250_waveform.analyzer.nnls_deconv` in
[`database/daq_config.json`](../../../database/daq_config.json) (the
historical name `nnls_deconv` is retained for compatibility; the
underlying solver is no longer NNLS):

```jsonc
"nnls_deconv": {
    "enabled":              false,
    "template_file":        "waveform/pulse_templates_024177.json",
    "apply_to_all_peaks":   false,
    "tau_r_range_ns":       [0.5, 15.0],
    "tau_f_range_ns":       [2.0, 100.0],
    "shape_window_factor":  1.5,
    "t0_window_ns":         8.0,
    "amp_max_factor":       2.0,
    "pre_samples":          8,
    "post_samples":         40
}
```

The shipped default is `enabled: false` while the per-type templates
remain under study. `enabled` and `apply_to_all_peaks` gate the
runtime path; `tau_*_range_ns`, `shape_window_factor`, `t0_window_ns`,
and `amp_max_factor` constrain the LM search; `pre_samples` and
`post_samples` set the per-peak integral window. `template_file` is
resolved against `PRAD2_DATABASE_DIR`.

**Wiring.** Every code path that uses `WaveAnalyzer` picks
deconvolution up automatically when (a) `daq_config` enables it,
(b) a `PulseTemplateStore` has been loaded, and (c) the analyzer is
bound to the store and the current channel key:

```cpp
// once at setup
fdec::PulseTemplateStore store;
if (cfg.wave_cfg.nnls_deconv.enabled
    && !cfg.wave_cfg.nnls_deconv.template_file.empty()) {
    store.LoadFromFile(db_dir + "/" + cfg.wave_cfg.nnls_deconv.template_file,
                       cfg.wave_cfg.nnls_deconv);
}
fdec::WaveAnalyzer ana(cfg.wave_cfg);
ana.SetTemplateStore(&store);

// per channel inside the existing loop
ana.SetChannelKey(roc.tag, slot, channel);
ana.Analyze(samples, n, wres);   // auto-deconv when conditions met
```

Per-peak deconvolved values are written into the same `peak.height`
and `peak.integral` fields the local-maxima path uses; consumers
distinguish the two paths via `Q_PEAK_DECONVOLVED`.

For visualising deconvolution on individual events, see
[`deconv_pileup_demo.py`](../../../analysis/pyscripts/deconv_pileup_demo.py).
The Python binding `wave_ana.deconvolve()` runs whenever a valid
template is supplied, regardless of the production `enabled` flag.

### Worked example — synthetic 3-pulse pile-up

The same input as
[fig 7](#crowded-windows-and-pile-up-flagging) — three PbWO₄-like
pulses at samples 20/35/50 with truth heights 800/600/350 ADC on a
146 ADC baseline. The script
[`scripts/plot_wave_deconv.py`](scripts/plot_wave_deconv.py) calls
`fdec::WaveAnalyzer::Deconvolve` through the bindings; every plotted
amplitude originates from the C++ solver.

![synth-deconv](plots/fig8_deconv_synth.png)

- **Open triangles** — `WaveAnalyzer.peaks[k].height` from the
  pre-deconv local-maxima path. Tail bleed-through biases peak 1 (the
  trailing shoulder of pulse 0 lifts its baseline) and peak 2
  (twice-piled).
- **Filled circles** — `DeconvOutput.height[k] = a_k · T_max(τ_r_k,
  τ_f_k)` from the LM fit, initialised on the per-type template
  `(τ_r = 10 ns, τ_f = 48 ns)`.
- **Black crosses** — the synthetic input heights (truth).
- **Orange curve** — `Σ a_k · T(t − t0_k; τ_r_k, τ_f_k)` evaluated at
  the design-matrix sample times; passes through every filled
  circle by construction.

Numerical comparison (script stdout):

| Peak | Truth | WA height | Deconv height | Δ(WA) | Δ(deconv) |
|---:|---:|---:|---:|---:|---:|
| 0 | 800 | ≈ 800 | ≈ 800 | small | small |
| 1 | 600 | biased high | ≈ 600 | large | small |
| 2 | 350 | biased high | ≈ 350 | large | small |

Deconvolved heights agree with truth to within ≪ 1 % on the leading
peak and a few percent on the buried peaks; the local-maxima method
shows tens of percent bias on peaks 1 and 2.

### Worked example — real PbWO₄ pile-up event

The same script with `--evio <path>` scans a production EVIO file for
the first physics event with a `Q_PEAK_PILED` peak on a channel whose
module type is covered by the per-type template store and whose LM
deconvolution converges:

```bash
cd docs/technical_notes/waveform_analysis
python scripts/plot_wave_deconv.py --evio /data/evio/data/prad_024202/prad_024202.evio.00000
```

![real-deconv](plots/fig9_deconv_real.png)

The legend matches fig 8. Truth markers are absent in real data; the
metric of interest is the magnitude of the WA → deconv shift on the
trailing peak, which signals pile-up bias on the local-maxima
estimate. The checked-in figure was generated from run 024202; rerun
the script against any local EVIO file to regenerate against your own
data.

### Parameter sensitivity

Two parameters visibly change the analyzer's output on the example
trace:

![params](plots/fig4_soft_parameters.png)

**Pedestal `ped_flatness` × `ped_max_iter`.** All 30 baseline samples
enter pass 1; samples deviating from the running mean by more than
`max(rms, ped_flatness) = 1.0` are dropped, and the procedure
repeats. The band collapses onto the dominant cluster (146/147 ADC),
rejecting outliers (143/144/145/150). `ped_flatness` floors the
band to prevent over-tight clipping on quiet channels.

> The demo runs the procedure on the raw samples for plot
> readability; the C++ runs it on the smoothed buffer, accounting
> for the 0.9 ADC difference in the converged mean. The
> kept/rejected pattern is identical.

**`int_tail_ratio`.** Integration walks outward from the peak and
stops when the pedsub waveform falls below `r · peak_height`:

| `int_tail_ratio` | Window | Samples | Integral |
|---:|:---:|---:|---:|
| 0.20 | [30, 41] | 12 | 8376 |
| 0.10 (default) | [30, 47] | 18 | 9477 |
| 0.05 | [30, 57] | 28 | 10332 |

The default `0.10` captures the prompt peak plus the first ≈ 70 ns
of the tail. Smaller values recover more tail energy at the cost of
sensitivity to baseline drift and downstream pulses.

**Smoothing — `smooth_order`.** Smoothing is invisible at the scale
of a 1247 ADC peak; it matters on small-signal channels where
per-sample fluctuations are comparable to the pulse height. The
figure below uses a small ≈ 24 ADC bump on a baseline with ±3 ADC
noise:

![smoothing](plots/fig5_smoothing.png)

| `smooth_order` | Spurious local maxima above +2 ADC | Smoothed peak height |
|---:|---:|---:|
| 1 (raw)       | 6 | 169 |
| 2 (default)   | 3 | 166 |
| 4             | 1 | 162 |

`smooth_order = 2` collapses the zig-zag without measurable peak
attenuation. `smooth_order = 4` removes essentially all baseline
structure but clips the peak by ≈ 7 ADC and is appropriate only for
very low-S/N channels.

The remaining parameters affect this trace only marginally:

- `peak_nsigma · rms = 5 · 0.45 = 2.25 ADC` lies below
  `min_peak_height = 10`, so the floor wins on this trace; on noisier
  channels (`rms ≳ 2 ADC`) the `peak_nsigma` rule activates.
- `min_peak_ratio` activates only when peaks share an integration
  range.

## Firmware emulator — `Fadc250FwAnalyzer`

The firmware analyzer reproduces the on-board pipeline so that
offline analysis can be cross-checked against firmware-reported
values without re-running the DAQ. The full algorithm specification
appears in
[`docs/clas_fadc/FADC250_algorithms.md`](../../clas_fadc/FADC250_algorithms.md);
this section walks through the parameters.

### Parameters

Parameters live under the `fadc250_waveform.firmware` block in
[`database/daq_config.json`](../../../database/daq_config.json).
`NSB` and `NSA` are in nanoseconds, floored to whole 4 ns samples by
the analyzer; remaining fields are unitless or in ADC counts.

| Field | Unit | Role |
|---|---|---|
| `TET`                    | ADC counts | Trigger Energy Threshold above pedestal; pulse rejected if `Vpeak − Vmin ≤ TET`. |
| `NSB`                    | ns         | Window before threshold crossing (Mode 2 integral); floored to whole samples. |
| `NSA`                    | ns         | Window after threshold crossing; same flooring as `NSB`. |
| `NPEAK` (= `MAX_PULSES`) | —          | Max pulses kept per channel per readout (1–4). |
| `NSAT`                   | samples    | Required consecutive-above-`TET` samples after `Tcross`; rejects single-sample spikes. `NSAT = 1` reproduces the legacy single-sample crossing. |
| `NPED`                   | samples    | Number of leading samples summed for the `Vnoise` estimate. |
| `MAXPED`                 | ADC counts | Online outlier rejection on the `Vnoise` sum; `0` disables. |
| `CLK_NS`                 | ns         | Sample period (4 ns at 250 MHz). |

Run-config defaults:

```jsonc
"fadc250_waveform": {
    "firmware": {
        "TET": 10.0, "NSB": 8, "NSA": 128, "NPEAK": 1,
        "NSAT": 4, "NPED": 3, "MAXPED": 1, "CLK_NS": 4.0
    }
}
```

### Pipeline

Step-by-step on the example waveform:

**1. Pedestal estimate (`Vnoise`).** Mean of the first `NPED = 3`
samples with `MAXPED = 1` outlier filter (drop any sample whose
deviation from the running mean exceeds 1 ADC). For the trace:
`(146 + 147 + 144) / 3 = 145.67`; sample 1 (147) is filtered;
refined mean is `145.0`.

**2. Pulse search.** `Vmin = Vnoise`. The buffer is scanned starting
at sample `NPED`; the first pulse is detected as soon as a sample
exceeds `Vnoise` and rises monotonically to a local maximum.

**3. Acceptance.** `Vpeak = 1393`; pedsub height is `1393 − 146 =
1247 ≫ TET = 10` → accepted.

**4. `Tcross`.** First leading-edge sample whose pedsub value exceeds
`TET`: sample 30, since `637 − 146 = 491 > 10`.

**5. NSAT gate.** With `NSAT = 4`, samples 30, 31, 32, 33 must all
exceed `TET`. They are (491, 1221, 1247, 1093) → accepted.
`NSAT = 1` reduces this gate to a no-op.

**6. TDC — `Va`, bracket, fine time.**

```
Va  = Vmin + (Vpeak − Vmin) / 2 = 146 + (1393 − 146)/2 = 769.5
```

Find the bracket on the rising edge: smallest `k` with `s[k] ≥ Va`.
Here `s[30] = 637 < 769.5`, `s[31] = 1367 ≥ 769.5` → `k = 31`,
giving `Vba = 637`, `Vaa = 1367`. Fine time:

```
fine = round( (Va − Vba) / (Vaa − Vba) · 64 )
     = round( (769.5 − 637) / (1367 − 637) · 64 )
     = round( 0.1815 · 64 ) = 12

coarse     = k − 1 = 30
time_units = coarse · 64 + fine = 1932    (LSB = 62.5 ps)
time_ns    = time_units · CLK_NS / 64 = 120.75 ns
```

In fig 2 (left), the dot-dash `Va` line crosses the rising edge
between the diamond `Vba` and square `Vaa` markers; the fine-time
arrow points from the `Vba` sample to the interpolated zero
crossing.

**7. Mode-2 integral.** Window
`[cross − NSB_s, cross + NSA_s]` with `NSB_s = NSB / 4 = 2`,
`NSA_s = NSA / 4 = 32`, giving `[28, 62]` (35 samples = 140 ns).
The integrand is `s' = max(0, s − Vnoise)`:

```
Σ s'[28..62] = 10589  (pedsub ADC·sample)
```

The shaded band in fig 2 (right) is this sum.

![firmware](plots/fig2_firmware_analysis.png)

**8. Quality bitmask.** `0x00 = Q_DAQ_GOOD` for the example. Set bits
indicate:

| Bit | Flag | Condition |
|---|---|---|
| `1<<0` | `Q_DAQ_PEAK_AT_BOUNDARY` | peak at the last sample |
| `1<<1` | `Q_DAQ_NSB_TRUNCATED`   | `cross − NSB_s < 0`; window clipped |
| `1<<2` | `Q_DAQ_NSA_TRUNCATED`   | `cross + NSA_s ≥ N`; window clipped |
| `1<<3` | `Q_DAQ_VA_OUT_OF_RANGE` | `Va` not bracketed on the rising edge |

## Side-by-side comparison

| Field | Waveform (`WaveAnalyzer`) | Firmware (`Fadc250FwAnalyzer`) |
|---|---|---|
| Pedestal | 145.61 ± 0.45 (28/30 samples; median+MAD seed, σ-clip) | 145.0 (3 samples; `MAXPED` filter) |
| Time | 126.6 ns (peak vertex, sub-sample) | 120.75 ns (mid-amplitude rising edge) |
| Height | 1247 ADC | 1247 ADC |
| Integral window | [29, 48] (20 samples, tail-driven, inclusive) | [28, 62] (35 samples, fixed `NSB/NSA`) |
| Integral | 9600 | 10589 |

Time observables differ by construction: `WaveAnalyzer.peak.time` is
the peak-vertex time (sub-sample, ≪ 1 ns precision for clean peaks);
the firmware time is mid-amplitude on the rising edge (62.5 ps LSB).
The firmware's wider window captures more of the slow scintillation
tail; with `NSA = 128 ns` the integral window stops at sample 62,
excluding samples 63–99.

## Reproducing the plots

The figures and numerical values in this note are produced by
scripts that drive the C++ analyzers via the `prad2py` Python
bindings; no parallel Python re-implementation can drift from the
production code. Build with `-DBUILD_PYTHON=ON` so that
`import prad2py` resolves, then run from this directory:

```bash
cd docs/technical_notes/waveform_analysis
python scripts/plot_wave_analysis.py                   # figs 1–7
python scripts/plot_wave_deconv.py [--evio <path>]     # figs 8–9
```

The only Python re-implementations in the scripts are pedagogical:
a small σ-clip helper used by the robustness figure to demonstrate
the simple-mean seed (the binding only exposes the median + MAD
seed), and the two-tau evaluation
`(1 − exp(−t/τ_r))·exp(−t/τ_f)` for the model-trace overlay in
figs 8–9. Heights, integrals, and convergence flags all originate
from the C++ solver.

## See also

- [`docs/clas_fadc/FADC250_algorithms.md`](../../clas_fadc/FADC250_algorithms.md)
  — full firmware algorithm specification with manual cross-references.
- [`prad2dec/include/WaveAnalyzer.h`](../../../prad2dec/include/WaveAnalyzer.h),
  [`WaveAnalyzer.cpp`](../../../prad2dec/src/WaveAnalyzer.cpp) — C++ source.
- [`prad2dec/include/Fadc250FwAnalyzer.h`](../../../prad2dec/include/Fadc250FwAnalyzer.h),
  [`Fadc250FwAnalyzer.cpp`](../../../prad2dec/src/Fadc250FwAnalyzer.cpp) — C++ source.
- [`prad2dec/include/Fadc250Data.h`](../../../prad2dec/include/Fadc250Data.h)
  — `Peak` / `Pedestal` structs and `Q_PEAK_*` / `Q_PED_*` / `Q_DAQ_*` flag definitions.
- [`docs/REPLAYED_DATA.md`](../../REPLAYED_DATA.md) — replay-tree branch layout.
