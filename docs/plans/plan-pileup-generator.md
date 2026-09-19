# Plan — `pileup_generator.py`

**Status:** draft (iteration 2, pre-approval — planning-critic feedback incorporated)
**Owner:** Jingyi Zhou (pile-up deconvolution study, PRad-II PbWO4 HyCal modules)
**Roadmap link:** item 2 from `jinst.md` — *"Pile-up generator. Poisson-random + correlated (cluster-coincident) pile-up. Configurable noise model: Gaussian + measured pedestal autocorrelation. Quote rate scaling vs PRad-II luminosity."*

> **Scope disclaimer.** v0.1 covers only a *piece* of roadmap item 2: a **deterministic single-channel response-surface generator** on a `(ΔT, ratio)` grid. Poisson-random pile-up timing, cluster-coincident cross-channel pile-up, and pedestal-autocorrelated noise are **explicitly deferred** (see §13). v0.1 produces the inputs for the "amplitude error vs ΔT vs ratio" plot; it does **not** produce the "Poisson-random pile-up ROC vs luminosity" plot.

**Target path:** `analysis/pyscripts/pileup_generator.py`
**Scope covered:** first-version, PbWO4-only, 2-pulse (base + 1 injected) synthetic events with an explicit no-extra-noise default (base pulse already carries real detector noise).
**Research inputs used:**
- `docs/technical_notes/waveform_analysis/wave_analysis.md` (`WaveAnalyzer` C++ API + `Q_PEAK_*` flags)
- `docs/technical_notes/waveform_analysis/fit_pulse_template.md` (§2b physics-selection defaults: `--height-min 500`, `--t0-min 25.0`; per-material `(τ_r, τ_f)` medians; template JSON schema)
- `analysis/pyscripts/fit_pulse_template.py` (base-event selection precedent: `_common.setup_pipeline`, `WaveAnalyzer.analyze`, `WaveAnalyzer.fit_pulse_shape`, clean-pulse cuts, pedestal accounting, multi-input EVIO discovery at lines 661–684)
- `analysis/pyscripts/deconv_pileup_demo.py` (downstream consumption pattern: `PulseTemplateStore.load_from_file`, `wave_ana.deconvolve`, per-peak `amplitude / height / t0_ns / tau_r_ns / tau_f_ns` outputs; `height[k] = amplitude[k] * T_max`)
- `prad2dec/include/PulseTemplateStore.h` (API: `type_template("PbWO4")` for per-type lookup)
- `analysis/pyscripts/output/pulse_templates_025308_h500_t0cut.json` (`_by_type["PbWO4"]`: `τ_r = 2.24 ns`, `τ_f = 24.0 ns`, `t0 = 26.3 ns`)
- `database/daq_config.json` (`clk_mhz = 250` → `4 ns/sample`; standard FADC250 window ≈ 200 samples = 800 ns)

**Open questions blocking approval:** none technical — a small number of design defaults (listed under §14 *Design decisions to confirm*) are called out for explicit sign-off before implementation.

---

## 1. Purpose recap

Produce a bag of synthetic pile-up waveforms where **every injected pulse has known truth** (`t0`, `height`, template shape). The downstream user runs `dec.WaveAnalyzer.deconvolve()` on each synthetic event and compares recovered heights and onsets to truth to produce the headline paper plot: amplitude error vs. ΔT vs. amplitude ratio.

The generator does *not* run any deconvolution itself, does *not* modify C++ code, and does *not* produce ROC curves. Those are separate deliverables.

Each synthetic event is (default noise mode):

```
synthetic[i] = base_waveform[i]        (real EVIO, real pedestal, real noise, real single pulse)
             + inject_pulse            (template shape, chosen height, chosen onset t0)
             (no additional readout noise added — base already carries it)
```

Under `--noise-model gaussian-iid` (opt-in stress test only), an extra `N(0, σ)` iid readout-noise term is added on top; that mode is documented as artificially increasing baseline RMS.

## 1a. Time and amplitude conventions (READ FIRST)

Two conventions matter for the downstream comparison against `dec_out`. Both are chosen to match what the C++ deconvolver reports.

### 1a.1 Time — peak time vs. template onset

The two-τ template is parameterised by an **onset** `t0`:

```
T(t; t0, τ_r, τ_f) = (1 - exp(-(t-t0)/τ_r)) · exp(-(t-t0)/τ_f) / T_max     for t > t0
```

The template *peak* sits at a fixed positive offset past the onset:

```
t_peak_offset = τ_r · log((τ_r + τ_f) / τ_r)
              ≈ 2.24 · log(26.24 / 2.24)  ≈ 5.5 ns    (PbWO4)
```

`WaveAnalyzer.analyze` returns `pk.time` = the sample-position of the **peak** (interpolated), *not* the onset. The C++ deconvolver `wave_ana.deconvolve` reports `dec_out.t0_ns` = the fitted **template onset**. These are two different quantities:

```
pk.time  ≈  fitted_onset_t0  +  t_peak_offset
```

**Definition of ΔT in this generator.** ΔT is **peak-to-peak separation** — the interval between the peak of the base pulse and the peak of the injected pulse. This is the physically intuitive quantity users draw on a scope. It equals the onset-to-onset separation because both pulses share the same `(τ_r, τ_f)` and therefore the same `t_peak_offset`.

**Storage.** Each synthetic event stores BOTH:
- `truth_peak_times_ns[:, k]` — peak time of pulse k (base and injected), for comparison to `pk.time` from `WaveAnalyzer.analyze`.
- `truth_onset_t0_ns[:, k]` — template onset of pulse k, for comparison to `dec_out.t0_ns` from the deconvolver.

Relationship: `truth_peak_times_ns = truth_onset_t0_ns + t_peak_offset` (stored in metadata for reference).

**Downstream rule.** When benchmarking the deconvolver's timing, compare `dec_out.t0_ns` against **`truth_onset_t0_ns`**, *not* against `truth_peak_times_ns`.

### 1a.2 Amplitude — peak height vs. model coefficient

The C++ deconvolver reports two related amplitude quantities per fitted peak:
- `dec_out.amplitude[k]` — the unnormalized model coefficient `a_k` fitted to the *raw* two-τ shape (before dividing by `T_max`).
- `dec_out.height[k]` — the physical peak height in ADC, equal to `a_k · T_max`.

This generator injects a **peak-normalised** pulse: `inject_wave = h_inject · T(t; t0, τ_r, τ_f)` where `T` already includes the `/T_max` normalization, so `h_inject` **is** the peak height in ADC.

**Storage.** Truth uses the **height** convention:
- `truth_heights_adc[:, k]` — physical peak height in ADC of pulse k. This is the primary truth column.
- `truth_model_amplitudes[:, k]` = `truth_heights_adc[:, k] / T_max` — stored as a convenience for anyone who wants to compare directly against `dec_out.amplitude`.

**Downstream rule.** Compare `dec_out.height` against **`truth_heights_adc`**. If you compare `dec_out.amplitude` against something, use `truth_model_amplitudes` — not `truth_heights_adc`.

## 2. File structure — top-level flow

Single script, ~600–750 LOC (see §12), organised as:

```
pileup_generator.py
├── module docstring (usage, physics motivation, output schema, §1a conventions)
├── constants: DEFAULT_DT_GRID_NS, DEFAULT_RATIO_GRID, DEFAULT_MATERIAL, MAX_ADC,
│              DEFAULT_HEIGHT_MIN=500, DEFAULT_T0_MIN=25.0, DEFAULT_CHI2_MAX
├── dataclasses:
│     BaseEvent      — one accepted base waveform + WaveAnalyzer + fit-shape diagnostics
│     GridPoint      — one (dt_ns, ratio) coordinate with its output buffer
│     TemplateShape  — resolved (tau_r, tau_f, t0_template, T_max, t_peak_offset)
├── template loading:
│     load_template_shape(json_path, material) -> TemplateShape
│       (computes T_max and t_peak_offset once, from tau_r/tau_f medians)
├── base-event selection:
│     iterate_base_events(pipeline, args, shape) -> Generator[BaseEvent]
│       (reuses _common.setup_pipeline + WaveAnalyzer.analyze + WaveAnalyzer.fit_pulse_shape;
│        applies clean-pulse cuts, fit-quality gate, and residual hidden-pileup veto)
├── injection kernel:
│     eval_template(t_ns, t0_ns, tau_r, tau_f, T_max) -> np.ndarray  (float, peak = 1)
│     inject_pulse(base_samples_u16, ped, rms, base_height, base_peak_time_ns,
│                  base_onset_t0_ns, dt_ns, ratio, shape, clk_ns,
│                  noise_model, rng) -> (u16 waveform, truth)
├── scanning driver:
│     scan_grid(base_events, dt_grid, ratio_grid, n_per_config, ...)
│       — paired-base design: each accepted BaseEvent is reused across ALL
│         grid cells, injecting a different (dt, ratio) per cell.
├── serialisation:
│     write_grid_point(out_dir, dt_ns, ratio, waveforms, truth, metadata_dict)
│       — writes .npz + sidecar .json (no allow_pickle needed)
├── progress reporter (mirrors fit_pulse_template.py:_emit_progress)
├── validate_output(npz_path, wave_ana=None)
│     — schema-only if wave_ana is None; full peak-detect smoke test if given
└── main() — argparse + orchestration
```

The main loop is:

1. Multi-input EVIO discovery — mirror `fit_pulse_template.py` lines 661–684 **exactly**:
   - Collect `args.evio_paths` (nargs='+').
   - For each, `C.discover_split_files(inp)`; dedupe order-preserving.
   - `SystemExit` if empty.
   - Call `setup_pipeline(evio_path=args.evio_paths[0], ...)` then override `p.evio_files = all_files`.
2. `load_template_shape(--template, --material)` → precompute `T_max` and `t_peak_offset`.
3. Build `GridPoint` objects for every `(dt_ns, ratio)` cell in the requested grid.
4. Walk EVIO with `iterate_base_events`; for each accepted `BaseEvent`, hand it to `scan_grid`, which reuses that single base across every `(dt, ratio)` cell (see §7).
5. Stop iterating EVIO once every cell has ≥ `n_per_config` events (or EVIO exhausts — then warn).
6. Write one `.npz` + one sidecar `.json` per grid cell via `write_grid_point`; write a top-level `manifest.json`.
7. Run `validate_output` on each written file (schema-only by default; add `wave_ana=` for the smoke test at the end of `main`).

## 3. Data model

### 3.1 In-memory `BaseEvent`

```python
@dataclass
class BaseEvent:
    samples: np.ndarray        # uint16, shape (n_samples,) — raw ADC, unmodified
    ped_mean: float            # from wave_ana.analyze
    ped_rms:  float            # noise σ; drives injected-pulse noise σ if --noise-model gaussian-iid
    base_height: float         # pk.height (ADC, pedsub) of the single found peak
    base_peak_time_ns: float   # pk.time (ns) — WaveAnalyzer's peak-position estimate
    base_onset_t0_ns: float    # fitted template onset from WaveAnalyzer.fit_pulse_shape
    base_peak_sample: int      # pk.pos (integer sample index)
    base_chi2_per_dof: float   # from fit_pulse_shape — used for the fit-quality gate
    base_fit_converged: bool   # from fit_pulse_shape
    # Provenance — one channel, one physics event, per accepted base
    base_roc_tag: int
    base_slot: int
    base_chan: int
    base_channel_id: str       # "<roc_tag>_<slot>_<chan>"
    base_channel_name: str     # e.g. "W123"
    module_type: str           # "PbWO4" (filtered by --material)
    base_source_file: str      # basename of the EVIO split
    base_physics_event_index: int   # physics-event counter within the run (for reproducibility)
    base_run_number: int
```

### 3.2 Per-grid-cell output `.npz`

Filename: `pileup_<material>_dt<DDD>_ratio<RRRR>.npz`
- `DDD` = ΔT in ns, zero-padded to 3 digits (`016`, `048`, …).
- `RRRR` = amplitude ratio × 1000, zero-padded to 4 digits (`0500`, `1000`, `2000`) — integer form keeps filenames sortable and free of `.`.

Example: `pileup_PbWO4_dt016_ratio0500.npz`. A sidecar `pileup_PbWO4_dt016_ratio0500.json` sits next to it (§3.4).

Arrays inside (all fixed shape per file so downstream `np.load(..., allow_pickle=False)['waveforms']` needs no ragged handling and no pickle):

| key                        | dtype   | shape           | meaning |
|----------------------------|---------|-----------------|---------|
| `waveforms`                | uint16  | (N, n_samples)  | synthetic ADC samples, ready for `wave_ana.analyze` |
| `truth_heights_adc`        | float32 | (N, 2)          | peak height in ADC. col 0 = injected, col 1 = base_reference. Compare to `dec_out.height`. |
| `truth_model_amplitudes`   | float32 | (N, 2)          | `truth_heights_adc / T_max` — convenience for comparing to `dec_out.amplitude`. |
| `truth_onset_t0_ns`        | float32 | (N, 2)          | template-onset time (ns). col 0 = injected, col 1 = base_reference. Compare to `dec_out.t0_ns`. |
| `truth_peak_times_ns`      | float32 | (N, 2)          | peak-position time (ns). col 0 = injected, col 1 = base_reference. Compare to `pk.time` from `analyze`. |
| `truth_source`             | `<U9`   | (2,)            | `["injected", "base_ref"]` — labels the two columns explicitly. |
| `pedestals`                | float32 | (N, 2)          | col 0 = ped mean, col 1 = ped rms |
| `base_channel_names`       | `<U8`   | (N,)            | channel name of the source base event |
| `base_channel_ids`         | `<U16`  | (N,)            | `"<roc>_<slot>_<chan>"` |
| `base_roc_tags`            | int32   | (N,)            | source-base ROC tag |
| `base_slots`               | int32   | (N,)            | source-base slot |
| `base_chans`               | int32   | (N,)            | source-base channel |
| `base_run_numbers`         | int32   | (N,)            | source-base run |
| `base_source_files`        | `<U64`  | (N,)            | source-base EVIO split basename |
| `base_physics_event_index` | int32   | (N,)            | source-base physics-event index |
| `base_chi2_per_dof`        | float32 | (N,)            | fit diagnostic recorded at base-selection time |
| `meta_dt_ns`               | float32 | ()              | mirror of sidecar `dt_ns` |
| `meta_ratio`               | float32 | ()              | mirror of sidecar `amplitude_ratio` |
| `meta_tau_r_ns`            | float32 | ()              | mirror of sidecar `tau_r_ns` |
| `meta_tau_f_ns`            | float32 | ()              | mirror of sidecar `tau_f_ns` |
| `meta_T_max`               | float32 | ()              | mirror of sidecar `T_max` |
| `meta_t_peak_offset_ns`    | float32 | ()              | mirror of sidecar `t_peak_offset_ns` |
| `meta_n_events`            | int32   | ()              | mirror of sidecar `n_events` |

**Column-order convention (documented in module docstring and enforced by `truth_source`):** column 0 = **injected** (exact truth); column 1 = **base_reference** (measured observables from WaveAnalyzer/fit_pulse_shape — accurate but not exact truth). See §7a for why this distinction matters.

**No pickle.** All arrays above use numeric or fixed-width string dtypes; downstream code can call `np.load(path, allow_pickle=False)`. Rich metadata lives in the sidecar JSON (§3.4).

### 3.3 Metadata format — sidecar JSON (no pickle)

For each `.npz`, `write_grid_point` also writes a sidecar `.json` with the same basename. This avoids `allow_pickle=True` on load and makes metadata inspectable with `jq`:

```json
{
    "material":          "PbWO4",
    "dt_ns":             16.0,
    "amplitude_ratio":   0.5,
    "tau_r_ns":          2.2356,
    "tau_f_ns":          24.009,
    "T_max":             0.7273,
    "t_peak_offset_ns":  5.51,
    "template_t0_ns":    26.279,
    "template_source":   ".../pulse_templates_025308_h500_t0cut.json",
    "n_events":          50,
    "n_samples":         200,
    "clk_ns":            4.0,
    "noise_model":       "none",
    "noise_sigma_source": "n/a (base pulse carries the readout noise)",
    "seed":              12345,
    "generator_version": "pileup_generator.py v0.1",
    "generated_utc":     "2026-09-17T…",
    "base_event_source": {
        "runs":              [25308],
        "n_evio_splits":     12,
        "cut_height_min":    500.0,
        "cut_height_max":    2000.0,
        "cut_height_rms_mult": 10.0,
        "cut_t0_min_ns":     25.0,
        "cut_chi2_per_dof_max": 5.0,
        "residual_veto_sigma": 5.0
    },
    "column_convention": {
        "truth_columns":   ["injected", "base_ref"],
        "notes":           "col 0 = exact injection truth; col 1 = base_reference observables from WaveAnalyzer + fit_pulse_shape (accurate but not exact)."
    }
}
```

Alternative accepted: if the user prefers a single-file output, store the same JSON as a **string scalar** inside the `.npz` under key `metadata_json` — `str(np.load(...)['metadata_json'])` decodes with `json.loads(...)` and still needs no pickle. The plan writes a sidecar `.json` by default; the string-scalar alternative can be enabled with `--metadata-in-npz` if desired.

### 3.4 Top-level `manifest.json`

Written once next to the `.npz` files:

```json
{
  "generator_version": "pileup_generator.py v0.1",
  "generated_utc":     "…",
  "material":          "PbWO4",
  "template_source":   "…/pulse_templates_025308_h500_t0cut.json",
  "clk_ns":            4.0,
  "n_samples":         200,
  "dt_grid_ns":        [8, 12, 16, …, 128],
  "ratio_grid":        [0.1, 0.2, 0.3, 0.5, 1.0, 2.0],
  "n_per_config":      50,
  "noise_model":       "none",
  "T_max":             0.7273,
  "t_peak_offset_ns":  5.51,
  "files": [
    {"dt_ns": 8,  "ratio": 0.1, "path": "pileup_PbWO4_dt008_ratio0100.npz",
     "sidecar": "pileup_PbWO4_dt008_ratio0100.json", "n_written": 50},
    …
  ],
  "seed":              12345,
  "cli":               "<argv joined>"
}
```

## 4. Base-event selection

An EVIO channel-event is accepted as a base iff **all** of:

1. `wave_ana.analyze(samples)` returns exactly one peak (`len(peaks) == 1`).
2. `pk.quality == 0` — no `Q_PEAK_PILED`, no other soft-analyzer quality flag set.
3. `not pk.overflow`.
4. `args.height_min ≤ pk.height ≤ args.height_max` — the physics regime. **Default `500 ≤ h ≤ 2000` ADC**, matching the deployed physics-template selection (`fit_pulse_template.md` §2b, `--height-min 500`).
5. `pk.height ≥ args.height_rms_mult * ped_rms` — matches `fit_pulse_template.py`.
6. Module type matches `--material` (default `PbWO4`) via `p.hycal.module_by_daq(crate, s, c).type.name`.
7. **Fit-quality gate (new).** Run `wave_ana.fit_pulse_shape(...)` on the base pulse window. Require:
   - fit converged (`fit_result.converged` true),
   - fitted onset `t0_ns ≥ args.t0_min` (**default 25.0 ns**, matching `fit_pulse_template.md` §2b `--t0-min 25.0` — Population A physics window),
   - `t0_ns` inside the readout window (`0 < t0_ns < n_samples * clk_ns`),
   - `chi2_per_dof ≤ args.chi2_max` (default 5.0 — loose enough to admit real pulses, tight enough to reject deformed ones).
   Store `base_chi2_per_dof` and `base_fit_converged` on the `BaseEvent`.
8. **Residual hidden-pileup veto (new).** Subtract the fitted single-pulse model from the base samples. Outside a `±window_ns` bracket around the main peak (default `±3·τ_f ≈ 72 ns`), scan for residual samples exceeding `args.residual_veto_sigma * ped_rms` (**default 5·σ**). If any such sample exists, reject — a "clean" single-peak event with a residual bump is a hidden pile-up that would contaminate the truth.
9. **Injection headroom cut.** Given the ΔT grid, the largest ΔT is `dt_max_ns`. The base peak sits at sample `pk.pos`. The injected peak will land near sample `pk.pos + round(dt_max_ns / clk_ns)`. We require that sample index plus a `post_pad` (default 40 samples ≈ 5 fall-time constants for PbWO4) to stay strictly less than `samples.shape[0]`. Base events that fail this are rejected — never truncated. Checked once against `dt_max_ns`, not per grid cell (each base is reused across all cells; see §7).
10. **Saturation rejection (not clipping).** For the largest requested ratio `ratio_max`, the peak of the composite would be near `ped_mean + base_height + ratio_max * base_height`. If that exceeds `MAX_ADC - safety` (default `MAX_ADC = 4095`, `safety = 50`), reject the base. **Never clip** — a clipped injected pulse is truth we no longer possess. Validation later re-asserts `waveforms.max() < MAX_ADC - safety_margin` (§10).

Cuts 1–5 mirror the deployed physics-template selection; cuts 7–8 tighten the base further so we only inject on top of pulses whose shape and cleanliness we've explicitly validated on this event. Cuts 9–10 preserve the "known truth" property.

Each cut increments a named counter (`n_reject_multi_peak`, `n_reject_quality`, …); the `[done]` line prints the full rejection histogram.

## 5. Template evaluation and injection

### 5.1 Formula

The template is the same `two_tau_unit` used in `fit_pulse_template.py`:

```
T(t; t0, τ_r, τ_f) = (1 - exp(-(t - t0)/τ_r)) · exp(-(t - t0)/τ_f) / T_max     for t > t0
                    0                                                            otherwise

T_max          = (1 - u) · u^(τ_r / τ_f),   u = τ_r / (τ_r + τ_f)
t_peak_offset  = τ_r · log((τ_r + τ_f) / τ_r)
```

Peak-normalised to 1 by construction. `T_max` and `t_peak_offset` are computed **once** during `load_template_shape` and cached on the `TemplateShape` dataclass so the hot loop doesn't recompute them.

### 5.2 Injection

```python
def inject_pulse(base_samples_u16, ped_mean, ped_rms,
                 base_height, base_peak_time_ns, base_onset_t0_ns,
                 dt_ns, ratio, shape, clk_ns, noise_model, rng):
    n = base_samples_u16.shape[0]
    t_ns = np.arange(n) * clk_ns

    # ΔT is peak-to-peak (§1a.1). Because both pulses share τ_r/τ_f, the
    # onset-to-onset separation equals the peak-to-peak separation.
    inject_peak_time_ns = base_peak_time_ns + dt_ns
    inject_onset_t0_ns  = inject_peak_time_ns - shape.t_peak_offset_ns
    # Height convention (§1a.2): inject_height IS the peak height in ADC.
    inject_height       = ratio * base_height

    # Unit-peak template on the readout grid, using material-median (τ_r, τ_f).
    pulse_unit = eval_template(t_ns, inject_onset_t0_ns,
                               shape.tau_r_ns, shape.tau_f_ns, shape.T_max)

    # Compose in float.
    composite_f = (base_samples_u16.astype(np.float32)
                   + inject_height * pulse_unit.astype(np.float32))

    if noise_model == "gaussian-iid":
        # Explicit opt-in stress test — see §6.  Doubles the readout-noise
        # contribution (the base already carries real detector noise).
        noise = rng.normal(0.0, ped_rms, size=n).astype(np.float32)
        composite_f += noise
    elif noise_model == "none":
        pass
    else:
        raise ValueError(f"unknown noise_model={noise_model!r}")

    # Saturation is a bug at this point — rejected upstream (§4 cut 10).
    # Assert not clip, so a violated invariant fails loudly.
    assert composite_f.max() < MAX_ADC, (
        f"saturation invariant violated: max={composite_f.max():.1f}; "
        f"§4 cut 10 should have rejected this base")
    composite_f = np.clip(composite_f, 0.0, float(MAX_ADC - 1))  # only for the 0-floor
    composite_u16 = np.rint(composite_f).astype(np.uint16)

    truth = {
        # injected — exact truth
        "inject_height_adc":    float(inject_height),
        "inject_model_amp":     float(inject_height / shape.T_max),
        "inject_onset_t0_ns":   float(inject_onset_t0_ns),
        "inject_peak_time_ns":  float(inject_peak_time_ns),
        # base_reference — measured, not exact
        "base_height_adc":      float(base_height),
        "base_model_amp":       float(base_height / shape.T_max),
        "base_onset_t0_ns":     float(base_onset_t0_ns),
        "base_peak_time_ns":    float(base_peak_time_ns),
    }
    return composite_u16, truth
```

**t0 origin.** Both `*_peak_time_ns` and `*_onset_t0_ns` are measured relative to the readout window start (sample 0). `base_peak_time_ns` comes from `wave_ana.analyze` (`pk.time`). `base_onset_t0_ns` comes from `wave_ana.fit_pulse_shape` (the fit's own t0). The injected pulse's peak is placed at `base_peak_time_ns + dt_ns`, so ΔT is exactly peak-to-peak by construction.

**Amplitude convention (relative vs absolute).** We scan in **relative** (`ratio = h_inject / h_base`). Rationale:
- Cleaner scan: the physics plot (height error vs. ΔT vs. ratio) is naturally 2D in `(ΔT, ratio)`.
- Robust across channels: base events span 500–2000 ADC; using a ratio removes that variability from the scan axis.
- Absolute injection can trivially be recovered post-hoc from `truth_heights_adc[:, 0]`.
- If the user later wants absolute-injection mode, a `--absolute-inject-adc` flag can be added without changing the data model.

## 6. Noise handling

**Default: `--noise-model none`.**

- The base event already carries real, correlated FADC noise + real pedestal drift + all channel-specific pathologies. That's the realistic readout noise for this event.
- Adding independent Gaussian samples on top would **double-count** the electronics noise — the deployed physics analysis operates on waveforms with exactly one readout-noise contribution per sample, and the synthetic events should too.
- The injected pulse is a mathematical template with no intrinsic noise; the noise on the injected-pulse *region* of the composite is inherited from the base samples underneath it.

**Opt-in stress-test: `--noise-model gaussian-iid`.**

- Adds an independent `N(0, ped_rms)` sample per timestep on top of the composite.
- Explicitly documented in the module docstring and in `metadata["noise_model"]` as *"artificially increases the baseline RMS relative to real data; use only for benchmarking the deconvolver's noise tolerance beyond nominal."*
- The `rng` is a per-`GridPoint` `np.random.default_rng` derived from a per-cell seed (see below), so a run with the same `--seed` and `--noise-model` is bit-reproducible even in gaussian-iid mode.

**Future extension (not v0.1): `--noise-model pedestal-acf`.**

- The `inject_pulse` `noise_model` argument is a string keyed to a dispatch table; a future PR adds `"pedestal-acf"` that generates coloured noise from a measured pedestal autocorrelation, without changing the injection kernel signature or the output schema. Noted in §13.

**Seeding (both modes).** Global `--seed`, then per-cell derived seeds via `SeedSequence(args.seed).spawn(len(grid_points))`. `--seed` is recorded in metadata; each cell's spawned SeedSequence is derived deterministically from `(args.seed, cell_index)`. In `none` mode the seed is still consumed (harmlessly) to keep the CLI shape uniform.

## 7. Scanning logic — paired-base design

The grid is the Cartesian product `dt_grid × ratio_grid`. For default grids of length 10 × 6 = 60 cells × 50 events = 3000 synthetic events; comfortably in-memory (≈ 3000 × 200 × 2 B = 1.2 MB of waveforms).

**Base-event allocation — paired-base.** Each accepted `BaseEvent` is reused **once per grid cell** — one base event produces 60 synthetic events (one per `(dt, ratio)` cell), all sharing the same underlying base pulse but with different `(dt, ratio)` injections. Consequences:
- The generator needs `n_per_config` accepted bases in total (not `n_per_config × n_cells`).
- Amortises the EVIO decoding + fit-quality-gate cost.
- Gives paired statistics: at fixed `ratio`, the ΔT sweep uses the same underlying base, so any per-base variability cancels in the error-vs-ΔT plot.
- The saturation rejection (§4 cut 10) uses the largest requested ratio, so any base that survives is safe for *every* cell.

Pseudocode:

```python
def scan_grid(base_iter, grid_points, n_per_config, shape, clk_ns,
              noise_model, rng_map):
    n_bases_done = 0
    for base in base_iter:
        if n_bases_done >= n_per_config:
            break
        for gp in grid_points:
            wf, truth = inject_pulse(base.samples, base.ped_mean, base.ped_rms,
                                     base.base_height, base.base_peak_time_ns,
                                     base.base_onset_t0_ns,
                                     gp.dt_ns, gp.ratio, shape, clk_ns,
                                     noise_model, rng_map[gp.key])
            gp.append(wf, truth, base)
        n_bases_done += 1
```

Each grid cell ends up with exactly `n_per_config` events, in the same base-event order (row *i* in every file comes from the same source base). If EVIO exhausts before `n_per_config` is reached, emit a warning to stderr and write partially-filled `.npz` files (with `metadata["n_events"]` updated) — better to hand the user something than nothing.

## 7a. Injected truth vs. base reference — output-field naming

The two pulses in each synthetic event are **not epistemically equivalent**:

- **Injected pulse (column 0).** Constructed from a closed-form template with chosen `(inject_height, inject_onset_t0)`. Its truth is **exact by construction** — no measurement error. This is the primary benchmark: how well does the deconvolver recover injected height and onset?
- **Base pulse (column 1).** Comes from real EVIO. Its "truth" is really *reference observables* — `pk.height` from `WaveAnalyzer.analyze` and `t0_ns` from `WaveAnalyzer.fit_pulse_shape`. These are the deployed-analysis best estimates but they carry finite measurement error. Downstream should treat them as high-quality reference values, not ground truth.

**Naming.** Every output field uses one of two prefixes:
- `truth_*` for injected-pulse quantities (exact truth).
- `base_reference_*` for base-pulse observables when used as reference (accurate but measured).

Where a single array holds both (e.g. `truth_heights_adc` shape `(N, 2)`), the column order is fixed at `[injected, base_ref]`, documented in the `truth_source` array (§3.2) and in `metadata["column_convention"]`. The **primary** benchmark uses column 0 only.

## 8. Progress reporting

Copy the pattern from `fit_pulse_template.py:_emit_progress`. Every `--progress-every` (default 500) accepted base events, print one line to stdout:

```
[progress] file 3/12 prad_025308.evio.00002  bases_kept=1250/2400  cells_full=17/60  rate=42 base/s  elapsed=29.8s
```

- `bases_kept` = accepted base events / **`n_per_config`** (the paired-base target — one accepted base fills a row in every cell).
- `cells_full` = grid cells that have reached `n_per_config`. In paired-base mode, this is 0/N until `n_bases_done == n_per_config`, then jumps to N/N in a single tick — that's expected. `bases_kept` is the meaningful progress bar.
- `rate` = accepted-base rate, not synthetic-event rate.
- `elapsed` = wall clock since `main()` started.

Additional log lines (once each):
- `[setup]` — template path, material, `(τ_r, τ_f, T_max, t_peak_offset)`, clk_ns, grid dims, target = `n_per_config`.
- `[file N/M] <path>` for each EVIO split opened, mirroring `fit_pulse_template.py`.
- `[write] <path> (N events) sidecar=<path>.json` for each `.npz` written.
- `[done] target=X kept=Y wrote=Z files elapsed=…` followed by the rejection-reason histogram.

## 9. Error handling

| Condition | Behaviour |
|-----------|-----------|
| `--template` JSON missing | `SystemExit(1)` with the path in the message. |
| Template JSON has no `_by_type[--material]` entry | `SystemExit(1)` listing the materials that *are* present in `_by_type`. |
| Template entry has NaN/None for `tau_r_ns.median` or `tau_f_ns.median` | `SystemExit(1)` — the shape is unusable. |
| `--dt-grid` or `--ratio-grid` parse failure | `argparse.ArgumentTypeError` from the parser callable. |
| Any `dt_grid` value ≤ 0 or `ratio_grid` value ≤ 0 | rejected in the parser callable. |
| No EVIO splits discovered across all `--evio_paths` | `SystemExit` with the input list (mirrors `fit_pulse_template.py`). |
| EVIO input has zero usable base events | write no `.npz`, write an empty manifest with `"files": []`, exit code 2 with a stderr warning. |
| Some grid cells fill, others don't | write everything that filled; write partial cells with fewer events + updated `metadata["n_events"]`; warn on stderr listing under-filled cells. |
| `open_auto(evio_path) != Status.success` on some split | warn to stderr, skip that split, continue (mirrors `fit_pulse_template.py`). |
| `KeyboardInterrupt` | write out whatever grid cells have data, print `[interrupted — partial output]`, exit non-zero. |

## 10. Validation / sanity checks

**Two modes.**

```python
def validate_output(npz_path: Path, wave_ana: "dec.WaveAnalyzer | None" = None) -> None:
    """Schema-only if wave_ana is None; adds peak-detect smoke test if provided."""
```

**Schema-only checks (always run, no C++ deps):**

1. `waveforms.dtype == np.uint16` and `waveforms.shape == (n_events, n_samples)`.
2. `truth_heights_adc.shape == (n_events, 2)`, all finite, all `> 0`.
3. `truth_onset_t0_ns.shape == (n_events, 2)`; `truth_peak_times_ns.shape == (n_events, 2)`.
4. `truth_peak_times_ns[:, 0] - truth_peak_times_ns[:, 1]` all within `1e-6` of the sidecar's `dt_ns` (round-trip on ΔT, peak-to-peak).
5. `truth_heights_adc[:, 0] / truth_heights_adc[:, 1]` all within `1e-6` of the sidecar's `amplitude_ratio` (round-trip on ratio).
6. `truth_peak_times_ns - truth_onset_t0_ns` all within `1e-6` of `meta_t_peak_offset_ns` (round-trip on the peak-vs-onset relationship).
7. `truth_model_amplitudes * meta_T_max` matches `truth_heights_adc` within `1e-4` relative (round-trip on amplitude convention).
8. **Saturation invariant.** `waveforms.max() < MAX_ADC - safety_margin` (default `4095 - 50 = 4045`). If any waveform is at or near `MAX_ADC`, an `AssertionError` fires with the offending row.
9. `np.load(npz_path, allow_pickle=False)` succeeds — the schema is pickle-free.
10. Sidecar JSON exists (or the `metadata_json` string scalar is present); JSON parses; declared `n_events` matches array shape.

**With `wave_ana` (opt-in, called at end of `main`):**

11. For a random 5% subsample of the file's events: rerun `wave_ana.analyze(waveforms[i])` and assert *at least 1 peak is detected*. Purposely lax: at small ΔT the soft analyzer *should* merge the two peaks into one `Q_PEAK_PILED` — exactly the regime the downstream deconv targets.

Validation failures are `AssertionError` with a descriptive message; they should never trigger in normal operation.

**`--self-test`** uses **schema-only mode** (no `wave_ana` needed): it builds a synthetic base waveform from the same template, runs the full inject → serialise → reload → schema-validate cycle on a 2×2 grid in a `tempfile.TemporaryDirectory`, and exits 0 if all checks pass. This means CI can run `--self-test` without needing DAQ config, HyCal map, or EVIO files.

## 11. CLI

Following `fit_pulse_template.py`'s conventions (including multi-input EVIO):

```
python pileup_generator.py <evio_path> [<evio_path> ...] \
    --template PATH             # required
    --out-dir DIR               # required
    [--material PbWO4]          # default PbWO4
    [--n-per-config 50]
    [--dt-grid 8,12,16,20,24,32,48,64,96,128]   # comma-sep ns
    [--ratio-grid 0.1,0.2,0.3,0.5,1.0,2.0]      # comma-sep
    [--height-min 500] [--height-max 2000]      # §2b physics-population defaults
    [--height-rms-mult 10.0]
    [--t0-min 25.0]                             # §2b Population A cut
    [--chi2-max 5.0]                            # fit-quality gate
    [--residual-veto-sigma 5.0]                 # hidden-pileup veto
    [--post-pad 40]             # samples of headroom past the last injected peak
    [--noise-model none]        # {none, gaussian-iid}; default: none
    [--seed 12345]
    [--max-events 0]            # cap on physics events scanned (0 = all)
    [--progress-every 500]
    [--daq-config PATH] [--hc-map-file PATH]
    [--force]                   # allow overwriting non-empty --out-dir
    [--metadata-in-npz]         # embed metadata JSON string inside .npz (no sidecar)
    [--self-test]               # schema-only in-memory smoke test; exit 0/nonzero
```

The two grid args use a custom argparse type:

```python
def _parse_float_grid(s: str) -> list[float]:
    vals = [float(x) for x in s.split(",") if x.strip()]
    if not vals or any(v <= 0 for v in vals):
        raise argparse.ArgumentTypeError(
            f"grid must be a non-empty comma-separated list of positive floats; got {s!r}")
    return sorted(set(vals))
```

Deduplication + sorting is intentional.

## 12. Rough size estimate

| Block | LOC |
|-------|-----|
| Module docstring (usage, §1a conventions, output schema, non-goals) | 90 |
| Imports + constants | 30 |
| Dataclasses (`BaseEvent`, `GridPoint`, `TemplateShape`) | 55 |
| `load_template_shape` + T_max / t_peak_offset precompute + diagnostics | 55 |
| `iterate_base_events` (EVIO loop + 10 cuts including fit-quality gate + residual veto) | 180 |
| `eval_template` + `inject_pulse` (with noise-model dispatch) | 80 |
| `scan_grid` + `GridPoint.append` | 45 |
| `write_grid_point` (.npz + sidecar .json) + manifest writer | 75 |
| `validate_output` (schema-only + optional wave_ana smoke test) + `--self-test` | 100 |
| Progress + logging | 30 |
| `main()` (argparse + orchestration + multi-input EVIO discovery) | 110 |
| **Total** | **~850** |

Slightly larger than iteration 1's estimate because of the fit-quality gate, residual-veto pass, and the split truth conventions. Still comfortably under `fit_pulse_template.py` (1080 LOC).

## 13. Explicit non-goals for v0.1

Called out in the module docstring and here so the review focuses on what v0.1 *does* cover:

1. **Coloured / autocorrelated noise (pedestal-ACF).** Roadmap item; the `inject_pulse` `noise_model` argument dispatches by string, so a future PR adds `"pedestal-acf"` without touching the injection kernel or output schema.
2. **Cluster-coincident correlated pile-up** (multiple channels injected together with a shared shower ancestor). v0.1 is a single-channel generator; a `cluster_pileup_generator.py` companion is the natural follow-up.
3. **Poisson-random pile-up timing.** Roadmap item; v0.1 produces a deterministic `(ΔT, ratio)` response surface, not Poisson-distributed pile-up. A separate driver script for the "Poisson-random pile-up ROC vs luminosity" plot is out of scope here.
4. **Luminosity-scaling / rate-quote table.** Belongs in the ROC analysis / write-up, not the generator.
5. **More than 2 pulses per event.** v0.1 injects exactly 1 pulse on top of the base.
6. **Per-channel `(τ_r, τ_f)` overrides.** v0.1 always uses the `_by_type[--material]` median.
7. **Amplitude scan in absolute ADC.** Ratio is the v0.1 scan axis; `--absolute-inject-adc` is a possible extension.
8. **Any C++ modification.** Everything the generator needs is already exposed by `prad2py`.
9. **Any deconvolution.** The generator produces inputs for `wave_ana.deconvolve()`; running it is downstream.
10. **Plotting.** The generator writes `.npz` + sidecar `.json` + `manifest.json` only.

**Bottom line on roadmap scope.** v0.1 satisfies the deterministic-response-surface *piece* of roadmap item 2, and produces the "amplitude error vs ΔT vs ratio" plot inputs. The Poisson-random-timing, cluster-coincident, and pedestal-ACF pieces of that same roadmap item are deferred to follow-up deliverables.

## 14. Design decisions to confirm before implementation

Everything below is a defaulted choice where genuine alternatives exist — the plan proceeds with the recommendation, but each is worth explicit sign-off:

| # | Decision | Recommendation | Alternative |
|---|----------|----------------|-------------|
| A | Amplitude convention | **Peak height (`truth_heights_adc`) as primary, model amplitude as convenience** — matches `dec_out.height` for the primary comparison | Only store model amplitude; force users to compute heights themselves |
| B | Base-event reuse | **Paired-base: one base reused across all 60 cells** (§7) — paired stats, EVIO-cost amortised | Draw fresh base per cell (uncorrelated stats, ~60× more EVIO reads, breaks paired-difference plots) |
| C | Default noise model | **`none`** (§6) — base already carries real readout noise; no double-counting | `gaussian-iid` — artificially inflates baseline RMS |
| D | Output layout | **One `.npz` per grid cell + sidecar `.json` per cell + top-level `manifest.json`** (§3.2, §3.4) | One big `.npz` with all events + a grid-index column; `--metadata-in-npz` toggle for a single-file variant |
| E | Base cut envelope | **`--height-min 500`, `--t0-min 25.0`, `--chi2-max 5.0`, residual-veto `5σ`** (§4) — matches `fit_pulse_template.md` §2b physics selection | Looser cuts (`--height-min 300`, no fit-quality gate) — admits more bases but contaminates truth with off-shape / hidden-pileup events |

Items removed from this table since iteration 1 (they don't need sign-off — the plan just does the right thing):
- Saturation guard on/off: guard is **always on**; clipping breaks truth, so there is no "off" that preserves the deliverable.
- Global `--seed` for reproducibility: no reasonable alternative — non-reproducibility is a bug.
- Overwrite policy: `--force` guard is standard defensive practice; not a design choice.
- Filename ratio encoding (`× 1000` zero-padded): a naming convention, not a design decision — sortability is a hard requirement.
- `allow_pickle=True`: no longer relevant — the sidecar-JSON format eliminates the need for pickle entirely.

## 15. Downstream consumption pattern (what the plan is designed to make easy)

The downstream user's benchmarking loop should read:

```python
import json
import numpy as np
from pathlib import Path
from prad2py import dec

# Load one grid cell — no pickle needed.
npz_path = Path("output/pileup_synth/pileup_PbWO4_dt016_ratio0500.npz")
data = np.load(npz_path, allow_pickle=False)
waveforms          = data["waveforms"]                # (N, n_samples), uint16
truth_heights      = data["truth_heights_adc"]        # (N, 2)  col 0 = injected, col 1 = base_ref
truth_onset_t0     = data["truth_onset_t0_ns"]        # (N, 2)  compare to dec_out.t0_ns
truth_peak_times   = data["truth_peak_times_ns"]      # (N, 2)  compare to pk.time from analyze
truth_model_amps   = data["truth_model_amplitudes"]   # (N, 2)  compare to dec_out.amplitude

# Sidecar metadata — plain JSON, jq-inspectable.
meta = json.loads(npz_path.with_suffix(".json").read_text())
dt_ns   = meta["dt_ns"]
ratio   = meta["amplitude_ratio"]
T_max   = meta["T_max"]

# Set up analyzer + template store once (identical to deconv_pileup_demo.py).
cfg   = dec.load_daq_config()
wcfg  = dec.WaveConfig(cfg.wave_cfg)
store = dec.PulseTemplateStore()
store.load_from_file(tmpl_path, wcfg)
wa    = dec.WaveAnalyzer(wcfg)

# v0.1 is PbWO4-only, so a single per-type template covers every event in the file.
tmpl = store.type_template(meta["material"])          # correct API — NOT lookup_by_type
assert tmpl is not None, f"no template loaded for material {meta['material']!r}"

# Per-event: analyze → deconvolve → compare to injected truth (column 0).
for i in range(waveforms.shape[0]):
    samples = waveforms[i]
    wres    = wa.analyze_result(samples)
    out     = wa.deconvolve(samples, wres, tmpl)
    # Correct comparisons — see §1a:
    #   dec_out.height   ↔  truth_heights_adc[i, 0]     (injected peak height in ADC)
    #   dec_out.t0_ns    ↔  truth_onset_t0_ns[i, 0]      (fitted template onset)
    #   dec_out.amplitude ↔ truth_model_amps[i, 0]       (unnormalized model coeff)
    #
    # Do NOT compare dec_out.t0_ns to truth_peak_times_ns — those differ by
    # t_peak_offset (~5.5 ns for PbWO4).
```

Three "make it easy" points the plan takes care of:

1. **Correct template lookup API.** `store.type_template(meta["material"])` — the real per-type accessor in `prad2dec/include/PulseTemplateStore.h`. (Iteration-1 draft used `store.lookup_by_type(...)`, which does not exist.)
2. **Correct amplitude comparison.** `dec_out.height` vs. `truth_heights_adc[:, 0]`. `dec_out.amplitude` (the unnormalized model coefficient) vs. `truth_model_amplitudes[:, 0]`. Mixing these differs by a factor of `T_max` and would silently look like a bias.
3. **Correct time comparison.** `dec_out.t0_ns` vs. `truth_onset_t0_ns[:, 0]`. Comparing against `truth_peak_times_ns[:, 0]` would produce a spurious ~5.5 ns timing bias.

Plus: fixed shape per file (every waveform has the same `n_samples`; every truth row is `(2,)`), no ragged arrays, no pickle.

## Todo

### Phase 1 — Skeleton and template loading

- [ ] Create `analysis/pyscripts/pileup_generator.py` with module docstring quoting §1, §1a (both conventions), §5.1 formula, §3 data model, and the non-goals from §13.
- [ ] Import `_common as C`, `_common.dec`, `numpy`, `argparse`, `json`, `time`, `pathlib`, `dataclasses`, `tempfile`.
- [ ] Define constants: `DEFAULT_DT_GRID_NS`, `DEFAULT_RATIO_GRID`, `DEFAULT_MATERIAL = "PbWO4"`, `MAX_ADC = 4095`, `SATURATION_SAFETY = 50`, `DEFAULT_POST_PAD = 40`, `DEFAULT_HEIGHT_MIN = 500`, `DEFAULT_HEIGHT_MAX = 2000`, `DEFAULT_T0_MIN_NS = 25.0`, `DEFAULT_CHI2_MAX = 5.0`, `DEFAULT_RESIDUAL_VETO_SIGMA = 5.0`.
- [ ] Define dataclasses `TemplateShape` (with cached `T_max` and `t_peak_offset_ns`), `BaseEvent` (with fit-quality diagnostics + full provenance per §3.1), `GridPoint`.
- [ ] Implement `load_template_shape(json_path, material) -> TemplateShape`, computing `T_max = (1-u)·u^(τ_r/τ_f)` and `t_peak_offset_ns = τ_r · log((τ_r+τ_f)/τ_r)`. Error paths: missing file, missing material, NaN medians.
- [ ] Wire a template-load self-check into `--self-test`: assert PbWO4 medians ≈ 2.24 / 24.0 within ±0.5 ns; assert `t_peak_offset ≈ 5.5 ns` and `T_max ≈ 0.727` within tolerance.

### Phase 2 — Base-event iteration (10 cuts including fit-quality + residual veto)

- [ ] Implement multi-input EVIO discovery mirroring `fit_pulse_template.py` lines 661–684 **exactly**: collect `args.evio_paths`, run `C.discover_split_files`, dedupe order-preserving, `SystemExit` on empty, call `C.setup_pipeline(evio_path=args.evio_paths[0], ...)`, then override `p.evio_files = all_files`.
- [ ] Implement `iterate_base_events(pipeline, args, shape) -> Generator[BaseEvent]` copying the EVIO/`ch.read`/`ch.scan`/`ch.decode_event` loop shape from `fit_pulse_template.py` (~lines 731–867).
- [ ] Apply the 10 cuts from §4 in order. Each cut increments a named counter (`n_reject_multi_peak`, `n_reject_quality`, `n_reject_overflow`, `n_reject_height`, `n_reject_snr`, `n_reject_material`, `n_reject_fit_converged`, `n_reject_fit_t0`, `n_reject_fit_chi2`, `n_reject_residual_veto`, `n_reject_headroom`, `n_reject_saturation`) so `[done]` prints the full histogram.
- [ ] Call `wave_ana.fit_pulse_shape(...)` for cut 7 (fit-quality gate). Store `chi2_per_dof`, `converged`, fitted `t0_ns` on the `BaseEvent`.
- [ ] Implement the residual hidden-pileup veto for cut 8: reconstruct fitted single-pulse waveform, subtract from samples, scan outside `±3·τ_f` around the main peak for residual samples exceeding `residual_veto_sigma * ped_rms`.
- [ ] Compute the injection-headroom cut against `max(dt_grid) / clk_ns + post_pad`; cache `dt_max_samples` once.
- [ ] Materialise `BaseEvent.samples` as a copy so the underlying C++-owned buffer is not aliased into the accumulator.
- [ ] Record full provenance on each accepted base: `roc_tag`, `slot`, `chan`, `channel_id`, `channel_name`, `source_file` (basename), `physics_event_index`, `run_number`.
- [ ] Test on a small slice of run 025308 (`--max-events 500`) and print rejection-reason histogram.

### Phase 3 — Injection kernel + noise-model dispatch

- [ ] Implement `eval_template(t_ns, t0_ns, tau_r, tau_f, T_max)` per §5.1. Vectorised; must handle `t0 > t_ns.max()` (return zeros).
- [ ] Implement `inject_pulse(base_samples_u16, ped_mean, ped_rms, base_height, base_peak_time_ns, base_onset_t0_ns, dt_ns, ratio, shape, clk_ns, noise_model, rng)` per §5.2 — returns `(uint16 waveform, truth_dict)` with all 8 truth keys.
- [ ] `noise_model` dispatch: `"none"` → no extra noise; `"gaussian-iid"` → add `N(0, ped_rms)` iid; unknown → `ValueError`.
- [ ] Add inline asserts: shape of output == shape of input, dtype == `uint16`, `truth["inject_peak_time_ns"] == base_peak_time_ns + dt_ns` exactly, `truth["inject_onset_t0_ns"] == truth["inject_peak_time_ns"] - shape.t_peak_offset_ns` exactly, saturation invariant (`composite_f.max() < MAX_ADC`).
- [ ] Unit test in `--self-test`: build a synthetic base (unit-peak template scaled to `h=1000` at `peak_time=100 ns`, on a 200-sample uint16 grid), inject at `ΔT=16, ratio=0.5`, verify the composite has a local maximum near sample `100/4 + 16/4 = 29` past the base peak with height `≈ 500 ADC`. Repeat under `--noise-model gaussian-iid` to smoke-test that branch.

### Phase 4 — Scanning driver + RNG plumbing

- [ ] Implement `GridPoint` with `.append(waveform, truth, base)` that grows pre-allocated numpy buffers (`n_per_config × n_samples` upfront; separate buffers for each truth column and each provenance column) to avoid list-append re-copies.
- [ ] Implement `scan_grid(base_iter, grid_points, n_per_config, shape, clk_ns, noise_model, rng_map)` per §7 — paired-base loop, stop when `n_bases_done == n_per_config`.
- [ ] Build `rng_map` using `np.random.SeedSequence(args.seed).spawn(len(grid_points))` — each cell gets an independent, reproducible stream. Consumed even in `--noise-model none` mode.
- [ ] Add a `full` predicate on `GridPoint` for consistency with the progress reporter (`cells_full`).

### Phase 5 — Serialisation (npz + sidecar JSON) + manifest

- [ ] Implement `write_grid_point(out_dir, gp, template_shape, meta_common)` per §3.2/§3.3. Use `np.savez` (not `_compressed`) and write **all** arrays listed in §3.2 (waveforms, both truth conventions, both time conventions, `truth_source`, pedestals, full base-provenance columns, `meta_*` scalars).
- [ ] Write sidecar `<basename>.json` per grid cell with the schema in §3.3 (including `column_convention`).
- [ ] If `--metadata-in-npz` is set, also embed the JSON string under `metadata_json` key inside the `.npz` (no pickle — pure fixed-width string scalar).
- [ ] Implement `write_manifest(out_dir, meta_common, grid_points, cli_argv)` per §3.4, listing every produced file with its sidecar path and `n_written`.
- [ ] Implement the `--force` guard on non-empty `--out-dir`.

### Phase 6 — Progress + logging + error handling

- [ ] Copy `_emit_progress` shape from `fit_pulse_template.py`; adapt fields per §8. In the paired-base design, `bases_kept` is the progress bar; `cells_full` jumps 0 → N in a single tick and that's expected — document this in the log-line comment.
- [ ] Wrap the main scan loop in `try / except KeyboardInterrupt` for partial-output preservation (§9).
- [ ] Add stderr warnings + non-zero exit code (2) for the "no usable base events" case.
- [ ] Add stderr warnings for under-filled cells at end of run.
- [ ] Print rejection-reason histogram in the `[done]` line.

### Phase 7 — Validation + self-test

- [ ] Implement `validate_output(npz_path, wave_ana=None)` with the 10 schema-only assertions from §10 plus the optional `wave_ana` smoke test (§10 item 11). Enforce `np.load(..., allow_pickle=False)`.
- [ ] Implement `--self-test` (§10 last paragraph): build a synthetic in-memory base, run the full inject/serialise/reload/schema-validate cycle on a 2×2 grid in a `tempfile.TemporaryDirectory`, print `[self-test] OK`, exit 0. **No `wave_ana` in self-test** — schema-only, so no DAQ config / HyCal map / EVIO needed.
- [ ] Verify `--self-test` runs standalone: `python pileup_generator.py --self-test` completes with exit 0.

### Phase 8 — CLI + integration test

- [ ] Wire `argparse` per §11 including the `_parse_float_grid` type callable, `--noise-model {none,gaussian-iid}`, `--t0-min`, `--chi2-max`, `--residual-veto-sigma`, `--metadata-in-npz`, and multi-input EVIO (`nargs='+'`).
- [ ] Run end-to-end on a small slice: `python pileup_generator.py <run025308_evio> --template output/pulse_templates_025308_h500_t0cut.json --out-dir output/pileup_synth_test --n-per-config 5 --dt-grid 16,32 --ratio-grid 0.5,1.0 --max-events 5000`. Confirm 4 `.npz` + 4 sidecar `.json` + 1 `manifest.json` written; `validate_output` (schema-only) passes on all.
- [ ] Round-trip test with `wave_ana`: `np.load` one output file, run the smoke path in `validate_output(path, wave_ana=wa)`, confirm at least one peak is detected in the 5% subsample.
- [ ] Round-trip test for the `--noise-model gaussian-iid` branch on the same small slice; confirm `metadata["noise_model"] == "gaussian-iid"`.

### Phase 9 — Documentation cross-links

- [ ] Add a two-paragraph note to `docs/technical_notes/waveform_analysis/wave_analysis.md` (or a new short doc) pointing at `pileup_generator.py` as the truth source for the pile-up benchmark and calling out the §1a conventions (peak vs onset, height vs model amplitude). **This documentation task is handed off to `doc-writer`; not in-scope for the script implementer.**
- [ ] Mention `pileup_generator.py` in `analysis/pyscripts/`'s pattern-of-usage reference if such a README exists (check during Phase 8).

### Non-goals reminder (do NOT do in this deliverable)

- Do NOT implement the ROC / benchmark analysis script.
- Do NOT implement Poisson-random pile-up timing (roadmap follow-up).
- Do NOT implement cluster-coincident cross-channel pile-up (roadmap follow-up).
- Do NOT implement pedestal-ACF noise (roadmap follow-up; hook exists via `noise_model` dispatch).
- Do NOT modify any C++ code.
- Do NOT commit or push.
- Do NOT add plotting.
