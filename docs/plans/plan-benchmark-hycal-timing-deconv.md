# Plan — `benchmark_hycal_timing_deconv.py`

**Status:** draft (iteration 2, pre-approval — planning-critic blocking
fixes applied: deconv success-state set, HyCalClusterConfig pattern,
`--require-trigger` ownership, `pk.integral` fallback)
**Scope:** HyCal timing-gated clustering benchmark that uses the C++
template-fit path (`WaveAnalyzer.deconvolve`) for pulse amplitudes and
the standard C++ CFD time (`WaveResult.peaks[k].time`) for timing,
compared across a 2×2 grid of {legacy, gated} × {largest-E,
smallest-|Δt|} neighbour-selection policies on the same event stream.
Python-only. No C++ changes. No pybind11 changes.
**Hard constraint honoured:** the literal matched-filter → template-fit
→ CFD pipeline is NOT in Phase 1 — see §3. Phase 1 uses the pipeline
that is already callable from Python today (§9).

## Research inputs used

- `analysis/pyscripts/study_hycal_timing.py` lines 150–210 (per-event
  decode loop, `is_hycal()` channel filter, `add_hit` + `collect_neighbor_timing`
  pattern, TSV 14-column schema).
- `analysis/pyscripts/benchmark_hycal_timing.py` lines 1–60, 370–420
  (legacy vs. gated dual-cluster pattern, `HAVE_SCIPY` idiom, Gaussian
  fit + CB extension, calibration workflow, `--signal-min` / `--e-range`
  / `--bin-width` / `--calibrate` CLI).
- `python/bind_dec.cpp` lines 240–295 (`PulseTemplateStore.type_template`),
  361–397 (`DeconvOutput` fields: `state, n, chi2_per_dof, amplitude[],
  height[], integral[], t0_ns[], tau_r_ns[], tau_f_ns[]`), 542–582
  (`WaveAnalyzer.analyze_result` + `deconvolve` signatures), 612–617
  (`Q_DECONV_*` constants exposed at module scope).
- `python/bind_det.cpp` lines 486–491 (`ModuleType` enum: `PbGlass`,
  `PbWO4`, `LMS`, `Veto`, `Unknown`), 529–572 (`Module` bindings:
  `name`, `id`, `index`, `type`, `is_pwo4()`, `is_glass()`, `is_hycal()`,
  `energize()`).
- `prad2det/include/HyCalSystem.h` lines 27, 170–174 (`enum class
  ModuleType : int8_t { PbGlass=0, PbWO4=1, LMS=2, Veto=3, Unknown=-1 }`;
  `is_hycal()` returns `is_pwo4() || is_glass()` — i.e. **PbWO4 and
  PbGlass only; LMS and Veto excluded**).
- `prad2det/src/HyCalCluster.cpp` lines 197–209 (BFS neighbour pulse
  selection: always largest-energy, no runtime switch), 784–794
  (`CollectNeighborTiming` emits one row per (seed-pulse, neighbor-pulse)
  pair — smallest-|Δt| grouping is pure Python-side post-processing).
- `prad2dec/src/WaveAnalyzer.cpp` lines 518–562 (CFD baked into private
  `findPeaks`, hardcoded `cfd_fraction = 0.5f`, result stored on
  `peak.time`), 832–835 (explicit `wa.deconvolve(samples, wres, tmpl)`
  bypasses the `nnls_deconv.enabled` production gate).
- `analysis/pyscripts/benchmark_matched_filter.py` lines 205, 246, 276
  (matched-filter prototype — Python-only, `np.correlate`-based, **not
  in Phase 1**; see §3 non-goals and §16 Phase 2 handoff).
- `analysis/pyscripts/fit_pulse_template.py` line 946
  (`mtype = mod.type.name` — canonical enum→string conversion for
  `type_template(...)` lookup).
- `analysis/pyscripts/_common.py` (CLI common args, `setup_pipeline`,
  `open_table_writer`).
- `database/daq_config.json` line 117 (`nnls_deconv.enabled: false` in
  production — motivates Phase 1).
- Template file `analysis/pyscripts/output/pulse_templates_025308_h50_RawSum_max20modules.json`
  (schema matches `PulseTemplateStore::LoadFromFile`; `_by_type` has
  PbGlass, PbWO4, Veto — no LMS; not an issue because `is_hycal()` has
  already excluded LMS upstream).

---

## 1. Purpose

Combine the per-event FADC-decode + HyCal-cluster-hit flow of
`study_hycal_timing.py` with the parallel-clusterer comparison and
σ_E/E Gaussian/Crystal-Ball fits of `benchmark_hycal_timing.py`, but
replace the raw `WaveAnalyzer.analyze` peak extraction with the
template-fit path: `wa.analyze_result(samples)` for CFD-timed
candidate detection, then `wa.deconvolve(samples, wres, tmpl)` for
per-pulse LM template-fit amplitudes. Compare four clustering
configurations on the same event stream:

| config | `seed_time_window` | neighbour-selection policy |
|---|---|---|
| **legacy, largest-E**    | ≤ 0 (disabled) | one largest-integral peak per module |
| **legacy, smallest-Δt**  | ≤ 0 (disabled) | one closest-to-seed-time peak per module |
| **gated, largest-E**     | `--window` ns (default 4) | all in-window peaks; C++ BFS picks largest-energy |
| **gated, smallest-Δt**   | `--window` ns (default 4) | Python-side pre-filter to one closest-to-seed-time peak per module, then C++ BFS runs its default largest-E logic on a trivial 1-pulse-per-module input |

The headline question Phase 1 answers: does using template-fit
amplitudes + CFD times + smallest-Δt neighbour selection restore the
expected `σ_E/E(gated) ≤ σ_E/E(legacy)` ordering, which is currently
inverted in `benchmark_hycal_timing.py`'s baseline on run 25320?

## 2. Inputs

| input                        | source                                                                                                      | note                                                                                            |
|------------------------------|-------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------|
| EVIO files                   | positional `EVIO_PATH` (file or glob) via `_common.setup_pipeline`                                           | same convention as `study_hycal_timing.py` and `benchmark_hycal_timing.py`.                     |
| HyCal map                    | `--hc-map-file`                                                                                              | passed through to `_common.setup_pipeline`.                                                     |
| HyCal calibration            | `--hc-calib-file`                                                                                            | `calibration_factor_2p2_fadc.json` for 2.2 GeV e-p; required for `mod.energize()` to return MeV.|
| DAQ config                   | `--daq-config`                                                                                               | `daq_config.json`; drives `WaveConfig` + `clk_mhz`.                                             |
| Pulse-template store         | `--template-file`                                                                                            | default `analysis/pyscripts/output/pulse_templates_025308_h50_RawSum_max20modules.json`.         |
| Trigger filter               | `--require-trigger`                                                                                          | default 0 (physics unbiased, same as `study_hycal_timing.py`).                                  |

The template file provides per-type templates for PbGlass, PbWO4,
Veto. LMS is absent from `_by_type` in this file — not an issue
because the channel filter (§12) drops LMS before any pulse
extraction.

## 3. Non-goals / explicitly out of scope

- **Matched filter in the detection loop.** The Python prototype in
  `benchmark_matched_filter.py` is NOT wired in. Candidate peaks come
  from `wa.analyze_result` (which internally runs the standard
  `findPeaks` local-maxima + prominence pipeline). Enabling the
  literal matched-filter → template-fit → CFD pipeline is Phase 2 —
  see `plan_waveform.md` §5. Reason: `WaveResult.npeaks` and
  `.peaks` are read-only in the pybind11 binding
  (`python/bind_dec.cpp:366`), so there is no supported path to
  construct a `WaveResult` from MF candidates and feed it to
  `deconvolve`. Making this work needs a binding change.
- **Isolated CFD call.** CFD timing is baked inside the private
  `WaveAnalyzer::findPeaks` at `prad2dec/src/WaveAnalyzer.cpp:518–562`
  with a hardcoded `cfd_fraction = 0.5f`. It is not callable in
  isolation from Python. We consume it indirectly via
  `wres.peaks[k].time`. Exposing a public `compute_cfd_time(...)`
  helper is Phase 2.
- **C++ changes.** No edits to `WaveAnalyzer`, `HyCalCluster`,
  `PulseTemplateStore`, or any header. No pybind11 changes.
- **Changing the C++ cluster BFS pulse-selection policy to
  smallest-Δt.** Would require editing
  `HyCalCluster::grow_island` at
  `prad2det/src/HyCalCluster.cpp:197–209`. Phase 3 — deferred until
  Phase 1 data shows this is unambiguously better.
- **Explicit V1–V4 veto-counter handling.** Phase 1 processes only
  HyCal physics calorimeter modules (PbWO4 + PbGlass) via
  `mod.is_hycal()`, which also incidentally excludes LMS and Veto
  modules (`HyCalSystem.h:174`). No additional V1–V4 veto-specific
  logic is implemented. This matches the convention used by sibling
  scripts `study_hycal_timing.py:168` and
  `benchmark_hycal_timing.py:379`. If a future downstream
  track-matching stage requires V-counter inhibit windows, that is a
  separate script — not Phase 1.
- **Writing or editing any production `daq_config.json` /
  `nnls_deconv.enabled` flag.** Phase 1 is a measurement, not a
  configuration change. Any recommendation to flip
  `nnls_deconv.enabled = true` is documented in `plan_waveform.md`
  Phase 3 and does not happen inside this script.

## 4. File location and naming

**Target path:** `analysis/pyscripts/benchmark_hycal_timing_deconv.py`

Single file. Consistent with sibling scripts
(`benchmark_hycal_timing.py`, `study_hycal_timing.py`,
`benchmark_matched_filter.py`, `benchmark_deconv.py`).

Output stem is `<out_path>` (positional), with the four artifact files
suffixed `_dt.tsv`, `_spectrum.png`, `_summary.csv`,
`_calibrated.png` (last only with `--calibrate`). Same stem-suffix
convention as `benchmark_hycal_timing.py`.

## 5. Dependencies

Hard: `numpy`, `matplotlib`, `prad2py` (via the `dec`, `det` modules).
Soft: `scipy.optimize.curve_fit` for the Gaussian / Crystal-Ball
fits, via the standard `HAVE_SCIPY` try/except idiom at the top of
`benchmark_hycal_timing.py:44–48`. If scipy is absent the script
still produces the TSV + PNG + summary CSV, but with NaN fit
parameters and a stderr warning (same policy as
`benchmark_hycal_timing.py` and `module_dt_diagnostic.py`).

No change to any requirements file. The module docstring lists the
deps explicitly.

## 6. CLI

Positional + common args come from `_common.add_common_args`
(`_common.py:364–392`): `EVIO_PATH`, `OUT_PATH`, `--hc-map-file`,
`--hc-calib-file`, `--daq-config`, `--run-num`, `--max-events`,
`--csv`, `--no-header`. `--require-trigger` is NOT in
`add_common_args` and must be declared by this script (see the
`event-level cuts` block below).

Script-specific flags:

```
input / pipeline:
  --template-file PATH          Per-type pulse-template JSON consumed by
                                PulseTemplateStore.load_from_file.
                                Default: analysis/pyscripts/output/
                                  pulse_templates_025308_h50_RawSum_max20modules.json

timing-gate configs (all four always run; this flag sets the gated-config window):
  --window W                    seed_time_window in ns for the two gated configs
                                (default 4.0; same as benchmark_hycal_timing.py).

pulse pre-filter (same semantics as study_hycal_timing.py):
  --pre-window LO HI            keep peaks with LO < pk.time < HI in ns
                                (default 100 200).

event-level cuts:
  --require-trigger MASK        Trigger-bits gate declared by this script
                                (NOT via _common.add_common_args).  Accepts
                                hex (0x...) or decimal.  Default 0 — accept
                                any physics trigger (unbiased sample).
                                Rationale: study_hycal_timing.py defaults
                                to C.PHYSICS_TRIGGER_BITS (SSP0 strict
                                equality) for its diagnostic TSV, but for
                                a σ_E/E benchmark we want the broadest
                                physics sample to minimise statistical
                                uncertainty on the fitted width.  Pass an
                                explicit mask (e.g. 0x1 for the matching
                                SSP0 bit) to reproduce the strict-equality
                                behaviour.
  --max-fired-modules N         skip events where more than N HyCal channels
                                produced ≥ 1 peak passing --pre-window
                                (default 50; backstop against LMS flashes that
                                illuminate physics channels + beam-splash).

CollectNeighborTiming:
  --max-qdist Q                 max quantized (dx_q,dy_q) distance (default 5.0;
                                passed to collect_neighbor_timing).

spectrum selection:
  --signal-min E_MeV            one-cluster cut (default 500 MeV; tuned for 2.2 GeV
                                e-p vs benchmark_hycal_timing.py's 3000 for 3.5 GeV).
  --e-range LO HI               energy histogram range in MeV (default 500 2500;
                                also tuned for 2.2 GeV e-p).
  --bin-width BW                histogram bin width in MeV (default 20).

calibration (identical semantics to benchmark_hycal_timing.py):
  --calibrate                   per-seed-module gain refinement + CB re-fit.
  --target-peak E_MeV           calibration target peak position (default 2239.0;
                                e-p elastic at 2.2 GeV).
  --inner-radius R_mm           radius cut for per-module gain fit (default 200.0).
  --min-seed-energy E_MeV       seed-energy cut for per-module gain fit (default 200.0).
  --min-per-module N            minimum events per module for per-module gain fit
                                (default 30).
  --exclude-rowcol R_LO R_HI C_LO C_HI
                                optional exclude region (default None; mirrors
                                benchmark_hycal_timing.py's exclude list).
  --cal-method                  "median" or "cb-peak" (default "median").
```

All defaults overridable; `_validate_args` enforces (`e_range[0] <
e_range[1]`, `pre_window[0] < pre_window[1]`, `--window > 0.0`,
`--bin-width > 0.0`).

## 7. Module structure (single file)

```
benchmark_hycal_timing_deconv.py
├── module docstring (deps, usage, output artifacts, Phase 1 scope note)
├── imports (stdlib: argparse, csv, math, sys, time, dataclasses, pathlib;
│            third-party: numpy, matplotlib; prad2py: dec, det;
│            optional: scipy.optimize.curve_fit via HAVE_SCIPY)
├── constants + dataclasses
│   ├── CONFIG_NAMES = ("legacy_largestE", "legacy_smallestDt",
│   │                   "gated_largestE",  "gated_smallestDt")
│   ├── TSV_COLUMNS  (identical 14-column schema as study_hycal_timing.py)
│   ├── @dataclass PeakRecord
│   │     module_index: int
│   │     module_id:    int
│   │     module_name:  str
│   │     module_type:  str                  # mod.type.name
│   │     energy_MeV:   float                # from out.integral[k] via mod.energize
│   │                                        # (both Q_DECONV_APPLIED and
│   │                                        # Q_DECONV_FALLBACK_GLOBAL produce
│   │                                        # valid out.integral[]; see §14)
│   │     time_ns:      float                # wres.peaks[k].time (CFD-timed)
│   │     deconv_state: int                  # out.state (one of Q_DECONV_*)
│   │     chi2_per_dof: float                # out.chi2_per_dof
│   │     used_fallback: bool                # True iff deconv state is neither
│   │                                        # Q_DECONV_APPLIED nor
│   │                                        # Q_DECONV_FALLBACK_GLOBAL, or
│   │                                        # out.n != wres.npeaks — in which
│   │                                        # case pk.integral was used
│   ├── @dataclass ConfigAccumulator
│   │     name:          str                 # one of CONFIG_NAMES
│   │     hist:          np.ndarray[int64]   # length = N_bins over --e-range
│   │     n_events:      int                 # events contributing a 1-cluster hit
│   │     energies:      list[float]         # per-event cluster-E in MeV, retained
│   │                                        # only with --calibrate (mirrors
│   │                                        # benchmark_hycal_timing.py:327)
│   │     cal_events:    list[tuple]         # per-event (seed_module_id, x, y,
│   │                                        # row, col, seed_energy_MeV,
│   │                                        # cluster_energy_MeV) — populated
│   │                                        # only with --calibrate; see
│   │                                        # benchmark_hycal_timing.py:328
│   │     fit_mu:        float = nan
│   │     fit_sigma:     float = nan
│   │     fit_mu_err:    float = nan
│   │     fit_sigma_err: float = nan
│   │     cb_mu, cb_sigma, cb_alpha, cb_n: float = nan   # populated only with --calibrate
│   │     (per-config deconv-health counters removed in iteration 2 —
│   │      counters are now run-level; see §14 and RunContext below)
│   └── @dataclass RunContext
│         e_lo, e_hi, bin_width,             # histogram grid
│         window, pre_window,                # timing cuts
│         require_trigger: int,              # from --require-trigger (default 0)
│         max_fired_modules, max_qdist,      # event-level cuts
│         signal_min,                        # single-cluster-E threshold
│         calibrate: bool, target_peak, inner_radius,
│         min_seed_energy, min_per_module,
│         exclude_rowcol, cal_method,
│         scipy_available: bool,
│         template_file: Path,
│         e_edges: np.ndarray,               # derived: np.arange(e_lo, e_hi+bw, bw)
│         # ---- run-level mutable tallies, populated during the event loop ----
│         # (see §14 for deconv-health semantics, §13 for fired-module cut)
│         n_events_total:            int = 0    # events seen after trigger filter
│         n_events_fired_module_cut: int = 0    # events skipped by --max-fired-modules
│         n_peaks_deconv_per_channel: int = 0   # out.state == Q_DECONV_APPLIED
│         n_peaks_deconv_global:      int = 0   # out.state == Q_DECONV_FALLBACK_GLOBAL
│         n_peaks_deconv_failed:      int = 0   # anything else, OR out.n != wres.npeaks
│         chi2_sum:                   float = 0.0  # running sum across both success states
│         chi2_count:                 int   = 0
├── load_templates(path, wave_cfg) -> (store, dict[str, PulseTemplate])
│     Build the store once, pre-fetch per-type templates via
│     store.type_template("PbWO4"), store.type_template("PbGlass").
│     Return the dict so the per-channel inner loop skips the store
│     lookup on every pulse.
│     Missing-file / parse-error policy (see §15.3): the LoadFromFile
│     call returns False on any such failure; load_templates catches
│     that, logs one WARN line, and returns ({}, {}).  Does NOT
│     raise.  Downstream, store_by_type.get(mod.type.name) returns
│     None for every channel, which falls through to the pk.integral
│     fallback branch in §9.
│     If the file loads but a per-type template is missing, log one
│     WARN line and set that entry to None; channels of that type
│     take the same pk.integral fallback path.
├── extract_peaks(wa, store_by_type, cd, mod, ctx) -> list[PeakRecord]
│     The per-channel pipeline per §9.  Returns zero or more
│     PeakRecord entries.  Isolated so --self-test can exercise it.
├── pre_filter_smallest_dt(peaks, seed_time) -> list[PeakRecord]
│     Group peaks by module_index, pick argmin |time_ns - seed_time| per
│     module.  See §10 for the seed-time definition.
├── run_configs(peaks, ctx) -> dict[str, Cluster] | None
│     Build four HyCalCluster instances, feed each the appropriate
│     filtered peak stream, call form_clusters().  Returns the four
│     cluster lists keyed by CONFIG_NAMES, or None if peaks is empty.
│     See §11 for the four-config assembly.
├── accumulate_spectrum(clusters_by_cfg, accs, ctx)
│     For each config: if exactly one cluster passes --signal-min, add
│     to accs[cfg].hist and bump accs[cfg].n_events.
├── fit_gaussian_per_config(acc, ctx)
│     Standard gauss(x; A, mu, sigma) curve_fit on acc.hist with
│     --signal-min as a left cut-off (same window logic as
│     benchmark_hycal_timing.py:220–234).  Populates fit_* fields.
├── fit_crystal_ball_per_config(acc, ctx)
│     Only called under --calibrate.  Mirrors benchmark_hycal_timing.py
│     crystal_ball() + bounds exactly.
├── calibrate_per_module(peaks_by_seed_mod, ctx) -> dict[int, float]
│     Per-seed-module gain refinement, same algorithm as
│     benchmark_hycal_timing.py's calibration routine.  Returns a
│     module_index → gain-correction-factor dict.
├── plot_spectrum(accs, ctx, out_path)
│     Four overlaid histograms (one per config) with per-config
│     Gaussian fit curves; legend shows σ_E and σ_E/E for each
│     config at the fitted μ.  Format consistent with
│     benchmark_hycal_timing.py's existing spectrum PNG.
├── plot_calibrated(accs, ctx, out_path)
│     Only under --calibrate.  Same layout as plot_spectrum but
│     overlaid with CB fits; mirrors benchmark_hycal_timing.py's
│     <out>_calibrated.png.
├── write_dt_tsv(all_rows, path)
│     14-column TSV identical to study_hycal_timing.py.  One row per
│     (seed-pulse, neighbor-pulse) pair produced by
│     CollectNeighborTiming.  Shared across all four configs
│     (CollectNeighborTiming runs once per event on the raw peak
│     stream before any config-specific pre-filter).
├── write_summary_csv(accs, ctx, path)
│     One row per config.  Schema in §8.3.  Comment header lines
│     (starting with `#`) carry: EVIO path, template-file path +
│     sha1, run-num, trigger-mask, window, pre-window, signal-min,
│     e-range, bin-width, calibration flags, scipy availability,
│     timestamp.
├── _validate_args(ns)
└── main(argv=None) -> int
      1. parse args, validate
      2. setup_pipeline (via _common.setup_pipeline) → EVIO stream,
         HyCalSystem, WaveAnalyzer, crate_map; build four HyCalCluster
         instances
      3. load_templates(args.template_file, p.wave_cfg)
      4. per-event loop (§9) — accumulate dt_rows + hits_by_cfg
      5. per-config: fit_gaussian_per_config; if --calibrate also
         calibrate_per_module + fit_crystal_ball_per_config
      6. write all four artifacts in a stable order:
            <out>_dt.tsv
            <out>_spectrum.png
            <out>_summary.csv
            <out>_calibrated.png    (only with --calibrate)
      7. print one-line summary + return 0
```

## 8. Output artifacts

### 8.1 `<out>_dt.tsv`

TSV, 14-column schema identical to `study_hycal_timing.py`:

```
event_num, trigger_bits,
seed_module, seed_id, seed_time, seed_energy,
neighbor_module, neighbor_id, neighbor_time, neighbor_energy,
dt, dx_q, dy_q, dist_q
```

One row per (seed, neighbour-pulse) pair returned by
`CollectNeighborTiming`. `CollectNeighborTiming` runs once per event
on the raw post-pre-filter peak stream (before the per-config
pre-filter), so this TSV captures the full picture; downstream
analysis can group by (seed_module, neighbor_module) and apply
smallest-|Δt| in post-processing to compare the picking policy
without re-running the pipeline.

### 8.2 `<out>_spectrum.png`

Single-cluster energy spectrum for all four configs overlaid, with
per-config Gaussian fits. Legend shows `σ_E (MeV)` and
`σ_E/E (%)` at each fitted μ. Format consistent with
`benchmark_hycal_timing.py`'s existing spectrum PNG — same colour
pairing (legend-mode-to-colour map): `legacy_largestE` and
`legacy_smallestDt` in cool tones; `gated_largestE` and
`gated_smallestDt` in warm tones. Vertical reference line at
`--target-peak` under `--calibrate`.

### 8.3 `<out>_summary.csv`

One row per config, in this fixed column order:

| column                      | meaning                                                                                                     |
|-----------------------------|-------------------------------------------------------------------------------------------------------------|
| `config_name`               | one of `legacy_largestE`, `legacy_smallestDt`, `gated_largestE`, `gated_smallestDt`                          |
| `n_events`                  | #events contributing a one-cluster hit above `--signal-min`                                                  |
| `n_in_fit_window`           | #events in the Gaussian fit window (post-signal cut, within `--e-range`)                                     |
| `mu_MeV`                    | Gaussian μ (NaN under `scipy-missing`)                                                                      |
| `sigma_MeV`                 | Gaussian σ                                                                                                   |
| `sigma_over_E_pct`          | 100 × σ / μ                                                                                                  |
| `mu_calibrated`             | Gaussian μ post-per-module gain refinement (NaN unless `--calibrate`)                                       |
| `sigma_calibrated`          | Gaussian σ post-calibration (NaN unless `--calibrate`)                                                      |
| `sigma_CB`                  | Crystal-Ball σ post-calibration (NaN unless `--calibrate`)                                                  |
| `sigma_over_E_CB`           | 100 × σ_CB / μ_calibrated (NaN unless `--calibrate`)                                                        |

**Deconv-health counters are run-level**, not per-config, and are
emitted as `#`-prefixed preamble lines at the top of the CSV
(not as columns). See §14.3 for the exact preamble format. Rationale:
deconv runs on the shared raw peak list before any config-specific
selection, so the three counters would be identical across the four
rows if emitted as columns.

Full preamble layout (lines starting with `#`):

- EVIO path, template-file path + sha1, run-num, trigger-mask
  (`--require-trigger` value), `--window`, `--pre-window`,
  `--signal-min`, `--e-range`, `--bin-width`, `--max-fired-modules`,
  `--max-qdist`, `--cal-method`, `--target-peak`, scipy
  availability, timestamp.
- Run-level deconv-health block per §14.3:
  `n_peaks_deconv_per_channel`, `n_peaks_deconv_global`,
  `n_peaks_deconv_failed`, `deconv_failed_fraction`, `mean_chi2_per_dof`.
- Run-level event-count block per §13:
  `n_events_total`, `n_events_fired_module_cut`,
  `fired_module_cut_fraction`.

### 8.4 `<out>_calibrated.png` (`--calibrate` only)

Same four-config overlay as `<out>_spectrum.png` but with the
per-module-gain-corrected spectrum and CB fits overlaid. Mirrors
`benchmark_hycal_timing.py` layout one-for-one.

## 9. Pipeline per event

Executed inside the EVIO loop, mirroring
`study_hycal_timing.py:140–210` with the deconv insertion. In
order:

1. **Trigger filter.** Skip events whose `trigger_bits` do not match
   `--require-trigger`. Default `0` keeps physics unbiased.
2. **Collect peaks per channel** (one pass over all ROCs/slots/channels):
   ```python
   peaks: list[PeakRecord] = []
   fired_module_set: set[int] = set()
   for ri in range(fadc_evt.nrocs):
       roc = fadc_evt.roc(ri)
       if not roc.present: continue
       crate = p.crate_map.get(roc.tag)
       if crate is None: continue
       for s in roc.present_slots():
           slot = roc.slot(s)
           for c in slot.present_channels():
               mod = p.hycal.module_by_daq(crate, s, c)
               if mod is None or not mod.is_hycal():
                   continue                            # drops LMS + Veto; see §12
               cd = slot.channel(c)
               if cd.nsamples <= 0:
                   continue
               wres = p.wave_ana.analyze_result(cd.samples)
               if wres.npeaks == 0:
                   continue

               tmpl = store_by_type.get(mod.type.name)   # may be None if template missing
               out  = (p.wave_ana.deconvolve(cd.samples, wres, tmpl)
                       if tmpl is not None else None)
               # Deconv success states (see §14 for the convention):
               #   Q_DECONV_APPLIED         — per-channel template fit converged
               #   Q_DECONV_FALLBACK_GLOBAL — per-type aggregate template fit
               #                              converged (THIS is the expected
               #                              success state for Phase 1 because
               #                              load_from_file loads per-type
               #                              templates only; see §14)
               # Both populate out.integral[] / out.height[] identically.
               used_deconv = (
                   out is not None
                   and out.state in (dec.Q_DECONV_APPLIED,
                                     dec.Q_DECONV_FALLBACK_GLOBAL)
                   and out.n == wres.npeaks
               )

               # Run-level deconv-health tally (§14).  Done once per channel-
               # per-event, BEFORE the per-peak pre-window cut, so the counts
               # reflect what the deconvolver did on the raw peak set.
               if out is not None and out.n == wres.npeaks:
                   if out.state == dec.Q_DECONV_APPLIED:
                       ctx.n_peaks_deconv_per_channel += out.n
                       ctx.chi2_sum   += out.chi2_per_dof * out.n
                       ctx.chi2_count += out.n
                   elif out.state == dec.Q_DECONV_FALLBACK_GLOBAL:
                       ctx.n_peaks_deconv_global += out.n
                       ctx.chi2_sum   += out.chi2_per_dof * out.n
                       ctx.chi2_count += out.n
                   else:
                       ctx.n_peaks_deconv_failed += out.n
               else:
                   # tmpl was None, or out.n != wres.npeaks — count every raw
                   # peak as a deconv failure for health-reporting purposes.
                   ctx.n_peaks_deconv_failed += wres.npeaks

               for k in range(wres.npeaks):
                   pk = wres.peaks[k]
                   if pk.time <= ctx.pre_window[0]: continue
                   if pk.time >= ctx.pre_window[1]: continue
                   if used_deconv:
                       adc = out.integral[k]
                   else:
                       adc = pk.integral                 # pk.integral fallback — same
                                                         # field sibling scripts use
                                                         # (study_hycal_timing.py:180,
                                                         # benchmark_hycal_timing.py:399)
                   peaks.append(PeakRecord(
                       module_index = mod.index,
                       module_id    = mod.id,
                       module_name  = mod.name,
                       module_type  = mod.type.name,
                       energy_MeV   = mod.energize(adc),
                       time_ns      = float(pk.time),
                       deconv_state = (out.state if out is not None
                                       else dec.Q_DECONV_NOT_RUN),
                       chi2_per_dof = (out.chi2_per_dof if used_deconv
                                       else float("nan")),
                       used_fallback = not used_deconv,
                   ))
                   fired_module_set.add(mod.index)
   ```
   **Why the fallback.** Explicitly calling `wa.deconvolve` bypasses
   the production `nnls_deconv.enabled` gate
   (`WaveAnalyzer.cpp:832–835`), so the call always runs. But the
   per-pulse LM fit may still fail (returning `Q_DECONV_LM_NOT_CONVERGED`,
   `Q_DECONV_NO_TEMPLATE`, `Q_DECONV_BAD_TEMPLATE`, `Q_DECONV_FAILED`,
   or `Q_DECONV_DISABLED`); when it does we use `wres.peaks[k].integral`
   from the C++ `findPeaks` pipeline instead. The sibling scripts
   (`study_hycal_timing.py:180`, `benchmark_hycal_timing.py:399`) use
   this same `pk.integral` field as the amplitude input to
   `mod.energize(...)`, so the fallback branch is numerically identical
   to the sibling scripts' baseline. The deconv-health counters (§14)
   record how often each state fired.
3. **Fired-module count cut.**
   ```python
   if len(fired_module_set) > ctx.max_fired_modules:
       continue                                       # skip event entirely
   ```
   Applied AFTER all channels have been scanned and pulse-pre-filtered,
   BEFORE `CollectNeighborTiming` or any clusterer call. See §13 for
   ordering rationale.
4. **CollectNeighborTiming on the raw peak stream.** Push every
   surviving `PeakRecord` into a single scratch `HyCalCluster`,
   call `collect_neighbor_timing(args.max_qdist)`, emit one TSV row
   per returned `SeedNeighborTiming`. This is the diagnostic TSV
   (§8.1); it is independent of the four-config spectrum.
5. **Build and fill the four configs.** See §11.
6. **Spectrum selection + calibration record.** For each config: if
   exactly one cluster passes `--signal-min`, add that cluster's
   energy to `accs[cfg].hist` and bump `accs[cfg].n_events`. Under
   `--calibrate`, additionally:
   - `accs[cfg].energies.append(cluster_energy_MeV)` — bulk energy
     list used for the per-config Gaussian + CB fit.
   - `accs[cfg].cal_events.append((seed_module_id, seed_x_mm,
     seed_y_mm, seed_row, seed_col, seed_energy_MeV,
     cluster_energy_MeV))` — per-event tuple consumed by
     `calibrate_per_module(...)` for the per-seed-module gain
     refinement. Tuple layout mirrors
     `benchmark_hycal_timing.py:325–328` one-for-one so the
     downstream calibration helper can be copied without
     adaptation.
   Both lists stay empty when `--calibrate` is off — the memory
   cost is opt-in.
7. **Progress line** every 5000 events (same cadence as
   `study_hycal_timing.py`).

## 10. Smallest-|Δt| selection algorithm

Used twice, with slightly different seed-time anchors:

- **Diagnostic TSV (§8.1):** no Python-side smallest-|Δt| picking
  happens before writing. `CollectNeighborTiming` emits one row per
  (seed, neighbor-pulse) pair, and the downstream consumer applies
  smallest-|Δt| in post-processing (just like it would for
  `study_hycal_timing.py`'s TSV). Keeps the two scripts' TSV schema
  identical.

- **Python pre-filter for `gated_smallestDt` (§11):** run inside
  `pre_filter_smallest_dt(peaks, seed_time)`.
  1. **Seed-time definition.** Use the single largest-amplitude peak
     across all surviving modules as the Python-side seed-time anchor,
     with a deterministic tie-break so the choice is stable across
     reruns on the same input:
     ```python
     seed = max(peaks,
                key=lambda p: (p.energy_MeV,
                               -p.module_index,
                               -p.time_ns))
     seed_time = seed.time_ns
     ```
     Primary key is `energy_MeV` (physics-motivated); secondary
     `-module_index` prefers lower module indices on exact ties
     (reproducible across the HyCal geometry); tertiary `-time_ns`
     prefers earlier peaks on the still-rare double tie. This
     mirrors what the C++ global seed-finder does
     (`HyCalCluster::FindSeeds` picks the global energy maximum
     first) while making the Python-side policy deterministic.
  2. **Per-module smallest-|Δt| pick.** Group `peaks` by
     `module_index`; for each module keep the single peak with
     minimum `(|time_ns - seed_time|, -energy_MeV)` (primary:
     closeness in time; secondary: larger energy on an exact Δt
     tie).
  3. Return the filtered list, ready for `add_hit` on the
     `gated_smallestDt` clusterer.

Behavioural edge cases the implementer should keep in mind (not
asserted in-code; a `--self-test` harness is out of Phase 1 scope):
zero peaks → empty list, no exception; one peak per module →
unchanged list; two peaks on one module, one clearly closer to the
global largest-amplitude peak → the closer one is kept.

## 11. Four-config clustering comparison

Build four `HyCalCluster` instances once per job (reused across
events; call `.clear()` at the top of each event). `seed_time_window`
is a field on `HyCalClusterConfig`, not a direct attribute of
`HyCalCluster` — it is set via `cfg.seed_time_window = ...` on a
freshly-constructed config that mirrors the pipeline's base config,
then installed with `cl.set_config(cfg)`. This is the exact pattern
used in `benchmark_hycal_timing.py:309–322`:

```python
# Four clusterers sharing the HyCal geometry; differ only in
# seed_time_window.  Base config is sourced from the pipeline's
# default clusterer so the min_* / split / log-weight knobs match
# whatever daq_config.json / HyCalClusterConfig defaults were loaded.
base_cfg = p.hc_clusterer.get_config()

def _make_cluster(seed_window: float) -> det.HyCalCluster:
    cfg = det.HyCalClusterConfig()
    for f in ("min_module_energy", "min_center_energy",
              "min_cluster_energy", "min_cluster_size",
              "corner_conn", "split_iter", "least_split",
              "log_weight_thres"):
        setattr(cfg, f, getattr(base_cfg, f))
    cfg.seed_time_window = seed_window
    cl = det.HyCalCluster(p.hycal)
    cl.set_config(cfg)
    return cl

cl_legacy_bigE  = _make_cluster(-1.0)
cl_legacy_smDt  = _make_cluster(-1.0)
cl_gated_bigE   = _make_cluster(ctx.window)
cl_gated_smDt   = _make_cluster(ctx.window)
```

The legacy configs use `seed_time_window = -1.0` (sentinel for "no
timing gate"); the gated configs use `ctx.window` (default 4.0 ns).
All four share the same base `min_*` / split / log-weight settings
so differences in the output spectra are attributable only to the
timing-gate × neighbour-selection factor under study.

Per-event feed:

| config               | input stream                                                              | why                                                                                                                                                                                                                                 |
|----------------------|---------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|
| `legacy_largestE`    | `max_energy_per_module(peaks)`                                            | upstream largest-integral-per-module reduction, same as the legacy `HyCalCluster` behaviour when `seed_time_window ≤ 0`.                                                                                                            |
| `legacy_smallestDt`  | `pre_filter_smallest_dt(peaks, seed_time)`                                | upstream smallest-|Δt|-per-module reduction. Seed time is the global largest-amplitude peak time (§10).                                                                                                                             |
| `gated_largestE`     | all `peaks` (every surviving pulse pushed as a separate `ModuleHit`)      | C++ BFS at `HyCalCluster.cpp:197–209` applies the seed-anchored timing coincidence cut using its default largest-energy neighbour-selection policy. This is the config that matches `benchmark_hycal_timing.py`'s existing "gated". |
| `gated_smallestDt`   | `pre_filter_smallest_dt(peaks, seed_time)`                                | C++ BFS still runs its default largest-E logic, but there is only one pulse per module in the input, so the policy is **effectively smallest-|Δt|**. See the ambiguity note below.                                                 |

**Ambiguity resolution (committed).** The C++ cluster BFS at
`HyCalCluster.cpp:197–209` always does largest-energy pulse selection
and has no runtime switch. The user-brief laid out two options:

- (a) Pre-select pulses in Python before `add_hit` so the C++ clusterer
  only sees one candidate per module (policy flip happens on the Python
  side).
- (b) Call `CollectNeighborTiming` directly, do smallest-|Δt| in Python,
  then rebuild clusters by hand bypassing `FormClusters`.

**Phase 1 uses option (a).** Option (b) rewrites the clusterer in
Python and loses the C++ BFS's shower-shape / dead-module / sector-grid
logic — not worth the fidelity hit for a measurement.

Call `form_clusters()` on each of the four instances after the hit
stream has been pushed, then forward the four cluster lists to
`accumulate_spectrum(...)`.

**Known approximation (not fixed in Phase 1).** The Python-side
smallest-|Δt| pre-filter uses a single global seed anchor per event
(the largest-amplitude peak across all surviving modules; §10). In
multi-cluster events (which the single-cluster `--signal-min` cut
is designed to reject but which can still influence the pre-filter
before that cut applies), "closest-to-event-largest-pulse" is only
an approximation of the physically ideal "closest-to-each-cluster-
seed" policy. The headline single-cluster spectrum is dominated by
events with a well-defined global seed, so the approximation is
expected to be adequate for the Phase 1 question — but it is a known
limitation to call out when interpreting the result (§15.2). A
per-cluster seed anchor would require re-running the clusterer in
Python, which was ruled out by the §11 ambiguity resolution.

## 12. Channel-type scope — `is_hycal()` convention

Phase 1 processes only HyCal physics calorimeter modules (PbWO4 +
PbGlass) via `mod.is_hycal()`, which also incidentally excludes LMS
and Veto modules. This matches the convention used by sibling
scripts `study_hycal_timing.py:168` and
`benchmark_hycal_timing.py:379`. From `prad2det/include/HyCalSystem.h:174`:

```cpp
bool is_hycal() const { return is_pwo4() || is_glass(); }
```

Which returns true iff the module's `ModuleType` is `PbWO4` or
`PbGlass`. LMS and Veto are excluded by construction — no new
predicate is needed; reuse the same idiom. The channel filter in §9
step 2 is literally the sibling scripts' line, copy-pasted:

```python
if mod is None or not mod.is_hycal():
    continue
```

Verified facts about the Python binding (`python/bind_det.cpp`):
- `mod.type` returns a `det.ModuleType` enum value (not a string).
  Canonical enum→string conversion is `mod.type.name` (used elsewhere
  in `fit_pulse_template.py:946`, `pileup_generator.py:572`).
- `mod.is_pwo4()`, `mod.is_glass()`, `mod.is_hycal()` are all bound
  helpers (`bind_det.cpp:554–556`). No `is_lms()` / `is_veto()`
  helper exists, which is fine — the complement of `is_hycal()`
  already covers both.
- The template-store lookup `store.type_template("PbWO4" | "PbGlass")`
  takes a string. Pass `mod.type.name` directly.

The template file `pulse_templates_025308_h50_RawSum_max20modules.json`
is missing an `LMS` entry in `_by_type`, which is **not an issue**:
because `is_hycal()` has already dropped LMS channels at the top of
the loop, `store.type_template("LMS")` is never called.

**One edge case to log.** If a `PbGlass` or `PbWO4` channel has
`mod.type.name` returning anything other than `"PbGlass"` or
`"PbWO4"` (should not happen; defensive), the `store_by_type.get(...)`
returns `None` and the pulse falls into the `pk.integral` fallback
branch (§9 step 2). The run-level `n_peaks_deconv_failed` counter
(§14.2) absorbs these, and the script emits one WARN line per
unique unknown type-name over the whole job.

## 13. Fired-module count cut

Ordering: this cut sits **after** per-channel pulse extraction and
**before** any cluster call (step 3 of §9). This ordering is
intentional:

- **Why after pulse extraction**, not before (e.g. channel-counting
  the ADC-non-zero channels in the raw EVIO payload): we want the
  cut to track the number of channels that produced a *real-looking
  pulse* passing the physics pre-window, not the number of channels
  that happened to have a saturated pedestal or a tiny ringing
  artifact. The post-pre-filter count is a cleaner proxy for
  "something physics-like fired".

- **Why before CollectNeighborTiming + the four clusterers**: both
  are O(N_hits) or worse, so skipping the event here saves work and
  avoids polluting the TSV + the energy spectra with LMS-flash
  events.

Rationale for the default `--max-fired-modules 50`:

- A clean 2.2 GeV e-p shower is typically ~10–30 fired HyCal
  modules (seed + first/second shower rings).
- An LMS laser flash illuminates physics channels simultaneously
  with the LMS reference channels. Even though the LMS reference
  channels themselves are already dropped by `is_hycal()`, the
  flash-induced response on surrounding PbWO4 / PbGlass channels is
  not — hence the count can balloon.
- Beam splash / high-rate pile-up also inflates the count.
- 50 is permissive enough to accept all normal e-p + Møller events
  with headroom, and tight enough to catch obvious LMS-flash and
  splash contamination. Expose as a CLI knob so the threshold can be
  tuned without a code change.

The cut skips the entire event (does not emit TSV rows, does not
contribute to any spectrum histogram). Both the count and the
fraction are logged to the summary CSV comment-header preamble
(§8.3) as:

```
# n_events_total=<N>
# n_events_fired_module_cut=<M>
# fired_module_cut_fraction=<M/N>
```

Report the fraction to two decimal places so a user eyeballing the
preamble sees immediately whether the backstop is firing on noise
(e.g. 0.01) or whether it is reshaping the dataset (e.g. 0.15).
Both counters live on `RunContext` (§7).

## 14. Deconv health reporting

### 14.1 Success-state convention

Per `prad2dec/src/WaveAnalyzer.cpp:1121–1122` and
`prad2dec/include/Fadc250Data.h:326–335`, the LM deconvolver returns
two distinct **success** states:

- `Q_DECONV_APPLIED` — the LM fit converged against a **per-channel**
  template (resolved by the store's `(roc_tag, slot, channel) →
  module_type` lookup). The C++ comment flags this as "currently
  unused" because the shipped templates are loaded as per-type
  aggregates, not per-channel entries.
- `Q_DECONV_FALLBACK_GLOBAL` — the LM fit converged against a
  **per-type aggregate** template obtained via `type_template("PbWO4")`
  / `type_template("PbGlass")`. **This is the expected success state
  for Phase 1**, because the user-specified template file
  (`pulse_templates_025308_h50_RawSum_max20modules.json`) populates
  only the `_by_type` block.

Both states produce valid `out.integral[]` / `out.height[]` /
`out.t0_ns[]` / `out.tau_r_ns[]` / `out.tau_f_ns[]` arrays indexed
0..`out.n-1`; downstream code reads them identically. Any other
state (`Q_DECONV_NOT_RUN`, `Q_DECONV_DISABLED`, `Q_DECONV_NO_TEMPLATE`,
`Q_DECONV_BAD_TEMPLATE`, `Q_DECONV_LM_NOT_CONVERGED`, `Q_DECONV_FAILED`),
or an `out.n != wres.npeaks` mismatch, is a **failure** and the
`pk.integral` fallback is used (§9).

### 14.2 Counters — run-level, not per-config

Deconv runs on the shared raw-peak list produced by
`extract_peaks(...)` BEFORE any config-specific pulse selection, so
per-config health counters would just duplicate the same run-wide
numbers. The counters live on `RunContext` (one set per run) and are
reported in the summary CSV **comment-header preamble**, not as
columns on the per-config rows.

Three run-level counters:

- `n_peaks_deconv_per_channel` — peaks where `out.state == Q_DECONV_APPLIED`.
  Expected to be **zero** in Phase 1 (no per-channel templates loaded).
  A non-zero value means a later template file added per-channel
  entries — worth flagging in the §15 interpretation.
- `n_peaks_deconv_global` — peaks where `out.state == Q_DECONV_FALLBACK_GLOBAL`.
  **Expected to dominate** in Phase 1 (per-type aggregate templates
  only).
- `n_peaks_deconv_failed` — peaks where the state is neither of the
  above, OR `out.n != wres.npeaks`, OR `tmpl is None` (template
  missing / parse failure). In this case `pk.integral` was used
  (see §9). Expected to be < a few % for a well-matched template.

Plus one derived quantity:

- `mean_chi2_per_dof` — arithmetic mean of `out.chi2_per_dof`
  weighted by peak count across both success states
  (`Q_DECONV_APPLIED` + `Q_DECONV_FALLBACK_GLOBAL`); failures do not
  contribute.

### 14.3 Reporting

Written as `#`-prefixed comment lines at the top of
`<out>_summary.csv` (§8.3 preamble) in a fixed format:

```
# n_peaks_deconv_per_channel=XXX
# n_peaks_deconv_global=YYY
# n_peaks_deconv_failed=ZZZ
# deconv_failed_fraction=0.FF
# mean_chi2_per_dof=X.XX
# n_events_total=NNN
# n_events_fired_module_cut=MMM
# fired_module_cut_fraction=0.GG   (see §13; this item is also run-level)
```

Also echoed once to stderr at end-of-run:

```
[deconv health] per_channel=<n>  global=<m>  failed=<k>
                mean_chi2/dof=<x.xx>  failed_frac=<pp.p>%
```

A `deconv_failed_fraction` above ~20% is a red flag that the
template file is poorly matched to the data run (wrong detector
state, wrong trigger mix, template extracted on a different
calibration epoch). This is diagnostic signal only — the script
does not stop or re-fit.

## 15. Validation — smoke test

### 15.1 One-file smoke test

```bash
python3 analysis/pyscripts/benchmark_hycal_timing_deconv.py \
    "$HOME/work/PRad/data/evio/prad_025320.evio.00099" \
    output/deconv_025320_smoke \
    --template-file analysis/pyscripts/output/pulse_templates_025308_h50_RawSum_max20modules.json \
    --max-events 50000 \
    --hc-calib-file "$PRAD2_DATABASE_DIR/calibration/calibration_factor_2p2_fadc.json" \
    --daq-config "$DAQ_CONFIG" \
    --hc-map-file "$HC_MAP"
```

Expected on success:
- Exit 0.
- Four files written:
  `output/deconv_025320_smoke_dt.tsv`,
  `output/deconv_025320_smoke_spectrum.png`,
  `output/deconv_025320_smoke_summary.csv`,
  (and `output/deconv_025320_smoke_calibrated.png` if `--calibrate`
  was passed).
- `_spectrum.png` shows four overlaid single-cluster energy
  distributions over 500–2500 MeV with four Gaussian fit curves and
  a legend showing σ_E for each config.
- `_summary.csv` has a `#`-prefixed preamble containing provenance
  + the run-level deconv-health block (§14.3) + the run-level
  event-count block (§13), followed by exactly four data rows (one
  per config) with non-NaN `mu_MeV`, `sigma_MeV`,
  `sigma_over_E_pct`.
- `_dt.tsv` first row header matches the 14-column schema; the body
  has one row per (seed-pulse, neighbor-pulse) pair.
- Runtime ~15 min at `--max-events 50000` on a single file
  (`wa.deconvolve` adds ~2× overhead per channel vs
  `benchmark_hycal_timing.py`'s raw-peak path). If the run exceeds
  ~30 min this is a signal that the deconv is taking the LM-not-
  converged path far more often than expected — rerun with `-v` and
  inspect the stderr deconv-health summary.

### 15.2 Result interpretation

- Preamble sanity check first: confirm
  `n_peaks_deconv_global` dominates and
  `n_peaks_deconv_per_channel == 0` (Phase 1 convention; see §14.1).
  A non-zero per-channel count means a later template file added
  per-channel entries — not wrong, but worth noting. A dominant
  `n_peaks_deconv_failed` means the deconv path is effectively off
  and the measurement collapses to the raw-integral baseline.
- `gated_smallestDt` σ_E/E ≤ `gated_largestE` σ_E/E — expected if
  closer-in-time pulses are better coincidence candidates
  (physically motivated). If this ordering is inverted, something
  is wrong with the Python-side seed-time anchor or the pre-filter.
  Caveat: the smallest-|Δt| policy uses a single global seed anchor
  per event (§11 "Known approximation"); multi-cluster-heavy
  samples will blur this comparison.
- Both `gated_*` σ_E/E ≤ both `legacy_*` σ_E/E — this is the key
  Phase 1 question. If true, template-fit + CFD cured the
  amplitude-dependent timing walk and `plan_waveform.md` Phase 3
  can recommend turning `nnls_deconv.enabled` on in production. If
  false (and the preamble sanity check above passed), Phase 2
  (matched filter in the detection loop) is next.
- `mean_chi2_per_dof` ≈ 1 — if much larger, the template file is
  off and this invalidates the amplitude comparison; try a
  different template file before drawing conclusions.

### 15.3 Degraded-mode checks

- Rerun on the same inputs with scipy uninstalled; confirm the TSV +
  PNG + CSV are still produced, with NaN fit columns in the summary
  CSV and a stderr `[WARN] scipy not available` line.
- Rerun with a nonsense `--template-file some/missing/file.json`.
  Behavior spec: `PulseTemplateStore.load_from_file(...)` returns
  False on file-not-found / parse-error / empty-contents;
  `load_templates(...)` (§7) catches this, logs one WARN line to
  stderr, and returns `({}, {})`. The per-channel loop in §9 then
  sees `tmpl = store_by_type.get(mod.type.name) == None` for every
  channel and falls through to the `pk.integral` fallback branch.
  **The script does not raise.** On a successful smoke-run with a
  missing template file: the run-level
  `n_peaks_deconv_failed` preamble counter (§14.3) should equal the
  total raw-peak count (every peak counted a failure); the
  `n_peaks_deconv_global` counter should be zero; the gated vs.
  legacy spectrum comparison still runs — just reduced to the
  baseline `benchmark_hycal_timing.py` configuration (raw peak
  integrals) with no template-fit improvement.

## 16. Deferred items

Captured here so they are not lost from the Phase 1 scope boundary;
tracked long-form in `plan_waveform.md`.

- **Explicit V1–V4 veto-counter handling** (e.g. in a downstream
  track-matching stage that needs V-counter inhibit windows).
  Phase 1 already drops all non-PbWO4/PbGlass modules via
  `mod.is_hycal()` (§12), so no veto-specific logic is needed inside
  this script. If a future downstream stage requires V-counter
  handling it is a separate script — Phase 3.
- **Literal matched-filter → template-fit → CFD pipeline.** Needs
  binding changes to make `WaveResult.npeaks` / `.peaks` writable
  (or an alternative "construct a WaveResult from candidate times"
  API), plus a public `compute_cfd_time(...)` C++ method. Full plan
  in `plan_waveform.md` §5 (Phase 2).
- **Change C++ cluster BFS to smallest-Δt policy.** Only worth doing
  if the Phase 1 2×2 grid shows `gated_smallestDt` is unambiguously
  better than `gated_largestE` — Phase 3.
- **Flipping `nnls_deconv.enabled = true` in
  `database/daq_config.json`.** Downstream production-config change
  recommended by `plan_waveform.md` Phase 3 if Phase 1 confirms the
  template-fit path cures the timing-gate regression. Not touched by
  this script.
- **Per-module-type `--window`.** PbWO4 timing is tighter than
  PbGlass; may want two separate `--window` values keyed by
  `mod.type.name`. Not in v1; a single global window is a correct
  starting point.

---

## Todo

Phase 1 implementation of `benchmark_hycal_timing_deconv.py`,
sequenced for incremental validation. Each item is specific enough
that an implementer can execute without re-reading the full plan.

### Phase A — scaffolding and I/O

- [ ] Create `analysis/pyscripts/benchmark_hycal_timing_deconv.py`
      with the module docstring, dependency list, and the Phase 1
      scope note from §3.
- [ ] Copy the `HAVE_SCIPY` try/except idiom and `gauss(...)` /
      `crystal_ball(...)` helpers from `benchmark_hycal_timing.py:44–60`.
- [ ] Define `CONFIG_NAMES`, `TSV_COLUMNS`, and the three dataclasses
      (`PeakRecord`, `ConfigAccumulator`, `RunContext`) per §7.
- [ ] Implement `_validate_args(ns)` with the bound checks from §6
      (window > 0, bin-width > 0, pre-window[0] < pre-window[1],
      e-range[0] < e-range[1]).
- [ ] Wire `_common.add_common_args` + the script-specific flags from
      §6 into `argparse`. Note: `--require-trigger` is NOT in
      `add_common_args` — declare it on the script's own argparse
      (hex-or-decimal parser, default 0 = unbiased physics sample).

### Phase B — pipeline bring-up

- [ ] Implement `load_templates(path, wave_cfg)` per §7: build the
      store, pre-fetch per-type templates into a dict keyed by
      type-name string.
- [ ] Implement `extract_peaks(wa, store_by_type, cd, mod, ctx)` per
      §9 step 2, treating BOTH `Q_DECONV_APPLIED` and
      `Q_DECONV_FALLBACK_GLOBAL` as deconv-success states, with the
      `out.n == wres.npeaks` sanity check, the `pk.integral` fallback
      branch, and the three run-level deconv-health counters
      (`n_peaks_deconv_per_channel`, `n_peaks_deconv_global`,
      `n_peaks_deconv_failed`) updated on `ctx`.
- [ ] Implement `pre_filter_smallest_dt(peaks, seed_time)` per §10
      including the three unit cases.
- [ ] Implement `max_energy_per_module(peaks)` helper for the
      `legacy_largestE` input stream.

### Phase C — clustering and spectrum accumulation

- [ ] Build the four `HyCalCluster` instances per §11 via the
      `HyCalClusterConfig` → `set_config(cfg)` pattern (NOT a direct
      `.seed_time_window = ...` assignment on the cluster object);
      base the config on `p.hc_clusterer.get_config()` to inherit
      all `min_*` / split / log-weight knobs; call `.clear()` at the
      top of each event.
- [ ] Implement `run_configs(peaks, ctx)` per §11 table: feed each
      clusterer the correct stream, call `form_clusters()`, return
      the four cluster lists.
- [ ] Implement `accumulate_spectrum(clusters_by_cfg, accs, ctx)`
      per §9 step 6, populating `accs[cfg].energies` and
      `accs[cfg].cal_events` only when `ctx.calibrate` is true.

### Phase D — event loop

- [ ] Write the main per-event loop per §9: trigger filter → extract
      peaks → fired-module count cut → CollectNeighborTiming →
      run_configs → accumulate_spectrum → progress-line cadence of
      5000 events.
- [ ] Push diagnostic TSV rows per event (one row per
      `SeedNeighborTiming` from `CollectNeighborTiming` on the raw
      peak stream, before any config pre-filter).

### Phase E — fits and output artifacts

- [ ] Implement `fit_gaussian_per_config(acc, ctx)` mirroring
      `benchmark_hycal_timing.py:220–234` fit window logic.
- [ ] Implement `write_dt_tsv(all_rows, path)` per §8.1 (14-column
      schema).
- [ ] Implement `plot_spectrum(accs, ctx, out_path)` per §8.2 with
      the four-config legend and σ_E / σ_E/E annotations.
- [ ] Implement `write_summary_csv(accs, ctx, path)` per §8.3,
      including the `#`-prefixed preamble with provenance, the
      run-level deconv-health block (§14.3), and the run-level
      event-count block (§13).

### Phase F — calibration path (`--calibrate`)

- [ ] Implement `calibrate_per_module(peaks_by_seed_mod, ctx)`
      mirroring `benchmark_hycal_timing.py` per-module gain
      refinement.
- [ ] Implement `fit_crystal_ball_per_config(acc, ctx)` mirroring
      `benchmark_hycal_timing.py` CB parameterisation.
- [ ] Implement `plot_calibrated(accs, ctx, out_path)` per §8.4.

### Phase G — validation

- [ ] Run the §15.1 smoke test on `prad_025320.evio.00099` with
      `--max-events 50000`; confirm the four artifacts and expected
      shape of the spectrum PNG.
- [ ] Run the §15.3 degraded-mode checks (scipy missing;
      template-file missing).
- [ ] Record the four `sigma_over_E_pct` values and the ordering
      against `benchmark_hycal_timing.py`'s baseline in
      `analysis/pyscripts/Commands.md` and `plan_waveform.md` §8
      change log.
- [ ] Update `plan_waveform.md` §2 "Current state" to mark Phase 1
      as implemented, and §4 "Phase 1" with the measured outcome
      (which of the three predicted branches fired).

### Phase H — follow-up decisions (planner artifacts, not code)

- [ ] If `gated_*` ≤ `legacy_*`: open `plan_waveform.md` Phase 3
      item for the `nnls_deconv.enabled = true` production config
      change; capture the measured σ_E/E improvement as the
      justification.
- [ ] If `gated_*` > `legacy_*` persists: open `plan_waveform.md`
      Phase 2 work — plan the binding changes needed to put the
      matched filter in the detection loop.
- [ ] If `gated_smallestDt` < `gated_largestE` consistently and
      significantly: open `plan_waveform.md` Phase 3 item for the
      C++ cluster BFS smallest-Δt policy change at
      `HyCalCluster.cpp:197–209`.
