# Plan — `module_dt_diagnostic.py`

**Status:** draft (iteration 10, pre-approval — §7.b flag vocabulary deferral, Phase 8 precondition deferral)
**Scope:** per-module ΔT = t_neighbour − t_seed diagnostic for HyCal, Python-only,
consumes the existing TSV produced by `analysis/pyscripts/study_hycal_timing.py`.
**Hard constraint honored:** no C++ changes; no modifications to `study_hycal_timing.py`
or `HyCalCluster`.

## Research inputs used

- `analysis/pyscripts/study_hycal_timing.py` (TSV producer; invocation + columns)
- `prad2det/src/HyCalCluster.cpp` lines 725–797 (`CollectNeighborTiming` semantics)
- `prad2det/include/HyCalCluster.h` lines 378–399 (`SeedNeighborTiming` struct)
- `prad2det/include/HyCalTimeCalib.h` lines 1–53 (per-module offset JSON schema)
- `database/hycal_time_offsets/time_calib.json` (keyed by `name` and/or `module_id`)
- `analysis/pyscripts/_common.py` lines 364–392, 417–429 (CLI + I/O conventions)
- `analysis/pyscripts/benchmark_hycal_timing.py` lines 44–235 (scipy Gaussian fit
  precedent; `HAVE_SCIPY` graceful-fallback pattern)
- `analysis/pyscripts/pileup_generator.py` lines 107–119 (dataclass+argparse style)
- `scripts/requirements.txt` and `calibration/requirements.txt` (project baseline
  `matplotlib`, `numpy`; neither lists `scipy` or `pandas`)
- `docs/research-study-hycal-timing-time-offset.md` (end-to-end
  confirmation that the TSV's `seed_time`, `neighbor_time`, and `dt`
  columns are raw / uncorrected; see §2.a Confirmation note)

---

## 1. Purpose

Given a run's per-pair timing TSV, flag HyCal modules whose observed `dt`
distribution against neighbouring seeds is inconsistent with a real
seed-anchored coincidence. Three failure modes we want to see:

1. **No peak** — ΔT distribution is essentially flat (pure accidentals): the
   module either never participates in real coincidences or its timing is
   randomized (dead/misconfigured timing).
2. **Offset peak** — ΔT distribution has a Gaussian peak but μ is far from 0:
   indicates an uncorrected per-module delay (missing/wrong row in the
   `hycal_time_offsets/*.json` table).
3. **Anomalously wide peak** — σ much larger than the detector baseline:
   indicates noisy/jittery timing (bad channel, cable issue).

The headline deliverable is a 2-D heatmap (y = module ID, x = ΔT, color =
counts) in which these three signatures jump out visually, plus a per-module
Gaussian + flat-background fit that quantifies them for downstream
flagging.

Accidentals-rate cross-check from singles rates
(`R_seed · R_module · 2·Δt_window`) is **out of scope** — see §12.

## 2. Input TSV — confirmed schema

Produced by `analysis/pyscripts/study_hycal_timing.py` (header at lines 100–105,
written by `C.open_table_writer` from `_common.py:417`).

**Columns (TSV by default; `--csv` switches to CSV):**

```
event_num, trigger_bits,
seed_module, seed_id, seed_time, seed_energy,
neighbor_module, neighbor_id, neighbor_time, neighbor_energy,
dt, dx_q, dy_q, dist_q
```

Key facts the implementer must respect:

- `*_module` is the HyCal internal index (`HyCalSystem::module(i).index`);
  `*_id` is the PrimEx ID (G: 1–576 for PbWO4, W: 1001–2152 for PbGlass).
  **Group and label modules by `neighbor_id`** (user-facing, matches the
  JSON offset table and other PRad tooling), not by the internal index.
- `dt`, `seed_time`, `neighbor_time` are ns; `*_energy` are MeV.
- `dx_q`, `dy_q`, `dist_q` are quantized distances in module units.
- One row per (seed, neighbour) pair, after the TSV producer has already:
  - filtered to physics-trigger events,
  - selected seeds by `min_center_energy` and the ±5 ns seed-claim window
    (`HyCalCluster.cpp:759–774`),
  - restricted neighbours to |dx_q|, |dy_q| ≤ `--max-qdist` (default 5).
- **No timing cut was applied by the producer** (that is the whole point of
  the study tool), so the ΔT distribution we read is the raw landscape.

### 2.a Offset-correction assumption — IMPORTANT

Looking at `study_hycal_timing.py:179–181`:

```python
p.hc_clusterer.add_hit(mod.index,
                        mod.energize(pk.integral),
                        float(pk.time))     # <-- raw peak time, from WaveAnalyzer
```

The script pushes the **raw** `pk.time` into the clusterer, and
`CollectNeighborTiming` just differences them (`HyCalCluster.cpp:789`:
`s.dt = hits_[k].time - st`). `Module::time_offset` is **never subtracted**
anywhere in this code path.

**Therefore the `dt` column in the TSV is uncorrected.** The μ this
diagnostic measures per neighbour module m is approximately

```
μ_m  ≈  (offset_m − offset_seed)   averaged over the seed-module mix
```

and only becomes "offset residual after calibration" if the calibration
has been baked into the TSV, which **is not the case today**.

Downstream implications the script must handle:

1. **Default mode: report raw μ.** No pretense that μ is an offset residual.
2. **`--time-calib <path.json>` mode (optional):** load the same JSON schema
   used by `HyCalTimeCalib` (keys `default`, `modules[].name` /
   `modules[].module_id` / `modules[].offset_ns`; see
   `prad2det/include/HyCalTimeCalib.h:37–53` and
   `database/hycal_time_offsets/time_calib.json`) and apply
   `dt_corr = dt − (offset_neighbour − offset_seed)` **in Python, at
   analysis time**, before histogramming. In this mode μ becomes the
   post-calibration residual — the quantity that should be ~0 for a
   well-calibrated module.
3. Always emit the mode (raw vs. offset-corrected, with the JSON path and
   SHA-1 of the file for provenance) into the summary CSV header comment
   and the markdown report.

### 2.a.1 Confirmation (researcher trace)

The "raw / uncorrected `dt`" finding is confirmed end-to-end by
`docs/research-study-hycal-timing-time-offset.md`, which traced
`pk.time` through the following path and found **no subtraction of
`Module::time_offset` anywhere on it**:

```
WaveAnalyzer::Analyze
  → bind_dec.cpp `analyze` lambda
  → study_hycal_timing.py add_hit loop
  → HyCalCluster::AddHit
  → HyCalCluster::CollectNeighborTiming
  → SeedNeighborTiming pybind11 binding
  → _common.open_table_writer TSV writer
```

Three mechanical facts from that trace are worth citing explicitly so
nobody re-opens the question:

- **Sign convention.**  Calibrated time = raw time − `Module::time_offset`.
  Documented at `prad2det/include/HyCalTimeCalib.h:8` and applied in
  the one C++ consumer that does use it,
  `prad2dec/src/WaveAnalyzer.cpp:678–684` (the subtraction line).  Our
  Python `apply_offsets` must therefore compute
  `dt_corr = dt_raw − (offset_neighbour − offset_seed)` — same sign,
  propagated through the (neighbour − seed) difference.
- **`WaveAnalyzer::Analyze` offset parameter is always zero in the
  Python path.**  The C++ signature accepts an optional
  `time_offset = 0.f`, but the pybind11 binding at
  `python/bind_dec.cpp:407–422` never forwards one, so every peak
  returned to Python has `pk.time` relative to the raw sample grid.
  No hidden offset gets applied inside `WaveAnalyzer` for us.
- **`fdec::Module::time_offset` is not exposed to Python.**  The
  `det.HyCalSystem` / `det.Module` binding at
  `python/bind_det.cpp:531–572` does not surface `time_offset` as a
  readable attribute.  **This is why `load_time_calib` must parse the
  JSON file directly instead of reading the per-module offset off the
  bound C++ object.**  Do not propose "just read `Module.time_offset`
  from the binding" in a future revision — it isn't there, and
  exposing it would be a C++/binding change, which is out of scope
  (see §3).

## 3. Non-goals / explicitly out of scope

- Any C++ change (`HyCalCluster`, `HyCalTimeCalib`, `PipelineBuilder`, etc.).
- Any change to `study_hycal_timing.py`.
- Re-running the EVIO decode / clusterer from this script. We consume the TSV.
- Singles-rate accidentals model `R_seed · R_module · 2Δt`. The script
  reports the data-driven background `B` from the fit and (optionally) a
  sideband-window cross-check; converting these to an expected rate from
  independent singles rates is a separable follow-up.
- Writing a new calibration JSON. If this diagnostic surfaces per-module
  offsets, producing a new `hycal_time_offsets/*.json` is a different task.

## 4. Location and naming

**Target path:** `analysis/pyscripts/module_dt_diagnostic.py`

Sibling convention (verified against the `analysis/pyscripts/` listing):

- `benchmark_*.py` and `study_*.py` run the pipeline; `plot_*.py` consume
  outputs. `module_dt_diagnostic.py` is a **consume-and-plot** script, so
  `plot_module_dt.py` would also fit — but "diagnostic" better advertises
  the fitting + flagging payload. Keep `module_dt_diagnostic.py`.
- No `_common.py` dependency is strictly required (we don't need
  `setup_pipeline`), but we WILL import `_common` for `open_table_writer`
  (`_common.py:417`) if we need to emit a secondary table and so that
  future CLI-arg sharing is cheap. Final call: **do not import `_common`**;
  this script is deliberately pipeline-free so it can run on any machine
  that has the TSV without needing the C++ build. Standalone imports only.

## 5. Dependencies

Baseline project deps (from `scripts/requirements.txt` and
`calibration/requirements.txt`): `matplotlib`, `numpy`. The repo does not
pin a `requirements.txt` for `analysis/pyscripts/` and does not ship a
`pyproject.toml` / `environment.yml`.

- **Must use:** `numpy`, `matplotlib` (incl. `matplotlib.backends.backend_pdf.PdfPages`).
- **Soft-require `scipy`:** needed for `scipy.optimize.curve_fit` (same
  Gaussian fit pattern as `benchmark_hycal_timing.py:44–48, 225–234`).
  Follow that script's `HAVE_SCIPY` try/except idiom:
  - If scipy is present, do the fits.
  - If absent, run in a **"plots-only, no fits" degraded mode**: still
    produce the per-module histograms, the 2-D heatmap, the sideband
    counts, and a summary CSV in which every row has
    `fit_status = "scipy-missing"` and NaN fit parameters.  Emit a
    single clear warning line to stderr
    (`"[WARN] scipy not available — running in plots-only mode"`) and do
    not crash.  The markdown report must call out the degraded mode in
    its metadata block.
  - **Where the warning is emitted.**  The `HAVE_SCIPY` try/except at
    module scope sets the flag silently; the warning line is printed
    from the start of `main()`, not at import time.  This keeps
    `python -c "import module_dt_diagnostic; module_dt_diagnostic._self_test()"`
    and any future reuse of the module as a library from spamming stderr
    just because the user imported it.  `--self-test` also calls `main`
    and therefore prints the warning once if scipy is missing — correct
    behaviour, since the synthetic Gaussian case in §11.a cannot pass
    without scipy.
- **No pandas.** Even though two sibling scripts mention pandas in
  comments, no script in `analysis/pyscripts/` actually imports it. Keep
  the dep surface flat. Use the stdlib `csv` module for reading and
  writing.

No change to any requirements file is needed. The script's module
docstring should list the deps explicitly so users on a slim env know
what to `pip install`.

## 6. CLI

Modelled on `C.add_common_args` (`_common.py:364–392`) in style but
with its own positional set because this script does not need an EVIO
path. Use `argparse` with `RawDescriptionHelpFormatter`.

```
module_dt_diagnostic.py INPUT_TSV OUTPUT_DIR [options]

positional:
  INPUT_TSV                 TSV/CSV from study_hycal_timing.py
  OUTPUT_DIR                directory to write plots + summary into
                            (created if missing)

input shape:
  --csv                     parse INPUT_TSV as CSV (default: TSV, matching
                            study_hycal_timing.py's default)
  --no-header               TSV/CSV has no header row; use the fixed column
                            order from §2 (mirrors _common.add_common_args)

grouping / filtering:
  --module-id ID [ID ...]   restrict to these neighbour PrimEx IDs
                            (default: all IDs that appear in the TSV)
  --min-seed-energy MEV     drop rows with seed_energy < MEV
                            (default 0 — accept all seeds already selected
                            upstream by min_center_energy)
  --max-dist-q F            drop rows with dist_q > F (default inf;
                            upstream already caps at --max-qdist)
  --min-entries N           modules with fewer than N rows are flagged
                            "low-stats" and skipped from fitting
                            (default 200)

ΔT windowing:
  --dt-window LO HI         histogram range in ns (default -40 40)
  --dt-bin F                bin width in ns (default 0.5)
  --fit-window F            half-width (ns) around the histogram argmax
                            used for the initial Gaussian seed window
                            (default 8.0)
  --sideband LO HI          off-time sideband range in ns for the
                            background cross-check (default -40 -20)

time calibration (optional, §2.a):
  --time-calib PATH         apply per-module offsets from this JSON
                            (same schema as HyCalTimeCalib)
                            Default: none → μ is raw (uncorrected).

flagging thresholds:
  --sig-threshold F         flag "no peak" if significance < F
                            (default 5.0)
  --mu-threshold F          flag "offset peak" if |μ| > F ns
                            (default 2.0)
  --sigma-threshold F       flag "wide peak" if σ > F ns (default 3.0)

streaming / performance (see §7.c):
  --chunk-rows N            rows per read chunk in the streaming loader
                            (default 1_000_000).  Performance knob only;
                            does not affect numerical output.
  --progress-every N        print a progress line every N million rows
                            read (default 10).  The line reports rows
                            read, rows kept after filters, elapsed wall
                            time, and ETA extrapolated from bytes
                            consumed vs. file size.  Set to 0 to silence.
  --max-rows N              hard cap on total rows read from the input
                            TSV (default 0 = no cap).  Useful for smoke
                            tests on huge TSVs; see §11.b validation.

output control:
  --run-id STR              appended to output filenames; default = sniffed
                            from the input filename ("run_<N>" if a
                            prad_<N>.* substring is found, else timestamp)
  --tag STR                 extra tag appended after run-id (default "")
  --no-pdf                  skip the per-module PDF (useful for CI runs)
  -v, --verbose             log per-module fit outcomes to stderr
```

Validation of CLI args happens in a dedicated `_validate_args(ns)` function
raising `SystemExit(2)` on bad combinations (e.g. `--sideband` inside
`--dt-window ± fit-window`).

**Scipy-missing degraded mode (see §5).**  When `scipy` is unavailable the
script runs with every fit-dependent threshold rendered meaningless: all
rows carry `fit_status = "scipy-missing"` and `--sig-threshold`,
`--mu-threshold`, `--sigma-threshold` collapse to no-ops.  Only
`--min-entries` still produces meaningful flags (`low-stats`).  Plots and
sideband counts are still produced; the warning emitted at startup is
the one source of truth for the implementer and the user.

## 7. Script structure (single file)

### 7.a Fit convention (READ FIRST)

The script fits a **density** model, not a counts-per-bin model.  This
single convention is used everywhere — histogramming, fitting,
significance, sideband check, summary CSV, per-module PDF — so that `A`,
`B`, `S`, and the sideband background density all have consistent units.

- Histogram ΔT with bin width `dt_bin` to produce raw per-bin counts
  `n_i`.
- Convert to **density** before fitting:
  `y_i = n_i / dt_bin`             (units: counts/ns)
  `σ_y_i = max(sqrt(n_i), 1) / dt_bin`  (floor the 0-count bin error at
  1 count so `curve_fit` has a non-zero weight everywhere).
- Fit `f(x) = A · exp(-½((x-μ)/σ)²) + B` to `(x_i, y_i, σ_y_i)` with
  `sigma=σ_y_i, absolute_sigma=True`.  `A` and `B` are now in counts/ns.
- Signal integral in ±3σ:
  `S = A · σ · √(2π) · erf(3/√2)`   (counts; closed form)
- Background counts in the same ±3σ window:
  `B_w = B · 6σ`   (counts)
- Significance: `S / sqrt(max(B_w, 1))`.
- Sideband density: `bg_per_ns = n_sideband / (sideband_hi − sideband_lo)`
  (counts/ns, directly comparable to the fit's `B`).
- Summary CSV column names reflect the units: `A_per_ns`, `B_per_ns`,
  `sideband_bg_per_ns`, `signal_counts`.

### 7.b Module layout

```
module_dt_diagnostic.py
├── module docstring (deps, usage, output artifacts)
├── imports (stdlib: argparse, csv, json, math, hashlib, sys, os, time,
│            dataclasses, pathlib, datetime, typing;
│            third-party: numpy, matplotlib, matplotlib.backends.backend_pdf.PdfPages;
│            optional: scipy.optimize.curve_fit via HAVE_SCIPY flag)
├── constants + dataclasses
│   ├── TSV_COLUMNS       (list[str], frozen schema from §2)
│   ├── @dataclass ModuleStats
│   │     Fields are the authoritative union of §8.a (summary-CSV
│   │     schema), §9.b (fit-status / primary-fit fields), and §14.d
│   │     (secondary-peak extension).  Do NOT re-enumerate here — all
│   │     downstream readers (write_summary_csv, write_markdown_report,
│   │     plot_per_module_pdf, plot_hycal_spatial_map) read §8.a /
│   │     §8.b / §14.d directly, and keeping two field lists in sync
│   │     has already caused one iteration's worth of drift.
│   └── @dataclass RunContext
│         input_path, output_dir, run_id, tag,
│         dt_lo, dt_hi, dt_bin,
│         fit_window, sideband_lo, sideband_hi,
│         sig_thr, mu_thr, sigma_thr, min_entries,
│         module_id_filter: set[int] | None,     # from --module-id
│         min_seed_energy: float,                 # from --min-seed-energy
│         max_dist_q: float,                      # from --max-dist-q
│         time_calib_path, time_calib_sha1,
│         offsets: dict[int, float],              # keyed by PrimEx id
│         default_offset: float,
│         scipy_available: bool,
│         chunk_rows: int,                        # from --chunk-rows
│         progress_every: int,                    # from --progress-every
│         max_rows: int,                          # from --max-rows (0 = no cap)
│         dt_edges: np.ndarray,                   # derived: np.arange(dt_lo,
│                                                 #   dt_hi + dt_bin, dt_bin);
│                                                 #   shared across every
│                                                 #   accumulator so fit_module
│                                                 #   and stream_tsv agree
│         sideband_edges: np.ndarray              # derived: analogous over
│                                                 #   [sideband_lo, sideband_hi]
├── load_time_calib(path) -> (offsets, default, sha1)
│     - parse JSON matching HyCalTimeCalib.h:37–53 schema
│     - resolution:  modules[].module_id (preferred)  OR  modules[].name
│       (strip leading G/W, map to PrimEx id: G→+0 offset in [1,576],
│        W→+1000 offset in [1001,2152]).  Pure helper, exercised by the
│        --self-test block (see §11).
│     - "default" value applied in the apply step to any id not present.
│     - compute SHA-1 of the file bytes for provenance.
├── @dataclass ModuleAccumulator            # see §7.c
│     module_id: int
│     dt_hist:       np.ndarray[int64]  # length = N_dt_bins
│     sideband_hist: np.ndarray[int64]  # length = N_sideband_bins
│     n_entries:     int                # total post-filter rows that
│                                       # fell inside [dt_lo, dt_hi]
│     n_total_seen:  int                # rows (post-filter) for this
│                                       # module BEFORE the dt-window
│                                       # cut, for diagnostics
├── stream_tsv(path, csv_mode, no_header, ctx) -> dict[int, ModuleAccumulator]
│     Streaming replacement for the old load_tsv + filter_rows +
│     apply_offsets + group_by_module stack.  See §7.c for the full
│     memory/throughput model.  High-level loop:
│
│       1. Open path with the stdlib `csv.reader` (TSV delimiter by
│          default; CSV when --csv).  If --no-header is NOT set, read
│          one header line and validate it against REQUIRED_COLUMNS;
│          build column-index lookups.  If --no-header IS set, use the
│          fixed TSV_COLUMNS order.
│       2. Pre-compute edges arrays once:
│            dt_edges       = np.arange(ctx.dt_lo,      ctx.dt_hi      + ctx.dt_bin, ctx.dt_bin)
│            sideband_edges = np.arange(ctx.sideband_lo, ctx.sideband_hi + ctx.dt_bin, ctx.dt_bin)
│       3. Pre-size per-module accumulators lazily (dict keyed by
│          neighbor_id; create on first sight).  Each accumulator holds
│          two int64 histograms + two counters.
│       4. Buffer rows into a chunk of size ctx.chunk_rows using plain
│          Python lists, then vectorise:
│            cols = {
│              "seed_id":           np.fromiter(..., dtype=np.int32),
│              "neighbor_id":       np.fromiter(..., dtype=np.int32),
│              "dt":                np.fromiter(..., dtype=np.float32),
│              "seed_energy":       np.fromiter(..., dtype=np.float32),
│              "dist_q":            np.fromiter(..., dtype=np.float32),
│            }
│          (seed_time, neighbor_time, neighbor_energy, event_num,
│          trigger_bits, seed_module, neighbor_module, dx_q, dy_q are
│          read but discarded — only the columns in REQUIRED_COLUMNS
│          minus the ones we don't need are kept.  Save on allocation.)
│       5. Apply the four row filters as per-chunk numpy masks, in the
│          SAME committed order as the old filter_rows step:
│            a. seed_energy >= ctx.min_seed_energy
│            b. dist_q      <= ctx.max_dist_q
│            c. neighbor_id ∈ ctx.module_id_filter  (if not None)
│            d. ctx.dt_lo <= dt_raw <= ctx.dt_hi    (on RAW dt, before offsets)
│       6. apply_offsets vectorised on the surviving chunk (see below).
│          dt_corr replaces the dt array for histogramming only; the
│          raw dt is no longer needed after step 5(d).
│       7. For every unique neighbor_id in the chunk, slice the chunk,
│          np.histogram the chunk's dt_corr slice against dt_edges and
│          sideband_edges, and add the resulting count vectors into the
│          accumulator (`acc.dt_hist += ...`, `acc.sideband_hist +=
│          ...`).  Update n_entries and n_total_seen.
│       8. Every ctx.progress_every × 1_000_000 rows read (when > 0),
│          print the progress line described in §6 to stderr.
│       9. Stop when EOF, when ctx.max_rows > 0 and that cap is reached,
│          or when a hard upper bound on module ids (unexpected garbage
│          id) triggers a loud error.
│      10. Return the dict of accumulators.
│
│     Per-chunk row counts before/after each mask are logged under -v.
├── apply_offsets_chunk(seed_ids, neighbor_ids, dt, offsets, default) -> np.ndarray
│     Vectorised chunk-level offset correction:
│       off_s = np.array([offsets.get(int(s), default) for s in seed_ids],
│                        dtype=np.float32)
│       off_n = np.array([offsets.get(int(n), default) for n in neighbor_ids],
│                        dtype=np.float32)
│       return dt - (off_n - off_s)
│     If ctx.offsets is empty, callers skip this and pass dt through
│     unchanged.  Dict lookups per row are acceptable at chunk_rows =
│     1M (sub-second per chunk on CPython); if profiling later shows
│     this is a bottleneck, replace with a precomputed np.ndarray[int32
│     -> float32] lookup table sized to max(primex_id)+1.
│     (NOTE: offset application happens inside stream_tsv step 6; this
│     helper exists so --self-test can exercise it directly, see §11.a.)
├── gaussian_plus_const(x, A, mu, sigma, B)
│     A * exp(-0.5*((x-mu)/sigma)**2) + B           # density
├── fit_module(acc, ctx) -> ModuleStats
│     Inputs: a ModuleAccumulator (contains dt_hist, n_entries) and the
│     RunContext.  The accumulator already carries the binned histogram,
│     so the fit skips the histogram step and goes straight to the
│     density conversion in §9.a step 3.  See §9 for the full procedure,
│     including:
│       - density conversion before curve_fit (y = dt_hist / dt_bin),
│       - the seeded A0 = max(0, max(y) - B0) initial amplitude,
│       - curve_fit(..., sigma=σ_y, absolute_sigma=True, maxfev=5000),
│       - the five possible fit_status values,
│       - the at-bound tolerance `max(0.5 · dt_bin, 1e-6)` applied to
│         both μ and σ (bin-sized, NOT a fraction of the bound span —
│         see §9.a step 7).
├── sideband_check(acc, ctx, mu0) -> (counts, bg_per_ns, overlap)
│     Inputs: ModuleAccumulator (carries sideband_hist, pre-binned over
│     [sideband_lo, sideband_hi] during the streaming pass) + mu0.
│     - if [ctx.sideband_lo, ctx.sideband_hi] overlaps
│       [mu0 - ctx.fit_window, mu0 + ctx.fit_window]:
│         return (NaN, NaN, overlap=True)
│         (caller appends "sideband-overlap" to the flag list; see §10.a)
│     - else counts = acc.sideband_hist.sum()
│     - divide by window width in ns → bg_per_ns
│     - this is a background DENSITY (counts/ns), directly comparable to
│       the fit's B_per_ns.  Comparison lives in the summary table.
│     - `mu0` is the fitted μ when fit_status == "ok", otherwise the
│       bin-centre argmax of acc.dt_hist (so overlap is still checked
│       in degraded modes).
├── flag_module(stats, ctx) -> list[str]
│     Flags, in fixed order (a module can collect several).  Full
│     vocabulary per §9.b (authoritative — currently six flags
│     including the §14 `secondary-peak`; do NOT re-enumerate here).
│       "low-stats"       if n_entries < min_entries   (n_entries is
│                                                       the post-filter,
│                                                       post-dt-window
│                                                       count accumulated
│                                                       by stream_tsv —
│                                                       see §7.b, §7.c)
│       "no-peak"         if fit_status != "ok" OR significance < sig_thr
│       "offset-peak"     if fit_status == "ok" and |mu|  > mu_thr
│       "wide-peak"       if fit_status == "ok" and sigma > sigma_thr
│       "sideband-overlap" appended by sideband_check (see §10.a);
│                          informational — does not imply "no-peak".
├── plot_heatmap(accumulators, ctx, path_stem)
│     - ids = sorted(accumulators.keys())
│     - counts2d[iy, :] = accumulators[id].dt_hist          (int64 → float)
│     - imshow with aspect='auto', norm=LogNorm(vmin=1), origin='lower',
│       extent=[dt_lo, dt_hi, 0, N_modules]; label y-axis with the actual
│       PrimEx ids (yticks at midpoints of visible bands; thin out if
│       >50 modules so labels remain legible — e.g. every 25th tick).
│     - vertical line at μ=0 for reference
│     - colorbar, title "Per-module ΔT (run <id>)", saved as PNG + PDF
│     - Memory: pure restacking of the per-module histograms, no extra
│       passes over the TSV.
├── plot_per_module_pdf(accumulators, stats_by_id, ctx, path)
│     - PdfPages, one page per module id (sorted)
│     - top: density histogram (acc.dt_hist / dt_bin; counts/ns) + fit
│       curve overlay + vertical lines at μ, μ±3σ
│     - bottom (text box): N, μ ± dμ, σ ± dσ, B_per_ns, S, significance,
│       flags, fit_status
│     - modules with fit_status != "ok" still get a page with the
│       histogram and a "NO FIT (<status>)" annotation
├── write_summary_csv(stats_list, ctx, path)
│     Writes columns per §8.a (authoritative schema).  Comment header
│     line(s) (starting with #) carry: input path, time-calib mode +
│     sha1, dt window, fit window, sideband window, thresholds, filter
│     settings, scipy availability, timestamp.
├── write_markdown_report(stats_list, ctx, path)
│     - run metadata block (incl. scipy-missing warning when applicable)
│     - counts by flag category
│     - per-flag section listing flagged modules; writes per-module
│       lines per §8.b format (authoritative), sorted so the worst
│       offenders are at the top
│     - closing note on raw-vs-calibrated μ interpretation (§2.a)
├── sniff_run_id(input_path) -> str
│     reuse the regex from _common._RUN_PAT style (local copy; do not import)
├── _self_test() -> int                             # see §11, --self-test
└── main(argv=None) -> int
      1. parse args, validate (and handle --self-test early return)
      2. ensure OUTPUT_DIR exists (Path.mkdir(parents=True, exist_ok=True))
      3. if --time-calib given: load_time_calib → fill ctx.offsets +
         ctx.default_offset + ctx.time_calib_sha1
      4. stream_tsv(path, ..., ctx) -> dict[int, ModuleAccumulator]
         This is the ONLY pass over the TSV.  It performs, in order:
         read-chunk → filter (seed_energy → dist_q → module_id →
         dt-window on RAW dt) → apply_offsets_chunk (no-op if empty) →
         per-id np.histogram into dt_hist and sideband_hist → accumulate.
      5. for each accumulator: fit_module (uses acc.dt_hist directly),
         sideband_check (uses acc.sideband_hist directly), flag_module.
      6. write artifacts in a stable order:
         <OUTPUT_DIR>/
           module_dt_summary_<run>[_<tag>].csv
           module_dt_heatmap_<run>[_<tag>].png
           module_dt_heatmap_<run>[_<tag>].pdf
           module_dt_fits_<run>[_<tag>].pdf       (skipped if --no-pdf)
           module_dt_flags_<run>[_<tag>].md
      7. print one-line summary to stdout + return 0 (or non-zero if
         zero rows survived filtering — a clear failure signal for CI)
```

### 7.c Streaming / memory model

The production target is a single-run TSV at ~324 M rows, ~25–40 GB on
disk (confirmed on run 25320).  The script must never materialise the
whole TSV in memory; it does one streaming pass and accumulates
pre-binned per-module histograms, which are all the downstream stages
need.

**Why histograms are sufficient.**  Every downstream consumer (the
Gaussian + flat-background fit in §9, the 2-D heatmap, the sideband
density check in §10) is a function of the per-module ΔT histogram
alone.  Raw rows carry no information the histogram does not.  Once
the TSV has been streamed once, the fit, heatmap, PDF, and summary
CSV are all O(N_modules), not O(N_rows).

**Memory budget (back-of-envelope).**

| component                   | size                                                                      |
|-----------------------------|---------------------------------------------------------------------------|
| per-module `dt_hist`         | `N_modules × N_dt_bins × 8 B` = 2100 × 160 × 8 ≈ 2.7 MB                   |
| per-module `sideband_hist`   | `N_modules × N_sideband_bins × 8 B` = 2100 × 40 × 8 ≈ 0.7 MB              |
| per-module counters          | `N_modules × 2 × 8 B` ≈ 34 KB                                             |
| chunk row buffer (vectorised)| `chunk_rows × ~5 kept cols × 4–8 B` ≈ 20–40 MB at `--chunk-rows 1_000_000`|
| offsets dict + lookup helper | O(N_modules) ≈ 100 KB                                                     |

Total steady-state footprint is **well under 100 MB** at the default
`--chunk-rows`.  The §11.b validation step pins a 500 MB soft ceiling
as a sanity bound; anything approaching that indicates a leak or a
chunk-buffer bug.

**Throughput and progress reporting.**  The streaming loader emits a
progress line every `--progress-every` million rows (default 10) with:
rows read, rows kept after filters, elapsed wall time, and an ETA
extrapolated from `(bytes_consumed / total_bytes) × elapsed`.  Set
`--progress-every 0` to silence (useful in CI).

**Correctness invariants the streaming loader must preserve.**

- Filter order remains `seed_energy → dist_q → module_id → dt_window`
  (committed in §7.b), applied per chunk.
- `--dt-window` is still applied to **raw** `dt`, before offsets
  (committed in §7.b with the intent note; same reasoning applies in
  the streaming path — no change in semantics).
- **`n_entries` semantics (pinned).**  `n_entries = int(acc.dt_hist.sum())`
  — the integer sum of the per-module ΔT histogram.  This explicitly
  excludes rows whose offset-corrected `dt_corr` fell outside
  `[ctx.dt_lo, ctx.dt_hi]` after `apply_offsets_chunk`, because such
  rows do not land in any bin of `acc.dt_hist`.  `stream_tsv` maintains
  `acc.n_entries` as a running integer so downstream code never has to
  re-sum, but the definition is the histogram sum.  This is what the
  fit sees, so it is the right denominator for `--min-entries` and for
  significance interpretation.
- `n_total_seen` is a diagnostic-only counter (post-filter, pre-dt-window
  — i.e. after `seed_energy → dist_q → module_id` masks but before the
  raw-dt-window mask).  The difference `n_total_seen − n_entries` is
  the "lost to raw-dt-window and/or offset-induced range exits" count.
  Large values indicate either (a) a mis-sized `--dt-window` for the
  run, or (b) a wrong `--time-calib` JSON that is pushing data out of
  the window — useful signal when diagnosing a bad offset table.
- `apply_offsets_chunk` is a pure function of (seed_id, neighbor_id,
  dt, offsets, default) and has no cross-chunk state, so the chunk
  size never affects numerical output — only throughput.

## 8. Output artifacts (summary)

| File | Format | Content |
|---|---|---|
| `module_dt_summary_<run>.csv` | CSV | one row per neighbour module id with μ, σ, B, S, significance, fit_status, flags, plus `n_entries` / `n_total_seen` for data-exclusion diagnostics |
| `module_dt_heatmap_<run>.{png,pdf}` | raster + vector | 2-D (module_id vs ΔT) heatmap — headline plot |
| `module_dt_fits_<run>.pdf` | multi-page PDF | per-module ΔT histogram + fit overlay (one page per module) |
| `module_dt_flags_<run>.md` | Markdown | flagged modules grouped by flag reason, sorted worst-first; each line shows `n_entries / n_total_seen` so a reader can see how much data was clipped |

All four share a `<run>[_<tag>]` suffix so multiple diagnostic runs can
coexist in the same directory.

### 8.a Summary CSV column list

In this order (ties `n_entries` and `n_total_seen` next to each other
so a reader can eyeball the clipped fraction):

| column                | meaning                                                                                                 |
|-----------------------|---------------------------------------------------------------------------------------------------------|
| `module_id`           | PrimEx neighbour id (G: 1–576, W: 1001–2152).                                                           |
| `n_entries`           | Rows in `acc.dt_hist` — post-filter, post-raw-dt-window, **post-offset-correction** (what the fit saw). |
| `n_total_seen`        | Rows passing seed_energy + dist_q + module_id filters, **before** the dt-window mask; `n_entries` is the subset that landed in `dt_hist` after offsets.  Difference `n_total_seen − n_entries` is the "lost to dt-window / offset" count (§7.c). |
| `mu_ns`, `mu_err_ns`  | Fitted peak position and 1-σ error (NaN if `fit_status != "ok"`).                                       |
| `sigma_ns`, `sigma_err_ns` | Fitted peak width and 1-σ error.                                                                   |
| `A_per_ns`            | Fitted Gaussian density amplitude (counts/ns) — see §7.a fit convention.                                |
| `B_per_ns`            | Fitted flat-background density (counts/ns).                                                             |
| `signal_counts`       | `S = A · σ · √(2π) · erf(3/√2)` — closed-form ±3σ signal integral.                                     |
| `significance`        | `S / √max(B·6σ, 1)`.                                                                                    |
| `fit_status`          | One of {`ok`, `low-stats`, `scipy-missing`, `fit-failed`, `at-bound`} per §9.b.                         |
| `sideband_bg_per_ns`  | Density in the configured sideband (counts/ns); NaN if `sideband-overlap` was flagged (§10.a).          |
| `sideband_counts`     | Integer entries in the sideband; NaN under `sideband-overlap`.                                          |
| `secondary_mu`        | Fitted μ₂ of the residual-detected secondary peak in ns (§14.a); NaN when no candidate or detection skipped. |
| `secondary_sigma`     | Fitted σ₂ of the secondary peak in ns; NaN otherwise.                                                   |
| `secondary_amp_per_ns`| Fitted amplitude of the secondary Gaussian (counts/ns); NaN otherwise.                                  |
| `secondary_significance` | `S₂ / √max(B·6σ₂, 1)` where `B` is the primary fit's background density (§14.a); NaN otherwise.      |
| `flags`               | Semicolon-separated subset of {`low-stats`, `no-peak`, `offset-peak`, `wide-peak`, `sideband-overlap`, `secondary-peak`}. |

The CSV comment header (lines starting with `#`) carries input path,
time-calib mode + sha1, dt window, fit window, sideband window, all
thresholds (incl. `--secondary-sig-threshold` and
`--secondary-search-window`), filter settings, scipy availability,
hycal-geometry file path when loaded, and generation timestamp.

### 8.b Markdown report per-module line format

Each flagged module appears under its flag heading as one line:

```
- <module_name> (id <module_id>)  n=<n_entries>/<n_total_seen>
  μ=<mu_ns> ± <mu_err_ns> ns  σ=<sigma_ns> ± <sigma_err_ns> ns
  sig=<significance>  [flags]
```

`n_entries/n_total_seen` makes it visually obvious when a flag is
driven by heavy clipping rather than by the fit itself.
`<module_name>` is the PrimEx `G<n>` / `W<n>` string resolved from the
PrimEx id (independent of the hycal-geometry file — it is a pure
`G/W + numeric suffix` mapping; see §14.b).  If the geometry file is
loaded the sector tag (`Top`/`Bottom`/`Left`/`Right`/`Center`) is
appended in parentheses after the name.

**For modules flagged `secondary-peak`**, append a second indented
line with the secondary-fit summary:

```
    secondary μ=<secondary_mu> ± <secondary_mu_err> ns
      σ=<secondary_sigma> ± <secondary_sigma_err> ns
      sig=<secondary_significance>
```

**Top-of-report section** (added before the per-flag sections): see
§14.c "Modules with secondary peak" — a dedicated roll-up section
listing every `secondary-peak`-flagged module as one line sorted by
`secondary_significance` descending.  This section is the direct
answer to the "print out the module name" requirement.

## 9. Fitting approach — mathematical details

Model (density; see §7.a): `f(x) = A · exp(-½((x-μ)/σ)²) + B`
where `A`, `B` are in counts/ns, consistent with `y_i = n_i / dt_bin`.

Rationale:

- A single Gaussian + flat background is the minimum model that separates
  "real coincidence peak" from "accidentals pedestal". It is also what
  `benchmark_hycal_timing.py` uses for its spectrum fits (same codebase
  precedent — reuse that style).
- Fit on the local window `[μ₀ − fit_window, μ₀ + fit_window]` where μ₀
  is the bin-centre argmax of the full-range density histogram. This
  avoids pulling σ toward the tails and matches
  `benchmark_hycal_timing.py:220–222`.
- Integrated signal `S = A · σ · √(2π) · erf(3/√2) ≈ A · σ · 2.5066 · 0.9973`
  (counts) — closed form, no numerical integration.
- Background in the ±3σ signal window: `B_w = B · 6σ` (counts).
- Significance: `S / √max(B_w, 1)` (floor to avoid div-by-zero when
  background is empty).

### 9.a Procedure inside `fit_module(acc: ModuleAccumulator, ctx: RunContext) -> ModuleStats`

The accumulator is pre-populated by `stream_tsv` (§7.b, §7.c):
`acc.dt_hist` and `acc.sideband_hist` are already the correct
pre-binned histograms over `ctx.dt_edges` and `ctx.sideband_edges`.
`fit_module` never sees raw ΔT rows and never calls `np.histogram`.

The number of entries that fed into the fit is
`n_entries = int(acc.dt_hist.sum())` — identical to the field
`acc.n_entries` maintained by `stream_tsv`.  Use `acc.n_entries`
everywhere for consistency; the sum form is only mentioned here as the
definition.  See the §7.c "n_entries semantics (pinned)" note for why
this count excludes rows whose offset-corrected `dt_corr` fell outside
`[ctx.dt_lo, ctx.dt_hi]`.

1. **Low-stats short-circuit.** If `acc.n_entries < ctx.min_entries`,
   return a `ModuleStats` with `fit_status = "low-stats"`, NaN fit
   parameters, `n_entries = acc.n_entries`, and
   `n_total_seen = acc.n_total_seen`.  Skip to the sideband step and
   return.
2. **Scipy-missing short-circuit.** If `not HAVE_SCIPY`, return
   `ModuleStats` with `fit_status = "scipy-missing"`, NaN fit
   parameters.  The plot layer reads `acc.dt_hist` directly, so no
   plot is lost.
3. **Density view of the histogram.**
   - `counts  = acc.dt_hist`                              # int64, pre-binned
   - `centres = 0.5 * (ctx.dt_edges[:-1] + ctx.dt_edges[1:])`
   - `y       = counts / ctx.dt_bin`                      # counts/ns
   - `y_err   = np.maximum(np.sqrt(counts), 1.0) / ctx.dt_bin`
4. **Seed parameters.**
   - `mu0 = centres[np.argmax(counts)]`
   - `B0 = median(y[ outer-fifth bins ])`   (bins outside the central
     60% window — these are the tails that approximate the pedestal)
   - `A0 = max(0.0, y.max() - B0)`           (do NOT double-count B0)
   - `sigma0 = 1.5`  ns
5. **Fit window mask.**
   `mask = (centres >= mu0 - ctx.fit_window) & (centres <= mu0 + ctx.fit_window)`
   — require at least 5 bins with `counts > 0` inside the mask; else
   return `fit_status = "low-stats"` (bin-level low stats inside the
   fit window specifically).
6. **curve_fit call.**
   ```
   bounds = (
       [0.0,             ctx.dt_lo,    0.1,              0.0],
       [np.inf,           ctx.dt_hi,    ctx.fit_window,   np.inf],
   )
   popt, pcov = curve_fit(
       gaussian_plus_const,
       centres[mask], y[mask],
       p0 = [A0, mu0, sigma0, B0],
       sigma = y_err[mask], absolute_sigma = True,
       bounds = bounds, maxfev = 5000,
   )
   ```
   - On `RuntimeError` or `OptimizeWarning`-promoted exception →
     `fit_status = "fit-failed"`, NaN params, return.
   - `mu_err, sigma_err = sqrt(pcov[1,1]), sqrt(pcov[2,2])` with NaN
     if `pcov` contains inf (singular covariance).
7. **At-bound detection.**  Both tolerances are tied to the bin size
   (small absolute values), not to the bound spans — otherwise a σ
   bound of `[0.1, fit_window=8.0]` gives a half-span tolerance of ~4 ns
   that would mark every realistic HyCal timing width (σ ≈ 1 ns) as
   `at-bound`.  Use:
   - `tol_mu    = max(0.5 * ctx.dt_bin, 1e-6)`
   - `tol_sigma = max(0.5 * ctx.dt_bin, 1e-6)`

   Then:
   - `|mu - ctx.dt_lo| < tol_mu` or `|mu - ctx.dt_hi| < tol_mu`
     → `at-bound`.
   - `|sigma - 0.1| < tol_sigma` or `|ctx.fit_window - sigma| < tol_sigma`
     → `at-bound`.
   - Keep the fitted params in `ModuleStats` so the user can inspect
     them; the flag is informational.
8. **Derived quantities.**
   - `S = A * sigma * sqrt(2π) * erf(3/√2)`
   - `B_w = B * 6 * sigma`
   - `significance = S / sqrt(max(B_w, 1.0))`
9. **Return.**  `fit_status = "ok"` unless an earlier branch set
   something else.

### 9.b Vocabulary — the full set of `fit_status` values

Exactly these five, no others:

| value          | meaning                                                      |
|----------------|--------------------------------------------------------------|
| `ok`           | fit converged, μ and σ inside bounds (within tolerance).     |
| `low-stats`    | `n_entries < min_entries` OR too few non-empty bins in the fit window. |
| `scipy-missing`| scipy not importable; fit was never attempted (see §5).      |
| `fit-failed`   | `curve_fit` raised, or returned a singular covariance.       |
| `at-bound`     | fit converged but μ or σ sits at a bound within tolerance.   |

(The earlier draft's `unresolved` was redundant with `at-bound` applied
to σ; it has been dropped and merged into `at-bound`.)

In every non-`ok` case the sideband step still runs so the summary CSV
and the per-module PDF page both carry something useful.

**Flag vocabulary addendum.**  The §14 extension introduces one new
informational flag, `secondary-peak`, into the fixed flag vocabulary
(§7.b `flag_module`).  It is appended by the secondary-peak detection
step in §14.a; it never changes `fit_status`.  Full flag set is now
{`low-stats`, `no-peak`, `offset-peak`, `wide-peak`, `sideband-overlap`,
`secondary-peak`}.

## 10. Sideband cross-check

Pick a window far from the expected peak (default `-40` to `-20` ns, 20 ns
wide) and count entries. Divide by window width to get a per-ns
background density directly comparable to the fit's `B`.

Reported in summary CSV columns `sideband_bg_per_ns` and `sideband_counts`.
Not used for flagging in v1 — users eyeball `sideband_bg_per_ns ≈ B` in
the summary CSV. (If a module shows `B` very different from the
sideband density, the fit model is wrong for that module.)

### 10.a Sideband-overlap check

The sideband window must not overlap the fit window, but the fit
window is **centered on the per-module data-driven argmax μ₀** — not on
anything the CLI knows up front. So the check cannot run at argparse
time; `_validate_args` only enforces coarse CLI sanity (e.g.
`sideband_lo < sideband_hi`, sideband inside `[dt_lo, dt_hi]`).

Per-module runtime check inside `sideband_check(acc, ctx, mu0)`
(the function reads `acc.sideband_hist` directly and takes μ₀ from
the fit stage, or the bin-centre argmax of `acc.dt_hist` when the
fit was skipped):

- If `[ctx.sideband_lo, ctx.sideband_hi]` overlaps
  `[mu0 - ctx.fit_window, mu0 + ctx.fit_window]`, do NOT compute a
  background density — set `sideband_bg_per_ns = NaN`,
  `sideband_counts = NaN`, and append the flag `sideband-overlap` to
  the module's flag list.
- Otherwise count entries in the sideband window as normal.

`sideband-overlap` joins the fixed flag vocabulary from §7.b
`flag_module`; it is informational (does not force `no-peak`) because
the fit itself is unaffected — only the data-driven background
cross-check is unavailable.

## 11. Validation — what the implementer must do after writing the script

### 11.a Built-in `--self-test`

Promote the earlier `__main__`-guarded block to a proper
`--self-test` CLI flag.  When passed, `main` ignores the positional
arguments, runs the two synthetic cases below, prints `PASS` / `FAIL`
for each, and exits with status 0 iff both pass.

1. **Synthetic Gaussian + flat background → `ModuleAccumulator`.**
   Build a `RunContext` with the §13.a defaults, then build a
   `ModuleAccumulator` whose `dt_hist` and `sideband_hist` carry
   known-truth samples:
   - draw `N_sig = 5000` samples from `N(μ=0.3 ns, σ=1.1 ns)` and
     `N_bg = 2000` uniform samples on `[-40, 40]` ns (background
     density `2000 / 80 = 25 counts/ns`),
   - concatenate, then `np.histogram` once into `acc.dt_hist` against
     `ctx.dt_edges` (do NOT keep the raw samples; this models exactly
     what `stream_tsv` produces),
   - likewise histogram the subset in `[ctx.sideband_lo, ctx.sideband_hi]`
     into `acc.sideband_hist` against `ctx.sideband_edges`,
   - set `acc.n_entries = int(acc.dt_hist.sum())` and
     `acc.n_total_seen = N_sig + N_bg` (lets the test exercise the
     n_total_seen plumbing too).

   Then call `fit_module(acc, ctx)` and assert:
   - `fit_status == "ok"`
   - `|mu   − 0.3| < 3 · mu_err`      and `|mu   − 0.3| < 0.1 ns`
   - `|sigma − 1.1| < 3 · sigma_err`  and `|sigma − 1.1| < 0.1 ns`
   - `|B_per_ns − 25| / 25 < 0.2`     (within 20% of the truth
     background density — this is the key check that **catches the
     §7.a normalization issue**)
   - `significance > 10`
   - `n_entries == int(acc.dt_hist.sum())`  (pins the §7.c definition)

2. **Minimal time-calib schema parse.**  Build an in-memory JSON string
   with `default = 0.5` and two entries — one keyed by `name = "W4"`
   (resolve to id 1004) and one keyed by `module_id = 42` with no
   `name`.  Write to a tempfile, call `load_time_calib`, assert both
   ids are present with the right offsets.  Then call
   `apply_offsets_chunk` directly with
   `seed_ids = np.array([1004, 42])`,
   `neighbor_ids = np.array([42, 1004])`,
   `dt = np.array([10.0, 10.0])`; assert the two corrected values are
   `dt_corr = 10 - (off_nbr - off_seed)` with the expected signs
   (i.e. one row is `10 - (off_42 - off_1004)` and the other is
   `10 - (off_1004 - off_42)`, so they must sum to 20.0).

The `--self-test` block must have **zero** dependency on an input TSV
or an output directory and should run in under a second so CI can
exercise it on every commit.

### 11.b Validation on real data

1. **Smoke test on an existing TSV.** Any TSV produced by
   `study_hycal_timing.py` (e.g. the file referenced in that script's
   docstring on line 37: `hcal_dt_023867.tsv`). Expected:
   - exits 0, prints a summary line like
     `"<N_modules> modules analysed, <k> flagged (no-peak=…, offset=…, wide=…, low-stats=…)"`.
   - `module_dt_summary_<run>.csv` has exactly one row per distinct
     `neighbor_id` seen in the TSV **after** the stream_tsv filter stack.
2. **Spot-check 2–3 modules.** Open `module_dt_heatmap_<run>.png`, pick
   two modules with a visible dark band near ΔT ≈ 0 and one that looks
   flat; confirm their individual pages in `module_dt_fits_<run>.pdf`
   match the heatmap slice (same peak position, same width).
3. **Fit-success case.** Confirm a high-stats central-crystal module
   (e.g. one of the W-series modules in the dense core) has
   `fit_status = "ok"`, significance ≫ 5, σ in the 0.5–2 ns range,
   and nonzero `mu_err`, `sigma_err` in the summary CSV.
4. **Graceful-failure case.** Confirm a module with only ~10 rows in the
   TSV lands with `fit_status = "low-stats"` and `flags` containing
   `low-stats`, and that the script did not crash.
5. **`--time-calib` round-trip on real data.** Rerun with
   `--time-calib database/hycal_time_offsets/time_calib.json`; confirm
   the resulting μ distribution across modules narrows toward zero
   relative to the uncorrected run (this is the whole point of the
   calibration table; if it doesn't narrow, either our schema parsing
   is wrong or the table was built against a different
   `study_hycal_timing.py` invocation).
6. **Scipy-missing degraded mode.** In a venv without scipy, rerun the
   smoke test; confirm the stderr warning fires, all rows carry
   `fit_status = "scipy-missing"`, and the heatmap + per-module PDF
   are still produced.
7. **CI-friendly rerun.** `--no-pdf` run completes in under a few
   seconds on the sample TSV and produces the three non-PDF artifacts.
8. **Large-TSV streaming sanity (production target).**  On a real
   production-size TSV (e.g. `output/hcal_dt_025320.tsv` at ~324 M
   rows / ~25–40 GB), confirm:
   - the script completes without OOM;
   - peak RSS stays well under a soft 500 MB ceiling (expected
     steady-state footprint is under 100 MB per §7.c);
   - the progress line fires at the `--progress-every` cadence with a
     monotonically shrinking ETA;
   - **chunking does not change results.**  Pick 3–5 modules.  Run
     the script twice on the same TSV with `--max-rows 10_000_000`,
     once at `--chunk-rows 100_000` and once at `--chunk-rows 2_000_000`.
     Assert the per-module `dt_hist` and `sideband_hist` arrays and
     `n_entries` are bitwise-identical between the two runs for every
     selected module id.  This is the off-by-one catch for the
     streaming reader.
   - **Reference-on-subset cross-check.**  Pick one high-stats module.
     Independently tally its row count and a coarse 10-bin dt
     histogram from the first 10 M rows of the TSV using a one-off
     awk / head+csv.reader snippet in a shell session (not committed
     code).  Compare against the script's output under
     `--max-rows 10_000_000 --module-id <that_id>`; both the row count
     and the 10-bin histogram must match exactly.

## 12. Follow-ups explicitly deferred

- **Observed, out of scope — stale C++ comment at
  `prad2dec/src/WaveAnalyzer.cpp:678–681`.**  The comment says "add the
  caller's offset" but the code subtracts it (correct per
  `HyCalTimeCalib.h:8`; the comment is wrong).  Pure comment fix on
  the C++ side; captured here so the finding from
  `docs/research-study-hycal-timing-time-offset.md` is not lost.  Not
  this plan's job.
- **Singles-rate accidentals.** Combine per-module live-time rates with
  the ΔT fit's `B` to produce an *expected* accidentals pedestal from
  first principles. Needs a per-module rate source that this TSV does
  not carry; separable task.
- **Automated writing of a new `hycal_time_offsets/*.json`.** The μ
  column of the summary CSV is the raw input, but the actual decision
  of which modules to re-calibrate and the regularization against the
  reference module (`W495` in the shipped table) is a judgment call
  that belongs in a dedicated calibration script.
- **Per-(seed, neighbour) 2-D μ map.** This diagnostic marginalizes
  over seed. If a specific (seed, neighbour) pair is suspect we would
  need a second script that keys on both IDs; out of scope here.
- **Cross-run comparison / trending.** The summary CSV is designed to
  make this easy (stable schema, run-id in filename), but the trending
  tool itself is a separate deliverable.

## 13. Committed defaults and remaining open questions

### 13.a Defaults committed for v1

The following defaults are now locked in — the implementer should bake
them into the CLI without further consultation.  All remain overridable
via their respective CLI flags.

| flag                | default        | note                                         |
|---------------------|----------------|----------------------------------------------|
| `--dt-window`       | `-40 40` ns    | diagnostic choice; wide enough to show both the peak and a clean sideband at `--dt-bin 0.5`.  Independent of `study_hycal_timing.py`'s `--pre-window`, which actually defaults to **None** (no pre-window — see `study_hycal_timing.py:65–70`), so the TSV's native ΔT range is whatever the WaveAnalyzer produced. |
| `--dt-bin`          | `0.5` ns       | matches typical HyCal timing resolution (~1 ns), ~160 bins total. |
| `--fit-window`      | `8.0` ns       | ±8 ns around argmax — comfortably wider than any plausible σ. |
| `--sideband`        | `-40 -20` ns   | 20 ns wide, well outside `--fit-window` of a near-zero peak. |
| `--min-entries`     | `200`          | low enough to fit most PbGlass outer modules, high enough for curve_fit to converge. |
| `--sig-threshold`   | `5.0`          | 5σ "no-peak" threshold. |
| `--mu-threshold`    | `2.0` ns       | "offset-peak" threshold. |
| `--sigma-threshold` | `3.0` ns       | "wide-peak" threshold. |
| `--secondary-search-window` | `5.0 25.0` ns | ΔT region in which to look for a residual secondary peak (§14.a). |
| `--secondary-sig-threshold` | `3.0`         | Minimum `S₂ / √(B·6σ₂)` to raise the `secondary-peak` flag (§14.a). |
| `--hycal-geometry`  | `$PRAD2_DATABASE_DIR/hycal_map.json`, else `<repo>/database/hycal_map.json` | Source for the §14.b 2-D spatial map.  Resolution mirrors `plot_template_2d_map.py:628–637`.  Missing file → warn + skip spatial map, everything else continues (§14.b). |

### 13.b Still open (do not block v1)

1. **Per-module-type thresholds.** If a target spec exists for PbWO4 vs
   PbGlass timing resolution (likely tighter σ for PbWO4), we may want
   two separate `--sigma-threshold`s keyed by `module_id < 1000` vs
   `>= 1001`. Not in v1; a single global threshold is a correct starting
   point.
2. **Per-module-type roll-up.** Easy to add (module id range split in
   `write_markdown_report`). Not in v1 unless you request it.

## 14. Secondary-peak detection and HyCal spatial map (extension)

Motivation: a production run on 25320 shows that some HyCal modules
have a clean primary ΔT peak at μ ≈ 0, σ ≈ 1 ns AND a secondary peak
near ΔT ≈ 15 ns.  This extension flags those modules, extracts the
secondary peak parameters, and visualises their spatial distribution
on the HyCal face so a user can see at a glance whether the affected
modules cluster (one physics cause) or scatter (multiple).

The extension is additive.  It does **not** change the primary fit,
the streaming architecture, or any existing output.  Everything here
runs after `fit_module` + `sideband_check` have populated the primary
`ModuleStats`.

### 14.a Secondary-peak detection

Signature: `detect_secondary_peak(acc: ModuleAccumulator, stats: ModuleStats, ctx: RunContext) -> None`
(mutates `stats` in place with the four `secondary_*` fields and
appends `secondary-peak` to `stats.flags` when the significance clears
`ctx.secondary_sig_thr`).

**Preconditions.**  Run only when ALL of the following hold:

1. `stats.fit_status == "ok"` (the stricter of {ok, at-bound} — an
   at-bound primary may be sitting on top of what IS the secondary,
   and we do not want to double-count).
2. `|stats.mu| <= 2 * ctx.mu_thr`     (default: `|μ_primary| <= 4 ns`).
3. `stats.sigma <= 1.5 * ctx.sigma_thr` (default: `σ_primary <= 4.5 ns`).

Rationale for 2 and 3: a broad or offset primary fit will leave a
badly-shaped residual (one-sided bump from an off-centre Gaussian, or
long tails from an over-wide one) that `find_peaks` happily picks up
as a false secondary.  The guardrail keeps secondary detection
restricted to modules where the primary is already clean enough that
the residual means what we think it means.

If any precondition fails, set all six `secondary_*` fields to NaN
and return without touching `stats.flags`.  Rationale for the "no
primary" case: we need the primary `(A, μ, σ, B)` to build the
residual; without a good primary fit the residual is meaningless.

**Procedure.**

1. **Build the residual density.**
   - `y        = acc.dt_hist / ctx.dt_bin`                 # counts/ns
   - `y_fit    = gaussian_plus_const(centres, A, mu, sigma, B)`
     (using the primary-fit parameters already in `stats`)
   - `residual = y - y_fit`

2. **Restrict to the secondary search window.**
   - `mask = (centres >= ctx.secondary_lo) & (centres <= ctx.secondary_hi)`
   - If `mask` captures fewer than 5 bins → set all `secondary_*` to
     NaN and return (silently; this is not a flag).
   - Semantics: `ctx.secondary_lo` / `ctx.secondary_hi` are in
     **absolute corrected ΔT coordinates** (same coordinate system as
     `acc.dt_hist`'s bin edges — i.e. post-offset when `--time-calib`
     is used, raw otherwise), NOT relative to `stats.mu`.  This is
     already how the code naturally works because `centres` is derived
     from `ctx.dt_edges`; the note is here so the implementer does not
     "helpfully" shift by `stats.mu`.

3. **Candidate search with `scipy.signal.find_peaks`.**
   - `local_B_counts = max(stats.B_per_ns * ctx.dt_bin, 1.0)`
   - `prominence    = max(3.0 * sqrt(local_B_counts), 5.0)`      # in count units
   - `prominence_density = prominence / ctx.dt_bin`              # convert to counts/ns
     (because `residual` is in counts/ns, the prominence threshold
     must be in the same units — this is a §7.a-style normalization
     cross-check waiting to happen if the implementer forgets)
   - `peaks, props = find_peaks(residual[mask], prominence=prominence_density)`
   - If `len(peaks) == 0` → set `secondary_*` to NaN, return.
   - Pick the tallest candidate:
     `tallest = peaks[np.argmax(residual[mask][peaks])]`
   - `mu2_seed = centres[mask][tallest]`

4. **Fit a single Gaussian to the residual around the candidate.**
   Use the SAME density + `absolute_sigma=True` + `maxfev=5000`
   convention as the primary fit (§9).  Model is a bare Gaussian
   (no constant term — the primary B has already been subtracted
   into the residual):
   `g(x; A2, μ2, σ2) = A2 · exp(-½((x-μ2)/σ2)²)`
   Seeds and bounds:
   - `A2_seed    = residual[mask][tallest]`                        (counts/ns)
   - `sigma2_seed = max(2.0 * ctx.dt_bin, 0.5)`                     (ns)
   - bounds `(low, high)` = `([0.0, ctx.secondary_lo, ctx.dt_bin], [np.inf, ctx.secondary_hi, 5.0])`
   - fit on `centres[mask]` → `residual[mask]` with
     `y_err_mask = sqrt(acc.dt_hist[mask] + local_B_counts) / ctx.dt_bin`
     (approximate Poisson error on the pre-subtraction counts, in
     counts/ns to match `residual`; `local_B_counts` from step 3
     already carries a `max(..., 1.0)` floor so empty bins have
     non-zero weight.  The subtraction of a fit does not add stat
     uncertainty here because we treat the primary-fit params as
     fixed.)

   On `RuntimeError` or singular covariance → set `secondary_*` to
   NaN, return.

5. **Derived quantities and flagging.**
   - `S2 = A2 * sigma2 * sqrt(2π) * erf(3/√2)`                     (counts)
   - `B_w2 = stats.B_per_ns * 6.0 * sigma2`                         (counts)
   - `secondary_significance = S2 / sqrt(max(B_w2, 1.0))`
   - Populate `stats.secondary_mu`, `stats.secondary_sigma`,
     `stats.secondary_amp_per_ns = A2`,
     `stats.secondary_significance`.
   - Also store `secondary_mu_err`, `secondary_sigma_err` on `stats`
     for the markdown report (§8.b) — fields added to `ModuleStats`
     in §14.d.
   - If `secondary_significance >= ctx.secondary_sig_thr`, append
     `secondary-peak` to `stats.flags`.  `fit_status` is NOT touched.

**Non-goals carried from the user's amendment brief.**  No
double-Gaussian primary fit; no tertiary peak search; no automatic
exclusion of secondary-peak modules from any downstream flow.  The
flag is informational.

### 14.b HyCal geometry loader

Signature: `load_hycal_geometry(path: Path) -> dict[int, ModuleGeom]`
where `ModuleGeom` is:

```
@dataclass
class ModuleGeom:
    module_id:   int        # PrimEx id (G: 1-576, W: 1001-2152)
    name:        str        # "G<n>" or "W<n>"
    x_mm:        float      # module centre x on the HyCal face (mm)
    y_mm:        float      # module centre y on the HyCal face (mm)
    size_x_mm:   float      # module face width  (mm) — "sx" in JSON
    size_y_mm:   float      # module face height (mm) — "sy" in JSON
    module_type: str        # "PbWO4" | "PbGlass" | "LMS" | "Veto" | "Unknown"
    sector:      str        # "Top" | "Bottom" | "Left" | "Right" | "Center"
                            #   — the "geo.sec" field in the JSON
```

**JSON schema** (confirmed by reading
`database/hycal_map.json:1–30` directly and the C++ parser at
`prad2det/src/HyCalSystem.cpp:115–134`; no ambiguity, no research gap):

```
[
  {
    "n":   "<name>",         # top-level: module name string
    "t":   "<type>",         # top-level: "PbWO4" | "PbGlass" | "LMS" | "Veto"
    "geo": {
      "sx":  <float, mm>,    # module face width
      "sy":  <float, mm>,    # module face height
      "x":   <float, mm>,    # module centre x on HyCal face
      "y":   <float, mm>,    # module centre y on HyCal face
      "sec": "<sector>",     # "Top" | "Bottom" | "Left" | "Right" | "Center"
      "row": <int>,          # 1-based row within sector (not needed for us)
      "col": <int>           # 1-based col within sector (not needed for us)
    },
    "daq": { ... }           # optional; we never read it
  },
  ...
]
```

**PrimEx id derivation** (do not import `prad2py` — pure Python, same
rule as the C++ `name_to_id` in `HyCalSystem.cpp:40–47`):
- `name` starts with `G` / `g` → id = `int(name[1:])`                 (range 1–576)
- `name` starts with `W` / `w` → id = `int(name[1:]) + 1000`          (range 1001–2152)
- anything else → skip the module with a `[WARN]` line.

**Resolution of the default `--hycal-geometry` path** mirrors the
convention established by `plot_template_2d_map.py:628–637`:
`$PRAD2_DATABASE_DIR/hycal_map.json` when the env var is set,
otherwise `<script_dir>/../../database/hycal_map.json`.  Both are
pre-resolved (`Path.resolve()`) at argparse time so later `cwd`
changes do not break anything.

**Missing-file policy.**  If `--hycal-geometry` points at a path that
does not exist, `load_hycal_geometry` returns `{}` (an empty dict) and
prints one `[WARN]` line to stderr.  Downstream, `plot_hycal_spatial_map`
sees the empty dict and skips the entire spatial-map output while
leaving every other artifact in place.  `--no-spatial-map` has the
same effect but silently.

**Malformed-JSON policy.**  Same warn-and-skip behaviour applies to
unreadable or ill-structured JSON files:

- JSON parse error (`json.JSONDecodeError`) → return `{}`, print
  `[WARN] load_hycal_geometry: cannot parse <path>: <error>`.
- Top-level value is not a list → return `{}`, print
  `[WARN] load_hycal_geometry: <path> is not a JSON array`.
- Individual module entry missing required keys (`n`, `t`, `geo`, or
  any of `geo.sx` / `geo.sy` / `geo.x` / `geo.y`) → skip that entry
  with `[WARN] load_hycal_geometry: skipping entry <idx>: missing
  key <name>`; keep going on the rest of the file.
- Individual module entry has a name that does not match the
  `G<n>` / `W<n>` PrimEx-id rule → skip that entry with the existing
  warn line from the "PrimEx id derivation" paragraph above.

The script never crashes on a bad geometry file.  The downstream
spatial-map step treats an empty dict identically to a missing file.

**Beam-hole handling.**  The PRad target hole lives at the centre of
the HyCal face.  `hycal_map.json` simply omits the modules that would
sit there — there is no explicit "hole" entry.  The spatial map
therefore gets a hole "for free" because the plotter only draws
modules present in the geometry dict.  No special-case code required.

### 14.c Spatial map figure and markdown roll-up

Signature:
`plot_hycal_spatial_map(stats_by_id: dict[int, ModuleStats], geometry: dict[int, ModuleGeom], ctx: RunContext, path_stem: Path) -> None`.

**Figure layout.**  One matplotlib `Figure` with three subplots in a
single row (`fig, axes = plt.subplots(1, 3, figsize=(18, 6))`),
sharing the HyCal footprint via `ax.set_aspect("equal")` and equal
`xlim`/`ylim` computed from `max(|x|+sx/2, |y|+sy/2)` across all
loaded modules.

- **(a) Binary flag map (axes[0]).**  For every module in `geometry`:
  - if the module id is in `stats_by_id` AND `"secondary-peak"` is in
    `stats_by_id[id].flags` → red filled square at `(x_mm, y_mm)`
    sized `(size_x_mm, size_y_mm)`.
  - else if the id is in `stats_by_id` → grey filled square (same
    size).
  - else → empty outline (black edge, no fill) at the module's slot
    (module was present in geometry but excluded from analysis, e.g.
    by `--module-id`).
  Title: "Modules with `secondary-peak` flag (red)".

- **(b) Secondary μ₂ map (axes[1]).**  For modules with a
  non-NaN `secondary_mu`: colour from the `viridis` colormap scaled
  over `[ctx.secondary_lo, ctx.secondary_hi]`.  Modules analysed but
  without a secondary peak → grey.  Modules not analysed → empty
  outline.  Add a colorbar labelled "secondary μ₂ (ns)".

- **(c) Secondary significance map (axes[2]).**  Same geometry, but
  colour from `plasma` scaled over `[0, max(10.0,
  np.nanmax(secondary_significance))]`.  Grey / empty outline rules
  identical to (b).  Colorbar label: "S₂ / √(B·6σ₂)".

**All-NaN empty case.**  If every module's `secondary_significance`
is NaN (no secondary peaks detected run-wide), `np.nanmax` would
raise a `RuntimeWarning` and return NaN, which breaks the colormap
scaling.  Guard against this: before computing the max, filter to
non-NaN values; if the filtered array is empty, short-circuit panel
(c) by filling every module with the grey "analysed-but-no-secondary"
colour and overlaying a centred text annotation "No secondary peaks
detected".  Panels (a) and (b) naturally render all-grey in this
case (panel (a) has no red squares because no `secondary-peak` flags
exist; panel (b) has no coloured squares because every
`secondary_mu` is NaN) — no special-case code needed for them.  The
figure is still produced so the user sees an explicit "nothing
detected" result rather than a missing artifact.

**Rectangle drawing.**  Use `matplotlib.patches.Rectangle` with
`xy=(x_mm - size_x_mm/2, y_mm - size_y_mm/2)` and the module's
`(size_x_mm, size_y_mm)` so PbWO4 (sx=sy≈20.5 mm) and PbGlass
(sx=sy=38.15 mm) draw at their correct physical scale.  This gives
the "dot size varies with module type" behaviour the user asked for
without any manual bookkeeping.

**Output files.**  Both `.png` and `.pdf`, written to
`<OUTPUT_DIR>/hycal_secondary_peak_map_<run>[_<tag>].{png,pdf}`.

**Markdown report section** (new, inserted at the top of
`<OUTPUT_DIR>/module_dt_flags_<run>[_<tag>].md`, before any per-flag
section):

```
## Modules with secondary peak

Sorted by secondary_significance (descending).  Threshold:
`--secondary-sig-threshold = <value>`; search window:
`--secondary-search-window = <lo> <hi>` ns.

| name  | id   | μ₂ ± err (ns)   | σ₂ ± err (ns)  | sig   |
|-------|------|-----------------|----------------|-------|
| W1234 | 2234 | 15.1 ± 0.3      | 1.1 ± 0.2      | 12.4  |
| ...
```

If no modules carry the `secondary-peak` flag, write "No modules
exceeded the secondary-peak significance threshold." in place of the
table.

### 14.d `ModuleStats` and `RunContext` additions

Added to `@dataclass ModuleStats` (default NaN for all floats):

```
secondary_mu:          float = float("nan")
secondary_mu_err:      float = float("nan")
secondary_sigma:       float = float("nan")
secondary_sigma_err:   float = float("nan")
secondary_amp_per_ns:  float = float("nan")
secondary_significance:float = float("nan")
```

Added to `@dataclass RunContext`:

```
secondary_lo:        float          # from --secondary-search-window[0]
secondary_hi:        float          # from --secondary-search-window[1]
secondary_sig_thr:   float          # from --secondary-sig-threshold
hycal_geometry_path: Path | None    # from --hycal-geometry (None = fall through)
spatial_map_enabled: bool           # False iff --no-spatial-map OR geometry empty
```

### 14.e CLI flags (summary; defaults locked in §13.a)

```
--secondary-search-window LO HI
                          ΔT region (ns) to search for a residual
                          secondary peak (default 5.0 25.0).

--secondary-sig-threshold T
                          Minimum S₂ / √(B·6σ₂) for the
                          secondary-peak flag (default 3.0).

--hycal-geometry PATH     Path to hycal_map.json (default:
                          $PRAD2_DATABASE_DIR/hycal_map.json, else
                          <repo>/database/hycal_map.json; mirrors
                          plot_template_2d_map.py:628–637).  Missing
                          file → warn + skip spatial map.

--no-spatial-map          Skip the HyCal spatial-map output.
                          Useful for CI, for iteration on other
                          outputs, or when geometry is unavailable.
```

`_validate_args` additions:
- `--secondary-search-window LO HI` must satisfy `LO < HI` and must
  sit inside `[ctx.dt_lo, ctx.dt_hi]`.
- `--secondary-sig-threshold` must be non-negative.
- `--hycal-geometry` is NOT required to exist at argparse time
  (missing-file is a warn-and-skip, not a hard failure).

### 14.f Scipy-missing interaction

`scipy.signal.find_peaks` and the secondary `curve_fit` call both
require scipy.  Under the §5 degraded mode:
- `detect_secondary_peak` becomes a no-op (sets all four `secondary_*`
  fields to NaN; no flag appended).
- The spatial map (§14.c) still renders panel (a) using only
  `stats.flags`, which will be empty of `secondary-peak` under the
  degraded mode — resulting in an all-grey panel.  Panels (b) and (c)
  are empty of coloured modules too; the figure is still produced so
  the user sees a trivial "nothing detected" result rather than a
  missing artifact.

### 14.g Main-flow integration point

`main()` step ordering (extends §7.b main() list; existing step
numbers are unchanged):

```
after step 5 (per-accumulator fit_module + sideband_check +
flag_module), insert:

  5a. for each accumulator: detect_secondary_peak(acc, stats, ctx).

before step 6 (write artifacts), load geometry once:

  5b. geometry = load_hycal_geometry(ctx.hycal_geometry_path)
      if ctx.spatial_map_enabled is still True AND geometry is non-empty.

step 6 artifact list gains:
      hycal_secondary_peak_map_<run>[_<tag>].{png,pdf}
      (skipped when ctx.spatial_map_enabled is False or geometry empty)
```

### 14.h Validation hooks (extend §11)

**§11.a self-test — new case 3** (added alongside the existing two):

- Build a `ModuleAccumulator` whose `dt_hist` carries three
  components: primary Gaussian `A1=1000 counts/ns, μ1=0 ns, σ1=1 ns`;
  secondary Gaussian `A2=200 counts/ns, μ2=15 ns, σ2=1 ns`; flat
  background `B=25 counts/ns`.  Draw ~N_sig1 ≈ 5000, N_sig2 ≈ 1000,
  N_bg ≈ 2000 samples (counts; the exact numbers follow from the
  amplitudes given the default `--dt-window` and `--dt-bin`), histogram
  once into `acc.dt_hist` against `ctx.dt_edges`.
- Call `fit_module(acc, ctx)` then `detect_secondary_peak(acc, stats, ctx)`.
- Assert:
  - primary `fit_status == "ok"`, `|mu − 0| < 0.1`, `|sigma − 1| < 0.1`;
  - `|secondary_mu − 15| < 0.5 ns`,
  - `secondary_significance > 10`,
  - `"secondary-peak" in stats.flags`.

**§11.b real-data validation — new steps**:

9. **False-positive sanity on clean modules.**  Identify two or three
   PbWO4 central-crystal modules from the first production run whose
   per-module page shows only a single clean primary peak (no visible
   residual structure in the ±15 ns neighbourhood).  Confirm the
   script does NOT raise `secondary-peak` on those modules at the
   default threshold.
10. **Spatial-map visual check.**  Confirm
    `hycal_secondary_peak_map_<run>.png` renders:
    - module squares at physically plausible positions (overall
      octagonal / rectangular HyCal footprint, PbWO4 core smaller
      than PbGlass ring),
    - a central beam-hole gap where `hycal_map.json` omits modules,
    - panel (a) coloured consistently with the `flags` column of
      `module_dt_summary_<run>.csv`.
11. **Missing-geometry degraded path.**  Rerun with
    `--hycal-geometry /tmp/does_not_exist.json`; confirm the script
    prints a single `[WARN]` line, skips the spatial map, and
    completes successfully with all other artifacts intact.

---

## Todo

### Phase 1 — scaffolding

- [ ] Create `analysis/pyscripts/module_dt_diagnostic.py` with module
      docstring (deps, usage synopsis, output artifact list from §8,
      pointer to this plan file, call-out of the §7.a density fit
      convention and the §5 scipy-missing degraded mode).
- [ ] Add top-level constants: `TSV_COLUMNS` (frozen schema from §2),
      `REQUIRED_COLUMNS = {"seed_id", "neighbor_id", "seed_time",
      "neighbor_time", "dt", "seed_energy", "neighbor_energy",
      "dist_q"}`.
- [ ] Add `HAVE_SCIPY` try/except around `from scipy.optimize import curve_fit`
      mirroring `benchmark_hycal_timing.py:44–48`.  The try/except at
      module scope sets the flag SILENTLY — do NOT print anything at
      import time (importing the module for reuse should not spam
      stderr).  The user-facing warning line from §5 is printed from
      the start of `main()` (Phase 5 Todo) when `not HAVE_SCIPY`.
- [ ] Define `@dataclass ModuleStats` per §7.b (fields include
      `mu_err`, `sigma_err`, `A_per_ns`, `B_per_ns`, `n_total_seen`).
- [ ] Define `@dataclass RunContext` per §7.b (fields include
      `module_id_filter`, `min_seed_energy`, `max_dist_q`,
      `scipy_available`, `default_offset`, the streaming knobs
      `chunk_rows`, `progress_every`, `max_rows`, AND the derived
      `dt_edges` / `sideband_edges` arrays computed once from
      `(dt_lo, dt_hi, dt_bin)` and `(sideband_lo, sideband_hi, dt_bin)`
      so `stream_tsv`, `fit_module`, and `sideband_check` all share
      the same bin partition).
- [ ] Define `@dataclass ModuleAccumulator` per §7.b (fields
      `module_id`, `dt_hist: np.ndarray[int64]`,
      `sideband_hist: np.ndarray[int64]`, `n_entries`, `n_total_seen`).
      Add a factory `ModuleAccumulator.empty(module_id, ctx)` that
      sizes `dt_hist = np.zeros(len(ctx.dt_edges) - 1, dtype=np.int64)`
      and `sideband_hist = np.zeros(len(ctx.sideband_edges) - 1, dtype=np.int64)`.

### Phase 2 — streaming I/O

- [ ] Implement `load_time_calib(path)`:
  - parse JSON per `HyCalTimeCalib.h:37–53` schema;
  - resolve each entry to a PrimEx id via `module_id` when present,
    otherwise strip `G`/`W` prefix from `name` and add `0` or `1000`;
  - compute SHA-1 of the file bytes (hashlib) for provenance;
  - return `(offsets: dict[int, float], default: float, sha1: str)`.
- [ ] Implement `apply_offsets_chunk(seed_ids, neighbor_ids, dt, offsets, default)`
      returning `dt_corr = dt − (off_neighbour − off_seed)` as a
      vectorised numpy op per §7.b.  No-op (return `dt` unchanged) when
      `offsets` is empty.  Exercised directly by `--self-test` (§11.a).
- [ ] Implement `stream_tsv(path, csv_mode, no_header, ctx) ->
      dict[int, ModuleAccumulator]` per §7.b and §7.c.  Checklist:
  - stdlib `csv.reader` (TSV delimiter by default; `--csv` switches
    to CSV); validate the header against `REQUIRED_COLUMNS` unless
    `--no-header` was passed;
  - buffer up to `ctx.chunk_rows` raw rows in Python lists, then
    convert to numpy arrays (int32 ids, float32 numeric cols) — keep
    only the columns actually used downstream (seed_id, neighbor_id,
    dt, seed_energy, dist_q);
  - apply the four filter masks in the committed order
    (seed_energy → dist_q → module_id → dt-window on RAW dt);
  - update per-module `n_total_seen` AFTER the first three masks and
    BEFORE the dt-window mask, so n_total_seen reflects post-filter
    pre-window counts (see §7.c correctness invariants);
  - call `apply_offsets_chunk` once per chunk (no-op when
    `ctx.offsets` empty);
  - for every unique `neighbor_id` in the surviving chunk, slice,
    `np.histogram` against the pre-computed `dt_edges` and
    `sideband_edges`, and `+=` into `acc.dt_hist` and
    `acc.sideband_hist`; update `acc.n_entries`;
  - every `ctx.progress_every × 1_000_000` rows read (skip when
    `progress_every == 0`) print the progress line from §6 (rows read,
    rows kept, elapsed, ETA from bytes consumed vs. file size);
  - honor `ctx.max_rows` as a hard cap on rows read from the input;
  - fail loudly (SystemExit with a non-zero code) on a required
    column missing, on a row that cannot be parsed, or on an
    unexpectedly huge module id (sanity bound 10_000).

### Phase 3 — analysis (on accumulators)

Note: filtering and per-chunk grouping now live inside `stream_tsv`
(Phase 2); there is no separate `filter_rows` or `group_by_module`
function.  The Phase 3 analysis stages consume `ModuleAccumulator`s
only, never raw rows.

- [ ] Implement `gaussian_plus_const(x, A, mu, sigma, B)` as the
      density model from §7.a.
- [ ] Implement `fit_module(acc, ctx)` per §9.a, producing exactly the
      five `fit_status` values from §9.b (`ok`, `low-stats`,
      `scipy-missing`, `fit-failed`, `at-bound`).  Input is a
      `ModuleAccumulator`; use `acc.dt_hist` directly (no re-histogram),
      compute `y = acc.dt_hist / ctx.dt_bin`, `y_err = max(sqrt(dt_hist), 1) / dt_bin`.
      Use `A0 = max(0.0, y.max() − B0)`, `curve_fit(..., maxfev=5000,
      sigma=y_err, absolute_sigma=True)`, μ bounds `[dt_lo, dt_hi]`
      (not symmetric), and the at-bound tolerance
      `max(0.5 · dt_bin, 1e-6)` applied to BOTH μ and σ (bin-sized
      absolute tolerance, NOT a fraction of the bound span — the old
      `0.5 · (upper − lower)` σ tolerance would mark realistic HyCal
      widths as at-bound; see §9.a step 7).  Capture `mu_err`,
      `sigma_err` from `sqrt(diag(pcov))` with NaN fallback on singular
      covariance.  `n_entries = int(acc.dt_hist.sum())`.
- [ ] Implement `sideband_check(acc, ctx, mu0)` returning
      `(counts, counts / window_width_ns, overlap_bool)`.  Input is a
      `ModuleAccumulator`; `counts = int(acc.sideband_hist.sum())`.
      Units must match `B_per_ns` from the fit (counts/ns) — this is
      the normalization cross-check from §7.a.  Per §10.a: if the
      sideband window overlaps `[mu0 ± fit_window]`, return
      `(NaN, NaN, True)` and let the caller append `sideband-overlap`
      to the flag list.  `mu0` is the fitted μ when `fit_status == "ok"`,
      otherwise the bin-centre argmax of `acc.dt_hist`.
- [ ] Implement `flag_module(stats, ctx)` applying flags per §7.b in a
      fixed order (low-stats → no-peak → offset-peak → wide-peak →
      sideband-overlap → secondary-peak).  Suppress
      `offset-peak`/`wide-peak` when `fit_status != "ok"`.
      `sideband-overlap` is appended by `sideband_check` directly;
      `secondary-peak` is appended by `detect_secondary_peak` (§14.a /
      Phase 8).  flag_module preserves both in the fixed order and
      never re-derives them.

### Phase 4 — outputs

- [ ] Implement `plot_heatmap(accumulators, ctx, path_stem)` writing
      both PNG and PDF.  Stack `acc.dt_hist` rows into a 2-D array
      (log-norm color scale, `vmin=1`; vertical line at ΔT=0).  Thin
      out y-tick labels if `len(ids) > 50`.  No per-row re-histogram.
- [ ] Implement `plot_per_module_pdf(accumulators, stats_by_id, ctx, path)`
      using `matplotlib.backends.backend_pdf.PdfPages`, one page per
      module (sorted by id). Plot the DENSITY histogram
      (`acc.dt_hist / dt_bin`; counts/ns, matching the fit's units) +
      fit curve + μ, μ±3σ guidelines + text box annotated with
      `μ ± dμ`, `σ ± dσ`, `B_per_ns`, `S`, significance, `fit_status`,
      flags.  Honor `--no-pdf`.
- [ ] Implement `write_summary_csv(stats_list, ctx, path)` with the
      `#`-prefixed metadata comment header described in §7.b and §8.a
      (incl. filter settings, scipy availability, hycal-geometry path,
      secondary-peak thresholds).  Columns per §8.a (authoritative
      schema, including the four `secondary_*` columns).
- [ ] Implement `write_markdown_report(stats_list, ctx, path)` grouped by
      flag reason, sorted worst-first per group; include the
      scipy-missing warning in the metadata block when applicable.
      Each flagged-module line follows the §8.b format (name + id +
      n_entries/n_total_seen + fit summary + flags), and for modules
      flagged `secondary-peak` append the §8.b "secondary μ=... σ=...
      sig=..." indented line.  Insert the §14.c "Modules with
      secondary peak" roll-up table at the top of the report (before
      any per-flag section).

### Phase 5 — CLI and glue

- [ ] Implement `sniff_run_id(input_path)` with the `prad_<run>` regex
      (local; do not import `_common`).
- [ ] Implement `_validate_args(ns)` catching only CLI-time sanity
      issues (the sideband-vs-fit-window overlap is per-module and
      data-driven, so it is checked at sideband_check time per §10.a,
      NOT here).  Catch:
        - `sideband_lo < sideband_hi` and sideband inside `[dt_lo, dt_hi]`,
        - `dt_bin >= (dt_hi − dt_lo) / 10`,
        - negative thresholds,
        - non-existent `--time-calib` file,
        - `--dt-window` with `lo >= hi`,
        - per §14.e: `--secondary-search-window LO HI` with `LO < HI`
          AND inside `[dt_lo, dt_hi]`;
          `--secondary-sig-threshold >= 0`.
          `--hycal-geometry` is NOT required to exist (warn-and-skip
          per §14.b).
- [ ] Add `--self-test` CLI flag; wire it to the §11.a `_self_test()`
      function so `main` short-circuits before any TSV I/O.
- [ ] Wire the three streaming CLI flags: `--chunk-rows` (int,
      default 1_000_000), `--progress-every` (int, default 10;
      `0 = silent`), `--max-rows` (int, default 0 = no cap).  Validate
      all three in `_validate_args` (non-negative; `chunk_rows >= 1`).
- [ ] Wire the §14.e secondary-peak CLI flags:
      `--secondary-search-window LO HI` (default `5.0 25.0` ns),
      `--secondary-sig-threshold` (default `3.0`),
      `--hycal-geometry PATH` (default resolved via
      `plot_template_2d_map.py:628–637`-style helper), `--no-spatial-map`
      switch.  Resolve the default geometry path at argparse time with
      `Path.resolve()`.
- [ ] Implement `main(argv=None)` wiring per §7.b main() step list
      augmented by §14.g (parse → optional load_time_calib →
      stream_tsv → per-accumulator fit/sideband/flag →
      detect_secondary_peak (step 5a) → optional load_hycal_geometry
      (step 5b) → write artifacts incl. hycal_secondary_peak_map).
      Return non-zero only when zero rows survive filtering (clear CI
      failure signal).  Immediately after argparse succeeds and
      before any heavy I/O, if `not HAVE_SCIPY` emit the §5 one-line
      warning to stderr.
- [ ] Add `if __name__ == "__main__": sys.exit(main())`.

### Phase 6 — validation (per §11)

- [ ] Implement and run `--self-test` (§11.a):
  - synthetic Gaussian + flat background with known truth, verify
    recovered μ, σ, `B_per_ns` within tolerance (this is the
    normalization check that catches §7.a mistakes);
  - minimal time-calib JSON with one `name`-keyed and one
    `module_id`-keyed entry, verify both resolution paths and the sign
    of correction (`dt_corr = dt − (off_nbr − off_seed)`);
  - §14.h case 3: synthetic primary (μ1=0, σ1=1) + secondary (μ2=15,
    σ2=1) + flat background; assert primary `fit_status == "ok"`,
    `|secondary_mu − 15| < 0.5 ns`, `secondary_significance > 10`,
    `"secondary-peak" in stats.flags`.
- [ ] Smoke test on an existing TSV; confirm one summary-CSV row per
      distinct `neighbor_id` that survives `stream_tsv`'s filter stack.
- [ ] Open the heatmap PNG and spot-check 2–3 module slices against
      their per-module pages in the fits PDF.
- [ ] Confirm one high-stats central-crystal module gets
      `fit_status = "ok"` with sensible μ, σ, significance, and nonzero
      `mu_err`, `sigma_err`.
- [ ] Confirm a sparse module gets `fit_status = "low-stats"` with no crash.
- [ ] Rerun with `--time-calib database/hycal_time_offsets/time_calib.json`
      and confirm the μ distribution across modules narrows toward 0.
- [ ] In a venv without scipy, rerun the smoke test and confirm the
      degraded mode: stderr warning fires, all rows carry
      `fit_status = "scipy-missing"`, plots are still produced.
- [ ] Rerun with `--no-pdf` and confirm sub-few-second completion for CI use.
- [ ] **Large-TSV streaming sanity** (per §11.b step 8).  Run on a
      production TSV (e.g. `output/hcal_dt_025320.tsv`, ~324 M rows):
      confirm no OOM, peak RSS stays well under 500 MB, progress lines
      fire at the `--progress-every` cadence with monotonically
      shrinking ETA.
- [ ] **Chunk-size invariance** (per §11.b step 8).  Pick 3–5 modules,
      run twice with `--max-rows 10_000_000` at `--chunk-rows 100_000`
      and `--chunk-rows 2_000_000`; assert the per-module `dt_hist`,
      `sideband_hist`, and `n_entries` match bitwise for every selected
      id.  This catches off-by-one errors in `stream_tsv`.
- [ ] **Reference-on-subset cross-check** (per §11.b step 8).  Tally
      row count and a coarse 10-bin dt histogram for one high-stats
      module from the first 10 M rows via an ad-hoc shell/awk snippet;
      compare against the script's output at
      `--max-rows 10_000_000 --module-id <that_id>`.  Both must match
      exactly.
- [ ] §14.h step 9: **false-positive sanity on clean modules**.  Pick
      2–3 PbWO4 central-crystal modules whose per-module page shows
      only a clean primary; confirm `secondary-peak` is NOT raised at
      the default threshold.
- [ ] §14.h step 10: **spatial-map visual check**.  Open
      `hycal_secondary_peak_map_<run>.png`; confirm the HyCal
      footprint, PbWO4-vs-PbGlass square sizing, central beam-hole
      gap, and panel (a) consistency with the summary CSV `flags`
      column.
- [ ] §14.h step 11: **missing-geometry degraded path**.  Rerun with
      `--hycal-geometry /tmp/does_not_exist.json`; confirm one
      `[WARN]` line, spatial map skipped, other artifacts intact.

### Phase 8 — secondary-peak detection + HyCal spatial map (§14)

Implementation items for the §14 extension.  Depends on Phase 3 (the
primary `fit_module` must be in place first) and reuses the §5
`HAVE_SCIPY` degraded-mode gate.

**On Phase 8's placement.**  Phase 8 is listed last among code phases
(before Phase 7 documentation) purely for logical grouping of the §14
extension — all §14 code lives in one block of Todos so an implementer
reads the extension coherently.  It is **not** a staged-delivery
escape hatch.  v1 ships with §14 integrated:
- Phase 5 wires `detect_secondary_peak` into `main()` step 5a.
- Phase 6 includes the §11.a case-3 secondary-peak self-test and the
  §14.h step 9–11 real-data validations.
- `--no-spatial-map` only skips the spatial-map output; secondary
  detection itself always runs so the CSV / markdown / per-module PDF
  stay consistent across invocations.

There is deliberately NO `--no-secondary-peak` flag.  The user asked
for this diagnostic; it is part of v1.

- [ ] Extend `@dataclass ModuleStats` with the six new fields from
      §14.d (`secondary_mu`, `secondary_mu_err`, `secondary_sigma`,
      `secondary_sigma_err`, `secondary_amp_per_ns`,
      `secondary_significance`), defaulting to NaN.
- [ ] Extend `@dataclass RunContext` with the five new fields from
      §14.d (`secondary_lo`, `secondary_hi`, `secondary_sig_thr`,
      `hycal_geometry_path`, `spatial_map_enabled`).
- [ ] Define `@dataclass ModuleGeom` per §14.b.
- [ ] Implement `load_hycal_geometry(path)` per §14.b: parse the JSON
      array, derive PrimEx ids via the `G/W + numeric suffix` rule
      (do NOT import prad2py), skip + warn on unparseable entries,
      return `{}` on missing file.  Resolution of the default path
      mirrors `plot_template_2d_map.py:628–637`.
- [ ] Implement `detect_secondary_peak(acc, stats, ctx)` per §14.a:
      preconditions per §14.a (strict ok + bounded μ + bounded σ;
      do NOT re-enumerate here), build the residual density, run
      `scipy.signal.find_peaks` with
      `prominence_density = max(3·sqrt(local_B_counts), 5) / dt_bin`
      (keep the per-ns normalization in step with §7.a), fit a bare
      Gaussian with `maxfev=5000` and `absolute_sigma=True`, derive
      `S2 / sqrt(max(B·6σ2, 1))`, append `secondary-peak` to
      `stats.flags` when the significance clears
      `ctx.secondary_sig_thr`.  No scipy → no-op (NaN fields, no flag).
- [ ] Implement `plot_hycal_spatial_map(stats_by_id, geometry, ctx, path_stem)`
      per §14.c: three subplots sharing the HyCal footprint, red /
      grey / outline-only convention, PbWO4 vs PbGlass sized via
      `(size_x_mm, size_y_mm)` rectangles, colorbars labeled per
      §14.c.  Write both PNG and PDF.  Early-return (and do nothing)
      when `geometry` is empty or `ctx.spatial_map_enabled` is False.
- [ ] Insert the §14.c "Modules with secondary peak" roll-up at the
      top of the markdown report (handled inside
      `write_markdown_report`; link back to §8.b for the per-line
      format with the indented secondary-peak follow-up line).
- [ ] Thread `detect_secondary_peak` into the `main()` per-accumulator
      loop as step 5a (§14.g).  Thread `load_hycal_geometry` +
      `plot_hycal_spatial_map` as step 5b and an artifact in step 6.

### Phase 7 — documentation (lightweight)

- [ ] *(optional)* Add a one-line entry to `analysis/pyscripts/Commands.md`
      under a new "Timing diagnostics" section pointing at this script
      and this plan file.  Skip if `Commands.md` conventions have
      drifted since this plan was written; defer to a doc-writer pass.
- [ ] If this plan surfaces changes during implementation that affect
      the deferred items (§12), capture them as a short "Follow-ups"
      section at the bottom of the script's module docstring so the
      next person has a map.
