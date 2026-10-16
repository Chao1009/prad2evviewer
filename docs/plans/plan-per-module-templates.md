# Plan — per-module pulse templates in two Python benchmarks

**Status:** draft (iteration 2, pre-approval — planning-critic blocking
fixes applied: missing `json` import, `--template` → `--template-file`
in recipes, per-channel τ range gate matching the C++ store)
**Scope:** add an opt-in per-module (per-channel) template / matched-filter
kernel path to two existing Python benchmark scripts, with automatic
fallback to the per-type aggregate when a channel has no usable per-channel
fit. Default behaviour (per-type aggregate) is preserved exactly.
**Hard constraints:**
- Python-only. No C++, no pybind11, no binding changes.
- Modify only the two existing scripts listed in §4; do NOT create new
  scripts. Do NOT touch `fit_pulse_template.py`.
- Downstream benchmark comparisons (Tasks A and B, §7–§8) are documented
  recipes only — this plan does NOT execute them.

## Research inputs used

- `python/bind_dec.cpp` lines 267–295 (`PulseTemplateStore.lookup` /
  `type_template` bindings — both return `py::none()` on miss), 318–332
  (`PulseTemplate` binding — exposes **only** `tau_r_ns`, `tau_f_ns`,
  `is_global`; Python constructor accepts `(tau_r_ns, tau_f_ns,
  is_global=False)`; **no `is_valid` predicate**), 559–582
  (`WaveAnalyzer.deconvolve(samples, wres, tmpl)` accepts any
  `PulseTemplate`).
- `prad2dec/src/PulseTemplateStore.cpp` lines 86–108 (JSON loader stores
  only `channel_id → module_type` map; **per-channel τ values from the
  JSON are intentionally discarded** — see comment at line 87–88),
  148–165 (`Lookup(roc, slot, ch)` returns the **per-type aggregate**
  for that channel's module type, not a per-channel template).
- `prad2dec/src/WaveAnalyzer.cpp` lines 806–821, 847–854 (deconv reads
  only `tmpl.tau_r_ns` / `tmpl.tau_f_ns` and reconstructs templates
  analytically on every LM iteration — `tmpl.grid` is **not** read;
  Python-constructed `PulseTemplate` works correctly), 1121 (`is_global`
  only toggles the returned deconv state between `Q_DECONV_APPLIED` and
  `Q_DECONV_FALLBACK_GLOBAL` — it has no effect on the fit itself).
- `prad2dec/include/Fadc250Data.h` line 272, 298 (`is_global` kept for
  backwards compatibility; set to `true` on every category aggregate).
- `analysis/pyscripts/output/pulse_templates_025308_h50_RawSum_max20modules.json`
  (schema — per-channel entries keyed by module name `"G129"`, `"W1234"`;
  each entry has `channel_id` as `"<roc>_<slot>_<channel>"`,
  `module_type`, `tau_r_ns.median`, `tau_f_ns.median`, `t0_ns.median`,
  `peak_amp_adc.median`, `chi2_per_dof`, `good_fit` bool,
  `n_pulses_used`; `_by_type` aggregate block alongside).
- `analysis/pyscripts/benchmark_matched_filter.py` lines 126–181
  (`_load_template_params` — currently reads only `_by_type[material]`),
  205–239 (`build_kernel(tau_r, tau_f, clk_ns, kernel_samples)`),
  449–579 (`_process_npz` event loop), 814–833 (CLI).
- `analysis/pyscripts/benchmark_hycal_timing_deconv.py` lines 199–201
  (`RunContext` deconv counters), 209–242 (`load_templates` — current
  per-type-only path), 248–315 (`extract_peaks` — per-channel loop with
  `store_by_type.get(mod.type.name)`), 697–782 (`write_summary_csv`
  preamble with `# n_peaks_deconv_*` lines).
- `analysis/pyscripts/pileup_generator.py` lines 60–75, 300–325
  (synthetic `.npz` schema already carries per-event `base_roc_tags`,
  `base_slots`, `base_chans`, `base_channel_ids` — so the matched-filter
  benchmark **already has (roc, slot, ch) per event** and can look up
  per-channel templates without any additional metadata work).
- `docs/plans/plan-benchmark-hycal-timing-deconv.md` (section
  numbering / preamble-line / counter conventions mirrored below).

---

## 1. Purpose

The user has verified that τ_r and τ_f vary substantially channel-to-channel
across HyCal (consistent with cable-length, PMT, and ADC-channel
differences). The current benchmarks always use the per-material aggregate
template, which is a known approximation. Both scripts need the same
treatment:

- `benchmark_matched_filter.py` builds one MF kernel per material and reuses
  it for every channel of that material. A per-channel kernel built from the
  channel's own fitted (τ_r, τ_f) should match the pulse shape better and
  therefore raise the matched-filter SNR / detection efficiency.
- `benchmark_hycal_timing_deconv.py` passes one per-material `PulseTemplate`
  to `wa.deconvolve` for every channel of that type. A per-channel template
  should lower the per-pulse χ²/dof and reduce `n_peaks_deconv_failed`.

The template JSON produced by `fit_pulse_template.py` already contains the
per-channel medians, so the data is in hand. This plan wires the Python
scripts to consume them, with graceful fallback for channels whose fit was
not reliable.

## 2. Inputs

### 2.1 Template file schema (already produced by `fit_pulse_template.py`)

Top-level JSON object. Keys that do **not** start with `_` are per-channel
entries (one per module); keys starting with `_` are metadata blocks.
Relevant fields:

- Per-channel entry (e.g. `"G129"`, `"W1234"`):
  - `channel_id` — string `"<roc>_<slot>_<channel>"`, e.g. `"4_14_12"`.
  - `module_type` — one of `"PbWO4"`, `"PbGlass"`, `"LMS"`, `"Veto"`.
  - `tau_r_ns.median` — float (ns), may be NaN if fit failed.
  - `tau_f_ns.median` — float (ns), may be NaN if fit failed.
  - `good_fit` — bool. **False** when the per-channel fit was unreliable
    (too few clean pulses, bad residuals, etc.) — these channels must
    fall back to the per-type aggregate.
  - `n_pulses_used`, `chi2_per_dof`, `t0_ns.median`, `peak_amp_adc.median`
    — informational, not consumed here.
- `_by_type` block:
  - Per-material entries (`"PbWO4"`, `"PbGlass"`, `"LMS"`, `"Veto"`) each
    with `tau_r_ns.median`, `tau_f_ns.median` — these are the aggregates
    used by the current scripts, and remain the default / fallback.

### 2.2 Why we read the JSON directly (not via `PulseTemplateStore.lookup`)

**Critical finding from research:** `PulseTemplateStore::Lookup(roc, slot,
ch)` in C++ does **not** return a per-channel template. It returns the
per-type aggregate for that channel's module type (see
`prad2dec/src/PulseTemplateStore.cpp` line 86–108 comment: *"per-channel
τ values are intentionally ignored"*). Both `lookup()` and
`type_template()` return the same object today.

Therefore, to actually consume per-channel τ_r / τ_f, both scripts must
read the JSON directly and construct `PulseTemplate(tau_r_ns, tau_f_ns)`
instances in Python. This is safe because:

- `WaveAnalyzer.deconvolve` reads only `tmpl.tau_r_ns` and `tmpl.tau_f_ns`
  and reconstructs templates analytically each LM iteration
  (`prad2dec/src/WaveAnalyzer.cpp` lines 806–821, 847); the
  `PulseTemplate::grid` field is never read. A Python-constructed
  `PulseTemplate(tr, tf)` with a zero-filled grid works correctly.
- The matched-filter kernel builder (`build_kernel` in
  `benchmark_matched_filter.py`) already takes raw `(tau_r, tau_f)`
  floats, not a `PulseTemplate` object, so it needs no bridging at all.

The alternative — extending the C++ loader to produce real per-channel
`PulseTemplate` objects — is explicitly out of scope (no C++ changes).

### 2.3 Inputs the scripts already have in hand

- `benchmark_matched_filter.py`: synthetic `.npz` events from
  `pileup_generator.py` already carry per-event `base_roc_tags`,
  `base_slots`, `base_chans`, `base_channel_ids` arrays (see
  `pileup_generator.py` §3.2 output schema). No metadata plumbing needed.
- `benchmark_hycal_timing_deconv.py`: event loop already iterates
  `(ri, s, c)` and resolves `mod = p.hycal.module_by_daq(crate, s, c)`
  (lines 1117–1130). `roc.tag` is `roc.tag` of the current ROC; the trio
  `(roc.tag, s, c)` is the lookup key.

## 3. Non-goals

- No C++ changes. No pybind11 binding changes. No `is_valid` predicate
  on `PulseTemplate` (does not exist today — see §10 Q1).
- No new scripts. No changes to `fit_pulse_template.py`.
- No changes to default script behaviour. The new flags default to **off**;
  without the flag, both scripts behave bit-for-bit as they do today.
- No implementation of Task A or Task B benchmark comparisons — those are
  documented recipes (§7, §8) for the user to run manually.
- No change to the C++ `PulseTemplateStore` call pattern inside
  `load_templates()` in the deconv benchmark — `store` and the per-type
  dict are still built the same way and still serve as the fallback.

## 4. Files to modify

| File | Scope |
|---|---|
| `analysis/pyscripts/benchmark_matched_filter.py` | Work item (a): per-module MF kernel + fallback. |
| `analysis/pyscripts/benchmark_hycal_timing_deconv.py` | Work item (c): per-module template + fallback. |

No other files are touched. No new files are created.

---

## 5. Work item (a) — per-module matched filter

### 5.1 Current behaviour

`_load_template_params(manifest, args.template)` (lines 126–181) reads
only the `_by_type[material]` block and returns a single `(tau_r, tau_f)`
pair. `main()` calls `build_kernel(tau_r, tau_f, clk_ns, kernel_samples)`
(line 859) to build one global MF kernel + kernel_norm. The same kernel
is applied to every event in `_process_npz` (line 528:
`matched_filter(pedsub, kernel, kernel_norm)`).

### 5.2 New behaviour

A new helper reads the full template JSON and returns a dict
`{(roc_tag, slot, chan): (tau_r_ns, tau_f_ns)}` built only from
per-channel entries that pass **all** of the following gates:

1. `good_fit` is strictly `True`.
2. `channel_id` is present and parses as `"<roc>_<slot>_<chan>"`.
3. `tau_r_ns.median` and `tau_f_ns.median` are both finite and positive.
4. **Both τ values are within the C++ deconvolver's accepted range**
   (`tau_r_min_ns ≤ tau_r ≤ tau_r_max_ns` and same for τ_f), using the
   bounds read from `WaveConfig.nnls_deconv` — identical to what
   `PulseTemplateStore::LoadFromFile` enforces at
   `prad2dec/src/PulseTemplateStore.cpp` lines 116–124. See §5.6 for
   why this check is mandatory.

Entries that fail any gate are **skipped silently**; the per-event
lookup then misses and the fallback (per-type aggregate) is used. The
helper therefore needs the `WaveConfig` (to read the τ bounds); the
signature becomes:

```python
def _load_per_channel_taus(
    path: Path,
    wcfg: "dec.WaveConfig",
) -> Dict[Tuple[int, int, int], Tuple[float, float]]: ...
```

The per-event path in `_process_npz` grows a cache
`kernel_cache: Dict[Tuple[int, int, int], Tuple[np.ndarray, float]]`
keyed by `(roc_tag, slot, chan)`. For each event:

```python
# only when per_module_mode is True
key = (int(base_roc_tags[i]), int(base_slots[i]), int(base_chans[i]))
entry = kernel_cache.get(key)
if entry is None:
    tr_tf = per_channel_taus.get(key)
    if tr_tf is not None:
        tr, tf = tr_tf
        k, kn = build_kernel(tr, tf, clk_ns, kernel_samples)
        entry = (k, kn, False)       # (kernel, kernel_norm, is_fallback)
        ctx_counters["per_module"] += 1
    else:
        entry = (fallback_kernel, fallback_kernel_norm, True)
        ctx_counters["fallback"] += 1
    kernel_cache[key] = entry
kernel_i, kernel_norm_i, _ = entry
mf_output = matched_filter(pedsub, kernel_i, kernel_norm_i)
```

When `--per-module-kernels` is **off**, this entire branch is bypassed and
the single global `(kernel, kernel_norm)` from `main()` is used exactly
as today.

Counters are incremented **once per channel the first time it is seen**
(cache miss), not per event. That matches the semantics the user asked
for: "channels where per-module kernel was used" / "channels where
fallback was triggered". Documented explicitly in the summary preamble.

### 5.3 CLI flag and help text

New argument added to the `argparse.ArgumentParser` block in `main()`
(after `--kernel-samples`, before `--mf-nsigma`):

```python
ap.add_argument(
    "--per-module-kernels", action="store_true",
    help=(
        "Build the matched-filter kernel per (roc,slot,chan) from the "
        "per-channel tau_r/tau_f medians in the template JSON, with "
        "fallback to the per-type aggregate for channels whose "
        "per-channel fit is missing or marked good_fit=false.  Default: "
        "off (single per-type-aggregate kernel for every event, exact "
        "current behaviour).  Requires a template JSON that contains "
        "per-channel entries, not just a '_by_type' block.  Kernels are "
        "cached per channel; two end-of-run counters "
        "(n_channels_per_module, n_channels_fallback) are added to "
        "mf_benchmark.json and the stdout summary."))
```

The help text is deliberately verbose because this is the only user-
facing documentation for the flag's semantics.

### 5.4 Caching strategy for per-channel kernels

Channels fire many times per run (one event per `.npz` row, hundreds of
thousands of rows across all cells × ratios). `build_kernel` is cheap
(one 30-element `np.exp` chain) but doing it per event would waste work.
A simple dict keyed by `(roc_tag, slot, chan)` is sufficient:

- Memory bound: ≤ ~1800 HyCal channels × ~500 bytes per kernel ≈ < 1 MB.
- Lookup: dict `.get` — O(1), amortised.
- Lifetime: cache is created once in `main()` and passed into
  `_process_npz` for every file. Cache survives across .npz files so
  channels that appear in multiple files only build once.

The cache stores the resolved kernel for both the per-module and fallback
cases (fallback entries point at the shared global `(fallback_kernel,
fallback_kernel_norm)` tuple — same NumPy array object, no copy, so
memory does not blow up with fallbacks).

**Rebuild on clk_ns correction.** `main()` already rebuilds the global
kernel if the DAQ-config `clk_ns` disagrees with the manifest value
(lines 870–878). The per-module cache must be built from the authoritative
`clk_ns` only, so cache construction happens **after** the DAQ-vs-manifest
clk_ns reconciliation. Implementation detail: `per_channel_taus` dict is
loaded once up front (just the τ dict, not kernels); the cache is a
separate dict populated lazily during `_process_npz` using the final
`clk_ns`.

### 5.5 Reporting counters

Three new counters attached to the `global_summary` dict and printed to
stdout:

- `n_channels_per_module` — distinct channels for which a per-channel
  kernel was built and used.
- `n_channels_fallback` — distinct channels for which the per-type
  aggregate kernel was used instead (per-channel fit missing or
  `good_fit=false` or τ out of `nnls_deconv` range).
- `n_perchannel_rejected_tau_range` — JSON-load-time counter (not
  per-event): per-channel entries that passed `good_fit=true` and
  finiteness but failed the τ-range gate (§5.6). Logged once at
  start-up; written into `mf_params` for auditability. Pure
  informational; a channel rejected here is identical to one that was
  never in the JSON — both route through the fallback counter.

Added to `global_summary` / `mf_params` in `main()`:

```python
if args.per_module_kernels:
    global_summary["n_channels_per_module"] = ctx_counters["per_module"]
    global_summary["n_channels_fallback"]   = ctx_counters["fallback"]
    mf_params["per_module_kernels"]              = True
    mf_params["n_perchannel_rejected_tau_range"] = n_rejected_tau_range
    print(f"[summary] per-module kernels: "
          f"{ctx_counters['per_module']} channels   "
          f"fallback (per-type aggregate): "
          f"{ctx_counters['fallback']} channels   "
          f"(rejected tau-range at load: "
          f"{n_rejected_tau_range} channels)", flush=True)
else:
    mf_params["per_module_kernels"] = False
```

When the flag is off, the two run-time counters are not added (keeps the
JSON schema backwards-compatible for existing consumers); only the
`per_module_kernels: false` boolean is written.

### 5.6 Per-channel τ-range gate (why it is mandatory)

`WaveAnalyzer.deconvolve` at `prad2dec/src/WaveAnalyzer.cpp` lines
847–854 applies its own τ-range gate against `WaveConfig.nnls_deconv`
and returns `Q_DECONV_BAD_TEMPLATE` if either τ is out of range. The
matched-filter script does not call `wa.deconvolve` so this failure
mode is a soft issue here — a per-channel kernel with absurd τ values
would still *build*, just with poor pulse match. However, the deconv
benchmark (§6) absolutely needs the gate: a per-channel entry with
positive-but-out-of-range τ (e.g. `tau_r = 0.1 ns` when
`tau_r_min_ns = 0.5`) would pass a naive finite/positive-only check,
reach `wa.deconvolve`, come back as `Q_DECONV_BAD_TEMPLATE`, and
increment `n_peaks_deconv_failed` — bypassing the fallback entirely.

Even though the matched-filter path tolerates out-of-range τ, we apply
the same gate in both scripts so that:

- the two scripts agree on which per-channel entries count as "usable"
  (lets us compare per-module coverage between runs);
- future refactoring that threads MF kernels through `wa.deconvolve`
  cannot regress silently;
- the semantics of `n_channels_per_module` match the semantics of
  `n_peaks_per_module_template` in §6.

Implementation: the τ bounds are read from `wcfg.nnls_deconv` where
`wcfg = _dec.WaveConfig(cfg.wave_cfg)` (already built at
`benchmark_matched_filter.py:193` and `benchmark_hycal_timing_deconv.py`
inside `_common.setup_pipeline`). The per-channel loader applies the
four comparisons inline during the JSON walk; rejection count is
returned alongside the dict.

---

## 6. Work item (c) — per-module template in deconv benchmark

### 6.1 Current behaviour

`load_templates(path, wave_cfg)` (lines 209–242) builds a
`PulseTemplateStore`, calls `store.type_template("PbWO4")` and
`store.type_template("PbGlass")`, and returns a dict `store_by_type`
keyed by module-type name. `extract_peaks(wa, store_by_type, cd, mod,
ctx)` (lines 248–315) resolves `tmpl = store_by_type.get(mod.type.name)`
and passes it to `wa.deconvolve(cd.samples, wres, tmpl)`.

Deconv health counters on `RunContext` (lines 199–201):
`n_peaks_deconv_per_channel` (`Q_DECONV_APPLIED`),
`n_peaks_deconv_global` (`Q_DECONV_FALLBACK_GLOBAL`),
`n_peaks_deconv_failed`.

### 6.2 New behaviour

Add a side-car dict `tmpl_by_channel: Dict[Tuple[int,int,int],
PulseTemplate]` populated from the same template JSON that
`load_templates` already consumes. Build it inside a new helper

```python
def load_per_channel_templates(
    path: Path,
    wave_cfg,                       # dec.WaveConfig; same type-tolerant
                                     # handling as load_templates
) -> Tuple[Dict[Tuple[int,int,int], "dec.PulseTemplate"], int]: ...
```

that:

1. Opens the JSON once with `json.load`. **Requires `import json` at
   module top (currently absent — the existing script has never needed
   `json` because `PulseTemplateStore` parses the JSON in C++).**
2. Normalises `wave_cfg` to `dec.WaveConfig` (same `isinstance` check as
   `load_templates`) and reads τ bounds:
   `tr_lo, tr_hi = wcfg.nnls_deconv.tau_r_min_ns, …max_ns`; same for τ_f.
3. Iterates top-level keys that do not start with `_`.
4. For each per-channel entry, **skip** if **any** of these fail:
   - `good_fit` is not strictly `True`;
   - `channel_id` missing or does not parse as `"<roc>_<slot>_<ch>"`;
   - `tau_r_ns.median` or `tau_f_ns.median` is missing / non-finite /
     non-positive;
   - **τ values are out of `nnls_deconv` range** (see §6.6) — count
     these in a dedicated local tally `n_rejected_tau_range` so the
     caller can report them separately from "entry not in JSON".
5. For survivors, construct `tmpl = dec.PulseTemplate(float(tr),
   float(tf), is_global=True)` and store under `(roc, slot, ch)`.
6. Return `(tmpl_by_channel, n_rejected_tau_range)`.

`is_global` is set to `True` so the returned deconv state remains
`Q_DECONV_FALLBACK_GLOBAL` and the existing three deconv-state counters
retain their current (approximate) meaning — see §6.4 for the reporting
tradeoff this creates. The Python-constructed template is valid for
`wa.deconvolve` because the C++ deconvolver reads only `tau_r_ns` /
`tau_f_ns` (§2.2).

`main()` grows:

```python
tmpl_by_channel: Dict[Tuple[int,int,int], object] = {}
n_rejected_tau_range = 0
if args.per_module_templates:
    tmpl_by_channel, n_rejected_tau_range = load_per_channel_templates(
        ctx.template_file, p.cfg.wave_cfg)
    print(f"[setup] per-module templates: "
          f"{len(tmpl_by_channel)} channels loaded "
          f"(from {ctx.template_file.name}; "
          f"{n_rejected_tau_range} rejected for tau out of "
          f"nnls_deconv range)", flush=True)
```

The event-loop call to `extract_peaks` grows one new argument:

```python
ch_peaks = extract_peaks(
    p.wave_ana, store_by_type, tmpl_by_channel,
    roc.tag, s, c,                # new: raw DAQ coords for lookup
    cd, mod, ctx)
```

`extract_peaks` is modified to prefer the per-channel template when
available, falling back to the per-type dict:

```python
def extract_peaks(wa, store_by_type, tmpl_by_channel,
                  roc_tag, slot, chan,
                  cd, mod, ctx):
    wres = wa.analyze_result(cd.samples)
    if wres.npeaks == 0:
        return []

    tmpl = tmpl_by_channel.get((int(roc_tag), int(slot), int(chan)))
    used_per_module = tmpl is not None
    if tmpl is None:
        tmpl = store_by_type.get(mod.type.name)   # per-type fallback

    out = (wa.deconvolve(cd.samples, wres, tmpl)
           if tmpl is not None else None)

    # ... existing used_deconv / counter logic unchanged ...

    # NEW: per-pulse per-module vs type-fallback counters
    if out is not None and out.n == wres.npeaks and used_deconv:
        if used_per_module:
            ctx.n_peaks_per_module_template += out.n
        else:
            ctx.n_peaks_type_fallback += out.n
    # (deconv-failed peaks are counted in n_peaks_deconv_failed as today;
    #  no need to split them further.)
```

**Null-template fallback semantics.** If a channel has neither a
per-channel entry (not in `tmpl_by_channel`) **nor** a per-type template
(e.g. unexpected module type), `tmpl` remains `None` and the existing
`pk.integral` fallback branch runs — behaviour unchanged from today.

### 6.3 CLI flag and help text

New argument added in the `argparse.ArgumentParser` block in `main()`:

```python
ap.add_argument(
    "--per-module-templates", action="store_true",
    help=(
        "Use the per-channel tau_r/tau_f medians from the template JSON "
        "for wa.deconvolve(), falling back to the per-type aggregate for "
        "any channel whose per-channel fit is missing or marked "
        "good_fit=false.  Default: off (per-type aggregate for every "
        "channel of a given material, exact current behaviour).  "
        "Requires a template JSON with per-channel entries, not just "
        "a '_by_type' block.  Adds two per-pulse counters "
        "(n_peaks_per_module_template, n_peaks_type_fallback) to the "
        "summary CSV preamble, orthogonal to the three existing "
        "n_peaks_deconv_* state counters."))
```

### 6.4 Compatibility with existing deconv-state counters

The three existing counters —
`n_peaks_deconv_per_channel` (= `Q_DECONV_APPLIED`),
`n_peaks_deconv_global` (= `Q_DECONV_FALLBACK_GLOBAL`),
`n_peaks_deconv_failed` — remain meaningful **and keep their current
semantics**. Reasoning:

- `wa.deconvolve` returns `Q_DECONV_FALLBACK_GLOBAL` whenever
  `tmpl.is_global == true` (`prad2dec/src/WaveAnalyzer.cpp` line 1121).
- The per-channel templates we build in §6.2 are constructed with
  `is_global=True` (because they are still category-level in the C++
  state model — the `is_global` flag describes the deconvolver's view
  of the template's aggregation level, not the Python layer's choice of
  per-channel vs. per-type).
- Therefore, when `--per-module-templates` is set **and** the deconv
  succeeds, the pulse continues to land in `n_peaks_deconv_global`
  exactly as it does today. The new per-module vs. fallback distinction
  lives in two **orthogonal** counters:
  - `n_peaks_per_module_template` — peaks where the input template came
    from the per-channel JSON entry.
  - `n_peaks_type_fallback` — peaks where the input template came from
    the per-type `_by_type` aggregate.

Setting `is_global=False` on the per-channel templates would re-route
successful deconvs into `n_peaks_deconv_per_channel`, which would silently
change the meaning of an existing counter and break any downstream
consumer that compares baseline vs. per-module runs on that number. The
chosen design keeps the state counters invariant across the two modes
and surfaces the per-module/fallback split in two new counters only.

**Deliberate reporting-compatibility tradeoff.** Marking a
Python-constructed per-channel template as `is_global=True` is a
*semantic abuse* of the flag as documented in `prad2dec/include/
Fadc250Data.h` lines 272–275 and `python/bind_dec.cpp` lines 319–322
(both state that `is_global=True` denotes a per-type aggregate). The
tradeoff is intentional and bounded:

- the three existing deconv-state counters
  (`n_peaks_deconv_per_channel`, `n_peaks_deconv_global`,
  `n_peaks_deconv_failed`) remain compatible across baseline and
  per-module runs, so direct `grep` / `diff` on summary CSVs continues
  to work;
- the cost is that `n_peaks_deconv_global` now counts **both**
  per-type aggregates **and** Python-constructed per-channel templates
  — the C++ state label is no longer a faithful description of the
  Python-level template source;
- the two **new** counters (`n_peaks_per_module_template`,
  `n_peaks_type_fallback`) are therefore the **authoritative** source
  for the per-channel-vs-per-type distinction. Any downstream analysis
  that cares about that split must consume the new counters, not the
  old state counters.

If a future change exposes a dedicated per-channel state constant
(e.g. `Q_DECONV_APPLIED_PER_CHANNEL`), the Python side should flip
`is_global=False` and the two new counters become redundant. Until
then, the two sets of counters together carry the full picture.

**Invariant (optional cross-check):** when the flag is on and the
deconv run succeeds, each peak increments exactly one of
{`n_peaks_per_module_template`, `n_peaks_type_fallback`} **and** exactly
one of {`n_peaks_deconv_per_channel`, `n_peaks_deconv_global`}. Sum of
the first pair should equal `n_peaks_deconv_per_channel +
n_peaks_deconv_global` within the run (modulo pulses where `tmpl` was
`None`, which contribute to neither new counter and increment
`n_peaks_deconv_failed`).

### 6.5 Reporting counters

Three new fields on `RunContext` (added to the `@dataclass` block near
lines 199–201):

```python
n_peaks_per_module_template:   int = 0
n_peaks_type_fallback:         int = 0
n_perchannel_rejected_tau_range: int = 0   # JSON-load-time, not per-event
```

Three new lines in the `write_summary_csv` preamble (added alongside the
existing `# n_peaks_deconv_*` lines near line 734):

```python
f"# n_peaks_per_module_template={ctx.n_peaks_per_module_template}",
f"# n_peaks_type_fallback={ctx.n_peaks_type_fallback}",
f"# n_perchannel_rejected_tau_range={ctx.n_perchannel_rejected_tau_range}",
```

Two new lines in the end-of-run `[deconv health]` stderr summary (around
line 1270):

```python
sys.stderr.write(
    f"[per-module]    per_module={ctx.n_peaks_per_module_template}  "
    f"type_fallback={ctx.n_peaks_type_fallback}  "
    f"(rejected tau-range at load: "
    f"{ctx.n_perchannel_rejected_tau_range} channels)\n"
)
```

Preamble is only emitted when `--per-module-templates` is on? **No** —
all three counters are always written so the CSV schema is stable across
runs. When the flag is off, all three are `0`.

`main()` populates `ctx.n_perchannel_rejected_tau_range` from the second
return value of `load_per_channel_templates` right after the loader call.

### 6.6 Per-channel τ-range gate

Same rationale as §5.6, with teeth here: `wa.deconvolve` at
`prad2dec/src/WaveAnalyzer.cpp` lines 847–854 returns
`Q_DECONV_BAD_TEMPLATE` whenever either τ falls outside
`nnls_deconv.tau_r_{min,max}_ns` or `nnls_deconv.tau_f_{min,max}_ns`.
Without a pre-filter, an out-of-range per-channel entry would:

1. populate `tmpl_by_channel`;
2. be preferred over the per-type fallback in `extract_peaks`;
3. come back from `wa.deconvolve` with `state ==
   Q_DECONV_BAD_TEMPLATE`;
4. miss both deconv-success branches and increment
   `n_peaks_deconv_failed` for every peak on that channel for the
   entire run — with no fallback ever attempted.

Applying the identical four comparisons that `PulseTemplateStore::
LoadFromFile` uses (`PulseTemplateStore.cpp` lines 116–124) at JSON-load
time keeps per-channel coverage honest: a channel whose per-channel τ
would be rejected by `wa.deconvolve` is excluded from `tmpl_by_channel`,
`extract_peaks` sees a cache miss, and the per-type fallback is used —
exactly as if the per-channel entry had `good_fit=false`. The rejection
count is surfaced separately (§6.5) so that an unusually noisy
`fit_pulse_template.py` run can be diagnosed from the counter rather
than from a mysterious drop in `n_peaks_per_module_template`.

---

## 7. Downstream Task A — deconv χ²/dof + amplitude recovery comparison (recipe only)

**Not implemented by this plan.** The user runs these two commands manually
on the same input, then diffs the two resulting `*_summary.csv` files.

### 7.1 Recipe

```bash
# Baseline (per-type aggregate, current default)
python3 analysis/pyscripts/benchmark_hycal_timing_deconv.py \
    /path/to/evio_dir \
    --template-file analysis/pyscripts/output/pulse_templates_025308_h50_RawSum_max20modules.json \
    --out-path output/deconv_run25320_baseline \
    --max-events 50000 \
    --run-num 25320

# Per-module (new flag)
python3 analysis/pyscripts/benchmark_hycal_timing_deconv.py \
    /path/to/evio_dir \
    --template-file analysis/pyscripts/output/pulse_templates_025308_h50_RawSum_max20modules.json \
    --out-path output/deconv_run25320_permodule \
    --max-events 50000 \
    --run-num 25320 \
    --per-module-templates
```

### 7.2 Comparison targets (read from each `*_summary.csv` preamble)

| Metric | Baseline | Per-module | Expectation |
|---|---|---|---|
| `# mean_chi2_per_dof=…` | X_base | X_pm | X_pm < X_base if per-module helps |
| `# n_peaks_deconv_failed=…` | F_base | F_pm | F_pm ≤ F_base |
| `# n_peaks_per_module_template=…` | 0 | P | — |
| `# n_peaks_type_fallback=…` | 0 | T | — |
| `# deconv_failed_fraction=…` | — | — | should drop or stay flat |

**Diagnostic sanity check.** If `T >> P` in the per-module run, most
channels fell back to the per-type aggregate — meaning the per-channel
fits were mostly `good_fit=false` for this detector. The two runs are
then effectively comparing baseline-vs-baseline on those channels and
no χ²/dof improvement should be expected. Record this ratio before
drawing conclusions.

## 8. Downstream Task B — σ_E/E comparison via `--calibrate` (recipe only)

**Not implemented by this plan.** Same two runs as Task A but add
`--calibrate` to both commands. Compare `sigma_over_E_CB` across the four
clustering configs × two template modes = 8 numbers.

### 8.1 Recipe

```bash
# Baseline + calibration
python3 analysis/pyscripts/benchmark_hycal_timing_deconv.py \
    /path/to/evio_dir \
    --template-file analysis/pyscripts/output/pulse_templates_025308_h50_RawSum_max20modules.json \
    --out-path output/deconv_run25320_baseline_cal \
    --max-events 50000 --run-num 25320 \
    --calibrate

# Per-module + calibration
python3 analysis/pyscripts/benchmark_hycal_timing_deconv.py \
    /path/to/evio_dir \
    --template-file analysis/pyscripts/output/pulse_templates_025308_h50_RawSum_max20modules.json \
    --out-path output/deconv_run25320_permodule_cal \
    --max-events 50000 --run-num 25320 \
    --calibrate \
    --per-module-templates
```

### 8.2 Comparison target

The `sigma_over_E_CB` column in each row of the summary CSV, for each of
`legacy_largestE`, `legacy_smallestDt`, `gated_largestE`,
`gated_smallestDt`. Headline question: does per-module template beat
per-type for **any** of the four clustering configs? If yes, which one(s)?

### 8.3 Both tasks are user-driven, not scripted

No new wrapper script is created. The user runs these commands manually
and inspects the summary CSVs by hand (or with the existing `diff` /
`csvkit` tooling they already use).

---

## 9. Validation (for both script extensions)

Mechanical validation only — no runs on full datasets.

1. **Byte-compile.** `python3 -m py_compile` on both modified scripts
   must succeed with zero output.
2. **CLI help.** `python3 analysis/pyscripts/benchmark_matched_filter.py
   --help` prints `--per-module-kernels`; `python3
   analysis/pyscripts/benchmark_hycal_timing_deconv.py --help` prints
   `--per-module-templates`. Each help block includes the "requires
   per-channel entries" and "fallback for good_fit=false" clauses from
   §5.3 / §6.3.
3. **Import check.** From `analysis/pyscripts/`:
   `python3 -c "import benchmark_matched_filter; import
   benchmark_hycal_timing_deconv; print('ok')"` must print `ok` without
   raising (deferred prad2py import keeps this safe even without
   bindings on the host).
4. **Smoke test — matched filter.** Point at a small synthetic dir
   (one `.npz` file is enough) and run with `--per-module-kernels`
   enabled, using a template file known to contain per-channel entries
   covering the test channels (e.g.
   `analysis/pyscripts/output/pulse_templates_025308_h50_RawSum_max20modules.json`).
   Confirm `n_channels_per_module` and `n_channels_fallback` both
   appear in `mf_benchmark.json["global"]`, that their sum is `≤` the
   number of distinct channels in that dir's manifest, **and that
   `n_channels_per_module > 0`**. The strict-greater check is
   non-negotiable: an all-fallback run (either a template file without
   per-channel entries, or a bug in the JSON walker that drops every
   entry) would otherwise silently pass validation while delivering
   zero real per-module coverage. The `mf_params.n_perchannel_rejected_tau_range`
   value should also appear and be sanity-checked (expect 0 or a small
   number — a large value signals a mis-matched `nnls_deconv` τ-range).
5. **Smoke test — deconv benchmark.** Run on one EVIO file with
   `--max-events 1000 --per-module-templates` and the same per-channel-
   covered template file. Confirm the output `*_summary.csv` preamble
   contains the three new preamble lines
   (`# n_peaks_per_module_template=…`, `# n_peaks_type_fallback=…`,
   `# n_perchannel_rejected_tau_range=…`) **and that
   `n_peaks_per_module_template > 0`**. Same silent-fallback trap as
   step 4.

Explicitly out of scope: no runs on the full 324 M-row TSV, no full
EVIO slices. Those are Task A / Task B territory (§7, §8) and belong
to the user.

---

## 10. Open questions

### Q1. `PulseTemplate.is_valid` or equivalent null-check predicate

**Resolved.** `PulseTemplate` exposes exactly three Python-visible
attributes: `tau_r_ns`, `tau_f_ns`, `is_global`
(`python/bind_dec.cpp` lines 330–332). There is **no `is_valid`
predicate** and no `__bool__`. The null-check must therefore be based on
the dict `.get()` returning `None`:

```python
tmpl = tmpl_by_channel.get((roc, slot, ch))   # None if not loaded
if tmpl is None:
    tmpl = store_by_type.get(mod.type.name)   # per-type fallback
```

Same pattern for `PulseTemplateStore.lookup` returning `py::none()`
(`python/bind_dec.cpp` line 272) — Python sees `None` and we test with
`is None`. No further researcher follow-up required.

### Q2. Should the per-channel JSON reader be shared across both scripts?

Both scripts need essentially the same JSON walker. Options:

- **(a) Duplicate.** Each script carries its own private helper
  (`_load_per_channel_taus` in `benchmark_matched_filter.py`,
  `load_per_channel_templates` in `benchmark_hycal_timing_deconv.py`).
  ~30 lines of duplication; no cross-file coupling introduced.
- **(b) Add to `_common.py`.** One shared helper
  `load_per_channel_template_entries(path)` returning the raw dict; each
  script wraps it with its own type-conversion (floats for MF, `dec.PulseTemplate`
  for deconv).

**Recommendation: (a).** Scope of this plan is "modify two existing
scripts, no new scripts" — adding to `_common.py` is technically allowed
but adds coupling for ~30 lines of code. If future work needs the helper
elsewhere, promote then.

Open for user to override.

### Q3. Should `--per-module-kernels` reject templates where `_by_type` is missing?

If the user passes `--per-module-kernels` and the template JSON has no
`_by_type` block (so the fallback is unavailable), channels with no
per-channel fit would have no kernel at all. Options:

- **(a) Hard-fail at startup** if `_by_type[material]` is missing — safe,
  explicit, aligns with how `_load_template_params` already fails today.
- **(b) Silently drop events for channels with no template** — bad,
  silent data loss.

**Recommendation: (a) — fail fast.** The existing
`_load_template_params` already raises `SystemExit` if the material is
not in `_by_type`; `--per-module-kernels` inherits that behaviour by
calling the existing function to resolve the fallback (τ_r, τ_f) before
building the kernel cache. No new error path needed.

### Q4. Should the plan emit a warning if per-channel coverage is low?

If the per-channel JSON covers (say) < 50 % of channels seen in the run,
the per-module benchmark is effectively a baseline-plus-noise run. Should
we print a `[WARN]` at end-of-run when
`n_channels_per_module / (n_channels_per_module + n_channels_fallback)
< 0.5`?

**Recommendation: no warning this iteration.** The counters are visible
in the summary; adding thresholds now would bake in a magic number. User
can inspect the ratio from the counters directly. Open for override.

---

## 11. Change log entry

- Add one entry to each script's module docstring under a "Changelog"
  or "History" section (if absent, add one):
  - `benchmark_matched_filter.py`:
    *"Added `--per-module-kernels` (default off): build the matched-
    filter kernel per (roc, slot, chan) from per-channel τ medians in
    the template JSON, with fallback to the per-type aggregate for
    channels whose per-channel fit is missing or marked
    `good_fit=false`. Two new counters
    (`n_channels_per_module`, `n_channels_fallback`) added to
    `mf_benchmark.json["global"]`."*
  - `benchmark_hycal_timing_deconv.py`:
    *"Added `--per-module-templates` (default off): pass per-channel
    `PulseTemplate(τ_r, τ_f)` built directly from the template JSON to
    `wa.deconvolve`, with fallback to the per-type aggregate for
    channels whose per-channel fit is missing or `good_fit=false`. Two
    new per-pulse counters (`n_peaks_per_module_template`,
    `n_peaks_type_fallback`) added to the summary CSV preamble; the
    three existing `n_peaks_deconv_*` state counters keep their current
    semantics (per-channel template is marked `is_global=True` so
    successful deconvs continue to land in `n_peaks_deconv_global`)."*

No separate `CHANGELOG.md` is updated — these scripts do not have one;
the docstring entry is the record of record.

---

## Todo

### Phase 0 — plan review (this file)

- [ ] User reviews §10 open questions and accepts recommendations (or
      annotates overrides).
- [ ] User confirms the "no new scripts, Python-only" interpretation of
      §2.2 (reading the JSON directly in Python, bypassing the C++
      loader for per-channel τ values).

### Phase 1 — `benchmark_matched_filter.py` edits

- [ ] Add module-level helper
      `_load_per_channel_taus(path: Path, wcfg: "dec.WaveConfig") ->
      Tuple[Dict[Tuple[int,int,int], Tuple[float, float]], int]`
      that walks the template JSON (`json` is already imported in this
      file), skips `_`-prefixed keys, requires `good_fit is True`,
      parses `channel_id`, validates τ medians finite and positive,
      **and applies the four `nnls_deconv` τ-range comparisons from
      §5.6**. Returns `(per_channel_taus_dict, n_rejected_tau_range)`.
- [ ] Add `--per-module-kernels` to `argparse` in `main()` with the help
      text from §5.3.
- [ ] In `main()`, after the DAQ-vs-manifest clk_ns reconciliation,
      load per-channel τ dict when the flag is set (passing `wcfg`
      built at line 193); build an empty kernel cache
      `Dict[Tuple[int,int,int], Tuple[np.ndarray, float, bool]]`.
- [ ] Add two counters `ctx_counters = {"per_module": 0, "fallback": 0}`
      scoped to `main()` and passed into `_process_npz`.
- [ ] Extend `_process_npz` signature with
      `per_module_mode: bool`, `per_channel_taus: Optional[Dict]`,
      `kernel_cache: Optional[Dict]`, `ctx_counters: Optional[Dict]`,
      `fallback_kernel: np.ndarray`, `fallback_kernel_norm: float`.
      When `per_module_mode` is false, use the existing single-kernel
      path untouched.
- [ ] When `per_module_mode` is true, load per-event
      `(roc_tag, slot, chan)` from `base_roc_tags[i]`, `base_slots[i]`,
      `base_chans[i]` arrays; look up / lazily build the per-channel
      kernel; increment the appropriate counter on cache miss.
- [ ] Extend `global_summary` dict in `main()` with the two run-time
      counters (only when flag is set) and `mf_params` with
      `per_module_kernels` (always, bool) and
      `n_perchannel_rejected_tau_range` (only when flag is set).
- [ ] Add stdout `[summary] per-module kernels: …  (rejected tau-range
      at load: …)` print when flag is set, per §5.5.
- [ ] Add Changelog entry at top of module docstring (§11).

### Phase 2 — `benchmark_hycal_timing_deconv.py` edits

- [ ] **Add `import json` to the import block** (currently absent from
      the file — the existing script delegates JSON parsing to
      `PulseTemplateStore` in C++, so `json` has never been imported).
      Insert alphabetically between `import hashlib` and `import math`
      at the current line 43–44.
- [ ] Add module-level helper
      `load_per_channel_templates(path: Path, wave_cfg) ->
      Tuple[Dict[Tuple[int,int,int], dec.PulseTemplate], int]`
      matching §6.2 filter rules; applies the four `nnls_deconv`
      τ-range comparisons from §6.6; constructs `dec.PulseTemplate(tr,
      tf, is_global=True)` per surviving entry. Returns
      `(tmpl_by_channel, n_rejected_tau_range)`.
- [ ] Add `--per-module-templates` to `argparse` in `main()` with the
      help text from §6.3.
- [ ] Add three fields to `RunContext` dataclass:
      `n_peaks_per_module_template: int = 0`,
      `n_peaks_type_fallback: int = 0`,
      `n_perchannel_rejected_tau_range: int = 0`.
- [ ] In `main()`, build `tmpl_by_channel` only when flag is set (pass
      `p.cfg.wave_cfg`), else leave empty dict and
      `n_rejected_tau_range = 0`. Assign the returned rejection count
      to `ctx.n_perchannel_rejected_tau_range`. Log `[setup] per-module
      templates: N channels loaded (… rejected for tau out of
      nnls_deconv range)`.
- [ ] Modify `extract_peaks` signature to accept
      `tmpl_by_channel`, `roc_tag`, `slot`, `chan`; resolve per-channel
      template first, fall back to `store_by_type.get(mod.type.name)`;
      set a local `used_per_module` bool.
- [ ] Inside `extract_peaks`, after the existing `used_deconv`
      block, increment `ctx.n_peaks_per_module_template` or
      `ctx.n_peaks_type_fallback` per pulse when `used_deconv` and
      `out.n == wres.npeaks`.
- [ ] Update the event-loop call site (around line 1144–1145) to pass
      the new arguments: `extract_peaks(p.wave_ana, store_by_type,
      tmpl_by_channel, roc.tag, s, c, cd, mod, ctx)`.
- [ ] Add three `# n_peaks_*` / `# n_perchannel_*` lines to the CSV
      preamble in `write_summary_csv` (always written; §6.5).
- [ ] Add `[per-module] per_module=… type_fallback=… (rejected
      tau-range at load: …)` line to the end-of-run stderr summary.
- [ ] Add Changelog entry at top of module docstring (§11).

### Phase 3 — Validation

- [ ] `python3 -m py_compile analysis/pyscripts/benchmark_matched_filter.py
      analysis/pyscripts/benchmark_hycal_timing_deconv.py` returns 0.
- [ ] `python3 analysis/pyscripts/benchmark_matched_filter.py --help`
      shows `--per-module-kernels` with the §5.3 wording.
- [ ] `python3 analysis/pyscripts/benchmark_hycal_timing_deconv.py --help`
      shows `--per-module-templates` with the §6.3 wording.
- [ ] From `analysis/pyscripts/`:
      `python3 -c "import benchmark_matched_filter,
      benchmark_hycal_timing_deconv; print('ok')"` → prints `ok`.
- [ ] Smoke-test MF: tiny synth dir with a per-channel-covered template
      (e.g. `pulse_templates_025308_h50_RawSum_max20modules.json`);
      confirm both counters appear in `mf_benchmark.json["global"]`
      **and `n_channels_per_module > 0`** (strict). Also confirm
      `mf_params.n_perchannel_rejected_tau_range` is present.
- [ ] Smoke-test deconv: `--max-events 1000`, same per-channel-covered
      template; confirm all three new preamble lines appear and
      **`n_peaks_per_module_template > 0`** (strict).

### Phase 4 — Hand-off to user for Tasks A / B

- [ ] User runs the two Task A commands from §7.1 on their chosen EVIO
      slice and compares the preamble metrics (§7.2).
- [ ] User runs the two Task B commands from §8.1 and compares
      `sigma_over_E_CB` across the 4 configs × 2 modes.
- [ ] (Optional) User reports back whether per-module helps for any
      config; results may motivate a follow-up plan.
