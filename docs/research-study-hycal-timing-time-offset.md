# Research — Are `study_hycal_timing.py` TSV times corrected for `Module::time_offset`?

**Scope.** Resolve a disagreement between two prior investigations about whether the `seed_time`, `neighbor_time`, and `dt` columns written by `analysis/pyscripts/study_hycal_timing.py` are post-calibration (per-module `Module::time_offset` already subtracted) or raw peak times out of `WaveAnalyzer`. The TSV is produced by tracing the `pk.time` value from `WaveAnalyzer::Analyze` through the pybind11 binding, through `HyCalCluster::AddHit` and `CollectNeighborTiming`, back into Python via the `SeedNeighborTiming` struct, and into the TSV row writer.

This note walks that path end-to-end and quotes every site where `time_offset` could in principle be applied.

---

## Answer

**The TSV columns `seed_time`, `neighbor_time`, and `dt` are RAW. `Module::time_offset` is loaded upstream by `PipelineBuilder` into `fdec::Module::time_offset`, but it is never subtracted anywhere on the path `WaveAnalyzer → add_hit → CollectNeighborTiming → SeedNeighborTiming → TSV`. The script passes `pk.time` directly to `HyCalCluster::AddHit`, and `AddHit` stores it unchanged. High confidence.**

The prior investigation recorded in `docs/plans/plan-module-dt-diagnostic.md` §2.a is correct: the Python-side offset application the plan describes is needed.

---

## Evidence trail

### 1. `WaveAnalyzer::Analyze` — raw peak time, with an optional per-call offset that defaults to zero

`prad2dec/include/WaveAnalyzer.h:165`
```cpp
void Analyze(const uint16_t *samples, int nsamples, WaveResult &result, float time_offset = 0.f) const;
```

The analyzer has a `time_offset` parameter, but its default is `0.f`. Peak times are set in `findPeaks` at `prad2dec/src/WaveAnalyzer.cpp:587` as the pickoff sample converted to nanoseconds:

```cpp
p.time     = t_pickoff * 1e3f / cfg.clk_mhz;  // ns
```

(where `t_pickoff` is a sub-sample position inside the waveform, in sample units). This is a *waveform-local* time in ns from sample 0 — it has no notion of a module or a calibration offset.

The caller-supplied `time_offset` is applied at the very end of `Analyze` in `prad2dec/src/WaveAnalyzer.cpp:678-684`:

```cpp
// ── Time offset: add the caller's offset to each peak's time.  This
// is a convenience for the caller to avoid having to loop over the
// peaks themselves. The offset is in ns, will add the value when you
// call Analyze for each waveform.
for (int p = 0; p < result.npeaks; ++p)
    result.peaks[p].time -= time_offset;
```

So Analyze *can* subtract a per-call offset, but only if the caller passes one. The comment above the loop is slightly misleading ("add the caller's offset") — the arithmetic in the code is subtraction, which matches the `calib_time = raw_time − mod.time_offset` convention documented in `prad2det/include/HyCalTimeCalib.h:8`.

### 2. pybind11 `WaveAnalyzer.analyze` binding — passes no `time_offset` (uses the C++ default = 0)

`python/bind_dec.cpp:407-422`
```cpp
.def("analyze",
    [](const fdec::WaveAnalyzer &self, py::array_t<uint16_t> samples) {
        py::buffer_info buf = samples.request();
        if (buf.ndim != 1)
            throw py::value_error("samples must be a 1-D uint16 array");
        fdec::WaveResult res;
        {
            py::gil_scoped_release rel;
            self.Analyze(static_cast<const uint16_t*>(buf.ptr),
                         static_cast<int>(buf.shape[0]), res);
        }
        py::list peaks;
        for (int i = 0; i < res.npeaks; ++i)
            peaks.append(res.peaks[i]);
        return py::make_tuple(res.ped.mean, res.ped.rms, peaks);
    },
    py::arg("samples"),
    "Analyze one channel's FADC waveform (uint16 numpy array). ...")
```

The lambda calls `self.Analyze(ptr, n, res)` with three arguments only — no `time_offset`. The C++ default `0.f` applies, so **the peak times returned to Python are in waveform-local ns from sample 0, with zero offset subtracted**. This is also true for the sibling `.def("analyze_full", ...)` at `bind_dec.cpp:428-448` (same signature, same three-argument call).

Nothing on the Python side ever opens a second path into `Analyze`; `WaveAnalyzer.analyze` and `WaveAnalyzer.analyze_full` are the only exposed entry points that run the full pipeline.

### 3. `study_hycal_timing.py` loop — pass-through, no offset applied in Python

`analysis/pyscripts/study_hycal_timing.py:173-181`
```python
_, _, peaks = p.wave_ana.analyze(cd.samples)
for pk in peaks:
    if pre_lo is not None and pk.time <= pre_lo:
        continue
    if pre_hi is not None and pk.time >= pre_hi:
        continue
    p.hc_clusterer.add_hit(mod.index,
                            mod.energize(pk.integral),
                            float(pk.time))
```

The script compares `pk.time` against the (optional) pre-window and otherwise forwards `float(pk.time)` directly into `add_hit`. No arithmetic touching `mod.time_offset` or any loaded calibration dict appears anywhere in the file (`grep "time_offset|time_calib" analysis/pyscripts/` returns zero matches across the whole directory). The pre-window ("`--pre-window LO HI`") is a trigger-latency gate, not a calibration — it is applied symmetrically to `pk.time` and would need the same thought experiment regardless of whether the time was calibrated.

Also note: `fdec::Module` is bound in `python/bind_det.cpp:531-572` and does **not** expose a `.time_offset` attribute at all (compare the exhaustive list of `.def_readonly`s there — `cal_factor`, `cal_base_energy`, `cal_non_linear_1/2` are exposed; `time_offset` is omitted). So even if the Python script wanted to read `mod.time_offset`, it currently cannot.

### 4. `HyCalCluster::AddHit` — stores the incoming `time` verbatim

`prad2det/src/HyCalCluster.cpp:72-78`
```cpp
void HyCalCluster::AddHit(int module_index, float energy, float time)
{
    if (module_index < 0 || module_index >= sys_.module_count()) return;
    if (test_bit(sys_.module(module_index).flag, kDeadModule)) return;
    if (energy > config_.min_module_energy)
        hits_.push_back({module_index, energy, time});
}
```

`AddHit` does not look at `sys_.module(module_index).time_offset`. The `time` parameter is pushed straight into `hits_` as `ModuleHit{module_index, energy, time}`. A `grep "time_offset"` across the entire `prad2det/src/HyCalCluster.cpp` returns zero matches — nothing in the clusterer ever consults the per-module offset.

### 5. `HyCalCluster::CollectNeighborTiming` — differences stored `hits_[k].time`

`prad2det/src/HyCalCluster.cpp:770,787-789`
```cpp
const float st = hits_[seed].time;
...
s.seed_time       = st;
s.neighbor_time   = hits_[k].time;
s.dt              = hits_[k].time - st;
```

Both `st` and `hits_[k].time` come straight from the `hits_` buffer populated by `AddHit`. Because `AddHit` stored the raw `pk.time` without touching `time_offset` (step 4), `seed_time`, `neighbor_time`, and `dt` here are all in the same uncorrected sample-based ns units.

Nothing in `CollectNeighborTiming` subtracts `module.time_offset` either (confirmed by the `grep "time_offset"` in step 4).

### 6. pybind11 return path for `SeedNeighborTiming` — no transformation

`python/bind_det.cpp:773-785` binds the struct with pure `.def_readonly`s:
```cpp
py::class_<fdec::HyCalCluster::SeedNeighborTiming>(m, "SeedNeighborTiming",
    ...
    .def_readonly("seed_module",     &fdec::HyCalCluster::SeedNeighborTiming::seed_module)
    .def_readonly("neighbor_module", &fdec::HyCalCluster::SeedNeighborTiming::neighbor_module)
    .def_readonly("seed_time",       &fdec::HyCalCluster::SeedNeighborTiming::seed_time)
    .def_readonly("neighbor_time",   &fdec::HyCalCluster::SeedNeighborTiming::neighbor_time)
    .def_readonly("dt",              &fdec::HyCalCluster::SeedNeighborTiming::dt)
    ...);
```

And `collect_neighbor_timing` at `bind_det.cpp:835-848` is a straight wrapper around `CollectNeighborTiming` that returns the vector of rows unmodified:
```cpp
.def("collect_neighbor_timing",
    [](const fdec::HyCalCluster &self, double max_quantized_dist) {
        std::vector<fdec::HyCalCluster::SeedNeighborTiming> out;
        { py::gil_scoped_release rel;
          self.CollectNeighborTiming(out, max_quantized_dist); }
        return out;
    }, ...);
```

The script writes those fields verbatim at `study_hycal_timing.py:197-205` with only format specifiers applied — no arithmetic. So the TSV field values are byte-identical to the C++ `SeedNeighborTiming::seed_time` / `neighbor_time` / `dt`.

### 7. Cross-check against analysis tools that *do* apply `time_offset` explicitly

Three C++ analysis tools demonstrate the project's convention of **the caller doing the subtraction before `AddHit` / before differencing**. This is strong independent evidence that `AddHit` and the clusterer do not already apply it (otherwise these tools would double-subtract).

**`analysis/tools/time_recon_check.cpp:362,370`** — fetches `time_offset` from the module and subtracts it from the raw peak time **before** calling `AddHit`:
```cpp
// timing offset for this module
float time_offset = mod->time_offset;
...
for (int p = 0; p < ev->npeaks[j]; ++p) {
    float peak_time = ev->peak_time[j][p] - time_offset; // apply module time offset
    if (peak_time <= hc_win.lo) continue;
    if (peak_time >= hc_win.hi) continue;
    float adc = ev->peak_integral[j][p] * gain;
    float energy = static_cast<float>(mod->energize(adc));
    clusterer.AddHit(mod->index, energy, peak_time);
}
```

**`analysis/tools/time_recon_check.cpp:414,429`** — same pattern a second time in the neighbor-dt loop (the "`dt = peak_time - hits[0].time`" differences). This also tells us that `hits[0].time` coming out of `ReconstructHits` is likewise *not* offset-corrected downstream (the subtraction is intentionally applied only on the per-peak side here).

**`analysis/tools/hycal_module_time_calib.cpp:274-275,286`**:
```cpp
float seed_time = hits[0].time;
if (validation) seed_time -= seed_mod->time_offset;
...
if (validation) event_times[mod->index] = ev->peak_time[j][0] - mod->time_offset;
else event_times[mod->index] = ev->peak_time[j][0];
```

The validation branch subtracts `time_offset` on both the seed cluster time and every per-peak time; the non-validation branch (building fresh calibrations) does not. The fact that the subtraction is gated on a boolean here confirms there is no implicit subtraction happening elsewhere.

**`analysis/tools/hycal_shower_profile.cpp` (lines 234, 242, 288, 306, 326)** — matches the same convention (`peak_time = ev->peak_time[j][p] - mod->time_offset` before any use). I spot-checked these to make sure the pattern is consistent; it is.

The three tools together confirm: **the convention project-wide is `calib_time = raw_time − mod.time_offset`, and that subtraction is the caller's responsibility.** Neither `WaveAnalyzer`, the Python binding, `AddHit`, nor `CollectNeighborTiming` subtracts it. The `WaveAnalyzer::Analyze` optional `time_offset` parameter is the only "built-in" path, and it is never invoked from Python.

### 8. `PipelineBuilder` loads the offsets but does not apply them to the `study_hycal_timing.py` data path

`prad2det/src/PipelineBuilder.cpp:357-389` loads `database/hycal_time_offsets/<file>.json` (via `prad2::LoadHyCalTimeCalib`) into `fdec::Module::time_offset`. The comment at `:358-359` is explicit:
```cpp
// Applies offsets directly into HyCalSystem::Module::time_offset so
// downstream code can use calib_time = raw_time - mod.time_offset.
```

So the pre-flight observation ("`PipelineBuilder` loads `time_calib.json`") is factually correct, but it is the *loading* step. The *application* step is the responsibility of the consumer, and the consumer in question — `study_hycal_timing.py` — never performs it (step 3) because:

1. The Python binding of `fdec::Module` does not expose `time_offset` (`python/bind_det.cpp:531-572`), so the script couldn't read it from the module even if it tried.
2. The Python binding of `WaveAnalyzer::Analyze` ignores the `time_offset` parameter (step 2).
3. The script itself does not implement any offset arithmetic (step 3).

The two prior investigations are therefore not in genuine conflict. Both facts are true independently:
- Offsets ARE loaded by `PipelineBuilder` (pre-flight was right).
- Offsets are NOT applied along the `study_hycal_timing.py` data path (the plan was right).

---

## Implication for the diagnostic script

The diagnostic plan at `docs/plans/plan-module-dt-diagnostic.md` §2.a — which calls for Python-side application of `time_offset` via a `--time-calib` option — is correct in substance. Specific confirmations:

- **JSON schema expected by the plan.** The file `database/hycal_time_offsets/time_calib.json` is consumed by `prad2::LoadHyCalTimeCalib` in `prad2det/include/HyCalTimeCalib.h:72-135`. The accepted schema is:
  ```json
  {
    "default": 0.0,
    "modules": [
      {"name": "W735", "offset_ns": -0.74},
      {"name": "W736", "offset_ns":  1.21}
    ]
  }
  ```
  Keyed primarily by `"name"` (`HyCalTimeCalib.h:53`), with `"module_id"` (PrimEx id) accepted as a fallback (`:117-119`). Any Python-side reader should mirror this precedence to stay consistent with the C++ loader. Entries without a matching module are silently skipped (`:120-122`).

- **Where to apply the subtraction in the diagnostic.** The cleanest insertion point is in `study_hycal_timing.py` right before line 179's `add_hit` call: read a per-module offset dict keyed on either `mod.name` or `mod.id`, and change the call to something like:
  ```python
  t_calib = float(pk.time) - offsets.get(mod.name, offsets_default)
  p.hc_clusterer.add_hit(mod.index, mod.energize(pk.integral), t_calib)
  ```
  Doing it before `add_hit` (rather than post-processing the TSV) means the `dt` column, the pre-window check, and `collect_neighbor_timing`'s own seed-claim window (5 ns, `HyCalCluster.cpp:759`) all see calibrated times — matching how `time_recon_check.cpp:370` structures its pipeline.

- **Alternative, equivalent fix inside the C++ path (optional, not required by the plan).** Expose `fdec::Module::time_offset` as a read-only `py::property` in `bind_det.cpp` and either (a) have the script loop and subtract it in Python using `p.hycal.module(mod.index).time_offset`, or (b) start passing it through `WaveAnalyzer::Analyze`'s existing `time_offset` parameter by extending the Python binding of `.def("analyze", ...)` to accept an optional `time_offset` kwarg. Any of these would eliminate the need for the script to re-parse `time_calib.json` on its own.

- **Residual ambiguity.** None material. The one footnote worth logging: the comment in `WaveAnalyzer.cpp:678-681` says the loop "adds" the offset while the code subtracts it. That's a comment bug, not a code bug — the subtraction matches the documented `calib_time = raw_time − mod.time_offset` convention from `HyCalTimeCalib.h:8` and the three analysis tools in step 7. Any future fix should correct the comment, not the arithmetic.

---

## Confidence and residual ambiguity

**Confidence: high.** Everything above is confirmed by direct reading of source:

- `WaveAnalyzer::Analyze` default `time_offset = 0.f` — confirmed in header (`WaveAnalyzer.h:165`) and `.cpp` offset-subtraction loop (`WaveAnalyzer.cpp:682-683`).
- The Python binding calls `Analyze` with three args, no `time_offset` — confirmed by direct read of `bind_dec.cpp:407-422` and `:428-448`.
- `HyCalCluster::AddHit` stores `time` verbatim — confirmed at `HyCalCluster.cpp:72-78`; corroborated by a project-wide `grep "time_offset"` returning zero matches in `HyCalCluster.cpp`.
- `CollectNeighborTiming` differences stored `hits_[k].time` with no additional arithmetic — confirmed at `HyCalCluster.cpp:770,787-789`.
- `fdec::Module` has no `time_offset` Python attribute — confirmed by exhaustive scan of the `py::class_<fdec::Module>` block at `bind_det.cpp:531-572`.
- `study_hycal_timing.py` and no other file under `analysis/pyscripts/` mentions `time_offset` or `time_calib` — confirmed by `grep -r` across the directory.
- The convention in sibling C++ analysis tools (`time_recon_check.cpp`, `hycal_module_time_calib.cpp`, `hycal_shower_profile.cpp`) is for the caller to subtract `mod->time_offset` explicitly before `AddHit` / before differencing — confirmed by direct read of the cited lines.
- `PipelineBuilder` does load the JSON into `Module::time_offset` — confirmed at `PipelineBuilder.cpp:357-389`.

The only inferences (not direct code reads) are: (a) the three cited analysis tools would be double-subtracting if `AddHit` applied the offset, so therefore `AddHit` must not apply it — this is airtight given the explicit subtractions on the raw-peak side; and (b) `hycal_shower_profile.cpp` lines were spot-checked rather than quoted in full here. Both inferences are low-risk.

No partial / conditional cases: the behaviour does not depend on run config, recon config, or any JSON file content. The offset is loaded if the file exists (per `PipelineBuilder.cpp:362-391`), but regardless of whether it is loaded, the `study_hycal_timing.py` TSV columns do not reflect it.
