# HyCal Time Calibration Explained (Code, Method, and Offset Usage)

This document focuses on two implementations in the current repository:

- `prad2det/include/HyCalTimeCalib.h`
- `analysis/tools/hycal_module_time_calib.cpp`

It answers four core questions:

1. What this HyCal time calibration is solving.
2. The mathematical principle behind it.
3. How the code computes per-module offsets step by step.
4. How to apply each module offset in replay and other analysis tools.

---

## 1. Problem Definition

Each HyCal crystal module has a quasi-constant time-zero bias (from electronics chain differences, channel-dependent delays, etc.).

If you compare events using raw pulse time directly (`peak_time` or `cl_time`), modules show systematic shifts. This leads to:

- Neighbor-pair time-difference distributions not centered at 0.
- Worse cluster-level timing resolution.
- May group noise/background pulse hits into clusters.
- Less stable physics selection based on timing windows.

So the target is to determine one constant per module, `offset_i`, and use:

$$
 t_{calib} = t_{raw} - offset_i
$$

This is the core purpose of the calibration.

---

## 2. Core Idea (Graph Model + Constrained Solve)

The key idea in `hycal_module_time_calib.cpp` is:

- Treat each module as a graph node.
- For each neighboring module pair $(i,j)$, measure:
  $$
  \Delta t_{ij} = t_i - t_j
  $$
- Accumulate events, fit each pair histogram with a Gaussian, and get mean $\mu_{ij}$ (and uncertainty).
- Convert to linear relations:
  $$
  offset_i - offset_j \approx \mu_{ij}
  $$
- Choose one reference module `ref`, set `offset_ref = 0`, and solve all other offsets with weighted least squares.

### 2.1 Why This Equation Holds

The code fills uncorrected (raw) time differences:

$$
\Delta t_{ij}^{raw} = t_i^{raw} - t_j^{raw}
$$

If the underlying physical time is approximately the same:

$$
 t_i^{raw} = t_{phys} + offset_i,\quad
 t_j^{raw} = t_{phys} + offset_j
$$

Then:

$$
\Delta t_{ij}^{raw} \approx offset_i - offset_j
$$

That is exactly the equation being solved.

---

## 3. `HyCalTimeCalib.h`: Loading and Injecting Offset Tables

`LoadHyCalTimeCalib(path, hycal, def_off)` behaves in four stages:

1. **Initialize all modules to a default**
   - Set every `hycal.module(i).time_offset = def_off`.

2. **Read JSON**
   - If file open or JSON parse fails, keep uniform defaults and print a warning.

3. **Apply file-level `default`**
   - If JSON has a numeric `default`, overwrite all modules with that value.

4. **Apply per-module overrides**
   - Iterate `modules` array, resolve module by `name` (preferred) or `module_id`.
   - Set that module’s `time_offset = offset_ns`.
   - Unknown modules are counted and skipped.

It returns a summary: final default and number of overrides.

### 3.1 JSON Format

Recommended format:

```json
{
  "default": 0.0,
  "modules": [
    {"name": "W735", "offset_ns": -0.74},
    {"name": "W736", "offset_ns": 1.21}
  ]
}
```

Resolution order per module:

- If listed in `modules`: use that `offset_ns`.
- Else if file has `default`: use `default`.
- Else: use caller-provided `def_off` (typically 0.0).

---

## 4. `hycal_module_time_calib.cpp`: Step-by-Step Computation

## 4.1 Input and Initialization

- Collect input `_raw.root` files.
- Build `TChain("events")` and bind `RawEventData` branches.
- Load run-dependent energy calibration and initialize HyCal.
- Load an existing time-offset file:
  - If `-V` is not given, default is `database/hycal_time_offsets/test.json`.
  - If `-V` is given, use that file.
- `-v` enables validation mode (explained in 4.4).

## 4.2 Build Histograms for Neighbor Pairs

Only PWO4 modules are considered.

For each module and each PWO4 neighbor, define a unique pair key:

$$
key = (\min(id_i,id_j),\max(id_i,id_j))
$$

Create one `TH1F` per pair (`[-10, 10] ns`, 200 bins) to accumulate per-event pair time differences.

## 4.3 Event Selection and Hit Reconstruction

Per-event flow:

1. Trigger and basic quality cuts:
   - Require `TBIT_sum`.
   - Require `nch <= 70`.

2. Build module hits:
   - Keep PWO4 modules only.
   - In run-config time window, pick best peak (largest integral).
   - Build cluster input with gain factor and energy calibration.

3. Form clusters and select single-cluster Mott-like events:
   - `hits.size() == 1`
   - `nblocks > 3`
   - `energy > 0.8 * Ebeam`

This provides a clean, timing-consistent sample for pair fitting.

## 4.4 Filling Pair $\Delta t$: Validation vs Non-Validation

The code builds `event_times[module_index]` as:

- **Default mode** (non-validation):
  $$
  t_i = t_i^{raw}
  $$
- **Validation mode** (`-v`):
  $$
  t_i = t_i^{raw} - offset_i
  $$
  where `offset_i` is `mod->time_offset` loaded from JSON.

Then for each neighboring pair in the same event, fill:

$$
\Delta t_{ij} = t_i - t_j
$$

This is the only place where “apply existing offset or not” is toggled.

- Default mode: derive new offsets.
- Validation mode: verify whether an existing offset table narrows pair distributions and pulls means toward 0.

## 4.5 Gaussian Fit Per Pair

For each pair histogram:

- Require `entries > 100`.
- Try fit near peak first; fallback to full range if needed.
- Record:
  - `mean = mu` (equation RHS)
  - `sigma`
  - `mean_error` (with floor `0.01 ns` to avoid infinite weight)
  - `chi2/ndf`

Only valid fitted pairs are used in the network solve.

## 4.6 Connected Graph and Reference Module

- Build adjacency graph from valid pairs.
- Current hard-coded reference: `reference_module = 1495`.
- BFS from reference defines `connected` component.

Only modules in this connected component can be solved in the same relative-offset system.

## 4.7 Weighted Least-Squares Solve

Unknowns: offsets of all non-reference modules in the connected component.

Each valid pair contributes:

$$
offset_{id1} - offset_{id2} = \mu
$$

with weight:

$$
w = \frac{1}{\sigma_\mu^2}
$$

where $\sigma_\mu$ comes from `mean_error` (floored at 0.01 ns).

Build normal equations `normal * x = rhs`, then solve with `TDecompSVD`:

- `solution(i)` gives module offset.
- If covariance is available, `sqrt(cov(i,i))` is `offset_error_ns`.

## 4.8 Outputs: JSON and ROOT

### JSON Output

One entry per PWO4 module, including:

- `name`, `module_id`
- `offset_ns`, `offset_error_ns`
- `connected` (in solved network)
- `resolved` (successfully solved)

Unresolved modules are currently written as `offset_ns = 0` with `resolved = false`.

### ROOT Output

- `TH2Poly`: HyCal offset map (display clamped to `[-8, 8] ns` by default).
- All neighbor-pair `h_dt_*` histograms for fit-quality inspection.

---

## 5. Where and How Module Offsets Are Applied

Distinguish two contexts: inside the calibration tool itself, and in general reconstruction/analysis flows.

## 5.1 Inside This Calibration Tool

In `hycal_module_time_calib.cpp`, application happens when building `event_times`:

- Validation OFF: raw peak times.
- Validation ON: `ev.peak_time - mod->time_offset`.

So the same algorithm supports two input definitions:

- Raw-time mode for solving offsets.
- Corrected-time mode for validating offsets.

## 5.2 Generic Reconstruction/Analysis Usage

`HyCalTimeCalib.h` writes offsets directly into `fdec::Module::time_offset`. Any place that already has a module pointer can use:

$$
 t_{calib} = t_{raw} - mod->time_offset
$$

Recommended event-loop convention:

1. Get `mod = hycal.module_by_id(module_id)`.
2. Read `raw_time = ev.peak_time[...]`.
3. Compute `calib_time = raw_time - mod->time_offset`.
4. Use `calib_time` consistently for timing cuts, pair differences, and cluster timing studies.

---

## 6. Applying Offsets in Replay and Other Tools

This section focuses on how calibrated offsets become active in real workflows.

The PRad-II raw-data-replay don't apply this time calibration. The recon-data-replay tool (`prad2ana_replay_recon`) can automatically apply offsets if the correct recon config is provided. Other analysis tools can also use the same `HyCalTimeCalib.h` interface to load and apply offsets.

## 6.1 Automatic Application in Replay (PRad-II Recon)

PRad-II reconstruction goes through `PipelineBuilder`, which reads HyCal time-calibration settings from recon config:

- In `database/reconstruction_config.json`:
  - `hycal.time_calib` (enable/disable)
  - `hycal.time_calib_file` (for example, `hycal_time_offsets/test.json`)
- In `prad2det/src/PipelineBuilder.cpp`, `LoadHyCalTimeCalib(...)` is called during build and writes values into each `Module::time_offset`.

Therefore `analysis/tools/replay_recon.cpp` does not need a separate manual offset-load step. Pass the correct recon config and the offsets are applied automatically.

Minimal config snippet (`database/reconstruction_config.json`):

```json
{
  "hycal": {
    "time_calib": true,
    "time_calib_file": "hycal_time_offsets/runXXXXX.json"
  }
}
```

Example replay recon run:

```bash
prad2ana_replay_recon /path/to/evio_or_dir \
  -c database/reconstruction_config.json \
  -o /path/to/recon_output
```

When `time_calib = true`, downstream replay logic can consistently use:

$$
t_{calib} = t_{raw} - mod->time_offset
$$

If `time_calib = false`, `PipelineBuilder` resets offsets to 0 (equivalent to no calibration).

## 6.2 Behavior in Replay (PRad-I Path)

The PRad-I branch does not use the PRad-II RF/offset pipeline. In current code, its offset table is kept uniform at 0 (no-op).

This means:

- PRad-I does not automatically consume `hycal.time_calib_file` in the same way.
- Compatibility interfaces exist, but no PRad-II-style per-module timing-offset correction is effectively applied.

## 6.3 Manual Application in Other Analysis Tools (Recommended Template)

For standalone tools (including `hycal_module_time_calib.cpp` itself or future timing diagnostics), use this standard pattern:

1. Initialize `HyCalSystem` (map + energy calibration).
2. Call `LoadHyCalTimeCalib(offset_json, hycal, 0.f)`.
3. In the event loop, whenever `mod` and `raw_time` are available, compute `calib_time` uniformly.

Reference code:

```cpp
fdec::HyCalSystem hycal;
hycal.Init(db_dir + "/hycal_map.json");
hycal.LoadCalibration(db_dir + "/" + gRunConfig.energy_calib_file);
prad2::LoadHyCalTimeCalib(db_dir + "/hycal_time_offsets/runXXXXX.json", hycal, 0.f);

const auto *mod = hycal.module_by_id(ev.module_id[j]);
if (!mod) continue;

const float raw_time = ev.peak_time[j][0];
const float calib_time = raw_time - mod->time_offset;
```

In practice, keep variable names explicit (`raw_time` vs `calib_time`) to avoid mixing corrected and uncorrected quantities in cuts or histograms.

## 6.4 Two Offset-Usage Modes in This Calibration Tool

`analysis/tools/hycal_module_time_calib.cpp` already supports both roles:

- Default mode: do not subtract offsets; infer a new table from raw times.
- `-v` + `-V <json>`: subtract existing offsets first, then evaluate whether pair means move closer to 0 and widths become smaller.

This is the mechanism that separates “offset production” from “offset validation”.

---

## 7. Practical Commands and Workflow

Assume executable: `prad2ana_hycal_module_time_calib`.

## 7.1 Produce a New Offset Table (No Validation)

```bash
prad2ana_hycal_module_time_calib /path/to/raw_or_dir \
  -o database/hycal_time_offsets/runXXXXX
```

Outputs:

- `database/hycal_time_offsets/runXXXXX.json`
- `database/hycal_time_offsets/runXXXXX.root`

## 7.2 Validate the Produced Table

```bash
prad2ana_hycal_module_time_calib /path/to/raw_or_dir \
  -o analysis/validation/runXXXXX_val \
  -v -V database/hycal_time_offsets/runXXXXX.json
```

What to check:

- Pair-histogram means should move closer to 0.
- RMS/sigma should become smaller.
- `Valid pair results` count and network connectivity should remain stable.

## 7.3 Enable This Table in Replay Recon

1. Place `runXXXXX.json` under `database/hycal_time_offsets/`.
2. Set in `database/reconstruction_config.json`:
   - `hycal.time_calib = true`
   - `hycal.time_calib_file = "hycal_time_offsets/runXXXXX.json"`
3. Re-run `prad2ana_replay_recon`.

---

## 8. Common Caveats

1. `reference_module` is currently hard-coded to `1495`.
   - If that module is poorly connected in the selected dataset, coverage can drop.

2. Only `is_pwo4()` modules are solved.
   - Non-PWO4 modules are outside this offset solve.

3. A pair must have `entries > 100` to be fitted.
   - Low statistics reduce solvable-module coverage.

4. `mean_error` has a floor of `0.01 ns`.
   - This prevents a single edge from getting unrealistically large weight.

5. Unresolved modules are written as `offset_ns = 0` and `resolved = false`.
   - Downstream quality control should use `resolved/connected` flags.

---

## 9. One-Sentence Summary

HyCal time calibration here is: build clean neighbor-pair timing constraints, convert them into a reference-anchored linear network, solve per-module constant offsets with weighted least squares, and then apply uniformly as `t_calib = t_raw - offset_module` across replay and downstream timing analysis.
