# HyCal Physics Calibration — `physics_calib` Algorithm and Usage Guide

Source code: [`analysis/tools/physics_calib.cpp`](../../analysis/tools/physics_calib.cpp)

This document explains how `physics_calib` calibrates the 1156 HyCal PbWO4
modules with elastic electron-proton (`ep`) events. It covers the input data and
configuration, event selection, 5×5 energy sum, elastic-peak fit, iterative
update, output interpretation, and a practical running procedure.

---

## Contents

1. [Calibration Goal and Workflow](#1-calibration-goal-and-workflow)
2. [Input Data and Configuration](#2-input-data-and-configuration)
3. [Event Reconstruction and Selection](#3-event-reconstruction-and-selection)
4. [Constructing the 5×5 Energy](#4-constructing-the-55-energy)
5. [Per-Module Fit and Calibration Update](#5-per-module-fit-and-calibration-update)
6. [Multithreaded Processing](#6-multithreaded-processing)
7. [Output Files and Interpretation](#7-output-files-and-interpretation)
8. [Building and Running](#8-building-and-running)
9. [Recommended Iteration and Validation Procedure](#9-recommended-iteration-and-validation-procedure)
10. [Implementation Details and Caveats](#10-implementation-details-and-caveats)

---

## 1. Calibration Goal and Workflow

### 1.1 Quantity being calibrated

For module $i$, let its calibration constant in iteration $n$ be
$C_i^{(n)}$. The raw pulse integral is first multiplied by the event-by-event
gain correction and then converted to energy with the current calibration:

$$
A_{i,\mathrm{corr}}=A_{i,\mathrm{integral}}g_i,
\qquad
E_i=\operatorname{energize}_i(A_{i,\mathrm{corr}};C_i^{(n)}).
$$

Here, $g_i$ comes from `hycal.gain_factor` in the raw ROOT file and corrects
time-dependent gain drift. The physics calibration constant $C_i^{(n)}$ is the
quantity updated by this program. These two corrections serve different
purposes and should not be confused.

For clean single-cluster events centered on module $i$, the program builds a
5×5 energy spectrum. If the fitted elastic peak is $E_{i,\mathrm{fit}}$ and the
kinematic expectation is $E_{i,\mathrm{expect}}$, the full correction would be

$$
R_i^{\mathrm{raw}}=
\frac{E_{i,\mathrm{expect}}}{E_{i,\mathrm{fit}}}.
$$

Only 70% of the deviation is applied to reduce overshoot and oscillation:

$$
R_i=1+0.7\left(R_i^{\mathrm{raw}}-1\right),
\qquad R_i\in[0.5,2.0],
$$

$$
C_i^{(n+1)}=C_i^{(n)}R_i.
$$

| Fit result | Interpretation | Update direction |
|---|---|---|
| $E_\mathrm{fit}<E_\mathrm{expect}$ | Reconstructed energy is too low | $R>1$; increase the factor |
| $E_\mathrm{fit}>E_\mathrm{expect}$ | Reconstructed energy is too high | $R<1$; decrease the factor |
| $E_\mathrm{fit}=E_\mathrm{expect}$ | Peak is aligned | $R=1$; no change |

### 1.2 Overall data flow

```text
raw ROOT files (events tree)
        |
        | load run info, detector map, reconstruction config,
        | time cuts/offsets, and current calibration constants
        v
per-file reconstruction in worker threads
        |
        | sum trigger -> PbWO4 peaks -> HyCal clustering
        | -> one clean cluster -> 5x5 energy
        v
per-thread module spectra
        |
        v
single-thread histogram merge
        |
        | elastic ep expected energy + Gaussian fit
        v
calib_factor_iterN.json + calib_result_iterN.{json,root}
        |
        v
iteration N+1 automatically uses calib_factor_iterN.json
```

---

## 2. Input Data and Configuration

### 2.1 Input ROOT files

The program reads the `events` TTree from raw ROOT files. Positional inputs may
be files, directories, or a mixture of both. For a directory input, only
regular files whose names contain `_raw.root` are collected, in sorted order.
A directly specified file is not filtered by name, but it must contain an
`events` tree.

The following branches are enabled:

```text
event_num, trigger_type, trigger_bits
hycal.nch, hycal.module_id, hycal.module_type
hycal.gain_factor, hycal.npeaks
hycal.peak_height, hycal.peak_time, hycal.peak_integral
```

The run number is extracted from the first input filename. It determines the
global run configuration and output directory. All inputs should therefore
belong to the same run; the current code does not reject mixed-run inputs.

### 2.2 Database directory

The database directory is resolved in this order:

1. `PRAD2_DATABASE_DIR` environment variable;
2. `../share/prad2evviewer/database`;
3. the build-time `DATABASE_DIR` value.

For a nonstandard installation or a build-tree run, set it explicitly:

```bash
export PRAD2_DATABASE_DIR=/path/to/prad2evviewer/database
```

The relevant database inputs include:

- `hycal_map.json`: module IDs, positions, sizes, types, and geometry;
- `reconstruction_config.json`: clustering parameters, including
  `seed_time_window`;
- `runinfo/general.json`: beam energy, HyCal $z$ position, dead modules, and
  other run-dependent settings;
- run-dependent HyCal time cuts and time offsets loaded by `PipelineBuilder`;
- the default iteration-1 calibration seed,
  `calibration/calibration_factor_3p5_June7.json`.

### 2.3 Calibration file used by each iteration

| Iteration | Input calibration | Output calibration |
|---|---|---|
| `-i 1` | `-c <seed.json>`, or the database default | `calib_factor_iter1.json` |
| `-i N`, $N>1$ | `calib_factor_iter(N-1).json` in the same run output directory | `calib_factor_iterN.json` |

The `-c` option affects iteration 1 only. At higher iterations it is ignored
with a warning, and the preceding iteration is used.

---

## 3. Event Reconstruction and Selection

Each worker calls `ProcessRawFiles` for one input file. `PipelineBuilder`
constructs the HyCal system and loads the calibration, clustering settings,
time cuts, and time offsets appropriate for that file.

### 3.1 Basic event selection

| Cut | Code condition | Purpose |
|---|---|---|
| Sum trigger | `(trigger_bits & TBIT_sum) != 0` | Select total-energy physics triggers |
| Channel multiplicity | `nch <= 100` | Reject high-occupancy or unclean events |

Only PbWO4 modules are used. Pb-glass modules and channels not found in the map
are skipped.

### 3.2 Constructing module hits from waveform peaks

Each peak time is corrected with the module offset:

$$
t_\mathrm{corr}=t_\mathrm{peak}-t_\mathrm{offset,module}.
$$

Only peaks strictly inside the run-dependent module time window are retained:

$$
t_\mathrm{lo}<t_\mathrm{corr}<t_\mathrm{hi}.
$$

Processing then depends on `seed_time_window`.

#### Multi-pulse mode: `seed_time_window > 0`

- Every in-window peak is sent to the clusterer and stored in `valid_peaks`.
- The clusterer applies seed-anchored timing consistency.
- During the 5×5 sum, the code again requires
  $|t_\mathrm{peak}-t_\mathrm{cluster}|\leq\texttt{seed_time_window}$.

#### Legacy mode: `seed_time_window <= 0`

- Only the largest-integral in-window peak is retained for each module.
- Each module contributes at most one hit.
- No additional seed-time coincidence cut is applied during the 5×5 sum.

The variable named `bestHeight` is compared against `peak_integral`; legacy
mode therefore selects the largest integral, not the largest `peak_height`.

### 3.3 Cluster-level selection

After `FormClusters()` and `ReconstructHits()`, the event must satisfy:

1. exactly one reconstructed cluster: `hits.size() == 1`;
2. at least four blocks: `hits[0].nblocks > 3`;
3. a PbWO4 seed/center module;
4. no `kDeadModule` flag on the cluster;
5. a reconstructed position near the seed-crystal center.

The normalized local coordinates are

$$
x_d=\frac{x_\mathrm{hit}-x_\mathrm{module}}{\mathrm{size}_x},
\qquad
y_d=\frac{y_\mathrm{hit}-y_\mathrm{module}}{\mathrm{size}_y},
$$

with the requirement

$$
|x_d|<0.3,\qquad |y_d|<0.3.
$$

This selects the central $60\%\times60\%$ of the seed crystal, reducing the
position dependence of shower leakage.

### 3.4 Additional transition-module cut

For a cluster carrying `kTransition`, only the side facing the HyCal center is
accepted:

- right side ($x>300$ mm): require $x_d<0$;
- left side ($x<-300$ mm): require $x_d>0$;
- top side ($y>300$ mm): require $y_d<0$;
- bottom side ($y<-300$ mm): require $y_d>0$.

At a corner, both applicable conditions must pass. This suppresses showers on
the outer side of transition regions, where leakage and material-boundary
effects are larger.

---

## 4. Constructing the 5×5 Energy

### 4.1 Geometric window

The code uses a fixed PbWO4 pitch:

$$
p=20.75\ \mathrm{mm}.
$$

Around the cluster center module, it retains `valid_peaks` satisfying

$$
|x_j-x_c|\leq2.5p,\qquad |y_j-y_c|\leq2.5p.
$$

The half-width is 51.875 mm, corresponding to offsets $-2,-1,0,1,2$
in both directions on a regular PbWO4 grid. All geometrically and temporally
accepted peak energies are summed:

$$
E_{5\times5}=\sum_{j\in5\times5}E_j.
$$

This sum is built from `valid_peaks`, rather than directly using
`hits[0].energy`, and is filled into the spectrum of the cluster center module.

### 4.2 Energy layers

For each peak, the position difference is divided by the pitch and rounded:

$$
n_x=\operatorname{round}(\Delta x/p),\qquad
n_y=\operatorname{round}(\Delta y/p),
$$

and the layer is the Chebyshev distance:

$$
L=\max(|n_x|,|n_y|).
$$

| Layer | Region | Modules in a complete window |
|---|---|---:|
| 0 | Center crystal | 1 |
| 1 | Central 3×3 excluding the center | 8 |
| 2 | Outer ring of the 5×5 | 16 |

The corresponding sums are stored as `center_energy`,
`second_layer_energy`, and `third_layer_energy`.

### 4.3 Center-energy fraction cut

The event must satisfy

$$
f_\mathrm{center}=
\frac{E_\mathrm{center}}{E_\mathrm{cluster}}>0.6.
$$

The denominator is the clusterer result `hits[0].energy`, not
$E_{5\times5}$. The diagnostic layer fractions use the same denominator:

$$
f_2=\frac{E_\mathrm{layer2}}{E_\mathrm{cluster}},
\qquad
f_3=\frac{E_\mathrm{layer3}}{E_\mathrm{cluster}}.
$$

The current ROOT histogram titles describe the latter two as
`E_layer/E_5x5`, but the implementation divides by `hits[0].energy`. Interpret
the histograms according to the implementation.

---

## 5. Per-Module Fit and Calibration Update

After merging all worker histograms, the main thread processes W1–W1156.

### 5.1 Minimum statistics and dead modules

- A spectrum with fewer than 100 entries is skipped completely: it is not fit
  and does not appear in `calib_result_iterN.json`.
- A run-config dead module is not updated.
- A dead module appears in the result JSON only if its spectrum first reaches
  100 entries; its numeric fields are zero and `fit_good=false`.
- `kDeadNeighbor` is recorded as `is_deadNeighbor`, but does not currently
  prevent fitting or updating.

The output calibration table is written from the complete calibration already
loaded into `HyCalSystem`. Skipped modules therefore normally retain their
input factors. The result JSON is a diagnostic subset, not a complete
calibration table.

### 5.2 Expected elastic `ep` energy

The representative angle is computed from the module center, not the
event-by-event reconstructed position:

$$
\theta_i=\tan^{-1}\left(
\frac{\sqrt{x_i^2+y_i^2}}{z_\mathrm{HyCal}}
\right).
$$

The module coordinates come from `hycal_map.json`; $z_\mathrm{HyCal}$ and
$E_\mathrm{beam}$ come from the run configuration. The elastic energy is

$$
E'=\frac{E_\mathrm{beam}M_p}
{M_p+E_\mathrm{beam}(1-\cos\theta)},
$$

followed by the approximate material energy-loss correction implemented in
`PhysicsTools::EnergyLoss`:

$$
E_\mathrm{expect}=E'-\Delta E_\mathrm{material}.
$$

Beam energy, HyCal position, module geometry, and the energy-loss model can
therefore all propagate into the final calibration.

### 5.3 Gaussian peak fit

`PhysicsTools::fitGaus(h, expected_peak)` performs the following steps:

1. require at least 100 entries;
2. find the highest local maximum in
   $[0.8E_\mathrm{expect},1.2E_\mathrm{expect}]$;
3. fall back to the global maximum if no local peak is found there;
4. find the boundaries where the spectrum drops below 40% of the peak height;
5. initialize and run a Gaussian fit over that interval;
6. return `(peak, sigma, chi2/ndf)`, or `(0,0,0)` for an invalid fit.

The program defines `fit_good` as

$$
E_\mathrm{fit}>0,\qquad \sigma>0,
$$

$$
\sigma<2\times0.03E_\mathrm{fit}
\sqrt{\frac{1000}{E_\mathrm{fit}}},
\qquad \chi^2/\mathrm{NDF}<10.
$$

The $0.03/\sqrt{E/1000}$ term is the approximate relative-resolution scale
used by the code; the fit cut permits twice that width.

### 5.4 Update rule

For a positive fitted peak:

```text
raw_ratio  = expected_peak / peak
ratio      = 1 + 0.7 * (raw_ratio - 1)
ratio      = clamp(ratio, 0.5, 2.0)
new_factor = old_factor * ratio
```

The module calibration base energy is also set to `expected_peak`.

Importantly, `fit_good` is currently a diagnostic flag, **not an update gate**.
A positive peak is used even when its sigma or chi-square fails the quality
criteria. If the fit returns `peak <= 0`, the code substitutes the expected
peak, giving ratio 1 and leaving the factor unchanged.

---

## 6. Multithreaded Processing

The default is four threads, limited to

$$
1\leq N_\mathrm{threads}\leq N_\mathrm{files}.
$$

Files are processed in rounds. Each thread handles one complete ROOT file per
round and reuses its `HistResult` in later rounds, so its histograms accumulate
across files. After all rounds, the main thread merges the 1156 module spectra,
global diagnostic histograms, and event counts.

Fitting, factor updates, and output writing are single-threaded. Increasing
`-j` accelerates raw-file reconstruction only, and the effective thread count
cannot exceed the number of input files.

---

## 7. Output Files and Interpretation

All outputs are written under

```text
<output_dir>/Physics_calib/run<run_number>/
```

For example, `-o ./calib_output` with run 23453 produces
`./calib_output/Physics_calib/run23453/`.

### 7.1 `calib_factor_iterN.json`

This is the complete calibration table written by
`HyCalSystem::PrintCalibConstants`. It is the actual input to iteration $N+1$
and the file to preserve or deploy as the final calibration.

Each array entry currently has this form:

```json
{
  "name": "W735",
  "factor": 0.1532,
  "base_energy": 3478.6,
  "nl1": 0.0,
  "nl2": 0.0
}
```

| Field | Meaning | Interpretation |
|---|---|---|
| `name` | Module name such as `W735` or `G1` | Only W1–W1156 are updated; other modules remain in the complete table |
| `factor` | Integral-to-energy calibration factor | The value used in the next iteration |
| `base_energy` | Calibration reference energy | Set to the elastic expected peak for a processed W module |
| `nl1` | First nonlinearity parameter | Not refit; inherited from the input table |
| `nl2` | Second nonlinearity parameter | Not refit; inherited from the input table |

For a normally updated module,

$$
\texttt{factor}_{N}
=\texttt{factor}_{N-1}\times\texttt{ratio}_{N}.
$$

Low-statistics, dead, and completely failed-fit modules normally keep their old
factor. A failed fit executes a no-change update with ratio 1; a dead module is
skipped before updating.

### 7.2 `calib_result_iterN.json`

This diagnostic array contains one object per processed module:

| Field | Meaning |
|---|---|
| `module_id` | Numeric ID; W1 corresponds to 1001 |
| `old_factor` | Input factor for this iteration |
| `new_factor` | Updated factor |
| `ratio` | Damped and clamped correction factor |
| `peak` | Fitted peak; replaced by the expectation after a complete fit failure |
| `expected_peak` | Elastic expectation at the module-center angle |
| `sigma` | Gaussian sigma |
| `chi2/ndf` | Gaussian-fit $\chi^2/\mathrm{NDF}$ |
| `fit_good` | Whether the code's fit-quality criteria pass |
| `is_dead` | Run-config dead-module flag |
| `is_deadNeighbor` | Dead-neighbor flag |

Only modules with at least 100 spectrum entries appear. Do not use this file
alone to determine whether all 1156 modules are present in the final table.

Example:

```json
{
  "module_id": 1735,
  "old_factor": 0.1500,
  "new_factor": 0.15525,
  "ratio": 1.035,
  "peak": 3300.0,
  "expected_peak": 3465.0,
  "sigma": 72.0,
  "chi2/ndf": 1.4,
  "fit_good": true,
  "is_dead": false,
  "is_deadNeighbor": false
}
```

The fitted peak is 5% below expectation. The undamped correction is
$3465/3300=1.05$; after damping it is $1+0.7\times0.05=1.035$, so the factor
increases from 0.1500 to 0.15525.

| Result pattern | Likely meaning | Recommended action |
|---|---|---|
| `ratio` near 1 and `fit_good=true` | Module is approaching convergence | Confirm stability across iterations |
| `ratio > 1` | Fitted peak is low | New factor increases |
| `ratio < 1` | Fitted peak is high | New factor decreases |
| `ratio` exactly 0.5 or 2.0 | Raw correction hit the clamp | Inspect the spectrum, seed, and data quality |
| `fit_good=false` but `peak > 0` | Peak exists but width or chi-square fails | The code still updates; inspect manually |
| `peak == expected_peak` and `sigma == 0` | Typical complete-fit fallback | Ratio 1 does not mean a perfect fit |
| `is_dead=true` with zero numeric fields | Dead module was skipped | Do not treat zero as a new valid factor |
| `is_deadNeighbor=true` | Leakage into a dead region is possible | Inspect peak shape and iteration stability |
| Module absent from result JSON | Usually fewer than 100 entries | Check its ROOT spectrum; factor normally remains unchanged |

### 7.3 `calib_result_iterN.root`

| ROOT object | Contents |
|---|---|
| `h1_E_mod_<id>_merged` | 5×5 seed spectrum for every W module |
| `h2_energy_theta_merged` | Selected $E_{5\times5}$ versus reconstructed $\theta$ |
| `hit_pos_merged` | Selected reconstructed hit positions |
| `h_E_1cl_merged` | Inclusive selected single-cluster 5×5 energy |
| `h_center_energy_fraction` | $E_\mathrm{center}/E_\mathrm{cluster}$ |
| `h_center_energy` | Center-crystal energy |
| `h_2nd_energy_fraction` | Actually $E_\mathrm{layer2}/E_\mathrm{cluster}$ in the code |
| `h_3rd_energy_fraction` | Actually $E_\mathrm{layer3}/E_\mathrm{cluster}$ in the code |
| `h_fit_peak_energy` | Fitted or fallback peaks across modules |
| `h_fit_peak_ratio` | Damped correction-factor distribution |
| `h_fit_peak_chi2ndf` | Module-fit $\chi^2/\mathrm{NDF}$ distribution |
| `h_fit_peak_sigma` | Module-fit sigma distribution |

#### Per-module input spectra

`h1_E_mod_<id>_merged` is the most direct evidence for whether a module update
is trustworthy. Check that the elastic peak is clear and near its expectation,
that the fit did not select a background peak, and that the spectrum has enough
entries. Double peaks, truncation, strong tails, and unusual edge or
dead-neighbor shapes require manual review. The JSON `peak`, `sigma`, and
`chi2/ndf` all come from these spectra.

#### Global event-selection diagnostics

- `h2_energy_theta_merged`: elastic events should form a continuous kinematic
  band. A global offset suggests energy-scale, beam-energy, or geometry issues;
  a local distortion is more likely module-specific.
- `hit_pos_merged`: reveals coverage, dead-region holes, and strong spatial
  selection nonuniformities.
- `h_E_1cl_merged`: shows the inclusive selected spectrum but cannot replace
  individual module spectra.
- `h_center_energy_fraction`: selected entries begin above 0.6. A large pileup
  at the threshold indicates sensitivity to this cut.
- The center and layer histograms show shower-energy sharing and help diagnose
  leakage near edges and dead regions.

#### Cross-module fit summaries

- `h_fit_peak_energy` need not be narrow because the expected energy changes
  with module angle.
- `h_fit_peak_ratio` is the main convergence summary. It should narrow toward 1
  over iterations. Pileups at 0.5 or 2.0 indicate widespread clamping.
- A long `h_fit_peak_chi2ndf` tail or many entries above 10 means a Gaussian
  poorly describes many spectra.
- Large sigma may indicate background, double peaks, leakage, or peak
  misidentification; unusually small sigma may indicate a fit over too few bins.

### 7.4 Relationship among the three outputs

```text
calib_result_iterN.root
  per-module spectrum h1_E_mod_<id>_merged
             |
             | Gaussian fit + expected ep energy
             v
calib_result_iterN.json
  old_factor, peak, expected_peak, ratio, quality flags
             |
             | new_factor = old_factor * ratio
             v
calib_factor_iterN.json
  complete calibration table used by iteration N+1
```

To investigate an abnormal module:

1. find its `module_id`, ratio, and flags in the result JSON;
2. inspect the corresponding ROOT module spectrum;
3. find `W<module_id-1000>` in the factor JSON and verify the final factor;
4. compare adjacent iterations to distinguish convergence, oscillation, low
   statistics, and wrong-peak selection.

---

## 8. Building and Running

### 8.1 Build

The tool requires ROOT and the analysis targets:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release -DBUILD_ANALYSIS=ON
cmake --build build --target physics_calib -j4
```

The build-tree executable is `build/bin/prad2ana_physics_calib`. For an
installed build, use `prad2ana_physics_calib` from `PATH`.

### 8.2 Command line

```text
prad2ana_physics_calib <input_raw.root|dir> [more files/dirs...]
    [-i iteration] [-o output_dir] [-c seed_calib.json] [-j num_threads]
```

| Argument | Default | Description |
|---|---|---|
| Positional input | Required | Raw ROOT file or directory containing `*_raw.root`; multiple allowed |
| `-i <N>` | `1` | Iteration number, at least 1 |
| `-o <dir>` | `.` | Output base; `Physics_calib/runN` is appended |
| `-c <json>` | Database seed | Initial calibration for iteration 1 only |
| `-j <N>` | `4` | Reconstruction workers, limited to the number of files |

The no-input usage string in the source still says `calib_5by5`; this is a
stale message. The actual executable is `prad2ana_physics_calib`.

### 8.3 Examples

```bash
export PRAD2_DATABASE_DIR=/path/to/prad2evviewer/database

./build/bin/prad2ana_physics_calib \
    /data/run23453/raw_root/ \
    -o ./calib_output \
    -i 1 \
    -j 8
```

With an explicit iteration-1 seed:

```bash
./build/bin/prad2ana_physics_calib \
    /data/run23453/prad_023453_0000_raw.root \
    /data/run23453/prad_023453_0001_raw.root \
    -o ./calib_output \
    -i 1 \
    -c ./seed_calibration.json \
    -j 2
```

Iteration 2 automatically reads
`./calib_output/Physics_calib/run23453/calib_factor_iter1.json`:

```bash
./build/bin/prad2ana_physics_calib \
    /data/run23453/raw_root/ \
    -o ./calib_output \
    -i 2 \
    -j 8
```

### 8.4 Iteration loop

```bash
#!/usr/bin/env bash
set -euo pipefail

EXE=./build/bin/prad2ana_physics_calib
RAW_DIR=/data/run23453/raw_root
OUT_DIR=./calib_output
THREADS=8

for ITER in 1 2 3 4 5; do
    "$EXE" "$RAW_DIR" -o "$OUT_DIR" -i "$ITER" -j "$THREADS"
done
```

Use the same input sample, output base, and run at every iteration. Otherwise
the preceding factor file may not be found, or changing statistics may obscure
the convergence behavior.

---

## 9. Recommended Iteration and Validation Procedure

### 9.1 Before running

1. Confirm that all raw ROOT files belong to one run and contain the required
   `events` branches.
2. Confirm that `PRAD2_DATABASE_DIR` matches the data-taking period.
3. Check `Ebeam`, `hycal_z`, and dead modules for the run in
   `runinfo/general.json`.
4. Confirm that gain corrections are present in `hycal.gain_factor`.
5. Validate the time cuts, time offsets, and reconstruction configuration.
6. Keep one consistent output path across iterations.

### 9.2 After each iteration

List suspicious modules:

```bash
jq '.[] | select((.fit_good == false) or (.ratio < 0.9) or (.ratio > 1.1))' \
    calib_output/Physics_calib/run23453/calib_result_iter1.json
```

Count good and bad fits:

```bash
jq '{total: length, good: map(select(.fit_good)) | length,
     bad: map(select(.fit_good == false)) | length}' \
    calib_output/Physics_calib/run23453/calib_result_iter1.json
```

In ROOT, inspect the ratio, chi-square, and sigma summaries; the energy-angle
band and hit coverage; and individual spectra for failed fits, clamped ratios,
dead neighbors, and transition modules.

### 9.3 Convergence criteria

Do not judge convergence only by successful program exit. Useful criteria are:

- most module ratios are close to 1;
- relative factor changes between consecutive iterations are small;
- the elastic peak no longer has a significant position-dependent bias;
- failed modules are stable and have been reviewed individually;
- few or no modules remain at the 0.5 or 2.0 clamp.

Because each iteration applies only 70% of the correction, multiple iterations
are normally required. The necessary count depends on the seed, statistics,
and background and should not be fixed mechanically.

---

## 10. Implementation Details and Caveats

1. **PbWO4 only:** the program calibrates W1–W1156; Pb-glass is outside this
   workflow.

2. **Raw ROOT input only:** the tree must be `events`, not `recon`; EVIO and
   reconstructed ROOT files are not valid inputs.

3. **The first file determines global run state:** each worker infers its own
   run for the pipeline, but mixed runs are not rejected.

4. **No event-count option:** an unused `max_events` variable exists, but there
   is no `-n` option. Every entry in every input file is processed.

5. **A failed input file does not abort the full job:** the worker reports
   `FAILED`, other files continue, and outputs may still be written. Review the
   log for failed files.

6. **Calibration loading is not explicitly validated by this tool:** verify
   that the seed or preceding JSON exists, parses correctly, and contains the
   expected modules.

7. **`fit_good` does not prevent an update:** any positive fitted peak produces
   a correction even if its width or chi-square is unacceptable.

8. **Low-statistics modules are absent from the result JSON:** modules below
   100 entries retain their old factor without an explicit low-stat status.

9. **Dead-neighbor is currently diagnostic only:** these modules are updated
   normally despite possible shower leakage into a dead crystal.

10. **The 5×5 uses a fixed 20.75 mm pitch:** it is a coordinate window rather
    than an explicit neighbor graph; edge and transition geometry is handled
    only by the geometric and inner-side cuts.

11. **Layer-fraction titles differ from the implementation:** the code divides
    layer 2 and layer 3 by cluster energy, not 5×5 energy.

12. **The merged hit-position title has an inconsistent unit:** the filled
    HyCal coordinates are used elsewhere as mm, although the merged histogram
    title says `X (cm);Y (cm)`.

In summary, this calibration selects clean single-cluster elastic-`ep` showers,
assigns each 5×5 spectrum to its seed PbWO4 module, aligns the fitted peak with
the kinematic expectation, and iteratively updates the module factor with a
damped multiplicative correction. A usable final calibration requires checking
peak shape, fit quality, dead and edge regions, statistical coverage, and
convergence across iterations.

---

## 11. Physics Calibration Viewer User Guide

`physics_calib_viewer.py` is a desktop viewer for inspecting the per-module
spectra and fit results produced by `physics_calib.cpp`. It is intended for
diagnosis and controlled manual refitting after an automated calibration pass.
The viewer can update the result and factor JSON files, so make a backup before
editing production calibration output.

### 11.1 Starting the viewer

From the source checkout, run the viewer with the `Physics_calib` directory as
the optional positional argument:

```bash
python3 scripts/physics_calib_viewer.py build/Physics_calib
```

The histogram source can be selected at startup or from the top toolbar:

```bash
python3 scripts/physics_calib_viewer.py build/Physics_calib --hist-mode 5by5
python3 scripts/physics_calib_viewer.py build/Physics_calib --hist-mode island
```

`5by5` reads the module spectra from the `modules_5by5` ROOT directory;
`island` reads them from `modules_island`. The run and iteration selectors are
populated from directories containing `calib_result_iterN.json`. The selected
iteration also loads its matching `calib_factor_iterN.json` and
`calib_result_iterN.root`.

The ROOT file is loaded in a background worker. While it is loading, the map
and JSON result information may already be visible, but a module histogram or
global ROOT plot may remain empty until loading finishes.

### 11.2 Main window layout

The window has four functional regions:

1. **Top toolbar**: open a `Physics_calib` directory, select a run and
   iteration, and choose the `5by5` or `island` histogram source.
2. **Map mode row and range control**: select the quantity painted on the
   HyCal map and adjust its color range, palette, or log mapping.
3. **HyCal map**: inspect spatial patterns, select modules, preview outer
   layers, and perform batch factor operations.
4. **Module detail and global diagnostics**: inspect one module's energy
   spectrum, refit it, and view run-wide ROOT diagnostic histograms.

The top controls are used as follows:

- **Open Physics_calib**: open a directory chooser and load another
  `Physics_calib` output directory. The directory must contain run subdirectories
  such as `run25320`.
- **Run**: choose which run directory to inspect. Changing it repopulates the
  iteration selector.
- **Iteration**: choose the matching `calib_result_iterN.json`, factor JSON, and
  ROOT file. Changing it clears temporary map selections and histogram display
  state, then starts a new ROOT load.
- **Histogram**: choose `5by5` or `island`. Changing it reloads the module
  spectra from the corresponding ROOT directory and resets the current module's
  temporary rebin state.
- **Map mode buttons**: select the quantity painted on the left map. Switching
  to `Has data` or `fit_good` fixes the range to `[0, 1]`; the other map modes
  recompute an automatic continuous range.
- **Map range controls**: use the minimum/maximum edits to set a manual range,
  `Auto` to recompute it, and `Log` to toggle logarithmic color mapping. These
  controls affect the map display only and do not alter calibration files.

### 11.3 Map modes and what to check

The map buttons show one value per module when that value exists in the result
JSON. Modules absent from the result JSON are displayed as no data. This is
important because the factor JSON contains a complete calibration table, while
the result JSON only contains modules with enough processed statistics.

#### `module flag`

This discrete map encodes the module status from the result JSON with a fixed
range of `[0, 2]`:

| Value | Status |
|---:|---|
| `0` | Normal module |
| `1` | Dead-neighbor module |
| `2` | Dead module |

Dead has priority over dead-neighbor if both flags are present. Use this map
to locate geometric regions affected by dead hardware and to compare them with
large sigma, poor chi-square, or large relative peak error. The flag is a
diagnostic classification; it does not by itself change the factor.

#### `Has data`

This is a binary map: `1` means the module has an entry in
`calib_result_iterN.json`, and `0` means it does not. The map range is fixed to
`[0, 1]`.

Use it to check coverage. A missing module normally means that its spectrum did
not reach the producer's minimum statistics threshold. It does not necessarily
mean that the module is dead or absent from the detector. Click a no-data
module to confirm that no result entry exists, then inspect its ROOT spectrum
if available.

#### `fit_good`

This is also a binary map with a fixed range of `[0, 1]`. A value of `1` means
the result passed the producer/viewer fit-quality flag; `0` means it did not.

Look for clusters of bad fits near detector edges, transition regions, or
dead-neighbor regions. A bad fit is a review signal, not proof that the factor
must be restored. Inspect the module histogram, fitted width, chi-square, and
peak location before changing it.

#### `chi2/ndf`

This map displays the Gaussian fit chi-square per degree of freedom. Its range
is automatically fitted to the values currently present.

Large isolated values usually indicate background, a second peak, a poor fit
window, or a module with an unusual shower shape. A broad spatial region of
large values can indicate a selection or calibration problem rather than one
bad channel.

#### `sigma`

This map displays the fitted Gaussian width in MeV. Use it to find unusually
wide peaks, which can indicate background, leakage, double peaks, or an
incorrect peak selection. Very small widths should also be checked because
they can result from too few useful bins or an overly narrow fit window.

#### `delta E`

This map displays the signed peak displacement:

$$
\Delta E = peak - expected\_peak.
$$

Positive values mean the fitted energy is above the kinematic expectation;
negative values mean it is below. Spatially coherent positive or negative
regions are useful for spotting geometry or calibration trends. This map uses
an automatic continuous range, so zero is not necessarily at the center of the
color palette.

#### `|delta E/expected|`

This map displays the absolute relative peak error:

$$
\left|\frac{\Delta E}{expected\_peak}\right|
= \frac{|peak - expected\_peak|}{expected\_peak}.
$$

It is useful for ranking modules by deviation without regard to the direction
of the error. Low values indicate alignment with the expected peak; high values
should be checked against `fit_good`, `sigma`, and the module spectrum.

#### `ratio`

This map displays the damped correction ratio stored in the result JSON. Values
near `1` indicate little correction. Values above `1` increase the factor;
values below `1` decrease it. Values at `0.5` or `2.0` indicate that the
producer clamp was reached and deserve direct inspection.

### 11.4 Map controls and markers

The range control below the map modes provides the following interactions:

- **Minimum and maximum range fields**: manually set the color scale bounds.
- **Auto**: fit the scale to the current map values. The viewer also applies
  automatic ranges when switching continuous map modes.
- **Log**: use logarithmic color mapping when the selected quantity is strictly
  positive. Do not use it to interpret signed `delta E` values.
- **Palette/color bar**: use the shared map palette controls to change the
  visual contrast.
- **Map zoom/pan**: use the map's existing zoom/pan interaction to inspect a
  detector region.

The map uses three independent marker styles:

- **Red circle**: a module whose factor or manual fit has been applied during
  the current viewer session.
- **Blue circle**: a module currently selected in multi-select mode.
- **Yellow dashed circle**: a temporary outer-layer preview.

Markers are visual state only; the yellow preview does not modify any JSON
file. Click `Preview outer` a second time to clear only the yellow preview.
Red and blue markers remain unchanged.

### 11.5 Factor and outer-layer buttons

The controls below the map perform the following operations:

- **Multi-select**: enable or disable module toggling by map clicks. Selected
  modules are shown with blue circles.
- **Factor**: enter a positive factor for the selected modules.
- **Apply selected**: write that factor to all selected module entries in
  `calib_factor_iterN.json`. After a successful write, the selected modules are
  marked red and multi-select mode exits.
- **Restore selected old_factor**: restore the producer `old_factor` from the
  result JSON for every selected module. It updates the factor JSON only, then
  marks the modules red and exits multi-select mode.
- **Outer layers**: choose how many outer layers to include. The spinbox has
  visible up/down buttons and changes by one layer per step.
- **Square / Circle**: choose the definition of an outer layer.
  `Square` uses the row/column distance from the PbWO4 array boundary;
  `Circle` uses radial distance from the geometric center of the PbWO4 array.
- **Restore outer W**: write the selected outer W modules' producer
  `old_factor` values to `calib_factor_iterN.json`.
- **Rebin N** next to the outer controls: choose the number of adjacent bins
  to combine when applying outer-layer rebinning.
- **Preview outer**: show the currently selected outer modules with yellow
  dashed circles. Press it again to clear only that preview.
- **Rebin outer**: apply the outer rebin setting to the selected outer module
  histograms for display. It does not modify the ROOT file.

The outer-layer selection is recalculated from the current shape and layer
count when `Preview outer`, `Restore outer W`, or `Rebin outer` is pressed.

### 11.6 Global diagnostics buttons

The `Global diagnostics` panel contains a row of buttons. Each button selects
one ROOT histogram loaded from `calib_result_iterN.root`:

- **Energy vs theta**: two-dimensional selected cluster energy versus angle.
  Look for the expected kinematic band and discontinuities.
- **Hit position**: two-dimensional reconstructed hit-position coverage.
  Look for holes, edge effects, and dead-region patterns.
- **One cluster energy**: inclusive selected single-cluster energy spectrum.
  Use it as a run-level overview, not as a replacement for module spectra.
- **Peak ratio**: distribution of module correction ratios. A wide distribution
  or pileup at the clamp limits suggests unstable fits or broad calibration
  variation.
- **Fit chi2/ndf**: distribution of fit quality values. A long high-value tail
  indicates spectra that require review.
- **Fit sigma**: distribution of fitted widths. Outliers should be compared
  with their individual module spectra.

The one-dimensional plots automatically adjust their x-axis to the non-zero
data range. The two-dimensional plots use a white background and logarithmic
z scaling so sparse and dense regions can be seen together.

For all global plots, drag a rectangle with the **left mouse button** to zoom
into that region. Drag with the **right mouse button** to zoom out around the
selected region. Scrolling is not required. Click **Reset zoom** to restore the
automatic limits for the current diagnostic plot.

### 11.7 Module histogram controls

Click a module on the map to populate the module detail panel. The information
block shows the module ID, producer peak, expected peak, sigma, chi-square,
ratio, old factor, new factor, and current factor.

Below the histogram:

- **Histogram rebin**: combine adjacent energy bins for the current module.
  The setting is stored separately for each module during the session.
  Rebinning affects display and refitting; it does not rewrite the ROOT file.
- **Restore old factor**: write only the current module's producer
  `old_factor` to `calib_factor_iterN.json`.
- **Range xmin/xmax**: optionally enter explicit fit boundaries in MeV.
- **Run fit**: fit the current histogram using the selected range, or use the
  automatic peak-window logic when no range is entered.
- **Apply fit + save**: write the accepted fit to both JSON outputs and keep
  the manual fit line visible on the plot.

### 11.8 Manual refit procedure

Use the following procedure for a module that needs review:

1. Select the correct run and iteration at the top of the window.
2. Choose `5by5` or `island`, depending on which spectrum should be reviewed.
3. Select a map mode such as `fit_good`, `chi2/ndf`, `delta E`, or
   `|delta E/expected|` to locate suspicious modules.
4. Click the module. Confirm the module information and wait for its histogram
   to load.
5. Set `Histogram rebin` if the spectrum needs a coarser binning. The fit uses
   the rebinned spectrum.
6. Drag across the desired peak region on the histogram. The viewer immediately
   runs a fit using that range. Alternatively, enter `xmin` and `xmax` and
   press `Run fit`.
7. Inspect the manual fit line and the displayed peak, sigma, chi-square,
   damped ratio, and new factor. Compare the manual result with the producer
   fit line and the expected peak.
8. Press `Apply fit + save` only after the fit window and peak are credible.
9. Confirm that the map updates immediately. In particular, `ratio`, `delta E`,
   `|delta E/expected|`, `sigma`, and `fit_good` should reflect the new result.

The ratio and factor use the same update rule as the producer:

$$
ratio = \operatorname{clamp}\left(
1+0.7\left(\frac{expected\_peak}{peak}-1\right), 0.5, 2.0\right),
$$

$$
new\_factor = old\_factor \times ratio.
$$

The `old_factor` used here is taken from the result row and is deliberately
not replaced by the current factor. The current fit updates these fields in
the original `calib_result_iterN.json`:

```text
peak, sigma, chi2/ndf, ratio, new_factor, fit_good
```

The original `old_factor` remains unchanged. The calculated factor is written
to the matching module in `calib_factor_iterN.json`. No
`calib_result_iterN.manual.json` sidecar is created.

If writing either JSON file fails, the viewer reports the failure and does not
claim that the fit was saved. Keep the original files backed up when doing
manual calibration work.

### 11.9 Recommended viewer review order

For a new iteration, a practical review sequence is:

1. Check `Has data` for coverage.
2. Check `fit_good` for broad regions of failed fits.
3. Check `ratio` and `|delta E/expected|` for modules requiring correction.
4. Check `chi2/ndf` and `sigma` for shape or fit-quality problems.
5. Use the global `Energy vs theta` and `Hit position` plots for run-wide
   selection or geometry problems.
6. Inspect suspicious modules individually and refit only when the histogram
   supports a different peak or fit window.
7. Save the manual change, verify the map refresh, and record which modules
   were changed before proceeding to the next iteration.
