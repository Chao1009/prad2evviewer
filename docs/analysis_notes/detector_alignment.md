# Detector Position Alignment — `det_align` Guide and Algorithm

Source code: [`analysis/tools/det_align.cpp`](../../analysis/tools/det_align.cpp)

---

## Table of Contents

1. [Physics Background](#1-physics-background)
2. [Detector Layout and Coordinate System](#2-detector-layout-and-coordinate-system)
3. [Algorithm Overview](#3-algorithm-overview)
4. [Step-by-Step Algorithm](#4-step-by-step-algorithm)
   - [4.1 Event Selection](#41-event-selection)
   - [4.2 Coordinate Transformation](#42-coordinate-transformation)
   - [4.3 Filling Alignment Histograms](#43-filling-alignment-histograms)
   - [4.4 Peak Extraction](#44-peak-extraction)
   - [4.5 Updating Detector Position Parameters](#45-updating-detector-position-parameters)
5. [Iterative Convergence](#5-iterative-convergence)
6. [Output Files](#6-output-files)
7. [Running the Program](#7-running-the-program)
8. [Command-Line Options](#8-command-line-options)
9. [Typical Workflow](#9-typical-workflow)
10. [Known Limitations and Caveats](#10-known-limitations-and-caveats)

---

## 1. Physics Background

Detector alignment uses two types of scattering processes:

**Møller scattering** ($e^- e^- \to e^- e^-$): beam electrons scatter off target electrons, producing two outgoing electrons with similar energies. In the center-of-mass frame the two electrons are back-to-back; in the lab frame the sum of their polar angles is uniquely determined by the beam energy. Møller events provide:
- The **Møller center** (intersection of the two hit-connecting lines from two Møller events) → determines beam position and detector transverse offset ($x$, $y$)
- The **vertex $z$ coordinate** (reconstructed from the transverse radii of the two hits in a single Møller event combined with beam energy) → constrains the detector longitudinal position ($z$)
- The **azimuthal angle difference** ($\Delta\phi$) between the same particle measured in GEM and HyCal within one Møller event → determines the GEM roll around the beam axis ($\phi$ alignment)

**Mott scattering** ($e^- p \to e^- p$): elastic electron-proton scattering producing a single outgoing electron. Mott events provide:
- Residuals between GEM chambers and between GEM and HyCal → used for relative GEM-to-GEM alignment

---

## 2. Detector Layout and Coordinate System

The system contains five detectors:

```
Beam direction →

  GEM2 (upstream left)    GEM0 (downstream left)
  GEM3 (upstream right)   GEM1 (downstream right)
                                          HyCal (calorimeter)
```

- **GEM numbering**: 0 = downstream left, 1 = downstream right, 2 = upstream left, 3 = upstream right
- **Coordinate system**: $z$ along the beam direction; facing the beam, $+x$ is to the right and $+y$ is upward (right-handed system, origin at the target)
- **Positions in config files**: stored as target-relative offsets. `LoadRunConfig` automatically subtracts the target offset on read, returning relative coordinates. In `main`, the target offset is added back before analysis to restore absolute coordinates.

---

## 3. Algorithm Overview

`det_align` is an **iterative** alignment program. Each run (one iteration) performs:

```
input run_config  →  transform events  →  fill histograms
→  fit peaks  →  compute corrections  →  output new run_config
```

The output `run_config_iter{N}.json` of iteration $N$ is fed as input to iteration $N+1$, and the corrections converge to zero over successive iterations.

---

## 4. Step-by-Step Algorithm

### 4.1 Event Selection

The program reads **reconstructed ROOT files** produced by `replay_recon` (filenames containing `_recon`, TTree named `recon`), which contain HyCal cluster reconstruction results and GEM matching information.

**Basic requirements**:
- Trigger type must be `TBIT_sum` (total-energy trigger)
- Number of HyCal clusters is 1 or 2, each cluster containing at least 3 crystals

**Møller events** (2 clusters) must satisfy:

| Condition | Description |
|-----------|-------------|
| $\theta_{1,2} > 0.65°$ | both electrons within HyCal acceptance |
| $\|\Delta\phi - 180°\| < 8°$ | azimuthal angles approximately back-to-back |
| $\|E_1 + E_2 - E_\text{beam}\| < 3\sigma$ | total energy conservation |
| $\|E_i - E_i^\text{expect}(\theta_i)\| < 3\sigma_i$ | each energy consistent with Møller kinematics |

**Mott events** (1 cluster + GEM match) must satisfy:

| Condition | Description |
|-----------|-------------|
| $\theta > 0.85°$ | within HyCal acceptance |
| $\|E - E_\text{expect}^\text{ep}(\theta)\| < 3\sigma$ | energy consistent with elastic $ep$ kinematics |
| exactly 1 matched hit per GEM chamber | no ambiguity |

### 4.2 Coordinate Transformation

After event selection, each HyCal cluster and GEM hit is re-expressed in the coordinate system of the current iteration's alignment parameters:

```
raw coordinates
  → ApplyToLocal(gRunConfig)    # transform to detector local frame using general.json
  → ApplyToLab(in_run_config)   # transform back to lab frame using current iteration config
```

This allows re-aligning hit positions without re-running the full replay.

### 4.3 Filling Alignment Histograms

#### Møller events: HyCal alignment

**Beam position / HyCal transverse alignment ($x$, $y$)**

For each Møller event, it is paired with the three preceding events. For each pair, the two hits of each Møller event define a line in the $xy$ plane; the **intersection** of these two lines is the Møller center:

$$\text{line}_i: y = a_i x + b_i,\quad a_i = \frac{y_i^{(1)} - y_i^{(2)}}{x_i^{(1)} - x_i^{(2)}},\quad b_i = y_i^{(1)} - a_i x_i^{(1)}$$

$$x_\text{center} = \frac{b_2 - b_1}{a_1 - a_2}, \quad y_\text{center} = a_1 x_\text{center} + b_1$$

The peak of the Møller center distribution gives the beam position in the detector coordinate system. If the beam is centered, the peak sits at 0; any offset is the correction $\Delta x$, $\Delta y$ to apply.

**Vertex $z$ alignment**

From a **single** Møller event, the transverse radii $R_1$, $R_2$ of the two hits at the detector face, combined with Møller kinematics, give the vertex $z$ coordinate:

$$R_i = \sqrt{x_i^2 + y_i^2}, \qquad z_\text{vertex} = \sqrt{\frac{(E_\text{beam} + m_e)\, R_1 R_2}{2\, m_e}}$$

The peak position directly constrains the $z$ coordinate of the scattering target relative to each detector.

**Upstream GEM beam center measurement**

Using `GEMup_moller` (Møller events reconstructed from the upstream GEM), the same Møller center method yields the beam position in the upstream GEM coordinate system.

#### Møller events: GEM alignment

For each GEM chamber $j$, the same Møller center and vertex $z$ calculations are performed using that chamber's hit coordinates:
- `h1_gem_CenterX[j]`: Møller center $x$ distribution measured by GEM $j$
- `h1_gem_CenterY[j]`: Møller center $y$ distribution measured by GEM $j$
- `h1_gem_Zdistance[j]`: vertex $z$ distribution reconstructed from GEM $j$

**GEM roll alignment ($\phi$ difference)**

For each particle ($k = 0, 1$) in the same Møller event, the azimuthal angle measured by HyCal and by GEM $j$ are compared:

$$\Delta\phi_j = \phi^\text{HyCal}_k - \phi^\text{GEM\_j}_k, \quad k=0,1$$

If GEM $j$'s coordinate system is rotated around the beam axis relative to HyCal (roll), the peak of $\Delta\phi$ will shift away from 0, and the shift gives the roll angle error.

#### Mott events: relative GEM alignment

Each GEM hit is projected to the HyCal $z$ plane:

$$x_\text{proj} = x_\text{GEM} \cdot \frac{z_\text{HyCal}}{z_\text{GEM}}$$

Residuals between chamber pairs are then computed:

| Histogram | Meaning |
|-----------|---------|
| `h1_deltaX/Y_gem_up` | GEM2 − GEM3 (same upstream layer; all 4 chambers matched) |
| `h1_deltaX/Y_gem_down` | GEM0 − GEM1 (same downstream layer; only $y > 0$ side) |
| `h1_deltaX/Y_gem_layer_left` | GEM0 − GEM2 (left-side inter-layer comparison) |
| `h1_deltaX/Y_gem_layer_right` | GEM1 − GEM3 (right-side inter-layer comparison) |
| `h1_deltaX/Y_gem_hycal[d]` | GEM $d$ − HyCal (each chamber vs. HyCal) |

### 4.4 Peak Extraction

The function `extract_peak(TH1F*)` extracts the peak center from a histogram:

1. Find the maximum bin
2. Use linear interpolation at 70% of the peak height to determine the left and right boundaries (FWHM-like window)
3. Perform a Gaussian fit within the window: $f(x) = A \exp\!\left(-\frac{(x-\mu)^2}{2\sigma^2}\right)$
4. Return the fitted mean $\mu$ (falls back to the maximum bin center if the fit fails)

### 4.5 Updating Detector Position Parameters

After peak extraction, the `RunConfig` is updated in the following order:

#### Step 0: Restore absolute coordinates

Config files store target-relative coordinates; the target offset is added back first:

```cpp
in_run_config.gem_x[det] += in_run_config.target_x;
in_run_config.gem_y[det] += in_run_config.target_y;
in_run_config.gem_z[det] += in_run_config.target_z;
```

#### Step 1: Beam position / target position ($x$, $y$)

The HyCal Møller center peak $(\Delta x, \Delta y)$ is the beam offset from the HyCal center:

```
target_x += peak(moller_HC_x)
target_y += peak(moller_HC_y)
```

If the target position has not yet been initialized (initial value is 0), it is set directly; otherwise the correction is added cumulatively.

#### Step 2: Detector $z$ position

The Møller vertex $z$ peak directly replaces the $z$ coordinate in the config:

```
hycal_z  = peak(moller_HC_z)
gem_z[d] = peak(moller_GEM{d}_z)
```

#### Step 3: GEM transverse position ($x$, $y$)

Corrections are applied at two levels (each using 50% of the measured value to prevent over-correction oscillation):

**3a. Global GEM offset relative to HyCal** (applied only when the beam position has converged, i.e. $|\Delta x|, |\Delta y| < 0.1$ mm):

```
gem_x[d] -= 0.5 × peak(moller_GEM{d}_x)
gem_y[d] -= 0.5 × peak(moller_GEM{d}_y)
```

**3b. Relative alignment between GEM chambers** (GEM3, the most upstream chamber, serves as the reference):

```
Align GEM2 to GEM3:  gem_x[2] -= 0.5 × peak(Upstream_GEM_dx_d2-d3)
Align GEM1 to GEM3:  gem_x[1] -= 0.5 × peak(GEM_layer_right_dx_d1-d3)
Align GEM0 to GEM2:  gem_x[0] -= 0.5 × peak(GEM_layer_left_dx_d0-d2)
(same for y)
```

> **Note**: Aligning GEM0 to GEM1 is currently commented out in the code, as the alignment chain GEM0→GEM2→GEM3 already covers it indirectly.

#### Step 4: GEM roll (rotation around $z$ axis)

```
gem_tilt_z[d] -= peak(moller_phi_diff_{d})
```

> **TODO**: Corrections for GEM tilt around the $x$ and $y$ axes (pitch/yaw) are not yet implemented.

---

## 5. Iterative Convergence

Each iteration produces:
- `run_config_iter{N}.json`: updated detector positions, used as input for the next iteration
- `alignment_summary_iter{N}.txt`: measured alignment parameter values for this iteration
- `alignment_histograms_iter{N}.root`: histograms and convergence graphs

**Convergence criterion**: inspect the `TGraph` objects under `Convergence/` in the ROOT file. When the correction (Delta) stabilizes near zero, iterations can stop. Typically 10–20 iterations are sufficient.

Sources for convergence graphs:
- HyCal $x$, $y$: `moller_HC_x/y` values read directly from `alignment_summary_iter{N}.txt` (each measurement is the residual offset)
- HyCal $z$: difference in `moller_HC_z` between consecutive iterations
- GEM $x$, $y$: absolute positions computed from consecutive `run_config_iter{N}.json` files; their difference is the actual correction applied each iteration
- GEM $z$, $\phi$: same as HyCal — consecutive differences or direct values from the summary

---

## 6. Output Files

Each run (iteration $N$) produces the following files:

### `run_config_iter{N}.json`

Updated run configuration containing corrected:
- `target` position (beam position)
- `hycal.position` (HyCal $x$, $y$, $z$)
- `gem[d].position` (each GEM chamber's $x$, $y$, $z$)
- `gem[d].tilting` (GEM roll angles)

### `alignment_summary_iter{N}.txt`

Measured alignment quantities for this iteration, formatted as:

```
Detector Alignment Parameters
==================================================
         moller_HC_x:  0.123 mm
         moller_HC_y: -0.045 mm
         moller_HC_z: 6234.567 mm
    moller_GEM0_x:  0.234 mm
    ...
```

### `alignment_histograms_iter{N}.root`

ROOT file containing the following directories:

| Directory | Contents |
|-----------|----------|
| `E_vs_Angle/` | Energy vs. polar angle 2D histograms (HyCal + each GEM) |
| `HyCal/` | HyCal hit positions, Møller center $x/y$, vertex $z$ |
| `GEM/` | GEM hit positions, Møller center, vertex $z$, $\phi$ difference |
| `Delta_GEMs/` | GEM inter-chamber residual histograms ($\Delta x$, $\Delta y$) |
| `Delta_GEM_vs_HyCal/` | GEM vs. HyCal residuals |
| `Convergence/` | Parameter convergence curves vs. iteration number (`TCanvas` + `TGraph`) |

---

## 7. Running the Program

### Build

```bash
cd build
make -j4 det_align
# executable at build/bin/prad2ana_det_align
```

### Basic usage

```bash
./bin/prad2ana_det_align <recon_data_dir_or_file> -o <output_dir> -i <iteration>
```

### Full example

```bash
# Iteration 1 (starting from general.json)
./bin/prad2ana_det_align /data/X17/recon/ -o ./align_output/ -i 1 -j 12

# Iteration 2 (automatically uses run_config_iter1.json)
./bin/prad2ana_det_align /data/X17/recon/ -o ./align_output/ -i 2 -j 12

# ... repeat until convergence

# Iteration N with a custom input config
./bin/prad2ana_det_align /data/X17/recon/ -o ./align_output/ -i 10 -j 12 \
    -c ./align_output/run_config_iter9.json
```

### Script for looped iterations

```bash
#!/bin/bash
RECON_DIR=/data/X17/recon/
OUT_DIR=./align_output/
EXE=./bin/prad2ana_det_align

for i in $(seq 1 20); do
    echo "===== Iteration $i ====="
    $EXE $RECON_DIR -o $OUT_DIR -i $i -j 12
done
```

---

## 8. Command-Line Options

| Option | Type | Default | Description |
|--------|------|---------|-------------|
| `<data>` | positional | required | reconstructed ROOT file(s) or directory (multiple allowed) |
| `-o <dir>` | string | required | output directory (created automatically) |
| `-i <N>` | int | 1 | iteration number; controls input/output file naming |
| `-j <N>` | int | 4 | number of parallel threads |
| `-c <file>` | string | `general.json` | input config file for this iteration (overrides auto-detection) |
| `-r <base>` | string | `<output_dir>/run_config.json` | base path for output config; written as `<base>_iter{N}.json` |
| `-f <N>` | int | −1 (all) | maximum number of files to process |
| `-n <N>` | int | −1 (all) | maximum number of events per file |

**Automatic input config resolution**:
- `-i 1`: uses `database/runinfo/general.json` (or the file given by `-c`)
- `-i N (N > 1)`: automatically uses `<base>_iter{N-1}.json`

---

## 9. Typical Workflow

### Checking convergence

After running several iterations, open the last ROOT output file:

```bash
root -l align_output/alignment_histograms_iter20.root
```

```cpp
// Browse interactively
new TBrowser

// Or draw directly
TFile f("alignment_histograms_iter20.root");
f.cd("Convergence");
HyCal_Convergence->Draw();   // HC_x, HC_y, HC_z
GEM0_Convergence->Draw();    // GEM0: x, y, z, phi
```

Convergence is reached when the Delta values no longer change appreciably (below ~0.1 mm).

### Using the final configuration

The converged `run_config_iter{N}.json` can be used as input for subsequent physics analysis replays:

```bash
replay_recon ... -c align_output/run_config_iter20.json
```

---

## 10. Known Limitations and Caveats

1. **Coordinate system convention**: detector coordinates in config files are target-relative; internally the program works in absolute coordinates, and subtracts the target offset when writing the output. Be careful when editing config files manually.

2. **50% damping factor**: GEM transverse corrections are multiplied by 0.5 to prevent over-correction oscillations, at the cost of requiring more iterations to converge.

3. **Prerequisite for GEM $x$, $y$ alignment**: the global transverse correction of each GEM relative to HyCal is only applied once the beam position has converged (`moller_HC_x/y` peak $< 0.1$ mm); before that, a global GEM translation is ill-defined.

4. **Missing iter 1 data point in GEM convergence graphs**: GEM $x$, $y$ convergence graphs are built by comparing consecutive `run_config_iterN.json` files. When running with iter $N > 1$, the path of the initial config used for iter 1 (`general.json`) is not recorded, so the iter 1 data point is absent; iter 2 and beyond appear normally.

5. **GEM tilt not yet implemented**: corrections for GEM rotation around the $x$ and $y$ axes (pitch/yaw) are not implemented (marked `// TODO` in the code).

6. **Parallel file processing**: files are processed concurrently by multiple threads, each writing to its own result slot. The main thread merges results in input-file order after all threads complete, ensuring deterministic output.
