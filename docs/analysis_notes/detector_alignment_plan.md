# PRad-II / X17 GEM Detector Precise Alignment Implementation Plan

## Executive Summary

This document outlines the detailed implementation strategy for the precision spatial alignment of the Gas Electron Multiplier (GEM) detector system for the **PRad-II / X17** experiment at Jefferson Lab (2.2 GeV electron beam, target located at $Z = 0$, GEM detectors positioned approximately 6 meters downstream of the target). 

The GEM detector system consists of **two layers** (Layer 1 and Layer 2) positioned sequentially along the beamline. Each layer is composed of **two separate chambers** (Left Chamber and Right Chamber), resulting in a total of **4 GEM chambers**. The Left and Right chambers overlap slightly in the central region around the beam pipe to minimize dead area and provide cross-calibration constraints.

The objective of this plan is to refine the GEM detector spatial alignment from the current coarse precision (~1 mm) down to the intrinsic detector spatial resolution limit (~72 $\mu$m). To achieve this, the alignment process utilizes pure $e^- e^-$ Møller elastic scattering kinematics and geometrical constraints, entirely bypassing reliance on calorimeter (HyCal) energy resolution limits.

---

## Geometry and Alignment Degrees of Freedom

The GEM detector system features 4 individual chambers:
* **Layer 1 Left Chamber (L1L)**
* **Layer 1 Right Chamber (L1R)**
* **Layer 2 Left Chamber (L2L)**
* **Layer 2 Right Chamber (L2R)**

Each chamber possesses **6 degrees of freedom (DoF)**:
1. **Translations:** $\Delta X$ (Horizontal), $\Delta Y$ (Vertical), $\Delta Z$ (Longitudinal / Beam direction)
2. **Rotations:** * **Pitch** ($\Delta\alpha$ / $\theta_x$): Rotation around the X-axis
   * **Yaw** ($\Delta\beta$ / $\theta_y$): Rotation around the Y-axis
   * **Roll** ($\Delta\gamma$ / $\theta_z$): Rotation around the Z-axis (Beam axis)

Total alignment parameters: $4 \text{ chambers} \times 6 \text{ DoF} = 24 \text{ global/local parameters}$.

---

## Alignment Strategy & Workflow Overview

The alignment follows a **"Global-First, Local-Second, Global-Fine"** sequential methodology using a standard 1D histogram slicing and Gaussian peak-fitting workflow (*Slice $\rightarrow$ Histogram $\rightarrow$ Fit Peak $\rightarrow$ Extract Parameter $\rightarrow$ Apply Correction*).

| Step | Scope | Target Parameters | Physical / Kinematic Constraint | Key Observables & Minimization Objective |
| :--- | :--- | :--- | :--- | :--- |
| **Step 1** | Global | $\Delta X_{global}, \Delta Y_{global}$, Global Roll ($\Delta\gamma$) | Møller scattering transverse spatial symmetry & coplanarity | Møller center aligned to $(0,0)$; $\Delta\phi$ peak centered at $180^\circ$ |
| **Step 2** | Global | Global $\Delta Z_{global}$ | Møller two-body elastic kinematics vertex reconstruction | Reconstructed Møller Vertex $Z$ peak centered at $Z = 0$ |
| **Step 3** | Local / Relative | Intra-layer & Inter-layer relative $\Delta X, \Delta Y, \Delta Z$, Pitch, Yaw | Overlap region hit coincidence & straight-line track projection | Residual slopes / parabolic coefficients vs. spatial coordinates driven to zero |
| **Step 4** | Global Fine | Global Pitch ($\Delta\alpha_{global}$), Global Yaw ($\Delta\beta_{global}$) | Møller polar angle projection dependence on azimuthal angle | Modulation amplitude of Vertex $Z$ vs. $\phi$ driven to zero |
| **Step 5** | System-wide | Parameter Convergence | Inter-parameter coupling minimization | Iterative cycle through Steps 1–4 until parameter variation $\Delta p < \text{Tolerance}$ |

---

## Detailed Step-by-Step Execution Plan

### Step 1: Global X, Y Translation and Global Roll Alignment

This step treats the entire 4-chamber GEM assembly as a rigid body to align its geometric center to the physics beamline $(0,0)$ and eliminate overall rotation around the beam axis.

#### 1.1 Global X and Y Translation ($\Delta X_{global}, \Delta Y_{global}$)
* **Event Selection:** Select clean, coincident double-electron Møller scattering events ($e^- e^- \rightarrow e^- e^-$).
* **Observable:** Reconstruct the geometric midpoint between the two Møller electron hit positions on the GEM plane (**Møller Center**):
  $$X_{center} = \frac{X_1 + X_2}{2}, \quad Y_{center} = \frac{Y_1 + Y_2}{2}$$
* **Histogram & Fitting:** 1. Fill 1D histograms for $X_{center}$ and $Y_{center}$.
  2. Perform Gaussian fits to extract the peak centers $X_{peak}$ and $Y_{peak}$.
* **Correction:** Apply global translation offsets $\Delta X_{global} = -X_{peak}$ and $\Delta Y_{global} = -Y_{peak}$ to align the Møller center to $(0,0)$.

#### 1.2 Global Roll Angle Alignment ($\Delta\gamma_{global}$)
* **Event Selection:** Coincident Møller scattering events.
* **Observable:** Azimuthal angle difference between the two outgoing Møller electrons:
  $$\Delta\phi = \vert{}\phi_1 - \phi_2\vert{}$$
  where $\phi_i = \arctan2(Y_i, X_i)$.
* **Histogram & Fitting:** 1. Fill a 1D histogram of $\Delta\phi$.
  2. Perform a Gaussian fit to extract the peak position $\Delta\phi_{peak}$.
* **Correction:** Rotate the entire GEM system around the Z-axis by $\Delta\gamma_{global}$ until $\Delta\phi_{peak} \equiv 180^\circ$ ($3.14159 \text{ rad}$).

---

### Step 2: Global Z Absolute Translation Alignment ($\Delta Z_{global}$)

This step determines the absolute physical distance along the longitudinal axis from the target center ($Z = 0$) to the GEM detector plane.

* **Event Selection:** Well-reconstructed Møller events.
* **Kinematic Principle:** For 2.2 GeV Møller scattering, the polar angles $\theta_1$ and $\theta_2$ measured by the GEM chambers uniquely determine the interaction vertex position $Z_{vertex}$ via two-body relativistic kinematics:
  $$\tan\theta_1 \cdot \tan\theta_2 = \frac{2 m_e}{E_0}$$
  where $E_0 = 2.2 \text{ GeV}$ is the beam energy and $m_e$ is the electron mass.
* **Observable:** Reconstructed longitudinal vertex position $Z_{vertex}$.
* **Histogram & Fitting:** 1. Fill a 1D histogram of $Z_{vertex}$.
  2. Fit the distribution with a Gaussian function (or Gaussian + polynomial background) to extract the peak $Z_{peak}$.
* **Correction:** Adjust the global Z position offset $\Delta Z_{global} = -Z_{peak}$ to anchor the GEM system to $Z = 0$.

---

### Step 3: Local Chamber Relative Alignment (Intra-Layer & Inter-Layer)

This step breaks the rigid-body assumption and aligns individual chambers relative to each other to ensure intra-layer coplanarity and inter-layer parallelism.

#### 3.1 Intra-Layer Alignment (Left Chamber vs. Right Chamber)
* **Region of Interest:** The central **Overlap Region** where Left and Right chambers overlap near the beam pipe.
* **Methodology (Y-Slicing for Relative Pitch):**
  1. Select tracks passing through the overlap region leaving hits in both Left and Right chambers within Layer 1 (or Layer 2).
  2. Divide the overlap region along the **Y-axis** into $N$ equal slice bins ($Y_1, Y_2, \dots, Y_N$).
  3. In each Y-slice bin $k$, plot the position residual histogram:
     $$\Delta Y_k = Y_{Left, k} - Y_{Right, k} \quad (\text{or } \Delta X_k = X_{Left, k} - X_{Right, k})$$
  4. Perform Gaussian fitting on each slice to extract the peak residual value $\text{Peak}(\Delta Y)_k$.
  5. Construct a profile plot of **$\text{Peak}(\Delta Y)$ vs. $Y$**.
* **Diagnostic & Parameter Extraction:**
  * If Left and Right chambers have relative Pitch ($\Delta\alpha_{Left} - \Delta\alpha_{Right} \neq 0$), the profile plot displays a linear slope:
    $$\Delta Y(Y) = a_0 + a_1 \cdot Y$$
  * Slope $a_1$ directly yields the relative Pitch angle.
* **Methodology (X-Slicing for Relative Yaw):**
  * Slicing along the X-axis in the overlap region yields relative Yaw ($\Delta\beta_{Left} - \Delta\beta_{Right}$).
* **Correction:** Adjust individual chamber Pitch and Yaw parameters until the profile slopes $a_1 \rightarrow 0$ and offsets $a_0 \rightarrow 0$.

#### 3.2 Inter-Layer Alignment (Layer 1 vs. Layer 2)
* **Region of Interest:** Full active acceptance of the GEM detector.
* **Methodology (Straight-Line Track Projection):**
  1. Assume tracks originate from the target center $(0,0,0)$ and travel in straight lines through Layer 1 and Layer 2 (magnetic field-free region).
  2. Project the trajectory from Layer 1 hit position $(X_1, Y_1, Z_1)$ to Layer 2 nominal position $Z_2$:
     $$X_{proj} = X_1 \cdot \frac{Z_2}{Z_1}, \quad Y_{proj} = Y_1 \cdot \frac{Z_2}{Z_1}$$
  3. Calculate projection residuals on Layer 2:
     $$\delta X = X_2 - X_{proj}, \quad \delta Y = Y_2 - Y_{proj}$$
  4. Slice Layer 2 along the **X-axis** into $M$ bins. In each bin, fit the 1D histogram of $\delta X$ to obtain $\text{Peak}(\delta X)_m$.
  5. Plot **$\text{Peak}(\delta X)$ vs. $X$**.
* **Diagnostic & Parameter Extraction:**
  * Relative Yaw ($\Delta\beta_{rel}$) between Layer 1 and Layer 2 produces a characteristic **quadratic (parabolic) trend** due to projection geometric expansion:
    $$\delta X(X) \approx c_0 + c_1 X + c_2 X^2$$
  * The quadratic coefficient $c_2$ is proportional to the relative Yaw angle.
  * Similarly, slicing along **Y-axis** and plotting **$\text{Peak}(\delta Y)$ vs. $Y$** extracts relative Pitch via quadratic term fitting.
* **Correction:** Rotate Layer 2 relative to Layer 1 until the parabolic coefficients $c_2 \rightarrow 0$.

---

### Step 4: Global Pitch and Yaw Fine Alignment ($\Delta\alpha_{global}, \Delta\beta_{global}$)

With individual chambers aligned internally and planar parallelized, Step 4 corrects any residual global inclination of the entire GEM plane vector relative to the beam line.

* **Event Selection:** Coincident Møller events.
* **Physical Effect of Global Tilt:**
  * If the GEM plane is tilted relative to the Z-axis, one side of the detector is closer to the target than the nominal $Z$, while the opposite side is farther.
  * Reconstructed polar angles $\theta$ become artificially enlarged on the closer side and compressed on the farther side.
  * Consequently, the reconstructed Møller interaction vertex $Z_{vertex}$ exhibits a periodic **sinusoidal modulation** as a function of the electron azimuthal angle $\phi$.

* **Methodology ($\phi$-Slicing):**
  1. Divide the full azimuthal coverage $\phi \in [0^\circ, 360^\circ]$ into $K$ uniform bins (e.g., $15^\circ$ per bin).
  2. In each $\phi$-slice bin $k$, fill a 1D histogram of the reconstructed Møller $Z_{vertex}$.
  3. Perform a Gaussian fit to extract the peak vertex position $\text{Peak}(Z_{vertex})_k$.
  4. Construct a profile plot of **$\text{Peak}(Z_{vertex})$ vs. $\phi$**.
* **Fitting & Interpretation:**
  * Fit the profile graph with a sine function:
    $$Z_{vertex}(\phi) = Z_0 + A \cdot \sin(\phi + \phi_0)$$
  * **Amplitude $A$:** Directly proportional to the global tilt magnitude $\sqrt{\Delta\alpha_{global}^2 + \Delta\beta_{global}^2}$.
  * **Phase $\phi_0$:** Determines the tilt direction (Phase at $0^\circ / 180^\circ$ corresponds to Global Yaw $\Delta\beta$; Phase at $90^\circ / 270^\circ$ corresponds to Global Pitch $\Delta\alpha$).
* **Correction:** Apply global rotations $\Delta\alpha_{global}$ and $\Delta\beta_{global}$ to flatten the sine wave modulation into a constant horizontal line at $Z_{vertex} = 0$.

---

### Step 5: Iterative Cycle & Convergence Criteria

Due to weak non-linear geometric coupling between spatial translation and angular orientation parameters (e.g., modifying global Pitch induces second-order shifts in outer boundary X/Y positions), steps 1 through 4 must be executed in an iterative loop.

```text
       +-----------------------------------------------+
       | Start Iteration Loop (Iteration N = 1, 2, ...) |
       +-----------------------------------------------+
                               |
                               v
       +-----------------------------------------------+
       | Step 1: Global X, Y Translation & Global Roll |
       +-----------------------------------------------+
                               |
                               v
       +-----------------------------------------------+
       | Step 2: Global Z Absolute Translation         |
       +-----------------------------------------------+
                               |
                               v
       +-----------------------------------------------+
       | Step 3: Local Chamber Relative Alignment      |
       |         (3.1 Intra-layer / 3.2 Inter-layer)   |
       +-----------------------------------------------+
                               |
                               v
       +-----------------------------------------------+
       | Step 4: Global Pitch & Yaw Fine Alignment     |
       +-----------------------------------------------+
                               |
                               v
       +-----------------------------------------------+
       | Check Convergence Criteria:                   |
       | - |Delta X, Y, Z| < 10 um                       |
       | - |Delta alpha, beta, gamma| < 0.1 mrad       |
       | - Profile residuals flattened                 |
       +-----------------------------------------------+
               /                               \
        [No]  /                                 \ [Yes]
             v                                   v
   Apply Corrections to                     ALIGNMENT
   Reconstruction Parameters                 COMPLETE
   & Rerun Tracking                         Save Geometry DB