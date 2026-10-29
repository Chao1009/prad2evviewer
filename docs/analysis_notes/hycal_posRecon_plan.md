# HyCal Position Reconstruction Research Plan

## 1. Objectives and Research Sequence

Compare HyCal cluster position reconstruction methods while separating position bias, energy response, and tracking uncertainties. Identify an algorithm suitable for offline reconstruction with a feasible path to online use.

Proceed in the following order, using the same evaluation protocol at each stage and comparing against the best method from the preceding stage:

1. **Log-weighted centroid baseline**: reproduce existing results and scan the weighting parameter and position reconstruction window.
2. **Log-weighted centroid with two-dimensional position corrections**: reduce intra-module S-curve bias.
3. **Shower template fitting**: study boundaries, dead channels, mixed materials, and overlapping showers.
4. **Machine learning regression**: evaluate additional gains using the datasets and benchmarks established in the preceding stages.

Use the maximum-energy module center and the linear energy-weighted centroid as controls in Stage 1, and neighboring-module energy ratios / eta methods as auxiliary controls in Stage 2. Neither a more uniform position density nor a smoother energy-response bowl is sufficient on its own to establish better position reconstruction.

This document is a research plan; it does not imply that the algorithms, scan interfaces, or output formats described below are already implemented. Do not replace the production default algorithm during the research phase.

## 2. Current Implementation and Reusable Components

| Component | Current role and research use |
|---|---|
| [HyCalCluster.cpp](../../prad2det/src/HyCalCluster.cpp) | `reconstruct_pos` computes a log-weighted centroid using the seed and a geometrically defined 3x3 neighborhood. Weights use the total cluster energy as the denominator. A leakage-corrected position, when available, overrides the initial position |
| [HyCalCluster.h](../../prad2det/include/HyCalCluster.h) | `ClusterConfig::log_weight_thres` defaults to 3.6; provides interfaces for clusters, module energies, and shower profiles |
| [MatchingTools.h](../../analysis/include/MatchingTools.h) | HyCal/GEM coordinate transformation, projection, and matching interfaces; `ClusterToLab` places the cluster at shower depth by default |
| [hycal_cluster_density.cpp](../../analysis/tools/hycal_cluster_density.cpp) | Reuse existing position distributions, residual slicing, and peak extraction after checking their selections and correction definitions |
| [hycal_energy_bowl_shape.cpp](../../analysis/tools/hycal_energy_bowl_shape.cpp) | Currently bins energy histograms by reconstructed HyCal position, supporting 5x5 and 10x10 grids. Also outputs intra-module HyCal/GEM position distributions, but does not yet compare energy responses binned by both position estimates |
| [REPLAYED_DATA.md](../REPLAYED_DATA.md) | Check input ROOT branches, coordinate definitions, and matching indices |
| [analysis/README.md](../../analysis/README.md) | Reuse the build and execution workflows for replay, filtering, and analysis tools |

The seed module ID is not a continuous reconstructed position and does not necessarily identify the module traversed by the particle. At interfaces between modules of different sizes, first reproduce the existing normalized-distance treatment rather than assuming it is equivalent to an ordinary centroid in physical coordinates.

## 3. Prerequisites: Common Samples and Reference Positions

### 3.1 Data and Event Selection

- [ ] Select runs with reliable GEM information, stable gains, and known geometry. Record beam energy, triggers, calibration versions, and bad-channel maps.
- [ ] Start with isolated single showers in the PbWO4 interior and valid, unambiguous GEM references. Then extend to PbGlass, material interfaces, the central hole, outer edges, and regions near dead channels.
- [ ] Treat the existing bowl tool's single-cluster requirement, minimum module count, energy window, and annular acceptance as selections to evaluate, not as inherently unbiased criteria.
- [ ] Check how module-count, energy, and matching-window cuts affect intra-module position distributions. Save a cut-flow for each selection.
- [ ] Split training, validation, and final test samples by run / independent acquisition period, preventing the same event or duplicated source data from appearing in multiple sets. Do not use the final test set for tuning.
- [ ] Determine statistics in each energy and spatial region. Mark low-statistics bins invalid or merge them using predefined rules; do not fill missing estimates with zero.

### 3.2 GEM References and Coordinate Conventions

- [ ] Define the primary evaluation quantity as the position at the effective shower depth, and separately provide results projected to the HyCal front face. Do not mix the two reference planes.
- [ ] First reproduce the current projection chain between `cl_x/cl_y/cl_z` and GEM, avoiding double application of the shower-depth correction.
- [ ] Distinguish a straight-line projection from the target origin through a single GEM hit from a track fitted to multiple GEM measurements. The existing `GetProjection` implements the former and must not be described as an independently measured track direction.
- [ ] Where the data allow, use GEM tracks without a HyCal position constraint. Otherwise, state the target-position assumption and its uncertainty.
- [ ] Determine global translations, rotations, and depth-dependent biases using an independent sample, then freeze the geometry parameters. Do not allow local position corrections to absorb global alignment errors.
- [ ] Check validity flags and cluster indices when matching. Missing GEM information must not become a valid zero coordinate. Count unmatched and ambiguous events.
- [ ] Compare algorithms using common reference samples and sufficiently broad matching criteria. Scan the matching window to ensure that narrow residual cuts do not artificially improve the measured resolution.

GEM provides a reference with measurement uncertainty, not an exact truth position. Consider subtracting tracking and multiple-scattering contributions from the residual variance only when errors are approximately independent, positions refer to the same projection plane, and the covariance model is valid. Otherwise, report only the HyCal-GEM residual resolution. Photon position performance requires simulation truth or an independent photon sample; charged-particle GEM samples are not a direct substitute.

### 3.3 Reconstruction Inputs and Reproducibility

- [ ] Check whether existing recon files retain the per-module energies, times, cluster membership, and split fractions required for recomputation. Final `cl_x/cl_y` values alone are insufficient for scanning position algorithms.
- [ ] If inputs are insufficient, reuse the existing replay chain on raw data to export research inputs or calculate candidate positions during replay. Do not infer module energies from the total energy.
- [ ] In the first comparison, freeze clustering, seeds, membership, and energy sharing, replacing only the position estimator. Subsequently test end-to-end effects coupled to cluster splitting / leakage iterations separately.
- [ ] Save uncorrected positions, all energy correction factors, and configurations. Label positions before and after leakage correction separately.
- [ ] Record the code version, input manifest, configuration, sample split, random seed, and output paths for every run.

## 4. Common Metrics and Stage Acceptance Criteria

Define residuals on the same reference plane:

$$
\Delta x=x_{\mathrm{HyCal}}-x_{\mathrm{GEM,proj}},\qquad
\Delta y=y_{\mathrm{HyCal}}-y_{\mathrm{GEM,proj}}.
$$

Define local coordinates normalized by module width:

$$
u=(x-x_{\mathrm{module}})/d_x,\qquad
v=(y-y_{\mathrm{module}})/d_y.
$$

**Prefer binning diagnostic plots by the GEM reference position**, while retaining versions binned by reconstructed position to identify migration and selection-induced biases. Record the seed module and the module containing the projected GEM position separately; do not require them to be identical.

Report at least the following for each algorithm:

| Metric | Definition / checks |
|---|---|
| Position bias | Mean residual and uncertainty in each bin; provide the median or main-peak center as auxiliary estimates |
| Robust resolution | `sigma68 = (q84 - q16) / 2`, reported separately for x/y |
| Fitted resolution | Core Gaussian sigma, fit status, and fit range using fixed fitting rules |
| Tails | Fraction exceeding a residual threshold established from the baseline and then frozen; do not redefine the threshold for each algorithm |
| Two-dimensional correlations | Residual covariance, two-dimensional bias maps, and continuity across module boundaries |
| Efficiency and failure rate | Total candidate events, valid positions, matching ambiguities, fit failures, and fallback counts |
| Resources | Median and p95 time per cluster, and peak memory, using identical hardware and inputs |

Group results by energy, material, intra-module position, incidence angle, distance to boundaries, module count, and run. Report results under both actual illumination weights and uniform spatial weights so that global averages do not hide local degradation.

**Suggested initial acceptance targets**: on the validation set, relative to the best method from the preceding stage, reduce the spatial bias RMS by at least 20% or sigma68 by at least 5%. In the main valid regions, sigma68 must not degrade by more than 2%, and the absolute increases in the fixed-threshold tail fraction and failure rate must not exceed 0.5 percentage points. Define spatial bias RMS as the square root of the fixed-weight mean of the squared mean residuals in valid local bins, reported separately for x/y.

These numbers are research targets, not established performance results. At the end of Stage 1, freeze the final thresholds, valid-bin requirements, and resource budget based on statistical sensitivity before examining the final test set. Use paired-event comparisons and run-grouped bootstrap estimates of uncertainty; do not claim statistically insignificant changes as improvements. Algorithms intended only for specific regions, such as boundaries, may use predeclared regional targets but must report their coverage.

## 5. Stage 1: Log-Weighted Centroid Baseline and Parameter Scans

### Method

$$
w_i=\max(0,W_0+\ln(E_i/E_{\mathrm{ref}})),\qquad
\hat{\boldsymbol r}=\frac{\sum_i w_i\boldsymbol r_i}{\sum_i w_i}.
$$

The expression above is the standard form for regions with equal-sized modules. Reproduction of the current implementation must also cover its geometric distances, neighborhood selection, positive-weight counts, and fallback behavior when no weight remains.

### Tasks

- [ ] Reproduce positions event by event under the current configuration, identifying differences between defaults and the actual run configuration.
- [ ] Add the maximum-energy module center and linear energy-weighted centroid as controls. Optionally test a small number of `E_i^alpha` weights without an unbounded parameter search.
- [ ] Start with `W0 = 3.0, 3.2, ..., 5.4`, then perform a finer scan near the validation-set optimum. This is an initial range, not an assumption that the optimum lies within it.
- [ ] First fix the current 3x3 window and total-cluster-energy denominator, scanning only W0.
- [ ] Then compare 3x3, 5x5, and full-cluster windows, explicitly defining membership and boundaries.
- [ ] Compare total cluster energy and window energy as separate denominator conventions, retuning W0 for each.
- [ ] Check how W0 affects `npos`, single-module dominance, low-energy thresholds, and the handling of nonpositive energies.
- [ ] Test a constant W0 first. Introduce simple energy / material dependence only if clear energy-dependent biases remain, avoiding per-module overfitting.
- [ ] Perform separate leakage-on/off ablation tests so that the final position override does not obscure the weight-scan effects.

### Deliverables and Completion Criteria

- A default-configuration reproduction report, including event-by-event differences and numerical tolerances.
- W0 / window / denominator scan results, local residual maps, energy dependence, and resource costs.
- A frozen best log-weighted centroid baseline, with training, validation, and test manifests shared by subsequent stages.
- If no change significantly outperforms the existing configuration, retain it as the baseline rather than forcing a new parameter choice.

## 6. Stage 2: Log-Weighted Centroid with Two-Dimensional Position Corrections

### Method

Use the best Stage 1 position as input and calibrate corrections that can be evaluated during actual reconstruction:

$$
x_{\mathrm{corr}}=x_{\mathrm{raw}}+\Delta x(u_{\mathrm{raw}},v_{\mathrm{raw}},E,\ldots),
$$

$$
y_{\mathrm{corr}}=y_{\mathrm{raw}}+\Delta y(u_{\mathrm{raw}},v_{\mathrm{raw}},E,\ldots).
$$

Projected GEM positions may serve as training labels, but application-time inputs must not depend on unknown truth positions or GEM labels. Bias curves plotted against reference positions cannot be used directly as correction tables in reconstructed coordinates without an inversion.

### Tasks

- [ ] Start with independent one-dimensional S-curve corrections for x/y as simple controls.
- [ ] Then compare two-dimensional lookup tables and smooth two-dimensional parameterizations, focusing on x-y coupling near module corners.
- [ ] Initially share corrections by material and geometric region. Introduce per-module terms only where statistics are sufficient and differences are established.
- [ ] Define the validity domain, interpolation rules, and low-statistics handling. Outside the domain, retain the baseline and record an explicit status rather than extrapolating without constraints.
- [ ] Check continuity across seed changes and module boundaries. Do not simply clip genuine cross-module events back into the seed module.
- [ ] Add left/right and top/bottom energy asymmetries or eta mappings as auxiliary controls, checking for insufficient neighboring signals near module centers.
- [ ] Treat density / cumulative-distribution corrections as auxiliary studies only. Validate them using reference illumination and acceptance models rather than assuming a uniform true distribution.
- [ ] Test generalization on independent runs, energy intervals, and modules excluded from calibration.

### Cross-Checks with the Energy-Response Bowl

- [ ] For the same events with valid GEM references, construct 5x5 and 10x10 energy grids using the uncorrected HyCal, corrected HyCal, and projected GEM positions separately.
- [ ] Compare module and grid-cell migration matrices, per-cell energy peaks, fit status, and statistics.
- [ ] Initially freeze the energy estimator and sample selection to isolate position-binning effects. Recalibrate position-dependent energy corrections afterward.
- [ ] Independently validate position and bowl energy calibrations. Check the stability of position-energy iterations to avoid apparent improvements caused by compensating errors.

### Deliverables and Completion Criteria

- Independent test comparisons of one- and two-dimensional corrections, with versioned correction parameters and validity domains.
- Bowl comparisons and migration reports for all three coordinate-based binnings.
- Meet the common acceptance criteria without significant boundary discontinuities; otherwise retain the log-weighted centroid baseline.

## 7. Stage 3: Shower Template Fitting

### Method

Model the expected response in each module:

$$
\mu_i=E_0 F_i(x_0,y_0,E_0,\theta,\ldots),
$$

Fit position and energy by minimizing a negative log-likelihood or

$$
\chi^2=(\boldsymbol E-\boldsymbol\mu)^T
V^{-1}(\boldsymbol E-\boldsymbol\mu)
$$

The template should describe the response integrated over the finite module area, rather than only the value of a radial function at the module center.

### Tasks

- [ ] Audit the variables, units, normalization, material coverage, and energy coverage of the existing `IClusterProfile` / `Geant4Profile`. Do not assume that profiles used for splitting are sufficient for precision position fitting.
- [ ] Start with closure tests on single showers in the PbWO4 interior, then extend to PbGlass and irregular boundaries.
- [ ] Use Stage 2 positions as initial values. Specify whether the energy parameter is free, which modules enter the fit, and the template's validity range.
- [ ] Begin with interpretable noise and shower-fluctuation models and test inter-module correlations. Do not impose a Poisson energy likelihood where the signals do not follow simple Poisson statistics.
- [ ] Distinguish below-threshold channels, unread channels, dead channels, and geometrically missing modules; do not treat them all as precisely measured zero energy.
- [ ] Let the model handle missing response at dead channels and boundaries directly, avoiding double counting with existing leakage energy corrections.
- [ ] Save convergence status, objective-function values, uncertainty estimates, and parameter-boundary flags. On failure, use an explicitly flagged baseline fallback and include it in failure rates and overall performance.
- [ ] After validating single-shower fitting, study joint multi-shower fits, separately evaluating positions, separation efficiency, false-splitting rate, and energy sharing.
- [ ] Evaluate systematic uncertainties from simulation/data differences, gain perturbations, and template choices.

### Deliverables and Completion Criteria

- Template generation and validation workflows, validity ranges, fit-quality results, and failure diagnostics.
- Regional gains and timing reports relative to Stage 2, including fallback events rather than only converged fits.
- A decision, based on measured gains, to apply the method globally, enable it only in special regions, or retain it as an offline cross-check.

## 8. Stage 4: Machine Learning Regression

### Tasks

- [ ] Start with a small BDT / MLP using local energy fractions, total energy, module geometry, material, and valid-channel masks.
- [ ] Distinguish dead channels, missing channels, and measured zero energy. Use the actual geometry rather than forcing mixed-size modules into an equally spaced image.
- [ ] Predict the position relative to the seed or the residual correction to the Stage 2 result. Clearly define the target plane and label uncertainties.
- [ ] Evaluate CNN / GNN models only if simpler models show clear limitations; model complexity is not an objective.
- [ ] Use GEM labels for data and truth labels for simulation, evaluating them separately. Quantify tracking-label noise, particle-type dependence, and simulation domain shifts.
- [ ] Exclude GEM positions, label-derived variables, and event identity information from inference features. If incidence angle is used, specify a source that is actually available at inference time.
- [ ] Use held-out runs, energies, and modules to prevent random event splits from hiding gain and geometry dependence.
- [ ] Test dead channels, gain drift, noise, low energies, boundaries, and out-of-distribution samples. Define how unsupported inputs are identified and handled with status-bearing fallbacks.
- [ ] If the model predicts uncertainties, check coverage on an independent test set. Record the model version, training configuration, feature units, and inference resources.

### Deliverables and Completion Criteria

- Reproducible training configurations and frozen models, with a common test table comparing Stages 1 through 3.
- Explicit model validity domains, systematic uncertainties, and deployment costs.
- Consider replacing simpler methods only if gains are significant, stable, and interpretable; otherwise retain the model as a research control.

## 9. Final Integration and Validation

Follow the dependency order "common evaluation inputs -> Stage 1 -> Stage 2 -> Stage 3 -> Stage 4." Freeze the best candidate and comparison results after each stage. Completing the research sequence does not require deploying a complex algorithm.

- [ ] Enable new methods through explicit algorithm configuration while keeping the old default reproducible. Prefer existing geometry, projection, profile, and calibration interfaces.
- [ ] Check how positions enter splitting, leakage, energy-bias corrections, shower depth, and GEM matching. Report fixed-membership algorithm comparisons separately from full replay results.
- [ ] Verify consistency between offline replay and online reconstruction, including coordinates, correction factors, and algorithm metadata in ROOT outputs.
- [ ] Add targeted tests for symmetric showers, known offsets, a single valid module, zero weights, threshold boundaries, dead channels, material interfaces, seed changes, fit failures, and fallbacks.
- [ ] Verify agreement with previous results when new corrections are disabled. When enabled, validate finite outputs, domains, units, position continuity, and status propagation.
- [ ] Perform the final end-to-end comparison on independent test runs and update directly related configuration and analysis documentation.

Final deliverables: a common metrics table for all algorithms; the best configuration, calibration, or model; input manifests for reproducibility; local bias and resolution maps; bowl cross-check results; failure and systematic-uncertainty reports; and an explicit recommendation on whether to replace the production default.
