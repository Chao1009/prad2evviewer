# HyCal Cluster Density and S-Shape Corrections

This document describes how cluster-density parameters are loaded and how the position and energy corrections are applied during HyCal reconstruction in the PRad1 software "PRadAnalyzer".

## Overview

The implementation contains two related corrections:

- **Position density correction**: removes the S-shaped bias in the reconstructed cluster position.
- **S-shape energy correction**: removes the position-dependent bias in the reconstructed cluster energy for elastic ep and Moller ee events.

The processing flow is:

```text
Configuration files
    -> PRadHyCalSystem loads both correction profiles
    -> PRadClusterDensity stores the parameters
    -> PRadHyCalReconstructor selects the profile set for the run
    -> Each reconstructed cluster calls CorrectBias()
    -> Cluster position and energy are updated
```

The main implementation is in:

- [`PRadClusterDensity.cpp`](../lib/prana/src/PRadClusterDensity.cpp)
- [`PRadClusterDensity.h`](../include/PRadClusterDensity.h)
- [`PRadHyCalReconstructor.cpp`](../lib/prana/src/PRadHyCalReconstructor.cpp)
- [`PRadHyCalSystem.cpp`](../lib/prana/src/PRadHyCalSystem.cpp)

## Configuration

The profile paths are configured in [`config/hycal.conf`](../config/hycal.conf):

```text
Density Profile [Set_1GeV] = ${DB_DIR}/density_params/set_1GeV.dat
Density Profile [Set_2GeV] = ${DB_DIR}/density_params/set_2GeV.dat

S-shape Energy Profile [Set_1GeV] = ${DB_DIR}/s_energy_params/ecorrect_1GeV.dat
S-shape Energy Profile [Set_2GeV] = ${DB_DIR}/s_energy_params/ecorrect_2GeV.dat
```

The correction switches are defined in [`config/hycal_cluster.conf`](../config/hycal_cluster.conf):

```text
Density Correction = true
S-shape Energy Correction = true
```

These switches become `config.den_corr` and `config.sene_corr` in the reconstructor.

## Loading the Parameters

During HyCal system initialization, `PRadHyCalSystem` loops over all `PRadClusterDensity::SetEnum` values and calls `LoadDensityParams()` once for each set. The reconstructor forwards the call to:

```cpp
PRadClusterDensity::Load(int set, const std::string& position_path,
                         const std::string& energy_path)
```

Each set is stored in one element of:

```cpp
ParamsSet psets[Max_SetEnums];
```

A `ParamsSet` contains:

```cpp
float beam_energy;
std::vector<float> energy_range;
std::vector<Params> ppars;
std::unordered_map<int, Params> epars_ep;
std::unordered_map<int, Params> epars_ee;
```

`Load()` reads the two files independently:

1. `processPosPars()` reads the position-density file.
2. `processEnePars()` reads the S-shape energy file.

## Position-Density File

The position files are [`database/density_params/set_1GeV.dat`](../database/density_params/set_1GeV.dat) and [`database/density_params/set_2GeV.dat`](../database/density_params/set_2GeV.dat).

Their header contains values such as:

```text
number_of_params = 4
geometry_groups = 112
energy_groups = 5
beam_energy = 1097
energy_range = 0, 250, 750, 1500
```

The fields mean:

- `geometry_groups`: number of detector geometry groups.
- `energy_groups`: number of position-parameter energy groups.
- `beam_energy`: beam energy used when identifying elastic and Moller events.
- `energy_range`: boundaries for ordinary energy bins.
- `number_of_params`: number of position coefficients read by the implementation for each coordinate parameter set.

The in-memory position container is:

```cpp
struct Params {
    std::vector<float> x, y;
};
```

The file contains one block for x-position parameters followed by one block for y-position parameters. The parameters are stored in a flat vector whose logical index is:

```text
position_parameter_index = geometry_index * energy_group_count + energy_index
```

The loader checks that:

```text
energy_groups == energy_range.size() + 1
```

The extra group is needed because the implementation reserves two special groups after the ordinary energy ranges: one for elastic ep events and one for Moller ee events. For example, with three range boundaries, the file declares five groups: three ordinary groups, one ep group, and one ee group.

## S-Shape Energy File

The energy files are [`database/s_energy_params/ecorrect_1GeV.dat`](../database/s_energy_params/ecorrect_1GeV.dat) and [`database/s_energy_params/ecorrect_2GeV.dat`](../database/s_energy_params/ecorrect_2GeV.dat).

They contain two sections:

```text
EP_PARAMS
<module name> <8 coefficients>
...

EE_PARAMS
<module name> <8 coefficients>
...
```

The loader stores the two sections separately:

- `EP_PARAMS` -> `epars_ep`
- `EE_PARAMS` -> `epars_ee`

The key is the HyCal module id obtained from the module name. Unlike position parameters, energy parameters are looked up by the cluster-center module id, not by a geometry-group index.

## Selecting the Active Profile Set

Both the 1 GeV and 2 GeV parameter sets are loaded, but only one is active at reconstruction time. The current selection is based on the run number in [`PRadHyCalSystem.cpp`](../lib/prana/src/PRadHyCalSystem.cpp):

```cpp
if(run < 1362)
    recon.ChooseDensitySet(PRadClusterDensity::Set_1GeV);
else
    recon.ChooseDensitySet(PRadClusterDensity::Set_2GeV);
```

Therefore:

```text
run < 1362  -> Set_1GeV
run >= 1362 -> Set_2GeV
```

The set is **not** selected separately for each cluster based on the cluster energy. Cluster energy is only used later to select an energy group within the already-selected set.

## Reconstruction and Correction Order

For each cluster, `PRadHyCalReconstructor::ReconstructHits()` performs the following sequence:

```text
FormCluster()
    -> CheckCluster()
    -> LeakCorr()
    -> Cluster2Hit()
```

Inside `Cluster2Hit()` the code first:

1. Creates a `HyCalHit` from the cluster.
2. Reconstructs the cluster position.
3. Applies the non-linearity energy correction when enabled.
4. Applies the shower-depth correction to `z`.
5. Calculates energy and position resolutions.
6. Calls `density.CorrectBias()`.

Thus, density correction operates on the already reconstructed `HyCalHit`, after leakage and non-linearity corrections.

## Energy-Group Lookup

`getEnergyIndex()` receives:

```text
cluster energy
polar angle
6 * hit.sig_ene
```

The lookup priority is:

1. Check whether the energy is consistent with the elastic ep energy.
2. Check whether the energy is consistent with the Moller ee energy.
3. Otherwise, find the ordinary energy range from `energy_range`.

The special event energies are calculated from the selected set's `beam_energy`:

```text
E_ep = mott_energy(beam_energy, theta)
E_ee = moller_energy(beam_energy, theta)
```

The tolerance is relative and is compared against `6 * hit.sig_ene`.

If no range matches, `getEnergyIndex()` returns `-1`, and no position correction is applied.

## Geometry-Group Lookup

`getGeometryIndex()` maps the cluster-center module into a shared detector geometry group.

- PbWO4 modules are grouped using their row and column, approximately in 2x2 regions with symmetry handling.
- PbGlass modules are grouped using the detector's 30x30 layout, approximately in 3x3 regions with boundary-region handling.

The geometry grouping reflects how the density correction study was performed: neighboring modules share the same fitted position parameters.

## Position Correction

When position correction is enabled and the hit does not already have `kDenCorr`, `CorrectBias()` computes:

```text
dx = (hit.x - center.x) / center.module_size_x
dy = (hit.y - center.y) / center.module_size_y
```

It then selects:

```text
geometry_index = getGeometryIndex(cluster_center)
energy_index   = getEnergyIndex(hit.E, theta, 6 * hit.sig_ene)
parameter_index = geometry_index * 5 + energy_index
```

The factor `5` is the configured number of energy groups in the current profiles.

The x and y biases are calculated independently with the same polynomial form and different coefficients:

```text
bias(d)
  = d * (d^2 - 0.25)
    * c0
    * (d^4 + c1*d^2 + c2)
    * (d^2 - c3)
```

The position is updated by adding the calculated bias back to the reconstructed position:

```text
hit.x += bias_x(dx) * center.module_size_x
hit.y += bias_y(dy) * center.module_size_y
```

After a successful position correction, the code sets `kDenCorr` so the same hit cannot be corrected twice.

## S-Shape Energy Correction

Energy correction is attempted when enabled and the hit does not already have `kSEneCorr`.

Only the two special energy groups are eligible:

- The elastic ep group uses `epars_ep`.
- The Moller ee group uses `epars_ee`.

The module lookup uses the cluster-center module id. If no matching module parameters are found, no energy correction is applied.

The normalized coordinates are calculated after the position correction block:

```text
dx = (hit.x - center.x) / center.module_size_x
dy = (hit.y - center.y) / center.module_size_y
```

For eight energy coefficients `c0` through `c7`, the correction is:

```text
denominator = 1
             + c1*dx^2
             + c2*dy^2
             + c3*dx^2*dy^2
             + c4*dx^4
             + c5*dy^4
             + c6*dx
             + c7*dy

E_corrected = E0 / c0 / denominator
E_Scorr     = E_corrected - E0
```

The hit is updated as follows:

```cpp
hit.E_Scorr = GetEneBias(parameters, dx, dy, hit.E);
hit.E += hit.E_Scorr;
```

`E_Scorr` stores the additive correction amount, not a multiplicative correction factor. After a successful correction, `kSEneCorr` is set.

## Correction Flags

The correction status is stored in `HyCalHit::flag`:

```cpp
kDenCorr   // position density correction has been applied
kSEneCorr  // S-shape energy correction has been applied
```

The two flags are independent. Position and energy correction can therefore be enabled or disabled separately, and each correction is protected against repeated application.

## Compact End-to-End Summary

```text
1. Load both position and energy files for Set_1GeV and Set_2GeV.
2. Select Set_1GeV or Set_2GeV from the run number.
3. Form a cluster and apply leakage correction.
4. Reconstruct the initial position and energy.
5. Determine an energy group from energy, angle, and resolution.
6. Determine a geometry group from the cluster-center module.
7. Apply the position bias to hit.x and hit.y.
8. If the event is identified as ep or ee, look up module-specific energy parameters.
9. Use the corrected x/y position to calculate and add the energy correction.
10. Set flags to prevent repeated correction.
```

## PRad2 Integration Notes

The PRad2 implementation uses the same PRad1 formulas and parameter files, but
the active correction set is selected from the beam energy resolved by
`LoadRunConfig()` from `database/runinfo/general.json`. The selection does not
use run number directly.

The switches live under `hycal` in `database/reconstruction_config.json`:

```json
"density_correction": false,
"s_shape_energy_correction": false
```

The available file mapping is also configured there. Current behavior is:

```text
beam energy < 1.0 GeV -> density_params/set_1GeV.dat + s_energy_params/ecorrect_1GeV.dat
2.0 GeV < beam energy < 3.0 GeV -> density_params/set_2GeV.dat + s_energy_params/ecorrect_2GeV.dat
beam energy > 3.0 GeV -> interface reserved; skip correction until files are provided
```

Missing files are non-fatal. The pipeline logs a warning and disables these two
corrections for that run.
