# prad2det

Static library for PRad-II detector reconstruction. Consumes decoded
front-end data from `prad2dec` and produces physics-level outputs
(per-module hits, strip clusters, 2-D GEM hits).

## Components

| Header | Role |
|---|---|
| `HyCalSystem` | HyCal geometry, DAQ map, calibration constants, sector-grid neighbour lookup, energy/position resolution model. |
| `HyCalCluster` | Island clustering with log-weighted center-of-gravity, profile-based shower split, and shower-depth correction on top of a `HyCalSystem`. See [`docs/technical_notes/hycal_clustering/hycal_clustering.md`](../docs/technical_notes/hycal_clustering/hycal_clustering.md). |
| `GemSystem` | GEM detector hierarchy and per-event chain: pedestal subtraction, sorting common-mode, zero-suppression, and APV → strip mapping. Full-readout vs online-zero-suppressed mode is auto-detected per APV. |
| `GemCluster` | Per-plane strip clustering (group + split + charge-weighted position) followed by X/Y matching (Cartesian-with-cuts or ADC-sorted). See [`docs/technical_notes/gem_clustering/gem_clustering.md`](../docs/technical_notes/gem_clustering/gem_clustering.md). |
| `GemPedestal` | Per-strip GEM pedestal accumulator; writes the APV-block text file `GemSystem::LoadPedestals` reads. |
| `DetectorTransform` | 3×3 rotation + translation for the detector → lab-frame transform; the rotation matrix is cached and invalidated through the `set(...)` mutators. |
| `RunInfoConfig` | Header-only loader for the run-period geometry / calibration JSON (`runinfo`). Picks the largest `run_number` ≤ the requested run. |
| `PipelineBuilder` | Fluent helper that wires up an entire reconstruction pipeline (HyCal + GEM + transforms + matching parameters) from the standard config files. |
| `EventData`, `EventData_io` | Shared per-event data structures and ROOT-tree branch helpers used by both the live monitor and the offline replay. |

## PipelineBuilder

Three callers — analysis scripts, the live server, and the Python
bindings — used to hand-wire the same `Init` → `LoadCalibration` →
`LoadPedestals` → `LoadCommonModeRange` → `SetReconConfigs` sequence.
Forgetting any step (most painfully the GEM crate remap derived from
`daq_cfg.roc_tags`) silently dropped data. `PipelineBuilder`
consolidates the wiring so it is impossible to omit:

```cpp
#include "PipelineBuilder.h"

auto p = prad2::PipelineBuilder()
    .set_run_number_from_evio(evio_path)
    .build();                              // throws on hard failures

fdec::HyCalCluster hc(p.hycal);
hc.SetConfig(p.hycal_cluster_cfg);
gem::GemCluster   gem_cl;
fdec::WaveAnalyzer wa(p.daq_cfg.wave_cfg);
// p.hycal, p.gem, p.hycal_transform, p.gem_transforms[…] are ready.
```

Boundary: the builder owns *detectors-ready* — initialised + calibrated
`HyCalSystem` + `GemSystem`, prepared `DetectorTransform`s, the
matching σ parameters, and the resolved input paths. It does **not**
own per-event scratch (`HyCalCluster`, `GemCluster`, `WaveAnalyzer`)
or server-only concerns (monitor config, trigger filters, EPICS,
livetime, histograms); those stay in the caller.

## Usage

```cpp
// GEM — per-event pipeline
gem::GemSystem gsys;
gsys.Init("database/gem_map.json");
gsys.LoadPedestals("gem_ped.txt");      // required only for full-readout data
gem::GemCluster gcl;

for each event {
    gsys.Clear();
    gsys.ProcessEvent(ssp_event_data);  // from prad2dec
    gsys.Reconstruct(gcl);
    for (const auto &h : gsys.GetAllHits()) {
        // h.x, h.y, h.det_id, h.x_charge, h.y_charge, ...
    }
}

// HyCal — feed per-module energies, cluster, read back
fdec::HyCalSystem hsys;
hsys.Init("database/hycal_map.json");
hsys.LoadCalibration("database/hycal_calib.json");
fdec::HyCalCluster hcl(hsys);

hcl.Clear();
for (const auto &[idx, E] : hits) hcl.AddHit(idx, E, time);
hcl.FormClusters();
std::vector<fdec::ClusterHit> out;
hcl.ReconstructHits(out);
```

## Dependencies

- [prad2dec](../prad2dec/README.md) — provides the SSP / FADC event-data
  types and the `DaqConfig` consumed by `PipelineBuilder`.
- [nlohmann/json](https://github.com/nlohmann/json) — fetched
  automatically; used for all JSON-backed configuration loaders.
