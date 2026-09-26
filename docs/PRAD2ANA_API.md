# prad2ana — C++ API Reference

`libprad2ana.a` — the offline replay and physics-analysis library.
Depends on [`prad2dec`](PRAD2DEC_API.md), [`prad2det`](PRAD2DET_API.md),
and ROOT 6.0+. All public symbols live in namespace `analysis`, except
the gain-factor helpers under `prad2`.

For a high-level walkthrough of the executables built on top of this
library (`prad2ana_replay_rawdata`, `prad2ana_replay_recon`,
`prad2ana_epCalib`, …) see [`analysis/README.md`](../analysis/README.md).
This document is the symbol reference for callers that link directly
against `libprad2ana.a` (analysis tools, ACLiC scripts in
`analysis/scripts/`, downstream user code).

| Header | Public symbols |
|---|---|
| [`Replay.h`](#replayh) | `analysis::Replay`, type aliases `EventVars`, `EventVars_Recon`, `FillPeaksFromWaveforms`, `ReconstructGemStrips` |
| [`ConfigSetup.h`](#configsetuph) | `analysis::RunConfig` (alias to `prad2::RunConfig`), `gRunConfig`, `LabTransforms`, `BuildLabTransforms`, `ApplyToLab`, `ApplyToLocal`, `ApplyToHyCal`, `get_run_str`, `get_run_int` |
| [`PhysicsTools.h`](#physicstoolsh) | `analysis::PhysicsTools`, `GEMHit`, `HCHit`, `DataPoint`, `MollerEvent`, `MollerData`, `kPbWO4Pitch`, `InHyCalRing` |
| [`MatchingTools.h`](#matchingtoolsh) | `analysis::MatchingTools`, `MatchHit`, `MatchHit_perChamber`, `MatchFlag`, `ProjectHit`, `GetProjection*`, `GetProjectionHits` |
| [`gain_factor.h`](#gain_factorh) | `prad2::GainFactor`, `GainFactorTable`, `GainCorrTable`, `FindGainFactorFile`, `LoadGainFactors`, `ComputeGainCorrection` |
| [`SlowControl.h`](#slowcontrolh) | `analysis::ScalerRow`, `EpicsRow`, `LoadScalerRows`, `LoadEpicsRows`, `SelectDscPair`, `SortByEvent`, `DeltaLivetime`, `ChargeSums` |
| [`ToolUtils.h`](#toolutilsh) | `analysis::ParseIntOption`, `ParseFloatOption`, `InvalidOption*`, `Is*Name`, `ExpandInputPath`, `CollectInputs`, `RunCommand`, `InitRootThreading`, `TreeEntries`, `DistributeEventBudget`, `RunFilesInRounds`, `ParallelFor`, `RunReplayPool`, `HistList`, `Book`, `AddAll`, `MergeTopLevelHistograms` |

---

## `Replay.h`

`analysis::Replay` — converts raw DAQ data (EVIO) into ROOT trees.
Three replay entry points:

1. **Raw replay** (`Process`) — per-channel waveform/peak data,
   `events` tree (`prad2::RawEventData`).
2. **Recon replay** (`ProcessWithRecon`) — full reconstruction
   (HyCal clusters, GEM hits, HyCal↔GEM matches), `recon` tree
   (`prad2::ReconEventData`).
3. **X17 recon replay** (`ProcessWithReconX17`) — cluster-trigger
   reconstruction with the X17 blind-sample selection, also written to
   a `recon` tree.

All replay entry points also write the side trees `scalers` (DSC2) and `epics`
(0x001F text banks) — see [`prad2det`](PRAD2DET_API.md#eventdata_ioh).

### Type aliases

```cpp
using EventVars       = prad2::RawEventData;
using EventVars_Recon = prad2::ReconEventData;
```

### Construction and configuration

```cpp
Replay r;
r.LoadDaqConfig(json_path);          // delegates to evc::load_daq_config
r.LoadHyCalMap(hycal_map_json_path); // populates the internal fdec::HyCalSystem used by moduleName/moduleType/moduleID
```

Calling `LoadHyCalMap` is strongly recommended — without it every
channel is unknown (`MOD_UNKNOWN`, `module_id` `-1`) and gets dropped.

### Channel introspection

```cpp
std::string moduleName(int roc, int slot, int ch) const;
prad2::ModuleType moduleType(int roc, int slot, int ch) const;
int moduleID(int roc, int slot, int ch) const;
```

`moduleID` returns the globally-unique encoding documented in
[`EventData.h`](PRAD2DET_API.md#eventdatah) (G: 1..1156, W: 1001..2152,
Veto: 3001..3004, LMS: 3100..3103). Returns `-1` for unknown channels.

### Replay drivers

```cpp
bool Process(const std::string &input_evio,
             const std::string &output_root,
             RunConfig &gRunConfig,
             const std::string &db_dir,
             int max_events = -1,
             bool write_peaks = false,
             const std::string &daq_config_file = "");

bool ProcessWithRecon(const std::string &input_evio,
                      const std::string &output_root,
                      RunConfig &gRunConfig,
                      const std::string &db_dir,
                      const std::string &daq_config_file = "",
                      const std::string &gem_ped_file    = "",
                      float zerosup_override            = 0.f,
                      bool prad1                        = false);

bool ProcessWithReconX17(const std::string &input_evio,
                         const std::string &output_root,
                         RunConfig &gRunConfig,
                         const std::string &db_dir,
                         const std::string &daq_config_file = "",
                         const std::string &gem_ped_file    = "",
                         float zerosup_override            = 0.f);
```

`max_events <= 0` ⇒ process every event. `write_peaks=true` adds the
optional soft-analyzer + firmware-mode peak branches to the `events`
tree (`hycal.peak_*`, `hycal.daq_peak_*`).

`ProcessWithRecon` runs the full pipeline via
[`prad2::PipelineBuilder`](PRAD2DET_API.md#pipelinebuilderh):
HyCal clustering, GEM strip clustering, HyCal↔GEM matching with the
runinfo σ parameters; `gem_ped_file` overrides the per-run pedestal,
`zerosup_override > 0` overrides the GEM ZS threshold, `prad1=true`
selects the legacy ADC1881M readout path.

`ProcessWithReconX17` uses the same detector pipeline and output data
type, but accepts events carrying `TBIT_1cl`, `TBIT_2cl`, `TBIT_3cl`,
`TBIT_lms`, or `TBIT_alpha`. It keeps every 1-/2-cluster event and only
the deterministic 10% of 3-cluster events for which
`event_num % 10 == 8`.

### Re-processing raw replay trees

```cpp
void FillPeaksFromWaveforms(prad2::RawEventData &ev, const fdec::HyCalSystem &hycal,
                            const fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

void ReconstructGemStrips(const prad2::RawEventData &ev, const gem::GemSystem &gem_sys,
                          gem::GemCluster &clusterer, std::vector<gem::GEMHit> &hits,
                          std::vector<std::array<std::vector<gem::StripCluster>, 2>>
                              *plane_clusters = nullptr);
```

`FillPeaksFromWaveforms` re-derives `npeaks` and `peak_height/time/integral`
of every PbWO4 channel from its stored samples, for raw trees written
without the peak branches (no module time offset is applied).
`ReconstructGemStrips` re-runs GEM clustering and X/Y matching, with
`gem_sys`'s per-detector configs, on the strip hits stored in a raw tree
(pedestal, common mode and zero suppression already applied): `hits`
receives the 2D hits of all detectors in detector order, `plane_clusters`
(when given) the kept clusters as `[det][0 = X, 1 = Y]`.  Strips of an
unknown detector or plane are skipped.

---

## `ConfigSetup.h`

Analysis-side helpers around `prad2::RunConfig`. Header-only.

### Re-exports

```cpp
using RunConfig = ::prad2::RunConfig;
using ::prad2::LoadRunConfig;
using ::prad2::WriteRunConfig;

inline RunConfig gRunConfig;   // single-run global, multi-run code should use locals
```

See the [`prad2det`](PRAD2DET_API.md#runinfoconfigh) reference for
`RunConfig` fields and `LoadRunConfig` selection rules.

### Lab-frame transforms

```cpp
struct LabTransforms {
    DetectorTransform                hycal;
    std::array<DetectorTransform, 4> gem;
};

LabTransforms BuildLabTransforms(const RunConfig &geo = gRunConfig);

template <typename Hit>
void ApplyToLab(const DetectorTransform &xform, Hit &h);   // in-place lab = R*[h.x,h.y,h.z] + t
template <typename Hit>
void ApplyToLocal(const DetectorTransform &xform, Hit &h); // in-place inverse (labToLocal)
template <typename Hit>
void ApplyToHyCal(Hit &h, const RunConfig &geo = gRunConfig);   // h.x += target_x, h.y += target_y
```

`BuildLabTransforms` populates each `DetectorTransform` via `set(...)`,
so the rotation matrices are precomputed before `toLab` is called per
hit. Build once per run and reuse. `ApplyToHyCal` moves a hit into the
HyCal coordinate system; project it to the HyCal surface first.

### Run-number filename parsers

```cpp
std::string get_run_str(const std::string &file_name);   // "unknown" on failure
int         get_run_int(const std::string &file_name);   // -1 on failure
```

Both use [`prad2::run_number_from_path`](PRAD2DEC_API.md#eviofilesh):
the first `prad_<digits>` or `run_<digits>` in the file name
(case-insensitive, directory part ignored).

---

## `PhysicsTools.h`

`analysis::PhysicsTools` — owns the per-module energy histograms,
Moller geometry, and the kinematic calculations used by the calibration
and matching tools.

### Per-event types

`GEMHit { x, y, z, det_id }` — `det_id ∈ 0..3` for GEM1..GEM4 (5 = unset).

`HCHit { x, y, z, energy, center_id, flag }`.

`DataPoint { x, y, z, E }` — used inside `MollerEvent`.

```cpp
typedef std::pair<DataPoint, DataPoint> MollerEvent;
typedef std::vector<MollerEvent>        MollerData;
```

### Fiducial ring (free functions)

```cpp
inline constexpr double kPbWO4Pitch = 20.75;   // mm, nominal PbWO4 module size
bool InHyCalRing(double x, double y, double inner, double outer,
                 double pitch = kPbWO4Pitch);
```

`InHyCalRing` is a square-annulus cut: outside the inner square of
half-width `inner*pitch` and inside the outer square of half-width
`outer*pitch`. It is frame-agnostic; the caller picks lab, HyCal or
module-centre coordinates.

### Construction

```cpp
explicit PhysicsTools(fdec::HyCalSystem &hycal);
```

The reference must outlive the `PhysicsTools` instance. Histograms are
constructed in the body of the constructor (one `TH1F` per HyCal module
plus the 2-D / Moller histograms).

### Per-module energy histograms

```cpp
void  FillModuleEnergy(int module_id, float energy);
TH1F *GetModuleEnergyHist(int module_id) const;

void  FillEnergyVsModule(int module_id, float energy);
TH2F *GetEnergyVsModuleHist() const;

void  FillEnergyVsTheta(float theta_deg, float energy);
TH2F *GetEnergyVsThetaHist() const;

void  FillNeventsModuleMap();    // populate from filled per-module hists
TH2F *GetNeventsModuleMapHist() const;

static TH2Poly *MakeModuleMap(const fdec::HyCalSystem &hycal, const char *name,
                              const char *title, double half_range,
                              std::vector<int> &bin_by_index);
```

`MakeModuleMap` builds a `TH2Poly` with one rectangular bin per PbWO4
module and axes ±`half_range` (mm), created with `new` in the current
directory (caller owns it). `bin_by_index` is resized to
`module_count()` and holds the `TH2Poly` bin of each PbWO4 module by
module index, `-1` for every other module.

### Moller-event histograms

```cpp
void  FillMollerPhiDiff(float phi_diff);
void  FillMollerXY    (float x, float y);
void  FillMollerZ     (float z);
TH1F *GetMollerPhiDiffHist() const;
TH1F *GetMollerXHist()       const;
TH1F *GetMollerYHist()       const;
TH1F *GetMollerZHist()       const;
```

### Resolution / peak fits

```cpp
std::array<float, 3> FitPeakResolution(int module_id) const;   // {peak, sigma, chi2}

static std::array<double, 5> fitGaus(TH1F *h, float expectPeak = 0.f,
                                     bool withError = false);
static std::array<double, 5> fitCrystalBall(TH1F *h, float expectPeak = 0.f,
                                            float alpha = 1.5f, float n = 5.0f,
                                            bool withError = false);
static std::array<double, 5> fitPeak(TH1F *h, float expectPeak = 0.f, bool withError = false,
                                     bool useCrystalBall = false,
                                     float alpha = 0.5f, float n = 5.0f);
```

The static fits return `{mean, sigma, chi2/ndf, mean_error, sigma_error}`;
the two errors are zero unless `withError` is true.

### Kinematics (static)

```cpp
static constexpr float kProtonMass   = 938.272f;     // MeV
static constexpr float kElectronMass = 0.51099895f;  // MeV

static float ExpectedEnergy(float theta_deg, float Ebeam,
                            const std::string &type);   // "ep" or "ee"
static float EnergyLoss   (float theta_deg, float E);   // target + windows
static bool  HitP4(float x, float y, float z, float E, float m, TLorentzVector &p4);
static bool  isMoller_kinematic(float theta_deg1, float energy1,
                                float theta_deg2, float energy2,
                                float EBeam, float resolution);
```

`HitP4` builds the four-momentum of a particle of mass `m` and energy `E`
(MeV) emitted from the target (origin) towards the hit at `(x, y, z)`; it
returns `false` and zeroes `p4` when `E < m` or the hit is at the origin.
`isMoller_kinematic` is the elastic e-e check for Moller selection: the
energy sum within 5σ of `EBeam` and each energy within 3.5σ of its
expected value, σ = `resolution` · E / √(E in GeV).

### Moller geometry (static)

```cpp
static std::array<float, 2> GetMollerCenter(const MollerEvent &e1, const MollerEvent &e2);
static float GetMollerZdistance(const MollerEvent &e, float Ebeam);
static float GetMollerPhiDiff  (const MollerEvent &e1);   // ≈ 180° for elastic ee
static bool  isBackToBack      (const MollerEvent &e, float max_dev_deg);
static float GetPhiAngle       (float x, float y);
static float GetThetaAngle     (float x, float y, float z);   // degrees, seen from the target
```

`isBackToBack` is `|GetMollerPhiDiff(e)| < max_dev_deg`.

---

## `MatchingTools.h`

`analysis::MatchingTools` — HyCal cluster ↔ GEM hit matching with the
projection-plus-cut algorithm ported from
`PRadAnalyzer/PRadDetMatch.cpp`.

### Matching flag enum

```cpp
enum MatchFlag : uint32_t {
    kGEM1Match = 0,  // bit 0
    kGEM2Match = 1,
    kGEM3Match = 2,
    kGEM4Match = 3,
};
```

### Projection helpers

```cpp
struct ProjectHit { float x_proj, y_proj, z_proj; };

ProjectHit GetProjectionHits(float x, float y, float z, float projection_z);
void GetProjection(HCHit &hc,                       float projection_z);
void GetProjection(std::vector<HCHit> &hc,          float projection_z);
void GetProjection(GEMHit &gem,                     float projection_z);
void GetProjection(std::vector<GEMHit> &gem,        float projection_z);
```

`GetProjection` updates `(x, y, z)` of each hit to project a straight
line from the target through the hit to the requested `projection_z`.

### Match outputs

`MatchHit` — one HyCal cluster paired with all candidate GEM hits per
detector, plus the chosen "best" pair:

```cpp
class MatchHit {
public:
    HCHit                hycal_hit;
    std::vector<GEMHit>  gem1_hits, gem2_hits, gem3_hits, gem4_hits;

    GEMHit               gem[2];      // best-matched upstream and downstream
    uint32_t             mflag = 0;   // OR of MatchFlag bits
    uint16_t             hycal_idx = 0;

    MatchHit(const HCHit &, std::vector<GEMHit> &g1, std::vector<GEMHit> &g2,
             const std::vector<GEMHit> &g3, const std::vector<GEMHit> &g4);
};
```

`MatchHit_perChamber` — the per-chamber variant; stores the best match
per detector as a flat `[det_id][x/y/z]` array for analyses that don't
collapse to one upstream/downstream pair:

```cpp
class MatchHit_perChamber {
public:
    HCHit       hycal_hit;
    float       gem_hits[4][3] = {};
    uint32_t    mflag = 0;
    uint16_t    hycal_idx = 0;

    explicit MatchHit_perChamber(const HCHit &);
};
```

### `MatchingTools` methods

```cpp
explicit MatchingTools(int postMatchMethod = 1);

std::vector<MatchHit> Match(
    std::vector<HCHit> &hycalHits,
    const std::vector<GEMHit> &gem1, const std::vector<GEMHit> &gem2,
    const std::vector<GEMHit> &gem3, const std::vector<GEMHit> &gem4) const;

std::vector<MatchHit_perChamber> MatchPerChamber(
    std::vector<HCHit> &hycalHits,
    const std::vector<GEMHit> &gem1, const std::vector<GEMHit> &gem2,
    const std::vector<GEMHit> &gem3, const std::vector<GEMHit> &gem4) const;

void SetMatchRange     (float range);   // mm; default 15
void SetSquareSelection(bool sq);       // true = square cut, false = circular
```

`postMatchMethod` is typically supplied from
`prad2::Pipeline::match_method` (loaded by `PipelineBuilder` from
`reconstruction_config.json:matching.match_method`).
Both modes are branches of `PostMatch`: `1` keeps, per upstream/downstream
GEM pair, the hit closest to the HyCal cluster; any other value picks the
upstream/downstream GEM-hit pair with minimum ΔR between the two hits.

---

## `gain_factor.h`

Header-only loaders for the per-module LMS gain-factor database. Lives
under namespace `prad2` (not `analysis`) so other libraries can share it
without dragging ROOT.

### File format

`<dir>/prad_XXXXXX_LMS.dat` — whitespace-delimited:

```
Name  lms_peak  lms_sigma  lms_chi2/ndf  g1  g2  g3
```

Only `W*` and `G*` lines are read; LMS header rows and other prefixes
are silently skipped.

### Selection rule

Same as `LoadRunConfig`: `run_num >= 0` ⇒ largest run ≤ requested;
`run_num < 0` ⇒ latest. If no file satisfies `run ≤ run_num`,
`FindGainFactorFile` falls back to the nearest available file with a
warning.

### Types

```cpp
struct GainFactor { float g[3] = {0, 0, 0}; };   // g1, g2, g3

struct GainFactorTable {
    static constexpr int MAX_W = 1157;   // W1..W1156
    static constexpr int MAX_G = 901;    // G1..G900
    GainFactor w[MAX_W];
    GainFactor g[MAX_G];
    int  run_number = -1;
    bool loaded     = false;
};

struct GainCorrTable {
    struct Entry {
        float corr[3] = {1.f, 1.f, 1.f};   // ref.g[j] / cur.g[j]
        float avg     = 1.f;               // mean over non-zero corr[]
    };
    Entry w[GainFactorTable::MAX_W];
    Entry g[GainFactorTable::MAX_G];
    int   ref_run = -1;
    int   cur_run = -1;
};
```

### Functions

```cpp
std::string     prad2::FindGainFactorFile(const std::string &dir, int run_num);
GainFactorTable prad2::LoadGainFactors   (const std::string &dir, int run_num);

GainCorrTable prad2::ComputeGainCorrection(const GainFactorTable &ref_tbl,
                                           const GainFactorTable &cur_tbl);
GainCorrTable prad2::ComputeGainCorrection(const std::string &dir,
                                           int cur_run, int ref_run);
```

`new_adc2mev = old_adc2mev * corr.w[id].avg` is the typical applier
(see top of header for the worked example).

---

## `SlowControl.h`

Slow-event rows of replayed ROOT files and the live-charge integration
shared by `prad2ana_replay_filter` and `prad2ana_live_charge`.

`ScalerRow` is a `prad2::RawScalerData` plus `bool good`;
`EpicsRow { event_number, ti_ticks, unix_time, sync_counter, run_number,
good, updates }` holds one `epics` row, `updates` mapping channel name to
the readings the row carries.

```cpp
bool LoadScalerRows(const std::vector<std::string> &files, std::vector<ScalerRow> &out,
                    const char *tag, bool *has_good = nullptr);
bool LoadEpicsRows (const std::vector<std::string> &files, std::vector<EpicsRow> &out,
                    const char *tag, bool *has_good = nullptr);

std::pair<uint32_t, uint32_t> SelectDscPair(const prad2::RawScalerData &row,
                                            const std::string &source, int channel);
template <class Row> std::vector<size_t> SortByEvent(const std::vector<Row> &rows);
std::vector<double> DeltaLivetime(const std::vector<ScalerRow> &rows,
                                  const std::vector<size_t> &order,
                                  const std::string &source, int channel,
                                  double scale = 1.0);
```

The loaders append every file's `scalers` / `epics` rows in file order
(`good` is replay_filter's per-row verdict, `true` without that branch)
and return `false` when a file cannot be opened.  `DeltaLivetime` gives
each row's slice-local livetime, Δgated / Δungated since the previous row
in event order × `scale` (`-1` where undefined), using
[`dsc::delta_live_ratio`](PRAD2DEC_API.md#dscdatah) with the first row's
predecessor at (0, 0).

`ChargeSums` accumulates Σ lf_b · Δt · ½(I_a + I_b) over adjacent
checkpoints with `AddPair(ticks_a, ticks_b, live_fraction_b, current_a,
current_b, good)`: gated sums (`value_nC`, `live_seconds`,
`real_seconds`, pair counters) for `good` pairs, `ungated_*` sums for all
pairs with valid data; `operator+=` merges two sums.

---

## `ToolUtils.h`

Helpers shared by the command-line tools in `analysis/tools/`.

- **Options** — `ParseIntOption` / `ParseFloatOption` (the whole argument
  must be a number in range; the value is untouched on failure),
  `InvalidOptionValue(flag, value)` and `InvalidOption(argv, optind,
  optopt)` (print the getopt error, return exit code 2).
- **Inputs** — name filters `IsEvioName`, `IsRootName`, `IsRawRootName`,
  `IsReconRootName`, `IsLmsRootName`; `ExpandInputPath(path, keep)`
  expands a directory to its regular files passing `keep(name)`, sorted;
  `CollectInputs(argc, argv, first, keep, max_files = -1)` does this for
  every argument in order.
- **Subprocesses** — `RunCommand(argv)` runs a program without a shell
  and returns its exit status (127 when it could not be started).
- **Threads** — `InitRootThreading()` (call first in `main`),
  `TreeEntries(path, tree)`, `DistributeEventBudget(files, tree,
  max_events)` (per-file entry limits), `RunFilesInRounds(files,
  nthreads, job, after_round)`, `ParallelFor(n, n_threads, fn)` and
  `RunReplayPool(inputs, n_threads, daq_config, daq_map, output_for, run,
  ok)` (one `Replay` per thread).
- **Histograms** — `HistList`, `Book<H>(reg, args...)` (detached
  histogram registered in booking order), `AddAll(dst, src)` (pairwise
  merge of bundles booked the same way) and
  `MergeTopLevelHistograms(inputs, output)` (sum by name of the top-level
  histograms of several files).

---

## Build / link

```cmake
target_link_libraries(your_tool PRIVATE prad2ana)
# transitively pulls in: prad2dec, prad2det, nlohmann_json,
#                        ROOT::{Core, Tree, RIO, Hist, Graf, Gpad, Spectrum}
```

`libprad2ana.a` is built with PIC so it can also be linked into ACLiC-
built shared objects (the `analysis/scripts/*.C` macros in installed
mode).

## Dependencies

- [`prad2dec`](PRAD2DEC_API.md) — EVIO reader, decoders, soft + firmware
  waveform analyzers.
- [`prad2det`](PRAD2DET_API.md) — `HyCalSystem`, `GemSystem`,
  `PipelineBuilder`, `RunConfig`, replay tree schema.
- ROOT 6.0+ — `TFile`, `TTree`, `TH1F`, `TH2F`, `TF1`, `TSpectrum`.
- [nlohmann/json](https://github.com/nlohmann/json) — fetched
  automatically by the top-level CMake.
