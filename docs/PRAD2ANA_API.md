# prad2ana — C++ API Reference

`libprad2ana.a` — the offline replay and physics-analysis library.
Depends on [`prad2dec`](PRAD2DEC_API.md), [`prad2det`](PRAD2DET_API.md),
and ROOT 6.0+. All public symbols live in namespace `analysis`, except
the gain-factor helpers under `prad2`.

For a high-level walkthrough of the executables built on top of this
library (`prad2ana_replay_rawdata`, `prad2ana_replay_recon`,
`prad2ana_replay_gainCorr`, …) see [`analysis/README.md`](../analysis/README.md).
This document is the symbol reference for callers that link directly
against `libprad2ana.a` (analysis tools, ACLiC scripts in
`analysis/scripts/`, downstream user code).

| Header | Public symbols |
|---|---|
| [`Replay.h`](#replayh) | `analysis::Replay`, type aliases `EventVars`, `EventVars_Recon`, `LMSEventVars`, `FillPeaksFromWaveforms`, `ReconstructGemStrips` |
| [`ConfigSetup.h`](#configsetuph) | `analysis::RunConfig` (alias to `prad2::RunConfig`), `gRunConfig`, `LabTransforms`, `BuildLabTransforms`, `ApplyToLab`, `ApplyToLocal`, `ApplyToHyCal`, `get_run_str`, `get_run_int` |
| [`PhysicsTools.h`](#physicstoolsh) | `analysis::PhysicsTools`, `GEMHit`, `HCHit`, `DataPoint`, `MollerEvent`, `MollerData`, `kPbWO4PitchX`, `kPbWO4PitchY`, `InHyCalRing` |
| [`MatchingTools.h`](#matchingtoolsh) | `analysis::MatchingTools`, `MatchHit`, `MatchHit_perChamber`, `MatchFlag`, `ProjectHit`, `GetProjection*`, `GetProjectionHits`, `ClusterToLab`, `GemHitToLab`, `FillReconMatches`, `MatchReconEvent` |
| [`gain_factor.h`](#gain_factorh) | `prad2::RefGainDir`, `GainCorrDir`, `RefGainFactor`, `RefGainTable`, `FindRunFile`, `FindRefGainFile`, `LoadRefGainFile`, `LoadRefGain`, `GainCorrTable`, `FindGainCorrRootFile`, `GainCorrTimeSeries`, `LoadGainCorrTimeSeries`, `FitResult`, `gain_hist_fitter` |
| [`GainCorrCompute.h`](#gaincorrcomputeh) | `analysis::GainBatch`, `GainPlotConfig`, `GainPlotStore`, `kGain*` constants, `MakeGainHist`, `MakeLMSOutputFile`, `GainCorrOutputPath`, `ReplayLMSFiles`, `EnsureGainCorr`, `SetupGainBranches`, `FlushGainBatch`, `ComputeGainCorrections` |
| [`SlowControl.h`](#slowcontrolh) | `analysis::ScalerRow`, `EpicsRow`, `LoadScalerRows`, `LoadEpicsRows`, `SelectDscPair`, `SortByEvent`, `DeltaLivetime`, `ChargeSums` |
| [`ToolUtils.h`](#toolutilsh) | `analysis::ParseIntOption`, `ParseFloatOption`, `InvalidOption*`, `Is*Name`, `ExpandInputPath`, `CollectInputs`, `RunCommand`, `InitRootThreading`, `TreeEntries`, `DistributeEventBudget`, `RunFilesInRounds`, `ParallelFor`, `RunReplayPool`, `HistList`, `Book`, `AddAll`, `MergeTopLevelHistograms` |

---

## `Replay.h`

`analysis::Replay` — converts raw DAQ data (EVIO) into ROOT trees.
Four replay entry points:

1. **Raw replay** (`Process`) — per-channel waveform/peak data,
   `events` tree (`prad2::RawEventData`).
2. **Recon replay** (`ProcessWithRecon`) — full reconstruction
   (HyCal clusters, GEM hits, HyCal↔GEM matches), `recon` tree
   (`prad2::ReconEventData`). The X17 modes are flags of this call.
3. **Raw → recon replay** (`ProcessRaw2Recon`) — the same reconstruction
   starting from a raw replay ROOT file (`events` tree) instead of EVIO.
4. **LMS gain replay** (`Process_LMSgainFactor`) — LMS and alpha events
   only, `lms_gain` tree (`prad2::LMSEventData`), input of the gain
   correction (see [`GainCorrCompute.h`](#gaincorrcomputeh)).

All replay entry points also write the side trees `scalers` (DSC2) and `epics`
(0x001F text banks) — see [`prad2det`](PRAD2DET_API.md#eventdata_ioh);
`ProcessRaw2Recon` copies them from its input file.

### Type aliases

```cpp
using EventVars       = prad2::RawEventData;
using EventVars_Recon = prad2::ReconEventData;
using LMSEventVars    = prad2::LMSEventData;
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
std::tuple<int, int, int> moduleLocation(int module_id) const;   // (crate, slot, ch)
```

`moduleID` returns the globally-unique encoding documented in
[`EventData.h`](PRAD2DET_API.md#eventdatah) (G: 1..900, W: 1001..2156,
Veto: 3001..3004, LMS: 3100..3103). Returns `-1` for unknown channels;
`moduleLocation` is the reverse lookup and returns `{-1, -1, -1}` for an
unknown `module_id`.

### Replay drivers

```cpp
bool Process(const std::string &input_evio,
             const std::string &output_root,
             RunConfig &gRunConfig,
             const std::string &db_dir,
             const std::string &recon_config_file,
             int max_events = -1,
             bool write_peaks = false,
             const std::string &daq_config_file = "",
             const float zerosup_override = 5.f,
             bool Ecalib = false,
             bool noWaveform = false);

bool ProcessWithRecon(const std::string &input_evio,
                      const std::string &output_root,
                      RunConfig &gRunConfig,
                      const std::string &db_dir,
                      const std::string &recon_config_file,
                      const std::string &daq_config_file = "",
                      const std::string &gem_ped_file    = "",
                      float zerosup_override            = 5.f,
                      bool prad1 = false, bool x17 = false, bool x17_blind = false,
                      bool random = false, bool gem_hit = false,
                      bool shower_depth_corr = true);

bool ProcessRaw2Recon(const std::string &input_raw,
                      const std::string &output_root,
                      RunConfig &gRunConfig,
                      const std::string &db_dir,
                      const std::string &recon_config_file,
                      const std::string &daq_config_file = "",
                      const std::string &gem_ped_file    = "",
                      bool x17 = false, bool x17_blind = false,
                      bool random = false, bool gem_hit = false,
                      bool shower_depth_corr = true);

bool Process_LMSgainFactor(const std::string &input_evio,
                           const std::string &output_root,
                           const std::string &db_dir,
                           const std::string &daq_config_file);
```

`recon_config_file` is the `reconstruction_config.json` handed to
[`prad2::PipelineBuilder`](PRAD2DET_API.md#pipelinebuilderh).
`max_events <= 0` ⇒ process every event. `write_peaks=true` adds the
optional soft-analyzer + firmware-mode peak branches to the `events`
tree (`hycal.peak_*`, `hycal.daq_peak_*`). `zerosup_override >= 0`
overrides the GEM ZS threshold (σ). `Ecalib=true` keeps only
`TBIT_sum` events, skips the GEM data and leaves out the `event_num`,
`timestamp` and waveform-sample branches; `noWaveform=true` leaves out
only the waveform-sample branches.

`ProcessWithRecon` runs the full pipeline via
`PipelineBuilder`: HyCal clustering, GEM strip clustering, HyCal↔GEM
matching with the runinfo σ parameters; `gem_ped_file` overrides the
per-run pedestal, `prad1=true` selects the legacy ADC1881M readout path.
Events are kept when they carry `TBIT_sum`, `TBIT_lms`, `TBIT_alpha` or
`TBIT_3cl` (`random=true` keeps every event as a sum trigger); with
`x17 && x17_blind` a `TBIT_3cl` event is only kept when
`event_num % 10 == 8`. `x17=true` also leaves the veto/LMS summary
branches out of the `recon` tree, `gem_hit=true` adds the GEM hit
branches (with their quality variables) and the `gem_cl_*` block, and `shower_depth_corr=false` turns off the shower-depth
correction of the cluster positions.

`ProcessRaw2Recon` applies the same selection and reconstruction to the
`events` tree of a raw replay file, re-running GEM clustering on the
stored strip hits (`ReconstructGemStrips`).

All three HyCal replays apply the per-event gain correction from
`prad2::LoadGainCorrTimeSeries` (see [`gain_factor.h`](#gain_factorh)).

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
inline constexpr double kPbWO4PitchX = 20.77;   // mm, nominal PbWO4 module size
inline constexpr double kPbWO4PitchY = 20.75;
bool InHyCalRing(double x, double y, double inner, double outer,
                 double pitch_x = kPbWO4PitchX, double pitch_y = kPbWO4PitchY);
```

`InHyCalRing` is a square-annulus cut: outside the inner rectangle of
half-widths `inner*pitch_x` × `inner*pitch_y` and inside the outer
rectangle of half-widths `outer*pitch_x` × `outer*pitch_y`. It is frame-agnostic; the caller picks lab, HyCal or
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
TH1F *GetMollerPhiDiffHist() const;
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

```cpp
HCHit  ClusterToLab(const DetectorTransform &hycal_xform, const fdec::ClusterHit &cluster,
                    bool at_shower_depth = true);
GEMHit GemHitToLab (const std::array<DetectorTransform, 4> &gem_xforms, const gem::GEMHit &hit);
```

`ClusterToLab` places the cluster at its shower-max depth behind the
HyCal face (on the face, z = 0, when `at_shower_depth` is false) before
the transform; `GemHitToLab` uses `gem_xforms[hit.det_id]` (a hit with
`det_id` outside 0..3 keeps its plane coordinates).

### Match outputs

`MatchHit` — one HyCal cluster paired with all candidate GEM hits per
detector, plus the chosen "best" pair:

```cpp
class MatchHit {
public:
    HCHit                hycal_hit;
    std::vector<GEMHit>  gem1_hits, gem2_hits, gem3_hits, gem4_hits;

    GEMHit               gem[2];      // [0] downstream GEM1/GEM2, [1] upstream GEM3/GEM4
    uint32_t             mflag = 0;   // OR of MatchFlag bits
    uint8_t              hycal_idx = 0;

    MatchHit(const HCHit &, const std::vector<GEMHit> &g1, const std::vector<GEMHit> &g2,
             const std::vector<GEMHit> &g3, const std::vector<GEMHit> &g4);
};
```

`MatchHit_perChamber` — the per-chamber variant; keeps, per GEM chamber,
the hits inside the matching window (closest first) for analyses that
don't collapse to one upstream/downstream pair:

```cpp
class MatchHit_perChamber {
public:
    HCHit               hycal_hit;
    std::vector<GEMHit> gem_hits[4];   // [0..3] = GEM1..GEM4
    uint32_t            mflag = 0;
    uint8_t             hycal_idx = 0;

    MatchHit_perChamber(const HCHit &);
};
```

### `MatchingTools` methods

```cpp
explicit MatchingTools(int postMatchMethod = 1);

std::vector<MatchHit> Match(
    const std::vector<HCHit> &hycalHits,
    const std::vector<GEMHit> &gem1, const std::vector<GEMHit> &gem2,
    const std::vector<GEMHit> &gem3, const std::vector<GEMHit> &gem4) const;

std::vector<MatchHit_perChamber> MatchPerChamber(
    const std::vector<HCHit> &hycalHits,
    const std::vector<GEMHit> &gem1, const std::vector<GEMHit> &gem2,
    const std::vector<GEMHit> &gem3, const std::vector<GEMHit> &gem4) const;

void SetMatchRange     (float range);   // mm; default 15
void SetSquareSelection(bool sq);       // true = square cut, false = circular
void SetEnergyDependent(bool ed);       // default false
void SetMatchSigma     (float sigma);   // mm; default 1, used when energy-dependent
void Configure(const prad2::RunConfig &cfg);   // all four from cfg.matching_*
```

```cpp
void FillReconMatches(prad2::ReconEventData &ev, const std::vector<MatchHit> &matched,
                      const std::vector<MatchHit_perChamber> &matched_per_chamber);
void MatchReconEvent (prad2::ReconEventData &ev, const MatchingTools &matching);
```

`FillReconMatches` writes the match results into the recon-tree matching
fields (clusters indexed by `hycal_idx`); `MatchReconEvent` matches the
event's lab-frame `cl_*` clusters to its GEM hits (`det_id` 0..3) and
fills them with `FillReconMatches`.

`postMatchMethod` is typically supplied from
`prad2::Pipeline::match_method` (loaded by `PipelineBuilder` from
`reconstruction_config.json:matching.match_method`).
Both modes are branches of `PostMatch`: `1` keeps, per upstream/downstream
GEM pair, the hit closest to the HyCal cluster; any other value picks the
upstream/downstream GEM-hit pair with minimum ΔR between the two hits.

---

## `gain_factor.h`

Header-only loaders for the per-module LMS reference gains and the
time-dependent gain correction. Lives under namespace `prad2` (not
`analysis`); it uses ROOT (`TFile`/`TTree`/`TH1F`/`TF1`) for the
correction files and the peak fitter.

### Files

| Kind | Directory | File |
|---|---|---|
| Reference gain (from `prad2ana_refGain_produce`) | `RefGainDir(db)` = `<db>/gain_factor/ref_gain/` | `prad_XXXXXX_LMS.dat` |
| Gain correction (from `prad2ana_replay_gainCorr`) | `GainCorrDir(db)` = `<db>/gain_factor/gain_correction/` | `prad_XXXXXX_gain_corr.root`, TTree `gain_corr`, one entry per LMS batch |

The `.dat` file is whitespace-delimited after one header line:

```
Name  lms_peak  lms_sigma  lms_chi2/ndf  g1  g2  g3
```

Only `W*` and `G*` lines are read; LMS rows and other prefixes are
skipped. The `gain_corr` files carry PbWO4 (W) modules only. Writers
always use `GainCorrDir(<db>)`; readers use `RunConfig::gain_data_dir`,
which defaults to the same directory.

### Selection rule

Both kinds of file are looked up by **exact** run number:
`FindRunFile(dir, run, suffix)` returns `dir/prad_<run as %06d><suffix>`
when that regular file exists, otherwise an empty string (no fallback to
another run; `run` outside 0..999999 never matches).

### Types

```cpp
struct RefGainFactor { float g[3] = {0, 0, 0}; };   // g1, g2, g3

struct RefGainTable {
    static constexpr int MAX_W = 1157;   // W1..W1156
    static constexpr int MAX_G = 901;    // G1..G900
    RefGainFactor w[MAX_W];
    RefGainFactor g[MAX_G];
    int  run_number = -1;
    bool loaded     = false;
};

struct GainCorrTable {
    struct Entry {
        float corr[3] = {1.f, 1.f, 1.f};   // ref.g[j] / cur.g[j]
        float avg     = 1.f;               // mean of the three, a failed fit counting as 1
    };
    Entry w[MAX_W];
    Entry g[MAX_G];
    int   ref_run = -1;
    int   cur_run = -1;
    float ModuleGain(int id, int ref_pmt = 23) const;
};

struct GainCorrTimeSeries {
    struct Batch {
        int           event_num_start = 0;
        int           event_num_end   = 0;
        uint32_t      unix_time       = 0;   // batch midpoint, 0 if unknown
        GainCorrTable corr;
    };
    std::vector<Batch> batches;   // ascending event_num_start
    int  run_num = -1;
    bool loaded  = false;
    const GainCorrTable &GetCorr(int event_num) const noexcept;
};

struct FitResult { float mean = 0, sigma = 0, chi2pndf = 0; };
```

### Functions

```cpp
std::string  prad2::RefGainDir (const std::string &db_dir);
std::string  prad2::GainCorrDir(const std::string &db_dir);
std::string  prad2::FindRunFile(const std::string &dir, int run_num, const char *suffix);

std::string  prad2::FindRefGainFile(const std::string &dir, int run_num);   // "_LMS.dat"
RefGainTable prad2::LoadRefGainFile(const std::string &path);
RefGainTable prad2::LoadRefGain    (const std::string &dir, int run_num);

std::string        prad2::FindGainCorrRootFile(const std::string &dir, int run_num); // "_gain_corr.root"
GainCorrTimeSeries prad2::LoadGainCorrTimeSeries(const RunConfig &run_cfg, int run_num,
                                                 bool use_precomputed = false);

FitResult prad2::gain_hist_fitter(TH1F *h, const float &fac);
```

The loaders return an unloaded (default) table/series on any failure.
`LoadGainCorrTimeSeries` reads `prad_<run_num>_gain_corr.root` from
`run_cfg.gain_data_dir`. By default each batch's correction is
`ref / gain_W[module][j]`, where `ref` is the reference run's
(`run_cfg.gain_ref_run`) own `gain_W` averaged over all of its batches;
`use_precomputed=true` uses the stored `gain_corr_W` branch instead.
Missing or non-positive values give 1. `GetCorr(event_num)` returns the
last batch with `event_num_start <= event_num` (the first batch for
earlier events, identity when there are no batches); after loading it is
safe to call from several threads.

Apply with `new_adc2mev = old_adc2mev * corr.ModuleGain(module_id, run_cfg.gain_ref_pmt)`.
For PbWO4 (`module_id = 1000 + W number`), `ref_pmt` selects the PMTs whose
corrections are averaged: `1`, `2`, `3`, `12`, `13`, `23` (default), or `123`.
PbGlass uses its stored `avg`; other IDs return identity. Invalid selections
throw `std::invalid_argument`.

`gain_hist_fitter` fits a Gaussian to the region of `h` above
`fac` × maximum (at least 5 bins around the peak); all fields are zero
when the fit fails.

---

## `GainCorrCompute.h`

LMS/alpha gain-correction batch calculation shared by
`prad2ana_replay_gainCorr` and `prad2ana_online_gain_monitor_base`;
`prad2ana_replay_rawdata` / `prad2ana_replay_recon` call `EnsureGainCorr`
(namespace `analysis`).

Constants: `kGainNW = 1156`, `kGainNLMS = 3`, `kGainLMSIDBase` (LMS1..3 =
3101..3103), `kGainWIDBase` (PbWO4 `module_id` = W number + 1000),
histogram binning `kGainHistBins = 600` over
[`kGainHistMin`, `kGainHistMax`] = [0, 15000], and the default LMS events
per batch `kGainBatchSizeOffline = 4000`, `kGainBatchSizeOnline = 1000`.

`GainBatch` holds one `gain_corr` entry (`batch_id`,
`event_num_start/end`, `n_lms_events`, `n_alpha_events`, `ref_run`,
`unix_time`, `refPMT_ratio[3]`, `gain_W`, `gain_W_ref`, `gain_corr_W`
`[1156][3]`, `fit_mean_ref_lms/alpha[3]`, `fit_mean_W_lms[1156]`).
`GainPlotConfig { enabled, max_hists, w_ids }` and `GainPlotStore`
(owns per-batch clones of the reference and selected module histograms)
are for optional plotting.

```cpp
TH1F       *MakeGainHist(const char *name);                  // detached, caller owns
std::string MakeLMSOutputFile(const std::string &evio_path); // <name without .evio>_lms.root
std::string GainCorrOutputPath(const std::string &db_dir, int run);

int  ReplayLMSFiles(const std::vector<std::string> &evio_files, const std::string &out_dir,
                    int num_threads, const std::string &daq_config,
                    const std::string &daq_map, const std::string &db_dir,
                    std::vector<std::string> *produced = nullptr);
void EnsureGainCorr(const std::vector<std::string> &evio_files, const std::string &db_dir,
                    const std::string &daq_config, const std::string &daq_map,
                    int num_threads);

void SetupGainBranches(TTree *tree, GainBatch &b);
void FlushGainBatch(GainBatch &b, TTree *tree, TH1F *mod_lms[kGainNW],
                    TH1F *ref_lms[kGainNLMS], TH1F *ref_alpha[kGainNLMS],
                    const prad2::RefGainTable &ref_tbl);
bool ComputeGainCorrections(const std::vector<std::string> &lms_files,
                            const std::string &gain_out, int batch_size,
                            int ref_run_num, const prad2::RefGainTable &ref_tbl,
                            const GainPlotConfig *plot_cfg = nullptr,
                            GainPlotStore *plot_store = nullptr);
```

`GainCorrOutputPath` is `GainCorrDir(db_dir)/prad_<run>_gain_corr.root`
(creating the directory). `ReplayLMSFiles` runs `Process_LMSgainFactor`
on every EVIO file into `out_dir`, returns the number of failures and
appends the written files to `*produced`. `EnsureGainCorr` runs
`prad2ana_replay_gainCorr` for every run among `evio_files` that has no
gain_corr file in `GainCorrDir(db_dir)` yet (failures are only reported;
that run then gets identity correction). `ComputeGainCorrections` reads
the `lms_gain` trees of `lms_files`, and every `batch_size` LMS events
`FlushGainBatch` fits the histograms and writes one `gain_corr` entry to
`gain_out`: `refPMT_ratio[j]` = LMS / alpha peak of reference PMT j,
`gain_W` = module LMS peak / `refPMT_ratio`, `gain_W_ref` from `ref_tbl`,
`gain_corr_W = gain_W_ref / gain_W`; a trailing batch with fewer than
`batch_size` LMS events is not written. `unix_time` is anchored on the EPICS
rows of the input files, or on the scaler rows when there are none.

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
