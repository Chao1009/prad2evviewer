#pragma once
//=============================================================================
// GemSystem.h — GEM detector system for PRad-II
//
// Manages the GEM detector hierarchy, DAQ↔detector channel mapping,
// pedestal subtraction, common mode correction, zero suppression,
// and strip hit collection.
//
// Usage:
//   gem::GemSystem sys;
//   sys.Init("gem_map.json");
//   sys.LoadPedestals("gem_ped.dat");
//
//   // per-event:
//   sys.Clear();
//   sys.ProcessEvent(ssp_evt);        // decoded SSP data
//   sys.Reconstruct(cluster);         // clustering + 2D hits
//   auto &hits = sys.GetHits(det_id); // reconstructed GEM hits
//=============================================================================

#include <cstdint>
#include <limits>
#include <map>
#include <string>
#include <vector>
#include <array>
#include <unordered_map>

#include "SspData.h"

namespace gem
{

// --- data structures --------------------------------------------------------

struct StripHit {
    int32_t strip       = 0;   // plane-wise strip number
    float   charge      = 0.f; // max charge across time samples
    short   max_timebin = 0;   // time sample with max charge
    float   position    = 0.f; // physical position in mm
    bool    cross_talk  = false;
    std::vector<float> ts_adc;  // all time sample ADC values
};

struct StripCluster {
    float   position     = 0.f; // charge-weighted position (mm)
    float   peak_charge  = 0.f; // highest strip charge in cluster
    float   total_charge = 0.f; // sum of all strip charges
    short   max_timebin  = -1;
    // CRITICAL: must be value-initialized — filterClusters() reads this
    // field, and a garbage value drops clusters non-deterministically.
    bool    cross_talk   = false;
    std::vector<StripHit> hits;

    // --- SBS-style (mpd_gem_view_ssp Cuts) quality variables ----------------
    // Filled for every cluster by GemCluster::reconstructCluster.  Seed =
    // first strip with the maximum `charge` (strict >).  NaN = undefined;
    // keep the default initializers (see cross_talk above).
    //   seed_time      StripMeanTime of the seed strip (ns)
    //   seed_peak_adc  max ts_adc of the seed strip (SBS "seed strip peak ADC";
    //                  equals peak_charge whenever charge = max sample, as on
    //                  every PRad path — a halved valley strip is never the seed)
    //   seed_sum_adc   sum of ts_adc of the seed strip, all samples
    //   max_strip_dt   max over non-seed strips of |StripMeanTime(i) - seed_time|
    //                  (ns); NaN for single-strip clusters
    //   min_ts_corr    min over non-seed strips of TimeSampleCorrelation(seed, i);
    //                  NaN for single-strip clusters
    float   seed_time     = std::numeric_limits<float>::quiet_NaN();
    float   seed_peak_adc = std::numeric_limits<float>::quiet_NaN();
    float   seed_sum_adc  = std::numeric_limits<float>::quiet_NaN();
    float   max_strip_dt  = std::numeric_limits<float>::quiet_NaN();
    float   min_ts_corr   = std::numeric_limits<float>::quiet_NaN();
};

struct GEMHit {
    float x = 0.f, y = 0.f, z = 0.f;
    int   det_id = -1;
    float x_charge = 0.f, y_charge = 0.f;
    float x_peak   = 0.f, y_peak   = 0.f;
    short x_max_timebin = 0, y_max_timebin = 0;
    int   x_size = 0, y_size = 0;

    // --- SBS-style quality of the X/Y cluster pair ---------------------------
    // Filled by GemCluster in both match modes, independent of whether any
    // cut is enabled.  NaN = undefined.
    //   x_time, y_time   seed mean times of the X / Y cluster (ns)
    //   time_diff        x_time - y_time (signed, ns)
    //   adc_asym         (x_peak - y_peak) / (x_peak + y_peak), signed; NaN if
    //                    the sum <= 0.  |adc_asym| is what match_adc_asymmetry cuts
    //   x/y_max_strip_dt, x/y_min_ts_corr   copied from the X / Y StripCluster
    float x_time    = std::numeric_limits<float>::quiet_NaN();
    float y_time    = std::numeric_limits<float>::quiet_NaN();
    float time_diff = std::numeric_limits<float>::quiet_NaN();
    float adc_asym  = std::numeric_limits<float>::quiet_NaN();
    float x_max_strip_dt = std::numeric_limits<float>::quiet_NaN();
    float y_max_strip_dt = std::numeric_limits<float>::quiet_NaN();
    float x_min_ts_corr  = std::numeric_limits<float>::quiet_NaN();
    float y_min_ts_corr  = std::numeric_limits<float>::quiet_NaN();
};

// --- APV pedestal -----------------------------------------------------------

struct ApvPedestal {
    float offset = 0.f;
    float noise  = 5000.f;     // large default → no hits until calibrated
};

// --- configuration ----------------------------------------------------------

struct ApvConfig {
    // DAQ address
    int crate_id    = -1;
    int mpd_id      = -1;
    int adc_ch      = -1;

    // Detector mapping
    int det_id      = -1;      // detector index
    int plane_type  = -1;      // 0=X, 1=Y
    int orient      = 0;       // 0 or 1 (strip reversal)
    int plane_index = -1;      // APV position on plane
    int det_pos     = 0;       // detector position in layer

    // Strip mapping parameters (APV channel → physical strip)
    int  pin_rotate  = 0;    // rotated connector pins (e.g. 16 for pos-11 near beam hole)
    int  shared_pos  = -1;   // effective plane position (-1 = use plane_index)
    bool hybrid_board = true; // hybrid board pin conversion (MPD electronics)
    // match: half-strip intersection constraint for beam hole region.
    //   "" = full strip, "+Y" = above hole, "-Y" = below hole.
    std::string match;

    // Pedestals (per-strip)
    ApvPedestal pedestal[ssp::APV_STRIP_SIZE];

    // Common-mode range from LoadCommonModeRange; not used by the
    // reconstruction (processApv uses the sorting estimator).
    float cm_range_min = 0.f;
    float cm_range_max = 5000.f;
};

struct PlaneConfig {
    int   type      = -1;       // 0=X, 1=Y
    float size      = 0.f;      // mm
    int   n_apvs    = 0;        // number of APVs on this plane
    float pitch     = 0.4f;     // strip pitch (mm)
};

struct DetectorConfig {
    std::string name;
    int    id         = -1;
    std::string type  = "PRADGEM";  // detector type for strip mapping
    PlaneConfig planes[2];          // [0]=X, [1]=Y
};

// --- strip-mapping pipeline (pure, stateless) -------------------------------
//
// Maps one APV25 channel index to the plane-wide strip number using the
// 6-step pipeline shared with the original PRadAnalyzer / mpd_gem_view_ssp
// code.  Parameters mirror ApvConfig so off-line analyses (including the
// Python bindings) can reproduce the layout without instantiating a
// GemSystem.  GemSystem::buildStripMap() delegates to this function —
// single source of truth.
//
//   ch             0..apv_channels-1
//   plane_index    APV position on the plane (0..n_apvs-1)
//   orient         strip-reversal flag (0 or 1)
//   pin_rotate     rotated connector pins (e.g. 16 for pos-11 near beam hole)
//   shared_pos     effective plane position (-1 = use plane_index)
//   hybrid_board   true → apply the MPD hybrid-board pin conversion (step 2)
//   apv_channels   N = channels per APV chip (default 128)
//   readout_center default 32; combined with pin_rotate to form the readout
//                  offset; readout_off<=0 skips step 3
int MapStrip(int ch, int plane_index, int orient,
             int  pin_rotate    = 0,
             int  shared_pos    = -1,
             bool hybrid_board  = true,
             int  apv_channels  = ssp::APV_STRIP_SIZE,
             int  readout_center = 32);

// Convenience: compute plane-wide strip numbers for every channel of an APV.
// Equivalent to calling MapStrip(ch, ...) for ch in [0, apv_channels).
std::vector<int> MapApvStrips(int plane_index, int orient,
                              int  pin_rotate    = 0,
                              int  shared_pos    = -1,
                              bool hybrid_board  = true,
                              int  apv_channels  = ssp::APV_STRIP_SIZE,
                              int  readout_center = 32);

// --- reconstruction config (per-detector knobs for GemCluster) --------------
//
// Lives here (rather than in GemCluster.h) so GemSystem can store one entry
// per detector by value.  Defaults follow the mpd_gem_view_ssp reconstruction
// chain, except that charac_dists is empty (cross-talk flagging off); the
// production distances come from reconstruction_config.json
// (gem.default.charac_dists).
struct ClusterConfig {
    int   min_cluster_hits  = 1;
    int   max_cluster_hits  = 20;
    int   consecutive_thres = 1;    // max gap between consecutive strips
    float split_thres       = 14.f; // charge valley depth for splitting
    float cross_talk_width  = 2.f;  // mm
    float cross_talk_peak_ratio_max = 0.40f; // weak/strong peak ratio upper bound
    std::vector<float> charac_dists;// cross-talk characteristic distances (mm); empty = off

    // XY matching mode: 0 = ADC-sorted 1:1, 1 = full Cartesian with cuts
    int   match_mode          = 1;
    // XY matching cuts (mode 1 only)
    float match_adc_asymmetry = 0.8f;  // max |Qx-Qy|/(Qx+Qy), <0 to disable
    float match_time_diff     = 50.f;  // max |mean_t_x - mean_t_y| in ns, <0 to disable
    float ts_period           = 25.f;  // ns per time sample (also used for strip mean times)

    // --- SBS-style (mpd_gem_view_ssp Cuts) quality cuts; all off by default ---
    // strip level: applied while grouping strips into clusters.  A failing strip
    // is dropped and terminates the current run of consecutive strips (SBS
    // IsGoodStrip semantics).  Strips without a positive sample fail an active
    // time window; empty ts_adc fails an active unimodal cut.
    float strip_time_min = -std::numeric_limits<float>::infinity();  // ns, inclusive strip mean-time window
    float strip_time_max =  std::numeric_limits<float>::infinity();  // (active if either finite)
    bool  strip_unimodal = false;   // SBS "use concave shape cut for strip"
    // cluster level: applied in filterClusters (after cross-talk flagging).
    // NaN quality values (e.g. single-strip cluster) pass.
    float seed_min_peak_adc    = 0.f;   // <= 0 disables (SBS: 30)
    float seed_min_sum_adc     = 0.f;   // <= 0 disables (SBS: 60)
    float strip_time_agreement = -1.f;  // ns; reject if max_strip_dt > value; < 0 disables (SBS: 50)
    float strip_ts_corr_min    = -1.f;  // reject if min_ts_corr < value; <= -1 disables (SBS: 0.7, never applied there)
};

// --- GemSystem class --------------------------------------------------------

class GemCluster;

class GemSystem
{
public:
    GemSystem();
    ~GemSystem();

    // --- initialization -----------------------------------------------------
    void Init(const std::string &map_file);

    // Load per-strip pedestals from the upstream APV-block text format
    // (one "APV crate slot fiber adc" header followed by 128 strip lines:
    // "<strip> <offset> <noise>"). The "slot" field is hardware metadata
    // and is ignored for matching — APVs are keyed by (crate, fiber, adc).
    //
    // crate_remap maps file-side hardware crate IDs to the logical crate
    // IDs used by gem_map.json (e.g. 146 -> 1, 147 -> 2).  Empty map = no
    // remap (file IDs are used directly).
    void LoadPedestals(const std::string &ped_file,
                       const std::map<int, int> &crate_remap = {});

    // Load per-APV common-mode range from the upstream text format:
    // "<crate> <slot> <fiber> <adc> <cm_min> <cm_max>" (one APV per line,
    // slot is ignored).  See LoadPedestals for crate_remap semantics.
    void LoadCommonModeRange(const std::string &cm_file,
                             const std::map<int, int> &crate_remap = {});

    // --- reconstruction config ---------------------------------------------
    // Set per-detector clustering / XY-matching parameters.  Application
    // layer (app_state_init.cpp / Replay) supplies one ClusterConfig per
    // detector after parsing reconstruction_config.json.  cfgs.size() should
    // equal GetNDetectors(); longer vectors are truncated, shorter ones padded
    // with library defaults.  Reconstruct() applies entry [d] before
    // clustering each detector.
    void SetReconConfigs(std::vector<ClusterConfig> cfgs);
    const std::vector<ClusterConfig>& GetReconConfigs() const { return per_det_cfgs_; }

    // --- per-event processing -----------------------------------------------
    void Clear();
    void ProcessEvent(const ssp::SspEventData &evt);
    void Reconstruct(GemCluster &clusterer);

    // --- accessors ----------------------------------------------------------
    int GetNDetectors() const { return static_cast<int>(detectors_.size()); }
    const std::vector<DetectorConfig>& GetDetectors() const { return detectors_; }

    const std::vector<StripHit>& GetPlaneHits(int det, int plane) const;
    const std::vector<StripCluster>& GetPlaneClusters(int det, int plane) const;
    const std::vector<GEMHit>& GetHits(int det) const;
    const std::vector<GEMHit>& GetAllHits() const { return all_hits_; }

    // DAQ→APV index lookup (O(1))
    int FindApvIndex(int crate, int mpd, int adc) const;

    // APV config access (for diagnostics/serialization)
    int GetNApvs() const { return static_cast<int>(apvs_.size()); }
    const ApvConfig& GetApvConfig(int idx) const { return apvs_[idx]; }

    // Compute beam hole X offset from detector center (mm), using match APV strip positions.
    // Returns 0 if no match APVs found. Uses detector 0 as reference (all identical).
    float GetHoleXOffset() const;

    // Active strip extent for a (det_id, plane) in detector-local coords (mm).
    // Returns the position of the minimum and maximum mapped strip number — the
    // tight bounding box of the real readout, which can be smaller than the
    // PlaneConfig.size bbox when split APVs (shared_pos) reuse strip numbers
    // and leave one APV's worth of bbox unused on the inner-edge side.
    // Returns (-size/2, +size/2) when no strips are mapped.
    std::pair<float, float> GetActiveExtent(int det_id, int plane) const;

    // Per-APV zero-suppression results (valid after ProcessEvent)
    bool  IsChannelHit(int apv_idx, int ch) const { return apv_work_[apv_idx].hit_pos[ch]; }
    bool  HasApvZsHits(int apv_idx) const {
        for (int ch = 0; ch < ssp::APV_STRIP_SIZE; ++ch)
            if (apv_work_[apv_idx].hit_pos[ch]) return true;
        return false;
    }
    float GetProcessedAdc(int apv_idx, int ch, int ts) const {
        return apv_work_[apv_idx].raw[ts * ssp::APV_STRIP_SIZE + ch];
    }

    // Configuration
    float GetCommonModeThreshold() const { return common_thres_; }
    float GetZeroSupThreshold()    const { return zerosup_thres_; }
    float GetCrossTalkThreshold()  const { return crosstalk_thres_; }
    bool  GetRejectFirstTimebin()  const { return reject_first_timebin_; }
    bool  GetRejectLastTimebin()   const { return reject_last_timebin_; }
    float GetMinPeakAdc()          const { return min_peak_adc_; }
    float GetMinSumAdc()           const { return min_sum_adc_; }
    void  SetCommonModeThreshold(float v) { common_thres_ = v; }
    void  SetZeroSupThreshold(float v)    { zerosup_thres_ = v; }
    void  SetRejectFirstTimebin(bool v)   { reject_first_timebin_ = v; }
    void  SetRejectLastTimebin(bool v)    { reject_last_timebin_ = v; }
    void  SetMinPeakAdc(float v)          { min_peak_adc_ = v; }
    void  SetMinSumAdc(float v)           { min_sum_adc_ = v; }

private:
    // --- per-APV processing -------------------------------------------------
    void processApv(int apv_idx, const ssp::ApvData &data);
    static float commonModeSorting(const float *buf, int size);
    void collectHits(int apv_idx);

    // --- strip mapping ------------------------------------------------------
    void buildStripMap(int apv_idx);
    // {min, max} mapped strip over the APVs of detectors_[det] on `plane`
    // (match APVs only if match_only); min > max when there are none.
    std::pair<int, int> stripRange(int det, int plane, bool match_only) const;

    // --- detector hierarchy -------------------------------------------------
    std::vector<DetectorConfig> detectors_;
    std::vector<ApvConfig> apvs_;
    std::unordered_map<uint64_t, int> apv_map_;  // prad2::pack_daq_key(crate,mpd,adc) → apv index

    // --- per-APV working data (pre-allocated) -------------------------------
    static constexpr int NUM_HIGH_STRIPS  = 20;  // for sorting CM algorithm

    struct ApvWorkData {
        float raw[ssp::APV_STRIP_SIZE * ssp::SSP_TIME_SAMPLES];
        bool  hit_pos[ssp::APV_STRIP_SIZE];
        int   strip_map[ssp::APV_STRIP_SIZE];    // APV channel → plane strip
    };
    std::vector<ApvWorkData> apv_work_;

    // --- per-plane data (hits + clusters) -----------------------------------
    struct PlaneData {
        std::vector<StripHit> hits;
        std::vector<StripCluster> clusters;
    };
    // plane_data_[det_id][plane_type]
    std::vector<std::array<PlaneData, 2>> plane_data_;

    // --- per-detector reconstructed hits ------------------------------------
    std::vector<std::vector<GEMHit>> det_hits_;
    std::vector<GEMHit> all_hits_;

    // --- global APV parameters -----------------------------------------------
    int   apv_channels_     = ssp::APV_STRIP_SIZE;  // channels per APV chip
    int   readout_center_   = 32;      // default readout mapping center

    // --- thresholds.  zerosup_thres_ runs on both paths (full readout
    // and online-ZS, see processApv).  common_thres_ and crosstalk_thres_
    // are stored but not used by the reconstruction.
    float common_thres_     = 20.f;
    float zerosup_thres_    = 5.f;
    float crosstalk_thres_  = 8.f;

    // --- strip-level cuts (from mpd_gem_view_ssp) -------------------------
    bool  reject_first_timebin_ = true;   // reject if peak at first time bin
    bool  reject_last_timebin_  = true;   // reject if peak at last time bin
    float min_peak_adc_         = 0.f;    // min peak ADC per strip (0=disabled)
    float min_sum_adc_          = 0.f;    // min sum ADC per strip (0=disabled)

    // --- per-detector reconstruction config (clustering + XY matching) ----
    // Sized to detectors_.size() during Init().  Library-default ClusterConfig
    // until SetReconConfigs() supplies parsed values.
    std::vector<ClusterConfig> per_det_cfgs_;
};

} // namespace gem
