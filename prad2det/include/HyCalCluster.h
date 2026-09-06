#pragma once
//=============================================================================
// HyCalCluster.h — Island clustering algorithm for HyCal
//
// Ported from PRadIslandCluster / PRadHyCalReconstructor (PRadAnalyzer).
// Operates on index-based hits referencing HyCalSystem modules.
//
// Usage:
//   HyCalSystem sys;
//   sys.Init("hycal_map.json");
//
//   HyCalCluster clusterer(sys);
//   // per-event:
//   clusterer.Clear();
//   clusterer.AddHit(module_index, energy, time);
//   clusterer.FormClusters();
//   for (auto &cl : clusterer.GetClusters()) { ... }
//=============================================================================

#include "HyCalSystem.h"
#include <array>
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>
#include <cmath>

namespace fdec
{

struct IClusterProfile;

// --- configuration ----------------------------------------------------------
struct ClusterConfig {
    // hit thresholds
    float min_module_energy  = 0.f;       // MeV, minimum hit energy
    float min_center_energy  = 10.f;      // MeV, minimum seed energy
    float min_cluster_energy = 50.f;      // MeV, minimum total cluster energy
    int   min_cluster_size   = 1;         // minimum number of hits in cluster

    // island algorithm
    bool  corner_conn        = false;     // include diagonal neighbors in grouping
    int   split_iter         = 6;         // iterations for fraction refinement
    float least_split        = 0.01f;     // minimum fraction to keep a split hit

    // position reconstruction
    float log_weight_thres   = 3.6f;      // W = max(0, thres + ln(E_i/E_tot))

    // energy correction
    bool  non_linear_corr    = true;      // apply per-module energy non-linearity correction
    std::shared_ptr<const IClusterProfile> profile;
    bool  leakage_correction = false;
    int   leakage_iterations = 6;
    float least_leakage_fraction = 0.01f;
    float max_leakage_fraction = 0.30f;
    float leakage_convergence_rel = 0.001f;

    // --- multi-pulse / timing coincidence ------------------------------------
    // Waveform data can produce more than one pulse per module per event.
    // When `seed_time_window > 0`, AddHit() may be called multiple times for
    // the same module (once per pulse); FormClusters() then:
    //   * sorts pulses by energy descending;
    //   * picks the largest unconsumed pulse satisfying min_center_energy as
    //     a cluster seed;
    //   * grows an island via BFS where each neighbour module contributes
    //     ONLY its LARGEST-energy unconsumed pulse whose time lies within
    //     ±seed_time_window of the seed (largest-amplitude is more reliable
    //     than closest-in-time — small pulses have noisy peak times);
    //   * marks every contributing pulse as consumed and repeats from the
    //     next-largest unconsumed pulse — so unmatched pulses stay eligible
    //     to seed additional clusters at different timings within the same
    //     event.
    // When `seed_time_window <= 0` (default), the time field is ignored and
    // the legacy single-pulse-per-module behaviour applies — callers that
    // only ever push one hit per module see no change.
    float seed_time_window   = -1.f;      // ns, ≤ 0 disables timing gating
};

// --- per-event hit ----------------------------------------------------------
struct ModuleHit {
    int   index;        // module index in HyCalSystem
    float energy;       // calibrated energy (MeV)
    float time;         // ADC peaking time (ns)
};

// --- cluster result ---------------------------------------------------------
struct ModuleCluster {
    ModuleHit              center;     // seed module (highest energy local max)
    std::vector<ModuleHit> hits;       // all hits (energy may be split)
    float                  energy = 0.f;
    float                  energy_square = 0.f;
    float                  leakage = 0.f;
    uint32_t               flag   = 0;
    bool                   has_leakage_position = false;
    float                  leakage_x = 0.f;
    float                  leakage_y = 0.f;
    int                    leakage_npos = 0;

    void add_hit(const ModuleHit &h)
    {
        hits.push_back(h);
        energy += h.energy;
    }
};

// --- reconstructed hit position ---------------------------------------------
struct ClusterHit {
    int   center_id;    // PrimEx ID of center module
    float x, y;         // reconstructed position (mm)
    float energy;       // total cluster energy (MeV)
    float time;         // ADC peaking time (ns) of center module
    int   nblocks;      // number of modules in cluster
    int   npos;         // number of modules used in position reconstruction
    uint32_t flag;      // cluster flags
    float linear_corr;    // linearity correction factor (E_corr / E_meas)
    float energy_square;  // raw module-energy sum in the 5x5 area around center
    float leakage;         // energy leakage correction (MeV)
};

// Shower-max depth into the calorimeter face for an EM shower of energy `E`
// (MeV) in the module with PrimEx id `center_id`.  Discriminates W/G by the
// PWO_ID0 boundary.  Returns 0 for E ≤ 0.  Units: mm.
//
//   t = X0 · (ln(E/Ec) − Cf)
//   Cf  = 0.5  for photon-induced showers
//   PWO4:    X0 = 8.6 mm, Ec = 1.1  MeV
//   PbGlass: X0 = 26.7 mm, Ec = 2.84 MeV
//
// (Same formula as the legacy analysis::PhysicsTools::GetShowerDepth, moved
// here so prad2det owns it and the python binding can expose it.)
float shower_depth(int center_id, float energy_mev);

// --- cluster profile (energy sharing lookup) --------------------------------
struct ProfileValue {
    float frac = 0.f;
    float err  = 0.f;
};

// Abstract interface — users can plug in their own profile data.
// Default implementation uses a simple analytical approximation.
struct IClusterProfile {
    virtual ~IClusterProfile() = default;
    virtual ProfileValue GetFractionValue(ModuleType type, float dist,
                                  float energy) const = 0;

    // Returns the fraction of energy at quantized distance `dist` for a cluster
    // of total energy `energy` (MeV) on a module of given type.
    float GetFraction(ModuleType type, float dist, float energy) const
    {
        return GetFractionValue(type, dist, energy).frac;
    }
};

// Simple analytical profile (exponential falloff in Moliere radius units)
struct SimpleProfile : public IClusterProfile {
    ProfileValue GetFractionValue(ModuleType type, float dist,
                          float /*energy*/) const override
    {
        // approximate transverse shower profile in quantized distance units
        // PbWO4 Moliere radius ~20mm ≈ module size, PbGlass ~38mm ≈ module size
        // so quantized distance ~1 corresponds to ~1 Moliere radius
        if (dist < 0.01f) return {0.78f, 0.f}; // center module: ~78% of energy
        float sigma = (type == ModuleType::PbWO4) ? 0.36f : 0.40f;
        return {0.78f * std::exp(-dist * dist / (2.f * sigma * sigma)), 0.f};
    }
};

// read Geant4 simulated shower profile from database
// get the fraction of energy at a given module corresponding to 
// the quantized distance "dist" from cluster center to the module center.
struct Geant4Profile : public IClusterProfile {
    Geant4Profile() = default;

    explicit Geant4Profile(const std::string &path)
    {
        Load(ModuleType::PbWO4, path);
    }

    Geant4Profile(const std::string &pwo_path, const std::string &glass_path)
    {
        Load(ModuleType::PbWO4, pwo_path);
        Load(ModuleType::PbGlass, glass_path);
    }

    bool Load(const std::string &path)
    {
        return Load(ModuleType::PbWO4, path);
    }

    bool Load(ModuleType type, const std::string &path)
    {
        const auto type_index = static_cast<int>(type);
        if (type_index < 0 || type_index >= static_cast<int>(profiles_.size()))
            return false;

        std::ifstream input(path);
        if (!input) return false;

        auto &profile = profiles_[static_cast<size_t>(type_index)];
        std::string line;
        float min_energy, max_energy, energy_step;
        float max_distance, distance_step;
        bool header_found = false;
        while (std::getline(input, line)) {
            if (line.empty() || line.find_first_not_of(" \t") == std::string::npos ||
                line[line.find_first_not_of(" \t")] == '#')
                continue;
            for (char &character : line)
                if (character == ',') character = ' ';
            std::istringstream header(line);
            if (!(header >> min_energy >> max_energy >> energy_step
                        >> max_distance >> distance_step))
                return false;
            header_found = true;
            break;
        }

        if (!header_found || energy_step <= 0.f || distance_step <= 0.f ||
            max_energy < min_energy || max_distance < 0.f)
            return false;

        const int energy_count = static_cast<int>((max_energy - min_energy) /
                                                   energy_step) + 1;
        const int distance_count = static_cast<int>(max_distance /
                                                     distance_step) + 1;
        if (energy_count <= 0 || distance_count <= 0) return false;

        profile.values.assign(static_cast<size_t>(energy_count * distance_count), {});
        profile.min_energy = min_energy;
        profile.max_energy = max_energy;
        profile.energy_step = energy_step;
        profile.max_distance = max_distance;
        profile.distance_step = distance_step;
        profile.energy_count = energy_count;
        profile.distance_count = distance_count;

        int energy_index, distance_index;
        float fraction, error;
        while (std::getline(input, line)) {
            if (line.empty() || line.find_first_not_of(" \t") == std::string::npos ||
                line[line.find_first_not_of(" \t")] == '#')
                continue;
            std::istringstream row(line);
            if (!(row >> energy_index >> distance_index >> fraction >> error))
                continue;
            if (energy_index < 0 || energy_index >= energy_count ||
                distance_index < 0 || distance_index >= distance_count)
                continue;
            profile.values[static_cast<size_t>(energy_index * distance_count +
                                               distance_index)] = {fraction, error};
        }

        profile.loaded = true;
        return true;
    }

    ProfileValue GetFractionValue(ModuleType type, float dist, float energy) const override
    {
        const auto type_index = static_cast<int>(type);
        if (type_index < 0 || type_index >= static_cast<int>(profiles_.size()))
            return SimpleProfile{}.GetFractionValue(type, dist, energy);
        const auto &profile = profiles_[static_cast<size_t>(type_index)];
        if (!profile.loaded)
            return SimpleProfile{}.GetFractionValue(type, dist, energy);
        if (dist < 0.f || dist >= profile.max_distance)
            return {};

        const int distance_index = static_cast<int>(dist /
                                                    profile.distance_step + 0.5f);
        const float normalized_energy = (energy - profile.min_energy) /
                                        profile.energy_step;
        int energy_index = static_cast<int>(normalized_energy);
        if (energy_index < 0) energy_index = 0;
        if (energy_index + 1 >= profile.energy_count)
            energy_index = profile.energy_count - 1;

        const auto value = [&profile, distance_index](int index) {
            return profile.values[static_cast<size_t>(index * profile.distance_count +
                                                      distance_index)];
        };
        if (energy_index == profile.energy_count - 1)
            return value(energy_index);
        const float remainder = normalized_energy -
                                static_cast<float>(energy_index);
        if (remainder < 0.05f) return value(energy_index);
        if (remainder > 0.95f) return value(energy_index + 1);
        const auto lower = value(energy_index);
        const auto upper = value(energy_index + 1);
        return {lower.frac * (1.f - remainder) + upper.frac * remainder,
            lower.err * (1.f - remainder) + upper.err * remainder};
    }

private:
    struct Profile {
        std::vector<ProfileValue> values;
        float min_energy = 0.f;
        float max_energy = 0.f;
        float energy_step = 0.f;
        float max_distance = 0.f;
        float distance_step = 0.f;
        int energy_count = 0;
        int distance_count = 0;
        bool loaded = false;
    };

    std::array<Profile, 4> profiles_;
};


// --- split container (static, reused across calls) --------------------------
static constexpr int SPLIT_MAX_HITS   = 100;
static constexpr int SPLIT_MAX_MAXIMA = 10;

struct SplitContainer {
    float frac[SPLIT_MAX_HITS][SPLIT_MAX_MAXIMA];
    float total[SPLIT_MAX_HITS];

    void sum_frac(int nhits, int nmax)
    {
        for (int i = 0; i < nhits; ++i) {
            total[i] = 0.f;
            for (int j = 0; j < nmax; ++j)
                total[i] += frac[i][j];
        }
    }

    float norm_frac(int imax, int ihit) const
    {
        return (total[ihit] > 0.f) ? frac[ihit][imax] / total[ihit] : 0.f;
    }
};

// --- main clustering class --------------------------------------------------
class HyCalCluster
{
public:
    explicit HyCalCluster(const HyCalSystem &sys);
    ~HyCalCluster();

    // non-copyable
    HyCalCluster(const HyCalCluster &) = delete;
    HyCalCluster &operator=(const HyCalCluster &) = delete;

    // set configuration
    void SetConfig(const ClusterConfig &cfg)
    {
        config_ = cfg;
        if (config_.profile)
            profile_ = config_.profile;
    }
    const ClusterConfig &GetConfig() const   { return config_; }

    // set cluster profile; the shared profile remains valid for all clusterers
    void SetProfile(std::shared_ptr<const IClusterProfile> prof);
    // compatibility overload; transfers ownership of the raw pointer
    void SetProfile(IClusterProfile *prof);

    // --- per-event interface ------------------------------------------------
    void Clear();
    // Push one pulse for `module_index`.  May be called multiple times for
    // the same module — each call records a separate ModuleHit.  See
    // ClusterConfig::seed_time_window for how multi-pulse input is grouped.
    void AddHit(int module_index, float energy, float time);
    void FormClusters();
    void ReconstructHits(std::vector<ClusterHit> &out) const;

    // Reconstruct and return paired (cluster, hit) for clusters passing thresholds.
    // This avoids fragile parallel iteration between GetClusters() and ReconstructHits().
    struct RecoResult {
        const ModuleCluster *cluster;
        ClusterHit hit;
    };
    void ReconstructMatched(std::vector<RecoResult> &out) const;

    // access results
    const std::vector<ModuleCluster> &GetClusters() const { return clusters_; }

    // --- timing-coincidence study tool --------------------------------------
    // For each event, identify seed candidates (the largest pulse passing
    // min_center_energy that hasn't already been claimed by a previous
    // seed in this scan) and emit one row per neighbouring pulse within
    // `max_quantized_dist` of the seed module — WITHOUT applying any
    // timing cut and WITHOUT consuming neighbour pulses across seeds.
    // Use this to histogram dt vs. spatial distance / energy on real
    // data and pick a value for ClusterConfig::seed_time_window before
    // turning the production cut on.
    struct SeedNeighborTiming {
        int    seed_module;     // HyCalSystem module index of the seed
        int    neighbor_module; // HyCalSystem module index of the neighbour
        float  seed_time;       // ns
        float  neighbor_time;   // ns
        float  dt;              // neighbor_time − seed_time (ns)
        float  seed_energy;     // MeV
        float  neighbor_energy; // MeV
        double dx_q;            // quantized distance from seed (module units)
        double dy_q;            // quantized distance from seed (module units)
    };
    void CollectNeighborTiming(std::vector<SeedNeighborTiming> &out,
                               double max_quantized_dist = 5.0) const;

private:
    // island algorithm steps
    void group_hits();
    void grow_island(int seed_idx, int group_id, std::vector<int> &group);
    void split_cluster(const std::vector<int> &group);
    std::vector<int> find_maxima(const std::vector<int> &group) const;
    void split_hits(const std::vector<int> &maxima,
                    const std::vector<int> &group);
    float calculate_energy_square(const ModuleHit &center) const;
    void eval_fraction(const std::vector<int> &maxima,
                       const std::vector<int> &group,
                       SplitContainer &split) const;

    struct LeakagePoint {
        float x = 0.f;
        float y = 0.f;
        float energy = 0.f;
        int   npos = 0;
    };

    struct LeakageHit {
        double x = 0.;
        double y = 0.;
        double dx = 0.;
        double dy = 0.;
        int sector = -1;
        ModuleType type = ModuleType::PbWO4;
        float energy = 0.f;
    };

    void apply_leakage_correction(ModuleCluster &cl) const;
    LeakagePoint reconstruct_leakage_position(const ModuleCluster &cl,
                                              const std::vector<LeakageHit> &leaks,
                                              float total_energy) const;
    double eval_cluster_profile(const LeakagePoint &pos,
                                const ModuleCluster &cl) const;

    // position reconstruction
    ClusterHit reconstruct_pos(const ModuleCluster &cl) const;
    float get_weight(float E, float E_total) const;

    // profile helper
    float get_profile_frac(const ModuleHit &center, const ModuleHit &hit) const;
    ProfileValue get_profile_value_at(float cx, float cy, float cE,
                                      double mx, double my, int msector,
                                      ModuleType type) const;
    ProfileValue get_pwo_profile_value_at(float cx, float cy, float cE,
                                          double mx, double my) const;
    float get_profile_frac_at(float cx, float cy, float cE,
                              const ModuleHit &hit) const;

    const HyCalSystem     &sys_;
    ClusterConfig          config_;
    std::shared_ptr<const IClusterProfile> profile_;

    // per-event data
    std::vector<ModuleHit>              hits_;
    std::vector<std::vector<int>>       groups_;        // groups of hit indices
    std::vector<ModuleCluster>          clusters_;
    std::vector<std::vector<int>>       mod_to_hits_;   // module_index → hit indices
    std::vector<int>                    hit_group_id_;  // hit_index → group_id
    std::vector<bool>                   consumed_;      // hit_index → seeded/claimed
};

} // namespace fdec
