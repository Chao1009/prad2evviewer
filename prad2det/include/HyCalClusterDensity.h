#pragma once
//=============================================================================
// HyCalClusterDensity.h — PRad1 cluster density and S-shape corrections
//=============================================================================

#include "HyCalCluster.h"

#include <array>
#include <string>
#include <unordered_map>
#include <vector>

namespace fdec
{

enum class DensitySet : int {
    Set_1GeV = 0,
    Set_2GeV = 1,
    Set_Above3GeV = 2,
    MaxSets  = 3
};

class HyCalClusterDensity
{
public:
    bool Load(DensitySet set, const std::string &position_path,
              const std::string &energy_path);
    bool ChooseSet(DensitySet set);
    bool HasSet(DensitySet set) const;
    DensitySet CurrentSet() const { return current_set_; }

    void CorrectBias(const Module &center, ClusterHit &hit,
                     float hycal_z, float relative_energy_window,
                     bool position_correction, bool energy_correction) const;

    static bool SelectSetForBeamEnergy(float beam_energy_mev, DensitySet &set);
    static const char *SetName(DensitySet set);
    static std::string trim(std::string text);

private:
    static constexpr int NB_POSITION_PARS = 4;
    static constexpr int NB_ENERGY_PARS = 8;

    struct Params {
        std::vector<float> x;
        std::vector<float> y;
    };

    struct ParamsSet {
        float beam_energy = 0.f;
        int energy_groups = 0;
        std::vector<float> energy_range;
        std::vector<Params> ppars;
        std::unordered_map<int, Params> epars_ep;
        std::unordered_map<int, Params> epars_ee;
        bool position_loaded = false;
        bool energy_loaded = false;

        void Resize(int geometry_groups, int groups, int position_params);
        bool loaded() const { return position_loaded || energy_loaded; }
    };

    bool process_position_params(const std::string &path, ParamsSet &set) const;
    bool process_energy_params(const std::string &path, ParamsSet &set) const;

    int get_geometry_index(const Module &center) const;
    int get_energy_index(const ParamsSet &set, float energy,
                         float theta, float relative_window) const;

    static float get_position_bias(const std::vector<float> &pars, float d);
    static float get_energy_bias(const std::vector<float> &pars,
                                 float dx, float dy, float energy);
    static float moller_energy(float beam_energy, float theta);
    static float mott_energy(float beam_energy, float theta);

    std::array<ParamsSet, static_cast<size_t>(DensitySet::MaxSets)> sets_;
    DensitySet current_set_ = DensitySet::Set_1GeV;
};

} // namespace fdec
