//=============================================================================
// HyCalClusterDensity.cpp — PRad1 cluster density and S-shape corrections
//=============================================================================

#include "HyCalClusterDensity.h"

#include <algorithm>
#include <cmath>
#include <cctype>
#include <fstream>
#include <iostream>
#include <sstream>

namespace fdec
{

namespace
{
constexpr float kDegToRad = 3.14159265358979323846f / 180.f;
constexpr float kElectronMassMeV = 0.51099895f;
constexpr float kProtonMassMeV = 938.27208816f;

bool parse_key_value(const std::string &line, std::string &key, std::string &value)
{
    const auto eq = line.find('=');
    if (eq == std::string::npos) return false;
    key = HyCalClusterDensity::trim(line.substr(0, eq));
    value = HyCalClusterDensity::trim(line.substr(eq + 1));
    std::replace(key.begin(), key.end(), ' ', '_');
    std::transform(key.begin(), key.end(), key.begin(),
                   [](unsigned char c) { return static_cast<char>(std::tolower(c)); });
    return true;
}

std::vector<float> parse_float_list(std::string text)
{
    for (char &c : text)
        if (c == ',') c = ' ';
    std::istringstream input(text);
    std::vector<float> values;
    float value = 0.f;
    while (input >> value)
        values.push_back(value);
    return values;
}

bool valid_set_index(DensitySet set)
{
    const int idx = static_cast<int>(set);
    return idx >= 0 && idx < static_cast<int>(DensitySet::MaxSets);
}
} // namespace

void HyCalClusterDensity::ParamsSet::Resize(int geometry_groups, int groups,
                                            int position_params)
{
    energy_groups = groups;
    ppars.clear();
    ppars.resize(static_cast<size_t>(geometry_groups * groups));
    for (auto &params : ppars) {
        params.x.assign(static_cast<size_t>(position_params), 0.f);
        params.y.assign(static_cast<size_t>(position_params), 0.f);
    }
}

bool HyCalClusterDensity::Load(DensitySet set, const std::string &position_path,
                               const std::string &energy_path)
{
    if (!valid_set_index(set)) return false;
    auto &params = sets_[static_cast<size_t>(set)];
    params = ParamsSet{};

    bool ok = true;
    if (!position_path.empty())
        ok = process_position_params(position_path, params) && ok;
    if (!energy_path.empty())
        ok = process_energy_params(energy_path, params) && ok;
    return ok && params.loaded();
}

bool HyCalClusterDensity::ChooseSet(DensitySet set)
{
    if (!HasSet(set)) return false;
    current_set_ = set;
    return true;
}

bool HyCalClusterDensity::HasSet(DensitySet set) const
{
    return valid_set_index(set) && sets_[static_cast<size_t>(set)].loaded();
}

void HyCalClusterDensity::CorrectBias(const Module &center, ClusterHit &hit,
                                      float hycal_z, float relative_energy_window,
                                      bool position_correction,
                                      bool energy_correction) const
{
    if (!position_correction && !energy_correction) return;
    if (!HasSet(current_set_)) return;

    const auto &params = sets_[static_cast<size_t>(current_set_)];
    if (hit.energy <= 0.f) return;

    const float theta = std::atan2(std::sqrt(hit.x * hit.x + hit.y * hit.y),
                                   hycal_z);
    const int energy_index = get_energy_index(params, hit.energy, theta,
                                              relative_energy_window);

    if (position_correction && !test_bit(hit.flag, kDenCorr) && params.position_loaded) {
        const int geometry_index = get_geometry_index(center);
        const int idx = energy_index + geometry_index * params.energy_groups;
        if (energy_index >= 0 && geometry_index >= 0 &&
            idx >= 0 && idx < static_cast<int>(params.ppars.size())) {
            const float dx = (hit.x - static_cast<float>(center.x)) /
                             static_cast<float>(center.size_x);
            const float dy = (hit.y - static_cast<float>(center.y)) /
                             static_cast<float>(center.size_y);
            const auto &pars = params.ppars[static_cast<size_t>(idx)];
            hit.x += get_position_bias(pars.x, dx) * static_cast<float>(center.size_x);
            hit.y += get_position_bias(pars.y, dy) * static_cast<float>(center.size_y);
            set_bit(hit.flag, kDenCorr);
        }
    }

    if (energy_correction && !test_bit(hit.flag, kSEneCorr) && params.energy_loaded) {
        const int ep_index = static_cast<int>(params.energy_range.size()) - 1;
        const int ee_index = ep_index + 1;
        const auto *energy_params = static_cast<const Params *>(nullptr);
        if (energy_index == ep_index) {
            auto it = params.epars_ep.find(center.id);
            if (it != params.epars_ep.end()) energy_params = &it->second;
        } else if (energy_index == ee_index) {
            auto it = params.epars_ee.find(center.id);
            if (it != params.epars_ee.end()) energy_params = &it->second;
        }

        if (energy_params && energy_params->x.size() >= NB_ENERGY_PARS) {
            const float dx = (hit.x - static_cast<float>(center.x)) /
                             static_cast<float>(center.size_x);
            const float dy = (hit.y - static_cast<float>(center.y)) /
                             static_cast<float>(center.size_y);
            hit.energy_s_shape_corr = get_energy_bias(energy_params->x, dx, dy,
                                                      hit.energy);
            hit.energy += hit.energy_s_shape_corr;
            set_bit(hit.flag, kSEneCorr);
        }
    }
}

bool HyCalClusterDensity::SelectSetForBeamEnergy(float beam_energy_mev,
                                                 DensitySet &set)
{
    if (beam_energy_mev < 1000.f) {
        set = DensitySet::Set_1GeV;
        return true;
    }
    if (beam_energy_mev > 2000.f && beam_energy_mev < 3000.f) {
        set = DensitySet::Set_2GeV;
        return true;
    }
    if (beam_energy_mev > 3000.f) {
        set = DensitySet::Set_Above3GeV;
        return true;
    }
    return false;
}

const char *HyCalClusterDensity::SetName(DensitySet set)
{
    switch (set) {
        case DensitySet::Set_1GeV: return "1GeV";
        case DensitySet::Set_2GeV: return "2GeV";
        case DensitySet::Set_Above3GeV: return ">3GeV";
        case DensitySet::MaxSets: break;
    }
    return "unknown";
}

bool HyCalClusterDensity::process_position_params(const std::string &path,
                                                  ParamsSet &set) const
{
    std::ifstream input(path);
    if (!input) return false;

    int number_of_params = 0;
    int geometry_groups = 0;
    int energy_groups = 0;
    std::string line;
    while (std::getline(input, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (line == "params") break;

        std::string key, value;
        if (!parse_key_value(line, key, value)) continue;
        if (key == "number_of_params") number_of_params = std::stoi(value);
        else if (key == "geometry_groups") geometry_groups = std::stoi(value);
        else if (key == "energy_groups") energy_groups = std::stoi(value);
        else if (key == "beam_energy") set.beam_energy = std::stof(value);
        else if (key == "energy_range") set.energy_range = parse_float_list(value);
    }

    if (number_of_params <= 0 || geometry_groups <= 0 || energy_groups <= 0)
        return false;
    if (energy_groups != static_cast<int>(set.energy_range.size()) + 1)
        return false;

    set.Resize(geometry_groups, energy_groups, number_of_params);

    size_t index = 0;
    const size_t size = set.ppars.size();
    while (std::getline(input, line) && index < 2 * size) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        std::istringstream row(line);
        std::string label;
        row >> label;
        auto &target = (index < size) ? set.ppars[index].x
                                      : set.ppars[index - size].y;
        for (int i = 0; i < number_of_params && row; ++i)
            row >> target[static_cast<size_t>(i)];
        ++index;
    }

    set.position_loaded = (index == 2 * size);
    return set.position_loaded;
}

bool HyCalClusterDensity::process_energy_params(const std::string &path,
                                                ParamsSet &set) const
{
    std::ifstream input(path);
    if (!input) return false;

    bool ep_section = true;
    bool have_section = false;
    std::string line;
    while (std::getline(input, line)) {
        line = trim(line);
        if (line.empty() || line[0] == '#') continue;
        if (line == "EP_PARAMS") {
            ep_section = true;
            have_section = true;
            continue;
        }
        if (line == "EE_PARAMS") {
            ep_section = false;
            have_section = true;
            continue;
        }
        if (!have_section) continue;

        std::istringstream row(line);
        std::string name;
        Params params;
        params.x.assign(NB_ENERGY_PARS, 0.f);
        row >> name;
        for (int i = 0; i < NB_ENERGY_PARS && row; ++i)
            row >> params.x[static_cast<size_t>(i)];
        const int id = HyCalSystem::name_to_id(name);
        if (id < 0) continue;
        if (ep_section) set.epars_ep[id] = std::move(params);
        else set.epars_ee[id] = std::move(params);
    }

    set.energy_loaded = !set.epars_ep.empty() || !set.epars_ee.empty();
    return set.energy_loaded;
}

int HyCalClusterDensity::get_geometry_index(const Module &center) const
{
    if (center.type == ModuleType::PbWO4) {
        const int row = center.row - 1;
        const int col = center.column - 1;
        if (row < 0 || col < 0) return -1;
        if (row < 17 && col < 17) return 32 + row / 2 * 9 + col / 2;
        if (row < 17) return 32 + row / 2 * 9 + (33 - col) / 2;
        if (col < 17) return 32 + (33 - row) / 2 * 9 + col / 2;
        return 32 + (33 - row) / 2 * 9 + (33 - col) / 2;
    }

    if (center.type == ModuleType::PbGlass) {
        const int row = (center.id - 1) / 30;
        const int col = (center.id - 1) % 30;
        if (row < 6 && col < 24) return row / 3 * 8 + col / 3;
        if (row < 24 && col >= 24) return 16 + (col - 24) / 3 * 8 + row / 3;
        if (row >= 24 && col >= 6) return (29 - row) / 3 * 8 + (29 - col) / 3;
        if (row >= 6 && col < 6) return 16 + (5 - col) / 3 * 8 + (29 - row) / 3;
    }

    return -1;
}

int HyCalClusterDensity::get_energy_index(const ParamsSet &set, float energy,
                                          float theta, float relative_window) const
{
    if (energy <= 0.f || set.energy_range.empty()) return -1;

    const float ep_energy = mott_energy(set.beam_energy, theta);
    if (ep_energy > 0.f && std::abs(energy / ep_energy - 1.f) < relative_window)
        return static_cast<int>(set.energy_range.size()) - 1;

    const float ee_energy = moller_energy(set.beam_energy, theta);
    if (ee_energy > 0.f && std::abs(energy / ee_energy - 1.f) < relative_window)
        return static_cast<int>(set.energy_range.size());

    for (int i = 0; i < static_cast<int>(set.energy_range.size()) - 1; ++i) {
        if (energy > set.energy_range[static_cast<size_t>(i)] &&
            energy <= set.energy_range[static_cast<size_t>(i + 1)])
            return i;
    }

    return -1;
}

float HyCalClusterDensity::get_position_bias(const std::vector<float> &pars,
                                             float d)
{
    if (pars.size() < NB_POSITION_PARS) return 0.f;
    const float d2 = d * d;
    const float d4 = d2 * d2;
    return d * (d2 - 0.25f) * pars[0] *
           (d4 + pars[1] * d2 + pars[2]) * (d2 - pars[3]);
}

float HyCalClusterDensity::get_energy_bias(const std::vector<float> &pars,
                                           float dx, float dy, float energy)
{
    if (pars.size() < NB_ENERGY_PARS || pars[0] == 0.f) return 0.f;
    const float dx2 = dx * dx;
    const float dx4 = dx2 * dx2;
    const float dy2 = dy * dy;
    const float dy4 = dy2 * dy2;
    const float denom = 1.f + pars[1] * dx2 + pars[2] * dy2 +
                        pars[3] * dx2 * dy2 + pars[4] * dx4 +
                        pars[5] * dy4 + pars[6] * dx + pars[7] * dy;
    if (denom == 0.f) return 0.f;
    return energy * (1.f / pars[0] / denom - 1.f);
}

float HyCalClusterDensity::moller_energy(float beam_energy, float theta)
{
    const float a = (beam_energy - kElectronMassMeV) /
                    (beam_energy + kElectronMassMeV);
    const float cth = std::cos(theta);
    const float denom = 1.f - a * cth * cth;
    if (denom == 0.f) return 0.f;
    return kElectronMassMeV * (1.f + a * cth * cth) / denom;
}

float HyCalClusterDensity::mott_energy(float beam_energy, float theta)
{
    const float sin_half = std::sin(theta / 2.f);
    return beam_energy / (1.f + 2.f * beam_energy * sin_half * sin_half /
                                    kProtonMassMeV);
}

std::string HyCalClusterDensity::trim(std::string text)
{
    const auto first = text.find_first_not_of(" \t\r\n");
    if (first == std::string::npos) return {};
    const auto last = text.find_last_not_of(" \t\r\n");
    return text.substr(first, last - first + 1);
}

} // namespace fdec
