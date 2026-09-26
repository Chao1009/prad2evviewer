#pragma once
//=============================================================================
// HyCalEnergyBias.h — position-dependent HyCal cluster-energy correction.
//
// Each input JSON contains one 5x5 bias map per PbWO4 module:
//   "W1": { "y0": { "x0": bias, ... }, ... }
// where bias = E_rec / E_expected - 1.  The correction factor is therefore
// 1 / (1 + bias).  Missing or invalid entries use zero bias.
//=============================================================================

#include "HyCalSystem.h"
#include "JsonUtil.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

namespace fdec {

struct HyCalEnergyBiasSet {
    float nominal_mev;
    const char *file_prefix;
};

inline HyCalEnergyBiasSet SelectHyCalEnergyBiasSet(float beam_energy)
{
    static constexpr std::array<HyCalEnergyBiasSet, 3> sets = {{
        {700.f,  "0p7GeV"},
        {2200.f, "2p2GeV"},
        {3500.f, "3p5GeV"}
    }};
    return *std::min_element(sets.begin(), sets.end(),
        [beam_energy](const HyCalEnergyBiasSet &left,
                      const HyCalEnergyBiasSet &right) {
            return std::fabs(beam_energy - left.nominal_mev) <
                   std::fabs(beam_energy - right.nominal_mev);
        });
}

struct HyCalEnergyBias {
    static constexpr int GRID_SIZE = 5;
    using Grid = std::array<float, GRID_SIZE * GRID_SIZE>;

    float beam_energy = 0.f;
    std::vector<Grid> ee;
    std::vector<Grid> ep;
    int ee_cells_loaded = 0;
    int ep_cells_loaded = 0;

    float ep_threshold() const
    {
        if (beam_energy <= 0.f) return 0.f;
        return beam_energy - 4.f * 0.03f * beam_energy /
                             std::sqrt(beam_energy / 1000.f);
    }

    bool is_ep(float cluster_energy) const
    {
        return cluster_energy > ep_threshold();
    }

    // Grid cell of a hit at (x, y) on `module`; the bias-table producer
    // (hycal_energy_bowl_shape) bins with the same rule.  False for a
    // non-finite position or a module without size.
    static bool cell(const Module &module, float x, float y,
                     int &column, int &row)
    {
        if (!std::isfinite(x) || !std::isfinite(y) ||
            module.size_x <= 0. || module.size_y <= 0.)
            return false;
        const auto [local_x, local_y] = module.cell_offset(x, y);
        column = std::clamp(
            static_cast<int>(std::floor((local_x + 0.5f) * GRID_SIZE)),
            0, GRID_SIZE - 1);
        row = std::clamp(
            static_cast<int>(std::floor((local_y + 0.5f) * GRID_SIZE)),
            0, GRID_SIZE - 1);
        return true;
    }

    float bias(const Module &module, float x, float y,
               float cluster_energy) const
    {
        int column, row;
        if (!module.is_pwo4() || module.index < 0 ||
            !std::isfinite(cluster_energy) || !cell(module, x, y, column, row))
            return 0.f;

        const auto &table = is_ep(cluster_energy) ? ep : ee;
        if (module.index >= static_cast<int>(table.size())) return 0.f;
        return table[module.index][row * GRID_SIZE + column];
    }

    float correction_factor(const Module &module, float x, float y,
                            float cluster_energy) const
    {
        const float denominator = 1.f + bias(module, x, y, cluster_energy);
        if (!std::isfinite(denominator) || denominator <= 0.f) return 1.f;
        return 1.f / denominator;
    }
};

namespace detail {

inline int LoadHyCalEnergyBiasFile(const std::string &path,
                                   const HyCalSystem &hycal,
                                   std::vector<HyCalEnergyBias::Grid> &table)
{
    table.assign(hycal.module_count(), {});
    if (path.empty()) return 0;

    nlohmann::json j;
    std::string err;
    if (!prad2::read_json_file(path, j, &err)) {
        std::cerr << "Warning: " << err << ", using zero bias.\n";
        return 0;
    }
    if (!j.is_object()) {
        std::cerr << "Warning: HyCal energy-bias file " << path
                  << " is not a JSON object, using zero bias.\n";
        return 0;
    }

    int loaded = 0;
    for (int module_index = 0; module_index < hycal.module_count(); ++module_index) {
        const auto &module = hycal.module(module_index);
        if (!module.is_pwo4()) continue;
        auto module_it = j.find(module.name);
        if (module_it == j.end() || !module_it->is_object()) continue;

        for (int row = 0; row < HyCalEnergyBias::GRID_SIZE; ++row) {
            auto row_it = module_it->find("y" + std::to_string(row));
            if (row_it == module_it->end() || !row_it->is_object()) continue;
            for (int column = 0; column < HyCalEnergyBias::GRID_SIZE; ++column) {
                auto cell_it = row_it->find("x" + std::to_string(column));
                if (cell_it == row_it->end() || !cell_it->is_number()) continue;
                const float value = cell_it->get<float>();
                if (!std::isfinite(value)) continue;
                table[module_index][row * HyCalEnergyBias::GRID_SIZE + column] = value;
                ++loaded;
            }
        }
    }
    return loaded;
}

} // namespace detail

inline std::shared_ptr<const HyCalEnergyBias> LoadHyCalEnergyBias(
    const std::string &ee_path, const std::string &ep_path,
    const HyCalSystem &hycal, float beam_energy)
{
    auto result = std::make_shared<HyCalEnergyBias>();
    result->beam_energy = beam_energy;
    result->ee_cells_loaded = detail::LoadHyCalEnergyBiasFile(
        ee_path, hycal, result->ee);
    result->ep_cells_loaded = detail::LoadHyCalEnergyBiasFile(
        ep_path, hycal, result->ep);
    return result;
}

} // namespace fdec
