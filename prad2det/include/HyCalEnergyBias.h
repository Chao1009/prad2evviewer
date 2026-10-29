#pragma once
//=============================================================================
// HyCalEnergyBias.h — position-dependent HyCal cluster-energy correction.
//
// Each input JSON contains one 5x5 or 10x10 bias map per PbWO4 module:
//   "W1": { "y0": { "x0": bias, ... }, ... }
// where bias = E_rec / E_expected - 1.  The correction factor is therefore
// 1 / (1 + bias).  The optional top-level "grid_size" selects the map
// resolution; files without it retain the legacy 5x5 format.
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
    static constexpr int MAX_GRID_SIZE = 10;
    using Grid = std::array<float, MAX_GRID_SIZE * MAX_GRID_SIZE>;

    float beam_energy = 0.f;
    std::vector<Grid> ee;
    std::vector<Grid> ep;
    int ee_grid_size = GRID_SIZE;
    int ep_grid_size = GRID_SIZE;
    int ee_cells_loaded = 0;
    int ep_cells_loaded = 0;

    float ep_threshold() const
    {
        if (beam_energy <= 0.f) return 0.f;
        return beam_energy - 3.5f * 0.03f * beam_energy /
                             std::sqrt(beam_energy / 1000.f);
    }

    bool is_ep(float cluster_energy) const
    {
        return cluster_energy > ep_threshold();
    }

    // Grid cell of a hit at (x, y) on `module`; the bias-table producer
    // (hycal_energy_bowl_shape) bins with the same rule. False for an invalid
    // position, unsupported grid size, or module without size.
    static bool cell(const Module &module, float x, float y,
                     int grid_size, int &column, int &row)
    {
        if (!std::isfinite(x) || !std::isfinite(y) ||
            module.size_x <= 0. || module.size_y <= 0. ||
            (grid_size != GRID_SIZE && grid_size != MAX_GRID_SIZE))
            return false;
        const auto [local_x, local_y] = module.cell_offset(x, y);
        column = std::clamp(
            static_cast<int>(std::floor((local_x + 0.5f) * grid_size)),
            0, grid_size - 1);
        row = std::clamp(
            static_cast<int>(std::floor((local_y + 0.5f) * grid_size)),
            0, grid_size - 1);
        return true;
    }

    static bool cell(const Module &module, float x, float y,
                     int &column, int &row)
    {
        return cell(module, x, y, GRID_SIZE, column, row);
    }

    float bias(const Module &module, float x, float y,
               float cluster_energy) const
    {
        if (!module.is_pwo4() || module.index < 0 ||
            !std::isfinite(cluster_energy))
            return 0.f;

        const bool use_ep = is_ep(cluster_energy);
        const auto &table = use_ep ? ep : ee;
        const int grid_size = use_ep ? ep_grid_size : ee_grid_size;
        int column, row;
        if (!cell(module, x, y, grid_size, column, row))
            return 0.f;
        if (module.index >= static_cast<int>(table.size())) return 0.f;
        return table[module.index][row * grid_size + column];
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
                                   std::vector<HyCalEnergyBias::Grid> &table,
                                   int &grid_size)
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

    grid_size = HyCalEnergyBias::GRID_SIZE;
    const auto grid_size_it = j.find("grid_size");
    if (grid_size_it != j.end()) {
        if (!grid_size_it->is_number_integer()) {
            std::cerr << "Warning: HyCal energy-bias file " << path
                      << " has an invalid grid_size, using zero bias.\n";
            return 0;
        }
        grid_size = grid_size_it->get<int>();
        if (grid_size != HyCalEnergyBias::GRID_SIZE &&
            grid_size != HyCalEnergyBias::MAX_GRID_SIZE) {
            std::cerr << "Warning: unsupported grid_size " << grid_size
                      << " in HyCal energy-bias file " << path
                      << ", using zero bias.\n";
            grid_size = HyCalEnergyBias::GRID_SIZE;
            return 0;
        }
    }

    int loaded = 0;
    for (int module_index = 0; module_index < hycal.module_count(); ++module_index) {
        const auto &module = hycal.module(module_index);
        if (!module.is_pwo4()) continue;
        auto module_it = j.find(module.name);
        if (module_it == j.end() || !module_it->is_object()) continue;

        for (int row = 0; row < grid_size; ++row) {
            auto row_it = module_it->find("y" + std::to_string(row));
            if (row_it == module_it->end() || !row_it->is_object()) continue;
            for (int column = 0; column < grid_size; ++column) {
                auto cell_it = row_it->find("x" + std::to_string(column));
                if (cell_it == row_it->end() || !cell_it->is_number()) continue;
                const float value = cell_it->get<float>();
                if (!std::isfinite(value)) continue;
                table[module_index][row * grid_size + column] = value;
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
        ee_path, hycal, result->ee, result->ee_grid_size);
    result->ep_cells_loaded = detail::LoadHyCalEnergyBiasFile(
        ep_path, hycal, result->ep, result->ep_grid_size);
    return result;
}

} // namespace fdec
