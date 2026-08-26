#pragma once
//=============================================================================
// HyCalDeadModules.h - apply per-run dead-module flags to HyCal.
//
// Dead modules are selected by LoadRunConfig from the independent
// `hycal_dead_modules` table in runinfo/general.json. This helper marks each
// listed module with kDeadModule and all of its side/corner neighbors with
// kDeadNeighbor. Existing dead flags are cleared before applying the list.
//=============================================================================

#include "HyCalSystem.h"

#include <string>
#include <vector>

namespace prad2 {

struct HyCalDeadModuleSummary {
    int n_dead = 0;
    int n_dead_neighbors = 0;
    int n_unknown = 0;
};

inline HyCalDeadModuleSummary ApplyHyCalDeadModules(
    const std::vector<std::string> &dead_module_names,
    fdec::HyCalSystem &hycal)
{
    HyCalDeadModuleSummary summary;
    const int n_modules = hycal.module_count();
    std::vector<bool> is_dead(n_modules, false);

    for (int i = 0; i < n_modules; ++i) {
        auto &mod = hycal.module(i);
        fdec::clear_bit(mod.flag, fdec::kDeadModule);
        fdec::clear_bit(mod.flag, fdec::kDeadNeighbor);
        fdec::clear_bit(mod.flag, fdec::kLeakage);

        mod.ClearVirtNeighbors();
    }

    std::vector<int> dead_indices;
    dead_indices.reserve(dead_module_names.size());
    constexpr int w_grid_size = 34;
    constexpr int beam_hole_min = 16;
    constexpr int beam_hole_max = 17;

    // Mark all the edge modules as having potential leakage
    // Set all the LG modules as dead initially since we do not expect them to be active
    for (int i = 0; i < n_modules; ++i) {
        auto &mod = hycal.module(i);
        if (mod.is_pwo4()) {
            // The physical beam hole is the missing 2x2 block at rows/columns
            // 16-17.  A center within two cells on both axes has at least one
            // of those four virtual cells in its 5x5 leakage window.
            if (mod.row >= beam_hole_min - 2 &&
                mod.row <= beam_hole_max + 2 &&
                mod.column >= beam_hole_min - 2 &&
                mod.column <= beam_hole_max + 2)
                fdec::set_bit(mod.flag, fdec::kLeakage);
            if (mod.row <= 1 || mod.row >= 32 || mod.column <= 1 || mod.column >= 32)
                fdec::set_bit(mod.flag, fdec::kLeakage);
        }

        if (!mod.is_glass()) continue;

        is_dead[i] = true;
        dead_indices.push_back(i);
        fdec::set_bit(mod.flag, fdec::kDeadModule);
        ++summary.n_dead;
    }

    for (const auto &name : dead_module_names) {
        const auto *mod = hycal.module_by_name(name);
        if (!mod || !mod->is_hycal()) {
            ++summary.n_unknown;
            continue;
        }
        if (is_dead[mod->index]) continue;

        is_dead[mod->index] = true;
        dead_indices.push_back(mod->index);
        fdec::set_bit(hycal.module(mod->index).flag, fdec::kDeadModule);
        ++summary.n_dead;
    }

    for (int dead_index : dead_indices) {
        hycal.for_each_neighbor(dead_index, true, [&](int neighbor_index) {
            if (neighbor_index < 0 || neighbor_index >= n_modules ||
                is_dead[neighbor_index]) {
                return;
            }

            auto &neighbor = hycal.module(neighbor_index);
            if (!fdec::test_bit(neighbor.flag, fdec::kDeadNeighbor)) {
                fdec::set_bit(neighbor.flag, fdec::kDeadNeighbor);
                ++summary.n_dead_neighbors;
            }
        });

        const auto &dead_module = hycal.module(dead_index);
        if (!dead_module.is_pwo4()) continue;
        for (int i = 0; i < n_modules; ++i) {
            const auto &module = hycal.module(i);
            if (!module.is_pwo4()) continue;

            double dx, dy;
            hycal.qdist(dead_module, module, dx, dy);
            if (std::abs(dx) < 2.51 && std::abs(dy) < 2.51)
                fdec::set_bit(hycal.module(i).flag, fdec::kLeakage);
        }
    }

    for (int owner_index = 0; owner_index < n_modules; ++owner_index) {
        auto &owner = hycal.module(owner_index);
        if (!owner.is_pwo4() ||
            !fdec::test_bit(owner.flag, fdec::kLeakage))
            continue;

        // add virtual neighbors for dead modules within two layers
        for (int dead_index : dead_indices) {
            const auto &dead = hycal.module(dead_index);
            if (!dead.is_pwo4() || dead_index == owner_index) continue;

            double dx, dy;
            hycal.qdist(owner, dead, dx, dy);
            if (std::abs(dx) < 2.51 && std::abs(dy) < 2.51) {
                owner.AddVirtNeighbor({dead.row, dead.column, dead.index,
                                       dead.x, dead.y, dx, dy,
                                       fdec::ModuleType::PbWO4});
            }
        }
        // Add virtual neighbors within two layers outside the array boundary.
        for (int row = owner.row - 2; row <= owner.row + 2; ++row) {
            for (int column = owner.column - 2;
                 column <= owner.column + 2; ++column) {
                const bool inside = row >= 0 && row < w_grid_size &&
                                    column >= 0 && column < w_grid_size;
                const bool within_two_layers = row >= -2 && row <= 35 &&
                                               column >= -2 && column <= 35;
                if (inside || !within_two_layers) continue;

                owner.AddVirtNeighbor({row, column, -1,
                                       owner.x + (column - owner.column) * owner.size_x,
                                       owner.y - (row - owner.row) * owner.size_y,
                                       static_cast<double>(column - owner.column),
                                       static_cast<double>(owner.row - row),
                                       fdec::ModuleType::PbWO4});
            }
        }

        // Add the physical 2x2 beam-hole cells (rows/columns 16-17) that lie
        // inside this owner's 5x5 window.
        for (int row = beam_hole_min; row <= beam_hole_max; ++row) {
            for (int column = beam_hole_min;
                 column <= beam_hole_max; ++column) {
                if (std::abs(row - owner.row) > 2 ||
                    std::abs(column - owner.column) > 2)
                    continue;
                owner.AddVirtNeighbor({row, column, -1,
                                       owner.x + (column - owner.column) * owner.size_x,
                                       owner.y - (row - owner.row) * owner.size_y,
                                       static_cast<double>(column - owner.column),
                                       static_cast<double>(owner.row - row),
                                       fdec::ModuleType::PbWO4});
            }
        }
    }

    return summary;
}

} // namespace prad2
