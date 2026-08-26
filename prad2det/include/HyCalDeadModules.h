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
    }

    std::vector<int> dead_indices;
    dead_indices.reserve(dead_module_names.size());
    // Set all the LG modules as dead initially since we do not expect them to be active
    for (int i = 0; i < n_modules; ++i) {
        const auto &mod = hycal.module(i);
        if (!mod.is_glass()) continue;

        is_dead[i] = true;
        dead_indices.push_back(i);
        fdec::set_bit(hycal.module(i).flag, fdec::kDeadModule);
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
        for (int i = 0; i < n_modules; ++i) {
            const auto &module = hycal.module(i);
            if (!module.is_hycal()) continue;

            double dx, dy;
            hycal.qdist(dead_module, module, dx, dy);
            if (std::abs(dx) < 2.51 && std::abs(dy) < 2.51)
                fdec::set_bit(hycal.module(i).flag, fdec::kLeakage);
        }
    }

    return summary;
}

} // namespace prad2
