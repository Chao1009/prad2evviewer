#pragma once
//=============================================================================
// HyCalRfOffsets.h — per-module HyCal→RF time offset table.
//
// Sized to hycal.module_count() so the per-event inner loop can index by
// `mod->index` — `at(idx)` for direct access, `apply(idx, dt)` for the
// common "subtract offset and re-fold" call.
//
// File format (under database/hycal_rf_offsets/, e.g. 24340.json; row lookup
// and "default" fallback per HyCalModuleTable.h):
//
//     {
//       "default": 0.0,
//       "modules": [
//         {"name": "G123", "offset_ns":  1.21},
//         {"name": "W735", "offset_ns": -0.74}
//       ]
//     }
//
// Offsets live on (−T_RF/2, T_RF/2] —
// anything outside that range can't be recovered from RF alignment alone
// (see docs/analysis_notes/rf_time_reconstruction_plan.md).
//=============================================================================

#include "HyCalModuleTable.h"
#include "RfTime.h"          // RF_PERIOD_NS, FoldRfDelta

#include <iostream>
#include <string>
#include <vector>

namespace prad2 {

struct HyCalRfOffsets {
    float default_off = 0.f;           // ns
    std::vector<float> off;            // per-module-index; size = module_count()
    int   n_overrides = 0;             // count of per-module rows applied
    int   n_unknown   = 0;             // rows naming no known module

    // Offset for a given module index (default if out of range).
    float at(int module_index) const {
        if (module_index < 0 || module_index >= static_cast<int>(off.size()))
            return default_off;
        return off[module_index];
    }

    // Apply per-module offset to an already-folded Δt and re-fold so the
    // result stays on (−T_RF/2, T_RF/2].  NaN passes through unchanged.
    float apply(int module_index, float folded_dt) const {
        return FoldRfDelta(folded_dt - at(module_index));
    }
};

// Uniform table at `def_off` when `path` is empty or unreadable.  Warnings
// go to *log (if non-null).
inline HyCalRfOffsets LoadHyCalRfOffsets(const std::string &path,
                                         const fdec::HyCalSystem &hycal,
                                         float def_off = 0.f,
                                         std::ostream *log = &std::cerr)
{
    HyCalRfOffsets table;
    table.default_off = def_off;
    const auto stats = LoadHyCalModuleTable(path, hycal, "offset_ns",
        ParseHyCalTableNumber, table.default_off, table.off, log);
    table.n_overrides = stats.n_overrides;
    table.n_unknown   = stats.n_unknown;
    return table;
}

} // namespace prad2
