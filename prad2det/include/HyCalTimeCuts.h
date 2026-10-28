#pragma once
//=============================================================================
// HyCalTimeCuts.h — per-module HyCal peak-time window with a default
// fallback.
//
// Sized to hycal.module_count() so the per-event inner loop can index it by
// `mod->index` — `at(idx)` for direct access, `in_window(idx, t)` for the
// common "is this peak inside the cut" call.
//
// File format (under database/hycal_time_cut/, e.g. cut_2p1.json; row lookup
// and "default" fallback per HyCalModuleTable.h):
//
//     {
//       "default": [lo_ns, hi_ns],
//       "modules": [
//         {"name": "G123", "window": [lo, hi]},
//         {"name": "W735", "window": [lo, hi]}
//       ]
//     }
//=============================================================================

#include "HyCalModuleTable.h"

#include <iostream>
#include <string>
#include <vector>

namespace prad2 {

struct HyCalTimeCuts {
    struct Window { float lo, hi; };   // ns

    float default_lo = 100.f;          // ns
    float default_hi = 200.f;          // ns
    std::vector<Window> win;           // per-module-index; size = module_count()
    int   n_overrides = 0;             // count of per-module rows applied
    int   n_unknown   = 0;             // rows naming no known module

    Window at(int module_index) const {
        if (module_index < 0 || module_index >= static_cast<int>(win.size()))
            return {default_lo, default_hi};
        return win[module_index];
    }

    bool in_window(int module_index, float t) const {
        const Window w = at(module_index);
        return t > w.lo && t < w.hi;
    }
};

// Uniform table at [def_lo, def_hi] when `path` is empty or unreadable.
// `def_lo`/`def_hi` come from RunConfig::hc_time_win_lo/_hi (i.e. the
// runinfo `time_cuts.hc_time_window`).  Windows that are not two numbers
// are skipped.  Warnings go to *log (if non-null).
inline HyCalTimeCuts LoadHyCalTimeCuts(const std::string &path,
                                       const fdec::HyCalSystem &hycal,
                                       float def_lo, float def_hi,
                                       std::ostream *log = &std::cerr)
{
    HyCalTimeCuts cuts;
    HyCalTimeCuts::Window def{def_lo, def_hi};
    const auto stats = LoadHyCalModuleTable(path, hycal, "window",
        [](const nlohmann::json &node, HyCalTimeCuts::Window &w) {
            if (!node.is_array() || node.size() < 2 ||
                !node[0].is_number() || !node[1].is_number())
                return false;
            w = {node[0].get<float>(), node[1].get<float>()};
            return true;
        },
        def, cuts.win, log);
    cuts.default_lo  = def.lo;
    cuts.default_hi  = def.hi;
    cuts.n_overrides = stats.n_overrides;
    cuts.n_unknown   = stats.n_unknown;
    return cuts;
}

} // namespace prad2
