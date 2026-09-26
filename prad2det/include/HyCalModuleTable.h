#pragma once
//=============================================================================
// HyCalModuleTable.h — loader for per-module HyCal JSON tables
//
//     {
//       "default": <value>,
//       "modules": [
//         {"name": "W735", "<key>": <value>},
//         {"module_id": 1735, "<key>": <value>}
//       ]
//     }
//
// Used by HyCalTimeCalib.h and HyCalRfOffsets.h (key "offset_ns", scalar)
// and HyCalTimeCuts.h (key "window", [lo, hi]).  A row resolves its module
// by "name", then by "module_id" (PrimEx ID); rows without a valid value are
// skipped, rows naming no known module are counted as unknown.
//=============================================================================

#include "HyCalSystem.h"
#include "JsonUtil.h"

#include <nlohmann/json.hpp>

#include <ostream>
#include <string>
#include <vector>

namespace prad2 {

struct HyCalModuleTableStats {
    int n_overrides = 0;   // per-module rows applied
    int n_unknown   = 0;   // rows naming no known module
};

// Fill `values` (one entry per module index) from the table at `path`.
// Every module starts at `def`; a valid "default" in the file replaces `def`
// (returned through it) before the per-module rows apply.  `parse(node, v)`
// reads one value and returns false to reject it.  An empty path loads
// nothing; an unreadable file warns on *log (when non-null) and keeps `def`.
template <class T, class Parse>
inline HyCalModuleTableStats LoadHyCalModuleTable(
    const std::string &path, const fdec::HyCalSystem &hycal,
    const char *value_key, Parse parse, T &def, std::vector<T> &values,
    std::ostream *log)
{
    HyCalModuleTableStats stats;
    values.assign(hycal.module_count(), def);
    if (path.empty()) return stats;

    nlohmann::json j;
    std::string err;
    if (!read_json_file(path, j, &err)) {
        if (log) *log << "Warning: " << err << ", using the default for every HyCal module.\n";
        return stats;
    }

    T v{};
    const auto dit = j.find("default");
    if (dit != j.end() && parse(*dit, v)) {
        def = v;
        values.assign(values.size(), def);
    }

    const auto mods = j.find("modules");
    if (mods == j.end() || !mods->is_array()) return stats;
    for (const auto &m : *mods) {
        const auto vit = m.find(value_key);
        if (vit == m.end() || !parse(*vit, v)) continue;

        const fdec::Module *mod = nullptr;
        const auto name = m.find("name");
        if (name != m.end() && name->is_string())
            mod = hycal.module_by_name(name->get<std::string>());
        const auto id = m.find("module_id");
        if (!mod && id != m.end() && id->is_number_integer())
            mod = hycal.module_by_id(id->get<int>());
        if (!mod) { ++stats.n_unknown; continue; }

        values[mod->index] = v;
        ++stats.n_overrides;
    }
    return stats;
}

// Parse functor for scalar tables.
inline bool ParseHyCalTableNumber(const nlohmann::json &node, float &out)
{
    if (!node.is_number()) return false;
    out = node.get<float>();
    return true;
}

} // namespace prad2
