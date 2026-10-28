#include "PulseTemplateStore.h"
#include "DaqKey.h"
#include "InstallPaths.h"
#include "JsonUtil.h"

#include <nlohmann/json.hpp>
#include <iostream>
#include <cmath>

using namespace fdec;
using nlohmann::json;

namespace {

// Read a {"median": <num>, "mad": <num>} sub-object's median field.
// The fitter writes NaN when no pulses contributed; we propagate NaN.
double read_median(const json &j, const char *key)
{
    if (!j.contains(key)) return std::nan("");
    const auto &sub = j[key];
    if (!sub.is_object() || !sub.contains("median")) return std::nan("");
    const auto &m = sub["median"];
    if (m.is_number()) return m.get<double>();
    return std::nan("");
}

} // anon

void PulseTemplateStore::Clear()
{
    channel_type_.clear();
    by_type_.clear();
    valid_ = false;
}

bool PulseTemplateStore::LoadFromFile(const std::string &path,
                                      const WaveConfig &cfg)
{
    Clear();

    json j;
    std::string err;
    if (!prad2::read_json_file(path, j, &err)) {
        std::cerr << "[PulseTemplateStore] WARN: " << err
                  << " — deconv will fall back to non-deconv mode.\n";
        return false;
    }

    const auto &dcfg = cfg.nnls_deconv;
    const float tr_lo = dcfg.tau_r_min_ns, tr_hi = dcfg.tau_r_max_ns;
    const float tf_lo = dcfg.tau_f_min_ns, tf_hi = dcfg.tau_f_max_ns;

    // ---- channel_id → module_type lookup --------------------------------
    // We only need each channel's category; per-channel τ values are
    // intentionally ignored (the deconvolver uses the per-type aggregate).
    for (auto it = j.begin(); it != j.end(); ++it) {
        const std::string &name = it.key();
        if (!name.empty() && name[0] == '_') continue;     // skip _meta etc.
        const auto &rec = it.value();
        if (!rec.is_object())                  continue;

        if (!rec.contains("channel_id") || !rec["channel_id"].is_string())
            continue;
        int roc = -1, slot = -1, ch = -1;
        if (!ParseChannelKey(rec["channel_id"].get<std::string>(),
                             roc, slot, ch))
            continue;

        if (!rec.contains("module_type") || !rec["module_type"].is_string())
            continue;
        std::string type_name = rec["module_type"].get<std::string>();
        if (type_name.empty() || type_name == "Unknown")    continue;

        channel_type_.emplace(prad2::pack_daq_key(roc, slot, ch),
                              std::move(type_name));
    }

    // ---- per-type templates from `_by_type` block ----------------------
    if (j.contains("_by_type") && j["_by_type"].is_object()) {
        for (auto bt = j["_by_type"].begin(); bt != j["_by_type"].end(); ++bt) {
            const std::string &type_name = bt.key();
            const auto &rec = bt.value();
            if (!rec.is_object()) continue;
            const double tr = read_median(rec, "tau_r_ns");
            const double tf = read_median(rec, "tau_f_ns");
            if (!std::isfinite(tr) || !std::isfinite(tf)) continue;
            if (tr < tr_lo || tr > tr_hi || tf < tf_lo || tf > tf_hi) continue;
            PulseTemplate t{};
            t.tau_r_ns    = static_cast<float>(tr);
            t.tau_f_ns    = static_cast<float>(tf);
            t.is_global   = true;     // a category aggregate, not a per-channel fit
            by_type_.emplace(type_name, t);
        }
    }

    if (by_type_.empty()) {
        std::cerr << "[PulseTemplateStore] WARN: " << path
                  << " yielded no per-type templates"
                  << " — deconv will fall back to non-deconv mode.\n";
        return false;
    }

    valid_ = true;
    std::cerr << "[PulseTemplateStore] loaded " << path
              << ": " << channel_type_.size() << " channels typed,"
              << " per-type:";
    for (const auto &kv : by_type_)
        std::cerr << " " << kv.first
                  << "(τ_r=" << kv.second.tau_r_ns
                  << ",τ_f=" << kv.second.tau_f_ns << ")";
    std::cerr << "\n";
    return true;
}

bool PulseTemplateStore::LoadFromConfig(const WaveConfig &cfg,
                                        const std::string &db_dir)
{
    Clear();
    const auto &file = cfg.nnls_deconv.template_file;
    if (!cfg.nnls_deconv.enabled || file.empty()) return false;
    return LoadFromFile(prad2::resolve_db_path(file, db_dir), cfg);
}

const PulseTemplate *
PulseTemplateStore::type_template(const std::string &type_name) const
{
    auto it = by_type_.find(type_name);
    return (it != by_type_.end()) ? &it->second : nullptr;
}

const PulseTemplate *
PulseTemplateStore::Lookup(int roc_tag, int slot, int channel) const
{
    if (!valid_) return nullptr;

    auto it = channel_type_.find(prad2::pack_daq_key(roc_tag, slot, channel));
    if (it == channel_type_.end()) return nullptr;

    auto bt = by_type_.find(it->second);
    return (bt != by_type_.end()) ? &bt->second : nullptr;
}
