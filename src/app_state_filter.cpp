#include "app_state.h"
#include "WaveAnalyzer.h"
#include "HyCalCluster.h"

using json = nlohmann::json;

void AppState::resolveFilterKeys()
{
    filter_wf_keys.clear();
    for (auto &name : waveform_filter.modules) {
        const auto *mod = hycal.module_by_name(name);
        if (!mod || mod->daq.crate < 0) continue;
        for (auto &[roc_tag, crate_id] : roc_to_crate) {
            if (crate_id == mod->daq.crate) {
                filter_wf_keys.insert(fdec::ChannelKey(roc_tag, mod->daq.slot, mod->daq.channel));
                break;
            }
        }
    }

    filter_cl_includes.clear();
    for (auto &name : cluster_filter.includes_modules) {
        const auto *mod = hycal.module_by_name(name);
        if (mod) filter_cl_includes.insert(mod->index);
    }

    filter_cl_centers.clear();
    for (auto &name : cluster_filter.center_modules) {
        const auto *mod = hycal.module_by_name(name);
        if (mod) filter_cl_centers.insert(mod->id);
    }
}

std::string AppState::loadFilter(const json &j)
{
    trigger_type_filter = {};
    waveform_filter = {};
    cluster_filter = {};

    if (j.contains("trigger_type")) {
        auto &tt = j["trigger_type"];
        if (tt.contains("enable")) trigger_type_filter.enable = tt["enable"];
        if (tt.contains("accept") && tt["accept"].is_array()) {
            for (auto &v : tt["accept"])
                trigger_type_filter.accept.push_back(static_cast<uint8_t>(v.get<int>()));
        }
    }

    if (j.contains("waveform")) {
        auto &w = j["waveform"];
        auto &f = waveform_filter;
        if (w.contains("enable"))       f.enable       = w["enable"];
        if (w.contains("n_peaks_min"))  f.n_peaks_min  = w["n_peaks_min"];
        if (w.contains("n_peaks_max"))  f.n_peaks_max  = w["n_peaks_max"];
        auto bound = [&](const char *key, std::optional<float> &b) {
            if (w.contains(key)) b = w[key].get<float>();
        };
        bound("time_min", f.peak.time_min);
        bound("time_max", f.peak.time_max);
        bound("integral_min", f.peak.integral_min);
        bound("integral_max", f.peak.integral_max);
        bound("height_min", f.peak.height_min);
        bound("height_max", f.peak.height_max);
        if (w.contains("modules"))
            for (auto &m : w["modules"])
                f.modules.push_back(m.get<std::string>());
    }

    if (j.contains("clustering")) {
        auto &c = j["clustering"];
        auto &f = cluster_filter;
        if (c.contains("enable"))           f.enable       = c["enable"];
        if (c.contains("n_min"))            f.n_min        = c["n_min"];
        if (c.contains("n_max"))            f.n_max        = c["n_max"];
        if (c.contains("energy_min"))       f.energy_min   = c["energy_min"];
        if (c.contains("energy_max"))       f.energy_max   = c["energy_max"];
        if (c.contains("size_min"))         f.size_min     = c["size_min"];
        if (c.contains("size_max"))         f.size_max     = c["size_max"];
        if (c.contains("includes_min"))     f.includes_min = c["includes_min"];
        if (c.contains("includes_modules"))
            for (auto &m : c["includes_modules"])
                f.includes_modules.push_back(m.get<std::string>());
        if (c.contains("center_modules"))
            for (auto &m : c["center_modules"])
                f.center_modules.push_back(m.get<std::string>());
    }

    resolveFilterKeys();

    std::cerr << "Filter loaded: waveform "
              << (waveform_filter.enable ? "ON" : "off")
              << " (" << filter_wf_keys.size() << " modules)"
              << ", cluster "
              << (cluster_filter.enable ? "ON" : "off") << "\n";
    return "";
}

void AppState::unloadFilter()
{
    trigger_type_filter = {};
    waveform_filter = {};
    cluster_filter = {};
    filter_wf_keys.clear();
    filter_cl_includes.clear();
    filter_cl_centers.clear();
    std::cerr << "Filters unloaded\n";
}

json AppState::filterToJson() const
{
    json r;
    r["active"] = filterActive();
    {
        auto &f = waveform_filter;
        json w;
        w["enable"] = f.enable;
        if (!f.modules.empty()) w["modules"] = f.modules;
        w["n_peaks_min"] = f.n_peaks_min;
        w["n_peaks_max"] = f.n_peaks_max;
        auto bound = [&](const char *key, const std::optional<float> &b) {
            if (b) w[key] = *b;
        };
        bound("time_min", f.peak.time_min);
        bound("time_max", f.peak.time_max);
        bound("integral_min", f.peak.integral_min);
        bound("integral_max", f.peak.integral_max);
        bound("height_min", f.peak.height_min);
        bound("height_max", f.peak.height_max);
        r["waveform"] = w;
    }
    {
        auto &f = cluster_filter;
        json c;
        c["enable"] = f.enable;
        c["n_min"] = f.n_min;
        c["n_max"] = f.n_max;
        if (f.energy_min > 0)    c["energy_min"] = f.energy_min;
        if (f.energy_max < 1e20f) c["energy_max"] = f.energy_max;
        c["size_min"] = f.size_min;
        c["size_max"] = f.size_max;
        if (!f.includes_modules.empty()) {
            c["includes_modules"] = f.includes_modules;
            c["includes_min"] = f.includes_min;
        }
        if (!f.center_modules.empty()) c["center_modules"] = f.center_modules;
        r["clustering"] = c;
    }
    return r;
}

bool AppState::evaluateFilter(fdec::EventData &event,
                              ssp::SspEventData *ssp) const
{
    if (!filterActive()) return true;

    // --- trigger type filter (fast, check first) ---
    if (trigger_type_filter.enable && !trigger_type_filter(event.info.trigger_type))
        return false;

    // --- waveform filter ---
    if (waveform_filter.enable) {
        if (filter_wf_keys.empty()) return false;  // enabled but no modules resolved

        const bool is_adc1881m = (daq_cfg.adc_format == "adc1881m");
        bool any_module_pass = false;
        auto ana = makeAnalyzer();
        fdec::WaveResult wres;

        for (int r = 0; r < event.nrocs && !any_module_pass; ++r) {
            auto &roc = event.rocs[r];
            if (!roc.present) continue;
            fdec::ForEachChannel(roc, [&](int s, int c, const fdec::ChannelData &cd) {
                if (any_module_pass) return;
                if (!filter_wf_keys.count(fdec::ChannelKey(roc.tag, s, c))) return;

                if (is_adc1881m) return;  // no peak analysis for ADC1881M

                ana.SetChannelKey(roc.tag, s, c);
                ana.Analyze(cd.samples, cd.nsamples, wres);
                int n_qual = 0;
                for (int p = 0; p < wres.npeaks; ++p)
                    if (waveform_filter.peak(wres.peaks[p])) n_qual++;
                if (n_qual >= waveform_filter.n_peaks_min &&
                    n_qual <= waveform_filter.n_peaks_max)
                    any_module_pass = true;
            });
        }
        if (!any_module_pass) return false;
    }

    // --- cluster filter ---
    if (cluster_filter.enable) {
        // build clusters from this event (local clusterer, no side effects)
        fdec::HyCalCluster clusterer(hycal);
        clusterer.SetConfig(cluster_cfg);

        auto ana = makeAnalyzer();
        fdec::WaveResult wres;
        feedClusterEvent(event, ana, wres, clusterer);

        clusterer.FormClusters();
        std::vector<fdec::HyCalCluster::RecoResult> results;
        clusterer.ReconstructMatched(results);

        int n_qual = 0;
        for (auto &rr : results) {
            auto &hit = rr.hit;
            auto *cl = rr.cluster;

            if (hit.energy < cluster_filter.energy_min ||
                hit.energy > cluster_filter.energy_max) continue;
            if (hit.nblocks < cluster_filter.size_min ||
                hit.nblocks > cluster_filter.size_max) continue;
            if (!filter_cl_includes.empty()) {
                int count = 0;
                for (auto &h : cl->hits)
                    if (filter_cl_includes.count(h.index)) count++;
                if (count < cluster_filter.includes_min) continue;
            }
            if (!filter_cl_centers.empty()) {
                if (!filter_cl_centers.count(hit.center_id)) continue;
            }
            n_qual++;
        }
        if (n_qual < cluster_filter.n_min || n_qual > cluster_filter.n_max)
            return false;
    }

    return true;
}
