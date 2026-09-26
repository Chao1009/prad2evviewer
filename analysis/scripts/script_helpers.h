#pragma once
//============================================================================
// script_helpers.h — small utility functions shared between analysis ACLiC
// scripts (gem_hycal_matching.C, plot_hits_at_hycal.C, …).
//
// Why a header instead of `static` helpers per-script:
//   Cling shares its dictionary scope across all ACLiC-loaded .C files in
//   the same ROOT session.  `static` / anonymous-namespace helpers in two
//   scripts therefore collide with `redefinition of …` errors at the
//   second `.L`.  Marking the helpers `inline` here gives them weak
//   external linkage so the dict-payload merge accepts them.
//
//   Each script just `#include "script_helpers.h"` and the symbols are
//   shared.  Add a new helper here whenever a second script needs it;
//   keep one-script-only helpers private to that script.
//============================================================================

#include "EvChannel.h"
#include "EventData.h"         // prad2::TBIT_sum
#include "Fadc250Data.h"
#include "HyCalCluster.h"
#include "HyCalSystem.h"
#include "PipelineBuilder.h"
#include "SspData.h"
#include "WaveAnalyzer.h"

#include <TString.h>           // Printf() — line-flushed message output

#include <exception>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

// Strip the extension off a path so "out.pdf" becomes "out".  Used by
// scripts that derive a sibling .root output from a user-supplied
// canvas filename.  Leaves the directory alone.
inline std::string strip_extension(const std::string &p)
{
    auto dot = p.find_last_of('.');
    auto slash = p.find_last_of("/\\");
    if (dot == std::string::npos) return p;
    if (slash != std::string::npos && dot < slash) return p;
    return p.substr(0, dot);
}

// Build a script's detector pipeline (DAQ config, runinfo, HyCal, GEM,
// transforms) from its file arguments, nullptr or "" meaning auto-discover
// from runinfo.  A run_num <= 0 is sniffed from the EVIO file name.  Prints
// "[ERROR] <what>" and returns false when the builder throws.
inline bool build_script_pipeline(prad2::Pipeline &pipeline,
                                  const char *evio_path, int run_num,
                                  const char *daq_config, const char *hc_calib_file,
                                  const char *gem_ped_file, const char *gem_cm_file,
                                  const char *hc_map_file, const char *gem_map_file)
{
    const auto arg = [](const char *s) { return std::string(s ? s : ""); };
    try {
        pipeline = prad2::PipelineBuilder()
            .set_daq_config(arg(daq_config))
            .set_hycal_calib(arg(hc_calib_file))
            .set_gem_pedestal(arg(gem_ped_file))
            .set_gem_common_mode(arg(gem_cm_file))
            .set_hycal_map(arg(hc_map_file))
            .set_gem_map(arg(gem_map_file))
            .set_run_number(run_num > 0 ? run_num : -1)
            .set_run_number_from_evio(arg(evio_path))
            .build();
    } catch (const std::exception &e) {
        Printf("[ERROR] %s", e.what());
        return false;
    }
    return true;
}

// HyCal clusters of one event, the scripts' quick way: for every HyCal
// channel the largest-height peak with t_lo < time < t_hi, energized without
// gain correction or time offsets and fed at time 0 (no multi-pulse mode).
// crate_map is DaqConfig::roc_crate_map() (ROC bank tag -> logical crate).
inline std::vector<fdec::ClusterHit>
reconstruct_hycal_event(const fdec::EventData &fadc,
                        const std::unordered_map<uint32_t, int> &crate_map,
                        const fdec::HyCalSystem &hycal, fdec::WaveAnalyzer &ana,
                        fdec::HyCalCluster &clusterer, float t_lo, float t_hi)
{
    clusterer.Clear();
    fdec::WaveResult wres;
    for (int r = 0; r < fadc.nrocs; ++r) {
        const auto &roc = fadc.rocs[r];
        if (!roc.present) continue;
        const auto cit = crate_map.find(roc.tag);
        if (cit == crate_map.end()) continue;     // not in roc_tags
        const int crate = cit->second;
        for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
            const auto &slot = roc.slots[s];
            if (!slot.present) continue;
            for (int c = 0; c < fdec::MAX_CHANNELS; ++c) {
                if (!(slot.channel_mask & (1ull << c))) continue;
                const auto *mod = hycal.module_by_daq(crate, s, c);
                if (!mod || !mod->is_hycal()) continue;
                const auto &cd = slot.channels[c];
                if (cd.nsamples <= 0) continue;
                ana.Analyze(cd.samples, cd.nsamples, wres);
                int   best = -1;
                float best_h = -1.f;
                for (int p = 0; p < wres.npeaks; ++p) {
                    const auto &pk = wres.peaks[p];
                    if (pk.time > t_lo && pk.time < t_hi && pk.height > best_h) {
                        best_h = pk.height; best = p;
                    }
                }
                if (best < 0) continue;
                float energy = static_cast<float>(mod->energize(wres.peaks[best].integral));
                clusterer.AddHit(mod->index, energy, 0.f);
            }
        }
    }
    clusterer.FormClusters();
    std::vector<fdec::ClusterHit> hits;
    clusterer.ReconstructHits(hits);
    return hits;
}

// Counters of for_each_physics_event.
struct EvioScanStats {
    long n_files_open = 0;   // EVIO files opened
    long n_read       = 0;   // EVIO records read
    long n_phys       = 0;   // physics events decoded
    long n_kept       = 0;   // of those, passed to the callback
};

// Decode the physics events of evio_files in order, printing "[file i/N]"
// per file and "[progress]" every 5000 physics events.  Only events with
// trigger_bits exactly 0x100 (the production physics trigger; LMS, alpha,
// cosmic etc. are skipped) are passed to on_event(fadc, ssp).  Stops, without
// closing the current file, once n_phys reaches max_events (> 0), counting
// every physics event so the scanned extent does not depend on the cut.
template <class OnEvent>
inline void for_each_physics_event(const evc::DaqConfig &cfg,
                                   const std::vector<std::string> &evio_files,
                                   long max_events, EvioScanStats &stats, OnEvent &&on_event)
{
    evc::EvChannel ch;
    ch.SetConfig(cfg);
    // Heap-allocated: the decoder structs hold large fixed-size sample
    // arrays, and on the stack they overflow the guard page at function
    // entry (the SEGV comes before any line of the body runs).
    auto fadc_evt = std::make_unique<fdec::EventData>();
    auto ssp_evt  = std::make_unique<ssp::SspEventData>();

    for (const auto &path : evio_files) {
        if (ch.OpenAuto(path) != evc::status::success) {
            Printf("[WARN] skip (cannot open): %s", path.c_str());
            continue;
        }
        ++stats.n_files_open;
        Printf("[file %ld/%zu] %s", stats.n_files_open, evio_files.size(), path.c_str());

        while (ch.Read() == evc::status::success) {
            ++stats.n_read;
            if (!ch.Scan()) continue;
            if (ch.GetEventType() != evc::EventType::Physics) continue;

            for (int i = 0; i < ch.GetNEvents(); ++i) {
                ssp_evt->clear();
                if (!ch.DecodeEvent(i, *fadc_evt, ssp_evt.get())) continue;
                ++stats.n_phys;
                if (fadc_evt->info.trigger_bits == prad2::TBIT_sum) {
                    ++stats.n_kept;
                    on_event(*fadc_evt, *ssp_evt);
                }
                if (max_events > 0 && stats.n_phys >= max_events) return;
            }
            if (stats.n_phys > 0 && stats.n_phys % 5000 == 0)
                Printf("[progress] %ld physics events", stats.n_phys);
        }
        ch.Close();
    }
}

// The scan lines that open a script's summary.
inline void print_scan_summary(const EvioScanStats &stats, size_t n_files)
{
    Printf("--- summary ---");
    Printf("  EVIO files opened     : %ld / %zu", stats.n_files_open, n_files);
    Printf("  EVIO records          : %ld", stats.n_read);
    Printf("  physics events        : %ld", stats.n_phys);
    Printf("  passed trig cut 0x100 : %ld", stats.n_kept);
}
