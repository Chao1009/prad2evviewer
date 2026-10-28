// gem/gem_dump.cpp — GEM data diagnostic tool
//
// Decodes SSP/MPD banks from EVIO data, optionally runs GEM reconstruction,
// and prints diagnostic output at each stage of the pipeline.
// Modes and options: see usage() (gem_dump -h).

#include "EvChannel.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "load_daq_config.h"
#include "GemSystem.h"
#include "GemPedestal.h"
#include "GemCluster.h"
#include "GemEventJson.h"
#include "PipelineBuilder.h"
#include "InstallPaths.h"

#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <set>
#include <algorithm>
#include <cstdlib>
#include <filesystem>
#include <getopt.h>
#include <memory>

using namespace evc;

// Event filter for evdump
//
// Single condition:    field=val[,val]:min_dets
// AND conditions:      field=val+field=val:min_dets
//
// Examples:
//   pos=10,11:3                 pos 10 or 11 (any plane) in >=3 dets
//   plane=X+pos=10,11:3         X-plane pos 10 or 11 in >=3 dets
//   match=+Y,-Y:3               beam-hole APVs in >=3 dets
struct FilterCondition {
    std::string field;
    std::vector<std::string> values;
};

struct EvdumpFilter {
    std::vector<FilterCondition> conditions;  // all must match (AND)
    int min_dets = 1;
};

static bool parseEvdumpFilter(const std::string &expr, EvdumpFilter &filt)
{
    // split at last ':'
    auto colon = expr.rfind(':');
    std::string lhs = (colon != std::string::npos) ? expr.substr(0, colon) : expr;
    filt.min_dets = (colon != std::string::npos) ? std::atoi(expr.substr(colon + 1).c_str()) : 1;

    // split conditions by '+'
    filt.conditions.clear();
    size_t start = 0;
    while (start < lhs.size()) {
        auto plus = lhs.find('+', start);
        if (plus == std::string::npos) plus = lhs.size();
        std::string term = lhs.substr(start, plus - start);
        start = plus + 1;

        auto eq = term.find('=');
        if (eq == std::string::npos) {
            std::cerr << "Error: filter term must be field=val[,val]\n";
            return false;
        }
        FilterCondition cond;
        cond.field = term.substr(0, eq);
        std::string vals = term.substr(eq + 1);
        size_t vs = 0;
        while (vs < vals.size()) {
            auto comma = vals.find(',', vs);
            if (comma == std::string::npos) comma = vals.size();
            cond.values.push_back(vals.substr(vs, comma - vs));
            vs = comma + 1;
        }
        if (cond.values.empty()) {
            std::cerr << "Error: empty value in filter term '" << term << "'\n";
            return false;
        }
        filt.conditions.push_back(std::move(cond));
    }
    return !filt.conditions.empty();
}

static std::string getApvField(const gem::ApvConfig &cfg, const std::string &field)
{
    if (field == "pos")    return std::to_string(cfg.plane_index);
    if (field == "plane")  return (cfg.plane_type == 0) ? "X" : "Y";
    if (field == "match")  return cfg.match;
    if (field == "orient") return std::to_string(cfg.orient);
    if (field == "det")    return std::to_string(cfg.det_id);
    return "";
}

static bool apvMatchesFilter(const gem::ApvConfig &cfg, const EvdumpFilter &filt)
{
    // all APV-level conditions must match (skip detector-level like "clusters")
    for (auto &cond : filt.conditions) {
        if (cond.field == "clusters" || cond.field == "cluster") continue;
        std::string val = getApvField(cfg, cond.field);
        bool found = false;
        for (auto &v : cond.values)
            if (v == val) { found = true; break; }
        if (!found) return false;
    }
    return true;
}

static bool checkEvdumpFilter(const EvdumpFilter &filt,
                               gem::GemSystem &sys,
                               const ssp::SspEventData &ssp,
                               int phys_ev = 0)
{
    // separate APV-level conditions from detector-level conditions
    bool has_apv_conds = false;
    bool has_cluster_cond = false;
    int cluster_min = 0;
    for (auto &cond : filt.conditions) {
        if (cond.field == "clusters" || cond.field == "cluster") {
            has_cluster_cond = true;
            cluster_min = cond.values.empty() ? 1 : std::atoi(cond.values[0].c_str());
        } else {
            has_apv_conds = true;
        }
    }

    std::set<int> matching_dets;

    // APV-based matching
    if (has_apv_conds) {
        ssp.forEachApv([&](const ssp::MpdData &mpd, int a, const ssp::ApvData &) {
            int idx = sys.FindApvIndex(mpd.crate_id, mpd.mpd_id, a);
            if (idx < 0) return;
            if (!sys.HasApvZsHits(idx)) return;

            auto &cfg = sys.GetApvConfig(idx);
            if (apvMatchesFilter(cfg, filt))
                matching_dets.insert(cfg.det_id);
        });
    }

    // cluster count matching
    if (has_cluster_cond) {
        std::set<int> cluster_dets;
        for (int d = 0; d < sys.GetNDetectors(); ++d) {
            int nc = static_cast<int>(sys.GetPlaneClusters(d, 0).size()
                                    + sys.GetPlaneClusters(d, 1).size());
            if (nc >= cluster_min)
                cluster_dets.insert(d);
        }
        if (has_apv_conds) {
            // AND: keep only dets that pass both APV and cluster conditions
            std::set<int> both;
            for (auto d : matching_dets)
                if (cluster_dets.count(d)) both.insert(d);
            matching_dets = both;
        } else {
            matching_dets = cluster_dets;
        }
    }

    bool pass = static_cast<int>(matching_dets.size()) >= filt.min_dets;
    if (pass || (phys_ev > 0 && phys_ev % 200 == 0)) {
        std::cerr << "  ev " << phys_ev << ": " << matching_dets.size() << " det(s) matched";
        if (!matching_dets.empty()) {
            std::cerr << " [";
            for (auto d : matching_dets) std::cerr << " GEM" << d;
            std::cerr << " ]";
        }
        std::cerr << (pass ? " PASS\n" : "\n");
    }
    return pass;
}

// Mode: raw — dump decoded SSP APV data
static void dumpRawSsp(const ssp::SspEventData &ssp, int phys_ev)
{
    std::cout << "--- Event " << phys_ev << ": "
              << ssp.nmpds << " MPD(s) ---\n";

    ssp.forEachApv([](const ssp::MpdData &mpd, int a, const ssp::ApvData &apv) {
        std::cout << "  MPD crate=" << mpd.crate_id
                  << " fiber=" << mpd.mpd_id
                  << " APV=" << a
                  << " strips=" << apv.nstrips;
        if (apv.has_online_cm) {
            std::cout << " online_cm=[";
            for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t) {
                if (t) std::cout << ",";
                std::cout << apv.online_cm[t];
            }
            std::cout << "]";
        }
        std::cout << "\n";

        // print first 8 strips with data (or all if ≤ 8)
        int printed = 0;
        for (int s = 0; s < ssp::APV_STRIP_SIZE && printed < 8; ++s) {
            if (!apv.hasStrip(s)) continue;
            printed++;
            std::cout << "    ch[" << std::setw(3) << s << "] =";
            for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t)
                std::cout << " " << std::setw(6) << apv.strips[s][t];
            std::cout << "\n";
        }
        if (apv.nstrips > 8)
            std::cout << "    ... (" << apv.nstrips - 8 << " more strips)\n";
    });
    std::cout << "\n";
}

// Mode: hits — show strip hits after GemSystem processing
static void dumpHits(gem::GemSystem &sys, int phys_ev)
{
    int total_hits = 0;
    for (int d = 0; d < sys.GetNDetectors(); ++d)
        for (int p = 0; p < 2; ++p)
            total_hits += (int)sys.GetPlaneHits(d, p).size();

    std::cout << "--- Event " << phys_ev << ": "
              << total_hits << " strip hit(s) ---\n";

    for (int d = 0; d < sys.GetNDetectors(); ++d) {
        auto &det = sys.GetDetectors()[d];
        for (int p = 0; p < 2; ++p) {
            auto &hits = sys.GetPlaneHits(d, p);
            if (hits.empty()) continue;

            std::cout << "  " << det.name << " " << (p == 0 ? "X" : "Y")
                      << ": " << hits.size() << " hit(s)\n";

            int show = std::min((int)hits.size(), 12);
            for (int i = 0; i < show; ++i) {
                auto &h = hits[i];
                std::cout << "    strip=" << std::setw(4) << h.strip
                          << " pos=" << std::fixed << std::setprecision(2)
                          << std::setw(8) << h.position << "mm"
                          << " charge=" << std::setw(8) << std::setprecision(1)
                          << h.charge
                          << " tbin=" << h.max_timebin;
                if (h.cross_talk) std::cout << " [xtalk]";
                if (!h.ts_adc.empty()) {
                    std::cout << "  ts=[";
                    for (size_t t = 0; t < h.ts_adc.size(); ++t) {
                        if (t) std::cout << ",";
                        std::cout << std::setprecision(0) << h.ts_adc[t];
                    }
                    std::cout << "]";
                }
                std::cout << "\n";
            }
            if ((int)hits.size() > show)
                std::cout << "    ... (" << hits.size() - show << " more)\n";
        }
    }
    std::cout << "\n";
}

// Mode: clusters — show clusters and 2D reconstructed hits
static void dumpClusters(gem::GemSystem &sys, int phys_ev)
{
    auto &all_hits = sys.GetAllHits();
    std::cout << "--- Event " << phys_ev << ": "
              << all_hits.size() << " reconstructed 2D hit(s) ---\n";

    for (int d = 0; d < sys.GetNDetectors(); ++d) {
        auto &det = sys.GetDetectors()[d];

        for (int p = 0; p < 2; ++p) {
            auto &clusters = sys.GetPlaneClusters(d, p);
            if (clusters.empty()) continue;

            std::cout << "  " << det.name << " " << (p == 0 ? "X" : "Y")
                      << ": " << clusters.size() << " cluster(s)\n";

            int show = std::min((int)clusters.size(), 8);
            for (int i = 0; i < show; ++i) {
                auto &cl = clusters[i];
                std::cout << "    pos=" << std::fixed << std::setprecision(2)
                          << std::setw(8) << cl.position << "mm"
                          << " peak=" << std::setprecision(1) << std::setw(8) << cl.peak_charge
                          << " total=" << std::setw(8) << cl.total_charge
                          << " size=" << cl.hits.size()
                          << " tbin=" << cl.max_timebin;
                if (cl.cross_talk) std::cout << " [xtalk]";
                std::cout << "\n";
            }
            if ((int)clusters.size() > show)
                std::cout << "    ... (" << clusters.size() - show << " more)\n";
        }

        auto &hits2d = sys.GetHits(d);
        if (!hits2d.empty()) {
            std::cout << "  " << det.name << " 2D hits: " << hits2d.size() << "\n";
            int show = std::min((int)hits2d.size(), 8);
            for (int i = 0; i < show; ++i) {
                auto &h = hits2d[i];
                std::cout << "    (" << std::fixed << std::setprecision(2)
                          << std::setw(8) << h.x << ", "
                          << std::setw(8) << h.y << ") mm"
                          << "  Qx=" << std::setprecision(0) << h.x_charge
                          << " Qy=" << h.y_charge
                          << " Nx=" << h.x_size << " Ny=" << h.y_size
                          << "\n";
            }
            if ((int)hits2d.size() > show)
                std::cout << "    ... (" << hits2d.size() - show << " more)\n";
        }
    }
    std::cout << "\n";
}

// Mode: evdump — dump single event data to JSON file
static int dumpEventJson(const ssp::SspEventData &ssp,
                         gem::GemSystem &sys,
                         int phys_ev,
                         int32_t trigger_num,
                         uint32_t trigger_bits,
                         const std::string &output_file,
                         bool include_raw = false)
{
    using json = nlohmann::json;

    json root;
    root["event_number"]  = phys_ev;
    root["trigger_number"] = trigger_num;
    root["trigger_bits"]   = trigger_bits;

    // --- raw APV data (optional, before pedestal/CM/zero-sup) ---
    if (include_raw) {
        json raw_arr = json::array();
        ssp.forEachApv([&](const ssp::MpdData &mpd, int a, const ssp::ApvData &apv) {
            json aj;
            aj["crate"] = mpd.crate_id;
            aj["mpd"]   = mpd.mpd_id;
            aj["adc"]   = a;

            int idx = sys.FindApvIndex(mpd.crate_id, mpd.mpd_id, a);
            if (idx >= 0) {
                auto &cfg = sys.GetApvConfig(idx);
                aj["det"]   = cfg.det_id;
                aj["plane"] = (cfg.plane_type == 0) ? "X" : "Y";
                aj["pos"]   = cfg.plane_index;
            }

            json ch_obj = json::object();
            for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s) {
                if (!apv.hasStrip(s)) continue;
                json samples = json::array();
                for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t)
                    samples.push_back(apv.strips[s][t]);
                ch_obj[std::to_string(s)] = samples;
            }
            aj["channels"] = ch_obj;
            raw_arr.push_back(aj);
        });
        root["raw_apvs"] = raw_arr;
    }

    root["zs_apvs"]   = gem::ZsApvsToJson(sys, true, &ssp);
    root["detectors"] = gem::DetectorsToJson(sys, true);

    // binary mode to avoid BOM/encoding issues on Windows
    std::ofstream of(output_file, std::ios::binary);
    if (!of.is_open()) {
        std::cerr << "Error: cannot write " << output_file << "\n";
        return 1;
    }
    of << root.dump(2) << "\n";
    of.close();
    std::cerr << "Written: " << output_file << "\n";
    return 0;
}

// Summary accumulator
struct EventStats {
    int nmpds       = 0;
    int napvs       = 0;
    int nstrips     = 0;
    int nhits_x     = 0;
    int nhits_y     = 0;
    int nclusters_x = 0;
    int nclusters_y = 0;
    int nhits_2d    = 0;
};

static void accumulateStats(const ssp::SspEventData &ssp,
                            gem::GemSystem *sys,
                            EventStats &st)
{
    for (int m = 0; m < ssp.nmpds; ++m) {
        if (!ssp.mpds[m].present) continue;
        st.nmpds++;
        for (int a = 0; a < ssp::MAX_APVS_PER_MPD; ++a) {
            if (!ssp.mpds[m].apvs[a].present) continue;
            st.napvs++;
            st.nstrips += ssp.mpds[m].apvs[a].nstrips;
        }
    }

    if (sys) {
        for (int d = 0; d < sys->GetNDetectors(); ++d) {
            st.nhits_x     += (int)sys->GetPlaneHits(d, 0).size();
            st.nhits_y     += (int)sys->GetPlaneHits(d, 1).size();
            st.nclusters_x += (int)sys->GetPlaneClusters(d, 0).size();
            st.nclusters_y += (int)sys->GetPlaneClusters(d, 1).size();
            st.nhits_2d    += (int)sys->GetHits(d).size();
        }
    }
}

static void usage(const char *prog)
{
    std::cerr
        << "GEM data diagnostic tool\n\n"
        << "Usage:\n"
        << "  " << prog << " <evio_file> [-D <daq_config.json>] [options]\n\n"
        << "Modes (-m, --mode <mode>; default: summary):\n"
        << "  -m raw        Dump raw SSP-decoded APV data\n"
        << "  -m hits       Strip hits after pedestal/CM/zero-sup\n"
        << "  -m clusters   Full reconstruction: clusters + 2D hits\n"
        << "  -m summary    Per-event statistics table\n"
        << "  -m ped        Compute per-strip pedestals → output file\n"
        << "  -m evdump     Dump event(s) with 2D hits to JSON (-n K for first K)\n\n"
        << "Options (short / long):\n"
        << "  -D, --daq-config <file>    DAQ config (auto-searches daq_config.json if omitted)\n"
        << "  -G, --gem-map <file>       GEM map file (default: gem_map.json next to the DAQ config)\n"
        << "  -P, --gem-ped <file>       GEM pedestal file (default: per-run file from database runinfo)\n"
        << "  -o, --output <file>        Output file (default: gem_ped.txt for ped, gem_event.json for evdump)\n"
        << "  -n, --num-events <N>       Max physics events (0 = no cap).  Default:\n"
        << "                             10 for summary/hits/raw/clusters,\n"
        << "                             0 (all) for ped,  1 for evdump.\n"
        << "  -t, --trigger-bit <bit>    Trigger bit filter (-1=all, default)\n"
        << "  -e, --event <N>            Dump only physics event N (1-based)\n"
        << "  -z, --zero-sup <sigma>     Override zero-suppression threshold (default: from gem_map)\n"
        << "  -R, --raw                  Include raw APV data in evdump output\n"
        << "  -f, --filter <expr>        Event filter for evdump: field=val[+field=val]:min_dets\n"
        << "                APV fields: pos, plane (X/Y), match (+Y/-Y), orient, det\n"
        << "                Detector field: clusters=N (>=N clusters per det)\n"
        << "                Use + to AND conditions. Examples:\n"
        << "                  -f plane=X+pos=10,11:3   X-plane pos 10/11 in >=3 dets\n"
        << "                  -f match=+Y,-Y:3         beam-hole APVs in >=3 dets\n"
        << "                  -f clusters=1:3           >=1 cluster in >=3 dets\n"
        << "                  -f clusters=2+match=+Y:2  >=2 clusters AND match APV in >=2 dets\n"
        << "  -h, --help                 Show this help\n";
}

int main(int argc, char *argv[])
{
    if (argc < 2) { usage(argv[0]); return 1; }

    std::string daq_config_file;
    std::string gem_map_file;
    std::string gem_ped_file;
    std::string output_file = "gem_ped.txt";
    std::string mode = "summary";
    std::string filter_expr;
    // -1 = -n not given (per-mode default applied below); 0 = no cap.
    int max_events  = -1;
    int trigger_bit = -1;   // -1 = accept all
    int target_event = 0;   // 0 = disabled
    float zerosup_override = -1.f;  // <0 = use gem_map default
    bool include_raw = false;

    // Long-option names are kept in sync with the Python GEM tools.
    static struct option long_opts[] = {
        {"daq-config",    required_argument, nullptr, 'D'},
        {"gem-map",       required_argument, nullptr, 'G'},
        {"gem-ped",       required_argument, nullptr, 'P'},
        {"output",        required_argument, nullptr, 'o'},
        {"mode",          required_argument, nullptr, 'm'},
        {"num-events",    required_argument, nullptr, 'n'},
        {"trigger-bit",   required_argument, nullptr, 't'},
        {"event",         required_argument, nullptr, 'e'},
        {"filter",        required_argument, nullptr, 'f'},
        {"raw",           no_argument,       nullptr, 'R'},
        {"zero-sup",      required_argument, nullptr, 'z'},
        {"help",          no_argument,       nullptr, 'h'},
        {nullptr, 0, nullptr, 0},
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "D:G:P:o:m:n:t:e:z:f:Rh",
                              long_opts, nullptr)) != -1) {
        switch (opt) {
        case 'D': daq_config_file = optarg; break;
        case 'G': gem_map_file = optarg; break;
        case 'P': gem_ped_file = optarg; break;
        case 'o': output_file = optarg; break;
        case 'm': mode = optarg; break;
        case 'n': max_events = std::atoi(optarg); break;
        case 't': trigger_bit = std::atoi(optarg); break;
        case 'e': target_event = std::atoi(optarg); max_events = 0; break;
        case 'f': filter_expr = optarg; break;
        case 'R': include_raw = true; break;
        case 'z': zerosup_override = std::atof(optarg); break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }
    if (optind >= argc) { usage(argv[0]); return 1; }
    std::string evio_file = argv[optind];

    // Per-mode default for max_events when the user didn't pass -n.
    //   ped        → 0 = all events (pedestals accumulate across run)
    //   evdump     → 1 = dump one matching event (overridden below)
    //   everything → 10 (keep terminal tables short)
    if (max_events < 0) {
        if (mode == "ped")     max_events = 0;
        else if (mode == "evdump") max_events = 1;
        else                   max_events = 10;
    }

    // evdump mode: dump first N events passing filter.
    //   -n K (K>0): dump up to K matching events
    //   -n 0      : dump every matching event
    //   -e N      : dump only physics event N, ignoring filter
    //   -f <expr> : custom APV/cluster filter (default: require 2D hits)
    int evdump_limit = 1;
    int evdump_count = 0;
    EvdumpFilter evdump_filter;
    bool has_evdump_filter = false;
    if (mode == "evdump") {
        if (target_event == 0) {
            evdump_limit = max_events;    // 0 means "all matching"
            max_events = 0;               // let evdump_limit control the loop
        }
        if (output_file == "gem_ped.txt")
            output_file = "gem_event.json";
        if (!filter_expr.empty()) {
            if (!parseEvdumpFilter(filter_expr, evdump_filter))
                return 1;
            has_evdump_filter = true;
            std::cerr << "Filter   : ";
            for (size_t c = 0; c < evdump_filter.conditions.size(); ++c) {
                if (c) std::cerr << " + ";
                auto &cond = evdump_filter.conditions[c];
                std::cerr << cond.field << "=";
                for (size_t i = 0; i < cond.values.size(); ++i)
                    std::cerr << (i ? "," : "") << cond.values[i];
            }
            std::cerr << " in >=" << evdump_filter.min_dets << " dets\n";
        }
    }

    bool need_gem = (mode == "hits" || mode == "clusters" || mode == "summary" || mode == "evdump");
    bool need_cluster = (mode == "clusters" || mode == "summary" || mode == "evdump");

    if (daq_config_file.empty())
        daq_config_file = prad2::find_database_file("daq_config.json");

    DaqConfig daq_cfg;
    if (daq_config_file.empty() || !load_daq_config(daq_config_file, daq_cfg)) {
        std::cerr << "Error: failed to load DAQ config"
                  << (daq_config_file.empty() ? " (not found)" : ": " + daq_config_file)
                  << "\n";
        return 1;
    }
    std::cerr << "DAQ config: " << daq_config_file
              << " (adc_format=" << daq_cfg.adc_format << ")\n";

    // GEM setup as in production; the DAQ config directory is the database
    // dir (gem_map.json, runinfo pedestal/CM files, reconstruction cuts).
    std::unique_ptr<prad2::Pipeline> pipe;
    gem::GemSystem *gem_sys = nullptr;
    std::unique_ptr<gem::GemCluster> gem_clusterer;

    if (need_gem) {
        namespace fs = std::filesystem;
        prad2::PipelineBuilder b;
        b.set_loaded_daq_config(daq_cfg)
         .set_daq_config(daq_config_file)
         .set_database_dir(fs::absolute(daq_config_file).parent_path().string())
         .set_run_number_from_evio(evio_file);
        // The builder resolves relative paths against the database dir, not the CWD.
        if (!gem_map_file.empty()) b.set_gem_map(fs::absolute(gem_map_file).string());
        if (!gem_ped_file.empty()) b.set_gem_pedestal(fs::absolute(gem_ped_file).string());
        try {
            pipe = std::make_unique<prad2::Pipeline>(b.build());
        } catch (const std::exception &e) {
            std::cerr << "Error: " << e.what() << "\n";
            return 1;
        }
        gem_sys = &pipe->gem;

        if (zerosup_override >= 0.f) {
            gem_sys->SetZeroSupThreshold(zerosup_override);
            std::cerr << "Zero-sup : " << zerosup_override << " sigma (override)\n";
        }

        if (need_cluster)
            gem_clusterer = std::make_unique<gem::GemCluster>();
    }

    EvChannel ch;
    ch.SetConfig(daq_cfg);
    if (ch.OpenAuto(evio_file) != status::success) {
        std::cerr << "Error: cannot open " << evio_file << "\n";
        return 1;
    }
    std::cerr << "File     : " << evio_file << "\n\n";

    uint32_t trigger_mask = 0;
    if (trigger_bit >= 0) {
        trigger_mask = 1u << trigger_bit;
        std::cerr << "Trigger  : bit " << trigger_bit
                  << " (mask 0x" << std::hex << trigger_mask << std::dec << ")\n";
    }

    if (mode == "summary") {
        std::cout << std::setw(6) << "ev#"
                  << std::setw(10) << "trigger#"
                  << std::setw(6) << "MPDs"
                  << std::setw(6) << "APVs"
                  << std::setw(8) << "strips"
                  << std::setw(8) << "hits_X"
                  << std::setw(8) << "hits_Y"
                  << std::setw(8) << "clus_X"
                  << std::setw(8) << "clus_Y"
                  << std::setw(8) << "2D_hits"
                  << "\n";
        std::cout << std::string(76, '-') << "\n";
    }

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto ssp_ptr   = std::make_unique<ssp::SspEventData>();
    auto &event    = *event_ptr;
    auto &ssp_evt  = *ssp_ptr;

    int phys_count = 0;
    int ssp_events = 0;
    int ped_events = 0;

    EventStats totals;
    gem::GemPedestal ped_accum;

    // Warn once if full-readout data arrives without a pedestal file: the
    // zero-suppressed reconstruction would produce no hits.  Only relevant
    // in modes that build the pipeline.
    bool ped_warning_emitted = false;
    bool ped_warning_applies = pipe && pipe->gem_pedestal_path.empty();

    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != EventType::Physics) continue;

        for (int i = 0; i < ch.GetNEvents(); ++i) {
            ssp_evt.clear();
            if (!ch.DecodeEvent(i, event, &ssp_evt)) continue;
            phys_count++;

            if (ped_warning_applies && !ped_warning_emitted && ssp_evt.nmpds > 0) {
                if (ssp_evt.hasFullReadout()) {
                    std::cerr
                        << "\n"
                        << "*** WARNING: full-readout GEM data detected, no pedestal file loaded. ***\n"
                        << "    Zero suppression will use the default noise value and produce\n"
                        << "    no hits.  Downstream output will appear empty.  Fix:\n"
                        << "      gem_dump -m ped <run.evio> -o gem_ped.txt\n"
                        << "      gem_dump -m " << mode << " <file.evio> -P gem_ped.txt ...\n\n";
                }
                // Latch either way: only the first event with SSP data is
                // checked; later events are assumed consistent with it.
                ped_warning_emitted = true;
            }

            if (trigger_mask && !(event.info.trigger_bits & trigger_mask))
                continue;

            if (target_event > 0 && phys_count != target_event)
                continue;

            if (ssp_evt.nmpds == 0) {
                if (target_event > 0) {
                    std::cout << "Event " << phys_count << ": no SSP data\n";
                    goto done;
                }
                continue;
            }

            ssp_events++;

            if (gem_sys) {
                gem_sys->Clear();
                gem_sys->ProcessEvent(ssp_evt);
                if (gem_clusterer)
                    gem_sys->Reconstruct(*gem_clusterer);
            }

            if (mode == "raw" || mode == "hits" || mode == "clusters")
                std::cout << "[physics event " << phys_count
                          << " trigger#=" << event.info.trigger_number
                          << " bits=0x" << std::hex << event.info.trigger_bits
                          << std::dec << "]\n";

            if (mode == "ped") {
                if (ped_accum.Accumulate(ssp_evt) > 0 && ++ped_events % 1000 == 0)
                    std::cerr << "  " << ped_events << " events...\r" << std::flush;
            }
            else if (mode == "raw")      dumpRawSsp(ssp_evt, phys_count);
            else if (mode == "hits")     dumpHits(*gem_sys, phys_count);
            else if (mode == "clusters") dumpClusters(*gem_sys, phys_count);
            else if (mode == "evdump") {
                // apply filter (unless -e targets a specific event)
                if (target_event == 0) {
                    if (has_evdump_filter) {
                        if (!checkEvdumpFilter(evdump_filter, *gem_sys, ssp_evt, phys_count))
                            continue;
                    } else {
                        if (gem_sys->GetAllHits().empty())
                            continue;
                    }
                }

                // multi-event: append event number to filename
                std::string out = output_file;
                if (target_event == 0 && evdump_limit != 1) {
                    auto dot = output_file.rfind('.');
                    std::string suffix = "_" + std::to_string(phys_count);
                    if (dot != std::string::npos)
                        out = output_file.substr(0, dot) + suffix + output_file.substr(dot);
                    else
                        out = output_file + suffix;
                }

                if (dumpEventJson(ssp_evt, *gem_sys, phys_count,
                                  event.info.trigger_number,
                                  event.info.trigger_bits, out,
                                  include_raw) != 0)
                    return 1;

                if (evdump_limit > 0 && ++evdump_count >= evdump_limit)
                    goto done;
            }
            else { // summary
                EventStats st;
                accumulateStats(ssp_evt, gem_sys, st);

                std::cout << std::setw(6) << phys_count
                          << std::setw(10) << event.info.trigger_number
                          << std::setw(6) << st.nmpds
                          << std::setw(6) << st.napvs
                          << std::setw(8) << st.nstrips
                          << std::setw(8) << st.nhits_x
                          << std::setw(8) << st.nhits_y
                          << std::setw(8) << st.nclusters_x
                          << std::setw(8) << st.nclusters_y
                          << std::setw(8) << st.nhits_2d
                          << "\n";

                totals.nmpds       += st.nmpds;
                totals.napvs       += st.napvs;
                totals.nstrips     += st.nstrips;
                totals.nhits_x     += st.nhits_x;
                totals.nhits_y     += st.nhits_y;
                totals.nclusters_x += st.nclusters_x;
                totals.nclusters_y += st.nclusters_y;
                totals.nhits_2d    += st.nhits_2d;
            }

            if (target_event > 0)
                goto done;
            if (max_events > 0 && ssp_events >= max_events)
                goto done;
        }
        if (max_events > 0 && ssp_events >= max_events)
            break;
    }

done:
    ch.Close();

    if (mode == "ped") {
        if (ped_accum.NumStrips() == 0) {
            std::cerr << "No full-readout APVs found - pedestals require full-readout events\n";
            return 1;
        }
        std::cerr << "Computed pedestals from " << ped_events << " events, "
                  << ped_accum.NumStrips() << " strips\n";
        int napvs = ped_accum.Write(output_file);
        if (napvs < 0) return 1;
        std::cerr << "Written: " << output_file << " (" << napvs << " APVs, "
                  << ped_accum.NumStrips() << " strips)\n";
        return 0;
    }

    if (mode == "summary" && ssp_events > 0) {
        std::cout << std::string(76, '-') << "\n";
        std::cout << "Totals: " << ssp_events << " events with SSP data"
                  << " (of " << phys_count << " physics events)\n";
        std::cout << "  Avg per event:"
                  << " MPDs=" << std::fixed << std::setprecision(1)
                  << (float)totals.nmpds / ssp_events
                  << " APVs=" << (float)totals.napvs / ssp_events
                  << " strips=" << (float)totals.nstrips / ssp_events
                  << "\n";
        if (gem_sys) {
            std::cout << "  Avg per event:"
                      << " hits_X=" << (float)totals.nhits_x / ssp_events
                      << " hits_Y=" << (float)totals.nhits_y / ssp_events
                      << " clus_X=" << (float)totals.nclusters_x / ssp_events
                      << " clus_Y=" << (float)totals.nclusters_y / ssp_events
                      << " 2D=" << (float)totals.nhits_2d / ssp_events
                      << "\n";
        }
    }
    else if (ssp_events == 0) {
        std::cerr << "No events with SSP data found in " << phys_count
                  << " physics events.\n";
    }

    std::cerr << "Done: parsed " << phys_count << " physics events, "
              << ssp_events << " with SSP data.\n";
    return 0;
}
