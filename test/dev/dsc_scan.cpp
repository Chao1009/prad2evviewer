// dsc_scan — explore DSC2 scaler banks (0xE115) in EVIO data to identify
// the right (crate, slot, channel) for live-time extraction.
//
// For every 0xE115 bank found, walks up to its parent ROC bank to identify
// which crate it belongs to, parses the per-slot DSC2 layout, and
// reports per-channel gated/ungated counts and the implied live time.
//
// On a typical PRad-II run only a single DSC2 module is read out, but the
// physics-trigger and reference-clock channels both offer a livetime.  This
// tool prints both so the user can pick.
//
// Accepts a single file, a base name (auto-discovers .00000, .00001, ...),
// or a directory (processes all .evio* files sorted by name).

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "Dsc2Decoder.h"
#include "InstallPaths.h"
#include "evio_inputs.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <map>
#include <getopt.h>

using namespace evc;

// Layout in DscData.h, parsed by dsc::Dsc2Decoder::ParsePayload.
static constexpr uint32_t DSC2_BANK_TAG = 0xE115;

static void usage(const char *prog)
{
    std::cerr << "Usage: " << prog << " <input> [-D daq_config.json] [-N max_events] [--all]\n";
}

int main(int argc, char *argv[])
{
    std::string input;
    std::string dcfg = prad2::database_dir() + "/daq_config.json";
    int n_events_max = 0;          // 0 = unlimited
    bool dump_all = false;         // print every DSC2 sighting (not just first/last)

    static struct option lopts[] = {
        {"all", no_argument, nullptr, 'a'},
        {nullptr, 0, nullptr, 0}
    };
    int opt;
    while ((opt = getopt_long(argc, argv, "D:N:ah", lopts, nullptr)) != -1) {
        switch (opt) {
        case 'D': dcfg = optarg; break;
        case 'N': n_events_max = std::atoi(optarg); break;
        case 'a': dump_all = true; break;
        default:
            usage(argv[0]);
            return opt == 'h' ? 0 : 1;
        }
    }
    if (optind < argc) input = argv[optind];
    if (input.empty()) {
        usage(argv[0]);
        return 1;
    }

    DaqConfig cfg;
    if (!load_daq_config(dcfg, cfg)) {
        std::cerr << "Failed to load DAQ config: " << dcfg << "\n";
        return 1;
    }

    // ROC tag → name lookup, for nicer reporting
    std::map<uint32_t, std::string> roc_name;
    std::map<uint32_t, int> roc_crate;
    for (auto &r : cfg.roc_tags) { roc_name[r.tag] = r.name; roc_crate[r.tag] = r.crate; }

    auto files = devtools::discover_evio_inputs(input);
    if (files.empty()) {
        std::cerr << "No EVIO files found for: " << input << "\n";
        return 1;
    }
    std::cerr << "Scanning " << files.size() << " file(s):\n";
    for (auto &f : files) std::cerr << "  " << f << "\n";

    EvChannel ch;
    ch.SetConfig(cfg);

    // Per (parent_tag, slot) keep the most-recent snapshot we saw and a count
    // of how many SYNC/Physics events carried it.
    struct Latest { dsc::DscEventData last{}; uint64_t hits{0}; uint32_t event_tag{0}; };
    std::map<std::pair<uint32_t,int>, Latest> seen;

    // For the very first sighting we also remember the values, to compute
    // delta = last - first → trigger live time over the whole scanned span.
    std::map<std::pair<uint32_t,int>, dsc::DscEventData> first_seen;

    uint64_t scanned = 0, with_dsc = 0;
    uint64_t sync_count = 0;
    uint32_t run_number = 0;

    for (auto &file : files) {
        if (ch.OpenAuto(file) != status::success) {
            std::cerr << "warn: cannot open " << file << ", skipping\n"; continue;
        }
        while (ch.Read() == status::success) {
            if (n_events_max > 0 && (int)scanned >= n_events_max) break;
            if (!ch.Scan()) continue;

            ++scanned;
            auto et = ch.GetEventType();
            if (et != EventType::Sync && et != EventType::Physics) continue;

            // O(1) lookup of every 0xE115 node in this event.
            const auto &idxs = ch.NodesForTag(DSC2_BANK_TAG);
            if (idxs.empty()) continue;
            ++with_dsc;
            if (et == EventType::Sync) ++sync_count;

            const auto &nodes = ch.GetNodes();
            for (int idx : idxs) {
                const EvNode &node = nodes[idx];
                if (node.data_words == 0) continue;
                uint32_t parent_tag = 0;
                if (node.parent >= 0 && node.parent < (int)nodes.size())
                    parent_tag = nodes[node.parent].tag;

                dsc::DscEventData s;
                // Probe the decoder's layouts plus offsets 3 and 5 to catch shifted wrappers.
                if (!dsc::Dsc2Decoder::ParsePayload(ch.GetData(node), node.data_words, s,
                                                    {0, 2, 3, 5})) {
                    if (with_dsc <= 2) {
                        const uint32_t *p = ch.GetData(node);
                        std::cerr << "DEBUG ev#" << scanned
                                  << " parent=0x" << std::hex << std::setw(4)
                                  << std::setfill('0') << parent_tag
                                  << "  bank=0x" << std::setw(4) << node.tag
                                  << "  words=" << std::dec << std::setfill(' ')
                                  << node.data_words << " — could not parse, first 8w:";
                        for (size_t i = 0; i < node.data_words && i < 8; ++i)
                            std::cerr << " 0x" << std::hex << std::setw(8)
                                      << std::setfill('0') << p[i];
                        std::cerr << std::dec << std::setfill(' ') << "\n";
                    }
                    continue;
                }
                auto key = std::make_pair(parent_tag, s.slot);
                auto &L  = seen[key];
                L.last = s;
                L.hits++;
                L.event_tag = ch.GetEvHeader().tag;
                if (first_seen.find(key) == first_seen.end())
                    first_seen[key] = s;
                if (dump_all) {
                    std::cout << "ev#" << scanned
                              << "  parent=0x" << std::hex << std::setw(4)
                              << std::setfill('0') << parent_tag
                              << std::setfill(' ') << std::dec
                              << "  slot=" << s.slot
                              << "  ref_g=" << s.ref_gated
                              << "  ref_u=" << s.ref_ungated
                              << "\n";
                }
            }

            if (run_number == 0 && ch.Sync().run_number != 0)
                run_number = ch.Sync().run_number;
        }
        ch.Close();
        if (n_events_max > 0 && (int)scanned >= n_events_max) break;
    }

    std::cout << "\n=== Scan summary ===\n";
    if (run_number) std::cout << "Run number       : " << run_number << "\n";
    std::cout << "Events scanned   : " << scanned   << "\n";
    std::cout << "Events with DSC2 : " << with_dsc  << "\n";
    std::cout << "SYNC w/ DSC2     : " << sync_count << "\n";
    std::cout << "Unique (crate,slot) DSC2 modules: " << seen.size() << "\n";

    if (seen.empty()) {
        std::cout << "\nNo 0xE115 banks were found.  Check the bank tag in daq_config.json.\n";
        return 0;
    }

    auto pct = [](double x) { std::ostringstream o; o << std::fixed
                                                      << std::setprecision(2) << x; return o.str(); };

    for (auto &[key, L] : seen) {
        uint32_t parent = key.first;
        int slot = key.second;
        auto rn = roc_name.find(parent);
        auto rc = roc_crate.find(parent);

        const dsc::DscEventData &s  = L.last;
        const dsc::DscEventData &s0 = first_seen[key];

        std::cout << "\n--- DSC2 module @ ROC tag 0x" << std::hex << std::setw(4)
                  << std::setfill('0') << parent << std::setfill(' ') << std::dec
                  << "  slot " << slot << "  (payload offset=" << s.offset << ")";
        if (rn != roc_name.end())  std::cout << "  (" << rn->second << ")";
        if (rc != roc_crate.end()) std::cout << "  crate=" << rc->second;
        std::cout << "  hits=" << L.hits << "\n";

        // cumulative (since GO) and scan-span (last - first sighting) ratios
        uint64_t dref_g = (uint64_t)s.ref_gated   - (uint64_t)s0.ref_gated;
        uint64_t dref_u = (uint64_t)s.ref_ungated - (uint64_t)s0.ref_ungated;
        double  rg_u    = (s.ref_ungated > 0) ? (double)s.ref_gated / (double)s.ref_ungated : -1;
        double  rg_d    = (dref_u > 0)        ? (double)dref_g / (double)dref_u            : -1;

        std::cout << "  Ref pair      cum: g=" << std::setw(10) << s.ref_gated
                  << "  u=" << std::setw(10) << s.ref_ungated
                  << "  g/u=" << pct(rg_u * 100) << "%  1-g/u=" << pct((1 - rg_u) * 100) << "%\n";
        if (dref_u > 0)
            std::cout << "                scan: dg=" << std::setw(10) << dref_g
                      << " du=" << std::setw(10) << dref_u
                      << "  g/u=" << pct(rg_d * 100) << "%  1-g/u=" << pct((1 - rg_d) * 100) << "%\n";

        // per-channel TRG / TDC tables — show only channels with any activity
        auto print_table = [&](const char *label,
                               const uint32_t *gated, const uint32_t *ungated,
                               const uint32_t *gated0, const uint32_t *ungated0) {
            bool any = false;
            for (int c = 0; c < dsc::DSC2_NCH; ++c)
                if (ungated[c] != 0) { any = true; break; }
            if (!any) return;

            std::cout << "\n  " << label << " channel scaler counts (cumulative since run start):\n"
                      << "   ch  |  gated         ungated         g/u(%)    1-g/u(%) | "
                      << "scan: dg          du              g/u(%)    1-g/u(%)\n";
            for (int c = 0; c < dsc::DSC2_NCH; ++c) {
                if (ungated[c] == 0) continue;
                double r_c  = (double)gated[c] / (double)ungated[c];
                uint64_t dg = (uint64_t)gated[c]   - (uint64_t)gated0[c];
                uint64_t du = (uint64_t)ungated[c] - (uint64_t)ungated0[c];
                double r_d  = (du > 0) ? (double)dg / (double)du : -1;
                std::cout << "   " << std::setw(2) << c
                          << "  | " << std::setw(12) << gated[c]
                          << "   " << std::setw(13) << ungated[c]
                          << "   " << std::setw(8) << pct(r_c * 100)
                          << "  " << std::setw(8) << pct((1 - r_c) * 100)
                          << "  | " << std::setw(12) << dg
                          << "   " << std::setw(13) << du;
                if (r_d >= 0) std::cout << "   " << std::setw(8) << pct(r_d * 100)
                                        << "  " << std::setw(8) << pct((1 - r_d) * 100);
                std::cout << "\n";
            }
        };

        print_table("TRG", s.trg_gated, s.trg_ungated, s0.trg_gated, s0.trg_ungated);
        print_table("TDC", s.tdc_gated, s.tdc_ungated, s0.tdc_gated, s0.tdc_ungated);
    }

    std::cout << "\n=== Recommendation ===\n"
        << "  • The DSC2 scaler bank lives at parent ROC tag 0x0027 (TI master) — fixed for this DAQ.\n"
        << "  • The slot value comes from the JLab BLKHDR (bits 26:22) printed above; record it as is.\n"
        << "  • Live-time convention: in PRad-II run 024246 the (gated, ungated) ratios are ~0.99,\n"
        << "    which means gated counts LIVE time (gate enabled while NOT busy).  Use the\n"
        << "    formula  live = gated/ungated, NOT 1 - gated/ungated.\n"
        << "  • For per-trigger live time, pick a TRG channel whose ungated count matches the\n"
        << "    expected trigger rate (column 'du' / scan duration).  Channel 2 looks active here.\n"
        << "  • Update database/daq_config.json:\n"
        << "      \"dsc_scaler\": { \"bank_tag\": \"0xE115\", \"slot\": <slot>, \"source\": \"ref|trg|tdc\", \"channel\": <c> }\n";
    return 0;
}
