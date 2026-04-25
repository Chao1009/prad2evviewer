//============================================================================
// gem_raw_dump.C — minimal example: read an EVIO file, find the GEM raw
// banks, and print the first few 32-bit words from each bank per event.
//
// What it shows
// -------------
// • How to open an EVIO file with EvChannel and the project's DaqConfig.
// • How to walk the per-event bank tree without invoking any decoder.
// • How to pull the raw uint32_t buffer for a specific bank tag (handy
//   when you want to test custom firmware-format parsing or feed raw
//   words to your own analysis).
//
// The GEM data lives in the SSP/MPD-style banks named in
// daq_config.json under "bank_tags.ssp_raw" — typically 0xE10C and
// 0x0DE9.  We use cfg.is_ssp_bank(tag) so the script picks up whatever
// the config says, and stays correct after future reconfigurations.
//
// Usage
// -----
//   cd build
//   root -l ../analysis/scripts/rootlogon.C
//   .x ../analysis/scripts/gem_raw_dump.C+( \
//       "/data/stage6/prad_023867/prad_023867.evio.00000", 5, 8)
//
//   args: evio_path, max_events (0 = all), n_words_to_print
//============================================================================

#include "EvChannel.h"
#include "DaqConfig.h"
#include "EvStruct.h"
#include "load_daq_config.h"

#include <TString.h>
#include <TSystem.h>

#include <algorithm>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <iostream>
#include <string>

using namespace evc;

//-----------------------------------------------------------------------------
// One ROC's GEM bank: print header line + a row of hex words.
//-----------------------------------------------------------------------------
static void dump_bank(const EvChannel &ch, const EvNode &node,
                      int n_words_show)
{
    const uint32_t *data = ch.GetData(node);
    const size_t    nw   = node.data_words;
    const size_t    show = (n_words_show > 0)
        ? std::min<size_t>(n_words_show, nw)
        : nw;

    // Tag, num, total words.  The tag tells you which bank format to
    // expect (0xE10C / 0x0DE9 are SSP-encoded MPD frames).  data_words
    // is the payload length AFTER the bank header (so a 12-word bank
    // here means 12 × 4 = 48 raw bytes of MPD/APV samples).
    std::printf("    bank tag=0x%04X num=%u depth=%d words=%zu\n",
                node.tag, node.num, node.depth, nw);

    if (nw == 0) return;
    std::printf("      data:");
    for (size_t i = 0; i < show; ++i)
        std::printf(" %08X", data[i]);
    if (show < nw) std::printf(" ... (+%zu)", nw - show);
    std::printf("\n");
}

//-----------------------------------------------------------------------------
// Entry point — runnable via `.x gem_raw_dump.C+(...)` after rootlogon.
//
// The full version takes 4 explicit args (no defaults).  Cling's
// default-arg synthesis SEGVs when an interpreter call leaves trailing
// `const char* = ...` defaults to be filled in alongside int/long
// defaults, so we expose convenience overloads for the common arg
// counts and let them delegate to the full version.
//-----------------------------------------------------------------------------
int gem_raw_dump(const char *evio_path,
                 long        max_events,
                 int         n_words_show,
                 const char *daq_config);

int gem_raw_dump(const char *evio_path)
{ return gem_raw_dump(evio_path, 5L, 8, ""); }
int gem_raw_dump(const char *evio_path, long max_events)
{ return gem_raw_dump(evio_path, max_events, 8, ""); }
int gem_raw_dump(const char *evio_path, long max_events, int n_words_show)
{ return gem_raw_dump(evio_path, max_events, n_words_show, ""); }

int gem_raw_dump(const char *evio_path,
                 long        max_events,
                 int         n_words_show,
                 const char *daq_config)
{
    //---- load DAQ config ----------------------------------------------------
    // PRAD2_DATABASE_DIR is set by rootlogon.C; honour an explicit override
    // if the caller wants to point at a different config (e.g. legacy run).
    // Default is an empty string rather than nullptr — cling SEGVs when
    // marshalling mixed-type defaults if any of them is `nullptr`.
    std::string cfg_path = (daq_config && *daq_config) ? daq_config : "";
    if (cfg_path.empty()) {
        const char *db = std::getenv("PRAD2_DATABASE_DIR");
        cfg_path = std::string(db ? db : "database") + "/daq_config.json";
    }
    DaqConfig cfg;
    if (!load_daq_config(cfg_path, cfg)) {
        std::cerr << "ERROR: cannot load DAQ config from " << cfg_path << "\n";
        return 1;
    }
    std::cout << "DAQ config : " << cfg_path << "\n";
    std::cout << "GEM (SSP) bank tags:";
    for (auto t : cfg.ssp_bank_tags) std::printf(" 0x%04X", t);
    std::cout << "\n";

    //---- open evio ----------------------------------------------------------
    EvChannel ch;
    ch.SetConfig(cfg);
    if (ch.OpenAuto(evio_path) != status::success) {
        std::cerr << "ERROR: cannot open " << evio_path << "\n";
        return 1;
    }
    std::cout << "reading    : " << evio_path << "\n";
    if (max_events > 0)
        std::cout << "limit      : " << max_events << " physics event(s)\n";
    std::cout << std::string(70, '-') << "\n";

    //---- single pass --------------------------------------------------------
    long n_read = 0, n_phys = 0, n_with_gem = 0;
    long total_gem_banks = 0;

    while (ch.Read() == status::success) {
        ++n_read;

        // Scan() parses the event header + child banks into a flat tree.
        // It does NOT decode any payloads (cheap; only header fields are
        // touched).  Skip Scan-failures (malformed events).
        if (!ch.Scan()) continue;

        // Only physics events carry GEM banks; SYNC / EPICS / control
        // events have nothing for us.
        if (ch.GetEventType() != EventType::Physics) continue;
        ++n_phys;

        // Find every GEM bank in this event by asking the config which
        // tags are SSP-style.  NodesForTag returns indices into
        // GetNodes() (O(1) — populated by Scan()).
        const auto &nodes = ch.GetNodes();
        bool printed_header = false;
        long banks_in_event = 0;

        for (auto tag : cfg.ssp_bank_tags) {
            for (int idx : ch.NodesForTag(tag)) {
                const EvNode &node = nodes[idx];
                if (node.data_words == 0) continue;   // empty placeholder

                if (!printed_header) {
                    auto hdr = ch.GetEvHeader();
                    std::printf("Event #%ld  tag=0x%04X  len=%uw\n",
                                n_phys, hdr.tag, hdr.length);
                    printed_header = true;
                }
                // node.parent is the index of the parent ROC bank — its
                // tag identifies which crate this GEM bank came from,
                // useful when more than one GEM crate is active.
                if (node.parent >= 0) {
                    std::printf("  ROC tag=0x%04X\n", nodes[node.parent].tag);
                }
                dump_bank(ch, node, n_words_show);
                ++banks_in_event;
                ++total_gem_banks;
            }
        }
        if (banks_in_event > 0) ++n_with_gem;

        if (max_events > 0 && n_phys >= max_events) break;
    }

    //---- summary ------------------------------------------------------------
    std::cout << std::string(70, '-') << "\n";
    std::cout << "records read              : " << n_read       << "\n";
    std::cout << "physics events            : " << n_phys       << "\n";
    std::cout << "physics events w/ GEM data: " << n_with_gem   << "\n";
    std::cout << "total GEM bank instances  : " << total_gem_banks << "\n";
    ch.Close();
    return 0;
}
