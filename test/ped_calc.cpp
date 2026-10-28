// test/ped_calc.cpp — compute per-channel pedestals from EVIO data
//
// Iterates all events matching a given trigger bit, accumulates raw ADC
// values per (crate, slot, channel), and outputs mean + rms as JSON.
// Usage and options: see usage().

#include "EvChannel.h"
#include "Fadc250Data.h"
#include "load_daq_config.h"
#include "InstallPaths.h"
#include "RunningStats.h"

#include <nlohmann/json.hpp>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <cmath>
#include <cstdlib>
#include <memory>

using namespace evc;
using json = nlohmann::json;

struct ChannelAccum : prad2::RunningStats {
    int crate = 0, slot = 0, channel = 0;
};

static void usage(const char *prog)
{
    std::cerr
        << "Compute per-channel pedestals from EVIO data\n\n"
        << "Usage:\n"
        << "  " << prog << " <evio_file> [-D <daq_config.json>] [options]\n\n"
        << "Options:\n"
        << "  -D <file>   DAQ configuration (auto-searches daq_config.json if omitted)\n"
        << "  -t <bit>    Trigger bit to select (default: 3 = 0x08)\n"
        << "  -o <file>   Output JSON file (default: pedestals_out.json)\n"
        << "  -n <N>      Max events to process (default: all)\n\n"
        << "Trigger bits (PRad-I):\n"
        << "  0 = PHYS_LeadGlassSum (0x01)\n"
        << "  1 = PHYS_TotalSum     (0x02)\n"
        << "  2 = LMS_Led           (0x04)\n"
        << "  3 = LMS_Alpha         (0x08)\n"
        << "PRad-II bits (24 = LMS, 25 = Alpha, ...): see database/trigger_bits.json\n";
}

int main(int argc, char *argv[])
{
    if (argc < 2) { usage(argv[0]); return 1; }

    std::string evio_file;
    std::string daq_config_file;
    std::string output_file = "pedestals_out.json";
    int trigger_bit = 3;  // PRad-I LMS_Alpha (0x08)
    int max_events  = 0;  // 0 = all

    for (int i = 1; i < argc; ++i) {
        std::string a = argv[i];
        if (a == "-D" && i + 1 < argc) { daq_config_file = argv[++i]; }
        else if (a == "-t" && i + 1 < argc) { trigger_bit = std::atoi(argv[++i]); }
        else if (a == "-o" && i + 1 < argc) { output_file = argv[++i]; }
        else if (a == "-n" && i + 1 < argc) { max_events = std::atoi(argv[++i]); }
        else if (evio_file.empty()) { evio_file = a; }
    }

    if (evio_file.empty()) { usage(argv[0]); return 1; }

    uint32_t trigger_mask = 1u << trigger_bit;
    std::cerr << "Input    : " << evio_file << "\n"
              << "Output   : " << output_file << "\n"
              << "Trigger  : bit " << trigger_bit << " (mask 0x"
              << std::hex << trigger_mask << std::dec << ")\n";

    if (daq_config_file.empty())
        daq_config_file = prad2::find_database_file("daq_config.json");

    DaqConfig daq_cfg;
    if (daq_config_file.empty() || !load_daq_config(daq_config_file, daq_cfg)) {
        std::cerr << "Error: failed to load DAQ config"
                  << (daq_config_file.empty() ? " (not found)" : ": " + daq_config_file)
                  << "\n";
        return 1;
    }
    std::cerr << "DAQ cfg  : " << daq_config_file
              << " (adc_format=" << daq_cfg.adc_format << ")\n";

    const auto roc_to_crate = daq_cfg.roc_crate_map(true);

    EvChannel ch;
    ch.SetConfig(daq_cfg);
    if (ch.OpenAuto(evio_file) != status::success) {
        std::cerr << "Error: cannot open " << evio_file << "\n";
        return 1;
    }

    // accumulate: key = (crate, slot, channel)
    std::map<uint64_t, ChannelAccum> accum;

    auto event_ptr = std::make_unique<fdec::EventData>();
    auto &event = *event_ptr;
    int total = 0, selected = 0;

    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != EventType::Physics) continue;

        for (int i = 0; i < ch.GetNEvents(); ++i) {
            if (!ch.DecodeEvent(i, event)) continue;
            total++;

            if (!(event.info.trigger_bits & trigger_mask)) continue;
            selected++;

            // accumulate raw ADC values (no pedestal subtraction)
            for (int r = 0; r < event.nrocs; ++r) {
                auto &roc = event.rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                if (cit == roc_to_crate.end()) continue;
                int crate = cit->second;

                fdec::ForEachChannel(roc, [&](int s, int c, const fdec::ChannelData &cd) {
                    auto &acc = accum[prad2::pack_daq_key(crate, s, c)];
                    if (acc.count == 0) {
                        acc.crate   = crate;
                        acc.slot    = s;
                        acc.channel = c;
                    }
                    acc.add(cd.samples[0]);
                });
            }

            if (max_events > 0 && selected >= max_events) break;
        }
        if (max_events > 0 && selected >= max_events) break;

        if (total % 5000 == 0)
            std::cerr << "  " << total << " events scanned, "
                      << selected << " selected...\r" << std::flush;
    }
    ch.Close();

    std::cerr << "\nDone: " << total << " physics events, "
              << selected << " matched trigger bit " << trigger_bit
              << ", " << accum.size() << " channels\n";

    json out = json::array();
    for (auto &[k, acc] : accum) {
        out.push_back({
            {"crate",   acc.crate},
            {"slot",    acc.slot},
            {"channel", acc.channel},
            {"mean",    std::round(acc.mean() * 1000.) / 1000.},
            {"rms",     std::round(acc.rms() * 10000.) / 10000.},
            {"count",   acc.count},
        });
    }

    std::ofstream of(output_file);
    if (!of.is_open()) {
        std::cerr << "Error: cannot write " << output_file << "\n";
        return 1;
    }
    of << out.dump(2) << "\n";
    of.close();

    std::cerr << "Written: " << output_file << " (" << out.size() << " channels)\n";
    return 0;
}
