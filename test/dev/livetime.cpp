// livetime — calculate DAQ live time from EVIO files
//
// Two independent methods:
//   1. DSC2 scalers: gated/ungated reference counts from the 0xe115 bank
//      (gated counts while NOT busy; see DscData.h).
//   2. Pulser counting: accepted 100 Hz pulser events vs expected from
//      elapsed time.
//
// Accepts a single file, a base name (auto-discovers .00000, .00001, ...),
// or a directory (processes all .evio* files sorted by name).

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "Fadc250Data.h"
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

static const char *tag_name(uint32_t tag)
{
    switch (tag) {
    case 0x00A9: return "SSP RawSum";
    case 0x00B0: return "Pulser 100Hz";
    case 0x00B9: return "LMS";
    case 0x00BA: return "Alpha";
    case 0x00BC: return "Master OR";
    case 0x00FA: return "Cluster";
    default:     return nullptr;
    }
}

static void usage(const char *prog)
{
    std::cerr << "Usage: " << prog
              << " <input> [-D daq_config.json] [-f freq_hz] [-t interval_sec]\n"
              << "  <input>: file, base name (finds .00000 .00001 ...), or directory\n";
}

int main(int argc, char *argv[])
{
    std::string input;
    std::string daq_config_file = prad2::database_dir() + "/daq_config.json";
    double pulser_freq = 100.0;   // Hz
    double report_interval = 10.0; // seconds

    int opt;
    while ((opt = getopt(argc, argv, "D:f:t:h")) != -1) {
        switch (opt) {
        case 'D': daq_config_file = optarg; break;
        case 'f': pulser_freq = std::atof(optarg); break;
        case 't': report_interval = std::atof(optarg); break;
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
    if (!load_daq_config(daq_config_file, cfg)) {
        std::cerr << "Failed to load DAQ config: " << daq_config_file << "\n";
        return 1;
    }

    auto files = devtools::discover_evio_inputs(input);
    if (files.empty()) {
        std::cerr << "No EVIO files found for: " << input << "\n";
        return 1;
    }
    std::cerr << "Processing " << files.size() << " file(s):\n";
    for (auto &f : files) std::cerr << "  " << f << "\n";

    // ---- state ----
    static constexpr uint32_t PULSER_TAG = 0x00B0;
    static constexpr uint32_t DSC2_BANK_TAG = 0xE115;

    auto event = std::make_unique<fdec::EventData>();
    std::map<uint32_t, uint64_t> tag_counts;
    uint64_t first_ts = 0, last_ts = 0;
    uint64_t total_physics = 0;
    uint32_t run_number = 0;
    uint32_t unix_start = 0, unix_end = 0;
    int sync_count = 0;

    // DSC2 cumulative (latest SYNC)
    dsc::DscEventData dsc2_data;

    // periodic reporting state
    double next_report = report_interval;

    auto print_header = []() {
        std::cout << std::left
                  << std::setw(10) << "Time(s)"
                  << std::setw(12) << "LT_DSC2(%)"
                  << std::setw(14) << "LT_Pulser(%)"
                  << std::setw(10) << "Physics"
                  << std::setw(10) << "Pulser"
                  << std::setw(8)  << "Syncs"
                  << std::setw(0)  << "File"
                  << "\n" << std::string(78, '-') << "\n";
    };

    std::string current_file;
    auto print_row = [&](double elapsed) {
        // cumulative DSC2 live time (latest SYNC), reference pair gated/ungated
        double lt_dsc2 = -1;
        if (dsc2_data.present && dsc2_data.ref_ungated > 0) {
            lt_dsc2 = static_cast<double>(dsc2_data.ref_gated)
                    / dsc2_data.ref_ungated * 100.0;
        }
        // cumulative pulser live time
        double lt_pulser = -1;
        uint64_t pulser_total = 0;
        auto it = tag_counts.find(PULSER_TAG);
        if (it != tag_counts.end()) pulser_total = it->second;
        double expected = pulser_freq * elapsed;
        if (expected > 0)
            lt_pulser = 100.0 * static_cast<double>(pulser_total) / expected;

        auto slash = current_file.rfind('/');
        if (slash == std::string::npos) slash = current_file.rfind('\\');
        std::string fname = (slash != std::string::npos) ? current_file.substr(slash + 1) : current_file;

        std::cout << std::fixed
                  << std::setw(10) << std::setprecision(1) << elapsed;
        if (lt_dsc2 >= 0)
            std::cout << std::setw(12) << std::setprecision(2) << lt_dsc2;
        else
            std::cout << std::setw(12) << "--";
        if (lt_pulser >= 0)
            std::cout << std::setw(14) << std::setprecision(2) << lt_pulser;
        else
            std::cout << std::setw(14) << "--";
        std::cout << std::setw(10) << total_physics
                  << std::setw(10) << pulser_total
                  << std::setw(8)  << sync_count
                  << fname
                  << "\n";
    };

    print_header();

    // ---- process files ----
    EvChannel ch;
    ch.SetConfig(cfg);

    for (auto &file : files) {
        current_file = file;
        if (ch.OpenAuto(file) != status::success) {
            std::cerr << "Warning: cannot open " << file << ", skipping\n";
            continue;
        }

        while (ch.Read() == status::success) {
            if (!ch.Scan()) continue;
            auto evtype = ch.GetEventType();

            if (evtype == EventType::Prestart || evtype == EventType::Go) {
                uint32_t ct = ch.Sync().unix_time;
                if (ct != 0 && unix_start == 0) unix_start = ct;
            }
            if (evtype == EventType::End) {
                uint32_t ct = ch.Sync().unix_time;
                if (ct != 0) unix_end = ct;
            }

            // DSC2 scalers
            if (evtype == EventType::Sync || evtype == EventType::Physics) {
                const EvNode *dsc2_node = ch.FindFirstByTag(DSC2_BANK_TAG);
                dsc::DscEventData parsed;
                if (dsc2_node && dsc::Dsc2Decoder::ParsePayload(ch.GetData(*dsc2_node),
                                                               dsc2_node->data_words, parsed)) {
                    dsc2_data = parsed;
                    if (evtype == EventType::Sync) sync_count++;
                }
            }

            if (evtype != EventType::Physics) continue;

            for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
                event->clear();
                if (!ch.DecodeEvent(ie, *event)) continue;

                tag_counts[event->info.event_tag]++;
                total_physics++;

                if (run_number == 0 && event->info.run_number != 0)
                    run_number = event->info.run_number;

                uint64_t ts = event->info.timestamp;
                if (ts != 0) {
                    if (first_ts == 0) first_ts = ts;
                    last_ts = ts;

                    double elapsed = static_cast<double>(ts - first_ts) * fdec::TI_TICK_SEC;
                    while (elapsed >= next_report) {
                        print_row(next_report);
                        next_report += report_interval;
                    }
                }
            }
        }
        ch.Close();
    }

    // ---- final report ----
    double elapsed_ti = (first_ts != 0 && last_ts > first_ts)
        ? static_cast<double>(last_ts - first_ts) * fdec::TI_TICK_SEC : 0.0;
    double elapsed_unix = (unix_start != 0 && unix_end > unix_start)
        ? static_cast<double>(unix_end - unix_start) : 0.0;
    double elapsed = (elapsed_ti > 0) ? elapsed_ti : elapsed_unix;

    // final row (partial interval)
    if (elapsed > 0)
        print_row(elapsed);

    // ---- summary ----
    std::cout << "\n=== Summary ===\n";
    if (run_number != 0)
        std::cout << "Run number     : " << run_number << "\n";
    std::cout << "Files          : " << files.size() << "\n";
    std::cout << "Total physics  : " << total_physics << "\n";
    std::cout << "SYNC events    : " << sync_count << "\n";
    if (elapsed_ti > 0)
        std::cout << std::fixed << std::setprecision(2)
                  << "Elapsed (TI)   : " << elapsed_ti << " sec\n";
    if (elapsed_unix > 0)
        std::cout << std::fixed << std::setprecision(0)
                  << "Elapsed (unix) : " << elapsed_unix << " sec\n";

    std::cout << "\n--- Trigger counts ---\n";
    std::cout << std::left << std::setw(10) << "Tag"
              << std::setw(20) << "Name"
              << std::right << std::setw(10) << "Count" << "\n";
    std::cout << std::string(40, '-') << "\n";
    for (auto &[tag, count] : tag_counts) {
        const char *name = tag_name(tag);
        char hex[16];
        snprintf(hex, sizeof(hex), "0x%04X", tag);
        std::cout << std::left << std::setw(10) << hex
                  << std::setw(20) << (name ? name : "unknown")
                  << std::right << std::setw(10) << count << "\n";
    }

    // DSC2 summary
    std::cout << "\n--- DSC2 scaler live time (cumulative) ---\n";
    if (!dsc2_data.present) {
        std::cout << "(no DSC2 scaler bank 0xE115 found)\n";
    } else {
        const auto &s = dsc2_data;
        double lt = (s.ref_ungated > 0)
            ? static_cast<double>(s.ref_gated) / s.ref_ungated * 100.0
            : 0.0;
        std::cout << std::fixed << std::setprecision(2);
        std::cout << "  DSC2 slot " << s.slot
                  << ": ref_gated=" << s.ref_gated
                  << "  ref_ungated=" << s.ref_ungated
                  << "  live=" << std::setprecision(3) << lt << "%\n";

        bool any = false;
        for (int c = 0; c < dsc::DSC2_NCH; ++c) {
            if (s.trg_ungated[c] == 0) continue;
            if (!any) {
                std::cout << "    TRG ch  gated(live)  ungated(total)  live%\n";
                any = true;
            }
            double cl = static_cast<double>(s.trg_gated[c])
                      / s.trg_ungated[c] * 100.0;
            std::cout << "      " << std::setw(2) << c
                      << std::setw(13) << s.trg_gated[c]
                      << std::setw(16) << s.trg_ungated[c]
                      << std::setw(9) << std::setprecision(2) << cl << "\n";
        }
    }

    // pulser summary
    uint64_t pulser_count = 0;
    auto pit = tag_counts.find(PULSER_TAG);
    if (pit != tag_counts.end()) pulser_count = pit->second;
    double expected_pulser = pulser_freq * elapsed;
    double lt_pulser = (expected_pulser > 0)
        ? 100.0 * static_cast<double>(pulser_count) / expected_pulser : 0.0;

    std::cout << "\n--- Pulser counting live time ---\n";
    std::cout << std::fixed << std::setprecision(2);
    std::cout << "  Pulser freq    : " << pulser_freq << " Hz\n";
    std::cout << "  Pulser accepted: " << pulser_count << "\n";
    if (expected_pulser > 0)
        std::cout << "  Pulser expected: " << std::setprecision(0) << expected_pulser << "\n";
    std::cout << std::setprecision(2);
    if (lt_pulser > 0)
        std::cout << "  Live time      : " << lt_pulser << " %\n";
    else
        std::cout << "  Live time      : N/A\n";

    return 0;
}
