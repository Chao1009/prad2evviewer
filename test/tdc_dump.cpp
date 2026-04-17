// tdc_dump — extract V1190 TDC hits (bank 0xE107) from EVIO files
//
// Usage:
//   tdc_dump <input.evio> [-o out.csv] [-b out.bin] [-n max_events]
//            [-D daq_config.json] [--roc 0x008E]
//
// Outputs one record per TDC hit.  Two output formats:
//   CSV (default, to stdout or -o file):
//     event_num,trigger_bits,roc_tag,slot,channel,edge,tdc
//   Binary (-b out.bin):
//     16-byte magic header "PRAD2_TDC_HITS_1"
//     uint32 record_count
//     then N records of 16 bytes each:
//       uint32 event_num, uint32 trigger_bits,
//       uint16 roc_tag,   uint8 slot, uint8 channel_edge, uint32 tdc
//     channel_edge packs:  bit 7 = edge, bits 6:0 = channel
//
// The binary format is what scripts/tdc_viewer.py expects; CSV is provided
// for ad-hoc inspection.

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"
#include "TdcData.h"

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <cstring>
#include <getopt.h>

using namespace evc;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace {

struct BinHit {
    uint32_t event_num;
    uint32_t trigger_bits;
    uint16_t roc_tag;
    uint8_t  slot;
    uint8_t  channel_edge;   // bit 7 = edge, bits 6:0 = channel
    uint32_t tdc;
};
static_assert(sizeof(BinHit) == 16, "BinHit must be 16 bytes");

constexpr const char *BIN_MAGIC = "PRAD2_TDC_HITS_1";

void print_usage(const char *argv0)
{
    std::cerr <<
        "Usage: " << argv0 << " <input.evio> [options]\n"
        "  -o <file>       CSV output file (default: stdout)\n"
        "  -b <file>       Binary output file (for tdc_viewer.py)\n"
        "  -n <N>          Stop after N physics events (0 = all)\n"
        "  -D <json>       DAQ config file (default: $PRAD2_DATABASE_DIR/daq_config.json)\n"
        "  --roc 0xTAG     Only keep hits from parent ROC with this tag\n"
        "  -q              Quiet (no progress)\n"
        "  -h              Help\n";
}

} // namespace

int main(int argc, char *argv[])
{
    std::string input, daq_config_file;
    std::string csv_path, bin_path;
    int max_events = 0;
    int roc_filter = -1;
    bool quiet = false;

    std::string db_dir = DATABASE_DIR;
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR")) db_dir = env;
    daq_config_file = db_dir + "/daq_config.json";

    static struct option long_opts[] = {
        {"roc",    required_argument, 0, 'R'},
        {0, 0, 0, 0}
    };

    int opt;
    while ((opt = getopt_long(argc, argv, "o:b:n:D:qh", long_opts, nullptr)) != -1) {
        switch (opt) {
        case 'o': csv_path = optarg;                                   break;
        case 'b': bin_path = optarg;                                   break;
        case 'n': max_events = std::atoi(optarg);                      break;
        case 'D': daq_config_file = optarg;                            break;
        case 'R': roc_filter = std::stoi(optarg, nullptr, 0);          break;
        case 'q': quiet = true;                                        break;
        case 'h': print_usage(argv[0]); return 0;
        default:  print_usage(argv[0]); return 1;
        }
    }
    if (optind < argc) input = argv[optind];
    if (input.empty()) { print_usage(argv[0]); return 1; }

    DaqConfig cfg;
    if (!load_daq_config(daq_config_file, cfg)) {
        std::cerr << "Failed to load DAQ config: " << daq_config_file << "\n";
        return 1;
    }
    if (cfg.tdc_bank_tag == 0) {
        std::cerr << "Warning: daq_config.bank_tags.tdc is not set; "
                  << "no 0xE107 banks will be decoded.\n";
    }

    EvChannel ch;
    ch.SetConfig(cfg);
    if (ch.Open(input) != status::success) {
        std::cerr << "Cannot open: " << input << "\n";
        return 1;
    }

    // --- output setup ------------------------------------------------------
    std::ofstream csv_file;
    std::ostream *csv_out = &std::cout;
    if (!csv_path.empty()) {
        csv_file.open(csv_path);
        if (!csv_file) {
            std::cerr << "Cannot open CSV output: " << csv_path << "\n";
            return 1;
        }
        csv_out = &csv_file;
    }
    const bool write_csv = !bin_path.empty() ? !csv_path.empty() : true;

    std::ofstream bin_file;
    if (!bin_path.empty()) {
        bin_file.open(bin_path, std::ios::binary);
        if (!bin_file) {
            std::cerr << "Cannot open binary output: " << bin_path << "\n";
            return 1;
        }
        bin_file.write(BIN_MAGIC, 16);
        uint32_t placeholder = 0;
        bin_file.write(reinterpret_cast<const char*>(&placeholder), 4);
    }

    if (write_csv)
        (*csv_out) << "event_num,trigger_bits,roc_tag,slot,channel,edge,tdc\n";

    // --- decode loop -------------------------------------------------------
    auto event   = std::make_unique<fdec::EventData>();
    auto tdc_evt = std::make_unique<tdc::TdcEventData>();

    uint64_t n_events = 0;
    uint64_t n_hits   = 0;

    while (ch.Read() == status::success) {
        if (!ch.Scan()) continue;
        if (ch.GetEventType() != EventType::Physics) continue;

        int nsub = ch.GetNEvents();
        for (int i = 0; i < nsub; ++i) {
            event->clear();
            tdc_evt->clear();
            if (!ch.DecodeEvent(i, *event, nullptr, nullptr, tdc_evt.get())) {
                // Event info is still populated even if no FADC data; keep going.
            }

            const uint32_t evnum  = event->info.event_number;
            const uint32_t tbits  = event->info.trigger_bits;

            for (int h = 0; h < tdc_evt->n_hits; ++h) {
                const tdc::TdcHit &hit = tdc_evt->hits[h];
                if (roc_filter >= 0 && static_cast<int>(hit.roc_tag) != roc_filter)
                    continue;

                if (write_csv) {
                    (*csv_out) << evnum            << ','
                               << tbits            << ','
                               << hit.roc_tag      << ','
                               << int(hit.slot)    << ','
                               << int(hit.channel) << ','
                               << int(hit.edge)    << ','
                               << hit.value        << '\n';
                }
                if (bin_file.is_open()) {
                    BinHit b{};
                    b.event_num    = evnum;
                    b.trigger_bits = tbits;
                    b.roc_tag      = static_cast<uint16_t>(hit.roc_tag);
                    b.slot         = hit.slot;
                    b.channel_edge = static_cast<uint8_t>(
                        ((hit.edge & 0x1) << 7) | (hit.channel & 0x7F));
                    b.tdc          = hit.value;
                    bin_file.write(reinterpret_cast<const char*>(&b), sizeof(b));
                }
                ++n_hits;
            }

            ++n_events;
            if (!quiet && (n_events % 50000) == 0)
                std::cerr << "  processed " << n_events << " events, "
                          << n_hits << " hits\n";
            if (max_events > 0 && static_cast<int64_t>(n_events) >= max_events)
                goto done;
        }
    }

done:
    if (bin_file.is_open()) {
        uint32_t count = static_cast<uint32_t>(n_hits);
        bin_file.seekp(16, std::ios::beg);
        bin_file.write(reinterpret_cast<const char*>(&count), 4);
    }

    if (!quiet)
        std::cerr << "Done: " << n_events << " events, " << n_hits << " hits\n";

    return 0;
}
