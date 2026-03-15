// src/test_main.cpp
// Basic smoke-test for the evc library.
// Usage:
//   evc_test <evio_file>            -- read from file
//   evc_test --et <ip> <port> <file> <station>  -- read from ET system

#include "EvChannel.h"
#include "EtChannel.h"
#include "Fadc250Decoder.h"
#include <iostream>
#include <string>
#include <cstdlib>

// FADC250 composite bank tag
static constexpr uint32_t TAG_FADC250_COMPOSITE = 0xe126;

static void usage(const char *prog)
{
    std::cerr << "Usage:\n"
              << "  " << prog << " <evio_file>\n"
              << "  " << prog << " --et <ip> <port> <et_file> <station>\n";
}

// Decode and print composite banks from the current event
static void decodeComposites(evc::EvChannel &ch, const fdec::Fadc250Decoder &decoder, int event_num)
{
    auto &composites = ch.GetCompositeInfos();
    if (composites.empty()) return;

    for (auto &ci : composites) {
        size_t nbytes = 0;
        const uint8_t *data = ch.GetCompositeData(ci, nbytes);

        auto hits = decoder.DecodeComposite(data, nbytes);

        std::cout << "  Event " << event_num
                  << "  ROC=" << ci.roc
                  << "  bank=0x" << std::hex << ci.tag << std::dec
                  << "  composite_hits=" << hits.size() << "\n";

        for (auto &hit : hits) {
            std::cout << "    slot=" << (int)hit.slot
                      << "  ch=" << (int)hit.channel
                      << "  nSamples=" << hit.samples.size();
            // print first few samples as a preview
            if (!hit.samples.empty()) {
                std::cout << "  samples=[";
                size_t nshow = std::min<size_t>(hit.samples.size(), 6);
                for (size_t i = 0; i < nshow; ++i) {
                    if (i) std::cout << ",";
                    std::cout << hit.samples[i];
                }
                if (hit.samples.size() > nshow) std::cout << ",...";
                std::cout << "]";
            }
            std::cout << "\n";
        }
    }
}

// ---- file mode -----------------------------------------------------------
static int testFile(const std::string &path)
{
    evc::EvChannel ch;
    fdec::Fadc250Decoder decoder(250.0);

    if (ch.Open(path) != evc::status::success) {
        std::cerr << "Failed to open: " << path << "\n";
        return 1;
    }

    // banks to scan: include the FADC250 composite tag
    std::vector<uint32_t> banks_of_interest = { TAG_FADC250_COMPOSITE };

    int nevents = 0;
    int ncomposite_events = 0;
    evc::status st;
    while ((st = ch.Read()) == evc::status::success) {
        ++nevents;
        auto hdr = ch.GetEvHeader();

        std::cout << "Event " << nevents
                  << "  tag=" << hdr.tag
                  << "  type=" << hdr.type
                  << "  length=" << hdr.length;

        // Try to scan the event for data banks
        // Physics events typically have tag >= 0xFF00 range, or specific detector tags.
        // We try scanning all events - non-physics events will just fail gracefully.
        bool scanned = false;
        try {
            scanned = ch.ScanBanks(banks_of_interest);
        } catch (...) {
            // Not a scannable event (e.g. control event), skip
        }

        if (scanned && !ch.GetCompositeInfos().empty()) {
            ++ncomposite_events;
            std::cout << "  [composite: " << ch.GetCompositeInfos().size() << " bank(s)]";
        }
        std::cout << "\n";

        // Decode composite data for this event
        if (scanned) {
            decodeComposites(ch, decoder, nevents);
        }
    }

    std::cout << "\nDone. Read " << nevents << " event(s), "
              << ncomposite_events << " with composite bank(s). "
              << "Final status: " << static_cast<int>(st) << "\n";
    ch.Close();
    return 0;
}

// ---- ET mode -------------------------------------------------------------
static int testET(const std::string &ip, int port,
                  const std::string &et_file, const std::string &station)
{
    evc::EtChannel ch;
    fdec::Fadc250Decoder decoder(250.0);

    if (ch.Connect(ip, port, et_file) != evc::status::success) {
        std::cerr << "Failed to connect to ET at " << ip << ":" << port << "\n";
        return 1;
    }
    if (ch.Open(station) != evc::status::success) {
        std::cerr << "Failed to open station: " << station << "\n";
        ch.Disconnect();
        return 1;
    }

    std::vector<uint32_t> banks_of_interest = { TAG_FADC250_COMPOSITE };

    int nevents = 0, max_events = 20;
    evc::status st;
    while (nevents < max_events && (st = ch.Read()) != evc::status::failure) {
        if (st == evc::status::empty) continue;
        ++nevents;
        auto hdr = ch.GetEvHeader();
        std::cout << "ET Event " << nevents
                  << "  tag=" << hdr.tag
                  << "  length=" << hdr.length << "\n";

        bool scanned = false;
        try {
            scanned = ch.ScanBanks(banks_of_interest);
        } catch (...) {}

        if (scanned) {
            decodeComposites(ch, decoder, nevents);
        }
    }

    std::cout << "Done. Read " << nevents << " event(s).\n";
    ch.Disconnect();
    return 0;
}

// --------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    if (argc < 2) { usage(argv[0]); return 1; }

    std::string first = argv[1];

    if (first == "--et") {
        if (argc < 6) { usage(argv[0]); return 1; }
        return testET(argv[2], std::atoi(argv[3]), argv[4], argv[5]);
    }

    return testFile(first);
}
