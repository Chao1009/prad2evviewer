// dsc_scan — explore DSC2 scaler banks (0xE115) in EVIO data to identify
// the right (crate, slot, channel) for live-time extraction.
//
// For every 0xE115 bank found, walks up to its parent ROC bank to identify
// which crate it belongs to, parses the per-slot 67-word DSC2 layout, and
// reports per-channel gated/ungated counts and the implied live time.
//
// On a typical PRad-II run only a single DSC2 module is read out, but the
// physics-trigger and reference-clock channels both offer a livetime.  This
// tool prints both so the user can pick.
//
// Usage: dsc_scan <input> [-D daq_config.json] [-N n_events] [--all]

#include "EvChannel.h"
#include "DaqConfig.h"
#include "load_daq_config.h"

#include <iostream>
#include <iomanip>
#include <string>
#include <cstring>
#include <cstdlib>
#include <vector>
#include <map>
#include <algorithm>
#include <getopt.h>

#ifndef _WIN32
#include <dirent.h>
#include <sys/stat.h>
#endif

using namespace evc;

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

// ---- DSC2 layout ----------------------------------------------------------
// In PRad-II run 024246, the 0xE115 bank is wrapped by a JLab-style
// SSP/VTP block: 3 prefix words (BLKHDR, EVTHDR, TRGTIME-like), then the
// 67-word DSC2 scaler payload (no 0xDCA00000 header word — that appears to
// be replaced by the BLKHDR), then 2 trailer words (FILLER, BLKTLR).
//
// 67-word DSC2 payload, per slot:
//   [0]      placeholder/header (was 0xDCA00000|(slot<<8)|rflag)
//   [1..16]  TRG  Grp1 (gated/busy)  — 16 channels
//   [17..32] TDC  Grp1 (gated/busy)  — 16 channels
//   [33..48] TRG  Grp2 (ungated)     — 16 channels
//   [49..64] TDC  Grp2 (ungated)     — 16 channels
//   [65]     Ref  Grp1 (gated/busy)  — 125 MHz clock
//   [66]     Ref  Grp2 (ungated)     — 125 MHz clock
//
// We accept either layout (with or without the 3-word prefix) by trying
// offsets 0 and 3, picking whichever places the ref words at sensible
// values (ref_ungated > ref_gated, both non-zero).
static constexpr uint32_t DSC2_BANK_TAG    = 0xE115;
static constexpr int      DSC2_NCH         = 16;
static constexpr double   DSC2_REF_FREQ    = 125.0e6;
static constexpr int      DSC2_PAYLOAD_W   = 67;

struct Dsc2Slot {
    uint32_t slot{0};
    uint32_t trg_gated[DSC2_NCH]{};
    uint32_t tdc_gated[DSC2_NCH]{};
    uint32_t trg_ungated[DSC2_NCH]{};
    uint32_t tdc_ungated[DSC2_NCH]{};
    uint32_t ref_gated{0};
    uint32_t ref_ungated{0};
    int      offset{0};   // payload offset inside the bank
};

static bool fill_slot_at(const uint32_t *data, size_t nwords, size_t off, Dsc2Slot &s)
{
    if (off + DSC2_PAYLOAD_W > nwords) return false;
    const uint32_t *p = &data[off + 1];           // skip the [0] header word
    std::memcpy(s.trg_gated,   p,      DSC2_NCH * 4);
    std::memcpy(s.tdc_gated,   p + 16, DSC2_NCH * 4);
    std::memcpy(s.trg_ungated, p + 32, DSC2_NCH * 4);
    std::memcpy(s.tdc_ungated, p + 48, DSC2_NCH * 4);
    s.ref_gated   = p[64];
    s.ref_ungated = p[65];
    s.offset      = (int)off;

    // Plausibility: ref_ungated must be non-zero and ≥ ref_gated.
    return s.ref_ungated > 0 && s.ref_ungated >= s.ref_gated;
}

static std::vector<Dsc2Slot> parse_dsc2(const uint32_t *data, size_t nwords)
{
    std::vector<Dsc2Slot> out;

    // Probe candidate offsets (header layouts).  The DSC2-firmware "rflag=0xFF"
    // legacy form puts a 0xDCA0 magic at offset 0; the "rflag=1" form (used in
    // PRad-II run 024246) wraps the payload in a JLab-style block: BLKHDR,
    // EVTHDR, TRGTIME, then payload[0]=placeholder, payload[1..16]=TRG_g, ...
    // i.e. the payload's "header" word lands at bank index 2.
    static const size_t kProbeOffsets[] = {0, 2, 3, 5};

    Dsc2Slot s{};
    for (size_t off : kProbeOffsets) {
        if (off + DSC2_PAYLOAD_W > nwords) continue;
        uint32_t hdr = data[off];
        if ((hdr & 0xFFFF0000u) == 0xDCA00000u) {
            // legacy DSC2 magic header
            if (!fill_slot_at(data, nwords, off, s)) continue;
            s.slot = (hdr >> 8) & 0xFF;
            out.push_back(s);
            return out;
        }
        if (off >= 1 && nwords >= 1 && (data[0] >> 27) == 0x10) {
            // JLab BLKHDR-wrapped form — slot is in BLKHDR bits 26:22.
            if (!fill_slot_at(data, nwords, off, s)) continue;
            s.slot = (data[0] >> 22) & 0x1F;
            out.push_back(s);
            return out;
        }
    }
    return out;
}

// ---- file discovery (same as livetime.cpp) ---------------------------------

static bool is_regular_file(const std::string &p)
{ struct stat st; return stat(p.c_str(), &st) == 0 && S_ISREG(st.st_mode); }
static bool is_directory(const std::string &p)
{ struct stat st; return stat(p.c_str(), &st) == 0 && S_ISDIR(st.st_mode); }

static std::vector<std::string> list_dir(const std::string &dir)
{
    std::vector<std::string> e;
    DIR *d = opendir(dir.c_str()); if (!d) return e;
    while (auto *en = readdir(d)) {
        if (en->d_name[0] == '.') continue;
        std::string f = dir + "/" + en->d_name;
        if (is_regular_file(f)) e.push_back(f);
    }
    closedir(d);
    return e;
}

static std::vector<std::string> discover(const std::string &path)
{
    std::vector<std::string> files;
    if (is_regular_file(path)) {
        auto dot = path.rfind('.');
        if (dot != std::string::npos) {
            std::string suf = path.substr(dot + 1);
            bool split = !suf.empty() && suf.find_first_not_of("0123456789") == std::string::npos;
            if (split) {
                std::string base = path.substr(0, dot);
                auto sl = base.rfind('/');
                std::string dir = sl != std::string::npos ? base.substr(0, sl) : ".";
                std::string bn  = sl != std::string::npos ? base.substr(sl + 1) : base;
                for (auto &f : list_dir(dir)) {
                    auto fsl = f.rfind('/');
                    std::string fn = fsl != std::string::npos ? f.substr(fsl + 1) : f;
                    if (fn.size() > bn.size() && fn.substr(0, bn.size() + 1) == bn + ".")
                        files.push_back(f);
                }
                std::sort(files.begin(), files.end());
                if (!files.empty()) return files;
            }
        }
        files.push_back(path); return files;
    }
    if (is_directory(path)) {
        for (auto &f : list_dir(path))
            if (f.find(".evio") != std::string::npos) files.push_back(f);
        std::sort(files.begin(), files.end());
    }
    return files;
}

// ---- main ------------------------------------------------------------------

int main(int argc, char *argv[])
{
    std::string input, dcfg;
    int n_events_max = 0;          // 0 = unlimited
    bool dump_all = false;         // print every DSC2 sighting (not just first/last)

    std::string db = DATABASE_DIR;
    if (auto *e = std::getenv("PRAD2_DATABASE_DIR")) db = e;
    dcfg = db + "/daq_config.json";

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
            std::cerr << "Usage: " << argv[0]
                      << " <input> [-D daq_config.json] [-N max_events] [--all]\n";
            return opt == 'h' ? 0 : 1;
        }
    }
    if (optind < argc) input = argv[optind];
    if (input.empty()) {
        std::cerr << "Usage: " << argv[0] << " <input> [-D daq_config.json] [-N max_events] [--all]\n";
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

    auto files = discover(input);
    if (files.empty()) {
        std::cerr << "No EVIO files found for: " << input << "\n";
        return 1;
    }
    std::cerr << "Scanning " << files.size() << " file(s):\n";
    for (auto &f : files) std::cerr << "  " << f << "\n";

    EvChannel ch;
    ch.SetConfig(cfg);

    // Per (parent_tag, slot) keep the most-recent Dsc2Slot we saw and a count
    // of how many SYNC/Physics events carried it.
    struct Latest { Dsc2Slot last{}; uint64_t hits{0}; uint32_t event_tag{0}; };
    std::map<std::pair<uint32_t,int>, Latest> seen;

    // For the very first sighting we also remember the values, to compute
    // delta = last - first → trigger live time over the whole scanned span.
    std::map<std::pair<uint32_t,int>, Dsc2Slot> first_seen;

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

            // record EPICS sightings too — useful for slow-tree design.
            if (et == EventType::Epics) {
                static int epics_dumped = 0;
                if (epics_dumped < 3) {
                    ++epics_dumped;
                    const auto &si = ch.Sync();
                    std::cerr << "EPICS ev#" << scanned
                              << "  top tag=0x" << std::hex << ch.GetEvHeader().tag
                              << "  sync.unix_time=" << std::dec << si.unix_time
                              << "  sync.sync_counter=" << si.sync_counter
                              << "  sync.run_number="  << si.run_number
                              << "  sync.event_tag=0x" << std::hex << si.event_tag
                              << std::dec << "\n";
                }
            }

            if (et != EventType::Sync && et != EventType::Physics) continue;

            // O(1) lookup of every 0xE115 node in this event.
            const auto &idxs = ch.NodesForTag(DSC2_BANK_TAG);
            if (idxs.empty()) continue;
            ++with_dsc;
            if (et == EventType::Sync) ++sync_count;

            // Per-DSC2-event timestamp dump (first 3 sightings).
            static int dsc_dumped = 0;
            if (dsc_dumped < 3) {
                ++dsc_dumped;
                const auto &si = ch.Sync();
                ch.SelectEvent(0);
                const auto &info = ch.Info();
                std::cerr << "DSC ev#" << scanned
                          << "  top tag=0x" << std::hex << ch.GetEvHeader().tag
                          << std::dec
                          << "  ev_num=" << info.event_number
                          << "  TI_ticks=" << info.timestamp
                          << "  sync.unix=" << si.unix_time
                          << "  sync.cnt=" << si.sync_counter
                          << "  sync.run=" << si.run_number
                          << "\n";
            }

            const auto &nodes = ch.GetNodes();
            for (int idx : idxs) {
                const EvNode &node = nodes[idx];
                if (node.data_words == 0) continue;
                uint32_t parent_tag = 0;
                if (node.parent >= 0 && node.parent < (int)nodes.size())
                    parent_tag = nodes[node.parent].tag;

                auto slots = parse_dsc2(ch.GetData(node), node.data_words);
                if (slots.empty() && with_dsc <= 2) {
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
                for (auto &s : slots) {
                    auto key = std::make_pair(parent_tag, (int)s.slot);
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

        const Dsc2Slot &s  = L.last;
        const Dsc2Slot &s0 = first_seen[key];

        std::cout << "\n--- DSC2 module @ ROC tag 0x" << std::hex << std::setw(4)
                  << std::setfill('0') << parent << std::setfill(' ') << std::dec
                  << "  slot " << slot << "  (payload offset=" << s.offset << ")";
        if (rn != roc_name.end())  std::cout << "  (" << rn->second << ")";
        if (rc != roc_crate.end()) std::cout << "  crate=" << rc->second;
        std::cout << "  hits=" << L.hits << "\n";

        // cumulative (since first sighting in scan) ratio
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
            for (int c = 0; c < DSC2_NCH; ++c)
                if (ungated[c] != 0) { any = true; break; }
            if (!any) return;

            std::cout << "\n  " << label << " channel scaler counts (cumulative since run start):\n"
                      << "   ch  |  gated         ungated         g/u(%)    1-g/u(%) | "
                      << "scan: dg          du              g/u(%)    1-g/u(%)\n";
            for (int c = 0; c < DSC2_NCH; ++c) {
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

    // Recommendation
    std::cout << "\n=== Recommendation ===\n"
        << "  • The DSC2 scaler bank lives at parent ROC tag 0x0027 (TI master) — fixed for this DAQ.\n"
        << "  • The slot value comes from the JLab BLKHDR (bits 26:22) printed above; record it as is.\n"
        << "  • Live-time convention: in PRad-II run 024246 the (gated, ungated) ratios are ~0.99,\n"
        << "    which means gated counts LIVE time (gate enabled while NOT busy).  Use the\n"
        << "    formula  live = gated/ungated, NOT 1 - gated/ungated.  If you keep the existing\n"
        << "    livetime.cpp 1-g/u formula, the reported live time will be the dead-time fraction.\n"
        << "  • For per-trigger live time, pick a TRG channel whose ungated count matches the\n"
        << "    expected trigger rate (column 'du' / scan duration).  Channel 2 looks active here.\n"
        << "  • Update database/daq_config.json:\n"
        << "      \"dsc_scaler\": { \"bank_tag\": \"0xE115\", \"slot\": <slot>, \"source\": \"ref|trg|tdc\", \"channel\": <c> }\n"
        << "    and verify that AppState::processDscBank uses the live-time convention you want.\n";
    return 0;
}
