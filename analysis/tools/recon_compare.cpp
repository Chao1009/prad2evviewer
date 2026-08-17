//=============================================================================
// recon_compare -- compare event-level reconstruction counts from two recon
// ROOT files, matched by event_num.
//=============================================================================

#include "EventData.h"
#include "EventData_io.h"

#include <TFile.h>
#include <TTree.h>

#include <algorithm>
#include <cstdlib>
#include <iostream>
#include <string>
#include <unordered_map>

namespace {

struct EventSummary {
    int n_clusters = 0;
    int n_gem_hits = 0;
    int match_num = 0;
    uint32_t match_flags[prad2::kMaxClusters] = {};
};

using EventMap = std::unordered_map<int, EventSummary>;

bool read_recon(const std::string &path, EventMap &events)
{
    TFile file(path.c_str(), "READ");
    if (file.IsZombie()) {
        std::cerr << "recon_compare: cannot open " << path << "\n";
        return false;
    }

    auto *tree = dynamic_cast<TTree *>(file.Get("recon"));
    if (!tree) {
        std::cerr << "recon_compare: " << path << " has no 'recon' tree\n";
        return false;
    }

    prad2::ReconEventData event;
    prad2::SetReconReadBranches(tree, event);
    const auto entries = tree->GetEntries();
    events.reserve(static_cast<size_t>(entries));

    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        EventSummary summary;
        summary.n_clusters = event.n_clusters;
        summary.n_gem_hits = event.n_gem_hits;
        summary.match_num = event.matchNum;
        const int n_clusters = std::clamp(event.n_clusters, 0, prad2::kMaxClusters);
        std::copy_n(event.matchFlag, n_clusters, summary.match_flags);
        events[event.event_num] = summary;
    }
    return true;
}

bool same_match_flags(const EventSummary &left, const EventSummary &right)
{
    const int count = std::min({left.n_clusters, right.n_clusters,
                                prad2::kMaxClusters});
    for (int index = 0; index < count; ++index) {
        if (left.match_flags[index] != right.match_flags[index]) return false;
    }
    return true;
}

void print_usage(const char *program)
{
    std::cerr << "Usage: " << program
              << " <evio_recon.root> <raw_recon.root> [-n max_differences]\n";
}

} // namespace

int main(int argc, char *argv[])
{
    if (argc < 3) {
        print_usage(argv[0]);
        return 1;
    }

    int max_differences = 20;
    if (argc == 5 && std::string(argv[3]) == "-n") {
        max_differences = std::max(0, std::atoi(argv[4]));
    } else if (argc != 3) {
        print_usage(argv[0]);
        return 1;
    }

    EventMap evio_events;
    EventMap raw_events;
    if (!read_recon(argv[1], evio_events) || !read_recon(argv[2], raw_events))
        return 1;

    int missing_from_raw = 0;
    int missing_from_evio = 0;
    int same = 0;
    int gem_difference = 0;
    int match_difference = 0;
    int hycal_difference = 0;
    int printed = 0;

    for (const auto &[event_num, evio] : evio_events) {
        auto raw_it = raw_events.find(event_num);
        if (raw_it == raw_events.end()) {
            ++missing_from_raw;
            continue;
        }

        const auto &raw = raw_it->second;
        const bool same_hycal = evio.n_clusters == raw.n_clusters;
        const bool same_gem = evio.n_gem_hits == raw.n_gem_hits;
        const bool same_match = evio.match_num == raw.match_num
                                && same_match_flags(evio, raw);
        if (same_hycal && same_gem && same_match) {
            ++same;
            continue;
        }
        if (!same_hycal) ++hycal_difference;
        if (!same_gem) ++gem_difference;
        if (!same_match) ++match_difference;

        if (printed++ < max_differences) {
            std::cout << "event " << event_num
                      << ": evio(cl=" << evio.n_clusters
                      << ", gem=" << evio.n_gem_hits
                      << ", match=" << evio.match_num
                      << ") raw(cl=" << raw.n_clusters
                      << ", gem=" << raw.n_gem_hits
                      << ", match=" << raw.match_num << ")\n";
        }
    }

    for (const auto &[event_num, raw] : raw_events) {
        (void)raw;
        if (evio_events.find(event_num) == evio_events.end())
            ++missing_from_evio;
    }

    std::cout << "EVIO events: " << evio_events.size()
              << "  raw events: " << raw_events.size() << "\n"
              << "identical: " << same
              << "  missing from raw: " << missing_from_raw
              << "  missing from EVIO: " << missing_from_evio << "\n"
              << "HyCal differences: " << hycal_difference
              << "  2D GEM-hit differences: " << gem_difference
              << "  matching differences: " << match_difference << "\n";
    return 0;
}
