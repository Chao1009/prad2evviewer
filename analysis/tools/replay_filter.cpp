//=============================================================================
// replay_filter.cpp — slow-control filter for replayed ROOT files
//
// Reads one or more replayed ROOT files (raw or recon), applies user-defined
// cuts on the slow streams (DSC2 livetime + EPICS values), and writes an
// output ROOT file containing:
//   * the events / recon tree with only the physics events bracketed by
//     two adjacent "good" slow-event checkpoints,
//   * the scalers and epics trees concatenated from every input file (no
//     filtering — they are small and useful as run-wide context),
// plus a JSON report with one entry per (cut-channel, slow-event) point so
// downstream tools can plot per-channel value traces with pass/fail status
// and the cut acceptance band.
//
// Cuts JSON schema:
//   {
//     "livetime": {
//       "source":  "ref",            // "ref" | "trg" | "tdc"
//       "channel": 0,                // ignored for "ref"
//       "abs":     { "min": 90, "max": 100 },
//       "rel_rms": 3
//     },
//     "epics": {
//       "<channel_name>": { "abs": {...}, "rel_rms": 3 },
//       ...
//     }
//   }
//
// `rel_rms: N` accepts points within N · σ̂ of the channel's median, where
// σ̂ = 1.4826 · MAD (median absolute deviation).  MAD is robust to heavy
// outliers (one bad reading does not pull the centre or width).
//
// Optional "split" block — classify each event by a slow-control PV level
// (e.g. target cell pressure) and write the "full" and "empty" subsets to
// separate output files.  Every checkpoint is labelled by its PV reading:
//   PV >= full  → full-target        PV <= empty → empty-target
//   empty < PV < full (the ramp), or no reading  → dropped (lands in neither)
// so events taken while the cell is filling/emptying are excluded from both.
//   "split": {
//     "channel": "TGT:PRad:Cell_P",  // EPICS channel to watch
//     "full":  500,                  // PV >= this  → full-target file
//     "empty": 5,                    // PV <= this  → empty-target file
//     "guard_checkpoints": 0,        // optional ± margin dropped at each state edge
//     "labels": ["full", "empty"]    // output suffixes (defaults: full / empty)
//   }
// Classification is per-checkpoint and stateless, so any number of full<->empty
// transitions in one run are handled, and a run that only ever shows one state
// produces only that one file (a pure full run → just <stem>_full.root).  Each
// side still gets every other configured cut and its own report + charge.
//
// Output ROOT file(s):
//   * events / recon — same schema as input, only kept events
//   * scalers / epics — concatenated from every input plus an extra
//     `good` boolean branch per row reflecting that checkpoint's
//     overall verdict (all cuts pass); with split on, `good` also requires
//     the row to belong to that file's side, so live_charge on a side file
//     reproduces that side's post-cut charge directly.
// With split on, one file <stem>_<label>.root + matching report is written per
// target state actually seen (one or two), instead of the single default output.
// With several inputs, -o is a directory: one filtered file per (input, side),
// prad_<run>_epics.root with the run's slow trees, and one run-level report.
// JSON report: see the write phase in the source for the full layout (a
// "split" block is added per side when run-splitting is active).
//=============================================================================

#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"     // analysis::get_run_int
#include "SlowControl.h"
#include "ToolUtils.h"

#include <TFile.h>
#include <TTree.h>

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <mutex>
#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include <getopt.h>

using json = nlohmann::json;

namespace {

// ── Cut configuration ────────────────────────────────────────────────────

struct AbsCut {
    bool   has_min = false;
    double min_val = 0;
    bool   has_max = false;
    double max_val = 0;
};

struct ChannelCut {
    AbsCut       abs;
    bool         has_rel_rms = false;
    double       rel_rms_n   = 0;
    // Optional: condition this channel's robust median/MAD on points where
    // every named gating channel's cut passed.  Use when the channel of
    // interest is physically meaningful only in a particular regime —
    // e.g. beam position is only meaningful when current is above some
    // floor.  Accepts either a single string or an array of strings in
    // the JSON; multiple gates are ANDed (a row's value is included in
    // the stats only if ALL listed gating channels passed).  All listed
    // channels must also be configured.  One level of gating is
    // supported (the gating channels themselves use ungated stats).
    // EPICS-on-EPICS only.
    std::vector<std::string> gated_by;

    // Robust statistics (filled in phase 2 if has_rel_rms).  `n_used` is the
    // input point count (MAD doesn't iterate-and-drop); `n_clipped` is points
    // outside [center − N·sigma, center + N·sigma], reported for traceability.
    bool   stats_valid = false;
    double robust_center = 0;   // median
    double robust_sigma  = 0;   // 1.4826 * MAD
    double mad           = 0;
    int    n_used        = 0;
    int    n_clipped     = 0;
};

struct LivetimeCut {
    bool        enabled = false;
    std::string source  = "ref";
    int         channel = 0;
    ChannelCut  cut;
};

// Live-charge integration over kept slow-event intervals.  Disabled if no
// `charge` block is present in the cut JSON; otherwise sums
//   Σ live_fraction × Δt × ½(I_i + I_{i+1})
// over each adjacent pair of accepted checkpoints, where live_fraction
// is the slice-local DSC2 livetime at the right endpoint and I is the
// configured beam-current EPICS channel.  Output units are
// (beam_current_unit · seconds).
struct ChargeCfg {
    bool        enabled              = false;
    std::string beam_current_channel;   // EPICS channel name to read
};

// Run-splitting by a slow-control PV level (the "split" block, see the file
// header): each checkpoint is classified by its forward-filled `channel` reading.
struct SplitCfg {
    bool        enabled  = false;
    std::string channel;                       // EPICS channel to watch

    // Level thresholds (require full > empty).  A checkpoint is full when its
    // PV >= full_thresh, empty when PV <= empty_thresh, dropped otherwise.
    double full_thresh  = 0.0;
    double empty_thresh = 0.0;

    // Optional extra margin: also drop the ± guard_checkpoints points on either
    // side of any state edge (0 = rely on the dead zone alone).
    int    guard_checkpoints = 0;

    // Output suffixes for the full / empty files (and their reports).
    std::string label_full  = "full";
    std::string label_empty = "empty";
};

struct CutConfig {
    LivetimeCut                       livetime;
    std::map<std::string, ChannelCut> epics;
    ChargeCfg                         charge;
    SplitCfg                          split;
    json                              raw;            // echoed in the report
};

void parse_abs(const json &j, AbsCut &a)
{
    if (!j.is_object()) return;
    if (j.contains("min") && j["min"].is_number()) {
        a.has_min = true; a.min_val = j["min"].get<double>();
    }
    if (j.contains("max") && j["max"].is_number()) {
        a.has_max = true; a.max_val = j["max"].get<double>();
    }
}

void parse_channel_cut(const json &j, ChannelCut &c)
{
    if (j.contains("abs")) parse_abs(j["abs"], c.abs);
    if (j.contains("rel_rms") && j["rel_rms"].is_number()) {
        c.has_rel_rms = true;
        c.rel_rms_n   = j["rel_rms"].get<double>();
    }
    if (j.contains("gated_by")) {
        const auto &g = j["gated_by"];
        if (g.is_string()) {
            c.gated_by.push_back(g.get<std::string>());
        } else if (g.is_array()) {
            for (const auto &item : g) {
                if (item.is_string()) c.gated_by.push_back(item.get<std::string>());
            }
        }
    }
}

bool load_cuts(const std::string &path, CutConfig &cfg)
{
    std::ifstream f(path);
    if (!f) {
        std::cerr << "replay_filter: cannot open cuts file: " << path << "\n";
        return false;
    }
    json j;
    try {
        j = json::parse(f, nullptr, true, /*allow_comments=*/true);
    } catch (json::parse_error &e) {
        std::cerr << "replay_filter: cuts JSON parse error: " << e.what() << "\n";
        return false;
    }
    cfg.raw = j;

    if (j.contains("livetime")) {
        const auto &lj = j["livetime"];
        cfg.livetime.enabled = true;
        if (lj.contains("source"))  cfg.livetime.source  = lj["source"].get<std::string>();
        if (lj.contains("channel")) cfg.livetime.channel = lj["channel"].get<int>();
        parse_channel_cut(lj, cfg.livetime.cut);
    }

    if (j.contains("epics") && j["epics"].is_object()) {
        for (auto it = j["epics"].begin(); it != j["epics"].end(); ++it) {
            ChannelCut c;
            parse_channel_cut(it.value(), c);
            cfg.epics[it.key()] = c;
        }
    }

    // Charge integration is opt-in: the cut JSON must name the EPICS
    // channel that carries the beam current.  Anything else (e.g. units)
    // is documented in the report so downstream consumers can convert.
    if (j.contains("charge") && j["charge"].is_object()) {
        const auto &cj = j["charge"];
        if (cj.contains("beam_current") && cj["beam_current"].is_string()) {
            cfg.charge.enabled = true;
            cfg.charge.beam_current_channel = cj["beam_current"].get<std::string>();
        }
    }

    // Run-splitting by PV level (opt-in: needs a `channel` plus `full` and
    // `empty` thresholds, with full > empty).
    if (j.contains("split") && j["split"].is_object()) {
        const auto &sj = j["split"];
        auto &s = cfg.split;
        if (sj.contains("channel") && sj["channel"].is_string())
            s.channel = sj["channel"].get<std::string>();
        const bool has_full  = sj.contains("full")  && sj["full"].is_number();
        const bool has_empty = sj.contains("empty") && sj["empty"].is_number();
        if (has_full)  s.full_thresh  = sj["full"].get<double>();
        if (has_empty) s.empty_thresh = sj["empty"].get<double>();
        if (sj.contains("guard_checkpoints") && sj["guard_checkpoints"].is_number())
            s.guard_checkpoints = std::max(0, sj["guard_checkpoints"].get<int>());
        if (sj.contains("labels") && sj["labels"].is_array()
            && sj["labels"].size() >= 2) {
            s.label_full  = sj["labels"][0].get<std::string>();
            s.label_empty = sj["labels"][1].get<std::string>();
        }
        // Need a channel and a well-ordered (full > empty) threshold pair.
        const bool ok = !s.channel.empty() && has_full && has_empty
                        && s.full_thresh > s.empty_thresh;
        s.enabled = ok;
        if (!s.channel.empty() && !ok) {
            if (!has_full || !has_empty)
                std::cerr << "replay_filter: split needs both `full` and `empty` "
                             "level thresholds — split disabled\n";
            else
                std::cerr << "replay_filter: split `full` (" << s.full_thresh
                          << ") must be > `empty` (" << s.empty_thresh
                          << ") — split disabled\n";
        }
    }
    return true;
}

// ── Robust statistics: median + MAD ──────────────────────────────────────
// Uses median absolute deviation (Hampel 1974).  More robust to heavy
// outliers than iterative sigma clipping: a single bad reading shifts the
// median negligibly and inflates MAD only via its own contribution.  The
// 1.4826 factor makes σ̂ a consistent estimator of stddev for a normal
// distribution, so cut thresholds in `rel_rms: N` keep their intuitive
// "N standard deviations" meaning.

struct RobustStats {
    double median    = 0;
    double mad       = 0;        // raw MAD
    double sigma     = 0;        // 1.4826 * MAD
    int    n_used    = 0;
    int    n_clipped = 0;        // points outside the cut band (informational)
};

double median_inplace(std::vector<double> &xs)
{
    if (xs.empty()) return 0;
    size_t n = xs.size();
    auto mid = xs.begin() + n / 2;
    std::nth_element(xs.begin(), mid, xs.end());
    double m = *mid;
    if ((n & 1u) == 0) {
        // even count: average mid with the largest of the lower half
        auto max_lo = std::max_element(xs.begin(), mid);
        m = 0.5 * (m + *max_lo);
    }
    return m;
}

RobustStats robust_mad(const std::vector<double> &xs, double n_sigma_for_clip_count)
{
    RobustStats r;
    if (xs.empty()) return r;
    r.n_used = static_cast<int>(xs.size());

    std::vector<double> tmp = xs;
    r.median = median_inplace(tmp);

    std::vector<double> dev;
    dev.reserve(xs.size());
    for (double x : xs) dev.push_back(std::fabs(x - r.median));
    r.mad   = median_inplace(dev);
    r.sigma = 1.4826 * r.mad;

    if (r.sigma > 0 && n_sigma_for_clip_count > 0) {
        for (double x : xs)
            if (std::fabs(x - r.median) > n_sigma_for_clip_count * r.sigma)
                ++r.n_clipped;
    }
    return r;
}

// ── Tree readers ─────────────────────────────────────────────────────────

// Pre-scan the events/recon tree across all input files and build a
// lookup event_num → ti_ticks (the 48-bit TI timestamp).  Only `event_num`
// and `timestamp` branches are activated, so this is fast even on
// millions-of-event runs.  Used for the time anchor and as the TI-tick
// fallback for EPICS rows from replays without ti_ticks_at_arrival.
bool build_evn_to_ticks(const std::vector<std::string> &files,
                        const std::string              &tree_name,
                        std::map<int32_t, int64_t>     &out)
{
    int       event_num = 0;
    long long timestamp = 0;
    for (const auto &path : files) {
        std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
        if (!f || f->IsZombie()) {
            std::cerr << "replay_filter: cannot open " << path << "\n";
            return false;
        }
        TTree *t = dynamic_cast<TTree *>(f->Get(tree_name.c_str()));
        if (!t) continue;

        t->SetBranchStatus("*", 0);
        if (t->GetBranch("event_num")) {
            t->SetBranchStatus("event_num", 1);
            t->SetBranchAddress("event_num", &event_num);
        } else {
            std::cerr << "replay_filter: '" << tree_name
                      << "' has no event_num branch in " << path << "\n";
            return false;
        }
        if (t->GetBranch("timestamp")) {
            t->SetBranchStatus("timestamp", 1);
            t->SetBranchAddress("timestamp", &timestamp);
        } else {
            std::cerr << "replay_filter: '" << tree_name
                      << "' has no timestamp branch in " << path << "\n";
            return false;
        }

        Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            // First-write-wins: in case of duplicate event_num across files
            // (shouldn't happen for a single run), keep the earliest tick.
            out.emplace(event_num, static_cast<int64_t>(timestamp));
        }
    }
    return true;
}

// ── Cut evaluation ───────────────────────────────────────────────────────

bool eval_channel_cut(const ChannelCut &c, double value)
{
    if (c.abs.has_min && !(value >= c.abs.min_val)) return false;
    if (c.abs.has_max && !(value <= c.abs.max_val)) return false;
    if (c.has_rel_rms && c.stats_valid && c.robust_sigma > 0) {
        if (std::fabs(value - c.robust_center) > c.rel_rms_n * c.robust_sigma)
            return false;
    }
    return true;
}

// ── Report point ─────────────────────────────────────────────────────────
//
// Each report point carries:
//   * `associated_evn` — the physics event the slow row is anchored to.
//     Scaler rows: their own event_number (the SYNC physics event whose
//     readout included the scaler bank).  EPICS rows: event_number_at_-
//     arrival (the most recent physics event seen at the time of the
//     EPICS event).  Both are integer keys into the events/recon tree.
//   * `associated_timestamp` (relative seconds) — the TI 48-bit tick carried
//     by the slow row (scalers: the SYNC event's own timestamp; EPICS:
//     ti_ticks_at_arrival, or the events-tree tick of associated_evn for
//     replays without that branch), in seconds since the earliest TI tick
//     seen in the run.  null when no tick is known (e.g.
//     event_number_at_arrival = -1 — EPICS arrived before any physics).
//   * `unix_time` (absolute Unix seconds) — from the 0xE112 HEAD bank.
//     Native for EPICS rows.  Explicitly null for scaler rows: the
//     scaler's cached unix_time can be a SYNC interval old, and emitting
//     it would invite mis-alignment.  Charts that need absolute time
//     should plot associated_timestamp and use any EPICS unix_time as
//     the absolute anchor (one EPICS pin is enough for the whole run).
struct ReportPoint {
    std::string channel;
    bool        pass;
    int32_t     event_number;     // = associated_evn
    bool        has_assoc_t;
    double      assoc_t_rel;      // seconds since the run's earliest event
    bool        has_unix_time;
    int64_t     unix_time;
    double      value;            // NaN ⇒ value not yet seen
    int         side = 0;         // split side this checkpoint belongs to;
                                  // -1 = dropped (ramp/guard). Backfilled in
                                  // phase 5, stays 0 when split is off.
};

// ── Main pipeline ────────────────────────────────────────────────────────

bool detect_event_tree(const std::string &path, std::string &name)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) return false;
    if (f->Get("events")) { name = "events"; return true; }
    if (f->Get("recon"))  { name = "recon";  return true; }
    return false;
}

std::string insert_before_root(const std::string &path, const std::string &suffix)
{
    auto p = std::filesystem::path(path);
    std::string name = p.filename().string();
    const std::string ext = ".root";
    if (name.size() >= ext.size() &&
        name.compare(name.size() - ext.size(), ext.size(), ext) == 0)
        name.insert(name.size() - ext.size(), suffix);
    else
        name += suffix;
    return name;
}

std::string filtered_output_path(const std::string &out_dir,
                                 const std::string &input_path,
                                 const std::string &label = "")
{
    auto p = std::filesystem::path(input_path);
    const std::string name = p.filename().string();
    const std::string marker = "_recon_";
    const std::string ext = ".root";
    const auto marker_pos = name.rfind(marker);
    if (marker_pos != std::string::npos &&
        name.size() > marker_pos + marker.size() + ext.size() &&
        name.compare(name.size() - ext.size(), ext.size(), ext) == 0) {
        const std::string index = name.substr(
            marker_pos + marker.size(),
            name.size() - marker_pos - marker.size() - ext.size());
        const std::string output_name =
            name.substr(0, marker_pos) + "_filter" +
            (label.empty() ? std::string() : "_" + label) +
            "_" + index + ext;
        return (std::filesystem::path(out_dir) / output_name).string();
    }

    const std::string unmerged_suffix = "_recon.root";
    if (name.size() > unmerged_suffix.size() &&
        name.compare(name.size() - unmerged_suffix.size(),
                     unmerged_suffix.size(), unmerged_suffix) == 0) {
        const std::string stem = name.substr(0, name.size() - unmerged_suffix.size());
        const auto index_pos = stem.rfind('.');
        if (index_pos != std::string::npos && index_pos + 1 < stem.size()) {
            std::string prefix = stem.substr(0, index_pos);
            if (prefix.size() >= 5 &&
                prefix.compare(prefix.size() - 5, 5, ".evio") == 0)
                prefix.resize(prefix.size() - 5);
            const std::string output_name =
                prefix + "_filter" +
                (label.empty() ? std::string() : "_" + label) +
                "_" + stem.substr(index_pos + 1) + ext;
            return (std::filesystem::path(out_dir) / output_name).string();
        }
    }

    std::string suffix = "_filter" + (label.empty() ? std::string() : "_" + label);
    return (std::filesystem::path(out_dir) / insert_before_root(input_path, suffix)).string();
}

std::string run_epics_path(const std::string &out_dir, int run_number)
{
    char name[64];
    std::snprintf(name, sizeof(name), "prad_%06d_epics.root", run_number);
    return (std::filesystem::path(out_dir) / name).string();
}

std::string run_report_path(const std::string &out_dir, int run_number)
{
    char name[80];
    std::snprintf(name, sizeof(name), "prad_%06d_filter_report.json", run_number);
    return (std::filesystem::path(out_dir) / name).string();
}

int run(const std::vector<std::string> &input_files,
        const std::string &output_path,
        const std::string &cuts_path,
        const std::string &report_path,
        int run_number_override,
        int num_threads)
{
    CutConfig cuts;
    if (!load_cuts(cuts_path, cuts)) return 1;

    // ---------- Phase 1: load slow streams into memory ----------
    // Rows stay in input-file order: the phase-7 re-reads iterate the files
    // in that order and look up each row's verdict by its load-order index.
    std::vector<analysis::ScalerRow> scalers;
    std::vector<analysis::EpicsRow>  epics_rows;
    if (!analysis::LoadScalerRows(input_files, scalers, "replay_filter"))   return 1;
    if (!analysis::LoadEpicsRows (input_files, epics_rows, "replay_filter")) return 1;
    std::cerr << "replay_filter: loaded " << scalers.size() << " scaler + "
              << epics_rows.size() << " epics rows from "
              << input_files.size() << " file(s)\n";

    int run_number = run_number_override;
    if (run_number < 0) run_number = analysis::get_run_int(input_files.front());
    if (run_number < 0 && !scalers.empty())    run_number = (int)scalers.front().run_number;
    if (run_number < 0 && !epics_rows.empty()) run_number = (int)epics_rows.front().run_number;

    struct FileInfo {
        std::string path;
        size_t scaler_offset = 0;
        size_t epics_offset = 0;
        Long64_t scaler_entries = 0;
        Long64_t epics_entries = 0;
    };
    std::vector<FileInfo> file_info;
    file_info.reserve(input_files.size());
    size_t scaler_offset = 0, epics_offset = 0;
    for (const auto &path : input_files) {
        FileInfo info;
        info.path = path;
        info.scaler_offset = scaler_offset;
        info.epics_offset = epics_offset;
        info.scaler_entries = analysis::TreeEntries(path, "scalers");
        info.epics_entries = analysis::TreeEntries(path, "epics");
        scaler_offset += static_cast<size_t>(info.scaler_entries);
        epics_offset += static_cast<size_t>(info.epics_entries);
        file_info.push_back(std::move(info));
    }

    // Sort scalers once and precompute delta livetime per row.  Cuts evaluate
    // and report against the slice-local live fraction (Δgated / Δungated),
    // not the run-cumulative ratio cached on each row, which would dilute a
    // recent dropout behind minutes of good livetime.
    auto sc_order = analysis::SortByEvent(scalers);
    std::vector<double> delta_live_pct;
    if (cuts.livetime.enabled)
        delta_live_pct = analysis::DeltaLivetime(scalers, sc_order, cuts.livetime.source,
                                                 cuts.livetime.channel, 100.0);

    // ---------- Phase 2: robust stats for rel_rms cuts ----------
    // Ungated channels first, then gated ones (see ChannelCut::gated_by);
    // gating chains are intentionally not supported to keep the JSON
    // unambiguous.
    auto fill_stats = [&](ChannelCut &c, const std::vector<double> &xs) {
        auto rs = robust_mad(xs, c.rel_rms_n);
        c.stats_valid   = (rs.n_used > 1) && rs.sigma > 0;
        c.robust_center = rs.median;
        c.robust_sigma  = rs.sigma;
        c.mad           = rs.mad;
        c.n_used        = rs.n_used;
        c.n_clipped     = rs.n_clipped;
    };

    // Walk EPICS rows in event-number order so forward-fill of the gating
    // channel reflects the actual time sequence.
    const auto ep_order = analysis::SortByEvent(epics_rows);
    auto values_of = [&](const std::string &channel) {
        std::vector<double> xs;
        for (size_t oi : ep_order) {
            auto it = epics_rows[oi].updates.find(channel);
            if (it != epics_rows[oi].updates.end()) xs.push_back(it->second);
        }
        return xs;
    };

    // 1. livetime (independent of EPICS).
    if (cuts.livetime.enabled && cuts.livetime.cut.has_rel_rms) {
        std::vector<double> xs;
        xs.reserve(scalers.size());
        for (size_t i = 0; i < scalers.size(); ++i) {
            double v = delta_live_pct[i];
            if (v >= 0) xs.push_back(v);
        }
        fill_stats(cuts.livetime.cut, xs);
    }

    // 2. ungated EPICS channels (stats from all observed values).
    for (auto &kv : cuts.epics) {
        if (!kv.second.has_rel_rms) continue;
        if (!kv.second.gated_by.empty()) continue;
        fill_stats(kv.second, values_of(kv.first));
    }

    // 3. gated EPICS channels (stats from rows where every gate's cut passed).
    for (auto &kv : cuts.epics) {
        if (!kv.second.has_rel_rms) continue;
        if (kv.second.gated_by.empty()) continue;

        // Resolve every gating channel.  If any is missing, fall back to
        // ungated stats (and log) — partial gating would be misleading.
        std::vector<const ChannelCut *> gates;
        gates.reserve(kv.second.gated_by.size());
        bool gates_ok = true;
        for (const auto &gname : kv.second.gated_by) {
            auto it = cuts.epics.find(gname);
            if (it == cuts.epics.end()) {
                std::cerr << "replay_filter: channel '" << kv.first
                          << "' is gated_by '" << gname
                          << "' which is not configured — falling back to ungated stats\n";
                gates_ok = false;
                break;
            }
            if (!it->second.gated_by.empty()) {
                std::cerr << "replay_filter: channel '" << kv.first
                          << "' gated_by '" << gname
                          << "' which is itself gated — chains not supported, "
                             "ignoring the inner gating\n";
            }
            gates.push_back(&it->second);
        }
        if (!gates_ok) {
            fill_stats(kv.second, values_of(kv.first));
            continue;
        }

        std::vector<double> xs;
        std::map<std::string, double> cur_eps;        // forward-fill across rows
        for (size_t oi : ep_order) {
            const auto &row = epics_rows[oi];
            for (const auto &up : row.updates) cur_eps[up.first] = up.second;

            auto val_it = row.updates.find(kv.first);
            if (val_it == row.updates.end()) continue;

            bool all_pass = true;
            for (size_t gi = 0; gi < gates.size(); ++gi) {
                auto gv_it = cur_eps.find(kv.second.gated_by[gi]);
                if (gv_it == cur_eps.end()
                    || !eval_channel_cut(*gates[gi], gv_it->second)) {
                    all_pass = false;
                    break;
                }
            }
            if (!all_pass) continue;
            xs.push_back(val_it->second);
        }
        fill_stats(kv.second, xs);
    }

    // ---------- Phase 3: anchor for relative associated_timestamp ----------
    // Slow rows carry their own TI tick (scalers: the SYNC event's
    // info.timestamp; EPICS: ti_ticks_at_arrival, captured at decode time so
    // it is independent of whether the anchor event was written to the events
    // tree).  The events tree is still detected (phase 7 copies it) and
    // pre-scanned: its ticks feed the anchor and are the fallback for EPICS
    // rows from replays without that branch.
    std::string ev_tree_name;
    if (!detect_event_tree(input_files.front(), ev_tree_name)) {
        std::cerr << "replay_filter: no events/recon tree in "
                  << input_files.front() << "\n";
        return 1;
    }
    std::map<int32_t, int64_t> evn_to_ticks;
    if (!build_evn_to_ticks(input_files, ev_tree_name, evn_to_ticks)) return 1;

    // Anchor = smallest TI tick across every source we have.  Considering
    // the slow rows (not just the events tree) keeps anchor monotonicity
    // when the events tree skips early events (e.g. trigger filter).
    int64_t  ti_anchor    = 0;
    bool     anchor_set   = false;
    auto consider_tick = [&](int64_t t) {
        if (t <= 0) return;
        if (!anchor_set || t < ti_anchor) { ti_anchor = t; anchor_set = true; }
    };
    for (const auto &kv : evn_to_ticks) consider_tick(kv.second);
    for (const auto &s  : scalers)      consider_tick(s.ti_ticks);
    for (const auto &e  : epics_rows)   consider_tick(e.ti_ticks);

    // ---------- Phase 4: walk merged timeline, mark good/bad ----------
    // Iterate via the sc_order / ep_order permutations so the parallel
    // verdict vectors stay aligned with the load-order vectors (used in
    // phase 7).
    std::vector<bool> scaler_verdict(scalers.size(),  false);
    std::vector<bool> epics_verdict (epics_rows.size(), false);

    struct Checkpoint {
        int32_t event_number;
        int64_t unix_time;
        int64_t ti_ticks;        // 0 if unknown — pair contributes no charge
        double  live_fraction;   // [0, 1] from cur_lt/100, NaN if unset
        double  beam_current;    // forward-filled, NaN if not seen yet
        bool    overall_pass;
        double  split_pv;        // forward-filled split channel, NaN if unset
        bool    is_scaler;       // true ⇒ orig indexes `scalers`, else `epics`
        size_t  orig;            // load-order index into the source vector
    };
    std::vector<Checkpoint>  timeline;
    std::vector<ReportPoint> report_points;

    // Forward-fill state for cut evaluation only.
    double                            cur_lt    = -1.0;   // % livetime
    std::map<std::string, double>     cur_eps;
    int64_t                           last_unix = 0;

    size_t i_sc = 0, i_ep = 0;
    while (i_sc < sc_order.size() || i_ep < ep_order.size()) {
        const bool take_sc =
            (i_sc < sc_order.size()) &&
            (i_ep >= ep_order.size() ||
             scalers[sc_order[i_sc]].event_number <=
             epics_rows[ep_order[i_ep]].event_number);

        int32_t cp_evn   = 0;
        int64_t cp_unix  = 0;
        int64_t cp_ticks = 0;        // TI tick captured on the slow row itself
        size_t  orig     = 0;
        bool    is_sc    = take_sc;

        bool emit_unix = false;     // true only for EPICS rows
        if (take_sc) {
            orig = sc_order[i_sc++];
            const auto &s = scalers[orig];
            cp_evn   = s.event_number;
            cp_ticks = s.ti_ticks;   // SYNC event's own info.timestamp
            // Slice-local live fraction (Δgated / Δungated), see
            // DeltaLivetime.  The first row's predecessor is (0, 0).
            cur_lt = cuts.livetime.enabled ? delta_live_pct[orig] : -1.0;
            // Scaler's cached unix_time is intentionally ignored (see
            // ReportPoint::unix_time).
            cp_unix = last_unix;
        } else {
            orig = ep_order[i_ep++];
            const auto &e = epics_rows[orig];
            cp_evn   = e.event_number;
            cp_ticks = e.ti_ticks;   // ti_ticks_at_arrival, captured at decode
            for (const auto &kv : e.updates) cur_eps[kv.first] = kv.second;
            if (e.unix_time > 0) last_unix = e.unix_time;
            cp_unix    = last_unix;
            emit_unix  = (e.unix_time > 0);
            // Replays without ti_ticks_at_arrival: join on
            // event_number_at_arrival.  Misses when the anchor event itself
            // was filtered out at replay time.
            if (cp_ticks <= 0 && cp_evn >= 0) {
                auto it = evn_to_ticks.find(cp_evn);
                if (it != evn_to_ticks.end()) cp_ticks = it->second;
            }
        }

        bool   pt_has_t = (cp_ticks > 0) && anchor_set;
        double pt_t     = pt_has_t
                          ? (cp_ticks - ti_anchor) * fdec::TI_TICK_SEC : 0.0;
        bool   pt_has_unix = emit_unix;
        int64_t pt_unix    = emit_unix ? (int64_t)last_unix : 0;

        // Per-channel report points with forward-filled values (dense traces
        // for plotting).
        if (cuts.livetime.enabled) {
            bool has  = (cur_lt >= 0);
            double v  = has ? cur_lt : std::numeric_limits<double>::quiet_NaN();
            bool pass = has && eval_channel_cut(cuts.livetime.cut, cur_lt);
            report_points.push_back({"livetime", pass, cp_evn,
                                     pt_has_t, pt_t,
                                     pt_has_unix, pt_unix, v});
        }
        for (const auto &kv : cuts.epics) {
            auto   it   = cur_eps.find(kv.first);
            bool   has  = (it != cur_eps.end());
            double v    = has ? it->second
                              : std::numeric_limits<double>::quiet_NaN();
            bool   pass = has && eval_channel_cut(kv.second, v);
            report_points.push_back({"epics:" + kv.first, pass, cp_evn,
                                     pt_has_t, pt_t,
                                     pt_has_unix, pt_unix, v});
        }

        // Overall verdict at this checkpoint = AND of every configured cut.
        // Channels that haven't reported yet count as "fail" — the user's
        // spec says events bracketed by an undefined endpoint are dropped.
        bool overall = true;
        if (cuts.livetime.enabled) {
            if (cur_lt < 0 || !eval_channel_cut(cuts.livetime.cut, cur_lt))
                overall = false;
        }
        for (const auto &kv : cuts.epics) {
            auto it = cur_eps.find(kv.first);
            if (it == cur_eps.end() || !eval_channel_cut(kv.second, it->second)) {
                overall = false;
            }
        }
        // Live fraction at this checkpoint (forward-filled %, scaled to
        // [0, 1]); beam current pulled from the configured EPICS channel
        // also via forward-fill.  Missing values stay NaN so the charge
        // integration knows to skip the surrounding pair.
        const double cp_live_fraction = (cur_lt >= 0)
            ? cur_lt * 0.01 : std::numeric_limits<double>::quiet_NaN();
        double cp_current = std::numeric_limits<double>::quiet_NaN();
        if (cuts.charge.enabled) {
            auto it = cur_eps.find(cuts.charge.beam_current_channel);
            if (it != cur_eps.end()) cp_current = it->second;
        }
        // Split PV (forward-filled like the others) so the transition scan
        // below sees a dense trace even on checkpoints carrying no update.
        double cp_split_pv = std::numeric_limits<double>::quiet_NaN();
        if (cuts.split.enabled) {
            auto it = cur_eps.find(cuts.split.channel);
            if (it != cur_eps.end()) cp_split_pv = it->second;
        }

        timeline.push_back({cp_evn, cp_unix, cp_ticks,
                            cp_live_fraction, cp_current, overall,
                            cp_split_pv, is_sc, orig});
        if (is_sc) scaler_verdict[orig] = overall;
        else       epics_verdict [orig] = overall;
    }

    // ---------- Phase 5: classify each checkpoint by PV level ----------
    // side 0 = full (label_full), 1 = empty (label_empty), -1 = dropped
    // (ramp / guard margin).  With split off, side[] is 0 everywhere so the
    // keep/charge/output code below is one path that runs for the present
    // side(s).  `split_active` mirrors cuts.split.enabled but degrades to false
    // if the channel never reports / never reaches either level.
    bool split_active = cuts.split.enabled;
    std::vector<int> side(timeline.size(), 0);
    std::vector<int> scaler_side(scalers.size(), 0);
    std::vector<int> epics_side (epics_rows.size(), 0);

    int  n_state_transitions = 0;
    json transitions = json::array();   // {evn, from, to, timestamp} per edge

    if (split_active) {
        const auto &S = cuts.split;
        // After phase 4's forward-fill the only NaNs are the head before the
        // channel's first report; back-fill that head with the first reading
        // (the run's starting state) so a pure run keeps its leading events.
        double first_val = std::numeric_limits<double>::quiet_NaN();
        for (const auto &cp : timeline)
            if (std::isfinite(cp.split_pv)) { first_val = cp.split_pv; break; }

        if (!std::isfinite(first_val)) {
            std::cerr << "replay_filter: split channel '" << S.channel
                      << "' never reported — cannot classify, writing the single "
                         "unsplit output instead\n";
            split_active = false;
        } else {
            // Raw per-checkpoint level classification.
            for (size_t k = 0; k < timeline.size(); ++k) {
                double v = timeline[k].split_pv;
                if (!std::isfinite(v)) v = first_val;          // back-fill head
                side[k] = (v >= S.full_thresh)  ? 0
                        : (v <= S.empty_thresh) ? 1
                        :                         -1;          // ramp / dead zone
            }
            // Optional extra margin: drop ± guard_checkpoints around any edge
            // between two adjacent, different, non-dropped states (a transition
            // faster than the checkpoint spacing leaves no dead-zone point
            // between the two levels).  Measured against the raw labels.
            if (S.guard_checkpoints > 0) {
                const std::vector<int> raw = side;
                const int N = S.guard_checkpoints;
                for (size_t k = 1; k < raw.size(); ++k)
                    if (raw[k] >= 0 && raw[k - 1] >= 0 && raw[k] != raw[k - 1])
                        for (int j = std::max(0, (int)k - N);
                             j < std::min((int)raw.size(), (int)k + N); ++j)
                            side[j] = -1;
            }
            // Count full<->empty transitions (using the last non-dropped state,
            // so a dead-zone gap is one transition not two) and record each
            // edge at the first checkpoint of the new state.
            int prev_state = -1;
            for (size_t k = 0; k < timeline.size(); ++k) {
                const int sd = side[k];
                if (sd < 0) continue;
                if (prev_state >= 0 && sd != prev_state) {
                    ++n_state_transitions;
                    const auto &cp = timeline[k];
                    const bool has_t = (cp.ti_ticks > 0 && anchor_set);
                    transitions.push_back({
                        {"evn",       cp.event_number},
                        {"from",      prev_state == 0 ? "full" : "empty"},
                        {"to",        sd == 0 ? "full" : "empty"},
                        {"timestamp", has_t
                            ? json((cp.ti_ticks - ti_anchor) * fdec::TI_TICK_SEC)
                            : json(nullptr)},
                    });
                }
                prev_state = sd;
            }
            // No checkpoint reached either level (PV sat in the dead zone the
            // whole run, or thresholds don't bracket the data) → nothing to
            // split.  Degrade to a single unsplit output rather than emit
            // empty files.
            int nf = 0, ne = 0;
            for (int sd : side) { if (sd == 0) ++nf; else if (sd == 1) ++ne; }
            if (nf == 0 && ne == 0) {
                std::cerr << "replay_filter: split channel '" << S.channel
                          << "' never reached full (>=" << S.full_thresh
                          << ") or empty (<=" << S.empty_thresh << ") — no split, "
                             "writing the single unsplit output instead\n";
                std::fill(side.begin(), side.end(), 0);
                split_active = false;
            }
        }

        // Back-map to load-order rows + tag report points with their side.
        if (split_active) {
            for (size_t k = 0; k < timeline.size(); ++k) {
                if (timeline[k].is_scaler) scaler_side[timeline[k].orig] = side[k];
                else                       epics_side [timeline[k].orig] = side[k];
            }
            const size_t ppc = (cuts.livetime.enabled ? 1u : 0u) + cuts.epics.size();
            if (ppc > 0 && report_points.size() == timeline.size() * ppc)
                for (size_t k = 0; k < timeline.size(); ++k)
                    for (size_t j = 0; j < ppc; ++j)
                        report_points[k * ppc + j].side = side[k];
        }
    }

    // ---------- Phase 6: build per-side keep-intervals (lo, hi] + charge ----------
    // A pair (cp_{i-1}, cp_i) belongs to a side only when both endpoints carry
    // the same non-guard label; pairs that straddle the transition or touch
    // the guard contribute to neither output (those events are dropped).  The
    // charge integration is bucketed per side:
    //   * gated — both endpoints overall_pass.  Canonical post-cut number;
    //     matches the events written to that side's file.
    //   * ungated — every valid-data pair on side s regardless of the cut
    //     verdict, so users see how much charge the cuts dropped.
    std::vector<std::pair<int32_t, int32_t>> keep[2];
    std::vector<std::pair<int32_t, int32_t>> span[2];   // ungated per-side ranges
    analysis::ChargeSums charge[2];
    for (size_t i = 1; i < timeline.size(); ++i) {
        const auto &a = timeline[i - 1];
        const auto &b = timeline[i];
        const int sa = side[i - 1], sb = side[i];
        const int ps = (sa >= 0 && sa == sb) ? sa : -1;   // pair's side, or none
        if (ps < 0) continue;                              // straddles / guard
        span[ps].emplace_back(a.event_number, b.event_number);  // ungated range
        const bool good_pair = (a.overall_pass && b.overall_pass);
        if (good_pair)
            keep[ps].emplace_back(a.event_number, b.event_number);
        if (!cuts.charge.enabled) continue;
        charge[ps].AddPair(a.ti_ticks, b.ti_ticks, b.live_fraction,
                           a.beam_current, b.beam_current, good_pair);
    }
    auto in_intervals = [](const std::vector<std::pair<int32_t, int32_t>> &iv,
                           int32_t evn) -> bool {
        if (iv.empty()) return false;
        auto it = std::upper_bound(
            iv.begin(), iv.end(), evn,
            [](int32_t e, const std::pair<int32_t, int32_t> &p) { return e < p.first; });
        if (it == iv.begin()) return false;
        --it;
        return evn > it->first && evn <= it->second;
    };
    auto is_kept = [&](int s, int32_t evn) { return in_intervals(keep[s], evn); };
    auto in_span = [&](int s, int32_t evn) { return in_intervals(span[s], evn); };

    // ---------- Phase 7: write the output(s) — ROOT file(s) + JSON report ----------
    // One input: one ROOT file + report per present side (split off: side 0 →
    // output_path).  Several inputs: one ROOT file per (input, present side),
    // the run's unfiltered slow trees and one run-level report.  A row's
    // `good` flag in a side's slow trees is its overall verdict AND-ed with
    // "belongs to this side", so running live_charge on a side file
    // reproduces that side's post-cut charge directly.
    const bool is_recon = (ev_tree_name == "recon");

    // Checkpoints per state (n_cp_side[0]=full, [1]=empty) and dropped
    // (ramp / dead zone / margin) — reported in the split block, and the
    // present-state check (>0) that drives which files get written.
    int n_cp_side[2] = {0, 0}, n_cp_guard = 0;
    for (int sd : side) { if (sd < 0) ++n_cp_guard; else ++n_cp_side[sd]; }

    auto side_label = [&](int s) -> const std::string & {
        return s == 0 ? cuts.split.label_full : cuts.split.label_empty;
    };

    // Concatenate the scalers / epics / runinfo trees of the given inputs into
    // out.  scalers and epics get a `good` branch (the checkpoint's phase-4
    // verdict, restricted to side s when restrict_side) so downstream tools
    // can colour the traces by pass/fail without recomputing.
    auto write_slow_trees = [&](TFile &out, const std::vector<size_t> &file_indices,
                                int s, bool restrict_side) {
        {
            prad2::RawScalerData sc;
            bool good = false;
            out.cd();
            TTree *out_sc = new TTree("scalers", "PRad2 DSC2 scaler readouts (concatenated)");
            prad2::SetScalerWriteBranches(out_sc, sc);
            out_sc->Branch("good", &good, "good/O");
            for (size_t fi : file_indices) {
                std::unique_ptr<TFile> f(TFile::Open(file_info[fi].path.c_str(), "READ"));
                TTree *t = f ? dynamic_cast<TTree *>(f->Get("scalers")) : nullptr;
                if (!t) continue;
                prad2::SetScalerReadBranches(t, sc);
                size_t seq = file_info[fi].scaler_offset;
                Long64_t n = t->GetEntries();
                for (Long64_t i = 0; i < n; ++i, ++seq) {
                    t->GetEntry(i);
                    good = (seq < scaler_verdict.size()) ? scaler_verdict[seq] : false;
                    if (good && restrict_side)
                        good = (seq < scaler_side.size() && scaler_side[seq] == s);
                    out.cd();
                    out_sc->Fill();
                }
            }
            out.cd();
            out_sc->Write();
        }
        // Rows of inputs without ti_ticks_at_arrival get it from the
        // events-tree lookup, so the output is always self-contained.
        {
            prad2::RawEpicsData ep;
            prad2::EpicsVectorBindings ep_vecs;
            bool good = false;
            out.cd();
            TTree *out_ep = new TTree("epics", "PRad2 EPICS slow control (concatenated)");
            prad2::SetEpicsWriteBranches(out_ep, ep);
            out_ep->Branch("good", &good, "good/O");
            for (size_t fi : file_indices) {
                std::unique_ptr<TFile> f(TFile::Open(file_info[fi].path.c_str(), "READ"));
                TTree *t = f ? dynamic_cast<TTree *>(f->Get("epics")) : nullptr;
                if (!t) continue;
                prad2::SetEpicsReadBranches(t, ep);
                prad2::BindEpicsVectorBranches(t, ep, ep_vecs);
                size_t seq = file_info[fi].epics_offset;
                Long64_t n = t->GetEntries();
                for (Long64_t i = 0; i < n; ++i, ++seq) {
                    ep.ti_ticks_at_arrival = 0;
                    t->GetEntry(i);
                    if (ep.ti_ticks_at_arrival <= 0 && ep.event_number_at_arrival >= 0) {
                        auto eit = evn_to_ticks.find(ep.event_number_at_arrival);
                        if (eit != evn_to_ticks.end()) ep.ti_ticks_at_arrival = eit->second;
                    }
                    good = (seq < epics_verdict.size()) ? epics_verdict[seq] : false;
                    if (good && restrict_side)
                        good = (seq < epics_side.size() && epics_side[seq] == s);
                    out.cd();
                    out_ep->Fill();
                }
            }
            out.cd();
            out_ep->Write();
        }
        // runinfo (one row per CODA control event, including the DAQ-config
        // text on PRESTART) is run-scoped metadata and is not filtered.
        {
            prad2::RawRunInfo ri;
            std::string *sp = &ri.daq_config;
            out.cd();
            TTree *out_ri = new TTree("runinfo", "PRad2 control events / DAQ config (concatenated)");
            prad2::SetRunInfoWriteBranches(out_ri, ri);
            for (size_t fi : file_indices) {
                std::unique_ptr<TFile> f(TFile::Open(file_info[fi].path.c_str(), "READ"));
                TTree *t = f ? dynamic_cast<TTree *>(f->Get("runinfo")) : nullptr;
                if (!t) continue;
                prad2::SetRunInfoReadBranches(t, ri);
                t->SetBranchAddress("daq_config", &sp);
                Long64_t n = t->GetEntries();
                for (Long64_t i = 0; i < n; ++i) {
                    ri.daq_config.clear();
                    t->GetEntry(i);
                    out.cd();
                    out_ri->Fill();
                }
            }
            out.cd();
            out_ri->Write();
        }
    };

    // Write the events of input fi kept on side s, plus that input's slow
    // trees, to the new file outp.  The Set*WriteBranches helpers keep the
    // schema of a replay.  ok is false when outp cannot be created.
    struct WriteStats { int64_t n_in = 0, n_out = 0; bool ok = false; };
    auto write_file_side = [&](size_t fi, int s, const std::string &outp) -> WriteStats {
        WriteStats stats;
        std::unique_ptr<TFile> out(TFile::Open(outp.c_str(), "RECREATE"));
        if (!out || out->IsZombie()) {
            std::cerr << "replay_filter: cannot create " << outp << "\n";
            return stats;
        }
        std::unique_ptr<TFile> f(TFile::Open(file_info[fi].path.c_str(), "READ"));
        TTree *t = f ? dynamic_cast<TTree *>(f->Get(ev_tree_name.c_str())) : nullptr;

        // Copy the entries of t kept on side s to out_ev, whose branches point
        // into ev; reset(ev) runs before every read.
        auto copy_kept = [&](TTree *out_ev, auto &ev, auto reset) {
            if (t) {
                Long64_t n = t->GetEntries();
                if (!split_active) stats.n_in += n;   // split off: every event counts
                for (Long64_t i = 0; i < n; ++i) {
                    reset(ev);
                    t->GetEntry(i);
                    if (split_active && in_span(s, ev.event_num)) ++stats.n_in;
                    if (is_kept(s, ev.event_num)) {
                        out->cd();
                        out_ev->Fill();
                        ++stats.n_out;
                    }
                }
            }
            out->cd();
            out_ev->Write();
            out_ev->ResetBranchAddresses();
        };

        if (!is_recon) {
            auto ev = std::make_unique<prad2::RawEventData>();
            prad2::RawReadStatus status;
            prad2::RawVectorBindings vb;
            if (t) {
                status = prad2::SetRawReadBranches(t, *ev);
                prad2::BindRawVectorBranches(t, *ev, vb);
            }
            out->cd();
            TTree *out_ev = new TTree("events", "PRad2 filtered replay (raw)");
            prad2::SetRawWriteBranches(out_ev, *ev, status.has_peaks);
            copy_kept(out_ev, *ev, [](prad2::RawEventData &e) { e.clear_banks(); });
        } else {
            auto ev = std::make_unique<prad2::ReconEventData>();
            prad2::ReconMatchVectorBindings match_bind;
            prad2::RawVectorBindings vb;
            if (t) {
                prad2::SetReconReadBranches(t, *ev);
                prad2::BindReconMatchVectorBranches(t, *ev, match_bind);
                prad2::BindRawVectorBranches(t, *ev, vb);
            }
            out->cd();
            TTree *out_ev = new TTree("recon", "PRad2 filtered replay (recon)");
            prad2::SetReconWriteBranches(out_ev, *ev, false); // not x17_mode
            copy_kept(out_ev, *ev, [](prad2::ReconEventData &e) { e.clear(); });
        }
        write_slow_trees(*out, {fi}, s, split_active);
        out->Close();
        stats.ok = true;
        return stats;
    };

    // JSON report of side s, or of the whole run for s < 0 (every checkpoint,
    // per-channel counts and keep intervals of all sides, charge summed).
    auto build_report = [&](int s, int64_t n_in, int64_t n_out) -> json {
        const bool this_side_only = split_active && s >= 0;
        auto to_json_or_null = [](bool valid, double v) -> json {
            return valid ? json(v) : json(nullptr);
        };
        auto stats_for = [&](const ChannelCut &c) -> json {
            json j = {
                {"abs_min", c.abs.has_min ? json(c.abs.min_val) : json(nullptr)},
                {"abs_max", c.abs.has_max ? json(c.abs.max_val) : json(nullptr)},
            };
            if (c.has_rel_rms) {
                j["rel_rms"]       = c.rel_rms_n;
                j["robust_center"] = to_json_or_null(c.stats_valid, c.robust_center);
                j["robust_sigma"]  = to_json_or_null(c.stats_valid, c.robust_sigma);
                j["mad"]           = to_json_or_null(c.stats_valid, c.mad);
                // n_used is the count *after* gating (if any) — useful for
                // sanity-checking that the gating restriction left enough data
                // to compute meaningful stats.
                j["n_used"]        = c.n_used;
                j["n_clipped"]     = c.n_clipped;
                // Always an array — even single-gate cases — so downstream
                // tools can iterate without checking type.
                if (!c.gated_by.empty()) j["gated_by"] = c.gated_by;
            }
            return j;
        };

        json report;
        report["run_number"]    = run_number;
        report["input_files"]   = input_files;
        report["cuts"]          = cuts.raw;
        report["robust_method"] = "mad";   // 1.4826 * MAD as σ̂

        json stats = json::object();
        if (cuts.livetime.enabled) {
            json ls = stats_for(cuts.livetime.cut);
            ls["source"]  = cuts.livetime.source;
            ls["channel"] = cuts.livetime.channel;
            stats["livetime"] = std::move(ls);
        }
        for (const auto &kv : cuts.epics)
            stats["epics:" + kv.first] = stats_for(kv.second);
        report["stats"] = std::move(stats);

        int n_pass_cp = 0, n_fail_cp = 0;
        for (size_t k = 0; k < timeline.size(); ++k)
            if (!this_side_only || side[k] == s)
                (timeline[k].overall_pass ? n_pass_cp : n_fail_cp)++;
        const int n_slow = n_pass_cp + n_fail_cp;

        json keep_intervals = json::array();
        for (int ks = 0; ks < 2; ++ks)
            if (s < 0 || ks == s)
                for (const auto &p : keep[ks]) keep_intervals.push_back({p.first, p.second});

        // Per-channel breakdown — number of slow-event checkpoints where this
        // channel's cut accepted vs rejected the value.  Helps the user see
        // immediately which cut is doing the rejecting (e.g. "beam current
        // killed 80% of points, livetime barely matters").
        std::map<std::string, std::pair<int, int>> per_channel;   // ch → {pass, fail}
        for (const auto &p : report_points) {
            if (this_side_only && p.side != s) continue;
            auto &c = per_channel[p.channel];
            if (p.pass) ++c.first; else ++c.second;
        }
        json per_channel_json = json::object();
        for (const auto &kv : per_channel) {
            int pass = kv.second.first, fail = kv.second.second;
            int tot  = pass + fail;
            per_channel_json[kv.first] = {
                {"n_pass",    pass},
                {"n_fail",    fail},
                {"pass_rate", tot > 0 ? double(pass) / double(tot) : 0.0},
            };
        }

        report["summary"] = {
            {"n_slow_events",      n_slow},
            {"n_slow_pass",        n_pass_cp},
            {"n_slow_reject",      n_fail_cp},
            {"slow_pass_rate",     n_slow > 0 ? double(n_pass_cp) / double(n_slow) : 0.0},
            {"n_physics_in",       n_in},
            {"n_physics_pass",     n_out},
            {"n_physics_reject",   n_in - n_out},
            {"physics_pass_rate",  n_in > 0 ? double(n_out) / double(n_in) : 0.0},
            // Keep-interval count (each is a (lo, hi] range of accepted events).
            {"n_keep_intervals",   (int)keep_intervals.size()},
            {"per_channel",        per_channel_json},
        };
        report["keep_intervals"] = std::move(keep_intervals);

        // Split metadata: the level thresholds, how many checkpoints fell in
        // each state vs the dropped ramp, and every full<->empty transition
        // seen, so neither side's report needs the other's.
        if (split_active) {
            json split = {
                {"enabled",               true},
                {"channel",               cuts.split.channel},
                {"full_threshold",        cuts.split.full_thresh},
                {"empty_threshold",       cuts.split.empty_thresh},
                {"guard_checkpoints",     cuts.split.guard_checkpoints},
                {"n_checkpoints_full",    n_cp_side[0]},
                {"n_checkpoints_empty",   n_cp_side[1]},
                {"n_checkpoints_dropped", n_cp_guard},
                {"n_state_transitions",   n_state_transitions},
                {"transitions",           transitions},
                {"pure_run",              n_state_transitions == 0},
            };
            if (s >= 0) {
                split["side"]               = s;
                split["state"]              = s == 0 ? "full" : "empty";
                split["label"]              = side_label(s);
                split["n_checkpoints_this"] = n_cp_side[s];
            } else {
                split["labels"] = json::array({cuts.split.label_full, cuts.split.label_empty});
            }
            report["split"] = std::move(split);
        }

        // Live-charge integration over kept intervals.  Units: assume the
        // configured EPICS beam-current channel publishes in nA (true for
        // hallb_IPM2C21A_CUR and the other Hall B IPM scalers), so
        // value = Σ live_fraction · Δt · I  ⇒  nA · s = nC.  Also emit the
        // accumulated live time so the average current is recoverable.
        // value_nC / live_seconds / real_seconds are the gated sums and
        // ungated_* the same over every valid-data pair (see phase 6).
        if (cuts.charge.enabled) {
            analysis::ChargeSums c = charge[s < 0 ? 0 : s];
            if (s < 0) c += charge[1];
            report["live_charge"] = {
                {"value_nC",                   c.value_nC},
                {"unit",                       "nC"},
                {"beam_current_channel",       cuts.charge.beam_current_channel},
                {"beam_current_unit",          "nA"},
                {"live_seconds",               c.live_seconds},
                {"real_seconds",               c.real_seconds},
                {"ungated_value_nC",           c.ungated_value_nC},
                {"ungated_live_seconds",       c.ungated_live_seconds},
                {"ungated_real_seconds",       c.ungated_real_seconds},
                {"n_pairs_integrated",         c.n_pairs_integrated},
                {"n_pairs_skipped",            c.n_pairs_skipped},
                {"n_ungated_pairs_integrated", c.n_ungated_pairs_integrated},
                {"n_ungated_pairs_skipped",    c.n_ungated_pairs_skipped},
            };
        }

        json pts = json::array();
        pts.get_ptr<json::array_t *>()->reserve(report_points.size());
        for (const auto &p : report_points) {
            pts.push_back({
                {"channel",              p.channel},
                {"status",               p.pass ? "pass" : "fail"},
                {"associated_evn",       p.event_number},
                {"associated_timestamp", p.has_assoc_t ? json(p.assoc_t_rel) : json(nullptr)},
                {"unix_time",            p.has_unix_time ? json(p.unix_time) : json(nullptr)},
                {"value",                std::isnan(p.value) ? json(nullptr) : json(p.value)},
            });
        }
        report["points"] = std::move(pts);
        return report;
    };

    auto write_report = [](const json &report, const std::string &path) {
        std::ofstream of(path);
        if (!of) {
            std::cerr << "replay_filter: cannot write " << path << "\n";
            return false;
        }
        of << report.dump(2) << "\n";
        return true;
    };

    auto fmt_pct = [](double r) {
        std::ostringstream o;
        o << std::fixed << std::setprecision(2) << (r * 100.0) << "%";
        return o.str();
    };

    auto print_summary = [&](const json &report) {
        const json &sm = report.at("summary");
        std::cerr << "  slow events  : " << sm.at("n_slow_events").get<int>()
                  << "  pass=" << sm.at("n_slow_pass").get<int>()
                  << "  reject=" << sm.at("n_slow_reject").get<int>()
                  << "  rate=" << fmt_pct(sm.at("slow_pass_rate").get<double>())
                  << "\n";
        std::cerr << "  keep intervals: " << sm.at("n_keep_intervals").get<int>() << "\n";
        std::cerr << "  physics      : in=" << sm.at("n_physics_in").get<int64_t>()
                  << "  pass="  << sm.at("n_physics_pass").get<int64_t>()
                  << "  reject=" << sm.at("n_physics_reject").get<int64_t>()
                  << "  rate=" << fmt_pct(sm.at("physics_pass_rate").get<double>())
                  << "\n";
        if (report.contains("live_charge")) {
            const json &lc = report.at("live_charge");
            std::cerr << "  live charge  : " << std::fixed << std::setprecision(3)
                      << lc.at("value_nC").get<double>() << " nC over "
                      << lc.at("live_seconds").get<double>()
                      << " s live" << std::defaultfloat << "\n";
        }
        const json &per_channel = sm.at("per_channel");
        if (!per_channel.empty()) {
            std::cerr << "  per-channel reject:\n";
            for (auto it = per_channel.begin(); it != per_channel.end(); ++it) {
                const int pass = it->at("n_pass").get<int>();
                const int fail = it->at("n_fail").get<int>();
                std::cerr << "    " << it.key() << ": " << fail << " / " << (pass + fail)
                          << " (" << fmt_pct(double(fail) / std::max(1, pass + fail))
                          << ")\n";
            }
        }
    };

    if (input_files.size() > 1) {
        std::filesystem::create_directories(output_path);

        std::vector<std::pair<size_t, int>> jobs;   // (input file, side)
        for (size_t fi = 0; fi < input_files.size(); ++fi) {
            if (!split_active) jobs.emplace_back(fi, 0);
            else {
                if (n_cp_side[0] > 0) jobs.emplace_back(fi, 0);
                if (n_cp_side[1] > 0) jobs.emplace_back(fi, 1);
            }
        }

        std::vector<std::string> output_files(jobs.size());
        std::vector<WriteStats>  job_stats(jobs.size());
        std::mutex log_mtx;
        analysis::ParallelFor(jobs.size(), num_threads, [&](size_t ji, int) {
            const auto [fi, s] = jobs[ji];
            output_files[ji] = filtered_output_path(output_path, file_info[fi].path,
                                                    split_active ? side_label(s) : "");
            job_stats[ji] = write_file_side(fi, s, output_files[ji]);
            if (!job_stats[ji].ok) return;
            std::lock_guard<std::mutex> lk(log_mtx);
            std::cerr << "replay_filter: output ROOT     " << output_files[ji] << "\n";
        });

        int64_t n_in_total = 0, n_pass_phys = 0;
        for (const auto &st : job_stats) {
            if (!st.ok) return 1;
            n_in_total  += st.n_in;
            n_pass_phys += st.n_out;
        }

        const std::string slow_out = run_epics_path(output_path, run_number);
        {
            std::unique_ptr<TFile> out(TFile::Open(slow_out.c_str(), "RECREATE"));
            if (!out || out->IsZombie()) {
                std::cerr << "replay_filter: cannot create " << slow_out << "\n";
                return 1;
            }
            std::vector<size_t> all_files(input_files.size());
            for (size_t i = 0; i < all_files.size(); ++i) all_files[i] = i;
            write_slow_trees(*out, all_files, 0, false);
            out->Close();
        }
        std::cerr << "replay_filter: run slow ROOT  " << slow_out << "\n";

        json report = build_report(-1, n_in_total, n_pass_phys);
        report["output_files"]     = output_files;
        report["slow_output_file"] = slow_out;
        if (!write_report(report, report_path)) return 1;
        std::cerr << "replay_filter: report written to " << report_path << "\n";
        print_summary(report);
        return 0;
    }

    // Side s of the single input: ROOT file outp + report repp.
    auto write_side = [&](int s, const std::string &outp, const std::string &repp) -> int {
        const WriteStats st = write_file_side(0, s, outp);
        if (!st.ok) return 1;
        json report = build_report(s, st.n_in, st.n_out);
        report["output_file"] = outp;
        if (!write_report(report, repp)) return 1;
        const std::string tag = split_active ? "[" + side_label(s) + "] " : "";
        std::cerr << "replay_filter: " << tag << "report written to " << repp << "\n";
        std::cerr << "replay_filter: " << tag << "output ROOT     " << outp << "\n";
        print_summary(report);
        return 0;
    };

    // ---------- Dispatch: single output, or one file per present target state ----------
    // with_label suffixes the file/report stems with a side label.  Handles the
    // compound ".report.json" suffix so a report becomes
    // "<stem>_<label>.report.json" rather than "<stem>.report_<label>.json".
    auto with_label = [](const std::string &path, const std::string &label) {
        static const std::vector<std::string> compound = {".report.json"};
        for (const auto &suf : compound) {
            if (path.size() > suf.size() &&
                path.compare(path.size() - suf.size(), suf.size(), suf) == 0)
                return path.substr(0, path.size() - suf.size())
                       + "_" + label + suf;
        }
        auto dot = path.rfind('.');
        return (dot == std::string::npos)
            ? path + "_" + label
            : path.substr(0, dot) + "_" + label + path.substr(dot);
    };

    if (!split_active)        // split off or degraded → the single default output
        return write_side(0, output_path, report_path);

    // Level split: one labelled file per target state that actually occurred.
    const char *names[2] = {"full", "empty"};
    std::cerr << "replay_filter: split on '" << cuts.split.channel
              << "' by level (full >= " << cuts.split.full_thresh
              << ", empty <= " << cuts.split.empty_thresh << "): "
              << n_cp_side[0] << " full + " << n_cp_side[1] << " empty checkpoints, "
              << n_cp_guard << " dropped, " << n_state_transitions
              << " transition(s)\n";
    int rc = 0;
    for (int s = 0; s <= 1 && rc == 0; ++s) {
        if (n_cp_side[s] == 0) {
            std::cerr << "replay_filter: no " << names[s]
                      << "-target checkpoints — skipping " << names[s] << " file\n";
            continue;
        }
        rc = write_side(s, with_label(output_path, side_label(s)),
                           with_label(report_path, side_label(s)));
    }
    return rc;
}

void usage(const char *prog)
{
    std::cerr <<
        "Usage: " << prog << " <input.root> [more.root ...]\n"
        "       -o <output.root|output_dir>  -c <cuts.json> [-j <report.json>]\n"
        "       [-r <run_num>] [-t threads] [-h]\n"
        "\n"
        "Filters replayed ROOT files by slow-control cuts (DSC2 livetime\n"
        "+ EPICS).  Writes a single ROOT file with the kept physics events\n"
        "and the full scaler/epics streams concatenated, plus a JSON report\n"
        "with per-(cut, slow-event) pass/fail status for chart plotting.\n"
        "With multiple inputs, -o is an output directory; one filtered ROOT\n"
        "is written per input as prad_<run>_filter_<index>.root,\n"
        "plus prad_<run>_epics.root and one run-level JSON report.\n"
        "\n"
        "With a \"split\" block events are instead classified by a PV level\n"
        "(e.g. target cell pressure) and the full/empty subsets written to\n"
        "separate files <stem>_<full>.root / <stem>_<empty>.root.  Single-input\n"
        "mode writes side reports; multi-input mode writes one run-level report.\n"
        "PV>=full and PV<=empty\n"
        "select the two states; the in-between ramp lands in neither.  Only\n"
        "the states that occur are written, so a pure run yields one file.\n"
        "\n"
        "Cut JSON example:\n"
        "  {\n"
        "    \"livetime\": { \"source\": \"ref\", \"abs\": { \"min\": 90 } },\n"
        "    \"epics\": {\n"
        "      \"hallb_IPM2C21A_CUR\":  { \"abs\": { \"min\": 3 } },\n"
        "      \"hallb_IPM2C21A_XPOS\": { \"rel_rms\": 3 },\n"
        "      \"hallb_IPM2C21A_YPOS\": { \"rel_rms\": 3 }\n"
        "    },\n"
        "    \"split\": {\n"
        "      \"channel\": \"TGT:PRad:Cell_P\", \"full\": 500, \"empty\": 5,\n"
        "      \"guard_checkpoints\": 0, \"labels\": [\"full\", \"empty\"]\n"
        "    }\n"
        "  }\n"
        "  (split needs full > empty; guard_checkpoints drops an extra +/- N\n"
        "   points at each state edge, beyond the dead zone.)\n";
}

} // anonymous namespace

int main(int argc, char *argv[])
{
    std::vector<std::string> inputs;
    std::string output, cuts_path, report_path;
    int run_override = -1;
    int num_threads = 1;

    analysis::InitRootThreading();

    int opt;
    while ((opt = getopt(argc, argv, "o:c:j:r:t:h")) != -1) {
        switch (opt) {
        case 'o': output       = optarg;           break;
        case 'c': cuts_path    = optarg;           break;
        case 'j': report_path  = optarg;           break;
        case 'r': run_override = std::atoi(optarg); break;
        case 't': num_threads  = std::atoi(optarg); break;
        case 'h': usage(argv[0]); return 0;
        default:  usage(argv[0]); return 1;
        }
    }
    for (int i = optind; i < argc; ++i) inputs.push_back(argv[i]);

    if (inputs.empty() || output.empty() || cuts_path.empty()) {
        usage(argv[0]);
        return 1;
    }
    num_threads = std::max(1, num_threads);
    if (report_path.empty()) {
        if (inputs.size() > 1) {
            int rn = run_override >= 0 ? run_override : analysis::get_run_int(inputs.front());
            report_path = run_report_path(output, rn);
        } else {
            auto dot = output.rfind('.');
            report_path = (dot == std::string::npos)
                ? output + ".report.json"
                : output.substr(0, dot) + ".report.json";
        }
    }
    return run(inputs, output, cuts_path, report_path, run_override, num_threads);
}
