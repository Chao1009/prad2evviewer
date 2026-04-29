#pragma once
//============================================================================
// script_helpers.h — small utility functions shared between analysis ACLiC
// scripts (gem_hycal_matching.C, plot_hits_at_hycal.C, …).
//
// Why a header instead of `static` helpers per-script:
//   Cling shares its dictionary scope across all ACLiC-loaded .C files in
//   the same ROOT session.  `static` / anonymous-namespace helpers in two
//   scripts therefore collide with `redefinition of …` errors at the
//   second `.L`.  Marking the helpers `inline` here gives them weak
//   external linkage so the dict-payload merge accepts them.
//
//   Each script just `#include "script_helpers.h"` and the symbols are
//   shared.  Add a new helper here whenever a second script needs it;
//   keep one-script-only helpers private to that script.
//============================================================================

#include "DaqConfig.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cstdio>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <map>
#include <regex>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

// Resolve a possibly-relative database path to an absolute one using
// PRAD2_DATABASE_DIR.  Empty / already-absolute paths pass through.
inline std::string resolve_db_path(const std::string &p)
{
    if (p.empty()) return p;
    if (p[0] == '/' || p[0] == '\\') return p;
    if (p.size() >= 2 && p[1] == ':') return p;       // Windows drive letter
    const char *db = std::getenv("PRAD2_DATABASE_DIR");
    if (!db) return p;
    return std::string(db) + "/" + p;
}

// Sniff the run number out of a path like "prad_NNNNNN.evio.*".  Returns
// -1 if no plausible match is found, in which case LoadRunConfig falls
// back to the largest known runinfo entry (with a warning).
inline int extract_run_number_from_path(const std::string &path)
{
    static const std::regex pat(R"((?:prad|run)_0*(\d+))",
                                std::regex_constants::icase);
    std::smatch m;
    if (std::regex_search(path, m, pat)) {
        try { return std::stoi(m[1].str()); } catch (...) {}
    }
    return -1;
}

// Read database/reconstruction_config.json (under PRAD2_DATABASE_DIR or
// ./database) and return the resolved runinfo path, or "" if the pointer
// is missing / malformed.
inline std::string discover_runinfo_path()
{
    const char *db = std::getenv("PRAD2_DATABASE_DIR");
    std::string db_dir = db ? db : "database";
    std::ifstream f(db_dir + "/reconstruction_config.json");
    if (!f) return {};
    auto j = nlohmann::json::parse(f, nullptr, false, true);
    if (j.is_discarded() || !j.contains("runinfo") || !j["runinfo"].is_string())
        return {};
    return resolve_db_path(j["runinfo"].get<std::string>());
}

// EVIO bank-tag → logical-crate index (every ROC type).  Mirrors the
// roc_to_crate map in src/app_state_init.cpp / analysis/src/Replay.cpp.
inline std::map<int, int> build_full_crate_remap(const evc::DaqConfig &cfg)
{
    std::map<int, int> remap;
    for (const auto &re : cfg.roc_tags)
        remap[(int)re.tag] = re.crate;
    return remap;
}

// Same shape but only GEM ROCs — for GemSystem::LoadPedestals().
inline std::map<int, int> build_gem_crate_remap(const evc::DaqConfig &cfg)
{
    std::map<int, int> remap;
    for (const auto &re : cfg.roc_tags)
        if (re.type == "gem") remap[(int)re.tag] = re.crate;
    return remap;
}

// Resolve an EVIO input path to the list of files to process.
//
// Modes (chosen by the input path):
//   * Glob mode — path contains `*` (e.g. `.../prad_023881.evio.*`):
//       enumerate every sibling `prad_<run>.evio.<digits>` in the
//       enclosing directory, sort by suffix, and warn (to stderr) about
//       any gaps in the suffix sequence — including suffixes < the
//       lowest one found, since splits are expected to start at .00000.
//   * Directory mode — path is a directory:
//       same enumeration as glob mode, sniffing the run number from the
//       directory's name.
//   * Single-file mode — anything else:
//       return just `{ any_path }` unchanged.  Use this to process one
//       specific split (e.g. for debugging a single segment).
//
// File pattern: `prad_<run>.evio.<digits>`.  The run number in the name
// can be unpadded (`prad_1234.evio.0`) or zero-padded to any width
// (`prad_023881.evio.00000`); both forms are accepted on either side.
inline std::vector<std::string>
discover_split_files(const std::string &any_path)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    fs::path p(any_path);

    const bool wants_glob = (any_path.find('*') != std::string::npos);
    const bool is_dir     = fs::is_directory(p, ec);

    // Single-file mode: pass through unchanged.
    if (!wants_glob && !is_dir) return { any_path };

    // Discovery mode: figure out the search dir + run number.
    fs::path dir;
    int run = -1;
    if (is_dir) {
        dir = p;
        run = extract_run_number_from_path(p.filename().string());
    } else {
        // Glob: strip the glob suffix, work in the parent directory.
        dir = p.parent_path();
        if (dir.empty()) dir = ".";
        run = extract_run_number_from_path(p.filename().string());
        if (run < 0)
            run = extract_run_number_from_path(dir.filename().string());
    }
    if (run < 0 || !fs::is_directory(dir, ec)) {
        std::fprintf(stderr,
            "[WARN] discover_split_files: cannot resolve run/dir from '%s' — "
            "passing through as a single file.\n", any_path.c_str());
        return { any_path };
    }

    std::regex pat("^prad_0*" + std::to_string(run) + R"(\.evio\.(\d+)$)",
                   std::regex_constants::icase);

    // Collect (suffix_int, full_path) so we can sort numerically and detect
    // gaps in one pass.
    std::vector<std::pair<int, std::string>> matched;
    for (const auto &entry : fs::directory_iterator(dir, ec)) {
        std::string name = entry.path().filename().string();
        std::smatch m;
        if (std::regex_match(name, m, pat)) {
            try {
                matched.emplace_back(std::stoi(m[1].str()),
                                     entry.path().string());
            } catch (...) {}
        }
    }
    std::sort(matched.begin(), matched.end());

    // Gap warning: expected sequence is .00000, .00001, ..., contiguous.
    // Report missing suffixes from 0 to the highest found (so the user
    // notices both internal gaps AND a missing-from-the-start situation).
    if (!matched.empty()) {
        int last = matched.back().first;
        std::vector<int> missing;
        size_t k = 0;
        for (int i = 0; i <= last; ++i) {
            if (k < matched.size() && matched[k].first == i) { ++k; continue; }
            missing.push_back(i);
        }
        if (!missing.empty()) {
            std::fprintf(stderr,
                "[WARN] split-file gaps in run %d (found %zu file(s), "
                "max suffix .%05d): missing",
                run, matched.size(), last);
            for (int i : missing) std::fprintf(stderr, " .%05d", i);
            std::fprintf(stderr, "\n");
        }
    }

    std::vector<std::string> out;
    out.reserve(matched.size());
    for (auto &pr : matched) out.push_back(std::move(pr.second));
    if (out.empty()) {
        std::fprintf(stderr,
            "[WARN] discover_split_files: no files matched 'prad_%d.evio.*' "
            "in %s\n", run, dir.string().c_str());
        return { any_path };
    }
    return out;
}

// Strip the extension off a path so "out.pdf" becomes "out".  Used by
// scripts that derive a sibling .root output from a user-supplied
// canvas filename.  Leaves the directory alone.
inline std::string strip_extension(const std::string &p)
{
    auto dot = p.find_last_of('.');
    auto slash = p.find_last_of("/\\");
    if (dot == std::string::npos) return p;
    if (slash != std::string::npos && dot < slash) return p;
    return p.substr(0, dot);
}
