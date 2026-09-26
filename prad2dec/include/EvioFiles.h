#pragma once
//=============================================================================
// EvioFiles.h — run number and split-file discovery for prad_<run>.evio.<N>
//
// Header-only so ACLiC scripts can use it as well as compiled tools.
//=============================================================================

#include <algorithm>
#include <cstdio>
#include <filesystem>
#include <regex>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

namespace prad2
{

// Run number from a name such as "prad_023881.evio.00000" or "run_1234"
// (first "prad_"/"run_" followed by digits, case-insensitive; the directory
// part of the path is ignored).  -1 if there is none.
inline int run_number_from_path(const std::string &path)
{
    const auto slash = path.find_last_of("/\\");
    const std::string base = (slash == std::string::npos) ? path : path.substr(slash + 1);
    static const std::regex pat(R"((?:prad|run)_0*(\d+))", std::regex_constants::icase);
    std::smatch m;
    if (std::regex_search(base, m, pat)) {
        try { return std::stoi(m[1].str()); } catch (...) {}
    }
    return -1;
}

// Resolve an EVIO input path to the list of split files to process.
//
//   * Glob (path contains '*', e.g. ".../prad_023881.evio.*") or directory:
//       every prad_<run>.evio.<digits> in the directory, with the run taken
//       from the file name (glob) or the directory name.  The glob itself is
//       not matched; all splits of the run are returned.
//   * With expand_siblings, an existing split file ("prad_023881.evio.00003")
//       or a base name that is not itself a file ("prad_023881.evio") is
//       treated like a glob in its directory.
//   * Anything else (including other files named after a run, such as
//       "prad_023881.00000_raw.root") is returned unchanged as { path }.
//
// The run in file names may be unpadded or zero-padded.  Results are sorted
// by numeric suffix, and missing suffixes from .00000 to the highest found
// are reported on stderr.  When the run or directory cannot be resolved, or
// nothing matches, a warning is printed and the list is empty.
inline std::vector<std::string>
discover_split_files(const std::string &path, bool expand_siblings = false)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    const fs::path p(path);

    const bool wants_glob = (path.find('*') != std::string::npos);
    const bool is_dir     = fs::is_directory(p, ec);
    bool expand = false;
    if (expand_siblings && !wants_glob && !is_dir) {
        static const std::regex split_or_base(R"(^prad_0*\d+\.evio(\.\d+)?$)",
                                              std::regex_constants::icase);
        const std::string name = p.filename().string();
        std::smatch m;
        if (std::regex_match(name, m, split_or_base)) {
            const bool is_file = fs::is_regular_file(p, ec);
            expand = m[1].matched ? is_file : !is_file;
        }
    }
    if (!wants_glob && !is_dir && !expand)
        return { path };

    fs::path dir;
    int run = -1;
    if (is_dir) {
        dir = p;
        run = run_number_from_path(p.filename().string());
    } else {
        dir = p.parent_path();
        if (dir.empty()) dir = ".";
        run = run_number_from_path(p.filename().string());
        if (run < 0)
            run = run_number_from_path(dir.filename().string());
    }
    if (run < 0 || !fs::is_directory(dir, ec)) {
        std::fprintf(stderr,
            "[WARN] discover_split_files: cannot resolve run/dir from '%s'\n",
            path.c_str());
        return {};
    }

    const std::regex pat("^prad_0*" + std::to_string(run) + R"(\.evio\.(\d+)$)",
                         std::regex_constants::icase);
    std::vector<std::pair<int, std::string>> matched;
    for (const auto &entry : fs::directory_iterator(dir, ec)) {
        const std::string name = entry.path().filename().string();
        std::smatch m;
        if (std::regex_match(name, m, pat)) {
            try {
                matched.emplace_back(std::stoi(m[1].str()), entry.path().string());
            } catch (...) {}
        }
    }
    std::sort(matched.begin(), matched.end());

    if (matched.empty()) {
        std::fprintf(stderr,
            "[WARN] discover_split_files: no files matched 'prad_%d.evio.*' "
            "in %s\n", run, dir.string().c_str());
        return {};
    }

    std::vector<int> missing;
    const int last = matched.back().first;
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

    std::vector<std::string> out;
    out.reserve(matched.size());
    for (auto &pr : matched) out.push_back(std::move(pr.second));
    return out;
}

} // namespace prad2
