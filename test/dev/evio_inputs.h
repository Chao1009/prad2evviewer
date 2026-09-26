#pragma once
//=============================================================================
// evio_inputs.h — EVIO input discovery shared by the test/dev scaler tools
// (livetime, dsc_scan).  Unlike prad2::discover_split_files it does not
// require prad_<run> names.
//=============================================================================

#include <algorithm>
#include <filesystem>
#include <string>
#include <system_error>
#include <vector>

namespace devtools {

// Regular, non-hidden files in dir whose name starts with "<prefix>.".
inline std::vector<std::string> list_with_prefix(const std::filesystem::path &dir,
                                                 const std::string &prefix)
{
    namespace fs = std::filesystem;
    std::vector<std::string> out;
    std::error_code ec;
    const std::string lead = prefix + ".";
    for (const auto &e : fs::directory_iterator(dir, ec)) {
        const std::string name = e.path().filename().string();
        if (name.empty() || name[0] == '.') continue;
        if (name.compare(0, lead.size(), lead) != 0) continue;
        if (fs::is_regular_file(e.path(), ec)) out.push_back(e.path().string());
    }
    std::sort(out.begin(), out.end());
    return out;
}

// Files to process for a user-supplied input:
//   * split file "x.evio.00003"  -> every "x.evio.*" next to it
//   * other existing file         -> itself
//   * base name "x.evio" (absent) -> every "x.evio.*" in its directory
//   * directory                   -> every file whose path contains ".evio"
// Sorted by path; empty if nothing matches.
inline std::vector<std::string> discover_evio_inputs(const std::string &path)
{
    namespace fs = std::filesystem;
    std::error_code ec;
    const fs::path p(path);
    fs::path dir = p.parent_path();
    if (dir.empty()) dir = ".";

    if (fs::is_regular_file(p, ec)) {
        const std::string name = p.filename().string();
        const auto dot = name.rfind('.');
        if (dot != std::string::npos && dot + 1 < name.size() &&
            name.find_first_not_of("0123456789", dot + 1) == std::string::npos) {
            auto files = list_with_prefix(dir, name.substr(0, dot));
            if (!files.empty()) return files;
        }
        return { path };
    }

    if (!fs::is_directory(p, ec))
        return list_with_prefix(dir, p.filename().string());

    std::vector<std::string> files;
    for (const auto &e : fs::directory_iterator(p, ec)) {
        const std::string name = e.path().filename().string();
        if (name.empty() || name[0] == '.') continue;
        if (e.path().string().find(".evio") == std::string::npos) continue;
        if (fs::is_regular_file(e.path(), ec)) files.push_back(e.path().string());
    }
    std::sort(files.begin(), files.end());
    return files;
}

} // namespace devtools
