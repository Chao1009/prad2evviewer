#pragma once
//=============================================================================
// JsonUtil.h — nlohmann::json helpers shared by the JSON config loaders.
//
// Header-only; like load_daq_config.h it needs nlohmann/json on the
// include path, so include it only where nlohmann is already used.
//=============================================================================

#include <nlohmann/json.hpp>
#include <cstddef>
#include <fstream>
#include <string>

namespace prad2 {

// Read and parse a JSON file; // and /* */ comments are allowed.  Never
// throws: on failure `out` is untouched, false is returned and, if given,
// *err is set to "cannot open <path>" or "cannot parse <path>: <reason>".
inline bool read_json_file(const std::string &path, nlohmann::json &out,
                           std::string *err = nullptr)
{
    std::ifstream f(path);
    if (!f.is_open()) {
        if (err) *err = "cannot open " + path;
        return false;
    }
    try {
        out = nlohmann::json::parse(f, nullptr, true, /*ignore_comments=*/true);
    } catch (const nlohmann::json::exception &e) {
        if (err) *err = "cannot parse " + path + ": " + e.what();
        return false;
    }
    return true;
}

// Assign arr[0], arr[1], ... to out... in order via get<T>(); extra
// elements are ignored.  Returns false, leaving out... untouched, unless
// arr is an array with at least sizeof...(out) elements.
template <class... T>
inline bool read_json_elements(const nlohmann::json &arr, T &...out)
{
    if (!arr.is_array() || arr.size() < sizeof...(T)) return false;
    std::size_t i = 0;
    ((out = arr[i++].template get<T>()), ...);
    return true;
}

// read_json_elements(j[key], out...); false when j has no such key.
template <class... T>
inline bool read_json_array(const nlohmann::json &j, const char *key, T &...out)
{
    const auto it = j.find(key);
    return it != j.end() && read_json_elements(*it, out...);
}

} // namespace prad2
