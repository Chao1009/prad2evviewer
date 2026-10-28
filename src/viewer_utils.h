#pragma once

#include "DetectorTransform.h"
#include "Fadc250Data.h"
#include "WaveAnalyzer.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cctype>
#include <fstream>
#include <string>
#include <string_view>
#include <vector>
#include <cmath>

// --- TI timestamp conversion ------------------------------------------------
// Safe (now − base) seconds for uint64 TI ticks: a plain unsigned difference
// with now < base wraps to ~2^64 × 4 ns ≈ 7.4e10 s.  now == 0 means "no anchor
// yet" and returns 0; now < base (ET reconnect, out-of-order events) returns
// an honest negative.
inline double ti_delta_sec(uint64_t now, uint64_t base) {
    using fdec::TI_TICK_SEC;
    if (now == 0) return 0.0;
    if (base == 0) return static_cast<double>(now) * TI_TICK_SEC;
    return (now >= base)
        ?  static_cast<double>(now  - base) * TI_TICK_SEC
        : -static_cast<double>(base - now ) * TI_TICK_SEC;
}

// --- file I/O helpers -------------------------------------------------------
inline std::string readFile(const std::string &path) {
    std::ifstream f(path);
    if (!f) return "";
    return {std::istreambuf_iterator<char>(f), {}};
}

inline std::string findFile(const std::string &name, const std::string &base) {
    { std::ifstream f(name); if (f.good()) return name; }
    std::string p = base + "/" + name;
    { std::ifstream f(p); if (f.good()) return p; }
    return "";
}

inline std::string contentType(const std::string &path) {
    if (path.size() >= 5 && path.substr(path.size()-5) == ".html") return "text/html; charset=utf-8";
    if (path.size() >= 4 && path.substr(path.size()-4) == ".css")  return "text/css; charset=utf-8";
    if (path.size() >= 3 && path.substr(path.size()-3) == ".js")   return "application/javascript; charset=utf-8";
    return "application/octet-stream";
}

// --- URL query helpers ------------------------------------------------------
// Percent-decode s.  A '%' not followed by two hex digits is kept as is.
// '+' means a space in a query string but not in a path segment.
inline std::string urlDecode(std::string_view s, bool plus_as_space = true)
{
    auto hex = [](char c) {
        return std::isdigit((unsigned char)c) ? c - '0'
                                              : std::tolower((unsigned char)c) - 'a' + 10;
    };
    std::string out;
    out.reserve(s.size());
    for (size_t i = 0; i < s.size(); ++i) {
        if (s[i] == '%' && i + 2 < s.size()
            && std::isxdigit((unsigned char)s[i + 1])
            && std::isxdigit((unsigned char)s[i + 2])) {
            out += (char)(hex(s[i + 1]) * 16 + hex(s[i + 2]));
            i += 2;
        } else if (s[i] == '+' && plus_as_space) {
            out += ' ';
        } else {
            out += s[i];
        }
    }
    return out;
}

// Decoded values of every `key=value` pair in the query part of uri.
inline std::vector<std::string> queryValues(std::string_view uri, std::string_view key)
{
    std::vector<std::string> values;
    auto q = uri.find('?');
    if (q == std::string_view::npos) return values;
    std::string_view query = uri.substr(q + 1);
    for (size_t pos = 0; pos < query.size();) {
        size_t amp = query.find('&', pos);
        if (amp == std::string_view::npos) amp = query.size();
        std::string_view kv = query.substr(pos, amp - pos);
        if (kv.size() > key.size() && kv.compare(0, key.size(), key) == 0
            && kv[key.size()] == '=')
            values.push_back(urlDecode(kv.substr(key.size() + 1)));
        pos = amp + 1;
    }
    return values;
}

// First decoded value of `key` in the query part of uri, or def.
inline std::string queryValue(std::string_view uri, std::string_view key,
                              std::string def = {})
{
    auto values = queryValues(uri, key);
    return values.empty() ? def : values.front();
}

// --- Histogram --------------------------------------------------------------
// Uniform binning from min in steps of step; configured in
// monitor_config.json as {min, max, step} (or with key prefixes such as
// angle_min / x_min when one block holds several axes).  With T = int the
// bounds stay integers (fractional JSON values truncate) and nbins() floors
// (max - min) / step.
template <typename T>
struct BasicHistAxis {
    T min = 0, max = 1, step = 1;
    int nbins() const { return std::max(1, (int)std::ceil((max - min) / step)); }
    // Reads <pfx>min / <pfx>max / <pfx>step where present.
    void parse(const nlohmann::json &j, const std::string &pfx = "") {
        if (j.contains(pfx + "min"))  min  = j[pfx + "min"];
        if (j.contains(pfx + "max"))  max  = j[pfx + "max"];
        if (j.contains(pfx + "step")) step = j[pfx + "step"];
    }
    // Writes <pfx>min / <pfx>max / <pfx>step into out.
    void toJson(nlohmann::json &out, const std::string &pfx) const {
        out[pfx + "min"] = min; out[pfx + "max"] = max; out[pfx + "step"] = step;
    }
    nlohmann::json toJson() const { return {{"min", min}, {"max", max}, {"step", step}}; }
};
using HistAxis    = BasicHistAxis<float>;
using IntHistAxis = BasicHistAxis<int>;

struct Histogram {
    int underflow = 0, overflow = 0;
    std::vector<int> bins;
    void init(int n) { bins.assign(n, 0); underflow = overflow = 0; }
    template <typename T>
    void init(const BasicHistAxis<T> &a) { init(a.nbins()); }
    void fill(float v, float bmin, float bstep) {
        if (v < bmin) { ++underflow; return; }
        int b = (int)((v - bmin) / bstep);
        if (b >= (int)bins.size()) { ++overflow; return; }
        ++bins[b];
    }
    template <typename T>
    void fill(float v, const BasicHistAxis<T> &a) { fill(v, a.min, a.step); }
    void clear() { std::fill(bins.begin(), bins.end(), 0); underflow = overflow = 0; }
};

struct Histogram2D {
    int nx = 0, ny = 0;
    std::vector<int> bins;  // row-major: bins[iy*nx + ix]
    void init(int nx_, int ny_) { nx = nx_; ny = ny_; bins.assign(nx * ny, 0); }
    void init(const HistAxis &ax, const HistAxis &ay) { init(ax.nbins(), ay.nbins()); }
    void fill(float vx, float vy, float xmin, float xstep, float ymin, float ystep) {
        int ix = (int)((vx - xmin) / xstep);
        int iy = (int)((vy - ymin) / ystep);
        if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) return;
        bins[iy * nx + ix]++;
    }
    void fill(float vx, float vy, const HistAxis &ax, const HistAxis &ay) {
        fill(vx, vy, ax.min, ax.step, ay.min, ay.step);
    }
    void clear() { std::fill(bins.begin(), bins.end(), 0); }
};

// --- Histogram config -------------------------------------------------------
// Pure binning info — peak-detection thresholds live in WaveConfig
// (daq_config.json fadc250_waveform.analyzer); per-tab Waveform-Tab cuts
// live in PeakFilter (monitor_config.json waveform.filter).
struct HistConfig {
    HistAxis integral{0.f, 20000.f, 100.f};
    HistAxis time{0.f, 400.f, 4.f};
    HistAxis height{0.f, 4000.f, 10.f};
};

// --- LMS entry --------------------------------------------------------------
struct LmsEntry {
    double time_sec;    // seconds since first LMS event (from TI timestamp)
    float  integral;    // peak integral within timing cut (or raw ADC for ADC1881M)
};

// --- Peak extraction helpers ------------------------------------------------
// Peaks in `wres` are already gated by the analyzer's
// max(peak_nsigma·ped.rms, min_peak_height) detection cut, so no extra
// height filter is needed in these picks.

// Best peak integral within a time window. Returns -1 if no peak qualifies.
inline float bestPeakInWindow(const fdec::WaveResult &wres,
                               float time_min, float time_max)
{
    float best = -1;
    for (int p = 0; p < wres.npeaks; ++p) {
        auto &pk = wres.peaks[p];
        if (pk.time >= time_min && pk.time <= time_max)
            if (pk.integral > best) best = pk.integral;
    }
    return best;
}

// Best peak integral across all detected peaks (no time cut) — the
// single-pulse clustering input.
inline float bestPeak(const fdec::WaveResult &wres)
{
    float best = -1;
    for (int p = 0; p < wres.npeaks; ++p)
        if (wres.peaks[p].integral > best) best = wres.peaks[p].integral;
    return best;
}
