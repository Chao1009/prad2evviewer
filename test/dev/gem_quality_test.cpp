//=============================================================================
// gem_quality_test.cpp — unit test for the SBS-style GEM quality variables
// and cuts in gem::GemCluster.  No data files: strips are hand-built; only
// the config-plumbing part reads the repo database (DATABASE_DIR), like
// hycal_energy_bias_test.
//
// Covers the three public helpers (StripMeanTime, TimeSampleCorrelation,
// IsUnimodalPulse), the StripCluster / GEMHit quality fields, every new
// ClusterConfig cut on its own, the unchanged mode-1 X/Y time cut, and the
// reconstruction_config.json parsing + [GEMCFG] log via PipelineBuilder.
//=============================================================================

#include "GemCluster.h"
#include "GemSystem.h"
#include "PipelineBuilder.h"

#include <nlohmann/json.hpp>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <random>
#include <sstream>
#include <string>
#include <unistd.h>
#include <utility>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;
using nlohmann::json;

namespace {

int failures = 0;

void check(bool condition, const std::string &message)
{
    if (condition) return;
    std::cerr << "FAIL: " << message << '\n';
    ++failures;
}

bool close_to(float actual, float expected, float tolerance = 1.e-4f)
{
    return std::fabs(actual - expected) <= tolerance;
}

bool same_bits(float a, float b)
{
    uint32_t ua, ub;
    std::memcpy(&ua, &a, sizeof ua);
    std::memcpy(&ub, &b, sizeof ub);
    return ua == ub;
}

// Pre-change seed mean time (GemCluster.cpp at f1a652f), kept verbatim as the
// reference for the bit-identity check of StripMeanTime.
float legacy_mean_time(const std::vector<float> &ts, float ts_period)
{
    if (ts.empty()) return -1.f;
    float sum_wt = 0.f, sum_w = 0.f;
    for (size_t i = 0; i < ts.size(); ++i) {
        float w = ts[i];
        if (w > 0.f) {
            sum_wt += w * static_cast<float>(i + 1) * ts_period;
            sum_w  += w;
        }
    }
    return (sum_w > 0.f) ? sum_wt / sum_w : -1.f;
}

// StripHit as GemSystem::collectHits builds it: charge = first max sample,
// max_timebin = its index, position = strip * 0.4 mm.
gem::StripHit make_strip(int strip, const std::vector<float> &ts)
{
    gem::StripHit h;
    h.strip = strip;
    h.ts_adc = ts;
    float mx = -1e9f;
    for (size_t i = 0; i < ts.size(); ++i)
        if (ts[i] > mx) { mx = ts[i]; h.max_timebin = static_cast<short>(i); }
    h.charge = mx;
    h.position = 0.4f * static_cast<float>(strip);
    return h;
}

// Pulse templates (6 samples, both strictly unimodal).  PROMPT: mean time
// 83.3 ns; LATE: 139.5 ns (peak in the last sample), r(PROMPT, LATE) = -0.55.
const std::vector<float> PROMPT = {0.1f, 0.5f, 1.0f, 0.7f, 0.3f, 0.1f};
const std::vector<float> LATE   = {0.005f, 0.01f, 0.02f, 0.05f, 0.4f, 1.0f};

std::vector<float> scaled(const std::vector<float> &shape, float amp)
{
    std::vector<float> v(shape);
    for (auto &x : v) x *= amp;
    return v;
}

// Summary of a cluster list: (first strip, size), sorted.
std::vector<std::pair<int, int>> summary(const std::vector<gem::StripCluster> &cls)
{
    std::vector<std::pair<int, int>> s;
    for (auto &c : cls)
        s.emplace_back(c.hits.empty() ? -1 : c.hits.front().strip,
                       static_cast<int>(c.hits.size()));
    std::sort(s.begin(), s.end());
    return s;
}

const gem::StripCluster *find_cluster(const std::vector<gem::StripCluster> &cls,
                                      int first_strip)
{
    for (auto &c : cls)
        if (!c.hits.empty() && c.hits.front().strip == first_strip) return &c;
    return nullptr;
}

// A hand-built plane with six clusters, each carrying one "bad" feature:
//   A 100-104  strip 101 has a tail bump (not unimodal), seed 102
//   B 200-203  seed 200 (index 0), strip 203 is LATE — after the seed, the
//              case SBS's early-stopping loop never checks (being late, it
//              is also poorly correlated with the seed)
//   C 300-302  strip 302 oscillates: compatible mean time, poor correlation
//   D 400      single strip, LATE (out of a [25,110] ns window)
//   E 500-502  seed peak 25 < 30 (sum 67.5 >= 60)
//   F 600-601  seed peak 40 >= 30 but sum 35 < 60 (negative samples count)
std::vector<gem::StripHit> make_plane()
{
    std::vector<gem::StripHit> h;
    // A
    h.push_back(make_strip(100, scaled(PROMPT, 150.f)));
    h.push_back(make_strip(101, {30.f, 150.f, 300.f, 150.f, 30.f, 60.f}));
    h.push_back(make_strip(102, scaled(PROMPT, 400.f)));
    h.push_back(make_strip(103, scaled(PROMPT, 250.f)));
    h.push_back(make_strip(104, scaled(PROMPT, 100.f)));
    // B
    h.push_back(make_strip(200, scaled(PROMPT, 500.f)));
    h.push_back(make_strip(201, scaled(PROMPT, 400.f)));
    h.push_back(make_strip(202, scaled(PROMPT, 300.f)));
    h.push_back(make_strip(203, scaled(LATE,   250.f)));
    // C
    h.push_back(make_strip(300, scaled(PROMPT, 200.f)));
    h.push_back(make_strip(301, scaled(PROMPT, 350.f)));
    h.push_back(make_strip(302, {150.f, 20.f, 150.f, 20.f, 150.f, 20.f}));
    // D
    h.push_back(make_strip(400, scaled(LATE, 300.f)));
    // E
    h.push_back(make_strip(500, scaled(PROMPT, 20.f)));
    h.push_back(make_strip(501, scaled(PROMPT, 25.f)));
    h.push_back(make_strip(502, scaled(PROMPT, 20.f)));
    // F
    h.push_back(make_strip(600, {-10.f, 5.f, 40.f, 10.f, -4.f, -6.f}));
    h.push_back(make_strip(601, {-10.f, 5.f, 30.f, 10.f, -4.f, -6.f}));
    // shuffle so FormClusters has to sort
    std::reverse(h.begin(), h.end());
    return h;
}

// Library-default config with cross-talk removal explicitly off.
gem::ClusterConfig base_config()
{
    gem::ClusterConfig cfg;
    cfg.charac_dists.clear();
    return cfg;
}

std::vector<gem::StripCluster> cluster_plane(const gem::ClusterConfig &cfg,
                                             std::vector<gem::StripHit> *hits_out = nullptr)
{
    gem::GemCluster gc;
    gc.SetConfig(cfg);
    auto hits = make_plane();
    std::vector<gem::StripCluster> cls;
    gc.FormClusters(hits, cls);
    if (hits_out) *hits_out = hits;
    return cls;
}

bool clusters_identical(const std::vector<gem::StripCluster> &a,
                        const std::vector<gem::StripCluster> &b)
{
    if (a.size() != b.size()) return false;
    for (size_t i = 0; i < a.size(); ++i) {
        const auto &x = a[i], &y = b[i];
        if (!same_bits(x.position, y.position) || !same_bits(x.peak_charge, y.peak_charge) ||
            !same_bits(x.total_charge, y.total_charge) || x.max_timebin != y.max_timebin ||
            x.cross_talk != y.cross_talk || x.hits.size() != y.hits.size() ||
            !same_bits(x.seed_time, y.seed_time) || !same_bits(x.seed_peak_adc, y.seed_peak_adc) ||
            !same_bits(x.seed_sum_adc, y.seed_sum_adc) ||
            !same_bits(x.max_strip_dt, y.max_strip_dt) || !same_bits(x.min_ts_corr, y.min_ts_corr))
            return false;
        for (size_t k = 0; k < x.hits.size(); ++k)
            if (x.hits[k].strip != y.hits[k].strip ||
                !same_bits(x.hits[k].charge, y.hits[k].charge)) return false;
    }
    return true;
}

// Builds a pipeline from `recon_path` and returns the [GEMCFG] d0 line.
std::string gemcfg_line(const std::string &recon_path,
                        std::vector<gem::ClusterConfig> *cfgs_out)
{
    std::ostringstream log;
    auto p = prad2::PipelineBuilder()
        .set_database_dir(DATABASE_DIR)
        .set_recon_config(recon_path)
        .set_run_number(24560)
        .set_log_stream(&log)
        .set_log_pedestal_checksum(false)
        .build();
    if (cfgs_out) *cfgs_out = p.gem.GetReconConfigs();
    std::istringstream in(log.str());
    std::string line;
    while (std::getline(in, line))
        if (line.rfind("[GEMCFG] d0", 0) == 0) return line;
    return {};
}

bool knobs_disabled(const gem::ClusterConfig &c)
{
    return std::isinf(c.strip_time_min) && c.strip_time_min < 0.f &&
           std::isinf(c.strip_time_max) && c.strip_time_max > 0.f &&
           !c.strip_unimodal && c.seed_min_peak_adc == 0.f &&
           c.seed_min_sum_adc == 0.f && c.strip_time_agreement == -1.f &&
           c.strip_ts_corr_min == -1.f;
}

} // namespace

int main()
{
    // ---------------------------------------------------------------------
    // StripMeanTime
    // ---------------------------------------------------------------------
    check(close_to(gem::StripMeanTime({0, 100, 200, 100, 0, 0}, 25.f), 75.f),
          "StripMeanTime {0,100,200,100,0,0} @25 = 75 ns");
    check(close_to(gem::StripMeanTime({0, 100, 200, 100, 0, 0}), 75.f),
          "StripMeanTime default ts_period is 25 ns");
    check(close_to(gem::StripMeanTime({0, 100, 200, 100, 0, 0}, 10.f), 30.f),
          "StripMeanTime scales with ts_period");
    check(close_to(gem::StripMeanTime({-50, 100, 200, 100, -30, -20}, 25.f), 75.f),
          "StripMeanTime ignores negative samples");
    check(close_to(gem::StripMeanTime({100, 0, 0, 0, 0, 0}, 25.f), 25.f),
          "StripMeanTime sample 0 maps to 1*ts_period");
    check(std::isnan(gem::StripMeanTime({0, -1, -2, 0, -3, 0}, 25.f)),
          "StripMeanTime with no positive sample is NaN");
    check(std::isnan(gem::StripMeanTime({}, 25.f)), "StripMeanTime of empty is NaN");
    {
        // bit-identity with the pre-change seed mean time
        std::mt19937 rng(12345);
        std::uniform_real_distribution<float> adc(-60.f, 1500.f);
        std::uniform_int_distribution<int> len(0, 9);
        int mismatches = 0;
        for (int n = 0; n < 20000; ++n) {
            std::vector<float> ts(len(rng));
            for (auto &x : ts) x = adc(rng);
            if (n % 7 == 0) for (auto &x : ts) x = -std::fabs(x);  // no positive sample
            float t_new = gem::StripMeanTime(ts, 25.f);
            float t_old = legacy_mean_time(ts, 25.f);
            bool ok = (t_old < 0.f) ? std::isnan(t_new) : same_bits(t_new, t_old);
            if (!ok) ++mismatches;
        }
        check(mismatches == 0, "StripMeanTime is bit-identical to the legacy seed mean time");
    }

    // ---------------------------------------------------------------------
    // TimeSampleCorrelation
    // ---------------------------------------------------------------------
    const std::vector<float> pulse = {10, 50, 100, 60, 20, 5};
    std::vector<float> reversed(pulse.rbegin(), pulse.rend());
    check(close_to(gem::TimeSampleCorrelation(pulse, pulse), 1.f, 1e-6f),
          "TimeSampleCorrelation identical = 1");
    check(close_to(gem::TimeSampleCorrelation(pulse, scaled(pulse, 3.7f)), 1.f, 1e-6f),
          "TimeSampleCorrelation scaled copy = 1");
    check(close_to(gem::TimeSampleCorrelation(pulse, scaled(pulse, -2.f)), -1.f, 1e-6f),
          "TimeSampleCorrelation negated copy = -1");
    check(gem::TimeSampleCorrelation(pulse, reversed) < 0.99f,
          "TimeSampleCorrelation reversed pulse < 1");
    check(std::isnan(gem::TimeSampleCorrelation(pulse, {7, 7, 7, 7, 7, 7})),
          "TimeSampleCorrelation flat = NaN");
    check(std::isnan(gem::TimeSampleCorrelation(pulse, {10, 50, 100, 60, 20})),
          "TimeSampleCorrelation size mismatch = NaN");
    check(std::isnan(gem::TimeSampleCorrelation({1.f}, {2.f})),
          "TimeSampleCorrelation size < 2 = NaN");

    // ---------------------------------------------------------------------
    // IsUnimodalPulse
    // ---------------------------------------------------------------------
    check(gem::IsUnimodalPulse(pulse), "IsUnimodalPulse rising/falling = true");
    check(!gem::IsUnimodalPulse({0, 100, 200, 100, 20, 25}), "IsUnimodalPulse tail bump = false");
    check(!gem::IsUnimodalPulse({10, 50, 100, 100, 20, 5}), "IsUnimodalPulse flat pair at peak = false");
    check(!gem::IsUnimodalPulse({10, 50, 50, 100, 20, 5}), "IsUnimodalPulse flat pair on rise = false");
    check(gem::IsUnimodalPulse({200, 100, 50, 20, 10, 5}), "IsUnimodalPulse peak at sample 0 = true");
    check(gem::IsUnimodalPulse({1, 2, 3, 4, 5, 6}), "IsUnimodalPulse peak at last sample = true");
    check(!gem::IsUnimodalPulse({}), "IsUnimodalPulse empty = false");

    // ---------------------------------------------------------------------
    // Default config: no new cut acts, quality fields filled
    // ---------------------------------------------------------------------
    const std::vector<std::pair<int, int>> all6 =
        {{100, 5}, {200, 4}, {300, 3}, {400, 1}, {500, 3}, {600, 2}};
    std::vector<gem::StripHit> hits_after;
    auto def = cluster_plane(base_config(), &hits_after);
    check(summary(def) == all6, "default config keeps all six clusters");
    check(hits_after.size() == 18, "FormClusters keeps every plane hit");
    {
        // explicitly-disabled values (every "off" spelling) change nothing
        gem::ClusterConfig off = base_config();
        off.strip_time_min = -std::numeric_limits<float>::infinity();
        off.strip_time_max =  std::numeric_limits<float>::infinity();
        off.strip_unimodal = false;
        off.seed_min_peak_adc = -5.f;
        off.seed_min_sum_adc = 0.f;
        off.strip_time_agreement = -0.5f;
        off.strip_ts_corr_min = -1.5f;
        check(clusters_identical(def, cluster_plane(off)),
              "explicitly disabled cuts give bit-identical clusters");
    }
    if (const auto *a = find_cluster(def, 100)) {
        check(close_to(a->seed_time, gem::StripMeanTime(scaled(PROMPT, 400.f))),
              "A seed_time is the mean time of strip 102");
        check(close_to(a->seed_peak_adc, 400.f), "A seed_peak_adc = 400");
        check(close_to(a->seed_sum_adc, 1080.f, 1e-2f), "A seed_sum_adc = 2.7*400");
        check(std::isfinite(a->max_strip_dt) && std::isfinite(a->min_ts_corr),
              "multi-strip cluster has finite dt / corr");
    } else check(false, "cluster A present");
    if (const auto *b = find_cluster(def, 200)) {
        float dt_late = gem::StripMeanTime(scaled(LATE, 250.f)) -
                        gem::StripMeanTime(scaled(PROMPT, 500.f));
        check(close_to(b->max_strip_dt, dt_late, 1e-3f),
              "B max_strip_dt sees the late strip AFTER the seed");
        check(b->max_strip_dt > 50.f, "B max_strip_dt > 50 ns");
    } else check(false, "cluster B present");
    if (const auto *c = find_cluster(def, 300)) {
        check(c->min_ts_corr < 0.7f, "C min_ts_corr < 0.7 (oscillating strip)");
        check(c->max_strip_dt < 10.f, "C oscillating strip has a compatible mean time");
    } else check(false, "cluster C present");
    if (const auto *d = find_cluster(def, 400)) {
        check(std::isfinite(d->seed_time), "single-strip seed_time finite");
        check(std::isnan(d->max_strip_dt), "single-strip max_strip_dt NaN");
        check(std::isnan(d->min_ts_corr), "single-strip min_ts_corr NaN");
    } else check(false, "cluster D present");
    if (const auto *f = find_cluster(def, 600)) {
        check(close_to(f->seed_peak_adc, 40.f), "F seed_peak_adc = 40");
        check(close_to(f->seed_sum_adc, 35.f), "F seed_sum_adc includes negative samples");
    } else check(false, "cluster F present");

    // ---------------------------------------------------------------------
    // Each new cut on its own
    // ---------------------------------------------------------------------
    {
        gem::ClusterConfig cfg = base_config();
        cfg.strip_unimodal = true;
        std::vector<gem::StripHit> hits;
        auto cls = cluster_plane(cfg, &hits);
        // A splits around strip 101; C loses the oscillating strip 302
        const std::vector<std::pair<int, int>> expect =
            {{100, 1}, {102, 3}, {200, 4}, {300, 2}, {400, 1}, {500, 3}, {600, 2}};
        check(summary(cls) == expect, "strip_unimodal drops strip 101 / 302 and splits A");
        check(hits.size() == 18, "strip cut does not erase plane hits");
    }
    {
        gem::ClusterConfig cfg = base_config();
        cfg.strip_time_min = 25.f;
        cfg.strip_time_max = 110.f;
        // late strips 203 (end of B) and 400 (D) fail the window
        const std::vector<std::pair<int, int>> expect =
            {{100, 5}, {200, 3}, {300, 3}, {500, 3}, {600, 2}};
        check(summary(cluster_plane(cfg)) == expect, "strip_mean_time_range drops late strips");

        // one-sided window, and a strip without positive sample fails it
        cfg.strip_time_min = -std::numeric_limits<float>::infinity();
        gem::GemCluster gc;
        gc.SetConfig(cfg);
        std::vector<gem::StripHit> hits = {
            make_strip(10, scaled(PROMPT, 100.f)),
            make_strip(11, {-5.f, -3.f, 0.f, -1.f, -2.f, -4.f}),
            make_strip(12, scaled(PROMPT, 90.f))};
        hits[1].charge = 50.f;
        std::vector<gem::StripCluster> cls;
        gc.FormClusters(hits, cls);
        check(summary(cls) == std::vector<std::pair<int, int>>{{10, 1}, {12, 1}},
              "strip without positive sample fails an active time window");
    }
    {
        gem::ClusterConfig cfg = base_config();
        cfg.seed_min_peak_adc = 30.f;
        const std::vector<std::pair<int, int>> expect =
            {{100, 5}, {200, 4}, {300, 3}, {400, 1}, {600, 2}};
        check(summary(cluster_plane(cfg)) == expect, "seed_min_peak_adc rejects E");
    }
    {
        gem::ClusterConfig cfg = base_config();
        cfg.seed_min_sum_adc = 60.f;
        const std::vector<std::pair<int, int>> expect =
            {{100, 5}, {200, 4}, {300, 3}, {400, 1}, {500, 3}};
        check(summary(cluster_plane(cfg)) == expect, "seed_min_sum_adc rejects F");
    }
    {
        gem::ClusterConfig cfg = base_config();
        cfg.strip_time_agreement = 50.f;
        const std::vector<std::pair<int, int>> expect =
            {{100, 5}, {300, 3}, {400, 1}, {500, 3}, {600, 2}};
        check(summary(cluster_plane(cfg)) == expect,
              "strip_time_agreement rejects B (late strip after the seed)");
    }
    {
        gem::ClusterConfig cfg = base_config();
        cfg.strip_ts_corr_min = 0.7f;
        const std::vector<std::pair<int, int>> expect =
            {{100, 5}, {400, 1}, {500, 3}, {600, 2}};
        auto cls = cluster_plane(cfg);
        check(summary(cls) == expect, "strip_ts_corr_min rejects B and C");
        check(find_cluster(cls, 400) != nullptr, "NaN min_ts_corr (single strip) passes");
    }

    // ---------------------------------------------------------------------
    // makeHit fields (both match modes) and the mode-1 time cut
    // ---------------------------------------------------------------------
    {
        gem::GemCluster gc;
        gem::ClusterConfig cfg = base_config();
        gc.SetConfig(cfg);
        std::vector<gem::StripHit> xh = {make_strip(10, scaled(PROMPT, 100.f)),
                                         make_strip(11, scaled(PROMPT, 300.f)),
                                         make_strip(12, scaled(PROMPT, 150.f))};
        // Y pulse ~8 ns later than X (inside match_time_diff)
        const std::vector<float> y_shape = {0.05f, 0.3f, 0.8f, 1.0f, 0.4f, 0.1f};
        std::vector<gem::StripHit> yh = {make_strip(50, scaled(y_shape, 200.f)),
                                         make_strip(51, scaled(y_shape, 120.f))};
        std::vector<gem::StripCluster> xc, yc;
        gc.FormClusters(xh, xc);
        gc.FormClusters(yh, yc);
        for (int mode = 0; mode < 2; ++mode) {
            cfg.match_mode = mode;
            gc.SetConfig(cfg);
            std::vector<gem::GEMHit> out;
            gc.CartesianReconstruct(xc, yc, out, 2);
            const std::string tag = " (mode " + std::to_string(mode) + ")";
            check(out.size() == 1, "one X/Y pair" + tag);
            if (out.size() != 1) continue;
            const auto &h = out[0];
            check(same_bits(h.x_time, xc[0].seed_time) && same_bits(h.y_time, yc[0].seed_time),
                  "x_time / y_time copy the seed times" + tag);
            check(same_bits(h.time_diff, h.x_time - h.y_time) && h.time_diff < -5.f,
                  "time_diff = x_time - y_time" + tag);
            check(close_to(h.adc_asym, (300.f - 200.f) / 500.f), "adc_asym signed (+0.2)" + tag);
            check(same_bits(h.x_max_strip_dt, xc[0].max_strip_dt) &&
                  same_bits(h.y_max_strip_dt, yc[0].max_strip_dt) &&
                  same_bits(h.x_min_ts_corr, xc[0].min_ts_corr) &&
                  same_bits(h.y_min_ts_corr, yc[0].min_ts_corr),
                  "cluster quality copied onto GEMHit" + tag);
            // swapped planes flip the sign of both signed quantities
            std::vector<gem::GEMHit> sw;
            gc.CartesianReconstruct(yc, xc, sw, 2);
            check(sw.size() == 1 && close_to(sw[0].adc_asym, -0.2f) &&
                  sw[0].time_diff == -h.time_diff,
                  "adc_asym / time_diff are signed" + tag);
        }
    }
    {
        gem::GemCluster gc;
        gem::ClusterConfig cfg = base_config();
        cfg.match_mode = 1;
        cfg.match_adc_asymmetry = 0.8f;
        cfg.match_time_diff = 50.f;
        gc.SetConfig(cfg);

        std::vector<gem::StripHit> xh = {make_strip(10, scaled(PROMPT, 300.f))};
        std::vector<gem::StripHit> yh_prompt = {make_strip(50, scaled(PROMPT, 280.f))};
        std::vector<gem::StripHit> yh_late   = {make_strip(50, scaled(LATE, 280.f))};
        std::vector<gem::StripHit> yh_undef  = {make_strip(50, {-5.f, -3.f, 0.f, -1.f, -2.f, -4.f})};
        yh_undef[0].charge = 280.f;   // no positive sample, but a valid charge
        std::vector<gem::StripCluster> xc, yp, yl, yu;
        gc.FormClusters(xh, xc);
        gc.FormClusters(yh_prompt, yp);
        gc.FormClusters(yh_late, yl);
        gc.FormClusters(yh_undef, yu);

        std::vector<gem::GEMHit> out;
        gc.CartesianReconstruct(xc, yp, out, 0);
        check(out.size() == 1, "mode 1: in-time pair accepted");
        gc.CartesianReconstruct(xc, yl, out, 0);
        check(std::fabs(xc[0].seed_time - yl[0].seed_time) > 50.f && out.empty(),
              "mode 1: pair with |dt| > match_time_diff rejected");
        gc.CartesianReconstruct(xc, yu, out, 0);
        check(out.size() == 1 && std::isnan(out[0].y_time) && std::isnan(out[0].time_diff),
              "mode 1: seed without positive sample is not rejected by the time cut");

        cfg.match_time_diff = -1.f;   // time cut off → late pair accepted
        gc.SetConfig(cfg);
        gc.CartesianReconstruct(xc, yl, out, 0);
        check(out.size() == 1 && out[0].time_diff < -50.f,
              "mode 1: time cut disabled keeps the late pair with its time_diff");

        // asymmetry sum <= 0 → NaN adc_asym, cut skipped as before
        std::vector<gem::StripCluster> zx(1), zy(1);
        zx[0].peak_charge = 0.f;
        zy[0].peak_charge = 0.f;
        gc.CartesianReconstruct(zx, zy, out, 0);
        check(out.size() == 1 && std::isnan(out[0].adc_asym) && std::isnan(out[0].x_time),
              "hand-built clusters: NaN adc_asym / times, pair kept");
    }

    // ---------------------------------------------------------------------
    // Config plumbing: shipped configs parse with every new cut disabled
    // ---------------------------------------------------------------------
    const std::string expect_log =
        " strip_t=[-inf,inf] unimodal=0 seed_peak=0 seed_sum=0 strip_dt=-1 ts_corr=-1";
    for (const char *name : {"reconstruction_config.json", "reconstruction_config_x17.json"}) {
        std::vector<gem::ClusterConfig> cfgs;
        std::string line;
        try {
            line = gemcfg_line(std::string(DATABASE_DIR) + "/" + name, &cfgs);
        } catch (const std::exception &e) {
            check(false, std::string(name) + ": PipelineBuilder threw: " + e.what());
            continue;
        }
        std::cout << name << ": " << line << '\n';
        check(line.size() >= expect_log.size() &&
              line.compare(line.size() - expect_log.size(), expect_log.size(), expect_log) == 0,
              std::string(name) + ": [GEMCFG] shows the new cuts disabled");
        check(!cfgs.empty(), std::string(name) + ": per-detector configs installed");
        for (auto &c : cfgs)
            check(knobs_disabled(c), std::string(name) + ": new knobs at library defaults");
    }
    {
        // explicit values + per-detector override; bad shapes are ignored
        const fs::path dir = fs::temp_directory_path() /
            ("prad2_gem_quality_" + std::to_string(static_cast<long long>(::getpid())));
        fs::create_directories(dir);
        const fs::path path = dir / "recon.json";
        json root;
        auto &d = root["gem"]["default"];
        d["strip_mean_time_range"] = json::array({25.0, 150.0});
        d["strip_unimodal_shape"]  = "yes";
        d["seed_min_peak_adc"]     = 30.0;
        d["seed_min_sum_adc"]      = 60;      // integer JSON value
        d["strip_time_agreement"]  = 50.0;
        d["strip_ts_corr_min"]     = 0.7;
        root["gem"]["1"]["strip_mean_time_range"] = json::array();
        root["gem"]["1"]["strip_unimodal_shape"]  = false;
        root["gem"]["2"]["strip_mean_time_range"] = nullptr;
        root["gem"]["3"]["strip_mean_time_range"] = json::array({10.0});
        std::ofstream(path) << root.dump(2) << '\n';

        std::vector<gem::ClusterConfig> cfgs;
        std::string line;
        try {
            line = gemcfg_line(path.string(), &cfgs);
        } catch (const std::exception &e) {
            check(false, std::string("fixture: PipelineBuilder threw: ") + e.what());
        }
        std::cout << "fixture: " << line << '\n';
        check(line.find(" strip_t=[25,150] unimodal=1 seed_peak=30 seed_sum=60"
                        " strip_dt=50 ts_corr=0.7") != std::string::npos,
              "fixture: [GEMCFG] d0 shows the parsed values");
        if (cfgs.size() >= 4) {
            const auto &c0 = cfgs[0];
            check(c0.strip_time_min == 25.f && c0.strip_time_max == 150.f && c0.strip_unimodal &&
                  c0.seed_min_peak_adc == 30.f && c0.seed_min_sum_adc == 60.f &&
                  c0.strip_time_agreement == 50.f && close_to(c0.strip_ts_corr_min, 0.7f, 1e-6f),
                  "fixture: default block parsed");
            check(std::isinf(cfgs[1].strip_time_min) && std::isinf(cfgs[1].strip_time_max) &&
                  !cfgs[1].strip_unimodal && cfgs[1].seed_min_peak_adc == 30.f,
                  "fixture: det 1 override [] disables the window, keeps the rest");
            check(std::isinf(cfgs[2].strip_time_min) && std::isinf(cfgs[2].strip_time_max),
                  "fixture: det 2 null disables the window");
            check(cfgs[3].strip_time_min == 25.f && cfgs[3].strip_time_max == 150.f,
                  "fixture: det 3 malformed range ignored");
        } else check(false, "fixture: four detector configs");
        fs::remove_all(dir);
    }

    if (failures != 0) {
        std::cerr << failures << " GEM quality test(s) failed\n";
        return 1;
    }
    std::cout << "GEM quality tests passed\n";
    return 0;
}
