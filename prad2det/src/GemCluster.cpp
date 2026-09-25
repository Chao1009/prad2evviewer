//============================================================================//
// GEM Cluster Reconstruction                                                 //
//                                                                            //
// Ported from mpd_gem_view_ssp GEMCluster                                    //
// Original authors: Xinzhan Bai, Kondo Gnanvo, Chao Peng                     //
//                                                                            //
// Also computes the SBS-style (mpd_gem_view_ssp Cuts) quality variables —    //
// seed mean time / peak / sum, seed-vs-strip time spread and time-sample     //
// correlation, X/Y time difference and ADC asymmetry — and applies the       //
// matching optional cuts (all off by default).  SBS bugs deliberately not    //
// ported: int-truncated sums / abs(), and the seed-vs-strip loop that stops  //
// at the seed index.                                                         //
//============================================================================//

#include "GemCluster.h"
#include <algorithm>
#include <cmath>
#include <limits>

using namespace gem;

//=============================================================================
// Strip time-sample helpers (public, see GemCluster.h)
//=============================================================================

float gem::StripMeanTime(const std::vector<float> &ts_adc, float ts_period)
{
    // positive samples only; float accumulators — must stay bit-identical to
    // the historical seed mean time used by the X/Y time cut
    float sum_wt = 0.f, sum_w = 0.f;
    for (size_t i = 0; i < ts_adc.size(); ++i) {
        float w = ts_adc[i];
        if (w > 0.f) {
            sum_wt += w * static_cast<float>(i + 1) * ts_period;
            sum_w  += w;
        }
    }
    return (sum_w > 0.f) ? sum_wt / sum_w
                         : std::numeric_limits<float>::quiet_NaN();
}

float gem::TimeSampleCorrelation(const std::vector<float> &a,
                                 const std::vector<float> &b)
{
    const size_t n = a.size();
    if (n < 2 || b.size() != n)
        return std::numeric_limits<float>::quiet_NaN();

    double mean_a = 0., mean_b = 0.;
    for (size_t i = 0; i < n; ++i) {
        mean_a += a[i];
        mean_b += b[i];
    }
    mean_a /= static_cast<double>(n);
    mean_b /= static_cast<double>(n);

    double s_ab = 0., s_aa = 0., s_bb = 0.;
    for (size_t i = 0; i < n; ++i) {
        const double da = a[i] - mean_a, db = b[i] - mean_b;
        s_ab += da * db;
        s_aa += da * da;
        s_bb += db * db;
    }
    // zero variance (flat samples) or NaN input → undefined
    if (!(s_aa > 0.) || !(s_bb > 0.))
        return std::numeric_limits<float>::quiet_NaN();

    // clamp rounding excursions just outside [-1, 1]
    const double r = s_ab / std::sqrt(s_aa * s_bb);
    return static_cast<float>(std::max(-1., std::min(1., r)));
}

bool gem::IsUnimodalPulse(const std::vector<float> &ts_adc)
{
    if (ts_adc.empty()) return false;

    // first maximum (strict >), as SBS Cuts::__get_max_timebin
    size_t max_bin = 0;
    for (size_t i = 1; i < ts_adc.size(); ++i)
        if (ts_adc[i] > ts_adc[max_bin]) max_bin = i;

    // strictly rising up to the maximum ...
    for (size_t i = 0; i < max_bin; ++i)
        if (!(ts_adc[i] < ts_adc[i + 1])) return false;
    // ... and strictly falling after it (SBS's extra f'' < 0 check at the
    // peak is implied by the two strict comparisons)
    for (size_t i = max_bin; i + 1 < ts_adc.size(); ++i)
        if (!(ts_adc[i] > ts_adc[i + 1])) return false;
    return true;
}

//=============================================================================
// Construction / destruction
//=============================================================================

GemCluster::GemCluster()
{
    // Default characteristic distances for cross-talk identification (mm)
    cfg_.charac_dists = {6.4f, 17.6f, 24.4f, 24.8f, 25.2f, 25.6f,
                         26.0f, 26.4f, 26.8f, 33.6f, 44.8f};
}

GemCluster::~GemCluster() = default;

//=============================================================================
// FormClusters — main entry point
//=============================================================================

void GemCluster::FormClusters(std::vector<StripHit> &hits,
                              std::vector<StripCluster> &clusters) const
{
    clusters.clear();
    if (hits.empty()) return;

    // group consecutive hits → preliminary clusters (with splitting)
    groupHits(hits, clusters);

    // reconstruct cluster position (+ SBS-style quality variables)
    for (auto &cluster : clusters)
        reconstructCluster(cluster);

    // mark cross-talk clusters
    setCrossTalk(clusters);

    // filter out bad clusters
    filterClusters(clusters);
}

//=============================================================================
// groupHits — sort by strip, then cluster consecutive strips
//=============================================================================

void GemCluster::groupHits(std::vector<StripHit> &hits,
                           std::vector<StripCluster> &clusters) const
{
    // sort by strip number
    std::sort(hits.begin(), hits.end(),
              [](const StripHit &a, const StripHit &b) {
                  return a.strip < b.strip;
              });

    // SBS-style strip cuts (IsGoodStrip).  A failing strip closes the open
    // run and is left out of every cluster; it is not erased from `hits`.
    // With no strip cut active the loop is the plain consecutive grouping.
    const bool strip_cuts = cfg_.strip_unimodal ||
                            std::isfinite(cfg_.strip_time_min) ||
                            std::isfinite(cfg_.strip_time_max);

    // cluster consecutive hits
    auto cbeg = hits.begin();
    for (auto it = hits.begin(); it != hits.end(); ++it) {
        if (strip_cuts && !isGoodStrip(*it)) {
            if (cbeg != it) splitCluster(cbeg, it, cfg_.split_thres, clusters);
            cbeg = it + 1;
            continue;
        }
        auto it_n = it + 1;
        if (it_n == hits.end() ||
            it_n->strip - it->strip > cfg_.consecutive_thres)
        {
            splitCluster(cbeg, it_n, cfg_.split_thres, clusters);
            cbeg = it_n;
        }
    }
}

//=============================================================================
// isGoodStrip — SBS-style strip-level cuts (only called when one is active)
//=============================================================================

bool GemCluster::isGoodStrip(const StripHit &hit) const
{
    // pulse shape: strictly rising then strictly falling (empty fails)
    if (cfg_.strip_unimodal && !IsUnimodalPulse(hit.ts_adc))
        return false;

    // mean-time window, inclusive; NaN (no positive sample) fails
    if (std::isfinite(cfg_.strip_time_min) || std::isfinite(cfg_.strip_time_max)) {
        float t = StripMeanTime(hit.ts_adc, cfg_.ts_period);
        if (!(t >= cfg_.strip_time_min && t <= cfg_.strip_time_max))
            return false;
    }
    return true;
}

//=============================================================================
// splitCluster — recursively split at local charge minima (valleys)
//=============================================================================

void GemCluster::splitCluster(std::vector<StripHit>::iterator beg,
                              std::vector<StripHit>::iterator end,
                              float thres,
                              std::vector<StripCluster> &clusters) const
{
    auto size = end - beg;
    if (size <= 0) return;

    // Don't split clusters smaller than 3 strips
    if (size < 3) {
        StripCluster cl;
        cl.hits.assign(beg, end);
        clusters.push_back(std::move(cl));
        return;
    }

    // Find the first local minimum (valley)
    bool descending = false, extremum = false;
    auto minimum = beg;

    for (auto it = beg, it_n = beg + 1; it_n != end; ++it, ++it_n) {
        if (descending) {
            if (it->charge < minimum->charge)
                minimum = it;
            // ascending trend confirms valley
            if (it_n->charge - it->charge > thres) {
                extremum = true;
                break;
            }
        } else {
            // descending trend — potential valley ahead
            if (it->charge - it_n->charge > thres) {
                descending = true;
                minimum = it_n;
            }
        }
    }

    if (extremum) {
        // halve the charge of the overlap strip
        minimum->charge /= 2.f;

        // left sub-cluster
        StripCluster cl;
        cl.hits.assign(beg, minimum);
        clusters.push_back(std::move(cl));

        // recurse on right portion
        splitCluster(minimum, end, thres, clusters);
    } else {
        StripCluster cl;
        cl.hits.assign(beg, end);
        clusters.push_back(std::move(cl));
    }
}

//=============================================================================
// reconstructCluster — charge-weighted position + SBS-style quality variables
//=============================================================================

void GemCluster::reconstructCluster(StripCluster &cluster) const
{
    if (cluster.hits.empty()) return;

    cluster.total_charge = 0.f;
    cluster.peak_charge  = 0.f;
    cluster.max_timebin  = -1;
    float weight_pos = 0.f;

    for (auto &hit : cluster.hits) {
        if (hit.charge > cluster.peak_charge) {
            cluster.peak_charge = hit.charge;
            cluster.max_timebin = hit.max_timebin;
        }
        cluster.total_charge += hit.charge;
        weight_pos += hit.position * hit.charge;
    }

    if (cluster.total_charge > 0.f)
        cluster.position = weight_pos / cluster.total_charge;

    // --- SBS-style quality variables -------------------------------------
    // seed = first strip with the maximum charge (strict >), the rule the
    // X/Y time cut has always used
    const float nan = std::numeric_limits<float>::quiet_NaN();
    const StripHit *seed = &cluster.hits.front();
    for (auto &hit : cluster.hits)
        if (hit.charge > seed->charge) seed = &hit;

    cluster.seed_time     = StripMeanTime(seed->ts_adc, cfg_.ts_period);
    cluster.seed_peak_adc = nan;
    cluster.seed_sum_adc  = nan;
    if (!seed->ts_adc.empty()) {
        float peak = seed->ts_adc.front(), sum = 0.f;
        for (float v : seed->ts_adc) {
            if (v > peak) peak = v;
            sum += v;
        }
        cluster.seed_peak_adc = peak;
        cluster.seed_sum_adc  = sum;
    }

    // seed vs every other strip (all i != seed, unlike SBS which stops at
    // the seed index); non-finite per-strip values are skipped
    cluster.max_strip_dt = nan;
    cluster.min_ts_corr  = nan;
    for (auto &hit : cluster.hits) {
        if (&hit == seed) continue;
        float dt = std::abs(StripMeanTime(hit.ts_adc, cfg_.ts_period) -
                            cluster.seed_time);
        if (std::isfinite(dt) &&
            (std::isnan(cluster.max_strip_dt) || dt > cluster.max_strip_dt))
            cluster.max_strip_dt = dt;
        float r = TimeSampleCorrelation(seed->ts_adc, hit.ts_adc);
        if (std::isfinite(r) &&
            (std::isnan(cluster.min_ts_corr) || r < cluster.min_ts_corr))
            cluster.min_ts_corr = r;
    }
}

//=============================================================================
// setCrossTalk — mark clusters at characteristic cross-talk distances
//=============================================================================

namespace {

// Check if all hits in a cluster are cross-talk strips
inline bool isPureCrossTalk(const StripCluster &cl)
{
    for (auto &hit : cl.hits)
        if (!hit.cross_talk) return false;
    return true;
}

// Check if cluster is at a characteristic CT distance from any later cluster
// and has a peak charge ratio below the threshold.
inline bool atCTDistance(std::vector<StripCluster>::iterator it,
                        std::vector<StripCluster>::iterator end,
                        float width,
                        const std::vector<float> &charac,
                        float peak_ratio_max)
{
    for (auto itn = it + 1; itn != end; ++itn) {
        bool dist_match = false;
        float delta = std::abs(it->position - itn->position);
        for (float dist : charac) {
            if (delta > dist - width && delta < dist + width) {
                dist_match = true;
                break;
            }
        }
        if (!dist_match || itn->peak_charge <= 0.f) continue;

        const float peak_ratio = it->peak_charge / itn->peak_charge;
        if (peak_ratio < peak_ratio_max)
            return true;
    }
    return false;
}

} // anonymous namespace

// // Below is the original implementation of setCrossTalk().

// void GemCluster::setCrossTalk(std::vector<StripCluster> &clusters) const
// {
//     if (cfg_.charac_dists.empty()) return;

//     // sort by peak charge ascending (check weakest clusters first)
//     std::sort(clusters.begin(), clusters.end(),
//               [](const StripCluster &a, const StripCluster &b) {
//                   return a.peak_charge < b.peak_charge;
//               });

//     for (auto it = clusters.begin(); it != clusters.end(); ++it) {
//         if (!isPureCrossTalk(*it)) continue;
//         it->cross_talk = atCTDistance(it, clusters.end(),
//                                      cfg_.cross_talk_width, cfg_.charac_dists,
//                                      cfg_.cross_talk_peak_ratio_max);
//     }
// }

// Below is the modified implementation of setCrossTalk() that looks at the charge ratio between the weak cluster and the strong cluster

void GemCluster::setCrossTalk(std::vector<StripCluster> &clusters) const
{
    if (cfg_.charac_dists.empty()) return;

    // Keep the original sorting: weakest to strongest
    std::sort(clusters.begin(), clusters.end(),
              [](const StripCluster &a, const StripCluster &b) {
                  return a.peak_charge < b.peak_charge;
              });

    for (auto it = clusters.begin(); it != clusters.end(); ++it) {
        // candidate cluster must already be made entirely of cross-talk-like strips
        // and has a peak charge ratio below the threshold.
        if (!isPureCrossTalk(*it)) continue;
        it->cross_talk = atCTDistance(it, clusters.end(),
                                     cfg_.cross_talk_width, cfg_.charac_dists,
                                     cfg_.cross_talk_peak_ratio_max);
    }
}


//=============================================================================
// filterClusters — remove bad clusters
//=============================================================================

void GemCluster::filterClusters(std::vector<StripCluster> &clusters) const
{
    clusters.erase(
        std::remove_if(clusters.begin(), clusters.end(),
            [this](const StripCluster &cl) {
                // bad size
                int sz = static_cast<int>(cl.hits.size());
                if (sz < cfg_.min_cluster_hits || sz > cfg_.max_cluster_hits)
                    return true;
                // cross-talk
                if (cl.cross_talk)
                    return true;
                // SBS-style quality cuts, each only when enabled; NaN
                // (undefined, e.g. single-strip cluster) passes
                if (cfg_.seed_min_peak_adc > 0.f &&
                    cl.seed_peak_adc < cfg_.seed_min_peak_adc)
                    return true;
                if (cfg_.seed_min_sum_adc > 0.f &&
                    cl.seed_sum_adc < cfg_.seed_min_sum_adc)
                    return true;
                if (cfg_.strip_time_agreement >= 0.f &&
                    cl.max_strip_dt > cfg_.strip_time_agreement)
                    return true;
                if (cfg_.strip_ts_corr_min > -1.f &&
                    cl.min_ts_corr < cfg_.strip_ts_corr_min)
                    return true;
                return false;
            }),
        clusters.end());
}

//=============================================================================
// CartesianReconstruct — match X and Y clusters to form 2D hits
//
// Mode 0 (ADC-sorted): sort by peak charge, pair 1:1 by rank
// Mode 1 (Cartesian):  all X×Y combinations with optional cuts:
//   - ADC asymmetry: |Qx_peak - Qy_peak| / (Qx_peak + Qy_peak) <= threshold
//   - Timing:        |mean_time_x_seed - mean_time_y_seed| <= threshold
//                    (StripCluster::seed_time, filled by reconstructCluster)
// Both modes record the X/Y seed times, time difference, signed ADC
// asymmetry and the cluster quality variables on every GEMHit.
//=============================================================================

static GEMHit makeHit(const StripCluster &xc, const StripCluster &yc,
                       int det_id)
{
    GEMHit hit;
    hit.x = xc.position;
    hit.y = yc.position;
    hit.z = 0.f;
    hit.det_id = det_id;
    hit.x_charge = xc.total_charge;
    hit.y_charge = yc.total_charge;
    hit.x_peak   = xc.peak_charge;
    hit.y_peak   = yc.peak_charge;
    hit.x_max_timebin = xc.max_timebin;
    hit.y_max_timebin = yc.max_timebin;
    hit.x_size = static_cast<int>(xc.hits.size());
    hit.y_size = static_cast<int>(yc.hits.size());

    // SBS-style X/Y quality (NaN = undefined)
    hit.x_time    = xc.seed_time;
    hit.y_time    = yc.seed_time;
    hit.time_diff = xc.seed_time - yc.seed_time;
    float sum = xc.peak_charge + yc.peak_charge;
    hit.adc_asym  = (sum > 0.f) ? (xc.peak_charge - yc.peak_charge) / sum
                                : std::numeric_limits<float>::quiet_NaN();
    hit.x_max_strip_dt = xc.max_strip_dt;
    hit.y_max_strip_dt = yc.max_strip_dt;
    hit.x_min_ts_corr  = xc.min_ts_corr;
    hit.y_min_ts_corr  = yc.min_ts_corr;
    return hit;
}

void GemCluster::CartesianReconstruct(
    const std::vector<StripCluster> &x_clusters,
    const std::vector<StripCluster> &y_clusters,
    std::vector<GEMHit> &container,
    int det_id) const
{
    container.clear();

    // Mode 0: ADC-sorted 1:1 matching
    if (cfg_.match_mode == 0) {
        std::vector<StripCluster> xc = x_clusters;
        std::vector<StripCluster> yc = y_clusters;
        std::sort(xc.begin(), xc.end(),
                  [](const StripCluster &a, const StripCluster &b) {
                      return a.peak_charge > b.peak_charge;
                  });
        std::sort(yc.begin(), yc.end(),
                  [](const StripCluster &a, const StripCluster &b) {
                      return a.peak_charge > b.peak_charge;
                  });
        size_t npairs = std::min(xc.size(), yc.size());
        for (size_t i = 0; i < npairs; ++i)
            container.push_back(makeHit(xc[i], yc[i], det_id));
        return;
    }

    // Mode 1: full Cartesian product with cuts
    const float adc_asym_cut = cfg_.match_adc_asymmetry;
    const float time_cut     = cfg_.match_time_diff;

    for (auto &xc : x_clusters) {
        for (auto &yc : y_clusters) {
            // ADC asymmetry cut
            if (adc_asym_cut >= 0.f) {
                float sum = xc.peak_charge + yc.peak_charge;
                if (sum > 0.f) {
                    float asym = std::abs(xc.peak_charge - yc.peak_charge) / sum;
                    if (asym > adc_asym_cut) continue;
                }
            }

            // timing asymmetry cut (seed mean times; skipped when either
            // is undefined, i.e. the seed has no positive sample)
            if (time_cut >= 0.f &&
                std::isfinite(xc.seed_time) && std::isfinite(yc.seed_time) &&
                std::abs(xc.seed_time - yc.seed_time) > time_cut)
                continue;

            container.push_back(makeHit(xc, yc, det_id));
        }
    }
}
