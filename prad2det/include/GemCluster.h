#pragma once
//=============================================================================
// GemCluster.h — GEM strip clustering and 2D hit reconstruction
//
// Algorithms ported from mpd_gem_view_ssp GEMCluster:
//   1. Group consecutive strip hits (optional SBS-style strip cuts)
//   2. Split multi-peak clusters at local minima
//   3. Charge-weighted position reconstruction + SBS-style quality variables
//      (seed mean time / peak / sum, seed-vs-strip time spread and
//      time-sample correlation)
//   4. Cross-talk identification by characteristic distance
//   5. Cluster filter (size, cross-talk, optional SBS-style quality cuts)
//   6. Cartesian X-Y cluster matching → 2D hits (X/Y time + ADC asymmetry
//      recorded on every GEMHit)
// The SBS-style cuts are all off by default (see ClusterConfig).
//=============================================================================

#include "GemSystem.h"   // ClusterConfig + StripHit/StripCluster/GEMHit
#include <vector>

namespace gem
{

// --- strip time-sample helpers (public so offline code can reuse them) ------

// ADC-weighted mean time of one strip (ns):
//   t = sum_{i: a_i>0} a_i*(i+1)*ts_period / sum_{i: a_i>0} a_i
// Only positive samples contribute (SBS includes negative samples).  Sample 0
// maps to 1*ts_period (SBS convention).  NaN if there is no positive sample.
// This is the seed mean time used by the X/Y time cut.
float StripMeanTime(const std::vector<float> &ts_adc, float ts_period = 25.f);

// Pearson correlation coefficient between two time-sample vectors (SBS "time
// sample correlation coefficient").  NaN if the sizes differ, size < 2, or
// either variance is zero.
float TimeSampleCorrelation(const std::vector<float> &a,
                            const std::vector<float> &b);

// SBS "concave shape" strip cut (Cuts::is_concave_shape): samples strictly
// rising up to the first maximum and strictly falling after it.  A peak in
// the first or last sample passes.  Empty => false.
bool IsUnimodalPulse(const std::vector<float> &ts_adc);

class GemCluster
{
public:
    GemCluster();
    ~GemCluster();

    void SetConfig(const ClusterConfig &cfg) { cfg_ = cfg; }
    const ClusterConfig &GetConfig() const   { return cfg_; }

    // Form 1D strip clusters from a list of strip hits.
    // hits will be sorted by strip number.
    void FormClusters(std::vector<StripHit> &hits,
                      std::vector<StripCluster> &clusters) const;

    // Match X and Y clusters to form 2D hits via Cartesian product.
    void CartesianReconstruct(const std::vector<StripCluster> &x_clusters,
                              const std::vector<StripCluster> &y_clusters,
                              std::vector<GEMHit> &hits,
                              int det_id) const;

private:
    // Group consecutive hits into preliminary clusters
    void groupHits(std::vector<StripHit> &hits,
                   std::vector<StripCluster> &clusters) const;

    // SBS-style strip-level cuts (time window, unimodal shape)
    bool isGoodStrip(const StripHit &hit) const;

    // Recursively split clusters at local charge minima
    void splitCluster(std::vector<StripHit>::iterator beg,
                      std::vector<StripHit>::iterator end,
                      float thres,
                      std::vector<StripCluster> &clusters) const;

    // Compute charge-weighted position and quality variables for a cluster
    void reconstructCluster(StripCluster &cluster) const;

    // Mark cross-talk clusters by characteristic distance
    void setCrossTalk(std::vector<StripCluster> &clusters) const;

    // Filter out bad clusters
    void filterClusters(std::vector<StripCluster> &clusters) const;

    ClusterConfig cfg_;
};

} // namespace gem
