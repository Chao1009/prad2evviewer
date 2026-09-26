// MatchingTools.cpp — tools for matching HyCal clusters to GEM hits
//=============================================================================
// Adapted from PRadAnalyzer/PRadDetMatch.cpp.
// Usage:
//   MatchingTools matcher;
//   auto matches = matcher.Match(hycalHits, gem1Hits, gem2Hits, gem3Hits, gem4Hits);
//   for (auto &m : matches) { 
//       m.hycal_hit is the HyCal cluster
//       m.gem[0] / m.gem[1] are the chosen GEM1/GEM2 and GEM3/GEM4 hits
//       m.gem1_hits, m.gem2_hits, m.gem3_hits, m.gem4_hits are the candidates in each plane (sorted by distance)
//       m.mflag has the MatchFlag bits of the planes of m.gem[0] / m.gem[1]
//       m.hycal_idx is the index of the cluster in the original vector
//   }
//=============================================================================


#include "MatchingTools.h"
#include "ConfigSetup.h"
#include "EventData.h"
#include "GemSystem.h"
#include "HyCalCluster.h"
#include <algorithm>
#include <set>

namespace analysis {

MatchingTools::MatchingTools(int postMatchMethod)
{
    postMatchMethod_ = postMatchMethod;
}

void MatchingTools::Configure(const prad2::RunConfig &cfg)
{
    SetMatchRange(cfg.matching_radius);
    SetSquareSelection(cfg.matching_use_square);
    SetEnergyDependent(cfg.matching_energy_dependent);
    SetMatchSigma(cfg.matching_sigma);
}

// --- Projection and lab-frame conversion ------------------------------------

ProjectHit GetProjectionHits(float x, float y, float z,
                                            float projection_z)
{   
    // simple linear projection from (x,y,z) to (x_proj, y_proj, projection_z)
    // in target and beam center coordinates
    float scale = projection_z / z;
    return ProjectHit(x * scale, y * scale, projection_z);
}

HCHit ClusterToLab(const DetectorTransform &hycal_xform, const fdec::ClusterHit &cluster,
                   bool at_shower_depth)
{
    HCHit hit = {cluster.x, cluster.y,
                 at_shower_depth ? fdec::shower_depth(cluster.center_id, cluster.energy) : 0.f,
                 cluster.energy, static_cast<uint16_t>(cluster.center_id), cluster.flag};
    ApplyToLab(hycal_xform, hit);
    return hit;
}

GEMHit GemHitToLab(const std::array<DetectorTransform, 4> &gem_xforms, const gem::GEMHit &hit)
{
    GEMHit lab = {hit.x, hit.y, 0.f, static_cast<uint8_t>(hit.det_id)};
    if (hit.det_id >= 0 && hit.det_id < 4)
        ApplyToLab(gem_xforms[hit.det_id], lab);
    return lab;
}

// Distance between HyCal cluster and GEM hit after projecting GEM to HyCal z
float MatchingTools::ProjectionDistance(const analysis::HCHit &h,
                                       const analysis::GEMHit &g) const
{
    ProjectHit proj = GetProjectionHits(g.x, g.y, g.z, h.z);
    float dx = h.x - proj.x_proj;
    float dy = h.y - proj.y_proj;
    return std::sqrt(dx * dx + dy * dy);
}

// Distance between two GEM hits projected to a common reference z
float MatchingTools::ProjectionDistance(const analysis::GEMHit &g1,
                                       const analysis::GEMHit &g2,
                                       float ref_z) const
{
    ProjectHit p1 = GetProjectionHits(g1.x, g1.y, g1.z, ref_z);
    ProjectHit p2 = GetProjectionHits(g2.x, g2.y, g2.z, ref_z);
    float dx = p1.x_proj - p2.x_proj;
    float dy = p1.y_proj - p2.y_proj;
    return std::sqrt(dx * dx + dy * dy);
}

// Pre-match: check if a GEM hit falls within the matching window of a cluster
bool MatchingTools::PreMatch(const analysis::HCHit &hycal,
                             const analysis::GEMHit &gem) const
{
    ProjectHit proj = GetProjectionHits(gem.x, gem.y, gem.z, hycal.z);
    float dx = std::fabs(hycal.x - proj.x_proj);
    float dy = std::fabs(hycal.y - proj.y_proj);

    float match_range = matchRange_;
    if (energyDependent_ && hycal.energy > 0.f) {
        match_range = matchSigma_ * (2.7924f / std::sqrt(hycal.energy / 1000.f)
                                    -0.2709f / (hycal.energy / 1000.f) - 0.0295f );
    }

    if (squareSel_) {
        return (dx <= match_range) && (dy <= match_range);
    } else {
        return (dx * dx + dy * dy) <= match_range * match_range;
    }
}

// Post-match: sort candidates per plane, pick the GEM1/GEM2 and GEM3/GEM4 hits,
// set flags
void MatchingTools::PostMatch(MatchHit &h) const
{
    if( (h.gem1_hits.empty() && h.gem2_hits.empty() ) ||
        (h.gem3_hits.empty() && h.gem4_hits.empty() ) )
        return; // require at least one match in both upstream and downstream pairs

    // sort each plane's candidates by projection distance (closest first)
    auto by_dist = [this, &h](const analysis::GEMHit &a, const analysis::GEMHit &b) {
        return ProjectionDistance(h.hycal_hit, a) < ProjectionDistance(h.hycal_hit, b);
    };
    for (auto *plane : {&h.gem1_hits, &h.gem2_hits, &h.gem3_hits, &h.gem4_hits})
        std::sort(plane->begin(), plane->end(), by_dist);

    if (postMatchMethod_ == 1) {
        // in each pair, the plane front() closest to the HyCal cluster
        auto closest = [&](const std::vector<analysis::GEMHit> &a,
                           const std::vector<analysis::GEMHit> &b) {
            float best_d = 1e9f;
            analysis::GEMHit best{};
            for (const auto *plane : {&a, &b}) {
                if (plane->empty()) continue;
                float d = ProjectionDistance(h.hycal_hit, plane->front());
                if (d < best_d) {
                    best_d = d;
                    best = plane->front();
                }
            }
            return best;
        };
        h.gem[0] = closest(h.gem1_hits, h.gem2_hits);
        h.gem[1] = closest(h.gem3_hits, h.gem4_hits);
    } else {
        //start from the closest candidate on upstream GEM planes (gem4 -> gem3)
        //for each hits on upstream GEM planes, loop over the sorted candidates on downstream GEM planes to find the best match
        //get a straightline from target to the upstream GEM hit, and project it to the downstream GEM planes to find the closest match
        //save the deltaR between the projected upstream hit and the downstream hit, finnaly find the best match by deltaR
        float best_deltaR = std::numeric_limits<float>::max();
        analysis::GEMHit best_gem_down{}, best_gem_up{};

        auto check_pair = [&](const analysis::GEMHit &up, const analysis::GEMHit &down) {
            float deltaR = ProjectionDistance(up, down, down.z);
            if (deltaR < best_deltaR) {
                best_deltaR = deltaR;
                best_gem_up = up;
                best_gem_down = down;
            }
        };

        // scan upstream candidates in priority order: GEM4 first, then GEM3
        for (const auto *ups : {&h.gem4_hits, &h.gem3_hits})
            for (const auto &up : *ups) {
                for (const auto &down : h.gem1_hits) check_pair(up, down);
                for (const auto &down : h.gem2_hits) check_pair(up, down);
            }

        if (best_deltaR == std::numeric_limits<float>::max()) return;

        h.gem[0] = best_gem_down;
        h.gem[1] = best_gem_up;
    }

    // set match flag for downstream pair
    if (h.gem[0].det_id == 0) fdec::set_bit(h.mflag, kGEM1Match);
    else if (h.gem[0].det_id == 1) fdec::set_bit(h.mflag, kGEM2Match);

    // set match flag for upstream pair
    if (h.gem[1].det_id == 2) fdec::set_bit(h.mflag, kGEM3Match);
    else if (h.gem[1].det_id == 3) fdec::set_bit(h.mflag, kGEM4Match);
}

// --- Main matching, adapted from PRadDetMatch::Match for 4 GEM planes -------

// comparator so GEMHit can be stored in std::set (identity by position)
struct GEMHitCmp {
    bool operator()(const analysis::GEMHit &a, const analysis::GEMHit &b) const
    {
        if (a.z != b.z) return a.z < b.z;
        if (a.x != b.x) return a.x < b.x;
        return a.y < b.y;
    }
};

std::vector<MatchHit> MatchingTools::Match(
    const std::vector<analysis::HCHit> &hycalHits,
    const std::vector<analysis::GEMHit> &gem1,
    const std::vector<analysis::GEMHit> &gem2,
    const std::vector<analysis::GEMHit> &gem3,
    const std::vector<analysis::GEMHit> &gem4) const
{
    const std::vector<analysis::GEMHit> *planes[4] = {&gem1, &gem2, &gem3, &gem4};
    std::vector<MatchHit> result;

    // keep track of GEM hits already claimed (higher-E cluster gets priority);
    // hycalHits are expected sorted by energy already (ReconstructHits order)
    std::set<analysis::GEMHit, GEMHitCmp> used[4];

    for (size_t i = 0; i < hycalHits.size(); ++i) {
        const auto &hit = hycalHits[i];

        // collect candidates in each GEM plane
        std::vector<analysis::GEMHit> cand[4];
        for (int d = 0; d < 4; ++d)
            for (const auto &g : *planes[d])
                if (PreMatch(hit, g) && used[d].find(g) == used[d].end())
                    cand[d].push_back(g);

        // skip if no candidates in any plane in one of the pairs (upstream or downstream)
        if ((cand[0].empty() && cand[1].empty()) || (cand[2].empty() && cand[3].empty()))
            continue;

        result.emplace_back(hit, cand[0], cand[1], cand[2], cand[3]);
        MatchHit &mhit = result.back();
        mhit.hycal_idx = static_cast<uint8_t>(i);

        // resolve best match and set flags
        PostMatch(mhit);

        // mark the closest candidate (front) of each flagged plane as used
        const std::vector<analysis::GEMHit> *sorted[4] =
            {&mhit.gem1_hits, &mhit.gem2_hits, &mhit.gem3_hits, &mhit.gem4_hits};
        for (int d = 0; d < 4; ++d)
            if (fdec::test_bit(mhit.mflag, d)) used[d].insert(sorted[d]->front());
    }

    return result;
}

// MatchPerChamber — for each HyCal hit, independently collect the matching
// GEM hits in each of the 4 chambers. No "used" exclusion across clusters.
// gem_hits[d] = all chamber-d hits inside the matching window, closest first.
// mflag bit d is set if chamber d has a match (bit 0 = GEM1, ..., bit 3 = GEM4).
std::vector<MatchHit_perChamber> MatchingTools::MatchPerChamber(
    const std::vector<analysis::HCHit> &hycalHits,
    const std::vector<analysis::GEMHit> &gem1,
    const std::vector<analysis::GEMHit> &gem2,
    const std::vector<analysis::GEMHit> &gem3,
    const std::vector<analysis::GEMHit> &gem4) const
{
    const std::vector<analysis::GEMHit> *planes[4] = {&gem1, &gem2, &gem3, &gem4};

    std::vector<MatchHit_perChamber> result;
    result.reserve(hycalHits.size());

    for (size_t i = 0; i < hycalHits.size(); ++i) {
        const auto &hit = hycalHits[i];
        MatchHit_perChamber mhit(hit);
        mhit.hycal_idx = static_cast<uint8_t>(i);

        auto by_dist = [this, &hit](const analysis::GEMHit &a, const analysis::GEMHit &b) {
            return ProjectionDistance(hit, a) < ProjectionDistance(hit, b);
        };
        for (int d = 0; d < 4; ++d) {
            for (const auto &g : *planes[d])
                if (PreMatch(hit, g)) mhit.gem_hits[d].push_back(g);
            std::sort(mhit.gem_hits[d].begin(), mhit.gem_hits[d].end(), by_dist);
            // flag bit d = GEM(d+1) match (MatchFlag)
            if (!mhit.gem_hits[d].empty()) fdec::set_bit(mhit.mflag, d);
        }

        result.push_back(mhit);
    }

    return result;
}

// --- Recon-tree fill --------------------------------------------------------

void FillReconMatches(prad2::ReconEventData &ev,
                      const std::vector<MatchHit> &matched,
                      const std::vector<MatchHit_perChamber> &matched_per_chamber)
{
    ev.clear_match_lists();
    for (const auto &m : matched_per_chamber) {
        const int cl_idx = m.hycal_idx;
        if (cl_idx < 0 || cl_idx >= ev.n_clusters) continue;
        for (int d = 0; d < 4; ++d)
            for (const auto &g : m.gem_hits[d])
                ev.add_match(cl_idx, d, g.x, g.y, g.z);
        ev.matchFlag[cl_idx] = m.mflag;
    }

    ev.matchNum = std::min(static_cast<int>(matched.size()), prad2::kMaxClusters);
    for (int i = 0; i < ev.matchNum; ++i) {
        const auto &m = matched[i];
        ev.mHit_E[i] = m.hycal_hit.energy;
        ev.mHit_x[i] = m.hycal_hit.x;
        ev.mHit_y[i] = m.hycal_hit.y;
        ev.mHit_z[i] = m.hycal_hit.z;
        for (int k = 0; k < 2; ++k) {
            ev.mHit_gx[i][k]  = m.gem[k].x;
            ev.mHit_gy[i][k]  = m.gem[k].y;
            ev.mHit_gz[i][k]  = m.gem[k].z;
            ev.mHit_gid[i][k] = m.gem[k].det_id;
        }
        ev.mHit_cl_index[i] = m.hycal_idx;
    }
}

void MatchReconEvent(prad2::ReconEventData &ev, const MatchingTools &matching)
{
    std::vector<HCHit> hc_hits;
    std::vector<GEMHit> gem_hits[4];
    hc_hits.reserve(ev.n_clusters);
    for (int i = 0; i < ev.n_clusters; ++i)
        hc_hits.push_back({ev.cl_x[i], ev.cl_y[i], ev.cl_z[i],
                           ev.cl_energy[i], ev.cl_center[i], ev.cl_flag[i]});
    for (int i = 0; i < ev.n_gem_hits; ++i) {
        const int d = ev.det_id[i];
        if (d < 0 || d >= 4) continue;
        gem_hits[d].push_back({ev.gem_x[i], ev.gem_y[i], ev.gem_z[i], ev.det_id[i]});
    }

    FillReconMatches(ev,
        matching.Match(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]),
        matching.MatchPerChamber(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]));
}

} // namespace analysis