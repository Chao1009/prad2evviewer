#pragma once
//=============================================================================
// MatchingTools.h — detector matching tools for PRad2
//
// adapted from PRadAnalyzer/PRadDetMatch.cpp
//=============================================================================

#include "PhysicsTools.h"
#include "DetectorTransform.h"
#include <array>
#include <vector>
#include <cstdint>
#include <cmath>
#include <limits>
#include <type_traits>

namespace fdec  { struct ClusterHit; }
namespace gem   { struct GEMHit; }
namespace prad2 { struct RunConfig; struct ReconEventData; }

namespace analysis {

// --- matching flag bit positions --------------------------------------------
enum MatchFlag : uint32_t {
    kGEM1Match = 0,
    kGEM2Match = 1,
    kGEM3Match = 2,
    kGEM4Match = 3,
};

struct ProjectHit
{
    float x_proj;
    float y_proj;
    float z_proj;

    ProjectHit(float x, float y, float z) : x_proj(x), y_proj(y), z_proj(z) {};
};

ProjectHit GetProjectionHits(float x, float y, float z, float projection_z);

// Move a lab-frame hit along the line from the target (origin) onto the
// plane z = projection_z.
template <typename Hit>
inline void GetProjection(Hit &hit, float projection_z)
{
    static_assert(std::is_same_v<Hit, HCHit> || std::is_same_v<Hit, GEMHit>,
                  "GetProjection expects an HCHit or a GEMHit");
    const ProjectHit p = GetProjectionHits(hit.x, hit.y, hit.z, projection_z);
    hit.x = p.x_proj;
    hit.y = p.y_proj;
    hit.z = p.z_proj;
}

template <typename Hit>
inline void GetProjection(std::vector<Hit> &hits, float projection_z)
{
    for (auto &hit : hits) GetProjection(hit, projection_z);
}

// Detector-frame HyCal cluster -> lab-frame HCHit.  The cluster is placed at
// its shower-max depth behind the HyCal face (z = 0 on the face when
// at_shower_depth is false) before the transform.
HCHit ClusterToLab(const DetectorTransform &hycal_xform, const fdec::ClusterHit &cluster,
                   bool at_shower_depth = true);

// Detector-plane GEM hit -> lab-frame GEMHit using gem_xforms[hit.det_id];
// a hit with det_id outside 0..3 keeps its plane coordinates.
GEMHit GemHitToLab(const std::array<DetectorTransform, 4> &gem_xforms, const gem::GEMHit &hit);

class MatchHit
{
    public:
        analysis::HCHit hycal_hit;
        std::vector<analysis::GEMHit> gem1_hits;
        std::vector<analysis::GEMHit> gem2_hits;
        std::vector<analysis::GEMHit> gem3_hits;
        std::vector<analysis::GEMHit> gem4_hits;

        MatchHit(const analysis::HCHit &hycal_hit, const std::vector<analysis::GEMHit> &g1, const std::vector<analysis::GEMHit> &g2,
                 const std::vector<analysis::GEMHit> &g3, const std::vector<analysis::GEMHit> &g4)
            : hycal_hit(hycal_hit), gem1_hits(g1), gem2_hits(g2), gem3_hits(g3), gem4_hits(g4) {}

        // chosen GEM hits: [0] downstream pair GEM1/GEM2 (det_id 0/1),
        //                  [1] upstream pair GEM3/GEM4 (det_id 2/3)
        analysis::GEMHit gem[2];
        uint32_t   mflag = 0;     //E.g. 0101 matching flags (see MatchFlag enum)
        uint8_t    hycal_idx = 0; // index into original hycal vector
};

class MatchHit_perChamber
{
    public:
        analysis::HCHit hycal_hit;

        MatchHit_perChamber(const analysis::HCHit &hycal_hit)
            : hycal_hit(hycal_hit){}

        // hits inside the matching window per GEM chamber ([0..3] = GEM1..GEM4), closest first
        std::vector<analysis::GEMHit> gem_hits[4];
        uint32_t   mflag = 0;     // matching flag
        uint8_t    hycal_idx = 0; // index into original hycal vector
};

class MatchingTools
{
public:
    explicit MatchingTools(int postMatchMethod = 1);

    std::vector<MatchHit> Match(const std::vector<analysis::HCHit> &hycalHits,
                            const std::vector<analysis::GEMHit> &gem1_hits,
                            const std::vector<analysis::GEMHit> &gem2_hits,
                            const std::vector<analysis::GEMHit> &gem3_hits,
                            const std::vector<analysis::GEMHit> &gem4_hits) const;

    std::vector<MatchHit_perChamber> MatchPerChamber(const std::vector<analysis::HCHit> &hycalHits,
                            const std::vector<analysis::GEMHit> &gem1_hits,
                            const std::vector<analysis::GEMHit> &gem2_hits,
                            const std::vector<analysis::GEMHit> &gem3_hits,
                            const std::vector<analysis::GEMHit> &gem4_hits) const;

    void SetMatchRange(float range)    { matchRange_ = range; }
    void SetSquareSelection(bool sq)   { squareSel_ = sq; }
    void SetEnergyDependent(bool ed)   { energyDependent_ = ed; }
    void SetMatchSigma(float sigma)    { matchSigma_ = sigma; }
    // all four from the run's matching_* settings
    void Configure(const prad2::RunConfig &cfg);

private:
    float matchRange_ = 15.f;   // mm, spatial matching window
    bool  squareSel_  = true;   // true = square window, false = circular

    float matchSigma_ = 1.f;   // mm, used for energy-dependent matching
    bool  energyDependent_ = false;
    int   postMatchMethod_ = 1;

    float ProjectionDistance(const analysis::HCHit &h, const analysis::GEMHit &g) const;
    float ProjectionDistance(const analysis::GEMHit &g1, const analysis::GEMHit &g2, float ref_z) const;
    bool  PreMatch(const analysis::HCHit &hycal, const analysis::GEMHit &gem) const;
    void  PostMatch(MatchHit &h) const;
};

// Write Match / MatchPerChamber results into the recon-tree matching fields
// (match lists, matchFlag, matchNum, mHit_*); clusters are indexed by hycal_idx.
void FillReconMatches(prad2::ReconEventData &ev,
                      const std::vector<MatchHit> &matched,
                      const std::vector<MatchHit_perChamber> &matched_per_chamber);

// Match the event's lab-frame clusters (cl_*) to its GEM hits (gem_*, det_id
// 0..3 only) and fill the results with FillReconMatches.
void MatchReconEvent(prad2::ReconEventData &ev, const MatchingTools &matching);

} // namespace analysis
