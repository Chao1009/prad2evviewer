//============================================================================
// gem_hycal_matching.C — full HyCal+GEM reconstruction with straight-line
// matching, dumping per-match cluster + strip-level info into a ROOT tree.
//
// Pipeline per physics event:
//
//   EvChannel.Read() → DecodeEvent() → FADC + SSP buffers
//                   → HyCalSystem(WaveAnalyzer) → HyCalCluster.FormClusters()
//                   → GemSystem.ProcessEvent() (pedestal + CM + ZS)
//                   → GemSystem.Reconstruct(GemCluster) (1D + 2D matching)
//                   → coord transform to lab (via runinfo geometry)
//                   → straight-line match HyCal cluster → each GEM plane
//                   → for each match, look up the X & Y constituent clusters
//                     by (position, total_charge) on the plane-cluster lists
//                   → TTree.Fill()
//
// Matching geometry: a line from the target through the HyCal cluster
// centroid (lab frame, with shower-depth applied to z) is intersected with
// each GEM plane.  A GEM 2D hit is a match if the 2D residual at the GEM
// plane is within N·σ_total, where
//
//   σ_hc_face = HyCalSystem::PositionResolution(E)   [mm at HyCal face]
//             = sqrt((A/sqrt(E_GeV))^2 + (B/E_GeV)^2 + C^2)
//   σ_hc@gem  = σ_hc_face · (z_gem / z_hc)           [scaled by line geometry]
//   σ_gem     = gem_pos_res[det_id]                  [mm, per-detector]
//   σ_total   = sqrt(σ_hc@gem² + σ_gem²)
//
// (A, B, C) and gem_pos_res are read from
// reconstruction_config.json:matching via the PipelineBuilder.
//
// Defaults: N = 3 (3-sigma).  The actual residual is stored in the tree so
// downstream cuts can be tightened or loosened without re-running.
//
// Best-match rule: HyCal cluster is the baseline.  For each (HC cluster,
// GEM detector) pair we keep at most ONE matched GEM hit — the candidate
// with the smallest 2D residual that's still ≤ N·σ_total.  A GEM hit
// CAN be the best match for multiple HC clusters (no exclusivity on the
// GEM side — that would require Hungarian assignment).  The Python
// counterpart in analysis/pyscripts/gem_hycal_matching.py applies the
// same rule.
//
// Trigger filter: only events with `trigger_bits == 0x100` (production
// physics trigger) are reconstructed and written.  Everything else
// (LMS / Alpha / cosmic / etc.) is skipped.  The summary lines report
// raw physics count vs. kept count.
//
// Multi-file mode is selected by the input path:
//   * `/data/.../prad_023881.evio.*`  → glob: enumerate every sibling
//     `prad_023881.evio.<digits>`, process them in suffix order, and
//     warn (to stderr) about any gap in the suffix sequence (including
//     missing from .00000).
//   * `/data/prad_023881/`            → directory: same enumeration,
//     run number sniffed from the directory name.
//   * `/data/.../prad_023881.evio.00000` → single specific split file.
// Use the glob form for a full-run replay; use the explicit suffix when
// you want to debug or re-process a single segment.
//
// Pedestals, common-mode files, and HyCal calibration are auto-discovered
// from database/reconstruction_config.json -> runinfo (matches the live monitor).  Pass
// "" (empty string) for any of the file args to use the discovered
// defaults; pass an explicit path to override.
//
// Tree (one entry per physics event):
//   event_num, trigger_bits                         scalar
//
//   ncl                                              scalar
//   hc_x, hc_y, hc_z, hc_energy, hc_center,
//   hc_nblocks, hc_flag                              vector<>(ncl)
//   hc_sigma                                         vector<>(ncl) [mm at HC face]
//
//   nmatch                                           scalar
//   m_hc_idx                                         vector<>(nmatch)  HC cluster index
//   m_det                                            vector<>(nmatch)  GEM detector 0..3
//   m_gem_x, m_gem_y, m_gem_z                        vector<>(nmatch)  lab frame
//   m_gem_x_charge, m_gem_y_charge                   vector<>(nmatch)
//   m_gem_x_size, m_gem_y_size                       vector<>(nmatch)  strip count
//   m_proj_x, m_proj_y                               vector<>(nmatch)  HC line @ GEM
//   m_residual, m_sigma_total                        vector<>(nmatch)  matching geometry
//
//   m_xcl_position, m_xcl_total, m_xcl_peak,
//   m_xcl_max_tb                                     vector<>(nmatch)
//   m_xcl_first, m_xcl_nstrips                       vector<>(nmatch)  slice into strip arrays
//   m_ycl_position, m_ycl_total, m_ycl_peak,
//   m_ycl_max_tb                                     vector<>(nmatch)
//   m_ycl_first, m_ycl_nstrips                       vector<>(nmatch)
//
//   nstrips                                          scalar
//   s_match_idx, s_plane (0=X, 1=Y)                  vector<>(nstrips)
//   s_strip, s_position, s_charge, s_max_tb,
//   s_cross_talk                                     vector<>(nstrips)
//   s_ts0 .. s_ts5                                   vector<>(nstrips)  full 6-sample waveform
//
// Usage
// -----
//   cd build
//   root -l ../analysis/scripts/rootlogon.C
//
//   # full run (glob — warns about any missing split):
//   .x ../analysis/scripts/gem_hycal_matching.C+( \
//       "/data/stage6/prad_023867/prad_023867.evio.*", \
//       "match_023867.root")
//
//   # single split (debugging):
//   .x ../analysis/scripts/gem_hycal_matching.C+( \
//       "/data/stage6/prad_023867/prad_023867.evio.00000", \
//       "match_023867_seg0.root")
//============================================================================

#include "DaqConfig.h"
#include "EvioFiles.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "WaveAnalyzer.h"

#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "GemSystem.h"
#include "GemCluster.h"
#include "RunInfoConfig.h"

#include "PhysicsTools.h"
#include "ConfigSetup.h"      // analysis::gRunConfig
#include "MatchingTools.h"    // analysis::ClusterToLab, GemHitToLab
#include "PipelineBuilder.h"  // prad2::Pipeline
#include "script_helpers.h"   // build_script_pipeline, for_each_physics_event

#include <TFile.h>
#include <TTree.h>
#include <TString.h>           // Printf() — line-flushed message output

#include <chrono>
#include <cmath>
#include <cstdint>
#include <string>
#include <vector>

namespace {

// Find the StripCluster in `clusters` whose (position, total_charge) match a
// 2D GEMHit's x or y component.  Returns nullptr if none — should be rare,
// since CartesianReconstruct copies these values verbatim.
static const gem::StripCluster *
find_constituent(const std::vector<gem::StripCluster> &clusters,
                 float position, float total_charge)
{
    constexpr float kPosTol    = 1e-3f;   // mm
    constexpr float kChargeTol = 1e-2f;   // ADC counts
    for (const auto &c : clusters) {
        if (std::abs(c.position - position) <= kPosTol &&
            std::abs(c.total_charge - total_charge) <= kChargeTol)
            return &c;
    }
    return nullptr;
}

// ----- event vars (TTree-bound) ---------------------------------------------

struct EventVars {
    int      event_num = 0;
    uint32_t trigger_bits = 0;

    int                 ncl = 0;
    std::vector<float>  hc_x, hc_y, hc_z;
    std::vector<float>  hc_energy;
    std::vector<short>  hc_center;
    std::vector<int>    hc_nblocks;
    std::vector<unsigned int> hc_flag;
    std::vector<float>  hc_sigma;     // σ at HyCal face

    int                 nmatch = 0;
    std::vector<int>    m_hc_idx, m_det;
    std::vector<float>  m_gem_x, m_gem_y, m_gem_z;
    std::vector<float>  m_gem_x_charge, m_gem_y_charge;
    std::vector<int>    m_gem_x_size,   m_gem_y_size;
    std::vector<float>  m_proj_x, m_proj_y;
    std::vector<float>  m_residual, m_sigma_total;
    std::vector<float>  m_xcl_position, m_xcl_total, m_xcl_peak;
    std::vector<short>  m_xcl_max_tb;
    std::vector<int>    m_xcl_first,   m_xcl_nstrips;
    std::vector<float>  m_ycl_position, m_ycl_total, m_ycl_peak;
    std::vector<short>  m_ycl_max_tb;
    std::vector<int>    m_ycl_first,   m_ycl_nstrips;

    int                 nstrips = 0;
    std::vector<int>    s_match_idx;
    std::vector<int>    s_plane;        // 0=X, 1=Y
    std::vector<int>    s_strip;
    std::vector<float>  s_position;
    std::vector<float>  s_charge;
    std::vector<short>  s_max_tb;
    std::vector<bool>   s_cross_talk;
    std::vector<float>  s_ts[6];

    void clear()
    {
        ncl = 0;
        hc_x.clear(); hc_y.clear(); hc_z.clear();
        hc_energy.clear(); hc_center.clear();
        hc_nblocks.clear(); hc_flag.clear(); hc_sigma.clear();

        nmatch = 0;
        m_hc_idx.clear(); m_det.clear();
        m_gem_x.clear(); m_gem_y.clear(); m_gem_z.clear();
        m_gem_x_charge.clear(); m_gem_y_charge.clear();
        m_gem_x_size.clear();   m_gem_y_size.clear();
        m_proj_x.clear(); m_proj_y.clear();
        m_residual.clear(); m_sigma_total.clear();
        m_xcl_position.clear(); m_xcl_total.clear(); m_xcl_peak.clear();
        m_xcl_max_tb.clear();
        m_xcl_first.clear(); m_xcl_nstrips.clear();
        m_ycl_position.clear(); m_ycl_total.clear(); m_ycl_peak.clear();
        m_ycl_max_tb.clear();
        m_ycl_first.clear(); m_ycl_nstrips.clear();

        nstrips = 0;
        s_match_idx.clear(); s_plane.clear();
        s_strip.clear(); s_position.clear();
        s_charge.clear(); s_max_tb.clear(); s_cross_talk.clear();
        for (int t = 0; t < 6; ++t) s_ts[t].clear();
    }
};

static void bind_branches(TTree *tree, EventVars &ev)
{
    tree->Branch("event_num",    &ev.event_num,    "event_num/I");
    tree->Branch("trigger_bits", &ev.trigger_bits, "trigger_bits/i");

    tree->Branch("ncl",        &ev.ncl, "ncl/I");
    tree->Branch("hc_x",       &ev.hc_x);
    tree->Branch("hc_y",       &ev.hc_y);
    tree->Branch("hc_z",       &ev.hc_z);
    tree->Branch("hc_energy",  &ev.hc_energy);
    tree->Branch("hc_center",  &ev.hc_center);
    tree->Branch("hc_nblocks", &ev.hc_nblocks);
    tree->Branch("hc_flag",    &ev.hc_flag);
    tree->Branch("hc_sigma",   &ev.hc_sigma);

    tree->Branch("nmatch",         &ev.nmatch, "nmatch/I");
    tree->Branch("m_hc_idx",       &ev.m_hc_idx);
    tree->Branch("m_det",          &ev.m_det);
    tree->Branch("m_gem_x",        &ev.m_gem_x);
    tree->Branch("m_gem_y",        &ev.m_gem_y);
    tree->Branch("m_gem_z",        &ev.m_gem_z);
    tree->Branch("m_gem_x_charge", &ev.m_gem_x_charge);
    tree->Branch("m_gem_y_charge", &ev.m_gem_y_charge);
    tree->Branch("m_gem_x_size",   &ev.m_gem_x_size);
    tree->Branch("m_gem_y_size",   &ev.m_gem_y_size);
    tree->Branch("m_proj_x",       &ev.m_proj_x);
    tree->Branch("m_proj_y",       &ev.m_proj_y);
    tree->Branch("m_residual",     &ev.m_residual);
    tree->Branch("m_sigma_total",  &ev.m_sigma_total);
    tree->Branch("m_xcl_position", &ev.m_xcl_position);
    tree->Branch("m_xcl_total",    &ev.m_xcl_total);
    tree->Branch("m_xcl_peak",     &ev.m_xcl_peak);
    tree->Branch("m_xcl_max_tb",   &ev.m_xcl_max_tb);
    tree->Branch("m_xcl_first",    &ev.m_xcl_first);
    tree->Branch("m_xcl_nstrips",  &ev.m_xcl_nstrips);
    tree->Branch("m_ycl_position", &ev.m_ycl_position);
    tree->Branch("m_ycl_total",    &ev.m_ycl_total);
    tree->Branch("m_ycl_peak",     &ev.m_ycl_peak);
    tree->Branch("m_ycl_max_tb",   &ev.m_ycl_max_tb);
    tree->Branch("m_ycl_first",    &ev.m_ycl_first);
    tree->Branch("m_ycl_nstrips",  &ev.m_ycl_nstrips);

    tree->Branch("nstrips",      &ev.nstrips, "nstrips/I");
    tree->Branch("s_match_idx",  &ev.s_match_idx);
    tree->Branch("s_plane",      &ev.s_plane);
    tree->Branch("s_strip",      &ev.s_strip);
    tree->Branch("s_position",   &ev.s_position);
    tree->Branch("s_charge",     &ev.s_charge);
    tree->Branch("s_max_tb",     &ev.s_max_tb);
    tree->Branch("s_cross_talk", &ev.s_cross_talk);
    for (int t = 0; t < 6; ++t)
        tree->Branch(TString::Format("s_ts%d", t), &ev.s_ts[t]);
}

// Append every strip in `cl` to the per-event arrays under match index `mi`
// and plane `plane` (0=X, 1=Y).  Returns (first_index, nstrips) so the
// match-row can record its slice.
static std::pair<int,int>
append_strips(EventVars &ev, int mi, int plane,
              const gem::StripCluster *cl)
{
    if (!cl) return {0, 0};
    const int first = ev.nstrips;
    for (const auto &h : cl->hits) {
        ev.s_match_idx.push_back(mi);
        ev.s_plane.push_back(plane);
        ev.s_strip.push_back(h.strip);
        ev.s_position.push_back(h.position);
        ev.s_charge.push_back(h.charge);
        ev.s_max_tb.push_back(h.max_timebin);
        ev.s_cross_talk.push_back(h.cross_talk);
        for (int t = 0; t < 6; ++t) {
            float v = (t < (int)h.ts_adc.size()) ? h.ts_adc[t] : 0.f;
            ev.s_ts[t].push_back(v);
        }
        ++ev.nstrips;
    }
    return {first, static_cast<int>(cl->hits.size())};
}

} // anonymous namespace

// Full 11-arg version, no defaults (defined below).
int gem_hycal_matching(const char *evio_path,
                         const char *out_path,
                         const char *gem_ped_file,
                         const char *gem_cm_file,
                         const char *hc_calib_file,
                         long        max_events,
                         int         run_num,
                         float       match_nsigma,
                         const char *daq_config,
                         const char *gem_map_file,
                         const char *hc_map_file);

// Convenience overloads forwarding to the full version with "" path
// defaults (auto-discovery via runinfo).  Overloads instead of default
// arguments: cling mis-marshals many mixed-type (`const char*` with
// long/int/float) defaults and SEGVs at the call site.
int gem_hycal_matching(const char *evio_path, const char *out_path)
{
    return gem_hycal_matching(evio_path, out_path,
                                "", "", "", 0L, -1, 3.0f, "", "", "");
}
int gem_hycal_matching(const char *evio_path, const char *out_path,
                         long max_events)
{
    return gem_hycal_matching(evio_path, out_path,
                                "", "", "", max_events, -1, 3.0f, "", "", "");
}
int gem_hycal_matching(const char *evio_path, const char *out_path,
                         long max_events, int run_num)
{
    return gem_hycal_matching(evio_path, out_path,
                                "", "", "", max_events, run_num, 3.0f, "", "", "");
}
int gem_hycal_matching(const char *evio_path, const char *out_path,
                         long max_events, int run_num, float match_nsigma)
{
    return gem_hycal_matching(evio_path, out_path,
                                "", "", "", max_events, run_num, match_nsigma,
                                "", "", "");
}

// Entry point — full version
int gem_hycal_matching(const char *evio_path,
                         const char *out_path,
                         const char *gem_ped_file,
                         const char *gem_cm_file,
                         const char *hc_calib_file,
                         long        max_events,
                         int         run_num,
                         float       match_nsigma,
                         const char *daq_config,
                         const char *gem_map_file,
                         const char *hc_map_file)
{
    Printf("[gem_hycal_matching] entered: evio=%s out=%s",
           evio_path ? evio_path : "(null)",
           out_path  ? out_path  : "(null)");

    //---- detector pipeline (DAQ config, runinfo, HyCal, GEM) ----------------
    prad2::Pipeline pipeline;
    if (!build_script_pipeline(pipeline, evio_path, run_num, daq_config,
                               hc_calib_file, gem_ped_file, gem_cm_file,
                               hc_map_file, gem_map_file))
        return 1;

    auto &cfg            = pipeline.daq_cfg;
    const auto crate_map = cfg.roc_crate_map();
    auto &geo            = pipeline.run_cfg;
    auto &hycal          = pipeline.hycal;
    auto &gem_sys        = pipeline.gem;
    auto &gem_pos_res    = pipeline.gem_pos_res;
    auto &hycal_xform    = pipeline.hycal_transform;
    auto &gem_xforms     = pipeline.gem_transforms;

    // Mirror into gRunConfig for any out-of-tree code that reads the
    // analysis-side global.
    analysis::gRunConfig = geo;

    // HyCal cluster + GEM cluster + WaveAnalyzer all hold per-event scratch,
    // so they live outside the pipeline (one per event loop, not one per
    // run).  HyCalCluster needs the wired-up HyCalSystem reference.
    fdec::HyCalCluster hc_clusterer(hycal);
    hc_clusterer.SetConfig(pipeline.hycal_cluster_cfg);
    gem::GemCluster gem_clusterer;

    //---- EVIO discovery -----------------------------------------------------
    auto evio_files = prad2::discover_split_files(evio_path ? evio_path : "");
    if (evio_files.empty()) {
        Printf("[ERROR] no EVIO files found for %s", evio_path ? evio_path : "(null)");
        return 1;
    }
    Printf("[setup] EVIO       : %zu split file(s) for input %s",
           evio_files.size(), evio_path ? evio_path : "(null)");
    for (const auto &f : evio_files) Printf("           %s", f.c_str());
    Printf("[setup] Match cut  : %.2f · sigma_total", match_nsigma);

    //---- ROOT output --------------------------------------------------------
    TFile fout(out_path, "RECREATE");
    if (fout.IsZombie()) {
        Printf("[ERROR] cannot create %s", out_path);
        return 1;
    }
    TTree *tree = new TTree("match",
        "HyCal+GEM clusters with straight-line matches and constituent strip waveforms");
    EventVars ev;
    bind_branches(tree, ev);

    //---- event loop ---------------------------------------------------------
    auto t0 = std::chrono::steady_clock::now();
    fdec::WaveAnalyzer ana;

    long n_filled = 0;
    long total_clusters = 0, total_matches = 0, total_strips = 0;
    long total_gem_2d   = 0;
    long gem_2d_per_det[4] = {0, 0, 0, 0};

    EvioScanStats scan;
    for_each_physics_event(cfg, evio_files, max_events, scan,
                           [&](const fdec::EventData &fadc_evt, const ssp::SspEventData &ssp_evt) {
        ev.clear();
        ev.event_num    = static_cast<int>(fadc_evt.info.event_number);
        ev.trigger_bits = fadc_evt.info.trigger_bits;

        // ---------- HyCal: waveform → energy → clusters ----------
        const auto hc_hits_raw = reconstruct_hycal_event(
            fadc_evt, crate_map, hycal, ana, hc_clusterer,
            geo.hc_time_win_lo, geo.hc_time_win_hi);

        // detector frame → lab/target-centered frame, at shower depth
        std::vector<analysis::HCHit> hc_hits;
        hc_hits.reserve(hc_hits_raw.size());
        for (const auto &h : hc_hits_raw)
            hc_hits.push_back(analysis::ClusterToLab(hycal_xform, h));

        // ---------- GEM: pedestal → CM → ZS → 1D + 2D ----------
        gem_sys.Clear();
        gem_sys.ProcessEvent(ssp_evt);
        gem_sys.Reconstruct(gem_clusterer);

        // Per-detector lab-frame hit lists for matching.
        std::vector<analysis::GEMHit> gem_lab[4];
        for (int d = 0; d < gem_sys.GetNDetectors() && d < 4; ++d) {
            const auto &raw = gem_sys.GetHits(d);
            gem_2d_per_det[d] += static_cast<long>(raw.size());
            total_gem_2d     += static_cast<long>(raw.size());
            for (const auto &h : raw)
                gem_lab[d].push_back(analysis::GemHitToLab(gem_xforms, h));
        }

        // ---------- record HyCal clusters in tree ----------
        ev.ncl = static_cast<int>(hc_hits.size());
        for (int k = 0; k < ev.ncl; ++k) {
            const auto &h = hc_hits[k];
            ev.hc_x.push_back(h.x);
            ev.hc_y.push_back(h.y);
            ev.hc_z.push_back(h.z);
            ev.hc_energy.push_back(h.energy);
            ev.hc_center.push_back(h.center_id);
            ev.hc_nblocks.push_back(hc_hits_raw[k].nblocks);
            ev.hc_flag.push_back(h.flag);
            // σ_HC at the HyCal face (mm) — see HyCalSystem::PositionResolution
            ev.hc_sigma.push_back(hycal.PositionResolution(h.energy));
        }

        // ---------- best-match per HC cluster × GEM detector ------------
        // Geometry and best-match rule: see the file header.
        for (int k = 0; k < ev.ncl; ++k) {
            const auto &h = hc_hits[k];
            if (h.z <= 0.f) continue;
            const float sigma_face = ev.hc_sigma[k];

            for (int d = 0; d < 4; ++d) {
                const auto &gl = gem_lab[d];
                if (gl.empty()) continue;
                // GEM z (lab) is the same for every hit on plane d —
                // just read the first one we have.
                const float z_gem = gl.front().z;
                if (z_gem <= 0.f) continue;
                const float scale  = z_gem / h.z;
                const float proj_x = h.x * scale;
                const float proj_y = h.y * scale;
                const float sig_hc_at_gem = sigma_face * scale;
                const float sig_gem = (d < (int)gem_pos_res.size())
                                        ? gem_pos_res[d] : 0.1f;
                const float sig_total = std::sqrt(
                    sig_hc_at_gem * sig_hc_at_gem + sig_gem * sig_gem);
                const float cut = match_nsigma * sig_total;

                // Find the closest GEM hit on detector d — must be
                // within `cut`.  best_dr starts at cut so any candidate
                // outside the window is automatically rejected.
                int   best_gi = -1;
                float best_dr = cut;
                for (size_t gi = 0; gi < gl.size(); ++gi) {
                    const float dx = gl[gi].x - proj_x;
                    const float dy = gl[gi].y - proj_y;
                    const float dr = std::sqrt(dx*dx + dy*dy);
                    if (dr <= best_dr) {
                        best_dr = dr;
                        best_gi = static_cast<int>(gi);
                    }
                }
                if (best_gi < 0) continue;

                const auto &g       = gl[best_gi];
                const auto &raw_g   = gem_sys.GetHits(d)[best_gi];
                const auto *xc = find_constituent(
                    gem_sys.GetPlaneClusters(d, 0),
                    raw_g.x, raw_g.x_charge);
                const auto *yc = find_constituent(
                    gem_sys.GetPlaneClusters(d, 1),
                    raw_g.y, raw_g.y_charge);

                const int mi = ev.nmatch;
                ev.m_hc_idx.push_back(k);
                ev.m_det.push_back(d);
                ev.m_gem_x.push_back(g.x);
                ev.m_gem_y.push_back(g.y);
                ev.m_gem_z.push_back(g.z);
                ev.m_gem_x_charge.push_back(raw_g.x_charge);
                ev.m_gem_y_charge.push_back(raw_g.y_charge);
                ev.m_gem_x_size.push_back(raw_g.x_size);
                ev.m_gem_y_size.push_back(raw_g.y_size);
                ev.m_proj_x.push_back(proj_x);
                ev.m_proj_y.push_back(proj_y);
                ev.m_residual.push_back(best_dr);
                ev.m_sigma_total.push_back(sig_total);

                if (xc) {
                    ev.m_xcl_position.push_back(xc->position);
                    ev.m_xcl_total.push_back(xc->total_charge);
                    ev.m_xcl_peak.push_back(xc->peak_charge);
                    ev.m_xcl_max_tb.push_back(xc->max_timebin);
                } else {
                    ev.m_xcl_position.push_back(0.f);
                    ev.m_xcl_total.push_back(0.f);
                    ev.m_xcl_peak.push_back(0.f);
                    ev.m_xcl_max_tb.push_back(-1);
                }
                if (yc) {
                    ev.m_ycl_position.push_back(yc->position);
                    ev.m_ycl_total.push_back(yc->total_charge);
                    ev.m_ycl_peak.push_back(yc->peak_charge);
                    ev.m_ycl_max_tb.push_back(yc->max_timebin);
                } else {
                    ev.m_ycl_position.push_back(0.f);
                    ev.m_ycl_total.push_back(0.f);
                    ev.m_ycl_peak.push_back(0.f);
                    ev.m_ycl_max_tb.push_back(-1);
                }

                auto xs = append_strips(ev, mi, 0, xc);
                ev.m_xcl_first.push_back(xs.first);
                ev.m_xcl_nstrips.push_back(xs.second);
                auto ys = append_strips(ev, mi, 1, yc);
                ev.m_ycl_first.push_back(ys.first);
                ev.m_ycl_nstrips.push_back(ys.second);

                ++ev.nmatch;
            }
        }

        tree->Fill();
        ++n_filled;
        total_clusters += ev.ncl;
        total_matches  += ev.nmatch;
        total_strips   += ev.nstrips;
    });

    auto t1 = std::chrono::steady_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();

    fout.cd();
    tree->Write();
    fout.Close();

    print_scan_summary(scan, evio_files.size());
    Printf("  tree entries written  : %ld", n_filled);
    Printf("  total HyCal clusters  : %ld", total_clusters);
    Printf("  total GEM 2D hits     : %ld  (det0=%ld det1=%ld det2=%ld det3=%ld)",
           total_gem_2d,
           gem_2d_per_det[0], gem_2d_per_det[1],
           gem_2d_per_det[2], gem_2d_per_det[3]);
    Printf("  total matches         : %ld", total_matches);
    Printf("  total strip rows      : %ld", total_strips);
    Printf("  elapsed (s)           : %.2f", secs);
    Printf("  wrote                 : %s", out_path);
    return 0;
}
