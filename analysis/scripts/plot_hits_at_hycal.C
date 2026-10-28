//============================================================================
// plot_hits_at_hycal.C — 2D occupancy of GEM hits (projected to HyCal
// surface) and HyCal cluster centroids (already on HyCal surface), drawn
// side-by-side in the lab / target-centered, beam-aligned frame at z =
// hycal_z.
//
// Pipeline per physics event:
//   EvChannel.Read()  → DecodeEvent() → FADC + SSP buffers
//                     → HyCal: WaveAnalyzer → energize → HyCalCluster
//                     → GEM:   GemSystem.ProcessEvent → Reconstruct
//                     → coord transform to lab (per-detector tilt + offset
//                       via prad2det's DetectorTransform)
//                     → GEM hits: GetProjection(hits, hycal_z) — straight
//                       line from target through (x,y,z) to z=hycal_z
//                     → fill the two TH2F occupancy maps
//
// Both plots share x/y range and binning.  GEM hits from all four
// detectors are combined into the left histogram; HyCal cluster centroids
// (one entry per cluster) populate the right histogram.
//
// Trigger filter: only events with `trigger_bits == 0x100` (production
// physics trigger) contribute.  Everything else (LMS / Alpha / cosmic /
// etc.) is skipped.
//
// Multi-file mode is selected by the input path:
//   * `/data/.../prad_023881.evio.*`  → glob: enumerate every sibling
//     `prad_023881.evio.<digits>`, fold them all into the same two
//     histograms, and warn (to stderr) about any gap in the suffix
//     sequence (including missing from .00000).
//   * `/data/prad_023881/`            → directory: same enumeration,
//     run number sniffed from the directory name.
//   * `/data/.../prad_023881.evio.00000` → single specific split file.
//
// Usage
// -----
//   cd build
//   root -l ../analysis/scripts/rootlogon.C
//
//   # full run (glob — warns about any missing split):
//   .x ../analysis/scripts/plot_hits_at_hycal.C+( \
//       "/data/stage6/prad_023867/prad_023867.evio.*", \
//       "hits_at_hycal.pdf")
//
//   # single split (debugging):
//   .x ../analysis/scripts/plot_hits_at_hycal.C+( \
//       "/data/stage6/prad_023867/prad_023867.evio.00000", \
//       "hits_at_hycal_seg0.pdf")
//
//   args (full): evio_path, out_path, max_events, run_num,
//                gem_ped_file, gem_cm_file, hc_calib_file,
//                daq_config, gem_map_file, hc_map_file
//   - out_path  : PDF/PNG/etc. for the canvas; an accompanying .root
//                 file alongside it stores both TH2Fs for re-plotting.
//   - max_events: 0 = all
//   - run_num   : -1 = sniff from EVIO basename (prad_NNNNNN.evio.*)
//   - all "_file"/"daq_config" args: "" = auto-discover via runinfo
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
#include "MatchingTools.h"    // ClusterToLab, GemHitToLab, GetProjection
#include "PipelineBuilder.h"  // prad2::Pipeline
#include "script_helpers.h"   // build_script_pipeline, for_each_physics_event, strip_extension

#include <TCanvas.h>
#include <TFile.h>
#include <TH2F.h>
#include <TStyle.h>
#include <TString.h>

#include <chrono>
#include <memory>
#include <string>
#include <vector>

// Full 10-arg version + convenience overloads (cling default-arg
// marshalling SEGVs for mixed-type signatures).
int plot_hits_at_hycal(const char *evio_path,
                       const char *out_path,
                       long        max_events,
                       int         run_num,
                       const char *gem_ped_file,
                       const char *gem_cm_file,
                       const char *hc_calib_file,
                       const char *daq_config,
                       const char *gem_map_file,
                       const char *hc_map_file);

int plot_hits_at_hycal(const char *evio_path, const char *out_path)
{
    return plot_hits_at_hycal(evio_path, out_path,
                              0L, -1, "", "", "", "", "", "");
}
int plot_hits_at_hycal(const char *evio_path, const char *out_path,
                       long max_events)
{
    return plot_hits_at_hycal(evio_path, out_path,
                              max_events, -1, "", "", "", "", "", "");
}
int plot_hits_at_hycal(const char *evio_path, const char *out_path,
                       long max_events, int run_num)
{
    return plot_hits_at_hycal(evio_path, out_path,
                              max_events, run_num, "", "", "", "", "", "");
}

// Entry point — full version
int plot_hits_at_hycal(const char *evio_path,
                       const char *out_path,
                       long        max_events,
                       int         run_num,
                       const char *gem_ped_file,
                       const char *gem_cm_file,
                       const char *hc_calib_file,
                       const char *daq_config,
                       const char *gem_map_file,
                       const char *hc_map_file)
{
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
    auto &hycal_xform    = pipeline.hycal_transform;
    auto &gem_xforms     = pipeline.gem_transforms;
    analysis::gRunConfig = geo;

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

    //---- histograms ---------------------------------------------------------
    // Lab frame (target-centered, beam-aligned) at z = hycal_z.
    // Range +/-650 mm covers PRad-II HyCal LG ring (~580 mm to outer edge);
    // 5 mm bins give a clean occupancy map without being too noisy.
    constexpr float kRange = 650.f;
    constexpr int   kBins  = 260;        // 5 mm bins
    auto h_gem = std::make_unique<TH2F>(
        "h_gem_at_hycal",
        TString::Format(
            "GEM hits projected to HyCal surface (z = %.0f mm);x (mm);y (mm)",
            geo.hycal_z),
        kBins, -kRange, kRange, kBins, -kRange, kRange);
    auto h_hc = std::make_unique<TH2F>(
        "h_hycal",
        TString::Format(
            "HyCal cluster centroids on HyCal surface (z = %.0f mm);x (mm);y (mm)",
            geo.hycal_z),
        kBins, -kRange, kRange, kBins, -kRange, kRange);

    //---- event loop ---------------------------------------------------------
    auto t0 = std::chrono::steady_clock::now();
    fdec::WaveAnalyzer ana;
    long n_hc_clusters = 0, n_gem_hits = 0;

    EvioScanStats scan;
    for_each_physics_event(cfg, evio_files, max_events, scan,
                           [&](const fdec::EventData &fadc_evt, const ssp::SspEventData &ssp_evt) {
        // ---------- HyCal: waveform → energy → clusters ----------
        const auto hc_raw = reconstruct_hycal_event(
            fadc_evt, crate_map, hycal, ana, hc_clusterer, 100.f, 200.f);

        // Build HCHit list with z = 0 (no shower depth) so transform
        // lands them at exactly z = hycal_z — i.e. on the HyCal face.
        std::vector<analysis::HCHit> hc_hits;
        hc_hits.reserve(hc_raw.size());
        for (const auto &h : hc_raw)
            hc_hits.push_back(analysis::ClusterToLab(hycal_xform, h, false));

        for (const auto &h : hc_hits) h_hc->Fill(h.x, h.y);
        n_hc_clusters += hc_hits.size();

        // ---------- GEM: pedestal → CM → ZS → 1D + 2D ----------
        gem_sys.Clear();
        gem_sys.ProcessEvent(ssp_evt);
        gem_sys.Reconstruct(gem_clusterer);

        // Per-detector lab-frame hit lists.  GetHits(d) returns local
        // plane hits (x, y, z=0); rotate + transform per-detector,
        // then project the line target->hit onto z = hycal_z.
        for (int d = 0; d < gem_sys.GetNDetectors() && d < 4; ++d) {
            const auto &raw = gem_sys.GetHits(d);
            if (raw.empty()) continue;
            std::vector<analysis::GEMHit> lab;
            lab.reserve(raw.size());
            for (const auto &h : raw) lab.push_back(analysis::GemHitToLab(gem_xforms, h));
            analysis::GetProjection(lab, geo.hycal_z);
            for (const auto &g : lab) h_gem->Fill(g.x, g.y);
            n_gem_hits += lab.size();
        }
    });
    auto t1 = std::chrono::steady_clock::now();
    double secs = std::chrono::duration<double>(t1 - t0).count();

    //---- draw + save --------------------------------------------------------
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    gStyle->SetNumberContours(99);

    TCanvas c("c_hits_at_hycal", "Hits at HyCal surface", 1600, 720);
    c.Divide(2, 1, 0.005, 0.005);

    c.cd(1);
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.10);
    h_gem->Draw("COLZ");

    c.cd(2);
    gPad->SetRightMargin(0.13);
    gPad->SetLeftMargin(0.10);
    h_hc->Draw("COLZ");

    c.SaveAs(out_path);

    // Sibling .root file with the two histograms for re-plotting.
    std::string root_out = strip_extension(out_path) + ".root";
    TFile fout(root_out.c_str(), "RECREATE");
    if (!fout.IsZombie()) {
        h_gem->Write();
        h_hc->Write();
        c.Write();
        fout.Close();
        Printf("[setup] Saved hists: %s", root_out.c_str());
    }

    print_scan_summary(scan, evio_files.size());
    Printf("  HyCal clusters total  : %ld  (avg %.2f / kept event)",
           n_hc_clusters,
           scan.n_kept ? double(n_hc_clusters) / scan.n_kept : 0.0);
    Printf("  GEM hits total (4 det): %ld  (avg %.2f / kept event)",
           n_gem_hits,
           scan.n_kept ? double(n_gem_hits) / scan.n_kept : 0.0);
    Printf("  elapsed (s)           : %.2f", secs);
    Printf("  wrote canvas          : %s", out_path);
    return 0;
}
