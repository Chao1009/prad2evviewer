//Check the leakage correction effect
// Choose module W566 as seed module, in 5 by 5 region no dead modules,
// then set 567 as dead module, then is 568 as dead module
// then check the energy beafore and after the dead module, and check the leakage correction effect

const int seed_id = 1565;
const int dead_id = 1566;
float seed_energy, dead_energy;

#include "Replay.h"
#include "PhysicsTools.h"
#include "MatchingTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "GemSystem.h"
#include "WaveAnalyzer.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "load_daq_config.h"
#include "RunInfoConfig.h"
#include "gain_factor.h"
#include "PulseTemplateStore.h"
#include "PipelineBuilder.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2Poly.h>
#include <TChain.h>
#include <TFitResult.h>
#include <TFile.h>
#include <TText.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <filesystem>
#include <vector>
#include <memory>
#include <algorithm>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <cmath>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>

#include <nlohmann/json.hpp>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;

using EventVars = prad2::RawEventData;
using namespace analysis;

// ── File collection helper ───────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find("_raw.root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

// ── Main ─────────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR")) db_dir = env;

    // ── Argument parsing ─────────────────────────────────────────────────────
    std::string output_path_name, daq_config_file, recon_config_file, gem_ped_file;
    int  max_events  = -1;
    int  num_threads = 4;
    int  num_files   = -1;

    int opt;
    while ((opt = getopt(argc, argv, "o:n:f:j:")) != -1) {
        switch (opt) {
            case 'o': output_path_name = optarg; break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'f': num_files        = std::atoi(optarg); break;
            case 'j': num_threads     = std::atoi(optarg); break;
        }
    }

    // Collect all input files
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectRootFiles(argv[i]);
        if (num_files > 0) {
            int remaining = num_files - static_cast<int>(root_files.size());
            if (remaining <= 0) break;
            int take = std::min(remaining, static_cast<int>(f.size()));
            root_files.insert(root_files.end(), f.begin(), f.begin() + take);
            if (static_cast<int>(root_files.size()) >= num_files) break;
        } else {
            root_files.insert(root_files.end(), f.begin(), f.end());
        }
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: hycal_shower_profile <input_raw.root|dir> [more...] "
                     "[-o ./output_name(no extension)] [-n max_events] [-f nfiles] [-j threads]\n";
        return 1;
    }

    if (output_path_name.empty()) {
        std::cerr << "No output prefix provided. Please pass -o <output_prefix>.\n";
        return 1;
    }

    TChain tree("events");
    for (const auto &file : root_files) {
        tree.Add(file.c_str());
    }

    auto ev = std::make_unique<EventVars>();
    prad2::SetRawReadBranches(&tree, *ev);
    const bool has_waveform = tree.GetBranch("hycal.samples") != nullptr;

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
    recon_config_file = db_dir + "/reconstruction_config.json";
    daq_config_file = db_dir + "/daq_config.json";

    evc::DaqConfig daq_cfg;

    // Detectors: PRad-II flows through PipelineBuilder so the wiring stays in
    // one place (see prad2det/include/PipelineBuilder.h).
    fdec::HyCalSystem                 hycal;
    gem::GemSystem                    gem_sys;
    fdec::ClusterConfig               cluster_cfg;
    prad2::HyCalTimeCuts              hc_time_cuts;
    prad2::HyCalRfOffsets             hc_rf_offsets;
    DetectorTransform                 hycal_transform;
    std::array<DetectorTransform, 4>  gem_transforms;
    std::unordered_map<int, int>      roc_to_crate;
    int                               match_method = 1;

    prad2::Pipeline pipeline = prad2::PipelineBuilder()
        .set_database_dir(db_dir)
        .set_recon_config(recon_config_file)
        .set_daq_config(daq_config_file)
        .set_gem_pedestal(gem_ped_file)     // empty falls back to RunConfig default
        .set_run_number(run_num)
        .set_log_stream(&std::cerr)
        .build();

    daq_cfg          = std::move(pipeline.daq_cfg);
    hycal            = std::move(pipeline.hycal);
    gem_sys          = std::move(pipeline.gem);
    cluster_cfg      = pipeline.hycal_cluster_cfg;
    hc_time_cuts     = std::move(pipeline.hycal_time_cuts);
    hc_rf_offsets    = std::move(pipeline.hycal_rf_offsets);
    hycal_transform  = pipeline.hycal_transform;
    gem_transforms   = pipeline.gem_transforms;
    match_method     = pipeline.match_method;

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);
    gem::GemCluster      gem_clusterer;
    MatchingTools        matching(match_method);

    //initialize tools for cluster reconstruction
    fdec::WaveAnalyzer ana(daq_cfg.wave_cfg);
    fdec::PulseTemplateStore template_store;
    if (daq_cfg.wave_cfg.nnls_deconv.enabled
        && !daq_cfg.wave_cfg.nnls_deconv.template_file.empty()) {
        template_store.LoadFromFile(
            db_dir + "/" + daq_cfg.wave_cfg.nnls_deconv.template_file,
            daq_cfg.wave_cfg);
    }
    ana.SetTemplateStore(&template_store);
    fdec::WaveResult wres;

    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(
        gRunConfig.gain_data_dir + "/gain_correction", run_num);

    // create histograms you want to fill for shower profile analysis
    TH1F *h1_cluster_energy = new TH1F("h1_cluster_energy", "Cluster Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_seed_fraction = new TH1F("h1_seed_fraction", "Seed Module Fraction;Fraction;Counts", 100, 0, 1);
    TH1F *h1_neighbor_fraction = new TH1F("h1_neighbor_fraction", "Fraction of Energy in Neighbor Module;Fraction;Counts", 100, 0, 1);
    TH2F *h2_pos = new TH2F("h2_pos_live", "Hit Position;X_d[20.75mm];Y_d[20.77mm]", 40, -1, 1, 40, -1, 1);

    // Here loop over the events in the TChain, read channels data, reconstruct clusters, and fill the histograms
    long long nentries = tree.GetEntries();
    for (long long i = 0; i < nentries; ++i) {
        tree.GetEntry(i);
        if (i >= max_events && max_events > 0) break;
        if (i % 10000 == 0) std::cout << "Processed " << i << " / " << nentries << " entries.\r" << std::flush;

        // assume you are selecting Mott-like events with a single cluster in the HyCal
        if ((ev->trigger_bits & prad2::TBIT_sum) == 0) continue;
        if (ev->nch > 70) continue; // channel numbers too high, likely not a clean event

        // Reconstruct clusters for this event.
        clusterer.Clear();

        dead_energy = 0.f;
        seed_energy = 0.f;

        // Per-event gain correction (time-series lookup by event number).
        const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(ev->event_num));

        for (int j = 0; j < ev->nch; ++j) {
            const auto *mod = hycal.module_by_id(ev->module_id[j]);
            if (!mod || !mod->is_pwo4()) continue;

            // Per-ID gain correction: average of three LMS channels.
            const float gain = (mod->id > 1000)
                ? (gain_corr.w[mod->id - 1000].corr[1] + gain_corr.w[mod->id - 1000].corr[2]) / 2.0f
                : 1.0f; // default gain factor
            // timing offset for this module
            float time_offset = mod->time_offset;

            const auto hc_win = hc_time_cuts.at(mod->index);
            // Multi-pulse mode: push every peak inside the trigger
            // window into the clusterer; the seed-anchored timing
            // coincidence cut is applied inside HyCalCluster, 
            // maximum 4 ns difference between pulses in a same cluster.
            if (has_waveform) 
            {
                ana.Analyze(ev->samples[j], ev->nsamples[j], wres, time_offset);
                for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                    const auto &pk = wres.peaks[p];
                    if (pk.time <= hc_win.lo) continue;
                    if (pk.time >= hc_win.hi) continue;
                    float adc = pk.integral * gain;
                    float energy = static_cast<float>(mod->energize(adc));
                    clusterer.AddHit(mod->index, energy, pk.time);
                    if (mod->id == dead_id) {
                        dead_energy = energy; // store the energy of the dead module for later analysis
                    }
                    if (mod->id == seed_id) {
                        seed_energy = energy; // store the energy of the seed module for later analysis
                    }
                }
            }
            else
            {
                for (int p = 0; p < ev->npeaks[j]; ++p) {
                    float peak_time = ev->peak_time[j][p] - time_offset; // apply module time offset
                    if (peak_time <= hc_win.lo) continue;
                    if (peak_time >= hc_win.hi) continue;
                    float adc = ev->peak_integral[j][p] * gain;
                    float energy = static_cast<float>(mod->energize(adc));
                    clusterer.AddHit(mod->index, energy, peak_time);
                }
            }
        }
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);

        // select single cluster Mott events
        //std::cout << hits.size() << " " << hits[0].nblocks << " " << hits[0].energy << " " << std::atan2(std::sqrt(hits[0].x * hits[0].x + hits[0].y * hits[0].y), gRunConfig.hycal_z) * 180.0 / M_PI << std::endl;
        if (hits.size() != 2 || hits[0].nblocks < 3 || hits[1].nblocks < 3) continue;
        if (hits[0].center_id != seed_id && hits[1].center_id != seed_id) continue;
        if (fabs(hits[0].energy + hits[1].energy - gRunConfig.Ebeam) > 3. * 0.033 * std::sqrt(gRunConfig.Ebeam * 1000.)) continue;
        float hc_x = hits[0].x, hc_y = hits[0].y, hc_z = gRunConfig.hycal_z;
        float theta = std::atan2(std::sqrt(hc_x * hc_x + hc_y * hc_y), hc_z) * 180.0 / M_PI;

        float phi1 = std::atan2(hits[0].y, hits[0].x) * 180.0 / M_PI;
        float phi2 = std::atan2(hits[1].y, hits[1].x) * 180.0 / M_PI;
        float dphi = std::abs(std::fabs(phi1 - phi2) - 180.f);
        if (dphi > 10.f) continue; // require back-to-back clusters
        

        const auto *seed_mod = hycal.module_by_id(seed_id);
        if (!seed_mod || !seed_mod->is_pwo4()) continue;

        if (hits[1].center_id == seed_id) {
            std::swap(hits[0], hits[1]);
        }

        // require hit to be in central 3x3 of a 5x5 grid in single central module (|xd|,|yd| < 0.3)
        float xd = (hits[0].x - (float)seed_mod->x) / (float)seed_mod->size_x;
        float yd = (hits[0].y - (float)seed_mod->y) / (float)seed_mod->size_y;
        h2_pos->Fill(xd, yd);
        if (std::abs(xd) >= 0.3f || std::abs(yd) >= 0.3f) continue;

        float dead_fraction = dead_energy / hits[0].energy;
        if (cluster_cfg.leakage_correction) {
            dead_fraction = hits[0].leakage / hits[0].energy; // use leakage energy if correction is applied
        }
        float seed_fraction = seed_energy / hits[0].energy;

        h1_cluster_energy->Fill(hits[0].energy);
        h1_seed_fraction->Fill(seed_fraction);
        h1_neighbor_fraction->Fill(dead_fraction);
    }

    // Save the histograms to a root file
    TFile *output_file = new TFile((output_path_name + ".root").c_str(), "RECREATE");
    h1_cluster_energy->Write();
    h1_seed_fraction->Write();
    h1_neighbor_fraction->Write();
    h2_pos->Write();
    output_file->Close();

}