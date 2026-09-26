// hycal_shower_profile.cpp
// Reads raw replay files (prad_<run>.*_raw.root from replay_rawdata), rebuilds
// HyCal clusters event by event, and compares the seed and neighbour module
// energies of single-cluster Mott events, whose matched GEM hit lies near the
// seed module centre, with the shower-profile projection of the cluster.  The
// peak branches (replay_rawdata -p) are used when present; otherwise the peaks
// are re-derived from the waveform samples.

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
#include "RunInfoConfig.h"
#include "gain_factor.h"
#include "PipelineBuilder.h"
#include "ToolUtils.h"

#include <TFile.h>
#include <TH1F.h>
#include <TChain.h>

#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <vector>
#include <memory>
#include <cstdio>
#include <mutex>
#include <cmath>

using EventVars = prad2::RawEventData;
using namespace analysis;

int main(int argc, char *argv[])
{
    std::string db_dir = prad2::database_dir();

    std::string output_path_name, daq_config_file, recon_config_file, gem_ped_file;
    int  max_events  = -1;
    int  num_threads = 4;
    int  num_files   = -1;
    bool worker_mode = false;

    int opt;
    while ((opt = getopt(argc, argv, "o:n:f:j:w")) != -1) {
        switch (opt) {
            case 'o': output_path_name = optarg; break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'f': num_files        = std::atoi(optarg); break;
            case 'j': num_threads     = std::atoi(optarg); break;
            case 'w': worker_mode      = true; break;
        }
    }

    std::vector<std::string> root_files =
        CollectInputs(argc, argv, optind, IsRawRootName, num_files);
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

    // ROOT global state and TChain branch addresses are not thread-safe: run
    // one isolated worker process per input file and merge their histograms.
    if (!worker_mode && root_files.size() > 1) {
        std::vector<std::string> worker_outputs(root_files.size());
        std::mutex io_mutex;
        ParallelFor(root_files.size(), num_threads, [&](size_t i, int) {
            const std::string worker_prefix =
                output_path_name + ".worker_" + std::to_string(i);
            worker_outputs[i] = worker_prefix + ".root";
            std::vector<std::string> args{argv[0], "-w", "-o", worker_prefix};
            if (max_events > 0) args.insert(args.end(), {"-n", std::to_string(max_events)});
            args.insert(args.end(), {"-j", "1", root_files[i]});
            const int rc = RunCommand(args);
            std::lock_guard<std::mutex> lock(io_mutex);
            if (rc != 0)
                std::cerr << "Worker failed for " << root_files[i]
                          << " (exit code " << rc << ")\n";
        });
        if (!MergeTopLevelHistograms(worker_outputs, output_path_name + ".root")) return 1;
        for (const auto &f : worker_outputs) std::remove(f.c_str());
        return 0;
    }

    TChain tree("events");
    for (const auto &file : root_files) {
        tree.Add(file.c_str());
    }

    auto ev = std::make_unique<EventVars>();
    prad2::SetRawReadBranches(&tree, *ev);
    const bool has_waveform = tree.GetBranch("hycal.samples") != nullptr;
    const bool has_peaks    = tree.GetBranch("hycal.npeaks") != nullptr;

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
    recon_config_file = db_dir + "/reconstruction_config.json";
    daq_config_file = db_dir + "/daq_config.json";

    // Detectors: PRad-II flows through PipelineBuilder so the wiring stays in
    // one place (see prad2det/include/PipelineBuilder.h).
    prad2::Pipeline pipeline = prad2::PipelineBuilder()
        .set_database_dir(db_dir)
        .set_recon_config(recon_config_file)
        .set_daq_config(daq_config_file)
        .set_gem_pedestal(gem_ped_file)     // empty falls back to RunConfig default
        .set_run_number(run_num)
        .set_log_stream(&std::cerr)
        .build();
    const auto &hycal        = pipeline.hycal;
    const auto &gem_sys      = pipeline.gem;
    const auto &hc_time_cuts = pipeline.hycal_time_cuts;
    const auto &hc_xform     = pipeline.hycal_transform;
    const auto &g_xform      = pipeline.gem_transforms;

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(pipeline.hycal_cluster_cfg);
    gem::GemCluster      gem_clusterer;
    MatchingTools        matching(pipeline.match_method);
    matching.Configure(gRunConfig);

    fdec::WaveAnalyzer ana(pipeline.daq_cfg.wave_cfg);
    fdec::WaveResult wres;

    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(gRunConfig, run_num);

    TH1F *h1_cluster_energy = new TH1F("h1_cluster_energy", "Cluster Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_seed_energy = new TH1F("h1_seed_energy", "Seed Module Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_neighbor_energy = new TH1F("h1_neighbor_energy", "Neighbor Module Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_seed_project = new TH1F("h1_seed_project", "Seed Module Projected Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_neighbor_project = new TH1F("h1_neighbor_project", "Neighbor Module Projected Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH2F *h2_pos = new TH2F("h2_pos_live", "Hit Position;X_d[20.75mm];Y_d[20.77mm]", 40, -1, 1, 40, -1, 1);
    TH1F *h1_seed_fraction = new TH1F("h1_seed_fraction", "Seed Module Energy Fraction;Fraction;Counts", 100, 0, 1);
    TH1F *h1_seed_fraction_proj = new TH1F("h1_seed_fraction_proj", "Seed Module Energy Fraction Projected;Fraction;Counts", 100, 0, 1);
    TH1F *h1_neighbor_fraction = new TH1F("h1_neighbor_fraction", "Neighbor Module Energy Fraction;Fraction;Counts", 100, 0, 1);
    TH1F *h1_neighbor_fraction_proj = new TH1F("h1_neighbor_fraction_proj", "Neighbor Module Energy Fraction Projected;Fraction;Counts", 100, 0, 1);
    TH1F *h1_neighbor_energy2 = new TH1F("h1_neighbor_energy2", "Neighbor Module 2 Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_neighbor_project2 = new TH1F("h1_neighbor_project2", "Neighbor Module 2 Projected Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_neighbor_fraction2 = new TH1F("h1_neighbor_fraction2", "Neighbor Module 2 Energy Fraction;Fraction;Counts", 100, 0, 1);
    TH1F *h1_neighbor_fraction_proj2 = new TH1F("h1_neighbor_fraction_proj2", "Neighbor Module 2 Energy Fraction Projected;Fraction;Counts", 100, 0, 1);

    long long nentries = tree.GetEntries();
    for (long long i = 0; i < nentries; ++i) {
        tree.GetEntry(i);
        if (i >= max_events && max_events > 0) break;
        if (i % 10000 == 0) std::cout << "Processed " << i << " / " << nentries << " entries.\r" << std::flush;

        // assume you are selecting Mott-like events with a single cluster in the HyCal
        if ((ev->trigger_bits & prad2::TBIT_sum) == 0) continue;
        if (ev->nch > 70) continue; // channel numbers too high, likely not a clean event

        clusterer.Clear();

        // Per-event gain correction (time-series lookup by event number).
        const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(ev->event_num));

        // raw trees written without peak branches: re-derive the peaks from the samples
        if (has_waveform && !has_peaks) FillPeaksFromWaveforms(*ev, hycal, ana, wres);

        for (int j = 0; j < ev->nch; ++j) {
            const auto *mod = hycal.module_by_id(ev->module_id[j]);
            if (!mod || !mod->is_pwo4()) continue;

            const float gain = gain_corr.ModuleGain(mod->id);
            float time_offset = mod->time_offset;

            const auto hc_win = hc_time_cuts.at(mod->index);
            // Multi-pulse mode: push every peak inside the trigger
            // window into the clusterer; the seed-anchored timing
            // coincidence cut (ClusterConfig::seed_time_window) is applied
            // inside HyCalCluster.
            for (int p = 0; p < ev->npeaks[j]; ++p) {
                float peak_time = ev->peak_time[j][p] - time_offset; // apply module time offset
                if (peak_time <= hc_win.lo) continue;
                if (peak_time >= hc_win.hi) continue;
                float adc = ev->peak_integral[j][p] * gain;
                float energy = static_cast<float>(mod->energize(adc));
                clusterer.AddHit(mod->index, energy, peak_time);
            }
        }
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);

        // No reconstructed cluster means there is no HyCal position to use
        // for the event selection below.
        if (hits.empty()) continue;

        HCHit hc_hit;
        GEMHit g_hit;
        hc_hit.x = hits[0].x;
        hc_hit.y = hits[0].y;
        hc_hit.z = gRunConfig.hycal_z;

        // Event Selection only use HyCal
        // select single cluster Mott events
        if (hits.size() != 1 || hits[0].nblocks < 3) continue;
        float theta = PhysicsTools::GetThetaAngle(hc_hit.x, hc_hit.y, hc_hit.z);
        if (std::abs(hits[0].energy - gRunConfig.Ebeam) > 3.f * hycal.EnergyResolution(gRunConfig.Ebeam)) continue;
        if (theta < 0.7 ) continue;
        if (fdec::test_bit(hits[0].flag, fdec::kSplit)) continue; // skip clusters with split hits

        int seed_id = hits[0].center_id;
        const auto *seed_mod = hycal.module_by_id(seed_id);
        if (!seed_mod || !seed_mod->is_pwo4()) continue;
        const auto *neighbor_mod = hycal.module_by_id(seed_mod->id + 1);
        if (!neighbor_mod) continue;
        const auto *neighbor_mod2 = hycal.module_by_id(seed_mod->id + 2);
        if (!neighbor_mod2) continue;

        if (!InHyCalRing(seed_mod->x, seed_mod->y, 4.0, 14.)) continue;

        // Here reconstruct the GEM parts do matching with HyCal clusters
        if (gem_sys.GetNDetectors() > 0) {
            std::vector<gem::GEMHit> all_gem_hits;
            ReconstructGemStrips(*ev, gem_sys, gem_clusterer, all_gem_hits);

            std::vector<HCHit> hc_hits;
            std::vector<GEMHit> gem_hits[4];
            for (const auto &hit : hits) {
                HCHit local_hit = ClusterToLab(hc_xform, hit);
                GetProjection(local_hit, gRunConfig.hycal_z);
                hc_hits.push_back(local_hit);
            }
            for (const auto &hit : all_gem_hits)
                if (hit.det_id >= 0 && hit.det_id < 4)
                    gem_hits[hit.det_id].push_back(GemHitToLab(g_xform, hit));
            const auto matched_hits = matching.Match(
                hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]);
            if (matched_hits.empty()) continue;
            // use the upstream GEM3/GEM4 pair hit (gem[1], det_id 2/3)
            if (matched_hits[0].gem[1].det_id < 2
                || matched_hits[0].gem[1].det_id > 3)
                continue;
            g_hit.x = matched_hits[0].gem[1].x;
            g_hit.y = matched_hits[0].gem[1].y;
            g_hit.z = matched_hits[0].gem[1].z;
        }

        // projection of the GEM hit onto the HyCal plane
        if (g_hit.z == 0.f) continue;
        GetProjection(g_hit, gRunConfig.hycal_z);

        //move back to HyCal coordinate system
        ApplyToHyCal(hc_hit, gRunConfig);
        ApplyToHyCal(g_hit, gRunConfig);

        // require the GEM hit near the seed module centre (|xd|,|yd| < 0.2 module sizes)
        const auto [xd, yd] = seed_mod->cell_offset<float>(g_hit.x, g_hit.y);
        h2_pos->Fill(xd, yd);
        if (std::abs(xd) >= 0.2f || std::abs(yd) >= 0.2f) continue;

        float seed_energy = 0.f;
        float neighbor_energy = 0.f;
        float neighbor_energy2 = 0.f;

        for (const auto &cluster : clusterer.GetClusters()) {
            if (cluster.center.index == seed_mod->index) {
                seed_energy = cluster.center.energy;
                for (const auto &hit : cluster.hits){
                    if (hit.index == neighbor_mod->index)
                        neighbor_energy = hit.energy;
                    if (hit.index == neighbor_mod2->index)
                        neighbor_energy2 = hit.energy;
                }
            }
        }

        // projected energy for seed and neighbor modules from the clusterer's shower profile
        const auto projected_energy = [&](const fdec::Module *mod) {
            return hits[0].energy * clusterer.ProfileFractionAt(
                hc_hit.x, hc_hit.y, hits[0].energy, mod->index);
        };
        const float seed_proj      = projected_energy(seed_mod);
        const float neighbor_proj  = projected_energy(neighbor_mod);
        const float neighbor_proj2 = projected_energy(neighbor_mod2);
        h1_seed_project->Fill(seed_proj);
        h1_neighbor_project->Fill(neighbor_proj);
        h1_neighbor_project2->Fill(neighbor_proj2);

        h1_cluster_energy->Fill(hits[0].energy);
        h1_seed_energy->Fill(seed_energy);
        h1_neighbor_energy->Fill(neighbor_energy);
        h1_neighbor_energy2->Fill(neighbor_energy2);
        h1_seed_fraction->Fill(seed_energy / hits[0].energy);
        h1_seed_fraction_proj->Fill(seed_proj / hits[0].energy);
        h1_neighbor_fraction->Fill(neighbor_energy / hits[0].energy);
        h1_neighbor_fraction_proj->Fill(neighbor_proj / hits[0].energy);
        h1_neighbor_fraction2->Fill(neighbor_energy2 / hits[0].energy);
        h1_neighbor_fraction_proj2->Fill(neighbor_proj2 / hits[0].energy);
    }

    TFile *output_file = new TFile((output_path_name + ".root").c_str(), "RECREATE");
    h1_cluster_energy->Write();
    h1_seed_energy->Write();
    h1_neighbor_energy->Write();
    h1_seed_project->Write();
    h1_neighbor_project->Write();
    h1_seed_fraction->Write();
    h1_seed_fraction_proj->Write();
    h1_neighbor_fraction->Write();
    h1_neighbor_fraction_proj->Write();
    h1_neighbor_energy2->Write();
    h1_neighbor_project2->Write();
    h1_neighbor_fraction2->Write();
    h1_neighbor_fraction_proj2->Write();
    h2_pos->Write();
    output_file->Close();

}