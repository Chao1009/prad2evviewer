// Check the leakage correction effect.  Seed module W565 (id 1565), neighbour
// W566 (id 1566): for elastic two-cluster events with the hit near the seed
// centre, compare the measured seed/neighbour module energies with the
// shower-profile projection of the reconstructed cluster.

const int seed_id = 1565;
const int neighbor_id = 1566;
float seed_energy, neighbor_energy;

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
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
#include <utility>
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

    int opt;
    while ((opt = getopt(argc, argv, "o:n:f:j:")) != -1) {
        switch (opt) {
            case 'o': output_path_name = optarg; break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'f': num_files        = std::atoi(optarg); break;
            case 'j': num_threads     = std::atoi(optarg); break;
        }
    }

    std::vector<std::string> root_files = CollectInputs(argc, argv, optind, IsRawRootName, num_files);
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: leakage_correction_check <input_raw.root|dir> [more...] "
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
    const bool has_peaks    = tree.GetBranch("hycal.npeaks") != nullptr;

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
    recon_config_file = db_dir + "/reconstruction_config.json";
    daq_config_file = db_dir + "/daq_config.json";

    prad2::Pipeline pipeline = prad2::PipelineBuilder()
        .set_database_dir(db_dir)
        .set_recon_config(recon_config_file)
        .set_daq_config(daq_config_file)
        .set_gem_pedestal(gem_ped_file)     // empty falls back to RunConfig default
        .set_run_number(run_num)
        .set_log_stream(&std::cerr)
        .build();

    const auto &hycal        = pipeline.hycal;
    const auto &hc_time_cuts = pipeline.hycal_time_cuts;

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(pipeline.hycal_cluster_cfg);

    fdec::WaveAnalyzer ana(pipeline.daq_cfg.wave_cfg);
    fdec::WaveResult wres;

    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(gRunConfig, run_num);

    TH1F *h1_cluster_energy = new TH1F("h1_cluster_energy", "Cluster Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_seed_energy = new TH1F("h1_seed_energy", "Seed Module Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_neighbor_energy = new TH1F("h1_neighbor_energy", "Neighbor Module Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_seed_project = new TH1F("h1_seed_project", "Seed Module Projected Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_neighbor_project = new TH1F("h1_neighbor_project", "Neighbor Module Projected Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH2F *h2_pos = new TH2F("h2_pos_live", "Hit Position;X_d[20.75mm];Y_d[20.77mm]", 40, -1, 1, 40, -1, 1);

    long long nentries = tree.GetEntries();
    for (long long i = 0; i < nentries; ++i) {
        tree.GetEntry(i);
        if (i >= max_events && max_events > 0) break;
        if (i % 10000 == 0) std::cout << "Processed " << i << " / " << nentries << " entries.\r" << std::flush;

        if ((ev->trigger_bits & prad2::TBIT_sum) == 0) continue;
        if (ev->nch > 70) continue; // channel numbers too high, likely not a clean event

        clusterer.Clear();

        seed_energy = 0.f;
        neighbor_energy = 0.f;

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
                    if (mod->id == neighbor_id) {
                        neighbor_energy = energy;
                    }
                    if (mod->id == seed_id) {
                        seed_energy = energy;
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
                    if (mod->id == neighbor_id) {
                        neighbor_energy = energy;
                    }
                    if (mod->id == seed_id) {
                        seed_energy = energy;
                    }
                }
            }
        }
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);

        // two-cluster elastic events, one of them seeded at seed_id
        if (hits.size() != 2 || hits[0].nblocks < 3 || hits[1].nblocks < 3) continue;
        if (hits[0].center_id != seed_id && hits[1].center_id != seed_id) continue;
        if (fabs(hits[0].energy + hits[1].energy - gRunConfig.Ebeam) > 3.f * hycal.EnergyResolution(gRunConfig.Ebeam)) continue;

        float phi1 = std::atan2(hits[0].y, hits[0].x) * 180.0 / M_PI;
        float phi2 = std::atan2(hits[1].y, hits[1].x) * 180.0 / M_PI;
        float dphi = std::abs(std::fabs(phi1 - phi2) - 180.f);
        if (dphi > 10.f) continue; // require back-to-back clusters
        

        const auto *seed_mod = hycal.module_by_id(seed_id);
        const auto *neighbor_mod = hycal.module_by_id(neighbor_id);
        if (!seed_mod || !neighbor_mod) continue;

        if (hits[1].center_id == seed_id) {
            std::swap(hits[0], hits[1]);
        }
        if (hits[0].energy < 1600. || hits[0].energy > 1850.) continue;

        if (fdec::test_bit(hits[0].flag, fdec::kSplit)) continue; // skip clusters with split hits

        // require the hit near the seed module centre (|xd|,|yd| < 0.2 module sizes)
        const auto [xd, yd] = seed_mod->cell_offset<float>(hits[0].x, hits[0].y);
        h2_pos->Fill(xd, yd);
        if (std::abs(xd) >= 0.2f || std::abs(yd) >= 0.2f) continue;

        // projected energy for seed and neighbor modules from the clusterer's shower profile
        const auto projected_energy = [&](const fdec::Module *mod) {
            return hits[0].energy * clusterer.ProfileFractionAt(hits[0].x, hits[0].y, hits[0].energy, mod->index);
        };
        h1_seed_project->Fill(projected_energy(seed_mod));
        h1_neighbor_project->Fill(projected_energy(neighbor_mod));

        h1_cluster_energy->Fill(hits[0].energy);
        h1_seed_energy->Fill(seed_energy);
        h1_neighbor_energy->Fill(neighbor_energy);
    }

    TFile *output_file = new TFile((output_path_name + ".root").c_str(), "RECREATE");
    h1_cluster_energy->Write();
    h1_seed_energy->Write();
    h1_neighbor_energy->Write();
    h1_seed_project->Write();
    h1_neighbor_project->Write();
    h2_pos->Write();
    output_file->Close();

}