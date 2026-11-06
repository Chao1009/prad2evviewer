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
#include <TProfile2D.h>

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

// ── Main ─────────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    std::string db_dir = prad2::database_dir();

    // ── Argument parsing ─────────────────────────────────────────────────────
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

    // create histograms you want to fill for shower profile analysis
    float ratioMin = 0.0;
    float ratioMax = 1.0;
    float ratioBins = (ratioMax - ratioMin)*1000.0;

    float eMin = 0;
    float eMax = 4000;
    float energyBins = (eMax - eMin);

    TH1F *h1_cluster_energy = new TH1F("h1_cluster_energy", "e-p Cluster Energy;Energy [MeV];Counts", energyBins, eMin, eMax);
    TH1F *h1_seed_energy = new TH1F("h1_seed_energy", "e-p Seed Module Energy;Energy [MeV];Counts", energyBins, eMin, eMax);
    
    TH1F *h1_seedOverClust = new TH1F("h1_seedOverClust", "e-p Seed Energy over CLuster Energy; Energy Ratio; Counts;", ratioBins, ratioMin, ratioMax);

    TH1F *h1_seedNFirstOverClust = new TH1F("h1_seedNFirstOverClust", "e-p Seed + First Layer Energy over Cluster Energy;Energy Ratio; Counts;", ratioBins, ratioMin, ratioMax);
    TH1F *h1_seedNFirstNSecOverClust = new TH1F("h1_seedNFirstNSecOverClust", "e-p Seed + First and Second Layer Energy over Cluster Energy;Energy Ratio; Counts;", ratioBins, ratioMin, ratioMax);
    
    TH1F *h1_FirstNSec_NoSeed_OverClust = new TH1F("h1_FirstNSec_NoSeed_OverClust", "e-p No Seed, First and Second Layer Energy over Cluster Energy;Energy Ratio; Counts;", ratioBins, ratioMin, ratioMax);
    TH1F *h1_First_NoSeed_OverClust = new TH1F("h1_First_NoSeed_OverClust", "e-p No Seed, First Layer Energy over Cluster Energy;Energy Ratio; Counts;", ratioBins, ratioMin, ratioMax);
    
    TH2F *h2_lowSeedE_XYPos = new TH2F("h2_lowSeedE_XYPos", 
        "Seed Module Energy < 0.5 of Total XY Position;X [mm];Y [mm]", 
        800, -400, 400, 800, -400, 400);
    
    TH2F *h2_cluster_vs_seed_energy = new TH2F("h2_cluster_vs_seed_energy", 
        "Cluster Energy vs Seed Module Energy;Seed Energy [MeV];Cluster Energy [MeV]", 
        4000, 0, 4000, 4000, 0, 4000);

    TH1F *h1_blocks_1st_layer = new TH1F("h1_blocks_1st_layer", "Number of Fired Modules in 1st Layer;Count;Events", 21, -0.5, 20.5);
    TH1F *h1_blocks_2nd_layer = new TH1F("h1_blocks_2nd_layer", "Number of Fired Modules in 2nd Layer;Count;Events", 21, -0.5, 20.5);
    TH1F *h1_time_diff_seed = new TH1F("h1_time_diff_seed", "Time Difference with Seed Module;#Delta t [ns];Events", 200, -10, 10);

    float max_x = 10.5;
    int xBins = 11;
    float minE = 50;
    float maxE = gRunConfig.Ebeam;
    int eBins = (int) (gRunConfig.Ebeam/1000);

    TString labels[4] = {"_right", "_left", "_top", "_bot"};
    TString sub[4] = {"R", "L", "T", "B"};
    TString ax[4] = {"x", "x", "y", "y"};
    TProfile2D *hProf2D[4];

    for(int q = 0; q < 4; q++){
        hProf2D[q] = new TProfile2D("hProf2D"+labels[q], 
            "2D Profile Ratio Map (E_{c}/E_{"+sub[q]+"}) for each Cluster Energy;Impact Position from Centroid "+ax[q]+" [mm];Cl_{E} [MeV];#LT E_{C} / E_{"+sub[q]+"} #GT",
            xBins, -1*max_x, max_x,   // Position bins
            eBins, minE, maxE   // Energy bins
        );
    }

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

        // Per-event gain correction (time-series lookup by event number).
        const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(ev->event_num));

        // raw trees written without peak branches: re-derive the peaks from the samples
        if (has_waveform && !has_peaks) FillPeaksFromWaveforms(*ev, hycal, ana, wres);

        for (int j = 0; j < ev->nch; ++j) {
            const auto *mod = hycal.module_by_id(ev->module_id[j]);
            if (!mod || !mod->is_pwo4()) continue;

            const float gain = gain_corr.ModuleGain(mod->id, gRunConfig.gain_ref_pmt);
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

        // Here could also reconstruct the GEM parts do matching with HyCal clusters if needed
        // TODO: Add GEM cluster reconstruction and matching if needed

        // select single cluster Mott events - Not Necessary anymore for getting profile at various energies.
        //if (hits.size() != 1 || hits[m].nblocks < 3) continue;

        for(int m = 0; m < hits.size(); m++){
            if(hits[m].nblocks < 3) continue;

            float hc_x = hits[m].x, hc_y = hits[m].y, hc_z = gRunConfig.hycal_z;
            if (!InHyCalRing(hc_x, hc_y, 3.0, 15.)) continue;

            float theta = std::atan2(std::sqrt(hc_x * hc_x + hc_y * hc_y), hc_z) * 180.0 / M_PI;
            if (std::abs(hits[m].energy - gRunConfig.Ebeam) > 3. * 0.03 * std::sqrt(gRunConfig.Ebeam * 1000.)) continue;
            if (theta < 0.7 ) continue; 

            int seed_id = hits[m].center_id;
            const auto *seed_mod = hycal.module_by_id(seed_id);
            if (!seed_mod || !seed_mod->is_pwo4()) continue;
            float sx = seed_mod->x, sy = seed_mod->y; //Get the centroid position of the seed module.

            float dx = hc_x - sx, dy = hc_y - sy;

            //the time window cut is already applied in the clusterer
            float dt_max = pipeline.hycal_cluster_cfg.seed_time_window;

            // counts for fired modules in the 1st and 2nd layer of neighbors
            int n1st_layer = 0;
            int n2nd_layer = 0;

            bool seed_not_clean = false;

            // loop over the hits in the cluster
            // Here go to the event by event analysis, add more analysis code or selection .etc here
            float energyTot = 0;
            float energyCent = 0;
            float energyCent_First = 0;
            float energyFirst = 0;
            float energyFirst_Sec = 0;
            float energySec = 0;
            bool lowSeed = false;

            float top = 0;
            float left = 0;
            float left_mid = 0;
            float mid_row = 0;
            float mid_col = 0;
            float bot = 0;
            float right = 0;
            float right_mid = 0;
            for (int j = 0; j < ev->nch; ++j) {
                const auto *mod = hycal.module_by_id(ev->module_id[j]);
                if (!mod || !mod->is_pwo4()) continue;

                // Per-ID gain correction: average of three LMS channels.
                const float gain = gain_corr.ModuleGain(mod->id, gRunConfig.gain_ref_pmt);
                // timing offset for this module
                float time_offset = mod->time_offset;
                const int drow = mod->row - seed_mod->row;
                const int dcolumn = mod->column - seed_mod->column;
                // select the seed module and its energy (require only one peak in seed module)
                if (mod->id == seed_id) {
                    if (ev->npeaks[j] != 1) {
                        seed_not_clean = true;
                        break;
                    }
                    float seed_energy = static_cast<float>(mod->energize(ev->peak_integral[j][0] * gain));
                    energyTot += seed_energy;
                    energyCent += seed_energy;
                    energyCent_First += seed_energy;
                    h1_seed_energy->Fill(seed_energy);
                    h2_cluster_vs_seed_energy->Fill(seed_energy, hits[m].energy);

                    if(seed_energy/hits[m].energy < 0.5){
                        lowSeed = true;
                    }
                }
                // select the 1st layer of neighboring modules
                if (std::abs(drow) <= 1 && std::abs(dcolumn) <= 1 && mod->id != seed_mod->id) {
                    // check if this module has a peak in seed time window
                    for (int p = 0; p < ev->npeaks[j]; ++p) {
                        float peak_time = ev->peak_time[j][p] - time_offset; // apply module time offset
                        float dt = peak_time - hits[m].time;
                        h1_time_diff_seed->Fill(dt);
                        if (std::abs(dt) < dt_max) {
                            n1st_layer++;
                            float energy = static_cast<float>(mod->energize(ev->peak_integral[j][p] * gain));
                            energyTot += energy;
                            energyCent_First += energy;
                            energyFirst += energy;
                            energyFirst_Sec += energy;
                                
                            // fill the histogram for shower profile analysis

                            // Classify the first-ring crystals by row and column.
                            if (dcolumn == 1 && std::abs(drow) <= 1) {
                                right += energy;
                                if (drow == 0) right_mid += energy;
                            }
                            if (dcolumn == -1 && std::abs(drow) <= 1) {
                                left += energy;
                                if (drow == 0) left_mid += energy;
                            }
                            if (drow == -1 && std::abs(dcolumn) <= 1) {
                                top += energy;
                            }
                            if (drow == 1 && std::abs(dcolumn) <= 1) {
                                bot += energy;
                            }
                            if (dcolumn == 0 && std::abs(drow) == 1) {
                                mid_col += energy;
                            }
                            if (drow == 0 && std::abs(dcolumn) == 1) {
                                mid_row += energy;
                            }
                            break; // only count one peak per module
                        }
                    }
                }
                // select the 2nd layer of neighboring modules
                if (std::abs(drow) <= 2 && std::abs(dcolumn) <= 2
                    && (std::abs(drow) == 2 || std::abs(dcolumn) == 2) && mod->id != seed_mod->id) {
                    
                    // check if this module has a peak in seed time window
                    for (int p = 0; p < ev->npeaks[j]; ++p) {
                        float peak_time = ev->peak_time[j][p] - time_offset; // apply module time offset
                        float dt = peak_time - hits[m].time;
                        h1_time_diff_seed->Fill(dt);
                        if (std::abs(dt) < dt_max) {
                            n2nd_layer++;
                            float energy = static_cast<float>(mod->energize(ev->peak_integral[j][p] * gain));
                            energyTot += energy;
                            energySec += energy;
                            energyFirst_Sec += energy;
                            // fill the histogram for shower profile analysis

                            break; // only count one peak per module
                        }
                    }
                }
            }
        
            // fill the histograms for the number of fired modules in the 1st and 2nd layers
            if (!seed_not_clean) {  
                if (right > 0) hProf2D[0]->Fill(dx, hits[m].energy, (energyCent+mid_col)/right);
                if (left > 0) hProf2D[1]->Fill(dx, hits[m].energy, (energyCent+mid_col)/left);
                if (top > 0) hProf2D[2]->Fill(dy, hits[m].energy, (energyCent+mid_row)/top);
                if (bot > 0) hProf2D[3]->Fill(dy, hits[m].energy, (energyCent+mid_row)/bot);

                h1_blocks_1st_layer->Fill(n1st_layer);
                h1_blocks_2nd_layer->Fill(n2nd_layer);
                h1_cluster_energy->Fill(hits[m].energy);

                h1_seedOverClust->Fill(energyCent/hits[m].energy);
                h1_seedNFirstOverClust->Fill(energyCent_First/hits[m].energy);
                h1_seedNFirstNSecOverClust->Fill(energyTot/hits[m].energy);

                h1_First_NoSeed_OverClust->Fill(energyFirst/hits[m].energy);
                h1_FirstNSec_NoSeed_OverClust->Fill(energyFirst_Sec/hits[m].energy);

                if(lowSeed){h2_lowSeedE_XYPos->Fill(hits[m].x, hits[m].y);}
            }
        }
    }

    // Save the histograms to a root file
    TFile *output_file = new TFile((output_path_name + ".root").c_str(), "RECREATE");
    for(int g = 0; g < 4; g++){
        hProf2D[g]->Write();
    }
    h1_cluster_energy->Write();
    h1_seed_energy->Write();
    h1_seedOverClust->Write();
    h1_seedNFirstOverClust->Write();
    h1_seedNFirstNSecOverClust->Write();
    h1_FirstNSec_NoSeed_OverClust->Write();
    h1_First_NoSeed_OverClust->Write();
    h2_lowSeedE_XYPos->Write();
    h2_cluster_vs_seed_energy->Write();
    h1_blocks_1st_layer->Write();
    h1_blocks_2nd_layer->Write();
    h1_time_diff_seed->Write();

    output_file->Close();

    return 0;
}