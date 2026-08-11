
//=============================================================================
// det_align — detector position calibration via Møller scattering
//
// Process: Read reconstructed root files, judge if it's 
// Moller events, analyze and fill histgroams, find out the detector alignment
//
// Usage:
//   det_align <recon_root_file_or_dir> [more files/dirs...]
//             -o output_dir
//=============================================================================

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"
#include "MatchingTools.h"
#include "PipelineBuilder.h"
#include "PulseTemplateStore.h"
#include "gain_factor.h"

#include <TClass.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

using EventVars_Recon = prad2::ReconEventData;

struct EventWithMoller {
    int event_num = -1;
    MollerEvent HC_moller;                    // Møller event reconstructed from HyCal
    MollerEvent GEMup_moller;                 // Møller event reconstructed from upstream GEM
    MollerEvent GEMdown_moller;               // Møller event reconstructed from downstream GEM
    MollerEvent GEM_moller[4] = {};           // Møller events on same chamber from GEM1, GEM2, GEM3, GEM4
};
struct EventWithMott {
    int event_num = -1;
    float HC_x = 0, HC_y = 0, HC_z = 0;
    float GEM_x[4] = {}, GEM_y[4] = {}, GEM_z[4] = {};
    bool match[4] = {false, false, false, false};
};

std::vector<EventWithMoller> AllMollerEvents;
std::vector<EventWithMott> AllMottEvents;

// ── forward declarations ──────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path);

static double extract_peak(TH1F *hist);

static bool ProcessFile(const std::string &input_root,
                        const std::string &db_dir,
                        const RunConfig &gRunConfig, const RunConfig &in_run_config,
                        std::vector<EventWithMoller> &all_moller_events, std::vector<EventWithMott> &all_mott_events,
                        bool show_progress = false,
                        int max_events = -1);


// ── Main ──────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");

    std::string  output_dir, run_config_in, run_config_out_path;
    int  max_files   = -1;
    int  num_threads = 4;
    int  iter        = 1;
    int  max_events     = -1;

    int opt;
    while ((opt = getopt(argc, argv, "o:f:j:c:r:i:n:")) != -1) {
        switch (opt) {
            case 'o': output_dir       = optarg; break;
            case 'f': max_files        = std::atoi(optarg); break;
            case 'j': num_threads      = std::atoi(optarg); break;
            case 'c': run_config_in       = optarg; break;
            case 'r': run_config_out_path       = optarg; break;
            case 'i': iter             = std::atoi(optarg); break;
            case 'n': max_events         = std::atoi(optarg); break;
            default: return 1;
        }
    }

    // collect input replay_recon files (files, directories, or mixed)
    std::vector<std::string> recon_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectRootFiles(argv[i]);
        recon_files.insert(recon_files.end(), f.begin(), f.end());
    }

    if (recon_files.empty() || output_dir.empty()) {
        std::cerr <<
            "Usage: det_calib <recon_file_or_dir> [more files/dirs...] -o output_dir\n"
            "       [-f max_files] [-j threads] [-i iteration] [-n max_events]\n"
            "       [-c run_config_in.json] [-r run_config_base_path.json]\n";
        return 1;
    }

    int num_files = static_cast<int>(recon_files.size());
    if (max_files > 0) num_files = std::min(num_files, max_files);
    num_threads = std::max(1, std::min(num_threads, num_files));

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    std::string run_config = db_dir + "/runinfo/general.json";

    int run_num = get_run_int(recon_files[0]);
    gRunConfig = LoadRunConfig(run_config, run_num);
    
    // Set up config paths with automatic iteration naming
    // If -r specified, use it as base path; otherwise use output_dir/run_config.json
    std::string config_base_path;
    if (!run_config_out_path.empty()) {
        config_base_path = run_config_out_path;
    } else {
        config_base_path = (fs::path(output_dir) / "run_config.json").string();
    }
    
    // Generate input config path for this iteration
    if (iter == 1) {
        if (run_config_in.empty()) {
            run_config_in = run_config;  // use general.json by default
        }
    } else {
        // Use previous iteration's output config as input
        fs::path base(config_base_path);
        std::string stem = base.stem().string();
        std::string ext = base.extension().string();
        if (ext.empty()) ext = ".json";  // Default to .json if no extension
        std::string prev_iter_config = (base.parent_path() / 
                                        (stem + "_iter" + std::to_string(iter - 1) + ext)).string();
        run_config_in = prev_iter_config;
        std::cout << "Iteration " << iter << ": using output from iteration " << (iter - 1) << "\n";
        std::cout << "  Input:  " << run_config_in << "\n";
    }

    RunConfig in_run_config = LoadRunConfig(run_config_in, run_num);

    // Each worker writes only to the result slot belonging to its replay_recon file.
    // This avoids locking the (potentially large) Moller vectors while events
    // are being reconstructed.  The main thread merges the slots after all
    // workers have finished.
    std::vector<std::vector<EventWithMoller>> events_per_file(num_files);
    std::vector<std::vector<EventWithMott>> events_per_file_mott(num_files);
    std::vector<char> processed_ok(num_files, 0);
    std::atomic<int> next_file{0};
    std::atomic<int> errors{0};
    std::mutex io_mtx;

    auto worker = [&]() {
        while (true) {
            const int idx = next_file.fetch_add(1);
            if (idx >= num_files) break;

            const bool ok = ProcessFile(recon_files[idx],
                                        db_dir,
                                        gRunConfig, in_run_config,
                                        events_per_file[idx],
                                        events_per_file_mott[idx],
                                        idx == 0, max_events);  // show progress only for first file
            processed_ok[idx] = ok ? 1 : 0;

            std::lock_guard<std::mutex> lock(io_mtx);
            if (ok) {
                std::cout << "  [" << (idx + 1) << "/" << num_files << "] "
                          << recon_files[idx] << ": "
                          << events_per_file[idx].size() << " HyCal Mollers\n";
            } else {
                ++errors;
                std::cerr << "  [" << (idx + 1) << "/" << num_files
                          << "] FAILED: " << recon_files[idx] << "\n";
            }
        }
    };

    std::cout << "Processing " << num_files << " Recon_ROOT files with "
              << num_threads << " threads\n";
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i)
        threads.emplace_back(worker);
    for (auto &thread : threads)
        thread.join();

    // Aggregate in input-file order, so results are deterministic regardless
    // of the order in which worker threads completed.
    std::vector<EventWithMoller> all_events;
    std::vector<EventWithMott> all_mott_events;
    for (int i = 0; i < num_files; ++i) {
        if (processed_ok[i]) {
            all_events.insert(all_events.end(),
                              events_per_file[i].begin(),
                              events_per_file[i].end());
            all_mott_events.insert(all_mott_events.end(),
                                   events_per_file_mott[i].begin(),
                                   events_per_file_mott[i].end());
        }
    }

    std::cout << "Finished processing all files. Total HyCal Mollers: "
              << all_events.size() << "\n";

    // Create the histogram for saving the results
    // Hits position histograms
    TH2F *h2_hycal_hits = new TH2F("h2_hycal_hits", "HyCal Hits", 720, -360, 360, 720, -360, 360);
    TH2F *h2_gem_hits[4];
    for (int i = 0; i < 4; ++i)
        h2_gem_hits[i] = new TH2F(Form("h2_gem%d_hits", i), Form("GEM%d Hits", i), 720, -360, 360, 720, -360, 360);
    // Energy vs angle histograms
    TH2F *h2_hycal_energy_vs_angle = new TH2F("h2_hycal_energy_vs_angle", "HyCal Energy vs Angle", 35*2, 0.5, 4, 230*4, 0, 2300);
    TH2F *h2_gem_energy_vs_angle[4];
    for (int i = 0; i < 4; ++i)
        h2_gem_energy_vs_angle[i] = new TH2F(Form("h2_gem%d_energy_vs_angle", i), Form("GEM%d Energy vs Angle", i), 35*2, 0.5, 4, 230*4, 0, 2300);
    // Moller center position histograms
    TH2F *h2_hycal_Center = new TH2F("h2_hycal_mollerCenterX", "HyCal Moller Center X", 720, -360, 360, 720, -360, 360);
    TH2F *h2_gem_Center[4];
    for (int i = 0; i < 4; ++i)
        h2_gem_Center[i] = new TH2F(Form("h2_gem%d_mollerCenter", i), Form("GEM%d Moller Center", i), 720, -360, 360, 720, -360, 360);
    // Moller center X, Y, and vertex Z distance histograms
    TH1F *h1_hycal_CenterX   = new TH1F("h1_hycal_mollerCenterX",   "HyCal Moller Center X",   800, -50, 50);
    TH1F *h1_hycal_CenterY   = new TH1F("h1_hycal_mollerCenterY",   "HyCal Moller Center Y",   800, -50, 50);
    TH1F *h1_hycal_Zdistance = new TH1F("h1_hycal_mollerZdistance", "HyCal Moller Z distance", 1200, 5000, 8000);
    TH1F *h1_gem_CenterX[4], *h1_gem_CenterY[4], *h1_gem_Zdistance[4];
    for (int i = 0; i < 4; ++i) {
        h1_gem_CenterX[i]   = new TH1F(Form("h1_gem%d_mollerCenterX", i),   Form("GEM%d Moller Center X", i),   1600, -50, 50);
        h1_gem_CenterY[i]   = new TH1F(Form("h1_gem%d_mollerCenterY", i),   Form("GEM%d Moller Center Y", i),   200, -50, 50);
        h1_gem_Zdistance[i] = new TH1F(Form("h1_gem%d_mollerZdistance", i), Form("GEM%d Moller Z distance", i), 1200, 5000, 8000);
    }
    // Beam center measurement histograms
    TH2F *h2_beam_Center = new TH2F("h2_beam_Center", "Beam Center", 800, -20, 20, 800, -20, 20);
    TH1F *h1_beam_CenterX = new TH1F("h1_beam_CenterX", "Beam Center X", 800, -20, 20);
    TH1F *h1_beam_CenterY = new TH1F("h1_beam_CenterY", "Beam Center Y", 800, -20, 20);
    // Phi difference between HyCal and GEMs
    TH1F *h1_phi_diff_hycal_gem[4];
    for (int i = 0; i < 4; ++i)
        h1_phi_diff_hycal_gem[i] = new TH1F(Form("h1_phi_diff_hycal_gem%d", i), Form("Phi Difference HyCal-GEM%d", i), 400, -10, 10);
    // internal layer detector alignment histograms(2 chambers in the same layer, overlap region)
    TH1F *h1_deltaX_gem_up = new TH1F("h1_deltaX_gem_up", "Delta X GEM Up", 1000, -20, 20);
    TH1F *h1_deltaX_gem_down = new TH1F("h1_deltaX_gem_down", "Delta X GEM Down", 1000, -20, 20);
    TH1F *h1_deltaY_gem_up = new TH1F("h1_deltaY_gem_up", "Delta Y GEM Up", 1000, -20, 20);
    TH1F *h1_deltaY_gem_down = new TH1F("h1_deltaY_gem_down", "Delta Y GEM Down", 1000, -20, 20);
    // 2 layers GEM alignment histograms
    // check the residuals deltaX and deltaY between the two GEM layers
    TH1F *h1_deltaX_gem_layer_left = new TH1F("h1_deltaX_gem_layer_left", "Delta X GEM Layer Left", 1000, -20, 20);
    TH1F *h1_deltaY_gem_layer_left = new TH1F("h1_deltaY_gem_layer_left", "Delta Y GEM Layer Left", 1000, -20, 20);
    TH1F *h1_deltaX_gem_layer_right = new TH1F("h1_deltaX_gem_layer_right", "Delta X GEM Layer Right", 1000, -20, 20);
    TH1F *h1_deltaY_gem_layer_right = new TH1F("h1_deltaY_gem_layer_right", "Delta Y GEM Layer Right", 1000, -20, 20);
    TH2F *h2_deltaX_vs_deltaY_gem_left = new TH2F("h2_deltaX_vs_deltaY_gem_left", "Delta X vs Delta Y GEM Layer Left", 1000, -20, 20, 1000, -20, 20);
    TH2F *h2_deltaX_vs_deltaY_gem_right = new TH2F("h2_deltaX_vs_deltaY_gem_right", "Delta X vs Delta Y GEM Layer Right", 1000, -20, 20, 1000, -20, 20);
    // check the residuals deltaX and deltaY between each GEM layers and HyCal
    TH1F *h1_deltaX_gem_hycal[4];
    TH1F *h1_deltaY_gem_hycal[4];
    for (int i = 0; i < 4; ++i) {
        h1_deltaX_gem_hycal[i] = new TH1F(Form("h1_deltaX_gem_hycal%d", i), Form("Delta X GEM%d vs HyCal", i), 1000, -20, 20);
        h1_deltaY_gem_hycal[i] = new TH1F(Form("h1_deltaY_gem_hycal%d", i), Form("Delta Y GEM%d vs HyCal", i), 1000, -20, 20);
    }
    // ----- Up to here are just translational alignment histograms, next will be global roll histograms (pitch, yaw, roll)

    
    // analysis and filling of histograms will be done here
    for (int i = 0; i < all_events.size(); ++i) {
        EventWithMoller &thisEvent = all_events[i];
        // 1. HyCal moller events
        float hc_x[2] = {thisEvent.HC_moller.first.x, thisEvent.HC_moller.second.x};
        float hc_y[2] = {thisEvent.HC_moller.first.y, thisEvent.HC_moller.second.y};
        float hc_z[2] = {thisEvent.HC_moller.first.z, thisEvent.HC_moller.second.z};
        float hc_energy[2] = {thisEvent.HC_moller.first.E, thisEvent.HC_moller.second.E};
        float hc_z_distance = PhysicsTools::GetMollerZdistance(thisEvent.HC_moller, gRunConfig.Ebeam);
        float hc_theta[2] = {static_cast<float>(std::atan2(std::sqrt(hc_x[0]*hc_x[0] + hc_y[0]*hc_y[0]), hc_z[0]) * 180.0 / M_PI),
                             static_cast<float>(std::atan2(std::sqrt(hc_x[1]*hc_x[1] + hc_y[1]*hc_y[1]), hc_z[1]) * 180.0 / M_PI)};
        float hc_phi[2] = {static_cast<float>(std::atan2(hc_y[0], hc_x[0]) * 180.0 / M_PI), static_cast<float>(std::atan2(hc_y[1], hc_x[1]) * 180.0 / M_PI)};
        h2_hycal_hits->Fill(hc_x[0], hc_y[0]);
        h2_hycal_hits->Fill(hc_x[1], hc_y[1]);
        h2_hycal_energy_vs_angle->Fill(hc_theta[0], hc_energy[0]);
        h2_hycal_energy_vs_angle->Fill(hc_theta[1], hc_energy[1]);
        h1_hycal_Zdistance->Fill(hc_z_distance);
        if(i >= 3) {
            auto center = PhysicsTools::GetMollerCenter(all_events[i-1].HC_moller, thisEvent.HC_moller);
            if(center[0] != 0 || center[1] != 0) {
                h2_hycal_Center->Fill(center[0], center[1]);
                h1_hycal_CenterX->Fill(center[0]); h1_hycal_CenterY->Fill(center[1]);
            }
            center = PhysicsTools::GetMollerCenter(all_events[i-2].HC_moller, thisEvent.HC_moller);
            if(center[0] != 0 || center[1] != 0) {
                h2_hycal_Center->Fill(center[0], center[1]);
                h1_hycal_CenterX->Fill(center[0]); h1_hycal_CenterY->Fill(center[1]);
            }
            center = PhysicsTools::GetMollerCenter(all_events[i-3].HC_moller, thisEvent.HC_moller);
            if(center[0] != 0 || center[1] != 0) {
                h2_hycal_Center->Fill(center[0], center[1]);
                h1_hycal_CenterX->Fill(center[0]); h1_hycal_CenterY->Fill(center[1]);
            }
        }

        // 2. Use upstream GEMs to measure the beam center
        if (i >= 3){
            auto beam_center = PhysicsTools::GetMollerCenter(all_events[i-1].GEMup_moller, thisEvent.GEMup_moller);
            if(beam_center[0] != 0 || beam_center[1] != 0) {
                h2_beam_Center->Fill(beam_center[0], beam_center[1]);
                h1_beam_CenterX->Fill(beam_center[0]); h1_beam_CenterY->Fill(beam_center[1]);
            }
            beam_center = PhysicsTools::GetMollerCenter(all_events[i-2].GEMup_moller, thisEvent.GEMup_moller);
            if(beam_center[0] != 0 || beam_center[1] != 0) {
                h2_beam_Center->Fill(beam_center[0], beam_center[1]);
                h1_beam_CenterX->Fill(beam_center[0]); h1_beam_CenterY->Fill(beam_center[1]);
            }
            beam_center = PhysicsTools::GetMollerCenter(all_events[i-3].GEMup_moller, thisEvent.GEMup_moller);
            if(beam_center[0] != 0 || beam_center[1] != 0) {
                h2_beam_Center->Fill(beam_center[0], beam_center[1]);
                h1_beam_CenterX->Fill(beam_center[0]); h1_beam_CenterY->Fill(beam_center[1]);
            }
        }

        // 3. per GEM chamber moller events
        float gem_x[4][2], gem_y[4][2], gem_z[4][2], gem_energy[4][2];
        float gem_theta[4][2], gem_phi[4][2], gem_z_distance[4];
        for(int j = 0; j < 4; j++) {
            // Skip if GEM_moller[j] is not valid (no energy)
            if (thisEvent.GEM_moller[j].first.E <= 0 || thisEvent.GEM_moller[j].second.E <= 0) continue;
            
            gem_x[j][0] = thisEvent.GEM_moller[j].first.x;
            gem_x[j][1] = thisEvent.GEM_moller[j].second.x;
            gem_y[j][0] = thisEvent.GEM_moller[j].first.y;
            gem_y[j][1] = thisEvent.GEM_moller[j].second.y;
            gem_z[j][0] = thisEvent.GEM_moller[j].first.z;
            gem_z[j][1] = thisEvent.GEM_moller[j].second.z;
            gem_energy[j][0] = thisEvent.GEM_moller[j].first.E;
            gem_energy[j][1] = thisEvent.GEM_moller[j].second.E;
            gem_theta[j][0] = static_cast<float>(std::atan2(std::sqrt(gem_x[j][0]*gem_x[j][0] + gem_y[j][0]*gem_y[j][0]), gem_z[j][0]) * 180.0 / M_PI);
            gem_theta[j][1] = static_cast<float>(std::atan2(std::sqrt(gem_x[j][1]*gem_x[j][1] + gem_y[j][1]*gem_y[j][1]), gem_z[j][1]) * 180.0 / M_PI);
            gem_phi[j][0] = static_cast<float>(std::atan2(gem_y[j][0], gem_x[j][0]) * 180.0 / M_PI);
            gem_phi[j][1] = static_cast<float>(std::atan2(gem_y[j][1], gem_x[j][1]) * 180.0 / M_PI);
            gem_z_distance[j] = PhysicsTools::GetMollerZdistance(thisEvent.GEM_moller[j], gRunConfig.Ebeam);
            // Fill histograms for GEM chamber j
            h2_gem_hits[j]->Fill(gem_x[j][0], gem_y[j][0]);
            h2_gem_hits[j]->Fill(gem_x[j][1], gem_y[j][1]);
            h1_gem_Zdistance[j]->Fill(gem_z_distance[j]);
            h2_gem_energy_vs_angle[j]->Fill(gem_theta[j][0], gem_energy[j][0]);
            h2_gem_energy_vs_angle[j]->Fill(gem_theta[j][1], gem_energy[j][1]);
            // Moller center positions - find up to 3 previous events with valid GEM data for this chamber
            int count = 0;
            for (int k = i - 1; k >= 0 && count < 3; --k) {
                if (all_events[k].GEM_moller[j].first.E > 0 && all_events[k].GEM_moller[j].second.E > 0) {
                    auto center = PhysicsTools::GetMollerCenter(all_events[k].GEM_moller[j], thisEvent.GEM_moller[j]);
                    if(center[0] != 0 || center[1] != 0) {
                        h2_gem_Center[j]->Fill(center[0], center[1]);
                        h1_gem_CenterX[j]->Fill(center[0]); 
                        h1_gem_CenterY[j]->Fill(center[1]);
                    }
                    ++count;
                }
            }
            // Phi difference between HyCal and GEM
            for (int k = 0; k < 2; ++k) {
                float phi_diff = hc_phi[k] - gem_phi[j][k];
                if (phi_diff > 180.0) phi_diff -= 360.0;
                else if (phi_diff < -180.0) phi_diff += 360.0;
                h1_phi_diff_hycal_gem[j]->Fill(phi_diff);
            }
        }
        
    }

    // 4. Check the residuals between 2 layers of GEM chambers
    //    and check the residuals between GEM layers and HyCal
    // Calculate and fill residuals between the two layers of GEM chambers for each Mott event
    for (auto& event : all_mott_events) {
        if (event.match[0] && event.match[1] && event.match[2] && event.match[3]) {
            for (int det = 0; det < 4; ++det) {
                float scale = event.HC_z / event.GEM_z[det];
                event.GEM_x[det] *= scale;
                event.GEM_y[det] *= scale;
            }
            h1_deltaX_gem_up->Fill(event.GEM_x[2] - event.GEM_x[3]);
            h1_deltaY_gem_up->Fill(event.GEM_y[2] - event.GEM_y[3]);
            if(event.GEM_y[1] > 0.) {
                h1_deltaX_gem_down->Fill(event.GEM_x[0] - event.GEM_x[1]);
                h1_deltaY_gem_down->Fill(event.GEM_y[0] - event.GEM_y[1]);
            }
        }
        if (event.match[0] && event.match[2]) {
            for (int det : {0, 2}) {
                float scale = event.HC_z / event.GEM_z[det];
                event.GEM_x[det] *= scale;
                event.GEM_y[det] *= scale;
                h1_deltaX_gem_hycal[det]->Fill(event.GEM_x[det] - event.HC_x);
                h1_deltaY_gem_hycal[det]->Fill(event.GEM_y[det] - event.HC_y);
            }
            h1_deltaX_gem_layer_left->Fill(event.GEM_x[0] - event.GEM_x[2]);
            h1_deltaY_gem_layer_left->Fill(event.GEM_y[0] - event.GEM_y[2]);
            h2_deltaX_vs_deltaY_gem_left->Fill(event.GEM_x[0] - event.GEM_x[2], event.GEM_y[0] - event.GEM_y[2]);
        }
        if (event.match[1] && event.match[3]) {
            for (int det : {1, 3}) {
                float scale = event.HC_z / event.GEM_z[det];
                event.GEM_x[det] *= scale;
                event.GEM_y[det] *= scale;
                h1_deltaX_gem_hycal[det]->Fill(event.GEM_x[det] - event.HC_x);
                h1_deltaY_gem_hycal[det]->Fill(event.GEM_y[det] - event.HC_y);
            }
            h1_deltaX_gem_layer_right->Fill(event.GEM_x[1] - event.GEM_x[3]);
            h1_deltaY_gem_layer_right->Fill(event.GEM_y[1] - event.GEM_y[3]);
            h2_deltaX_vs_deltaY_gem_right->Fill(event.GEM_x[1] - event.GEM_x[3], event.GEM_y[1] - event.GEM_y[3]);
        }
    }

    // fit peak and resolve alignment parameters
    fs::create_directories(output_dir);
    gROOT->SetBatch(kTRUE);
    
    std::map<std::string, double> alignment_params;
    
    // Fit HyCal parameters
    alignment_params["moller_HC_x"] = extract_peak(h1_hycal_CenterX);
    alignment_params["moller_HC_y"] = extract_peak(h1_hycal_CenterY);
    alignment_params["moller_HC_z"] = extract_peak(h1_hycal_Zdistance);
    
    // Fit GEM and phi parameters
    for (int det = 0; det < 4; ++det) {
        alignment_params["moller_GEM" + std::to_string(det) + "_x"] = extract_peak(h1_gem_CenterX[det]);
        alignment_params["moller_GEM" + std::to_string(det) + "_y"] = extract_peak(h1_gem_CenterY[det]);
        alignment_params["moller_GEM" + std::to_string(det) + "_z"] = extract_peak(h1_gem_Zdistance[det]);
        alignment_params["moller_phi_diff_" + std::to_string(det)] = extract_peak(h1_phi_diff_hycal_gem[det]);
    }

    //  internal layer alignment parameters (same GEM layer left and right)
    alignment_params["Upstream_GEM_dx_d2-d3"] = extract_peak(h1_deltaX_gem_up);
    alignment_params["Upstream_GEM_dy_d2-d3"] = extract_peak(h1_deltaY_gem_up);
    alignment_params["Downstream_GEM_dx_d0-d1"] = extract_peak(h1_deltaX_gem_down);
    alignment_params["Downstream_GEM_dy_d0-d1"] = extract_peak(h1_deltaY_gem_down);
    // alignment parameters for 2 layers of GEM detectors (upstream and downstream)
    alignment_params["GEM_layer_left_dx_d0-d2"] = extract_peak(h1_deltaX_gem_layer_left);
    alignment_params["GEM_layer_left_dy_d0-d2"] = extract_peak(h1_deltaY_gem_layer_left);
    alignment_params["GEM_layer_right_dx_d1-d3"] = extract_peak(h1_deltaX_gem_layer_right);
    alignment_params["GEM_layer_right_dy_d1-d3"] = extract_peak(h1_deltaY_gem_layer_right);
    // Alignment parameters for GEM vs HyCal residuals
    for (int det = 0; det < 4; ++det) {
        alignment_params["deltaX_gem_hycal_" + std::to_string(det)] = extract_peak(h1_deltaX_gem_hycal[det]);
        alignment_params["deltaY_gem_hycal_" + std::to_string(det)] = extract_peak(h1_deltaY_gem_hycal[det]);
    }

    // Beam center measurement with upstream GEM
    alignment_params["beam_center_x"] = extract_peak(h1_beam_CenterX);
    alignment_params["beam_center_y"] = extract_peak(h1_beam_CenterY);
    std::cout << " ----- ***** ----- ***** ------ ***** ----- ******" << std::endl;
    std::cout << "The coordinates in this configuration: " << std::endl;
    std::cout << "Beam center (X, Y): " << alignment_params["beam_center_x"] << ", " 
              << alignment_params["beam_center_y"] << std::endl;
    std::cout << "The coordinates of HyCal center: " << std::endl;
    std::cout << "Beam center (X, Y): " << alignment_params["beam_center_x"] + in_run_config.target_x << ", "
              << alignment_params["beam_center_y"] + in_run_config.target_y << std::endl;
    // Target center Z measurement with upstream GEM
    alignment_params["target_center_z_gem3"] = in_run_config.gem_z[3] - alignment_params["moller_GEM3_z"];
    alignment_params["target_center_z_gem0"] = in_run_config.gem_z[0] - alignment_params["moller_GEM0_z"];
    alignment_params["target_center_z_gem1"] = in_run_config.gem_z[1] - alignment_params["moller_GEM1_z"];
    alignment_params["target_center_z_gem2"] = in_run_config.gem_z[2] - alignment_params["moller_GEM2_z"];
    std::cout << "Target center Z (GEM3): " << alignment_params["target_center_z_gem3"] << std::endl;
    std::cout << "Target center Z (GEM2): " << alignment_params["target_center_z_gem2"] << std::endl;
    std::cout << "Target center Z (GEM1): " << alignment_params["target_center_z_gem1"] << std::endl;
    std::cout << "Target center Z (GEM0): " << alignment_params["target_center_z_gem0"] << std::endl;
    std::cout << " ----- ***** ----- ***** ------ ***** ----- ******" << std::endl;

    // resolve the alignment parameters based on the extracted peak positions

    // the detector position was shifted by the target position when read from the configuration file
    // So we need to account for this shift, move the detector positions back by the target position
    in_run_config.hycal_x += in_run_config.target_x;
    in_run_config.hycal_y += in_run_config.target_y;
    in_run_config.hycal_z += in_run_config.target_z;
    for (int det = 0; det < 4; ++det) {
        in_run_config.gem_z[det] += in_run_config.target_z;
        in_run_config.gem_x[det] += in_run_config.target_x;
        in_run_config.gem_y[det] += in_run_config.target_y;
    }

    // 1. Extract beam position from HyCal X/Y alignment parameters
    float beam_x = alignment_params["moller_HC_x"];
    float beam_y = alignment_params["moller_HC_y"];
    // set target position based on beam position if not already set
    if(in_run_config.target_x == 0 && in_run_config.target_y == 0) {
        in_run_config.target_x = beam_x;
        in_run_config.target_y = beam_y;
    }
    // otherwise, shift the target position based on this measured beam position
    else {
        in_run_config.target_x += beam_x;
        in_run_config.target_y += beam_y;
    }

    // 2. Extract Detector Z position from HyCal/GEMs Z alignment parameter
    in_run_config.hycal_z = alignment_params["HC_z"];
    in_run_config.gem_z[0] = alignment_params["moller_GEM0_z"];
    in_run_config.gem_z[1] = alignment_params["moller_GEM1_z"];
    in_run_config.gem_z[2] = alignment_params["moller_GEM2_z"];
    in_run_config.gem_z[3] = alignment_params["moller_GEM3_z"];

    // 3. Extract GEMs X/Y alignment parameters(alignment to HyCal center)
    // should firstly complete the beam position extraction
    if(beam_x < 0.1 && beam_y < 0.1) {
        for (int det = 0; det < 4; ++det) {
            in_run_config.gem_x[det] -= 0.5 * alignment_params["moller_GEM" + std::to_string(det) + "_x"];
            in_run_config.gem_y[det] -= 0.5 * alignment_params["moller_GEM" + std::to_string(det) + "_y"];
        }
    }
    // 3.1 make other 3 GEMs alignment to GEM3 (the most upstream GEM)
    // 3.1.1 Align GEM2 to GEM3
    in_run_config.gem_x[2] -= 0.5 * alignment_params["Upstream_GEM_dx_d2-d3"];
    in_run_config.gem_y[2] -= 0.5 * alignment_params["Upstream_GEM_dy_d2-d3"];
    // 3.1.2 Align GEM0 to GEM1
    //in_run_config.gem_x[0] -= 0.5 * alignment_params["Downstream_GEM_dx_d0-d1"];
    //in_run_config.gem_y[0] -= 0.5 * alignment_params["Downstream_GEM_dy_d0-d1"];
    // 3.1.3 Align GEM1 to GEM3
    in_run_config.gem_x[1] -= 0.5 * alignment_params["GEM_layer_right_dx_d1-d3"];
    in_run_config.gem_y[1] -= 0.5 * alignment_params["GEM_layer_right_dy_d1-d3"];
    // 3.1.4 Align GEM0 to GEM2
    in_run_config.gem_x[0] -= 0.5 * alignment_params["GEM_layer_left_dx_d0-d2"];
    in_run_config.gem_y[0] -= 0.5 * alignment_params["GEM_layer_left_dy_d0-d2"];

    // 4. Extract GEMs phi alignment(roll rotate around Z axis) parameters (alignment to HyCal coordinate system)
    // TODO: don't sure which direction is positive for the phi rotation
    for (int det = 0; det < 4; ++det) {
        in_run_config.gem_tilt_z[det] -= alignment_params["phi_diff_" + std::to_string(det)];
    }

    //5. TODO: Extract GEMs rotation around X/Y axes (tilt) parameters

    // output the resolved new run_config file
    // Automatically generate output filename with iteration number
    fs::path base(config_base_path);
    std::string stem = base.stem().string();
    std::string ext = base.extension().string();
    if (ext.empty()) ext = ".json";  // Default to .json if no extension
    std::string actual_output_path = (base.parent_path() / 
                                     (stem + "_iter" + std::to_string(iter) + ext)).string();
    std::cout << "Iteration " << iter << " output: " << actual_output_path << "\n";
    WriteRunConfig(actual_output_path, run_num, in_run_config);
    
    // output results and summary
    const std::string summary_path = (fs::path(output_dir) / 
                                     ("alignment_summary_iter" + std::to_string(iter) + ".txt")).string();
    std::ofstream summary(summary_path);
    summary << "Detector Alignment Parameters\n";
    summary << std::string(50, '=') << "\n";
    for (const auto &[key, val] : alignment_params) {
        summary << std::setw(20) << key << ": " << std::fixed << std::setprecision(3) << val << " mm\n";
    }
    summary.close();
    std::cout << "Alignment summary saved to " << summary_path << "\n";

    // Build convergence graphs by reading all iteration summary files
    std::vector<TGraph*> convergence_graphs;
    TGraph *g_hc_x = nullptr;
    TGraph *g_hc_y = nullptr;
    TGraph *g_hc_z = nullptr;
    TGraph *g_gem[4][4];  // [chamber][param: 0=x, 1=y, 2=z, 3=phi]
    std::memset(g_gem, 0, sizeof(g_gem));
    {
        auto parse_summary = [](const std::string &path) -> std::map<std::string, double> {
            std::map<std::string, double> params;
            std::ifstream f(path);
            std::string line;
            while (std::getline(f, line)) {
                auto colon = line.rfind(':');
                if (colon == std::string::npos) continue;
                std::string key = line.substr(0, colon);
                key.erase(0, key.find_first_not_of(" \t"));
                key.erase(key.find_last_not_of(" \t") + 1);
                std::string val_str = line.substr(colon + 1);
                auto mm_pos = val_str.find("mm");
                if (mm_pos != std::string::npos) val_str = val_str.substr(0, mm_pos);
                try { params[key] = std::stod(val_str); } catch (...) {}
            }
            return params;
        };

        std::vector<std::map<std::string, double>> iter_params(iter + 1);
        for (int it = 1; it <= iter; ++it) {
            std::string p = (fs::path(output_dir) /
                            ("alignment_summary_iter" + std::to_string(it) + ".txt")).string();
            iter_params[it] = parse_summary(p);
        }

        // zdiff=true: plot change between consecutive iters (for absolute Z positions)
        auto make_graph = [&](const std::string &key, bool zdiff) -> TGraph* {
            std::vector<double> xs, ys;
            for (int it = 1; it <= iter; ++it) {
                if (!iter_params[it].count(key)) continue;
                double val = iter_params[it][key];
                if (zdiff) {
                    double prev = (it > 1 && iter_params[it-1].count(key))
                                  ? iter_params[it-1][key] : val;
                    val = val - prev;
                }
                xs.push_back(it);
                ys.push_back(val);
            }
            if (xs.empty()) return nullptr;
            auto *g = new TGraph((int)xs.size(), xs.data(), ys.data());
            g->SetName(("g_conv_" + key).c_str());
            std::string ytitle = zdiff ? "Delta (mm)" : "Value (mm)";
            g->SetTitle((key + " convergence;Iteration;" + ytitle).c_str());
            g->SetMarkerStyle(20);
            g->SetMarkerSize(1.2);
            return g;
        };

        // HyCal: x and y are direct adjustment values; z is iteration-to-iteration delta
        g_hc_x = make_graph("HC_x", false);
        g_hc_y = make_graph("HC_y", false);
        g_hc_z = make_graph("HC_z", true);
        
        if (g_hc_x) convergence_graphs.push_back(g_hc_x);
        if (g_hc_y) convergence_graphs.push_back(g_hc_y);
        if (g_hc_z) convergence_graphs.push_back(g_hc_z);
        
        // GEM x/y: actual applied correction derived from consecutive run_config files
        auto build_gem_pos_graph = [&](int det, int coord) -> TGraph* {
            std::string pfx = "GEM" + std::to_string(det);
            std::vector<double> xs, ys;
            for (int it = 1; it <= iter; ++it) {
                fs::path bp(config_base_path);
                std::string st = bp.stem().string(), ex = bp.extension().string();
                if (ex.empty()) ex = ".json";
                std::string p_curr = (bp.parent_path() / (st + "_iter" + std::to_string(it) + ex)).string();
                std::string p_prev;
                if (it == 1)
                    p_prev = (iter == 1) ? run_config_in : "";  // iter0 only known when running iter1
                else
                    p_prev = (bp.parent_path() / (st + "_iter" + std::to_string(it - 1) + ex)).string();
                if (p_prev.empty() || !fs::exists(p_prev) || !fs::exists(p_curr)) continue;
                // suppress LoadRunConfig's per-call stderr log
                std::streambuf *cerr_buf = std::cerr.rdbuf(nullptr);
                RunConfig cp = LoadRunConfig(p_prev, run_num);
                RunConfig cc = LoadRunConfig(p_curr, run_num);
                std::cerr.rdbuf(cerr_buf);
                double prev_abs = (coord == 0) ? (cp.gem_x[det] + cp.target_x) : (cp.gem_y[det] + cp.target_y);
                double curr_abs = (coord == 0) ? (cc.gem_x[det] + cc.target_x) : (cc.gem_y[det] + cc.target_y);
                xs.push_back(it);
                ys.push_back(prev_abs - curr_abs);
            }
            if (xs.empty()) return nullptr;
            std::string axis = (coord == 0) ? "X" : "Y";
            auto *g = new TGraph((int)xs.size(), xs.data(), ys.data());
            g->SetName(("g_conv_" + pfx + "_d" + (coord == 0 ? "x" : "y")).c_str());
            g->SetTitle((pfx + " " + axis + " correction;Iteration;Delta (mm)").c_str());
            g->SetMarkerStyle(20);
            g->SetMarkerSize(1.2);
            return g;
        };

        // GEM 0-3: x, y from run_config position deltas; z diff; phi direct
        for (int det = 0; det < 4; ++det) {
            std::string pfx = "GEM" + std::to_string(det);
            g_gem[det][0] = build_gem_pos_graph(det, 0);
            g_gem[det][1] = build_gem_pos_graph(det, 1);
            g_gem[det][2] = make_graph(pfx + "_z", true);
            g_gem[det][3] = make_graph("phi_diff_" + std::to_string(det), false);

            if (g_gem[det][0]) convergence_graphs.push_back(g_gem[det][0]);
            if (g_gem[det][1]) convergence_graphs.push_back(g_gem[det][1]);
            if (g_gem[det][2]) convergence_graphs.push_back(g_gem[det][2]);
            if (g_gem[det][3]) convergence_graphs.push_back(g_gem[det][3]);
        }
    }

    // write histograms to output file
    const std::string output_root = (fs::path(output_dir) / 
                                    ("alignment_histograms_iter" + std::to_string(iter) + ".root")).string();
    std::unique_ptr<TFile> out_file(TFile::Open(output_root.c_str(), "RECREATE"));
    if (out_file && !out_file->IsZombie()) {
        out_file->cd();
        out_file->mkdir("E_vs_Angle");
        out_file->cd("E_vs_Angle");
        h2_hycal_energy_vs_angle->Write();
        h2_gem_energy_vs_angle[0]->Write();
        h2_gem_energy_vs_angle[1]->Write();
        h2_gem_energy_vs_angle[2]->Write();
        h2_gem_energy_vs_angle[3]->Write();
        out_file->cd();
        out_file->mkdir("HyCal");
        out_file->cd("HyCal");
        h2_hycal_hits->Write();
        h1_hycal_CenterX->Write();
        h1_hycal_CenterY->Write();
        h1_hycal_Zdistance->Write();
        out_file->cd();
        out_file->mkdir("GEM");
        out_file->cd("GEM");
        for (int i = 0; i < 4; ++i) {
            h2_gem_hits[i]->Write();
            h1_gem_CenterX[i]->Write();
            h1_gem_CenterY[i]->Write();
            h1_gem_Zdistance[i]->Write();
            h1_phi_diff_hycal_gem[i]->Write();
        }
        out_file->cd();
        out_file->mkdir("Delta_GEMs");
        out_file->cd("Delta_GEMs");
        h1_deltaX_gem_up->Write();
        h1_deltaY_gem_up->Write();
        h1_deltaX_gem_down->Write();
        h1_deltaY_gem_down->Write();
        h1_deltaX_gem_layer_left->Write();
        h1_deltaY_gem_layer_left->Write();
        h1_deltaX_gem_layer_right->Write();
        h1_deltaY_gem_layer_right->Write();
        h2_deltaX_vs_deltaY_gem_left->Write();
        h2_deltaX_vs_deltaY_gem_right->Write();
        out_file->mkdir("Delta_GEM_vs_HyCal");
        out_file->cd("Delta_GEM_vs_HyCal");
        for (int i = 0; i < 4; ++i) {
            h1_deltaX_gem_hycal[i]->Write();
            h1_deltaY_gem_hycal[i]->Write();
        }
        out_file->cd();
        out_file->mkdir("Convergence");
        out_file->cd("Convergence");
        
        // Create a canvas for HyCal convergence plots
        if (g_hc_x && g_hc_y && g_hc_z) {
            auto *canvas_hc = new TCanvas("HyCal_Convergence", "HyCal Convergence Parameters", 400, 1200);
            canvas_hc->Divide(1, 3);
            
            // Plot HC_x
            canvas_hc->cd(1);
            g_hc_x->Draw("APL");
            g_hc_x->GetXaxis()->SetTitle("Iteration");
            g_hc_x->GetYaxis()->SetTitle("Value (mm)");
            
            // Plot HC_y
            canvas_hc->cd(2);
            g_hc_y->Draw("APL");
            g_hc_y->GetXaxis()->SetTitle("Iteration");
            g_hc_y->GetYaxis()->SetTitle("Value (mm)");
            
            // Plot HC_z
            canvas_hc->cd(3);
            g_hc_z->Draw("APL");
            g_hc_z->GetXaxis()->SetTitle("Iteration");
            g_hc_z->GetYaxis()->SetTitle("Delta (mm)");
            
            canvas_hc->Write();
        }
        
        // Create canvases for GEM chambers (2x2 layout for each: x, y, z, phi)
        for (int det = 0; det < 4; ++det) {
            if (g_gem[det][0] && g_gem[det][1] && g_gem[det][2] && g_gem[det][3]) {
                auto *canvas_gem = new TCanvas(
                    Form("GEM%d_Convergence", det), 
                    Form("GEM%d Convergence Parameters", det), 
                    400, 1600
                );
                canvas_gem->Divide(1, 4);
                
                // Plot GEM_x
                canvas_gem->cd(1);
                g_gem[det][0]->Draw("APL");
                g_gem[det][0]->GetXaxis()->SetTitle("Iteration");
                g_gem[det][0]->GetYaxis()->SetTitle("Delta (mm)");
                
                // Plot GEM_y
                canvas_gem->cd(2);
                g_gem[det][1]->Draw("APL");
                g_gem[det][1]->GetXaxis()->SetTitle("Iteration");
                g_gem[det][1]->GetYaxis()->SetTitle("Delta (mm)");
                
                // Plot GEM_z
                canvas_gem->cd(3);
                g_gem[det][2]->Draw("APL");
                g_gem[det][2]->GetXaxis()->SetTitle("Iteration");
                g_gem[det][2]->GetYaxis()->SetTitle("Delta (mm)");
                
                // Plot phi_diff
                canvas_gem->cd(4);
                g_gem[det][3]->Draw("APL");
                g_gem[det][3]->GetXaxis()->SetTitle("Iteration");
                g_gem[det][3]->GetYaxis()->SetTitle("Value (deg)");
                
                canvas_gem->Write();
            }
        }
        
        for (auto *g : convergence_graphs)
            if (g) g->Write();
        out_file->Close();
        std::cout << "Histograms saved to " << output_root << "\n";
    }
}

// ── Helper: Extract peak center from histogram via FWHM method ──────────
static double extract_peak(TH1F *hist)
{
    if (!hist || hist->GetEntries() <= 0) return 0.0;
    
    const int max_bin = hist->GetMaximumBin();
    const double peak = hist->GetBinContent(max_bin);
    if (peak <= 0.0) return 0.0;
    
    const double half_peak = 0.7 * peak;
    const int nbins = hist->GetNbinsX();
    
    // Find left and right edges at 70% of peak
    int left_bin = max_bin, right_bin = max_bin;
    while (left_bin > 1 && hist->GetBinContent(left_bin) >= half_peak) --left_bin;
    while (right_bin < nbins && hist->GetBinContent(right_bin) >= half_peak) ++right_bin;
    
    // Linear interpolation for precise FWHM edges
    auto crossing = [hist](int bin0, int bin1) {
        if (std::abs(hist->GetBinContent(bin1) - hist->GetBinContent(bin0)) < 1e-12)
            return hist->GetXaxis()->GetBinCenter(bin0);
        double frac = (0.7 * hist->GetBinContent(hist->GetMaximumBin()) - hist->GetBinContent(bin0)) /
                      (hist->GetBinContent(bin1) - hist->GetBinContent(bin0));
        frac = std::max(0.0, std::min(1.0, frac));
        return hist->GetXaxis()->GetBinCenter(bin0) + frac * hist->GetXaxis()->GetBinWidth(bin0);
    };
    
    double fit_lo = crossing(left_bin, left_bin + 1);
    double fit_hi = crossing(right_bin - 1, right_bin);
    double sigma = (fit_hi - fit_lo) / 2.354820045;
    
    // Gaussian fit for refined peak position
    TF1 fit_func("gaus_fit", "gaus", fit_lo, fit_hi);
    fit_func.SetParameters(peak, hist->GetXaxis()->GetBinCenter(max_bin), std::max(sigma, hist->GetXaxis()->GetBinWidth(max_bin)));
    
    int fit_status = hist->Fit(&fit_func, "QR");  // Draw fit curve on histogram
    if (fit_status == 0 && std::isfinite(fit_func.GetParameter(1)))
        return fit_func.GetParameter(1);  // return fitted center
    
    return hist->GetXaxis()->GetBinCenter(max_bin);  // fallback to bin center
}

static bool ProcessFile(const std::string &input_root, 
                        const std::string &db_dir,
                        const RunConfig &gRunConfig, const RunConfig &in_run_config,
                        std::vector<EventWithMoller> &all_moller_events, std::vector<EventWithMott> &all_mott_events,
                        bool show_progress, int max_events)
{
    // Implement the file processing logic here

    std::unique_ptr<TFile> f(TFile::Open(input_root.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << input_root << "\n";
        return false;
    }
    TTree *tree = dynamic_cast<TTree *>(f->Get("recon"));
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in " << input_root << "\n";
        return false;
    }
    
    if (show_progress) {
        std::cout << "Reading file: " << fs::path(input_root).filename() << "\n";
    }

    DetectorTransform                 hycal_transform;
    std::array<DetectorTransform, 4>  gem_transforms;

    hycal_transform.set(
        gRunConfig.hycal_x, gRunConfig.hycal_y, gRunConfig.hycal_z,
        gRunConfig.hycal_tilt_x, gRunConfig.hycal_tilt_y, gRunConfig.hycal_tilt_z);
    for (int d = 0; d < 4; ++d) {
        gem_transforms[d].set(
            gRunConfig.gem_x[d], gRunConfig.gem_y[d], gRunConfig.gem_z[d],
            gRunConfig.gem_tilt_x[d], gRunConfig.gem_tilt_y[d], gRunConfig.gem_tilt_z[d]);
    }

    DetectorTransform in_hycal_transform;
    std::array<DetectorTransform, 4>  in_gem_transforms;

    in_hycal_transform.set(
        in_run_config.hycal_x, in_run_config.hycal_y, in_run_config.hycal_z,
        in_run_config.hycal_tilt_x, in_run_config.hycal_tilt_y, in_run_config.hycal_tilt_z);
    for (int d = 0; d < 4; ++d) {
        in_gem_transforms[d].set(
            in_run_config.gem_x[d], in_run_config.gem_y[d], in_run_config.gem_z[d],
            in_run_config.gem_tilt_x[d], in_run_config.gem_tilt_y[d], in_run_config.gem_tilt_z[d]);
    }

    MatchingTools matching(2);
    matching.SetMatchRange(in_run_config.matching_radius);
    matching.SetSquareSelection(in_run_config.matching_use_square);
    matching.SetEnergyDependent(in_run_config.matching_energy_dependent);
    matching.SetMatchSigma(in_run_config.matching_sigma);

    EventVars_Recon ev;
    prad2::SetReconReadBranches(tree, ev);
    prad2::ReconMatchVectorBindings match_bindings;
    prad2::BindReconMatchVectorBranches(tree, ev, match_bindings);
    Long64_t n = tree->GetEntries();

    for (Long64_t i = 0; i < n; ++i) {
        tree->GetEntry(i);

        if (max_events > 0 && i >= max_events) break;
        
        // Progress display
        if (show_progress && (i + 1) % 1000 == 0) {
            int percent = static_cast<int>(100.0 * (i + 1) / n);
            std::cout << "  [" << std::setw(3) << percent << "%] Event " 
                      << std::setw(8) << (i + 1) << " / " << n << "\r";
            std::cout.flush();
        }
        
        // Process the event here

        // trigger selection
        bool is_sum      = (ev.trigger_bits & prad2::TBIT_sum) != 0;
        if (!is_sum) continue;

        // select events and analyze
        if (ev.n_clusters != 2 && ev.n_clusters != 1) continue;
        if (ev.cl_nblocks[0] < 3 || (ev.n_clusters == 2 && ev.cl_nblocks[1] < 3)) continue;

        // store all HyCal hits and GEM hits
        // transform hits from lab to detector coordinates
        // then transform them to the detector coordinates for the input run_config
        std::vector<HCHit> hycal_hits;
        std::vector<GEMHit> gem_hits[4];
        for (int i = 0; i < ev.n_clusters; ++i) {
            HCHit hit{
                ev.cl_x[i], ev.cl_y[i], ev.cl_z[i],
                ev.cl_energy[i], ev.cl_center[i], ev.cl_flag[i]
            };
            ApplyToLocal(hycal_transform, hit);
            ApplyToLab(in_hycal_transform, hit);
            hycal_hits.push_back(hit);
            ev.cl_x[i] = hit.x;
            ev.cl_y[i] = hit.y;
            ev.cl_z[i] = hit.z;
        }
        for (int d = 0; d < 4; ++d) {
            for (int j = 0; j < ev.n_gem_hits; ++j) {
                if (ev.det_id[j] != d) continue;
                GEMHit hit{
                    ev.gem_x[j], ev.gem_y[j], ev.gem_z[j], ev.det_id[j]
                };
                ApplyToLocal(gem_transforms[d], hit);
                ApplyToLab(in_gem_transforms[d], hit);
                gem_hits[d].push_back(hit);
            }
        }

        // do the matching between HyCal hits and GEM hits
        std::vector<MatchHit> matched_hits = matching.Match(hycal_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]);
        std::vector<MatchHit_perChamber> matched_hits_chamber = matching.MatchPerChamber(hycal_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]); 
        ev.clear_match_lists();
        for(int i = 0; i < matched_hits_chamber.size(); i++){
            auto &m = matched_hits_chamber[i];
            int cl_idx = m.hycal_idx;
            if( cl_idx != i) std::cerr << "Warning: cluster index mismatch in matched_hits_chamber: " << cl_idx << " vs " << i << "\n";
            for(int j = 0; j < 4; j++){
                for (const auto &gh : m.gem_hits[j]) {
                    ev.add_match(i, j, gh.x, gh.y, gh.z);
                }
            }
            ev.matchFlag[i] = m.mflag;
        }

        ev.matchNum = std::min((int)matched_hits.size(), prad2::kMaxClusters);
        for (int i = 0; i < ev.matchNum; i++){
            // save the matched GEM hit (must 2 matchings) info in mHit_ arrays for quick check
            ev.mHit_E[i] = matched_hits[i].hycal_hit.energy;
            ev.mHit_x[i] = matched_hits[i].hycal_hit.x;
            ev.mHit_y[i] = matched_hits[i].hycal_hit.y;
            ev.mHit_z[i] = matched_hits[i].hycal_hit.z;
            for(int j = 0; j < 2; j++) {
                ev.mHit_gx[i][j] =  matched_hits[i].gem[j].x;
                ev.mHit_gy[i][j] =  matched_hits[i].gem[j].y;
                ev.mHit_gz[i][j] =  matched_hits[i].gem[j].z;
                ev.mHit_gid[i][j] = matched_hits[i].gem[j].det_id; // placeholder for GEM hit ID if needed
            }
            ev.mHit_cl_index[i] = matched_hits[i].hycal_idx;
        }
        if (ev.matchNum == 1 && ev.n_clusters == 1) {
            // select single matched events for Mott electron
            float theta = std::atan2(std::sqrt(ev.cl_y[0]*ev.cl_y[0] + ev.cl_x[0]*ev.cl_x[0]), ev.cl_z[0]) * 180.0 / M_PI;
            float E = ev.cl_energy[0];
            float expectE = PhysicsTools::ExpectedEnergy(theta, gRunConfig.Ebeam, "ep");
            float sigma = 0.03 * std::sqrt(expectE * 1000.0);
            if (std::abs(E - expectE) < 3.0 * sigma && theta > 0.85) {
                // mark this event as a good Mott electron candidate
                EventWithMott thisEvent;
                thisEvent.HC_x = ev.cl_x[0];
                thisEvent.HC_y = ev.cl_y[0];
                thisEvent.HC_z = ev.cl_z[0];
                // count the number of matched GEM hits for this Mott electron candidate
                int count[4] = {0, 0, 0, 0};
                for(auto did : ev.match_det_id) {
                    count[did]++;
                }
                for(int did = 0; did < 4; did ++){
                    if ((ev.matchFlag[0] & (1u << did)) != 0 && count[did] == 1) {
                        float x, y, z;
                        if (!ev.first_match(0, did, x, y, z)) continue;
                        thisEvent.GEM_x[did] = x;
                        thisEvent.GEM_y[did] = y;
                        thisEvent.GEM_z[did] = z;
                        thisEvent.match[did] = true;
                    }
                }
                all_mott_events.push_back(thisEvent);
            }
        }

        if (ev.matchNum == 2 && ev.n_clusters == 2) {
            // use HyCal to judge and select Moller events should be good enough
            float theta1 = std::atan2(std::sqrt(ev.cl_y[0]*ev.cl_y[0] + ev.cl_x[0]*ev.cl_x[0]), ev.cl_z[0]) * 180.0 / M_PI;
            float theta2 = std::atan2(std::sqrt(ev.cl_y[1]*ev.cl_y[1] + ev.cl_x[1]*ev.cl_x[1]), ev.cl_z[1]) * 180.0 / M_PI;
            float phi1 = std::atan2(ev.cl_y[0], ev.cl_x[0]) * 180.0 / M_PI;
            float phi2 = std::atan2(ev.cl_y[1], ev.cl_x[1]) * 180.0 / M_PI;
            float E1 = ev.cl_energy[0];
            float E2 = ev.cl_energy[1];
            float expectE1 = PhysicsTools::ExpectedEnergy(theta1, gRunConfig.Ebeam, "ee");
            float expectE2 = PhysicsTools::ExpectedEnergy(theta2, gRunConfig.Ebeam, "ee");
            float sigma1 = 0.03 * std::sqrt(expectE1 * 1000.0);
            float sigma2 = 0.03 * std::sqrt(expectE2 * 1000.0);
            float sigma_sum = std::sqrt(sigma1*sigma1 + sigma2*sigma2);

            if (theta1 < 0.65 || theta2 < 0.65) continue;
            if (std::abs(phi1 - phi2) -180.0 > 8.0) continue;
            if (std::abs(E1 + E2 - gRunConfig.Ebeam) > 3.0 * sigma_sum) continue;
            if (std::abs(E1 - expectE1) > 3.0 * sigma1) continue;
            if (std::abs(E2 - expectE2) > 3.0 * sigma2) continue;

            EventWithMoller thisEvent;
            MollerEvent m_hycal, m_gemUp, m_gemDown;
            std::array<MollerEvent, 4> m_gem;
            // Initialize all m_gem entries to zero-energy Moller events
            for (int d = 0; d < 4; ++d) {
                m_gem[d] = MollerEvent(DataPoint(0, 0, 0, 0), DataPoint(0, 0, 0, 0));
            }

            m_hycal = MollerEvent(
                DataPoint(ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]),
                DataPoint(ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1]));
            /*
            for(int did = 0; did <= 1; did ++){
                if (ev.mHit_gid[0][0] == did
                    && ev.mHit_gid[1][0] == did) {
                    float x0, y0, z0, x1, y1, z1;
                    x0 = ev.mHit_gx[0][0];
                    y0 = ev.mHit_gy[0][0];
                    z0 = ev.mHit_gz[0][0];
                    x1 = ev.mHit_gx[1][0];
                    y1 = ev.mHit_gy[1][0];
                    z1 = ev.mHit_gz[1][0];
                    m_gem[did] = MollerEvent(
                        DataPoint(x0, y0, z0, ev.mHit_E[0]),
                        DataPoint(x1, y1, z1, ev.mHit_E[1]));
                }
            }
            for(int did = 2; did <= 3; did ++){
                if (ev.mHit_gid[0][1] == did
                    && ev.mHit_gid[1][1] == did) {
                    float x0, y0, z0, x1, y1, z1;
                    x0 = ev.mHit_gx[0][1];
                    y0 = ev.mHit_gy[0][1];
                    z0 = ev.mHit_gz[0][1];
                    x1 = ev.mHit_gx[1][1];
                    y1 = ev.mHit_gy[1][1];
                    z1 = ev.mHit_gz[1][1];
                    m_gem[did] = MollerEvent(
                        DataPoint(x0, y0, z0, ev.mHit_E[0]),
                        DataPoint(x1, y1, z1, ev.mHit_E[1]));
                }
            }*/
            for(int did = 0; did < 4; did ++){
                if (((ev.matchFlag[0] & (1u << did)) != 0)
                    && ((ev.matchFlag[1] & (1u << did)) != 0)) {
                    float x0, y0, z0, x1, y1, z1;
                    if (!ev.first_match(0, did, x0, y0, z0)) continue;
                    if (!ev.first_match(1, did, x1, y1, z1)) continue;
                    m_gem[did] = MollerEvent(
                        DataPoint(x0, y0, z0, ev.cl_energy[0]),
                        DataPoint(x1, y1, z1, ev.cl_energy[1]));
                }
            }
            m_gemUp = MollerEvent(
                DataPoint(ev.mHit_gx[0][1], ev.mHit_gy[0][1], ev.mHit_gz[0][1], ev.mHit_E[0]),
                DataPoint(ev.mHit_gx[1][1], ev.mHit_gy[1][1], ev.mHit_gz[1][1], ev.mHit_E[1]));
            m_gemDown = MollerEvent(
                DataPoint(ev.mHit_gx[0][0], ev.mHit_gy[0][0], ev.mHit_gz[0][0], ev.mHit_E[0]),
                DataPoint(ev.mHit_gx[1][0], ev.mHit_gy[1][0], ev.mHit_gz[1][0], ev.mHit_E[1]));
            
            thisEvent.event_num = ev.event_num;
            thisEvent.HC_moller = m_hycal;
            for(int did = 0; did < 4; did++){
                thisEvent.GEM_moller[did] = m_gem[did];
            }
            thisEvent.GEMup_moller = m_gemUp;
            thisEvent.GEMdown_moller = m_gemDown;
            all_moller_events.push_back(thisEvent);
        }
    }
    if (show_progress) {
        std::cout << "\n";  // newline after progress display
    }
    return true;
}

// ── Helpers ──────────────────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            std::string name = entry.path().filename().string();
            if (entry.is_regular_file() &&
                name.find("_recon") != std::string::npos &&
                name.size() >= 5 && name.compare(name.size() - 5, 5, ".root") == 0) {
                files.push_back(entry.path().string());
            }
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}


