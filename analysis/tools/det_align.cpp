
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

std::vector<EventWithMoller> AllMollerEvents;

// ── forward declarations ──────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path);

static double extract_peak(TH1F *hist);

static bool ProcessFile(const std::string &input_root,
                        const std::string &db_dir,
                        const RunConfig &gRunConfig, const RunConfig &in_run_config,
                        std::vector<EventWithMoller> &all_moller_events);


// ── Main ──────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");

    std::string  output_dir, run_config_in, run_config_out;
    int  max_files   = -1;
    int  num_threads = 4;

    int opt;
    while ((opt = getopt(argc, argv, "o:f:j:c:r:")) != -1) {
        switch (opt) {
            case 'o': output_dir       = optarg; break;
            case 'f': max_files        = std::atoi(optarg); break;
            case 'j': num_threads      = std::atoi(optarg); break;
            case 'c': run_config_in       = optarg; break;
            case 'r': run_config_out       = optarg; break;
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
            "       [-f max_files] [-j threads] \n"
            "       [-c run_config_in.json] [-r run_config_out.json]\n";
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

    RunConfig in_run_config = LoadRunConfig(run_config_in, run_num);

    // Each worker writes only to the result slot belonging to its replay_recon file.
    // This avoids locking the (potentially large) Moller vectors while events
    // are being reconstructed.  The main thread merges the slots after all
    // workers have finished.
    std::vector<std::vector<EventWithMoller>> events_per_file(num_files);
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
                                        events_per_file[idx]);
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
    for (int i = 0; i < num_files; ++i) {
        if (processed_ok[i]) {
            all_events.insert(all_events.end(),
                              events_per_file[i].begin(),
                              events_per_file[i].end());
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
    TH1F *h1_hycal_CenterX   = new TH1F("h1_hycal_mollerCenterX",   "HyCal Moller Center X",   400, -50, 50);
    TH1F *h1_hycal_CenterY   = new TH1F("h1_hycal_mollerCenterY",   "HyCal Moller Center Y",   400, -50, 50);
    TH1F *h1_hycal_Zdistance = new TH1F("h1_hycal_mollerZdistance", "HyCal Moller Z distance", 4000, 0, 10000);
    TH1F *h1_gem_CenterX[4], *h1_gem_CenterY[4], *h1_gem_Zdistance[4];
    for (int i = 0; i < 4; ++i) {
        h1_gem_CenterX[i]   = new TH1F(Form("h1_gem%d_mollerCenterX", i),   Form("GEM%d Moller Center X", i),   800, -50, 50);
        h1_gem_CenterY[i]   = new TH1F(Form("h1_gem%d_mollerCenterY", i),   Form("GEM%d Moller Center Y", i),   200, -50, 50);
        h1_gem_Zdistance[i] = new TH1F(Form("h1_gem%d_mollerZdistance", i), Form("GEM%d Moller Z distance", i), 4000, 0, 10000);
    }
    // Phi difference between HyCal and GEMs
    TH1F *h1_phi_diff_hycal_gem[4];
    for (int i = 0; i < 4; ++i)
        h1_phi_diff_hycal_gem[i] = new TH1F(Form("h1_phi_diff_hycal_gem%d", i), Form("Phi Difference HyCal-GEM%d", i), 200, -10, 10);
    // internal layer detector alignment histograms(2 chambers in the same layer, overlap region)
    TH1F *h1_deltaX_gem_up = new TH1F("h1_deltaX_gem_up", "Delta X GEM Up", 400, -20, 20);
    TH1F *h1_deltaX_gem_down = new TH1F("h1_deltaX_gem_down", "Delta X GEM Down", 400, -20, 20);
    TH1F *h1_deltaY_gem_up = new TH1F("h1_deltaY_gem_up", "Delta Y GEM Up", 400, -20, 20);
    TH1F *h1_deltaY_gem_down = new TH1F("h1_deltaY_gem_down", "Delta Y GEM Down", 400, -20, 20);
    TH1F *h1_deltaPhi_gem_up = new TH1F("h1_deltaPhi_gem_up", "Delta Phi GEM Up", 400, -20, 20);
    TH1F *h1_deltaPhi_gem_down = new TH1F("h1_deltaPhi_gem_down", "Delta Phi GEM Down", 200, -10, 10);
    // 2 layers GEM alignment histograms
    TH1F *h1_deltaX_gem_layer = new TH1F("h1_deltaX_gem_layer", "Delta X GEM Layer", 400, -20, 20);
    TH1F *h1_deltaY_gem_layer = new TH1F("h1_deltaY_gem_layer", "Delta Y GEM Layer", 400, -20, 20);
    TH1F *h1_deltaPhi_gem_layer = new TH1F("h1_deltaPhi_gem_layer", "Delta Phi GEM Layer", 200, -10, 10);
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

        // 2. per GEM chamber moller events
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

    // fit peak and resolve alignment parameters
    fs::create_directories(output_dir);
    gROOT->SetBatch(kTRUE);
    
    std::map<std::string, double> alignment_params;
    
    // Fit HyCal parameters
    alignment_params["HC_x"] = extract_peak(h1_hycal_CenterX);
    alignment_params["HC_y"] = extract_peak(h1_hycal_CenterY);
    alignment_params["HC_z"] = extract_peak(h1_hycal_Zdistance);
    
    // Fit GEM and phi parameters
    for (int det = 0; det < 4; ++det) {
        alignment_params["GEM" + std::to_string(det) + "_x"] = extract_peak(h1_gem_CenterX[det]);
        alignment_params["GEM" + std::to_string(det) + "_y"] = extract_peak(h1_gem_CenterY[det]);
        alignment_params["GEM" + std::to_string(det) + "_z"] = extract_peak(h1_gem_Zdistance[det]);
        alignment_params["phi_diff_" + std::to_string(det)] = extract_peak(h1_phi_diff_hycal_gem[det]);
    }
    
    // output results and summary
    const std::string summary_path = fs::path(output_dir) / "alignment_summary.txt";
    std::ofstream summary(summary_path);
    summary << "Detector Alignment Parameters\n";
    summary << std::string(50, '=') << "\n";
    for (const auto &[key, val] : alignment_params) {
        summary << std::setw(20) << key << ": " << std::fixed << std::setprecision(3) << val << " mm\n";
    }
    summary.close();
    std::cout << "Alignment summary saved to " << summary_path << "\n";
    
    // write histograms to output file
    const std::string output_root = fs::path(output_dir) / "alignment_histograms.root";
    std::unique_ptr<TFile> out_file(TFile::Open(output_root.c_str(), "RECREATE"));
    if (out_file && !out_file->IsZombie()) {
        h2_hycal_energy_vs_angle->Write();
        h2_gem_energy_vs_angle[0]->Write();
        h2_gem_energy_vs_angle[1]->Write();
        h2_gem_energy_vs_angle[2]->Write();
        h2_gem_energy_vs_angle[3]->Write();
        h2_hycal_hits->Write();
        h1_hycal_CenterX->Write();
        h1_hycal_CenterY->Write();
        h1_hycal_Zdistance->Write();
        for (int i = 0; i < 4; ++i) {
            h2_gem_hits[i]->Write();
            h1_gem_CenterX[i]->Write();
            h1_gem_CenterY[i]->Write();
            h1_gem_Zdistance[i]->Write();
            h1_phi_diff_hycal_gem[i]->Write();
        }
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
    
    const double half_peak = 0.8 * peak;
    const int nbins = hist->GetNbinsX();
    
    // Find left and right edges at 80% of peak
    int left_bin = max_bin, right_bin = max_bin;
    while (left_bin > 1 && hist->GetBinContent(left_bin) >= half_peak) --left_bin;
    while (right_bin < nbins && hist->GetBinContent(right_bin) >= half_peak) ++right_bin;
    
    // Linear interpolation for precise FWHM edges
    auto crossing = [hist](int bin0, int bin1) {
        if (std::abs(hist->GetBinContent(bin1) - hist->GetBinContent(bin0)) < 1e-12)
            return hist->GetXaxis()->GetBinCenter(bin0);
        double frac = (0.8 * hist->GetBinContent(hist->GetMaximumBin()) - hist->GetBinContent(bin0)) /
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
                        std::vector<EventWithMoller> &all_moller_events)
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
        // Process the event here

        // trigger selection
        bool is_sum      = (ev.trigger_bits & prad2::TBIT_sum) != 0;
        if (!is_sum) continue;

        // select events and analyze
        if (ev.n_clusters != 2) continue;
        if (ev.cl_nblocks[0] < 3 || ev.cl_nblocks[1] < 3) continue;

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

        if (ev.matchNum != 2) continue;

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

        if (theta1 < 0.8 || theta2 < 0.8) continue;
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
            DataPoint(ev.mHit_gx[0][1], ev.mHit_gy[0][1], ev.mHit_gz[0][1], ev.cl_energy[0]),
            DataPoint(ev.mHit_gx[1][1], ev.mHit_gy[1][1], ev.mHit_gz[1][1], ev.cl_energy[1]));
        m_gemDown = MollerEvent(
            DataPoint(ev.mHit_gx[0][0], ev.mHit_gy[0][0], ev.mHit_gz[0][0], ev.cl_energy[0]),
            DataPoint(ev.mHit_gx[1][0], ev.mHit_gy[1][0], ev.mHit_gz[1][0], ev.cl_energy[1]));
        
        thisEvent.event_num = ev.event_num;
        thisEvent.HC_moller = m_hycal;
        for(int did = 0; did < 4; did++){
            thisEvent.GEM_moller[did] = m_gem[did];
        }
        thisEvent.GEMup_moller = m_gemUp;
        thisEvent.GEMdown_moller = m_gemDown;
        all_moller_events.push_back(thisEvent);
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


