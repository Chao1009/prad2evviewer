// physics_calib.cpp: multi-threaded HyCal physics calibration
// Raw ROOT files are processed in rounds, with one file assigned to each
// worker thread. Each worker reconstructs HyCal clusters and accumulates a
// 5x5 energy spectrum for every PbWO4 module. The per-thread histograms are
// merged before the module fits and calibration update are performed.
//
// Events are selected from the sum trigger by requiring one reconstructed
// cluster with at least four blocks, a hit near the center of its seed crystal,
// and at least 60% of the cluster energy in that crystal. Dead seed modules
// from the run configuration are rejected. For transition modules, only hits
// on the inner side of the crystal are retained. Time-compatible PbWO4 hits in
// a 5x5 window are then summed into the seed module's energy spectrum.
//
// For each non-dead module, the reconstructed elastic e-p peak is fitted near
// the expected energy calculated from the run beam energy and detector
// geometry. The expected/fitted peak ratio is damped to 70% of the full
// correction, limited to [0.5, 2.0], and applied to the current calibration
// constant. Later iterations use the preceding result. Dead modules are not
// calibrated, and dead/dead-neighbor flags are recorded in the fit-result JSON.
//=============================================================================
//
// Usage: physics_calib <input_raw.root|dir> [more files/dirs...]
//                      [-i iteration] [-o output_dir]
//                      [-c seed_calib.json] [-j num_threads]
//   - input_raw.root|dir: input ROOT file or directory containing *_raw.root
//   - iteration: calibration iteration, starting from 1 (default: 1)
//   - output_dir: base output directory (default: current directory)
//   - seed_calib.json: input calibration for iteration 1
//                      (default: database calibration seed)
//   - num_threads: number of worker threads (default: 4)
//
// Outputs are written under output_dir/Physics_calib/run<run_number>/:
//   - calib_factor_iterN.json: updated HyCal calibration constants
//   - calib_result_iterN.json: per-module fit and correction results
//   - calib_result_iterN.root: merged spectra and diagnostic histograms
//=============================================================================

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "WaveAnalyzer.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "load_daq_config.h"
#include "RunInfoConfig.h"
#include "gain_factor.h"
#include "PipelineBuilder.h"

#include <TFile.h>
#include <TTree.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TROOT.h>
#include <TClass.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <filesystem>
#include <vector>
#include <thread>
#include <mutex>
#include <memory>
#include <algorithm>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;

using EventVars = prad2::RawEventData;
using namespace analysis;

// ── Per-thread accumulated results ──────────────────────────────────────────
struct HistResult {
    std::vector<std::unique_ptr<TH1F>>       h1_E_modules;   // indexed by module ID (index 0 is module W1)
    std::unique_ptr<TH2F>                    h2_energy_theta;
    std::unique_ptr<TH2F>                    hit_pos;
    std::unique_ptr<TH1F>                    h_E_1cl;
    std::unique_ptr<TH1F>                    h_center_energy_fraction;
    std::unique_ptr<TH1F>                    h_center_energy;
    std::unique_ptr<TH1F>                    h_2nd_energy_fraction;
    std::unique_ptr<TH1F>                    h_3rd_energy_fraction;
    std::unique_ptr<TH1F>                    h_fit_peak_energy;
    std::unique_ptr<TH1F>                    h_fit_peak_ratio;
    std::unique_ptr<TH1F>                    h_fit_peak_chi2ndf;
    std::unique_ptr<TH1F>                    h_fit_peak_sigma;
    long long                                events_processed = 0;
};

const int angle_bins = 50; const double angle_min = 0., angle_max = 5.;
const int energy_bins = 500; const double energy_min = 0., energy_max = 5000.;
const int pos_bins = 720; const double pos_min = -360., pos_max = 360.;
const int center_energy_fraction_bins = 100; const double center_energy_fraction_min = 0., center_energy_fraction_max = 1.;
const int fit_ratio_bins = 200; const double fit_ratio_min = 0., fit_ratio_max = 2.;
const int fit_chi2ndf_bins = 200; const double fit_chi2ndf_min = 0., fit_chi2ndf_max = 50.;
const int fit_sigma_bins = 200; const double fit_sigma_min = 0., fit_sigma_max = 200.;

bool ProcessRawFiles (const std::string &input_raw, RunConfig &gRunConfig, 
                      const std::string &db_dir, const std::string &recon_config_file,
                      const std::string &calib_file, HistResult *res);

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
    // ROOT multi-thread safety (must be called before any ROOT object creation)
    ROOT::EnableThreadSafety();
    // Force dictionary loading in main thread
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");
    TClass::GetClass("TH1F");
    TClass::GetClass("TH2F");

    // ── Argument parsing ─────────────────────────────────────────────────────
    std::string output_path, daq_config_file, seed_calib_file;
    int  iteration   = 1;
    int  max_events  = -1;
    int  num_threads = 4;

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR")) db_dir = env;

    std::string recon_config_file = db_dir + "/reconstruction_config.json";

    int opt;
    while ((opt = getopt(argc, argv, "i:o:j:c:")) != -1) {
        switch (opt) {
            case 'i': iteration        = std::atoi(optarg); break;
            case 'o': output_path = optarg; break;
            case 'c': seed_calib_file  = optarg; break;
            case 'j': num_threads      = std::atoi(optarg); break;
        }
    }

    // Collect all input files
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: calib_5by5 <input_raw.root|dir> [more...] "
                     "[-i iter] [-o output_dir] [-c seed_calib.json] [-j threads]\n";
        return 1;
    }

    // ── Run number / output paths ─────────────────────────────────────────────
    int run_num = get_run_int(root_files[0]);
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
    std::string run_str = "run" + std::to_string(run_num);
    if (output_path.empty()) output_path = ".";
    std::string run_out_dir = output_path + "/Physics_calib/" + run_str;
    fs::create_directories(run_out_dir);
    std::cerr << "Output directory: " << run_out_dir << "\n";

    std::string input_calib_file, output_calib_file, output_root_file, output_json_file;
    if (iteration == 1)
        input_calib_file = !seed_calib_file.empty()
            ? seed_calib_file
            : db_dir + "/calibration/calibration_factor_3p5_June7.json";
    else if (iteration > 1)
        input_calib_file = run_out_dir + Form("/calib_factor_iter%d.json", iteration - 1);
    else {
        std::cerr << "Invalid iteration number: " << iteration << ". Must be >= 1.\n";
        return 1;
    }
    if (iteration > 1 && !seed_calib_file.empty())
        std::cerr << "Warning: -c " << seed_calib_file
                  << " is ignored for iteration > 1; using " << input_calib_file << "\n";
    output_calib_file = run_out_dir + Form("/calib_factor_iter%d.json", iteration);

    output_root_file = run_out_dir + Form("/calib_result_iter%d.root", iteration);
    output_json_file = run_out_dir + Form("/calib_result_iter%d.json", iteration);

    // ── Thread count ─────────────────────────────────────────────────────────
    int n_files    = static_cast<int>(root_files.size());
    num_threads    = std::max(1, std::min(num_threads, n_files));
    int num_rounds = (n_files + num_threads - 1) / num_threads;
    std::cout << "Processing " << n_files << " file(s) with "
              << num_threads << " thread(s), " << num_rounds << " round(s)\n";

    // ── Initialize per-thread results (once, reused across rounds) ───────────
    std::vector<std::unique_ptr<HistResult>> results(num_threads);
    std::mutex io_mtx;

    for (int tid = 0; tid < num_threads; ++tid) {
        auto res = std::make_unique<HistResult>();

        fdec::HyCalSystem hycal;
        hycal.Init(db_dir + "/hycal_map.json");

        res->h1_E_modules.resize(1156);
        for (int i = 0; i < 1156; ++i) {
            const int mod_id = i + 1000 + 1; // module IDs start at 1001(W1)
            res->h1_E_modules[i] = std::make_unique<TH1F>(
                Form("h1_E_mod_%d_tid%d", mod_id, tid),
                Form("Module W%d cluster energy;E (MeV);Counts", mod_id-1000),
                energy_bins, energy_min, energy_max);
            res->h1_E_modules[i]->SetDirectory(nullptr);
        }

        res->h2_energy_theta = std::make_unique<TH2F>(
            Form("h2_energy_theta_tid%d", tid),
            "Energy vs Theta;Theta (deg);Energy (MeV)",
            angle_bins, angle_min, angle_max, energy_bins, energy_min, energy_max);
        res->h2_energy_theta->SetDirectory(nullptr);

        res->hit_pos = std::make_unique<TH2F>(
            Form("hit_pos_tid%d", tid),
            "One-cluster hit positions;hycal X (mm);hycal Y (mm)",
            pos_bins, pos_min, pos_max, pos_bins, pos_min, pos_max);
        res->hit_pos->SetDirectory(nullptr);

        res->h_E_1cl = std::make_unique<TH1F>(
            Form("h_E_1cl_tid%d", tid),
            "Single-cluster energy;E (MeV);Counts",
            energy_bins, energy_min, energy_max);
        res->h_E_1cl->SetDirectory(nullptr);

        res->h_center_energy_fraction = std::make_unique<TH1F>(
            Form("h_center_energy_fraction_tid%d", tid),
            "Center energy fraction;E_{center}/E_{cluster};Counts",
            center_energy_fraction_bins,
            center_energy_fraction_min, center_energy_fraction_max);
        res->h_center_energy_fraction->SetDirectory(nullptr);

        res->h_center_energy = std::make_unique<TH1F>(
            Form("h_center_energy_tid%d", tid),
            "Center module energy;E_{center} (MeV);Counts",
            energy_bins, energy_min, energy_max);
        res->h_center_energy->SetDirectory(nullptr);

        res->h_2nd_energy_fraction = std::make_unique<TH1F>(
            Form("h_2nd_energy_fraction_tid%d", tid),
            "Second 5x5 energy layer fraction;E_{layer 2}/E_{5x5};Counts",
            center_energy_fraction_bins,
            center_energy_fraction_min, center_energy_fraction_max);
        res->h_2nd_energy_fraction->SetDirectory(nullptr);

        res->h_3rd_energy_fraction = std::make_unique<TH1F>(
            Form("h_3rd_energy_fraction_tid%d", tid),
            "Third 5x5 energy layer fraction;E_{layer 3}/E_{5x5};Counts",
            center_energy_fraction_bins,
            center_energy_fraction_min, center_energy_fraction_max);
        res->h_3rd_energy_fraction->SetDirectory(nullptr);

        res->h_fit_peak_energy = std::make_unique<TH1F>(
            Form("h_fit_peak_energy_tid%d", tid),
            "Fitted peak energy;E_{peak} (MeV);Modules",
            energy_bins, energy_min, energy_max);
        res->h_fit_peak_energy->SetDirectory(nullptr);

        res->h_fit_peak_ratio = std::make_unique<TH1F>(
            Form("h_fit_peak_ratio_tid%d", tid),
            "Calibration ratio;E_{expected}/E_{peak};Modules",
            fit_ratio_bins, fit_ratio_min, fit_ratio_max);
        res->h_fit_peak_ratio->SetDirectory(nullptr);

        res->h_fit_peak_chi2ndf = std::make_unique<TH1F>(
            Form("h_fit_peak_chi2ndf_tid%d", tid),
            "Peak-fit #chi^{2}/NDF;#chi^{2}/NDF;Modules",
            fit_chi2ndf_bins, fit_chi2ndf_min, fit_chi2ndf_max);
        res->h_fit_peak_chi2ndf->SetDirectory(nullptr);

        res->h_fit_peak_sigma = std::make_unique<TH1F>(
            Form("h_fit_peak_sigma_tid%d", tid),
            "Fitted peak sigma;#sigma (MeV);Modules",
            fit_sigma_bins, fit_sigma_min, fit_sigma_max);
        res->h_fit_peak_sigma->SetDirectory(nullptr);

        results[tid] = std::move(res);
    }

    // ── Process files in rounds: num_threads files per round, 1 file/thread ──
    //  Each round does a local work pass, then the main thread can merge later.
    for (int round = 0; round < num_rounds; ++round) {
        int round_start        = round * num_threads;
        int round_end          = std::min(round_start + num_threads, n_files);
        int threads_this_round = round_end - round_start;

        std::cout << "\nRound " << (round + 1) << "/" << num_rounds
                  << ": files [" << round_start << ", " << round_end - 1 << "]\n";

        std::vector<std::thread> threads;
        threads.reserve(threads_this_round);

        for (int t = 0; t < threads_this_round; ++t) {
            threads.emplace_back([&, t, round]() {
                int fi = round * num_threads + t;
                auto *res = results[t].get();
                (void)res;

                bool ok = ProcessRawFiles(root_files[fi], gRunConfig,
                    db_dir, recon_config_file, input_calib_file, res);
                {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cout << "[thread " << t << "] file " << fi
                              << " / " << (n_files - 1) << ": "
                              << root_files[fi] << " -> "
                              << (ok ? "OK" : "FAILED") << "\n";
                }
            });
        }

        for (auto &th : threads) {
            if (th.joinable()) th.join();
        }
    }

    // ── Merge histograms (single-threaded) ────────────────────────────────────
    std::cout << "\nAll rounds finished. Merging histograms...\n";
    HistResult merged_result;
    // Initialize merged histograms
    merged_result.h1_E_modules.resize(1156);
    for (int i = 0; i < 1156; ++i) {
        const int mod_id = i + 1000 + 1; // module IDs start at 1001(W1)
        merged_result.h1_E_modules[i] = std::make_unique<TH1F>(
            Form("h1_E_mod_%d_merged", mod_id),
            Form("Module W%d cluster energy;E (MeV);Counts", mod_id-1000),
            energy_bins, energy_min, energy_max);
        merged_result.h1_E_modules[i]->SetDirectory(nullptr);
    }
    merged_result.h2_energy_theta = std::make_unique<TH2F>(
        "h2_energy_theta_merged",
        "Cluster energy vs theta;#theta (deg);E (MeV)",
        angle_bins, angle_min, angle_max,
        energy_bins, energy_min, energy_max);
    merged_result.h2_energy_theta->SetDirectory(nullptr);

    merged_result.hit_pos = std::make_unique<TH2F>(
        "hit_pos_merged",
        "Hit position;X (cm);Y (cm)",
        pos_bins, pos_min, pos_max,
        pos_bins, pos_min, pos_max);
    merged_result.hit_pos->SetDirectory(nullptr);

    merged_result.h_E_1cl = std::make_unique<TH1F>(
        "h_E_1cl_merged",
        "Single-cluster energy;E (MeV);Counts",
        energy_bins, energy_min, energy_max);
    merged_result.h_E_1cl->SetDirectory(nullptr);

    merged_result.h_center_energy_fraction = std::make_unique<TH1F>(
        "h_center_energy_fraction",
        "Center energy fraction;E_{center}/E_{cluster};Counts",
        center_energy_fraction_bins,
        center_energy_fraction_min, center_energy_fraction_max);
    merged_result.h_center_energy_fraction->SetDirectory(nullptr);

    merged_result.h_center_energy = std::make_unique<TH1F>(
        "h_center_energy",
        "Center module energy;E_{center} (MeV);Counts",
        energy_bins, energy_min, energy_max);
    merged_result.h_center_energy->SetDirectory(nullptr);

    merged_result.h_2nd_energy_fraction = std::make_unique<TH1F>(
        "h_2nd_energy_fraction",
        "Second 5x5 energy layer fraction;E_{layer 2}/E_{5x5};Counts",
        center_energy_fraction_bins,
        center_energy_fraction_min, center_energy_fraction_max);
    merged_result.h_2nd_energy_fraction->SetDirectory(nullptr);

    merged_result.h_3rd_energy_fraction = std::make_unique<TH1F>(
        "h_3rd_energy_fraction",
        "Third 5x5 energy layer fraction;E_{layer 3}/E_{5x5};Counts",
        center_energy_fraction_bins,
        center_energy_fraction_min, center_energy_fraction_max);
    merged_result.h_3rd_energy_fraction->SetDirectory(nullptr);

    merged_result.h_fit_peak_energy = std::make_unique<TH1F>(
        "h_fit_peak_energy",
        "Fitted peak energy;E_{peak} (MeV);Modules",
        energy_bins, energy_min, energy_max);
    merged_result.h_fit_peak_energy->SetDirectory(nullptr);

    merged_result.h_fit_peak_ratio = std::make_unique<TH1F>(
        "h_fit_peak_ratio",
        "Calibration ratio;E_{expected}/E_{peak};Modules",
        fit_ratio_bins, fit_ratio_min, fit_ratio_max);
    merged_result.h_fit_peak_ratio->SetDirectory(nullptr);

    merged_result.h_fit_peak_chi2ndf = std::make_unique<TH1F>(
        "h_fit_peak_chi2ndf",
        "Peak-fit #chi^{2}/NDF;#chi^{2}/NDF;Modules",
        fit_chi2ndf_bins, fit_chi2ndf_min, fit_chi2ndf_max);
    merged_result.h_fit_peak_chi2ndf->SetDirectory(nullptr);

    merged_result.h_fit_peak_sigma = std::make_unique<TH1F>(
        "h_fit_peak_sigma",
        "Fitted peak sigma;#sigma (MeV);Modules",
        fit_sigma_bins, fit_sigma_min, fit_sigma_max);
    merged_result.h_fit_peak_sigma->SetDirectory(nullptr);
    merged_result.events_processed = 0;

    for (int tid = 0; tid < num_threads; ++tid) {
        auto *res = results[tid].get();
        if (!res) continue;

        // Merge per-module energy histograms
        for (int i = 0; i < 1156; ++i) {
            if (res->h1_E_modules[i]) {
                TH1F *main_h = merged_result.h1_E_modules[i].get();
                main_h->Add(res->h1_E_modules[i].get());
            }
        }

        // Merge E-vs-theta 2D histogram
        if (res->h2_energy_theta) {
            TH2F *main_etheta = merged_result.h2_energy_theta.get();
            main_etheta->Add(res->h2_energy_theta.get());
        }

        // Merge hit position histogram
        if (res->hit_pos) {
            TH2F *main_hitpos = merged_result.hit_pos.get();
            main_hitpos->Add(res->hit_pos.get());
        }

        // Merge single-cluster energy histogram
        if (res->h_E_1cl) {
            TH1F *main_hE1cl = merged_result.h_E_1cl.get();
            main_hE1cl->Add(res->h_E_1cl.get());
        }
        if (res->h_center_energy_fraction) {
            merged_result.h_center_energy_fraction->Add(
                res->h_center_energy_fraction.get());
        }
        if (res->h_center_energy) {
            merged_result.h_center_energy->Add(res->h_center_energy.get());
        }
        if (res->h_2nd_energy_fraction) {
            merged_result.h_2nd_energy_fraction->Add(
                res->h_2nd_energy_fraction.get());
        }
        if (res->h_3rd_energy_fraction) {
            merged_result.h_3rd_energy_fraction->Add(
                res->h_3rd_energy_fraction.get());
        }
        if (res->h_fit_peak_energy) {
            merged_result.h_fit_peak_energy->Add(res->h_fit_peak_energy.get());
        }
        if (res->h_fit_peak_ratio) {
            merged_result.h_fit_peak_ratio->Add(res->h_fit_peak_ratio.get());
        }
        if (res->h_fit_peak_chi2ndf) {
            merged_result.h_fit_peak_chi2ndf->Add(res->h_fit_peak_chi2ndf.get());
        }
        if (res->h_fit_peak_sigma) {
            merged_result.h_fit_peak_sigma->Add(res->h_fit_peak_sigma.get());
        }
        merged_result.events_processed += res->events_processed;
    }

    // rsolve new calibration constants from the histograms of each module's energy distribution
    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");
    hycal.LoadCalibration(input_calib_file);
    prad2::ApplyHyCalDeadModules(gRunConfig.hycal_dead_modules, hycal);
    analysis::PhysicsTools physics(hycal);

    // save the new calibration results, later we can write them into a JSON file
    struct CalibrationResult {
        int module_id;
        float old_calib_factor;
        float new_calib_factor;
        float fit_ratio; // ratio of expected peak to fitted peak
        float fit_peak;
        float expected_peak;
        // fitting quality metrics
        float fit_sigma;
        float fit_chi2ndf;
        bool fit_good; // true if the fit is considered good
        bool is_dead; // true if the module is dead (read from RunConfig)
        bool is_deadNeighbor; // true if the module is in a 3 by 3 region of dead modules
    };
    std::vector<CalibrationResult> calib_results;

    int n_calibrated = 0;
    for (int i = 0; i < 1156; ++i) {
        if (!merged_result.h1_E_modules[i]) continue;
        TH1F *h = merged_result.h1_E_modules[i].get();
        if (h->GetEntries() < 100) continue; // skip modules with too few entries

        int mod_id = i + 1000 + 1; // module IDs start at 1001(W1)
        auto mod = hycal.module_by_id(mod_id);
        if (!mod) continue;

        bool is_dead = fdec::test_bit(mod->flag, fdec::kDeadModule);
        bool is_deadNeighbor = fdec::test_bit(mod->flag, fdec::kDeadNeighbor);
        if (is_dead) {
            calib_results.push_back({mod_id, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, false, is_dead, is_deadNeighbor});
            continue; // skip dead modules
        }

        float theta_deg = std::atan(std::sqrt(mod->x * mod->x + mod->y * mod->y)
                                    / gRunConfig.hycal_z) * 180.f / 3.14159265f;
        float expected_peak = analysis::PhysicsTools::ExpectedEnergy(theta_deg, gRunConfig.Ebeam, "ep");

        auto [peak, sigma, chi2] = physics.fitGaus(h, expected_peak);
        bool fit_good = (peak > 0 && sigma > 0 && sigma < 2. * 0.03*peak/std::sqrt(peak/1000.f) && chi2 < 10.f);
        if (!fit_good) {
            std::cout << "Check!!! Module W" << (mod_id - 1000)
                 << ": fit failed (peak=" << peak
                 << ", sigma=" << sigma
                 << ", chi2/ndf=" << chi2 << ")\n";
        }
        if (peak <= 0) peak = expected_peak; // fallback to expected if fit failed
        float ratio         = expected_peak / peak; 
        // apply a conservative factor to avoid over-correction
        ratio = (ratio - 1.f) * 0.7f + 1.f; // apply a conservative factor to avoid over-correction
        if(ratio < 0.5) ratio = 0.5;
        if(ratio > 2.0) ratio = 2.0;

        merged_result.h_fit_peak_energy->Fill(peak);
        merged_result.h_fit_peak_ratio->Fill(ratio);
        merged_result.h_fit_peak_chi2ndf->Fill(chi2);
        merged_result.h_fit_peak_sigma->Fill(sigma);

        double current_factor = hycal.GetCalibConstant(mod_id);
        double new_factor = current_factor * ratio;
        hycal.SetCalibConstant(mod_id, new_factor);
        hycal.SetCalibBaseEnergy(mod_id, expected_peak);
        n_calibrated++;

        calib_results.push_back({mod_id, 
            static_cast<float>(current_factor), 
            static_cast<float>(new_factor),
            ratio, static_cast<float>(peak), expected_peak,
            static_cast<float>(sigma), static_cast<float>(chi2), fit_good,
            is_dead, is_deadNeighbor});
    }

    // Write calibration results to a json file
    hycal.PrintCalibConstants(output_calib_file);
    std::ofstream json_out(output_json_file);
    if (json_out.is_open()) {
        json_out << "[\n";
        for (size_t i = 0; i < calib_results.size(); ++i) {
            const auto &res = calib_results[i];
            json_out << "  {"
                     << "\"module_id\": " << res.module_id << ", "
                     << "\"old_factor\": " << res.old_calib_factor << ", "
                     << "\"new_factor\": " << res.new_calib_factor << ", "
                     << "\"ratio\": " << res.fit_ratio << ", "
                     << "\"peak\": " << res.fit_peak << ", "
                     << "\"expected_peak\": " << res.expected_peak << ", "
                     << "\"sigma\": " << res.fit_sigma << ", "
                     << "\"chi2/ndf\": " << res.fit_chi2ndf << ", "
                     << "\"fit_good\": " << (res.fit_good ? "true" : "false") << ", "
                     << "\"is_dead\": " << (res.is_dead ? "true" : "false") << ", "
                     << "\"is_deadNeighbor\": " << (res.is_deadNeighbor ? "true" : "false")
                     << "}" << (i + 1 < calib_results.size() ? "," : "") << "\n";
        }
        json_out << "]\n";
        json_out.close();
    }

    // ── Save merged histograms to output ROOT file ───────────────────────────
    TFile *outfile = TFile::Open(output_root_file.c_str(), "RECREATE");
    for (int i = 0; i < 1156; ++i) {
        if (merged_result.h1_E_modules[i]) merged_result.h1_E_modules[i]->Write();
    }
    if (merged_result.h2_energy_theta) merged_result.h2_energy_theta->Write();
    if (merged_result.hit_pos) merged_result.hit_pos->Write();
    if (merged_result.h_E_1cl) merged_result.h_E_1cl->Write();
    if (merged_result.h_center_energy_fraction) merged_result.h_center_energy_fraction->Write();
    if (merged_result.h_center_energy) merged_result.h_center_energy->Write();
    if (merged_result.h_2nd_energy_fraction) merged_result.h_2nd_energy_fraction->Write();
    if (merged_result.h_3rd_energy_fraction) merged_result.h_3rd_energy_fraction->Write();
    if (merged_result.h_fit_peak_energy) merged_result.h_fit_peak_energy->Write();
    if (merged_result.h_fit_peak_ratio) merged_result.h_fit_peak_ratio->Write();
    if (merged_result.h_fit_peak_chi2ndf) merged_result.h_fit_peak_chi2ndf->Write();
    if (merged_result.h_fit_peak_sigma) merged_result.h_fit_peak_sigma->Write();

    outfile->Close();
    delete outfile;
}


bool ProcessRawFiles (const std::string &input_raw, RunConfig &gRunConfig, 
                      const std::string &db_dir, const std::string &recon_config_file,
                      const std::string &calib_file, HistResult *res)
{
    // Detectors: PRad-II flows through PipelineBuilder so the wiring stays in
    // one place (see prad2det/include/PipelineBuilder.h).
    fdec::HyCalSystem                 hycal;
    fdec::ClusterConfig               cluster_cfg;
    prad2::HyCalTimeCuts              hc_time_cuts;
    prad2::HyCalRfOffsets             hc_rf_offsets;

    int run_num = get_run_int(input_raw);

    prad2::Pipeline pipeline = prad2::PipelineBuilder()
        .set_hycal_calib(calib_file)
        .set_recon_config(recon_config_file)
        .set_database_dir(db_dir)
        .set_daq_config("")
        .set_hycal_map("") // empty falls back to defaults
        .set_gem_map("") // empty falls back to defaults
        .set_run_number_from_evio(input_raw)
        .set_log_stream(&std::cerr)
        .build();

    hycal            = std::move(pipeline.hycal);
    cluster_cfg      = pipeline.hycal_cluster_cfg;
    hc_time_cuts     = std::move(pipeline.hycal_time_cuts);
    hc_rf_offsets    = std::move(pipeline.hycal_rf_offsets);

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);

    // set up raw read branches for the input tree
    TFile *infile = TFile::Open(input_raw.c_str(), "READ");
    if (!infile || !infile->IsOpen()) {
        std::cerr << "Replay: cannot open " << input_raw << "\n";
        return false;
    }
    TTree *tree_in = dynamic_cast<TTree *>(infile->Get("events"));
    if (!tree_in) {
        std::cerr << "Replay: input raw file has no 'events' tree\n";
        return false;
    }
    tree_in->SetBranchStatus("*", 0);
    auto enable_branch = [tree_in](const char *name) {
        if (tree_in->GetBranch(name)) tree_in->SetBranchStatus(name, 1);
    };

    for (const char *name : {
             "event_num", "trigger_type", "trigger_bits",
             "hycal.nch", "hycal.module_id", "hycal.module_type",
             "hycal.gain_factor",
             "hycal.npeaks",
             "hycal.peak_height",
             "hycal.peak_time",
             "hycal.peak_integral"}) {
        enable_branch(name);
    }

    auto in = std::make_unique<EventVars>();
    prad2::SetRawReadBranches(tree_in, *in);

    // loop over events in the input file
    long long nentries = tree_in->GetEntries();
    for (long long i = 0; i < nentries; ++i) {
        tree_in->GetEntry(i);
        if ((in->trigger_bits & prad2::TBIT_sum) == 0) continue;

        if (in->nch > 100) continue; // too many hits, likely not a clean event

        clusterer.Clear();

        struct PeakInfo {
            int module_id;
            float mod_x;
            float mod_y;
            float energy;
            float time;
            int peak_n;
        };
        std::vector<PeakInfo> valid_peaks;

        for (int j = 0; j < in->nch; ++j) {
            const auto *mod = hycal.module_by_id(in->module_id[j]);
            if (!mod || !mod->is_pwo4()) continue;

            // Per-ID gain correction: average of 2 LMS channels(LMS 2 and 3, 1 is not used).
            float gain = in->gain_factor[j];

            // timing offset for this module
            float time_offset = mod->time_offset;

            auto hc_win = hc_time_cuts.at(mod->index);

            if (cluster_cfg.seed_time_window > 0.f) {
                // Multi-pulse mode: push every peak inside the trigger
                // window into the clusterer; the seed-anchored timing
                // coincidence cut is applied inside HyCalCluster.
                for (int p = 0; p < in->npeaks[j] && p < fdec::MAX_PEAKS; ++p) {
                    float peak_time = in->peak_time[j][p] - time_offset;
                    if (peak_time <= hc_win.lo) continue;
                    if (peak_time >= hc_win.hi) continue;
                    float adc = in->peak_integral[j][p] * gain;
                    float energy = static_cast<float>(mod->energize(adc));
                    clusterer.AddHit(mod->index, energy, peak_time);
                    valid_peaks.push_back({static_cast<int>(in->module_id[j]),
                        static_cast<float>(mod->x), static_cast<float>(mod->y),
                        energy, peak_time, in->npeaks[j]});
                }
            } else {
                // Legacy: pick the largest in-window peak as the single
                // module hit, time field unused downstream.
                int bestIdx = -1;
                float bestHeight = -1.f;
                for (int p = 0; p < in->npeaks[j] && p < fdec::MAX_PEAKS; ++p) {
                    float peak_time = in->peak_time[j][p] - time_offset;
                    if (peak_time <= hc_win.lo) continue;
                    if (peak_time >= hc_win.hi) continue;
                    if (in->peak_integral[j][p] > bestHeight) {
                        bestHeight = in->peak_integral[j][p];
                        bestIdx = p;
                    }
                }
                if (bestIdx < 0) continue;
                float adc = in->peak_integral[j][bestIdx] * gain;
                float energy = static_cast<float>(mod->energize(adc));
                clusterer.AddHit(mod->index, energy, in->peak_time[j][bestIdx] - time_offset);
                valid_peaks.push_back({static_cast<int>(in->module_id[j]),
                    static_cast<float>(mod->x), static_cast<float>(mod->y),
                    energy, in->peak_time[j][bestIdx] - time_offset, in->npeaks[j]});
            }
        }
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);

        if (hits.size() != 1) continue; // only keep single-cluster events
        if (hits[0].nblocks <= 3) continue; // require cluster to be at least 4 blocks (5x5) for this calibration
        
        auto *mod = hycal.module_by_id(hits[0].center_id);
        if (!mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals

        if (fdec::test_bit(hits[0].flag, fdec::kDeadModule)) continue; // skip clusters with dead modules

        // require hit to be in central 3x3 of a 5x5 grid in single central module (|xd|,|yd| < 0.3)
        float xd = (hits[0].x - (float)mod->x) / (float)mod->size_x;
        float yd = (hits[0].y - (float)mod->y) / (float)mod->size_y;
        if (std::abs(xd) >= 0.3f || std::abs(yd) >= 0.3f) continue;
        if (fdec::test_bit(hits[0].flag, fdec::kTransition)) {
            if (hits[0].x >  300.0 && xd >= 0.0f) continue; // only keep hits on the inner side for transition modules
            if (hits[0].x < -300.0 && xd <= 0.0f) continue;
            if (hits[0].y >  300.0 && yd >= 0.0f) continue;
            if (hits[0].y < -300.0 && yd <= 0.0f) continue;
        }

        // 5×5 energy sum: select modules whose center lies within
        // ±2 crystal pitches (20.75 mm) in both x and y from center
        constexpr float crystal_pitch = 20.75f;
        constexpr float half_win = 2.5f * crystal_pitch;
        float E5x5 = 0.f, center_energy = 0.f;
        float second_layer_energy = 0.f, third_layer_energy = 0.f;
        for (const auto &peak : valid_peaks) {
            float dx = peak.mod_x - mod->x;
            float dy = peak.mod_y - mod->y;
            if (std::abs(dx) <= half_win && std::abs(dy) <= half_win) {
                if (cluster_cfg.seed_time_window > 0.f &&
                    std::abs(peak.time - hits[0].time) > cluster_cfg.seed_time_window) {
                    continue; // skip peaks outside the seed time window
                }
                E5x5 += peak.energy;
                int x_offset = static_cast<int>(std::lround(dx / crystal_pitch));
                int y_offset = static_cast<int>(std::lround(dy / crystal_pitch));
                int layer = std::max(std::abs(x_offset), std::abs(y_offset));
                if (layer == 0) {
                    center_energy += peak.energy;
                } else if (layer == 1) {
                    second_layer_energy += peak.energy;
                } else if (layer == 2) {
                    third_layer_energy += peak.energy;
                }
            }
        }
        // require center module to have at least 60% of cluster energy
        if (hits[0].energy <= 0.f || E5x5 <= 0.f) continue;
        float center_energy_fraction = center_energy / hits[0].energy;
        if (center_energy_fraction < 0.6f) continue;

        float second_energy_fraction = second_layer_energy / hits[0].energy;
        float third_energy_fraction = third_layer_energy / hits[0].energy;

        res->h1_E_modules[mod->id-1001]->Fill(E5x5);
        float theta = std::atan2(std::sqrt(hits[0].x*hits[0].x + hits[0].y*hits[0].y), gRunConfig.hycal_z) * 180.0f / M_PI;
        res->h2_energy_theta->Fill(theta, E5x5);
        res->hit_pos->Fill(hits[0].x, hits[0].y);
        res->h_E_1cl->Fill(E5x5);
        res->h_center_energy_fraction->Fill(center_energy_fraction);
        res->h_center_energy->Fill(center_energy);
        res->h_2nd_energy_fraction->Fill(second_energy_fraction);
        res->h_3rd_energy_fraction->Fill(third_energy_fraction);
        res->events_processed++;
    }
    infile->Close();
    delete infile;
    return true;
}
