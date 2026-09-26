// physics_calib.cpp: multi-threaded HyCal physics calibration
// Raw ROOT files are processed in rounds, with one file assigned to each
// worker thread. Each worker reconstructs HyCal clusters and accumulates a
// 5x5 energy spectrum for every PbWO4 module. The per-thread histograms are
// merged before the module fits and calibration update are performed.
//
// Events are selected from the sum trigger by requiring one reconstructed
// cluster with at least three blocks, seeded by a PbWO4 module and flagged
// neither dead nor split. Unless -a is given, clusters between 2 and 16 module
// pitches from the beam must hit the central |xd|,|yd| < 0.3 of the seed
// crystal. For transition modules, only hits on the inner side of the crystal
// are retained. The cluster's 5x5 energy (energy_square) fills the seed
// module's spectrum.
//
// For each non-dead module with at least 40 entries, the reconstructed elastic
// e-p peak is fitted near the expected energy calculated from the run beam
// energy and detector geometry. The expected/fitted peak ratio is damped to
// 85% of the full correction, limited to [0.5, 2.0], and applied to the current
// calibration constant. Later iterations use the preceding result. Dead
// modules are not calibrated, and dead/dead-neighbor flags are recorded in the
// fit-result JSON.
//=============================================================================
//
// Usage: physics_calib <input_raw.root|dir> [more files/dirs...]
//                      [-i iteration] [-o output_dir]
//                      [-c seed_calib.json] [-j num_threads]
//                      [-f gaus|crystalball] [-a]
//   - input_raw.root|dir: input ROOT file or directory containing *_raw.root
//   - iteration: calibration iteration, starting from 1 (default: 1)
//   - output_dir: base output directory (default: current directory)
//   - seed_calib.json: input calibration for iteration 1
//                      (default: the run's calibration file in runinfo)
//   - num_threads: number of worker threads (default: 4)
//   - -f: peak-fit function, gaus or crystalball (default: gaus)
//   - -a: disable the central-region cut
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
#include "RunInfoConfig.h"
#include "gain_factor.h"
#include "PipelineBuilder.h"
#include "ToolUtils.h"

#include <TFile.h>
#include <TTree.h>

#include <iostream>
#include <fstream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <filesystem>
#include <vector>
#include <thread>
#include <mutex>
#include <memory>
#include <algorithm>

namespace fs = std::filesystem;

using EventVars = prad2::RawEventData;
using namespace analysis;

// ── Per-thread accumulated results ──────────────────────────────────────────
struct HistResult {
    HistList                                 all;   // every histogram below, in booking order
    std::vector<std::unique_ptr<TH1F>>       h1_E_modules;   // indexed by module ID (index 0 is module W1)
    std::vector<std::unique_ptr<TH1F>>       h1_E_modules_island;  // same as h1_E_modules but energy from island clustering
    std::unique_ptr<TH2F>                    h2_energy_theta;
    std::unique_ptr<TH2F>                    hit_pos;
    std::unique_ptr<TH1F>                    h_E_1cl;
    std::unique_ptr<TH1F>                    h_center_energy_fraction;
    std::unique_ptr<TH1F>                    h_center_energy;
    std::unique_ptr<TH1F>                    h_fit_peak_energy;
    std::unique_ptr<TH1F>                    h_fit_peak_ratio;
    std::unique_ptr<TH1F>                    h_fit_peak_chi2ndf;
    std::unique_ptr<TH1F>                    h_fit_peak_sigma;
    long long                                events_processed = 0;
    // more histograms here but do not go to the physics_calib_viewer
    std::unique_ptr<TH1F>                    h1_E_1cl_island;
    std::unique_ptr<TH1F>                    h1_E_1cl_square;
    std::unique_ptr<TH1F>                    h1_dE_1cl; // E_island - E_square
    std::unique_ptr<TH2F>                    h2_cl_module_occupancy; // the occupancy of neighboring modules in one cluster
};

const int angle_bins = 50; const double angle_min = 0., angle_max = 5.;
const int energy_bins = 500; const double energy_min = 0., energy_max = 5000.;
const int pos_bins = 720; const double pos_min = -360., pos_max = 360.;
const int center_energy_fraction_bins = 100; const double center_energy_fraction_min = 0., center_energy_fraction_max = 1.;
const int fit_ratio_bins = 200; const double fit_ratio_min = 0., fit_ratio_max = 2.;
const int fit_chi2ndf_bins = 200; const double fit_chi2ndf_min = 0., fit_chi2ndf_max = 50.;
const int fit_sigma_bins = 200; const double fit_sigma_min = 0., fit_sigma_max = 200.;
const int extra_energy_bins = 500; const double extra_energy_min = 0., extra_energy_max = 5000.;
const int extra_denergy_bins = 1000; const double extra_denergy_min = -100., extra_denergy_max = 100.;
const int cl_module_occupancy_bins = 9; const double cl_module_occupancy_min = -4.5, cl_module_occupancy_max = 4.5;

bool ProcessRawFiles (const std::string &input_raw, const RunConfig &run_cfg,
                      const std::string &db_dir, const std::string &recon_config_file,
                      const std::string &calib_file, HistResult *res, bool central_region);

// One histogram set.  The merged set (sfx "") is written under these names,
// which physics_calib_viewer.py looks up; per-thread sets only need unique names.
static std::unique_ptr<HistResult> makeHistResult(const std::string &sfx)
{
    auto res = std::make_unique<HistResult>();
    const char *s = sfx.c_str();

    res->h1_E_modules.resize(1156);
    res->h1_E_modules_island.resize(1156);
    for (int i = 0; i < 1156; ++i) {
        const int mod_id = i + 1000 + 1; // module IDs start at 1001(W1)
        res->h1_E_modules[i] = Book<TH1F>(res->all,
            Form("h1_E_mod_%d_merged%s", mod_id, s),
            Form("Module W%d cluster energy;E (MeV);Counts", mod_id-1000),
            energy_bins, energy_min, energy_max);
        res->h1_E_modules_island[i] = Book<TH1F>(res->all,
            Form("h1_E_mod_%d_island_merged%s", mod_id, s),
            Form("Module W%d island cluster energy;E (MeV);Counts", mod_id-1000),
            energy_bins, energy_min, energy_max);
    }
    res->h2_energy_theta = Book<TH2F>(res->all,
        Form("h2_energy_theta_merged%s", s),
        "Cluster energy vs theta;#theta (deg);E (MeV)",
        angle_bins, angle_min, angle_max,
        energy_bins, energy_min, energy_max);
    res->hit_pos = Book<TH2F>(res->all,
        Form("hit_pos_merged%s", s),
        "Hit position;X (mm);Y (mm)",
        pos_bins, pos_min, pos_max,
        pos_bins, pos_min, pos_max);
    res->h_E_1cl = Book<TH1F>(res->all,
        Form("h_E_1cl_merged%s", s),
        "Single-cluster energy;E (MeV);Counts",
        energy_bins, energy_min, energy_max);
    res->h2_cl_module_occupancy = Book<TH2F>(res->all,
        Form("h2_cl_module_occupancy%s", s),
        "Modules in selected cluster;#Delta x (module);#Delta y (module)",
        cl_module_occupancy_bins, cl_module_occupancy_min, cl_module_occupancy_max,
        cl_module_occupancy_bins, cl_module_occupancy_min, cl_module_occupancy_max);
    res->h_center_energy_fraction = Book<TH1F>(res->all,
        Form("h_center_energy_fraction%s", s),
        "Center energy fraction;E_{center}/E_{cluster};Counts",
        center_energy_fraction_bins,
        center_energy_fraction_min, center_energy_fraction_max);
    res->h_center_energy = Book<TH1F>(res->all,
        Form("h_center_energy%s", s),
        "Center module energy;E_{center} (MeV);Counts",
        energy_bins, energy_min, energy_max);
    res->h_fit_peak_energy = Book<TH1F>(res->all,
        Form("h_fit_peak_energy%s", s),
        "Fitted peak energy;E_{peak} (MeV);Modules",
        energy_bins, energy_min, energy_max);
    res->h_fit_peak_ratio = Book<TH1F>(res->all,
        Form("h_fit_peak_ratio%s", s),
        "Calibration ratio;E_{expected}/E_{peak};Modules",
        fit_ratio_bins, fit_ratio_min, fit_ratio_max);
    res->h_fit_peak_chi2ndf = Book<TH1F>(res->all,
        Form("h_fit_peak_chi2ndf%s", s),
        "Peak-fit #chi^{2}/NDF;#chi^{2}/NDF;Modules",
        fit_chi2ndf_bins, fit_chi2ndf_min, fit_chi2ndf_max);
    res->h_fit_peak_sigma = Book<TH1F>(res->all,
        Form("h_fit_peak_sigma%s", s),
        "Fitted peak sigma;#sigma (MeV);Modules",
        fit_sigma_bins, fit_sigma_min, fit_sigma_max);
    res->h1_E_1cl_island = Book<TH1F>(res->all,
        Form("h1_E_1cl_island%s", s),
        "Single-cluster island energy;E_{island} (MeV);Counts",
        extra_energy_bins, extra_energy_min, extra_energy_max);
    res->h1_E_1cl_square = Book<TH1F>(res->all,
        Form("h1_E_1cl_square%s", s),
        "Single-cluster square energy;E_{square} (MeV);Counts",
        extra_energy_bins, extra_energy_min, extra_energy_max);
    res->h1_dE_1cl = Book<TH1F>(res->all,
        Form("h1_dE_1cl%s", s),
        "Island minus square energy;E_{island}-E_{square} (MeV);Counts",
        extra_denergy_bins, extra_denergy_min, extra_denergy_max);
    return res;
}

int main(int argc, char *argv[])
{
    analysis::InitRootThreading();

    // ── Argument parsing ─────────────────────────────────────────────────────
    std::string output_path, seed_calib_file;
    int  iteration   = 1;
    int  num_threads = 4;
    bool use_crystal_ball = false;
    bool central_region = true;

    std::string db_dir = prad2::database_dir();

    std::string recon_config_file = db_dir + "/reconstruction_config.json";

    int opt;
    while ((opt = getopt(argc, argv, "i:o:j:c:f:a")) != -1) {
        switch (opt) {
            case 'i': iteration        = std::atoi(optarg); break;
            case 'o': output_path = optarg; break;
            case 'c': seed_calib_file  = optarg; break;
            case 'j': num_threads      = std::atoi(optarg); break;
            case 'a': central_region = false; break;
            case 'f': {
                std::string mode = optarg;
                std::transform(mode.begin(), mode.end(), mode.begin(), [](unsigned char ch) {
                    return std::tolower(ch);
                });
                if (mode == "gaus" || mode == "g") {
                    use_crystal_ball = false;
                } else if (mode == "crystalball" || mode == "crystal_ball" || mode == "cb" || mode == "c") {
                    use_crystal_ball = true;
                } else {
                    std::cerr << "Unknown fit mode: '" << optarg
                              << "'. Use 'gaus' or 'crystalball'.\n";
                    return 1;
                }
                break;
            }
        }
    }

    std::vector<std::string> root_files = CollectInputs(argc, argv, optind, IsRawRootName);
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: physics_calib <input_raw.root|dir> [more...] "
                     "[-i iter] [-o output_dir] [-c seed_calib.json] [-j threads] "
                     "[-f gaus|crystalball] [-a(central_region = false)]\n";
        return 1;
    }

    // ── Run number / output paths ─────────────────────────────────────────────
    int run_num = get_run_int(root_files[0]);
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
    std::string run_str = "run" + std::to_string(run_num);
    if (output_path.empty()) output_path = ".";
    std::string run_out_dir = output_path + "/Physics_calib/" + run_str;
    fs::create_directories(run_out_dir);
    // Make absolute: PipelineBuilder treats a relative calib path as
    // relative to db_dir, which would misresolve calib_factor_iterN.json.
    run_out_dir = fs::absolute(run_out_dir).lexically_normal().string();
    std::cerr << "Output directory: " << run_out_dir << "\n";

    std::string input_calib_file, output_calib_file, output_root_file, output_json_file;
    if (iteration == 1)
        input_calib_file = fs::absolute(!seed_calib_file.empty()
            ? fs::path(seed_calib_file)
            : fs::path(db_dir) / gRunConfig.energy_calib_file).lexically_normal().string();
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

    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");
    if (hycal.LoadCalibration(input_calib_file) <= 0) {
        std::cerr << "Cannot load calibration " << input_calib_file << "\n";
        return 1;
    }

    // ── Thread count ─────────────────────────────────────────────────────────
    int n_files    = static_cast<int>(root_files.size());
    num_threads    = std::max(1, std::min(num_threads, n_files));
    int num_rounds = (n_files + num_threads - 1) / num_threads;
    std::cout << "Processing " << n_files << " file(s) with "
              << num_threads << " thread(s), " << num_rounds << " round(s)\n";

    // ── Initialize per-thread results (once, reused across rounds) ───────────
    std::vector<std::unique_ptr<HistResult>> results(num_threads);
    std::mutex io_mtx;

    for (int tid = 0; tid < num_threads; ++tid)
        results[tid] = makeHistResult(Form("_tid%d", tid));

    // ── Process files in rounds: num_threads files per round, 1 file/thread ──
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

                bool ok = ProcessRawFiles(root_files[fi], gRunConfig,
                    db_dir, recon_config_file, input_calib_file, res, central_region);
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
    auto merged = makeHistResult("");
    for (const auto &res : results) {
        AddAll(merged->all, res->all);
        merged->events_processed += res->events_processed;
    }

    // resolve new calibration constants from the histograms of each module's energy distribution
    prad2::ApplyHyCalDeadModules(gRunConfig.hycal_dead_modules, hycal);
    analysis::PhysicsTools physics(hycal);

    // per-module results, written to calib_result_iterN.json
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
    int n_good_fit = 0;
    for (int i = 0; i < 1156; ++i) {
        if (!merged->h1_E_modules[i]) continue;
        TH1F *h = merged->h1_E_modules[i].get();
        if (h->GetEntries() < 40) continue;

        int mod_id = i + 1000 + 1;
        auto mod = hycal.module_by_id(mod_id);
        if (!mod) continue;

        bool is_dead = fdec::test_bit(mod->flag, fdec::kDeadModule);
        bool is_deadNeighbor = fdec::test_bit(mod->flag, fdec::kDeadNeighbor);
        if (is_dead) {
            calib_results.push_back({mod_id, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, 0.f, false, is_dead, is_deadNeighbor});
            continue;
        }

        float theta_deg = std::atan(std::sqrt(mod->x * mod->x + mod->y * mod->y)
                                    / gRunConfig.hycal_z) * 180.f / 3.14159265f;
        float expected_peak = analysis::PhysicsTools::ExpectedEnergy(theta_deg, gRunConfig.Ebeam, "ep");

        auto fit_result = physics.fitPeak(h, expected_peak, false, use_crystal_ball);
        float peak = static_cast<float>(fit_result[0]);
        float sigma = static_cast<float>(fit_result[1]);
        float chi2 = static_cast<float>(fit_result[2]);
        float expected_sigma = 0.03f*peak/std::sqrt(peak/1000.f);
        bool fit_good = (peak > 0 && sigma > 0.5f * expected_sigma && sigma < 1.5f * expected_sigma && chi2 < 2.5f);
        if (!fit_good) {
            std::cout << "Check!!! Module W" << (mod_id - 1000)
                 << ": fit failed (peak=" << peak
                 << ", sigma=" << sigma
                 << ", chi2/ndf=" << chi2 << ")\n";
        } else {
            n_good_fit++;
        }
        if (peak <= 0) peak = expected_peak; // fallback to expected if fit failed
        float ratio         = expected_peak / peak; 
        ratio = (ratio - 1.f) * 0.85f + 1.f; // apply a conservative factor to avoid over-correction
        if(ratio < 0.5) ratio = 0.5;
        if(ratio > 2.0) ratio = 2.0;

        merged->h_fit_peak_energy->Fill(peak);
        merged->h_fit_peak_ratio->Fill(ratio);
        merged->h_fit_peak_chi2ndf->Fill(chi2);
        merged->h_fit_peak_sigma->Fill(sigma);

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

    std::cout << "Calibration iteration " << iteration << " completed. "
              << n_calibrated << " modules calibrated. " << n_good_fit << " fits were good.\n";

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
    outfile->cd();
    outfile->mkdir("modules_island");
    outfile->cd("modules_island");
    for (int i = 0; i < 1156; ++i) {
        if (merged->h1_E_modules_island[i]) merged->h1_E_modules_island[i]->Write();
    }
    outfile->cd();
    outfile->mkdir("modules_5by5");
    outfile->cd("modules_5by5");
    for (int i = 0; i < 1156; ++i) {
        if (merged->h1_E_modules[i]) merged->h1_E_modules[i]->Write();
    }
    outfile->cd();
    if (merged->h2_energy_theta) merged->h2_energy_theta->Write();
    if (merged->hit_pos) merged->hit_pos->Write();
    if (merged->h_E_1cl) merged->h_E_1cl->Write();
    if (merged->h_center_energy_fraction) merged->h_center_energy_fraction->Write();
    if (merged->h_center_energy) merged->h_center_energy->Write();
    if (merged->h_fit_peak_energy) merged->h_fit_peak_energy->Write();
    if (merged->h_fit_peak_ratio) merged->h_fit_peak_ratio->Write();
    if (merged->h_fit_peak_chi2ndf) merged->h_fit_peak_chi2ndf->Write();
    if (merged->h_fit_peak_sigma) merged->h_fit_peak_sigma->Write();
    if (merged->h1_E_1cl_island) merged->h1_E_1cl_island->Write();
    if (merged->h1_E_1cl_square) merged->h1_E_1cl_square->Write();
    if (merged->h1_dE_1cl) merged->h1_dE_1cl->Write();
    if (merged->h2_cl_module_occupancy) {
        int entries = merged->h2_cl_module_occupancy->GetBinContent(cl_module_occupancy_bins/2+1, cl_module_occupancy_bins/2+1);
        if (entries > 0) merged->h2_cl_module_occupancy->Scale(100.0 / entries);
        merged->h2_cl_module_occupancy->Write();
    }

    outfile->Close();
    delete outfile;
}


bool ProcessRawFiles (const std::string &input_raw, const RunConfig &run_cfg,
                      const std::string &db_dir, const std::string &recon_config_file,
                      const std::string &calib_file, HistResult *res, bool central_region)
{   
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
    const auto &hycal        = pipeline.hycal;
    const auto &cluster_cfg  = pipeline.hycal_cluster_cfg;
    const auto &hc_time_cuts = pipeline.hycal_time_cuts;

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);

    fdec::WaveAnalyzer ana(pipeline.daq_cfg.wave_cfg);
    fdec::WaveResult wres;

    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(run_cfg, run_num);

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
    const bool has_waveform = tree_in->GetBranch("hycal.samples") != nullptr;
    const bool has_peaks    = tree_in->GetBranch("hycal.npeaks") != nullptr;
    std::cout << "Input file: " << input_raw
              << " (waveform branches: " << (has_waveform ? "yes" : "no")
              << ", peak branches: " << (has_peaks ? "yes" : "no") << ")\n";
    tree_in->SetBranchStatus("*", 0);
    auto enable_branch = [tree_in](const char *name) {
        if (tree_in->GetBranch(name)) tree_in->SetBranchStatus(name, 1);
    };

    for (const char *name : {
             "event_num", "trigger_type", "trigger_bits",
             "hycal.nch", "hycal.module_id", "hycal.module_type",
             "hycal.gain_factor",
             "hycal.nsamples", "hycal.samples",
             "hycal.npeaks",
             "hycal.peak_height",
             "hycal.peak_time",
             "hycal.peak_integral"}) {
        enable_branch(name);
    }

    auto in = std::make_unique<EventVars>();
    prad2::SetRawReadBranches(tree_in, *in);

    long long nentries = tree_in->GetEntries();
    for (long long i = 0; i < nentries; ++i) {
        tree_in->GetEntry(i);
        if ((in->trigger_bits & prad2::TBIT_sum) == 0) continue;

        if (in->nch > 100) continue; // too many hits, likely not a clean event

        clusterer.Clear();

        // raw trees written without peak branches: re-derive the peaks from the samples
        if (has_waveform && !has_peaks) FillPeaksFromWaveforms(*in, hycal, ana, wres);

        // Per-event gain correction (time-series lookup by event number).
        const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(in->event_num));

        for (int j = 0; j < in->nch; ++j) {
            const auto *mod = hycal.module_by_id(in->module_id[j]);
            if (!mod || !mod->is_pwo4()) continue;

            float gain = gain_corr.ModuleGain(mod->id);
            if (gain <= 0.f || gain == 1.f) gain = in->gain_factor[j];

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
            }
        }
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);

        if (hits.size() != 1) continue; // only keep single-cluster events
        if (hits[0].nblocks < 3) continue;
        
        auto *mod = hycal.module_by_id(hits[0].center_id);
        if (!mod || !mod->is_pwo4()) continue;

        if (fdec::test_bit(hits[0].flag, fdec::kDeadModule)) continue;
        if (fdec::test_bit(hits[0].flag, fdec::kSplit)) continue;

        // require hit to be in central 3x3 of a 5x5 grid in single central module (|xd|,|yd| < 0.3)
        const auto [xd, yd] = mod->cell_offset<float>(hits[0].x, hits[0].y);
        if ((std::abs(xd) >= 0.3f || std::abs(yd) >= 0.3f) && central_region
            && InHyCalRing(hits[0].x, hits[0].y, 2.0, 16.)) continue;
        if (fdec::test_bit(hits[0].flag, fdec::kTransition)) {
            if (hits[0].x >  300.0 && xd >= 0.0f) continue; // only keep hits on the inner side for transition modules
            if (hits[0].x < -300.0 && xd <= 0.0f) continue;
            if (hits[0].y >  300.0 && yd >= 0.0f) continue;
            if (hits[0].y < -300.0 && yd <= 0.0f) continue;
        }

        float center_energy = 0.f;
        const fdec::ModuleCluster *selected_cluster = nullptr;
        for (const auto &cluster : clusterer.GetClusters()) {
            if (cluster.center.index == mod->index) {
                center_energy = cluster.center.energy;
                selected_cluster = &cluster;
                break;
            }
        }

        if (hits[0].energy <= 0.f || hits[0].energy_square <= 0.f) continue;
        float center_energy_fraction = center_energy / hits[0].energy;

        res->h1_E_modules[mod->id-1001]->Fill(hits[0].energy_square);
        res->h1_E_modules_island[mod->id-1001]->Fill(hits[0].energy);
        float theta = std::atan2(std::sqrt(hits[0].x*hits[0].x + hits[0].y*hits[0].y), run_cfg.hycal_z) * 180.0f / M_PI;
        res->h2_energy_theta->Fill(theta, hits[0].energy_square);
        res->hit_pos->Fill(hits[0].x, hits[0].y);
        res->h_E_1cl->Fill(hits[0].energy_square);
        res->h_center_energy_fraction->Fill(center_energy_fraction);
        res->h_center_energy->Fill(center_energy);
        res->events_processed++;

        if (InHyCalRing(mod->x, mod->y, 4.0, 14.)) {
            res->h1_E_1cl_island->Fill(hits[0].energy);
            res->h1_E_1cl_square->Fill(hits[0].energy_square);
            res->h1_dE_1cl->Fill(hits[0].energy - hits[0].energy_square);
            if (selected_cluster) {
                const auto &center_mod = hycal.module(selected_cluster->center.index);
                for (const auto &cluster_hit : selected_cluster->hits) {
                    const auto &hit_mod = hycal.module(cluster_hit.index);
                    double dx = 0., dy = 0.;
                    hycal.qdist(center_mod, hit_mod, dx, dy);
                    res->h2_cl_module_occupancy->Fill(dx, dy);
                }
            }
        }
    }
    infile->Close();
    delete infile;
    return true;
}
