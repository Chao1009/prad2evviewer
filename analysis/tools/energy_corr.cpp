// energy_corr.cpp : to correct the reconstructed energy non-uniformity depending on
// the position of the cluster within the calorimeter modules.
//
// Measure the reconstructed-energy response at different positions within each
// module and compare its fitted peak with the expected energy. Each module is divided into a grid,
// First try to a grid of 5 by 5, the map is below, could be saved into a 2D array of 1D hist
// with one reconstructed-energy histogram for each cell in the 5x5 grid.
// column   0   1   2   3   4
// row     +---+---+---+---+---+     beam top  ^
//  0      |   |   |   |   |   |               |
//  1      |   |   |   |   |   |
//  2      |   |   |   |   |   |     beam right ->
//  3      |   |   |   |   |   |
//  4      |   |   |   |   |   |
//         +---+---+---+---+---+
//
// Usage:
//   energy_corr [-a 3.5GeV_file1.root ...] [-b 2.2GeV_file1.root ...] [-c 0.7GeV_file1.root ...]
//               -o output.root [-n max_events]
//
// At least one of -a / -b / -c must be provided.  Each energy group is processed
// in a separate thread (up to 3 threads running in parallel).  Results are written
// into per-energy subdirectories inside the output ROOT file:
//   3p5GeV/fit_results, 3p5GeV/h_energy_*
//   2p2GeV/fit_results, 2p2GeV/h_energy_*
//   0p7GeV/fit_results, 0p7GeV/h_energy_*

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TMarker.h>
#include <TLegend.h>
#include <TROOT.h>
#include <TClass.h>
#include <TDirectory.h>

#include <nlohmann/json.hpp>

#include <fstream>
#include <iostream>
#include <array>
#include <map>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <thread>
#include <mutex>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;

using EventVars_Recon = prad2::ReconEventData;

constexpr int kGridSize = 5;
using HistGrid = std::array<std::array<TH1F*, kGridSize>, kGridSize>;
using EnergyHistMap = std::map<int, HistGrid>;

struct ExpectedEnergies {
    double ep = 0.;
    double ee = 0.;
};
using ExpectedEnergyGrid = std::array<std::array<ExpectedEnergies, kGridSize>, kGridSize>;
using ExpectedEnergyMap = std::map<int, ExpectedEnergyGrid>;

// Per-cell fit outcome stored in memory for JSON export.
struct CellFitData {
    double expected_ep = 0.;
    double expected_ee = 0.;
    double correction  = 1.;   // raw (uncapped) value; 0 means not fitted
    bool   fit_valid   = false;
};
using FitDataGrid   = std::array<std::array<CellFitData, kGridSize>, kGridSize>;
using FitDataByName = std::map<std::string, FitDataGrid>; // module_name → 5×5

// Per-energy-run state: histograms, expected energies, and metadata.
struct EnergyRun {
    std::string label;           // human-readable tag, e.g. "3.5GeV"
    float beam_energy;           // beam energy in MeV
    std::string subdir;          // subdirectory inside the output ROOT file
    std::vector<std::string> input_paths;
    EnergyHistMap energy_hists;
    ExpectedEnergyMap expected_energies;
    FitDataByName fit_data;      // populated by write_energy_run(), used for JSON
};

// ---------------------------------------------------------------------------
// Forward declarations
// ---------------------------------------------------------------------------
void init_energy_run(EnergyRun &run, const fdec::HyCalSystem &hycal);
void write_energy_run(EnergyRun &run, const fdec::HyCalSystem &hycal, TDirectory *dir);
void write_corr_json(const std::vector<EnergyRun> &runs,
                     const fdec::HyCalSystem &hycal,
                     const std::string &path);
long long process_event(TTree *tree, const EventVars_Recon &ev,
                        const fdec::HyCalSystem &hycal,
                        EnergyHistMap &energy_hists, float Ebeam,
                        int max_events = -1, const std::string &label = "",
                        std::mutex *io_mtx = nullptr);

float resolution = 0.035; // pre-defined energy resolution

constexpr float E3p5 = 3485.41f; // Energy for 3.5 GeV beam, MeV
constexpr float E2p2 = 2239.51f; // Energy for 2.2 GeV beam, MeV
constexpr float E0p7 =  728.9f;  // Energy for 0.7 GeV beam, MeV
constexpr float kHyCalZ = 6269.f; // HyCal distance from the target, mm

bool Vetoed(float cl_time, float sci_time, float sci_int)
{
    const float time_shift    = 35.f;   // ns
    const float time_window   =  7.f;   // ns
    const float int_threshold = 2000.f; // arbitrary units
    return (fabs(cl_time - sci_time - time_shift) < time_window)
           && (sci_int > int_threshold);
}

// ---------------------------------------------------------------------------
int main(int argc, char *argv[])
{
    // --- parse command line arguments ---
    std::vector<std::string> inputs_a, inputs_b, inputs_c;
    std::string output_path;
    std::string json_path;
    int max_events = -1;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-a" || arg == "-b" || arg == "-c") {
            auto &inputs = (arg == "-a") ? inputs_a
                         : (arg == "-b") ? inputs_b
                                         : inputs_c;
            if (++i >= argc || argv[i][0] == '-') {
                std::cerr << arg << " requires at least one input ROOT file.\n";
                return 1;
            }
            do {
                inputs.emplace_back(argv[i]);
                ++i;
            } while (i < argc && argv[i][0] != '-');
            --i;
        } else if (arg == "-o") {
            if (++i >= argc) {
                std::cerr << "-o requires an output ROOT file.\n";
                return 1;
            }
            output_path = argv[i];
        } else if (arg == "-j") {
            if (++i >= argc) {
                std::cerr << "-j requires an output JSON file.\n";
                return 1;
            }
            json_path = argv[i];
        } else if (arg == "-n") {
            if (++i >= argc) {
                std::cerr << "-n requires a maximum event count.\n";
                return 1;
            }
            max_events = std::atoi(argv[i]);
            if (max_events <= 0) {
                std::cerr << "-n must be a positive integer.\n";
                return 1;
            }
        } else {
            std::cerr << "Unknown argument: " << arg << "\n";
            return 1;
        }
    }

    if (inputs_a.empty() && inputs_b.empty() && inputs_c.empty()) {
        std::cerr << "At least one of -a / -b / -c must be specified.\n"
                  << "Usage: " << argv[0]
                  << " [-a 3.5GeV_files...] [-b 2.2GeV_files...]"
                     " [-c 0.7GeV_files...] -o output.root [-j output.json] [-n max_events]\n";
        return 1;
    }
    if (output_path.empty()) {
        std::cerr << "-o output.root is required.\n";
        return 1;
    }

    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- init detector system (shared read-only across all threads) ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");

    // --- build per-energy runs ---
    std::vector<EnergyRun> runs;
    runs.reserve(3);
    if (!inputs_a.empty())
        runs.push_back({"3.5GeV", E3p5, "3p5GeV", std::move(inputs_a), {}, {}});
    if (!inputs_b.empty())
        runs.push_back({"2.2GeV", E2p2, "2p2GeV", std::move(inputs_b), {}, {}});
    if (!inputs_c.empty())
        runs.push_back({"0.7GeV", E0p7, "0p7GeV", std::move(inputs_c), {}, {}});

    for (auto &run : runs)
        init_energy_run(run, hycal);

    // --- enable ROOT thread safety before launching threads ---
    ROOT::EnableThreadSafety();

    // --- launch one thread per energy group ---
    std::mutex io_mtx;
    std::vector<std::thread> threads;
    threads.reserve(runs.size());

    for (auto &run : runs) {
        threads.emplace_back([&run, &hycal, max_events, &io_mtx]() {
            TChain chain("recon");
            for (const auto &f : run.input_paths) {
                if (chain.Add(f.c_str()) == 0) {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cerr << "[" << run.label << "] Failed to add: " << f << "\n";
                }
            }
            if (chain.GetEntries() <= 0) {
                std::lock_guard<std::mutex> lk(io_mtx);
                std::cerr << "[" << run.label << "] No entries found – skipping.\n";
                return;
            }
            EventVars_Recon ev;
            prad2::SetReconReadBranches(&chain, ev);
            process_event(&chain, ev, hycal, run.energy_hists, run.beam_energy,
                          max_events, run.label, &io_mtx);
        });
    }

    for (auto &t : threads)
        t.join();

    // --- write output ROOT file ---
    TFile output_file(output_path.c_str(), "RECREATE");
    if (output_file.IsZombie()) {
        std::cerr << "Failed to create output ROOT file: " << output_path << '\n';
        return 1;
    }

    for (auto &run : runs) {
        TDirectory *dir = output_file.mkdir(run.subdir.c_str());
        write_energy_run(run, hycal, dir);
    }

    output_file.Close();

    if (!json_path.empty())
        write_corr_json(runs, hycal, json_path);

    return 0;
}

// ---------------------------------------------------------------------------
// init_energy_run – allocate histograms and compute expected energies for one
// energy configuration.  Histogram x-axis: 10 MeV/bin, range [0, 1.25 * Ebeam].
// ---------------------------------------------------------------------------
void init_energy_run(EnergyRun &run, const fdec::HyCalSystem &hycal)
{
    const float xmax     = run.beam_energy * 1.25f;
    const float bin_width = (run.beam_energy < 1000.f) ? 5.f : 10.f;
    const int   nbins    = std::max(50, static_cast<int>(xmax / bin_width));

    for (int i = 0; i < hycal.module_count(); ++i) {
        const auto &module = hycal.module(i);
        if (!module.is_pwo4()) continue;
        // exclude outmost layer modules and modules under the absorber
        if (fabs(module.x) > 16.f * 20.77f || fabs(module.y) > 16.f * 20.75f) continue;
        if (fabs(module.x) <  2.f * 20.77f && fabs(module.y) <  2.f * 20.75f) continue;

        auto &grid = run.energy_hists[module.id];
        for (int row = 0; row < kGridSize; ++row) {
            for (int col = 0; col < kGridSize; ++col) {
                const std::string name = "h_energy_" + module.name + "_r"
                    + std::to_string(row) + "_c" + std::to_string(col);
                const std::string title = module.name + " energy distribution, cell ("
                    + std::to_string(row) + ", " + std::to_string(col)
                    + ");E_{reconstructed} (MeV);Counts";
                grid[row][col] = new TH1F(name.c_str(), title.c_str(), nbins, 0., xmax);
                grid[row][col]->SetDirectory(nullptr);

                // Geometric centre of this cell (row 0 = top, col 0 = left).
                const double cell_x = module.x
                    + (static_cast<double>(col) + 0.5) / kGridSize * module.size_x
                    - 0.5 * module.size_x;
                const double cell_y = module.y
                    + 0.5 * module.size_y
                    - (static_cast<double>(row) + 0.5) / kGridSize * module.size_y;
                const double theta = std::atan2(std::hypot(cell_x, cell_y), kHyCalZ)
                    * 180. / M_PI;
                auto &expected = run.expected_energies[module.id][row][col];
                expected.ep = PhysicsTools::ExpectedEnergy(theta, run.beam_energy, "ep");
                expected.ee = PhysicsTools::ExpectedEnergy(theta, run.beam_energy, "ee");
            }
        }
    }
}

// ---------------------------------------------------------------------------
// write_energy_run – fit each cell's histogram and write fit_results TTree
// plus all histograms into *dir*.
// ---------------------------------------------------------------------------
void write_energy_run(EnergyRun &run, const fdec::HyCalSystem &hycal, TDirectory *dir)
{
    dir->cd();

    constexpr Long64_t kMinFitEntries = 30;

    TTree fit_results("fit_results", "Per-cell reconstructed-energy Gaussian fit results");
    fit_results.SetDirectory(dir);

    int         module_id          = 0;
    std::string module_name;
    int         row                = 0;
    int         col                = 0;
    Long64_t    entries            = 0;
    double      histogram_mean     = 0.;
    double      expected_energy_ep = 0.;
    double      expected_energy_ee = 0.;
    double      fit_mean           = 0.;
    double      fit_mean_error     = 0.;
    double      fit_sigma          = 0.;
    double      fit_sigma_error    = 0.;
    double      chi2_ndf           = 0.;
    double      correction         = 0.;
    double      correction_error   = 0.;
    int         root_fit_status    = -1;
    bool        fit_valid          = false;

    fit_results.Branch("module_id",          &module_id);
    fit_results.Branch("module_name",        &module_name);
    fit_results.Branch("row",                &row);
    fit_results.Branch("col",                &col);
    fit_results.Branch("entries",            &entries);
    fit_results.Branch("histogram_mean",     &histogram_mean);
    fit_results.Branch("expected_energy_ep", &expected_energy_ep);
    fit_results.Branch("expected_energy_ee", &expected_energy_ee);
    fit_results.Branch("fit_mean",           &fit_mean);
    fit_results.Branch("fit_mean_error",     &fit_mean_error);
    fit_results.Branch("fit_sigma",          &fit_sigma);
    fit_results.Branch("fit_sigma_error",    &fit_sigma_error);
    fit_results.Branch("chi2_ndf",           &chi2_ndf);
    fit_results.Branch("correction",         &correction);
    fit_results.Branch("correction_error",   &correction_error);
    fit_results.Branch("root_fit_status",    &root_fit_status);
    fit_results.Branch("fit_valid",          &fit_valid);

    Long64_t successful_fits = 0;

    for (auto &[id, grid] : run.energy_hists) {
        module_id = id;
        const auto module = hycal.module_by_id(module_id);
        module_name = module ? module->name : "";

        for (row = 0; row < kGridSize; ++row) {
            for (col = 0; col < kGridSize; ++col) {
                TH1F *hist = grid[row][col];
                entries        = static_cast<Long64_t>(hist->GetEntries());
                histogram_mean = hist->GetMean();
                expected_energy_ep = 0.;
                expected_energy_ee = 0.;
                fit_mean           = 0.;
                fit_mean_error     = 0.;
                fit_sigma          = 0.;
                fit_sigma_error    = 0.;
                chi2_ndf           = 0.;
                correction         = 0.;
                correction_error   = 0.;
                root_fit_status    = -1;
                fit_valid          = false;

                const auto expected_it = run.expected_energies.find(module_id);
                if (expected_it != run.expected_energies.end()) {
                    const auto &expected = expected_it->second[row][col];
                    expected_energy_ep = expected.ep;
                    expected_energy_ee = expected.ee;
                }

                if (entries >= kMinFitEntries && expected_energy_ep > 0.) {
                    const double expected_sigma = resolution * expected_energy_ep
                        / std::sqrt(expected_energy_ep / 1000.);
                    const double fit_lo = std::fmax(hist->GetXaxis()->GetXmin(),
                                                    expected_energy_ep - 2. * expected_sigma);
                    const double fit_hi = std::fmin(hist->GetXaxis()->GetXmax(),
                                                    expected_energy_ep + 2. * expected_sigma);
                    const std::string fit_name = "f_energy_" + module_name + "_r"
                        + std::to_string(row) + "_c" + std::to_string(col);
                    TF1 gaussian(fit_name.c_str(), "gaus", fit_lo, fit_hi);
                    gaussian.SetParameters(hist->GetMaximum(), expected_energy_ep, expected_sigma);
                    gaussian.SetParLimits(1, fit_lo, fit_hi);
                    gaussian.SetParLimits(2, 0.25 * hist->GetBinWidth(1),
                                          5. * expected_sigma);

                    // "N" keeps the short-lived TF1 out of the histogram ownership list.
                    root_fit_status = hist->Fit(&gaussian, "RQL0N");
                    const double mean        = gaussian.GetParameter(1);
                    const double mean_error  = gaussian.GetParError(1);
                    const double sigma       = std::fabs(gaussian.GetParameter(2));
                    const double sigma_error = gaussian.GetParError(2);

                    fit_valid = root_fit_status == 0
                        && std::isfinite(mean)        && mean > 0.
                        && std::isfinite(mean_error)  && mean_error >= 0.
                        && std::isfinite(sigma)       && sigma > 0.
                        && std::isfinite(sigma_error) && sigma_error >= 0.;
                    if (fit_valid) {
                        fit_mean         = mean;
                        fit_mean_error   = mean_error;
                        fit_sigma        = sigma;
                        fit_sigma_error  = sigma_error;
                        chi2_ndf         = gaussian.GetNDF() > 0
                            ? gaussian.GetChisquare() / gaussian.GetNDF() : 0.;
                        correction       = expected_energy_ep / fit_mean;
                        correction_error = expected_energy_ep * fit_mean_error
                            / (fit_mean * fit_mean);
                        ++successful_fits;
                    }
                }

                // Store fit outcome for JSON export.
                auto &cd       = run.fit_data[module_name][row][col];
                cd.expected_ep = expected_energy_ep;
                cd.expected_ee = expected_energy_ee;
                cd.correction  = correction;
                cd.fit_valid   = fit_valid;

                fit_results.Fill();
            }
        }
    }

    std::cerr << "[" << run.label << "] Fitted " << successful_fits
              << "/" << fit_results.GetEntries() << " module cells successfully.\n";

    dir->cd();
    fit_results.Write();

    for (auto &[mid, grid] : run.energy_hists) {
        (void)mid;
        for (auto &r : grid)
            for (TH1F *hist : r)
                hist->Write();
    }
}

// ---------------------------------------------------------------------------
// process_event – fill per-cell energy histograms from a TTree/TChain.
// ---------------------------------------------------------------------------
long long process_event(TTree *tree, const EventVars_Recon &ev,
                        const fdec::HyCalSystem &hycal,
                        EnergyHistMap &energy_hists, float Ebeam,
                        int max_events, const std::string &label,
                        std::mutex *io_mtx)
{
    auto log_msg = [&](const std::string &msg, bool flush = false) {
        if (io_mtx) {
            std::lock_guard<std::mutex> lk(*io_mtx);
            std::cerr << msg;
            if (flush) std::cerr << std::flush;
        } else {
            std::cerr << msg;
            if (flush) std::cerr << std::flush;
        }
    };

    long long n_accepted = 0;
    const Long64_t n_entries = tree->GetEntries();

    for (Long64_t i = 0; i < n_entries; ++i) {
        if (max_events > 0 && i >= max_events) {
            log_msg("[" + label + "] Reached max events limit: "
                    + std::to_string(max_events) + "\n");
            break;
        }
        tree->GetEntry(i);
        if (i % 10000 == 0) {
            log_msg("[" + label + "] Processing event " + std::to_string(i)
                    + "/" + std::to_string(n_entries) + "\r", true);
        }
        if ((ev.trigger_bits & prad2::TBIT_sum) == 0) continue;
        ++n_accepted;

        for (int j = 0; j < ev.n_clusters; j++) {
            int mod_id = ev.cl_center[j];
            if (ev.cl_nblocks[j] <= 3) continue;
            auto mod = hycal.module_by_id(mod_id);
            if (!mod || !mod->is_pwo4()) continue;
            auto hist_it = energy_hists.find(mod_id);
            if (hist_it == energy_hists.end()) continue;

            const float mod_x      = static_cast<float>(mod->x);
            const float mod_y      = static_cast<float>(mod->y);
            const float mod_size_x = static_cast<float>(mod->size_x);
            const float mod_size_y = static_cast<float>(mod->size_y);

            const float c_x = ev.cl_x[j];
            const float c_y = ev.cl_y[j];
            const float c_z = ev.cl_z[j];

            const float xd = (c_x - mod_x) / mod_size_x;
            const float yd = (c_y - mod_y) / mod_size_y;

            const float theta  = std::atan2(std::sqrt(c_x*c_x + c_y*c_y), c_z)
                                 * 180.f / M_PI;
            const float energy = ev.cl_energy[j];

            bool veto = false;
            float sci_time, sci_int;
            for (int k = 0; k < ev.veto_nch; k++) {
                for (int p = 0; p < ev.veto_npeaks[k]; p++) {
                    sci_time = ev.veto_peak_time[k][p];
                    sci_int  = ev.veto_peak_integral[k][p];
                    veto = Vetoed(ev.cl_time[j], sci_time, sci_int);
                    if (veto) break;
                }
                if (veto) break;
            }
            if (theta > 1.3f) veto = false;
            if (veto && energy > 600.f && Ebeam < 1000.f) continue;

            const int col_idx = static_cast<int>(std::floor((xd + 0.5f) * kGridSize));
            const int row_idx = static_cast<int>(std::floor((0.5f - yd) * kGridSize));
            if (row_idx < 0 || row_idx >= kGridSize) continue;
            if (col_idx < 0 || col_idx >= kGridSize) continue;

            hist_it->second[row_idx][col_idx]->Fill(energy);
        }
    }

    log_msg("[" + label + "] Done. Accepted events: " + std::to_string(n_accepted) + "\n");
    return n_accepted;
}

// ---------------------------------------------------------------------------
// write_corr_json – export per-cell correction factors to a JSON file.
//
// Output guarantees:
//   • Every eligible module is written, sorted by module ID (ascending).
//   • All three beam-energy subdirs (3p5GeV / 2p2GeV / 0p7GeV) are always
//     present even if that energy was not processed: expected energies are
//     computed from detector geometry, corrections default to 1.0.
//   • 5×5 grids are written as compact rows: [v0, v1, v2, v3, v4]
//
// Correction caps:
//   fit not valid / no entries → 1.0
//   correction > 1.08          → 1.08
//   correction < 0.92          → 0.92
// ---------------------------------------------------------------------------
void write_corr_json(const std::vector<EnergyRun> &runs,
                     const fdec::HyCalSystem &hycal,
                     const std::string &path)
{
    constexpr double kCapHi = 1.08;
    constexpr double kCapLo = 0.92;

    auto cap = [](double corr, bool valid) -> double {
        constexpr double hi = 1.08, lo = 0.92;
        if (!valid || corr <= 0. || !std::isfinite(corr)) return 1.;
        return std::max(lo, std::min(hi, corr));
    };

    // Number formatter: removes trailing zeros, always keeps ≥ 1 decimal.
    auto fmt = [](double v, int max_dec) -> std::string {
        char buf[32];
        std::snprintf(buf, sizeof(buf), "%.*f", max_dec, v);
        std::string s = buf;
        const auto dot = s.find('.');
        if (dot != std::string::npos) {
            auto last = s.find_last_not_of('0');
            if (last == dot) last = dot + 1;
            s.resize(last + 1);
        }
        return s;
    };

    // Standard beam-energy configs – always written in this fixed order.
    struct EConfig { std::string subdir; float beam_energy; };
    const std::vector<EConfig> kEnergies = {
        {"3p5GeV", E3p5},
        {"2p2GeV", E2p2},
        {"0p7GeV", E0p7},
    };

    // Lookup: subdir → run fit_data.
    std::map<std::string, const FitDataByName*> run_map;
    for (const auto &run : runs)
        run_map[run.subdir] = &run.fit_data;

    // Build name→id map from hycal, then id→name (sorted) from fit_data.
    std::map<std::string, int> name_to_id;
    for (int i = 0; i < hycal.module_count(); ++i) {
        const auto &m = hycal.module(i);
        name_to_id[std::string(m.name)] = m.id;
    }
    std::map<int, std::string> id_to_name;   // integer-sorted by module ID
    for (const auto &run : runs) {
        for (const auto &[mod_name, ignored] : run.fit_data) {
            auto it = name_to_id.find(mod_name);
            if (it != name_to_id.end())
                id_to_name[it->second] = mod_name;
        }
    }

    // Compute expected ep/ee for a single cell from detector geometry.
    auto cell_exp = [&](const auto &m, int r, int c, float Ebeam)
        -> std::pair<double, double>
    {
        const double cx = m.x + (c + 0.5) / kGridSize * m.size_x - 0.5 * m.size_x;
        const double cy = m.y + 0.5 * m.size_y - (r + 0.5) / kGridSize * m.size_y;
        const double theta = std::atan2(std::hypot(cx, cy), kHyCalZ) * 180. / M_PI;
        return {PhysicsTools::ExpectedEnergy(theta, Ebeam, "ep"),
                PhysicsTools::ExpectedEnergy(theta, Ebeam, "ee")};
    };

    std::ofstream out(path);
    if (!out) {
        std::cerr << "Cannot write JSON: " << path << "\n";
        return;
    }

    using Mat5 = std::array<std::array<double, kGridSize>, kGridSize>;

    // Write a 5×5 matrix with compact rows at fixed nesting depth:
    //   key at 8-space indent → rows at 10-space indent → closing ] at 8 spaces.
    auto write_mat5 = [&](const Mat5 &mat, int dec) {
        out << "[\n";
        for (int r = 0; r < kGridSize; ++r) {
            out << "          [";
            for (int c = 0; c < kGridSize; ++c) {
                if (c > 0) out << ", ";
                out << fmt(mat[r][c], dec);
            }
            out << "]";
            if (r < kGridSize - 1) out << ",";
            out << "\n";
        }
        out << "        ]";
    };

    out << "{\n";
    bool first_mod = true;

    for (const auto &[mod_id, mod_name] : id_to_name) {
        if (!first_mod) out << ",\n";
        first_mod = false;
        out << "  \"" << mod_name << "\": {\n";

        const auto *mod_ptr = hycal.module_by_id(mod_id);

        bool first_e = true;
        for (const auto &ec : kEnergies) {
            if (!first_e) out << ",\n";
            first_e = false;
            out << "    \"" << ec.subdir << "\": {\n";

            Mat5 ep_exp{}, ee_exp{}, corr_mat{};
            for (auto &row : corr_mat) row.fill(1.0);   // default: identity

            // Try to use stored fit results for this energy.
            bool has_fit = false;
            auto sd_it = run_map.find(ec.subdir);
            if (sd_it != run_map.end()) {
                auto mod_it = sd_it->second->find(mod_name);
                if (mod_it != sd_it->second->end()) {
                    for (int r = 0; r < kGridSize; ++r) {
                        for (int c = 0; c < kGridSize; ++c) {
                            const auto &cd = mod_it->second[r][c];
                            ep_exp[r][c]   = cd.expected_ep;
                            ee_exp[r][c]   = cd.expected_ee;
                            corr_mat[r][c] = cap(cd.correction, cd.fit_valid);
                        }
                    }
                    has_fit = true;
                }
            }
            // Fall back to geometry-derived expected energies for missing runs.
            if (!has_fit && mod_ptr) {
                for (int r = 0; r < kGridSize; ++r) {
                    for (int c = 0; c < kGridSize; ++c) {
                        auto [ep, ee] = cell_exp(*mod_ptr, r, c, ec.beam_energy);
                        ep_exp[r][c] = ep;
                        ee_exp[r][c] = ee;
                        // corr_mat already set to 1.0
                    }
                }
            }

            // Round values before writing.
            for (int r = 0; r < kGridSize; ++r) {
                for (int c = 0; c < kGridSize; ++c) {
                    ep_exp[r][c]   = std::round(ep_exp[r][c]   * 10.)    / 10.;
                    ee_exp[r][c]   = std::round(ee_exp[r][c]   * 10.)    / 10.;
                    corr_mat[r][c] = std::round(corr_mat[r][c] * 10000.) / 10000.;
                }
            }

            // ep section
            out << "      \"ep\": {\n";
            out << "        \"correction\": ";      write_mat5(corr_mat, 4); out << ",\n";
            out << "        \"expected_energy\": "; write_mat5(ep_exp,   1); out << "\n";
            out << "      },\n";

            // ee section
            out << "      \"ee\": {\n";
            out << "        \"correction\": ";      write_mat5(corr_mat, 4); out << ",\n";
            out << "        \"expected_energy\": "; write_mat5(ee_exp,   1); out << "\n";
            out << "      }\n";

            out << "    }";
        }
        out << "\n  }";
    }
    out << "\n}\n";

    std::cerr << "Written correction JSON (" << id_to_name.size()
              << " modules) to " << path << "\n";
}
