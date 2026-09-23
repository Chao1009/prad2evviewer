// energy_corr.cpp : to correct the reconstructed energy non-uniformity depending on
// the position of the cluster within the calorimeter modules.
//
// Measure the reconstructed-energy response at different positions within each
// module and compare its fitted peak with the expected energy. Each module is divided into a grid,
// First try to a grid of 5 by 5, the map is below, could be saved into a 2D array of 1D hist
// with one reconstructed-energy histogram for each cell in the 5x5 grid.
// The histograms use h1_energy_grid[module][column][row]. Columns increase
// with the local HyCal X coordinate, and rows increase with local HyCal Y.
// column   0   1   2   3   4
// row     +---+---+---+---+---+     beam top (+Y) ^
//  4      |   |   |   |   |   |                   |
//  3      |   |   |   |   |   |
//  2      |   |   |   |   |   |     beam right (+X) ->
//  1      |   |   |   |   |   |
//  0      |   |   |   |   |   |
//         +---+---+---+---+---+

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "MatchingTools.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TF1.h>
#include <TF2.h>
#include <TGraphErrors.h>
#include <TKey.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPad.h>
#include <TROOT.h>
#include <TClass.h>
#include <TLorentzVector.h>

#include <nlohmann/json.hpp>

#include <iostream>
#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <atomic>
#include <cstdio>
#include <future>
#include <map>
#include <memory>
#include <mutex>
#include <limits>
#include <thread>
#include <getopt.h>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;

//anagles bin edges
const int Nbins = 33;
const Double_t binEdge[Nbins+1] = {
    0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.775, 0.800, 0.825, 0.850,
    0.875, 0.900, 0.940, 0.975, 1.014, 1.057, 1.105, 1.157, 1.211, 1.270,
    1.338, 1.417, 1.514, 1.634, 1.787, 2.000, 2.213, 2.492, 2.792, 3.092,
    3.392, 3.692, 3.992, 4.292
};
const int energy_bins = 350; const double energy_min = 500., energy_max = 4000.;
const int grids = 5;

// Modules map that want to draw
// W456 W457 W458 W459 W460 W461 W462 W463
// W490 W491 W492 W493 W494 W495 W496 W497
// W524 W525                     W530 W531
// W558 W559                     W564 W565
// W592 W593                     W598 W599
// W626 W627                     W632 W633
// W660 W661 W662 W663 W664 W665 W666 W667
// W694 W695 W696 W697 W698 W699 W700 W701
const std::array<int, 48> module_numbers = {
    456, 457, 458, 459, 460, 461, 462, 463,
    490, 491, 492, 493, 494, 495, 496, 497,
    524, 525, 530, 531,
    558, 559, 564, 565,
    592, 593, 598, 599,
    626, 627, 632, 633,
    660, 661, 662, 663, 664, 665, 666, 667,
    694, 695, 696, 697, 698, 699, 700, 701
};
const std::array<int, 48> module_canvas_pads = {
     1,  2,  3,  4,  5,  6,  7,  8,
     9, 10, 11, 12, 13, 14, 15, 16,
    17, 18, 23, 24,
    25, 26, 31, 32,
    33, 34, 39, 40,
    41, 42, 47, 48,
    49, 50, 51, 52, 53, 54, 55, 56,
    57, 58, 59, 60, 61, 62, 63, 64
};
const int module_count = module_numbers.size();

struct HistResult {
    std::unique_ptr<TH2F> h2_hit_module_hycal;
    std::unique_ptr<TH2F> h2_hit_module_gem;
    std::unique_ptr<TH1F> h1_energy_grid[module_count][grids][grids];
    std::unique_ptr<TH1F> h1_energy_grid_allModule[1156][grids][grids];
    long long events_processed = 0;
};

struct SharedFillLocks {
    std::array<std::mutex, 1156> module;
    std::mutex hit_maps;
    std::atomic<long long> events_processed{0};
};

static std::unique_ptr<HistResult> makeHistResult(const std::string &suffix);
static bool processRootFile(const std::string &input_file, const RunConfig &run_config,
                            const std::string &db_dir, long long max_events,
                            HistResult *result, SharedFillLocks *fill_locks);

static std::vector<std::string> collectRootFiles(const std::string &path);

static std::string shell_quote(const std::string &value)
{
    std::string quoted = "'";
    for (char ch : value) {
        if (ch == '\'') quoted += "'\\''";
        else quoted += ch;
    }
    return quoted + "'";
}

static std::string outputFileName(const std::string &output_name, bool corr = false)
{
    const fs::path output_path(output_name);
    const std::string file_name = output_path.filename().string()
        + (corr ? ".corr" : "") + ".root";
    return (output_path.parent_path() / file_name).string();
}

static std::string outputJsonFileName(const std::string &output_name, bool corr = false)
{
    const fs::path output_path(output_name);
    const std::string file_name = output_path.filename().string()
        + (corr ? ".corr" : "") + ".json";
    return (output_path.parent_path() / file_name).string();
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
                name.size() >= 5 && name.compare(name.size() - 5, 5, ".root") == 0)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}


bool inHyCal(float xmm, float ymm) {
    const float module = 20.75; // mm
    return (fabs(xmm) > module * 2.0 || fabs(ymm) > module * 2.0)
        && (fabs(xmm) < module * 16. && fabs(ymm) < module * 16.);
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
    std::string output_name;
    int  max_events  = -1;
    int  num_threads = 4;
    int  num_files   = -1;
    bool worker_mode = false;
    bool corr = false;

    static option long_options[] = {
        {"corr", no_argument, nullptr, 'c'},
        {nullptr, 0, nullptr, 0}
    };
    int opt;
    while ((opt = getopt_long_only(argc, argv, "o:n:f:j:",
                                   long_options, nullptr)) != -1) {
        switch (opt) {
            case 'o': output_name = optarg; break;
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
        std::cerr << "Usage: hycal_Erecon_check <input_recon.root|dir> [more...] "
                 "-o <output_name> [-n max_events] [-f nfiles] [-j threads] [-corr]\n";
        return 1;
    }

    if (output_name.empty()) {
        std::cerr << "No output prefix provided. Please pass -o <output_prefix>.\n";
        return 1;
    }
    ROOT::EnableThreadSafety();
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TH1F");
    TClass::GetClass("TH2F");

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);

    std::vector<long long> file_limits(root_files.size(), -1);
    if (max_events >= 0) {
        long long remaining = max_events;
        for (size_t i = 0; i < root_files.size(); ++i) {
            TFile input_file(root_files[i].c_str(), "READ");
            auto *tree = input_file.IsOpen()
                ? dynamic_cast<TTree *>(input_file.Get("recon")) : nullptr;
            const long long entries = tree ? tree->GetEntries() : 0;
            file_limits[i] = std::max(0LL, std::min(entries, remaining));
            remaining -= file_limits[i];
        }
    }

    auto merged = makeHistResult("");
    SharedFillLocks fill_locks;
    const int threads_count = std::max(1, std::min(num_threads,
        static_cast<int>(root_files.size())));
    const int rounds = (static_cast<int>(root_files.size()) + threads_count - 1)
        / threads_count;
    std::mutex io_mutex;
    std::cout << "Processing " << root_files.size() << " file(s) with "
              << threads_count << " thread(s), " << rounds << " round(s)\n";

    for (int round = 0; round < rounds; ++round) {
        const int first = round * threads_count;
        const int last = std::min(first + threads_count,
                                  static_cast<int>(root_files.size()));
        std::vector<std::thread> workers;
        workers.reserve(last - first);

        for (int file_index = first; file_index < last; ++file_index) {
            workers.emplace_back([&, file_index, first]() {
                const long long limit = max_events >= 0 ? file_limits[file_index] : -1;
                const bool ok = processRootFile(root_files[file_index], gRunConfig,
                                                db_dir, limit, merged.get(), &fill_locks);
                std::lock_guard<std::mutex> lock(io_mutex);
                std::cout << "[worker " << (file_index - first) << "] file "
                          << file_index << " / " << (root_files.size() - 1)
                          << ": " << root_files[file_index] << " -> "
                          << (ok ? "OK" : "FAILED") << "\n";
            });
        }
        for (auto &worker : workers) worker.join();
    }
    merged->events_processed = fill_locks.events_processed.load();

    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");
    analysis::PhysicsTools physics(hycal);

    std::array<std::unique_ptr<TH2F>, module_count> bowl_shapes;
    std::array<double, module_count> expected_energies{};
    for (int m = 0; m < module_count; ++m) {
        const int module_number = module_numbers[m];
        bowl_shapes[m] = std::make_unique<TH2F>(
            Form("h2_bowl_shape_W%d", module_number),
            Form("Bowl Shape W%d;X;Y;E_{recon}/E_{expect}", module_number),
            grids, 0, grids, grids, 0, grids);
        const auto *module = hycal.module_by_id(1000 + module_number);
        if (!module) {
            std::cerr << "Cannot find HyCal module W" << module_number << "\n";
            return 1;
        }
        const double angle = std::atan2(
            std::sqrt(module->x * module->x + module->y * module->y),
            gRunConfig.hycal_z);
        expected_energies[m] = analysis::PhysicsTools::ExpectedEnergy(
            angle, gRunConfig.Ebeam, "ep");

        for (int i = 0; i < grids; ++i) {
            for (int j = 0; j < grids; ++j) {
                auto fit = physics.fitPeak(
                    merged->h1_energy_grid[m][i][j].get(),
                    static_cast<float>(expected_energies[m]), true);
                // If the fit failed, use the mean of the histogram as the energy.
                // or the entries are too few to perform a reliable fit.
                if (merged->h1_energy_grid[m][i][j]->GetEntries() < 200) fit[0] = 0;
                const double energy = fit[0] != 0
                    ? fit[0] : merged->h1_energy_grid[m][i][j]->GetMean();
                bowl_shapes[m]->SetBinContent(
                    i + 1, j + 1, energy / expected_energies[m]);
            }
        }
    }

    TCanvas bowl_modules("bowl_modules", "HyCal Module Bowl Shapes", 3200, 3200);
    bowl_modules.Divide(8, 8, 0.002, 0.002);
    for (int m = 0; m < module_count; ++m) {
        bowl_modules.cd(module_canvas_pads[m]);
        gPad->SetLeftMargin(0.12);
        gPad->SetRightMargin(0.14);
        gPad->SetBottomMargin(0.12);
        bowl_shapes[m]->SetStats(0);
        bowl_shapes[m]->SetMinimum(0.98);
        bowl_shapes[m]->SetMaximum(1.02);
        bowl_shapes[m]->GetXaxis()->SetTitleSize(0.06);
        bowl_shapes[m]->GetYaxis()->SetTitleSize(0.06);
        bowl_shapes[m]->GetZaxis()->SetTitleSize(0.05);
        bowl_shapes[m]->GetXaxis()->SetLabelSize(0.05);
        bowl_shapes[m]->GetYaxis()->SetLabelSize(0.05);
        bowl_shapes[m]->GetZaxis()->SetLabelSize(0.045);
        bowl_shapes[m]->Draw("COLZ");
    }
    bowl_modules.Update();

    // fit the grids histograms for each module
    float energy_bias[1156][grids][grids] = {{{0}}};
    int fit_count = 0, mean_count = 0, analyzed_module_count = 0;
    for (int m = 0; m < 1156; ++m) {
        const auto *module = hycal.module_by_id(1001 + m);
        if (!module) continue;
        if (std::fabs(module->x) < 20.75 * 2.0 && std::fabs(module->y) < 20.75 * 2.0) continue;
        if (std::fabs(module->x) > 20.75 * 16.0 || std::fabs(module->y) > 20.75 * 16.0) continue;
        ++analyzed_module_count;
        const double angle = std::atan2(std::sqrt(module->x * module->x + module->y * module->y), gRunConfig.hycal_z);
        const double expected_energy = analysis::PhysicsTools::ExpectedEnergy(angle, gRunConfig.Ebeam, "ep");
        for (int i = 0; i < grids; ++i) {
            for (int j = 0; j < grids; ++j) {
                if (merged->h1_energy_grid_allModule[m][i][j]->GetEntries() < 200) continue;
                auto fit = physics.fitPeak(
                    merged->h1_energy_grid_allModule[m][i][j].get(),
                    static_cast<float>(expected_energy), true);
                if (merged->h1_energy_grid_allModule[m][i][j]->GetEntries() < 400) fit[0] = 0;
                const double energy = fit[0] != 0
                    ? (fit_count++, fit[0]) : (mean_count++, merged->h1_energy_grid_allModule[m][i][j]->GetMean());
                energy_bias[m][i][j] = std::clamp(energy / expected_energy - 1.0, -0.03, 0.03);
            }
        }
    }
    std::cout << "Fit count: " << fit_count << ", Mean count: " << mean_count << ", Module count: " << analyzed_module_count << ", Grids: " << analyzed_module_count * grids * grids << "\n";

    // output energy bias for each module and grid to a json file to put in database
    // structure of the JSON file:
    // {
    //     "W555": {
    //         "y4": { "x0": bias, "x1": bias, "x2": bias, "x3": bias, "x4": bias}
    //         "y3": { "x0": bias, "x1": bias, "x2": bias, "x3": bias, "x4": bias}
    //         "y2": { "x0": bias, "x1": bias, "x2": bias, "x3": bias, "x4": bias}
    //         "y1": { "x0": bias, "x1": bias, "x2": bias, "x3": bias, "x4": bias}
    //         "y0": { "x0": bias, "x1": bias, "x2": bias, "x3": bias, "x4": bias}
    //     }
    // }
    const std::string output_json_name = outputJsonFileName(output_name, corr);
    std::ofstream json_output(output_json_name);
    if (!json_output) {
        std::cerr << "Cannot create output file " << output_json_name << "\n";
        return 1;
    }
    json_output << std::fixed << std::setprecision(6) << "{\n";
    bool first_module = true;
    for (int m = 0; m < 1156; ++m) {
        if (!hycal.module_by_id(1001 + m)) continue;
        if (!first_module) json_output << ",\n";
        first_module = false;
        json_output << "  \"W" << m + 1 << "\": {\n";
        for (int row = grids - 1; row >= 0; --row) {
            json_output << "    \"y" << row << "\": {";
            for (int col = 0; col < grids; ++col) {
                if (col > 0) json_output << ",";
                json_output << " \"x" << col << "\": "
                            << energy_bias[m][col][row];
            }
            json_output << " }" << (row > 0 ? "," : "") << "\n";
        }
        json_output << "  }";
    }
    json_output << "\n}\n";
    json_output.close();
    std::cout << "Wrote energy-bias JSON to " << output_json_name << "\n";

    // TH1::Fit draws into the active pad, which is pad 64 after the loop above.
    // Restore W701 after all fits before persisting the combined canvas.
    bowl_modules.cd(module_canvas_pads.back());
    gPad->Clear();
    bowl_shapes.back()->Draw("COLZ");
    bowl_modules.Modified();
    bowl_modules.Update();

    const std::string output_file_name = outputFileName(output_name, corr);
    TFile output_file(output_file_name.c_str(), "RECREATE");
    if (output_file.IsZombie()) {
        std::cerr << "Cannot create output file " << output_file_name << "\n";
        return 1;
    }
    merged->h2_hit_module_hycal->Write();
    merged->h2_hit_module_gem->Write();

    for (int m = 0; m < module_count; ++m) {
        const std::string directory = Form("energy_grids_W%d", module_numbers[m]);
        output_file.mkdir(directory.c_str());
        output_file.cd(directory.c_str());
        for (int i = 0; i < grids; ++i) {
            for (int j = 0; j < grids; ++j) {
                merged->h1_energy_grid[m][i][j]->Write();
            }
        }
    }
    output_file.cd();
    for (const auto &bowl_shape : bowl_shapes) bowl_shape->Write();
    bowl_modules.Write("bowl_shapes_8x8");
    output_file.mkdir("energy_grids_allModule");
    output_file.cd("energy_grids_allModule");
    for (int m = 0; m < 1156; ++m) {
        for (int i = 0; i < grids; ++i) {
            for (int j = 0; j < grids; ++j) {
                merged->h1_energy_grid_allModule[m][i][j]->Write();
            }
        }
    }
    output_file.Close();

}

static std::unique_ptr<HistResult> makeHistResult(const std::string &suffix)
{
    auto result = std::make_unique<HistResult>();
    const std::string name_suffix = suffix.empty() ? "" : "_" + suffix;
    result->h2_hit_module_hycal = std::make_unique<TH2F>(
        Form("h2_hit_module_hycal%s", name_suffix.c_str()),
        "HyCal Hit Distribution (Module);(X_{hycal}-X_{cell center})/d_{cell size};(Y_{hycal}-Y_{cell center})/d_{cell size}",
        100, -0.5, 0.5, 100, -0.5, 0.5);
    result->h2_hit_module_gem = std::make_unique<TH2F>(
        Form("h2_hit_module_gem%s", name_suffix.c_str()),
        "GEM Hit Distribution (Module);(X_{gem}-X_{cell center})/d_{cell size};(Y_{gem}-Y_{cell center})/d_{cell size}",
        100, -0.5, 0.5, 100, -0.5, 0.5);
    for (int m = 0; m < module_count; ++m) {
        for (int i = 0; i < grids; ++i) {
            for (int j = 0; j < grids; ++j) {
                result->h1_energy_grid[m][i][j] = std::make_unique<TH1F>(
                    Form("h1_energy_grid_W%d_%d_%d%s",
                         module_numbers[m], i, j, name_suffix.c_str()),
                    Form("Energy Grid W%d;E_{recon} [MeV];Counts",
                         module_numbers[m]),
                    energy_bins, energy_min, energy_max);
            }
        }
    }
    for (int m = 0; m < 1156; ++m) {
        for (int i = 0; i < grids; ++i) {
            for (int j = 0; j < grids; ++j) {
                result->h1_energy_grid_allModule[m][i][j] = std::make_unique<TH1F>(
                    Form("h1_energy_grid_allModule_W%d_%d_%d%s",
                        m + 1, i, j, name_suffix.c_str()),
                    Form("Energy Grid All Module W%d;E_{recon} [MeV];Counts",
                        m + 1),
                    energy_bins, energy_min, energy_max);
            }
        }
    }
    return result;
}

static bool processRootFile(const std::string &input_file, const RunConfig &run_config,
                            const std::string &db_dir, long long max_events,
                            HistResult *result, SharedFillLocks *fill_locks)
{
    TFile input(input_file.c_str(), "READ");
    if (input.IsZombie() || !input.IsOpen()) {
        std::cerr << "Cannot open input file " << input_file << "\n";
        return false;
    }
    auto *tree = dynamic_cast<TTree *>(input.Get("recon"));
    if (!tree) {
        std::cerr << "Input file has no 'recon' tree: " << input_file << "\n";
        return false;
    }

    EventVars_Recon event;
    prad2::SetReconReadBranches(tree, event);
    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");

    Long64_t entries = tree->GetEntries();
    if (max_events >= 0) entries = std::min(entries, max_events);
    long long events_processed = 0;
    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        if ((event.trigger_bits & prad2::TBIT_sum) == 0) continue;
        if (event.n_clusters != 1) continue;
        //if (event.matchNum != 1) continue;
        if (event.cl_nblocks[0] < 3) continue;
        if (fabs(event.cl_energy[0] - run_config.Ebeam) > 3.0 * 0.03 * std::sqrt(run_config.Ebeam * 1000.)) continue;

        HCHit hc_hit;
        GEMHit gem_hit;
        hc_hit.x = event.cl_x[0];
        hc_hit.y = event.cl_y[0];
        hc_hit.z = event.cl_z[0];
        hc_hit.energy = event.cl_energy[0];
        if (event.matchNum == 1){
            gem_hit.x = event.mHit_gx[0][0];
            gem_hit.y = event.mHit_gy[0][0];
            gem_hit.z = event.mHit_gz[0][0];
        }

        if (gem_hit.z != 0.f) {
            const float scale = hc_hit.z / gem_hit.z;
            gem_hit.x *= scale;
            gem_hit.y *= scale;
            gem_hit.z *= scale;
        }
        ApplyToHyCal(gem_hit, run_config);
        ApplyToHyCal(hc_hit, run_config);
        if (!inHyCal(hc_hit.x, hc_hit.y)) continue;

        const auto *mod = hycal.module_by_id(event.cl_center[0]);
        if (!mod) continue;
        float xd_hycal = (hc_hit.x - mod->x) / mod->size_x;
        float yd_hycal = (hc_hit.y - mod->y) / mod->size_y;

        float xd_gem = (gem_hit.x - mod->x) / mod->size_x;
        float yd_gem = (gem_hit.y - mod->y) / mod->size_y;

        // Fill the selected module's energy-grid histogram.
        int col = static_cast<int>((xd_hycal + 0.5f) * grids);
        int row = static_cast<int>((yd_hycal + 0.5f) * grids);
        if (col < 0) col = 0;
        if (col >= grids) col = grids - 1;
        if (row < 0) row = 0;
        if (row >= grids) row = grids - 1;
        const int all_module_index = mod->id - 1001;
        {
            std::lock_guard<std::mutex> lock(fill_locks->module[all_module_index]);
            const auto module_it = std::find(
                module_numbers.begin(), module_numbers.end(), mod->id - 1000);
            if (module_it != module_numbers.end()) {
                const int module_index = std::distance(module_numbers.begin(), module_it);
                result->h1_energy_grid[module_index][col][row]->Fill(hc_hit.energy);
            }
            result->h1_energy_grid_allModule[all_module_index][col][row]->Fill(hc_hit.energy);
        }

        {
            std::lock_guard<std::mutex> lock(fill_locks->hit_maps);
            result->h2_hit_module_hycal->Fill(xd_hycal, yd_hycal);
            if (event.matchNum == 1) result->h2_hit_module_gem->Fill(xd_gem, yd_gem);
        }
        ++events_processed;
    }
    fill_locks->events_processed.fetch_add(events_processed, std::memory_order_relaxed);
    return true;
}