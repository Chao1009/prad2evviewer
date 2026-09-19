// energy_corr.cpp : to correct the reconstructed energy non-uniformity depending on
// the position of the cluster within the calorimeter modules.
//
// Measure the reconstructed-energy response at different positions within each
// module and compare its fitted peak with the expected energy. Each module is divided into a grid,
// First try to a grid of 5 by 5, the map is below, could be saved into a 2D array of 1D hist
// with one reconstructed-energy histogram for each cell in the 5x5 grid.
// The hists would be created for each cell in the 5x5 grid
// for example hist[5][5], the first index is the row and the second index is the column.
// column   0   1   2   3   4
// row     +---+---+---+---+---+     beam top  ^
//  0      |   |   |   |   |   |               |
//  1      |   |   |   |   |   |
//  2      |   |   |   |   |   |     beam right ->
//  3      |   |   |   |   |   |
//  4      |   |   |   |   |   |
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

#include <iostream>
#include <array>
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
const int energy_bins = 500; const double energy_min = 0., energy_max = 5000.;
const int grids = 5;

struct HistResult {
    std::unique_ptr<TH2F> h2_hit_module_hycal;
    std::unique_ptr<TH2F> h2_hit_module_gem;
    std::unique_ptr<TH1F> h1_energy_grid_W521[grids][grids];
    std::unique_ptr<TH1F> h1_energy_grid_W522[grids][grids];
    std::unique_ptr<TH1F> h1_energy_grid_W523[grids][grids];
    std::unique_ptr<TH1F> h1_energy_grid_W633[grids][grids];
    std::unique_ptr<TH1F> h1_energy_grid_W634[grids][grids];
    std::unique_ptr<TH1F> h1_energy_grid_W635[grids][grids];
    std::vector<std::unique_ptr<TH1F>> h1_energy_grid[grids][grids];
    long long events_processed = 0;
};

static std::unique_ptr<HistResult> makeHistResult(const std::string &suffix);
static bool processRootFile(const std::string &input_file, const RunConfig &run_config,
                            const std::string &db_dir, long long max_events,
                            HistResult *result);

static void mergeHistResult(HistResult &destination, const HistResult &source)
{
    destination.h2_hit_module_hycal->Add(source.h2_hit_module_hycal.get());
    destination.h2_hit_module_gem->Add(source.h2_hit_module_gem.get());
    for (int i = 0; i < grids; ++i) {
        for (int j = 0; j < grids; ++j) {
            destination.h1_energy_grid_W521[i][j]->Add(source.h1_energy_grid_W521[i][j].get());
            destination.h1_energy_grid_W522[i][j]->Add(source.h1_energy_grid_W522[i][j].get());
            destination.h1_energy_grid_W523[i][j]->Add(source.h1_energy_grid_W523[i][j].get());
            destination.h1_energy_grid_W633[i][j]->Add(source.h1_energy_grid_W633[i][j].get());
            destination.h1_energy_grid_W634[i][j]->Add(source.h1_energy_grid_W634[i][j].get());
            destination.h1_energy_grid_W635[i][j]->Add(source.h1_energy_grid_W635[i][j].get());
        }
    }
    for (int m = 0; m < 1156; ++m) {
        for (int i = 0; i < grids; ++i) {
            for (int j = 0; j < grids; ++j) {
                destination.h1_energy_grid[m][i][j]->Add(source.h1_energy_grid[m][i][j].get());
            }
        }   
    }
    destination.events_processed += source.events_processed;
}

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
        && (fabs(xmm) < module * 15. && fabs(ymm) < module * 15.);
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
        std::vector<std::unique_ptr<HistResult>> results(last - first);
        std::vector<std::thread> workers;
        workers.reserve(last - first);

        for (int file_index = first; file_index < last; ++file_index) {
            workers.emplace_back([&, file_index, first]() {
                auto result = makeHistResult(Form("worker%d", file_index - first));
                const long long limit = max_events >= 0 ? file_limits[file_index] : -1;
                const bool ok = processRootFile(root_files[file_index], gRunConfig,
                                                db_dir, limit, result.get());
                results[file_index - first] = std::move(result);
                std::lock_guard<std::mutex> lock(io_mutex);
                std::cout << "[worker " << (file_index - first) << "] file "
                          << file_index << " / " << (root_files.size() - 1)
                          << ": " << root_files[file_index] << " -> "
                          << (ok ? "OK" : "FAILED") << "\n";
            });
        }
        for (auto &worker : workers) worker.join();
        for (const auto &result : results) {
            if (result) mergeHistResult(*merged, *result);
        }
    }

    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");
    analysis::PhysicsTools physics(hycal);

    // Fit or get the mean energy for each grid cell, fill the bowl shape histograms
    TH2F *h2_bowl_shape_W521 = new TH2F("h2_bowl_shape_W521", "Bowl Shape W521;X;Y;E_{recon}/E_{expect}", grids, 0, grids, grids, 0, grids);
    TH2F *h2_bowl_shape_W522 = new TH2F("h2_bowl_shape_W522", "Bowl Shape W522;X;Y;E_{recon}/E_{expect}", grids, 0, grids, grids, 0, grids);
    TH2F *h2_bowl_shape_W523 = new TH2F("h2_bowl_shape_W523", "Bowl Shape W523;X;Y;E_{recon}/E_{expect}", grids, 0, grids, grids, 0, grids);
    TH2F *h2_bowl_shape_W633 = new TH2F("h2_bowl_shape_W633", "Bowl Shape W633;X;Y;E_{recon}/E_{expect}", grids, 0, grids, grids, 0, grids);
    TH2F *h2_bowl_shape_W634 = new TH2F("h2_bowl_shape_W634", "Bowl Shape W634;X;Y;E_{recon}/E_{expect}", grids, 0, grids, grids, 0, grids);
    TH2F *h2_bowl_shape_W635 = new TH2F("h2_bowl_shape_W635", "Bowl Shape W635;X;Y;E_{recon}/E_{expect}", grids, 0, grids, grids, 0, grids);

    double E_mean_W521[grids][grids], E_mean_W522[grids][grids], E_mean_W523[grids][grids];
    double E_mean_W633[grids][grids], E_mean_W634[grids][grids], E_mean_W635[grids][grids];
    memset(E_mean_W521, 0, sizeof(E_mean_W521));
    memset(E_mean_W522, 0, sizeof(E_mean_W522));
    memset(E_mean_W523, 0, sizeof(E_mean_W523));
    memset(E_mean_W633, 0, sizeof(E_mean_W633));
    memset(E_mean_W634, 0, sizeof(E_mean_W634));
    memset(E_mean_W635, 0, sizeof(E_mean_W635));
    const auto &mod_W521 = hycal.module_by_id(1521);
    double angle_W521 = std::atan2(std::sqrt(mod_W521->x * mod_W521->x + mod_W521->y * mod_W521->y), gRunConfig.hycal_z);
    double expected_energy_W521 = analysis::PhysicsTools::ExpectedEnergy(angle_W521, gRunConfig.Ebeam, "ep");
    const auto &mod_W522 = hycal.module_by_id(1522);
    double angle_W522 = std::atan2(std::sqrt(mod_W522->x * mod_W522->x + mod_W522->y * mod_W522->y), gRunConfig.hycal_z);
    double expected_energy_W522 = analysis::PhysicsTools::ExpectedEnergy(angle_W522, gRunConfig.Ebeam, "ep");
    const auto &mod_W523 = hycal.module_by_id(1523);
    double angle_W523 = std::atan2(std::sqrt(mod_W523->x * mod_W523->x + mod_W523->y * mod_W523->y), gRunConfig.hycal_z);
    double expected_energy_W523 = analysis::PhysicsTools::ExpectedEnergy(angle_W523, gRunConfig.Ebeam, "ep");
    const auto &mod_W633 = hycal.module_by_id(1633);
    double angle_W633 = std::atan2(std::sqrt(mod_W633->x * mod_W633->x + mod_W633->y * mod_W633->y), gRunConfig.hycal_z);
    double expected_energy_W633 = analysis::PhysicsTools::ExpectedEnergy(angle_W633, gRunConfig.Ebeam, "ep");
    const auto &mod_W634 = hycal.module_by_id(1634);
    double angle_W634 = std::atan2(std::sqrt(mod_W634->x * mod_W634->x + mod_W634->y * mod_W634->y), gRunConfig.hycal_z);
    double expected_energy_W634 = analysis::PhysicsTools::ExpectedEnergy(angle_W634, gRunConfig.Ebeam, "ep");
    const auto &mod_W635 = hycal.module_by_id(1635);
    double angle_W635 = std::atan2(std::sqrt(mod_W635->x * mod_W635->x + mod_W635->y * mod_W635->y), gRunConfig.hycal_z);
    double expected_energy_W635 = analysis::PhysicsTools::ExpectedEnergy(angle_W635, gRunConfig.Ebeam, "ep");
    for (int i = 0; i < grids; ++i) {
        for (int j = 0; j < grids; ++j) {
            const auto fit_W521 = physics.fitPeak(merged->h1_energy_grid_W521[i][j].get(), static_cast<float>(expected_energy_W521), true);
            if (fit_W521[0] != 0) E_mean_W521[i][j] = fit_W521[0]/expected_energy_W521;
            else E_mean_W521[i][j] = merged->h1_energy_grid_W521[i][j]->GetMean()/expected_energy_W521;
            const auto fit_W522 = physics.fitPeak(merged->h1_energy_grid_W522[i][j].get(), static_cast<float>(expected_energy_W522), true);
            if (fit_W522[0] != 0) E_mean_W522[i][j] = fit_W522[0]/expected_energy_W522;
            else E_mean_W522[i][j] = merged->h1_energy_grid_W522[i][j]->GetMean()/expected_energy_W522;
            const auto fit_W523 = physics.fitPeak(merged->h1_energy_grid_W523[i][j].get(), static_cast<float>(expected_energy_W523), true);
            if (fit_W523[0] != 0) E_mean_W523[i][j] = fit_W523[0]/expected_energy_W523;
            else E_mean_W523[i][j] = merged->h1_energy_grid_W523[i][j]->GetMean()/expected_energy_W523;
            const auto fit_W633 = physics.fitPeak(merged->h1_energy_grid_W633[i][j].get(), static_cast<float>(expected_energy_W633), true);
            if (fit_W633[0] != 0) E_mean_W633[i][j] = fit_W633[0]/expected_energy_W633;
            else E_mean_W633[i][j] = merged->h1_energy_grid_W633[i][j]->GetMean()/expected_energy_W633;
            const auto fit_W634 = physics.fitPeak(merged->h1_energy_grid_W634[i][j].get(), static_cast<float>(expected_energy_W634), true);
            if (fit_W634[0] != 0) E_mean_W634[i][j] = fit_W634[0]/expected_energy_W634;
            else E_mean_W634[i][j] = merged->h1_energy_grid_W634[i][j]->GetMean()/expected_energy_W634;
            const auto fit_W635 = physics.fitPeak(merged->h1_energy_grid_W635[i][j].get(), static_cast<float>(expected_energy_W635), true);
            if (fit_W635[0] != 0) E_mean_W635[i][j] = fit_W635[0]/expected_energy_W635;
            else E_mean_W635[i][j] = merged->h1_energy_grid_W635[i][j]->GetMean()/expected_energy_W635;
        }
    }

    for (int i = 0; i < grids; ++i) {
        for (int j = 0; j < grids; ++j) {
            h2_bowl_shape_W521->SetBinContent(i+1, j+1, E_mean_W521[i][j]);
            h2_bowl_shape_W522->SetBinContent(i+1, j+1, E_mean_W522[i][j]);
            h2_bowl_shape_W523->SetBinContent(i+1, j+1, E_mean_W523[i][j]);
            h2_bowl_shape_W633->SetBinContent(i+1, j+1, E_mean_W633[i][j]);
            h2_bowl_shape_W634->SetBinContent(i+1, j+1, E_mean_W634[i][j]);
            h2_bowl_shape_W635->SetBinContent(i+1, j+1, E_mean_W635[i][j]);
        }
    }
    TCanvas bowl_W521("bowl_W521", "Bowl Shape of W521", 800, 600);
    bowl_W521.SetLeftMargin(0.10);
    bowl_W521.SetBottomMargin(0.10);
    bowl_W521.SetRightMargin(0.05);
    gStyle->SetPalette(kRainBow);
    h2_bowl_shape_W521->SetStats(0);
    h2_bowl_shape_W521->GetXaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W521->GetYaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W521->GetZaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W521->GetXaxis()->SetTitleOffset(1.8);
    h2_bowl_shape_W521->GetYaxis()->SetTitleOffset(2.0);
    h2_bowl_shape_W521->GetZaxis()->SetTitleOffset(1.3);
    h2_bowl_shape_W521->GetXaxis()->CenterTitle();
    h2_bowl_shape_W521->GetYaxis()->CenterTitle();
    h2_bowl_shape_W521->GetZaxis()->CenterTitle();
    h2_bowl_shape_W521->SetMinimum(0.99);
    h2_bowl_shape_W521->SetMaximum(1.02);
    h2_bowl_shape_W521->Draw("LEGO2Z");  // 3D colored blocks with Z-palette
    gPad->SetTheta(30.);
    gPad->SetPhi(40.);
    gPad->Update();

    TCanvas bowl_W522("bowl_W522", "Bowl Shape of W522", 800, 600);
    bowl_W522.SetLeftMargin(0.10);
    bowl_W522.SetBottomMargin(0.10);
    bowl_W522.SetRightMargin(0.05);
    gStyle->SetPalette(kRainBow);
    h2_bowl_shape_W522->SetStats(0);
    h2_bowl_shape_W522->GetXaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W522->GetYaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W522->GetZaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W522->GetXaxis()->SetTitleOffset(1.8);
    h2_bowl_shape_W522->GetYaxis()->SetTitleOffset(2.0);
    h2_bowl_shape_W522->GetZaxis()->SetTitleOffset(1.3);
    h2_bowl_shape_W522->GetXaxis()->CenterTitle();
    h2_bowl_shape_W522->GetYaxis()->CenterTitle();
    h2_bowl_shape_W522->GetZaxis()->CenterTitle();
    h2_bowl_shape_W522->SetMinimum(0.99);
    h2_bowl_shape_W522->SetMaximum(1.02);
    h2_bowl_shape_W522->Draw("LEGO2Z");  // 3D colored blocks with Z-palette
    gPad->SetTheta(30.);
    gPad->SetPhi(40.);
    gPad->Update();

    TCanvas bowl_W523("bowl_W523", "Bowl Shape of W523", 800, 600);
    bowl_W523.SetLeftMargin(0.10);
    bowl_W523.SetBottomMargin(0.10);
    bowl_W523.SetRightMargin(0.05);
    gStyle->SetPalette(kRainBow);
    h2_bowl_shape_W523->SetStats(0);
    h2_bowl_shape_W523->GetXaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W523->GetYaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W523->GetZaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W523->GetXaxis()->SetTitleOffset(1.8);
    h2_bowl_shape_W523->GetYaxis()->SetTitleOffset(2.0);
    h2_bowl_shape_W523->GetZaxis()->SetTitleOffset(1.3);
    h2_bowl_shape_W523->GetXaxis()->CenterTitle();
    h2_bowl_shape_W523->GetYaxis()->CenterTitle();
    h2_bowl_shape_W523->GetZaxis()->CenterTitle();
    h2_bowl_shape_W523->SetMinimum(0.99);
    h2_bowl_shape_W523->SetMaximum(1.02);
    h2_bowl_shape_W523->Draw("LEGO2Z");  // 3D colored blocks with Z-palette
    gPad->SetTheta(30.);
    gPad->SetPhi(40.);
    gPad->Update();

    TCanvas bowl_W633("bowl_W633", "Bowl Shape of W633", 800, 600);
    bowl_W633.SetLeftMargin(0.10);
    bowl_W633.SetBottomMargin(0.10);
    bowl_W633.SetRightMargin(0.05);
    gStyle->SetPalette(kRainBow);
    h2_bowl_shape_W633->SetStats(0);
    h2_bowl_shape_W633->GetXaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W633->GetYaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W633->GetZaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W633->GetXaxis()->SetTitleOffset(1.8);
    h2_bowl_shape_W633->GetYaxis()->SetTitleOffset(2.0);
    h2_bowl_shape_W633->GetZaxis()->SetTitleOffset(1.3);
    h2_bowl_shape_W633->GetXaxis()->CenterTitle();
    h2_bowl_shape_W633->GetYaxis()->CenterTitle();
    h2_bowl_shape_W633->GetZaxis()->CenterTitle();
    h2_bowl_shape_W633->SetMinimum(0.99);
    h2_bowl_shape_W633->SetMaximum(1.02);
    h2_bowl_shape_W633->Draw("LEGO2Z");  // 3D colored blocks with Z-palette
    gPad->SetTheta(30.);
    gPad->SetPhi(40.);
    gPad->Update();

    TCanvas bowl_W634("bowl_W634", "Bowl Shape of W634", 800, 600);
    bowl_W634.SetLeftMargin(0.10);
    bowl_W634.SetBottomMargin(0.10);
    bowl_W634.SetRightMargin(0.05);
    gStyle->SetPalette(kRainBow);
    h2_bowl_shape_W634->SetStats(0);
    h2_bowl_shape_W634->GetXaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W634->GetYaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W634->GetZaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W634->GetXaxis()->SetTitleOffset(1.8);
    h2_bowl_shape_W634->GetYaxis()->SetTitleOffset(2.0);
    h2_bowl_shape_W634->GetZaxis()->SetTitleOffset(1.3);
    h2_bowl_shape_W634->GetXaxis()->CenterTitle();
    h2_bowl_shape_W634->GetYaxis()->CenterTitle();
    h2_bowl_shape_W634->GetZaxis()->CenterTitle();
    h2_bowl_shape_W634->SetMinimum(0.99);
    h2_bowl_shape_W634->SetMaximum(1.02);
    h2_bowl_shape_W634->Draw("LEGO2Z");  // 3D colored blocks with Z-palette
    gPad->SetTheta(30.);
    gPad->SetPhi(40.);
    gPad->Update();

    TCanvas bowl_W635("bowl_W635", "Bowl Shape of W635", 800, 600);
    bowl_W635.SetLeftMargin(0.10);
    bowl_W635.SetBottomMargin(0.10);
    bowl_W635.SetRightMargin(0.05);
    gStyle->SetPalette(kRainBow);
    h2_bowl_shape_W635->SetStats(0);
    h2_bowl_shape_W635->GetXaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W635->GetYaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W635->GetZaxis()->SetTitleSize(0.04);
    h2_bowl_shape_W635->GetXaxis()->SetTitleOffset(1.8);
    h2_bowl_shape_W635->GetYaxis()->SetTitleOffset(2.0);
    h2_bowl_shape_W635->GetZaxis()->SetTitleOffset(1.3);
    h2_bowl_shape_W635->GetXaxis()->CenterTitle();
    h2_bowl_shape_W635->GetYaxis()->CenterTitle();
    h2_bowl_shape_W635->GetZaxis()->CenterTitle();
    h2_bowl_shape_W635->SetMinimum(0.99);
    h2_bowl_shape_W635->SetMaximum(1.02);
    h2_bowl_shape_W635->Draw("LEGO2Z");  // 3D colored blocks with Z-palette
    gPad->SetTheta(30.);
    gPad->SetPhi(40.);
    gPad->Update();

    const std::string output_file_name = outputFileName(output_name, corr);
    TFile output_file(output_file_name.c_str(), "RECREATE");
    if (output_file.IsZombie()) {
        std::cerr << "Cannot create output file " << output_file_name << "\n";
        return 1;
    }
    merged->h2_hit_module_hycal->Write();
    merged->h2_hit_module_gem->Write();
    output_file.cd();
    output_file.mkdir("energy_grids_W521");
    output_file.mkdir("energy_grids_W522");
    output_file.mkdir("energy_grids_W523");
    output_file.mkdir("energy_grids_W633");
    output_file.mkdir("energy_grids_W634");
    output_file.mkdir("energy_grids_W635");
    for (int i = 0; i < grids; ++i) {
        for (int j = 0; j < grids; ++j) {
            output_file.cd("energy_grids_W521");
            merged->h1_energy_grid_W521[i][j]->Write();
            output_file.cd("energy_grids_W522");
            merged->h1_energy_grid_W522[i][j]->Write();
            output_file.cd("energy_grids_W523");
            merged->h1_energy_grid_W523[i][j]->Write();
            output_file.cd("energy_grids_W633");
            merged->h1_energy_grid_W633[i][j]->Write();
            output_file.cd("energy_grids_W634");
            merged->h1_energy_grid_W634[i][j]->Write();
            output_file.cd("energy_grids_W635");
            merged->h1_energy_grid_W635[i][j]->Write();
        }
    }
    output_file.cd();
    h2_bowl_shape_W521->Write("h2_bowl_shape_W521");
    h2_bowl_shape_W522->Write("h2_bowl_shape_W522");
    h2_bowl_shape_W523->Write("h2_bowl_shape_W523");
    h2_bowl_shape_W633->Write("h2_bowl_shape_W633");
    h2_bowl_shape_W634->Write("h2_bowl_shape_W634");
    h2_bowl_shape_W635->Write("h2_bowl_shape_W635");
    bowl_W521.Write("bowl_shape_W521");
    bowl_W522.Write("bowl_shape_W522");
    bowl_W523.Write("bowl_shape_W523");
    bowl_W633.Write("bowl_shape_W633");
    bowl_W634.Write("bowl_shape_W634");
    bowl_W635.Write("bowl_shape_W635");
    output_file.cd();
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
    for (int i = 0; i < grids; ++i) {
        for (int j = 0; j < grids; ++j) {
                    result->h1_energy_grid_W521[i][j] = std::make_unique<TH1F>(
                Form("h1_energy_grid_W521_%d_%d%s", i, j, name_suffix.c_str()),
                "Energy Grid W521;E_{recon} [MeV];Counts",
                energy_bins, energy_min, energy_max);
            result->h1_energy_grid_W522[i][j] = std::make_unique<TH1F>(
                Form("h1_energy_grid_W522_%d_%d%s", i, j, name_suffix.c_str()),
                "Energy Grid W522;E_{recon} [MeV];Counts",
                energy_bins, energy_min, energy_max);
            result->h1_energy_grid_W523[i][j] = std::make_unique<TH1F>(
                Form("h1_energy_grid_W523_%d_%d%s", i, j, name_suffix.c_str()),
                "Energy Grid W523;E_{recon} [MeV];Counts",
                energy_bins, energy_min, energy_max);
            result->h1_energy_grid_W633[i][j] = std::make_unique<TH1F>(
                Form("h1_energy_grid_W633_%d_%d%s", i, j, name_suffix.c_str()),
                "Energy Grid W633;E_{recon} [MeV];Counts",
                energy_bins, energy_min, energy_max);
            result->h1_energy_grid_W634[i][j] = std::make_unique<TH1F>(
                Form("h1_energy_grid_W634_%d_%d%s", i, j, name_suffix.c_str()),
                "Energy Grid W634;E_{recon} [MeV];Counts",
                energy_bins, energy_min, energy_max);
            result->h1_energy_grid_W635[i][j] = std::make_unique<TH1F>(
                Form("h1_energy_grid_W635_%d_%d%s", i, j, name_suffix.c_str()),
                "Energy Grid W635;E_{recon} [MeV];Counts",
                energy_bins, energy_min, energy_max);
        }
    }
    return result;
}

static bool processRootFile(const std::string &input_file, const RunConfig &run_config,
                            const std::string &db_dir, long long max_events,
                            HistResult *result)
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
    for (Long64_t entry = 0; entry < entries; ++entry) {
        tree->GetEntry(entry);
        if ((event.trigger_bits & prad2::TBIT_sum) == 0) continue;
        if (event.n_clusters != 1) continue;
        if (event.matchNum != 1) continue;
        if (event.cl_nblocks[0] < 3) continue;
        if (fabs(event.cl_energy[0] - run_config.Ebeam) > 3.0 * 0.03 * std::sqrt(run_config.Ebeam * 1000.)) continue;

        HCHit hc_hit;
        GEMHit gem_hit;
        hc_hit.x = event.cl_x[0];
        hc_hit.y = event.cl_y[0];
        hc_hit.z = event.cl_z[0];
        hc_hit.energy = event.cl_energy[0];
        gem_hit.x = event.mHit_gx[0][0];
        gem_hit.y = event.mHit_gy[0][0];
        gem_hit.z = event.mHit_gz[0][0];

        if (gem_hit.z == 0.f) continue;
        const float scale = hc_hit.z / gem_hit.z;
        gem_hit.x *= scale;
        gem_hit.y *= scale;
        gem_hit.z *= scale;
        ApplyToHyCal(gem_hit, run_config);
        ApplyToHyCal(hc_hit, run_config);
        if (!inHyCal(hc_hit.x, hc_hit.y)) continue;

        const auto *mod = hycal.module_by_id(event.cl_center[0]);
        if (!mod) continue;
        const float dx = hc_hit.x - gem_hit.x;
        const float dy = hc_hit.y - gem_hit.y;
        float xd_hycal = (hc_hit.x - mod->x) / mod->size_x;
        float yd_hycal = (hc_hit.y - mod->y) / mod->size_y;
        //if (xd_hycal < -0.5f) xd_hycal += 1.0f;
        //if (xd_hycal >  0.5f) xd_hycal -= 1.0f;
        //if (yd_hycal < -0.5f) yd_hycal += 1.0f;
        //if (yd_hycal >  0.5f) yd_hycal -= 1.0f;
        float xd_gem = (gem_hit.x - mod->x) / mod->size_x;
        float yd_gem = (gem_hit.y - mod->y) / mod->size_y;
        //if (xd_gem < -0.5f) xd_gem += 1.0f;
        //if (xd_gem >  0.5f) xd_gem -= 1.0f;
        //if (yd_gem < -0.5f) yd_gem += 1.0f;
        //if (yd_gem >  0.5f) yd_gem -= 1.0f;

        // Fill the energy grid histograms for the W521-W526 modules
        int col = static_cast<int>((xd_hycal + 0.5f) * grids);
        int row = static_cast<int>((yd_hycal + 0.5f) * grids);
        if (col < 0) col = 0;
        if (col >= grids) col = grids - 1;
        if (row < 0) row = 0;
        if (row >= grids) row = grids - 1;
        if (mod->id == 1521) result->h1_energy_grid_W521[col][row]->Fill(hc_hit.energy);
        if (mod->id == 1522) result->h1_energy_grid_W522[col][row]->Fill(hc_hit.energy);
        if (mod->id == 1523) result->h1_energy_grid_W523[col][row]->Fill(hc_hit.energy);
        if (mod->id == 1633) result->h1_energy_grid_W633[col][row]->Fill(hc_hit.energy);
        if (mod->id == 1634) result->h1_energy_grid_W634[col][row]->Fill(hc_hit.energy);
        if (mod->id == 1635) result->h1_energy_grid_W635[col][row]->Fill(hc_hit.energy);

        result->h2_hit_module_hycal->Fill(xd_hycal, yd_hycal);
        result->h2_hit_module_gem->Fill(xd_gem, yd_gem);
        ++result->events_processed;
    }
    return true;
}