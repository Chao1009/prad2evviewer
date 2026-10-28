// hycal_energy_bowl_shape.cpp : to correct the reconstructed energy non-uniformity depending on
// the position of the cluster within the calorimeter modules.
//
// Measure the reconstructed-energy response at different positions within each
// module and compare its fitted peak with the expected energy. Each module is divided into
// a 5x5 grid with one reconstructed-energy histogram per cell (map below).
// The histograms use h1_energy_grid_allModule[module][column][row]. Columns increase
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
#include "HyCalEnergyBias.h"
#include "MatchingTools.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"
#include "ToolUtils.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TCanvas.h>

#include <iostream>
#include <array>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <atomic>
#include <memory>
#include <mutex>
#include <getopt.h>

using namespace analysis;

using EventVars_Recon = prad2::ReconEventData;

const int energy_bins = 350; const double energy_min = 500., energy_max = 4000.;
const int grids = fdec::HyCalEnergyBias::GRID_SIZE;

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

int main(int argc, char *argv[])
{
    std::string db_dir = prad2::database_dir();

    // ── Argument parsing ─────────────────────────────────────────────────────
    std::string output_name;
    int  max_events  = -1;
    int  num_threads = 4;
    int  num_files   = -1;

    // getopt_long_only reports an unknown word option such as '-xyz' as one
    // unrecognized option instead of splitting it into short options.
    static option long_options[] = {
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

    std::vector<std::string> root_files =
        CollectInputs(argc, argv, optind, IsReconRootName, num_files);
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: hycal_energy_bowl_shape <input_recon.root|dir> [more...] "
                 "-o <output_name> [-n max_events] [-f nfiles] [-j threads]\n";
        return 1;
    }

    if (output_name.empty()) {
        std::cerr << "No output prefix provided. Please pass -o <output_prefix>.\n";
        return 1;
    }
    InitRootThreading();

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);

    const auto file_limits = DistributeEventBudget(root_files, "recon", max_events);

    auto merged = makeHistResult("");
    SharedFillLocks fill_locks;
    RunFilesInRounds(root_files, num_threads, [&](int idx, int) {
        return processRootFile(root_files[idx], gRunConfig, db_dir, file_limits[idx],
                               merged.get(), &fill_locks);
    });
    merged->events_processed = fill_locks.events_processed.load();

    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");

    // Copies of the drawn modules' grids for energy_grids_W<N>/, so the
    // display fits below attach their functions to the copies only.
    std::unique_ptr<TH1F> energy_grids[module_count][grids][grids];
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
                auto &grid = energy_grids[m][i][j];
                grid.reset(static_cast<TH1F *>(
                    merged->h1_energy_grid_allModule[module_number - 1][i][j]->Clone(
                        Form("h1_energy_grid_W%d_%d_%d", module_number, i, j))));
                grid->SetTitle(Form("Energy Grid W%d;E_{recon} [MeV];Counts", module_number));
                auto fit = analysis::PhysicsTools::fitPeak(
                    grid.get(), static_cast<float>(expected_energies[m]), true);
                // Use the histogram mean when the fit failed or there are too few entries.
                if (grid->GetEntries() < 200) fit[0] = 0;
                const double energy = fit[0] != 0 ? fit[0] : grid->GetMean();
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
        if (!InHyCalRing(module->x, module->y, 2.0, 16.)) continue;
        ++analyzed_module_count;
        const double angle = std::atan2(
            std::sqrt(module->x * module->x + module->y * module->y),
            gRunConfig.hycal_z);
        const double expected_energy = analysis::PhysicsTools::ExpectedEnergy(
            angle, gRunConfig.Ebeam, "ep");
        for (int i = 0; i < grids; ++i) {
            for (int j = 0; j < grids; ++j) {
                if (merged->h1_energy_grid_allModule[m][i][j]->GetEntries() < 200) continue;
                auto fit = analysis::PhysicsTools::fitPeak(
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
    const std::string output_json_name = output_name + ".json";
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

    const std::string output_file_name = output_name + ".root";
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
                energy_grids[m][i][j]->Write();
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
        if (event.cl_nblocks[0] < 3) continue;
        if (fabs(event.cl_energy[0] - run_config.Ebeam) > 3.0 * 0.03 * std::sqrt(run_config.Ebeam * 1000.)) continue;

        HCHit hc_hit;
        GEMHit gem_hit;
        hc_hit.x = event.cl_x[0];
        hc_hit.y = event.cl_y[0];
        hc_hit.z = event.cl_z[0];
        hc_hit.energy = event.cl_energy[0];
        if (event.matchNum == 1){
            // [0][0]: the matched hit of the downstream GEM pair (GEM1/GEM2)
            gem_hit.x = event.mHit_gx[0][0];
            gem_hit.y = event.mHit_gy[0][0];
            gem_hit.z = event.mHit_gz[0][0];
        }

        if (gem_hit.z != 0.f) GetProjection(gem_hit, hc_hit.z);
        ApplyToHyCal(gem_hit, run_config);
        ApplyToHyCal(hc_hit, run_config);
        if (!InHyCalRing(hc_hit.x, hc_hit.y, 2.0, 16.)) continue;

        const auto *mod = hycal.module_by_id(event.cl_center[0]);
        if (!mod) continue;
        const auto [xd_hycal, yd_hycal] = mod->cell_offset(hc_hit.x, hc_hit.y);
        const auto [xd_gem, yd_gem] = mod->cell_offset(gem_hit.x, gem_hit.y);

        // Fill the selected module's energy-grid histogram.
        int col, row;
        if (!fdec::HyCalEnergyBias::cell(*mod, hc_hit.x, hc_hit.y, col, row)) continue;
        const int all_module_index = mod->id - 1001;
        {
            std::lock_guard<std::mutex> lock(fill_locks->module[all_module_index]);
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