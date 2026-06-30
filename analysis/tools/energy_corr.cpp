// energy_corr.cpp : to correct the reconstructed energy non-uniformity depending on
// the position of the cluster within the calorimeter modules.
//
// Measure the ratio of the reconstructed energy to the expected energy for different
// positions within each module. For example, each module can be divided into a grid,
// First try to a grid of 5 by 5, the map is below, could be saved into a 2D array of 1D hist
// for each module, each cell in the 5x5 grid, you can fill a 1D histogram with the energy ratio for that cell.
// column   0   1   2   3   4
// row     +---+---+---+---+---+     beam top  ^
//  0      |   |   |   |   |   |               |
//  1      |   |   |   |   |   |
//  2      |   |   |   |   |   |     beam right ->
//  3      |   |   |   |   |   |
//  4      |   |   |   |   |   |
//         +---+---+---+---+---+
//
// Plan:
// 1. Loop over events in the TTree.
// 2. For each event, determine the module and the position within the module.
// 3. Need a histogram map
//    This map will hold the 1D histograms for each module and each cell in the 5x5 grid.
//    The key can be the module number, and the value can be 2D array of 25 histograms (for the 5x5 grid).
// 3. Fill the corresponding 1D histogram for the energy ratio in the 5x5 grid cell.
// 4. After processing all events, fit the histograms to get the correction factors to correct the energy non-uniformity.

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

long long process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal,
    EnergyHistMap &energy_hists, PhysicsTools &physics, float Ebeam, int max_events = -1,
    const std::string &label = "", std::mutex *io_mtx = nullptr);

float resolution = 0.035; // pre-defined energy resolution

float E3p5 = 3485.41f; // Energy for 3.5 GeV beam
float E2p2 = 2239.51f; // Energy for 2.2 GeV beam
float E0p7 = 728.9f;  // Energy for 0.7 GeV beam

bool Vetoed(float cl_time, float sci_time, float sci_int){
    // Simple veto logic: if the cluster time is within a certain window of the scintillator time, and the scintillator signal is above a threshold, we consider it a vetoed event.
    const float time_shift = 35.f; // ns
    const float time_window = 7.f; // ns
    const float int_threshold = 2000.f; // arbitrary units
    return (fabs(cl_time - sci_time - time_shift) < time_window) && (sci_int > int_threshold);
}

int main(int argc, char *argv[]){

    // --- parse command line arguments ---
    std::vector<std::string> input_paths;
    std::string output_path;
    int max_events = -1;

    for (int i = 1; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "-i") {
            if (++i >= argc || argv[i][0] == '-') {
                std::cerr << "-i requires at least one input ROOT file.\n";
                return 1;
            }
            do {
                input_paths.emplace_back(argv[i]);
                ++i;
            } while (i < argc && argv[i][0] != '-');
            --i;
        } else if (arg == "-o") {
            if (++i >= argc) {
                std::cerr << "-o requires an output ROOT file.\n";
                return 1;
            }
            output_path = argv[i];
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

    if (input_paths.empty() || output_path.empty()) {
        std::cerr << "Usage: " << argv[0]
                  << " -i input1.root [input2.root ...] -o output.root [-n max_events]\n";
        return 1;
    }

    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);

    // One 5x5 histogram grid for each PbWO4 module.
    EnergyHistMap energy_hists;
    for (int i = 0; i < hycal.module_count(); ++i) {
        const auto &module = hycal.module(i);
        if (!module.is_pwo4()) continue;
        // exclude outmost layer modules and under absorber modules
        if (fabs(module.x) > 16.f * 20.77f || fabs(module.y) > 16.f * 20.75f) continue;
        if (fabs(module.x) <  2.f * 20.77f && fabs(module.y) <  2.f * 20.75f) continue;

        auto &grid = energy_hists[module.id];
        for (int row = 0; row < kGridSize; ++row) {
            for (int col = 0; col < kGridSize; ++col) {
                const std::string name = "h_energy_ratio_" + module.name + "_r"
                    + std::to_string(row) + "_c" + std::to_string(col);
                const std::string title = module.name + " energy ratio, cell ("
                    + std::to_string(row) + ", " + std::to_string(col)
                    + ");E_{reconstructed}/E_{expected};Counts";
                grid[row][col] = new TH1F(name.c_str(), title.c_str(), 80, 0.8, 1.2);
                grid[row][col]->SetDirectory(nullptr);
            }
        }
    }

    TChain tree("recon");
    for (const auto &file : input_paths) {
        if (tree.Add(file.c_str()) == 0) {
            std::cerr << "Failed to add input ROOT file: " << file << '\n';
            return 1;
        }
    }
    if (tree.GetEntries() <= 0) {
        std::cerr << "No entries found in input ROOT files.\n";
        return 1;
    }

    EventVars_Recon ev;
    prad2::SetReconReadBranches(&tree, ev);

    process_event(&tree, ev, hycal, energy_hists, physics, E3p5, max_events, "3.5GeV");

    TFile output_file(output_path.c_str(), "RECREATE");
    if (output_file.IsZombie()) {
        std::cerr << "Failed to create output ROOT file: " << output_path << '\n';
        return 1;
    }
    for (auto &[module_id, grid] : energy_hists) {
        (void)module_id;
        for (auto &row : grid)
            for (TH1F *hist : row)
                hist->Write();
    }
    output_file.Close();
}

long long process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal,
    EnergyHistMap &energy_hists, PhysicsTools &physics, float Ebeam, int max_events,
    const std::string &label, std::mutex *io_mtx)
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
    for (int i = 0; i < tree->GetEntries(); i++) {
        // entry-count cap; checked before the trigger cut so -n N stops at entry N
        if (max_events > 0 && i >= max_events) {
            log_msg("[" + label + "] Reached max events limit: "
                    + std::to_string(max_events) + "\n");
            break;
        }
        tree->GetEntry(i);
        if( i % 10000 == 0) {
            log_msg("[" + label + "] Processing event " + std::to_string(i)
                    + "/" + std::to_string(tree->GetEntries()) + "\r", true);
        }
        if ((ev.trigger_bits & prad2::TBIT_sum) == 0) continue;
        n_accepted++;

        if (ev.n_clusters != 1) continue;

        int mod_id = ev.cl_center[0];
        if (ev.cl_nblocks[0] <= 3) continue;
        auto mod = hycal.module_by_id(mod_id);
        if ( !mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals
        auto hist_it = energy_hists.find(mod_id);
        if (hist_it == energy_hists.end()) continue; // excluded during histogram setup

        float mod_x = (float)mod->x;
        float mod_y = (float)mod->y;
        float mod_size_x = (float)mod->size_x;
        float mod_size_y = (float)mod->size_y;

        float c_x = ev.cl_x[0], c_y = ev.cl_y[0], c_z = ev.cl_z[0];

        // hit in a single module grid
        float xd = (c_x - mod_x) / mod_size_x;
        float yd = (c_y - mod_y) / mod_size_y;

        float theta = std::atan2(std::sqrt(c_x*c_x + c_y*c_y), c_z) * 180.f / M_PI;
        float energy = ev.cl_energy[0];

        bool veto = false;
            float sci_time, sci_int;
            for(int k = 0; k < ev.veto_nch; k++){
                for(int p = 0; p < ev.veto_npeaks[k]; p++){
                    sci_time = ev.veto_peak_time[k][p];
                    sci_int = ev.veto_peak_integral[k][p];
                    veto = Vetoed(ev.cl_time[0], sci_time, sci_int);
                    if(veto) break;
                }
                if(veto) break;
            }
            if(theta > 1.3) veto = false;

        if(veto && energy > 600. && Ebeam < 1000.f) continue;

        const int col = static_cast<int>(std::floor((xd + 0.5f) * kGridSize));
        const int row = static_cast<int>(std::floor((0.5f - yd) * kGridSize));
        if (row < 0 || row >= kGridSize || col < 0 || col >= kGridSize) continue;

        const float expected_energy = PhysicsTools::ExpectedEnergy(theta, Ebeam, "ep");
        if (expected_energy <= 0.f) continue;
        hist_it->second[row][col]->Fill(energy / expected_energy);
    }
    return n_accepted;
}
