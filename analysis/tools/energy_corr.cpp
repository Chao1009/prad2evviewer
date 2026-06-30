
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
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <thread>
#include <mutex>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

using EventVars_Recon = prad2::ReconEventData;

long long process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal,
    std::map<int, TH1F*> &energy_hists, PhysicsTools &physics, float Ebeam, int max_events = -1,
    const std::string &label = "", std::mutex *io_mtx = nullptr);

float resolution = 0.035; // pre-defined energy resolution

float E3p5 = 3485.41f; // Energy for 3.5 GeV beam
float E2p2 = 2239.51f; // Energy for 2.2 GeV beam
float E0p7 = 728.9f;  // Energy for 0.7 GeV beam

static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find("_recon.root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

bool Vetoed(float cl_time, float sci_time, float sci_int){
    // Simple veto logic: if the cluster time is within a certain window of the scintillator time, and the scintillator signal is above a threshold, we consider it a vetoed event.
    const float time_shift = 35.f; // ns
    const float time_window = 7.f; // ns
    const float int_threshold = 2000.f; // arbitrary units
    return (fabs(cl_time - sci_time - time_shift) < time_window) && (sci_int > int_threshold);
}

int main(int argc, char *argv[]){

    // --- parse command line arguments ---


    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);
}

long long process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal,
    std::map<int, TH1F*> &energy_hists, PhysicsTools &physics, float Ebeam, int max_events,
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

        for( int j = 0; j < ev.n_clusters; j++) {
            int mod_id = ev.cl_center[j];
            if (ev.cl_nblocks[j] <= 3) continue;
            auto mod = hycal.module_by_id(mod_id);
            if ( !mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals

            float mod_x = (float)mod->x;
            float mod_y = (float)mod->y;
            float mod_size_x = (float)mod->size_x;
            float mod_size_y = (float)mod->size_y;

            float c_x = ev.cl_x[j], c_y = ev.cl_y[j], c_z = ev.cl_z[j];

            // hit in a single module grid
            float xd = (c_x - mod_x) / mod_size_x;
            float yd = (c_y - mod_y) / mod_size_y;

            float theta = std::atan2(std::sqrt(c_x*c_x + c_y*c_y), c_z) * 180.f / M_PI;
            float energy = ev.cl_energy[j];

            bool veto = false;
                float sci_time, sci_int;
                for(int k = 0; k < ev.veto_nch; k++){
                    for(int p = 0; p < ev.veto_npeaks[k]; p++){
                        sci_time = ev.veto_peak_time[k][p];
                        sci_int = ev.veto_peak_integral[k][p];
                        veto = Vetoed(ev.cl_time[j], sci_time, sci_int);
                        if(veto) break;
                    }
                    if(veto) break;
                }
                if(theta > 1.3) veto = false;

            if(veto && energy > 600. && Ebeam < 1000.f) continue;

            energy_hists[mod_id]->Fill(energy);
        }
    }
    return n_accepted;
}
