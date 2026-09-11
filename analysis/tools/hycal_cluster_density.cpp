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
#include <TKey.h>
#include <TLatex.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TROOT.h>
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
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;

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

static std::string outputFileName(const std::string &output_prefix)
{
    return output_prefix + "hycal_cluster_density.root";
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
    return (fabs(xmm) > module * 2.2 || fabs(ymm) > module * 2.2)
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
    std::string output_path_name, daq_config_file, recon_config_file, gem_ped_file;
    int  max_events  = -1;
    int  num_threads = 4;
    int  num_files   = -1;
    bool worker_mode = false;

    int opt;
    while ((opt = getopt(argc, argv, "o:n:f:j:w")) != -1) {
        switch (opt) {
            case 'o': output_path_name = optarg; break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'f': num_files        = std::atoi(optarg); break;
            case 'j': num_threads     = std::atoi(optarg); break;
            case 'w': worker_mode      = true; break;
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
        std::cerr << "Usage: hycal_shower_profile <input_recon.root|dir> [more...] "
                     "[-o ./output_name(no extension)] [-n max_events] [-f nfiles] [-j threads]\n";
        return 1;
    }

    if (output_path_name.empty()) {
        std::cerr << "No output prefix provided. Please pass -o <output_prefix>.\n";
        return 1;
    }

    // ROOT TTree branch addresses are not safe to share across worker threads.
    // Match hycal_shower_profile: run one isolated worker process per input
    // file, then merge the worker histograms in this parent process.
    if (!worker_mode && root_files.size() > 1) {
        num_threads = std::max(1, std::min(num_threads,
                                           static_cast<int>(root_files.size())));
        const std::string executable = shell_quote(argv[0]);
        std::vector<std::string> worker_outputs(root_files.size());
        std::atomic<size_t> next_file{0};
        std::mutex worker_mutex;
        std::vector<std::future<void>> workers;

        auto run_worker = [&]() {
            while (true) {
                const size_t file_index = next_file.fetch_add(1);
                if (file_index >= root_files.size()) return;

                const std::string worker_prefix =
                    output_path_name + ".worker_" + std::to_string(file_index);
                worker_outputs[file_index] = outputFileName(worker_prefix);
                std::string command = executable + " -w -o "
                    + shell_quote(worker_prefix);
                if (max_events > 0)
                    command += " -n " + std::to_string(max_events);
                command += " -j 1 " + shell_quote(root_files[file_index]);

                const int status = std::system(command.c_str());
                std::lock_guard<std::mutex> lock(worker_mutex);
                if (status != 0)
                    std::cerr << "Worker failed for " << root_files[file_index]
                              << " (status " << status << ")\n";
            }
        };
        for (int i = 0; i < num_threads; ++i)
            workers.push_back(std::async(std::launch::async, run_worker));
        for (auto &worker : workers) worker.get();

        const std::string merged_output = outputFileName(output_path_name);
        TFile *merged = TFile::Open(merged_output.c_str(), "RECREATE");
        if (!merged || merged->IsZombie()) {
            std::cerr << "Cannot create merged output " << merged_output << "\n";
            return 1;
        }

        std::map<std::string, std::unique_ptr<TH1>> merged_histograms;
        for (const auto &worker_output : worker_outputs) {
            TFile *input = TFile::Open(worker_output.c_str(), "READ");
            if (!input || input->IsZombie()) {
                if (input) delete input;
                continue;
            }
            TIter keys(input->GetListOfKeys());
            while (auto *key = dynamic_cast<TKey *>(keys())) {
                TObject *object = key->ReadObj();
                auto *hist = dynamic_cast<TH1 *>(object);
                if (!hist) {
                    delete object;
                    continue;
                }

                const std::string name = hist->GetName();
                auto it = merged_histograms.find(name);
                if (it == merged_histograms.end()) {
                    std::unique_ptr<TH1> copy(
                        dynamic_cast<TH1 *>(hist->Clone(name.c_str())));
                    if (copy) {
                        copy->SetDirectory(nullptr);
                        merged_histograms.emplace(name, std::move(copy));
                    }
                } else {
                    it->second->Add(hist);
                }
                delete object;
            }
            input->Close();
            delete input;
        }

        merged->cd();
        for (auto &[name, hist] : merged_histograms) {
            hist->SetDirectory(merged);
            hist->Write(name.c_str(), TObject::kOverwrite);
            hist->SetDirectory(nullptr);
        }
        merged->Close();
        delete merged;
        for (const auto &worker_output : worker_outputs)
            std::remove(worker_output.c_str());
        return 0;
    }

    TChain tree("recon");
    for (const auto &file : root_files) {
        tree.Add(file.c_str());
    }
    EventVars_Recon ev;
    prad2::SetReconReadBranches(&tree, ev);

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(db_dir + "/hycal_map.json");

    // Histograms
    TH2F *h2_hit_hycal = new TH2F("h2_hit_hycal", "HyCal Hit Distribution;(X_{hycal}-X_{cell center})/d_{cell size};(Y_{hycal}-Y_{cell center})/d_{cell size}", 200, -1.0, 1.0, 200, -1.0, 1.0);
    TH2F *h2_hit_gem = new TH2F("h2_hit_gem", "GEM Hit Distribution;(X_{gem}-X_{cell center})/d_{cell size};(Y_{gem}-Y_{cell center})/d_{cell size}", 200, -1.0, 1.0, 200, -1.0, 1.0);
    TH2F *h2_dist_dx_xd = new TH2F("h2_dist_dx_xd", "dx vs xd_hycal;relative x to cell center;x_hycal-x_gem [mm]", 200, -1.0, 1.0, 400, -5.0, 5.0);
    TH2F *h2_dist_dy_yd = new TH2F("h2_dist_dy_yd", "dy vs yd_hycal;relative y to cell center;y_hycal-y_gem [mm]", 200, -1.0, 1.0, 400, -5.0, 5.0);
    TH1F *h1_density_hycal_xd = new TH1F("h1_density_hycal_xd", "Density of HyCal xd", 200, -1.0, 1.0);
    TH1F *h1_density_hycal_yd = new TH1F("h1_density_hycal_yd", "Density of HyCal yd", 200, -1.0, 1.0);
    TH1F *h1_density_gem_xd = new TH1F("h1_density_gem_xd", "Density of GEM xd", 200, -1.0, 1.0);
    TH1F *h1_density_gem_yd = new TH1F("h1_density_gem_yd", "Density of GEM yd", 200, -1.0, 1.0);

    Long64_t n = tree.GetEntries();
    if (max_events >= 0 && max_events < n) n = max_events;

    for (Long64_t i = 0; i < n; i++) {
        tree.GetEntry(i);
        if (i % 10000 == 0) std::cout << "Processed " << i << " / " << n << " entries.\r" << std::flush;

        // trigger selection
        bool is_sum      = (ev.trigger_bits & prad2::TBIT_sum) != 0;
        if (!is_sum) continue;

        //Event selection, single cluster e-p events, no "kSplit" flag
        if (ev.n_clusters != 1 || ev.matchNum != 1) continue;
        if (ev.cl_nblocks[0] < 3) continue;
        if (fdec::test_bit(ev.cl_flag[0], fdec::kSplit)) continue;
        if (std::fabs(ev.cl_energy[0] - gRunConfig.Ebeam) > 3.0f * 0.03f * std::sqrt(gRunConfig.Ebeam * 1000.f)) continue;

        HCHit hc_hit;
        GEMHit g_hit;

        hc_hit.x = ev.cl_x[0];
        hc_hit.y = ev.cl_y[0];
        hc_hit.z = ev.cl_z[0];
        hc_hit.energy = ev.cl_energy[0];

        g_hit.x = ev.mHit_gx[0][0];
        g_hit.y = ev.mHit_gy[0][0];
        g_hit.z = ev.mHit_gz[0][0];

        float scale = hc_hit.z / g_hit.z;
        g_hit.x *= scale;
        g_hit.y *= scale;
        g_hit.z *= scale;

        ApplyToHyCal(g_hit, gRunConfig);
        ApplyToHyCal(hc_hit, gRunConfig);

        // only look at one module first, W566
        const auto &mod = hycal.module_by_id(1567+34);
        if ( !( ev.cl_x[0] < mod->x + mod->size_x / 2. && ev.cl_x[0] > mod->x - mod->size_x / 2. &&
                ev.cl_y[0] < mod->y + mod->size_y / 2. && ev.cl_y[0] > mod->y - mod->size_y / 2. ) ) continue;

        float dx = hc_hit.x - g_hit.x;
        float dy = hc_hit.y - g_hit.y;
        
        float xd_hycal = (hc_hit.x - mod->x) / mod->size_x;
        float yd_hycal = (hc_hit.y - mod->y) / mod->size_y;

        float xd_gem = (g_hit.x - mod->x) / mod->size_x;
        float yd_gem = (g_hit.y - mod->y) / mod->size_y;

        h2_hit_hycal->Fill(xd_hycal, yd_hycal);
        h2_hit_gem->Fill(xd_gem, yd_gem);
        h2_dist_dx_xd->Fill(xd_hycal, dx);
        h2_dist_dy_yd->Fill(yd_hycal, dy);
        h1_density_hycal_xd->Fill(xd_hycal);
        h1_density_hycal_yd->Fill(yd_hycal);
        h1_density_gem_xd->Fill(xd_gem);
        h1_density_gem_yd->Fill(yd_gem);

    }

    TFile output_file(outputFileName(output_path_name).c_str(), "RECREATE");
    h2_hit_hycal->Write();
    h2_hit_gem->Write();
    h2_dist_dx_xd->Write();
    h2_dist_dy_yd->Write();
    h1_density_hycal_xd->Write();
    h1_density_hycal_yd->Write();
    h1_density_gem_xd->Write();
    h1_density_gem_yd->Write();
    output_file.Close();

}