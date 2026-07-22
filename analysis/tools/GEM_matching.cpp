// A quick check tool to test the matching result between HyCal clusters and GEM hits in the replay output
// Usage:
//   matching <input_recon.root|dir> [more files...] [-o out.root] [-n max_events]
//   -o  output ROOT file (default: input filename with _matching.root suffix)
//   -n  max events to process (default: all)
// Example:
//   matching recon.root -o recon_matching.root -n 10000
//   matching recon_dir/ recon.root...  -n 100000

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "MatchingTools.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>

#include <iostream>
#include <string>
#include <vector>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <algorithm>
#include <unistd.h>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output = "matching_result.root";
    
    int max_events = -1, nfiles = -1;
    int opt;
    while ((opt = getopt(argc, argv, "o:n:f:")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'n':
                if (!optarg) {
                    std::cerr << "Option -n requires an argument.\n";
                    return 1;
                }
                max_events = std::atoi(optarg);
                break;
            case 'f':
                if (!optarg) {
                    std::cerr << "Option -f requires an argument.\n";
                    return 1;
                }
                nfiles = std::atoi(optarg);
                break;
            default:
                std::cerr << "Usage: GEM_matching <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-f nfiles]\n";
                return 1;
        }
    }
    // collect input files (can be files, directories, or mixed)
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; i++) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-f nfiles]\n";
        return 1;
    }

    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- load run configuration from JSON ---
    int run_num = get_run_int(root_files[0]);
    gRunConfig = LoadRunConfig(dbDir + "/runinfo/general.json", run_num);
    
    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);
    MatchingTools matching;

    // --- setup TChain and branches ---
    TChain *chain = new TChain("recon");
    int n_added_files = 0;
    for (const auto &f : root_files) {
        chain->Add(f.c_str());
        std::cerr << "Added file: " << f << "\n";
        n_added_files++;
        if (nfiles > 0 && n_added_files >= nfiles) break;
    }
    TTree *tree = chain;
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in input files\n";
        return 1;
    }

    EventVars_Recon ev;
    prad2::SetReconReadBranches(tree, ev);
    prad2::ReconMatchVectorBindings match_bindings;
    prad2::BindReconMatchVectorBranches(tree, ev, match_bindings);

    //setup histograms here
    int bin_num = 400, bin_lo = -20, bin_hi = 20;
    // select Mott events from Raw Sum trigger bits
    TH1F *h1_deltaX_gem = new TH1F("h1_deltaX_gem", "Delta X Between GEMs;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_gem = new TH1F("h1_deltaY_gem", "Delta Y Between GEMs;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_gem = new TH2F("h2_deltaXY_gem", "Delta X vs Delta Y Between GEMs;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaX_hycal = new TH1F("h1_deltaX_hycal", "Delta X Between GEMs and HyCal;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_hycal = new TH1F("h1_deltaY_hycal", "Delta Y Between GEMs and HyCal;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_hycal = new TH2F("h2_deltaXY_hycal", "Delta X vs Delta Y Between GEMs and HyCal;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_Nhits_matched[4];
    for (int i = 0; i < 4; i++) {
        h1_Nhits_matched[i] = new TH1F(Form("h1_Nhits_matched_%d", i), Form("Number of Hits in Matching Radius for GEM %d;N_{hits};Entries", i), 30, 0, 30);
    }
    bin_num = 600; bin_lo = -30; bin_hi = 30;
    // select 3cl events from 3-cluster trigger bits
    TH1F *h1_deltaX_gem_100MeV = new TH1F("h1_deltaX_gem_100MeV", "Delta X Between GEMs for 100 MeV;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_gem_100MeV = new TH1F("h1_deltaY_gem_100MeV", "Delta Y Between GEMs for 100 MeV;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_gem_100MeV = new TH2F("h2_deltaXY_gem_100MeV", "Delta X vs Delta Y Between GEMs for 100 MeV;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaX_hycal_100MeV = new TH1F("h1_deltaX_hycal_100MeV", "Delta X Between GEMs and HyCal for 100 MeV;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_hycal_100MeV = new TH1F("h1_deltaY_hycal_100MeV", "Delta Y Between GEMs and HyCal for 100 MeV;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_hycal_100MeV = new TH2F("h2_deltaXY_hycal_100MeV", "Delta X vs Delta Y Between GEMs and HyCal for 100 MeV;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaX_gem_300MeV = new TH1F("h1_deltaX_gem_300MeV", "Delta X Between GEMs for 300 MeV;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_gem_300MeV = new TH1F("h1_deltaY_gem_300MeV", "Delta Y Between GEMs for 300 MeV;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_gem_300MeV = new TH2F("h2_deltaXY_gem_300MeV", "Delta X vs Delta Y Between GEMs for 300 MeV;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaX_hycal_300MeV = new TH1F("h1_deltaX_hycal_300MeV", "Delta X Between GEMs and HyCal for 300 MeV;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_hycal_300MeV = new TH1F("h1_deltaY_hycal_300MeV", "Delta Y Between GEMs and HyCal for 300 MeV;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_hycal_300MeV = new TH2F("h2_deltaXY_hycal_300MeV", "Delta X vs Delta Y Between GEMs and HyCal for 300 MeV;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaX_gem_500MeV = new TH1F("h1_deltaX_gem_500MeV", "Delta X Between GEMs for 500 MeV;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_gem_500MeV = new TH1F("h1_deltaY_gem_500MeV", "Delta Y Between GEMs for 500 MeV;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_gem_500MeV = new TH2F("h2_deltaXY_gem_500MeV", "Delta X vs Delta Y Between GEMs for 500 MeV;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaX_hycal_500MeV = new TH1F("h1_deltaX_hycal_500MeV", "Delta X Between GEMs and HyCal for 500 MeV;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_hycal_500MeV = new TH1F("h1_deltaY_hycal_500MeV", "Delta Y Between GEMs and HyCal for 500 MeV;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_hycal_500MeV = new TH2F("h2_deltaXY_hycal_500MeV", "Delta X vs Delta Y Between GEMs and HyCal for 500 MeV;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaX_gem_1000MeV = new TH1F("h1_deltaX_gem_1000MeV", "Delta X Between GEMs for 1000 MeV;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_gem_1000MeV = new TH1F("h1_deltaY_gem_1000MeV", "Delta Y Between GEMs for 1000 MeV;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_gem_1000MeV = new TH2F("h2_deltaXY_gem_1000MeV", "Delta X vs Delta Y Between GEMs for 1000 MeV;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaX_hycal_1000MeV = new TH1F("h1_deltaX_hycal_1000MeV", "Delta X Between GEMs and HyCal for 1000 MeV;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_hycal_1000MeV = new TH1F("h1_deltaY_hycal_1000MeV", "Delta Y Between GEMs and HyCal for 1000 MeV;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_hycal_1000MeV = new TH2F("h2_deltaXY_hycal_1000MeV", "Delta X vs Delta Y Between GEMs and HyCal for 1000 MeV;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaX_gem_2000MeV = new TH1F("h1_deltaX_gem_2000MeV", "Delta X Between GEMs for 2000 MeV;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_gem_2000MeV = new TH1F("h1_deltaY_gem_2000MeV", "Delta Y Between GEMs for 2000 MeV;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_gem_2000MeV = new TH2F("h2_deltaXY_gem_2000MeV", "Delta X vs Delta Y Between GEMs for 2000 MeV;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaX_hycal_2000MeV = new TH1F("h1_deltaX_hycal_2000MeV", "Delta X Between GEMs and HyCal for 2000 MeV;#DeltaX [mm];Entries", bin_num, bin_lo, bin_hi);
    TH1F *h1_deltaY_hycal_2000MeV = new TH1F("h1_deltaY_hycal_2000MeV", "Delta Y Between GEMs and HyCal for 2000 MeV;#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi);
    TH2F *h2_deltaXY_hycal_2000MeV = new TH2F("h2_deltaXY_hycal_2000MeV", "Delta X vs Delta Y Between GEMs and HyCal for 2000 MeV;#DeltaX [mm];#DeltaY [mm];Entries", bin_num, bin_lo, bin_hi, bin_num, bin_lo, bin_hi);
    TH1F *h1_Nhits_matched_100MeV[4], *h1_Nhits_matched_300MeV[4], *h1_Nhits_matched_500MeV[4], *h1_Nhits_matched_1000MeV[4], *h1_Nhits_matched_2000MeV[4];
    for (int i = 0; i < 4; i++) {
        h1_Nhits_matched_100MeV[i] = new TH1F(Form("h1_Nhits_matched_100MeV_%d", i), Form("Number of Hits in Matching Radius for GEM %d at 100 MeV;N_{hits};Entries", i), 30, 0, 30);
        h1_Nhits_matched_300MeV[i] = new TH1F(Form("h1_Nhits_matched_300MeV_%d", i), Form("Number of Hits in Matching Radius for GEM %d at 300 MeV;N_{hits};Entries", i), 30, 0, 30);
        h1_Nhits_matched_500MeV[i] = new TH1F(Form("h1_Nhits_matched_500MeV_%d", i), Form("Number of Hits in Matching Radius for GEM %d at 500 MeV;N_{hits};Entries", i), 30, 0, 30);
        h1_Nhits_matched_1000MeV[i] = new TH1F(Form("h1_Nhits_matched_1000MeV_%d", i), Form("Number of Hits in Matching Radius for GEM %d at 1000 MeV;N_{hits};Entries", i), 30, 0, 30);
        h1_Nhits_matched_2000MeV[i] = new TH1F(Form("h1_Nhits_matched_2000MeV_%d", i), Form("Number of Hits in Matching Radius for GEM %d at 2000 MeV;N_{hits};Entries", i), 30, 0, 30);
    }

    // event loop
    Long64_t nentries = tree->GetEntries();
    if (max_events > 0 && max_events < nentries) nentries = max_events;
    for (Long64_t i = 0; i < nentries; i++) {
        if (i % 10000 == 0) {
            std::cerr << "Processing event " << i << " / " << nentries << "\r" << std::flush;
        }
        tree->GetEntry(i);

        // trigger selection
        bool is_3cluster = (ev.trigger_bits & prad2::TBIT_3cl) != 0;
        bool is_sum      = (ev.trigger_bits & prad2::TBIT_sum) != 0;
        if (!is_3cluster && !is_sum) continue;

        if (is_sum) {
            // Process sum trigger events
            if(ev.n_clusters == 1 && ev.matchNum == 1 && ev.cl_nblocks[0] > 1 &&
                (ev.cl_energy[0] - gRunConfig.Ebeam) < 3. * 0.033 * sqrt(gRunConfig.Ebeam * 1000.)) 
            {
                HCHit hc_hit; // the HyCal cluster 
                GEMHit gem_hit_match[2]; // the best matched GEM hits, 0 for upstream, 1 for downstream
                
                hc_hit = {ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]};
                gem_hit_match[0] = {ev.mHit_gx[0][1], ev.mHit_gy[0][1], ev.mHit_gz[0][1], ev.mHit_gid[0][1]};
                gem_hit_match[1] = {ev.mHit_gx[0][0], ev.mHit_gy[0][0], ev.mHit_gz[0][0], ev.mHit_gid[0][0]};

                GetProjection(gem_hit_match[0], gem_hit_match[1].z);
                float dx_gem = gem_hit_match[1].x - gem_hit_match[0].x;
                float dy_gem = gem_hit_match[1].y - gem_hit_match[0].y;
                h1_deltaX_gem->Fill(dx_gem);
                h1_deltaY_gem->Fill(dy_gem);
                h2_deltaXY_gem->Fill(dx_gem, dy_gem);

                GetProjection(gem_hit_match[0], hc_hit.z);
                float dx_hc = hc_hit.x - gem_hit_match[0].x;
                float dy_hc = hc_hit.y - gem_hit_match[0].y;
                h1_deltaX_hycal->Fill(dx_hc);
                h1_deltaY_hycal->Fill(dy_hc);
                h2_deltaXY_hycal->Fill(dx_hc, dy_hc);

                int N_matched[4] = {0, 0, 0, 0};
                int cl_idx = ev.mHit_cl_index[0];
                const size_t n_match = ev.match_cl_idx.size();
                for (size_t i = 0; i < n_match; ++i) {
                    const int det = static_cast<int>(ev.match_det_id[i]);
                    if (det >= 0 && det < 4) {
                        N_matched[det]++;
                    }
                }

                for (int j = 0; j < 4; ++j) {
                    if(N_matched[j] > 0) h1_Nhits_matched[j]->Fill(N_matched[j]);
                }
            }
        }
        if(is_3cluster){
            for(int j = 0; j < ev.matchNum; ++j) {
                // Process each of the clusters
                const int cl_idx = static_cast<int>(ev.mHit_cl_index[j]);
                if (cl_idx < 0 || cl_idx >= ev.n_clusters || cl_idx >= prad2::kMaxClusters) continue;
                if(ev.cl_nblocks[cl_idx] <= 1) continue; // Skip clusters with 1 or fewer blocks
                HCHit hc_hit; // the HyCal cluster 
                GEMHit gem_hit_match[2]; // the best matched GEM hits, 0 for upstream, 1 for downstream

                hc_hit = {ev.cl_x[cl_idx], ev.cl_y[cl_idx], ev.cl_z[cl_idx], ev.cl_energy[cl_idx]}; // Initialize the HyCal cluster hit
                gem_hit_match[0] = {ev.mHit_gx[j][1], ev.mHit_gy[j][1], ev.mHit_gz[j][1], ev.mHit_gid[j][1]};
                gem_hit_match[1] = {ev.mHit_gx[j][0], ev.mHit_gy[j][0], ev.mHit_gz[j][0], ev.mHit_gid[j][0]};

                GetProjection(gem_hit_match[0], gem_hit_match[1].z);
                float dx_gem = gem_hit_match[1].x - gem_hit_match[0].x;
                float dy_gem = gem_hit_match[1].y - gem_hit_match[0].y;
                GetProjection(gem_hit_match[0], hc_hit.z);
                float dx_hc = hc_hit.x - gem_hit_match[0].x;
                float dy_hc = hc_hit.y - gem_hit_match[0].y;
                int N_matched[4] = {0, 0, 0, 0};
                const size_t n_match = ev.match_cl_idx.size();
                for (size_t i = 0; i < n_match; ++i) {
                    const int det = static_cast<int>(ev.match_det_id[i]);
                    const int cl = static_cast<int>(ev.match_cl_idx[i]);
                    if (det >= 0 && det < 4 && cl == cl_idx) {
                        N_matched[det]++;
                    }
                }
                
                if(fabs(ev.cl_energy[cl_idx] - 2000.) < 3. * 0.033 * sqrt(2000. * 1000.)) {
                    h1_deltaX_gem_2000MeV->Fill(dx_gem);
                    h1_deltaY_gem_2000MeV->Fill(dy_gem);
                    h2_deltaXY_gem_2000MeV->Fill(dx_gem, dy_gem);
                    h2_deltaXY_hycal_2000MeV->Fill(dx_hc, dy_hc);
                    h1_deltaX_hycal_2000MeV->Fill(dx_hc);
                    h1_deltaY_hycal_2000MeV->Fill(dy_hc);
                    for (int k = 0; k < 4; ++k) if(N_matched[k] > 0) h1_Nhits_matched_2000MeV[k]->Fill(N_matched[k]);
                }
                if(fabs(ev.cl_energy[cl_idx] - 1000.) < 3. * 0.033 * sqrt(1000. * 1000.)) {
                    h1_deltaX_gem_1000MeV->Fill(dx_gem);
                    h1_deltaY_gem_1000MeV->Fill(dy_gem);
                    h2_deltaXY_gem_1000MeV->Fill(dx_gem, dy_gem);
                    h1_deltaX_hycal_1000MeV->Fill(dx_hc);
                    h1_deltaY_hycal_1000MeV->Fill(dy_hc);
                    h2_deltaXY_hycal_1000MeV->Fill(dx_hc, dy_hc);
                    for (int k = 0; k < 4; ++k) if(N_matched[k] > 0) h1_Nhits_matched_1000MeV[k]->Fill(N_matched[k]);
                }
                if(fabs(ev.cl_energy[cl_idx] - 500.) < 3. * 0.033 * sqrt(500. * 1000.)) {
                    h1_deltaX_gem_500MeV->Fill(dx_gem);
                    h1_deltaY_gem_500MeV->Fill(dy_gem);
                    h2_deltaXY_gem_500MeV->Fill(dx_gem, dy_gem);
                    h1_deltaX_hycal_500MeV->Fill(dx_hc);
                    h1_deltaY_hycal_500MeV->Fill(dy_hc);
                    h2_deltaXY_hycal_500MeV->Fill(dx_hc, dy_hc);
                    for (int k = 0; k < 4; ++k) if(N_matched[k] > 0) h1_Nhits_matched_500MeV[k]->Fill(N_matched[k]);
                }
                if(fabs(ev.cl_energy[cl_idx] - 300.) < 3. * 0.033 * sqrt(300. * 1000.)) {
                    h1_deltaX_gem_300MeV->Fill(dx_gem);
                    h1_deltaY_gem_300MeV->Fill(dy_gem);
                    h2_deltaXY_gem_300MeV->Fill(dx_gem, dy_gem);
                    h1_deltaX_hycal_300MeV->Fill(dx_hc);
                    h1_deltaY_hycal_300MeV->Fill(dy_hc);
                    h2_deltaXY_hycal_300MeV->Fill(dx_hc, dy_hc);
                    for (int k = 0; k < 4; ++k) if(N_matched[k] > 0) h1_Nhits_matched_300MeV[k]->Fill(N_matched[k]);
                }
                if(fabs(ev.cl_energy[cl_idx] - 100.) < 3. * 0.033 * sqrt(100. * 1000.)) {
                    h1_deltaX_gem_100MeV->Fill(dx_gem);
                    h1_deltaY_gem_100MeV->Fill(dy_gem);
                    h2_deltaXY_gem_100MeV->Fill(dx_gem, dy_gem);
                    h1_deltaX_hycal_100MeV->Fill(dx_hc);
                    h1_deltaY_hycal_100MeV->Fill(dy_hc);
                    h2_deltaXY_hycal_100MeV->Fill(dx_hc, dy_hc);
                    for (int k = 0; k < 4; ++k) if(N_matched[k] > 0) h1_Nhits_matched_100MeV[k]->Fill(N_matched[k]);
                }
            }
        }
    }

    // End of event loop
    // Output results to the file
    TFile *output_file = TFile::Open(output.c_str(), "RECREATE");
    output_file->cd();
    h1_deltaX_gem->Write();
    h1_deltaY_gem->Write();
    h2_deltaXY_gem->Write();
    h1_deltaX_hycal->Write();
    h1_deltaY_hycal->Write();
    h2_deltaXY_hycal->Write();
    for (int j = 0; j < 4; ++j) {
        h1_Nhits_matched[j]->Write();
    }
    output_file->mkdir("energy_bins");
    output_file->cd("energy_bins");
    h1_deltaX_gem_2000MeV->Write();
    h1_deltaY_gem_2000MeV->Write();
    h2_deltaXY_gem_2000MeV->Write();
    h1_deltaX_hycal_2000MeV->Write();
    h1_deltaY_hycal_2000MeV->Write();
    h2_deltaXY_hycal_2000MeV->Write();
    h1_deltaX_gem_1000MeV->Write();
    h1_deltaY_gem_1000MeV->Write();
    h2_deltaXY_gem_1000MeV->Write();
    h1_deltaX_hycal_1000MeV->Write();
    h1_deltaY_hycal_1000MeV->Write();
    h2_deltaXY_hycal_1000MeV->Write();
    h1_deltaX_gem_500MeV->Write();
    h1_deltaY_gem_500MeV->Write();
    h2_deltaXY_gem_500MeV->Write();
    h1_deltaX_hycal_500MeV->Write();
    h1_deltaY_hycal_500MeV->Write();
    h2_deltaXY_hycal_500MeV->Write();
    h1_deltaX_gem_300MeV->Write();
    h1_deltaY_gem_300MeV->Write();
    h2_deltaXY_gem_300MeV->Write();
    h1_deltaX_hycal_300MeV->Write();
    h1_deltaY_hycal_300MeV->Write();
    h2_deltaXY_hycal_300MeV->Write();
    h1_deltaX_gem_100MeV->Write();
    h1_deltaY_gem_100MeV->Write();
    h2_deltaXY_gem_100MeV->Write();
    h1_deltaX_hycal_100MeV->Write();
    h1_deltaY_hycal_100MeV->Write();
    h2_deltaXY_hycal_100MeV->Write();
    for (int k = 0; k < 4; ++k) {
        h1_Nhits_matched_2000MeV[k]->Write();
        h1_Nhits_matched_1000MeV[k]->Write();
        h1_Nhits_matched_500MeV[k]->Write();
        h1_Nhits_matched_300MeV[k]->Write();
        h1_Nhits_matched_100MeV[k]->Write();
    }

    // Print a concise summary so users can verify content without opening ROOT.
    for (int j = 0; j < 4; ++j) {
        int nonzero_bins = 0;
        for (int b = 1; b <= h1_Nhits_matched[j]->GetNbinsX(); ++b) {
            if (h1_Nhits_matched[j]->GetBinContent(b) > 0) nonzero_bins++;
        }
        std::cerr << "Summary h1_Nhits_matched_" << j
                  << ": entries=" << h1_Nhits_matched[j]->GetEntries()
                  << ", mean=" << h1_Nhits_matched[j]->GetMean()
                  << ", nonzero_bins=" << nonzero_bins
                  << "\n";
    }

    output_file->Close();
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
