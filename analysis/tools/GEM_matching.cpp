// A quick check tool to test the matching result between HyCal clusters and GEM hits in the replay output
// Usage:
//   GEM_matching <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-f nfiles] [-j threads]
//   -o  output ROOT file (default: matching_result.root)
//   -n  max events to process across all files (default: all)
//   -f  max number of input files to add (default: all)
//   -j  number of worker threads (default: 4)

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "MatchingTools.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"

#include <TClass.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TROOT.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <thread>
#include <unistd.h>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

using EventVars_Recon = prad2::ReconEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);

namespace {

struct HistSet {
    static constexpr int kNGem = 4;
    static constexpr int kNEnergy = 5;

    std::unique_ptr<TH1F> h1_deltaX_gem;
    std::unique_ptr<TH1F> h1_deltaY_gem;
    std::unique_ptr<TH2F> h2_deltaXY_gem;
    std::unique_ptr<TH1F> h1_deltaX_hycal;
    std::unique_ptr<TH1F> h1_deltaY_hycal;
    std::unique_ptr<TH2F> h2_deltaXY_hycal;
    std::array<std::unique_ptr<TH1F>, kNGem> h1_Nhits_matched;

    std::array<std::unique_ptr<TH1F>, kNEnergy> h1_deltaX_gem_E;
    std::array<std::unique_ptr<TH1F>, kNEnergy> h1_deltaY_gem_E;
    std::array<std::unique_ptr<TH2F>, kNEnergy> h2_deltaXY_gem_E;
    std::array<std::unique_ptr<TH1F>, kNEnergy> h1_deltaX_hycal_E;
    std::array<std::unique_ptr<TH1F>, kNEnergy> h1_deltaY_hycal_E;
    std::array<std::unique_ptr<TH2F>, kNEnergy> h2_deltaXY_hycal_E;
    std::array<std::array<std::unique_ptr<TH1F>, kNGem>, kNEnergy> h1_Nhits_matched_E;

    explicit HistSet(const std::string &suffix)
    {
        auto make_h1 = [&suffix](const std::string &name,
                                 const std::string &title,
                                 int bins,
                                 double lo,
                                 double hi) {
            auto h = std::make_unique<TH1F>((name + suffix).c_str(), title.c_str(), bins, lo, hi);
            h->SetDirectory(nullptr);
            return h;
        };
        auto make_h2 = [&suffix](const std::string &name,
                                 const std::string &title,
                                 int bins_x,
                                 double lo_x,
                                 double hi_x,
                                 int bins_y,
                                 double lo_y,
                                 double hi_y) {
            auto h = std::make_unique<TH2F>((name + suffix).c_str(), title.c_str(), bins_x, lo_x, hi_x, bins_y, lo_y, hi_y);
            h->SetDirectory(nullptr);
            return h;
        };

        h1_deltaX_gem = make_h1("h1_deltaX_gem", "Delta X Between GEMs;#DeltaX [mm];Entries", 600, -30, 30);
        h1_deltaY_gem = make_h1("h1_deltaY_gem", "Delta Y Between GEMs;#DeltaY [mm];Entries", 600, -30, 30);
        h2_deltaXY_gem = make_h2("h2_deltaXY_gem", "Delta X vs Delta Y Between GEMs;#DeltaX [mm];#DeltaY [mm];Entries", 600, -30, 30, 600, -30, 30);
        h1_deltaX_hycal = make_h1("h1_deltaX_hycal", "Delta X Between GEMs and HyCal;#DeltaX [mm];Entries", 600, -30, 30);
        h1_deltaY_hycal = make_h1("h1_deltaY_hycal", "Delta Y Between GEMs and HyCal;#DeltaY [mm];Entries", 600, -30, 30);
        h2_deltaXY_hycal = make_h2("h2_deltaXY_hycal", "Delta X vs Delta Y Between GEMs and HyCal;#DeltaX [mm];#DeltaY [mm];Entries", 600, -30, 30, 600, -30, 30);
        for (int i = 0; i < kNGem; ++i) {
            h1_Nhits_matched[i] = make_h1(
                Form("h1_Nhits_matched_%d", i),
                Form("Number of Hits in Matching Radius for GEM %d;N_{hits};Entries", i),
                30, 0, 30);
        }

        const std::array<int, kNEnergy> energies = {100, 300, 500, 1000, 1500};
        for (int ie = 0; ie < kNEnergy; ++ie) {
            const int E = energies[ie];
            const std::string tag = Form("_%dMeV", E);
            h1_deltaX_gem_E[ie] = make_h1(
                "h1_deltaX_gem" + tag,
                Form("Delta X Between GEMs for %d MeV;#DeltaX [mm];Entries", E),
                600, -30, 30);
            h1_deltaY_gem_E[ie] = make_h1(
                "h1_deltaY_gem" + tag,
                Form("Delta Y Between GEMs for %d MeV;#DeltaY [mm];Entries", E),
                600, -30, 30);
            h2_deltaXY_gem_E[ie] = make_h2(
                "h2_deltaXY_gem" + tag,
                Form("Delta X vs Delta Y Between GEMs for %d MeV;#DeltaX [mm];#DeltaY [mm];Entries", E),
                600, -30, 30, 600, -30, 30);
            h1_deltaX_hycal_E[ie] = make_h1(
                "h1_deltaX_hycal" + tag,
                Form("Delta X Between GEMs and HyCal for %d MeV;#DeltaX [mm];Entries", E),
                600, -30, 30);
            h1_deltaY_hycal_E[ie] = make_h1(
                "h1_deltaY_hycal" + tag,
                Form("Delta Y Between GEMs and HyCal for %d MeV;#DeltaY [mm];Entries", E),
                600, -30, 30);
            h2_deltaXY_hycal_E[ie] = make_h2(
                "h2_deltaXY_hycal" + tag,
                Form("Delta X vs Delta Y Between GEMs and HyCal for %d MeV;#DeltaX [mm];#DeltaY [mm];Entries", E),
                600, -30, 30, 600, -30, 30);

            for (int i = 0; i < kNGem; ++i) {
                h1_Nhits_matched_E[ie][i] = make_h1(
                    Form("h1_Nhits_matched_%dMeV_%d", E, i),
                    Form("Number of Hits in Matching Radius for GEM %d at %d MeV;N_{hits};Entries", i, E),
                    30, 0, 30);
            }
        }
    }

    void mergeFrom(const HistSet &other)
    {
        h1_deltaX_gem->Add(other.h1_deltaX_gem.get());
        h1_deltaY_gem->Add(other.h1_deltaY_gem.get());
        h2_deltaXY_gem->Add(other.h2_deltaXY_gem.get());
        h1_deltaX_hycal->Add(other.h1_deltaX_hycal.get());
        h1_deltaY_hycal->Add(other.h1_deltaY_hycal.get());
        h2_deltaXY_hycal->Add(other.h2_deltaXY_hycal.get());
        for (int i = 0; i < kNGem; ++i) {
            h1_Nhits_matched[i]->Add(other.h1_Nhits_matched[i].get());
        }

        for (int ie = 0; ie < kNEnergy; ++ie) {
            h1_deltaX_gem_E[ie]->Add(other.h1_deltaX_gem_E[ie].get());
            h1_deltaY_gem_E[ie]->Add(other.h1_deltaY_gem_E[ie].get());
            h2_deltaXY_gem_E[ie]->Add(other.h2_deltaXY_gem_E[ie].get());
            h1_deltaX_hycal_E[ie]->Add(other.h1_deltaX_hycal_E[ie].get());
            h1_deltaY_hycal_E[ie]->Add(other.h1_deltaY_hycal_E[ie].get());
            h2_deltaXY_hycal_E[ie]->Add(other.h2_deltaXY_hycal_E[ie].get());
            for (int i = 0; i < kNGem; ++i) {
                h1_Nhits_matched_E[ie][i]->Add(other.h1_Nhits_matched_E[ie][i].get());
            }
        }
    }

    void writeTo(TFile *output_file)
    {
        output_file->cd();
        h1_deltaX_gem->Write("h1_deltaX_gem");
        h1_deltaY_gem->Write("h1_deltaY_gem");
        h2_deltaXY_gem->Write("h2_deltaXY_gem");
        h1_deltaX_hycal->Write("h1_deltaX_hycal");
        h1_deltaY_hycal->Write("h1_deltaY_hycal");
        h2_deltaXY_hycal->Write("h2_deltaXY_hycal");
        for (int i = 0; i < kNGem; ++i) {
            h1_Nhits_matched[i]->Write(Form("h1_Nhits_matched_%d", i));
        }

        output_file->mkdir("energy_bins");
        output_file->cd("energy_bins");

        const std::array<int, kNEnergy> energies = {100, 300, 500, 1000, 1500};
        for (int ie = kNEnergy - 1; ie >= 0; --ie) {
            const int E = energies[ie];
            h1_deltaX_gem_E[ie]->Write(Form("h1_deltaX_gem_%dMeV", E));
            h1_deltaY_gem_E[ie]->Write(Form("h1_deltaY_gem_%dMeV", E));
            h2_deltaXY_gem_E[ie]->Write(Form("h2_deltaXY_gem_%dMeV", E));
            h1_deltaX_hycal_E[ie]->Write(Form("h1_deltaX_hycal_%dMeV", E));
            h1_deltaY_hycal_E[ie]->Write(Form("h1_deltaY_hycal_%dMeV", E));
            h2_deltaXY_hycal_E[ie]->Write(Form("h2_deltaXY_hycal_%dMeV", E));
        }

        for (int i = 0; i < kNGem; ++i) {
            h1_Nhits_matched_E[4][i]->Write(Form("h1_Nhits_matched_1500MeV_%d", i));
            h1_Nhits_matched_E[3][i]->Write(Form("h1_Nhits_matched_1000MeV_%d", i));
            h1_Nhits_matched_E[2][i]->Write(Form("h1_Nhits_matched_500MeV_%d", i));
            h1_Nhits_matched_E[1][i]->Write(Form("h1_Nhits_matched_300MeV_%d", i));
            h1_Nhits_matched_E[0][i]->Write(Form("h1_Nhits_matched_100MeV_%d", i));
        }
    }
};

static bool inEnergyWindow(float energy, float center)
{
    return std::fabs(energy - center) < 3.f * 0.033f * std::sqrt(center * 1000.f);
}

static void processTree(
    TTree *tree,
    EventVars_Recon &ev,
    HistSet &h,
    std::atomic<long long> &global_events,
    long long max_events,
    std::mutex &io_mtx,
    const std::string &label)
{
    const Long64_t nentries = tree->GetEntries();

    for (Long64_t ie = 0; ie < nentries; ++ie) {
        if (max_events > 0) {
            const long long ticket = global_events.fetch_add(1, std::memory_order_relaxed);
            if (ticket >= max_events) break;
        }

        if (ie % 10000 == 0) {
            std::lock_guard<std::mutex> lk(io_mtx);
            std::cerr << "[" << label << "] Processing event " << ie << " / " << nentries << "\r" << std::flush;
        }

        tree->GetEntry(ie);

        const bool is_3cluster = (ev.trigger_bits & prad2::TBIT_3cl) != 0;
        const bool is_sum = (ev.trigger_bits & prad2::TBIT_sum) != 0;
        if (!is_3cluster && !is_sum) continue;

        if (is_sum) {
            if (ev.n_clusters == 1 && ev.matchNum == 1 && ev.cl_nblocks[0] > 1 &&
                (ev.cl_energy[0] - gRunConfig.Ebeam) < 3.f * 0.033f * std::sqrt(gRunConfig.Ebeam * 1000.f)) {
                HCHit hc_hit = {ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]};
                GEMHit gem_hit_match[2];
                gem_hit_match[0] = {ev.mHit_gx[0][1], ev.mHit_gy[0][1], ev.mHit_gz[0][1], ev.mHit_gid[0][1]};
                gem_hit_match[1] = {ev.mHit_gx[0][0], ev.mHit_gy[0][0], ev.mHit_gz[0][0], ev.mHit_gid[0][0]};

                GetProjection(gem_hit_match[0], gem_hit_match[1].z);
                const float dx_gem = gem_hit_match[1].x - gem_hit_match[0].x;
                const float dy_gem = gem_hit_match[1].y - gem_hit_match[0].y;
                h.h1_deltaX_gem->Fill(dx_gem);
                h.h1_deltaY_gem->Fill(dy_gem);
                h.h2_deltaXY_gem->Fill(dx_gem, dy_gem);

                GetProjection(gem_hit_match[0], hc_hit.z);
                const float dx_hc = hc_hit.x - gem_hit_match[0].x;
                const float dy_hc = hc_hit.y - gem_hit_match[0].y;
                h.h1_deltaX_hycal->Fill(dx_hc);
                h.h1_deltaY_hycal->Fill(dy_hc);
                h.h2_deltaXY_hycal->Fill(dx_hc, dy_hc);

                int N_matched[4] = {0, 0, 0, 0};
                const size_t n_match = ev.match_cl_idx.size();
                for (size_t im = 0; im < n_match; ++im) {
                    const int det = static_cast<int>(ev.match_det_id[im]);
                    if (det >= 0 && det < 4) {
                        N_matched[det]++;
                    }
                }
                for (int d = 0; d < 4; ++d) {
                    if (N_matched[d] > 0) h.h1_Nhits_matched[d]->Fill(N_matched[d]);
                }
            }
        }

        if (is_3cluster && ev.matchNum == 3) {
            float x[3], y[3], z[3], E[3];
            for (int k = 0; k < 3; ++k) {
                x[k] = ev.mHit_x[k];
                y[k] = ev.mHit_y[k];
                z[k] = ev.mHit_z[k];
                E[k] = ev.mHit_E[k];
            }

            auto get_pt = [](float x0, float y0, float z0, float energy) {
                constexpr float electron_mass = 0.51099895f;
                const float norm = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
                if (norm <= 0.f || energy < electron_mass) {
                    return std::pair<float, float>{0.f, 0.f};
                }
                const float p = std::sqrt(std::max(0.f, energy * energy - electron_mass * electron_mass));
                return std::pair<float, float>{p * x0 / norm, p * y0 / norm};
            };

            const auto [px1, py1] = get_pt(x[0], y[0], z[0], E[0]);
            const auto [px2, py2] = get_pt(x[1], y[1], z[1], E[1]);
            const auto [px3, py3] = get_pt(x[2], y[2], z[2], E[2]);
            const float ptx = px1 + px2 + px3;
            const float pty = py1 + py2 + py3;

            if (std::fabs(E[0] + E[1] + E[2] - gRunConfig.Ebeam) > 250.f || std::sqrt(ptx * ptx + pty * pty) > 5.f) continue;
            if (E[0] < 70.f || E[1] < 70.f || E[2] < 70.f || E[0] > 1800.f || E[1] > 1800.f || E[2] > 1800.f) continue;

            for (int j = 0; j < ev.matchNum; ++j) {
                const int cl_idx = static_cast<int>(ev.mHit_cl_index[j]);
                if (cl_idx < 0 || cl_idx >= ev.n_clusters || cl_idx >= prad2::kMaxClusters) continue;
                if (ev.cl_nblocks[cl_idx] <= 1) continue;

                HCHit hc_hit = {ev.cl_x[cl_idx], ev.cl_y[cl_idx], ev.cl_z[cl_idx], ev.cl_energy[cl_idx]};
                GEMHit gem_hit_match[2];
                gem_hit_match[0] = {ev.mHit_gx[j][1], ev.mHit_gy[j][1], ev.mHit_gz[j][1], ev.mHit_gid[j][1]};
                gem_hit_match[1] = {ev.mHit_gx[j][0], ev.mHit_gy[j][0], ev.mHit_gz[j][0], ev.mHit_gid[j][0]};

                GetProjection(gem_hit_match[0], gem_hit_match[1].z);
                const float dx_gem = gem_hit_match[1].x - gem_hit_match[0].x;
                const float dy_gem = gem_hit_match[1].y - gem_hit_match[0].y;
                GetProjection(gem_hit_match[0], hc_hit.z);
                const float dx_hc = hc_hit.x - gem_hit_match[0].x;
                const float dy_hc = hc_hit.y - gem_hit_match[0].y;

                int N_matched[4] = {0, 0, 0, 0};
                const size_t n_match = ev.match_cl_idx.size();
                for (size_t im = 0; im < n_match; ++im) {
                    const int det = static_cast<int>(ev.match_det_id[im]);
                    const int cl = static_cast<int>(ev.match_cl_idx[im]);
                    if (det >= 0 && det < 4 && cl == cl_idx) {
                        N_matched[det]++;
                    }
                }

                const std::array<float, HistSet::kNEnergy> centers = {100.f, 300.f, 500.f, 1000.f, 1500.f};
                for (int ieh = 0; ieh < HistSet::kNEnergy; ++ieh) {
                    if (!inEnergyWindow(ev.cl_energy[cl_idx], centers[ieh])) continue;
                    h.h1_deltaX_gem_E[ieh]->Fill(dx_gem);
                    h.h1_deltaY_gem_E[ieh]->Fill(dy_gem);
                    h.h2_deltaXY_gem_E[ieh]->Fill(dx_gem, dy_gem);
                    h.h1_deltaX_hycal_E[ieh]->Fill(dx_hc);
                    h.h1_deltaY_hycal_E[ieh]->Fill(dy_hc);
                    h.h2_deltaXY_hycal_E[ieh]->Fill(dx_hc, dy_hc);
                    for (int d = 0; d < 4; ++d) {
                        if (N_matched[d] > 0) h.h1_Nhits_matched_E[ieh][d]->Fill(N_matched[d]);
                    }
                }
            }
        }
    }
}

} // namespace

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    ROOT::EnableThreadSafety();
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");
    TClass::GetClass("TH1F");
    TClass::GetClass("TH2F");

    std::string output = "matching_result.root";

    int max_events = -1;
    int nfiles = -1;
    int num_threads = 4;

    int opt;
    while ((opt = getopt(argc, argv, "o:n:f:j:")) != -1) {
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
            case 'j':
                if (!optarg) {
                    std::cerr << "Option -j requires an argument.\n";
                    return 1;
                }
                num_threads = std::max(1, std::atoi(optarg));
                break;
            default:
                std::cerr << "Usage: GEM_matching <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-f nfiles] [-j threads]\n";
                return 1;
        }
    }

    std::vector<std::string> root_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }

    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: GEM_matching <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-f nfiles] [-j threads]\n";
        return 1;
    }

    if (nfiles > 0 && static_cast<int>(root_files.size()) > nfiles) {
        root_files.resize(nfiles);
    }

    num_threads = std::min(num_threads, static_cast<int>(root_files.size()));
    num_threads = std::max(1, num_threads);

    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    const int run_num = get_run_int(root_files[0]);
    gRunConfig = LoadRunConfig(dbDir + "/runinfo/general.json", run_num);

    std::cerr << "Processing " << root_files.size() << " file(s) with " << num_threads << " thread(s)\n";

    std::atomic<size_t> next_file{0};
    std::atomic<long long> global_events{0};
    std::mutex io_mtx;

    std::vector<std::unique_ptr<HistSet>> local_hists;
    local_hists.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        local_hists.push_back(std::make_unique<HistSet>(Form("_t%d", i)));
    }

    std::vector<std::thread> workers;
    workers.reserve(num_threads);

    for (int tid = 0; tid < num_threads; ++tid) {
        workers.emplace_back([&, tid]() {
            EventVars_Recon ev;

            while (true) {
                if (max_events > 0 && global_events.load(std::memory_order_relaxed) >= max_events) break;

                const size_t idx = next_file.fetch_add(1, std::memory_order_relaxed);
                if (idx >= root_files.size()) break;

                const std::string &file = root_files[idx];
                {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cerr << "[thread " << tid << "] Added file: " << file << "\n";
                }

                std::unique_ptr<TFile> in_file(TFile::Open(file.c_str(), "READ"));
                if (!in_file || in_file->IsZombie()) {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cerr << "[thread " << tid << "] Cannot open file: " << file << "\n";
                    continue;
                }

                TTree *tree = dynamic_cast<TTree *>(in_file->Get("recon"));
                if (!tree) {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cerr << "[thread " << tid << "] Cannot find TTree 'recon' in: " << file << "\n";
                    continue;
                }

                prad2::SetReconReadBranches(tree, ev);
                prad2::ReconMatchVectorBindings match_bindings;
                prad2::BindReconMatchVectorBranches(tree, ev, match_bindings);

                processTree(tree, ev, *local_hists[tid], global_events, max_events, io_mtx, Form("thread %d", tid));
            }
        });
    }

    for (auto &t : workers) {
        t.join();
    }

    std::cerr << "\nProcessed events: " << global_events.load() << "\n";

    HistSet merged("");
    for (int i = 0; i < num_threads; ++i) {
        merged.mergeFrom(*local_hists[i]);
    }

    std::unique_ptr<TFile> output_file(TFile::Open(output.c_str(), "RECREATE"));
    if (!output_file || output_file->IsZombie()) {
        std::cerr << "Cannot create output file: " << output << "\n";
        return 1;
    }

    merged.writeTo(output_file.get());
    output_file->Close();

    return 0;
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
                name.size() >= 5 && name.compare(name.size() - 5, 5, ".root") == 0) {
                files.push_back(entry.path().string());
            }
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}
