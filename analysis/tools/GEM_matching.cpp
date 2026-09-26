// A quick check tool to test the matching result between HyCal clusters and GEM hits in the replay output
// Usage:
//   GEM_matching <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-f nfiles] [-j threads]
//   -o  output ROOT file (default: matching_result.root)
//   -n  max events to process across all files (default: all)
//   -f  max number of input files to add (default: all)
//   -j  number of worker threads (default: 4)

#include "PhysicsTools.h"
#include "MatchingTools.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"
#include "ToolUtils.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TTree.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <memory>
#include <mutex>
#include <string>
#include <unistd.h>
#include <vector>

using namespace analysis;

using EventVars_Recon = prad2::ReconEventData;

namespace {

// GEM-GEM and GEM-HyCal residuals of one matched cluster (mm).
struct Residuals {
    float dx_gem, dy_gem, dx_hc, dy_hc;
};

struct HistSet {
    static constexpr int kNGem = 4;
    static constexpr int kNEnergy = 5;
    static constexpr std::array<int, kNEnergy> kEnergies = {100, 300, 500, 1000, 1500};  // MeV

    // Residuals and number of GEM hits within the matching radius, for all
    // clusters (E = 0) or for the clusters around E MeV.
    struct ResidualHists {
        std::unique_ptr<TH1F> h1_deltaX_gem;
        std::unique_ptr<TH1F> h1_deltaY_gem;
        std::unique_ptr<TH2F> h2_deltaXY_gem;
        std::unique_ptr<TH1F> h1_deltaX_hycal;
        std::unique_ptr<TH1F> h1_deltaY_hycal;
        std::unique_ptr<TH2F> h2_deltaXY_hycal;
        std::array<std::unique_ptr<TH1F>, kNGem> h1_Nhits_matched;

        void book(HistList &reg, int E, const std::string &suffix)
        {
            const std::string tag = E > 0 ? Form("_%dMeV", E) : "";
            const std::string for_E = E > 0 ? Form(" for %d MeV", E) : "";
            const std::string at_E = E > 0 ? Form(" at %d MeV", E) : "";
            auto make_h1 = [&](const std::string &name, const std::string &title) {
                return Book<TH1F>(reg, (name + tag + suffix).c_str(), title.c_str(), 600, -30., 30.);
            };
            auto make_h2 = [&](const std::string &name, const std::string &title) {
                return Book<TH2F>(reg, (name + tag + suffix).c_str(), title.c_str(), 600, -30., 30., 600, -30., 30.);
            };

            h1_deltaX_gem = make_h1("h1_deltaX_gem", "Delta X Between GEMs" + for_E + ";#DeltaX [mm];Entries");
            h1_deltaY_gem = make_h1("h1_deltaY_gem", "Delta Y Between GEMs" + for_E + ";#DeltaY [mm];Entries");
            h2_deltaXY_gem = make_h2("h2_deltaXY_gem", "Delta X vs Delta Y Between GEMs" + for_E + ";#DeltaX [mm];#DeltaY [mm];Entries");
            h1_deltaX_hycal = make_h1("h1_deltaX_hycal", "Delta X Between GEMs and HyCal" + for_E + ";#DeltaX [mm];Entries");
            h1_deltaY_hycal = make_h1("h1_deltaY_hycal", "Delta Y Between GEMs and HyCal" + for_E + ";#DeltaY [mm];Entries");
            h2_deltaXY_hycal = make_h2("h2_deltaXY_hycal", "Delta X vs Delta Y Between GEMs and HyCal" + for_E + ";#DeltaX [mm];#DeltaY [mm];Entries");
            for (int i = 0; i < kNGem; ++i) {
                h1_Nhits_matched[i] = Book<TH1F>(reg,
                    Form("h1_Nhits_matched%s_%d%s", tag.c_str(), i, suffix.c_str()),
                    Form("Number of Hits in Matching Radius for GEM %d%s;N_{hits};Entries", i, at_E.c_str()),
                    30, 0., 30.);
            }
        }

        void fill(const Residuals &d, const std::array<int, kNGem> &n_matched)
        {
            h1_deltaX_gem->Fill(d.dx_gem);
            h1_deltaY_gem->Fill(d.dy_gem);
            h2_deltaXY_gem->Fill(d.dx_gem, d.dy_gem);
            h1_deltaX_hycal->Fill(d.dx_hc);
            h1_deltaY_hycal->Fill(d.dy_hc);
            h2_deltaXY_hycal->Fill(d.dx_hc, d.dy_hc);
            for (int i = 0; i < kNGem; ++i) {
                if (n_matched[i] > 0) h1_Nhits_matched[i]->Fill(n_matched[i]);
            }
        }

        void writeResiduals() const
        {
            h1_deltaX_gem->Write();
            h1_deltaY_gem->Write();
            h2_deltaXY_gem->Write();
            h1_deltaX_hycal->Write();
            h1_deltaY_hycal->Write();
            h2_deltaXY_hycal->Write();
        }
    };

    HistList all;
    ResidualHists total;
    std::array<ResidualHists, kNEnergy> byE;

    explicit HistSet(const std::string &suffix)
    {
        total.book(all, 0, suffix);
        for (int ie = 0; ie < kNEnergy; ++ie) byE[ie].book(all, kEnergies[ie], suffix);
    }

    void mergeFrom(const HistSet &other) { AddAll(all, other.all); }

    // Keys are the histogram names, so write the set booked with no suffix.
    void writeTo(TFile *output_file) const
    {
        output_file->cd();
        total.writeResiduals();
        for (const auto &h : total.h1_Nhits_matched) h->Write();

        output_file->mkdir("energy_bins");
        output_file->cd("energy_bins");
        for (int ie = kNEnergy - 1; ie >= 0; --ie) byE[ie].writeResiduals();
        for (int i = 0; i < kNGem; ++i) {
            for (int ie = kNEnergy - 1; ie >= 0; --ie) byE[ie].h1_Nhits_matched[i]->Write();
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
    // Matched cluster j (HyCal cluster cl_idx): its upstream GEM hit (mHit_g*[j][1])
    // projected onto its downstream GEM hit ([j][0]) and onto the cluster.
    const auto residuals = [&ev](int j, int cl_idx) {
        GEMHit gem = {ev.mHit_gx[j][1], ev.mHit_gy[j][1], ev.mHit_gz[j][1]};
        Residuals d;
        GetProjection(gem, ev.mHit_gz[j][0]);
        d.dx_gem = ev.mHit_gx[j][0] - gem.x;
        d.dy_gem = ev.mHit_gy[j][0] - gem.y;
        GetProjection(gem, ev.cl_z[cl_idx]);
        d.dx_hc = ev.cl_x[cl_idx] - gem.x;
        d.dy_hc = ev.cl_y[cl_idx] - gem.y;
        return d;
    };

    // Number of stored GEM matches of cluster cl_idx, per GEM.
    const auto count_matched = [&ev](int cl_idx) {
        std::array<int, HistSet::kNGem> n{};
        for (size_t im = 0; im < ev.match_cl_idx.size(); ++im) {
            const int det = static_cast<int>(ev.match_det_id[im]);
            if (det >= 0 && det < HistSet::kNGem && static_cast<int>(ev.match_cl_idx[im]) == cl_idx) {
                n[det]++;
            }
        }
        return n;
    };

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

        if (is_sum && ev.n_clusters == 1 && ev.matchNum == 1 && ev.cl_nblocks[0] > 1 &&
            (ev.cl_energy[0] - gRunConfig.Ebeam) < 3.f * 0.033f * std::sqrt(gRunConfig.Ebeam * 1000.f)) {
            h.total.fill(residuals(0, 0), count_matched(0));
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
                const float norm = std::sqrt(x0 * x0 + y0 * y0 + z0 * z0);
                if (norm <= 0.f || energy < PhysicsTools::kElectronMass) {
                    return std::pair<float, float>{0.f, 0.f};
                }
                const float p = std::sqrt(std::max(0.f, energy * energy - PhysicsTools::kElectronMass * PhysicsTools::kElectronMass));
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

                const Residuals d = residuals(j, cl_idx);
                const auto n_matched = count_matched(cl_idx);
                for (int ieh = 0; ieh < HistSet::kNEnergy; ++ieh) {
                    if (inEnergyWindow(ev.cl_energy[cl_idx], HistSet::kEnergies[ieh])) h.byE[ieh].fill(d, n_matched);
                }
            }
        }
    }
}

} // namespace

int main(int argc, char *argv[])
{
    InitRootThreading();

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

    std::vector<std::string> root_files = CollectInputs(argc, argv, optind, IsReconRootName, nfiles);

    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: GEM_matching <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-f nfiles] [-j threads]\n";
        return 1;
    }

    num_threads = std::min(num_threads, static_cast<int>(root_files.size()));
    num_threads = std::max(1, num_threads);

    std::string dbDir = prad2::database_dir();

    const int run_num = get_run_int(root_files[0]);
    gRunConfig = LoadRunConfig(dbDir + "/runinfo/general.json", run_num);

    std::cerr << "Processing " << root_files.size() << " file(s) with " << num_threads << " thread(s)\n";

    std::atomic<long long> global_events{0};
    std::mutex io_mtx;

    std::vector<EventVars_Recon> events(num_threads);
    std::vector<std::unique_ptr<HistSet>> local_hists;
    local_hists.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i) {
        local_hists.push_back(std::make_unique<HistSet>(Form("_t%d", i)));
    }

    ParallelFor(root_files.size(), num_threads, [&](size_t idx, int tid) {
        if (max_events > 0 && global_events.load(std::memory_order_relaxed) >= max_events) return;

        const std::string &file = root_files[idx];
        {
            std::lock_guard<std::mutex> lk(io_mtx);
            std::cerr << "[thread " << tid << "] Added file: " << file << "\n";
        }

        std::unique_ptr<TFile> in_file(TFile::Open(file.c_str(), "READ"));
        if (!in_file || in_file->IsZombie()) {
            std::lock_guard<std::mutex> lk(io_mtx);
            std::cerr << "[thread " << tid << "] Cannot open file: " << file << "\n";
            return;
        }

        TTree *tree = dynamic_cast<TTree *>(in_file->Get("recon"));
        if (!tree) {
            std::lock_guard<std::mutex> lk(io_mtx);
            std::cerr << "[thread " << tid << "] Cannot find TTree 'recon' in: " << file << "\n";
            return;
        }

        EventVars_Recon &ev = events[tid];
        prad2::SetReconReadBranches(tree, ev);
        prad2::ReconMatchVectorBindings match_bindings;
        prad2::BindReconMatchVectorBranches(tree, ev, match_bindings);

        processTree(tree, ev, *local_hists[tid], global_events, max_events, io_mtx, Form("thread %d", tid));
    });

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
