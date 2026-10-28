// quick_check.cpp — quick physics check of replayed recon files
//
// Reads reconstructed ROOT tree (output of replay_recon), runs physics
// analysis using PhysicsTools, and saves histograms to an output ROOT file.
// Usage:
//   quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-j threads]
//   -o  output ROOT file (default: input filename with _quick_check.root suffix)
//   -n  max events to process (default: all)
//   -j  number of input-file worker threads (default: 4)
// Example:
//   quick_check recon.root -o recon_check.root -n 10000
//   quick_check recon_dir/ recon.root...  -n 100000

#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"
#include "ToolUtils.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
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
#include <memory>
#include <mutex>
#include <limits>
#include <thread>
#include <unistd.h>

using namespace analysis;
namespace fs = std::filesystem;

using EventVars_Recon = prad2::ReconEventData;

static std::string makeDefaultOutput(const std::string &input_path);

static bool inHyCal(float x, float y) { return InHyCalRing(x, y, 2.2, 16.); }

static float electronPairInvariantMass(
    float x1, float y1, float z1, float E1,
    float x2, float y2, float z2, float E2)
{
    constexpr float electron_mass = PhysicsTools::kElectronMass;
    const float position[2][3] = {{x1, y1, z1}, {x2, y2, z2}};
    const float energy[2] = {E1, E2};
    float momentum[2][3] = {};

    for (int i = 0; i < 2; ++i) {
        const float norm = std::sqrt(
            position[i][0] * position[i][0]
            + position[i][1] * position[i][1]
            + position[i][2] * position[i][2]);
        if (energy[i] < electron_mass || norm <= 0.f)
            return std::numeric_limits<float>::quiet_NaN();
        const float p = std::sqrt(std::max(
            0.f, energy[i] * energy[i] - electron_mass * electron_mass));
        momentum[i][0] = p * position[i][0] / norm;
        momentum[i][1] = p * position[i][1] / norm;
        momentum[i][2] = p * position[i][2] / norm;
    }

    const float total_energy = energy[0] + energy[1];
    const float px = momentum[0][0] + momentum[1][0];
    const float py = momentum[0][1] + momentum[1][1];
    const float pz = momentum[0][2] + momentum[1][2];
    const float mass2 = total_energy * total_energy - px * px - py * py - pz * pz;
    return std::sqrt(std::max(0.f, mass2));
}

const int Nbins = 33;
const float binEdge[Nbins+1] = {
    0.500, 0.550, 0.600, 0.650, 0.700, 0.750, 0.775, 0.800, 0.825, 0.850,
    0.875, 0.900, 0.940, 0.975, 1.014, 1.057, 1.105, 1.157, 1.211, 1.270,
    1.338, 1.417, 1.514, 1.634, 1.787, 2.000, 2.213, 2.492, 2.792, 3.092,
    3.392, 3.692, 3.992, 4.292
};

struct QuickResult {
    std::unique_ptr<PhysicsTools> physics;
    HistList all;   // every histogram below except h_ep_ee_ratio, merged by AddAll
    std::unique_ptr<TH2F> hit_pos;
    std::unique_ptr<TH1F> h_1cl;
    std::unique_ptr<TH1F> h_2cl;
    std::unique_ptr<TH1F> h_all;
    std::unique_ptr<TH1F> h_tot;
    std::unique_ptr<TH2F> h2_energy_theta_ep_ee;
    MollerData mollers;
    MollerData mollers_hc;
    Long64_t processed = 0;

    std::unique_ptr<TH2F> h2_ep_hits;
    std::unique_ptr<TH2F> h2_ee_hits;
    std::unique_ptr<TH2F> h2_ep_E_angle;
    std::unique_ptr<TH2F> h2_ee_E_angle;

    std::unique_ptr<TH1F> h_ep_yield;
    std::unique_ptr<TH1F> h_ee_yield;
    std::unique_ptr<TH1F> h_ep_ee_ratio;
    std::unique_ptr<TH1F> h_ee_tDiff;

    std::unique_ptr<TH1F> h_ee_center_x;
    std::unique_ptr<TH1F> h_ee_center_y;
    std::unique_ptr<TH1F> h_ee_vertex_z;

    std::unique_ptr<TH2F> h2_ep_hits_hc;
    std::unique_ptr<TH2F> h2_ee_hits_hc;
    std::unique_ptr<TH2F> h2_ep_E_angle_hc;
    std::unique_ptr<TH2F> h2_ee_E_angle_hc;

    std::unique_ptr<TH1F> h_ee_center_x_hc;
    std::unique_ptr<TH1F> h_ee_center_y_hc;
    std::unique_ptr<TH1F> h_ee_vertex_z_hc;

    std::unique_ptr<TH1F> h_ee_invariant_mass;

    // For X17
    // gamma decay channel (e gamma gamma)
    // step by step cuts histograms
    std::unique_ptr<TH1F> h_gamma_totalE[5];
    std::unique_ptr<TH2F> h2_gamma_hits[5];
    std::unique_ptr<TH1F> h_gamma_E[5];
    std::unique_ptr<TH1F> h_gamma_E_gamma[5];
    std::unique_ptr<TH1F> h_gamma_E_electron[5];
    std::unique_ptr<TH2F> h2_gamma_E_gamma_vs_E_electron[5];
    std::unique_ptr<TH2F> h2_gamma_E_gamma_vs_E_gamma[5];
    std::unique_ptr<TH2F> h2_gamma_E_angle_gamma[5];
    std::unique_ptr<TH2F> h2_gamma_E_angle_electron[5];
    std::unique_ptr<TH1F> h_gamma_ptx[5];
    std::unique_ptr<TH1F> h_gamma_pty[5];
    std::unique_ptr<TH2F> h2_gamma_Pt[5];
    std::unique_ptr<TH1F> h_gamma_tDiff[5];
    std::unique_ptr<TH1F> h_gamma_dphi[5];
    std::unique_ptr<TH1F> h_gamma_mass[5];

    // e+/e- decay channel (3 clusters)
    std::unique_ptr<TH1F> h_3cl_cluster_num;
    std::unique_ptr<TH1F> h_3cl_cluster_num_cut_cl;
    std::unique_ptr<TH1F> h_3cl_tDiff_raw;
    std::unique_ptr<TH1F> h_3cl_cluster_num_cut_cl_t;
    // step by step cuts histograms
    std::unique_ptr<TH1F> h_3cl_totalE[6];
    std::unique_ptr<TH2F> h2_3cl_hits[6];
    std::unique_ptr<TH1F> h_3cl_E[6];
    std::unique_ptr<TH1F> h_3cl_yield[6];
    std::unique_ptr<TH1F> h_3cl_ptx[6];
    std::unique_ptr<TH1F> h_3cl_pty[6];
    std::unique_ptr<TH2F> h2_3cl_Pt[6];
    std::unique_ptr<TH1F> h_3cl_tDiff[6];
    std::unique_ptr<TH1F> h_3cl_dphi[6];
    std::unique_ptr<TH1F> h_3cl_vertexZ[6];
    std::unique_ptr<TH2F> h2_3cl_E_angle[6];
    std::unique_ptr<TH1F> h_3cl_mass[6];
    std::unique_ptr<TH1F> h_3cl_mass_1comb[6];
    std::unique_ptr<TH1F> h_3cl_mass_2comb[6];
};

static std::unique_ptr<QuickResult> makeResult(fdec::HyCalSystem &hycal)
{
    auto r = std::make_unique<QuickResult>();
    r->physics = std::make_unique<PhysicsTools>(hycal);
    r->hit_pos = Book<TH2F>(r->all, "hit_pos",
        "Hit positions;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h_1cl = Book<TH1F>(r->all, "one_cluster_energy",
        "Single-cluster energy;E (MeV);Counts", 4000, 0, 4000);
    r->h_2cl = Book<TH1F>(r->all, "two_cluster_energy",
        "Two-cluster energy;E (MeV);Counts", 4000, 0, 4000);
    r->h_all = Book<TH1F>(r->all, "clusters_energy",
        "All clusters;E (MeV);Counts", 4000, 0, 4000);
    r->h_tot = Book<TH1F>(r->all, "total_energy",
        "Total energy per event;E (MeV);Counts", 4000, 0, 4000);
    r->h2_energy_theta_ep_ee = Book<TH2F>(r->all, "energy_vs_theta",
        "Energy vs Theta(1 cluster);Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);

    r->h2_ep_hits = Book<TH2F>(r->all, "ep_hits",
        "EP Hit positions;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ee_hits = Book<TH2F>(r->all, "ee_hits",
        "EE Hit positions;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ep_E_angle = Book<TH2F>(r->all, "ep_E_angle",
        "EP Energy vs Angle;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);
    r->h2_ee_E_angle = Book<TH2F>(r->all, "ee_E_angle",
        "EE Energy vs Angle;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);

    r->h_ep_yield = Book<TH1F>(r->all, "ep_yield",
        "EP Yield;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_ee_yield = Book<TH1F>(r->all, "ee_yield",
        "EE Yield;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_ep_ee_ratio = std::make_unique<TH1F>("ep_ee_ratio",
        "EP/EE Yield Ratio;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_ee_tDiff = Book<TH1F>(r->all, "ee_tDiff",
        "EE Time Difference;Time Difference (ns);Counts", 400, -10, 10);

    r->h_ee_center_x = Book<TH1F>(r->all, "ee_center_x",
        "EE Center X;X (mm);Counts", 800, -20, 20);
    r->h_ee_center_y = Book<TH1F>(r->all, "ee_center_y",
        "EE Center Y;Y (mm);Counts", 800, -20, 20);
    r->h_ee_vertex_z = Book<TH1F>(r->all, "ee_vertex_z",
        "EE Vertex Z;Z (mm);Counts", 8000, 5000, 9000);

    r->h2_ep_hits_hc = Book<TH2F>(r->all, "ep_hits_hc",
        "EP Hit positions hycal;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ee_hits_hc = Book<TH2F>(r->all, "ee_hits_hc",
        "EE Hit positions hycal;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ep_E_angle_hc = Book<TH2F>(r->all, "ep_E_angle_hc",
        "EP Energy vs Angle hycal;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);
    r->h2_ee_E_angle_hc = Book<TH2F>(r->all, "ee_E_angle_hc",
        "EE Energy vs Angle hycal;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);

    r->h_ee_center_x_hc = Book<TH1F>(r->all, "ee_center_x_hc",
        "EE Center X hycal;X (mm);Counts", 800, -20, 20);
    r->h_ee_center_y_hc = Book<TH1F>(r->all, "ee_center_y_hc",
        "EE Center Y hycal;Y (mm);Counts", 800, -20, 20);
    r->h_ee_vertex_z_hc = Book<TH1F>(r->all, "ee_vertex_z_hc",
        "EE Vertex Z hycal;Z (mm);Counts", 8000, 5000, 9000);

    r->h_ee_invariant_mass = Book<TH1F>(r->all, "ee_invariant_mass",
        "EE Invariant Mass;Mass (MeV);Counts", 400, 0, 100);

    // X17 gamma decay channel
    for (int i = 0; i < 5; i++) {
        r->h_gamma_totalE[i] = Book<TH1F>(r->all, Form("gamma_totalE_step_%d", i),
            Form("Gamma Channel Total Energy - Step %d;Total Energy (MeV);Counts", i), 2500, 0, 2500);
        r->h2_gamma_hits[i] = Book<TH2F>(r->all, Form("gamma_hits_step_%d", i),
            Form("Gamma Channel Hit positions hycal - Step %d;X (mm);Y (mm)", i), 720, -360, 360, 720, -360, 360);
        r->h_gamma_E[i] = Book<TH1F>(r->all, Form("gamma_E_step_%d", i),
            Form("Gamma Channel Energy - Step %d;Energy (MeV);Counts", i), 2500, 0, 2500);
        r->h_gamma_E_gamma[i] = Book<TH1F>(r->all, Form("gamma_E_gamma_step_%d", i),
            Form("Gamma Channel Energy Gamma - Step %d;Energy (MeV);Counts", i), 2500, 0, 2500);
        r->h_gamma_E_electron[i] = Book<TH1F>(r->all, Form("gamma_E_electron_step_%d", i),
            Form("Gamma Channel Energy Electron - Step %d;Energy (MeV);Counts", i), 2500, 0, 2500);
        r->h2_gamma_E_gamma_vs_E_electron[i] = Book<TH2F>(r->all, Form("gamma_E_gamma_vs_E_electron_step_%d", i),
            Form("Gamma Channel E_gamma vs E_electron - Step %d;E_electron (MeV);E_gamma (MeV)", i), 2500, 0, 2500, 2500, 0, 2500);
        r->h2_gamma_E_gamma_vs_E_gamma[i] = Book<TH2F>(r->all, Form("gamma_E_gamma_vs_E_gamma_step_%d", i),
            Form("Gamma Channel E_gamma vs E_gamma - Step %d;E_gamma (MeV);E_gamma (MeV)", i), 2500, 0, 2500, 2500, 0, 2500);
        r->h2_gamma_E_angle_gamma[i] = Book<TH2F>(r->all, Form("gamma_E_angle_gamma_step_%d", i),
            Form("Gamma Channel E_gamma vs Angle Gamma - Step %d;Theta (deg);E_gamma (MeV)", i), 80, 0, 4, 2500, 0, 2500);
        r->h2_gamma_E_angle_electron[i] = Book<TH2F>(r->all, Form("gamma_E_angle_electron_step_%d", i),
            Form("Gamma Channel E_electron vs Angle Electron - Step %d;Theta (deg);E_electron (MeV)", i), 80, 0, 4, 2500, 0, 2500);
        r->h_gamma_ptx[i] = Book<TH1F>(r->all, Form("gamma_ptx_step_%d", i),
            Form("Gamma Channel Ptx - Step %d;Ptx (MeV/c);Counts", i), 200, -50, 50);
        r->h_gamma_pty[i] = Book<TH1F>(r->all, Form("gamma_pty_step_%d", i),
            Form("Gamma Channel Pty - Step %d;Pty (MeV/c);Counts", i), 200, -50, 50);
        r->h2_gamma_Pt[i] = Book<TH2F>(r->all, Form("gamma_Pt_step_%d", i),
            Form("Gamma Channel Pt - Step %d;Ptx (MeV/c);Pty (MeV/c)", i), 400, -50, 50, 400, -50, 50);
        r->h_gamma_tDiff[i] = Book<TH1F>(r->all, Form("gamma_tDiff_step_%d", i),
            Form("Gamma Channel Time Difference - Step %d;#Deltat (ns);Counts", i), 240*2, -12, 12);
        r->h_gamma_dphi[i] = Book<TH1F>(r->all, Form("gamma_dphi_step_%d", i),
            Form("Gamma Channel Delta Phi - Step %d;#Delta#phi (rad);Counts", i), 360*3, 0, 360);
        r->h_gamma_mass[i] = Book<TH1F>(r->all, Form("gamma_mass_step_%d", i),
            Form("Gamma Channel Invariant Mass - Step %d;Mass (MeV/c^{2});Counts", i), 1000, 0, 100);
    }

    // use gem matching to cut the 3-cluster events
    r->h_3cl_cluster_num = Book<TH1F>(r->all, "3cl_cluster_gem_num",
        "GEM-matched Cluster Number;Number of Clusters;Counts", 20, 0, 20);
    r->h_3cl_cluster_num_cut_cl = Book<TH1F>(r->all, "3cl_cluster_num_cut_cl",
        "Candidate Number after Cluster Quality Cuts;Number of Clusters;Counts", 20, 0, 20);
    r->h_3cl_tDiff_raw = Book<TH1F>(r->all, "3cl_tDiff_raw",
        "Time Difference to Leading-E Cluster;#Deltat (ns);Counts", 400, -10, 10);
    r->h_3cl_cluster_num_cut_cl_t = Book<TH1F>(r->all, "3cl_cluster_num_cut_cl_t",
        "Candidate Number after Cluster+Timing Cuts;Number of Clusters;Counts", 20, 0, 20);
    // step by step cuts histograms
    for(int i = 0; i < 6; i++){
        r->h2_3cl_hits[i] = Book<TH2F>(r->all, Form("3cl_hits_step_%d", i),
            Form("3-Cluster Hit positions on hycal - Step %d;X (mm);Y (mm)", i), 720, -360, 360, 720, -360, 360);
        r->h2_3cl_E_angle[i] = Book<TH2F>(r->all, Form("3cl_E_angle_step_%d", i),
            Form("3-Cluster Energy vs Angle- Step %d;Theta (deg);Energy (MeV)", i), 80, 0, 4, 2500, 0, 2500);
        r->h_3cl_E[i] = Book<TH1F>(r->all, Form("3cl_E_step_%d", i),
            Form("3-Cluster Energy - Step %d;Energy (MeV);Counts", i), 2500, 0, 2500);
        r->h_3cl_totalE[i] = Book<TH1F>(r->all, Form("3cl_totalE_step_%d", i),
            Form("3-Cluster Total Energy - Step %d;Total Energy (MeV);Counts", i), 2500, 0, 2500);
        r->h_3cl_yield[i] = Book<TH1F>(r->all, Form("3cl_yield_step_%d", i),
            Form("3-Cluster Yield - Step %d;Scattering Angle (deg);Counts", i), Nbins, binEdge);
        r->h_3cl_ptx[i] = Book<TH1F>(r->all, Form("3cl_ptx_step_%d", i),
            Form("3-Cluster Ptx - Step %d;Ptx (MeV);Counts", i), 200, -50, 50);
        r->h_3cl_pty[i] = Book<TH1F>(r->all, Form("3cl_pty_step_%d", i),
            Form("3-Cluster Pty - Step %d;Pty (MeV);Counts", i), 200, -50, 50);
        r->h_3cl_tDiff[i] = Book<TH1F>(r->all, Form("3cl_tDiff_step_%d", i),
            Form("3-Cluster Time Difference - Step %d;Time Difference (ns);Counts", i), 240*2, -12, 12);
        r->h_3cl_dphi[i] = Book<TH1F>(r->all, Form("3cl_dphi_step_%d", i),
            Form("3-Cluster Phi Difference - Step %d;#Delta#phi (deg);Counts", i), 360*3, 0, 360);
        r->h2_3cl_Pt[i] = Book<TH2F>(r->all, Form("3cl_Pt_step_%d", i),
            Form("3-Cluster Pt hycal - Step %d;Ptx (MeV);Pty (MeV);Counts", i), 400, -50, 50, 400, -50, 50);
        r->h_3cl_vertexZ[i] = Book<TH1F>(r->all, Form("3cl_VertexZ_step_%d", i),
            Form("3-Cluster Vertex Z - Step %d;Vertex Z (mm);Counts", i), 1100, -3500, 7500);
        r->h_3cl_mass[i] = Book<TH1F>(r->all, Form("3cl_mass_step_%d", i),
            Form("3-Cluster Inv. Mass - Step %d;Inv. Mass (MeV);Counts", i), 1000, 0, 100);
        r->h_3cl_mass_1comb[i] = Book<TH1F>(r->all, Form("3cl_mass_1comb_step_%d", i),
            Form("3-Cluster Inv. Mass - Step %d, 1 Combination;Inv. Mass (MeV);Counts", i), 1000, 0, 100);
        r->h_3cl_mass_2comb[i] = Book<TH1F>(r->all, Form("3cl_mass_2comb_step_%d", i),
            Form("3-Cluster Inv. Mass - Step %d, 2 Combinations;Inv. Mass (MeV);Counts", i), 1000, 0, 100);
    }

    return r;
}

// Keeps the last three Moller pairs in buf and fills the centres that m forms
// with the earlier ones.
static void fillMollerCenters(MollerData &buf, const MollerEvent &m, TH1F *hx, TH1F *hy)
{
    buf.push_back(m);
    if (buf.size() > 3) buf.erase(buf.begin());
    for (size_t k = 2; k <= buf.size(); ++k) {
        const auto c = PhysicsTools::GetMollerCenter(buf[buf.size() - k], m);
        hx->Fill(c[0]);
        hy->Fill(c[1]);
    }
}

static bool processFile(const std::string &path,
                        Long64_t max_entries,
                        float Ebeam,
                        QuickResult &out)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) {
        std::cerr << "Cannot open " << path << "\n";
        return false;
    }
    TTree *tree = dynamic_cast<TTree *>(f->Get("recon"));
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in " << path << "\n";
        return false;
    }

    EventVars_Recon ev;
    prad2::SetReconReadBranches(tree, ev);
    Long64_t n = tree->GetEntries();
    if (max_entries >= 0 && max_entries < n) n = max_entries;

    auto &physics = *out.physics;
    for (Long64_t i = 0; i < n; i++) {
        tree->GetEntry(i);

        // trigger selection
        bool is_3cluster = (ev.trigger_bits & prad2::TBIT_3cl) != 0;
        bool is_sum      = (ev.trigger_bits & prad2::TBIT_sum) != 0;

        if (!is_3cluster && !is_sum) continue;

        if(is_sum){
            for (int j = 0; j < ev.n_clusters; j++) {
                float r = std::sqrt(ev.cl_x[j]*ev.cl_x[j] + ev.cl_y[j]*ev.cl_y[j]);
                float theta = std::atan(r / ev.cl_z[j]) * 180.f / M_PI;

                physics.FillEnergyVsModule(ev.cl_center[j], ev.cl_energy[j]);
                out.hit_pos->Fill(ev.cl_x[j], ev.cl_y[j]);
                out.h_all->Fill(ev.cl_energy[j]);

                if (ev.cl_nblocks[j] > 1 && inHyCal(ev.cl_x[j], ev.cl_y[j])) {
                    physics.FillEnergyVsTheta(theta, ev.cl_energy[j]);
                }
            }
            out.h_tot->Fill(ev.total_energy);

            if (ev.n_clusters == 1) {
                physics.FillModuleEnergy(ev.cl_center[0], ev.cl_energy[0]);
                out.h_1cl->Fill(ev.cl_energy[0]);
                out.h2_energy_theta_ep_ee->Fill(
                    std::atan(std::sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]) / ev.cl_z[0]) * 180.f / M_PI,
                    ev.cl_energy[0]);
                out.h2_ep_hits_hc->Fill(ev.cl_x[0], ev.cl_y[0]);
                out.h2_ep_E_angle_hc->Fill(
                    std::atan(std::sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]) / ev.cl_z[0]) * 180.f / M_PI,
                    ev.cl_energy[0]);
            }

            if (ev.n_clusters == 2 && inHyCal(ev.cl_x[0], ev.cl_y[0]) && inHyCal(ev.cl_x[1], ev.cl_y[1])) {
                out.h_2cl->Fill(ev.cl_energy[0]);
                out.h_2cl->Fill(ev.cl_energy[1]);

                float Epair = ev.cl_energy[0] + ev.cl_energy[1];
                float sigma = Ebeam * 0.033f / std::sqrt(Ebeam / 1000.f);
                if (std::abs(Epair - Ebeam) < 3. * sigma) {
                    MollerEvent mp(
                        {ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]},
                        {ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1]});
                    physics.FillMollerPhiDiff(physics.GetMollerPhiDiff(mp));
                    if (PhysicsTools::GetMollerPhiDiff(mp) < 10.f) {
                        out.h2_ee_hits_hc->Fill(ev.cl_x[0], ev.cl_y[0]);
                        out.h2_ee_hits_hc->Fill(ev.cl_x[1], ev.cl_y[1]);
                        float t1 = std::atan2(std::sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]), ev.cl_z[0]) * 180.f / M_PI;
                        float t2 = std::atan2(std::sqrt(ev.cl_x[1]*ev.cl_x[1] + ev.cl_y[1]*ev.cl_y[1]), ev.cl_z[1]) * 180.f / M_PI;
                        out.h2_ee_E_angle_hc->Fill(t1, ev.cl_energy[0]);
                        out.h2_ee_E_angle_hc->Fill(t2, ev.cl_energy[1]);
                        fillMollerCenters(out.mollers_hc, mp, out.h_ee_center_x_hc.get(), out.h_ee_center_y_hc.get());
                        float vertex = physics.GetMollerZdistance(mp, Ebeam);
                        out.h_ee_vertex_z_hc->Fill(vertex);

                        const float invariant_mass = electronPairInvariantMass(
                            ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0],
                            ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1]);
                        if (std::isfinite(invariant_mass))
                            out.h_ee_invariant_mass->Fill(invariant_mass);
                    }
                }
            }

            // loop over GEM matched hits find e-p events
            const int n_match = std::clamp(ev.matchNum, 0, prad2::kMaxClusters);
            for (int j = 0; j < n_match; ++j) {
                float x = ev.mHit_gx[j][1];
                float y = ev.mHit_gy[j][1];
                float z = ev.mHit_gz[j][1];
                float E = ev.mHit_E[j];
                float scale = ev.mHit_z[j] / z;
                x *= scale;
                y *= scale;
                z *= scale;

                if (!inHyCal(x, y)) continue;

                float theta = std::atan(std::sqrt(x*x + y*y) / z) * 180.f / M_PI;
                float expectE = physics.ExpectedEnergy(theta, Ebeam, "ep");
                if (fabs(E - expectE) < 3.f * expectE * 0.035f / std::sqrt(E/1000.f)) {
                    out.h2_ep_hits->Fill(x, y);
                    out.h2_ep_E_angle->Fill(theta, E);
                    out.h_ep_yield->Fill(theta);
                }
            }

            //select GEM matched Moller events 
            if (ev.matchNum == 2) {
                // 0 is downstream, 1 is upstream
                float x[2] = {ev.mHit_gx[0][1], ev.mHit_gx[1][1]};
                float y[2] = {ev.mHit_gy[0][1], ev.mHit_gy[1][1]};
                float z[2] = {ev.mHit_gz[0][1], ev.mHit_gz[1][1]};
                float E[2] = {ev.mHit_E[0], ev.mHit_E[1]};
                int idx[2] = {ev.mHit_cl_index[0], ev.mHit_cl_index[1]};
                if (idx[0] < 0 || idx[1] < 0
                    || idx[0] >= ev.n_clusters || idx[1] >= ev.n_clusters
                    || idx[0] >= prad2::kMaxClusters || idx[1] >= prad2::kMaxClusters)
                    continue;
                float time[2] = {ev.cl_time[idx[0]], ev.cl_time[idx[1]]};
                int mod_id[2] = {ev.cl_center[idx[0]], ev.cl_center[idx[1]]};
                float scale[2] = {ev.mHit_z[0] / z[0], ev.mHit_z[1] / z[1]};
                for (int j = 0; j < 2; j++) {
                    x[j] *= scale[j];
                    y[j] *= scale[j];
                    z[j] *= scale[j];
                }
                if (!inHyCal(x[0], y[0]) || !inHyCal(x[1], y[1])) continue;
                float theta[2] = {
                    std::atan(std::sqrt(x[0]*x[0] + y[0]*y[0]) / z[0]) * 180.f / static_cast<float>(M_PI),
                    std::atan(std::sqrt(x[1]*x[1] + y[1]*y[1]) / z[1]) * 180.f / static_cast<float>(M_PI)
                };

                MollerEvent mev({x[0], y[0], z[0], E[0]}, {x[1], y[1], z[1], E[1]});
                if (PhysicsTools::isMoller_kinematic(theta[0], E[0], theta[1], E[1], Ebeam, 0.033f)
                    && PhysicsTools::isBackToBack(mev, 10.f))
                {
                    out.h2_ee_hits->Fill(x[0], y[0]);
                    out.h2_ee_hits->Fill(x[1], y[1]);
                    out.h2_ee_E_angle->Fill(theta[0], E[0]);
                    out.h2_ee_E_angle->Fill(theta[1], E[1]);
                    out.h_ee_yield->Fill(theta[0]);
                    out.h_ee_yield->Fill(theta[1]);
                    float vertex = physics.GetMollerZdistance(mev, Ebeam);
                    out.h_ee_vertex_z->Fill(vertex);
                    float delta_time = time[0] - time[1];
                    if (mod_id[0] < mod_id[1]) delta_time = -delta_time;
                    out.h_ee_tDiff->Fill(delta_time);
                    fillMollerCenters(out.mollers, mev, out.h_ee_center_x.get(), out.h_ee_center_y.get());
                }
            }
        }

        // x17 trigger selection
        if (is_3cluster) {
            //try to find the gamma decay channel,
            //firstly try on the clean events(only 3 clusters on HyCal)
            if(ev.n_clusters == 3 && ev.matchNum == 1) {
                int e_idx = ev.mHit_cl_index[0];
                float E_e = ev.cl_energy[e_idx];
                float x_e = ev.mHit_gx[0][1];
                float y_e = ev.mHit_gy[0][1];
                float z_e = ev.mHit_gz[0][1];
                float t_e = ev.cl_time[e_idx];

                float scale = ev.cl_z[e_idx] / z_e;
                x_e *= scale;
                y_e *= scale;
                z_e *= scale;

                float E_g[2], x_g[2], y_g[2], z_g[2], t_g[2];
                int gamma_idx = 0;
                bool noMatch[2] = {true, true};
                for (int j = 0; j < 3; j++) {
                    if ( j == e_idx) continue;
                    E_g[gamma_idx] = ev.cl_energy[j];
                    x_g[gamma_idx] = ev.cl_x[j];
                    y_g[gamma_idx] = ev.cl_y[j];
                    z_g[gamma_idx] = ev.cl_z[j];
                    t_g[gamma_idx] = ev.cl_time[j];
                    if (ev.matchFlag[j] != 0) noMatch[gamma_idx] = false;
                    gamma_idx++;
                }
                bool bothNoMatch = noMatch[0] && noMatch[1];

                float Sigma_e = 0.035f * std::sqrt(E_e * 1000.f);
                float Sigma[2] = {0.035f * std::sqrt(E_g[0] * 1000.f),
                                  0.035f * std::sqrt(E_g[1] * 1000.f)};
                float totalSigma = std::sqrt(Sigma_e*Sigma_e + Sigma[0]*Sigma[0] + Sigma[1]*Sigma[1]);
                float totalE = E_e + E_g[0] + E_g[1];

                float theta_g[2], theta_e;
                theta_e = std::atan2(std::sqrt(x_e*x_e + y_e*y_e), z_e) * 180.f / static_cast<float>(M_PI);
                theta_g[0] = std::atan2(std::sqrt(x_g[0]*x_g[0] + y_g[0]*y_g[0]), z_g[0]) * 180.f / static_cast<float>(M_PI);
                theta_g[1] = std::atan2(std::sqrt(x_g[1]*x_g[1] + y_g[1]*y_g[1]), z_g[1]) * 180.f / static_cast<float>(M_PI);

                float tDiff = std::max({std::fabs(t_e - t_g[0]), std::fabs(t_e - t_g[1]), std::fabs(t_g[0] - t_g[1])});
                float dt[2] = {t_g[0] - t_e, t_g[1] - t_e};

                // 4-momentum calculation for each single hit and gamma pair hits
                TLorentzVector p_e, p_g[2], p_pair;
                PhysicsTools::HitP4(x_e, y_e, z_e, E_e, PhysicsTools::kElectronMass, p_e);
                for (int i = 0; i < 2; ++i)
                    PhysicsTools::HitP4(x_g[i], y_g[i], z_g[i], E_g[i], 0.f, p_g[i]);
                p_pair = p_g[0] + p_g[1];

                // Pt x and Pt y calculation using TLorentzVector
                float ptx = p_e.Px() + p_g[0].Px() + p_g[1].Px();
                float pty = p_e.Py() + p_g[0].Py() + p_g[1].Py();

                float mass = (p_g[0] + p_g[1]).M();

                // get the azimuthal angles for each single hit and the pair of gamma hits
                float phi_e = std::atan2(p_e.Py(), p_e.Px()) * 180.f / static_cast<float>(M_PI);
                float phi_pair = std::atan2(p_pair.Py(), p_pair.Px()) * 180.f / static_cast<float>(M_PI);

                // Phi difference for the pair of gamma hits and the electron hit
                float dphi = std::fabs(phi_pair - phi_e);

                bool nblocks_ok = (ev.cl_nblocks[0] > 1 && ev.cl_nblocks[1] > 1 && ev.cl_nblocks[2] > 1);
                bool totalE_pass = std::fabs(totalE - Ebeam) < 4. * totalSigma;
                bool Pt_pass = std::sqrt(ptx * ptx + pty * pty) < 5.0f;
                bool pos_pass = inHyCal(x_e, y_e) && inHyCal(x_g[0], y_g[0]) && inHyCal(x_g[1], y_g[1]);
                bool time_pass = tDiff < 2.f;
                bool clusterE_pass = E_e > 70.f && E_g[0] > 70.f && E_g[1] > 70.f && E_e < 0.75 * Ebeam && E_g[0] < 0.75 * Ebeam && E_g[1] < 0.75 * Ebeam;
                bool dphi_pass = std::fabs(dphi - 180.f) < 10.0f;

                // cut step s applies the cuts of steps 0..s
                const bool step_cut[5] = {
                    bothNoMatch && nblocks_ok && pos_pass && clusterE_pass,  // 0: cluster quality, acceptance, cluster energy
                    time_pass,    // 1: timing
                    totalE_pass,  // 2: total energy
                    Pt_pass,      // 3: Pt
                    dphi_pass     // 4: azimuthal angle
                };
                for (int s = 0; s < 5 && step_cut[s]; ++s) {
                    out.h_gamma_totalE[s]->Fill(totalE);
                    out.h2_gamma_hits[s]->Fill(x_g[0], y_g[0]);
                    out.h2_gamma_hits[s]->Fill(x_g[1], y_g[1]);
                    out.h2_gamma_hits[s]->Fill(x_e, y_e);
                    out.h_gamma_E[s]->Fill(E_g[0]);
                    out.h_gamma_E[s]->Fill(E_g[1]);
                    out.h_gamma_E[s]->Fill(E_e);
                    out.h_gamma_E_gamma[s]->Fill(E_g[0]);
                    out.h_gamma_E_gamma[s]->Fill(E_g[1]);
                    out.h_gamma_E_electron[s]->Fill(E_e);
                    out.h2_gamma_E_gamma_vs_E_electron[s]->Fill(E_g[0], E_e);
                    out.h2_gamma_E_gamma_vs_E_electron[s]->Fill(E_g[1], E_e);
                    out.h2_gamma_E_gamma_vs_E_gamma[s]->Fill(E_g[0], E_g[1]);
                    out.h2_gamma_E_angle_gamma[s]->Fill(theta_g[0], E_g[0]);
                    out.h2_gamma_E_angle_gamma[s]->Fill(theta_g[1], E_g[1]);
                    out.h2_gamma_E_angle_electron[s]->Fill(theta_e, E_e);
                    out.h_gamma_ptx[s]->Fill(ptx);
                    out.h_gamma_pty[s]->Fill(pty);
                    out.h2_gamma_Pt[s]->Fill(ptx, pty);
                    out.h_gamma_tDiff[s]->Fill(dt[0]);
                    out.h_gamma_tDiff[s]->Fill(dt[1]);
                    out.h_gamma_dphi[s]->Fill(dphi);
                    out.h_gamma_mass[s]->Fill(mass);
                }
            }

            //loop over all clusters for GEM matching

            if (ev.matchNum == 0) continue;

            struct Hits{
                float xu, yu, zu; // upstream
                float xd, yd, zd; // downstream
                float x, y, z; // projected to HyCal plane
                float E, t; // cluster energy and time
            };

            std::vector<Hits> hits_candidate;
            out.h_3cl_cluster_num->Fill(ev.matchNum);
            const int n_match_3cl = std::clamp(ev.matchNum, 0, prad2::kMaxClusters);
            for (int j = 0; j < n_match_3cl; ++j) {
                int idx = ev.mHit_cl_index[j];
                if (idx < 0 || idx >= ev.n_clusters || idx >= prad2::kMaxClusters) continue;
                if(ev.cl_nblocks[idx] < 2) continue;
                if(fdec::test_bit(ev.cl_flag[idx], fdec::kInnerBound)) continue;
                if(fdec::test_bit(ev.cl_flag[idx], fdec::kOuterBound)) continue;
                if(ev.cl_energy[idx] < 70.f || ev.cl_energy[idx] > 0.75 * Ebeam) continue;

                Hits hit{};
                hit.xu = ev.mHit_gx[j][1];
                hit.yu = ev.mHit_gy[j][1];
                hit.zu = ev.mHit_gz[j][1];
                hit.xd = ev.mHit_gx[j][0];
                hit.yd = ev.mHit_gy[j][0];
                hit.zd = ev.mHit_gz[j][0];
                float scale = ev.mHit_z[j] / hit.zu;
                hit.x = hit.xu * scale;
                hit.y = hit.yu * scale;
                hit.z = hit.zu * scale;
                hit.E = ev.mHit_E[j];
                hit.t = ev.cl_time[idx];
                hits_candidate.push_back(hit);
            }

            std::vector<Hits> hits;
            std::sort(hits_candidate.begin(), hits_candidate.end(),
                      [](const Hits &a, const Hits &b) { return a.E > b.E; });
            if (hits_candidate.empty()) continue;
            out.h_3cl_cluster_num_cut_cl->Fill(hits_candidate.size());
            hits.push_back(hits_candidate[0]);
            // loop over the candidates to check the timing correlation,
            // should +-2ns around the highest energy matched cluster which is the 1st in the candidates vector
            for (int j = 1; j < hits_candidate.size(); j++){
                out.h_3cl_tDiff_raw->Fill(hits_candidate[j].t - hits_candidate[0].t);
                if (fabs(hits_candidate[j].t - hits_candidate[0].t) < 2.0f)
                    hits.push_back(hits_candidate[j]);
            }
            out.h_3cl_cluster_num_cut_cl_t->Fill(hits.size());

            if(hits.size() == 3){
                float dt[2] = {hits[1].t - hits[0].t, hits[2].t - hits[0].t};
                float totalE = hits[0].E + hits[1].E + hits[2].E;
                float Sigma[3] = {0.035f * std::sqrt(hits[0].E * 1000.f),
                                  0.035f * std::sqrt(hits[1].E * 1000.f),
                                  0.035f * std::sqrt(hits[2].E * 1000.f)};
                float totalSigma = std::sqrt(Sigma[0]*Sigma[0] + Sigma[1]*Sigma[1] + Sigma[2]*Sigma[2]);
                float theta[3] = {std::atan2(std::sqrt(hits[0].x * hits[0].x + hits[0].y * hits[0].y), hits[0].z) * 180.f / static_cast<float>(M_PI),
                                  std::atan2(std::sqrt(hits[1].x * hits[1].x + hits[1].y * hits[1].y), hits[1].z) * 180.f / static_cast<float>(M_PI),
                                  std::atan2(std::sqrt(hits[2].x * hits[2].x + hits[2].y * hits[2].y), hits[2].z) * 180.f / static_cast<float>(M_PI)};

                // 4-momentum calculation for each single hit and each pair of hits
                TLorentzVector p[3], p12, p02, p01;
                for (int k = 0; k < 3; ++k)
                    PhysicsTools::HitP4(hits[k].x, hits[k].y, hits[k].z, hits[k].E,
                                        PhysicsTools::kElectronMass, p[k]);

                p12 = p[1] + p[2];
                p02 = p[0] + p[2];
                p01 = p[0] + p[1];

                // get the azimuthal angles for each single hit and each pair of hits
                float phi[3] = {
                    std::atan2(hits[0].y, hits[0].x) * 180.f / static_cast<float>(M_PI),
                    std::atan2(hits[1].y, hits[1].x) * 180.f / static_cast<float>(M_PI),
                    std::atan2(hits[2].y, hits[2].x) * 180.f / static_cast<float>(M_PI)
                };
                float phi_pair[3] = {
                    static_cast<float>(std::atan2(p12.Y(), p12.X()) * 180.f / M_PI),
                    static_cast<float>(std::atan2(p02.Y(), p02.X()) * 180.f / M_PI),
                    static_cast<float>(std::atan2(p01.Y(), p01.X()) * 180.f / M_PI)
                };

                // Phi difference for each combination (the pair of hits and the remaining single hit)
                float dphi[3] = {
                    std::fabs(phi_pair[0] - phi[0]),
                    std::fabs(phi_pair[1] - phi[1]),
                    std::fabs(phi_pair[2] - phi[2])
                };
                bool dphi_pass[3] = {
                    fabs(dphi[0] - 180.0f) < 10.0,
                    fabs(dphi[1] - 180.0f) < 10.0,
                    fabs(dphi[2] - 180.0f) < 10.0
                };

                // Pt x and Pt y calculation using TLorentzVector
                float ptx = p[0].Px() + p[1].Px() + p[2].Px();
                float pty = p[0].Py() + p[1].Py() + p[2].Py();

                // Invariant mass calculation using TLorentzVector for each pair of hits
                float mass[3];
                mass[0] = p12.M();
                mass[1] = p02.M();
                mass[2] = p01.M();

                // judge if it's Moller event plus an accidental smaller energy cluster
                // For the 2 combinations of biggest energy hit with the ohter 2 hits,
                // check (total energy, enrgy angle correlation, coplanarity) to identify Moller events
                bool Moller_event = false;
                for (int j = 1; j < 3; ++j) {
                    if (fabs(hits[0].E + hits[j].E - Ebeam) < 3. * sqrt(Sigma[0]*Sigma[0] + Sigma[j]*Sigma[j])) {
                        if (fabs(fabs(phi[0] - phi[j]) - 180.0) < 10.0) {
                            float expectE1 = physics.ExpectedEnergy(theta[0], Ebeam, "ee");
                            float expectE2 = physics.ExpectedEnergy(theta[j], Ebeam, "ee");
                            if (fabs(hits[0].E - expectE1) < 3. * 0.035 * sqrt(expectE1 * 1000.0) && 
                                fabs(hits[j].E - expectE2) < 3. * 0.035 * sqrt(expectE2 * 1000.0)) {
                                Moller_event = true;
                            }
                        }
                    }
                }

                // Vertex Z calculation for 3-cluster events
                // Use radial projection to estimate vertex Z position
                float vertexZ[3];
                for (int j = 0; j < 3; ++j) {
                    float r1 = std::sqrt(hits[j].xu * hits[j].xu + hits[j].yu * hits[j].yu);
                    float z1 = hits[j].zu;
                    float r2 = std::sqrt(hits[j].xd * hits[j].xd + hits[j].yd * hits[j].yd);
                    float z2 = hits[j].zd;

                    // project to r = 0 to estimate vertex Z position
                    if (r2 <= r1) vertexZ[j] = -9999.;
                    else vertexZ[j] = z1 - r1 * (z2 - z1) / (r2 - r1);
                }

                bool totalE_pass = std::fabs(totalE - Ebeam) < 4. * totalSigma;
                bool Pt_pass = std::sqrt(ptx * ptx + pty * pty) < 5.0f;
                bool pos_pass = inHyCal(hits[0].x, hits[0].y) && inHyCal(hits[1].x, hits[1].y) && inHyCal(hits[2].x, hits[2].y);
                bool anzimuthal_pass = dphi_pass[0] || dphi_pass[1] || dphi_pass[2];
                bool vertexZ_pass = std::fabs(vertexZ[0]) < 2000.0 && std::fabs(vertexZ[1]) < 2000.0 && std::fabs(vertexZ[2]) < 2000.0;

                // cut step s applies the cuts of steps 0..s
                const bool step_cut[6] = {
                    pos_pass,         // 0: acceptance position
                    totalE_pass,      // 1: total energy
                    Pt_pass,          // 2: Pt
                    anzimuthal_pass,  // 3: azimuthal angle
                    !Moller_event,    // 4: Moller events
                    vertexZ_pass      // 5: vertex Z
                };
                for (int s = 0; s < 6 && step_cut[s]; ++s) {
                    // from step 3 on, the mass and dphi plots keep only back-to-back combinations
                    const bool all_comb = s < 3;
                    out.h_3cl_totalE[s]->Fill(totalE);
                    out.h_3cl_tDiff[s]->Fill(dt[0]);
                    out.h_3cl_tDiff[s]->Fill(dt[1]);
                    out.h_3cl_ptx[s]->Fill(ptx);
                    out.h_3cl_pty[s]->Fill(pty);
                    out.h2_3cl_Pt[s]->Fill(ptx, pty);
                    for (int j = 0; j < 3; ++j) {
                        out.h_3cl_E[s]->Fill(hits[j].E);
                        out.h2_3cl_hits[s]->Fill(hits[j].x, hits[j].y);
                        out.h2_3cl_E_angle[s]->Fill(theta[j], hits[j].E);
                        out.h_3cl_yield[s]->Fill(theta[j]);
                        out.h_3cl_vertexZ[s]->Fill(vertexZ[j]);
                        if (all_comb || dphi_pass[j]) {
                            out.h_3cl_mass[s]->Fill(mass[j]);
                            out.h_3cl_dphi[s]->Fill(dphi[j]);
                        }
                    }
                    if (all_comb || dphi_pass[2]) out.h_3cl_mass_1comb[s]->Fill(mass[2]);
                    if (all_comb || dphi_pass[0]) out.h_3cl_mass_2comb[s]->Fill(mass[0]);
                    if (all_comb || dphi_pass[1]) out.h_3cl_mass_2comb[s]->Fill(mass[1]);
                }
            }
        }
    }
    out.processed += n;
    return true;
}

static void mergeResult(QuickResult &dst, const QuickResult &src, fdec::HyCalSystem &hycal)
{
    AddAll(dst.all, src.all);
    dst.physics->GetEnergyVsModuleHist()->Add(src.physics->GetEnergyVsModuleHist());
    dst.physics->GetEnergyVsThetaHist()->Add(src.physics->GetEnergyVsThetaHist());
    dst.physics->GetMollerPhiDiffHist()->Add(src.physics->GetMollerPhiDiffHist());
    for (int i = 0; i < hycal.module_count(); ++i) {
        int module_id = hycal.module(i).id;
        TH1F *d = dst.physics->GetModuleEnergyHist(module_id);
        TH1F *s = src.physics->GetModuleEnergyHist(module_id);
        if (d && s) d->Add(s);
    }
    dst.processed += src.processed;
}

// ── Main ─────────────────────────────────────────────────────────────────

int main(int argc, char *argv[])
{
    std::string output;
    float Ebeam = 2108.f;
    int max_events = -1;
    int num_threads = 4;
    int opt;
    while ((opt = getopt(argc, argv, "o:n:j:")) != -1) {
        switch (opt) {
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'j': num_threads = std::atoi(optarg); break;
        }
    }
    // collect input files (can be files, directories, or mixed)
    const std::vector<std::string> root_files = CollectInputs(argc, argv, optind, IsReconRootName);
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-j threads]\n";
        return 1;
    }
    num_threads = std::max(1, std::min(num_threads, static_cast<int>(root_files.size())));
    analysis::InitRootThreading();
    TH1::AddDirectory(kFALSE);

    std::string dbDir = prad2::database_dir();

    // --- load run config: assign run_id and Ebeam from gRunConfig ---
    const int run_id = analysis::get_run_int(root_files[0]);
    gRunConfig = analysis::LoadRunConfig(dbDir + "/runinfo/general.json", run_id);
    Ebeam = gRunConfig.Ebeam > 0.f ? gRunConfig.Ebeam : Ebeam;

    std::cout << "Processing run " << run_id << " with Ebeam = " << Ebeam << " MeV\n";

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    std::cout << "Processing " << root_files.size() << " file(s) with "
              << num_threads << " thread(s)\n";

    const auto file_limits = DistributeEventBudget(root_files, "recon", max_events > 0 ? max_events : -1);

    auto merged = makeResult(hycal);
    std::atomic<size_t> next_file{0};
    std::atomic<int> errors{0};
    std::mutex io_mtx;
    std::mutex merge_mtx;
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int t = 0; t < num_threads; ++t) {
        threads.emplace_back([&]() {
            auto res = makeResult(hycal);
            while (true) {
                size_t idx = next_file.fetch_add(1);
                if (idx >= root_files.size()) break;
                {
                    std::lock_guard<std::mutex> lk(io_mtx);
                    std::cerr << "Processing file [" << (idx + 1) << "/"
                              << root_files.size() << "]: " << root_files[idx] << "\n";
                }
                if (!processFile(root_files[idx], file_limits[idx], Ebeam, *res)) {
                    ++errors;
                    continue;
                }
            }
            {
                std::lock_guard<std::mutex> lk(merge_mtx);
                mergeResult(*merged, *res, hycal);
            }
        });
    }
    for (auto &t : threads) t.join();
    if (errors > 0) return 1;
    for (int i = 1; i <= merged->h_ep_ee_ratio->GetNbinsX(); i++) {
        double ep = merged->h_ep_yield->GetBinContent(i);
        double ee = merged->h_ee_yield->GetBinContent(i);
        merged->h_ep_ee_ratio->SetBinContent(i, ee > 0. ? ep / ee : 0.);
    }
    PhysicsTools &physics = *merged->physics;
    TH2F *hit_pos = merged->hit_pos.get();
    TH1F *h_1cl = merged->h_1cl.get();
    TH1F *h_2cl = merged->h_2cl.get();
    TH1F *h_all = merged->h_all.get();
    TH1F *h_tot = merged->h_tot.get();
    TH2F *h2_energy_theta_ep_ee = merged->h2_energy_theta_ep_ee.get();

    TString outName = output;
    if (outName.IsNull())
        outName = makeDefaultOutput(root_files[0]);
    TFile outfile(outName, "RECREATE");

    // --- write output ---
    outfile.cd();
    hit_pos->Write();

    merged->h2_ep_hits->Write();
    merged->h2_ee_hits->Write();
    merged->h_ee_center_x->Write();
    merged->h_ee_center_y->Write();
    merged->h_ee_vertex_z->Write();
    merged->h2_ep_hits_hc->Write();
    merged->h2_ee_hits_hc->Write();
    merged->h_ee_center_x_hc->Write();
    merged->h_ee_center_y_hc->Write();
    merged->h_ee_vertex_z_hc->Write();

    outfile.mkdir("energy_plots"); outfile.cd("energy_plots");
    if (physics.GetEnergyVsModuleHist()) physics.GetEnergyVsModuleHist()->Write();
    if (physics.GetEnergyVsThetaHist())  physics.GetEnergyVsThetaHist()->Write();
    h_1cl->Write(); h_2cl->Write(); h_all->Write(); h_tot->Write();
    h2_energy_theta_ep_ee->Write();
    merged->h2_ep_E_angle->Write();
    merged->h2_ee_E_angle->Write();
    merged->h2_ep_E_angle_hc->Write();
    merged->h2_ee_E_angle_hc->Write();

    outfile.cd();
    outfile.mkdir("physics_yields"); outfile.cd("physics_yields");
    merged->h_ep_yield->Write();
    merged->h_ee_yield->Write();
    merged->h_ep_ee_ratio->Write();

    outfile.cd();
    outfile.mkdir("moller_analysis"); outfile.cd("moller_analysis");
    if (physics.GetMollerPhiDiffHist()) physics.GetMollerPhiDiffHist()->Write();
    if (physics.GetMollerXHist()) physics.GetMollerXHist()->Write();
    if (physics.GetMollerYHist()) physics.GetMollerYHist()->Write();
    if (physics.GetMollerZHist()) physics.GetMollerZHist()->Write();
    merged->h_ee_tDiff->Write();
    merged->h_ee_invariant_mass->Write();

    for (int i = 0; i < 5; i++){
        outfile.cd();
        outfile.mkdir(Form("x17_gamma_cut_steps_%d", i)); outfile.cd(Form("x17_gamma_cut_steps_%d", i));
        merged->h_gamma_totalE[i]->Write();
        merged->h2_gamma_hits[i]->Write();
        merged->h_gamma_E[i]->Write();
        merged->h_gamma_E_gamma[i]->Write();
        merged->h_gamma_E_electron[i]->Write();
        merged->h2_gamma_E_gamma_vs_E_electron[i]->Write();
        merged->h2_gamma_E_gamma_vs_E_gamma[i]->Write();
        merged->h2_gamma_E_angle_gamma[i]->Write();
        merged->h2_gamma_E_angle_electron[i]->Write();
        merged->h_gamma_ptx[i]->Write();
        merged->h_gamma_pty[i]->Write();
        merged->h2_gamma_Pt[i]->Write();
        merged->h_gamma_tDiff[i]->Write();
        merged->h_gamma_dphi[i]->Write();
        merged->h_gamma_mass[i]->Write();
    }

    outfile.cd();
    outfile.mkdir("x17_gem"); outfile.cd("x17_gem");
    merged->h_3cl_cluster_num->Write();
    merged->h_3cl_cluster_num_cut_cl->Write();
    merged->h_3cl_tDiff_raw->Write();
    merged->h_3cl_cluster_num_cut_cl_t->Write();
    for (int i = 0; i < 6; ++i) {
        outfile.cd();
        outfile.mkdir(Form("x17_gem_cut_steps_%d", i)); outfile.cd(Form("x17_gem_cut_steps_%d", i));
        merged->h_3cl_totalE[i]->Write();
        merged->h2_3cl_hits[i]->Write();
        merged->h2_3cl_E_angle[i]->Write();
        merged->h_3cl_E[i]->Write();
        merged->h_3cl_yield[i]->Write();
        merged->h_3cl_mass[i]->Write();
        merged->h_3cl_ptx[i]->Write();
        merged->h_3cl_pty[i]->Write();
        merged->h2_3cl_Pt[i]->Write();
        merged->h_3cl_tDiff[i]->Write();
        merged->h_3cl_dphi[i]->Write();
        merged->h_3cl_vertexZ[i]->Write();
        merged->h_3cl_mass_1comb[i]->Write();
        merged->h_3cl_mass_2comb[i]->Write();
    }

    outfile.mkdir("module_energy"); outfile.cd("module_energy");
    for (int i = 0; i < hycal.module_count(); i++) {
        int module_id = hycal.module(i).id;
        TH1F *h = physics.GetModuleEnergyHist(module_id);
        if (h && h->GetEntries() > 0) h->Write();
    }

    outfile.Close();

    std::cerr << "Result saved -> " << outName.Data() << "\n";
}

// ── Helpers ──────────────────────────────────────────────────────────────
static std::string makeDefaultOutput(const std::string &input_path)
{
    fs::path p(input_path);
    std::string name = p.filename().string();
    const std::string ext = ".root";
    if (name.size() >= ext.size() &&
        name.compare(name.size() - ext.size(), ext.size(), ext) == 0) {
        name.insert(name.size() - ext.size(), "_quick_check");
    } else {
        name += "_quick_check.root";
    }
    return (p.parent_path() / name).string();
}
