// quick_check.C — ROOT script version 
//
// Reads reconstructed ROOT tree (output of replay_recon), runs physics
// analysis using PhysicsTools and MatchingTools from prad2det, and saves
// histograms to an output ROOT file.
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
#include "MatchingTools.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "ConfigSetup.h"

#include <TFile.h>
#include <TTree.h>
#include <TH1F.h>
#include <TH2F.h>
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

constexpr float M_ELECTRON = 0.5109989461; // MeV

// Aliases for the shared replay data structures
using EventVars_Recon = prad2::ReconEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);
static std::string makeDefaultOutput(const std::string &input_path);

bool inHyCal(float xmm, float ymm) {
    const float module = 20.75; // mm
    return (fabs(xmm) > module * 2.2 || fabs(ymm) > module * 2.2)
        && (fabs(xmm) < module * 16. && fabs(ymm) < module * 16.);
}

static float electronPairInvariantMass(
    float x1, float y1, float z1, float E1,
    float x2, float y2, float z2, float E2, bool is_gamma = false)
{
    float electron_mass = M_ELECTRON; // MeV
    if (is_gamma) electron_mass = 0.f;
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
    std::unique_ptr<TH1F> h_gamma_totalE;
    std::unique_ptr<TH2F> h2_gamma_hits;
    std::unique_ptr<TH2F> h2_gamma_E_angle;
    std::unique_ptr<TH1F> h_gamma_E;
    std::unique_ptr<TH1F> h_gamma_yield;
    std::unique_ptr<TH1F> h_gamma_mass;
    std::unique_ptr<TH1F> h_gamma_ptx;
    std::unique_ptr<TH1F> h_gamma_pty;
    std::unique_ptr<TH2F> h2_gamma_Pt;
    std::unique_ptr<TH1F> h_gamma_tDiff;

    std::unique_ptr<TH1F> h_gamma_totalE_pass;
    std::unique_ptr<TH2F> h2_gamma_hits_pass;
    std::unique_ptr<TH2F> h2_gamma_E_angle_pass;
    std::unique_ptr<TH1F> h_gamma_E_pass;
    std::unique_ptr<TH1F> h_gamma_yield_pass;
    std::unique_ptr<TH1F> h_gamma_mass_pass;
    std::unique_ptr<TH1F> h_gamma_ptx_pass;
    std::unique_ptr<TH1F> h_gamma_pty_pass;
    std::unique_ptr<TH2F> h2_gamma_Pt_pass;
    std::unique_ptr<TH1F> h_gamma_tDiff_pass;
    std::unique_ptr<TH1F> h_gamma_E_gamma_pass;
    std::unique_ptr<TH1F> h_gamma_E_electron_pass;
    std::unique_ptr<TH2F> h2_gamma_E_gamma_vs_E_electron_pass;
    std::unique_ptr<TH2F> h2_gamma_Etheta_gamma_pass;
    std::unique_ptr<TH2F> h2_gamma_Etheta_electron_pass;

    std::unique_ptr<TH1F> h_3cl_cluster_gem_num;
    std::unique_ptr<TH1F> h_3cl_totalE_gem;
    std::unique_ptr<TH2F> h2_3cl_hits_gem;
    std::unique_ptr<TH2F> h2_3cl_E_angle_gem;
    std::unique_ptr<TH1F> h_3cl_E_gem;
    std::unique_ptr<TH1F> h_3cl_yield_gem;
    std::unique_ptr<TH1F> h_3cl_mass_gem;
    std::unique_ptr<TH1F> h_3cl_ptx_gem;
    std::unique_ptr<TH1F> h_3cl_pty_gem;
    std::unique_ptr<TH2F> h2_3cl_Pt_gem;
    std::unique_ptr<TH1F> h_3cl_tDiff_gem;
    std::unique_ptr<TH1F> h_3cl_dphi_gem;

    std::unique_ptr<TH1F> h_3cl_totalE_gem_cut;
    std::unique_ptr<TH2F> h2_3cl_hits_gem_cut;
    std::unique_ptr<TH2F> h2_3cl_E_angle_gem_cut;
    std::unique_ptr<TH1F> h_3cl_E_gem_cut;
    std::unique_ptr<TH1F> h_3cl_yield_gem_cut;
    std::unique_ptr<TH1F> h_3cl_mass_gem_cut;
    std::unique_ptr<TH1F> h_3cl_ptx_gem_cut;
    std::unique_ptr<TH1F> h_3cl_pty_gem_cut;
    std::unique_ptr<TH2F> h2_3cl_Pt_gem_cut;
    std::unique_ptr<TH1F> h_3cl_tDiff_gem_cut;
    std::unique_ptr<TH1F> h_3cl_dphi_gem_cut;
};

static void detach(TH1 *h)
{
    if (h) h->SetDirectory(nullptr);
}

static std::unique_ptr<QuickResult> makeResult(fdec::HyCalSystem &hycal)
{
    auto r = std::make_unique<QuickResult>();
    r->physics = std::make_unique<PhysicsTools>(hycal);
    r->hit_pos = std::make_unique<TH2F>("hit_pos",
        "Hit positions;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h_1cl = std::make_unique<TH1F>("one_cluster_energy",
        "Single-cluster energy;E (MeV);Counts", 4000, 0, 4000);
    r->h_2cl = std::make_unique<TH1F>("two_cluster_energy",
        "Two-cluster energy;E (MeV);Counts", 4000, 0, 4000);
    r->h_all = std::make_unique<TH1F>("clusters_energy",
        "All clusters;E (MeV);Counts", 4000, 0, 4000);
    r->h_tot = std::make_unique<TH1F>("total_energy",
        "Total energy per event;E (MeV);Counts", 4000, 0, 4000);
    r->h2_energy_theta_ep_ee = std::make_unique<TH2F>("energy_vs_theta",
        "Energy vs Theta(1 cluster);Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);

    r->h2_ep_hits = std::make_unique<TH2F>("ep_hits",
        "EP Hit positions;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ee_hits = std::make_unique<TH2F>("ee_hits",
        "EE Hit positions;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ep_E_angle = std::make_unique<TH2F>("ep_E_angle",
        "EP Energy vs Angle;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);
    r->h2_ee_E_angle = std::make_unique<TH2F>("ee_E_angle",
        "EE Energy vs Angle;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);

    r->h_ep_yield = std::make_unique<TH1F>("ep_yield",
        "EP Yield;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_ee_yield = std::make_unique<TH1F>("ee_yield",
        "EE Yield;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_ep_ee_ratio = std::make_unique<TH1F>("ep_ee_ratio",
        "EP/EE Yield Ratio;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_ee_tDiff = std::make_unique<TH1F>("ee_tDiff",
        "EE Time Difference;Time Difference (ns);Counts", 400, -10, 10);

    r->h_ee_center_x = std::make_unique<TH1F>("ee_center_x",
        "EE Center X;X (mm);Counts", 800, -20, 20);
    r->h_ee_center_y = std::make_unique<TH1F>("ee_center_y",
        "EE Center Y;Y (mm);Counts", 800, -20, 20);
    r->h_ee_vertex_z = std::make_unique<TH1F>("ee_vertex_z",
        "EE Vertex Z;Z (mm);Counts", 8000, 5000, 9000);

    r->h2_ep_hits_hc = std::make_unique<TH2F>("ep_hits_hc",
        "EP Hit positions hycal;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ee_hits_hc = std::make_unique<TH2F>("ee_hits_hc",
        "EE Hit positions hycal;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_ep_E_angle_hc = std::make_unique<TH2F>("ep_E_angle_hc",
        "EP Energy vs Angle hycal;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);
    r->h2_ee_E_angle_hc = std::make_unique<TH2F>("ee_E_angle_hc",
        "EE Energy vs Angle hycal;Theta (deg);Energy (MeV)", 160, 0, 8, 7500, 0, 5000);

    r->h_ee_center_x_hc = std::make_unique<TH1F>("ee_center_x_hc",
        "EE Center X hycal;X (mm);Counts", 800, -20, 20);
    r->h_ee_center_y_hc = std::make_unique<TH1F>("ee_center_y_hc",
        "EE Center Y hycal;Y (mm);Counts", 800, -20, 20);
    r->h_ee_vertex_z_hc = std::make_unique<TH1F>("ee_vertex_z_hc",
        "EE Vertex Z hycal;Z (mm);Counts", 8000, 5000, 9000);

    r->h_ee_invariant_mass = std::make_unique<TH1F>("ee_invariant_mass",
        "EE Invariant Mass;Mass (MeV);Counts", 400, 0, 100);

    //For X17
    r->h2_gamma_hits = std::make_unique<TH2F>("gamma_hits",
        "Gamma Channel Hit positions hycal;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_gamma_E_angle = std::make_unique<TH2F>("gamma_E_angle",
        "Gamma Channel Energy vs Angle hycal;Theta (deg);Energy (MeV)", 80, 0, 4, 7500, 0, 5000);
    r->h_gamma_E = std::make_unique<TH1F>("gamma_E",
        "Gamma Channel Cluster Energy;Energy (MeV);Counts", 3750, 0, 2500);
    r->h_gamma_totalE = std::make_unique<TH1F>("gamma_totalE",
        "Gamma Channel Total Energy;Total Energy (MeV);Counts", 7500, 0, 5000);
    r->h_gamma_yield = std::make_unique<TH1F>("gamma_yield",
        "Gamma Channel Yield;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_gamma_mass = std::make_unique<TH1F>("gamma_mass",
        "Gamma Pair Inv. Mass;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_gamma_ptx = std::make_unique<TH1F>("gamma_ptx",
        "Gamma Channel Ptx;Ptx (MeV);Counts", 200, -50, 50);
    r->h_gamma_pty = std::make_unique<TH1F>("gamma_pty",
        "Gamma Channel Pty;Pty (MeV);Counts", 200, -50, 50);
    r->h2_gamma_Pt = std::make_unique<TH2F>("gamma_Pt",
        "Gamma Channel Pt hycal;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);
    r->h_gamma_tDiff = std::make_unique<TH1F>("gamma_tDiff",
        "Gamma Channel Time Difference;Time Difference (ns);Counts", 200, 0, 20);

    r->h2_gamma_hits_pass = std::make_unique<TH2F>("gamma_hits_pass",
        "Gamma Channel Hit positions hycal - Pass;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_gamma_E_angle_pass = std::make_unique<TH2F>("gamma_E_angle_pass",
        "Gamma Channel Energy vs Angle hycal - Pass;Theta (deg);Energy (MeV)", 80, 0, 4, 7500, 0, 5000);
    r->h_gamma_E_pass = std::make_unique<TH1F>("gamma_E_pass",
        "Gamma Channel Cluster Energy - Pass;Energy (MeV);Counts", 3750, 0, 2500);
    r->h_gamma_totalE_pass = std::make_unique<TH1F>("gamma_totalE_pass",
        "Gamma Channel Total Energy - Pass;Total Energy (MeV);Counts", 7500, 0, 5000);
    r->h_gamma_yield_pass = std::make_unique<TH1F>("gamma_yield_pass",
        "Gamma Channel Yield - Pass;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_gamma_mass_pass = std::make_unique<TH1F>("gamma_mass_pass",
        "Gamma Pair Inv. Mass - Pass;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_gamma_ptx_pass = std::make_unique<TH1F>("gamma_ptx_pass",
        "Gamma Channel Ptx - Pass;Ptx (MeV);Counts", 200, -50, 50);
    r->h_gamma_pty_pass = std::make_unique<TH1F>("gamma_pty_pass",
        "Gamma Channel Pty - Pass;Pty (MeV);Counts", 200, -50, 50);
    r->h2_gamma_Pt_pass = std::make_unique<TH2F>("gamma_Pt_pass",
        "Gamma Channel Pt hycal - Pass;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);
    r->h_gamma_tDiff_pass = std::make_unique<TH1F>("gamma_tDiff_pass",
        "Gamma Channel Time Difference - Pass;Time Difference (ns);Counts", 200, 0, 20);
    r->h_gamma_E_gamma_pass = std::make_unique<TH1F>("gamma_E_gamma_pass",
        "Gamma Cluster Energy - Pass;Energy (MeV);Counts", 3750, 0, 2500);
    r->h_gamma_E_electron_pass = std::make_unique<TH1F>("gamma_E_electron_pass",
        "Electron Cluster Energy - Pass;Energy (MeV);Counts", 3750, 0, 2500);
    r->h2_gamma_E_gamma_vs_E_electron_pass = std::make_unique<TH2F>("gamma_E_gamma_vs_E_electron_pass",
        "Gamma vs Electron Cluster Energy - Pass;Electron Energy (MeV);Gamma Energy (MeV)", 7500, 0, 5000, 7500, 0, 5000);
    r->h2_gamma_Etheta_gamma_pass = std::make_unique<TH2F>("gamma_Etheta_gamma_pass",
        "Gamma Cluster Energy vs Angle - Pass;Theta (deg);Energy (MeV)", 80, 0, 4, 7500, 0, 5000);
    r->h2_gamma_Etheta_electron_pass = std::make_unique<TH2F>("gamma_Etheta_electron_pass",
        "Electron Cluster Energy vs Angle - Pass;Theta (deg);Energy (MeV)", 80, 0, 4, 7500, 0, 5000);

    // use gem matching to cut the 3-cluster events
    r->h_3cl_cluster_gem_num = std::make_unique<TH1F>( "3cl_cluster_gem_num",
        "GEM-matched Cluster Number;Number of Clusters;Counts", 20, 0, 20);
    r->h2_3cl_hits_gem = std::make_unique<TH2F>("3cl_hits_gem",
        "3-Cluster Hit positions hycal with GEM matching;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_3cl_E_angle_gem = std::make_unique<TH2F>("3cl_E_angle_gem",
        "3-Cluster Energy vs Angle hycal with GEM matching;Theta (deg);Energy (MeV)", 80, 0, 4, 7500, 0, 5000);
    r->h_3cl_E_gem = std::make_unique<TH1F>("3cl_E_gem",
        "3-Cluster Energy with GEM matching;Energy (MeV);Counts", 7500, 0, 5000);
    r->h_3cl_totalE_gem = std::make_unique<TH1F>("3cl_totalE_gem",
        "3-Cluster Total Energy with GEM matching;Total Energy (MeV);Counts", 7500, 0, 5000);
    r->h_3cl_yield_gem = std::make_unique<TH1F>("3cl_yield_gem",
        "3-Cluster Yield with GEM matching;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_3cl_mass_gem = std::make_unique<TH1F>("3cl_mass_gem",
        "3-Cluster Inv. Mass with GEM matching;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_3cl_ptx_gem = std::make_unique<TH1F>("3cl_ptx_gem",
        "3-Cluster Ptx with GEM matching;Ptx (MeV);Counts", 200, -50, 50);
    r->h_3cl_pty_gem = std::make_unique<TH1F>("3cl_pty_gem",
        "3-Cluster Pty with GEM matching;Pty (MeV);Counts", 200, -50, 50);
    r->h_3cl_tDiff_gem = std::make_unique<TH1F>("3cl_tDiff_gem",
        "3-Cluster Time Difference with GEM matching;Time Difference (ns);Counts", 200, 0, 20);
    r->h_3cl_dphi_gem = std::make_unique<TH1F>("3cl_dphi_gem",
        "3-Cluster Phi Difference with GEM matching;#Delta#phi (deg);Counts", 360*3, 0, 360);
    r->h2_3cl_Pt_gem = std::make_unique<TH2F>("3cl_Pt_gem",
        "3-Cluster Pt hycal with GEM matching;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);

    r->h2_3cl_hits_gem_cut = std::make_unique<TH2F>("3cl_hits_gem_cut",
        "3-Cluster Hit positions with GEM matching - Cut;X (mm);Y (mm)", 720, -360, 360, 720, -360, 360);
    r->h2_3cl_E_angle_gem_cut = std::make_unique<TH2F>("3cl_E_angle_gem_cut",
        "3-Cluster Energy vs Angle with GEM matching - Cut;Theta (deg);Energy (MeV)", 80, 0, 4, 7500, 0, 5000);
    r->h_3cl_E_gem_cut = std::make_unique<TH1F>("3cl_E_gem_cut",
        "3-Cluster Energy with GEM matching - Cut;Energy (MeV);Counts", 7500, 0, 5000);
    r->h_3cl_totalE_gem_cut = std::make_unique<TH1F>("3cl_totalE_gem_cut",
        "3-Cluster Total Energy with GEM matching - Cut;Total Energy (MeV);Counts", 7500, 0, 5000);
    r->h_3cl_yield_gem_cut = std::make_unique<TH1F>("3cl_yield_gem_cut",
        "3-Cluster Yield with GEM matching - Cut;Scattering Angle (deg);Counts", Nbins, binEdge);
    r->h_3cl_mass_gem_cut = std::make_unique<TH1F>("3cl_mass_gem_cut",
        "3-Cluster Inv. Mass with GEM matching - Cut;Inv. Mass (MeV);Counts", 1000, 0, 100);
    r->h_3cl_ptx_gem_cut = std::make_unique<TH1F>("3cl_ptx_gem_cut",
        "3-Cluster Ptx with GEM matching - Cut;Ptx (MeV);Counts", 200, -50, 50);
    r->h_3cl_pty_gem_cut = std::make_unique<TH1F>("3cl_pty_gem_cut",
        "3-Cluster Pty with GEM matching - Cut;Pty (MeV);Counts", 200, -50, 50);
    r->h2_3cl_Pt_gem_cut = std::make_unique<TH2F>("3cl_Pt_gem_cut",
        "3-Cluster Pt hycal with GEM matching - Cut;Ptx (MeV);Pty (MeV)", 400, -50, 50, 400, -50, 50);
    r->h_3cl_tDiff_gem_cut = std::make_unique<TH1F>("3cl_tDiff_gem_cut",
        "3-Cluster Time Difference with GEM matching - Cut;Time Difference (ns);Counts", 200, 0, 20);
    r->h_3cl_dphi_gem_cut = std::make_unique<TH1F>("3cl_dphi_gem_cut",
        "3-Cluster Phi Difference with GEM matching - Cut;#Delta#phi (deg);Counts", 360, 0, 360);

    detach(r->hit_pos.get());
    detach(r->h_1cl.get());
    detach(r->h_2cl.get());
    detach(r->h_all.get());
    detach(r->h_tot.get());
    detach(r->h2_energy_theta_ep_ee.get());
    detach(r->h2_ep_hits.get());
    detach(r->h2_ee_hits.get());
    detach(r->h2_ep_E_angle.get());
    detach(r->h2_ee_E_angle.get());
    detach(r->h_ep_yield.get());
    detach(r->h_ee_yield.get());
    detach(r->h_ep_ee_ratio.get());
    detach(r->h_ee_tDiff.get());
    detach(r->h_ee_center_x.get());
    detach(r->h_ee_center_y.get());
    detach(r->h_ee_vertex_z.get());
    detach(r->h2_ep_hits_hc.get());
    detach(r->h2_ee_hits_hc.get());
    detach(r->h2_ep_E_angle_hc.get());
    detach(r->h2_ee_E_angle_hc.get());
    detach(r->h_ee_center_x_hc.get());
    detach(r->h_ee_center_y_hc.get());
    detach(r->h_ee_vertex_z_hc.get());

    detach(r->h_ee_invariant_mass.get());

    detach(r->h2_gamma_hits.get());
    detach(r->h2_gamma_E_angle.get());
    detach(r->h_gamma_E.get());
    detach(r->h_gamma_totalE.get());
    detach(r->h_gamma_yield.get());
    detach(r->h_gamma_mass.get());
    detach(r->h_gamma_ptx.get());
    detach(r->h_gamma_pty.get());
    detach(r->h2_gamma_Pt.get());
    detach(r->h_gamma_tDiff.get());

    detach(r->h2_gamma_hits_pass.get());
    detach(r->h2_gamma_E_angle_pass.get());
    detach(r->h_gamma_E_pass.get());
    detach(r->h_gamma_totalE_pass.get());
    detach(r->h_gamma_yield_pass.get());
    detach(r->h_gamma_mass_pass.get());
    detach(r->h_gamma_ptx_pass.get());
    detach(r->h_gamma_pty_pass.get());
    detach(r->h2_gamma_Pt_pass.get());
    detach(r->h_gamma_tDiff_pass.get());
    detach(r->h_gamma_E_gamma_pass.get());
    detach(r->h_gamma_E_electron_pass.get());
    detach(r->h2_gamma_E_gamma_vs_E_electron_pass.get());
    detach(r->h2_gamma_Etheta_gamma_pass.get());
    detach(r->h2_gamma_Etheta_electron_pass.get());
    
    detach(r->h_3cl_cluster_gem_num.get());
    detach(r->h2_3cl_hits_gem.get());
    detach(r->h2_3cl_E_angle_gem.get());
    detach(r->h_3cl_E_gem.get());
    detach(r->h_3cl_totalE_gem.get());
    detach(r->h_3cl_yield_gem.get());
    detach(r->h_3cl_mass_gem.get());
    detach(r->h_3cl_ptx_gem.get());
    detach(r->h_3cl_pty_gem.get());
    detach(r->h2_3cl_Pt_gem.get());
    detach(r->h_3cl_tDiff_gem.get());
    detach(r->h_3cl_dphi_gem.get());
    detach(r->h2_3cl_hits_gem_cut.get());
    detach(r->h2_3cl_E_angle_gem_cut.get());
    detach(r->h_3cl_E_gem_cut.get());
    detach(r->h_3cl_totalE_gem_cut.get());
    detach(r->h_3cl_yield_gem_cut.get());
    detach(r->h_3cl_mass_gem_cut.get());
    detach(r->h_3cl_ptx_gem_cut.get());
    detach(r->h_3cl_pty_gem_cut.get());
    detach(r->h2_3cl_Pt_gem_cut.get());
    detach(r->h_3cl_tDiff_gem_cut.get());
    detach(r->h_3cl_dphi_gem_cut.get());
    return r;
}

static Long64_t reconEntries(const std::string &path)
{
    std::unique_ptr<TFile> f(TFile::Open(path.c_str(), "READ"));
    if (!f || f->IsZombie()) return 0;
    TTree *t = dynamic_cast<TTree *>(f->Get("recon"));
    return t ? t->GetEntries() : 0;
}

static bool processFile(const std::string &path,
                        Long64_t max_entries,
                        fdec::HyCalSystem &hycal,
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
                    if(physics.GetMollerPhiDiff(mp) < 10.f) {
                        out.h2_ee_hits_hc->Fill(ev.cl_x[0], ev.cl_y[0]);
                        out.h2_ee_hits_hc->Fill(ev.cl_x[1], ev.cl_y[1]);
                        float t1 = std::atan2(std::sqrt(ev.cl_x[0]*ev.cl_x[0] + ev.cl_y[0]*ev.cl_y[0]), ev.cl_z[0]) * 180.f / M_PI;
                        float t2 = std::atan2(std::sqrt(ev.cl_x[1]*ev.cl_x[1] + ev.cl_y[1]*ev.cl_y[1]), ev.cl_z[1]) * 180.f / M_PI;
                        out.h2_ee_E_angle_hc->Fill(t1, ev.cl_energy[0]);
                        out.h2_ee_E_angle_hc->Fill(t2, ev.cl_energy[1]);
                        out.mollers_hc.push_back(mp);
                        if (out.mollers_hc.size() > 3) out.mollers_hc.erase(out.mollers_hc.begin());

                        if (out.mollers_hc.size() > 1) {
                            auto center = physics.GetMollerCenter(out.mollers_hc[out.mollers_hc.size() - 2], mp);
                            out.h_ee_center_x_hc->Fill(center[0]);
                            out.h_ee_center_y_hc->Fill(center[1]);
                            if (out.mollers_hc.size() > 2) {
                                auto center2 = physics.GetMollerCenter(out.mollers_hc[out.mollers_hc.size() - 3], mp);
                                out.h_ee_center_x_hc->Fill(center2[0]);
                                out.h_ee_center_y_hc->Fill(center2[1]);
                            }
                        }
                        float vertex = physics.GetMollerZdistance(mp, Ebeam);
                        out.h_ee_vertex_z_hc->Fill(vertex);

                        // Calculate invariant mass
                        const float invariant_mass = electronPairInvariantMass(
                            ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0],
                            ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1]);
                        if (std::isfinite(invariant_mass))
                            out.h_ee_invariant_mass->Fill(invariant_mass);
                    }
                }
            }

            //loop over GEM matched hits find e-p events
            for(int j = 0; j < ev.matchNum; j++) {
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
                if (physics.isMoller_kinematic(theta[0], E[0], theta[1], E[1], Ebeam, 0.033f)
                    && fabs(physics.GetMollerPhiDiff(mev)) < 10.f)
                {
                    out.mollers.push_back(mev);
                    if (out.mollers.size() > 3) out.mollers.erase(out.mollers.begin());
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
                    if (out.mollers.size() > 1) {
                        auto center = physics.GetMollerCenter(out.mollers[out.mollers.size() - 2], mev);
                        out.h_ee_center_x->Fill(center[0]);
                        out.h_ee_center_y->Fill(center[1]);
                    }
                    if (out.mollers.size() > 2) {
                        auto center2 = physics.GetMollerCenter(out.mollers[out.mollers.size() - 3], mev);
                        out.h_ee_center_x->Fill(center2[0]);
                        out.h_ee_center_y->Fill(center2[1]);
                    }
                }
            }
        }

        // x17 trigger selection
        if (is_3cluster) {
            //try to find the gamma decay channel,
            //firstly try on the clean events(only 3 clusters on HyCal)
            if(ev.n_clusters == 3) {
                // code to analyze 3-cluster events for gamma decay channel goes here
                float x[3], y[3], z[3], E[3], t[3], theta[3];
                x[0] = ev.cl_x[0]; y[0] = ev.cl_y[0]; z[0] = ev.cl_z[0];
                x[1] = ev.cl_x[1]; y[1] = ev.cl_y[1]; z[1] = ev.cl_z[1];
                x[2] = ev.cl_x[2]; y[2] = ev.cl_y[2]; z[2] = ev.cl_z[2];
                E[0] = ev.cl_energy[0]; E[1] = ev.cl_energy[1]; E[2] = ev.cl_energy[2];

                t[0] = ev.cl_time[0]; t[1] = ev.cl_time[1]; t[2] = ev.cl_time[2];
                float tDiff = std::max({std::fabs(t[0] - t[1]), std::fabs(t[0] - t[2]), std::fabs(t[1] - t[2])});

                bool nblocks_ok = true;
                if(ev.cl_nblocks[0] <= 2 || ev.cl_nblocks[1] <= 2 || ev.cl_nblocks[2] <= 2)
                    nblocks_ok = false;

                theta[0] = std::atan2(std::sqrt(x[0]*x[0] + y[0]*y[0]), z[0]) * 180.f / M_PI;
                theta[1] = std::atan2(std::sqrt(x[1]*x[1] + y[1]*y[1]), z[1]) * 180.f / M_PI;
                theta[2] = std::atan2(std::sqrt(x[2]*x[2] + y[2]*y[2]), z[2]) * 180.f / M_PI;

                // check the GEM matching, the 2 gamma clusters should not have matching on both of 2 layers
                // the one electron cluster should have matching on both layers
                bool gem_ok = true;
                int gamma_count = 0, electron_count = 0;
                int electron_index = -1, gamma_index1 = -1, gamma_index2 = -1;
                for (int i = 0; i < 3; i++) {
                    bool has_downstream = ev.matchFlag[i] & 1u << 0 || ev.matchFlag[i] & 1u << 1;
                    bool has_upstream = ev.matchFlag[i] & 1u << 2 || ev.matchFlag[i] & 1u << 3;
                    if (has_downstream && has_upstream) {
                        electron_index = i;
                        electron_count++;
                    } else if (!has_downstream && !has_upstream) {
                        if (gamma_index1 == -1) gamma_index1 = i;
                        else gamma_index2 = i;
                        gamma_count++;
                    } else {
                        gem_ok = false;
                    }
                }
                if (gamma_count != 2 || electron_count != 1) gem_ok = false;

                if (gem_ok) {
                    // Pt x and Pt y calculation
                    auto get_pt = [](float x, float y, float z, float energy, bool is_gamma) {
                        float mass = M_ELECTRON;
                        if (is_gamma) mass = 0.f;
                        const float norm = std::sqrt(x*x + y*y + z*z);
                        if (norm <= 0.f || energy < mass)
                            return std::pair<float, float>{0.f, 0.f};
                        const float p = std::sqrt(std::max(
                            0.f, energy*energy - mass*mass));
                        return std::pair<float, float>{p*x/norm, p*y/norm};
                    };
                    const auto [pxe, pye] = get_pt(x[electron_index], y[electron_index], z[electron_index], E[electron_index], false);
                    const auto [pxgamma1, pygamma1] = get_pt(x[gamma_index1], y[gamma_index1], z[gamma_index1], E[gamma_index1], true);
                    const auto [pxgamma2, pygamma2] = get_pt(x[gamma_index2], y[gamma_index2], z[gamma_index2], E[gamma_index2], true);
                    const float ptx = pxe + pxgamma1 + pxgamma2;
                    const float pty = pye + pygamma1 + pygamma2;

                    float mass = electronPairInvariantMass(x[gamma_index1], y[gamma_index1], z[gamma_index1], E[gamma_index1],
                        x[gamma_index2], y[gamma_index2], z[gamma_index2], E[gamma_index2], true);

                    bool time_cut = false, totalE_cut = false, Pt_cut = false, clusterE_cut = false, inHyCal_cut = false;

                    if(tDiff < 2.f)
                        time_cut = true;
                    if(E[gamma_index1] + E[gamma_index2] + E[electron_index] < Ebeam + 250.f && E[gamma_index1] + E[gamma_index2] + E[electron_index] > 0.8 * Ebeam)
                        totalE_cut = true;
                    if(std::sqrt(ptx*ptx + pty*pty) < 5.f)
                        Pt_cut = true;
                    if(inHyCal(x[gamma_index1], y[gamma_index1]) && inHyCal(x[gamma_index2], y[gamma_index2]) && inHyCal(x[electron_index], y[electron_index]))
                        inHyCal_cut = true;
                    if(E[gamma_index1] > 70.f && E[gamma_index2] > 70.f && E[electron_index] > 70.f && E[gamma_index1] < 0.75 * Ebeam && E[gamma_index2] < 0.75 * Ebeam && E[electron_index] < 0.75 * Ebeam)
                        clusterE_cut = true;
                        
                    out.h_gamma_totalE->Fill(E[gamma_index1] + E[gamma_index2] + E[electron_index]);
                    out.h_gamma_ptx->Fill(ptx);
                    out.h_gamma_pty->Fill(pty);
                    out.h2_gamma_Pt->Fill(ptx, pty);
                    out.h_gamma_tDiff->Fill(tDiff);
                    out.h_gamma_E->Fill(E[gamma_index1]);
                    out.h_gamma_E->Fill(E[gamma_index2]);
                    out.h_gamma_E->Fill(E[electron_index]);
                    out.h2_gamma_hits->Fill(x[gamma_index1], y[gamma_index1]);
                    out.h2_gamma_hits->Fill(x[gamma_index2], y[gamma_index2]);
                    out.h2_gamma_hits->Fill(x[electron_index], y[electron_index]);
                    out.h2_gamma_E_angle->Fill(theta[gamma_index1], E[gamma_index1]);
                    out.h2_gamma_E_angle->Fill(theta[gamma_index2], E[gamma_index2]);
                    out.h2_gamma_E_angle->Fill(theta[electron_index], E[electron_index]);
                    out.h_gamma_yield->Fill(theta[gamma_index1]);
                    out.h_gamma_yield->Fill(theta[gamma_index2]);
                    out.h_gamma_yield->Fill(theta[electron_index]);
                    if (std::isfinite(mass)) out.h_gamma_mass->Fill(mass);
                    
                    if(nblocks_ok && totalE_cut && time_cut && Pt_cut && inHyCal_cut && clusterE_cut) {
                        out.h_gamma_totalE_pass->Fill(E[gamma_index1] + E[gamma_index2] + E[electron_index]);
                        out.h_gamma_ptx_pass->Fill(ptx);
                        out.h_gamma_pty_pass->Fill(pty);
                        out.h2_gamma_Pt_pass->Fill(ptx, pty);
                        out.h_gamma_tDiff_pass->Fill(tDiff);
                        out.h_gamma_E_pass->Fill(E[gamma_index1]);
                        out.h_gamma_E_pass->Fill(E[gamma_index2]);
                        out.h_gamma_E_pass->Fill(E[electron_index]);
                        out.h2_gamma_hits_pass->Fill(x[gamma_index1], y[gamma_index1]);
                        out.h2_gamma_hits_pass->Fill(x[gamma_index2], y[gamma_index2]);
                        out.h2_gamma_hits_pass->Fill(x[electron_index], y[electron_index]);
                        out.h2_gamma_E_angle_pass->Fill(theta[gamma_index1], E[gamma_index1]);
                        out.h2_gamma_E_angle_pass->Fill(theta[gamma_index2], E[gamma_index2]);
                        out.h2_gamma_E_angle_pass->Fill(theta[electron_index], E[electron_index]);
                        out.h_gamma_yield_pass->Fill(theta[gamma_index1]);
                        out.h_gamma_yield_pass->Fill(theta[gamma_index2]);
                        out.h_gamma_yield_pass->Fill(theta[electron_index]);
                        out.h_gamma_E_gamma_pass->Fill(E[gamma_index1]);
                        out.h_gamma_E_gamma_pass->Fill(E[gamma_index2]);
                        out.h_gamma_E_electron_pass->Fill(E[electron_index]);
                        out.h2_gamma_E_gamma_vs_E_electron_pass->Fill(E[electron_index], E[gamma_index1]);
                        out.h2_gamma_E_gamma_vs_E_electron_pass->Fill(E[electron_index], E[gamma_index2]);
                        out.h2_gamma_Etheta_gamma_pass->Fill(theta[gamma_index1], E[gamma_index1]);
                        out.h2_gamma_Etheta_gamma_pass->Fill(theta[gamma_index2], E[gamma_index2]);
                        out.h2_gamma_Etheta_electron_pass->Fill(theta[electron_index], E[electron_index]);
                        if (std::isfinite(mass)) out.h_gamma_mass_pass->Fill(mass);
                    }
                }
            }

            //loop over all clusters for GEM matching
            struct Hits{
                float xu, yu, zu; // upstream
                float xd, yd, zd; // downstream
                float x, y, z; // projected to HyCal plane
                float E, t; // cluster energy and time
            };

            std::vector<Hits> hits_candidate;
            out.h_3cl_cluster_gem_num->Fill(ev.matchNum);
            for (int j = 0; j < ev.matchNum; j++) {
                int idx = ev.mHit_cl_index[j];
                if(ev.cl_nblocks[idx] <= 2) continue;
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
            hits.push_back(hits_candidate[0]);
            // loop over the candidates to check the timing correlation,
            // should +-2ns around the highest energy matched cluster which is the 1st in the candidates vector
            for (int j = 1; j < hits_candidate.size(); j++){
                if (fabs(hits_candidate[j].t - hits_candidate[0].t) < 2.0f)
                    hits.push_back(hits_candidate[j]);
            }

            if(hits.size() == 3){
                float dt[2] = {hits[1].t - hits[0].t, hits[2].t - hits[0].t};
                float totalE = hits[0].E + hits[1].E + hits[2].E;
                float theta[3] = {std::atan2(std::sqrt(hits[0].x * hits[0].x + hits[0].y * hits[0].y), hits[0].z) * 180.f / M_PI,
                                  std::atan2(std::sqrt(hits[1].x * hits[1].x + hits[1].y * hits[1].y), hits[1].z) * 180.f / M_PI,
                                  std::atan2(std::sqrt(hits[2].x * hits[2].x + hits[2].y * hits[2].y), hits[2].z) * 180.f / M_PI};

                // 4-momentum calculation for each single hit and each pair of hits
                TLorentzVector p[3], p12, p02, p01;
                for (int k = 0; k < 3; ++k) {
                    const float norm = std::sqrt(
                        hits[k].x * hits[k].x + hits[k].y * hits[k].y + hits[k].z * hits[k].z);
                    const float p_mag = std::sqrt(hits[k].E * hits[k].E - M_ELECTRON * M_ELECTRON);
                    const float ux = hits[k].x / norm;
                    const float uy = hits[k].y / norm;
                    const float uz = hits[k].z / norm;
                    p[k] = TLorentzVector(p_mag * ux, p_mag * uy, p_mag * uz, hits[k].E);
                }

                p12 = p[1] + p[2];
                p02 = p[0] + p[2];
                p01 = p[0] + p[1];

                // get the azimuthal angles for each single hit and each pair of hits
                float phi[3] = {
                    std::atan2(hits[0].y, hits[0].x) * 180.f / M_PI,
                    std::atan2(hits[1].y, hits[1].x) * 180.f / M_PI,
                    std::atan2(hits[2].y, hits[2].x) * 180.f / M_PI
                };
                float phi_pair[3] = {
                    std::atan2(p12.Y(), p12.X()) * 180.f / M_PI,
                    std::atan2(p02.Y(), p02.X()) * 180.f / M_PI,
                    std::atan2(p01.Y(), p01.X()) * 180.f / M_PI
                };

                // Phi difference for each combination (the pair of hits and the remaining single hit)
                float dphi[3] = {
                    std::fabs(phi_pair[0] - phi[0]),
                    std::fabs(phi_pair[1] - phi[1]),
                    std::fabs(phi_pair[2] - phi[2])
                };

                // Pt x and Pt y calculation using TLorentzVector
                float ptx = p[0].Px() + p[1].Px() + p[2].Px();
                float pty = p[0].Py() + p[1].Py() + p[2].Py();

                // Invariant mass calculation using TLorentzVector for each pair of hits
                float mass[3];
                mass[0] = p12.M();
                mass[1] = p02.M();
                mass[2] = p01.M();

                bool totalE_pass = std::fabs(totalE - Ebeam) < 3. * 0.035 * sqrt(Ebeam * 1000.f);
                bool Pt_pass = std::sqrt(ptx * ptx + pty * pty) < 5.0f;
                bool pos_pass = inHyCal(hits[0].x, hits[0].y) && inHyCal(hits[1].x, hits[1].y) && inHyCal(hits[2].x, hits[2].y);

                out.h_3cl_totalE_gem->Fill(totalE);
                out.h_3cl_tDiff_gem->Fill(dt[0]);
                out.h_3cl_tDiff_gem->Fill(dt[1]);
                out.h_3cl_ptx_gem->Fill(ptx);
                out.h_3cl_pty_gem->Fill(pty);
                out.h2_3cl_Pt_gem->Fill(ptx, pty);
                for (int j = 0; j < 3; ++j) {
                    out.h_3cl_E_gem->Fill(hits[j].E);
                    out.h2_3cl_hits_gem->Fill(hits[j].x, hits[j].y);
                    out.h2_3cl_E_angle_gem->Fill(theta[j], hits[j].E);
                    out.h_3cl_yield_gem->Fill(theta[j]);
                    out.h_3cl_mass_gem->Fill(mass[j]);
                    out.h_3cl_dphi_gem->Fill(dphi[j]);
                }
                if (totalE_pass && Pt_pass && pos_pass) {
                    // Fill histograms for events passing all three cuts
                    out.h_3cl_totalE_gem_cut->Fill(totalE);
                    out.h_3cl_ptx_gem_cut->Fill(ptx);
                    out.h_3cl_pty_gem_cut->Fill(pty);
                    out.h2_3cl_Pt_gem_cut->Fill(ptx, pty);
                    out.h_3cl_tDiff_gem_cut->Fill(dt[0]);
                    out.h_3cl_tDiff_gem_cut->Fill(dt[1]);
                    for (int j = 0; j < 3; ++j) {
                        out.h_3cl_E_gem_cut->Fill(hits[j].E);
                        out.h2_3cl_hits_gem_cut->Fill(hits[j].x, hits[j].y);
                        out.h2_3cl_E_angle_gem_cut->Fill(theta[j], hits[j].E);
                        out.h_3cl_yield_gem_cut->Fill(theta[j]);
                        out.h_3cl_mass_gem_cut->Fill(mass[j]);
                        out.h_3cl_dphi_gem_cut->Fill(dphi[j]);
                    }
                }
            }
        }
    }
    out.processed += n;
    return true;
}

static void mergeResult(QuickResult &dst, const QuickResult &src, fdec::HyCalSystem &hycal)
{
    dst.hit_pos->Add(src.hit_pos.get());
    dst.h_1cl->Add(src.h_1cl.get());
    dst.h_2cl->Add(src.h_2cl.get());
    dst.h_all->Add(src.h_all.get());
    dst.h_tot->Add(src.h_tot.get());
    dst.h2_energy_theta_ep_ee->Add(src.h2_energy_theta_ep_ee.get());
    dst.h2_ep_hits->Add(src.h2_ep_hits.get());
    dst.h2_ee_hits->Add(src.h2_ee_hits.get());
    dst.h2_ep_E_angle->Add(src.h2_ep_E_angle.get());
    dst.h2_ee_E_angle->Add(src.h2_ee_E_angle.get());
    dst.h_ep_yield->Add(src.h_ep_yield.get());
    dst.h_ee_yield->Add(src.h_ee_yield.get());
    dst.h_ee_center_x->Add(src.h_ee_center_x.get());
    dst.h_ee_center_y->Add(src.h_ee_center_y.get());
    dst.h_ee_vertex_z->Add(src.h_ee_vertex_z.get());
    dst.h_ee_tDiff->Add(src.h_ee_tDiff.get());
    dst.h2_ep_hits_hc->Add(src.h2_ep_hits_hc.get());
    dst.h2_ee_hits_hc->Add(src.h2_ee_hits_hc.get());
    dst.h2_ep_E_angle_hc->Add(src.h2_ep_E_angle_hc.get());
    dst.h2_ee_E_angle_hc->Add(src.h2_ee_E_angle_hc.get());
    dst.h_ee_center_x_hc->Add(src.h_ee_center_x_hc.get());
    dst.h_ee_center_y_hc->Add(src.h_ee_center_y_hc.get());
    dst.h_ee_vertex_z_hc->Add(src.h_ee_vertex_z_hc.get());
    dst.h_ee_invariant_mass->Add(src.h_ee_invariant_mass.get());

    // X17 three-cluster histograms.
    dst.h_gamma_totalE->Add(src.h_gamma_totalE.get());
    dst.h2_gamma_hits->Add(src.h2_gamma_hits.get());
    dst.h2_gamma_E_angle->Add(src.h2_gamma_E_angle.get());
    dst.h_gamma_E->Add(src.h_gamma_E.get());
    dst.h_gamma_yield->Add(src.h_gamma_yield.get());
    dst.h_gamma_mass->Add(src.h_gamma_mass.get());
    dst.h_gamma_ptx->Add(src.h_gamma_ptx.get());
    dst.h_gamma_pty->Add(src.h_gamma_pty.get());
    dst.h2_gamma_Pt->Add(src.h2_gamma_Pt.get());
    dst.h_gamma_tDiff->Add(src.h_gamma_tDiff.get());

    dst.h_gamma_totalE_pass->Add(src.h_gamma_totalE_pass.get());
    dst.h2_gamma_hits_pass->Add(src.h2_gamma_hits_pass.get());
    dst.h2_gamma_E_angle_pass->Add(src.h2_gamma_E_angle_pass.get());
    dst.h_gamma_E_pass->Add(src.h_gamma_E_pass.get());
    dst.h_gamma_yield_pass->Add(src.h_gamma_yield_pass.get());
    dst.h_gamma_mass_pass->Add(src.h_gamma_mass_pass.get());
    dst.h_gamma_ptx_pass->Add(src.h_gamma_ptx_pass.get());
    dst.h_gamma_pty_pass->Add(src.h_gamma_pty_pass.get());
    dst.h2_gamma_Pt_pass->Add(src.h2_gamma_Pt_pass.get());
    dst.h_gamma_tDiff_pass->Add(src.h_gamma_tDiff_pass.get());
    dst.h_gamma_E_gamma_pass->Add(src.h_gamma_E_gamma_pass.get());
    dst.h_gamma_E_electron_pass->Add(src.h_gamma_E_electron_pass.get());
    dst.h2_gamma_E_gamma_vs_E_electron_pass->Add(src.h2_gamma_E_gamma_vs_E_electron_pass.get());
    dst.h2_gamma_Etheta_gamma_pass->Add(src.h2_gamma_Etheta_gamma_pass.get());
    dst.h2_gamma_Etheta_electron_pass->Add(src.h2_gamma_Etheta_electron_pass.get());

    // X17 three-cluster histograms with GEM matching.
    dst.h_3cl_cluster_gem_num->Add(src.h_3cl_cluster_gem_num.get());
    dst.h_3cl_totalE_gem->Add(src.h_3cl_totalE_gem.get());
    dst.h2_3cl_hits_gem->Add(src.h2_3cl_hits_gem.get());
    dst.h2_3cl_E_angle_gem->Add(src.h2_3cl_E_angle_gem.get());
    dst.h_3cl_E_gem->Add(src.h_3cl_E_gem.get());
    dst.h_3cl_yield_gem->Add(src.h_3cl_yield_gem.get());
    dst.h_3cl_mass_gem->Add(src.h_3cl_mass_gem.get());
    dst.h_3cl_ptx_gem->Add(src.h_3cl_ptx_gem.get());
    dst.h_3cl_pty_gem->Add(src.h_3cl_pty_gem.get());
    dst.h2_3cl_Pt_gem->Add(src.h2_3cl_Pt_gem.get());
    dst.h_3cl_tDiff_gem->Add(src.h_3cl_tDiff_gem.get());
    dst.h_3cl_dphi_gem->Add(src.h_3cl_dphi_gem.get());

    dst.h_3cl_totalE_gem_cut->Add(src.h_3cl_totalE_gem_cut.get());
    dst.h2_3cl_hits_gem_cut->Add(src.h2_3cl_hits_gem_cut.get());
    dst.h2_3cl_E_angle_gem_cut->Add(src.h2_3cl_E_angle_gem_cut.get());
    dst.h_3cl_E_gem_cut->Add(src.h_3cl_E_gem_cut.get());
    dst.h_3cl_yield_gem_cut->Add(src.h_3cl_yield_gem_cut.get());
    dst.h_3cl_mass_gem_cut->Add(src.h_3cl_mass_gem_cut.get());
    dst.h_3cl_ptx_gem_cut->Add(src.h_3cl_ptx_gem_cut.get());
    dst.h_3cl_pty_gem_cut->Add(src.h_3cl_pty_gem_cut.get());
    dst.h2_3cl_Pt_gem_cut->Add(src.h2_3cl_Pt_gem_cut.get());
    dst.h_3cl_tDiff_gem_cut->Add(src.h_3cl_tDiff_gem_cut.get());
    dst.h_3cl_dphi_gem_cut->Add(src.h_3cl_dphi_gem_cut.get());
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
    int run_id = 12345;
    
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
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; i++) {
        auto f = collectRootFiles(argv[i]);
        root_files.insert(root_files.end(), f.begin(), f.end());
    }
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: quick_check <input_recon.root|dir> [more files...] [-o out.root] [-n max_events] [-j threads]\n";
        return 1;
    }
    num_threads = std::max(1, std::min(num_threads, static_cast<int>(root_files.size())));
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);

    // --- database path ---
    std::string dbDir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);

    // --- load run config: assign run_id and Ebeam from gRunConfig ---
    run_id = analysis::get_run_int(root_files[0]);
    gRunConfig = analysis::LoadRunConfig(dbDir + "/runinfo/general.json", run_id);
    Ebeam = gRunConfig.Ebeam > 0.f ? gRunConfig.Ebeam : Ebeam;

    std::cout << "Processing run " << run_id << " with Ebeam = " << Ebeam << " MeV\n";

    // --- init detector system ---
    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    std::cout << "Processing " << root_files.size() << " file(s) with "
              << num_threads << " thread(s)\n";

    std::vector<Long64_t> file_limits(root_files.size(), -1);
    if (max_events > 0) {
        Long64_t remaining = max_events;
        for (size_t i = 0; i < root_files.size(); ++i) {
            Long64_t n = reconEntries(root_files[i]);
            file_limits[i] = std::min(n, remaining);
            remaining -= file_limits[i];
            if (remaining <= 0) {
                for (size_t j = i + 1; j < root_files.size(); ++j)
                    file_limits[j] = 0;
                break;
            }
        }
    }

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
                if (!processFile(root_files[idx], file_limits[idx], hycal, Ebeam, *res)) {
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

    outfile.cd();
    outfile.mkdir("x17_gamma"); outfile.cd("x17_gamma");
    merged->h_gamma_totalE->Write();
    merged->h2_gamma_hits->Write();
    merged->h2_gamma_E_angle->Write();
    merged->h_gamma_E->Write();
    merged->h_gamma_yield->Write();
    merged->h_gamma_mass->Write();
    merged->h_gamma_ptx->Write();
    merged->h_gamma_pty->Write();
    merged->h2_gamma_Pt->Write();
    merged->h_gamma_tDiff->Write();

    outfile.cd();
    outfile.mkdir("x17_gamma_pass"); outfile.cd("x17_gamma_pass");
    merged->h_gamma_totalE_pass->Write();
    merged->h2_gamma_hits_pass->Write();
    merged->h2_gamma_E_angle_pass->Write();
    merged->h_gamma_E_pass->Write();
    merged->h_gamma_yield_pass->Write();
    merged->h_gamma_mass_pass->Write();
    merged->h_gamma_ptx_pass->Write();
    merged->h_gamma_pty_pass->Write();
    merged->h2_gamma_Pt_pass->Write();
    merged->h_gamma_tDiff_pass->Write();
    merged->h_gamma_E_gamma_pass->Write();
    merged->h_gamma_E_electron_pass->Write();
    merged->h2_gamma_E_gamma_vs_E_electron_pass->Write();
    merged->h2_gamma_Etheta_gamma_pass->Write();
    merged->h2_gamma_Etheta_electron_pass->Write();

    outfile.cd();
    outfile.mkdir("x17_gem"); outfile.cd("x17_gem");
    merged->h_3cl_cluster_gem_num->Write();
    merged->h_3cl_totalE_gem->Write();
    merged->h2_3cl_hits_gem->Write();
    merged->h2_3cl_E_angle_gem->Write();
    merged->h_3cl_E_gem->Write();
    merged->h_3cl_yield_gem->Write();
    merged->h_3cl_mass_gem->Write();
    merged->h_3cl_ptx_gem->Write();
    merged->h_3cl_pty_gem->Write();
    merged->h2_3cl_Pt_gem->Write();
    merged->h_3cl_tDiff_gem->Write();
    merged->h_3cl_dphi_gem->Write();

    outfile.cd();
    outfile.mkdir("x17_gem_cut"); outfile.cd("x17_gem_cut");
    merged->h_3cl_totalE_gem_cut->Write();
    merged->h2_3cl_hits_gem_cut->Write();
    merged->h2_3cl_E_angle_gem_cut->Write();
    merged->h_3cl_E_gem_cut->Write();
    merged->h_3cl_yield_gem_cut->Write();
    merged->h_3cl_mass_gem_cut->Write();
    merged->h_3cl_ptx_gem_cut->Write();
    merged->h_3cl_pty_gem_cut->Write();
    merged->h2_3cl_Pt_gem_cut->Write();
    merged->h_3cl_tDiff_gem_cut->Write();
    merged->h_3cl_dphi_gem_cut->Write();

    outfile.cd("moller_analysis");

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
