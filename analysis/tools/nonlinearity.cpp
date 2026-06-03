
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
#include <TCanvas.h>
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

using EventVars_Recon = prad2::ReconEventData;

static std::vector<std::string> collectRootFiles(const std::string &path);
void process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal, 
    std::map<int, TH1F*> &energy_hists, PhysicsTools &physics, float Ebeam);

float resolution = 0.035; // pre-defined energy resolution

int main(int argc, char *argv[]){

    std::string output, input_3p5, input_0p7;
    
    int max_events = -1;
    int opt;
    while ((opt = getopt(argc, argv, "a:b:o:n:")) != -1) {
        switch (opt) {
            case 'a': input_3p5 = optarg; break;
            case 'b': input_0p7 = optarg; break;
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
        }
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

    // Energy histogram for each crystal
    std::map<int, TH1F*> energy_hists_3p5;
    std::map<int, TH1F*> energy_hists_0p7;
    for (int i = 0; i < hycal.module_count(); ++i) {
        const auto &m = hycal.module(i);
        if (!m.is_pwo4()) continue;
        std::string hname = "h_energy_" + m.name;
        std::string htitle = "Energy " + m.name + ";E (MeV);Counts";
        energy_hists_3p5[m.id] = new TH1F((hname + "_3p5").c_str(), (htitle + " (3.5)").c_str(), 400, 0, 4000);
        energy_hists_0p7[m.id] = new TH1F((hname + "_0p7").c_str(), (htitle + " (0.7)").c_str(), 400, 0, 4000);
    }

    // --- setup TChain and branches ---
    TChain *chain = new TChain("recon");
        chain->Add(input_3p5.c_str());
        std::cerr << "Added file: " << input_3p5 << "\n";
    TTree *tree = chain;
    if (!tree) {
        std::cerr << "Cannot find TTree 'recon' in input files\n";
        return 1;
    }

    EventVars_Recon ev;
    prad2::SetReconReadBranches(tree, ev);

    process_event(tree, ev, hycal, energy_hists_3p5, physics, 3485.f);

    // --- repeat for 0.7 GeV data ---
    TChain *chain2 = new TChain("recon");
        chain2->Add(input_0p7.c_str());
        std::cerr << "Added file: " << input_0p7 << "\n";
    TTree *tree2 = chain2;
    if (!tree2) {
        std::cerr << "Cannot find TTree 'recon' in input files\n";
        return 1;
    }

    EventVars_Recon ev2;
    prad2::SetReconReadBranches(tree2, ev2);
    process_event(tree2, ev2, hycal, energy_hists_0p7, physics, 729.f);

    // calculate non-linearity module by module and save to output file
    TFile outFile(output.empty() ? "nonlinearity_results.root" : output.c_str(), "RECREATE");
    for (int i = 0; i < hycal.module_count(); i++) {
        const auto &mod = hycal.module(i);
        int mod_id = mod.id;
        if (!mod.is_pwo4()) continue; // only look at PbWO4 crystals

        auto it_3p5 = energy_hists_3p5.find(mod_id);
        if (it_3p5 == energy_hists_3p5.end() || !it_3p5->second) continue;
        auto hist_3p5 = it_3p5->second;
        auto it_0p7 = energy_hists_0p7.find(mod_id);
        if (it_0p7 == energy_hists_0p7.end() || !it_0p7->second) continue;
        auto hist_0p7 = it_0p7->second;

        float x = mod.x, y = mod.y, z = 6270.f;
        float theta = std::atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;
        float e_p_exp_3p5 = physics.ExpectedEnergy(theta, 3485.f, "ep");
        float e_e_exp_3p5 = physics.ExpectedEnergy(theta, 3485.f, "ee");
        float e_p_exp_0p7 = physics.ExpectedEnergy(theta, 729.f, "ep");
        float e_e_exp_0p7 = physics.ExpectedEnergy(theta, 729.f, "ee");

        float sigma_ep_3p5 = resolution * e_p_exp_3p5 / sqrt(e_p_exp_3p5/1000.f);
        float sigma_ee_3p5 = resolution * e_e_exp_3p5 / sqrt(e_e_exp_3p5/1000.f);
        float sigma_ep_0p7 = resolution * e_p_exp_0p7 / sqrt(e_p_exp_0p7/1000.f);
        float sigma_ee_0p7 = resolution * e_e_exp_0p7 / sqrt(e_e_exp_0p7/1000.f);
        // Fit the energy spectrum with a Gaussian to find the peak position
        TF1 gaus("gaus", "gaus", e_p_exp_3p5 - 2.*sigma_ep_3p5, e_p_exp_3p5 + 2.*sigma_ep_3p5);
        hist_3p5->Fit(&gaus,"RQ0");
        float peak_ep_3p5 = gaus.GetParameter(1);
        
        gaus.SetRange(e_e_exp_3p5 - 2.*sigma_ee_3p5, e_e_exp_3p5 + 2.*sigma_ee_3p5);
        hist_3p5->Fit(&gaus,"RQ0");
        float peak_ee_3p5 = gaus.GetParameter(1);

        gaus.SetRange(e_p_exp_0p7 - 2.*sigma_ep_0p7, e_p_exp_0p7 + 2.*sigma_ep_0p7);
        hist_0p7->Fit(&gaus,"RQ0");
        float peak_ep_0p7 = gaus.GetParameter(1);

        gaus.SetRange(e_e_exp_0p7 - 2.*sigma_ee_0p7, e_e_exp_0p7 + 2.*sigma_ee_0p7);
        hist_0p7->Fit(&gaus,"RQ0");
        float peak_ee_0p7 = gaus.GetParameter(1);

        // make a canvas, measured E vs expected E, and save to output file
        TCanvas *c = new TCanvas(Form("c_mod%d", mod_id), Form("Module %d Non-linearity", mod_id), 800, 600);
        c->SetGrid();
        TGraph *g = new TGraph(4);
        g->SetPoint(0, e_p_exp_3p5, peak_ep_3p5);
        g->SetPoint(1, e_e_exp_3p5, peak_ee_3p5);
        g->SetPoint(2, e_p_exp_0p7, peak_ep_0p7);
        g->SetPoint(3, e_e_exp_0p7, peak_ee_0p7);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.5);
        g->SetTitle(Form("Module %d Non-linearity;Expected Energy (MeV);Measured Peak Position (MeV)", mod_id));
        g->Draw("AP");

        // perfect linearity reference line (y = x)
        double xmin = g->GetXaxis()->GetXmin();
        double xmax = g->GetXaxis()->GetXmax();
        TLine *ref = new TLine(xmin, xmin, xmax, xmax);
        ref->SetLineColor(kRed);
        ref->SetLineStyle(2);
        ref->Draw();

        c->Write();
        delete ref;
        delete c;
        delete g;
    }

    std::cout << "Results saved to " << outFile.GetName() << "\n";
    outFile.Close();

}

void process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal, 
    std::map<int, TH1F*> &energy_hists, PhysicsTools &physics, float Ebeam)
{   
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if( i % 1000 == 0) {
            std::cerr << "Processing event " << i << "/" << tree->GetEntries() << "\n";
        }

        //e-p events selection
        if ( ev.n_clusters == 1) {
            int mod_id = ev.cl_center[0];
            if ( ev.cl_nblocks[0] < 4) continue;
            auto mod = hycal.module_by_id(mod_id);
            if ( !mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals
            float x = mod->x, y = mod->y, z = ev.cl_z[0];
            float theta = std::atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;
            float Eexp = physics.ExpectedEnergy(theta, Ebeam, "ep");
            float energy = ev.cl_energy[0];
            if(std::abs(energy - Eexp) < 3. * resolution * Eexp / sqrt(Eexp/1000.f)) {
                energy_hists[mod_id]->Fill(energy);
            }
        }
        else if ( ev.n_clusters == 2) {
            int mod_id1 = ev.cl_center[0];
            int mod_id2 = ev.cl_center[1];
            float energy1 = ev.cl_energy[0];
            float energy2 = ev.cl_energy[1];
            auto mod1 = hycal.module_by_id(mod_id1);
            auto mod2 = hycal.module_by_id(mod_id2);
            if ( !mod1 || !mod1->is_pwo4() || !mod2 || !mod2->is_pwo4()) continue; // only look at PbWO4 crystals
            float x1 = mod1->x, y1 = mod1->y, z1 = ev.cl_z[0];
            float x2 = mod2->x, y2 = mod2->y, z2 = ev.cl_z[1];
            float theta1 = std::atan2(std::sqrt(x1*x1 + y1*y1), z1) * 180.f / M_PI;
            float theta2 = std::atan2(std::sqrt(x2*x2 + y2*y2), z2) * 180.f / M_PI;
            float phi1 = physics.GetPhiAngle(x1, y1);
            float phi2 = physics.GetPhiAngle(x2, y2);
            float Eexp1 = physics.ExpectedEnergy(theta1, Ebeam, "ee");
            float Eexp2 = physics.ExpectedEnergy(theta2, Ebeam, "ee");

            if(std::abs(energy1 + energy2 - Ebeam) > 4. * resolution * Ebeam / sqrt(Ebeam/1000.f)) continue;
            if(std::abs(energy1 - Eexp1) > 3. * resolution * Eexp1 / sqrt(Eexp1/1000.f)) continue;
            if(std::abs(energy2 - Eexp2) > 3. * resolution * Eexp2 / sqrt(Eexp2/1000.f)) continue;
            if ( std::abs(fabs(phi1 - phi2) - 180.f) > 10.f ) continue; // expect back-to-back in phi

            // if passed all the cuts, fill the 2D histogram for energy
            energy_hists[mod_id1]->Fill(energy1);
            energy_hists[mod_id2]->Fill(energy2);
        }
    }
}