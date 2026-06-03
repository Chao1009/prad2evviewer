
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
#include <TSpectrum.h>
#include <TGraph.h>
#include <TLine.h>
#include <TLatex.h>
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
    std::map<int, TH1F*> &energy_hists, PhysicsTools &physics, float Ebeam, int max_events = -1);

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

    process_event(tree, ev, hycal, energy_hists_3p5, physics, 3485.f, max_events);

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
    process_event(tree2, ev2, hycal, energy_hists_0p7, physics, 729.f, max_events);

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
        // Use TSpectrum to find peak nearest to expected energy, then fit Gaussian
        auto findPeak = [](TH1F *h, double Eexp, double sigma) -> float {
            TSpectrum spec(5);
            spec.Search(h, 3, "nobackground nodraw", 0.05);
            int nfound = spec.GetNPeaks();
            double center = Eexp;
            if (nfound > 0) {
                const double *xpeaks = spec.GetPositionX();
                double best_dist = 9999999.;
                for (int k = 0; k < nfound; ++k) {
                    if(xpeaks[k] < Eexp - 3*sigma || xpeaks[k] > Eexp + 3*sigma) continue;
                    double dist = std::abs(xpeaks[k] - Eexp);
                    if (dist < best_dist) { best_dist = dist; center = xpeaks[k]; }
                }
            }
            TF1 gfit("_gfit_", "gaus", center - 1.5*sigma, center + 1.5*sigma);
            gfit.SetParameters(h->GetBinContent(h->FindBin(center)), center, sigma);
            h->Fit(&gfit, "RQ0");
            return static_cast<float>(gfit.GetParameter(1));
        };

        float peak_ep_3p5 = findPeak(hist_3p5, e_p_exp_3p5, sigma_ep_3p5);
        float peak_ee_3p5 = findPeak(hist_3p5, e_e_exp_3p5, sigma_ee_3p5);
        float peak_ep_0p7 = findPeak(hist_0p7, e_p_exp_0p7, sigma_ep_0p7);
        float peak_ee_0p7 = findPeak(hist_0p7, e_e_exp_0p7, sigma_ee_0p7);

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

        //fitting use 1st order polynomial, anchored at the 3.5 GeV e-p point
        // f(x) = peak_ep_3p5 + slope * (x - e_p_exp_3p5)
        TF1 *fitLine = new TF1("fitLine",
            [](double *x, double *p){ return p[2] + p[0] * (x[0] - p[1]); },
            xmin, xmax, 3);
        fitLine->SetParameter(0, 1.0);
        fitLine->FixParameter(1, e_p_exp_3p5);
        fitLine->FixParameter(2, peak_ep_3p5);
        fitLine->SetLineColor(kBlue);
        fitLine->SetLineWidth(2);
        g->Fit(fitLine, "Q");
        fitLine->Draw("same");

        // draw fit formula and result on the canvas
        double slope = fitLine->GetParameter(0);
        double slope_err = fitLine->GetParError(0);
        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.035);
        tex->DrawLatex(0.15, 0.82, Form("f(x) = peak_{ep,3.5} + slope #times (x - E_{ep,3.5})"));
        tex->DrawLatex(0.15, 0.76, Form("slope = %.4f #pm %.4f", slope, slope_err));
        tex->DrawLatex(0.15, 0.70, Form("E_{ep,3.5} = %.1f MeV,  peak = %.1f MeV",
                                         (double)e_p_exp_3p5, (double)peak_ep_3p5));

        c->Write();
        delete tex;
        delete fitLine;
        delete ref;
        delete c;
        delete g;
    }

    std::cout << "Results saved to " << outFile.GetName() << "\n";
    outFile.Close();

}

void process_event( TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal, 
    std::map<int, TH1F*> &energy_hists, PhysicsTools &physics, float Ebeam, int max_events)
{   
    for (int i = 0; i < tree->GetEntries(); i++) {
        tree->GetEntry(i);
        if( i % 1000 == 0) {
            std::cerr << "Processing event " << i << "/" << tree->GetEntries() << "\r" << std::flush;
        }
        if (max_events > 0 && i >= max_events) {
            std::cerr << "Reached max events limit: " << max_events << "\n";
            break;
        }

        if ( ev.n_clusters == 1) {
            int mod_id = ev.cl_center[0];
            if ( ev.cl_nblocks[0] < 4) continue;
            auto mod = hycal.module_by_id(mod_id);
            if ( !mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals
            float x = mod->x, y = mod->y, z = ev.cl_z[0];
            float theta = std::atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;
            float Eexp = physics.ExpectedEnergy(theta, Ebeam, "ep");
            float energy = ev.cl_energy[0];
            if(std::abs(energy - Eexp) < 4. * resolution * Eexp / sqrt(Eexp/1000.f)) {
                energy_hists[mod_id]->Fill(energy);
            }
        }
        for( int j = 0; j < ev.n_clusters; j++) {
            int mod_id = ev.cl_center[j];
            if (ev.cl_nblocks[j] < 3) continue;
            auto mod = hycal.module_by_id(mod_id);
            if ( !mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals
            float x = mod->x, y = mod->y, z = ev.cl_z[j];
            float theta = std::atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;
            
            float Eexp = physics.ExpectedEnergy(theta, Ebeam, "ee");

            if(Ebeam - Eexp < 3. * Ebeam*resolution/sqrt(Ebeam/1000.f)) continue; // skip if expected energy is too close to beam energy (likely e-p events)

            float energy = ev.cl_energy[j];

            if(std::abs(energy - Eexp) < 5. * resolution * Eexp / sqrt(Eexp/1000.f))
                energy_hists[mod_id]->Fill(energy);
        }
    }
}