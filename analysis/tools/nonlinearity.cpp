
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
    std::string pngDir = "module_hists";
    
    int max_events = -1;
    int opt;
    while ((opt = getopt(argc, argv, "a:b:o:n:p:")) != -1) {
        switch (opt) {
            case 'a': input_3p5 = optarg; break;
            case 'b': input_0p7 = optarg; break;
            case 'o': output = optarg; break;
            case 'n': max_events = std::atoi(optarg); break;
            case 'p': pngDir = optarg; break;
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

    TH1F *h_energy_peak_3p5 = new TH1F("h_energy_peak_3p5", "Energy Peak Distribution;Energy (MeV);Counts", 4000, 0, 4000);

    // calculate non-linearity module by module and save to output file
    gSystem->mkdir(pngDir.c_str(), true);
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

        // Find peak with TSpectrum, fit Gaussian, draw fit on current pad, return peak center.
        // Must be called after hist->Draw() with the target pad cd()'d.
        int _fit_uid = mod_id * 10;
        auto fitPeakAndDraw = [&_fit_uid](TH1F *h, double Eexp, double sigma, int color) -> float {
            ++_fit_uid;
            TSpectrum spec(5);
            spec.Search(h, 4, "nobackground nodraw", 0.03);
            int nfound = spec.GetNPeaks();
            int n_in_window = 0;
            double center = 0.;
            if (nfound > 0) {
                const double *xpeaks = spec.GetPositionX();
                for (int k = 0; k < nfound; ++k) {
                    if (xpeaks[k] >= Eexp - 4*sigma && xpeaks[k] <= Eexp + 4*sigma) {
                        ++n_in_window;
                        center = xpeaks[k];
                    }
                }
            }
            if (n_in_window != 1) return 0.f;
            TF1 *gfit = new TF1(Form("_gfit_%d", _fit_uid), "gaus",
                center - 2.*sigma, center + 2.*sigma);
            gfit->SetParameters(h->GetBinContent(h->FindBin(center)), center, sigma);
            h->Fit(gfit, "RQ0", "", center - 2.*sigma, center + 2.*sigma);
            gfit->SetLineColor(color);
            gfit->SetLineWidth(2);
            gfit->Draw("same");
            return static_cast<float>(gfit->GetParameter(1));
        };

        // --- Draw both beam-energy histograms on a two-pad canvas, save PNG ---
        TCanvas *ch = new TCanvas(Form("ch_mod_W%d", mod_id-1000),
            Form("Module W%d Histograms", mod_id-1000), 800, 800);
        ch->Divide(1, 2, 0, 0);

        // --- top pad: 3.5 GeV ---
        ch->cd(1);
        gPad->SetBottomMargin(0.005);
        gPad->SetTopMargin(0.10);
        gPad->SetLeftMargin(0.12);
        hist_3p5->GetXaxis()->SetLabelSize(0);
        hist_3p5->GetXaxis()->SetTitleSize(0);
        hist_3p5->SetTitle(Form("Module W%d;  ;Counts", mod_id-1000));
        hist_3p5->SetLineColor(kBlack);
        hist_3p5->SetLineWidth(2);
        hist_3p5->SetStats(0);
        hist_3p5->Draw("HIST");
        float peak_ep_3p5 = fitPeakAndDraw(hist_3p5, e_p_exp_3p5, sigma_ep_3p5, kRed);
        float peak_ee_3p5 = fitPeakAndDraw(hist_3p5, e_e_exp_3p5, sigma_ee_3p5, kBlue);
        {
            TLatex lat;
            lat.SetNDC(); lat.SetTextSize(0.050);
            lat.SetTextColor(kRed);
            lat.DrawLatex(0.50, 0.86, Form("e-p: exp=%.0f  meas=%s",
                (double)e_p_exp_3p5, peak_ep_3p5 > 0.f ? Form("%.0f MeV", (double)peak_ep_3p5) : "N/A"));
            lat.SetTextColor(kBlue);
            lat.DrawLatex(0.50, 0.78, Form("e-e: exp=%.0f  meas=%s",
                (double)e_e_exp_3p5, peak_ee_3p5 > 0.f ? Form("%.0f MeV", (double)peak_ee_3p5) : "N/A"));
            lat.SetTextColor(kBlack);
            lat.DrawLatex(0.15, 0.86, "E_{beam} = 3.5 GeV");
        }

        // --- bottom pad: 0.7 GeV ---
        ch->cd(2);
        gPad->SetTopMargin(0.005);
        gPad->SetBottomMargin(0.14);
        gPad->SetLeftMargin(0.12);
        hist_0p7->SetTitle(";Energy (MeV);Counts");
        hist_0p7->SetLineColor(kBlack);
        hist_0p7->SetLineWidth(2);
        hist_0p7->SetStats(0);
        hist_0p7->Draw("HIST");
        float peak_ep_0p7 = fitPeakAndDraw(hist_0p7, e_p_exp_0p7, sigma_ep_0p7, kRed);
        float peak_ee_0p7 = 0.f;
        if(729. - e_e_exp_0p7 > 6. * 729.*resolution/sqrt(729./1000.f))
            peak_ee_0p7 = fitPeakAndDraw(hist_0p7, e_e_exp_0p7, sigma_ee_0p7, kBlue);
        {
            TLatex lat;
            lat.SetNDC(); lat.SetTextSize(0.050);
            lat.SetTextColor(kRed);
            lat.DrawLatex(0.50, 0.86, Form("e-p: exp=%.0f  meas=%s",
                (double)e_p_exp_0p7, peak_ep_0p7 > 0.f ? Form("%.0f MeV", (double)peak_ep_0p7) : "N/A"));
            lat.SetTextColor(kBlue);
            lat.DrawLatex(0.50, 0.78, Form("e-e: exp=%.0f  meas=%s",
                (double)e_e_exp_0p7, peak_ee_0p7 > 0.f ? Form("%.0f MeV", (double)peak_ee_0p7) : "N/A"));
            lat.SetTextColor(kBlack);
            lat.DrawLatex(0.15, 0.86, "E_{beam} = 0.7 GeV");
        }

        ch->SaveAs(Form("%s/mod_W%d.png", pngDir.c_str(), mod_id-1000));
        delete ch;

        // if the anchor point (3.5 GeV e-p) has no clean peak, skip this module
        if (peak_ep_3p5 == 0.f) continue;

        h_energy_peak_3p5->Fill(peak_ep_3p5);

        // make a canvas, measured E vs expected E; only add points with valid peaks
        TCanvas *c = new TCanvas(Form("c_mod_W%d", mod_id-1000), Form("Module W%d Non-linearity", mod_id-1000), 1400, 800);
        c->SetGrid();
        TGraph *g = new TGraph();
        int np = 0;
        auto addPoint = [&](double Eexp, float peak) {
            if (peak != 0.f) g->SetPoint(np++, Eexp, peak);
        };
        addPoint(e_p_exp_3p5, peak_ep_3p5);
        addPoint(e_e_exp_3p5, peak_ee_3p5);
        addPoint(e_p_exp_0p7, peak_ep_0p7);
        if(729. - e_e_exp_0p7 > 6. * 729.*resolution/sqrt(729./1000.f)) 
            addPoint(e_e_exp_0p7, peak_ee_0p7);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.5);
        g->SetTitle(Form("Module W%d Non-linearity;Expected Energy (MeV);Measured Peak Position (MeV)", mod_id-1000));
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
    h_energy_peak_3p5->Write();
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
        if ( ev.trigger_bits & ( 1 << 8) == 0) continue;
        if (max_events > 0 && i >= max_events) {
            std::cerr << "Reached max events limit: " << max_events << "\n";
            break;
        }

        if ( ev.n_clusters == 1) {
            int mod_id = ev.cl_center[0];
            if ( ev.cl_nblocks[0] < 3) continue;
            auto mod = hycal.module_by_id(mod_id);
            if ( !mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals
            float x = ev.cl_x[0], y = ev.cl_y[0], z = ev.cl_z[0];
            float theta = std::atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;
            float Eexp = physics.ExpectedEnergy(theta, Ebeam, "ep");
            float energy = ev.cl_energy[0];
            if(Ebeam < 1000. && theta < 1.1) continue;
            if(std::abs(energy - Eexp) < 5. * resolution * Eexp / sqrt(Eexp/1000.f)) {
                energy_hists[mod_id]->Fill(energy);
            }
        }
        if ( ev.n_clusters == 2) {
            int mod_id_1 = ev.cl_center[0];
            int mod_id_2 = ev.cl_center[1];
            auto mod1 = hycal.module_by_id(mod_id_1);
            auto mod2 = hycal.module_by_id(mod_id_2);
            if ( !mod1 || !mod1->is_pwo4() || !mod2 || !mod2->is_pwo4()) continue; // only look at PbWO4 crystals
            float x1 = ev.cl_x[0], y1 = ev.cl_y[0], z1 = ev.cl_z[0];
            float x2 = ev.cl_x[1], y2 = ev.cl_y[1], z2 = ev.cl_z[1];
            float theta1 = std::atan2(std::sqrt(x1*x1 + y1*y1), z1) * 180.f / M_PI;
            float theta2 = std::atan2(std::sqrt(x2*x2 + y2*y2), z2) * 180.f / M_PI;
            float Eexp1 = physics.ExpectedEnergy(theta1, Ebeam, "ee");
            float Eexp2 = physics.ExpectedEnergy(theta2, Ebeam, "ee");
            float energy1 = ev.cl_energy[0];
            float energy2 = ev.cl_energy[1];
            float phi1 = physics.GetPhiAngle(x1, y1);
            float phi2 = physics.GetPhiAngle(x2, y2);
            float phi_diff = fabs(phi1 - phi2) - 180.f;
            if (phi_diff < -10.f || phi_diff > 10.f) continue;
            if ()
            if(std::abs(energy1 + energy2 - Ebeam) < 4. * resolution * Ebeam / sqrt(Ebeam/1000.f)) {
                if(std::abs(energy1 - Eexp1) < 4. * resolution * Eexp1 / sqrt(Eexp1/1000.f) &&
                   std::abs(energy2 - Eexp2) < 4. * resolution * Eexp2 / sqrt(Eexp2/1000.f)) {
                    energy_hists[mod_id_1]->Fill(energy1);
                    energy_hists[mod_id_2]->Fill(energy2);
                }
            }
        }

        /*for( int j = 0; j < ev.n_clusters; j++) {
            int mod_id = ev.cl_center[j];
            if (ev.cl_nblocks[j] < 3) continue;
            auto mod = hycal.module_by_id(mod_id);
            if ( !mod || !mod->is_pwo4()) continue; // only look at PbWO4 crystals
            float x = mod->x, y = mod->y, z = ev.cl_z[j];
            float theta = std::atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;

            if (theta < 3.3) continue;
            
            float Eexp = physics.ExpectedEnergy(theta, Ebeam, "ee");

            if(Ebeam - Eexp < 6. * Ebeam*resolution/sqrt(Ebeam/1000.f)) continue; // skip if expected energy is too close to beam energy (likely e-p events)

            float energy = ev.cl_energy[j];

            if(std::abs(energy - Eexp) < 5. * resolution * Eexp / sqrt(Eexp/1000.f))
                energy_hists[mod_id]->Fill(energy);
        }*/
    }
}