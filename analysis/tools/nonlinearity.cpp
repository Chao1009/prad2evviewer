
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "ToolUtils.h"

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

#include <iostream>
#include <string>
#include <vector>
#include <array>
#include <map>
#include <utility>
#include <cmath>
#include <cstdlib>
#include <algorithm>
#include <mutex>

using namespace analysis;

using EventVars_Recon = prad2::ReconEventData;

// returns the number of events passing the sum-trigger selection
long long process_event( bool use_GEM, TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal,
    std::map<int, TH1F*> &energy_hists, float Ebeam, int max_events = -1,
    const std::string &label = "", std::mutex *io_mtx = nullptr);

float resolution = 0.035; // pre-defined energy resolution

// One beam-energy data set (-a/-b/-c): its inputs, per-module energy
// histograms and output directory.
struct BeamSet {
    const char *tag;    // histogram/directory suffix
    const char *label;  // beam energy in GeV, for titles and log lines
    float Ebeam;        // MeV
    float min_ee_sep;   // fit the e-e peak only if it is this far (MeV) from the e-p peak
    std::vector<std::string> inputs;
    std::map<int, TH1F*> hists;
    long long n_sum = 0;
    TDirectory *dir = nullptr;
};
constexpr int kNBeams = 3;

bool Vetoed(float cl_time, float sci_time, float sci_int){
    // Simple veto logic: if the cluster time is within a certain window of the scintillator time, and the scintillator signal is above a threshold, we consider it a vetoed event.
    const float time_shift = 35.f; // ns
    const float time_window = 7.f; // ns
    const float int_threshold = 2000.f; // arbitrary units
    return (fabs(cl_time - sci_time - time_shift) < time_window) && (sci_int > int_threshold);
}

int main(int argc, char *argv[]){
    analysis::InitRootThreading();

    std::array<BeamSet, kNBeams> beams{{
        {"3p5", "3.5", 3485.41f, 0.f},
        {"2p2", "2.2", 2239.51f, 0.f},
        {"0p7", "0.7", 728.9f, 170.f},
    }};

    std::string output;
    std::string pngDir = "module_hists";

    int max_events = -1;
    bool use_GEM = false;
    {
        std::vector<std::string> *cur = nullptr;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "-a")       { cur = &beams[0].inputs; }
            else if (arg == "-b")  { cur = &beams[1].inputs; }
            else if (arg == "-c")  { cur = &beams[2].inputs; }
            else if (arg == "-g")  { cur = nullptr; use_GEM = true; }
            else if (arg == "-o")  { cur = nullptr; if (++i < argc) output = argv[i]; }
            else if (arg == "-n")  { cur = nullptr; if (++i < argc) max_events = std::atoi(argv[i]); }
            else if (arg == "-p")  { cur = nullptr; if (++i < argc) pngDir = argv[i]; }
            else if (!arg.empty() && arg[0] != '-' && cur) { cur->push_back(arg); }
        }
    }
    

    std::string dbDir = prad2::database_dir();

    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);

    // Energy histogram for each crystal
    for (int i = 0; i < hycal.module_count(); ++i) {
        const auto &m = hycal.module(i);
        if (!m.is_pwo4()) continue;
        std::string hname = "h_energy_" + m.name;
        std::string htitle = "Energy " + m.name + ";E (MeV);Counts";
        for (auto &beam : beams) {
            TH1F *h = new TH1F((hname + "_" + beam.tag).c_str(),
                               (htitle + " (" + beam.label + ")").c_str(), 400, 0, 4000);
            h->SetDirectory(nullptr);
            beam.hists[m.id] = h;
        }
    }

    std::mutex io_mtx;

    auto run_energy = [&](BeamSet &beam) {
        const std::string label = std::string(beam.label) + " GeV";
        TChain chain("recon");
        for (const auto &f : beam.inputs) {
            chain.Add(f.c_str());
            std::lock_guard<std::mutex> lk(io_mtx);
            std::cerr << "[" << label << "] Added file: " << f << "\n";
        }

        EventVars_Recon ev;
        prad2::SetReconReadBranches(&chain, ev);

        {
            std::lock_guard<std::mutex> lk(io_mtx);
            std::cerr << "[" << label << "] Processing "
                      << chain.GetEntries() << " event(s)\n";
        }
        beam.n_sum = process_event(use_GEM, &chain, ev, hycal, beam.hists,
                                   beam.Ebeam, max_events, label, &io_mtx);
        {
            std::lock_guard<std::mutex> lk(io_mtx);
            std::cerr << "[" << label << "] Selected " << beam.n_sum
                      << " sum-trigger event(s)\n";
        }
    };
    ParallelFor(beams.size(), kNBeams, [&](size_t b, int) { run_energy(beams[b]); });

    // a file with no sum-trigger events would silently produce an all-default
    // calibration (every module skipped); refuse to write output in that case
    if (std::any_of(beams.begin(), beams.end(), [](const BeamSet &b) { return b.n_sum <= 0; })) {
        std::cerr << "No sum-trigger events selected (";
        for (int b = 0; b < kNBeams; ++b)
            std::cerr << (b ? ", " : "") << beams[b].label << " GeV: " << beams[b].n_sum;
        std::cerr << "); check trigger_bits in the inputs. Aborting before writing calibration.\n";
        return 1;
    }

    std::string calib_file = dbDir + "/" + "calibration/calibration_factor_3p5_June1.json";
    int nmatched = hycal.LoadCalibration(calib_file);
    if (nmatched >= 0)
        std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";

    TH1F *h_energy_peak_3p5 = new TH1F("h_energy_peak_3p5", "Energy Peak Distribution;Energy (MeV);Counts", 4000, 0, 4000);

    // calculate non-linearity module by module and save to output file
    gSystem->mkdir(pngDir.c_str(), true);
    TFile outFile(output.empty() ? "nonlinearity_results.root" : output.c_str(), "RECREATE");
    for (auto &beam : beams) beam.dir = outFile.mkdir(Form("energy_%sGeV", beam.tag));
    TDirectory *dir_lin = outFile.mkdir("linearity");
    for (int i = 0; i < hycal.module_count(); i++) {
        const auto &mod = hycal.module(i);
        int mod_id = mod.id;
        if (!mod.is_pwo4()) continue;

        float x = mod.x, y = mod.y, z = 6275.f;
        float theta = std::atan2(std::sqrt(x*x + y*y), z) * 180.f / M_PI;

        int _fit_uid = mod_id * 10;

        // Peak fit in [Eexp±6σ]: the weighted centroid seeds a two-stage
        // Gaussian fit (mean±2σ, then mean±1σ); returns the final fitted mean.
        auto fitPeakAndDraw = [&](TH1F *h, double Eexp, double sigma, int color) -> float {
            ++_fit_uid;
            int b0 = std::max(1, h->FindBin(Eexp - 6.*sigma));
            int b1 = std::min(h->GetNbinsX(), h->FindBin(Eexp + 6.*sigma));
            double ypad_min = 0;
            double ypad = h->GetMaximum() * 1.2;
            // search-window: dashed, full y-axis height
            TLine *wl = new TLine(Eexp - 6.*sigma, ypad_min, Eexp - 6.*sigma, ypad);
            wl->SetLineColor(6); wl->SetLineStyle(7); wl->SetLineWidth(2); wl->Draw();
            TLine *wr = new TLine(Eexp + 6.*sigma, ypad_min, Eexp + 6.*sigma, ypad);
            wr->SetLineColor(6); wr->SetLineStyle(7); wr->SetLineWidth(2); wr->Draw();
            // weighted mean within search window
            double wsum = 0., wpos = 0.;
            for (int ib = b0; ib <= b1; ++ib) {
                double c = h->GetBinContent(ib);
                if (c > 0.) { wsum += c; wpos += c * h->GetBinCenter(ib); }
            }
            if (wsum <= 0.) return 0.f;
            double mean = wpos / wsum;
            // first Gaussian fit within [mean ± 2σ]
            auto estimateSig = [](double E) -> double {
                return (E > 0.) ? E * 0.035 / std::sqrt(E / 1000.) : 1.;
            };
            double sig = estimateSig(mean);
            {
                TF1 g1(Form("_fpk_g1_%d", _fit_uid), "gaus", mean - 2.*sig, mean + 2.*sig);
                g1.SetParameters(h->GetMaximum(), mean, sig);
                h->Fit(&g1, "RQ0");
                mean = g1.GetParameter(1);
            }
            // final Gaussian fit within [mean ± 1σ]
            sig = estimateSig(mean);
            TF1 *g = new TF1(Form("_gfit_%d", _fit_uid), "gaus", mean - sig, mean + sig);
            g->SetParameters(h->GetMaximum(), mean, sig);
            h->Fit(g, "RQ0");
            mean = g->GetParameter(1);
            double fit_amp = g->GetParameter(0);
            g->SetLineColor(color); g->SetLineWidth(3); g->Draw("same");
            // draw a marker at the fitted mean, at Gaussian peak height
            TMarker *mk = new TMarker(mean, fit_amp, 29); // star symbol
            mk->SetMarkerColor(color); mk->SetMarkerSize(2.5); mk->Draw();
            return static_cast<float>(mean);
        };

        // --- one pad per beam energy (3.5, 2.2, 0.7 GeV from top), save PNG ---
        TCanvas *ch = new TCanvas(Form("ch_mod_W%d", mod_id-1000),
            Form("Module W%d Histograms", mod_id-1000), 800, 1200);
        ch->Divide(1, kNBeams, 0, 0);

        float exp_ep[kNBeams], exp_ee[kNBeams], peak_ep[kNBeams] = {}, peak_ee[kNBeams] = {};
        for (int b = 0; b < kNBeams; ++b) {
            const BeamSet &beam = beams[b];
            const bool top = (b == 0), bottom = (b == kNBeams - 1);
            TH1F *h = beam.hists.at(mod_id);
            exp_ep[b] = physics.ExpectedEnergy(theta, beam.Ebeam, "ep");
            exp_ee[b] = physics.ExpectedEnergy(theta, beam.Ebeam, "ee");
            float sigma_ep = resolution * exp_ep[b] / sqrt(exp_ep[b]/1000.f);
            float sigma_ee = resolution * exp_ee[b] / sqrt(exp_ee[b]/1000.f);

            ch->cd(b + 1);
            gPad->SetTopMargin(top ? 0.10 : 0.005);
            gPad->SetBottomMargin(bottom ? 0.14 : 0.005);
            gPad->SetLeftMargin(0.12);
            if (!bottom) {
                h->GetXaxis()->SetLabelSize(0);
                h->GetXaxis()->SetTitleSize(0);
            }
            h->SetTitle(top ? Form("Module W%d;  ;Counts", mod_id-1000)
                            : bottom ? ";Energy (MeV);Counts" : ";  ;Counts");
            h->SetLineColor(kBlack);
            h->SetLineWidth(2);
            h->SetStats(0);
            h->Draw("HIST");
            peak_ep[b] = fitPeakAndDraw(h, exp_ep[b], sigma_ep, kRed);
            if (std::abs(exp_ep[b] - exp_ee[b]) > beam.min_ee_sep)
                peak_ee[b] = fitPeakAndDraw(h, exp_ee[b], sigma_ee, kBlue);

            TLatex lat;
            lat.SetNDC(); lat.SetTextSize(0.050);
            lat.SetTextColor(kRed);
            lat.DrawLatex(0.50, 0.86, Form("e-p: exp=%.0f  meas=%s",
                (double)exp_ep[b], peak_ep[b] > 0.f ? Form("%.0f MeV", (double)peak_ep[b]) : "N/A"));
            lat.SetTextColor(kBlue);
            lat.DrawLatex(0.50, 0.78, Form("e-e: exp=%.0f  meas=%s",
                (double)exp_ee[b], peak_ee[b] > 0.f ? Form("%.0f MeV", (double)peak_ee[b]) : "N/A"));
            lat.SetTextColor(kBlack);
            lat.DrawLatex(0.15, 0.86, Form("E_{beam} = %s GeV", beam.label));
        }

        ch->SaveAs(Form("%s/mod_W%d.png", pngDir.c_str(), mod_id-1000));
        delete ch;

        // if the anchor point (3.5 GeV e-p) has no clean peak, skip this module
        if (peak_ep[0] == 0.f) continue;
        const float E_base = exp_ep[0];

        h_energy_peak_3p5->Fill(peak_ep[0]);

        // (E_exp, E_rec) of every valid peak, the anchor first
        std::vector<std::pair<double, float>> points;
        for (int b = 0; b < kNBeams; ++b) {
            if (peak_ep[b] != 0.f) points.emplace_back(exp_ep[b], peak_ep[b]);
            if (peak_ee[b] != 0.f) points.emplace_back(exp_ee[b], peak_ee[b]);
        }

        // make a canvas, E_rec/E_exp vs E_rec
        TCanvas *c = new TCanvas(Form("c_mod_W%d", mod_id-1000), Form("Module W%d Non-linearity", mod_id-1000), 1400, 800);
        c->SetGrid();
        TGraph *g = new TGraph();
        for (const auto &[Eexp, peak] : points) g->SetPoint(g->GetN(), peak, peak / Eexp);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(1.5);
        g->SetTitle(Form("Module W%d Non-linearity;E_{rec} (MeV);E_{rec}/E_{exp}", mod_id-1000));
        g->Draw("AP");
        g->GetYaxis()->SetRangeUser(0.85, 1.15);
        g->Draw("AP");

        // perfect linearity reference line (y = 1)
        double xmin = g->GetXaxis()->GetXmin();
        double xmax = g->GetXaxis()->GetXmax();
        TLine *ref = new TLine(xmin, 1.0, xmax, 1.0);
        ref->SetLineColor(kRed);
        ref->SetLineStyle(2);
        ref->SetLineWidth(2);
        ref->Draw();

        // 1st order fit: E_rec/E_exp = 1 + nl1 * (E_rec - E_base)/1000
        TF1 *fitLine = new TF1("fitLine",
            [](double *x, double *p){ return 1.0 + p[0] * (x[0] - p[1])/1000.0; },
            xmin, xmax, 2);
        fitLine->SetParameter(0, 0.01);
        fitLine->FixParameter(1, E_base);
        fitLine->SetLineColor(kBlue);
        fitLine->SetLineWidth(2);
        g->Fit(fitLine, "RQ0");
        fitLine->Draw("same");

        double nl        = fitLine->GetParameter(0);
        double nl_err    = fitLine->GetParError(0);
        double chi2      = fitLine->GetChisquare();
        int    ndf       = fitLine->GetNDF();

        // 2nd order fit: E_rec/E_exp = 1 + nl1*(E_rec-E_base)/1000 + nl2*((E_rec-E_base)/1000)^2
        TF1 *fitLine2 = new TF1("fitLine2",
            [](double *x, double *p){
                double t = (x[0] - p[2]) / 1000.0;
                return 1.0 + p[0] * t + p[1] * t * t;
            },
            xmin, xmax, 3);
        fitLine2->SetParameter(0, nl);
        fitLine2->SetParameter(1, 0.0);
        fitLine2->FixParameter(2, E_base);
        fitLine2->SetLineColor(kMagenta+1);
        fitLine2->SetLineWidth(2);
        fitLine2->SetLineStyle(7);
        g->Fit(fitLine2, "RQ0+");
        fitLine2->Draw("same");

        double nl2_1     = fitLine2->GetParameter(0);
        double nl2_1_err = fitLine2->GetParError(0);
        double nl2_2     = fitLine2->GetParameter(1);
        double nl2_2_err = fitLine2->GetParError(1);
        double chi2_2    = fitLine2->GetChisquare();
        int    ndf_2     = fitLine2->GetNDF();

        TLatex *tex = new TLatex();
        tex->SetNDC();
        tex->SetTextSize(0.030);
        tex->SetTextColor(kBlue);
        tex->DrawLatex(0.15, 0.88, "1st: E_{rec}/E_{exp} = 1 + nl_{1} #times (E_{rec}-E_{base})/1000");
        tex->DrawLatex(0.15, 0.83, Form("nl_{1} = %.4f #pm %.4f, #chi^{2}/ndf = %.2f/%d", nl, nl_err, chi2, ndf));
        tex->SetTextColor(kMagenta+1);
        tex->DrawLatex(0.15, 0.78, "2nd: + nl_{2} #times ((E_{rec}-E_{base})/1000)^{2}");
        tex->DrawLatex(0.15, 0.73, Form("nl_{1} = %.4f #pm %.4f, nl_{2} = %.4f #pm %.4f", nl2_1, nl2_1_err, nl2_2, nl2_2_err));
        tex->DrawLatex(0.15, 0.68, Form("#chi^{2}/ndf = %.2f/%d,  E_{base} = %.1f MeV", chi2_2, ndf_2, (double)E_base));
        tex->SetTextColor(kBlack);

        // corrected points using 1st order: E_corr = E_rec / (1 + nl*(E_rec-E_base)/1000)
        TGraph *gCorr = new TGraph();
        gCorr->SetMarkerStyle(24);
        gCorr->SetMarkerSize(1.0);
        gCorr->SetMarkerColor(kGreen+2);
        for (const auto &[Eexp, peak] : points) {
            double denom = 1.0 + nl * (peak - E_base)/1000.0;
            double E_corr = (denom != 0.) ? peak / denom : peak;
            gCorr->SetPoint(gCorr->GetN(), E_corr, E_corr / Eexp);
        }
        gCorr->Draw("P same");

        // corrected points using 2nd order fit
        TGraph *gCorr2 = new TGraph();
        gCorr2->SetMarkerStyle(25);
        gCorr2->SetMarkerSize(1.0);
        gCorr2->SetMarkerColor(kOrange+2);
        for (const auto &[Eexp, peak] : points) {
            double t = (peak - E_base) / 1000.0;
            double denom = 1.0 + nl2_1 * t + nl2_2 * t * t;
            double E_corr = (denom != 0.) ? peak / denom : peak;
            gCorr2->SetPoint(gCorr2->GetN(), E_corr, E_corr / Eexp);
        }
        gCorr2->Draw("P same");

        TLegend *leg = new TLegend(0.60, 0.12, 0.92, 0.52);
        leg->SetBorderSize(1);
        leg->SetTextSize(0.026);
        leg->AddEntry(g,        "Measured points",         "p");
        leg->AddEntry(fitLine,  "1st order fit",           "l");
        leg->AddEntry(fitLine2, "2nd order fit",           "l");
        leg->AddEntry(ref,      "Perfect linearity (y=1)", "l");
        leg->AddEntry(gCorr,    "Corrected points (1st)",  "p");
        leg->AddEntry(gCorr2,   "Corrected points (2nd)",  "p");
        leg->Draw();

        dir_lin->cd();
        c->Write();
        outFile.cd();
        delete leg;
        delete gCorr2;
        delete gCorr;
        delete tex;
        delete fitLine2;
        delete fitLine;
        delete ref;
        delete c;
        delete g;

        // outermost ring or absorber/beam-hole region: no non-linearity correction
        if (fdec::test_bit(mod.flag, fdec::kTransition) || fdec::test_bit(mod.flag, fdec::kInnerBound))
            nl = 0.0;
        hycal.SetCalibNonLinearity(mod_id, nl);
    }
    for (auto &beam : beams) {
        beam.dir->cd();
        for (auto &[id, h] : beam.hists) if (h) h->Write();
    }
    outFile.cd();
    h_energy_peak_3p5->Write();
    std::cout << "Results saved to " << outFile.GetName() << "\n";
    outFile.Close();

    hycal.PrintCalibConstants("new_calibration_NonLinearity.json");

}

long long process_event(bool use_GEM, TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal,
    std::map<int, TH1F*> &energy_hists, float Ebeam, int max_events,
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
        if( i % 1000 == 0) {
            log_msg("[" + label + "] Processing event " + std::to_string(i)
                    + "/" + std::to_string(tree->GetEntries()) + "\r", true);
        }
        if ((ev.trigger_bits & prad2::TBIT_sum) == 0) continue;
        n_accepted++;

        for( int j = 0; j < ev.n_clusters; j++) {
            int mod_id = ev.cl_center[j];
            if (ev.cl_nblocks[j] < 4) continue;
            auto mod = hycal.module_by_id(mod_id);
            if ( !mod || !mod->is_pwo4()) continue;

            float c_x, c_y, c_z;
            if(!use_GEM){
                c_x = ev.cl_x[j];
                c_y = ev.cl_y[j];
                c_z = ev.cl_z[j];
            }
            else{
                bool match[4] = {false, false, false, false};
                for(int d = 0; d < 4; d++){
                    if(ev.matchFlag[j] & 1 << d) match[d] = true;
                }
                if( (match[0] || match[1]) && (match[2] || match[3]) ){
                    if(match[0]) ev.first_match(j, 0, c_x, c_y, c_z);
                    else ev.first_match(j, 1, c_x, c_y, c_z);
                    //projection onto the HyCal module plane
                    float scale = 6275.f / c_z;
                    c_x *= scale;
                    c_y *= scale;
                    c_z = 6275.f;
                }
                else{
                    c_x = -999.f;
                    c_y = -999.f;
                    c_z = -999.f;
                }
                
            }

            // require the hit near the seed module centre (|xd|,|yd| < 0.2 module sizes)
            const auto [xd, yd] = mod->cell_offset<float>(c_x, c_y);
            if (std::abs(xd) >= 0.2f || std::abs(yd) >= 0.2f) continue;

            float theta = std::atan2(std::sqrt(c_x*c_x + c_y*c_y), 6275.f) * 180.f / M_PI;
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
