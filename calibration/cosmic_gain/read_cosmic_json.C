#include "cosmic_common.h"

void read_cosmic_json(){

    std::string json_file = "cosmic_modules_run570.json";
    std::vector<CosmicEntry> modules;
    if (!read_cosmic_modules(json_file, modules)) {
        std::cerr << "Cannot open json file" << std::endl;
        return;
    }

    TH1D *h_ph_mean = new TH1D("h_ph_mean", ("Peak Height Mean per Module "+json_file+"; Peak Height ADC; Counts").c_str(), 50, 0, 100);
    TH1D *h_pi_mean = new TH1D("h_pi_mean", ("Peak Integral Mean per Module "+json_file+"; Integral ADC; Counts").c_str(), 50, 0, 500);
    for (const auto &m : modules) {
        if (m.count > 0) {
            h_ph_mean->Fill(m.ph_mean);
            h_pi_mean->Fill(m.pi_mean);
        }
    }

    TCanvas *c_peak = new TCanvas("c_peak", "Cosmic Peak per Module", 1400, 500);
    c_peak->Divide(2, 1);

    c_peak->cd(1);
    h_ph_mean->SetLineColor(kBlue);
    h_ph_mean->SetLineWidth(2);
    h_ph_mean->Draw();
    TLine *line_ph = new TLine(35, 0, 35, h_ph_mean->GetMaximum());
    line_ph->SetLineColor(kRed);
    line_ph->SetLineWidth(2);
    line_ph->SetLineStyle(2);
    line_ph->Draw();

    c_peak->cd(2);
    h_pi_mean->SetLineColor(kRed);
    h_pi_mean->SetLineWidth(2);
    h_pi_mean->Draw();
    TLine *line_pi = new TLine(250, 0, 250, h_pi_mean->GetMaximum());
    line_pi->SetLineColor(kBlue);
    line_pi->SetLineWidth(2);
    line_pi->SetLineStyle(2);
    line_pi->Draw();

}
