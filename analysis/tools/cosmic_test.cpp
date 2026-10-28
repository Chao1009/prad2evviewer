// cosmic_test.cpp: per-module HyCal cosmic and LMS peak spectra.
//
// Usage: cosmic_test -r <run> -n <file_count> [-j existing_json] [-L]
//
// Reads the `events` tree of prad_023<run>.000NN_raw.root, NN = 0..file_count-1,
// in the current directory.  Cosmic events (5-70 channels) and LMS events
// (more than 900 channels) fill single-peak integral and height spectra of
// every PbWO4 module and of the lead-glass ring; each spectrum is fitted with
// a Gaussian around its maximum (plots in ./fit_canvas/).  -L analyses the
// lead-glass ring alone, from single-channel events.
//
// Outputs: cosmic_run[_LG]_<run>.root, cosmic_peak_<run>.dat,
// cosmic_eventNum_<run>.dat, lms_run_<run>.json (not with -L) and
// cosmic_modules_run<run>.json, or with -j this run's entries appended to
// every module of existing_json.

#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"

#include <TCanvas.h>
#include <TChain.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <unistd.h>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

using EventVars = prad2::RawEventData;

namespace {

constexpr int kNW = 1156;  // PbWO4 modules W1..W1156, module id 1000 + W

// Lead-glass ring around the PbWO4 region, in output order.
constexpr int kNLG = 76;
constexpr int kLGModuleIds[kNLG] = {
    156, 157, 158, 159, 160, 161, 162, 163, 164, 165, 166, 167, 168, 169, 170, 171, 172, 173, 174,
    186, 216, 246, 276, 306, 336, 366, 396, 426, 456, 486, 516, 546, 576, 606, 636, 666, 696, 726,
    175, 205, 235, 265, 295, 325, 355, 385, 415, 445, 475, 505, 535, 565, 595, 625, 655, 685, 715,
    727, 728, 729, 730, 731, 732, 733, 734, 735, 736, 737, 738, 739, 740, 741, 742, 743, 744, 745
};

int lgIndex(int module_id)
{
    for (int k = 0; k < kNLG; ++k)
        if (kLGModuleIds[k] == module_id) return k;
    return -1;
}

enum Spectrum { kIntegral, kHeight, kLmsIntegral, kLmsHeight, kNSpectra };
const char *const kHistPrefix[kNSpectra] = {"peak_", "peakHeight_", "lms_", "lms_height_"};
const char *const kPlotPrefix[kNSpectra] = {"", "peakHeight_", "LMS_", "LMS_height_"};

struct PeakFit { float mean, sigma; };

struct Module {
    int id;      // HyCal module id
    int number;  // number in the module name: W number or G id
    TH1F *h[kNSpectra];
    PeakFit fit[kNSpectra];
};

struct ModuleGroup {
    char key;           // name prefix, 'W' or 'G'
    std::string tag;    // histogram / plot name tag
    const char *dir;    // output directory of the cosmic spectra
    bool shoulder_cut;  // see fitPeak; applied to the cosmic spectra
    TH2F *hit_map;
    std::vector<Module> modules;
};

void addModule(ModuleGroup &g, int id, int number)
{
    const int nbins[kNSpectra] = {80, 80, 80*4, 80*4};
    const double hi[kNSpectra] = {800, 200, 800*4, 200*4};
    Module m{id, number, {}, {}};
    for (int s = 0; s < kNSpectra; ++s) {
        const std::string name = kHistPrefix[s] + g.tag + "module_" + std::to_string(number);
        m.h[s] = new TH1F(name.c_str(), name.c_str(), nbins[s], 0, hi[s]);
    }
    g.modules.push_back(m);
}

// Gaussian fit in [0.7, 1.5] x the centre of the maximum bin, drawn to png;
// {0.1, 1e5} for an empty histogram or a failed fit.  shoulder_cut first
// empties bins 1-8 when either of the two bins below the maximum is under
// 50% / 40% of it.
PeakFit fitPeak(TH1F *h, bool shoulder_cut, const std::string &png)
{
    if (h->GetEntries() <= 0) return {0.1f, 1e5f};
    if (shoulder_cut) {
        const double top = h->GetBinContent(h->GetMaximumBin());
        const float max_l1 = h->GetBinContent(h->GetMaximumBin() - 1);
        const float max_l2 = h->GetBinContent(h->GetMaximumBin() - 2);
        if (max_l1 < 0.5 * top || max_l2 < 0.4 * top)
            for (int b = 1; b <= 8; b++) h->SetBinContent(b, 0);
    }
    const float max = h->GetBinCenter(h->GetMaximumBin());
    h->Fit("gaus", "Q", "r", max*0.7, max*1.5);
    TF1 *fit = h->GetFunction("gaus");
    if (!fit) return {0.1f, 1e5f};
    const PeakFit result{static_cast<float>(fit->GetParameter(1)),
                         static_cast<float>(fit->GetParameter(2))};
    TCanvas *c = new TCanvas();
    h->Draw();
    fit->Draw("same");
    c->SaveAs(png.c_str());
    delete c;
    return result;
}

void fitGroup(ModuleGroup &g, Spectrum s)
{
    const bool shoulder_cut = g.shoulder_cut && (s == kIntegral || s == kHeight);
    for (auto &m : g.modules)
        m.fit[s] = fitPeak(m.h[s], shoulder_cut, std::string("./fit_canvas/fit_") + kPlotPrefix[s]
                           + g.tag + "module_" + std::to_string(m.number) + ".png");
}

} // namespace

int main(int argc, char *argv[])
{
    std::string in_json;
    int run_number = -1, file_number = -1;
    bool lg_only = false;
    int opt;
    while ((opt = getopt(argc, argv, "r:n:j:L")) != -1) {
        switch (opt) {
            case 'r': run_number  = std::atoi(optarg); break;
            case 'n': file_number = std::atoi(optarg); break;
            case 'j': in_json     = optarg; break;
            case 'L': lg_only     = true; break;
            default:
                std::cerr << "Usage: " << argv[0] << " [-r run_number] [-n file_number] [-j existing_json] [-L]\n";
                return 1;
        }
    }

    if (run_number < 0 || file_number < 0) {
        std::cerr << "Usage: " << argv[0]
                  << " -r <run_number> -n <file_count> [-j existing_json] [-L]\n";
        return 1;
    }

    TChain *tree = new TChain("events");
    for (int i = 0; i <= file_number - 1; i++) {
        std::string filename = Form("prad_023%d.000%02d_raw.root", run_number, i);
        tree->Add(filename.c_str());
    }
    auto ev = std::make_unique<EventVars>();
    prad2::SetRawReadBranches(tree, *ev);

    fdec::HyCalSystem hycal;
    std::string db_dir = prad2::database_dir();
    hycal.Init(db_dir + "/hycal_map.json");

    ModuleGroup pwo{'W', "", "peak_histograms", true, nullptr, {}};
    ModuleGroup lg{'G', "LG_", "peak_histograms_LG", false, nullptr, {}};
    if (!lg_only)
        for (int n = 1; n <= kNW; n++) addModule(pwo, 1000 + n, n);
    for (int id : kLGModuleIds) addModule(lg, id, id);
    if (!lg_only)
        pwo.hit_map = new TH2F("cosmic_eventNum", "Cosmic Event Number", 34, -17.*20.75, 17.*20.75, 34, -17.*20.75, 17.*20.75);
    lg.hit_map = new TH2F("cosmic_eventNum_LG", "Cosmic Event Number for LG Modules", 34, -17.*38.15, 17.*38.15, 34, -17.*38.15, 17.*38.15);
    std::vector<ModuleGroup *> groups;
    if (!lg_only) groups.push_back(&pwo);
    groups.push_back(&lg);

    // Cosmic hits per module id: every hit with a peak, or with -L only the
    // single-peak lead-glass hits.
    int event_num_module[3000] = {};

    int nentries = tree->GetEntries();
    for (int i = 0; i < nentries; i++) {
        tree->GetEntry(i);
        std::cout << "Event " << ev->event_num << ": nch = " << ev->nch << "\r" << std::flush;
        const bool lms    = !lg_only && ev->nch > 900;
        const bool cosmic = lg_only ? ev->nch <= 1 : (ev->nch > 4 && ev->nch <= 70);
        if (!lms && !cosmic) continue;
        for (int j = 0; j < ev->nch; j++) {
            const auto *mod = hycal.module_by_id(ev->module_id[j]);
            if (!mod || !mod->is_hycal()) continue;
            if (ev->npeaks[j] <= 0) continue;
            if (cosmic && !lg_only) event_num_module[mod->id]++;
            if (ev->npeaks[j] != 1) continue;

            const int w = mod->id - 1000;
            const bool is_pwo = !lg_only && w >= 1 && w <= kNW;
            const int k = lgIndex(mod->id);
            if (!is_pwo && k < 0) continue;
            if (lg_only) event_num_module[mod->id]++;
            ModuleGroup &g = is_pwo ? pwo : lg;
            Module &m = g.modules[is_pwo ? w - 1 : k];
            if (lms) {
                m.h[kLmsIntegral]->Fill(ev->peak_integral[j][0]);
                m.h[kLmsHeight]->Fill(ev->peak_height[j][0]);
            } else {
                m.h[kIntegral]->Fill(ev->peak_integral[j][0]);
                m.h[kHeight]->Fill(ev->peak_height[j][0]);
                g.hit_map->Fill(mod->x, mod->y);
            }
        }
    }

    TFile outfile(Form("cosmic_run_%s%d.root", lg_only ? "LG_" : "", run_number), "RECREATE");
    outfile.cd();
    for (auto *g : groups) {
        outfile.mkdir(g->dir)->cd();
        for (const auto &m : g->modules) {
            if (m.h[kIntegral]->GetEntries() > 0) m.h[kIntegral]->Write();
            if (m.h[kHeight]->GetEntries() > 0) m.h[kHeight]->Write();
        }
        outfile.cd();
    }
    for (auto *g : groups)
        if (g->hit_map->GetEntries() > 0) g->hit_map->Write();

    // The spectra are written unfitted above.  Keep this fit order: each fit
    // starts from the parameter errors the previous one left in ROOT's global
    // "gaus" function.
    for (auto *g : groups) {
        fitGroup(*g, kIntegral);
        fitGroup(*g, kHeight);
    }
    for (Spectrum s : {kLmsIntegral, kLmsHeight})
        for (auto *g : groups) fitGroup(*g, s);

    TH1F *peak_module = new TH1F("peak_module", "Peak Integral by Module", 100, 0, 500);
    TH1F *rms_module = new TH1F("rms_module", "RMS of Peak Integral by Module", 100, 0, 400);
    for (auto *g : groups)
        for (const auto &m : g->modules) {
            peak_module->Fill(m.fit[kIntegral].mean);
            rms_module->Fill(m.fit[kIntegral].sigma);
        }
    peak_module->Write();
    rms_module->Write();

    std::ofstream csv_out(Form("cosmic_peak_%d.dat", run_number));
    csv_out << "ModuleID  PeakIntegral  RMS\n";
    for (auto *g : groups)
        for (const auto &m : g->modules)
            csv_out << g->key << m.number << "  " << m.fit[kIntegral].mean << "  " << m.fit[kIntegral].sigma << "\n";
    csv_out.close();

    std::ofstream rate_out(Form("cosmic_eventNum_%d.dat", run_number));
    rate_out << "ModuleID  EventCount\n";
    for (auto *g : groups)
        for (const auto &m : g->modules)
            rate_out << g->key << m.number << "  " << event_num_module[m.id] << "\n";
    rate_out.close();

    // ── JSON output ───────────────────────────────────────────────────────
    if (run_number > 0) {
        auto make_entry = [&](const Module &m) -> std::string {
            char buf[512];
            std::snprintf(buf, sizeof(buf),
                "{\"run\": %d, \"peak_height_mean\": %g"
                ", \"peak_height_sigma\": %g"
                ", \"peak_height_diff\": %g"
                ", \"peak_integral_mean\": %g"
                ", \"peak_integral_sigma\": %g"
                ", \"peak_integral_diff\": %g"
                ", \"count\": %d}",
                run_number,
                m.fit[kHeight].mean, m.fit[kHeight].sigma, m.fit[kHeight].mean - 35.,
                m.fit[kIntegral].mean, m.fit[kIntegral].sigma, m.fit[kIntegral].mean - 250.,
                event_num_module[m.id]);
            return std::string(buf);
        };

        if (!in_json.empty()) {
            // ── Append mode: add this run to each module's array ──────────
            std::ifstream fin(in_json);
            if (!fin) {
                std::cerr << "Cannot open input JSON: " << in_json << "\n";
            } else {
                std::vector<std::string> lines;
                std::string line;
                while (std::getline(fin, line)) lines.push_back(line);
                fin.close();

                // module lines look like `  "W12": [{...}, ...{...}],`
                for (auto &l : lines) {
                    char key = 0;
                    int number = 0;
                    const auto pos = l.rfind("}]");
                    if (pos == std::string::npos
                        || std::sscanf(l.c_str(), " \"%c%d\"", &key, &number) != 2)
                        continue;
                    for (auto *g : groups) {
                        if (g->key != key) continue;
                        for (const auto &m : g->modules)
                            if (m.number == number) l.insert(pos + 1, ", " + make_entry(m));
                    }
                }

                std::ofstream fout(in_json);
                for (auto &l : lines) fout << l << "\n";
                fout.close();
                std::cerr << "Appended run " << run_number << " to " << in_json << "\n";
            }
        } else {
            // ── Create mode: write new JSON ───────────────────────────────
            std::string out_path = Form("cosmic_modules_run%d.json", run_number);
            std::ofstream json_out(out_path);
            json_out << "{\n";
            for (auto *g : groups)
                for (size_t i = 0; i < g->modules.size(); i++) {
                    json_out << "  \"" << g->key << g->modules[i].number << "\": ["
                             << make_entry(g->modules[i]) << "]";
                    if (i + 1 < g->modules.size()) json_out << ",";
                    json_out << "\n";
                }
            json_out << "}\n";
            json_out.close();
            std::cerr << "JSON written to " << out_path << "\n";
        }
    }

    if (!lg_only) {
        std::ofstream lms_json_out(Form("lms_run_%d.json", run_number));
        lms_json_out << "{\n";
        for (auto *g : groups)
            for (size_t i = 0; i < g->modules.size(); i++) {
                const auto &f = g->modules[i].fit[kLmsIntegral];
                lms_json_out << "  \"" << g->key << g->modules[i].number << "\": {\"run\": " << run_number
                             << ", \"lms_peak\": " << f.mean << ", \"lms_rms\": " << f.sigma << "}";
                if (i + 1 < g->modules.size()) lms_json_out << ",";
                lms_json_out << "\n";
            }
        lms_json_out << "}\n";
        lms_json_out.close();
    }

    outfile.Close();

    return 0;
}
