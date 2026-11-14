
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
#include <TF1.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TLine.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TString.h>
#include <TSystem.h>
#include <TChain.h>
#include <TMarker.h>
#include <TLegend.h>
#include <TNamed.h>

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
long long process_event(TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal,
    std::map<int, TH1F*> &energy_hists, TH2F *energy_vs_theta,
    float Ebeam, const RunConfig &run_config, int max_events = -1,
    const std::string &label = "", std::mutex *io_mtx = nullptr);

// One beam-energy data set (-a/-b/-c): its inputs, per-module energy
// histograms and output directory.
struct BeamSet {
    const char *tag;    // histogram/directory suffix
    const char *label;  // beam energy in GeV, for titles and log lines
    float Ebeam;        // MeV
    float min_ee_sep;   // fit the e-e peak only if it is this far (MeV) from the e-p peak
    std::vector<std::string> inputs;
    RunConfig run_config;
    std::map<int, TH1F*> hists;
    TH2F *energy_vs_theta = nullptr;
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
    std::string seed_calib_file;
    int base_beam_index = 1;

    int max_events = -1;
    int num_threads = kNBeams;
    auto print_usage = [&]() {
        std::cerr << "Usage: nonlinearity -a <recon.root|dir> -b <recon.root|dir> "
                     "-c <recon.root|dir> [more files/dirs...] "
                     "[-o output.root] [-n max_events] [-p png_dir] "
                     "[-s seed_calibration.json] [-B 3.5|2.2|0.7] [-j cores]\n"
                     "  -a  3.5 GeV recon ROOT files or directories (repeatable)\n"
                     "  -b  2.2 GeV recon ROOT files or directories (repeatable)\n"
                     "  -c  0.7 GeV recon ROOT files or directories (repeatable)\n"
                     "  -s  seed calibration JSON (default: built-in calibration)\n"
                     "  -B  base beam energy for non-linearity fit (default: 2.2 GeV)\n"
                     "  -j  number of worker threads (default: 3; may exceed beam count)\n";
    };
    if (argc == 1) {
        print_usage();
        return 0;
    }

    auto require_value = [&](int &i, const std::string &option) -> bool {
        if (i + 1 < argc && argv[i + 1][0] != '-') return true;
        std::cerr << "Option " << option << " requires a value.\n";
        print_usage();
        return false;
    };
    {
        std::vector<std::string> *cur = nullptr;
        for (int i = 1; i < argc; ++i) {
            std::string arg = argv[i];
            if (arg == "-h")       { print_usage(); return 0; }
            else if (arg == "-a")  { cur = &beams[0].inputs; }
            else if (arg == "-b")  { cur = &beams[1].inputs; }
            else if (arg == "-c")  { cur = &beams[2].inputs; }
            else if (arg == "-o")  { cur = nullptr; if (!require_value(i, arg)) return 1; output = argv[++i]; }
            else if (arg == "-n")  { cur = nullptr; if (!require_value(i, arg)) return 1; max_events = std::atoi(argv[++i]); }
            else if (arg == "-p")  { cur = nullptr; if (!require_value(i, arg)) return 1; pngDir = argv[++i]; }
            else if (arg == "-j")  { cur = nullptr; if (!require_value(i, arg)) return 1; num_threads = std::atoi(argv[++i]); }
            else if (arg == "-s")  {
                cur = nullptr;
                if (!require_value(i, arg)) return 1;
                seed_calib_file = argv[++i];
            }
            else if (arg == "-B") {
                cur = nullptr;
                if (!require_value(i, arg)) return 1;
                const std::string base = argv[++i];
                if (base == "3.5" || base == "3p5") base_beam_index = 0;
                else if (base == "2.2" || base == "2p2") base_beam_index = 1;
                else if (base == "0.7" || base == "0p7") base_beam_index = 2;
                else {
                    std::cerr << "Invalid base beam energy '" << base
                              << "'. Choose 3.5, 2.2, or 0.7 GeV.\n";
                    print_usage();
                    return 1;
                }
            }
            else if (!arg.empty() && arg[0] != '-' && cur) {
                auto inputs = ExpandInputPath(arg, IsReconRootName);
                cur->insert(cur->end(), inputs.begin(), inputs.end());
            } else {
                std::cerr << "Unexpected argument: " << arg << "\n";
                print_usage();
                return 1;
            }
        }
    }

    for (const auto &beam : beams) {
        if (beam.inputs.empty()) {
            std::cerr << "No input files specified for " << beam.label << " GeV.\n";
            print_usage();
            return 1;
        }
    }
    if (num_threads <= 0) {
        std::cerr << "Thread count must be greater than zero.\n";
        return 1;
    }
    std::string dbDir = prad2::database_dir();

    for (auto &beam : beams) {
        if (beam.inputs.empty()) continue;
        const int run_num = get_run_int(beam.inputs.front());
        beam.run_config = LoadRunConfig(dbDir + "/runinfo/general.json", run_num);
    }

    fdec::HyCalSystem hycal;
    hycal.Init(dbDir + "/hycal_map.json");
    PhysicsTools physics(hycal);

    // Energy histogram for each crystal
    for (int i = 0; i < hycal.module_count(); ++i) {
        const auto &m = hycal.module(i);
        if (!m.is_pwo4()) continue;
        if (!InHyCalRing(m.x, m.y, 2.0, 16.)) continue;
        std::string hname = "h_energy_" + m.name;
        std::string htitle = "Energy " + m.name + ";E (MeV);Counts";
        for (auto &beam : beams) {
            TH1F *h = new TH1F((hname + "_" + beam.tag).c_str(),
                               (htitle + " (" + beam.label + ")").c_str(), 400, 0, 4000);
            h->SetDirectory(nullptr);
            beam.hists[m.id] = h;
        }
    }
    for (auto &beam : beams) {
        beam.energy_vs_theta = new TH2F(
            Form("h2_energy_vs_theta_%s", beam.tag),
            Form("Cluster Energy vs Scattering Angle (%s GeV);Scattering angle (deg);Cluster energy (MeV)",
                 beam.label),
            45, 0., 4.5, 4000, 0., 4000.);
        beam.energy_vs_theta->SetDirectory(nullptr);
    }

    struct InputJob {
        size_t beam_index;
        std::string path;
        Long64_t max_events;
    };
    std::vector<InputJob> jobs;
    for (size_t b = 0; b < beams.size(); ++b) {
        const auto &inputs = beams[b].inputs;
        if (max_events > 0) {
            const auto limits = DistributeEventBudget(inputs, "recon", max_events);
            for (size_t i = 0; i < inputs.size(); ++i)
                if (limits[i] > 0) jobs.push_back({b, inputs[i], limits[i]});
        } else {
            for (const auto &path : inputs) jobs.push_back({b, path, -1});
        }
    }
    if (jobs.empty()) {
        std::cerr << "No input ROOT files contain events to process.\n";
        return 1;
    }

    struct WorkerData {
        std::array<std::map<int, TH1F *>, kNBeams> hists;
        std::array<TH2F *, kNBeams> energy_vs_theta{};
        std::array<long long, kNBeams> n_sum{};
    };
    const int worker_count = std::min(num_threads, static_cast<int>(jobs.size()));
    std::vector<WorkerData> workers(static_cast<size_t>(worker_count));
    for (int w = 0; w < worker_count; ++w) {
        for (size_t b = 0; b < beams.size(); ++b) {
            for (const auto &[id, source] : beams[b].hists) {
                auto *hist = dynamic_cast<TH1F *>(
                    source->Clone(Form("%s_worker%d", source->GetName(), w)));
                if (!hist) {
                    std::cerr << "Failed to create worker histogram for "
                              << source->GetName() << "\n";
                    return 1;
                }
                hist->SetDirectory(nullptr);
                hist->Reset();
                workers[w].hists[b].emplace(id, hist);
            }
            auto *source = beams[b].energy_vs_theta;
            workers[w].energy_vs_theta[b] = dynamic_cast<TH2F *>(
                source->Clone(Form("%s_worker%d", source->GetName(), w)));
            if (!workers[w].energy_vs_theta[b]) {
                std::cerr << "Failed to create worker histogram for "
                          << source->GetName() << "\n";
                return 1;
            }
            workers[w].energy_vs_theta[b]->SetDirectory(nullptr);
            workers[w].energy_vs_theta[b]->Reset();
        }
    }

    std::mutex io_mtx;
    std::cerr << "Processing " << jobs.size() << " input file(s) with "
              << worker_count << " worker thread(s)\n";
    ParallelFor(jobs.size(), worker_count, [&](size_t job_index, int worker_index) {
        const auto &job = jobs[job_index];
        auto &beam = beams[job.beam_index];
        auto &worker = workers[static_cast<size_t>(worker_index)];
        const std::string label = std::string(beam.label) + " GeV: " + job.path;

        TChain chain("recon");
        if (chain.Add(job.path.c_str()) == 0) {
            std::lock_guard<std::mutex> lk(io_mtx);
            std::cerr << "Failed to add input ROOT file: " << job.path << "\n";
            return;
        }
        EventVars_Recon ev;
        prad2::SetReconReadBranches(&chain, ev);
        const Long64_t entries = chain.GetEntries();
        {
            std::lock_guard<std::mutex> lk(io_mtx);
            std::cerr << "[" << label << "] Processing " << entries << " event(s)\n";
        }
        worker.n_sum[job.beam_index] += process_event(
            &chain, ev, hycal, worker.hists[job.beam_index],
            worker.energy_vs_theta[job.beam_index], beam.Ebeam, beam.run_config,
            static_cast<int>(job.max_events), label, &io_mtx);
    });

    for (size_t b = 0; b < beams.size(); ++b) {
        for (auto &worker : workers) {
            beams[b].n_sum += worker.n_sum[b];
            for (const auto &[id, source] : worker.hists[b])
                beams[b].hists.at(id)->Add(source);
            beams[b].energy_vs_theta->Add(worker.energy_vs_theta[b]);
        }
        std::cerr << "[" << beams[b].label << " GeV] Selected "
                  << beams[b].n_sum << " sum-trigger event(s)\n";
    }

    // a file with no sum-trigger events would silently produce an all-default
    // calibration (every module skipped); refuse to write output in that case
    if (std::any_of(beams.begin(), beams.end(), [](const BeamSet &b) { return b.n_sum <= 0; })) {
        std::cerr << "No sum-trigger events selected (";
        for (int b = 0; b < kNBeams; ++b)
            std::cerr << (b ? ", " : "") << beams[b].label << " GeV: " << beams[b].n_sum;
        std::cerr << "); check trigger_bits in the inputs. Aborting before writing calibration.\n";
        return 1;
    }

    std::string calib_file = seed_calib_file.empty()
        ? dbDir + "/calibration/calibration_factor_3p5_June1.json"
        : seed_calib_file;
    int nmatched = hycal.LoadCalibration(calib_file);
    if (nmatched <= 0) {
        std::cerr << "Failed to load seed calibration from " << calib_file << "\n";
        return 1;
    }
    std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";

    TH1F *h_energy_peak_3p5 = new TH1F("h_energy_peak_3p5", "Energy Peak Distribution;Energy (MeV);Counts", 4000, 0, 4000);

    // calculate non-linearity module by module and save to output file
    gSystem->mkdir(pngDir.c_str(), true);
    TFile outFile(output.empty() ? "nonlinearity_results.root" : output.c_str(), "RECREATE");
    TNamed base_beam_energy("base_beam_energy_GeV", beams[base_beam_index].label);
    base_beam_energy.Write();
    for (auto &beam : beams) beam.dir = outFile.mkdir(Form("energy_%sGeV", beam.tag));
    TDirectory *dir_lin = outFile.mkdir("linearity");
    for (int i = 0; i < hycal.module_count(); i++) {
        const auto &mod = hycal.module(i);
        int mod_id = mod.id;
        if (!mod.is_pwo4()) continue;

        float x = mod.x, y = mod.y;
        if (!InHyCalRing(x, y, 2.0, 16.)) continue;

        // Fit and draw a peak using the same adaptive Gaussian fit as physics_calib.
        auto fitPeakAndDraw = [&](TH1F *h, double Eexp, int color,
                                  double &mean_error) -> float {
            double ypad_min = 0;
            double ypad = h->GetMaximum() * 1.2;
            // Candidate peak search range used by PhysicsTools::fitPeak.
            TLine *wl = new TLine(0.8 * Eexp, ypad_min, 0.8 * Eexp, ypad);
            wl->SetLineColor(6); wl->SetLineStyle(7); wl->SetLineWidth(2); wl->Draw();
            TLine *wr = new TLine(1.2 * Eexp, ypad_min, 1.2 * Eexp, ypad);
            wr->SetLineColor(6); wr->SetLineStyle(7); wr->SetLineWidth(2); wr->Draw();
            auto fit = physics.fitPeak(h, static_cast<float>(Eexp), true);
            const double mean = fit[0];
            if (mean <= 0.) return 0.f;
            mean_error = fit[3] > 0. ? fit[3] : 1.0;

            TF1 *g = h->GetFunction("_fg_");
            if (g) {
                g->SetLineColor(color);
                g->SetLineWidth(3);
                g->Draw("same");
            }
            // draw a marker at the fitted mean, at Gaussian peak height
            TMarker *mk = new TMarker(mean, g ? g->GetParameter(0) : h->GetMaximum(), 29); // star symbol
            mk->SetMarkerColor(color); mk->SetMarkerSize(2.5); mk->Draw();
            return static_cast<float>(mean);
        };

        // --- one pad per beam energy (3.5, 2.2, 0.7 GeV from top), save PNG ---
        TCanvas *ch = new TCanvas(Form("ch_mod_W%d", mod_id-1000),
            Form("Module W%d Histograms", mod_id-1000), 800, 1200);
        ch->Divide(1, kNBeams, 0, 0);

        float exp_ep[kNBeams], exp_ee[kNBeams], peak_ep[kNBeams] = {}, peak_ee[kNBeams] = {};
        double peak_ep_err[kNBeams] = {}, peak_ee_err[kNBeams] = {};
        for (int b = 0; b < kNBeams; ++b) {
            const BeamSet &beam = beams[b];
            const bool top = (b == 0), bottom = (b == kNBeams - 1);
            TH1F *h = beam.hists.at(mod_id);
            float theta = std::atan2(std::sqrt(x*x + y*y), beam.run_config.hycal_z) * 180.f / M_PI;
            exp_ep[b] = physics.ExpectedEnergy(theta, beam.Ebeam, "ep");
            exp_ee[b] = physics.ExpectedEnergy(theta, beam.Ebeam, "ee");

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
            peak_ep[b] = fitPeakAndDraw(h, exp_ep[b], kRed, peak_ep_err[b]);
            if (std::abs(exp_ep[b] - exp_ee[b]) > beam.min_ee_sep)
                peak_ee[b] = fitPeakAndDraw(h, exp_ee[b], kBlue, peak_ee_err[b]);

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

        // The selected base-energy e-p peak anchors the non-linearity fit.
        if (peak_ep[base_beam_index] == 0.f) continue;
        const float E_base = exp_ep[base_beam_index];

        h_energy_peak_3p5->Fill(peak_ep[0]);

        struct FitPoint {
            double E_exp;
            double E_rec;
            double E_rec_error;
            double weight;
        };
        // Give every peak at the selected base beam twice the statistical weight.
        std::vector<FitPoint> points;
        for (int b = 0; b < kNBeams; ++b) {
            const double weight = (b == base_beam_index) ? 2.0 : 1.0;
            if (peak_ep[b] != 0.f)
                points.push_back({exp_ep[b], peak_ep[b], peak_ep_err[b], weight});
            if (peak_ee[b] != 0.f)
                points.push_back({exp_ee[b], peak_ee[b], peak_ee_err[b], weight});
        }

        // make a canvas, E_rec/E_exp vs E_rec
        TCanvas *c = new TCanvas(Form("c_mod_W%d", mod_id-1000), Form("Module W%d Non-linearity", mod_id-1000), 1400, 800);
        c->SetGrid();
        TGraphErrors *g = new TGraphErrors();
        for (const auto &point : points) {
            const int index = g->GetN();
            const double ratio_error = point.E_rec_error / point.E_exp
                                     / std::sqrt(point.weight);
            g->SetPoint(index, point.E_rec, point.E_rec / point.E_exp);
            g->SetPointError(index, 0.0, ratio_error);
        }
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
        constexpr double kNlFitEnergyMax = 3500.0;
        constexpr double kNlFitMaxDeviation = 0.2;
        const double max_abs_t = std::max(
            std::abs(static_cast<double>(E_base)),
            std::abs(kNlFitEnergyMax - static_cast<double>(E_base))) / 1000.0;
        const double nl1_limit = kNlFitMaxDeviation / max_abs_t;
        TF1 *fitLine = new TF1("fitLine",
            [](double *x, double *p){ return 1.0 + p[0] * (x[0] - p[1])/1000.0; },
            xmin, xmax, 2);
        fitLine->SetParameter(0, 0.01);
        fitLine->SetParLimits(0, -nl1_limit, nl1_limit);
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
        const double nl2_nl1_limit = (kNlFitMaxDeviation / 2.0) / max_abs_t;
        const double nl2_nl2_limit = (kNlFitMaxDeviation / 2.0) / (max_abs_t * max_abs_t);
        fitLine2->SetParameter(0, std::clamp(nl, -nl2_nl1_limit, nl2_nl1_limit));
        fitLine2->SetParameter(1, 0.0);
        fitLine2->SetParLimits(0, -nl2_nl1_limit, nl2_nl1_limit);
        fitLine2->SetParLimits(1, -nl2_nl2_limit, nl2_nl2_limit);
        fitLine2->FixParameter(2, E_base);
        fitLine2->SetLineColor(kMagenta+1);
        fitLine2->SetLineWidth(2);
        fitLine2->SetLineStyle(7);
        const int fit_status_2 = g->Fit(fitLine2, "RQ0+");
        fitLine2->Draw("same");

        double nl2_1     = (fit_status_2 == 0 && points.size() >= 3) ? fitLine2->GetParameter(0) : 0.;
        double nl2_1_err = (fit_status_2 == 0 && points.size() >= 3) ? fitLine2->GetParError(0) : 0.;
        double nl2_2     = (fit_status_2 == 0 && points.size() >= 3) ? fitLine2->GetParameter(1) : 0.;
        double nl2_2_err = (fit_status_2 == 0 && points.size() >= 3) ? fitLine2->GetParError(1) : 0.;
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
        for (const auto &point : points) {
            const double Eexp = point.E_exp;
            const double peak = point.E_rec;
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
        for (const auto &point : points) {
            const double Eexp = point.E_exp;
            const double peak = point.E_rec;
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
        const bool force_zero = fdec::test_bit(mod.flag, fdec::kTransition)
                             || fdec::test_bit(mod.flag, fdec::kInnerBound);
        if (force_zero) {
            nl2_1 = 0.0;
            nl2_2 = 0.0;
        }
        hycal.SetCalibBaseEnergy(mod_id, E_base);
        hycal.SetCalibNonLinearity(mod_id, nl2_1, nl2_2);
    }
    for (auto &beam : beams) {
        beam.dir->cd();
        TNamed hycal_z("hycal_z_mm", Form("%.9g", beam.run_config.hycal_z));
        hycal_z.Write();
        for (auto &[id, h] : beam.hists) if (h) h->Write();
        if (beam.energy_vs_theta) beam.energy_vs_theta->Write();
    }
    outFile.cd();
    h_energy_peak_3p5->Write();
    std::cout << "Results saved to " << outFile.GetName() << "\n";
    outFile.Close();

    hycal.PrintCalibConstants("new_calibration_NonLinearity.json");

}

long long process_event(TTree *tree, const EventVars_Recon &ev, const fdec::HyCalSystem &hycal,
    std::map<int, TH1F*> &energy_hists, TH2F *energy_vs_theta,
    float Ebeam, const RunConfig &run_config, int max_events,
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
        if( i % 10000 == 0) {
            log_msg("[" + label + "] Processing event " + std::to_string(i)
                    + "/" + std::to_string(tree->GetEntries()) + "\r", true);
        }
        if ((ev.trigger_bits & prad2::TBIT_sum) == 0) continue;
        if (ev.n_clusters > 2) continue;
        if (ev.n_clusters == 2) {
            const MollerEvent moller(
                {ev.cl_x[0], ev.cl_y[0], ev.cl_z[0], ev.cl_energy[0]},
                {ev.cl_x[1], ev.cl_y[1], ev.cl_z[1], ev.cl_energy[1]});
            const float theta1 = PhysicsTools::GetThetaAngle(
                moller.first.x, moller.first.y, moller.first.z);
            const float theta2 = PhysicsTools::GetThetaAngle(
                moller.second.x, moller.second.y, moller.second.z);
            const bool is_moller =
                std::isfinite(theta1) && std::isfinite(theta2)
                && PhysicsTools::isMoller_kinematic(
                    theta1, moller.first.E, theta2, moller.second.E, Ebeam, 0.03f * 2.0f)
                && PhysicsTools::isBackToBack(moller, 8.f);
            if (!is_moller) continue;
        }
        n_accepted++;

        for( int j = 0; j < ev.n_clusters; j++) {
            int mod_id = ev.cl_center[j];
            if (ev.cl_nblocks[j] < 3) continue;
            auto mod = hycal.module_by_id(mod_id);
            if ( !mod || !mod->is_pwo4()) continue;
            if (!InHyCalRing(mod->x, mod->y, 2.0, 16.)) continue;

            HCHit hc_hit;
            hc_hit.x = ev.cl_x[j];
            hc_hit.y = ev.cl_y[j];
            hc_hit.z = ev.cl_z[j];

            // get the energy ratio of measured energy to expected energy calculated with beam position correction
            // then multiply by the expected energy calculated with module center position correction
            // to produce a corrected energy value that is independent of the beam position
            float theta = std::atan2(std::sqrt(hc_hit.x*hc_hit.x + hc_hit.y*hc_hit.y), run_config.hycal_z) * 180.f / M_PI;
            float expect_energy = PhysicsTools::ExpectedEnergy(theta, Ebeam, "ee");
            float theta_moduleCenter = std::atan2(std::sqrt(mod->x*mod->x + mod->y*mod->y), run_config.hycal_z) * 180.f / M_PI;
            float expect_energy_moduleCenter = PhysicsTools::ExpectedEnergy(theta_moduleCenter, Ebeam, "ee");
            float energy = ev.cl_energy[j] / expect_energy * expect_energy_moduleCenter;

            // Move back to HyCal center coordinates
            hc_hit.x -= ev.hycal_center_x;
            hc_hit.y -= ev.hycal_center_y;

            // require the hit near the seed module centre (|xd|,|yd| < 0.1 module sizes)
            const auto [xd, yd] = mod->cell_offset<float>(hc_hit.x, hc_hit.y);
            if (InHyCalRing(mod->x, mod->y, 2.0, 12.) && (std::abs(xd) >= 0.15f || std::abs(yd) >= 0.15f)) continue;
            if ((std::abs(xd) >= 0.3f || std::abs(yd) >= 0.3f)) continue;

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
            if(theta > 1.1) veto = false;

            if(veto && energy > 600. && Ebeam < 1000.f) continue;

            energy_hists[mod_id]->Fill(energy);
            energy_vs_theta->Fill(theta, energy);
        }
    }
    return n_accepted;
}
