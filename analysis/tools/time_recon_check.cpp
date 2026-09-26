// time_recon_check.cpp
// Time-reconstruction inspection tool for HyCal waveforms.
//
// Input:
//   raw replay ROOT files (events tree), usually files matching *_raw.root,
//   produced by replay_rawdata with peak analysis enabled.
//
// Workflow:
//   - Reconstruct HyCal clusters event-by-event.
//   - Select single-cluster Mott-like events.
//   - Collect representative waveform categories (small/big, in-time/out-time,
//     veto-module small/big).
//   - Re-run WaveAnalyzer on saved waveforms and draw per-waveform overlays
//     (fit curve + reconstructed timing marker).
//
// Output:
//   - One ROOT file: <output_prefix>.root
//   - Five PDF files:
//       <output_prefix>_wave_small_inTime.pdf
//       <output_prefix>_wave_small_outTime.pdf
//       <output_prefix>_wave_big.pdf
//       <output_prefix>_wave_veto_big.pdf
//       <output_prefix>_wave_veto_small.pdf

#include "ConfigSetup.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "WaveAnalyzer.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "RunInfoConfig.h"
#include "gain_factor.h"
#include "PipelineBuilder.h"
#include "ToolUtils.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2Poly.h>
#include <TChain.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TMarker.h>

#include <iostream>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <vector>
#include <array>
#include <memory>
#include <optional>
#include <algorithm>
#include <cmath>

using EventVars = prad2::RawEventData;
using namespace analysis;

struct StoredWaveform {
    int nsamples = 0;
    float dt = 0.0f;
    int peak_height = 0;
    std::array<uint16_t, fdec::MAX_SAMPLES> samples{};
};

// Up to 100 saved waveforms of one category, each also kept as histogram
// h1_<suffix><i> (first 100 samples).  dir is the ROOT directory created for it.
struct WaveCategory {
    const char *suffix, *title, *dir;
    std::array<TH1F *, 100> hists{};
    std::array<StoredWaveform, 100> waves{};
    int count = 0;
};

// Log-normal fit parameters of the first peak for one module type.
struct FitParamHists {
    TH1F *mu, *sigma, *chi2_ndf;
    TH2F *mu_vs_sigma, *mu_vs_height, *sigma_vs_height;

    FitParamHists(const char *tag, const char *label)
        : mu(new TH1F(Form("h1_mu_%s", tag), Form("Fitted Mu for %s;Mu [ADC];Counts", label), 100, 0, 2)),
          sigma(new TH1F(Form("h1_sigma_%s", tag), Form("Fitted Sigma for %s;Sigma [ADC];Counts", label), 100, 0, 1)),
          chi2_ndf(new TH1F(Form("h1_chi2_ndf_%s", tag), Form("Fitted Chi2/NDF for %s;Chi2/NDF [ADC];Counts", label), 100, 0, 15)),
          mu_vs_sigma(new TH2F(Form("h2_mu_vs_sigma_%s", tag), Form("Mu vs Sigma for %s;Mu [ADC];Sigma [ADC]", label), 100, 0, 2, 100, 0, 1)),
          mu_vs_height(new TH2F(Form("h2_mu_vs_height_%s", tag), Form("Mu vs Height for %s;Mu [ADC];Height [ADC]", label), 100, 0, 2, 1000, 0, 3000)),
          sigma_vs_height(new TH2F(Form("h2_sigma_vs_height_%s", tag), Form("Sigma vs Height for %s;Sigma [ADC];Height [ADC]", label), 100, 0, 1, 1000, 0, 3000))
    {}

    void Fill(const fdec::LogNormalFitResult &fit, float height) const
    {
        mu->Fill(fit.mu);
        sigma->Fill(fit.sigma);
        chi2_ndf->Fill(fit.chi2_per_dof/fit.A);
        mu_vs_sigma->Fill(fit.mu, fit.sigma);
        mu_vs_height->Fill(fit.mu, height);
        sigma_vs_height->Fill(fit.sigma, height);
    }

    void Write() const
    {
        mu->Write();
        sigma->Write();
        chi2_ndf->Write();
        mu_vs_sigma->Write();
        mu_vs_height->Write();
        sigma_vs_height->Write();
    }
};

int main(int argc, char *argv[])
{
    std::string db_dir = prad2::database_dir();

    // ── Argument parsing ─────────────────────────────────────────────────────
    std::string output_path_name, daq_config_file, recon_config_file, gem_ped_file;
    int  max_events  = -1;
    int  num_threads = 4;
    int  num_files   = -1;

    int opt;
    while ((opt = getopt(argc, argv, "o:n:f:j:")) != -1) {
        switch (opt) {
            case 'o': output_path_name = optarg; break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'f': num_files        = std::atoi(optarg); break;
            case 'j': num_threads     = std::atoi(optarg); break;
        }
    }

    std::vector<std::string> root_files = CollectInputs(argc, argv, optind, IsRawRootName, num_files);
    if (root_files.empty()) {
        std::cerr << "No input files specified.\n";
        std::cerr << "Usage: time_recon_check <input_raw.root|dir> [more...] "
                     "[-o ./output_name(no extension)] [-n max_events] [-f nfiles] [-j threads]\n";
        return 1;
    }

    if (output_path_name.empty()) {
        std::cerr << "No output prefix provided. Please pass -o <output_prefix>.\n";
        return 1;
    }

    TChain tree("events");
    for (const auto &file : root_files) {
        tree.Add(file.c_str());
    }

    auto ev = std::make_unique<EventVars>();
    prad2::SetRawReadBranches(&tree, *ev);

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
    recon_config_file = db_dir + "/reconstruction_config.json";
    daq_config_file = db_dir + "/daq_config.json";

    prad2::Pipeline pipeline = prad2::PipelineBuilder()
        .set_database_dir(db_dir)
        .set_recon_config(recon_config_file)
        .set_daq_config(daq_config_file)
        .set_gem_pedestal(gem_ped_file)     // empty falls back to RunConfig default
        .set_run_number(run_num)
        .set_log_stream(&std::cerr)
        .build();

    const auto &hycal        = pipeline.hycal;
    auto &cluster_cfg        = pipeline.hycal_cluster_cfg;
    const auto &hc_time_cuts = pipeline.hycal_time_cuts;

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);
    cluster_cfg.seed_time_window = 0;

    fdec::WaveAnalyzer ana(pipeline.daq_cfg.wave_cfg);
    fdec::WaveResult wres;

    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(gRunConfig, run_num);

    TH1F *h1_cluster_energy = new TH1F("h1_cluster_energy", "Cluster Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_time_diff_seed = new TH1F("h1_time_diff_seed", "Time Difference with Seed Module;#Delta t [ns];Events", 500, -25, 25);
    TH2F *h2_height_vs_dtime = new TH2F("h2_height_vs_dtime", "Peak Height vs delta Time;delta Time [ns];Height [ADC]", 500, -25, 25, 1000, 0, 3000);
    // Fitting parameters for the waveforms
    FitParamHists fit_crystal("crystal", "Crystal");
    FitParamHists fit_veto("veto", "Veto");

    // TH2Poly map for per-module out-of-time ratio.
    std::vector<int> module_bin_by_index;
    TH2Poly *h2_ratio_out4ns_per_crystal = PhysicsTools::MakeModuleMap(
        hycal, "h2_ratio_out4ns_per_crystal", "Neighbor - Seed Ratio |#Deltat| > 4 ns;X;Y", 350.,
        module_bin_by_index);
    h2_ratio_out4ns_per_crystal->SetContour(255);
    h2_ratio_out4ns_per_crystal->SetStats(false);
    h2_ratio_out4ns_per_crystal->SetMinimum(0);
    h2_ratio_out4ns_per_crystal->SetMaximum(0.15);
    h2_ratio_out4ns_per_crystal->SetOption("colz");
    std::vector<int> module_total_counts(hycal.module_count(), 0);
    std::vector<int> module_out4ns_counts(hycal.module_count(), 0);

    // Save the waveforms for later analysis (in PDF and ROOT write order)
    std::array<WaveCategory, 5> wave_cats{{
        {"wave_small_inTime",  "Small In-Time Waveform",  "waveforms_small_inTime"},
        {"wave_small_outTime", "Small Out-Time Waveform", "waveforms_small_outTime"},
        {"wave_big",           "Big Waveform",            "waveforms_big"},
        {"wave_veto_big",      "Veto Big Waveform",       "waveforms_veto_big"},
        {"wave_veto_small",    "Veto Small Waveform",     "waveforms_veto_small"},
    }};
    auto &[wave_small_inTime, wave_small_outTime, wave_big, wave_veto_big, wave_veto_small] = wave_cats;
    for (auto &cat : wave_cats)
        for (int i = 0; i < 100; ++i)
            cat.hists[i] = new TH1F(Form("h1_%s%d", cat.suffix, i), Form("%s%d;Samples;ADC", cat.title, i), 100, -0.5, 99.5);

    auto store_wave = [&](WaveCategory &cat, int j, float dt) {
        if (cat.count >= 100) return;
        auto &wf = cat.waves[cat.count];
        const int ns = std::min<int>(static_cast<int>(ev->nsamples[j]), fdec::MAX_SAMPLES);
        wf.nsamples = ns;
        wf.dt = dt;
        wf.peak_height = static_cast<int>(ev->peak_height[j][0]);
        for (int s = 0; s < ns; ++s) {
            wf.samples[s] = ev->samples[j][s];
            if (s < 100) {
                cat.hists[cat.count]->SetBinContent(s + 1, ev->samples[j][s]);
            }
        }
        ++cat.count;
    };

    long long nentries = tree.GetEntries();
    for (long long i = 0; i < nentries; ++i) {
        tree.GetEntry(i);
        if (i >= max_events && max_events > 0) break;
        if (i % 10000 == 0) std::cout << "Processed " << i << " / " << nentries << " entries.\r" << std::flush;

        if ((ev->trigger_bits & prad2::TBIT_sum) == 0) continue;
        if (ev->nch > 70) continue; // channel numbers too high, likely not a clean event

        clusterer.Clear();

        // Per-event gain correction (time-series lookup by event number).
        const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(ev->event_num));

        for (int j = 0; j < ev->nch; ++j) {

            if (ev->module_id[j] > prad2::kVetoIdBase && ev->module_id[j] <= prad2::kVetoIdBase + 4
                && ev->npeaks[j] == 1) {
                if (ev->peak_height[j][0] > 10) {
                    ana.Analyze(ev->samples[j], ev->nsamples[j], wres);
                    fit_veto.Fill(wres.peaks_fit[0], ev->peak_height[j][0]);
                }
                if (ev->peak_height[j][0] < 50) store_wave(wave_veto_small, j, 0.0f);
                if (ev->peak_height[j][0] >= 100) store_wave(wave_veto_big, j, 0.0f);
            }

            const auto *mod = hycal.module_by_id(ev->module_id[j]);
            if (!mod) continue;
            if (!mod->is_pwo4()) continue;

            const float gain = gain_corr.ModuleGain(mod->id);
            float time_offset = mod->time_offset;

            const auto hc_win = hc_time_cuts.at(mod->index);
            // Multi-pulse mode: push every peak inside the trigger
            // window into the clusterer; the seed-anchored timing
            // coincidence cut is applied inside HyCalCluster.
            for (int p = 0; p < ev->npeaks[j]; ++p) {
                float peak_time = ev->peak_time[j][p] - time_offset;
                if (peak_time <= hc_win.lo) continue;
                if (peak_time >= hc_win.hi) continue;
                float adc = ev->peak_integral[j][p] * gain;
                float energy = static_cast<float>(mod->energize(adc));
                clusterer.AddHit(mod->index, energy, peak_time);
            }
        }
        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);

        // select single cluster Mott events
        if (hits.size() != 1 || hits[0].nblocks < 3) continue;
        float hc_x = hits[0].x, hc_y = hits[0].y, hc_z = gRunConfig.hycal_z;
        float theta = PhysicsTools::GetThetaAngle(hc_x, hc_y, hc_z);
        if (std::abs(hits[0].energy - gRunConfig.Ebeam) > 3.f * hycal.EnergyResolution(gRunConfig.Ebeam)) continue;
        if (theta < 0.7 ) continue; 

        int seed_id = hits[0].center_id;
        const auto *seed_mod = hycal.module_by_id(seed_id);
        if (!seed_mod || !seed_mod->is_pwo4()) continue;

        //the time window cut is already applied in the clusterer
        float dt_max = cluster_cfg.seed_time_window;
        if (dt_max <= 0 ) dt_max = 100.0f;

        bool seed_not_clean = false;

        for (int j = 0; j < ev->nch; ++j) {
            const auto *mod = hycal.module_by_id(ev->module_id[j]);
            if (!mod || !mod->is_pwo4()) continue;

            float time_offset = mod->time_offset;

            // require only one peak in the seed module
            if (mod->id == seed_id && ev->npeaks[j] != 1) {
                seed_not_clean = true;
                break;
            }
            // select the 1st layer of neighboring modules
            if (std::abs(mod->x - seed_mod->x) < mod->size_x * 1.5f &&
                std::abs(mod->y - seed_mod->y) < mod->size_y * 1.5f && mod->id != seed_mod->id
                && ev->npeaks[j] == 1 && ev->module_id[j] > 0) 
            {
                float peak_time = ev->peak_time[j][0] - time_offset;
                float dt = peak_time - hits[0].time;
                const int mod_index = mod->index;
                if (mod_index >= 0 && mod_index < static_cast<int>(module_total_counts.size())) {
                    module_total_counts[mod_index]++;
                    if (dt > 4.0f) {
                        module_out4ns_counts[mod_index]++;
                    }
                }
                h1_time_diff_seed->Fill(dt);
                h2_height_vs_dtime->Fill(dt, ev->peak_height[j][0]);
                if (dt > 4.0 && ev->peak_height[j][0] < 30) store_wave(wave_small_outTime, j, dt);
                if (std::abs(dt) < dt_max && ev->peak_height[j][0] < 30) store_wave(wave_small_inTime, j, dt);
                if (ev->peak_height[j][0] > 10) {
                    ana.Analyze(ev->samples[j], ev->nsamples[j], wres);
                    fit_crystal.Fill(wres.peaks_fit[0], ev->peak_height[j][0]);
                    store_wave(wave_big, j, dt);
                }
            }
        }
        if (!seed_not_clean) {
            h1_cluster_energy->Fill(hits[0].energy);
        }
    }

    for (int m = 0; m < hycal.module_count(); ++m) {
        const int bin_id = module_bin_by_index[m];
        if (bin_id < 0) continue;

        const int total_count = module_total_counts[m];
        if (total_count <= 0) continue;

        h2_ratio_out4ns_per_crystal->SetBinContent(
            bin_id, static_cast<double>(module_out4ns_counts[m]) / total_count);
    }

    // Draw waveform pages to PDFs: one waveform per page with WaveAnalyzer re-fit overlay.
    auto draw_waveforms_pdf = [&](const WaveCategory &cat) {
        const std::string pdf_name = output_path_name + "_" + cat.suffix + ".pdf";
        TCanvas c_wave(Form("c_wave_%s", cat.suffix), cat.title, 900, 650);
        c_wave.Print((pdf_name + "[").c_str());

        for (int i = 0; i < cat.count; ++i) {
            auto *h = cat.hists[i];
            const auto &wf = cat.waves[i];
            if (!h || wf.nsamples <= 0) continue;

            std::optional<TF1> fit_curve;
            std::optional<TMarker> reco_point;

            h->SetTitle(Form("%s %d;Samples;ADC", cat.title, i));
            h->SetStats(0);
            h->SetLineColor(kBlue + 1);
            h->SetLineWidth(1);
            h->SetLineStyle(2);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(0.6);
            h->SetMarkerColor(kBlue + 1);
            h->GetXaxis()->SetRangeUser(20.0, 70.0);
            h->Draw("hist l");
            h->Draw("P same");

            float smoothed[fdec::MAX_SAMPLES];
            ana.smooth(wf.samples.data(), wf.nsamples, smoothed);
            // draw the smoothed waveform on top of the original histogram
            std::vector<double> x_vals(wf.nsamples);
            std::vector<double> y_vals(wf.nsamples);
            for (int s = 0; s < wf.nsamples; ++s) {
                x_vals[s] = static_cast<double>(s);
                y_vals[s] = static_cast<double>(smoothed[s]);
            }
            TGraph g_smoothed(wf.nsamples, x_vals.data(), y_vals.data());
            g_smoothed.SetLineColor(kGreen + 1);
            g_smoothed.SetLineWidth(2);
            g_smoothed.Draw("L same");

            fdec::WaveResult wave_res;
            ana.Analyze(wf.samples.data(), wf.nsamples, wave_res, 0.0f);

            int best_peak = -1;
            float best_h = -1.0f;
            for (int k = 0; k < wave_res.npeaks; ++k) {
                if (wave_res.peaks[k].height > best_h) {
                    best_h = wave_res.peaks[k].height;
                    best_peak = k;
                }
            }

            bool fit_ok = false;
            fdec::LogNormalFitResult ln_fit{};
            int fit_left = 0;
            int fit_right = 0;
            if (best_peak >= 0) {
                ln_fit = wave_res.peaks_fit[best_peak];
                fit_ok = ln_fit.ok;
                fit_left = std::max(0, wave_res.peaks[best_peak].left - 4);
                fit_right = std::min(wf.nsamples - 1, wave_res.peaks[best_peak].pos + 4);
            }

            const float ped_mean = wave_res.ped.mean;
            const auto ln_value = [ped_mean, ln_fit](float sample) {
                return fdec::WaveAnalyzer::log_normal_pulse_value(
                    sample, ped_mean, ln_fit.A, ln_fit.t0, ln_fit.mu, ln_fit.sigma);
            };

            if (fit_ok) {
                fit_curve.emplace(Form("f_wave_fit_%s_%d", cat.suffix, i),
                    [=](double *x, double *) {
                        return static_cast<double>(ln_value(static_cast<float>(x[0])));
                    },
                    static_cast<double>(fit_left), static_cast<double>(fit_right), 0);
                fit_curve->SetLineColor(kRed + 1);
                fit_curve->SetLineWidth(3);
                fit_curve->SetNpx(500);
                fit_curve->Draw("same");
            }

            const float clk_ns = ana.cfg.clk_ns();
            float reco_time_ns = -9999.0f;
            float reco_sample = -1.0f;
            const char *reco_algo = "unknown";
            if (best_peak >= 0) {
                reco_time_ns = wave_res.peaks[best_peak].time;
                reco_sample = reco_time_ns / clk_ns;
                switch (wave_res.peaks[best_peak].time_algo) {
                    case fdec::T_PICKOFF_FIT_CFD:
                        reco_algo = "fit cfd";
                        break;
                    case fdec::T_PICKOFF_LINEAR_CFD:
                        reco_algo = "linear cfd";
                        break;
                    case fdec::T_PICKOFF_PEAKING_SUBSAMPLE:
                        reco_algo = "peaking subsample";
                        break;
                    default:
                        reco_algo = "unknown";
                        break;
                }
            }
            if (reco_sample >= 0.0f && reco_sample <= static_cast<float>(wf.nsamples)) {
                float reco_y = 0.0f;
                if (fit_ok) {
                    reco_y = ln_value(reco_sample);
                } else {
                    // Fall back to the waveform value if the fit is unavailable.
                    const int i0 = std::max(0, std::min(wf.nsamples - 1, static_cast<int>(std::floor(reco_sample))));
                    const int i1 = std::max(0, std::min(wf.nsamples - 1, i0 + 1));
                    const float y0 = static_cast<float>(wf.samples[i0]);
                    const float y1 = static_cast<float>(wf.samples[i1]);
                    const float frac = std::max(0.0f, std::min(1.0f, reco_sample - static_cast<float>(i0)));
                    reco_y = y0 + (y1 - y0) * frac;
                }

                reco_point.emplace(reco_sample, reco_y, 20);
                reco_point->SetMarkerColor(kGreen + 2);
                reco_point->SetMarkerSize(1.2f);
                reco_point->Draw("same");
            }

            TLatex label;
            label.SetNDC();
            label.SetTextSize(0.030f);
            label.SetTextAlign(31);
            label.DrawLatex(0.88f, 0.92f, Form("dt = %.2f ns, peak height = %d ADC", wf.dt, wf.peak_height));
            label.DrawLatex(0.88f, 0.88f, Form("WaveAnalyzer peaks = %d, fit ok = %d", wave_res.npeaks, fit_ok ? 1 : 0));
            if (best_peak >= 0) {
                label.DrawLatex(0.88f, 0.84f, Form("reco time = %.2f ns (sample %.2f), algo: %s",
                                                   reco_time_ns, reco_sample, reco_algo));
            }
            if (fit_ok) {
                label.DrawLatex(0.88f, 0.80f, Form("logn fit: t0=%.2f, mu=%.2f, sigma=%.2f, chi2/ndf=%.3f",
                                                   ln_fit.t0, ln_fit.mu, ln_fit.sigma, ln_fit.chi2_per_dof));
            }

            c_wave.Print(pdf_name.c_str());
        }
        c_wave.Print((pdf_name + "]").c_str());
    };

    for (const auto &cat : wave_cats) draw_waveforms_pdf(cat);

    TFile *output_file = new TFile((output_path_name + ".root").c_str(), "RECREATE");
    h1_cluster_energy->Write();
    h1_time_diff_seed->Write();
    h2_height_vs_dtime->Write();
    h2_ratio_out4ns_per_crystal->Write();
    fit_crystal.Write();
    fit_veto.Write();
    for (const auto &cat : wave_cats) {
        output_file->cd();
        output_file->mkdir(cat.dir);
        for (int i = 0; i < cat.count; ++i) cat.hists[i]->Write();
    }
    output_file->Close();

}