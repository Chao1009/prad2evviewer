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

#include "Replay.h"
#include "PhysicsTools.h"
#include "MatchingTools.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"
#include "GemSystem.h"
#include "WaveAnalyzer.h"
#include "EventData.h"
#include "EventData_io.h"
#include "InstallPaths.h"
#include "load_daq_config.h"
#include "RunInfoConfig.h"
#include "gain_factor.h"
#include "PulseTemplateStore.h"
#include "PipelineBuilder.h"

#include <TFile.h>
#include <TH1F.h>
#include <TH2Poly.h>
#include <TChain.h>
#include <TFitResult.h>
#include <TFile.h>
#include <TText.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TGraph.h>
#include <TLatex.h>
#include <TLine.h>
#include <TMarker.h>

#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstdlib>
#include <getopt.h>
#include <filesystem>
#include <vector>
#include <array>
#include <memory>
#include <optional>
#include <algorithm>
#include <map>
#include <queue>
#include <set>
#include <utility>
#include <cmath>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompSVD.h>

#include <nlohmann/json.hpp>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

namespace fs = std::filesystem;

using EventVars = prad2::RawEventData;
using namespace analysis;

// ── File collection helper ───────────────────────────────────────────────────
static std::vector<std::string> collectRootFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &entry : fs::directory_iterator(path)) {
            if (entry.is_regular_file() &&
                entry.path().filename().string().find("_raw.root") != std::string::npos)
                files.push_back(entry.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

// ── Main ─────────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    if (const char *env = std::getenv("PRAD2_DATABASE_DIR")) db_dir = env;

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

    // Collect all input files
    std::vector<std::string> root_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectRootFiles(argv[i]);
        if (num_files > 0) {
            int remaining = num_files - static_cast<int>(root_files.size());
            if (remaining <= 0) break;
            int take = std::min(remaining, static_cast<int>(f.size()));
            root_files.insert(root_files.end(), f.begin(), f.begin() + take);
            if (static_cast<int>(root_files.size()) >= num_files) break;
        } else {
            root_files.insert(root_files.end(), f.begin(), f.end());
        }
    }
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

    EventVars ev;
    prad2::SetRawReadBranches(&tree, ev);

    int run_num = get_run_int(root_files.front());
    gRunConfig = LoadRunConfig(db_dir + "/runinfo/general.json", run_num);
    recon_config_file = db_dir + "/reconstruction_config.json";
    daq_config_file = db_dir + "/daq_config.json";

    evc::DaqConfig daq_cfg;

    // Detectors: PRad-II flows through PipelineBuilder so the wiring stays in
    // one place (see prad2det/include/PipelineBuilder.h).
    fdec::HyCalSystem                 hycal;
    gem::GemSystem                    gem_sys;
    fdec::ClusterConfig               cluster_cfg;
    prad2::HyCalTimeCuts              hc_time_cuts;
    prad2::HyCalRfOffsets             hc_rf_offsets;
    DetectorTransform                 hycal_transform;
    std::array<DetectorTransform, 4>  gem_transforms;
    std::unordered_map<int, int>      roc_to_crate;
    int                               match_method = 1;

    prad2::Pipeline pipeline = prad2::PipelineBuilder()
        .set_database_dir(db_dir)
        .set_recon_config(recon_config_file)
        .set_daq_config(daq_config_file)
        .set_gem_pedestal(gem_ped_file)     // empty falls back to RunConfig default
        .set_run_number(run_num)
        .set_log_stream(&std::cerr)
        .build();

    daq_cfg          = std::move(pipeline.daq_cfg);
    hycal            = std::move(pipeline.hycal);
    gem_sys          = std::move(pipeline.gem);
    cluster_cfg      = pipeline.hycal_cluster_cfg;
    hc_time_cuts     = std::move(pipeline.hycal_time_cuts);
    hc_rf_offsets    = std::move(pipeline.hycal_rf_offsets);
    hycal_transform  = pipeline.hycal_transform;
    gem_transforms   = pipeline.gem_transforms;
    match_method     = pipeline.match_method;

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);
    cluster_cfg.seed_time_window = 0;
    gem::GemCluster      gem_clusterer;
    MatchingTools        matching(match_method);

    //initialize tools for cluster reconstruction
    fdec::WaveAnalyzer ana(daq_cfg.wave_cfg);
    fdec::PulseTemplateStore template_store;
    if (daq_cfg.wave_cfg.nnls_deconv.enabled
        && !daq_cfg.wave_cfg.nnls_deconv.template_file.empty()) {
        template_store.LoadFromFile(
            db_dir + "/" + daq_cfg.wave_cfg.nnls_deconv.template_file,
            daq_cfg.wave_cfg);
    }
    ana.SetTemplateStore(&template_store);
    fdec::WaveResult wres;

    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(
        gRunConfig.gain_data_dir + "/gain_correction", run_num);

    // create histograms you want to fill for shower profile analysis
    TH1F *h1_cluster_energy = new TH1F("h1_cluster_energy", "Cluster Energy;Energy [MeV];Counts", 4000, 0, 4000);
    TH1F *h1_time_diff_seed = new TH1F("h1_time_diff_seed", "Time Difference with Seed Module;#Delta t [ns];Events", 500, -25, 25);
    TH2F *h2_height_vs_dtime = new TH2F("h2_height_vs_dtime", "Peak Height vs delta Time;delta Time [ns];Height [ADC]", 500, -25, 25, 1000, 0, 3000);
    // Save the waveforms for later analysis
    TH1F *h1_wave_small_inTime[100], *h1_wave_small_outTime[100];
    TH1F *h1_wave_big[100];
    TH1F *h1_wave_veto_big[100], *h1_wave_veto_small[100];
    // Fitting parameters for the waveforms
    TH1F *h1_mu_crystal = new TH1F("h1_mu_crystal", "Fitted Mu for Crystal;Mu [ADC];Counts", 100, 0, 4);
    TH1F *h1_sigma_crystal = new TH1F("h1_sigma_crystal", "Fitted Sigma for Crystal;Sigma [ADC];Counts", 100, 0, 4);
    TH1F *h1_chi2_ndf_crystal = new TH1F("h1_chi2_ndf_crystal", "Fitted Chi2/NDF for Crystal;Chi2/NDF [ADC];Counts", 100, 0, 50);
    TH1F *h1_mu_veto = new TH1F("h1_mu_veto", "Fitted Mu for Veto;Mu [ADC];Counts", 100, 0, 4);
    TH1F *h1_sigma_veto = new TH1F("h1_sigma_veto", "Fitted Sigma for Veto;Sigma [ADC];Counts", 100, 0, 4);
    TH1F *h1_chi2_ndf_veto = new TH1F("h1_chi2_ndf_veto", "Fitted Chi2/NDF for Veto;Chi2/NDF [ADC];Counts", 100, 0, 50);

    for (int i = 0; i < 100; ++i) {
        h1_wave_small_inTime[i] = new TH1F(Form("h1_wave_small_inTime%d", i), Form("Small In-Time Waveform%d;Samples;ADC", i), 100, -0.5, 99.5);
        h1_wave_small_outTime[i] = new TH1F(Form("h1_wave_small_outTime%d", i), Form("Small Out-Time Waveform%d;Samples;ADC", i), 100, -0.5, 99.5);
        h1_wave_big[i] = new TH1F(Form("h1_wave_big%d", i), Form("Big Waveform%d;Samples;ADC", i), 100, -0.5, 99.5);
        h1_wave_veto_big[i] = new TH1F(Form("h1_wave_veto_big%d", i), Form("Veto Big Waveform%d;Samples;ADC", i), 100, -0.5, 99.5);
        h1_wave_veto_small[i] = new TH1F(Form("h1_wave_veto_small%d", i), Form("Veto Small Waveform%d;Samples;ADC", i), 100, -0.5, 99.5);
    }
    int wave_count_small_inTime = 0;
    int wave_count_small_outTime = 0;
    int wave_count_big = 0;
    int wave_count_veto_big = 0;
    int wave_count_veto_small = 0;

    struct StoredWaveform {
        int nsamples = 0;
        float dt = 0.0f;
        int peak_height = 0;
        std::array<uint16_t, fdec::MAX_SAMPLES> samples{};
    };
    std::array<StoredWaveform, 100> stored_wave_small_inTime;
    std::array<StoredWaveform, 100> stored_wave_small_outTime;
    std::array<StoredWaveform, 100> stored_wave_big;
    std::array<StoredWaveform, 100> stored_wave_veto_big;
    std::array<StoredWaveform, 100> stored_wave_veto_small;

    // Open an ouput PDF file to save the histograms

    // Here loop over the events in the TChain, read channels data, reconstruct clusters, and fill the histograms
    long long nentries = tree.GetEntries();
    for (long long i = 0; i < nentries; ++i) {
        tree.GetEntry(i);
        if (i >= max_events && max_events > 0) break;
        if (i % 10000 == 0) std::cout << "Processed " << i << " / " << nentries << " entries.\r" << std::flush;

        // assume you are selecting Mott-like events with a single cluster in the HyCal
        if ((ev.trigger_bits & prad2::TBIT_sum) == 0) continue;
        if (ev.nch > 70) continue; // channel numbers too high, likely not a clean event

        // Reconstruct clusters for this event.
        clusterer.Clear();

        // Per-event gain correction (time-series lookup by event number).
        const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(ev.event_num));

        for (int j = 0; j < ev.nch; ++j) {

            if (ev.module_id[j] >= 3001 && ev.module_id[j] <= 3004 && ev.npeaks[j] == 1) {
                if (ev.peak_height[j][0] > 100 && ev.npeaks[j] == 1) {
                    ana.Analyze(ev.samples[j], ev.nsamples[j], wres);
                    h1_mu_veto->Fill(wres.peaks_fit[0].mu);
                    h1_sigma_veto->Fill(wres.peaks_fit[0].sigma);
                    h1_chi2_ndf_veto->Fill(wres.peaks_fit[0].chi2_per_dof/wres.peaks_fit[0].A);
                }
                // save the waveforms for veto modules
                if (wave_count_veto_small < 100 && ev.peak_height[j][0] < 50) {
                    auto &wf = stored_wave_veto_small[wave_count_veto_small];
                    const int ns = std::min<int>(static_cast<int>(ev.nsamples[j]), fdec::MAX_SAMPLES);
                    wf.nsamples = ns;
                    wf.dt = 0.0f;
                    wf.peak_height = static_cast<int>(ev.peak_height[j][0]);
                    for (int s = 0; s < ns; ++s) {
                        wf.samples[s] = ev.samples[j][s];
                        if (s < 100) {
                            h1_wave_veto_small[wave_count_veto_small]->SetBinContent(s + 1, ev.samples[j][s]);
                        }
                    }
                    wave_count_veto_small++;
                }
                if (wave_count_veto_big < 100 && ev.peak_height[j][0] >= 100) {
                    auto &wf = stored_wave_veto_big[wave_count_veto_big];
                    const int ns = std::min<int>(static_cast<int>(ev.nsamples[j]), fdec::MAX_SAMPLES);
                    wf.nsamples = ns;
                    wf.dt = 0.0f;
                    wf.peak_height = static_cast<int>(ev.peak_height[j][0]);
                    for (int s = 0; s < ns; ++s) {
                        wf.samples[s] = ev.samples[j][s];
                        if (s < 100) {
                            h1_wave_veto_big[wave_count_veto_big]->SetBinContent(s + 1, ev.samples[j][s]);
                        }
                    }
                    wave_count_veto_big++;
                }
            }

            const auto *mod = hycal.module_by_id(ev.module_id[j]);
            if (!mod) continue;
            if (!mod->is_pwo4()) continue;

            // Per-ID gain correction: average of three LMS channels.
            const float gain = (mod->id > 1000)
                ? (gain_corr.w[mod->id - 1000].corr[1] + gain_corr.w[mod->id - 1000].corr[2]) / 2.0f
                : 1.0f; // default gain factor
            // timing offset for this module
            float time_offset = mod->time_offset;

            const auto hc_win = hc_time_cuts.at(mod->index);
            // Multi-pulse mode: push every peak inside the trigger
            // window into the clusterer; the seed-anchored timing
            // coincidence cut is applied inside HyCalCluster, 
            // maximum 2 ns difference between pulses in a same cluster.
            for (int p = 0; p < ev.npeaks[j]; ++p) {
                float peak_time = ev.peak_time[j][p] - time_offset; // apply module time offset
                if (peak_time <= hc_win.lo) continue;
                if (peak_time >= hc_win.hi) continue;
                float adc = ev.peak_integral[j][p] * gain;
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
        float theta = std::atan2(std::sqrt(hc_x * hc_x + hc_y * hc_y), hc_z) * 180.0 / M_PI;
        if (std::abs(hits[0].energy - gRunConfig.Ebeam) > 3. * 0.033 * std::sqrt(gRunConfig.Ebeam * 1000.)) continue;
        if (theta < 0.7 ) continue; 

        int seed_id = hits[0].center_id;
        const auto *seed_mod = hycal.module_by_id(seed_id);
        if (!seed_mod || !seed_mod->is_pwo4()) continue;

        //the time window cut is already applied in the clusterer
        float dt_max = cluster_cfg.seed_time_window;
        if (dt_max <= 0 ) dt_max = 100.0f;

        // counts for fired modules in the 1st and 2nd layer of neighbors
        int n1st_layer = 0;
        int n2nd_layer = 0;

        bool seed_not_clean = false;

        // loop over the hits in the cluster
        // Here go to the event by event analysis, add more analysis code or selection .etc here
        for (int j = 0; j < ev.nch; ++j) {
            const auto *mod = hycal.module_by_id(ev.module_id[j]);
            if (!mod || !mod->is_pwo4()) continue;

            // Per-ID gain correction: average of three LMS channels.
            const float gain = (mod->id > 1000)
                ? (gain_corr.w[mod->id - 1000].corr[1] + gain_corr.w[mod->id - 1000].corr[2]) / 2.0f
                : 1.0f; // default gain factor
            // timing offset for this module
            float time_offset = mod->time_offset;

            // select the seed module and its energy (require only one peak in seed module)
            if (mod->id == seed_id) {
                if (ev.npeaks[j] != 1) {
                    seed_not_clean = true;
                    break;
                }
                float seed_energy = static_cast<float>(mod->energize(ev.peak_integral[j][0] * gain));
            }
            // select the 1st layer of neighboring modules
            if (std::abs(mod->x - seed_mod->x) < mod->size_x * 1.5f &&
                std::abs(mod->y - seed_mod->y) < mod->size_y * 1.5f && mod->id != seed_mod->id
                && ev.npeaks[j] == 1) 
            {
                float peak_time = ev.peak_time[j][0] - time_offset;
                float dt = peak_time - hits[0].time;
                h1_time_diff_seed->Fill(dt);
                h2_height_vs_dtime->Fill(dt, ev.peak_height[j][0]);
                if (dt > 6.0 && ev.peak_height[j][0] < 30) {
                    if (wave_count_small_outTime < 100) {
                        // Store the waveforms
                        auto &wf = stored_wave_small_outTime[wave_count_small_outTime];
                        const int ns = std::min<int>(static_cast<int>(ev.nsamples[j]), fdec::MAX_SAMPLES);
                        wf.nsamples = ns;
                        wf.dt = dt;
                        wf.peak_height = ev.peak_height[j][0];
                        for (int s = 0; s < ns; ++s) {
                            wf.samples[s] = ev.samples[j][s];
                            if (s < 100) {
                                h1_wave_small_outTime[wave_count_small_outTime]->SetBinContent(s + 1, ev.samples[j][s]);
                            }
                        }
                    }
                    wave_count_small_outTime++;
                }
                if (std::abs(dt) < dt_max && ev.peak_height[j][0] < 30) {
                    if (wave_count_small_inTime < 100) {
                        // Store the waveforms
                        auto &wf = stored_wave_small_inTime[wave_count_small_inTime];
                        const int ns = std::min<int>(static_cast<int>(ev.nsamples[j]), fdec::MAX_SAMPLES);
                        wf.nsamples = ns;
                        wf.dt = dt;
                        wf.peak_height = ev.peak_height[j][0];
                        for (int s = 0; s < ns; ++s) {
                            wf.samples[s] = ev.samples[j][s];
                            if (s < 100) {
                                h1_wave_small_inTime[wave_count_small_inTime]->SetBinContent(s + 1, ev.samples[j][s]);
                            }
                        }
                    }
                    wave_count_small_inTime++;
                }
                if (ev.peak_height[j][0] > 100) {
                    ana.Analyze(ev.samples[j], ev.nsamples[j], wres);
                    h1_mu_crystal->Fill(wres.peaks_fit[0].mu);
                    h1_sigma_crystal->Fill(wres.peaks_fit[0].sigma);
                    h1_chi2_ndf_crystal->Fill(wres.peaks_fit[0].chi2_per_dof/wres.peaks_fit[0].A);
                    if (wave_count_big < 100) {
                        // Store the waveforms
                        auto &wf = stored_wave_big[wave_count_big];
                        const int ns = std::min<int>(static_cast<int>(ev.nsamples[j]), fdec::MAX_SAMPLES);
                        wf.nsamples = ns;
                        wf.dt = dt;
                        wf.peak_height = ev.peak_height[j][0];
                        for (int s = 0; s < ns; ++s) {
                            wf.samples[s] = ev.samples[j][s];
                            if (s < 100) {
                                h1_wave_big[wave_count_big]->SetBinContent(s + 1, ev.samples[j][s]);
                            }
                        }
                    }
                    wave_count_big++;
                }
            }
        }
        // fill the histograms for the number of fired modules in the 1st and 2nd layers
        if (!seed_not_clean) {
            h1_cluster_energy->Fill(hits[0].energy);
        }
    }

    // Draw waveform pages to PDFs: one waveform per page with WaveAnalyzer re-fit overlay.
    auto draw_waveforms_pdf = [&](const std::string &suffix,
                                  const std::string &title_prefix,
                                  TH1F *const hist_arr[100],
                                  const std::array<StoredWaveform, 100> &wave_arr,
                                  int wave_count) {
        const std::string pdf_name = output_path_name + "_" + suffix + ".pdf";
        TCanvas c_wave(Form("c_wave_%s", suffix.c_str()), title_prefix.c_str(), 900, 650);
        c_wave.Print((pdf_name + "[").c_str());

        const int nwave = std::min(wave_count, 100);
        for (int i = 0; i < nwave; ++i) {
            auto *h = hist_arr[i];
            const auto &wf = wave_arr[i];
            if (!h || wf.nsamples <= 0) continue;

            std::optional<TF1> fit_curve;
            std::optional<TMarker> reco_point;

            h->SetTitle(Form("%s %d;Samples;ADC", title_prefix.c_str(), i));
            h->SetStats(0);
            h->SetLineColor(kBlue + 1);
            h->SetLineWidth(1);
            h->SetMarkerStyle(20);
            h->SetMarkerSize(0.6);
            h->SetMarkerColor(kBlue + 1);
            h->GetXaxis()->SetRangeUser(20.0, 70.0);
            h->Draw("P");

            // stack-allocated scratch buffer for smoothed waveform
            float smoothed[100];
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

            if (fit_ok) {
                const float ped_mean = wave_res.ped.mean;
                const float A = ln_fit.A;
                const float t0 = ln_fit.t0;
                const float mu = ln_fit.mu;
                const float sigma = ln_fit.sigma;

                fit_curve.emplace(Form("f_wave_fit_%s_%d", suffix.c_str(), i),
                    [=](double *x, double *) {
                        const float sample = static_cast<float>(x[0]);
                        return static_cast<double>(
                            fdec::WaveAnalyzer::log_normal_pulse_value(
                                sample, ped_mean, A, t0, mu, sigma));
                    },
                    static_cast<double>(fit_left), static_cast<double>(fit_right), 0);
                fit_curve->SetLineColor(kRed + 1);
                fit_curve->SetLineWidth(3);
                fit_curve->SetNpx(500);
                fit_curve->Draw("same");
            }

            const float clk_ns = (ana.cfg.clk_mhz > 0.0f) ? (1000.0f / ana.cfg.clk_mhz) : 4.0f;
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
                    const float ped_mean = wave_res.ped.mean;
                    const float A = ln_fit.A;
                    const float t0 = ln_fit.t0;
                    const float mu = ln_fit.mu;
                    const float sigma = ln_fit.sigma;
                    reco_y = fdec::WaveAnalyzer::log_normal_pulse_value(
                        reco_sample, ped_mean, A, t0, mu, sigma);
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

    draw_waveforms_pdf("wave_small_inTime", "Small In-Time Waveform", h1_wave_small_inTime,
                       stored_wave_small_inTime, wave_count_small_inTime);
    draw_waveforms_pdf("wave_small_outTime", "Small Out-Time Waveform", h1_wave_small_outTime,
                       stored_wave_small_outTime, wave_count_small_outTime);
    draw_waveforms_pdf("wave_big", "Big Waveform", h1_wave_big,
                       stored_wave_big, wave_count_big);
    draw_waveforms_pdf("wave_veto_big", "Veto Big Waveform", h1_wave_veto_big,
                       stored_wave_veto_big, wave_count_veto_big);
    draw_waveforms_pdf("wave_veto_small", "Veto Small Waveform", h1_wave_veto_small,
                       stored_wave_veto_small, wave_count_veto_small);

    // Save the histograms to a root file
    TFile *output_file = new TFile((output_path_name + ".root").c_str(), "RECREATE");
    h1_cluster_energy->Write();
    h1_time_diff_seed->Write();
    h2_height_vs_dtime->Write();
    h1_mu_crystal->Write();
    h1_sigma_crystal->Write();
    h1_chi2_ndf_crystal->Write();
    h1_mu_veto->Write();
    h1_sigma_veto->Write();
    h1_chi2_ndf_veto->Write();
    output_file->cd();
    output_file->mkdir("waveforms_small_inTime");
    for (int i = 0; i < std::min(wave_count_small_inTime, 100); ++i) h1_wave_small_inTime[i]->Write();
    output_file->cd();
    output_file->mkdir("waveforms_small_outTime");
    for (int i = 0; i < std::min(wave_count_small_outTime, 100); ++i) h1_wave_small_outTime[i]->Write();
    output_file->cd();
    output_file->mkdir("waveforms_big");
    for (int i = 0; i < std::min(wave_count_big, 100); ++i) h1_wave_big[i]->Write();
    output_file->cd();
    output_file->mkdir("waveforms_veto_big");
    for (int i = 0; i < std::min(wave_count_veto_big, 100); ++i) h1_wave_veto_big[i]->Write();
    output_file->cd();
    output_file->mkdir("waveforms_veto_small");
    for (int i = 0; i < std::min(wave_count_veto_small, 100); ++i) h1_wave_veto_small[i]->Write();
    output_file->Close();

}