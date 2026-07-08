
//=============================================================================
// det_calib — detector position calibration via Møller scattering
//
// Process: Read evio files, decode and reconstruct, judge if it's 
// Moller events, analyze and fill histgroams, find out the detector alignment
//
// Usage:
//   det_calib <evio_file_or_dir> [more files/dirs...]
//             -o output_dir
//=============================================================================

#include "Replay.h"
#include "PhysicsTools.h"
#include "HyCalSystem.h"
#include "EventData.h"
#include "EventData_io.h"
#include "ConfigSetup.h"
#include "InstallPaths.h"
#include "MatchingTools.h"
#include "PipelineBuilder.h"
#include "PulseTemplateStore.h"
#include "gain_factor.h"

#include <TClass.h>
#include <TROOT.h>
#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TString.h>
#include <TSystem.h>
#include <TLatex.h>
#include <TCanvas.h>
#include <TF1.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdlib>
#include <filesystem>
#include <getopt.h>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <unordered_map>
#include <vector>

#ifndef DATABASE_DIR
#define DATABASE_DIR "."
#endif

using namespace analysis;
namespace fs = std::filesystem;

using EventVars_Recon = prad2::ReconEventData;

// ── forward declarations ──────────────────────────────────────────────────
static std::vector<std::string> collectEvioFiles(const std::string &path);
static bool ProcessEVIO(const std::string &input_evio,
                        RunConfig &run_config,
                        const std::string &db_dir,
                        const std::string &run_config_file,
                        const std::string &daq_config_file,
                        const std::string &hycal_map_file,
                        const std::string &gem_ped_file,
                        float zerosup_override,
                        int max_events,
                        MollerData &hycal_mollers,
                        std::array<MollerData, 4> &gem_mollers);

// ── file helpers ──────────────────────────────────────────────────────────
static std::vector<std::string> collectEvioFiles(const std::string &path)
{
    std::vector<std::string> files;
    if (fs::is_directory(path)) {
        for (auto &e : fs::directory_iterator(path)) {
            if (e.is_regular_file() &&
                e.path().filename().string().find(".evio") != std::string::npos)
                files.push_back(e.path().string());
        }
        std::sort(files.begin(), files.end());
    } else {
        files.push_back(path);
    }
    return files;
}

// ── Main ──────────────────────────────────────────────────────────────────
int main(int argc, char *argv[])
{
    ROOT::EnableThreadSafety();
    TH1::AddDirectory(kFALSE);
    TClass::GetClass("TTree");
    TClass::GetClass("TFile");
    TClass::GetClass("TBranch");

    std::string daq_config, daq_map, gem_ped_file, output_dir, run_config;
    float zerosup_override = 0.f;
    int  max_files   = -1;
    int  num_threads = 4;
    int  max_events  = -1;

    int opt;
    while ((opt = getopt(argc, argv, "o:f:n:j:c:C:d:g:z:")) != -1) {
        switch (opt) {
            case 'o': output_dir       = optarg; break;
            case 'f': max_files        = std::atoi(optarg); break;
            case 'n': max_events       = std::atoi(optarg); break;
            case 'j': num_threads      = std::atoi(optarg); break;
            case 'c': run_config       = optarg; break;
            case 'C': daq_config       = optarg; break;
            case 'd': daq_map          = optarg; break;
            case 'g': gem_ped_file     = optarg; break;
            case 'z': zerosup_override = std::atof(optarg); break;
            default: return 1;
        }
    }

    // collect input EVIO files (files, directories, or mixed)
    std::vector<std::string> evio_files;
    for (int i = optind; i < argc; ++i) {
        auto f = collectEvioFiles(argv[i]);
        evio_files.insert(evio_files.end(), f.begin(), f.end());
    }

    if (evio_files.empty() || output_dir.empty()) {
        std::cerr <<
            "Usage: det_calib <evio_file_or_dir> [more files/dirs...] -o output_dir\n"
            "       [-f max_files] [-j threads] [-n max_events]\n"
            "       [-C daq_config.json] [-d hycal_map.json]\n"
            "       [-c run_config.json] [-g gem_pedestal.json]\n";
        return 1;
    }

    int num_files = static_cast<int>(evio_files.size());
    if (max_files > 0) num_files = std::min(num_files, max_files);
    num_threads = std::max(1, std::min(num_threads, num_files));

    std::string db_dir = prad2::resolve_data_dir(
        "PRAD2_DATABASE_DIR",
        {"../share/prad2evviewer/database"},
        DATABASE_DIR);
    if (daq_config.empty()) daq_config = db_dir + "/daq_config.json";
    if (daq_map.empty())    daq_map    = db_dir + "/hycal_map.json";
    if (run_config.empty()) run_config = db_dir + "/runinfo/general.json";

    int run_num = get_run_int(evio_files[0]);
    gRunConfig = LoadRunConfig(run_config, run_num);

    // Each worker writes only to the result slot belonging to its EVIO file.
    // This avoids locking the (potentially large) Moller vectors while events
    // are being reconstructed.  The main thread merges the slots after all
    // workers have finished.
    std::vector<MollerData> hycal_results(num_files);
    std::vector<std::array<MollerData, 4>> gem_results(num_files);
    std::vector<char> processed_ok(num_files, 0);
    std::atomic<int> next_file{0};
    std::atomic<int> errors{0};
    std::mutex io_mtx;

    auto worker = [&]() {
        while (true) {
            const int idx = next_file.fetch_add(1);
            if (idx >= num_files) break;

            RunConfig file_config = gRunConfig;
            const bool ok = ProcessEVIO(
                evio_files[idx], file_config, db_dir, run_config, daq_config,
                daq_map, gem_ped_file, zerosup_override, max_events,
                hycal_results[idx], gem_results[idx]);
            processed_ok[idx] = ok ? 1 : 0;

            std::lock_guard<std::mutex> lock(io_mtx);
            if (ok) {
                std::cout << "  [" << (idx + 1) << "/" << num_files << "] "
                          << evio_files[idx] << ": "
                          << hycal_results[idx].size() << " HyCal Mollers\n";
            } else {
                ++errors;
                std::cerr << "  [" << (idx + 1) << "/" << num_files
                          << "] FAILED: " << evio_files[idx] << "\n";
            }
        }
    };

    std::cout << "Processing " << num_files << " EVIO files with "
              << num_threads << " threads\n";
    std::vector<std::thread> threads;
    threads.reserve(num_threads);
    for (int i = 0; i < num_threads; ++i)
        threads.emplace_back(worker);
    for (auto &thread : threads)
        thread.join();

    // Aggregate in input-file order, so results are deterministic regardless
    // of the order in which worker threads completed.
    MollerData hycal_mollers;
    std::array<MollerData, 4> gem_mollers;
    size_t total_hycal = 0;
    std::array<size_t, 4> total_gem{};
    for (int i = 0; i < num_files; ++i) {
        if (!processed_ok[i]) continue;
        total_hycal += hycal_results[i].size();
        for (size_t det = 0; det < gem_mollers.size(); ++det)
            total_gem[det] += gem_results[i][det].size();
    }
    hycal_mollers.reserve(total_hycal);
    for (size_t det = 0; det < gem_mollers.size(); ++det)
        gem_mollers[det].reserve(total_gem[det]);

    for (int i = 0; i < num_files; ++i) {
        if (!processed_ok[i]) continue;
        hycal_mollers.insert(hycal_mollers.end(),
                             hycal_results[i].begin(), hycal_results[i].end());
        for (size_t det = 0; det < gem_mollers.size(); ++det) {
            gem_mollers[det].insert(gem_mollers[det].end(),
                                    gem_results[i][det].begin(),
                                    gem_results[i][det].end());
        }
    }

    std::cout << "Collected " << hycal_mollers.size()
              << " HyCal Moller events\n";
    for (size_t det = 0; det < gem_mollers.size(); ++det)
        std::cout << "  GEM " << det << ": " << gem_mollers[det].size()
                  << " Moller events\n";

    // TODO: analyze the combined hycal_mollers and gem_mollers data.
    return errors.load() == 0 ? 0 : 1;
}

static bool ProcessEVIO (const std::string &input_evio, RunConfig &gRunConfig,
                                const std::string &db_dir,
                                const std::string &run_config_file,
                                const std::string &daq_config_file,
                                const std::string &hycal_map_file,
                                const std::string &gem_ped_file,
                                const float zerosup_override,
                                const int max_events,
                                MollerData &hycal_mollers,
                                std::array<MollerData, 4> &gem_mollers)
{
    fdec::HyCalSystem                 hycal;
    gem::GemSystem                    gem_sys;
    fdec::ClusterConfig               cluster_cfg;
    prad2::HyCalTimeCuts              hc_time_cuts;
    DetectorTransform                 hycal_transform;
    std::array<DetectorTransform, 4>  gem_transforms;
    std::unordered_map<int, int>      roc_to_crate;

    prad2::Pipeline pipeline;
    try {
        pipeline = prad2::PipelineBuilder()
            .set_database_dir(db_dir)
            .set_daq_config(daq_config_file)
            .set_runinfo(run_config_file)
            .set_hycal_map(hycal_map_file)
            .set_gem_pedestal(gem_ped_file)
            .set_run_number_from_evio(input_evio)
            .set_log_stream(&std::cerr)
            .build();
    } catch (const std::exception &e) {
        std::cerr << "Replay: setup failed for " << input_evio
                  << ": " << e.what() << "\n";
        return false;
    }

    evc::DaqConfig daq_cfg = std::move(pipeline.daq_cfg);
    gRunConfig      = pipeline.run_cfg;
    hycal            = std::move(pipeline.hycal);
    gem_sys          = std::move(pipeline.gem);
    cluster_cfg      = pipeline.hycal_cluster_cfg;
    hc_time_cuts     = std::move(pipeline.hycal_time_cuts);
    hycal_transform  = pipeline.hycal_transform;
    gem_transforms   = pipeline.gem_transforms;

    // ROC→crate map from the same DAQ config the builder consumed.
    for (const auto &re : daq_cfg.roc_tags) {
        if (re.crate < 0) continue;
        if (!re.type.empty() && re.type != "roc" && re.type != "gem") continue;
        roc_to_crate[re.tag] = re.crate;
    }

    if (zerosup_override > 0.f)
        gem_sys.SetZeroSupThreshold(zerosup_override);

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);
    gem::GemCluster      gem_clusterer;
    MatchingTools        matching;
    matching.SetMatchRange(gRunConfig.matching_radius);
    matching.SetSquareSelection(gRunConfig.matching_use_square);
    matching.SetEnergyDependent(gRunConfig.matching_energy_dependent);
    matching.SetMatchSigma(gRunConfig.matching_sigma);
    PhysicsTools physics(hycal);
    //open EVIO file and output ROOT file
    evc::EvChannel ch;
    ch.SetConfig(daq_cfg);

    if (ch.OpenAuto(input_evio) != evc::status::success) {
        std::cerr << "Replay: cannot open " << input_evio << "\n";
        return false;
    }

    auto ev = std::make_unique<EventVars_Recon>();

    //initialize tools for event decoder and cluster reconstruction
    auto event = std::make_unique<fdec::EventData>();
    auto ssp_evt = std::make_unique<ssp::SspEventData>();
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

    int total = 0;
    int processed_events = 0;

    int run_num = get_run_int(input_evio);
    std::string gain_data_dir = gRunConfig.gain_data_dir;
    if (gain_data_dir.empty())
        gain_data_dir = db_dir + "/gain_factor";
    else if (!fs::path(gain_data_dir).is_absolute())
        gain_data_dir = db_dir + "/" + gain_data_dir;
    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(
        gain_data_dir + "/gain_correction", run_num);

    // Per-detector lab transforms — set up by either branch of the detector
    // wiring above (PipelineBuilder for PRad-II, BuildLabTransforms for PRad-1).
    const auto &hc_xform = hycal_transform;
    const auto &g_xform  = gem_transforms;

    while ((max_events <= 0 || processed_events < max_events)
           && ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            if (max_events > 0 && processed_events >= max_events) break;
            event->clear();
            ssp_evt->clear();
            clusterer.Clear();
            if (!ch.DecodeEvent(ie, *event, ssp_evt.get())) continue;
            ++processed_events;

            *ev = EventVars_Recon{};
            ev->event_num    = event->info.event_number;
            ev->trigger_bits = event->info.trigger_bits;

            bool is_sum      = (ev->trigger_bits & prad2::TBIT_sum)   != 0;
            if (!is_sum) continue;

            // Per-event gain correction (time-series lookup by event number).
            const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(ev->event_num));

            // decode FADC250 and reconstruct HyCal data
            int nch = 0;
            for (int r = 0; r < event->nrocs; ++r) {
                auto &roc = event->rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                if (cit == roc_to_crate.end()) continue;
                int crate = cit->second;
                for (int s = 0; s < fdec::MAX_SLOTS; ++s) {
                    if (!roc.slots[s].present) continue;
                    for (int c = 0; c < 64; ++c) { //should be 16, a bigger number to adapt PRad1 data
                        if (!(roc.slots[s].channel_mask & (1ull << c))) continue;
                        auto &cd = roc.slots[s].channels[c];
                        if (cd.nsamples <= 0) continue;

                        const auto *mod = hycal.module_by_daq(crate, s, c);
                        if (!mod || !mod->is_hycal()) continue;
                        // Per-ID gain correction: average of LMS 2/3 channels.
                        const float gain = (mod->id > 1000)
                            ? (gain_corr.w[mod->id - 1000].corr[1] + gain_corr.w[mod->id - 1000].corr[2]) / 2.0f
                            : gain_corr.g[mod->id].avg;

                        ana.SetChannelKey(roc.tag, s, c);
                        ana.Analyze(cd.samples, cd.nsamples, wres);
                        if (wres.npeaks <= 0) continue;

                        const auto hc_win = hc_time_cuts.at(mod->index);
                        if (cluster_cfg.seed_time_window > 0.f) {
                            // Multi-pulse mode: push every peak inside the trigger
                            // window into the clusterer; the seed-anchored timing
                            // coincidence cut is applied inside HyCalCluster.
                            for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                const auto &pk = wres.peaks[p];
                                if (pk.time <= hc_win.lo) continue;
                                if (pk.time >= hc_win.hi) continue;
                                float adc = pk.integral * gain;
                                float energy = static_cast<float>(mod->energize(adc));
                                clusterer.AddHit(mod->index, energy, pk.time);
                                ev->total_energy += energy;
                                nch++;
                            }
                        } else {
                            // Legacy: pick the largest in-window peak as the single
                            // module hit, time field unused downstream.
                            int bestIdx = -1;
                            float bestHeight = -1.f;
                            for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                const auto &pk = wres.peaks[p];
                                if (pk.time > hc_win.lo &&
                                    pk.time < hc_win.hi &&
                                    pk.height > bestHeight) {
                                    bestHeight = pk.height;
                                    bestIdx = p;
                                }
                            }
                            if (bestIdx < 0) continue;
                            float adc = wres.peaks[bestIdx].integral * gain;
                            float energy = static_cast<float>(mod->energize(adc));
                            clusterer.AddHit(mod->index, energy, wres.peaks[bestIdx].time);
                            ev->total_energy += energy;
                            nch++;
                        }
                    }
                }
            }
            if(nch > 200) continue; // too many hits, likely noise, skip the event

            clusterer.FormClusters();
            std::vector<fdec::ClusterHit> hits;
            clusterer.ReconstructHits(hits);
            //HyCal event reconstrued
            ev->n_clusters = std::min((int)hits.size(), prad2::kMaxClusters);
            for (int i = 0; i < ev->n_clusters; ++i) {
                ev->cl_nblocks[i] = hits[i].nblocks;
                ev->cl_time[i]    = hits[i].time;
                //transform the cluster positions to the lab coordinate
                HCHit local_hit = {hits[i].x, hits[i].y, fdec::shower_depth(hits[i].center_id, hits[i].energy),
                    hits[i].energy, static_cast<uint16_t>(hits[i].center_id), hits[i].flag};
                analysis::ApplyToLab(hc_xform, local_hit);
                GetProjection(local_hit, gRunConfig.hycal_z);
                ev->cl_x[i] = local_hit.x;
                ev->cl_y[i] = local_hit.y;
                ev->cl_z[i] = local_hit.z;
                ev->cl_energy[i] = local_hit.energy;
                ev->cl_linear_corr[i] = hits[i].linear_corr;
                ev->cl_center[i] = local_hit.center_id;
                ev->cl_flag[i] = local_hit.flag;
            }

            //decode GEM data and reconstruct GEM hits
            if(gem_sys.GetNDetectors() > 0){
                gem_sys.Clear();
                gem_sys.ProcessEvent(*ssp_evt);
                gem_sys.Reconstruct(gem_clusterer);
                auto &all_hits = gem_sys.GetAllHits();
                ev->n_gem_hits = 0;
                for (const auto &h : all_hits) {
                    if (ev->n_gem_hits >= prad2::kMaxGemHits) break;
                    if (h.det_id < 0 || h.det_id >= 4) continue;
                    const int i = ev->n_gem_hits++;
                    ev->det_id[i] = h.det_id;
                    ev->gem_x_charge[i] = h.x_charge;
                    ev->gem_y_charge[i] = h.y_charge;
                    ev->gem_x_peak[i] = h.x_peak;
                    ev->gem_y_peak[i] = h.y_peak;
                    ev->gem_x_size[i] = h.x_size;
                    ev->gem_y_size[i] = h.y_size;
                    ev->gem_x_mTbin[i] = h.x_max_timebin;
                    ev->gem_y_mTbin[i] = h.y_max_timebin;
                    //transform the GEM hit positions to the lab coordinate
                    GEMHit local_hit = {h.x, h.y, 0.f, static_cast<uint8_t>(h.det_id)};
                    int d = local_hit.det_id;
                    if (d >= 0 && d < 4) {
                        analysis::ApplyToLab(g_xform[d], local_hit);
                    }
                    ev->gem_x[i] = local_hit.x;
                    ev->gem_y[i] = local_hit.y;
                    ev->gem_z[i] = local_hit.z;
                }

                // Perform matching between HyCal clusters and GEM hits
                //store all the hits on HyCal and GEMs in this event
                std::vector<HCHit> hc_hits;
                std::vector<GEMHit> gem_hits[4]; // separate vector for each GEM
                for (int i = 0; i < ev->n_clusters; ++i)
                    hc_hits.push_back({ev->cl_x[i], ev->cl_y[i], ev->cl_z[i], ev->cl_energy[i], ev->cl_center[i], ev->cl_flag[i]});
                for (int i = 0; i < ev->n_gem_hits; ++i) {
                    const int det_id = ev->det_id[i];
                    if (det_id < 0 || det_id >= 4) continue;
                    gem_hits[det_id].push_back(GEMHit{
                        ev->gem_x[i], ev->gem_y[i], ev->gem_z[i],
                        static_cast<uint8_t>(det_id)});
                }
                
                // already transform to the coordinates

                std::vector<MatchHit> matched_hits = matching.Match(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]);
                std::vector<MatchHit_perChamber> matched_hits_chamber = matching.MatchPerChamber(hc_hits, gem_hits[0], gem_hits[1], gem_hits[2], gem_hits[3]); 
                
                for (const auto &m : matched_hits_chamber) {
                    const int cl_idx = m.hycal_idx;
                    if (cl_idx < 0 || cl_idx >= prad2::kMaxClusters) continue;
                    for(int j = 0; j < 4; j++){
                        ev->matchGEMx[cl_idx][j] = m.gem_hits[j][0];
                        ev->matchGEMy[cl_idx][j] = m.gem_hits[j][1];
                        ev->matchGEMz[cl_idx][j] = m.gem_hits[j][2];
                    }
                    ev->matchFlag[cl_idx] = m.mflag;
                }

                ev->matchNum = std::min((int)matched_hits.size(), prad2::kMaxClusters);
                for (int i = 0; i < ev->matchNum; i++){
                    // save the matched GEM hit (must 2 matchings) info in mHit_ arrays for quick check
                    ev->mHit_E[i] = matched_hits[i].hycal_hit.energy;
                    ev->mHit_x[i] = matched_hits[i].hycal_hit.x;
                    ev->mHit_y[i] = matched_hits[i].hycal_hit.y;
                    ev->mHit_z[i] = matched_hits[i].hycal_hit.z;
                    for(int j = 0; j < 2; j++) {
                        ev->mHit_gx[i][j] =  matched_hits[i].gem[j].x;
                        ev->mHit_gy[i][j] =  matched_hits[i].gem[j].y;
                        ev->mHit_gz[i][j] =  matched_hits[i].gem[j].z;
                        ev->mHit_gid[i][j] = matched_hits[i].gem[j].det_id; // placeholder for GEM hit ID if needed
                    }
                }
            }

            // select events and analyze
            if (ev->n_clusters != 2) continue;
            if (ev->matchNum != 2) continue;
            if (ev->cl_nblocks[0] < 3 || ev->cl_nblocks[1] < 3) continue;
            if (std::fabs(ev->cl_x[0]) < 20.75f * 2.5f && std::fabs(ev->cl_y[0]) < 20.75f * 2.5f) continue;
            if (std::fabs(ev->cl_x[1]) < 20.75f * 2.5f && std::fabs(ev->cl_y[1]) < 20.75f * 2.5f) continue;
            if (std::fabs(ev->cl_x[0]) > 20.75f * 16.f || std::fabs(ev->cl_y[0]) > 20.75f * 16.f) continue;
            if (std::fabs(ev->cl_x[1]) > 20.75f * 16.f || std::fabs(ev->cl_y[1]) > 20.75f * 16.f) continue;

            MollerEvent m_hc(
                {ev->cl_x[0], ev->cl_y[0], ev->cl_z[0], ev->cl_energy[0]},
                {ev->cl_x[1], ev->cl_y[1], ev->cl_z[1], ev->cl_energy[1]});

            constexpr float kRadToDeg = 57.29577951308232f;
            float theta1 = std::atan2(std::hypot(m_hc.first.x, m_hc.first.y), m_hc.first.z) * kRadToDeg;
            float theta2 = std::atan2(std::hypot(m_hc.second.x, m_hc.second.y), m_hc.second.z) * kRadToDeg;
            if (!physics.isMoller_kinematic(theta1, m_hc.first.E,
                                            theta2, m_hc.second.E,
                                            gRunConfig.Ebeam, 0.035f))
                continue;
            if (std::fabs(physics.GetMollerPhiDiff(m_hc)) > 6.f) continue;

            //add some scattering angle cuts for 0.7GeV
            if(gRunConfig.Ebeam > 0.f && gRunConfig.Ebeam < 1000.f) {
                if(theta1 < 1.5f || theta2 < 1.5f) continue;
            }

            hycal_mollers.push_back(m_hc);

            for(int did = 0; did < 4; did ++){
                if (((ev->matchFlag[0] & (1u << did)) != 0)
                    && ((ev->matchFlag[1] & (1u << did)) != 0)) {
                    gem_mollers[did].push_back(MollerEvent(
                        {ev->matchGEMx[0][did], ev->matchGEMy[0][did], ev->matchGEMz[0][did], ev->cl_energy[0]},
                        {ev->matchGEMx[1][did], ev->matchGEMy[1][did], ev->matchGEMz[1][did], ev->cl_energy[1]}));
                }
            }

            total++;
            if (total % 1000 == 0)
                std::cerr << "\rFind: " << total << " moller events on HyCal" << std::flush;
        }
    }
    std::cerr << "\rReplay: " << total << " moller events reconstructed on HyCal\n";

    return true;
}
