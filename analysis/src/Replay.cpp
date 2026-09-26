//=============================================================================
// Replay.cpp — EVIO to ROOT tree conversion
//=============================================================================

#include "Replay.h"
#include "DaqConfig.h"
#include "EventData_io.h"
#include "PulseTemplateStore.h"
#include "HyCalSystem.h"
#include "GemSystem.h"
#include "HyCalCluster.h"
#include "GemCluster.h"
#include "MatchingTools.h"
#include "ConfigSetup.h"
#include "PipelineBuilder.h"
#include "RfTime.h"
#include "TdcDecoder.h"
#include "VtpDecoder.h"
#include "gain_factor.h"

#include <iostream>

namespace analysis {

namespace {

// GEM quality (SBS-style) recon-tree fillers.  Shared by ProcessWithRecon
// (EVIO) and ProcessRaw2Recon (raw ROOT) so both paths write identical
// values; the branches themselves are only booked with -gem_hit.

void fillGemHitQA(EventVars_Recon &ev, int i, const gem::GEMHit &h)
{
    ev.gem_x_time[i]     = h.x_time;
    ev.gem_y_time[i]     = h.y_time;
    ev.gem_xy_dt[i]      = h.time_diff;
    ev.gem_xy_asym[i]    = h.adc_asym;
    ev.gem_x_max_sdt[i]  = h.x_max_strip_dt;
    ev.gem_y_max_sdt[i]  = h.y_max_strip_dt;
    ev.gem_x_min_corr[i] = h.x_min_ts_corr;
    ev.gem_y_min_corr[i] = h.y_min_ts_corr;
}

// Append the (already filtered) 1D clusters of one detector plane to the
// per-cluster block.  Stops silently at kMaxGemClusters, like the hits.
// size / mTbin use the same narrowing as gem_x_size / gem_x_mTbin.
void appendGemClusters(EventVars_Recon &ev, int det, int plane,
                       const std::vector<gem::StripCluster> &cls)
{
    for (const auto &c : cls) {
        if (ev.n_gem_cl >= prad2::kMaxGemClusters) return;
        const int k = ev.n_gem_cl++;
        ev.gem_cl_det[k]       = static_cast<uint8_t>(det);
        ev.gem_cl_plane[k]     = static_cast<uint8_t>(plane);
        ev.gem_cl_size[k]      = static_cast<uint8_t>(c.hits.size());
        ev.gem_cl_mTbin[k]     = static_cast<uint8_t>(c.max_timebin);
        ev.gem_cl_pos[k]       = c.position;
        ev.gem_cl_peak[k]      = c.peak_charge;
        ev.gem_cl_charge[k]    = c.total_charge;
        ev.gem_cl_time[k]      = c.seed_time;
        ev.gem_cl_seed_peak[k] = c.seed_peak_adc;
        ev.gem_cl_seed_sum[k]  = c.seed_sum_adc;
        ev.gem_cl_max_sdt[k]   = c.max_strip_dt;
        ev.gem_cl_min_corr[k]  = c.min_ts_corr;
    }
}

// Slow-control side trees written next to the event tree: one DSC2 scaler
// row per SYNC physics event, one EPICS row per EPICS event and one row per
// CODA control event.  The PRESTART row carries the long DAQ-config text
// (0xE10E STRING bank); GO and END are recorded too so analysis can recover
// the run start/end time even when no PRESTART is in the input.  The rows
// are filled by prad2dec accessors; see EventData_io.h for the format and
// the join-by-event_number scheme.  Create with make_unique: the trees keep
// the addresses of the rows.
struct SideTrees {
    TTree *scalers = new TTree("scalers", "PRad2 DSC2 scaler readouts");
    TTree *epics   = new TTree("epics",   "PRad2 EPICS slow control");
    TTree *runinfo = new TTree("runinfo", "PRad2 control events / DAQ config");
    prad2::RawScalerData sc_row;
    prad2::RawEpicsData  ep_row;
    prad2::RawRunInfo    ri_row;

    SideTrees()
    {
        prad2::SetScalerWriteBranches (scalers, sc_row);
        prad2::SetEpicsWriteBranches  (epics,   ep_row);
        prad2::SetRunInfoWriteBranches(runinfo, ri_row);
    }
    SideTrees(const SideTrees &) = delete;
    SideTrees &operator=(const SideTrees &) = delete;

    // Record a control or EPICS event.  Returns false only for a physics
    // event, which is left to the caller.
    bool record(const evc::EvChannel &ch)
    {
        const auto et = ch.GetEventType();
        if (et == evc::EventType::Prestart ||
            et == evc::EventType::Go       ||
            et == evc::EventType::End)
        {
            // Scan() refreshed ch.Sync() with this event's run number /
            // unix_time / run_type; the config text only ships on PRESTART.
            std::string cfg_text;
            if (et == evc::EventType::Prestart)
                cfg_text = ch.ExtractDaqConfigText();
            prad2::FillRunInfoRow(ch.Sync(), cfg_text, ri_row);
            runinfo->Fill();
        } else if (et == evc::EventType::Epics) {
            const auto &rec = ch.Epics();
            if (rec.present) {
                prad2::FillEpicsRow(rec, ep_row);
                epics->Fill();
            }
        }
        return et != evc::EventType::Physics;
    }

    // DSC2 lives at the CODA-event level (one bank per Read() covering all
    // its sub-events), but the carrying sub-event has a unique event_number:
    // call this for the first decoded sub-event, before any per-event cut,
    // so the rows track SYNC arrivals 1:1.
    void recordScalers(const evc::EvChannel &ch, const fdec::EventInfo &info,
                       const evc::DaqConfig::DscScaler &cfg)
    {
        const auto &dsc = ch.Dsc();
        if (!dsc.present) return;
        prad2::FillScalerRow(dsc, ch.Sync(), info, cfg, sc_row);
        scalers->Fill();
    }

    void write()
    {
        scalers->Write();
        epics->Write();
        runinfo->Write();
    }
};

// Every leaf bank with one of `tags` in the current read group (the ones
// ch.Vtp() / ch.Tdc() decode), kept verbatim as a flat triple (parent ROC
// tag or 0, word count, concatenated words) so offline tools can re-decode
// it without the EVIO file and without a ROOT dictionary for nested STL
// collections.  VTP: up to ~9 banks (7 HyCal + 2 GEM VTPs; PRAD_CLUSTER
// TAG_EXP 0x1CC, TRIGGER 0x1D, ...).  TDC: V1190/V1290 banks (0x40 "rf": the
// divided CEBAF RF on slot 16, ch 0 + ch 8; see RawEventData.tdc_words and
// TdcDecoder.h).
struct BankSnapshot {
    std::vector<uint32_t> roc_tags, nwords, words;
};

BankSnapshot snapshotBanks(const evc::EvChannel &ch, const std::vector<uint32_t> &tags)
{
    BankSnapshot snap;
    for (uint32_t tag : tags) {
        ch.ForEachLeafBank(tag, [&](const evc::EvNode &n, uint32_t roc_tag) {
            const uint32_t *p = ch.GetData(n);
            snap.roc_tags.push_back(roc_tag);
            snap.nwords.push_back(static_cast<uint32_t>(n.data_words));
            snap.words.insert(snap.words.end(), p, p + n.data_words);
        });
    }
    return snap;
}

// Words of the first `tag` bank of the current read group (empty if none),
// e.g. the one 0xE10C SSP trigger bank shared by all its sub-events.
std::vector<uint32_t> snapshotFirstBank(const evc::EvChannel &ch, uint32_t tag)
{
    std::vector<uint32_t> words;
    if (const auto *n = ch.FindFirstByTag(tag)) {
        const uint32_t *p = ch.GetData(*n);
        words.assign(p, p + n->data_words);
    }
    return words;
}

// Trigger selection of the recon replays.
// TODO: use config-driven trigger filter (monitor_config.json "physics"
// section accept_trigger_bits/reject_trigger_bits) instead of hardcoded bits.
struct TriggerSel {
    bool sum = false, lms = false, alpha = false, cl3 = false;
    bool keep = false;
};

// Keeps sum (every event with `random`), LMS, alpha and 3-cluster triggers.
// For X17 blind analysis an event with the 3-cluster bit is only kept when
// its event number ends in 8 (lucky number decided by students); the
// raw-sum events stay for calibration and monitoring.
TriggerSel selectTrigger(uint32_t bits, int event_num, bool random, bool x17, bool x17_blind)
{
    TriggerSel t;
    t.sum   = random || (bits & prad2::TBIT_sum) != 0;
    t.lms   = (bits & prad2::TBIT_lms)   != 0;
    t.alpha = (bits & prad2::TBIT_alpha) != 0;
    t.cl3   = (bits & prad2::TBIT_3cl)   != 0;
    t.keep  = (t.sum || t.lms || t.alpha || t.cl3)
           && !(x17 && x17_blind && t.cl3 && event_num % 10 != 8);
    return t;
}

// VTP PRAD_CLUSTER online trigger clusters; center is the HyCal module id.
void fillVtpClusters(EventVars_Recon &ev, const vtp::VtpEventData &vtp_event)
{
    for (int i = 0; i < vtp_event.n_prad_clusters && ev.vtp_cl_n < vtp::MAX_PRAD_CLUSTERS; ++i) {
        const auto &cl = vtp_event.prad_clusters[i];
        ev.vtp_cl_time[ev.vtp_cl_n]   = cl.time;
        ev.vtp_cl_energy[ev.vtp_cl_n] = cl.energy;
        ev.vtp_cl_center[ev.vtp_cl_n] = cl.hycal_id();
        ev.vtp_cl_blocks[ev.vtp_cl_n] = cl.nhits;
        ++ev.vtp_cl_n;
    }
}

// Decode the RF reference from a TDC bank snapshot into the rf_n_a/_b and
// rf_ns_a/_b (leading-edge ns) branches; returns it for the per-cluster
// cl_dt_rf.
tdc::RfTimeData fillRf(EventVars_Recon &ev, const std::vector<uint32_t> &roc_tags,
                       const std::vector<uint32_t> &nwords,
                       const std::vector<uint32_t> &words)
{
    tdc::RfTimeData rf;
    tdc::RfTimeDecoder::DecodeReplay(roc_tags, nwords, words, rf);
    ev.rf_n_a = static_cast<uint8_t>(rf.n_a);
    ev.rf_n_b = static_cast<uint8_t>(rf.n_b);
    std::copy(rf.ns_a, rf.ns_a + rf.n_a, ev.rf_ns_a);
    std::copy(rf.ns_b, rf.ns_b + rf.n_b, ev.rf_ns_b);
    return rf;
}

// HyCal clusters (at most kMaxClusters): lab-frame position at shower depth
// projected to hycal_z, and the RF time difference.
void fillClusters(EventVars_Recon &ev, const std::vector<fdec::ClusterHit> &hits,
                  const DetectorTransform &hycal_xform, float hycal_z,
                  const tdc::RfTimeData &rf, const fdec::HyCalSystem &hycal,
                  const prad2::HyCalRfOffsets &rf_offsets)
{
    ev.n_clusters = std::min((int)hits.size(), prad2::kMaxClusters);
    for (int i = 0; i < ev.n_clusters; ++i) {
        const auto &h = hits[i];
        ev.cl_nblocks[i] = h.nblocks;
        ev.cl_npos[i]    = h.npos;
        ev.cl_time[i]    = h.time;
        HCHit lab = ClusterToLab(hycal_xform, h);
        GetProjection(lab, hycal_z);
        ev.cl_x[i]           = lab.x;
        ev.cl_y[i]           = lab.y;
        ev.cl_z[i]           = lab.z;
        ev.cl_energy[i]      = lab.energy;
        ev.cl_linear_corr[i] = h.linear_corr;
        ev.cl_bias_corr[i]   = h.bias_corr;
        ev.cl_center[i]      = lab.center_id;
        ev.cl_flag[i]        = lab.flag;

        // Per-cluster RF Δt — fold (cl_time − nearest_a) onto
        // (−T_RF/2, T_RF/2], then subtract per-module offset and
        // re-fold.  NaN when rf has no ch-A hits this event
        // (apply() preserves NaN through both steps).
        const float dt0 = prad2::ClusterDeltaRf(h.time, rf);
        const auto *mod = hycal.module_by_id(h.center_id);
        ev.cl_dt_rf[i] = rf_offsets.apply(mod ? mod->index : -1, dt0);
    }
}

// GEM 2D hits (at most kMaxGemHits): charge, size and timing QA, and the
// lab-frame position.
void fillGemHits(EventVars_Recon &ev, const std::vector<gem::GEMHit> &hits,
                 const std::array<DetectorTransform, 4> &gem_xforms)
{
    ev.n_gem_hits = std::min((int)hits.size(), prad2::kMaxGemHits);
    for (int i = 0; i < ev.n_gem_hits; ++i) {
        const auto &h = hits[i];
        ev.det_id[i]       = h.det_id;
        ev.gem_x_charge[i] = h.x_charge;
        ev.gem_y_charge[i] = h.y_charge;
        ev.gem_x_peak[i]   = h.x_peak;
        ev.gem_y_peak[i]   = h.y_peak;
        ev.gem_x_size[i]   = h.x_size;
        ev.gem_y_size[i]   = h.y_size;
        ev.gem_x_mTbin[i]  = h.x_max_timebin;
        ev.gem_y_mTbin[i]  = h.y_max_timebin;
        fillGemHitQA(ev, i, h);
        const GEMHit lab = GemHitToLab(gem_xforms, h);
        ev.gem_x[i] = lab.x;
        ev.gem_y[i] = lab.y;
        ev.gem_z[i] = lab.z;
    }
}

// hycal_map.json category -> on-disk module_type code.
prad2::ModuleType toModuleType(fdec::ModuleType t)
{
    switch (t) {
        case fdec::ModuleType::PbGlass: return prad2::MOD_PbGlass;
        case fdec::ModuleType::PbWO4:   return prad2::MOD_PbWO4;
        case fdec::ModuleType::Veto:    return prad2::MOD_VETO;
        case fdec::ModuleType::LMS:     return prad2::MOD_LMS;
        default:                        return prad2::MOD_UNKNOWN;
    }
}

} // anonymous namespace

void Replay::LoadHyCalMap(const std::string &json_path)
{
    if (!hycal_map_.Init(json_path)) {
        std::cerr << "Replay: cannot load HyCal map: " << json_path << "\n";
        return;
    }
    int n_daq = 0;
    for (int i = 0; i < hycal_map_.module_count(); ++i)
        if (hycal_map_.module(i).daq.crate >= 0) ++n_daq;
    std::cerr << "Replay: loaded " << hycal_map_.module_count()
              << " modules (" << n_daq << " with daq) from "
              << json_path << "\n";
}

std::string Replay::moduleName(int roc, int slot, int ch) const
{
    const auto *m = hycal_map_.module_by_daq(roc, slot, ch);
    return m ? m->name : "";
}

prad2::ModuleType Replay::moduleType(int roc, int slot, int ch) const
{
    const auto *m = hycal_map_.module_by_daq(roc, slot, ch);
    return m ? toModuleType(m->type) : prad2::MOD_UNKNOWN;
}

int Replay::moduleID(int roc, int slot, int ch) const
{
    // Globally-unique ID encoding — see RawEventData docs.  The numeric
    // ranges are deliberately disjoint so HyCalSystem::module_by_id(...)
    // returns nullptr for Veto/LMS, letting existing HyCal consumers
    // skip them via their existing nullptr / is_hycal() checks.
    const auto *m = hycal_map_.module_by_daq(roc, slot, ch);
    if (!m) return -1;
    const std::string &name = m->name;
    switch (m->type) {
        case fdec::ModuleType::PbGlass:
        case fdec::ModuleType::PbWO4:
            return m->id;

        case fdec::ModuleType::Veto:
            if (name.size() >= 2 && name[0] == 'V')
                try { return prad2::kVetoIdBase + std::stoi(name.substr(1)); } catch (...) {}
            return -1;

        case fdec::ModuleType::LMS:
            if (name == "LMSPin") return prad2::kLmsIdBase;
            if (name.rfind("LMS", 0) == 0 && name.size() == 4 && name[3] >= '1' && name[3] <= '9')
                return prad2::kLmsIdBase + (name[3] - '0');
            return -1;

        default:
            return -1;
    }
}

std::tuple<int, int, int> Replay::moduleLocation(int module_id) const
{
    using prad2::kVetoIdBase;
    using prad2::kLmsIdBase;
    const fdec::Module *m =
          module_id <  kVetoIdBase ? hycal_map_.module_by_id(module_id)
        : module_id <  kLmsIdBase  ? hycal_map_.module_by_name("V" + std::to_string(module_id - kVetoIdBase))
        : module_id == kLmsIdBase  ? hycal_map_.module_by_name("LMSPin")
        : hycal_map_.module_by_name("LMS" + std::to_string(module_id - kLmsIdBase));
    if (!m || m->daq.crate < 0) return {-1, -1, -1};
    return {m->daq.crate, m->daq.slot, m->daq.channel};
}

bool Replay::Process(const std::string &input_evio, const std::string &output_root, RunConfig &gRunConfig,
                     const std::string &db_dir, const std::string &recon_config_file,
                     int max_events, bool write_peaks , const std::string &daq_config_file,
                     const float zerosup_override,bool Ecalib, bool noWaveform)
{   
    // Detectors flow through PipelineBuilder so the wiring stays in one place
    // (see prad2det/include/PipelineBuilder.h).  daq_cfg_ moves through the
    // builder (which then attaches map paths) and comes back populated with
    // everything the per-event loop needs.
    std::string hycal_map_override = daq_cfg_.hycal_map_file;
    std::string gem_map_override   = daq_cfg_.gem_map_file;

    prad2::Pipeline pipeline = prad2::PipelineBuilder()
        .set_recon_config(recon_config_file)
        .set_database_dir(db_dir)
        .set_loaded_daq_config(std::move(daq_cfg_))
        .set_daq_config(daq_config_file)        // logging only
        .set_hycal_map(std::move(hycal_map_override))
        .set_gem_map(std::move(gem_map_override))
        .set_gem_pedestal("")         // empty falls back to RunConfig default
        .set_run_number_from_evio(input_evio)
        .set_log_stream(&std::cerr)
        .build();

    daq_cfg_ = std::move(pipeline.daq_cfg);
    auto &gem_sys = pipeline.gem;
    const auto roc_to_crate = daq_cfg_.roc_crate_map(true);

    if (zerosup_override >= 0.f) {
        gem_sys.SetZeroSupThreshold(zerosup_override);
        std::cerr << "Zero-sup : " << zerosup_override << " sigma (override)\n";
    }

    evc::EvChannel ch;
    ch.SetConfig(daq_cfg_);

    if (ch.OpenAuto(input_evio) != evc::status::success) {
        std::cerr << "Replay: cannot open " << input_evio << "\n";
        return false;
    }

    TFile *outfile = TFile::Open(output_root.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Replay: cannot create " << output_root << "\n";
        return false;
    }

    TTree *tree = new TTree("events", "PRad2 replay data");
    auto ev = std::make_unique<EventVars>();
    prad2::SetRawWriteBranches(tree, *ev, write_peaks, Ecalib, noWaveform);

    auto side = std::make_unique<SideTrees>();

    auto event = std::make_unique<fdec::EventData>();
    auto ssp_evt = std::make_unique<ssp::SspEventData>();
    fdec::WaveAnalyzer ana(daq_cfg_.wave_cfg);
    fdec::PulseTemplateStore template_store;
    template_store.LoadFromConfig(daq_cfg_.wave_cfg, db_dir);
    ana.SetTemplateStore(&template_store);
    fdec::WaveResult wres;
    // Firmware-mode emulator (FADC250 Modes 1/2/3).  Configured from the
    // optional "fadc250_waveform.firmware" block in daq_config.json — defaults
    // are safe for DAQ signal studies but should be overridden to match the
    // actual run's TET/NSB/NSA/MAX_PULSES if comparing to firmware output.
    fdec::Fadc250FwAnalyzer fw_ana(daq_cfg_.fadc250_fw);
    fdec::DaqWaveResult dwres;
    int total = 0;

    int run_num = get_run_int(input_evio);
    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(gRunConfig, run_num);

    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;

        if (side->record(ch)) continue;

        // Raw banks of this read group, shared by all its sub-events: the
        // 0xE10C SSP trigger bank and every VTP and TDC bank.
        const auto ssp_raw = snapshotFirstBank(ch, 0xE10C);
        const auto vtp_banks = snapshotBanks(ch, ch.GetVtpTags());
        const auto tdc_banks = snapshotBanks(ch, ch.GetTdcTags());

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            event->clear();
            ssp_evt->clear();
            if (!ch.DecodeEvent(ie, *event, ssp_evt.get())) continue;
            if (max_events > 0 && total >= max_events) break;

            if (ie == 0) side->recordScalers(ch, event->info, daq_cfg_.dsc_scaler);

            ev->clear();
            ev->event_num    = event->info.event_number;
            ev->trigger_type = event->info.trigger_type;
            ev->trigger_bits      = event->info.trigger_bits;
            ev->timestamp    = event->info.timestamp;
            ev->ssp_raw      = ssp_raw;
            ev->vtp_roc_tags = vtp_banks.roc_tags;
            ev->vtp_nwords   = vtp_banks.nwords;
            ev->vtp_words    = vtp_banks.words;
            ev->tdc_roc_tags = tdc_banks.roc_tags;
            ev->tdc_nwords   = tdc_banks.nwords;
            ev->tdc_words    = tdc_banks.words;

            if (Ecalib) {
                bool is_sum = (ev->trigger_bits & prad2::TBIT_sum) != 0;
                if (!is_sum) continue;
            }

            // Per-event gain correction (time-series lookup by event number).
            const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(ev->event_num));

            // Decode FADC250 data — single pass over all channels (HyCal +
            // Veto + LMS).  Type dispatch comes from hycal_map.json's "t"
            // field, not module-name prefix; module_type[nch] records the
            // category.
            int nch = 0;
            for (int r = 0; r < event->nrocs; ++r) {
                auto &roc = event->rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                int crate = (cit == roc_to_crate.end()) ? (int)roc.tag : cit->second;
                fdec::ForEachChannel(roc, [&](int s, int c, const fdec::ChannelData &cd) {
                    if (nch >= prad2::kMaxChannels) return;

                    int  mod_id   = moduleID(crate, s, c);
                    auto mod_type = moduleType(crate, s, c);
                    // Drop channels with no DAQ-map / module-info entry —
                    // we have no way to interpret them downstream.
                    if (mod_id < 0) return;

                    ev->module_id[nch]   = static_cast<uint16_t>(mod_id);
                    ev->module_type[nch] = static_cast<uint8_t>(mod_type);
                    ev->nsamples[nch]    = static_cast<uint8_t>(cd.nsamples);
                    for (int i = 0; i < cd.nsamples && i < fdec::MAX_SAMPLES; ++i)
                        ev->samples[nch][i] = cd.samples[i];

                    // Gain correction applies to HyCal modules only; Veto /
                    // LMS get 1.
                    ev->gain_factor[nch] = gain_corr.ModuleGain(mod_id);

                    if (write_peaks) {
                        // Soft analyzer drives both peaks AND the
                        // pedestal estimate that the firmware analyzer
                        // consumes — only run it when its output is
                        // being written.
                        ana.SetChannelKey(roc.tag, s, c);
                        ana.Analyze(cd.samples, cd.nsamples, wres);
                        ev->ped_mean[nch]    = wres.ped.mean;
                        ev->ped_rms[nch]     = wres.ped.rms;
                        ev->ped_nused[nch]   = wres.ped.nused;
                        ev->ped_quality[nch] = wres.ped.quality;
                        ev->ped_slope[nch]   = wres.ped.slope;
                        ev->npeaks[nch]   = static_cast<uint8_t>(wres.npeaks);
                        for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; p++) {
                            ev->peak_height[nch][p]   = wres.peaks[p].height;
                            ev->peak_time[nch][p]     = wres.peaks[p].time;
                            ev->peak_integral[nch][p] = wres.peaks[p].integral;
                            ev->peak_quality[nch][p]  = wres.peaks[p].quality;
                        }
                        fw_ana.Analyze(cd.samples, cd.nsamples, wres.ped.mean, dwres);
                        ev->daq_npeaks[nch] = static_cast<uint8_t>(dwres.npeaks);
                        for (int p = 0; p < dwres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                            const auto &dp = dwres.peaks[p];
                            ev->daq_peak_vp[nch][p]       = dp.vpeak;
                            ev->daq_peak_integral[nch][p] = dp.integral;
                            ev->daq_peak_time[nch][p]     = dp.time_ns;
                            ev->daq_peak_cross[nch][p]    = dp.cross_sample;
                            ev->daq_peak_pos[nch][p]      = dp.peak_sample;
                            ev->daq_peak_coarse[nch][p]   = dp.coarse;
                            ev->daq_peak_fine[nch][p]     = dp.fine;
                            ev->daq_peak_quality[nch][p]  = dp.quality;
                        }
                    }
                    nch++;
                }, 16);
            }
            ev->nch = nch;

            // decode and process GEM SSP data via GemSystem
            if (!Ecalib) {
                gem_sys.Clear();
                gem_sys.ProcessEvent(*ssp_evt);
                int gem_ch = 0;
                for (int d = 0; d < gem_sys.GetNDetectors(); ++d) {
                    for (int p = 0; p < 2; ++p) {
                        for (const auto &h : gem_sys.GetPlaneHits(d, p)) {
                            if (gem_ch >= prad2::kMaxGemStrips) break;
                            ev->gem_det[gem_ch]   = static_cast<uint8_t>(d);
                            ev->gem_plane[gem_ch] = static_cast<uint8_t>(p);
                            ev->gem_strip[gem_ch] = h.strip;
                            ev->gem_charge[gem_ch] = h.charge;
                            ev->gem_max_tb[gem_ch] = h.max_timebin;
                            ev->gem_pos[gem_ch]   = h.position;
                            ev->gem_xtalk[gem_ch] = h.cross_talk ? 1u : 0u;
                            for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t)
                                ev->gem_ts_adc[gem_ch][t] = t < (int)h.ts_adc.size() ? h.ts_adc[t] : 0.f;
                            ++gem_ch;
                        }
                    }
                }
                ev->gem_nch = gem_ch;
            } // end of if (!Ecalib)
            tree->Fill();
            total++;

            if (total % 10000 == 0)
                std::cerr << "\rReplay: " << total << " events processed" << std::flush;
        }
        if (max_events > 0 && total >= max_events) break;
    }

    std::cerr << "\rReplay: " << total << " events written to " << output_root << "\n";
    outfile->cd();
    tree->Write();
    side->write();
    delete outfile;
    return true;
}

bool Replay::ProcessWithRecon(const std::string &input_evio, const std::string &output_root, RunConfig &gRunConfig,
                                const std::string &db_dir, const std::string &recon_config_file,
                                const std::string &daq_config_file, const std::string &gem_ped_file,
                                const float zerosup_override, bool prad1, bool x17, bool x17_blind, bool random, bool gem_hit)
{
    // Similar to Process(), but with HyCal reconstruction and GEM hit reconstruction
    // before filling the ROOT tree.
    // The main differences are:
    // - After decoding, we run the HyCal clusterer to reconstruct clusters and hits.
    // - We also run the GemSystem reconstruction to get GEM hits.
    // - We fill a different TTree with reconstructed quantities instead of raw data.

    // Detectors: PRad-II flows through PipelineBuilder so the wiring stays in
    // one place (see prad2det/include/PipelineBuilder.h).  PRad-1 keeps its
    // hand-wired path because the builder is PRad-II-shaped (no GEM, different
    // hycal map, ADC1881M pedestals).
    fdec::HyCalSystem                 hycal;
    gem::GemSystem                    gem_sys;
    fdec::ClusterConfig               cluster_cfg;
    prad2::HyCalTimeCuts              hc_time_cuts;
    prad2::HyCalRfOffsets             hc_rf_offsets;
    DetectorTransform                 hycal_transform;
    std::array<DetectorTransform, 4>  gem_transforms;
    std::unordered_map<uint32_t, int> roc_to_crate;
    int                               match_method = 1;

    if (prad1) {
        // Legacy PRad-1 setup — no GEM, ADC1881M pedestals.
        std::string hycal_map_file = db_dir + "/prad1/prad_hycal_map.json";
        hycal.Init(hycal_map_file);
        evc::load_pedestals(db_dir + "/prad1/adc1881m_pedestals.json", daq_cfg_);

        std::string calib_file = db_dir + "/" + gRunConfig.energy_calib_file;
        int nmatched = hycal.LoadCalibration(calib_file);
        if (nmatched >= 0)
            std::cerr << "Calibration: " << calib_file << " (" << nmatched << " modules)\n";

        // Every roc_tags entry of the loaded DAQ config (no data-ROC filter).
        if (!daq_config_file.empty()) roc_to_crate = daq_cfg_.roc_crate_map();

        // PRad-1 transforms come from the externally-loaded gRunConfig (the
        // builder owns this for PRad-II).
        auto t = analysis::BuildLabTransforms(gRunConfig);
        hycal_transform = t.hycal;
        gem_transforms  = t.gem;

        // PRad-1 has no per-module time-cut file; build a uniform table so
        // the per-event loop uses the same call as PRad-II.
        hc_time_cuts = prad2::LoadHyCalTimeCuts(
            "", hycal,
            gRunConfig.hc_time_win_lo, gRunConfig.hc_time_win_hi);
        // PRad-1 had no RF readout — keep the offset table as uniform 0 so
        // the per-event apply() call is a no-op (and cl_dt_rf branch stays
        // NaN since rf_n_a == 0 for every event).
        hc_rf_offsets = prad2::LoadHyCalRfOffsets("", hycal, 0.f);
    } else {
        std::string hycal_map_override = daq_cfg_.hycal_map_file;
        std::string gem_map_override   = daq_cfg_.gem_map_file;

        prad2::Pipeline pipeline = prad2::PipelineBuilder()
            .set_recon_config(recon_config_file)
            .set_database_dir(db_dir)
            .set_loaded_daq_config(std::move(daq_cfg_))
            .set_daq_config(daq_config_file)        // logging only
            .set_hycal_map(std::move(hycal_map_override))
            .set_gem_map(std::move(gem_map_override))
            .set_gem_pedestal(gem_ped_file)         // empty falls back to RunConfig default
            .set_run_number_from_evio(input_evio)
            .set_log_stream(&std::cerr)
            .build();

        daq_cfg_         = std::move(pipeline.daq_cfg);
        hycal            = std::move(pipeline.hycal);
        gem_sys          = std::move(pipeline.gem);
        cluster_cfg      = pipeline.hycal_cluster_cfg;
        hc_time_cuts     = std::move(pipeline.hycal_time_cuts);
        hc_rf_offsets    = std::move(pipeline.hycal_rf_offsets);
        hycal_transform  = pipeline.hycal_transform;
        gem_transforms   = pipeline.gem_transforms;
        match_method     = pipeline.match_method;
        roc_to_crate     = daq_cfg_.roc_crate_map(true);

        if (zerosup_override >= 0.f) {
            gem_sys.SetZeroSupThreshold(zerosup_override);
            std::cerr << "Zero-sup : " << zerosup_override << " sigma (override)\n";
        }
    }

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);
    gem::GemCluster      gem_clusterer;
    MatchingTools        matching(match_method);
    matching.Configure(gRunConfig);
    evc::EvChannel ch;
    ch.SetConfig(daq_cfg_);

    if (ch.OpenAuto(input_evio) != evc::status::success) {
        std::cerr << "Replay: cannot open " << input_evio << "\n";
        return false;
    }

    TFile *outfile = TFile::Open(output_root.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Replay: cannot create " << output_root << "\n";
        return false;
    }

    TTree *tree = new TTree("recon", "PRad2 replay reconstruction");
    auto ev = std::make_unique<EventVars_Recon>();
    prad2::SetReconWriteBranches(tree, *ev, x17, gem_hit);

    auto side = std::make_unique<SideTrees>();

    auto event = std::make_unique<fdec::EventData>();
    auto ssp_evt = std::make_unique<ssp::SspEventData>();
    fdec::WaveAnalyzer ana(daq_cfg_.wave_cfg);
    fdec::PulseTemplateStore template_store;
    template_store.LoadFromConfig(daq_cfg_.wave_cfg, db_dir);
    ana.SetTemplateStore(&template_store);
    fdec::WaveResult wres;

    int run_num = get_run_int(input_evio);
    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(gRunConfig, run_num);

    // Per-detector lab transforms — set up by either branch of the detector
    // wiring above (PipelineBuilder for PRad-II, BuildLabTransforms for PRad-1).
    const auto &hc_xform = hycal_transform;
    const auto &g_xform  = gem_transforms;
    
    int total = 0;
    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;

        if (side->record(ch)) continue;

        // Raw banks of this read group (see Process()); the TDC banks give
        // the per-cluster RF Δt.
        const auto ssp_raw = snapshotFirstBank(ch, 0xE10C);
        const auto vtp_banks = snapshotBanks(ch, ch.GetVtpTags());
        const auto tdc_banks = snapshotBanks(ch, ch.GetTdcTags());

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            event->clear();
            ssp_evt->clear();
            clusterer.Clear();
            if (!ch.DecodeEvent(ie, *event, ssp_evt.get())) continue;

            if (ie == 0) side->recordScalers(ch, event->info, daq_cfg_.dsc_scaler);

            ev->clear();
            ev->event_num    = event->info.event_number;
            ev->trigger_type = event->info.trigger_type;
            ev->trigger_bits = event->info.trigger_bits;
            ev->timestamp    = event->info.timestamp;
            ev->ssp_raw      = ssp_raw;
            ev->vtp_roc_tags = vtp_banks.roc_tags;
            ev->vtp_nwords   = vtp_banks.nwords;
            ev->vtp_words    = vtp_banks.words;

            fillVtpClusters(*ev, ch.Vtp());

            const auto rf = fillRf(*ev, tdc_banks.roc_tags, tdc_banks.nwords, tdc_banks.words);

            const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(ev->event_num));

            const auto trig = selectTrigger(ev->trigger_bits, ev->event_num,
                                            random, x17, x17_blind);
            if (!trig.keep) continue;

            // decode FADC250 and reconstruct HyCal data
            int veto_nch = 0;
            int lms_nch = 0;
            int nch = 0;
            for (int r = 0; r < event->nrocs; ++r) {
                auto &roc = event->rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                if (cit == roc_to_crate.end()) continue;
                int crate = cit->second;
                // All 64 channels per slot, not the FADC250's 16: PRad-1 ADC1881M.
                fdec::ForEachChannel(roc, [&](int s, int c, const fdec::ChannelData &cd) {
                    std::string mod_name = moduleName(crate, s, c);
                    if(mod_name.empty()) return;
                    const auto *mod = hycal.module_by_daq(crate, s, c);

                    if(trig.lms || trig.alpha) {
                        if(mod_name[0] == 'L'){
                            if(mod_name.length() != 4) return;
                            if(lms_nch >= 4) return; // guard against overflow
                            if(mod_name[3] == 'P') ev->lms_id[lms_nch] = 0;
                            else ev->lms_id[lms_nch] = mod_name[3] - '0';
                            ana.SetChannelKey(roc.tag, s, c);
                            ana.Analyze(cd.samples, cd.nsamples, wres);
                            ev->lms_npeaks[lms_nch] = wres.npeaks;
                            if(wres.npeaks <= 0) return;
                            for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                ev->lms_peak_height[lms_nch][p] = wres.peaks[p].height;
                                ev->lms_peak_integral[lms_nch][p] = wres.peaks[p].integral;
                                ev->lms_peak_time[lms_nch][p] = wres.peaks[p].time;
                            }
                            lms_nch++;
                        }
                        else return;
                    }

                    if((trig.sum || trig.cl3) && !trig.lms) {
                        if(mod_name[0] == 'V'){
                            if(mod_name.length() != 2) return;
                            if(veto_nch >= 4) return; // guard against overflow
                            // "V1".."V4" → 1..4
                            ev->veto_id[veto_nch] = mod_name[1] - '0';
                            ana.SetChannelKey(roc.tag, s, c);
                            ana.Analyze(cd.samples, cd.nsamples, wres);
                            ev->veto_npeaks[veto_nch] = wres.npeaks;
                            if(wres.npeaks <= 0) return;
                            for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                ev->veto_peak_height[veto_nch][p] = wres.peaks[p].height;
                                ev->veto_peak_integral[veto_nch][p] = wres.peaks[p].integral;
                                ev->veto_peak_time[veto_nch][p] = wres.peaks[p].time;
                            }
                            veto_nch++;
                        }
                        else{
                            if (!mod || !mod->is_hycal()) return;
                            const float gain = gain_corr.ModuleGain(mod->id);

                            if (prad1 == true) {
                                float adc = cd.samples[0] * 0.543f * gain; //0.543 for prad1 run1308, correct to 1.1GeV
                                float energy = static_cast<float>(mod->energize(adc));
                                clusterer.AddHit(mod->index, energy, 0.f);
                                ev->total_energy += energy;
                                nch++;
                                return;
                            }

                            ana.SetChannelKey(roc.tag, s, c);
                            ana.Analyze(cd.samples, cd.nsamples, wres, mod->time_offset);
                            if (wres.npeaks <= 0) return;

                            auto hc_win = hc_time_cuts.at(mod->index);
                            if (random) {hc_win.lo = 0; hc_win.hi = 400;}
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
                                // module hit.
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
                                if (bestIdx < 0) return;
                                float adc = wres.peaks[bestIdx].integral * gain;
                                float energy = static_cast<float>(mod->energize(adc));
                                clusterer.AddHit(mod->index, energy, wres.peaks[bestIdx].time);
                                ev->total_energy += energy;
                                nch++;
                            }
                        }
                    }
                });
            }
            ev->veto_nch = veto_nch;
            ev->lms_nch = lms_nch;
            if(nch > 1000) continue; // too many hits, likely noise, skip the event

            clusterer.FormClusters();
            std::vector<fdec::ClusterHit> hits;
            clusterer.ReconstructHits(hits);
            fillClusters(*ev, hits, hc_xform, gRunConfig.hycal_z, rf, hycal, hc_rf_offsets);

            //decode GEM data and reconstruct GEM hits
        if(!prad1 && gem_sys.GetNDetectors() > 0){
            gem_sys.Clear();
            gem_sys.ProcessEvent(*ssp_evt);
            gem_sys.Reconstruct(gem_clusterer);
            fillGemHits(*ev, gem_sys.GetAllHits(), g_xform);
            // per-cluster QA block: det 0 X, det 0 Y, det 1 X, ...
            for (int d = 0; d < gem_sys.GetNDetectors(); ++d)
                for (int p = 0; p < 2; ++p)
                    appendGemClusters(*ev, d, p, gem_sys.GetPlaneClusters(d, p));

            MatchReconEvent(*ev, matching);

        }
            tree->Fill();
            total++;
            if (total % 10000 == 0)
                std::cerr << "\rReplay: " << total << " events processed" << std::flush;
        }
    }
    std::cerr << "\rReplay: " << total << " events reconstructed -> " << output_root << "\n";
    outfile->cd();
    tree->Write();
    side->write();
    delete outfile;

    return true;
}

bool Replay::ProcessRaw2Recon(const std::string &input_raw, const std::string &output_root, RunConfig &gRunConfig,
                                const std::string &db_dir, const std::string &recon_config_file,
                                const std::string &daq_config_file, const std::string &gem_ped_file,
                                bool x17, bool x17_blind, bool random, bool gem_hit)
{
    // Similar to ProcessWithRecon(), with HyCal reconstruction and GEM hit reconstruction
    // before filling the ROOT tree. But unlike ProcessWithRecon(), it starts from raw root files 
    // input rather than EVIO events.
    // The main differences are:
    // - Read decoded data, run the HyCal clusterer to reconstruct clusters and hits.
    // - We also run the GemSystem reconstruction to get GEM hits, from strip-level data.
    // - We fill a different TTree with reconstructed quantities instead of raw data.

    int run_num = get_run_int(input_raw);

    std::string hycal_map_override = daq_cfg_.hycal_map_file;
    std::string gem_map_override   = daq_cfg_.gem_map_file;

    prad2::Pipeline pipeline = prad2::PipelineBuilder()
        .set_recon_config(recon_config_file)
        .set_database_dir(db_dir)
        .set_loaded_daq_config(std::move(daq_cfg_))
        .set_daq_config(daq_config_file)        // logging only
        .set_hycal_map(std::move(hycal_map_override))
        .set_gem_map(std::move(gem_map_override))
        .set_gem_pedestal(gem_ped_file)         // empty falls back to RunConfig default
        .set_run_number(run_num)
        .set_log_stream(&std::cerr)
        .build();

    daq_cfg_ = std::move(pipeline.daq_cfg);
    const auto &hycal         = pipeline.hycal;
    const auto &gem_sys       = pipeline.gem;
    const auto &cluster_cfg   = pipeline.hycal_cluster_cfg;
    const auto &hc_time_cuts  = pipeline.hycal_time_cuts;
    const auto &hc_rf_offsets = pipeline.hycal_rf_offsets;
    const auto &hc_xform      = pipeline.hycal_transform;
    const auto &g_xform       = pipeline.gem_transforms;

    fdec::HyCalCluster   clusterer(hycal);
    clusterer.SetConfig(cluster_cfg);
    gem::GemCluster      gem_clusterer;
    MatchingTools        matching(pipeline.match_method);
    matching.Configure(gRunConfig);

    TFile *outfile = TFile::Open(output_root.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Replay: cannot create " << output_root << "\n";
        return false;
    }

    TFile *infile = TFile::Open(input_raw.c_str(), "READ");
    if (!infile || !infile->IsOpen()) {
        std::cerr << "Replay: cannot open " << input_raw << "\n";
        return false;
    }
    TTree *tree_in = dynamic_cast<TTree *>(infile->Get("events"));
    if (!tree_in) {
        std::cerr << "Replay: input raw file has no 'events' tree\n";
        return false;
    }
    if (tree_in->GetBranch("gem.nch") &&
        (!tree_in->GetBranch("gem.det") ||
         !tree_in->GetBranch("gem.plane") ||
         !tree_in->GetBranch("gem.charge") ||
         !tree_in->GetBranch("gem.max_tb") ||
         !tree_in->GetBranch("gem.pos") ||
         !tree_in->GetBranch("gem.xtalk") ||
         !tree_in->GetBranch("gem.ts_adc"))) {
        std::cerr << "Replay: input raw file uses the legacy GEM schema; "
                  << "regenerate it with the current replay_rawdata before "
                  << "running raw-to-recon\n";
        return false;
    }

    // A raw replay tree carries waveform-analysis diagnostics that this
    // reconstruction path never consumes.  Disable everything first so
    // GetEntry() only decompresses data used below or copied to recon output.
    const bool has_waveform = tree_in->GetBranch("hycal.samples") != nullptr;
    tree_in->SetBranchStatus("*", 0);
    auto enable_branch = [tree_in](const char *name) {
        if (tree_in->GetBranch(name)) tree_in->SetBranchStatus(name, 1);
    };
    for (const char *name : {
             "event_num", "trigger_type", "trigger_bits", "timestamp",
             "hycal.nch", "hycal.module_id", "hycal.module_type",
             "hycal.gain_factor",
             "gem.nch", "gem.det", "gem.plane", "gem.strip",
             "gem.charge", "gem.max_tb", "gem.pos", "gem.xtalk",
             "gem.ts_adc",
             "ssp_raw",
             "vtp_roc_tags", "vtp_nwords", "vtp_words",
             "tdc_roc_tags", "tdc_nwords", "tdc_words"}) {
        enable_branch(name);
    }
    if (has_waveform) {
        enable_branch("hycal.nsamples");
        enable_branch("hycal.samples");
    } else {
        enable_branch("hycal.npeaks");
        enable_branch("hycal.peak_height");
        enable_branch("hycal.peak_time");
        enable_branch("hycal.peak_integral");
    }

    auto in = std::make_unique<EventVars>();
    prad2::SetRawReadBranches(tree_in, *in);
    prad2::RawVectorBindings raw_vecs;
    prad2::BindRawVectorBranches(tree_in, *in, raw_vecs);

    // All newly-created trees must belong to the output file, not the input
    // file most recently opened above.
    outfile->cd();

    TTree *tree = new TTree("recon", "PRad2 replay reconstruction");
    auto ev = std::make_unique<EventVars_Recon>();
    prad2::SetReconWriteBranches(tree, *ev, x17, gem_hit);

    // Side trees already exist in replay_raw output.  Copy them verbatim so
    // their event-number join semantics survive raw-to-recon conversion.
    auto copy_side_tree = [infile, outfile](const char *name) -> TTree * {
        auto *source = dynamic_cast<TTree *>(infile->Get(name));
        if (!source) {
            std::cerr << "Replay: input raw file has no '" << name
                      << "' tree; skipping it\n";
            return nullptr;
        }
        outfile->cd();
        return source->CloneTree(-1, "fast");
    };
    TTree *scalers_tree = copy_side_tree("scalers");
    TTree *epics_tree   = copy_side_tree("epics");
    TTree *runinfo_tree = copy_side_tree("runinfo");

    fdec::WaveAnalyzer ana(daq_cfg_.wave_cfg);
    fdec::PulseTemplateStore template_store;
    template_store.LoadFromConfig(daq_cfg_.wave_cfg, db_dir);
    ana.SetTemplateStore(&template_store);
    fdec::WaveResult wres;

    auto gain_corr_ts = prad2::LoadGainCorrTimeSeries(gRunConfig, run_num);

    int total = 0;
    long long nentries = tree_in->GetEntries();
    for (long long i = 0; i < nentries; ++i) {
        in->clear_banks();
        tree_in->GetEntry(i);
        if (i % 1000 == 0) std::cout << "Processed " << i << " / " << nentries << " entries.\r" << std::flush;

        ev->clear();
        ev->event_num    = in->event_num;
        ev->trigger_type = in->trigger_type;
        ev->trigger_bits = in->trigger_bits;
        ev->timestamp    = in->timestamp;
        ev->ssp_raw      = in->ssp_raw;
        ev->vtp_roc_tags = in->vtp_roc_tags;
        ev->vtp_nwords   = in->vtp_nwords;
        ev->vtp_words    = in->vtp_words;

        clusterer.Clear();

        // Re-decode PRAD_CLUSTER records from the flat VTP bank snapshot.
        vtp::VtpEventData vtp_event;
        if (!vtp::VtpDecoder::DecodeReplay(in->vtp_roc_tags, in->vtp_nwords,
                                           in->vtp_words, vtp_event))
            std::cerr << "Replay: malformed VTP snapshot in event "
                      << ev->event_num << "\n";
        fillVtpClusters(*ev, vtp_event);

        const auto rf = fillRf(*ev, in->tdc_roc_tags, in->tdc_nwords, in->tdc_words);

        const auto &gain_corr = gain_corr_ts.GetCorr(static_cast<int>(ev->event_num));

        const auto trig = selectTrigger(ev->trigger_bits, ev->event_num,
                                        random, x17, x17_blind);
        if (!trig.keep) continue;

        int lms_nch = 0;
        if (trig.lms || trig.alpha) {
            for (int j = 0; j < in->nch; ++j) {
                if (in->module_type[j] != prad2::MOD_LMS) continue;
                if(lms_nch >= 4) continue; // guard against overflow
                ev->lms_id[lms_nch] = in->module_id[j] - prad2::kLmsIdBase;
                // if has waveform data, reanalyze the waveforms,
                // or just use the existing peak information.
                if (has_waveform)
                {
                    int crate, slot, ch;
                    std::tie(crate, slot, ch) = moduleLocation(in->module_id[j]);
                    ana.SetChannelKey(crate, slot, ch);
                    ana.Analyze(in->samples[j], in->nsamples[j], wres);
                    ev->lms_npeaks[lms_nch] = wres.npeaks;
                    if(wres.npeaks <= 0) continue;
                    for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                        ev->lms_peak_height[lms_nch][p] = wres.peaks[p].height;
                        ev->lms_peak_integral[lms_nch][p] = wres.peaks[p].integral;
                        ev->lms_peak_time[lms_nch][p] = wres.peaks[p].time;
                    }
                }
                else 
                {
                    ev->lms_npeaks[lms_nch] = in->npeaks[j];
                    if(in->npeaks[j] <= 0) continue;
                    for (int p = 0; p < in->npeaks[j] && p < fdec::MAX_PEAKS; ++p) {
                        ev->lms_peak_height[lms_nch][p] = in->peak_height[j][p];
                        ev->lms_peak_integral[lms_nch][p] = in->peak_integral[j][p];
                        ev->lms_peak_time[lms_nch][p] = in->peak_time[j][p];
                    }
                }
                lms_nch++;
            }
        }

        int veto_nch = 0;
        int nch = 0;
        if ((trig.sum || trig.cl3) && !trig.lms) {
            for (int j = 0; j < in->nch; ++j) {
                if (in->module_type[j] == prad2::MOD_VETO) {
                    if(veto_nch >= 4) continue; // guard against overflow
                    ev->veto_id[veto_nch] = in->module_id[j] - prad2::kVetoIdBase;
                    if (has_waveform)
                    {
                        int crate, slot, ch;
                        std::tie(crate, slot, ch) = moduleLocation(in->module_id[j]);
                        ana.SetChannelKey(crate, slot, ch);
                        ana.Analyze(in->samples[j], in->nsamples[j], wres);
                        ev->veto_npeaks[veto_nch] = wres.npeaks;
                        if(wres.npeaks <= 0) continue;
                        for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                            ev->veto_peak_height[veto_nch][p] = wres.peaks[p].height;
                            ev->veto_peak_integral[veto_nch][p] = wres.peaks[p].integral;
                            ev->veto_peak_time[veto_nch][p] = wres.peaks[p].time;
                        }
                    }
                    else
                    {
                        ev->veto_npeaks[veto_nch] = in->npeaks[j];
                        for (int p = 0; p < in->npeaks[j] && p < fdec::MAX_PEAKS; ++p) {
                            ev->veto_peak_height[veto_nch][p] = in->peak_height[j][p];
                            ev->veto_peak_integral[veto_nch][p] = in->peak_integral[j][p];
                            ev->veto_peak_time[veto_nch][p] = in->peak_time[j][p];
                        }
                    }
                    veto_nch++;
                }
                else {
                    const auto *mod = hycal.module_by_id(in->module_id[j]);
                    if (!mod || !mod->is_pwo4()) continue;
                    int crate, slot, ch;
                    std::tie(crate, slot, ch) = moduleLocation(in->module_id[j]);

                    float time_offset = mod->time_offset;

                    float gain = gain_corr.ModuleGain(mod->id);
                    if (gain <= 0.f || gain == 1.0f) gain = in->gain_factor[j];
                    
                    if (has_waveform) {
                        ana.SetChannelKey(crate, slot, ch);
                        ana.Analyze(in->samples[j], in->nsamples[j], wres, time_offset);
                        if (wres.npeaks <= 0) continue;
                    }

                    auto hc_win = hc_time_cuts.at(mod->index);
                    if (random) { hc_win.lo = 0; hc_win.hi = 400; }
                    if (cluster_cfg.seed_time_window > 0.f) {
                        // Multi-pulse mode, as in ProcessWithRecon().
                        if (has_waveform)
                        {
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
                        }
                        else
                        {
                            for (int p = 0; p < in->npeaks[j] && p < fdec::MAX_PEAKS; ++p) {
                                float peak_time = in->peak_time[j][p] - time_offset;
                                if (peak_time <= hc_win.lo) continue;
                                if (peak_time >= hc_win.hi) continue;
                                float adc = in->peak_integral[j][p] * gain;
                                float energy = static_cast<float>(mod->energize(adc));
                                clusterer.AddHit(mod->index, energy, peak_time);
                                ev->total_energy += energy;
                                nch++;
                            }
                        }
                    } else {
                        // Legacy: pick the largest in-window peak as the single
                        // module hit.
                        int bestIdx = -1;
                        float bestHeight = -1.f;
                        if (has_waveform)
                        {
                            for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; ++p) {
                                const auto &pk = wres.peaks[p];
                                if (pk.time > hc_win.lo &&
                                    pk.time < hc_win.hi &&
                                    pk.height > bestHeight) {
                                    bestHeight = pk.height;
                                    bestIdx = p;
                                }
                            }
                        }
                        else
                        {
                            for (int p = 0; p < in->npeaks[j] && p < fdec::MAX_PEAKS; ++p) {
                                float peak_time = in->peak_time[j][p] - time_offset;
                                if (peak_time <= hc_win.lo) continue;
                                if (peak_time >= hc_win.hi) continue;
                                if (in->peak_integral[j][p] > bestHeight) {
                                    bestHeight = in->peak_integral[j][p];
                                    bestIdx = p;
                                }
                            }
                        }
                        if (bestIdx < 0) continue;
                        float adc = (has_waveform ? wres.peaks[bestIdx].integral : in->peak_integral[j][bestIdx]) * gain;
                        float energy = static_cast<float>(mod->energize(adc));
                        clusterer.AddHit(mod->index, energy, (has_waveform ? wres.peaks[bestIdx].time : in->peak_time[j][bestIdx] - time_offset));
                        ev->total_energy += energy;
                        nch++;
                    }
                }
            }
        }
        ev->veto_nch = veto_nch;
        ev->lms_nch = lms_nch;
        if(nch > 1000) continue; // too many hits, likely noise, skip the event

        clusterer.FormClusters();
        std::vector<fdec::ClusterHit> hits;
        clusterer.ReconstructHits(hits);
        fillClusters(*ev, hits, hc_xform, gRunConfig.hycal_z, rf, hycal, hc_rf_offsets);
        // Reconstruct GEM hits from the strip-level data saved on the raw
        // replay tree.  This starts after pedestal/common-mode/ZS, so it
        // resumes at the same clustering and X/Y matching stage used by
        // ProcessWithRecon().
        if (gem_sys.GetNDetectors() > 0) {
            std::vector<gem::GEMHit> all_gem_hits;
            std::vector<std::array<std::vector<gem::StripCluster>, 2>> plane_clusters;
            ReconstructGemStrips(*in, gem_sys, gem_clusterer, all_gem_hits, &plane_clusters);
            fillGemHits(*ev, all_gem_hits, g_xform);
            // per-cluster QA block, same order as ProcessWithRecon
            for (int d = 0; d < gem_sys.GetNDetectors(); ++d)
                for (int p = 0; p < 2; ++p)
                    appendGemClusters(*ev, d, p, plane_clusters[d][p]);

            MatchReconEvent(*ev, matching);
        }
        tree->Fill();
        ++total;
    }

    std::cerr << "\rReplay: " << total << " events reconstructed -> "
              << output_root << "\n";
    outfile->cd();
    tree->Write();
    if (scalers_tree) scalers_tree->Write();
    if (epics_tree) epics_tree->Write();
    if (runinfo_tree) runinfo_tree->Write();
    delete infile;
    delete outfile;
    return true;
}

bool Replay::Process_LMSgainFactor(const std::string &input_evio, const std::string &output_root,
                     const std::string &db_dir, const std::string &daq_config_file)
{
    // ROC tag → crate over every roc_tags entry of the loaded DAQ config
    std::unordered_map<uint32_t, int> roc_to_crate;
    if (!daq_config_file.empty()) {
        std::cout << "Loading DAQ config from " << daq_config_file << "\n";
        roc_to_crate = daq_cfg_.roc_crate_map();
    }
    else {
        std::cerr << "No DAQ config file provided, ROC tag to crate mapping will be unavailable.\n";
    }

    evc::EvChannel ch;
    ch.SetConfig(daq_cfg_);

    auto ev = std::make_unique<LMSEventVars>();

    if (ch.OpenAuto(input_evio) != evc::status::success) {
        std::cerr << "Replay: cannot open " << input_evio << "\n";
        return false;
    }

    TFile *outfile = TFile::Open(output_root.c_str(), "RECREATE");
    if (!outfile || !outfile->IsOpen()) {
        std::cerr << "Replay: cannot create " << output_root << "\n";
        return false;
    }

    TTree *tree = new TTree("lms_gain", "LMS gain factor calculation");
    prad2::SetLMSWriteBranches(tree, *ev);

    auto side = std::make_unique<SideTrees>();

    auto event = std::make_unique<fdec::EventData>();
    fdec::WaveAnalyzer ana(daq_cfg_.wave_cfg);
    fdec::PulseTemplateStore template_store;
    template_store.LoadFromConfig(daq_cfg_.wave_cfg, db_dir);
    ana.SetTemplateStore(&template_store);
    fdec::WaveResult wres;
    
    int total = 0;

    int run_num = get_run_int(input_evio);

    while (ch.Read() == evc::status::success) {
        if (!ch.Scan()) continue;

        if (side->record(ch)) continue;

        for (int ie = 0; ie < ch.GetNEvents(); ++ie) {
            event->clear();
            if (!ch.DecodeEvent(ie, *event, nullptr)) continue;

            if (ie == 0) side->recordScalers(ch, event->info, daq_cfg_.dsc_scaler);

            ev->clear();
            ev->event_num    = event->info.event_number;
            ev->trigger_type = event->info.trigger_type;
            ev->trigger_bits      = event->info.trigger_bits;
            ev->timestamp    = event->info.timestamp;

            bool trig_lms = (ev->trigger_bits & prad2::TBIT_lms) != 0;
            bool trig_alpha = (ev->trigger_bits & prad2::TBIT_alpha) != 0;
            if (!trig_lms && !trig_alpha) continue;

            int nch = 0;
            for (int r = 0; r < event->nrocs; ++r) {
                auto &roc = event->rocs[r];
                if (!roc.present) continue;
                auto cit = roc_to_crate.find(roc.tag);
                int crate = (cit == roc_to_crate.end()) ? (int)roc.tag : cit->second;
                fdec::ForEachChannel(roc, [&](int s, int c, const fdec::ChannelData &cd) {
                    if (nch >= prad2::kMaxChannels) return;

                    int  mod_id   = moduleID(crate, s, c);
                    auto mod_type = moduleType(crate, s, c);
                    if (mod_id < 0) return;

                    ev->module_id[nch] = mod_id;
                    ev->module_type[nch] = mod_type;

                    ana.SetChannelKey(roc.tag, s, c);
                    ana.Analyze(cd.samples, cd.nsamples, wres);
                    ev->npeaks[nch]   = static_cast<uint8_t>(wres.npeaks);
                    for (int p = 0; p < wres.npeaks && p < fdec::MAX_PEAKS; p++) {
                        ev->peak_height[nch][p]   = wres.peaks[p].height;
                        ev->peak_time[nch][p]     = wres.peaks[p].time;
                        ev->peak_integral[nch][p] = wres.peaks[p].integral;
                    }
                    nch++;
                }, 16);
            }
            ev->nch = nch;

            //Because the LMS and alpha trigger_bits can not believe, 
            // we need to use "nch" to seperate the LMS and alpha events
            bool is_lms = (trig_lms && ev->nch > 1000);
            bool is_alpha = (trig_alpha && ev->nch < 50);
            
            if(is_lms) ev->event_type = 0; // LMS event
            if(is_alpha) ev->event_type = 1; // alpha event

            if(is_lms || is_alpha) tree->Fill();

            total++;
        }
    }

    std::cerr << "\rReplay: " << total << " events written to " << output_root << "\n";
    outfile->cd();
    tree->Write();
    side->write();
    delete outfile;
    return true;
}

void FillPeaksFromWaveforms(prad2::RawEventData &ev, const fdec::HyCalSystem &hycal,
                            const fdec::WaveAnalyzer &ana, fdec::WaveResult &wres)
{
    for (int j = 0; j < ev.nch; ++j) {
        const auto *mod = hycal.module_by_id(ev.module_id[j]);
        if (!mod || !mod->is_pwo4()) continue;

        ana.Analyze(ev.samples[j], ev.nsamples[j], wres);
        ev.npeaks[j] = std::min(wres.npeaks, fdec::MAX_PEAKS);
        for (int p = 0; p < ev.npeaks[j]; ++p) {
            const auto &pk = wres.peaks[p];
            ev.peak_height[j][p]   = pk.height;
            ev.peak_time[j][p]     = pk.time;
            ev.peak_integral[j][p] = pk.integral;
        }
    }
}

void ReconstructGemStrips(const prad2::RawEventData &ev, const gem::GemSystem &gem_sys,
                          gem::GemCluster &clusterer, std::vector<gem::GEMHit> &hits,
                          std::vector<std::array<std::vector<gem::StripCluster>, 2>>
                              *plane_clusters)
{
    const int n_det = gem_sys.GetNDetectors();
    std::vector<std::array<std::vector<gem::StripHit>, 2>> plane_hits(n_det);
    const int n_strips = std::min(ev.gem_nch, prad2::kMaxGemStrips);
    for (int i = 0; i < n_strips; ++i) {
        const int det = ev.gem_det[i];
        const int plane = ev.gem_plane[i];
        if (det < 0 || det >= n_det || plane < 0 || plane > 1) continue;

        gem::StripHit hit;
        hit.strip       = ev.gem_strip[i];
        hit.charge      = ev.gem_charge[i];
        hit.max_timebin = ev.gem_max_tb[i];
        hit.position    = ev.gem_pos[i];
        hit.cross_talk  = ev.gem_xtalk[i] != 0;
        hit.ts_adc.assign(ev.gem_ts_adc[i], ev.gem_ts_adc[i] + ssp::SSP_TIME_SAMPLES);
        plane_hits[det][plane].push_back(std::move(hit));
    }

    hits.clear();
    if (plane_clusters) plane_clusters->assign(n_det, {});
    const auto &cfgs = gem_sys.GetReconConfigs();
    std::array<std::vector<gem::StripCluster>, 2> clusters;
    std::vector<gem::GEMHit> det_hits;
    for (int det = 0; det < n_det; ++det) {
        clusterer.SetConfig(cfgs[det]);
        for (int p = 0; p < 2; ++p)
            clusterer.FormClusters(plane_hits[det][p], clusters[p]);
        clusterer.CartesianReconstruct(clusters[0], clusters[1], det_hits, det);
        hits.insert(hits.end(), det_hits.begin(), det_hits.end());
        if (plane_clusters) (*plane_clusters)[det] = std::move(clusters);
    }
}

} // namespace analysis
