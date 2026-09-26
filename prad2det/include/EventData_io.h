#pragma once
//=============================================================================
// EventData_io.h — TTree branch I/O helpers for RawEventData / ReconEventData
//
// Header-only.  prad2det's library does NOT compile this file; it is included
// only by ROOT-aware consumers (analysis tools, the viewer's root data
// source, sim2replay).  The library still has no link-time ROOT dependency.
//
// All helpers are inline.  Writers set up the branch list selected by
// their flags.  Readers use TTree::GetBranch to skip any branch that's
// missing on disk, so they read older files that pre-date later branches.
//
// Usage:
//   #include "EventData_io.h"
//   prad2::SetRawWriteBranches(tree, ev, /*with_peaks=*/true);
//   ...
//   auto status = prad2::SetRawReadBranches(tree, ev);
//   if (status.has_peaks) { ... }
//
// Single source of truth for the replay tree schema.  Match this against
// `docs/REPLAYED_DATA.md` (or vice versa) when the layout changes.
//=============================================================================

#include "EventData.h"

#include "DscData.h"
#include "EpicsData.h"
#include "SyncData.h"
#include "Fadc250Data.h"
#include "DaqConfig.h"

#include <TTree.h>
#include <TString.h>   // for Form()

#include <cstring>
#include <type_traits>

namespace prad2 {

// Binds branch addresses, skipping branches the tree does not have.
// std::vector / std::string branches need a held pointer that outlives every
// GetEntry: bind(name, obj, holder) points `holder` at `obj` and binds it.
struct BranchBinder {
    TTree *tree;

    void operator()(const char *name, void *addr) const
    {
        if (tree->GetBranch(name)) tree->SetBranchAddress(name, addr);
    }

    template <class T>
    void operator()(const char *name, T &obj, T *&holder) const
    {
        holder = &obj;
        if (tree->GetBranch(name)) tree->SetBranchAddress(name, &holder);
    }
};

// ── Reader status ────────────────────────────────────────────────────────
struct RawReadStatus {
    bool has_peaks      = false;   // soft-analyzer: ped_mean/rms + peak_*
    bool has_daq_peaks  = false;   // firmware-mode: daq_peak_*
    bool has_gem        = false;   // gem.*
    bool has_ssp_raw    = false;   // ssp_raw vector
    bool has_vtp_raw    = false;   // vtp_roc_tags + vtp_nwords + vtp_words
    bool has_tdc_raw    = false;   // tdc_roc_tags + tdc_nwords + tdc_words
};

struct ReconReadStatus {
    bool has_match_num  = false;   // mHit_* quick-access arrays
    bool has_per_cl_match = false; // matchFlag / match_cl_idx / match_det_id / match_gem_*
    bool has_veto       = false;
    bool has_lms        = false;
    bool has_ssp_raw    = false;
    bool has_vtp_raw    = false;   // vtp_roc_tags + vtp_nwords + vtp_words
    bool has_vtp_cl     = false;   // vtp_cl_n + vtp_cl_time + vtp_cl_energy + vtp_cl_center + vtp_cl_blocks
    bool has_rf         = false;   // rf_n_a/b + rf_ns_a/b + cl_dt_rf
    bool has_gem_hits   = false;   // n_gem_hits + det_id + gem_* (replay_recon -gem_hit)
    bool has_gem_qa     = false;   // per-hit gem_x/y_time, gem_xy_dt/asym, gem_x/y_max_sdt, gem_x/y_min_corr
    bool has_gem_cl     = false;   // n_gem_cl + gem_cl_* per-cluster block
};

// Holder for ROOT vector-branch pointer-to-pointer binding.
// Keep one instance alive for as long as GetEntry is called.
struct ReconMatchVectorBindings {
    std::vector<uint16_t> *match_cl_idx = nullptr;
    std::vector<uint8_t>  *match_det_id = nullptr;
    std::vector<float>    *match_gem_x  = nullptr;
    std::vector<float>    *match_gem_y  = nullptr;
    std::vector<float>    *match_gem_z  = nullptr;
};

inline void BindReconMatchVectorBranches(TTree *tree,
                                         ReconEventData &ev,
                                         ReconMatchVectorBindings &b)
{
    const BranchBinder bind{tree};
    bind("match_cl_idx", ev.match_cl_idx, b.match_cl_idx);
    bind("match_det_id", ev.match_det_id, b.match_det_id);
    bind("match_gem_x",  ev.match_gem_x,  b.match_gem_x);
    bind("match_gem_y",  ev.match_gem_y,  b.match_gem_y);
    bind("match_gem_z",  ev.match_gem_z,  b.match_gem_z);
}

// Same holder pattern for the raw bank-word vectors (ssp_raw, vtp_*, tdc_*).
struct RawVectorBindings {
    std::vector<uint32_t> *ssp_raw      = nullptr;
    std::vector<uint32_t> *vtp_roc_tags = nullptr;
    std::vector<uint32_t> *vtp_nwords   = nullptr;
    std::vector<uint32_t> *vtp_words    = nullptr;
    std::vector<uint32_t> *tdc_roc_tags = nullptr;
    std::vector<uint32_t> *tdc_nwords   = nullptr;
    std::vector<uint32_t> *tdc_words    = nullptr;
};

// Event is RawEventData or ReconEventData.  Recon trees carry no tdc_*
// (RF is decoded into rf_*), so those holders stay null there.
template <class Event>
inline void BindRawVectorBranches(TTree *tree, Event &ev, RawVectorBindings &b)
{
    const BranchBinder bind{tree};
    bind("ssp_raw",      ev.ssp_raw,      b.ssp_raw);
    bind("vtp_roc_tags", ev.vtp_roc_tags, b.vtp_roc_tags);
    bind("vtp_nwords",   ev.vtp_nwords,   b.vtp_nwords);
    bind("vtp_words",    ev.vtp_words,    b.vtp_words);
    if constexpr (std::is_same_v<Event, RawEventData>) {
        bind("tdc_roc_tags", ev.tdc_roc_tags, b.tdc_roc_tags);
        bind("tdc_nwords",   ev.tdc_nwords,   b.tdc_nwords);
        bind("tdc_words",    ev.tdc_words,    b.tdc_words);
    }
}

// ── Raw "events" tree — write ────────────────────────────────────────────
inline void SetRawWriteBranches(TTree *tree, RawEventData &ev, bool with_peaks, bool Ecalib = false, bool noWaveform = false)
{
    if(!Ecalib) {
        tree->Branch("event_num",    &ev.event_num,    "event_num/I");
        tree->Branch("timestamp",    &ev.timestamp,    "timestamp/L");
    }
    tree->Branch("trigger_type", &ev.trigger_type, "trigger_type/b");
    tree->Branch("trigger_bits", &ev.trigger_bits, "trigger_bits/i");
    // Unified FADC250 channel array (HyCal + Veto + LMS); see RawEventData.
    tree->Branch("hycal.nch",         &ev.nch,         "hycal.nch/I");
    tree->Branch("hycal.module_id",   ev.module_id,    "hycal.module_id[hycal.nch]/s");
    tree->Branch("hycal.module_type", ev.module_type,  "hycal.module_type[hycal.nch]/b");
    if (!Ecalib && !noWaveform){
        tree->Branch("hycal.nsamples",    ev.nsamples,     "hycal.nsamples[hycal.nch]/b");
        tree->Branch("hycal.samples",     ev.samples,
                    Form("hycal.samples[hycal.nch][%d]/s", fdec::MAX_SAMPLES));
    }
    tree->Branch("hycal.gain_factor", ev.gain_factor,  "hycal.gain_factor[hycal.nch]/F");

    if (with_peaks) {
        // ped_* / peak_* are soft-analyzer products — only meaningful when
        // the analyzer ran.
        tree->Branch("hycal.ped_mean",      ev.ped_mean,       "hycal.ped_mean[hycal.nch]/F");
        tree->Branch("hycal.ped_rms",       ev.ped_rms,        "hycal.ped_rms[hycal.nch]/F");
        tree->Branch("hycal.ped_nused",     ev.ped_nused,      "hycal.ped_nused[hycal.nch]/b");
        tree->Branch("hycal.ped_quality",   ev.ped_quality,    "hycal.ped_quality[hycal.nch]/b");
        tree->Branch("hycal.ped_slope",     ev.ped_slope,      "hycal.ped_slope[hycal.nch]/F");
        tree->Branch("hycal.npeaks",        ev.npeaks,         "hycal.npeaks[hycal.nch]/b");
        tree->Branch("hycal.peak_height",   ev.peak_height,
                     Form("hycal.peak_height[hycal.nch][%d]/F",   fdec::MAX_PEAKS));
        tree->Branch("hycal.peak_time",     ev.peak_time,
                     Form("hycal.peak_time[hycal.nch][%d]/F",     fdec::MAX_PEAKS));
        tree->Branch("hycal.peak_integral", ev.peak_integral,
                     Form("hycal.peak_integral[hycal.nch][%d]/F", fdec::MAX_PEAKS));
        tree->Branch("hycal.peak_quality",  ev.peak_quality,
                     Form("hycal.peak_quality[hycal.nch][%d]/b",  fdec::MAX_PEAKS));

        // Firmware-mode (FADC250 Modes 1/2/3) emulation peaks.
        tree->Branch("hycal.daq_npeaks",        ev.daq_npeaks,    "hycal.daq_npeaks[hycal.nch]/b");
        tree->Branch("hycal.daq_peak_vp",       ev.daq_peak_vp,
                     Form("hycal.daq_peak_vp[hycal.nch][%d]/F",       fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_integral", ev.daq_peak_integral,
                     Form("hycal.daq_peak_integral[hycal.nch][%d]/F", fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_time",     ev.daq_peak_time,
                     Form("hycal.daq_peak_time[hycal.nch][%d]/F",     fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_cross",    ev.daq_peak_cross,
                     Form("hycal.daq_peak_cross[hycal.nch][%d]/I",    fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_pos",      ev.daq_peak_pos,
                     Form("hycal.daq_peak_pos[hycal.nch][%d]/I",      fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_coarse",   ev.daq_peak_coarse,
                     Form("hycal.daq_peak_coarse[hycal.nch][%d]/I",   fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_fine",     ev.daq_peak_fine,
                     Form("hycal.daq_peak_fine[hycal.nch][%d]/I",     fdec::MAX_PEAKS));
        tree->Branch("hycal.daq_peak_quality",  ev.daq_peak_quality,
                     Form("hycal.daq_peak_quality[hycal.nch][%d]/b",  fdec::MAX_PEAKS));
    }

    // GEM strip hit data (processed by GemSystem).
    if (!Ecalib){
        tree->Branch("gem.nch",    &ev.gem_nch,    "gem.nch/I");
        tree->Branch("gem.det",     ev.gem_det,    "gem.det[gem.nch]/b");
        tree->Branch("gem.plane",   ev.gem_plane,  "gem.plane[gem.nch]/b");
        tree->Branch("gem.strip",   ev.gem_strip,  "gem.strip[gem.nch]/I");
        tree->Branch("gem.charge",  ev.gem_charge, "gem.charge[gem.nch]/F");
        tree->Branch("gem.max_tb",  ev.gem_max_tb, "gem.max_tb[gem.nch]/S");
        tree->Branch("gem.pos",     ev.gem_pos,    "gem.pos[gem.nch]/F");
        tree->Branch("gem.xtalk",   ev.gem_xtalk,  "gem.xtalk[gem.nch]/b");
        tree->Branch("gem.ts_adc", ev.gem_ts_adc,
                        Form("gem.ts_adc[gem.nch][%d]/F", ssp::SSP_TIME_SAMPLES));

        // Raw 0xE10C SSP trigger bank words.
        tree->Branch("ssp_raw", &ev.ssp_raw);

        // Raw 0xE122 VTP bank words.  Flat triple of parallel vectors so ROOT
        // can serialize without a custom dictionary — see RawEventData for
        // the offset/decoding convention.
        tree->Branch("vtp_roc_tags", &ev.vtp_roc_tags);
        tree->Branch("vtp_nwords",   &ev.vtp_nwords);
        tree->Branch("vtp_words",    &ev.vtp_words);

        // Raw 0xE107 TDC bank words (RF reference + tagger).  Same layout
        // as the vtp_* triple — see RawEventData for the per-hit bit fields.
        tree->Branch("tdc_roc_tags", &ev.tdc_roc_tags);
        tree->Branch("tdc_nwords",   &ev.tdc_nwords);
        tree->Branch("tdc_words",    &ev.tdc_words);
    }
}

// ── Raw "events" tree — read ─────────────────────────────────────────────
// Binds the addresses of every branch that exists on `tree`; reports back
// which optional groups are present.
inline RawReadStatus SetRawReadBranches(TTree *tree, RawEventData &ev)
{
    RawReadStatus s;
    const BranchBinder bind{tree};

    bind("event_num",    &ev.event_num);
    bind("trigger_type", &ev.trigger_type);
    bind("trigger_bits", &ev.trigger_bits);
    bind("timestamp",    &ev.timestamp);

    bind("hycal.nch",         &ev.nch);
    bind("hycal.module_id",   ev.module_id);
    bind("hycal.module_type", ev.module_type);
    bind("hycal.nsamples",    ev.nsamples);
    bind("hycal.samples",     ev.samples);
    bind("hycal.gain_factor", ev.gain_factor);

    s.has_peaks = (tree->GetBranch("hycal.npeaks") != nullptr);
    if (s.has_peaks) {
        bind("hycal.ped_mean",      ev.ped_mean);
        bind("hycal.ped_rms",       ev.ped_rms);
        // ped_nused / ped_quality / ped_slope / peak_quality were added in
        // Mar-2026; files replayed earlier lack them.
        bind("hycal.ped_nused",     ev.ped_nused);
        bind("hycal.ped_quality",   ev.ped_quality);
        bind("hycal.ped_slope",     ev.ped_slope);
        bind("hycal.npeaks",        ev.npeaks);
        bind("hycal.peak_height",   ev.peak_height);
        bind("hycal.peak_time",     ev.peak_time);
        bind("hycal.peak_integral", ev.peak_integral);
        bind("hycal.peak_quality",  ev.peak_quality);
    }

    s.has_daq_peaks = (tree->GetBranch("hycal.daq_npeaks") != nullptr);
    if (s.has_daq_peaks) {
        bind("hycal.daq_npeaks",        ev.daq_npeaks);
        bind("hycal.daq_peak_vp",       ev.daq_peak_vp);
        bind("hycal.daq_peak_integral", ev.daq_peak_integral);
        bind("hycal.daq_peak_time",     ev.daq_peak_time);
        bind("hycal.daq_peak_cross",    ev.daq_peak_cross);
        bind("hycal.daq_peak_pos",      ev.daq_peak_pos);
        bind("hycal.daq_peak_coarse",   ev.daq_peak_coarse);
        bind("hycal.daq_peak_fine",     ev.daq_peak_fine);
        bind("hycal.daq_peak_quality",  ev.daq_peak_quality);
    }

    s.has_gem = (tree->GetBranch("gem.nch") != nullptr);
    if (s.has_gem) {
        bind("gem.nch",    &ev.gem_nch);
        bind("gem.det",     ev.gem_det);
        bind("gem.plane",   ev.gem_plane);
        bind("gem.strip",   ev.gem_strip);
        bind("gem.charge",  ev.gem_charge);
        bind("gem.max_tb",  ev.gem_max_tb);
        bind("gem.pos",     ev.gem_pos);
        bind("gem.xtalk",   ev.gem_xtalk);
        bind("gem.ts_adc",  ev.gem_ts_adc);
    }

    // ssp_raw / vtp_* / tdc_* are std::vector branches, bound separately by
    // BindRawVectorBranches.
    s.has_ssp_raw = (tree->GetBranch("ssp_raw") != nullptr);
    s.has_vtp_raw = (tree->GetBranch("vtp_words")    != nullptr)
                    && (tree->GetBranch("vtp_nwords")   != nullptr)
                    && (tree->GetBranch("vtp_roc_tags") != nullptr);
    s.has_tdc_raw = (tree->GetBranch("tdc_words")    != nullptr)
                    && (tree->GetBranch("tdc_nwords")   != nullptr)
                    && (tree->GetBranch("tdc_roc_tags") != nullptr);

    return s;
}

// ── Recon tree — write ───────────────────────────────────────────────────
inline void SetReconWriteBranches(TTree *tree, ReconEventData &ev, bool x17_mode, bool with_gem_hits = false)
{
    tree->Branch("event_num",    &ev.event_num,    "event_num/I");
    tree->Branch("trigger_type", &ev.trigger_type, "trigger_type/b");
    tree->Branch("trigger_bits", &ev.trigger_bits, "trigger_bits/i");
    tree->Branch("timestamp",    &ev.timestamp,    "timestamp/L");
    tree->Branch("total_energy", &ev.total_energy, "total_energy/F");

    // HyCal cluster branches (lab frame: target/beam-centred).
    tree->Branch("n_clusters", &ev.n_clusters, "n_clusters/I");
    tree->Branch("cl_x",       ev.cl_x,        "cl_x[n_clusters]/F");
    tree->Branch("cl_y",       ev.cl_y,        "cl_y[n_clusters]/F");
    tree->Branch("cl_z",       ev.cl_z,        "cl_z[n_clusters]/F");
    tree->Branch("cl_energy",  ev.cl_energy,   "cl_energy[n_clusters]/F");
    tree->Branch("cl_linear_corr", ev.cl_linear_corr, "cl_linear_corr[n_clusters]/F");
    tree->Branch("cl_bias_corr", ev.cl_bias_corr, "cl_bias_corr[n_clusters]/F");
    tree->Branch("cl_nblocks", ev.cl_nblocks,  "cl_nblocks[n_clusters]/b");
    tree->Branch("cl_npos",    ev.cl_npos,     "cl_npos[n_clusters]/b");
    tree->Branch("cl_center",  ev.cl_center,   "cl_center[n_clusters]/s");
    tree->Branch("cl_time",    ev.cl_time,     "cl_time[n_clusters]/F");
    tree->Branch("cl_flag",    ev.cl_flag,     "cl_flag[n_clusters]/i");

    // Flattened HyCal↔GEM matches across all clusters/chambers.
    tree->Branch("matchFlag", ev.matchFlag, "matchFlag[n_clusters]/i");
    tree->Branch("match_cl_idx", &ev.match_cl_idx);
    tree->Branch("match_det_id", &ev.match_det_id);
    tree->Branch("match_gem_x",  &ev.match_gem_x);
    tree->Branch("match_gem_y",  &ev.match_gem_y);
    tree->Branch("match_gem_z",  &ev.match_gem_z);

    // Quick-access matched pairs (clusters with ≥2 GEMs matched).
    tree->Branch("match_num", &ev.matchNum, "match_num/I");
    tree->Branch("mHit_E",  ev.mHit_E,  "mHit_E[match_num]/F");
    tree->Branch("mHit_x",  ev.mHit_x,  "mHit_x[match_num]/F");
    tree->Branch("mHit_y",  ev.mHit_y,  "mHit_y[match_num]/F");
    tree->Branch("mHit_z",  ev.mHit_z,  "mHit_z[match_num]/F");
    tree->Branch("mHit_gx", ev.mHit_gx, "mHit_gx[match_num][2]/F");
    tree->Branch("mHit_gy", ev.mHit_gy, "mHit_gy[match_num][2]/F");
    tree->Branch("mHit_gz", ev.mHit_gz, "mHit_gz[match_num][2]/F");
    tree->Branch("mHit_gid", ev.mHit_gid, "mHit_gid[match_num][2]/b");
    tree->Branch("mHit_cl_index", ev.mHit_cl_index, "mHit_cl_index[match_num]/b");

    // GEM hits (lab frame, per-detector plane).
    if (with_gem_hits) {
        tree->Branch("n_gem_hits",   &ev.n_gem_hits,   "n_gem_hits/I");
        tree->Branch("det_id",       ev.det_id,        "det_id[n_gem_hits]/b");
        tree->Branch("gem_x",        ev.gem_x,         "gem_x[n_gem_hits]/F");
        tree->Branch("gem_y",        ev.gem_y,         "gem_y[n_gem_hits]/F");
        tree->Branch("gem_z",        ev.gem_z,         "gem_z[n_gem_hits]/F");
        tree->Branch("gem_x_charge", ev.gem_x_charge,  "gem_x_charge[n_gem_hits]/F");
        tree->Branch("gem_y_charge", ev.gem_y_charge,  "gem_y_charge[n_gem_hits]/F");
        tree->Branch("gem_x_peak",   ev.gem_x_peak,    "gem_x_peak[n_gem_hits]/F");
        tree->Branch("gem_y_peak",   ev.gem_y_peak,    "gem_y_peak[n_gem_hits]/F");
        tree->Branch("gem_x_size",   ev.gem_x_size,    "gem_x_size[n_gem_hits]/b");
        tree->Branch("gem_y_size",   ev.gem_y_size,    "gem_y_size[n_gem_hits]/b");
        tree->Branch("gem_x_mTbin",  ev.gem_x_mTbin,   "gem_x_mTbin[n_gem_hits]/b");
        tree->Branch("gem_y_mTbin",  ev.gem_y_mTbin,   "gem_y_mTbin[n_gem_hits]/b");

        // Per-hit SBS-style quality variables (see gem::GEMHit; NaN = undefined).
        tree->Branch("gem_x_time",     ev.gem_x_time,     "gem_x_time[n_gem_hits]/F");
        tree->Branch("gem_y_time",     ev.gem_y_time,     "gem_y_time[n_gem_hits]/F");
        tree->Branch("gem_xy_dt",      ev.gem_xy_dt,      "gem_xy_dt[n_gem_hits]/F");
        tree->Branch("gem_xy_asym",    ev.gem_xy_asym,    "gem_xy_asym[n_gem_hits]/F");
        tree->Branch("gem_x_max_sdt",  ev.gem_x_max_sdt,  "gem_x_max_sdt[n_gem_hits]/F");
        tree->Branch("gem_y_max_sdt",  ev.gem_y_max_sdt,  "gem_y_max_sdt[n_gem_hits]/F");
        tree->Branch("gem_x_min_corr", ev.gem_x_min_corr, "gem_x_min_corr[n_gem_hits]/F");
        tree->Branch("gem_y_min_corr", ev.gem_y_min_corr, "gem_y_min_corr[n_gem_hits]/F");

        // GEM 1D clusters after filtering (detector-local position).
        tree->Branch("n_gem_cl",         &ev.n_gem_cl,        "n_gem_cl/I");
        tree->Branch("gem_cl_det",       ev.gem_cl_det,       "gem_cl_det[n_gem_cl]/b");
        tree->Branch("gem_cl_plane",     ev.gem_cl_plane,     "gem_cl_plane[n_gem_cl]/b");
        tree->Branch("gem_cl_size",      ev.gem_cl_size,      "gem_cl_size[n_gem_cl]/b");
        tree->Branch("gem_cl_mTbin",     ev.gem_cl_mTbin,     "gem_cl_mTbin[n_gem_cl]/b");
        tree->Branch("gem_cl_pos",       ev.gem_cl_pos,       "gem_cl_pos[n_gem_cl]/F");
        tree->Branch("gem_cl_peak",      ev.gem_cl_peak,      "gem_cl_peak[n_gem_cl]/F");
        tree->Branch("gem_cl_charge",    ev.gem_cl_charge,    "gem_cl_charge[n_gem_cl]/F");
        tree->Branch("gem_cl_time",      ev.gem_cl_time,      "gem_cl_time[n_gem_cl]/F");
        tree->Branch("gem_cl_seed_peak", ev.gem_cl_seed_peak, "gem_cl_seed_peak[n_gem_cl]/F");
        tree->Branch("gem_cl_seed_sum",  ev.gem_cl_seed_sum,  "gem_cl_seed_sum[n_gem_cl]/F");
        tree->Branch("gem_cl_max_sdt",   ev.gem_cl_max_sdt,   "gem_cl_max_sdt[n_gem_cl]/F");
        tree->Branch("gem_cl_min_corr",  ev.gem_cl_min_corr,  "gem_cl_min_corr[n_gem_cl]/F");
    }
    // Veto + LMS soft-peak summaries.
    if(!x17_mode) {
        tree->Branch("veto_nch",         &ev.veto_nch,         "veto_nch/I");
        tree->Branch("veto_id",          ev.veto_id,           "veto_id[veto_nch]/b");
        tree->Branch("veto_npeaks",      ev.veto_npeaks,       "veto_npeaks[veto_nch]/I");
        tree->Branch("veto_peak_time",   ev.veto_peak_time,
                    Form("veto_peak_time[veto_nch][%d]/F",     fdec::MAX_PEAKS));
        tree->Branch("veto_peak_height", ev.veto_peak_height,
                    Form("veto_peak_height[veto_nch][%d]/F",   fdec::MAX_PEAKS));
        tree->Branch("veto_peak_integral", ev.veto_peak_integral,
                    Form("veto_peak_integral[veto_nch][%d]/F", fdec::MAX_PEAKS));
    }

    tree->Branch("lms_nch",         &ev.lms_nch,         "lms_nch/I");
    tree->Branch("lms_id",          ev.lms_id,           "lms_id[lms_nch]/b");
    tree->Branch("lms_npeaks",      ev.lms_npeaks,       "lms_npeaks[lms_nch]/I");
    tree->Branch("lms_peak_time",   ev.lms_peak_time,
                 Form("lms_peak_time[lms_nch][%d]/F",     fdec::MAX_PEAKS));
    tree->Branch("lms_peak_height", ev.lms_peak_height,
                 Form("lms_peak_height[lms_nch][%d]/F",   fdec::MAX_PEAKS));
    tree->Branch("lms_peak_integral", ev.lms_peak_integral,
                 Form("lms_peak_integral[lms_nch][%d]/F", fdec::MAX_PEAKS));

    // Raw 0xE10C SSP trigger bank words.
    tree->Branch("ssp_raw", &ev.ssp_raw);

    // Raw 0xE122 VTP bank words — same flat triple as the raw events tree.
    tree->Branch("vtp_roc_tags", &ev.vtp_roc_tags);
    tree->Branch("vtp_nwords",   &ev.vtp_nwords);
    tree->Branch("vtp_words",    &ev.vtp_words);

    // VTP PRAD_CLUSTER online trigger data
    tree->Branch("vtp_cl_n",      &ev.vtp_cl_n,      "vtp_cl_n/I");
    tree->Branch("vtp_cl_time",   ev.vtp_cl_time,    "vtp_cl_time[vtp_cl_n]/s");
    tree->Branch("vtp_cl_energy", ev.vtp_cl_energy,  "vtp_cl_energy[vtp_cl_n]/s");
    tree->Branch("vtp_cl_center", ev.vtp_cl_center,  "vtp_cl_center[vtp_cl_n]/s");
    tree->Branch("vtp_cl_blocks", ev.vtp_cl_blocks,  "vtp_cl_blocks[vtp_cl_n]/b");

    // RF reference + per-cluster folded Δt.  See ReconEventData and
    // prad2det/include/RfTime.h for the folding rule.
    tree->Branch("rf_n_a",   &ev.rf_n_a,   "rf_n_a/b");
    tree->Branch("rf_n_b",   &ev.rf_n_b,   "rf_n_b/b");
    tree->Branch("rf_ns_a",  ev.rf_ns_a,
                 Form("rf_ns_a[%d]/F", tdc::RfTimeData::MAX_HITS_PER_CH));
    tree->Branch("rf_ns_b",  ev.rf_ns_b,
                 Form("rf_ns_b[%d]/F", tdc::RfTimeData::MAX_HITS_PER_CH));
    tree->Branch("cl_dt_rf", ev.cl_dt_rf, "cl_dt_rf[n_clusters]/F");
}

// ── Recon tree — read ────────────────────────────────────────────────────
inline ReconReadStatus SetReconReadBranches(TTree *tree, ReconEventData &ev)
{
    ReconReadStatus s;
    const BranchBinder bind{tree};

    // Defaults first, so branches missing from older files read the writer's
    // defaults (corrections 1, cl_dt_rf / GEM QA NaN, counts 0) instead of
    // stale or zero values.
    ev.clear();
    ev.fill_gem_qa_nan();

    bind("event_num",    &ev.event_num);
    bind("trigger_type", &ev.trigger_type);
    bind("trigger_bits", &ev.trigger_bits);
    bind("timestamp",    &ev.timestamp);
    bind("total_energy", &ev.total_energy);

    bind("n_clusters", &ev.n_clusters);
    bind("cl_x",       ev.cl_x);
    bind("cl_y",       ev.cl_y);
    bind("cl_z",       ev.cl_z);
    bind("cl_energy",  ev.cl_energy);
    bind("cl_linear_corr", ev.cl_linear_corr);
    bind("cl_bias_corr", ev.cl_bias_corr);
    bind("cl_nblocks", ev.cl_nblocks);
    bind("cl_npos",    ev.cl_npos);
    bind("cl_center",  ev.cl_center);
    bind("cl_time",    ev.cl_time);
    bind("cl_flag",    ev.cl_flag);

    s.has_per_cl_match = (tree->GetBranch("matchFlag") != nullptr);
    if (s.has_per_cl_match) {
        bind("matchFlag", ev.matchFlag);
        // match_* vectors: BindReconMatchVectorBranches.
    }

    s.has_match_num = (tree->GetBranch("match_num") != nullptr);
    if (s.has_match_num) {
        bind("match_num", &ev.matchNum);
        bind("mHit_E",  ev.mHit_E);
        bind("mHit_x",  ev.mHit_x);
        bind("mHit_y",  ev.mHit_y);
        bind("mHit_z",  ev.mHit_z);
        bind("mHit_gx", ev.mHit_gx);
        bind("mHit_gy", ev.mHit_gy);
        bind("mHit_gz", ev.mHit_gz);
        bind("mHit_gid", ev.mHit_gid);
        bind("mHit_cl_index", ev.mHit_cl_index);
    }

    bind("n_gem_hits",   &ev.n_gem_hits);
    bind("det_id",       ev.det_id);
    bind("gem_x",        ev.gem_x);
    bind("gem_y",        ev.gem_y);
    bind("gem_z",        ev.gem_z);
    bind("gem_x_charge", ev.gem_x_charge);
    bind("gem_y_charge", ev.gem_y_charge);
    bind("gem_x_peak",   ev.gem_x_peak);
    bind("gem_y_peak",   ev.gem_y_peak);
    bind("gem_x_size",   ev.gem_x_size);
    bind("gem_y_size",   ev.gem_y_size);
    bind("gem_x_mTbin",  ev.gem_x_mTbin);
    bind("gem_y_mTbin",  ev.gem_y_mTbin);
    s.has_gem_hits = (tree->GetBranch("n_gem_hits") != nullptr);

    s.has_gem_qa = (tree->GetBranch("gem_xy_dt") != nullptr);
    s.has_gem_cl = (tree->GetBranch("n_gem_cl") != nullptr);
    bind("gem_x_time",       ev.gem_x_time);
    bind("gem_y_time",       ev.gem_y_time);
    bind("gem_xy_dt",        ev.gem_xy_dt);
    bind("gem_xy_asym",      ev.gem_xy_asym);
    bind("gem_x_max_sdt",    ev.gem_x_max_sdt);
    bind("gem_y_max_sdt",    ev.gem_y_max_sdt);
    bind("gem_x_min_corr",   ev.gem_x_min_corr);
    bind("gem_y_min_corr",   ev.gem_y_min_corr);
    bind("n_gem_cl",         &ev.n_gem_cl);
    bind("gem_cl_det",       ev.gem_cl_det);
    bind("gem_cl_plane",     ev.gem_cl_plane);
    bind("gem_cl_size",      ev.gem_cl_size);
    bind("gem_cl_mTbin",     ev.gem_cl_mTbin);
    bind("gem_cl_pos",       ev.gem_cl_pos);
    bind("gem_cl_peak",      ev.gem_cl_peak);
    bind("gem_cl_charge",    ev.gem_cl_charge);
    bind("gem_cl_time",      ev.gem_cl_time);
    bind("gem_cl_seed_peak", ev.gem_cl_seed_peak);
    bind("gem_cl_seed_sum",  ev.gem_cl_seed_sum);
    bind("gem_cl_max_sdt",   ev.gem_cl_max_sdt);
    bind("gem_cl_min_corr",  ev.gem_cl_min_corr);

    s.has_veto = (tree->GetBranch("veto_nch") != nullptr);
    if (s.has_veto) {
        bind("veto_nch",          &ev.veto_nch);
        bind("veto_id",           ev.veto_id);
        bind("veto_npeaks",       ev.veto_npeaks);
        bind("veto_peak_time",    ev.veto_peak_time);
        bind("veto_peak_integral", ev.veto_peak_integral);
        bind("veto_peak_height",  ev.veto_peak_height);
    }

    s.has_lms = (tree->GetBranch("lms_nch") != nullptr);
    if (s.has_lms) {
        bind("lms_nch",          &ev.lms_nch);
        bind("lms_id",           ev.lms_id);
        bind("lms_npeaks",       ev.lms_npeaks);
        bind("lms_peak_time",    ev.lms_peak_time);
        bind("lms_peak_integral", ev.lms_peak_integral);
        bind("lms_peak_height",  ev.lms_peak_height);
    }

    // ssp_raw / vtp_* vectors are bound by BindRawVectorBranches.
    s.has_ssp_raw = (tree->GetBranch("ssp_raw") != nullptr);
    s.has_vtp_raw = (tree->GetBranch("vtp_words")    != nullptr)
                    && (tree->GetBranch("vtp_nwords")   != nullptr)
                    && (tree->GetBranch("vtp_roc_tags") != nullptr);

    // VTP PRAD_CLUSTER online trigger data — recon files replayed after 2026-06.
    s.has_vtp_cl = (tree->GetBranch("vtp_cl_n") != nullptr);
    if (s.has_vtp_cl) {
        bind("vtp_cl_n",      &ev.vtp_cl_n);
        bind("vtp_cl_time",   ev.vtp_cl_time);
        bind("vtp_cl_energy", ev.vtp_cl_energy);
        bind("vtp_cl_center", ev.vtp_cl_center);
        bind("vtp_cl_blocks", ev.vtp_cl_blocks);
    }

    // RF branches — recon files replayed after 2026-05.
    s.has_rf = (tree->GetBranch("rf_n_a") != nullptr);
    if (s.has_rf) {
        bind("rf_n_a",   &ev.rf_n_a);
        bind("rf_n_b",   &ev.rf_n_b);
        bind("rf_ns_a",  ev.rf_ns_a);
        bind("rf_ns_b",  ev.rf_ns_b);
        bind("cl_dt_rf", ev.cl_dt_rf);
    }

    return s;
}

// ── Scaler tree ("scalers") — write / read ───────────────────────────────
// Branch names mirror the field names so analysis can reference them
// directly.  Per-channel TRG/TDC arrays use a fixed [kDscChannels] length —
// the underlying DSC2 always emits 16 channels regardless of how many are
// hooked up — so no count branch is needed.
inline void SetScalerWriteBranches(TTree *tree, RawScalerData &sc)
{
    tree->Branch("event_number", &sc.event_number, "event_number/I");
    tree->Branch("ti_ticks",     &sc.ti_ticks,     "ti_ticks/L");
    tree->Branch("unix_time",    &sc.unix_time,    "unix_time/i");
    tree->Branch("sync_counter", &sc.sync_counter, "sync_counter/i");
    tree->Branch("run_number",   &sc.run_number,   "run_number/i");
    tree->Branch("trigger_type", &sc.trigger_type, "trigger_type/b");

    tree->Branch("slot",         &sc.slot,         "slot/I");
    tree->Branch("gated",        &sc.gated,        "gated/i");
    tree->Branch("ungated",      &sc.ungated,      "ungated/i");
    tree->Branch("live_ratio",   &sc.live_ratio,   "live_ratio/F");
    tree->Branch("source",       &sc.source,       "source/b");   // 0=ref 1=trg 2=tdc
    tree->Branch("channel",      &sc.channel,      "channel/b");

    tree->Branch("ref_gated",    &sc.ref_gated,    "ref_gated/i");
    tree->Branch("ref_ungated",  &sc.ref_ungated,  "ref_ungated/i");
    tree->Branch("trg_gated",    sc.trg_gated,
                 Form("trg_gated[%d]/i",   kDscChannels));
    tree->Branch("trg_ungated",  sc.trg_ungated,
                 Form("trg_ungated[%d]/i", kDscChannels));
    tree->Branch("tdc_gated",    sc.tdc_gated,
                 Form("tdc_gated[%d]/i",   kDscChannels));
    tree->Branch("tdc_ungated",  sc.tdc_ungated,
                 Form("tdc_ungated[%d]/i", kDscChannels));
}

inline void SetScalerReadBranches(TTree *tree, RawScalerData &sc)
{
    const BranchBinder bind{tree};
    bind("event_number", &sc.event_number);
    bind("ti_ticks",     &sc.ti_ticks);
    bind("unix_time",    &sc.unix_time);
    bind("sync_counter", &sc.sync_counter);
    bind("run_number",   &sc.run_number);
    bind("trigger_type", &sc.trigger_type);
    bind("slot",         &sc.slot);
    bind("gated",        &sc.gated);
    bind("ungated",      &sc.ungated);
    bind("live_ratio",   &sc.live_ratio);
    bind("source",       &sc.source);
    bind("channel",      &sc.channel);
    bind("ref_gated",    &sc.ref_gated);
    bind("ref_ungated",  &sc.ref_ungated);
    bind("trg_gated",    sc.trg_gated);
    bind("trg_ungated",  sc.trg_ungated);
    bind("tdc_gated",    sc.tdc_gated);
    bind("tdc_ungated",  sc.tdc_ungated);
}

// ── EPICS tree ("epics") — write / read ──────────────────────────────────
// `channel` and `value` are std::vector branches: readers that need them
// call BindEpicsVectorBranches after SetEpicsReadBranches.
inline void SetEpicsWriteBranches(TTree *tree, RawEpicsData &ep)
{
    tree->Branch("event_number_at_arrival", &ep.event_number_at_arrival,
                 "event_number_at_arrival/I");
    tree->Branch("ti_ticks_at_arrival", &ep.ti_ticks_at_arrival,
                 "ti_ticks_at_arrival/L");
    tree->Branch("unix_time",    &ep.unix_time,    "unix_time/i");
    tree->Branch("sync_counter", &ep.sync_counter, "sync_counter/i");
    tree->Branch("run_number",   &ep.run_number,   "run_number/i");
    tree->Branch("channel",      &ep.channel);
    tree->Branch("value",        &ep.value);
}

// ── Fillers — copy decoded prad2dec records into the side-tree structs ──
// All conversion / convention details (gated→live etc.) stay in prad2dec.
inline void FillScalerRow(const dsc::DscEventData &dsc,
                          const psync::SyncInfo &sync,
                          const fdec::EventInfo &info,
                          const evc::DaqConfig::DscScaler &cfg,
                          RawScalerData &out)
{
    using DSrc = evc::DaqConfig::DscScaler::Source;
    out = RawScalerData{};
    out.event_number = info.event_number;
    out.ti_ticks     = static_cast<long long>(info.timestamp);
    out.unix_time    = sync.unix_time;
    out.sync_counter = sync.sync_counter;
    out.run_number   = sync.run_number;
    out.trigger_type = info.trigger_type;

    out.slot       = dsc.slot;
    out.gated      = dsc.gated;
    out.ungated    = dsc.ungated;
    out.live_ratio = static_cast<float>(dsc.live_ratio());
    out.source     = (cfg.source == DSrc::Ref) ? 0
                   : (cfg.source == DSrc::Trg) ? 1 : 2;
    out.channel    = static_cast<uint8_t>(cfg.channel);

    out.ref_gated   = dsc.ref_gated;
    out.ref_ungated = dsc.ref_ungated;
    static_assert(kDscChannels == dsc::DSC2_NCH,
                  "RawScalerData/DscEventData channel count mismatch");
    std::memcpy(out.trg_gated,   dsc.trg_gated,   kDscChannels * sizeof(uint32_t));
    std::memcpy(out.trg_ungated, dsc.trg_ungated, kDscChannels * sizeof(uint32_t));
    std::memcpy(out.tdc_gated,   dsc.tdc_gated,   kDscChannels * sizeof(uint32_t));
    std::memcpy(out.tdc_ungated, dsc.tdc_ungated, kDscChannels * sizeof(uint32_t));
}

inline void FillEpicsRow(const epics::EpicsRecord &rec, RawEpicsData &out)
{
    out = RawEpicsData{};
    out.event_number_at_arrival = rec.event_number_at_arrival;
    out.ti_ticks_at_arrival     = static_cast<long long>(rec.timestamp_at_arrival);
    out.unix_time               = rec.unix_time;
    out.sync_counter            = rec.sync_counter;
    out.run_number              = rec.run_number;
    out.channel                 = rec.channel;
    out.value                   = rec.value;
}

inline void SetEpicsReadBranches(TTree *tree, RawEpicsData &ep)
{
    const BranchBinder bind{tree};
    bind("event_number_at_arrival", &ep.event_number_at_arrival);
    bind("ti_ticks_at_arrival",     &ep.ti_ticks_at_arrival);
    bind("unix_time",    &ep.unix_time);
    bind("sync_counter", &ep.sync_counter);
    bind("run_number",   &ep.run_number);
}

// Same holder pattern for the EPICS vector branches.
struct EpicsVectorBindings {
    std::vector<std::string> *channel = nullptr;
    std::vector<double>      *value   = nullptr;
};

inline void BindEpicsVectorBranches(TTree *tree, RawEpicsData &ep,
                                    EpicsVectorBindings &b)
{
    const BranchBinder bind{tree};
    bind("channel", ep.channel, b.channel);
    bind("value",   ep.value,   b.value);
}

// ── Run-info tree ("runinfo") — write / read ─────────────────────────────
// One row per CODA control event (PRESTART / GO / END); the PRESTART row
// is the one carrying the long `daq_config` text.  The string branch is
// stored directly (ROOT serializes std::string fine); readers that need it
// bind a held pointer: BranchBinder{tree}("daq_config", ri.daq_config, sp).
inline void SetRunInfoWriteBranches(TTree *tree, RawRunInfo &ri)
{
    tree->Branch("run_number", &ri.run_number, "run_number/i");
    tree->Branch("unix_time",  &ri.unix_time,  "unix_time/i");
    tree->Branch("run_type",   &ri.run_type,   "run_type/b");
    tree->Branch("event_tag",  &ri.event_tag,  "event_tag/b");
    tree->Branch("daq_config", &ri.daq_config);
}

inline void SetRunInfoReadBranches(TTree *tree, RawRunInfo &ri)
{
    const BranchBinder bind{tree};
    bind("run_number", &ri.run_number);
    bind("unix_time",  &ri.unix_time);
    bind("run_type",   &ri.run_type);
    bind("event_tag",  &ri.event_tag);
}

inline void FillRunInfoRow(const psync::SyncInfo &sync,
                           const std::string     &daq_config_text,
                           RawRunInfo            &out)
{
    out = RawRunInfo{};
    out.run_number = sync.run_number;
    out.unix_time  = sync.unix_time;
    out.run_type   = sync.run_type;
    out.event_tag  = static_cast<uint8_t>(sync.event_tag & 0xFFu);
    out.daq_config = daq_config_text;
}

// ── LMS gain-correction tree ("gain") — write / read ─────────────────────
// Only peak-level data is kept (no waveform samples): HyCal channels are
// stored with the same hycal.* prefix used in the events tree so that
// existing channel-lookup helpers work unchanged (LMS PMTs share the same
// flat array, module_type MOD_LMS).
inline void SetLMSWriteBranches(TTree *tree, LMSEventData &ev)
{
    tree->Branch("event_num",    &ev.event_num,    "event_num/I");
    tree->Branch("trigger_type", &ev.trigger_type, "trigger_type/b");
    tree->Branch("trigger_bits", &ev.trigger_bits, "trigger_bits/i");
    tree->Branch("timestamp",    &ev.timestamp,    "timestamp/L");
    tree->Branch("event_type",   &ev.event_type,   "event_type/I");

    tree->Branch("hycal.nch",         &ev.nch,        "hycal.nch/I");
    tree->Branch("hycal.module_id",   ev.module_id,   "hycal.module_id[hycal.nch]/s");
    tree->Branch("hycal.module_type", ev.module_type, "hycal.module_type[hycal.nch]/b");
    tree->Branch("hycal.npeaks",      ev.npeaks,      "hycal.npeaks[hycal.nch]/b");
    tree->Branch("hycal.peak_height",   ev.peak_height,
                 Form("hycal.peak_height[hycal.nch][%d]/F",   fdec::MAX_PEAKS));
    tree->Branch("hycal.peak_time",     ev.peak_time,
                 Form("hycal.peak_time[hycal.nch][%d]/F",     fdec::MAX_PEAKS));
    tree->Branch("hycal.peak_integral", ev.peak_integral,
                 Form("hycal.peak_integral[hycal.nch][%d]/F", fdec::MAX_PEAKS));
}

inline void SetLMSReadBranches(TTree *tree, LMSEventData &ev)
{
    const BranchBinder bind{tree};

    bind("event_num",    &ev.event_num);
    bind("trigger_type", &ev.trigger_type);
    bind("trigger_bits", &ev.trigger_bits);
    bind("timestamp",    &ev.timestamp);
    bind("event_type",   &ev.event_type);

    bind("hycal.nch",         &ev.nch);
    bind("hycal.module_id",   ev.module_id);
    bind("hycal.module_type", ev.module_type);
    bind("hycal.npeaks",      ev.npeaks);
    bind("hycal.peak_height",   ev.peak_height);
    bind("hycal.peak_time",     ev.peak_time);
    bind("hycal.peak_integral", ev.peak_integral);
}

} // namespace prad2
