#pragma once
//=============================================================================
// GainCorrCompute.h — shared LMS/alpha gain-correction batch calculation
//=============================================================================

#include <TH1F.h>
#include <TTree.h>

#include "EventData.h"
#include "HyCalSystem.h"
#include "gain_factor.h"

#include <cstdint>
#include <map>
#include <set>
#include <string>
#include <vector>

namespace analysis {

static constexpr int   kGainNW        = 1156;
static constexpr int   kGainNLMS      = 3;
static constexpr int   kGainLMSIDBase = prad2::kLmsIdBase + 1;  // LMS1..3 = 3101..3103
static constexpr int   kGainWIDBase   = fdec::PWO_ID0;          // PbWO4 module_id = W-id + 1000
static constexpr int   kGainHistBins  = 600;
static constexpr float kGainHistMin   = 0.f;
static constexpr float kGainHistMax   = 15000.f;

// Default LMS events per gain_corr batch: offline replay vs. online monitor.
static constexpr int   kGainBatchSizeOffline = 4000;
static constexpr int   kGainBatchSizeOnline  = 1000;

struct GainBatch {
    int   batch_id         = 0;
    int   event_num_start  = 0;
    int   event_num_end    = 0;
    int   n_lms_events     = 0;
    int   n_alpha_events   = 0;
    int   ref_run          = 0;
    uint32_t unix_time     = 0;  // absolute Unix seconds at batch midpoint

    float refPMT_ratio       [kGainNLMS]            = {};
    float gain_W             [kGainNW][kGainNLMS]   = {};
    float gain_W_ref         [kGainNW][kGainNLMS]   = {};
    float gain_corr_W        [kGainNW][kGainNLMS]   = {};
    float fit_mean_ref_lms   [kGainNLMS]            = {};
    float fit_mean_ref_alpha [kGainNLMS]            = {};
    float fit_mean_W_lms     [kGainNW]              = {};
};

struct GainPlotConfig {
    bool          enabled   = false;
    int           max_hists = 10;
    std::set<int> w_ids;
};

struct GainPlotStore {
    std::vector<TH1F*> ref_lms  [kGainNLMS];
    std::vector<TH1F*> ref_alpha[kGainNLMS];
    std::map<int, std::vector<TH1F*>> mod_w;  // key = W-id (1-based)

    ~GainPlotStore();
};

// LMS/alpha amplitude histogram: kGainHistBins bins over [kGainHistMin,
// kGainHistMax], detached from any directory (the caller owns it).
TH1F *MakeGainHist(const char *name);

std::string MakeLMSOutputFile(const std::string &evio_path);

// prad2::GainCorrDir(db_dir) + "/prad_RUN_gain_corr.root"; creates the
// directory (warning only on failure).
std::string GainCorrOutputPath(const std::string &db_dir, int run);

// Replay every EVIO file through Process_LMSgainFactor into
// out_dir/MakeLMSOutputFile(file) with up to num_threads workers, logging
// "[i/N]" progress.  Returns the number of failed files; the outputs that were
// written are appended to *produced in input order.
int ReplayLMSFiles(const std::vector<std::string> &evio_files,
                   const std::string              &out_dir,
                   int                             num_threads,
                   const std::string              &daq_config,
                   const std::string              &daq_map,
                   const std::string              &db_dir,
                   std::vector<std::string>       *produced = nullptr);

// For each run among evio_files that has no gain_corr file in
// prad2::GainCorrDir(db_dir) yet, run prad2ana_replay_gainCorr (installed next
// to this executable) on that run's files.  Failures are only reported; the
// replay then proceeds with identity gain correction for that run.
void EnsureGainCorr(const std::vector<std::string> &evio_files,
                    const std::string              &db_dir,
                    const std::string              &daq_config,
                    const std::string              &daq_map,
                    int                             num_threads);

void SetupGainBranches(TTree *tree, GainBatch &b);

void FlushGainBatch(GainBatch &b, TTree *tree,
                    TH1F *mod_lms[kGainNW],
                    TH1F *ref_lms[kGainNLMS],
                    TH1F *ref_alpha[kGainNLMS],
                    const prad2::RefGainTable &ref_tbl);

bool ComputeGainCorrections(const std::vector<std::string> &lms_files,
                            const std::string             &gain_out,
                            int                            batch_size,
                            int                            ref_run_num,
                            const prad2::RefGainTable     &ref_tbl,
                            const GainPlotConfig          *plot_cfg = nullptr,
                            GainPlotStore                 *plot_store = nullptr);

} // namespace analysis
