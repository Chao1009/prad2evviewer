#pragma once
//=============================================================================
// gain_factor.h — per-module LMS gain factors and time-dependent corrections
//
// ── Reference gain (.dat files, produced by refGain_produce) ─────────────────
//   Directory: RefGainDir(<db>) = <db>/gain_factor/ref_gain/
//   File:      prad_XXXXXX_LMS.dat  (7 columns: Name lms_peak lms_sigma
//              lms_chi2/ndf g1 g2 g3)
//   W/G module lines only; LMS header lines ignored.
//
//   auto tbl = prad2::LoadRefGain(prad2::RefGainDir(db), run_num);
//   float g1_W1 = tbl.w[1].g[0];
//
// ── Time-dependent correction (.root files, produced by replay_gainCorr) ──────
//   Directory: GainCorrDir(<db>) = <db>/gain_factor/gain_correction/
//   File:      prad_XXXXXX_gain_corr.root  (TTree "gain_corr", one entry/batch)
//   Only W (PbWO4) modules; G (PbGlass) corrections not stored in these files.
//   Writers always use GainCorrDir(<db>); readers use RunConfig::gain_data_dir,
//   which defaults to the same directory (general.json gain_factor.data_dir
//   only redirects the reader).  The correction methods are described at
//   LoadGainCorrTimeSeries().
//
//   // One-time setup (single-threaded):
//   auto ts = prad2::LoadGainCorrTimeSeries(gRunConfig, run_num);
//   // Per-event lookup (read-only → safe from multiple threads after init):
//   const auto& corr = ts.GetCorr(event_num);
//   new_adc2mev = old_adc2mev * corr.ModuleGain(module_id);
//
// Both kinds of file are looked up by exact run number (FindRunFile).
//=============================================================================

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <regex>
#include <string>
#include <vector>

#include <TROOT.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TTree.h>

#include "RunInfoConfig.h"

namespace prad2 {

inline std::string RefGainDir(const std::string &db_dir)
{
    return db_dir + "/gain_factor/ref_gain";
}

inline std::string GainCorrDir(const std::string &db_dir)
{
    return db_dir + "/gain_factor/gain_correction";
}

// Single module's three gain factors (g1, g2, g3 from the LMS/alpha fit).
// Zero-initialised by default so missing entries are safe to use.
struct RefGainFactor {
    float g[3] = {0.f, 0.f, 0.f};  // g[0]=g1, g[1]=g2, g[2]=g3
};

// Full table for one run.  Arrays indexed by module numeric ID; index 0 is
// unused (modules are 1-based in the dat file).
struct RefGainTable {
    static constexpr int MAX_W = 1157;  // W1 .. W1156
    static constexpr int MAX_G = 901;   // G1 .. G900

    RefGainFactor w[MAX_W];  // PWO crystal modules
    RefGainFactor g[MAX_G];  // PbGlass modules
    int  run_number = -1;
    bool loaded     = false;
};

// Path of the regular file dir/prad_<run_num as %06d><suffix>, or an empty
// string when there is none (run_num outside 0..999999 never matches).
inline std::string FindRunFile(const std::string &dir, int run_num, const char *suffix)
{
    if (dir.empty() || run_num < 0 || run_num > 999999) return {};
    char name[64];
    std::snprintf(name, sizeof(name), "prad_%06d%s", run_num, suffix);
    const std::filesystem::path path = std::filesystem::path(dir) / name;
    std::error_code ec;
    return std::filesystem::is_regular_file(path, ec) ? path.string() : std::string{};
}

inline std::string FindRefGainFile(const std::string &dir, int run_num)
{
    return FindRunFile(dir, run_num, "_LMS.dat");
}

// Load reference gain factors from one explicit .dat file.
// On any failure returns a default-constructed (unloaded) table.
inline RefGainTable LoadRefGainFile(const std::string &path)
{
    RefGainTable tbl;

    std::ifstream f(path);
    if (!f) {
        std::cerr << "Warning: cannot open gain factor file " << path << "\n";
        return tbl;
    }

    // Parse the embedded run number from the file name.
    {
        static const std::regex pat(R"(prad_(\d{6})_LMS\.dat)");
        std::smatch m;
        std::string fname = std::filesystem::path(path).filename().string();
        if (std::regex_search(fname, m, pat))
            tbl.run_number = std::stoi(m[1].str());
    }

    // Skip the header line (contains string column names, not numeric data).
    { std::string header; std::getline(f, header); }

    std::string name;
    float col2, col3, col4, g1, g2, g3;
    while (f >> name >> col2 >> col3 >> col4 >> g1 >> g2 >> g3) {
        // W/G module lines only; LMS lines and anything else are skipped.
        if (name.empty() || (name[0] != 'W' && name[0] != 'G')) continue;
        const bool is_w = name[0] == 'W';
        const int id = std::stoi(name.substr(1));
        if (id >= 1 && id < (is_w ? RefGainTable::MAX_W : RefGainTable::MAX_G))
            (is_w ? tbl.w : tbl.g)[id] = {{g1, g2, g3}};
    }

    tbl.loaded = true;
    std::cerr << "RefGain: loaded run " << tbl.run_number
              << " from " << path << "\n";
    return tbl;
}

// Load reference gain factors for the given run number.
// On any failure returns a default-constructed (unloaded) table.
inline RefGainTable LoadRefGain(const std::string &dir, int run_num)
{
    std::string path = FindRefGainFile(dir, run_num);
    if (path.empty()) {
        std::cerr << "Warning: no gain factor file found in " << dir
                  << " for run " << run_num << "\n";
        return RefGainTable{};
    }
    return LoadRefGainFile(path);
}

// Per-module gain correction factor: correction[id] = g_ref / g_current.
// Applying this to the current ADC->MeV scale compensates for gain drift.
// A value of 1.0 means no correction needed; > 1.0 means gain dropped.
struct GainCorrTable {
    static constexpr int MAX_W = RefGainTable::MAX_W;
    static constexpr int MAX_G = RefGainTable::MAX_G;

    // correction[id][j] = ref.g[j] / cur.g[j]  (j = 0,1,2 for g1,g2,g3)
    // avg[id]           = mean of the three per-LMS corrections, a failed
    //                     fit counting as 1
    struct Entry {
        float corr[3] = {1.f, 1.f, 1.f};
        float avg      = 1.f;
    };

    Entry w[MAX_W];
    Entry g[MAX_G];

    int ref_run = -1;
    int cur_run = -1;

    // Gain factor for a module id (PbGlass: G number, PbWO4: 1000 + W
    // number).  PbWO4 uses the mean of the LMS2 and LMS3 corrections (LMS1
    // is not used), PbGlass its avg; any other id (Veto, LMS) gets 1.
    float ModuleGain(int id) const
    {
        if (id > 1000 && id - 1000 < MAX_W)
            return (w[id - 1000].corr[1] + w[id - 1000].corr[2]) / 2.0f;
        if (id > 0 && id < MAX_G) return g[id].avg;
        return 1.f;
    }
};

// ── Time-dependent gain correction from replay_gainCorr ROOT output ───────────

inline std::string FindGainCorrRootFile(const std::string &dir, int run_num)
{
    return FindRunFile(dir, run_num, "_gain_corr.root");
}

// A time series of gain correction tables loaded from a replay_gainCorr ROOT
// file (one entry per LMS batch, including its midpoint Unix time when the
// source file supplied a usable EPICS/scaler anchor).  After construction the
// object is read-only and therefore safe to access concurrently from multiple
// threads.
struct GainCorrTimeSeries {
    struct Batch {
        int           event_num_start = 0;
        int           event_num_end   = 0;
        uint32_t      unix_time       = 0;
        GainCorrTable corr;
    };

    std::vector<Batch> batches;   // sorted ascending by event_num_start
    int  run_num = -1;
    bool loaded  = false;

    // Return the correction table for the given event number.
    // Finds the last batch whose event_num_start <= event_num.
    // Falls back to the first batch if event_num precedes all batches.
    const GainCorrTable &GetCorr(int event_num) const noexcept
    {
        // C++11 guarantees thread-safe initialisation of function-scope statics.
        static const GainCorrTable kIdentity{};
        if (batches.empty()) return kIdentity;
        for (auto it = batches.rbegin(); it != batches.rend(); ++it)
            if (event_num >= it->event_num_start) return it->corr;
        return batches.front().corr;
    }
};

// Load the gain correction time series from the replay_gainCorr ROOT output.
//
// run_cfg         : RunConfig for this run — supplies gain_data_dir (the
//                    directory holding the gain_corr files) and gain_ref_run
//                    (the reference run number to divide into; the ref_run
//                    stored in the current run's file is unreliable and not
//                    read).
// run_num         : run whose prad_XXXXXX_gain_corr.root is loaded
// use_precomputed : false (default) — compute the correction at load time as
//                    ref_run.gain_W / cur_run.gain_W, where ref_run.gain_W is
//                    averaged over every batch in the reference run's own
//                    gain_corr.root (one average per module/LMS channel);
//                    the gain_corr_W branch is ignored.  true — use the
//                    pre-computed gain_corr_W branch instead (legacy
//                    behaviour: ref_tbl.g from the .dat file / cur_run's
//                    gain_W, computed once by replay_gainCorr).
//
// Thread safety
//   Call once from a single thread during setup.
//   The returned GainCorrTimeSeries may then be shared across threads for
//   read-only access via GetCorr().  Requires ROOT::EnableThreadSafety() to
//   have been called before spawning worker threads.
inline GainCorrTimeSeries LoadGainCorrTimeSeries(const RunConfig &run_cfg,
                                                  int              run_num,
                                                  bool             use_precomputed = false)
{
    const std::string corr_dir = run_cfg.gain_data_dir;
    const int         ref_run  = run_cfg.gain_ref_run;

    GainCorrTimeSeries ts;
    ts.run_num = run_num;

    std::string path = FindGainCorrRootFile(corr_dir, run_num);
    if (path.empty()) {
        std::cerr << "Warning: no gain_corr root file in " << corr_dir
                  << " for run " << run_num << "\n";
        return ts;
    }

    TFile *f = TFile::Open(path.c_str(), "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Warning: cannot open " << path << "\n";
        delete f;
        return ts;
    }

    TTree *tree = nullptr;
    f->GetObject("gain_corr", tree);
    if (!tree) {
        std::cerr << "Warning: no 'gain_corr' tree in " << path << "\n";
        delete f;
        return ts;
    }

    // Written with gain_corr_W[N_W][N_LMS] or gain_W[N_W][N_LMS], N_W = 1156, N_LMS = 3.
    static constexpr int kNW   = GainCorrTable::MAX_W - 1;  // 1156
    static constexpr int kNLMS = 3;

    int      ev_start = 0, ev_end = 0;
    uint32_t unix_time = 0;
    float corr_W[kNW][kNLMS];
    float gain_W[kNW][kNLMS];

    tree->SetBranchAddress("event_num_start", &ev_start);
    tree->SetBranchAddress("event_num_end",   &ev_end);
    if (use_precomputed)
        tree->SetBranchAddress("gain_corr_W", corr_W);
    else
        tree->SetBranchAddress("gain_W",      gain_W);
    const bool has_unix_time = tree->GetBranch("unix_time") != nullptr;
    if (has_unix_time)
        tree->SetBranchAddress("unix_time", &unix_time);

    // Default method: average the reference run's own gain_W over all
    // of its batches, then divide it into every current-run batch below.
    float ref_gain_W[kNW][kNLMS];
    bool  have_ref_gain = false;
    if (!use_precomputed) {
        if (ref_run < 0) {
            std::cerr << "Warning: RunConfig has no valid gain_ref_run\n";
        } else {
            std::string ref_path = FindGainCorrRootFile(corr_dir, ref_run);
            if (ref_path.empty()) {
                std::cerr << "Warning: no gain_corr root file in " << corr_dir
                          << " for ref run " << ref_run << "\n";
            } else {
                TFile *ref_f = TFile::Open(ref_path.c_str(), "READ");
                if (!ref_f || ref_f->IsZombie()) {
                    std::cerr << "Warning: cannot open ref-run gain_corr file "
                              << ref_path << "\n";
                } else {
                    TTree *ref_tree = nullptr;
                    ref_f->GetObject("gain_corr", ref_tree);
                    const Long64_t n_ref_entries = ref_tree ? ref_tree->GetEntries() : 0;
                    if (n_ref_entries == 0) {
                        std::cerr << "Warning: no 'gain_corr' entries in "
                                  << ref_path << "\n";
                    } else {
                        float ref_batch_W[kNW][kNLMS];
                        float ref_sum[kNW][kNLMS] = {};
                        int   ref_cnt[kNW][kNLMS] = {};
                        ref_tree->SetBranchAddress("gain_W", ref_batch_W);
                        for (Long64_t rk = 0; rk < n_ref_entries; ++rk) {
                            ref_tree->GetEntry(rk);
                            for (int wi = 0; wi < kNW; ++wi) {
                                for (int j = 0; j < kNLMS; ++j) {
                                    if (ref_batch_W[wi][j] <= 0.f) continue;
                                    ref_sum[wi][j] += ref_batch_W[wi][j];
                                    ++ref_cnt[wi][j];
                                }
                            }
                        }
                        for (int wi = 0; wi < kNW; ++wi)
                            for (int j = 0; j < kNLMS; ++j)
                                ref_gain_W[wi][j] = (ref_cnt[wi][j] > 0)
                                    ? ref_sum[wi][j] / ref_cnt[wi][j] : 0.f;
                        have_ref_gain = true;
                    }
                }
                delete ref_f;
            }
        }
    }

    const Long64_t nentries = tree->GetEntries();
    ts.batches.reserve(static_cast<size_t>(nentries));

    for (Long64_t ie = 0; ie < nentries; ++ie) {
        tree->GetEntry(ie);

        GainCorrTimeSeries::Batch b;
        b.event_num_start = ev_start;
        b.event_num_end   = ev_end;
        b.unix_time       = has_unix_time ? unix_time : 0;
        b.corr.cur_run    = run_num;
        b.corr.ref_run    = ref_run;

        for (int wi = 0; wi < kNW; ++wi) {
            float sum = 0.f;
            for (int j = 0; j < kNLMS; ++j) {
                float v;
                if (use_precomputed) {
                    // 0 in the ROOT file signals a failed fit — treat as identity.
                    v = (corr_W[wi][j] > 0.f) ? corr_W[wi][j] : 1.f;
                } else {
                    v = (have_ref_gain && ref_gain_W[wi][j] > 0.f && gain_W[wi][j] > 0.f)
                        ? ref_gain_W[wi][j] / gain_W[wi][j] : 1.f;
                }
                b.corr.w[wi + 1].corr[j] = v;
                sum += v;
            }
            b.corr.w[wi + 1].avg = sum / static_cast<float>(kNLMS);
        }
        ts.batches.push_back(std::move(b));
    }

    delete f;   // TFile owns TTree; both are freed here

    // Ensure ascending order (defensive; writer should already sort).
    std::sort(ts.batches.begin(), ts.batches.end(),
              [](const GainCorrTimeSeries::Batch &a,
                 const GainCorrTimeSeries::Batch &b) {
                  return a.event_num_start < b.event_num_start;
              });

    ts.loaded = true;
    std::cerr << "GainCorrTS: run " << run_num << ": "
              << ts.batches.size() << " batches from " << path
              << (use_precomputed ? " (precomputed gain_corr_W)"
                                  : " (ref-run gain_W ratio)") << "\n";
    return ts;
}

// Gaussian fit of an LMS/alpha amplitude peak (all zero when the fit failed).
struct FitResult{
    float mean     = 0.;
    float sigma    = 0.;
    float chi2pndf = 0.;
};

inline FitResult gain_hist_fitter(TH1F* h, const float & fac)
{
    FitResult result;

    if (!h) {
        return result;
    }

    const int nBins = h->GetNbinsX();
    if (nBins <= 0) {
        return result;
    }

    const int maxBin = h->GetMaximumBin();
    const double maxContent = h->GetBinContent(maxBin);

    if (maxContent <= 0.) {
        return result;
    }

    const double threshold = fac * maxContent;

    int leftBin  = maxBin;
    int rightBin = maxBin;

    // Find first bin to the left below threshold (fac * max)
    for (int ibin = maxBin; ibin >= 1; --ibin) {
        if (h->GetBinContent(ibin) < threshold) {
            leftBin = ibin;
            break;
        }
        if (ibin == 1) {
            leftBin = 1;
        }
    }

    // Find first bin to the right below threshold
    for (int ibin = maxBin; ibin <= nBins; ++ibin) {
        if (h->GetBinContent(ibin) < threshold) {
            rightBin = ibin;
            break;
        }
        if (ibin == nBins) {
            rightBin = nBins;
        }
    }

    // Keep the fit range inside the above-threshold region: shift inward by
    // one bin when the threshold-crossing bin itself is below threshold.
    if (leftBin < maxBin && h->GetBinContent(leftBin) < threshold) {
        leftBin++;
    }
    if (rightBin > maxBin && h->GetBinContent(rightBin) < threshold) {
        rightBin--;
    }

    // Ensure at least 5 bins around the maximum are included so that a
    // 3-parameter Gaussian fit is always well-constrained, even when the
    // peak is very narrow (only a few filled bins).
    if (leftBin  > maxBin - 2) leftBin  = maxBin - 2;
    if (rightBin < maxBin + 2) rightBin = maxBin + 2;

    if (leftBin < 1) leftBin = 1;
    if (rightBin > nBins) rightBin = nBins;
    if (leftBin >= rightBin) {
        return result;
    }
    
    const double xLow  = h->GetXaxis()->GetBinLowEdge(leftBin);
    const double xHigh = h->GetXaxis()->GetBinUpEdge(rightBin);

    const double peakX = h->GetXaxis()->GetBinCenter(maxBin);

    // A reasonable initial sigma guess from fit window width
    double sigmaGuess = 0.5 * (xHigh - xLow) / 2.0;
    if (sigmaGuess <= 0.) {
        sigmaGuess = h->GetRMS();
    }
    if (sigmaGuess <= 0.) {
        sigmaGuess = h->GetBinWidth(maxBin);
    }

    std::string fitName = std::string(h->GetName()) + "_gaus_fit";
    TF1 * gausFit = new TF1(fitName.c_str(), "gaus", xLow, xHigh);
    gausFit->SetParameters(maxContent, peakX, sigmaGuess);

    // R = fit in function range, Q = quiet, N = do not store function in histogram
    // "N" is required to avoid double-ownership: new TF1 is registered in gROOT's
    // global list; without "N", Fit() also stores it in the histogram's list,
    // causing a double-free during ROOT cleanup at program exit.
    int fitStatus = h->Fit(gausFit, "RQN");

    if (fitStatus != 0) {
        delete gausFit;
        return result;
    }

    const float mean  = static_cast<float>(gausFit->GetParameter(1));
    const float sigma = static_cast<float>(gausFit->GetParameter(2));
    if (!std::isfinite(mean) || !std::isfinite(sigma) || sigma <= 0.f) {
        delete gausFit;
        return result;
    }
    result.mean  = mean;
    result.sigma = sigma;

    const double ndf = gausFit->GetNDF();
    if (ndf > 0) {
        result.chi2pndf = static_cast<float>(gausFit->GetChisquare() / ndf);
    }

    // Transfer ownership from gROOT's global list to the histogram so the fit
    // curve is saved in the output file and there is only one owner (no double-free).
    gROOT->GetListOfFunctions()->Remove(gausFit);
    h->GetListOfFunctions()->Add(gausFit);

    return result;
}

} // namespace prad2
