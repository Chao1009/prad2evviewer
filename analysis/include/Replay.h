#pragma once
//=============================================================================
// Replay.h — convert raw DAQ data (EVIO) to ROOT trees
//
// Decodes EVIO events and writes per-channel waveform/peak data to a TTree.
// Depends on prad2dec (decoder) and ROOT (TFile/TTree).
//=============================================================================

#include "EvChannel.h"
#include "EventData.h"
#include "WaveAnalyzer.h"
#include "Fadc250FwAnalyzer.h"
#include "DaqConfig.h"
#include "ConfigSetup.h"
#include "load_daq_config.h"
#include "HyCalSystem.h"
#include "HyCalCluster.h"

#include <TFile.h>
#include <TTree.h>

#include <array>
#include <memory>
#include <string>
#include <tuple>
#include <vector>

namespace gem { class GemSystem; class GemCluster; struct GEMHit; struct StripCluster; }

namespace analysis {

using EventVars       = prad2::RawEventData;
using EventVars_Recon = prad2::ReconEventData;
using LMSEventVars    = prad2::LMSEventData;

class Replay
{
public:
    Replay() = default;

    // Load DAQ configuration (event tags, ADC format, etc.).
    void LoadDaqConfig(const std::string &json_path) { evc::load_daq_config(json_path, daq_cfg_); }

    // Load the merged HyCal map that the (crate, slot, ch) lookups below use.
    // The "t" field of each record ("PbGlass" / "PbWO4" / "Veto" / "LMS") is
    // the single source of truth for category dispatch; entries without a
    // "daq" block are not reachable by DAQ address.  Without a map every
    // channel is unknown (MOD_UNKNOWN, module_id -1) and gets dropped.
    void LoadHyCalMap(const std::string &json_path);

    // Module name of the channel, "" if it is not in the map.
    std::string moduleName(int roc, int slot, int ch) const;
    // prad2::ModuleType of the channel, MOD_UNKNOWN if it is not in the map.
    prad2::ModuleType moduleType(int roc, int slot, int ch) const;
    // Globally-unique module_id, see RawEventData in EventData.h; -1 if unknown.
    int moduleID(int roc, int slot, int ch) const;
    // Reverse lookup: returns the (crate, slot, ch) tuple for a given module_id, or {-1,-1,-1}.
    std::tuple<int, int, int> moduleLocation(int module_id) const;

    // Convert an EVIO file to a ROOT file with a TTree.
    // max_events <= 0 means process all. write_peaks adds peak branches.
    bool Process(const std::string &input_evio, const std::string &output_root, RunConfig &gRunConfig,
                 const std::string &db_dir, const std::string &recon_config_file,
                 int max_events = -1, bool write_peaks = false, const std::string &daq_config_file = "",
                 const float zerosup_override = 5.f, bool Ecalib = false, bool noWaveform = false);

    bool ProcessWithRecon(const std::string &input_evio, const std::string &output_root, RunConfig &gRunConfig,
                            const std::string &db_dir, const std::string &recon_config_file,
                            const std::string &daq_config_file = "",
                            const std::string &gem_ped_file = "", float zerosup_override = 5.f,
                            bool prad1 = false, bool x17 = false, bool x17_blind = false, bool random = false, bool gem_hit = false);
    
    bool ProcessRaw2Recon(const std::string &input_raw, const std::string &output_root, RunConfig &gRunConfig,
                            const std::string &db_dir, const std::string &recon_config_file,
                            const std::string &daq_config_file = "",
                            const std::string &gem_ped_file = "",
                            bool x17 = false, bool x17_blind = false, bool random = false, bool gem_hit = false);
    
    bool Process_LMSgainFactor(const std::string &input_evio, const std::string &output_root,
                                const std::string &db_dir, const std::string &daq_config_file);

private:
    // Module lookups of moduleName/Type/ID/Location, filled by LoadHyCalMap().
    fdec::HyCalSystem hycal_map_;
    evc::DaqConfig daq_cfg_;
};

// ── Re-processing raw replay trees ──────────────────────────────────────────

// Re-derive npeaks and peak_height/time/integral of every PbWO4 channel from
// its stored samples, for raw trees written without the peak branches.  No
// module time offset is applied.
void FillPeaksFromWaveforms(prad2::RawEventData &ev, const fdec::HyCalSystem &hycal,
                            const fdec::WaveAnalyzer &ana, fdec::WaveResult &wres);

// Re-run GEM clustering and X/Y matching, with gem_sys's per-detector
// configs, on the strip hits stored in a raw tree (pedestal, common mode and
// zero suppression already applied).  hits receives the 2D hits of all
// detectors in detector order; plane_clusters, when given, the kept clusters
// as [det][0 = X, 1 = Y].  Strips of an unknown detector or plane are skipped.
void ReconstructGemStrips(const prad2::RawEventData &ev, const gem::GemSystem &gem_sys,
                          gem::GemCluster &clusterer, std::vector<gem::GEMHit> &hits,
                          std::vector<std::array<std::vector<gem::StripCluster>, 2>>
                              *plane_clusters = nullptr);

} // namespace analysis
