#pragma once
//=============================================================================
// GemPedestal.h — per-strip pedestal accumulator for GEM (APV25 front-end).
//
// Algorithm (per event, per APV):
//   1) Per time sample — sort 128 strip ADCs, drop the top CM_DISCARD and
//      bottom CM_DISCARD, average the middle 72 → common_mode[ts].
//      Subtract common_mode[ts] from every strip's value.
//   2) Per strip — average the CM-corrected values across the 6 time
//      samples → contribution.  Accumulate into per-strip mean/RMS.
//
// After all events, Write() writes the text format GemSystem::LoadPedestals
// reads: per APV an "APV <crate> -1 <mpd> <adc>" header (the MPD slot is not
// known here and LoadPedestals ignores it) followed by 128
// "<strip> <offset> <noise>" lines, offset = mean and noise = RMS.
//
// Usage:
//   gem::GemPedestal ped;
//   while (read_next_event(ssp_evt))
//       ped.Accumulate(ssp_evt);        // returns the number of APVs folded
//   if (ped.NumStrips() > 0) ped.Write("gem_ped.txt");
//=============================================================================

#include <memory>
#include <string>

namespace ssp { struct SspEventData; }

namespace gem {

class GemPedestal
{
public:
    GemPedestal();
    ~GemPedestal();

    GemPedestal(const GemPedestal &)            = delete;
    GemPedestal &operator=(const GemPedestal &) = delete;

    void Clear();

    // Fold one event's SSP data into the running accumulators and return
    // the number of APVs folded.  Only full-readout APVs (nstrips == 128)
    // contribute: online-ZS strips are already pedestal/CM-subtracted by
    // the firmware.  Returns 0 for a pure online-ZS event.
    int Accumulate(const ssp::SspEventData &evt);

    // Number of APVs that received at least one contribution.
    int NumApvs() const;
    // Number of strips (across all APVs) with at least one contribution.
    int NumStrips() const;

    // Write the accumulated mean/RMS (offset to 0.001, noise to 0.0001).
    // Returns the number of APVs written, or a negative value on I/O failure.
    int Write(const std::string &output_path) const;

private:
    struct Impl;
    std::unique_ptr<Impl> impl_;
};

} // namespace gem
