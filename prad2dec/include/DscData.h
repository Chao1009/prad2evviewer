#pragma once
//=============================================================================
// DscData.h — DSC2 scaler bank decoded record
//
// One DSC2 readout = one of these.  In PRad-II run 024246 the bank is emitted
// in physics events that carry the SYNC flag (~every 1–2 s) inside the TI
// master crate (parent ROC tag 0x0027).  Counts ACCUMULATE from the GO
// transition; subtract two consecutive snapshots to get a windowed live time.
//
// Bank format — two physical layouts handled by Dsc2Decoder:
//
//   Legacy 67-word: payload at offset 0, slot in 0xDCA0... magic header word.
//   PRad-II rflag=1 (run 024246): 72 words = 3-word JLab BLKHDR/EVTHDR/TRGTIME
//     prefix, 67-word payload at offset 2, 2-word FILLER/BLKTLR trailer; slot
//     lives in BLKHDR bits 26:22.
//
// 67-word payload, per slot:
//   [0]      header (placeholder in the rflag=1 form)
//   [1..16]  TRG Grp1 (gated)   — 16 channels
//   [17..32] TDC Grp1 (gated)   — 16 channels
//   [33..48] TRG Grp2 (ungated) — 16 channels
//   [49..64] TDC Grp2 (ungated) — 16 channels
//   [65]     Ref Grp1 (gated)
//   [66]     Ref Grp2 (ungated)
//
// Convention in this DAQ: Group A (gated) is enabled while NOT busy, so it
// counts LIVE time and live_fraction = gated/ungated.
//=============================================================================

#include <cstdint>
#include <cstring>

namespace dsc
{

static constexpr int DSC2_NCH = 16;

struct DscEventData {
    bool     present  = false;       // false ⇒ no DSC2 bank in the event
    int      slot     = -1;          // physical DSC2 slot (BLKHDR or 0xDCA0 hdr)
    int      offset   = 0;           // payload offset inside the bank (0 or 2)

    // Selected source counters — what the configured (source, channel) picks.
    // gated counts during live; live_ratio = gated / ungated.
    uint32_t gated    = 0;
    uint32_t ungated  = 0;

    // Full per-channel + reference counters, kept for diagnostics and for
    // writers that want to record more than just the configured pick.
    uint32_t trg_gated[DSC2_NCH]   = {};
    uint32_t tdc_gated[DSC2_NCH]   = {};
    uint32_t trg_ungated[DSC2_NCH] = {};
    uint32_t tdc_ungated[DSC2_NCH] = {};
    uint32_t ref_gated   = 0;
    uint32_t ref_ungated = 0;

    void clear() { *this = DscEventData{}; }

    // Live fraction in [0,1] from the selected (gated, ungated) pair.
    // Returns -1 if ungated is zero (not yet populated).
    double live_ratio() const
    {
        return ungated > 0 ? (double)gated / (double)ungated : -1.0;
    }
};

} // namespace dsc
