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
//     prefix, 64-word counter array + 5-word trailer at offset 2; slot lives
//     in BLKHDR bits 26:22.
//
// 64 per-channel counters + a reference pair, organised as 4 groups of
// 16 input connectors followed by 2 reference words — group 1 + group 2
// are gated, group 3 + group 4 are ungated, and groups 1/2 (resp. 3/4)
// report nearly identical counts (the firmware writes each pair from
// the same hardware counter).  Indices below use the bank-relative
// numbering after the 3-word header skip:
//
//   [3..18]  group 1 (gated)   — 16 channels — `trg_gated[]`
//   [19..34] group 2 (gated)   — 16 channels — `tdc_gated[]`   (≈ group 1)
//   [35..50] group 3 (ungated) — 16 channels — `trg_ungated[]`
//   [51..66] group 4 (ungated) — 16 channels — `tdc_ungated[]` (≈ group 3)
//   [67]     reference (gated)                — `ref_gated`
//   [68]     reference (ungated)              — `ref_ungated`
//   [69..71] block trailer (FILLER + BLKTLR + tail) — not decoded
//
// Convention: gated counters are enabled while NOT busy, so they count
// LIVE time and live_fraction = gated/ungated.  Faraday cup readout is
// typically on channel 0, channel 1 carries an auxiliary signal, and
// the remaining 14 connectors are usually unused (cabling varies by
// run — verify with a few non-zero channels in the live data).
//
// Recommended primary source for livetime is `trg` (group 1 + group 3,
// per-channel).  The `ref` pair tracks the same livetime via the bank-
// level reference inputs and gives a similar fraction; it is a useful
// cross-check but not the canonical pick.  The `tdc_*` field names
// predate the spec where group 2/4 are duplicates of group 1/3 (not a
// separate TDC source) — kept for back-compat.
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
