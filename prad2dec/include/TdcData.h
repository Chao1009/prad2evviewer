#pragma once
//=============================================================================
// TdcData.h — pre-allocated event data for V1190 TDC readout (0xE107)
//
// The 0xE107 bank is emitted by rol2.c after reformatting the raw v1190/v1290
// hardware stream. Each 32-bit word is a single TDC hit packed as:
//
//   bits 31:27 — slot       (5 bits, 0-31)
//   bit  26    — edge       (0 = leading, 1 = trailing)
//   bits 25:19 — channel    (7 bits, 0-127)
//   bits 18:00 — TDC value  (19 bits)
//
// For V1190 boards, rol2.c left-shifts the native 19-bit value by 2 to match
// the v1290's 25 ps LSB, truncating the result back to 19 bits. V1290 data
// enters the 0xE107 stream already in this layout with a different native
// channel width (5 bits) — the XML dictionary documents only the V1190 form
// used in PRad-II.
//
// Tagger readout: 1 VME crate with 3 V1190 boards × 128 channels = 384 ch.
//=============================================================================

#include <cstdint>
#include <cstddef>

namespace tdc
{

// --- capacity limits --------------------------------------------------------
// 3 boards × 128 channels × worst-case multi-hit burst × 2 edges.
// 4096 comfortably covers a saturated tagger event without heap growth.
static constexpr int MAX_TDC_HITS = 4096;

// Optional per-slot/channel index for fast lookup.
static constexpr int MAX_TDC_SLOTS    = 32;   // V1190 slot field is 5 bits
static constexpr int MAX_TDC_CHANNELS = 128;  // V1190 has 128 channels/board

// --- one TDC hit ------------------------------------------------------------
struct TdcHit
{
    uint32_t roc_tag;   // parent ROC bank tag (e.g. 0x008E for the tagger crate)
    uint8_t  slot;      // 5-bit slot, V1190 board position in the VME crate
    uint8_t  channel;   // 7-bit channel, 0-127
    uint8_t  edge;      // 0 = leading, 1 = trailing
    uint32_t value;     // 19-bit TDC value (LSB = 25 ps after rol2 shift)
};

// --- full event data --------------------------------------------------------
struct TdcEventData
{
    int    n_hits = 0;
    TdcHit hits[MAX_TDC_HITS];

    void clear() { n_hits = 0; }

    // Count hits for a given slot (linear scan — fine for tagger-sized events).
    int countSlot(uint8_t slot) const
    {
        int n = 0;
        for (int i = 0; i < n_hits; ++i)
            if (hits[i].slot == slot) ++n;
        return n;
    }

    // Count hits for a given (slot, channel).
    int countChannel(uint8_t slot, uint8_t channel) const
    {
        int n = 0;
        for (int i = 0; i < n_hits; ++i)
            if (hits[i].slot == slot && hits[i].channel == channel) ++n;
        return n;
    }
};

} // namespace tdc
