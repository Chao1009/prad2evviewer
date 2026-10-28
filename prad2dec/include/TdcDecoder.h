#pragma once
//=============================================================================
// TdcDecoder.h — decode 0xE107 V1190 TDC Data banks
//
// Hit-word layout: see TdcData.h.  Produced by rol2.c from the raw 0xE10B
// v1190/v1290 hardware stream.
//
// rol2.c strips the v1190 global header, TDC headers, TDC EOB markers and
// global trailer before writing 0xE107, so no framing words are present in
// the payload — it is a flat array of hits.
//=============================================================================

#include "TdcData.h"

#include <vector>

namespace tdc
{

class TdcDecoder
{
public:
    // Decode one ROC's 0xE107 payload. Appended hits carry roc_tag for
    // identification across multi-crate events.
    //
    // Returns the number of hits appended. Stops early when TdcEventData
    // reaches MAX_TDC_HITS.
    static int DecodeRoc(const uint32_t *data, size_t nwords,
                         uint32_t roc_tag, TdcEventData &evt);

    // Decode the replay tree's flat-triple-of-vectors representation
    // (ev.tdc_roc_tags / ev.tdc_nwords / ev.tdc_words from
    // RawEventData / ReconEventData) into TdcEventData. `evt` is cleared
    // first. Returns the total number of hits appended.
    //
    // This is the canonical way for analysis code to access the TDC bank
    // from a replayed ROOT file — no bit shifts in user code.
    static int DecodeReplay(const std::vector<uint32_t> &roc_tags,
                            const std::vector<uint32_t> &nwords,
                            const std::vector<uint32_t> &words,
                            TdcEventData &evt);
};

// --- RF time extractor -----------------------------------------------------
// Two static convenience entry points for getting RfTimeData (per-event,
// per-channel ns arrays for the PRad-II RF reference). Both clear `out`
// first, drop trailing-edge hits, and silently truncate at
// RfTimeData::MAX_HITS_PER_CH (the observed max is ~6).
class RfTimeDecoder
{
public:
    // Extract RF hits from an already-decoded TdcEventData (e.g. produced
    // by EvChannel::Tdc() during live decoding).
    static void Extract(const TdcEventData &all, RfTimeData &out,
                        uint32_t rf_roc_tag = RF_ROC_TAG,
                        uint8_t  rf_slot    = RF_SLOT,
                        uint8_t  rf_ch_a    = RF_CH_A,
                        uint8_t  rf_ch_b    = RF_CH_B);

    // Decode directly from the replay tree's flat-triple representation.
    // Equivalent to DecodeReplay() into a temporary TdcEventData followed
    // by Extract(), but without the heavyweight intermediate buffer.
    static void DecodeReplay(const std::vector<uint32_t> &roc_tags,
                             const std::vector<uint32_t> &nwords,
                             const std::vector<uint32_t> &words,
                             RfTimeData &out,
                             uint32_t rf_roc_tag = RF_ROC_TAG,
                             uint8_t  rf_slot    = RF_SLOT,
                             uint8_t  rf_ch_a    = RF_CH_A,
                             uint8_t  rf_ch_b    = RF_CH_B);
};

} // namespace tdc
