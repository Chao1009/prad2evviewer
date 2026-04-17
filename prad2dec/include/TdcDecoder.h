#pragma once
//=============================================================================
// TdcDecoder.h — decode 0xE107 V1190 TDC Data banks
//
// Format (from docs/rols/clonbanks_20260406.xml, produced by rol2.c from the
// raw 0xE10B v1190/v1290 hardware stream):
//
//   Each 32-bit word is one TDC hit:
//     bits 31:27  SLOT     (0-31)
//     bit  26     EDGE     (0 = leading, 1 = trailing)
//     bits 25:19  CH       (0-127)
//     bits 18:00  TDC      (19 bits, LSB = 25 ps for v1190/v1290 after rol2)
//
// rol2.c strips the v1190 global header, TDC headers, TDC EOB markers and
// global trailer before writing 0xE107, so no framing words are present in
// the payload — it is a flat array of hits.
//=============================================================================

#include "TdcData.h"

namespace tdc
{

class TdcDecoder
{
public:
    // Decode one ROC's 0xE107 payload. Appended hits carry roc_tag for
    // identification across multi-crate events (we currently expect only
    // 0x008E for the tagger crate).
    //
    // Returns the number of hits appended. Stops early when TdcEventData
    // reaches MAX_TDC_HITS.
    static int DecodeRoc(const uint32_t *data, size_t nwords,
                         uint32_t roc_tag, TdcEventData &evt);
};

} // namespace tdc
