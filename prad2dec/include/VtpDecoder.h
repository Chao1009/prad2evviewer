#pragma once
//=============================================================================
// VtpDecoder.h — decode 0xE122 VTP Hardware Data banks
//
// Format (from docs/rols/clonbanks_20260406.xml):
//   Each record starts on a "defining" word where bit 31 = 1 and
//   bits[31:27] select the record type. Continuation words have bit 31 = 0
//   and are consumed by the active record type.
//
// Record types we parse:
//   0x10 BLKHDR      — slot, block#, block_level
//   0x11 BLKTLR      — slot, nwords   (closes the current block)
//   0x12 EVTHDR      — event_number
//   0x13 TRGTIME     — 48-bit trigger time (2 words)
//   0x14 EC_PEAK     — ECAL trigger peak (2 words)
//   0x15 EC_CLUSTER  — ECAL trigger cluster (2 words)
//   0x1C TAG_EXP     — tag-expanded record; 9-bit tag 0x1CC is PRAD_CLUSTER
//                      (HyCal trigger cluster, 2 words), other tags skipped
//   0x1D TRIGGER     — trigger time + pattern (2 words, skipped)
//   0x1E DNV         — data not valid, skipped
//   0x1F FILLER      — filler, skipped
// Other types (HTCC/FT/FTOF/CTOF/CND/PCU clusters) are stepped over but not
// stored — they are CLAS12/HPS, not PRad-II.
//=============================================================================

#include "VtpData.h"

#include <vector>

namespace vtp
{

class VtpDecoder
{
public:
    // Decode one ROC's 0xE122 payload. Entries are appended to evt.
    // roc_tag is the parent ROC bank tag, used to annotate records.
    // Returns the number of EC_PEAK, EC_CLUSTER and PRAD_CLUSTER records
    // appended, or 0 for stub banks (block header/trailer only).
    static int DecodeRoc(const uint32_t *data, size_t nwords,
                         uint32_t roc_tag, VtpEventData &evt);

    // Decode the replay tree's flat bank list (ev.vtp_roc_tags /
    // ev.vtp_nwords / ev.vtp_words) into `evt`, which is cleared first.
    // Returns false if the three vectors are inconsistent (see
    // evc::ForEachFlatBank); banks before the inconsistency are kept.
    static bool DecodeReplay(const std::vector<uint32_t> &roc_tags,
                             const std::vector<uint32_t> &nwords,
                             const std::vector<uint32_t> &words,
                             VtpEventData &evt);
};

} // namespace vtp
