#pragma once
//=============================================================================
// Dsc2Decoder.h — decode 0xE115 DSC2 scaler banks
//
// Single entry point: DecodeBank() probes the supported bank layouts (legacy
// 0xDCA0-magic vs PRad-II rflag=1 BLKHDR-wrapped), then unpacks all 16 TRG +
// 16 TDC + 2 Ref counters into DscEventData and applies the (source, channel)
// selection from the DSC2 config to populate the .gated / .ungated pair.
//
// The bank-format and live-time conventions are documented in DscData.h.
//=============================================================================

#include <initializer_list>
#include "DscData.h"
#include "DaqConfig.h"

namespace dsc
{

class Dsc2Decoder
{
public:
    // Parse one DSC2 0xE115 bank.  Returns true (and fills out.present) when
    // the data matches a known layout AND the slot matches cfg.slot; returns
    // false otherwise, with out cleared unless the config is disabled
    // (cfg.slot < 0 or bank_tag < 0), in which case out is left untouched.
    static bool DecodeBank(const uint32_t *data, size_t nwords,
                           const evc::DaqConfig::DscScaler &cfg,
                           DscEventData &out);

    // Lower-level: parse the bank without applying a (source, channel)
    // selection.  Returns true if a payload was found, populating slot/offset
    // and the per-channel + ref arrays.  out.gated / out.ungated stay at 0.
    // out is cleared first, even when nothing parses, so parse into a
    // temporary to keep a previous snapshot.  Useful for diagnostic tools
    // that want the full counter set.
    static bool ParsePayload(const uint32_t *data, size_t nwords,
                             DscEventData &out);

    // Same, probing the given payload offsets instead of the decoder's {0, 2}.
    static bool ParsePayload(const uint32_t *data, size_t nwords,
                             DscEventData &out,
                             std::initializer_list<size_t> probe_offsets);
};

} // namespace dsc
