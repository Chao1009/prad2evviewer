#pragma once
//=============================================================================
// DaqKey.h — pack a three-level DAQ address into one map key
//
// (crate, slot, channel) for FADC/ADC1881M, (crate, mpd, apv) for SSP/MPD,
// (roc_tag, slot, channel) for pulse templates.  Each field keeps its low
// 16 bits:  key = a << 32 | b << 16 | c.
//=============================================================================

#include <cstdint>

namespace prad2
{

constexpr uint64_t pack_daq_key(int a, int b, int c)
{
    return (static_cast<uint64_t>(static_cast<uint16_t>(a)) << 32) |
           (static_cast<uint64_t>(static_cast<uint16_t>(b)) << 16) |
           static_cast<uint64_t>(static_cast<uint16_t>(c));
}

} // namespace prad2
