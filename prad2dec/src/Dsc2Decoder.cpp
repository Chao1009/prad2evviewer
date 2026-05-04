//=============================================================================
// Dsc2Decoder.cpp — see Dsc2Decoder.h / DscData.h for layout + convention.
//=============================================================================

#include "Dsc2Decoder.h"

#include <cstring>

namespace dsc
{

namespace {

constexpr int      kPayloadW = 67;            // 67-word DSC2 payload
constexpr uint32_t kHdrMask  = 0xFFFF0000u;
constexpr uint32_t kHdrId    = 0xDCA00000u;   // legacy DSC2 header magic
constexpr uint32_t kBlkHdr   = 0x10u;         // JLab block-header type

// Copy the 16+16+16+16+2 counters out of the payload at `off`, validate that
// the ref pair looks sane (ungated ≥ gated, ungated > 0).  Slot is left
// untouched — caller fills it from whichever header word it parsed.
bool fill_counters(const uint32_t *data, size_t nwords, size_t off, DscEventData &s)
{
    if (off + (size_t)kPayloadW > nwords) return false;
    const uint32_t *p = &data[off + 1];
    std::memcpy(s.trg_gated,   p,      DSC2_NCH * sizeof(uint32_t));
    std::memcpy(s.tdc_gated,   p + 16, DSC2_NCH * sizeof(uint32_t));
    std::memcpy(s.trg_ungated, p + 32, DSC2_NCH * sizeof(uint32_t));
    std::memcpy(s.tdc_ungated, p + 48, DSC2_NCH * sizeof(uint32_t));
    s.ref_gated   = p[64];
    s.ref_ungated = p[65];
    s.offset      = static_cast<int>(off);
    return s.ref_ungated > 0 && s.ref_ungated >= s.ref_gated;
}

} // namespace

bool Dsc2Decoder::ParsePayload(const uint32_t *data, size_t nwords, DscEventData &out)
{
    out.clear();
    if (data == nullptr || nwords == 0) return false;

    // Probe the two layouts we know about.  The first one whose ref pair
    // looks sane wins.
    static const size_t kProbes[] = {0, 2};
    for (size_t off : kProbes) {
        if (off + (size_t)kPayloadW > nwords) continue;
        const uint32_t hdr = data[off];

        if ((hdr & kHdrMask) == kHdrId) {
            if (!fill_counters(data, nwords, off, out)) continue;
            out.slot    = (hdr >> 8) & 0xFF;
            out.present = true;
            return true;
        }

        if (off >= 1 && (data[0] >> 27) == kBlkHdr) {
            if (!fill_counters(data, nwords, off, out)) continue;
            out.slot    = (data[0] >> 22) & 0x1F;
            out.present = true;
            return true;
        }
    }
    return false;
}

bool Dsc2Decoder::DecodeBank(const uint32_t *data, size_t nwords,
                             const evc::DaqConfig::DscScaler &cfg,
                             DscEventData &out)
{
    if (!cfg.enabled()) return false;
    if (!ParsePayload(data, nwords, out)) return false;
    if (out.slot != cfg.slot) { out.clear(); return false; }

    using DSrc = evc::DaqConfig::DscScaler::Source;
    switch (cfg.source) {
    case DSrc::Ref:
        out.gated   = out.ref_gated;
        out.ungated = out.ref_ungated;
        break;
    case DSrc::Trg:
        if (cfg.channel < 0 || cfg.channel >= DSC2_NCH) {
            out.clear();
            return false;
        }
        out.gated   = out.trg_gated[cfg.channel];
        out.ungated = out.trg_ungated[cfg.channel];
        break;
    case DSrc::Tdc:
        if (cfg.channel < 0 || cfg.channel >= DSC2_NCH) {
            out.clear();
            return false;
        }
        out.gated   = out.tdc_gated[cfg.channel];
        out.ungated = out.tdc_ungated[cfg.channel];
        break;
    }
    return true;
}

} // namespace dsc
