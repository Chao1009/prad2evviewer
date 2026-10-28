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

// Copy the 16+16+16+16+2 counters out of the payload at `off` (the caller
// bounds-checks it), validate that the ref pair looks sane (ungated ≥ gated,
// ungated > 0).  Slot is left untouched — caller fills it from whichever
// header word it parsed.
bool fill_counters(const uint32_t *data, size_t off, DscEventData &s)
{
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
    // The two layouts we know about: 0 = legacy 0xDCA0 header, 2 = BLKHDR-wrapped.
    return ParsePayload(data, nwords, out, {0, 2});
}

bool Dsc2Decoder::ParsePayload(const uint32_t *data, size_t nwords, DscEventData &out,
                               std::initializer_list<size_t> probe_offsets)
{
    out.clear();
    if (data == nullptr || nwords == 0) return false;

    // The first offset whose ref pair looks sane wins.
    for (size_t off : probe_offsets) {
        if (off + (size_t)kPayloadW > nwords) continue;
        const uint32_t hdr = data[off];

        int slot;
        if ((hdr & kHdrMask) == kHdrId)
            slot = (hdr >> 8) & 0xFF;
        else if (off >= 1 && (data[0] >> 27) == kBlkHdr)
            slot = (data[0] >> 22) & 0x1F;
        else
            continue;

        if (!fill_counters(data, off, out)) continue;
        out.slot    = slot;
        out.present = true;
        return true;
    }
    out.clear();
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
    case DSrc::Tdc: {
        if (cfg.channel < 0 || cfg.channel >= DSC2_NCH) {
            out.clear();
            return false;
        }
        const bool trg = cfg.source == DSrc::Trg;
        out.gated   = (trg ? out.trg_gated   : out.tdc_gated)[cfg.channel];
        out.ungated = (trg ? out.trg_ungated : out.tdc_ungated)[cfg.channel];
        break;
    }
    }
    return true;
}

} // namespace dsc
