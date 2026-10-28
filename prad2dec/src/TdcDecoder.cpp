#include "TdcDecoder.h"
#include "EvStruct.h"

#include <cmath>
#include <limits>

using namespace tdc;

namespace
{

TdcHit decode_hit(uint32_t w, uint32_t roc_tag)
{
    TdcHit h;
    h.roc_tag = roc_tag;
    h.slot    = static_cast<uint8_t>((w >> 27) & 0x1F);
    h.edge    = static_cast<uint8_t>((w >> 26) & 0x1);
    h.channel = static_cast<uint8_t>((w >> 19) & 0x7F);
    h.value   = w & 0x7FFFF;
    return h;
}

void push_rf(RfTimeData &out, const TdcHit &h,
             uint8_t rf_slot, uint8_t rf_ch_a, uint8_t rf_ch_b)
{
    if (h.slot != rf_slot || h.edge != 0) return;   // leading edges only
    float t_ns = static_cast<float>(h.value * TDC_LSB_NS);
    if (h.channel == rf_ch_a) {
        if (out.n_a < RfTimeData::MAX_HITS_PER_CH)
            out.ns_a[out.n_a++] = t_ns;
    } else if (h.channel == rf_ch_b) {
        if (out.n_b < RfTimeData::MAX_HITS_PER_CH)
            out.ns_b[out.n_b++] = t_ns;
    }
}

} // namespace

int TdcDecoder::DecodeRoc(const uint32_t *data, size_t nwords,
                          uint32_t roc_tag, TdcEventData &evt)
{
    int appended = 0;
    for (size_t i = 0; i < nwords; ++i) {
        if (evt.n_hits >= MAX_TDC_HITS) break;
        evt.hits[evt.n_hits++] = decode_hit(data[i], roc_tag);
        ++appended;
    }
    return appended;
}

int TdcDecoder::DecodeReplay(const std::vector<uint32_t> &roc_tags,
                             const std::vector<uint32_t> &nwords,
                             const std::vector<uint32_t> &words,
                             TdcEventData &evt)
{
    evt.clear();
    int total = 0;
    evc::ForEachFlatBank(roc_tags, nwords, words,
        [&](uint32_t roc, const uint32_t *d, size_t n) {
            total += DecodeRoc(d, n, roc, evt);
        });
    return total;
}

// --- RfTimeData --------------------------------------------------------------

float RfTimeData::nearest(const float *v, int n, float t_ref)
{
    if (n <= 0) return std::numeric_limits<float>::quiet_NaN();
    int best_i = 0;
    float best_dt = std::abs(v[0] - t_ref);
    for (int i = 1; i < n; ++i) {
        float dt = std::abs(v[i] - t_ref);
        if (dt < best_dt) { best_dt = dt; best_i = i; }
    }
    return v[best_i];
}

// --- RfTimeDecoder -----------------------------------------------------------

void RfTimeDecoder::Extract(const TdcEventData &all, RfTimeData &out,
                            uint32_t rf_roc_tag, uint8_t rf_slot,
                            uint8_t rf_ch_a,    uint8_t rf_ch_b)
{
    out.clear();
    for (int i = 0; i < all.n_hits; ++i) {
        const TdcHit &h = all.hits[i];
        if (h.roc_tag == rf_roc_tag)
            push_rf(out, h, rf_slot, rf_ch_a, rf_ch_b);
    }
}

void RfTimeDecoder::DecodeReplay(const std::vector<uint32_t> &roc_tags,
                                 const std::vector<uint32_t> &nwords,
                                 const std::vector<uint32_t> &words,
                                 RfTimeData &out,
                                 uint32_t rf_roc_tag, uint8_t rf_slot,
                                 uint8_t rf_ch_a,    uint8_t rf_ch_b)
{
    out.clear();
    evc::ForEachFlatBank(roc_tags, nwords, words,
        [&](uint32_t roc, const uint32_t *d, size_t n) {
            if (roc != rf_roc_tag) return;
            for (size_t k = 0; k < n; ++k)
                push_rf(out, decode_hit(d[k], roc), rf_slot, rf_ch_a, rf_ch_b);
        });
}
