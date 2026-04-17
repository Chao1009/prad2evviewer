#include "TdcDecoder.h"

using namespace tdc;

int TdcDecoder::DecodeRoc(const uint32_t *data, size_t nwords,
                          uint32_t roc_tag, TdcEventData &evt)
{
    int appended = 0;
    for (size_t i = 0; i < nwords; ++i) {
        if (evt.n_hits >= MAX_TDC_HITS) break;
        uint32_t w = data[i];
        TdcHit &h = evt.hits[evt.n_hits++];
        h.roc_tag = roc_tag;
        h.slot    = static_cast<uint8_t>((w >> 27) & 0x1F);
        h.edge    = static_cast<uint8_t>((w >> 26) & 0x1);
        h.channel = static_cast<uint8_t>((w >> 19) & 0x7F);
        h.value   = w & 0x7FFFF;
        ++appended;
    }
    return appended;
}
