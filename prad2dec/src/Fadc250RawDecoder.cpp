#include "Fadc250RawDecoder.h"
#include <iostream>

using namespace fdec;

//=============================================================================
// FADC250 hardware-format words (JLab FADC250 firmware / rol1.c)
//
// Defining words have bit 31 = 1 and the type in bits 30:27, so
// (w >> 27) = 0x10 + type:
//   0x10: Block Header    — slot 26:22
//   0x11: Block Trailer
//   0x12: Event Header    — trigger# 21:0
//   0x13: Trigger Time    — time_low 23:0; continuation: time_high 23:0
//   0x14: Window Raw Data — channel 26:23, window width 11:0,
//         then (width+1)/2 sample words, 2 ADC samples each
//   0x1E: Data Not Valid, 0x1F: Filler — skip
//
// Continuation words have bit 31 = 0.  Sample words:
//   Bits 28:16 = ADC sample i   (13 bits, 12-bit value + valid flag at bit 29)
//   Bits 12:0  = ADC sample i+1 (13 bits, 12-bit value + valid flag at bit 13)
//=============================================================================

namespace {

enum : uint32_t {
    T_BLOCK_HEADER  = 0x10,
    T_BLOCK_TRAILER = 0x11,
    T_EVENT_HEADER  = 0x12,
    T_TRIGGER_TIME  = 0x13,
    T_WINDOW_RAW    = 0x14,
    T_FILLER        = 0x1F,
};

inline uint32_t type_tag(uint32_t w) { return (w >> 27) & 0x1F; }

// Block Header
inline uint32_t bh_slot(uint32_t w)     { return (w >> 22) & 0x1F; }

// Event Header
inline uint32_t eh_trigger(uint32_t w)  { return w & 0x003FFFFF; }

// Trigger Time
inline uint32_t tt_time(uint32_t w)     { return w & 0x00FFFFFF; }

// Window Raw Data header
inline uint32_t wr_channel(uint32_t w)  { return (w >> 23) & 0x0F; }
inline uint32_t wr_width(uint32_t w)    { return w & 0x0FFF; }

// Sample extraction (13-bit ADC values, 2 per word)
inline uint16_t sample_hi(uint32_t w)   { return (w >> 16) & 0x1FFF; }
inline uint16_t sample_lo(uint32_t w)   { return w & 0x1FFF; }

} // anonymous namespace

int Fadc250RawDecoder::DecodeRoc(const uint32_t *data, size_t nwords, RocData &roc)
{
    if (!data || nwords == 0) return 0;

    int nslots = 0;
    SlotData *sd = nullptr;

    for (size_t i = 0; i < nwords; ++i) {
        uint32_t w = data[i];
        uint32_t tt = type_tag(w);

        switch (tt) {

        case T_BLOCK_HEADER: {
            uint32_t slot_id = bh_slot(w);
            if (slot_id >= MAX_SLOTS) {
                std::cerr << "Fadc250RawDecoder: slot_id=" << slot_id
                          << " >= MAX_SLOTS\n";
                sd = nullptr;
                break;
            }
            sd = &roc.slots[slot_id];
            sd->present = true;
            sd->nchannels = 0;
            sd->channel_mask = 0;
            nslots++;
            break;
        }

        case T_BLOCK_TRAILER: {
            sd = nullptr;
            break;
        }

        case T_EVENT_HEADER: {
            if (!sd) break;
            sd->trigger = static_cast<int32_t>(eh_trigger(w));
            break;
        }

        case T_TRIGGER_TIME: {
            if (!sd) break;
            uint64_t time_low = tt_time(w);
            // Continuation word has bit 31 = 0 and carries high 24 bits
            if (i + 1 < nwords && (data[i + 1] >> 31) == 0) {
                ++i;
                uint64_t time_high = tt_time(data[i]);
                sd->timestamp = static_cast<int64_t>(time_low | (time_high << 24));
            } else {
                sd->timestamp = static_cast<int64_t>(time_low);
            }
            break;
        }

        case T_WINDOW_RAW: {
            if (!sd) break;
            uint32_t ch = wr_channel(w);
            uint32_t width = wr_width(w);
            if (ch >= MAX_CHANNELS) break;

            ChannelData &cd = sd->channels[ch];
            uint32_t nsamp = 0;
            uint32_t max_samp = (width < static_cast<uint32_t>(MAX_SAMPLES))
                                ? width : static_cast<uint32_t>(MAX_SAMPLES);

            // Read (width+1)/2 continuation words, each packing 2 samples
            uint32_t nwords_expected = (width + 1) / 2;
            for (uint32_t j = 0; j < nwords_expected && i + 1 < nwords; ++j) {
                ++i;
                uint32_t sw = data[i];
                if (nsamp < max_samp) cd.samples[nsamp++] = sample_hi(sw);
                if (nsamp < max_samp) cd.samples[nsamp++] = sample_lo(sw);
            }

            cd.nsamples = static_cast<int>(nsamp);
            sd->channel_mask |= (1ull << ch);
            sd->nchannels++;
            break;
        }

        case T_FILLER:
        default:
            break;
        }
    }

    roc.nslots = nslots;
    return nslots;
}
