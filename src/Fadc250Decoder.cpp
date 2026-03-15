//
// A simple decoder to process the JLab FADC250 data
// Reference: https://www.jlab.org/Hall-B/ftof/manuals/FADC250UsersManual.pdf
//
// Author: Chao Peng
// Date: 2020/08/22
//
// Updated: 2025 - Added composite data type (tag 0xe126) decoding
//

#include "Fadc250Decoder.h"

using namespace fdec;

#define SET_BIT(n,i)  ( (n) |= (1ULL << i) )
#define TEST_BIT(n,i)  ( (bool)( n & (1ULL << i) ) )

inline void print_word(uint32_t word)
{
    std::cout << "0x" << std::hex << std::setw(8) << std::setfill('0') << word << std::dec << "\n";
}

inline Fadc250Data &get_channel(Fadc250Event &ev, uint32_t ch)
{
    return ev.channels[ch];
}

template<class Container>
inline uint32_t fill_in_words(const uint32_t *buf, size_t beg, Container &raw_data, size_t max_words = -1)
{
    uint32_t nwords = 0;
    for (uint32_t i = beg + 1; raw_data.size() < max_words; ++i, ++nwords) {
        auto data = buf[i];
        // finished
        if ((data & 0x80000000) && nwords > 0) {
            return nwords;
        }

        if (!(data & 0x20000000)) {
            raw_data.push_back((data >> 16) & 0x1FFF);
        }
        if (!(data & 0x2000)) {
            raw_data.push_back((data & 0x1FFF));
        }
    }
    return nwords;
}


Fadc250Decoder::Fadc250Decoder(double clk)
: _clk(clk)
{
    // place holder
}

// a help structure to save peak infos
struct PeakBuffer {
    uint32_t height = 0., integral = 0., time = 0.;
    bool in_data = false;
};

// ===================================================================
//  Legacy decoder: raw FADC250 module words
// ===================================================================
void Fadc250Decoder::DecodeEvent(Fadc250Event &res, const uint32_t *buf, size_t buflen)
const
{
    res.Clear();

    // sanity check
    if (!buflen) {
        return;
    }

    auto header = buf[0];
    if (!(header & 0x80000000) || ((header >> 27) & 0xF) != EventHeader) {
        std::cout << "Fadc250Decoder Error: incorrect event header:";
        print_word(buf[0]);
        return;
    }

    res.number = (header & 0x3FFFFF);
    std::vector<std::vector<PeakBuffer>> peak_buffers(res.channels.size());
    uint32_t type = FillerWord;

    for (size_t iw = 1; iw < buflen; ++iw) {
        uint32_t data = buf[iw];

        // new type word, update the current type
        bool new_type = (data & 0x80000000);
        if (new_type) {
            type = (data >> 27) & 0xF;
            SET_BIT(res.mode, type);
        }

        switch (type) {
        // trigger timing, might be multiple timing words
        case TriggerTime:
            res.time.push_back(data & 0xFFFFFF);
            break;
        // window raw data
        case WindowRawData:
            if (new_type) {
                // get channel and window size
                uint32_t ch = (data >> 23) & 0xF;
                size_t nwords= (data & 0xFFF);
                auto &raw_data = get_channel(res, ch).raw;
                raw_data.clear();
                iw += fill_in_words(buf, iw, raw_data, nwords);
            } else {
                std::cout << "Fadc250Decoder Error: unexpected window raw data word. ";
                print_word(data);
            }
            break;
        // pulse raw data, TODO: currently unsupported
        case PulseRawData:
            std::cout << "Fadc250Decoder Warning: unsupported data mode: pulse raw data, skip it" << std::endl;
            break;
        // pulse integral
        case PulseIntegral:
            {
                uint32_t ch = (data >> 23) & 0xF;
                uint32_t pulse_num = (data >> 21) & 0x3;
                if (peak_buffers[ch].size() < pulse_num + 1) {
                    peak_buffers[ch].resize(4);
                }
                peak_buffers[ch][pulse_num].integral = data & 0x7FFFF;
                peak_buffers[ch][pulse_num].in_data = true;
            }
            break;
        case PulseTime:
            {
                uint32_t ch = (data >> 23) & 0xF;
                uint32_t pulse_num = (data >> 21) & 0x3;
                auto &chan = get_channel(res, ch);
                if (peak_buffers[ch].size() < pulse_num + 1) {
                    peak_buffers[ch].resize(4);
                }
                peak_buffers[ch][pulse_num].time = data & 0xFFFF;
                peak_buffers[ch][pulse_num].in_data = true;
            }
            break;
        case Scaler:
            break;
        case InvalidData:
        case FillerWord:
            break;
        default:
            std::cout << "Error: unexpected data type " << type << " in header processing. ";
            print_word(data);
            return;
        }
    }

    // fill peak buffers to result
    for (size_t i = 0; i < peak_buffers.size(); ++i) {
        for (auto &peak : peak_buffers[i]) {
            if (!peak.in_data) {
                continue;
            }
            res.channels[i].peaks.emplace_back(static_cast<double>(peak.height),
                                               static_cast<double>(peak.integral),
                                               static_cast<double>(peak.time)*15.625/_clk);
        }
    }

    return;
}


// ===================================================================
//  Composite decoder (tag 0xe126)
//
//  Format string: c,m(c,ms)
//
//  Byte stream layout (big-endian):
//    Repeating until end of data:
//      [1 byte]  slot number               (c)
//      [1 byte]  nChannels fired            (m)
//        For each channel:
//          [1 byte]  channel number         (c)
//          [2 bytes] nSamples (BE uint16)   (m)
//          For each sample:
//            [2 bytes] sample value (BE int16)  (s)
//      [padding to 4-byte boundary]
// ===================================================================

static inline uint16_t read_be16(const uint8_t *p)
{
    return (static_cast<uint16_t>(p[0]) << 8) | static_cast<uint16_t>(p[1]);
}

std::vector<CompositeHit>
Fadc250Decoder::DecodeComposite(const uint8_t *data, size_t nbytes) const
{
    std::vector<CompositeHit> hits;

    if (!data || nbytes == 0) {
        return hits;
    }

    size_t pos = 0;

    while (pos < nbytes) {

        // --- slot number (c = 1 byte) ---
        if (pos >= nbytes) break;
        uint8_t slot = data[pos++];

        // --- number of channels fired (m = 1 byte) ---
        if (pos >= nbytes) {
            std::cerr << "Fadc250Decoder::DecodeComposite: truncated at nChannels for slot "
                      << (int)slot << "\n";
            break;
        }
        uint8_t nChannels = data[pos++];

        // --- loop over channels ---
        for (uint8_t ich = 0; ich < nChannels; ++ich) {

            if (pos >= nbytes) {
                std::cerr << "Fadc250Decoder::DecodeComposite: truncated at channel number, slot "
                          << (int)slot << "\n";
                break;
            }
            uint8_t channel = data[pos++];

            if (pos + 2 > nbytes) {
                std::cerr << "Fadc250Decoder::DecodeComposite: truncated at nSamples, slot "
                          << (int)slot << " ch " << (int)channel << "\n";
                break;
            }
            uint16_t nSamples = read_be16(&data[pos]);
            pos += 2;

            CompositeHit hit;
            hit.slot = slot;
            hit.channel = channel;
            hit.samples.reserve(nSamples);

            if (pos + 2 * static_cast<size_t>(nSamples) > nbytes) {
                std::cerr << "Fadc250Decoder::DecodeComposite: truncated sample data, slot "
                          << (int)slot << " ch " << (int)channel
                          << " expected " << nSamples << " samples\n";
                nSamples = static_cast<uint16_t>((nbytes - pos) / 2);
            }

            for (uint16_t is = 0; is < nSamples; ++is) {
                uint16_t sample = read_be16(&data[pos]);
                pos += 2;
                hit.samples.push_back(sample);
            }

            hits.push_back(std::move(hit));
        }

        // Pad to 4-byte boundary after each slot block
        size_t remainder = pos % 4;
        if (remainder != 0) {
            pos += (4 - remainder);
        }
    }

    return hits;
}


std::map<uint8_t, Fadc250Event>
Fadc250Decoder::DecodeCompositeToEvents(const uint8_t *data, size_t nbytes, size_t nchans) const
{
    std::map<uint8_t, Fadc250Event> result;

    auto hits = DecodeComposite(data, nbytes);

    for (auto &hit : hits) {
        auto it = result.find(hit.slot);
        if (it == result.end()) {
            auto ins = result.emplace(hit.slot, Fadc250Event(0, static_cast<uint32_t>(nchans)));
            it = ins.first;
        }

        auto &evt = it->second;
        if (hit.channel >= evt.channels.size()) {
            evt.channels.resize(hit.channel + 1);
        }

        auto &ch = evt.channels[hit.channel];
        ch.raw.clear();
        ch.raw.reserve(hit.samples.size());
        for (auto s : hit.samples) {
            ch.raw.push_back(static_cast<uint32_t>(s));
        }
    }

    return result;
}
