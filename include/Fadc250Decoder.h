#pragma once

//
// A simple decoder to process the JLab FADC250 data
// Reference: https://www.jlab.org/Hall-B/ftof/manuals/FADC250UsersManual.pdf
//
// Author: Chao Peng
// Date: 2020/08/22
//
// Updated: 2025 - Added composite data type (tag 0xe126) decoding
//

#include <iostream>
#include <iomanip>
#include <map>
#include <vector>
#include <cstring>
#include "Fadc250Data.h"


namespace fdec
{

// data type
enum Fadc250Type {
    BlockHeader = 0,
    BlockTrailer = 1,
    EventHeader = 2,
    TriggerTime = 3,
    WindowRawData = 4,
    WindowSum = 5,
    PulseRawData = 6,
    PulseIntegral = 7,
    PulseTime = 8,
    StreamingRawData = 9,
    // 10 - 11 user defined
    Scaler = 12,
    EventTrailer = 13,
    InvalidData = 14,
    FillerWord = 15,
};

class Fadc250Event
{
public:
    uint32_t number, mode;
    std::vector<uint32_t> time;
    std::vector<Fadc250Data> channels;

    Fadc250Event(uint32_t n = 0, uint32_t nch = 16)
        : number(n), mode(0)
    {
        channels.resize(nch);
    }

    void Clear()
    {
        mode = 0;
        time.clear();
        for (auto &ch : channels) { ch.Clear(); }
    }
};

// Result of decoding a composite bank: one entry per (slot, channel)
struct CompositeHit
{
    uint8_t slot;
    uint8_t channel;
    std::vector<uint16_t> samples;
};

class Fadc250Decoder
{
public:
    Fadc250Decoder(double clk = 250.);

    // ---------------------------------------------------------------
    //  Legacy: decode raw FADC250 module words (entangled/disentangled)
    // ---------------------------------------------------------------
    void DecodeEvent(Fadc250Event &event, const uint32_t *buf, size_t len) const;
    inline Fadc250Event DecodeEvent(const uint32_t *buf, size_t len, size_t nchans = 16) const
    {
        Fadc250Event evt;
        evt.channels.resize(nchans);
        DecodeEvent(evt, buf, len);
        return evt;
    }

    // ---------------------------------------------------------------
    //  Composite (tag 0xe126): "FADC250 Window Raw Data (mode 1 packed/compressed)"
    //
    //  Format string: c,m(c,ms)
    //      c     -> slot number       (uint8)
    //      m     -> nChannelsFired     (uint8, loop count)
    //        c   -> channel number    (uint8)
    //        m   -> nSamples          (uint16, loop count)
    //        s   -> sample values     (int16 each)
    //
    //  `data` / `nbytes` point to the raw data payload of the inner
    //  BANK inside the composite envelope (after the TagSegment + format
    //  string + inner Bank header have been stripped).
    //
    //  Returns a vector of CompositeHit, one per (slot, channel).
    // ---------------------------------------------------------------
    std::vector<CompositeHit> DecodeComposite(const uint8_t *data, size_t nbytes) const;

    // Convenience: decode composite and fill Fadc250Event per slot.
    std::map<uint8_t, Fadc250Event> DecodeCompositeToEvents(const uint8_t *data, size_t nbytes,
                                                             size_t nchans = 16) const;

private:
    double _clk;
};

}; // namespace fdec
