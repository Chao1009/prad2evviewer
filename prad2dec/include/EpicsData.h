#pragma once
//=============================================================================
// EpicsData.h — single-event EPICS slow-control record
//
// EPICS events (top-level tag 0x001F) carry a string bank (0xE114) of
// "value channel_name\n" lines plus a 0xE112 HEAD bank with the absolute
// unix_time and a monotonic sync_counter (see SyncData.h).
//
// `EpicsRecord` is the per-event POD: parser output (sparse parallel
// name+value vectors) plus the absolute-time / event-anchor metadata.
// It is what EvChannel::Epics() returns and what offline tools write into
// the "epics" TTree (one row per EPICS event).
//
// For the run-scoped accumulator (channel registry, persistent values
// across snapshots, O(log N) lookup by event_number) see EpicsStore.h.
//=============================================================================

#include <cstdint>
#include <string>
#include <vector>

namespace epics
{

struct EpicsRecord {
    bool        present                 = false;
    uint32_t    unix_time               = 0;   // 0xE112 HEAD d[3]; absolute
    uint32_t    sync_counter            = 0;   // 0xE112 HEAD d[2]
    uint32_t    run_number              = 0;   // 0xE112 HEAD d[1]
    int32_t     event_number_at_arrival = -1;  // physics event_number seen
                                               // most recently before this
                                               // EPICS event (-1 if none yet)
    uint64_t    timestamp_at_arrival    = 0;   // TI 48-bit tick of the same
                                               // physics event (0 if none
                                               // yet); see EvChannel::
                                               // GetLastPhysicsTimestamp()

    // Channel readings — parallel arrays so consumers can dump them straight
    // into a TTree without per-row std::pair overhead.  Both must be the
    // same length on every populated record.
    std::vector<std::string> channel;
    std::vector<double>      value;

    void clear()
    {
        present                 = false;
        unix_time               = 0;
        sync_counter            = 0;
        run_number              = 0;
        event_number_at_arrival = -1;
        timestamp_at_arrival    = 0;
        channel.clear();
        value.clear();
    }
};

// Parse an EPICS text payload (one line per channel, format
// "value channel_name") into the channel/value arrays of `out`.  Trims
// whitespace and skips empty / unparsable lines.  Returns the number of
// (channel, value) pairs produced.
int ParseEpicsText(const std::string &text, EpicsRecord &out);

} // namespace epics
