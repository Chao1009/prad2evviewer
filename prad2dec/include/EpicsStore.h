#pragma once
//=============================================================================
// EpicsStore.h — run-scoped EPICS snapshot accumulator with channel registry
//
// Companion to EpicsData.h.  EpicsRecord is the per-event POD; EpicsStore
// is the multi-event store the monitor server queries to answer
// "what's the most recent value of channel X at event N?"
//
// What it adds beyond a list of EpicsRecords:
//   * channel registry — a stable string→int id assigned the first time a
//     channel is seen, so later snapshots store dense `vector<float>`
//     instead of repeating channel names per record;
//   * value persistence — every Feed() carries forward the previous
//     snapshot's values, so slow channels (which only update on a subset
//     of EPICS events) keep their last-known reading;
//   * O(log N) lookup by event_number (snapshots are kept in arrival
//     order, which matches monotonic event_number).
//
// Both EpicsStore::Feed() and EvChannel::Epics() route raw text through
// `epics::ParseEpicsText` (EpicsData.h) — there is one parser.
//
// Usage:
//   EpicsStore epics;
//   // in event loop:
//   if (ch.GetEventType() == EventType::Epics)
//       epics.Feed(event_number, timestamp, ch.ExtractEpicsText());
//   // later, for any physics event:
//   float beam_current;
//   if (epics.GetValue(event_number, "beam_current", beam_current)) { ... }
//=============================================================================

#include <string>
#include <vector>
#include <unordered_map>
#include <cstdint>

namespace epics
{

class EpicsStore
{
public:
    EpicsStore() = default;

    // --- feeding data -------------------------------------------------------

    // Parse raw EPICS text and store a snapshot.
    // Text format: one "value  channel_name" pair per line.
    // event_number: the trigger/event number at the time of this EPICS update.
    // timestamp: 48-bit TI timestamp (0 if unavailable).
    void Feed(int32_t event_number, uint64_t timestamp, const std::string &text);

    // --- querying -----------------------------------------------------------

    // Get the most recent value of a channel at or before the given event number.
    // Returns true if found, false if channel unknown or no snapshot before this event.
    bool GetValue(int32_t event_number, const std::string &channel, float &value) const;

    struct Snapshot {
        int32_t              event_number;
        uint64_t             timestamp;
        std::vector<float>   values;     // indexed by channel id
    };

    // Most recent snapshot at or before event_number, or nullptr if none.
    // Use GetChannelId() to map names to indices into its values.
    const Snapshot *FindSnapshot(int32_t event_number) const;

    // --- channel info -------------------------------------------------------

    int  GetChannelCount() const { return static_cast<int>(channel_names_.size()); }
    int  GetChannelId(const std::string &name) const;
    const std::string &GetChannelName(int id) const { return channel_names_[id]; }

    const std::vector<std::string> &GetChannelNames() const { return channel_names_; }

    int  GetSnapshotCount() const { return static_cast<int>(snapshots_.size()); }

    // direct snapshot access by index (0 = oldest)
    const Snapshot &GetSnapshot(int index) const { return snapshots_[index]; }

    // trim oldest snapshots to keep at most max_count
    void Trim(int max_count);

    // --- reset --------------------------------------------------------------
    void Clear();

private:
    int get_or_create_channel(const std::string &name);

    std::vector<std::string>                   channel_names_;
    std::unordered_map<std::string, int>       channel_map_;    // name → id
    std::vector<Snapshot>                      snapshots_;      // sorted by event_number
};

} // namespace epics
