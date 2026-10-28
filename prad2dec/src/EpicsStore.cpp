//=============================================================================
// EpicsStore.cpp — EPICS snapshot accumulation and lookup (see EpicsStore.h)
//=============================================================================

#include "EpicsStore.h"
#include "EpicsData.h"

#include <algorithm>

namespace epics
{

// --- Channel management -----------------------------------------------------

int EpicsStore::get_or_create_channel(const std::string &name)
{
    auto it = channel_map_.find(name);
    if (it != channel_map_.end())
        return it->second;

    int id = static_cast<int>(channel_names_.size());
    channel_names_.push_back(name);
    channel_map_[name] = id;
    return id;
}

int EpicsStore::GetChannelId(const std::string &name) const
{
    auto it = channel_map_.find(name);
    return (it != channel_map_.end()) ? it->second : -1;
}

// --- Feed — parse text and store snapshot -----------------------------------

void EpicsStore::Feed(int32_t event_number, uint64_t timestamp, const std::string &text)
{
    if (text.empty()) return;

    // Parse via the shared decoder (also used by EvChannel::Epics() and the
    // offline replay's epics tree).
    EpicsRecord rec;
    ParseEpicsText(text, rec);
    if (rec.channel.empty()) return;

    // Carry forward previous channel values — slow EPICS channels only
    // update on a subset of events, and consumers expect the last-known
    // reading to persist.
    Snapshot snap;
    snap.event_number = event_number;
    snap.timestamp    = timestamp;
    if (!snapshots_.empty())
        snap.values = snapshots_.back().values;

    for (size_t i = 0; i < rec.channel.size(); ++i) {
        int id = get_or_create_channel(rec.channel[i]);
        if (id >= static_cast<int>(snap.values.size()))
            snap.values.resize(id + 1, 0.f);
        snap.values[id] = static_cast<float>(rec.value[i]);
    }

    if (snap.values.size() < channel_names_.size())
        snap.values.resize(channel_names_.size(), 0.f);

    snapshots_.push_back(std::move(snap));
}

// --- Lookup — find most recent snapshot at or before event_number -----------

const EpicsStore::Snapshot *EpicsStore::FindSnapshot(int32_t event_number) const
{
    if (snapshots_.empty()) return nullptr;

    // snapshots are in insertion order (monotonically increasing event numbers)
    // binary search for the last snapshot with event_number <= query
    auto it = std::upper_bound(
        snapshots_.begin(), snapshots_.end(), event_number,
        [](int32_t ev, const Snapshot &s) { return ev < s.event_number; }
    );

    if (it == snapshots_.begin()) return nullptr;
    --it;
    return &(*it);
}

bool EpicsStore::GetValue(int32_t event_number, const std::string &channel, float &value) const
{
    int id = GetChannelId(channel);
    if (id < 0) return false;

    const Snapshot *snap = FindSnapshot(event_number);
    if (!snap) return false;
    if (id >= static_cast<int>(snap->values.size())) return false;

    value = snap->values[id];
    return true;
}

// --- Trim / Clear -----------------------------------------------------------

void EpicsStore::Trim(int max_count)
{
    if (max_count <= 0 || static_cast<int>(snapshots_.size()) <= max_count) return;
    int to_remove = static_cast<int>(snapshots_.size()) - max_count;
    snapshots_.erase(snapshots_.begin(), snapshots_.begin() + to_remove);
}

void EpicsStore::Clear()
{
    channel_names_.clear();
    channel_map_.clear();
    snapshots_.clear();
}

} // namespace epics
