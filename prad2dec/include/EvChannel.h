#pragma once
//=============================================================================
// EvChannel.h — read evio events, scan bank tree, lazily decode per product
//
// Usage:
//   EvChannel ch;
//   ch.SetConfig(cfg);             // builds per-product tag lists
//   ch.OpenAuto("file.evio");      // RA when supported, sequential otherwise
//   while (ch.Read() == status::success) {
//       if (!ch.Scan()) continue;
//       if (ch.GetEventType() != EventType::Physics) continue;
//       for (int i = 0; i < ch.GetNEvents(); ++i) {
//           ch.SelectEvent(i);
//           const auto &info = ch.Info();   // always cheap
//           const auto &fadc = ch.Fadc();   // decoded on first call
//           // multiple requests for the same product reuse the cached result
//       }
//   }
//=============================================================================

#include "EvStruct.h"
#include "Fadc250Data.h"
#include "SspData.h"
#include "VtpData.h"
#include "TdcData.h"
#include "SyncData.h"
#include "DscData.h"
#include "EpicsData.h"
#include "DaqConfig.h"
#include <string>
#include <vector>
#include <unordered_map>
#include <memory>

namespace evc {

enum class status : int { failure = -1, success = 1, incomplete = 2, empty = 3, eof = 4 };

class EvChannel
{
public:
    EvChannel(size_t buflen = 1024 * 2000);
    virtual ~EvChannel() { Close(); }
    EvChannel(const EvChannel &) = delete;
    EvChannel &operator=(const EvChannel &) = delete;

    // --- configuration ------------------------------------------------------
    // Stores the config and precomputes per-product tag lists used by the
    // lazy accessors.  Every tag set in the bank-tag fields but missing from
    // data_banks (e.g. JSON without a bank_structure section) gets a default
    // entry, so older configs keep working.
    void SetConfig(const DaqConfig &cfg);
    const DaqConfig &GetConfig() const   { return config; }

    // Explicit sequential evio open (evio "r" mode).  Pairs symmetrically
    // with OpenRandomAccess().  For most callers OpenAuto() below is the
    // right choice — it tries random-access first and falls back to this.
    virtual status OpenSequential(const std::string &path);
    virtual void   Close();
    virtual status Read();

    // --- random-access mode (evio "ra") -------------------------------------
    // Opens the file with evio's random-access mode: the file is mmap'd and
    // an event-pointer table is built during Open.  Use ReadEventByIndex()
    // afterwards; Read() is for sequential mode only.  The two modes are
    // mutually exclusive — calling OpenRandomAccess on an already-open
    // handle closes it first.
    virtual status OpenRandomAccess(const std::string &path);

    // Total number of events in the random-access table (0 if not opened in
    // random-access mode).  Each "event" here is one evio event (one block
    // in CODA built-trigger streams); use Scan() + GetNEvents() afterwards
    // to iterate physics sub-events within the block.
    int GetRandomAccessEventCount() const { return ra_count; }

    // Copy the event at the given 0-based index into the internal buffer, so
    // Scan() / SelectEvent() / Info() / Fadc() / ... work identically to the
    // sequential path.  Returns status::failure if the index is out of range
    // or the handle wasn't opened in random-access mode.
    status ReadEventByIndex(int evio_event_index);

    // --- convenience: open with RA when possible, else sequential ---------
    // Tries OpenRandomAccess first; if that fails (e.g. evio-4 files whose
    // blocks lack the optional index array), falls back to OpenSequential().
    // After success, IsRandomAccess() tells the caller which mode was
    // actually selected so they can dispatch between ReadEventByIndex() and
    // Read().
    virtual status OpenAuto(const std::string &path);

    // True iff the current handle was opened in random-access mode.
    bool IsRandomAccess() const { return ra_count > 0; }

    // --- scan the current event into a flat tree ----------------------------
    // Rebuilds nodes[] and the tag index; invalidates the per-product cache.
    bool Scan();

    // --- event type (valid after Scan) --------------------------------------
    EventType GetEventType() const { return evtype; }

    // --- tree accessors -----------------------------------------------------
    BankHeader                  GetEvHeader() const { return BankHeader(&buffer[0]); }
    const std::vector<EvNode>  &GetNodes()    const { return nodes; }
    const EvNode               &GetChild(const EvNode &n, size_t i) const { return nodes[n.child_first + i]; }
    std::vector<const EvNode*>  FindByTag(uint32_t tag) const;

    // Find first node with given tag (no allocation). Returns nullptr if not found.
    const EvNode *FindFirstByTag(uint32_t tag) const;

    // O(1) lookup of every node index carrying a given tag in the current
    // event (populated by Scan).  Empty span if the tag is not present.
    const std::vector<int> &NodesForTag(uint32_t tag) const;

    const uint32_t *GetData(const EvNode &n) const { return &buffer[n.data_begin]; }
    const uint8_t  *GetBytes(const EvNode &n) const
    { return reinterpret_cast<const uint8_t*>(&buffer[n.data_begin]); }
    size_t GetDataBytes(const EvNode &n) const { return n.data_words * sizeof(uint32_t); }
    const uint8_t *GetCompositePayload(const EvNode &n, size_t &nbytes) const;

    // Text of a CHARSTAR8/CHAR8 bank up to its first NUL.  EVIO pads string
    // banks with NULs and may pack several strings; only the first is read.
    std::string GetString(const EvNode &n) const;

    // Visit every non-empty UINT32 bank carrying `tag` that lies inside the
    // buffer and is not part of a composite bank, in scan order, as
    // fn(const EvNode &bank, uint32_t roc_tag); roc_tag is the parent bank's
    // tag (0 at top level).
    template <class Fn>
    void ForEachLeafBank(uint32_t tag, Fn &&fn) const
    {
        for (int ni : NodesForTag(tag)) {
            const EvNode &n = nodes[ni];
            if (n.type != DATA_UINT32 || n.data_words == 0) continue;
            if (n.data_begin + n.data_words > buffer.size()) continue;
            if (n.parent >= 0 && nodes[n.parent].type == DATA_COMPOSITE) continue;
            fn(n, n.parent >= 0 ? nodes[n.parent].tag : 0u);
        }
    }

    uint32_t       *GetRawBuffer()       { return buffer.data(); }
    const uint32_t *GetRawBuffer() const { return buffer.data(); }

    // Number of events in this block (from the physics event header num field).
    // For single-event mode this is 1. For multi-event blocks this is M.
    int GetNEvents() const { return nevents; }

    // Bank tags decoded by Tdc() / Vtp(), in decode order (from SetConfig).
    const std::vector<uint32_t> &GetTdcTags() const { return tdc_tags; }
    const std::vector<uint32_t> &GetVtpTags() const { return vtp_tags; }

    // --- lazy data-product accessors ----------------------------------------
    //
    // Choose the sub-event index subsequent accessor calls refer to.  Clears
    // the product cache if the index changed; for PRad-II single-event data,
    // pass 0 (the default after Scan()).  Safe to call repeatedly.
    void SelectEvent(int i) const;

    // Each accessor decodes on first call after SelectEvent(), then returns
    // a cached reference.  References are invalidated by the next Read(),
    // Scan(), or SelectEvent() call.
    const fdec::EventInfo    &Info() const;  // always cheap
    const fdec::EventData    &Fadc() const;  // FADC250 + ADC1881M waveforms
    const ssp::SspEventData  &Gem()  const;  // SSP/MPD GEM strips
    const tdc::TdcEventData  &Tdc()  const;  // V1190 timing hits
    const vtp::VtpEventData  &Vtp()  const;  // VTP trigger clusters (see VtpData.h)

    // Absolute-time / run-state snapshot.  Unlike the other accessors this
    // one's result PERSISTS across events — it's only refreshed when Scan()
    // parses a SYNC/EPICS event (0xE112 HEAD bank) or a control event
    // (PRESTART/GO/END with a 3-word payload).  Intervening physics events
    // see the most-recent snapshot, giving them an absolute anchor for
    // their 48-bit TI timestamp delta.  To detect "a new SYNC arrived",
    // compare `Sync().sync_counter` to your last-seen value — it's a
    // monotonic counter for 0xE112 banks and stays 0 for control events
    // (distinguish those via `event_tag`).
    const psync::SyncInfo &Sync() const;

    // DSC2 scaler readout (0xE115 inside the TI master crate).  Populated
    // on physics events that carry the SYNC flag — at most once per scan.
    // Returns a record with `present == false` when no DSC2 bank is in this
    // event (or the configured slot/source doesn't match); see DscData.h.
    const dsc::DscEventData &Dsc() const;

    // EPICS slow-control record.  Populated on EPICS events (top-level
    // tag 0x001F) and stamped with absolute unix_time + sync_counter from
    // the 0xE112 HEAD bank plus the most-recent physics event_number seen
    // before this EPICS event.  Returns a record with `present == false`
    // on any other event type.
    const epics::EpicsRecord &Epics() const;

    // Most-recent physics event_number observed since this channel was
    // opened.  Refreshed automatically when a Physics event is scanned and
    // its TI/trigger-bank info has been decoded (Info() / DecodeEvent()).
    // Used to stamp slow events (EPICS, DSC2 SYNCs) so analysis can join
    // them to the physics tree by integer key.  Returns -1 before any
    // physics event has been seen.
    int32_t GetLastPhysicsEventNumber() const { return last_physics_event_number_; }

    // Companion to GetLastPhysicsEventNumber: TI 48-bit timestamp of the
    // same physics event.  Captured unconditionally at decode time so a
    // downstream filter that drops the event from the output tree does
    // not desynchronise slow rows from the timeline.  Returns 0 before
    // any physics event has been seen.
    uint64_t GetLastPhysicsTimestamp() const { return last_physics_timestamp_; }

    // --- legacy API (compat, writes directly to caller-owned structs) -------
    //
    // Preserves the original semantics: populate the caller's structs without
    // touching the lazy cache.  Kept so existing consumers compile unchanged;
    // new code should use SelectEvent() + Info()/Fadc()/Gem()/Tdc()/Vtp().
    bool DecodeEvent(int i, fdec::EventData &evt,
                     ssp::SspEventData *ssp_evt = nullptr,
                     vtp::VtpEventData *vtp_evt = nullptr,
                     tdc::TdcEventData *tdc_evt = nullptr) const;

    // --- EPICS extraction (call when GetEventType() == Epics) ---------------

    // Extract raw EPICS text from the current event buffer.
    // Returns the text payload (lines of "value  channel_name").
    // Returns empty string if no EPICS data found.
    std::string ExtractEpicsText() const;

    // --- DAQ config extraction (call when GetEventType() == Prestart) -------
    //
    // PRESTART events carry a 0xE10E STRING bank that holds the full
    // concatenated DAQ configuration text used by the run (TI / DSC2 /
    // FADC / TDC / SSP / VTP / TS settings, per-channel pedestals & gains,
    // trigger masks, etc.).  Returns the text payload (newline-separated
    // lines as written by CODA), or empty if the bank is absent.
    std::string ExtractDaqConfigText() const;

    // debug
    void PrintTree(std::ostream &os) const;

protected:
    DaqConfig config;
    int fHandle;
    std::vector<uint32_t> buffer;
    std::vector<EvNode>   nodes;
    int nevents = 0;
    EventType evtype = EventType::Unknown;
    int ra_count = 0;   // event count from evGetRandomAccessTable, 0 if sequential
    int ra_pos   = 0;   // next event index consumed by Read() in RA mode

    // tag → every node index in the current event that carries it.
    // Rebuilt by Scan(); consulted by the lazy accessors to avoid re-scanning.
    std::unordered_map<uint32_t, std::vector<int>> tag_index;

    // Per-product tag lists derived from config.data_banks at SetConfig().
    std::vector<uint32_t> fadc_tags;
    std::vector<uint32_t> gem_tags;
    std::vector<uint32_t> tdc_tags;
    std::vector<uint32_t> vtp_tags;

    // --- product cache (populated by Info/Fadc/Gem/Tdc/Vtp) -----------------
    // Cleared on Read/Scan/SelectEvent(i != current).  Heap-allocated on first
    // use — each product struct is sized for the worst-case event (EventData
    // alone is ~8 MB), so keeping them inline would blow the stack when
    // callers declare `EvChannel ch;` as a local.  Marked mutable so the
    // accessors stay const-callable — the cache is an implementation detail.
    mutable int  cached_event_idx = -1;
    mutable bool info_ready  = false;
    mutable bool fadc_ready  = false;
    mutable bool gem_ready   = false;
    mutable bool tdc_ready   = false;
    mutable bool vtp_ready   = false;
    mutable bool dsc_ready   = false;
    mutable bool epics_ready = false;
    mutable std::unique_ptr<fdec::EventData>   cache_fadc;   // .info also serves Info()
    mutable std::unique_ptr<ssp::SspEventData> cache_gem;
    mutable std::unique_ptr<tdc::TdcEventData> cache_tdc;
    mutable std::unique_ptr<vtp::VtpEventData> cache_vtp;
    mutable std::unique_ptr<dsc::DscEventData> cache_dsc;
    mutable std::unique_ptr<epics::EpicsRecord> cache_epics;

    // See GetLastPhysicsEventNumber() / GetLastPhysicsTimestamp().
    mutable int32_t  last_physics_event_number_ = -1;
    mutable uint64_t last_physics_timestamp_    = 0;

    // Persistent across Scan() — refreshed only when a SYNC/EPICS or control
    // event is scanned, otherwise carries the most recent snapshot.  The
    // `sync_decoded_this_event_` flag guards against re-decoding 0xE112 /
    // control payload more than once per event when Sync() is called
    // repeatedly, without clobbering the snapshot on physics events.
    mutable psync::SyncInfo last_sync_info_;
    mutable bool           sync_decoded_this_event_ = false;

    // --- per-bank decoders (shared by legacy and lazy paths) ----------------
    void decodeTriggerInfo(const EvNode &node, fdec::EventInfo &info) const;
    void decodeTIBank(const EvNode &node, fdec::EventInfo &info) const;
    void decodeRunInfo(const EvNode &node, fdec::EventInfo &info) const;

    // --- per-product dispatchers (write into caller-supplied structs) -------
    void decodeInfoInto(fdec::EventInfo &info) const;
    void decodeFadcInto(fdec::EventData &evt) const;    // fills evt.info too
    // Returns the total APV count across every SSP/MPD bank decoded — used by
    // the legacy DecodeEvent compat wrapper to preserve "true iff data found".
    int  decodeGemInto (ssp::SspEventData &ssp) const;
    void decodeTdcInto (tdc::TdcEventData &tdc) const;
    void decodeVtpInto (vtp::VtpEventData &vtp) const;
    // Reads 0xE112 HEAD (SYNC/EPICS events) or the first-child UINT32 bank of
    // a control event (PRESTART/GO/END) and merges into `out`, leaving prior
    // fields untouched when the current event contributes nothing.  Returns
    // true if `out` was updated.
    bool decodeSyncInto(psync::SyncInfo &out) const;
    // DSC2 0xE115 scaler bank, configured slot/source/channel from
    // config.dsc_scaler.  Leaves `out.present` false when there is no DSC2
    // bank or the slot doesn't match the configured one.
    void decodeDscInto(dsc::DscEventData &out) const;
    // EPICS slow-control bank.  Stamps with the latest Sync() unix_time /
    // sync_counter / run_number plus last_physics_event_number_.  Leaves
    // `out.present` false on non-EPICS events.
    void decodeEpicsInto(epics::EpicsRecord &out) const;

    // Invalidate all product cache flags.
    void clearCache() const;

    // Close any open handle, then evOpen `path` in evio `mode` ("r" or "ra")
    // with the random-access cursor reset.  Returns the evio status code.
    int openHandle(const std::string &path, const char *mode);

    size_t scanBank      (size_t off, int depth, int parent);
    // SEGMENT or TAGSEGMENT (1-word header): Hdr is SegmentHeader or
    // TagSegmentHeader.
    template <class Hdr>
    size_t scanSegment   (size_t off, int depth, int parent);
    void   scanChildren  (size_t off, size_t nwords, uint32_t ptype, int depth, int pidx);
};

} // namespace evc
