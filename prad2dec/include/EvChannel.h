#pragma once
//=============================================================================
// EvChannel.h — read evio events, scan bank tree, decode FADC data
//
// Usage:
//   EvChannel ch;
//   ch.Open("file.evio");
//   while (ch.Read() == status::success) {
//       if (!ch.Scan()) continue;
//       auto hdr = ch.GetEvHeader();
//       if (hdr.tag != 0xfe) continue;       // skip non-physics
//
//       int nevt = ch.GetNEvents();           // events in this block
//       for (int i = 0; i < nevt; ++i) {
//           ch.DecodeEvent(i, event);         // fills pre-allocated EventData
//           // ... use event.rocs[r].slots[s].channels[c].samples[] ...
//       }
//   }
//=============================================================================

#include "EvStruct.h"
#include "Fadc250Data.h"
#include <string>
#include <vector>

namespace evc {

enum class status : int { failure = -1, success = 1, incomplete = 2, empty = 3, eof = 4 };

class EvChannel
{
public:
    EvChannel(size_t buflen = 1024 * 2000);
    virtual ~EvChannel() { Close(); }
    EvChannel(const EvChannel &) = delete;
    EvChannel &operator=(const EvChannel &) = delete;

    virtual status Open(const std::string &path);
    virtual void   Close();
    virtual status Read();

    // --- scan the current event into a flat tree ----------------------------
    bool Scan();

    // --- tree accessors -----------------------------------------------------
    BankHeader                  GetEvHeader() const { return BankHeader(&buffer[0]); }
    const std::vector<EvNode>  &GetNodes()    const { return nodes; }
    const EvNode               &GetChild(const EvNode &n, size_t i) const { return nodes[n.child_first + i]; }
    std::vector<const EvNode*>  FindByTag(uint32_t tag) const;

    const uint32_t *GetData(const EvNode &n) const { return &buffer[n.data_begin]; }
    const uint8_t  *GetBytes(const EvNode &n) const
    { return reinterpret_cast<const uint8_t*>(&buffer[n.data_begin]); }
    size_t GetDataBytes(const EvNode &n) const { return n.data_words * sizeof(uint32_t); }
    const uint8_t *GetCompositePayload(const EvNode &n, size_t &nbytes) const;

    uint32_t       *GetRawBuffer()       { return buffer.data(); }
    const uint32_t *GetRawBuffer() const { return buffer.data(); }

    // --- event-by-event access (after Scan) ---------------------------------

    // Number of events in this block (from the physics event header num field).
    // For single-event mode this is 1. For multi-event blocks this is M.
    int GetNEvents() const { return nevents; }

    // Decode the i-th event (0-based) into the pre-allocated EventData.
    // For single-event blocks (nevents=1), i must be 0.
    // Returns true on success.
    bool DecodeEvent(int i, fdec::EventData &evt) const;

    // debug
    void PrintTree(std::ostream &os) const;

protected:
    int fHandle;
    std::vector<uint32_t> buffer;
    std::vector<EvNode>   nodes;
    int nevents = 0;

    size_t scanBank      (size_t off, int depth, int parent);
    size_t scanSegment   (size_t off, int depth, int parent);
    size_t scanTagSegment(size_t off, int depth, int parent);
    void   scanChildren  (size_t off, size_t nwords, uint32_t ptype, int depth, int pidx);
};

} // namespace evc
