//=============================================================================
// Class EvChannel                                                           ||
// Read event from CODA evio file, it can also scan data banks to locate     ||
// event buffers                                                             ||
//                                                                           ||
// Developer:                                                                ||
// Chao Peng                                                                 ||
// 09/07/2020                                                                ||
//                                                                           ||
// Updated: 2025 - Added composite bank scanning support                     ||
//=============================================================================
#pragma once

#include "EvStruct.h"
#include <iostream>
#include <string>
#include <vector>
#include <exception>
#include <unordered_map>


namespace evc {

// status enum
enum class status : int
{
    failure = -1,
    success = 1,
    incomplete = 2,
    empty = 3,
    eof = 4,
};

// buffer address
struct BufferAddress
{
    uint32_t roc, bank, slot;

    BufferAddress(uint32_t r = 0, uint32_t b = 0, uint32_t s = 0) : roc(r), bank(b), slot(s) {}
    bool operator== (const BufferAddress &a) const { return (roc == a.roc) && (bank == a.bank) && (slot == a.slot); }
};

// buffer info
struct BufferInfo
{
    uint32_t iword, len;

    BufferInfo(const uint32_t i = 0, uint32_t l = 0) : iword(i), len(l) {}
};

// Info for a composite bank: points to the data payload bytes within the event buffer
struct CompositeInfo
{
    uint32_t roc;           // ROC tag from the parent bank
    uint32_t tag;           // bank tag (e.g. 0xe126)
    uint32_t data_iword;    // word index into the event buffer where the data payload starts
    uint32_t data_nwords;   // number of 32-bit words of data payload

    CompositeInfo(uint32_t r = 0, uint32_t t = 0, uint32_t i = 0, uint32_t n = 0)
        : roc(r), tag(t), data_iword(i), data_nwords(n) {}
};

class BufferHash
{
public:
    size_t operator()(const BufferAddress& a) const { return ((a.roc << 24) | (a.bank << 8) | (a.slot)); }
};

class EvChannel
{
public:
    EvChannel(size_t buflen = 1024*2000);
    virtual ~EvChannel() { Close(); }

    EvChannel(const EvChannel &)  = delete;
    void operator =(const EvChannel &)  = delete;

    virtual status Open(const std::string &path);
    virtual void Close();
    virtual status Read();

    bool ScanBanks(const std::vector<uint32_t> &banks);
    bool Scan() { return ScanBanks({}); }

    uint32_t *GetRawBuffer() { return &buffer[0]; }
    const uint32_t *GetRawBuffer() const { return &buffer[0]; }

    std::vector<uint32_t> &GetRawBufferVec() { return buffer; }
    const std::vector<uint32_t> &GetRawBufferVec() const { return buffer; }

    BankHeader GetEvHeader() const { return BankHeader(&buffer[0]); }
    const std::unordered_map<BufferAddress, std::vector<BufferInfo>, BufferHash> &GetEvBuffers() const
    {
        return buffer_info;
    }
    const std::vector<BufferInfo> &GetEvBuffer(uint32_t roc, uint32_t bank, uint32_t slot) const
    {
        auto it = buffer_info.find(BufferAddress{roc, bank, slot});
        if (it != buffer_info.end()) {
            return it->second;
        }
        std::string error = "No data found for ROC " + std::to_string(roc)
                          + ", bank " + std::to_string(bank)
                          + ", slot " + std::to_string(slot);
        throw std::runtime_error(error);
    }

    const uint32_t *GetEvBuffer(uint32_t roc, uint32_t bank, uint32_t slot, uint32_t blk, size_t &len) const
    {
        auto it = buffer_info.find(BufferAddress{roc, bank, slot});
        if (it != buffer_info.end()) {
            auto &info = it->second[blk];
            len = info.len;
            return &buffer[info.iword];
        }
        std::string error = "No data found for ROC " + std::to_string(roc)
                          + ", bank " + std::to_string(bank)
                          + ", slot " + std::to_string(slot)
                          + ", block_level " + std::to_string(blk);
        throw std::runtime_error(error);
    }

    // ----- Composite bank access -----

    // Get all composite banks found in the last ScanBanks() call
    const std::vector<CompositeInfo> &GetCompositeInfos() const { return composite_info; }

    // Get the byte pointer and byte count for a composite bank's data payload
    // (for feeding into Fadc250Decoder::DecodeComposite)
    const uint8_t *GetCompositeData(const CompositeInfo &ci, size_t &nbytes) const
    {
        nbytes = ci.data_nwords * sizeof(uint32_t);
        return reinterpret_cast<const uint8_t*>(&buffer[ci.data_iword]);
    }

protected:
    size_t scanTriggerBank(const uint32_t *buf, size_t gindex);
    size_t scanRocBank(const uint32_t *buf, size_t gindex, const std::vector<uint32_t> &banks);
    void scanDataBank(const uint32_t *buf, size_t buflen, uint32_t roc, uint32_t bank, size_t gindex);
    void scanCompositeBank(const uint32_t *buf, size_t buflen, uint32_t roc, uint32_t bank, size_t gindex);

    int fHandle;
    std::vector<uint32_t> buffer;
    std::unordered_map<BufferAddress, std::vector<BufferInfo>, BufferHash> buffer_info;
    std::vector<CompositeInfo> composite_info;
};

} // namespace evc

// operator for address output
static std::ostream &operator << (std::ostream &os, const evc::BufferAddress &a)
{
    return os << "(" << a.roc << ", " << a.bank << ", " << a.slot << ")";
}
