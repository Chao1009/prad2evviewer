#include "EvChannel.h"
#include "evio.h"
#include <cstring>
#include <iostream>
#include <iomanip>
#include <exception>
#include <algorithm>

using namespace evc;


// convert evio status to the enum
static inline status evio_status (int code)
{
    switch (code) {
    case S_SUCCESS:
        return status::success;
    case EOF:
    case S_EVFILE_UNXPTDEOF:
        return status::eof;
    case S_EVFILE_TRUNC:
        return status::incomplete;
    default:
        return status::failure;
    }
}

EvChannel::EvChannel(size_t buflen)
: fHandle(-1)
{
    buffer.resize(buflen);
}


status EvChannel::Open(const std::string &path)
{
    if (fHandle > 0) {
        Close();
    }
    char *cpath = strdup(path.c_str()), *copt = strdup("r");
    int status = evOpen(cpath, copt, &fHandle);
    free(cpath); free(copt);
    return evio_status(status);
}

void EvChannel::Close()
{
    evClose(fHandle);
    fHandle = -1;
}

status EvChannel::Read()
{
    return evio_status(evRead(fHandle, &buffer[0], buffer.size()));
}

bool EvChannel::ScanBanks(const std::vector<uint32_t> &banks)
{
    buffer_info.clear();
    composite_info.clear();

    auto evh = BankHeader(&buffer[0]);
    // skip the header
    size_t iword = BankHeader::size();

    // sanity checks
    if (evh.length > buffer.size()) {
        std::cout << "Ev Channel Error: Incomplete or corrupted event: event length = " << evh.length
                  << ", while buffer size is only " << buffer.size() << std::endl;
        return false;
    }

    if (evh.type != DATA_BANK) {
        std::cout << "Ev Channel Error: Expected DATA_BANK at the begining of an event, but got "
                  << evh.type << std::endl;
        return false;
    }

    // scan event, first one is the trigger bank
    try {
        iword += scanTriggerBank(&buffer[iword], iword);

        // scan ROC banks
        while (iword < evh.length + 1) {
            iword += scanRocBank(&buffer[iword], iword, banks);
        }

    } catch (std::exception const& e) {
        std::cerr << e.what() << std::endl;
        return false;
    }

    return true;
}

// scan trigger bank
size_t EvChannel::scanTriggerBank(const uint32_t *buf, size_t /* gindex */)
{
    auto header = BankHeader(buf);
    size_t iword = BankHeader::size();
    // sanity check
    if (header.type != DATA_SEGMENT) {
        throw(std::runtime_error("unexpected data type for trigger bank: " + std::to_string(header.type)));
    }

    // loop over the bank
    while (iword < header.length + 1) {
        SegmentHeader seg(buf + iword);
        // TODO, utilize the trigger bank info
        switch (seg.type) {
        // time stamp segment
        case DATA_ULONG64:
        // event type segment
        case DATA_USHORT16:
        // ROC segment
        case DATA_UINT32:
            break;
        // unexpected segment
        default:
            throw(std::runtime_error("unexpected segment in trigger bank: " + std::to_string(seg.type)));
        }
        // pass this segment
        iword += seg.num + 1;
    }

    return iword;
}

size_t EvChannel::scanRocBank(const uint32_t *buf, size_t gindex, const std::vector<uint32_t> &banks)
{
    auto header = BankHeader(buf);
    size_t iword = BankHeader::size();

    switch (header.type) {
    // bank of banks
    case DATA_BANK:
        break;
    case DATA_UINT32:
        throw(std::runtime_error("uint32_t data in ROC bank is not supported yet"));
    default:
        throw(std::runtime_error("unexpected data type in ROC bank: " + std::to_string(header.type)));
    }

    while (iword < header.length + 1) {
        auto bh = BankHeader(buf + iword);
        iword += BankHeader::size();

        // check if this bank is of interest
        bool interested = banks.empty() || (std::find(banks.begin(), banks.end(), bh.tag) != banks.end());

        if (interested) {
            if (bh.type == DATA_COMPOSITE) {
                // ---- Composite bank (e.g. tag 0xe126) ----
                scanCompositeBank(&buf[iword], bh.length - 1, header.tag, bh.tag, gindex + iword);
            } else {
                // ---- Legacy raw data bank ----
                scanDataBank(&buf[iword], bh.length - 1, header.tag, bh.tag, gindex + iword);
            }
        }

        iword += bh.length - 1;
    }

    return iword;
}

void EvChannel::scanDataBank(const uint32_t *buf, size_t buflen, uint32_t roc, uint32_t bank, size_t gindex)
{
    uint32_t slot, type, iev = 0;
    std::vector<BufferInfo> event_buffers;
    // scan the data bank
    for (size_t iword = 0; iword < buflen; ++iword) {
        // not a defininition word
        if (!(buf[iword] & 0x80000000)) { continue; }

        type = (buf[iword] >> 27) & 0xF;

        switch (type) {
        case BLOCK_HEADER:
            {
                BlockHeader blk(buf + iword);
                event_buffers.clear();
                slot = blk.slot;
            }
            break;
        case BLOCK_TRAILER:
            {
                BlockTrailer blk(buf + iword);
                if (slot != blk.slot) {
                    std::string mes = "warning: unmatched slot between block header (" + std::to_string(slot)
                                    + ") and trailer (" + std::to_string(blk.slot) + "), skip an event for roc "
                                    + std::to_string(roc) + " bank " + std::to_string(bank);
                    throw(std::runtime_error(mes));
                }
                if (event_buffers.size()) {
                    event_buffers.back().len = iword - event_buffers.back().len;
                    buffer_info[BufferAddress(roc, bank, slot)] = event_buffers;
                }
            }
            break;
        case EVENT_HEADER:
            {
                EventHeader evt(buf + iword);
                if (slot != evt.slot) {
                    std::string mes = "warning: unmatched slot between block header (" + std::to_string(slot)
                                    + ") and event header (" + std::to_string(evt.slot) + "), skip an event for roc "
                                    + std::to_string(roc) + " bank " + std::to_string(bank);
                    throw(std::runtime_error(mes));
                }
                if (event_buffers.size()) {
                    event_buffers.back().len = iword - event_buffers.back().len;
                }
                event_buffers.emplace_back(gindex + iword, iword);
            }
            break;
        // skip other headers
        default:
            break;
        }
    }
}

// ===================================================================
//  Composite bank scanner
//
//  Composite evio structure (e.g. tag 0xe126):
//
//    [TagSegment header: 1 word]   tag:12 | type:4 | length:16
//    [format string: ts.length words]     e.g. "c,m(c,ms)\0" padded
//    [Bank header: 2 words]               length | tag:16 | pad:2 | type:6 | num:8
//    [data payload: bank.length-1 words]  raw bytes in composite format
//
//  We parse the envelope to locate the data payload, then store a
//  CompositeInfo so the caller can feed the bytes to a decoder
//  (e.g. Fadc250Decoder::DecodeComposite).
// ===================================================================
void EvChannel::scanCompositeBank(const uint32_t *buf, size_t buflen, uint32_t roc, uint32_t bank, size_t gindex)
{
    if (buflen < 4) {
        std::cerr << "EvChannel Warning: composite bank too short (" << buflen << " words) in roc "
                  << roc << " bank 0x" << std::hex << bank << std::dec << "\n";
        return;
    }

    CompositeHeader ch(buf);

    // sanity: make sure the data payload fits
    if (ch.data_offset + ch.data_nwords > buflen) {
        std::cerr << "EvChannel Warning: composite data payload overflows bank boundary in roc "
                  << roc << " bank 0x" << std::hex << bank << std::dec
                  << " (data_offset=" << ch.data_offset << " data_nwords=" << ch.data_nwords
                  << " buflen=" << buflen << ")\n";
        return;
    }

    composite_info.emplace_back(roc, bank, static_cast<uint32_t>(gindex + ch.data_offset), static_cast<uint32_t>(ch.data_nwords));
}
