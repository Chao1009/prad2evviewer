//=============================================================================
// EvStruct                                                                  ||
// Basic information about CODA evio file format                             ||
//                                                                           ||
// Developer:                                                                ||
// Chao Peng                                                                 ||
// 09/07/2020                                                                ||
//                                                                           ||
// Updated: 2025 - Added TagSegmentHeader and CompositeHeader for            ||
//          composite bank support (tag 0xe126)                              ||
//=============================================================================
#pragma once

#include <cstdint>
#include <cstddef>


namespace evc
{

// evio bank data type
enum DataType {
    DATA_UNKNOWN32    =  (0x0),
    DATA_UINT32       =  (0x1),
    DATA_FLOAT32      =  (0x2),
    DATA_CHARSTAR8    =  (0x3),
    DATA_SHORT16      =  (0x4),
    DATA_USHORT16     =  (0x5),
    DATA_CHAR8        =  (0x6),
    DATA_UCHAR8       =  (0x7),
    DATA_DOUBLE64     =  (0x8),
    DATA_LONG64       =  (0x9),
    DATA_ULONG64      =  (0xa),
    DATA_INT32        =  (0xb),
    DATA_TAGSEGMENT   =  (0xc),
    DATA_ALSOSEGMENT  =  (0xd),
    DATA_ALSOBANK     =  (0xe),
    DATA_COMPOSITE    =  (0xf),
    DATA_BANK         =  (0x10),
    DATA_SEGMENT      =  (0x20)
};

// data word definitions
enum WordDefinition {
    BLOCK_HEADER = 0,
    BLOCK_TRAILER = 1,
    EVENT_HEADER = 2,
};

/* 32 bit bank header structure
 * ----------------------------------
 * |          length:32             |
 * ----------------------------------
 * |   tag:16   |:2| type:6 | num:8 |
 * ----------------------------------
 */
struct BankHeader
{
    uint32_t length, num, type, tag;

    BankHeader() : length(0) {}
    BankHeader(const uint32_t *buf)
    {
        length = buf[0];
        uint32_t word = buf[1];
        tag = (word >> 16) & 0xFFFF;
        type = (word >> 8) & 0x3F;
        num = (word & 0xFF);
    }

    static size_t size() { return 2; }
};

struct SegmentHeader
{
    uint32_t num, type, tag;

    SegmentHeader() {}
    SegmentHeader(const uint32_t *buf)
    {
        uint32_t word = buf[0];
        tag = (word >> 24) & 0xFF;
        type = (word >> 16) & 0x3F;
        num = (word & 0xFFFF);
    }

    static size_t size() { return 1; }
};

/* TagSegment header (used inside composite data):
 * ----------------------------------
 * | tag:12 | type:4 | length:16    |
 * ----------------------------------
 * Note: no padding field in tagsegments.
 */
struct TagSegmentHeader
{
    uint32_t length, type, tag;

    TagSegmentHeader() : length(0), type(0), tag(0) {}
    TagSegmentHeader(const uint32_t *buf)
    {
        uint32_t word = buf[0];
        tag    = (word >> 20) & 0xFFF;
        type   = (word >> 16) & 0xF;
        length = (word & 0xFFFF);
    }

    static size_t size() { return 1; }
};

/* Composite data envelope:
 *   [TagSegment header]   -- contains format string (type = DATA_CHARSTAR8)
 *   [format string words] -- padded to 32-bit boundary
 *   [Bank header]         -- contains the actual data payload
 *   [data payload words]
 *
 * CompositeHeader parses both the tagsegment and inner bank headers
 * from a pointer to the start of the composite data (i.e. the first
 * word after the outer bank header whose type == DATA_COMPOSITE).
 */
struct CompositeHeader
{
    TagSegmentHeader ts_header;    // tagsegment with format string
    BankHeader       data_header;  // inner bank with actual data

    // offsets in 32-bit words from the start of the composite data
    size_t format_offset;    // where the format string bytes start
    size_t data_offset;      // where the data payload starts (after inner bank header)
    size_t data_nwords;      // number of payload words

    CompositeHeader() : format_offset(0), data_offset(0), data_nwords(0) {}
    CompositeHeader(const uint32_t *buf)
    {
        // parse tagsegment header
        ts_header = TagSegmentHeader(buf);
        format_offset = TagSegmentHeader::size();

        // skip past the format string (ts_header.length words)
        size_t inner_bank_start = TagSegmentHeader::size() + ts_header.length;

        // parse the inner bank header
        data_header = BankHeader(buf + inner_bank_start);
        data_offset = inner_bank_start + BankHeader::size();
        data_nwords = data_header.length - 1; // length is exclusive of the second header word
    }

    // total size of the composite envelope + data in 32-bit words
    size_t total_words() const
    {
        return TagSegmentHeader::size() + ts_header.length + BankHeader::size() + data_nwords;
    }
};

struct BlockHeader
{
    bool valid;
    uint32_t nevents, number, module, slot;

    BlockHeader() : valid(false) {}
    BlockHeader(const uint32_t *buf)
    {
        uint32_t word = buf[0];
        valid = (word & 0x80000000) && (((word >> 27) & 0xF) == BLOCK_HEADER);
        slot = (word >> 22) & 0x1F;
        module = (word >> 18) & 0xF;
        number = (word >> 8) & 0x3FF;
        nevents = (word & 0xFF);
    }

    static size_t size() { return 1; }
};

struct BlockTrailer
{
    bool valid;
    uint32_t nwords, slot;

    BlockTrailer() : valid(false) {}
    BlockTrailer(const uint32_t *buf)
    {
        uint32_t word = buf[0];
        valid = (word & 0x80000000) && (((word >> 27) & 0xF) == BLOCK_TRAILER);
        slot = (word >> 22) & 0x1F;
        nwords = (word & 0x3FFFFF);
    }

    static size_t size() { return 1; }
};

struct EventHeader
{
    bool valid;
    uint32_t number, slot;

    EventHeader() : valid(false) {}
    EventHeader(const uint32_t *buf)
    {
        uint32_t word = buf[0];
        valid = (word & 0x80000000) && (((word >> 27) & 0xF) == EVENT_HEADER);
        slot = (word >> 22) & 0x1F;
        number = (word & 0x3FFFFF);
    }

    static size_t size() { return 1; }
};

} // namespace evc
