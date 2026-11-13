//=============================================================================
// test_evio_scan.cpp — EvChannel::Scan() on hand-built evio events.
//
//   T1: well-formed nested event — node table, child counts, GetChild()
//       through siblings that have descendants, composite children
//   T2: child container whose length overruns the event → rejected
//   T3: leaf whose length overruns the event → rejected
//   T4: child overrunning its parent but inside the buffer (stale words after
//       the event, as in file mode) → rejected
//   T5: top-level length 0xFFFFFFFF (32-bit wrap) → rejected
//
// Build with -fsanitize=address to check that the malformed cases are
// rejected without reading past the buffer.
//
// Build target: prad2dec_test_evio_scan  (CTest: prad2dec.evio_scan).
//=============================================================================

#include "EvChannel.h"

#include <cstdio>
#include <vector>

using namespace evc;

namespace {
int g_total    = 0;
int g_failures = 0;

#define EXPECT_EQ(actual, expected) do {                                    \
    ++g_total;                                                              \
    auto _a = (actual); auto _e = (expected);                               \
    if (!(_a == _e)) {                                                      \
        ++g_failures;                                                       \
        std::fprintf(stderr,                                                \
            "  FAIL  %s:%d  %s == %s  (got %lld, want %lld)\n",             \
            __FILE__, __LINE__, #actual, #expected,                         \
            static_cast<long long>(_a), static_cast<long long>(_e));        \
    }                                                                       \
} while (0)

// Exposes the event buffer so a test can load hand-built words.
class TestChannel : public EvChannel
{
public:
    TestChannel() : EvChannel(0) {}
    void Load(const std::vector<uint32_t> &words) { buffer = words; }
};

uint32_t bank_word(uint32_t tag, uint32_t type, uint32_t num = 0)
{
    return (tag << 16) | (type << 8) | num;
}

// T1 layout (word index: content)
//   0  event   BANK  tag 0xFF50  len 21
//   2   ROC1   BANK  tag 0x0080  len 7
//   4    leaf  UINT32 tag 0xE10A len 2   (1 data word)
//   7    leaf  UINT32 tag 0xE101 len 2   (1 data word)
//  10   comp   COMPOSITE tag 0xE102 len 5
//  12    tagseg CHAR8 tag 1 len 1        (1 data word)
//  14    bank  UINT32 tag 0 len 1        (0 data words)
//  16   ROC2   BANK of SEGMENTs  tag 0x0082  len 5
//  18    seg-container SEGMENT tag 3 len 3
//  19     seg UINT32 tag 4 len 2         (2 data words)
//  22 (end)
std::vector<uint32_t> nested_event()
{
    std::vector<uint32_t> w = {
        21, bank_word(0xFF50, DATA_BANK, 1),
        7,  bank_word(0x0080, DATA_BANK),
        2,  bank_word(0xE10A, DATA_UINT32), 0x11,
        2,  bank_word(0xE101, DATA_UINT32), 0x22,
        5,  bank_word(0xE102, DATA_COMPOSITE),
        (1u << 20) | (DATA_CHAR8 << 16) | 1, 0x00004E28,
        1,  bank_word(0, DATA_UINT32),
        5,  bank_word(0x0082, DATA_SEGMENT),
        (3u << 24) | (DATA_SEGMENT << 16) | 3,
        (4u << 24) | (DATA_UINT32 << 16) | 2, 0x33, 0x44,
    };
    return w;
}

void test_nested()
{
    std::printf("T1: well-formed nested event\n");
    TestChannel ch;
    DaqConfig cfg;
    ch.SetConfig(cfg);
    ch.Load(nested_event());
    EXPECT_EQ(ch.Scan(), true);

    const auto &nodes = ch.GetNodes();
    EXPECT_EQ(nodes.size(), size_t(10));
    if (nodes.size() != 10) return;

    const EvNode &ev = nodes[0];
    EXPECT_EQ(ev.child_count, size_t(3));
    EXPECT_EQ(ch.GetChild(ev, 0).tag, 0x0080u);
    EXPECT_EQ(ch.GetChild(ev, 1).tag, 0xE102u);
    EXPECT_EQ(ch.GetChild(ev, 2).tag, 0x0082u);
    for (size_t i = 0; i < ev.child_count; ++i)
        EXPECT_EQ(ch.GetChild(ev, i).depth, 1);

    const EvNode &roc1 = ch.GetChild(ev, 0);
    EXPECT_EQ(roc1.child_count, size_t(2));
    EXPECT_EQ(ch.GetChild(roc1, 0).tag, 0xE10Au);
    EXPECT_EQ(ch.GetChild(roc1, 1).tag, 0xE101u);
    EXPECT_EQ(ch.GetChild(roc1, 1).data_words, size_t(1));
    EXPECT_EQ(ch.GetData(ch.GetChild(roc1, 1))[0], 0x22u);

    const EvNode &comp = ch.GetChild(ev, 1);
    EXPECT_EQ(comp.child_count, size_t(2));
    EXPECT_EQ(ch.GetChild(comp, 0).type, uint32_t(DATA_CHAR8));
    EXPECT_EQ(ch.GetChild(comp, 1).data_begin, size_t(16));

    const EvNode &roc2 = ch.GetChild(ev, 2);
    EXPECT_EQ(roc2.child_count, size_t(1));
    const EvNode &segc = ch.GetChild(roc2, 0);
    EXPECT_EQ(segc.child_count, size_t(1));
    EXPECT_EQ(ch.GetChild(segc, 0).tag, 4u);
    EXPECT_EQ(ch.GetChild(segc, 0).data_words, size_t(2));
    EXPECT_EQ(ch.GetData(ch.GetChild(segc, 0))[1], 0x44u);

    // last children terminate the sibling chain
    EXPECT_EQ(roc2.next_sibling, -1);
    EXPECT_EQ(ch.GetChild(roc1, 1).next_sibling, -1);
}

void expect_rejected(const char *name, const std::vector<uint32_t> &words)
{
    std::printf("%s\n", name);
    TestChannel ch;
    DaqConfig cfg;
    ch.SetConfig(cfg);
    ch.Load(words);
    EXPECT_EQ(ch.Scan(), false);
    EXPECT_EQ(ch.GetNodes().size(), size_t(0));
}

void test_container_overrun()
{
    // 16-word buffer, well-formed top level, child container claiming 1000
    std::vector<uint32_t> w(16, 0);
    w[0] = 15; w[1] = bank_word(0xFF50, DATA_BANK, 1);
    w[2] = 1000; w[3] = bank_word(0x0080, DATA_BANK);
    w[4] = 2; w[5] = bank_word(0xE10A, DATA_UINT32);
    expect_rejected("T2: child container overruns the event", w);
}

void test_leaf_overrun()
{
    std::vector<uint32_t> w(16, 0);
    w[0] = 15; w[1] = bank_word(0xFF50, DATA_BANK, 1);
    w[2] = 100000; w[3] = bank_word(0xE10A, DATA_UINT32);
    expect_rejected("T3: leaf overruns the event", w);
}

void test_parent_overrun_in_buffer()
{
    // event of 8 words in a 64-word buffer; ROC claims 4 words but holds a
    // leaf claiming 10 — inside the buffer, outside its parent
    std::vector<uint32_t> w(64, 0);
    w[0] = 7; w[1] = bank_word(0xFF50, DATA_BANK, 1);
    w[2] = 5; w[3] = bank_word(0x0080, DATA_BANK);
    w[4] = 10; w[5] = bank_word(0xE10A, DATA_UINT32);
    expect_rejected("T4: child overruns its parent inside the buffer", w);

    // same, with the ROC itself overrunning the event
    std::vector<uint32_t> v(64, 0);
    v[0] = 7; v[1] = bank_word(0xFF50, DATA_BANK, 1);
    v[2] = 20; v[3] = bank_word(0x0080, DATA_BANK);
    expect_rejected("T4b: ROC overruns the event inside the buffer", v);
}

void test_length_wrap()
{
    std::vector<uint32_t> w(16, 0);
    w[0] = 0xFFFFFFFFu; w[1] = bank_word(0xFF50, DATA_BANK, 1);
    expect_rejected("T5: top-level length 0xFFFFFFFF", w);

    std::vector<uint32_t> v(16, 0);
    v[0] = 15; v[1] = bank_word(0xFF50, DATA_BANK, 1);
    v[2] = 0xFFFFFFFFu; v[3] = bank_word(0xE10A, DATA_UINT32);
    expect_rejected("T5b: child length 0xFFFFFFFF", v);
}
} // namespace

int main()
{
    test_nested();
    test_container_overrun();
    test_leaf_overrun();
    test_parent_overrun_in_buffer();
    test_length_wrap();

    std::printf("\n%d checks, %d failures\n", g_total, g_failures);
    return g_failures == 0 ? 0 : 1;
}
