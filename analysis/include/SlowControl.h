#pragma once
//=============================================================================
// SlowControl.h — slow-event rows of replayed ROOT files (DSC2 `scalers`,
// `epics`), slice-local livetime and live-charge integration
//=============================================================================

#include "EventData.h"

#include <algorithm>
#include <cstdint>
#include <map>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

namespace analysis {

// One `scalers` row.  good is replay_filter's per-row verdict, true when the
// tree has no `good` branch.
struct ScalerRow : prad2::RawScalerData {
    bool good = true;
};

// One `epics` row.  ti_ticks is ti_ticks_at_arrival, 0 on files written
// before that branch existed; updates holds the readings carried by the row.
struct EpicsRow {
    int32_t  event_number = 0;   // event_number_at_arrival
    int64_t  ti_ticks     = 0;
    int64_t  unix_time    = 0;
    uint32_t sync_counter = 0;
    uint32_t run_number   = 0;
    bool     good         = true;
    std::map<std::string, double> updates;
};

// Append the rows of every file's `scalers` / `epics` tree, in file order;
// files without the tree are skipped.  A file that cannot be opened prints
// "<tag>: cannot open <path>" and returns false.  *has_good is set to true
// when any tree carries the `good` branch.
bool LoadScalerRows(const std::vector<std::string> &files, std::vector<ScalerRow> &out,
                    const char *tag, bool *has_good = nullptr);
bool LoadEpicsRows(const std::vector<std::string> &files, std::vector<EpicsRow> &out,
                   const char *tag, bool *has_good = nullptr);

// (gated, ungated) DSC2 counters of source "ref", "trg" or "tdc" (channel
// clamped to 0..15, unused for ref); {0, 0} for any other source.
std::pair<uint32_t, uint32_t> SelectDscPair(const prad2::RawScalerData &row,
                                            const std::string &source, int channel);

// Row indices in event_number order.
template <class Row>
std::vector<size_t> SortByEvent(const std::vector<Row> &rows)
{
    std::vector<size_t> idx(rows.size());
    std::iota(idx.begin(), idx.end(), size_t{0});
    std::sort(idx.begin(), idx.end(), [&](size_t a, size_t b) {
        return rows[a].event_number < rows[b].event_number;
    });
    return idx;
}

// Slice-local livetime of each row (indexed like rows, `order` from
// SortByEvent): Δgated / Δungated since the previous row, times scale
// (1 = fraction, 100 = percent).  -1 where it is undefined (ungated did not
// advance or gated advanced more).  The counters are cumulative, so the
// first row's delta is taken from (0, 0), i.e. it covers everything counted
// before it; a counter moving backward (DSC2 reset) also rebases to (0, 0).
std::vector<double> DeltaLivetime(const std::vector<ScalerRow> &rows,
                                  const std::vector<size_t> &order,
                                  const std::string &source, int channel,
                                  double scale = 1.0);

// Live charge Σ lf_b · Δt · ½(I_a + I_b) over adjacent checkpoints (a, b),
// with Δt from the TI ticks and lf_b the live fraction at b.  Beam currents
// in nA give nC.  Every pair with valid data enters the ungated sums, pairs
// added with good = true also the gated ones.  A pair with a missing or
// non-increasing tick, an unknown or negative lf_b, or an unknown current is
// counted as skipped.
struct ChargeSums {
    double  value_nC           = 0.0;
    double  live_seconds       = 0.0;
    double  real_seconds       = 0.0;
    int64_t n_pairs_integrated = 0;
    int64_t n_pairs_skipped    = 0;

    double  ungated_value_nC           = 0.0;
    double  ungated_live_seconds       = 0.0;
    double  ungated_real_seconds       = 0.0;
    int64_t n_ungated_pairs_integrated = 0;
    int64_t n_ungated_pairs_skipped    = 0;

    void AddPair(int64_t ticks_a, int64_t ticks_b, double live_fraction_b,
                 double current_a, double current_b, bool good);
    ChargeSums &operator+=(const ChargeSums &o);
};

} // namespace analysis
