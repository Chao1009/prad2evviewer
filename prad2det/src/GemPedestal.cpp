#include "GemPedestal.h"
#include "SspData.h"
#include "DaqKey.h"
#include "RunningStats.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>

namespace gem {

namespace {

// Pack (crate, mpd, apv, strip) → uint64 key — fits in a flat std::map without
// needing a custom hash for a 4-tuple.  key >> 16 is the per-APV key.
inline uint64_t packKey(int crate, int mpd, int apv, int strip)
{
    return (prad2::pack_daq_key(crate, mpd, apv) << 16) |
           static_cast<uint16_t>(strip);
}

// Drop top/bottom N of 128 when computing the per-time-sample common mode.
constexpr int CM_DISCARD = 28;

}   // namespace

// --- Impl (pImpl) ------------------------------------------------------------
struct GemPedestal::Impl {
    std::map<uint64_t, prad2::RunningStats> accum;
};

GemPedestal::GemPedestal() : impl_(std::make_unique<Impl>()) {}
GemPedestal::~GemPedestal() = default;

void GemPedestal::Clear() { impl_->accum.clear(); }

int GemPedestal::NumStrips() const { return static_cast<int>(impl_->accum.size()); }

int GemPedestal::NumApvs() const
{
    std::map<uint64_t, int> apvs;
    for (const auto &[key, _] : impl_->accum) {
        uint64_t akey = (key >> 16);
        apvs[akey] += 1;
    }
    return static_cast<int>(apvs.size());
}

int GemPedestal::Accumulate(const ssp::SspEventData &evt)
{
    int    nfolded = 0;
    float  sorted[ssp::APV_STRIP_SIZE];
    float  cm_corrected[ssp::APV_STRIP_SIZE][ssp::SSP_TIME_SAMPLES];

    evt.forEachApv([&](const ssp::MpdData &mpd, int a, const ssp::ApvData &apv) {
        if (!apv.isFullReadout()) return;
        ++nfolded;

        // Per-time-sample common-mode subtraction.
        for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t) {
            for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s)
                sorted[s] = float(apv.strips[s][t]);

            std::sort(sorted, sorted + ssp::APV_STRIP_SIZE);

            double cm_sum  = 0.;
            int    cm_count = ssp::APV_STRIP_SIZE - 2 * CM_DISCARD;
            for (int s = CM_DISCARD; s < ssp::APV_STRIP_SIZE - CM_DISCARD; ++s)
                cm_sum += sorted[s];
            const double cm = cm_sum / cm_count;

            for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s)
                cm_corrected[s][t] = float(apv.strips[s][t]) - cm;
        }

        // Per-strip average over the 6 time samples, then accumulate.
        for (int s = 0; s < ssp::APV_STRIP_SIZE; ++s) {
            double avg = 0.;
            for (int t = 0; t < ssp::SSP_TIME_SAMPLES; ++t)
                avg += cm_corrected[s][t];
            avg /= ssp::SSP_TIME_SAMPLES;

            impl_->accum[packKey(mpd.crate_id, mpd.mpd_id, a, s)].add(avg);
        }
    });
    return nfolded;
}

int GemPedestal::Write(const std::string &output_path) const
{
    std::ofstream of(output_path);
    if (!of.is_open()) {
        std::cerr << "GemPedestal::Write: cannot write " << output_path << "\n";
        return -1;
    }

    // accum is ordered by (crate, mpd, apv, strip): each APV's strips are
    // contiguous and ascending.  The slot field is written as -1.
    of << std::fixed;
    int napvs = 0;
    uint64_t cur_apv = UINT64_MAX;
    for (const auto &[key, acc] : impl_->accum) {
        if ((key >> 16) != cur_apv) {
            cur_apv = key >> 16;
            ++napvs;
            of << "APV " << ((key >> 48) & 0xFFFF) << " -1 "
               << ((key >> 32) & 0xFFFF) << " " << ((key >> 16) & 0xFFFF) << "\n";
        }
        of << (key & 0xFFFF)
           << " " << std::setprecision(3) << std::round(acc.mean() *  1000.) /  1000.
           << " " << std::setprecision(4) << std::round(acc.rms()  * 10000.) / 10000. << "\n";
    }
    of.close();
    return of ? napvs : -1;
}

} // namespace gem
