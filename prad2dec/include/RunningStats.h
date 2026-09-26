#pragma once
// RunningStats.h — running mean / population RMS accumulator.

#include <cmath>

namespace prad2 {

struct RunningStats {
    double sum   = 0.;
    double sum2  = 0.;
    int    count = 0;

    void   add(double v) { sum += v; sum2 += v * v; ++count; }
    double mean()  const { return count > 0 ? sum / count : 0.; }
    // 0 for fewer than two entries.
    double rms()   const {
        if (count < 2) return 0.;
        double m = mean();
        double var = sum2 / count - m * m;
        return var > 0 ? std::sqrt(var) : 0.;
    }
};

} // namespace prad2
