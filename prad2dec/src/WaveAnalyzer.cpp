#include "WaveAnalyzer.h"

using namespace fdec;

// --- triangular-kernel smoothing (your SmoothSpectrum, zero-alloc) ----------
void WaveAnalyzer::smooth(const uint16_t *raw, int n, float *buf) const
{
    int res = cfg.resolution;
    if (res <= 1) {
        for (int i = 0; i < n; ++i) buf[i] = raw[i];
        return;
    }
    for (int i = 0; i < n; ++i) {
        float val = raw[i];
        float wsum = 1.0f;
        for (int j = 1; j < res; ++j) {
            if (j > i || i + j >= n) continue;
            float w = 1.0f - j / static_cast<float>(res + 1);
            val  += w * (raw[i - j] + raw[i + j]);
            wsum += 2.0f * w;
        }
        buf[i] = val / wsum;
    }
}

// --- iterative pedestal with outlier rejection (your CalcPedestal idea) -----
void WaveAnalyzer::findPedestal(const float *buf, int n, Pedestal &ped) const
{
    int nped = std::min(cfg.ped_nsamples, n);
    if (nped <= 0) { ped = {0, 0}; return; }

    // copy pedestal samples to scratch (on stack, max 200 samples)
    float scratch[MAX_SAMPLES];
    for (int i = 0; i < nped; ++i) scratch[i] = buf[i];

    // compute initial mean + rms
    auto calc = [](const float *s, int m, float &mean, float &rms) {
        float sum = 0, sum2 = 0;
        for (int i = 0; i < m; ++i) { sum += s[i]; sum2 += s[i] * s[i]; }
        mean = sum / m;
        float var = sum2 / m - mean * mean;
        rms = (var > 0) ? std::sqrt(var) : 0;
    };

    float mean, rms;
    calc(scratch, nped, mean, rms);

    // iterative outlier rejection: remove samples > 1σ from mean
    for (int iter = 0; iter < cfg.ped_max_iter; ++iter) {
        int count = 0;
        for (int i = 0; i < nped; ++i) {
            if (std::abs(scratch[i] - mean) < std::max(rms, cfg.ped_flatness))
                scratch[count++] = scratch[i];
        }
        if (count == nped || count < 5) break;
        nped = count;
        calc(scratch, nped, mean, rms);
    }

    ped.mean = mean;
    ped.rms  = rms;
}

// --- local-maxima peak search (your SearchMaxima approach, zero-alloc) ------
void WaveAnalyzer::findPeaks(const uint16_t *raw, const float *buf, int n,
                             float ped_mean, float ped_rms, float thr,
                             WaveResult &result) const
{
    result.npeaks = 0;
    if (n < 3) return;

    // trend: +1 rising, -1 falling, 0 flat
    auto trend = [](float a, float b) -> int {
        float d = a - b;
        return (std::abs(d) < 0.1f) ? 0 : (d > 0 ? 1 : -1);
    };

    for (int i = 1; i < n - 1 && result.npeaks < MAX_PEAKS; ++i) {
        int tr1 = trend(buf[i], buf[i - 1]);  // +1 if buf[i] > left
        int tr2 = trend(buf[i], buf[i + 1]);  // +1 if buf[i] > right

        // local maximum: higher than (or equal to) both neighbors, with at least one strict
        if (tr1 * tr2 < 0 || (tr1 == 0 && tr2 == 0)) continue;

        // handle flat plateau: if flat on the right side, walk to end of plateau
        // and use the center as the peak position
        int flat_end = i;
        if (tr2 == 0) {
            while (flat_end < n - 1 && trend(buf[flat_end], buf[flat_end + 1]) == 0)
                ++flat_end;
            // plateau must fall on the right to be a real maximum
            if (flat_end >= n - 1 || trend(buf[flat_end], buf[flat_end + 1]) <= 0)
                continue;
        }
        int peak_pos = (i + flat_end) / 2;

        // expand peak range: walk left while rising, walk right while falling/flat
        int left = i, right = flat_end;
        while (left > 0 && trend(buf[left], buf[left - 1]) > 0)
            --left;
        while (right < n - 1 && trend(buf[right], buf[right + 1]) >= 0)
            ++right;

        // estimate local baseline from edges (handles peaks on a slope)
        int span = right - left;
        if (span <= 0) continue;
        float base = (buf[left] * (right - peak_pos) + buf[right] * (peak_pos - left))
                   / static_cast<float>(span);

        // height above local baseline on smoothed data
        float smooth_height = buf[peak_pos] - base;
        if (smooth_height < thr) { i = right; continue; }

        // height above pedestal
        float ped_height = buf[peak_pos] - ped_mean;
        if (ped_height < thr || ped_height < 3.0f * ped_rms) { i = right; continue; }

        // --- integrate: walk outward from peak, stop at baseline crossing ---
        float integral = buf[peak_pos] - ped_mean;
        int int_left = peak_pos, int_right = peak_pos;

        for (int j = peak_pos - 1; j >= left; --j) {
            float v = buf[j] - ped_mean;
            if (v < ped_rms || v * ped_height < 0) { int_left = j; break; }
            integral += v;
            int_left = j;
        }
        for (int j = peak_pos + 1; j <= right; ++j) {
            float v = buf[j] - ped_mean;
            if (v < ped_rms || v * ped_height < 0) { int_right = j; break; }
            integral += v;
            int_right = j;
        }

        // --- correct peak position: find max in raw samples near smoothed peak ---
        int raw_pos = peak_pos;
        float raw_height = raw[peak_pos] - ped_mean;
        int search = std::max(1, cfg.resolution) + (flat_end - i) / 2;  // widen for plateaus
        for (int j = 1; j <= search; ++j) {
            if (peak_pos - j >= 0) {
                float h = raw[peak_pos - j] - ped_mean;
                if (h > raw_height) { raw_height = h; raw_pos = peak_pos - j; }
            }
            if (peak_pos + j < n) {
                float h = raw[peak_pos + j] - ped_mean;
                if (h > raw_height) { raw_height = h; raw_pos = peak_pos + j; }
            }
        }

        // --- fill peak ---
        Peak &p = result.peaks[result.npeaks];
        p.pos      = raw_pos;
        p.left     = int_left;
        p.right    = int_right;
        p.height   = raw_height;
        p.integral = integral;
        p.time     = raw_pos * 1e3f / cfg.clk_mhz;  // ns
        p.overflow = (raw[raw_pos] >= cfg.overflow);
        result.npeaks++;

        // skip past this peak's range to avoid double-counting
        i = right;
    }
}

// --- main entry point -------------------------------------------------------
void WaveAnalyzer::Analyze(const uint16_t *samples, int nsamples, WaveResult &result) const
{
    result.clear();
    if (!samples || nsamples <= 0 || nsamples > MAX_SAMPLES) return;

    // stack-allocated scratch buffer for smoothed waveform
    float buf[MAX_SAMPLES];
    smooth(samples, nsamples, buf);

    findPedestal(buf, nsamples, result.ped);

    float thr = std::max(cfg.threshold * result.ped.rms, cfg.min_threshold);
    findPeaks(samples, buf, nsamples, result.ped.mean, result.ped.rms, thr, result);
}
