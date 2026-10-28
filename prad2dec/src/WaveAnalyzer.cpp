#include "WaveAnalyzer.h"
#include "PulseTemplateStore.h"

#include <limits>

using namespace fdec;

namespace {
inline float cholesky(const float *M, float *L, int K);
inline void chol_solve(const float *L, const float *b, float *x, int K);

// Levenberg-Marquardt damping schedule shared by every fit in this file.
constexpr float LM_LAMBDA0    = 1.0e-3f;
constexpr float LM_LAMBDA_UP  = 10.0f;
constexpr float LM_LAMBDA_DN  = 10.0f;
constexpr float LM_LAMBDA_MAX = 1.0e10f;

// One damped Gauss-Newton step: solve (JᵀJ + λI) δ = Jᵀr for npar ≤ NPMAX
// parameters, with J[j * MAX_SAMPLES + i] = ∂r_i/∂p_j over n samples.
// Weighted multiplies every JᵀJ and Jᵀr term by w.  Returns false when the
// damped matrix is not positive definite (the caller raises λ and retries).
template <int NPMAX, bool Weighted = false>
inline bool lm_step(const float *J, const float *r, int n, int npar,
                    float lambda, float *delta, float w = 1.0f)
{
    float A[NPMAX * NPMAX], L[NPMAX * NPMAX], g[NPMAX];
    for (int j1 = 0; j1 < npar; ++j1) {
        g[j1] = 0.0f;
        for (int j2 = 0; j2 <= j1; ++j2) A[j1 * npar + j2] = 0.0f;
    }
    for (int i = 0; i < n; ++i) {
        for (int j1 = 0; j1 < npar; ++j1) {
            const float jj1 = J[j1 * MAX_SAMPLES + i];
            if constexpr (Weighted) g[j1] += jj1 * r[i] * w;
            else                    g[j1] += jj1 * r[i];
            for (int j2 = 0; j2 <= j1; ++j2) {
                const float jj2 = J[j2 * MAX_SAMPLES + i];
                if constexpr (Weighted) A[j1 * npar + j2] += jj1 * jj2 * w;
                else                    A[j1 * npar + j2] += jj1 * jj2;
            }
        }
    }
    // cholesky() reads only the lower triangle.
    for (int j = 0; j < npar; ++j) A[j * npar + j] += lambda;
    if (cholesky(A, L, npar) < 0.0f) return false;
    chol_solve(L, g, delta, npar);
    return true;
}

// Last j in [max(lo, 0), hi) where the pedsub signal rises through thr,
// v[j] < thr <= v[j+1] with v = buf - ped, by more than min_dv; -1 if none.
inline int last_rising_crossing(const float *buf, int lo, int hi, float ped,
                                float thr, float min_dv)
{
    int last = -1;
    for (int j = std::max(lo, 0); j < hi; ++j) {
        const float v0 = buf[j]     - ped;
        const float v1 = buf[j + 1] - ped;
        if (v0 < thr && v1 >= thr && v1 - v0 > min_dv) last = j;
    }
    return last;
}
}

float WaveAnalyzer::log_normal_pulse_value(float t, float ped_mean, float A,
                                           float t0, float mu, float sigma)
{
    if (!(sigma > 0.0f) || t <= t0) return ped_mean;
    const float dt = t - t0;
    const float z = (std::log(dt) - mu) / sigma;
    return ped_mean + A * std::exp(-0.5f * z * z);
}

float WaveAnalyzer::log_normal_cfd_time_sample(float t0, float mu, float sigma,
                                               float cfd_fraction)
{
    if (!(cfd_fraction > 0.0f) || !(cfd_fraction < 1.0f))
        return std::numeric_limits<float>::quiet_NaN();
    const float root_term = std::sqrt(-2.0f * std::log(cfd_fraction));
    return t0 + std::exp(mu - sigma * root_term);
}

LogNormalFitResult WaveAnalyzer::fit_log_normal_cfd(const float *buf, int nsamples,
                                                    int fit_left_bound, int raw_pos,
                                                    float cfd_fraction,
                                                    float ped_mean, float ped_rms)
{
    LogNormalFitResult res;
    if (!buf || nsamples <= 0 || raw_pos < 0 || raw_pos >= nsamples) return res;
    if (!(cfd_fraction > 0.0f) || !(cfd_fraction < 1.0f)) return res;

    const int fit_start = std::max(0, fit_left_bound - 4);
    const int fit_stop  = std::min(nsamples - 1, raw_pos + 4);
    const int nfit = fit_stop - fit_start + 1;
    if (nfit < 8) return res;

    // Only fit when the rising edge actually crosses the CFD level.
    const float v_thr = cfd_fraction * (buf[raw_pos] - ped_mean);
    if (last_rising_crossing(buf, fit_left_bound - 2, raw_pos, ped_mean, v_thr, 0.0f) < 0)
        return res;

    constexpr int NPAR = 4;
    float t[MAX_SAMPLES];
    float y[MAX_SAMPLES];
    for (int i = 0; i < nfit; ++i) {
        const int idx = fit_start + i;
        t[i] = static_cast<float>(idx);
        y[i] = buf[idx];
    }

    const float peak_time = static_cast<float>(raw_pos);
    const float A_guess = std::max(buf[raw_pos] - ped_mean, 1.0f);

    float p[NPAR] = {
        A_guess,
        std::max(0.0f, fit_left_bound - 0.5f),
        1.0f,
        0.5f,
    };
    float p_lo[NPAR] = {
        0.0f,
        std::max(0.0f, fit_left_bound - 1.5f),
        0.5f,
        0.1f,
    };
    float p_hi[NPAR] = {
        3.0f * A_guess,
        peak_time - 1.0e-3f,
        2.0f,
        0.9f,
    };
    for (int j = 0; j < NPAR; ++j)
        p[j] = std::clamp(p[j], p_lo[j], p_hi[j]);

    const float sigma_fit = std::max(ped_rms, 1.0f);
    const float inv_sigma2 = 1.0f / (sigma_fit * sigma_fit);
    const int dof = std::max(1, nfit - NPAR);

    auto eval_chi2 = [&](const float *params) -> float {
        if (!(params[0] > 0.0f) || !(params[3] > 0.0f) || !(params[1] < peak_time))
            return std::numeric_limits<float>::infinity();
        float chi2 = 0.0f;
        for (int i = 0; i < nfit; ++i) {
            const float r = y[i] - log_normal_pulse_value(
                t[i], ped_mean, params[0], params[1], params[2], params[3]);
            chi2 += r * r;
        }
        return (chi2 * inv_sigma2) / static_cast<float>(dof);
    };

    float chi2 = eval_chi2(p);
    if (!std::isfinite(chi2)) return res;

    float p_best[NPAR];
    for (int j = 0; j < NPAR; ++j) p_best[j] = p[j];
    float chi2_best = chi2;

    constexpr int   MAX_ITER = 60;
    constexpr float STEP_TOL = 1.0e-4f;

    float lambda = LM_LAMBDA0;

    float r[MAX_SAMPLES];
    float J[NPAR * MAX_SAMPLES];
    float delta[NPAR], p_new[NPAR];

    for (int iter = 0; iter < MAX_ITER; ++iter) {
        for (int i = 0; i < nfit; ++i) {
            r[i] = y[i] - log_normal_pulse_value(t[i], ped_mean, p[0], p[1], p[2], p[3]);
        }

        const float h_A     = std::max(1.0e-3f * std::max(std::abs(p[0]), 1.0f), 1.0e-3f);
        const float h_t0    = 1.0e-3f;
        const float h_mu    = std::max(1.0e-3f * std::max(std::abs(p[2]), 1.0f), 1.0e-4f);
        const float h_sigma = std::max(1.0e-3f * std::max(std::abs(p[3]), 1.0f), 1.0e-4f);
        const float h[NPAR] = {h_A, h_t0, h_mu, h_sigma};

        for (int j = 0; j < NPAR; ++j) {
            float p_fd[NPAR];
            for (int k = 0; k < NPAR; ++k) p_fd[k] = p[k];
            p_fd[j] = std::min(p_hi[j], p[j] + h[j]);
            const float step = p_fd[j] - p[j];
            if (!(step > 0.0f)) {
                for (int i = 0; i < nfit; ++i) J[j * MAX_SAMPLES + i] = 0.0f;
                continue;
            }
            for (int i = 0; i < nfit; ++i) {
                const float f0 = log_normal_pulse_value(t[i], ped_mean, p[0],    p[1],    p[2],    p[3]);
                const float fp = log_normal_pulse_value(t[i], ped_mean, p_fd[0], p_fd[1], p_fd[2], p_fd[3]);
                J[j * MAX_SAMPLES + i] = -(fp - f0) / step;
            }
        }

        if (!lm_step<NPAR, true>(J, r, nfit, NPAR, lambda, delta, inv_sigma2)) {
            lambda *= LM_LAMBDA_UP;
            if (lambda > LM_LAMBDA_MAX) break;
            continue;
        }

        float dpar = 0.0f;
        for (int j = 0; j < NPAR; ++j) {
            p_new[j] = std::clamp(p[j] - delta[j], p_lo[j], p_hi[j]);
            const float scale = std::max(std::abs(p[j]), 1.0e-3f);
            dpar += std::abs(p_new[j] - p[j]) / scale;
        }

        const float chi2_new = eval_chi2(p_new);
        if (std::isfinite(chi2_new) && chi2_new < chi2_best) {
            chi2_best = chi2_new;
            for (int j = 0; j < NPAR; ++j) p_best[j] = p_new[j];
        }

        if (std::isfinite(chi2_new) && chi2_new < chi2) {
            for (int j = 0; j < NPAR; ++j) p[j] = p_new[j];
            chi2 = chi2_new;
            lambda = std::max(lambda / LM_LAMBDA_DN, 1.0e-12f);
            if (dpar < STEP_TOL) break;
        } else {
            lambda *= LM_LAMBDA_UP;
            if (lambda > LM_LAMBDA_MAX) break;
        }
    }

    const float t_cfd_sample = log_normal_cfd_time_sample(
        p_best[1], p_best[2], p_best[3], cfd_fraction);
    const float t_peak_sample = p_best[1] + std::exp(p_best[2]);
    const float fit_lo_sample = static_cast<float>(fit_start);
    const float fit_hi_sample = static_cast<float>(raw_pos);
    if (!std::isfinite(t_cfd_sample) || !std::isfinite(t_peak_sample)) return res;
    if (!(t_cfd_sample >= fit_lo_sample && t_cfd_sample <= fit_hi_sample)) return res;
    if (!(t_peak_sample >= fit_lo_sample && t_peak_sample <= static_cast<float>(fit_stop) + 1.0f)) return res;

    res.ok = true;
    res.t_cfd_sample = t_cfd_sample;
    res.chi2_per_dof = chi2_best;
    res.A = p_best[0];
    res.t0 = p_best[1];
    res.mu = p_best[2];
    res.sigma = p_best[3];
    return res;
}

// --- triangular-kernel smoothing (zero-alloc) --------------------------------
void WaveAnalyzer::smooth(const uint16_t *raw, int n, float *buf) const
{
    int res = cfg.smooth_order;
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

// --- iterative pedestal with median/MAD bootstrap + outlier rejection -------
//
// Median + MAD (×1.4826) seed is robust against ≤50% contamination — a
// previous-event tail or early ringing in the leading window biases the
// simple-mean seed badly, which then loosens the σ-clip band and the
// iteration can converge on a contaminated baseline.  Median-seeded σ-clip
// recovers the right baseline immediately and matches the simple-mean
// behaviour on clean baselines.
void WaveAnalyzer::findPedestal(const float *buf, int start, int nped,
                                Pedestal &ped) const
{
    ped = {};
    if (nped <= 0) return;

    // Copy the window plus original sample indices (needed for slope below,
    // since the survivor set after σ-clip is a subset of the window).
    float scratch[MAX_SAMPLES];
    int   orig_idx[MAX_SAMPLES];
    for (int i = 0; i < nped; ++i) {
        scratch[i]  = buf[start + i];
        orig_idx[i] = start + i;
    }
    int active = nped;

    // ── Median + MAD bootstrap.
    float sorted[MAX_SAMPLES];
    auto median = [&]() {
        std::sort(sorted, sorted + nped);
        return (nped % 2 == 1) ? sorted[nped / 2]
                               : 0.5f * (sorted[nped / 2 - 1] + sorted[nped / 2]);
    };
    for (int i = 0; i < nped; ++i) sorted[i] = scratch[i];
    float mean = median();
    for (int i = 0; i < nped; ++i) sorted[i] = std::abs(scratch[i] - mean);
    const float mad = median();
    float rms = mad * 1.4826f;        // MAD → σ for normally-distributed noise

    // ── Iterative σ-clip from the robust seed.  scratch / orig_idx track
    // surviving samples in lock-step so we can compute slope on the actual
    // survivor set (not on samples that pass the final band post-hoc).
    bool converged = false;
    for (int iter = 0; iter < cfg.ped_max_iter; ++iter) {
        const float band = std::max(rms, cfg.ped_flatness);
        int count = 0;
        for (int i = 0; i < active; ++i) {
            if (std::abs(scratch[i] - mean) < band) {
                scratch[count]  = scratch[i];
                orig_idx[count] = orig_idx[i];
                ++count;
            }
        }
        if (count == active) { converged = true; break; }
        if (count < 5) {
            ped.quality |= Q_PED_TOO_FEW_SAMPLES;
            active = count;
            break;     // keep prior mean/rms — too few survivors to refit
        }
        active = count;
        float sum = 0, sum2 = 0;
        for (int i = 0; i < active; ++i) { sum += scratch[i]; sum2 += scratch[i] * scratch[i]; }
        mean = sum / active;
        const float var = sum2 / active - mean * mean;
        rms = (var > 0) ? std::sqrt(var) : 0;
    }
    if (!converged && !(ped.quality & Q_PED_TOO_FEW_SAMPLES))
        ped.quality |= Q_PED_NOT_CONVERGED;
    if (rms < cfg.ped_flatness)
        ped.quality |= Q_PED_FLOOR_ACTIVE;

    // ── Linear least-squares slope on the survivors (ADC/sample).  Catches
    // baseline drift / pulse-tail contamination that the σ-clip alone can
    // hide (e.g. a slow tail tilts every sample similarly so none of them
    // register as outliers).
    float slope = 0.0f;
    if (active >= 2) {
        double sx = 0, sy = 0;
        for (int i = 0; i < active; ++i) { sx += orig_idx[i]; sy += scratch[i]; }
        const double xbar = sx / active, ybar = sy / active;
        double sxy = 0, sxx = 0;
        for (int i = 0; i < active; ++i) {
            const double dx = orig_idx[i] - xbar;
            sxy += dx * (scratch[i] - ybar);
            sxx += dx * dx;
        }
        if (sxx > 0) slope = static_cast<float>(sxy / sxx);
    }

    ped.mean  = mean;
    ped.rms   = rms;
    ped.nused = static_cast<uint8_t>(active < 255 ? active : 255);
    ped.slope = slope;
}

// --- local-maxima peak search (zero-alloc) ----------------------------------
void WaveAnalyzer::findPeaks(const uint16_t *raw, const float *buf, int n,
                             float ped_mean, float ped_rms, float thr,
                             WaveResult &result) const
{
    result.npeaks = 0;
    if (n < 3) return;

    // track peak-finding ranges (left/right) separately from integration bounds
    int pk_range[MAX_PEAKS][2];  // [i][0]=left, [i][1]=right

    // Trend: +1 rising, -1 falling, 0 flat.  The flat-tolerance scales with
    // the pedestal RMS — a hardcoded 0.1 ADC threshold treats noise-level
    // wiggles as "rising/falling" on quiet channels (ped_rms ~ 0.5), which
    // splits real plateaus into spurious mini-peaks.  Floor at 0.1 keeps
    // behaviour reasonable on raw integer ADC.
    const float trend_tol = std::max(0.1f, 0.5f * ped_rms);
    auto trend = [trend_tol](float a, float b) -> int {
        const float d = a - b;
        return (std::abs(d) < trend_tol) ? 0 : (d > 0 ? 1 : -1);
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
        if (ped_height < thr) { i = right; continue; }

        // --- integrate: walk outward from peak, stop at baseline or tail cutoff ---
        // Termination requires N = cfg.tail_break_n consecutive sub-threshold
        // samples, so a single noise dip in the tail does not truncate the
        // integral early.  Below-threshold samples seen during a not-yet-
        // confirmed run are held in `pending` and either committed (on
        // recovery) or discarded (when the run reaches N).
        //
        // int_left / int_right are INCLUSIVE bounds: they're only advanced
        // when a sample is actually added to `integral`, so they always
        // point to the outermost above-threshold sample on each side.
        float integral = buf[peak_pos] - ped_mean;
        const float tail_cut = ped_height * cfg.int_tail_ratio;
        const int   N_break  = std::max(1, cfg.tail_break_n);
        int int_left = peak_pos, int_right = peak_pos;

        auto is_below = [&](float v) {
            return v < tail_cut || v < ped_rms || v * ped_height < 0;
        };

        // Walk from peak_pos in direction step (±1) within [left, right].
        auto walk = [&](int step, int &bound) {
            int   below_run = 0;
            float pending   = 0.0f;
            for (int j = peak_pos + step; j >= left && j <= right; j += step) {
                const float v = buf[j] - ped_mean;
                if (is_below(v)) {
                    ++below_run;
                    pending += v;
                    if (below_run >= N_break) break;
                } else {
                    integral += pending + v;
                    pending = 0.0f;
                    below_run = 0;
                    bound = j;
                }
            }
        };
        walk(-1, int_left);
        walk(+1, int_right);

        // --- correct peak position: find max in raw samples near smoothed peak ---
        int raw_pos = peak_pos;
        float raw_height = raw[peak_pos] - ped_mean;
        int search = std::max(1, cfg.smooth_order) + (flat_end - i) / 2;  // widen for plateaus
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

        // --- reject if overlapping a previous peak and local height too small ---
        // Use peak-finding range (left/right), not integration bounds, for overlap test.
        // smooth_height is the height above the line connecting left/right edges,
        // i.e., how much this peak rises above the tail it sits on.
        bool rejected = false;
        for (int k = 0; k < result.npeaks; ++k) {
            if (left <= pk_range[k][1] && right >= pk_range[k][0]) {
                if (smooth_height < result.peaks[k].height * cfg.min_peak_ratio) {
                    rejected = true;
                    break;
                }
            }
        }
        if (rejected) { i = right; continue; }

        // --- quadratic peak-time interpolation ---
        // Fit y = a x² + b x + c through the 3 raw samples around raw_pos;
        // the parabola vertex sits at δ = (h[-1] - h[+1]) / (2·(h[-1] - 2·h[0] + h[+1]))
        // relative to raw_pos.  Lifts the time resolution from 4 ns
        // (sample-quantised) to ≪ 1 ns for clean peaks.  Guarded by:
        //   - raw_pos not at the buffer edge,
        //   - denom < 0 (real concave-down max — flat plateaus and numerical
        //     noise have denom ≥ 0 and skip interpolation),
        //   - δ clamped to ±1 sample for robustness.
        float t_subsample = 0.0f;
        if (raw_pos > 0 && raw_pos < n - 1) {
            const float h_minus = raw[raw_pos - 1];
            const float h_zero  = raw[raw_pos];
            const float h_plus  = raw[raw_pos + 1];
            const float denom = h_minus - 2.0f * h_zero + h_plus;
            if (denom < -1e-3f) {
                const float delta = 0.5f * (h_minus - h_plus) / denom;
                t_subsample = std::max(-1.0f, std::min(1.0f, delta));
            }
        }

        // --- linear interpolation constant fraction discrimination ---
        // For large pulses, use a simple digital CFD time on the raw
        // pedsub samples; if crossing search or interpolation fails,
        // fall back to quadratic peak-time interpolation above.
        float t_pickoff = raw_pos + t_subsample - 1.5f;  // 50% cfd is ~1.5 samples before peak
        const float cfd_fraction = 0.5f;
        uint8_t time_algo = T_PICKOFF_PEAKING_SUBSAMPLE;
        if (raw_height > 5.0f * ped_rms) {
            const float v_thr = cfd_fraction * raw_height;

            // Scan the rising edge from int_left - 1 up to the peak position.
            // Keep the last valid crossing so the pickoff stays nearest
            // to the peak if small pre-rise oscillations exist.
            const int jc = last_rising_crossing(buf, int_left - 1, raw_pos, ped_mean,
                                                v_thr, 1.5f * ped_rms);
            if (jc >= 0) {
                const float v0 = buf[jc]     - ped_mean;
                const float v1 = buf[jc + 1] - ped_mean;
                const float frac = (v_thr - v0) / (v1 - v0);
                t_pickoff = jc + std::max(0.0f, std::min(1.0f, frac));
                time_algo = T_PICKOFF_LINEAR_CFD;
            } else {
                t_pickoff = raw_pos + t_subsample;
            }
        }

        // Refine the leading-edge time with the local Log-Normal CFD
        // idea. If the bounded fit becomes
        // non-physical or unstable, keep the simpler CFD / quadratic fallback.
        LogNormalFitResult fit_info;
        if (raw_height > 10.0f * ped_rms) {
            fit_info = fit_log_normal_cfd(
                buf, n, int_left, raw_pos, cfd_fraction, ped_mean, ped_rms);
            if (fit_info.ok) {
                t_pickoff = fit_info.t_cfd_sample;
                time_algo = T_PICKOFF_FIT_CFD;
            }
        }

        // --- pile-up detection ---
        // Flag this peak (and the matching previously-found peak) when
        // their integration windows touch or overlap within
        // cfg.peak_pileup_gap samples — diagnostic for downstream cuts on
        // isolated vs piled-up pulses.
        uint8_t my_quality = Q_PEAK_GOOD;
        const int gap = std::max(1, cfg.peak_pileup_gap);
        for (int k = 0; k < result.npeaks; ++k) {
            const Peak &prev = result.peaks[k];
            if (int_left  <= prev.right + gap &&
                int_right >= prev.left  - gap) {
                result.peaks[k].quality |= Q_PEAK_PILED;
                my_quality |= Q_PEAK_PILED;
            }
        }

        // --- fill peak ---
        Peak &p = result.peaks[result.npeaks];
        p.pos      = raw_pos;
        p.left     = int_left;
        p.right    = int_right;
        p.height   = raw_height;
        p.integral = integral;
        p.time     = t_pickoff * 1e3f / cfg.clk_mhz;  // ns
        p.overflow = (raw[raw_pos] >= cfg.overflow);
        p.quality  = my_quality;
        p.time_algo = time_algo;
        result.peaks_fit[result.npeaks] = fit_info;
        pk_range[result.npeaks][0] = left;
        pk_range[result.npeaks][1] = right;
        result.npeaks++;

        // skip past this peak's range to avoid double-counting
        i = right;
    }
}

// --- main entry point -------------------------------------------------------
void WaveAnalyzer::Analyze(const uint16_t *samples, int nsamples, WaveResult &result, float time_offset) const
{
    result.clear();
    if (!samples || nsamples <= 0 || nsamples > MAX_SAMPLES) return;

    // stack-allocated scratch buffer for smoothed waveform
    float buf[MAX_SAMPLES];
    smooth(samples, nsamples, buf);

    auto window_overflow = [&](int wstart, int wlen) -> bool {
        const uint16_t ovr = cfg.overflow;
        for (int i = wstart; i < wstart + wlen; ++i)
            if (samples[i] >= ovr) return true;
        return false;
    };

    const int nped_window = std::min(cfg.ped_nsamples, nsamples);

    // ── Leading-window pedestal estimate.
    Pedestal P_lead;
    findPedestal(buf, 0, nped_window, P_lead);
    if (window_overflow(0, nped_window))
        P_lead.quality |= Q_PED_OVERFLOW;

    // ── Adaptive: if the leading window looks suspicious (didn't converge,
    // lost > 50% of samples to rejection, or hit overflow), try the
    // trailing window — only if the two don't overlap.  Pick whichever
    // has the lower RMS (with nused as tiebreaker); flag the choice with
    // Q_PED_TRAILING_WINDOW.
    Pedestal P_use         = P_lead;
    int      ped_win_start = 0;
    const bool lead_suspicious =
        (P_lead.quality & (Q_PED_NOT_CONVERGED |
                           Q_PED_TOO_FEW_SAMPLES |
                           Q_PED_OVERFLOW))
        || (P_lead.nused * 2 < nped_window);

    if (lead_suspicious && nsamples >= 2 * nped_window) {
        const int trail_start = nsamples - nped_window;
        Pedestal P_trail;
        findPedestal(buf, trail_start, nped_window, P_trail);
        if (window_overflow(trail_start, nped_window))
            P_trail.quality |= Q_PED_OVERFLOW;
        const bool trail_better =
            (P_trail.rms < P_lead.rms) ||
            (P_trail.rms == P_lead.rms && P_trail.nused > P_lead.nused);
        if (trail_better) {
            P_use         = P_trail;
            P_use.quality |= Q_PED_TRAILING_WINDOW;
            ped_win_start = trail_start;
        }
    }
    result.ped = P_use;

    // ── Peak finding uses the chosen pedestal.
    const float thr = std::max(cfg.peak_nsigma * result.ped.rms, cfg.min_peak_height);
    findPeaks(samples, buf, nsamples, result.ped.mean, result.ped.rms, thr, result);

    // ── Post-hoc: was a real pulse inside the pedestal window we used?
    // Diagnostic for downstream filters — doesn't influence the estimate
    // (the median+MAD seed already absorbs single-pulse contamination on
    // most channels), but lets analyses optionally cut on clean events.
    const int ped_win_end = ped_win_start + nped_window;
    for (int p = 0; p < result.npeaks; ++p) {
        const int pos = result.peaks[p].pos;
        if (pos >= ped_win_start && pos < ped_win_end) {
            result.ped.quality |= Q_PED_PULSE_IN_WINDOW;
            break;
        }
    }

    // ── Pile-up deconvolution.  Silent no-op unless the caller bound a
    // template store and the current channel key.
    applyAutoDeconv(samples, nsamples, result);

    // ── Subtract the caller's time offset (ns) from every peak time.
    for (int p = 0; p < result.npeaks; ++p)
        result.peaks[p].time -= time_offset;
}

void WaveAnalyzer::applyAutoDeconv(const uint16_t *samples, int nsamples,
                                   WaveResult &result) const
{
    if (!cfg.nnls_deconv.enabled)               return;
    if (template_store_ == nullptr)             return;
    if (ck_roc_ < 0 || ck_slot_ < 0 || ck_chan_ < 0) return;
    if (result.npeaks <= 0)                     return;

    // Cheap gate: skip clean events unless config asks otherwise.  This
    // is the dominant cost saving — most channels see no pile-up and
    // we'd otherwise pay a deconvolution fit on every event.
    if (!cfg.nnls_deconv.apply_to_all_peaks) {
        bool any_piled = false;
        for (int k = 0; k < result.npeaks; ++k) {
            if (result.peaks[k].quality & Q_PEAK_PILED) {
                any_piled = true; break;
            }
        }
        if (!any_piled) return;
    }

    const PulseTemplate *tmpl = template_store_->Lookup(
        ck_roc_, ck_slot_, ck_chan_);
    if (tmpl == nullptr) return;

    DeconvOutput out;
    Deconvolve(samples, nsamples, result, *tmpl, out);

    // On success: replace each peak's height/integral with the deconv
    // values and mark Q_PEAK_DECONVOLVED.  Failure paths
    // (Q_DECONV_BAD_TEMPLATE / Q_DECONV_LM_NOT_CONVERGED) leave the peaks as-is
    // so downstream code falls back to WaveAnalyzer's tail-cutoff values.
    if (out.state == Q_DECONV_APPLIED || out.state == Q_DECONV_FALLBACK_GLOBAL) {
        const int K = (out.n < result.npeaks) ? out.n : result.npeaks;
        for (int k = 0; k < K; ++k) {
            result.peaks[k].height    = out.height[k];
            result.peaks[k].integral  = out.integral[k];
            result.peaks[k].quality  |= Q_PEAK_DECONVOLVED;
        }
    }
}

// === Per-pulse-fit pile-up deconvolution ====================================
//
// Given pedsub samples b[0..n-1] and K peak times τ_1..τ_K from the
// WaveAnalyzer, we fit a 4K-parameter model
//
//   model(t_i) = Σ_k a_k · T(t_i; t0_k, τ_r_k, τ_f_k) / T_max(τ_r_k, τ_f_k)
//
// via Levenberg-Marquardt with the channel template providing the
// initial guess (a_k = WA peak.height / T_max(template), t0_k from
// peak.time, (τ_r, τ_f) = template) and tight bounds around it
// (cfg.nnls_deconv.shape_window_factor, t0_window_ns, amp_max_factor).
//
// Algorithm: classical LM with forward-difference Jacobian.  The
// shape-param FD is per-peak (each FD perturbation only recomputes
// the k-th template column, not the whole model — analytic block
// structure of the Jacobian).

namespace {

// Lower-triangular Cholesky factorisation of a KxK SPD matrix M (row-major,
// L is also row-major).  Returns the smallest pivot squared (= L[k,k]²) so
// callers can do a conditioning check.  Returns -1 on failure (M not SPD).
inline float cholesky(const float *M, float *L, int K)
{
    float min_pivot_sq = std::numeric_limits<float>::infinity();
    for (int i = 0; i < K; ++i) {
        for (int j = 0; j <= i; ++j) {
            float sum = M[i * K + j];
            for (int k = 0; k < j; ++k)
                sum -= L[i * K + k] * L[j * K + k];
            if (i == j) {
                if (sum <= 0.0f) return -1.0f;
                if (sum < min_pivot_sq) min_pivot_sq = sum;
                L[i * K + j] = std::sqrt(sum);
            } else {
                L[i * K + j] = sum / L[j * K + j];
            }
        }
        // Zero the strict upper triangle so chol_solve is well-defined
        for (int j = i + 1; j < K; ++j) L[i * K + j] = 0.0f;
    }
    return min_pivot_sq;
}

// Solve L L^T x = b using a precomputed Cholesky factor.  Sized to
// accommodate the per-pulse-fit deconvolver's 4·MAX_PEAKS-parameter
// Hessian (so K can be up to 32, not just MAX_PEAKS=8).
inline void chol_solve(const float *L, const float *b, float *x, int K)
{
    float y[4 * MAX_PEAKS];
    // Forward: L y = b
    for (int i = 0; i < K; ++i) {
        float s = b[i];
        for (int k = 0; k < i; ++k) s -= L[i * K + k] * y[k];
        y[i] = s / L[i * K + i];
    }
    // Back: L^T x = y
    for (int i = K - 1; i >= 0; --i) {
        float s = y[i];
        for (int k = i + 1; k < K; ++k) s -= L[k * K + i] * x[k];
        x[i] = s / L[i * K + i];
    }
}

// Two-tau pulse shape (1 − e^(−dt/τ_r)) · e^(−dt/τ_f), dt ns after onset.
inline float two_tau_shape(float dt, float tr, float tf)
{
    return (1.0f - std::exp(-dt / tr)) * std::exp(-dt / tf);
}

// Peak value of two_tau_shape, (1 − u) · u^(τ_r/τ_f) with u = τ_r/(τ_r + τ_f).
inline float two_tau_max(float tr, float tf)
{
    const float u = tr / (tr + tf);
    return (1.0f - u) * std::pow(u, tr / tf);
}

// Closed-form template peak position offset (ns) and peak value.
//   t_peak = τ_r · ln((τ_r + τ_f) / τ_r)
inline void template_peak(float tr, float tf, float &t_off, float &t_max)
{
    const float u = tr / (tr + tf);
    t_off = tr * std::log(1.0f / u);
    t_max = two_tau_max(tr, tf);
}

// Unit-amplitude template at sample times t_i = i·clk_ns, with template
// onset at t0 (ns).  Writes n values into `col`.  Two exp() per sample;
// (τ_r, τ_f) vary per peak per LM iteration, so columns are always
// computed analytically.
inline void template_column(float *col, int n, float clk_ns,
                            float t0, float tr, float tf)
{
    for (int i = 0; i < n; ++i) {
        const float t = i * clk_ns;
        col[i] = (t <= t0) ? 0.0f : two_tau_shape(t - t0, tr, tf);
    }
}

} // namespace (anonymous)

void WaveAnalyzer::Deconvolve(const uint16_t *samples, int nsamples,
                              const WaveResult &wres,
                              const PulseTemplate &tmpl,
                              DeconvOutput &dec_out) const
{
    dec_out.clear();

    // Explicit API: always runs when given valid inputs and a usable
    // template.  The cfg.nnls_deconv.enabled gate only governs the auto
    // path (applyAutoDeconv inside Analyze) so the Python diagnostic
    // can compute deconv values without flipping the production switch.
    if (!samples || nsamples <= 0 || nsamples > MAX_SAMPLES) {
        dec_out.state = Q_DECONV_NOT_RUN;
        return;
    }
    const int K = wres.npeaks;
    if (K <= 0 || K > MAX_PEAKS) {
        dec_out.state = Q_DECONV_NOT_RUN;
        return;
    }

    const auto &dcfg = cfg.nnls_deconv;
    const float tr_tmpl = tmpl.tau_r_ns;
    const float tf_tmpl = tmpl.tau_f_ns;
    if (!(tr_tmpl >= dcfg.tau_r_min_ns && tr_tmpl <= dcfg.tau_r_max_ns) ||
        !(tf_tmpl >= dcfg.tau_f_min_ns && tf_tmpl <= dcfg.tau_f_max_ns) ||
        !(tr_tmpl > 0.0f) || !(tf_tmpl > 0.0f)) {
        dec_out.state = Q_DECONV_BAD_TEMPLATE;
        return;
    }

    const float clk_ns = cfg.clk_ns();
    const float ped    = wres.ped.mean;
    const float sigma  = std::max(wres.ped.rms, 1.0f);
    const float inv_sigma2 = 1.0f / (sigma * sigma);

    // Pedsub data.
    float b[MAX_SAMPLES];
    for (int i = 0; i < nsamples; ++i)
        b[i] = static_cast<float>(samples[i]) - ped;

    // Initial param vector laid out as 4K floats with the per-peak block
    //   [a_k, t0_k, τ_r_k, τ_f_k]
    // at index 4*k.  Bounds in p_lo / p_hi follow the same layout.
    constexpr int NPP = 4;                          // params per peak
    constexpr int NPMAX = NPP * MAX_PEAKS;          // 32
    const int npar = NPP * K;

    float t_off_tmpl, t_max_tmpl;
    template_peak(tr_tmpl, tf_tmpl, t_off_tmpl, t_max_tmpl);
    const float inv_t_max_tmpl = 1.0f / t_max_tmpl;

    const float fac = std::max(dcfg.shape_window_factor, 1.001f);
    const float t0_win = std::max(dcfg.t0_window_ns, 0.0f);

    float p[NPMAX]   = {0};
    float p_lo[NPMAX] = {0};
    float p_hi[NPMAX] = {0};
    float peak_t[MAX_PEAKS];

    for (int k = 0; k < K; ++k) {
        const float t_pk = wres.peaks[k].time;           // ns from sample 0
        peak_t[k] = t_pk;
        const float h_wa = std::max(wres.peaks[k].height, 1.0f);
        const float a_init  = h_wa * inv_t_max_tmpl;     // height / T_max
        const float t0_init = t_pk - t_off_tmpl;

        p[NPP*k + 0] = a_init;
        p[NPP*k + 1] = t0_init;
        p[NPP*k + 2] = tr_tmpl;
        p[NPP*k + 3] = tf_tmpl;

        p_lo[NPP*k + 0] = 0.0f;
        p_hi[NPP*k + 0] = std::max(dcfg.amp_max_factor, 1.001f) * h_wa
                         * inv_t_max_tmpl;
        p_lo[NPP*k + 1] = t0_init - t0_win;
        p_hi[NPP*k + 1] = t0_init + t0_win;
        p_lo[NPP*k + 2] = tr_tmpl / fac;
        p_hi[NPP*k + 2] = tr_tmpl * fac;
        p_lo[NPP*k + 3] = tf_tmpl / fac;
        p_hi[NPP*k + 3] = tf_tmpl * fac;
    }

    // Helpers to (re)compute the K template columns at given params and
    // evaluate the residual chi².
    auto compute_M = [&](const float *params, float *M) {
        for (int k = 0; k < K; ++k) {
            template_column(&M[k * MAX_SAMPLES], nsamples, clk_ns,
                            params[NPP*k + 1],
                            params[NPP*k + 2],
                            params[NPP*k + 3]);
        }
    };
    auto eval_chi2 = [&](const float *params, const float *M) -> float {
        float c = 0.0f;
        for (int i = 0; i < nsamples; ++i) {
            float fit = 0.0f;
            for (int k = 0; k < K; ++k)
                fit += params[NPP*k + 0] * M[k * MAX_SAMPLES + i];
            const float rr = (b[i] - fit);
            c += rr * rr;
        }
        return c * inv_sigma2;
    };

    // Initial M and chi².
    float M[MAX_PEAKS * MAX_SAMPLES];
    compute_M(p, M);
    float chi2 = eval_chi2(p, M);

    // Best-seen tracking (return whatever we found, scipy-style).
    float chi2_best = chi2;
    float p_best[NPMAX];
    for (int j = 0; j < npar; ++j) p_best[j] = p[j];

    // LM hyperparameters.
    constexpr int   MAX_ITER  = 100;
    constexpr float TOL_PARAM = 1e-5f;

    float lambda = LM_LAMBDA0;
    int iter = 0;
    bool any_accepted = false;

    // Scratch buffers reused across iterations.
    float r_vec[MAX_SAMPLES];
    float J[NPMAX * MAX_SAMPLES];   // J[j * MAX_SAMPLES + i] = ∂r_i/∂p_j
    float delta[NPMAX];
    float p_new[NPMAX];
    float M_new[MAX_PEAKS * MAX_SAMPLES];
    float col_perturbed[MAX_SAMPLES];

    for (; iter < MAX_ITER; ++iter) {
        // Residuals at current point.
        for (int i = 0; i < nsamples; ++i) {
            float fit = 0.0f;
            for (int k = 0; k < K; ++k)
                fit += p[NPP*k + 0] * M[k * MAX_SAMPLES + i];
            r_vec[i] = b[i] - fit;
        }

        // Jacobian — block structure: only the k-th template column
        // depends on (t0_k, τ_r_k, τ_f_k); all peaks contribute to a_k
        // through the analytic ∂model/∂a_k = -T_k.
        for (int k = 0; k < K; ++k) {
            const float a_k  = p[NPP*k + 0];
            const float t0_k = p[NPP*k + 1];
            const float tr_k = p[NPP*k + 2];
            const float tf_k = p[NPP*k + 3];

            // ∂r/∂a_k = -T_k(t_i)  (analytic — reuse M).
            for (int i = 0; i < nsamples; ++i)
                J[(NPP*k + 0) * MAX_SAMPLES + i] = -M[k * MAX_SAMPLES + i];

            const float h_t0 = std::max(1e-3f * clk_ns, 1e-6f);
            const float h_tr = std::max(1e-3f * tr_k,    1e-6f);
            const float h_tf = std::max(1e-3f * tf_k,    1e-6f);

            // ∂r/∂t0_k = -a_k · ∂T_k/∂t0_k
            template_column(col_perturbed, nsamples, clk_ns,
                            t0_k + h_t0, tr_k, tf_k);
            for (int i = 0; i < nsamples; ++i) {
                const float dT = (col_perturbed[i] - M[k * MAX_SAMPLES + i]) / h_t0;
                J[(NPP*k + 1) * MAX_SAMPLES + i] = -a_k * dT;
            }
            template_column(col_perturbed, nsamples, clk_ns,
                            t0_k, tr_k + h_tr, tf_k);
            for (int i = 0; i < nsamples; ++i) {
                const float dT = (col_perturbed[i] - M[k * MAX_SAMPLES + i]) / h_tr;
                J[(NPP*k + 2) * MAX_SAMPLES + i] = -a_k * dT;
            }
            template_column(col_perturbed, nsamples, clk_ns,
                            t0_k, tr_k, tf_k + h_tf);
            for (int i = 0; i < nsamples; ++i) {
                const float dT = (col_perturbed[i] - M[k * MAX_SAMPLES + i]) / h_tf;
                J[(NPP*k + 3) * MAX_SAMPLES + i] = -a_k * dT;
            }
        }

        if (!lm_step<NPMAX>(J, r_vec, nsamples, npar, lambda, delta)) {
            // Indefinite — back off and try a tighter LM step.
            lambda *= LM_LAMBDA_UP;
            if (lambda > LM_LAMBDA_MAX) break;
            continue;
        }

        // Trial point — same sign convention as FitPulseShape:
        //   p_new = p − δ, then clamp to bounds.
        float dpar = 0.0f;
        for (int j = 0; j < npar; ++j) {
            p_new[j] = std::clamp(p[j] - delta[j], p_lo[j], p_hi[j]);
            const float scale = std::max(std::abs(p[j]), 1e-3f);
            dpar += std::abs(p_new[j] - p[j]) / scale;
        }

        compute_M(p_new, M_new);
        const float chi2_new = eval_chi2(p_new, M_new);

        if (std::isfinite(chi2_new) && chi2_new < chi2_best) {
            chi2_best = chi2_new;
            for (int j = 0; j < npar; ++j) p_best[j] = p_new[j];
        }

        if (chi2_new < chi2) {
            // Accept.
            for (int j = 0; j < npar; ++j) p[j] = p_new[j];
            for (int idx = 0; idx < K * MAX_SAMPLES; ++idx) M[idx] = M_new[idx];
            chi2 = chi2_new;
            lambda = std::max(lambda / LM_LAMBDA_DN, 1e-12f);
            any_accepted = true;
            if (dpar < TOL_PARAM) { ++iter; break; }
        } else {
            lambda *= LM_LAMBDA_UP;
            if (lambda > LM_LAMBDA_MAX) break;
        }
    }

    if (!any_accepted) {
        dec_out.state = Q_DECONV_LM_NOT_CONVERGED;
        return;
    }

    // Use the best-seen params.  Recompute M one more time so the
    // integral window sums use the converged shape.
    for (int j = 0; j < npar; ++j) p[j] = p_best[j];
    compute_M(p, M);

    dec_out.n = K;
    for (int k = 0; k < K; ++k) {
        const float a_k  = p[NPP*k + 0];
        const float t0_k = p[NPP*k + 1];
        const float tr_k = p[NPP*k + 2];
        const float tf_k = p[NPP*k + 3];

        if (!std::isfinite(a_k) || !std::isfinite(t0_k) ||
            !std::isfinite(tr_k) || !std::isfinite(tf_k)) {
            dec_out.state = Q_DECONV_LM_NOT_CONVERGED;
            return;
        }

        float t_off_k, t_max_k;
        template_peak(tr_k, tf_k, t_off_k, t_max_k);

        dec_out.amplitude[k] = a_k;
        dec_out.height[k]    = a_k * t_max_k;
        dec_out.t0_ns[k]     = t0_k;
        dec_out.tau_r_ns[k]  = tr_k;
        dec_out.tau_f_ns[k]  = tf_k;

        // Per-peak integral: window centred on the WaveAnalyzer peak
        // index (consistent with how downstream consumers slice raw data),
        // summed against the per-peak fitted template.
        const int i_pk = static_cast<int>(std::lround(peak_t[k] / clk_ns));
        const int lo   = std::max(0,        i_pk - dcfg.pre_samples);
        const int hi   = std::min(nsamples, i_pk + dcfg.post_samples + 1);
        float sum = 0.0f;
        for (int i = lo; i < hi; ++i) sum += M[k * MAX_SAMPLES + i];
        dec_out.integral[k] = a_k * sum;
    }

    const int dof = std::max(1, nsamples - npar);
    dec_out.chi2_per_dof = chi2_best / static_cast<float>(dof);
    dec_out.state = tmpl.is_global ? Q_DECONV_FALLBACK_GLOBAL
                                   : Q_DECONV_APPLIED;
}

// === Per-pulse shape fit (Levenberg-Marquardt on normalised two-tau model) ===
//
// Jacobian via forward finite differences (FD) on the model itself —
// analytic partials are correct but tedious and bring no real speed-up
// here since the dominant cost is the exp() calls inside the model
// evaluator and FD reuses those.

namespace {

// Unit-amplitude pulse models for fit_two_tau.  q = (t0, τ_r, τ_f[, p]) with
// times in ns; unit() is the model at time t times t_max_inv = 1 / t_max(q).
struct TwoTau {
    static constexpr int NP = 3;
    static float t_max(const float *q) { return two_tau_max(q[1], q[2]); }
    static float unit(float t, const float *q, float t_max_inv)
    {
        if (t <= q[0]) return 0.0f;
        return two_tau_shape(t - q[0], q[1], q[2]) * t_max_inv;
    }
};

// Two-tau with rise-edge exponent p (p = 1 is TwoTau up to rounding):
//   T(t)  = [1 − exp(−(t−t0)/τ_r)]^p · exp(−(t−t0)/τ_f)
//   T_max = u_peak^p · (τ_r/(τ_r + p·τ_f))^(τ_r/τ_f),  u_peak = p·τ_f/(τ_r + p·τ_f)
struct TwoTauP {
    static constexpr int NP = 4;
    static float t_max(const float *q)
    {
        const float tr = q[1], tf = q[2], p = q[3];
        if (!(tr > 0.0f) || !(tf > 0.0f) || !(p > 0.0f)) return 1.0f;
        const float denom = tr + p * tf;
        return std::pow(p * tf / denom, p) * std::pow(tr / denom, tr / tf);
    }
    static float unit(float t, const float *q, float t_max_inv)
    {
        if (t <= q[0]) return 0.0f;
        const float dt = t - q[0];
        const float u  = 1.0f - std::exp(-dt / q[1]);
        if (u <= 0.0f) return 0.0f;
        return std::pow(u, q[3]) * std::exp(-dt / q[2]) * t_max_inv;
    }
};

// Sum of squared residuals of model M against the normalised pulse y.
template <class M>
float two_tau_chi2(const float *y, int n, float clk_ns, const float *q)
{
    const float t_max_inv = 1.0f / M::t_max(q);
    float s = 0.0f;
    for (int i = 0; i < n; ++i) {
        const float r = y[i] - M::unit(i * clk_ns, q, t_max_inv);
        s += r * r;
    }
    return s;
}

// Levenberg-Marquardt fit of model M to one slice after pedsub and
// normalisation to the raw peak height (see WaveAnalyzer::FitPulseShape).
// peak_amp is set as soon as the slice validates; q[0..M::NP), chi2_per_dof
// and n_iter only when the fit succeeds.
template <class M>
bool fit_two_tau(const uint16_t *slice, int nslice, int peak_idx,
                 float ped, float ped_rms, float clk_ns, float model_err_floor,
                 float *q, float &peak_amp, float &chi2_per_dof, int &n_iter)
{
    constexpr int   NP        = M::NP;
    constexpr int   MAX_ITER  = 100;       // ~scipy curve_fit default budget
    constexpr float TOL_PARAM = 1e-5f;     // relative param step in clk_ns units

    if (!slice || nslice < 8 || nslice > MAX_SAMPLES) return false;
    if (peak_idx < 0 || peak_idx >= nslice) return false;

    const float amp = static_cast<float>(slice[peak_idx]) - ped;
    if (amp <= 0.0f) return false;
    peak_amp = amp;

    float y[MAX_SAMPLES];
    const float inv_amp = 1.0f / amp;
    for (int i = 0; i < nslice; ++i)
        y[i] = (static_cast<float>(slice[i]) - ped) * inv_amp;

    // σ on the normalised pulse: relative noise, floored at the model-
    // error scale so χ²/dof stays sane on high-amplitude pulses.
    const float sigma_noise = std::max(ped_rms, 1.0f) * inv_amp;
    const float sigma       = std::max(sigma_noise, model_err_floor);
    const float w_inv2      = 1.0f / (sigma * sigma);

    // Initial guesses and bounds for (t0, τ_r, τ_f, p).  p starts at the
    // sigmoidal 2.0; p < 1 (sharp rise) is allowed but unusual.
    const float init[4] = {peak_idx * clk_ns - 2.0f * clk_ns,
                           1.0f * clk_ns, 5.0f * clk_ns, 2.0f};
    const float lo[4]   = {-2.0f * clk_ns, 0.2f * clk_ns, 1.0f * clk_ns, 0.3f};
    const float hi[4]   = {(nslice - 1) * clk_ns, 10.0f * clk_ns,
                           80.0f * clk_ns, 10.0f};

    float p[NP], p_new[NP], p_best[NP];
    std::copy(init, init + NP, p);
    std::copy(init, init + NP, p_best);
    float lambda = LM_LAMBDA0;
    float chi2   = two_tau_chi2<M>(y, nslice, clk_ns, p);

    // The best point seen is returned unconditionally (scipy curve_fit
    // semantics); requiring a strictly improving step would report ~12% of
    // fits whose initial guess is already near-optimal as failures.
    float chi2_best = chi2;

    int iter = 0;
    for (; iter < MAX_ITER; ++iter) {
        float r[MAX_SAMPLES];
        const float t_max_inv = 1.0f / M::t_max(p);
        for (int i = 0; i < nslice; ++i)
            r[i] = y[i] - M::unit(i * clk_ns, p, t_max_inv);

        // Forward-difference Jacobian; T_max is re-evaluated at each
        // perturbed point (it does not depend on t0).
        float J[NP * MAX_SAMPLES];
        for (int j = 0; j < NP; ++j) {
            float ph[NP];
            std::copy(p, p + NP, ph);
            const float h = std::max(1e-3f * (j == 0 ? clk_ns : p[j]), 1e-6f);
            ph[j] = p[j] + h;
            const float tmi_h = 1.0f / M::t_max(ph);
            for (int i = 0; i < nslice; ++i) {
                const float t  = i * clk_ns;
                const float f0 = M::unit(t, p,  t_max_inv);
                const float fp = M::unit(t, ph, tmi_h);
                J[j * MAX_SAMPLES + i] = -(fp - f0) / h;
            }
        }

        // The uniform weight w_inv2 factors out of the normal equations.
        float delta[NP];
        if (!lm_step<NP>(J, r, nslice, NP, lambda, delta)) {
            lambda *= LM_LAMBDA_UP;
            if (lambda > LM_LAMBDA_MAX) break;
            continue;
        }

        // J = −∂model/∂p, so the Gauss-Newton step is p − δ.
        for (int j = 0; j < NP; ++j)
            p_new[j] = std::clamp(p[j] - delta[j], lo[j], hi[j]);
        const float chi2_new = two_tau_chi2<M>(y, nslice, clk_ns, p_new);

        if (std::isfinite(chi2_new) && chi2_new < chi2_best) {
            chi2_best = chi2_new;
            std::copy(p_new, p_new + NP, p_best);
        }

        if (chi2_new < chi2) {
            float dpar = 0.0f;   // times in clk_ns units, p relative to itself
            for (int j = 0; j < NP; ++j)
                dpar += std::abs(delta[j]) / (j < 3 ? clk_ns : std::max(p[j], 0.1f));
            std::copy(p_new, p_new + NP, p);
            chi2   = chi2_new;
            lambda = std::max(lambda / LM_LAMBDA_DN, 1.0e-12f);
            if (dpar < TOL_PARAM) { ++iter; break; }
        } else {
            lambda *= LM_LAMBDA_UP;
            if (lambda > LM_LAMBDA_MAX) break;
        }
    }

    for (int j = 0; j < NP; ++j)
        if (!std::isfinite(p_best[j])) return false;

    std::copy(p_best, p_best + NP, q);
    chi2_per_dof = (chi2_best * w_inv2) / static_cast<float>(std::max(1, nslice - NP));
    n_iter = iter;
    return true;
}

} // anon

WaveAnalyzer::PulseFitResult
WaveAnalyzer::FitPulseShape(const uint16_t *slice, int nslice,
                            int peak_idx_in_slice,
                            float ped, float ped_rms,
                            float clk_ns,
                            float model_err_floor)
{
    PulseFitResult res{};
    float q[TwoTau::NP];
    res.ok = fit_two_tau<TwoTau>(slice, nslice, peak_idx_in_slice, ped, ped_rms,
                                 clk_ns, model_err_floor, q, res.peak_amp,
                                 res.chi2_per_dof, res.n_iter);
    if (res.ok) {
        res.t0_ns    = q[0];
        res.tau_r_ns = q[1];
        res.tau_f_ns = q[2];
    }
    return res;
}

WaveAnalyzer::PulseFitTwoTauPResult
WaveAnalyzer::FitPulseShapeTwoTauP(const uint16_t *slice, int nslice,
                                    int peak_idx_in_slice,
                                    float ped, float ped_rms,
                                    float clk_ns,
                                    float model_err_floor)
{
    PulseFitTwoTauPResult res{};
    float q[TwoTauP::NP];
    res.ok = fit_two_tau<TwoTauP>(slice, nslice, peak_idx_in_slice, ped, ped_rms,
                                  clk_ns, model_err_floor, q, res.peak_amp,
                                  res.chi2_per_dof, res.n_iter);
    if (res.ok) {
        res.t0_ns    = q[0];
        res.tau_r_ns = q[1];
        res.tau_f_ns = q[2];
        res.p        = q[3];
    }
    return res;
}
