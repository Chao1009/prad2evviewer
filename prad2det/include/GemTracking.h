#pragma once
//=============================================================================
// GemTracking.h — straight-line track primitives (HyCal + GEM hits, lab frame)
//
// Header-only and templated on the scalar type: the online GEM efficiency
// monitor instantiates them in float, the Python bindings (gem_eff_audit.py)
// in double.  The fit sums are double in both cases.
//=============================================================================

#include <cmath>

namespace gem
{

template <class T>
struct TrackLine {
    T ax = 0, bx = 0;          // x(z) = ax + bx*z
    T ay = 0, by = 0;          // y(z) = ay + by*z
    T chi2_per_dof = 0;
};

// Line through two points.  |z2 - z1| < 1e-6 gives the flat line through
// point 1; callers pair points at different z (e.g. HyCal and a GEM).
template <class T>
inline TrackLine<T> SeedLine(T x1, T y1, T z1, T x2, T y2, T z2)
{
    TrackLine<T> L;
    T dz = z2 - z1;
    if (std::abs(dz) < T(1e-6)) { L.ax = x1; L.ay = y1; return L; }
    L.bx = (x2 - x1) / dz;  L.ax = x1 - L.bx * z1;
    L.by = (y2 - y1) / dz;  L.ay = y1 - L.by * z1;
    return L;
}

// Weighted least-squares fit v(z) = a + b*z.  False when singular.
template <class T>
inline bool FitWeightedLine1D(int N, const T *z, const T *v, const T *w,
                              double &a, double &b)
{
    double Sw = 0, Sz = 0, Szz = 0, Sv = 0, Svz = 0;
    for (int i = 0; i < N; ++i) {
        double wi = w[i];
        Sw  += wi;
        Sz  += wi * z[i];
        Szz += wi * z[i] * z[i];
        Sv  += wi * v[i];
        Svz += wi * v[i] * z[i];
    }
    double D = Sw * Szz - Sz * Sz;
    if (std::abs(D) < 1e-9) return false;
    b = (Sw * Svz - Sz * Sv) / D;
    a = (Sv - b * Sz) / Sw;
    return true;
}

// Independent weighted fits in (z, x) with weights wx and in (z, y) with wy;
// pass wx twice when sigma_x = sigma_y.  Distinct weights handle anisotropic
// points, e.g. a target constraint whose sigma_z couples into x and y through
// different slopes.  chi2_per_dof uses dof = 2N - 4 (0 when dof <= 0).
// False, with `out` untouched, when N < 2 or either fit is singular.
template <class T>
inline bool FitWeightedLine(int N, const T *z, const T *x, const T *y,
                            const T *wx, const T *wy, TrackLine<T> &out)
{
    if (N < 2) return false;
    double ax, bx, ay, by;
    if (!FitWeightedLine1D(N, z, x, wx, ax, bx) ||
        !FitWeightedLine1D(N, z, y, wy, ay, by))
        return false;
    out.ax = static_cast<T>(ax); out.bx = static_cast<T>(bx);
    out.ay = static_cast<T>(ay); out.by = static_cast<T>(by);
    int dof = 2 * N - 4;
    if (dof > 0) {
        double chi2 = 0;
        for (int i = 0; i < N; ++i) {
            double dxp = (ax + bx * z[i]) - x[i];
            double dyp = (ay + by * z[i]) - y[i];
            chi2 += wx[i] * dxp * dxp + wy[i] * dyp * dyp;
        }
        out.chi2_per_dof = static_cast<T>(chi2 / dof);
    } else {
        out.chi2_per_dof = 0;
    }
    return true;
}

} // namespace gem
