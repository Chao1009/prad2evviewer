#!/usr/bin/env python3
"""mf_pulse.py — matched-filter pulse reconstruction for HyCal PbWO4 FADC waveforms.

Pure NumPy (no prad2 import).  Used by benchmark_hycal_timing_test_20261021.py
with --pulse-method mf; see plan_20261022.md §4 and its "Decisions (user,
2026-10-22)" section.

Model F4:  p(t) = (1 - exp(-t/tau_r)) exp(-t/tau_f) / T_max   (t > 0, else 0),
unit peak.  Per-module (tau_r, tau_f) come from output/pulse_templates_h50.dat
(25320 column); invalid / missing PbWO4 modules use the PbWO4 aggregate.

Reconstruction (MFReconstructor.run), batched over any number of channels:
  1. baseline = clipped mean of samples 0..ped_nsamples-1, sigma = 1.4826 MAD
  2. candidate channel: max 3-sample moving average > max(3 sigma, 3 ADC)
  3. matched filter over onset t0 = 4 j - phi (phi = 0..3 ns), L = 24 samples
  4. accept the best (j, phi) if A >= max(nsigma * sigma_A, A_min); a secondary
     inside an accepted peak's kernel span also needs A >= r_model * A_parent
  5. parabola on the 1 ns score grid + 3 Gauss-Newton steps on (a, t0)
  6. subtract, repeat up to kmax peaks per channel
  7. joint K x K amplitude fit (positivity drop), one (a, t0) refinement pass
     for K >= 2, chi2/dof with a 0.03 model-error floor
  8. clipped samples (>= overflow) are masked (weight 0); channels whose peak
     window has < 6 unclipped samples or < 2 on the leading edge are flagged
     for a WaveAnalyzer fallback (done by the caller).

Banks: TemplateBank (two-tau, above) or SampledTemplateBank (per-module
sampled shapes from a spline_templates_v1 .npz, via build_bank on a ".npz"
path).  Both provide shape(m, t) / shape_and_slope(m, t), K, t50, t_pk, I_tab,
iint, so MFReconstructor does not care which one it holds.

Energy: E_units = A * I_int[m](phase) * kappa [* f(A)], where I_int is the WaveAnalyzer
walk-out integral (WaveAnalyzer.cpp:391-462, smoothing :235) of the noiseless
sample-aligned unit template, in the high-amplitude limit (cut = int_tail_ratio
* smoothed peak, tail_break_n consecutive sub-cut samples end the walk).
The sampled bank (iint_mode "continuous", default) uses the phase average of that
walk-out over the continuous template (I_cont = I_mean; no sub-sample phase
dependence; continuous_walk_integral is kept as the I_cont_walk diagnostic) and may
carry a global amplitude factor f(A) (ampcorr keys of the .npz); see SampledTemplateBank.

    python3 mf_pulse.py --self-test
"""

from __future__ import annotations

import hashlib
import json
import math
import os
import sys

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

DT = 4.0                       # ns per sample
TAU_R_RANGE = (1.0, 6.0)       # validity of a per-module pair (ns)
TAU_F_RANGE = (15.0, 35.0)
NPHASE = 4                     # kernel onset phases (1 ns grid)
NIPH = 256                     # I_int phase grid points per 4 ns sample (+1 wrap)
MODEL_FLOOR = 0.03
LW = 32                        # fit window: samples [ceil(t0/4)-3, +LW)
LEAD = 3
MIN_SEP_NS = 12.0              # minimum onset separation of two pulses in a channel
SIGMA_MIN = 0.3                # ADC floor on the raw sigma
_TBL_N, _TBL_I0 = 64, 12       # I_int table waveform length / first-sample index


# ---------------------------------------------------------------- template --
def tmax_of(tr, tf):
    u = tr / (tr + tf)
    return (1.0 - u) * u ** (tr / tf)


def template(t, tr, tf, tm):
    """Unit-peak two-tau pulse at times t (ns) relative to onset (broadcasts)."""
    tp = np.maximum(t, 0.0)
    return (1.0 - np.exp(-tp / tr)) * np.exp(-tp / tf) / tm


def template_and_slope(t, tr, tf, tm):
    """(s, ds/dt) of the template; ds/dt = 0 for t <= 0."""
    tp = np.maximum(t, 0.0)
    e1 = np.exp(-tp / tr)
    e2 = np.exp(-tp / tf)
    s = (1.0 - e1) * e2 / tm
    ds = np.where(t > 0.0, (e1 / tr - (1.0 - e1) / tf) * e2 / tm, 0.0)
    return s, ds


def smooth_tri(x, order):
    """Triangular smoothing of WaveAnalyzer::smooth (edges left unsmoothed)."""
    x = np.asarray(x, dtype=np.float64)
    if order <= 1:
        return x.copy()
    n = x.shape[-1]
    val = x.copy()
    wsum = np.ones(n)
    for j in range(1, order):
        if 2 * j >= n:
            break
        w = 1.0 - j / float(order + 1)
        val[..., j:n - j] += w * (x[..., :n - 2 * j] + x[..., 2 * j:])
        wsum[j:n - j] += 2.0 * w
    return val / wsum


def walk_integral(v, ratio, nbreak, floor=None):
    """WaveAnalyzer outward walk on unimodal smoothed rows v[R, N] (peak > 0).

    A sample is 'below' if v < max(ratio * peak, floor) or v < 0; nbreak
    consecutive below samples end the walk and the pending ones are discarded
    (WaveAnalyzer.cpp:391-462).  Returns the integral of the kept samples.
    """
    R, N = v.shape
    ar = np.arange(R)
    pk = np.argmax(v, axis=1)
    h = v[ar, pk]
    cut = ratio * h
    if floor is not None:
        cut = np.maximum(cut, floor)
    integral = h.copy()
    nb = max(1, int(nbreak))
    for side in (-1, 1):
        run = np.zeros(R, dtype=np.int64)
        pend = np.zeros(R)
        active = np.ones(R, dtype=bool)
        for d in range(1, N):
            pos = pk + side * d
            inr = active & (pos >= 0) & (pos < N)
            if not inr.any():
                break
            vv = v[ar, np.clip(pos, 0, N - 1)]
            below = (vv < cut) | (vv < 0.0)
            mb = inr & below
            run = np.where(mb, run + 1, run)
            pend = np.where(mb, pend + vv, pend)
            active &= ~(mb & (run >= nb))
            ma = inr & ~below
            integral += np.where(ma, pend + vv, 0.0)
            pend = np.where(ma, 0.0, pend)
            run = np.where(ma, 0, run)
    return integral


def unit_integral(tr, tf, tm, u, smooth_order, ratio, nbreak, floor=None):
    """I_int of unit templates sampled at tau_i = 4 (i - i0) + u.

    tr, tf, tm, u (and floor, in unit-peak units) are 1-D arrays of equal length.
    """
    tau = DT * (np.arange(_TBL_N) - _TBL_I0)[None, :] + np.asarray(u)[:, None]
    s = template(tau, np.asarray(tr)[:, None], np.asarray(tf)[:, None],
                 np.asarray(tm)[:, None])
    return walk_integral(smooth_tri(s, smooth_order), ratio, nbreak, floor)


def unit_integral_fn(shape_fn, u, smooth_order, ratio, nbreak, floor=None):
    """unit_integral for any shape: shape_fn maps times [R, N] (ns) to values [R, N].

    The caller closes over the row -> module mapping; u (and floor) are 1-D, length R.
    """
    tau = DT * (np.arange(_TBL_N) - _TBL_I0)[None, :] + np.asarray(u)[:, None]
    return walk_integral(smooth_tri(shape_fn(tau), smooth_order), ratio, nbreak, floor)


def _cont_smooth(shape_fn, smooth_order, dt, t_max, nrow):
    """Continuous triangle-smoothed template on a fine grid: S[nrow, n].

    S(t) = sum_j w_j s(t + DT j) / sum_j w_j, j = -(order-1)..(order-1),
    w_j = 1 - |j|/(order+1) (smooth_tri); t = -DT (order-1) + dt i up to t_max.
    s is evaluated once on an extended grid (DT / dt must be an integer).
    """
    q = int(round(DT / dt))
    if abs(q * dt - DT) > 1e-9:
        raise ValueError("dt must divide the 4 ns sample")
    order = max(int(smooth_order), 1)
    k = q * (order - 1)                       # index shift of the outermost term
    n = int(round((t_max + DT * (order - 1)) / dt)) + 1
    t = dt * (np.arange(n + 2 * k) - k) - DT * (order - 1)
    s = np.asarray(shape_fn(np.broadcast_to(t, (nrow, t.size))), dtype=np.float64)
    if order == 1:
        return s
    S = s[:, k:k + n].copy()
    wsum = 1.0
    for j in range(1, order):
        w = 1.0 - j / float(order + 1)
        S += w * (s[:, k - q * j:k - q * j + n] + s[:, k + q * j:k + q * j + n])
        wsum += 2.0 * w
    return S / wsum


def _cont_walk(S, ratio, nbreak, floor, dt):
    """Walk-out integral of the continuous smoothed rows S[R, n] -> (I, peak).

    Continuous analogue of walk_integral: cut = max(ratio * peak, floor);
    the kept region is the connected set around the peak where S >= cut, in
    which a below-cut stretch (or S < 0) shorter than (nbreak - 0.5) DT is
    bridged (kept, below-cut values included) and a stretch of at least that
    length ends the walk on that side; a trailing stretch with nothing above
    the cut after it is dropped.  The discrete rule ends the walk on nbreak
    consecutive below-cut samples, i.e. on a stretch holding >= nbreak samples;
    a stretch of length L holds nbreak samples with probability 0 for
    L <= (nbreak - 1) DT and 1 for L >= nbreak DT, 1/2 at (nbreak - 0.5) DT,
    so this is the phase average of that rule.  I = integral of S over the
    kept region / DT (trapezoid on the dt grid): the phase average of the
    discrete sample sum.  A cut above the peak gives I = peak (one sample).
    """
    R, n = S.shape
    ar = np.arange(R)
    pk = np.argmax(S, axis=1)
    h = S[ar, pk]
    cut = ratio * h
    if floor is not None:
        cut = np.maximum(cut, floor)
    below = (S < cut[:, None]) | (S < 0.0)
    Lp = max(1, int(round((max(1, int(nbreak)) - 0.5) * DT / dt)))
    idx = np.arange(n, dtype=np.int32)[None, :]

    def end(bel, p):
        la = np.maximum.accumulate(np.where(bel, np.int32(-1), idx), axis=1)
        hit = bel & ((idx - la) >= Lp) & (idx > p[:, None])
        i1 = np.where(hit.any(axis=1), hit.argmax(axis=1), n - 1)
        return np.clip(la[ar, i1], 0, n - 1)

    r = end(below, pk)
    l = (n - 1) - end(below[:, ::-1], (n - 1) - pk)
    cs = np.cumsum(S, axis=1)
    I = (cs[ar, r] - cs[ar, l] + 0.5 * (S[ar, l] - S[ar, r])) * (dt / DT)
    return np.where(cut > h, h, I), h


def continuous_walk_integral(shape_fn, smooth_order, ratio, nbreak, floor=None,
                             dt=0.01, t_max=400.0, nrow=None):
    """Walk-out I_int on the continuous smoothed template (no sub-sample phase).

    shape_fn maps times [R, N] (ns after onset) to values [R, N] as for
    unit_integral_fn; R = len(floor) if floor (unit-peak units, per row) is
    given, else nrow.  Returns I[R]: the phase average of unit_integral_fn
    (see _cont_walk for the continuous walk rule).  Keep R modest (~32): the
    working arrays are R x (t_max + 4 (order-1)) / dt.
    """
    R = len(floor) if floor is not None else nrow
    if R is None:
        raise ValueError("give floor or nrow")
    S = _cont_smooth(shape_fn, smooth_order, dt, t_max, R)
    fl = None if floor is None else np.asarray(floor, dtype=np.float64)
    return _cont_walk(S, ratio, nbreak, fl, dt)[0]


# -------------------------------------------------------------- h50 loader --
AGG_FALLBACK = (2.499476432800293, 24.406997680664062)


def table_columns(path):
    """Run columns of a dump_pbwo4_table.py table, in header order."""
    cols = []
    with open(path) as fh:
        for line in fh:
            if line.startswith("# module"):
                for c in line[2:].split("\t"):
                    c = c.strip()
                    if c.startswith("tau_r_mean_"):
                        cols.append(c[len("tau_r_mean_"):])
                break
    return cols


def resolve_column(path, column="auto", prefer="25320"):
    """`column` as given, or `prefer` if the table has it, else its last column.

    "sampled" for a spline_templates_v1 .npz (no run columns).
    """
    if str(path).endswith(".npz"):
        return "sampled"
    if column not in (None, "auto"):
        return str(column)
    cols = table_columns(path)
    if not cols:
        raise ValueError(f"{path}: no tau_r_mean_<run> columns")
    return prefer if prefer in cols else cols[-1]


def load_h50(path, column="25320"):
    """Read a per-module (tau_r, tau_f) table -> ({name: (tau_r, tau_f)}, sha1, json).

    The `column` run's median is used.  Only the h50 table (columns 25308 and
    25320, column 25320) has a secondary column: a NaN pair there takes the
    25308 value; any other table falls back to its own aggregate.  Rows failing
    the validity rule (1<=tau_r<=6, 15<=tau_f<=35 ns) are dropped, as are modules the
    run's JSON does not record as good_fit or with n_pulses_used >= min_pulses (the
    table's own parameter block).  `json` is the
    aggregate file named by the table's own `# Run #<column>: <json>` header line.
    """
    with open(path, "rb") as fh:
        raw = fh.read()
    sha1 = hashlib.sha1(raw).hexdigest()
    column = str(column)
    lines = raw.decode().splitlines()
    jpath, min_pulses, runs = None, 50, []
    for line in lines:
        if line.startswith(f"# Run #{column}:"):
            jpath = line.split(":", 1)[1].strip()
        t = line[1:].split()
        if t[:1] == ["Parameter"]:                 # "Parameter  Run A  Run B"
            runs = [x for x in t[1:] if x != "Run"]
        elif t[:1] == ["min_pulses"] and column in runs:
            try:
                min_pulses = int(t[1 + runs.index(column)])
            except (ValueError, IndexError):
                pass
    if jpath is None:
        raise ValueError(f"{path}: no '# Run #{column}: <json>' header line")
    if not os.path.isfile(jpath):
        alt = os.path.join(os.path.dirname(os.path.abspath(path)),
                           os.path.basename(jpath))
        if os.path.isfile(alt):
            jpath = alt
    # a column value is trusted only if the run's JSON has the channel as
    # good_fit or with n_pulses_used >= min_pulses (the table has no such filter)
    with open(jpath) as fh:
        jch = json.load(fh)
    trusted = {k for k, e in jch.items()
               if isinstance(e, dict) and (e.get("good_fit")
                                           or e.get("n_pulses_used", 0) >= min_pulses)}
    cols = None
    rows = {}
    other = None
    for line in lines:
        if line.startswith("# module"):
            cols = [c.strip() for c in line[2:].split("\t")]
            if f"tau_r_mean_{column}" not in cols:
                raise ValueError(f"{path}: no column {column} (have "
                                 f"{table_columns(path)})")
            if column == "25320" and "tau_r_mean_25308" in cols:
                other = "25308"
            continue
        if not line.strip() or line.startswith("#"):
            continue
        if cols is None:
            raise ValueError(f"{path}: data before the column header")
        f = line.split()
        d = dict(zip(cols, f))
        vals = []
        for k in ("tau_r_mean_", "tau_f_mean_"):
            v = float(d[k + column])
            if not math.isfinite(v) and other:
                v = float(d[k + other])
            elif f[0] not in trusted:
                v = float("nan")
            vals.append(v)
        tr, tf = vals
        if (math.isfinite(tr) and math.isfinite(tf)
                and TAU_R_RANGE[0] <= tr <= TAU_R_RANGE[1]
                and TAU_F_RANGE[0] <= tf <= TAU_F_RANGE[1]):
            rows[f[0]] = (tr, tf)
    return rows, sha1, jpath


def load_aggregate(json_path):
    """PbWO4 aggregate (tau_r, tau_f) from a template JSON's `_by_type`."""
    with open(json_path) as fh:
        d = json.load(fh)["_by_type"]["PbWO4"]
    return (float(d["tau_r_ns"]["median"]), float(d["tau_f_ns"]["median"]))


# ------------------------------------------------------------ template bank --
class TemplateBank:
    """Per-module kernel bank, t50 and I_int tables (indexed by module index)."""

    kind = "two_tau"

    def __init__(self, tau_r, tau_f, pwo4, *, smooth_order, int_tail_ratio,
                 tail_break_n, L=24):
        self.names = None              # build_bank sets the module names
        self.fallback = []             # names on the aggregate (build_bank)
        self.T_end = 400.0             # ns, unused: the analytic shape has no end
        self.tau_r = np.asarray(tau_r, dtype=np.float64)
        self.tau_f = np.asarray(tau_f, dtype=np.float64)
        self.pwo4 = np.asarray(pwo4, dtype=bool)
        nm = self.tau_r.size
        self.L = L
        self.smooth_order = int(smooth_order)
        self.int_tail_ratio = float(int_tail_ratio)
        self.tail_break_n = int(tail_break_n)
        # non-PbWO4 rows hold the aggregate so the arrays stay finite
        self.tau_r = np.where(self.pwo4, self.tau_r, AGG_FALLBACK[0])
        self.tau_f = np.where(self.pwo4, self.tau_f, AGG_FALLBACK[1])
        self.tmx = tmax_of(self.tau_r, self.tau_f)
        self.t_pk = self.tau_r * np.log1p(self.tau_f / self.tau_r)
        # kernel bank K[m, phi, l] = p(4 l + phi)
        tl = (DT * np.arange(L)[None, None, :]
              + np.arange(NPHASE)[None, :, None])
        K = template(tl, self.tau_r[:, None, None], self.tau_f[:, None, None],
                     self.tmx[:, None, None])
        self.K = K.astype(np.float32)
        self.K2 = (K * K).astype(np.float32)
        self.Kn2 = (K * K).sum(axis=2).astype(np.float32)          # [nm, 4]
        self.Kinv = (1.0 / np.sqrt(np.maximum(self.Kn2, 1e-12))).astype(np.float32)
        # leading-edge 50% time relative to onset
        g = np.linspace(0.0, 1.0, 4001)[None, :] * self.t_pk[:, None]
        sg = template(g, self.tau_r[:, None], self.tau_f[:, None], self.tmx[:, None])
        i1 = np.argmax(sg >= 0.5, axis=1)
        i0 = np.maximum(i1 - 1, 0)
        ar = np.arange(nm)
        s0, s1 = sg[ar, i0], sg[ar, i1]
        fr = np.where(s1 > s0, (0.5 - s0) / np.maximum(s1 - s0, 1e-300), 0.0)
        self.t50 = g[ar, i0] + fr * (g[ar, i1] - g[ar, i0])
        self.build_iint(int_tail_ratio)

    def shape(self, m, t, dtype=np.float32):
        """Unit-peak shape at module index m, t ns after onset (m, t broadcast).

        The module parameters are cast to `dtype` (float32 as the reconstructor
        always did; float64 for the synthetic generator).
        """
        return template(t, self.tau_r[m].astype(dtype), self.tau_f[m].astype(dtype),
                        self.tmx[m].astype(dtype))

    def shape_and_slope(self, m, t, dtype=np.float32):
        """(shape, d shape/dt); both 0 for t <= 0."""
        return template_and_slope(t, self.tau_r[m].astype(dtype),
                                  self.tau_f[m].astype(dtype),
                                  self.tmx[m].astype(dtype))

    def build_iint(self, ratio):
        """I_int[m, phase-grid] (NIPH+1 points over one 4 ns sample)."""
        self.int_tail_ratio = float(ratio)
        nm = self.tau_r.size
        u = DT * np.arange(NIPH + 1) / NIPH
        tr = np.repeat(self.tau_r, u.size)
        tf = np.repeat(self.tau_f, u.size)
        tm = np.repeat(self.tmx, u.size)
        uu = np.tile(u, nm)
        tab = unit_integral(tr, tf, tm, uu, self.smooth_order, self.int_tail_ratio,
                            self.tail_break_n).reshape(nm, u.size)
        self.I_tab = tab.astype(np.float32)
        self.I_mean = tab[:, :NIPH].mean(axis=1)
        self.area = np.array([
            template(DT * np.arange(400), a, b, c).sum()
            for a, b, c in zip(self.tau_r, self.tau_f, self.tmx)])

    def iint(self, m, t0):
        """I_int at module index m for onset t0 (ns): sample-aligned phase."""
        u = np.mod(-np.asarray(t0, dtype=np.float64), DT) * (NIPH / DT)
        k = np.minimum(np.floor(u).astype(np.int64), NIPH - 1)
        fr = u - k
        return (self.I_tab[m, k] * (1.0 - fr) + self.I_tab[m, k + 1] * fr)

    def iint_floor(self, m, t0, amp, ped_rms):
        """I_int with the ped_rms floor (v < max(ratio*h, ped_rms)), per pulse."""
        u = np.mod(-np.asarray(t0, dtype=np.float64), DT)
        fl = np.asarray(ped_rms, dtype=np.float64) / np.maximum(amp, 1e-6)
        return unit_integral(self.tau_r[m], self.tau_f[m], self.tmx[m], u,
                             self.smooth_order, self.int_tail_ratio,
                             self.tail_break_n, floor=fl)


_LUT_DT = 0.05                  # ns, lookup step of the sampled bank


class _ClippedSpline:
    """BSpline callable that holds its end value beyond the base interval."""

    def __init__(self, spl):
        self.spl = spl
        self.tmax = float(spl.t[-spl.k - 1])

    def __call__(self, x):
        return self.spl(np.minimum(x, self.tmax))

    def derivative(self):
        return _ClippedSpline(self.spl.derivative())


_CONT_DT, _CONT_TMAX, _CONT_CHUNK = 0.01, 400.0, 32   # continuous walk-out grid


class SampledTemplateBank:
    """TemplateBank interface on per-module sampled shapes (spline_templates_v1).

    shape_rows: one entry per module (bank order) -- a callable f(t) of ns after
    onset (unit peak; a `.derivative()` method gives the slope, else it is
    differenced) or an array of values on t = 0, 0.05, ..., T_end.  Tables
    lut/dlut[nm, n] (float32, step 0.05 ns) hold value and slope; shape() is
    their linear interpolation, 0 outside (0, T_end).  Non-PbWO4 rows hold
    whatever shape the caller gives them (from_npz: the average).

    Energy integral: iint_mode "continuous" (default) takes I_cont[nm] = I_mean, the
    phase average of the discrete walk-out over the continuous template (no
    sub-sample phase dependence); "phase" keeps the 256-phase table I_tab.
    I_tab and the strictly continuous walk I_cont_walk are kept as diagnostics.
    ampcorr = (log A, f) (optional, from_npz) is a global amplitude-dependent
    integral factor applied by pulse_energy_units through f_amp(amp).  f may also
    be a per-module table [n_modules, len(log A)] (rows = this bank's modules; from_npz
    maps the file's per-module ampcorr_f to them); f_amp(amp, midx) then takes the
    module of every pulse.
    """

    kind = "sampled"

    def __init__(self, shape_rows, names, pwo4, *, smooth_order, int_tail_ratio,
                 tail_break_n, L=24, T_end=160.0, fallback=(), meta=None,
                 iint_mode="continuous", ampcorr=None, ampcorr_meta=None):
        if iint_mode not in ("continuous", "phase"):
            raise ValueError(f"iint_mode must be 'continuous' or 'phase', got {iint_mode!r}")
        self.iint_mode = iint_mode
        if ampcorr is not None:
            la, fa = (np.asarray(x, dtype=np.float64) for x in ampcorr)
            per_module = fa.ndim == 2
            if (la.ndim != 1 or la.size < 1
                    or (la.shape != fa.shape and not (
                        per_module and fa.shape == (len(names), la.size)))
                    or not (np.isfinite(la).all() and np.isfinite(fa).all())
                    or (np.diff(la) <= 0).any()):
                raise ValueError("ampcorr: need equal-length 1-D finite arrays (f may be "
                                 "[n_modules, len(log A)]), log A strictly increasing")
            ampcorr = (la, fa)
        self.ampcorr = ampcorr
        self.ampcorr_meta = dict(ampcorr_meta or {})
        self.names = list(names)
        self.pwo4 = np.asarray(pwo4, dtype=bool)
        nm = len(self.names)
        if len(shape_rows) != nm or self.pwo4.size != nm:
            raise ValueError("shape_rows, names and pwo4 must have equal length")
        self.L = L
        self.smooth_order = int(smooth_order)
        self.int_tail_ratio = float(int_tail_ratio)
        self.tail_break_n = int(tail_break_n)
        self.T_end = float(T_end)
        self.fallback = list(fallback)
        self.n_fallback = len(self.fallback)
        self.n_file = int(self.pwo4.sum()) - self.n_fallback   # rows with own shape
        self.meta = dict(meta or {})
        self.run_type = self.meta.get("run_type")
        n = self.n_lut = int(round(self.T_end / _LUT_DT)) + 1
        tg = _LUT_DT * np.arange(n)
        self.lut = np.empty((nm, n), dtype=np.float32)
        self.dlut = np.empty((nm, n), dtype=np.float32)
        cache = {}                     # rows sharing one object (the average)
        for i, f in enumerate(shape_rows):
            if id(f) not in cache:
                if callable(f):
                    v = np.asarray(f(tg), dtype=np.float64)
                    d = (f.derivative()(tg) if hasattr(f, "derivative")
                         else np.gradient(v, _LUT_DT))
                else:
                    v = np.asarray(f, dtype=np.float64)
                    if v.shape != (n,):
                        raise ValueError(f"row {self.names[i]}: need {n} samples")
                    d = np.gradient(v, _LUT_DT)
                cache[id(f)] = (v.astype(np.float32),
                                np.asarray(d, dtype=np.float64).astype(np.float32))
            self.lut[i], self.dlut[i] = cache[id(f)]
        self._lut = self.lut.reshape(-1)
        self._dlut = self.dlut.reshape(-1)
        mi = np.arange(nm)
        # kernel bank K[m, phi, l] = p(4 l + phi)
        tl = (DT * np.arange(L)[None, None, :]
              + np.arange(NPHASE)[None, :, None])
        K = self.shape(mi[:, None, None], tl, dtype=np.float64)
        self.K = K.astype(np.float32)
        self.K2 = (K * K).astype(np.float32)
        self.Kn2 = (K * K).sum(axis=2).astype(np.float32)          # [nm, 4]
        self.Kinv = (1.0 / np.sqrt(np.maximum(self.Kn2, 1e-12))).astype(np.float32)
        self.t_pk = np.argmax(self.lut, axis=1) * _LUT_DT
        # leading-edge 50% time relative to onset (grid search as TemplateBank)
        g = np.linspace(0.0, 1.0, 4001)[None, :] * self.t_pk[:, None]
        sg = self.shape(mi[:, None], g, dtype=np.float64)
        i1 = np.argmax(sg >= 0.5, axis=1)
        i0 = np.maximum(i1 - 1, 0)
        s0, s1 = sg[mi, i0], sg[mi, i1]
        fr = np.where(s1 > s0, (0.5 - s0) / np.maximum(s1 - s0, 1e-300), 0.0)
        self.t50 = g[mi, i0] + fr * (g[mi, i1] - g[mi, i0])
        self.build_iint(int_tail_ratio)

    def _locate(self, m, t, dtype):
        t = np.asarray(t, dtype=dtype)
        x = t * dtype(1.0 / _LUT_DT)
        i = np.clip(np.floor(x), 0, self.n_lut - 2)
        f = x - i
        base = np.asarray(m, dtype=np.int64) * self.n_lut + i.astype(np.int64)
        return base, f, (t > 0) & (t < self.T_end)

    def _interp(self, tab, base, f, inside, dtype):
        v = tab[base].astype(dtype) * (1 - f) + tab[base + 1].astype(dtype) * f
        return np.where(inside, v, dtype(0))

    def shape(self, m, t, dtype=np.float32):
        """Unit-peak shape at module index m, t ns after onset (m, t broadcast).

        dtype (np.float32 / np.float64) is the arithmetic of the interpolation.
        """
        base, f, inside = self._locate(m, t, dtype)
        return self._interp(self._lut, base, f, inside, dtype)

    def shape_and_slope(self, m, t, dtype=np.float32):
        """(shape, d shape/dt); both 0 for t <= 0 and t >= T_end."""
        base, f, inside = self._locate(m, t, dtype)
        return (self._interp(self._lut, base, f, inside, dtype),
                self._interp(self._dlut, base, f, inside, dtype))

    def build_iint(self, ratio):
        """I_int[m, phase-grid] (NIPH+1 points over one 4 ns sample)."""
        self.int_tail_ratio = float(ratio)
        nm = len(self.names)
        u = DT * np.arange(NIPH + 1) / NIPH
        tab = np.empty((nm, u.size))
        for a in range(0, nm, 128):                  # bounded memory
            mrow = np.repeat(np.arange(a, min(a + 128, nm)), u.size)
            tab[a:a + 128] = unit_integral_fn(
                lambda tau: self.shape(mrow[:, None], tau, dtype=np.float64),
                np.tile(u, mrow.size // u.size), self.smooth_order,
                self.int_tail_ratio, self.tail_break_n).reshape(-1, u.size)
        self.I_tab = tab.astype(np.float32)
        self.I_mean = tab[:, :NIPH].mean(axis=1)
        self.area = self.shape(np.arange(nm)[:, None], DT * np.arange(400)[None, :],
                               dtype=np.float64).sum(axis=1)
        # continuous walk-out integral; rows with an identical shape computed once
        _, first, inv = np.unique(self.lut, axis=0, return_index=True, return_inverse=True)
        Iu, hu = np.empty(len(first)), np.empty(len(first))
        for a in range(0, len(first), _CONT_CHUNK):
            mu = np.arange(a, min(a + _CONT_CHUNK, len(first)))
            mrow = first[mu]
            S = _cont_smooth(lambda tau: self.shape(mrow[:, None], tau, dtype=np.float64),
                             self.smooth_order,
                             _CONT_DT, _CONT_TMAX, mu.size)
            Iu[mu], hu[mu] = _cont_walk(S, self.int_tail_ratio, self.tail_break_n,
                                        None, _CONT_DT)
        # I_cont_walk: walk-out on the continuous smoothed template (cut from the true
        # smoothed peak); it sits 0.1-0.7 % below the phase average of the discrete
        # rule because the discrete cut uses the sample peak.  The energy integral
        # I_cont is the phase average of the exact discrete walk-out (I_mean): no
        # sub-sample phase dependence and equal to the C++ rule averaged over phase.
        self.I_cont_walk = Iu[inv.reshape(-1)]
        self.I_cont = np.asarray(self.I_mean, dtype=np.float64).copy()
        self.h_cont = hu[inv.reshape(-1)]          # peak of the smoothed unit template

    def iint(self, m, t0):
        """I_int at module index m for onset t0 (ns); continuous mode: no t0 dependence."""
        if self.iint_mode == "phase":
            return TemplateBank.iint(self, m, t0)
        return self.I_cont[np.asarray(m, dtype=np.int64)] + np.zeros(np.shape(t0))

    def iint_floor(self, m, t0, amp, ped_rms):
        """I_int with the ped_rms floor (v < max(ratio*h, ped_rms)), per pulse.

        Continuous mode: pulses whose floor ped_rms/amp stays below ratio * (smoothed
        peak) have the unfloored cut and take I_cont; the others are walked on the
        continuous template with that floor (slower, ~1 ms per pulse).
        """
        mm = np.asarray(m, dtype=np.int64)
        fl = np.asarray(ped_rms, dtype=np.float64) / np.maximum(amp, 1e-6)
        if self.iint_mode == "phase":
            u = np.mod(-np.asarray(t0, dtype=np.float64), DT)
            return unit_integral_fn(
                lambda tau: self.shape(mm[:, None], tau, dtype=np.float64), u,
                self.smooth_order, self.int_tail_ratio, self.tail_break_n, floor=fl)
        mm, fl = np.broadcast_arrays(mm, fl)
        shp = mm.shape
        mm, fl = mm.reshape(-1), fl.reshape(-1)
        out = self.I_cont[mm].copy()
        need = np.flatnonzero(fl > self.int_tail_ratio * self.h_cont[mm])
        for a in range(0, need.size, _CONT_CHUNK):
            sel = need[a:a + _CONT_CHUNK]
            um, back = np.unique(mm[sel], return_inverse=True)
            S = _cont_smooth(lambda tau: self.shape(um[:, None], tau, dtype=np.float64),
                             self.smooth_order, _CONT_DT, _CONT_TMAX, um.size)
            # floored walk, rescaled by I_cont / I_cont_walk so floored and unfloored
            # pulses share the I_cont normalisation
            out[sel] = (_cont_walk(S[back.reshape(-1)], self.int_tail_ratio,
                                   self.tail_break_n, fl[sel], _CONT_DT)[0]
                        * self.I_cont[mm[sel]] / self.I_cont_walk[mm[sel]])
        return out.reshape(shp)

    def f_amp(self, amp, midx=None):
        """Amplitude-dependent integral factor f(A): log-linear in A between the
        ampcorr knots, constant outside; 1 if the bank carries none.  A per-module
        table (2-D f) needs `midx`, the bank module of every pulse; its rows are
        interpolated like np.interp does the 1-D table."""
        amp = np.asarray(amp, dtype=np.float64)
        if self.ampcorr is None:
            return np.ones(amp.shape)
        if self.ampcorr[1].ndim == 1:
            return np.interp(np.log(np.maximum(amp, 1.0)), *self.ampcorr)
        if midx is None:
            raise ValueError("f_amp: the bank carries a per-module f(A) table; "
                             "pass the module index of every pulse")
        la, F = self.ampcorr
        x = np.log(np.maximum(amp, 1.0))
        m = np.broadcast_to(np.asarray(midx), x.shape)
        j = np.clip(np.searchsorted(la, x, side="right") - 1, 0, la.size - 2) \
            if la.size > 1 else np.zeros(x.shape, dtype=np.intp)
        f0 = F[m, j]
        if la.size == 1:
            return f0
        f1 = F[m, j + 1]
        out = (f1 - f0) / (la[j + 1] - la[j]) * (x - la[j]) + f0
        return np.where(x <= la[0], F[m, 0], np.where(x >= la[-1], F[m, -1], out))

    def ampcorr_str(self):
        """Short description of ampcorr (None if absent)."""
        if self.ampcorr is None:
            return None
        la, fa = self.ampcorr
        return (f"f(A) {'per-module ' if fa.ndim == 2 else ''}{la.size} points, "
                f"A {math.exp(la[0]):.0f}-{math.exp(la[-1]):.0f} ADC, "
                f"f {fa.min():.4f}-{fa.max():.4f}")

    @classmethod
    def from_npz(cls, path, names, pwo4, *, smooth_order, int_tail_ratio,
                 tail_break_n, L=24, iint_mode="continuous"):
        """Bank from a spline_templates_v1 .npz, rows matched to `names` by name.

        Each file row is a unit-peak cubic B-spline on the onset-shifted axis.  A
        PbWO4 name missing from the file or with valid == 0 takes the average
        spline; those and the file's own fallback == 2 rows are listed in
        `fallback`.  Non-PbWO4 rows hold the average.  Optional keys ampcorr_logA,
        ampcorr_f (both or neither; ln A in ADC, multiplicative factors) and
        ampcorr_meta (JSON string) give bank.ampcorr and bank.f_amp.  A 2-D
        ampcorr_f [file rows, len(ampcorr_logA)] (`ampcorr --per-module-f`) is a
        per-module table: bank row = the file row of the same name; names without a
        file row (non-PbWO4, absent) take ampcorr_f_rt (the run-type f on the same
        grid) if the file has it, else the median over the file rows.  The band
        arrays of the file (lo / mid / top, or the generalised band_edges /
        shape_band / r_band / n_band of `build --bands`) are diagnostics and are
        not read.
        """
        from scipy.interpolate import BSpline
        z = np.load(path, allow_pickle=False)
        k = int(z["spline_k"])
        T_end = float(z["T_end"])

        def spline(t, c):
            s = _ClippedSpline(BSpline(np.asarray(t, dtype=np.float64),
                                       np.asarray(c, dtype=np.float64), k))
            if s.tmax < T_end - 1e-6:
                raise ValueError(f"{path}: spline ends at {s.tmax} ns < T_end {T_end}")
            return s

        avg_t, avg_c = z["avg_spline_t"], z["avg_spline_c"]
        avg = spline(avg_t[np.isfinite(avg_t)], avg_c[np.isfinite(avg_c)])
        row = {str(nm): r for r, nm in enumerate(z["names"])}
        valid, fbs = z["valid"], z["fallback"]
        nt, nc = z["spline_nt"], z["spline_nc"]
        rows, fallback = [], []
        for nm, is_p in zip(names, pwo4):
            r = row.get(nm) if is_p else None
            if r is not None and valid[r] > 0:
                rows.append(spline(z["spline_t"][r, :nt[r]], z["spline_c"][r, :nc[r]]))
                if fbs[r] == 2:
                    fallback.append(nm)
            else:
                rows.append(avg)
                if is_p:
                    fallback.append(nm)
        meta = json.loads(str(z["meta"])) if "meta" in z.files else {}
        has = ("ampcorr_logA" in z.files, "ampcorr_f" in z.files)
        if has[0] != has[1]:
            raise ValueError(f"{path}: ampcorr_logA and ampcorr_f must come together")
        ampcorr = (z["ampcorr_logA"], z["ampcorr_f"]) if has[0] else None
        if has[0] and ampcorr[1].ndim == 2:
            zf = ampcorr[1]
            dflt = (z["ampcorr_f_rt"] if "ampcorr_f_rt" in z.files
                    else np.median(zf, axis=0))
            ampcorr = (ampcorr[0], np.array(
                [zf[row[nm]] if is_p and nm in row else dflt
                 for nm, is_p in zip(names, pwo4)]))
        acm = json.loads(str(z["ampcorr_meta"])) if "ampcorr_meta" in z.files else {}
        return cls(rows, names, pwo4, smooth_order=smooth_order,
                   int_tail_ratio=int_tail_ratio, tail_break_n=tail_break_n, L=L,
                   T_end=T_end, fallback=fallback, meta=meta, iint_mode=iint_mode,
                   ampcorr=ampcorr, ampcorr_meta=acm)


def build_bank(names, pwo4, table_path, wave_cfg, *, column="25320",
               int_tail_ratio=None, L=24, iint_mode="continuous"):
    """TemplateBank for a geometry (names/pwo4 indexed by module index).

    Modules missing from the table or failing the validity rule get the PbWO4
    aggregate of the JSON named in the table header.  Returns (bank, info);
    info has sha1, n_module, fallback (names), column, agg, agg_json, kind.
    A table_path ending in ".npz" (spline_templates_v1) gives a
    SampledTemplateBank instead (column "sampled", agg None; info also has
    n_fallback, run_type, iint_mode and ampcorr (description or None));
    iint_mode applies to the .npz bank only.
    """
    ratio = (wave_cfg["int_tail_ratio"] if int_tail_ratio is None
             else int_tail_ratio)
    if str(table_path).endswith(".npz"):
        with open(table_path, "rb") as fh:
            sha1 = hashlib.sha1(fh.read()).hexdigest()
        bank = SampledTemplateBank.from_npz(
            table_path, names, pwo4, smooth_order=wave_cfg["smooth_order"],
            int_tail_ratio=ratio, tail_break_n=wave_cfg["tail_break_n"], L=L,
            iint_mode=iint_mode)
        info = dict(sha1=sha1, n_module=bank.n_file, fallback=list(bank.fallback),
                    column="sampled", agg=None, agg_json=None,
                    int_tail_ratio=float(ratio),
                    tail_break_n=int(wave_cfg["tail_break_n"]), kind="sampled",
                    n_fallback=bank.n_fallback, run_type=bank.run_type,
                    iint_mode=bank.iint_mode, ampcorr=bank.ampcorr_str())
        return bank, info
    rows, sha1, agg_json = load_h50(table_path, column)
    agg = load_aggregate(agg_json)
    tr = np.full(len(names), agg[0])
    tf = np.full(len(names), agg[1])
    fallback, n_file = [], 0
    for i, (nm, is_p) in enumerate(zip(names, pwo4)):
        if not is_p:
            continue
        if nm in rows:
            tr[i], tf[i] = rows[nm]
            n_file += 1
        else:
            fallback.append(nm)
    bank = TemplateBank(tr, tf, pwo4, smooth_order=wave_cfg["smooth_order"],
                        int_tail_ratio=ratio,
                        tail_break_n=wave_cfg["tail_break_n"], L=L)
    bank.names = list(names)
    bank.fallback = fallback
    info = dict(sha1=sha1, n_module=n_file, fallback=fallback, column=column,
                agg=agg, agg_json=agg_json, int_tail_ratio=float(ratio),
                tail_break_n=int(wave_cfg["tail_break_n"]), kind="two_tau")
    return bank, info


# ---------------------------------------------------------- reconstruction --
class MFResult:
    """Output of MFReconstructor.run (channel-level and pulse-level arrays)."""
    __slots__ = ("ped", "sigma", "clipped", "cand", "npk", "wa_fb",
                 "p_ch", "p_t0", "p_amp", "p_chi2", "p_time", "p_onset",
                 "stats")


class MFReconstructor:
    def __init__(self, bank, *, nsigma=5.0, min_amp=5.0, kmax=4, r_model=0.05,
                 ped_nsamples=30, jmax=None, cand_nsigma=3.0):
        self.b = bank
        self.nsigma, self.min_amp, self.kmax = nsigma, min_amp, int(kmax)
        self.r_model = r_model
        self.ped_n = int(ped_nsamples)
        self.jmax = jmax
        self.cand_ns = cand_nsigma
        self.L = bank.L

    # -- pieces -----------------------------------------------------------
    def _pedestal(self, Y):
        x = Y[:, :self.ped_n]
        med = np.median(x, axis=1)
        mad = np.median(np.abs(x - med[:, None]), axis=1)
        sig = np.maximum(1.4826 * mad, SIGMA_MIN).astype(np.float32)
        band = np.maximum(2.0 * sig, 1.0)[:, None]
        ok = np.abs(x - med[:, None]) <= band
        cnt = np.maximum(ok.sum(axis=1), 1)
        ped = ((x * ok).sum(axis=1) / cnt).astype(np.float32)
        return ped, sig

    def _J(self, ns):
        J = ns - self.L + 1
        if self.jmax is not None:
            J = min(J, self.jmax + 1)
        return J

    def _corr(self, r, Kb):
        J = self._J(r.shape[1])
        Wv = sliding_window_view(r, self.L, axis=1)[:, :J]
        return np.matmul(Wv, Kb.transpose(0, 2, 1))             # [n, J, 4]

    def _basis(self, m, t0, ns):
        """Full-length unit model s(4 i - t0) for peaks (arrays of equal length)."""
        b = self.b
        tau = DT * np.arange(ns, dtype=np.float32)[None, :] - t0[:, None]
        return b.shape(m[:, None], tau).astype(np.float32)

    def _window(self, t0, ns):
        """(sample idx [n, LW], valid [n, LW]) fit window of onset t0."""
        w0 = np.ceil(t0 / DT).astype(np.int64) - LEAD
        idx = w0[:, None] + np.arange(LW)[None, :]
        return np.clip(idx, 0, ns - 1), (idx >= 0) & (idx < ns), idx

    def amp_at(self, y, m, t0, w=None, sigma=None):
        """Fixed-onset LS template amplitude (the K = 1 closed form of _group).

        y [n, ns] pedestal-subtracted, m [n] module indices, t0 [n] onsets in ns
        (same convention as MFResult.p_t0), w [n, ns] optional sample weights.
        Returns (A [n], sigma_A [n] or None); A = <y, s>/<s, s> over the fit
        window of t0, no threshold, no positivity cut; sigma_A = sigma/sqrt(<s, s>).
        """
        b = self.b
        y = np.asarray(y, dtype=np.float32)
        m = np.asarray(m, dtype=np.int64)
        t0 = np.asarray(t0, dtype=np.float32)
        ic, valid, idx = self._window(t0, y.shape[1])
        ww = valid.astype(np.float32)
        if w is not None:
            ww = ww * np.take_along_axis(w, ic, axis=1)
        yy = np.take_along_axis(y, ic, axis=1)
        s = b.shape(m[:, None], DT * idx.astype(np.float32) - t0[:, None])
        Gss = (ww * s * s).sum(1)
        A = (ww * s * yy).sum(1) / np.maximum(Gss, 1e-20)
        sA = None if sigma is None else (
            np.asarray(sigma, dtype=np.float32) / np.sqrt(np.maximum(Gss, 1e-20)))
        return A, sA

    def _gn(self, r, w, m, t0, a, nit=3):
        """3 Gauss-Newton steps on (a, t0) for one peak per row of r."""
        b = self.b
        ns = r.shape[1]
        ic, valid, idx = self._window(t0, ns)
        rr = np.take_along_axis(r, ic, axis=1)
        ww = valid.astype(np.float32)
        if w is not None:
            ww = ww * np.take_along_axis(w, ic, axis=1)
        mc = m[:, None]
        xt = DT * idx.astype(np.float32)
        t0 = t0.astype(np.float32).copy()
        a = a.astype(np.float32).copy()
        for _ in range(nit):
            s, ds = b.shape_and_slope(mc, xt - t0[:, None])
            g = -a[:, None] * ds
            res = rr - a[:, None] * s
            Gss = (ww * s * s).sum(1)
            Gsg = (ww * s * g).sum(1)
            Ggg = (ww * g * g).sum(1)
            bs = (ww * s * res).sum(1)
            bg = (ww * g * res).sum(1)
            det = Gss * Ggg - Gsg * Gsg
            good = det > 1e-6 * np.maximum(Gss * Ggg, 1e-20)
            dsafe = np.where(good, det, 1.0)
            da = np.where(good, (bs * Ggg - bg * Gsg) / dsafe,
                          bs / np.maximum(Gss, 1e-20))
            dt = np.where(good, (Gss * bg - Gsg * bs) / dsafe, 0.0)
            dt = np.clip(dt, -2.0, 2.0)
            a = np.maximum(a + da, 0.2 * a)
            t0 = t0 + dt
        return t0, a

    def _solve_joint(self, y, w, m, t0, a, valid, ns):
        """Joint amplitudes for padded peaks [c, kmax]; drops non-positive ones."""
        c, kmax = t0.shape
        valid = valid.copy()
        wm = np.zeros((c, ns), dtype=np.float32)
        ii = np.arange(ns)[None, :]
        S = np.zeros((c, kmax, ns), dtype=np.float32)
        for k in range(kmax):
            vk = valid[:, k]
            if not vk.any():
                continue
            rows = np.nonzero(vk)[0]
            S[rows, k] = self._basis(m[rows], t0[rows, k], ns)
            w0 = np.ceil(t0[rows, k] / DT).astype(np.int64) - LEAD
            inw = (ii >= w0[:, None]) & (ii < w0[:, None] + LW)
            wm[rows] = np.maximum(wm[rows], inw.astype(np.float32))
        if w is not None:
            wm = wm * w
        eye = np.eye(kmax, dtype=np.float32)[None]
        anew = np.zeros((c, kmax), dtype=np.float32)
        for _ in range(kmax):
            Sv = S * valid[:, :, None]
            Sw = Sv * wm[:, None, :]
            G = np.matmul(Sw, Sv.transpose(0, 2, 1))
            G = G + (~valid)[:, :, None] * eye + 1e-6 * eye * G.max(axis=(1, 2), keepdims=True)
            bvec = np.matmul(Sw, y[:, :, None])[:, :, 0]
            anew = np.linalg.solve(G, bvec[:, :, None])[:, :, 0].astype(np.float32)
            anew = np.where(valid, anew, 0.0)
            neg = valid & (anew <= 0.0)
            if not neg.any():
                break
            worst = np.argmin(np.where(valid, anew, np.inf), axis=1)
            has = neg.any(axis=1)
            valid[np.nonzero(has)[0], worst[has]] = False
        return anew, valid, S, wm

    def _refine_multi(self, y, w, m, t0, a, valid, ns, nit=4):
        """Joint amplitude solve, joint Gauss-Newton on all (a_k, t0_k), solve."""
        c, kmax = t0.shape
        t0 = t0.astype(np.float32).copy()
        a, valid, S, wmask = self._solve_joint(y, w, m, t0, a, valid, ns)
        b = self.b
        mc = m[:, None]
        xt = DT * np.arange(ns, dtype=np.float32)[None, :]
        eye = np.eye(2 * kmax, dtype=np.float32)[None]
        redo = np.ones(c, dtype=bool)        # stage 2 only for channels that dropped
        for stage in range(2):
            upd = valid & redo[:, None]
            for _ in range(nit):
                Sk = np.zeros((c, kmax, ns), dtype=np.float32)
                Gk = np.zeros((c, kmax, ns), dtype=np.float32)
                for k in range(kmax):
                    if not valid[:, k].any():
                        continue
                    s_, ds_ = b.shape_and_slope(mc, xt - t0[:, k, None])
                    vk = valid[:, k, None]
                    Sk[:, k] = s_ * vk
                    Gk[:, k] = -a[:, k, None] * ds_ * vk
                J = np.concatenate([Sk, Gk], axis=1)
                Jw = J * wmask[:, None, :]
                r = (y - np.matmul(a[:, None, :], Sk)[:, 0, :]) * wmask
                N = np.matmul(Jw, J.transpose(0, 2, 1))
                dg = np.diagonal(N, axis1=1, axis2=2)
                inv = np.concatenate([~valid, ~valid], axis=1) | (dg <= 1e-12)
                N = N + eye * (1e-3 * dg[:, :, None] + inv[:, :, None])
                rhs = np.matmul(Jw, r[:, :, None])[:, :, 0]
                sol = np.linalg.solve(N, rhs[:, :, None])[:, :, 0]
                da = np.where(upd, sol[:, :kmax], 0.0)
                dt = np.where(upd, np.clip(sol[:, kmax:], -2.0, 2.0), 0.0)
                a = np.where(valid, np.maximum(a + da, 0.2 * a), 0.0).astype(np.float32)
                t0 = (t0 + dt).astype(np.float32)
            # pulses closer than MIN_SEP_NS are one pulse: drop the smaller, refit
            tsep = np.abs(t0[:, :, None] - t0[:, None, :])
            both = valid[:, :, None] & valid[:, None, :] & ~np.eye(kmax, dtype=bool)[None]
            close = both & (tsep < MIN_SEP_NS)
            if not close.any():
                break
            redo = close.any(axis=(1, 2))
            for ci in np.nonzero(redo)[0]:
                while True:
                    pairs = np.argwhere(close[ci] & valid[ci][:, None] & valid[ci][None, :])
                    if pairs.size == 0:
                        break
                    i, j = pairs[0]
                    valid[ci, i if a[ci, i] < a[ci, j] else j] = False
            a = np.where(valid, a, 0.0).astype(np.float32)
        a, valid, S, wmask = self._solve_joint(y, w, m, t0, a, valid, ns)
        return t0, a, valid, S, wmask

    # -- group reconstruction (all channels share the weights mode) ------------
    def _group(self, y, w, m, sig):
        """y [n,ns] baseline-subtracted (clipped samples zeroed if w), w or None."""
        b = self.b
        n, ns = y.shape
        kmax, L = self.kmax, self.L
        J = self._J(ns)
        res = y.copy()
        pt0 = np.zeros((n, kmax), dtype=np.float32)
        pa = np.zeros((n, kmax), dtype=np.float32)
        pv = np.zeros((n, kmax), dtype=bool)
        active = np.ones(n, dtype=bool)
        for k in range(kmax):
            ia = np.nonzero(active)[0]
            na = ia.size
            if na == 0:
                break
            r = res[ia]
            mk = m[ia]
            Kb = b.K[mk]
            C = self._corr(r, Kb)
            if w is None:
                S = C * b.Kinv[mk][:, None, :]
            else:
                D = self._corr(w[ia], b.K2[mk])
                S = np.where(D >= 0.05 * b.Kn2[mk][:, None, :],
                             C / np.sqrt(np.maximum(D, 1e-12)), -np.inf)
            Sg = np.ascontiguousarray(S[:, :, ::-1]).reshape(na, 4 * J)
            if k > 0:
                # no new pulse within MIN_SEP_NS of an accepted one (both sides)
                tg = np.arange(4 * J, dtype=np.float32)[None, :] - 3.0
                near = np.zeros((na, 4 * J), dtype=bool)
                for kk in range(k):
                    near |= (pv[ia, kk][:, None]
                             & (np.abs(tg - pt0[ia, kk][:, None]) < MIN_SEP_NS))
                Sg = np.where(near, -np.inf, Sg)
            q = np.argmax(Sg, axis=1)
            ar = np.arange(na)
            j = q // 4
            ph = 3 - (q % 4)
            Cb = C[ar, j, ph]
            if w is None:
                Db = b.Kn2[mk, ph]
            else:
                Db = D[ar, j, ph]
            Db = np.maximum(Db, 1e-12)
            A = Cb / Db
            sA = sig[ia] / np.sqrt(Db)
            t0 = (q - 3).astype(np.float32)
            # parabola on the 1 ns score grid
            qm = np.clip(q - 1, 0, 4 * J - 1)
            qp = np.clip(q + 1, 0, 4 * J - 1)
            s0, sm, sp = Sg[ar, q], Sg[ar, qm], Sg[ar, qp]
            with np.errstate(invalid="ignore"):          # -inf neighbours
                den = sm - 2.0 * s0 + sp
                okp = (q > 0) & (q < 4 * J - 1) & (den < 0) & np.isfinite(den)
                delta = np.where(okp, 0.5 * (sm - sp) / np.where(okp, den, 1.0),
                                 0.0)
            t0 = t0 + np.clip(delta, -1.0, 1.0).astype(np.float32)
            # acceptance
            ok = np.isfinite(s0) & (A >= np.maximum(self.nsigma * sA, self.min_amp))
            if k > 0:
                ok &= ~(pv[ia] & (np.abs(t0[:, None] - pt0[ia]) < MIN_SEP_NS)).any(axis=1)
                cov = (pv[ia] & (t0[:, None] >= pt0[ia])
                       & (t0[:, None] < pt0[ia] + DT * L))
                parent = np.where(cov, pa[ia], 0.0).max(axis=1)
                ok &= (A >= self.r_model * parent)
            active[:] = False
            io = np.nonzero(ok)[0]
            if io.size == 0:
                break
            g = ia[io]
            wg = None if w is None else w[g]
            t0r, ar_ = self._gn(res[g], wg, m[g], t0[io], A[io])
            pos = ar_ > 0
            g, io = g[pos], io[pos]
            t0r, ar_ = t0r[pos], ar_[pos]
            if g.size == 0:
                break
            pt0[g, k], pa[g, k], pv[g, k] = t0r, ar_, True
            if k == 0:
                res[g] -= ar_[:, None] * self._basis(m[g], t0r, ns)
                if w is not None:
                    res[g] *= w[g]
            else:
                # back-fit every accepted peak of these channels with the new one
                wg = None if w is None else w[g]
                tt, aa, vv, Sg_, _ = self._refine_multi(
                    y[g], wg, m[g], pt0[g].copy(), pa[g].copy(), pv[g].copy(), ns)
                pt0[g], pa[g], pv[g] = tt, aa, vv
                res[g] = y[g] - np.matmul(aa[:, None, :], Sg_)[:, 0, :]
                if w is not None:
                    res[g] *= w[g]
            active[g] = True

        npk = pv.sum(axis=1)
        # ---- joint amplitude fit -------------------------------------------
        n_joint = int((npk >= 2).sum())
        chi2 = np.zeros(n, dtype=np.float32)
        have = npk > 0
        # K = 1: closed form on the fit window
        one = np.nonzero(npk == 1)[0]
        if one.size:
            k1 = pv[one].argmax(axis=1)        # the surviving slot
            t1 = pt0[one, k1]
            ic, valid, idx = self._window(t1, ns)
            ww = valid.astype(np.float32)
            if w is not None:
                ww = ww * np.take_along_axis(w[one], ic, axis=1)
            yy = np.take_along_axis(y[one], ic, axis=1)
            s = b.shape(m[one][:, None], DT * idx.astype(np.float32) - t1[:, None])
            Gss = (ww * s * s).sum(1)
            a1 = (ww * s * yy).sum(1) / np.maximum(Gss, 1e-20)
            keep = a1 > 0
            pa[one, k1] = np.where(keep, a1, 0.0)
            pv[one[~keep], k1[~keep]] = False
            mdl = a1[:, None] * s
            var = sig[one][:, None] ** 2 + (MODEL_FLOOR * mdl) ** 2
            x2 = (ww * (yy - mdl) ** 2 / var).sum(1)
            dof = np.maximum(ww.sum(1) - 2.0, 1.0)
            chi2[one] = x2 / dof
        # K >= 2: joint solve, refinement pass, joint solve
        mul = np.nonzero(npk >= 2)[0]
        if mul.size:
            ym, mm, sm_ = y[mul], m[mul], sig[mul]
            wm_ = None if w is None else w[mul]
            t0m, am, vm = pt0[mul].copy(), pa[mul].copy(), pv[mul].copy()
            # the back-fit of the last accepted peak already refined (a, t0) of
            # every peak of these channels; re-solve the amplitudes on the final
            # windows (positivity drop) and evaluate chi2
            am, vm, S, wmask = self._solve_joint(ym, wm_, mm, t0m, am, vm, ns)
            mdl = np.matmul(am[:, None, :], S)[:, 0, :]
            var = sm_[:, None] ** 2 + (MODEL_FLOOR * mdl) ** 2
            x2 = (wmask * (ym - mdl) ** 2 / var).sum(1)
            dof = np.maximum(wmask.sum(1) - 2.0 * vm.sum(1), 1.0)
            chi2[mul] = x2 / dof
            pt0[mul], pa[mul], pv[mul] = t0m, am, vm
        return pt0, pa, pv, chi2, n_joint

    # -- public ------------------------------------------------------------
    def run(self, Y, midx, overflow=4095):
        """Reconstruct channels Y[n, ns] (uint16 or float) of module indices midx."""
        b = self.b
        Y = np.asarray(Y)
        n, ns = Y.shape
        midx = np.asarray(midx, dtype=np.int64)
        res = MFResult()
        Yf = Y.astype(np.float32)
        ped, sig = self._pedestal(Yf)
        y = Yf - ped[:, None]
        clipped = (Yf >= overflow).any(axis=1)
        ma = (y[:, :-2] + y[:, 1:-1] + y[:, 2:]) * (1.0 / 3.0)
        cand = ma.max(axis=1) > np.maximum(self.cand_ns * sig, 3.0)
        res.ped, res.sigma, res.clipped, res.cand = ped, sig, clipped, cand
        res.npk = np.zeros(n, dtype=np.int64)
        res.wa_fb = np.zeros(n, dtype=bool)
        pch, pt0_, pamp, pchi, = [], [], [], []
        n_joint = 0
        sat_fb = 0
        for masked in (False, True):
            sel = np.nonzero(cand & (clipped == masked))[0]
            if sel.size == 0:
                continue
            ys = y[sel]
            ms = midx[sel]
            if masked:
                w = (Yf[sel] < overflow).astype(np.float32)
                ys = ys * w
            else:
                w = None
            t0g, ag, vg, chig, nj = self._group(ys, w, ms, sig[sel])
            n_joint += nj
            if masked:
                bad = np.zeros(sel.size, dtype=bool)
                ii = np.arange(ns)[None, :]
                for k in range(self.kmax):
                    rows = np.nonzero(vg[:, k])[0]
                    if rows.size == 0:
                        continue
                    t0k = t0g[rows, k]
                    w0 = np.ceil(t0k / DT).astype(np.int64) - LEAD
                    inw = (ii >= w0[:, None]) & (ii < w0[:, None] + LW)
                    n_ok = (inw * w[rows]).sum(1)
                    tt = DT * ii - t0k[:, None]
                    lead = ((tt > -DT) & (tt <= b.t_pk[ms[rows]][:, None]))
                    n_lead = (lead * w[rows]).sum(1)
                    bad[rows] |= (n_ok < 6) | (n_lead < 2)
                bad |= ~vg.any(axis=1)          # clipped channel with no pulse
                res.wa_fb[sel[bad]] = True
                sat_fb += int(bad.sum())
                vg = vg & ~bad[:, None]
            for k in range(self.kmax):
                rows = np.nonzero(vg[:, k])[0]
                if rows.size == 0:
                    continue
                pch.append(sel[rows])
                pt0_.append(t0g[rows, k])
                pamp.append(ag[rows, k])
                pchi.append(chig[rows])
        if pch:
            p_ch = np.concatenate(pch)
            p_t0 = np.concatenate(pt0_).astype(np.float64)
            p_amp = np.concatenate(pamp).astype(np.float64)
            p_chi2 = np.concatenate(pchi).astype(np.float64)
            o = np.lexsort((p_t0, p_ch))
            p_ch, p_t0, p_amp, p_chi2 = p_ch[o], p_t0[o], p_amp[o], p_chi2[o]
        else:
            p_ch = np.zeros(0, dtype=np.int64)
            p_t0 = p_amp = p_chi2 = np.zeros(0)
        np.add.at(res.npk, p_ch, 1)
        res.p_ch, res.p_t0, res.p_amp, res.p_chi2 = p_ch, p_t0, p_amp, p_chi2
        res.p_time = p_t0 + b.t50[midx[p_ch]]
        res.p_onset = p_t0 / DT
        n_multi = int((res.npk >= 2).sum())
        res.stats = dict(n_cand=int(cand.sum()), n_peaks=int(p_ch.size),
                         n_multi=n_multi, n_joint=int(n_joint),
                         n_sat=int(clipped.sum()), n_sat_fallback=int(sat_fb))
        return res


def pulse_energy_units(bank, midx_of_pulse, t0, amp, kappa=1.0, ped_rms=None):
    """A * I_int(phase) * kappa; with ped_rms (per pulse) the ped floor is applied.

    A bank carrying ampcorr (sampled) also multiplies I_int by f(A).
    """
    if ped_rms is None:
        I = bank.iint(midx_of_pulse, t0)
    else:
        I = bank.iint_floor(midx_of_pulse, t0, amp, ped_rms)
    if getattr(bank, "ampcorr", None) is not None:
        I = I * (bank.f_amp(amp) if bank.ampcorr[1].ndim == 1
                 else bank.f_amp(amp, midx_of_pulse))
    return amp * I * kappa


# ---------------------------------------------------------------- self test --
def _synth(bank, m, pulses, ns=100, ped=200.0, noise=0.0, rng=None, clip=None):
    """pulses: list of (onset ns, amplitude).  Returns float waveform."""
    t = DT * np.arange(ns)
    y = np.full(ns, ped, dtype=np.float64)
    for t0, a in pulses:
        y += a * bank.shape(m, t - t0, dtype=np.float64)
    if noise > 0:
        y += rng.normal(0.0, noise, ns)
    if clip is not None:
        y = np.minimum(y, clip)
    return y


def self_test(verbose=True):
    rng = np.random.default_rng(12345)
    wcfg = dict(smooth_order=2, int_tail_ratio=0.1, tail_break_n=2)
    nm = 3
    taus = [(2.499, 24.41), (2.0, 20.0), (3.2, 28.0)]
    bank = TemplateBank([t[0] for t in taus], [t[1] for t in taus], [True] * nm,
                        smooth_order=2, int_tail_ratio=0.1, tail_break_n=2)
    rec = MFReconstructor(bank)
    ok_all = True
    print(f"template peak check: max p = "
          f"{template(np.linspace(0, 100, 100001), 2.499, 24.41, tmax_of(2.499, 24.41)).max():.6f}")

    def report(tag, passed, txt):
        nonlocal ok_all
        ok_all &= bool(passed)
        print(f"[{'PASS' if passed else 'FAIL'}] {tag}: {txt}")

    # ---- single pulses
    N = 600
    mi = rng.integers(0, nm, N)
    amp = np.exp(rng.uniform(np.log(20), np.log(3000), N))
    t0 = rng.uniform(80.0, 160.0, N)
    Y = np.stack([_synth(bank, mi[i], [(t0[i], amp[i])]) for i in range(N)]
                 ).astype(np.float32)
    r = rec.run(Y, mi)
    found = r.npk == 1
    dA = np.full(N, np.nan)
    dt = np.full(N, np.nan)
    for c, a_, t_ in zip(r.p_ch, r.p_amp, r.p_t0):
        dA[c] = abs(a_ - amp[c]) / amp[c]
        dt[c] = abs(t_ - t0[c])
    report("single", found.all() and np.nanmax(dA) < 2e-3 and np.nanmax(dt) < 0.05,
           f"N={N} found={int(found.sum())} max|dA|/A={np.nanmax(dA)*100:.4f}% "
           f"max|dt|={np.nanmax(dt):.4f} ns")
    tt = r.p_time - (t0[r.p_ch] + bank.t50[mi[r.p_ch]])
    print(f"        reported-time check (t - (t0+t50)): max|.|={np.abs(tt).max():.4f} ns")

    # ---- pairs
    cells = []
    for dT in (16, 20, 24, 32, 48, 64):
        for ratio in (0.2, 0.5, 1.0, 2.0):
            cells.append((dT, ratio))
    worst = 0.0
    allfound = True
    lines = []
    for dT, ratio in cells:
        n = 40
        mi = rng.integers(0, nm, n)
        a1 = np.exp(rng.uniform(np.log(200), np.log(3500.0 / (1.0 + ratio)), n))
        a2 = a1 * ratio
        ta = rng.uniform(90.0, 130.0, n)
        tb = ta + dT
        Y = np.stack([_synth(bank, mi[i], [(ta[i], a1[i]), (tb[i], a2[i])])
                      for i in range(n)]).astype(np.float32)
        r = rec.run(Y, mi)
        nf = 0
        eA = []
        eT = []
        for c in range(n):
            sel = np.nonzero(r.p_ch == c)[0]
            tc, ac = r.p_t0[sel], r.p_amp[sel]
            if sel.size < 2:
                continue
            ia = np.argmin(np.abs(tc - ta[c]))
            ib = np.argmin(np.abs(tc - tb[c]))
            if ia == ib:
                continue
            if abs(tc[ia] - ta[c]) < 3 and abs(tc[ib] - tb[c]) < 3:
                nf += 1
                eA.append(max(abs(ac[ia] - a1[c]) / a1[c], abs(ac[ib] - a2[c]) / a2[c]))
                eT.append(max(abs(tc[ia] - ta[c]), abs(tc[ib] - tb[c])))
        mx = max(eA) if eA else float('nan')
        lines.append(f"dT={dT:>2} ratio={ratio:<3} both found {nf}/{n} "
                     f"max|dA|/A={mx*100:6.3f}% max|dt|={max(eT) if eT else float('nan'):.3f} ns")
        if dT >= 24:
            allfound &= (nf == n)
            worst = max(worst, mx if eA else 1.0)
    report("pairs (dT>=24)", allfound and worst < 0.01,
           f"all found={allfound} worst max|dA|/A={worst*100:.3f}%")
    if verbose:
        for l in lines:
            print("        " + l)

    # ---- clipped
    n = 200
    mi = rng.integers(0, nm, n)
    amp = np.exp(rng.uniform(np.log(300), np.log(3000), n))
    t0 = rng.uniform(80.0, 160.0, n)
    Y = np.empty((n, 100), dtype=np.float32)
    clv = np.empty(n)
    for i in range(n):
        pk = 200.0 + amp[i]
        clv[i] = 200.0 + 0.6 * amp[i]
        Y[i] = _synth(bank, mi[i], [(t0[i], amp[i])], clip=clv[i])
    errs, nfb, nfound = [], 0, 0
    for i in range(n):
        r = rec.run(Y[i:i + 1], mi[i:i + 1], overflow=clv[i] - 1e-3)
        if r.wa_fb[0]:
            nfb += 1
            continue
        if r.p_ch.size == 1:
            nfound += 1
            errs.append(abs(r.p_amp[0] - amp[i]) / amp[i])
    report("clipped 0.6x", nfound + nfb == n and max(errs) < 0.01,
           f"N={n} fit={nfound} wa_fallback={nfb} max|dA|/A={max(errs)*100:.4f}% "
           f"median={np.median(errs)*100:.4f}%")

    # ---- fixed-onset amplitude amp_at (P3-1)
    n = 400
    mi = rng.integers(0, nm, n)
    amp = np.exp(rng.uniform(np.log(20), np.log(3000), n))
    t0 = rng.uniform(80.0, 160.0, n)
    Y = np.stack([_synth(bank, mi[i], [(t0[i], amp[i])]) for i in range(n)]
                 ).astype(np.float32)
    r = rec.run(Y, mi)
    one = np.nonzero(r.npk == 1)[0]
    sel = np.searchsorted(r.p_ch, one)
    yp = Y - r.ped[:, None]
    A, sA = rec.amp_at(yp[one], mi[one], r.p_t0[sel].astype(np.float32),
                       sigma=r.sigma[one])
    dr = np.max(np.abs(A - r.p_amp[sel]) / r.p_amp[sel])
    report("amp_at == K=1 closed form", one.size == n and dr <= 1e-5,
           f"N={one.size} max rel diff={dr:.2e} (<=1e-5) sigma_A[0]={sA[0]:.4f}")

    n = 400
    mi = rng.integers(0, nm, n)
    amp = np.exp(rng.uniform(np.log(5), np.log(3000), n))
    t0 = rng.uniform(80.0, 160.0, n)
    Y = np.stack([_synth(bank, mi[i], [(t0[i], amp[i])], ped=0.0)
                  for i in range(n)]).astype(np.float32)
    A, _ = rec.amp_at(Y, mi, t0)
    e0 = np.max(np.abs(A / amp - 1.0))
    report("amp_at dt=0", e0 <= 2e-3,
           f"N={n} A 5-3000 ADC max|A/A_true-1|={e0*100:.5f}% (<=0.2%)")

    # loss curve at a wrong fixed onset, unit template of the first bank module
    f7 = {0.5: 0.2, 1.0: 1.0, 1.5: 2.6, 2.0: 5.0, 3.0: 8.0}
    n = 2000
    t0 = rng.uniform(80.0, 160.0, n)
    Y = np.stack([_synth(bank, 0, [(t0[i], 1.0)], ped=0.0)
                  for i in range(n)]).astype(np.float32)
    m0 = np.zeros(n, dtype=np.int64)
    print("        F7 loss curve (module 0, 2.499/24.41 ns), mean A/A_true-1 over "
          f"{n} random phases:")
    for d, f in f7.items():
        ap, an = rec.amp_at(Y, m0, t0 + d)[0] - 1.0, rec.amp_at(Y, m0, t0 - d)[0] - 1.0
        print(f"        err={d:+.1f} ns: mean {ap.mean()*100:+7.3f}% worst {ap.min()*100:+7.3f}%"
              f" | err={-d:+.1f} ns: mean {an.mean()*100:+7.3f}% worst {an.min()*100:+7.3f}%"
              f" | F7 loss {f:.1f}%")
    print("[info] F7 comparison is informational (not a gate): see plan_20261023 F7")

    # ---- noisy sanity (not a gate)
    n = 400
    mi = rng.integers(0, nm, n)
    amp = np.exp(rng.uniform(np.log(30), np.log(1500), n))
    t0 = rng.uniform(90.0, 150.0, n)
    Y = np.stack([_synth(bank, mi[i], [(t0[i], amp[i])], noise=1.5, rng=rng)
                  for i in range(n)]).astype(np.float32)
    r = rec.run(Y, mi)
    one = np.nonzero(r.npk == 1)[0]
    dA = np.zeros(n) + np.nan
    dT_ = np.zeros(n) + np.nan
    for c, a_, t_ in zip(r.p_ch, r.p_amp, r.p_t0):
        dA[c], dT_[c] = (a_ - amp[c]), (t_ - t0[c])
    big = amp > 150
    print(f"[info] noise 1.5 ADC: found 1 pulse in {one.size}/{n}; extra pulses "
          f"{int((r.npk > 1).sum())}; A>150: sigma(dA)/A="
          f"{np.nanstd(dA[big]/amp[big])*100:.2f}% sigma(dt)={np.nanstd(dT_[big]):.3f} ns; "
          f"A<60: miss={int(((amp<60)&(r.npk==0)).sum())}/{int((amp<60).sum())}")

    # ---- I_int regression of the aggregate two-tau (recorded from the pre-sampled code)
    agg_bank = TemplateBank([AGG_FALLBACK[0]], [AGG_FALLBACK[1]], [True], **wcfg)
    ref_mean, ref_00 = 7.212860678845598, 7.027331829071045
    e_mean = abs(agg_bank.I_mean[0] / ref_mean - 1.0)
    e_00 = abs(float(agg_bank.I_tab[0, 0]) / ref_00 - 1.0)
    report("I_int regression (aggregate)", e_mean < 1e-9 and e_00 < 1e-9,
           f"I_mean={agg_bank.I_mean[0]:.12f} (ref {ref_mean:.12f}, rel {e_mean:.1e}) "
           f"I_tab[0,0]={float(agg_bank.I_tab[0, 0]):.9f} (ref {ref_00:.9f}, rel {e_00:.1e})")

    # ---- sampled bank built from the analytic shapes
    import tempfile
    from scipy.interpolate import make_interp_spline
    rng2 = np.random.default_rng(777)
    tg = np.arange(0.0, 400.0 + 1e-9, 0.1)       # T_end 400: analytic tail kept
    spl = [make_interp_spline(tg, template(tg, tr_, tf_, tmax_of(tr_, tf_)), k=3)
           for tr_, tf_ in taus]
    sb = SampledTemplateBank(spl, [f"M{i}" for i in range(nm)], [True] * nm,
                             smooth_order=2, int_tail_ratio=0.1, tail_break_n=2,
                             T_end=400.0)
    dK = np.abs(sb.K - bank.K).max()
    dt50 = np.abs(sb.t50 - bank.t50).max()
    rI = np.abs(sb.I_tab - bank.I_tab) / bank.I_tab
    nflip = int((rI > 1e-4).sum())             # walk-out cut flips on a sample at the cut
    dIm = np.abs(sb.I_mean / bank.I_mean - 1.0).max()
    report("sampled == analytic (K, t50, I_tab)",
           dK < 1e-4 and dt50 < 1e-3 and np.median(rI) < 1e-4 and nflip <= 3
           and dIm < 1e-3,
           f"max|dK|={dK:.2e} max|dt50|={dt50:.2e} ns (<=1e-3) median rel|dI_tab|="
           f"{np.median(rI):.1e} (<=1e-4), {nflip}/{rI.size} phases differ >1e-4 "
           f"(cut flips, <=3), max rel dI_mean={dIm:.1e} (<=1e-3)")

    # ---- continuous walk-out integral vs the phase average of the discrete table
    rc = np.abs(sb.I_cont_walk / sb.I_mean - 1.0)
    trv, tfv, tmv = np.array([t[0] for t in taus]), np.array([t[1] for t in taus]), bank.tmx
    I_an = continuous_walk_integral(
        lambda tau: template(tau, trv[:, None], tfv[:, None], tmv[:, None]), 2, 0.1, 2,
        nrow=nm)
    rc_an = np.abs(I_an / bank.I_mean - 1.0)
    pk2 = bank.iint(np.zeros(5, dtype=int), np.linspace(100, 104, 5))   # two-tau still phase-dependent
    sb_ph = SampledTemplateBank(spl, [f"M{i}" for i in range(nm)], [True] * nm,
                                smooth_order=2, int_tail_ratio=0.1, tail_break_n=2,
                                T_end=400.0, iint_mode="phase")
    t0s = np.linspace(100.0, 104.0, 64, endpoint=False)
    Ic_ = sb.iint(np.full(64, 1), t0s)
    Ip_ = sb_ph.iint(np.full(64, 1), t0s)
    flc = sb.iint_floor(np.full(64, 1), t0s, np.full(64, 400.0), np.full(64, 2.0))
    fl_hi = 0.2                                   # floor above ratio * peak: the floor path
    flc2 = sb.iint_floor(np.full(64, 2), t0s, np.full(64, 5.0), np.full(64, 5.0 * fl_hi))
    flp2 = sb_ph.iint_floor(np.full(64, 2), t0s, np.full(64, 5.0), np.full(64, 5.0 * fl_hi))
    rf = abs(flc2.mean() / flp2.mean() - 1.0)
    report("continuous I_int vs phase average (I_mean)",
           rc.max() < 2e-3 and rc_an.max() < 2e-3 and np.ptp(pk2) > 0.01
           and sb.iint_mode == "continuous" and sb_ph.iint_mode == "phase"
           and np.all(Ic_ == sb.I_cont[1]) and Ic_.shape == t0s.shape
           and np.ptp(Ip_) > 0.01 and np.array_equal(sb_ph.I_tab, sb.I_tab)
           and np.all(flc == sb.I_cont[1]) and np.ptp(flc2) == 0.0 and rf < 5e-3
           and flc2[0] < sb.I_cont[2],
           f"max rel|I_cont/I_mean-1| = {rc.max():.2e} (sampled), {rc_an.max():.2e} (analytic) "
           f"(<=2e-3); phase-table p-p {np.ptp(Ip_) / Ip_.mean() * 100:.1f}% vs continuous 0; "
           f"floored (floor {fl_hi}) cont/phase-mean-1 = {rf:.1e} (<=5e-3), "
           f"I {flc2[0]:.3f} < unfloored {sb.I_cont[2]:.3f}")

    def single_gate(tag, sbank):
        n = 600
        mi_ = rng2.integers(0, nm, n)
        amp_ = np.exp(rng2.uniform(np.log(50), np.log(3000), n))   # A>=50: the analytic gate itself reaches ~0.2-0.3 % near A~24 (onset-kink GN difference)
        t0_ = rng2.uniform(80.0, 160.0, n)
        Y_ = np.stack([_synth(bank, mi_[i], [(t0_[i], amp_[i])])
                       for i in range(n)]).astype(np.float32)
        r_ = MFReconstructor(sbank).run(Y_, mi_)
        dA_, dt_ = np.full(n, np.nan), np.full(n, np.nan)
        for c, a_, t_ in zip(r_.p_ch, r_.p_amp, r_.p_t0):
            dA_[c], dt_[c] = abs(a_ - amp_[c]) / amp_[c], abs(t_ - t0_[c])
        report(tag, (r_.npk == 1).all() and np.nanmax(dA_) < 2e-3 and np.nanmax(dt_) < 0.05,
               f"N={n} found={int((r_.npk == 1).sum())} max|dA|/A={np.nanmax(dA_)*100:.4f}% "
               f"max|dt|={np.nanmax(dt_):.4f} ns (analytic waveforms)")

    single_gate("sampled single", sb)

    # ---- round trip through the spline_templates_v1 .npz and build_bank
    def pad(rows_):
        w = max(len(r_) for r_ in rows_)
        return np.array([np.r_[r_, np.full(w - len(r_), np.nan)] for r_ in rows_])

    with tempfile.TemporaryDirectory() as td:
        path = os.path.join(td, "spline_test.npz")
        np.savez(path, names=np.array(["A", "B", "C"]), valid=np.ones(3, dtype=int),
                 fallback=np.zeros(3, dtype=int), T_end=400.0, spline_k=3,
                 spline_t=pad([x.t for x in spl]), spline_c=pad([x.c for x in spl]),
                 spline_nt=np.array([len(x.t) for x in spl]),
                 spline_nc=np.array([len(x.c) for x in spl]),
                 avg_spline_t=spl[0].t, avg_spline_c=spl[0].c,
                 meta=json.dumps(dict(run_type="selftest")))
        gn, gp = ["A", "B", "C", "D", "E"], [True, True, True, True, False]
        fb_bank, info = build_bank(gn, gp, path, wcfg)
        ph_bank, _ = build_bank(gn, gp, path, wcfg, iint_mode="phase")
        sel_ = [0, 1, 2]
        dK2 = np.abs(fb_bank.K[sel_] - bank.K).max()
        dI2 = np.median(np.abs(ph_bank.I_tab[sel_] - bank.I_tab) / bank.I_tab)
        dav = np.abs(fb_bank.K[3] - bank.K[0]).max()
        okb = (info["column"] == "sampled" and info["kind"] == "sampled"
               and info["fallback"] == ["D"] and info["n_module"] == 3
               and info["n_fallback"] == 1 and info["run_type"] == "selftest"
               and info["agg"] is None and info["agg_json"] is None
               and info["iint_mode"] == "continuous" and info["ampcorr"] is None
               and fb_bank.ampcorr is None
               and len(info["sha1"]) == 40 and fb_bank.fallback == ["D"]
               and resolve_column(path) == "sampled")
        report("build_bank(.npz)", okb and dK2 < 1e-4 and dI2 < 1e-4 and dav < 1e-6
               and np.isfinite(fb_bank.I_mean).all(),
               f"max|dK|={dK2:.2e} median rel|dI_tab|={dI2:.2e} fallback D==avg: {dav:.1e} "
               f"info={ {k_: info[k_] for k_ in ('n_module', 'n_fallback', 'run_type', 'column')} }")
        single_gate("sampled single (npz round trip)", fb_bank)
        pe = pulse_energy_units(ph_bank, np.array([0, 1]), np.array([100.0, 101.3]),
                                np.array([500.0, 800.0]), ped_rms=np.array([1.0, 2.0]))
        pe0 = pulse_energy_units(bank, np.array([0, 1]), np.array([100.0, 101.3]),
                                 np.array([500.0, 800.0]), ped_rms=np.array([1.0, 2.0]))
        report("sampled pulse_energy_units (ped floor)", np.abs(pe / pe0 - 1).max() < 1e-3,
               f"max rel diff to analytic = {np.abs(pe / pe0 - 1).max():.2e}")
        pc = pulse_energy_units(fb_bank, np.array([0, 1]), np.array([100.0, 101.3]),
                                np.array([500.0, 800.0]), ped_rms=np.array([1.0, 2.0]))
        pc2 = pulse_energy_units(fb_bank, np.array([0, 1]), np.array([100.0, 101.3]),
                                 np.array([500.0, 800.0]))
        pm_ = np.array([500.0, 800.0]) * bank.I_mean[:2]
        report("sampled pulse_energy_units (continuous)",
               np.array_equal(pc, pc2) and np.abs(pc / pm_ - 1).max() < 2e-3,
               f"max rel diff to A*I_mean = {np.abs(pc / pm_ - 1).max():.2e}")
        pb = pulse_energy_units(bank, np.array([0, 1]), np.array([100.0, 101.3]),
                                np.array([500.0, 800.0]), 1.01)
        report("two-tau pulse_energy_units unchanged (no ampcorr)",
               np.array_equal(pb, np.array([500.0, 800.0])
                              * bank.iint(np.array([0, 1]), np.array([100.0, 101.3])) * 1.01)
               and not hasattr(bank, "ampcorr"), "A * I_int * kappa exactly")

        # ---- amplitude correction f(A) from the .npz keys
        kn_la, kn_f = np.log([200.0, 1000.0]), np.array([1.01, 0.99])
        path2 = os.path.join(td, "spline_test_ampcorr.npz")
        np.savez(path2, names=np.array(["A", "B", "C"]), valid=np.ones(3, dtype=int),
                 fallback=np.zeros(3, dtype=int), T_end=400.0, spline_k=3,
                 spline_t=pad([x.t for x in spl]), spline_c=pad([x.c for x in spl]),
                 spline_nt=np.array([len(x.t) for x in spl]),
                 spline_nc=np.array([len(x.c) for x in spl]),
                 avg_spline_t=spl[0].t, avg_spline_c=spl[0].c,
                 meta=json.dumps(dict(run_type="selftest")),
                 ampcorr_logA=kn_la, ampcorr_f=kn_f,
                 ampcorr_meta=json.dumps(dict(scope="global", note="selftest")))
        ac_bank, ac_info = build_bank(gn, gp, path2, wcfg)
        Aa = np.array([100.0, 200.0, 447.0, 1000.0, 3000.0])
        ma = np.array([0, 1, 2, 0, 1])
        ta = np.linspace(100, 103, 5)
        got = pulse_energy_units(ac_bank, ma, ta, Aa, 1.01) / pulse_energy_units(fb_bank, ma, ta, Aa, 1.01)
        got_p = (pulse_energy_units(ac_bank, ma, ta, Aa, ped_rms=np.full(5, 1.0))
                 / pulse_energy_units(fb_bank, ma, ta, Aa, ped_rms=np.full(5, 1.0)))
        fx = np.interp(np.log(np.maximum(Aa, 1.0)), kn_la, kn_f)
        f447 = 1.01 + (0.99 - 1.01) * (np.log(447.0) - kn_la[0]) / (kn_la[1] - kn_la[0])
        report("amplitude correction f(A)",
               np.abs(got / fx - 1).max() < 1e-12 and np.abs(got_p / fx - 1).max() < 1e-12
               and abs(fx[2] - f447) < 1e-12 and fx[0] == 1.01 and fx[1] == 1.01
               and fx[3] == 0.99 and fx[4] == 0.99
               and np.array_equal(ac_bank.ampcorr[0], kn_la)
               and np.array_equal(ac_bank.ampcorr[1], kn_f)
               and ac_bank.ampcorr_meta == dict(scope="global", note="selftest")
               and isinstance(ac_info["ampcorr"], str) and ac_info["iint_mode"] == "continuous"
               and np.array_equal(ac_bank.f_amp(Aa), fx)
               and np.array_equal(fb_bank.f_amp(Aa), np.ones(5)),
               f"f at A={Aa.tolist()}: {np.round(got, 6).tolist()} (exact interp, "
               f"max dev {np.abs(got / fx - 1).max():.1e}); info ampcorr: {ac_info['ampcorr']}")

        # ---- the generalised band layout of spline_template_extract `build --bands`
        # (band_edges / shape_band / r_band / n_band / A_med_band) is ignored by the bank
        path3 = os.path.join(td, "spline_test_bands.npz")
        np.savez(path3, **dict(np.load(path2), band_edges=np.array([100.0, 250.0, 500.0]),
                               shape_band=np.full((2, 3, 5), np.nan, dtype=np.float32),
                               r_band=np.ones((2, 3)), n_band=np.zeros((2, 3), dtype=np.int32),
                               A_med_band=np.full((2, 3), np.nan)))
        bd_bank, bd_info = build_bank(gn, gp, path3, wcfg)
        report("band-layout npz (band_edges, shape_band, ...) loads like the three-band one",
               np.array_equal(bd_bank.ampcorr[0], kn_la) and np.array_equal(bd_bank.ampcorr[1], kn_f)
               and np.array_equal(bd_bank.f_amp(Aa), fx)
               and np.array_equal(pulse_energy_units(bd_bank, ma, ta, Aa, 1.01),
                                  pulse_energy_units(ac_bank, ma, ta, Aa, 1.01)),
               f"info ampcorr: {bd_info['ampcorr']}")

        # ---- per-module f(A) table [file rows, grid] (`ampcorr --per-module-f`)
        pm_la = np.log([200.0, 400.0, 1000.0])
        pm_f = np.interp(pm_la, kn_la, kn_f)
        path4 = os.path.join(td, "spline_test_pm.npz")
        np.savez(path4, **dict(np.load(path2), ampcorr_logA=pm_la,
                               ampcorr_f=np.tile(pm_f, (3, 1)),
                               ampcorr_meta=json.dumps(dict(scope="per_run_type",
                                                            per_module=True))))
        pm_bank, pm_info = build_bank(gn, gp, path4, wcfg)
        Ab = np.array([100.0, 200.0, 300.0, 447.0, 1000.0, 3000.0, 600.0])
        mb = np.array([0, 1, 2, 3, 0, 1, 2])
        tb = np.linspace(100, 103, Ab.size)
        pm_e = pulse_energy_units(pm_bank, mb, tb, Ab, 1.01)
        d_e = pulse_energy_units(ac_bank, mb, tb, Ab, 1.01)
        pm_ep = pulse_energy_units(pm_bank, mb, tb, Ab, ped_rms=np.full(Ab.size, 1.0))
        d_ep = pulse_energy_units(ac_bank, mb, tb, Ab, ped_rms=np.full(Ab.size, 1.0))
        # rows differ: module B x 1.02, file default row (D) x 0.97, A and C unchanged
        pm_f2 = np.tile(pm_f, (3, 1))
        pm_f2[1] *= 1.02
        path5 = os.path.join(td, "spline_test_pm2.npz")
        np.savez(path5, **dict(np.load(path4), ampcorr_f=pm_f2, ampcorr_f_rt=pm_f * 0.97))
        pm2_bank, _ = build_bank(gn, gp, path5, wcfg)
        pm2_e = pulse_energy_units(pm2_bank, mb, tb, Ab, 1.01)
        want = d_e * np.array([1.0, 1.02, 1.0, 0.97, 1.0, 1.02, 1.0])
        report("per-module f(A) table",
               np.array_equal(pm_e, d_e) and np.array_equal(pm_ep, d_ep)
               and pm_bank.ampcorr[1].shape == (5, 3)
               and isinstance(pm_info["ampcorr"], str) and "per-module" in pm_info["ampcorr"]
               and np.abs(pm2_e / want - 1).max() < 1e-12
               and np.array_equal(pm_bank.f_amp(Ab, mb),
                                  np.interp(np.log(np.maximum(Ab, 1.0)), pm_la, pm_f)),
               f"equal rows == 1-D table: max |dE/E| {np.abs(pm_e / d_e - 1).max():.1e} "
               f"(ped floor {np.abs(pm_ep / d_ep - 1).max():.1e}); distinct rows "
               f"max dev {np.abs(pm2_e / want - 1).max():.1e}; info ampcorr: {pm_info['ampcorr']}")

    print("SELF-TEST", "PASS" if ok_all else "FAIL")
    return ok_all


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        sys.exit(0 if self_test() else 1)
    print(__doc__)
