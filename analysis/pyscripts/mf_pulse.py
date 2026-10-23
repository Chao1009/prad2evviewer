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

Energy: E_units = A * I_int[m](phase) * kappa, where I_int is the WaveAnalyzer
walk-out integral (WaveAnalyzer.cpp:391-462, smoothing :235) of the noiseless
sample-aligned unit template, in the high-amplitude limit (cut = int_tail_ratio
* smoothed peak, tail_break_n consecutive sub-cut samples end the walk).

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
    """`column` as given, or `prefer` if the table has it, else its last column."""
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

    def __init__(self, tau_r, tau_f, pwo4, *, smooth_order, int_tail_ratio,
                 tail_break_n, L=24):
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


def build_bank(names, pwo4, table_path, wave_cfg, *, column="25320",
               int_tail_ratio=None, L=24):
    """TemplateBank for a geometry (names/pwo4 indexed by module index).

    Modules missing from the table or failing the validity rule get the PbWO4
    aggregate of the JSON named in the table header.  Returns (bank, info);
    info has sha1, n_module, fallback (names), column, agg, agg_json.
    """
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
    ratio = (wave_cfg["int_tail_ratio"] if int_tail_ratio is None
             else int_tail_ratio)
    bank = TemplateBank(tr, tf, pwo4, smooth_order=wave_cfg["smooth_order"],
                        int_tail_ratio=ratio,
                        tail_break_n=wave_cfg["tail_break_n"], L=L)
    info = dict(sha1=sha1, n_module=n_file, fallback=fallback, column=column,
                agg=agg, agg_json=agg_json, int_tail_ratio=float(ratio),
                tail_break_n=int(wave_cfg["tail_break_n"]))
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
        return template(tau, b.tau_r[m][:, None].astype(np.float32),
                        b.tau_f[m][:, None].astype(np.float32),
                        b.tmx[m][:, None].astype(np.float32)).astype(np.float32)

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
        s = template(DT * idx.astype(np.float32) - t0[:, None],
                     b.tau_r[m][:, None].astype(np.float32),
                     b.tau_f[m][:, None].astype(np.float32),
                     b.tmx[m][:, None].astype(np.float32))
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
        tr = b.tau_r[m][:, None].astype(np.float32)
        tf = b.tau_f[m][:, None].astype(np.float32)
        tm = b.tmx[m][:, None].astype(np.float32)
        xt = DT * idx.astype(np.float32)
        t0 = t0.astype(np.float32).copy()
        a = a.astype(np.float32).copy()
        for _ in range(nit):
            s, ds = template_and_slope(xt - t0[:, None], tr, tf, tm)
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
        tr = b.tau_r[m][:, None].astype(np.float32)
        tf = b.tau_f[m][:, None].astype(np.float32)
        tm = b.tmx[m][:, None].astype(np.float32)
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
                    s_, ds_ = template_and_slope(xt - t0[:, k, None], tr, tf, tm)
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
            s = template(DT * idx.astype(np.float32) - t1[:, None],
                         b.tau_r[m[one]][:, None].astype(np.float32),
                         b.tau_f[m[one]][:, None].astype(np.float32),
                         b.tmx[m[one]][:, None].astype(np.float32))
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
    """A * I_int(phase) * kappa; with ped_rms (per pulse) the ped floor is applied."""
    if ped_rms is None:
        I = bank.iint(midx_of_pulse, t0)
    else:
        I = bank.iint_floor(midx_of_pulse, t0, amp, ped_rms)
    return amp * I * kappa


# ---------------------------------------------------------------- self test --
def _synth(bank, m, pulses, ns=100, ped=200.0, noise=0.0, rng=None, clip=None):
    """pulses: list of (onset ns, amplitude).  Returns float waveform."""
    t = DT * np.arange(ns)
    y = np.full(ns, ped, dtype=np.float64)
    for t0, a in pulses:
        y += a * template(t - t0, bank.tau_r[m], bank.tau_f[m], bank.tmx[m])
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
    print("SELF-TEST", "PASS" if ok_all else "FAIL")
    return ok_all


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        sys.exit(0 if self_test() else 1)
    print(__doc__)
