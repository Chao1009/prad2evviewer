#!/usr/bin/env python3
"""mf_ampshape_diag.py -- option (ii) offline feasibility diagnostic (plan_20261025 section 4.4).

Does an amplitude-indexed template remove the residual amplitude dependence of
Q_template / Q_rawsum that the run-type factor f(A) of option (i) leaves?  On the stored
waveforms of the study-pulse dumps (<dump>/pulses.npz: wave, wave_row, wave_ped, wave_sigma)
every selected pulse is refit twice:

  (a) module template + f(A)           E = A * I_int(module shape) * f(A)         (the bank as run)
  (b) shape interpolated linearly in ln A between the two neighbouring band templates of the
      module (shape_band of the bank), E = A * I_int(interpolated shape), NO f(A).

Pulses: stored-wave pulses (single mf pulse, unclipped, A >= 50) of PbWO4 modules that also have
one unclipped single-peak (npk == 1) rawsum row of the same (event, midx) with an onset inside
[3, ns-4] and a fit window inside the waveform; at most --max-per-module per module and run
(seeded random choice over the run's splits).  Q_raw = rawsum E_units * 4 ns, Q = E_units / kappa
* 4 ns, as in mf_fit_study intraw.

Fit (both arms, one pulse at a time, vectorised): pedestal and sigma as the mf (first
ped_nsamples samples, 2-sigma clipped mean; sigma = the dump's wave_sigma); fit window of the
mf, LW = 32 samples from ceil(t0/4) - 3 with t0 the dump's onset, held fixed; unweighted LSQ
(the mf's K = 1 closed form) for amplitude and onset by Gauss-Newton on (A, t0), 6 steps from the
dump's (A, t0), the amplitude finally in closed form at the final t0.  chi2/dof =
sum (y - A s)^2 / (sigma^2 + (0.03 A s)^2) over the window / (n_window - 2), as the mf.  In (b)
the shape depends on the fitted A: 3 passes, the interpolation weight recomputed from the
previous pass's A.  time = t0 + t50 of the arm's own shape (50 % leading-edge point).

Band shapes: knot of band b = median over the modules with a usable band-b shape (finite shape
and A_med_band, band present for >= MIN_BAND_MODULES modules) of A_med_band[b]; interpolation
is linear in ln A between the module's two neighbouring knots, the outer band outside the
outer knots; a module lacking a band uses its module template for that side (so arm (b) there
degrades towards the module template without f(A)).  I_int of any shape is the bank's energy
integral: the phase average of the WaveAnalyzer walk-out over the continuous template
(mf_pulse.SampledTemplateBank.I_cont), evaluated on tables of the interpolation weight (WGRID
steps per knot interval) and read by linear interpolation in the weight.  "cov" rows restrict to
pulses whose two interpolation sides are both true band shapes.

Output: <out>/ampshape_diag.tsv, <out>/ampshape_diag.md, <out>/ampshape_ratio_<run>.png.

Usage (MPLBACKEND=Agg; needs setup_env.sh for prad2py):
  python3 mf_ampshape_diag.py \\
      --dumps "output/mf_fit_study_20261025/dump/fa/*/pulses.npz" \\
      --templates "output/templates_spline_20261025/{runtype}/templates_{runtype}.npz" \\
      --out output/templates_spline_20261025/ampshape_diag/
"""
from __future__ import annotations

import argparse
import glob
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import mf_pulse as M
from mf_fit_study import bin_med_se, flatness

NAN = float("nan")
EDGES = np.array([100, 250, 500, 1000, 2000, 4000], float)
MIN_BAND_MODULES = 10          # a band is a knot only if this many modules have a usable shape
WGRID = 8                      # I_int tables: interpolation weight steps per knot interval
NIT = 6                        # Gauss-Newton steps per pass
NPASS_B = 3                    # arm (b): passes of the A-dependent shape
CHUNK = 20000
DT, LW, LEAD = M.DT, M.LW, M.LEAD
LUT_DT = M._LUT_DT
MODEL_FLOOR = M.MODEL_FLOOR


# ------------------------------------------------------------------ selection --
def load_run(paths, max_per_module, seed):
    """Selected pulses of one run (dump splits concatenated).  Returns dict of arrays."""
    parts, names, kappa, run = [], None, None, None
    for ip, p in enumerate(sorted(paths)):
        z = np.load(p)
        nm = [str(x) for x in z["module_names"]]
        if names is None:
            names, kappa, run = nm, float(z["kappa"]), int(z["run"])
        elif nm != names or int(z["run"]) != run:
            sys.exit(f"FATAL: {p}: module list or run differs between the dumps of one run")
        ev, mi = z["event"].astype(np.int64), z["midx"].astype(np.int64)
        key = (ev << 12) | mi
        raw = np.nonzero(z["method"] == 0)[0]
        rk = key[raw]
        uk, first, cnt = np.unique(rk, return_index=True, return_counts=True)
        r1 = raw[first]                       # the row of keys with exactly one rawsum row
        ok_raw = (cnt == 1) & (z["npk"][r1] == 1) & ~z["clipped"][r1] & (z["E_units"][r1] > 0) \
            & (z["onset"][r1] >= 3)
        wr = z["wave_row"]
        Y = z["wave"]
        ns = Y.shape[1]
        j = np.clip(np.searchsorted(uk, key[wr]), 0, max(uk.size - 1, 0))
        hit = (uk[j] == key[wr]) & ok_raw[j]
        hit &= (z["onset"][r1[j]] <= ns - 4)
        hit &= (z["method"][wr] == 1) & (z["npk"][wr] == 1) & ~z["clipped"][wr] \
            & ~z["wa_fallback"][wr] & (z["amp"][wr] >= 50.0)
        hit &= np.array([n.startswith("W") for n in names])[mi[wr]]     # PbWO4
        t0 = z["onset"][wr] * DT
        w0 = np.ceil(t0 / DT).astype(np.int64) - LEAD
        hit &= (w0 >= 0) & (w0 + LW <= ns)
        k = np.nonzero(hit)[0]
        rr = wr[k]
        parts.append(dict(
            wave=Y[k], sigma=z["wave_sigma"][k].astype(np.float64), midx=mi[rr],
            t0=t0[k], amp=z["amp"][rr], time=z["time_ns"][rr], chi2=z["chi2dof"][rr],
            E=z["E_units"][rr] / kappa * DT, raw=z["E_units"][r1[j[k]]] * DT,
            split=np.full(k.size, ip)))
    S = {k: np.concatenate([q[k] for q in parts]) for k in parts[0]}
    # cap per module, seeded random choice
    rng = np.random.default_rng(seed)
    order = rng.permutation(len(S["midx"]))
    o2 = order[np.argsort(S["midx"][order], kind="stable")]
    ms = S["midx"][o2]
    rank = np.arange(len(o2)) - np.searchsorted(ms, ms)
    keep = np.sort(o2[rank < max_per_module])
    S = {k: v[keep] for k, v in S.items()}
    S.update(names=names, kappa=kappa, run=run)
    return S


# ------------------------------------------------------------------ band stack --
def resample_rows(tg, X, T_end):
    """Rows of X on tg (0.1 ns, cubic) -> the bank LUT grid 0, 0.05, ..., T_end (float32)."""
    from scipy.interpolate import CubicSpline
    n = int(round(T_end / LUT_DT)) + 1
    t = LUT_DT * np.arange(n)
    return CubicSpline(tg, X, axis=1)(t).astype(np.float32)


def iint_rows(rows, names_u, tail, T_end):
    """I_cont (the bank's energy integral) of unit-peak LUT rows [nu, n_lut]."""
    b = M.SampledTemplateBank(list(rows), names_u, np.ones(len(names_u), bool), T_end=T_end, **tail)
    return b.I_cont.copy()


class Stack:
    """Shapes of the compact module set `used`: stack index 0 = module template (the bank's own
    LUT and slope), b + 1 = band b's template, or the module template where the module has none."""

    def __init__(self, bank, zt, used, tail):
        self.used = used
        self.nu = len(used)
        self.n_lut, self.T_end = bank.n_lut, bank.T_end
        names_u = [bank.names[m] for m in used]
        row = {str(n): r for r, n in enumerate(zt["names"])}
        rr = np.array([row[n] for n in names_u])
        sb = zt["shape_band"][:, rr].astype(np.float64)             # [nb, nu, nt]
        am = zt["A_med_band"][:, rr].astype(np.float64)
        nb = sb.shape[0]
        usable = np.isfinite(sb).all(axis=2) & np.isfinite(am)
        enough = usable.sum(axis=1) >= MIN_BAND_MODULES
        usable &= enough[:, None]
        self.usable = usable
        self.kb = np.nonzero(enough)[0]                              # bands that are knots
        self.knot = np.array([np.median(am[b][usable[b]]) for b in self.kb])
        if (np.diff(self.knot) <= 0).any():
            sys.exit(f"FATAL: band knots not increasing: {self.knot}")
        self.edges = zt["band_edges"]
        nt = zt["T_end"]
        lut = np.empty((nb + 1, self.nu, self.n_lut), np.float32)
        dlut = np.empty_like(lut)
        lut[0], dlut[0] = bank.lut[used], bank.dlut[used]
        for b in range(nb):
            if b in self.kb:
                v = np.where(usable[b][:, None], sb[b], 0.0)
                lut[b + 1] = resample_rows(zt["t_grid"], v, float(nt))
                dlut[b + 1] = np.gradient(lut[b + 1], LUT_DT, axis=1)
            else:
                lut[b + 1], dlut[b + 1] = lut[0], dlut[0]
                continue
            bad = ~usable[b]
            lut[b + 1][bad], dlut[b + 1][bad] = lut[0][bad], dlut[0][bad]
        self.lut, self.dlut = lut.reshape(nb + 1, -1), dlut.reshape(nb + 1, -1)
        self.Iband = np.empty((nb + 1, self.nu))                     # I_cont of every stack entry
        self.Iband[0] = bank.I_cont[used]
        Ig = {}
        for b in self.kb:
            Ig[b] = iint_rows(lut[b + 1], names_u, tail, self.T_end)
        # grid-resampled module template vs the bank's (spline) integral: sanity
        self.I_resamp_mod = iint_rows(resample_rows(zt["t_grid"], zt["shape"][rr].astype(np.float64),
                                                    float(nt)), names_u, tail, self.T_end) / self.Iband[0]
        for b in range(nb):
            self.Iband[b + 1] = Ig.get(b, self.Iband[0])
        # tables over the interpolation weight, per knot interval
        nk = len(self.kb)
        self.nint = max(nk - 1, 1)
        self.T = np.empty((self.nint, WGRID + 1, self.nu))
        for k in range(self.nint):
            b0, b1 = self.kb[k], self.kb[min(k + 1, nk - 1)]
            self.T[k, 0], self.T[k, WGRID] = self.Iband[b0 + 1], self.Iband[b1 + 1]
            for iw in range(1, WGRID):
                w = iw / WGRID
                self.T[k, iw] = iint_rows((1 - w) * lut[b0 + 1] + w * lut[b1 + 1], names_u, tail, self.T_end) \
                    if nk > 1 else self.Iband[b0 + 1]

    def weights(self, mi_c, amp):
        """Per pulse (interval k, lower/upper stack index, weight w) of the ln A interpolation."""
        la, lk = np.log(np.maximum(amp, 1.0)), np.log(self.knot)
        nk = len(self.kb)
        if nk == 1:
            k = np.zeros(len(amp), np.int64)
            w = np.zeros(len(amp))
            iL = iH = np.full(len(amp), self.kb[0] + 1)
            return k, iL, iH, w
        k = np.clip(np.searchsorted(lk, la) - 1, 0, nk - 2)
        w = np.clip((la - lk[k]) / (lk[k + 1] - lk[k]), 0.0, 1.0)
        return k, self.kb[k] + 1, self.kb[k + 1] + 1, w

    def iint(self, k, w, mi_c):
        x = w * WGRID
        i0 = np.minimum(np.floor(x).astype(np.int64), WGRID - 1)
        f = x - i0
        return self.T[k, i0, mi_c] * (1 - f) + self.T[k, i0 + 1, mi_c] * f

    def both_real(self, iL, iH, mi_c):
        return (iL > 0) & (iH > 0) & self.usable[np.maximum(iL - 1, 0), mi_c] \
            & self.usable[np.maximum(iH - 1, 0), mi_c]

    def eval(self, iL, iH, w, mi_c, tau, slope=True):
        """Interpolated unit shape (and slope) at module rows mi_c (compact), times tau [n, nw]."""
        x = tau * (1.0 / LUT_DT)
        i = np.clip(np.floor(x), 0, self.n_lut - 2)
        f = x - i
        inside = (tau > 0) & (tau < self.T_end)
        base = (mi_c[:, None].astype(np.int64) * self.n_lut + i.astype(np.int64))
        out = []
        for tab in ((self.lut, self.dlut) if slope else (self.lut,)):
            v = []
            for idx in (iL, iH):
                v.append(tab[idx[:, None], base] * (1 - f) + tab[idx[:, None], base + 1] * f)
            out.append(np.where(inside, (1 - w[:, None]) * v[0] + w[:, None] * v[1], 0.0))
        return out if slope else out[0]


def t50_of(stk, iL, iH, w, mi_c, nsamp=2400):
    """50 % leading-edge time of the interpolated shapes (first crossing on the LUT grid)."""
    out = np.empty(len(w))
    for a in range(0, len(w), 4000):
        sl = slice(a, a + 4000)
        base = mi_c[sl, None].astype(np.int64) * stk.n_lut + np.arange(nsamp)[None, :]
        s = (1 - w[sl, None]) * stk.lut[iL[sl, None], base] + w[sl, None] * stk.lut[iH[sl, None], base]
        i1 = np.argmax(s >= 0.5, axis=1)
        i0 = np.maximum(i1 - 1, 0)
        ar = np.arange(s.shape[0])
        s0, s1 = s[ar, i0], s[ar, i1]
        fr = np.where(s1 > s0, (0.5 - s0) / np.maximum(s1 - s0, 1e-30), 0.0)
        out[sl] = (i0 + fr * (i1 - i0)) * LUT_DT
    return out


# ------------------------------------------------------------------------- fit --
def fit_arm(stk, iL, iH, w, mi_c, y, sig, idx, valid, t0, a0):
    """Gauss-Newton on (A, t0) in a fixed window; unweighted LSQ.  Returns t0, A, chi2dof."""
    xt = DT * idx.astype(np.float64)
    yy = np.take_along_axis(y, np.clip(idx, 0, y.shape[1] - 1), axis=1)
    ww = valid.astype(np.float64)
    t0, a = t0.copy(), a0.copy()
    for _ in range(NIT):
        s, ds = stk.eval(iL, iH, w, mi_c, xt - t0[:, None])
        g = -a[:, None] * ds
        res = yy - a[:, None] * s
        Gss, Gsg, Ggg = (ww * s * s).sum(1), (ww * s * g).sum(1), (ww * g * g).sum(1)
        bs, bg = (ww * s * res).sum(1), (ww * g * res).sum(1)
        det = Gss * Ggg - Gsg * Gsg
        good = det > 1e-6 * np.maximum(Gss * Ggg, 1e-20)
        dsafe = np.where(good, det, 1.0)
        da = np.where(good, (bs * Ggg - bg * Gsg) / dsafe, bs / np.maximum(Gss, 1e-20))
        dt = np.where(good, (Gss * bg - Gsg * bs) / dsafe, 0.0)
        a = np.maximum(a + da, 0.2 * a)
        t0 = t0 + np.clip(dt, -2.0, 2.0)
    s = stk.eval(iL, iH, w, mi_c, xt - t0[:, None], slope=False)
    a = (ww * s * yy).sum(1) / np.maximum((ww * s * s).sum(1), 1e-20)
    mdl = a[:, None] * s
    var = sig[:, None] ** 2 + (MODEL_FLOOR * mdl) ** 2
    chi2 = (ww * (yy - mdl) ** 2 / var).sum(1) / np.maximum(ww.sum(1) - 2.0, 1.0)
    return t0, a, chi2


def refit(S, bank, zt, tail, ped_nsamples):
    """Arm (a) and (b) refits of the selected pulses of one run."""
    used = np.unique(S["midx"])
    cmap = np.full(len(bank.names), -1, np.int64)
    cmap[used] = np.arange(len(used))
    mi_c = cmap[S["midx"]]
    stk = Stack(bank, zt, used, tail)
    print(f"[ampshape] run {S['run']}: bands as knots (ADC) "
          + ", ".join(f"{int(stk.edges[b])}-{int(stk.edges[b + 1])}: A={k:.0f} "
                      f"({int(stk.usable[b].sum())} modules)" for b, k in zip(stk.kb, stk.knot)))
    rec = M.MFReconstructor(bank, ped_nsamples=ped_nsamples)
    n = len(mi_c)
    out = {k: np.empty(n) for k in ("Aa", "Ab", "ta", "tb", "c2a", "c2b", "Qa", "Qb", "Ia", "Ib",
                                    "t50a", "t50b", "t0a", "t0b", "cov")}
    for c0 in range(0, n, CHUNK):
        sl = slice(c0, min(c0 + CHUNK, n))
        Yf = S["wave"][sl].astype(np.float32)
        ped, _ = rec._pedestal(Yf)
        y = (Yf - ped[:, None]).astype(np.float64)
        sig, m_c, m_b = S["sigma"][sl], mi_c[sl], S["midx"][sl]
        t0d = S["t0"][sl]
        w0 = np.ceil(t0d / DT).astype(np.int64) - LEAD
        idx = w0[:, None] + np.arange(LW)[None, :]
        valid = (idx >= 0) & (idx < y.shape[1])
        zero = np.zeros(len(m_c), np.int64)
        wz = np.zeros(len(m_c))
        # (a) module template (stack 0) + f(A)
        A0, _ = rec.amp_at(y.astype(np.float32), m_b, t0d.astype(np.float32))
        A0 = np.asarray(A0, np.float64)
        ta, Aa, c2a = fit_arm(stk, zero, zero, wz, m_c, y, sig, idx, valid, t0d, A0)
        Ia = stk.Iband[0][m_c] * bank.f_amp(Aa, m_b)
        # (b) interpolated band shapes, no f(A)
        tb, Ab = ta, Aa
        for _ in range(NPASS_B):
            k, iL, iH, w = stk.weights(m_c, Ab)
            tb, Ab, c2b = fit_arm(stk, iL, iH, w, m_c, y, sig, idx, valid, tb, Ab)
        k, iL, iH, w = stk.weights(m_c, Ab)
        Ib = stk.iint(k, w, m_c)
        t50a = t50_of(stk, zero, zero, wz, m_c)
        t50b = t50_of(stk, iL, iH, w, m_c)
        for key, v in (("Aa", Aa), ("Ab", Ab), ("ta", ta + t50a), ("tb", tb + t50b), ("c2a", c2a),
                       ("c2b", c2b), ("Qa", Aa * Ia * DT),
                       ("Qb", Ab * Ib * DT), ("Ia", Ia), ("Ib", Ib), ("t50a", t50a),
                       ("t50b", t50b), ("t0a", ta), ("t0b", tb),
                       ("cov", stk.both_real(iL, iH, m_c).astype(float))):
            out[key][sl] = v
    out["I_resamp_mod"] = stk.I_resamp_mod
    out["knots"] = [(int(stk.edges[b]), int(stk.edges[b + 1]), float(k), int(stk.usable[b].sum()))
                    for b, k in zip(stk.kb, stk.knot)]
    return out


# ------------------------------------------------------------------- reporting --
def fmt(x, p=4):
    return "nan" if not np.isfinite(x) else f"{x:.{p}f}"


def summarise(run, rt, S, R, tsv, md, plots_dir):
    A = S["amp"]
    bi = np.digitize(A, EDGES) - 1
    labels = [f"{lo:g}-{hi:g}" for lo, hi in zip(EDGES[:-1], EDGES[1:])]
    qa, qb = R["Qa"] / S["raw"], R["Qb"] / S["raw"]
    dt = R["tb"] - R["ta"]
    cov = R["cov"] > 0
    md.append(f"## run {run:06d} ({rt})")
    md.append("")
    md.append(f"{len(A)} pulses over {len(np.unique(S['midx']))} modules.  Band knots (ADC, "
              "median A_med, modules): " + "; ".join(f"{a}-{b}: {k:.0f} ({n})"
                                                       for a, b, k, n in R["knots"]) + ".")
    qd = R["Qa"] / S["E"]
    md.append(f"Sanity (a) vs the dump's mf row: Q_a/Q_dump median {np.median(qd):.5f}, "
              f"1-99 % {np.percentile(qd, 1):.5f}-{np.percentile(qd, 99):.5f}; A_a/A_dump median "
              f"{np.median(R['Aa'] / A):.5f}; t_a - t_dump median "
              f"{np.median(R['ta'] - S['time']):.4f} ns (|.| 99 % "
              f"{np.percentile(np.abs(R['ta'] - S['time']), 99):.3f}); chi2/dof_a / chi2/dof_dump "
              f"median {np.median(R['c2a'] / np.maximum(S['chi2'], 1e-9)):.4f}.  Grid-resampled "
              f"module shape I_int / bank (spline) I_int: median {np.median(R['I_resamp_mod']):.5f}, "
              f"min-max {R['I_resamp_mod'].min():.5f}-{R['I_resamp_mod'].max():.5f}.")
    md.append("")
    res = {}
    for pop, sel in (("all", np.ones(len(A), bool)), ("cov", cov)):
        st = {"a": [], "b": []}
        md.append(f"Population {pop} ({'every pulse' if pop == 'all' else 'both interpolation sides are true band templates'}; "
                  f"{int(sel.sum())} pulses).")
        md.append("")
        md.append("| A bin (ADC) | n | cov % | Q_a/Q_raw +- SE | Q_b/Q_raw +- SE | chi2/dof a | "
                  "chi2/dof b | median t_b - t_a (ns) |")
        md.append("|---|---|---|---|---|---|---|---|")
        for k, lab in enumerate(labels):
            m = sel & (bi == k)
            sa, sb_ = bin_med_se(qa[m]), bin_med_se(qb[m])
            st["a"].append(sa)
            st["b"].append(sb_)
            c2a = float(np.median(R["c2a"][m])) if m.any() else NAN
            c2b = float(np.median(R["c2b"][m])) if m.any() else NAN
            d = float(np.median(dt[m])) if m.any() else NAN
            cf = 100.0 * float(cov[bi == k].mean()) if (bi == k).any() else NAN
            tsv.append((run, rt, pop, lab, int(m.sum()), cf, sa[1], sa[2], sb_[1], sb_[2], c2a, c2b, d))
            md.append(f"| {lab} | {int(m.sum())} | {cf:.0f} | "
                      f"{fmt(sa[1])} +- {fmt(sa[2])} | {fmt(sb_[1])} +- {fmt(sb_[2])} | {fmt(c2a, 3)} | "
                      f"{fmt(c2b, 3)} | {fmt(d, 3)} |")
        fa, fb = flatness(st["a"], EDGES), flatness(st["b"], EDGES)
        res[pop] = (fa, fb)
        md.append("")
        for arm, f in (("a", fa), ("b", fb)):
            md.append(f"- arm ({arm}) flatness over the bins >= 100 ADC: "
                      + ("n/a (< 2 bins)" if f is None else
                         f"{'FLAT' if f[0] else 'NOT flat'}, max pairwise pull {f[2]:.2f} SE-units "
                         f"(rule <= 2), spread {f[1]:.3f} %"))
        if fa and fb:
            md.append(f"- pass (b flat where a is not): "
                      f"{'YES' if (fb[0] and not fa[0]) else 'NO'}; b spread "
                      f"{'smaller' if fb[1] < fa[1] else 'not smaller'} than a "
                      f"({fb[1]:.3f} vs {fa[1]:.3f} %)")
        md.append("")
    plot_run(run, rt, A, qa, qb, cov, plots_dir)
    return res


def plot_run(run, rt, A, qa, qb, cov, out_dir):
    edges = np.geomspace(50, 4000, 14)
    xc = np.sqrt(edges[:-1] * edges[1:])
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2), sharey=True)
    for x, (ttl, sel) in zip(ax, (("all pulses", np.ones(len(A), bool)),
                                  ("both sides true band templates", cov))):
        for q, lab, col in ((qa, "(a) module + f(A)", "#1f77b4"), (qb, "(b) ln A band interp.", "#d95f02")):
            st = [bin_med_se(q[sel & (A >= lo) & (A < hi)]) for lo, hi in zip(edges[:-1], edges[1:])]
            x.errorbar(xc, [s[1] for s in st], [s[2] for s in st], fmt="o-", ms=3, lw=1, color=col, label=lab)
        for e in (100, 250, 500, 1000, 2000):
            x.axvline(e, color="gray", lw=.5, ls=":")
        x.set(xscale="log", xlabel="A (ADC)", title=ttl)
    ax[0].set_ylabel("median Q_template / Q_raw (+- SE)")
    ax[0].legend(fontsize=8)
    fig.suptitle(f"run {run:06d} ({rt}): option (ii) offline refit", fontsize=9)
    fig.tight_layout()
    fig.savefig(out_dir / f"ampshape_ratio_{run:06d}.png", dpi=110)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dumps", required=True, help="glob of study-pulse dumps (pulses.npz)")
    ap.add_argument("--templates", required=True,
                    help="bank npz pattern with {runtype} (the dump run's run type)")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--max-per-module", type=int, default=200,
                    help="pulses per module and run (default 200)")
    ap.add_argument("--seed", type=int, default=20261025)
    args = ap.parse_args()
    paths = sorted(glob.glob(args.dumps))
    if not paths:
        sys.exit(f"no dumps match {args.dumps}")
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    from prad2py import dec
    from gate_scan_20261022 import runtype_slug
    wcfg = dec.load_daq_config().wave_cfg
    tail = dict(smooth_order=int(wcfg.smooth_order), int_tail_ratio=float(wcfg.int_tail_ratio),
                tail_break_n=int(wcfg.tail_break_n))
    byrun = {}
    for p in paths:
        byrun.setdefault(int(np.load(p)["run"]), []).append(p)
    tsv, verdict = [], []
    md = ["# option (ii) offline feasibility: amplitude-indexed templates "
                            "(plan_20261025 section 4.4)", "",
                            "(a) module template + f(A); (b) shape interpolated linearly in ln A "
                            "between neighbouring band templates, no f(A).  Q/Q_raw per acceptance bin "
                            "(median +- SE, SE = 1.2533 * 1.4826 MAD / sqrt(n)); flatness: every pairwise "
                            "bin difference <= 2 hypot(SE).  "
                            "See the docstring of mf_ampshape_diag.py for the fit and the band rules.", ""]
    for run in sorted(byrun):
        rt = runtype_slug(run)
        S = load_run(byrun[run], args.max_per_module, args.seed)
        tp = Path(args.templates.replace("{runtype}", rt))
        if not tp.is_file():
            tp = HERE / tp
        zt = np.load(tp, allow_pickle=False)
        names = S["names"]
        pwo4 = np.array([n.startswith("W") for n in names])
        bank = M.SampledTemplateBank.from_npz(str(tp), names, pwo4, **tail)
        print(f"[ampshape] run {run}: {len(S['midx'])} pulses from {len(byrun[run])} splits, "
              f"bank {tp}")
        R = refit(S, bank, zt, tail, int(wcfg.ped_nsamples))
        verdict.append((run, rt, summarise(run, rt, S, R, tsv, md, out)))
    hdr = ("run", "runtype", "pop", "A_bin", "n", "cov_pct", "Qa_med", "Qa_se", "Qb_med", "Qb_se",
           "chi2dof_a", "chi2dof_b", "dt_b_minus_a_ns")
    with open(out / "ampshape_diag.tsv", "w") as f:
        f.write("\t".join(hdr) + "\n")
        for r in tsv:
            f.write("\t".join(f"{x:.6g}" if isinstance(x, float) else str(x) for x in r) + "\n")
    md += ["## Verdict", "", "| run | pop | a flat | a spread % | b flat | b spread % | pass | b smaller |",
           "|---|---|---|---|---|---|---|---|"]
    for run, rt, res in verdict:
        for pop, (fa, fb) in res.items():
            if fa is None or fb is None:
                continue
            md.append(f"| {run:06d} {rt} | {pop} | {'yes' if fa[0] else 'no'} | {fa[1]:.3f} | "
                      f"{'yes' if fb[0] else 'no'} | {fb[1]:.3f} | "
                      f"{'YES' if (fb[0] and not fa[0]) else 'no'} | {'yes' if fb[1] < fa[1] else 'no'} |")
    (out / "ampshape_diag.md").write_text("\n".join(md) + "\n")
    print(f"[ampshape] wrote {out}/ampshape_diag.tsv, ampshape_diag.md, ampshape_ratio_*.png")


if __name__ == "__main__":
    main()
