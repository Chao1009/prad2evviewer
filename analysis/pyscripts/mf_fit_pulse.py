#!/usr/bin/env python3
"""mf_fit_pulse.py -- per-pulse free-tau fit of the two-tau HyCal pulse (batched LM).

Companion of mf_pulse.py (pure NumPy, float64).  The matched filter identifies a
single, unclipped pulse and gives (A, t0); this module refits (A, t0, tau_r, tau_f)
of the same pulse function (mf_pulse.template, unit peak) on the same fixed sample
window as the mf K=1 closed form, seeded at the mf (A, t0) and the module's template
(tau_r, tau_f), constrained by the template's per-module spread sigma = 1.4826 * MAD
(load_tau_sigma):

  mode "bounded": tau within template +- nsig*sigma (projected LM, hard box)
  mode "prior"  : Gaussian prior ((tau - tau0)/sigma)^2 added to the chi2, only a
                  wide safety box (tau_r <= 15 ns, 5 <= tau_f <= 80 ns); ending on
                  it -> status 3

Model f_i = A s(4 i - t0; tau_r, tau_f), window idx = ceil(t0_seed/4) - LEAD +
arange(LW) (fixed at the seed).  Weights are fixed per pulse (option `weight`):
  "uniform" (default): one variance for the whole window, as the h50 extraction
      (WaveAnalyzer fit_pulse_shape): v = max(max(sigma, 1), floor * a_seed)^2
  "seedmodel": v_i = sigma^2 + (floor * a_seed * s_seed_i)^2 (seed-model shape)
Objective sum (y-f)^2/v (+ prior, in the same chi2 units).  chi2dof, chi2dof_seed
and the kink-escape triggers always use the seed-model variance (the mf definition),
whatever `weight` is, so chi2 stays comparable with the mf p_chi2.

Levenberg-Marquardt (N + lambda diag(N)), per-row lambda, batched.  A is eliminated
(variable projection): every trial point re-solves A as the weighted LS optimum for the
other parameters.  A step is accepted only if the objective decreases; lambda follows
the gain ratio (actual / model decrease: x0.3 above 0.75, x3 below 0.25, x5 on a
rejected step, floor 1e-9).  Steps: |dt0| <= 2 ns, no tau changes by more than 30%
(the whole step is scaled), box clip with bound parameters frozen (also when the coupled
step would push them out).  A row converges when the model-predicted remaining decrease
(2 g.u - u^T N u at the undamped box-projected Gauss-Newton step u; negative = not
converged) is <= 1e-6 * objective (+1e-10), or when a
rejected step is already < 1e-4 ns in t0, tau_r, tau_f, or when lambda exceeds 1e8
(stalled).  Not converged after max_iter (default 60) -> status 1.
The model has a kink in t0 at every sample time (zero for t <= t0), which traps LM
next to the kink; converged rows near one are restarted across it ("kink escape",
see FreeTauFitter._lm) and keep the lower objective.

Uncertainties (per pulse, noise only): sandwich covariance
C = H^-1 (J^T V^-1 sigma^2 V^-1 J) H^-1, H = J^T V^-1 J (+ prior diag, prior arm);
sigw_* is the as-fit sqrt(diag((J^T V^-1 J)^-1)) (likelihood only).  For rows sitting
on a bound of the bounded arm these are the unconstrained local values.  The pedestal
is not a fit parameter: its estimation error is not in sig_*.

Energy: E = A * I_int(tau_r_fit, tau_f_fit, phase) * kappa, I_int the WaveAnalyzer
walk-out unit integral of the fitted shape at u = mod(-t0_fit, 4 ns) (A-independent,
same convention as TemplateBank.iint).

Failure (status != 0): amp/t0/tau/I_int are the seed values (amp = fixed-template LS
amplitude at the seed t0, I_int = bank.iint(m, t0)); sig_*, sigw_*, chi2dof are NaN
and at_* are False.  status: 0 ok, 1 not converged at max_iter, 2 non-finite /
singular / unusable seed, 3 ended on the prior arm's safety box.

    python3 mf_fit_pulse.py --self-test
"""

from __future__ import annotations

import math
import os
import sys
import time

import numpy as np

import mf_pulse
from mf_pulse import DT, LEAD, LW, template, tmax_of

TAU_R_SAFE_HI = 15.0               # prior arm safety box (tau_r lower edge = tau_r_min)
TAU_F_SAFE = (5.0, 80.0)
TOL_STEP = 1e-4                    # ns
TOL_OBJ = 1e-6                     # predicted remaining decrease / objective
LAM_STALL = 1e8
TAU_REL = 0.3                      # max relative tau change per LM step
BOUND_EPS = 1e-6                   # ns, "at the bound"
KINK_NSIG = 3.0                    # kink-escape trigger (see FreeTauFitter._lm)
ST_OK, ST_NOCONV, ST_BAD, ST_BOX = 0, 1, 2, 3


# ----------------------------------------------------------------- sigmas --
def load_tau_sigma(table_path, column, names, pwo4):
    """Per-module (sig_r, sig_f) = 1.4826*MAD from tau_{r,f}_sigma_<column>.

    Same rule as mf_pulse.load_h50: for column 25320 a NaN mean pair takes the
    25308 column (its sigma too); modules that load_h50 does not accept (not in the
    table, untrusted, out of range) or without a finite positive sigma, and non
    PbWO4 modules, get the median over the modules that have one.
    Returns (sig_r[nm], sig_f[nm], info); info n_fallback counts every module that
    got the median (non-PbWO4 ones included).
    """
    column = str(column)
    rows, _, _ = mf_pulse.load_h50(table_path, column)
    cols, tab = None, {}
    with open(table_path) as fh:
        for line in fh:
            if line.startswith("# module"):
                cols = [c.strip() for c in line[2:].split("\t")]
                continue
            if cols is None or not line.strip() or line.startswith("#"):
                continue
            f = line.split()
            tab[f[0]] = dict(zip(cols, f))
    other = "25308" if (column == "25320" and cols and "tau_r_mean_25308" in cols) else None
    nm = len(names)
    sr = np.full(nm, np.nan)
    sf = np.full(nm, np.nan)
    for i, (name, is_p) in enumerate(zip(names, pwo4)):
        if not is_p or name not in rows or name not in tab:
            continue
        d = tab[name]
        col = column
        if other and not (math.isfinite(float(d["tau_r_mean_" + column]))
                          and math.isfinite(float(d["tau_f_mean_" + column]))):
            col = other
        a, b = float(d["tau_r_sigma_" + col]), float(d["tau_f_sigma_" + col])
        if math.isfinite(a) and math.isfinite(b) and a > 0 and b > 0:
            sr[i], sf[i] = 1.4826 * a, 1.4826 * b
    have = np.isfinite(sr) & np.isfinite(sf)
    if not have.any():
        raise ValueError(f"{table_path}: no finite tau sigma in column {column}")
    msr, msf = float(np.median(sr[have])), float(np.median(sf[have]))
    sr = np.where(have, sr, msr)
    sf = np.where(have, sf, msf)
    info = dict(n_module=int(have.sum()), n_fallback=int(nm - have.sum()),
                median_sig_r=msr, median_sig_f=msf)
    return sr, sf, info


# ----------------------------------------------------------------- result --
class FitResult:
    """Plain attribute container of the per-pulse fit outputs (arrays of length n)."""

    def __init__(self, n):
        for k in ("amp", "t0", "tau_r", "tau_f", "sig_tau_r", "sig_tau_f", "sig_amp",
                  "sigw_tau_r", "sigw_tau_f", "chi2dof", "chi2dof_seed", "I_int"):
            setattr(self, k, np.full(n, np.nan))
        self.niter = np.zeros(n, dtype=np.int64)
        self.status = np.zeros(n, dtype=np.int8)
        for k in ("at_lo_r", "at_hi_r", "at_lo_f", "at_hi_f"):
            setattr(self, k, np.zeros(n, dtype=bool))


# ---------------------------------------------------------- model and jacobian --
def _model_jac(P, xt):
    """f [k, LW] and J [k, LW, 4] (A, t0, tau_r, tau_f) for P [k, 4], sample times xt."""
    A, t0, tr, tf = P[:, 0:1], P[:, 1:2], P[:, 2:3], P[:, 3:4]
    x = np.maximum(xt - t0, 0.0)
    pos = (xt - t0) > 0.0
    e1 = np.exp(-x / tr)
    e2 = np.exp(-x / tf)
    tm = tmax_of(tr, tf)
    tpk = tr * np.log1p(tf / tr)
    dlr = -np.exp(-tpk / tr) * tpk / (tr * tr) * np.exp(-tpk / tf) / tm   # dln(tm)/dtau_r
    dlf = tpk / (tf * tf)                                                # dln(tm)/dtau_f
    g = (1.0 - e1) * e2
    s = g / tm
    J = np.empty(P.shape[:1] + xt.shape[1:] + (4,))
    J[..., 0] = s
    J[..., 1] = np.where(pos, -A * (e1 / tr - (1.0 - e1) / tf) * e2 / tm, 0.0)
    J[..., 2] = A * (-e1 * x / (tr * tr) * e2 / tm - s * dlr)
    J[..., 3] = A * (g * x / (tf * tf) / tm - s * dlf)
    return A * s, J


def _solve(M, g):
    """Batched solve; exactly singular rows fall back to a per-row try (NaN)."""
    try:
        return np.linalg.solve(M, g[..., None])[..., 0]
    except np.linalg.LinAlgError:
        out = np.zeros_like(g)
        for i in range(M.shape[0]):
            try:
                out[i] = np.linalg.solve(M[i], g[i])
            except np.linalg.LinAlgError:
                out[i] = np.nan
        return out


def _inv(M):
    try:
        return np.linalg.inv(M)
    except np.linalg.LinAlgError:
        out = np.full_like(M, np.nan)
        for i in range(M.shape[0]):
            try:
                out[i] = np.linalg.inv(M[i])
            except np.linalg.LinAlgError:
                pass
        return out


# ------------------------------------------------------------------ fitter --
class FreeTauFitter:
    def __init__(self, bank, sig_r, sig_f, *, mode, nsig=3.0, max_iter=60,
                 floor=mf_pulse.MODEL_FLOOR, tau_r_min=0.25, weight="uniform"):
        if mode not in ("bounded", "prior"):
            raise ValueError("mode must be 'bounded' or 'prior'")
        if weight not in ("uniform", "seedmodel"):
            raise ValueError("weight must be 'uniform' or 'seedmodel'")
        self.weight = weight
        self.b = bank
        self.sig_r = np.asarray(sig_r, dtype=np.float64)
        self.sig_f = np.asarray(sig_f, dtype=np.float64)
        self.mode, self.nsig = mode, float(nsig)
        self.max_iter, self.floor = int(max_iter), float(floor)
        self.tau_r_min = float(tau_r_min)

    # -- helpers --------------------------------------------------------------
    def _box(self, m):
        """(tau_r lo, tau_r hi, tau_f lo, tau_f hi) per row: the bounded-arm box or
        the prior arm's safety box."""
        b, rm = self.b, self.tau_r_min
        tr0, tf0 = b.tau_r[m], b.tau_f[m]
        if self.mode == "bounded":
            sr, sf = self.nsig * self.sig_r[m], self.nsig * self.sig_f[m]
            return (np.maximum(tr0 - sr, rm), tr0 + sr,
                    np.maximum(tf0 - sf, 2.0 * rm), tf0 + sf)
        k = np.ones_like(tr0)
        return rm * k, TAU_R_SAFE_HI * k, TAU_F_SAFE[0] * k, TAU_F_SAFE[1] * k

    def _iint(self, tr, tf, t0):
        b = self.b
        return mf_pulse.unit_integral(tr, tf, tmax_of(tr, tf), np.mod(-t0, DT),
                                      b.smooth_order, b.int_tail_ratio, b.tail_break_n)

    def energy_units(self, fr, kappa=1.0):
        return fr.amp * fr.I_int * kappa

    # -- fit ------------------------------------------------------------------
    def fit(self, y, m, t0, a, sigma):
        b = self.b
        y = np.asarray(y, dtype=np.float64)
        n, ns = y.shape
        m = np.asarray(m, dtype=np.int64)
        t0 = np.asarray(t0, dtype=np.float64)
        a = np.asarray(a, dtype=np.float64)
        sigma = np.maximum(np.asarray(sigma, dtype=np.float64), 1e-6)
        fr = FitResult(n)
        if n == 0:
            return fr
        tr0, tf0, tm0 = b.tau_r[m], b.tau_f[m], b.tmx[m]

        # fixed window (same as mf), seed model, fixed-template LS amplitude, chi2 seed
        w0 = np.ceil(t0.astype(np.float32) / DT).astype(np.int64) - LEAD
        idx = w0[:, None] + np.arange(LW)[None, :]
        valid = (idx >= 0) & (idx < ns)
        yy = np.take_along_axis(y, np.clip(idx, 0, ns - 1), axis=1)
        xt = DT * idx.astype(np.float64)
        nval = valid.sum(1)
        s0 = template(xt - t0[:, None], tr0[:, None], tf0[:, None], tm0[:, None])
        with np.errstate(all="ignore"):
            a_ls = (valid * s0 * yy).sum(1) / np.maximum((valid * s0 * s0).sum(1), 1e-20)
            mdl = a_ls[:, None] * s0
            v2 = sigma[:, None] ** 2 + (self.floor * mdl) ** 2
            fr.chi2dof_seed = ((valid * (yy - mdl) ** 2 / v2).sum(1)
                               / np.maximum(nval - 2.0, 1.0))
        # fit weights: fixed at the seed model of the given amplitude
        v = sigma[:, None] ** 2 + (self.floor * np.maximum(a, 0.0)[:, None] * s0) ** 2
        wq_rep = valid / np.sqrt(v)          # seed-model weights: chi2dof, kink triggers
        if self.weight == "seedmodel":
            wq = wq_rep                                             # 0 off-window
        else:
            # one variance per pulse, as the h50 extraction (WaveAnalyzer
            # fit_pulse_shape): max(max(sigma, 1), floor * A_seed)^2
            vu = np.maximum(np.maximum(sigma, 1.0),
                            self.floor * np.maximum(a, 0.0)) ** 2
            wq = valid / np.sqrt(vu)[:, None]

        usable = ((a > 0) & np.isfinite(a) & np.isfinite(t0) & (nval >= 5)
                  & np.isfinite(yy).all(1) & np.isfinite(a_ls))
        rows = np.nonzero(usable)[0]
        status = np.full(n, ST_BAD, dtype=np.int8)
        res = {}
        if rows.size:
            res = self._lm(yy[rows], xt[rows], wq[rows], wq_rep[rows], np.stack(
                [a[rows], t0[rows], tr0[rows], tf0[rows]], axis=1), m[rows],
                sigma[rows], nval[rows])
            status[rows] = res["status"]
            fr.niter[rows] = res["niter"]

        # outputs: fitted on success, seed otherwise
        ok = status == ST_OK
        fr.amp[:], fr.t0[:], fr.tau_r[:], fr.tau_f[:] = a_ls, t0, tr0, tf0
        fr.I_int[:] = b.iint(m, t0)
        if ok.any():
            ro = ok[rows]
            io = rows[ro]
            P = res["P"][ro]
            fr.amp[io], fr.t0[io], fr.tau_r[io], fr.tau_f[io] = P.T
            fr.I_int[io] = self._iint(P[:, 2], P[:, 3], P[:, 1])
            for k in ("sig_amp", "sig_tau_r", "sig_tau_f", "sigw_tau_r", "sigw_tau_f"):
                getattr(fr, k)[io] = res[k][ro]
            fr.chi2dof[io] = res["chi2dof"][ro]
            for k in ("at_lo_r", "at_hi_r", "at_lo_f", "at_hi_f"):
                getattr(fr, k)[io] = res[k][ro]
        fr.status = status
        return fr

    def _core(self, yy, xt, wq, P, lo, hi, pw, P0):
        """Projected LM on rows (P updated in place). Returns obj, niter, done."""
        k = P.shape[0]
        eye = np.eye(4)[None]

        def objective(f, Pq, sel):
            r = wq[sel] * (yy[sel] - f)
            return r, (r * r).sum(1) + (pw[sel] * (Pq - P0[sel]) ** 2).sum(1)

        def model_opt(Pq, sel):
            """Model and Jacobian at Pq with A replaced by its weighted LS optimum for
            the other parameters (variable projection: the zig-zag of A against the
            shape parameters in large-residual fits is gone)."""
            Pu = Pq.copy()
            Pu[:, 0] = 1.0
            f1, J1 = _model_jac(Pu, xt[sel])
            w2 = wq[sel] ** 2
            A = np.maximum((w2 * f1 * yy[sel]).sum(1)
                           / np.maximum((w2 * f1 * f1).sum(1), 1e-300), 1e-3)
            Pq[:, 0] = A
            J1[..., 1:] *= A[:, None, None]
            return A[:, None] * f1, J1

        with np.errstate(all="ignore"):
            f, J = model_opt(P, np.arange(k))
        r, obj = objective(f, P, np.arange(k))
        Jw = J * wq[:, :, None]
        lam = np.full(k, 1e-3)
        done = np.zeros(k, dtype=bool)
        niter = np.zeros(k, dtype=np.int64)
        for _ in range(self.max_iter):
            act = np.nonzero(~done)[0]
            if act.size == 0:
                break
            Pa, la, ha = P[act], lo[act], hi[act]
            Ja, ra = Jw[act], r[act]
            N = np.matmul(Ja.transpose(0, 2, 1), Ja)
            g = np.einsum("kli,kl->ki", Ja, ra) - pw[act] * (Pa - P0[act])
            N = N + pw[act][:, :, None] * eye
            dg = np.maximum(np.diagonal(N, axis1=1, axis2=2), 1e-10)
            # parameters on a bound with the gradient pushing outward are frozen
            frz = ((Pa <= la + 1e-12) & (g < 0.0)) | ((Pa >= ha - 1e-12) & (g > 0.0))
            fm = (~frz).astype(np.float64)
            Md = N + lam[act][:, None, None] * dg[:, :, None] * eye
            for _pass in range(3):
                M = Md * fm[:, :, None] * fm[:, None, :] + (1.0 - fm)[:, :, None] * eye
                d = _solve(M, g * fm)
                # a bound parameter whose coupled step still points outward is frozen
                out = (((Pa <= la + 1e-12) & (d < 0.0)) | ((Pa >= ha - 1e-12) & (d > 0.0))) \
                    & (fm > 0.0)
                if not out.any():
                    break
                fm = fm * (~out)
            # predicted remaining decrease of the objective: quadratic model at the
            # undamped (Gauss-Newton) projected step, 2 g.u - u^T N u.  Components
            # clipped by the box are held at the bound and the rest re-solved; a
            # negative value means no model descent: not converged.
            M0 = N + 1e-9 * dg[:, :, None] * eye
            fm0 = (M0 * fm[:, :, None] * fm[:, None, :]
                   + (1.0 - fm)[:, :, None] * eye)
            u1 = _solve(fm0, g * fm)
            u = np.clip(Pa + u1, la, ha) - Pa
            cl = (np.abs(u - u1) > 1e-12) & (fm > 0)
            if cl.any():
                fm2 = fm * (~cl)
                uc = u * cl
                rhs = (g - np.einsum("kij,kj->ki", N, uc)) * fm2
                m2 = (M0 * fm2[:, :, None] * fm2[:, None, :]
                      + (1.0 - fm2)[:, :, None] * eye)
                u = np.clip(Pa + _solve(m2, rhs) * fm2 + uc, la, ha) - Pa
            pred = (2.0 * g * u).sum(1) - np.einsum("ki,kij,kj->k", u, N, u)
            d[:, 1] = np.clip(d[:, 1], -2.0, 2.0)
            # trust region: scale the step so no tau changes by more than TAU_REL
            rel = np.max(np.abs(d[:, 2:]) / Pa[:, 2:], axis=1)
            d = d * np.minimum(1.0, TAU_REL / np.maximum(rel, 1e-300))[:, None]
            Pt = np.clip(Pa + d, la, ha)
            with np.errstate(all="ignore"):
                ft, Jt = model_opt(Pt, act)
                rt, ot = objective(ft, Pt, act)
            oa = obj[act]
            acc = np.isfinite(ot) & (ot < oa)
            step = Pt - Pa
            small = (np.abs(step[:, 1:]) < TOL_STEP).all(1)
            onew = np.where(acc, ot, oa)
            conv = ((pred >= 0.0) & (pred <= TOL_OBJ * onew + 1e-10)) | (small & ~acc)
            ia = act[acc]
            P[ia], r[ia], Jw[ia], obj[ia] = (Pt[acc], rt[acc],
                                             Jt[acc] * wq[ia][:, :, None], ot[acc])
            # damping from the gain ratio (actual / model decrease of the step taken):
            # good -> less damping, poor (zig-zag overshoot) -> more, rejected -> x5
            pstep = (2.0 * g * step).sum(1) - np.einsum("ki,kij,kj->k", step, N, step)
            with np.errstate(all="ignore"):
                rho = (oa - ot) / np.maximum(pstep, 1e-300)
            fac = np.where(acc, np.where(rho > 0.75, 0.3, np.where(rho < 0.25, 3.0, 1.0)),
                           5.0)
            lam[act] = np.clip(lam[act] * fac, 1e-9, None)
            niter[act] += 1
            done[act[conv]] = True
            done |= lam > LAM_STALL
        return obj, niter, done

    def _lm(self, yy, xt, wq, wq_rep, P, m, sigma, nval):
        """Batched projected LM on the rows that can be fitted (+ kink escape)."""
        k = P.shape[0]
        lo = np.full((k, 4), -np.inf)
        hi = np.full((k, 4), np.inf)
        lo[:, 2], hi[:, 2], lo[:, 3], hi[:, 3] = self._box(m)
        P0 = P.copy()
        P[:, 2:] = np.clip(P[:, 2:], lo[:, 2:], hi[:, 2:])
        pw = np.zeros((k, 4))
        if self.mode == "prior":
            pw[:, 2] = 1.0 / self.sig_r[m] ** 2
            pw[:, 3] = 1.0 / self.sig_f[m] ** 2
        eye = np.eye(4)[None]
        obj, niter, done = self._core(yy, xt, wq, P, lo, hi, pw, P0)

        # Kink escape.  The model is zero for x = t - t0 <= 0 and rises linearly at
        # x = 0+, so f(t0) has a kink at every sample time and LM can converge to a
        # local minimum next to one: (a) the sample nearest to onset on the inactive
        # side (x <= 0) has data > KINK_NSIG sigma above the zero model, the true onset
        # is earlier (the step cannot cross, J = 0 there); (b) the first active sample
        # (x > 0) has a model value < KINK_NSIG sigma, its curvature freezes the step
        # and the onset may be later.  Such rows restart just across the kink (A, tau
        # kept) and keep the result if the objective is lower.
        ar = np.arange(k)
        for _ in range(1):
            gained = False
            for side in (-1, 1):
                x = xt - P[:, 1:2]
                if side < 0:
                    xs = np.where((x <= 0.0) & (wq_rep > 0.0), -x, np.inf)   # distance
                else:
                    xs = np.where((x > 0.0) & (wq_rep > 0.0), x, np.inf)
                kk = np.argmin(xs, axis=1)
                xk = xs[ar, kk]
                with np.errstate(all="ignore"):
                    if side < 0:
                        trig = yy[ar, kk] * wq_rep[ar, kk] > KINK_NSIG
                    else:
                        fk = P[:, 0] * template(xk, P[:, 2], P[:, 3],
                                                tmax_of(P[:, 2], P[:, 3]))
                        trig = fk * wq_rep[ar, kk] < KINK_NSIG
                cand = np.nonzero(done & np.isfinite(xk) & trig)[0]
                if cand.size == 0:
                    continue
                P2 = P[cand].copy()
                P2[:, 1] += side * (xk[cand] + 1e-3)       # just across the kink
                o2, n2, d2 = self._core(yy[cand], xt[cand], wq[cand], P2, lo[cand],
                                        hi[cand], pw[cand], P0[cand])
                better = np.isfinite(o2) & (o2 < obj[cand] - 1e-9 - 1e-9 * obj[cand])
                ib = cand[better]
                P[ib], obj[ib], done[ib] = P2[better], o2[better], d2[better]
                niter[cand] += n2
                gained |= bool(better.any())
            if not gained:
                break
        status = np.where(done, ST_OK, ST_NOCONV).astype(np.int8)

        # final-point diagnostics
        with np.errstate(all="ignore"):
            f, J = _model_jac(P, xt)
            Jw = J * wq[:, :, None]
            H0 = np.matmul(Jw.transpose(0, 2, 1), Jw)
            Hp = H0 + pw[:, :, None] * eye
            Hi0, Hi = _inv(H0), _inv(Hp)
            Jm = Jw * wq[:, :, None]
            meat = sigma[:, None, None] ** 2 * np.matmul(Jm.transpose(0, 2, 1), Jm)
            C = np.matmul(np.matmul(Hi, meat), Hi)
            dC = np.sqrt(np.diagonal(C, axis1=1, axis2=2))
            d0 = np.sqrt(np.diagonal(Hi0, axis1=1, axis2=2))
            rr = wq_rep * (yy - f)
            chi2dof = (rr * rr).sum(1) / np.maximum(nval - 4.0, 1.0)
        out = dict(P=P, niter=niter, chi2dof=chi2dof,
                   sig_amp=dC[:, 0], sig_tau_r=dC[:, 2], sig_tau_f=dC[:, 3],
                   sigw_tau_r=d0[:, 2], sigw_tau_f=d0[:, 3])
        bad = ~(np.isfinite(P).all(1) & np.isfinite(dC).all(1) & np.isfinite(d0).all(1)
                & np.isfinite(chi2dof))
        status[bad & (status == ST_OK)] = ST_BAD
        at = dict(at_lo_r=P[:, 2] <= lo[:, 2] + BOUND_EPS,
                  at_hi_r=P[:, 2] >= hi[:, 2] - BOUND_EPS,
                  at_lo_f=P[:, 3] <= lo[:, 3] + BOUND_EPS,
                  at_hi_f=P[:, 3] >= hi[:, 3] - BOUND_EPS)
        if self.mode == "prior":
            box = np.zeros(k, dtype=bool)
            for v_ in at.values():
                box |= v_
            status[box & (status == ST_OK)] = ST_BOX
            at = {kk_: np.zeros(k, dtype=bool) for kk_ in at}
        out.update(at)
        out["status"] = status
        return out


# ---------------------------------------------------------------- self test --
def _wave(t0, amp, tr, tf, ns=100, noise=0.0, rng=None, ped=0.0):
    t = DT * np.arange(ns)[None, :]
    y = ped + amp[:, None] * template(t - t0[:, None], tr[:, None], tf[:, None],
                                      tmax_of(tr, tf)[:, None])
    if noise > 0:
        y = y + rng.normal(0.0, noise, y.shape)
    return y


def _seeds(rec, Y, mi):
    """Run the mf; keep rows with exactly one pulse.
    Returns (rows, m, t0, a, sigma, ped, p_chi2)."""
    r = rec.run(Y.astype(np.float32), mi)
    rows = np.nonzero(r.npk == 1)[0]
    sel = np.searchsorted(r.p_ch, rows)
    return (rows, mi[rows], r.p_t0[sel], r.p_amp[sel], r.sigma[rows].astype(np.float64),
            r.ped[rows].astype(np.float64), r.p_chi2[sel])


def _self_test_w(weight):
    FTF = lambda *a, **k: FreeTauFitter(*a, weight=weight, **k)
    rng = np.random.default_rng(20261023)
    taus = [(2.499, 24.41), (2.0, 20.0), (3.2, 28.0)]
    nm = 3
    bank = mf_pulse.TemplateBank([t[0] for t in taus], [t[1] for t in taus], [True] * nm,
                                 smooth_order=2, int_tail_ratio=0.1, tail_break_n=2)
    rec = mf_pulse.MFReconstructor(bank)
    sr = np.array([0.6, 0.7, 0.8])
    sf = np.array([1.5, 1.2, 2.0])
    ok_all = True

    def report(tag, passed, txt):
        nonlocal ok_all
        ok_all &= bool(passed)
        print(f"[{'PASS' if passed else 'FAIL'}] {tag}: {txt}")

    # (1) Jacobian vs central differences
    k = 400
    P = np.stack([rng.uniform(50, 2000, k), rng.uniform(90, 150, k),
                  rng.uniform(0.8, 5.0, k), rng.uniform(12, 40, k)], axis=1)
    xt = DT * (np.arange(LW)[None, :] + np.ceil(P[:, 1:2] / DT) - LEAD)
    _, J = _model_jac(P, xt)
    worst = 0.0
    for j in range(4):
        h = 1e-6 * np.maximum(np.abs(P[:, j]), 1.0)
        Pp, Pm = P.copy(), P.copy()
        Pp[:, j] += h
        Pm[:, j] -= h
        fd = (_model_jac(Pp, xt)[0] - _model_jac(Pm, xt)[0]) / (2 * h)[:, None]
        scale = np.abs(J[:, :, j]).max(1, keepdims=True)
        worst = max(worst, float((np.abs(fd - J[:, :, j]) / scale).max()))
    report("jacobian", worst < 1e-5, f"max rel err vs central FD = {worst:.2e} (<1e-5)")

    wide_r, wide_f = sr * 100.0, sf * 100.0
    fb = FTF(bank, sr, sf, mode="bounded")
    fp = FTF(bank, sr, sf, mode="prior")
    fpw = FTF(bank, wide_r, wide_f, mode="prior")

    # (2) noise-free, true tau = template
    n = 300
    mi = rng.integers(0, nm, n)
    amp = np.exp(rng.uniform(np.log(100), np.log(2000), n))
    t0 = rng.uniform(80, 160, n)
    tr_t, tf_t = bank.tau_r[mi], bank.tau_f[mi]
    Y = _wave(t0, amp, tr_t, tf_t, ped=200.0)
    rows, m_, t0s, as_, sg, ped, _ = _seeds(rec, Y, mi)
    yf = _wave(t0, amp, tr_t, tf_t)[rows]
    for tag, fit in (("bounded", fb), ("prior", fp)):
        fr = fit.fit(yf, m_, t0s, as_, sg)
        eA = np.abs(fr.amp / amp[rows] - 1).max()
        eT = np.abs(fr.t0 - t0[rows]).max()
        er = np.abs(fr.tau_r / tr_t[rows] - 1).max()
        ef = np.abs(fr.tau_f / tf_t[rows] - 1).max()
        ed = np.abs(fr.I_int / bank.iint(m_, fr.t0) - 1)
        report(f"noise-free template tau, {tag}",
               rows.size == n and (fr.status == 0).all() and max(eA, er, ef) < 1e-4
               and eT < 1e-3 and fr.chi2dof.max() < 1e-4 and np.median(ed) < 1e-4,
               f"N={rows.size} max rel A {eA:.1e} tau_r {er:.1e} tau_f {ef:.1e}; "
               f"max|dt0| {eT:.1e} ns; max chi2/dof {fr.chi2dof.max():.1e}; "
               f"I_int vs bank.iint median {np.median(ed):.1e} 90% "
               f"{np.percentile(ed, 90):.1e} max {ed.max():.1e} (256-pt phase table, "
               f"I_int steps in phase); iter mean {fr.niter.mean():.1f} max {fr.niter.max()}")

    # (3) noise-free, true tau shifted inside the box
    n = 300
    mi = np.zeros(n, dtype=np.int64)
    trs, tfs = np.full(n, 2.0), np.full(n, 26.5)
    amp = np.exp(rng.uniform(np.log(100), np.log(2000), n))
    t0 = rng.uniform(80, 160, n)
    Y = _wave(t0, amp, trs, tfs, ped=200.0)
    rows, m_, t0s, as_, sg, ped, _ = _seeds(rec, Y, mi)
    yf = _wave(t0, amp, trs, tfs)[rows]
    Itrue = mf_pulse.unit_integral(trs[rows], tfs[rows], tmax_of(trs[rows], tfs[rows]),
                                   np.mod(-t0[rows], DT), 2, 0.1, 2)
    for tag, fit in (("bounded", fb), ("prior(wide)", fpw)):
        fr = fit.fit(yf, m_, t0s, as_, sg)
        e = max(np.abs(fr.tau_r / 2.0 - 1).max(), np.abs(fr.tau_f / 26.5 - 1).max(),
                np.abs(fr.amp / amp[rows] - 1).max())
        eE = np.abs(fit.energy_units(fr, 1.3) / (1.3 * amp[rows] * Itrue) - 1).max()
        report(f"noise-free shifted tau, {tag}",
               (fr.status == 0).all() and e < 1e-3 and eE < 1e-4,
               f"N={rows.size} max rel err (A, tau) {e:.1e}; energy/kappa vs "
               f"A_true*I_true {eE:.1e}; iter mean {fr.niter.mean():.1f}")
    Efix = as_ * bank.iint(m_, t0s) / (amp[rows] * Itrue) - 1.0
    print(f"        fixed-template K=1 amp*iint vs truth: mean {Efix.mean()*100:+.2f}% "
          f"min {Efix.min()*100:+.2f}% max {Efix.max()*100:+.2f}%  (tau 2.0/26.5 vs "
          f"template {bank.tau_r[0]:.3f}/{bank.tau_f[0]:.2f})")

    # (4) noisy pulses, pulls (floor 0 so v is the true noise)
    n = 4000
    mi = rng.integers(0, nm, n)
    amp = np.exp(rng.uniform(np.log(30), np.log(2000), n))
    t0 = rng.uniform(80, 160, n)
    trt = bank.tau_r[mi] * rng.uniform(0.85, 1.15, n)
    tft = bank.tau_f[mi] * rng.uniform(0.9, 1.1, n)
    Yn = _wave(t0, amp, trt, tft, noise=1.5, rng=rng, ped=200.0)
    rows, m_, t0s, as_, sg, ped, chi_mf = _seeds(rec, Yn, mi)
    yn = Yn[rows] - ped[:, None]
    big = amp[rows] > 300
    bb = FTF(bank, sr, sf, mode="bounded", nsig=10.0, floor=0.0)
    pp = FTF(bank, wide_r, wide_f, mode="prior", floor=0.0)
    yx = Yn[rows] - 200.0                       # exact pedestal
    sgt = np.full(rows.size, 1.5)               # true noise
    for tag, fit in (("bounded(10 sig box)", bb), ("prior(wide)", pp)):
        for lab, yin, sin in (("exact ped+sigma", yx, sgt), ("est. ped+sigma", yn, sg)):
            fr = fit.fit(yin, m_, t0s, as_, sin)
            okb = big & (fr.status == 0)
            pr = (fr.tau_r - trt[rows])[okb] / fr.sig_tau_r[okb]
            pf = (fr.tau_f - tft[rows])[okb] / fr.sig_tau_f[okb]
            pa = (fr.amp - amp[rows])[okb] / fr.sig_amp[okb]
            txt = (f"N(A>300)={int(okb.sum())}/{int(big.sum())} SD pull tau_r {pr.std():.3f} "
                   f"tau_f {pf.std():.3f} (A {pa.std():.3f}); mean {pr.mean():+.3f} "
                   f"{pf.mean():+.3f}; frac|pull|>4 {np.mean(np.abs(pr) > 4)*100:.2f}% "
                   f"{np.mean(np.abs(pf) > 4)*100:.2f}%; status!=0 {int((fr.status != 0).sum())}; "
                   f"chi2/dof mean {fr.chi2dof[okb].mean():.3f}")
            if lab.startswith("exact"):
                report(f"noisy pulls, {tag}, {lab}",
                       0.85 <= pr.std() <= 1.15 and 0.85 <= pf.std() <= 1.15, txt)
            else:
                print(f"[info] noisy pulls, {tag}, {lab}: {txt}")
    # chi2dof_seed equals the mf K=1 chi2 (floor 0.03)
    fr = fb.fit(yn, m_, t0s, as_, sg)
    dchi = np.abs(fr.chi2dof_seed - chi_mf).max()
    report("chi2dof_seed == mf p_chi2", dchi < 1e-3, f"max abs diff {dchi:.2e}")

    # (5) truth outside the box (noise-free): bounded arm sits on the bound
    n = 100
    mi = np.zeros(n, dtype=np.int64)
    amp = np.exp(rng.uniform(np.log(200), np.log(2000), n))
    t0 = rng.uniform(80, 160, n)
    for tag, trs_, tfs_, key, which in (
            ("tau_r high", 5.0, 24.41, "at_hi_r", 1), ("tau_f low", 2.499, 18.5, "at_lo_f", 2)):
        trs, tfs = np.full(n, trs_), np.full(n, tfs_)
        Y = _wave(t0, amp, trs, tfs, ped=200.0)
        rows, m_, t0s, as_, sg, ped, _ = _seeds(rec, Y, mi)
        fr = fb.fit(_wave(t0, amp, trs, tfs)[rows], m_, t0s, as_, sg)
        lo_r, hi_r, lo_f, hi_f = fb._box(m_)
        bound = hi_r if which == 1 else lo_f
        val = fr.tau_r if which == 1 else fr.tau_f
        good = fr.status == 0
        flag = getattr(fr, key)
        report(f"truth outside box, {tag}",
               good.mean() > 0.95 and flag[good].all()
               and np.abs(val[good] - bound[good]).max() < 1e-6,
               f"N={rows.size} status0 {int(good.sum())}; {key} on all of them "
               f"{bool(flag[good].all())}; max|tau-bound| {np.abs(val[good]-bound[good]).max():.1e}; "
               f"status counts {np.bincount(fr.status, minlength=4).tolist()}")

    # (6) failure path
    n = 40
    mi = np.zeros(n, dtype=np.int64)
    amp = np.full(n, 400.0)
    t0 = rng.uniform(80, 160, n)
    Y = _wave(t0, amp, np.full(n, 2.0), np.full(n, 26.5), noise=1.5, rng=rng, ped=200.0)
    rows, m_, t0s, as_, sg, ped, _ = _seeds(rec, Y, mi)
    yn = Y[rows] - ped[:, None]
    f1 = FTF(bank, sr, sf, mode="bounded", max_iter=1)
    fr = f1.fit(yn, m_, t0s, as_, sg)
    bad = fr.status == 1
    seedok = (np.allclose(fr.tau_r[bad], bank.tau_r[m_][bad]) and
              np.allclose(fr.tau_f[bad], bank.tau_f[m_][bad]) and
              np.allclose(fr.t0[bad], t0s[bad]) and
              np.allclose(fr.I_int[bad], bank.iint(m_, t0s)[bad]) and
              np.allclose(fr.amp[bad], as_[bad], rtol=1e-4))
    z = np.zeros((2, 100))                      # zero waveform and a non-finite one
    z[1, 50] = np.nan
    frz = f1.fit(z, np.zeros(2, dtype=np.int64), np.array([100.0, 100.0]),
                 np.array([0.0, 50.0]), np.array([1.0, 1.0]))
    report("failure path", bad.sum() > 0 and seedok and (frz.status == 2).all()
           and np.isfinite(frz.I_int).all(),
           f"max_iter=1: {int(bad.sum())}/{rows.size} status 1 with seed outputs "
           f"({seedok}); zero/NaN waveform status {frz.status.tolist()}")

    # prior arm safety box: truth tau_f = 100 ns > 80 ns -> status 3 (direct seeds)
    n = 20
    t0 = rng.uniform(80, 160, n)
    amp = np.full(n, 1000.0)
    yb = _wave(t0, amp, np.full(n, 2.5), np.full(n, 100.0))
    fr = fpw.fit(yb, np.zeros(n, dtype=np.int64), t0, amp, np.full(n, 0.3))
    report("prior safety box", (fr.status != 0).all() and (fr.status == 3).any()
           and np.isfinite(fr.I_int).all() and np.allclose(fr.tau_f, bank.tau_f[0]),
           f"truth tau_f=100 ns: status counts {np.bincount(fr.status, minlength=4).tolist()}"
           f" (none ok; 3 = on the box, 1 = still creeping at max_iter; seed outputs)")

    # sigma table (real table, if present)
    tp = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output",
                      "pulse_templates_h50.dat")
    if os.path.isfile(tp):
        try:
            names = [l.split()[0] for l in open(tp) if l.strip() and not l.startswith("#")]
            pwo = np.array([True] * len(names) + [False])
            r_, f_, info = load_tau_sigma(tp, "25320", names + ["GX"], pwo)
            report("load_tau_sigma", np.isfinite(r_).all() and np.isfinite(f_).all()
                   and (r_ > 0).all() and r_[-1] == info["median_sig_r"], f"{info}")
        except (OSError, ValueError) as e:
            print(f"[info] load_tau_sigma skipped: {e}")

    # timing: 4000 noisy pulses, template tau, A 30-2000
    n = 4000
    mi4 = rng.integers(0, nm, n)
    Yn = _wave(rng.uniform(80, 160, n), np.exp(rng.uniform(np.log(30), np.log(2000), n)),
               bank.tau_r[mi4], bank.tau_f[mi4], noise=1.5, rng=rng, ped=200.0)
    rows, m_, t0s, as_, sg, ped, _ = _seeds(rec, Yn, mi4)
    yn = Yn[rows] - ped[:, None]
    for tag, fit in (("bounded", fb), ("prior", fp)):
        fit.fit(yn, m_, t0s, as_, sg)
        best = 1e9
        for _ in range(3):
            t = time.perf_counter()
            fr = fit.fit(yn, m_, t0s, as_, sg)
            best = min(best, time.perf_counter() - t)
        it = fr.niter[fr.niter > 0]
        print(f"[time] {tag}: n={rows.size} {best*1e6/rows.size:.1f} us/pulse "
              f"(full fit incl. I_int), mean iter {it.mean():.1f} max {it.max()}, "
              f"status counts {np.bincount(fr.status, minlength=4).tolist()}")
    t = time.perf_counter()
    fb._iint(fr.tau_r, fr.tau_f, fr.t0)
    print(f"[time] I_int alone: {(time.perf_counter()-t)*1e6/rows.size:.1f} us/pulse")

    return ok_all


def _invariant_test():
    """Tau pinned (bounded arm, nsig = 1e-6), uniform weights, production floor:
    the fitted energy A * I_int must equal the mf one, a_mf * I_int(t0_mf), where
    the mf (A, t0) is the same LS solution.  The fit moves t0 (up to ~2 ns on noisy
    low-amplitude pulses) and A with it, along a valley of constant energy, so the
    reference is a_mf times the exact unit integral at the mf t0 (not the 256-point
    table) and rows where I_int jumps with the phase (|I_int(t0_fit)/I_int(t0_mf) - 1|
    >= 1e-4) are dropped.  seedmodel weights are printed as info only."""
    rng = np.random.default_rng(424242)
    taus = [(2.499, 24.41), (2.0, 20.0), (3.2, 28.0)]
    bank = mf_pulse.TemplateBank([t[0] for t in taus], [t[1] for t in taus], [True] * 3,
                                 smooth_order=2, int_tail_ratio=0.1, tail_break_n=2)
    rec = mf_pulse.MFReconstructor(bank)
    sr, sf = np.array([0.6, 0.7, 0.8]), np.array([1.5, 1.2, 2.0])
    n = 1500
    mi = rng.integers(0, 3, n)
    amp = np.exp(rng.uniform(np.log(5), np.log(3000), n))
    t0 = rng.uniform(80, 160, n)
    Y = _wave(t0, amp, bank.tau_r[mi], bank.tau_f[mi], noise=1.5, rng=rng, ped=200.0)
    rows, m_, t0s, as_, sg, ped, _ = _seeds(rec, Y, mi)
    yn = Y[rows] - ped[:, None]
    i0 = mf_pulse.unit_integral(bank.tau_r[m_], bank.tau_f[m_], bank.tmx[m_],
                                np.mod(-t0s, DT), 2, 0.1, 2)
    out = {}
    for w in ("uniform", "seedmodel"):
        fr = FreeTauFitter(bank, sr, sf, mode="bounded", nsig=1e-6, weight=w) \
            .fit(yn, m_, t0s, as_, sg)
        away = (fr.status == 0) & (np.abs(fr.I_int / i0 - 1) < 1e-4)
        r = np.abs(fr.amp * fr.I_int / (as_ * i0) - 1)[away]
        out[w] = (r.max(), np.percentile(r, 99), int(away.sum()), int((fr.status == 0).sum()),
                  np.abs(fr.t0 - t0s)[fr.status == 0].max(),
                  np.abs(fr.amp / as_ - 1)[fr.status == 0].max())
    mx, p99, na, nok, dt, dA = out["uniform"]
    passed = mx < 2e-4 and p99 < 3e-5
    print(f"[{'PASS' if passed else 'FAIL'}] energy invariant, tau pinned, uniform "
          f"weights: N={rows.size} status0 {nok} away-from-jumps {na}; max|A*I_int/"
          f"(a_mf*I_int(t0_mf))-1| = {mx:.2e} (<2e-4), 99% {p99:.1e} (<3e-5); fit moves "
          f"t0 by up to {dt:.2f} ns and A by up to {dA*100:.0f}% (valley of constant E)")
    mx, p99, na, nok, dt, dA = out["seedmodel"]
    print(f"[info] same, weight=seedmodel (not gated): away-from-jumps {na}; max {mx:.2e}, "
          f"99% {p99:.1e}")
    return passed


def self_test():
    ok = bool(_invariant_test())
    for w in ("uniform", "seedmodel"):
        print(f"--- weight={w}")
        ok &= bool(_self_test_w(w))
    print("SELF-TEST", "PASS" if ok else "FAIL")
    return ok


if __name__ == "__main__":
    if "--self-test" in sys.argv:
        sys.exit(0 if self_test() else 1)
    print(__doc__)
