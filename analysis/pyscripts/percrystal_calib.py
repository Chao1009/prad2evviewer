#!/usr/bin/env python3
"""percrystal_calib.py — per-crystal (member-level) HyCal gain solver.

Implements plan_20261023.md section 4.2.  One gain g_i per PbWO4 module that
appears as a cluster member, fitted so that the gain-corrected cluster energy
sum_i g_i E_ki of elastic events matches the kinematic E'_k:

    chi2(g) = sum_k w (sum_i g_i E_ki - E'_k)^2 + sum_i lambda_i (g_i - g0)^2

k runs over the pool events that lie inside the elastic window of their seed
module (windows from the benchmark's fit_elastic_window_gains, re-evaluated on
the current E_cal each iteration), i over the cluster members.

Deviation from the plan text (lead decision): the calibration files carry a
+8-10 % global offset (mu_global/E' ~ 1.09), so the prior centre is
g0 = eprime_mean / mu_global of the FIRST window pass (g = 1), not 1.  The
prior strength is lambda_i = D_min = 30 (0.64 E'_mean)^2 w with
w = 1/sigma_glob^2 of the first pass (constant per run, also used as the
event weight).  prior_frac_i = lambda_i / (D_ii + lambda_i).  A module that is
never a member of a selected event keeps g = g0 (also when the gains are
applied to any pool event).

Iteration: start g = 1 (or an initial gain table); each pass evaluates
E_cal = sum_i g_i E_ki, re-runs the window machinery on it (first N_REFIT = 3
passes), solves, and stops when |dg_i| < max(1e-3, 0.5 sigma_g,i) for every
free module (prior_frac < 0.5; the sigma term is a statistical floor for the
low-statistics modules) or after 6 passes.  Deviation: the per-module window fits
(Crystal-Ball, window edges) are discontinuous in E_cal; on a 60k-event split
a few low-statistics modules flip between two window states and g jumps by
about 1 sigma_g from pass to pass although the bulk is converged to 1e-5
(median |dg|), so a full refit never reaches 1e-3.  After the last refit pass
the window edges are therefore frozen (in E_cal units) and only the event
membership is re-evaluated with the new E_cal, which is a contraction.
`rows` of the result are those of the last refit pass.
fixed_windows=True (fit --fixed-windows, added after the 2026-10-23 scans):
the event set is selected once with the g = 1 windows in pre-refinement
energy, as in the seed-indexed refinement, so the result does not depend on
the start (P1-6 max |dg| 0 vs 4e-3); on 025320.00050 the gains differ from
the default by median 8.5e-4, RMS 2.6e-3, sigma/E +0.004 pp (Gauss).

Module-level import is numpy-only (no prad2py).  The offline CLI imports the
benchmark lazily to obtain fit_elastic_window_gains.

Usage
  percrystal_calib.py --self-test
  percrystal_calib.py fit --dumps GLOB [GLOB ...] --splits even|odd|all
                          --out DIR [--exclude-rowcol 15 20 15 20]
  percrystal_calib.py stability --gains GLOB [GLOB ...] [--path legacy]
  percrystal_calib.py compare --gains TSV [--window TSV] [--path legacy]
"""

from __future__ import annotations

import argparse
import functools
import glob
import hashlib
import json
import math
import os
import re
import sys

import numpy as np

SCHEMA = 1
SEED_FRAC = 0.64            # mean seed energy fraction (plan section 4.2)
PRIOR_MIN_EVENTS = 30       # events defining D_min (= --min-per-module)
FREE_FRAC = 0.5             # prior_frac below this = "free" module
DG_TOL = 1e-3
FIXED_WINDOWS = False        # True: windows and event set of the g = 1 (seed-indexed) pass
TOL_SIGMA = 0.5             # statistical floor of the per-module dg threshold
MAX_ITER = 6
N_REFIT = 3                 # iterations that re-run the window fits
GAIN_COLS = ("path", "module", "index", "row", "col", "g", "sigma_g", "n_seed",
             "n_member", "D_ii", "prior_frac", "g_seedindexed", "used")


def file_sha1(path):
    h = hashlib.sha1()
    with open(path, "rb") as f:
        for blk in iter(lambda: f.read(1 << 20), b""):
            h.update(blk)
    return h.hexdigest()


# ---------------------------------------------------------------------------
# Solver
# ---------------------------------------------------------------------------
def _normal_equations(ptr, midx_c, mE, ev_tgt, n_mod, chunk_elems=2_000_000):
    """A^T A (without weights) and A^T t for events given as CSR members.
    ptr: int64 (n_ev+1), midx_c: compact module index per member, mE: member
    energies, ev_tgt: target per event.  Built per member-count group with
    bincount over the flattened outer products (no Python loop over events)."""
    n_ev = len(ptr) - 1
    counts = np.diff(ptr)
    N = np.zeros(n_mod * n_mod)
    for c in np.unique(counts):
        evs = np.nonzero(counts == c)[0]
        step = max(1, chunk_elems // int(c * c))
        for a in range(0, len(evs), step):
            ev = evs[a:a + step]
            pos = ptr[ev][:, None] + np.arange(c)[None, :]
            m2 = midx_c[pos]
            e2 = mE[pos]
            flat = (m2[:, :, None] * n_mod + m2[:, None, :]).ravel()
            val = (e2[:, :, None] * e2[:, None, :]).ravel()
            N += np.bincount(flat, weights=val, minlength=n_mod * n_mod)
    ev_of = np.repeat(np.arange(n_ev), counts)
    b = np.bincount(midx_c, weights=mE * ev_tgt[ev_of], minlength=n_mod)
    return N.reshape(n_mod, n_mod), b


def solve(events, window_fn, *, seed_to_idx=None, init=None,
          max_iter=MAX_ITER, n_refit=N_REFIT, tol=DG_TOL,
          tol_sigma=TOL_SIGMA, free_frac=FREE_FRAC, fixed_windows=FIXED_WINDOWS,
          prior_min_events=PRIOR_MIN_EVENTS, seed_frac=SEED_FRAC,
          verbose=False):
    """Iterate window selection and the regularised least-squares solve.

    events      list of event tuples (seed_key, x, y, row, col, seed_E,
                cluster_E, f_eff, members_idx, members_E, ...); the pool
                (already outside the rectangle etc.); members_idx are module
                indices.
    window_fn   callable(list of 8-tuples with cluster_E = current E_cal) ->
                (gains, raw, corr, n_mod, rows, sigma_glob, eprime_mean,
                mu_global), i.e. the benchmark's fit_elastic_window_gains
                (functools.partial of it).  Only rows (dicts with module,
                win_lo, win_hi, used, target[, gain]), sigma_glob,
                eprime_mean, mu_global are used.
    seed_to_idx dict seed_key -> module index (None: the key is the index).
    init        None (start g = 1) or dict module index -> g (modules not
                listed start at 1); g0 and w always come from the g = 1 pass.

    Returns a dict (None if the window pass fails), see the keys at the end.
    """
    n = len(events)
    if n == 0:
        return None
    nan = float('nan')
    keys = np.array([e[0] for e in events])
    tuples8 = [e[:8] for e in events]
    ce0 = np.array([e[6] for e in events], dtype=float)
    mi_list = [np.asarray(e[8]) for e in events]
    counts = np.fromiter((len(m) for m in mi_list), np.int64, n)
    ptr = np.zeros(n + 1, np.int64)
    np.cumsum(counts, out=ptr[1:])
    midx = np.concatenate(mi_list).astype(np.int64)
    mE = np.concatenate([np.asarray(e[9], dtype=float) for e in events])
    evid = np.repeat(np.arange(n), counts)
    nmax = int(midx.max()) + 1
    if seed_to_idx is None:
        seed_idx = keys.astype(np.int64)
    else:
        seed_idx = np.array([seed_to_idx[k] for k in keys], dtype=np.int64)
        nmax = max(nmax, int(seed_idx.max()) + 1)
    ukeys, kinv = np.unique(keys, return_inverse=True)

    def ecal_of(gfull):
        return np.bincount(evid, weights=mE * gfull[midx], minlength=n)

    def window_pass(ecal):
        evs = [t[:6] + (c,) + t[7:] for t, c in zip(tuples8, ecal)]
        return window_fn(evs)

    # Pass 0 at g = 1: prior centre, weight, seed-indexed reference.
    out0 = window_pass(ce0)
    rows0, sg0, ep_mean, mu_g0 = out0[4], out0[5], out0[6], out0[7]
    if not (math.isfinite(sg0) and math.isfinite(mu_g0) and sg0 > 0):
        return None
    g0 = ep_mean / mu_g0
    w = 1.0 / (sg0 * sg0)
    lam = prior_min_events * (seed_frac * ep_mean) ** 2 * w
    gsi = {r['module']: r.get('gain', nan) for r in rows0 if r['used']}

    gprev = np.ones(nmax)
    if init is not None:
        for i, v in init.items():
            if i < nmax:
                gprev[i] = v
    history, converged, it = [], False, 0
    res = None
    kpos = {k: i for i, k in enumerate(ukeys.tolist())}
    while True:
        it += 1
        ecal = ecal_of(gprev)
        if it <= (1 if fixed_windows else n_refit):
            out = out0 if (it == 1 and (init is None or fixed_windows)) else window_pass(ecal)
            rows, sg, ep_it, mu_g = out[4], out[5], out[6], out[7]
            if not (math.isfinite(sg) and math.isfinite(mu_g)):
                return None
            lo_k = np.full(len(ukeys), np.nan)
            hi_k = np.full(len(ukeys), np.nan)
            tg_k = np.full(len(ukeys), np.nan)
            used_k = np.zeros(len(ukeys), bool)
            for r in rows:
                i = kpos.get(r['module'])
                if i is None:
                    continue
                tg_k[i] = r['target']
                if r['used']:
                    used_k[i] = True
                    lo_k[i], hi_k[i] = r['win_lo'], r['win_hi']
        # else: windows frozen at the last refit pass (energy edges in E_cal
        # units); only the membership is re-evaluated with the new E_cal.
        esel = ce0 if fixed_windows else ecal
        sel = used_k[kinv] & (esel >= lo_k[kinv]) & (esel <= hi_k[kinv])
        n_sel = int(sel.sum())
        if n_sel == 0:
            return None
        msel = sel[evid]
        uniq, inv_m = np.unique(midx[msel], return_inverse=True)
        M = len(uniq)
        sel_ids = np.nonzero(sel)[0]
        ptr_s = np.zeros(n_sel + 1, np.int64)
        np.cumsum(counts[sel_ids], out=ptr_s[1:])
        E_s = mE[msel]
        tgt_s = tg_k[kinv][sel_ids]
        Nm, b = _normal_equations(ptr_s, inv_m, E_s, tgt_s, M)
        Nm *= w
        b *= w
        D = np.diag(Nm).copy()
        A = Nm + lam * np.eye(M)
        gsol = np.linalg.solve(A, b + lam * g0)
        # formal errors at w = 1/sg0^2, rescaled to the current-pass width sg
        # (sg0 is the uncalibrated g = 1 width, ~1.75x the calibrated one)
        sigma_g = np.sqrt(np.diag(np.linalg.inv(A))) * (sg / sg0)
        pf = lam / (D + lam)
        ev_new = np.bincount(np.repeat(np.arange(n_sel), np.diff(ptr_s)),
                             weights=E_s * gsol[inv_m], minlength=n_sel)
        chi2 = float(w * np.sum((ev_new - tgt_s) ** 2)
                     + lam * np.sum((gsol - g0) ** 2))
        gnew = np.full(nmax, g0)
        gnew[uniq] = gsol
        free = pf < free_frac
        dgv = np.abs(gnew[uniq][free] - gprev[uniq][free])
        dg = float(dgv.max()) if free.any() else nan
        # events crossing a window edge move a low-statistics module by up to
        # ~0.5 sigma_g from pass to pass: the threshold has a statistical floor
        ok_conv = bool(free.any() and (dgv < np.maximum(
            tol, tol_sigma * sigma_g[free])).all())
        history.append(dg)
        if verbose:
            print(f"    [percrystal] iter {it}: n_sel {n_sel}  modules {M}  "
                  f"free {int(free.sum())}  max|dg| {dg:.5f}", flush=True)
        seed_sel = seed_idx[sel_ids]
        pos = np.searchsorted(uniq, seed_sel)
        okp = (pos < M) & (uniq[np.minimum(pos, M - 1)] == seed_sel)
        n_seed = np.bincount(pos[okp], minlength=M)
        n_member = np.bincount(inv_m, minlength=M)
        gprev = gnew
        res = dict(uniq=uniq, g=gsol, sigma_g=sigma_g, D=D, prior_frac=pf,
                   n_seed=n_seed, n_member=n_member, n_sel=n_sel, chi2=chi2,
                   rows=rows, sigma_glob=sg, eprime_mean=ep_it,
                   mu_global=mu_g, sel=sel, used_keys=ukeys[used_k])
        if ok_conv:
            converged = True
            break
        if it >= max_iter:
            break

    uniq, pf = res['uniq'], res['prior_frac']
    seed_of_mod = {int(i): k for i, k in
                   zip(seed_idx.tolist(), keys.tolist())}   # idx -> seed key
    key_to_idx = {k: i for i, k in seed_of_mod.items()}
    used_idx = {key_to_idx[k] for k in res['used_keys'].tolist()}
    n_free = int((pf < free_frac).sum())
    gseed = np.array([gsi.get(seed_of_mod.get(int(i)), nan) for i in uniq])
    res.update(
        gfull=gprev, g0=g0, w=w, lam=lam, iterations=it, converged=converged,
        history=history, n_free=n_free, n_prior=len(uniq) - n_free,
        sigma_glob_first=sg0, mu_global_first=mu_g0,
        eprime_mean_first=ep_mean, ecal=ecal_of(gprev),
        g_seedindexed=gseed, used_idx=used_idx, seed_idx=seed_idx,
        n_events=n)
    return res


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
def gain_rows(res, modinfo):
    """Per-module table rows (dicts) of a solve() result.  modinfo: module
    index -> (name, row, col)."""
    out = []
    for j, i in enumerate(res['uniq'].tolist()):
        name, row, col = modinfo.get(i, (str(i), -1, -1))
        out.append(dict(module=name, index=i, row=row, col=col,
                        g=res['g'][j], sigma_g=res['sigma_g'][j],
                        n_seed=int(res['n_seed'][j]),
                        n_member=int(res['n_member'][j]),
                        D_ii=res['D'][j], prior_frac=res['prior_frac'][j],
                        g_seedindexed=res['g_seedindexed'][j],
                        used=int(i in res['used_idx'])))
    return out


def write_gain_table(path, results, modinfo, extra_header=()):
    """results: {path_label: solve() result}.  One file, `path` column."""
    with open(path, "w") as f:
        f.write("# percrystal per-module gain table (percrystal_calib.py, "
                f"schema {SCHEMA})\n")
        for line in extra_header:
            f.write(f"# {line}\n")
        for lab, res in results.items():
            if res is None:
                f.write(f"# path={lab} failed\n")
                continue
            f.write(f"# path={lab} g0={res['g0']:.6f} "
                    f"iterations={res['iterations']} "
                    f"converged={int(res['converged'])} "
                    f"n_free={res['n_free']} n_prior={res['n_prior']} "
                    f"sigma_glob={res['sigma_glob']:.3f} "
                    f"eprime_mean={res['eprime_mean']:.3f} "
                    f"mu_global={res['mu_global']:.3f} "
                    f"sigma_glob_first={res['sigma_glob_first']:.3f} "
                    f"mu_global_first={res['mu_global_first']:.3f} "
                    f"n_sel_events={res['n_sel']}\n")
        f.write("\t".join(GAIN_COLS) + "\n")
        for lab, res in results.items():
            if res is None:
                continue
            for r in gain_rows(res, modinfo):
                f.write(f"{lab}\t{r['module']}\t{r['index']}\t{r['row']}\t"
                        f"{r['col']}\t{r['g']:.6f}\t{r['sigma_g']:.6f}\t"
                        f"{r['n_seed']}\t{r['n_member']}\t{r['D_ii']:.6g}\t"
                        f"{r['prior_frac']:.4f}\t{r['g_seedindexed']:.6f}\t"
                        f"{r['used']}\n")


def _read_tsv(path):
    """-> (header comment lines, column names, list of row dicts of str)."""
    hdr, cols, rows = [], None, []
    with open(path) as f:
        for line in f:
            line = line.rstrip("\n")
            if not line.strip():
                continue
            if line.startswith("#"):
                hdr.append(line)
            elif cols is None:
                cols = line.split("\t")
            else:
                rows.append(dict(zip(cols, line.split("\t"))))
    return hdr, cols, rows


def _header_kv(hdr, path_label):
    """key=value pairs of the '# path=<label> ...' header line."""
    for line in hdr:
        toks = line[1:].split()
        kv = dict(t.split("=", 1) for t in toks if "=" in t)
        if kv.get("path") == path_label:
            return kv
    return {}


def read_gain_table(path, path_label=None):
    """-> (dict module index -> g, g0).  A table with several paths yields
    `path_label` (default 'legacy' if present, else the first path).  g0 comes
    from the '# path=<label> g0=..' header line."""
    hdr, cols, rows = _read_tsv(path)
    if cols is None or "g" not in cols:
        raise ValueError(f"{path}: not a percrystal gain table")
    if "path" in cols:
        labs = []
        for r in rows:
            if r["path"] not in labs:
                labs.append(r["path"])
        lab = path_label or ("legacy" if "legacy" in labs else labs[0])
        rows = [r for r in rows if r["path"] == lab]
    else:
        lab = path_label
    kv = _header_kv(hdr, lab)
    if "g0" not in kv:
        raise ValueError(f"{path}: no g0 in the header for path {lab!r}")
    gains = {int(r["index"]): float(r["g"]) for r in rows
             if math.isfinite(float(r["g"]))}
    return gains, float(kv["g0"])


# ---------------------------------------------------------------------------
# Self-test (P1-2): synthetic 300-module detector, no prad2py
# ---------------------------------------------------------------------------
def _synthetic(rng, n_clusters=50000, nx=20, ny=15, pitch=20.77, g_mean=0.91,
               g_rms=0.02, smear=0.02):
    n_mod = nx * ny
    gx = (np.arange(n_mod) % nx - (nx - 1) / 2.0) * pitch
    gy = (np.arange(n_mod) // nx - (ny - 1) / 2.0) * pitch
    row = np.arange(n_mod) // nx
    col = np.arange(n_mod) % nx
    rad = np.hypot(gx, gy)
    g_true = g_mean * (1.0 + rng.normal(0.0, g_rms, n_mod))
    eprime = 2200.0 - 2.0 * rad               # simple angle model (MeV)
    # elastic yield falls steeply with radius
    wgt = np.exp(-rad / 60.0)
    wgt /= wgt.sum()
    seeds = rng.choice(n_mod, n_clusters, p=wgt)
    frac = np.clip(rng.normal(0.644, 0.124, n_clusters), 0.3, 0.95)
    etot = eprime[seeds] * (1.0 + smear * rng.standard_normal(n_clusters))

    def nb(m, r):
        out = []
        for dy in range(-r, r + 1):
            for dx in range(-r, r + 1):
                if (dx or dy) and 0 <= col[m] + dx < nx and 0 <= row[m] + dy < ny:
                    out.append(m + dy * nx + dx)
        return out
    events = []
    for k in range(n_clusters):
        s = int(seeds[k])
        use5 = rng.random() < 0.2
        nbs = nb(s, 2 if use5 else 1)
        idx = [s] + nbs
        wn = np.array([3.0 if (abs(col[m] - col[s]) + abs(row[m] - row[s]) == 1)
                       else (1.0 if max(abs(col[m] - col[s]),
                                        abs(row[m] - row[s])) == 1 else 0.25)
                       for m in nbs])
        wn = wn * rng.gamma(4.0, 0.25, len(wn))
        share = np.concatenate(([frac[k]], (1.0 - frac[k]) * wn / wn.sum()))
        e_true = share * etot[k]
        idx = np.array(idx, np.int32)
        e_meas = e_true / g_true[idx]
        keep = e_meas > 1.0
        idx, e_meas = idx[keep], e_meas[keep]
        events.append((s, gx[s], gy[s], int(row[s]), int(col[s]),
                       float(e_meas[0]), float(e_meas.sum()), 1.0,
                       idx, e_meas, gx[s], gy[s], k))
    return events, g_true, eprime, (gx, gy, row, col)


def _standin_window_fn(eprime_of, nsig=2.0, min_n=30):
    """Per-seed median +- nsig*sigma (sigma = pooled relative robust width);
    used if >= min_n in window.  Same return layout as the benchmark's."""
    def fn(events):
        sid = np.array([e[0] for e in events])
        ce = np.array([e[6] for e in events])
        ep = np.array([eprime_of[s] for s in sid])
        rel = ce / ep
        med_all = float(np.median(rel))
        sig_rel = 1.4826 * float(np.median(np.abs(rel - med_all)))
        ep_mean = float(np.mean([eprime_of[s] for s in np.unique(sid)]))
        rows = []
        for s in np.unique(sid):
            v = ce[sid == s]
            m = float(np.median(v))
            lo, hi = m - nsig * sig_rel * m, m + nsig * sig_rel * m
            nin = int(np.count_nonzero((v >= lo) & (v <= hi)))
            used = nin >= min_n
            rows.append(dict(module=int(s), win_lo=lo, win_hi=hi, used=int(used),
                             target=eprime_of[s], n_in_window=nin,
                             gain=(eprime_of[s] / float(np.median(
                                 v[(v >= lo) & (v <= hi)])) if used else float('nan'))))
        return ({}, None, None, 0, rows, sig_rel * med_all * ep_mean,
                ep_mean, med_all * ep_mean)
    return fn


def self_test(seed=7, verbose=True):
    rng = np.random.default_rng(seed)
    events, g_true, eprime, _ = _synthetic(rng)
    eprime_of = {i: float(v) for i, v in enumerate(eprime)}
    wfn = _standin_window_fn(eprime_of)
    import time
    t0 = time.time()
    res = solve(events, wfn, verbose=verbose)
    dt = time.time() - t0
    if res is None:
        print("[self-test] FAIL: solve returned None")
        return 1
    uniq, g, sg, pf = res['uniq'], res['g'], res['sigma_g'], res['prior_frac']
    D, lam, g0 = res['D'], res['lam'], res['g0']
    dmin = lam
    gt = g_true[uniq]
    big = D > 4.0 * dmin
    rms_big = float(np.sqrt(np.mean((g[big] / gt[big] - 1.0) ** 2)))
    mean_big = float(np.mean(g[big] / gt[big] - 1.0))
    # The prior shrinks every module by prior_frac towards g0 (~10 % for
    # D_ii ~ 4-10 D_min, i.e. 0.2 % RMS at a 2 % gain spread); the pass rule
    # compares with the lambda-weighted expectation g_exp.
    gexp_big = (pf[big] * g0 + (1.0 - pf[big]) * gt[big]) / gt[big] - 1.0
    rms_exp = float(np.sqrt(np.mean((g[big] / gt[big] - 1.0 - gexp_big) ** 2)))
    pri = pf > FREE_FRAC
    ok_prior = True
    z_rms = float('nan')
    if pri.any():
        gexp = pf[pri] * g0 + (1.0 - pf[pri]) * gt[pri]
        z = (g[pri] - gexp) / sg[pri]
        z_rms = float(np.sqrt(np.mean(z ** 2)))
        ok_prior = z_rms < 2.0
    # init independence at a tight tolerance: g = 1 start vs true-gain start
    resa = solve(events, wfn, tol=1e-5, tol_sigma=0.0, max_iter=12, n_refit=12)
    resb = solve(events, wfn, tol=1e-5, tol_sigma=0.0, max_iter=12, n_refit=12,
                 init={int(i): float(g_true[i]) for i in range(len(g_true))})
    dinit = float('nan')
    if resa is not None and resb is not None:
        dinit = float(np.max(np.abs(resa['g'] - resb['g'])
                             [resa['prior_frac'] < FREE_FRAC]))
    print(f"[self-test] {len(events)} clusters, {len(uniq)} modules solved in "
          f"{dt:.1f} s; iterations {res['iterations']} "
          f"converged {res['converged']}; n_free {res['n_free']} "
          f"n_prior {res['n_prior']}; g0 {g0:.4f} (true mean "
          f"{g_true.mean():.4f})")
    print(f"[self-test] D_ii > 4 D_min: {int(big.sum())} modules, "
          f"RMS(g/g_true - 1) {rms_big*100:.3f} % (raw, includes the prior "
          f"shrinkage), mean {mean_big*100:+.3f} %; RMS vs lambda-weighted "
          f"expectation {rms_exp*100:.3f} % (true spread "
          f"{np.std(g_true/g_true.mean())*100:.2f} %)")
    print(f"[self-test] prior-dominated: {int(pri.sum())} modules, "
          f"RMS z vs lambda-weighted expectation "
          f"{z_rms:.2f}")
    print(f"[self-test] init independence (tol 1e-5, g=1 vs true-gain "
          f"start): max|dg| {dinit:.2e} (reported, not gated: window "
          f"re-selection moves g by ~1e-3)")
    ok = (rms_exp < 0.003 and ok_prior and res['converged']
          and res['iterations'] <= 4)
    print("[self-test] " + ("PASS" if ok else "FAIL"))
    return 0 if ok else 1


# ---------------------------------------------------------------------------
# Offline CLI
# ---------------------------------------------------------------------------
def _expand(globs):
    out = []
    for g in globs:
        out.extend(sorted(glob.glob(g)))
    return sorted(dict.fromkeys(out))


def _bench():
    """Lazy import of the benchmark module (needs prad2py)."""
    import importlib
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    return importlib.import_module("benchmark_hycal_timing_test_20261021")


def load_dump(path):
    """-> (header dict, list of event tuples, modinfo).  Tuples: the 8-tuple
    of the benchmark with seed key = seed module index, then members_idx,
    members_E, x_cl, y_cl, evnum."""
    with np.load(path, allow_pickle=False) as z:
        d = {k: z[k] for k in z.files}      # NpzFile re-reads on every access
    hdr = json.loads(str(d['header']))
    ptr = d['ptr']
    midx, E = d['midx'], d['E'].astype(np.float64)
    cols = [d[k].tolist() for k in ('seed_idx', 'x_seed', 'y_seed', 'row',
                                    'col', 'E_seed', 'E_cluster', 'f_lms',
                                    'x_cl', 'y_cl', 'evnum')]
    ev = []
    for i, (si, xs, ys, ro, co, es, ec, fl, xc, yc, en) in enumerate(zip(*cols)):
        a, b = ptr[i], ptr[i + 1]
        ev.append((si, xs, ys, ro, co, es, ec, fl, midx[a:b], E[a:b], xc, yc,
                   en))
    modinfo = {int(i): (str(nm), int(r), int(c)) for i, nm, r, c in
               zip(d['mod_index'], d['mod_name'], d['mod_row'], d['mod_col'])}
    return hdr, ev, modinfo


def window_fn_from_header(hdr, bench, exclude_rowcol):
    h = hdr
    excl = exclude_rowcol if exclude_rowcol is not None else \
        (tuple(h['exclude_rowcol']) if h.get('exclude_rowcol') else None)
    return functools.partial(
        bench.fit_elastic_window_gains, ebeam=h['ebeam'], hycal_z=h['hycal_z'],
        m_target=h['m_target'], escale=h['kS'], kE=h['kE'],
        inner_r=h['inner_radius'], min_seed_E=h['min_seed_energy'],
        min_per_mod=h['min_per_module'], exclude_rowcol=excl,
        nsigma=h['window_nsigma'], bin_width=h['bin_width']), excl


def cmd_fit(args):
    bench = _bench()
    files = _expand(args.dumps)
    if not files:
        print("[ERROR] no dump matches", file=sys.stderr)
        return 2
    byrun = {}
    for f in files:
        d = np.load(f, allow_pickle=False)
        h = json.loads(str(d['header']))
        byrun.setdefault(int(h['run']), []).append((int(h['split']), f))
    os.makedirs(args.out, exist_ok=True)
    for run, lst in sorted(byrun.items()):
        lst.sort()
        # even/odd parity must match the harness (position in the run's full
        # split list): abort if a split dir of this run has no dump
        sdirs = {os.path.dirname(os.path.dirname(os.path.abspath(f)))
                 for _, f in lst}
        n_dirs = sum(len(glob.glob(os.path.join(d, "split_*"))) for d in sdirs)
        if args.splits != 'all' and n_dirs != len(lst):
            print(f"[ERROR] run {run}: {len(lst)} dumps but {n_dirs} split dirs; "
                  f"even/odd parity would not match the harness", file=sys.stderr)
            return 2
        if args.splits == 'even':
            sel = lst[0::2]
        elif args.splits == 'odd':
            sel = lst[1::2]
        else:
            sel = lst
        if not sel:
            print(f"[WARN] run {run}: no splits for --splits {args.splits}")
            continue
        events, modinfo, hdr0 = [], {}, None
        for sp, f in sel:
            hdr, ev, mi = load_dump(f)
            if hdr0 is None:
                hdr0 = hdr
            events.extend(ev)
            modinfo.update(mi)
        wfn, excl = window_fn_from_header(hdr0, bench, args.exclude_rowcol and
                                          tuple(args.exclude_rowcol))
        pool = bench._select_pool(events, inner_r=hdr0['inner_radius'],
                                  min_seed_E=hdr0['min_seed_energy'],
                                  exclude_rowcol=excl)
        print(f"[fit] run {run}: splits {[s for s, _ in sel]} "
              f"({args.splits}), {len(events)} dump events -> pool {len(pool)} "
              f"after rectangle {excl}", flush=True)
        res = solve(pool, wfn, verbose=True, max_iter=args.max_iter,
                    fixed_windows=args.fixed_windows)
        if res is None:
            print(f"[ERROR] run {run}: solver failed (window pass)")
            continue
        tag = args.splits
        base = os.path.join(args.out, f"gains_{run:06d}_{tag}")
        calsha = hdr0.get('calibration_sha1')
        write_gain_table(base + ".tsv", {'legacy': res}, modinfo,
                         extra_header=[
                             f"run={run} splits={','.join(str(s) for s, _ in sel)} "
                             f"half={tag} n_events={len(pool)} "
                             f"exclude_rowcol={excl} rowcol_origin=1 "
                             f"fixed_windows={int(args.fixed_windows)} "
                             f"max_iter={args.max_iter}",
                             f"calibration_file={hdr0.get('calibration_file')} "
                             f"sha1={calsha}"])
        prov = dict(schema=SCHEMA, run=run, half=tag,
                    splits=[s for s, _ in sel], dumps=[f for _, f in sel],
                    n_events=len(pool), exclude_rowcol=excl, rowcol_origin=1,
                    fixed_windows=bool(args.fixed_windows),
                    max_iter=args.max_iter, argv=sys.argv,
                    iterations=res['iterations'], converged=res['converged'],
                    history=res['history'], n_free=res['n_free'],
                    n_prior=res['n_prior'], g0=res['g0'],
                    sigma_glob_first=res['sigma_glob_first'],
                    n_sel_events=res['n_sel'], calibration_sha1=calsha,
                    dump_header=hdr0)
        with open(base + ".json", "w") as f:
            json.dump(prov, f, indent=1, default=str)
        g = res['g']
        free = res['prior_frac'] < FREE_FRAC
        print(f"[fit] run {run}: iterations {res['iterations']} converged "
              f"{res['converged']} n_free {res['n_free']} n_prior "
              f"{res['n_prior']} g0 {res['g0']:.4f}; g(free) mean "
              f"{g[free].mean():.4f} sd {g[free].std():.4f}; wrote {base}.tsv")
    return 0


_RUN_SPLIT = re.compile(r"scan_0*(\d+)[/_]split_0*(\d+)")


def _job_key(path):
    m = _RUN_SPLIT.search(path)
    return (int(m.group(1)), int(m.group(2))) if m else None


def cmd_stability(args):
    files = _expand(args.gains)
    jobs = {}
    for f in files:
        key = _job_key(os.path.abspath(f))
        if key is None:
            print(f"[WARN] cannot read run/split from {f}; skipped")
            continue
        hdr, cols, rows = _read_tsv(f)
        rows = [r for r in rows if r.get("path", args.path) == args.path]
        win = {}
        wpath = re.sub(r"_gains\.tsv$", "_window.tsv", f)
        if wpath != f and os.path.isfile(wpath):
            _, _, wr = _read_tsv(wpath)
            win = {r["module"]: float(r["gain"]) for r in wr
                   if r.get("path") == args.path and r.get("used") == "1"}
        d = {}
        for r in rows:
            d[int(r["index"])] = dict(
                g=float(r["g"]), pf=float(r["prior_frac"]),
                gsi=float(r["g_seedindexed"]), win=win.get(r["module"],
                                                           float('nan')))
        jobs[key] = d
    if not jobs:
        print("[ERROR] no usable gain tables", file=sys.stderr)
        return 2
    runs = sorted({k[0] for k in jobs})
    print(f"[stability] {len(jobs)} jobs, runs "
          + ", ".join(f"{r}({sum(1 for k in jobs if k[0] == r)})" for r in runs))
    common = set.intersection(*[set(d) for d in jobs.values()])
    common = {i for i in common
              if all(d[i]['pf'] < FREE_FRAC for d in jobs.values())}
    print(f"  modules with prior_frac < {FREE_FRAC} in all jobs: {len(common)}")
    mods = sorted(common)

    def stats(field, label):
        per_run, within, nsp = [], [], []
        for r in runs:
            ks = [k for k in jobs if k[0] == r]
            arr = np.array([[jobs[k][i][field] for i in mods] for k in ks])
            with np.errstate(all='ignore'):
                per_run.append(np.nanmean(arr, axis=0))
                nsp.append(len(ks))
                within.append((len(ks), np.nanvar(arr, axis=0, ddof=1)
                               if len(ks) > 1 else None))
        good = np.isfinite(np.vstack(per_run)).all(axis=0)
        pr = np.vstack(per_run)[:, good]
        if len(runs) < 2 or not good.any():
            print(f"  {label:<22} n/a (modules {int(good.sum())})")
            return
        sd_runs = pr.std(axis=0, ddof=1)
        num = sum((n - 1) * v[good] for n, v in within if v is not None)
        den = sum((n - 1) for n, v in within if v is not None)
        if den > 0:
            sd_split = np.sqrt(num / den)
            exp_mean = sd_split * math.sqrt(np.mean([1.0 / n for n in nsp]))
            s = (f"split-to-split SD mean {sd_split.mean():.4f} median "
                 f"{np.median(sd_split):.4f}; expected on a run mean mean "
                 f"{exp_mean.mean():.4f} median {np.median(exp_mean):.4f}")
        else:
            s = "split-to-split SD n/a (one split per run)"
        print(f"  {label:<22} run-to-run SD of run means: mean "
              f"{sd_runs.mean():.4f} median {np.median(sd_runs):.4f} "
              f"({int(good.sum())} modules, {len(runs)} runs); {s}")
    stats('g', 'percrystal g')
    stats('gsi', 'seed-indexed (table)')
    if any(np.isfinite(jobs[k][i]['win']) for k in jobs for i in mods[:50]):
        stats('win', 'window.tsv gain')
    return 0


def cmd_compare(args):
    hdr, cols, rows = _read_tsv(args.gains)
    rows = [r for r in rows if r.get("path", args.path) == args.path]
    if args.window:
        _, _, wr = _read_tsv(args.window)
        gm = {r["module"]: float(r["gain"]) for r in wr
              if r.get("path") == args.path and r.get("used") == "1"}
        src = args.window
    else:
        gm = {r["module"]: float(r["g_seedindexed"]) for r in rows
              if math.isfinite(float(r["g_seedindexed"]))}
        src = "g_seedindexed column"
    names, gi, gmv, pfv, sgv = [], [], [], [], []
    for r in rows:
        if r["module"] in gm and math.isfinite(gm[r["module"]]):
            names.append(r["module"])
            gi.append(float(r["g"]))
            gmv.append(gm[r["module"]])
            pfv.append(float(r["prior_frac"]))
            sgv.append(float(r["sigma_g"]))
    gi, gmv, pfv, sgv = map(np.array, (gi, gmv, pfv, sgv))
    d = gi - gmv
    print(f"[compare] {args.gains}  vs  {src}  (path {args.path})")
    if len(d) == 0:
        print("  no common modules")
        return 1
    print(f"  N = {len(d)}  RMS(g_i - g_M) = {np.sqrt(np.mean(d**2)):.5f}  "
          f"median = {np.median(d):+.5f}  mean = {d.mean():+.5f}  "
          f"corr = {np.corrcoef(gi, gmv)[0, 1]:.4f}")
    big = np.nonzero(np.abs(d) > 0.02)[0]
    print(f"  modules with |g_i - g_M| > 2 %: {len(big)}")
    for j in big:
        print(f"    {names[j]:<7} g_i {gi[j]:.4f} +- {sgv[j]:.4f}  "
              f"g_M {gmv[j]:.4f}  diff {d[j]*100:+.2f} %  prior_frac "
              f"{pfv[j]:.2f}")
    return 0


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    if "--self-test" in argv:
        return self_test()
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    f = sub.add_parser("fit", help="pooled per-crystal fit of cluster dumps "
                                   "(one or several runs, each pooled alone)")
    f.add_argument("--dumps", nargs="+", required=True, metavar="GLOB")
    f.add_argument("--splits", choices=['even', 'odd', 'all'], default='all',
                   help="0-based position in the run's sorted split list "
                        "(default all).")
    f.add_argument("--out", required=True, metavar="DIR")
    f.add_argument("--max-iter", type=int, default=MAX_ITER)
    f.add_argument("--fixed-windows", action="store_true",
                   help="select the events once, with the windows of the g = 1 (seed-indexed) pass in pre-refinement energy; the solve is then independent of the starting gains (P1-6).  Default: windows refit in passes 1-3 on E_cal (the algorithm of the 2026-10-23 scans)")
    f.add_argument("--exclude-rowcol", type=int, nargs=4, default=None,
                   metavar=("R_LO", "R_HI", "C_LO", "C_HI"),
                   help="Rectangle to drop from the pool, 1-based "
                        "inclusive rows/cols (default: the dump header's "
                        "value; the dumps of output/scan_20261023 carry 15 "
                        "20 15 20, now read 1-based, i.e. centred).")
    s = sub.add_parser("stability", help="run-to-run vs split-to-split gain "
                                         "spread from per-job gain tables")
    s.add_argument("--gains", nargs="+", required=True, metavar="GLOB")
    s.add_argument("--path", default="legacy")
    c = sub.add_parser("compare", help="g_i vs seed-indexed g_M of one job")
    c.add_argument("--gains", required=True, metavar="TSV")
    c.add_argument("--window", default=None, metavar="TSV",
                   help="_window.tsv of a seed-indexed (default elastic-window) "
                        "job on the same events, not of a percrystal job (its "
                        "gain column is the residual of the last pass); "
                        "default: the g_seedindexed column of --gains.")
    c.add_argument("--path", default="legacy")
    args = ap.parse_args(argv)
    return {'fit': cmd_fit, 'stability': cmd_stability,
            'compare': cmd_compare}[args.cmd](args)


if __name__ == "__main__":
    sys.exit(main())
