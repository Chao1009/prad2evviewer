#!/usr/bin/env python3
"""template_amp_binned.py -- does the PbWO4 pulse shape depend on pulse amplitude?

Question: the best-50 spline templates are built from the top-50 mf_energy pulses of a module,
so the pool amplitude mix differs between run types (2.24 GeV H2: 12 % of the pool below
600 ADC, Ta (X17): 45 %, 0.73 GeV H2: most).  If the shape depends on amplitude, that mix
alone moves the templates.  Does it, by how much, and does it explain (a) the per-module
Ta - H2 scatter and (b) the 0.73 GeV bank being narrower?

1. Per-bin templates.  Candidates of the 2.24 GeV H2 select pools (the pools the bank is built
   from; tier >= 1, band b0 [300, 2500] ADC; optionally one run) are split by amplitude
   (default 300-600 / 600-1200 / 1200-2500).  A module is used if every bin has
   n_pool + n_ho pulses.  Per module and bin, n_pool pulses are drawn at random (seeded)
   and n_ho more from the same bin are the hold-out of the lambda selection; the template is
   spline_template_extract.process_module on exactly these pulses, with the options of the
   shipped bank (npz meta.argv, loaded as in template_doublet_audit.load_ctx).
2. Noise.  The same build on two disjoint draws of one bin (the null, equal N) gives the
   split-half noise of a bin difference; sigma = 1.4826 MAD over modules.
3. Metrics of the peak-normalised shapes: rise 10-90, w99, FWHM, t_pk (parabolic) and max|df|
   (template_doublet_audit.shape_metrics / maxdiff).
4. f(A) check.  r_b / r_mod of the bin templates, as spline_template_extract._fa_ratios does
   (noise-free integer waveforms of the bin template at 16 phases, fitted with the module
   template by mf_pulse), with the analytic discrete walk-out integral H * I_cont in place of
   the C++ integral; validated on the shipped lo / mid / top band templates against the
   stored r_lo / r_mid / r_top and the ampcorr knots.
5. Ta - H2 and 0.73 - 2.24 H2: per-module shape differences of the shipped banks against the
   pool amplitude mix (mean ln A of the shipped pools), predicted with the amplitude slope.
6. Reproduction: process_module on a module's whole candidate set and on the rebuilt pool
   with the shipped hold-out must give the shipped template (max|df| = 0).

    MPLBACKEND=Agg python3 template_amp_binned.py                    # all 5 runs, n 40 + 12
    MPLBACKEND=Agg python3 template_amp_binned.py --run 25320 --n-pool 25 --n-ho 8 \\
        --out output/templates_spline_best50_20261026/amp_binned/run25320

Output: output/templates_spline_best50_20261026/amp_binned/ (report.md, modules.tsv,
ta_h2_modules.tsv, amp_binned.pdf / .png).  numpy / scipy / matplotlib.
"""
from __future__ import annotations

import argparse
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import spearmanr, theilslopes

import template_doublet_audit as A
import spline_template_extract as SE

HERE = Path(__file__).resolve().parent
B50 = A.B50
SLUG = "2p24GeV_H2"
SLUG_TA = "2p24GeV_Ta_X17"
SLUG_073 = "0p73GeV_H2"
TG = SE.T_GRID
MET = ["rise", "w99", "fwhm", "tpk"]
UNIT = {"rise": "ns", "w99": "ns", "fwhm": "ns", "tpk": "ns", "maxd": "%"}
_CTX = None
_Z = None


# ----------------------------------------------------------------- helpers --
def rsig(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if x.size > 2 else np.nan


def pct(x, q=(16, 50, 84)):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return np.percentile(x, q) if x.size else np.full(len(q), np.nan)


def med_se(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return (float(np.median(x)), float(1.2533 * rsig(x) / np.sqrt(x.size))) if x.size > 2 else (np.nan, np.nan)


def fmt(m, lo, hi, nd=2, sign=True):
    f = f"{{:+.{nd}f}}" if sign else f"{{:.{nd}f}}"
    return f"{f.format(m)} [{f.format(lo)}, {f.format(hi)}]"


# ------------------------------------------------------------------ builds --
def build_pool(mi, pool, ho):
    """process_module on `pool` with `ho` as the hold-out of the lambda selection."""
    ctx = _CTX
    hold = ctx["hold"]
    pool = np.asarray(pool, np.int64)
    ho = np.asarray(ho, np.int64)
    idx = np.unique(np.concatenate([pool, ho]))
    saved = hold[idx].copy()
    hold[pool] = False
    hold[ho] = True
    try:
        return SE.process_module(ctx, mi, idx)
    finally:
        hold[idx] = saved


def _summ(r):
    kg = r.get("kept_gidx")
    return dict(fb=int(r["fallback"]), shape=r.get("shape"), lam=float(r.get("lam_final", np.nan)),
                n_pool=int(r["n_pool"]), n_kept=None if kg is None else int(kg.size),
                spl=(r.get("t"), r.get("c")))


def _worker_init():
    SE.BAND_MIN = 10 ** 9           # skip the lo / hi / mid / top band builds (shape unaffected)


def _task(t):
    kind, key, mi, pool, ho = t
    if kind == "repro_full":
        r = SE.process_module(_CTX, mi, A.module_idx(_CTX, mi))
    else:
        r = build_pool(mi, pool, ho)
    return key, _summ(r)


def _fa_task(t):
    return _fa_one(*t)


# --------------------------------------------------------------- f(A) proxy --
_LUT_N = None


def _lut_row(shape):
    import mf_pulse
    tg = mf_pulse._LUT_DT * np.arange(int(round(SE.T_END / mf_pulse._LUT_DT)) + 1)
    return np.interp(tg, TG, np.nan_to_num(np.asarray(shape, np.float64)), left=0.0, right=0.0)


def _fa_one(mod_shape, bins, H_mod):
    """r_b / r_mod of the bin templates `bins` = [(shape, H)] against the module template.
    As _fa_ratios, with Q = H * I_cont[row] (analytic) instead of the C++ integral."""
    import mf_pulse
    rows = [_lut_row(mod_shape)] + [_lut_row(s) for s, _ in bins]
    nm = len(rows)
    bank = mf_pulse.SampledTemplateBank(rows, [f"s{i}" for i in range(nm)], [True] * nm,
                                        smooth_order=2, int_tail_ratio=0.1, tail_break_n=2,
                                        T_end=SE.T_END)
    mf = mf_pulse.MFReconstructor(bank, ped_nsamples=30, **SE.FA_MF)
    phi = SE.DT_NS * np.arange(SE.FA_NPH) / SE.FA_NPH
    out = []
    for k, H in [(0, H_mod)] + [(i + 1, h) for i, (_, h) in enumerate(bins)]:
        t0 = 140.0 + phi
        S = bank.shape(np.full((SE.FA_NPH, 1), k), SE.DT_NS * np.arange(SE.NS_EXPECT)[None, :] - t0[:, None],
                       dtype=np.float64)
        Y = np.rint(200.0 + H * S).astype(np.uint16)
        res = mf.run(Y, np.zeros(SE.FA_NPH, np.int64))
        E = np.full(SE.FA_NPH, np.nan)
        one = np.nonzero(res.npk == 1)[0]
        if one.size:
            pos = np.searchsorted(res.p_ch, one)
            E[one] = mf_pulse.pulse_energy_units(bank, np.zeros(one.size, np.int64), res.p_t0[pos],
                                                 res.p_amp[pos], kappa=1.0)
        r = H * bank.I_cont[k] / E
        out.append(float(np.nanmean(r)) if np.isfinite(r).any() else np.nan)
    return np.array(out[1:]) / out[0]


# --------------------------------------------------------------------- main --
def parse(argv):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run", default="all", help="candidate run number of the select pools, or 'all' [all]")
    p.add_argument("--edges", default="300,600,1200,2500", help="amplitude bin edges, ADC [300,600,1200,2500]")
    p.add_argument("--n-pool", type=int, default=40, help="pulses per bin template [40]")
    p.add_argument("--n-ho", type=int, default=12, help="same-bin hold-out pulses of the lambda selection [12]")
    p.add_argument("--keep-shipped", action="store_true",
                   help="do not exclude the shipped bank's pool and hold-out pulses from the bins (default: excluded, so the "
                        "bins are independent of the shipped 2.24 H2 template they are compared with)")
    p.add_argument("--workers", type=int, default=8, help="processes (<= 8) [8]")
    p.add_argument("--seed", type=int, default=20261028)
    p.add_argument("--quick", type=int, default=0, help="use every k-th selected module (test)")
    p.add_argument("--n-fa-valid", type=int, default=120, help="modules of the f(A) proxy validation [120]")
    p.add_argument("--n-repro", type=int, default=4, help="modules of the reproduction check [4]")
    p.add_argument("--out", default=str(B50 / "amp_binned"), help="output directory")
    return p.parse_args(argv)


def pool_amps(slug):
    """Per module: shipped pool (top-50 mf_energy) amplitude summary, from the audit's own rebuild."""
    ctx, z = A.load_ctx(slug)
    D = ctx["D"]
    fbk = np.asarray(z["fallback"])
    out = {}
    for r in range(len(fbk)):
        if fbk[r] in (0, 1):
            p, _ = A.module_pool(ctx, int(z["module_index"][r]), int(fbk[r]))
            a = D["amp"][p].astype(np.float64)
            out[str(z["names"][r])] = dict(n=int(p.size), lnA=float(np.mean(np.log(a))), med=float(np.median(a)),
                                          flo=float(np.mean(a < 600)), a=a)
    zz = {k: np.asarray(z[k]) for k in ("names", "fallback", "shape", "shape_avg", "module_index")}
    return out, zz


def main(argv=None):
    global _CTX
    a = parse(argv)
    t_start = time.time()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    edges = [float(x) for x in a.edges.split(",")]
    nb = len(edges) - 1
    bl = [f"{edges[i]:.0f}-{edges[i + 1]:.0f}" for i in range(nb)]
    workers = min(int(a.workers), 8)
    npl, nho = a.n_pool, a.n_ho
    blk = npl + nho

    # ---------------------------------------------------------------- setup
    ctx, z = A.load_ctx(SLUG)
    _CTX = ctx
    D = ctx["D"]
    names = [str(x) for x in z["names"]]
    fbk = np.asarray(z["fallback"])
    mid_of_row = np.asarray(z["module_index"])
    shape_ship = np.asarray(z["shape"], np.float64)
    am = json.loads(str(z["ampcorr_meta"]))
    ok = (D["tier"] >= 1) & D["band_b0"] & np.isfinite(D["mf_energy"])
    if a.run != "all":
        ok &= D["run"] == int(a.run)
    amp = D["amp"].astype(np.float64)
    binof = np.digitize(amp, edges[1:-1]) if nb > 1 else np.zeros(amp.size, int)
    inrange = (amp >= edges[0]) & (amp <= edges[-1])
    ok &= inrange
    print(f"[setup] {D['midx'].size} candidates, {int(ok.sum())} eligible "
          f"(run {a.run}, {edges[0]:.0f}-{edges[-1]:.0f} ADC); bins {bl}; N {npl} + {nho}; "
          f"shipped pool + hold-out {'kept' if a.keep_shipped else 'excluded'}", flush=True)

    sel = []                              # module rows with every bin populated
    rng_all = {}
    for r in range(len(names)):
        if fbk[r] != 0:
            continue
        mi = int(mid_of_row[r])
        idx = A.module_idx(ctx, mi)
        idx = idx[ok[idx]]
        if not a.keep_shipped:
            pool_s, idx_all = A.module_pool(ctx, mi, 0)
            idx = np.setdiff1d(idx, np.concatenate([pool_s, idx_all[ctx["hold"][idx_all]]]))
        pools, epools = [], []
        good = True
        rng = np.random.default_rng([a.seed, mi])
        for b in range(nb):
            ib = idx[binof[idx] == b]
            if ib.size < blk:
                good = False
                break
            perm = rng.permutation(ib)
            nblk = min(int(ib.size // blk), 2)             # a second disjoint block = the null draw
            pools.append([perm[k * blk:(k + 1) * blk] for k in range(nblk)] + [ib.size])
            epools.append(ib[np.argsort(-D["mf_energy"][ib], kind="stable")][:blk])     # the bank's rule, within the bin
        if good:
            sel.append((r, mi, pools, epools))
    if a.quick:
        sel = sel[::a.quick]
    print(f"[setup] {len(sel)} modules with >= {blk} pulses in every bin; "
          f"null draws per bin: " + ", ".join(f"{bl[b]} {sum(len(p[b]) == 3 for _, _, p, _e in sel)}" for b in range(nb)),
          flush=True)

    tasks = []
    for r, mi, pools, epools in sel:
        for b in range(nb):
            for k in range(len(pools[b]) - 1):
                blkidx = pools[b][k]
                tasks.append((("bin" if k == 0 else "null", r, b), ("bin", (k, r, b), mi, blkidx[:npl], blkidx[npl:])))
            tasks.append((("ener", r, b), ("bin", (2, r, b), mi, epools[b][:npl], epools[b][npl:])))
    # reproduction (shipped pool + shipped hold-out; the build's own whole-module call)
    rng = np.random.default_rng(a.seed)
    rows0 = [r for r in range(len(names)) if fbk[r] == 0]
    spot = [int(x) for x in rng.choice(rows0, min(a.n_repro, len(rows0)), replace=False)]
    for r in spot:
        mi = int(mid_of_row[r])
        pool, idx = A.module_pool(ctx, mi, 0)
        h = idx[ctx["hold"][idx]]
        tasks.append((("repro_pool", r), ("repro_pool", ("rp", r), mi, pool, h)))
        tasks.append((("repro_full", r), ("repro_full", ("rf", r), mi, None, None)))
    print(f"[build] {len(tasks)} process_module calls on {workers} workers", flush=True)
    t0 = time.time()
    pw = mp.get_context("fork").Pool(workers, initializer=_worker_init)
    try:
        res = pw.map(_task, [t[1] for t in tasks], chunksize=1)
    finally:
        pw.close()
        pw.join()
    print(f"[build] done in {time.time() - t0:.0f} s", flush=True)
    R = {k: v for k, v in res}

    # ---------------------------------------------------------- reproduction
    repro = []
    for r in spot:
        sh = shape_ship[r]
        d_pool = float(np.nanmax(np.abs(np.asarray(R[("rp", r)]["shape"], np.float64) - sh)))
        d_full = float(np.nanmax(np.abs(np.asarray(R[("rf", r)]["shape"], np.float64) - sh)))
        repro.append((names[r], d_pool, d_full))
    print("[repro] " + "; ".join(f"{n}: pool {dp:.1e} full {df:.1e}" for n, dp, df in repro), flush=True)

    # --------------------------------------------------------------- metrics
    def mets(shape):
        if shape is None:
            return np.full(5, np.nan)
        s = np.asarray(shape, np.float64)
        return np.append(A.shape_metrics(s), 0.0)

    rec = []                              # per module
    for r, mi, pools, epools in sel:
        d = dict(row=r, name=names[r], midx=mi, ok=True, m=[], A=[], lam=[], fb=[], npool=[], lnA=[], null=[], null_shape=[],
                 shape=[], nin=[], mE=[], AE=[], lnAE=[], lamE=[], shapeE=[])
        for b in range(nb):
            o = R[(0, r, b)]
            if o["shape"] is None or o["fb"] != 0:
                d["ok"] = False
                break
            ai = amp[pools[b][0][:npl]]
            d["m"].append(mets(o["shape"])[:4])
            d["A"].append(float(np.median(ai)))
            d["lnA"].append(float(np.mean(np.log(ai))))
            d["lam"].append(o["lam"])
            d["fb"].append(o["fb"])
            d["npool"].append(o["n_pool"])
            d["shape"].append(np.asarray(o["shape"], np.float32))
            d["nin"].append(pools[b][-1])
            oe = R[(2, r, b)]
            if oe["shape"] is None or oe["fb"] != 0:
                d["ok"] = False
                break
            ae = amp[epools[b][:npl]]
            d["mE"].append(mets(oe["shape"])[:4])
            d["AE"].append(float(np.median(ae)))
            d["lnAE"].append(float(np.mean(np.log(ae))))
            d["lamE"].append(oe["lam"])
            d["shapeE"].append(np.asarray(oe["shape"], np.float32))
            nk = (1, r, b)
            if nk in R and R[nk]["shape"] is not None and R[nk]["fb"] == 0:
                d["null"].append(np.append(mets(R[nk]["shape"])[:4] - d["m"][-1], A.maxdiff(R[nk]["shape"], o["shape"])))
            else:
                d["null"].append(None)
        if d["ok"]:
            rec.append(d)
    nrec = len(rec)
    print(f"[metrics] {nrec} of {len(sel)} modules with all bins fallback 0", flush=True)
    M = np.array([d["m"] for d in rec])                   # [mod, bin, 4]
    Ab = np.array([d["A"] for d in rec])
    lnAb = np.array([d["lnA"] for d in rec])
    LAM = np.array([d["lam"] for d in rec])
    ME = np.array([d["mE"] for d in rec])
    AbE = np.array([d["AE"] for d in rec])
    lnAbE = np.array([d["lnAE"] for d in rec])
    LAME = np.array([d["lamE"] for d in rec])
    lnA_med = np.log(np.median(Ab, axis=0))               # bin axis (median amplitude over modules)
    # null sigmas per bin (difference of two disjoint draws)
    sig_null = np.full((nb, 5), np.nan)
    n_null = np.zeros(nb, int)
    null_all = []
    for b in range(nb):
        x = np.array([d["null"][b] for d in rec if d["null"][b] is not None])
        n_null[b] = len(x)
        null_all.append(x)
        if len(x) > 5:
            sig_null[b] = [rsig(x[:, k]) for k in range(5)]
    sig_fill = np.nanmean(sig_null, axis=0)               # bins without null draw: mean of the others
    sig_b = np.array([sig_null[b] if np.isfinite(sig_null[b]).all() else sig_fill for b in range(nb)])

    pairs = [(0, 1), (1, 2), (0, 2)] if nb == 3 else [(0, nb - 1)]
    rows_tab = []
    pair_res = {}
    for (i, j) in pairs:
        dM = M[:, j, :] - M[:, i, :]                      # high - low amplitude
        dmax = np.array([A.maxdiff(rec[m]["shape"][j], rec[m]["shape"][i]) for m in range(nrec)])
        sp = np.sqrt((sig_b[i, :4] ** 2 + sig_b[j, :4] ** 2) / 2.0)      # sigma of a difference of two N-pulse templates
        for k, mn in enumerate(MET):
            x = dM[:, k]
            lo_, m_, hi_ = pct(x)
            _, se = med_se(x)
            z2 = x / sp[k]
            n2p, n2m = int(np.sum(z2 > 2)), int(np.sum(z2 < -2))
            rows_tab.append((f"{bl[j]} - {bl[i]}", mn, m_, lo_, hi_, se, sp[k], n2p, n2m, int(np.isfinite(x).sum()), (hi_ - lo_) / 2))
        mx_lo, mx_m, mx_hi = pct(dmax)
        nulmax = [np.median(null_all[b][:, 4]) if len(null_all[b]) else np.nan for b in (i, j)]
        rows_tab.append((f"{bl[j]} - {bl[i]}", "maxd", mx_m, mx_lo, mx_hi, med_se(dmax)[1], np.nanmean(nulmax), 0, 0,
                         nrec, np.nan))
        pair_res[(i, j)] = dict(dM=dM, dmax=dmax, sp=sp)

    # per-module amplitude slope: metric vs ln A (all bins), least squares per module
    slopes = np.full((nrec, 4), np.nan)
    for m in range(nrec):
        x = lnAb[m]
        for k in range(4):
            slopes[m, k] = np.polyfit(x, M[m, :, k], 1)[0] if np.isfinite(M[m, :, k]).all() else np.nan
    slope_med = np.array([med_se(slopes[:, k]) for k in range(4)])          # [metric, (median, se)]
    slope_pct = np.array([pct(slopes[:, k]) for k in range(4)])
    # null of the slope: 3-point slope of independent draws is not available; use lo-hi slope sigma
    i0, j0 = pairs[-1]
    dlnA_pair = np.median(lnAb[:, j0] - lnAb[:, i0])
    slope_lohi = np.array([med_se((M[:, j0, k] - M[:, i0, k]) / (lnAb[:, j0] - lnAb[:, i0])) for k in range(4)])

    # bank-like selection (top-N mf_energy within the bin): slopes, lo -> hi differences, and E - R per bin
    slopes_E = np.full((nrec, 4), np.nan)
    for m in range(nrec):
        for k in range(4):
            slopes_E[m, k] = np.polyfit(lnAbE[m], ME[m, :, k], 1)[0] if np.isfinite(ME[m, :, k]).all() else np.nan
    slope_med_E = np.array([med_se(slopes_E[:, k]) for k in range(4)])
    slope_pct_E = np.array([pct(slopes_E[:, k]) for k in range(4)])
    dME = ME[:, j0, :] - ME[:, i0, :]
    dmaxE = np.array([A.maxdiff(rec[m]["shapeE"][j0], rec[m]["shapeE"][i0]) for m in range(nrec)])
    spE = pair_res[(i0, j0)]["sp"]
    rowsE = []
    for k, mn in enumerate(MET):
        x = dME[:, k]
        lo_, m_, hi_ = pct(x)
        rowsE.append((f"{bl[j0]} - {bl[i0]} (energy-top)", mn, m_, lo_, hi_, med_se(x)[1], spE[k], int(np.sum(x / spE[k] > 2)),
                      int(np.sum(x / spE[k] < -2)), int(np.isfinite(x).sum()), (hi_ - lo_) / 2))
    mx_lo, mx_m, mx_hi = pct(dmaxE)
    rowsE.append((f"{bl[j0]} - {bl[i0]} (energy-top)", "maxd", mx_m, mx_lo, mx_hi, med_se(dmaxE)[1],
                  float(np.nanmean([np.median(null_all[b][:, 4]) for b in (i0, j0)])), 0, 0, nrec, np.nan))
    dER = ME - M                                          # energy-top minus random, same bin and N
    er_tab = []
    for b in range(nb):
        row = [bl[b], float(np.median(AbE[:, b])), float(np.median(Ab[:, b]))]
        for k in range(3):
            lo_, m_, hi_ = pct(dER[:, b, k])
            row += [m_, lo_, hi_, med_se(dER[:, b, k])[1]]
        er_tab.append(row)
    lam_E = [(bl[b], float(np.mean(LAME[:, b] <= 10)), float(np.mean(LAME[:, b] == 100)), float(np.mean(LAME[:, b] == 1000)),
              float(np.mean(LAME[:, b] >= 1e4))) for b in range(nb)]

    # --------------------------------------------------------- lambda / bins
    lam_tab = []
    for b in range(nb):
        v = LAM[:, b]
        lam_tab.append((bl[b], float(np.median(Ab[:, b])), float(np.mean(v <= 10)), float(np.mean(v == 100)),
                        float(np.mean(v == 1000)), float(np.mean(v >= 1e4))))

    # ------------------------------------------------- f(A): bin proxy + validation
    print("[fa] r_b / r_mod proxy", flush=True)
    t0 = time.time()
    H_b0 = float(am["H_b0"]) if "H_b0" in am else float(np.median(np.concatenate([D["amp"][A.module_idx(ctx, int(mid_of_row[d["row"]]))] for d in rec[:50]])))
    fa_tasks = []
    for d in rec:
        fa_tasks.append((shape_ship[d["row"]], [(d["shape"][b], float(np.median(Ab[:, b]))) for b in range(nb)], H_b0))
    for d in rec:
        fa_tasks.append((shape_ship[d["row"]], [(d["shapeE"][b], float(np.median(AbE[:, b]))) for b in range(nb)], H_b0))
    # validation on the shipped band templates (lo / mid / top)
    zf = np.load(B50 / SLUG / f"templates_{SLUG}.npz", allow_pickle=False)
    knots = am["knots"]
    kA = [knots[k]["A"] for k in ("lo", "mid", "top")]
    shp_b = [np.asarray(zf["shape_lo"], np.float64), np.asarray(zf["shape_mid"], np.float64),
             np.asarray(zf["shape_top"], np.float64)]
    okb = [r for r in range(len(names)) if fbk[r] == 0 and all(np.isfinite(s[r]).all() for s in shp_b)
           and np.isfinite(zf["r_mod"][r]) and all(np.isfinite(zf[k][r]) for k in ("r_lo", "r_mid", "r_top"))]
    vr = np.random.default_rng(a.seed + 1)
    vrows = [int(x) for x in vr.choice(okb, min(a.n_fa_valid, len(okb)), replace=False)]
    for r in vrows:
        fa_tasks.append((shape_ship[r], [(s[r], h) for s, h in zip(shp_b, kA)], H_b0))
    pw = mp.get_context("fork").Pool(workers)
    try:
        fa_out = pw.map(_fa_task, fa_tasks, chunksize=2)
    finally:
        pw.close()
        pw.join()
    g_bins = np.array(fa_out[:nrec])                       # [mod, bin] = r_bin / r_mod (proxy)
    g_binsE = np.array(fa_out[nrec:2 * nrec])
    g_val = np.array(fa_out[2 * nrec:])                    # [mod, 3]
    print(f"[fa] {len(fa_tasks)} banks in {time.time() - t0:.0f} s", flush=True)
    rc = np.array([[zf["r_lo"][r], zf["r_mid"][r], zf["r_top"][r]] for r in vrows]) / np.array([zf["r_mod"][r] for r in vrows])[:, None]
    val_ratio = {
        "lo": med_se(g_val[:, 0]), "mid_lo": med_se(g_val[:, 1] / g_val[:, 0]), "top_mid": med_se(g_val[:, 2] / g_val[:, 1]),
        "cpp_lo": med_se(rc[:, 0]), "cpp_mid_lo": med_se(rc[:, 1] / rc[:, 0]), "cpp_top_mid": med_se(rc[:, 2] / rc[:, 1]),
        "tab_lo": 1.0, "tab_mid_lo": knots["mid"]["f"] / knots["lo"]["f"], "tab_top_mid": knots["top"]["f"] / knots["mid"]["f"],
        "diff_lo": med_se(g_val[:, 0] - rc[:, 0]), "diff_mid": med_se(g_val[:, 1] - rc[:, 1]), "diff_top": med_se(g_val[:, 2] - rc[:, 2]),
        "n": len(vrows)}
    la, fa_ = np.asarray(am["logA"]), np.asarray(am["f"])
    fA_at = lambda A_: float(np.interp(np.log(A_), la, fa_))
    Abm = np.median(Ab, axis=0)
    fbin = [med_se(g_bins[:, b]) for b in range(nb)]
    fa_tab = []
    for b in range(nb):
        fa_tab.append((bl[b], float(Abm[b]), fbin[b][0], fbin[b][1], fA_at(Abm[b])))
    AbmE = np.median(AbE, axis=0)
    fbinE = [med_se(g_binsE[:, b]) for b in range(nb)]
    fa_tabE = [(bl[b], float(AbmE[b]), fbinE[b][0], fbinE[b][1], fA_at(AbmE[b])) for b in range(nb)]
    fstep = []
    for (i, j) in pairs:
        s_ = med_se(g_bins[:, j] / g_bins[:, i])
        s2 = med_se(g_binsE[:, j] / g_binsE[:, i])
        fstep.append((f"{bl[j]} / {bl[i]}", s_[0], s_[1], fA_at(Abm[j]) / fA_at(Abm[i]), s2[0], s2[1],
                      fA_at(AbmE[j]) / fA_at(AbmE[i])))

    # where does the integral change come from: unit-peak area of the core (0-50 ns) and of the tail (50 ns - end)
    core = (TG > 0) & (TG <= 50)
    tail = TG > 50
    def area_ratio(key, i, j):
        out_ = []
        for m_ in range(nrec):
            a_ = np.nan_to_num(np.asarray(rec[m_][key][i], np.float64))
            b_ = np.nan_to_num(np.asarray(rec[m_][key][j], np.float64))
            out_.append([np.trapezoid(b_[msk], TG[msk]) / np.trapezoid(a_[msk], TG[msk]) for msk in (core, tail)]
                        + [np.trapezoid(b_[TG > 0], TG[TG > 0]) / np.trapezoid(a_[TG > 0], TG[TG > 0])])
        out_ = np.array(out_)
        return [med_se(out_[:, c_]) for c_ in range(3)]
    area_tab = dict(R=area_ratio("shape", i0, j0), E=area_ratio("shapeE", i0, j0))

    # ----------------------------------------- Ta - H2, 0.73 - 2.24 H2 (shipped banks)
    print("[bank] Ta - H2 and 0.73 - 2.24 H2", flush=True)
    _CTX = None
    del ctx, D
    pa_h2, zh2 = pool_amps(SLUG)
    pa_ta, zta = pool_amps(SLUG_TA)
    pa_07, z07 = pool_amps(SLUG_073)

    def bank_mets(zz):
        nm = [str(x) for x in zz["names"]]
        return {n: A.shape_metrics(zz["shape"][i]) for i, n in enumerate(nm) if zz["fallback"][i] <= 1}

    mh2, mta, m07 = bank_mets(zh2), bank_mets(zta), bank_mets(z07)
    avg = {k: A.shape_metrics(zz["shape_avg"]) for k, zz in (("2.24 H2", zh2), ("Ta", zta), ("0.73 H2", z07))}

    slope_of = {d["name"]: (slopes[m], slopes_E[m]) for m, d in enumerate(rec)}

    def compare(mx, pax, my, pay, tag):
        """x - y per common module: observed shape differences vs pool mix differences."""
        nm = [n for n in mx if n in my and n in pax and n in pay]
        dobs = np.array([mx[n] - my[n] for n in nm])                            # [mod, 4]
        dln = np.array([pax[n]["lnA"] - pay[n]["lnA"] for n in nm])
        dlo = np.array([pax[n]["flo"] - pay[n]["flo"] for n in nm])
        flo_x = np.array([pax[n]["flo"] for n in nm])
        medA_x = np.array([pax[n]["med"] for n in nm])
        out = dict(tag=tag, names=nm, dobs=dobs, dln=dln, dlo=dlo, flo=flo_x, medA=medA_x,
                   corr={}, pred={}, fit={})
        for k, mn in enumerate(MET):
            x = dobs[:, k]
            fin = np.isfinite(x)
            out["corr"][mn] = dict(
                lnA=spearmanr(dln[fin], x[fin]), lo=spearmanr(dlo[fin], x[fin]), flo=spearmanr(flo_x[fin], x[fin]))
            pr = slope_med[k, 0] * dln
            prE = slope_med_E[k, 0] * dln
            out["pred"][mn] = dict(predE=pct(prE), residE=pct(x - prE), pred=pct(pr), obs=pct(x), pred_med_se=med_se(pr), obs_med_se=med_se(x),
                                   resid=pct(x - pr), resid_sig=rsig(x - pr), obs_sig=rsig(x))
            ts = theilslopes(x[fin], dln[fin])
            out["fit"][mn] = dict(slope=float(ts[0]), lo=float(ts[2]), hi=float(ts[3]))
            # per-module prediction: the module's own amplitude slope x its d<ln A>
            cm = [i for i, n in enumerate(nm) if n in slope_of and np.isfinite(x[i])]
            pm = {}
            for tg_, kk in (("R", 0), ("E", 1)):
                xs = x[cm]
                ps = np.array([slope_of[nm[i]][kk][k] * dln[i] for i in cm])
                g = np.isfinite(ps)
                rho = spearmanr(ps[g], xs[g])
                tsl = theilslopes(xs[g], ps[g])
                pm[tg_] = dict(n=int(g.sum()), rho=float(rho[0]), p=float(rho[1]), slope=float(tsl[0]), lo=float(tsl[2]),
                               hi=float(tsl[3]), sd_pred=rsig(ps[g]), sd_obs=rsig(xs[g]))
            out.setdefault("permod", {})[mn] = pm
        return out

    cmp_ta = compare(mta, pa_ta, mh2, pa_h2, "Ta - 2.24 H2")
    cmp_07 = compare(m07, pa_07, mh2, pa_h2, "0.73 - 2.24 H2")
    # run-type averages: pool mix = mean over modules of the pool mean ln A
    mean_lnA = {k: float(np.mean([v["lnA"] for v in pa.values()])) for k, pa in
                (("2.24 H2", pa_h2), ("Ta", pa_ta), ("0.73 H2", pa_07))}
    mean_flo = {k: float(np.mean(np.concatenate([v["a"] < 600 for v in pa.values()]))) for k, pa in
                (("2.24 H2", pa_h2), ("Ta", pa_ta), ("0.73 H2", pa_07))}
    med_amp = {k: float(np.median(np.concatenate([v["a"] for v in pa.values()]))) for k, pa in
               (("2.24 H2", pa_h2), ("Ta", pa_ta), ("0.73 H2", pa_07))}

    # --------------------------------------------- bank band templates (lo / mid / top)
    bmet = {}
    for nm_, s_ in zip(("lo", "mid", "top"), shp_b):
        bmet[nm_] = np.array([A.shape_metrics(s_[r]) if (fbk[r] == 0 and np.isfinite(s_[r]).all()) else np.full(4, np.nan)
                              for r in range(len(names))])
    bank_band = {}
    for (p, q) in (("lo", "mid"), ("mid", "top"), ("lo", "top")):
        d = bmet[q] - bmet[p]
        bank_band[(p, q)] = [pct(d[:, k]) for k in range(4)] + [int(np.isfinite(d[:, 0]).sum())]
    lnk = np.log(kA)
    bank_slope = {mn: med_se((bmet["top"][:, k] - bmet["lo"][:, k]) / (lnk[2] - lnk[0])) for k, mn in enumerate(MET)}

    # ---------------------------------------------------------------- outputs
    write_all(a, out, dict(
        bl=bl, edges=edges, nb=nb, npl=npl, nho=nho, run=a.run, keep=a.keep_shipped, nsel=len(sel), nrec=nrec, rec=rec, M=M, Ab=Ab, lnAb=lnAb,
        rows_tab=rows_tab, sig_b=sig_b, sig_null=sig_null, n_null=n_null, slope_med=slope_med, slope_pct=slope_pct,
        slope_lohi=slope_lohi, dlnA_pair=dlnA_pair, lam_tab=lam_tab, repro=repro, fa_tab=fa_tab, fstep=fstep,
        val=val_ratio, cmp_ta=cmp_ta, cmp_07=cmp_07, mean_lnA=mean_lnA, mean_flo=mean_flo, med_amp=med_amp,
        avg=avg, bank_band=bank_band, bank_slope=bank_slope, kA=kA, pairs=pairs, pair_res=pair_res, g_bins=g_bins,
        null_all=null_all, slope_med_E=slope_med_E, slope_pct_E=slope_pct_E, rowsE=rowsE, er_tab=er_tab, lam_E=lam_E,
        fa_tabE=fa_tabE, area_tab=area_tab, AbE=AbE, lnAbE=lnAbE, ME=ME, slopes_E=slopes_E, g_binsE=g_binsE, pa_ta=pa_ta, pa_h2=pa_h2, pa_07=pa_07, slopes=slopes, lnA_med=lnA_med, knots=knots,
        wall=time.time() - t_start, workers=workers, argv=list(sys.argv)))
    print(f"[done] {time.time() - t_start:.0f} s -> {out}", flush=True)


# ------------------------------------------------------------------ writers --
def write_all(a, out, S):
    bl, nb = S["bl"], S["nb"]
    rec, M = S["rec"], S["M"]
    L = []
    P = L.append
    P(f"# Amplitude dependence of the PbWO4 pulse shape (2.24 GeV H2 select pools)\n")
    P(f"Generated by `template_amp_binned.py` ({time.strftime('%Y-%m-%d %H:%M')}), {S['wall']:.0f} s wall, {S['workers']} workers. "
      f"Command: `{' '.join(S['argv'])}`\n")
    P("## Setup\n")
    P(f"- Pools: 2.24 GeV H2 select candidates (tier >= 1, band b0 {S['edges'][0]:.0f}-{S['edges'][-1]:.0f} ADC), run `{S['run']}`. "
      f"Bins {', '.join(bl)} ADC; the shipped bank's pool and hold-out pulses are {'kept' if S['keep'] else 'excluded'}. Per module and bin: {S['npl']} random pulses (seeded) build the template with "
      f"`process_module` (options of the shipped bank), {S['nho']} more of the same bin are the lambda hold-out.")
    P(f"- Modules: {S['nrec']} (of {S['nsel']} with >= {S['npl'] + S['nho']} pulses in every bin; all bins fallback 0). "
      f"Median amplitude of the bin pulses: " + ", ".join(f"{bl[b]}: {np.median(S['Ab'][:, b]):.0f}" for b in range(nb)) + " ADC.")
    P(f"- Noise: the same build on two disjoint draws of one bin. Modules with a null draw per bin: "
      + ", ".join(f"{bl[b]} {S['n_null'][b]}" for b in range(nb)) + ". sigma of one bin-pair difference (null sigma of a difference of two "
      f"{S['npl']}-pulse templates, robust MAD) rise / w99 / FWHM / t_pk (ns): "
      + "; ".join(f"{bl[b]}: " + " / ".join(f"{S['sig_b'][b, k]:.2f}" for k in range(4)) for b in range(nb)) + ".")
    P("- Reproduction (max|df| of the unit-peak shape vs the shipped template; rebuilt pool with shipped hold-out / the build's whole-module call): "
      + "; ".join(f"{n} {dp:.1e} / {df:.1e}" for n, dp, df in S["repro"]) + ".\n")
    P("## 1. Bin-to-bin shape differences (high - low amplitude), median [16, 84 %]\n")
    P("| pair | metric | median [16, 84 %] | SE(median) | noise sigma | excess scatter | n > +2sigma | n < -2sigma | n |")
    P("|---|---|---|---|---:|---:|---:|---:|---:|")
    for (pr, mn, m_, lo_, hi_, se, sg, np_, nm_, n, _f) in S["rows_tab"] + S["rowsE"]:
        nd = 2 if mn != "maxd" else 1
        lab = {"rise": "rise 10-90", "w99": "w99", "fwhm": "FWHM", "tpk": "t_pk", "maxd": "max abs df (%)"}[mn]
        sgs = f"{sg:.2f}" + (" (null median)" if mn == "maxd" else "")
        exc = "-" if mn == "maxd" else f"{np.sqrt(max(_f ** 2 - sg ** 2, 0.0)):.2f}"
        P(f"| {pr} | {lab} | {fmt(m_, lo_, hi_, nd, sign=(mn != 'maxd'))} | {se:.3f} | {sgs} | {exc} | "
          f"{np_ if mn != 'maxd' else '-'} | {nm_ if mn != 'maxd' else '-'} | {n} |")
    P("\nExcess scatter = sqrt(((p84 - p16) / 2)^2 - noise sigma^2), ns (the module-dependent part). Expected from noise alone: 2.3 % of modules on each side of +-2 sigma. `(energy-top)` rows: pulses chosen as the bank does, "
      "the top-N mf_energy of the bin (sigma from the random-draw null).\n")
    P("Energy-top minus random draw, same bin and N, median [16, 84 %] (SE): rise / w99 / FWHM (ns)\n")
    P("| bin | median A (energy-top / random) | rise | w99 | FWHM |")
    P("|---|---|---|---|---|")
    for r in S["er_tab"]:
        P(f"| {r[0]} | {r[1]:.0f} / {r[2]:.0f} | " + " | ".join(f"{fmt(r[3 + 4 * k], r[4 + 4 * k], r[5 + 4 * k])} ({r[6 + 4 * k]:.3f})" for k in range(3)) + " |")
    P("")
    P("## 2. Amplitude slope (per ln A, per-module least squares over the bins)\n")
    P("| metric | random draw: median slope [16, 84 %] (SE) | lo->hi slope (median of per-module ratio) | energy-top: median slope [16, 84 %] (SE) | bank band templates (lo 100-400 / top 1000-2500), slope | n modules |")
    P("|---|---|---|---|---|---:|")
    for k, mn in enumerate(MET):
        bs = S["bank_slope"][mn]
        P(f"| {mn} | {fmt(*S['slope_pct'][k][[1, 0, 2]], nd=3)} ({S['slope_med'][k, 1]:.3f}) | "
          f"{S['slope_lohi'][k, 0]:+.3f} +- {S['slope_lohi'][k, 1]:.3f} | {fmt(*S['slope_pct_E'][k][[1, 0, 2]], nd=3)} ({S['slope_med_E'][k, 1]:.3f}) | "
          f"{bs[0]:+.3f} +- {bs[1]:.3f} | {int(np.isfinite(S['slopes'][:, k]).sum())} |")
    P("\nSlopes in ns per unit ln A.")
    P(f"\nBin median ln A: " + ", ".join(f"{bl[b]}: {S['lnA_med'][b]:.2f}" for b in range(nb)) + f"; lo->hi dlnA {S['dlnA_pair']:.2f}.\n")
    P("Bank band templates (shipped `shape_lo/mid/top`, <= 50 random pulses each, band_shape; medians of the knots A "
      + "/".join(f"{x:.0f}" for x in S["kA"]) + " ADC), difference median [16, 84 %] (rise / w99 / FWHM / t_pk, ns; n):\n")
    for (p, q), v in S["bank_band"].items():
        P(f"- {q} - {p}: " + " / ".join(fmt(*v[k][[1, 0, 2]], nd=2) for k in range(4)) + f"  (n {v[4]})")
    P("\n## 3. Hold-out lambda by bin (fraction of modules)\n")
    P("| bin | median A | lambda <= 10 | 100 | 1e3 | >= 1e4 | energy-top: <= 10 / 100 / 1e3 / >= 1e4 |")
    P("|---|---:|---:|---:|---:|---:|---|")
    for r, e in zip(S["lam_tab"], S["lam_E"]):
        P(f"| {r[0]} | {r[1]:.0f} | {r[2]:.2f} | {r[3]:.2f} | {r[4]:.2f} | {r[5]:.2f} | {e[1]:.2f} / {e[2]:.2f} / {e[3]:.2f} / {e[4]:.2f} |")
    P("\n## 4. Does ampcorr f(A) cover it?\n")
    v = S["val"]
    P(f"Proxy validation on {v['n']} shipped band-template sets (lo 100-400 / mid / top): r_b / r_mod from the stored C++ r_lo, r_mid, r_top: "
      f"lo {v['cpp_lo'][0]:.4f}, mid/lo {v['cpp_mid_lo'][0]:.4f}, top/mid {v['cpp_top_mid'][0]:.4f}; analytic-integral proxy "
      f"lo {v['lo'][0]:.4f}, mid/lo {v['mid_lo'][0]:.4f}, top/mid {v['top_mid'][0]:.4f}; ampcorr table steps "
      f"{v['tab_mid_lo']:.4f}, {v['tab_top_mid']:.4f}. Per-module proxy - C++: lo {v['diff_lo'][0]:+.4f}, mid {v['diff_mid'][0]:+.4f}, top {v['diff_top'][0]:+.4f}.\n")
    P("| bin | median A (random / energy-top) | proxy r_bin / r_mod, random | proxy, energy-top | table f(A) at the two A |")
    P("|---|---|---|---|---|")
    for r, e in zip(S["fa_tab"], S["fa_tabE"]):
        P(f"| {r[0]} | {r[1]:.0f} / {e[1]:.0f} | {r[2]:.4f} +- {r[3]:.4f} | {e[2]:.4f} +- {e[3]:.4f} | {r[4]:.4f} / {e[4]:.4f} |")
    P("\n| step | proxy ratio, random | table f ratio at the bin A | proxy ratio, energy-top | table f ratio at the bin A |\n|---|---|---:|---|---:|")
    for r in S["fstep"]:
        P(f"| {r[0]} | {r[1]:.4f} +- {r[2]:.4f} | {r[3]:.4f} | {r[4]:.4f} +- {r[5]:.4f} | {r[6]:.4f} |")
    at = S["area_tab"]
    P(f"\nUnit-peak template area, {bl[S['pairs'][-1][1]]} / {bl[S['pairs'][-1][0]]}, median +- SE over modules (core 0-50 ns / tail 50 ns-end / total): "
      f"random draws {at['R'][0][0]:.4f} +- {at['R'][0][1]:.4f} / {at['R'][1][0]:.4f} +- {at['R'][1][1]:.4f} / {at['R'][2][0]:.4f} +- {at['R'][2][1]:.4f}; "
      f"energy-top {at['E'][0][0]:.4f} +- {at['E'][0][1]:.4f} / {at['E'][1][0]:.4f} +- {at['E'][1][1]:.4f} / {at['E'][2][0]:.4f} +- {at['E'][2][1]:.4f}.")
    P("\n## 5. Ta and 0.73 GeV against 2.24 GeV H2 (shipped best-50 banks)\n")
    P(f"Pool mix (shipped pools, modules pooled): mean ln A / median A / fraction < 600 ADC: "
      + "; ".join(f"{k} {S['mean_lnA'][k]:.2f} / {S['med_amp'][k]:.0f} / {S['mean_flo'][k]:.2f}" for k in S["mean_lnA"]) + ".\n")
    for c in (S["cmp_ta"], S["cmp_07"]):
        P(f"### {c['tag']} ({len(c['names'])} modules)\n")
        P("| metric | observed median [16, 84 %] | predicted, random-draw slope | predicted, energy-top slope | observed - predicted (random) | Spearman vs d<ln A> | vs d(frac < 600) | vs frac < 600 of the first bank | Theil-Sen slope of obs on d<ln A> [95 %] vs bins (random / energy-top) |")
        P("|---|---|---|---|---|---|---|---|---|")
        for k, mn in enumerate(MET[:3]):
            pr, co, ft = c["pred"][mn], c["corr"][mn], c["fit"][mn]
            P(f"| {mn} | {fmt(*pr['obs'][[1, 0, 2]])} | {fmt(*pr['pred'][[1, 0, 2]])} | {fmt(*pr['predE'][[1, 0, 2]])} | {fmt(*pr['resid'][[1, 0, 2]])} | "
              f"rho {co['lnA'][0]:+.2f} (p {co['lnA'][1]:.1e}) | rho {co['lo'][0]:+.2f} (p {co['lo'][1]:.1e}) | "
              f"rho {co['flo'][0]:+.2f} (p {co['flo'][1]:.1e}) | {ft['slope']:+.3f} [{ft['lo']:+.3f}, {ft['hi']:+.3f}] vs {S['slope_med'][k, 0]:+.3f} / {S['slope_med_E'][k, 0]:+.3f} |")
        P("")
    P("Per-module prediction (each module's own amplitude slope x its d<ln A>) against the observed difference, modules with both: "
      "Spearman rho (p), Theil-Sen slope of observed on predicted [95 %] (1 = fully explained), robust sigma of predicted / observed (ns), for random-draw / energy-top slopes.\n")
    P("| comparison | metric | n | random: rho (p) | slope [95 %] | sigma pred / obs | energy-top: rho (p) | slope [95 %] | sigma pred / obs |")
    P("|---|---|---:|---|---|---|---|---|---|")
    for c in (S["cmp_ta"], S["cmp_07"]):
        for mn in MET[:3]:
            r_, e_ = c["permod"][mn]["R"], c["permod"][mn]["E"]
            P(f"| {c['tag']} | {mn} | {r_['n']} | {r_['rho']:+.2f} ({r_['p']:.1e}) | {r_['slope']:+.2f} [{r_['lo']:+.2f}, {r_['hi']:+.2f}] | "
              f"{r_['sd_pred']:.2f} / {r_['sd_obs']:.2f} | {e_['rho']:+.2f} ({e_['p']:.1e}) | {e_['slope']:+.2f} [{e_['lo']:+.2f}, {e_['hi']:+.2f}] | "
              f"{e_['sd_pred']:.2f} / {e_['sd_obs']:.2f} |")
    P("")
    P("Run-type averages (`shape_avg`) rise / w99 / FWHM (ns): " + "; ".join(
        f"{k}: {v_[0]:.2f} / {v_[1]:.2f} / {v_[2]:.2f}" for k, v_ in S["avg"].items()) + ".\n")
    d73 = S["avg"]["0.73 H2"] - S["avg"]["2.24 H2"]
    dTa = S["avg"]["Ta"] - S["avg"]["2.24 H2"]
    dl73 = S["mean_lnA"]["0.73 H2"] - S["mean_lnA"]["2.24 H2"]
    dlTa = S["mean_lnA"]["Ta"] - S["mean_lnA"]["2.24 H2"]
    P("Run-type average offsets vs 2.24 H2 (rise / w99 / FWHM, ns): observed 0.73 " + " / ".join(f"{x:+.2f}" for x in d73[:3])
      + "; predicted (median slope x d<ln A> of the pool means, d<ln A> = " + f"{dl73:+.2f}): "
      + " / ".join(f"{S['slope_med'][k, 0] * dl73:+.2f}" for k in range(3)) + " (energy-top slope: "
      + " / ".join(f"{S['slope_med_E'][k, 0] * dl73:+.2f}" for k in range(3))
      + "). Ta observed " + " / ".join(f"{x:+.2f}" for x in dTa[:3]) + f"; predicted (d<ln A> = {dlTa:+.2f}): "
      + " / ".join(f"{S['slope_med'][k, 0] * dlTa:+.2f}" for k in range(3)) + " (energy-top: "
      + " / ".join(f"{S['slope_med_E'][k, 0] * dlTa:+.2f}" for k in range(3)) + ").\n")
    (out / "report.md").write_text("\n".join(L) + "\n")

    # per-module TSV
    hdr = ["module", "midx"]
    for b in range(nb):
        hdr += [f"A_{bl[b]}", f"lam_{bl[b]}"] + [f"{m}_{bl[b]}" for m in MET] + [f"E_{m}_{bl[b]}" for m in MET]
    hdr += [f"slope_{m}" for m in MET] + [f"slopeE_{m}" for m in MET] + [f"rproxy_{bl[b]}" for b in range(nb)] \
        + [f"rproxyE_{bl[b]}" for b in range(nb)]
    hdr += [f"null_d{m}_{bl[b]}" for b in range(nb) for m in MET + ["maxd"]]
    with open(out / "modules.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for m, d in enumerate(rec):
            row = [d["name"], d["midx"]]
            for b in range(nb):
                row += [d["A"][b], d["lam"][b]] + list(d["m"][b]) + list(d["mE"][b])
            row += list(S["slopes"][m]) + list(S["slopes_E"][m]) + list(S["g_bins"][m]) + list(S["g_binsE"][m])
            for b in range(nb):
                row += list(d["null"][b]) if d["null"][b] is not None else [np.nan] * 5
            fh.write("\t".join(f"{x:.6g}" if isinstance(x, (float, np.floating)) else str(x) for x in row) + "\n")
    for tag, c in (("ta_h2", S["cmp_ta"]), ("h073_h2", S["cmp_07"])):
        with open(out / f"{tag}_modules.tsv", "w") as fh:
            fh.write("\t".join(["module", "d_rise", "d_w99", "d_fwhm", "d_tpk", "d_meanlnA", "d_frac_lt600", "frac_lt600_first", "medA_first"]) + "\n")
            for i, n in enumerate(c["names"]):
                fh.write("\t".join([n] + [f"{x:.6g}" for x in list(c["dobs"][i]) + [c["dln"][i], c["dlo"][i], c["flo"][i], c["medA"][i]]]) + "\n")
    figure(out, S)


def figure(out, S):
    bl, nb, M, rec = S["bl"], S["nb"], S["M"], S["rec"]
    fig, ax = plt.subplots(3, 3, figsize=(16.5, 13.5))
    col = ["#4c78a8", "#f58518", "#54a24b"]
    ci, cj = S["pairs"][-1]
    # (a) mean shape difference hi - lo
    ax0 = ax[0, 0]
    sh = np.array([[d["shape"][b] for b in range(nb)] for d in rec], np.float64)      # [mod, bin, t]
    sh = sh / np.nanmax(sh, axis=2, keepdims=True)
    win = (SE.T_GRID > -3) & (SE.T_GRID < 45)
    t = SE.T_GRID[win]
    for (i, j), c in zip(S["pairs"], ["#f58518", "#54a24b", "#e45756"] if len(S["pairs"]) == 3 else ["#e45756"]):
        d = sh[:, j, :][:, win] - sh[:, i, :][:, win]
        m = np.nanmedian(d, axis=0) * 100
        lo, hi = np.nanpercentile(d, [16, 84], axis=0) * 100
        ax0.plot(t, m, color=c, label=f"{bl[j]} - {bl[i]}")
        if (i, j) == (ci, cj):
            ax0.fill_between(t, lo, hi, color=c, alpha=0.15, lw=0)
    ax0.axhline(0, color="k", lw=0.5)
    ax0.set(xlabel="t after onset (ns)", ylabel="normalised shape difference (%)", title="(a) median over modules, band 16-84 % (hi - lo)")
    ax0.legend(fontsize=8)
    # (b) metric vs ln A (relative to lowest bin)
    ax1 = ax[0, 1]
    xs = S["lnA_med"]
    for k, (mn, c) in enumerate(zip(MET[:3], ["#4c78a8", "#e45756", "#54a24b"])):
        y = M[:, :, k] - M[:, [0], k]
        mm = np.nanmedian(y, axis=0)
        se = 1.2533 * np.array([rsig(y[:, b]) for b in range(nb)]) / np.sqrt(M.shape[0])
        ax1.errorbar(np.exp(xs), mm, se, marker="o", color=c, label=mn)
    ax1.axhline(0, color="k", lw=0.5)
    ax1.set_xscale("log")
    ax1.set(xlabel="bin median amplitude (ADC)", ylabel="metric - lowest bin (ns)", title="(b) amplitude trend, median +- SE")
    ax1.legend(fontsize=8)
    # (c) max|df|
    ax2 = ax[0, 2]
    dm = S["pair_res"][(ci, cj)]["dmax"]
    nl = S["null_all"][ci][:, 4] if len(S["null_all"][ci]) else np.array([])
    bins = np.linspace(0, max(np.nanpercentile(dm, 99), 1), 30)
    ax2.hist(dm, bins, histtype="step", color="#e45756", label=f"{bl[cj]} vs {bl[ci]}")
    if nl.size:
        ax2.hist(nl, bins, histtype="step", color="grey", label=f"split halves of {bl[ci]}")
    ax2.set(xlabel="max abs df (%)", ylabel="modules", title="(c) max |df|, normalised shapes")
    ax2.legend(fontsize=8)
    # (d-f) metric differences obs vs null
    for k, (mn, a_) in enumerate(zip(MET[:3], ax[1])):
        x = S["pair_res"][(ci, cj)]["dM"][:, k]
        sp = S["pair_res"][(ci, cj)]["sp"][k]
        lim = 4.5 * max(sp, np.nanstd(x) * 0.5)
        bins = np.linspace(-lim, lim, 41)
        a_.hist(x, bins, color="#e45756", alpha=0.6, label=f"{bl[cj]} - {bl[ci]} (n {np.isfinite(x).sum()})")
        nn = S["null_all"][ci]
        if len(nn):
            sc = sp / S["sig_b"][ci, k]
            a_.hist(nn[:, k] * sc, bins, histtype="step", color="k", label="split-half noise (same N)")
        a_.axvline(np.nanmedian(x), color="#e45756", lw=1)
        a_.axvline(0, color="k", lw=0.5)
        a_.set(xlabel=f"delta {mn} (ns)", ylabel="modules", title=f"({'def'[k]}) {mn}: median {np.nanmedian(x):+.2f}, sigma_noise {sp:.2f}")
        a_.legend(fontsize=8)
    # (g) f(A)
    ax7 = ax[2, 0]
    Abm = np.array([r[1] for r in S["fa_tab"]])
    fm = np.array([r[2] for r in S["fa_tab"]])
    fe = np.array([r[3] for r in S["fa_tab"]])
    ft = np.array([r[4] for r in S["fa_tab"]])
    norm = fm[-1] / ft[-1]
    ax7.errorbar(Abm, fm / norm, fe, marker="o", color="#4c78a8", label="bin templates (proxy r_bin / r_mod), scaled to f(A) at top")
    xg = np.linspace(np.log(Abm[0]), np.log(Abm[-1]), 50)
    ax7.plot(np.exp(xg), np.interp(xg, np.log([S["knots"][k]["A"] for k in ("lo", "mid", "top")]),
                                   [S["knots"][k]["f"] for k in ("lo", "mid", "top")]),
             color="k", label="ampcorr f(A) (knot interpolation)")
    ax7.set_xscale("log")
    ax7.set(xlabel="amplitude (ADC)", ylabel="integral factor", title="(g) integral factor: shape-implied vs f(A)")
    ax7.legend(fontsize=7)
    # (h, i) Ta and 0.73
    for a_, c, lab in ((ax[2, 1], S["cmp_ta"], "Ta - 2.24 H2"), (ax[2, 2], S["cmp_07"], "0.73 - 2.24 H2")):
        x, y = c["dln"], c["dobs"][:, 0]
        fin = np.isfinite(y)
        a_.scatter(x[fin], y[fin], s=5, alpha=0.4, color="#4c78a8")
        xx = np.linspace(np.nanpercentile(x, 1), np.nanpercentile(x, 99), 10)
        a_.plot(xx, S["slope_med"][0, 0] * xx, color="#e45756", label=f"bins slope {S['slope_med'][0, 0]:+.2f} ns / ln A")
        ft_ = c["fit"]["rise"]
        a_.plot(xx, np.median(y[fin] - ft_["slope"] * x[fin]) + ft_["slope"] * xx, color="k", label=f"Theil-Sen {ft_['slope']:+.2f}")
        a_.axhline(0, color="k", lw=0.4)
        co = c["corr"]["rise"]["lnA"]
        a_.set(xlabel="d<ln A> of the shipped pools", ylabel="d rise (ns)", ylim=(-1.2, 1.2),
               title=f"({'hi'[0 if 'Ta' in lab else 1]}) {lab}: rho {co[0]:+.2f} (p {co[1]:.0e})")
        a_.legend(fontsize=8)
    fig.suptitle(f"PbWO4 shape vs amplitude, 2.24 GeV H2 select pools (run {S['run']}, {S['nrec']} modules, N {S['npl']} per bin)", y=0.995)
    fig.tight_layout()
    fig.savefig(out / "amp_binned.pdf")
    fig.savefig(out / "amp_binned.png", dpi=110)
    plt.close(fig)


if __name__ == "__main__":
    main()
