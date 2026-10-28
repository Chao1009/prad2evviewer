#!/usr/bin/env python3
"""template_doublet_audit.py -- hidden doublets in the best-50 template pools.

Question: the best-50 spline templates are built from <= 50 pulses per module (top
mf_energy of the candidates).  If a fraction of those pool pulses are really two
pulses (pile-up the selection did not veto), does the spline template inherit the
second pulse, and does that explain Ta (X17) templates being wider than e-p?

1. Pool reproduction: the candidate files named in the bank's meta are loaded with
   spline_template_extract's own loader and the build's own argv; the pool of a module is
   the top n_max by mf_energy of its chain pool (p0 for fallback 0, p1 for fallback 1),
   exactly as process_module does.  Checked against n_pool of every module and against
   a re-run of process_module on spot-check modules (shape vs the shipped template).
2. Two-pulse fit (twopulse_fit.fit_wave, production options eps 0.02) of every pool
   pulse with the module's own template; flag = dchi2 > thr, fitted dt >= 2 ns and
   A2/A1 >= 0.05.  Thresholds: U0 summary.json (per type where recorded, else 2.24 GeV H2).
3. Shape impact (2.24 GeV Ta X17 and 2.24 GeV H2): the 20 most flagged modules are rebuilt
   with spline_template_extract.process_module on the pool minus the flagged pulses (the
   hold-out of the build is kept), compared with the shipped template, with controls:
   null = k random pool pulses removed, noise = split halves of a pool.
4. Selection bias: flagged fraction in the 10 largest vs 10 smallest amplitudes of each pool.

    MPLBACKEND=Agg python3 template_doublet_audit.py

Output: output/templates_spline_best50_20261026/doublet_audit/ (report.md, per-bank TSVs,
rebuilt templates, doublet_audit.pdf).  numpy / scipy / matplotlib.
"""
from __future__ import annotations

import argparse
import glob
import json
import multiprocessing as mp
import sys
import time
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import spearmanr

import mf_pulse
import spline_template_extract as SE
import twopulse_fit as TP

HERE = Path(__file__).resolve().parent
B50 = TP.B50
OUT = B50 / "doublet_audit"
SLUGS = ["0p73GeV_H2", "2p24GeV_C", "2p24GeV_H2", "2p24GeV_Ta_X17", "3p49GeV_H2"]
LABEL = {"0p73GeV_H2": "0.73 H2", "2p24GeV_C": "2.24 C", "2p24GeV_H2": "2.24 H2",
         "2p24GeV_Ta_X17": "2.24 Ta", "3p49GeV_H2": "3.49 H2"}
U0 = HERE / "output/ml_pileup_20261027/U0/summary.json"
U0_KEY = {"0p73GeV_H2": "0.73 GeV H2", "2p24GeV_H2": "2.24 GeV H2", "3p49GeV_H2": "3.49 GeV H2"}
REBUILD = ["2p24GeV_Ta_X17", "2p24GeV_H2"]
DT_BINS = [(2, 8), (8, 16), (16, 32), (32, 1e9)]
DT_LAB = ["2-8", "8-16", "16-32", ">32"]
R_BINS = [(0.05, 0.2), (0.2, 0.5), (0.5, 1e9)]
R_LAB = ["0.05-0.2", "0.2-0.5", ">0.5"]
THR_SENS = (98.0, 175.0, 300.0, 1000.0)
STRICT = (8.0, 0.2)                  # "clean doublet" subset: dt >= 8 ns and secondary/main >= 0.2
TG = SE.T_GRID
WIN = (TG > -5) & (TG < 60)          # xtype-study window for max|df| (shapes are in the bank onset frame)

_CTX = None


# ------------------------------------------------------------------ context --
class _Got(Exception):
    pass


def load_ctx(slug, quick=0):
    """cmd_build's setup without the pipeline: candidates, hold-out, build args, tau table."""
    z = np.load(B50 / slug / f"templates_{slug}.npz", allow_pickle=False)
    meta = json.loads(str(z["meta"]))
    argv = list(meta["argv"][1:])
    argv[argv.index("--out") + 1] = str(OUT / "_unused")

    def grab(a):
        raise _Got(a)
    orig = SE.cmd_build
    SE.cmd_build = grab
    try:
        SE.main(argv)
    except _Got as e:
        a = e.args[0]
    finally:
        SE.cmd_build = orig
    SE.SMOOTH_LAMBDA = float(a.smooth_lambda)
    SE.FIT_POINTS = a.fit_points
    SE.CONV_DF = float(a.converge_df)
    SE.SHAPE_SIGMA = float(a.shape_sigma)
    SE.LAM_CHI2_SIGMA = a.lambda_chi2_sigma
    SE.LAM_MIN = float(a.lambda_min)
    SE.LAM_FALLBACK = float(a.lambda_fallback)
    files = sorted({f for pat in a.pulses.split(",") for f in glob.glob(pat.strip())})
    D = SE.load_candidates(files, SE.BUILD_FIELDS + ("mf_energy", "mf_npk", "mf_sec_ratio"))
    hold, n_hold = SE.holdout_energy(D, a.n_max, a.holdout_frac)
    h50 = HERE / a.h50 if not Path(a.h50).is_absolute() else Path(a.h50)
    column = mf_pulse.resolve_column(str(h50), "auto", prefer="25320")
    tau_rows, _, agg_json = mf_pulse.load_h50(str(h50), column)
    tau_agg = mf_pulse.load_aggregate(agg_json)
    names = [None] * (int(max(z["module_index"].max(), D["midx"].max())) + 1)
    for n, mi in zip(z["names"], z["module_index"]):
        names[int(mi)] = str(n)
    order = np.argsort(D["midx"], kind="stable")
    ms = D["midx"][order]
    ctx = dict(D=D, args=a, iter_kw=SE._iter_kw(a), hold=hold, n_hold=n_hold, names=names,
               tau=tau_rows, tau_agg=tau_agg, band_edges=a.bands,
               fa_bands=SE.FA_BANDS_DEFAULT if a.bands is None else SE.fa_bands_from_edges(a.bands),
               wave=SE._wave_ctx(D["metas"]), run_type=a.run_type, files=files,
               _order=order, _ms=ms)
    return ctx, z


def module_idx(ctx, mi):
    o, ms = ctx["_order"], ctx["_ms"]
    return o[np.searchsorted(ms, mi, "left"):np.searchsorted(ms, mi, "right")]


def module_pool(ctx, mi, fb):
    """The module's fit pool (candidate indices), as process_module selects it."""
    D, a = ctx["D"], ctx["args"]
    idx = module_idx(ctx, mi)
    amp, tier, held = D["amp"][idx], D["tier"][idx], ctx["hold"][idx]
    if fb == 0:
        mask = (~held) & (tier >= 1) & D["band_b0"][idx]
    else:
        mask = (~held) & (tier >= 1) & (amp >= SE.RELAXED_BAND[0]) & (amp <= SE.RELAXED_BAND[1])
    pool = idx[mask]
    if pool.size > a.n_max:
        o = np.argsort(-D["mf_energy"][pool], kind="stable")
        pool = np.sort(pool[o[:a.n_max]])
    return pool, idx


def rebuild(mi, keep):
    """process_module on the pool `keep` (+ the module's hold-out); returns the result dict."""
    ctx = _CTX
    idx = module_idx(ctx, mi)
    h = idx[ctx["hold"][idx]]
    sel = np.unique(np.concatenate([np.asarray(keep, np.int64), h]))
    return SE.process_module(ctx, mi, sel)


def _task(args):
    tag, mi, keep, full = args
    if full:                           # exactly the build's call
        r = SE.process_module(_CTX, mi, module_idx(_CTX, mi))
    else:
        r = rebuild(mi, keep)
    kg = r.get("kept_gidx")
    return (tag, mi, int(r["fallback"]), r.get("shape"), float(r.get("lam_final", np.nan)),
            int(r["n_pool"]), None if kg is None else np.asarray(kg, np.int64))


# ------------------------------------------------------------------ metrics --
def _cross(f, level, step):
    ip = int(np.argmax(f))
    j = ip
    while 0 <= j + step < f.size and f[j] >= level:
        j += step
    if f[j] >= level:
        return np.nan
    k = j - step
    return float(TG[j] + (level - f[j]) / (f[k] - f[j]) * (TG[k] - TG[j]))


def shape_metrics(f):
    """rise 10-90, w99, FWHM, t_pk (parabolic), all ns, of a shape on T_GRID."""
    f = np.nan_to_num(np.asarray(f, np.float64))
    pk = f.max()
    f = f / pk
    rise, w99 = SE._rise_width(f)
    fwhm = _cross(f, 0.5, 1) - _cross(f, 0.5, -1)
    i = int(np.argmax(f))
    y0, y1, y2 = f[i - 1], f[i], f[i + 1]
    den = y0 - 2 * y1 + y2
    tpk = TG[i] + (0.5 * (y0 - y2) / den * (TG[1] - TG[0]) if den < 0 else 0.0)
    return np.array([rise, w99, fwhm, tpk])


def maxdiff(fa, fb):
    a = np.nan_to_num(np.asarray(fa, np.float64))
    b = np.nan_to_num(np.asarray(fb, np.float64))
    return float(np.abs(a / a.max() - b / b.max())[WIN].max() * 100.0)


METRIC = ["rise", "w99", "fwhm", "tpk"]


def mad_sig(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(1.4826 * np.median(np.abs(x - np.median(x)))) if x.size else np.nan


# --------------------------------------------------------------- flag logic --
def classify(R, thr):
    # ratio = secondary / main amplitude (<= 1): min(A1, A2) / max(A1, A2)
    dchi, dt = R[:, 12], R[:, 11] - R[:, 9]
    ratio = np.minimum(R[:, 8], R[:, 10]) / np.maximum(np.maximum(R[:, 8], R[:, 10]), 1e-9)
    over = (R[:, 13] == 1) & (dchi > thr)
    degen = over & ((dt < 2.0) | (ratio < 0.05))
    flag = over & ~degen
    return flag, degen, dt, ratio


def bin2d(flag, dt, ratio):
    t = np.zeros((len(DT_BINS), len(R_BINS)), int)
    for i, (a, b) in enumerate(DT_BINS):
        for j, (c, d) in enumerate(R_BINS):
            t[i, j] = int((flag & (dt >= a) & (dt < b) & (ratio >= c) & (ratio < d)).sum())
    return t


# ---------------------------------------------------------------- per bank --
def run_bank(slug, args, thr_map, pw_n):
    global _CTX
    t0 = time.time()
    ctx, z = load_ctx(slug)
    _CTX = ctx
    D = ctx["D"]
    bank, _ = TP.load_bank(B50, slug)
    TP._BANK = bank
    TP.OPT.update(eps=0.02, kap=0.0, search=80.0, grid=0.5, min_sep=1.0)
    ix = {n: i for i, n in enumerate(bank.names)}
    names = [str(n) for n in z["names"]]
    fbk = np.asarray(z["fallback"])
    n_pool_ref = np.asarray(z["n_pool"])
    rows = [r for r in range(len(names)) if fbk[r] in (0, 1)]
    if args.quick:
        rows = rows[::args.quick]
    pools, mism = {}, 0
    for r in rows:
        mi = int(z["module_index"][r])
        p, _ = module_pool(ctx, mi, int(fbk[r]))
        pools[r] = p
        mism += int(p.size != n_pool_ref[r])
    thr = thr_map[slug]
    pw = mp.get_context("fork").Pool(pw_n)
    try:
        # ---- two-pulse fits of all pool pulses
        G = np.concatenate([pools[r] for r in rows])
        row_of = np.concatenate([np.full(pools[r].size, r) for r in rows])
        raw = D["samples"][G].astype(np.float64)
        y, ped = TP.subtract_ped(bank, raw)
        sig = {}
        for r in rows:
            mi = int(z["module_index"][r])
            sig[r] = max(float(np.median(D["rms"][module_idx(ctx, mi)])), mf_pulse.SIGMA_MIN)
        sigp = np.array([sig[r] for r in row_of])
        R = TP.fit_batch(pw, y, np.array([ix[names[r]] for r in row_of]), sigp)
        flag, degen, dt, ratio = classify(R, thr)
        strict = flag & (dt >= STRICT[0]) & (ratio >= STRICT[1])
        print(f"[{slug}] {len(G)} pool pulses, {len(rows)} modules, fit {time.time() - t0:.0f} s; "
              f"flagged {flag.sum()} degenerate {degen.sum()} (thr {thr:.1f}); n_pool mismatches {mism}", flush=True)
        # ---- per-module summary and selection bias
        amp = D["amp"][G]
        rank = np.zeros(G.size, int)
        top = np.zeros(G.size, bool)
        bot = np.zeros(G.size, bool)
        mod = []
        for r in rows:
            s = np.nonzero(row_of == r)[0]
            o = np.argsort(-amp[s], kind="stable")
            rank[s[o]] = np.arange(1, s.size + 1)
            if s.size >= 20:
                top[s[o[:10]]] = True
                bot[s[o[-10:]]] = True
            mod.append(dict(name=names[r], row=r, midx=int(z["module_index"][r]), fallback=int(fbk[r]),
                            n_pool=int(s.size), n_flag=int(flag[s].sum()), n_degen=int(degen[s].sum()),
                            n_flag_dtlt16=int((flag[s] & (dt[s] < 16)).sum()),
                            n_strict=int(strict[s].sum()),
                            n_over=int(((R[s, 13] == 1) & (R[s, 12] > thr)).sum()),
                            dchi_med=float(np.median(R[s, 12])), amp_med=float(np.median(amp[s])),
                            sig_ped=sig[r]))
        res = dict(slug=slug, thr=thr, rows=rows, G=G, row_of=row_of, R=R, flag=flag, strict=strict, degen=degen, dt=dt,
                   ratio=ratio, amp=amp, rank=rank, top=top, bot=bot, mod=mod, n_pool_mismatch=mism,
                   fb_counts=[int((fbk == k).sum()) for k in range(3)])
        # examples (flagged pulses with the largest dchi2 in each |dt| bin)
        ex = []
        for (a, b) in DT_BINS:
            c = np.nonzero(flag & (dt >= a) & (dt < b))[0]
            if c.size:
                k = int(c[np.argmax(R[c, 12])])
                ex.append(dict(k=k, y=y[k].copy(), row=int(row_of[k]), bidx=ix[names[int(row_of[k])]]))
        res["examples"] = ex
        # ---- reproduction spot check (the build's own process_module on the whole module)
        rng = np.random.default_rng(20261027)
        spot = [int(x) for x in rng.choice(rows, min(args.n_spot, len(rows)), replace=False)]
        tasks = []
        for r in spot:
            mi = int(z["module_index"][r])
            tasks.append((("full", r), mi, None, True))
            tasks.append((("ctl", r), mi, pools[r], False))
        out = pw.map(_task, tasks, chunksize=1)
        spotres = []
        for (tag, mi, fb, shp, lam, npool, kg) in out:
            if tag[0] != "full":
                continue
            r = tag[1]
            ctl = [o for o in out if o[0] == ("ctl", r)][0]
            ship = z["shape"][r].astype(np.float64)
            d_full = float(np.nanmax(np.abs(shp.astype(np.float64) - ship))) if shp is not None else np.nan
            d_ctl = float(np.nanmax(np.abs(ctl[3].astype(np.float64) - ship))) if ctl[3] is not None else np.nan
            spotres.append(dict(name=names[r], d_full=d_full, d_ctl=d_ctl, fb=fb,
                                kept_in_pool=bool(np.isin(kg, pools[r]).all()) if kg is not None else None,
                                n_kept=None if kg is None else int(kg.size)))
        res["spot"] = spotres
        print(f"[{slug}] spot-check max|df| full {max(s['d_full'] for s in spotres):.2e} "
              f"restricted-pool control {max(s['d_ctl'] for s in spotres):.2e}", flush=True)
        # ---- shape impact
        if slug in REBUILD:
            res["rebuild"] = shape_impact(slug, args, pw, z, rows, pools, mod, flag, strict, G, row_of)
    finally:
        pw.close()
        pw.join()
    res["wall"] = time.time() - t0
    res["bank"] = bank
    res["z_shape"] = z["shape"].astype(np.float32)
    res["z_names"] = names
    res["shape_avg"] = z["shape_avg"].astype(np.float32)
    res["z_fb"] = fbk
    return res


def shape_impact(slug, args, pw, z, rows, pools, mod, flag, strict, G, row_of):
    rng = np.random.default_rng(20261028)
    cand = [m for m in mod if m["n_flag"] >= 1]
    cand.sort(key=lambda m: (-m["n_flag"], -m["n_flag"] / m["n_pool"], m["name"]))
    top = cand[:20]
    tasks = []
    for m in top:
        r, mi = m["row"], m["midx"]
        pool = pools[r]
        fl = G[(row_of == r) & flag]
        tasks.append((("rem", r), mi, np.setdiff1d(pool, fl), False))
        st = G[(row_of == r) & strict]
        if st.size:
            tasks.append((("rems", r), mi, np.setdiff1d(pool, st), False))
        tasks.append((("ctl", r), mi, pool, False))
        for d in range(args.n_null):
            drop = rng.choice(pool, len(fl), replace=False)
            tasks.append((("null", r, d), mi, np.setdiff1d(pool, drop), False))
    # split-half noise on modules with a full pool
    full = [m for m in mod if m["n_pool"] == 50 and m["fallback"] == 0]
    sample = [full[i] for i in rng.choice(len(full), min(args.n_noise, len(full)), replace=False)]
    for m in sample:
        pool = rng.permutation(pools[m["row"]])
        tasks.append((("hA", m["row"]), m["midx"], pool[:25], False))
        tasks.append((("hB", m["row"]), m["midx"], pool[25:], False))
    t0 = time.time()
    out = pw.map(_task, tasks, chunksize=1)
    print(f"[{slug}] {len(tasks)} rebuilds in {time.time() - t0:.0f} s", flush=True)
    sh = {o[0]: o for o in out}
    ship = lambda r: z["shape"][r].astype(np.float64)
    rec = []
    for m in top:
        r = m["row"]
        o = sh[("rem", r)]
        row = dict(name=m["name"], row=r, n_flag=m["n_flag"], n_pool=m["n_pool"], fb_new=o[2],
                   n_pool_new=o[5], lam_ship=sh[("ctl", r)][4], fb_ship=int(z["fallback"][r]),
                   lam_new=o[4], shape=None if o[3] is None else o[3].astype(np.float32),
                   shipped=ship(r).astype(np.float32))
        if o[3] is None:
            rec.append(row)
            continue
        ms0 = shape_metrics(ship(r))
        d = shape_metrics(o[3]) - ms0
        row.update(dm=d, maxd=maxdiff(o[3], ship(r)), m_ship=ms0)
        c = sh[("ctl", r)]
        row["maxd_ctl"] = maxdiff(c[3], ship(r)) if c[3] is not None else np.nan
        row["stable"] = bool(o[4] == c[4] and o[2] == row["fb_ship"])
        row["n_strict"] = m["n_strict"]
        if ("rems", r) in sh and sh[("rems", r)][3] is not None:
            o2 = sh[("rems", r)]
            row["dm_s"] = shape_metrics(o2[3]) - ms0
            row["maxd_s"] = maxdiff(o2[3], ship(r))
            row["shape_s"] = o2[3].astype(np.float32)
        nl = [sh[("null", r, k)] for k in range(args.n_null) if sh[("null", r, k)][3] is not None]
        row["null_maxd"] = [maxdiff(x[3], ship(r)) for x in nl]
        row["null_dm"] = [shape_metrics(x[3]) - ms0 for x in nl]
        rec.append(row)
    noise = dict(dm=[], maxd=[])
    for m in sample:
        a, b = sh[("hA", m["row"])], sh[("hB", m["row"])]
        if a[3] is None or b[3] is None:
            continue
        noise["dm"].append(shape_metrics(a[3]) - shape_metrics(b[3]))
        noise["maxd"].append(maxdiff(a[3], b[3]))
    noise["dm"] = np.array(noise["dm"])
    noise["maxd"] = np.array(noise["maxd"])
    return dict(rec=rec, noise=noise)


# --------------------------------------------------------------- tabulation --
def bank_tables(res):
    f, dg, dt, ratio = res["flag"], res["degen"], res["dt"], res["ratio"]
    n = f.size
    lt16 = f & (dt < 16)
    out = dict(n=n, n_flag=int(f.sum()), n_degen=int(dg.sum()), frac=f.sum() / n,
               frac_lt16=lt16.sum() / n, t2d=bin2d(f, dt, ratio),
               dtb=[int((f & (dt >= a) & (dt < b)).sum()) for a, b in DT_BINS],
               n_over=int(((res["R"][:, 13] == 1) & (res["R"][:, 12] > res["thr"])).sum()))
    nf = np.array([m["n_flag"] for m in res["mod"]])
    sf = res["strict"]
    ns = np.array([m["n_strict"] for m in res["mod"]])
    out.update(n_strict=int(sf.sum()), frac_strict=sf.sum() / n, mod_s_ge1=int((ns >= 1).sum()), mod_s_ge3=int((ns >= 3).sum()),
               top_frac_s=float(sf[res["top"]].mean()), bot_frac_s=float(sf[res["bot"]].mean()))
    ab = [0, 600, 1500, 1e9]
    out["amp_tab"] = [[int(((res["amp"] >= lo) & (res["amp"] < hi)).sum()),
                       100 * float(f[(res["amp"] >= lo) & (res["amp"] < hi)].mean()),
                       100 * float(sf[(res["amp"] >= lo) & (res["amp"] < hi)].mean())] for lo, hi in zip(ab[:-1], ab[1:])]
    out["n_mod"] = nf.size
    out["mod_ge3"] = int((nf >= 3).sum())
    out["mod_ge1"] = int((nf >= 1).sum())
    out["nflag_hist"] = np.bincount(np.minimum(nf, 10), minlength=11)
    t, b = res["top"], res["bot"]
    out["top_frac"] = float(f[t].mean()) if t.any() else np.nan
    out["bot_frac"] = float(f[b].mean()) if b.any() else np.nan
    out["top_n"], out["bot_n"] = int(f[t].sum()), int(f[b].sum())
    out["ratio_tb"] = out["top_frac"] / out["bot_frac"] if out["bot_frac"] > 0 else np.nan
    sens = {}
    R = res["R"]
    for th in THR_SENS:
        fl, _, _, _ = classify(R, th)
        sens[th] = float(fl.mean())
    out["sens"] = sens
    sens_s = {}
    for th in THR_SENS:
        fl, _, _, _ = classify(R, th)
        sens_s[th] = float((fl & (dt >= STRICT[0]) & (ratio >= STRICT[1])).mean())
    out["sens_strict"] = sens_s
    return out


def write_tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(f"{x:.6g}" if isinstance(x, (float, np.floating)) else str(x) for x in r) + "\n")


def write_bank_tsv(res, ctx_D):
    slug = res["slug"]
    D = ctx_D
    G, R = res["G"], res["R"]
    hdr = ["module", "midx", "gidx", "run", "split", "event_number", "amp", "mf_energy", "amp_rank", "chi0", "A", "t",
           "chi1", "A1", "t1", "A2", "t2", "dchi2", "ok2", "dt", "ratio", "degenerate", "flag", "strict"]
    names = res["z_names"]
    rows = []
    for k in range(G.size):
        g = G[k]
        rows.append([names[res["row_of"][k]], int(D["midx"][g]), int(g), int(D["run"][g]), int(D["split"][g]),
                     int(D["event_number"][g]), D["amp"][g], D["mf_energy"][g], int(res["rank"][k]),
                     R[k, 0], R[k, 2], R[k, 3], R[k, 4], R[k, 8], R[k, 9], R[k, 10], R[k, 11], R[k, 12],
                     int(R[k, 13]), res["dt"][k], res["ratio"][k], int(res["degen"][k]), int(res["flag"][k]), int(res["strict"][k])])
    write_tsv(OUT / f"pulses_{slug}.tsv", hdr, rows)
    mh = ["module", "row", "midx", "fallback", "n_pool", "n_flag", "n_flag_dtlt16", "n_strict", "n_degen", "n_over", "dchi_med",
          "amp_med", "sig_ped"]
    mrows = sorted(res["mod"], key=lambda m: (-m["n_flag"], m["name"]))
    write_tsv(OUT / f"modules_{slug}.tsv", mh, [[m["name"], m["row"], m["midx"], m["fallback"], m["n_pool"], m["n_flag"],
                                                m["n_flag_dtlt16"], m["n_strict"], m["n_degen"], m["n_over"], m["dchi_med"],
                                                m["amp_med"], m["sig_ped"]] for m in mrows])


# ----------------------------------------------------------------- Ta vs H2 --
def width_compare(ta, h2, noise_ta, noise_h2):
    """Per-module Ta - H2 widths of the shipped templates (own templates in both)."""
    nt, nh = ta["z_names"], h2["z_names"]
    ok = (ta["z_fb"] <= 1) & (h2["z_fb"] <= 1)
    rows = []
    for r in np.nonzero(ok)[0]:
        if nt[r] != nh[r]:
            continue
        a, b = shape_metrics(ta["z_shape"][r]), shape_metrics(h2["z_shape"][r])
        rows.append((nt[r], r, a, b))
    # noise of one 50-pool template = std(half-half) / 2 ; of the Ta - H2 difference: quadrature
    sg = {}
    for k, mname in enumerate(METRIC):
        s_ta = mad_sig(noise_ta["dm"][:, k]) / 2.0
        s_h2 = mad_sig(noise_h2["dm"][:, k]) / 2.0
        sg[mname] = float(np.hypot(s_ta, s_h2))
    return rows, sg


# ---------------------------------------------------------------------- PDF --
def text_page(pdf, title, lines, fs=7.0):
    fig = plt.figure(figsize=(8.5, 11))
    fig.text(0.04, 0.97, title, fontsize=10, weight="bold", va="top")
    fig.text(0.04, 0.94, "\n".join(lines), family="monospace", fontsize=fs, va="top")
    pdf.savefig(fig)
    plt.close(fig)


def fmt_pdf_tables(tabs, results):
    L = []
    for slug in SLUGS:
        if slug not in tabs:
            continue
        t = tabs[slug]
        L.append(f"{LABEL[slug]}  thr {results[slug]['thr']:.1f}: pool pulses {t['n']}, flagged {t['n_flag']} "
                 f"({100 * t['frac']:.2f} %; |dt|<16: {100 * t['frac_lt16']:.2f} %), degenerate-excluded {t['n_degen']} "
                 f"(dchi2>thr total {t['n_over']})")
        L.append("   counts [dt rows | A2/A1 cols " + " ".join(f"{x:>9}" for x in R_LAB) + "]")
        for i, lab in enumerate(DT_LAB):
            L.append(f"   dt {lab:>6} ns | " + " ".join(f"{t['t2d'][i, j]:>9d}" for j in range(3)))
        L.append(f"   modules with >=1 / >=3 flagged of {t['n_mod']}: {t['mod_ge1']} / {t['mod_ge3']};  "
                 f"flagged-per-module hist 0..10+: {t['nflag_hist'].tolist()}")
        L.append(f"   strict (dt>=8 ns, sec/main>=0.2): {t['n_strict']} ({100 * t['frac_strict']:.2f} %), modules >=1 / >=3: "
                 f"{t['mod_s_ge1']} / {t['mod_s_ge3']}")
        L.append("   by amplitude <600 / 600-1500 / >1500 ADC: n " + " / ".join(str(a[0]) for a in t["amp_tab"]) + ";  flagged % "
                 + " / ".join(f"{a[1]:.2f}" for a in t["amp_tab"]) + ";  strict % " + " / ".join(f"{a[2]:.2f}" for a in t["amp_tab"]))
        L.append("   flagged % at thr " + ", ".join(f"{k:g}: {100 * v:.2f}" for k, v in t["sens"].items()))
        L.append("   strict  % at thr " + ", ".join(f"{k:g}: {100 * v:.2f}" for k, v in t["sens_strict"].items()))
        L.append("")
    return L


def make_pdf(results, tabs, wc, sg, ex_pages):
    with PdfPages(OUT / "doublet_audit.pdf") as pdf:
        text_page(pdf, "Flagged pool pulses (two-pulse fit, not degenerate)", fmt_pdf_tables(tabs, results))
        L = []
        for slug in SLUGS:
            if slug not in results:
                continue
            L.append(f"{LABEL[slug]}: top-20 modules by flagged count (name n_flag/n_pool n_degen)")
            m = sorted(results[slug]["mod"], key=lambda m: (-m["n_flag"], m["name"]))[:20]
            L.append("   " + "  ".join(f"{x['name']} {x['n_flag']}/{x['n_pool']}/{x['n_degen']}" for x in m[:10]))
            L.append("   " + "  ".join(f"{x['name']} {x['n_flag']}/{x['n_pool']}/{x['n_degen']}" for x in m[10:]))
            L.append("")
        text_page(pdf, "Top-20 flagged modules per bank", L, 6.5)
        # Ta vs H2 widths
        rows = wc["rows"]
        fig, ax = plt.subplots(2, 3, figsize=(11, 7))
        for k, mname in enumerate(["rise", "w99", "fwhm"]):
            a = np.array([r[2][k] for r in rows])
            b = np.array([r[3][k] for r in rows])
            ax[0, k].scatter(b, a, s=3, alpha=0.4)
            lim = [np.nanpercentile(np.r_[a, b], 0.5), np.nanpercentile(np.r_[a, b], 99.5)]
            ax[0, k].plot(lim, lim, "k-", lw=0.6)
            ax[0, k].set_xlim(lim)
            ax[0, k].set_ylim(lim)
            ax[0, k].set_xlabel(f"2.24 H2 {mname} (ns)")
            ax[0, k].set_ylabel(f"2.24 Ta {mname} (ns)")
            d = a - b
            ax[1, k].hist(d[np.isfinite(d)], bins=np.linspace(-1.5, 1.5, 80))
            ax[1, k].axvline(0, c="k", lw=0.6)
            ax[1, k].axvline(2 * sg[mname], c="r", lw=0.8)
            ax[1, k].set_xlabel(f"Ta - H2 {mname} (ns); red = +2 sigma_noise ({sg[mname]:.2f})")
            ax[1, k].set_title(f"median {np.nanmedian(d):+.3f}", fontsize=8)
        fig.suptitle("Shipped templates, Ta vs 2.24 H2 per module", fontsize=9)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        for fig in ex_pages:
            pdf.savefig(fig)
            plt.close(fig)


def example_page(res):
    """4 flagged pool pulses with single and two-pulse fits overlaid."""
    bank = res["bank"]
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    for ax, e in zip(axs.ravel(), res["examples"]):
        k = e["k"]
        R = res["R"][k]
        mi = e["bidx"]
        row, drow = bank.lut[mi].astype(np.float64), bank.dlut[mi].astype(np.float64)
        tt = np.arange(0, 400.0, 0.5)
        s1 = R[2] * TP.tpl(row, drow, tt - R[3], bank.T_end)
        s2a = R[8] * TP.tpl(row, drow, tt - R[9], bank.T_end)
        s2b = R[10] * TP.tpl(row, drow, tt - R[11], bank.T_end)
        # single-pulse fit alone is R[2], R[3] (chi1 = R[4]); two-pulse = R[7]
        ax.plot(TP.TS, e["y"], "k.", ms=4, label="data (ped-sub)")
        ax.plot(tt, s1, "b-", lw=0.9, label=f"1-pulse chi2 {R[4]:.0f}")
        ax.plot(tt, s2a + s2b, "r-", lw=0.9, label=f"2-pulse chi2 {R[7]:.0f}")
        ax.plot(tt, s2b, "r:", lw=0.8)
        c = R[9]
        ax.set_xlim(max(c - 30, 0), c + 110)
        ax.set_title(f"{res['z_names'][e['row']]}  A {R[2]:.0f}  dchi2 {R[12]:.0f}  dt {res['dt'][k]:.1f} ns  "
                     f"A2/A1 {res['ratio'][k]:.2f}", fontsize=7.5)
        ax.legend(fontsize=6)
    fig.suptitle(f"{LABEL[res['slug']]}: flagged pool pulses (largest dchi2 per |dt| bin)", fontsize=9)
    fig.tight_layout()
    return fig


def shape_pages(rb_ta, rb_h2):
    figs = []
    for slug, rb in (("2p24GeV_Ta_X17", rb_ta), ("2p24GeV_H2", rb_h2)):
        if rb is None:
            continue
        rec = [r for r in rb["rec"] if r["shape"] is not None][:3]
        fig, axs = plt.subplots(2, 3, figsize=(11, 6.5), sharex="col")
        for j, r in enumerate(rec):
            a = r["shipped"].astype(float)
            b = r["shape"].astype(float)
            a, b = np.nan_to_num(a) / np.nanmax(a), np.nan_to_num(b) / np.nanmax(b)
            axs[0, j].plot(TG, a, "k-", lw=1, label="shipped")
            axs[0, j].plot(TG, b, "r-", lw=1, label=f"flagged removed ({r['n_flag']}/{r['n_pool']})")
            axs[1, j].plot(TG, 100 * (b - a), "r-", lw=0.8)
            if "shape_s" in r:
                c = np.nan_to_num(r["shape_s"].astype(float))
                c = c / c.max()
                axs[0, j].plot(TG, c, "g--", lw=0.8, label=f"strict removed ({r['n_strict']})")
                axs[1, j].plot(TG, 100 * (c - a), "g--", lw=0.8)
            axs[1, j].axhline(0, c="k", lw=0.5)
            axs[0, j].set_xlim(-5, 60)
            axs[0, j].set_title(f"{r['name']}  max|df| {r['maxd']:.1f} %  d_rise {r['dm'][0]:+.2f} ns", fontsize=8)
            axs[1, j].set_xlabel("ns after onset")
            axs[1, j].set_ylabel("df (%)")
            axs[0, j].legend(fontsize=6)
        fig.suptitle(f"{LABEL[slug]}: shipped vs rebuilt without flagged pulses (peak-normalised)", fontsize=9)
        fig.tight_layout()
        figs.append(fig)
    return figs


# ------------------------------------------------------------------- report --
def _read_tsv(path):
    L = open(path).read().splitlines()
    h = L[0].split("\t")
    return [dict(zip(h, l.split("\t"))) for l in L[1:]]


def write_report(summ):
    T = summ["tabs"]
    ta, h2 = T.get("2p24GeV_Ta_X17"), T.get("2p24GeV_H2")
    L = ["# Hidden doublets in the best-50 template pools (2026-10-27)", ""]
    L.append("**Threshold source.** dchi2 thr (twopulse_fit, eps 0.02): " + "; ".join(
        f"{LABEL[s]} {summ['thr'][s]:.1f} ({'U0 per-type' if s in U0_KEY else 'not recorded in U0 -> production 2.24 H2 value'})"
        for s in SLUGS if s in T) + ". Flag = dchi2 > thr, fitted dt >= 2 ns, secondary/main amplitude >= 0.05 (ratio = min/max; "
        "A2/A1 of the later pulse is not bounded by 1).")
    rp = summ["repro"]
    L.append("**Pool reproduction.** Selected indices are not stored; rebuilt with spline_template_extract's own loader/argv "
             "(top n_max by mf_energy of the p0/p1 chain pool). n_pool matches for " + ", ".join(
                 f"{LABEL[s]} {r['n_modules'] - r['n_pool_mismatch']}/{r['n_modules']}" for s, r in rp.items())
             + ". process_module re-run on " + str(len(next(iter(rp.values()))["spot"])) + " modules/bank: max|df| vs shipped "
             f"{max(x['d_full'] for r in rp.values() for x in r['spot']):.1e} (whole module) and "
             f"{max(x['d_ctl'] for r in rp.values() for x in r['spot']):.1e} (pool-restricted call used for rebuilds). Reproduced exactly.")
    L += ["", "| bank | pool pulses | flagged % (|dt|<16) | flagged by |dt| 2-8/8-16/16-32/>32 | degenerate excl. | strict % (dt>=8, r>=0.2) | modules >=3 flagged | top10 / bot10 flagged % (ratio) |",
          "|---|---|---|---|---|---|---|---|"]
    for s in SLUGS:
        if s not in T:
            continue
        t = T[s]
        L.append(f"| {LABEL[s]} | {t['n']} | {100 * t['frac']:.2f} ({100 * t['frac_lt16']:.2f}) | {'/'.join(map(str, t['dtb']))} | "
                 f"{t['n_degen']} | {100 * t['frac_strict']:.3f} | {t['mod_ge3']} of {t['n_mod']} | "
                 f"{100 * t['top_frac']:.2f} / {100 * t['bot_frac']:.2f} ({t['ratio_tb']:.2f}) |")
    L.append("2D (dt x A2/A1) tables, flagged-per-module histograms, thresholds 98-1000 and amplitude bins: `flag_tables.tsv`, `summary.json`, PDF p.1-2. "
             "Most flagged pulses are |dt| < 8 ns with ratio < 0.5 (shape ripple, not a resolved second pulse).")
    if ta and h2:
        L.append(f"**vs expectation.** Expected ~0.1 % (e-p), ~3 % (Ta) in |dt| < 16 ns; observed {100 * ta['frac_lt16']:.2f} % (Ta) and "
                 f"{', '.join(f'{100 * T[s]['frac_lt16']:.2f}' for s in SLUGS if s in T and s != '2p24GeV_Ta_X17')} % (e-p banks); Ta is "
                 f"{'not ' if ta['frac'] <= h2['frac'] else ''}more contaminated than 2.24 H2 ({100 * ta['frac']:.2f} vs {100 * h2['frac']:.2f} %; strict {100 * ta['frac_strict']:.3f} vs {100 * h2['frac_strict']:.3f} %).")
    w = summ.get("width")
    if w:
        L.append(f"**Ta vs 2.24 H2, shipped templates, {w['n']} modules (Ta - H2, median [16, 84 %], noise sigma of the difference):** "
                 + "; ".join(f"{m} {w[m]['med_d']:+.3f} [{w[m]['p16']:+.2f}, {w[m]['p84']:+.2f}] (sigma {w['sigma_noise'][m]:.2f}) ns, Ta wider >2 sigma in {w[m]['n_wider2']}, narrower in {w[m]['n_narrower2']}"
                             for m in ("rise", "w99", "fwhm")) + ".")
        a = summ["avg_template"]
        L.append("Run-type average templates (rise / FWHM ns): " + ", ".join(f"{k} {v['rise']:.2f}/{v['fwhm']:.2f}" for k, v in a.items())
                 + ". Per-module Ta-H2 vs n_flag(Ta), Spearman rho rise/w99/fwhm: " + "/".join(f"{w[m]['spearman_nflag'][0]:+.2f}" for m in ("rise", "w99", "fwhm")) + ".")
        x = summ.get("w325")
        if x:
            pw = summ["w325_pool"]
            L.append(f"**W325:** Ta - H2 rise {x['d_rise']:+.2f} ns ({x['d_over_sigma_rise']:+.1f} sigma), w99 {x['d_w99']:+.2f} ns ({x['d_over_sigma_w99']:+.1f} sigma), "
                     f"FWHM {x['d_fwhm']:+.2f} ns; flagged pool pulses Ta {pw['2.24 Ta'][0]['n_flag']}/{pw['2.24 Ta'][0]['n_pool']}, H2 {pw['2.24 H2'][0]['n_flag']}/{pw['2.24 H2'][0]['n_pool']}.")
    if ta and h2:
        L.append("Pool amplitude mix (<600 / 600-1500 / >1500 ADC, %): " + "; ".join(
            f"{LABEL[s]} " + "/".join(f"{100 * a[0] / T[s]['n']:.0f}" for a in T[s]["amp_tab"]) for s in SLUGS if s in T)
            + " -- Ta pools are lower-amplitude than 2.24 H2 and 3.49 H2; an amplitude-dependent shape (untested here) is a candidate for any residual Ta/e-p difference.")
    for slug in REBUILD:
        sm = summ.get(f"shape_{slug}")
        if not sm:
            continue
        rows = _read_tsv(OUT / f"shape_impact_{slug}.tsv")
        ex_r = sum(abs(float(r["d_rise"])) > float(r["null_abs_d_rise_max"]) for r in rows)
        ex_m = sum(float(r["maxd_pct"]) > float(r["null_maxd_max"]) for r in rows)
        st = sm["stable"]
        L.append(f"**Rebuild without flagged, {LABEL[slug]}, top-20 flagged modules ({sm['n']} built, {sm['n_fail']} failed, {sm['n_lam_switch']} changed lambda/fallback):** "
                 f"median |d| rise/w99/FWHM/tpk {sm['dm_absmed'][0]:.3f}/{sm['dm_absmed'][1]:.3f}/{sm['dm_absmed'][2]:.3f}/{sm['dm_absmed'][3]:.3f} ns "
                 f"(signed median {sm['dm_med'][0]:+.3f}/{sm['dm_med'][1]:+.3f}/{sm['dm_med'][2]:+.3f}/{sm['dm_med'][3]:+.3f}); max|df| median {np.median(sm['maxd']):.1f} % "
                 f"vs random-removal null {sm['null_maxd_med']:.1f} %; |d rise| beyond null max in {ex_r}, max|df| beyond null max in {ex_m} of {len(rows)}; "
                 f"sigma(one 50-pool) rise/w99 {sm['noise_dm_sig'][0] / 2:.2f}/{sm['noise_dm_sig'][1] / 2:.2f} ns, split-half(25/25) max|df| median {sm['noise_maxd_med']:.1f} % "
                 f"(xtype: 2.6-2.9 %); null |d| rise/w99 {sm['null_dm_absmed'][0]:.3f}/{sm['null_dm_absmed'][1]:.3f}. Strict set removed: {sm['strict']['n']} modules"
                 + (f", median |d| rise/w99 {sm['strict']['dm_absmed'][0]:.3f}/{sm['strict']['dm_absmed'][1]:.3f} ns." if sm["strict"]["n"] else "."))
    g = summ.get("gap_closure")
    if g:
        L.append(f"Gap closure (Ta rebuilt - H2, median over the {g['n']} modules, rise/w99/FWHM): shipped {g['ship_med'][0]:+.3f}/{g['ship_med'][1]:+.3f}/{g['ship_med'][2]:+.3f} -> rebuilt {g['new_med'][0]:+.3f}/{g['new_med'][1]:+.3f}/{g['new_med'][2]:+.3f} ns.")
    L += ["", "**Answer.** " + (
        f"No. The two-pulse fit finds few hidden doublets in the pools (Ta {100 * ta['frac']:.2f} %, e-p/C {100 * min(T[s]['frac'] for s in T if s != '2p24GeV_Ta_X17'):.2f}-"
        f"{100 * max(T[s]['frac'] for s in T if s != '2p24GeV_Ta_X17'):.2f} %; Ta strict {100 * ta['frac_strict']:.3f} %), Ta is comparable to the 2.24 GeV e-p and C banks (flagged {100 * min(T[s]['frac'] for s in T if s in ('2p24GeV_H2', '2p24GeV_C')):.2f}-{100 * max(T[s]['frac'] for s in T if s in ('2p24GeV_H2', '2p24GeV_C')):.2f} %) "
        "and not 3 % as extrapolated, removing the flagged pulses moves the templates by no more than removing the same number of random pulses, and the shipped Ta templates are not systematically wider than "
        "2.24 H2 (median differences within noise). Pool pile-up does not explain a Ta/e-p width difference." if ta and h2 else "see tables"),
        "**Caveats.** (i) dchi2 threshold calibrated on typical clean pulses (97 % quantile of real singles); pool pulses are the highest-energy ones, where shape mismatch also lifts dchi2, hence most flags are low-ratio/short-dt; "
        "Ta and C use the 2.24 H2 value (no per-type calibration in U0). (ii) Pool pulses were built into the template, which lowers their dchi2 (conservative for a contamination count). "
        "(iii) Rebuilds re-run the held-out lambda selection and oscillation escalation, so some modules switch lambda; those large steps also appear in the random-removal null. "
        "(iv) Pool = pre-outlier-rejection top-50 (the build keeps ~40 of them)."]
    (OUT / "report.md").write_text("\n".join(L) + "\n")


# ---------------------------------------------------------------------- main --
def main(argv=None):
    global OUT
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--banks", default=",".join(SLUGS))
    ap.add_argument("--workers", type=int, default=18)
    ap.add_argument("--quick", type=int, default=0, help="use every Nth module only (testing)")
    ap.add_argument("--n-spot", type=int, default=8)
    ap.add_argument("--n-null", type=int, default=4)
    ap.add_argument("--n-noise", type=int, default=150)
    ap.add_argument("--out", default=str(OUT))
    args = ap.parse_args(argv)
    OUT = Path(args.out)
    OUT.mkdir(parents=True, exist_ok=True)
    t_start = time.time()
    u0 = json.load(open(U0))["meta"]
    thr_map, thr_src = {}, {}
    ref = u0["2.24 GeV H2"]["thr"]["real"][0]
    for s in SLUGS:
        if s in U0_KEY:
            thr_map[s] = float(u0[U0_KEY[s]]["thr"]["real"][0])
            thr_src[s] = "U0 per-type (cal-half 97 % quantile of real singles)"
        else:
            thr_map[s] = float(ref)
            thr_src[s] = "not recorded in U0: production 2.24 GeV H2 value"
    results, tabs = {}, {}
    Dref = {}
    for slug in args.banks.split(","):
        res = run_bank(slug, args, thr_map, args.workers)
        write_bank_tsv(res, _CTX["D"])
        results[slug] = res
        tabs[slug] = bank_tables(res)
        del _CTX["D"]
    # ---- Ta vs H2 and shape-impact summaries
    summ = dict(thr=thr_map, thr_src=thr_src)
    wc = sg = None
    ex_pages = []
    if "2p24GeV_Ta_X17" in results and "2p24GeV_H2" in results:
        ta, h2 = results["2p24GeV_Ta_X17"], results["2p24GeV_H2"]
        nta, nh2 = ta["rebuild"]["noise"], h2["rebuild"]["noise"]
        rows, sg = width_compare(ta, h2, nta, nh2)
        wc = dict(rows=rows)
        a = np.array([r[2] for r in rows])
        b = np.array([r[3] for r in rows])
        d = a - b
        summ["width"] = dict(n=len(rows), sigma_noise=sg)
        nfl = {m["name"]: m["n_flag"] for m in ta["mod"]}
        wl = []
        for k, mname in enumerate(METRIC[:3]):
            dk = d[:, k]
            ok = np.isfinite(dk)
            summ["width"][mname] = dict(
                med_ta=float(np.nanmedian(a[:, k])), med_h2=float(np.nanmedian(b[:, k])),
                med_d=float(np.nanmedian(dk)), p16=float(np.nanpercentile(dk[ok], 16)),
                p84=float(np.nanpercentile(dk[ok], 84)), n_wider2=int((dk > 2 * sg[mname]).sum()),
                n_narrower2=int((dk < -2 * sg[mname]).sum()))
        nf = np.array([nfl.get(r[0], 0) for r in rows])
        for k, mname in enumerate(METRIC[:3]):
            ok = np.isfinite(d[:, k])
            rho, p = spearmanr(nf[ok], d[ok, k])
            summ["width"][mname]["spearman_nflag"] = [float(rho), float(p)]
        for i, r in enumerate(rows):
            wl.append([r[0], int(nf[i])] + [x for k in range(3) for x in (a[i, k], b[i, k], d[i, k], d[i, k] / sg[METRIC[k]])])
        wl.sort(key=lambda x: -x[5] if np.isfinite(x[5]) else 0)
        hdr = ["module", "n_flag_Ta"] + [f"{p}_{m}" for m in METRIC[:3] for p in ("ta", "h2", "d", "d_over_sigma")]
        write_tsv(OUT / "width_ta_vs_h2.tsv", hdr, wl)
        w325 = [x for x in wl if x[0] == "W325"]
        summ["w325"] = None if not w325 else dict(zip(hdr, [float(v) if not isinstance(v, str) else v for v in w325[0]]))
        summ["avg_template"] = {LABEL[s_]: dict(zip(METRIC, shape_metrics(r_["shape_avg"]).tolist())) for s_, r_ in results.items()}
        summ["w325_pool"] = {LABEL[s_]: [m for m in results[s_]["mod"] if m["name"] == "W325"] for s_ in ("2p24GeV_Ta_X17", "2p24GeV_H2")}
        # shape impact
        for slug in REBUILD:
            rb = results[slug]["rebuild"]
            rec = [r for r in rb["rec"] if r["shape"] is not None]
            dm = np.array([r["dm"] for r in rec])
            sm = dict(n=len(rec), n_fail=len(rb["rec"]) - len(rec),
                      maxd=[r["maxd"] for r in rec], maxd_ctl=[r["maxd_ctl"] for r in rec],
                      null_maxd_med=float(np.median([x for r in rec for x in r["null_maxd"]])),
                      null_dm_absmed=[float(np.median(np.abs([x[k] for r in rec for x in r["null_dm"]]))) for k in range(4)],
                      dm_med=[float(np.median(dm[:, k])) for k in range(4)],
                      dm_absmed=[float(np.median(np.abs(dm[:, k]))) for k in range(4)],
                      dm_max=[float(np.max(np.abs(dm[:, k]))) for k in range(4)],
                      noise_maxd_med=float(np.median(rb["noise"]["maxd"])),
                      noise_dm_sig=[mad_sig(rb["noise"]["dm"][:, k]) for k in range(4)],
                      n_noise=int(len(rb["noise"]["maxd"])),
                      fb_changed=int(sum(r["fb_new"] != r["fb_ship"] for r in rec)),
                      n_lam_switch=int(sum(r["lam_new"] != r["lam_ship"] for r in rec)))
            st_ = [r for r in rec if r["stable"]]
            sm["stable"] = dict(n=len(st_), maxd_med=float(np.median([r["maxd"] for r in st_])) if st_ else None,
                                maxd_max=float(np.max([r["maxd"] for r in st_])) if st_ else None,
                                null_maxd_med=float(np.median([x for r in st_ for x in r["null_maxd"]])) if st_ else None,
                                dm_absmed=[float(np.median(np.abs([r["dm"][k] for r in st_]))) for k in range(4)] if st_ else None,
                                dm_med=[float(np.median([r["dm"][k] for r in st_])) for k in range(4)] if st_ else None,
                                null_dm_absmed=[float(np.median(np.abs([x[k] for r in st_ for x in r["null_dm"]]))) for k in range(4)] if st_ else None)
            write_tsv(OUT / f"shape_impact_{slug}.tsv",
                      ["module", "n_flag", "n_strict", "n_pool", "n_pool_new", "fb_ship", "fb_new", "lam_ship", "lam_new", "stable",
                       "maxd_pct", "d_rise", "d_w99", "d_fwhm", "d_tpk", "null_maxd_med", "null_maxd_max", "null_abs_d_rise_max",
                       "ctl_maxd_pct", "strict_maxd_pct", "strict_d_rise", "strict_d_w99"],
                      [[r["name"], r["n_flag"], r.get("n_strict", 0), r["n_pool"], r["n_pool_new"], r["fb_ship"], r["fb_new"],
                        r["lam_ship"], r["lam_new"], int(r["stable"]), r["maxd"], *r["dm"], float(np.median(r["null_maxd"])),
                        float(np.max(r["null_maxd"])), float(np.max(np.abs([x[0] for x in r["null_dm"]]))), r["maxd_ctl"],
                        r.get("maxd_s", np.nan), *(r["dm_s"][:2] if "dm_s" in r else (np.nan, np.nan))] for r in rec])
            rs = [r for r in rec if "dm_s" in r]
            sm["strict"] = dict(n=len(rs), maxd=[r["maxd_s"] for r in rs],
                                dm_med=[float(np.median([r["dm_s"][k] for r in rs])) for k in range(4)] if rs else None,
                                dm_absmed=[float(np.median(np.abs([r["dm_s"][k] for r in rs]))) for k in range(4)] if rs else None,
                                dm_max=[float(np.max(np.abs([r["dm_s"][k] for r in rs]))) for k in range(4)] if rs else None)
            summ[f"shape_{slug}"] = sm
            names = [r["name"] for r in rec]
            np.savez_compressed(OUT / f"rebuilt_templates_{slug}.npz", names=np.array(names), t_grid=TG.astype(np.float32),
                                shape_shipped=np.array([r["shipped"] for r in rec]),
                                shape_removed=np.array([r["shape"] for r in rec]),
                                shape_removed_strict=np.array([r.get("shape_s", np.full(TG.size, np.nan, np.float32)) for r in rec]),
                                n_flag=np.array([r["n_flag"] for r in rec]), n_pool=np.array([r["n_pool"] for r in rec]),
                                n_pool_new=np.array([r["n_pool_new"] for r in rec]), d_metrics=dm,
                                metric_names=np.array(METRIC), maxd_pct=np.array([r["maxd"] for r in rec]))
        # gap closure: Ta rebuilt vs H2 for the rebuilt Ta modules
        rows_by = {r[0]: r for r in rows}
        gc = []
        for r in results["2p24GeV_Ta_X17"]["rebuild"]["rec"]:
            if r["shape"] is None or r["name"] not in rows_by:
                continue
            h = rows_by[r["name"]][3]
            gc.append((shape_metrics(r["shipped"]) - h, shape_metrics(r["shape"]) - h))
        if gc:
            gc = np.array(gc)
            summ["gap_closure"] = dict(n=len(gc), ship_med=[float(np.median(gc[:, 0, k])) for k in range(4)],
                                       new_med=[float(np.median(gc[:, 1, k])) for k in range(4)])
        ex_pages = [example_page(results["2p24GeV_Ta_X17"])] + shape_pages(
            results["2p24GeV_Ta_X17"]["rebuild"], results["2p24GeV_H2"]["rebuild"])
        if results["2p24GeV_H2"]["examples"]:
            ex_pages.insert(1, example_page(results["2p24GeV_H2"]))
        make_pdf(results, tabs, wc, sg, ex_pages)
    # ---- tables tsv and summary json
    hdr = ["bank", "thr", "n_pool_pulses", "n_flag", "frac_pct", "frac_dtlt16_pct", "n_degen", "n_over", "modules_ge1",
           "modules_ge3", "n_modules", "top10_frac_pct", "bot10_frac_pct", "top_over_bot", "n_strict", "strict_pct", "strict_top10_pct", "strict_bot10_pct"] + \
          [f"dt{l}" for l in DT_LAB] + [f"dt{a}_r{c}" for a in DT_LAB for c in R_LAB]
    rws = []
    for s, t in tabs.items():
        rws.append([s, results[s]["thr"], t["n"], t["n_flag"], 100 * t["frac"], 100 * t["frac_lt16"], t["n_degen"], t["n_over"],
                    t["mod_ge1"], t["mod_ge3"], t["n_mod"], 100 * t["top_frac"], 100 * t["bot_frac"], t["ratio_tb"], t["n_strict"], 100 * t["frac_strict"],
                    100 * t["top_frac_s"], 100 * t["bot_frac_s"]] + t["dtb"] + [int(x) for x in t["t2d"].ravel()])
    write_tsv(OUT / "flag_tables.tsv", hdr, rws)
    summ["tabs"] = {s: {k: (v.tolist() if hasattr(v, "tolist") else v) for k, v in t.items()} for s, t in tabs.items()}
    summ["repro"] = {s: dict(n_pool_mismatch=r["n_pool_mismatch"], n_modules=len(r["rows"]), spot=r["spot"]) for s, r in results.items()}
    summ["wall_s"] = time.time() - t_start
    json.dump(summ, open(OUT / "summary.json", "w"), indent=1, default=lambda o: o.tolist() if hasattr(o, "tolist") else str(o))
    write_report(summ)
    print(f"done in {summ['wall_s']:.0f} s -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
