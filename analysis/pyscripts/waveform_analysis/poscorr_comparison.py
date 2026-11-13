#!/usr/bin/env python3
"""poscorr_comparison.py -- scratch driver: position-correction variants for the PbWO4 energy-resolution benchmark.

    python3 waveform_analysis/poscorr_comparison.py <cache dir>

Imports energy_resolution_benchmark as a module (no scan; reads output/energy_resolution/scan/*/*.npz).  Setting fixed for every
variant: gains on, Moller rejection off, non-linearity off (= the benchmark's arm gp).  Writes energy_resolution/poscorr_comparison.tsv
and .pdf.  Step 1 clusters every split once (legacy and gate 12 ns) into the cache; step 2 gains per (run, training set); step 3
position-correction variants per (run type, training set); step 4 evaluation per split.
"""
import math
import multiprocessing as mp
import pickle
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import energy_resolution_benchmark as bm   # noqa: E402

NW = 6
GATE = 12.0
SCAN = HERE.parent / "output" / "energy_resolution" / "scan"
OUT = HERE / "energy_resolution"
RTS = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old", "2p24GeV_Ta_X17"]
BORROW = {"2p24GeV_Ta_X17": "2p24GeV_H2"}
CURRENT = {"0p73GeV_H2": "map12", "2p24GeV_H2": "hybrid3000", "3p49GeV_H2": "hybrid3000", "3p49GeV_H2_old": "hybrid3000",
           "2p24GeV_Ta_X17": "map12"}
VARIANTS = ["none", "map8", "map12", "map16", "map24", "hybrid1000", "hybrid3000", "hybrid10000", "even4_global", "even6_global",
            "even4_pm", "even6_pm"]
BINS = {"map8": 8, "map12": 12, "map16": 16, "map24": 24}
OWN_N = {"hybrid1000": 1000, "hybrid3000": 3000, "hybrid10000": 10000}
TERMS = {"even4": 3, "even6": 6}          # leading columns of (1, u^2, v^2, u^2 v^2, u^4, v^4, u, v) kept
PATHS = ("legacy", f"g{GATE:g}")
CACHE = Path(sys.argv[1] if len(sys.argv) > 1 else "/tmp/sup-9d445b-pc.3e5S/cache")


# --- step 1: clustering cache --------------------------------------------------------------------------------------------
def cache_job(a):
    fn, split = a
    out = CACHE / f"{Path(fn).stem}_{rt_of(fn)}_{split}.pkl"
    if out.exists():
        return str(out)
    z = np.load(fn, allow_pickle=False)
    run, rt = int(z["run"]), str(z["run_type"])
    cls, geo = bm.clusterers(run, str(z["calib"]), (GATE,))
    kE, kS = round(geo["ebeam"] / bm.E_REF, 6), round(math.sqrt(geo["ebeam"] / bm.E_REF), 6)
    evs, h_ev = np.flatnonzero(z["ev_split"] == split), z["h_ev"]
    lo, hi = np.searchsorted(h_ev, [evs[0], evs[-1] + 1]) if evs.size else (0, 0)
    h = {q: z[q][lo:hi] for q in ("h_mod", "h_t", "h_e", "h_k", "h_extra")} | {"h_ev": h_ev[lo:hi]}
    t0 = time.monotonic()
    base, _ = bm.cluster_events(h, cls, geo, bm.SIGNAL_MIN * kE, bm.SEED_MIN * kE)
    with open(out, "wb") as f:
        pickle.dump(dict(run=run, rt=rt, split=split, geo=geo, kE=kE, kS=kS, pools={p: base[p] for p in PATHS}), f)
    print(f"  [cache] {out.name}: {time.monotonic() - t0:.0f} s", flush=True)
    return str(out)


def rt_of(fn):
    return Path(fn).parent.name


def load(p):
    with open(p, "rb") as f:
        return pickle.load(f)


# --- step 2: gains ---------------------------------------------------------------------------------------------------------
def gain_job(a):
    files, rt = a
    ds = [load(f) for f in files]
    geo, kE, kS = ds[0]["geo"], ds[0]["kE"], ds[0]["kS"]
    pool = bm.nl_view(bm.cat([d["pools"]["legacy"] for d in ds]), 0)
    g = bm.solve_gains(pool, geo, bm.pool_mass(rt), kE, kS)
    return g["g"] if g is not None else None


def ecal_of(pool, g):
    counts = np.array([m.size for m in pool["mi"]], np.int64)
    return np.bincount(np.repeat(np.arange(counts.size), counts), weights=np.concatenate(pool["me"]) * g[np.concatenate(pool["mi"]).astype(np.int64)],
                       minlength=counts.size)


# --- step 3: position-correction fits ----------------------------------------------------------------------------------------
def fit_poly(u, v, r, k):
    """bm.fit_module with the first k of the eight terms (k = 8: identical)."""
    cen = -bm.PC_CLAMP + (np.arange(bm.PC_GRID) + 0.5) * 2 * bm.PC_CLAMP / bm.PC_GRID
    iu, iv, inb = bm.bin_index(u, v, bm.PC_GRID, bm.PC_CLAMP)
    pts = []
    for b, idx in bm.groups(np.where(inb, iu * bm.PC_GRID + iv, -1)).items():
        t = bm.trunc_mean(r[idx]) if b >= 0 and idx.size >= bm.PC_GRID_MIN else None
        if t:
            pts.append((cen[b // bm.PC_GRID], cen[b % bm.PC_GRID], t[0], t[1] / math.sqrt(t[2])))
    if len(pts) < bm.PC_GRID_BINS:
        return None
    uc, vc, Rv, se = np.array(pts).T
    X = np.stack([np.ones_like(uc), uc ** 2, vc ** 2, uc ** 2 * vc ** 2, uc ** 4, vc ** 4, uc, vc], axis=-1)[:, :k]
    a = np.linalg.lstsq(X / se[:, None], Rv / se, rcond=None)[0]
    c = np.zeros(8)
    c[0], c[1:k] = a[0], a[1:] / a[0]
    c[0] /= bm.poly_r(c, u, v).mean()
    g = np.linspace(-bm.PC_CLAMP, bm.PC_CLAMP, 13)
    Rg = bm.poly_r(c, *np.meshgrid(g, g))
    ok = np.isfinite(c).all() and bm.PC_ACCEPT[0] <= Rg.min() and Rg.max() <= bm.PC_ACCEPT[1]
    return c if ok else None


def pc_job(a):
    """All variants of one (run type, training set): {variant: (spec, diagnostics)}; spec = dict(R, own, glob) or None (R = 1)."""
    rt, files_by_run, g_by_run = a
    u, v, s, e = [], [], [], []
    for run in sorted(files_by_run):
        ds = [load(f) for f in files_by_run[run]]
        geo, kS = ds[0]["geo"], ds[0]["kS"]
        pool = bm.cat([d["pools"]["legacy"] for d in ds])
        uu, vv = bm.uv(pool, geo)
        u.append(uu), v.append(vv), s.append(pool["seed"].astype(np.int64)), e.append(ecal_of(pool, g_by_run[run]))
    u, v, s, e = map(np.concatenate, (u, v, s, e))
    ep = bm.e_elastic(bm.theta_deg(geo["x"], geo["y"], geo["hz"]), geo["ebeam"], bm.pool_mass(rt))
    r = e / ep[s]
    tref = float(np.mean(ep[np.unique(s)]))              # as fit_pc: before the r window
    m = np.isfinite(u) & np.isfinite(v) & np.isfinite(r) & (r >= bm.PC_R[0]) & (r <= bm.PC_R[1])
    u, v, r, s = u[m], v[m], r[m], s[m]
    out = {"none": (None, dict())}
    maps = {}
    for name, B in BINS.items():
        bm.PC_BINS = B
        R, n, filled = bm.fit_map(u, v, r, bm.PC_FIT[rt], tref, kS)
        maps[B] = R
        out[name] = (dict(R=R, own={}, glob=None), dict(own=0, fitted=int((~filled).sum()), total=B * B, glob="NA", own_failed=0))
    bm.PC_BINS = 12
    mod_idx = bm.groups(s)
    own8 = {mod: (bm.fit_module(u[idx], v[idx], r[idx]), idx.size) for mod, idx in mod_idx.items() if idx.size >= 1000}
    base = out["map12"][1]
    for name, N in OWN_N.items():
        own = {mod: c for mod, (c, n) in own8.items() if n >= N and c is not None}
        failed = sum(1 for _, (c, n) in own8.items() if n >= N and c is None)
        out[name] = (dict(R=maps[12], own=own, glob=None), base | dict(own=len(own), own_failed=failed))
    for form, k in TERMS.items():
        c = fit_poly(u, v, r, k)
        out[f"{form}_global"] = (dict(R=None, own={}, glob=c), dict(own=0, fitted=0, total=0, glob="pass" if c is not None else "FAIL", own_failed=0))
        own = {}
        failed = 0
        for mod, idx in mod_idx.items():
            if idx.size >= 3000:
                cc = fit_poly(u[idx], v[idx], r[idx], k)
                if cc is None:
                    failed += 1
                else:
                    own[mod] = cc
        out[f"{form}_pm"] = (dict(R=maps[12], own=own, glob=None), base | dict(own=len(own), own_failed=failed))
    out["_n_train"] = (None, dict(n=int(m.sum())))
    return out


# --- step 4: evaluation --------------------------------------------------------------------------------------------------------
def map_eval(G, u, v):
    B = G.shape[0]
    cu = -bm.PC_RANGE + (np.arange(B) + 0.5) * 2 * bm.PC_RANGE / B
    uc, vc = np.clip(u, cu[0], cu[-1]), np.clip(v, cu[0], cu[-1])
    iu = np.clip(np.searchsorted(cu, uc, side="right") - 1, 0, B - 2)
    iv = np.clip(np.searchsorted(cu, vc, side="right") - 1, 0, B - 2)
    tu, tv = (uc - cu[iu]) / (cu[iu + 1] - cu[iu]), (vc - cu[iv]) / (cu[iv + 1] - cu[iv])
    return ((1 - tu) * (1 - tv) * G[iu, iv] + tu * (1 - tv) * G[iu + 1, iv] + (1 - tu) * tv * G[iu, iv + 1]
            + tu * tv * G[iu + 1, iv + 1])


def r_eval(spec, pool, geo):
    u, v = bm.uv(pool, geo)
    if spec is None:
        return np.ones(u.size)
    if spec["glob"] is not None or spec["R"] is None:
        return bm.poly_r(spec["glob"], u, v) if spec["glob"] is not None else np.ones(u.size)
    R = map_eval(spec["R"], u, v)
    for mod, c in spec["own"].items():
        m = pool["seed"] == mod
        R[m] = bm.poly_r(c, u[m], v[m])
    return R


def eval_job(a):
    """One split: {(sample, variant, path): (res_cb, res_g2, mu_ep)}; a = (file, {sample: (gains, {variant: spec})})."""
    fn, corr = a
    d = load(fn)
    geo, kE, kS, mass = d["geo"], d["kE"], d["kS"], bm.pool_mass(d["rt"])
    out = {}
    for sample, (g, specs) in corr.items():
        ecal = {p: ecal_of(d["pools"][p], g) for p in PATHS}
        memo = {}
        for var, spec in specs.items():
            key = ("map", spec["R"].shape[0]) if spec is not None and spec["R"] is not None and spec["glob"] is None and not spec["own"] else (var,)
            if key not in memo:
                pools = {p: bm.nl_view(d["pools"][p], 0) | {"e": ecal[p] / r_eval(spec, d["pools"][p], geo)} for p in PATHS}
                res = bm.evaluate(pools, geo, mass, kE, kS)[0]
                memo[key] = {p: (res[p]["res_cb"], res[p]["res_g2"], res[p]["mu_ep"]) for p in PATHS}
            for p in PATHS:
                out[sample, var, p] = memo[key][p]
    print(f"  [eval] {Path(fn).name}", flush=True)
    return out


def se(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(x.std(ddof=1) / math.sqrt(x.size)) if x.size > 1 else float("nan")


def mean(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return float(x.mean()) if x.size else float("nan")


def main():
    t_start = time.monotonic()
    CACHE.mkdir(parents=True, exist_ok=True)
    info = {}                                           # rt -> {run: [(split, file, cache path)]}
    jobs = []
    for rt in RTS:
        for fn in sorted((SCAN / rt).glob("*.npz")):
            with np.load(fn, allow_pickle=False) as z:
                run, splits = int(z["run"]), [int(x) for x in z["splits"]]
            for sp in splits:
                jobs.append((str(fn), sp))
                info.setdefault(rt, {}).setdefault(run, []).append((sp, str(CACHE / f"{fn.stem}_{rt}_{sp}.pkl")))
    ctx = mp.get_context("fork")
    with ctx.Pool(NW) as pool:
        pool.map(cache_job, jobs, chunksize=1)
        t1 = time.monotonic()
        print(f"[cache] done {t1 - t_start:.0f} s", flush=True)
        # ordinals, training sets
        ordinal = {rt: {run: {sp: i for i, (sp, _) in enumerate(sorted(l))} for run, l in runs.items()} for rt, runs in info.items()}
        train = {}                                      # (rt, tset) -> {run: [cache files]}
        for rt, runs in info.items():
            for ts in ("all", "even", "odd"):
                train[rt, ts] = {run: [f for sp, f in l if ts == "all" or ordinal[rt][run][sp] % 2 == (ts == "odd")] for run, l in runs.items()}
        gk = [(rt, ts, run) for (rt, ts), d in train.items() for run in d]
        gres = dict(zip(gk, pool.map(gain_job, [(train[rt, ts][run], rt) for rt, ts, run in gk], chunksize=1)))
        failed = {k: v for k, v in gres.items() if v is None}
        for k in failed:
            print(f"[WARN] gain solve failed {k}, gains 1", flush=True)
        ones = np.ones(load(next(iter(info[RTS[0]].values()))[0][1])["geo"]["x"].size)
        gains = {k: (v if v is not None else ones) for k, v in gres.items()}
        t2 = time.monotonic()
        print(f"[gains] done {t2 - t1:.0f} s", flush=True)
        pk = [(rt, ts) for rt in RTS if rt not in BORROW for ts in ("all", "even", "odd")]
        pres = dict(zip(pk, pool.map(pc_job, [(rt, train[rt, ts], {r: gains[rt, ts, r] for r in train[rt, ts]}) for rt, ts in pk], chunksize=1)))
        for k in pk:
            pres[k].pop("_n_train")
        t3 = time.monotonic()
        print(f"[poscorr] done {t3 - t2:.0f} s", flush=True)
        # evaluation jobs: per split, in-sample and out-of-sample
        ej, keys = [], []
        for rt in RTS:
            src = BORROW.get(rt, rt)
            for run, l in info[rt].items():
                for sp, f in l:
                    oth = "odd" if ordinal[rt][run][sp] % 2 == 0 else "even"
                    corr = {}
                    for sample, ts in (("in", "all"), ("out", oth)):
                        pts = "all" if rt in BORROW else ts
                        corr[sample] = (gains[rt, ts, run], {v: pres[src, pts][v][0] for v in VARIANTS})
                    ej.append((f, corr))
                    keys.append((rt, run, sp, "even" if oth == "odd" else "odd"))
        eres = pool.map(eval_job, ej, chunksize=1)
        t4 = time.monotonic()
        print(f"[eval] done {t4 - t3:.0f} s", flush=True)
    pickle.dump(dict(keys=keys, eres=eres, pres={k: {v: d for v, (_, d) in x.items()} for k, x in pres.items()}, failed=list(failed)),
                open(CACHE / "results.pkl", "wb"))
    wall = time.monotonic() - t_start
    write_outputs(keys, eres, pres, failed, wall)


def write_outputs(keys, eres, pres, failed, wall):
    vals = {}                                           # (rt, var, path, sample) -> array [3 quantities, n splits]
    half = {}                                           # split index -> half label (the half the split belongs to)
    for rt in RTS:
        idx = [i for i, k in enumerate(keys) if k[0] == rt]
        half[rt] = np.array([keys[i][3] for i in idx])
        for var in VARIANTS:
            for p in PATHS:
                for sample in ("in", "out"):
                    vals[rt, var, p, sample] = np.array([eres[i][sample, var, p] for i in idx]).T
    rows = []
    for rt in RTS:
        cur = CURRENT[rt]
        for var in VARIANTS:
            for p in PATHS:
                for sample in ("in", "out"):
                    a = vals[rt, var, p, sample]
                    c = vals[rt, cur, p, sample]
                    gap = vals[rt, var, p, "in"][0] - vals[rt, var, p, "out"][0]
                    dg = lambda ts: pres[(BORROW.get(rt, rt), "all" if rt in BORROW else ts)][var][1]    # noqa: E731
                    ds = [dg("all")] if sample == "in" else [dg("even"), dg("odd")]
                    jn = lambda k: "|".join(str(d.get(k, "NA")) for d in ds)   # noqa: E731
                    gf = [sum(1 for f in failed if f[0] == rt and f[1] == ts) for ts in (("all",) if sample == "in" else ("even", "odd"))]
                    ho = half[rt] == "odd"
                    rows.append([rt, var, p, sample, a.shape[1], *[f"{x:.5g}" for i in range(3) for x in (mean(a[i]), se(a[i]))],
                                 *([f"{mean(a[0][ho]):.5g}", f"{mean(a[0][~ho]):.5g}"] if sample == "out" else ["NA", "NA"]),
                                 f"{mean(a[0] - c[0]):.5g}", f"{se(a[0] - c[0]):.5g}",
                                 *([f"{mean(gap):.5g}", f"{se(gap):.5g}"] if sample == "out" else ["NA", "NA"]),
                                 jn("own"), jn("own_failed"), jn("fitted"), jn("total"), jn("glob"), "|".join(map(str, gf))])
    hdr = ("run_type", "variant", "path", "sample", "n_splits", "res_cb", "res_cb_se", "res_g2", "res_g2_se", "mu_ep", "mu_ep_se",
           "res_cb_odd_half_splits", "res_cb_even_half_splits", "d_cb", "d_cb_se", "gap_cb", "gap_cb_se", "own_modules",
           "own_fits_failed", "map_bins_fitted", "map_bins_total", "global_fit", "gain_solves_failed")
    comment = ("# poscorr_comparison.py: gains on, Moller rejection off, non-linearity off (arm gp); res_cb, res_g2 (sigma/E, %) and mu_ep (%) are means "
               "and standard errors over splits.  sample in = gains and position correction fitted on all splits of the run type, evaluated "
               "on every split; sample out = even/odd cross fit (within each run ordinal split numbers 0, 1, 2, ...; fit on one half, evaluate "
               "the other half; every split has one out-of-sample value; mean +- standard error over all splits).  res_cb_odd_half_splits / "
               "res_cb_even_half_splits (out rows): mean over the odd-ordinal / even-ordinal splits (fit on the even / odd half).  d_cb, d_cb_se: "
               "paired per-split res_cb difference to the current choice of the run type at the same sample and path "
               "(0 on the current row itself; current = hybrid3000 for 2.24 and 3.49 GeV hydrogen, map12 for 0.73 GeV and X17).  gap_cb, gap_cb_se "
               "(out rows): paired per-split in minus out res_cb (negative: out-of-sample worse).  Splits whose fit is NaN (3p49GeV_H2 run 24937: "
               "gain solve fails, as in the benchmark report) count in n_splits but not in means.  Diagnostics: own_modules = seed modules with an own polynomial, own_fits_failed = modules "
               "above the count threshold whose fit failed, map_bins_fitted / total, global_fit pass or FAIL (FAIL: R = 1); out rows: "
               "even-half|odd-half training.  2p24GeV_Ta_X17: gains are its own; the position correction of every variant is the in-sample fit of the "
               "same variant on 2p24GeV_H2 (including its own polynomials), for both samples; the current code borrows only the map, so "
               f"X17 map12 is the current production.  Total wall time of the driver {wall:.0f} s.")
    OUT.mkdir(exist_ok=True)
    bm.write_tsv(OUT / "poscorr_comparison.tsv", hdr, rows, comment)
    write_pdf(vals)
    print(f"[done] wall {wall:.0f} s", flush=True)


LAB = {"none": "none", "map8": "bilinear map\n8 x 8", "map12": "bilinear map\n12 x 12", "map16": "bilinear map\n16 x 16",
       "map24": "bilinear map\n24 x 24", "hybrid1000": "hybrid, own polynomial\nfrom 1000 clusters",
       "hybrid3000": "hybrid, own polynomial\nfrom 3000 clusters", "hybrid10000": "hybrid, own polynomial\nfrom 10000 clusters",
       "even4_global": "even4, one for\nthe run type", "even6_global": "even6, one for\nthe run type", "even4_pm": "even4,\nper module",
       "even6_pm": "even6,\nper module"}
TITLE = {"0p73GeV_H2": "0.73 GeV, hydrogen", "2p24GeV_H2": "2.24 GeV, hydrogen", "3p49GeV_H2": "3.49 GeV, hydrogen",
         "3p49GeV_H2_old": "3.49 GeV, hydrogen, previous runs", "2p24GeV_Ta_X17": "2.24 GeV, tantalum (X17 runs),\ncorrection from 2.24 GeV hydrogen"}


def write_pdf(vals):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams.update({
        "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
        "xtick.labelsize": 12, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
        "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
        "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
    x = np.arange(len(VARIANTS))
    with PdfPages(OUT / "poscorr_comparison.pdf") as pdf:
        for p, name in zip(PATHS, ("legacy clustering", f"time-gated clustering, {GATE:g} ns")):
            fig, axs = plt.subplots(2, 3, figsize=(26, 15))
            for ax, rt in zip(axs.flat, RTS):
                for sample, mk, lab, off in (("in", "o", "in-sample", -0.12), ("out", "s", "out-of-sample", 0.12)):
                    m = [mean(vals[rt, v, p, sample][0]) for v in VARIANTS]
                    e = [se(vals[rt, v, p, sample][0]) for v in VARIANTS]
                    ax.errorbar(x + off, m, e, fmt=mk, capsize=3, ms=8, label=lab)
                ax.axvspan(VARIANTS.index(CURRENT[rt]) - 0.45, VARIANTS.index(CURRENT[rt]) + 0.45, color="0.85", zorder=0,
                           label="current choice")
                ax.set_xticks(x, [LAB[v] for v in VARIANTS], rotation=60, ha="right")
                ax.set_ylabel("Crystal Ball sigma/E (%)")
                ax.set_title(TITLE[rt], fontsize=15)
                ax.grid(alpha=0.3)
            ax = axs.flat[-1]
            ax.axis("off")
            h, l = axs.flat[0].get_legend_handles_labels()
            ax.legend(h, l, loc="center", frameon=False)
            ax.text(0.5, 0.2, "points: mean over splits, bars: standard error;\nout-of-sample: even/odd cross fit", ha="center", fontsize=14)
            fig.suptitle(f"Position-correction variants, {name}: gains on, Moller rejection and non-linearity off")
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)


if __name__ == "__main__":
    main()
