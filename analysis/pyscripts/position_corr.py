#!/usr/bin/env python3
"""position_corr.py — cluster-energy dependence on the impact position inside the seed crystal.

The map is built from a benchmark scan directory written with --write-cluster-tsv (columns
cluster_x_mm / cluster_y_mm, the C++ log-weighted centroid) laid out as
<scan>/scan_<run:06d>/split_<n:05d>/{result.tsv, result_window.tsv, result_clusters.tsv}.

  u = (cluster_x_mm - seed_x_mm) / size_x,  v = (cluster_y_mm - seed_y_mm) / size_y
  r = cluster_E_cal_MeV / E'(seed module)   (E' = eprime_MeV of result_window.tsv)

Pool = plot_module_resolution.py's pool (chosen path, seed r < --inner-radius, seed energy >=
header min_seed_energy, outside the 1-based --exclude-rowcol rectangle, finite cluster_E_cal).
Runs are grouped by driver run type (all seeds pooled); Ta (X17) is mapped and reported only.

Subcommands:
  map  --scan-dir DIR --out OUTDIR [--splits even|odd|all] [--path legacy|gated]
       [--bins 12 --range 0.6 --min-n 300] [--gauss-window 2sigma|fixed150]
       Per run type: bins x bins grid over u, v in [-range, range]; per bin the preferred fit
       (gate_scan R_PREF: gauss or cb, the benchmark's fit functions as plot_module_resolution
       uses them; Gaussian window per --gauss-window) of r in [0.85, 1.12] -> mu (= R), sigma, n, mu_err,
       and the median.  mu_err = bootstrap std of the refitted mu (--boot N refits of the
       bin's resampled events, default 30, seed 12345 + bin index; 0 = sigma/sqrt(n_in_fit)).
       --splits even|odd = the harness's definition (0-based position in each
       run's sorted stage-R / stage-RX pair list, gate_scan.r_pairs); all = every complete split
       (24386.00060 is only in 'all').  Writes map_<slug>_<set>.tsv/.png, sigma_<slug>_<set>.png,
       n_<slug>_<set>.png, rings_<slug>_<set>.tsv and updates map_summary.tsv (--path gated:
       OUTDIR gets the suffix _gated).
       `map --splits even|odd` also prints the pull RMS of (R_even - R_odd) / hypot(se) per run
       type when both maps exist in OUTDIR (error calibration; eo_check.tsv).
  fit  --map MAP_TSV --form bilinear|poly8|both --out JSON
       Writes the poscorr_v1 JSON (bilinear: binned R, sparse bins filled from populated
       neighbours, normalised to the population-weighted mean 1; poly8: weighted least squares
       of R = c0 (1 + c1 u^2 + c2 v^2 + c3 u^2 v^2 + c4 u^4 + c5 v^4 + c6 u + c7 v)).
       --form both writes JSON (bilinear) and <stem>_poly8.json.
  fit  --global --maps GLOB --out DIR
       Pools the per-type maps (Ta excluded) by inverse variance -> poscorr_global_<set>.json
       and consistency_global.tsv (verdict 'global' only if >= 95 % of populated bins agree
       within 2 sigma in every H2 and C run type, else 'per_run_type').
  fit  --per-module --scan-dir DIR --fallback FALLBACK_JSON --out OUT_JSON [--splits all|even|odd]
       [--slug SLUG] [--min-clusters 3000] [--clamp 0.6]
       Hybrid correction (poscorr_hybrid_v1; Study B form d): per seed module with >= --min-clusters
       training clusters (r in [0.85, 1.12], the `map` pool of the runs whose run-type slug is SLUG,
       default the fallback's slug) a poly8 fitted to 6x6-bin +-2 sigma truncated means (bins with
       >= 30 events, >= 10 bins), normalised to event mean 1, accepted if R in [0.9, 1.1] on a
       13x13 grid over [-clamp, clamp]^2; every other module uses the fallback poscorr_v1 JSON
       (embedded verbatim).  Also writes <stem>_modules.tsv (all seeding modules).

Usage:
  python3 position_corr.py map --scan-dir output/scan_20261025_fa --splits even --out output/poscorr_20261025/
  python3 position_corr.py fit --map output/poscorr_20261025/map_2p24GeV_H2_even.tsv --form both \\
      --out output/poscorr_20261025/2p24GeV_H2/poscorr_2p24GeV_H2_even.json
  python3 position_corr.py fit --global --maps "output/poscorr_20261025/map_*_all.tsv" \\
      --out output/poscorr_20261025/global/
  python3 position_corr.py fit --per-module --scan-dir output/scan_20261025_ref --splits even \\
      --fallback output/poscorr_20261025/2p24GeV_H2/poscorr_2p24GeV_H2_even.json --out hybrid_even.json
"""

from __future__ import annotations

import argparse
import concurrent.futures
import contextlib
import datetime
import hashlib
import glob
import io
import json
import math
import multiprocessing
import sys
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

_HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(_HERE))

import plot_module_resolution as R   # noqa: E402  (loader helpers; R.bench = the benchmark)
import gate_scan_20261022 as G       # noqa: E402  (r_pairs, R_PREF, run-type slug)

bench = R.bench
SCHEMA = "poscorr_v1"
R_LO, R_HI = 0.85, 1.12             # fit range of r = E_cal / E' (target_ref * [0.85, 1.12])
MIN_FIT = 30                        # fit a bin only with this many events in range
GATE_PCT = 0.3                      # RMS(R) minus statistical part, % (plan 5.2 gate)
REPORT_ONLY = ("2.24 GeV Ta (X17)",)
RING_BINS = 6


# ---- pool -------------------------------------------------------------------

def split_filter(scan_dir, run, which):
    """Set of split numbers of `run` for --splits even|odd (harness rule, G.r_pairs), None for
    all complete splits."""
    if which == "all":
        return None
    pairs = G.r_pairs("R") + G.r_pairs("RX")
    par = 0 if which == "even" else 1
    pos = {}
    keep = set()
    for r, s in pairs:      # 0-based position in each run's sorted pairs, as gate_scan cmd_run
        if r == run and pos.get(r, 0) % 2 == par:
            keep.add(s)
        pos[r] = pos.get(r, 0) + 1
    return keep


def load_pool(scan_dir, run, args, geo):
    """Pool events of a run: R.load_run's cuts, keeping the centroid.  Returns a dict of
    arrays (u, v, r, row, col, seed = seed module name) plus hdr and the (run, split) list, or None."""
    wanted = split_filter(scan_dir, run, args.splits)
    if wanted is not None and not wanted:
        return None
    splits = R.complete_splits(scan_dir / f"scan_{run:06d}", wanted)
    if not splits:
        return None
    cfg = "legacy" if args.path == "legacy" else "gated"
    wpath = "legacy" if args.path == "legacy" else "new"
    r2, (r0, r1, c0, c1) = args.inner_radius ** 2, args.exclude_rowcol
    out = {k: [] for k in ("u", "v", "r", "row", "col", "seed")}
    for d in splits:
        h = R.parse_header(d / "result.tsv")
        df = pd.read_csv(d / "result_clusters.tsv", sep="\t")
        if "cluster_x_mm" not in df.columns:
            sys.exit(f"{d}: result_clusters.tsv has no cluster_x_mm/cluster_y_mm columns "
                     f"(scan with --write-cluster-tsv from the centroid-column benchmark)")
        df = df[df.config_name == cfg]
        inrect = (df.seed_row + 1).between(r0, r1) & (df.seed_col + 1).between(c0, c1)  # 1-based
        df = df[(df.seed_x_mm ** 2 + df.seed_y_mm ** 2 < r2)
                & (df.seed_energy_MeV >= h["min_seed"]) & ~inrect
                & np.isfinite(df.cluster_E_cal_MeV)]
        w = pd.read_csv(d / "result_window.tsv", sep="\t", comment="#")
        eprime = w[w.path == wpath].set_index("module").eprime_MeV
        ep = df.seed_module_id.map(eprime).to_numpy(float)
        sx = df.seed_module_id.map(geo.sx).to_numpy(float)
        sy = df.seed_module_id.map(geo.sy).to_numpy(float)
        u = (df.cluster_x_mm.to_numpy() - df.seed_x_mm.to_numpy()) / sx
        v = (df.cluster_y_mm.to_numpy() - df.seed_y_mm.to_numpy()) / sy
        r = df.cluster_E_cal_MeV.to_numpy() / ep
        ok = np.isfinite(u) & np.isfinite(v) & np.isfinite(r) & (ep > 0)
        for k, a in (("u", u), ("v", v), ("r", r), ("row", df.seed_row.to_numpy() + 1),
                     ("col", df.seed_col.to_numpy() + 1), ("seed", df.seed_module_id.to_numpy())):
            out[k].append(a[ok])
    out = {k: np.concatenate(a) for k, a in out.items()}
    out["hdr"] = R.parse_header(splits[0] / "result.tsv")
    out["splits"] = [(run, int(d.name.split("_")[1])) for d in splits]
    return out


def pool_type(pools):
    pools = [p for p in pools if p]
    if not pools:
        return None
    d = {k: np.concatenate([p[k] for p in pools]) for k in ("u", "v", "r", "row", "col", "seed")}
    d["hdr"] = pools[0]["hdr"]
    d["splits"] = [s for p in pools for s in p["splits"]]
    return d


# ---- per-bin fits -----------------------------------------------------------

def fit_ratio(r, hdr, kind):
    """Preferred fit of r (peak near 1) via the benchmark functions on r * target_ref, as
    R.fit_one does; the Gaussian window follows hdr["gw"] (2sigma, the default, or fixed150,
    set by `map --gauss-window`).  Returns (mu, sigma, n_in_fit) in r units (NaN on failure)."""
    tref, bw, ks = hdr["target_ref"], hdr["bw"], hdr["kS"]
    x = r * tref
    lo, hi = tref * R_LO, tref * R_HI
    nan = float("nan")
    with contextlib.redirect_stderr(io.StringIO()), warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gmu, gsig, _, gn, edges, counts = bench.fit_peak(x, lo, hi, bw, ks)
        if kind == "gauss":
            if hdr.get("gw", "2sigma") == "2sigma":     # --gauss-window; CB seed below stays on the fixed fit
                g = bench.fit_gauss_reported(0.5 * (edges[:-1] + edges[1:]), counts, bw, ks, lo, hi,
                                             "2sigma", start=(gmu, gsig))
                return g[0] / tref, g[1] / tref, g[3]
            return gmu / tref, gsig / tref, gn
        cmu, csig, _, _, _, cn, *_ = bench.fit_peak_cb(
            x, lo, hi, bw, mu0_hint=None if math.isnan(gmu) else gmu, escale=ks)
    return (cmu / tref, csig / tref, cn) if np.isfinite(cmu) else (nan, nan, cn)


def boot_err(rb, hdr, kind, nboot, seed):
    """Bootstrap errors of the fitted (mu, sigma/mu): std over `nboot` refits (fit_ratio) of
    the bin's events resampled with replacement.  NaN if fewer than 5 refits succeed."""
    rs = np.random.default_rng(seed)
    mus, sos = [], []
    for _ in range(nboot):
        mu, sig, nf = fit_ratio(rb[rs.integers(0, len(rb), len(rb))], hdr, kind)
        if np.isfinite(mu) and np.isfinite(sig) and nf > 0:
            mus.append(mu)
            sos.append(sig / mu)
    return (float(np.std(mus, ddof=1)), float(np.std(sos, ddof=1))) if len(mus) >= 5 else (np.nan, np.nan)


def _boot_task(t):
    return boot_err(*t)


def fit_grid(u, v, r, hdr, kind, bins, rng, min_n, boot=0, pool=None):
    """bins x bins map over u, v in [-rng, rng]; arrays indexed [iu, iv].  boot > 0: mu_err
    and the sigma/mu error are bootstrap errors (seed 12345 + bin index) instead of
    sigma/sqrt(n_in_fit); only for bins with n >= min_n.  pool: optional executor map."""
    m = (r >= R_LO) & (r <= R_HI)
    u, v, r = u[m], v[m], r[m]
    w = 2 * rng / bins
    iu = np.floor((u + rng) / w).astype(int)
    iv = np.floor((v + rng) / w).astype(int)
    if bins % 2 == 0:       # exact 0 (single-module centroid) sits on a bin edge: alternate by event parity
        idx = np.arange(len(u))
        iu = np.where(u == 0, bins // 2 - 1 + idx % 2, iu)
        iv = np.where(v == 0, bins // 2 - 1 + (idx // 2) % 2, iv)
    inb = (iu >= 0) & (iu < bins) & (iv >= 0) & (iv < bins)
    flat = (iu * bins + iv)[inb]
    r = r[inb]
    order = np.argsort(flat, kind="stable")
    ids, start = np.unique(flat[order], return_index=True)
    groups = dict(zip(ids, np.split(order, start[1:])))
    g = {k: np.full((bins, bins), np.nan) for k in ("mu", "mu_err", "sigma", "median")}
    g["n"] = np.zeros((bins, bins), int)
    g["n_fit"] = np.zeros((bins, bins), int)
    flag = np.full((bins, bins), "empty", dtype=object)
    g["s_err_boot"] = np.full((bins, bins), np.nan)
    todo = []
    for b, idx in groups.items():
        i, j = divmod(int(b), bins)
        rb = r[idx]
        g["n"][i, j] = len(rb)
        g["median"][i, j] = float(np.median(rb))
        if len(rb) < MIN_FIT:
            flag[i, j] = "low_n"
            continue
        mu, sig, nfit = fit_ratio(rb, hdr, kind)
        g["n_fit"][i, j] = nfit
        g["mu"][i, j], g["sigma"][i, j] = mu, sig
        if np.isfinite(mu) and np.isfinite(sig) and nfit > 0:
            g["mu_err"][i, j] = sig / math.sqrt(nfit)
            flag[i, j] = "ok" if len(rb) >= min_n else "low_n"
            if boot > 0 and len(rb) >= min_n:
                todo.append(((i, j), (rb, hdr, kind, boot, 12345 + int(b))))
        else:
            flag[i, j] = "fit_fail"
    if todo:
        res = (pool.map(_boot_task, [t for _, t in todo], chunksize=4) if pool is not None
               else map(_boot_task, [t for _, t in todo]))
        for ((i, j), _), (e, se) in zip(todo, res):
            if np.isfinite(e):
                g["mu_err"][i, j] = e
                g["s_err_boot"][i, j] = se
    g["flag"] = flag
    g["sigma_over_mu"] = g["sigma"] / g["mu"]
    g["s_err"] = np.where(np.isfinite(g["s_err_boot"]), g["s_err_boot"],
                          g["sigma_over_mu"] / np.sqrt(2.0 * np.maximum(g["n_fit"], 1)))
    return g


def edges_of(bins, rng):
    return np.linspace(-rng, rng, bins + 1)


def populated(g):
    return (g["flag"] == "ok") & np.isfinite(g["mu"]) & np.isfinite(g["mu_err"])


def wstats(val, err, n, pop):
    """peak-to-peak, population-weighted RMS about the weighted mean, statistical part
    sqrt(<err^2>_w), and RMS with the statistical part removed (all as fractions)."""
    if not pop.any():
        return dict(ptp=np.nan, rms=np.nan, stat=np.nan, rms_corr=np.nan)
    x, e, w = val[pop], err[pop], n[pop].astype(float)
    mean = (w * x).sum() / w.sum()
    rms = math.sqrt((w * (x - mean) ** 2).sum() / w.sum())
    stat = math.sqrt((w * e ** 2).sum() / w.sum())
    return dict(ptp=float(x.max() - x.min()), rms=rms, stat=stat,
                rms_corr=math.sqrt(max(0.0, rms ** 2 - stat ** 2)))


# ---- plots ------------------------------------------------------------------

def draw_grid(a, rng, title, label, cmap, vmin, vmax, path):
    fig, ax = plt.subplots(figsize=(6.2, 5.4))
    cm = plt.get_cmap(cmap).copy()
    cm.set_bad("0.88")
    im = ax.imshow(np.ma.masked_invalid(a).T, origin="lower", extent=(-rng, rng, -rng, rng),
                   cmap=cm, vmin=vmin, vmax=vmax, aspect="equal")
    for x in (-0.5, 0.5):
        ax.axvline(x, color="k", lw=0.5, ls=":")
        ax.axhline(x, color="k", lw=0.5, ls=":")
    ax.set(xlabel="u = (x_cl - x_seed) / size_x", ylabel="v = (y_cl - y_seed) / size_y")
    fig.colorbar(im, ax=ax, shrink=0.85).set_label(label)
    ax.set_title(title, fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def make_plots(g, pop, rng, slug, rt, sset, fit, st, sst, outdir):
    tag = f"{rt} | {sset} | {fit}"
    mu = np.where(pop, g["mu"], np.nan)
    dev = 100 * (mu - 1)
    amp = max(0.5, float(np.nanmax(np.abs(dev))) if pop.any() else 0.5)
    draw_grid(dev, rng, f"{tag} | R - 1 (%) | RMS {100 * st['rms']:.2f}% stat {100 * st['stat']:.2f}%"
              f" ptp {100 * st['ptp']:.2f}%", "R - 1 (%)", "RdBu_r", -amp, amp,
              outdir / f"map_{slug}_{sset}.png")
    s = np.where(pop, 100 * g["sigma_over_mu"], np.nan)
    draw_grid(s, rng, f"{tag} | sigma/mu (%) | RMS {100 * sst['rms']:.2f}%", "sigma / mu (%)",
              "viridis", float(np.nanmin(s)) if pop.any() else 0, float(np.nanmax(s)) if pop.any() else 1,
              outdir / f"sigma_{slug}_{sset}.png")
    n = np.where(g["n"] > 0, g["n"], np.nan).astype(float)
    draw_grid(n, rng, f"{tag} | events per bin", "N (r in [0.85, 1.12])", "magma_r", 0,
              float(np.nanmax(n)) if np.isfinite(n).any() else 1, outdir / f"n_{slug}_{sset}.png")


# ---- tables -----------------------------------------------------------------

def grid_rows(g, edges, extra=None):
    rows = []
    nb = len(edges) - 1
    for i in range(nb):
        for j in range(nb):
            row = dict(iu=i, iv=j, u_lo=edges[i], u_hi=edges[i + 1], v_lo=edges[j], v_hi=edges[j + 1],
                       n=int(g["n"][i, j]), mu=g["mu"][i, j], mu_err=g["mu_err"][i, j],
                       sigma=g["sigma"][i, j], sigma_over_mu=g["sigma_over_mu"][i, j],
                       median=g["median"][i, j], flag=g["flag"][i, j], n_fit=int(g["n_fit"][i, j]))
            if extra:
                row.update({k: v[i, j] for k, v in extra.items()})
            rows.append(row)
    return rows


def write_map_tsv(path, rows, meta):
    with open(path, "w") as fh:
        for k, v in meta.items():
            fh.write(f"# {k}={v}\n")
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, float_format="%.6g", mode="a")


def read_map_tsv(path):
    meta = {}
    with open(path) as fh:
        for line in fh:
            if not line.startswith("#"):
                break
            k, _, v = line[1:].strip().partition("=")
            meta[k] = v
    df = pd.read_csv(path, sep="\t", comment="#")
    nb = int(meta["bins"])
    g = {}
    for k in ("n", "n_fit", "mu", "mu_err", "sigma", "sigma_over_mu", "median"):
        a = df[k].to_numpy().reshape(nb, nb)
        g[k] = a.astype(int) if k in ("n", "n_fit") else a.astype(float)
    g["flag"] = df["flag"].to_numpy().reshape(nb, nb).astype(object)
    return meta, g


def ring_masks(row, col):
    """Ring groups (1-based row/col): ring k = max(|row-17.5|, |col-17.5|) - 2.5; sides of
    ring 1 for non-corner modules; top = +y (row < 17.5), right = +x (col > 17.5)."""
    dr, dc = row - 17.5, col - 17.5
    k = np.rint(np.maximum(np.abs(dr), np.abs(dc)) - 2.5).astype(int)
    ar, ac = np.abs(dr), np.abs(dc)
    noncorner = (k == 1) & ((ar == 3.5) ^ (ac == 3.5))
    return {"all": np.ones(len(row), bool), "ring1": k == 1, "ring2": k == 2,
            "ring3-4": (k == 3) | (k == 4),
            "ring1_top": noncorner & (ar == 3.5) & (dr < 0),
            "ring1_bottom": noncorner & (ar == 3.5) & (dr > 0),
            "ring1_left": noncorner & (ac == 3.5) & (dc < 0),
            "ring1_right": noncorner & (ac == 3.5) & (dc > 0)}


def rings_table(d, kind, rng, min_n, path, boot=0, pool=None):
    """Map on a RING_BINS x RING_BINS grid per ring group; max |R_ring - R_all| in SE units
    (both maps normalised to their n-weighted mean over the bins populated in both; the group
    is a subset of 'all', so SE^2 = se_group^2 - se_all^2 = se_group^2 (1 - n_group/n_all) per bin;
    the bootstrap se_all is not used: its noise made the difference unstable)."""
    edges = edges_of(RING_BINS, rng)
    masks = ring_masks(d["row"], d["col"])
    grids = {k: fit_grid(d["u"][m], d["v"][m], d["r"][m], d["hdr"], kind, RING_BINS, rng, min_n, boot, pool)
             for k, m in masks.items()}
    def normed(g, both):
        w = g["n"][both].astype(float)
        nm = float((w * g["mu"][both]).sum() / w.sum())
        return g["mu"] / nm, g["mu_err"] / nm

    ga = grids["all"]
    pa = populated(ga)
    rows, summ = [], {}
    for k, g in grids.items():
        both = populated(g) & pa
        if both.any():
            mg, eg = normed(g, both)
            ma, ea = normed(ga, both)
            f = np.minimum(g["n"] / np.maximum(ga["n"], 1), 1.0)       # subset fraction per bin
            sd = np.maximum(eg * np.sqrt(1.0 - f), 1e-12)    # se_g^2 - se_all^2 with se^2 ~ 1/n
            z = np.where(both, np.abs(mg - ma) / sd, np.nan)
        else:
            mg = ma = z = np.full_like(g["mu"], np.nan)
        zmax = float(np.nanmax(z)) if both.any() else float("nan")
        summ[k] = (int(masks[k].sum()), int(both.sum()), zmax)
        for r in grid_rows(g, edges, {"mu_norm": mg, "mu_all_norm": ma, "dz": z}):
            rows.append(dict(group=k, n_group=int(masks[k].sum()), max_abs_dz=zmax, **r))
    pd.DataFrame(rows).to_csv(path, sep="\t", index=False, float_format="%.6g")
    return summ


def eo_check(out, slugs=None):
    """Empirical error calibration: pull of (R_even - R_odd) / hypot(se_even, se_odd) over the
    bins populated in both maps of OUTDIR, per run type; appends to eo_check.tsv."""
    rows = []
    for fe in sorted(Path(out).glob("map_*_even.tsv")):
        fo = fe.with_name(fe.name.replace("_even.tsv", "_odd.tsv"))
        if not fo.exists() or (slugs is not None and fe.name[4:-9] not in slugs):
            continue
        me, ge = read_map_tsv(fe)
        mo, go = read_map_tsv(fo)
        both = populated(ge) & populated(go)
        if not both.any():
            continue
        we = ge["n"][both].astype(float)
        wo = go["n"][both].astype(float)
        ne = float((we * ge["mu"][both]).sum() / we.sum())       # each normalised to its mean
        no = float((wo * go["mu"][both]).sum() / wo.sum())
        z = (ge["mu"][both] / ne - go["mu"][both] / no) / np.hypot(ge["mu_err"][both] / ne,
                                                                    go["mu_err"][both] / no)
        rows.append(dict(slug=me["slug"], run_type=me["run_type"], boot=me.get("boot", "0"),
                         n_bins=int(both.sum()), pull_rms=float(np.sqrt(np.mean(z ** 2))),
                         pull_mean=float(z.mean()), max_abs_pull=float(np.abs(z).max()),
                         frac_within_2=float((np.abs(z) <= 2).mean())))
        print(f"even vs odd [{me['slug']}]: {both.sum()} bins, pull RMS {rows[-1]['pull_rms']:.2f} "
              f"(mean {rows[-1]['pull_mean']:+.2f}, max {rows[-1]['max_abs_pull']:.1f}, "
              f"{100 * rows[-1]['frac_within_2']:.0f}% within 2) boot={rows[-1]['boot']}")
    if rows:
        pd.DataFrame(rows).to_csv(Path(out) / "eo_check.tsv", sep="\t", index=False, float_format="%.4f")
    return len(rows)


# ---- map command ------------------------------------------------------------

SUMMARY_COLS = ["slug", "set", "run_type", "path", "fit", "n_runs", "n_splits", "n_events",
                "n_populated_bins", "R_ptp_pct", "R_rms_pct", "R_stat_pct", "R_rms_stat_removed_pct",
                "s_ptp_pct", "s_rms_pct", "s_stat_pct", "s_rms_stat_removed_pct",
                "s_edge_pct", "s_centre_pct", "gate_below_0p3pct", "report_only"]


def update_summary(path, row):
    old = pd.read_csv(path, sep="\t", dtype={"set": str}) if path.exists() else pd.DataFrame(columns=SUMMARY_COLS)
    old = old[~((old.slug == row["slug"]) & (old["set"] == row["set"]) & (old.path == row["path"]))]
    new = pd.concat([old, pd.DataFrame([row])], ignore_index=True)[SUMMARY_COLS]
    new.to_csv(path, sep="\t", index=False, float_format="%.4f")


def cmd_map(a):
    scan_dir = Path(a.scan_dir)
    out = Path(a.out)
    if a.path == "gated":
        out = out.with_name(out.name + "_gated")
    out.mkdir(parents=True, exist_ok=True)
    runs = sorted(int(p.name.split("_")[1]) for p in scan_dir.glob("scan_*") if p.is_dir())
    types = {}
    for run in runs:
        rt = G._run_type(run)
        if rt in G.R_PREF:
            types.setdefault(rt, []).append(run)
        else:
            print(f"[skip] run {run}: run type {rt!r} has no preferred fit")
    edges = edges_of(a.bins, a.range)
    gates = {}
    pool = None
    if a.workers > 1 and a.boot > 0:     # fork: workers inherit the imported benchmark
        pool = concurrent.futures.ProcessPoolExecutor(
            max_workers=min(6, a.workers), mp_context=multiprocessing.get_context("fork"))
    slugs = []
    for rt in sorted(types, key=lambda t: list(G.R_PREF).index(t)):
        t0 = time.time()
        rs = types[rt]
        slug = G.runtype_slug(rs[0])
        geo = R.load_geometry(rs[0])
        d = pool_type([load_pool(scan_dir, run, a, geo) for run in rs])
        if d is None:
            print(f"[skip] {rt}: no complete split in set {a.splits}")
            continue
        d["hdr"]["gw"] = a.gauss_window
        slugs.append(slug)
        kind = G.R_PREF[rt]
        g = fit_grid(d["u"], d["v"], d["r"], d["hdr"], kind, a.bins, a.range, a.min_n, a.boot, pool)
        pop = populated(g)
        mu, n = g["mu"], g["n"]
        st = wstats(mu, g["mu_err"], n, pop)
        sst = wstats(g["sigma_over_mu"], g["s_err"], g["n_fit"], pop)
        cen = 0.5 * (edges[:-1] + edges[1:])
        rad = np.maximum(np.abs(cen)[:, None], np.abs(cen)[None, :])
        def wsig(sel):
            sel = sel & pop
            return float(np.average(g["sigma_over_mu"][sel], weights=g["n_fit"][sel])) if sel.any() else float("nan")
        s_edge, s_cen = wsig(rad > 0.4), wsig(rad < 0.2)
        ro = rt in REPORT_ONLY
        gate = bool(st["rms_corr"] * 100 < GATE_PCT)      # False for an empty map (NaN)
        gates[rt] = (gate, ro)
        n_ev = len(d["r"])
        splits = ",".join(f"{r}.{s}" for r, s in d["splits"])
        meta = dict(schema="poscorr_map_v1", scan_dir=scan_dir, run_type=rt, slug=slug, set=a.splits,
                    path=a.path, fit=kind, gauss_window=a.gauss_window, bins=a.bins, range=a.range, min_n=a.min_n,
                    boot=a.boot, inner_radius=a.inner_radius, exclude_rowcol=" ".join(map(str, a.exclude_rowcol)),
                    n_events=n_ev, n_splits=len(d["splits"]), target_ref=d["hdr"]["target_ref"],
                    splits=splits)
        write_map_tsv(out / f"map_{slug}_{a.splits}.tsv", grid_rows(g, edges), meta)
        make_plots(g, pop, a.range, slug, rt, a.splits, kind, st, sst, out)
        summ = rings_table(d, kind, a.range, a.min_n, out / f"rings_{slug}_{a.splits}.tsv", a.boot, pool)
        update_summary(out / "map_summary.tsv", {
            "slug": slug, "set": a.splits, "run_type": rt, "path": a.path, "fit": kind,
            "n_runs": len({s[0] for s in d["splits"]}), "n_splits": len(d["splits"]), "n_events": n_ev,
            "n_populated_bins": int(pop.sum()), "R_ptp_pct": 100 * st["ptp"], "R_rms_pct": 100 * st["rms"],
            "R_stat_pct": 100 * st["stat"], "R_rms_stat_removed_pct": 100 * st["rms_corr"],
            "s_ptp_pct": 100 * sst["ptp"], "s_rms_pct": 100 * sst["rms"], "s_stat_pct": 100 * sst["stat"],
            "s_rms_stat_removed_pct": 100 * sst["rms_corr"], "s_edge_pct": 100 * s_edge,
            "s_centre_pct": 100 * s_cen, "gate_below_0p3pct": gate if np.isfinite(st["rms_corr"]) else "n/a",
            "report_only": "report only" if ro else ""})
        print(f"{rt} [{a.splits}, {a.path}, {kind}]: {len(d['splits'])} splits, {n_ev} events, "
              f"{int(pop.sum())}/{a.bins ** 2} bins populated; R ptp {100 * st['ptp']:.2f}% "
              f"RMS {100 * st['rms']:.3f}% stat {100 * st['stat']:.3f}% -> stat removed "
              f"{100 * st['rms_corr']:.3f}%{' (report only)' if ro else ''}; sigma/mu RMS "
              f"{100 * sst['rms']:.3f}% (stat-removed {100 * sst['rms_corr']:.3f}%), edge "
              f"{100 * s_edge:.2f}% vs centre {100 * s_cen:.2f}%")
        print("  rings (n, bins, max|dR|/SE): " + "; ".join(
            f"{k} {v[0]}, {v[1]}, {v[2]:.1f}" for k, v in summ.items())
              + f"  [{time.time() - t0:.0f} s, boot {a.boot}]")
    if pool is not None:
        pool.shutdown()
    if a.splits in ("even", "odd"):
        eo_check(out, slugs)
    chk = {t: v[0] for t, v in gates.items() if not v[1]}
    if chk:
        print(f"GATE (RMS of R, statistical part removed, < {GATE_PCT}% in every H2 and C run type): "
              + ("PASS - skip fit/scan, go to the checkpoint" if all(chk.values())
                 else "FAIL - correction worth fitting") + "  " + str({t: bool(v) for t, v in chk.items()}))
    return 0


# ---- fit command ------------------------------------------------------------

def fill_neighbours(val, pop):
    """Fill unpopulated bins with the mean of the 8 neighbours that already have a value,
    simultaneous passes until the grid is full."""
    val = np.where(pop, val, np.nan)
    filled = np.zeros(val.shape, bool)
    nb = val.shape[0]
    if not pop.any():
        sys.exit("no populated bin: cannot build a correction")
    while np.isnan(val).any():
        new = val.copy()
        for i in range(val.shape[0]):
            for j in range(val.shape[1]):
                if np.isnan(val[i, j]):
                    blk = val[max(0, i - 1):i + 2, max(0, j - 1):j + 2]
                    if np.isfinite(blk).any():
                        new[i, j] = np.nanmean(blk)
                        filled[i, j] = True
        val = new
    return val, filled


def fill_grid(g):
    """Normalised, NaN-free (R, se, filled) from a map: bilinear table."""
    pop = populated(g)
    Rv, filled = fill_neighbours(g["mu"], pop)
    se = np.where(pop, g["mu_err"], np.nan)
    se_fill = float(np.nanmax(se)) if pop.any() else 1.0
    se = np.where(np.isfinite(se), se, se_fill)
    n = g["n"].astype(float)
    if n.sum() == 0:
        sys.exit("empty map")
    norm = float((n * Rv).sum() / n.sum())
    return Rv / norm, se / norm, filled, norm, pop


def poly_design(u, v):
    u, v = np.asarray(u, float), np.asarray(v, float)
    return np.stack([np.ones_like(u), u ** 2, v ** 2, u ** 2 * v ** 2, u ** 4, v ** 4, u, v], axis=-1)


def poly_fit(Rv, se, n, pop, uc, vc):
    """Weighted LSQ (1/se^2, populated bins) of R = sum a_i T_i(u, v); c0 = a0, ci = ai/a0,
    then c0 rescaled so that the n-weighted mean of the poly over all bins is 1.
    Returns (c, norm_poly, residual rms vs the map, chi2/ndf)."""
    U, V = np.meshgrid(uc, vc, indexing="ij")
    X = poly_design(U.ravel(), V.ravel())
    m = pop.ravel()
    w = 1.0 / se.ravel()[m]
    coef, *_ = np.linalg.lstsq(X[m] * w[:, None], Rv.ravel()[m] * w, rcond=None)
    a0 = coef[0]
    c = np.concatenate([[a0], coef[1:] / a0])
    pv = (X * np.concatenate([[1.0], c[1:]])).sum(axis=1) * c[0]
    nn = n.ravel().astype(float)
    norm_p = float((nn * pv).sum() / nn.sum())
    c[0] /= norm_p
    pv = pv / norm_p
    res = pv[m] - Rv.ravel()[m]
    wn = nn[m]
    rms = math.sqrt((wn * res ** 2).sum() / wn.sum())
    chi2 = float(((res / se.ravel()[m]) ** 2).sum() / max(1, m.sum() - 8))
    return c, norm_p, rms, chi2


def poly_dict(c):
    return {f"c{i}": float(c[i]) for i in range(8)}


def build_json(form, rt, slug, edges_u, edges_v, Rv, n, se, filled, norm, source, c=None, norm_p=None):
    gi = G.git_info()
    js = {"schema": SCHEMA, "form": form, "run_type": rt, "slug": slug,
          "u_edges": [float(x) for x in edges_u], "v_edges": [float(x) for x in edges_v],
          "R": Rv.tolist(), "n": n.astype(int).tolist(), "se": se.tolist(),
          "filled": filled.astype(int).tolist(), "norm": float(norm),
          "poly": poly_dict(c) if c is not None else None,
          "source": source, "created": datetime.datetime.now().isoformat(timespec="seconds"),
          "git": f"{gi['head']}{'+dirty' if gi['dirty'] else ''}" if gi["head"] else None}
    if norm_p is not None:
        js["poly_norm"] = float(norm_p)
    if not np.isfinite(np.asarray(Rv)).all():
        sys.exit("internal error: NaN in R")
    return js


def write_json(path, js):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as fh:
        json.dump(js, fh, indent=1)
    print(f"wrote {path}")


def source_of(meta, map_path):
    sp = [[int(x.split(".")[0]), int(x.split(".")[1])] for x in meta["splits"].split(",") if x]
    return {"scan_dir": meta["scan_dir"], "splits": sp, "fit": meta["fit"],
            "n_events": int(meta["n_events"]), "map_tsv": str(map_path), "set": meta["set"],
            "path": meta["path"]}


def cmd_fit_map(a):
    meta, g = read_map_tsv(a.map)
    nb, rng = int(meta["bins"]), float(meta["range"])
    edges = edges_of(nb, rng)
    cen = 0.5 * (edges[:-1] + edges[1:])
    Rv, se, filled, norm, pop = fill_grid(g)
    src = source_of(meta, a.map)
    out = Path(a.out)
    print(f"{meta['run_type']} [{meta['set']}]: {int(pop.sum())}/{nb * nb} populated bins, "
          f"{int(filled.sum())} filled, norm {norm:.6f}")
    if a.form in ("bilinear", "both"):
        write_json(out, build_json("bilinear", meta["run_type"], meta["slug"], edges, edges, Rv,
                                   g["n"], se, filled, norm, src))
    if a.form in ("poly8", "both"):
        c, norm_p, rms, chi2 = poly_fit(Rv, se, g["n"], pop, cen, cen)
        print(f"  poly8: c = {np.array2string(c, precision=4)}; residual RMS vs binned map "
              f"{100 * rms:.3f}% (n-weighted), chi2/ndf {chi2:.2f}")
        path = out if a.form == "poly8" else out.with_name(out.stem + "_poly8" + out.suffix)
        write_json(path, build_json("poly8", meta["run_type"], meta["slug"], edges, edges, Rv,
                                    g["n"], se, filled, norm, src, c=c, norm_p=norm_p))
    return 0


def cmd_fit_global(a):
    paths = sorted(glob.glob(a.maps))
    if not paths:
        sys.exit(f"no maps match {a.maps}")
    maps = []
    for p in paths:
        meta, g = read_map_tsv(p)
        maps.append((meta, g, p))
    sets = {m[0]["set"] for m in maps}
    if len(sets) != 1:
        sys.exit(f"maps of several split sets: {sorted(sets)}")
    sset = sets.pop()
    nb, rng = int(maps[0][0]["bins"]), float(maps[0][0]["range"])
    if any(int(m[0]["bins"]) != nb or float(m[0]["range"]) != rng for m in maps):
        sys.exit("maps differ in bins/range")
    edges = edges_of(nb, rng)
    cen = 0.5 * (edges[:-1] + edges[1:])
    per = {}
    for meta, g, p in maps:
        pop = populated(g)
        if not pop.any():
            continue
        w = g["n"][pop].astype(float)
        norm = float((w * g["mu"][pop]).sum() / w.sum())
        per[meta["run_type"]] = dict(meta=meta, g=g, pop=pop, R=g["mu"] / norm, se=g["mu_err"] / norm, path=p)
    pool_types = [t for t in per if t not in REPORT_ONLY]
    num = np.zeros((nb, nb))
    den = np.zeros((nb, nb))
    nsum = np.zeros((nb, nb))
    for t in pool_types:
        pt = per[t]
        wv = np.where(pt["pop"], 1.0 / np.where(pt["pop"], pt["se"], 1.0) ** 2, 0.0)
        num += wv * np.where(pt["pop"], pt["R"], 0.0)
        den += wv
        nsum += np.where(pt["pop"], pt["g"]["n"], 0)
    gpop = den > 0
    Rg = np.where(gpop, num / np.where(gpop, den, 1.0), np.nan)
    seg = np.where(gpop, 1.0 / np.sqrt(np.where(gpop, den, 1.0)), np.nan)
    Rf, filled = fill_neighbours(Rg, gpop)
    se_fill = float(np.nanmax(seg))
    seg_f = np.where(gpop, seg, se_fill)
    nn = nsum.astype(int)
    normg = float((nn * Rf).sum() / max(nn.sum(), 1))
    Rf, seg_f = Rf / normg, seg_f / normg
    Rg_n = Rg / normg
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    src = {"scan_dir": per[pool_types[0]]["meta"]["scan_dir"], "maps": [per[t]["path"] for t in pool_types],
           "run_types": pool_types, "set": sset, "fit": "per-type preferred",
           "splits": sorted({(int(x.split('.')[0]), int(x.split('.')[1])) for t in pool_types
                             for x in per[t]["meta"]["splits"].split(",") if x}),
           "n_events": int(sum(int(per[t]["meta"]["n_events"]) for t in pool_types)),
           "map_tsv": "global (inverse-variance pool)"}
    src["splits"] = [list(s) for s in src["splits"]]
    for form in (("bilinear", "poly8") if a.form == "both" else (a.form,)):
        kw = {}
        if form == "poly8":
            c, norm_p, rms, chi2 = poly_fit(Rf, seg_f, nn, gpop, cen, cen)
            print(f"global poly8: residual RMS vs pooled map {100 * rms:.3f}% chi2/ndf {chi2:.2f}")
            kw = dict(c=c, norm_p=norm_p)
        sfx = "" if form == "bilinear" else "_poly8"
        write_json(out / f"poscorr_global_{sset}{sfx}.json",
                   build_json(form, "global", "global", edges, edges, Rf, nn, seg_f, filled, normg, src, **kw))
    rows = []
    ok_all = True
    for t, pt in per.items():
        ro = t in REPORT_ONLY
        if ro:      # not pooled: compare with the full global
            both = pt["pop"] & gpop
            ref, sref = Rg_n, seg / normg
            Rt, st_ = pt["R"], pt["se"]
        else:       # leave-one-out global: the pool without this type
            wt = np.where(pt["pop"], 1.0 / np.where(pt["pop"], pt["se"], 1.0) ** 2, 0.0)
            dl = den - wt
            both = pt["pop"] & (dl > 1e-9 * np.maximum(den, 1e-300))
            ref = np.where(both, (num - wt * np.where(pt["pop"], pt["R"], 0.0)) / np.where(both, dl, 1.0), np.nan)
            sref = np.where(both, 1.0 / np.sqrt(np.where(both, dl, 1.0)), np.nan)
            Rt, st_ = pt["R"], pt["se"]
            if both.any():
                w = pt["g"]["n"][both].astype(float)
                nt = float((w * Rt[both]).sum() / w.sum())
                nr = float((w * ref[both]).sum() / w.sum())
                Rt, st_, ref, sref = Rt / nt, st_ / nt, ref / nr, sref / nr
        z = np.abs(Rt - ref) <= 2 * np.hypot(st_, sref)
        npass, nb_ = int((z & both).sum()), int(both.sum())
        frac = npass / nb_ if nb_ else float("nan")
        ro = t in REPORT_ONLY
        if not ro and not (frac >= 0.95):
            ok_all = False
        rows.append(dict(set=sset, run_type=t, slug=pt["meta"]["slug"], n_bins=nb_, n_pass=npass,
                         frac_pass=frac, pass95=bool(frac >= 0.95),
                         pooled_into_global=not ro, reference="full global" if ro else "leave-one-out",
                         report_only="report only" if ro else ""))
    verdict = "global" if ok_all else "per_run_type"
    for r in rows:
        r["verdict"] = verdict
    with open(out / "consistency_global.tsv", "w") as fh:
        fh.write("# reference: pooled types vs the leave-one-out global (pool without the type); "
                 "report-only types vs the full global; 2*hypot(se) test\n")
    pd.DataFrame(rows).to_csv(out / "consistency_global.tsv", sep="\t", index=False, float_format="%.4f",
                              mode="a")
    for r in rows:
        print(f"  {r['run_type']:<20s} {r['n_pass']}/{r['n_bins']} bins within 2 sigma ({100 * r['frac_pass']:.1f}%)"
              f"{'  [report only]' if r['report_only'] else ''}")
    print(f"VERDICT: {verdict} (leave-one-out; global requires >= 95% in every H2 and C run type; Ta excluded)")
    return 0


# ---- hybrid per-module poly8 (Study B form d) --------------------------------

HYB_SCHEMA = "poscorr_hybrid_v1"
HYB_RMIN, HYB_RMAX = 0.85, 1.12     # training r window (before the per-module cluster count)
HYB_GRID = 6                        # 6x6 bins over [-clamp, clamp]^2
HYB_MIN_BIN = 30                    # events per populated bin
HYB_MIN_BINS = 10                   # populated bins needed for the poly8 (8 parameters)
HYB_RANGE = (0.9, 1.1)              # accepted R range of a module fit over the clamp square
HYB_TRUNC_K = 2.0
_PHI = math.exp(-HYB_TRUNC_K ** 2 / 2) / math.sqrt(2 * math.pi)
_CTRUNC = 1 - 2 * HYB_TRUNC_K * _PHI / math.erf(HYB_TRUNC_K / math.sqrt(2))   # var of the +-2 sigma truncated Gaussian / sigma^2


def poly_basis(u, v):
    return np.stack([np.ones_like(u), u ** 2, v ** 2, u ** 2 * v ** 2, u ** 4, v ** 4, u, v], axis=-1)


def poly_R(c, u, v):
    """c0 (1 + c1 u^2 + c2 v^2 + c3 u^2 v^2 + c4 u^4 + c5 v^4 + c6 u + c7 v); c (..., 8)."""
    c = np.asarray(c, float)
    return c[..., 0] * (1 + c[..., 1] * u ** 2 + c[..., 2] * v ** 2 + c[..., 3] * u ** 2 * v ** 2
                        + c[..., 4] * u ** 4 + c[..., 5] * v ** 4 + c[..., 6] * u + c[..., 7] * v)


def poly_R_clamped(c, u, v, clamp):
    return poly_R(c, np.clip(u, -clamp, clamp), np.clip(v, -clamp, clamp))


def trunc_mean(r, k=HYB_TRUNC_K, max_iter=20):
    """Iterative +-2 sigma truncated mean: start at the median and 1.4826 MAD, then mean and the
    truncation-corrected sigma of the events in [mu-2s, mu+2s] until stable.  Returns
    (mu, sigma, n_in, converged) or None."""
    r = np.asarray(r, float)
    if len(r) < 5:
        return None
    mu = float(np.median(r))
    sig = 1.4826 * float(np.median(np.abs(r - mu)))
    if not sig > 0:
        return None
    n_in = 0
    conv = False
    for _ in range(max_iter):
        sel = np.abs(r - mu) <= k * sig
        n_in = int(sel.sum())
        if n_in < 5:
            return None
        mu2 = float(r[sel].mean())
        sig2 = float(r[sel].std(ddof=1)) / math.sqrt(_CTRUNC)
        done = abs(mu2 - mu) < 1e-7 and abs(sig2 - sig) < 1e-4 * sig
        mu, sig = mu2, sig2
        if done:
            conv = True
            break
    sel = np.abs(r - mu) <= k * sig
    return mu, sig, int(sel.sum()), conv


def bin_index(u, v, bins, rng):
    """fit_grid's binning (exact 0 alternates over the two central bins by event parity)."""
    w = 2 * rng / bins
    iu = np.floor((u + rng) / w).astype(int)
    iv = np.floor((v + rng) / w).astype(int)
    if bins % 2 == 0:
        idx = np.arange(len(u))
        iu = np.where(u == 0, bins // 2 - 1 + idx % 2, iu)
        iv = np.where(v == 0, bins // 2 - 1 + (idx // 2) % 2, iv)
    inb = (iu >= 0) & (iu < bins) & (iv >= 0) & (iv < bins)
    return iu, iv, inb


def grid_estimates(u, v, r, bins, rng, min_n):
    """Truncated-mean R and SE per bin (n >= min_n).  Returns centres (uc, vc), R, SE, n arrays
    for the populated bins, or None."""
    iu, iv, inb = bin_index(u, v, bins, rng)
    flat = np.where(inb, iu * bins + iv, -1)
    cen = -rng + (np.arange(bins) + 0.5) * 2 * rng / bins
    order = np.argsort(flat, kind="stable")
    ids, start = np.unique(flat[order], return_index=True)
    out = []
    for b, idx in zip(ids, np.split(order, start[1:])):
        if b < 0 or len(idx) < min_n:
            continue
        t = trunc_mean(r[idx])
        if t is None:
            continue
        mu, sig, nin, _ = t
        i, j = divmod(int(b), bins)
        out.append((cen[i], cen[j], mu, sig / math.sqrt(nin), len(idx)))
    if not out:
        return None
    a = np.array(out)
    return a[:, 0], a[:, 1], a[:, 2], a[:, 3], a[:, 4]


def fit_poly8_bins(uc, vc, Rv, se):
    X = poly_basis(uc, vc)
    w = 1.0 / se
    coef, *_ = np.linalg.lstsq(X * w[:, None], Rv * w, rcond=None)
    a0 = coef[0]
    return np.concatenate([[a0], coef[1:] / a0])


def normalise(p, u, v, clamp):
    """Rescale c0 so the event mean of R(u, v) over the given points is 1."""
    p = np.array(p, float)
    p[0] /= poly_R_clamped(p, u, v, clamp).mean()
    return p


def range_ok(p, clamp):
    g = np.linspace(-clamp, clamp, 13)
    U, V = np.meshgrid(g, g, indexing="ij")
    Rv = poly_R_clamped(p, U.ravel(), V.ravel(), clamp)
    return (bool(np.isfinite(Rv).all() and Rv.min() >= HYB_RANGE[0] and Rv.max() <= HYB_RANGE[1]),
            float(Rv.min()), float(Rv.max()))


def fit_events(u, v, r, clamp):
    """6x6 truncated-mean estimator + poly8 fit + event normalisation.  Returns (p, reason)."""
    g = grid_estimates(u, v, r, HYB_GRID, clamp, HYB_MIN_BIN)
    if g is None:
        return None, "no_bins"
    uc, vc, Rv, se, n = g
    if len(Rv) < HYB_MIN_BINS:
        return None, "few_bins"
    try:
        p = fit_poly8_bins(uc, vc, Rv, se)
    except Exception:
        return None, "fit_fail"
    if p is None or not np.isfinite(p).all():
        return None, "fit_fail"
    p = normalise(p, u, v, clamp)
    ok, lo, hi = range_ok(p, clamp)
    if not ok:
        return None, f"range({lo:.3f},{hi:.3f})"
    return p, "ok"


def lifts(p, clamp):
    """Edge lift (mean R at (+-0.5,0),(0,+-0.5)) - R(0,0) and corner lift (+-0.5,+-0.5) - R(0,0), in %."""
    p = np.asarray(p, float)
    e = [(0.5, 0), (-0.5, 0), (0, 0.5), (0, -0.5)]
    c = [(0.5, 0.5), (0.5, -0.5), (-0.5, 0.5), (-0.5, -0.5)]
    f = lambda pts: float(np.mean([poly_R_clamped(p, np.array(a), np.array(b), clamp) for a, b in pts]))
    r0 = float(poly_R_clamped(p, np.array(0.0), np.array(0.0), clamp))
    return 100 * (f(e) - r0), 100 * (f(c) - r0)


def cmd_fit_hybrid(a):
    if not (a.scan_dir and a.fallback):
        sys.exit("fit --per-module needs --scan-dir and --fallback")
    scan_dir = Path(a.scan_dir)
    fb_path = Path(a.fallback)
    fb_bytes = fb_path.read_bytes()
    fb = json.loads(fb_bytes)
    slug = a.slug or fb["slug"]
    runs = [run for run in sorted(int(p.name.split("_")[1]) for p in scan_dir.glob("scan_*") if p.is_dir())
            if G.runtype_slug(run) == slug]
    if not runs:
        sys.exit(f"no run of {scan_dir} has run-type slug {slug}")
    rt = G._run_type(runs[0])
    geo = R.load_geometry(runs[0])
    d = pool_type([load_pool(scan_dir, run, a, geo) for run in runs])
    if d is None:
        sys.exit(f"no complete split in set {a.splits} for {slug}")
    m = np.isfinite(d["u"]) & np.isfinite(d["v"]) & np.isfinite(d["r"]) & (d["r"] >= HYB_RMIN) & (d["r"] <= HYB_RMAX)
    u, v, r, seed = d["u"][m], d["v"][m], d["r"][m], d["seed"][m]
    n_ev = len(r)
    names, inv = np.unique(seed, return_inverse=True)
    order = np.argsort(inv, kind="stable")
    cnt = np.bincount(inv, minlength=len(names))
    bounds = np.concatenate([[0], np.cumsum(cnt)])
    mods, recs, failed = {}, [], {}
    n_own_ev = 0
    n_elig = 0
    for i, name in enumerate(names):
        ev = order[bounds[i]:bounds[i + 1]]
        mid = int(geo.module_id[name])
        rec = dict(module=name, id=mid, n=int(cnt[i]), own=0, reason=f"n<{a.min_clusters}")
        if cnt[i] >= a.min_clusters:
            n_elig += 1
            p, reason = fit_events(u[ev], v[ev], r[ev], a.clamp)
            rec["reason"] = reason
            if p is None:
                failed[name] = reason
            else:
                le, lc = lifts(p, a.clamp)
                rec.update(own=1, L_e=le, L_c=lc, **{f"c{k}": float(x) for k, x in enumerate(p)})
                mods[name] = {"id": mid, "n": int(cnt[i]), "c": [float(x) for x in p], "L_e": le, "L_c": lc}
                n_own_ev += int(cnt[i])
        recs.append(rec)
    nseed, n_own = len(names), len(mods)
    stats = {"seeding_modules": nseed, "eligible": n_elig, "own": n_own,
             "own_frac_modules": n_own / nseed, "own_frac_events": n_own_ev / n_ev,
             "n_events": n_ev, "failed": failed}
    gi = G.git_info()
    js = {"schema": HYB_SCHEMA, "form": "hybrid", "run_type": rt, "slug": slug,
          "module_form": "poly8", "clamp": float(a.clamp), "min_clusters": int(a.min_clusters),
          "fallback": fb, "fallback_path": str(a.fallback),
          "fallback_sha1": hashlib.sha1(fb_bytes).hexdigest(),
          "modules": dict(sorted(mods.items(), key=lambda kv: kv[1]["id"])), "stats": stats,
          "estimator": {"grid": HYB_GRID, "min_bin": HYB_MIN_BIN, "trunc_k": HYB_TRUNC_K,
                        "r_range": [HYB_RMIN, HYB_RMAX], "min_bins": HYB_MIN_BINS,
                        "range_check": list(HYB_RANGE), "norm": "module event mean"},
          "source": {"scan_dir": str(a.scan_dir), "splits": [list(s) for s in d["splits"]],
                     "set": a.splits, "path": a.path, "n_events": n_ev},
          "created": datetime.datetime.now().isoformat(timespec="seconds"),
          "git": f"{gi['head']}{'+dirty' if gi['dirty'] else ''}" if gi["head"] else None}
    out = Path(a.out)
    write_json(out, js)
    cols = ["module", "id", "n", "own", "reason", "L_e", "L_c"] + [f"c{k}" for k in range(8)]
    pd.DataFrame(recs).reindex(columns=cols).sort_values("id").to_csv(
        out.with_name(out.stem + "_modules.tsv"), sep="\t", index=False, float_format="%.6g")
    print(f"{slug} [{a.splits}]: seeding {nseed}, eligible (>= {a.min_clusters}) {n_elig}, own {n_own} "
          f"({100 * n_own / nseed:.1f}% of modules, {100 * n_own_ev / n_ev:.1f}% of {n_ev} events), "
          f"failed {len(failed)}" + (f" {failed}" if failed else ""))
    return 0


def cmd_fit(a):
    if a.per_module:
        return cmd_fit_hybrid(a)
    if a.glob:
        if not a.maps:
            sys.exit("fit --global needs --maps GLOB")
        return cmd_fit_global(a)
    if not a.map:
        sys.exit("fit needs --map MAP_TSV (or --global --maps GLOB)")
    return cmd_fit_map(a)


# ---- main -------------------------------------------------------------------

def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    m = sub.add_parser("map", help="position map per run type from a scan directory")
    m.add_argument("--scan-dir", required=True)
    m.add_argument("--out", required=True)
    m.add_argument("--splits", choices=("even", "odd", "all"), default="all",
                   help="harness definition: 0-based position in each run's sorted pair list")
    m.add_argument("--path", choices=("legacy", "gated"), default="legacy")
    m.add_argument("--gauss-window", choices=("2sigma", "fixed150"), default="2sigma",
                   help="window of the Gaussian bin fits: 2sigma (default) = iterative refit in "
                        "[mu-2s, mu+2s] (benchmark fit_gauss_2sigma_hist), fixed150 = the fixed "
                        "+-150*kS MeV window of the maps made before this option; recorded in the "
                        "map TSV metadata as gauss_window= (maps without it are fixed150)")
    m.add_argument("--bins", type=int, default=12)
    m.add_argument("--range", type=float, default=0.6, help="u, v in [-range, range]")
    m.add_argument("--min-n", type=int, default=300, help="bins with fewer events are flagged and not used")
    m.add_argument("--exclude-rowcol", type=int, nargs=4, default=(15, 20, 15, 20),
                   metavar=("ROW0", "ROW1", "COL0", "COL1"), help="excluded seed rectangle, 1-based inclusive")
    m.add_argument("--inner-radius", type=float, default=200.0, help="pool seed radius (mm)")
    m.add_argument("--boot", type=int, default=30,
                   help="bootstrap refits per bin for mu_err (0: sigma/sqrt(n_in_fit))")
    m.add_argument("--workers", type=int, default=1, help="processes for the bootstrap (capped at 6)")
    m.set_defaults(func=cmd_map)
    f = sub.add_parser("fit", help="poscorr_v1 JSON from a map (or the pooled global map), or a poscorr_hybrid_v1 (--per-module)")
    f.add_argument("--map", help="map_<slug>_<set>.tsv")
    f.add_argument("--form", choices=("bilinear", "poly8", "both"), default="bilinear")
    f.add_argument("--out", required=True, help="JSON path (DIR with --global)")
    f.add_argument("--global", dest="glob", action="store_true", help="pool the per-type maps")
    f.add_argument("--maps", help="with --global: glob of map_*_<set>.tsv")
    f.add_argument("--per-module", dest="per_module", action="store_true",
                   help="hybrid: per-module poly8 for seed modules with >= --min-clusters clusters, "
                        "--fallback elsewhere (poscorr_hybrid_v1 JSON + <stem>_modules.tsv)")
    f.add_argument("--scan-dir", help="with --per-module: scan directory (cluster TSVs)")
    f.add_argument("--fallback", help="with --per-module: poscorr_v1 JSON used outside the own fits")
    f.add_argument("--splits", choices=("even", "odd", "all"), default="all",
                   help="with --per-module: split set (as `map`)")
    f.add_argument("--slug", help="with --per-module: run-type slug of the runs used "
                                  "(default: the fallback's slug)")
    f.add_argument("--min-clusters", type=int, default=3000, help="with --per-module: own fit threshold")
    f.add_argument("--clamp", type=float, default=0.6, help="with --per-module: u, v clipped to +-clamp")
    f.add_argument("--path", choices=("legacy", "gated"), default="legacy")
    f.add_argument("--exclude-rowcol", type=int, nargs=4, default=(15, 20, 15, 20),
                   metavar=("ROW0", "ROW1", "COL0", "COL1"), help="excluded seed rectangle, 1-based inclusive")
    f.add_argument("--inner-radius", type=float, default=200.0, help="pool seed radius (mm)")
    f.set_defaults(func=cmd_fit)
    args = ap.parse_args(argv)
    return args.func(args)


if __name__ == "__main__":
    sys.exit(main())
