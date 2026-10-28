#!/usr/bin/env python3
"""cluster_regression.py - V0/V0b: cluster-level energy regression on top of E_R3 (plan_ml_clustering_20261027).

Run command (rebuild ~20 s, fit ~1 min per variant, oracle ~10 s on 16 cores; numpy/scipy only):

    MPLBACKEND=Agg python3 cluster_regression.py all --jobs 16
    # ablations used by report.md (each writes ridge_<features>_w<window>_*.tsv/json):
    for f in no_nmem nmem_only eabs_only ortho; do python3 cluster_regression.py fit --features $f; done
    for f in full no_nmem; do python3 cluster_regression.py fit --features $f --train-window 3; done
    python3 cluster_regression.py report

Stages (each also runs alone; outputs under output/ml_cluster_20261027/V0/):

  rebuild  per split: rebuild E_R3 from the R3 per-cluster dumps (pooled per-crystal gains, poscorr_20261025
           map, then the benchmark's own seed-indexed elastic-window step fit_elastic_window_gains), compute
           sigma_E/E with the benchmark's own fit functions (preferred fit per run type: CB for 3.49 H2 and
           2.24 C, Gauss +-2 sigma for 2.24 H2 and 0.73 H2), cache the per-cluster features, and write
           gate.tsv (rebuilt vs the R3 scan result.tsv, legacy path).  GATE: |diff| <= 0.005 pp on the 4
           reference splits of each e-p run type (best50_vs_r3_splits.tsv).
  fit      ridge (numpy, Huber-IRLS) on shape features, augment arm: y = E'/E_R3, E_corr = E_R3 * y_hat.
           Train on even split rank, evaluate on odd rank (0-based rank in the run's sorted split list), lambda
           by grouped CV inside the training set.  Shuffled-truth control: y permuted within seed module on
           the training set (E' is ~constant within a module, so permuting E' would be a no-op), refit with
           the same lambda.  Paired delta sigma/E = model - R3 per odd split, mean +- SE over splits.
  oracle   V0b: in-sample lookup median(y) per (seed module, u bin, v bin, E_seed/E_cl bin), min count; fit on
           all splits, evaluated in-sample on the odd splits.  Ceiling of the features.
  report   assemble report.md from the TSVs.

Features (no absolute energy): u, v (+ powers, products); seed module one-hot; E_seed/E_cl; 3x3 / 5x5 energy
fractions; ln n_members; second moments of the member energies; max-neighbour fraction and direction; seed ring
x (u, v) terms.  Member energies carry the per-crystal gains.  Training pool: calibrated-block clusters with
E_R3/E' in [0.85, 1.12]; evaluation pool: the untouched calibrated block.
"""
from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import math
import os
import re
import sys
import time
from multiprocessing import Pool
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

DUMPS = "output/scan_20261024/spl_percrystal/scan_{run:06d}/split_*/result_dump.npz"
EXTRA_DUMPS = ["output/percrystal_spl/dumps_extra/024386_00060_dump.npz"]
GAINS = "output/percrystal_spl/{run:06d}/gains_{run:06d}_all.tsv"
POSCORR = {"3.49 GeV H2": "output/poscorr_20261025/hybrid/3p49GeV_H2/poscorr_hybrid_3p49GeV_H2_all.json",
           "2.24 GeV H2": "output/poscorr_20261025/hybrid/2p24GeV_H2/poscorr_hybrid_2p24GeV_H2_all.json",
           "2.24 GeV C": "output/poscorr_20261025/hybrid/2p24GeV_C/poscorr_hybrid_2p24GeV_C_all.json",
           "0.73 GeV H2": "output/poscorr_20261025/0p73GeV_H2/poscorr_0p73GeV_H2_all.json"}
TYPES = ["3.49 GeV H2", "2.24 GeV H2", "0.73 GeV H2", "2.24 GeV C"]
PREF = {"3.49 GeV H2": "CB", "2.24 GeV H2": "G", "0.73 GeV H2": "G", "2.24 GeV C": "CB"}
RUNS = {"3.49 GeV H2": [24386, 24503, 24512], "0.73 GeV H2": [24590, 24603],
        "2.24 GeV H2": [25228, 25257, 25260, 25299, 25320], "2.24 GeV C": [25308]}
REF_SPLITS = {"3.49 GeV H2": (24503, [20, 40, 60, 80]), "0.73 GeV H2": (24590, [20, 40, 60, 80]),
              "2.24 GeV C": (25308, [100, 200, 300, 400]), "2.24 GeV H2": (25320, [100, 200, 300, 400])}
REF_TSV = "output/scan_20261026_best50/compare/best50_vs_r3_splits.tsv"
R3_DIR = "output/resolution_fit_20261026/r3_splits"
OUT = Path("output/ml_cluster_20261027/V0")
GATE_PP = 0.005
WIN = (0.85, 1.12)          # E_R3/E' training window (position_corr.py window)
SLUG = lambda t: t.replace(" ", "_").replace(".", "p")          # noqa: E731

_B = None


def bench():
    """The benchmark module (fit functions, fit_elastic_window_gains, poscorr_eval); read only."""
    global _B
    if _B is None:
        spec = importlib.util.spec_from_file_location("bench_v0", HERE / "benchmark_hycal_timing_test_20261021.py")
        _B = importlib.util.module_from_spec(spec)
        sys.modules["bench_v0"] = _B
        spec.loader.exec_module(_B)
    return _B


def run_type(run):
    from benchmark_hycal_timing_scan import _run_type
    return _run_type(run)


def split_sort_key(p):
    return p[1]


def list_splits():
    """{run: [(split, dump_path)]} sorted by split number, e-p run types only (no Ta)."""
    out = {}
    for t in TYPES:
        for run in RUNS[t]:
            ps = {}
            for f in sorted(Path(".").glob(DUMPS.format(run=run))):
                ps[int(re.search(r"split_(\d+)", str(f)).group(1))] = str(f)
            for f in EXTRA_DUMPS:
                m = re.match(r".*/(\d{6})_(\d{5})_dump", f)
                if m and int(m.group(1)) == run and Path(f).exists():
                    ps[int(m.group(2))] = f
            out[run] = sorted(ps.items())
    return out


def tasks_from_splits():
    """[(run, split, path, rank)] with the 0-based rank in the run's sorted split list."""
    tl = []
    for run, lst in list_splits().items():
        for rank, (sp, p) in enumerate(lst):
            tl.append((run, sp, p, rank))
    return tl


def cache_path(cache, run, sp):
    return Path(cache) / f"{run:06d}_{sp:05d}.npz"


# ---------------------------------------------------------------------------
# Estimator (the benchmark's own functions, as in main() of the benchmark)
# ---------------------------------------------------------------------------
def spectrum_sigma(E, kE, kS, bw, tref):
    """(G, CB, mu_G, mu_CB, n_in) of the post-calibration spectrum E: fit_peak -> fit_peak_cb ->
    fit_gauss_reported('2sigma') on [tref -+ 800 kE], as the benchmark reports sigma_E/E (fractions)."""
    B = bench()
    lo, hi = tref - 800.0 * kE, tref + 800.0 * kE
    nan = float("nan")
    try:
        mu, sg, A, _nf, edges, counts = B.fit_peak(E, lo, hi, bw, kS)
        cen = 0.5 * (edges[:-1] + edges[1:])
        cb = B.fit_peak_cb(E, lo, hi, bw, mu0_hint=mu, escale=kS)
        gw = B.fit_gauss_reported(cen, counts, bw, kS, lo, hi, "2sigma", start=(mu, sg))
        mG, sG = gw[0], gw[1]
        mC, sC = cb[0], cb[1]
        return (sG / mG if mG > 0 else nan, sC / mC if mC > 0 else nan, mG, mC, int(len(E)))
    except Exception as exc:        # noqa: BLE001
        sys.stderr.write(f"[WARN] spectrum_sigma failed: {exc}\n")
        return nan, nan, nan, nan, int(len(E))


# ---------------------------------------------------------------------------
# Rebuild
# ---------------------------------------------------------------------------
_PCD = {}
_GT = {}


def load_gains(run, nmod):
    if run not in _GT:
        import percrystal_calib as PCC
        gtab, g0 = PCC.read_gain_table(GAINS.format(run=run))
        garr = np.full(max(nmod, max(gtab) + 1), g0)
        for i, g in gtab.items():
            garr[i] = g
        _GT[run] = garr
    return _GT[run]


def load_poscorr(rt):
    if rt not in _PCD:
        _PCD[rt] = bench().poscorr_load(POSCORR[rt])
    return _PCD[rt]


def member_features(d, sel, garr):
    """Shape features of the selected clusters (sel: bool over the dump rows) from the member lists,
    with the per-crystal gains applied to the member energies.  Returns a dict of float arrays."""
    ptr, midx, Em = d["ptr"], d["midx"], d["E"]
    cnt = np.diff(ptr)[sel]
    starts = ptr[:-1][sel]
    nsel = len(cnt)
    ev = np.repeat(np.arange(nsel), cnt)
    off = np.cumsum(cnt) - cnt
    pos = np.arange(cnt.sum()) - np.repeat(off, cnt) + np.repeat(starts, cnt)
    mi = midx[pos]
    Eg = garr[mi] * Em[pos].astype(np.float64)
    names = d["mod_name"]
    isW = np.char.startswith(names, "W")[mi]
    mrow, mcol = d["mod_row"][mi], d["mod_col"][mi]
    srow, scol = d["row"][sel][ev], d["col"][sel][ev]
    dx = np.where(isW, mcol - scol, 99).astype(float)             # +x = +col
    dy = np.where(isW, -(mrow - srow), 99).astype(float)          # +y = -row
    isseed = mi == d["seed_idx"][sel][ev]
    tot = np.bincount(ev, Eg, nsel)
    fs = np.bincount(ev, np.where(isseed, Eg, 0.0), nsel) / tot
    in3 = isW & (np.abs(dx) <= 1) & (np.abs(dy) <= 1)
    in5 = isW & (np.abs(dx) <= 2) & (np.abs(dy) <= 2)
    f3 = np.bincount(ev, np.where(in3, Eg, 0.0), nsel) / tot
    f5 = np.bincount(ev, np.where(in5, Eg, 0.0), nsel) / tot
    dxf = np.where(isW, dx, 0.0)
    dyf = np.where(isW, dy, 0.0)
    cx = np.bincount(ev, Eg * dxf, nsel) / tot
    cy = np.bincount(ev, Eg * dyf, nsel) / tot
    Sxx = np.bincount(ev, Eg * (dxf - cx[ev]) ** 2, nsel) / tot
    Syy = np.bincount(ev, Eg * (dyf - cy[ev]) ** 2, nsel) / tot
    Sxy = np.bincount(ev, Eg * (dxf - cx[ev]) * (dyf - cy[ev]), nsel) / tot
    Ens = np.where(isseed, 0.0, Eg)
    mx = np.maximum.reduceat(Ens, off)
    hit = (Ens == mx[ev]) & (mx[ev] > 0)
    evh, ih = np.unique(ev[hit], return_index=True)
    dirx = np.zeros(nsel)
    diry = np.zeros(nsel)
    kk = np.nonzero(hit)[0][ih]
    dirx[evh] = np.clip(dxf[kk], -1, 1)
    diry[evh] = np.clip(dyf[kk], -1, 1)
    return dict(fs=fs, f3=f3, f5=f5, lnn=np.log(cnt.astype(float)), Sxx=Sxx, Syy=Syy, Sxy=Sxy,
                fmax=mx / tot, dirx=dirx, diry=diry)


def rebuild_one(task):
    run, sp, path, rank, cache = task
    B = bench()
    t0 = time.time()
    d = dict(np.load(path))
    hdr = json.loads(str(d["header"]))
    rt = run_type(run)
    n = len(d["evnum"])
    nmod = len(d["mod_index"])
    garr = load_gains(run, nmod)
    pcd = load_poscorr(rt)
    sidx = d["seed_idx"]
    names = d["mod_name"][sidx]
    sid = d["mod_id"][sidx]
    nev = np.diff(d["ptr"])
    ev = np.repeat(np.arange(n), nev)
    Eg_cl = np.bincount(ev, garr[d["midx"]] * d["E"].astype(np.float64), n)
    sx, sy = d["size_x"], d["size_y"]
    u = (d["x_cl"] - d["x_seed"]) / sx
    v = (d["y_cl"] - d["y_seed"]) / sy
    if pcd["form"] == "hybrid":
        R, _ = B.poscorr_eval(pcd, u, v, names)
    else:
        R, _ = B.poscorr_eval(pcd, u, v)
    ok = np.isfinite(u) & np.isfinite(v) & np.isfinite(R)
    Ep = np.where(ok, Eg_cl / np.where(ok, R, 1.0), Eg_cl)
    events = list(zip(sid.tolist(), d["x_seed"].tolist(), d["y_seed"].tolist(), d["row"].tolist(),
                      d["col"].tolist(), d["E_seed"].tolist(), Ep.tolist(), d["f_lms"].tolist()))
    kE, kS, bw = float(hdr["kE"]), float(hdr["kS"]), float(hdr["bin_width"])
    excl = tuple(hdr["exclude_rowcol"])
    ew_kw = dict(ebeam=float(hdr["ebeam"]), hycal_z=float(hdr["hycal_z"]), m_target=float(hdr["m_target"]),
                 escale=kS, kE=kE, inner_r=float(hdr["inner_radius"]), min_seed_E=float(hdr["min_seed_energy"]),
                 min_per_mod=int(hdr["min_per_module"]), exclude_rowcol=excl,
                 nsigma=float(hdr["window_nsigma"]), bin_width=bw, moller_nsigma=B.MOLLER_GAP_NSIGMA)
    gains, raw, corr, nmod_used, rows, sig_glob, ep_mean, mu_glob = B.fit_elastic_window_gains(events, **ew_kw)
    # pool mask (same cuts as _select_pool) and used-module mask
    r2 = d["x_seed"] ** 2 + d["y_seed"] ** 2
    pool = (r2 < ew_kw["inner_r"] ** 2) & (d["E_seed"] >= ew_kw["min_seed_E"])
    row1, col1 = d["row"] + 1, d["col"] + 1
    pool &= ~((excl[0] <= row1) & (row1 <= excl[1]) & (excl[2] <= col1) & (col1 <= excl[3]))
    gsid = np.array([gains.get(int(s), np.nan) for s in sid])
    sel = pool & np.isfinite(gsid)
    E3 = Ep * gsid
    chk = float(np.max(np.abs(E3[sel] - corr))) if sel.sum() == len(corr) and len(corr) else float("nan")
    if sel.sum() != len(corr):
        raise RuntimeError(f"{path}: block size mismatch {sel.sum()} vs {len(corr)}")
    G, CB, mG, mC, nin = spectrum_sigma(corr, kE, kS, bw, ep_mean)
    feat = member_features(d, sel, garr)
    ring = np.maximum(np.abs(d["row"][sel] - 16.5), np.abs(d["col"][sel] - 16.5)) - 0.5
    np.savez_compressed(
        cache_path(cache, run, sp),
        E=corr, eprime=d["eprime_MeV"][sel], u=u[sel], v=v[sel], seed_idx=sidx[sel].astype(np.int32),
        ring=ring.astype(np.float32), Eseed_cl=feat["fs"], f3=feat["f3"], f5=feat["f5"], lnn=feat["lnn"],
        Sxx=feat["Sxx"], Syy=feat["Syy"], Sxy=feat["Sxy"], fmax=feat["fmax"], dirx=feat["dirx"],
        diry=feat["diry"], sid=sid[sel].astype(np.int32), xs=d["x_seed"][sel], ys=d["y_seed"][sel],
        srow=d["row"][sel], scol=d["col"][sel], Eseed=d["E_seed"][sel], ewkw=np.array(json.dumps(ew_kw)),
        meta=np.array([kE, kS, bw, ep_mean, G, CB, mG, mC, nmod_used, len(corr), rank]))
    return dict(run=run, split=sp, rank=rank, type=rt, G=G, CB=CB, n_mod=nmod_used, n_ev=int(len(corr)),
                n_pool=int(pool.sum()), chk=chk, sec=time.time() - t0)


def r3_reference():
    """{(run, split): (G, CB)} in percent from the R3 scan result.tsv files (legacy), Gaussian refit with the
    2 sigma window from the stored histograms (compare_scan_dirs conventions)."""
    import compare_scan_dirs as C
    recs = C.load_scan(Path(R3_DIR), set())
    C.apply_gauss({"R3": recs}, "2sigma")
    out = {}
    for (run, sp), rec in recs.items():
        out[(run, sp)] = (C.sig(rec, "legacy", run, "G_post"), C.sig(rec, "legacy", run, "CB_post"))
    return out


def read_ref_tsv():
    """The 4 reference splits: {(run, split): (R3_G_post, R3_CB_post)} (percent, legacy)."""
    out = {}
    with open(REF_TSV) as f:
        for r in csv.DictReader(f, delimiter="\t"):
            if r["path"] != "legacy":
                continue
            out[(int(r["run"]), int(r["split"]))] = (float(r["R3_G_post"]), float(r["R3_CB_post"]))
    return out


def write_tsv(path, rows, cols):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(r[c] if isinstance(r[c], str) else f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c])
                              for c in cols) + "\n")


def cmd_rebuild(a):
    OUT.mkdir(parents=True, exist_ok=True)
    cache = Path(a.cache)
    cache.mkdir(parents=True, exist_ok=True)
    tl = [(r, s, p, k, str(cache)) for (r, s, p, k) in tasks_from_splits()]
    print(f"rebuild: {len(tl)} splits, {a.jobs} jobs", flush=True)
    t0 = time.time()
    with Pool(a.jobs) as pool:
        res = []
        for i, r in enumerate(pool.imap_unordered(rebuild_one, tl)):
            res.append(r)
            if (i + 1) % 10 == 0:
                print(f"  {i + 1}/{len(tl)}  {time.time() - t0:.0f} s", flush=True)
    res.sort(key=lambda r: (r["run"], r["split"]))
    ref_all = r3_reference()
    ref4 = read_ref_tsv()
    rows = []
    for r in res:
        key = (r["run"], r["split"])
        pref = PREF[r["type"]]
        mine = (r["CB"] if pref == "CB" else r["G"]) * 100
        rg = ref_all.get(key)
        refv = (rg[1] if pref == "CB" else rg[0]) if rg else float("nan")
        r4 = ref4.get(key)
        ref4v = (r4[1] if pref == "CB" else r4[0]) if r4 else float("nan")
        rows.append(dict(type=r["type"], run=r["run"], split=r["split"], rank=r["rank"], fit=pref,
                         rebuilt_pct=mine, r3_result_pct=refv, diff_pp=mine - refv,
                         r3_best50tsv_pct=ref4v, diff4_pp=mine - ref4v, ref_split=int(r4 is not None),
                         G_pct=r["G"] * 100, CB_pct=r["CB"] * 100, n_mod=r["n_mod"], n_ev=r["n_ev"],
                         n_pool=r["n_pool"], chk_MeV=r["chk"]))
    cols = ["type", "run", "split", "rank", "fit", "rebuilt_pct", "r3_result_pct", "diff_pp", "r3_best50tsv_pct",
            "diff4_pp", "ref_split", "G_pct", "CB_pct", "n_mod", "n_ev", "n_pool", "chk_MeV"]
    write_tsv(OUT / "gate.tsv", rows, cols)
    ref_rows = [r for r in rows if r["ref_split"]]
    worst4 = max(abs(r["diff4_pp"]) for r in ref_rows) if ref_rows else float("nan")
    fin = [abs(r["diff_pp"]) for r in rows if np.isfinite(r["diff_pp"])]
    print(f"GATE: {len(ref_rows)} reference splits, max |diff| = {worst4:.4f} pp (limit {GATE_PP}); "
          f"all {len(fin)} splits with an R3 value: max {max(fin):.4f}, mean {np.mean(fin):.4f} pp; "
          f"{time.time() - t0:.0f} s", flush=True)
    return worst4


# ---------------------------------------------------------------------------
# Features, ridge, oracle
# ---------------------------------------------------------------------------
NUMERIC = ["u", "v", "u2", "v2", "uv", "u2v2", "u4", "v4", "fs", "fs2", "f3", "f5", "lnn", "Sxx", "Syy", "Sxy",
           "fmax", "fmax_dx", "fmax_dy", "fs_u", "fs_v", "fs_u2", "fs_v2", "f3_u", "f3_v", "fmax_u_dx", "fmax_v_dy"]
RING_BINS = [1, 2, 3, 4, 5, 6]          # ring index clipped to 6; one-hot x (u, v, u^2, v^2)
ALL_NUM = NUMERIC + [f"ring{rb}_{t}" for rb in RING_BINS for t in ("u", "v", "u2", "v2")] + ["E", "E2"]


def feature_keep(feat):
    """Column indices of the numeric block for a feature set: full | no_nmem (drop ln n_members, an energy
    proxy) | nmem_only (ln n_members alone, the energy-proxy control)."""
    base = [i for i, n in enumerate(ALL_NUM) if n not in ("E", "E2")]
    if feat == "full":
        return base
    if feat in ("no_nmem", "ortho"):
        return [i for i in base if ALL_NUM[i] != "lnn"]
    if feat == "nmem_only":
        return [ALL_NUM.index("lnn")]
    if feat == "eabs_only":
        return [ALL_NUM.index("E"), ALL_NUM.index("E2")]
    raise ValueError(feat)


def load_cache(run, sp, cache):
    z = np.load(cache_path(cache, run, sp))
    return {k: z[k] for k in z.files}


def numeric_block(c):
    u, v = c["u"].astype(float), c["v"].astype(float)
    fs, f3, f5 = c["Eseed_cl"].astype(float), c["f3"].astype(float), c["f5"].astype(float)
    fm, dx, dy = c["fmax"].astype(float), c["dirx"].astype(float), c["diry"].astype(float)
    cols = [u, v, u * u, v * v, u * v, u * u * v * v, u ** 4, v ** 4, fs, fs * fs, f3, f5, c["lnn"].astype(float),
            c["Sxx"].astype(float), c["Syy"].astype(float), c["Sxy"].astype(float), fm, fm * dx, fm * dy,
            fs * u, fs * v, fs * u * u, fs * v * v, f3 * u, f3 * v, fm * u * dx, fm * v * dy]
    ring = np.clip(np.rint(c["ring"].astype(float)), 1, 6)
    for rb in RING_BINS:
        m = (ring == rb).astype(float)
        cols += [m * u, m * v, m * u * u, m * v * v]
    Ea = c["E"] / 1000.0                       # eabs_only control only: pure sculpting proxy, never in the real sets
    cols += [Ea, Ea * Ea]
    return np.column_stack(cols)


class Design:
    """Standardised numeric block + unpenalised-scale module one-hot (modules seen in training)."""

    def __init__(self, train_blocks, train_mods, keep, ortho=False):
        self.keep = keep
        self.nnum = len(keep)
        self.ortho = ortho
        Nall = np.vstack(train_blocks)
        self.eidx = [ALL_NUM.index("E"), ALL_NUM.index("E2")]
        if ortho:       # residualise the numeric features against (1, E_R3, E_R3^2) (training fit)
            self.emu, self.esd = Nall[:, self.eidx].mean(0), Nall[:, self.eidx].std(0)
            Z = np.column_stack([np.ones(len(Nall)), (Nall[:, self.eidx] - self.emu) / self.esd])
            self.coef = np.linalg.lstsq(Z, Nall[:, keep], rcond=None)[0]
        N = self._res(Nall)
        self.mu = N.mean(0)
        self.sd = N.std(0)
        self.sd[self.sd == 0] = 1.0
        self.mods = np.unique(np.concatenate(train_mods))
        self.col = {int(m): i for i, m in enumerate(self.mods)}

    def _res(self, N):
        out = N[:, self.keep]
        if self.ortho:
            Z = np.column_stack([np.ones(len(N)), (N[:, self.eidx] - self.emu) / self.esd])
            out = out - Z @ self.coef
        return out

    def X(self, N, mods):
        n = len(N)
        X = np.zeros((n, self.nnum + len(self.mods)), dtype=np.float32)
        X[:, :self.nnum] = (self._res(N) - self.mu) / self.sd
        idx = np.array([self.col.get(int(m), -1) for m in mods])
        ok = idx >= 0
        X[np.nonzero(ok)[0], self.nnum + idx[ok]] = 1.0
        return X


def stats(Xs, ts, ws):
    """Per-group sufficient statistics (X^T W X, X^T W t, sum w) in float64."""
    out = []
    for X, t, w in zip(Xs, ts, ws):
        Xw = X.astype(np.float64) * w[:, None]
        out.append((Xw.T @ X.astype(np.float64), Xw.T @ t, float(w.sum())))
    return out


def solve(stat_list, lam):
    A = sum(s[0] for s in stat_list)
    b = sum(s[1] for s in stat_list)
    n = sum(s[2] for s in stat_list)
    return np.linalg.solve(A / n + lam * np.eye(len(b)), b / n)


def huber_w(res, delta):
    a = np.abs(res)
    return np.where(a <= delta, 1.0, delta / np.maximum(a, 1e-12))


def huber_loss(res, delta):
    a = np.abs(res)
    return np.where(a <= delta, 0.5 * a * a, delta * (a - 0.5 * delta))


def fit_ridge(Xs, ts, lam, n_iter=3, lam0=1e-4):
    """Huber IRLS ridge on t = y - 1; initial weights from a plain ridge at lam0.  Returns (beta, delta, weights)."""
    ws = [np.ones(len(t)) for t in ts]
    beta = solve(stats(Xs, ts, ws), lam0)
    res = np.concatenate([t - X @ beta for X, t in zip(Xs, ts)])
    delta = 1.345 * 1.4826 * np.median(np.abs(res - np.median(res)))
    for _ in range(n_iter):
        ws = [huber_w(t - X @ beta, delta) for X, t in zip(Xs, ts)]
        beta = solve(stats(Xs, ts, ws), lam)
    return beta, delta, ws


def cv_lambda(Xs, ts, ws, delta, grid, k=5):
    """Grouped K-fold (groups = training splits) held-out Huber loss, weights fixed from the full fit."""
    st = stats(Xs, ts, ws)
    folds = [list(range(i, len(Xs), k)) for i in range(k)]
    folds = [f for f in folds if f]
    out = []
    for lam in grid:
        L = 0.0
        W = 0.0
        for f in folds:
            beta = solve([s for i, s in enumerate(st) if i not in f], lam)
            for i in f:
                r = ts[i] - Xs[i] @ beta
                L += float(huber_loss(r, delta).sum())
                W += len(r)
        out.append(L / W)
    return out


def type_data(rt, cache, tasks):
    """Per-split cached arrays for a run type, sorted, with the split rank parity."""
    items = []
    for run in RUNS[rt]:
        for (r, s, p, k) in tasks:
            if r == run:
                c = load_cache(r, s, cache)
                c["run"], c["split"], c["rank"] = r, s, k
                items.append(c)
    return items


def train_window(rt, k):
    """E_R3/E' training window: the plan's (0.85, 1.12) for k <= 0, else 1 -+ k sigma_E/E of the run type's R3
    preferred fit (mean over its splits, from gate.tsv)."""
    if k <= 0:
        return WIN
    g = [float(r["rebuilt_pct"]) / 100 for r in csv.DictReader(open(OUT / "gate.tsv"), delimiter="\t")
         if r["type"] == rt]
    s = float(np.mean(g))
    return (1 - k * s, 1 + k * s)


def tag_of(a):
    return f"{a.features}_w{'plan' if a.train_window <= 0 else format(a.train_window, 'g') + 'sig'}"


def target(c):
    return c["eprime"] / c["E"]


def eval_task(args):
    """Worker: sigma_E/E of every energy vector of one split.  A name ending in '_gm' gets the seed-indexed
    elastic-window step refitted on that energy vector (the benchmark's fit_elastic_window_gains on the cached
    calibrated-block clusters) before the spectrum fit; other names are fitted as they are."""
    run, sp, cpath, Es = args
    z = np.load(cpath)
    kE, kS, bw, tref = z["meta"][:4]
    out = {}
    kw = None
    for m, E in Es.items():
        if m.endswith("_gm"):
            if kw is None:
                kw = json.loads(str(z["ewkw"]))
                kw["exclude_rowcol"] = tuple(kw["exclude_rowcol"])
            ev = list(zip(z["sid"].tolist(), z["xs"].tolist(), z["ys"].tolist(), z["srow"].tolist(),
                          z["scol"].tolist(), z["Eseed"].tolist(), E.tolist(), [1.0] * len(E)))
            gains, _raw, corr, _n, _rows, _sg, ep, _mg = bench().fit_elastic_window_gains(ev, **kw)
            out[m] = spectrum_sigma(corr, kE, kS, bw, ep)
        else:
            out[m] = spectrum_sigma(E, kE, kS, bw, tref)
    return run, sp, out


def paired(rows, model, pref):
    """mean, SE, n of (model - R3) in pp over the split rows of one run type."""
    ix = 0 if pref == "G" else 1
    d = np.array([100.0 * (r["res"][model][ix] - r["res"]["R3"][ix]) for r in rows])
    d = d[np.isfinite(d)]
    n = len(d)
    return float(d.mean()), (float(d.std(ddof=1) / math.sqrt(n)) if n > 1 else float("nan")), n


def diag_stats(rows_E):
    """tail fraction (E/E' < 0.90), peak-region robust half-width of E/E' (percent), centre vs edge (|u| or |v| > 0.4)."""
    r = np.concatenate([x[0] / x[1] for x in rows_E])
    u = np.concatenate([x[2] for x in rows_E])
    v = np.concatenate([x[3] for x in rows_E])
    out = dict(tail=float(np.mean(r < 0.90)), med=float(np.median(r)))
    edge = (np.abs(u) > 0.4) | (np.abs(v) > 0.4)
    for lab, m in (("edge", edge), ("centre", ~edge)):
        x = r[m & (r > 0.94) & (r < 1.06)]
        q = np.quantile(x, [0.16, 0.84])
        out["w_" + lab] = float(50.0 * (q[1] - q[0]))
    return out


def cmd_fit(a):
    tasks = tasks_from_splits()
    OUT.mkdir(parents=True, exist_ok=True)
    grid = [10.0 ** e for e in (-6, -5, -4, -3, -2.5, -2, -1.5, -1, 0)]
    rng = np.random.default_rng(a.seed)
    rows_out, summary, lines = [], {}, []
    t0 = time.time()
    for rt in TYPES:
        items = type_data(rt, a.cache, tasks)
        tr = [c for c in items if c["rank"] % 2 == 0]
        te = [c for c in items if c["rank"] % 2 == 1]
        # training pool: window on E_R3/E'
        wlo, whi = train_window(rt, a.train_window)
        tr_sel = [(c, (c["E"] / c["eprime"] >= wlo) & (c["E"] / c["eprime"] <= whi)) for c in tr]
        des = Design([numeric_block(c)[s] for c, s in tr_sel], [c["seed_idx"][s] for c, s in tr_sel],
                     feature_keep(a.features), ortho=(a.features == "ortho"))
        Xs = [des.X(numeric_block(c)[s], c["seed_idx"][s]) for c, s in tr_sel]
        ts = [target(c)[s] - 1.0 for c, s in tr_sel]
        mods_tr = [c["seed_idx"][s] for c, s in tr_sel]
        beta, delta, ws = fit_ridge(Xs, ts, 1e-4)
        cvl = cv_lambda(Xs, ts, ws, delta, grid)
        lam = grid[int(np.argmin(cvl))]
        beta, delta, ws = fit_ridge(Xs, ts, lam)
        # shuffled truth: permute t within seed module over the pooled training set
        allm = np.concatenate(mods_tr)
        allt = np.concatenate(ts)
        perm = np.arange(len(allt))
        for m in np.unique(allm):
            ii = np.nonzero(allm == m)[0]
            perm[ii] = rng.permutation(ii)
        tsh = allt[perm]
        cuts = np.cumsum([0] + [len(t) for t in ts])
        ts_sh = [tsh[cuts[i]:cuts[i + 1]] for i in range(len(ts))]
        beta_sh, _d, _w = fit_ridge(Xs, ts_sh, lam)
        # evaluation on odd splits: untouched R3 pool; every model with and without the refitted g_M step
        jobs, tmp = [], []
        for c in te:
            Xe = des.X(numeric_block(c), c["seed_idx"])
            Em = {"R3": c["E"], "R3_gm": c["E"], "ridge": c["E"] * (1.0 + Xe @ beta),
                  "shuffle": c["E"] * (1.0 + Xe @ beta_sh)}
            Es = dict(Em)
            Es["ridge_gm"], Es["shuffle_gm"] = Em["ridge"], Em["shuffle"]
            jobs.append((c["run"], c["split"], cache_path(a.cache, c["run"], c["split"]), Es))
            tmp.append((c, Es))
        with Pool(a.jobs) as pool:
            res = pool.map(eval_task, jobs)
        rows = [dict(run=r, split=s_, res=x) for (r, s_, x) in res]
        pref = PREF[rt]
        ix = 0 if pref == "G" else 1
        for m in ("ridge", "shuffle", "ridge_gm", "shuffle_gm", "R3_gm"):
            summary[(rt, m)] = paired(rows, m, pref)
        for r in rows:
            cm = [c for c, _ in tmp if c["run"] == r["run"] and c["split"] == r["split"]][0]
            R = r["res"]
            row = dict(type=rt, run=r["run"], split=r["split"], rank=int(cm["rank"]), fit=pref, n_ev=R["R3"][4])
            for m in ("R3", "R3_gm", "ridge", "shuffle", "ridge_gm", "shuffle_gm"):
                row[m + "_pct"] = 100 * R[m][ix]
                if m != "R3":
                    row["d_" + m + "_pp"] = 100 * (R[m][ix] - R["R3"][ix])
            row["dG_ridge_gm_pp"] = 100 * (R["ridge_gm"][0] - R["R3"][0])
            row["dCB_ridge_gm_pp"] = 100 * (R["ridge_gm"][1] - R["R3"][1])
            rows_out.append(row)
        dg = diag_stats([(E["R3"], c["eprime"], c["u"], c["v"]) for c, E in tmp])
        dr = diag_stats([(E["ridge"], c["eprime"], c["u"], c["v"]) for c, E in tmp])
        ixm = 2 if pref == "G" else 3
        mushift = float(np.mean([100 * (r["res"]["ridge"][ixm] / r["res"]["R3"][ixm] - 1) for r in rows]))
        summary[(rt, "info")] = dict(lam=lam, cv=cvl, nfeat=des.nnum + len(des.mods), nnum=des.nnum, nmod=len(des.mods),
                                     window=[wlo, whi], ntrain=int(sum(len(t) for t in ts)), ntr_splits=len(tr),
                                     nte=len(te), delta=float(delta), diag_R3=dg, diag_ridge=dr, mushift=mushift)
        p_ = lambda k: f"{summary[(rt, k)][0]:+.4f}+-{summary[(rt, k)][1]:.4f}"      # noqa: E731
        print(f"{rt}: lambda={lam:g} nfeat={des.nnum + len(des.mods)} train={sum(len(t) for t in ts)} n={len(te)} | "
              f"ridge {p_('ridge')}  shuffle {p_('shuffle')} | +gM refit: ridge {p_('ridge_gm')}  shuffle "
              f"{p_('shuffle_gm')}  R3_gm {p_('R3_gm')}  {time.time() - t0:.0f} s", flush=True)
    cols = ["type", "run", "split", "rank", "fit", "n_ev"] + [x for m in ("R3", "R3_gm", "ridge", "shuffle", "ridge_gm", "shuffle_gm")
                                                           for x in ([m + "_pct"] + (["d_" + m + "_pp"] if m != "R3" else []))] \
        + ["dG_ridge_gm_pp", "dCB_ridge_gm_pp"]
    tag = tag_of(a)
    write_tsv(OUT / f"ridge_{tag}_splits.tsv", rows_out, cols)
    with open(OUT / f"ridge_{tag}_summary.json", "w") as f:
        json.dump({f"{k[0]}|{k[1]}": v for k, v in summary.items()}, f, indent=1, default=float)


ORC_EDGES = np.linspace(-0.5, 0.5, 6)
MIN_CELL = 30


def oracle_keys(c, fs_cuts, mod_to_i):
    iu = np.clip(np.digitize(c["u"], ORC_EDGES) - 1, 0, 4)
    iv = np.clip(np.digitize(c["v"], ORC_EDGES) - 1, 0, 4)
    ifs = np.digitize(c["Eseed_cl"], fs_cuts)
    return mod_to_i * 75 + iu * 15 + iv * 3 + ifs


def cmd_oracle(a):
    tasks = tasks_from_splits()
    summary, rows_out = {}, []
    t0 = time.time()
    for rt in TYPES:
        items = type_data(rt, a.cache, tasks)
        te = [c for c in items if c["rank"] % 2 == 1]
        allfs = np.concatenate([c["Eseed_cl"] for c in items])
        fs_cuts = np.quantile(allfs, [1 / 3, 2 / 3])
        mods = np.unique(np.concatenate([c["seed_idx"] for c in items]))
        m2i = {int(m): i for i, m in enumerate(mods)}
        key = np.concatenate([oracle_keys(c, fs_cuts, np.array([m2i[int(m)] for m in c["seed_idx"]])) for c in items])
        y = np.concatenate([target(c) for c in items])
        r = 1.0 / y
        wlo, whi = train_window(rt, a.train_window)
        w = (r >= wlo) & (r <= whi)
        kw, yw = key[w], y[w]
        order = np.argsort(kw, kind="stable")
        kw, yw = kw[order], yw[order]
        uk, st, cnt = np.unique(kw, return_index=True, return_counts=True)
        corr = {}
        for k_, s_, n_ in zip(uk, st, cnt):
            if n_ >= MIN_CELL:
                corr[int(k_)] = float(np.median(yw[s_:s_ + n_]))
        jobs = []
        cover = []
        for c in te:
            kk = oracle_keys(c, fs_cuts, np.array([m2i[int(m)] for m in c["seed_idx"]]))
            yh = np.array([corr.get(int(k_), 1.0) for k_ in kk])
            cover.append(np.mean([int(k_) in corr for k_ in kk]))
            jobs.append((c["run"], c["split"], cache_path(a.cache, c["run"], c["split"]),
                         {"R3": c["E"], "oracle": c["E"] * yh, "oracle_gm": c["E"] * yh}))
        with Pool(a.jobs) as pool:
            res = pool.map(eval_task, jobs)
        rows = [dict(run=r_, split=s, res=x) for (r_, s, x) in res]
        pref = PREF[rt]
        m, se, n = paired(rows, "oracle", pref)
        mg, seg, _n = paired(rows, "oracle_gm", pref)
        summary[rt] = dict(d=m, se=se, n=n, d_gm=mg, se_gm=seg, ncells=len(corr), cover=float(np.mean(cover)),
                           min_cell=MIN_CELL)
        ix = 0 if pref == "G" else 1
        for r_ in rows:
            rows_out.append(dict(type=rt, run=r_["run"], split=r_["split"], fit=pref,
                                 R3_pct=100 * r_["res"]["R3"][ix], oracle_pct=100 * r_["res"]["oracle"][ix],
                                 d_oracle_pp=100 * (r_["res"]["oracle"][ix] - r_["res"]["R3"][ix]),
                                 d_oracle_gm_pp=100 * (r_["res"]["oracle_gm"][ix] - r_["res"]["R3"][ix])))
        print(f"{rt}: oracle d={m:+.4f}+-{se:.4f} pp (n={n}); +gM refit {mg:+.4f}+-{seg:.4f}; cells={len(corr)} "
              f"coverage={np.mean(cover):.2f} {time.time() - t0:.0f} s", flush=True)
    tag = f"w{'plan' if a.train_window <= 0 else format(a.train_window, 'g') + 'sig'}"
    write_tsv(OUT / f"oracle_{tag}_splits.tsv", rows_out, ["type", "run", "split", "fit", "R3_pct", "oracle_pct", "d_oracle_pp", "d_oracle_gm_pp"])
    with open(OUT / f"oracle_{tag}_summary.json", "w") as f:
        json.dump(summary, f, indent=1)


# ---------------------------------------------------------------------------
# Report
# ---------------------------------------------------------------------------
VARIANTS = [("full_wplan", "full (plan features, plan window)"), ("no_nmem_wplan", "no ln n_members"),
            ("nmem_only_wplan", "ln n_members only (energy-proxy control)"),
            ("ortho_wplan", "no ln n_members, features residualised on E_R3, E_R3^2"),
            ("eabs_only_wplan", "E_R3 only (pure sculpting control)"),
            ("full_w3sig", "full, window 1 +- 3 sigma"), ("no_nmem_w3sig", "no ln n_members, window 1 +- 3 sigma")]


def cmd_report(a):
    gate = list(csv.DictReader(open(OUT / "gate.tsv"), delimiter="\t"))
    rs = {t: json.load(open(OUT / f"ridge_{t}_summary.json")) for t, _ in VARIANTS if (OUT / f"ridge_{t}_summary.json").exists()}
    orc = json.load(open(OUT / "oracle_wplan_summary.json"))
    fm = lambda v: f"{v[0]:+.3f}+-{v[1]:.3f}"          # noqa: E731
    L = ["# V0/V0b cluster-level energy regression (augment arm, ridge), 2026-10-27", ""]
    L.append("Gate: rebuilt E_R3 sigma_E/E (legacy, preferred fit) vs the R3 scan, 4 reference splits per run type, limit 0.005 pp:")
    L.append("")
    L.append("| run type | fit | R3 % | rebuilt % | max abs diff pp | all-split max / mean pp | |")
    L.append("|---|---|---|---|---|---|---|")
    for rt in TYPES:
        g = [r for r in gate if r["type"] == rt]
        r4 = [r for r in g if r["ref_split"] == "1"]
        mx = max(abs(float(r["diff4_pp"])) for r in r4)
        fin = [abs(float(r["diff_pp"])) for r in g if r["diff_pp"] != "nan"]
        L.append(f"| {rt} | {PREF[rt]} | {np.mean([float(r['r3_best50tsv_pct']) for r in r4]):.3f} | "
                 f"{np.mean([float(r['rebuilt_pct']) for r in r4]):.3f} | {mx:.4f} | {max(fin):.4f} / {np.mean(fin):.4f} | "
                 f"{'PASS' if mx <= GATE_PP else 'FAIL'} |")
    L.append("")
    L.append("Paired d sigma_E/E = model - R3 (pp, mean +- SE over odd-rank splits), plan features and window; "
             "'+gM' = seed-indexed step refitted on the model energies; shuffled = y permuted within seed module; oracle in-sample:")
    L.append("")
    L.append("| run type | n eval | ridge | ridge +gM | shuffled | shuffled +gM | oracle | oracle +gM | lambda | n feat |")
    L.append("|---|---|---|---|---|---|---|---|---|---|")
    p = rs["full_wplan"]
    for rt in TYPES:
        i, o = p[f"{rt}|info"], orc[rt]
        L.append(f"| {rt} | {p[f'{rt}|ridge'][2]} | {fm(p[f'{rt}|ridge'])} | {fm(p[f'{rt}|ridge_gm'])} | "
                 f"{fm(p[f'{rt}|shuffle'])} | {fm(p[f'{rt}|shuffle_gm'])} | {o['d']:+.3f}+-{o['se']:.3f} | "
                 f"{o['d_gm']:+.3f}+-{o['se_gm']:.3f} | {i['lam']:g} | {i['nfeat']} |")
    L.append("")
    L.append("Ablations, ridge +gM, d sigma_E/E (pp); 0.73 GeV fails its shuffled control in every variant (background-dominated pool):")
    L.append("")
    L.append("| variant | " + " | ".join(TYPES) + " |")
    L.append("|---|" + "---|" * len(TYPES))
    for t, lab in VARIANTS:
        if t in rs:
            L.append(f"| {lab} | " + " | ".join(fm(rs[t][f"{rt}|ridge_gm"]).replace("+nan+-nan", "collapsed") for rt in TYPES) + " |")
    L.append("")
    r349, o349 = p["3.49 GeV H2|ridge_gm"], orc["3.49 GeV H2"]
    o_ = rs.get("ortho_wplan", {}).get("3.49 GeV H2|ridge_gm", [float("nan")] * 3)
    L.append(f"Kill rule (3.49 GeV H2, d vs 0.02 pp): ridge {r349[0]:+.3f}, oracle {o349['d_gm']:+.3f} -> not triggered by the letter "
             "(both > 0.02 pp). But the gain is carried by ln n_members (an energy proxy) and vanishes "
             f"({o_[0]:+.3f}+-{o_[1]:.3f}) once the features are made orthogonal to E_R3: recommendation HOLD V1.")
    d3, dr = p["3.49 GeV H2|info"]["diag_R3"], p["3.49 GeV H2|info"]["diag_ridge"]
    L.append(f"Diagnostics 3.49 GeV (full, no gM): E/E'<0.90 tail {d3['tail']:.4f} -> {dr['tail']:.4f}; peak half-width centre "
             f"{d3['w_centre']:.2f} -> {dr['w_centre']:.2f} %, edge {d3['w_edge']:.2f} -> {dr['w_edge']:.2f} %.")
    L.append("Caveats: R3 gains/maps/g_M in-sample (favour R3 by 0.01-0.02 pp); shuffle permutes y within module (plan: E', a no-op); "
             "Huber-IRLS ridge on y-1; the orthogonalisation also removes legitimate leakage-fluctuation information.")
    (OUT / "report.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name in ("rebuild", "fit", "oracle", "report", "all"):
        p = sub.add_parser(name)
        p.add_argument("--jobs", type=int, default=16)
        p.add_argument("--cache", default=str(OUT / "cache"), help="per-split feature cache directory")
        p.add_argument("--model", default="ridge", choices=["ridge"])
        p.add_argument("--features", default="full", choices=["full", "no_nmem", "nmem_only", "eabs_only", "ortho"],
                       help="feature set (no_nmem drops ln n_members, an energy proxy; nmem_only is the proxy control)")
        p.add_argument("--train-window", type=float, default=0.0,
                       help="0: E_R3/E' in [0.85, 1.12] (plan); k>0: 1 -+ k sigma_E/E of the run type")
        p.add_argument("--seed", type=int, default=20261027, help="shuffle-control seed")
    a = ap.parse_args(argv)
    os.environ.setdefault("MPLBACKEND", "Agg")
    if a.cmd in ("rebuild", "all"):
        w = cmd_rebuild(a)
        if a.cmd == "all" and not (w <= GATE_PP):
            print("GATE FAILED: stopping", file=sys.stderr)
            return 2
    if a.cmd in ("fit", "all"):
        cmd_fit(a)
    if a.cmd in ("oracle", "all"):
        cmd_oracle(a)
    if a.cmd in ("report", "all"):
        cmd_report(a)
    return 0


if __name__ == "__main__":
    sys.exit(main())
