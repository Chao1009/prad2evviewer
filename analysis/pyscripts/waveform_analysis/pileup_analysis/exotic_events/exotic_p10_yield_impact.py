#!/usr/bin/env python3
"""exotic_p10_yield_impact.py -- P10 of the exotic-events impact packages: what a veto on flagged clusters costs in ep and Moller yield, and what the Ta inner-ring f1 rate follows.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p10_yield_impact.py [--jobs 8] [--out <dir>] [--cache <pkl>]

Inputs.  P7 cluster cache (output/pileup_rates/exotic_p7_cache.npz, read only; join, populations, elastic windows and Moller pairs exactly as exotic_p7_cluster_level.py:
derive / populations / moller_keys are imported, not copied), the exotic-scan npz (exotic_full: 19 Ta X17 splits; exotic_sample: H2 and the Ta 26138 splits), the benchmark scan
npz (re-clustered for the out-of-cluster f1 rows, same code path as P7 pass 1), RCDB TSV (notes/rcdb_runs_prad2_20261117.tsv) and the elog summary for run currents.
PbWO4 only, run 24937 excluded.  Clusters with E_g < 0.1 E_beam are dropped before the populations are built (the elastic and Moller populations need E_g > 0.22 E_beam; one Ta Moller pair loses a member to the cut and is not counted).

Veto definitions (per cluster; seed = largest member, categories are the exclusive f4 > f5 > f1 > f2 > f3 of P7):
  V1  seed flagged (any category)
  V2  any member flagged
  V3  saturated seed (f4 or f5) AND seed energy fraction > 0.95 (not a normal multi-module shower: the "fakes", 1-5-member clusters with a clipped seed; coordinator rule after the
      P7 saturated-seed recount, which found 5 fakes with 2-4 tiny neighbours that an n_members <= 2 rule misses)
  V3b saturated seed (f4 or f5) in a normal cluster (seed fraction <= 0.95): real showers with a clipped seed (high ADC/MeV module); NOT vetoed, listed separately
  V4  V3 + seed in f1 or f3
  (V3s = V3 with n_members == 1, shown as a sub-row.)
Loss = fraction of the population vetoed (Wilson 68 %).  Elastic window = P7's (+-3 sigma per ring); a Moller pair is lost when either member is vetoed (pair ring = inner member).

Outputs (--out): exotic_P10_yield_impact_20261118.pdf and .tsv.  Intermediate cache: output/pileup_rates/exotic_p10_cache.pkl (own file; the P7 cache is never written).
"""
from __future__ import annotations

import os
import re
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import csv
import time
import pickle
import argparse
import multiprocessing as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import exotic_p1_quantities as p1                       # noqa: E402
import exotic_p4_module_rates as p4                     # noqa: E402
import exotic_p7_cluster_level as p7                    # noqa: E402
from exotic_p1_quantities import RT_LABEL, write_tsv, EXCLUDE_RUNS     # noqa: E402
from exotic_p3_sample import text_page as _text_page    # noqa: E402
import textwrap                                         # noqa: E402


def text_page(pdf, plt, title, lines, size=9.5):
    """P3 text page with prose lines wrapped (tables keep their columns) and the font size reduced until width and height fit."""
    out = []
    for ln in lines:
        if len(ln) > 200 and not re.search(r"\S {3,}\S", ln[12:]):
            ind = len(ln) - len(ln.lstrip())
            out += textwrap.wrap(ln, 200, subsequent_indent=" " * (ind + 4)) or [""]
        else:
            out.append(ln)
    w = max((len(x) for x in out), default=1)
    sz = min(size, 1500.0 / (0.6 * w), 700.0 / (1.45 * max(len(out), 1)))
    _text_page(pdf, plt, title, out, size=sz)

STEM = "exotic_P10_yield_impact_20261118"
SCAN_ROOT = p7.SCAN_ROOT
BENCH_SCAN = p7.BENCH_SCAN
RADBINS = p4.RADBINS
EBINS = p7.EBINS
E_MIN = 0.1                                             # E_g / E_beam kept in the veto stage
VETOES = ("V1", "V2", "V3", "V3s", "V4")
VLAB = {"V1": "V1 seed flagged", "V2": "V2 any member flagged", "V3": "V3 sat. seed, non-shower", "V3s": "V3s sat. seed, 1 member", "V4": "V4 = V3 + f1/f3 seed",
        "V3b": "V3b sat. seed, normal shower (kept)"}
VCOL = {"V1": "tab:blue", "V2": "tab:cyan", "V3": "tab:red", "V3s": "tab:orange", "V4": "tab:purple", "V3b": "tab:green"}
NONSHOWER_F = 0.95
INNER_MM = RADBINS[1]
RCDB = p1.PYDIR / "notes" / "rcdb_runs_prad2_20261117.tsv"
# beam current (nA) per run: RCDB 'requested' value where present; blank in RCDB -> elog summary (flagged 'elog')
CURRENT = {24590: (23, "RCDB"), 24603: (25, "RCDB"), 24698: (400, "RCDB"), 24911: (700, "RCDB"), 25228: (320, "RCDB"), 25257: (320, "elog 06-15..23 320 nA"),
           25260: (320, "RCDB"), 25299: (320, "RCDB"), 25320: (320, "RCDB"), 25800: (3.0, "elog 07-05.. 3.0 nA (RCDB blank)"), 26138: (3.5, "RCDB"), 26230: (3.5, "RCDB")}
OWN_DT = 8.0                                            # ns: |hit time - (t0 + offset)| below this = the hit is the flagged row's own pulse
SEED_E_MEV = 30.0                                       # benchmark min_center_energy
TWIN = (100.0, 200.0)                                   # benchmark hit window (ns)


def wil(k, n):
    lo, hi = p4.wilson(k, n)
    return float(lo), float(hi)


def veto_masks(D):
    scat, nm, cnt = D["scat"].astype(np.int64), D["nmask"].astype(np.int64), D["cnt"].astype(np.int64)
    sat = (scat == 4) | (scat == 5)
    return sat, scat, nm, cnt


# ---------------------------------------------------------------------------------------------------------------- stage A: veto
def veto_stage(cache):
    c = np.load(cache, allow_pickle=True)
    metas = list(c["metas"])
    D0 = {q: c["d_" + q] for q in p7.REC + ("job",)}
    eb = np.array([m["ebeam"] for m in metas])[D0["job"]]
    sel = D0["eg"] > E_MIN * eb
    D = {q: v[sel] for q, v in D0.items()}
    del D0
    print(f"[veto] {sel.sum()} of {sel.size} clusters kept (E_g > {E_MIN} E_beam)", flush=True)
    with np.load(next((SCAN_ROOT / "exotic_sample").glob("run_*.npz"))) as z:
        names = [str(x) for x in z["mod_name"]]
    X = p7.derive(metas, D, names)
    dss = p7.datasets(metas)
    mol = p7.moller_keys()
    sat, scat, nm, cnt = veto_masks(D)
    nonsh = X["f_seed"] > NONSHOWER_F
    V = {"V1": scat > 0, "V2": (scat > 0) | (nm > 0), "V3": sat & nonsh, "V3s": sat & (cnt == 1), "V4": (sat & nonsh) | (scat == 1) | (scat == 3), "V3b": sat & ~nonsh}
    rows, extra = [], {}
    for ds in dss:
        key, lab, rt, src = ds
        P, fits, mst = p7.populations(metas, D, X, ds, mol)
        ms = [m for m in metas if m["rt"] == rt and m["src"] == src]
        n_ev = sum(m["n_events"] for m in ms)
        base = P["all"]
        el = P["elastic"]
        pops = {"hiE": base & (D["eg"] > 0.5 * X["eb"]), "elastic": el}
        for pop, pm in pops.items():
            for axn, nb in (("overall", 1), ("ring", len(RADBINS) - 1)):
                for b in range(nb):
                    m = pm if axn == "overall" else pm & (X["ring"] == b)
                    n = int(m.sum())
                    for v in VETOES + ("V3b",):
                        k = int((m & V[v]).sum())
                        lo, hi = wil(k, n)
                        rows.append(dict(table="veto", dataset=key, population=pop, axis=axn, bin="all" if axn == "overall" else f"{RADBINS[b]:g}-{RADBINS[b + 1]:g}",
                                         veto=v, n=n, k=k, frac=k / n if n else np.nan, lo=lo, hi=hi, n_events=n_ev, per_kevent=1000.0 * k / n_ev))
        # Moller pairs (Ta): either member vetoed; pair ring = inner member
        if "moller" in P:
            idx = np.flatnonzero(P["moller"])
            pk = D["job"][idx].astype(np.int64) * 10 ** 11 + D["evn"][idx].astype(np.int64)
            o = np.argsort(pk, kind="stable")
            idx, pk = idx[o], pk[o]
            first = np.r_[True, pk[1:] != pk[:-1]]
            gid = np.cumsum(first) - 1
            npair = gid[-1] + 1
            sizes = np.bincount(gid)
            ok = sizes == 2
            pring = np.full(npair, 99)
            np.minimum.at(pring, gid, X["ring"][idx])
            for v in VETOES + ("V3b",):
                pv = np.zeros(npair, bool)
                np.logical_or.at(pv, gid, V[v][idx])
                for axn, nb in (("overall", 1), ("ring", len(RADBINS) - 1)):
                    for b in range(nb):
                        mm = ok if axn == "overall" else ok & (pring == b)
                        n, k = int(mm.sum()), int((mm & pv).sum())
                        lo, hi = wil(k, n)
                        rows.append(dict(table="veto", dataset=key, population="moller_pair", axis=axn, bin="all" if axn == "overall" else f"{RADBINS[b]:g}-{RADBINS[b + 1]:g}",
                                         veto=v, n=n, k=k, frac=k / n if n else np.nan, lo=lo, hi=hi, n_events=n_ev, per_kevent=1000.0 * k / n_ev))
            extra[key + "_pairs"] = (int(npair), int(ok.sum()))
        # energy dependence: all clusters E_g > 0.1 E_beam
        eb_ = X["ebin"]
        for b in range(len(EBINS) - 1):
            m = base & (eb_ == b)
            n = int(m.sum())
            for v in VETOES + ("V3b",):
                k = int((m & V[v]).sum())
                kin = int((m & V[v] & el).sum())
                rows.append(dict(table="ebin", dataset=key, population="all", axis="ebin", bin=f"{EBINS[b]:g}-{EBINS[b + 1]:g}", veto=v, n=n, k=k, k_in_window=kin,
                                 frac=k / n if n else np.nan, n_events=n_ev, per_kevent=1000.0 * k / n_ev))
        # r histogram of single-cluster candidates (pre-window), vetoed overlay
        sg = base & X["single"]
        r = X["r"]
        bins = np.linspace(0.5, 1.5, 101)
        h = {"all": np.histogram(r[sg & np.isfinite(r)], bins)[0]}
        for v in VETOES:
            h[v] = np.histogram(r[sg & V[v] & np.isfinite(r)], bins)[0]
        extra[key + "_rhist"] = (bins, h, [(lab_, mu, s, n_, own) for lab_, mu, s, n_, own in fits])
        # saturated seeds in the elastic window: populations (a) and (b) per module
        sm = el & sat
        mods = {}
        for i in np.flatnonzero(sm):
            nmn = names[int(D["seed"][i])]
            e = mods.setdefault(nmn, [0, 0])
            e[0 if nonsh[i] else 1] += 1
        for nmn, (na, nb_) in mods.items():
            rows.append(dict(table="satmod", dataset=key, bin=nmn, n=na + nb_, k=na, k_in_window=nb_, veto="a_nonshower|b_shower"))
        rows.append(dict(table="scope", dataset=key, n=int(base.sum()), k=int(el.sum()), n_events=n_ev, bin=src))
    return rows, extra, metas


# ---------------------------------------------------------------------------------------------------------------- shared row-level helper
def rows_of(z):
    """Row-level arrays of an exotic-scan npz: PbWO4 flag, exclusive category (uniform thresholds), main pulse A and t, k, module, event index."""
    midx = z["ch_midx"].astype(np.int64)
    sel = z["mod_pwo4"][midx]
    k = z["ch_k"].astype(np.int64)
    a, t = z["ch_a"], z["ch_t"]
    main = np.argmax(np.where(np.isfinite(t), a, -1.0), axis=1)
    ar = np.arange(k.size)
    t0, a0 = t[ar, main].astype(np.float64), a[ar, main].astype(np.float64)
    with np.errstate(invalid="ignore"):
        trunc = (np.ceil(t0 / p4.DT) - p4.LEAD + p4.LW > p4.NS) & np.isfinite(t0)
    c = dict(db=z["ch_ex_db"].astype(np.float64), tail=z["ch_ex_tail"].astype(np.float64), w50=z["ch_w50"].astype(np.float64), k=k,
             nm=z["ch_ex_nm"].astype(np.int64), c2s=z["ch_ex_chi2s"].astype(np.float64), c2k=z["ch_ex_chi2k"].astype(np.float64),
             nclip=z["ch_nclip"].astype(np.int64), trunc=trunc, s4=z["ch_ex_f4"].astype(bool))
    F = p4.flags_of(c, p4.UNI)
    return dict(midx=midx, sel=sel, k=k, a0=a0, t0=t0, cat=F["cat"], ev=z["ch_ev"].astype(np.int64))


def ring_of(names):
    rad = p4.radius_of(names)
    return rad


# ---------------------------------------------------------------------------------------------------------------- stage B: occupancy
NKMAX = 200


def occ_work(job):
    src, fn = job
    z = np.load(fn)
    run, rt = int(z["run"]), str(z["run_type"])
    R = rows_of(z)
    rad = ring_of([str(x) for x in z["mod_name"]])
    inner_m = np.nan_to_num(rad[R["midx"]], nan=1e9) < INNER_MM
    evn = z["ev_number"].astype(np.int64)
    b = np.load(p7.bench_file(rt, run))
    bn, bs = b["ev_number"], b["ev_split"]
    o = np.argsort(bn)
    pos = np.searchsorted(bn[o], evn)
    ok = (pos < bn.size) & (bn[o][np.minimum(pos, bn.size - 1)] == evn)
    esplit = np.where(ok, bs[o][np.minimum(pos, bn.size - 1)], -1)
    out = []
    nev = evn.size
    sel = R["sel"]
    out_row = sel & ~inner_m
    nk_out = np.bincount(R["ev"][out_row], minlength=nev)
    n_in = np.bincount(R["ev"][sel & inner_m], minlength=nev)
    n_f1in = np.bincount(R["ev"][sel & inner_m & (R["cat"] == 1)], minlength=nev)
    n_f1 = np.bincount(R["ev"][sel & (R["cat"] == 1)], minlength=nev)
    n_all = np.bincount(R["ev"][sel], minlength=nev)
    ksum = np.bincount(R["ev"][sel], weights=R["k"][sel], minlength=nev)
    for s in sorted(set(esplit.tolist()) - {-1}):
        m = esplit == s
        d = dict(src=src, run=run, rt=rt, split=int(s), n_ev=int(m.sum()), rows=int(n_all[m].sum()), rows_in=int(n_in[m].sum()), f1_in=int(n_f1in[m].sum()), f1_all=int(n_f1[m].sum()),
                 mean_nk_out=float(nk_out[m].mean()), mean_nk_all=float(n_all[m].mean()), mean_k=float(ksum[m].sum() / max(n_all[m].sum(), 1)), mean_in=float(n_in[m].mean()))
        nk = np.minimum((n_all - n_f1in)[m], NKMAX)             # event occupancy: kept PbWO4 channels excluding the inner-ring f1 rows being counted
        d["h_ev"] = np.bincount(nk, minlength=NKMAX + 1)
        d["h_rows_in"] = np.bincount(nk, weights=n_in[m], minlength=NKMAX + 1)
        d["h_f1_in"] = np.bincount(nk, weights=n_f1in[m], minlength=NKMAX + 1)
        out.append(d)
    print(f"  [occ] {src} run {run} {rt}: {len(out)} splits", flush=True)
    return out


def occ_jobs():
    jobs, seen = [], set()
    for f in sorted((SCAN_ROOT / "exotic_full").glob("run_*_s*.npz")):
        with np.load(f) as z:
            if int(z["run"]) in EXCLUDE_RUNS:
                continue
            seen |= {(int(z["run"]), int(s)) for s in z["splits"]}
            jobs.append(("full", str(f)))
    for f in sorted((SCAN_ROOT / "exotic_sample").glob("run_*.npz")):
        with np.load(f) as z:
            if int(z["run"]) in EXCLUDE_RUNS:
                continue
            if all((int(z["run"]), int(s)) in seen for s in z["splits"]):
                continue
            jobs.append(("sample", str(f)))
    return jobs


# ---------------------------------------------------------------------------------------------------------------- stage C: out-of-cluster f1 rows (Ta full)
def ooc_work(job):
    fn, split = job
    erb = p7.erb_mod()
    t00 = time.monotonic()
    z = np.load(fn)
    run, rt = int(z["run"]), str(z["run_type"])
    R = rows_of(z)
    evn = z["ev_number"].astype(np.int64)
    key_all = evn[R["ev"]] * 4096 + R["midx"]
    flg = R["sel"] & (R["cat"] > 0)
    lut = {int(k_): (int(c), 0) for k_, c in zip(key_all[flg], R["cat"][flg])}
    fn_b = p7.bench_file(rt, run)
    b = np.load(fn_b)
    bn = b["ev_number"]
    evsel = np.flatnonzero((b["ev_split"] == split) & np.isin(bn, evn))
    gate = erb.CONTROL_GATE[rt]
    cls, geo = erb.clusterers(run, str(b["calib"]), (gate,))
    kE = round(geo["ebeam"] / erb.E_REF, 6)
    gains = p7.run_gains(rt, run, geo["name"])
    h = p7.split_hits(b, fn_b, split, evsel, geo, rt)
    ins_ev = np.isin(key_all // 4096, bn[evsel])
    flag_ev = set((key_all[flg & ins_ev] // 4096).tolist())
    _, _, used, _, _ = p7.cluster_split(h, cls[f"g{gate:g}"], geo, lut, flag_ev, gains, kE, bn)
    # hits
    hk = bn[h["h_ev"]].astype(np.int64) * 4096 + h["h_mod"].astype(np.int64)
    ho = np.argsort(hk, kind="stable")
    hks = hk[ho]
    # MeV per ADC per module from k == 1 rows with a matching hit
    sel = R["sel"] & ins_ev
    rk = key_all
    pos = np.searchsorted(hks, rk)
    posc = np.minimum(pos, hks.size - 1)
    has = sel & (hks[posc] == rk)
    hi_ = ho[posc]
    dt = h["h_t"][hi_].astype(float) - R["t0"]
    cal = has & (R["k"] == 1) & (np.abs(dt) < 8) & (R["a0"] > 50)
    ratio = h["h_e"][hi_][cal] / R["a0"][cal]
    cm = np.full(1735, np.nan)
    mm = R["midx"][cal]
    oo = np.argsort(mm)
    mm, ratio = mm[oo], ratio[oo]
    u, st = np.unique(mm, return_index=True)
    for i, m_ in enumerate(u):
        cm[m_] = np.median(ratio[st[i]:(st[i + 1] if i + 1 < len(u) else mm.size)])
    cm = np.where(np.isfinite(cm), cm, np.nanmedian(cm))
    f1 = sel & (R["cat"] == 1)
    idx = np.flatnonzero(f1)
    inc = np.array([int(rk[i]) in used for i in idx], bool)
    # neighbour search for rows outside clusters
    hev, hmod, he, ht = h["h_ev"], h["h_mod"], h["h_e"], h["h_t"].astype(float)
    pitch = 35.0
    gx, gy = np.asarray(geo["x"], float), np.asarray(geo["y"], float)
    nbE = np.zeros(idx.size)
    nbdt = np.full(idx.size, np.nan)
    ownE = np.full(idx.size, np.nan)
    ownT = np.full(idx.size, np.nan)
    bidx = {int(e): j for j, e in enumerate(bn)}
    for q, i in enumerate(idx):
        if inc[q]:
            continue
        ev_b = bidx[int(evn[R["ev"][i]])]
        lo, hi = np.searchsorted(hev, [ev_b, ev_b + 1])
        mods = hmod[lo:hi]
        m0 = int(R["midx"][i])
        mine = np.flatnonzero(mods == m0)
        if mine.size:
            j = mine[np.argmin(np.abs(ht[lo:hi][mine] - R["t0"][i]))]
            ownE[q], ownT[q] = he[lo + j], ht[lo + j]
        near = (np.hypot(gx[mods] - gx[m0], gy[mods] - gy[m0]) < pitch) & (mods != m0)
        if near.any():
            jj = np.flatnonzero(near)
            top = jj[np.argmax(he[lo:hi][jj])]
            nbE[q] = he[lo + top]
            ref = ownT[q] if np.isfinite(ownT[q]) else R["t0"][i]
            nbdt[q] = ht[lo + top] - ref
    rec = dict(run=np.full(idx.size, run), split=np.full(idx.size, split), midx=R["midx"][idx], a0=R["a0"][idx], t0=R["t0"][idx], k=R["k"][idx], inc=inc,
               eest=R["a0"][idx] * cm[R["midx"][idx]], ownE=ownE, ownT=ownT, nbE=nbE, nbdt=nbdt)
    print(f"  [ooc] run {run} split {split}: f1 rows {idx.size}, in cluster {int(inc.sum())} ({time.monotonic() - t00:.0f} s)", flush=True)
    return rec, dict(gate=gate, median_dt_cal=float(np.median(dt[cal])), n_cal=int(cal.sum()), names=list(geo["name"]), rad=ring_of(list(geo["name"])))


def ooc_jobs():
    jobs = []
    for f in sorted((SCAN_ROOT / "exotic_full").glob("run_*_s*.npz")):
        with np.load(f) as z:
            if int(z["run"]) in EXCLUDE_RUNS or str(z["run_type"]) != "2p24GeV_Ta_X17":
                continue
            jobs += [(str(f), int(s)) for s in z["splits"]]
    return jobs


def run_stages(a):
    cache = Path(a.cache)
    C = pickle.load(open(cache, "rb")) if cache.exists() else {}
    if "veto" not in C:
        rows, extra, metas = veto_stage(SCAN_ROOT / "exotic_p7_cache.npz")
        C["veto"] = (rows, extra, [{k: m[k] for k in ("src", "run", "rt", "split", "n_events", "ebeam")} for m in metas])
        pickle.dump(C, open(cache, "wb"))
    if "occ" not in C:
        with mp.get_context("fork").Pool(a.jobs) as pool:
            res = pool.map(occ_work, occ_jobs(), chunksize=1)
        C["occ"] = [d for r in res for d in r]
        pickle.dump(C, open(cache, "wb"))
    if "ooc" not in C:
        with mp.get_context("fork").Pool(a.jobs) as pool:
            res = pool.map(ooc_work, ooc_jobs(), chunksize=1)
        rec = {q: np.concatenate([r[0][q] for r in res]) for q in res[0][0]}
        C["ooc"] = (rec, [r[1] for r in res])
        pickle.dump(C, open(cache, "wb"))
    return C


# ---------------------------------------------------------------------------------------------------------------- analysis helpers
def line_fit(x, y, sig):
    """Weighted straight line y = a + b x.  Returns (a, b, sa, sb, chi2, dof); errors scaled by sqrt(chi2/dof) when chi2/dof > 1."""
    x, y, sig = (np.asarray(v, float) for v in (x, y, sig))
    w = 1.0 / sig ** 2
    S, Sx, Sy, Sxx, Sxy = w.sum(), (w * x).sum(), (w * y).sum(), (w * x * x).sum(), (w * x * y).sum()
    d = S * Sxx - Sx ** 2
    b, a = (S * Sxy - Sx * Sy) / d, (Sxx * Sy - Sx * Sxy) / d
    sb, sa = np.sqrt(S / d), np.sqrt(Sxx / d)
    chi2 = float((w * (y - a - b * x) ** 2).sum())
    dof = max(len(x) - 2, 1)
    sc = max(1.0, np.sqrt(chi2 / dof))
    return a, b, sa * sc, sb * sc, chi2, len(x) - 2


def occ_tables(occ):
    """Per-split rows, per-run rows and event-level binned rows."""
    sp = []
    for d in occ:
        r = d["f1_in"] / d["rows_in"]
        sp.append(dict(d, rate=r, err=max(np.sqrt(max(d["f1_in"], 1)) / d["rows_in"], 1e-9)))
    return sp


def ev_bins(occ, runs, edges):
    """Event-level f1 inner-ring rate vs event occupancy for the splits of the given runs: [(xmean, rate per inner row, err, n_f1, n_rows, f1 per event, err)]."""
    he = sum(d["h_ev"] for d in occ if d["run"] in runs)
    hr = sum(d["h_rows_in"] for d in occ if d["run"] in runs)
    hf = sum(d["h_f1_in"] for d in occ if d["run"] in runs)
    nk = np.arange(he.size)
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (nk >= lo) & (nk < hi)
        if he[m].sum() < 200 or hr[m].sum() < 2000:
            continue
        out.append(((he[m] * nk[m]).sum() / he[m].sum(), hf[m].sum() / hr[m].sum(), np.sqrt(max(hf[m].sum(), 1)) / hr[m].sum(), hf[m].sum(), hr[m].sum(),
                    hf[m].sum() / he[m].sum(), np.sqrt(max(hf[m].sum(), 1)) / he[m].sum()))
    return out


def fmt_fit(a, b, sa, sb, ref):
    return f"slope {b:.3g} +- {sb:.2g} per kept channel ({100 * b * ref / max(a + b * ref, 1e-12):.0f} % of the rate at the mean occupancy per +1 mean-occupancy unit x {ref:.1f})"


def pc(x, d=3):
    return f"{100 * x:.{d}f}" if np.isfinite(x) else "nan"


def get(rows, **kw):
    for r in rows:
        if all(r.get(k) == v for k, v in kw.items()):
            return r
    return None


DSN = {"0p73GeV_H2": "0.73 GeV H2", "2p24GeV_H2": "2.24 GeV H2", "3p49GeV_H2": "3.49 GeV H2", "2p24GeV_Ta_X17": "2.24 GeV Ta X17 (full)", "Ta_sample": "Ta X17 (sample)"}
MAIN = ("0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17")


def ci_s(r):
    if not r or not r["n"]:
        return "-"
    return f"{pc(r['frac'], 3)} [{pc(r['lo'], 3)},{pc(r['hi'], 3)}] ({r['k']}/{r['n']})"


def veto_tables(rows):
    L = ["Loss (% of the population vetoed, Wilson 68 % in brackets, k/n).  Elastic = P7 +-3 sigma window; Moller = selected pairs (either member vetoed); hiE = all clusters with E_g > 0.5 E_beam.", ""]
    for pop in ("elastic", "moller_pair", "hiE"):
        L.append(f"population: {pop}")
        L.append(f"{'dataset':<26}" + "".join(f"{v:<36}" for v in VETOES + ("V3b",)))
        for ds in MAIN + ("Ta_sample",):
            if get(rows, table="veto", dataset=ds, population=pop, axis="overall", veto="V1") is None:
                continue
            L.append(f"{DSN[ds]:<26}" + "".join(f"{ci_s(get(rows, table='veto', dataset=ds, population=pop, axis='overall', veto=v)):<36}" for v in VETOES + ("V3b",)))
        L.append("")
    return L


def ring_tables(rows, pop, vs=("V1", "V3", "V4")):
    L = [f"{pop}: loss % per ring (seed radius, mm), k/n"]
    L.append(f"{'dataset':<20}{'veto':<6}" + "".join(f"{f'{RADBINS[i]:g}-{RADBINS[i + 1]:g}':>22}" for i in range(len(RADBINS) - 1)))
    for ds in MAIN:
        if get(rows, table="veto", dataset=ds, population=pop, axis="ring", veto="V1") is None:
            continue
        for v in vs:
            cells = []
            for i in range(len(RADBINS) - 1):
                r = get(rows, table="veto", dataset=ds, population=pop, axis="ring", bin=f"{RADBINS[i]:g}-{RADBINS[i + 1]:g}", veto=v)
                cells.append(f"{pc(r['frac'], 2)} ({r['k']}/{r['n']})" if r and r["n"] else "-")
            L.append(f"{DSN[ds]:<20}{v:<6}" + "".join(f"{c:>22}" for c in cells))
    return L


def pages_veto(pdf, plt, rows, extra, R):
    text_page(pdf, plt, "Table 1: yield loss under each veto, per run type", veto_tables(rows), size=8.2)
    text_page(pdf, plt, "Table 2: loss per radius ring (ep elastic window, Moller pairs)", ring_tables(rows, "elastic") + [""] + ring_tables(rows, "moller_pair"), size=7.8)
    # loss vs ring
    fig, axs = plt.subplots(2, 3, figsize=(16, 9.5), sharex=True)
    for ax, (pop, ds) in zip(axs.ravel(), [("elastic", d) for d in MAIN] + [("moller_pair", "2p24GeV_Ta_X17")]):
        for j, v in enumerate(VETOES):
            for i in range(len(RADBINS) - 1):
                r = get(rows, table="veto", dataset=ds, population=pop, axis="ring", bin=f"{RADBINS[i]:g}-{RADBINS[i + 1]:g}", veto=v)
                if r and r["n"] >= 20:
                    x = i + 0.12 * (j - 2)
                    if r["k"] == 0:
                        ax.plot(x, 100 * r["hi"], "v", color=VCOL[v], ms=6)
                    else:
                        ax.errorbar(x, 100 * r["frac"], yerr=[[100 * (r["frac"] - r["lo"])], [100 * (r["hi"] - r["frac"])]], fmt="o", color=VCOL[v], ms=5, capsize=2)
            ax.plot([], [], "o", color=VCOL[v], label=VLAB[v])
        ax.set_yscale("log")
        ax.set_ylim(1e-3, 50)
        ax.set_xticks(range(len(RADBINS) - 1))
        ax.set_xticklabels([f"{RADBINS[i]:g}-{RADBINS[i + 1]:g}" for i in range(len(RADBINS) - 1)], fontsize=8)
        ax.set_title(f"{DSN[ds]}: {'ep elastic window' if pop == 'elastic' else 'Moller pairs'}", fontsize=10)
        ax.set_ylabel("fraction vetoed (%)")
        ax.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=7, loc="lower left")
    axs[1, 2].axis("off")
    for ax in axs[1]:
        ax.set_xlabel("seed ring (mm)")
    fig.suptitle("Yield loss per ring, 68 % Wilson; triangles = 0 vetoed, value = upper edge of the interval (fewer than 20 clusters in a ring are not drawn)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)
    # energy dependence
    fig, axs = plt.subplots(2, 2, figsize=(16, 9.5), sharex=True)
    eb = p7.EBINS
    xc = np.arange(len(eb) - 1)
    for ax, ds in zip(axs.ravel(), MAIN):
        for v in VETOES + ("V3b",):
            ys = [get(rows, table="ebin", dataset=ds, bin=f"{eb[i]:g}-{eb[i + 1]:g}", veto=v) for i in range(len(eb) - 1)]
            y = np.array([r["per_kevent"] if r else np.nan for r in ys])
            kin = np.array([1000.0 * r["k_in_window"] / r["n_events"] if r else np.nan for r in ys])
            ax.plot(xc, np.where(y > 0, y, np.nan), "o-", color=VCOL[v], label=VLAB[v], ms=4)
            if v in ("V3", "V4"):
                ax.plot(xc, np.where(kin > 0, kin, np.nan), "s--", color=VCOL[v], ms=4, mfc="none", label=VLAB[v] + ", in ep window")
        ax.set_yscale("log")
        ax.set_xticks(xc)
        ax.set_xticklabels([f"{eb[i]:g}-{eb[i + 1]:g}" if np.isfinite(eb[i + 1]) else f">{eb[i]:g}" for i in xc], fontsize=8)
        ax.set_title(DSN[ds], fontsize=10)
        ax.set_ylabel("vetoed clusters per 1000 events")
        ax.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=6.5)
    for ax in axs[1]:
        ax.set_xlabel("cluster energy E_g / E_beam")
    fig.suptitle("Veto-flagged clusters vs cluster energy (E_g > 0.1 E_beam); dashed = those inside the ep elastic window (the contamination a veto removes)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)
    # r histogram
    fig, axs = plt.subplots(2, 2, figsize=(16, 9.5))
    for ax, ds in zip(axs.ravel(), MAIN):
        bins, h, fits = extra[ds + "_rhist"]
        c = 0.5 * (bins[1:] + bins[:-1])
        ax.step(c, np.maximum(h["all"], 0.5), where="mid", color="k", label="single-cluster events (pre-window)")
        for v in ("V1", "V3", "V4"):
            ax.step(c, np.where(h[v] > 0, h[v], np.nan), where="mid", color=VCOL[v], label=VLAB[v])
        _, mu, s, _, _ = fits[0]
        ax.axvspan(mu - 3 * s, mu + 3 * s, color="gold", alpha=0.25, label="elastic window (all-ring fit)")
        ax.set_yscale("log")
        ax.set_title(DSN[ds], fontsize=10)
        ax.set_xlabel("E_g / E'(theta)")
        ax.set_ylabel("clusters")
        ax.grid(alpha=0.3)
    axs[0, 0].legend(fontsize=7)
    fig.suptitle("Where the vetoed single clusters sit relative to the elastic window (0.73 GeV: before the ep / Moller separation cut)", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def sat_lines(rows):
    L = ["Saturated (f4/f5) seeds inside the ep elastic window, split into (a) non-shower (seed fraction > 0.95: vetoed by V3) and (b) normal shower with a clipped seed (kept),",
         "per run type and module (seed module).  (b) recurring = high-ADC/MeV module (list for P8/P9).", ""]
    L.append(f"{'dataset':<24}{'elastic n':>10}{'(a) n':>8}{'(b) n':>8}   modules (a: count)  /  modules (b: count)")
    for ds in MAIN:
        sc = get(rows, table="scope", dataset=ds)
        ms = [r for r in rows if r["table"] == "satmod" and r["dataset"] == ds]
        a = sorted(((r["bin"], r["k"]) for r in ms if r["k"]), key=lambda t: -t[1])
        b = sorted(((r["bin"], r["k_in_window"]) for r in ms if r["k_in_window"]), key=lambda t: -t[1])
        L.append(f"{DSN[ds]:<24}{sc['k']:>10}{sum(x for _, x in a):>8}{sum(x for _, x in b):>8}   " + ", ".join(f"{m}:{n}" for m, n in a[:8]) + "  /  " + (", ".join(f"{m}:{n}" for m, n in b) or "none"))
    L += ["", "All clusters with E_g > 0.5 E_beam (V3 non-shower vs V3b shower, n):"]
    for ds in MAIN:
        r3, rb = get(rows, table="veto", dataset=ds, population="hiE", axis="overall", veto="V3"), get(rows, table="veto", dataset=ds, population="hiE", axis="overall", veto="V3b")
        L.append(f"  {DSN[ds]:<24} V3 {r3['k']:>6}   V3b {rb['k']:>6}   of {r3['n']}")
    return L


# ---------------------------------------------------------------------------------------------------------------- occupancy pages
def occ_page(pdf, plt, C):
    occ = occ_tables(C["occ"])
    ta = [d for d in occ if d["rt"] == "2p24GeV_Ta_X17"]
    h2 = [d for d in occ if d["rt"] != "2p24GeV_Ta_X17"]
    txt = []
    fig, axs = plt.subplots(2, 3, figsize=(16.5, 9.8))
    colr = {"25800": "tab:blue", "26230": "tab:red", "26138": "tab:green"}
    # (0,0) Ta split level vs mean kept channels per event
    ax = axs[0, 0]
    x, y, e = [np.array([d[q] for d in ta]) for q in ("mean_nk_all", "rate", "err")]
    for d in ta:
        ax.errorbar(d["mean_nk_all"], 100 * d["rate"], yerr=100 * d["err"], fmt="o", color=colr[str(d["run"])], ms=4)
    for r_, c_ in colr.items():
        ax.plot([], [], "o", color=c_, label=f"run {r_}")
    a, b, sa, sb, chi2, dof = line_fit(x, y, e)
    xs = np.linspace(x.min() - 0.1, x.max() + 0.1, 20)
    ax.plot(xs, 100 * (a + b * xs), "k-")
    ax.set_xlabel("mean kept PbWO4 channels per event")
    ax.set_ylabel("f1 rate, inner ring (0-120 mm) (%)")
    ax.set_title(f"Ta X17 splits: slope {100 * b:.3f} +- {100 * sb:.3f} %/channel", fontsize=10)
    ax.legend(fontsize=8)
    txt.append(("Ta split level, x = mean kept PbWO4 channels per event", a, b, sa, sb, chi2, dof, x.mean(), y.mean(), len(x)))
    # (0,1) Ta split level vs mean k
    ax = axs[0, 1]
    xk = np.array([d["mean_k"] for d in ta])
    for d in ta:
        ax.errorbar(d["mean_k"], 100 * d["rate"], yerr=100 * d["err"], fmt="o", color=colr[str(d["run"])], ms=4)
    a2, b2, sa2, sb2, chi22, dof2 = line_fit(xk, y, e)
    xs = np.linspace(xk.min() - 0.005, xk.max() + 0.005, 20)
    ax.plot(xs, 100 * (a2 + b2 * xs), "k-")
    ax.set_xlabel("mean fitted pulses per kept channel (k)")
    ax.set_ylabel("f1 rate, inner ring (%)")
    ax.set_title(f"Ta X17 splits: slope {100 * b2:.2f} +- {100 * sb2:.2f} %/unit k", fontsize=10)
    txt.append(("Ta split level, x = mean k per kept channel", a2, b2, sa2, sb2, chi22, dof2, xk.mean(), y.mean(), len(xk)))
    # (0,2) Ta event level, per run and pooled
    ax = axs[0, 2]
    edges = [0, 6, 12, 16, 20, 24, 28, 32, 36, 40, 45, 50, 60, 80, 201]
    for runs, lab, c_ in (({25800}, "25800 (3.0 nA)", "tab:blue"), ({26230, 26138}, "26230+26138 (3.5 nA)", "tab:red")):
        bn = ev_bins(C["occ"], runs, edges)
        xb, yb, eb_ = np.array([t[0] for t in bn]), np.array([t[1] for t in bn]), np.array([t[2] for t in bn])
        ax.errorbar(xb, 100 * yb, yerr=100 * eb_, fmt="o", color=c_, ms=4, label=lab)
        a3, b3, sa3, sb3, chi23, dof3 = line_fit(xb, yb, eb_)
        xs = np.linspace(xb.min(), xb.max(), 20)
        ax.plot(xs, 100 * (a3 + b3 * xs), "-", color=c_)
        txt.append((f"Ta event level {lab}, x = kept PbWO4 channels in the event (f1 rows removed)", a3, b3, sa3, sb3, chi23, dof3, float(np.average(xb, weights=[t[4] for t in bn])), yb.mean(), len(xb)))
    ax.set_xlabel("kept PbWO4 channels in the event")
    ax.set_ylabel("f1 rate, inner ring (%)")
    ax.set_title("Ta X17, event multiplicity bins", fontsize=10)
    ax.legend(fontsize=8)
    # (1,0) H2 split level
    ax = axs[1, 0]
    mk = {"0p73GeV_H2": "s", "2p24GeV_H2": "o", "3p49GeV_H2": "^"}
    cc = {"0p73GeV_H2": "tab:purple", "2p24GeV_H2": "tab:orange", "3p49GeV_H2": "tab:green"}
    for rt in mk:
        dd = [d for d in h2 if d["rt"] == rt]
        ax.errorbar([d["mean_nk_all"] for d in dd], [100 * d["rate"] for d in dd], yerr=[100 * d["err"] for d in dd], fmt=mk[rt], color=cc[rt], ms=5, label=RT_LABEL[rt])
    xh, yh, eh = [np.array([d[q] for d in h2]) for q in ("mean_nk_all", "rate", "err")]
    a4, b4, sa4, sb4, chi24, dof4 = line_fit(xh, yh, eh)
    xs = np.linspace(xh.min(), xh.max(), 20)
    ax.plot(xs, 100 * (a4 + b4 * xs), "k-")
    ax.set_xlabel("mean kept PbWO4 channels per event")
    ax.set_ylabel("f1 rate, inner ring (%)")
    ax.set_title(f"H2 splits (same scale as Ta in panel 1): slope {100 * b4:.4f} +- {100 * sb4:.4f} %/channel", fontsize=9)
    ax.legend(fontsize=8)
    txt.append(("H2 split level (all three energies), x = mean kept PbWO4 channels per event", a4, b4, sa4, sb4, chi24, dof4, xh.mean(), yh.mean(), len(xh)))
    # (1,1) H2 event level
    ax = axs[1, 1]
    for rt in mk:
        runs = {d["run"] for d in h2 if d["rt"] == rt}
        bn = ev_bins(C["occ"], runs, [0, 2, 4, 6, 8, 10, 13, 16, 20, 30, 201])
        xb, yb, eb_ = np.array([t[0] for t in bn]), np.array([t[1] for t in bn]), np.array([t[2] for t in bn])
        ax.errorbar(xb, 100 * yb, yerr=100 * eb_, fmt=mk[rt] + "-", color=cc[rt], ms=4, label=RT_LABEL[rt])
        a5, b5, sa5, sb5, chi25, dof5 = line_fit(xb, yb, eb_)
        txt.append((f"{RT_LABEL[rt]} event level, x = kept PbWO4 channels in the event", a5, b5, sa5, sb5, chi25, dof5, float(np.average(xb, weights=[t[4] for t in bn])), yb.mean(), len(xb)))
    ax.set_xlabel("kept PbWO4 channels in the event")
    ax.set_ylabel("f1 rate, inner ring (%)")
    ax.set_title("H2, event multiplicity bins", fontsize=10)
    ax.legend(fontsize=8)
    # (1,2) vs current, per run
    ax = axs[1, 2]
    runs = sorted({d["run"] for d in occ})
    cur_rows = []
    for rn in runs:
        dd = [d for d in occ if d["run"] == rn]
        f, n = sum(d["f1_in"] for d in dd), sum(d["rows_in"] for d in dd)
        lo, hi = wil(f, n)
        cur, src = CURRENT[rn]
        rt = dd[0]["rt"]
        cur_rows.append(dict(run=rn, rt=rt, current=cur, source=src, f1_in=f, rows_in=n, rate=f / n, lo=lo, hi=hi, mean_nk=float(np.mean([d["mean_nk_all"] for d in dd])), n_splits=len(dd)))
        ax.errorbar(cur, 100 * f / n, yerr=[[100 * (f / n - lo)], [100 * (hi - f / n)]], fmt="o" if rt.endswith("X17") else mk[rt], color="tab:red" if rt.endswith("X17") else cc[rt], ms=6)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("beam current (nA; RCDB request, elog where RCDB blank)")
    ax.set_ylabel("f1 rate, inner ring (%)")
    ax.set_title("per run (red: Ta X17, markers as left)", fontsize=10)
    ax.grid(alpha=0.3, which="both")
    for ax in axs.ravel()[:5]:
        ax.grid(alpha=0.3)
    fig.suptitle("P10: f1 inner-ring (0-120 mm) rate vs occupancy and beam current", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)
    # event level, f1 rows per event (not per inner row)
    fig, axs = plt.subplots(1, 2, figsize=(14, 5.6))
    for ax, groups, edges, ttl in ((axs[0], ((({25800}), "25800 (3.0 nA)", "tab:blue"), (({26230, 26138}), "26230+26138 (3.5 nA)", "tab:red")), [0, 6, 12, 16, 20, 24, 28, 32, 36, 40, 45, 50, 60, 80, 201], "Ta X17"),
                                   (axs[1], tuple(({d["run"] for d in h2 if d["rt"] == rt}, RT_LABEL[rt], cc[rt]) for rt in mk), [0, 2, 4, 6, 8, 10, 13, 16, 20, 30, 201], "H2")):
        for runs, lab, c_ in groups:
            bn = ev_bins(C["occ"], runs, edges)
            xb, yb, eb_ = (np.array([t[i] for t in bn]) for i in (0, 5, 6))
            ax.errorbar(xb, yb, yerr=eb_, fmt="o-", color=c_, ms=4, label=lab)
            a6, b6, sa6, sb6, chi26, dof6 = line_fit(xb, yb, eb_)
            txt.append((f"{lab} event level, f1 inner rows PER EVENT (not per row) vs kept channels", a6, b6, sa6, sb6, chi26, dof6, float(np.average(xb, weights=[t[4] for t in bn])), yb.mean(), len(xb)))
        ax.set_xlabel("kept PbWO4 channels in the event")
        ax.set_ylabel("inner-ring f1 rows per event")
        ax.set_title(ttl, fontsize=10)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=8)
    fig.suptitle("f1 rows per event vs in-event multiplicity: flat or falling means the f1 population is not driven by in-time activity of the same event", fontsize=11)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    pdf.savefig(fig)
    plt.close(fig)
    L = ["f1 rate (exclusive f1 / kept PbWO4 channels; a, b in fraction units, not %) in the inner ring, y = a + b x, weighted least squares (errors scaled by sqrt(chi2/dof) when > 1); n = points.",
         "Elasticity = d ln(rate) / d ln(x) at the mean = b x_mean / (a + b x_mean): 1 means rate proportional to x; intercept a / rate(x_mean) is the share of the rate that does not follow x.", ""]
    L.append(f"{'fit':<100}{'n':>4}{'a':>12}{'b +- sb':>28}{'chi2/dof':>10}{'elasticity':>12}")
    fits = []
    for lab, a, b, sa, sb, chi2, dof, xm, ym, n in txt:
        el = b * xm / max(a + b * xm, 1e-12)
        L.append(f"{lab:<100}{n:>4}{a:>12.4g}{f'{b:.4g} +- {sb:.2g}':>28}{chi2 / max(dof, 1):>10.2f}{el:>12.2f}")
        fits.append(dict(table="occ_fit", dataset=lab, n=n, a=a, b=b, sa=sa, sb=sb, chi2=chi2, dof=dof, elasticity=el))
    L += ["", f"{'run':<8}{'type':<18}{'current nA':>11}  {'source':<36}{'splits':>7}{'mean kept ch/ev':>17}{'f1 inner rows':>14}{'rate % [68 % CI]':>30}"]
    for r in cur_rows:
        L.append(f"{r['run']:<8}{RT_LABEL[r['rt']]:<18}{r['current']:>11g}  {r['source']:<36}{r['n_splits']:>7}{r['mean_nk']:>17.2f}{r['f1_in']:>14}{f'{100 * r['rate']:.4f} [{100 * r['lo']:.4f}, {100 * r['hi']:.4f}]':>30}")
    t25, t26 = [r for r in cur_rows if r["run"] == 25800][0], [r for r in cur_rows if r["run"] == 26230][0]
    L += ["", f"Ta current: RCDB gives 3.5 nA for 26138 and 26230 and is blank for 25800 (elog: 3.0 nA at 07-12).  Only two current levels exist, so no slope against current is fitted:",
          f"26230 / 25800: rate ratio {t26['rate'] / t25['rate']:.3f}, current ratio {t26['current'] / t25['current']:.3f}, mean kept channels ratio {t26['mean_nk'] / t25['mean_nk']:.3f}.",
          "H2 runs: 25257 is blank in RCDB (320 nA from the elog); 3.49 GeV 400 nA (24698) vs 700 nA (24911) is the only within-energy current lever."]
    return L, fits, cur_rows


# ---------------------------------------------------------------------------------------------------------------- out-of-cluster page
def ooc_page(pdf, plt, C):
    rec, info = C["ooc"]
    dtm = float(np.median([i["median_dt_cal"] for i in info]))
    gate = info[0]["gate"]
    inc = rec["inc"]
    o = ~inc
    tc = rec["t0"] + dtm
    anyhit = np.isfinite(rec["ownE"])
    hit = anyhit & (np.abs(rec["ownT"] - tc) < OWN_DT)          # the module's benchmark hit is this row's own pulse only if the times agree
    inw = (tc >= TWIN[0]) & (tc < TWIN[1])
    nb_hi = rec["nbE"] > np.where(hit, rec["ownE"], rec["eest"])
    n_o = int(o.sum())
    cls = [("no own-pulse hit, pulse before the window (t < 100 ns)", o & ~hit & (tc < TWIN[0])), ("no own-pulse hit, pulse after the window (t >= 200 ns)", o & ~hit & (tc >= TWIN[1])),
           ("no own-pulse hit, pulse inside the window", o & ~hit & inw), ("own-pulse hit, E < 30 MeV (below the seed threshold)", o & hit & (rec["ownE"] < SEED_E_MEV)),
           ("own-pulse hit, E >= 30 MeV, higher neighbour hit (not a local maximum)", o & hit & (rec["ownE"] >= SEED_E_MEV) & nb_hi),
           ("own-pulse hit, E >= 30 MeV, no higher neighbour (not clustered: cluster energy/size cuts)", o & hit & (rec["ownE"] >= SEED_E_MEV) & ~nb_hi)]
    L = [f"Ta X17 f1 rows (exclusive f1, PbWO4, 19 full-scan splits, 25800 + 26230): {inc.size} total, {int(inc.sum())} in a P7 cluster, {n_o} outside any cluster ({100 * n_o / inc.size:.1f} %).",
         f"NOTE: the P7 Log line says '191k of 558k outside'; the counts are the reverse, 191,223 are INSIDE a cluster and 366,283 are outside (P7 per-split n_rows_in_cluster_cat f1 / n_rows_flag_cat f1 = 0.36).",
         f"Time axis: scan pulse time t0 (largest fitted pulse) + {dtm:.1f} ns, the median offset between t0 and the benchmark hit time of k = 1 rows (n = {sum(i['n_cal'] for i in info)}); benchmark hit window {TWIN[0]:g}-{TWIN[1]:g} ns, cluster time gate {gate:g} ns, seed threshold {SEED_E_MEV:g} MeV.",
         "Energy: benchmark hit energy where a hit exists, else the fitted main amplitude x the module's median MeV/ADC (k = 1 rows matched to a hit).", "",
         f"{'class of an out-of-cluster f1 row':<92}{'n':>9}{'%':>7}{'A_main med ADC':>16}{'E med MeV':>11}"]
    ooc_rows = []
    for lab, m in cls:
        E = np.where(hit, rec["ownE"], rec["eest"])[m]
        L.append(f"{lab:<92}{int(m.sum()):>9}{100 * m.sum() / n_o:>7.1f}{np.median(rec['a0'][m]) if m.any() else np.nan:>16.1f}{np.median(E) if m.any() else np.nan:>11.1f}")
        ooc_rows.append(dict(table="ooc_class", bin=lab, n=int(m.sum()), frac=m.sum() / n_o))
    nh = o & ~hit
    L += ["", "Could they seed fake clusters?  (a) rows without an own-pulse hit (incl. rows whose module has a hit from another pulse): estimated E >= 30 MeV and no higher neighbour hit, if the benchmark trigger-window were widened by w on both sides:",
          f"{'w (ns)':>8}{'rows entering the window':>26}{'E_est >= 30 MeV':>18}{'... and no higher nb (would seed)':>36}"]
    for w in (0, 8, 16, 32, 64, 128):
        m = nh & (tc >= TWIN[0] - w) & (tc < TWIN[1] + w)
        s = m & (rec["eest"] >= SEED_E_MEV)
        s2 = s & ~(rec["nbE"] > rec["eest"])
        L.append(f"{w:>8}{int(m.sum()):>26}{int(s.sum()):>18}{int(s2.sum()):>36}")
        ooc_rows.append(dict(table="ooc_widen", bin=f"w={w}", n=int(m.sum()), k=int(s.sum()), k_in_window=int(s2.sum())))
    L += ["(b) rows with a hit of E >= 30 MeV and a higher neighbour hit: time difference to that neighbour, i.e. which cluster time gate would absorb them:", f"{'gate (ns)':>10}{'rows with |dt| < gate':>24}{'%':>8}"]
    m = o & hit & (rec["ownE"] >= SEED_E_MEV) & nb_hi
    for g in (4, 8, 16, 32, 64):
        k = int((m & (np.abs(rec["nbdt"]) < g)).sum())
        L.append(f"{g:>10}{k:>24}{100 * k / max(int(m.sum()), 1):>8.1f}")
        ooc_rows.append(dict(table="ooc_gate", bin=f"gate={g}", n=int(m.sum()), k=k))
    mi = inc & ((rec["t0"] + dtm < TWIN[0] - 5) | (rec["t0"] + dtm >= TWIN[1] + 5))
    sd = lambda w: int((nh & (tc >= TWIN[0] - w) & (tc < TWIN[1] + w) & (rec["eest"] >= SEED_E_MEV) & ~(rec["nbE"] > rec["eest"])).sum())   # noqa: E731
    w32, w64 = sd(32), sd(64)
    big = o & (np.where(hit, rec["ownE"], rec["eest"]) >= 500)
    L += ["", f"Energy reach: outside-cluster f1 rows with E (own-pulse hit or estimate) >= 500 MeV (the ep selection's SIGNAL_MIN): {int(big.sum())}; >= 100 MeV: {int((o & (np.where(hit, rec['ownE'], rec['eest']) >= 100)).sum())}.",
          f"In-cluster f1 rows: {int(mi.sum())} of {int(inc.sum())} ({100 * mi.sum() / inc.sum():.1f} %) have their largest (flagged) pulse outside the hit window by > 5 ns; the clustered hit of that module is then another, in-window pulse of the same record riding on the flagged pulse's tail.",
          "", f"Reading: the outside rows are mostly early pulses (t ~ 0, ~100 ns before the window) whose own pulse never became a hit; {int((o & anyhit & ~hit).sum())} of them share a module with a different in-window hit (another pulse).",
          f"Almost none reach the ep or Moller energy range (>= 500 MeV: see above); a trigger window wider by 32 ns would let {w32} of {n_o} rows seed 30+ MeV clusters, by 64 ns {w64}."]
    text_page(pdf, plt, "Out-of-cluster Ta f1 rows: what they are and whether they could seed clusters", L, size=7.6)
    fig, axs = plt.subplots(1, 3, figsize=(16.5, 5.6))
    ax = axs[0]
    bins = np.linspace(0, 300, 61)
    ax.hist(np.clip(rec["a0"][o], 0, 299), bins, histtype="step", density=True, color="tab:red", label=f"outside clusters ({n_o})")
    ax.hist(np.clip(rec["a0"][inc], 0, 299), bins, histtype="step", density=True, color="tab:blue", label=f"in a cluster ({int(inc.sum())})")
    ax.set_xlabel("fitted main amplitude (ADC)")
    ax.set_ylabel("density")
    ax.set_yscale("log")
    ax.legend(fontsize=8)
    ax = axs[1]
    bins = np.arange(-40, 380, 4)
    ax.hist(np.clip(tc[o], -39, 379), bins, histtype="step", color="tab:red", label="outside clusters")
    ax.hist(np.clip(tc[inc], -39, 379), bins, histtype="step", color="tab:blue", label="in a cluster")
    ax.axvspan(*TWIN, color="gold", alpha=0.25, label="benchmark hit window")
    ax.set_yscale("log")
    ax.set_xlabel("main pulse time t0 + offset (ns)")
    ax.legend(fontsize=8)
    ax = axs[2]
    E = np.where(hit, rec["ownE"], rec["eest"])
    bins = np.logspace(0, 3, 61)
    ax.hist(np.clip(E[o & hit], 1, 999), bins, histtype="step", color="tab:green", label="outside, own-pulse hit (hit E)")
    ax.hist(np.clip(rec["eest"][o & ~hit], 1, 999), bins, histtype="step", color="tab:red", label="outside, no own-pulse hit (E_est from A)")
    ax.axvline(SEED_E_MEV, color="k", ls="--", label="seed threshold 30 MeV")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("energy (MeV)")
    ax.legend(fontsize=7)
    for ax in axs:
        ax.grid(alpha=0.3)
    fig.suptitle("Ta X17 f1 rows outside any cluster: amplitude, time against the hit window, energy", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    pdf.savefig(fig)
    plt.close(fig)
    return L, ooc_rows, dict(n=inc.size, n_in=int(inc.sum()), n_out=n_o, cls=[(lab, int(m.sum())) for lab, m in cls], w0=sd(0), gate=gate, n500=int(big.sum()), mi=int(mi.sum()), medE=[float(np.median(np.where(hit, rec['ownE'], rec['eest'])[m])) if m.any() else np.nan for _, m in cls], other_hit=int((o & anyhit & ~hit).sum()), w32=w32, w64=w64)


# ---------------------------------------------------------------------------------------------------------------- summary
def sat_check(rows):
    p = HERE / "exotic_P7_saturated_seeds_20261118.tsv"
    n = {}
    if p.exists():
        for r in csv.DictReader(open(p), delimiter="\t"):
            n[r["dataset"]] = n.get(r["dataset"], 0) + 1
    mine = {ds: (get(rows, table="veto", dataset=ds, population="elastic", axis="overall", veto="V3")["k"], get(rows, table="veto", dataset=ds, population="elastic", axis="overall", veto="V3b")["k"]) for ds in ("2p24GeV_H2", "3p49GeV_H2")}
    return ("Cross-check against exotic_P7_saturated_seeds_20261118.tsv (saturated elastic seeds per dataset: " + ", ".join(f"{k} {v}" for k, v in n.items()) + "): this script V3 + V3b = "
            + ", ".join(f"{k} {a}+{b}" for k, (a, b) in mine.items()) + " (V3b = the W398 showers kept; the fakes with tiny neighbours are in V3).")


def summary_lines(rows, fits, cur_rows, ooc):
    g = lambda ds, pop, v, ax="overall", bn="all": get(rows, table="veto", dataset=ds, population=pop, axis=ax, bin=bn, veto=v)   # noqa: E731
    sc = lambda ds: get(rows, table="scope", dataset=ds)                                                                     # noqa: E731
    L = ["Loss = fraction of ep elastic-window clusters (Moller pairs) vetoed, Wilson 68 %, k/n.  V3 = saturated seed with seed energy fraction > 0.95 (the fakes); V4 = V3 + f1/f3 seeds; V2 = any member flagged.",
         "The vetoed saturated non-shower clusters are fakes by the user's decision, so for V3 the 'loss' is contamination removed, not real ep yield.",
         sat_check(rows), ""]
    for ds in MAIN:
        e = {v: g(ds, "elastic", v) for v in ("V1", "V2", "V3", "V4", "V3b")}
        n_ev = sc(ds)["n_events"]
        L.append(f"{DSN[ds]}  (elastic n = {e['V1']['n']}, events {n_ev}):")
        L.append(f"  V3 {ci_s(e['V3'])} %   V4 {ci_s(e['V4'])} %   V1 {ci_s(e['V1'])} %   V2 {ci_s(e['V2'])} %   shower with clipped seed (kept) {e['V3b']['k']}")
        if ds == "2p24GeV_Ta_X17":
            L.append(f"  Moller pairs: V3 {ci_s(g(ds, 'moller_pair', 'V3'))}   V4 {ci_s(g(ds, 'moller_pair', 'V4'))}   V2 {ci_s(g(ds, 'moller_pair', 'V2'))} %")
            r0 = g(ds, "elastic", "V4", "ring", f"{RADBINS[0]:g}-{RADBINS[1]:g}")
            L.append(f"  V4 elastic loss in the inner ring (0-120 mm) {ci_s(r0)} %; other rings <= {100 * max(get(rows, table='veto', dataset=ds, population='elastic', axis='ring', bin=f'{RADBINS[i]:g}-{RADBINS[i + 1]:g}', veto='V4')['frac'] for i in range(1, len(RADBINS) - 1)):.2f} % (Table 2)")
    L += ["", "Recommendation (veto vs resolution gain is decided at the gate after P9):",
          f"  H2 2.24 and 3.49 GeV: V3.  It removes {g('2p24GeV_H2', 'elastic', 'V3')['k']} (2.24) and {g('3p49GeV_H2', 'elastic', 'V3')['k']} (3.49) fake clusters from the elastic window at a loss of "
          f"{pc(g('2p24GeV_H2', 'elastic', 'V3')['frac'], 4)} and {pc(g('3p49GeV_H2', 'elastic', 'V3')['frac'], 4)} %, so no efficiency correction is needed (< 0.01 % : {100 * max(g('2p24GeV_H2', 'elastic', 'V3')['hi'], g('3p49GeV_H2', 'elastic', 'V3')['hi']):.3f} % upper edge).  Keep the {g('3p49GeV_H2', 'elastic', 'V3b')['k']} clipped-seed",
          "  showers (W398 at 3.49 GeV, a high-ADC/MeV module) in the sample and treat them as an energy-estimate issue (P8/P9).",
          f"  H2 0.73 GeV: no saturated seed in the window (0/{g('0p73GeV_H2', 'elastic', 'V3')['n']}); V3 is a no-op there.  V4 would cost {pc(g('0p73GeV_H2', 'elastic', 'V4')['frac'], 3)} %: not worth a veto without a P9 gain.",
          "  H2 V4 (adds f1/f3 seeds): 2.24 GeV costs "
          f"{pc(g('2p24GeV_H2', 'elastic', 'V4')['frac'], 3)} % (V3 {pc(g('2p24GeV_H2', 'elastic', 'V3')['frac'], 3)} %), 3.49 GeV {pc(g('3p49GeV_H2', 'elastic', 'V4')['frac'] - g('3p49GeV_H2', 'elastic', 'V3')['frac'], 3)} % extra; adopt only if P9 shows a resolution or peak-shift gain.",
          f"  Ta X17: V3 is a no-op (0 saturated seeds in {g('2p24GeV_Ta_X17', 'elastic', 'V3')['n']} elastic clusters).  V4 costs {pc(g('2p24GeV_Ta_X17', 'elastic', 'V4')['frac'], 3)} % of elastic and "
          f"{pc(g('2p24GeV_Ta_X17', 'moller_pair', 'V4')['frac'], 3)} % of Moller pairs, concentrated in the inner ring; if P9 shows a gain, correct the yield ring by ring with Table 2.",
          f"  V2 (any member flagged) is not viable: {pc(g('2p24GeV_Ta_X17', 'elastic', 'V2')['frac'], 1)} % (Ta elastic), {pc(g('2p24GeV_H2', 'elastic', 'V2')['frac'], 1)}-{pc(g('3p49GeV_H2', 'elastic', 'V2')['frac'], 1)} % (H2 2.24/3.49), "
          f"{pc(g('2p24GeV_Ta_X17', 'moller_pair', 'V2')['frac'], 1)} % of Ta Moller pairs.", ""]
    k1 = next(f for f in fits if f["dataset"].startswith("Ta split level, x = mean kept PbWO4"))
    kev = [f for f in fits if f["dataset"].startswith("Ta event level") and "PER EVENT" not in f["dataset"]]
    t25, t26 = [r for r in cur_rows if r["run"] == 25800][0], [r for r in cur_rows if r["run"] == 26230][0]
    h700, h400 = [r for r in cur_rows if r["run"] == 24911][0], [r for r in cur_rows if r["run"] == 24698][0]
    h73b, h73a = [r for r in cur_rows if r["run"] == 24603][0], [r for r in cur_rows if r["run"] == 24590][0]
    kpe = [f for f in fits if f["dataset"].startswith(("25800", "26230")) and "PER EVENT" in f["dataset"]]
    t38 = [r for r in cur_rows if r["run"] == 26138][0]
    L += ["Ta inner-ring f1 rate vs occupancy and current:",
          f"  Across runs the f1 rate follows the occupancy: slope {100 * k1['b']:.4f} +- {100 * k1['sb']:.4f} % per kept channel per event ({k1['n']} splits, elasticity {k1['elasticity']:.2f}, chi2/dof {k1['chi2'] / max(k1['dof'], 1):.0f}: run-level scatter above the binomial error, only 3 runs).",
          f"  26230 / 25800: rate x{t26['rate'] / t25['rate']:.2f}, current x{t26['current'] / t25['current']:.2f} (3.5 / 3.0 nA; 25800 from the elog, RCDB blank), kept channels x{t26['mean_nk'] / t25['mean_nk']:.2f}.  "
          f"26138 has the same 3.5 nA as 26230 but fewer kept channels ({t38['mean_nk']:.1f} vs {t26['mean_nk']:.1f}) and a lower rate ({100 * t38['rate']:.2f} vs {100 * t26['rate']:.2f} %): current and occupancy cannot be separated with three runs.",
          f"  Within a run, event multiplicity is only a partial proxy: f1 rows per event rise sub-proportionally (elasticity {kpe[0]['elasticity']:.2f} for 25800, {kpe[1]['elasticity']:.2f} for 26230+26138) and the rate per inner row falls "
          f"(elasticity {kev[0]['elasticity']:.2f}, {kev[1]['elasticity']:.2f}).",
          "  So f1 tracks run-level occupancy (elasticity ~1, 3 runs); current, occupancy, trigger version (25800: X17_v3; 26138/26230: X17_v4_RawSum950) and the 08-10 inner-layer HV change are confounded and cannot be separated here.",
          "  Consistent with descending tails of earlier pulses (their rate grows with the pulse rate), not proof: the within-run dependence on in-event multiplicity is only partial.",
          f"  H2 (event level, three energies): per-row elasticity {min(f['elasticity'] for f in fits if 'GeV H2 event level, x' in f['dataset']):.1f}-{max(f['elasticity'] for f in fits if 'GeV H2 event level, x' in f['dataset']):.1f} (chi2/dof {min(f['chi2'] / max(f['dof'], 1) for f in fits if 'GeV H2 event level, x' in f['dataset']):.0f}-{max(f['chi2'] / max(f['dof'], 1) for f in fits if 'GeV H2 event level, x' in f['dataset']):.0f}: not a straight line) - the H2 f1 rows rise faster than the multiplicity, a different mix from Ta.",
          f"  H2: rates {100 * min(r['rate'] for r in cur_rows if not r['rt'].endswith('X17')):.3f}-{100 * max(r['rate'] for r in cur_rows if not r['rt'].endswith('X17')):.3f} % (Ta {100 * t25['rate']:.2f}-{100 * t26['rate']:.2f} %); 3.49 GeV 700 vs 400 nA: rate x{h700['rate'] / h400['rate']:.2f} for current x1.75 "
          f"(single runs, statistical errors only); 0.73 GeV 24603 / 24590: x{h73b['rate'] / h73a['rate']:.2f} for current x{h73b['current'] / h73a['current']:.2f}, so another run condition (24603 has the 250 MeV threshold) matters there.", ""]
    L += ["Out-of-cluster Ta f1 rows:",
          f"  {ooc['n_out']} of {ooc['n']} f1 rows ({100 * ooc['n_out'] / ooc['n']:.0f} %) are outside any cluster (the P7 Log line has the counts reversed: 191k are inside).  Classes: "
          + "; ".join(f"{100 * n_ / ooc['n_out']:.1f} % {lab}, E med {me:.0f} MeV" for (lab, n_), me in zip(ooc["cls"], ooc["medE"]) if n_ > 0.001 * ooc["n_out"]) + ".",
          f"  {ooc['other_hit']} of them share a module with a different in-window hit that is a different pulse (not counted as their own hit).",
          f"  Practically none can enter the ep or Moller windows (E >= 500 MeV: {ooc['n500']} rows, {100 * ooc['n500'] / ooc['n_out']:.2f} %).  {ooc['w0']} in-window rows have E_est >= 30 MeV and no higher neighbour but did not become hits (the module's hit is often another pulse); a wider trigger window would, see the table (a trigger window wider by 32 / 64 ns: {ooc['w32']} / {ooc['w64']} rows would seed 30+ MeV clusters).  No action for the yields.",
          f"  {ooc['mi']} of {ooc['n_in']} in-cluster f1 rows are flagged on a pulse outside the hit window (mostly early) while the clustered hit is another pulse of the record: relevant for how P8 assigns the f1 correction."]
    return L


def report(C, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    out = Path(out)
    rows, extra, metas = C["veto"]
    with PdfPages(out / f"{STEM}.pdf") as pdf:
        head = [__doc__.split("Inputs.")[1].split("Outputs")[0].strip()]
        text_page(pdf, plt, "P10 impact on yields: definitions and data", "\n".join(head).split("\n"), size=8.4)
        pages_veto(pdf, plt, rows, extra, None)
        text_page(pdf, plt, "Saturated seeds in the elastic window: non-shower fakes vs real showers with a clipped seed", sat_lines(rows), size=8.2)
        Lo, fits, cur_rows = occ_page(pdf, plt, C)
        text_page(pdf, plt, "f1 inner-ring rate vs occupancy and beam current: fits and run table", Lo, size=7.6)
        Lc, ooc_rows, ooc_info = ooc_page(pdf, plt, C)
        S = summary_lines(rows, fits, cur_rows, ooc_info)
        text_page(pdf, plt, "Summary: recommended veto per run type (input to the gate after P9 / P10)", S, size=8.2)
    allrows = [dict(r) for r in rows if r["table"] in ("veto", "ebin", "satmod", "scope")]
    for f in fits:
        allrows.append(f)
    for r in cur_rows:
        allrows.append(dict(table="occ_run", dataset=f"{r['run']}", bin=r["rt"], n=r["rows_in"], k=r["f1_in"], frac=r["rate"], lo=r["lo"], hi=r["hi"], a=r["current"], b=r["mean_nk"]))
    for d in occ_tables(C["occ"]):
        allrows.append(dict(table="occ_split", dataset=f"{d['run']}", bin=f"{d['rt']}_s{d['split']}", n=d["rows_in"], k=d["f1_in"], frac=d["rate"], n_events=d["n_ev"], a=d["mean_nk_all"], b=d["mean_k"]))
    allrows += ooc_rows
    cols = ("table", "dataset", "population", "axis", "bin", "veto", "n", "k", "k_in_window", "frac", "lo", "hi", "n_events", "per_kevent", "a", "b", "sa", "sb", "chi2", "dof", "elasticity")
    write_tsv(out / f"{STEM}.tsv", allrows, cols)
    print("\n".join(S))
    print("\n".join(Lo))
    print("\n".join(Lc))


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--cache", default=str(SCAN_ROOT / "exotic_p10_cache.pkl"))
    ap.add_argument("--jobs", type=int, default=8)
    a = ap.parse_args()
    report(run_stages(a), a.out)
