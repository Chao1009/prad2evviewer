#!/usr/bin/env python3
"""exotic_p13_recluster.py -- P13 of the exotic-events impact packages: the P7 / P8 / P10 (ep) cluster-level studies recomputed under a channel-level DROP policy.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p13_recluster.py [--jobs 8] [--cache <npz>] [--tsv <path>] [--pdf <path>] [--check-p7]

Policy (user decision).  Any kept PbWO4 channel record of exclusive category f1..f5 (uniform P2 thresholds, exotic_p7_cluster_level.row_categories, f4 > f5 > f1 > f2 > f3)
is not used in the clustering or in the energy calibration: every benchmark hit of that (run, EVIO event number, module) is removed from the hit list BEFORE clustering.
No cluster veto, no correction.

Inputs.  Exactly P7's join (exotic_p7_cluster_level.work: exotic-scan split x benchmark-scan split, same hits, phase factor, time offsets, gated HyCalCluster at the run type's
CONTROL_GATE, control-v3 gains), imported helpers row_categories / bench_file / split_hits / run_gains / datasets / derive / populations / moller_keys / fit_window.
Data scope as P7 minus the Ta 'sample' cross-check: H2 run types from exotic_sample, Ta X17 from the 19 exotic_full npz.  PbWO4 only, run 24937 excluded.

Procedure (pass 1, per split).  The same benchmark hits are clustered twice: A = default (reproduces the P7 cache record for record), D = the hits minus every hit whose
(event number, module) has a flagged row.  Only events that lose at least one hit are reclustered for D ("changed" events); in every other event D = A.  D records are therefore
stored for changed events only (REC format, key 'job' as A); an A cluster in an unchanged event has fate 0, eg_D = eg_A, dmatch = -1.  Per A cluster (matched within the event):
fate 0 = survives with the same seed module; 1 = re-seeded (its seed hit was dropped; the largest not yet claimed D cluster seeded on one of its former member modules); 2 = vanished
(seed dropped, no D cluster seeded on a former member); 4 = merged into a neighbour (seed dropped, the only D clusters seeded on former members are matched to other A clusters, i.e.
its members ended up in a surviving neighbour); 3 = seed not dropped but no same-seed D cluster (counted; eg_D = the largest unclaimed D cluster seeded on a member, NaN if none).
Each D cluster is matched to at most one A cluster (fate 0 first, then fate 1 largest A cluster first, then fate 3).  S2 = (sum of eg_D over the distinct matched D clusters
- sum E_A) / sum E_A over the sample's A clusters: vanished, merged (their energy is inside the neighbour's D cluster, counted there) and unmatched fate 3 contribute 0 to the D sum.  Also eg_D, the removed
(flagged, gain-corrected) member energy, number of dropped members, centre dropped, the index of the matched D record, and one record per flagged member (cluster, category,
role seed = largest member / neighbour, gain-corrected energy).

Pass 2 (tables).  Populations (P7's populations(): all, hiE, elastic = exactly one cluster above SIGNAL_MIN, seed >= SEED_MIN, +-3 sigma per ring with the 0.73 GeV Moller veto /
separation, moller = Moller-pair clusters) are built on A.  P7 section: affected, seed lost (vanish / re-seed), removed energy fraction, fakes (V3: f4/f5 seed with seed fraction
> 0.95) removed.  P8 section (first order): flagged members dropped per category x role x ring (their full default energy), cluster dE/E of affected clusters, sample energy-scale
shift S1 = -sum E_flagged / sum E_cluster (ratio-estimator error), second order S2 = (sum E_D - sum E_A) / sum E_A, fraction of clusters whose seed is dropped.  P10 section (ep):
A elastic-window clusters lost under D (vanish; the event no longer passes the single-cluster selection; re-formed cluster outside A's window with A's per-ring mu / sigma,
no refit) and clusters entering the window (D elastic, not the match of an A elastic cluster), with the P10 V3 / V4 veto losses beside them
(exotic_P10_yield_impact_20261118.tsv).  The Moller-pair part of P10 is not done here.

Outputs.  Cache output/pileup_rates/exotic_p13_cache.npz (own file; the P7 / P8 / P10 caches are never written).  Standalone: --tsv long table, --pdf pages.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import csv
import time
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
from exotic_p1_quantities import RT_LABEL, write_tsv     # noqa: E402
from exotic_p3_sample import text_page                  # noqa: E402

STEM = "exotic_P13_drop_flagged_20261118"
CATS = p7.CATS
REC = p7.REC
RADBINS = p4.RADBINS
CACHE = p7.SCAN_ROOT / "exotic_p13_cache_v2.npz"
P10_TSV = HERE / "exotic_P10_yield_impact_20261118.tsv"
NONSHOWER_F = 0.95                                      # P10 V3: seed energy fraction above this = not a normal shower
FATE_LAB = {0: "same seed", 1: "re-seeded", 2: "vanished", 3: "seed kept, other D cluster / none", 4: "merged into a neighbour"}
SAMPLES = ("all", "hiE", "elastic", "moller")
FLOAT_Q = ("seed_e", "e", "eg", "x", "y", "es", "en")
INT16_Q = ("seed", "nbig", "cnt", "scat", "sraw", "nmask")
FATE_Q = ("fate", "egD", "erem", "ndrop", "cdrop", "dmatch", "chg")
MEM_Q = ("m_cl", "m_cat", "m_role", "m_e")


# ---------------------------------------------------------------------------------------------------------------- pass 1: per split
def clusters_ev(cl, idx, mod, e, t, pw, sig, gains):
    """PbWO4-seeded clusters of one event from the hits idx (same calls as p7.cluster_split): list of dicts (member energies em = {module: gain-corrected E}), and the LG-seeded count."""
    cl.clear()
    for i in idx:
        cl.add_hit(mod[i], e[i], t[i])
    cl.form_clusters()
    res = cl.reconstruct_matched()
    nbig = sum(r.hit.energy > sig for r in res)
    out, nlg = [], 0
    for r in res:
        seed = r.cluster.center.index
        if not pw[seed]:
            nlg += 1
            continue
        em = {}
        for x in r.cluster.hits:
            em[x.index] = em.get(x.index, 0.0) + x.energy * gains[x.index]
        out.append(dict(seed=seed, seed_e=r.cluster.center.energy, e=r.hit.energy, x=r.hit.x, y=r.hit.y, nbig=nbig, em=em, ndup=len(r.cluster.hits) - len(em)))
    return out, nlg


def rec_of(c, en_, lut, fl):
    """P7 REC tuple of a cluster dict c (same rules as p7.cluster_split), the flagged members [(module, category, is_seed role, energy)], the number of seed != largest member."""
    mi, me = list(c["em"]), np.fromiter(c["em"].values(), float, len(c["em"]))
    eg = float(me.sum())
    smax = mi[int(me.argmax())]
    scat = sraw = nmask = 0
    es = 0.0
    enb = [0.0] * 6
    fm = []
    if fl:
        for m, x in zip(mi, me):
            v = lut.get(en_ * 4096 + m)
            if v is None:
                continue
            cc, rw = v
            fm.append((m, cc, int(m == smax), float(x)))
            if m == smax:
                scat, sraw, es = cc, rw, float(x)
            else:
                nmask |= 1 << (cc - 1)
                enb[cc] += float(x)
    return (en_, 0, c["seed"], c["seed_e"], c["e"], eg, c["x"], c["y"], c["nbig"], len(mi), scat, sraw, nmask, es, enb[1:]), fm, int(smax != c["seed"])


def fates_ev(ca, fms, Dev):
    """[(fate, eg_D, matched D index)] of the A clusters ca (flagged members fms) of one event against its D clusters Dev; every D cluster is matched to at most one A cluster.
    1st pass: A clusters whose seed hit survived take the unclaimed D cluster with the same seed (nearest energy): fate 0.  2nd pass: A clusters whose seed was dropped, largest first,
    take the largest unclaimed D cluster seeded on one of their member modules: fate 1 (re-seeded); none unclaimed but such a D cluster exists (it belongs to a neighbour: the
    members were merged into it): fate 4; no such D cluster: fate 2 (vanished).  3rd pass: seed kept but no same-seed D cluster: the largest unclaimed D cluster seeded on a member
    (eg_D, matched) or none (eg_D = NaN), fate 3."""
    eg = [sum(a["em"].values()) for a in ca]
    egd = [sum(d["em"].values()) for d in Dev]
    drop = [{m for m, _, _, _ in fm} for fm in fms]
    out = [None] * len(ca)
    claimed = set()
    for i, a in enumerate(ca):
        if a["seed"] in drop[i]:
            continue
        cand = [j for j, d in enumerate(Dev) if d["seed"] == a["seed"] and j not in claimed]
        if cand:
            j = min(cand, key=lambda j: abs(egd[j] - eg[i]))
            claimed.add(j)
            out[i] = (0, float(egd[j]), j)
    for i in sorted((i for i in range(len(ca)) if out[i] is None and ca[i]["seed"] in drop[i]), key=lambda i: -eg[i]):
        allc = [j for j, d in enumerate(Dev) if d["seed"] in ca[i]["em"] and d["seed"] not in drop[i]]
        cand = [j for j in allc if j not in claimed]
        if cand:
            j = max(cand, key=lambda j: egd[j])
            claimed.add(j)
            out[i] = (1, float(egd[j]), j)
        else:
            out[i] = (4 if allc else 2, np.nan, -1)
    for i in range(len(ca)):
        if out[i] is None:
            cand = [j for j, d in enumerate(Dev) if d["seed"] in ca[i]["em"] and d["seed"] not in drop[i] and j not in claimed]
            if cand:
                j = max(cand, key=lambda j: egd[j])
                claimed.add(j)
                out[i] = (3, float(egd[j]), j)
            else:
                out[i] = (3, np.nan, -1)
    return out


def work(job):
    """One exotic-scan split joined to its benchmark-scan split, clustered with all hits (A) and with the flagged-row hits dropped (D)."""
    src, fn, split = job
    erb = p7.erb_mod()
    t0 = time.monotonic()
    z = np.load(fn)
    run, rt = int(z["run"]), str(z["run_type"])
    key, cat, raw, n_rows = p7.row_categories(z)
    fn_b = p7.bench_file(rt, run)
    b = np.load(fn_b)
    bn = b["ev_number"]
    evsel = np.flatnonzero((b["ev_split"] == split) & np.isin(bn, z["ev_number"]))
    gate = erb.CONTROL_GATE[rt]
    cls, geo = erb.clusterers(run, str(b["calib"]), (gate,))
    cl = cls[f"g{gate:g}"]
    kE = round(geo["ebeam"] / erb.E_REF, 6)
    gains = p7.run_gains(rt, run, geo["name"])
    h = p7.split_hits(b, fn_b, split, evsel, geo, rt)
    ins = np.isin(key // 4096, bn[evsel])
    lut = {int(k_): (int(c), int(r)) for k_, c, r in zip(key[ins], cat[ins], raw[ins])}
    flag_ev = set((key[ins] // 4096).tolist())
    ev = h["h_ev"]
    A = {q: [] for q in REC}
    D = {q: [] for q in REC}
    F = {q: [] for q in FATE_Q}
    M = {q: [] for q in MEM_Q}
    nlg = nmax = ndup = n_chg = n_hdrop = 0
    n_fate = np.zeros(5, np.int64)
    if ev.size:
        mod, e = h["h_mod"].tolist(), h["h_e"].tolist()
        t = (h["h_t"].astype(np.float64) - geo["toff"][h["h_mod"]]).tolist() if geo.get("toff") is not None else h["h_t"].tolist()
        keyh = bn[ev].astype(np.int64) * 4096 + h["h_mod"].astype(np.int64)
        dropped = np.isin(keyh, np.fromiter(lut, np.int64, len(lut))) if lut else np.zeros(ev.size, bool)
        lo = np.flatnonzero(np.r_[True, ev[1:] != ev[:-1]])
        hi = np.r_[lo[1:], ev.size]
        pw, sig = geo["pwo4"], erb.SIGNAL_MIN * kE
        for s, f in zip(lo.tolist(), hi.tolist()):
            en_ = int(bn[ev[s]])
            fl = en_ in flag_ev
            ca, nl = clusters_ev(cl, range(s, f), mod, e, t, pw, sig, gains)
            nlg += nl
            chg = bool(dropped[s:f].any())
            if chg:
                n_chg += 1
                n_hdrop += int(dropped[s:f].sum())
                keep = [i for i in range(s, f) if not dropped[i]]
                cd, _ = clusters_ev(cl, keep, mod, e, t, pw, sig, gains) if keep else ([], 0)
                d0 = len(D["e"])
                for c in cd:
                    r_, _, _ = rec_of(c, en_, lut, False)
                    for q, v in zip(REC, r_):
                        D[q].append(v)
            else:
                cd = []
            recs = [rec_of(c, en_, lut, fl) for c in ca]
            fts = fates_ev(ca, [x[1] for x in recs], cd) if chg else None
            for k_, (c, (r_, fm, nm)) in enumerate(zip(ca, recs)):
                nmax += nm
                ndup += c["ndup"]
                ia = len(A["e"])
                for q, v in zip(REC, r_):
                    A[q].append(v)
                for m_, cc, ro, x in fm:
                    M["m_cl"].append(ia)
                    M["m_cat"].append(cc)
                    M["m_role"].append(ro)
                    M["m_e"].append(x)
                if chg:
                    fa, egd, j = fts[k_]
                    dm = d0 + j if j >= 0 else -1
                else:
                    fa, egd, dm = 0, r_[5], -1
                n_fate[fa] += 1
                for q, v in zip(FATE_Q, (fa, egd, sum(x for *_, x in fm), len(fm), c["seed"] in {m_ for m_, *_ in fm}, dm, chg)):
                    F[q].append(v)
    n = len(A["e"])
    nd = len(D["e"])

    def pack(R, nn):
        d = {q: np.array(v, dtype=np.float32 if q in FLOAT_Q else None) for q, v in R.items()}
        d["en"] = d["en"].reshape(nn, 5) if nn else np.zeros((0, 5), np.float32)
        d["split"] = np.full(nn, split, np.int32)
        for q in INT16_Q:
            d[q] = d[q].astype(np.int16)
        return d

    dA, dD = pack(A, n), pack(D, nd)
    f_ = {"fate": np.array(F["fate"], np.int8), "egD": np.array(F["egD"], np.float32), "erem": np.array(F["erem"], np.float32), "ndrop": np.array(F["ndrop"], np.int16),
          "cdrop": np.array(F["cdrop"], bool), "dmatch": np.array(F["dmatch"], np.int64), "chg": np.array(F["chg"], bool)}
    m_ = {"m_cl": np.array(M["m_cl"], np.int64), "m_cat": np.array(M["m_cat"], np.int8), "m_role": np.array(M["m_role"], np.int8), "m_e": np.array(M["m_e"], np.float32)}
    ncat = {c: int(sum(cc == c for cc in M["m_cat"])) for c in CATS}
    meta = dict(src=src, run=run, rt=rt, split=split, fn=str(fn), n_events=int(evsel.size), n_bench_split=int((b["ev_split"] == split).sum()), n_exotic_events=int(z["ev_number"].size),
                n_missing_file=int((~np.isin(z["ev_number"], bn)).sum()), n_rows_flag=int(ins.sum()), n_lg_seed=nlg, n_seed_not_max=int(nmax), n_dup_hits=int(ndup),
                n_members_flag_cat=ncat, n_clusters=n, n_clusters_D=nd, n_changed_events=n_chg, n_hits_dropped=n_hdrop, n_fate=n_fate.tolist(), n_members=len(M["m_e"]),
                wall=time.monotonic() - t0, ebeam=float(geo["ebeam"]), hz=float(geo["hz"]), gate=gate)
    print(f"  [{src}] run {run} split {split}: {evsel.size} events, {n} A / {nd} D(changed) clusters, {n_chg} changed events, fates {n_fate.tolist()} ({meta['wall']:.0f} s)", flush=True)
    return meta, dA, dD, f_, m_


def jobs13():
    """P7's jobs without the Ta sample cross-check files."""
    return [j for j in p7.jobs_of() if not (j[0] == "sample" and str(np.load(j[1])["run_type"]) == "2p24GeV_Ta_X17")]


def collect(nproc, cache=CACHE):
    """(metas, A, D, fate): A = default clusters (P7 REC dict + 'job'), D = clusters of the changed events only (REC + 'job'), fate = per A cluster arrays (FATE_Q; dmatch indexes D)
    and per flagged member arrays (MEM_Q; m_cl indexes A).  Read from the cache when present."""
    if cache and Path(cache).exists():
        c = np.load(cache, allow_pickle=True)
        print(f"[cache] {cache}", flush=True)
        return (list(c["metas"]), {q: c["a_" + q] for q in REC + ("job",)}, {q: c["d_" + q] for q in REC + ("job",)},
                {q: c["f_" + q] for q in FATE_Q + MEM_Q})
    t0 = time.monotonic()
    jobs = jobs13()
    with mp.get_context("fork").Pool(nproc) as pool:
        res = pool.map(work, sorted(jobs, key=lambda j: -os.path.getsize(j[1])), chunksize=1)
    metas = [r[0] for r in res]
    A = {q: np.concatenate([r[1][q] for r in res]) for q in REC}
    D = {q: np.concatenate([r[2][q] for r in res]) for q in REC}
    A["job"] = np.concatenate([np.full(r[0]["n_clusters"], i, np.int32) for i, r in enumerate(res)])
    D["job"] = np.concatenate([np.full(r[0]["n_clusters_D"], i, np.int32) for i, r in enumerate(res)])
    ao = np.r_[0, np.cumsum([r[0]["n_clusters"] for r in res])]
    do = np.r_[0, np.cumsum([r[0]["n_clusters_D"] for r in res])]
    Fa = {q: np.concatenate([r[3][q] for r in res]) for q in FATE_Q}
    Fa["dmatch"] = np.concatenate([np.where(r[3]["dmatch"] >= 0, r[3]["dmatch"] + do[i], -1) for i, r in enumerate(res)])
    Fa.update({q: np.concatenate([r[4][q] + (ao[i] if q == "m_cl" else 0) for i, r in enumerate(res)]) for q in MEM_Q})
    print(f"[collect] {len(jobs)} splits, {A['e'].size} A clusters, {D['e'].size} D clusters (changed events), {time.monotonic() - t0:.0f} s", flush=True)
    if cache:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache, metas=np.array(metas, dtype=object), **{"a_" + q: v for q, v in A.items()}, **{"d_" + q: v for q, v in D.items()},
                            **{"f_" + q: v for q, v in Fa.items()})
    return metas, A, D, Fa


def check_p7(metas, A, p7cache=p7.SCAN_ROOT / "exotic_p7_cache.npz"):
    """Reproduction of the P7 cache by A: per split same count and same (evn, seed, eg); returns {n_splits, n_match_splits, n_clusters, n_equal, rate, bad}."""
    c = np.load(p7cache, allow_pickle=True)
    pm = list(c["metas"])
    pidx = {(m["fn"], m["split"]): i for i, m in enumerate(pm)}
    pj = c["d_job"]
    pe, ps, pg = c["d_evn"], c["d_seed"], c["d_eg"]
    n_s = n_ok = n_c = n_eq = 0
    bad = []
    for i, m in enumerate(metas):
        k = (m["fn"], m["split"])
        if k not in pidx:
            bad.append((k, "not in P7"))
            continue
        n_s += 1
        sa = A["job"] == i
        sp = pj == pidx[k]
        n_c += int(sa.sum())
        same = int(sa.sum()) == int(sp.sum()) and bool(np.all(A["evn"][sa] == pe[sp]) and np.all(A["seed"][sa] == ps[sp]) and np.array_equal(A["eg"][sa], pg[sp]))
        n_ok += same
        n_eq += int(sa.sum()) if same else 0
        if not same:
            bad.append((k, f"A {int(sa.sum())} vs P7 {int(sp.sum())}"))
    return dict(n_splits=n_s, n_match_splits=n_ok, n_clusters=n_c, n_equal=n_eq, rate=n_eq / max(n_c, 1), bad=bad)


# ---------------------------------------------------------------------------------------------------------------- pass 2: tables
def names_of():
    with np.load(next((p7.SCAN_ROOT / "exotic_sample").glob("run_*.npz"))) as z:
        return [str(x) for x in z["mod_name"]]


def wil(k, n):
    lo, hi = p4.wilson(k, n)
    return float(lo), float(hi)


def wf(name, k, n):
    """{name_n, name, name_lo, name_hi}: count, Wilson 68 % fraction of n."""
    lo, hi = wil(k, n)
    return {f"n_{name}": int(k), name: k / n if n else np.nan, f"{name}_lo": lo if n else np.nan, f"{name}_hi": hi if n else np.nan}


def qs(x, qq=(0.5, 0.84, 0.99)):
    return np.quantile(x, qq) if np.size(x) else (np.nan,) * len(qq)


def qerr(x, qq, nb=200, seed=1):
    """Bootstrap (nb resamples, fixed seed) standard errors of the quantiles qq of x."""
    x = np.asarray(x, float)
    if x.size < 2:
        return (np.nan,) * len(qq)
    rng = np.random.default_rng(seed)
    n = min(x.size, 20000)
    xs = x if x.size <= n else rng.choice(x, n, replace=False)
    bs = np.array([np.quantile(rng.choice(xs, xs.size), qq) for _ in range(nb)])
    return tuple(bs.std(0, ddof=1) * np.sqrt(xs.size / x.size) if x.size > n else bs.std(0, ddof=1))


def se(x):
    x = np.asarray(x, float)
    return float(x.std(ddof=1) / np.sqrt(x.size)) if x.size > 1 else np.nan


def ratio(num, den):
    """(sum num / sum den, ratio-estimator error) over independent clusters."""
    sd = float(den.sum())
    if not sd:
        return np.nan, np.nan
    r = float(num.sum()) / sd
    return r, float(np.sqrt(np.sum((num - r * den) ** 2))) / sd


def window_mask(fits, X, base, rt, ring_lim):
    """Elastic window of P7's populations() with FIXED mu / sigma (fits of A: [all rings, ring 0, ...])."""
    sing = base & X["single"]
    if rt in p7.MOLLER_VETO_RT:
        with np.errstate(invalid="ignore"):
            mid = 0.5 * (X["rm"] + 1.0)
            sing &= (mid < p7.MOLLER_SEP_MAX) & (X["r"] > mid)
    el = np.zeros(base.size, bool)
    for i in range(ring_lim):
        m = sing & (X["ring"] == i)
        el |= m & (np.abs(X["r"] - fits[i + 1][1]) < p7.WIN_NSIG * fits[i + 1][2])
    return el, sing


def p10_ref():
    """{(dataset, population, axis, bin, veto): (frac, lo, hi)} of the P10 TSV, elastic population."""
    out = {}
    if not P10_TSV.exists():
        return out
    with open(P10_TSV) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["table"] == "veto" and r["population"] == "elastic" and r["veto"] in ("V3", "V4"):
                out[(r["dataset"], r["axis"], r["bin"], r["veto"])] = tuple(float(r[q]) for q in ("frac", "lo", "hi"))
    return out


def tables(metas, A, D, fate):
    """Long table rows (table, section P7 | P8 | P10, run_type, sample, ...) of the drop policy."""
    names = names_of()
    XA = p7.derive(metas, A, names)
    XD = p7.derive(metas, D, names)
    mol = p7.moller_keys()
    nring = len(RADBINS) - 1
    ringl = [f"{RADBINS[i]:g}-{RADBINS[i + 1]:g}" for i in range(nring)]
    ref = p10_ref()
    fa, egD, erem, ndrop, cdrop, dm = (fate[q] for q in ("fate", "egD", "erem", "ndrop", "cdrop", "dmatch"))
    egA = A["eg"].astype(float)
    ef = np.where(egA > 0, erem / np.maximum(egA, 1e-9), 0.0)
    fake = ((A["scat"] == 4) | (A["scat"] == 5)) & (XA["f_seed"] > NONSHOWER_F)
    eflag = A["es"].astype(float) + A["en"].astype(float).sum(1)
    eD = np.where(np.isfinite(egD), egD, 0.0).astype(float)
    ring = XA["ring"]
    rows = []
    for ds in p7.datasets(metas):
        key, lab, rt, src = ds
        P, fits, mst = p7.populations(metas, A, XA, ds, mol)
        base_d = (XD["rt"] == rt) & (XD["src"] == src)
        for smp in SAMPLES:
            if smp not in P:
                continue
            pm = P[smp]
            ids = [("overall", "all", pm)]
            if smp == "elastic":
                ids += [("ring", ringl[i], pm & (ring == i)) for i in range(nring)]
            for ax, bn, m in ids:
                n = int(m.sum())
                cm = {"table": "p7_cluster", "run_type": rt, "sample": smp, "axis": ax, "bin": bn, "n": n}
                aff = m & (ndrop > 0)
                sl = m & cdrop
                van, res_, f3 = m & (fa == 2), m & (fa == 1), m & (fa == 3)
                q = qs(ef[aff])
                fk = m & fake
                r_ = np.where(egA[res_] > 0, egD[res_] / np.maximum(egA[res_], 1e-9), np.nan)
                rows.append(dict(cm, section="P7", **wf("affected", int(aff.sum()), n), **wf("seedlost", int(sl.sum()), n), **wf("vanish", int(van.sum()), n),
                                 **wf("reseed", int(res_.sum()), n), **wf("merged", int((m & (fa == 4)).sum()), n), n_fate3=int(f3.sum()), reseed_egratio_med=float(np.nanmedian(r_)) if r_.size else np.nan,
                                 ef_med=q[0], ef_q84=q[1], ef_q99=q[2], **dict(zip(("ef_med_err", "ef_q84_err", "ef_q99_err"), qerr(ef[aff], (0.5, 0.84, 0.99)))),
                                 ef_mean=float(ef[aff].mean()) if aff.any() else np.nan, ef_mean_se=se(ef[aff]), reseed_egratio_med_err=qerr(r_, (0.5,))[0], n_fake=int(fk.sum()),
                                 **wf("fake_removed", int((fk & ((fa == 1) | (fa == 2) | (fa == 4))).sum()), n),
                                 frac_fake_removed_of_fake=float((fk & ((fa == 1) | (fa == 2) | (fa == 4))).sum() / fk.sum()) if fk.sum() else np.nan))
                if ax == "overall":
                    # P8 section: sample-wide scale shifts and cluster dE/E of affected clusters
                    s1, s1e = ratio(-eflag[m], egA[m])
                    s2, s2e = ratio(eD[m] - egA[m], egA[m])
                    dE1 = -ef[aff]
                    ok2 = aff & np.isfinite(egD)
                    dE2 = egD[ok2] / np.maximum(egA[ok2], 1e-9) - 1.0
                    q1, q2 = qs(dE1, (0.16, 0.5, 0.84)), qs(dE2, (0.16, 0.5, 0.84))
                    rows.append(dict(table="p8_scale", section="P8", run_type=rt, sample=smp, n=n, S1=s1, S1_err=s1e, S2=s2, S2_err=s2e,
                                     S1_from_erem=ratio(-erem[m].astype(float), egA[m])[0], n_aff=int(aff.sum()), dE1_q16=q1[0], dE1_med=q1[1], dE1_q84=q1[2],
                                     n_surv_aff=int(ok2.sum()), dE2_q16=q2[0], dE2_med=q2[1], dE2_q84=q2[2],
                                     dE1_med_err=qerr(dE1, (0.5,))[0], dE1_mean=float(dE1.mean()) if dE1.size else np.nan, dE1_mean_se=se(dE1),
                                     dE2_med_err=qerr(dE2, (0.5,))[0], dE2_mean=float(dE2.mean()) if dE2.size else np.nan, dE2_mean_se=se(dE2), **wf("seeddrop", int(sl.sum()), n)))
        # P8: flagged members dropped per category x role x ring (A clusters of this dataset)
        mc = fate["m_cl"]
        for smp in ("all", "elastic"):
            msk = P[smp][mc]
            for c in CATS:
                for role, rv in (("seed", 1), ("neighbour", 0)):
                    sel = msk & (fate["m_cat"] == c) & (fate["m_role"] == rv)
                    for ax, bn, sm in [("overall", "all", sel)] + [("ring", ringl[i], sel & (ring[mc] == i)) for i in range(nring)]:
                        ee = fate["m_e"][sm].astype(float)
                        rows.append(dict(table="p8_member", section="P8", run_type=rt, sample=smp, category=f"f{c}", role=role, axis=ax, bin=bn, n_dropped=int(sm.sum()),
                                         e_med=float(np.median(ee)) if ee.size else np.nan, e_mean=float(ee.mean()) if ee.size else np.nan, e_mean_se=se(ee), e_med_err=qerr(ee, (0.5,))[0], e_sum=float(ee.sum())))
        # P10 (ep): A elastic-window clusters under D
        elA = P["elastic"]
        elD, singD = window_mask(fits, XD, base_d, rt, nring)
        ia = np.flatnonzero(elA)
        d = dm[ia]
        chg = fate["chg"][ia]
        dj = np.where(d >= 0, d, 0)
        has = (d >= 0) | ~chg                                   # D cluster exists (same event unchanged -> A itself)
        merged = ~has & (fa[ia] == 4)
        vanish = ~has & ~merged
        nosingle = has & chg & ~singD[dj]
        outwin = has & chg & singD[dj] & ~elD[dj]
        lost = vanish | merged | nosingle | outwin
        matched_el = np.zeros(elD.size, bool)
        matched_el[d[(d >= 0) & ~lost]] = True
        entry = elD & ~matched_el                               # D window clusters (changed events) that are not an A window cluster
        # D window clusters whose A match was lost cannot exist (lost = not in window); an A cluster that is non-elastic in A but matched to an elastic D is an entry
        for ax, bn, sm, ent in [("overall", "all", np.ones(ia.size, bool), np.ones(elD.size, bool))] + \
                [("ring", ringl[i], ring[ia] == i, XD["ring"] == i) for i in range(nring)]:
            n = int(sm.sum())
            ne = int((entry & ent & base_d).sum())
            nd_ = int(((elD & ent & base_d).sum())) + int((sm & ~chg).sum())      # unchanged-event A window clusters are D window clusters
            r = dict(table="p10_elastic_loss", section="P10", run_type=rt, sample="elastic", axis=ax, bin=bn, n_A=n, n_D=nd_, **wf("lost", int((lost & sm).sum()), n),
                     **wf("vanish", int((vanish & sm).sum()), n), **wf("merged", int((merged & sm).sum()), n), **wf("nosingle", int((nosingle & sm).sum()), n), **wf("outwin", int((outwin & sm).sum()), n),
                     **wf("enter", ne, nd_))
            for v in ("V3", "V4"):
                f = ref.get((rt, ax, "all" if ax == "overall" else bn, v))
                if f:
                    r.update({f"{v}_frac": f[0], f"{v}_lo": f[1], f"{v}_hi": f[2]})
            rows.append(r)
    return rows


# ---------------------------------------------------------------------------------------------------------------- pages
def pc(x, nd=3):
    return "    n/a" if x is None or not np.isfinite(x) else f"{100 * x:{nd + 4}.{nd}f}"


def pe(v, e):
    return "n/a" if v is None or not np.isfinite(v) else f"{100 * v:.3f}+-{100 * e:.3f}"


def ci(r, k):
    return f"{pc(r[k])} [{pc(r[k + '_lo'])},{pc(r[k + '_hi'])}]"


def cik(k, n):
    """Wilson 68 % fraction k / n as 'f [lo,hi]' (%)."""
    if not n:
        return "n/a"
    lo, hi = wil(int(k), int(n))
    return f"{pc(k / n)} [{pc(lo)},{pc(hi)}]"


def civ(r, v):
    """P10 veto fraction of row r ('V3' / 'V4') with its Wilson interval."""
    f = r.get(f"{v}_frac")
    return "n/a" if f is None or not np.isfinite(f) else f"{pc(f)} [{pc(r[v + '_lo'])},{pc(r[v + '_hi'])}]"


def coverage_lines(metas):
    L = ["Coverage: benchmark split events present in the exotic scan of the split (join = intersection)", f"{'run type':<18}{'splits':>7}{'bench events':>14}{'exotic events':>15}{'joined':>12}{'joined/bench %':>16}"]
    for rt in p1.TYPES:
        ms = [m for m in metas if m["rt"] == rt]
        if ms:
            nb, ne, nj = (sum(m[q] for m in ms) for q in ("n_bench_split", "n_exotic_events", "n_events"))
            L.append(f"{RT_LABEL[rt]:<18}{len(ms):>7d}{nb:>14d}{ne:>15d}{nj:>12d}{100 * nj / max(nb, 1):>16.1f}")
    return L


def num(r):
    """Row with numeric strings (a row read back from a TSV) converted: ints, floats, '' -> nan."""
    out = {}
    for k, v in r.items():
        if isinstance(v, str) and k not in ("table", "section", "run_type", "sample", "axis", "bin", "category", "role"):
            try:
                v = np.nan if v == "" else (int(v) if v.lstrip("-").isdigit() else float(v))
            except ValueError:
                pass
        out[k] = v
    return out


def ep_cost(rows):
    """{run type: (k, n, frac, lo, hi)}: overall ep elastic-window clusters of A lost under D (all reasons)."""
    out = {}
    for r in rows:
        if r.get("table") == "p10_elastic_loss" and r["axis"] == "overall":
            r = num(r)
            out[r["run_type"]] = (int(r["n_lost"]), int(r["n_A"]), float(r["lost"]), float(r["lost_lo"]), float(r["lost_hi"]))
    return out


def pages(pdf, plt, rows, section, metas=None):
    """Text pages of one section (P7 | P8 | P10); metas (optional) adds the coverage table to the P7 page."""
    R = [num(r) for r in rows if r["section"] == section and str(r.get("table", "")).startswith(("p7_", "p8_", "p10_"))]
    rts = [rt for rt in p1.TYPES if any(r["run_type"] == rt for r in R)]
    if section == "P7":
        L = ["Channel-level DROP of every f1-f5 flagged channel before clustering; clusters recomputed.  fractions of the sample, Wilson 68 % in [ ].", "",
             f"{'run type':<16}{'sample':<9}{'N':>9}  {'affected %':<24}{'seed lost %':<24}{'vanish %':<24}{'re-seed %':<24}{'merged %':<24}{'fate3':>6}{'eg_D/eg_A re-seed':>22}"]
        for r in [r for r in R if r["axis"] == "overall"]:
            L.append(f"{RT_LABEL[r['run_type']]:<16}{r['sample']:<9}{r['n']:>9d}  {ci(r, 'affected'):<24}{ci(r, 'seedlost'):<24}{ci(r, 'vanish'):<24}{ci(r, 'reseed'):<24}{ci(r, 'merged'):<24}{r['n_fate3']:>6d}{r['reseed_egratio_med']:>14.3f}+-{r['reseed_egratio_med_err']:.3f}")
        L += ["", "Removed energy per affected cluster / E_g (%): median, 84 %, 99 % (bootstrap errors); mean (SE); fakes (f4/f5 seed, seed fraction > 0.95) removed (vanish, re-seed or merged)", "",
              f"{'run type':<16}{'sample':<9}{'median+-err':>16}{'84 %+-err':>16}{'99 %+-err':>16}{'mean+-SE':>16}{'N fake':>9}{'fake removed':>14}  {'% of sample':<24}{'% of fakes':>26}"]
        for r in [r for r in R if r["axis"] == "overall"]:
            L.append(f"{RT_LABEL[r['run_type']]:<16}{r['sample']:<9}{pe(r['ef_med'], r['ef_med_err']):>16}{pe(r['ef_q84'], r['ef_q84_err']):>16}{pe(r['ef_q99'], r['ef_q99_err']):>16}{pe(r['ef_mean'], r['ef_mean_se']):>16}{r['n_fake']:>9d}{r['n_fake_removed']:>14d}  {ci(r, 'fake_removed'):<24}{cik(r['n_fake_removed'], r['n_fake']):>26}")
        text_page(pdf, plt, "P7 under drop: affected clusters, seed loss, energy removed, fakes", L, size=9)
        L = ["Elastic-window clusters per P4 ring (seed module radius, mm)", "",
             f"{'run type':<16}{'ring':<10}{'N':>8}  {'affected %':<24}{'seed lost %':<24}{'vanish %':<24}{'re-seed %':<24}{'merged %':<24}{'ef med %':>9}{'fake rm':>8}"]
        for r in [r for r in R if r["axis"] == "ring"]:
            L.append(f"{RT_LABEL[r['run_type']]:<16}{r['bin']:<10}{r['n']:>8d}  {ci(r, 'affected'):<24}{ci(r, 'seedlost'):<24}{ci(r, 'vanish'):<24}{ci(r, 'reseed'):<24}{ci(r, 'merged'):<24}{pc(r['ef_med']):>9}{r['n_fake_removed']:>8d}")
        if metas:
            L += [""] + coverage_lines(metas)
        text_page(pdf, plt, "P7 under drop: elastic window by ring; coverage", L, size=8.5)
    elif section == "P8":
        S = [r for r in R if r["table"] == "p8_scale"]
        L = ["First order: S1 = -sum E_flagged / sum E_cluster (ratio-estimator error); second order S2 = (sum E_D - sum E_A) / sum E_A (vanished = 0); in per mille", "",
             f"{'run type':<16}{'sample':<9}{'N':>9}{'S1 (permille)':>16}{'+-':>8}{'S2 (permille)':>16}{'+-':>8}  {'seed dropped %':<24}{'affected N':>11}"
             f"{'dE/E 1st: 16/50(+-)/84 %':>34}{'dE/E 2nd (survivors)':>34}{'mean 1st+-SE %':>18}"]
        for r in S:
            L.append(f"{RT_LABEL[r['run_type']]:<16}{r['sample']:<9}{r['n']:>9d}{1e3 * r['S1']:>16.3f}{1e3 * r['S1_err']:>8.3f}{1e3 * r['S2']:>16.3f}{1e3 * r['S2_err']:>8.3f}  "
                     f"{ci(r, 'seeddrop'):<24}{r['n_aff']:>11d}   {100 * r['dE1_q16']:6.1f}/{100 * r['dE1_med']:6.1f}(+-{100 * r['dE1_med_err']:.1f})/{100 * r['dE1_q84']:6.1f}   {100 * r['dE2_q16']:6.1f}/{100 * r['dE2_med']:6.1f}(+-{100 * r['dE2_med_err']:.1f})/{100 * r['dE2_q84']:6.1f}   {pe(r['dE1_mean'], r['dE1_mean_se']):>16}")
        text_page(pdf, plt, "P8 under drop: energy-scale shift and cluster dE/E of affected clusters", L, size=8.5)
        L = ["Flagged members dropped (all PbWO4-seeded clusters), their full gain-corrected default energy (MeV): N, median, mean", "",
             f"{'run type':<16}{'cat':<5}{'role':<11}{'N':>9}{'med+-err':>14}{'mean+-SE':>14}   ring N (mm: " + " ".join(f"{RADBINS[i]:g}-{RADBINS[i + 1]:g}" for i in range(len(RADBINS) - 1)) + ")"]
        for rt in rts:
            for c in CATS:
                for role in ("seed", "neighbour"):
                    r = next(r for r in R if r["table"] == "p8_member" and r["run_type"] == rt and r["sample"] == "all" and r["category"] == f"f{c}" and r["role"] == role and r["axis"] == "overall")
                    rg = [next(q for q in R if q["table"] == "p8_member" and q["run_type"] == rt and q["sample"] == "all" and q["category"] == f"f{c}" and q["role"] == role and q["bin"] == f"{RADBINS[i]:g}-{RADBINS[i + 1]:g}")["n_dropped"]
                          for i in range(len(RADBINS) - 1)]
                    L.append(f"{RT_LABEL[rt]:<16}f{c:<4}{role:<11}{r['n_dropped']:>9d}{r['e_med']:>8.1f}+-{r['e_med_err']:<5.1f}{r['e_mean']:>8.1f}+-{r['e_mean_se']:<5.1f}   " + " ".join(f"{v:>8d}" for v in rg))
        text_page(pdf, plt, "P8 under drop: flagged members removed per category and role", L, size=8.5)
    else:
        L = ["A ep elastic-window clusters lost when the flagged channels are dropped (A's per-ring mu / sigma, no refit): vanish = seed dropped, no re-seed; merged = members absorbed by a neighbour cluster; no-single = event no longer",
             "passes exactly-one-cluster / SIGNAL_MIN / SEED_MIN; out-window = re-formed cluster outside +-3 sigma.  entering = D window clusters that were not an A window cluster.  Wilson 68 % in [ ].", "",
             f"{'run type':<16}{'ring':<10}{'N_A':>8}  {'lost %':<24}{'vanish %':<24}{'merged %':<24}{'no-single %':<24}{'out-window %':<24}{'enter N':>8}  {'entering % of D':<24}{'P10 V3 %':>26}{'P10 V4 %':>26}"]
        for r in R:
            L.append(f"{RT_LABEL[r['run_type']]:<16}{r['bin']:<10}{r['n_A']:>8d}  {ci(r, 'lost'):<24}{ci(r, 'vanish'):<24}{ci(r, 'merged'):<24}{ci(r, 'nosingle'):<24}{ci(r, 'outwin'):<24}{r['n_enter']:>8d}  {ci(r, 'enter'):<24}"
                     f"{civ(r, 'V3'):>26}{civ(r, 'V4'):>26}")
        text_page(pdf, plt, "P10 under drop (ep elastic): window loss and entries, beside the P10 V3 / V4 vetoes", L, size=8.5)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--cache", default=str(CACHE))
    ap.add_argument("--tsv", default="")
    ap.add_argument("--pdf", default="")
    ap.add_argument("--check-p7", action="store_true", help="compare A with the P7 cache")
    a = ap.parse_args()
    metas, A, D, F = collect(a.jobs, a.cache)
    cargs = (metas, A, D, F)
    if a.check_p7:
        c = check_p7(metas, A)
        print(f"[check] A vs P7 cache: {c['n_match_splits']}/{c['n_splits']} splits identical, {c['n_equal']}/{c['n_clusters']} clusters ({100 * c['rate']:.4f} %); {c['bad'][:5]}", flush=True)
    rows = tables(*cargs)
    if a.tsv:
        cols = []
        for r in rows:
            cols += [k for k in r if k not in cols]
        write_tsv(a.tsv, rows, cols)
        print(f"[tsv] {a.tsv} ({len(rows)} rows)", flush=True)
    if a.pdf:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.backends.backend_pdf import PdfPages
        with PdfPages(a.pdf) as pdf:
            for s in ("P7", "P8", "P10"):
                pages(pdf, plt, rows, s, metas)
        print(f"[pdf] {a.pdf}", flush=True)


if __name__ == "__main__":
    main()
