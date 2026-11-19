#!/usr/bin/env python3
"""exotic_p7_cluster_level.py -- P7 of the exotic-events impact packages: the exclusive exotic categories joined to reconstructed clusters.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p7_cluster_level.py [--jobs 8] [--out <dir>] [--cache <npz>] [--plot-only]
                                                                                       [--examples 10] [--evio-dir <dir>]

Join.  Clusters are the benchmark's own (energy_resolution_benchmark.py, control v3): the hits of the benchmark scan file of the same run
(output/energy_resolution/scan/<run type>/run_<run>.npz: matched-filter energy, hit time, 100 < t < 200 ns) with the residual onset-phase factor
(PHASE_RESID_BY_RT) and the per-module time offsets, clustered by erb.clusterers' gated HyCalCluster at the run type's CONTROL_GATE (builder
configuration, min_center_energy 30 MeV, Geant4 profile, log_weight_thres per run type).  Every exotic-scan file split (pileup_rates scan,
PILEUP_EXOTIC=1) is also a benchmark-scan split, so each exotic-scan event is found in the benchmark scan by its EVIO event number, and each
cluster member module gets the exclusive category of its exotic-scan row (run, event number, module) from exotic_p4_module_rates.flags_of with
the uniform P2 thresholds (f4 > f5 > f1 > f2 > f3; f2 never on a clipped or truncated row); a member without a kept exotic row, or with an
unflagged one, is unflagged.  PbWO4-seeded clusters only (seed module = the clusterer's centre, used for the ring and the Moller match).  Role: seed = the largest
member by gain-corrected energy (the centre in all but 0.1-0.4 % of the clusters, counted on page 1), neighbour = any other member.  The gated path
makes every in-window pulse a hit, so a module with two in-time pulses can seed two clusters; both get its category.  The benchmark's report-side below-threshold virtual members and position correction are not applied (they carry no flag).

Energies.  Member energies are multiplied by the control-v3 per-run gains (v3_defaults/<run type>/gains.tsv, default arm, moller = 0,
nonlinear = 0; a module without a row takes the run's median gain): E_g = sum g_i E_i.  The flagged energy fraction of a cluster is the sum of
g_i E_i over flagged members / E_g.  Energy bins use E_g / E_beam.

Populations (per run type).  all = every PbWO4-seeded cluster; hiE = E_g > 0.5 E_beam; elastic = events with exactly one cluster above
SIGNAL_MIN kE (500 MeV at 2.24 GeV, the benchmark selection), seed energy >= SEED_MIN kE, and |r - mu| < 3 sigma, r = E_g / E'(theta), theta from
the cluster position (no beam offset), E' the elastic energy on H2 (Ta for X17) minus the energy loss (benchmark e_elastic); mu, sigma from an
iterated Gaussian within 2 sigma (benchmark gauss_2sigma) per run type and radius ring, at 0.73 GeV only where the ep and Moller peaks separate and after a Moller veto (MOLLER_VETO_RT) (the run-type fit where a ring has < MIN_WIN_FIT clusters);
moller = the clusters of the selected Moller pairs of moller_calib.select (control settings, gate 4 ns, 60 mm isolation) read from its candidate
cache (output/energy_resolution/moller_x17/pairs_g4.npz), Ta X17 only, matched by (run, event number, seed module).
Rings: the P4 radius bins of the seed module centre (exotic_p4_module_rates.RADBINS).

Data.  H2 run types: the 2-file sample (output/pileup_rates/exotic_sample, two splits per run).  Ta X17: the 19 full-scan npz
(output/pileup_rates/exotic_full) as the main source; the Ta sample files (2 splits of 25800, 26138, 26230) as a cross-check row.  The H2
full-scan pickles carry no event ids and cannot be joined.  Run 24937 excluded.

Outputs (--out): exotic_P7_cluster_level_20261118.pdf, .tsv (long table: source, run type, population, axis, bin, category, role, counts,
Wilson 68 % fraction, energy fractions, fake-cluster rates), _examples.tsv (the example events).  The pass-1 cluster records are an intermediate
cache (--cache, default output/pileup_rates/exotic_p7_cache.npz; delete it to rebuild).
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
sys.path.insert(0, str(HERE.parent.parent))             # waveform_analysis: energy_resolution_benchmark, moller_calib
sys.path.insert(0, str(HERE.parent))                    # pileup_analysis: pileup_rates
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import exotic_p1_quantities as p1                       # noqa: E402
import exotic_p4_module_rates as p4                     # noqa: E402
from exotic_p1_quantities import TYPES, RT_LABEL, write_tsv, EXCLUDE_RUNS     # noqa: E402
from exotic_p3_sample import text_page                  # noqa: E402

STEM = "exotic_P7_cluster_level_20261118"
CATS = (4, 5, 1, 2, 3)                                  # exclusive categories in priority order (p4.CAT_ORDER)
CNAME = {0: "any", 4: "f4 clipped-multi", 5: "f5 saturated single", 1: "f1 tail-baseline", 2: "f2 wide", 3: "f3 smooth-wide"}
CCOL = {4: "tab:red", 5: "tab:orange", 1: "tab:blue", 2: "tab:green", 3: "tab:purple"}
RADBINS = p4.RADBINS                                    # mm, seed module centre radius (P4 rings)
EBINS = (0.0, 0.1, 0.2, 0.35, 0.5, 0.7, 0.9, 1.1, np.inf)   # E_g / E_beam
WIN_NSIG = 3.0                                          # elastic window half width (sigma)
MIN_WIN_FIT = 500                                       # single clusters per ring for a ring-own elastic window fit
EF_LEVEL = 1e-3                                         # mean flagged energy fraction that matters for the calibration (0.1 %)
MOLLER_VETO_RT = ("0p73GeV_H2",)                     # elastic window: only angles where the ep and Moller peaks separate, i.e. the midpoint of E_moller(theta) / E' and 1 is
MOLLER_SEP_MAX = 0.90                                   # below this (~ mu - 3 sigma at 0.73 GeV, sigma ~ 3.3 %; seed radius >~ 110 mm), and above that midpoint (Moller veto), before the fit
MIN_CELL = 5                                            # flagged clusters a cell needs to be judged against EF_LEVEL
SCAN_ROOT = p1.PYDIR / "output" / "pileup_rates"
BENCH_SCAN = p1.PYDIR / "output" / "energy_resolution" / "scan"
POPS = ("all", "hiE", "elastic", "moller")
POP_LABEL = {"all": "all clusters", "hiE": "E_g > 0.5 E_beam", "elastic": "ep elastic window +-3 sigma", "moller": "Moller pairs (selected)"}
REC = ("evn", "split", "seed", "seed_e", "e", "eg", "x", "y", "nbig", "cnt", "scat", "sraw", "nmask", "es", "en")   # per cluster


def erb_mod():
    import energy_resolution_benchmark as erb
    return erb


# ---------------------------------------------------------------------------------------------------------------- pass 1: per split
def row_categories(z):
    """(key = event number * 4096 + module index, exclusive category) of the flagged PbWO4 rows of an exotic-scan npz (uniform thresholds),
    plus the raw-flag bit mask (bit i-1 = raw f_i) of the same rows and the number of kept PbWO4 rows."""
    midx = z["ch_midx"].astype(np.int64)
    sel = z["mod_pwo4"][midx]
    k = z["ch_k"].astype(np.int64)
    a, t = z["ch_a"], z["ch_t"]
    main = np.argmax(np.where(np.isfinite(t), a, -1.0), axis=1)
    t0 = t[np.arange(k.size), main].astype(np.float64)
    with np.errstate(invalid="ignore"):
        trunc = (np.ceil(t0 / p4.DT) - p4.LEAD + p4.LW > p4.NS) & np.isfinite(t0)
    c = dict(db=z["ch_ex_db"].astype(np.float64), tail=z["ch_ex_tail"].astype(np.float64), w50=z["ch_w50"].astype(np.float64), k=k,
             nm=z["ch_ex_nm"].astype(np.int64), c2s=z["ch_ex_chi2s"].astype(np.float64), c2k=z["ch_ex_chi2k"].astype(np.float64),
             nclip=z["ch_nclip"].astype(np.int64), trunc=trunc, s4=z["ch_ex_f4"].astype(bool))
    F = p4.flags_of(c, p4.UNI)
    flag = sel & (F["cat"] > 0)
    raw = sum(F[i].astype(np.int64) << (i - 1) for i in (1, 2, 3, 4, 5))
    key = z["ev_number"][z["ch_ev"][flag].astype(np.int64)].astype(np.int64) * 4096 + midx[flag]
    return key, F["cat"][flag].astype(np.int8), raw[flag].astype(np.int8), int(sel.sum())


_GAIN = {}


def run_gains(rt, run, names):
    """[n modules] control-v3 gain of each module for run (default arm), the run's median for a module without a row."""
    if run not in _GAIN:
        import moller_calib as mc
        g = {}
        with open(mc.REF_ROOT / rt / "gains.tsv") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                if int(r["run"]) == run and r["moller"] == "0" and r["nonlinear"] == "0":
                    g[r["module"]] = float(r["g"])
        if not g:
            sys.exit(f"[ERROR] no control-v3 gains for run {run} in {mc.REF_ROOT / rt / 'gains.tsv'}")
        med = float(np.median(list(g.values())))
        _GAIN[run] = np.array([g.get(str(n), med) for n in names])
    return _GAIN[run]


def bench_file(rt, run):
    return BENCH_SCAN / rt / f"run_{run}.npz"


def split_hits(b, fn_b, split, evsel, geo, rt):
    """Benchmark hits of the events evsel (indices into b['ev_number']) of one split, with the residual onset-phase factor."""
    erb = erb_mod()
    h_ev = b["h_ev"]
    lo, hi = np.searchsorted(h_ev, [evsel[0], evsel[-1] + 1]) if evsel.size else (0, 0)
    h = {q: b[q][lo:hi] for q in ("h_mod", "h_t", "h_e", "h_k", "h_extra")} | {"h_ev": h_ev[lo:hi]}
    keep = np.isin(h["h_ev"], evsel)
    h = {q: v[keep] for q, v in h.items()}
    if erb.PHASE_RESID_BY_RT.get(rt):
        erb.PHASE_RESID = str(erb.TABLE_DIR / erb.PHASE_RESID_BY_RT[rt])
        h = erb.phase_resid(h, b, geo, split, str(fn_b))
    return h


def cluster_split(h, cl, geo, lut, flag_ev, gains, kE, evnum):
    """Every cluster of the events of h (gated path cl).  Returns REC lists, the number of LG-seeded clusters dropped, the keys of the flagged
    rows that are cluster members, the number of clusters whose centre is not their largest (gain-corrected) member module, and the number of
    second hits of a module within one cluster (summed into that module)."""
    erb = erb_mod()
    out = {q: [] for q in REC}
    ev = h["h_ev"]
    if not ev.size:
        return out, 0, set(), 0, 0
    lo = np.flatnonzero(np.r_[True, ev[1:] != ev[:-1]])
    hi = np.r_[lo[1:], ev.size]
    mod, e = h["h_mod"].tolist(), h["h_e"].tolist()
    t = (h["h_t"].astype(np.float64) - geo["toff"][h["h_mod"]]).tolist() if geo.get("toff") is not None else h["h_t"].tolist()
    pw, sig = geo["pwo4"], erb.SIGNAL_MIN * kE
    nlg = nmax = ndup = 0
    used = set()
    for s, f in zip(lo.tolist(), hi.tolist()):
        cl.clear()
        for i in range(s, f):
            cl.add_hit(mod[i], e[i], t[i])
        cl.form_clusters()
        res = cl.reconstruct_matched()
        nbig = sum(r.hit.energy > sig for r in res)
        en_ = int(evnum[ev[s]])
        fl = en_ in flag_ev
        for r in res:
            seed = r.cluster.center.index
            if not pw[seed]:
                nlg += 1
                continue
            em = {}                                     # gain-corrected member energy per module (the gated path may hold two pulses of one module)
            for x in r.cluster.hits:
                em[x.index] = em.get(x.index, 0.0) + x.energy * gains[x.index]
            ndup += len(r.cluster.hits) - len(em)
            mi, me = list(em), np.fromiter(em.values(), float, len(em))
            eg = float(me.sum())
            smax = mi[int(me.argmax())]                 # role "seed" = the largest (gain-corrected) member module
            nmax += smax != seed
            scat = sraw = nmask = 0
            es = 0.0
            enb = [0.0] * 6
            if fl:
                for m, x in zip(mi, me):
                    v = lut.get(en_ * 4096 + m)
                    if v is None:
                        continue
                    c, rw = v
                    used.add(en_ * 4096 + m)
                    if m == smax:
                        scat, sraw, es = c, rw, float(x)
                    else:
                        nmask |= 1 << (c - 1)
                        enb[c] += float(x)
            for q, v in zip(REC, (en_, 0, seed, r.cluster.center.energy, r.hit.energy, eg, r.hit.x, r.hit.y, nbig, len(mi), scat, sraw, nmask, es, enb[1:])):
                out[q].append(v)
    return out, nlg, used, nmax, ndup


def work(job):
    """One exotic-scan split joined to its benchmark-scan split: per-cluster records (REC) and counts."""
    src, fn, split = job
    erb = erb_mod()
    t0 = time.monotonic()
    z = np.load(fn)
    run, rt = int(z["run"]), str(z["run_type"])
    key, cat, raw, n_rows = row_categories(z)
    fn_b = bench_file(rt, run)
    b = np.load(fn_b)
    bn = b["ev_number"]
    evsel = np.flatnonzero((b["ev_split"] == split) & np.isin(bn, z["ev_number"]))
    gate = erb.CONTROL_GATE[rt]
    cls, geo = erb.clusterers(run, str(b["calib"]), (gate,))
    kE = round(geo["ebeam"] / erb.E_REF, 6)
    gains = run_gains(rt, run, geo["name"])
    h = split_hits(b, fn_b, split, evsel, geo, rt)
    ins = np.isin(key // 4096, bn[evsel])
    lut = {int(k_): (int(c), int(r)) for k_, c, r in zip(key[ins], cat[ins], raw[ins])}
    flag_ev = set((key[ins] // 4096).tolist())
    out, nlg, used, nmax, ndup = cluster_split(h, cls[f"g{gate:g}"], geo, lut, flag_ev, gains, kE, bn)
    n = len(out["e"])
    d = {q: np.array(v, dtype=np.float32 if q in ("seed_e", "e", "eg", "x", "y", "es", "en") else None) for q, v in out.items()}
    d["en"] = d["en"].reshape(n, 5) if n else np.zeros((0, 5), np.float32)
    d["split"] = np.full(n, split, np.int32)
    for q in ("seed", "nbig", "cnt", "scat", "sraw", "nmask"):
        d[q] = d[q].astype(np.int16)
    meta = dict(src=src, run=run, rt=rt, split=split, fn=str(fn), n_events=int(evsel.size), n_missing_file=int((~np.isin(z["ev_number"], bn)).sum()),
                n_rows_flag=int(ins.sum()), n_lg_seed=nlg, n_seed_not_max=int(nmax), n_dup_hits=int(ndup),
                n_rows_flag_cat={c: int((cat[ins] == c).sum()) for c in CATS},
                n_rows_in_cluster_cat={c: int(sum(lut[k_][0] == c for k_ in used)) for c in CATS}, n_clusters=n, wall=time.monotonic() - t0, ebeam=float(geo["ebeam"]), hz=float(geo["hz"]), gate=gate)
    print(f"  [{src}] run {run} split {split}: {evsel.size} events, {n} clusters, {ins.sum()} flagged rows ({meta['wall']:.0f} s)", flush=True)
    return meta, d


def jobs_of():
    jobs = []
    for f in sorted((SCAN_ROOT / "exotic_sample").glob("run_*.npz")):
        with np.load(f) as z:
            if int(z["run"]) in EXCLUDE_RUNS:
                continue
            jobs += [("sample", str(f), int(s)) for s in z["splits"]]
    for f in sorted((SCAN_ROOT / "exotic_full").glob("run_*_s*.npz")):
        with np.load(f) as z:
            if int(z["run"]) in EXCLUDE_RUNS:
                continue
            jobs += [("full", str(f), int(s)) for s in z["splits"]]
    return jobs


def collect(jobs, nproc, cache):
    if cache and Path(cache).exists():
        c = np.load(cache, allow_pickle=True)
        print(f"[cache] {cache}", flush=True)
        return list(c["metas"]), {q: c["d_" + q] for q in REC + ("job",)}
    t0 = time.monotonic()
    with mp.get_context("fork").Pool(nproc) as pool:
        res = pool.map(work, sorted(jobs, key=lambda j: -os.path.getsize(j[1])), chunksize=1)
    metas = [r[0] for r in res]
    D = {q: np.concatenate([r[1][q] for r in res]) for q in REC}
    D["job"] = np.concatenate([np.full(r[0]["n_clusters"], i, np.int32) for i, r in enumerate(res)])
    print(f"[collect] {len(jobs)} splits, {D['e'].size} clusters, {time.monotonic() - t0:.0f} s", flush=True)
    if cache:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache, metas=np.array(metas, dtype=object), **{"d_" + q: v for q, v in D.items()})
    return metas, D


# ---------------------------------------------------------------------------------------------------------------- pass 2: populations
def datasets(metas):
    """[(key, label, run type, source)]: H2 run types from the sample, Ta X17 from the full scan (main), Ta sample (cross-check)."""
    have = {(m["rt"], m["src"]) for m in metas}
    out = [(rt, RT_LABEL[rt] + " (sample)", rt, "sample") for rt in TYPES if rt != "2p24GeV_Ta_X17" and (rt, "sample") in have]
    if ("2p24GeV_Ta_X17", "full") in have:
        out.append(("2p24GeV_Ta_X17", RT_LABEL["2p24GeV_Ta_X17"] + " (full)", "2p24GeV_Ta_X17", "full"))
    if ("2p24GeV_Ta_X17", "sample") in have:
        out.append(("Ta_sample", RT_LABEL["2p24GeV_Ta_X17"] + " (sample)", "2p24GeV_Ta_X17", "sample"))
    return out


def fit_window(r):
    """(mu, sigma, n) of r = E_g / E': iterated Gaussian within 2 sigma (benchmark gauss_2sigma, in per-mille units, 0.2 % bins on 0.6-1.4)
    started at the mode of the 1 %-smoothed histogram on 0.8-1.2 with sigma 3 %."""
    erb = erb_mod()
    v = r[(r > 0.6) & (r < 1.4)]
    if v.size < 50:
        return np.nan, np.nan, int(v.size)
    y, e = np.histogram(v, bins=400, range=(0.6, 1.4))
    c = 0.5 * (e[1:] + e[:-1])
    ys = np.convolve(y, np.ones(5) / 5.0, mode="same")
    inr = (c > 0.8) & (c < 1.2)
    mu, sg = float(c[inr][np.argmax(ys[inr])]), 0.03
    m, s = erb.gauss_2sigma(1000.0 * c, y.astype(float), 1000.0 * mu, 1000.0 * sg)
    return m / 1000.0, s / 1000.0, int(v.size)


def moller_keys():
    """{(run, event number, seed module)} of the clusters of the selected Moller pairs (moller_calib.select on its candidate cache, gate path
    g4, control settings), plus {key: cluster energy} for the cross-check and the cache path."""
    import moller_calib as mc
    erb = erb_mod()
    f = mc.CACHE_ROOT / "pairs_g4.npz"
    if not f.exists():
        return None, {}, str(f)
    p = np.load(f, allow_pickle=False)
    d = {q[3:]: p[q] for q in p.files if q.startswith("g4:")}
    sp = {q[6:]: p[q] for q in p.files if q.startswith("split:")}
    run0 = int(d["run"][0])
    with np.load(bench_file(mc.RT, run0)) as b:
        geo = erb.clusterers(run0, str(b["calib"]), (4.0,))[1]
    s = mc.select(d, geo, sp, "g4", mc.CL_MIN, mc.R_ISO)
    sel = np.repeat(s["sel"], 2)
    keys, E = set(), {}
    for run in np.unique(d["run"][sel]).tolist():
        with np.load(bench_file(mc.RT, run)) as b:
            bn = b["ev_number"]
        m = sel & (d["run"] == run)
        for ev, seed, e in zip(bn[d["ev"][m].astype(np.int64)].tolist(), d["seed"][m].tolist(), d["e"][m].tolist()):
            keys.add((run, ev, int(seed)))
            E[(run, ev, int(seed))] = e
    return keys, E, str(f)


def derive(metas, D, names):
    """Per-cluster derived arrays: run type, source, run, beam energy, ring, E bin, single, r, flagged energy fractions."""
    erb = erb_mod()
    J = D["job"]
    pick = lambda k, dt=None: np.array([m[k] for m in metas], dtype=dt)[J]      # noqa: E731
    X = dict(rt=pick("rt"), src=pick("src"), run=pick("run", np.int64), eb=pick("ebeam", float), hz=pick("hz", float))
    rad = p4.radius_of(names)
    seed = D["seed"].astype(np.int64)
    X["rad"] = rad[seed]
    X["ring"] = np.clip(np.digitize(X["rad"], RADBINS) - 1, 0, len(RADBINS) - 2)
    X["ebin"] = np.clip(np.digitize(D["eg"] / X["eb"], EBINS) - 1, 0, len(EBINS) - 2)
    kE = X["eb"] / erb.E_REF
    X["single"] = (D["nbig"] == 1) & (D["e"] > erb.SIGNAL_MIN * kE) & (D["seed_e"] >= erb.SEED_MIN * kE)
    th = erb.theta_deg(D["x"].astype(float), D["y"].astype(float), X["hz"])
    ep = np.zeros(th.size)
    for rt in np.unique(X["rt"]):
        m = X["rt"] == rt
        ep[m] = erb.e_elastic(th[m], X["eb"][m], erb.pool_mass(str(rt)))
    X["r"] = D["eg"] / np.where(ep > 0, ep, np.nan)
    X["rm"] = erb.e_moller(th, X["eb"]) / np.where(ep > 0, ep, np.nan)      # Moller E(theta) / E'
    eg = np.maximum(D["eg"].astype(float), 1e-9)
    X["f_seed"] = D["es"] / eg                                   # flagged seed energy fraction
    X["f_nb"] = D["en"].astype(float) / eg[:, None]              # [n, 5] neighbour flagged energy fraction per category f1..f5
    X["f_all"] = X["f_seed"] + X["f_nb"].sum(1)
    return X


def populations(metas, D, X, ds, mol):
    """{pop: bool mask} of dataset ds, the elastic-window fits [(ring label, mu, sigma, n, own fit)], Moller match stats."""
    key, _, rt, src = ds
    base = (X["rt"] == rt) & (X["src"] == src)
    P = {"all": base, "hiE": base & (D["eg"] > 0.5 * X["eb"])}
    sing = base & X["single"]
    if rt in MOLLER_VETO_RT:
        with np.errstate(invalid="ignore"):
            mid = 0.5 * (X["rm"] + 1.0)
            sing &= (mid < MOLLER_SEP_MAX) & (X["r"] > mid)
    mu0, s0, n0 = fit_window(X["r"][sing])
    fits = [("all rings", mu0, s0, n0, True)]
    el = np.zeros(base.size, bool)
    for i in range(len(RADBINS) - 1):
        m = sing & (X["ring"] == i)
        mu, s, n = fit_window(X["r"][m])
        own = n >= MIN_WIN_FIT and np.isfinite(mu)
        if not own:
            mu, s = mu0, s0
        fits.append((f"{RADBINS[i]:g}-{RADBINS[i + 1]:g} mm", mu, s, n, own))
        el |= m & (np.abs(X["r"] - mu) < WIN_NSIG * s)
    P["elastic"] = el
    mst = None
    if rt == "2p24GeV_Ta_X17" and mol[0] is not None:
        idx = np.flatnonzero(base)
        runs = X["run"][idx]
        k = list(zip(runs.tolist(), D["evn"][idx].tolist(), D["seed"][idx].astype(np.int64).tolist()))
        hit = np.array([q in mol[0] for q in k], bool)
        mm = np.zeros(base.size, bool)
        mm[idx[hit]] = True
        P["moller"] = mm
        # expected: selected Moller clusters whose (run, split) is in this dataset's joined splits
        splits = {(m["run"], m["split"]) for m in metas if m["rt"] == rt and m["src"] == src}
        evn_split = {}
        for j in idx[hit]:
            evn_split[(int(X["run"][j]), int(D["evn"][j]), int(D["seed"][j]))] = float(D["e"][j])
        de = np.array([evn_split[q] / mol[1][q] - 1.0 for q in evn_split])
        mst = dict(n_match=int(hit.sum()), splits=len(splits), de_max=float(np.max(np.abs(de))) if de.size else np.nan)
    return P, fits, mst


def poisson_wilson(k, z=1.0):
    """Wilson-score 68 % interval of a Poisson mean given k counts (finite at k = 0: [0, 1])."""
    c, h = k + 0.5 * z * z, z * np.sqrt(k + 0.25 * z * z)
    return max(c - h, 0.0), c + h


def frac_rows(lab, P, D, X, n_ev):
    """TSV rows (table frac) of one dataset: population x axis (overall, ring, ebin) x category (any, f4, f5, f1, f2, f3) x role (any, seed, neighbour)."""
    rows = []
    scat, nmask = D["scat"].astype(np.int64), D["nmask"].astype(np.int64)
    roles = {}
    for c in (0,) + CATS:
        if c == 0:
            s, nb = scat > 0, nmask > 0
            fs, fn = X["f_seed"], X["f_nb"].sum(1)
        else:
            s, nb = scat == c, (nmask >> (c - 1) & 1) > 0
            fs, fn = np.where(scat == c, X["f_seed"], 0.0), X["f_nb"][:, c - 1]
        roles[c] = {"any": (s | nb, fs + fn), "seed": (s, fs), "neighbour": (nb, fn)}
    axes = [("overall", [("all", None)])]
    axes.append(("ring", [(f"{RADBINS[i]:g}-{RADBINS[i + 1]:g}", ("ring", i)) for i in range(len(RADBINS) - 1)]))
    axes.append(("ebin", [(f"{EBINS[i]:g}-{EBINS[i + 1]:g}", ("ebin", i)) for i in range(len(EBINS) - 1)]))
    for pop, pm in P.items():
        for ax, bins in axes:
            for bl, cond in bins:
                m = pm if cond is None else pm & (X[cond[0]] == cond[1])
                n = int(m.sum())
                for c in (0,) + CATS:
                    for role, (flag, f) in roles[c].items():
                        k = int((m & flag).sum())
                        lo, hi = p4.wilson(k, n)
                        ff = f[m & flag]
                        rlo, rhi = poisson_wilson(k)
                        q = np.quantile(ff, (0.5, 0.84, 0.99)) if ff.size else (np.nan,) * 3
                        rows.append(dict(table="frac", dataset=lab, population=pop, axis=ax, bin=bl, category=CNAME[c], role=role, n_clusters=n, n_flag=k,
                                         frac=k / n if n else np.nan, frac_lo=float(lo) if n else np.nan, frac_hi=float(hi) if n else np.nan,
                                         ef_mean=float(f[m].sum() / n) if n else np.nan, ef_sem=float(np.std(f[m], ddof=1) / np.sqrt(n)) if n > 1 else np.nan, ef_med_flag=q[0], ef_q84_flag=q[1], ef_q99_flag=q[2],
                                         n_events=n_ev, per_kevent=1000.0 * k / n_ev if n_ev else np.nan,
                                         rate_lo68=1000.0 * rlo / n_ev if n_ev else np.nan, rate_hi68=1000.0 * rhi / n_ev if n_ev else np.nan))
    return rows


def get(rows, **kw):
    for r in rows:
        if all(r.get(k) == v for k, v in kw.items()):
            return r
    return None


# ---------------------------------------------------------------------------------------------------------------- examples
def pick_examples(P, D, X, metas, ds, n, rng):
    """n elastic-window clusters with a flagged seed, round robin over the seed categories present (random within a category)."""
    m = P["elastic"] & (D["scat"] > 0)
    pools = {c: list(rng.permutation(np.flatnonzero(m & (D["scat"] == c)))) for c in CATS}
    out = []
    while len(out) < n and any(pools.values()):
        for c in CATS:
            if pools[c] and len(out) < n:
                out.append(int(pools[c].pop()))
    return [dict(dataset=ds[0], rt=ds[2], run=int(X["run"][i]), split=int(D["split"][i]), evn=int(D["evn"][i]), seed=int(D["seed"][i]),
                 scat=int(D["scat"][i]), eg=float(D["eg"][i]), r=float(X["r"][i]), fn=metas[int(D["job"][i])]["fn"], ebeam=float(X["eb"][i]))
            for i in out]


_CTX = {}


def example_job(a):
    """Seed waveform (EVIO re-read, pileup_rates fit with the exotic scan's stored noise reference) and the cluster's members (benchmark hits
    reclustered as in pass 1) with the exotic category of every module of the event."""
    ex, evio_dir = a
    import event_display as ed
    erb = erb_mod()
    run, split, evn, seed = ex["run"], ex["split"], ex["evn"], ex["seed"]
    out = dict(ex)
    z = np.load(ex["fn"])
    key, cat, raw, _ = row_categories(z)
    sel = key // 4096 == evn
    out["mod_cat"] = {int(k % 4096): (int(c), int(r)) for k, c, r in zip(key[sel], cat[sel], raw[sel])}
    rt = ex["rt"]
    fn_b = bench_file(rt, run)
    b = np.load(fn_b)
    ie = int(np.flatnonzero(b["ev_number"] == evn)[0])
    gate = erb.CONTROL_GATE[rt]
    cls, geo = erb.clusterers(run, str(b["calib"]), (gate,))
    h = split_hits(b, fn_b, split, np.array([ie]), geo, rt)
    g = run_gains(rt, run, geo["name"])
    cl = cls[f"g{gate:g}"]
    cl.clear()
    t = h["h_t"].astype(np.float64) - (geo["toff"][h["h_mod"]] if geo.get("toff") is not None else 0.0)
    for m_, e_, t_ in zip(h["h_mod"].tolist(), h["h_e"].tolist(), t.tolist()):
        cl.add_hit(m_, e_, t_)
    cl.form_clusters()
    c = next((r for r in cl.reconstruct_matched() if r.cluster.center.index == seed), None)
    mem = {}
    for x in (c.cluster.hits if c is not None else []):
        mem[x.index] = mem.get(x.index, 0.0) + x.energy * g[x.index]
    out["members"] = list(mem.items())
    hm = {}
    for m_, e_ in zip(h["h_mod"].tolist(), h["h_e"].tolist()):
        hm[m_] = hm.get(m_, 0.0) + e_ * g[m_]
    out["hits"] = hm
    out["geo"] = {q: geo[q] for q in ("x", "y", "sx", "sy", "name", "pwo4")}
    # seed waveform
    if run not in _CTX:
        _CTX[run] = ed.Ctx(run)
    ctx = _CTX[run]
    f = Path(evio_dir) / f"prad_{run:06d}.evio.{split:05d}"
    out["wave"] = None
    fe = ed.find_event(ctx, f, evn) if f.exists() else None
    if fe is not None:
        j = list(z["splits"]).index(split)
        b_ref, sig_ref = z["mod_b_ref"][j].astype(float), z["mod_sig_ref"][j].astype(float)
        r = ed.fit_event(ctx, fe, b_ref, sig_ref)
        if seed in r["row"]:
            i = r["row"][seed]
            k = int(r["k"][i])
            out["wave"] = dict(y=r["Y"][i].astype(float) - r["ped"][i], comp=r["comp"][i][:k], k=k, chi2=float(r["chi2"][i]), a=r["a"][i][:k].tolist(),
                               t=r["t"][i][:k].tolist(), nclip=int(r["nclip"][i]))
    print(f"  [example] run {run} ev {evn} {geo['name'][seed]} f{ex['scat']}: wave {'ok' if out['wave'] else 'missing'}", flush=True)
    return out


# ---------------------------------------------------------------------------------------------------------------- report
def pct(x):
    return f"{max(100 * x, 0.0):.4f}" if np.isfinite(x) else "nan"          # fractions are >= 0 (Wilson lower edge rounding)


def ci(r):
    return f"{pct(r['frac'])} [{pct(r['frac_lo'])}, {pct(r['frac_hi'])}]" if r["n_clusters"] else "-"


def efs(r):
    """<fE> +- SEM in %."""
    return f"{pct(r['ef_mean'])} +- {100 * r['ef_sem']:.4f}" if r["n_clusters"] > 1 else "-"


def rate(r):
    """rate per 1000 events [Poisson-Wilson 68 %]."""
    return f"{r['per_kevent']:.3f} [{r['rate_lo68']:.3f}, {r['rate_hi68']:.3f}]"


def table_main(rows, dss, pops):
    g = lambda ds, pop, c, role: get(rows, dataset=ds[0], population=pop, axis="overall", category=CNAME[c], role=role)       # noqa: E731
    L = ["1a: clusters with >= 1 flagged module (any category): % [68 % CI]; <fE> +- SEM (% of E_g, mean over all clusters of the population)",
         f"{'dataset':<26}{'population':<12}{'N clusters':>10}  {'any role':<30}{'seed':<30}{'neighbour':<30}{'<fE> +- SEM':>18}{'fE|fl med/84/99 %':>20}"]
    for ds in dss:
        for pop in pops:
            r = g(ds, pop, 0, "any")
            if r is None:
                continue
            q = "/".join(f"{100 * r[k]:.1f}" if np.isfinite(r[k]) else "-" for k in ("ef_med_flag", "ef_q84_flag", "ef_q99_flag"))
            L.append(f"{ds[1]:<26}{pop:<12}{r['n_clusters']:>10d}  {ci(r):<30}{ci(g(ds, pop, 0, 'seed')):<30}{ci(g(ds, pop, 0, 'neighbour')):<30}{efs(r):>18}{q:>20}")
        L.append("")
    L.append("1b: clusters with >= 1 module of exclusive category c (any role): % [68 % CI]")
    L.append(f"{'dataset':<26}{'population':<12}" + "".join(f"{'f' + str(c):<26}" for c in CATS))
    for ds in dss:
        for pop in pops:
            if g(ds, pop, 0, "any") is None:
                continue
            L.append(f"{ds[1]:<26}{pop:<12}" + "".join(f"{ci(g(ds, pop, c, 'any')):<26}" for c in CATS))
        L.append("")
    return L


def table_roles(rows, dss, pop):
    L = [f"population: {POP_LABEL[pop]}.  Per category and role: fraction of clusters (%, Wilson 68 %), mean flagged energy fraction <fE> +- SEM (% of E_g, over all clusters of the population)",
         f"{'dataset':<26}{'category':<22}{'seed % [68 % CI]':<36}{'<fE> seed +- SEM':>20}   {'neighbour % [68 % CI]':<36}{'<fE> nb +- SEM':>20}"]
    for ds in dss:
        for c in CATS:
            s = get(rows, dataset=ds[0], population=pop, axis="overall", category=CNAME[c], role="seed")
            nb = get(rows, dataset=ds[0], population=pop, axis="overall", category=CNAME[c], role="neighbour")
            if s is None:
                continue
            L.append(f"{ds[1]:<26}{CNAME[c]:<22}{ci(s):<36}{efs(s):>20}   {ci(nb):<36}{efs(nb):>20}")
        L.append("")
    return L


def table_fake(rows, dss, sat):
    L = ["Fake-cluster candidates: clusters whose seed (largest member) is flagged.  Population: all PbWO4-seeded clusters.  n; % of clusters [Wilson 68 %]; rate per 1000 events [Poisson-Wilson 68 %]"]
    for c in (0,) + CATS:
        L += ["", f"{CNAME[c]}:", f"{'dataset':<26}{'events':>10}{'n':>8}  {'% of clusters [68 % CI]':<30}{'rate / 1000 events [68 % CI]':<30}"]
        for ds in dss:
            r = get(rows, dataset=ds[0], population="all", axis="overall", category=CNAME[c], role="seed")
            L.append(f"{ds[1]:<26}{r['n_events']:>10d}{r['n_flag']:>8d}  {ci(r):<30}{rate(r):<30}")
    L += ["", "Flagged seed (any category) vs E_g / E_beam: n; rate per 1000 events [68 %]:"]
    for ds in dss:
        L.append(f"{ds[1]}:")
        for i in range(len(EBINS) - 1):
            r = get(rows, dataset=ds[0], population="all", axis="ebin", bin=f"{EBINS[i]:g}-{EBINS[i + 1]:g}", category="any", role="seed")
            L.append(f"   {EBINS[i]:g}-{EBINS[i + 1]:g}".ljust(14) + f"{r['n_flag']:>8d}  {rate(r)}")
    for c in CATS:
        L += ["", f"seed {CNAME[c]} vs E_g / E_beam: n (rate per 1000 events [68 %]):"]
        for ds in dss:
            L.append(f"{ds[1]:<26}" + "".join(
                f"{(r := get(rows, dataset=ds[0], population='all', axis='ebin', bin=f'{EBINS[i]:g}-{EBINS[i + 1]:g}', category=CNAME[c], role='seed'))['n_flag']:>5d} ({rate(r)})  "
                for i in range(len(EBINS) - 1)))
    L += ["", "Saturated seed (f4 or f5) in the ep elastic window: n / N, % [68 % CI]; of these, clusters with <= 2 member modules (n, % of saturated-seed clusters [68 % CI]), median member modules",
          "(saturated-seed clusters / all elastic-window clusters), median seed energy share of the saturated-seed clusters"]
    for ds in dss:
        k, n, k2, mc, ma, fs = sat[ds[0]]
        lo, hi = p4.wilson(k, n)
        l2, h2 = p4.wilson(k2, k)
        L.append(f"{ds[1]:<26}{k:>6d} / {n:<9d} {pct(k / max(n, 1))} [{pct(float(lo))}, {pct(float(hi))}]   <= 2 modules: {k2:>4d}" + (f" = {pct(k2 / k)} [{pct(float(l2))}, {pct(float(h2))}]" if k else "") + "   "
                 f"median modules {mc:.1f} / {ma:.1f}   seed share {fs:.3f}")
    return L


def cells(rows, dss):
    """(matter, negligible-max) per dataset: cells (population elastic / moller overall and all x E bin; category x role) with <fE> >= EF_LEVEL."""
    out = {}
    for ds in dss:
        hit, rest = [], {"ef_mean": 0.0, "ef_sem": 0.0, "n_clusters": 2}
        for r in rows:
            if r["table"] != "frac" or r["dataset"] != ds[0] or r["category"] == "any" or r["role"] == "any":
                continue
            if not (r["population"] in ("elastic", "moller") and r["axis"] in ("overall", "ebin")):
                continue
            if r["n_flag"] < MIN_CELL or not np.isfinite(r["ef_mean"]):
                continue
            if r["ef_mean"] >= EF_LEVEL:
                hit.append(r)
            else:
                rest = r if r["ef_mean"] > rest["ef_mean"] else rest
        out[ds[0]] = (hit, rest)
    return out


def meaning_lines(rows, dss, sat):
    C = cells(rows, dss)
    L = ["Criterion: a cell (population x category x role [x E bin]) matters for the calibration at the 0.1 % level if the mean energy fraction carried by its",
         "flagged modules, <fE> over all clusters of the cell, is >= 0.1 %; P8's per-channel energy error (at most ~100 % of the flagged energy) times <fE> bounds the",
         "shift of the mean cluster energy, so cells below 0.1 % cannot move the overall energy scale by 0.1 % whatever P8 finds (<fE> is a population",
         "average: a single hot module can exceed it, which P8 / P9 judge per module).  Cells with fewer than",
         f"{MIN_CELL} flagged clusters are not judged.  Populations judged: the calibration samples, elastic and moller (overall and per E_g / E_beam bin).",
         "Flagged seeds outside them (e.g. f4 / f5 seeds at E_g > 1.1 E_beam) are rate / fake-cluster effects: Table 3.", ""]
    for ds in dss:
        hit, rest = C[ds[0]]
        L.append(f"{ds[1]}:")
        if not hit:
            L.append(f"  no cell reaches 0.1 % (largest <fE> {efs(rest)} %): negligible for the calibration.")
        grp = {}
        for r in hit:
            grp.setdefault((r["population"], r["category"], r["role"]), []).append(r)
        for (pop, cat, role), rs in sorted(grp.items(), key=lambda kv: -max(r["ef_mean"] for r in kv[1])):
            where = ", ".join(f"{r['bin']} {efs(r)} % ({r['n_flag']} of {r['n_clusters']})" for r in sorted(rs, key=lambda r: r["bin"])
                              ) if pop == "all" else f"{efs(rs[0])} % ({rs[0]['n_flag']} of {rs[0]['n_clusters']} clusters)"
            L.append(f"  matters: {POP_LABEL[pop]}, {cat}, {role}: " + ("<fE> per E_g / E_beam bin: " if pop == "all" else "<fE> ") + where)
        if hit:
            L.append(f"  every other cell < 0.1 % (largest {efs(rest)} %).")
        k, n, k2 = sat[ds[0]][:3]
        lo, hi = p4.wilson(k, n)
        L.append(f"  saturated seed in the elastic window: {k} of {n} ({pct(k / max(n, 1))} % [{pct(float(lo))}, {pct(float(hi))}]), {k2} of them with <= 2 member modules.")
        L.append("")
    return L


def frac_axis(ax, rows, ds, pop, axis, cat, labels, style, lab_pop):
    x = np.arange(len(labels))
    for role, col, dx in (("seed", "tab:red", -0.12), ("neighbour", "tab:blue", 0.12)):
        R = [get(rows, dataset=ds[0], population=pop, axis=axis, bin=b, category=CNAME[cat], role=role) for b in labels]
        f = np.array([r["frac"] if r and r["n_clusters"] else np.nan for r in R]) * 100
        lo = np.array([r["frac_lo"] if r and r["n_clusters"] else np.nan for r in R]) * 100
        hi = np.array([r["frac_hi"] if r and r["n_clusters"] else np.nan for r in R]) * 100
        pos = f > 0
        ax.errorbar(x[pos] + dx + (0.05 if style == "open" else 0), f[pos], yerr=[f[pos] - lo[pos], hi[pos] - f[pos]], fmt="o" if style == "solid" else "s",
                    color=col, mfc=col if style == "solid" else "white", ms=5, capsize=2, lw=1, label=f"{role}, {lab_pop}")
        z = (f == 0) & np.isfinite(hi)
        ax.plot(x[z] + dx, hi[z], "v", color=col, mfc="none", ms=6, label=None)


def axis_pages(pdf, plt, rows, dss, axis):
    labels = ([f"{RADBINS[i]:g}-{RADBINS[i + 1]:g}" for i in range(len(RADBINS) - 1)] if axis == "ring"
              else [f"{EBINS[i]:g}-{EBINS[i + 1]:g}" for i in range(len(EBINS) - 1)])
    for cat in (0,) + CATS:
        fig, axs = plt.subplots(2, 2, figsize=(15, 10.5))
        for ax, ds in zip(axs.ravel(), dss[:4]):
            frac_axis(ax, rows, ds, "all", axis, cat, labels, "solid", "all clusters")
            if axis == "ring":
                frac_axis(ax, rows, ds, "elastic", axis, cat, labels, "open", "elastic window")
            ax.set_yscale("log")
            ax.set_xticks(np.arange(len(labels)))
            ax.set_xticklabels(labels, fontsize=9)
            ax.set_xlabel("seed module radius (mm)" if axis == "ring" else "E_g / E_beam")
            ax.set_ylabel("clusters with a flagged module (%)")
            ax.set_title(f"{ds[1]}")
            ax.grid(alpha=0.3, which="both")
            ax.legend(fontsize=8, loc="best")
        fig.suptitle(f"{CNAME[cat]}: fraction of clusters with a flagged seed / neighbour vs {'seed radius ring (P4 rings)' if axis == 'ring' else 'cluster energy'}\n"
                     "Wilson 68 %; open triangle = upper edge of the interval where no cluster is flagged", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        pdf.savefig(fig)
        plt.close(fig)


def ef_page(pdf, plt, D, X, dss, Ps):
    fig, axs = plt.subplots(2, 2, figsize=(15, 10.5))
    bins = np.logspace(-4, 0, 41)
    for ax, ds in zip(axs.ravel(), dss[:4]):
        P = Ps[ds[0]]
        for pop, col in (("all", "0.3"), ("hiE", "tab:blue"), ("elastic", "tab:red"), ("moller", "tab:green")):
            if pop not in P:
                continue
            f = X["f_all"][P[pop] & (X["f_all"] > 0)]
            if f.size:
                ax.hist(np.clip(f, 1e-4, 1), bins=bins, histtype="step", color=col, lw=1.5, label=f"{POP_LABEL[pop]} (n = {f.size})")
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlabel("flagged energy fraction of the cluster, fE (clusters with >= 1 flagged module)")
        ax.set_ylabel("clusters")
        ax.set_title(ds[1])
        ax.legend(fontsize=8)
        ax.grid(alpha=0.3)
    fig.suptitle("Fraction of the cluster energy carried by flagged modules (all categories)", fontsize=13)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def example_pages(pdf, plt, exs):
    from matplotlib.patches import Rectangle
    from matplotlib.colors import LogNorm
    TT = 4.0 * np.arange(100)
    for k0 in range(0, len(exs), 2):
        fig, axs = plt.subplots(2, 2, figsize=(15, 10.5), gridspec_kw=dict(width_ratios=(1.25, 1)))
        for row, ex in zip(axs, exs[k0:k0 + 2]):
            g = ex["geo"]
            nm = g["name"][ex["seed"]]
            w = ex["wave"]
            ax = row[0]
            if w is not None:
                ax.plot(TT, w["y"], "k.", ms=3.5, label="data - baseline")
                if w["k"]:
                    ax.plot(TT, np.sum(w["comp"], 0), "r-", lw=1.4, label="fit (sum)")
                    for j in range(w["k"]):
                        ax.plot(TT, w["comp"][j], "--", lw=1.1, label=f"pulse {j + 1}: A {w['a'][j]:.0f}, t {w['t'][j]:.0f} ns")
                ax.axhline(0, color="0.6", lw=0.8)
                ax.set_title(f"{RT_LABEL[ex['rt']]}, run {ex['run']} ev {ex['evn']}: seed {nm}, {CNAME[ex['scat']]}\nk = {w['k']}, chi2/dof = {w['chi2']:.1f}, nclip = {w['nclip']}", fontsize=9)
                ax.legend(fontsize=7)
            else:
                ax.text(0.5, 0.5, "waveform not re-read (EVIO file or event missing)", ha="center", transform=ax.transAxes)
            ax.set_xlabel("time (ns)")
            ax.set_ylabel("ADC")
            ax = row[1]
            x0, y0 = g["x"][ex["seed"]], g["y"][ex["seed"]]
            mem = dict(ex["members"])
            emax = max(list(ex["hits"].values()) + [1.0])
            norm = LogNorm(vmin=5.0, vmax=emax)
            cmap = plt.get_cmap("viridis")
            near = [m for m in range(len(g["x"])) if abs(g["x"][m] - x0) < 3.6 * g["sx"][ex["seed"]] and abs(g["y"][m] - y0) < 3.6 * g["sy"][ex["seed"]]]
            for m in near:
                e = ex["hits"].get(m, 0.0)
                fc = cmap(norm(max(e, 5.0))) if e > 0 else "#f0f0f0"
                cc = ex["mod_cat"].get(m, (0, 0))[0]
                ax.add_patch(Rectangle((g["x"][m] - g["sx"][m] / 2, g["y"][m] - g["sy"][m] / 2), g["sx"][m], g["sy"][m], fc=fc,
                                       ec=CCOL[cc] if cc else ("k" if m in mem else "0.75"), lw=3.0 if cc else (1.2 if m in mem else 0.4)))
                if e > 0:
                    ax.text(g["x"][m], g["y"][m], f"{e:.0f}" + (f"\nf{cc}" if cc else ""), ha="center", va="center", fontsize=7,
                            color="white" if norm(max(e, 5.0)) < 0.6 else "black")
            ax.plot(x0, y0, "r+", ms=14)
            ax.set_xlim(x0 - 3.6 * g["sx"][ex["seed"]], x0 + 3.6 * g["sx"][ex["seed"]])
            ax.set_ylim(y0 - 3.6 * g["sy"][ex["seed"]], y0 + 3.6 * g["sy"][ex["seed"]])
            ax.set_aspect("equal")
            esum = sum(mem.values())
            fl = sum(e for m, e in mem.items() if ex["mod_cat"].get(m, (0, 0))[0])
            ax.set_title(f"cluster E_g {esum:.0f} MeV, E_g / E' = {ex['r']:.3f}, flagged {100 * fl / max(esum, 1e-9):.1f} % of E_g", fontsize=9)
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
        for row in axs[len(exs[k0:k0 + 2]):]:
            for ax in row:
                ax.axis("off")
        fig.suptitle("Left: seed waveform, EVIO re-read and pileup_rates fit.  Right: 7 x 7 modules, gain-corrected benchmark hits (MeV); cluster member = black edge,\n"
                     "flagged module = coloured edge (f4 red, f5 orange, f1 blue, f2 green, f3 purple); red cross = seed centre", fontsize=10)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        pdf.savefig(fig)
        plt.close(fig)


def sub(d, idx):
    return {q: v[idx] for q, v in d.items()}


COLS = ("table", "dataset", "population", "axis", "bin", "category", "role", "n_clusters", "n_flag", "frac", "frac_lo", "frac_hi", "ef_mean",
        "ef_med_flag", "ef_q84_flag", "ef_q99_flag", "n_events", "per_kevent", "mu", "sigma", "n_fit", "own_fit", "source", "n_splits",
        "ef_sem", "rate_lo68", "rate_hi68")


def report(metas, D, out, n_ex, evio_dir, nproc):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with np.load(next((SCAN_ROOT / "exotic_sample").glob("run_*.npz"))) as z:
        names = [str(x) for x in z["mod_name"]]
    X = derive(metas, D, names)
    dss = datasets(metas)
    mol = moller_keys()
    rows, Ps, fits, msts, sat, exs = [], {}, {}, {}, {}, []
    rng = np.random.default_rng(20261118)
    for ds in dss:
        P, fit, mst = populations(metas, D, X, ds, mol)
        Ps[ds[0]], fits[ds[0]], msts[ds[0]] = P, fit, mst
        ms = [m for m in metas if m["rt"] == ds[2] and m["src"] == ds[3]]
        n_ev = sum(m["n_events"] for m in ms)
        idx = np.flatnonzero(P["all"])
        rows += frac_rows(ds[0], {p: v[idx] for p, v in P.items()}, sub(D, idx), sub(X, idx), n_ev)
        for lab, mu, s, n, own in fit:
            rows.append(dict(table="window", dataset=ds[0], population="single", axis="ring", bin=lab, mu=mu, sigma=s, n_fit=n, own_fit=int(own)))
        el = P["elastic"]
        ms_ = el & np.isin(D["scat"], (4, 5))
        sat[ds[0]] = (int(ms_.sum()), int(el.sum()), int((ms_ & (D["cnt"] <= 2)).sum()),
                      float(np.median(D["cnt"][ms_])) if ms_.any() else np.nan, float(np.median(D["cnt"][el])) if el.any() else np.nan,
                      float(np.median(X["f_seed"][ms_])) if ms_.any() else np.nan)
        k_, n_ = sat[ds[0]][:2]
        sl, sh = p4.wilson(k_, n_)
        rows.append(dict(table="sat", dataset=ds[0], population="elastic", axis="overall", bin="all", category="f4+f5", role="seed", n_clusters=n_, n_flag=k_,
                         frac=k_ / n_ if n_ else np.nan, frac_lo=float(sl) if n_ else np.nan, frac_hi=float(sh) if n_ else np.nan))
        rows.append(dict(table="scope", dataset=ds[0], source=ds[3], n_splits=len(ms), n_events=n_ev, n_clusters=int(P["all"].sum())))
        if ds[0] != "Ta_sample" and n_ex:
            exs += pick_examples(P, D, X, metas, ds, n_ex, rng)
    write_tsv(out / f"{STEM}.tsv", rows, COLS)
    # examples
    exd = []
    if exs:
        t0 = time.monotonic()
        with mp.get_context("fork").Pool(nproc) as pool:
            exd = pool.map(example_job, [(e, evio_dir) for e in exs], chunksize=1)
        print(f"[examples] {len(exd)} events, {time.monotonic() - t0:.0f} s", flush=True)
        write_tsv(out / f"{STEM}_examples.tsv",
                  [dict(e, module=e["geo"]["name"][e["seed"]], category=CNAME[e["scat"]], wave=int(e["wave"] is not None), k=e["wave"]["k"] if e["wave"] else "",
                        n_members=len(e["members"])) for e in exd],
                  ("dataset", "run", "split", "evn", "module", "category", "eg", "r", "ebeam", "k", "wave", "n_members"))
    # text
    cov = []
    for ds in dss:
        ms = [m for m in metas if m["rt"] == ds[2] and m["src"] == ds[3]]
        nf = {c: sum(m["n_rows_flag_cat"][c] for m in ms) for c in CATS}
        nc = {c: sum(m["n_rows_in_cluster_cat"][c] for m in ms) for c in CATS}
        cov.append(f"  {ds[1]:<26} splits {len(ms):>3}, events {sum(m['n_events'] for m in ms):>9}, clusters {sum(m['n_clusters'] for m in ms):>9}, "
                   f"LG-seeded dropped {sum(m['n_lg_seed'] for m in ms)}, centre != largest member {sum(m['n_seed_not_max'] for m in ms)}, "
                   f"repeated-module hits {sum(m['n_dup_hits'] for m in ms)}, "
                   f"exotic events missing from the benchmark scan {sum({m['fn']: m['n_missing_file'] for m in ms}.values())}")
        cov.append("      flagged rows in a cluster / flagged rows: " + ", ".join(f"f{c} {nc[c]}/{nf[c]}" for c in CATS))
    win = []
    for ds in dss:
        win.append(f"  {ds[1]:<26}" + "  ".join(f"{lab}: {mu:.4f}+-{s:.4f} (n {n}{'' if own else ', run-type fit'})" for lab, mu, s, n, own in fits[ds[0]][:4]))
        win.append(f"  {'':<26}" + "  ".join(f"{lab}: {mu:.4f}+-{s:.4f} (n {n}{'' if own else ', run-type fit'})" for lab, mu, s, n, own in fits[ds[0]][4:]))
    mlines = []
    for ds in dss:
        if msts[ds[0]]:
            m = msts[ds[0]]
            mlines.append(f"  {ds[1]:<26} selected Moller clusters matched {m['n_match']} in {m['splits']} splits; max |E_P7 / E_moller_cache - 1| = {m['de_max']:.2e}")
    head = [f"Join: exotic-scan rows (pileup_rates, PILEUP_EXOTIC=1) -> benchmark-scan clusters by (run, EVIO event number, module); clusters = erb.clusterers' gated",
            f"HyCalCluster at CONTROL_GATE ({', '.join(f'{k} {v:g} ns' for k, v in erb_mod().CONTROL_GATE.items() if k in TYPES)}) on the benchmark hits with the M1 phase factor",
            "and time offsets (control v3), member energies x control-v3 per-run gains.  Exclusive categories f4 > f5 > f1 > f2 > f3, uniform P2 thresholds, PbWO4 only.",
            "Seed role = largest member (gain-corrected); neighbour = other member.  Populations: all PbWO4-seeded clusters; hiE: E_g > 0.5 E_beam; elastic: single-cluster events (one cluster",
            "> 500 kE MeV), seed >= 200 kE MeV, |E_g / E' - mu| < 3 sigma per P4 ring; moller: selected pairs of moller_calib.select (Ta only, gate 4 ns, 60 mm isolation).",
            "Sources: H2 = exotic_sample (2 splits per run; 0.73: 24590, 24603; 2.24: 25228, 25257, 25260, 25299, 25320; 3.49: 24698, 24911); Ta X17 main = the 19",
            "exotic_full npz (25800 s100-s600, 26230 s100-s1300); Ta sample (25800, 26138, 26230, 2 splits each) = cross-check (4 of its 6 splits are also in the full set).",
            f"Run 24937 excluded.  0.73 GeV elastic window: only angles where ep and Moller separate (midpoint of E_moller / E' and 1 < {MOLLER_SEP_MAX:g},",
            "seed radius >~ 110 mm; the inner merged region is left out) and E_g / E' above that midpoint, before the fit.", "",
            "Data scope and join coverage:"] + cov + ["", "Elastic window fits (E_g / E', mu +- sigma per P4 ring; window = mu +- 3 sigma):"] + win + (
            ["", "Moller sample:"] + mlines if mlines else [])
    pops = POPS
    with PdfPages(out / f"{STEM}.pdf") as pdf:
        text_page(pdf, plt, "P7 exotic categories at cluster level: definitions, data, join", head, size=8.6)
        text_page(pdf, plt, "Table 1: clusters with >= 1 flagged module, per run type and population",
                  ["fractions in % (Wilson 68 %); f<c> % = clusters with >= 1 module of category c (any role); <fE> = mean over all clusters of the flagged energy",
                   "fraction; fE|fl = median / 84 % / 99 % of fE over the clusters with >= 1 flagged module.  Source per row: dataset label (H2 sample, Ta full or sample)", ""]
                  + table_main(rows, dss, pops), size=8.2)
        for pop in ("all", "elastic", "moller"):
            text_page(pdf, plt, f"Table 2 ({pop}): category x role", table_roles(rows, [d for d in dss if pop != "moller" or d[2] == "2p24GeV_Ta_X17"], pop), size=8.4)
        text_page(pdf, plt, "Table 3: fake-cluster candidates (flagged seed) and saturated seeds", table_fake(rows, dss, sat), size=7.6)
        text_page(pdf, plt, "What this means for P8 / P9", meaning_lines(rows, dss, sat), size=8.6)
        axis_pages(pdf, plt, rows, dss, "ring")
        axis_pages(pdf, plt, rows, dss, "ebin")
        ef_page(pdf, plt, D, X, dss, Ps)
        if exd:
            text_page(pdf, plt, "Examples: flagged seed in the elastic window", [f"{n_ex} per run type (Ta: full scan), round robin over the seed categories present, random within a category (seed 20261118).",
                      "Left: seed waveform re-read from EVIO and refitted by pileup_rates (stored noise reference of the exotic scan); right: 7 x 7 module",
                      "neighbourhood, gain-corrected benchmark hit energies (MeV), cluster members black, flagged modules coloured by category.", ""]
                      + [f"  {RT_LABEL[e['rt']]:<18} run {e['run']} split {e['split']:>5} ev {e['evn']:>10}  seed {e['geo']['name'][e['seed']]:<6} {CNAME[e['scat']]:<20} E_g {e['eg']:>7.0f} MeV  E_g/E' {e['r']:.3f}"
                         for e in exd], size=8.4)
            example_pages(pdf, plt, exd)
    for L in (table_main(rows, dss, pops), table_fake(rows, dss, sat), meaning_lines(rows, dss, sat), head):
        print("\n".join(L))
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--cache", default=str(SCAN_ROOT / "exotic_p7_cache.npz"), help="pass-1 cluster records (reused if present; delete to rebuild)")
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--examples", type=int, default=10, help="example events per run type (0 = none)")
    ap.add_argument("--evio-dir", default=p1.EVIO_DIR)
    ap.add_argument("--only", default="", help="debug: substring of the exotic file names to process (no cache written)")
    a = ap.parse_args()
    jobs = jobs_of()
    if a.only:
        jobs = [j for j in jobs if a.only in Path(j[1]).name]
        a.cache = ""
    metas, D = collect(jobs, a.jobs, a.cache)
    report(metas, D, a.out, a.examples, a.evio_dir, a.jobs)


if __name__ == "__main__":
    main()
