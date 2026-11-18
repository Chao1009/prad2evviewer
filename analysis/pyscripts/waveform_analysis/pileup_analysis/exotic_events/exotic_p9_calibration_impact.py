#!/usr/bin/env python3
"""exotic_p9_calibration_impact.py -- P9 of the exotic-events impact packages: what vetoing or correcting the flagged clusters does to the calibration and the resolution.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p9_calibration_impact.py [--stage hooks|run|analyse|all] [--jobs 2] [--root <dir>] [--out <dir>]
                                                                                            [--nboot 200]

Arms (identical events per run type; PbWO4 only; run 24937 excluded):
  A   default control v3 (energy_resolution_benchmark.py report, ep calibration; Ta X17 also the Moller two-cluster calibration, moller_calib.py)
  V3  veto: clusters with a saturated seed (f4 or f5) and seed energy fraction > 0.95 dropped (P10 V3, the policy on 2.24 / 3.49 GeV H2)
  V4  veto variant: V3 + clusters with an f1 or f3 seed (P10 V4)
  C   corrected: each flagged member's raw energy x (1 + dE / E_def), dE the credible P8 alternative of its category (f1 sloped baseline everywhere, f3 stretched single
      in Ta only, f5 k = 1 extrapolation for mild clipping, nclip <= NCLIP_MILD overflow samples; f2, f4, deep f5 and H2 f3 keep the default); first order: clustering,
      seed, position, pool selection and the Moller candidate list are not redone.
Data: H2 run types on the splits of the 2-file exotic sample (output/pileup_rates/exotic_sample), Ta X17 on the 19 full-scan splits (output/pileup_rates/exotic_full;
25800 x 6, 26230 x 13), as P7 / P8.  The vetoes are P7's clusters (output/pileup_rates/exotic_p7_cache.npz: control gate, keyed by run, EVIO event number, seed module),
the corrections P8's records (output/pileup_rates/exotic_p8_cache.npz: keyed by run, event number, seed, member module).  0.73 GeV: every arm drops clusters seeded
below RMIN mm (user decision 2026-11-18: calibration sample at R >~ 107 mm, where the ep and Moller peaks separate).

Machinery.  Stage hooks writes one npz per run type and arm (the benchmark's opt-in EXOTIC_HOOK: splits, veto keys, rmin, member corrections); stage run calls
`energy_resolution_benchmark.py report --set EXOTIC_HOOK=... --out <root>/<arm>/...` per run type (ep: gains per run in and out of sample, position correction, ring
calibration, fits per split, exactly the control-v3 chain on these splits) and `report --calib moller` for Ta X17 (candidate cache under the arm's root, the
2.24 GeV H2 position correction borrowed from the control-v3 root for both Ta calibrations).  Stage analyse reads the reports (summary.tsv / splits.tsv, gains.tsv)
and the per-cluster dumps (p9_clusters.npz, p9_moller.npz): arm values are the benchmark's own (mean +- standard error over the splits); differences B - A and C - A
are paired: the same clusters (shared keys, used in both arms) refitted pooled with the benchmark's fit_path (ep: two-peak Crystal Ball, tails beyond -+3 sigma) or
moller_calib.cb_fit (Moller), error from NBOOT bootstrap resamples of the shared clusters (Moller: of the pairs); per P4 ring of the seed radius likewise.  Gains:
per run and module from each arm's gains.tsv (ep: pool solve; Moller: the per-run in-sample solve), dg = g_arm - g_A against sigma_g of A.

Outputs (--out): exotic_P9_calibration_impact_20261118.pdf and .tsv.  Intermediates under --root (default output/pileup_rates/exotic_p9: hooks, the arm reports, logs).
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import csv
import math
import shutil
import zlib
import argparse
import subprocess
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
from exotic_p1_quantities import RT_LABEL, EXCLUDE_RUNS # noqa: E402
from exotic_p3_sample import text_page                  # noqa: E402

STEM = "exotic_P9_calibration_impact_20261118"
RTS = ("0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17")
TA = "2p24GeV_Ta_X17"
ARMS = ("A", "V3", "V4", "C")
CMP = ("V3", "V4", "C")                                 # compared with A
ARM_LABEL = {"A": "A default", "V3": "B veto V3", "V4": "B' veto V4", "C": "C corrected"}
RMIN = {"0p73GeV_H2": 107.0}                            # mm, seed-centre radius cut of every arm (user decision 2026-11-18)
NCLIP_MILD = 3                                          # f5: k = 1 extrapolation credible up to this many overflow samples (P8 JUDGE: "mild clipping", W398 = 1)
NONSHOWER_F = 0.95                                      # V3: seed energy fraction above this (P10)
RADBINS = p4.RADBINS                                    # P4 rings (mm, seed module centre)
NBOOT = 200
MIN_RING = 300                                          # shared clusters a ring needs for a fit
TA_MODS = ("W564", "W559", "W493", "W663", "W593")      # Ta inner ring (P8)
H2_MODS = ("W398", "W457", "W613")                      # W398 high-gain (clipped real showers), W457 fakes, W613 f4
NAMED = TA_MODS + H2_MODS
PYDIR = p1.PYDIR
WA = PYDIR / "waveform_analysis"
SCAN_ROOT = PYDIR / "output" / "pileup_rates"
BENCH_SCAN = PYDIR / "output" / "energy_resolution" / "scan"
V3_ROOT = PYDIR / "output" / "energy_resolution" / "control_phase2_v3_20261114" / "v3_defaults"
P7_CACHE = SCAN_ROOT / "exotic_p7_cache.npz"
P8_CACHE = SCAN_ROOT / "exotic_p8_cache.npz"
P10_TSV = HERE / "exotic_P10_yield_impact_20261118.tsv"
BASELINE = {"0p73GeV_H2": "3.160", "2p24GeV_H2": "1.746", "3p49GeV_H2": "1.470", TA: "1.987 (ep), Moller pair sum 2.09 +- 0.04"}   # full-data control v3, oos CB sigma/E %


def erb_mod():
    import energy_resolution_benchmark as erb
    return erb


# ---------------------------------------------------------------------------------------------------------------- stage hooks
def splits_of():
    """{run type: [(run, split)]}: H2 from the exotic sample, Ta X17 from the full-scan npz (P7 / P8 scope)."""
    out = {rt: [] for rt in RTS}
    for f in sorted((SCAN_ROOT / "exotic_sample").glob("run_*.npz")):
        with np.load(f) as z:
            run, rt = int(z["run"]), str(z["run_type"])
            if run in EXCLUDE_RUNS or rt == TA:
                continue
            out[rt] += [(run, int(s)) for s in z["splits"]]
    for f in sorted((SCAN_ROOT / "exotic_full").glob("run_*_s*.npz")):
        with np.load(f) as z:
            if str(z["run_type"]) == TA and int(z["run"]) not in EXCLUDE_RUNS:
                out[TA] += [(int(z["run"]), int(s)) for s in z["splits"]]
    return {rt: sorted(v) for rt, v in out.items()}


def in_scope(rt, src):
    return (rt == TA and src == "full") or (rt != TA and src == "sample")


def veto_keys():
    """{run type: {"V3": [n, 3], "V4": [n, 3]}} (run, event number, seed) of P7's clusters (control gate, scope sources)."""
    c = np.load(P7_CACHE, allow_pickle=True)
    metas = list(c["metas"])
    J = c["d_job"]
    scat, es, eg = c["d_scat"].astype(np.int64), c["d_es"].astype(np.float64), c["d_eg"].astype(np.float64)
    sat = (scat == 4) | (scat == 5)
    v3 = sat & (es / np.maximum(eg, 1e-9) > NONSHOWER_F)
    v4 = v3 | (scat == 1) | (scat == 3)
    rt_j = np.array([m["rt"] for m in metas])
    ok_j = np.array([in_scope(m["rt"], m["src"]) for m in metas])
    run_j = np.array([m["run"] for m in metas], np.int64)
    evn, seed = c["d_evn"].astype(np.int64), c["d_seed"].astype(np.int64)
    out = {}
    for rt in RTS:
        m = ok_j[J] & (rt_j[J] == rt)
        out[rt] = {}
        for name, v in (("V3", v3), ("V4", v4)):
            i = np.flatnonzero(m & v)
            out[rt][name] = np.unique(np.stack([run_j[J[i]], evn[i], seed[i]], 1), axis=0) if i.size else np.zeros((0, 3), np.int64)
    return out


def credible_de(rt, cat, de, nclip):
    """dE (MeV, gain-corrected) of the credible P8 alternative per record; 0 where the default stays."""
    out = np.zeros(cat.size)
    f1 = cat == 1
    out[f1] = de["bs"][f1]
    if rt == TA:
        out[cat == 3] = de["st"][cat == 3]
    mild = (cat == 5) & (nclip <= NCLIP_MILD)
    out[mild] = de["k1"][mild]
    return np.nan_to_num(out, nan=0.0)


def corrections():
    """{run type: (corr [n, 4] (run, event number, seed, module), rho [n], counts)} from the P8 records."""
    c = np.load(P8_CACHE, allow_pickle=True)
    metas = list(c["metas"])
    J = c["r_job"].astype(np.int64)
    rt_j = np.array([m["rt"] for m in metas])
    ok_j = np.array([in_scope(m["rt"], m["src"]) for m in metas])
    run_j = np.array([m["run"] for m in metas], np.int64)
    ci = c["r_ci"].astype(np.int64)
    cat, nclip, edef = c["r_cat"].astype(np.int64), c["r_nclip"].astype(np.int64), c["r_e_def"].astype(np.float64)
    de = {v: c[f"r_de_{v}"].astype(np.float64) for v in ("bs", "st", "k1")}
    out = {}
    for rt in RTS:
        m = ok_j[J] & (rt_j[J] == rt)
        i = np.flatnonzero(m)
        d = credible_de(rt, cat[i], {k: v[i] for k, v in de.items()}, nclip[i])
        rho = np.where(edef[i] > 0, d / np.maximum(edef[i], 1e-9), 0.0)
        key = np.stack([run_j[J[i]], c["r_evn"][i].astype(np.int64), c["c_c_seed"][ci[i]].astype(np.int64), c["r_midx"][i].astype(np.int64)], 1)
        _, inv, cnt = np.unique(key, axis=0, return_inverse=True, return_counts=True)
        amb = cnt[inv.ravel()] > 1                      # two time-split clusters of one event with the same seed: the hook key cannot tell them apart -> default
        nz = np.flatnonzero((rho != 0) & ~amb)
        counts = {f"f{k}": int(((cat[i] == k) & (rho != 0)).sum()) for k in (1, 3, 5)} | dict(records=int(i.size), corrected=int(nz.size),
                                                                                           f5_deep=int(((cat[i] == 5) & (nclip[i] > NCLIP_MILD)).sum()),
                                                                                           ambiguous=int(((rho != 0) & amb).sum()))
        out[rt] = (key[nz], rho[nz], counts)
    return out


def hook_path(root, rt, arm):
    return Path(root) / "hooks" / f"{rt}_{arm}.npz"


def stage_hooks(root):
    sp = splits_of()
    vk = veto_keys()
    co = corrections()
    (Path(root) / "hooks").mkdir(parents=True, exist_ok=True)
    rows = []
    for rt in RTS:
        empty = np.zeros((0, 3), np.int64)
        for arm in ARMS:
            veto = vk[rt][arm] if arm in ("V3", "V4") else empty
            corr, rho = (co[rt][0], co[rt][1]) if arm == "C" else (np.zeros((0, 4), np.int64), np.zeros(0))
            np.savez(hook_path(root, rt, arm), splits=np.array(sp[rt], np.int64).reshape(-1, 2), veto=veto, rmin=np.float64(RMIN.get(rt, 0.0)),
                     corr=corr, corr_rho=rho)
            rows.append(dict(table="hook", run_type=rt, arm=arm, quantity="splits", value=len(sp[rt])))
            rows.append(dict(table="hook", run_type=rt, arm=arm, quantity="veto_keys", value=len(veto)))
            rows.append(dict(table="hook", run_type=rt, arm=arm, quantity="corr_rows", value=len(corr)))
        for k, v in co[rt][2].items():
            rows.append(dict(table="hook", run_type=rt, arm="C", quantity=f"p8_{k}", value=v))
        print(f"[hooks] {rt}: {len(sp[rt])} splits, V3 {len(vk[rt]['V3'])} / V4 {len(vk[rt]['V4'])} veto keys, C {co[rt][2]}", flush=True)
    write_rows(Path(root) / "hooks" / "hooks.tsv", rows)


# ---------------------------------------------------------------------------------------------------------------- stage run
def arm_root(root, rt, arm):
    """Report root of an arm: H2 run types share <root>/<arm>/ep; Ta X17 gets its own (<root>/<arm>/ta) holding the borrowed 2.24 GeV H2 correction."""
    return Path(root) / arm / ("ta" if rt == TA else "ep")


def run_cmds(root):
    py = sys.executable
    bench = str(WA / "energy_resolution_benchmark.py")
    cmds = []
    for arm in ARMS:
        for rt in RTS:
            hk = hook_path(root, rt, arm)
            out = arm_root(root, rt, arm)
            done = out / rt / "p9_clusters.npz"
            cmds.append((f"{rt}_{arm}", [py, bench, "report", "--set", f"EXOTIC_HOOK={str(hk)!r}", "--out", str(out), str(BENCH_SCAN / rt)], done, out))
        hk = hook_path(root, TA, arm)
        out = Path(root) / arm / "moller"
        cmds.append((f"moller_{arm}", [py, bench, "report", "--calib", "moller", "--set", f"EXOTIC_HOOK={str(hk)!r}", "--out", str(out), "--cache",
                                       str(out / "cache"), str(BENCH_SCAN / TA)], out / f"{TA}_moller" / "p9_moller.npz", out))
    return cmds


def run_one(c):
    name, cmd, done, out = c
    log = Path(c[3]).parents[1] / "logs" / f"{name}.log"
    log.parent.mkdir(parents=True, exist_ok=True)
    if done.exists():
        return name, 0, "exists"
    out.mkdir(parents=True, exist_ok=True)
    if name.startswith(TA):                             # the ep chain of Ta borrows the control-v3 2.24 GeV H2 position correction
        src, dst = V3_ROOT / "2p24GeV_H2", out / "2p24GeV_H2"
        dst.mkdir(parents=True, exist_ok=True)
        for f in ("poscorr_modules.tsv", "poscorr_map.tsv", "poscorr_extra.tsv"):
            if (src / f).exists():
                shutil.copy2(src / f, dst / f)
    with open(log, "w") as fh:
        r = subprocess.run(cmd, stdout=fh, stderr=subprocess.STDOUT, cwd=str(PYDIR))
    return name, r.returncode, str(log)


def stage_run(root, jobs):
    cmds = [c for c in run_cmds(root) if not c[2].exists()]
    cmds.sort(key=lambda c: (0 if (TA in c[0] or "moller" in c[0]) else 1))    # long Ta jobs first
    print(f"[run] {len(cmds)} reports, {jobs} at a time", flush=True)
    with mp.get_context("spawn").Pool(jobs) as pool:
        for name, rc, info in pool.imap_unordered(run_one, cmds):
            print(f"[run] {name}: rc {rc} ({info})", flush=True)
            if rc:
                raise SystemExit(f"[ERROR] {name} failed, see {info}")


# ---------------------------------------------------------------------------------------------------------------- stage analyse
def write_rows(path, rows):
    cols = []
    for r in rows:
        cols += [k for k in r if k not in cols]
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")


def fmt(v):
    if isinstance(v, (float, np.floating)):
        return "nan" if not math.isfinite(v) else f"{v:.6g}"
    return str(v)


def tsv(path):
    with open(path) as fh:
        return list(csv.DictReader((x for x in fh if not x.startswith("#")), delimiter="\t"))


def fnum(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")



# ---- ep: arm values, paired pooled differences
def gate_path(rt):
    return f"g{erb_mod().CONTROL_GATE[rt]:g}"


def ep_dir(root, rt, arm):
    return arm_root(root, rt, arm) / rt


def load_dump(root, rt, arm):
    z = np.load(ep_dir(root, rt, arm) / "p9_clusters.npz", allow_pickle=False)
    return {q: z[q] for q in z.files}


def keys_of(d):
    return np.stack([d["run"], d["split"], d["evn"], d["seed"]], 1).astype(np.int64)


def join(ka, kb):
    """Row indices (ia, ib) of the keys present in both [n, k] key arrays (each key unique within an array)."""
    a = {tuple(k): i for i, k in enumerate(ka.tolist())}
    ia, ib = [], []
    for j, k in enumerate(kb.tolist()):
        i = a.get(tuple(k))
        if i is not None:
            ia.append(i)
            ib.append(j)
    return np.array(ia, np.int64), np.array(ib, np.int64)


_G = {}                                                 # fork-shared data of the bootstrap workers


def ep_fit(rt, d, idx, e):
    """Pooled fit_path of the clusters idx of dump d with energies e[idx] (benchmark two-peak CB, tails with the run type's cut)."""
    erb = erb_mod()
    hz, eb = float(d["hz"]), float(d["ebeam"])
    seed = d["seed"][idx].astype(np.int64)
    us = np.unique(seed)
    ep = erb.e_elastic(erb.theta_deg(d["mod_x"][us], d["mod_y"][us], hz), eb, erb.pool_mass(rt))
    used = {int(m): float(v) for m, v in zip(us, ep)}
    kE, kS = round(eb / erb.E_REF, 6), round(math.sqrt(eb / erb.E_REF), 6)
    pool = dict(seed=seed, e=e[idx], x=d["x"][idx], y=d["y"][idx], pu=d["pu"][idx], clip=d["clip"][idx])
    geo = dict(hz=hz, ebeam=eb, x=d["mod_x"])
    ref = float(np.nanmean(d["_ref"][idx]))
    with np.errstate(all="ignore"):
        r, _ = erb.fit_path(pool, used, float(np.mean(ep)), ref, geo, kE, kS, False, erb.SEP_CUT_MEV[rt])
    return {q: float(r.get(q, np.nan)) for q in ("res_cb", "mu_ep", "tail_lo3", "tail_hi3")} | dict(n=int(idx.size))


EPQ = ("res_cb", "mu_ep", "tail_lo3", "tail_hi3")


def ep_boot_job(a):
    """Paired difference X - A on the shared clusters (sub = mask of them, e.g. a ring): point values and bootstrap sd."""
    rt, cmpn, prot, axis, b, nboot = a
    S = _G["ep"][rt, cmpn, prot]
    m = np.ones(S["iA"].size, bool) if axis == "overall" else S["ring"] == b
    iA, iX = S["iA"][m], S["iX"][m]
    if iA.size < MIN_RING:
        return a, None
    dA, dX, eA, eX = S["dA"], S["dX"], S["eA"], S["eX"]
    fA, fX = ep_fit(rt, dA, iA, eA), ep_fit(rt, dX, iX, eX)
    rng = np.random.default_rng(zlib.crc32(repr((rt, cmpn, prot, axis, b)).encode()))
    bs = []
    for _ in range(nboot):
        k = rng.integers(0, iA.size, iA.size)
        ga, gx = ep_fit(rt, dA, iA[k], eA), ep_fit(rt, dX, iX[k], eX)
        bs.append([gx[q] - ga[q] for q in EPQ])
    bs = np.array(bs)
    sd = {q: float(np.nanstd(bs[:, j], ddof=1)) for j, q in enumerate(EPQ)}
    return a, dict(A=fA, X=fX, sd=sd, nboot=int(np.isfinite(bs[:, 0]).sum()))


def ring_of(d, idx):
    return np.clip(np.digitize(np.hypot(d["mod_x"], d["mod_y"])[d["seed"][idx].astype(np.int64)], RADBINS) - 1, 0, len(RADBINS) - 2)


def prep_ep(root):
    """_G["ep"][rt, cmp, prot]: shared used clusters of A and the compared arm (in sample: e, used; out of sample: e_oos, used_oos)."""
    out = {}
    for rt in RTS:
        D = {arm: load_dump(root, rt, arm) for arm in ARMS}
        for arm, d in D.items():
            K = keys_of(d)
            if np.unique(K, axis=0).shape[0] != K.shape[0]:
                print(f"[WARN] {rt} {arm}: duplicate cluster keys in p9_clusters.npz", flush=True)
        for prot, eq, uq, rq in (("in", "e", "used", "ref"), ("oos", "e_oos", "used_oos", "ref_oos")):
            for c in CMP:
                dA, dX = D["A"], D[c]
                ua, ux = np.flatnonzero(dA[uq]), np.flatnonzero(dX[uq])
                ia, ix = join(keys_of(dA)[ua], keys_of(dX)[ux])
                iA, iX = ua[ia], ux[ix]
                dA2, dX2 = dict(dA, _ref=dA[rq]), dict(dX, _ref=dX[rq])
                out[rt, c, prot] = dict(dA=dA2, dX=dX2, eA=dA[eq], eX=dX[eq], iA=iA, iX=iX, ring=ring_of(dA, iA), nA=int(ua.size), nX=int(ux.size))
    return out


# ---- Moller (Ta X17)
def load_moller(root, arm):
    z = np.load(Path(root) / arm / "moller" / f"{TA}_moller" / "p9_moller.npz", allow_pickle=False)
    return {q: z[q] for q in z.files}


MQ = ("res", "peak", "tail_lo3", "tail_hi3")


def m_fit(x, T):
    import moller_calib as mc
    with np.errstate(all="ignore"):
        f = mc.cb_fit(x, T)
    mu, sg = f["mu"], f["sig"]
    lo, hi = mc.fit_scale(T)[3:5]
    inr = x[(x >= lo) & (x <= hi)]
    ok = math.isfinite(mu) and math.isfinite(sg) and inr.size
    return dict(res=f["res"], peak=100.0 * (mu / T - 1.0), tail_lo3=100.0 * float(np.mean(inr < mu - 3 * sg)) if ok else np.nan,
                tail_hi3=100.0 * float(np.mean(inr > mu + 3 * sg)) if ok else np.nan, n=int(x.size))


def m_values(d, eq, pairs, kind, ring=None):
    """x and T of the pairs (indices) for kind 'sum' (E1 + E2 vs Ts) or 'cluster' (E / E_pred x T, optionally one ring of the seed radius)."""
    e = d[eq]
    if kind == "sum":
        return e[2 * pairs] + e[2 * pairs + 1], float(d["Ts"])
    cl = np.r_[2 * pairs, 2 * pairs + 1]
    if ring is not None:
        rg = np.clip(np.digitize(np.hypot(d["x"], d["y"])[cl], RADBINS) - 1, 0, len(RADBINS) - 2)   # cluster position radius (seed centre not stored)
        cl = cl[rg == ring]
    return e[cl] / d["t"][cl] * float(d["T"]), float(d["T"])


def m_boot_job(a):
    cmpn, prot, kind, ring, nboot = a
    S = _G["m"][cmpn]
    eq = "e_in" if prot == "in" else "e_oos"
    pA, pX = S["pA"], S["pX"]
    xa, T = m_values(S["dA"], eq, pA, kind, ring)
    if xa.size < MIN_RING:
        return a, None
    fA, fX = m_fit(xa, T), m_fit(m_values(S["dX"], eq, pX, kind, ring)[0], T)
    rng = np.random.default_rng(zlib.crc32(repr(a[:4]).encode()))
    bs = []
    for _ in range(nboot):
        k = rng.integers(0, pA.size, pA.size)
        ga = m_fit(m_values(S["dA"], eq, pA[k], kind, ring)[0], T)
        gx = m_fit(m_values(S["dX"], eq, pX[k], kind, ring)[0], T)
        bs.append([gx[q] - ga[q] for q in MQ])
    bs = np.array(bs)
    return a, dict(A=fA, X=fX, sd={q: float(np.nanstd(bs[:, j], ddof=1)) for j, q in enumerate(MQ)}, nboot=int(np.isfinite(bs[:, 0]).sum()))


def prep_moller(root):
    D = {arm: load_moller(root, arm) for arm in ARMS}
    out = {}
    pk = lambda d: np.stack([d["run"][0::2], d["split"][0::2], d["evn"][0::2]], 1).astype(np.int64)   # noqa: E731
    for c in CMP:
        ia, ix = join(pk(D["A"]), pk(D[c]))
        out[c] = dict(dA=D["A"], dX=D[c], pA=ia, pX=ix, nA=D["A"]["run"].size // 2, nX=D[c]["run"].size // 2)
    return out


# ---- gains
def module_info(root):
    """Module names, x, y, ring, region code (benchmark module_regions) from the A dump of each run type and the clusterer geometry."""
    erb = erb_mod()
    d = load_dump(root, "2p24GeV_H2", "A")
    fn = sorted((BENCH_SCAN / "2p24GeV_H2").glob("run_*.npz"))[0]
    with np.load(fn, allow_pickle=False) as z:
        run, calib = int(z["run"]), str(z["calib"])
    _, geo = erb.clusterers(run, calib, (erb.CONTROL_GATE["2p24GeV_H2"],))
    assert list(geo["name"]) == list(d["mod_name"]), "module order differs between dump and geometry"
    rad = np.hypot(geo["x"], geo["y"])
    return dict(name=np.asarray(geo["name"]).astype(str), x=geo["x"], y=geo["y"], rad=rad, reg=erb.module_regions(geo), pwo4=geo["pwo4"],
                ring=np.clip(np.digitize(rad, RADBINS) - 1, 0, len(RADBINS) - 2), pool=erb.pool_modules(geo))


def gains_of(root, rt, arm, calib):
    """{(run, module): (g, sigma_g, n_seed, prior_frac, table)}; ep: pool solve (gains.tsv) and the ring-scoped solve (ring/gains.tsv) for modules not free in the pool
    solve; Moller: the per-run in-sample solves."""
    out = {}
    if calib == "moller":
        for r in tsv(Path(root) / arm / "moller" / f"{TA}_moller" / "gains.tsv"):
            run, k = r["solve"].split(":")
            if k == "in":
                out[int(run), r["module"]] = (fnum(r["g"]), fnum(r["sigma_g"]), fnum(r["n_seed"]), fnum(r["prior_frac"]), "moller")
        return out
    base = ep_dir(root, rt, arm)
    for tab, f in (("pool", base / "gains.tsv"), ("ring", base / "ring" / "gains.tsv")):
        if not f.exists():
            continue
        for r in tsv(f):
            if r["moller"] != "0" or r["nonlinear"] != "0":
                continue
            k = (int(r["run"]), r["module"])
            v = (fnum(r["g"]), fnum(r["sigma_g"]), fnum(r["n_seed"]), fnum(r["prior_frac"]), tab)
            if tab == "pool" or k not in out or out[k][3] >= 0.5 > v[3]:
                out[k] = v if tab == "pool" or k not in out else (v if out[k][3] >= 0.5 > v[3] else out[k])
    return out


def gain_rows(root, M):
    """Per (calibration, comparison, run, module): dg = g_X - g_A, sigma_g of A; free = prior fraction of A < 0.5 (the data set the gain)."""
    idx = {n: i for i, n in enumerate(M["name"])}
    rows = []
    for rt, calib in [(rt, "ep") for rt in RTS] + [(TA, "moller")]:
        G = {arm: gains_of(root, rt, arm, calib) for arm in ARMS}
        for c in CMP:
            for k, (g, sg, ns, pf, tab) in G["A"].items():
                if k not in G[c]:
                    continue
                i = idx.get(k[1])
                if i is None or not M["pwo4"][i]:
                    continue
                gx = G[c][k][0]
                rows.append(dict(table="gain", run_type=rt, calib=calib, comparison=f"{c}-A", run=k[0], module=k[1], ring=f"{RADBINS[M['ring'][i]]:g}-{RADBINS[M['ring'][i] + 1]:g}",
                                 pool=int(M["pool"][i]), solve=tab, g_A=g, sigma_g=sg, n_seed=ns, prior_frac=pf, free=int(pf < 0.5), dg=gx - g,
                                 z=(gx - g) / sg if sg > 0 else np.nan))
    return rows


# ---- benchmark's own arm values
def bench_rows(root):
    """Arm values as the benchmark quotes them (mean +- SE over the splits, control-gate path) and per-split paired differences."""
    rows = []
    for rt in RTS:
        gp = gate_path(rt)
        S, P = {}, {}
        for arm in ARMS:
            s = next(r for r in tsv(ep_dir(root, rt, arm) / "summary.tsv") if r["path"] == gp)
            S[arm] = s
            P[arm] = {(r["run"], r["split"]): r for r in tsv(ep_dir(root, rt, arm) / "splits.tsv") if r["path"] == gp and r["variant"] == "all"}
            for q in ("res_cb", "mu_ep", "tail", "res_cb_oos", "mu_ep_oos"):
                rows.append(dict(table="bench", run_type=rt, calib="ep", arm=arm, quantity=q, value=fnum(s[q]), se=fnum(s[f"{q}_se"]), n=fnum(s["n_used"]),
                                 n_splits=fnum(s["n_splits"])))
        for c in CMP:
            for q in ("res_cb", "mu_ep", "tail", "res_cb_oos", "mu_ep_oos"):
                d = np.array([fnum(P[c][k][q]) - fnum(P["A"][k][q]) for k in P["A"] if k in P[c]])
                d = d[np.isfinite(d)]
                rows.append(dict(table="bench_paired", run_type=rt, calib="ep", comparison=f"{c}-A", quantity=q, value=float(d.mean()) if d.size else np.nan,
                                 se=float(d.std(ddof=1) / math.sqrt(d.size)) if d.size > 1 else np.nan, n_splits=int(d.size)))
    for arm in ARMS:
        for r in tsv(Path(root) / arm / "moller" / f"{TA}_moller" / "summary.tsv"):
            if r["path"] == gate_path(TA) and r["quantity"] in ("sum", "cluster"):
                for q in ("res", "res_pooled", "peak"):
                    rows.append(dict(table="bench", run_type=TA, calib="moller", arm=arm, protocol=r["protocol"], quantity=f"{r['quantity']}:{q}", value=fnum(r[q]),
                                     se=fnum(r.get(f"{q}_se", "nan")), n=fnum(r["n"]), n_splits=fnum(r["n_splits"])))
    return rows


def boot_rows(root, nboot, nproc):
    _G["ep"] = prep_ep(root)
    _G["m"] = prep_moller(root)
    ep_jobs = [(rt, c, prot, "overall", -1, nboot) for rt in RTS for c in CMP for prot in ("in", "oos")]
    ep_jobs += [(rt, c, "oos", "ring", b, max(nboot // 2, 50)) for rt in RTS for c in CMP for b in range(len(RADBINS) - 1)]
    m_jobs = [(c, prot, kind, None, nboot) for c in CMP for prot in ("in", "oos") for kind in ("sum", "cluster")]
    m_jobs += [(c, "oos", "cluster", b, max(nboot // 2, 50)) for c in CMP for b in range(len(RADBINS) - 1)]
    with mp.get_context("fork").Pool(nproc) as pool:
        ep_res = pool.map(ep_boot_job, ep_jobs, chunksize=1)
        m_res = pool.map(m_boot_job, m_jobs, chunksize=1)
    rows = []
    for (rt, c, prot, axis, b, nb), r in ep_res:
        S = _G["ep"][rt, c, prot]
        base = dict(table="paired", run_type=rt, calib="ep", comparison=f"{c}-A", protocol=prot, axis=axis,
                    bin="all" if axis == "overall" else f"{RADBINS[b]:g}-{RADBINS[b + 1]:g}", n_A=S["nA"], n_X=S["nX"])
        if r is None:
            rows.append(base | dict(quantity="n_shared", value=int((S["ring"] == b).sum())))
            continue
        for q in EPQ:
            dv = r["X"][q] - r["A"][q]
            rows.append(base | dict(quantity=q, n_shared=r["A"]["n"], value_A=r["A"][q], value_X=r["X"][q], diff=dv, se=r["sd"][q],
                                    z=dv / r["sd"][q] if r["sd"][q] > 0 else np.nan, nboot=r["nboot"]))
    for (c, prot, kind, b, nb), r in m_res:
        S = _G["m"][c]
        base = dict(table="paired", run_type=TA, calib="moller", comparison=f"{c}-A", protocol=prot, axis="overall" if b is None else "ring",
                    bin=(kind if b is None else f"{RADBINS[b]:g}-{RADBINS[b + 1]:g}"), n_A=S["nA"], n_X=S["nX"])
        if r is None:
            rows.append(base | dict(quantity="n_shared", value=0))
            continue
        for q in MQ:
            dv = r["X"][q] - r["A"][q]
            rows.append(base | dict(quantity=q, n_shared=r["A"]["n"], value_A=r["A"][q], value_X=r["X"][q], diff=dv, se=r["sd"][q],
                                    z=dv / r["sd"][q] if r["sd"][q] > 0 else np.nan, nboot=r["nboot"]))
    return rows


def p10_losses():
    """{(dataset, population, veto): (k, n, frac %)} of P10's overall veto table."""
    out = {}
    for r in tsv(P10_TSV):
        if r["table"] == "veto" and r["axis"] == "overall" and r["veto"] in ("V3", "V4"):
            out[r["dataset"], r["population"], r["veto"]] = (int(fnum(r["k"])), int(fnum(r["n"])), 100.0 * fnum(r["frac"]))
    return out


def coverage_rows(root):
    """Per run type: how many C corrections and V3 / V4 keys reach the control-gate pool clusters of the A dump (and the used ones, in sample)."""
    rows = []
    for rt in RTS:
        d = load_dump(root, rt, "A")
        kd = set(zip(d["run"].tolist(), d["evn"].tolist(), d["seed"].tolist()))
        ku = set(zip(d["run"][d["used"]].tolist(), d["evn"][d["used"]].tolist(), d["seed"][d["used"]].tolist()))
        for arm, q in (("C", "corr"), ("V3", "veto"), ("V4", "veto")):
            with np.load(hook_path(root, rt, arm)) as h:
                k = set(map(tuple, h[q][:, :3].tolist()))
            rows.append(dict(table="coverage", run_type=rt, arm=arm, quantity="clusters_with_key", value=len(k)))
            rows.append(dict(table="coverage", run_type=rt, arm=arm, quantity="in_pool", value=len(k & kd)))
            rows.append(dict(table="coverage", run_type=rt, arm=arm, quantity="in_used", value=len(k & ku)))
    return rows


def stage_analyse(root, out, nboot, nproc=12):
    M = module_info(root)
    rows = tsv(Path(root) / "hooks" / "hooks.tsv")
    rows = [dict(r) for r in rows]
    rows += bench_rows(root)
    rows += boot_rows(root, nboot, nproc)
    rows += gain_rows(root, M)
    rows += coverage_rows(root)
    write_rows(out / f"{STEM}.tsv", rows)
    print(f"[tsv] {out / f'{STEM}.tsv'}: {len(rows)} rows", flush=True)
    report(rows, M, out, p10_losses())


def stage_report(root, out):
    """PDF (and the coverage rows) from an existing TSV, without the bootstrap."""
    rows = [dict(r) for r in tsv(out / f"{STEM}.tsv") if r.get("table") != "coverage"] + coverage_rows(root)
    write_rows(out / f"{STEM}.tsv", rows)
    report(rows, module_info(root), out, p10_losses())



# ---------------------------------------------------------------------------------------------------------------- report
QL = {"res_cb": "sigma/E CB (%)", "mu_ep": "peak/E' - 1 (%)", "tail_lo3": "tail < -3 sigma (%)", "tail_hi3": "tail > +3 sigma (%)",
      "res": "sigma/E CB (%)", "peak": "peak/E - 1 (%)"}
Z_SIG = 2.0                                             # |difference| / its paired error above this = beyond the statistical uncertainty


def get(rows, **kw):
    return [r for r in rows if all(str(r.get(k, "")) == str(v) for k, v in kw.items())]


def one(rows, **kw):
    g = get(rows, **kw)
    return g[0] if g else {}


def pm(v, e, nd=3):
    v, e = fnum(v), fnum(e)
    return f"{v:.{nd}f} +- {e:.{nd}f}" if math.isfinite(v) else "n/a"


def dz(r, nd=4):
    d, e = fnum(r.get("diff")), fnum(r.get("se"))
    if not math.isfinite(d):
        return "n/a"
    if d == 0 and e == 0:
        return f"{0:+.{nd}f} (identical)"
    z = d / e if e > 0 else float("nan")
    return f"{d:+.{nd}f} +- {e:.{nd}f} ({z:+.1f})"


def table_lines(rows):
    L = ["Arm values: the benchmark's own numbers on these splits (control-gate path, mean +- SE over the splits; in = in sample, oos = gains / position fitted on the other half).",
         "Differences: paired, pooled fit of the clusters used in both arms (same keys), bootstrap error; value (z = diff / error).  Full-data control v3 oos CB sigma/E for",
         "reference: " + ", ".join(f"{RT_LABEL[rt]} {BASELINE[rt]}" for rt in RTS) + " %.", ""]
    for rt in RTS:
        L.append(f"{RT_LABEL[rt]}  (splits {int(fnum(one(rows, table='hook', run_type=rt, arm='A', quantity='splits').get('value')))}"
                 + (f", seed radius >= {RMIN[rt]:g} mm" if rt in RMIN else "") + ")")
        L.append(f"  {'arm':<14}{'n used':>9}{'sigma/E in':>20}{'sigma/E oos':>20}{'mu in':>20}{'mu oos':>20}{'tail 3 sigma in':>20}")
        for arm in ARMS:
            b = lambda q: one(rows, table="bench", run_type=rt, calib="ep", arm=arm, quantity=q)   # noqa: E731
            L.append(f"  {ARM_LABEL[arm]:<14}{int(fnum(b('res_cb').get('n'))):>9}" + "".join(f"{pm(b(q).get('value'), b(q).get('se')):>20}"
                                                                                       for q in ("res_cb", "res_cb_oos", "mu_ep", "mu_ep_oos", "tail")))
        L.append(f"  {'paired diff':<14}{'n shared':>9}{'d sigma/E oos':>27}{'d mu oos':>27}{'d tail lo oos':>27}{'d tail hi oos':>27}{'d sigma/E in':>27}")
        for c in CMP:
            g = lambda q, prot="oos": one(rows, table="paired", run_type=rt, calib="ep", comparison=f"{c}-A", protocol=prot, axis="overall", quantity=q)   # noqa: E731
            L.append(f"  {ARM_LABEL[c] + ' - A':<14}{int(fnum(g('res_cb').get('n_shared'))):>9}" + "".join(f"{dz(g(q)):>27}" for q in EPQ)
                     + f"{dz(g('res_cb', 'in')):>27}")
        L.append("")
    L.append(f"{RT_LABEL[TA]}, Moller two-cluster calibration (pairs, gate path; sum = E1 + E2 vs Ebeam + m_e - losses, cluster = E / E_pred)")
    L.append(f"  {'arm':<14}{'pairs':>8}{'sum sigma/E oos':>22}{'pooled':>22}{'peak oos':>10}{'cluster sigma/E oos':>24}")
    for arm in ARMS:
        b = lambda q: one(rows, table="bench", calib="moller", arm=arm, protocol="out-of-sample", quantity=q)   # noqa: E731
        L.append(f"  {ARM_LABEL[arm]:<14}{int(fnum(b('sum:res').get('n'))):>8}{pm(b('sum:res').get('value'), b('sum:res').get('se')):>22}"
                 f"{pm(b('sum:res_pooled').get('value'), b('sum:res_pooled').get('se')):>22}{fnum(b('sum:peak').get('value')):>10.3f}"
                 f"{pm(b('cluster:res').get('value'), b('cluster:res').get('se')):>24}")
    L.append(f"  {'paired diff':<14}{'pairs':>8}{'d sum sigma/E oos':>27}{'d sum peak oos':>27}{'d sum tail lo':>27}{'d sum tail hi':>27}{'d cluster sigma/E oos':>27}")
    for c in CMP:
        g = lambda q, kind="sum": one(rows, table="paired", calib="moller", comparison=f"{c}-A", protocol="oos", axis="overall", bin=kind, quantity=q)   # noqa: E731
        L.append(f"  {ARM_LABEL[c] + ' - A':<14}{int(fnum(g('res').get('n_shared'))):>8}" + "".join(f"{dz(g(q)):>27}" for q in MQ) + f"{dz(g('res', 'cluster')):>27}")
    return L


def ring_lines(rows):
    L = ["Per P4 ring of the seed-module centre radius, out of sample, paired as above (rings with < %d shared clusters omitted; ep: the calibration pool reaches 200 mm," % MIN_RING,
         "so the 170-230 ring is 170-200 mm; Moller: the cluster position radius, E / E_pred x <E_pred>).", ""]
    for rt, calib in [(rt, "ep") for rt in RTS] + [(TA, "moller")]:
        q1, q2 = ("res_cb", "mu_ep") if calib == "ep" else ("res", "peak")
        L.append(f"{RT_LABEL[rt]} ({calib})")
        L.append(f"  {'ring (mm)':<11}{'n':>8}{'sigma/E A':>11}" + "".join(f"{'d sigma/E ' + c:>27}" for c in CMP) + f"{'mu A':>9}" + "".join(f"{'d mu ' + c:>27}" for c in CMP))
        for b in range(len(RADBINS) - 1):
            lab = f"{RADBINS[b]:g}-{RADBINS[b + 1]:g}"
            g = lambda c, q: one(rows, table="paired", run_type=rt, calib=calib, comparison=f"{c}-A", protocol="oos", axis="ring", bin=lab, quantity=q)   # noqa: E731
            r0 = g("C", q1)
            if not r0:
                continue
            L.append(f"  {lab:<11}{int(fnum(r0.get('n_shared'))):>8}{fnum(r0.get('value_A')):>11.3f}" + "".join(f"{dz(g(c, q1)):>27}" for c in CMP)
                     + f"{fnum(g('C', q2).get('value_A')):>9.3f}" + "".join(f"{dz(g(c, q2)):>27}" for c in CMP))
        L.append("")
    return L


def gain_summary(rows):
    """{(rt, calib, cmp): dict(n, n1, n2, worst row)} over the free (data-set) PbWO4 modules of each calibration."""
    out = {}
    for r in rows:
        if r.get("table") != "gain" or not int(fnum(r["free"])):
            continue
        k = (r["run_type"], r["calib"], r["comparison"])
        o = out.setdefault(k, dict(n=0, n1=0, n2=0, worst=None))
        z = abs(fnum(r["z"]))
        o["n"] += 1
        o["n1"] += int(z > 1)
        o["n2"] += int(z > Z_SIG)
        if o["worst"] is None or z > abs(fnum(o["worst"]["z"])):
            o["worst"] = r
    return out


def gain_lines(rows, top=12):
    gs = gain_summary(rows)
    L = ["Per-module gain change dg = g_X - g_A per run (ep: pool solve, ring-scoped solve for modules that only it sets; Moller: per-run in-sample solve), in units of",
         "sigma_g of A (the solve's statistical error).  Free = prior fraction < 0.5 (the data set the gain); the prior-held modules move with g0 only.", ""]
    L.append(f"  {'calibration':<26}{'comparison':<8}{'free modules x runs':>21}{'|dg| > 1 sigma':>16}{'|dg| > 2 sigma':>16}   largest |dg| / sigma_g")
    for (rt, calib, c), o in sorted(gs.items(), key=lambda x: (RTS.index(x[0][0]), x[0][1], x[0][2])):
        w = o["worst"]
        L.append(f"  {RT_LABEL[rt] + ' ' + calib:<26}{c:<8}{o['n']:>21}{o['n1']:>16}{o['n2']:>16}   {w['module']} run {w['run']} dg {fnum(w['dg']):+.5f} "
                 f"({100 * fnum(w['dg']) / fnum(w['g_A']):+.3f} %, {fnum(w['z']):+.2f} sigma)")
    L += ["", f"Largest |dg| / sigma_g (free modules, top {top} over every calibration and comparison):"]
    fr = sorted([r for r in rows if r.get("table") == "gain" and int(fnum(r["free"]))], key=lambda r: -abs(fnum(r["z"])))[:top]
    L.append(f"  {'calibration':<26}{'cmp':<6}{'run':>7}{'module':>8}{'ring':>10}{'g_A':>10}{'sigma_g':>10}{'dg':>11}{'dg %':>9}{'z':>8}{'n_seed':>8}")
    for r in fr:
        L.append(f"  {RT_LABEL[r['run_type']] + ' ' + r['calib']:<26}{r['comparison']:<6}{r['run']:>7}{r['module']:>8}{r['ring']:>10}{fnum(r['g_A']):>10.5f}"
                 f"{fnum(r['sigma_g']):>10.5f}{fnum(r['dg']):>+11.5f}{100 * fnum(r['dg']) / fnum(r['g_A']):>+9.3f}{fnum(r['z']):>+8.2f}{int(fnum(r['n_seed'])):>8}")
    L += ["", "Named modules (Ta inner ring W564 W559 W493 W663 W593; W398 high-gain, W457 fakes, W613 f4): dg / sigma_g per run, C - A | V4 - A | V3 - A  (* = prior-held)"]
    for rt, calib in [(rt, "ep") for rt in RTS] + [(TA, "moller")]:
        parts = []
        for mname in NAMED:
            rr = get(rows, table="gain", run_type=rt, calib=calib, module=mname)
            if not rr:
                continue
            runs = sorted({r["run"] for r in rr})
            z = lambda c, run: fnum(one(rr, comparison=f"{c}-A", run=run).get("z"))   # noqa: E731
            star = "" if all(int(fnum(one(rr, comparison="C-A", run=run).get("free", 0))) for run in runs) else "*"
            parts.append(f"{mname}{star} " + " ".join(f"{z('C', run):+.2f}|{z('V4', run):+.2f}|{z('V3', run):+.2f}" for run in runs))
        L.append(f"  {RT_LABEL[rt] + ' ' + calib:<26}" + ("; ".join(parts) if parts else "none in the solve"))
    return L


def gate_eval(rows):
    """{(rt, calib, cmp): dict(flags, text)}: does the arm change sigma/E, mu (oos, overall, |z| > Z_SIG) or any free module gain beyond Z_SIG sigma_g?"""
    gs = gain_summary(rows)
    out = {}
    for rt, calib in [(rt, "ep") for rt in RTS] + [(TA, "moller")]:
        for c in CMP:
            if calib == "ep":
                g = lambda q: one(rows, table="paired", run_type=rt, calib="ep", comparison=f"{c}-A", protocol="oos", axis="overall", quantity=q)   # noqa: E731
                rs, mu = g("res_cb"), g("mu_ep")
            else:
                g = lambda q: one(rows, table="paired", calib="moller", comparison=f"{c}-A", protocol="oos", axis="overall", bin="sum", quantity=q)   # noqa: E731
                rs, mu = g("res"), g("peak")
            o = gs.get((rt, calib, f"{c}-A"), dict(n=0, n1=0, n2=0, worst=None))
            zr, zm = fnum(rs.get("z")), fnum(mu.get("z"))
            out[rt, calib, c] = dict(res=abs(zr) > Z_SIG, mu=abs(zm) > Z_SIG, gain=o["n2"] > 0, n2=o["n2"], zr=zr, zm=zm, n1=o["n1"], n=o["n"],
                                     dres=fnum(rs.get("diff")), dmu=fnum(mu.get("diff")), worst=o["worst"])
    return out


def gate_lines(rows, p10):
    G = gate_eval(rows)
    L = [f"Gate test per run type: does the arm move sigma/E or the peak (out of sample, paired, |z| > {Z_SIG:g}) or any free module gain by more than {Z_SIG:g} sigma_g (same",
         "threshold; the 1-sigma counts are on the gain page)?  The verdict uses the out-of-sample protocol; in-sample rows are in the three-arm table and the TSV.",
         "P10 cost = elastic-window loss of the veto (Moller pairs for Ta Moller), Wilson 68 %, from exotic_P10_yield_impact_20261118.tsv.", ""]
    L.append(f"  {'calibration':<26}{'arm':<6}{'d sigma/E (z)':>16}{'d mu (z)':>16}{'gains > 2 sigma':>17}   {'changes?':<10}{'P10 cost':>24}")
    for (rt, calib, c), o in G.items():
        ds = "2p24GeV_Ta_X17" if rt == TA else rt
        cost = p10.get((ds, "moller_pair" if calib == "moller" else "elastic", c)) if c in ("V3", "V4") else None
        cs = f"{cost[2]:.4f} % ({cost[0]}/{cost[1]})" if cost else ("none (no veto)" if c == "C" else "n/a")
        ch = "yes" if (o["res"] or o["mu"] or o["gain"]) else "no"
        zs = lambda z, d: "(  --)" if d == 0 else f"({z:+.1f})"   # noqa: E731  identical arms
        L.append(f"  {RT_LABEL[rt] + ' ' + calib:<26}{c:<6}{o['dres']:>+9.4f} {zs(o['zr'], o['dres'])}{o['dmu']:>+9.4f} {zs(o['zm'], o['dmu'])}{o['n2']:>9} of {o['n']:<5}   {ch:<10}{cs:>24}")
    L += [""] + RECOMMEND
    return L


RECOMMEND = [   # judgement from the 2026-11-18 run (round 2, ambiguous P8 keys kept at the default; numbers: this page, the three-arm table, the gain page)
    "Reading: no arm moves sigma/E or the peak beyond 2 sigma of the paired error at any run type (largest |z| over the gated quantities, oos overall: 3.49 GeV",
    "C peak +0.0009 +- 0.0005 %, z 1.8; Ta Moller V4 sum peak +0.0055 %, z 1.4).  In sample the same holds (largest: 3.49 GeV V4 peak -0.0012 +- 0.0006 %, z -1.9).",
    "No free module gain moves by 2 sigma_g; at most 1 module-run per calibration passes 1 sigma_g (Ta Moller V4 W659 run 26230 -0.87 %, -1.98 sigma; 0.73 V4",
    "W567 -0.56 %, -1.73; Ta ep C W626 +0.24 %, +1.39; Ta ep V4 W566 +0.28 %, +1.28): refit noise from a few hundred removed or changed events, not a pattern.",
    "The C arm's energy-scale shift (P8: Ta elastic -0.017 %) is absorbed by the gain refit: the peak stays put (Ta ep -0.0004 +- 0.0047 %) and no gain moves by",
    "more than 0.4 % (at most 1.4 sigma_g).  The Ta inner-ring modules (W564 W559 W493 W663 W593) sit in the central block H that seeds neither Ta calibration: their",
    "gains are prior-held in both, so P9 cannot measure their dgain on Ta (on 2.24 / 3.49 GeV H2, where the ring calibration frees them, |dg| < 0.2 sigma_g).",
    "W398 is prior-held at 0.73 GeV and on Ta, free on some 2.24 / 3.49 GeV runs with dg ~ 0; its clipped showers remain a P8 energy-estimate issue.",
    "",
    "Recommendation for the gate (calibration and resolution side; yields from P10):",
    "  H2 2.24 / 3.49 GeV: keep V3 as adopted (it removes the fakes at 0.004 / 0.007 % elastic cost and is neutral for calibration and resolution);",
    "     do not add V4 (no gain, 0.008 % cost at 2.24); no in-place correction (C changes nothing measurable).  Ignore f1 / f3 / f2 / f5 otherwise.",
    "  H2 0.73 GeV: ignore - V3 has no key in the pool, V4 removes 10 used clusters with no effect, and no corrected cluster reaches the pool (0 of 146",
    "     keyed clusters; almost all seed below 107 mm, the rest fail the pool selection).",
    "  Ta X17 (ep and Moller): ignore - V3 is a no-op, V4 costs 0.29 % of elastic and 0.12 % of Moller pairs (inner ring 0.31 %) for no measurable change, and",
    "     the credible correction (f1 sloped baseline, f3 stretch) is absorbed by the gains.  The flags stay event markers.",
]
CAVEATS = [
    "First order: the C arm rescales the flagged members' raw energies inside the benchmark's clusters; clustering, seed choice, position, the one-cluster pool selection",
    "  and the Moller candidate list are those of the default.  The vetoes drop P7's clusters by (run, event number, seed module) on every path (legacy included,",
    "  where a legacy cluster with the same seed stands in for the gated one), so the gain solves (legacy pool) and the fits both lose them.",
    "  The legacy clusters (gain-solve input) get the gate-path rho = dE / E_def by the same key (first order); a key held by two P8 records (two time-split clusters of",
    "  one event with the same seed) is ambiguous and keeps the default.",
    "Scope: H2 run types on the 2-file exotic sample (4 / 10 / 4 splits at 0.73 / 2.24 / 3.49 GeV), Ta X17 on 19 splits of 25800 and 26230 (26138 not in the full",
    "  scan).  The calibrations are refit on these splits only (one split per run per out-of-sample half for H2), so the arm values sit above the full-data control v3",
    "  numbers; the comparison is between arms on identical inputs, not with the full-data chain.",
    "Paired errors: bootstrap of the shared clusters with the gains fixed per arm; they measure whether the energies of the same events moved, not the calibration's",
    "  run-to-run spread.  Gain changes are compared with the solve's own sigma_g of A.",
    "P8's credible dE sets a member to zero (rho = -1) when its pulse disappears in the alternative fit (f1 sloped baseline ~70 %, f3 stretch ~30 %): 81k of 367k Ta corrections.",
    "f5: k = 1 extrapolation only up to %d overflow samples (P8 'mild'); deeper saturation keeps the default in C (veto territory)." % NCLIP_MILD,
]


def fig_diffs(pdf, plt, rows):
    fig, axs = plt.subplots(2, 2, figsize=(16, 11))
    cols = {"V3": "tab:red", "V4": "tab:purple", "C": "tab:blue"}
    cal = [(rt, "ep") for rt in RTS] + [(TA, "moller")]
    for ax, (q_ep, q_m, title) in zip(axs[0], (("res_cb", "res", "d sigma/E (oos, percentage points)"), ("mu_ep", "peak", "d peak (oos, %)"))):
        for j, c in enumerate(CMP):
            xs, ys, es = [], [], []
            for i, (rt, calib) in enumerate(cal):
                r = one(rows, table="paired", run_type=rt, calib=calib, comparison=f"{c}-A", protocol="oos", axis="overall",
                        quantity=q_ep if calib == "ep" else q_m, **({} if calib == "ep" else {"bin": "sum"}))
                xs.append(i + (j - 1) * 0.2)
                ys.append(fnum(r.get("diff")))
                es.append(fnum(r.get("se")))
            ax.errorbar(xs, ys, es, fmt="o", color=cols[c], label=f"{c} - A", capsize=3)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(range(len(cal)))
        ax.set_xticklabels([RT_LABEL[rt] + ("\nMoller sum" if calib == "moller" else "\nep") for rt, calib in cal], fontsize=9)
        ax.set_title(title)
        ax.legend(fontsize=9)
    for ax, q in zip(axs[1], ("res_cb", "mu_ep")):
        for j, c in enumerate(CMP):
            for i, rt in enumerate(RTS):
                xs, ys, es = [], [], []
                for b in range(len(RADBINS) - 1):
                    r = one(rows, table="paired", run_type=rt, calib="ep", comparison=f"{c}-A", protocol="oos", axis="ring", bin=f"{RADBINS[b]:g}-{RADBINS[b + 1]:g}", quantity=q)
                    if r:
                        xs.append(i * 4 + b + (j - 1) * 0.2)
                        ys.append(fnum(r["diff"]))
                        es.append(fnum(r["se"]))
                ax.errorbar(xs, ys, es, fmt="o", ms=4, color=cols[c], label=f"{c} - A" if i == 0 else None, capsize=2)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks([i * 4 + 1 for i in range(len(RTS))])
        ax.set_xticklabels([RT_LABEL[rt] for rt in RTS], fontsize=9)
        ax.set_title(f"ep pool per P4 ring (0-120, 120-170, 170-200 mm), d {'sigma/E' if q == 'res_cb' else 'peak'} oos")
        ax.legend(fontsize=9)
    fig.suptitle("Paired differences to the default arm (bootstrap errors)", fontsize=15)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


def fig_maps(pdf, plt, rows, M):
    from matplotlib.colors import TwoSlopeNorm
    idx = {n: i for i, n in enumerate(M["name"])}
    for rt, calib in [(rt, "ep") for rt in RTS] + [(TA, "moller")]:
        fig, axs = plt.subplots(1, 3, figsize=(16, 6.2))
        for ax, c in zip(axs, CMP):
            rr = [r for r in get(rows, table="gain", run_type=rt, calib=calib, comparison=f"{c}-A")]
            agg = {}
            for r in rr:
                agg.setdefault(r["module"], []).append((fnum(r["z"]), int(fnum(r["free"]))))
            pw = np.flatnonzero(M["pwo4"])
            ax.scatter(M["x"][pw], M["y"][pw], s=9, marker="s", color="0.9")
            mods = [m for m in agg if m in idx]
            zz = np.array([np.mean([z for z, _ in agg[m]]) for m in mods])
            fr = np.array([any(f for _, f in agg[m]) for m in mods])
            ii = np.array([idx[m] for m in mods], int)
            lim = max(1.0, float(np.nanmax(np.abs(zz))) if zz.size else 1.0)
            sc = None
            if zz.size:
                sc = ax.scatter(M["x"][ii[fr]], M["y"][ii[fr]], c=zz[fr], s=12, marker="s", cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim))
                ax.scatter(M["x"][ii[~fr]], M["y"][ii[~fr]], c=zz[~fr], s=4, marker=".", cmap="RdBu_r", norm=TwoSlopeNorm(0, -lim, lim))
            for mname in NAMED:
                if mname in idx:
                    ax.annotate(mname, (M["x"][idx[mname]], M["y"][idx[mname]]), fontsize=6)
            ax.set_aspect("equal")
            ax.set_xlim(-560, 560)
            ax.set_ylim(-560, 560)
            n1 = int(sum(1 for r in rr if int(fnum(r["free"])) and abs(fnum(r["z"])) > 1))
            ax.set_title(f"{c} - A: mean dg / sigma_g over runs\n(squares: free modules; |z| > 1: {n1} module-runs)", fontsize=10)
            if sc is not None:
                fig.colorbar(sc, ax=ax, shrink=0.75)
        fig.suptitle(f"Per-module gain change, {RT_LABEL[rt]} ({calib} calibration)", fontsize=14)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)


def hook_lines(rows):
    L = ["Inputs per run type (hook files): splits, veto keys (P7 clusters in scope), member corrections (P8 records with a credible non-zero dE):", ""]
    L.append(f"  {'run type':<18}{'splits':>8}{'V3 keys':>9}{'V4 keys':>9}{'C rows':>9}{'  of P8 records':>17}{'f1':>9}{'f3':>8}{'f5 mild':>9}{'f5 deep (kept)':>16}{'ambiguous (kept)':>18}")
    for rt in RTS:
        g = lambda arm, q: int(fnum(one(rows, table="hook", run_type=rt, arm=arm, quantity=q).get("value")))   # noqa: E731
        L.append(f"  {RT_LABEL[rt]:<18}{g('A', 'splits'):>8}{g('V3', 'veto_keys'):>9}{g('V4', 'veto_keys'):>9}{g('C', 'corr_rows'):>9}{g('C', 'p8_records'):>17}"
                 f"{g('C', 'p8_f1'):>9}{g('C', 'p8_f3'):>8}{g('C', 'p8_f5'):>9}{g('C', 'p8_f5_deep'):>16}{g('C', 'p8_ambiguous'):>18}")
    L += ["", "Keys that reach the benchmark's control-gate pool (A dump: all pool clusters | used in sample) of the keyed clusters:"]
    for rt in RTS:
        g = lambda arm, q: int(fnum(one(rows, table="coverage", run_type=rt, arm=arm, quantity=q).get("value")))   # noqa: E731
        L.append(f"  {RT_LABEL[rt]:<18}" + "   ".join(f"{arm}: {g(arm, 'in_pool')} | {g(arm, 'in_used')} of {g(arm, 'clusters_with_key')}" for arm in ("V3", "V4", "C")))
    L += ["", "Clusters in each arm's dump (control gate, pool) and the shared sets of the paired comparisons (out of sample):"]
    for rt in RTS:
        r = {c: one(rows, table="paired", run_type=rt, calib="ep", comparison=f"{c}-A", protocol="oos", axis="overall", quantity="res_cb") for c in CMP}
        L.append(f"  {RT_LABEL[rt]:<18}A used {int(fnum(r['C'].get('n_A'))):>8}" + "".join(f"   {c} used {int(fnum(r[c].get('n_X'))):>8} shared {int(fnum(r[c].get('n_shared'))):>8}" for c in CMP))
    r = {c: one(rows, table="paired", calib="moller", comparison=f"{c}-A", protocol="oos", axis="overall", bin="sum", quantity="res") for c in CMP}
    L.append(f"  {'Ta X17 Moller':<18}A pairs {int(fnum(r['C'].get('n_A'))):>7}" + "".join(f"   {c} pairs {int(fnum(r[c].get('n_X'))):>7} shared {int(fnum(r[c].get('n_shared'))):>7}" for c in CMP))
    return L


def report(rows, M, out, p10):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    pdf_path = out / f"{STEM}.pdf"
    with PdfPages(pdf_path) as pdf:
        text_page(pdf, plt, "P9 exotic events: impact on calibration and resolution (arms A default, B veto V3 / V4, C corrected)",
                  [x.strip() for x in __doc__.split("Arms")[1].split("Outputs")[0].splitlines()] + [""] + hook_lines(rows), size=8.6)
        text_page(pdf, plt, "Three-arm table: sigma/E, peak, tails (ep per run type; Ta Moller)", table_lines(rows), size=7.4)
        fig_diffs(pdf, plt, rows)
        text_page(pdf, plt, "Per ring: paired differences of sigma/E and the peak", ring_lines(rows), size=7.6)
        fig_maps(pdf, plt, rows, M)
        text_page(pdf, plt, "Per-module gain changes: extremes and named modules", gain_lines(rows), size=7.6)
        text_page(pdf, plt, "Gate after P9 / P10", gate_lines(rows, p10), size=8.2)
        text_page(pdf, plt, "Caveats", CAVEATS, size=9.5)
    print(f"[pdf] {pdf_path}", flush=True)


# ---------------------------------------------------------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", default="all", choices=("hooks", "run", "analyse", "report", "all"))
    ap.add_argument("--jobs", type=int, default=2, help="reports at a time (each uses the benchmark's NPROC workers)")
    ap.add_argument("--root", default=str(SCAN_ROOT / "exotic_p9"))
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--nboot", type=int, default=NBOOT)
    a = ap.parse_args()
    if a.stage in ("hooks", "all"):
        stage_hooks(a.root)
    if a.stage in ("run", "all"):
        stage_run(a.root, a.jobs)
    if a.stage in ("analyse", "all"):
        stage_analyse(a.root, Path(a.out), a.nboot)
    if a.stage == "report":
        stage_report(a.root, Path(a.out))


if __name__ == "__main__":
    main()
