#!/usr/bin/env python3
"""moller_calib.py -- Moller two-cluster calibration and energy resolution of the 2p24GeV_Ta_X17 runs (two-cluster trigger).

    python3 waveform_analysis/energy_resolution_benchmark.py report --calib moller [--out <dir>] <scan dir> [<gate_ns>]
    python3 waveform_analysis/moller_calib.py report [--out <dir>] [--ref <v3 report root>] [--cache <dir>] [--cl-min none|50|100|200|300]
                                                     [--r-iso <mm>] [--combine-runs] <scan dir> [<gate_ns>]
    python3 waveform_analysis/moller_calib.py test          # unit tests of the Moller kinematics

The control-v3 ep path (`report` without --calib, or --calib ep) is not touched: this module imports the benchmark and reuses its
hit preprocessing (residual onset-phase factor, time offsets), clusterers, Crystal Ball fit, gain normal equations and convergence
test, and the borrowed 2.24 GeV H2 position correction of the control-v3 report root (--ref).

Chain (README: energy_resolution/2p24GeV_Ta_X17_moller/README.md)
  1. Clustering per split, legacy path and the gate path (default CONTROL_GATE of the run type, 4 ns): candidates = the two largest
     clusters of an event, both PbWO4-seeded, loosely back to back and loosely at Ebeam (cached under CACHE_ROOT).
  2. Selection: isolation, no other cluster above CL_STORE kE within R_ISO of either pair cluster (--r-iso; 0 = off); optional global
     third-cluster threshold, exactly two clusters above CL_MIN kE (--cl-min, off by default; 300 = the previous default); both seeds PbWO4 outside the central-hole block H and the PbWO4 / lead-glass
     transition layers T (module_regions) with seed energy >= SEED_MIN kE; energy sum (calibration-file energies) within ESUM_TOL of
     Ebeam; coplanar | |phi1 - phi2| - 180 | < DPHI_CUT about the beam centre; Moller angle relation
     |tan th1 tan th2 (gamma + 1) / 2 - 1| < PROD_CUT.  Beam centre (x0, y0) and vertex distance z per run from the selected pairs
     (centre: least-squares crossing of the pair lines; z: median of sqrt((gamma + 1) / 2 R1 R2), as PhysicsTools::GetMollerZdistance),
     iterated GEO_ITER times.
  3. E_pred(theta) = Moller E'(theta) at Ebeam minus the energy loss (benchmark e_moller), tan theta = R / z about (x0, y0).
  4. Gains: one per member module, sum_i g_i E_i of each selected cluster matched to its E_pred (both clusters of a pair), ridge prior to
     g0 as solve_gains (same lambda form, module windows, refit passes and convergence test; a module without its own window uses the
     global one); legacy path, one solve per run as the ep chain (`--combine-runs`: one over all runs); in sample (all splits of the
     run) and out of sample (even / odd split ordinal within the run, the other half's gains).
  5. Corrected E = sum g_i E_i / R(u, v) (borrowed position correction).  Resolution on each path, Crystal Ball as cb_single, per split
     (mean +- standard error) and pooled (bootstrap over the splits): (i) E / E_pred per cluster, (ii) E1 + E2 (angle independent; also
     the symmetric pairs, sqrt(2) sigma_sum / E_sum = single-cluster sigma / E), (iii) E / E_pred in N_EBINS bins of E_pred; compared
     with the H2 curve H2_CURVE.
"""
from __future__ import annotations

import hashlib
import json
import inspect
import math
import multiprocessing as mp
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import energy_resolution_benchmark as erb  # noqa: E402

RT = "2p24GeV_Ta_X17"
CL_MIN = None                          # MeV x kE or None: global third-cluster threshold, exactly two clusters above it (None = off; 300 = the previous default)
R_ISO = 60.0                           # mm: isolation, no other cluster above CL_STORE kE within this of either pair cluster (0 = off); README
COMPARE = ((None, 60.0), (None, 100.0), (300.0, 0.0), (100.0, 0.0), (50.0, 0.0))   # (CL_MIN, R_ISO) of the comparison table
RUN_DATE = {25800: "2026-07-12", 26138: "2026-08-06", 26230: "2026-08-12"}   # run start (notes/PbWO4_energy_resolution_20261028.md, elog summary 11-17)
CL_COUNT = (50.0, 100.0, 200.0, 300.0)  # MeV x kE: thresholds of the stored cluster-multiplicity histograms (CL_MIN must be one; 50 = any other cluster)
CL_STORE = 50.0                        # MeV x kE: clusters considered at all (builder min_cluster_energy)
DPHI_STORE, ESUM_STORE = 30.0, 0.30    # loose candidate prefilter (deg about the origin; fraction of Ebeam)
ESUM_TOL = 0.20                        # |E1 + E2 - Ebeam| < this x Ebeam (calibration-file energies)
DPHI_CUT = 10.0                        # deg: | |phi1 - phi2| - 180 | below this (about the beam centre)
PROD_CUT = 0.15                        # |tan th1 tan th2 (gamma + 1) / 2 - 1| below this (core sigma ~0.058: ~2.6 sigma)
GEO_ITER = 3                           # beam-centre / vertex-distance iterations
SYM_MAX = 0.15                         # symmetric pairs: |E1_pred - E2_pred| < this x (E1_pred + E2_pred)
SYM_Q = f"sum, symmetric pairs |E1 - E2| < {SYM_MAX:g} (E1 + E2)"
N_EBINS = 4                            # E_pred bins (quantiles of the selected gate-path clusters)
N_EXAMPLES = 20                        # example events in the PDF
NPROC = 4                              # worker processes (one split each)
REF_ROOT = erb.PYDIR / "output" / "energy_resolution" / "control_phase2_v3_20261114" / "v3_defaults"   # control-v3 report root
CACHE_ROOT = erb.PYDIR / "output" / "energy_resolution" / "moller_x17"                                 # clustered pairs per gate
OUT_ROOT = erb.OUT_ROOT                # report output: <out>/2p24GeV_Ta_X17_moller/
H2_CURVE = (0.9015, 2.05, 1.35)        # %: sigma/E = a (+) b / sqrt(E/GeV) (+) c / (E/GeV), energy_resolution/resolution_fit_20261115.tsv
PAIR_F = ("seed", "seed_e", "e", "er", "x", "y", "pu", "clip")   # per cluster; members as flat arrays + counts


# kinematics ==================================================================
def moller_e(theta_deg, ebeam):
    """Lab energy (MeV, total) of a Moller electron at polar angle theta (deg); ebeam = beam total energy (MeV); no energy loss."""
    c2 = np.cos(np.radians(theta_deg)) ** 2
    g = ebeam / erb.M_ELECTRON
    return erb.M_ELECTRON * ((g + 1.0) + (g - 1.0) * c2) / ((g + 1.0) - (g - 1.0) * c2)


def moller_partner(theta_deg, ebeam):
    """Polar angle (deg) of the partner: tan th1 tan th2 = 2 / (gamma + 1)."""
    g = ebeam / erb.M_ELECTRON
    return np.degrees(np.arctan(2.0 / ((g + 1.0) * np.tan(np.radians(theta_deg)))))


def run_tests():
    eb, m = 2239.51, erb.M_ELECTRON
    g = eb / m
    ts = math.degrees(math.atan(math.sqrt(2.0 / (g + 1.0))))          # symmetric angle
    e1 = float(moller_e(ts, eb))
    assert abs(float(moller_partner(ts, eb)) - ts) < 1e-9
    assert abs(e1 - (eb + m) / 2.0) < 1e-6 * eb, (e1, (eb + m) / 2)
    for th in (0.5, 0.8, 1.0, ts, 1.5, 2.0, 3.0):
        t2 = float(moller_partner(th, eb))
        s = float(moller_e(th, eb) + moller_e(t2, eb))
        assert abs(s - (eb + m)) < 1e-6 * eb, (th, s, eb + m)
        # the benchmark's e_moller is the same formula minus the energy loss
        assert abs(float(erb.e_moller(th, eb) + erb.energy_loss(th)) - float(moller_e(th, eb))) < 1e-6 * eb
    # momentum balance: p1 sin th1 = p2 sin th2, p1 cos th1 + p2 cos th2 = p_beam
    th = 0.9
    t2 = float(moller_partner(th, eb))
    p = lambda e: math.sqrt(e * e - m * m)  # noqa: E731
    p1, p2 = p(float(moller_e(th, eb))), p(float(moller_e(t2, eb)))
    r = math.radians
    assert abs(p1 * math.sin(r(th)) - p2 * math.sin(r(t2))) < 1e-6 * eb
    assert abs(p1 * math.cos(r(th)) + p2 * math.cos(r(t2)) - p(eb)) < 1e-6 * eb
    print(f"[test] Moller kinematics ok: symmetric angle {ts:.4f} deg, E = {e1:.3f} MeV = (Ebeam + m_e) / 2; "
          f"E1 + E2 = Ebeam + m_e and momentum balance at 7 angles")


# clustering ==================================================================
def pair_job(job):
    """One split: per path (legacy, gate), the two largest clusters of the events whose two largest clusters (above CL_STORE kE) are
    PbWO4-seeded, loosely back to back about the origin (DPHI_STORE) and loosely at Ebeam (ESUM_STORE, calibration-file energies),
    consecutive, with e3 = energy of the event's third-largest cluster (0 without one), dn / en = distance (mm) from the cluster to the nearest
    other cluster above CL_STORE kE and its energy (NaN / 0 without one); and the number of prefiltered events (hit sum
    above (1 - ESUM_TOL) Ebeam) by their number of clusters above each CL_COUNT threshold x kE (index 5 = 5 or more)."""
    fn, split, gate = job
    z = np.load(fn, allow_pickle=False)
    run = int(z["run"])
    cls, geo = erb.clusterers(run, str(z["calib"]), (gate,))
    kE, eb = round(geo["ebeam"] / erb.E_REF, 6), geo["ebeam"]
    evs, h_ev = np.flatnonzero(z["ev_split"] == split), z["h_ev"]
    lo, hi = (np.searchsorted(h_ev, [evs[0], evs[-1] + 1]) if evs.size else (0, 0))
    h = {q: z[q][lo:hi] for q in ("h_mod", "h_t", "h_e", "h_k", "h_extra")} | {"h_ev": h_ev[lo:hi]}
    if erb.PHASE_RESID:
        h = erb.phase_resid(h, z, geo, split, fn)
    if erb.EXOTIC_HOOK:                            # P13: dropped channels removed before clustering
        h = erb.exotic_drop(h, run, z["ev_number"], True)
    ev, pwo4, cstore = h["h_ev"], geo["pwo4"], CL_STORE * kE
    out = {p: {q: [] for q in PAIR_F + ("ev", "e3", "dn", "en", "mi", "me")} for p in cls}
    ncl = {p: np.zeros((len(CL_COUNT), 6), np.int64) for p in cls}
    if ev.size:
        lo_ = np.flatnonzero(np.r_[True, ev[1:] != ev[:-1]])
        hi_ = np.r_[lo_[1:], ev.size]
        keep = np.add.reduceat(h["h_e"].astype(np.float64), lo_) > (1.0 - ESUM_TOL) * eb   # the energy-sum cut needs this
        mod, e, k, extra = (h[q].tolist() for q in ("h_mod", "h_e", "h_k", "h_extra"))
        t = (h["h_t"].astype(np.float64) - geo["toff"][h["h_mod"]]).tolist() if geo.get("toff") is not None else h["h_t"].tolist()
        for s, f in zip(lo_[keep].tolist(), hi_[keep].tolist()):
            multi = {m for m, c in Counter(mod[s:f]).items() if c >= 2}
            wa = {m for m, kk in zip(mod[s:f], k[s:f]) if kk == 0}
            for p, cl in cls.items():
                cl.clear()
                for i in range(s, f):
                    if p != "legacy" or not extra[i]:
                        cl.add_hit(mod[i], e[i], t[i])
                cl.form_clusters()
                big = sorted((r for r in cl.reconstruct_matched() if r.hit.energy > cstore), key=lambda r: -r.hit.energy)
                for j, c in enumerate(CL_COUNT):
                    ncl[p][j, min(sum(r.hit.energy > c * kE for r in big), 5)] += 1
                if len(big) < 2 or not (pwo4[big[0].cluster.center.index] and pwo4[big[1].cluster.center.index]):
                    continue
                a, b = big[0].hit, big[1].hit
                if dphi_dev(a.x, a.y, b.x, b.y) >= DPHI_STORE:
                    continue
                mis = [np.array([x.index for x in r.cluster.hits], np.int64) for r in big[:2]]
                mes = [np.array([x.energy for x in r.cluster.hits]) for r in big[:2]]
                if abs(mes[0].sum() + mes[1].sum() - eb) >= ESUM_STORE * eb:
                    continue
                for r, mi, me in zip(big[:2], mis, mes):
                    seed = r.cluster.center.index
                    dd = [(math.hypot(o.hit.x - r.hit.x, o.hit.y - r.hit.y), o.hit.energy) for o in big[2:]]
                    dn, en = min(dd) if dd else (float("nan"), 0.0)          # nearest other cluster (isolation)
                    for q, v in zip(PAIR_F + ("ev", "e3", "dn", "en", "mi", "me"),
                                    (seed, r.cluster.center.energy, r.hit.energy, me.sum(), r.hit.x, r.hit.y,
                                     bool(multi.intersection(mi.tolist())), bool(seed in wa and pwo4[seed]), int(ev[s]),
                                     big[2].hit.energy if len(big) > 2 else 0.0, dn, en, mi, me)):
                        out[p][q].append(v)
    res = {}
    for p, d in out.items():
        n = len(d["e"])
        res[p] = {q: np.array(d[q]) for q in PAIR_F + ("ev", "e3", "dn", "en")} | dict(
            run=np.full(n, run), split=np.full(n, split), cnt=np.array([x.size for x in d["mi"]], np.int64),
            mi=np.concatenate(d["mi"]) if n else np.zeros(0, np.int64), me=np.concatenate(d["me"]) if n else np.zeros(0),
            ncl=ncl[p], n_events=np.int64(evs.size))
    print(f"  [split] run {run} split {split}: {evs.size} events, {len(out['legacy']['e']) // 2} legacy candidates", flush=True)
    return res


def dphi_dev(x1, y1, x2, y2):
    """| |phi1 - phi2| - 180 | in degrees, folded into [0, 180]."""
    d = np.abs(np.degrees(np.arctan2(y1, x1) - np.arctan2(y2, x2))) % 360.0
    return np.abs(np.minimum(d, 360.0 - d) - 180.0)


def resolve_report_globals():
    """The run-type-dependent benchmark switches that report() resolves before forking (control v3 for the Ta runs)."""
    if erb.BELOW_THR == "auto":
        erb.BELOW_THR = str(erb.TABLE_DIR / erb.BELOW_THR_BY_RT[RT]) if RT in erb.BELOW_THR_BY_RT else ""
    if erb.PHASE_RESID == "auto":
        erb.PHASE_RESID = str(erb.TABLE_DIR / erb.PHASE_RESID_BY_RT[RT]) if RT in erb.PHASE_RESID_BY_RT else ""
    if erb.PC_EXTRA == "auto":
        erb.PC_EXTRA = erb.PC_EXTRA_BY_RT.get(RT, "")
    unsupported = dict(BELOW_THR=erb.BELOW_THR, PC_EXTRA=erb.PC_EXTRA, ASYM_WINDOW=erb.ASYM_WINDOW, SECONDARY_RULE=erb.SECONDARY_RULE,
                       TRAIN_REGIONS=erb.TRAIN_REGIONS, ENERGY_DEF=erb.ENERGY_DEF != "integral", TIME_DEF=erb.TIME_DEF != "mf",
                       GAIN_MODE=erb.GAIN_MODE != "member", NONLINEAR=erb.NONLINEAR, MOLLER_REJECT=erb.MOLLER_REJECT)
    bad = [k for k, v in unsupported.items() if v]
    if bad:
        sys.exit(f"[ERROR] --calib moller supports the control-v3 Ta settings only; set: {', '.join(bad)}")


def collect(scan_dir, gate, cache):
    """Candidates of every split of the scan directory (pair_job), cached in cache/pairs_g<gate>.npz: {path: arrays} plus per split
    (run, split, n_events, ncl)."""
    f = Path(cache) / f"pairs_g{gate:g}.npz"
    files = sorted(Path(scan_dir).glob("*.npz"))
    jobs = []
    want = set(map(tuple, erb.exotic_hook()["splits"].tolist())) if erb.EXOTIC_HOOK else None   # P9: the hook's splits only
    if erb.EXOTIC_HOOK and "drop_calib_only" in erb.exotic_hook() and bool(erb.exotic_hook()["drop_calib_only"]):
        sys.exit("[ERROR] EXOTIC_HOOK drop_calib_only: not supported by --calib moller (one candidate set serves the solve and the evaluation)")
    for fn in files:
        with np.load(fn, allow_pickle=False) as z:
            if str(z["run_type"]) != RT:
                sys.exit(f"[ERROR] {fn}: run type {z['run_type']}, --calib moller is for {RT}")
            if int(z["run"]) in erb.EXCLUDE_RUNS:
                continue
            jobs += [(str(fn), int(s), gate) for s in z["splits"] if want is None or (int(z["run"]), int(s)) in want]
    files_id = []
    for fn in sorted({j[0] for j in jobs}):
        st = Path(fn).resolve().stat()
        with np.load(fn, allow_pickle=False) as z:
            files_id.append((str(Path(fn).resolve()), st.st_size, int(st.st_mtime), str(z["calib"])))
    key = repr(files_id) + repr(sorted((Path(j[0]).name, j[1]) for j in jobs)) + repr((CL_STORE, DPHI_STORE, ESUM_STORE, ESUM_TOL, CL_COUNT, erb.PHASE_RESID,
                                                                       erb.TIME_OFFSETS, erb.OVERRIDES,
                                                                       hashlib.sha1(inspect.getsource(pair_job).encode()).hexdigest()))
    if erb.EXOTIC_HOOK and "drop" in erb.exotic_hook():                       # P13: the dropped channels enter the candidates
        key += hashlib.sha1(np.ascontiguousarray(erb.exotic_hook()["drop"], np.int64).tobytes()).hexdigest()
    if f.exists():
        z = np.load(f, allow_pickle=False)
        if str(z["key"]) == key:
            print(f"[cache] {f}", flush=True)
            paths = [str(p) for p in z["paths"]]
            return {p: {q[len(p) + 1:]: z[q] for q in z.files if q.startswith(p + ":")} for p in paths + ["split"]}
    t0 = time.monotonic()
    with mp.get_context("fork").Pool(NPROC) as workers:
        res = workers.map(pair_job, jobs, chunksize=1)
    paths = list(res[0])
    out = {p: {q: np.concatenate([r[p][q] for r in res]) for q in res[0][p] if q not in ("ncl", "n_events")} for p in paths}
    out["split"] = dict(run=np.array([int(np.load(j[0])["run"]) for j in jobs]), split=np.array([j[1] for j in jobs]),
                        n_events=np.array([r[paths[0]]["n_events"] for r in res]),
                        **{f"ncl_{p}": np.array([r[p]["ncl"] for r in res]) for p in paths})
    Path(cache).mkdir(parents=True, exist_ok=True)
    np.savez_compressed(f, key=key, paths=np.array(paths), **{f"{p}:{q}": v for p, d in out.items() for q, v in d.items()})
    print(f"[collect] {len(jobs)} splits, {time.monotonic() - t0:.0f} s -> {f}", flush=True)
    return out


def exotic_pairs(D, scan_dir, geo):
    """P9 (erb.EXOTIC_HOOK) on the candidates D (per path, the two clusters of a pair consecutive): member corrections as erb.exotic_apply, and a
    pair is dropped when either cluster is vetoed (the seed-radius cut is not applied).  Returns D in the same form plus evn (EVIO event number)."""
    evn_of = {}
    for fn in sorted(Path(scan_dir).glob("*.npz")):
        with np.load(fn, allow_pickle=False) as z:
            evn_of[int(z["run"])] = z["ev_number"]
    out = {"split": D["split"]}
    for p, d in D.items():
        if p == "split":
            continue
        n = d["e"].size
        ptr = np.r_[0, np.cumsum(d["cnt"])]
        mi = erb.objects([d["mi"][ptr[i]:ptr[i + 1]] for i in range(n)])
        me = erb.objects([d["me"][ptr[i]:ptr[i + 1]] for i in range(n)])
        e, er, alive = d["e"].astype(np.float64), d["er"].astype(np.float64), np.ones(n, bool)
        evn = np.zeros(n, np.int64)
        for run in np.unique(d["run"]).tolist():
            m = np.flatnonzero(d["run"] == run)
            evn[m] = evn_of[run][d["ev"][m].astype(np.int64)]
            r = erb.exotic_apply({"x": dict(ev=d["ev"][m], seed=d["seed"][m], mi=mi[m], me=me[m], e=e[m], er=er[m], idx=m)}, run, evn_of[run], geo, rmin=False)["x"]
            alive[m] = np.isin(m, r["idx"])
            e[r["idx"]], er[r["idx"]] = r["e"], r["er"]
            me[r["idx"]] = r["me"]
        keep = np.repeat(alive[0::2] & alive[1::2], 2)
        q = {k: v[keep] for k, v in d.items() if k not in ("mi", "me") and np.ndim(v) and len(v) == n}
        q |= dict(e=e[keep].astype(d["e"].dtype), er=er[keep].astype(d["er"].dtype), evn=evn[keep],
                  mi=np.concatenate(list(mi[keep])) if keep.any() else d["mi"][:0], me=np.concatenate(list(me[keep])) if keep.any() else d["me"][:0])
        out[p] = q | {k: v for k, v in d.items() if k not in q and k not in ("mi", "me")}
        print(f"[exotic] {p}: {int((~keep).sum()) // 2} of {n // 2} candidate pairs vetoed", flush=True)
    return out


# selection and geometry ======================================================
def kin_centre_z(x1, y1, x2, y2, eb):
    """Beam centre (x0, y0) = least-squares point nearest to the lines through the two clusters of each pair, and the vertex distance
    z = median of sqrt((gamma + 1) / 2 R1 R2) about it (tan th1 tan th2 = 2 / (gamma + 1), tan th = R / z)."""
    ux, uy = x2 - x1, y2 - y1
    ln = np.maximum(np.hypot(ux, uy), 1e-9)
    nx, ny = -uy / ln, ux / ln
    x0, y0 = np.linalg.lstsq(np.c_[nx, ny], nx * x1 + ny * y1, rcond=None)[0]
    g = eb / erb.M_ELECTRON
    z = float(np.median(np.sqrt(0.5 * (g + 1.0) * np.hypot(x1 - x0, y1 - y0) * np.hypot(x2 - x0, y2 - y0))))
    return float(x0), float(y0), z


def select(d, geo, split_tab, path, cl_min, r_iso):
    """Cut flow and kinematics of the candidates d of one path (pairs = consecutive clusters).  Returns dict: per pair run, split, dphi
    (about the run's beam centre), esum = (E1 + E2) / Ebeam (calibration-file energies), prod = tan th1 tan th2 (gamma + 1) / 2; per
    cluster th (deg), t = E_pred (MeV); dmin = distance from the pair to the nearest other cluster above CL_STORE kE (mm, NaN without
    one); flow = [(cut, events after it)]; masks = {cut: bool per pair} (cumulative); sel = final mask; geo_run = {run: (x0, y0, z,
    pairs used)}.  cl_min (MeV or None): global third-cluster threshold (exactly two clusters above it); r_iso (mm, 0 = off):
    isolation, no other cluster above CL_STORE kE within r_iso of either pair cluster."""
    eb, kE = geo["ebeam"], geo["ebeam"] / erb.E_REF
    a, b = slice(0, None, 2), slice(1, None, 2)
    reg = erb.module_regions(geo)
    s1, s2 = d["seed"][a].astype(np.int64), d["seed"][b].astype(np.int64)
    fid_mod = (reg >= 0) & (reg <= 4)                      # PbWO4 seed outside H (central hole) and T (PbWO4 / lead-glass transition)
    run = d["run"][a]
    two = ((d["e3"][a] <= cl_min * kE) & (d["e"][b] > cl_min * kE) & (d["e"][a] > cl_min * kE) if cl_min is not None
           else np.ones(run.size, bool))
    dmin = np.fmin(d["dn"][a], d["dn"][b])
    iso = ~(dmin < r_iso) if r_iso > 0 else np.ones(run.size, bool)       # NaN (no other cluster) passes
    fid = fid_mod[s1] & fid_mod[s2] & (d["seed_e"][a] >= erb.SEED_MIN * kE) & (d["seed_e"][b] >= erb.SEED_MIN * kE)
    esum = (d["er"][a] + d["er"][b]) / eb
    pre = two & iso & fid & (np.abs(esum - 1.0) < ESUM_TOL)
    g = eb / erb.M_ELECTRON
    dphi, prod, th = np.full(run.size, np.nan), np.full(run.size, np.nan), np.full(d["e"].size, np.nan)
    geo_run = {}
    for r in np.unique(run).tolist():
        rr = run == r
        x0, y0, z = 0.0, 0.0, geo["hz"]
        for it in range(GEO_ITER + 1):
            x1, y1, x2, y2 = d["x"][a] - x0, d["y"][a] - y0, d["x"][b] - x0, d["y"][b] - y0
            dp = dphi_dev(x1, y1, x2, y2)
            pr = 0.5 * (g + 1.0) * np.hypot(x1, y1) * np.hypot(x2, y2) / z ** 2
            m = pre & rr & (dp < DPHI_CUT) & ((np.abs(pr - 1.0) < PROD_CUT) if it else True)
            if it == GEO_ITER:
                break
            x0, y0, z = kin_centre_z(d["x"][a][m], d["y"][a][m], d["x"][b][m], d["y"][b][m], eb)
        geo_run[r] = (x0, y0, z, int(m.sum()))
        dphi[rr], prod[rr] = dp[rr], pr[rr]
        rc = np.repeat(rr, 2)
        th[rc] = np.degrees(np.arctan(np.hypot(d["x"][rc] - x0, d["y"][rc] - y0) / z))
    t = erb.e_moller(th, eb)
    masks = {}
    ncl = split_tab[f"ncl_{path}"]
    flow = [("all events", int(split_tab["n_events"].sum())),
            (f"hit sum > {1 - ESUM_TOL:g} Ebeam (prefilter)", int(ncl[:, 0].sum()))]
    if cl_min is not None:
        flow.append((f"exactly two clusters above {cl_min:g} MeV", int(ncl[:, CL_COUNT.index(cl_min), 2].sum())))
    m = two.copy()
    cuts = [(f"two largest clusters PbWO4, loose: |dphi - 180| < {DPHI_STORE:g} deg about (0, 0), |E1 + E2 - Ebeam| < {ESUM_STORE:g} Ebeam", m)]
    if r_iso > 0:
        cuts.append((f"isolation: no other cluster above {CL_STORE:g} MeV within {r_iso:g} mm of either", iso))
    for name, c in cuts + [
                    ("fiducial: both seeds outside H and T, seed >= 200 MeV", fid),
                    (f"energy sum: |E1 + E2 - Ebeam| < {ESUM_TOL:g} Ebeam", np.abs(esum - 1.0) < ESUM_TOL),
                    (f"coplanar: |dphi - 180| < {DPHI_CUT:g} deg about the beam centre", dphi < DPHI_CUT),
                    (f"Moller angles: |tan th1 tan th2 (gamma + 1) / 2 - 1| < {PROD_CUT:g}", np.abs(prod - 1.0) < PROD_CUT)]:
        m = m & c
        masks[name] = m.copy()
        flow.append((name, int(m.sum())))
    return dict(run=run, split=d["split"][a], dphi=dphi, esum=esum, prod=prod, th=th, t=t, flow=flow, masks=masks, sel=m,
                pre=pre, two=two, iso=iso, fid=fid, dmin=dmin, geo_run=geo_run)


# gains =======================================================================
def sub(d, m):
    """Clusters of d selected by the bool mask m (per cluster), with their members."""
    mm = np.repeat(m, d["cnt"])
    return {q: v[m] for q, v in d.items() if q not in ("mi", "me") and np.ndim(v) and v.shape[0] == m.size} | dict(mi=d["mi"][mm], me=d["me"][mm])


def ecal(c, g):
    """sum_i g_i E_i per cluster of c."""
    return np.bincount(np.repeat(np.arange(c["cnt"].size), c["cnt"]), weights=c["me"] * g[c["mi"]], minlength=c["cnt"].size)


def fit_scale(T):
    """kE, kS, bin width and fit range of a spectrum around T MeV, the benchmark's convention at energy T."""
    kE = T / erb.E_REF
    kS = math.sqrt(kE)
    return kE, kS, erb.BIN * kS, T - erb.FIT_HALF * kE, T + erb.FIT_HALF * kE


def mwindows(x, seed, T):
    """Global peak (CB, else fixed Gaussian) of x = E / E_pred x T and per seed module its window mu_M +- WINDOW_NSIGMA sigma_M (module
    Crystal Ball on the prefit range with >= MIN_FIT_EVENTS clusters, else median and the global sigma), kept with >= MIN_PER_MODULE
    clusters inside (module_set without the elastic-Moller gap).  Returns (windows, sigma_glob, mu_global)."""
    kE, kS, bw, lo, hi = fit_scale(T)
    c, y = erb.hist(x, lo, hi, bw)
    mu_g, sig_g = erb.gauss_fixed(c, y, lo, hi, bw, kS)
    mu_cb, sig_cb = erb.cb_single(c, y, lo, hi, bw, kS, mu_g) if math.isfinite(mu_g) else erb.NAN2
    sg, mu = (sig_cb, mu_cb) if math.isfinite(sig_cb) and math.isfinite(mu_cb) else (sig_g, mu_g)
    win = {}
    if not (math.isfinite(sg) and math.isfinite(mu)):
        return win, sg, mu
    plo, phi = mu * (1.0 + erb.PREFIT_RANGE[0] / kS), mu * (1.0 + erb.PREFIT_RANGE[1] / kS)
    order = np.argsort(seed, kind="stable")
    ss, xs = seed[order], x[order]
    bnd = np.flatnonzero(np.r_[True, ss[1:] != ss[:-1], True])
    for i0, i1 in zip(bnd[:-1], bnd[1:]):
        vals = xs[i0:i1]
        sel = vals[(vals >= plo) & (vals <= phi)]
        if sel.size == 0:
            continue
        f = erb.module_cb(vals, plo, phi, mu, bw, sg, kS) if sel.size >= erb.MIN_FIT_EVENTS else None
        mu_m, sig_m = f if f else (float(np.median(sel)), sg)
        n_in = int(((vals >= mu_m - erb.WINDOW_NSIGMA * sig_m) & (vals <= mu_m + erb.WINDOW_NSIGMA * sig_m)).sum())
        if n_in >= erb.MIN_PER_MODULE:
            win[int(ss[i0])] = (mu_m - erb.WINDOW_NSIGMA * sig_m, mu_m + erb.WINDOW_NSIGMA * sig_m)
    return win, sg, mu


def solve_moller(c, t, n_all, T):
    """Per-member gains so that sum_i g_i E_i of each training cluster matches its E_pred t: solve_gains with per-cluster targets.
    Training clusters: x = E / t x T inside their seed module's window (mwindows; modules without one use the global window mu_global
    +- WINDOW_NSIGMA sigma_glob), windows refitted in the first GAIN_N_REFIT passes; prior lambda (g - g0)^2, g0 = T / mu_global,
    lambda = MIN_PER_MODULE (GAIN_SEED_FRAC T)^2 / sigma_glob^2 of the g = 1 pass; converged when |dg| < max(GAIN_DG_TOL,
    GAIN_TOL_SIGMA sigma_g) for every free module (prior fraction < GAIN_FREE_FRAC).  Returns dict(g, g0, rows, iterations,
    converged, n_free, hist, own = modules with their own window) or None."""
    seed = c["seed"].astype(np.int64)
    evid = np.repeat(np.arange(seed.size), c["cnt"])
    mi, me = c["mi"].astype(np.int64), c["me"].astype(float)
    win0, sg0, mu0 = mwindows(ecal(c, np.ones(n_all)) / t * T, seed, T)
    if not (math.isfinite(sg0) and math.isfinite(mu0) and sg0 > 0):
        return None
    g0, w = T / mu0, 1.0 / sg0 ** 2
    lam = erb.MIN_PER_MODULE * (erb.GAIN_SEED_FRAC * T) ** 2 * w
    g, converged, hist = np.ones(n_all), False, []
    for it in range(1, erb.GAIN_MAX_ITER + 1):
        x = ecal(c, g) / t * T
        if it <= erb.GAIN_N_REFIT:
            win, sg, mu = (win0, sg0, mu0) if it == 1 else mwindows(x, seed, T)
            if not (math.isfinite(sg) and math.isfinite(mu)):
                return None
            lo, hi = np.full(n_all, mu - erb.WINDOW_NSIGMA * sg), np.full(n_all, mu + erb.WINDOW_NSIGMA * sg)
            for m_, (a_, b_) in win.items():
                lo[m_], hi[m_] = a_, b_
        sel = (x >= lo[seed]) & (x <= hi[seed])
        if not sel.any():
            return None
        msel = sel[evid]
        uniq, inv = np.unique(mi[msel], return_inverse=True)
        Nm, b = erb.normal_equations(np.r_[0, np.cumsum(c["cnt"][sel])], inv, me[msel], t[sel], uniq.size)
        Nm, b = Nm * w, b * w
        D = np.diag(Nm).copy()
        A = Nm + lam * np.eye(uniq.size)
        gsol = np.linalg.solve(A, b + lam * g0)
        sigma_g = np.sqrt(np.diag(np.linalg.inv(A))) * (sg / sg0)
        pf = lam / (D + lam)
        gnew = np.full(n_all, g0)
        gnew[uniq] = gsol
        free = pf < erb.GAIN_FREE_FRAC
        dg = np.abs(gnew[uniq] - g[uniq])[free]
        converged = bool(free.any() and (dg < np.maximum(erb.GAIN_DG_TOL, erb.GAIN_TOL_SIGMA * sigma_g[free])).all())
        hist.append(dict(it=it, n_sel=int(sel.sum()), n_mod=int(uniq.size), n_free=int(free.sum()), n_win=len(win),
                         max_dg=float(dg.max()) if dg.size else float("nan"), med_dg=float(np.median(dg)) if dg.size else float("nan"),
                         mu=mu / T, sg=sg / mu))
        g = gnew
        if converged:
            break
    n_seed = np.bincount(np.searchsorted(uniq, seed[sel]), minlength=uniq.size)
    n_member = np.bincount(inv, minlength=uniq.size)
    rows = [dict(module=int(m_), g=gsol[j], sigma_g=sigma_g[j], n_seed=int(n_seed[j]), n_member=int(n_member[j]), prior_frac=pf[j])
            for j, m_ in enumerate(uniq.tolist())]
    return dict(g=g, g0=g0, rows=rows, iterations=it, converged=converged, n_free=int(free.sum()), hist=hist, own=set(win))


# fits and evaluation =========================================================
def h2_curve(e_mev):
    """H2 resolution curve sigma/E (%) at E (MeV)."""
    a, b, c = H2_CURVE
    e = np.asarray(e_mev, float) / 1000.0
    return np.sqrt(a ** 2 + b ** 2 / e + c ** 2 / e ** 2)


def cb_fit(x, T):
    """Crystal Ball of the values x (MeV) around T exactly as the benchmark's cb_single (fixed-window Gaussian start, window
    [mu0 - win, mu0 + win / 2], same start values and bounds; kE = T / E_REF): dict(mu, sig, res = 100 sig / mu, res_err from the
    covariance, n in range, p, c, y)."""
    nan = float("nan")
    kE, kS, bw, lo, hi = fit_scale(T)
    c, y = erb.hist(x, lo, hi, bw)
    out = dict(mu=nan, sig=nan, res=nan, res_err=nan, n=int(y.sum()), p=None, c=c, y=y)
    mu0, _ = erb.gauss_fixed(c, y, lo, hi, bw, kS)
    if not math.isfinite(mu0):
        return out
    win = max(400.0 * kS, 12.0 * bw)
    m = (c > mu0 - win) & (c < mu0 + 0.5 * win)
    if y[m].sum() < 50:
        return out
    p, cov = erb.fit_curve(erb.crystal_ball, c[m], y[m], [float(y[m].max()), mu0, max(60.0 * kS, 3 * bw), 1.5, 2.0],
                           [0, lo - 200.0 * kS, 1.0, 0.1, 1.01], [np.inf, hi + 200.0 * kS, win, 20.0, 50.0], maxfev=10000)
    if p is None:
        return out
    mu, sig = float(p[1]), abs(float(p[2]))
    J = np.array([-sig / mu ** 2, np.sign(p[2]) / mu])
    var = float(J @ np.asarray(cov)[1:3, 1:3] @ J) if np.ndim(cov) == 2 else nan
    return out | dict(mu=mu, sig=sig, res=100.0 * sig / mu, res_err=100.0 * math.sqrt(var) if var >= 0 else nan, p=p)


def corrected(c, g, pc, geo):
    """Corrected cluster energy sum g_i E_i / R(u, v) (gains, then the borrowed position correction)."""
    return ecal(c, g) / erb.pc_eval(pc, c, geo)


def prod_width(prod):
    """Core sigma of tan th1 tan th2 (gamma + 1) / 2: iterated Gaussian within 2 sigma (benchmark gauss_2sigma) on 1000 x prod."""
    c, y = erb.hist(1000.0 * prod, 700.0, 1300.0, 5.0)
    mu, sg = erb.gauss_2sigma(c, y, *erb.gauss_fixed(c, y, 700.0, 1300.0, 5.0, 0.5))
    return sg / 1000.0


def bootstrap_se(groups, fn, n_boot=100, seed=1):
    """Standard deviation of fn over resamples (with replacement) of the groups (one array of values per split); NaN fits skipped."""
    rng, v = np.random.default_rng(seed), []
    for _ in range(n_boot):
        pick = rng.integers(0, len(groups), len(groups))
        r = fn(np.concatenate([groups[i] for i in pick]))
        if math.isfinite(r):
            v.append(r)
    return float(np.std(v, ddof=1)) if len(v) > 1 else float("nan")


def results(x, T, keys, ref_res=float("nan"), extra=None):
    """Per-split CB of x around T (mean +- standard error over the splits keys) and the pooled fit (+- bootstrap over the splits)."""
    splits = sorted(set(keys))
    groups = [x[keys == k] for k in splits]
    per = [cb_fit(gx, T)["res"] for gx in groups]
    pooled = cb_fit(x, T)
    ms = erb.mean_se(per)
    return dict(E=T, n=int(x.size), n_splits=len(splits), n_fit_splits=int(np.isfinite(per).sum()), res=ms[0], res_se=ms[1],
                res_pooled=pooled["res"], res_pooled_se=bootstrap_se(groups, lambda v: cb_fit(v, T)["res"]),
                res_pooled_fit_err=pooled["res_err"], peak=100.0 * (pooled["mu"] / T - 1.0), h2=ref_res, fit=pooled, per=per) | (extra or {})


# report ======================================================================
def analyse(D, geo, sp, paths, gp, cl_min, r_iso, combine_runs, pc, parity, quick=False):
    """Selection, gain solves and fits of one configuration.  quick: only the gate path, out of sample, the sum and the symmetric pairs
    (comparison table).  Returns dict(S, C, solves, rows, keep, T, Ts, sig_p)."""
    n_all, eb = geo["x"].size, geo["ebeam"]
    S = {p: select(D[p], geo, sp, p, cl_min, r_iso) for p in paths}
    C = {}
    for p in paths:
        cm = np.repeat(S[p]["sel"], 2)
        C[p] = sub(D[p], cm) | dict(t=S[p]["t"][cm], th=S[p]["th"][cm])
        C[p]["par"] = np.array([parity[r, s] for r, s in zip(C[p]["run"].tolist(), C[p]["split"].tolist())], np.int64)
        C[p]["key"] = C[p]["run"].astype(np.int64) * 10000 + C[p]["split"]
    tr = C["legacy"]
    T = float(np.mean(tr["t"]))
    groups = {"all": np.ones(tr["t"].size, bool)} if combine_runs else {str(r): tr["run"] == r for r in np.unique(tr["run"]).tolist()}
    solves = {}                                     # "<run or all>:<in | 0 | 1>"; per run as the ep chain unless combine_runs
    for gname, gm in groups.items():
        for k in ("in", 0, 1):
            m = gm & ((tr["par"] == k) if k != "in" else True)
            solves[f"{gname}:{k}"] = solve_moller(sub(tr, m) | dict(t=tr["t"][m]), tr["t"][m], n_all, T)
    if any(v is None for v in solves.values()):
        sys.exit(f"[ERROR] Moller gain solve failed (cl_min {cl_min}, r_iso {r_iso:g})")
    if not quick:
        for k, v in solves.items():
            print(f"[gains] solve {k}: {len(v['rows'])} modules ({v['n_free']} free, {len(v['own'])} with own window), g0 {v['g0']:.4f}, "
                  f"{v['iterations']} passes, converged {v['converged']}", flush=True)
    sig_p = prod_width(S[gp]["prod"][S[gp]["sel"]])
    th_s = math.degrees(math.atan(math.sqrt(2.0 / (eb / erb.M_ELECTRON + 1.0))))
    Ts = eb + erb.M_ELECTRON - 2.0 * float(erb.energy_loss(th_s))      # E1 + E2 expected: Ebeam + m_e - losses (symmetric angle)
    edges = np.quantile(C[gp]["t"], np.linspace(0, 1, N_EBINS + 1))
    rows, keep = [], {}
    for p in ([gp] if quick else paths):
        c = C[p]
        gname = np.full(c["t"].size, "all") if combine_runs else c["run"].astype(str)
        e_in, e_oos = np.zeros(c["t"].size), np.zeros(c["t"].size)
        for gn in np.unique(gname).tolist():
            for h in (0, 1):
                m = (gname == gn) & (c["par"] == h)
                cm_ = sub(c, m)
                e_in[m] = corrected(cm_, solves[f"{gn}:in"]["g"], pc, geo)
                e_oos[m] = corrected(cm_, solves[f"{gn}:{1 - h}"]["g"], pc, geo)
        for prot, e in ((("out-of-sample", e_oos),) if quick else (("in-sample", e_in), ("out-of-sample", e_oos))):
            r = e / c["t"]
            res = {} if quick else {"cluster": results(r * T, T, c["key"], float(h2_curve(T)))}
            es = e[0::2] + e[1::2]                                       # the sum against a constant: angle independent
            s1, s2 = c["t"][0::2] * h2_curve(c["t"][0::2]) / 100.0, c["t"][1::2] * h2_curve(c["t"][1::2]) / 100.0
            res["sum"] = results(es, Ts, c["key"][0::2], 100.0 * math.sqrt(float(np.mean(s1 ** 2 + s2 ** 2))) / Ts)
            sym = np.abs(c["t"][0::2] - c["t"][1::2]) < SYM_MAX * (c["t"][0::2] + c["t"][1::2])   # sigma_cluster / E = sqrt(2) sigma_sum / E_sum
            rs = results(es[sym], Ts, c["key"][0::2][sym], 100.0 * math.sqrt(float(np.mean(s1[sym] ** 2 + s2[sym] ** 2))) / Ts)
            rs |= dict(res_cluster_equiv=math.sqrt(2.0) * rs["res"], res_cluster_equiv_se=math.sqrt(2.0) * rs["res_se"],
                       res_cluster_equiv_pooled=math.sqrt(2.0) * rs["res_pooled"],
                       res_cluster_equiv_pooled_se=math.sqrt(2.0) * rs["res_pooled_se"],
                       E_cluster=float(np.mean(np.r_[c["t"][0::2][sym], c["t"][1::2][sym]])))
            res[SYM_Q] = rs
            if not quick:
                ib = np.clip(np.searchsorted(edges, c["t"], side="right") - 1, 0, N_EBINS - 1)
                for b_ in range(N_EBINS):
                    mb = ib == b_
                    Tb = float(np.mean(c["t"][mb]))
                    res[f"cluster, E_pred bin {b_ + 1}"] = results(r[mb] * Tb, Tb, c["key"][mb], float(h2_curve(Tb)),
                                                                   dict(lo=edges[b_], hi=edges[b_ + 1]))
            for q, v in res.items():
                rows.append(dict(path=p, protocol=prot, quantity=q) | v)
            keep[p, prot] = dict(e=e, r=r, res=res)
    return dict(S=S, C=C, solves=solves, rows=rows, keep=keep, T=T, Ts=Ts, sig_p=sig_p)


def config_name(cl_min, r_iso):
    """Label of a selection configuration."""
    parts = ([f"isolation {r_iso:g} mm"] if r_iso > 0 else []) + ([f"exactly two above {cl_min:g} MeV"] if cl_min is not None else [])
    return ", ".join(parts) or "no third-cluster cut"


def run_breakdown(res, geo, gp, calib_of):
    """Per run (gate path, out of sample): pairs, sum sigma/E per split and pooled, and the effective gain 1 / (calibration factor x
    Moller gain, in-sample solve of the run) of the inner layers relative to the first run (layer = ring around the beam hole: 1-2 =
    the central block H, the HV-raised innermost two layers; 3 = the next ring), also the calibration-file part alone."""
    c, e = res["C"][gp], res["keep"][gp, "out-of-sample"]["e"]
    es, run = e[0::2] + e[1::2], c["run"][0::2]
    pw, row, col = geo["pwo4"], geo["row"].astype(float), geo["col"].astype(float)
    layer = np.where(pw, np.maximum(np.abs(row - 16.5), np.abs(col - 16.5)) - 0.5, -1)
    runs = sorted(calib_of)
    fac = {}
    for r in runs:
        f = {x["name"]: float(x["factor"]) for x in json.load(open(calib_of[r]))}
        fac[r] = np.array([f.get(str(n), np.nan) for n in geo["name"]])
    out = []
    for r in runs:
        m = run == r
        rr = results(es[m], res["Ts"], c["key"][0::2][m]) if m.any() else {}
        geff = lambda q: 1.0 / (fac[q] * res["solves"][f"{q}:in"]["g"]) if f"{q}:in" in res["solves"] else np.full(pw.size, np.nan)  # noqa: E731
        row_ = dict(run=r, date=RUN_DATE.get(r, ""), calib=Path(calib_of[r]).name, pairs=int(m.sum()), n_splits=rr.get("n_splits", 0),
                    res=rr.get("res", np.nan), res_se=rr.get("res_se", np.nan), res_pooled=rr.get("res_pooled", np.nan),
                    res_pooled_se=rr.get("res_pooled_se", np.nan))
        for lname, lm in (("L12", (layer >= 1) & (layer <= 2)), ("L3", layer == 3)):
            ok = lm & np.isfinite(fac[runs[0]]) & (fac[runs[0]] > 0) & np.isfinite(fac[r]) & (fac[r] > 0)
            row_[f"gain_{lname}"] = float(np.mean(geff(r)[ok] / geff(runs[0])[ok]))
            row_[f"calib_{lname}"] = float(np.mean(fac[runs[0]][ok] / fac[r][ok]))
            row_[f"n_{lname}"] = int(ok.sum())
        out.append(row_)
    return out


def report(scan_dir, gate=None, out_root=OUT_ROOT, ref_root=REF_ROOT, cache=CACHE_ROOT, cl_min=CL_MIN, r_iso=R_ISO, combine_runs=False):
    t0 = time.monotonic()
    if cl_min is not None and cl_min not in CL_COUNT:
        sys.exit(f"[ERROR] --cl-min {cl_min:g}: one of {', '.join(f'{c:g}' for c in CL_COUNT)} MeV, or none")
    resolve_report_globals()
    gate = float(gate if gate is not None else erb.CONTROL_GATE[RT])
    D = collect(scan_dir, gate, cache)
    sp = D["split"]
    paths, gp = ["legacy", f"g{gate:g}"], f"g{gate:g}"
    calib_of = {}
    for fn in sorted(Path(scan_dir).glob("*.npz")):
        with np.load(fn, allow_pickle=False) as z:
            if int(z["run"]) not in erb.EXCLUDE_RUNS:
                calib_of[int(z["run"])] = str(z["calib"])
    if erb.EXOTIC_HOOK:                            # P9: the runs of the hook's splits; veto and member corrections on the candidates
        calib_of = {r: c for r, c in calib_of.items() if r in set(sp["run"].tolist())}
    r0 = min(calib_of)
    _, geo = erb.clusterers(r0, calib_of[r0], (gate,))
    if erb.EXOTIC_HOOK:
        D = exotic_pairs(D, scan_dir, geo)
    pc_src = Path(ref_root) / erb.PC_SETTING[RT][0]
    pc = erb.read_pc(pc_src, geo["name"])[(0, 0)]
    suffix = ("" if cl_min == CL_MIN else f"_clmin{cl_min:g}" if cl_min is not None else "_clminnone") + \
             ("" if r_iso == R_ISO else f"_riso{r_iso:g}") + ("_combined" if combine_runs else "")
    out = Path(out_root) / f"{RT}_moller{suffix}"
    out.mkdir(parents=True, exist_ok=True)
    order = {}
    for r, s in zip(sp["run"].tolist(), sp["split"].tolist()):
        order.setdefault(r, []).append(s)
    parity = {(r, s): sorted(v).index(s) % 2 for r, v in order.items() for s in v}
    main = analyse(D, geo, sp, paths, gp, cl_min, r_iso, combine_runs, pc, parity)
    print(f"[angles] core sigma of tan th1 tan th2 (gamma + 1) / 2: {main['sig_p']:.4f}", flush=True)
    comp = []                                      # comparison of the third-cluster rules (gate path, out of sample)
    for cm_, ri in COMPARE:
        res = main if (cm_, ri) == (cl_min, r_iso) else analyse(D, geo, sp, paths, gp, cm_, ri, combine_runs, pc, parity, quick=True)
        get = lambda q: next(x for x in res["rows"] if x["path"] == gp and x["protocol"] == "out-of-sample" and x["quantity"] == q)  # noqa: E731
        own = sum(len(v["own"]) for k, v in res["solves"].items() if k.endswith(":in"))
        comp.append(dict(config=config_name(cm_, ri), cl_min=cm_, r_iso=ri, main=(cm_, ri) == (cl_min, r_iso), pairs=int(res["S"][gp]["sel"].sum()),
                         own=own, sum=get("sum"), sym=get(SYM_Q)))
        print(f"[compare] {comp[-1]['config']}: {comp[-1]['pairs']} pairs, sum {comp[-1]['sum']['res']:.3f} +- {comp[-1]['sum']['res_se']:.3f} "
              f"(pooled {comp[-1]['sum']['res_pooled']:.3f})", flush=True)
    if erb.EXOTIC_HOOK:                            # P9: the selected gate-path clusters (pairs consecutive), energies in and out of sample
        c = main["C"][gp]
        np.savez_compressed(out / "p9_moller.npz", T=main["T"], Ts=main["Ts"], **{q: c[q] for q in ("run", "split", "evn", "seed", "x", "y", "t")},
                            e_in=main["keep"][gp, "in-sample"]["e"], e_oos=main["keep"][gp, "out-of-sample"]["e"])
    per_run = run_breakdown(main, geo, gp, calib_of)
    S0 = select(D[gp], geo, sp, gp, None, 0.0)       # isolation study: every other cut, no third-cluster cut
    iso_study = dict(dmin=S0["dmin"][S0["sel"]], esum=(D[gp]["er"][0::2] + D[gp]["er"][1::2])[S0["sel"]] / geo["ebeam"])
    write_outputs(out, main["rows"], main["keep"], main["S"], main["C"], D, main["solves"], geo, paths, gp, main["T"], cl_min, r_iso, pc_src,
                  main["sig_p"], comp, per_run, iso_study)
    print(f"[done] {RT} Moller: {time.monotonic() - t0:.0f} s -> {out}", flush=True)


SUM_COLS = ("path", "protocol", "quantity", "E", "n", "n_splits", "n_fit_splits", "res", "res_se", "res_pooled", "res_pooled_se",
            "res_pooled_fit_err", "peak", "h2", "lo", "hi", "res_cluster_equiv", "res_cluster_equiv_se", "res_cluster_equiv_pooled",
            "res_cluster_equiv_pooled_se", "E_cluster")


def write_outputs(out, rows, keep, S, C, D, solves, geo, paths, gp, T, cl_min, r_iso, pc_src, sig_p, comp, per_run, iso_study):
    """summary.tsv, cutflow.tsv, geometry.tsv, gains.tsv, compare.tsv, runs.tsv and report.pdf."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    nan = float("nan")
    erb.write_tsv(out / "summary.tsv", SUM_COLS, [[r.get(q, nan) for q in SUM_COLS] for r in rows],
                  f"# {RT} Moller calibration: res = per-split Crystal Ball sigma / peak (%) mean, res_se its standard error over the splits; "
                  f"res_pooled = all splits in one fit, res_pooled_se bootstrap over the splits; E = mean E_pred (MeV; sum: Ebeam + m_e - "
                  f"2 x energy loss at the symmetric angle; the sum fitted is E1 + E2 itself); peak = CB peak / E - 1 (%); h2 = H2 curve {H2_CURVE} "
                  f"at E (sum: quadrature of the two clusters' curve values / E); core sigma of tan th1 tan th2 (gamma + 1) / 2: {sig_p:.4f}; lo, hi = E_pred bin edges (MeV); res_cluster_equiv(_se) = sqrt(2) res of "
                  f"the symmetric-pair sum = single-cluster sigma / E at E_cluster (MeV); selection: {config_name(cl_min, r_iso)}")
    erb.write_tsv(out / "cutflow.tsv", ("path", "step", "cut", "events"),
                  [[p, i, name, n] for p in paths for i, (name, n) in enumerate(S[p]["flow"])])
    erb.write_tsv(out / "geometry.tsv", ("path", "run", "x0_mm", "y0_mm", "z_mm", "hycal_z_mm", "pairs"),
                  [[p, r, *v[:3], geo["hz"], v[3]] for p in paths for r, v in S[p]["geo_run"].items()])
    erb.write_tsv(out / "gains.tsv", ("solve", "module", "g", "sigma_g", "n_seed", "n_member", "prior_frac", "own_window"),
                  [[k, str(geo["name"][x["module"]]), x["g"], x["sigma_g"], x["n_seed"], x["n_member"], x["prior_frac"], int(x["module"] in v["own"])]
                   for k, v in solves.items() for x in v["rows"]],
                  "# solve <run or all>:<in | 0 | 1>: in = all splits (in sample), 0 / 1 = even / odd split ordinal (out-of-sample gains of the other half); g0: "
                  + ", ".join(f"{k} {v['g0']:.5f}" for k, v in solves.items()))
    erb.write_tsv(out / "compare.tsv", ("config", "default", "pairs", "own_windows", "sum_res", "sum_res_se", "sum_n_fit_splits", "sum_res_pooled",
                                        "sum_res_pooled_se", "cluster_equiv", "cluster_equiv_se", "cluster_equiv_pooled", "cluster_equiv_pooled_se"),
                  [[x["config"], int(x["main"]), x["pairs"], x["own"], x["sum"]["res"], x["sum"]["res_se"], x["sum"]["n_fit_splits"],
                    x["sum"]["res_pooled"], x["sum"]["res_pooled_se"], x["sym"]["res_cluster_equiv"], x["sym"]["res_cluster_equiv_se"],
                    x["sym"]["res_cluster_equiv_pooled"], x["sym"]["res_cluster_equiv_pooled_se"]] for x in comp],
                  f"# third-cluster rules, gate path {gp}, out of sample, per-run gains: E1 + E2 sigma/E (%) per split (mean, s.e., fitted splits of 25) "
                  f"and pooled; cluster_equiv = sqrt(2) sigma_sum / E_sum of the symmetric pairs (one cluster at 1.12 GeV); own_windows summed over the runs' in-sample solves")
    rk = ("run", "date", "calib", "pairs", "n_splits", "res", "res_se", "res_pooled", "res_pooled_se", "gain_L12", "calib_L12", "n_L12",
          "gain_L3", "calib_L3", "n_L3")
    erb.write_tsv(out / "runs.tsv", rk, [[x[k] for k in rk] for x in per_run],
                  f"# per run, gate path {gp}, out of sample: E1 + E2 sigma/E (%) per split (mean, s.e.) and pooled; gain_L* = mean over the modules of "
                  f"layer L (rings around the beam hole: L12 = innermost two, the central block H; L3 = next ring) of 1 / (calibration factor x Moller "
                  f"gain of the run's in-sample solve), relative to the first run; calib_L* = the calibration-file part alone")
    plt.rcParams.update({"axes.labelsize": 13, "axes.titlesize": 13, "axes.titleweight": "bold", "xtick.labelsize": 11,
                         "ytick.labelsize": 11, "legend.fontsize": 10, "lines.linewidth": 1.6})
    eb = geo["ebeam"]
    s = S[gp]
    head = f"{RT}: Moller two-cluster calibration (gate path {gp[1:]} ns, out of sample unless stated; {config_name(cl_min, r_iso)})"
    with PdfPages(out / "report.pdf") as pdf:
        # page 1: cut flow, geometry, headline
        fig = plt.figure(figsize=(14, 9.5))
        fig.suptitle(head, fontsize=13, fontweight="bold")
        ax = fig.add_axes([0.03, 0.05, 0.94, 0.85])
        ax.axis("off")
        lines = ["Cut flow (events)                                                                                 legacy        " + gp]
        for (name, n1), (_, n2) in zip(S["legacy"]["flow"], S[gp]["flow"]):
            lines.append(f"  {name:<96s}{n1:>9d}  {n2:>9d}")
        lines += ["", "Beam centre and vertex distance per run (from the coplanar pairs; hycal_z %.1f mm)" % geo["hz"]]
        for p in paths:
            for r, v in S[p]["geo_run"].items():
                lines.append(f"  {p:>7s} run {r}: x0 {v[0]:+.2f} mm  y0 {v[1]:+.2f} mm  z {v[2]:.1f} mm  ({v[3]} pairs)")
        lines += ["", "Resolution, Crystal Ball sigma / peak (%): per-split mean +- s.e. | pooled +- bootstrap s.e. | H2 curve at E"]
        for r in rows:
            if r["protocol"] == "out-of-sample" and not r["quantity"].startswith("cluster, E_pred"):
                lines.append(f"  {r['path']:>7s} {r['quantity']:<26s} E {r['E'] / 1000:5.3f} GeV  n {r['n']:6d}: {r['res']:.3f} +- {r['res_se']:.3f} | "
                             f"{r['res_pooled']:.3f} +- {r['res_pooled_se']:.3f} | {r['h2']:.3f}")
        lines += ["", "Gain solves, <run>:<in = all splits | 0 / 1 = even / odd split ordinal, used out of sample on the other half>:"]
        lines += [f"  {k}: {v['iterations']} passes, converged {v['converged']}, {v['n_free']} free of {len(v['rows'])} modules, "
                  f"{len(v['own'])} with their own window, g0 {v['g0']:.4f}" for k, v in solves.items()]
        lines += [f"Position correction borrowed from {pc_src}",
                  f"Core sigma of tan th1 tan th2 (gamma + 1) / 2: {sig_p:.4f} (E_pred inherits this angle spread)"]
        ax.text(0, 1, "\n".join(lines), family="monospace", fontsize=8.2, va="top")
        pdf.savefig(fig)
        plt.close(fig)
        # page 2: selection distributions (gate path), each with the other cuts applied
        fig, axs = plt.subplots(2, 3, figsize=(16, 10))
        fig.suptitle(f"Selection distributions, path {gp} (each with all other cuts applied)", fontsize=14, fontweight="bold")
        base = s["two"] & s["fid"]
        es_ok, dp_ok, pr_ok = np.abs(s["esum"] - 1) < ESUM_TOL, s["dphi"] < DPHI_CUT, np.abs(s["prod"] - 1) < PROD_CUT
        a = axs[0, 0]
        a.hist(s["dphi"][base & es_ok & pr_ok], bins=np.arange(0, DPHI_STORE + 0.5, 0.5), histtype="step", color="k")
        a.axvline(DPHI_CUT, color="r", ls="--")
        a.set(xlabel="| |phi1 - phi2| - 180 | (deg)", ylabel="pairs", title="coplanarity")
        a = axs[0, 1]
        a.hist(s["esum"][base & dp_ok & pr_ok], bins=np.arange(1 - ESUM_STORE, 1 + ESUM_STORE + 1e-9, 0.005), histtype="step", color="k")
        for v in (1 - ESUM_TOL, 1 + ESUM_TOL):
            a.axvline(v, color="r", ls="--")
        a.set(xlabel="(E1 + E2) / Ebeam, calibration-file energies", ylabel="pairs", title="energy sum")
        a = axs[0, 2]
        a.hist(s["prod"][base & es_ok & dp_ok], bins=np.arange(0.4, 1.6, 0.01), histtype="step", color="k")
        for v in (1 - PROD_CUT, 1 + PROD_CUT):
            a.axvline(v, color="r", ls="--")
        a.set(xlabel="tan th1 tan th2 (gamma + 1) / 2", ylabel="pairs", title="Moller angle relation")
        a = axs[1, 0]
        th = s["th"]
        a.hist2d(th[0::2][base & es_ok & dp_ok], th[1::2][base & es_ok & dp_ok], bins=80, cmap="Greys", cmin=1)
        tt = np.linspace(0.4, 3.0, 200)
        a.plot(tt, moller_partner(tt, eb), "r-", label="tan th1 tan th2 = 2 / (gamma + 1)")
        a.set(xlabel="th1 (deg)", ylabel="th2 (deg)", title="th2 versus th1", xlim=(0.4, 3.0), ylim=(0.4, 3.0))
        a.legend(loc="upper right")
        a = axs[1, 1]
        sel2 = np.repeat(s["sel"], 2)
        a.plot(th[sel2], s["t"][sel2], ",", color="0.4")
        a.plot(tt, erb.e_moller(tt, eb), "r-", label="E_pred(th) = Moller E' - energy loss")
        a.set(xlabel="th (deg)", ylabel="E_pred (MeV)", title="E_pred versus th (selected clusters)", xlim=(0.4, 3.0))
        a.legend(loc="upper right")
        a = axs[1, 2]
        cm_ = cl_min if cl_min is not None else CL_STORE
        ncl = D["split"][f"ncl_{gp}"][:, CL_COUNT.index(cm_)].sum(0)
        a.bar(np.arange(6), ncl, color="0.5")
        a.set_xticks(np.arange(6), ["0", "1", "2", "3", "4", ">=5"])
        a.set(xlabel=f"clusters above {cm_:g} MeV", ylabel="prefiltered events", title="cluster multiplicity")
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        pdf.savefig(fig)
        plt.close(fig)
        # page 3: calibration convergence
        fig, axs = plt.subplots(2, 2, figsize=(14, 10))
        fig.suptitle("Moller gain solve (legacy path, both clusters of each selected pair)", fontsize=14, fontweight="bold")
        for k, v in solves.items():
            if not k.endswith(":in"):              # the half solves are in gains.tsv and on page 1
                continue
            it = [h["it"] for h in v["hist"]]
            axs[0, 0].semilogy(it, [h["max_dg"] for h in v["hist"]], "o-", label=f"{k}: max |dg|")
            axs[0, 0].semilogy(it, [h["med_dg"] for h in v["hist"]], "s--", label=f"{k}: median |dg|")
            axs[0, 1].plot(it, [h["n_sel"] for h in v["hist"]], "o-", label=f"{k}")
        axs[0, 0].axhline(erb.GAIN_DG_TOL, color="r", ls=":", label="GAIN_DG_TOL")
        axs[0, 0].set(xlabel="pass", ylabel="|g_new - g| (free modules)", title="gain change per pass")
        axs[0, 0].legend(fontsize=8)
        axs[0, 1].set(xlabel="pass", ylabel="training clusters in window", title="training sample per pass")
        axs[0, 1].legend()
        ins = [v for k, v in solves.items() if k.endswith(":in")]
        rin = [x for v in ins for x in v["rows"]]
        gg, pf, ns = (np.array([x[q] for x in rin]) for q in ("g", "prior_frac", "n_seed"))
        free = pf < erb.GAIN_FREE_FRAC
        axs[1, 0].hist([gg[free], gg[~free]], bins=60, stacked=True, label=[f"free ({free.sum()})", f"prior-dominated ({(~free).sum()})"])
        for v in ins:
            axs[1, 0].axvline(v["g0"], color="r", ls="--", lw=1)
        axs[1, 0].set(xlabel="gain (in-sample solves; dashed: g0)", ylabel="module rows", title="gains")
        axs[1, 0].legend()
        axs[1, 1].semilogx(np.maximum(ns, 0.5), gg, ".", color="0.3")
        axs[1, 1].set(xlabel="training clusters seeded in the module", ylabel="gain", title="gain versus seed count (in sample)")
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        pdf.savefig(fig)
        plt.close(fig)
        # pages 4-5: Crystal Ball fits (gate path, Moller gains, out of sample; in-sample dashed)
        def draw(a, q, title, xlabel):
            for prot, st in (("out-of-sample", "-"), ("in-sample", ":")):
                r = keep[gp, prot]["res"][q]
                f = r["fit"]
                a.step(f["c"], f["y"], where="mid", color="k" if st == "-" else "0.6", lw=1.2, label=f"{prot}")
                if f["p"] is not None:
                    xx = np.linspace(f["c"][0], f["c"][-1], 600)
                    a.plot(xx, erb.crystal_ball(xx, *f["p"]), "r" + st,
                           label=f"CB {prot}: sigma/E {f['res']:.2f} +- {f['res_err']:.2f} %")
            r = keep[gp, "out-of-sample"]["res"][q]
            a.set(xlabel=xlabel, ylabel="clusters" if q != "sum" else "pairs",
                  title=f"{title}\nper split {r['res']:.2f} +- {r['res_se']:.2f} %, H2 curve {r['h2']:.2f} %")
            a.legend(fontsize=8)
        fig, axs = plt.subplots(1, 2, figsize=(16, 7))
        fig.suptitle(f"Crystal Ball fits, path {gp}, Moller gains + position correction", fontsize=14, fontweight="bold")
        draw(axs[0], "cluster", "(i) per cluster, E / E_pred", f"E / E_pred x {T:.0f} MeV")
        draw(axs[1], "sum", "(ii) E1 + E2 (angle independent)", "E1 + E2 (MeV)")
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        pdf.savefig(fig)
        plt.close(fig)
        fig, axs = plt.subplots(2, 2, figsize=(15, 11))
        fig.suptitle(f"(iii) per E_pred bin, path {gp}", fontsize=14, fontweight="bold")
        for b, a in enumerate(axs.ravel()):
            q = f"cluster, E_pred bin {b + 1}"
            r = keep[gp, "out-of-sample"]["res"][q]
            draw(a, q, f"E_pred {r['lo']:.0f}-{r['hi']:.0f} MeV (mean {r['E']:.0f})", f"E / E_pred x {r['E']:.0f} MeV")
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        pdf.savefig(fig)
        plt.close(fig)
        # page 6: resolution versus energy
        fig, a = plt.subplots(figsize=(11, 7.5))
        ee = np.linspace(300, 3700, 300)
        a.plot(ee / 1000, h2_curve(ee), "k-", label="H2 curve 0.90 (+) 2.05 / sqrt(E) (+) 1.35 / E (points at 0.73, 2.24, 3.49 GeV)")
        for p, mk in ((gp, "o"), ("legacy", "s")):
            for prot, fill in (("out-of-sample", True), ("in-sample", False)):
                bins = [r for r in rows if r["path"] == p and r["protocol"] == prot and "bin" in r["quantity"]]
                a.errorbar([r["E"] / 1000 for r in bins], [r["res_pooled"] for r in bins], [r["res_pooled_se"] for r in bins], fmt=mk,
                           mfc=None if fill else "none", capsize=3, label=f"Moller per E_pred bin, {p}, {prot} (pooled)")
        a.errorbar([2.237], [1.987], [0.025], fmt="D", color="C3", capsize=3, label="ep, control v3, gate 4 ns, out of sample (1.987 %)")
        hc = next(r for r in rows if r["path"] == gp and r["protocol"] == "out-of-sample" and r["quantity"] == "cluster")
        a.errorbar([hc["E"] / 1000], [hc["res"]], [hc["res_se"]], fmt="*", ms=14, color="C2", capsize=3,
                   label=f"Moller per cluster, all bins, per-split mean ({gp}, out of sample)")
        hs = next(r for r in rows if r["path"] == gp and r["protocol"] == "out-of-sample"
                  and r["quantity"].startswith("sum, symmetric"))
        a.errorbar([hs["E_cluster"] / 1000], [hs["res_cluster_equiv"]], [hs["res_cluster_equiv_se"]], fmt="P", ms=12, color="C4", capsize=3,
                   label="Moller symmetric pairs: sqrt(2) sigma_sum / E_sum (angle independent)")
        a.set(xlabel="E (GeV)", ylabel="sigma / E (%)", ylim=(0, 8),
              title="PbWO4 resolution: Moller (Ta X17) and the H2 curve\n(E / E_pred points include the E_pred smearing from the measured angles)")
        a.legend(fontsize=8)
        a.grid(alpha=0.3)
        fig.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)
        # pages 7-8: example events
        c = C[gp]
        e_oos = keep[gp, "out-of-sample"]["e"]
        npair = c["t"].size // 2
        pick = np.sort(np.random.default_rng(7).choice(npair, min(N_EXAMPLES, npair), replace=False))
        ptr = np.r_[0, np.cumsum(c["cnt"])]
        isel = np.flatnonzero(s["sel"])
        for page in range(0, pick.size, 10):
            fig = plt.figure(figsize=(16, 10.5))
            fig.suptitle(f"Example events {page + 1}-{min(page + 10, pick.size)} (path {gp}; E = out-of-sample Moller gains + position "
                         f"correction; members: module MeV, raw)", fontsize=12, fontweight="bold")
            ax = fig.add_axes([0.01, 0.01, 0.98, 0.93])
            ax.axis("off")
            L = []
            for i in pick[page:page + 10]:
                j = isel[i]
                L.append(f"run {c['run'][2 * i]} event {c['ev'][2 * i]}: |dphi - 180| {s['dphi'][j]:.2f} deg, tan th1 tan th2 (g+1)/2 "
                         f"{s['prod'][j]:.3f}, (E1 + E2) / Ebeam {(e_oos[2 * i] + e_oos[2 * i + 1]) / eb:.4f}")
                for k in (2 * i, 2 * i + 1):
                    mm = np.argsort(-c["me"][ptr[k]:ptr[k + 1]])[:4]
                    mem = " ".join(f"{geo['name'][c['mi'][ptr[k] + q]]}:{c['me'][ptr[k] + q]:.0f}" for q in mm)
                    L.append(f"    seed {geo['name'][int(c['seed'][k])]:>6s} ({c['cnt'][k]:2d} members: {mem:<44s}) x {c['x'][k]:+7.1f} y "
                             f"{c['y'][k]:+7.1f} mm  th {c['th'][k]:.3f} deg  E_raw {c['er'][k]:7.1f}  E {e_oos[k]:7.1f}  E_pred "
                             f"{c['t'][k]:7.1f} MeV  E / E_pred {e_oos[k] / c['t'][k]:.3f}")
            ax.text(0, 1, "\n".join(L), family="monospace", fontsize=7.6, va="top")
            pdf.savefig(fig)
            plt.close(fig)
        # page 9: per split
        fig, axs = plt.subplots(1, 2, figsize=(16, 6))
        fig.suptitle(f"Per-split resolution, path {gp}, Moller gains", fontsize=14, fontweight="bold")
        for a, q in zip(axs, ("cluster", "sum")):
            for prot, mk in (("out-of-sample", "o"), ("in-sample", "s")):
                r = keep[gp, prot]["res"][q]
                a.plot(np.arange(len(r["per"])), r["per"], mk, mfc="none" if prot == "in-sample" else None,
                       label=f"{prot}: {r['res']:.2f} +- {r['res_se']:.2f} %")
            a.set(xlabel="split (run, then split order)", ylabel="CB sigma / E (%)", title=q)
            a.legend()
        fig.tight_layout(rect=(0, 0, 1, 0.93))
        pdf.savefig(fig)
        plt.close(fig)
        # page 10: isolation radius from the data
        fig, axs = plt.subplots(1, 2, figsize=(16, 6.5))
        fig.suptitle(f"Isolation: nearest other cluster (> {CL_STORE:g} MeV) to the pair, path {gp}, every other cut applied, no "
                     f"third-cluster cut", fontsize=13, fontweight="bold")
        dm, esr = iso_study["dmin"], iso_study["esum"]
        a = axs[0]
        a.hist(dm[np.isfinite(dm)], bins=np.arange(0, 400, 5), histtype="step", color="k")
        a.axvline(r_iso, color="r", ls="--", label=f"R_iso {r_iso:g} mm")
        a.set(xlabel="distance to the nearest other cluster (mm)", ylabel="pairs",
              title=f"{int(np.isfinite(dm).sum())} pairs with, {int((~np.isfinite(dm)).sum())} without another cluster")
        a.legend()
        a = axs[1]
        ed = np.array([0, 20, 30, 40, 50, 60, 70, 80, 100, 120, 150, 200, 300])
        xc, pk, md = [], [], []
        for lo, hi in zip(ed[:-1], ed[1:]):
            k = (dm >= lo) & (dm < hi)
            if k.sum() < 100:
                continue
            f = cb_fit(esr[k] * eb, eb)
            xc.append(0.5 * (lo + hi))
            pk.append(100.0 * (f["mu"] / eb - 1.0))
            md.append(100.0 * (float(np.median(esr[k])) - 1.0))
        k = ~np.isfinite(dm)
        a.plot(xc, pk, "o-", label="Crystal Ball peak of E1 + E2")
        a.plot(xc, md, "s--", label="median of E1 + E2")
        if k.sum() >= 100:
            a.axhline(100.0 * (cb_fit(esr[k] * eb, eb)["mu"] / eb - 1.0), color="C0", ls=":", label="peak, no other cluster")
            a.axhline(100.0 * (float(np.median(esr[k])) - 1.0), color="C1", ls=":", label="median, no other cluster")
        a.axvline(r_iso, color="r", ls="--")
        a.set(xlabel="distance to the nearest other cluster (mm)", ylabel="(E1 + E2) / Ebeam - 1 (%), calibration-file energies",
              title="pair-sum shift versus distance")
        a.legend(fontsize=9)
        fig.tight_layout(rect=(0, 0, 1, 0.92))
        pdf.savefig(fig)
        plt.close(fig)
        # page 11: comparison of the third-cluster rules and per-run breakdown
        fig = plt.figure(figsize=(15, 9))
        fig.suptitle("Third-cluster rules and runs (E1 + E2 sigma/E, gate path, out of sample, per-run gains)", fontsize=13, fontweight="bold")
        ax = fig.add_axes([0.02, 0.03, 0.96, 0.88])
        ax.axis("off")
        L = [f"{'selection':<34s}{'pairs':>8s}{'own win':>9s}   {'sum: per split (fitted)':<28s}{'pooled':<16s}{'one cluster 1.12 GeV: split; pooled':<34s}"]
        for x in comp:
            L.append(f"{x['config'] + (' *' if x['main'] else ''):<34s}{x['pairs']:>8d}{x['own']:>9d}   "
                     f"{x['sum']['res']:.2f} +- {x['sum']['res_se']:.2f} ({x['sum']['n_fit_splits']:>2d}){'':<11s}"
                     f"{x['sum']['res_pooled']:.2f} +- {x['sum']['res_pooled_se']:.2f}    "
                     f"{x['sym']['res_cluster_equiv']:.2f} +- {x['sym']['res_cluster_equiv_se']:.2f}; "
                     f"{x['sym']['res_cluster_equiv_pooled']:.2f} +- {x['sym']['res_cluster_equiv_pooled_se']:.2f}")
        L += ["  * = this report", "", f"Per run ({config_name(cl_min, r_iso)}):",
              f"{'run':<7s}{'date':<12s}{'calibration file':<20s}{'pairs':>7s}{'splits':>8s}   {'sum per split':<17s}{'pooled':<16s}"
              f"{'gain L1-2 / calib part':<26s}{'gain L3 / calib part':<22s}"]
        for x in per_run:
            L.append(f"{x['run']:<7d}{x['date']:<12s}{x['calib']:<20s}{x['pairs']:>7d}{x['n_splits']:>8d}   {x['res']:.2f} +- {x['res_se']:.2f}     "
                     f"{x['res_pooled']:.2f} +- {x['res_pooled_se']:.2f}    {x['gain_L12']:.3f} / {x['calib_L12']:.3f}{'':<11s}"
                     f"{x['gain_L3']:.3f} / {x['calib_L3']:.3f}")
        L += ["", "gain = mean over the layer's modules of 1 / (calibration factor x Moller gain), relative to the first run; calib part = the",
              "calibration-file ratio alone.  Layers 1-2 (central block H, 32 modules) are never Moller seeds: their Moller gains are mostly the prior.",
              "Elog: x17_v4 gains 08-08 (inner modules 50-70 % damage); inner-layer HV raised 08-10, partly reverted 08-14."]
        ax.text(0, 1, "\n".join(L), family="monospace", fontsize=8.6, va="top")
        pdf.savefig(fig)
        plt.close(fig)


def main(argv):
    a = list(argv)
    if a[:1] == ["test"]:
        run_tests()
        return
    if a[:1] != ["report"]:
        sys.exit(__doc__)
    a = a[1:]
    kw = {}
    for opt, key, conv in (("--out", "out_root", Path), ("--ref", "ref_root", Path), ("--cache", "cache", Path), ("--r-iso", "r_iso", float),
                           ("--cl-min", "cl_min", lambda v: None if v == "none" else float(v))):
        while opt in a:
            i = a.index(opt)
            kw[key] = conv(a[i + 1])
            del a[i:i + 2]
    if "--combine-runs" in a:                  # variant: one gain solve over all runs (the default solves per run, as the ep chain)
        a.remove("--combine-runs")
        kw["combine_runs"] = True
    if not a or len(a) > 2:
        sys.exit(__doc__)
    report(a[0], float(a[1]) if len(a) > 1 else None, **kw)


if __name__ == "__main__":
    main(erb.take_overrides(sys.argv[1:]))
