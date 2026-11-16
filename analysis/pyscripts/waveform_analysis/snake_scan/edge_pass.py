#!/usr/bin/env python3
"""
edge_pass.py - per-module endpoint edge of the 3x3 energy sum in dwell events (plan P14, D3, work package P14.2).

Used as `snake_scan.py edge --dwell SS/out/dwell_23701.tsv <evio files>`; exposes add_parser/run for the dispatcher.

  resolve_calib(run, db)            calibration file of the run (DB/runinfo/general.json chain) -> Path
  load_calib(path)                  {module name: factor W (MeV/ADC)}
  neighbours(map_path)              {module: [9 member names: the module and its PbWO4 3x3 neighbours]} (sorted by name)
  collect_events(evio, rows, map, calib=None, pipe=None, max_events=0) -> {module: dict(members, I (n,9), n_pileup)}
        Physics events (no LMS/alpha trigger bit) inside the event ranges of the 'dwell' rows of this file;
        I = WaveAnalyzer integral (python/bind_dec.cpp:412, Peak.integral) per member, 0 where the channel did not fire.
        Single event = the central module fired with a peak and no fired member has >= 2 peaks; the rest of the events
        with a central peak are n_pileup (a member that fired with no peak adds 0 and does not flag the event).
  make_hist(ev, W, extra_scale=None) -> dict(sum3, central, n_single, n_pileup, edges, members)
        10 MeV bins 0-4000; MeV = W_j x integral, with W_j divided by extra_scale[j] for neighbours (second pass).
  edge_pass(evio, dwell_rows, map, calib, cache_dir=None, hist_dir=None) -> {module: hist}
  fit_edge(hist) -> (E_edge, sigma, plateau, slope, err, window_spread)
        model (linear plateau) x (erfc step at E_edge, width sigma) + constant, fitted over [0.6, 1.3] x the
        half-plateau crossing; window_spread = (max - min)/mean of E_edge over the windows
        [0.5,1.2], [0.6,1.3], [0.7,1.4] (E_edge is the [0.6,1.3] one); err = fit error of E_edge (MeV).
  second_pass(...), write_edge_tsv(...)
"""
from __future__ import annotations

import argparse
import json
import os
import shutil
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

SS = Path(__file__).resolve().parent
if str(SS) not in sys.path:
    sys.path.insert(0, str(SS))
import snake_scan as SSN  # noqa: E402

E_BEAM = 2108.85
BIN, EMAX = 10.0, 4000.0
EDGES = np.arange(0.0, EMAX + BIN / 2, BIN)
WINDOWS = ((0.5, 1.2), (0.6, 1.3), (0.7, 1.4))
PITCH = 20.77
HIST_DIR = SS / "out" / "edge_hist"
COLS = ["module", "s", "s_err", "window_spread", "n_single", "n_pileup", "s_central", "s_pass1"]


# ---------------------------------------------------------------------------------------------------------------
# calibration, geometry
# ---------------------------------------------------------------------------------------------------------------

def resolve_calib(run: int, db: Path = SSN.DB) -> Path:
    cfg = json.loads((Path(db) / "runinfo" / "general.json").read_text())
    f = None
    for c in sorted(cfg["configurations"], key=lambda c: c.get("from_run", c.get("run_number", 0))):
        if c.get("from_run", c.get("run_number", 0)) <= run and "file" in c.get("calibration", {}):
            f = c["calibration"]["file"]
    if f is None:
        raise RuntimeError(f"no calibration file for run {run}")
    return Path(db) / f


def load_calib(path) -> dict:
    return {e["name"]: float(e["factor"]) for e in json.loads(Path(path).read_text())}


def neighbours(map_path=SSN.DERIVED_MAP) -> dict:
    mods = SSN.load_modules(map_path)
    out = {}
    for n, (x, y) in mods.items():
        mem = [m for m, (xx, yy) in mods.items() if abs(xx - x) < 1.5 * PITCH and abs(yy - y) < 1.5 * PITCH]
        out[n] = sorted(mem)
    return out


def _run_of(path) -> int:
    import re
    m = re.search(r"prad_0*(\d+)\.evio", str(path))
    return int(m.group(1)) if m else -1


# ---------------------------------------------------------------------------------------------------------------
# event pass
# ---------------------------------------------------------------------------------------------------------------

def collect_events(evio, rows, map=SSN.DERIVED_MAP, calib=None, pipe=None, max_events: int = 0) -> dict:
    from prad2py import dec
    fname = Path(evio).name
    ranges = sorted((r["ev_first"], r["ev_last"], r["module"]) for r in rows
                    if r["kind"] == "dwell" and r["file"] == fname)
    if not ranges:
        return {}
    nb = neighbours(map)
    pipe = pipe or SSN.open_pipeline(map, str(evio))
    hy, wa = pipe.hycal, pipe.wave_ana
    ch = SSN.open_evio(evio, pipe)
    cache = {}

    def modname(crate, slot, c):
        k = (crate, slot, c)
        if k not in cache:
            m = hy.module_by_daq(crate, slot, c)
            cache[k] = m.name if (m is not None and m.is_pwo4()) else None
        return cache[k]

    acc = {m: dict(members=nb[m], I=[], n_pileup=0) for _, _, m in ranges}
    ri_cur, nsc = 0, 0
    done = False
    while not done and ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            ch.select_event(i)
            inf = ch.info()
            evn = int(inf.event_number)
            while ri_cur < len(ranges) and evn > ranges[ri_cur][1]:
                ri_cur += 1
            if ri_cur >= len(ranges):
                done = True
                break
            if evn < ranges[ri_cur][0] or int(inf.trigger_bits) & (SSN.TBIT_LMS | SSN.TBIT_ALPHA):
                continue
            nsc += 1
            mod = ranges[ri_cur][2]
            members = acc[mod]["members"]
            pos = {n: k for k, n in enumerate(members)}
            vals = np.zeros(len(members), np.float32)
            pile, central_ok = False, False
            e = ch.fadc()
            for ri in range(e.nrocs):
                roc = e.roc(ri)
                if not roc.present:
                    continue
                crate = pipe.crate_map.get(roc.tag)
                if crate is None:
                    continue
                for s in roc.present_slots():
                    slot = roc.slot(s)
                    for c in slot.present_channels():
                        nm = modname(crate, s, c)
                        k = pos.get(nm)
                        if k is None:
                            continue
                        cd = slot.channel(c)
                        if cd.nsamples < 9:
                            continue
                        _, _, peaks = wa.analyze(np.asarray(cd.samples, np.uint16))
                        if len(peaks) >= 2:
                            pile = True
                        if peaks:
                            vals[k] = peaks[0].integral if len(peaks) == 1 else sum(p.integral for p in peaks)
                            if nm == mod:
                                central_ok = True
            if not central_ok:
                continue
            if pile:
                acc[mod]["n_pileup"] += 1
            else:
                acc[mod]["I"].append(vals)
            if max_events and nsc >= max_events:
                done = True
                break
    for m, a in acc.items():
        a["I"] = np.array(a["I"], np.float32).reshape(-1, len(a["members"]))
    return acc


def make_hist(ev, W, extra_scale=None) -> dict:
    members = list(ev["members"])
    w = np.array([W[n] for n in members], np.float64)
    if extra_scale:
        w = w / np.array([extra_scale.get(n, 1.0) for n in members])
    I = ev["I"].astype(np.float64)
    mev = I * w
    c = members.index(ev["module"]) if "module" in ev else None
    return dict(sum3=np.histogram(mev.sum(1), EDGES)[0], central=np.histogram(mev[:, c], EDGES)[0],
                n_single=len(I), n_pileup=int(ev["n_pileup"]), edges=EDGES, members=np.array(members))


def edge_pass(evio, dwell_rows, map=SSN.DERIVED_MAP, calib=None, cache_dir=None, hist_dir=None, max_events=0) -> dict:
    W = load_calib(calib or resolve_calib(_run_of(evio)))
    acc = collect_events(evio, dwell_rows, map, calib, max_events=max_events)
    out = {}
    for m, a in acc.items():
        a["module"] = m
        if cache_dir:
            np.savez_compressed(Path(cache_dir) / f"{Path(evio).name}_{m}.npz", I=a["I"],
                                members=np.array(a["members"]), n_pileup=a["n_pileup"])
        out[m] = make_hist(a, W)
        out[m]["events"] = a
        if hist_dir:
            save_hist(out[m], Path(hist_dir) / f"edge_{m}_{Path(evio).name}.npz")
    return out


def save_hist(h, path):
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, **{k: v for k, v in h.items() if k != "events"})


def add_hists(a, b):
    return dict(sum3=a["sum3"] + b["sum3"], central=a["central"] + b["central"],
                n_single=a["n_single"] + b["n_single"], n_pileup=a["n_pileup"] + b["n_pileup"],
                edges=a["edges"], members=a["members"])


# ---------------------------------------------------------------------------------------------------------------
# fit
# ---------------------------------------------------------------------------------------------------------------

def _model(E, p0, p1, E0, sig, c, Eref):
    from scipy.special import erfc
    return (p0 + p1 * (E - Eref) / Eref) * 0.5 * erfc((E - E0) / (np.sqrt(2.0) * sig)) + c


def _half_crossing(cnt, ctr):
    """(E_half, plateau): steepest fall of the smoothed spectrum at 0.7-1.4 x 2108.85 MeV -> plateau just below it -> half crossing."""
    sm = np.convolve(cnt.astype(float), np.ones(5) / 5.0, mode="same")
    d = sm[5:] - sm[:-5]
    ok = (ctr[:-5] > 0.7 * E_BEAM) & (ctr[:-5] < 1.4 * E_BEAM)      # edge seed within s = 0.7-1.4
    k = int(np.argmin(np.where(ok, d, np.inf)))
    Eg = ctr[k] + 25.0
    sel = (ctr > 0.7 * Eg) & (ctr < 0.9 * Eg)
    plateau = float(np.mean(sm[sel]))
    j = int(np.argmax((ctr > 0.9 * Eg) & (sm < 0.5 * plateau)))
    j0 = j - 1
    t = (0.5 * plateau - sm[j0]) / (sm[j] - sm[j0]) if sm[j] != sm[j0] else 0.0
    return float(ctr[j0] + t * BIN), plateau


def _fit_window(cnt, ctr, Eh, plateau, lo, hi):
    from scipy.optimize import curve_fit
    sel = (ctr >= lo * Eh) & (ctr <= hi * Eh)
    x, y = ctr[sel], cnt[sel].astype(float)
    sg = np.sqrt(np.maximum(y, 1.0))
    f = lambda E, p0, p1, E0, sig, c: _model(E, p0, p1, E0, sig, c, Eh)
    lb = [0.0, -np.inf, 0.7 * Eh, 10.0, 0.0]
    ub = [np.inf, np.inf, 1.3 * Eh, 0.15 * Eh, np.inf]
    p, cov = curve_fit(f, x, y, p0=[plateau, 0.0, Eh, 0.03 * Eh, 0.02 * plateau], sigma=sg, absolute_sigma=True,
                       bounds=(lb, ub), maxfev=20000)
    return p, np.sqrt(np.diag(cov))


def fit_edge(hist, which: str = "sum3"):
    cnt = np.asarray(hist[which])
    edges = np.asarray(hist["edges"])
    ctr = 0.5 * (edges[1:] + edges[:-1])
    Eh, plateau = _half_crossing(cnt, ctr)
    res = [_fit_window(cnt, ctr, Eh, plateau, lo, hi) for lo, hi in WINDOWS]
    E = np.array([r[0][2] for r in res])
    p, e = res[1]
    spread = float((E.max() - E.min()) / E.mean())
    return float(p[2]), float(abs(p[3])), float(p[0]), float(p[1] / Eh), float(e[2]), spread


# ---------------------------------------------------------------------------------------------------------------
# second pass, table
# ---------------------------------------------------------------------------------------------------------------

def _fit_all(hists):
    return {m: fit_edge(h) for m, h in hists.items()}


def second_pass(events: dict, W: dict, s1: dict) -> dict:
    """Re-histogram with W_j / s_j for the neighbours j of each module that have a dwell (the module itself keeps W)."""
    out = {}
    for m, ev in events.items():
        scale = {n: s1[n] for n in ev["members"] if n != m and n in s1}
        out[m] = make_hist(ev, W, scale)
    return out


def write_edge_tsv(path, rows, sources, map_path, calib, note=""):
    with open(path, "w") as f:
        f.write(f"# source: {' '.join(str(s) for s in sources)} ; map {map_path} ; calib {calib}\n")
        f.write("# s = E_edge/2108.85 of the 3x3 MeV sum (second pass: neighbours with a dwell use W_j/s_j); s_err = hypot(fit "
                "error, half window spread x s); window_spread = (max-min)/mean over 3 windows; s_central = same for the "
                "central module alone; s_pass1 = s before the neighbour pass" + note + "\n")
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(str(r[c]) if c in ("module", "n_single", "n_pileup") else f"{r[c]:.5f}" for c in COLS) + "\n")


def _worker(a):
    evio, rows, map_path, calib, cache_dir, hist_dir, mx = a
    t0 = time.time()
    h = edge_pass(evio, rows, map_path, calib, cache_dir=cache_dir, hist_dir=None, max_events=mx)
    for v in h.values():
        v.pop("events", None)
    return evio, h, time.time() - t0


def run(args):
    rows = SSN.read_dwell_tsv(args.dwell)
    files = [Path(f) for f in args.evio]
    calib = args.calib or resolve_calib(_run_of(files[0]))
    W = load_calib(calib)
    cache = Path(args.cache_dir) if args.cache_dir else Path(tempfile.mkdtemp(prefix="sup-edge-"))
    cache.mkdir(parents=True, exist_ok=True)
    hist_dir = Path(args.hist_dir)
    try:
        jobs = [(str(f), rows, args.map, str(calib), str(cache), None, args.max_events) for f in files
                if any(r["file"] == f.name and r["kind"] == "dwell" for r in rows)]
        if args.from_cache:
            jobs = []
        if args.jobs > 1:
            from multiprocessing import Pool
            with Pool(min(args.jobs, 4)) as p:
                res = p.map(_worker, jobs)
        else:
            res = [_worker(j) for j in jobs]
        ev = {}
        for evio, h, dt in res:
            print(f"[edge] {Path(evio).name}: {', '.join(f'{m} n1={v['n_single']} np={v['n_pileup']}' for m, v in h.items())} "
                  f"({dt:.0f} s)")
        for f in cache.glob("*.npz"):                   # combine the event caches per module
            d = np.load(f)
            m = f.name.split("_")[-1][:-4]
            a = ev.setdefault(m, dict(module=m, members=list(d["members"]), I=[], n_pileup=0))
            a["I"].append(d["I"])
            a["n_pileup"] += int(d["n_pileup"])
        for a in ev.values():
            a["I"] = np.concatenate(a["I"])
        hists = {m: make_hist(a, W) for m, a in ev.items()}      # modules with dwell rows in several files combine here
        fit1 = _fit_all(hists)
        s1 = {m: v[0] / E_BEAM for m, v in fit1.items()}
        hists2 = second_pass(ev, W, s1)
        fit2 = _fit_all(hists2)
        out = []
        for m in sorted(hists2, key=lambda n: int(n[1:])):
            E, sig, pl, sl, err, spread = fit2[m]
            sc = fit_edge(hists2[m], "central")[0] / E_BEAM
            s = E / E_BEAM
            h = hists2[m]
            save_hist(h, hist_dir / f"edge_{m}.npz")
            out.append(dict(module=m, s=s, s_err=float(np.hypot(err / E_BEAM, 0.5 * spread * s)), window_spread=spread,
                            n_single=h["n_single"], n_pileup=h["n_pileup"], s_central=sc, s_pass1=s1[m]))
            print(f"[edge] {m}: E_edge={E:.1f} MeV sigma={sig:.1f} s={s:.4f} (pass1 {s1[m]:.4f}, "
                  f"{100 * (s / s1[m] - 1):+.2f}%) spread={spread:.4f} s_central={sc:.4f}")
        write_edge_tsv(args.out, out, [str(f) for f in files], args.map, calib)
        print(f"[edge] wrote {args.out}")
    finally:
        if not args.cache_dir:
            shutil.rmtree(cache, ignore_errors=True)
    return 0


def add_parser(sub):
    p = sub.add_parser("edge", help="endpoint edge of the 3x3 sum per dwell module (D3)")
    p.add_argument("evio", nargs="+")
    p.add_argument("--dwell", required=True, help="dwell TSV (snake_scan.py path --out)")
    p.add_argument("--map", default=str(SSN.DERIVED_MAP))
    p.add_argument("--calib", default="", help="calibration JSON (default: the run's file of runinfo/general.json)")
    p.add_argument("--out", default=str(SS / "out" / "edge_23701.tsv"))
    p.add_argument("--hist-dir", default=str(HIST_DIR))
    p.add_argument("--cache-dir", default="", help="per-file/module event cache (npz); temporary when omitted")
    p.add_argument("--from-cache", action="store_true", help="skip the event pass; use the npz files of --cache-dir")
    p.add_argument("--max-events", type=int, default=0, help="stop after N dwell events per file (testing)")
    p.add_argument("-j", "--jobs", type=int, default=1, help="parallel files (max 4)")
    p.set_defaults(fn=run)
    return p


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    add_parser(ap.add_subparsers(dest="cmd", required=True))
    sys.exit(run(ap.parse_args()) or 0)
