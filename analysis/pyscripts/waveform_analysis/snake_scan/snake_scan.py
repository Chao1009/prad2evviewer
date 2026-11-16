#!/usr/bin/env python3
"""
snake_scan.py - snake-scan constants through LMS ratios (plan P14, notes/plan_p14_snake_scan_20261114.md).

Subcommands (D10):  map | path | edge | lms | transfer.
  map   : derive the HyCal map for runs <= 23849 (D1) -> SS/hycal_map_til_run023849.json
  path  : per-file dwell table from the occupancy argmax + EPICS ptrans (D2) and the ptrans->module affine map
  edge, lms, transfer : live in separate optional modules edge_pass.py, lms_pass.py, transfer.py (written by
          later work packages).  Convention: each module exposes
              add_parser(subparsers)   register its subcommand on the shared argparse subparsers
              run(args)                execute it (returns an int exit code or None)
          If the module file is absent the dispatcher prints "not implemented".

Environment: cd analysis/pyscripts && source setup_env.sh ; MPLBACKEND=Agg.
SS = this directory.  Reusable helpers (import with sys.path pointing at SS, `import snake_scan as SSN`):

  derive_map(hycal_map: Path, history: Path, out: Path) -> list[dict]
      hycal_map.json entries ({n,t,geo,daq}) with `daq` replaced by the history entry of the same name; history
      names missing from the map are appended (LMSP, t "LMS", geo continuing the LMS row); writes `out`.
  map_diff(hycal_map, derived) -> list[(name, old_daq|None, new_daq)]       changed or added entries
  open_pipeline(map_path: str = DERIVED_MAP, evio_path: str = "") -> _common.Pipeline
      _common.setup_pipeline(evio_path=..., hc_map_file=map_path): the C++ PipelineBuilder with
      set_hycal_map(map_path).  Run number sniffed from evio_path (runinfo / calibration of that run).
  open_evio(path, pipe) -> dec.EvChannel          opened and configured; loop `while ch.read()==dec.Status.success`
  parse_epics_text(text: str) -> dict[str, float]  EPICS payload "value  name" lines; a repeated name keeps the last
  scan_blocks(evio, map_path=DERIVED_MAP, block=5000, pipe=None) -> (blocks, epics)
      blocks: list of dict(ev_first, ev_last, ts_first, ts_last, n, counts{module: n_argmax}); epics: list of
      dict(idx, ts, vals) with idx = number of scanned physics events before the snapshot.
  dwell_table(evio, map, block=5000, ...) -> rows(file, ev_first, ev_last, module, share, n_events, ptrans_x,
      ptrans_y, ts_first, ts_last, kind)      kind = dwell | ramp (D2)
  fit_affine(rows, hycal_modules) -> dict     ptrans -> module position, per axis (x = a*ptrans_x + b), residuals
  ptrans_to_module(ptx, pty, affine, modules) -> (name, distance_mm)       nearest PbWO4 module
  write_dwell_tsv(rows, path, sources, map_path)   TSV with a `# source:` line
  read_dwell_tsv(path) -> rows
Event selection for the path scan: Physics events without LMS (bit 24) or alpha (bit 25) trigger bits; the argmax is
over PbWO4 modules of max(x - median(x[:8])).
"""
from __future__ import annotations

import argparse
import importlib
import json
import os
import sys
import time
from pathlib import Path

import numpy as np

SS = Path(__file__).resolve().parent
PYSCRIPTS = SS.parents[1]
DB = Path(os.environ.get("PRAD2_DATABASE_DIR", str(Path.home() / "work/PRad/prad2evviewer/database")))
DERIVED_MAP = SS / "hycal_map_til_run023849.json"
HISTORY = DB / "history" / "daq_map_20260424_til_run023849.json"
TBIT_LMS, TBIT_ALPHA = 1 << 24, 1 << 25
PEAK_HI = 300.0           # ADC; an event "has beam" when its largest PbWO4 peak sample exceeds this (W564 is weak)
SHARE_MIN = 0.90          # D2: dwell block purity
WEAK_FRAC = 0.5           # blocks with fewer than this fraction of beam events (peak > PEAK_HI) are "weak beam" and
SHARE_MIN_WEAK = 0.75     # use this purity (background modules take a larger fraction of the few hi events, W564)
MIN_DWELL_BLOCKS = 10     # a dwell lasts >= 50 k events (modules take 1.5-2 M)
MAX_PTRANS_SPAN = 2.0     # mm; EPICS position moving inside a candidate dwell -> ramp
if str(PYSCRIPTS) not in sys.path:
    sys.path.insert(0, str(PYSCRIPTS))


# ---------------------------------------------------------------------------------------------------------------
# D1: derived map
# ---------------------------------------------------------------------------------------------------------------

def derive_map(hycal_map: Path, history: Path, out: Path) -> list:
    entries = json.loads(Path(hycal_map).read_text())
    hist = {e["name"]: e for e in json.loads(Path(history).read_text())}
    names = set()
    for e in entries:
        names.add(e["n"])
        h = hist.get(e["n"])
        if h is not None:
            e["daq"] = {"crate": h["crate"], "slot": h["slot"], "channel": h["channel"]}
    lms = [e for e in entries if e["t"] == "LMS"]
    for nm, h in hist.items():
        if nm in names:
            continue
        if not nm.startswith("LMS"):
            raise ValueError(f"history entry {nm} is not in the map and not an LMS channel")
        ref = max(lms, key=lambda e: e["geo"]["x"])
        geo = dict(ref["geo"])
        geo["x"] = geo["x"] + 50.0          # LMSP: one more 50 mm step along the LMS row (no geometry is used)
        entries.append({"n": nm, "t": "LMS", "geo": geo,
                        "daq": {"crate": h["crate"], "slot": h["slot"], "channel": h["channel"]}})
    Path(out).write_text(json.dumps(entries, indent=1) + "\n")
    return entries


def map_diff(hycal_map: Path, derived) -> list:
    old = {e["n"]: e.get("daq") for e in json.loads(Path(hycal_map).read_text())}
    if not isinstance(derived, list):
        derived = json.loads(Path(derived).read_text())
    return [(e["n"], old.get(e["n"]), e["daq"]) for e in derived if old.get(e["n"]) != e["daq"]]


# ---------------------------------------------------------------------------------------------------------------
# pipeline / EVIO / EPICS helpers
# ---------------------------------------------------------------------------------------------------------------

def open_pipeline(map_path=DERIVED_MAP, evio_path: str = ""):
    import _common as C
    return C.setup_pipeline(evio_path=evio_path, hc_map_file=str(map_path))


def open_evio(path, pipe):
    from prad2py import dec
    ch = dec.EvChannel()
    ch.set_config(pipe.cfg)
    if ch.open_auto(str(path)) != dec.Status.success:
        raise RuntimeError(f"cannot open {path}")
    return ch


def parse_epics_text(text: str) -> dict:
    out = {}
    for ln in text.splitlines():
        parts = ln.split()
        if len(parts) < 2:
            continue
        try:
            out[parts[1]] = float(parts[0])
        except ValueError:
            continue
    return out


# ---------------------------------------------------------------------------------------------------------------
# D2: blocks, dwells
# ---------------------------------------------------------------------------------------------------------------

def scan_blocks(evio, map_path=DERIVED_MAP, block: int = 5000, pipe=None, max_events: int = 0, cache=""):
    """Per-block occupancy argmax counts.  `cache` (a JSON path) is read if present and written otherwise."""
    if cache and Path(cache).exists():
        d = json.loads(Path(cache).read_text())
        return d["blocks"], d["epics"]
    res = _scan_blocks(evio, map_path, block, pipe, max_events)
    if cache:
        Path(cache).write_text(json.dumps(dict(blocks=res[0], epics=res[1])))
    return res


def _scan_blocks(evio, map_path, block, pipe, max_events):
    from prad2py import dec
    pipe = pipe or open_pipeline(map_path, str(evio))
    hy = pipe.hycal
    ch = open_evio(evio, pipe)
    cache = {}                                  # (crate, slot, ch) -> module name or None (not PbWO4)

    def modname(crate, slot, c):
        k = (crate, slot, c)
        if k not in cache:
            m = hy.module_by_daq(crate, slot, c)
            cache[k] = m.name if (m is not None and m.is_pwo4()) else None
        return cache[k]

    blocks, epics = [], []
    cur = None
    nsc = 0

    def flush():
        nonlocal cur
        if cur is not None and cur["n"] > 0:
            blocks.append(cur)
        cur = None

    while ch.read() == dec.Status.success:
        if not ch.scan():
            continue
        et = ch.get_event_type()
        if et == dec.EventType.Epics:
            txt = ch.extract_epics_text()
            if txt:
                ts = blocks[-1]["ts_last"] if cur is None and blocks else (cur["ts_last"] if cur else 0)
                epics.append(dict(idx=nsc, ts=ts, vals=parse_epics_text(txt)))
            continue
        if et != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            ch.select_event(i)
            inf = ch.info()
            tb = int(inf.trigger_bits)
            if tb & (TBIT_LMS | TBIT_ALPHA):
                continue
            if cur is None:
                cur = dict(ev_first=int(inf.event_number), ts_first=int(inf.timestamp), n=0, n_hi=0,
                           counts={}, counts_hi={})
            cur["ev_last"] = int(inf.event_number)
            cur["ts_last"] = int(inf.timestamp)
            cur["n"] += 1
            nsc += 1
            e = ch.fadc()
            best, best_v = None, -1.0
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
                        if nm is None:
                            continue
                        cd = slot.channel(c)
                        if cd.nsamples < 9:
                            continue
                        x = np.asarray(cd.samples, np.float64)
                        v = x.max() - np.median(x[:8])
                        if v > best_v:
                            best, best_v = nm, v
            if best is not None:
                cur["counts"][best] = cur["counts"].get(best, 0) + 1
                if best_v >= PEAK_HI:
                    cur["n_hi"] += 1
                    cur["counts_hi"][best] = cur["counts_hi"].get(best, 0) + 1
            if cur["n"] >= block:
                flush()
            if max_events and nsc >= max_events:
                break
        if max_events and nsc >= max_events:
            break
    flush()
    # a short trailing block (< 1000 events) is folded into its predecessor
    if len(blocks) > 1 and blocks[-1]["n"] < 1000:
        a, b = blocks[-2], blocks.pop()
        a["ev_last"], a["ts_last"], a["n"] = b["ev_last"], b["ts_last"], a["n"] + b["n"]
        a["n_hi"] += b["n_hi"]
        for key in ("counts", "counts_hi"):
            for k, v in b[key].items():
                a[key][k] = a[key].get(k, 0) + v
    return blocks, epics


def _top(counts, n):
    if not counts:
        return "", 0.0
    k = max(counts, key=counts.get)
    return k, counts[k] / n if n else 0.0


def _ptrans(epics, i0, i1, last_before):
    """Median ptrans_x/y of the snapshots with scan index in [i0, i1) (else the last one before i0) and the
    max-min span of ptrans_x, ptrans_y."""
    sel = [e for e in epics if i0 <= e["idx"] < i1 and "ptrans_x" in e["vals"]]
    if not sel:
        sel = [last_before] if last_before else []
    if not sel:
        return float("nan"), float("nan"), 0.0
    px = np.array([e["vals"]["ptrans_x"] for e in sel])
    py = np.array([e["vals"]["ptrans_y"] for e in sel])
    return float(np.median(px)), float(np.median(py)), float(max(np.ptp(px), np.ptp(py)))


def blocks_to_rows(fname, blocks, epics, share_min=SHARE_MIN, min_blocks=MIN_DWELL_BLOCKS):
    # scan index at the start of each block
    starts = np.cumsum([0] + [b["n"] for b in blocks])
    lab = []
    for b in blocks:
        m, sh = _top(b["counts_hi"], b["n_hi"])
        if b["n_hi"] < WEAK_FRAC * b["n"]:
            sh = sh * share_min / SHARE_MIN_WEAK       # rescale so that the common test sh >= share_min applies
        lab.append((m, sh))
    # maximal runs of clean blocks with one argmax
    cand = []
    i = 0
    while i < len(blocks):
        m, sh = lab[i]
        if sh >= share_min:
            j = i
            while j + 1 < len(blocks) and lab[j + 1][0] == m and lab[j + 1][1] >= share_min:
                j += 1
            cand.append((i, j + 1, m))
            i = j + 1
        else:
            i += 1
    runs = []
    for a, b, m in cand:
        # ptrans is the transporter set point: the first snapshot differing from the run's first by more than
        # MAX_PTRANS_SPAN marks the start of the next move, blocks from there on are ramp
        snaps = [e for e in epics if starts[a] <= e["idx"] < starts[b] and "ptrans_x" in e["vals"]]
        if snaps:
            p0 = (snaps[0]["vals"]["ptrans_x"], snaps[0]["vals"]["ptrans_y"])
            for e in snaps[1:]:
                if max(abs(e["vals"]["ptrans_x"] - p0[0]), abs(e["vals"]["ptrans_y"] - p0[1])) > MAX_PTRANS_SPAN:
                    b = int(np.searchsorted(starts, e["idx"], side="right")) - 1
                    break
        if b - a >= min_blocks:
            runs.append((a, b, "dwell"))
    # group everything else into ramp runs
    allruns, k = [], 0
    dw = {a: (a, b) for a, b, _ in runs}
    while k < len(blocks):
        if k in dw:
            allruns.append((dw[k][0], dw[k][1], "dwell"))
            k = dw[k][1]
        else:
            j = k
            while j + 1 < len(blocks) and (j + 1) not in dw:
                j += 1
            allruns.append((k, j + 1, "ramp"))
            k = j + 1
    rows = []
    for a, b, kd in allruns:
        counts, n, nhi = {}, 0, 0
        for blk in blocks[a:b]:
            n += blk["n"]
            nhi += blk["n_hi"]
            for mm, c in blk["counts_hi"].items():
                counts[mm] = counts.get(mm, 0) + c
        mod, share = _top(counts, nhi)
        # a ramp row is named by its partner module: the one the beam leaves (ramp before a dwell, first block)
        # or reaches (ramp after a dwell, last block); share is that module's fraction over the row
        if kd == "ramp":
            if allruns and a > 0 and any(r[2] == "dwell" and r[1] == a for r in allruns):
                mod = _top(blocks[b - 1]["counts_hi"], blocks[b - 1]["n_hi"])[0] or mod
            elif any(r[2] == "dwell" and r[0] == b for r in allruns):
                mod = _top(blocks[a]["counts_hi"], blocks[a]["n_hi"])[0] or mod
            share = counts.get(mod, 0) / nhi if nhi else 0.0
        before = [e for e in epics if e["idx"] < starts[a] and "ptrans_x" in e["vals"]]
        px, py, _ = _ptrans(epics, starts[a], starts[b], before[-1] if before else None)
        rows.append(dict(file=os.path.basename(str(fname)), ev_first=blocks[a]["ev_first"],
                         ev_last=blocks[b - 1]["ev_last"], module=mod, share=share, n_events=n,
                         ptrans_x=px, ptrans_y=py, ts_first=blocks[a]["ts_first"],
                         ts_last=blocks[b - 1]["ts_last"], kind=kd))
    return rows


def dwell_table(evio, map=DERIVED_MAP, block: int = 5000, pipe=None, max_events: int = 0, cache=""):
    blocks, epics = scan_blocks(evio, map, block, pipe, max_events, cache)
    return blocks_to_rows(evio, blocks, epics)


COLS = ["file", "ev_first", "ev_last", "module", "share", "n_events", "ptrans_x", "ptrans_y",
        "ts_first", "ts_last", "kind"]


def write_dwell_tsv(rows, path, sources, map_path=DERIVED_MAP):
    with open(path, "w") as f:
        f.write(f"# source: {' '.join(str(s) for s in sources)} ; map {map_path}\n")
        f.write(f"# share = fraction of events with a PbWO4 peak sample > {PEAK_HI:.0f} ADC whose argmax is the module; "
                "n_events counts all physics events (no LMS/alpha) of the range; ts in TI ticks (4 ns)\n")
        f.write("\t".join(COLS) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[c]:.4f}" if c in ("share", "ptrans_x", "ptrans_y") else str(r[c])
                              for c in COLS) + "\n")


def read_dwell_tsv(path):
    rows = []
    for ln in Path(path).read_text().splitlines():
        if ln.startswith("#") or ln.startswith("file\t") or not ln.strip():
            continue
        v = ln.split("\t")
        r = dict(zip(COLS, v))
        for c in ("ev_first", "ev_last", "n_events", "ts_first", "ts_last"):
            r[c] = int(r[c])
        for c in ("share", "ptrans_x", "ptrans_y"):
            r[c] = float(r[c])
        rows.append(r)
    return rows


# ---------------------------------------------------------------------------------------------------------------
# EPICS ptrans -> module
# ---------------------------------------------------------------------------------------------------------------

def load_modules(map_path=DERIVED_MAP):
    """name -> (x, y) of the PbWO4 modules (map geometry, mm)."""
    return {e["n"]: (e["geo"]["x"], e["geo"]["y"]) for e in json.loads(Path(map_path).read_text())
            if e["t"] == "PbWO4"}


def fit_affine(rows, modules):
    """Per-axis linear map ptrans -> module centre fitted on the dwell rows (module position = a*ptrans + b).
    An axis whose ptrans does not vary (rms < 1 mm) is fitted with unit slope of the sign that fits better."""
    d = [r for r in rows if r["kind"] == "dwell" and r["module"] in modules and np.isfinite(r["ptrans_x"])]
    out = {"n": len(d)}
    for ax, key in (("x", "ptrans_x"), ("y", "ptrans_y")):
        p = np.array([r[key] for r in d])
        m = np.array([modules[r["module"]][0 if ax == "x" else 1] for r in d])
        if np.std(p) >= 1.0 and len(set(m.round(1))) >= 3:
            a, b = np.polyfit(p, m, 1)
            fixed = False
        else:
            best = None
            for a_ in (1.0, -1.0):
                b_ = float(np.mean(m - a_ * p))
                rs = float(np.std(m - a_ * p - b_))
                if best is None or rs < best[0]:
                    best = (rs, a_, b_)
            _, a, b = best
            fixed = True
        res = m - (a * p + b)
        out[ax] = dict(a=float(a), b=float(b), rms=float(np.sqrt(np.mean(res ** 2))),
                       maxabs=float(np.abs(res).max()), slope_fixed=fixed, ptrans_rms=float(np.std(p)),
                       residuals={r["module"]: float(x) for r, x in zip(d, res)})
    return out


def ptrans_to_module(ptx, pty, affine, modules):
    x = affine["x"]["a"] * ptx + affine["x"]["b"]
    y = affine["y"]["a"] * pty + affine["y"]["b"]
    best = min(modules, key=lambda k: (modules[k][0] - x) ** 2 + (modules[k][1] - y) ** 2)
    return best, float(np.hypot(modules[best][0] - x, modules[best][1] - y))


def write_affine(aff, rows, path, sources, map_path=DERIVED_MAP):
    modules = load_modules(map_path)
    L = [f"# source: {' '.join(str(s) for s in sources)} ; map {map_path}",
         f"# affine map of EPICS ptrans (mm) to HyCal module-centre position (mm), fitted on {aff['n']} clean dwells",
         "# module position = a * ptrans + b, per axis (axis-aligned; ptrans_y spans only rows 17-18)"]
    for ax in ("x", "y"):
        q = aff[ax]
        L.append(f"{ax}: a={q['a']:.4f} b={q['b']:.3f} mm  residual rms={q['rms']:.2f} mm  max|res|={q['maxabs']:.2f} mm"
                 f"  ptrans_rms={q['ptrans_rms']:.1f} mm  slope_fixed={q['slope_fixed']}")
    L.append("module\tptrans_x\tptrans_y\tres_x_mm\tres_y_mm\tnearest_module\tdist_mm")
    for r in rows:
        if r["kind"] != "dwell":
            continue
        nm, dist = ptrans_to_module(r["ptrans_x"], r["ptrans_y"], aff, modules)
        L.append(f"{r['module']}\t{r['ptrans_x']:.2f}\t{r['ptrans_y']:.2f}\t{aff['x']['residuals'].get(r['module'], float('nan')):.2f}"
                 f"\t{aff['y']['residuals'].get(r['module'], float('nan')):.2f}\t{nm}\t{dist:.2f}")
    verdict = "usable" if max(aff['x']['maxabs'], aff['y']['maxabs']) < 5.0 else "NOT usable (residual > 5 mm): use occupancy argmax"
    L.append(f"# verdict: ptrans {verdict}")
    Path(path).write_text("\n".join(L) + "\n")


# ---------------------------------------------------------------------------------------------------------------
# subcommands
# ---------------------------------------------------------------------------------------------------------------

def _path_worker(a):
    f, mp, block, maxev, cache = a
    t0 = time.time()
    rows = dwell_table(f, mp, block, max_events=maxev, cache=cache)
    print(f"[path] {os.path.basename(f)}: {len(rows)} rows ({time.time() - t0:.0f} s)", flush=True)
    return rows


def cmd_map(a):
    ents = derive_map(Path(a.hycal_map), Path(a.history), Path(a.out))
    d = map_diff(Path(a.hycal_map), ents)
    print(f"[map] wrote {a.out}: {len(ents)} entries, {len(d)} changed/added daq entries")
    for nm, o, n in d:
        print(f"  {nm}: {o} -> {n}")


def cmd_path(a):
    files = list(a.evio)
    cdir = Path(a.cache_dir) if a.cache_dir else None
    if cdir:
        cdir.mkdir(parents=True, exist_ok=True)
    jobs = [(f, a.map, a.block, a.max_events, str(cdir / (os.path.basename(f) + ".blocks.json")) if cdir else "")
            for f in files]
    if a.jobs > 1 and len(jobs) > 1:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(min(a.jobs, 4)) as pool:
            res = pool.map(_path_worker, jobs, chunksize=1)
    else:
        res = [_path_worker(j) for j in jobs]
    rows = [r for rr in res for r in rr]
    hdr = "\t".join(COLS)
    print(hdr)
    for r in rows:
        print("\t".join(f"{r[c]:.4f}" if c in ("share", "ptrans_x", "ptrans_y") else str(r[c]) for c in COLS))
    if a.out:
        write_dwell_tsv(rows, a.out, files, a.map)
    if a.affine_out:
        aff = fit_affine(rows, load_modules(a.map))
        write_affine(aff, rows, a.affine_out, files, a.map)
        print(f"[path] affine: x a={aff['x']['a']:.3f} b={aff['x']['b']:.1f} rms={aff['x']['rms']:.2f} mm; "
              f"y a={aff['y']['a']:.3f} b={aff['y']['b']:.1f} rms={aff['y']['rms']:.2f} mm")


OPTIONAL = {"edge": "edge_pass", "lms": "lms_pass", "transfer": "transfer"}


def _load_optional(name):
    """Import SS/<module>.py if it exists; None when absent (a later work package writes it)."""
    if not (SS / f"{OPTIONAL[name]}.py").exists():
        return None
    if str(SS) not in sys.path:
        sys.path.insert(0, str(SS))
    return importlib.import_module(OPTIONAL[name])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("map", help="derive the HyCal map for runs <= 23849 (D1)")
    p.add_argument("--hycal-map", default=str(DB / "hycal_map.json"))
    p.add_argument("--history", default=str(HISTORY))
    p.add_argument("--out", default=str(DERIVED_MAP))
    p.set_defaults(fn=cmd_map)
    p = sub.add_parser("path", help="dwell table of EVIO file(s) (D2)")
    p.add_argument("evio", nargs="+")
    p.add_argument("--map", default=str(DERIVED_MAP))
    p.add_argument("--block", type=int, default=5000)
    p.add_argument("--max-events", type=int, default=0, help="stop after N scanned events per file (testing)")
    p.add_argument("-j", "--jobs", type=int, default=1, help="parallel files (max 4)")
    p.add_argument("--cache-dir", default="", help="per-file block cache (JSON), reused when present")
    p.add_argument("--out", default="", help="dwell TSV")
    p.add_argument("--affine-out", default="", help="ptrans -> module affine fit text")
    p.set_defaults(fn=cmd_path)
    for name in OPTIONAL:
        mod = _load_optional(name)
        if mod is not None:
            mod.add_parser(sub)
        else:
            q = sub.add_parser(name, help=f"{name} (not implemented: {OPTIONAL[name]}.py absent)",
                               add_help=False)
            q.add_argument("rest", nargs=argparse.REMAINDER)
            q.set_defaults(fn=lambda a, n=name: print(f"{n}: not implemented ({OPTIONAL[n]}.py absent)") or 0)
    a = ap.parse_args(argv)
    if a.cmd in OPTIONAL:
        mod = _load_optional(a.cmd)
        return mod.run(a) if mod is not None else a.fn(a)
    return a.fn(a)


if __name__ == "__main__":
    sys.exit(main() or 0)
