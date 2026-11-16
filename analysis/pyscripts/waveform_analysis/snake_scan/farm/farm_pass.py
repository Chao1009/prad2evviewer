#!/usr/bin/env python3
"""
farm_pass.py - P14.7 per-file product of one EVIO split (single heavy decode loop), run by farm_job.sh.

    python farm_pass.py EVIO OUTDIR [--map PATH] [--calib PATH] [--max-events N] [--mode auto|scan|cosmic] [--force]

writes OUTDIR/<run>/<evio basename>.npz (atomic: written as a partial file, then renamed).

Mode scan (every run except 23744): ONE Python loop over the physics events (the expensive part, about 0.5 ms/event of
channel iteration; the C++ decode itself is about 2 s per file) produces, per file,
  * the dwell classifier input exactly as snake_scan._scan_blocks (5000-event blocks without LMS/alpha trigger bits,
    argmax of max(x - median(x[:8])) over PbWO4 channels, counts of events with the largest peak sample >= PEAK_HI,
    short trailing block folded) and the EPICS ptrans snapshots -> farm_merge.py applies SSN.blocks_to_rows unchanged;
  * the edge-pass inputs per event, around a *running* beam module: the 3x3 members of the module that is the mode of
    the argmax of the last GUESS_N events with a large peak (stable only: mode share >= GUESS_MIN_SHARE).  Inside a dwell
    this equals the dwell module of wave 1 (97 % of the events have it as argmax); only the first GUESS_WARM large-peak
    events of a file or dwell are lost.  Per event: WaveAnalyzer integral (sum of all peaks) and number of peaks of each
    member (edge_pass.collect_events convention).  Counters per (block, module): events with a central peak, wave-1 single
    events (no member with >= 2 peaks), central-single events (central has exactly 1 peak) -> n_single/n_pileup of every
    pile-up setting without the events below the energy cut.  The per-event arrays are kept only for blocks whose
    large-peak argmax share is >= KEEP_SHARE and for events of that block's module with a 3x3 MeV sum inside
    [KEEP_MEV_LO, KEEP_MEV_HI] (the edge fit never uses anything below 0.5 x half-plateau crossing >= 700 MeV) and with
    event_number % KEEP_EVERY == 0 (size: about 10 MB per dwell file at 2, 20 MB at 1; the fit needs only about 1e5 events);
  * the LMS/alpha/pedestal series by calling lms_pass.lms_pass unchanged (a second decode of about 2-5 s: no Python
    channel loop runs for the other events) -> keys with the prefix L_.
Mode cosmic (run 23744): cosmic_mip.scan_file unchanged (MIP hits of vertical tracks, LMS and alpha events).

npz keys, scan mode (all arrays; strings as 0-d unicode):
  meta_json                    file, run, map, calib, n_physics, n_scanned, seconds, version, constants
  mod_names[nM]                PbWO4 module names (ids of ev_ctr / cnt_ctr index this)
  b_evfirst, b_evlast, b_tsfirst, b_tslast [nb] int64; b_n, b_nhi [nb] int32; b_counts_hi_json   blocks (folded)
  epics_json                   [{idx, ts, vals:{ptrans_*}}]
  cnt[k, 5] int32              block index, module id, n_central_peak, n_single_w1, n_single_central
  ev_evn[m] uint32, ev_ctr[m] uint16, ev_I[m, 9] float32, ev_npk[m, 9] uint8      kept per-event arrays; columns = members of
                               edge_pass.neighbours(map)[module] (sorted by name), zero-padded to 9 (modules at the HyCal edge/hole)
  L_*                          lms_pass keys (lms_I, lms_ref, alpha_ref, lms_ts, alpha_ts, lms_nch, lms_pk_med, ped, mod_names, n_all)
npz keys, cosmic mode: see cosmic_mip.scan_file; hits_* (module, peak_adc, integral_adc, column_len), counts_json,
  trig_hi_json, lms_*, alpha_*, names, mod_x, mod_y, W.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import socket
import sys
import tempfile
import time
from collections import Counter, deque
from pathlib import Path

import numpy as np

FARM = Path(__file__).resolve().parent
SS = FARM.parent
PYSCRIPTS = SS.parents[1]
for p in (str(SS), str(PYSCRIPTS)):
    if p not in sys.path:
        sys.path.insert(0, p)
import snake_scan as SSN        # noqa: E402
import edge_pass as EP          # noqa: E402
import lms_pass as LP           # noqa: E402

VERSION = 1
COSMIC_RUN = 23744
KEEP_SHARE = 0.60       # block kept when its large-peak argmax share >= this (wave-1 dwell test: 0.90, weak beam 0.75)
KEEP_MEV_LO = 650.0     # 3x3 MeV (W of the run) window of the stored events
KEEP_MEV_HI = 4100.0
KEEP_EVERY = 2          # stored events: event_number % KEEP_EVERY == 0 (counters stay exact); 1 = all, about 2x the size
GUESS_N = 400           # window (events with a large peak) of the running beam module
GUESS_WARM = 100
GUESS_MIN_SHARE = 0.6
DEFAULT_MAP = FARM / "hycal_map_til_run023849.json"


def run_of(path) -> int:
    m = re.search(r"prad_0*(\d+)\.evio", str(path))
    if not m:
        raise ValueError(f"cannot read the run number from {path}")
    return int(m.group(1))


def _top(counts, n):
    return SSN._top(counts, n)


class RunningModule:
    """Mode of the argmax of the last `n` large-peak events; None until `warm` of them and while the mode share < min_share."""

    def __init__(self, n=GUESS_N, warm=GUESS_WARM, min_share=GUESS_MIN_SHARE):
        self.q, self.c, self.n, self.warm, self.min_share = deque(), Counter(), n, warm, min_share
        self.cur = None

    def push(self, m):
        self.q.append(m)
        self.c[m] += 1
        if len(self.q) > self.n:
            o = self.q.popleft()
            self.c[o] -= 1
            if not self.c[o]:
                del self.c[o]
        if len(self.q) < self.warm:
            self.cur = None
            return
        k = max(self.c, key=self.c.get)
        self.cur = k if self.c[k] >= self.min_share * len(self.q) else None


def scan_pass(evio, map_path, calib, max_events=0, keep_every=KEEP_EVERY) -> dict:
    from prad2py import dec
    t0 = time.time()
    pipe = SSN.open_pipeline(map_path, str(evio))
    hy, wa = pipe.hycal, pipe.wave_ana
    ch = SSN.open_evio(evio, pipe)
    nb = EP.neighbours(map_path)
    names = sorted(nb, key=lambda n: int(n[1:]))
    mid = {n: i for i, n in enumerate(names)}
    W = EP.load_calib(calib)
    wvec = {m: np.array([W[n] for n in nb[m]] + [0.0] * (9 - len(nb[m])), np.float64) for m in names}   # edge modules have < 9 members: padded
    posof = {m: {n: k for k, n in enumerate(nb[m])} for m in names}
    cidx = {m: nb[m].index(m) for m in names}
    cache = {}

    def modname(crate, slot, c):
        k = (crate, slot, c)
        if k not in cache:
            m = hy.module_by_daq(crate, slot, c)
            cache[k] = m.name if (m is not None and m.is_pwo4()) else None
        return cache[k]

    run = RunningModule()
    blocks, epics, cnt = [], [], {}
    buf = dict(evn=[], ctr=[], I=[], npk=[])
    kept = dict(evn=[], ctr=[], I=[], npk=[])
    cur = None
    nsc = n_phys = 0

    def flush():
        nonlocal cur
        if cur is not None and cur["n"] > 0:
            top, share = _top(cur["counts_hi"], cur["n_hi"])
            if cur["n_hi"] and share >= KEEP_SHARE and top in mid:
                t = mid[top]
                for j, c in enumerate(buf["ctr"]):
                    if c == t:
                        for k in kept:
                            kept[k].append(buf[k][j])
            blocks.append(cur)
        cur = None
        for v in buf.values():
            v.clear()

    done = False
    while not done and ch.read() == dec.Status.success:
        if not ch.scan():
            continue
        et = ch.get_event_type()
        if et == dec.EventType.Epics:
            txt = ch.extract_epics_text()
            if txt:
                ts = blocks[-1]["ts_last"] if cur is None and blocks else (cur["ts_last"] if cur else 0)
                vals = {k: v for k, v in SSN.parse_epics_text(txt).items() if k.startswith("ptrans")}
                epics.append(dict(idx=nsc, ts=ts, vals=vals))
            continue
        if et != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            ch.select_event(i)
            inf = ch.info()
            n_phys += 1
            if int(inf.trigger_bits) & (SSN.TBIT_LMS | SSN.TBIT_ALPHA):
                continue
            evn = int(inf.event_number)
            if cur is None:
                cur = dict(ev_first=evn, ts_first=int(inf.timestamp), n=0, n_hi=0, counts={}, counts_hi={})
            cur["ev_last"] = evn
            cur["ts_last"] = int(inf.timestamp)
            cur["n"] += 1
            nsc += 1
            g = run.cur
            if g is not None:
                pos = posof[g]
                vals = np.zeros(9, np.float32)
                npk = np.zeros(9, np.uint8)
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
                        if g is not None:
                            k = pos.get(nm)
                            if k is not None:
                                _, _, peaks = wa.analyze(np.asarray(cd.samples, np.uint16))
                                npk[k] = min(len(peaks), 255)
                                if peaks:
                                    vals[k] = peaks[0].integral if len(peaks) == 1 else sum(p.integral for p in peaks)
            if best is not None:
                cur["counts"][best] = cur["counts"].get(best, 0) + 1
                if best_v >= SSN.PEAK_HI:
                    cur["n_hi"] += 1
                    cur["counts_hi"][best] = cur["counts_hi"].get(best, 0) + 1
                    run.push(best)
            if g is not None and npk[cidx[g]] >= 1:
                key = (len(blocks), mid[g])
                row = cnt.setdefault(key, [0, 0, 0])
                row[0] += 1
                row[1] += int(not (npk >= 2).any())
                row[2] += int(npk[cidx[g]] == 1)
                mev = float((vals.astype(np.float64) * wvec[g]).sum())
                if KEEP_MEV_LO <= mev <= KEEP_MEV_HI and evn % keep_every == 0:
                    buf["evn"].append(evn); buf["ctr"].append(mid[g]); buf["I"].append(vals); buf["npk"].append(npk)
            if cur["n"] >= 5000:
                flush()
            if max_events and nsc >= max_events:
                done = True
                break
    flush()
    if len(blocks) > 1 and blocks[-1]["n"] < 1000:             # fold as snake_scan._scan_blocks
        a, b = blocks[-2], blocks.pop()
        a["ev_last"], a["ts_last"], a["n"] = b["ev_last"], b["ts_last"], a["n"] + b["n"]
        a["n_hi"] += b["n_hi"]
        for key in ("counts", "counts_hi"):
            for k, v in b[key].items():
                a[key][k] = a[key].get(k, 0) + v
        nb_ = len(blocks)                                       # the folded block now has index nb_ - 1
        folded = {}
        for (bi, mi), v in cnt.items():
            k2 = (min(bi, nb_ - 1), mi)
            r = folded.setdefault(k2, [0, 0, 0])
            for j in range(3):
                r[j] += v[j]
        cnt = folded
    out = dict(
        mod_names=np.array(names),
        b_evfirst=np.array([b["ev_first"] for b in blocks], np.int64), b_evlast=np.array([b["ev_last"] for b in blocks], np.int64),
        b_tsfirst=np.array([b["ts_first"] for b in blocks], np.int64), b_tslast=np.array([b["ts_last"] for b in blocks], np.int64),
        b_n=np.array([b["n"] for b in blocks], np.int32), b_nhi=np.array([b["n_hi"] for b in blocks], np.int32),
        b_counts_hi_json=np.array(json.dumps([b["counts_hi"] for b in blocks])),
        epics_json=np.array(json.dumps(epics)),
        cnt=np.array([[bi, mi] + v for (bi, mi), v in sorted(cnt.items())], np.int32).reshape(-1, 5),
        ev_evn=np.array(kept["evn"], np.uint32), ev_ctr=np.array(kept["ctr"], np.uint16),
        ev_I=np.array(kept["I"], np.float32).reshape(-1, 9), ev_npk=np.array(kept["npk"], np.uint8).reshape(-1, 9))
    return dict(arrays=out, n_physics=n_phys, n_scanned=nsc, seconds=time.time() - t0)


def lms_arrays(evio, map_path, tmp_dir) -> dict:
    p = LP.lms_pass(evio, str(map_path), cache=Path(tmp_dir) / "lms.npz", force=True)
    z = np.load(p, allow_pickle=False)
    return {"L_" + k: z[k] for k in z.files}


def cosmic_arrays(evio, map_path, max_events=0) -> dict:
    import cosmic_mip as CM
    r = CM.scan_file(evio, map_path, max_events)
    out = dict(hits_module=r["hits"]["module"], hits_peak_adc=r["hits"]["peak_adc"], hits_integral_adc=r["hits"]["integral_adc"],
               hits_column_len=r["hits"]["column_len"], counts_json=np.array(json.dumps(r["counts"])),
               trig_hi_json=np.array(json.dumps({str(k): v for k, v in r["trig_hi"].items()})),
               names=np.array(r["names"]), mod_x=r["mod_x"], mod_y=r["mod_y"], W=r["W"])
    nW = len(r["names"])
    for k in ("ts", "ev"):
        out["lms_" + k] = np.asarray(r["lms"][k], np.int64)
        out["alpha_" + k] = np.asarray(r["alpha"][k], np.int64)
    out["lms_I"] = np.asarray(r["lms"]["I"], np.float32).reshape(-1, nW)
    out["lms_ref"] = np.asarray(r["lms"]["ref"], np.float32).reshape(-1, 3)
    out["alpha_ref"] = np.asarray(r["alpha"]["ref"], np.float32).reshape(-1, 3)
    return dict(arrays=out, n_physics=r["counts"].get("n_all", 0), n_scanned=r["counts"].get("n_all", 0))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("evio")
    ap.add_argument("outdir")
    ap.add_argument("--map", default=str(DEFAULT_MAP))
    ap.add_argument("--calib", default="", help="calibration JSON (default: the run's file of runinfo/general.json)")
    ap.add_argument("--max-events", type=int, default=0, help="stop after N scanned events (testing)")
    ap.add_argument("--keep-every", type=int, default=KEEP_EVERY, help="store every N-th event (by event number) per-event; counters stay exact")
    ap.add_argument("--mode", choices=("auto", "scan", "cosmic"), default="auto")
    ap.add_argument("--force", action="store_true", help="overwrite an existing product")
    a = ap.parse_args(argv)
    evio = Path(a.evio)
    if not evio.exists():
        print(f"[farm_pass] ERROR: {evio} does not exist (not on /cache yet? run jcache.sh --check)", file=sys.stderr)
        return 2
    run = run_of(evio)
    out = Path(a.outdir) / str(run) / (evio.name + ".npz")
    if out.exists() and not a.force:
        print(f"[farm_pass] {out} exists, skipping (--force to redo)")
        return 0
    map_path = Path(a.map).resolve()          # the C++ side resolves a relative map under database/
    if not map_path.exists():
        print(f"[farm_pass] ERROR: map {map_path} not found", file=sys.stderr)
        return 2
    mode = a.mode if a.mode != "auto" else ("cosmic" if run == COSMIC_RUN else "scan")
    out.parent.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    calib = ""
    with tempfile.TemporaryDirectory(prefix="farm_pass_") as tmp:
        if mode == "scan":
            calib = a.calib or str(EP.resolve_calib(run))
            r = scan_pass(evio, map_path, calib, a.max_events, a.keep_every)
            arrays = r["arrays"]
            arrays.update(lms_arrays(evio, map_path, tmp))
        else:
            r = cosmic_arrays(evio, map_path, a.max_events)
            arrays = r["arrays"]
        meta = dict(file=evio.name, run=run, mode=mode, map=map_path.name, calib=calib, n_physics=r["n_physics"],
                    n_scanned=r["n_scanned"], seconds=round(time.time() - t0, 1), version=VERSION,
                    host=socket.gethostname(), max_events=a.max_events, keep_every=a.keep_every if mode == "scan" else 0,
                    constants=dict(KEEP_SHARE=KEEP_SHARE, KEEP_MEV_LO=KEEP_MEV_LO, KEEP_MEV_HI=KEEP_MEV_HI, GUESS_N=GUESS_N,
                                   GUESS_WARM=GUESS_WARM, GUESS_MIN_SHARE=GUESS_MIN_SHARE, PEAK_HI=SSN.PEAK_HI))
        arrays["meta_json"] = np.array(json.dumps(meta))
        part = out.with_name("." + out.name[:-4] + ".partial.npz")
        try:
            np.savez_compressed(part, **arrays)
            os.replace(part, out)
        finally:
            if part.exists():
                part.unlink()
    print(f"[farm_pass] {out} {out.stat().st_size / 1024:.0f} KB, mode {mode}, {meta['n_scanned']} events, "
          f"{meta['seconds']:.0f} s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
