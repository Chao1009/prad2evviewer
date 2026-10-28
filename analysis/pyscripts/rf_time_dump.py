#!/usr/bin/env python3
"""rf_time_dump.py - per-event CEBAF RF TDC hits (one split -> one npz), companion of module_timing_dump.py.

PURPOSE
    Records, for every physics-type event of an EVIO split, the event number, trigger bits,
    TI timestamp and the leading-edge times of the divided CEBAF RF reference (V1190 TDC,
    ROC 0x40, slot 16, channels 0 = A and 8 = B), so that plot_module_timing_rf.py can join
    them to the module timing dump by (run, split, event number).

DEFINITIONS
    RF hit      leading edge (edge == 0) of ROC 0x40 / slot 16 / channel 0 or 8, value * TDC_LSB_NS
                (prad2py.dec.TDC_LSB_NS = 23.436 ps), ns within the TDC trigger-matching window.
                Same selection as tdc::RfTimeDecoder::Extract, but without its 16-hit cap.  The
                recorded edges are the 4.008 ns RF divided by 32 (spacing ~128.25 ns).
    Events      every event of a physics-type record (no trigger or decode-quality filter);
                read with the lazy EvChannel accessors info() and tdc() only (no FADC decode).

OUTPUT  DIR/<run:06d>_<split:05d>.npz  (np.savez_compressed, written atomically)
    per event   ev_num int64, trig uint32, ts uint64 (TI timestamp, 4 ns ticks),
                rf_off int32 (n_ev + 1; hits of event i are rf_*[rf_off[i]:rf_off[i+1]], in
                decoder order)
    per hit     rf_t float32 ns, rf_ch int8 (0 = A, 8 = B)
    metadata    run, split, evio, tdc_lsb_ns, n_read, n_phys, n_hits_a, n_hits_b, elapsed_s,
                script_version

USAGE
    python3 rf_time_dump.py FILE [FILE ...] --out DIR
    python3 rf_time_dump.py --all --parallel 12 [--dump-dir DIR] [--out DIR]
    --all takes the split list from the module timing dump directory (every <run>_<split>.npz
    there, run 24578 skipped).  Existing outputs are skipped (resume).
"""
import argparse
import glob
import os
import re
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT_VERSION = "rf_time_dump 1.0 (20261025)"
EVIO_DIR = "/Users/L00338853/work/PRad/data/evio"
DEF_DUMP = os.path.join(HERE, "output", "module_timing_20261025", "dump")
DEF_OUT = os.path.join(HERE, "output", "module_timing_20261025", "rf")
SKIP_RUNS = (24578,)


def run_split(path):
    run = int(re.search(r"prad_0*(\d+)", os.path.basename(path)).group(1))
    return run, int(path.rsplit(".", 1)[1])


def process_split(path, out_dir):
    sys.path.insert(0, HERE)
    import _common as C
    from prad2py import dec

    run, split = run_split(path)
    out_path = os.path.join(out_dir, f"{run:06d}_{split:05d}.npz")
    t0 = time.monotonic()
    p = C.setup_pipeline(evio_path=path)
    roc, slot = int(dec.RF_ROC_TAG), int(dec.RF_SLOT)
    cha, chb = int(dec.RF_CH_A), int(dec.RF_CH_B)
    lsb = float(dec.TDC_LSB_NS)
    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    if ch.open_auto(path) != dec.Status.success:
        raise RuntimeError(f"cannot open {path}")
    ev_num, trig, ts, off, tt, cc = [], [], [], [0], [], []
    n_read = n = 0
    while ch.read() == dec.Status.success:
        n_read += 1
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            ch.select_event(i)
            info = ch.info()
            h = ch.tdc().hits_numpy
            k = ((h["roc_tag"] == roc) & (h["slot"] == slot) & (h["edge"] == 0)
                 & ((h["channel"] == cha) | (h["channel"] == chb)))
            h = h[k]
            ev_num.append(int(info.event_number))
            trig.append(int(info.trigger_bits))
            ts.append(int(info.timestamp))
            tt.append((h["value"] * lsb).astype(np.float32))
            cc.append(h["channel"].astype(np.int8))
            n += len(h)
            off.append(n)
    rf_t = np.concatenate(tt) if tt else np.zeros(0, np.float32)
    rf_ch = np.concatenate(cc) if cc else np.zeros(0, np.int8)
    elapsed = time.monotonic() - t0
    out = dict(ev_num=np.asarray(ev_num, np.int64), trig=np.asarray(trig, np.uint32),
               ts=np.asarray(ts, np.uint64), rf_off=np.asarray(off, np.int32), rf_t=rf_t, rf_ch=rf_ch,
               run=run, split=split, evio=path, tdc_lsb_ns=lsb, n_read=n_read, n_phys=len(ev_num),
               n_hits_a=int((rf_ch == cha).sum()), n_hits_b=int((rf_ch == chb).sum()),
               elapsed_s=elapsed, script_version=SCRIPT_VERSION)
    os.makedirs(out_dir, exist_ok=True)
    tmp = out_path + ".tmp"
    with open(tmp, "wb") as fh:
        np.savez_compressed(fh, **out)
    os.replace(tmp, out_path)
    n_ev = max(len(ev_num), 1)
    msg = (f"[done] run {run} split {split:05d}: {n_read} read, {len(ev_num)} physics, "
           f"RF A {out['n_hits_a'] / n_ev:.2f} / B {out['n_hits_b'] / n_ev:.2f} hits/ev, "
           f"{int((np.diff(off) > 0).sum())} ev with RF, {elapsed:.1f} s -> {out_path}")
    return msg


def _job(a):
    try:
        return process_split(*a)
    except Exception as e:  # report and continue with the other splits
        return f"[FAIL] {a[0]}: {e!r}"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="*", help="EVIO split files")
    ap.add_argument("--out", default=DEF_OUT, metavar="DIR")
    ap.add_argument("--all", action="store_true", help="every split present in --dump-dir")
    ap.add_argument("--dump-dir", default=DEF_DUMP, metavar="DIR")
    ap.add_argument("--parallel", type=int, default=1, metavar="N")
    a = ap.parse_args()
    files = list(a.files)
    if a.all:
        for fn in sorted(glob.glob(os.path.join(a.dump_dir, "*.npz"))):
            m = re.match(r"(\d{6})_(\d{5})\.npz$", os.path.basename(fn))
            if m and int(m.group(1)) not in SKIP_RUNS:
                files.append(os.path.join(EVIO_DIR, f"prad_{m.group(1)}.evio.{m.group(2)}"))
    if not files:
        ap.error("no input files (give FILE ... or --all)")
    todo = [f for f in files
            if not os.path.exists(os.path.join(a.out, "%06d_%05d.npz" % run_split(f)))]
    print(f"[driver] {len(files)} split(s), {len(files) - len(todo)} already done, {len(todo)} to run",
          flush=True)
    t0 = time.monotonic()
    if a.parallel <= 1:
        res = []
        for f in todo:
            res.append(_job((f, a.out)))
            print(res[-1], flush=True)
    else:
        import multiprocessing as mp
        with mp.get_context("spawn").Pool(a.parallel) as pool:
            res = []
            for r in pool.imap_unordered(_job, [(f, a.out) for f in todo]):
                print(r, flush=True)
                res.append(r)
    n_fail = sum(r.startswith("[FAIL]") for r in res)
    print(f"[driver] finished, {n_fail} failed, wall {time.monotonic() - t0:.1f} s", flush=True)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
