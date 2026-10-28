#!/usr/bin/env python3
"""module_timing_dump.py - per-event module / veto timing extraction (one split -> one npz).

PURPOSE
    Light first pass over EVIO splits that records, for every decoded physics-type
    event, the trigger bits, the veto-scintillator pulses (V1..V4) and, for every
    PbWO4 module with a pulse, its dominant pulse (time, energy).  A later plotting
    script builds t_module - t_veto and t_module - t_eventref per module from these
    arrays; nothing is plotted or cut here (no trigger filter, no time window).

DEFINITIONS
    PbWO4 pulse time   mf (matched-filter) reconstruction with the reference spline
                       templates (mf_pulse.MFReconstructor, nsigma 5, min_amp 5,
                       kmax 4, r_model 0.3 = the R2 config).  time = res.p_time =
                       fitted onset + template t50, i.e. the 50 % leading-edge time
                       (ns).  No per-module time offsets.
    PbWO4 pulse energy Eu = A * I_int(phase) * f(A) * kappa (kappa = 1.0), converted
                       with mod.energize(Eu).  NO LMS correction and NO per-crystal
                       gains: the energy (MeV) is only a selection variable.
    dominant pulse     the pulse with the largest Eu anywhere in the readout window
                       (first one in onset order on ties).  m_npk = number of pulses
                       in the channel; m_inwin = dominant time in the reference
                       pre-window 100 < t < 200 ns; m_fb = the channel fell back to
                       WaveAnalyzer peaks (mf rejected it, or nsamples differs from
                       the block's); then (time, Eu) = (peak.time, peak.integral).
                       Channels without any pulse get no row.
    veto pulse time    WaveAnalyzer CFD leading-edge pickoff (peak.time, ns); every
                       pulse of V1..V4 is stored with integral, height and overflow.

OUTPUT  DIR/<run:06d>_<split:05d>.npz  (written atomically)
    per event   ev_num int64, trig uint32            (n_ev entries)
    PbWO4       m_ev int32 (index into event arrays), m_idx int16 (p.hycal module
                index), m_t float32 ns, m_E float32 MeV, m_npk int8, m_inwin bool,
                m_fb bool
    veto        v_ev int32, v_ch int8 (1..4), v_t float32, v_int float32,
                v_h float32, v_ovf bool
    metadata    run, split, module_names (p.hycal index order), module_is_pwo4,
                templates, calib_file, n_read, n_phys, elapsed_s, script_version

USAGE
    python3 module_timing_dump.py FILE [FILE ...] --out DIR [--max-events N]
    python3 module_timing_dump.py --all --parallel 12 --out DIR [--runs 25320 ...]
        [--templates-root DIR]     (default output/templates_spline, DIR/<type>/templates_<type>.npz)
    Each file is processed alone (only that split).  Existing npz are skipped
    (resume).  With --parallel, each split runs as a subprocess, log in DIR/logs/.
    Run-type settings (templates, calibration file) come from the RUN_CFG table.
"""
import argparse
import glob
import os
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT_VERSION = "module_timing_dump 1.0 (20261025)"
EVIO_DIR = "/Users/L00338853/work/PRad/data/evio"
CALIB_DIR = "/Users/L00338853/work/PRad/prad2evviewer/database/calibration"
PRE_LO, PRE_HI = 100.0, 200.0
BLOCK_EVENTS = 256

TEMPLATES_ROOT = "output/templates_spline"     # --templates-root default (relative to pyscripts)

# run -> (template npz relative to the templates root, calibration file name)
_T = "{0}/templates_{0}.npz"
_TYPES = {
    "3p49": (_T.format("3p49GeV_H2"), "calibration_factor_fadc_3p5.json",
             (24386, 24503, 24512)),
    "2p24": (_T.format("2p24GeV_H2"), "calibration_factor_2p2_fadc.json",
             (25228, 25257, 25260, 25299, 25320)),
    "0p73": (_T.format("0p73GeV_H2"), "calibration_factor_fadc_0p7.json",
             (24590, 24603)),
    "x17": (_T.format("2p24GeV_Ta_X17"), "x17_25825_fadc.json",
            (25800, 26138, 26230)),
}
RUN_CFG = {r: (t, c) for t, c, runs in _TYPES.values() for r in runs}


def split_of(path):
    return int(path.rsplit(".", 1)[1])


def out_name(path, out_dir):
    import re
    run = int(re.search(r"prad_0*(\d+)", os.path.basename(path)).group(1))
    return run, os.path.join(out_dir, f"{run:06d}_{split_of(path):05d}.npz")


# ----------------------------------------------------------------------------
# one split
# ----------------------------------------------------------------------------
def process_split(path, out_dir, max_events=0, templates_root=None):
    import _common as C
    import mf_pulse
    from prad2py import dec

    run, out_path = out_name(path, out_dir)
    if run not in RUN_CFG:
        sys.exit(f"[ERROR] run {run} not in RUN_CFG")
    tmpl = os.path.join(templates_root or os.path.join(HERE, TEMPLATES_ROOT), RUN_CFG[run][0])
    calib = os.path.join(CALIB_DIR, RUN_CFG[run][1])
    t_start = time.monotonic()

    # setup_pipeline takes a single file path unchanged (no split expansion);
    # only that file is opened below.
    p = C.setup_pipeline(evio_path=path, hc_calib_file=calib)
    hy = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    names = np.array([m.name for m in hy])
    pwo4 = [bool(m.is_pwo4()) for m in hy]
    wcfg = p.wave_ana.cfg
    overflow = int(wcfg.overflow)
    bank, info = mf_pulse.build_bank(
        [m.name for m in hy], pwo4, tmpl,
        dict(smooth_order=int(wcfg.smooth_order),
             int_tail_ratio=float(wcfg.int_tail_ratio),
             tail_break_n=int(wcfg.tail_break_n)),
        column=mf_pulse.resolve_column(tmpl, None, prefer="25320"),
        int_tail_ratio=None)
    rec = mf_pulse.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4,
                                   r_model=0.3, ped_nsamples=int(wcfg.ped_nsamples))
    print(f"[setup] run {run} templates {tmpl} ({info['kind']}) calib {calib}",
          flush=True)

    # outputs (lists of per-block / per-row chunks)
    ev_num, trig = [], []
    M = {k: [] for k in ("ev", "idx", "t", "E", "npk", "inwin", "fb")}
    M_DT = dict(ev=np.int32, idx=np.int16, t=np.float32, E=np.float32,
                npk=np.int8, inwin=bool, fb=bool)
    M_chunks = []            # compact dtype arrays; M lists are packed per block

    def pack_m():
        if M["ev"]:
            M_chunks.append({k: np.asarray(v, dtype=M_DT[k]) for k, v in M.items()})
            for v in M.values():
                v.clear()
    V = {k: [] for k in ("ev", "ch", "t", "int", "h", "ovf")}

    def put_m(ev, idx, t, E, npk, fb):
        M["ev"].append(ev); M["idx"].append(idx); M["t"].append(t)
        M["E"].append(E); M["npk"].append(npk)
        M["inwin"].append(PRE_LO < t < PRE_HI); M["fb"].append(fb)

    def put_wa_dominant(ev, mod, peaks):      # WaveAnalyzer peaks -> dominant
        if not peaks:
            return
        best = peaks[0]
        for pk in peaks[1:]:
            if pk.integral > best.integral:
                best = pk
        put_m(ev, mod.index, best.time, mod.energize(best.integral),
              len(peaks), True)

    # block buffer of PbWO4 rows
    st = dict(Y=None, mid=None, n=0, ns=None)
    row_ev = []

    def flush():
        n = st["n"]
        if n == 0:
            return
        Y, mid = st["Y"][:n], st["mid"][:n]
        res = rec.run(Y, mid, overflow)
        pc = res.p_ch
        Eu = mf_pulse.pulse_energy_units(bank, mid[pc], res.p_t0, res.p_amp,
                                         1.0, None)
        tl = res.p_time
        starts = np.searchsorted(pc, np.arange(n + 1))
        # dominant pulse per row: group by channel, largest Eu first (stable)
        order = np.lexsort((-Eu, pc))
        first = order[starts[:-1][starts[1:] > starts[:-1]]] if len(pc) else []
        for j in first:
            r = int(pc[j])
            mod = hy[int(mid[r])]
            put_m(row_ev[r], mod.index, float(tl[j]),
                  mod.energize(float(Eu[j])), int(starts[r + 1] - starts[r]),
                  False)
        for r in np.nonzero(res.wa_fb)[0]:
            _, _, pks = p.wave_ana.analyze(Y[r])
            put_wa_dominant(row_ev[r], hy[int(mid[r])], list(pks))
        st["n"] = 0
        row_ev.clear()
        pack_m()

    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    n_read = n_phys = n_evt = 0
    if ch.open_auto(path) != dec.Status.success:
        sys.exit(f"[ERROR] cannot open {path}")
    done = False
    while not done and ch.read() == dec.Status.success:
        n_read += 1
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            e = d["event"]
            n_phys += 1
            ev = n_evt
            n_evt += 1
            ev_num.append(int(e.info.event_number))
            trig.append(int(e.info.trigger_bits))
            for ri in range(e.nrocs):
                roc = e.roc(ri)
                if not roc.present:
                    continue
                crate = p.crate_map.get(roc.tag)
                if crate is None:
                    continue
                for s in roc.present_slots():
                    slot = roc.slot(s)
                    for c in slot.present_channels():
                        mod = p.hycal.module_by_daq(crate, s, c)
                        if mod is None:
                            continue
                        cd = slot.channel(c)
                        if cd.nsamples <= 0:
                            continue
                        if mod.type.name == "Veto":
                            _, _, pks = p.wave_ana.analyze(cd.samples)
                            for pk in pks:
                                V["ev"].append(ev); V["ch"].append(int(mod.name[1:]))
                                V["t"].append(pk.time); V["int"].append(pk.integral)
                                V["h"].append(pk.height)
                                V["ovf"].append(bool(pk.overflow))
                            continue
                        if not (mod.is_hycal() and pwo4[mod.index]):
                            continue
                        if st["ns"] is None:
                            st["ns"] = cd.nsamples
                            cap = BLOCK_EVENTS * 256
                            st["Y"] = np.zeros((cap, st["ns"]), dtype=np.uint16)
                            st["mid"] = np.zeros(cap, dtype=np.int64)
                        if cd.nsamples != st["ns"]:      # odd length: WaveAnalyzer
                            _, _, pks = p.wave_ana.analyze(cd.samples)
                            put_wa_dominant(ev, mod, list(pks))
                            continue
                        r = st["n"]
                        if r >= len(st["mid"]):
                            st["Y"] = np.concatenate([st["Y"], np.zeros_like(st["Y"])])
                            st["mid"] = np.concatenate([st["mid"],
                                                        np.zeros_like(st["mid"])])
                        st["Y"][r] = cd.samples
                        st["mid"][r] = mod.index
                        row_ev.append(ev)
                        st["n"] = r + 1
            if n_evt % BLOCK_EVENTS == 0:
                flush()
            if max_events > 0 and n_phys >= max_events:
                done = True
                break
    flush()
    pack_m()
    elapsed = time.monotonic() - t_start

    # assemble (rows sorted by event; the odd-length/fallback rows interleave)
    Mc = {k: (np.concatenate([c[k] for c in M_chunks]) if M_chunks
              else np.zeros(0, dtype=M_DT[k])) for k in M_DT}
    m_ev = Mc["ev"]
    o = np.argsort(m_ev, kind="stable")
    ev_a = np.asarray(ev_num, dtype=np.int64)
    out = dict(
        ev_num=ev_a, trig=np.asarray(trig, dtype=np.uint32),
        m_ev=m_ev[o],
        m_idx=Mc["idx"][o],
        m_t=Mc["t"][o],
        m_E=Mc["E"][o],
        m_npk=Mc["npk"][o],
        m_inwin=Mc["inwin"][o],
        m_fb=Mc["fb"][o],
        v_ev=np.asarray(V["ev"], dtype=np.int32),
        v_ch=np.asarray(V["ch"], dtype=np.int8),
        v_t=np.asarray(V["t"], dtype=np.float32),
        v_int=np.asarray(V["int"], dtype=np.float32),
        v_h=np.asarray(V["h"], dtype=np.float32),
        v_ovf=np.asarray(V["ovf"], dtype=bool),
        run=run, split=split_of(path), module_names=names,
        module_is_pwo4=np.asarray(pwo4, dtype=bool),
        templates=tmpl, calib_file=calib, n_read=n_read, n_phys=n_phys,
        elapsed_s=elapsed, script_version=SCRIPT_VERSION)
    os.makedirs(out_dir, exist_ok=True)
    tmp = out_path + ".tmp"
    with open(tmp, "wb") as fh:
        np.savez(fh, **out)
    os.replace(tmp, out_path)
    n_veto_ev = len(np.unique(out["v_ev"]))
    print(f"[done] run {run} split {split_of(path):05d}: {n_read} read, "
          f"{n_phys} physics, {n_veto_ev} with veto pulse, "
          f"{len(out['m_ev'])} module rows, {elapsed:.1f} s "
          f"({n_phys / max(elapsed, 1e-9):.0f} ev/s) -> {out_path}", flush=True)


# ----------------------------------------------------------------------------
# driver
# ----------------------------------------------------------------------------
def all_splits(runs):
    files = []
    for r in sorted(runs):
        files += sorted(glob.glob(os.path.join(EVIO_DIR, f"prad_{r:06d}.evio.[0-9]*")))
    return files


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("files", nargs="*", help="EVIO split files (each processed alone)")
    ap.add_argument("--out", required=True, metavar="DIR")
    ap.add_argument("--max-events", type=int, default=0,
                    help="stop after N physics-type events per split (0 = all)")
    ap.add_argument("--all", action="store_true",
                    help="every local split of the RUN_CFG runs")
    ap.add_argument("--runs", type=int, nargs="+", metavar="RUN",
                    help="restrict --all to these runs")
    ap.add_argument("--parallel", type=int, default=1, metavar="N",
                    help="run N splits concurrently as subprocesses")
    ap.add_argument("--templates-root", default=TEMPLATES_ROOT, metavar="DIR",
                    help="spline template root (DIR/<type>/templates_<type>.npz)")
    a = ap.parse_args()
    troot = (os.path.join(HERE, TEMPLATES_ROOT) if a.templates_root == TEMPLATES_ROOT
             else os.path.abspath(a.templates_root))

    files = list(a.files)
    if a.all:
        runs = a.runs or list(RUN_CFG)
        files += all_splits(runs)
    elif a.runs:
        files = [f for f in files if out_name(f, a.out)[0] in a.runs]
    if not files:
        ap.error("no input files (give FILE ... or --all)")
    todo = [f for f in files if not os.path.exists(out_name(f, a.out)[1])]
    print(f"[driver] {len(files)} split(s), {len(files) - len(todo)} already done, "
          f"{len(todo)} to run", flush=True)

    if a.parallel <= 1 or len(todo) <= 1:
        for f in todo:
            process_split(f, a.out, a.max_events, troot)
        return

    log_dir = os.path.join(a.out, "logs")
    os.makedirs(log_dir, exist_ok=True)
    running, queue, n_fail = [], list(todo), 0
    while queue or running:
        while queue and len(running) < a.parallel:
            f = queue.pop(0)
            run, op = out_name(f, a.out)
            log = open(os.path.join(log_dir, os.path.basename(op)[:-4] + ".log"), "w")
            cmd = [sys.executable, os.path.abspath(__file__), f, "--out", a.out,
                   "--max-events", str(a.max_events), "--templates-root", troot]
            running.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT),
                            f, log))
        time.sleep(1.0)
        for item in list(running):
            pr, f, log = item
            if pr.poll() is None:
                continue
            running.remove(item)
            log.close()
            ok = pr.returncode == 0
            n_fail += not ok
            with open(log.name) as fh:
                last = [ln for ln in fh.read().splitlines() if ln.startswith("[done]")]
            print(last[-1] if ok and last else f"[FAIL rc={pr.returncode}] {f}",
                  flush=True)
    print(f"[driver] finished, {n_fail} failed", flush=True)
    sys.exit(1 if n_fail else 0)


if __name__ == "__main__":
    main()
