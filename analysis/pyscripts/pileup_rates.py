#!/usr/bin/env python3
"""pileup_rates.py - real-data matched-filter (MF) pile-up rates for PbWO4 (dump + analyse).

PURPOSE
    How often does the MF find more than one pulse in a PbWO4 channel of real physics
    events, per run type, per cluster, per module and as a function of the peak-time
    separation and the peak-height ratio of the extra pulse.  Two stages:
      dump     one EVIO split -> one npz with EVERY MF pulse of every PbWO4 channel and the
               legacy clusters of the event (needs prad2py, minutes per split);
      analyse  npz -> tables, module maps, Delta-t x ratio distributions and report.md
               (numpy / matplotlib only; reruns without the MF).

DEFINITIONS
    MF                 mf_pulse.MFReconstructor on the best-50 spline bank, nsigma 5, min_amp 5,
                       kmax 4, r_model 0.3 (the R2/R3 config).  p_time = fitted onset + t50 (ns),
                       A = p_amp (ADC), E = mod.energize(A * I_int(phase) * f(A)) (MeV, no LMS,
                       no per-crystal gain: only a selection variable).
    event              physics-type EVIO event with (trig & 0xf00) != 0 and (trig & 0x3000000) == 0
                       (SSP accept bits; LMS flasher / alpha pulser excluded).  Denominator.
    fired channel      PbWO4 channel with >= 1 MF pulse.  Channels whose MF fit was rejected
                       (WaveAnalyzer fallback, flag fb, A = peak height) or whose length differs
                       are EXCLUDED from every pile-up rate and counted separately (n_fb).
    pile-up channel    fired channel with >= 2 MF pulses anywhere in the readout window.
    pulse 1            the IN-TIME pulse: the largest-E pulse with 100 < t < 200 ns (the pulse the
                       legacy clusterer uses).  Every other pulse k of the channel is an extra
                       pulse: dt = t_k - t_1 (signed; tables use |dt|), ratio = A_k / A_1 (> 1
                       possible: the extra pulse can be larger than the in-time one).  Channels
                       without an in-window pulse have no pulse 1 (counted separately, no dt).
    near pile-up       an extra pulse with |dt| <= 60 ns of pulse 1 (the ones that can distort the
                       in-time pulse), as opposed to "any" (any second pulse in the window).
    cluster            det.HyCalCluster, legacy mode (one hit per module = pulse 1), config from
                       reconstruction_config.json, hits from the MF + WaveAnalyzer for lead glass;
                       reconstruct_matched(); E = cluster hit energy (MeV), kept if >= config
                       min_cluster_energy.  Cluster pile-up = any member channel / the seed channel
                       is a pile-up channel (any or near, as above).
    MF-incomplete      |dt| < 16 ns or ratio < 0.3: the secondary search (mf_pulse.py ~1054-1091:
                       12 ns MIN_SEP, 0.3*A_parent gate only after the parent onset, argmax then
                       stop) is biased low there (output/templates_wassim_20261027/why_hard_cells.md).
                       Cells carry a trailing '*' in the TSV/report.
    RF comb            fold dt to (-P/2, P/2], P = 4.0077 ns; on-comb = |fold| < 0.5 ns; random
                       expectation 1.0/4.0077 = 25 %.  dt is within one channel: offsets cancel.
    run groups         0p73, 2p24, 3p49 (H2, pooled over the runs of the type), c2p24 (run 25308),
                       Ta runs x17_25800 / x17_26138 / x17_26230 (never pooled).

OUTPUT
    dump     DIR/<run:06d>_<split:05d>.npz (atomic): per pulse pu_ev, pu_idx, pu_t, pu_A, pu_E, pu_fb;
             clusters cl_ev, cl_seed, cl_E, cl_x, cl_y, cl_nb and members cm_cl, cm_idx; n_sel,
             n_phys, module geometry/names, metadata.
    analyse  rates_by_type.tsv, rf_comb.tsv, cluster_rates.tsv, module_rates_<g>.tsv,
             module_map_<g>.png, dt_ratio_<g>.tsv, dt_ratio_distributions.pdf, report.md

USAGE
    python3 pileup_rates.py dump --sample --parallel 8 --out output/pileup_rates_20261027/dump
        [--max-events 40000]        (--sample = the built-in SAMPLE split list; or FILE ...)
    python3 pileup_rates.py analyse --dump-dir output/pileup_rates_20261027/dump
        --out output/pileup_rates_20261027
    Existing npz are skipped (resume).  With --parallel each split is a subprocess, log in
    DIR/logs/.
"""
import argparse
import glob
import os
import re
import subprocess
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
SCRIPT_VERSION = "pileup_rates 1.0 (20261027)"
EVIO_DIR = "/Users/L00338853/work/PRad/data/evio"
CALIB_DIR = "/Users/L00338853/work/PRad/prad2evviewer/database/calibration"
TEMPLATES_ROOT = "output/templates_spline_best50_20261026"
PRE_LO, PRE_HI = 100.0, 200.0
BLOCK_EVENTS = 256
TRIG_MASK, TRIG_EXCL = 0xF00, 0x3000000
RF_P = 4.0077
GRID_P = 4.000                          # FADC sample period (ns)
PAIR_BINS = [(16, 120), (120, 260), (260, 400)]   # |separation| bins of within-channel pulse pairs
NEAR_NS = 60.0
INC_DT, INC_RATIO = 16.0, 0.3          # MF-incomplete region

_T = "{0}/templates_{0}.npz"
# group -> (label, template slug, calibration file, runs)
GROUPS = {
    "0p73": ("0.73 GeV H2", "0p73GeV_H2", "calibration_factor_fadc_0p7.json", (24590, 24603)),
    "2p24": ("2.24 GeV H2", "2p24GeV_H2", "calibration_factor_2p2_fadc.json",
             (25228, 25257, 25260, 25299, 25320)),
    "3p49": ("3.49 GeV H2", "3p49GeV_H2", "calibration_factor_fadc_3p5.json", (24386, 24503, 24512)),
    "c2p24": ("2.24 GeV C", "2p24GeV_C", "calibration_factor_2p2_fadc.json", (25308,)),
    "x17_25800": ("Ta 25800", "2p24GeV_Ta_X17", "x17_25825_fadc.json", (25800,)),
    "x17_26138": ("Ta 26138", "2p24GeV_Ta_X17", "x17_25825_fadc.json", (26138,)),
    "x17_26230": ("Ta 26230", "2p24GeV_Ta_X17", "x17_25825_fadc.json", (26230,)),
}
ORDER = list(GROUPS)
RUN_GROUP = {r: g for g, v in GROUPS.items() for r in v[3]}
RUN_CFG = {r: (_T.format(v[1]), v[2]) for g, v in GROUPS.items() for r in v[3]}

# built-in sample: (run, split) pairs spread over the run(s) of each type
SAMPLE = [
    (24386, 20), (24503, 50), (24512, 30), (24386, 100),                  # 3.49 H2
    (25228, 50), (25260, 400), (25299, 700), (25320, 300),                # 2.24 H2
    (24590, 10), (24590, 80), (24603, 100), (24603, 250),                 # 0.73 H2
    (25308, 50), (25308, 150), (25308, 250), (25308, 400),                # 2.24 C
    (25800, 100), (25800, 400), (25800, 800),                             # Ta
    (26138, 99), (26138, 299), (26138, 399),
    (26230, 100), (26230, 600), (26230, 1300),
]


def split_of(path):
    return int(path.rsplit(".", 1)[1])


def out_name(path, out_dir):
    run = int(re.search(r"prad_0*(\d+)", os.path.basename(path)).group(1))
    return run, os.path.join(out_dir, f"{run:06d}_{split_of(path):05d}.npz")


# ----------------------------------------------------------------------------
# dump: one split
# ----------------------------------------------------------------------------
def process_split(path, out_dir, max_events=0, templates_root=None):
    import _common as C
    import mf_pulse
    from prad2py import dec, det

    run, out_path = out_name(path, out_dir)
    if run not in RUN_CFG:
        sys.exit(f"[ERROR] run {run} not in RUN_CFG")
    tmpl = os.path.join(templates_root or os.path.join(HERE, TEMPLATES_ROOT), RUN_CFG[run][0])
    calib = os.path.join(CALIB_DIR, RUN_CFG[run][1])
    t_start = time.monotonic()

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
    # legacy clusterer, config copied from the pipeline's (reconstruction_config.json)
    base_cfg = p.hc_clusterer.get_config()
    ccfg = det.HyCalClusterConfig()
    for f in ("min_module_energy", "min_center_energy", "min_cluster_energy",
              "min_cluster_size", "corner_conn", "split_iter", "least_split",
              "log_weight_thres"):
        setattr(ccfg, f, getattr(base_cfg, f))
    ccfg.seed_time_window = -1.0
    cl = det.HyCalCluster(p.hycal)
    cl.set_config(ccfg)
    print(f"[setup] run {run} templates {tmpl} ({info['kind']}) calib {calib} "
          f"min_cluster_energy {ccfg.min_cluster_energy}", flush=True)

    PU_DT = dict(ev=np.int32, idx=np.int16, t=np.float32, A=np.float32, E=np.float32, fb=bool)
    PU = {k: [] for k in PU_DT}
    CL = {k: [] for k in ("ev", "seed", "E", "x", "y", "nb")}
    CM = {"cl": [], "idx": []}
    ev_hits = {}            # event -> [(module index, E, t)]  (pulse-1 hit of every channel)
    n_cl = [0]
    st = dict(Y=None, mid=None, n=0, ns=None, blk_lo=0)
    row_ev = []
    cnt = dict(n_fb_rows=0, n_odd=0)

    def put_pulse(ev, idx, t, A, E, fb):
        PU["ev"].append(ev); PU["idx"].append(idx); PU["t"].append(t)
        PU["A"].append(A); PU["E"].append(E); PU["fb"].append(fb)

    def put_wa(ev, mod, peaks):          # WaveAnalyzer peaks: all stored (fb), best in window -> hit
        best, bi = None, -1.0
        for pk in peaks:
            E = mod.energize(pk.integral)
            if pwo4[mod.index]:
                put_pulse(ev, mod.index, float(pk.time), float(pk.height), E, True)
            if PRE_LO < pk.time < PRE_HI and pk.integral > bi:
                best, bi = (mod.index, E, float(pk.time)), pk.integral
        if best is not None:
            ev_hits.setdefault(ev, []).append(best)

    def cluster_block(ev_lo, ev_hi):
        for ev in range(ev_lo, ev_hi):
            hits = ev_hits.pop(ev, None)
            if not hits:
                continue
            cl.clear()
            for idx, E, t in hits:
                cl.add_hit(idx, E, t)
            cl.form_clusters()
            for r in cl.reconstruct_matched():
                row = n_cl[0]
                n_cl[0] += 1
                CL["ev"].append(ev)
                CL["seed"].append(p.hycal.module_by_id(r.hit.center_id).index)
                CL["E"].append(r.hit.energy); CL["x"].append(r.hit.x)
                CL["y"].append(r.hit.y); CL["nb"].append(r.hit.nblocks)
                for h in r.cluster.hits:
                    CM["cl"].append(row); CM["idx"].append(h.index)

    def flush(n_evt):
        n = st["n"]
        if n:
            Y, mid = st["Y"][:n], st["mid"][:n]
            res = rec.run(Y, mid, overflow)
            pc = res.p_ch
            Eu = mf_pulse.pulse_energy_units(bank, mid[pc], res.p_t0, res.p_amp, 1.0, None)
            tl, amp = res.p_time, res.p_amp
            fbrow = np.asarray(res.wa_fb, dtype=bool)
            best = {}
            for j in range(len(pc)):
                r = int(pc[j])
                if fbrow[r]:
                    continue
                mod = hy[int(mid[r])]
                E = mod.energize(float(Eu[j]))
                put_pulse(row_ev[r], mod.index, float(tl[j]), float(amp[j]), E, False)
                if PRE_LO < tl[j] < PRE_HI and (r not in best or E > best[r][1]):
                    best[r] = (mod.index, E, float(tl[j]))
            for r, h in best.items():
                ev_hits.setdefault(row_ev[r], []).append(h)
            for r in np.nonzero(fbrow)[0]:
                cnt["n_fb_rows"] += 1
                _, _, pks = p.wave_ana.analyze(Y[r])
                put_wa(row_ev[r], hy[int(mid[r])], list(pks))
            st["n"] = 0
            row_ev.clear()
        cluster_block(st["blk_lo"], n_evt)
        st["blk_lo"] = n_evt

    ev_num, trig = [], []
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
            tb = int(e.info.trigger_bits)
            if (tb & TRIG_MASK) and not (tb & TRIG_EXCL):
                ev = n_evt
                n_evt += 1
                ev_num.append(int(e.info.event_number))
                trig.append(tb)
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
                            if mod is None or not mod.is_hycal():
                                continue
                            cd = slot.channel(c)
                            if cd.nsamples <= 0:
                                continue
                            if not pwo4[mod.index]:               # lead glass: hit only
                                _, _, pks = p.wave_ana.analyze(cd.samples)
                                put_wa(ev, mod, list(pks))
                                continue
                            if st["ns"] is None:
                                st["ns"] = cd.nsamples
                                cap = BLOCK_EVENTS * 256
                                st["Y"] = np.zeros((cap, st["ns"]), dtype=np.uint16)
                                st["mid"] = np.zeros(cap, dtype=np.int64)
                            if cd.nsamples != st["ns"]:           # odd length: WaveAnalyzer
                                cnt["n_odd"] += 1
                                _, _, pks = p.wave_ana.analyze(cd.samples)
                                put_wa(ev, mod, list(pks))
                                continue
                            r = st["n"]
                            if r >= len(st["mid"]):
                                st["Y"] = np.concatenate([st["Y"], np.zeros_like(st["Y"])])
                                st["mid"] = np.concatenate([st["mid"], np.zeros_like(st["mid"])])
                            st["Y"][r] = cd.samples
                            st["mid"][r] = mod.index
                            row_ev.append(ev)
                            st["n"] = r + 1
                if n_evt % BLOCK_EVENTS == 0:
                    flush(n_evt)
            if max_events > 0 and n_phys >= max_events:
                done = True
                break
    flush(n_evt)
    elapsed = time.monotonic() - t_start

    PUa = {k: np.asarray(v, dtype=PU_DT[k]) for k, v in PU.items()}
    o = np.lexsort((PUa["t"], PUa["idx"], PUa["ev"]))
    out = {f"pu_{k}": v[o] for k, v in PUa.items()}
    out.update(
        cl_ev=np.asarray(CL["ev"], np.int32), cl_seed=np.asarray(CL["seed"], np.int16),
        cl_E=np.asarray(CL["E"], np.float32), cl_x=np.asarray(CL["x"], np.float32),
        cl_y=np.asarray(CL["y"], np.float32), cl_nb=np.asarray(CL["nb"], np.int16),
        cm_cl=np.asarray(CM["cl"], np.int32), cm_idx=np.asarray(CM["idx"], np.int16),
        ev_num=np.asarray(ev_num, np.int64), trig=np.asarray(trig, np.uint32),
        n_sel=n_evt, n_phys=n_phys, n_read=n_read, n_fb_rows=cnt["n_fb_rows"],
        n_odd=cnt["n_odd"], run=run, split=split_of(path), group=RUN_GROUP[run],
        module_names=names, module_is_pwo4=np.asarray(pwo4, dtype=bool),
        module_x=np.array([m.x for m in hy], np.float32), module_y=np.array([m.y for m in hy], np.float32),
        module_sx=np.array([m.size_x for m in hy], np.float32),
        module_sy=np.array([m.size_y for m in hy], np.float32),
        min_cluster_energy=float(ccfg.min_cluster_energy),
        templates=tmpl, calib_file=calib, elapsed_s=elapsed, script_version=SCRIPT_VERSION)
    os.makedirs(out_dir, exist_ok=True)
    tmp = out_path + ".tmp"
    with open(tmp, "wb") as fh:
        np.savez(fh, **out)
    os.replace(tmp, out_path)
    print(f"[done] run {run} split {split_of(path):05d}: {n_read} read, {n_phys} physics, "
          f"{n_evt} selected, {len(out['pu_ev'])} pulses, {len(out['cl_ev'])} clusters, "
          f"{elapsed:.1f} s -> {out_path}", flush=True)


def sample_files():
    fs = []
    for r, s in SAMPLE:
        f = os.path.join(EVIO_DIR, f"prad_{r:06d}.evio.{s:05d}")
        if not os.path.exists(f):
            sys.exit(f"[ERROR] missing {f}")
        fs.append(f)
    return fs


def dump_main(a):
    troot = (os.path.join(HERE, TEMPLATES_ROOT) if a.templates_root == TEMPLATES_ROOT
             else os.path.abspath(a.templates_root))
    files = list(a.files)
    if a.sample:
        files += sample_files()
    if not files:
        sys.exit("[ERROR] no input files (give FILE ... or --sample)")
    out = os.path.abspath(a.out)
    todo = [f for f in files if not os.path.exists(out_name(f, out)[1])]
    print(f"[driver] {len(files)} split(s), {len(files) - len(todo)} already done, "
          f"{len(todo)} to run", flush=True)
    if a.parallel <= 1 or len(todo) <= 1:
        for f in todo:
            process_split(f, out, a.max_events, troot)
        return
    log_dir = os.path.join(out, "logs")
    os.makedirs(log_dir, exist_ok=True)
    running, queue, n_fail = [], list(todo), 0
    while queue or running:
        while queue and len(running) < a.parallel:
            f = queue.pop(0)
            _, op = out_name(f, out)
            log = open(os.path.join(log_dir, os.path.basename(op)[:-4] + ".log"), "w")
            cmd = [sys.executable, os.path.abspath(__file__), "dump", f, "--out", out,
                   "--max-events", str(a.max_events), "--templates-root", troot]
            running.append((subprocess.Popen(cmd, stdout=log, stderr=subprocess.STDOUT), f, log))
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
            print(last[-1] if ok and last else f"[FAIL rc={pr.returncode}] {f}", flush=True)
    print(f"[driver] finished, {n_fail} failed", flush=True)
    sys.exit(1 if n_fail else 0)


# ----------------------------------------------------------------------------
# analyse
# ----------------------------------------------------------------------------
DT_EDGES = np.array(list(range(0, 61, 4)) + [100, 1e9], float)          # |dt| bins (ns)
DT_CAP = np.minimum(DT_EDGES, 240.0)      # plotting edges (>100 bin drawn to 240 ns; max observed |dt| ~238)
R_EDGES = np.array([0.0, 0.1, 0.3, 0.5, 1.0, 2.0, 1e9])                 # ratio bins
E_BINS = [(50, 200), (200, 500), (500, 1000), (1000, 2000), (2000, 1e9)]
RSEED_BINS = [(0, 100), (100, 150), (150, 200), (200, 300), (300, 400), (400, 1e9)]


def binom(k, n):
    """fraction and binomial (sqrt(p(1-p)/n)) error."""
    k, n = np.asarray(k, float), np.asarray(n, float)
    p = np.divide(k, n, out=np.full(np.shape(k), np.nan), where=n > 0)
    e = np.sqrt(np.clip(p * (1 - p), 0, None) / np.maximum(n, 1))
    return p, np.where(n > 0, e, np.nan)


def fmt(x, nd=3):
    return "nan" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"


def pct(k, n, nd=2):
    p, e = binom(k, n)
    return f"{100 * float(p):.{nd}f} +- {100 * float(e):.{nd}f}"


def fold(dt, P=RF_P):
    return np.mod(dt + P / 2, P) - P / 2


def load_group(files):
    """Concatenate the splits of one group into flat arrays (event ids made unique)."""
    D = dict(pu_ev=[], pu_idx=[], pu_t=[], pu_A=[], pu_E=[], pu_fb=[],
             cl_ev=[], cl_seed=[], cl_E=[], cl_x=[], cl_y=[], cm_cl=[], cm_idx=[])
    n_sel = n_phys = n_fb_rows = n_odd = 0
    ev_off = cl_off = 0
    meta = None
    splits = []
    for f in files:
        z = np.load(f, allow_pickle=True)
        if meta is None:
            meta = {k: z[k] for k in ("module_names", "module_is_pwo4", "module_x", "module_y",
                                      "module_sx", "module_sy")}
            meta["min_cluster_energy"] = float(z["min_cluster_energy"])
        for k in ("pu_idx", "pu_t", "pu_A", "pu_E", "pu_fb", "cl_seed", "cl_E", "cl_x", "cl_y",
                  "cm_idx"):
            D[k].append(z[k])
        D["pu_ev"].append(z["pu_ev"].astype(np.int64) + ev_off)
        D["cl_ev"].append(z["cl_ev"].astype(np.int64) + ev_off)
        D["cm_cl"].append(z["cm_cl"].astype(np.int64) + cl_off)
        ev_off += int(z["n_sel"])
        cl_off += len(z["cl_ev"])
        n_sel += int(z["n_sel"]); n_phys += int(z["n_phys"])
        n_fb_rows += int(z["n_fb_rows"]); n_odd += int(z["n_odd"])
        splits.append((int(z["run"]), int(z["split"]), int(z["n_sel"]), float(z["elapsed_s"])))
    D = {k: (np.concatenate(v) if v else np.zeros(0)) for k, v in D.items()}
    D.update(n_sel=n_sel, n_phys=n_phys, n_fb_rows=n_fb_rows, n_odd=n_odd, meta=meta, splits=splits)
    return D


def channel_table(D):
    """Per fired MF channel (fb channels excluded): ev, idx, npk, A1, E1, inwin flag, plus the
    extra pulses (dt, ratio, t, A) relative to the in-time pulse 1."""
    fb = D["pu_fb"].astype(bool)
    ev, idx = D["pu_ev"][~fb], D["pu_idx"][~fb].astype(np.int64)
    t, A, E = D["pu_t"][~fb].astype(float), D["pu_A"][~fb].astype(float), D["pu_E"][~fb].astype(float)
    key = ev * 4096 + idx
    o = np.argsort(key, kind="stable")
    ev, idx, t, A, E, key = ev[o], idx[o], t[o], A[o], E[o], key[o]
    uk, first, inv, npk = np.unique(key, return_index=True, return_inverse=True, return_counts=True)
    nch = len(uk)
    inw = (t > PRE_LO) & (t < PRE_HI)
    # pulse 1 = largest-E in-window pulse per channel
    Ew = np.where(inw, E, -np.inf)
    order = np.lexsort((-Ew, inv))
    p1 = order[first]                                   # first (largest Ew) pulse of each channel
    has1 = inw[p1]
    is1 = np.zeros(len(t), bool)
    is1[p1[has1]] = True
    A1 = np.where(has1, A[p1], np.nan)
    E1 = np.where(has1, E[p1], np.nan)
    t1 = np.where(has1, t[p1], np.nan)
    ex = ~is1
    ex_has = has1[inv] & ex
    dt = t[ex_has] - t1[inv][ex_has]
    ratio = A[ex_has] / A1[inv][ex_has]
    ch = dict(ev=ev[first], idx=idx[first], npk=npk, A1=A1, E1=E1, has1=has1, key=uk)
    xs = dict(ch=inv[ex_has], dt=dt, ratio=ratio, A1=A1[inv][ex_has], E1=E1[inv][ex_has])
    tmin = np.minimum.reduceat(t, first)                # earliest pulse of each channel
    sep = t - tmin[inv]
    ch["pair_sep"] = sep[sep > 0]                       # within-channel separation from the earliest pulse
    return ch, xs, nch


def near_flag(ch, xs):
    nch = len(ch["npk"])
    near = np.zeros(nch, bool)
    m = np.abs(xs["dt"]) <= NEAR_NS
    near[np.unique(xs["ch"][m])] = True
    return near


def analyse_group(D):
    ch, xs, nch = channel_table(D)
    near = near_flag(ch, xs)
    npk = ch["npk"]
    pu = npk >= 2
    n_extra = int((npk - 1).sum())
    a100 = ch["has1"] & (ch["A1"] >= 100)
    e100 = ch["has1"] & (ch["E1"] >= 100)
    ev_pu = len(np.unique(ch["ev"][pu]))
    ev_near = len(np.unique(ch["ev"][near]))
    r = dict(n_ev=D["n_sel"], n_phys=D["n_phys"], n_fired=nch, n_pu=int(pu.sum()), n_inwin=int(ch["has1"].sum()),
             n_near=int(near.sum()), n_a100=int(a100.sum()), n_a100_pu=int((a100 & pu).sum()),
             n_e100=int(e100.sum()), n_e100_pu=int((e100 & pu).sum()),
             ev_pu=ev_pu, ev_near=ev_near, n_extra=n_extra, n_fb_rows=D["n_fb_rows"], n_odd=D["n_odd"],
             n_a100_near=int((a100 & near).sum()), n_e100_near=int((e100 & near).sum()),
             n_before=int((xs["dt"] < 0).sum()), n_after=int((xs["dt"] > 0).sum()),
             n_extra_x=len(xs["dt"]))
    # RF comb
    f = np.abs(fold(xs["dt"]))
    fg = np.abs(fold(xs["dt"], GRID_P))
    for lab, m in (("all", np.ones(len(f), bool)), ("ge16", np.abs(xs["dt"]) >= 16)):
        r[f"rf_n_{lab}"] = int(m.sum())
        r[f"rf_on_{lab}"] = int((f[m] < 0.5).sum())
        r[f"grid_on_{lab}"] = int((fg[m] < 0.5).sum())
    # within-channel pulse pairs (separation from the channel's earliest pulse), no pulse-1 needed
    sp = ch["pair_sep"]
    for lo, hi in PAIR_BINS:
        m = (sp >= lo) & (sp < hi)
        r[f"pair_n_{lo}_{hi}"] = int(m.sum())
        r[f"pair_on_{lo}_{hi}"] = int((np.abs(fold(sp[m])) < 0.5).sum())
        r[f"pair_grid_{lo}_{hi}"] = int((np.abs(fold(sp[m], GRID_P)) < 0.5).sum())
    return ch, xs, near, r


def cluster_rates(D, ch, xs, near):
    """Per cluster pile-up flags (any/near) for members and seed, with E and seed r."""
    meta = D["meta"]
    nm = len(meta["module_names"])
    # channel lookup (ev, idx) -> channel row
    key = ch["ev"] * 4096 + ch["idx"]                    # already sorted (ascending)
    cl_ev, cl_seed = D["cl_ev"].astype(np.int64), D["cl_seed"].astype(np.int64)
    ncl = len(cl_ev)
    pu = ch["npk"] >= 2

    def flag(arr, evs, idxs):
        k = evs * 4096 + idxs
        pos = np.searchsorted(key, k)
        pos = np.clip(pos, 0, max(len(key) - 1, 0))
        ok = (key[pos] == k) if len(key) else np.zeros(len(k), bool)
        return np.where(ok, arr[pos], False) if len(key) else np.zeros(len(k), bool)

    out = {}
    for nmf, arr in (("any", pu), ("near", near)):
        s = flag(arr, cl_ev, cl_seed)
        mem_ev = cl_ev[D["cm_cl"]]
        mf = flag(arr, mem_ev, D["cm_idx"].astype(np.int64))
        anym = np.zeros(ncl, bool)
        np.logical_or.at(anym, D["cm_cl"], mf)
        out[f"seed_{nmf}"] = s
        out[f"mem_{nmf}"] = anym | s
    sx, sy = meta["module_x"][cl_seed], meta["module_y"][cl_seed]
    out["E"] = D["cl_E"].astype(float)
    out["rseed"] = np.hypot(sx, sy).astype(float)
    return out


# ---- 2-D table --------------------------------------------------------------
def dt_ratio_table(xs, n_chan_inwin, sel=None):
    s = np.ones(len(xs["dt"]), bool) if sel is None else sel
    H, _, _ = np.histogram2d(np.abs(xs["dt"][s]), xs["ratio"][s], bins=[DT_EDGES, R_EDGES])
    return H


def incomplete(i_dt, j_r):
    return DT_EDGES[i_dt + 1] <= INC_DT or R_EDGES[j_r + 1] <= INC_RATIO


def edge_label(lo, hi, nd=0):
    return f"{lo:g}-{hi:g}" if hi < 1e8 else f">{lo:g}"


# ---- module map -------------------------------------------------------------
def module_map(path, meta, frac, n, label, min_n):
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    pw = meta["module_is_pwo4"]
    fig, ax = plt.subplots(figsize=(8.2, 7.2))
    idx = np.nonzero(pw)[0]
    patches = [Rectangle((meta["module_x"][i] - meta["module_sx"][i] / 2,
                          meta["module_y"][i] - meta["module_sy"][i] / 2),
                         meta["module_sx"][i], meta["module_sy"][i]) for i in idx]
    val = np.where(n[idx] >= min_n, 100 * frac[idx], np.nan)
    pc = PatchCollection(patches, cmap="viridis", edgecolor="none")
    pc.set_array(np.ma.masked_invalid(val))
    ok = np.isfinite(val)
    vmax = np.nanpercentile(val, 99) if ok.any() else 1.0
    pc.set_clim(0, vmax)
    ax.add_collection(pc)
    ax.set_xlim(meta["module_x"][idx].min() - 20, meta["module_x"][idx].max() + 20)
    ax.set_ylim(meta["module_y"][idx].min() - 20, meta["module_y"][idx].max() + 20)
    ax.set_aspect("equal")
    ax.set_xlabel("x (mm)"); ax.set_ylabel("y (mm)")
    cb = fig.colorbar(pc, ax=ax, shrink=0.8)
    cb.set_label("fraction of fired channels with >= 2 MF pulses (%)")
    ax.set_title(f"{label}: PbWO4 pile-up fraction per module (N fired >= {min_n}; "
                 f"grey = fewer)\nMF-incomplete at |dt| < 16 ns / ratio < 0.3: rates are lower limits",
                 fontsize=9)
    ax.set_facecolor("#dddddd")
    fig.tight_layout()
    fig.savefig(path, dpi=130)
    plt.close(fig)


def analyse_main(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    out = a.out
    os.makedirs(out, exist_ok=True)
    files = sorted(glob.glob(os.path.join(a.dump_dir, "*.npz")))
    by_g = {g: [] for g in ORDER}
    for f in files:
        m = re.search(r"(\d{6})_(\d{5})\.npz$", f)
        by_g[RUN_GROUP[int(m.group(1))]].append(f)
    R, CLU, XS, CH, CLD, DAT = {}, {}, {}, {}, {}, {}
    for g in ORDER:
        if not by_g[g]:
            continue
        D = load_group(by_g[g])
        ch, xs, near, r = analyse_group(D)
        r["splits"] = ";".join(f"{ru}/{sp}" for ru, sp, _, _ in D["splits"])
        r["n_splits"] = len(D["splits"])
        R[g], XS[g], CH[g], DAT[g] = r, xs, ch, D
        CLD[g] = cluster_rates(D, ch, xs, near)
        print(f"[analyse] {g}: {len(D['splits'])} splits, {D['n_sel']} events, "
              f"{len(ch['npk'])} fired channels", flush=True)
    groups = list(R)

    # ---- rates_by_type.tsv / rf_comb.tsv
    cols = ["group", "label", "n_splits", "splits", "n_events", "n_fired_ch", "n_inwin_ch", "n_fb_rows",
            "ev_with_pu", "ev_with_pu_frac", "ev_with_pu_err", "ev_with_near", "ev_with_near_frac",
            "ev_with_near_err", "ch_pu", "ch_pu_frac", "ch_pu_err", "ch_near_frac", "ch_near_err",
            "ch_A100", "ch_A100_pu", "ch_A100_pu_frac", "ch_A100_pu_err",
            "ch_E100", "ch_E100_pu", "ch_E100_pu_frac", "ch_E100_pu_err",
            "extra_pulses", "extra_per_event", "extra_before", "extra_after", "extra_no_pulse1"]
    with open(os.path.join(out, "rates_by_type.tsv"), "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for g in groups:
            r = R[g]
            pe, ee = binom(r["ev_pu"], r["n_ev"])
            pn, en = binom(r["ev_near"], r["n_ev"])
            pc_, ec = binom(r["n_pu"], r["n_fired"])
            pcn, ecn = binom(r["n_near"], r["n_inwin"])
            p1, e1 = binom(r["n_a100_pu"], r["n_a100"])
            p2, e2 = binom(r["n_e100_pu"], r["n_e100"])
            n_nop1 = int(r["n_extra"] - r["n_extra_x"])
            row = [g, GROUPS[g][0], r["n_splits"], r["splits"], r["n_ev"], r["n_fired"], r["n_inwin"],
                   r["n_fb_rows"], r["ev_pu"], fmt(pe, 5), fmt(ee, 5), r["ev_near"], fmt(pn, 5), fmt(en, 5),
                   r["n_pu"], fmt(pc_, 5), fmt(ec, 5), fmt(pcn, 5), fmt(ecn, 5),
                   r["n_a100"], r["n_a100_pu"], fmt(p1, 5), fmt(e1, 5),
                   r["n_e100"], r["n_e100_pu"], fmt(p2, 5), fmt(e2, 5),
                   r["n_extra"], fmt(r["n_extra"] / max(r["n_ev"], 1), 4),
                   r["n_before"], r["n_after"], n_nop1]
            fh.write("\t".join(map(str, row)) + "\n")
    with open(os.path.join(out, "rf_comb.tsv"), "w") as fh:
        fh.write("# on_comb = |fold(sep, 4.0077 ns)| < 0.5 ns; on_grid = same for the 4.000 ns FADC sample grid; random = 1.0/period.\n"
                 "# pulse-1 dt rows: extra pulses relative to the in-time pulse (|dt| max ~236 ns, none beyond 260 ns).\n"
                 "# pair rows: every non-earliest MF pulse of a channel vs the channel's earliest pulse, |sep| bin in ns.\n")
        fh.write("group\tsubset\tn\tn_on_comb\ton_comb_frac\ton_comb_err\tn_on_grid\ton_grid_frac\ton_grid_err\trandom_expect\n")
        for g in groups:
            r = R[g]
            rows = [("pulse-1 dt, all |dt|", r["rf_n_all"], r["rf_on_all"], r["grid_on_all"]),
                    ("pulse-1 dt, |dt|>=16 ns", r["rf_n_ge16"], r["rf_on_ge16"], r["grid_on_ge16"])]
            rows += [(f"pair sep {lo}-{hi} ns", r[f"pair_n_{lo}_{hi}"], r[f"pair_on_{lo}_{hi}"],
                      r[f"pair_grid_{lo}_{hi}"]) for lo, hi in PAIR_BINS]
            for tag, n, k, kg in rows:
                p, e = binom(k, n)
                pg, eg = binom(kg, n)
                fh.write(f"{g}\t{tag}\t{n}\t{k}\t{fmt(p, 4)}\t{fmt(e, 4)}\t{kg}\t{fmt(pg, 4)}\t{fmt(eg, 4)}\t{1.0 / RF_P:.4f}\n")

    # ---- cluster rates
    with open(os.path.join(out, "cluster_rates.tsv"), "w") as fh:
        fh.write("group\tbin_kind\tbin\tn_clusters\tany_member_pu\tany_member_frac\tany_member_err\t"
                 "seed_pu\tseed_frac\tseed_err\tnear_member_frac\tnear_seed_frac\n")
        for g in groups:
            c = CLD[g]
            sets = [("all", "all", np.ones(len(c["E"]), bool))]
            sets += [("E_MeV", edge_label(lo, hi), (c["E"] >= lo) & (c["E"] < hi)) for lo, hi in E_BINS]
            sets += [("seed_r_mm", edge_label(lo, hi), (c["rseed"] >= lo) & (c["rseed"] < hi))
                     for lo, hi in RSEED_BINS]
            for kind, lab, m in sets:
                n = int(m.sum())
                km, ks = int(c["mem_any"][m].sum()), int(c["seed_any"][m].sum())
                pm, em = binom(km, n)
                ps, es = binom(ks, n)
                pmn, _ = binom(int(c["mem_near"][m].sum()), n)
                psn, _ = binom(int(c["seed_near"][m].sum()), n)
                fh.write(f"{g}\t{kind}\t{lab}\t{n}\t{km}\t{fmt(pm, 5)}\t{fmt(em, 5)}\t{ks}\t{fmt(ps, 5)}\t"
                         f"{fmt(es, 5)}\t{fmt(pmn, 5)}\t{fmt(psn, 5)}\n")

    # ---- modules
    MODR = {}
    for g in groups:
        ch, meta = CH[g], DAT[g]["meta"]
        nm = len(meta["module_names"])
        nf = np.bincount(ch["idx"], minlength=nm)
        npu = np.bincount(ch["idx"][ch["npk"] >= 2], minlength=nm)
        frac, err = binom(npu, nf)
        MODR[g] = (nf, npu, frac, err)
        with open(os.path.join(out, f"module_rates_{g}.tsv"), "w") as fh:
            fh.write("module\tn_fired\tn_ge2\tfrac\terr\n")
            for i in np.nonzero(meta["module_is_pwo4"] & (nf > 0))[0]:
                fh.write(f"{meta['module_names'][i]}\t{nf[i]}\t{npu[i]}\t{fmt(frac[i], 5)}\t{fmt(err[i], 5)}\n")
        module_map(os.path.join(out, f"module_map_{g}.png"), meta, np.nan_to_num(frac), nf,
                   GROUPS[g][0], a.map_min_fired)

    # ---- dt x ratio
    TAB = {}
    for g in groups:
        xs, r = XS[g], R[g]
        N1 = max(r["n_inwin"], 1)
        N1a = max(r["n_a100"], 1)
        H = dt_ratio_table(xs, N1)
        Ha = dt_ratio_table(xs, N1, sel=xs["A1"] >= 100)
        TAB[g] = (H, Ha, N1, N1a)
        with open(os.path.join(out, f"dt_ratio_{g}.tsv"), "w") as fh:
            fh.write(f"# group {g}; pulse 1 = in-time pulse (largest-E in 100-200 ns); normalised to the "
                     f"{r['n_inwin']} fired channels with a pulse 1 ({r['n_a100']} with A1>=100); "
                     "'*' = MF-incomplete cell (|dt|<16 ns or ratio<0.3)\n")
            fh.write("dt_lo\tdt_hi\tratio_lo\tratio_hi\tn\tprob\terr\tn_A1ge100\tprob_A1ge100\tflag\n")
            for i in range(len(DT_EDGES) - 1):
                for j in range(len(R_EDGES) - 1):
                    k = H[i, j]
                    p, e = binom(k, N1)
                    p2, _ = binom(Ha[i, j], N1a)
                    fh.write(f"{DT_EDGES[i]:g}\t{DT_EDGES[i + 1]:g}\t{R_EDGES[j]:g}\t{R_EDGES[j + 1]:g}\t"
                             f"{int(k)}\t{p:.3e}\t{e:.1e}\t{int(Ha[i, j])}\t{p2:.3e}\t"
                             f"{'*' if incomplete(i, j) else ''}\n")

    # ---- figure: marginals + 2-D per group
    cmap = plt.get_cmap("tab10")
    with PdfPages(os.path.join(out, "dt_ratio_distributions.pdf")) as pdf:
        fig, axs = plt.subplots(1, 2, figsize=(11, 4.4))
        for k, g in enumerate(groups):
            H, _, N1, _ = TAB[g]
            ls = "--" if g.startswith("x17") else "-"
            dtp = H.sum(1) / N1 / np.diff(DT_CAP)                # per ns, all bins (last bin drawn to DT_CAP)
            dtp = np.where(dtp > 0, dtp, np.nan)
            axs[0].stairs(dtp, DT_CAP, color=cmap(k), ls=ls, label=GROUPS[g][0])
            rp = H.sum(0) / N1
            rp = np.where(rp > 0, rp, np.nan)
            axs[1].step(range(len(rp)), rp, where="mid", color=cmap(k), ls=ls, label=GROUPS[g][0])
        axs[0].axvspan(0, INC_DT, color="grey", alpha=0.25)
        axs[0].text(INC_DT / 2, 0.97, "MF-incomplete", rotation=90, va="top", ha="center", fontsize=7,
                    transform=axs[0].get_xaxis_transform())
        axs[0].set_yscale("log")
        axs[0].set_xlabel("|dt| = |t_k - t_1| (ns); bins 4 ns to 60, 60-100, >100 (drawn to 240 ns, max observed ~238)")
        axs[0].set_ylabel("extra pulses per fired in-time channel per ns")
        axs[0].set_title("dt marginal (all bins; per ns)", fontsize=9)
        axs[0].legend(fontsize=7)
        rl = [edge_label(R_EDGES[j], R_EDGES[j + 1]) for j in range(len(R_EDGES) - 1)]
        axs[1].set_xticks(range(len(rl))); axs[1].set_xticklabels(rl)
        axs[1].axvspan(-0.5, 1.5, color="grey", alpha=0.25)
        axs[1].text(0.5, 0.97, "MF-incomplete", ha="center", va="top", fontsize=7,
                    transform=axs[1].get_xaxis_transform())
        axs[1].set_yscale("log")
        axs[1].set_xlabel("ratio A_k / A_1"); axs[1].set_ylabel("extra pulses per fired in-time channel")
        axs[1].set_title("ratio marginal (all |dt|)", fontsize=9)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        for g in groups:
            H, Ha, N1, N1a = TAB[g]
            P = H / N1
            fig, ax = plt.subplots(figsize=(9, 5.2))
            im = ax.imshow(np.where(P > 0, P, np.nan), origin="lower", aspect="auto",
                           norm=matplotlib.colors.LogNorm(vmin=max(np.nanmin(P[P > 0]), 1e-7) if (P > 0).any() else 1e-6,
                                                          vmax=max(P.max(), 1e-6)), cmap="viridis")
            for i in range(P.shape[0]):
                for j in range(P.shape[1]):
                    if incomplete(i, j):
                        ax.add_patch(matplotlib.patches.Rectangle((j - .5, i - .5), 1, 1, fill=False,
                                                                  hatch="///", edgecolor="grey", lw=0))
            ax.set_xticks(range(len(rl))); ax.set_xticklabels(rl)
            dl = [edge_label(DT_EDGES[i], DT_EDGES[i + 1]) for i in range(len(DT_EDGES) - 1)]
            ax.set_yticks(range(len(dl))); ax.set_yticklabels(dl, fontsize=7)
            ax.set_xlabel("ratio A_k / A_1"); ax.set_ylabel("|dt| (ns)")
            ax.set_title(f"{GROUPS[g][0]}: extra pulses per fired in-time channel (N = {N1}); "
                         "hatched = MF-incomplete", fontsize=9)
            fig.colorbar(im, ax=ax).set_label("probability per channel")
            fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

    write_report(a, R, CLD, TAB, MODR, DAT)


def tab3(H, N1):
    """3x3: |dt| (<16, 16-60, >60) x ratio (<0.3, 0.3-1, >1)."""
    d_edges = [0, 16, 60, 1e9]
    r_edges = [0, 0.3, 1.0, 1e9]
    # fold fine bins: DT bins edges are multiples of 4 up to 60, then 100, inf
    di = [(0, 4), (4, 15), (15, len(DT_EDGES) - 1)]     # bins [0,16), [16,60), [60,inf)
    rj = [(0, 2), (2, 4), (4, 6)]                       # ratio [0,0.3), [0.3,1), [1,inf)
    out = np.zeros((3, 3))
    for a_, (i0, i1) in enumerate(di):
        for b_, (j0, j1) in enumerate(rj):
            out[a_, b_] = H[i0:i1, j0:j1].sum() / N1
    return out


def write_report(a, R, CLD, TAB, MODR, DAT):
    out = a.out
    L = []
    L.append("# Real-data MF pile-up rates, PbWO4 (2026-10-27)")
    L.append("")
    L.append("Best-50 spline bank, MF nsigma 5 / min_amp 5 / kmax 4 / r_model 0.3. Event = physics trigger "
             "(0xf00, no LMS/alpha). Pulse 1 = in-time pulse (largest-E, 100<t<200 ns); extra pulse = any "
             "other MF pulse in the readout window. `*` = MF-incomplete (|dt|<16 ns or ratio<0.3).")
    L.append("")
    L.append("Table 1 (percent +- binomial): any = >= 2 MF pulses in the readout window; near = an extra pulse within 60 ns of pulse 1. "
             "Denominators: ev = selected events; ch any = fired MF channels (n_fired, >= 1 pulse); ch near = channels with a pulse 1 "
             "(n_inwin); A1/E1 columns = fired channels with a pulse 1 and A1 >= 100 ADC / E1 >= 100 MeV; clusters = E >= 50 MeV clusters.")
    L.append("")
    L.append("| group | splits | events | ev any (%) | ev near (%) | ch any (%) | ch near (%) | ch any, A1>=100 (%) | ch any, E1>=100 MeV (%) | extra/event | clu member any / near (%) | clu seed any / near (%) |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
    for g, r in R.items():
        c = CLD[g]
        n = len(c["E"])
        pp = lambda k, m, nd=3: f"{100 * k / max(m, 1):.{nd}f}"
        L.append(f"| {GROUPS[g][0]} | {r['n_splits']} | {r['n_ev']} | {pct(r['ev_pu'], r['n_ev'], 2)} | "
                 f"{pct(r['ev_near'], r['n_ev'], 2)} | {pct(r['n_pu'], r['n_fired'], 3)} | {pct(r['n_near'], r['n_inwin'], 3)} | "
                 f"{pct(r['n_a100_pu'], r['n_a100'], 3)} | {pct(r['n_e100_pu'], r['n_e100'], 3)} | "
                 f"{r['n_extra'] / max(r['n_ev'], 1):.3f} | {pp(c['mem_any'].sum(), n, 2)} / {pp(c['mem_near'].sum(), n, 2)} | "
                 f"{pp(c['seed_any'].sum(), n, 2)} / {pp(c['seed_near'].sum(), n, 2)} |")
    L.append("")
    L.append("Table 2: RF comb (extra pulses with |dt mod 4.0077| < 0.5 ns; random = 25 %) and the 3x3 per-channel probability (%) of an extra pulse, "
             "|dt| (<16 / 16-60 / >60 ns) x ratio A_k/A_1 (<0.3 / 0.3-1 / >1), normalised to the channels with a pulse 1 (n_inwin; the 3x3 sums extra pulses, a channel can count in several cells); "
             "`*` = MF-incomplete. Full 4 ns x 6 ratio bins in `dt_ratio_<g>.tsv`.")
    L.append("")
    L.append("| group | on-comb all (%) | on-comb &#124;dt&#124;>=16 (%) | extra before / after pulse 1 | <16 x <0.3* | <16 x 0.3-1* | <16 x >1* | 16-60 x <0.3* | 16-60 x 0.3-1 | 16-60 x >1 | >60 x <0.3* | >60 x 0.3-1 | >60 x >1 |")
    L.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for g, r in R.items():
        t = 100 * tab3(TAB[g][0], TAB[g][2])
        L.append(f"| {GROUPS[g][0]} | {pct(r['rf_on_all'], r['rf_n_all'], 1)} | "
                 f"{pct(r['rf_on_ge16'], r['rf_n_ge16'], 1)} | {r['n_before']} / {r['n_after']} | "
                 + " | ".join(f"{v:.3f}" for v in t.ravel()) + " |")
    L.append("")
    L.append("Top modules by pile-up fraction (N fired >= %d), per group: see `module_rates_<g>.tsv`, maps `module_map_<g>.png`." % a.min_fired)
    for g in R:
        nf, npu, fr, _ = MODR[g]
        pw = DAT[g]["meta"]["module_is_pwo4"]
        ok = np.nonzero(pw & (nf >= a.min_fired))[0]
        top = ok[np.argsort(-fr[ok])][:5]
        L.append(f"- {GROUPS[g][0]}: " + ", ".join(f"{DAT[g]['meta']['module_names'][i]} {100 * fr[i]:.1f}% ({nf[i]})" for i in top))
    L.append("")
    c2 = CLD["2p24"] if "2p24" in CLD else None
    if c2 is not None:
        def fr(m, k):
            return 100 * c2[k][m].sum() / max(m.sum(), 1)
        eb = " / ".join(f"{fr((c2['E'] >= lo) & (c2['E'] < hi), 'mem_any'):.1f}" for lo, hi in E_BINS)
        rb = " / ".join(f"{fr((c2['rseed'] >= lo) & (c2['rseed'] < hi), 'mem_any'):.1f}" for lo, hi in RSEED_BINS)
        L.append("Cluster dependence, 2.24 H2, any-member pile-up (%%): E 50-200 / 200-500 / 500-1000 / 1000-2000 / >2000 MeV = %s; "
                 "seed r 0-100 / 100-150 / 150-200 / 200-300 / 300-400 / >400 mm = %s (outer edge >400 mm is the highest bin, few clusters). "
                 "All bins in `cluster_rates.tsv`." % (eb, rb))
        L.append("")
    fb = ", ".join(f"{g} {r['n_fb_rows']}" for g, r in R.items())
    rng = lambda key, nkey, sel=None: (min(100 * R[g][key] / max(R[g][nkey], 1) for g in R if sel is None or sel(g)),
                                       max(100 * R[g][key] / max(R[g][nkey], 1) for g in R if sel is None or sel(g)))
    on_lo, on_hi = rng("rf_on_all", "rf_n_all")
    gr_lo, gr_hi = rng("grid_on_all", "rf_n_all")
    ta = lambda g: g.startswith("x17")
    p_on = rng("pair_on_260_400", "pair_n_260_400", ta)
    p_gr = rng("pair_grid_260_400", "pair_n_260_400", ta)
    nn = [R[g]["pair_n_260_400"] for g in R if not ta(g)]
    nh = (min(nn), max(nn))
    L.append("Caveats: (1) MF secondary search (mf_pulse.py ~1054-1091: 12 ns MIN_SEP, 0.3*A_parent gate only after the "
             "parent onset, argmax-then-stop) misses |dt|<16-20 ns and ratio<0.3 (output/templates_wassim_20261027/why_hard_cells.md): all "
             "rates are lower limits for true pile-up, and the 'any' rate includes out-of-time background pulses far from pulse 1. "
             "(2) Event caps: first %d physics events per split, early/mid/late splits (see `rates_by_type.tsv:splits`); errors binomial, "
             "channels within an event are correlated. (3) WaveAnalyzer-fallback rows excluded (n_fb_rows per group: %s). (4) Ta runs not pooled; "
             "H2 types pooled over their runs. (5) Cluster pile-up flags use MF channels only; lead-glass members have none; map colouring needs "
             "N fired >= %d, top-5 lists N fired >= %d. "
             "(6) RF comb: pulse-1 dt never exceeds ~236 ns. For |dt| < ~260 ns the 4.000 ns FADC sample grid lies within 0.5 ns of the 4.0077 ns comb, "
             "so the on-comb fraction of pulse-1 dt (%.0f-%.0f %% over the groups, random 25 %%) equals the 4.000 ns grid fraction (%.0f-%.0f %%) "
             "and is NOT separable from grid locking. The grid-decorrelated test is the within-channel pulse pairs (separation from the earliest pulse) at "
             "|sep| 260-400 ns (`rf_comb.tsv`): Ta on-comb %.0f-%.0f %% vs grid %.0f-%.0f %% (random 25 %%), i.e. a real but smaller RF excess than the "
             "pulse-1 numbers suggest; H2/C have only %d-%d such pairs per group, too few for the test."
             % (a.max_events_note, fb, a.map_min_fired, a.min_fired, on_lo, on_hi, gr_lo, gr_hi,
                p_on[0], p_on[1], p_gr[0], p_gr[1], nh[0], nh[1]))
    L.append("")
    L.append("Command: `python3 pileup_rates.py dump --sample --parallel 8 --out output/pileup_rates_20261027/dump --max-events %d`; "
             "`python3 pileup_rates.py analyse --dump-dir output/pileup_rates_20261027/dump --out output/pileup_rates_20261027`." % a.max_events_note)
    with open(os.path.join(out, "report.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    d = sub.add_parser("dump", help="EVIO split(s) -> npz")
    d.add_argument("files", nargs="*", help="EVIO split files (each processed alone)")
    d.add_argument("--out", required=True, metavar="DIR")
    d.add_argument("--sample", action="store_true", help="the built-in SAMPLE split list")
    d.add_argument("--max-events", type=int, default=0,
                   help="stop after N physics-type events per split (0 = all)")
    d.add_argument("--parallel", type=int, default=1, metavar="N")
    d.add_argument("--templates-root", default=TEMPLATES_ROOT, metavar="DIR")
    n = sub.add_parser("analyse", help="npz -> tables, maps, report")
    n.add_argument("--dump-dir", required=True)
    n.add_argument("--out", required=True)
    n.add_argument("--min-fired", type=int, default=2000,
                   help="min fired events for a module to enter the top-5 lists in report.md")
    n.add_argument("--map-min-fired", type=int, default=200,
                   help="min fired events for a module to be coloured in module_map_<g>.png")
    n.add_argument("--max-events-note", type=int, default=60000,
                   help="event cap used in the dump (report text only)")
    a = ap.parse_args()
    if a.cmd == "dump":
        dump_main(a)
    else:
        analyse_main(a)


if __name__ == "__main__":
    main()
