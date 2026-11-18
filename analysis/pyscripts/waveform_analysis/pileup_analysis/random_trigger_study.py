#!/usr/bin/env python3
"""random_trigger_study.py -- what the local run 24578 (random trigger, 0.7 GeV beam; elog 4498531) tells about pedestal, noise, occupancy and pulse shape.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/random_trigger_study.py scan <out.npz> <run> <evio file> [<evio file> ...] [--max-events N] [--ref-from <scan.npz>]
    python3 waveform_analysis/pileup_analysis/random_trigger_study.py report <scan dir> <out stem> [--ref-sigma <run_24590.npz>]

scan: one EVIO file per call (as many as given, each with its own noise reference).  Event classes by the trigger word: bits == 0 (class 0, the random
trigger of run 24578), bits == 0x100 (class 1, the production trigger; run 24590 has only this class), anything else only counted.  Events read with
select_event/info/fadc (no GEM decoding).  Waveforms: HyCal channels, dead modules left out (as pileup_rates.scan), blocks of BLOCK events.
The DEFAULT pileup_rates fit (MONO_EDGE, NOISE_FLAG off, SHAPE_SCALE 1: forced below), noise reference per file from the first N_REF_EVENTS events of the
classes 0 and 1 (pileup_rates.noise_reference; same recipe as the production scan; with --ref-from the fit takes sigma_ref and b_ref of another run's scan instead: the random trigger gives < 20 waveforms per module in 20000 events), bank of the run type 0p73GeV_H2 (run 24578 is mapped to it in memory only;
pileup_rates.py is not edited).  Fit only if the raw peak (max - median of the first 30 samples) >= PRECUT_ADC and the waveform is not bad.
Per waveform: baseline statistics of the first 30 and of all 100 samples, fit k, A[6], onset t[6], chi2/dof.  Sums of the residual profile of isolated PbWO4
pulses (k == 1, A 100-2000 ADC, no overflow) per class, amplitude bin and tau bin (T1 recipe, without the phase split).  Reservoirs of example waveforms.
report: PDF + TSV (see the report docstring text on the first page).  Cache under output/random_trigger_24578/ (not for git).

PbWO4 (W*) only everywhere, as the standing rule says: lead glass is read and fitted but dropped in the report.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
for _v, _d in (("PILEUP_MONO_EDGE", "0"), ("PILEUP_NOISE_FLAG", "0"), ("PILEUP_SHAPE_SCALE", "1")):
    os.environ[_v] = _d                                  # default fit behaviour only

import sys
import re
import time
import glob
import argparse
import collections
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "exotic_events"))
import numpy as np                                      # noqa: E402

RUN_TYPE_EXTRA = {24578: "0p73GeV_H2"}                  # in memory only: the 0.7 GeV bank for the random-trigger run
NS, DT = 100, 4.0
TT = DT * np.arange(NS)
KMAX = 6
CLS_BITS = {0: 0x0, 1: 0x100}
N_REF_EVENTS = 20000
WIDE_NCH = 300                                          # events with this many HyCal channels are array-wide (vetoed as in pileup_rates.VETO_NKEPT)
ISO_T = (20.0, 340.0)                                   # ns, onset range of the isolated pulses of the residual profile (fit window complete)
BLOCK = 256
AMP_EDGES = np.array([100.0, 200.0, 400.0, 800.0, 2000.0])
NAB = len(AMP_EDGES) - 1
TAU = np.arange(-20.0, 200.1, 4.0)
NT = TAU.size
EX_AMP_EDGES = np.array([10.0, 20.0, 50.0, 200.0, 1000.0, np.inf])     # example bins by fitted amplitude of the k == 1 pulse
EX_PER = dict(k1=6, k0=25, k2=25)                       # reservoir size per file (k1 per amplitude bin)
SEED = 20261117


def clock_s(ts):
    """Event timestamp (4 ns ticks of the 250 MHz TI clock) to seconds."""
    return np.asarray(ts, dtype=np.float64) * 4e-9


# ===================================================================== scan
def scan(out, evios, run, max_events, ref_from=None):
    import event_display as ed
    from prad2py import dec
    pr = ed.pr
    assert not pr.MONO_EDGE and not pr.NOISE_FLAG and pr.SHAPE_SCALE == 1.0
    for r, t in RUN_TYPE_EXTRA.items():
        pr.RUN_TYPE_OF.setdefault(r, t)
    t_start = time.monotonic()
    ctx = ed.Ctx(run)
    nm = len(ctx.names)
    bank_row = ctx.bank.row
    rng = np.random.default_rng(SEED + run + int(re.search(r"\.(\d+)$", evios[0]).group(1)))
    print(f"[setup] run {run} {ctx.rt}, {nm} modules, {len(evios)} files, max events {max_events or 'all'}", flush=True)

    fields = dict(fi=np.int8, cls=np.int8, ev=np.int32, mi=np.int16, ped=np.float32, m30=np.float32, s30=np.float32, mad30=np.float32, fm=np.float32, fs=np.float32,
                  raw=np.float32, ovf=np.bool_, bad=np.bool_, k=np.int8, a=np.float32, t=np.float32, chi2=np.float32, dped=np.float32, trunc=np.bool_)
    rec = {k: [] for k in fields}
    ev_ts, ev_cls, ev_nch, ev_num, ev_fi, ev_nwf = [], [], [], [], [], []
    bits_count = collections.Counter()
    G_n = np.zeros((2, NAB, NT)); G_sx = np.zeros_like(G_n); G_sxx = np.zeros_like(G_n); G_sv = np.zeros_like(G_n)
    ex = {}
    ref_n, ref_b, ref_sig = [], [], []

    def reservoir(name, new, per):
        old = ex.get(name)
        if old is not None:
            new = {k: np.concatenate([old[k], new[k]]) for k in new}
        o = np.argsort(new["key"])[:per]
        ex[name] = {k: v[o] for k, v in new.items()}

    def events(ch, limit):
        """Yield (cls, bits, info, event data or None) of the physics events of an open channel."""
        n = 0
        while ch.read() == dec.Status.success:
            if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
                continue
            for i in range(ch.get_n_events()):
                ch.select_event(i)
                inf = ch.info()
                bits = int(inf.trigger_bits)
                cls = 0 if bits == CLS_BITS[0] else 1 if bits == CLS_BITS[1] else -1
                n += 1
                yield cls, bits, inf, (ch.fadc() if cls >= 0 else None)
                if limit and n >= limit:
                    return

    def reference(f):
        X, M, rows, mis = [], [], [], []
        ch = ctx.open(f)
        n = 0
        for cls, bits, inf, fe in events(ch, N_REF_EVENTS):
            if cls < 0:
                continue
            for mi, smp in ctx.channels(fe):
                rows.append(smp)
                mis.append(mi)
            if len(rows) > 20000:
                Y, mm = np.stack(rows), np.array(mis)
                rows.clear(), mis.clear()
                good = ~pr.bad_waveform(Y)
                X.append(Y[good, :pr.PED_NSAMPLES]); M.append(mm[good])
        if rows:
            Y, mm = np.stack(rows), np.array(mis)
            good = ~pr.bad_waveform(Y)
            X.append(Y[good, :pr.PED_NSAMPLES]); M.append(mm[good])
        ch.close()
        return pr.noise_reference(np.concatenate(X), np.concatenate(M), nm)

    ext = None
    if ref_from:                                         # noise reference of another run (the random trigger leaves most modules without enough waveforms in 20000 events)
        z = np.load(ref_from, allow_pickle=False)
        okr = z["mod_n_ref"] >= pr.MIN_REF_WAVEFORMS
        with np.errstate(all="ignore"):
            sg = np.nanmedian(np.where(okr, z["mod_sig_ref"], np.nan), axis=0)
            bb = np.nanmedian(np.where(okr, z["mod_b_ref"], np.nan), axis=0)
        sg = np.where(np.isfinite(sg), sg, np.nanmedian(sg)).astype(np.float32)
        ext = (bb.astype(np.float32), sg)
        print(f"[reference] sigma_ref and b_ref taken from {Path(ref_from).name}: sigma_ref median {np.median(sg):.3f} ADC", flush=True)
    cur = dict(fi=0, sig_ref=None, b_ref=None)
    rows, evs, mis = [], [], []

    def flush():
        if not rows:
            return
        Y, ev, mi = np.stack(rows), np.array(evs), np.array(mis)
        rows.clear(), evs.clear(), mis.clear()
        n = len(mi)
        sig_ref, b_ref = cur["sig_ref"][mi], cur["b_ref"][mi]
        bad = pr.bad_waveform(Y)
        med30 = np.median(Y[:, :pr.PED_NSAMPLES], axis=1)
        sel = np.nonzero(~bad & (Y.max(axis=1) - med30 >= pr.PRECUT_ADC))[0]
        k = np.zeros(n, np.int8)
        A = np.zeros((n, KMAX), np.float32)
        T = np.full((n, KMAX), np.nan, np.float32)
        chi2 = np.full(n, np.nan, np.float32)
        dped = np.zeros(n, np.float32)
        trunc = np.zeros(n, bool)
        f = None
        if sel.size:
            f = pr.fit(ctx.mf, Y[sel], bank_row[mi[sel]], sig_ref[sel], b_ref[sel])
            f.pop("far")
            k[sel], A[sel], T[sel], chi2[sel], dped[sel], trunc[sel] = f["k"], f["a"], f["t"], f["chi2"], f["dped"], f["trunc_main"]
        x = Y[:, :pr.PED_NSAMPLES].astype(np.float64)
        ped = pr.pedestal(Y)
        ovf = (Y >= pr.OVERFLOW).any(axis=1)
        cls = np.array(ev_cls)[ev]
        vals = dict(fi=np.full(n, cur["fi"]), cls=cls, ev=ev, mi=mi, ped=ped, m30=x.mean(axis=1), s30=x.std(axis=1),
                    mad30=1.4826 * np.median(np.abs(x - np.median(x, axis=1)[:, None]), axis=1), fm=Y.mean(axis=1), fs=Y.std(axis=1),
                    raw=Y.max(axis=1) - ped, ovf=ovf, bad=bad, k=k, a=A, t=T, chi2=chi2, dped=dped, trunc=trunc)
        for key, v in vals.items():
            rec[key].append(np.asarray(v).astype(fields[key]))
        if f is None:
            return
        # isolated PbWO4 pulses: residual profile sums
        wide = np.array(ev_nch)[ev] >= WIDE_NCH                 # array-wide events (LED / pulser flashes): left out of the sums and the examples
        pw = ctx.pwo4[mi[sel]] & ~wide[sel]
        a1 = f["a"][:, 0]
        t1 = f["t"][:, 0]
        iso = (f["k"] == 1) & pw & ~ovf[sel] & (a1 >= AMP_EDGES[0]) & (a1 < AMP_EDGES[-1]) & (t1 >= ISO_T[0]) & (t1 <= ISO_T[1])
        idx = np.nonzero(iso)[0]
        if idx.size:
            valid = np.arange(KMAX)[None, :] < f["k"][idx][:, None]
            t0 = np.where(valid, np.nan_to_num(f["t"][idx]), 0.0)
            a = f["a"][idx]
            noise = np.sqrt(f["sig"][idx].astype(np.float64) ** 2 + 1.0 / 12.0)
            var = ctx.mf.variance(bank_row[mi[sel][idx]], t0.astype(np.float32), a, valid, noise, NS)
            res = Y[sel[idx]] - f["ped"][idx][:, None] - f["comp"][idx].sum(axis=1)
            A0 = a[:, 0].astype(np.float64)[:, None]
            xr, vr = res / A0, var / A0 ** 2
            ab = np.digitize(a[:, 0], AMP_EDGES) - 1
            gc = cls[sel[idx]]
            ok_c = gc >= 0
            tau = TT[None, :] - t0[:, :1]
            ib = np.floor((tau + 22.0) / 4.0).astype(np.int64)
            ok = (ib >= 0) & (ib < NT) & ok_c[:, None]
            g = (((gc * NAB + ab)[:, None]) * NT + ib)[ok]
            for arr, w in ((G_n, None), (G_sx, xr[ok]), (G_sxx, xr[ok] ** 2), (G_sv, vr[ok])):
                arr += np.bincount(g, weights=w, minlength=arr.size).reshape(arr.shape)
        # example reservoirs (class 0 PbWO4)
        s0 = pw & (cls[sel] == 0)
        am = f["a"].max(axis=1)
        picks = [("k0", np.nonzero(~bad & ~wide & (cls == 0) & ctx.pwo4[mi] & (k == 0))[0], None)]
        picks += [("k2", sel[np.nonzero(s0 & (f["k"] >= 2))[0]], None)]
        ab_ex = np.digitize(am, EX_AMP_EDGES) - 1
        picks += [(f"k1_{b}", sel[np.nonzero(s0 & (f["k"] == 1) & (ab_ex == b))[0]], None) for b in range(len(EX_AMP_EDGES) - 1)]
        pos_in_sel = {int(s): j for j, s in enumerate(sel)}
        for name, rows_, _ in picks:
            if rows_.size == 0:
                continue
            key = rng.random(rows_.size)
            keep = np.argsort(key)[:EX_PER["k0" if name == "k0" else "k2" if name == "k2" else "k1"]]
            rows_, key = rows_[keep], key[keep]
            new = dict(key=key, Y=Y[rows_], mi=mi[rows_], evn=np.array(ev_num)[ev[rows_]], fi=np.full(rows_.size, cur["fi"]), k=k[rows_], t=T[rows_], a=A[rows_], chi2=chi2[rows_],
                       ped=ped[rows_], sig=sig_ref[rows_].astype(np.float32), comp=np.zeros((rows_.size, KMAX, NS), np.float32))
            if name != "k0":
                jj = np.array([pos_in_sel[int(r)] for r in rows_])
                new["comp"] = f["comp"][jj].astype(np.float32)
                new["ped"] = f["ped"][jj].astype(np.float32)
            reservoir(name, new, EX_PER["k0" if name == "k0" else "k2" if name == "k2" else "k1"])

    n_seen = 0
    for fi, evio in enumerate(evios):
        flush()
        t_ref = time.monotonic()
        n_ref, b_ref, sig_ref = reference(evio)
        ref_n.append(n_ref), ref_b.append(b_ref), ref_sig.append(sig_ref)
        cur.update(fi=fi, sig_ref=sig_ref, b_ref=b_ref)
        if ext is not None:
            cur.update(b_ref=ext[0], sig_ref=ext[1])
        print(f"[reference] {Path(evio).name}: sigma_ref median {np.nanmedian(sig_ref[n_ref >= pr.MIN_REF_WAVEFORMS]):.3f} ADC, {int((n_ref < pr.MIN_REF_WAVEFORMS).sum())} modules below {pr.MIN_REF_WAVEFORMS} waveforms; {time.monotonic() - t_ref:.0f} s", flush=True)
        ch = ctx.open(evio)
        for cls, bits, inf, fe in events(ch, max_events):
            n_seen += 1
            bits_count[(fi, bits)] += 1
            if cls < 0:
                continue
            ev = len(ev_num)
            chans = ctx.channels(fe)
            ev_num.append(int(inf.event_number)); ev_ts.append(int(inf.timestamp)); ev_cls.append(cls); ev_nch.append(len(chans)); ev_fi.append(fi)
            for mi, smp in chans:
                rows.append(smp); evs.append(ev); mis.append(mi)
            if len(ev_num) % BLOCK == 0:
                flush()
            if n_seen % 20000 == 0:
                print(f"  [progress] {n_seen} events, {time.monotonic() - t_start:.0f} s", flush=True)
        ch.close()
    flush()
    r = {f"w_{k}": (np.concatenate(v) if v else np.zeros((0, KMAX) if k in ("a", "t") else 0, fields[k])) for k, v in rec.items()}
    exo = {f"ex_{n}_{k}": v for n, d in ex.items() for k, v in d.items()}
    bk = sorted(bits_count)
    os.makedirs(os.path.dirname(os.path.abspath(out)), exist_ok=True)
    np.savez_compressed(out, run=run, files=np.array([Path(e).name for e in evios]), n_seen=n_seen, wall_s=time.monotonic() - t_start,
                        ev_ts=np.array(ev_ts, np.int64), ev_cls=np.array(ev_cls, np.int8), ev_nch=np.array(ev_nch, np.int32), ev_num=np.array(ev_num, np.int64), ev_fi=np.array(ev_fi, np.int8),
                        bits_fi=np.array([b[0] for b in bk]), bits_val=np.array([b[1] for b in bk]), bits_n=np.array([bits_count[b] for b in bk]),
                        mod_name=ctx.names, mod_pwo4=ctx.pwo4, mod_n_ref=np.array(ref_n), mod_b_ref=np.array(ref_b), mod_sig_ref=np.array(ref_sig),
                        G_n=G_n, G_sx=G_sx, G_sxx=G_sxx, G_sv=G_sv, **r, **exo)
    print(f"[done] {out}: {n_seen} events, {len(ev_num)} in classes 0/1, {r['w_mi'].size} waveforms; {time.monotonic() - t_start:.0f} s", flush=True)


# ===================================================================== report
FID = (160.0, 360.0)          # ns, onset range of the rate: the late plateau of the onset density of the random pulses (the report shows the depleted windows below it)
FID_ANY_NS = 340.0            # ns, record length with a complete fit (onsets 20-360 ns): the time of the 'anywhere in the record' prediction
PLAT2 = (56.0, 120.0)         # ns, second plateau, before the depleted in-time window
DEPL = ((0.0, 48.0), (120.0, 160.0))   # ns, windows where the onset density is depleted
WIN170 = 170.0                # ns, overlap window of the predicted random pile-up (+-85 ns around the main pulse)
INTIME = (140.0, 156.0)       # ns, in-time onset range of the production trigger
THRS = (10.0, 20.0)           # ADC, amplitude thresholds of the rates and of the measured k >= 2
KEEP = 10.0                   # ADC, kept channel (pileup_rates.KEEP_ADC)
PF_T = 120.0                  # ns, first onset after the 30 pedestal samples
MIN_WF = 5                    # waveforms per module for a per-module statistic
RING_BANDS = [(1, 2), (3, 4), (5, 6), (7, 9), (10, 13), (14, 16)]
BLOCK_MIN = 10.0


def load_runs(d):
    """Concatenate the scan files of a directory: dict of waveform arrays (w_*, with a global event index and file number), event arrays, noise reference arrays, sums."""
    fs = sorted(glob.glob(str(Path(d) / "*.npz")))
    if not fs:
        sys.exit(f"[ERROR] no scan files in {d}")
    Z = [np.load(f, allow_pickle=False) for f in fs]
    out = {"files": [str(z["files"][0]) for z in Z], "names": Z[0]["mod_name"], "pwo4": Z[0]["mod_pwo4"]}
    off = np.cumsum([0] + [len(z["ev_ts"]) for z in Z])
    for key in Z[0].files:
        if key.startswith("w_"):
            out[key] = np.concatenate([z[key] for z in Z])
    out["w_ev"] = np.concatenate([z["w_ev"].astype(np.int64) + o for z, o in zip(Z, off[:-1])])
    out["w_gf"] = np.concatenate([np.full(z["w_mi"].size, i, np.int16) for i, z in enumerate(Z)])
    for key in ("ev_ts", "ev_cls", "ev_nch", "ev_num"):
        out[key] = np.concatenate([z[key] for z in Z])
    out["ev_gf"] = np.concatenate([np.full(len(z["ev_ts"]), i, np.int16) for i, z in enumerate(Z)])
    out["n_seen"] = np.array([int(z["n_seen"]) for z in Z])
    out["wall"] = np.array([float(z["wall_s"]) for z in Z])
    out["sig_ref"] = np.stack([z["mod_sig_ref"][0] for z in Z])
    out["b_ref"] = np.stack([z["mod_b_ref"][0] for z in Z])
    out["n_ref"] = np.stack([z["mod_n_ref"][0] for z in Z])
    bits = collections.Counter()
    for z in Z:
        for b, n in zip(z["bits_val"], z["bits_n"]):
            bits[int(b)] += int(n)
    out["bits"] = bits
    for key in ("G_n", "G_sx", "G_sxx", "G_sv"):
        out[key] = sum(z[key] for z in Z)
    exo = collections.defaultdict(list)
    for z in Z:
        for key in z.files:
            if key.startswith("ex_"):
                exo[key].append(z[key])
    out["ex"] = {k: np.concatenate(v) for k, v in exo.items()}
    return out


def first_onset(t):
    with np.errstate(all="ignore"):
        return np.nanmin(np.where(np.isfinite(t), t, np.inf), axis=1)


def ring_of(names, geo):
    """Chebyshev ring index (1 = innermost PbWO4 ring around the beam hole) per module; 0 for modules that are not PbWO4."""
    r = np.zeros(len(names), np.int64)
    for i, n in enumerate(names):
        g = geo.get(str(n))
        if g is not None and str(n).startswith("W"):
            r[i] = int(np.floor(max(abs(g["x"]), abs(g["y"])) / g["sx"]))
    return r


def report(scan_dir, stem, ref_sigma):
    import json
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle
    t_start = time.monotonic()
    D = load_runs(Path(scan_dir) / "scan")
    P = load_runs(Path(scan_dir) / "ref24590")
    names, pwo4 = D["names"], D["pwo4"]
    nm = len(names)
    geo = {m["n"]: m["geo"] for m in json.load(open(Path(os.environ["PRAD2_DATABASE_DIR"]) / "hycal_map.json"))}
    ring = ring_of(names, geo)
    band = np.zeros(nm, np.int64)
    for j, (a, b) in enumerate(RING_BANDS):
        band[(ring >= a) & (ring <= b)] = j + 1
    band_lbl = [f"ring {a}-{b}" for a, b in RING_BANDS]
    rows = []                                                       # TSV: section, item, key, value, err, n, note

    def add(section, item, key, value, err="", n="", note=""):
        rows.append((section, item, key, value, err, n, note))

    def mk(R, cls):
        """Waveform mask: class cls, PbWO4, not in an array-wide event; also returns per-class event count and vetoed events."""
        ev_ok = (R["ev_cls"] == cls) & (R["ev_nch"] < 300)
        m = (R["w_cls"] == cls) & R["pwo4"][R["w_mi"]] & ev_ok[R["w_ev"]]
        return m, int(ev_ok.sum()), int(((R["ev_cls"] == cls) & (R["ev_nch"] >= 300)).sum())
    for R in (D, P):
        R["pwo4"] = pwo4
    mR, nR, vR = mk(D, 0)                                           # random trigger, run 24578
    mQ, nQ, vQ = mk(D, 1)                                           # production-trigger events of run 24578
    mP, nP, vP = mk(P, 1)                                           # production run 24590
    notes = []

    def say(s=""):
        notes.append(s)
        print(s, flush=True)
    span = {}
    for tag, R in (("24578", D), ("24590", P)):
        span[tag] = (R["ev_ts"].min() * 4e-9, R["ev_ts"].max() * 4e-9)
    say(f"random_trigger_study report: run 24578, {len(D['files'])} EVIO splits ({', '.join(f[-5:] for f in D['files'])}), {int(D['n_seen'].sum())} events read (all events of each split); reference run 24590, {len(P['files'])} splits ({', '.join(f[-5:] for f in P['files'])}), {int(P['n_seen'].sum())} events (first 60000 of each split).")
    say("Event fraction used: 24578 100 % of the events of 10 local splits (splits 300, 600, ..., 3000 of the run: every 300th split, about 3 % of the run, ~11 min apart); 24590: 2 splits, first 60000 events each (~50 % of each split).")
    say(f"Waveform source: PbWO4 (W*) only, default fit, run-type bank 0p73GeV_H2 for both runs; noise reference: 24590 per file as the production scan, 24578 the median of the 24590 files' sigma_ref and b_ref per module (the random trigger gives < 20 waveforms per module in the first 20000 events); array-wide events (>= 300 channels) removed: 24578 {vR} random + {vQ} production-trigger events, 24590 {vP} events.")
    ntot = sum(D["bits"].values())
    say("")
    say("== 1. Trigger class and onset times ==")
    bits_txt = ", ".join(f"0x{b:x}: {n} ({100 * n / ntot:.2f} %)" for b, n in D["bits"].most_common(6))
    say(f"24578 trigger words (events): {bits_txt}")
    add("trigger", "24578", "n_events", ntot)
    for b, n in D["bits"].most_common():
        add("trigger", "24578", f"bits_0x{b:x}", n, "", ntot)
    ne_ch = D["ev_nch"][(D["ev_cls"] == 0) & (D["ev_nch"] < 300)]
    say(f"24578 random-trigger events (bits 0): {nR} used; events with no HyCal channel (all channels below the readout threshold): {100 * np.mean(ne_ch == 0):.2f} %; mean channels per event {ne_ch.mean():.3f}; production 24590: {100 * np.mean(P['ev_nch'][P['ev_cls'] == 1] == 0):.2f} % empty, mean {P['ev_nch'][(P['ev_cls'] == 1) & (P['ev_nch'] < 300)].mean():.2f} channels per event")
    dts = []
    for g in range(len(D["files"])):
        t = D["ev_ts"][(D["ev_gf"] == g)]
        dts.append(np.diff(t) * 4e-3)                               # us (4 ns ticks)
    dts = np.concatenate(dts)
    dts_all = dts
    dts = dts[dts < 1e4]
    say(f"24578 inter-event time (all classes, per split): median {np.median(dts):.1f} us, mean {dts.mean():.1f} us, i.e. event rate {1e3 / dts.mean():.1f} kHz; std/mean = {dts.std() / dts.mean():.2f} (1 for Poisson; gaps > 10 ms dropped: {int((dts_all >= 1e4).sum())})")
    add("trigger", "24578", "event_rate_kHz", 1e3 / dts.mean(), "", len(dts), "1/mean inter-event time, 4 ns clock ticks assumed")
    kmax = 4
    kdist = {}
    for tag, m, R in (("24578 random (bits 0)", mR, D), ("24578 production trigger (0x100)", mQ, D), ("24590 production (0x100)", mP, P)):
        k = R["w_k"][m & ~R["w_bad"]]
        kdist[tag] = np.bincount(np.minimum(k, kmax), minlength=kmax + 1) / max(k.size, 1)
        ncl = k.size
        say(f"pulses per stored waveform {tag}: n_wf {ncl}; k = 0 {100 * kdist[tag][0]:.2f} %, 1 {100 * kdist[tag][1]:.2f} %, 2 {100 * kdist[tag][2]:.2f} %, 3 {100 * kdist[tag][3]:.2f} %, >= 4 {100 * kdist[tag][4]:.2f} %")
        for j in range(kmax + 1):
            add("k_distribution", tag, f"k={j}" if j < kmax else f"k>={kmax}", kdist[tag][j], "", ncl)
    intime = {}
    for tag, m, R in (("24578 random (bits 0)", mR, D), ("24578 production trigger (0x100)", mQ, D), ("24590 production (0x100)", mP, P)):
        mm = m & ~R["w_bad"] & (R["w_k"] >= 1)
        t = R["w_t"][mm]
        a = R["w_a"][mm]
        inw = (t >= INTIME[0]) & (t <= INTIME[1]) & (a > 0)
        fw = inw.any(axis=1).mean() if t.shape[0] else float("nan")
        fp = inw.sum() / max((a > 0).sum(), 1)
        sel = (a >= 20) & np.isfinite(t)
        fp20 = (inw & (a >= 20)).sum() / max(sel.sum(), 1)
        intime[tag] = (fw, fp, fp20)
        say(f"in-time onset fraction ({INTIME[0]:.0f}-{INTIME[1]:.0f} ns) {tag}: waveforms with k >= 1 and a pulse in the range {100 * fw:.2f} % ({t.shape[0]} wf); pulses {100 * fp:.2f} % (A >= 20 ADC: {100 * fp20:.2f} % of {int(sel.sum())})")
        add("in_time", tag, "frac_waveforms_k>=1", fw, "", t.shape[0], f"onset {INTIME[0]:.0f}-{INTIME[1]:.0f} ns")
        add("in_time", tag, "frac_pulses", fp, "", int((a > 0).sum()))
        add("in_time", tag, "frac_pulses_A>=20", fp20, "", int(sel.sum()))
    flat_exp = (INTIME[1] - INTIME[0]) / FID_ANY_NS
    say(f"flat-density expectation for the range: {100 * flat_exp:.2f} % of the pulses (record length with a complete fit {FID_ANY_NS:.0f} ns)")
    # onset flatness
    edges = np.arange(0, 401, 16.0)
    hist_t = {}
    for tag, m, R in (("24578 random (bits 0)", mR, D), ("24590 production (0x100)", mP, P)):
        mm = m & ~R["w_bad"]
        t = R["w_t"][mm]
        a = R["w_a"][mm]
        ok = (a >= 20) & np.isfinite(t)
        hist_t[tag] = np.histogram(t[ok], edges)[0]
    h = hist_t["24578 random (bits 0)"]
    c = 0.5 * (edges[1:] + edges[:-1])
    fidm = (c >= FID[0]) & (c <= FID[1])

    def flat_stats(h_, lo, hi):
        m_ = (c >= lo) & (c <= hi)
        mean_ = h_[m_].mean()
        chi_ = ((h_[m_] - mean_) ** 2 / np.maximum(mean_, 1)).sum() / (m_.sum() - 1)
        sl_ = np.polyfit(c[m_], h_[m_] / mean_, 1)[0] * 100.0
        return mean_, chi_, sl_, int(h_[m_].sum()), int(m_.sum() - 1)
    # finer, amplitude-split densities (8 ns bins) for the depletion windows, from the vetoed PbWO4 random waveforms
    e8 = np.arange(0, 401, 8.0)
    c8 = 0.5 * (e8[1:] + e8[:-1])
    dens = {}
    for lo_a, hi_a in ((20, 50), (50, 200), (200, 1e9)):
        mm = mR & ~D["w_bad"]
        tt_, aa_ = D["w_t"][mm], D["w_a"][mm]
        ok_ = (aa_ >= lo_a) & (aa_ < hi_a) & np.isfinite(tt_)
        dens[(lo_a, hi_a)] = np.histogram(tt_[ok_], e8)[0]
    mean_h, chi, slope, nflat, ndf = flat_stats(h, *FID)
    _, chi2_, slope2, n2_, ndf2 = flat_stats(h, *PLAT2)
    say(f"onset flatness of random pulses with A >= 20 ADC (16 ns bins, array-wide events removed): late plateau {FID[0]:.0f}-{FID[1]:.0f} ns: chi2/ndf vs a constant {chi:.2f} (ndf {ndf}, n = {nflat}), relative slope {slope:+.2f} % per 100 ns; early plateau {PLAT2[0]:.0f}-{PLAT2[1]:.0f} ns: chi2/ndf {chi2_:.2f}, density relative to the late plateau {100 * h[(c >= PLAT2[0]) & (c <= PLAT2[1])].mean() / mean_h:.0f} %")
    add("onset_flat", "24578 random", "chi2_ndf_late_plateau", chi, "", ndf, f"A>=20 ADC, {FID[0]:.0f}-{FID[1]:.0f} ns, 16 ns bins")
    add("onset_flat", "24578 random", "slope_late_plateau_pct_per_100ns", slope, "", nflat)
    add("onset_flat", "24578 random", "chi2_ndf_early_plateau", chi2_, "", ndf2, f"{PLAT2[0]:.0f}-{PLAT2[1]:.0f} ns")
    for (lo_a, hi_a), hh in dens.items():
        pl = hh[(c8 >= FID[0]) & (c8 <= FID[1])].mean()
        pe = hh[(c8 >= PLAT2[0]) & (c8 <= PLAT2[1])].mean()
        parts = [f"{d[0]:.0f}-{d[1]:.0f} ns: {100 * hh[(c8 >= d[0]) & (c8 < d[1])].mean() / pl:.0f} %" for d in DEPL]
        lab = f"A {lo_a:.0f}-{hi_a:.0f}" if hi_a < 1e8 else f"A >= {lo_a:.0f}"
        say(f"   {lab} ADC ({int(hh.sum())} pulses): density relative to the late plateau: early plateau {100 * pe / pl:.0f} %, depleted windows " + ", ".join(parts))
        for d in DEPL:
            add("onset_depletion", lab, f"{d[0]:.0f}-{d[1]:.0f}ns_rel_to_late_plateau", float(hh[(c8 >= d[0]) & (c8 < d[1])].mean() / pl), "", int(hh.sum()))
        add("onset_depletion", lab, "early_plateau_rel_to_late_plateau", float(pe / pl), "", int(hh.sum()))
    say(f"   -> the random sample (trigger word 0) lacks pulses with onset 120-160 ns (they fire the production trigger, word 0x100) and, for the larger pulses, with onset 0-48 ns; the rate below uses the late plateau only. fraction of A >= 20 pulses with onset > 360 ns: {100 * h[c > FID[1]].sum() / h.sum():.1f} %")
    for j in range(len(c)):
        add("onset_hist", "24578 random", f"{edges[j]:.0f}-{edges[j + 1]:.0f}", int(h[j]), "", "", "A>=20 ADC pulses")
        add("onset_hist", "24590 production", f"{edges[j]:.0f}-{edges[j + 1]:.0f}", int(hist_t["24590 production (0x100)"][j]), "", "", "A>=20 ADC pulses")

    # ------------------------------------------------------------------ 2. pedestal and noise
    say("")
    say("== 2. Pedestal and noise (PbWO4, waveforms whose first 30 samples are pulse-free: k == 0 with max - baseline < 8 ADC (PRECUT_ADC) or first onset >= %.0f ns, not bad, no overflow) ==" % PF_T)

    def pf_mask(R, m):
        t1 = first_onset(R["w_t"])
        return m & ~R["w_bad"] & ~R["w_ovf"] & (((R["w_k"] == 0) & (R["w_raw"] < 8.0)) | (t1 >= PF_T))

    fac = np.sqrt(30.0 / 29.0)
    pfR, pfP = pf_mask(D, mR), pf_mask(P, mP)

    def per_module(R, m, vals_fn):
        """Per-module median of vals over waveforms of mask m; count."""
        mi = R["w_mi"][m]
        v = vals_fn(R)[m]
        o = np.argsort(mi, kind="stable")
        mi, v = mi[o], v[o]
        cnt = np.bincount(mi, minlength=nm)
        end = np.cumsum(cnt)
        med = np.full(nm, np.nan)
        sd = np.full(nm, np.nan)
        for j in np.nonzero(cnt >= MIN_WF)[0]:
            x = v[end[j] - cnt[j]:end[j]]
            med[j] = np.median(x)
            sd[j] = 1.4826 * np.median(np.abs(x - med[j]))
        return med, sd, cnt
    ped_m, ped_sd, cnt_pf = per_module(D, pfR, lambda R: R["w_ped"])
    rms_m, rms_sd, _ = per_module(D, pfR, lambda R: R["w_s30"] * fac)
    mad_m, _, _ = per_module(D, pfR, lambda R: R["w_mad30"])
    rmsP, _, cnt_pfP = per_module(P, pfP, lambda R: R["w_s30"] * fac)
    pedP, _, _ = per_module(P, pfP, lambda R: R["w_ped"])
    k0R = mR & ~D["w_bad"] & ~D["w_ovf"] & (D["w_k"] == 0)
    fs_m, _, cnt_k0 = per_module(D, k0R, lambda R: R["w_fs"])
    s30_k0, _, _ = per_module(D, k0R, lambda R: R["w_s30"] * fac)
    ref_files = np.load(ref_sigma, allow_pickle=False)
    sref_prod = np.nanmedian(np.where(ref_files["mod_n_ref"] >= 20, ref_files["mod_sig_ref"], np.nan), axis=0)
    bref_prod = np.nanmedian(np.where(ref_files["mod_n_ref"] >= 20, ref_files["mod_b_ref"], np.nan), axis=0)
    sref_own = np.nanmedian(np.where(D["n_ref"] >= 20, D["sig_ref"], np.nan), axis=0)
    sref_used = np.median(D["sig_ref"], axis=0)
    sref_prod_all = np.median(ref_files["mod_sig_ref"], axis=0)
    W = pwo4
    ok = W & np.isfinite(rms_m)
    q = lambda x, p: np.nanpercentile(x, p)       # noqa: E731
    say(f"pulse-free waveforms of 24578 random: {int(pfR.sum())} (of {int((mR & ~D['w_bad']).sum())} stored PbWO4 waveforms); per-module statistics for modules with >= {MIN_WF}: {int(ok.sum())} of {int(W.sum())} PbWO4 modules (median {int(np.median(cnt_pf[W]))} waveforms per module, 16-84 %: {q(cnt_pf[W], 16):.0f}-{q(cnt_pf[W], 84):.0f})")
    say(f"baseline rms of the first 30 samples (median over a module's waveforms x sqrt(30/29)): median over modules {np.nanmedian(rms_m[ok]):.3f} ADC, 16-84 %: {q(rms_m[ok], 16):.3f}-{q(rms_m[ok], 84):.3f}, 2.5-97.5 %: {q(rms_m[ok], 2.5):.3f}-{q(rms_m[ok], 97.5):.3f}; std over modules {np.nanstd(rms_m[ok]):.3f} ({100 * np.nanstd(rms_m[ok]) / np.nanmedian(rms_m[ok]):.1f} % of the median)")
    say(f"   MAD-based (1.4826 MAD of the 30 samples, median over waveforms): median {np.nanmedian(mad_m[ok]):.3f} ADC (quantised: the samples are integers, MAD = 1 ADC gives 1.483; use the rms)")
    say(f"baseline mean (clipped pedestal, median over a module's waveforms): median over modules {np.nanmedian(ped_m[ok]):.2f} ADC, 16-84 %: {q(ped_m[ok], 16):.2f}-{q(ped_m[ok], 84):.2f}, std over modules {np.nanstd(ped_m[ok]):.2f} ADC; within-module spread of the per-waveform baseline (1.4826 MAD): median {np.nanmedian(ped_sd[ok]):.3f} ADC")
    say(f"whole-record (100 samples) std of the waveforms with no pulse found (k == 0): median over modules with >= {MIN_WF} such waveforms {np.nanmedian(fs_m[W & np.isfinite(fs_m)]):.3f} ADC (n modules {int((W & np.isfinite(fs_m)).sum())}; {int(k0R.sum())} waveforms): not a noise measure, see the next line")
    k0Y = D["ex"]["ex_k0_Y"].astype(float)
    k0ped = np.median(k0Y[:, :30], axis=1)
    mad100 = 1.4826 * np.median(np.abs(k0Y - np.median(k0Y, axis=1)[:, None]), axis=1)
    late = np.argmax(k0Y, axis=1) >= 90
    rawk0 = k0Y.max(axis=1) - k0ped
    say(f"   the {len(k0Y)} random k == 0 example waveforms (reservoir): 1.4826 MAD of the 100 samples median {np.median(mad100):.3f} ADC (quantised, integer samples), std {np.median(k0Y.std(axis=1)):.3f} ADC; max sample - baseline median {np.median(rawk0):.1f} ADC (16-84 %: {np.percentile(rawk0, 16):.1f}-{np.percentile(rawk0, 84):.1f}); {100 * late.mean():.0f} % have their maximum in the last 10 samples (pulse onset at the end of the record, rejected by the fit): these are not pulse-free records; true pulse-free records are rare")
    add("noise", "24578 k0 reservoir", "mad100_median", float(np.median(mad100)), "", len(k0Y))
    add("noise", "24578 k0 reservoir", "frac_max_in_last_10_samples", float(late.mean()), "", len(k0Y))
    for tag, v in (("sigma_ref used by the fitter for 0p73GeV_H2 (run 24590 scans, median over files)", sref_prod_all), ("sigma_ref of the 24578 files' own reference (NOT used by the fit: < 20 waveforms per module in 20000 events; the fit uses the 24590 reference)", sref_used)):
        okk = ok & np.isfinite(v)
        r = rms_m[okk] / v[okk]
        say(f"vs {tag}: sigma_ref median {np.nanmedian(v[okk]):.3f} ADC; rms/sigma_ref median {np.median(r):.3f}, 16-84 %: {q(r, 16):.3f}-{q(r, 84):.3f}")
        add("noise", tag, "median_rms_over_sigma_ref", float(np.median(r)), "", int(okk.sum()))
    okp = ok & np.isfinite(rmsP)
    say(f"24590 production pulse-free first-30 rms (same recipe, modules with >= {MIN_WF}): median {np.nanmedian(rmsP[W & np.isfinite(rmsP)]):.3f} ADC ({int((W & np.isfinite(rmsP)).sum())} modules, {int(pfP.sum())} waveforms); per module 24578 / 24590 rms median {np.median(rms_m[okp] / rmsP[okp]):.3f}, 16-84 %: {q(rms_m[okp] / rmsP[okp], 16):.3f}-{q(rms_m[okp] / rmsP[okp], 84):.3f}")
    say(f"selection-bias check (random trigger, per-module medians of the first-30 rms x sqrt(30/29)): waveforms with a late pulse or a pulse-free record: {np.nanmedian(rms_m[ok]):.3f}; k == 0 waveforms only: {np.nanmedian(s30_k0[W & np.isfinite(s30_k0)]):.3f} ADC (n modules {int((W & np.isfinite(s30_k0)).sum())}); sigma_ref (files of 24578, from all stored waveforms of the first 20000 events): {np.nanmedian(sref_own[W]):.3f}")
    rel_err = 1.25 / np.sqrt(2.0 * 29.0 * np.maximum(cnt_pf, 1))
    intr = np.sqrt(max(np.nanvar(rms_m[ok] / np.nanmedian(rms_m[ok])) - np.mean(rel_err[ok] ** 2), 0.0))
    say(f"spread of the per-module rms over modules: {100 * np.nanstd(rms_m[ok]) / np.nanmedian(rms_m[ok]):.1f} % of the median; expected from the statistics of the per-module median ({MIN_WF}+ waveforms, 30 samples each) {100 * np.sqrt(np.mean(rel_err[ok] ** 2)):.1f} %; intrinsic module-to-module spread {100 * intr:.1f} %")
    add("noise", "24578", "rms_spread_pct_of_median_observed", float(100 * np.nanstd(rms_m[ok]) / np.nanmedian(rms_m[ok])), "", int(ok.sum()))
    add("noise", "24578", "rms_spread_pct_intrinsic", float(100 * intr), "", int(ok.sum()))
    add("noise", "24578", "rms30_median_over_modules", float(np.nanmedian(rms_m[ok])), "", int(ok.sum()))
    for tag, R, pf in (("24578 random", D, pfR), ("24590 production", P, pfP)):
        okw = pf & np.isfinite(sref_prod_all[R["w_mi"]])
        rw = R["w_s30"][okw] * fac / sref_prod_all[R["w_mi"][okw]]
        say(f"per waveform (all {int(okw.sum())} pulse-free waveforms of {tag}, every module): first-30 rms / sigma_ref(0p73GeV_H2) median {np.median(rw):.3f}, 16-84 %: {q(rw, 16):.3f}-{q(rw, 84):.3f}; baseline minus b_ref(24590): median {np.nanmedian(R['w_ped'][okw] - bref_prod[R['w_mi'][okw]]):+.2f} ADC (modules with a b_ref)")
        add("noise", tag, "waveform_rms30_over_sigma_ref_median", float(np.median(rw)), "", int(okw.sum()))
    okb = ok & np.isfinite(bref_prod)
    db = ped_m[okb] - bref_prod[okb]
    say(f"baseline mean 24578 minus b_ref of 24590 (median of the files' module medians): median {np.median(db):.2f} ADC, 16-84 %: {q(db, 16):.2f}-{q(db, 84):.2f} (n {int(okb.sum())} modules)")
    for i in np.nonzero(W)[0]:
        add("module", str(names[i]), "n_pulse_free", int(cnt_pf[i]))
        for key, v in (("baseline_mean", ped_m), ("baseline_spread_mad", ped_sd), ("rms30", rms_m), ("rms_mad30", mad_m), ("rms100_k0", fs_m), ("sigma_ref_0p73_prod", sref_prod_all), ("sigma_ref_24578", sref_used), ("b_ref_24590", bref_prod), ("rms30_24590", rmsP)):
            if np.isfinite(v[i]):
                add("module", str(names[i]), key, float(v[i]))
    # drift per file / 10-min block
    tmin = D["ev_ts"] * 4e-9 / 60.0
    t0min = tmin.min()
    drift = []
    ped_ref = np.where(np.isfinite(ped_m), ped_m, np.nan)
    rms_ref = np.where(np.isfinite(rms_m), rms_m, np.nan)
    for g in range(len(D["files"])):
        mg = pfR & (D["w_gf"] == g)
        mi = D["w_mi"][mg]
        dped = D["w_ped"][mg] - ped_ref[mi]
        drms = D["w_s30"][mg] * fac / rms_ref[mi]
        okm = np.isfinite(dped) & np.isfinite(drms)
        # per module medians within the file, then median over modules (>= 3 waveforms): robust to hot modules
        u, inv = np.unique(mi[okm], return_inverse=True)
        cnt = np.bincount(inv)
        sel = cnt >= 3
        md = np.array([np.median(dped[okm][inv == j]) for j in np.nonzero(sel)[0]])
        mr = np.array([np.median(drms[okm][inv == j]) for j in np.nonzero(sel)[0]])
        tg = tmin[D["ev_gf"] == g]
        blk = int(np.floor(((tg.mean() - t0min)) / BLOCK_MIN))
        sem_d = 1.2533 * 1.4826 * np.median(np.abs(md - np.median(md))) / np.sqrt(max(md.size, 1))
        sem_r = 1.2533 * 1.4826 * np.median(np.abs(mr - np.median(mr))) / np.sqrt(max(mr.size, 1))
        drift.append((g, float(tg.mean() - t0min), blk, int(mg.sum()), int(sel.sum()), float(np.median(md)), float(sem_d), float(np.median(mr)), float(sem_r)))
        add("drift", D["files"][g], "t_min_from_first", float(tg.mean() - t0min), "", int(mg.sum()), f"{BLOCK_MIN:.0f}-min block {blk}")
        add("drift", D["files"][g], "median_dBaseline_ADC", float(np.median(md)), float(sem_d), int(sel.sum()), "median over modules (>= 3 waveforms) of the per-module median minus the run median")
        add("drift", D["files"][g], "median_rms_ratio", float(np.median(mr)), float(sem_r), int(sel.sum()), "same, rms / run median rms")
    dd = np.array([d[5] for d in drift]); rr = np.array([d[7] for d in drift])
    say(f"drift over the run (one split per {BLOCK_MIN:.0f}-min block, {len(drift)} blocks, t = {drift[0][1]:.0f} ... {drift[-1][1]:.0f} min): baseline offset relative to the run median per file: min {dd.min():+.3f}, max {dd.max():+.3f} ADC, std {dd.std():.3f} ADC (typical sem {np.median([d[6] for d in drift]):.3f}); rms ratio min {rr.min():.3f}, max {rr.max():.3f}, std {rr.std():.3f}")

    # ------------------------------------------------------------------ 3. occupancy and pile-up
    say("")
    say("== 3. Occupancy and predicted random pile-up (24578 random trigger, PbWO4) ==")
    good = mR & ~D["w_bad"]
    a_all = D["w_a"]
    t_all = D["w_t"]
    wf_k1 = (D["w_k"] >= 1)
    wf_k2 = (D["w_k"] >= 2)
    nwf = int(good.sum())
    say(f"stored PbWO4 waveforms {nwf} in {nR} random events ({nwf / nR / W.sum() * 1e4:.3f} per 10^4 module-events); with >= 1 pulse (fit k >= 1, A >= 5 ADC): {100 * np.mean(wf_k1[good]):.2f} %; >= 2 pulses: {100 * np.mean(wf_k2[good]):.3f} % ({int(wf_k2[good].sum())})")
    expo = nR * (FID[1] - FID[0]) * 1e-3           # events x us per module
    lam = {}
    for thr in THRS:
        sel = good[:, None] & (a_all >= thr) & (t_all >= FID[0]) & (t_all <= FID[1])
        npul = np.bincount(np.repeat(D["w_mi"], KMAX)[sel.ravel()], minlength=nm)
        lam[thr] = npul / expo
        p1 = np.bincount(D["w_mi"][good & (np.where(sel, a_all, 0).max(axis=1) >= thr)], minlength=nm) / nR
        tot = npul[W].sum()
        say(f"A >= {thr:.0f} ADC, onset {FID[0]:.0f}-{FID[1]:.0f} ns: {int(tot)} pulses in {nR} events x {int(W.sum())} modules; per-module rate mean {tot / expo / W.sum() * 1e3:.4f} per ms (= {tot / expo / W.sum():.2e} per us), median over modules {np.median(lam[thr][W]) * 1e3:.4f} per ms; array total {tot / expo:.3f} pulses per us; P(>= 1 pulse with onset in the plateau window) mean over modules {np.mean(p1[W]):.2e}")
        add("occupancy", f"A>={thr:.0f}", "pulses_total", int(tot), "", nR, f"onset {FID[0]:.0f}-{FID[1]:.0f} ns")
        add("occupancy", f"A>={thr:.0f}", "rate_per_module_mean_per_us", float(tot / expo / W.sum()), float(np.sqrt(tot) / expo / W.sum()), int(tot))
        add("occupancy", f"A>={thr:.0f}", "rate_per_module_median_per_us", float(np.median(lam[thr][W])))
        add("occupancy", f"A>={thr:.0f}", "array_rate_per_us", float(tot / expo), float(np.sqrt(tot) / expo), int(tot))
        add("occupancy", f"A>={thr:.0f}", "P_ge1_per_module_event_mean", float(np.mean(p1[W])))
        for i in np.nonzero(W)[0]:
            add("module", str(names[i]), f"npulses_A>={thr:.0f}", int(npul[i]))
            add("module", str(names[i]), f"rate_per_us_A>={thr:.0f}", float(lam[thr][i]))
    a1 = good & (D["w_k"] >= 1)
    mA = a_all.max(axis=1)
    for thr in THRS:
        n1 = int((good & (mA >= thr)).sum()); n2 = int((good & (np.sort(a_all, axis=1)[:, -2] >= thr)).sum())
        say(f"waveform level, A >= {thr:.0f} ADC: P(>= 1 pulse) among {nwf} stored waveforms {100 * n1 / nwf:.2f} %, P(>= 2) {100 * n2 / nwf:.3f} % ({n2}); per random event and module P(>= 1) = {n1 / nR / W.sum():.2e}, P(>= 2) = {n2 / nR / W.sum():.2e}")
        add("occupancy", f"A>={thr:.0f}", "P_ge1_per_stored_waveform", n1 / nwf, "", nwf)
        add("occupancy", f"A>={thr:.0f}", "P_ge2_per_stored_waveform", n2 / nwf, "", nwf)
    # per band rates
    for thr in THRS:
        for j, bl in enumerate(band_lbl):
            mm = band == j + 1
            tot = lam[thr][mm].sum() * expo
            add("ring_rate", bl, f"rate_per_module_per_us_A>={thr:.0f}", float(tot / expo / mm.sum()), float(np.sqrt(max(tot, 1)) / expo / mm.sum()), int(tot), f"{int(mm.sum())} modules")
    hot = np.argsort(-lam[THRS[0]] * W)[:8]
    say("hottest modules (rate per us, A >= %.0f): " % THRS[0] + ", ".join(f"{names[i]} {lam[THRS[0]][i]:.3f}" for i in hot))
    # clustering: observed waveforms with >= 2 pulses (A >= thr, onset 20-360) vs independent (Poisson) pulses at the module's plateau rate
    for thr in THRS:
        sel2 = (a_all >= thr) & (t_all >= 20.0) & (t_all <= 360.0)
        n2o = int((good & (sel2.sum(axis=1) >= 2)).sum())
        mu = lam[thr] * FID_ANY_NS * 1e-3
        n2e = float((nR * (1.0 - np.exp(-mu) - mu * np.exp(-mu)))[W].sum())
        top = np.argsort(-np.bincount(D["w_mi"][good & (sel2.sum(axis=1) >= 2)], minlength=nm))[:5]
        say(f"clustering, A >= {thr:.0f} ADC: stored waveforms with >= 2 pulses (onset 20-360 ns): {n2o}; expected for independent pulses at each module's own plateau rate: {n2e:.2f}; ratio {n2o / max(n2e, 1e-9):.0f} (second pulses of the random sample are far from independent: clustered in a few modules); modules with most: " + ", ".join(f"{names[i]} {int(np.sum(good & (sel2.sum(axis=1) >= 2) & (D['w_mi'] == i)))}" for i in top))
        add("clustering", f"A>={thr:.0f}", "n_k2_observed", n2o, "", nR)
        add("clustering", f"A>={thr:.0f}", "n_k2_expected_poisson", n2e, "", nR)
    # measured production k >= 2 and prediction
    gP = mP & ~P["w_bad"]
    kept = gP & (P["w_a"].max(axis=1) >= KEEP)
    aP, tP = P["w_a"], P["w_t"]
    main = np.argmax(aP, axis=1)
    tm = tP[np.arange(len(main)), main]
    meas = {}
    wk = np.bincount(P["w_mi"][kept], minlength=nm).astype(float)                    # kept channels per module (production)
    res_rows = []
    for thr in THRS:
        oth = (aP >= thr) & (np.arange(KMAX)[None, :] != main[:, None])
        anyk = oth.any(axis=1)
        near = (oth & (np.abs(tP - tm[:, None]) < 0.5 * WIN170)).any(axis=1)
        for name, flag in (("any", anyk), ("170ns", near)):
            m_all = kept & flag
            cm = np.bincount(P["w_mi"][m_all], minlength=nm)
            meas[(thr, name)] = (cm, wk)
        T_fid = FID_ANY_NS * 1e-3
        for name, T in (("any", T_fid), ("170ns", WIN170 * 1e-3)):
            Pm = 1.0 - np.exp(-lam[thr] * T)
            pred = (wk * Pm)[W].sum() / wk[W].sum()
            cm, _ = meas[(thr, name)]
            mval = cm[W].sum() / wk[W].sum()
            err = np.sqrt(cm[W].sum()) / wk[W].sum()
            say(f"kept production channels (24590, A_main >= {KEEP:.0f} ADC, {int(wk[W].sum())} channels): other pulse A >= {thr:.0f} ADC {'anywhere in the record' if name == 'any' else 'within +-85 ns of the main pulse'}: measured {100 * mval:.3f} +- {100 * err:.3f} %; random-rate prediction (1 - exp(-rate x {T * 1e3:.0f} ns), module rates weighted by the production kept channels): {100 * pred:.4f} %; ratio predicted / measured {pred / mval:.4f}")
            add("pileup", f"A_other>={thr:.0f} {name}", "measured_k2_fraction_24590", float(mval), float(err), int(wk[W].sum()))
            add("pileup", f"A_other>={thr:.0f} {name}", "predicted_random_fraction", float(pred), "", int(wk[W].sum()))
            for j, bl in enumerate(band_lbl):
                mm = band == j + 1
                pr_ = (wk * Pm)[mm].sum() / max(wk[mm].sum(), 1)
                mv = cm[mm].sum() / max(wk[mm].sum(), 1)
                res_rows.append((thr, name, bl, mv, pr_, int(wk[mm].sum())))
                add("pileup_ring", f"A_other>={thr:.0f} {name}", bl, float(pr_), float(mv), int(wk[mm].sum()), "value = predicted random fraction; err column = measured k>=2 fraction (24590)")
    say("caveat: the random-trigger records are zero-suppressed, so a small pulse enters the rate only if it passes the readout threshold alone (spectrum rolls off at 10-13 ADC): the A >= 20 ADC comparison is the headline, A >= 10 ADC a lower bound of the prediction")
    k2all = np.sum(kdist["24590 production (0x100)"][2:])
    say(f"for reference, 24590 stored PbWO4 waveforms with k >= 2 (any amplitude >= 5 ADC): {100 * k2all:.2f} % (pileup_rates report, kept channels with A_main >= 10 ADC, run 24590 + 24603: 1.84 % PbWO4; the 5.0 % quoted in the brief was not found in the stored tables)")
    add("pileup", "24590 stored waveforms", "k>=2_fraction", float(k2all))

    # ------------------------------------------------------------------ 4. pulse shape and amplitude spectrum
    Gn, Gsx, Gsxx, Gsv = D["G_n"][0], D["G_sx"][0], D["G_sxx"][0], D["G_sv"][0]
    Pn, Psx, Psxx, Psv = P["G_n"][1], P["G_sx"][1], P["G_sxx"][1], P["G_sv"][1]
    say("")
    say("== 4. Pulse shape of isolated pulses (k == 1, PbWO4, no overflow) and amplitude spectrum ==")

    def stats(n, sx, sxx, sv):
        with np.errstate(invalid="ignore", divide="ignore"):
            mean = sx / n
            std = np.sqrt(np.maximum(sxx / n - mean ** 2, 0.0))
            return mean, std / np.sqrt(n), std, np.sqrt(sv / n)
    amp_lbl = [f"{int(AMP_EDGES[i])}-{int(AMP_EDGES[i + 1])}" for i in range(NAB)]
    inwin = (TAU >= -12.0) & (TAU <= 116.0)
    shape_txt = []
    for ab in range(NAB):
        for tag, n, sx, sxx, sv in (("24578 random", Gn[ab], Gsx[ab], Gsxx[ab], Gsv[ab]), ("24590 production", Pn[ab], Psx[ab], Psxx[ab], Psv[ab])):
            mean, sem, rms, asg = stats(n, sx, sxx, sv)
            okt = n > 100
            nwf_ = int(n[np.argmax(n)])
            coh = okt & (np.abs(mean) > 2 * sem)
            mx = np.nanmax(np.abs(mean[okt])) if okt.any() else float("nan")
            ratio = np.nanmedian((rms / asg)[okt & inwin]) if (okt & inwin).any() else float("nan")
            line = f"A {amp_lbl[ab]} {tag}: ~{nwf_} pulses; max |mean r/A| {100 * mx:.2f} %; bins with |mean| > 2 sem: {int(coh.sum())} of {int(okt.sum())}; rms / assumed sigma in the fit window median {ratio:.2f}"
            shape_txt.append(line)
            say(line)
            add("shape", f"{tag} A {amp_lbl[ab]}", "n_pulses_approx", nwf_)
            add("shape", f"{tag} A {amp_lbl[ab]}", "max_abs_mean_r_over_A", float(mx))
            add("shape", f"{tag} A {amp_lbl[ab]}", "rms_over_sigma_fitwindow_median", float(ratio))
            for j in range(NT):
                if n[j] > 0:
                    add("shape_profile", f"{tag} A {amp_lbl[ab]}", f"tau={TAU[j]:.0f}", float(mean[j]), float(sem[j]), int(n[j]), f"rms {rms[j]:.5f} assumed {asg[j]:.5f}")
    # amplitude spectrum of the random pulses (onset in the fiducial range)
    ae = np.logspace(0.5, 3.6, 32)
    sel = good[:, None] & (a_all > 0) & (t_all >= FID[0]) & (t_all <= FID[1])
    spec = np.histogram(a_all[sel], ae)[0]
    rate_a = spec / expo
    ac = np.sqrt(ae[1:] * ae[:-1])
    for j in range(len(ac)):
        add("amp_spectrum", "24578 random", f"{ae[j]:.1f}-{ae[j + 1]:.1f}", float(rate_a[j]), float(np.sqrt(max(spec[j], 1)) / expo), int(spec[j]), "array rate per us, all PbWO4 modules")
    peak = ac[np.argmax(spec)]
    say(f"amplitude spectrum of the random pulses (onset {FID[0]:.0f}-{FID[1]:.0f} ns): {int(spec.sum())} pulses, A median {np.median(a_all[sel]):.1f} ADC, fraction with A < 20: {100 * np.mean(a_all[sel] < 20):.1f} %, < 50: {100 * np.mean(a_all[sel] < 50):.1f} %, >= 100: {100 * np.mean(a_all[sel] >= 100):.1f} %, >= 500: {100 * np.mean(a_all[sel] >= 500):.2f} %; most populated bin {ae[np.argmax(spec)]:.0f}-{ae[np.argmax(spec) + 1]:.0f} ADC (= threshold roll-off region below)")

    # ------------------------------------------------------------------ examples
    ex = D["ex"]

    def pick(prefix, n, rng_):
        if f"{prefix}_key" not in ex:
            return np.array([], int)
        return np.argsort(ex[f"{prefix}_key"])[:n]
    k1_idx = []
    for b in range(len(EX_AMP_EDGES) - 1):
        pre = f"ex_k1_{b}"
        if f"{pre}_key" in ex:
            k1_idx += [(pre, j) for j in np.argsort(ex[f"{pre}_key"])[:4]]
    if len(k1_idx) < 20:                                           # fill from the bins that have more
        have = set(k1_idx)
        extra = []
        for b in range(len(EX_AMP_EDGES) - 1):
            pre = f"ex_k1_{b}"
            if f"{pre}_key" in ex:
                extra += [(pre, j) for j in np.argsort(ex[f"{pre}_key"])[4:] if (pre, j) not in have]
        k1_idx += extra[:20 - len(k1_idx)]
    k1_idx = k1_idx[:20]
    k0_idx = [("ex_k0", j) for j in pick("ex_k0", 10, None)]
    k2_idx = [("ex_k2", j) for j in pick("ex_k2", 10, None)]
    say(f"examples: {len(k1_idx)} single-pulse (4 per amplitude bin {EX_AMP_EDGES[:-1].astype(int).tolist()} ADC), {len(k0_idx)} pulse-free (k == 0, PbWO4, not bad), {len(k2_idx)} with k >= 2 (all from the random-trigger events of the 10 splits; available: k0 {len(ex.get('ex_k0_key', []))}, k2 {len(ex.get('ex_k2_key', []))})")
    for cat, lst in (("k>=1", k1_idx), ("k==0", k0_idx), ("k>=2", k2_idx)):
        for pre, j in lst:
            kk = int(ex[f"{pre}_k"][j])
            add("example", cat, f"{names[int(ex[f'{pre}_mi'][j])]} ev {int(ex[f'{pre}_evn'][j])}", kk, "", "",
                "t=" + ",".join(f"{x:.1f}" for x in ex[f"{pre}_t"][j][:kk]) + " A=" + ",".join(f"{x:.1f}" for x in ex[f"{pre}_a"][j][:kk]) + f" chi2={ex[f'{pre}_chi2'][j]:.2f}")

    # ------------------------------------------------------------------ figures
    plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": 0.3})
    rect = {}
    for i, nme in enumerate(names):
        g = geo.get(str(nme))
        if g is not None and pwo4[i]:
            rect[i] = Rectangle((g["x"] - g["sx"] / 2, g["y"] - g["sy"] / 2), g["sx"], g["sy"])

    def mod_map(ax, vals, title, vmin=None, vmax=None, log=False, label=""):
        from matplotlib.colors import LogNorm
        has = [i for i in rect if np.isfinite(vals[i]) and (vals[i] > 0 or not log)]
        no = [i for i in rect if i not in set(has)]
        ax.add_collection(PatchCollection([Rectangle(rect[i].get_xy(), rect[i].get_width(), rect[i].get_height()) for i in no], facecolor="0.85", edgecolor="0.7", lw=0.2))
        pc = PatchCollection([Rectangle(rect[i].get_xy(), rect[i].get_width(), rect[i].get_height()) for i in has], cmap="viridis", edgecolor="0.4", lw=0.15)
        pc.set_array(np.asarray(vals)[has])
        if log:
            pc.set_norm(LogNorm(vmin, vmax))
        else:
            pc.set_clim(vmin, vmax)
        ax.add_collection(pc)
        ax.set_xlim(-360, 360); ax.set_ylim(-360, 360); ax.set_aspect("equal"); ax.grid(False)
        ax.set_title(title, fontsize=8)
        plt.colorbar(pc, ax=ax, fraction=0.046, pad=0.02, label=label)

    def text_pages(pdf, title, lines, size=6.2, per_page=70):
        wrapped = [w for ln in lines for w in (textwrap.wrap(ln, 200, subsequent_indent="    ") or [""])]
        for p0 in range(0, len(wrapped), per_page):
            fig = plt.figure(figsize=(11, 8.5))
            fig.text(0.03, 0.97, f"{title} ({p0 // per_page + 1}/{(len(wrapped) - 1) // per_page + 1})", fontsize=10, weight="bold", va="top")
            for j, w in enumerate(wrapped[p0:p0 + per_page]):
                fig.text(0.03, 0.94 - 0.0128 * j, w, fontsize=size, va="top", family="monospace")
            pdf.savefig(fig); plt.close(fig)

    cols = ["#0072B2", "#D55E00", "#009E73", "#CC79A7"]
    pdf_path = Path(str(stem) + ".pdf")
    with PdfPages(pdf_path) as pdf:
        # summary pages
        text_pages(pdf, "Run 24578 random-trigger study: numbers", notes)
        # 1. trigger / k / onset
        fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
        a_ = ax[0, 0]
        bl = sorted(D["bits"].items(), key=lambda kv: -kv[1])[:8]
        a_.bar(range(len(bl)), [n for _, n in bl], color=cols[0]); a_.set_yscale("log")
        a_.set_xticks(range(len(bl))); a_.set_xticklabels([f"0x{b:x}" for b, _ in bl], rotation=45)
        a_.set_title("run 24578: events per trigger word (10 splits)"); a_.set_ylabel("events")
        a_ = ax[0, 1]
        a_.hist(dts, bins=np.logspace(0, 3.5, 60), color=cols[0]); a_.set_xscale("log"); a_.set_title("24578 inter-event time (all events, per split)"); a_.set_xlabel("us")
        a_ = ax[1, 0]
        w = 0.27
        for j, (tag, v) in enumerate(kdist.items()):
            a_.bar(np.arange(kmax + 1) + (j - 1) * w, 100 * v, w, label=tag, color=cols[j])
        a_.set_yscale("log"); a_.set_xticks(range(kmax + 1)); a_.set_xticklabels(["0", "1", "2", "3", ">=4"]); a_.set_xlabel("fitted pulses per stored waveform k"); a_.set_ylabel("% of stored PbWO4 waveforms"); a_.legend(fontsize=6)
        a_.set_title("pulses per waveform")
        a_ = ax[1, 1]
        for j, (tag, v) in enumerate(intime.items()):
            a_.bar(np.arange(3) + (j - 1) * w, [100 * x for x in v], w, label=tag, color=cols[j])
        a_.axhline(100 * flat_exp, color="k", ls="--", lw=0.8, label=f"flat density ({100 * flat_exp:.1f} %)")
        a_.set_xticks(range(3)); a_.set_xticklabels(["waveforms k>=1 with a pulse\nin 140-156 ns", "pulses in 140-156 ns", "pulses A>=20 in 140-156 ns"], fontsize=6)
        a_.set_ylabel("%"); a_.legend(fontsize=6); a_.set_title("in-time onset fraction")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        # onset times
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.5))
        for j, (tag, hh) in enumerate(hist_t.items()):
            ax[0].step(c, hh / hh.sum() / 16.0, where="mid", label=tag, color=cols[j])
            ax[1].step(c, hh / (hh[fidm].mean() if tag.startswith("24578") else hh.max()), where="mid", label=tag + (" / plateau" if tag.startswith("24578") else " / max"), color=cols[j])
        ax[0].set_yscale("log"); ax[0].set_xlabel("fitted onset (ns)"); ax[0].set_ylabel("pulses A >= 20 ADC per ns (normalised)"); ax[0].legend(fontsize=7)
        ax[1].axvspan(*FID, color="0.9"); ax[1].set_xlabel("fitted onset (ns)"); ax[1].set_ylabel("relative"); ax[1].legend(fontsize=7)
        ax[1].set_title(f"random: late plateau (grey) chi2/ndf {chi:.2f}; depleted at 0-48 and 120-160 ns (trigger bias)", fontsize=7)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        # 2. noise distributions
        fig, ax = plt.subplots(2, 3, figsize=(11, 7))
        okm = ok
        ax[0, 0].hist(rms_m[okm], bins=np.linspace(0.5, 3.0, 60), color=cols[0], alpha=0.8, label="24578 random")
        ax[0, 0].hist(rmsP[W & np.isfinite(rmsP)], bins=np.linspace(0.5, 3.0, 60), histtype="step", color=cols[1], label="24590 production (same recipe)")
        ax[0, 0].hist(sref_prod_all[W], bins=np.linspace(0.5, 3.0, 60), histtype="step", color="k", label="sigma_ref 0p73GeV_H2 (24590 scans)")
        ax[0, 0].set_xlabel("baseline rms, first 30 samples (ADC)"); ax[0, 0].set_ylabel("modules"); ax[0, 0].legend(fontsize=6)
        ax[0, 1].hist(rms_m[okm] / sref_prod_all[okm], bins=np.linspace(0.5, 2.5, 60), color=cols[0]); ax[0, 1].set_xlabel("rms / sigma_ref (0p73GeV_H2, 24590)"); ax[0, 1].set_ylabel("modules")
        ax[0, 2].hist(ped_m[okm], bins=60, color=cols[0]); ax[0, 2].set_xlabel("baseline mean (ADC)"); ax[0, 2].set_ylabel("modules")
        ax[1, 0].scatter(sref_prod_all[okm], rms_m[okm], s=3, color=cols[0]); ax[1, 0].plot([0.5, 3], [0.5, 3], "k--", lw=0.6); ax[1, 0].set_xlim(0.5, 3); ax[1, 0].set_ylim(0.5, 3)
        ax[1, 0].set_xlabel("sigma_ref 0p73GeV_H2 (ADC)"); ax[1, 0].set_ylabel("24578 pulse-free rms (ADC)")
        ax[1, 1].hist(ped_m[okb] - bref_prod[okb], bins=np.linspace(-5, 5, 60), color=cols[0]); ax[1, 1].set_xlabel("baseline 24578 - b_ref 24590 (ADC)"); ax[1, 1].set_ylabel("modules")
        ax[1, 2].hist(cnt_pf[W], bins=40, color=cols[0]); ax[1, 2].axvline(MIN_WF, color="k", ls="--"); ax[1, 2].set_xlabel("pulse-free waveforms per module"); ax[1, 2].set_ylabel("modules")
        fig.suptitle("Pedestal and noise per module (PbWO4; pulse-free first 30 samples)"); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        fig, ax = plt.subplots(1, 3, figsize=(11, 4.2))
        mod_map(ax[0], np.where(ok, rms_m, np.nan), "baseline rms 24578 (ADC)", 0.8, 2.0, label="ADC")
        mod_map(ax[1], np.where(ok & np.isfinite(sref_prod_all), rms_m / sref_prod_all, np.nan), "rms / sigma_ref(0p73GeV_H2)", 0.8, 1.6)
        mod_map(ax[2], np.where(ok, ped_m, np.nan), "baseline mean (ADC)", float(q(ped_m[ok], 2)), float(q(ped_m[ok], 98)), label="ADC")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        # drift
        fig, ax = plt.subplots(1, 3, figsize=(11, 3.8))
        tt = [d[1] for d in drift]
        ax[0].errorbar(tt, [d[5] for d in drift], [d[6] for d in drift], fmt="o", color=cols[0]); ax[0].axhline(0, color="k", lw=0.5); ax[0].set_xlabel("time since first split (min)"); ax[0].set_ylabel("baseline offset (ADC)")
        ax[1].errorbar(tt, [d[7] for d in drift], [d[8] for d in drift], fmt="o", color=cols[0]); ax[1].axhline(1, color="k", lw=0.5); ax[1].set_xlabel("time since first split (min)"); ax[1].set_ylabel("rms / run median")
        ax[2].bar(range(len(drift)), [d[3] for d in drift], color=cols[0]); ax[2].set_xticks(range(len(drift))); ax[2].set_xticklabels([f"{d[1]:.0f}" for d in drift], rotation=45); ax[2].set_ylabel("pulse-free waveforms"); ax[2].set_xlabel("minute")
        fig.suptitle(f"Drift over the run, one split per {BLOCK_MIN:.0f}-min block (median over modules of per-module medians)"); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        # 3. occupancy
        fig, ax = plt.subplots(2, 2, figsize=(11, 9))
        lm = np.where(W, lam[THRS[0]] * 1e3, np.nan)
        mod_map(ax[0, 0], np.where(lm > 0, lm, np.nan), f"pulse rate per module, A >= {THRS[0]:.0f} ADC (per ms), onset {FID[0]:.0f}-{FID[1]:.0f} ns", 1e-2, max(np.nanmax(lm), 1.0), log=True, label="per ms")
        pos = lm[W & (lm > 0)]
        ax[0, 1].hist(lam[THRS[0]][W] * 1e3, bins=np.concatenate([[0], np.logspace(-2, np.log10(max(pos.max(), 1.1)), 30)]), color=cols[0]); ax[0, 1].set_xscale("symlog", linthresh=1e-2); ax[0, 1].set_yscale("log")
        ax[0, 1].set_xlabel(f"rate per module (per ms), A >= {THRS[0]:.0f}"); ax[0, 1].set_ylabel("modules")
        x = np.arange(len(band_lbl))
        for j, thr in enumerate(THRS):
            rb = [lam[thr][band == b + 1].sum() / max((band == b + 1).sum(), 1) * 1e3 for b in range(len(band_lbl))]
            ax[1, 0].bar(x + (j - 0.5) * 0.35, rb, 0.35, label=f"A >= {thr:.0f} ADC", color=cols[j])
        ax[1, 0].set_xticks(x); ax[1, 0].set_xticklabels(band_lbl, rotation=30); ax[1, 0].set_ylabel("mean rate per module (per ms)"); ax[1, 0].legend(fontsize=7); ax[1, 0].set_title("by Chebyshev ring around the beam hole")
        ax[1, 1].axis("off")
        lines = []
        for thr, name, bl_, mv, pr_, nk in res_rows:
            if thr == THRS[0] and name == "170ns":
                lines.append(f"{bl_:9s} kept {nk:8d}  measured (24590, other pulse within +-85 ns, A>={thr:.0f}) {100 * mv:6.3f} %   predicted random {100 * pr_:8.5f} %")
        ax[1, 1].text(0, 1, "Predicted random pile-up vs measured, by ring\n(window 170 ns, A_other >= %.0f ADC)\n\n" % THRS[0] + "\n".join(lines), va="top", family="monospace", fontsize=6)
        fig.suptitle("Occupancy of the random-trigger run"); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        # pile-up bar chart
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        for pi, name in enumerate(("170ns", "any")):
            for j, thr in enumerate(THRS):
                mv = [r[3] for r in res_rows if r[0] == thr and r[1] == name]
                pv = [r[4] for r in res_rows if r[0] == thr and r[1] == name]
                ax[pi].bar(x + (j - 0.5) * 0.35, [100 * v for v in mv], 0.35, color=cols[j], label=f"measured 24590, A_other >= {thr:.0f}")
                ax[pi].plot(x + (j - 0.5) * 0.35, [100 * v for v in pv], "k_", ms=14, mew=2, label="predicted random" if j == 0 else None)
            ax[pi].set_yscale("log"); ax[pi].set_xticks(x); ax[pi].set_xticklabels(band_lbl, rotation=30); ax[pi].set_ylabel("fraction of kept channels (%)")
            ax[pi].set_title("another pulse within +-85 ns of the main" if name == "170ns" else "another pulse anywhere in the record"); ax[pi].legend(fontsize=6)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        # amplitude spectrum
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        ax[0].errorbar(ac, rate_a, np.sqrt(np.maximum(spec, 1)) / expo, fmt="o-", color=cols[0], ms=3)
        ax[0].set_xscale("log"); ax[0].set_yscale("log"); ax[0].set_xlabel("fitted amplitude (ADC)"); ax[0].set_ylabel("PbWO4 array rate per us per bin"); ax[0].set_title("amplitude spectrum of random pulses (32 log bins)")
        ax[1].errorbar(ac, rate_a / np.diff(ae), np.sqrt(np.maximum(spec, 1)) / expo / np.diff(ae), fmt="o-", color=cols[0], ms=3)
        ax[1].set_xscale("log"); ax[1].set_yscale("log"); ax[1].set_xlabel("fitted amplitude (ADC)"); ax[1].set_ylabel("dN/dA per us per ADC"); ax[1].set_title("per ADC")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        # residual profiles
        fig, ax = plt.subplots(2, 2, figsize=(11, 8), sharex=True)
        for ab in range(NAB):
            a_ = ax[ab // 2, ab % 2]
            for j, (tag, n, sx, sxx, sv) in enumerate((("24578 random", Gn[ab], Gsx[ab], Gsxx[ab], Gsv[ab]), ("24590 production", Pn[ab], Psx[ab], Psxx[ab], Psv[ab]))):
                mean, sem, rms, asg = stats(n, sx, sxx, sv)
                okt = n > 100
                a_.errorbar(TAU[okt], 100 * mean[okt], 100 * sem[okt], fmt="o-", ms=2.5, lw=0.8, color=cols[j], label=f"{tag} (max bin n = {int(n.max())})")
                if j == 0:
                    a_.fill_between(TAU[okt], -100 * asg[okt] / np.sqrt(n[okt]), 100 * asg[okt] / np.sqrt(n[okt]), color="0.85", label="+-assumed sigma / sqrt(n)")
            a_.axhline(0, color="k", lw=0.5); a_.set_title(f"A {amp_lbl[ab]} ADC"); a_.set_ylabel("mean r / A (%)"); a_.legend(fontsize=6)
            if ab >= 2:
                a_.set_xlabel("tau = t - t0 (ns)")
        fig.suptitle("Isolated single pulses: mean fit residual relative to the 0p73GeV_H2 template"); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        # rms / assumed
        fig, ax = plt.subplots(1, 2, figsize=(11, 4))
        for ab in range(NAB):
            for j, (tag, n, sx, sxx, sv) in enumerate((("24578 random", Gn[ab], Gsx[ab], Gsxx[ab], Gsv[ab]), ("24590 production", Pn[ab], Psx[ab], Psxx[ab], Psv[ab]))):
                mean, sem, rms, asg = stats(n, sx, sxx, sv)
                okt = n > 100
                ax[j].plot(TAU[okt], (rms / asg)[okt], "o-", ms=2.5, lw=0.8, label=f"A {amp_lbl[ab]}")
                ax[j].set_title(f"{tag}: rms of r/A over assumed sigma(t)"); ax[j].set_xlabel("tau (ns)"); ax[j].axhline(1, color="k", lw=0.5); ax[j].legend(fontsize=6)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        # examples
        def draw(pdf, lst, title, nr, nc):
            fig, axs = plt.subplots(nr, nc, figsize=(11, 8.5))
            for ax_, (pre, j) in zip(axs.ravel(), lst):
                Y = ex[f"{pre}_Y"][j].astype(float)
                ped = float(ex[f"{pre}_ped"][j])
                comp = ex[f"{pre}_comp"][j]
                kk = int(ex[f"{pre}_k"][j])
                ax_.plot(TT, Y, "o", ms=1.8, color="k")
                if kk > 0:
                    ax_.plot(TT, ped + comp.sum(axis=0), "-", color=cols[1], lw=1)
                    for ci in range(kk):
                        ax_.plot(TT, ped + comp[ci], "--", lw=0.6, color=cols[2])
                ax_.axhline(ped, color="0.5", lw=0.5)
                ax_.set_title(f"{names[int(ex[f'{pre}_mi'][j])]} ev {int(ex[f'{pre}_evn'][j])} k={kk} chi2 {ex[f'{pre}_chi2'][j]:.1f}\n" + ("A " + ",".join(f"{x:.0f}" for x in ex[f"{pre}_a"][j][:kk]) + " t " + ",".join(f"{x:.0f}" for x in ex[f"{pre}_t"][j][:kk]) if kk else f"no pulse; sd30 {Y[:30].std():.2f}"), fontsize=5.5)
                ax_.tick_params(labelsize=5)
            for ax_ in axs.ravel()[len(lst):]:
                ax_.axis("off")
            fig.suptitle(title, fontsize=9); fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        draw(pdf, k1_idx, "Random-trigger waveforms with exactly 1 fitted pulse (4 per amplitude bin 10-20, 20-50, 50-200, 200-1000, >1000 ADC; black samples, orange fit)", 5, 4)
        draw(pdf, k2_idx, "Random-trigger waveforms with >= 2 fitted pulses", 2, 5)
        draw(pdf, k0_idx, "Random-trigger waveforms without a fitted pulse (stored, k == 0)", 2, 5)
    tsv = Path(str(stem) + ".tsv")
    with open(tsv, "w") as fh:
        fh.write("section\titem\tkey\tvalue\terr\tn\tnote\n")
        for r in rows:
            fh.write("\t".join(f"{v:.6g}" if isinstance(v, (float, np.floating)) else str(v) for v in r) + "\n")
    print(f"[done] {pdf_path} {tsv}; {len(rows)} rows; {time.monotonic() - t_start:.0f} s", flush=True)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    s = sub.add_parser("scan")
    s.add_argument("out"); s.add_argument("run", type=int); s.add_argument("evio", nargs="+"); s.add_argument("--max-events", type=int, default=0)
    s.add_argument("--ref-from", default=None, help="scan npz whose mod_sig_ref / mod_b_ref (median over files, >= 20 waveforms) replace the file's own noise reference")
    r = sub.add_parser("report")
    r.add_argument("scan_dir"); r.add_argument("stem")
    r.add_argument("--ref-sigma", default=str(HERE.parent.parent / "output" / "pileup_rates" / "scan" / "run_24590.npz"))
    a = ap.parse_args()
    if a.cmd == "scan":
        scan(a.out, a.evio, a.run, a.max_events, a.ref_from)
    else:
        report(a.scan_dir, a.stem, a.ref_sigma)


if __name__ == "__main__":
    main()
