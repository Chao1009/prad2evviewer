#!/usr/bin/env python3
"""pileup_artefact_study.py -- the +12-40 ns secondary excess of the pile-up scan.

Plan 20261024 task B (Checkpoint 2 decisions).  The 119-split scan (spline banks,
kappa 1.0, E_dom > 50 MeV, r_model 0.10) shows an excess of accepted secondaries at
dt = +12..40 ns in e-p / C runs (the h50 two-tau bank puts it at +36..76 ns), i.e.
a template-dependent artefact.  This script characterises it.

Subcommands (outputs under --out, default output/pileup_20261024/artefact/):

  scan     one task per (run, split) (--jobs N): walks the EVIO split ONCE, runs the
           mf with the h50 bank (kappa 1.010) and the spline bank (kappa 1.0) on the
           SAME waveforms (scan settings: nsigma 5, min_amp 5, kmax 4, r_model 0.10,
           trigger 0x100, unclipped, E_dom > 50 MeV, early / edge exclusions as
           plot_pileup_examples.Scanner) and writes
             scan/sec_<run>_<split>.npz  EVERY accepted non-dominant pulse of every
                 channel with a dominant above threshold (both banks): dt, A and E
                 ratios, dchi2 (plot_pileup_examples.fit_terms), joint dominant, and
                 the dominant-alone fit (MFReconstructor kmax=1) amplitude / onset /
                 energy; plus per-module / per-amplitude-bin denominators
             scan/prof_<run>_<split>.npz  residual (data - single fit)/A profiles vs
                 time since onset (1 ns bins), per bank x population x amplitude bin
                 x fallback x radius, and per hot channel
             scan/raw_<run>_<split>.npz   raw waveforms + fit parameters of the
                 window-secondary channels, of the hot channels and of a random
                 clean sample (for the average residual of single channels)
  analyze  tables / histograms / profiles / correlations / energy effect from the
           scan outputs (items 1, 2, 3, 5), summary.json
  synth    item 4 on the pseudo-data cells (output/pileup_eff_20261024/<slug>/synth,
           spline bank, r_model 0.10): same dchi2 per accepted secondary + truth
  rules    item 4: acceptance rules evaluated on the recorded secondaries and on
           the synth rows; writes rules_*.tsv and the figures
  all      scan, synth, analyze, rules

  MPLBACKEND=Agg python3 pileup_artefact_study.py all --jobs 6
  MPLBACKEND=Agg python3 pileup_artefact_study.py scan --splits 24590:10 --max-events 3000 --out /tmp/x
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import os
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402

import _common as C                                   # noqa: E402
from prad2py import dec                               # noqa: E402
import mf_pulse                                       # noqa: E402
import plot_pileup_examples as ppe                    # noqa: E402
from benchmark_hycal_timing_scan import (             # noqa: E402
    _DATABASE_DIR, _DEFAULT_EVIO_DIR, _HC_CALIB_MAP)

HERE = Path(__file__).resolve().parent
OUT_DEFAULT = HERE / "output" / "pileup_20261024" / "artefact"
SPL_PAT = str(HERE / "output" / "templates_spline" / "{runtype}" / "templates_{runtype}.npz")
SYNTH_ROOT = HERE / "output" / "pileup_eff_20261024"
SPLITS_DEFAULT = ["25320:50", "25320:100", "24503:10", "24503:20", "24590:10", "24590:20"]
BANKS = ("s", "h")                       # spline, h50
BANK_NAME = {"s": "spline", "h": "h50"}
KAPPA = {"s": 1.0, "h": 1.010}
THR_MEV = 50.0
R_MODEL = 0.10
HOT = ("W335", "W1", "W1019", "W1020", "W934")
AMP_EDGES = (150.0, 400.0, 1000.0, 2500.0)          # ADC; bins 0..4 (4: > 2500)
AMP_LABELS = ("<150", "150-400", "400-1000", "1000-2500", ">2500")
R_EDGES = (200.0, 350.0)                            # mm; radius bins 0..2
R_LABELS = ("R<200", "200-350", "R>350")
POPS = ("clean", "weak", "win", "other", "honly")
# clean: one spline pulse and no h50 secondary above r_model*E_dom; weak: spline
# secondaries all below r_model*E_dom (not reported); win: a reported spline
# secondary in WIN; other: reported, none in WIN; honly: one spline pulse, a reported
# h50 secondary
WIN = (12.0, 40.0)
PROF_LO, PROF_HI = -20, 240                         # ns, 1 ns bins
NPB = PROF_HI - PROF_LO
RAW_CLEAN_PER_SPLIT = 2000
RAW_HOT_CAP = 4000
SAME_NS = 8.0            # the dominant-alone onset must lie this close to the joint one
SEC_COLS = ("event", "mod", "fb", "npk", "k", "n_rep",
            "E_dom", "A_dom", "t0_dom", "time_dom", "ped", "sig",
            "t0", "time", "dt", "A", "E", "E_ratio", "A_ratio", "dchi2",
            "chi2_full", "dof", "chi2dof_mf", "edge", "reported",
            "A1", "t01", "E1", "chi1", "has1", "same1", "npk_other")
FB_NAMES = ("fb0", "fb1", "fb2")


def parse_splits(spec):
    out = []
    for s in spec:
        r, sp = s.split(":")
        out.append((int(r), int(sp)))
    return out


def run_slug(run):
    from gate_scan_20261022 import runtype_slug
    return runtype_slug(run)


def amp_bin(a):
    return np.digitize(a, AMP_EDGES)


def rad_bin(r):
    return np.digitize(r, R_EDGES)


def bank_fallback(slug, names):
    """Spline fallback status (0 own, 1 relaxed, 2 run-type average) per module name."""
    z = np.load(SPL_PAT.format(runtype=slug), allow_pickle=False)
    d = {str(n): int(f) for n, f in zip(z["names"], z["fallback"])}
    return np.array([d.get(n, 2) for n in names], dtype=np.int8)


def bank_widths(bank, nm):
    """Per module: t_rise50, t_fall50 (ns after onset), FWHM, t_peak of the unit template."""
    t = np.arange(0.0, 160.0, 0.1)
    out = np.zeros((nm, 4))
    for m in range(nm):
        s = bank.shape(m, t, dtype=np.float64)
        ip = int(np.argmax(s))
        up = np.nonzero(s[:ip + 1] >= 0.5)[0]
        dn = np.nonzero(s[ip:] >= 0.5)[0]
        tr = t[up[0]] if up.size else np.nan
        tf = t[ip + dn[-1]] if dn.size else np.nan
        out[m] = (tr, tf, tf - tr, t[ip])
    return out


# ------------------------------------------------------------------ scan ----
class Study(ppe.Scanner):
    """ppe.Scanner's EVIO collection; flush() runs both banks and records everything."""

    def __init__(self, args, p, banks, recs, recs1, ebeam, slug):
        super().__init__(args, p, banks["s"], recs["s"], KAPPA["s"], ebeam)
        self.banks, self.recs, self.recs1 = banks, recs, recs1
        self.slug = slug
        nm = len(self.mods)
        self.fb_m = bank_fallback(slug, [m.name for m in self.mods])
        xy = np.array([[m.x, m.y] for m in self.mods], dtype=np.float64)
        self.rad_m = np.hypot(xy[:, 0], xy[:, 1])
        self.rb_m = rad_bin(self.rad_m)
        self.hot_idx = {n: i for i, m in enumerate(self.mods)
                        for n in HOT if m.name == n}
        self.hot_of_m = np.full(nm, -1, dtype=np.int64)
        for j, n in enumerate(HOT):
            if n in self.hot_idx:
                self.hot_of_m[self.hot_idx[n]] = j
        self.den = {b: np.zeros((nm, 5), dtype=np.int64) for b in BANKS}   # [mod, A bin]
        self.cnt = {b: dict(n_early=0, n_edge=0) for b in BANKS}
        self.rows = []
        self.meta = []
        self.sec = {b: {c: [] for c in SEC_COLS} for b in BANKS}
        self.sec_split = {b: [] for b in BANKS}
        ng = len(POPS) * 5 * 3 * 3
        self.prof = {b: np.zeros((ng, NPB, 3)) for b in BANKS}            # sum, sumsq, n
        self.prof_hot = {b: np.zeros((len(POPS), len(HOT), NPB, 3)) for b in BANKS}
        self.n_pop = np.zeros((len(POPS), 5), dtype=np.int64)             # [pop, A bin]
        self.raw = []
        self.rng = np.random.default_rng(12345)
        self.n_clean_seen = 0
        self.hot_raw_n = {}
        self._reset()

    # -- blocks ----------------------------------------------------------
    def flush(self):
        if not self.rows:
            return
        a = self.a
        Y = np.stack(self.rows)
        mid = np.array([m[3] for m in self.meta], dtype=np.int64)
        evs = np.array([m[1] for m in self.meta], dtype=np.int64)
        spl = np.array([m[0] for m in self.meta], dtype=np.int64)
        n, ns = Y.shape
        ovf = self.overflow
        res = {b: self.recs[b].run(Y, mid, ovf) for b in BANKS}
        clipped = res["s"].clipped
        if not (res["s"].p_ch.size and res["h"].p_ch.size):
            self._reset()
            return
        self.n_clip += int(clipped.sum())
        t_edge = mf_pulse.DT * (ns - self.rec.L - 2)
        t_early = mf_pulse.DT * self.rec.ped_n / 2.0
        info = {}
        for b in BANKS:
            r = res[b]
            pc = r.p_ch
            E = (mf_pulse.pulse_energy_units(self.banks[b], mid[pc], r.p_t0, r.p_amp,
                                             KAPPA[b]) * self.cal[mid[pc]])
            starts = np.searchsorted(pc, np.arange(n + 1))
            emax = np.zeros(n)
            np.maximum.at(emax, pc, E)
            tmin = np.full(n, np.inf)
            np.minimum.at(tmin, pc, r.p_t0)
            dom = np.full(n, -1, dtype=np.int64)
            isd = E == emax[pc]
            dom[pc[isd][::-1]] = np.nonzero(isd)[0][::-1]       # first of ties
            high = (emax > self.thr) & ~clipped
            early = high & (tmin < t_early)
            edge = high & ~early & (r.p_t0[np.maximum(dom, 0)] >= t_edge)
            ok = high & ~early & ~edge
            self.cnt[b]["n_early"] += int(early.sum())
            self.cnt[b]["n_edge"] += int(edge.sum())
            np.add.at(self.den[b], (mid[ok], amp_bin(r.p_amp[dom[ok]])), 1)
            # pulses other than the dominant above the secondary rule (any E_dom)
            rep = np.zeros(n, dtype=bool)
            sel = (E > R_MODEL * emax[pc]) & (np.arange(pc.size) != dom[pc]) & (r.p_t0 < t_edge)
            rep[pc[sel]] = True
            info[b] = dict(E=E, starts=starts, dom=dom, ok=ok, emax=emax, rep=rep)
        okany = info["s"]["ok"] | info["h"]["ok"]
        R = np.nonzero(okany)[0]
        res1 = {}
        a1 = {b: np.full(n, np.nan) for b in BANKS}
        t1 = {b: np.full(n, np.nan) for b in BANKS}
        c1 = {b: np.full(n, np.nan) for b in BANKS}
        e1 = {b: np.full(n, np.nan) for b in BANKS}
        if R.size:
            for b in BANKS:
                r1 = self.recs1[b].run(Y[R], mid[R], ovf)
                pcl = r1.p_ch
                if pcl.size:
                    g = R[pcl]
                    a1[b][g], t1[b][g], c1[b][g] = r1.p_amp, r1.p_t0, r1.p_chi2
                    e1[b][g] = (mf_pulse.pulse_energy_units(
                        self.banks[b], mid[g], r1.p_t0, r1.p_amp, KAPPA[b])
                        * self.cal[mid[g]])
        # ---- per channel: secondaries, populations ----------------------
        npk_s, npk_h = res["s"].npk, res["h"].npk
        label = np.full(n, -1, dtype=np.int64)       # population index (spline ok rows)
        fitinfo = {b: {} for b in BANKS}
        for b in BANKS:
            r, inf = res[b], info[b]
            bank = self.banks[b]
            for ch in np.nonzero(inf["ok"] & (r.npk >= 2))[0]:
                lo, hi = inf["starts"][ch], inf["starts"][ch + 1]
                kd = int(inf["dom"][ch] - lo)
                y_raw = Y[ch]
                ped, sig = float(r.ped[ch]), float(r.sigma[ch])
                y = y_raw.astype(np.float64) - ped
                w = np.ones(ns)
                t0s = r.p_t0[lo:hi].astype(np.float64)
                amps = r.p_amp[lo:hi].astype(np.float64)
                _, _, chi2, dof, dchi2 = ppe.fit_terms(y, w, sig, bank, int(mid[ch]),
                                                       t0s, amps)
                Er = inf["E"][lo:hi]
                tim = r.p_time[lo:hi].astype(np.float64)
                fitinfo[b][ch] = (kd, t0s, amps, Er, tim, chi2, dof, dchi2, ped, sig)
        for ch in np.nonzero(info["s"]["ok"])[0]:
            if ch in fitinfo["s"]:
                kd, t0s, amps, Er, tim, *_ = fitinfo["s"][ch]
                rep = [k for k in range(len(t0s)) if k != kd
                       and Er[k] / Er[kd] > a.secondary_frac and t0s[k] < t_edge]
                if rep:
                    dts = [tim[k] - tim[kd] for k in rep]
                    win = any(WIN[0] <= d < WIN[1] for d in dts)
                    label[ch] = POPS.index("win" if win else "other")
                else:
                    label[ch] = POPS.index("weak")
            else:
                label[ch] = POPS.index("honly" if info["h"]["rep"][ch] else "clean")
        # ---- secondary rows ---------------------------------------------
        for b in BANKS:
            r, inf = res[b], info[b]
            S = self.sec[b]
            for ch, (kd, t0s, amps, Er, tim, chi2, dof, dchi2, ped, sig) in fitinfo[b].items():
                reps = [k for k in range(len(t0s)) if k != kd
                        and Er[k] / Er[kd] > self.a.secondary_frac and t0s[k] < t_edge]
                m = int(mid[ch])
                has1 = np.isfinite(a1[b][ch])
                same = bool(has1 and abs(t1[b][ch] - t0s[kd]) < SAME_NS)
                for k in range(len(t0s)):
                    if k == kd:
                        continue
                    row = dict(
                        event=evs[ch], mod=m, fb=self.fb_m[m], npk=len(t0s), k=k,
                        n_rep=len(reps), E_dom=Er[kd], A_dom=amps[kd],
                        t0_dom=t0s[kd], time_dom=tim[kd], ped=ped, sig=sig,
                        t0=t0s[k], time=tim[k], dt=tim[k] - tim[kd], A=amps[k],
                        E=Er[k], E_ratio=Er[k] / Er[kd], A_ratio=amps[k] / amps[kd],
                        dchi2=dchi2[k], chi2_full=chi2, dof=dof,
                        chi2dof_mf=float(r.p_chi2[inf["starts"][ch]]),
                        edge=int(t0s[k] >= t_edge), reported=int(k in reps),
                        A1=a1[b][ch], t01=t1[b][ch], E1=e1[b][ch], chi1=c1[b][ch],
                        has1=int(has1), same1=int(same),
                        npk_other=int((npk_h if b == "s" else npk_s)[ch]))
                    for c in SEC_COLS:
                        S[c].append(row[c])
                    self.sec_split[b].append(spl[ch])
        # ---- residual profiles + raw store ------------------------------
        dom_s = info["s"]["dom"]
        t0d = np.where(dom_s >= 0, res["s"].p_t0[np.maximum(dom_s, 0)], np.nan)
        good = np.nonzero((label >= 0) & np.isfinite(a1["s"]) & np.isfinite(a1["h"])
                          & (np.abs(t1["s"] - t0d) < SAME_NS))[0]
        if good.size:
            self._accumulate(Y, mid, label, good, res, a1, t1)
        self._store_raw(Y, mid, evs, spl, label, good, res, info, a1, t1)
        self._reset()

    def _accumulate(self, Y, mid, label, good, res, a1, t1):
        ns = Y.shape[1]
        ts = mf_pulse.DT * np.arange(ns)
        As = a1["s"][good]
        ab = amp_bin(As)
        m = mid[good]
        gid = ((label[good] * 5 + ab) * 3 + self.fb_m[m]) * 3 + self.rb_m[m]
        np.add.at(self.n_pop, (label[good], ab), 1)
        hot = self.hot_of_m[m]
        for b in BANKS:
            A1, T1 = a1[b][good], t1[b][good]
            y = Y[good].astype(np.float64) - res[b].ped[good][:, None]
            tau = ts[None, :] - T1[:, None]
            mdl = A1[:, None] * self.banks[b].shape(m[:, None], tau, dtype=np.float64)
            rr = (y - mdl) / A1[:, None]
            ib = np.floor(tau - PROF_LO).astype(np.int64)
            sel = (ib >= 0) & (ib < NPB)
            g = np.broadcast_to(gid[:, None], tau.shape)[sel]
            bi = ib[sel]
            v = rr[sel]
            flat = g * NPB + bi
            P = self.prof[b].reshape(-1, 3)
            np.add.at(P[:, 0], flat, v)
            np.add.at(P[:, 1], flat, v * v)
            np.add.at(P[:, 2], flat, 1.0)
            hh = hot >= 0
            if hh.any():
                pp = np.broadcast_to((label[good] * len(HOT) + np.maximum(hot, 0))[:, None],
                                     tau.shape)
                selh = sel & hh[:, None]
                fl = pp[selh] * NPB + ib[selh]
                PH = self.prof_hot[b].reshape(-1, 3)
                vh = rr[selh]
                np.add.at(PH[:, 0], fl, vh)
                np.add.at(PH[:, 1], fl, vh * vh)
                np.add.at(PH[:, 2], fl, 1.0)

    def _store_raw(self, Y, mid, evs, spl, label, good, res, info, a1, t1):
        ns = Y.shape[1]
        gs = set(good.tolist())
        keep = []
        for ch in good:
            lab = label[ch]
            hot = self.hot_of_m[mid[ch]]
            if lab == POPS.index("win"):
                keep.append(ch)
            elif hot >= 0:
                k = self.hot_raw_n.get(int(hot), 0)
                if k < RAW_HOT_CAP:
                    self.hot_raw_n[int(hot)] = k + 1
                    keep.append(ch)
            elif lab == POPS.index("clean"):
                self.n_clean_seen += 1
                if self.rng.random() < RAW_CLEAN_PER_SPLIT / 900000.0:
                    keep.append(ch)
        for ch in keep:
            d = dict(raw=Y[ch].astype(np.uint16), mod=mid[ch], event=evs[ch],
                     label=label[ch])
            for b in BANKS:
                r, inf = res[b], info[b]
                lo, hi = inf["starts"][ch], inf["starts"][ch + 1]
                t0 = np.zeros(4)
                am = np.zeros(4)
                t0[:hi - lo] = r.p_t0[lo:hi]
                am[:hi - lo] = r.p_amp[lo:hi]
                d[f"{b}_t0"], d[f"{b}_amp"], d[f"{b}_ped"] = t0, am, float(r.ped[ch])
                d[f"{b}_sig"] = float(r.sigma[ch])
                d[f"{b}_A1"], d[f"{b}_t01"] = a1[b][ch], t1[b][ch]
                d[f"{b}_kd"] = int(inf["dom"][ch] - lo)
            self.raw.append(d)


def scan_task(job):
    a = SimpleNamespace(**job["a"])
    run, split, fpath = job["run"], job["split"], job["path"]
    sdir = Path(a.out) / "scan"
    t_start = time.monotonic()
    calib = _HC_CALIB_MAP.get(run)
    calib_path = str(_DATABASE_DIR / "calibration" / calib) if calib else ""
    p = C.setup_pipeline(evio_path=fpath, run_num=run, hc_calib_file=calib_path)
    ebeam = float(p.geo.Ebeam)
    slug = run_slug(run)
    mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    wcfg = p.wave_ana.cfg
    wdict = dict(smooth_order=int(wcfg.smooth_order),
                 int_tail_ratio=float(wcfg.int_tail_ratio),
                 tail_break_n=int(wcfg.tail_break_n))
    names, pwo4 = [m.name for m in mods], [bool(m.is_pwo4()) for m in mods]
    paths = {"s": SPL_PAT.format(runtype=slug)}
    hp, hcol, _ = ppe.resolve_templates("auto", run)
    paths["h"] = hp
    banks, infos = {}, {}
    for b, col in (("s", "sampled"), ("h", hcol)):
        banks[b], infos[b] = mf_pulse.build_bank(names, pwo4, paths[b], wdict, column=col)
    ped_n = int(wcfg.ped_nsamples)
    recs = {b: mf_pulse.MFReconstructor(banks[b], nsigma=5.0, min_amp=5.0, kmax=4,
                                        r_model=R_MODEL, ped_nsamples=ped_n)
            for b in BANKS}
    recs1 = {b: mf_pulse.MFReconstructor(banks[b], nsigma=5.0, min_amp=5.0, kmax=1,
                                         r_model=R_MODEL, ped_nsamples=ped_n)
             for b in BANKS}
    print(f"[setup] run {run} split {split:05d} {slug} E_beam {ebeam:.0f}  spline "
          f"{os.path.basename(paths['s'])} ({infos['s']['n_fallback']} fb2)  h50 "
          f"{os.path.basename(paths['h'])}:{hcol}", flush=True)
    sa = SimpleNamespace(dominant_min_mev=THR_MEV, dominant_frac=0.5,
                         secondary_frac=R_MODEL, include_clipped=False)
    sc = Study(sa, p, banks, recs, recs1, ebeam, slug)
    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    if ch.open_auto(str(fpath)) != dec.Status.success:
        raise RuntimeError(f"cannot open {fpath}")
    n_ev = n_blk = 0
    stop = False
    while not stop and ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            fe = d["event"]
            if int(fe.info.trigger_bits) != C.PHYSICS_TRIGGER_BITS:
                continue
            sc.collect(split, fe)
            n_ev += 1
            n_blk += 1
            if n_blk >= ppe.BLOCK:
                sc.flush()
                n_blk = 0
                if n_ev % (ppe.BLOCK * 40) == 0:
                    print(f"  run {run} split {split}: {n_ev} events, "
                          f"{time.monotonic() - t_start:.0f} s", flush=True)
            if a.max_events and n_ev >= a.max_events:
                stop = True
                break
    sc.flush()
    os.makedirs(sdir, exist_ok=True)
    tag = f"{run}_{split:05d}"
    cols = {}
    for b in BANKS:
        for c in SEC_COLS:
            cols[f"{b}_{c}"] = np.asarray(sc.sec[b][c])
        cols[f"{b}_den"] = sc.den[b]
        cols[f"{b}_n_early"] = sc.cnt[b]["n_early"]
        cols[f"{b}_n_edge"] = sc.cnt[b]["n_edge"]
        cols[f"{b}_sha1"] = np.array(infos[b]["sha1"])
    meta = dict(run=run, split=split, slug=slug, ebeam=ebeam, n_events=n_ev,
                n_clip=sc.n_clip, thr_mev=THR_MEV, r_model=R_MODEL, kappa=KAPPA,
                templates={b: paths[b] for b in BANKS}, hot=sc.hot_idx,
                t_edge=mf_pulse.DT * (sc.ns - sc.rec.L - 2),
                t_early=mf_pulse.DT * sc.rec.ped_n / 2.0, ns=int(sc.ns),
                elapsed=time.monotonic() - t_start, n_clean_seen=sc.n_clean_seen)
    np.savez_compressed(
        sdir / f"sec_{tag}.npz", meta=np.array(json.dumps(meta)),
        mod_name=np.array(names), fb_m=sc.fb_m, mod_x=np.array([m.x for m in mods]),
        mod_y=np.array([m.y for m in mods]), **cols)
    np.savez_compressed(
        sdir / f"prof_{tag}.npz", meta=np.array(json.dumps(meta)), prof_s=sc.prof["s"],
        prof_h=sc.prof["h"], hot_s=sc.prof_hot["s"], hot_h=sc.prof_hot["h"],
        n_pop=sc.n_pop, pops=np.array(POPS), hot=np.array(HOT))
    rw = {}
    if sc.raw:
        for k in sc.raw[0]:
            rw[k] = np.array([d[k] for d in sc.raw])
    np.savez_compressed(sdir / f"raw_{tag}.npz", meta=np.array(json.dumps(meta)), **rw)
    el = time.monotonic() - t_start
    ns_ = {b: int(sc.den[b].sum()) for b in BANKS}
    print(f"[done] run {run} split {split:05d}: {n_ev} events, den spline {ns_['s']} "
          f"h50 {ns_['h']}, secondaries spline {len(sc.sec['s']['dt'])} h50 "
          f"{len(sc.sec['h']['dt'])}, raw {len(sc.raw)}  ({el:.0f} s)", flush=True)
    return dict(run=run, split=split, n_events=n_ev, elapsed=el)


def cmd_scan(a):
    out = Path(a.out)
    (out / "scan").mkdir(parents=True, exist_ok=True)
    jobs = []
    for run, split in parse_splits(a.splits):
        fp = dict(ppe.local_splits(a.evio_dir, run)).get(split)
        if fp is None:
            sys.exit(f"[ERROR] no EVIO file for run {run} split {split}")
        if not a.force and (out / "scan" / f"raw_{run}_{split:05d}.npz").exists():
            print(f"[skip] {run} {split:05d} done (--force to redo)")
            continue
        jobs.append(dict(run=run, split=split, path=str(fp),
                         a=dict(out=str(out), max_events=a.max_events)))
    if a.jobs <= 1:
        for j in jobs:
            scan_task(j)
    else:
        with cf.ProcessPoolExecutor(max_workers=a.jobs) as ex:
            for r in ex.map(scan_task, jobs):
                pass


# --------------------------------------------------------------- analysis ---
import pandas as pd                                    # noqa: E402

RUNTYPE_LABEL = {"0p73GeV_H2": "0.73 GeV H2", "2p24GeV_H2": "2.24 GeV H2",
                 "3p49GeV_H2": "3.49 GeV H2"}
BANK_COLOR = {"s": "C0", "h": "C3"}
POP_COLOR = {"clean": "k", "weak": "C2", "win": "C0", "other": "C1", "honly": "C3"}


def load_scan(out):
    """Concatenate the scan outputs: ({bank: DataFrame}, den table, meta list, extras)."""
    sdir = Path(out) / "scan"
    files = sorted(sdir.glob("sec_*.npz"))
    if not files:
        sys.exit(f"[ERROR] no sec_*.npz under {sdir}; run the scan first")
    dfs = {b: [] for b in BANKS}
    dens, metas = [], []
    names = None
    for f in files:
        z = np.load(f, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        metas.append(meta)
        nm = z["mod_name"]
        names = nm
        for b in BANKS:
            d = {c: z[f"{b}_{c}"] for c in SEC_COLS}
            df = pd.DataFrame(d)
            df["run"], df["split"], df["slug"] = meta["run"], meta["split"], meta["slug"]
            df["name"] = nm[df["mod"].to_numpy()] if len(df) else []
            dfs[b].append(df)
        dens.append(dict(run=meta["run"], split=meta["split"], slug=meta["slug"],
                         fb=z["fb_m"], x=z["mod_x"], y=z["mod_y"],
                         den={b: z[f"{b}_den"] for b in BANKS}))
    out_df = {}
    for b in BANKS:
        df = pd.concat(dfs[b], ignore_index=True)
        mx = dens[0]["x"]
        my = dens[0]["y"]
        df["x"], df["y"] = mx[df["mod"].to_numpy()], my[df["mod"].to_numpy()]
        df["R"] = np.hypot(df["x"], df["y"])
        df["ab"] = amp_bin(df["A_dom"].to_numpy())
        df["ch"] = (df["run"].astype(np.int64) * 100000 + df["split"]) * 10**11 \
            + df["event"].astype(np.int64) * 2000 + df["mod"]
        df["inwin"] = (df["dt"] >= WIN[0]) & (df["dt"] < WIN[1])
        out_df[b] = df
    return out_df, dens, metas, names


def den_summary(dens, b, slug, by="fb"):
    """Denominator channel counts for a run type by fallback status (0/1/2) and A bin."""
    n = {f: 0 for f in (0, 1, 2)}
    nab = np.zeros(5, dtype=np.int64)
    for d in dens:
        if d["slug"] != slug:
            continue
        den = d["den"][b]
        for f in (0, 1, 2):
            n[f] += int(den[d["fb"] == f].sum())
        nab += den.sum(axis=0)
    return n, nab


def flat_expect(dt, lo_side=((-140.0, -20.0), (60.0, 180.0)), win=WIN):
    """Flat-in-dt expectation of the window count from the sidebands."""
    side = sum(((dt >= a) & (dt < b)).sum() for a, b in lo_side)
    width = sum(b - a for a, b in lo_side)
    return side / width * (win[1] - win[0])


def cmd_analyze(a):
    out = Path(a.out)
    dfs, dens, metas, names = load_scan(out)
    adir = out
    adir.mkdir(parents=True, exist_ok=True)
    slugs = sorted({m["slug"] for m in metas})
    for b in BANKS:                      # events with >= NB_BURST candidate channels
        d_ = dfs[b]
        rp = d_[(d_["reported"] == 1) & (d_["edge"] == 0)]
        d_["burst"] = False
        d_.loc[rp.index, "burst"] = burst_events(rp, NB_BURST)
    summ = dict(window_ns=list(WIN), burst_nb=NB_BURST, thr_mev=THR_MEV, r_model=R_MODEL, splits=[
        [m["run"], m["split"], m["n_events"]] for m in metas])
    # ------------------------------------------------ (1) tables + histograms
    for b in BANKS:
        df = dfs[b]
        cols = ["run", "split", "event", "name", "fb", "npk", "n_rep", "E_dom",
                "A_dom", "time_dom", "dt", "E_ratio", "A_ratio", "dchi2", "chi2_full",
                "dof", "reported", "edge", "E1", "A1", "same1", "R"]
        df[cols].to_csv(adir / f"secondaries_{BANK_NAME[b]}.tsv.gz", sep="\t",
                        index=False, float_format="%.4g", compression="gzip")
    rows = []
    for slug in slugs:
        for b in BANKS:
            df = dfs[b]
            dn, dab = den_summary(dens, b, slug)
            for veto in (False, True):
                for fbv in (None, 0, 1, 2):
                    sel = (df["slug"] == slug) & (df["reported"] == 1) & (df["edge"] == 0)
                    if veto:
                        sel &= ~df["burst"]
                    den_n = sum(dn.values()) if fbv is None else dn[fbv]
                    if fbv is not None:
                        sel &= df["fb"] == fbv
                    s = df[sel]
                    chs = s.drop_duplicates("ch")
                    chw = s[s["inwin"]].drop_duplicates("ch")
                    cho = s.groupby("ch")["inwin"].agg(lambda v: (~v).any()).sum() if len(s) else 0
                    n_w = int(s["inwin"].sum())
                    fe = flat_expect(s["dt"].to_numpy()) if len(s) else 0.0
                    rows.append(dict(
                        runtype=slug, bank=BANK_NAME[b], burst_veto=veto,
                        fb="all" if fbv is None else fbv,
                        den=den_n, n_sec=len(s), n_cand=len(chs), n_cand_win=len(chw),
                        n_cand_nonwin=int(cho), rate_pct=100.0 * len(chs) / max(den_n, 1),
                        rate_nonwin_pct=100.0 * cho / max(den_n, 1),
                        n_sec_win=n_w, win_frac=n_w / max(len(s), 1),
                        flat_expect_win=fe, excess_win=n_w - fe,
                        excess_frac=(n_w - fe) / max(len(s), 1)))
    tab = pd.DataFrame(rows)
    tab.to_csv(adir / "summary_by_runtype_bank_fb.tsv", sep="\t", index=False,
               float_format="%.4g")
    summ["by_runtype_bank_fb"] = tab.to_dict("records")
    # dt / ratio histograms
    dt_edges = np.arange(-160.0, 244.0, 4.0)
    r_edges = np.logspace(math.log10(0.02), math.log10(5.0), 33)
    hrows = []
    fig, axes = plt.subplots(len(slugs) * 3, 2, figsize=(11, 3.0 * 3 * len(slugs)),
                             squeeze=False)
    for i, slug in enumerate(slugs):
        for f in (0, 1, 2):
            ax_dt, ax_r = axes[3 * i + f]
            for b in BANKS:
                df = dfs[b]
                sel = ((df["slug"] == slug) & (df["fb"] == f) & (df["edge"] == 0))
                for tag, sub in (("reported", sel & (df["reported"] == 1)),
                                 ("all", sel)):
                    h, _ = np.histogram(df.loc[sub, "dt"], dt_edges)
                    hr, _ = np.histogram(df.loc[sub, "E_ratio"], r_edges)
                    for j in range(len(h)):
                        hrows.append((slug, BANK_NAME[b], f, tag, "dt", dt_edges[j],
                                      dt_edges[j + 1], int(h[j])))
                    for j in range(len(hr)):
                        hrows.append((slug, BANK_NAME[b], f, tag, "E_ratio", r_edges[j],
                                      r_edges[j + 1], int(hr[j])))
                    if tag == "reported":
                        ax_dt.stairs(h, dt_edges, color=BANK_COLOR[b],
                                     label=f"{BANK_NAME[b]} (n={int(h.sum())})")
                        ax_r.stairs(hr, r_edges, color=BANK_COLOR[b])
            ax_dt.axvspan(*WIN, color="0.9", zorder=0)
            ax_dt.set_title(f"{RUNTYPE_LABEL.get(slug, slug)}  fb{f}: reported secondaries "
                            f"(E/E_dom > {R_MODEL})", fontsize=9)
            ax_dt.set_xlabel("dt = t_sec - t_dom (ns)")
            ax_dt.set_ylabel("secondaries / 4 ns")
            ax_dt.legend(fontsize=7)
            ax_r.set_xscale("log")
            ax_r.set_xlabel("E_sec / E_dom")
    fig.tight_layout()
    fig.savefig(adir / "hist_dt_ratio_by_runtype_fb.png", dpi=110)
    plt.close(fig)
    pd.DataFrame(hrows, columns=["runtype", "bank", "fb", "sel", "var", "lo", "hi",
                                 "count"]).to_csv(adir / "hist_dt_ratio.tsv", sep="\t",
                                                  index=False, float_format="%.5g")
    # -------------------------------------------------- (2) residual profiles
    analyze_profiles(out, slugs, summ, dens)
    # -------------------------------------------------- (3) fallback + hot
    analyze_hot(out, dfs, dens, slugs, summ, names)
    # -------------------------------------------------- (5) energy effect
    analyze_energy(out, dfs, dens, slugs, summ)
    analyze_crossbank(out, dfs, summ)
    analyze_window_shape(out, dfs, dens, slugs, summ, names)
    update_summary(out, analyze=summ)
    print(f"[analyze] wrote {adir}/summary.json")


def eb(ax, x, y, yerr, fmt=None, lw=1.0, color=None, ls="-", label=None):
    """Line with a +-1 sigma band."""
    ln, = ax.plot(x, y, ls=ls, lw=lw, color=color, label=label)
    ax.fill_between(x, y - yerr, y + yerr, color=ln.get_color(), alpha=0.2, lw=0)


def update_summary(out, **kw):
    p = Path(out) / "summary.json"
    d = {}
    if p.exists():
        d = json.loads(p.read_text())
    d.update(kw)
    p.write_text(json.dumps(d, indent=1, default=float))


TAU_WIN = {"pre[-20,-4)": (-20, -4), "[0,12)": (0, 12), "[12,40)": (12, 40),
           "[40,80)": (40, 80), "[80,160)": (80, 160), "[160,240)": (160, 240)}


def _msn(P):
    """sums [..., 3] -> mean, sem, rms over the last-but-one axis already summed."""
    n = P[..., 2]
    with np.errstate(invalid="ignore", divide="ignore"):
        mean = P[..., 0] / n
        var = np.maximum(P[..., 1] / n - mean ** 2, 0.0)
        return mean, np.sqrt(var / n), np.sqrt(var), n


def load_prof(out):
    prof = {}
    for f in sorted((Path(out) / "scan").glob("prof_*.npz")):
        z = np.load(f, allow_pickle=False)
        meta = json.loads(str(z["meta"]))
        sl = meta["slug"]
        d = prof.setdefault(sl, dict(n=0, s=0.0, h=0.0, hs=0.0, hh=0.0, n_pop=0))
        d["s"] = d["s"] + z["prof_s"].reshape(len(POPS), 5, 3, 3, NPB, 3)
        d["h"] = d["h"] + z["prof_h"].reshape(len(POPS), 5, 3, 3, NPB, 3)
        d["hs"] = d["hs"] + z["hot_s"]
        d["hh"] = d["hh"] + z["hot_h"]
        d["n_pop"] = d["n_pop"] + z["n_pop"]
        d["n"] += 1
    return prof


def win_mean(P, lo, hi):
    """Sum a profile array [..., NPB, 3] over tau in [lo, hi) and return mean, sem."""
    i0, i1 = int(lo - PROF_LO), int(hi - PROF_LO)
    return _msn(P[..., i0:i1, :].sum(axis=-2))[:2]


def analyze_profiles(out, slugs, summ, dens):
    out = Path(out)
    prof = load_prof(out)
    tau = np.arange(PROF_LO, PROF_HI) + 0.5
    rows, wrows = [], []
    res = {}
    for slug in slugs:
        d = prof[slug]
        fig, axs = plt.subplots(3, 4, figsize=(17, 10), sharex=True)
        for ab in range(4):
            for row, (pops, ylim) in enumerate((
                    ((("s", "clean"), ("h", "clean"), ("s", "weak")), 2.0),
                    ((("s", "win"), ("s", "other"), ("h", "honly")), 25.0),
                    ((("s", "win"), ("h", "win")), 25.0))):
                ax = axs[row, ab]
                for b, pop in pops:
                    P = d[b][POPS.index(pop), ab].sum(axis=(0, 1))     # [NPB, 3]
                    nb = P[:, 2]
                    k = 4
                    Pk = P[: (NPB // k) * k].reshape(-1, k, 3).sum(axis=1)
                    mean, sem, _, n = _msn(Pk)
                    ok = n >= 5
                    eb(ax, tau[: (NPB // k) * k].reshape(-1, k).mean(1)[ok],
                       100 * mean[ok], 100 * sem[ok], lw=1.1,
                       color=POP_COLOR[pop], ls="-" if b == "s" else "--",
                       label=f"{BANK_NAME[b]} {pop} (n={int(nb.max())})")
                ax.axvspan(*WIN, color="0.92", zorder=0)
                ax.axhline(0, color="0.4", lw=0.6)
                ax.set_ylim(-ylim, ylim)
                ax.legend(fontsize=7)
                if row == 0:
                    ax.set_title(f"{RUNTYPE_LABEL.get(slug, slug)}  A {AMP_LABELS[ab]} ADC")
                if ab == 0:
                    ax.set_ylabel("mean (data - fit)/A  [% of peak]")
                if row == 2:
                    ax.set_xlabel("time since onset (ns), bank's own onset convention")
        fig.tight_layout()
        fig.savefig(out / f"resid_profile_{slug}.png", dpi=100)
        plt.close(fig)
        # tables
        for b in BANKS:
            for pi, pop in enumerate(POPS):
                for ab in range(5):
                    P = d[b][pi, ab].sum(axis=(0, 1))
                    mean, sem, rms, n = _msn(P)
                    for j in range(NPB):
                        if n[j] > 0:
                            rows.append((slug, BANK_NAME[b], pop, AMP_LABELS[ab],
                                         PROF_LO + j, 100 * mean[j], 100 * sem[j],
                                         100 * rms[j], int(n[j])))
        # windowed means by amp / fb / radius
        for b in BANKS:
            for pi, pop in enumerate(POPS):
                for ab in range(5):
                    for fbv in (-1, 0, 1, 2):
                        for rb in (-1, 0, 1, 2):
                            P = d[b][pi, ab]
                            P = P[:, :, :, :] if False else P
                            Pf = P if fbv < 0 else P[fbv:fbv + 1]
                            Pf = Pf.sum(axis=0) if fbv < 0 else Pf[0]
                            Pr = Pf.sum(axis=0) if rb < 0 else Pf[rb]
                            if Pr[..., 2].max() < 30:
                                continue
                            rec = [slug, BANK_NAME[b], pop, AMP_LABELS[ab],
                                   "all" if fbv < 0 else fbv,
                                   "all" if rb < 0 else R_LABELS[rb],
                                   int(Pr[..., 2].max())]
                            for k, (lo, hi) in TAU_WIN.items():
                                m, s = win_mean(Pr, lo, hi)
                                rec += [100 * float(m), 100 * float(s)]
                            wrows.append(rec)
        res[slug] = {pop: {int(n): int(d["n_pop"][i, ab]) for ab, n in enumerate(range(5))}
                     for i, pop in enumerate(POPS)}
    pd.DataFrame(rows, columns=["runtype", "bank", "pop", "A_bin", "tau_ns",
                                "mean_pct", "sem_pct", "rms_pct", "n"]).to_csv(
        out / "resid_profiles.tsv.gz", sep="\t", index=False, float_format="%.4g",
        compression="gzip")
    cols = ["runtype", "bank", "pop", "A_bin", "fb", "radius", "n_max"]
    for k in TAU_WIN:
        cols += [f"mean_{k}_pct", f"sem_{k}_pct"]
    wtab = pd.DataFrame(wrows, columns=cols)
    wtab.to_csv(out / "resid_window_means.tsv", sep="\t", index=False, float_format="%.4g")
    # largest 4 ns-binned |mean| residual of the clean population in [0, 120) ns and the
    # window / tail means by fallback status and radius
    bump = {}
    for slug in slugs:
        d = prof[slug]
        for b in BANKS:
            for ab in range(4):
                P = d[b][POPS.index("clean"), ab].sum(axis=(0, 1))
                Pk = P[: (NPB // 4) * 4].reshape(-1, 4, 3).sum(axis=1)
                mean, sem, _, n = _msn(Pk)
                x = np.arange(PROF_LO, PROF_LO + (NPB // 4) * 4, 4) + 2.0
                ok = (x >= 0) & (x < 120) & (n >= 50)
                j = int(np.argmax(np.where(ok, np.abs(mean), -1)))
                bump[f"{slug}|{BANK_NAME[b]}|A{AMP_LABELS[ab]}"] = dict(
                    max_abs_mean_pct=100 * float(abs(mean[j])), at_tau_ns=float(x[j]),
                    sem_pct=100 * float(sem[j]), mean_signed_pct=100 * float(mean[j]))
    summ["clean_max_abs_bin_mean"] = bump
    fbr = wtab[(wtab["pop"] == "clean") & (wtab["bank"] == "spline")
               & (wtab["A_bin"].isin(["<150", "150-400", "400-1000", "1000-2500"]))]
    summ["clean_spline_window_means_by_fb_radius"] = fbr[
        ["runtype", "A_bin", "fb", "radius", "n_max", "mean_[12,40)_pct", "sem_[12,40)_pct",
         "mean_[40,80)_pct", "mean_[80,160)_pct"]].to_dict("records")
    summ["population_counts"] = res
    # headline: clean-population window means per A bin and bank (all fb / radius)
    head = wtab[(wtab["pop"].isin(["clean", "win", "other"])) & (wtab["fb"] == "all")
                & (wtab["radius"] == "all")]
    summ["resid_window_means_headline"] = head.to_dict("records")


def _q(v, ps=(16, 50, 84)):
    v = np.asarray(v, dtype=float)
    return [float(np.percentile(v, p)) if len(v) else float("nan") for p in ps]


def hot_scan119(out_scan="output/pileup_20261024/scan"):
    """Hot-channel candidate statistics over ALL splits of the 119-split scan."""
    sdir = HERE / out_scan
    rows = []
    den_run = {}
    cands = []
    for f in sorted(sdir.glob("counts_*.npz")):
        z = np.load(f, allow_pickle=False)
        run = int(z["run"])
        names = z["mod_name"]
        d = den_run.setdefault(run, np.zeros(len(names), dtype=np.int64))
        d += z["den_m"]
        den_run["names"] = names
    for f in sorted(sdir.glob("cands_*.tsv")):
        df = pd.read_csv(f, sep="\t", comment="#")
        if len(df):
            cands.append(df)
    c = pd.concat(cands, ignore_index=True)
    names = list(den_run["names"])
    # first-secondary quantities (the candidate's first reported secondary)
    c["dt1"] = c["dt_ns"].astype(str).str.split(";").str[0].astype(float)
    c["r1"] = c["E_ratio"].astype(str).str.split(";").str[0].astype(float)
    c["dchi1"] = c["dchi2"].astype(str).str.split(";").str[0].astype(float)
    c["inwin"] = (c["dt1"] >= WIN[0]) & (c["dt1"] < WIN[1])
    return c, den_run, names


def analyze_hot(out, dfs, dens, slugs, summ, names):
    out = Path(out)
    df = dfs["s"]
    rep = df[(df["reported"] == 1) & (df["edge"] == 0)]
    # ---- fallback / radius / amplitude dependence of the window rate -----
    rows = []
    for slug in slugs:
        dd = [d for d in dens if d["slug"] == slug]
        R = np.hypot(dd[0]["x"], dd[0]["y"])
        rbm = rad_bin(R)
        fbm = dd[0]["fb"]
        denmat = sum(d["den"]["s"] for d in dd)                 # [mod, ab]
        s = rep[rep["slug"] == slug]
        sw = s[s["inwin"]].drop_duplicates("ch")
        so = s.groupby("ch").filter(lambda g: (~g["inwin"]).any()).drop_duplicates("ch")
        for fbv in (-1, 0, 1, 2):
            for rb in (-1, 0, 1, 2):
                for ab in (-1, 0, 1, 2, 3):
                    mm = np.ones(len(fbm), bool)
                    if fbv >= 0:
                        mm &= fbm == fbv
                    if rb >= 0:
                        mm &= rbm == rb
                    den = int(denmat[mm][:, slice(None) if ab < 0 else ab].sum())
                    if den < 200:
                        continue
                    def sel(t, fbv=fbv, rb=rb, ab=ab):
                        m_ = np.ones(len(t), dtype=bool)
                        if fbv >= 0:
                            m_ &= (t["fb"] == fbv).to_numpy()
                        if rb >= 0:
                            m_ &= rad_bin(t["R"].to_numpy()) == rb
                        if ab >= 0:
                            m_ &= (t["ab"] == ab).to_numpy()
                        return m_
                    nw = int(sel(sw).sum()) if len(sw) else 0
                    no = int(sel(so).sum()) if len(so) else 0
                    rows.append(dict(runtype=slug, fb="all" if fbv < 0 else fbv,
                                     radius="all" if rb < 0 else R_LABELS[rb],
                                     A_bin="all" if ab < 0 else AMP_LABELS[ab], den=den,
                                     n_win_cand=nw, n_nonwin_cand=no,
                                     rate_win_pct=100 * nw / den,
                                     rate_nonwin_pct=100 * no / den))
    t = pd.DataFrame(rows)
    t.to_csv(out / "window_rate_by_fb_radius_amp.tsv", sep="\t", index=False,
             float_format="%.4g")
    summ["window_rate_by_fb_radius_amp"] = t[
        (t["A_bin"] == "all")].to_dict("records")
    # ---- hot channels, my splits ---------------------------------------
    hrows = []
    for slug in slugs:
        dd = [d for d in dens if d["slug"] == slug]
        nm = list(names)
        denm = sum(d["den"]["s"] for d in dd).sum(axis=1)
        for h in HOT:
            if h not in nm:
                continue
            mi = nm.index(h)
            s = rep[(rep["slug"] == slug) & (rep["mod"] == mi)]
            if denm[mi] == 0:
                continue
            chs = s.drop_duplicates("ch")
            hrows.append(dict(runtype=slug, module=h, fb=int(dd[0]["fb"][mi]),
                              den=int(denm[mi]), n_cand=len(chs),
                              rate_pct=100 * len(chs) / denm[mi],
                              n_sec=len(s), n_win=int(s["inwin"].sum()),
                              dt_q16_50_84=_q(s["dt"]), Eratio_q16_50_84=_q(s["E_ratio"]),
                              A_dom_q16_50_84=_q(s["A_dom"])))
    pd.DataFrame(hrows).to_csv(out / "hot_channels_splits.tsv", sep="\t", index=False,
                               float_format="%.4g")
    summ["hot_channels_splits"] = hrows
    # share of the window candidates carried by the hot channels / by fb
    tot = {}
    for slug in slugs:
        s = rep[(rep["slug"] == slug) & rep["inwin"]].drop_duplicates("ch")
        hotm = [list(names).index(h) for h in HOT if h in list(names)]
        tot[slug] = dict(n_win_cand=len(s), n_hot=int(s["mod"].isin(hotm).sum()),
                         top10_modules_share=float(
                             s["mod"].value_counts().head(10).sum() / max(len(s), 1)),
                         n_modules=int(s["mod"].nunique()))
    summ["window_candidates_hot_share"] = tot
    # ---- hot channels over ALL splits of the 119-split scan ---------------
    c, den_run, nm119 = hot_scan119()
    rr = []
    for run in sorted(k for k in den_run if k != "names"):
        for h in HOT:
            if h not in nm119:
                continue
            mi = nm119.index(h)
            s = c[(c["run"] == run) & (c["module"] == h)]
            den = int(den_run[run][mi])
            if den == 0 and s.empty:
                continue
            rr.append(dict(run=run, module=h, den=den, n_cand=len(s),
                           rate_pct=100 * len(s) / max(den, 1),
                           dt_q16_50_84=_q(s["dt1"]), Eratio_q16_50_84=_q(s["r1"]),
                           A_dom_q16_50_84=_q(s["A_dominant_ADC"]),
                           dchi2_med=float(np.median(s["dchi1"])) if len(s) else np.nan,
                           frac_in_win=float(s["inwin"].mean()) if len(s) else np.nan))
    pd.DataFrame(rr).to_csv(out / "hot_channels_scan119.tsv", sep="\t", index=False,
                            float_format="%.4g")
    summ["hot_channels_scan119"] = rr
    # contribution of the hot channels to the window excess, all 119 splits
    cw = c[c["inwin"]]
    hs = {}
    for run in sorted(c["run"].unique()):
        a_ = cw[cw["run"] == run]
        n_all = int((c["run"] == run).sum())
        hs[int(run)] = dict(n_cand=n_all, n_win=len(a_),
                            n_win_hot=int(a_["module"].isin(HOT).sum()),
                            n_cand_hot=int(c[(c["run"] == run)]["module"].isin(HOT).sum()))
    summ["scan119_window_hot_share"] = hs
    # ---- W934 --------------------------------------------------------------
    w = c[c["module"] == "W934"]
    if len(w):
        ab = np.digitize(w["A_dominant_ADC"], AMP_EDGES)
        wtab = []
        for run in sorted(w["run"].unique()):
            for k in range(5):
                s = w[(w["run"] == run) & (ab == k)]
                if len(s) >= 5:
                    wtab.append(dict(run=int(run), A_bin=AMP_LABELS[k], n=len(s),
                                     dt_q16_50_84=_q(s["dt1"]), Eratio_q16_50_84=_q(s["r1"]),
                                     dchi2_med=float(np.median(s["dchi1"]))))
        pd.DataFrame(wtab).to_csv(out / "W934_dt_ratio_by_run_amp.tsv", sep="\t",
                                  index=False, float_format="%.4g")
        summ["W934"] = dict(
            n=len(w), dt_q16_50_84=_q(w["dt1"]), Eratio_q16_50_84=_q(w["r1"]),
            n_by_run={int(r): int((w["run"] == r).sum()) for r in w["run"].unique()},
            by_run_amp=wtab)
        fig, axs = plt.subplots(1, 3, figsize=(16, 4.4))
        for run in sorted(w["run"].unique()):
            s = w[w["run"] == run]
            axs[0].plot(s["A_dominant_ADC"], s["dt1"], ".", ms=3, label=f"{run} (n={len(s)})")
            axs[1].plot(s["A_dominant_ADC"], s["r1"], ".", ms=3)
        axs[0].set_xscale("log")
        axs[1].set_xscale("log")
        axs[0].set_xlabel("A_dom (ADC)")
        axs[0].set_ylabel("dt (ns)")
        axs[0].legend(fontsize=7)
        axs[1].set_xlabel("A_dom (ADC)")
        axs[1].set_ylabel("E_sec / E_dom")
        axs[0].set_title("W934: first secondary of every candidate (119-split scan)")
        prof = load_prof(out)
        tau = np.arange(PROF_LO, PROF_HI) + 0.5
        for slug, d in prof.items():
            if "W934" not in HOT:
                break
            for b in BANKS:
                P = d["hs" if b == "s" else "hh"][:, HOT.index("W934")].sum(axis=0)
                k = 4
                Pk = P[: (NPB // k) * k].reshape(-1, k, 3).sum(axis=1)
                mean, sem, _, n = _msn(Pk)
                ok = n >= 5
                if ok.sum():
                    eb(axs[2],tau[: (NPB // k) * k].reshape(-1, k).mean(1)[ok],
                                    100 * mean[ok], 100 * sem[ok], fmt="none", lw=0.9,
                                    color=BANK_COLOR[b], ls="-" if b == "s" else "--",
                                    label=f"{RUNTYPE_LABEL.get(slug, slug)} {BANK_NAME[b]}")
        axs[2].axvspan(*WIN, color="0.92", zorder=0)
        axs[2].axhline(0, color="0.4", lw=0.6)
        axs[2].set_xlabel("time since onset (ns)")
        axs[2].set_ylabel("mean (data - single fit)/A  [% of peak]")
        axs[2].set_title("W934 average residual, all pulses (my splits)")
        axs[2].legend(fontsize=7)
        fig.tight_layout()
        fig.savefig(out / "hot_W934.png", dpi=110)
        plt.close(fig)
    # ---- average residual of every hot channel (my splits) ----------------
    prof = load_prof(out)
    fig, axs = plt.subplots(len(HOT), 1, figsize=(9, 2.6 * len(HOT)), sharex=True)
    hr = {}
    for j, h in enumerate(HOT):
        ax = axs[j]
        for slug, d in prof.items():
            for b in BANKS:
                P = d["hs" if b == "s" else "hh"][:, j].sum(axis=0)
                if P[:, 2].max() < 30:
                    continue
                k = 4
                Pk = P[: (NPB // k) * k].reshape(-1, k, 3).sum(axis=1)
                mean, sem, _, n = _msn(Pk)
                ok = n >= 5
                eb(ax,tau[: (NPB // k) * k].reshape(-1, k).mean(1)[ok], 100 * mean[ok],
                            100 * sem[ok], fmt="none", lw=0.9, ls="-" if b == "s" else "--",
                            label=f"{RUNTYPE_LABEL.get(slug, slug)} {BANK_NAME[b]} "
                                  f"(n={int(P[:, 2].max())})")
                hr[f"{h}|{slug}|{BANK_NAME[b]}"] = dict(
                    n=int(P[:, 2].max()),
                    win=[float(x) for x in
                         (100 * win_mean(P, *WIN)[0], 100 * win_mean(P, *WIN)[1])])
        ax.axvspan(*WIN, color="0.92", zorder=0)
        ax.axhline(0, color="0.4", lw=0.6)
        ax.set_ylabel(f"{h}\nresid [% of peak]")
        ax.legend(fontsize=6, ncol=2)
    axs[-1].set_xlabel("time since onset (ns)")
    fig.tight_layout()
    fig.savefig(out / "hot_channels_residual.png", dpi=100)
    plt.close(fig)
    summ["hot_residual_window_mean_pct"] = hr


def analyze_energy(out, dfs, dens, slugs, summ):
    out = Path(out)
    rows = []
    res = {}
    fig, axs = plt.subplots(len(slugs), 4, figsize=(15, 3.2 * len(slugs)), squeeze=False)
    for i, slug in enumerate(slugs):
        for b in BANKS:
            df = dfs[b]
            dd = dens
            dn, dab = den_summary(dens, b, slug)
            den = sum(dn.values())
            s = df[(df["slug"] == slug) & (df["reported"] == 1) & (df["edge"] == 0)]
            for tag, sel in (("window", s["inwin"]), ("window_isolated", s["inwin"] & ~s["burst"]),
                             ("window_burst", s["inwin"] & s["burst"]), ("nonwindow", ~s["inwin"]),
                             ("any", s["inwin"] | ~s["inwin"])):
                t = s[sel]
                if tag == "nonwindow":      # channels with a reported non-window secondary
                    t = t[~t["ch"].isin(s[s["inwin"]]["ch"])]
                if tag == "window":
                    pass
                t = t.drop_duplicates("ch")
                aff = len(t)
                tt = t[(t["has1"] == 1) & (t["same1"] == 1)]
                dE = (tt["E_dom"] - tt["E1"]) / tt["E1"]
                for ab in range(-1, 5):
                    m = np.ones(len(tt), bool) if ab < 0 else (tt["ab"] == ab).to_numpy()
                    if m.sum() < 5:
                        continue
                    x = dE.to_numpy()[m]
                    nden = den if ab < 0 else int(dab[ab])
                    rows.append(dict(
                        runtype=slug, bank=BANK_NAME[b], sample=tag,
                        A_bin="all" if ab < 0 else AMP_LABELS[ab], n=int(m.sum()),
                        den=nden, frac_of_den_pct=100 * int(m.sum()) / max(nden, 1),
                        dE_over_E1_median_pct=100 * float(np.median(x)),
                        dE_q16_pct=100 * float(np.percentile(x, 16)),
                        dE_q84_pct=100 * float(np.percentile(x, 84)),
                        dE_mean_pct=100 * float(x.mean()),
                        frac_shift_gt1pct=float(np.mean(np.abs(x) > 0.01)),
                        frac_shift_gt5pct=float(np.mean(np.abs(x) > 0.05)),
                        sum_over_E1_median_pct=100 * float(np.median(
                            ((tt["E_dom"] + tt["E"]) / tt["E1"] - 1).to_numpy()[m]))
                        if tag != "any" else np.nan))
            if b == "s":
                t = s[s["inwin"]].drop_duplicates("ch")
                tt = t[(t["has1"] == 1) & (t["same1"] == 1)]
                for ab in range(4):
                    x = (100 * (tt["E_dom"] - tt["E1"]) / tt["E1"])[tt["ab"] == ab]
                    ax = axs[i, ab]
                    if len(x):
                        ax.hist(np.clip(x, -40, 15), bins=np.linspace(-40, 15, 56),
                                color="C0")
                    ax.axvline(0, color="k", lw=0.6)
                    ax.set_title(f"{RUNTYPE_LABEL.get(slug, slug)} A {AMP_LABELS[ab]} "
                                 f"(n={len(x)})", fontsize=9)
                    ax.set_xlabel("(E_dom joint - E_dom alone)/E_alone  [%]")
        res[slug] = dict(den=den)
    fig.tight_layout()
    fig.savefig(out / "energy_effect_window_secondary.png", dpi=100)
    plt.close(fig)
    t = pd.DataFrame(rows)
    t.to_csv(out / "energy_effect.tsv", sep="\t", index=False, float_format="%.4g")
    summ["energy_effect"] = t[(t["A_bin"] == "all")].to_dict("records")
    summ["energy_effect_by_A"] = t[t["sample"].isin(["window", "window_isolated"])
                                   & (t["bank"] == "spline")].to_dict("records")


# ------------------------------------------------------------ synth stage ---
SYNTH_SLUGS = ("0p73GeV_H2", "2p24GeV_C", "2p24GeV_H2", "2p24GeV_Ta_X17", "3p49GeV_H2")
TOL_NS = 4.0
NB_BURST = 10          # an event with >= this many candidate channels is a burst event


def synth_task(job):
    import benchmark_pileup_mf as bm
    slug = job["slug"]
    sdir = Path(job["synth"])
    tpath = SPL_PAT.format(runtype=slug)
    cells, man = bm.find_cells(sdir)
    first = np.load(cells[0][2], allow_pickle=False)
    run_num = int(first["base_run_numbers"][0])
    p = C.setup_pipeline(evio_path="", run_num=run_num,
                         daq_config=bm._manifest_cli_value(man, "--daq-config"),
                         hc_map_file=bm._manifest_cli_value(man, "--hc-map-file"))
    hy = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    names = [m.name for m in hy]
    ix = {n: i for i, n in enumerate(names)}
    wcfg = p.wave_ana.cfg
    bank, info = mf_pulse.build_bank(
        names, [bool(m.is_pwo4()) for m in hy], tpath,
        dict(smooth_order=int(wcfg.smooth_order), int_tail_ratio=float(wcfg.int_tail_ratio),
             tail_break_n=int(wcfg.tail_break_n)), column="sampled")
    rec = mf_pulse.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4, r_model=R_MODEL,
                                   ped_nsamples=int(wcfg.ped_nsamples))
    t_early = mf_pulse.DT * rec.ped_n / 2.0
    rows = []
    for dt_c, ratio_c, pth in cells:
        d = np.load(pth, allow_pickle=False)
        wf = d["waveforms"]
        n_ev, ns = wf.shape
        mid = np.array([ix[n] for n in d["base_channel_names"]], dtype=np.int64)
        tru = d["truth_onset_sampled_ns"].astype(np.float64)       # [n, 2] inj, base
        res = rec.run(wf, mid)
        pc = res.p_ch
        starts = np.searchsorted(pc, np.arange(n_ev + 1))
        t_edge = mf_pulse.DT * (ns - rec.L - 2)
        E = (mf_pulse.pulse_energy_units(bank, mid[pc], res.p_t0, res.p_amp, 1.0)
             if pc.size else np.zeros(0))
        for ev in range(n_ev):
            lo, hi = starts[ev], starts[ev + 1]
            if hi - lo == 0:
                rows.append(dict(slug=slug, dt_cell=dt_c, ratio_cell=ratio_c, ev=ev,
                                 npk=0))
                continue
            Er = E[lo:hi]
            kd = int(np.argmax(Er))
            t0s = res.p_t0[lo:hi].astype(np.float64)
            amps = res.p_amp[lo:hi].astype(np.float64)
            tim = res.p_time[lo:hi].astype(np.float64)
            ped, sig = float(res.ped[ev]), float(res.sigma[ev])
            if hi - lo >= 2:
                y = wf[ev].astype(np.float64) - ped
                _, _, chi2, dof, dch = ppe.fit_terms(y, np.ones(ns), sig, bank, int(mid[ev]),
                                                     t0s, amps)
            else:
                chi2 = dof = np.nan
                dch = [np.nan]
            base_ok = float(tru[ev, 1]) >= t_early
            for k in range(hi - lo):
                ni = abs(t0s[k] - tru[ev, 0]) <= TOL_NS
                nb = abs(t0s[k] - tru[ev, 1]) <= TOL_NS
                rows.append(dict(
                    slug=slug, dt_cell=dt_c, ratio_cell=ratio_c, ev=ev, mod=int(mid[ev]),
                    npk=hi - lo, k=k, dom=int(k == kd), t0=t0s[k], time=tim[k],
                    A=amps[k], E=Er[k], E_ratio=Er[k] / Er[kd], A_ratio=amps[k] / amps[kd],
                    dt=tim[k] - tim[kd], dchi2=dch[k] if k != kd or hi - lo == 1 else np.nan,
                    chi2_full=chi2, dof=dof, A_dom=amps[kd], sig=sig,
                    edge=int(t0s[k] >= t_edge), near_inj=int(ni), near_base=int(nb),
                    base_clean=int(base_ok), fb=int(d["inject_template_fallback"][ev]),
                    t_inj=tru[ev, 0], t_base=tru[ev, 1]))
    df = pd.DataFrame(rows)
    return slug, df, bank_widths(bank, len(names)), names


def cmd_synth(a):
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)
    jobs = [dict(slug=s, synth=str(SYNTH_ROOT / s / "synth")) for s in SYNTH_SLUGS
            if (SYNTH_ROOT / s / "synth").exists()]
    allp = []
    widths = {}
    if a.jobs <= 1:
        results = [synth_task(j) for j in jobs]
    else:
        with cf.ProcessPoolExecutor(max_workers=min(a.jobs, len(jobs))) as ex:
            results = list(ex.map(synth_task, jobs))
    for slug, df, wd, names in results:
        allp.append(df)
        widths[slug] = wd
        print(f"[synth] {slug}: {len(df)} pulse rows", flush=True)
    df = pd.concat(allp, ignore_index=True)
    df.to_csv(out / "synth_pulses.tsv.gz", sep="\t", index=False, float_format="%.5g",
              compression="gzip")
    np.savez_compressed(out / "bank_widths.npz", **{k: v for k, v in widths.items()},
                        names=np.array(results[0][3]))
    print(f"[synth] wrote {out}/synth_pulses.tsv.gz ({len(df)} rows)")


# ------------------------------------------------------------- rules stage ---
def rule_specs():
    """[(name, spec)] acceptance rules; a secondary is kept when rule_keep() is True."""
    R = [("none", dict(kind="none"))]
    for T in (25, 50, 100, 200, 500, 1000, 3000):
        R.append((f"dchi2>={T}", dict(kind="dchi2", T=T)))
    for T in (2, 5, 10, 20, 50, 100):
        R.append((f"dchi2/dof>={T}", dict(kind="dchi2_dof", T=T)))
    for f in (0.1, 0.3, 0.5, 0.7):
        R.append((f"dchi2/chi2_wo>={f}", dict(kind="dchi2_rel", f=f)))
    for rb in (0.17, 0.19, 0.25, 0.29, 0.49):
        for D in (24, 32, 40, 48, 60):
            R.append((f"dt>={D}|r>={rb}", dict(kind="dtabs", D=D, rb=rb)))
        for k in (0.75, 1.0, 1.25, 1.5):
            R.append((f"dt>={k}FWHM|r>={rb}", dict(kind="dtfwhm", k=k, rb=rb)))
        for k in (0.5, 0.75, 1.0):
            R.append((f"dt>={k}tfall50|r>={rb}", dict(kind="dtfall", k=k, rb=rb)))
    for D, rb, T in ((40, 0.19, 1000), (40, 0.25, 1000), (40, 0.29, 1000), (40, 0.29, 3000)):
        R.append((f"dt>={D}|r>={rb}|dchi2>={T}", dict(kind="combo", D=D, rb=rb, T=T)))
    return R


def rule_keep(d, spec):
    """d: DataFrame/dict of arrays dt, E_ratio, dchi2, chi2_full, dof, fwhm, tfall."""
    k = spec["kind"]
    n = len(d["dt"])
    if k == "none":
        return np.ones(n, dtype=bool)
    if k == "dchi2":
        return d["dchi2"] >= spec["T"]
    if k == "dchi2_dof":
        return d["dchi2"] / np.maximum(d["dof"], 1) >= spec["T"]
    if k == "dchi2_rel":
        return d["dchi2"] / np.maximum(d["chi2_full"] + d["dchi2"], 1e-9) >= spec["f"]
    pos = (d["dt"] >= 0)
    if k == "dtabs":
        rej = pos & (d["dt"] < spec["D"])
    elif k == "dtfwhm":
        rej = pos & (d["dt"] < spec["k"] * d["fwhm"])
    elif k == "dtfall":
        rej = pos & (d["dt"] < spec["k"] * d["tfall"])
    elif k == "combo":
        rej = pos & (d["dt"] < spec["D"]) & (d["dchi2"] < spec["T"])
    else:
        raise ValueError(k)
    return ~(rej & (d["E_ratio"] < spec["rb"]))


def _ch_rate(keep, chs, den):
    """Channel-level candidate count: channels with >= 1 kept secondary."""
    if den <= 0 or len(chs) == 0:
        return 0, 0.0
    u, inv = np.unique(chs, return_inverse=True)
    n = int((np.bincount(inv, weights=keep.astype(float)) > 0).sum())
    return n, 100.0 * n / den


def burst_events(rep, nb):
    """Boolean mask of the secondaries that belong to an event with >= nb candidate channels."""
    key = rep["run"].astype(np.int64) * 10**4 * 10**9 + rep["split"].astype(np.int64) * 10**9 + rep["event"].astype(np.int64)
    ch = rep.drop_duplicates("ch").assign(_k=lambda t: t["run"].astype(np.int64) * 10**13
                                          + t["split"].astype(np.int64) * 10**10 + t["event"].astype(np.int64))
    cnt = ch.groupby("_k").size()
    kk = rep["run"].astype(np.int64) * 10**13 + rep["split"].astype(np.int64) * 10**10 + rep["event"].astype(np.int64)
    return kk.map(cnt).fillna(0).to_numpy() >= nb


def cmd_rules(a):
    out = Path(a.out)
    dfs, dens, metas, names = load_scan(out)
    wz = np.load(out / "bank_widths.npz", allow_pickle=False)
    wnames = [str(n) for n in wz["names"]]
    assert wnames == [str(n) for n in names], "bank_widths.npz geometry differs from the scan"
    df = dfs["s"]
    rep = df[(df["reported"] == 1) & (df["edge"] == 0)].copy()
    fw = np.zeros(len(rep))
    tf = np.zeros(len(rep))
    for slug in rep["slug"].unique():
        m = (rep["slug"] == slug).to_numpy()
        w = wz[slug]
        fw[m] = w[rep.loc[m, "mod"].to_numpy(), 2]
        tf[m] = w[rep.loc[m, "mod"].to_numpy(), 1]
    rep["fwhm"], rep["tfall"] = fw, tf
    summ = dict(bank_widths_median={s: dict(t_rise50=float(np.nanmedian(wz[s][:, 0])),
                                           t_fall50=float(np.nanmedian(wz[s][:, 1])),
                                           fwhm=float(np.nanmedian(wz[s][:, 2])),
                                           t_peak=float(np.nanmedian(wz[s][:, 3])))
                                    for s in wz.files if s != "names"})
    dn = {s: den_summary(dens, "s", s)[0] for s in rep["slug"].unique()}
    den_tot = {s: sum(v.values()) for s, v in dn.items()}
    # bursts: events with many candidate channels
    bs = {nb: burst_events(rep, nb) for nb in (3, 5, 10, 20)}
    NB_ISO = NB_BURST
    brows = []
    for slug in sorted(rep["slug"].unique()):
        m = (rep["slug"] == slug).to_numpy()
        chs_all = rep.loc[m].drop_duplicates("ch")
        for nb in bs:
            mb = bs[nb][m]
            nbc = rep.loc[m].loc[mb].drop_duplicates("ch")
            win = rep.loc[m]["inwin"].to_numpy()
            brows.append(dict(runtype=slug, nb=nb, n_cand=len(chs_all), n_cand_in_bursts=len(nbc),
                              frac=len(nbc) / max(len(chs_all), 1),
                              n_win_sec=int(win.sum()), n_win_sec_in_bursts=int((win & mb).sum()),
                              n_burst_events=int(nbc.drop_duplicates(["run", "split", "event"]).shape[0])))
    pd.DataFrame(brows).to_csv(out / "burst_events.tsv", sep="\t", index=False, float_format="%.4g")
    bb = rep[bs[NB_ISO]]
    bev = []
    for (run, split, ev), g in bb.groupby(["run", "split", "event"]):
        g1 = g.drop_duplicates("ch")
        bev.append(dict(run=run, split=split, event=ev, n_cand_channels=len(g1),
                        n_win_sec=int(g["inwin"].sum()), dt_median=float(g["dt"].median()),
                        E_ratio_median=float(g["E_ratio"].median()),
                        time_dom_median=float(g["time_dom"].median()),
                        chi2dof_median=float((g["chi2_full"] / g["dof"]).median()),
                        E_dom_median=float(g["E_dom"].median())))
    pd.DataFrame(bev).to_csv(out / "burst_events_list.tsv", sep="\t", index=False, float_format="%.4g")
    summ["burst_events_nb10"] = bev
    summ["bursts"] = brows
    # --------------------------------------------------- real-data rules
    specs = rule_specs()
    rrows = []
    cols = ["dt", "E_ratio", "dchi2", "chi2_full", "dof", "fwhm", "tfall"]
    for burst_nb in (0, 5, NB_ISO):
        for slug in sorted(rep["slug"].unique()) + ["ALL"]:
            m = np.ones(len(rep), bool) if slug == "ALL" else (rep["slug"] == slug).to_numpy()
            if burst_nb:
                m = m & ~bs[burst_nb]
            sub = rep[m]
            d = {c: sub[c].to_numpy() for c in cols}
            win = sub["inwin"].to_numpy()
            dtv = sub["dt"].to_numpy()
            chs = sub["ch"].to_numpy()
            den = (sum(den_tot.values()) if slug == "ALL" else den_tot[slug])
            n_win0 = int(win.sum())
            fe0 = flat_expect(dtv)
            n_out0 = int((~win).sum())
            n0, r0 = _ch_rate(np.ones(len(sub), bool), chs, den)
            for name, spec in specs:
                keep = rule_keep(d, spec)
                nk, rk = _ch_rate(keep, chs, den)
                nw = int((keep & win).sum())
                rrows.append(dict(
                    rule=name, burst_veto_nb=burst_nb, runtype=slug, den=den,
                    n_win=n_win0, n_out=n_out0, flat_win=fe0,
                    removed_win=1 - nw / max(n_win0, 1),
                    removed_out=1 - int((keep & ~win).sum()) / max(n_out0, 1),
                    excess_surviving=(nw - fe0) / max(n_win0 - fe0, 1e-9),
                    n_cand_after=nk, rate_after_pct=rk, rate_before_pct=r0,
                    n_win_after=nw))
    rr = pd.DataFrame(rrows)
    rr.to_csv(out / "rules_real.tsv", sep="\t", index=False, float_format="%.4g")
    # --------------------------------------------------- pseudo-data rules
    sp = pd.read_csv(out / "synth_pulses.tsv.gz", sep="\t")
    sp = sp[sp["npk"] > 0].copy()
    fw = np.zeros(len(sp))
    tf = np.zeros(len(sp))
    for slug in sp["slug"].unique():
        m = (sp["slug"] == slug).to_numpy()
        w = wz[slug]
        fw[m] = w[sp.loc[m, "mod"].astype(int).to_numpy(), 2]
        tf[m] = w[sp.loc[m, "mod"].astype(int).to_numpy(), 1]
    sp["fwhm"], sp["tfall"] = fw, tf
    sp["reported"] = ((sp["dom"] == 1) | ((sp["E_ratio"] > R_MODEL) & (sp["edge"] == 0))).astype(int)
    sp = sp[sp["reported"] == 1].copy()
    sp["evid"] = (sp["slug"].astype("category").cat.codes.astype(np.int64) * 10**9
                  + (sp["dt_cell"] * 10).astype(np.int64) * 10**6
                  + (sp["ratio_cell"] * 1000).astype(np.int64) * 10**3 * 0 + sp["ev"].astype(np.int64))
    sp["evid"] = (sp["slug"].astype(str) + "|" + sp["dt_cell"].astype(str) + "|"
                  + sp["ratio_cell"].astype(str) + "|" + sp["ev"].astype(str))
    sp["evid"] = sp["evid"].astype("category").cat.codes
    nev = sp["evid"].max() + 1
    srows = []
    d_all = {c: sp[c].to_numpy() for c in cols}
    isdom = sp["dom"].to_numpy() == 1
    evid = sp["evid"].to_numpy()
    ni, nbse = sp["near_inj"].to_numpy(), sp["near_base"].to_numpy()
    cell = sp.groupby("evid")[["slug", "dt_cell", "ratio_cell"]].first()
    cell_idx = cell.index.to_numpy()
    for name, spec in specs:
        keep = isdom | rule_keep(d_all, spec)
        na = np.bincount(evid, weights=keep & (ni == 1), minlength=nev)
        nb_ = np.bincount(evid, weights=keep & (nbse == 1), minlength=nev)
        nu = np.bincount(evid, weights=keep & ((ni == 1) | (nbse == 1)), minlength=nev)
        match = (na >= 1) & (nb_ >= 1) & (nu >= 2)
        fake = np.bincount(evid, weights=keep & (ni == 0) & (nbse == 0), minlength=nev)
        nsec = np.bincount(evid, weights=keep & ~isdom, minlength=nev)
        t = pd.DataFrame(dict(slug=cell["slug"].to_numpy(), dt_cell=cell["dt_cell"].to_numpy(),
                              ratio_cell=cell["ratio_cell"].to_numpy(),
                              match=match[cell_idx].astype(float), fake=fake[cell_idx], nsec=nsec[cell_idx]))
        for grp, sel in (("ratio>=0.2&dt>=24", (t["ratio_cell"] >= 0.2) & (t["dt_cell"] >= 24)),
                         ("ratio>=0.3&dt>=24", (t["ratio_cell"] >= 0.3) & (t["dt_cell"] >= 24)),
                         ("ratio=0.2&dt24-32", (t["ratio_cell"] == 0.2) & (t["dt_cell"].isin([24, 32]))),
                         ("ratio=0.2&dt>=48", (t["ratio_cell"] == 0.2) & (t["dt_cell"] >= 48)),
                         ("ratio>=0.3&dt24-32", (t["ratio_cell"] >= 0.3) & (t["dt_cell"].isin([24, 32]))),
                         ("ratio=0.1&dt>=24", (t["ratio_cell"] == 0.1) & (t["dt_cell"] >= 24)),
                         ("dt<24", t["dt_cell"] < 24), ("all", t["dt_cell"] > 0)):
            for sl in ("ALL", *sorted(t["slug"].unique())):
                ss = sel if sl == "ALL" else sel & (t["slug"] == sl)
                if ss.sum() == 0:
                    continue
                srows.append(dict(rule=name, group=grp, slug=sl, n_events=int(ss.sum()),
                                  eff_match=float(t.loc[ss, "match"].mean()),
                                  fake_rate=float(t.loc[ss, "fake"].mean()),
                                  n_sec_per_event=float(t.loc[ss, "nsec"].mean())))
    sr = pd.DataFrame(srows)
    sr.to_csv(out / "rules_synth.tsv", sep="\t", index=False, float_format="%.5g")
    # --------------------------------------------------- joint summary
    base = sr[(sr["rule"] == "none")].set_index(["group", "slug"])
    rows = []
    for name, spec in specs:
        r1 = rr[(rr["rule"] == name) & (rr["burst_veto_nb"] == 0)].set_index("runtype")
        r2 = rr[(rr["rule"] == name) & (rr["burst_veto_nb"] == NB_ISO)].set_index("runtype")
        s1 = sr[sr["rule"] == name].set_index(["group", "slug"])
        row = dict(rule=name)
        for sl in ("0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "ALL"):
            row[f"removed_win_{sl}"] = r1.loc[sl, "removed_win"]
            row[f"removed_out_{sl}"] = r1.loc[sl, "removed_out"]
            row[f"rate_after_{sl}"] = r1.loc[sl, "rate_after_pct"]
        row["rate_before_ALL"] = r1.loc["ALL", "rate_before_pct"]
        row["removed_win_isolated_ALL"] = r2.loc["ALL", "removed_win"]
        row["removed_out_isolated_ALL"] = r2.loc["ALL", "removed_out"]
        row["rate_after_isolated_ALL"] = r2.loc["ALL", "rate_after_pct"]
        row["rate_before_isolated_ALL"] = r2.loc["ALL", "rate_before_pct"]
        for grp in ("ratio>=0.2&dt>=24", "ratio>=0.3&dt>=24", "ratio=0.2&dt24-32", "ratio=0.2&dt>=48",
                    "ratio>=0.3&dt24-32"):
            row[f"d_eff[{grp}]"] = s1.loc[(grp, "ALL"), "eff_match"] - base.loc[(grp, "ALL"), "eff_match"]
        row["eff_base[ratio>=0.2&dt>=24]"] = base.loc[("ratio>=0.2&dt>=24", "ALL"), "eff_match"]
        row["d_fake[all]"] = s1.loc[("all", "ALL"), "fake_rate"] - base.loc[("all", "ALL"), "fake_rate"]
        row["fake_base[all]"] = base.loc[("all", "ALL"), "fake_rate"]
        rows.append(row)
    js = pd.DataFrame(rows)
    js.to_csv(out / "rules_summary.tsv", sep="\t", index=False, float_format="%.4g")
    summ["rules_summary_file"] = "rules_summary.tsv"
    # figure: removal of the window excess vs eff loss
    fig, axs = plt.subplots(1, 2, figsize=(13, 5))
    for ax, grp in zip(axs, ("ratio>=0.2&dt>=24", "ratio>=0.3&dt>=24")):
        x = -100 * js[f"d_eff[{grp}]"]
        y = 100 * js["removed_win_ALL"]
        ax.scatter(x, y, s=12, c="0.6")
        for kind, c in (("dchi2", "C1"), ("dt>=", "C0")):
            mm = js["rule"].str.startswith(kind) if kind == "dchi2" else js["rule"].str.startswith("dt>=")
            ax.scatter(x[mm], y[mm], s=18, c=c, label=kind)
        for i in js.index:
            if js.loc[i, "rule"] in ("dt>=40|r>=0.19", "dt>=40|r>=0.25", "dt>=40|r>=0.29", "dt>=1.0FWHM|r>=0.29",
                                     "dt>=32|r>=0.29", "dchi2>=1000", "dt>=40|r>=0.49",
                                     "dt>=24|r>=0.29"):
                ax.annotate(js.loc[i, "rule"], (x[i], y[i]), fontsize=6)
        ax.set_xlabel(f"eff_match loss [percentage points], cells {grp}")
        ax.set_ylabel("window (+12..40 ns) secondaries removed [%]  (e-p splits)")
        ax.legend(fontsize=7)
        ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(out / "rules_tradeoff.png", dpi=110)
    plt.close(fig)
    update_summary(out, rules=summ)
    print(f"[rules] wrote rules_real.tsv, rules_synth.tsv, rules_summary.tsv ({len(specs)} rules)")


# ------------------------------------------------------- neighbours stage ---
NB_SPLITS = [(24590, 10), (25320, 50)]


def nb_task(job):
    """Re-read the EVIO events of the isolated window candidates and test whether the
    neighbouring modules of the same event show a pulse at the same time as the
    candidate's secondary (a second shower would; a channel-level feature would not)."""
    run, split, fpath, targets = job["run"], job["split"], job["path"], job["targets"]
    calib = _HC_CALIB_MAP.get(run)
    calib_path = str(_DATABASE_DIR / "calibration" / calib) if calib else ""
    p = C.setup_pipeline(evio_path=fpath, run_num=run, hc_calib_file=calib_path)
    mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    wcfg = p.wave_ana.cfg
    names, pwo4 = [m.name for m in mods], [bool(m.is_pwo4()) for m in mods]
    slug = run_slug(run)
    bank, _ = mf_pulse.build_bank(
        names, pwo4, SPL_PAT.format(runtype=slug),
        dict(smooth_order=int(wcfg.smooth_order), int_tail_ratio=float(wcfg.int_tail_ratio),
             tail_break_n=int(wcfg.tail_break_n)), column="sampled")
    rec = mf_pulse.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4, r_model=R_MODEL,
                                   ped_nsamples=int(wcfg.ped_nsamples))
    xy = np.array([[m.x, m.y, m.size_x] for m in mods])
    sa = SimpleNamespace(dominant_min_mev=THR_MEV, dominant_frac=0.5,
                         secondary_frac=R_MODEL, include_clipped=False)
    sc = ppe.Scanner(sa, p, bank, rec, 1.0, float(p.geo.Ebeam))
    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    ch.open_auto(str(fpath))
    tset = {t[0] for t in targets}
    evrows = {}
    seen = 0
    while ch.read() == dec.Status.success and seen < len(tset):
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            fe = d["event"]
            if int(fe.info.trigger_bits) != C.PHYSICS_TRIGGER_BITS:
                continue
            if int(fe.info.event_number) not in tset:
                continue
            sc._reset()
            sc.collect(split, fe)
            if not sc.rows:
                continue
            Y = np.stack(sc.rows)
            mid = np.array([m[3] for m in sc.meta])
            r = rec.run(Y, mid, sc.overflow)
            pul = {}
            for j, (c_, t0, am) in enumerate(zip(r.p_ch, r.p_t0, r.p_amp)):
                pul.setdefault(int(mid[c_]), []).append((float(t0), float(am)))
            evrows[int(fe.info.event_number)] = pul
            seen += 1
    offs = [-50, -40, -30, 30, 40, 50, 60, 70]
    out = []
    for ev, mi, t_dom, t_sec in targets:
        pul = evrows.get(ev)
        if not pul:
            continue
        near = [m for m in pul if m != mi and math.hypot(*(xy[m, :2] - xy[mi, :2])) < 1.6 * xy[mi, 2]]
        # neighbours with an in-time pulse (within 12 ns of the candidate's dominant)
        for m in near:
            dom = [(t, a_) for t, a_ in pul[m] if abs(t - t_dom) < 12]
            if not dom:
                continue
            a_nb = max(a_ for t, a_ in dom)
            rowx = dict(run=run, split=split, event=ev, mod=mi, nbr=m, A_nbr=a_nb)
            for off in [0] + offs:
                ts = t_sec + off
                rowx[f"sec_{off}"] = int(any(abs(t - ts) < 5.0 and abs(t - t_dom) >= 10
                                             for t, a_ in pul[m]))
            out.append(rowx)
    return out


def cmd_neighbours(a):
    out = Path(a.out)
    dfs, dens, metas, names = load_scan(out)
    df = dfs["s"]
    r = df[(df["reported"] == 1) & (df["edge"] == 0)]
    jobs = []
    for run, split in NB_SPLITS:
        s = r[(r["run"] == run) & (r["split"] == split) & r["inwin"]]
        nwin = s.groupby("event")["ch"].nunique()
        iso = s[s["event"].isin(nwin[nwin <= 2].index)].sort_values("A_dom", ascending=False)
        iso = iso.drop_duplicates("ch").sample(frac=1.0, random_state=1).head(a.nb_events)
        tg = [(int(x.event), int(x.mod), float(x.time_dom), float(x.time_dom + x.dt)) for x in iso.itertuples()]
        # t_dom and t_sec in the mf p_t0 (onset) convention: shift by the module's t50
        fp = dict(ppe.local_splits(a.evio_dir, run)).get(split)
        jobs.append(dict(run=run, split=split, path=str(fp), targets=[]))
        jobs[-1]["_tg"] = tg
    # onset-convention conversion needs t50: onset = time - t50 for both pulses, so the
    # difference is unchanged; dominant time is only used for the +-12 ns in-time test
    # (t50 ~ 8 ns shifts both identically, the neighbours' own t0 are onsets) -> use onsets
    for j in jobs:
        z = np.load(SPL_PAT.format(runtype=run_slug(j["run"])), allow_pickle=False)
        t50 = {str(n): float(v) for n, v in zip(z["names"], z["t50"])}
        nm = np.array([str(n) for n in names])
        j["targets"] = [(ev, mi, td - t50.get(nm[mi], 0.0), ts - t50.get(nm[mi], 0.0))
                        for ev, mi, td, ts in j.pop("_tg")]
    res = []
    with cf.ProcessPoolExecutor(max_workers=min(a.jobs, len(jobs))) as ex:
        for rr in ex.map(nb_task, jobs):
            res += rr
    t = pd.DataFrame(res)
    t.to_csv(out / "neighbours.tsv", sep="\t", index=False)
    summ = {}
    edges = (0, 30, 100, 300, 1e9)
    for run in sorted(t["run"].unique()):
        s = t[t["run"] == run]
        rec = dict(n_candidates=int(s.drop_duplicates(["event", "mod"]).shape[0]),
                   n_neighbours=len(s), by_A_nbr={})
        for lo, hi in zip(edges[:-1], edges[1:]):
            q = s[(s["A_nbr"] >= lo) & (s["A_nbr"] < hi)]
            if len(q) == 0:
                continue
            ctrl = [q[f"sec_{o}"].mean() for o in (-50, -40, -30, 30, 40, 50, 60, 70)]
            rec["by_A_nbr"][f"{lo:g}-{hi:g}"] = dict(
                n=len(q), frac_at_candidate_dt=float(q["sec_0"].mean()),
                frac_control=float(np.mean(ctrl)))
        summ[int(run)] = rec
        print(f"[neighbours] run {run}: {rec}")
    update_summary(out, neighbours=summ)


def analyze_window_shape(out, dfs, dens, slugs, summ, names):
    """Fixed-ratio / fixed-delay test of the window secondaries of the isolated (non-burst)
    events, and the per-module spread of the window rate."""
    df = dfs["s"]
    rep = df[(df["reported"] == 1) & (df["edge"] == 0) & df["inwin"] & ~df["burst"]].copy()
    rep["A_sec"] = rep["A_ratio"] * rep["A_dom"]
    rows = []
    for slug in slugs + ["ALL"]:
        sl = rep if slug == "ALL" else rep[rep["slug"] == slug]
        for ab in range(-1, 5):
            t = sl if ab < 0 else sl[sl["ab"] == ab]
            if len(t) < 10:
                continue
            rows.append(dict(runtype=slug, A_bin="all" if ab < 0 else AMP_LABELS[ab], n=len(t),
                             E_ratio_q16_50_84=_q(t["E_ratio"]), A_ratio_q16_50_84=_q(t["A_ratio"]),
                             dt_q16_50_84=_q(t["dt"]), A_sec_q16_50_84=_q(t["A_sec"]),
                             A_dom_q16_50_84=_q(t["A_dom"]),
                             chi2dof_q16_50_84=_q(t["chi2_full"] / t["dof"])))
    pd.DataFrame(rows).to_csv(out / "window_secondaries_isolated.tsv", sep="\t", index=False,
                              float_format="%.4g")
    summ["window_secondaries_isolated"] = rows
    # evidence figure: fixed ratio / fixed delay, burst events
    allrep = df[(df["reported"] == 1) & (df["edge"] == 0)]
    fig, axs = plt.subplots(1, 4, figsize=(19, 4.3))
    for slug, c in zip(slugs, ("C0", "C1", "C2")):
        t = allrep[(allrep["slug"] == slug) & ~allrep["burst"]]
        h, e = np.histogram(t["dt"], np.arange(-100, 160, 2))
        axs[0].stairs(h, e, color=c, label=RUNTYPE_LABEL.get(slug, slug))
    axs[0].axvspan(*WIN, color="0.92", zorder=0)
    axs[0].set_yscale("log")
    axs[0].set_xlabel("dt (ns)")
    axs[0].set_ylabel("reported secondaries / 2 ns (burst events removed)")
    axs[0].legend(fontsize=7)
    axs[1].hist2d(np.log10(rep["A_dom"]), rep["E_ratio"].clip(0, 1), bins=[40, 40],
                  range=[[1.4, 3.6], [0, 1]], cmap="viridis", cmin=1)
    axs[1].set_xlabel("log10 A_dom (ADC)")
    axs[1].set_ylabel("E_sec / E_dom")
    axs[1].set_title("window secondaries, isolated events")
    axs[2].loglog(rep["A_dom"], rep["A_sec"], ",", alpha=0.4)
    xx = np.logspace(1.7, 3.4, 10)
    axs[2].loglog(xx, 0.18 * xx, "r-", lw=0.8, label="A_sec = 0.18 A_dom")
    axs[2].set_xlabel("A_dom (ADC)")
    axs[2].set_ylabel("A_sec (ADC)")
    axs[2].legend(fontsize=7)
    for slug, c in zip(slugs, ("C0", "C1", "C2")):
        t = allrep[allrep["slug"] == slug].drop_duplicates("ch")
        n = t.groupby(["split", "event"]).size()
        h, e = np.histogram(n, np.logspace(0, 3.3, 24))
        axs[3].stairs(h, e, color=c, label=RUNTYPE_LABEL.get(slug, slug))
    axs[3].set_xscale("log")
    axs[3].set_yscale("log")
    axs[3].set_xlabel("candidate channels per event")
    axs[3].set_ylabel("events")
    axs[3].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(out / "window_secondaries_evidence.png", dpi=105)
    plt.close(fig)
    # per-module spread of the window-candidate rate (isolated events)
    nm = [str(n) for n in names]
    mrows = {}
    for slug in slugs:
        dd = [d for d in dens if d["slug"] == slug]
        denm = sum(d["den"]["s"] for d in dd).sum(axis=1)
        t = rep[rep["slug"] == slug].drop_duplicates("ch")
        cnt = np.bincount(t["mod"].to_numpy(), minlength=len(nm))
        ok = denm >= 200
        rate = 100.0 * cnt[ok] / denm[ok]
        top = np.argsort(-np.where(ok, cnt, -1))[:10]
        mrows[slug] = dict(n_modules=int(ok.sum()), rate_pct_q10_50_90_99=_q(rate, (10, 50, 90, 99)),
                           n_modules_rate_gt2pct=int((rate > 2).sum()),
                           top10_share=float(cnt[top].sum() / max(cnt.sum(), 1)),
                           top5=[(nm[i], int(cnt[i]), int(denm[i])) for i in top[:5]])
    summ["window_rate_per_module_isolated"] = mrows
    # fraction of channels per time-since-trigger: dominant onset of the isolated window sample
    summ["window_time_dom_q16_50_84_isolated"] = _q(rep["time_dom"])
    summ["nonwindow_time_dom_q16_50_84"] = _q(df[(df["reported"] == 1) & (df["edge"] == 0)
                                                 & ~df["inwin"] & ~df["burst"]]["time_dom"])


def analyze_crossbank(out, dfs, summ):
    """The same waveforms under both banks: where does the other bank put the secondary?"""
    res = {}
    rows = []
    for (b1, b2) in (("s", "h"), ("h", "s")):
        d1, d2 = dfs[b1], dfs[b2]
        r1 = d1[(d1["reported"] == 1) & (d1["edge"] == 0) & d1["inwin"] & ~d1["burst"]]
        c1 = r1.drop_duplicates("ch")
        r2 = d2[(d2["reported"] == 1) & (d2["edge"] == 0)]
        first2 = r2.sort_values("E_ratio", ascending=False).drop_duplicates("ch").set_index("ch")
        for slug in sorted(c1["slug"].unique()):
            cs = c1[c1["slug"] == slug]
            have = cs["ch"].isin(first2.index)
            other = first2.loc[cs.loc[have, "ch"]]
            res[f"{BANK_NAME[b1]}->{BANK_NAME[b2]}|{slug}"] = dict(
                n_window_channels=len(cs), frac_other_bank_has_reported=float(have.mean()),
                other_dt_q16_50_84=_q(other["dt"]), other_Eratio_q16_50_84=_q(other["E_ratio"]),
                other_dt_in_own_window=float(((other["dt"] >= WIN[0]) & (other["dt"] < WIN[1])).mean())
                if len(other) else float("nan"))
    summ["cross_bank_isolated_window_channels"] = res


# ------------------------------------------------------------------- main ---
def build_parser():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("cmd", choices=("scan", "analyze", "synth", "rules", "neighbours", "all"))
    ap.add_argument("--out", default=str(OUT_DEFAULT))
    ap.add_argument("--splits", nargs="+", default=SPLITS_DEFAULT,
                    help="run:split pairs for the scan [25320:50 25320:100 24503:10 "
                         "24503:20 24590:10 24590:20]")
    ap.add_argument("--jobs", type=int, default=6, help="processes (max 6 here)")
    ap.add_argument("--evio-dir", default=_DEFAULT_EVIO_DIR)
    ap.add_argument("--max-events", type=int, default=0)
    ap.add_argument("--force", action="store_true", help="redo finished scan splits")
    ap.add_argument("--nb-events", type=int, default=400,
                    help="neighbours: isolated window-candidate events per split [400]")
    return ap


def main(argv=None):
    a = build_parser().parse_args(argv)
    if a.cmd in ("scan", "all"):
        cmd_scan(a)
    if a.cmd in ("synth", "all"):
        cmd_synth(a)
    if a.cmd in ("analyze", "all"):
        cmd_analyze(a)
    if a.cmd in ("neighbours", "all"):
        cmd_neighbours(a)
    if a.cmd in ("rules", "all"):
        cmd_rules(a)


if __name__ == "__main__":
    main()
