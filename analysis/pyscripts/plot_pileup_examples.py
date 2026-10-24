#!/usr/bin/env python3
"""plot_pileup_examples.py -- single-module pile-up rates, heat maps and examples (matched filter).

Two stages.  Scan: one task per (run, split) walks the EVIO split, reconstructs
every PbWO4 waveform with the matched filter (mf_pulse.MFReconstructor, imported
unchanged) and keeps the waveforms with a clear secondary pulse:

  * >= 2 pulses accepted by the mf (5 sigma in MF units and >= 5 ADC, MIN_SEP
    12 ns; the secondary acceptance r_model is set to --secondary-frac, so a
    pulse inside a parent's kernel span needs A >= 0.10 A_parent),
  * dominant pulse (highest E) with E_dom > --dominant-min-mev MeV (when > 0)
    or else E_dom > --dominant-frac * E_beam,
  * at least one other pulse with E / E_dom > --secondary-frac.  Only those are
    reported as secondaries (more than one is allowed); weaker accepted pulses
    stay in the fit and are drawn in grey.

The denominator is the channels with such a dominant pulse; the numerator the
candidates among them, both counted per module.  E = A * I_int(phase) * kappa *
cal_factor (mf_pulse.pulse_energy_units + Module.energize; no LMS correction;
the sampled bank adds its amplitude factor f(A) inside pulse_energy_units).
Clipped waveforms (any sample at the FADC overflow) are skipped by default: their
mf amplitude is an extrapolation over the masked samples (--include-clipped keeps
the ones the mf does not hand to the WaveAnalyzer fallback).  Pulses with an mf
onset within 2 samples of the last kernel position (sample ns - L) are truncated
and mis-fitted: such a dominant drops the channel, such a secondary never counts.
A channel with any pulse onset before half the pedestal window (ped_nsamples/2
samples) is dropped: the pulse biases the pedestal, sigma and all amplitudes.
Per secondary, dchi2 = chi2(fit without it, other amplitudes re-solved at fixed
onsets, same samples) - chi2(full fit).

Pile-up fix (both off by default = the behaviour of the 20261023 scan; the
production rerun uses --burst-veto 10 --window-rule 40:0.16):
  --burst-veto N   an event with >= N candidate channels (channels with a
                   secondary under the rule above, in that event) is removed from
                   numerator and denominator for all its channels (burst events:
                   late dominant times, chi2/dof in the thousands);
  --window-rule DT_MAX:RATIO  a secondary with 0 <= dt < DT_MAX ns and E/E_dom <
                   RATIO is not a pile-up secondary (isolated second bumps at
                   dt ~ 23 ns, ratio ~ 0.18); the channel still qualifies through
                   another secondary.
The raw counts stay in the outputs (rate_raw); the fixed ones are (veto + rule),
and both are also given without the --hot-channels.  The candidate TSV lists all
raw candidates with a trailing column removed_by (burst / window / empty).

Report: per run rate_m = num_m / den_m per module with binomial errors, a heat
map (rate in %, modules with den_m < 20 grey; denominator panel; fixed counts
when the scan has them, hot channels outlined) and the first --n-examples KEPT
candidates of the run's first split drawn with the selected bank.

Outputs in --out:
  scan/counts_<run>_<split>.npz        per-module den_m, num_m, counters, geometry
  scan/cands_<run>_<split>.tsv         every candidate of the split
  scan/examples_<run>_<split>.npz      first-split run only: first 20 candidates
                                       (raw waveform, ped, sigma, fit components)
  pileup_rate_summary.tsv, pileup_rate_module_<run>.tsv, pileup_rate_<run>.txt,
  pileup_hot_channels.tsv
  pileup_heatmap_<run>.png, pileup_heatmaps_all_runs.pdf
  pileup_examples_all_runs.pdf         per run a cover page + the examples

  MPLBACKEND=Agg python3 plot_pileup_examples.py --runs 25320 --splits 50 --out DIR
  python3 plot_pileup_examples.py --jobs 12 --mf-templates \
      "output/templates_spline/{runtype}/templates_{runtype}.npz" \
      --mf-scale 1.0 --dominant-min-mev 50 --out output/pileup_20261024
  python3 plot_pileup_examples.py --report-only --out output/pileup_20261024
"""

from __future__ import annotations

import argparse
import concurrent.futures as cf
import json
import math
import os
import re
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                       # noqa: E402
from matplotlib.collections import PatchCollection    # noqa: E402
from matplotlib.colors import LogNorm                 # noqa: E402
from matplotlib.patches import Rectangle              # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402

import _common as C                                   # noqa: E402
from prad2py import dec                               # noqa: E402
import mf_pulse                                       # noqa: E402
from benchmark_hycal_timing_scan import (             # noqa: E402
    _DATABASE_DIR, _DEFAULT_EVIO_DIR, _EXCLUDED_RUNS, _HC_CALIB_MAP)
from benchmark_hycal_timing_test_20261021 import X17_RUNS  # noqa: E402

HERE = Path(__file__).resolve().parent
BLOCK = 256           # events per batched mf pass
ALL_RUNS = [24386, 24503, 24512, 24590, 24603, 25228, 25257, 25260, 25299,
            25308, 25320, 25800, 26138, 26230]
N_RAW = 20            # candidates of a run's first split stored with raw waveforms
MIN_DEN = 20          # heat map: modules with fewer denominator channels are grey
HOT_DEFAULT = "W335,W1,W1019,W1020,W934"
RAW_CAP_BURST = 3     # examples of a burst-vetoed scan: raw waveforms for 3 x N_RAW candidates


def parse_window_rule(txt):
    """'DT_MAX:RATIO' -> (dt_max_ns, ratio) or None."""
    if not txt:
        return None
    try:
        d, r = txt.split(":")
        return float(d), float(r)
    except ValueError:
        sys.exit(f"[ERROR] --window-rule expects DT_MAX:RATIO (e.g. 40:0.17), got {txt!r}")


# ------------------------------------------------------------------ helpers --
def local_splits(evio_dir, run):
    pat = re.compile(rf"^prad_0*{run}\.evio\.(\d+)$", re.IGNORECASE)
    out = []
    for f in sorted(Path(evio_dir).iterdir()):
        m = pat.match(f.name)
        if m:
            out.append((int(m.group(1)), f))
    return sorted(out)


def fit_terms(y, w, sig, bank, m, t0s, amps):
    """Model samples per pulse, fit mask, chi2 of the full fit, dchi2 per pulse.

    Same definitions as MFReconstructor (union of the pulses' LW windows, clipped
    samples masked, var = sigma^2 + (MODEL_FLOOR * model)^2).  dchi2[k] is the
    chi2 increase when pulse k is removed and the others are re-solved (linear
    LS at fixed onsets, non-positive amplitudes dropped and re-solved) on the
    same samples.
    """
    ns = y.size
    i = np.arange(ns)
    S = np.array([bank.shape(m, mf_pulse.DT * i - t0, dtype=np.float64)
                  for t0 in t0s])
    wm = np.zeros(ns)
    for t0 in t0s:
        w0 = int(math.ceil(t0 / mf_pulse.DT)) - mf_pulse.LEAD
        wm[max(w0, 0):max(w0 + mf_pulse.LW, 0)] = 1.0
    wm *= w

    def chi2(Ss, a):
        mdl = a @ Ss if len(a) else np.zeros(ns)
        var = sig ** 2 + (mf_pulse.MODEL_FLOOR * mdl) ** 2
        return float((wm * (y - mdl) ** 2 / var).sum())

    full = chi2(S, np.asarray(amps))
    dof = max(wm.sum() - 2.0 * len(t0s), 1.0)
    dchi2 = []
    for k in range(len(t0s)):
        keep = [j for j in range(len(t0s)) if j != k]
        while True:
            Sk = S[keep]
            G = (Sk * wm) @ Sk.T
            b = (Sk * wm) @ y
            a = np.linalg.solve(G + 1e-9 * np.eye(len(keep)) * max(G.max(), 1e-12), b)
            if len(keep) <= 1 or (a > 0).all():
                break
            del keep[int(np.argmin(a))]
        dchi2.append(chi2(Sk, np.maximum(a, 0.0)) - full)
    return S, wm, full, dof, dchi2


# --------------------------------------------------------------------- scan --
class Scanner:
    def __init__(self, args, p, bank, rec, kappa, ebeam, keep_raw=False):
        self.a, self.p, self.bank, self.rec = args, p, bank, rec
        self.kappa, self.ebeam, self.keep_raw = kappa, ebeam, keep_raw
        self.mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
        self.thr = (args.dominant_min_mev if args.dominant_min_mev > 0
                    else args.dominant_frac * ebeam)      # dominant E threshold, MeV
        self.den_m = np.zeros(len(self.mods), dtype=np.int64)
        self.num_m = np.zeros(len(self.mods), dtype=np.int64)
        self.pwo4 = np.array([bool(m.is_pwo4()) for m in self.mods])
        self.cal = np.array([m.energize(1.0) for m in self.mods])
        self.overflow = int(p.wave_ana.cfg.overflow)
        self.ns = None
        self.cands = []          # all candidates (dicts)
        self.n_den = 0           # channels with E_dom above the threshold
        self.n_fb = 0            # clipped channels handed to the WaveAnalyzer fallback
        self.n_clip = 0          # clipped channels skipped (no --include-clipped)
        self.n_edge = 0          # E_dom above threshold but dominant onset at the window end
        self.n_early = 0         # E_dom above threshold but a pulse inside the pedestal window
        self.n_events = 0
        self.burst_n = int(args.burst_veto)
        self.win = parse_window_rule(args.window_rule)
        self.raw_cap = N_RAW * (RAW_CAP_BURST if self.burst_n > 0 else 1)
        self.n_raw = 0
        self.ev_rec = []         # per denominator channel: (event, module, raw cand, fixed cand)
        self.den_fix_m = np.zeros(len(self.mods), dtype=np.int64)
        self.num_fix_m = np.zeros(len(self.mods), dtype=np.int64)
        self.n_burst_events = self.n_burst_den = self.n_burst_cand = 0
        self.n_win_removed = self.n_win_removed_burst = 0
        self._reset()

    def _reset(self):
        self.rows, self.meta = [], []    # samples, (split, ev, trig, mod index)

    def collect(self, split, fadc_evt):
        p = self.p
        ev = int(fadc_evt.info.event_number)
        trig = int(fadc_evt.info.trigger_bits)
        for ri in range(fadc_evt.nrocs):
            roc = fadc_evt.roc(ri)
            if not roc.present:
                continue
            crate = p.crate_map.get(roc.tag)
            if crate is None:
                continue
            for s in roc.present_slots():
                slot = roc.slot(s)
                for c in slot.present_channels():
                    mod = p.hycal.module_by_daq(crate, s, c)
                    if mod is None or not mod.is_hycal() or not self.pwo4[mod.index]:
                        continue
                    cd = slot.channel(c)
                    if cd.nsamples <= 0:
                        continue
                    if self.ns is None:
                        self.ns = cd.nsamples
                    if cd.nsamples != self.ns:
                        continue
                    self.rows.append(np.asarray(cd.samples, dtype=np.uint16))
                    self.meta.append((split, ev, trig, mod.index))

    def flush(self):
        if not self.rows:
            return
        a, bank = self.a, self.bank
        Y = np.stack(self.rows)
        mid = np.array([m[3] for m in self.meta], dtype=np.int64)
        n = Y.shape[0]
        res = self.rec.run(Y, mid, self.overflow)
        self.n_fb += int(res.wa_fb.sum())
        skip = res.wa_fb if a.include_clipped else res.clipped
        if not a.include_clipped:
            self.n_clip += int(res.clipped.sum())
        pc = res.p_ch
        if pc.size:
            Eu = mf_pulse.pulse_energy_units(bank, mid[pc], res.p_t0, res.p_amp,
                                             self.kappa)
            E = Eu * self.cal[mid[pc]]
            starts = np.searchsorted(pc, np.arange(n + 1))
            emax = np.zeros(n)
            np.maximum.at(emax, pc, E)
            high = np.nonzero((emax > self.thr) & ~skip)[0]
            # the mf cannot place an onset past sample ns - L: pulses at that
            # edge are truncated and mis-fitted, so they never count
            t_edge = mf_pulse.DT * (Y.shape[1] - self.rec.L - 2)
            t_early = mf_pulse.DT * self.rec.ped_n / 2.0
            for r in high:
                lo, hi = starts[r], starts[r + 1]
                Er = E[lo:hi]
                kd = int(np.argmax(Er))
                if res.p_t0[lo:hi].min() < t_early:
                    self.n_early += 1
                    continue
                if res.p_t0[lo + kd] >= t_edge:
                    self.n_edge += 1
                    continue
                self.n_den += 1
                self.den_m[self.meta[r][3]] += 1
                sec = [k for k in range(hi - lo)
                       if k != kd and Er[k] / Er[kd] > a.secondary_frac
                       and res.p_t0[lo + k] < t_edge]
                sk = sec
                if self.win is not None:
                    dmax, rmax = self.win
                    sk = [k for k in sec
                          if not (0.0 <= res.p_time[lo + k] - res.p_time[lo + kd] < dmax
                                  and Er[k] / Er[kd] < rmax)]
                self.ev_rec.append((self.meta[r][1], self.meta[r][3],
                                    int(bool(sec)), int(bool(sk))))
                if sec:
                    self._record(r, Y[r], res, lo, hi, Er, kd, sec, sk, t_edge)
        self._reset()

    def finalize(self):
        """Burst veto and fixed counts (needs all channels of each event)."""
        nm = len(self.mods)
        n = len(self.ev_rec)
        arr = np.array(self.ev_rec, dtype=np.int64).reshape(n, 4)
        ev, m, rc, fc = arr.T
        isb = np.zeros(n, dtype=bool)
        bev = set()
        if self.burst_n > 0 and n:
            u, inv = np.unique(ev, return_inverse=True)
            cnt = np.bincount(inv, weights=rc)
            isb = cnt[inv] >= self.burst_n
            bev = {int(e) for e in u[cnt >= self.burst_n]}
        keep = ~isb
        self.den_fix_m = np.bincount(m[keep], minlength=nm).astype(np.int64)
        self.num_fix_m = np.bincount(m[keep], weights=fc[keep],
                                     minlength=nm).astype(np.int64)
        self.n_burst_events = len(bev)
        self.n_burst_den = int(isb.sum())
        self.n_burst_cand = int(rc[isb].sum())
        lost = (rc == 1) & (fc == 0)
        self.n_win_removed = int((lost & keep).sum())
        self.n_win_removed_burst = int((lost & isb).sum())
        for c in self.cands:
            c["removed_by"] = ("burst" if c["event"] in bev
                               else "window" if not c["sec_keep"] else "")

    def _record(self, r, y_raw, res, lo, hi, Er, kd, sec, sk, t_edge):
        split, ev, trig, m = self.meta[r]
        ped, sig = float(res.ped[r]), float(res.sigma[r])
        y = y_raw.astype(np.float64) - ped
        w = (y_raw < self.overflow).astype(np.float64)
        t0s = res.p_t0[lo:hi].astype(np.float64)
        amps = res.p_amp[lo:hi].astype(np.float64)
        S, wm, chi2, dof, dchi2 = fit_terms(y, w, sig, self.bank, m, t0s, amps)
        self.num_m[m] += 1
        c = dict(
            split=split, event=ev, trig=trig, m=m, name=self.mods[m].name,
            mid=int(self.mods[m].id), x=float(self.mods[m].x),
            y=float(self.mods[m].y),
            raw=y_raw.copy(), ped=ped, sig=sig, clipped=bool(res.clipped[r]),
            t0=t0s, amp=amps, time=res.p_time[lo:hi].astype(np.float64),
            E=np.asarray(Er, dtype=np.float64), kd=kd, sec=sec,
            chi2dof=float(res.p_chi2[lo]), chi2=chi2, dof=dof, dchi2=dchi2,
            edge=[t >= t_edge for t in t0s], sec_keep=list(sk))
        if self.keep_raw and sk and self.n_raw < self.raw_cap:
            self.n_raw += 1
            tf = fine_grid(y_raw.size)
            c.update(w=w, S=S, wm=wm, comp=np.array(
                [amps[k] * self.bank.shape(m, tf - t0s[k], dtype=np.float64)
                 for k in range(len(t0s))]))
        self.cands.append(c)


# --------------------------------------------------------------------- plot --
def in_pre(t, a):
    return a.pre_window[0] <= t <= a.pre_window[1]


def fine_grid(ns):
    return np.linspace(0.0, mf_pulse.DT * (ns - 1), 4 * ns * 4)


def plot_candidate(pdf, c, a, run, ebeam, kind="two_tau"):
    """One page.  a: pre_window, secondary_frac; c carries w, S, wm, comp."""
    y_raw = c["raw"]
    ns = y_raw.size
    ts = mf_pulse.DT * np.arange(ns)
    y = y_raw - c["ped"]
    w, S, wm = c["w"], c["S"], c["wm"]
    model = c["amp"] @ S
    tf = fine_grid(ns)
    kd, sec = c["kd"], c["sec"]
    sec_rm = c.get("sec_rm", [])

    fig = plt.figure(figsize=(11.0, 8.5))
    gs = fig.add_gridspec(2, 2, width_ratios=[2.7, 1.6], height_ratios=[3.0, 1.0],
                          hspace=0.06, wspace=0.05, left=0.07, right=0.98,
                          top=0.92, bottom=0.08)
    ax = fig.add_subplot(gs[0, 0])
    axr = fig.add_subplot(gs[1, 0], sharex=ax)
    axt = fig.add_subplot(gs[:, 1])
    axt.axis("off")

    ax.axvspan(*a.pre_window, color="0.92", zorder=0, label="pre-window "
               f"{a.pre_window[0]:g}-{a.pre_window[1]:g} ns")
    ax.axhline(0.0, color="0.4", ls="--", lw=1.0,
               label=f"pedestal {c['ped']:.1f} ADC (sigma {c['sig']:.2f})")
    ok = w > 0
    ax.errorbar(ts[ok], y[ok], yerr=c["sig"], fmt="o", ms=3.0, color="k",
                elinewidth=0.6, label="data (ped-subtracted)", zorder=3)
    if (~ok).any():
        ax.plot(ts[~ok], y[~ok], "x", color="r", ms=6, label="clipped (masked)")
    colors = ["C0", "C1", "C2", "C3", "C4", "C5"]
    ci = 1
    for k in range(len(c["t0"])):
        comp = c["comp"][k]
        if k == kd:
            col, lab, ls = "C0", f"#{k} dominant", "-"
        elif k in sec:
            col, lab, ls = colors[ci % len(colors)], f"#{k} secondary", "-"
            ci += 1
        else:
            col, ls = "0.55", ":"
            lab = (f"#{k} (onset at window end)" if c["edge"][k]
                   else f"#{k} (window rule)" if k in sec_rm
                   else f"#{k} (< {a.secondary_frac:g} E_dom)")
        ax.plot(tf, comp, ls, color=col, lw=1.4, label=lab)
        ax.axvline(c["time"][k], color=col, lw=0.7, alpha=0.6)
    tot = c["comp"].sum(axis=0)
    ax.plot(tf, tot, "-", color="m", lw=1.0, alpha=0.8, label="sum of components")
    ax.set_ylabel("ADC - pedestal")
    ax.legend(fontsize=7.5, loc="upper right")
    ax.tick_params(labelbottom=False)
    ax.set_title(f"run {run}  split {c['split']:05d}  event {c['event']}  "
                 f"module {c['name']}", fontsize=12)

    resid = y - model
    fit = wm > 0
    axr.axhspan(-c["sig"], c["sig"], color="0.88")
    axr.axhline(0.0, color="0.4", lw=0.8)
    axr.plot(ts[fit], resid[fit], "o", ms=3.0, color="k", label="in fit window")
    axr.plot(ts[~fit & ok], resid[~fit & ok], "o", ms=2.5, mfc="none",
             color="0.5", label="outside fit window")
    axr.set_xlabel("time (ns)")
    axr.set_ylabel("data - fit")
    axr.legend(fontsize=7, loc="upper right")
    axr.set_xlim(ts[0] - 2, ts[-1] + 2)

    Ed = c["E"][kd]
    lines = [
        f"run {run}   split {c['split']:05d}",
        f"event {c['event']}",
        f"module {c['name']}" + ("   (clipped)" if c["clipped"] else ""),
        f"trigger bits 0x{c['trig']:x}",
        f"E_beam {ebeam:.1f} MeV",
        f"chi2/dof (mf joint fit) {c['chi2dof']:.2f}",
        f"chi2 {c['chi2']:.1f} / dof {c['dof']:.0f}",
        "",
        " #  role  t(ns)    A(ADC)  E(MeV)  E/Ed  A/Ad  pre  dchi2",
    ]
    Ad = c["amp"][kd]
    for k in range(len(c["t0"])):
        role = ("dom" if k == kd else "sec" if k in sec
                else "edge" if c["edge"][k] else "win" if k in sec_rm else "low")
        dch = "-" if k == kd else f"{c['dchi2'][k]:.0f}"
        lines.append(
            f"{k:2d}  {role:4s} {c['time'][k]:6.1f} {c['amp'][k]:8.1f} "
            f"{c['E'][k]:7.1f} {c['E'][k] / Ed:5.3f} {c['amp'][k] / Ad:5.3f} "
            f"{'yes' if in_pre(c['time'][k], a) else 'no':>4s} {dch:>6s}")
    lines += [
        "",
        "t = mf onset t0 + t50 (raw, no time offset)",
        f"pre = t in [{a.pre_window[0]:g}, {a.pre_window[1]:g}] ns",
        "dchi2 = chi2(without pulse) - chi2(full)",
        "",
        "E = A I_int(phase) kappa cal_factor" + (" f(A);" if kind == "sampled" else ";"),
        *(["fixed template: E/Ed = A/Ad up to f(A)"] if kind == "sampled" else
          ["fixed template: E/Ed = A/Ad up to the", "I_int phase dependence"]),
    ]
    axt.text(0.0, 1.0, "\n".join(lines), va="top", ha="left", family="monospace",
             fontsize=7.6, transform=axt.transAxes)
    pdf.savefig(fig)
    plt.close(fig)


# --------------------------------------------------------------- scan tasks --
def resolve_templates(spec, run):
    """(path, column, x17): auto = the h50 .dat as the 20261023 scan; otherwise a
    path or a pattern with {runtype} (gate_scan_20261022.runtype_slug)."""
    x17 = run in X17_RUNS
    if spec == "auto":
        path = str(HERE / "output" / ("pulse_templates_x17_combined_h50.dat" if x17
                                      else "pulse_templates_h50.dat"))
    elif "{runtype}" in spec:
        from gate_scan_20261022 import runtype_slug
        path = spec.replace("{runtype}", runtype_slug(run))
    else:
        path = spec
    if path.endswith(".npz"):
        return path, "sampled", x17
    return path, mf_pulse.resolve_column(path, "auto",
                                         prefer="27000" if x17 else "25320"), x17


def write_npz_examples(path, cands):
    d = dict(n=np.array(len(cands)))
    for i, c in enumerate(cands):
        for k in ("raw", "ped", "sig", "t0", "amp", "time", "E", "w", "S", "wm",
                  "comp", "dchi2", "chi2dof", "chi2", "dof", "kd", "m", "split",
                  "event", "trig", "clipped"):
            d[f"c{i}_{k}"] = np.asarray(c[k])
        sk = c.get("sec_keep", c["sec"])
        d[f"c{i}_sec"] = np.asarray(sk, dtype=np.int64)
        d[f"c{i}_secrm"] = np.asarray([k for k in c["sec"] if k not in sk],
                                      dtype=np.int64)
        d[f"c{i}_edge"] = np.asarray(c["edge"], dtype=bool)
        d[f"c{i}_name"] = np.array(c["name"])
    np.savez(path, **d)


def read_npz_examples(path):
    z = np.load(path)
    out = []
    for i in range(int(z["n"])):
        c = {k: z[f"c{i}_{k}"] for k in (
            "raw", "ped", "sig", "t0", "amp", "time", "E", "w", "S", "wm", "comp",
            "dchi2")}
        c.update(sec=[int(s) for s in z[f"c{i}_sec"]], edge=list(z[f"c{i}_edge"]),
                 sec_rm=([int(s) for s in z[f"c{i}_secrm"]]
                         if f"c{i}_secrm" in z.files else []),
                 name=str(z[f"c{i}_name"]), kd=int(z[f"c{i}_kd"]),
                 m=int(z[f"c{i}_m"]), split=int(z[f"c{i}_split"]),
                 event=int(z[f"c{i}_event"]), trig=int(z[f"c{i}_trig"]),
                 clipped=bool(z[f"c{i}_clipped"]),
                 ped=float(z[f"c{i}_ped"]), sig=float(z[f"c{i}_sig"]),
                 chi2dof=float(z[f"c{i}_chi2dof"]), chi2=float(z[f"c{i}_chi2"]),
                 dof=float(z[f"c{i}_dof"]))
        out.append(c)
    return out


def scan_task(job):
    """One (run, split): scan the EVIO file and write counts/cands (+ examples)."""
    a = argparse.Namespace(**job["a"])
    run, split, fpath, first = job["run"], job["split"], job["path"], job["first"]
    sdir = Path(a.out) / "scan"
    t_start = time.monotonic()
    calib = _HC_CALIB_MAP.get(run)
    calib_path = str(_DATABASE_DIR / "calibration" / calib) if calib else ""
    p = C.setup_pipeline(evio_path=fpath, run_num=run, hc_calib_file=calib_path)
    ebeam = float(p.geo.Ebeam)
    mf_path, column, _ = resolve_templates(a.mf_templates, run)
    mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    wcfg = p.wave_ana.cfg
    bank, info = mf_pulse.build_bank(
        [m.name for m in mods], [bool(m.is_pwo4()) for m in mods], mf_path,
        dict(smooth_order=int(wcfg.smooth_order),
             int_tail_ratio=float(wcfg.int_tail_ratio),
             tail_break_n=int(wcfg.tail_break_n)), column=column)
    rec = mf_pulse.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4,
                                   r_model=a.secondary_frac,
                                   ped_nsamples=int(wcfg.ped_nsamples))
    thr_txt = (f"{a.dominant_min_mev:g} MeV" if a.dominant_min_mev > 0
               else f"{a.dominant_frac:g} E_beam")
    print(f"[setup] run {run} split {split:05d}  E_beam {ebeam:.2f} MeV  calib "
          f"{calib_path or str(p.geo.energy_calib_file)}  templates {mf_path} "
          f"(column {column}, {info['kind']}, {len(info['fallback'])} fallback)  "
          f"kappa {a.mf_scale:g}  r_model {a.secondary_frac:g}  E_dom > {thr_txt}",
          flush=True)

    sc = Scanner(a, p, bank, rec, a.mf_scale, ebeam, keep_raw=first)
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
                continue                # SSP0 only, as the benchmark
            sc.collect(split, fe)
            n_ev += 1
            n_blk += 1
            if n_blk >= BLOCK:
                sc.flush()
                n_blk = 0
            if a.max_events and n_ev >= a.max_events:
                stop = True
                break
    sc.flush()
    sc.finalize()
    sc.n_events = n_ev
    n_num = len(sc.cands)
    assert int(sc.den_m.sum()) == sc.n_den and int(sc.num_m.sum()) == n_num
    elapsed = time.monotonic() - t_start

    os.makedirs(sdir, exist_ok=True)
    meta = dict(
        ebeam=ebeam, dominant_frac=a.dominant_frac,
        dominant_min_mev=a.dominant_min_mev, secondary_frac=a.secondary_frac,
        kappa=a.mf_scale, include_clipped=bool(a.include_clipped),
        pre_window=list(a.pre_window), max_events=a.max_events,
        burst_veto=int(a.burst_veto),
        window_rule=list(sc.win) if sc.win is not None else None,
        templates=mf_path, column=column, kind=info["kind"],
        sha1=info["sha1"], n_fallback=len(info["fallback"]),
        bank_run_type=info.get("run_type"),
        calib=os.path.basename(calib_path) or "runinfo default",
        early_ns=float(mf_pulse.DT * rec.ped_n / 2.0), min_sep_ns=mf_pulse.MIN_SEP_NS,
        thr_txt=thr_txt)
    np.savez(sdir / f"counts_{run}_{split:05d}.npz",
             run=run, split=split, n_events=n_ev, n_den=sc.n_den, n_num=n_num,
             n_early=sc.n_early, n_edge=sc.n_edge, n_clip=sc.n_clip, n_fb=sc.n_fb,
             den_m=sc.den_m, num_m=sc.num_m, elapsed=elapsed,
             den_fix_m=sc.den_fix_m, num_fix_m=sc.num_fix_m,
             n_burst_events=sc.n_burst_events, n_burst_den=sc.n_burst_den,
             n_burst_cand=sc.n_burst_cand, n_win_removed=sc.n_win_removed,
             n_win_removed_burst=sc.n_win_removed_burst,
             mod_name=np.array([m.name for m in mods]),
             mod_id=np.array([int(m.id) for m in mods]),
             mod_x=np.array([float(m.x) for m in mods]),
             mod_y=np.array([float(m.y) for m in mods]),
             mod_sx=np.array([float(m.size_x) for m in mods]),
             mod_sy=np.array([float(m.size_y) for m in mods]),
             mod_pwo4=sc.pwo4, meta=np.array(json.dumps(meta)))
    rate = n_num / sc.n_den if sc.n_den else float("nan")
    rate_line = (f"RATE run {run} split {split:05d}: {n_num} qualifying / {sc.n_den} "
                 f"{'' if a.include_clipped else 'unclipped '}PbWO4 channels with "
                 f"E_dom > {thr_txt} = {rate:.4%}  [{n_ev} events]")
    print(rate_line, flush=True)
    n_den_fix, n_num_fix = int(sc.den_fix_m.sum()), int(sc.num_fix_m.sum())
    fix_line = None
    if sc.burst_n > 0 or sc.win is not None:
        fr = n_num_fix / n_den_fix if n_den_fix else float("nan")
        fix_line = (f"FIXED run {run} split {split:05d}: {n_num_fix} / {n_den_fix} = "
                    f"{fr:.4%}  [burst veto {sc.burst_n}: {sc.n_burst_events} events, "
                    f"{sc.n_burst_den} channels, {sc.n_burst_cand} candidates; window "
                    f"rule {sc.win}: {sc.n_win_removed} candidates removed outside "
                    f"bursts, {sc.n_win_removed_burst} inside]")
        print(fix_line, flush=True)
    with open(sdir / f"cands_{run}_{split:05d}.tsv", "w") as fh:
        fh.write(f"# run {run}  split {split}  E_beam {ebeam:.2f} MeV  E_dom > "
                 f"{thr_txt}  secondary_frac {a.secondary_frac}  kappa "
                 f"{a.mf_scale}  templates {os.path.basename(mf_path)}:{column}\n")
        fh.write(f"# {rate_line}\n")
        if fix_line:
            fh.write(f"# {fix_line}\n")
        fh.write("# lists (';'-separated) are per reported secondary; dt = t_sec - t_dom\n")
        fh.write("run\tsplit\tevent\ttrigger_bits\tmodule\tn_pulses\tn_secondary\t"
                 "E_dominant_MeV\tt_dominant_ns\tA_dominant_ADC\tE_secondary_MeV\t"
                 "dt_ns\tE_ratio\tA_ratio\tdchi2\tsec_in_prewindow\tchi2dof\t"
                 "clipped\tplotted\tmodule_id\tx\ty\tremoved_by\n")
        n_kept_raw = 0
        for i, c in enumerate(sc.cands):
            plotted = 0
            if first and not c["removed_by"] and "comp" in c:
                plotted = int(n_kept_raw < a.n_examples)
                n_kept_raw += 1
            kd, sec = c["kd"], c["sec"]
            j = lambda vals: ";".join(vals)   # noqa: E731
            fh.write("\t".join(str(x) for x in (
                run, c["split"], c["event"], f"0x{c['trig']:x}", c["name"],
                len(c["t0"]), len(sec), f"{c['E'][kd]:.1f}",
                f"{c['time'][kd]:.2f}", f"{c['amp'][kd]:.1f}",
                j(f"{c['E'][k]:.1f}" for k in sec),
                j(f"{c['time'][k] - c['time'][kd]:.2f}" for k in sec),
                j(f"{c['E'][k] / c['E'][kd]:.4f}" for k in sec),
                j(f"{c['amp'][k] / c['amp'][kd]:.4f}" for k in sec),
                j(f"{c['dchi2'][k]:.1f}" for k in sec),
                j("1" if in_pre(c["time"][k], a) else "0" for k in sec),
                f"{c['chi2dof']:.3f}", int(c["clipped"]),
                plotted, c["mid"],
                f"{c['x']:.2f}", f"{c['y']:.2f}", c["removed_by"])) + "\n")
    if first:
        kept = [c for c in sc.cands if not c["removed_by"] and "comp" in c]
        write_npz_examples(sdir / f"examples_{run}_{split:05d}.npz", kept[:N_RAW])
    return dict(run=run, split=split, n_events=n_ev, n_den=sc.n_den, n_num=n_num,
                elapsed=elapsed)


# ------------------------------------------------------------------- report --
def load_run(out, run):
    """Sum the run's counts_<run>_*.npz; None when there are none."""
    files = sorted((Path(out) / "scan").glob(f"counts_{run}_*.npz"))
    if not files:
        return None
    r = None
    for f in files:
        z = np.load(f)
        if r is not None:
            m = json.loads(str(z["meta"]))
            if any(m.get(k) != r["meta"].get(k)
                   for k in ("sha1", "kappa", "burst_veto", "window_rule", "thr_txt")):
                print(f"[WARN] {f.name}: template sha1 / kappa / burst veto / window rule / threshold differ from "
                      f"{files[0].name}", flush=True)
        if r is None:
            r = dict(run=run, splits=[], n_events=0, n_den=0, n_num=0, n_early=0,
                     n_edge=0, n_clip=0, n_fb=0, elapsed=0.0, has_fix=True,
                     n_burst_events=0, n_burst_den=0, n_burst_cand=0,
                     n_win_removed=0, n_win_removed_burst=0,
                     den_m=np.zeros_like(z["den_m"]), num_m=np.zeros_like(z["num_m"]),
                     den_fix_m=np.zeros_like(z["den_m"]),
                     num_fix_m=np.zeros_like(z["num_m"]),
                     meta=json.loads(str(z["meta"])),
                     geo={k: z[k] for k in ("mod_name", "mod_id", "mod_x", "mod_y",
                                            "mod_sx", "mod_sy", "mod_pwo4")})
        r["splits"].append(int(z["split"]))
        for k in ("n_events", "n_den", "n_num", "n_early", "n_edge", "n_clip", "n_fb"):
            r[k] += int(z[k])
        r["elapsed"] += float(z["elapsed"])
        r["den_m"] += z["den_m"]
        r["num_m"] += z["num_m"]
        if "den_fix_m" in z.files:       # scans before the pile-up fix have no fixed counts
            r["den_fix_m"] += z["den_fix_m"]
            r["num_fix_m"] += z["num_fix_m"]
            for k in ("n_burst_events", "n_burst_den", "n_burst_cand",
                      "n_win_removed", "n_win_removed_burst"):
                r[k] += int(z[k])
        else:
            r["has_fix"] = False
    return r


def binom_err(num, den):
    num, den = np.asarray(num, float), np.asarray(den, float)
    with np.errstate(divide="ignore", invalid="ignore"):
        rt = num / den
        return rt, np.sqrt(rt * (1.0 - rt) / den)


def run_type_name(run):
    from benchmark_hycal_timing_scan import _run_type
    return _run_type(run)


def rate_line_of(r, ic):
    rate = r["n_num"] / r["n_den"] if r["n_den"] else float("nan")
    return (f"RATE run {r['run']}: {r['n_num']} qualifying / {r['n_den']} "
            f"{'' if ic else 'unclipped '}PbWO4 channels with E_dom > "
            f"{r['meta']['thr_txt']} = {rate:.4%}  [{r['n_events']} events, splits "
            f"{','.join(f'{s:05d}' for s in r['splits'])}]")


def fix_active(r):
    return r["has_fix"] and (r["meta"].get("burst_veto", 0) > 0
                             or r["meta"].get("window_rule") is not None)


def fixed_totals(r):
    return int(r["den_fix_m"].sum()), int(r["num_fix_m"].sum())


def fix_line_of(r):
    den, num = fixed_totals(r)
    rate = num / den if den else float("nan")
    return (f"FIXED run {r['run']}: {num} / {den} = {rate:.4%}  [burst veto "
            f"{r['meta'].get('burst_veto', 0)}: {r['n_burst_events']} events, "
            f"{r['n_burst_den']} channels, {r['n_burst_cand']} candidates removed; "
            f"window rule {fmt_rule(r['meta'].get('window_rule'))}: "
            f"{r['n_win_removed']} candidates removed outside bursts, "
            f"{r['n_win_removed_burst']} inside]")


def fmt_rule(w):
    return "off" if not w else f"dt < {w[0]:g} ns & E/E_dom < {w[1]:g}"


def hot_mask(g, hot):
    return np.array([str(n) in hot for n in g["mod_name"]]) & g["mod_pwo4"].astype(bool)


def sums(den, num, sel):
    d, n = int(den[sel].sum()), int(num[sel].sum())
    return d, n, (n / d if d else float("nan"))


def excl_text(r):
    return (f"excluded: early onset {r['n_early']}, dominant at window end "
            f"{r['n_edge']}, clipped {r['n_clip']}")


def heatmap_figure(r, hot=()):
    g = r["geo"]
    pw = g["mod_pwo4"].astype(bool)
    den, num = ((r["den_fix_m"], r["num_fix_m"]) if r["has_fix"]
                else (r["den_m"], r["num_m"]))
    hm = hot_mask(g, hot)
    rate, _ = binom_err(num, den)
    rate = rate * 100.0
    good = pw & (den >= MIN_DEN)
    vmax = float(np.percentile(rate[good], 98)) if good.any() else 1.0
    if not vmax > 0:
        vmax = float(rate[good].max()) if good.any() and rate[good].max() > 0 else 1.0

    def rects(sel):
        return [Rectangle((g["mod_x"][i] - g["mod_sx"][i] / 2.0,
                           g["mod_y"][i] - g["mod_sy"][i] / 2.0),
                          g["mod_sx"][i], g["mod_sy"][i]) for i in np.nonzero(sel)[0]]

    fig, (axl, axr) = plt.subplots(1, 2, figsize=(15.0, 7.6), layout="constrained")
    pc = PatchCollection(rects(good), cmap="viridis", linewidths=0)
    pc.set_array(rate[good])
    pc.set_clim(0.0, vmax)
    axl.add_collection(PatchCollection(rects(pw & ~good), facecolor="0.82", linewidths=0))
    axl.add_collection(pc)
    fig.colorbar(pc, ax=axl, fraction=0.046, pad=0.02,
                 label=f"pile-up rate (%); vmax = 98th percentile")
    withden = pw & (den > 0)
    pd = PatchCollection(rects(withden), cmap="plasma", linewidths=0,
                         norm=LogNorm(vmin=1.0, vmax=max(float(den[pw].max()), 2.0)))
    pd.set_array(den[withden].astype(float))
    axr.add_collection(PatchCollection(rects(pw & ~withden), facecolor="0.82",
                                       linewidths=0))
    axr.add_collection(pd)
    if hm.any():
        for ax in (axl, axr):
            ax.add_collection(PatchCollection(rects(hm), facecolor="none",
                                              edgecolor="red", linewidths=1.3))
    fig.colorbar(pd, ax=axr, fraction=0.046, pad=0.02, label="denominator channels")
    xl = (float((g["mod_x"] - g["mod_sx"] / 2)[pw].min()),
          float((g["mod_x"] + g["mod_sx"] / 2)[pw].max()))
    yl = (float((g["mod_y"] - g["mod_sy"] / 2)[pw].min()),
          float((g["mod_y"] + g["mod_sy"] / 2)[pw].max()))
    for ax, t in ((axl, f"rate in %  (grey: denominator < {MIN_DEN}"
                        + (f"; red outline: hot {','.join(hot)})" if hm.any() else ")")),
                  (axr, "denominator channels per module (log)")):
        ax.set_xlim(*xl)
        ax.set_ylim(*yl)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(t, fontsize=10)
    rate_g = r["n_num"] / r["n_den"] if r["n_den"] else float("nan")
    if r["has_fix"]:
        fd, fn = fixed_totals(r)
        fix_txt = (f"\nfixed {fn / fd if fd else float('nan'):.3%} "
                   f"(N_den {fd}, N_num {fn}; burst veto "
                   f"{r['meta'].get('burst_veto', 0)}: {r['n_burst_events']} events / "
                   f"{r['n_burst_den']} channels vetoed; window rule "
                   f"{fmt_rule(r['meta'].get('window_rule'))}: "
                   f"{r['n_win_removed']} cand. removed)")
    else:
        fix_txt = ""
    fig.suptitle(
        f"run {r['run']} ({run_type_name(r['run'])})  raw pile-up rate {rate_g:.3%}  "
        f"N_den {r['n_den']}  N_num {r['n_num']}{fix_txt}\n{len(r['splits'])} split(s), "
        f"{r['n_events']} events   {excl_text(r)}   E_dom > {r['meta']['thr_txt']}, "
        f"secondary > {r['meta']['secondary_frac']:g} E_dom, "
        f"templates {os.path.basename(r['meta']['templates'])}", fontsize=10.5)
    return fig


def cover_page(pdf, run, rate_line, info_lines):
    fig = plt.figure(figsize=(11.0, 8.5))
    fig.text(0.06, 0.92, f"Run {run}: single-module pile-up examples (matched filter)",
             fontsize=15, weight="bold")
    fig.text(0.06, 0.86, "\n".join(info_lines + ["", rate_line]), va="top",
             family="monospace", fontsize=9)
    pdf.savefig(fig)
    plt.close(fig)


def info_lines_of(r, n_plotted):
    mt = r["meta"]
    tmpl = (f"{os.path.basename(mt['templates'])} (sampled, run type "
            f"{mt['bank_run_type']}, {mt['n_fallback']} fallback modules)"
            if mt["kind"] == "sampled" else
            f"{os.path.basename(mt['templates'])} column {mt['column']}")
    return [
        f"run type {run_type_name(r['run'])}   E_beam {mt['ebeam']:.2f} MeV   calib {mt['calib']}",
        f"templates {tmpl}   kappa {mt['kappa']:g}",
        f"mf: nsigma 5, min_amp 5 ADC, kmax 4, r_model {mt['secondary_frac']:g}, "
        f"MIN_SEP {mt['min_sep_ns']:g} ns",
        f"selection: >= 2 mf pulses, E_dom > {mt['thr_txt']}, "
        f">= 1 secondary with E/E_dom > {mt['secondary_frac']:g}",
        f"trigger: bits == 0x100   pre-window {mt['pre_window'][0]:g}-{mt['pre_window'][1]:g} ns",
        f"splits scanned: {','.join(f'{s:05d}' for s in r['splits'])}   "
        f"events {r['n_events']}   max-events {mt['max_events']}",
        f"channels with E_dom above threshold and a pulse onset < {mt['early_ns']:g} ns "
        f"(pedestal window; not counted): {r['n_early']}",
        f"channels with E_dom above threshold but dominant onset at the window end "
        f"(not counted): {r['n_edge']}",
        (f"clipped channels skipped: {r['n_clip']}" if not mt["include_clipped"] else
         f"clipped channels kept (WaveAnalyzer fallback ones skipped: {r['n_fb']})"),
        f"plotted: {n_plotted} candidates (first kept ones of the run's first split"
        + ("; burst-vetoed and window-rule candidates not shown" if fix_active(r) else "")
        + f"); all candidates of all splits (raw): {r['n_num']}",
        *([f"pile-up fix: burst veto {mt.get('burst_veto', 0)} (events with >= N candidate "
           f"channels removed), window rule {fmt_rule(mt.get('window_rule'))}",
           f"  vetoed: {r['n_burst_events']} events, {r['n_burst_den']} channels, "
           f"{r['n_burst_cand']} candidates; window rule removed {r['n_win_removed']} "
           f"candidates outside bursts ({r['n_win_removed_burst']} inside)"]
          if fix_active(r) else []),
    ]


def report(a, runs):
    out = Path(a.out)
    hot = tuple(h for h in a.hot_channels.split(",") if h)
    rows, hot_rows, loaded = [], [], []
    for run in runs:
        r = load_run(out, run)
        if r is None:
            print(f"[WARN] run {run}: no counts in {out / 'scan'}; skipped", flush=True)
            continue
        loaded.append(r)
    if not loaded:
        sys.exit("[ERROR] nothing to report")
    with PdfPages(out / "pileup_heatmaps_all_runs.pdf") as hpdf, \
            PdfPages(out / "pileup_examples_all_runs.pdf") as epdf:
        for r in loaded:
            run, g = r["run"], r["geo"]
            ic = r["meta"]["include_clipped"]
            rl = rate_line_of(r, ic)
            fl = fix_line_of(r) if fix_active(r) else None
            print(rl, flush=True)
            if fl:
                print(fl, flush=True)
            (out / f"pileup_rate_{run}.txt").write_text(rl + "\n" + (fl + "\n" if fl else ""))
            rate, err = binom_err(r["n_num"], r["n_den"])
            pw = g["mod_pwo4"].astype(bool)
            hm = hot_mask(g, hot)
            nan = float("nan")
            if r["has_fix"]:
                fd, fn = fixed_totals(r)
                frate = fn / fd if fd else nan
                _, _, frate_nh = sums(r["den_fix_m"], r["num_fix_m"], pw & ~hm)
                fixed_cols = (fd, fn, f"{frate:.6g}", f"{frate_nh:.6g}",
                              r["n_burst_events"], r["n_burst_den"], r["n_burst_cand"],
                              r["n_win_removed"], r["n_win_removed_burst"])
            else:
                fixed_cols = (nan,) * 9
            dnh, nnh, rnh = sums(r["den_m"], r["num_m"], pw & ~hm)
            rows.append((run, run_type_name(run),
                         ",".join(f"{s:05d}" for s in r["splits"]), r["n_events"],
                         r["n_den"], r["n_num"], f"{float(rate):.6g}",
                         f"{float(err):.3g}", r["n_early"], r["n_edge"], r["n_clip"],
                         f"{float(rate):.6g}", dnh, nnh, f"{rnh:.6g}") + fixed_cols)
            rm, em = binom_err(r["num_m"], r["den_m"])
            rfm, _ = binom_err(r["num_fix_m"], r["den_fix_m"])
            with open(out / f"pileup_rate_module_{run}.tsv", "w") as fh:
                fh.write("module\tmodule_id\tx\ty\tden\tnum\trate\terr\t"
                         "den_fixed\tnum_fixed\trate_fixed\thot\n")
                for i in np.nonzero(pw)[0]:
                    fh.write(f"{g['mod_name'][i]}\t{g['mod_id'][i]}\t{g['mod_x'][i]:.2f}\t"
                             f"{g['mod_y'][i]:.2f}\t{r['den_m'][i]}\t{r['num_m'][i]}\t"
                             f"{rm[i]:.6g}\t{em[i]:.3g}\t{r['den_fix_m'][i]}\t"
                             f"{r['num_fix_m'][i]}\t{rfm[i]:.6g}\t{int(hm[i])}\n")
            for i in np.nonzero(hm)[0]:
                hot_rows.append((run, g["mod_name"][i], int(r["den_m"][i]),
                                 int(r["num_m"][i]), f"{rm[i]:.6g}",
                                 int(r["den_fix_m"][i]), int(r["num_fix_m"][i]),
                                 f"{rfm[i]:.6g}"))
            fig = heatmap_figure(r, hot)
            fig.savefig(out / f"pileup_heatmap_{run}.png", dpi=130)
            hpdf.savefig(fig)
            plt.close(fig)
            ex = sorted((out / "scan").glob(f"examples_{run}_*.npz"))
            cands = read_npz_examples(ex[0])[:a.n_examples] if ex else []
            cover_page(epdf, run, rl + ("\n" + fl if fl else ""),
                       info_lines_of(r, len(cands)))
            st = SimpleNamespace(pre_window=tuple(r["meta"]["pre_window"]),
                                 secondary_frac=r["meta"]["secondary_frac"])
            for c in cands:
                plot_candidate(epdf, c, st, run, r["meta"]["ebeam"], r["meta"]["kind"])
    with open(out / "pileup_rate_summary.tsv", "w") as fh:
        fh.write("run\trun_type\tsplits\tn_events\tN_den\tN_num\trate\terror\t"
                 "n_early\tn_edge\tn_clipped\trate_raw\tN_den_raw_nohot\t"
                 "N_num_raw_nohot\trate_raw_nohot\tN_den_fixed\tN_num_fixed\t"
                 "rate_fixed\trate_fixed_nohot\tn_burst_events\tn_burst_channels\t"
                 "n_burst_cands\tn_window_removed\tn_window_removed_in_burst\n")
        for row in rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    with open(out / "pileup_hot_channels.tsv", "w") as fh:
        fh.write("run\tmodule\tden\tnum\trate\tden_fixed\tnum_fixed\trate_fixed\n")
        for row in hot_rows:
            fh.write("\t".join(str(x) for x in row) + "\n")
    print(f"[out] {out}/pileup_rate_summary.tsv, pileup_heatmaps_all_runs.pdf, "
          f"pileup_examples_all_runs.pdf, pileup_heatmap_<run>.png, "
          f"pileup_rate_module_<run>.tsv, pileup_hot_channels.tsv", flush=True)


# --------------------------------------------------------------------- main --
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("--runs", type=int, nargs="+", default=None,
                    help="runs to scan (default: the 14 PbWO4 runs)")
    ap.add_argument("--run", type=int, default=None, help="alias: --runs RUN")
    ap.add_argument("--splits", nargs="+", default=["all"],
                    help="'all' (default) or split numbers; a run lacking a "
                         "requested split skips it (the first listed split of a "
                         "run supplies its examples)")
    ap.add_argument("--jobs", type=int, default=1,
                    help="processes, one task per (run, split)")
    ap.add_argument("--n-examples", type=int, default=10)
    ap.add_argument("--dominant-frac", type=float, default=0.5,
                    help="dominant pulse E > frac * E_beam (default 0.5)")
    ap.add_argument("--dominant-min-mev", type=float, default=0.0,
                    help="dominant pulse E > this many MeV (default 0 = off, "
                         "use --dominant-frac)")
    ap.add_argument("--secondary-frac", type=float, default=0.10,
                    help="secondary E / E_dom > frac; also the mf r_model "
                         "(default 0.10)")
    ap.add_argument("--mf-templates", default="auto",
                    help="auto (h50 .dat as before), a .dat/.npz path, or a "
                         "pattern with {runtype} (e.g. .../{runtype}/"
                         "templates_{runtype}.npz)")
    ap.add_argument("--out", default=str(HERE / "output" / "pileup_20261024"))
    ap.add_argument("--evio-dir", default=_DEFAULT_EVIO_DIR)
    ap.add_argument("--max-events", type=int, default=0,
                    help="cap on triggered events per split (0 = full split)")
    ap.add_argument("--pre-window", type=float, nargs=2, default=(100.0, 200.0))
    ap.add_argument("--mf-scale", type=float, default=1.010, help="kappa")
    ap.add_argument("--include-clipped", action="store_true",
                    help="keep clipped waveforms (default: skipped, mf amplitude "
                         "extrapolated over masked samples)")
    ap.add_argument("--burst-veto", type=int, default=0,
                    help="remove events with >= N candidate channels from numerator "
                         "and denominator (default 0 = off; production 10)")
    ap.add_argument("--window-rule", default=None, metavar="DT_MAX:RATIO",
                    help="a secondary with 0 <= dt < DT_MAX ns and E/E_dom < RATIO "
                         "is not pile-up (default off; production 40:0.16)")
    ap.add_argument("--hot-channels", default=HOT_DEFAULT,
                    help="report: modules listed separately, rate also given "
                         f"without them (default {HOT_DEFAULT})")
    ap.add_argument("--scan-only", action="store_true", help="skip the report stage")
    ap.add_argument("--report-only", action="store_true",
                    help="only the report stage, from the counts in --out/scan")
    a = ap.parse_args(argv)
    if a.scan_only and a.report_only:
        sys.exit("[ERROR] --scan-only and --report-only are exclusive")
    parse_window_rule(a.window_rule)
    runs = a.runs or ([a.run] if a.run else ALL_RUNS)
    os.makedirs(a.out, exist_ok=True)

    failed = []
    if not a.report_only:
        jobs = []
        want = None if "all" in a.splits else [int(s) for s in a.splits]
        for run in runs:
            if run in _EXCLUDED_RUNS:
                print(f"[WARN] run {run} is excluded from the scan scope "
                      f"({_EXCLUDED_RUNS[run]}); continuing as requested.", flush=True)
            have = dict(local_splits(a.evio_dir, run))
            if want is None:
                sel = sorted(have)
            else:
                sel = [s for s in want if s in have]
                if len(sel) < len(want):
                    print(f"[WARN] run {run}: splits not found locally: "
                          f"{[s for s in want if s not in have]}", flush=True)
            if not sel:
                print(f"[WARN] run {run}: no splits to scan in {a.evio_dir}", flush=True)
                continue
            for i, s in enumerate(sel):
                jobs.append(dict(a=vars(a), run=run, split=s, path=str(have[s]),
                                 first=(i == 0)))
        print(f"[scan] {len(jobs)} task(s), --jobs {a.jobs}", flush=True)
        t0 = time.monotonic()
        if a.jobs <= 1:
            results = []
            for jb in jobs:
                results.append(scan_task(jb))
        else:
            results = []
            with cf.ProcessPoolExecutor(max_workers=a.jobs) as ex:
                futs = {ex.submit(scan_task, jb): jb for jb in jobs}
                for f in cf.as_completed(futs):
                    jb = futs[f]
                    try:
                        results.append(f.result())
                    except Exception as e:                      # noqa: BLE001
                        failed.append((jb["run"], jb["split"], repr(e)))
                        print(f"[ERROR] run {jb['run']} split {jb['split']:05d}: {e!r}",
                              flush=True)
        for r in sorted(results, key=lambda r: (r["run"], r["split"])):
            print(f"[task] run {r['run']} split {r['split']:05d}: {r['n_events']} "
                  f"events, {r['n_num']} / {r['n_den']}, {r['elapsed']:.0f} s",
                  flush=True)
        print(f"[scan] wall {time.monotonic() - t0:.0f} s", flush=True)
    if not a.scan_only:
        report(a, runs)
    if failed:
        print(f"[ERROR] {len(failed)} task(s) failed: {failed}", flush=True)
        sys.exit(1)


if __name__ == "__main__":
    main()
