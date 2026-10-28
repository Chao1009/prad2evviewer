#!/usr/bin/env python3
"""twopulse_realdata.py -- U5: real-data fake rate and trade-off of the two-pulse chi2 fit.

Plan: plans/plan_ml_pileup_20261027.md (U5).  The production fit of twopulse_fit.py (eps 0.02,
kap 0, search 80 ns, min_sep 1 ns; best-50 bank of the run type, Ta = 2p24GeV_Ta_X17) is run on
the raw waveform of every PbWO4 channel that has an MF pulse 1 (largest-E MF pulse with
100 < t < 200 ns, non-fallback) in the first N selected events of each of the 25 pileup_rates
dump splits.  The dumps (output/pileup_rates_20261027/dump) hold only MF pulse lists, so the
same EVIO split is re-decoded with the same event selection as pileup_rates.process_split
(physics events, (trig & 0xf00) != 0 and (trig & 0x3000000) == 0) and the event numbers of the
first N selected events are asserted equal to the dump's ev_num.  The dump stays the
production-MF comparator.

Stages (scratch caches, resumable):
  decode   split -> raw samples of the wanted channels (one subprocess per split)
  fit      raw -> per-channel one/two-pulse fit table (dchi2, t1, A1, t2, A2, chi2_1, ...)
  report   caches + dumps -> TSVs, report numbers, u5_realdata.pdf (no refit)

    MPLBACKEND=Agg python3 twopulse_realdata.py all --workers 10
    MPLBACKEND=Agg python3 twopulse_realdata.py report          (--report-only: no decode / fit)
    MPLBACKEND=Agg python3 twopulse_realdata.py all --smoke      (1 split, few hundred events)

Output: output/ml_pileup_20261027/U5/ (threshold_scan, mf_agreement, rf_phase, operating_points,
fake_by_amplitude, dt_ratio_<group> TSVs, u5_realdata.pdf, summary.json).  Caches in SCRATCH
(raw/ fit/ mf/ u0_eff.pkl); `report` needs the U0 diagnostics state pickle (scratchpad/diag_state_full.pkl).
Fit options are fixed to the U0 production set; thresholds: U0 0.03 quantiles 107 / 98 / 175 (0.73 / 2.24 / 3.49 GeV H2),
98 for C and Ta.
numpy / scipy / matplotlib; prad2py for the decode stage only.
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import multiprocessing as mp
import os
import re
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

import pileup_rates as pr
import twopulse_fit as tf

HERE = Path(__file__).resolve().parent
OUT = HERE / "output/ml_pileup_20261027/U5"
DUMP = HERE / "output/pileup_rates_20261027/dump"
SCRATCH = Path("/private/tmp/claude-505/-Users-L00338853-work-PRad-prad2evviewer-analysis-pyscripts/"
               "b1de5a72-6ab3-4e2b-8b9b-77e955030d7b/scratchpad/U5")
OPT_PROD = dict(eps=0.02, kap=0.0, search=80.0, grid=0.5, min_sep=1.0)
NSAMP = 100


def log(*a):
    print(*a, flush=True)


# ------------------------------------------------------------ dump channels --
def dump_files():
    fs = sorted(glob.glob(str(DUMP / "*.npz")))
    out = []
    for f in fs:
        m = re.search(r"(\d{6})_(\d{5})\.npz$", f)
        out.append((int(m.group(1)), int(m.group(2)), f))
    return out


def channel_pulses(D):
    """Per MF channel (non-fallback pulses) of one dump: key = ev * 4096 + idx (sorted), pulse 1
    (largest-E pulse with 100 < t < 200 ns; same rule as pileup_rates.channel_table), and the
    pulse range [start, start + n) in the arrays returned (sorted by key, then t)."""
    fb = D["pu_fb"].astype(bool)
    ev, idx = D["pu_ev"][~fb].astype(np.int64), D["pu_idx"][~fb].astype(np.int64)
    t, A, E = (D[k][~fb].astype(float) for k in ("pu_t", "pu_A", "pu_E"))
    key = ev * 4096 + idx
    o = np.lexsort((t, key))
    ev, idx, t, A, E, key = ev[o], idx[o], t[o], A[o], E[o], key[o]
    uk, first, inv, npk = np.unique(key, return_index=True, return_inverse=True, return_counts=True)
    inw = (t > pr.PRE_LO) & (t < pr.PRE_HI)
    Ew = np.where(inw, E, -np.inf)
    order = np.lexsort((-Ew, inv))
    p1 = order[first]
    has1 = inw[p1]
    return dict(key=uk, ev=ev[first], idx=idx[first], first=first, npk=npk, has1=has1,
                p1=p1, t=t, A=A, E=E, t1=np.where(has1, t[p1], np.nan), A1=np.where(has1, A[p1], np.nan),
                E1=np.where(has1, E[p1], np.nan))


# ------------------------------------------------------------ banks / sigma --
_BANKS = {}
_SIG = {}


def get_bank(slug):
    if slug not in _BANKS:
        _BANKS[slug] = tf.load_bank(tf.B50, slug)[0]
    return _BANKS[slug]


def slug_of(template_path):
    return os.path.basename(os.path.dirname(str(template_path)))


def sig_mod(slug):
    """Module pedestal rms: median over the stored best-50 pulse rows (as twopulse_fit.Pool),
    all splits of the slug, floored at mf_pulse.SIGMA_MIN."""
    if slug in _SIG:
        return _SIG[slug]
    import mf_pulse
    cf = SCRATCH / f"sigmod_{slug}.json"
    if cf.exists():
        _SIG[slug] = json.loads(cf.read_text())
        return _SIG[slug]
    names, rms = [], []
    for f in sorted(glob.glob(str(tf.B50 / slug / "pulses" / "*.npz"))):
        z = np.load(f, allow_pickle=False)
        names.append(np.asarray(z["name"]))
        rms.append(np.asarray(z["rms"], float))
    nm, rm = np.concatenate(names), np.concatenate(rms)
    sg = {str(n): max(float(np.median(rm[nm == n])), mf_pulse.SIGMA_MIN) for n in np.unique(nm)}
    cf.write_text(json.dumps(sg))
    _SIG[slug] = sg
    return sg


# ------------------------------------------------------------------- decode --
def n_events_for(cp, n_sel, target, ok_idx_mask=None):
    """Smallest number of selected events whose in-bank pulse-1 channels reach `target`."""
    h = cp["has1"] & (ok_idx_mask[cp["idx"]] if ok_idx_mask is not None else True)
    cum = np.cumsum(np.bincount(cp["ev"][h], minlength=n_sel))
    k = int(np.searchsorted(cum, target)) + 1
    return min(k, n_sel), int(cum[min(k, n_sel) - 1])


def decode_split(run, split, dump_path, target, max_ev_cap=0):
    """Re-decode the split, asserting the dump's event numbers, and keep the raw samples of the
    wanted channels.  Returns the cache dict (also written to SCRATCH/raw)."""
    import _common as C
    from prad2py import dec
    z = np.load(dump_path, allow_pickle=True)
    names = np.asarray(z["module_names"])
    pwo4 = np.asarray(z["module_is_pwo4"], bool)
    tmpl = str(z["templates"])
    slug = slug_of(tmpl)
    bank = get_bank(slug)
    ix = {n: i for i, n in enumerate(bank.names)}
    fbset = set(bank.fallback)
    in_bank = np.array([(str(n) in ix) and (str(n) not in fbset) and bool(pw) for n, pw in zip(names, pwo4)])
    D = {k: z[k] for k in ("pu_ev", "pu_idx", "pu_t", "pu_A", "pu_E", "pu_fb")}
    cp = channel_pulses(D)
    n_sel = int(z["n_sel"])
    nev, n_exp = n_events_for(cp, n_sel, target, in_bank)
    if max_ev_cap:
        nev = min(nev, max_ev_cap)
    sel = cp["has1"] & (cp["ev"] < nev) & in_bank[cp["idx"]]
    want = {int(k): i for i, k in enumerate(cp["key"]) if sel[i]}
    rows = np.nonzero(sel)[0]
    pos = {int(cp["key"][i]): j for j, i in enumerate(rows)}
    Y = np.zeros((len(rows), NSAMP), np.uint16)
    got = np.zeros(len(rows), bool)
    p = C.setup_pipeline(evio_path=os.path.join(pr.EVIO_DIR, f"prad_{run:06d}.evio.{split:05d}"),
                         hc_calib_file=str(z["calib_file"]))
    path = os.path.join(pr.EVIO_DIR, f"prad_{run:06d}.evio.{split:05d}")
    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    if ch.open_auto(path) != dec.Status.success:
        sys.exit(f"[ERROR] cannot open {path}")
    ev_num = []
    n_evt = 0
    done = False
    t0 = time.time()
    while not done and ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            e = d["event"]
            tb = int(e.info.trigger_bits)
            if not ((tb & pr.TRIG_MASK) and not (tb & pr.TRIG_EXCL)):
                continue
            ev = n_evt
            n_evt += 1
            ev_num.append(int(e.info.event_number))
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
                        if mod is None or not mod.is_hycal() or not pwo4[mod.index]:
                            continue
                        j = pos.get(ev * 4096 + mod.index)
                        if j is None:
                            continue
                        cd = slot.channel(c)
                        assert cd.nsamples == NSAMP, (run, split, ev, mod.index, cd.nsamples)
                        Y[j] = cd.samples
                        got[j] = True
            if n_evt >= nev:
                done = True
                break
    ev_num = np.asarray(ev_num, np.int64)
    assert len(ev_num) == nev and np.array_equal(ev_num, np.asarray(z["ev_num"])[:nev]), \
        f"event alignment failed {run}/{split}: decoded {len(ev_num)} vs dump {nev}"
    assert got.all(), f"{(~got).sum()} wanted channels not found in {run}/{split}"
    res = dict(run=run, split=split, group=str(z["group"]), slug=slug, templates=tmpl, n_ev=nev,
               n_sel=n_sel, n_chan=len(rows), ev=cp["ev"][rows], idx=cp["idx"][rows],
               mf_t1=cp["t1"][rows], mf_A1=cp["A1"][rows], mf_E1=cp["E1"][rows],
               mf_first=cp["first"][rows], mf_npk=cp["npk"][rows], raw=Y, ev_num=ev_num,
               n_hasp1_all=int(cp["has1"][cp["ev"] < nev].sum()), n_not_in_bank=int((cp["has1"] & (cp["ev"] < nev)).sum() - len(rows)),
               decode_s=time.time() - t0)
    return res


def save_raw(res):
    d = SCRATCH / "raw"
    d.mkdir(parents=True, exist_ok=True)
    f = d / f"{res['run']:06d}_{res['split']:05d}.npz"
    np.savez(str(f) + ".tmp.npz", **res)
    os.replace(str(f) + ".tmp.npz", f)
    return f


# ---------------------------------------------------------------------- fit --
def _fit_worker(a):
    slug, Y, m, sigp = a
    b = _BANKS[slug]
    out = np.empty((len(m), tf.NFIT))
    for i in range(len(m)):
        mi = int(m[i])
        out[i] = tf.fit_wave(Y[i], b.lut[mi].astype(np.float64), b.dlut[mi].astype(np.float64),
                             b.T_end, b.t_pk[mi], sigp[i], OPT_PROD)
    return out


def prep_fit(raw_file):
    """Pedestal-subtracted waveforms, bank module index, sigma_ped of a raw cache."""
    z = np.load(raw_file, allow_pickle=True)
    slug = str(z["slug"])
    bank = get_bank(slug)
    zz = np.load(sorted(DUMP.glob(f"{int(z['run']):06d}_{int(z['split']):05d}.npz"))[0], allow_pickle=True)
    mn = [str(x) for x in zz["module_names"]]
    ix = {n: i for i, n in enumerate(bank.names)}
    midx = np.array([ix[mn[int(i)]] for i in z["idx"]], np.int64)
    sg = sig_mod(slug)
    sigp = np.array([sg[mn[int(i)]] for i in z["idx"]])
    y, ped = tf.subtract_ped(bank, z["raw"])
    return z, slug, y, ped, midx, sigp


def fit_split(raw_file, pool, chunk=40):
    z, slug, y, ped, midx, sigp = prep_fit(raw_file)
    n = len(midx)
    t0 = time.time()
    tasks = [(slug, y[i:i + chunk], midx[i:i + chunk], sigp[i:i + chunk]) for i in range(0, n, chunk)]
    R = np.concatenate(pool.map(_fit_worker, tasks, chunksize=1)) if tasks else np.zeros((0, tf.NFIT))
    t50 = get_bank(slug).t50[midx]
    out = dict(run=int(z["run"]), split=int(z["split"]), group=str(z["group"]), slug=slug,
               n_ev=int(z["n_ev"]), n_sel=int(z["n_sel"]), ev=z["ev"], idx=z["idx"], R=R, ped=ped, sigp=sigp,
               midx=midx, t50=t50, ymax=y.max(1), mf_t1=z["mf_t1"], mf_A1=z["mf_A1"], mf_E1=z["mf_E1"],
               mf_first=z["mf_first"], mf_npk=z["mf_npk"], fit_s=time.time() - t0,
               opt=np.array([OPT_PROD[k] for k in sorted(OPT_PROD)]))
    return out


# ============================================================================ report ==
ORDER = ["0p73", "2p24", "3p49", "c2p24", "x17_25800", "x17_26138", "x17_26230"]
LAB = {"0p73": "0.73 GeV H2", "2p24": "2.24 GeV H2", "3p49": "3.49 GeV H2", "c2p24": "2.24 GeV C",
       "x17_25800": "Ta 25800", "x17_26138": "Ta 26138", "x17_26230": "Ta 26230"}
EP = ("0p73", "2p24", "3p49")
TA = ("x17_25800", "x17_26138", "x17_26230")
U0_NAME = {"0p73": "0.73 GeV H2", "2p24": "2.24 GeV H2", "3p49": "3.49 GeV H2"}
PROD_THR = {"0p73": 107.24, "2p24": 97.87, "3p49": 174.74, "c2p24": 97.87,
            "x17_25800": 97.87, "x17_26138": 97.87, "x17_26230": 97.87}   # U0 0.03 quantile; C and Ta use 98
U0_01 = {"0p73": 240.46, "2p24": 648.74, "3p49": 947.88}                    # U0 0.1 % thresholds
THR_NAMED = [25, 50, 100, 175, 300, 500, 1000, 2000, 5000, 10000]
THR_FINE = sorted(set(int(round(x)) for x in np.geomspace(20, 20000, 49)) | set(THR_NAMED)
                  | {98, 107, 175, 240, 649, 948})
VETO_SEP, VETO_RATIO = 2.0, 0.05
DT_BINS = [(0, 2), (2, 8), (8, 16), (16, 60), (60, 80)]
R_BINS = [(0, 0.3), (0.3, 1.0), (1.0, 1e9)]
EFF_CELLS = [(8, 1.0), (8, 0.3), (12, 1.0), (12, 0.3)]
C_CLOSE = 16.0
REACH = 80.0


def wilson_err(k, n):
    p = k / max(n, 1)
    return math.sqrt(max(p * (1 - p), 0) / max(n, 1))


def fold(x, P):
    return np.mod(x + P / 2, P) - P / 2


def mf_cache(run, split, z):
    """MF pulse lists of the sampled channels (flat arrays, channel index = row of the fit cache)."""
    f = SCRATCH / "mf" / f"{run:06d}_{split:05d}.npz"
    if f.exists():
        return dict(np.load(f))
    zz = np.load(DUMP / f"{run:06d}_{split:05d}.npz", allow_pickle=True)
    D = {k: zz[k] for k in ("pu_ev", "pu_idx", "pu_t", "pu_A", "pu_E", "pu_fb")}
    cp = channel_pulses(D)
    first, npk = z["mf_first"], z["mf_npk"]
    ch = np.repeat(np.arange(len(first)), npk)
    off = np.arange(int(npk.sum())) - np.repeat(np.cumsum(npk) - npk, npk)
    rows = np.repeat(first, npk) + off
    assert np.allclose(cp["A1"][np.searchsorted(cp["key"], z["ev"].astype(np.int64) * 4096 + z["idx"])], z["mf_A1"])
    out = dict(ch=ch, t=cp["t"][rows], A=cp["A"][rows], E=cp["E"][rows])
    f.parent.mkdir(parents=True, exist_ok=True)
    np.savez(f, **out)
    return out


class Grp:
    """Concatenated fit table + MF extras of one run group, with derived pair variables."""

    def __init__(self, g, files):
        self.g = g
        Z = [np.load(f, allow_pickle=True) for f in files]
        self.splits = [(int(z["run"]), int(z["split"]), len(z["R"]), int(z["n_ev"])) for z in Z]
        self.n_outside = sum(int(np.load(SCRATCH / "raw" / f"{int(z['run']):06d}_{int(z['split']):05d}.npz", allow_pickle=True)["n_not_in_bank"]) for z in Z)
        self.R = np.concatenate([z["R"] for z in Z])
        self.t50 = np.concatenate([z["t50"] for z in Z])
        self.mf_t1 = np.concatenate([z["mf_t1"] for z in Z])
        self.mf_A1 = np.concatenate([z["mf_A1"] for z in Z])
        self.ymax = np.concatenate([z["ymax"] for z in Z])
        self.ped = np.concatenate([z["ped"] for z in Z])
        self.midx = np.concatenate([z["midx"] for z in Z])
        self.sigp = np.concatenate([z["sigp"] for z in Z])
        self.slug = str(Z[0]["slug"])
        self.sid = np.concatenate([np.full(len(z["R"]), int(z["run"]) * 100000 + int(z["split"])) for z in Z])
        self.row = np.concatenate([np.arange(len(z["R"])) for z in Z])
        self.ev = np.concatenate([z["ev"] for z in Z])
        self.idx = np.concatenate([z["idx"] for z in Z])
        self.n = len(self.R)
        self.n_ev = sum(s[3] for s in self.splits)
        ch, t, A, off = [], [], [], 0
        for z in Z:
            m = mf_cache(int(z["run"]), int(z["split"]), z)
            ch.append(m["ch"] + off)
            t.append(m["t"])
            A.append(m["A"])
            off += len(z["R"])
        ch, t, A = np.concatenate(ch), np.concatenate(t), np.concatenate(A)
        isp1 = (t == self.mf_t1[ch]) & (A == self.mf_A1[ch])
        ex = ~isp1
        # exactly one pulse 1 per channel
        assert np.bincount(ch[isp1], minlength=self.n).min() >= 1
        self.xch, self.xdt = ch[ex], t[ex] - self.mf_t1[ch[ex]]
        self.xratio = A[ex] / self.mf_A1[ch[ex]]
        self.xt = t[ex]
        R = self.R
        self.valid2 = R[:, 13] == 1
        self.dchi = R[:, 12]
        tt1, tt2 = R[:, 9] + self.t50, R[:, 11] + self.t50
        a1, a2 = R[:, 8], R[:, 10]
        first = np.abs(tt1 - self.mf_t1) <= np.abs(tt2 - self.mf_t1)
        self.pt = np.where(first, tt1, tt2)
        self.pa = np.where(first, a1, a2)
        self.ot = np.where(first, tt2, tt1)
        self.oa = np.where(first, a2, a1)
        self.dt = self.ot - self.pt
        self.ratio = self.oa / np.maximum(self.pa, 1e-9)
        self.on_bank = np.where(first, R[:, 11], R[:, 9])          # secondary onset on the bank axis
        self.p1_bank = np.where(first, R[:, 9], R[:, 11])
        with np.errstate(invalid="ignore"):
            self.degen = (np.abs(tt2 - tt1) < VETO_SEP) | (np.minimum(a1, a2) / np.maximum(np.maximum(a1, a2), 1e-9) < VETO_RATIO)
        self.a100 = self.mf_A1 >= 100

    def det(self, thr):
        return self.valid2 & ~self.degen & (self.dchi > thr)

    def veto(self, thr):
        return self.valid2 & self.degen & (self.dchi > thr)


def load_groups():
    by = {g: [] for g in ORDER}
    for f in sorted((SCRATCH / "fit").glob("*.npz")):
        by[str(np.load(f, allow_pickle=True)["group"])].append(f)
    return {g: Grp(g, by[g]) for g in ORDER if by[g]}


def u0_eff():
    """U0 synthetic efficiency curves (sorted dchi of matched events, n) from the diagnostics state pickle."""
    cf = SCRATCH / "u0_eff.pkl"
    if cf.exists():
        import pickle
        return pickle.load(open(cf, "rb"))
    import pickle
    from twopulse_diagnostics import match_t
    SS = pickle.load(open(SCRATCH.parent / "diag_state_full.pkl", "rb"))
    out = {}
    for rt, S in SS.items():
        evs = S["Rb"][S["pool"]["parity"] == 1, 12]
        out[(rt, "fake_eval")] = np.sort(evs)
        for k, c in S["cells"].items():
            if k[0] != "real2":
                continue
            m = match_t(c["R"], c["ev"]["t1"], c["ev"]["t2"])
            out[(rt, k[1], k[2])] = (np.sort(c["R"][m, 12]), len(m))
        out[(rt, "thr")] = S["thr"]
    pickle.dump(out, open(cf, "wb"))
    return out


def eff_at(U, rt, dt, r, thr):
    """U0 real2 efficiency at thr; the run type's own cell if the U0 state has it, else the 2.24 GeV H2 cell
    (U0 state: dt 8 for 0.73 / 3.49 GeV H2, the full grid for 2.24 GeV H2; C and Ta use 2.24 GeV H2)."""
    k = (rt, dt, r)
    if k not in U:
        k = ("2.24 GeV H2", dt, r)
    if k not in U:
        return np.nan
    d, n = U[k]
    return (len(d) - np.searchsorted(d, thr, side="right")) / n


def expected_close(G, msk=None):
    """Genuine close-pair rate per fired channel from the MF density of ratio > 0.3 extra pulses
    (both signs): bands 24-60 and 16-60 ns, scaled to |dt| < 16 (32 ns) and to |dt| < 80 (measured
    24/16-80 + band density x the rest)."""
    out = {}
    ad, r = np.abs(G.xdt), G.xratio
    sel = np.ones(len(ad), bool) if msk is None else msk[G.xch]
    n = G.n if msk is None else int(msk.sum())
    for tag, lo in (("a", 24.0), ("b", 16.0)):
        m = (ad > lo) & (ad < 60) & (r > 0.3) & sel
        dens = m.sum() / (n * 2 * (60 - lo))
        mw = (ad > lo) & (ad < REACH) & (r > 0.3) & sel
        out[tag] = dict(dens=dens, close=dens * 2 * C_CLOSE, win=mw.sum() / n + dens * 2 * lo, n_band=int(m.sum()))
    return out


def scan_rows(G):
    rows = []
    ad = np.abs(G.dt)
    for sample, msk in (("all", np.ones(G.n, bool)), ("A1ge100", G.a100)):
        n = int(msk.sum())
        for thr in sorted(set(THR_FINE) | {round(PROD_THR[G.g], 2)}):
            d = G.det(thr) & msk
            c = d & (ad < C_CLOSE)
            c3 = c & (G.ratio > 0.3)
            w = d & (ad <= REACH)
            w3 = w & (G.ratio > 0.3)
            rows.append(dict(group=G.g, sample=sample, thr=thr, n_ch=n, n_det=int(d.sum()), n_close=int(c.sum()),
                             n_neg=int((c & (G.dt < 0)).sum()), n_pos=int((c & (G.dt > 0)).sum()),
                             n_close_r03=int(c3.sum()), n_win=int(w.sum()), n_win_r03=int(w3.sum()),
                             n_veto=int((G.veto(thr) & msk).sum()), n2raw=int(((G.dchi > thr) & G.valid2 & msk).sum())))
    return rows


def rate(k, n):
    return k / max(n, 1)


def op_threshold(rows_g, exp_lo, key="n_close", sample="all"):
    """Lowest scanned threshold with fake < real by both estimators: det < 2 exp_lo (fake = det - exp) and
    the symmetric fake estimate (N(0<dt<16) - N(-16<dt<0)) / n < exp_lo, in the sample."""
    for r in sorted((r for r in rows_g if r["sample"] == sample), key=lambda r: r["thr"]):
        if rate(r[key], r["n_ch"]) < 2 * exp_lo and rate(r["n_pos"] - r["n_neg"], r["n_ch"]) < exp_lo:
            return r["thr"]
    return None


def table_cells(G, thr_list, msk):
    """Fit and MF extra-pulse counts per (|dt| bin, ratio bin) over the channels in msk."""
    ad_f, ad_m = np.abs(G.dt), np.abs(G.xdt)
    xm = msk[G.xch]
    out = []
    for (dl, dh) in DT_BINS + [(0, 16)]:
        for (rl, rh) in R_BINS + [(0, 1e9)]:
            mm = (ad_m >= dl) & (ad_m < dh) & (G.xratio >= rl) & (G.xratio < rh) & xm
            row = dict(dt_lo=dl, dt_hi=dh, r_lo=rl, r_hi=rh, mf=int(mm.sum()))
            for nm, thr in thr_list:
                d = G.det(thr) & msk & (ad_f >= dl) & (ad_f < dh) & (G.ratio >= rl) & (G.ratio < rh)
                row[nm] = int(d.sum())
            out.append(row)
    return out


def agreement(G, thr, msk):
    """(a) fit detections vs MF extra pulses in the MF-complete band, +-4 ns (channels in msk)."""
    d = G.det(thr) & msk
    ad = np.abs(G.dt)
    cand = d & (ad > 24) & (ad < REACH) & (G.ratio > 0.3)
    # MF extra within 4 ns of the fitted secondary
    near = np.abs(G.xt - G.ot[G.xch]) <= 4.0
    has_mf = np.zeros(G.n, bool)
    has_mf[G.xch[near]] = True
    k1, n1 = int((cand & has_mf).sum()), int(cand.sum())
    mx = (np.abs(G.xdt) > 24) & (np.abs(G.xdt) < 76) & (G.xratio > 0.3) & msk[G.xch]
    c = G.xch[mx]
    hit = d[c] & (np.abs(G.xt[mx] - G.ot[c]) <= 4.0)
    k2, n2 = int(hit.sum()), int(mx.sum())
    return dict(n_fit=n1, n_fit_in_mf=k1, n_mf=n2, n_mf_in_fit=k2)


def comb_rows(G, thr, fref):
    rows = []
    d = G.det(thr)
    ad = np.abs(G.dt)
    for lab, lo, hi in (("<16", 0, 16), ("16-80", 16, REACH + 1e-6), ("all", 0, REACH + 1e-6)):
        m = d & (ad >= lo) & (ad < hi)
        n = int(m.sum())
        x = G.dt[m]
        kon = int((np.abs(fold(x, pr.RF_P)) < 0.5).sum())
        kgr = int((np.abs(fold(x, pr.GRID_P)) < 0.5).sum())
        f = kon / max(n, 1)
        fl = 1 - (f - 0.25) / (fref - 0.25) if n else np.nan
        fle = math.sqrt(max(f * (1 - f), 0) / max(n, 1)) / (fref - 0.25) if n else np.nan
        ph2 = np.mod(G.on_bank[m], 4.0)
        ph1 = np.mod(G.p1_bank[m], 4.0)
        h2 = np.histogram(ph2, bins=8, range=(0, 4))[0]
        h1 = np.histogram(ph1, bins=8, range=(0, 4))[0]
        ex = max(n, 1) / 8
        rows.append(dict(group=G.g, thr=thr, subset=lab, n=n, on_comb=kon, f_on=f, f_on_err=math.sqrt(max(f * (1 - f), 0) / max(n, 1)),
                         on_grid=kgr, f_grid=kgr / max(n, 1), f_flat=fl, f_flat_err=fle, f_ref=fref,
                         phase2_chi2_ndf=float(((h2 - ex) ** 2 / ex).sum() / 7) if n else np.nan,
                         phase1_chi2_ndf=float(((h1 - ex) ** 2 / ex).sum() / 7) if n else np.nan,
                         phase2_hist=" ".join(map(str, h2))))
    return rows


def f_ref_ta(GG):
    """On-comb fraction of MF-complete pairs (24 < |dt| < 80, ratio > 0.3) over the Ta sample, pooled / per run."""
    k = n = 0
    per = {}
    for g in TA:
        if g not in GG:
            continue
        G = GG[g]
        m = (np.abs(G.xdt) > 24) & (np.abs(G.xdt) < REACH) & (G.xratio > 0.3)
        kk = int((np.abs(fold(G.xdt[m], pr.RF_P)) < 0.5).sum())
        per[g] = (kk, int(m.sum()))
        k += kk
        n += int(m.sum())
    return k / max(n, 1), n, per


def fmtp(k, n, nd=3, scale=100):
    return f"{scale * k / max(n, 1):.{nd}f}"


def write_tsv(path, rows, cols=None):
    cols = cols or list(rows[0].keys())
    with open(path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            fh.write("\t".join(f"{r[c]:.5g}" if isinstance(r[c], (float, np.floating)) else str(r[c]) for c in cols) + "\n")


def stage_report(a):
    import matplotlib
    matplotlib.use("Agg")
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    GG = load_groups()
    U = u0_eff()
    log(f"[report] loaded {len(GG)} groups, {sum(g.n for g in GG.values())} channels, {time.time() - t0:.0f} s")
    S = dict(groups={})
    EXP = {g: {"all": expected_close(G), "A1ge100": expected_close(G, G.a100)} for g, G in GG.items()}
    ROWS = []
    for g, G in GG.items():
        ROWS += scan_rows(G)
    # --- operating points: lowest thr with det(|dt|<16, all ratios) < 2 x lower expectation
    OP = {}
    for g, G in GG.items():
        e = EXP[g]["all"]
        lo, hi = min(e["a"]["close"], e["b"]["close"]), max(e["a"]["close"], e["b"]["close"])
        e1 = EXP[g]["A1ge100"]
        lo1, hi1 = min(e1["a"]["close"], e1["b"]["close"]), max(e1["a"]["close"], e1["b"]["close"])
        rg = [r for r in ROWS if r["group"] == g]
        OP[g] = dict(thr_all=op_threshold(rg, lo, "n_close"), thr_r03=op_threshold(rg, lo, "n_close_r03"),
                     thr_all_hi=op_threshold(rg, hi, "n_close"), exp_lo=lo, exp_hi=hi,
                     thr_a100=op_threshold(rg, lo1, "n_close", "A1ge100"), exp_lo_a100=lo1, exp_hi_a100=hi1)
    # --- threshold_scan.tsv with expectation and fake and U0 efficiency
    cols = ["group", "label", "sample", "thr", "n_ch", "n_det", "n_close", "n_neg", "n_pos", "n_close_r03", "n_win", "n_win_r03", "n_veto", "n2raw",
            "det_per_ch", "sym_genuine", "sym_fake", "sym_err", "close_per_ch", "close_err", "close_r03_per_ch", "win_per_ch", "win_r03_per_ch",
            "exp_close_lo", "exp_close_hi", "fake_close_hi", "fake_close_lo", "fake_close_r03_hi", "fake_close_r03_lo",
            "exp_win_a", "exp_win_b", "fake_win_r03_hi", "fake_win_r03_lo",
            "eff8_r1", "eff8_r0.3", "eff12_r1", "eff12_r0.3", "u0_fake_eval"]
    out_rows = []
    for r in ROWS:
        g = r["group"]
        e = EXP[g][r["sample"]]
        n = r["n_ch"]
        rt = U0_NAME.get(g, "2.24 GeV H2")
        elo = min(e["a"]["close"], e["b"]["close"])
        ehi = max(e["a"]["close"], e["b"]["close"])
        c = r["n_close"] / n
        c3 = r["n_close_r03"] / n
        w3 = r["n_win_r03"] / n
        row = dict(r, label=LAB[g], sym_genuine=2 * r["n_neg"] / n, sym_fake=(r["n_pos"] - r["n_neg"]) / n,
                   sym_err=math.sqrt(r["n_pos"] + r["n_neg"]) / n, det_per_ch=r["n_det"] / n, close_per_ch=c, close_err=wilson_err(r["n_close"], n),
                   close_r03_per_ch=c3, win_per_ch=r["n_win"] / n, win_r03_per_ch=w3,
                   exp_close_lo=elo, exp_close_hi=ehi, fake_close_hi=c - elo, fake_close_lo=c - ehi,
                   fake_close_r03_hi=c3 - elo, fake_close_r03_lo=c3 - ehi,
                   exp_win_a=e["a"]["win"], exp_win_b=e["b"]["win"],
                   fake_win_r03_hi=w3 - min(e["a"]["win"], e["b"]["win"]), fake_win_r03_lo=w3 - max(e["a"]["win"], e["b"]["win"]))
        ev = U[(rt, "fake_eval")]
        row["u0_fake_eval"] = (len(ev) - np.searchsorted(ev, r["thr"], side="right")) / len(ev) if g in EP else np.nan
        for (dt, rr), nm in zip(EFF_CELLS, ("eff8_r1", "eff8_r0.3", "eff12_r1", "eff12_r0.3")):
            row[nm] = eff_at(U, rt, dt, rr, r["thr"])
        out_rows.append(row)
    write_tsv(OUT / "threshold_scan.tsv", out_rows, cols)
    # --- agreement
    ag_rows = []
    for g, G in GG.items():
        for sample, msk in (("all", np.ones(G.n, bool)), ("A1ge100", G.a100)):
            for thr in sorted({round(PROD_THR[g], 2), 175, 300, 500, 1000, 2000}):
                ag = agreement(G, thr, msk)
                ag_rows.append(dict(group=g, label=LAB[g], sample=sample, thr=thr, **ag,
                                    frac_fit_in_mf=ag["n_fit_in_mf"] / max(ag["n_fit"], 1),
                                    frac_mf_in_fit=ag["n_mf_in_fit"] / max(ag["n_mf"], 1)))
    write_tsv(OUT / "mf_agreement.tsv", ag_rows)
    # --- rf comb
    fref, nref, per = f_ref_ta(GG)
    rf_rows = []
    for g, G in GG.items():
        thr_list = [round(PROD_THR[g], 2)]
        if OP[g]["thr_all"] and abs(OP[g]["thr_all"] - thr_list[0]) > 1:
            thr_list.append(OP[g]["thr_all"])
        for thr in thr_list:
            rf_rows += comb_rows(G, thr, fref)
    write_tsv(OUT / "rf_phase.tsv", rf_rows)
    # --- operating points
    ops = []
    for g, G in GG.items():
        rt = U0_NAME.get(g, "2.24 GeV H2")
        o = OP[g]
        for kind, thr in (("lowest_all", o["thr_all"]), ("lowest_all_A1ge100", o["thr_a100"]), ("lowest_r03", o["thr_r03"]),
                          ("production", round(PROD_THR[g], 2))):
            a1 = kind.endswith("A1ge100")
            msk = G.a100 if a1 else np.ones(G.n, bool)
            nn = int(msk.sum())
            elo, ehi = (o["exp_lo_a100"], o["exp_hi_a100"]) if a1 else (o["exp_lo"], o["exp_hi"])
            row = dict(group=g, label=LAB[g], kind=kind, thr=thr if thr is not None else "none", n_ch=nn,
                       exp_close_lo=elo, exp_close_hi=ehi)
            if thr is not None:
                d = G.det(thr) & msk
                ad = np.abs(G.dt)
                c = d & (ad < C_CLOSE)
                neg, pos = int((c & (G.dt < 0)).sum()), int((c & (G.dt > 0)).sum())
                row.update(close_per_ch=c.sum() / nn, close_r03_per_ch=(c & (G.ratio > 0.3)).sum() / nn,
                           fake_hi=c.sum() / nn - elo, fake_lo=c.sum() / nn - ehi,
                           sym_genuine=2 * neg / nn, sym_fake=(pos - neg) / nn, sym_err=math.sqrt(pos + neg) / nn)
                mf1660 = int(((np.abs(G.xdt) > 16) & (np.abs(G.xdt) < 60) & msk[G.xch]).sum()) / nn
                cfs = c.sum() / nn - (pos - neg) / nn
                row.update(close_fake_sub=cfs, mf_16_60_per_ch=mf1660, real_frac_lt16=cfs / (cfs + mf1660))
                for (dt, rr), nm in zip(EFF_CELLS, ("eff8_r1", "eff8_r0.3", "eff12_r1", "eff12_r0.3")):
                    row[nm] = eff_at(U, rt, dt, rr, thr)
            ops.append(row)
    write_tsv(OUT / "operating_points.tsv", ops, ["group", "label", "kind", "thr", "n_ch", "exp_close_lo", "exp_close_hi", "close_per_ch",
                                                  "close_r03_per_ch", "fake_hi", "fake_lo", "sym_genuine", "sym_fake", "sym_err",
                                                  "close_fake_sub", "mf_16_60_per_ch", "real_frac_lt16", "eff8_r1", "eff8_r0.3", "eff12_r1", "eff12_r0.3"])
    # --- dt x ratio tables (rescan): fit at the operating point and at production, MF on the same channels
    for g, G in GG.items():
        with open(OUT / f"dt_ratio_{g}.tsv", "w") as fh:
            fh.write(f"# {LAB[g]}: {G.n} fired PbWO4 channels with an MF pulse 1 ({int(G.a100.sum())} with A1>=100). Counts of extra pulses per cell, prob = count / channels of the sample; "
                     "fit = two-pulse fit detection (non-degenerate; the other pulse vs the fitted pulse 1; at most one extra pulse per channel) at the operating threshold "
                     "of the sample (lowest thr with detections(|dt|<16) < 2 x MF-extrapolated genuine rate) and at the production threshold; "
                     "mf = production-MF extra pulses on the same channels (|dt| < 80 ns only)\n")
            fh.write("sample\tthr_op\tdt_lo\tdt_hi\tratio_lo\tratio_hi\tn_fit_op\tprob_fit_op\tn_fit_prod\tprob_fit_prod\tn_mf\tprob_mf\tprob_mf_err\n")
            for sample, msk, key in (("all", np.ones(G.n, bool), "thr_all"), ("A1ge100", G.a100, "thr_a100")):
                thr_op = OP[g][key] if OP[g][key] is not None else max(THR_FINE)
                nn = int(msk.sum())
                for c in table_cells(G, [("fit_op", thr_op), ("fit_prod", round(PROD_THR[g], 2))], msk):
                    fh.write(f"{sample}\t{thr_op}\t{c['dt_lo']:g}\t{c['dt_hi']:g}\t{c['r_lo']:g}\t{c['r_hi']:g}\t{c['fit_op']}\t{c['fit_op'] / nn:.3e}\t"
                             f"{c['fit_prod']}\t{c['fit_prod'] / nn:.3e}\t{c['mf']}\t{c['mf'] / nn:.3e}\t{wilson_err(c['mf'], nn):.1e}\n")
    # --- amplitude dependence of the fake rate
    ab = [5, 20, 50, 100, 200, 500, 1000, 2000, 4100]
    amp_rows = []
    for g, G in GG.items():
        ad = np.abs(G.dt)
        for lo, hi in zip(ab[:-1], ab[1:]):
            m = (G.mf_A1 >= lo) & (G.mf_A1 < hi)
            if m.sum() < 200:
                continue
            row = dict(group=g, label=LAB[g], A1_lo=lo, A1_hi=hi, n_ch=int(m.sum()))
            row["thr_prod"] = round(PROD_THR[g], 2)
            for nm, thr in (("prod", round(PROD_THR[g], 2)), ("300", 300), ("500", 500), ("1000", 1000)):
                row[f"close_per_ch_{nm}"] = float((G.det(thr) & m & (ad < C_CLOSE)).sum() / m.sum())
                row[f"veto_per_ch_{nm}"] = float((G.veto(thr) & m).sum() / m.sum())
            row["median_dchi"] = float(np.median(G.dchi[m]))
            amp_rows.append(row)
    cols_a = []
    for r in amp_rows:
        for k in r:
            if k not in cols_a:
                cols_a.append(k)
    with open(OUT / "fake_by_amplitude.tsv", "w") as fh:
        fh.write("\t".join(cols_a) + "\n")
        for r in amp_rows:
            fh.write("\t".join(f"{r.get(c, ''):.5g}" if isinstance(r.get(c, ''), float) else str(r.get(c, '')) for c in cols_a) + "\n")
    off = {}
    for g, G in GG.items():
        m = G.valid2 & ~G.degen & (G.dchi > PROD_THR[g]) & G.a100
        o1 = (G.R[:, 3] + G.t50 - G.mf_t1)[G.a100]
        o2 = (G.pt - G.mf_t1)[m]
        off[g] = dict(single_median=float(np.median(o1)), single_p16=float(np.percentile(o1, 16)), single_p84=float(np.percentile(o1, 84)),
                      pair_p1_median=float(np.median(o2)), pair_p1_p16=float(np.percentile(o2, 16)), pair_p1_p84=float(np.percentile(o2, 84)),
                      ymax_over_A1=float(np.median(G.ymax[G.a100] / G.mf_A1[G.a100])))
    S.update(EXP=EXP, OP=OP, fref=fref, nref=nref, fref_per=per)
    json.dump(dict(EXP=EXP, OP=OP, fref=fref, nref=nref, fref_per=per, offsets=off,
                   groups={g: dict(n=G.n, n_ev=G.n_ev, n_a100=int(G.a100.sum()), n_outside_bank=G.n_outside, splits=G.splits) for g, G in GG.items()},
                   frac_agree=ag_rows), open(OUT / "summary.json", "w"), indent=1, default=float)
    figures(GG, U, EXP, OP, rf_rows, out_rows, fref)
    log(f"[report] done {time.time() - t0:.0f} s")


COL = {"0p73": "#0072B2", "2p24": "#D55E00", "3p49": "#009E73", "c2p24": "#CC79A7",
       "x17_25800": "#444444", "x17_26138": "#E69F00", "x17_26230": "#56B4E9"}


def thr_used(g, OP):
    t = OP[g]["thr_all"]
    return t if t is not None else max(THR_FINE)


def ex_panel(fig, gs_main, gs_res, G, i, thr, title):
    """Raw waveform, 1- and 2-pulse fit and residuals of channel i."""
    bank = get_bank(G.slug)
    m = int(G.midx[i])
    row, T_end = bank.lut[m].astype(np.float64), bank.T_end
    raw = np.load(SCRATCH / "raw" / f"{G.sid[i] // 100000:06d}_{G.sid[i] % 100000:05d}.npz", allow_pickle=True)["raw"][G.row[i]]
    y = raw.astype(float) - G.ped[i]
    R = G.R[i]
    tf_ = np.arange(0.0, 400.0, 0.5)
    sh = lambda t, A: A * tf.tpl(row, None, tf_ - t, T_end)
    shs = lambda t, A: A * tf.tpl(row, None, tf.TS - t, T_end)
    m1, m2 = shs(R[3], R[2]), shs(R[9], R[8]) + shs(R[11], R[10])
    ax = fig.add_subplot(gs_main)
    ar = fig.add_subplot(gs_res, sharex=ax)
    ax.plot(tf.TS, y, "o", color="k", ms=3, label="data")
    ax.plot(tf_, sh(R[3], R[2]), color="#0072B2", lw=1, label="1 pulse")
    ax.plot(tf_, sh(R[9], R[8]) + sh(R[11], R[10]), color="#D55E00", lw=1.2, label="2 pulses")
    ax.plot(tf_, sh(R[9], R[8]), color="#D55E00", lw=0.6, ls=":")
    ax.plot(tf_, sh(R[11], R[10]), color="#D55E00", lw=0.6, ls=":")
    ar.axhline(0, color="0.6", lw=0.5)
    ar.plot(tf.TS, y - m1, "o-", color="#0072B2", ms=2, lw=0.6)
    ar.plot(tf.TS, y - m2, "s-", color="#D55E00", ms=2, lw=0.6)
    lo = max(min(R[9], R[11]) - 40, 0)
    ax.set_xlim(lo, lo + 160)
    ax.set_title(f"{title}\n{LAB[G.g]} {G.sid[i] // 100000}/{G.sid[i] % 100000} ev {G.ev[i]} idx {G.idx[i]}"
                 f" dchi {G.dchi[i]:.0f}\ndt {G.dt[i]:+.1f} ns ratio {G.ratio[i]:.2f} A1 {G.pa[i]:.0f} MF npk {int((G.xch == i).sum()) + 1}", fontsize=6)
    ax.tick_params(labelsize=6)
    ar.tick_params(labelsize=6)
    import matplotlib.pyplot as plt
    plt.setp(ax.get_xticklabels(), visible=False)
    ar.set_xlabel("t (ns)", fontsize=6)
    ar.set_ylabel("res (ADC)", fontsize=6)
    ax.set_ylabel("ADC", fontsize=6)
    return ax


def pick(rng, mask, k):
    ii = np.nonzero(mask)[0]
    return list(rng.choice(ii, size=min(k, len(ii)), replace=False)) if len(ii) else []


def figures(GG, U, EXP, OP, rf_rows, scan_rows_, fref):
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    plt.rcParams.update({"font.size": 8, "axes.grid": True, "grid.alpha": 0.25})
    thr_g = np.array(sorted(set(THR_FINE)))
    with PdfPages(OUT / "u5_realdata.pdf") as pdf:
        # --- (1) detection rate vs threshold
        fig, axs = plt.subplots(2, 2, figsize=(11, 8))
        for ci, (grp, ttl) in enumerate((((("0p73", "2p24", "3p49", "c2p24")), "H2 and C"), (TA, "Ta (per run)"))):
            for ri, (key, lab) in enumerate((("close", "|dt| < 16 ns, any ratio"), ("win", "|dt| <= 80 ns (fit window), any ratio"))):
                ax = axs[ri, ci]
                for g in grp:
                    if g not in GG:
                        continue
                    G = GG[g]
                    ad = np.abs(G.dt)
                    msk = ad < C_CLOSE if key == "close" else ad <= REACH
                    y = np.array([(G.det(t) & msk).sum() / G.n for t in thr_g])
                    ax.plot(thr_g, np.where(y > 0, y, np.nan), "-o", ms=2.5, lw=1, color=COL[g], label=LAB[g])
                    e = EXP[g]["all"]
                    if key == "close":
                        ax.axhspan(min(e["a"]["close"], e["b"]["close"]), max(e["a"]["close"], e["b"]["close"]), color=COL[g], alpha=0.18, lw=0)
                    else:
                        ax.axhspan(min(e["a"]["win"], e["b"]["win"]), max(e["a"]["win"], e["b"]["win"]), color=COL[g], alpha=0.18, lw=0)
                    ax.plot([PROD_THR[g]], [np.interp(PROD_THR[g], thr_g, y)], marker="*", ms=9, color=COL[g], mec="k", ls="none")
                ax.set_xscale("log")
                ax.set_yscale("log")
                ax.set_xlabel("dchi2 threshold")
                ax.set_ylabel("detections per fired channel")
                ax.set_title(f"{ttl}: {lab}", fontsize=8)
                ax.legend(fontsize=6)
        fig.suptitle("Real-data two-pulse detections per fired PbWO4 channel (non-degenerate; star = production threshold).\n"
                     "Shaded: genuine expectation from the MF density (ratio > 0.3 only, bands 24-60 / 16-60 ns); detections above it are fakes or small-ratio pulses", fontsize=8)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        pdf.savefig(fig)
        plt.close(fig)
        # --- (2) dt distribution at the operating point
        fig, axs = plt.subplots(2, 4, figsize=(13, 6.5), sharex=True)
        edges = np.arange(-80, 80.1, 2.0)
        for ax, g in zip(axs.ravel(), ORDER):
            if g not in GG:
                ax.axis("off")
                continue
            G = GG[g]
            thr = thr_used(g, OP)
            d = G.det(thr)
            h = np.histogram(G.dt[d], bins=edges)[0] / G.n / 2.0
            hm = np.histogram(G.xdt, bins=edges)[0] / G.n / 2.0
            ax.stairs(np.where(h > 0, h, np.nan), edges, color=COL[g], lw=1.2, label=f"fit, thr {thr:g}")
            ax.stairs(np.where(hm > 0, hm, np.nan), edges, color="k", lw=0.8, ls="--", label="MF extras")
            ax.axvspan(-16, 16, color="0.5", alpha=0.12)
            ax.set_yscale("log")
            ax.set_title(f"{LAB[g]} ({G.n} ch)", fontsize=8)
            ax.legend(fontsize=6)
            ax.set_xlabel("dt = t_other - t_1 (ns)")
            ax.set_ylabel("per channel per ns")
        axs.ravel()[-1].axis("off")
        fig.suptitle("dt distribution of real two-pulse detections (2 ns bins) at the operating threshold vs production-MF extra pulses on the same channels", fontsize=8)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        pdf.savefig(fig)
        plt.close(fig)
        # --- (3) comb / grid phase
        fig, axs = plt.subplots(2, 2, figsize=(11, 8))
        ax = axs[0, 0]
        xs_ = np.arange(len(ORDER))
        for k, (lab, off, mk) in enumerate((("<16", -0.2, "o"), ("16-80", 0.0, "s"))):
            f_ = []
            e_ = []
            for g in ORDER:
                r = [x for x in rf_rows if x["group"] == g and x["subset"] == lab and abs(x["thr"] - thr_used(g, OP)) < 1e-6]
                r = r or [x for x in rf_rows if x["group"] == g and x["subset"] == lab]
                f_.append(r[0]["f_on"] if r else np.nan)
                e_.append(r[0]["f_on_err"] if r else np.nan)
            ax.errorbar(xs_ + off, f_, e_, fmt=mk, label=f"fit detections |dt| {lab}")
        mfv, mfe = [], []
        for g in ORDER:
            G = GG[g]
            m = (np.abs(G.xdt) > 24) & (np.abs(G.xdt) < REACH) & (G.xratio > 0.3)
            n = int(m.sum())
            k = int((np.abs(fold(G.xdt[m], pr.RF_P)) < 0.5).sum())
            mfv.append(k / max(n, 1))
            mfe.append(wilson_err(k, n))
        ax.errorbar(xs_ + 0.2, mfv, mfe, fmt="^", color="k", label="MF-complete pairs 24-80 ns")
        ax.axhline(0.25, color="r", lw=0.8, ls="--", label="random 25 %")
        ax.set_xticks(xs_)
        ax.set_xticklabels([LAB[g] for g in ORDER], rotation=30, fontsize=6)
        ax.set_ylabel("on-comb fraction (|fold| < 0.5 ns, P = 4.0077)")
        ax.legend(fontsize=6)
        ax.set_title("on-comb fraction", fontsize=8)
        ax = axs[0, 1]
        bins = np.linspace(-2, 2, 17)
        for g in ORDER:
            G = GG[g]
            thr = thr_used(g, OP)
            m = G.det(thr) & (np.abs(G.dt) >= 16) & (np.abs(G.dt) <= REACH)
            ax.stairs(np.histogram(fold(G.dt[m], pr.RF_P), bins=bins)[0] / max(m.sum(), 1) / (bins[1] - bins[0]), bins, color=COL[g], lw=1, label=LAB[g])
            ax.stairs(np.histogram(fold(G.dt[m], pr.GRID_P), bins=bins)[0] / max(m.sum(), 1) / (bins[1] - bins[0]), bins, color=COL[g], lw=0.8, ls=":")
            mc = G.det(thr) & (np.abs(G.dt) < 16)
            ax.stairs(np.histogram(fold(G.dt[mc], pr.RF_P), bins=bins)[0] / max(mc.sum(), 1) / (bins[1] - bins[0]), bins, color=COL[g], lw=0.7, ls="--")
        ax.axhline(1 / pr.RF_P, color="r", lw=0.8, ls="--")
        ax.set_xlabel("fold(dt, 4.0077) (solid) / fold(dt, 4.000) (dotted) (ns), 16 <= |dt| <= 80")
        ax.set_ylabel("per ns (normalised)")
        ax.legend(fontsize=6)
        ax.set_title("folded dt of detections (solid = RF period, dotted = grid, 16-80 ns; dashed = RF period, |dt| < 16)", fontsize=7)
        ax = axs[1, 0]
        pb = np.linspace(0, 4, 9)
        for g in ORDER:
            G = GG[g]
            thr = thr_used(g, OP)
            m = G.det(thr)
            h = np.histogram(np.mod(G.on_bank[m], 4.0), bins=pb)[0]
            ax.stairs(h / max(h.sum(), 1) * 8, pb, color=COL[g], lw=1, label=f"{LAB[g]} (n {int(m.sum())})")
            h1 = np.histogram(np.mod(G.p1_bank[m], 4.0), bins=pb)[0]
            ax.stairs(h1 / max(h1.sum(), 1) * 8, pb, color=COL[g], lw=0.8, ls=":")
        ax.axhline(1.0, color="r", lw=0.8, ls="--")
        ax.set_xlabel("secondary onset (bank axis) mod 4.000 ns")
        ax.set_ylabel("relative to flat")
        ax.set_title("secondary (solid) and pulse-1 (dotted) onset phase on the sample grid, detections at the operating threshold", fontsize=7)
        ax.legend(fontsize=6)
        ax = axs[1, 1]
        for k, lab in enumerate(("<16", "16-80")):
            f_, e_ = [], []
            for g in ORDER:
                r = [x for x in rf_rows if x["group"] == g and x["subset"] == lab and abs(x["thr"] - thr_used(g, OP)) < 1e-6]
                r = r or [x for x in rf_rows if x["group"] == g and x["subset"] == lab]
                f_.append(r[0]["f_flat"] if r else np.nan)
                e_.append(r[0]["f_flat_err"] if r else np.nan)
            ax.errorbar(xs_ + (k - 0.5) * 0.25, f_, e_, fmt="os"[k], label=f"|dt| {lab}")
        ax.axhline(0, color="0.5", lw=0.6)
        ax.axhline(1, color="r", lw=0.6, ls="--")
        ax.set_xticks(xs_)
        ax.set_xticklabels([LAB[g] for g in ORDER], rotation=30, fontsize=6)
        ax.set_ylabel(f"f_flat = 1 - (f_on - 0.25)/({fref:.3f} - 0.25)")
        ax.legend(fontsize=6)
        ax.set_title("flat-phase (RF-uncorrelated) fraction estimate", fontsize=8)
        fig.suptitle("RF comb vs sample grid for fitted close pairs (operating threshold; comb == grid at |dt| < 80 ns)", fontsize=8)
        fig.tight_layout(rect=(0, 0, 1, 0.96))
        pdf.savefig(fig)
        plt.close(fig)
        # --- (4) examples of detections, (5) examples of rejected / fake candidates
        rng = np.random.default_rng(20261028)
        ex_det, ex_fake = [], []
        for g in ("2p24", "0p73", "3p49", "x17_25800", "x17_26138", "x17_26230"):
            G = GG[g]
            thr = thr_used(g, OP)
            near_mf = np.zeros(G.n, bool)
            mm = np.abs(G.xt - G.ot[G.xch]) <= 20
            near_mf[G.xch[mm]] = True
            mask = G.det(thr) & (np.abs(G.dt) >= 4) & (np.abs(G.dt) < 16) & (G.ratio > 0.3) & (G.pa > 200) & (G.pa < 3000) & ~near_mf
            ex_det += [(G, i, f"detection |dt|<16 (thr {thr:g}, MF has no pulse there)") for i in pick(rng, mask, 1)]
        for g in ("2p24", "c2p24", "3p49"):
            G = GG[g]
            thr = round(PROD_THR[g], 2)
            mask = G.veto(thr) & (G.pa > 150)
            ex_fake += [(G, i, f"vetoed degenerate (|dt|<2 or ratio<0.05), thr {thr:g}") for i in pick(rng, mask, 1)]
        for g in ("2p24", "0p73", "3p49"):
            G = GG[g]
            thr = round(PROD_THR[g], 2)
            mask = G.det(thr) & (G.dt > 16) & (G.ratio < 0.12) & (G.pa > 150)
            ex_fake += [(G, i, f"e-p detection after main, ratio<0.12 (tail ripple?), thr {thr:g}") for i in pick(rng, mask, 1)]
        for name, ex in (("detections", ex_det), ("rejected / fake candidates", ex_fake)):
            fig = plt.figure(figsize=(12, 9))
            gs = fig.add_gridspec(4, 3, height_ratios=[3, 1, 3, 1], hspace=0.55, wspace=0.28)
            for k, (G, i, ttl) in enumerate(ex):
                r_, c_ = divmod(k, 3)
                ex_panel(fig, gs[2 * r_, c_], gs[2 * r_ + 1, c_], G, int(i), 0, ttl)
            fig.suptitle(f"Real-data examples: {name} (random picks, fixed seed). Black: data (ped-subtracted); blue: 1-pulse fit; orange: 2-pulse fit (dotted: components); lower: data - fit"
                         + ("\nnot classed genuine / fake by eye: at e-p rates (0.02-0.05 % close pairs) most such detections near the threshold are expected to be fakes" if name == "detections" else ""), fontsize=7)
            pdf.savefig(fig)
            plt.close(fig)


# ------------------------------------------------------------------- driver --
def target_for(group, args):
    return args.target_ta if group.startswith("x17") else args.target


def run_decode_cli(a):
    run, split = a.run, a.split
    dp = DUMP / f"{run:06d}_{split:05d}.npz"
    g = str(np.load(dp, allow_pickle=True)["group"])
    res = decode_split(run, split, str(dp), target_for(g, a), a.max_ev)
    f = save_raw(res)
    log(f"[decode] {run}/{split} {g}: {res['n_ev']} events, {res['n_chan']} channels "
        f"({res['n_not_in_bank']} pulse-1 channels outside the bank), alignment OK, {res['decode_s']:.0f} s -> {f}")


def stage_decode(a):
    jobs = [(r, s) for r, s, f in dump_files() if not a.smoke or (r, s) == (25228, 50)]
    todo = [(r, s) for r, s in jobs if a.force or not (SCRATCH / "raw" / f"{r:06d}_{s:05d}.npz").exists()]
    log(f"[decode] {len(jobs)} splits, {len(todo)} to do")
    (SCRATCH / "logs").mkdir(parents=True, exist_ok=True)
    running = []
    q = list(todo)
    while q or running:
        while q and len(running) < a.workers:
            r, s = q.pop(0)
            cmd = [sys.executable, os.path.abspath(__file__), "decode-one", "--run", str(r), "--split", str(s),
                   "--target", str(a.target), "--target-ta", str(a.target_ta)]
            if a.smoke:
                cmd += ["--max-ev", str(a.smoke_events)]
            lg = open(SCRATCH / "logs" / f"decode_{r:06d}_{s:05d}.log", "w")
            running.append((subprocess.Popen(cmd, stdout=lg, stderr=subprocess.STDOUT, cwd=str(HERE)), r, s, lg))
        time.sleep(1.0)
        for it in list(running):
            pr_, r, s, lg = it
            if pr_.poll() is None:
                continue
            running.remove(it)
            lg.close()
            txt = open(lg.name).read().splitlines()
            last = [x for x in txt if x.startswith("[decode]")]
            log(last[-1] if pr_.returncode == 0 and last else f"[FAIL rc={pr_.returncode}] {r}/{s}: " + " | ".join(txt[-3:]))
            if pr_.returncode:
                sys.exit(1)


def stage_fit(a):
    files = sorted((SCRATCH / "raw").glob("*.npz"))
    if a.smoke:
        files = [f for f in files if f.name.startswith("025228_00050")]
    (SCRATCH / "fit").mkdir(parents=True, exist_ok=True)
    todo = [f for f in files if a.force or not (SCRATCH / "fit" / f.name).exists()]
    log(f"[fit] {len(files)} raw caches, {len(todo)} to fit")
    for f in todo:                       # load every bank before the fork
        get_bank(str(np.load(f, allow_pickle=True)["slug"]))
    pool = mp.get_context("fork").Pool(a.workers)
    t0 = time.time()
    for f in todo:
        o = fit_split(f, pool)
        ff = SCRATCH / "fit" / f.name
        np.savez(str(ff) + ".tmp.npz", **o)
        os.replace(str(ff) + ".tmp.npz", ff)
        log(f"[fit] {f.name}: {len(o['R'])} channels, {o['fit_s']:.0f} s ({1e3 * o['fit_s'] * a.workers / max(len(o['R']), 1):.1f} ms/ch/worker); "
            f"total {time.time() - t0:.0f} s")
    pool.close()
    pool.join()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0], formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("stage", choices=["all", "decode", "fit", "report", "decode-one"])
    ap.add_argument("--workers", type=int, default=10)
    ap.add_argument("--target", type=int, default=30000, help="pulse-1 channels per H2 / C split")
    ap.add_argument("--target-ta", type=int, default=30000, help="pulse-1 channels per Ta split (3 per run)")
    ap.add_argument("--smoke", action="store_true")
    ap.add_argument("--smoke-events", type=int, default=300)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--run", type=int)
    ap.add_argument("--split", type=int)
    ap.add_argument("--max-ev", type=int, default=0)
    a = ap.parse_args(argv)
    if a.stage == "decode-one":
        return run_decode_cli(a)
    if a.stage in ("all", "decode"):
        stage_decode(a)
    if a.stage in ("all", "fit"):
        stage_fit(a)
    if a.stage in ("all", "report"):
        stage_report(a)


if __name__ == "__main__":
    main()
