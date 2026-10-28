#!/usr/bin/env python3
"""noise_diag.py -- FADC baseline-noise diagnostic (2026-10-28).

Question: is the HyCal PbWO4 FADC baseline noise white / uncorrelated (stored per-module sigma usable as the
diagonal chi2 sigma, ML denoising closed) or structured?  And why is the stored sig_ped (WaveAnalyzer pedestal
rms, median ~0.40 ADC) ~2.3x below the measured baseline sample std (~0.92 ADC,
output/ml_pileup_20261027/U0/overfit_diagnostics.md)?

Data.  The HyCal FADC readout is zero-suppressed: a physics event carries only ~16 (e-p) / ~43 (Ta) of the 1152
PbWO4 channels, every one of them with a pulse, so pulse-free channel-events do not exist.  Baseline sample =
the pre-pulse window (samples 0..29, the WaveAnalyzer pedestal window) of every present channel-event whose earliest
MF pulse (pileup_rates dump, fallback rows included) starts >= 132 ns (= sample 30 + 3 samples margin; MF onset =
t - bank t50); these are ~99 % of the channel-events (main pulses start at ~145 ns).  Late-only channels (earliest
onset >= 204 ns, ~1 %) give a 48-sample window.  The only events with ALL channels present are the LMS flasher events
((trig & 0x3000000) == 0x1000000, ~65 per split file; their pulse starts after sample ~36, window 0..29 clean); they
give complete boards for the same-sample correlation and the leave-one-out (LOO) common mode.  Physics selection =
pileup_rates.process_split ((trig & 0xf00) != 0, (trig & 0x3000000) == 0), event numbers asserted against the dump.

Stages (scratch caches under SCRATCH, resumable; one subprocess per split for decode):
  decode   split -> sparse physics rows (ev, col, raw[100]) + LMS events (all channels)
  analyse  per split: pedestal definitions, spectra, acf, fixed pattern, correlation, drift, common mode (physics rows;
           LMS events pooled over all decoded files)
  report   caches -> TSVs, noise_diag.pdf, numbers.json (report.md / notes are written by hand from these)

    MPLBACKEND=Agg python3 noise_diag.py all [--workers 6]
"""
from __future__ import annotations

import argparse
import json
import math
import multiprocessing as mp
import os
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
OUT = HERE / "output/noise_diag_20261028"
DUMP = HERE / "output/pileup_rates_20261027/dump"
HMAP = HERE.parents[1] / "database/hycal_map.json"
SCRATCH = Path("/private/tmp/claude-505/-Users-L00338853-work-PRad-prad2evviewer-analysis-pyscripts/"
               "b1de5a72-6ab3-4e2b-8b9b-77e955030d7b/scratchpad/noise_diag")
NS, NW, NL, NCH = 100, 30, 48, 1152
SPLITS = [("0p73", 24590, 10), ("2p24", 25228, 50), ("3p49", 24503, 50),
          ("x17_26230", 26230, 100), ("x17_25800", 25800, 100)]
LMS_EXTRA = [(24386, 20), (25320, 300), (24590, 80), (26230, 600), (25800, 400)]   # LMS events only
LAB = {"0p73": "0.73 GeV H2", "2p24": "2.24 GeV H2", "3p49": "3.49 GeV H2",
       "x17_26230": "Ta 26230", "x17_25800": "Ta 25800"}
FS = 250.0                                 # MHz sample rate
SM_ORDER = 2                               # database/daq_config.json wave smooth_order
ONSET30, ONSET48 = 12.0 + 4 * NW, 12.0 + 4 * NL   # ns: window end + 3 samples margin
CH_BLK = 300
NF30, NF48 = NW // 2 + 1, NL // 2 + 1
KACF = 12
MIN_FREE_B = 4                             # present clean channels in a board for a LOO estimate
LMS_TB = 0x1000000
SGN = (-1.0) ** np.arange(NS)               # +/- alternating (Nyquist) pattern, sample index parity
f32, f64 = np.float32, np.float64


def log(*a):
    print(*a, flush=True)


def key(run, split):
    return f"{run:06d}_{split:05d}"


# ------------------------------------------------------------------ WaveAnalyzer replica --
def smooth32(raw, order=SM_ORDER):
    """WaveAnalyzer::smooth in float32 (fma emulated); edges drop the missing taps (as C++)."""
    x = raw.astype(f32)
    n = x.shape[-1]
    val, wsum = x.copy(), np.ones(n, f32)
    for j in range(1, order):
        w = f32(1) - f32(j) / f32(order + 1)
        s = np.zeros_like(x)
        s[..., j:n - j] = x[..., :n - 2 * j] + x[..., 2 * j:]
        val = (f64(w) * s.astype(f64) + val.astype(f64)).astype(f32)
        ws = np.zeros(n, f32)
        ws[j:n - j] = 1
        wsum = (2.0 * f64(w) * ws + wsum.astype(f64)).astype(f32)
    return (val / wsum).astype(f32)


def wa_ped(buf, flat=1.0, max_iter=3):
    """WaveAnalyzer::findPedestal (median/MAD seed, sigma-clip band max(rms, flat), <= 3 iterations,
    float32 running sums as the C++).  buf [n, N] smoothed window; returns (mean, rms)."""
    n, N = buf.shape
    h = N // 2
    srt = np.sort(buf, 1)
    mean = ((srt[:, h - 1] + srt[:, h]) * f32(.5)).astype(f32)
    dev = np.sort(np.abs(buf - mean[:, None]), 1)
    rms = (((dev[:, h - 1] + dev[:, h]) * f32(.5)) * f32(1.4826)).astype(f32)
    active, done = np.ones(buf.shape, bool), np.zeros(n, bool)
    for _ in range(max_iter):
        band = np.maximum(rms, f32(flat))
        keep = active & (np.abs(buf - mean[:, None]) < band[:, None])
        cnt, acnt = keep.sum(1), active.sum(1)
        conv = cnt == acnt
        few = (~done) & (~conv) & (cnt < 5)
        upd = (~done) & (~conv) & (~few)
        active = np.where(done[:, None], active, keep)
        sm, sm2 = np.zeros(n, f32), np.zeros(n, f32)
        for i in range(N):
            b = np.where(keep[:, i], buf[:, i], f32(0))
            sm = (sm.astype(f64) + b).astype(f32)
            sm2 = (sm2.astype(f64) + b.astype(f64) * b.astype(f64)).astype(f32)
        k = np.maximum(cnt, 1).astype(f32)
        m2 = (sm / k).astype(f32)
        var = (sm2 / k - m2 * m2).astype(f32)
        r2 = np.where(var > 0, np.sqrt(np.maximum(var, f32(0))), f32(0)).astype(f32)
        mean, rms = np.where(upd, m2, mean), np.where(upd, r2, rms)
        done |= conv | few
    return mean, rms


def mad_rms(x):
    """1.4826 * MAD of x[n, N] (median of |x - median|), the WaveAnalyzer seed."""
    med = np.median(x, 1, keepdims=True)
    return (f32(1.4826) * np.median(np.abs(x - med), 1)).astype(f32)


# ------------------------------------------------------------------ geometry --
def hmap_info(names):
    m = {e["n"]: e["daq"] for e in json.load(open(HMAP))}
    cr = np.array([m[n]["crate"] for n in names])
    sl = np.array([m[n]["slot"] for n in names])
    chn = np.array([m[n]["channel"] for n in names])
    return cr, sl, chn, cr * 100 + sl


# ------------------------------------------------------------------ decode --
def decode_one(run, split, lms_only):
    import _common as C
    import pileup_rates as pr
    from prad2py import dec
    z = np.load(DUMP / f"{key(run, split)}.npz", allow_pickle=True)
    pwo4 = np.asarray(z["module_is_pwo4"], bool)
    names_all = np.asarray(z["module_names"])
    cols = np.flatnonzero(pwo4)
    colof = np.full(len(pwo4), -1, np.int64)
    colof[cols] = np.arange(len(cols))
    assert len(cols) == NCH
    n_sel = int(z["n_sel"])
    path = os.path.join(pr.EVIO_DIR, f"prad_{run:06d}.evio.{split:05d}")
    p = C.setup_pipeline(evio_path=path, hc_calib_file=str(z["calib_file"]))
    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    if ch.open_auto(path) != dec.Status.success:
        sys.exit(f"[ERROR] cannot open {path}")
    cache, ev_num = {}, []
    P_ev, P_col, P_raw = [], [], []
    L_raw, L_got, L_num = [], [], []
    n_evt, t0 = 0, time.time()
    while ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            e = d["event"]
            tb = int(e.info.trigger_bits)
            lms = (tb & 0x3000000) == LMS_TB
            phys = bool(tb & pr.TRIG_MASK) and not (tb & pr.TRIG_EXCL) and n_evt < n_sel and not lms_only
            if not (phys or lms):
                continue
            if lms:
                lr, lg = np.zeros((NCH, NS), np.uint16), np.zeros(NCH, bool)
            if phys:
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
                        k = (crate, s, c)
                        col = cache.get(k)
                        if col is None:
                            mod = p.hycal.module_by_daq(crate, s, c)
                            col = int(colof[mod.index]) if (mod is not None and mod.is_hycal()) else -1
                            cache[k] = col
                        if col < 0:
                            continue
                        cd = slot.channel(c)
                        assert cd.nsamples == NS
                        if phys:
                            P_ev.append(n_evt); P_col.append(col); P_raw.append(cd.samples)
                        if lms:
                            lr[col] = cd.samples; lg[col] = True
            if lms:
                L_raw.append(lr); L_got.append(lg); L_num.append(int(e.info.event_number))
            if phys:
                n_evt += 1
    d_ = SCRATCH / "raw"
    d_.mkdir(parents=True, exist_ok=True)
    if not lms_only:
        assert n_evt == n_sel and np.array_equal(np.asarray(ev_num), np.asarray(z["ev_num"])), f"alignment failed {run}/{split}"
        names = names_all[cols]
        cr, sl, chn, board = hmap_info(names)
        np.savez(d_ / f"{key(run, split)}_meta.npz", ev_num=np.asarray(ev_num, np.int64), cols=cols, names=names,
                 crate=cr, slot=sl, chan=chn, board=board, group=str(z["group"]), templates=str(z["templates"]),
                 run=run, split=split, n_sel=n_sel)
        np.savez(str(d_ / f"{key(run, split)}_phys.tmp.npz"), ev=np.asarray(P_ev, np.int32), col=np.asarray(P_col, np.int16),
                 raw=np.asarray(P_raw, np.uint16))
        os.replace(d_ / f"{key(run, split)}_phys.tmp.npz", d_ / f"{key(run, split)}_phys.npz")
    np.savez(str(d_ / f"{key(run, split)}_lms.tmp.npz"), raw=np.asarray(L_raw, np.uint16).reshape(-1, NCH, NS),
             got=np.asarray(L_got, bool).reshape(-1, NCH), ev_num=np.asarray(L_num, np.int64))
    os.replace(d_ / f"{key(run, split)}_lms.tmp.npz", d_ / f"{key(run, split)}_lms.npz")
    log(f"[decode] {run}/{split}: {len(ev_num)} selected events, {len(P_ev)} rows, {len(L_num)} LMS events, "
        f"{time.time() - t0:.0f} s")


def stage_decode(a):
    jobs = [(r, s, False) for _, r, s in SPLITS] + [(r, s, True) for r, s in LMS_EXTRA]
    todo = [j for j in jobs if a.force or not (SCRATCH / "raw" / f"{key(j[0], j[1])}_lms.npz").exists()]
    log(f"[decode] {len(jobs)} splits, {len(todo)} to do")
    procs = []
    (SCRATCH / "logs").mkdir(parents=True, exist_ok=True)
    for r, s, lo in todo:
        lg = open(SCRATCH / "logs" / f"decode_{key(r, s)}.log", "w")
        cmd = [sys.executable, os.path.abspath(__file__), "decode-one", "--run", str(r), "--split", str(s)]
        if lo:
            cmd.append("--lms-only")
        procs.append((subprocess.Popen(cmd, stdout=lg, stderr=subprocess.STDOUT, cwd=str(HERE)), r, s, lg))
    for pr_, r, s, lg in procs:
        rc = pr_.wait()
        lg.close()
        txt = open(lg.name).read().splitlines()
        log(([x for x in txt if x.startswith("[decode]")] or txt[-3:])[-1] if rc == 0 else f"[FAIL] {r}/{s}: " + " | ".join(txt[-3:]))
        if rc:
            sys.exit(1)


# ------------------------------------------------------------------ analysis accumulators --
def acf_sums(x, K):
    """x [m, N] mean-removed; sum over windows of r_k / (N-k), k < K."""
    N = x.shape[1]
    F = np.fft.rfft(x.astype(f64), n=2 * N, axis=1)
    r = np.fft.irfft(F.real ** 2 + F.imag ** 2, axis=1)[:, :K]
    return (r / (N - np.arange(K))).sum(0)


def bincount_cols(ch, M, n=NCH):
    return np.stack([np.bincount(ch, weights=M[:, j], minlength=n) for j in range(M.shape[1])], 1)


HIST_BINS = np.linspace(0, 3.0, 121)
Q = ("rstd0", "rmad", "sstd0", "smad", "war")
CMK = ("n", "ssq0", "ssq1", "wm", "mm", "inv")


class Acc:
    """Accumulates every baseline statistic over dense chunks X[c, NCH, NS] with window masks."""

    def __init__(self, board, res=False):
        self.res = res                      # True: subtract the per-window alternating (Nyquist) amplitude first
        self.board = board
        self.ub, self.bidx = np.unique(board, return_inverse=True)
        self.members = [np.flatnonzero(self.bidx == b) for b in range(len(self.ub))]
        self.qs = {q: [] for q in Q}
        self.hist = {q: np.zeros(len(HIST_BINS) - 1) for q in Q}
        self.ch, self.ev, self.mu, self.ssq = [], [], [], []
        self.S30, self.S48 = np.zeros((NCH, NF30)), np.zeros((NCH, NF48))
        self.n30, self.n48 = np.zeros(NCH), np.zeros(NCH)
        self.acf30, self.acf48 = np.zeros(NW), np.zeros(NL)
        self.pat, self.oe = np.zeros((NCH, NW)), np.zeros((NCH, 2))
        self.CM = {k: np.zeros(len(self.ub)) for k in CMK}
        self.SA, self.SB, self.NP = np.zeros((NCH, NCH)), np.zeros((NCH, NCH)), np.zeros((NCH, NCH))
        self.han = np.hanning(NL + 1)[:-1]
        self.ev0 = 0
        self.alt, self.altev, self.altch = [], [], []
        self.Pcm = np.zeros((3, NF30))       # periodograms of: LOO board mean, LOO residual, channel itself (boards with >= 5 clean channels)
        self.Pcm_n = 0.0
        self.mom = np.zeros(4)               # n, sum a, sum a^2, sum a^4 of the 30-sample alternating amplitude
        self.coh = np.zeros(4)               # event-pair numerator / denominator (all channels), same board
        self.per = np.zeros(6)               # persistence: n, ssq before, ssq after, sum ae*al, sum ae^2, sum al^2
        onehot = np.zeros((NCH, len(self.ub)), np.float32)
        onehot[np.arange(NCH), self.bidx] = 1
        self.onehot = onehot

    def add(self, X, free30, free48):
        c = X.shape[0]
        a_all = (X[:, :, :NW] * SGN[:NW]).mean(2)             # alternating amplitude (sum of SGN over 30 samples is 0)
        f30 = free30.astype(f32)
        if not self.res:
            a_f = a_all * f30
            n_e = f30.sum(1)
            S_e, Q_e = a_f.sum(1), (a_f ** 2).sum(1)
            self.coh[0] += float((S_e ** 2 - Q_e).sum()); self.coh[1] += float(((n_e - 1).clip(0) * Q_e).sum())
            Sb, Qb, nb = a_f @ self.onehot, (a_f ** 2) @ self.onehot, f30 @ self.onehot
            self.coh[2] += float((Sb ** 2 - Qb).sum()); self.coh[3] += float(((nb - 1).clip(0) * Qb).sum())
            av = a_all[free30].astype(f64)
            self.mom += [av.size, av.sum(), (av ** 2).sum(), (av ** 4).sum()]
            ce0, ch0 = np.nonzero(free30)
            self.alt.append(a_all[ce0, ch0].astype(f32)); self.altev.append((ce0 + self.ev0).astype(np.int32)); self.altch.append(ch0.astype(np.int16))
            cp_, hp_ = np.nonzero(free48)
            if len(cp_):
                dd = X[cp_, hp_][:, :NL].astype(f64)
                dd = dd - dd[:, :NW].mean(1, keepdims=True)
                ae = a_all[cp_, hp_].astype(f64)
                al = (dd[:, NW:NL] * SGN[NW:NL]).mean(1)
                rs = dd[:, NW:NL] - ae[:, None] * SGN[NW:NL]
                self.per += [len(cp_), (dd[:, NW:NL] ** 2).sum(), (rs ** 2).sum(), (ae * al).sum(), (ae ** 2).sum(), (al ** 2).sum()]
        else:
            X[:, :, :NW] -= (a_all * f30)[:, :, None] * SGN[:NW]
            X[:, :, NW:NL] -= (a_all * free48.astype(f32))[:, :, None] * SGN[NW:NL]
        ce, ch = np.nonzero(free30)
        Xw = X[ce, ch]
        w30 = Xw[:, :NW]
        fm = free30[:, :, None].astype(f32)
        gm = (X[:, :, :NW] * fm).sum((0, 2)) / np.maximum(fm.sum((0, 2)) * NW, 1)   # per-channel chunk mean (offset-kept centre)
        mu = w30.mean(1)
        d30 = w30 - mu[:, None]
        ssq = (d30 * d30).sum(1)
        s30 = smooth32(Xw)[:, :NW]
        _, war = wa_ped(s30)
        q = dict(rstd0=np.sqrt(ssq / NW), rmad=mad_rms(w30), sstd0=s30.std(1), smad=mad_rms(s30), war=war)
        for k in Q:
            self.qs[k].append(q[k].astype(f32))
            self.hist[k] += np.histogram(q[k], HIST_BINS)[0]
        self.ch.append(ch.astype(np.int16)); self.ev.append((ce + self.ev0).astype(np.int32))
        self.mu.append(mu); self.ssq.append(ssq)
        self.n30 += np.bincount(ch, minlength=NCH)
        F = np.fft.rfft(d30.astype(f64), axis=1)
        self.S30 += bincount_cols(ch, (F.real ** 2 + F.imag ** 2) / NW)
        self.acf30 += acf_sums(d30, NW)                         # per-window centred; debiased exactly in the report stage
        pw = d30[:, 0::2].mean(1) - d30[:, 1::2].mean(1)
        self.oe[:, 0] += np.bincount(ch, weights=pw, minlength=NCH)
        self.oe[:, 1] += np.bincount(ch, weights=pw * pw, minlength=NCH)
        self.pat += bincount_cols(ch, d30)
        # long window (late-only channels)
        ce2, ch2 = np.nonzero(free48)
        if len(ce2):
            d48 = X[ce2, ch2][:, :NL]
            d48 = d48 - d48.mean(1, keepdims=True)
            F = np.fft.rfft(d48.astype(f64) * self.han, axis=1)
            self.S48 += bincount_cols(ch2, (F.real ** 2 + F.imag ** 2) / (self.han ** 2).sum(), NCH)
            self.n48 += np.bincount(ch2, minlength=NCH)
            self.acf48 += acf_sums(d48, NL)
        # correlation, 30-sample windows
        Xc = X[:, :, :NW]
        A_ = ((Xc - Xc.mean(2, keepdims=True)) * fm).transpose(0, 2, 1).reshape(-1, NCH)
        B_ = ((Xc - gm[None, :, None]) * fm).transpose(0, 2, 1).reshape(-1, NCH)
        self.SA += (A_.T @ A_).astype(f64)
        self.SB += (B_.T @ B_).astype(f64)
        ff = free30.astype(f32)
        self.NP += (ff.T @ ff).astype(f64) * NW
        # leave-one-out board common mode, window 0..29, per-window mean removed
        for b, idx in enumerate(self.members):
            fb = free30[:, idx]
            nf = fb.sum(1)
            okev = nf >= MIN_FREE_B
            if not okev.any():
                continue
            Xb = X[:, idx, :NW]
            w = (Xb - Xb.mean(2, keepdims=True)) * fb[:, :, None]
            m = (w.sum(1)[:, None, :] - w) / np.maximum(nf - 1, 1)[:, None, None]
            val = (fb & okev[:, None])[:, :, None]
            w, m = w * val, m * val
            r = w - m
            big = (val[:, :, 0] & (nf >= 5)[:, None])
            if big.any():
                for q_, arr in enumerate((m, r, w)):
                    Fq = np.fft.rfft(arr[big].astype(f64), axis=1)
                    self.Pcm[q_] += ((Fq.real ** 2 + Fq.imag ** 2) / NW).sum(0)
                self.Pcm_n += int(big.sum())
            self.CM["n"][b] += val[:, :, 0].sum()
            self.CM["inv"][b] += (val[:, :, 0] / np.maximum(nf - 1, 1)[:, None]).sum()
            self.CM["ssq0"][b] += (w ** 2).sum()
            self.CM["ssq1"][b] += (r ** 2).sum()
            self.CM["wm"][b] += (w * m).sum()
            self.CM["mm"][b] += (m ** 2).sum()
        self.ev0 += c

    def result(self):
        chs, evs = np.concatenate(self.ch), np.concatenate(self.ev)
        qs = {k: np.concatenate(v) for k, v in self.qs.items()}
        mus, ssqs = np.concatenate(self.mu), np.concatenate(self.ssq)
        o = np.lexsort((evs, chs))
        chs, evs, mus, ssqs = chs[o], evs[o], mus[o], ssqs[o]
        qs = {k: v[o] for k, v in qs.items()}
        st = np.searchsorted(chs, np.arange(NCH + 1))
        med = {k: np.full(NCH, np.nan) for k in Q}
        drift = np.full((NCH, 2), np.nan)
        pooled = np.full(NCH, np.nan)
        for c in range(NCH):
            a, b = st[c], st[c + 1]
            if b - a < 20:
                continue
            for k in Q:
                med[k][c] = np.median(qs[k][a:b])
            pooled[c] = math.sqrt(ssqs[a:b].mean() / (NW - 1))
            m_ = mus[a:b].astype(f64)
            drift[c, 0] = m_.std(ddof=1)
            drift[c, 1] = math.sqrt(np.mean(np.diff(m_) ** 2) / 2)
        Np = np.maximum(self.NP, 1)
        corr = {}
        for tag, S in (("A", self.SA), ("B", self.SB)):
            cov = S / Np
            dg = np.sqrt(np.diag(cov))
            corr[tag] = (cov / np.outer(dg, dg)).astype(f32)
        n30 = np.maximum(self.n30, 1)
        pat = self.pat / n30[:, None]
        iu = np.triu_indices(NCH, 1)
        res = dict(nwin=self.n30, nwin48=self.n48, pooled=pooled, drift=drift, S30=self.S30, S48=self.S48,
                   acf30=self.acf30 / self.acf30[0], acf48=self.acf48 / self.acf48[0] if self.acf48[0] > 0 else self.acf48,
                   pattern=pat.astype(f32), oe_mean=self.oe[:, 0] / n30, oe_rms=np.sqrt(self.oe[:, 1] / n30),
                   corrA=corr["A"], corrB=corr["B"], npair=float(np.median(self.NP[iu])), np_mat=self.NP.astype(f32),
                   mom=self.mom, coh=self.coh, per=self.per, Pcm=self.Pcm / max(self.Pcm_n, 1), Pcm_n=self.Pcm_n,
                   npair_zero=float((self.NP[iu] == 0).mean()), ub=self.ub,
                   **{"med_" + k: v for k, v in med.items()}, **{"hist_" + k: v for k, v in self.hist.items()},
                   **{"cm_" + k: v for k, v in self.CM.items()})
        return res


def pulse_onsets(run, split, meta, nev):
    """Earliest MF-pulse onset (ns) per (event, col) from the pileup dump; fallback rows included."""
    import twopulse_realdata as trd
    z = np.load(DUMP / f"{key(run, split)}.npz", allow_pickle=True)
    names = [str(n) for n in z["module_names"]]
    bank = trd.get_bank(trd.slug_of(str(z["templates"])))
    ix = {n: i for i, n in enumerate(bank.names)}
    t50 = np.array([bank.t50[ix[n]] if n in ix else np.nan for n in names])
    t50 = np.where(np.isfinite(t50), t50, np.nanmedian(t50))
    colof = np.full(len(names), -1, np.int64)
    colof[meta["cols"]] = np.arange(NCH)
    ev, idx = z["pu_ev"].astype(np.int64), z["pu_idx"].astype(np.int64)
    on = z["pu_t"].astype(np.float64) - t50[idx]
    k = (ev < nev) & (colof[idx] >= 0)
    o = np.full((nev, NCH), np.inf)
    np.minimum.at(o, (ev[k], colof[idx[k]]), on[k])
    return o


RF_KEYS = ("Pcm", "pooled", "acf30", "acf48", "S30", "S48", "nwin", "nwin48", "corrA", "corrB", "cm_n", "cm_ssq0", "cm_ssq1", "cm_wm",
           "cm_mm", "cm_inv", "npair", "med_war", "med_sstd0", "med_rstd0", "drift", "ub")


def rf_phase_test(run, split, meta, acc):
    """Modulation of the 30-sample alternating amplitude with the RF phase (first channel-A hit of the event):
    returns {period: (|mean a exp(-i 2 pi t/T)| * 2, null expectation)} for T = RF/2 period, 2x, bunch period."""
    rf = np.load(HERE / "output/module_timing_20261025/rf" / f"{key(run, split)}.npz")
    r_ev, off = rf["ev_num"], rf["rf_off"].astype(np.int64)
    rt, rch = rf["rf_t"].astype(np.float64), rf["rf_ch"]
    ha = np.flatnonzero(rch == 0)
    eva = np.repeat(np.arange(len(r_ev)), np.diff(off))[ha]
    ta = rt[ha]
    uq, first = np.unique(eva, return_index=True)
    t_first = np.full(len(r_ev), np.nan)
    t_first[uq] = ta[first]
    en = meta["ev_num"]
    pos = np.minimum(np.searchsorted(r_ev, en), len(r_ev) - 1)
    t_ev = np.where(r_ev[pos] == en, t_first[pos], np.nan)
    a, ev = np.concatenate(acc.alt).astype(np.float64), np.concatenate(acc.altev)
    x = t_ev[ev]
    g = np.isfinite(x)
    out = {}
    for name, T in (("RF/2 (8.0154 ns)", 8.0154), ("RF/4 (16.0308 ns)", 16.0308), ("bunch (4.0077 ns)", 4.0077), ("clk/2 (8 ns)", 8.0)):
        z = (a[g] * np.exp(-2j * np.pi * x[g] / T)).mean()
        out[name] = (float(abs(z) * 2), float(a[g].std() / math.sqrt(g.sum()) * 2))
    return out


def analyse_split(run, split):
    t0 = time.time()
    d_ = SCRATCH / "raw"
    meta = dict(np.load(d_ / f"{key(run, split)}_meta.npz", allow_pickle=True))
    ph = np.load(d_ / f"{key(run, split)}_phys.npz")
    ev, col, raw = ph["ev"], ph["col"].astype(np.int64), ph["raw"]
    nev = int(meta["n_sel"])
    onset = pulse_onsets(run, split, meta, nev)
    acc, accr = Acc(meta["board"]), Acc(meta["board"], res=True)
    n_rows = n30 = n48 = 0
    for e0 in range(0, nev, CH_BLK):
        e1 = min(e0 + CH_BLK, nev)
        r0, r1 = np.searchsorted(ev, [e0, e1])
        X = np.zeros((e1 - e0, NCH, NS), f32)
        pres = np.zeros((e1 - e0, NCH), bool)
        X[ev[r0:r1] - e0, col[r0:r1]] = raw[r0:r1]
        pres[ev[r0:r1] - e0, col[r0:r1]] = True
        on = onset[e0:e1]
        f30, f48 = pres & (on >= ONSET30), pres & (on >= ONSET48)
        n_rows += int(pres.sum()); n30 += int(f30.sum()); n48 += int(f48.sum())
        acc.add(X, f30, f48)
        accr.add(X, f30, f48)
    res = acc.result()
    resr = accr.result()
    res.update({"r_" + k: resr[k] for k in RF_KEYS})
    res["rf_test"] = json.dumps(rf_phase_test(run, split, meta, acc))
    res.update(run=run, split=split, group=str(meta["group"]), names=meta["names"], board=meta["board"], chan=meta["chan"],
               crate=meta["crate"], n_events=nev, n_rows=n_rows, n_valid30=n30, n_valid48=n48, elapsed=time.time() - t0)
    f = SCRATCH / "analysis" / f"{key(run, split)}.npz"
    f.parent.mkdir(parents=True, exist_ok=True)
    np.savez(str(f) + ".tmp.npz", **res)
    os.replace(str(f) + ".tmp.npz", f)
    return f"[analyse] {run}/{split}: {nev} ev, {n_rows} rows, {n30} clean-30 ({100 * n30 / n_rows:.1f} %), {n48} clean-48, {time.time() - t0:.0f} s"


def analyse_lms():
    """LMS flasher events of every decoded file pooled (all channels present)."""
    t0 = time.time()
    d_ = SCRATCH / "raw"
    meta = dict(np.load(d_ / f"{key(*SPLITS[1][1:])}_meta.npz", allow_pickle=True))
    acc, accr = Acc(meta["board"]), Acc(meta["board"], res=True)
    n = 0
    mean_wave = np.zeros(NS)
    files = [(r, s) for _, r, s in SPLITS] + LMS_EXTRA
    for r, s in files:
        z = np.load(d_ / f"{key(r, s)}_lms.npz")
        raw, got = z["raw"], z["got"]
        for e0 in range(0, len(raw), 100):
            X = raw[e0:e0 + 100].astype(f32)
            g = got[e0:e0 + 100]
            mean_wave += ((X - X[:, :, :10].mean(2, keepdims=True)) * g[:, :, None]).sum((0, 1))
            acc.add(X, g, np.zeros_like(g))
            accr.add(X, g, np.zeros_like(g))
            n += len(X)
    res = acc.result()
    resr = accr.result()
    res.update({"r_" + k: resr[k] for k in RF_KEYS})
    res.update(n_lms=n, mean_wave=mean_wave / max(1, res["nwin"].sum()), names=meta["names"], board=meta["board"],
               chan=meta["chan"], crate=meta["crate"], elapsed=time.time() - t0)
    f = SCRATCH / "analysis" / "lms.npz"
    f.parent.mkdir(parents=True, exist_ok=True)
    np.savez(str(f) + ".tmp.npz", **res)
    os.replace(str(f) + ".tmp.npz", f)
    return f"[analyse] LMS: {n} events pooled, {time.time() - t0:.0f} s"


def _job(a):
    return analyse_lms() if a == "lms" else analyse_split(*a)


def stage_analyse(a):
    todo = [(r, s) for _, r, s in SPLITS if a.force or not (SCRATCH / "analysis" / f"{key(r, s)}.npz").exists()]
    if a.force or not (SCRATCH / "analysis" / "lms.npz").exists():
        todo = ["lms"] + todo
    log(f"[analyse] {len(todo)} to do")
    if todo:
        with mp.get_context("fork").Pool(min(a.workers, len(todo))) as pool:
            for m in pool.imap_unordered(_job, todo):
                log(m)


# ------------------------------------------------------------------ linear-fix test on real pulses --
FIX_SPLITS = ("2p24", "0p73", "3p49")


def _rf_phase(rf, ev_num_rows, t, P=None):
    """Nearest channel-A RF hit in the same event (as plot_module_timing_rf.prep); returns (phi, P, ok)."""
    import plot_module_timing_rf as RF
    r_ev, off = rf["ev_num"], rf["rf_off"].astype(np.int64)
    rf_t, rf_ch = rf["rf_t"].astype(np.float64), rf["rf_ch"]
    n_rf = len(r_ev)
    pos = np.minimum(np.searchsorted(r_ev, ev_num_rows), n_rf - 1)
    matched = r_ev[pos] == ev_num_rows
    ha = np.flatnonzero(rf_ch == 0)
    ev_a = np.repeat(np.arange(n_rf), np.diff(off))[ha]
    ta = rf_t[ha]
    o = np.lexsort((ta, ev_a))
    ev_a, ta = ev_a[o], ta[o]
    if P is None:
        dd = np.diff(ta)[ev_a[1:] == ev_a[:-1]]
        dd = dd[(dd > 120.0) & (dd < 136.0)]
        P = float(dd.mean() / RF.N_SUB) if dd.size else RF.P_NOM
    K = ev_a * RF.KEY_BIG + ta
    i = np.searchsorted(K, pos * RF.KEY_BIG + t)
    il, ih = np.clip(i - 1, 0, len(K) - 1), np.clip(i, 0, len(K) - 1)
    dl = np.where(ev_a[il] == pos, np.abs(ta[il] - t), np.inf)
    dh = np.where(ev_a[ih] == pos, np.abs(ta[ih] - t), np.inf)
    t_rf = np.where(dh < dl, ta[ih], ta[il])
    return RF.fold(t - t_rf, P), P, matched & np.isfinite(np.minimum(dl, dh))


def fix_split(run, split, workers):
    """Real-data test of the linear fix that IS possible with zero-suppressed readout: subtract the alternating
    (Nyquist) component a (-1)^s, a fitted on the channel's own pre-pulse window ('own') or taken as the mean of the
    other clean channels of the event ('common', >= 3 others), from the whole waveform; then MF time vs RF (E >= 100 MeV,
    plot_module_timing_rf estimator) and single-pulse two-pulse-fit dchi2 (OPT_PROD, stored sigma)."""
    import _common as C
    import mf_pulse
    import plot_module_timing_rf as RF
    import twopulse_fit as tf
    import twopulse_realdata as trd
    t0 = time.time()
    d_ = SCRATCH / "raw"
    meta = dict(np.load(d_ / f"{key(run, split)}_meta.npz", allow_pickle=True))
    ph = np.load(d_ / f"{key(run, split)}_phys.npz")
    pev, pcol, praw = ph["ev"].astype(np.int64), ph["col"].astype(np.int64), ph["raw"]
    nev = int(meta["n_sel"])
    onset = pulse_onsets(run, split, meta, nev)
    valid = onset[pev, pcol] >= ONSET30
    a_row = (praw[:, :NW].astype(f32) * SGN[:NW].astype(f32)).mean(1)
    S_ev = np.bincount(pev[valid], weights=a_row[valid], minlength=nev)
    n_ev = np.bincount(pev[valid], minlength=nev)
    pkey = pev * 2048 + pcol
    po = np.argsort(pkey)
    z = np.load(DUMP / f"{key(run, split)}.npz", allow_pickle=True)
    names_all = [str(n) for n in z["module_names"]]
    D = {k: z[k] for k in ("pu_ev", "pu_idx", "pu_t", "pu_A", "pu_E", "pu_fb")}
    cp = trd.channel_pulses(D)
    colof = np.full(len(names_all), -1, np.int64)
    colof[meta["cols"]] = np.arange(NCH)
    tmpl, slug = str(z["templates"]), trd.slug_of(str(z["templates"]))
    bank_t = trd.get_bank(slug)
    ix = {n: i for i, n in enumerate(bank_t.names)}
    fbset = set(bank_t.fallback)
    inbank = np.array([(n in ix) and (n not in fbset) for n in names_all])
    sel_t = cp["has1"] & (cp["E1"] >= 80.0)
    sel_s = cp["has1"] & (cp["npk"] == 1) & (cp["A1"] >= 100) & (cp["A1"] <= 4000) & inbank[cp["idx"]]
    rows = np.flatnonzero(sel_t | sel_s)
    r_ev, r_mod = cp["ev"][rows].astype(np.int64), cp["idx"][rows].astype(np.int64)
    r_col = colof[r_mod]
    k_ = r_ev * 2048 + r_col
    pos = np.minimum(np.searchsorted(pkey[po], k_), len(po) - 1)
    ok = (pkey[po][pos] == k_) & (r_col >= 0)
    rows, r_ev, r_mod, r_col, pos = rows[ok], r_ev[ok], r_mod[ok], r_col[ok], po[pos[ok]]
    is_t, is_s = sel_t[rows], sel_s[rows]
    v = valid[pos]
    own = a_row[pos]
    nother = n_ev[r_ev] - v.astype(int)
    com = np.where(nother >= 3, (S_ev[r_ev] - np.where(v, own, 0.0)) / np.maximum(nother, 1), np.nan)
    Y0 = praw[pos].astype(f32)
    use = v & np.isfinite(com)
    Yo = Y0 - own[:, None].astype(f32) * SGN.astype(f32)
    Yc = Y0 - np.where(use, com, 0.0)[:, None].astype(f32) * SGN.astype(f32)
    n = len(rows)
    log(f"[fix] {run}/{split}: {n} rows ({is_t.sum()} timing, {is_s.sum()} singles); valid {v.mean() * 100:.1f} %, "
        f"common-tone available {use.mean() * 100:.1f} %, {time.time() - t0:.0f} s")
    # ---- MF re-run
    p = C.setup_pipeline(evio_path=os.path.join(trd.pr.EVIO_DIR, f"prad_{run:06d}.evio.{split:05d}"), hc_calib_file=str(z["calib_file"]))
    hy = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    wcfg = p.wave_ana.cfg
    bank, _ = mf_pulse.build_bank([m.name for m in hy], [bool(m.is_pwo4()) for m in hy], tmpl,
                                  dict(smooth_order=int(wcfg.smooth_order), int_tail_ratio=float(wcfg.int_tail_ratio),
                                       tail_break_n=int(wcfg.tail_break_n)),
                                  column=mf_pulse.resolve_column(tmpl, None, prefer="25320"), int_tail_ratio=None)
    rec = mf_pulse.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4, r_model=0.3, ped_nsamples=int(wcfg.ped_nsamples))
    overflow = int(wcfg.overflow)
    ti = np.flatnonzero(is_t)

    def mf_dom(Y):
        t, E, okk = np.full(n, np.nan), np.full(n, np.nan), np.zeros(n, bool)
        for s0 in range(0, len(ti), 20000):
            ii = ti[s0:s0 + 20000]
            res = rec.run(Y[ii], r_mod[ii], overflow)
            pc = res.p_ch
            if len(pc) == 0:
                continue
            Eu = mf_pulse.pulse_energy_units(bank, r_mod[ii][pc], res.p_t0, res.p_amp, 1.0, None)
            order = np.lexsort((-Eu, pc))
            st = np.searchsorted(pc, np.arange(len(ii) + 1))
            has = st[1:] > st[:-1]
            fi = order[st[:-1][has]]
            r_ = np.flatnonzero(has)
            t[ii[r_]] = res.p_time[fi]
            E[ii[r_]] = [hy[int(r_mod[ii[q]])].energize(float(Eu[j])) for q, j in zip(r_, fi)]
            okk[ii[r_]] = ~res.wa_fb[r_]
        return t, E, okk

    V = {"raw": mf_dom(Y0), "own": mf_dom(Yo), "common": mf_dom(Yc)}
    log(f"[fix] {run}/{split}: MF done, {time.time() - t0:.0f} s")
    rf = dict(np.load(HERE / "output/module_timing_20261025/rf" / f"{key(run, split)}.npz"))
    ev_num_rows = meta["ev_num"][r_ev]
    base = is_t & use
    for k in V:
        base &= V[k][2] & np.isfinite(V[k][0])
    base &= V["raw"][1] >= 100.0
    phi, P = {}, None
    for k in V:
        phi[k], P, okp = _rf_phase(rf, ev_num_rows, np.where(base, V[k][0], 0.0), P)
        base &= okp
    NBINS = RF.NB
    cnt = np.bincount(r_col[base], minlength=NCH)
    mods = np.flatnonzero(cnt >= 100)
    tabt = []
    for m in mods:
        k = base & (r_col == m)
        row = [m, cnt[m]]
        for v_ in V:
            h = np.bincount(np.clip(np.floor((phi[v_][k] / P + .5) * NBINS).astype(np.int64), 0, NBINS - 1), minlength=NBINS)
            fr = RF.fit_unwrapped((h, P))
            row += [fr["sigma"] if fr["ok"] else np.nan, 1.4826 * RF.circ_mad(phi[v_][k], P)]
        tabt.append(row)
    tabt = np.array(tabt, float).reshape(-1, 8)       # mod, N, (fit, MAD) x (raw, own, common)
    pool = {}
    kk = base & (cnt[r_col] >= 20)
    for v_ in V:
        cen = np.zeros(NCH)
        for m in np.flatnonzero(cnt >= 20):
            cen[m] = RF.circ_median(phi[v_][base & (r_col == m)], P)
        pphi = RF.fold(phi[v_][kk] - cen[r_col[kk]], P)
        h = np.bincount(np.clip(np.floor((pphi / P + .5) * NBINS).astype(np.int64), 0, NBINS - 1), minlength=NBINS)
        fr = RF.fit_unwrapped((h, P))
        pool[v_] = (fr["sigma"] if fr["ok"] else np.nan, 1.4826 * float(np.median(np.abs(pphi))), int(kk.sum()))
    dt_rms = {k: float(1.4826 * np.median(np.abs(V[k][0][base] - V["raw"][0][base]))) for k in ("own", "common")}
    dt_max = {k: float(np.percentile(np.abs(V[k][0][base] - V["raw"][0][base]), 99)) for k in ("own", "common")}
    # ---- singles dchi2
    si = np.flatnonzero(is_s & (Y0.max(1) < 4095))
    sg = trd.sig_mod(slug)
    sigp = np.array([sg[names_all[int(r_mod[i])]] for i in si])
    midx = np.array([ix[names_all[int(r_mod[i])]] for i in si], np.int64)
    trd._BANKS[slug] = bank_t
    res_d = {}
    with mp.get_context("fork").Pool(workers) as pl:
        for tag, Y in (("raw", Y0), ("own", Yo), ("common", Yc)):
            y, _ = tf.subtract_ped(bank_t, Y[si])
            tasks = [(slug, y[i:i + 40], midx[i:i + 40], sigp[i:i + 40]) for i in range(0, len(si), 40)]
            R = np.concatenate(pl.map(trd._fit_worker, tasks, chunksize=1))
            res_d[tag] = np.where(R[:, 13] == 1, R[:, 12], 0.0)
    log(f"[fix] {run}/{split}: {len(si)} singles fitted, {time.time() - t0:.0f} s")
    out = dict(run=run, split=split, P=P, n_rows=n, n_t=int(base.sum()), tab_t=tabt, dt_rms_own=dt_rms["own"], dt_rms_common=dt_rms["common"], dt_p99_own=dt_max["own"], dt_p99_common=dt_max["common"],
               pool_raw=np.array(pool["raw"]), pool_own=np.array(pool["own"]), pool_common=np.array(pool["common"]), n_single=len(si),
               single_A=cp["A1"][rows][si], single_use=use[si], thr=trd.PROD_THR.get(str(z["group"]), 97.87), group=str(z["group"]),
               **{"dchi_" + k: v for k, v in res_d.items()}, a_own_rms=float(own[v].std()), elapsed=time.time() - t0)
    f = SCRATCH / "fix" / f"{key(run, split)}.npz"
    f.parent.mkdir(parents=True, exist_ok=True)
    np.savez(str(f) + ".tmp.npz", **out)
    os.replace(str(f) + ".tmp.npz", f)


def stage_fix(a):
    for g, r, s in SPLITS:
        if g in FIX_SPLITS and (a.force or not (SCRATCH / "fix" / f"{key(r, s)}.npz").exists()):
            fix_split(r, s, max(a.workers, 12))


# ------------------------------------------------------------------ report --
def load_all():
    A = {}
    for g, r, s in SPLITS:
        f = SCRATCH / "analysis" / f"{key(r, s)}.npz"
        if f.exists():
            A[g] = dict(np.load(f, allow_pickle=True))
    L = dict(np.load(SCRATCH / "analysis" / "lms.npz", allow_pickle=True))
    return A, L


def stored_sigma(r, s, names):
    import twopulse_realdata as trd
    z = np.load(DUMP / f"{key(r, s)}.npz", allow_pickle=True)
    sg = trd.sig_mod(trd.slug_of(str(z["templates"])))
    return np.array([sg.get(str(n), np.nan) for n in names])


def tsv(path, header, rows):
    with open(path, "w") as fh:
        fh.write("\t".join(header) + "\n")
        for r in rows:
            fh.write("\t".join(f"{x:.5g}" if isinstance(x, (float, np.floating)) else str(x) for x in r) + "\n")


def debias_acf(c_hat, N):
    """Exact inverse of the window-mean-removal bias.  c_hat[k] (k < N, any normalisation) are the pooled estimates
    sum_i x_i x_{i+k} / (N-k) of mean-removed windows of length N.  With P = I - 11^T/N and C(c) the Toeplitz matrix of the
    stationary acf c, E[c_hat_k] = sum_i (P C P)_{i,i+k} / (N-k) is linear in c; solve it (min-norm lstsq: the constant
    c_k = const is the null direction, a per-event offset, which mean-removed windows cannot see).  Returns rho = c / c_0."""
    from scipy.linalg import toeplitz
    P = np.eye(N) - 1.0 / N
    M = np.zeros((N, N))
    for j in range(N):
        e = np.zeros(N)
        e[j] = 1.0
        Q = P @ toeplitz(e) @ P
        M[:, j] = [np.trace(Q, offset=k) / (N - k) for k in range(N)]
    c = np.linalg.lstsq(M, np.asarray(c_hat[:N], float), rcond=1e-9)[0]
    return c / c[0]


def kernel_like(rho, kern=(2 / 3, 1.0, 2 / 3), N=NW):
    """sqrt of (expected per-window mean-removed variance of the smoothed buffer / that of the raw buffer) for the
    stationary acf rho (interior rows use the 3-tap kernel, sample 0 stays raw as in WaveAnalyzer::smooth)."""
    from scipy.linalg import toeplitz
    P = np.eye(N) - 1.0 / N
    r = np.zeros(N + 2)
    r[:len(rho)] = rho[:N + 2]
    C = toeplitz(r)
    k = np.asarray(kern) / np.sum(kern)
    S = np.zeros((N, N + 2))
    S[0, 0] = 1.0
    for i in range(1, N):
        S[i, i - 1:i + 2] = k
    num = np.trace(P @ S @ C @ S.T @ P)
    den = np.trace(P @ C[:N, :N] @ P)
    return math.sqrt(num / den)


def kernel_ratio(rho, kern):
    """sigma_smoothed / sigma_raw for stationary noise with acf rho[0..] (rho[0] = 1) and the smoothing kernel."""
    k = np.asarray(kern, float)
    k = k / k.sum()
    tot = sum(k[i] * k[j] * (rho[abs(i - j)] if abs(i - j) < len(rho) else 0.0) for i in range(len(k)) for j in range(len(k)))
    return math.sqrt(tot)


def toeplitz_cov(rho, n=NS):
    from scipy.linalg import toeplitz
    r = np.zeros(n)
    r[:min(n, len(rho))] = rho[:n]
    C = toeplitz(r)
    ev, V = np.linalg.eigh(C)
    return (V * np.clip(ev, 1e-3, None)) @ V.T


def template_numbers(rho, sig_raw, sig_stored, a_adc):
    """Noise on the template amplitude and onset for the measured covariance (acf rho, std sig_raw), unweighted LS fit
    (2.24 GeV H2 bank, 20 modules, 4 onset phases, medians): effective diagonal sigma of the amplitude (ls_eff) and onset
    sigmas (ns) at amplitude a_adc for a white sigma (stored / raw) and for the measured covariance (st_ls_cov)."""
    import twopulse_fit as tf
    import twopulse_realdata as trd
    bank = trd.get_bank("2p24GeV_H2")
    C = toeplitz_cov(rho)
    out = []
    for mi in np.linspace(0, len(bank.names) - 1, 20).astype(int):
        for ph in (0.0, 1.0, 2.0, 3.0):
            t0 = 4 * 34 - bank.t_pk[mi] + ph
            s, d = tf.tpl(bank.lut[mi].astype(f64), bank.dlut[mi].astype(f64), tf.TS - t0, bank.T_end, True)
            ss, dd = float(s @ s), float(d @ d)
            if ss <= 0 or dd <= 0:
                continue
            out.append((math.sqrt(float(s @ C @ s) / ss), math.sqrt(float(d @ C @ d) / dd), 1.0 / math.sqrt(dd)))
    o = np.nanmedian(np.array(out), 0)
    return dict(ls_eff=sig_raw * o[0], ls_fac=o[0], ls_fac_onset=o[1], st_white_stored=sig_stored / a_adc * o[2],
                st_white_raw=sig_raw / a_adc * o[2], st_ls_cov=sig_raw / a_adc * o[2] * o[1])


def pair_stats(x, pre, tag, nmin=300):
    """Mean / std / count of the pair correlations with >= nmin co-present samples: within board, other board of the
    same crate, other crate."""
    iu = np.triu_indices(NCH, 1)
    c = x[pre + "corr" + tag][iu].astype(f64)
    good = np.isfinite(c) & (x["np_mat"][iu] >= nmin)
    bd, cr = x["board"], x["crate"]
    sb = (bd[:, None] == bd[None, :])[iu]
    sc = (cr[:, None] == cr[None, :])[iu]
    out = []
    for m in (good & sb, good & ~sb & sc, good & ~sc):
        out += [float(c[m].mean()) if m.any() else np.nan, float(c[m].std()) if m.any() else np.nan, int(m.sum())]
    return out


def common_fraction(x, pre):
    """Per-board common-variance fraction f from the LOO regression R2 = f^2 / (f + (1-f) <1/(n-1)>)."""
    out = []
    for b in range(len(x["ub"])):
        n_ = x[pre + "cm_n"][b]
        if n_ < 100:
            continue
        s0, s1, wm, mm = x[pre + "cm_ssq0"][b], x[pre + "cm_ssq1"][b], x[pre + "cm_wm"][b], x[pre + "cm_mm"][b]
        r2, inv = wm * wm / (mm * s0), x[pre + "cm_inv"][b] / n_
        f = (r2 * (1 - inv) + math.sqrt((r2 * (1 - inv)) ** 2 + 4 * r2 * inv)) / 2
        out.append((int(x["ub"][b]), int(n_), math.sqrt(s0 / n_ / NW), math.sqrt(s1 / n_ / NW), 1 - s1 / s0, r2, wm / mm,
                    math.sqrt(1 + inv), f))
    return out


def stage_report(a):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    OUT.mkdir(parents=True, exist_ok=True)
    A, L = load_all()
    G = [g for g, _, _ in SPLITS if g in A]
    EP = [g for g in ("0p73", "2p24", "3p49") if g in A]
    TA = [g for g in G if g not in EP]
    R = {}
    ALL = dict(A)
    ALL["lms"] = L
    nm = lambda g: LAB.get(g, "LMS events")
    mean_ep = lambda f: np.mean([f(A[g]) for g in EP], 0)
    # ---- 1. sigma table and decomposition
    rows, dec, per = [], [], {}
    for g, r, s in SPLITS:
        if g not in A:
            continue
        x = A[g]
        st = stored_sigma(r, s, x["names"])
        per[g] = st
        for c in range(NCH):
            rows.append((LAB[g], x["names"][c], int(x["board"][c]), int(x["nwin"][c]), x["pooled"][c], x["r_pooled"][c], x["med_rstd0"][c],
                         x["med_rmad"][c], x["med_sstd0"][c], x["med_smad"][c], x["med_war"][c], st[c],
                         x["pooled"][c] / st[c], x["med_war"][c] / st[c]))
        mk = np.isfinite(st) & np.isfinite(x["med_war"]) & np.isfinite(x["pooled"])
        dec.append([LAB[g], int(x["n_valid30"]), *[float(np.median(x[k][mk])) for k in ("med_rstd0", "med_rmad", "med_sstd0", "med_smad", "med_war")],
                    float(np.median(st[mk])), float(np.median(x["pooled"][mk])), float(np.median(x["pooled"][mk] / st[mk])),
                    float(np.median(x["med_war"][mk] / st[mk])), float(np.median(x["r_pooled"][mk]))])
    tsv(OUT / "noise_sigma_channels.tsv",
        ["split", "module", "board", "n_windows", "sigma_raw_pooled_ddof1", "sigma_after_alt_removal", "med_raw_std_ddof0", "med_raw_1.4826MAD",
         "med_smoothed_std", "med_smoothed_1.4826MAD", "med_WA_ped_rms_replica", "stored_sig_ped", "ratio_raw_pooled/stored",
         "ratio_WAreplica/stored"], rows)
    tsv(OUT / "noise_sigma_decomposition.tsv",
        ["split", "n_clean_windows", "raw_std_ddof0", "raw_1.4826MAD", "smoothed_std", "smoothed_1.4826MAD", "WA_clip_rms_replica",
         "stored_sig_ped", "raw_pooled_ddof1", "ratio_raw_pooled/stored", "ratio_WAreplica/stored", "sigma_after_alt_removal"], dec)
    R["dec"] = dec
    rat = lambda gs: np.concatenate([A[g]["pooled"] / per[g] for g in gs])
    ra = rat(EP)
    ra = ra[np.isfinite(ra)]
    R["ratio_ep"] = [float(np.percentile(ra, p)) for p in (5, 50, 95)]
    ra = rat(G)
    ra = ra[np.isfinite(ra)]
    R["ratio_all"] = [float(np.percentile(ra, p)) for p in (5, 50, 95)]
    R["n_clean"] = {g: [int(A[g]["n_events"]), int(A[g]["n_rows"]), int(A[g]["n_valid30"]), int(A[g]["n_valid48"])] for g in G}
    R["n_lms"] = int(L["n_lms"])
    R["lms_wave"] = [float(v) for v in L["mean_wave"][28:44]]
    hr = np.sum([A[g]["hist_rmad"] for g in EP], 0)
    R["mad_quant_frac"] = float(hr[np.searchsorted(HIST_BINS, 1.4826) - 1] / hr.sum())
    # acf and kernel prediction
    rho_raw = mean_ep(lambda x: x["acf30"])                                # per-window mean-removed (biased), 30 lags
    rho = debias_acf(rho_raw, NW)                                           # stationary within-event acf, 30 lags
    rho_ta = debias_acf(np.mean([A[g]["acf30"] for g in TA], 0), NW)
    rho_l = debias_acf(L["acf30"], NW)
    rho_r = debias_acf(mean_ep(lambda x: x["r_acf30"]), NW)
    rho48_raw = mean_ep(lambda x: x["acf48"])
    rho48 = debias_acf(rho48_raw, NL)
    kern = [2 / 3, 1.0, 2 / 3]
    sr = float(np.median([d[2] for d in dec[:len(EP)]]))
    sm_meas = float(np.median([d[4] for d in dec[:len(EP)]]))
    kl = kernel_like(rho)
    R.update(acf_ep_raw=[float(v) for v in rho_raw[:12]], acf_ep=[float(v) for v in rho[:12]], acf_ep48=[float(v) for v in rho48[:12]],
             acf_ta=[float(v) for v in rho_ta[:6]], acf_lms=[float(v) for v in rho_l[:6]], acf_res=[float(v) for v in rho_r[:6]],
             acf_sum_raw=float(1 + 2 * sum(rho_raw[1:12])), acf_sum=float(1 + 2 * sum(rho[1:12])), acf_sum_full=float(1 + 2 * sum((1 - k / NW) * rho[k] for k in range(1, NW))),
             km_white=kernel_ratio([1.0], kern), km_acf=kernel_ratio(rho, kern), km_like=kl,
             km_meas=float(np.median([d[4] / d[2] for d in dec[:len(EP)]])), sig_raw_std_ep=sr, sig_smoothed_meas=sm_meas,
             sig_pred_smoothed=kernel_ratio(rho, kern) * sr, sig_pred_smoothed_like=kl * sr)
    tsv(OUT / "noise_acf.tsv", ["lag"] + ["acf30_raw_" + LAB[g] for g in G] + ["acf30_raw_epmean", "acf30_debiased_epmean", "acf30_alt_removed_debiased_epmean",
                                                                              "acf30_debiased_Ta", "acf48_late_only_raw_ep", "acf48_late_only_debiased_ep", "acf30_debiased_LMS"],
        [(k,) + tuple(A[g]["acf30"][k] for g in G) + (rho_raw[k], rho[k], rho_r[k], rho_ta[k], rho48_raw[k], rho48[k], rho_l[k]) for k in range(NW)])
    # ---- 2. spectra
    fr30, fr48 = np.arange(NF30) * FS / NW, np.arange(NF48) * FS / NL

    def nspec(S, n):
        s = S / np.maximum(n, 1)[:, None]
        ok = n >= 5
        return s[ok] / s[ok][:, 1:-1].mean(1, keepdims=True), np.flatnonzero(ok)

    spec30 = {g: nspec(ALL[g]["S30"], ALL[g]["nwin"]) for g in G + ["lms"]}
    spec30r = {g: nspec(ALL[g]["r_S30"], ALL[g]["r_nwin"]) for g in G + ["lms"]}
    spec48 = {g: nspec(A[g]["S48"], A[g]["nwin48"]) for g in G}
    ov30 = np.mean([spec30[g][0].mean(0) for g in EP], 0)
    ov30r = np.mean([spec30r[g][0].mean(0) for g in EP], 0)
    ov30t = np.mean([spec30[g][0].mean(0) for g in TA], 0)
    s48e, _ = nspec(sum(A[g]["S48"] for g in EP), sum(A[g]["nwin48"] for g in EP))
    s48t, _ = nspec(sum(A[g]["S48"] for g in TA), sum(A[g]["nwin48"] for g in TA))
    ov48, ov48t = s48e.mean(0), s48t.mean(0)
    w48, w30 = np.ones(NF48), np.ones(NF30)
    w48[[0, -1]] = 0.5
    w30[[0, -1]] = 0.5
    bd = A[G[0]]["board"]
    brd = []
    with open(OUT / "noise_spectra.tsv", "w") as fh:
        fh.write("# periodogram / channel band mean (white = 1). 30: samples 0..29 rect window, all clean channel-events; 48: late-only channels (Hann)\n")
        fh.write("window\tfreq_MHz\t" + "\t".join(LAB[g] for g in G) + "\tLMS\te-p mean\te-p alt-removed\tTa mean\n")
        for k, f_ in enumerate(fr30):
            fh.write(f"30\t{f_:.2f}\t" + "\t".join(f"{spec30[g][0][:, k].mean():.4f}" for g in G + ["lms"]) + f"\t{ov30[k]:.4f}\t{ov30r[k]:.4f}\t{ov30t[k]:.4f}\n")
        fh.write("window\tfreq_MHz\te-p pooled\tTa pooled\n")
        for k, f_ in enumerate(fr48):
            fh.write(f"48\t{f_:.2f}\t{ov48[k]:.4f}\t{ov48t[k]:.4f}\n")
        fh.write("# common-mode spectrum, boards with >= 5 clean channels (periodogram, ADC^2 per bin; e-p mean): freq, LOO board mean, LOO residual, channel\n")
        pc = mean_ep(lambda x: x["Pcm"])
        for k, f_ in enumerate(fr30):
            fh.write(f"CM\t{f_:.2f}\t{pc[0][k]:.4f}\t{pc[1][k]:.4f}\t{pc[2][k]:.4f}\n")
        fh.write("# per-board mean spectrum, 30-sample window, e-p physics rows\nboard\t" + "\t".join(f"{f_:.1f}" for f_ in fr30) + "\n")
        for b in np.unique(bd):
            vs = [spec30[g][0][np.isin(spec30[g][1], np.flatnonzero(bd == b))].mean(0) for g in EP if np.isin(spec30[g][1], np.flatnonzero(bd == b)).any()]
            if vs:
                brd.append(np.mean(vs, 0))
                fh.write(f"{b}\t" + "\t".join(f"{x:.4f}" for x in brd[-1]) + "\n")
    brd = np.array(brd)
    tv = lambda s_: float(1 - 0.5 * np.mean(np.abs(s_[1:-1] / s_[1:-1].mean() - 1)))
    tot = (ov30 * w30).sum()
    wref = np.ones(NF30)
    wref[0] = 0                      # white, mean-removed 30-sample window: DC bin empty, bins 1..15 flat
    mom = np.sum([A[g]["mom"] for g in EP], 0)
    pool2 = np.mean([np.nanmedian(A[g]["pooled"]) ** 2 * (NW - 1) / NW for g in EP])
    R.update(psd30_ep=[float(v) for v in ov30], psd30_ta=[float(v) for v in ov30t], psd30_res=[float(v) for v in ov30r], psd48_ep=[float(v) for v in ov48],
             psd48_ta=[float(v) for v in ov48t], flat30=float(np.exp(np.mean(np.log(ov30[1:-1]))) / ov30[1:-1].mean()), tv30=tv(ov30), tv30_res=tv(ov30r),
             hi30=float((ov30 * w30)[fr30 > 30].sum() / tot), hi30_white=float((wref * w30)[fr30 > 30].sum() / (wref * w30).sum()), lowf=float((ov30 * w30)[fr30 < 10].sum() / tot),
             lowf_ta=float((ov30t * w30)[fr30 < 10].sum() / (ov30t * w30).sum()), lowf_white=float((wref * w30)[fr30 < 10].sum() / (wref * w30).sum()), nyq_share=float(0.5 * ov30[-1] / tot * 1.0),
             nyq_peak=float(ov30[-1]), board_spread=float(np.std(brd[:, 1:-1] / brd[:, 1:-1].mean(1, keepdims=True), axis=0).max()),
             n48=[int(sum(A[g]["nwin48"].sum() for g in EP)), int(sum(A[g]["nwin48"].sum() for g in TA))],
             alt_rms=float(math.sqrt(mom[2] / mom[0])), alt_var_frac=float(mom[2] / mom[0] / pool2), alt_kurt=float(mom[3] / mom[0] / (mom[2] / mom[0]) ** 2),
             alt_mean=float(mom[1] / mom[0]), coh_all=float(np.mean([A[g]["coh"][0] / A[g]["coh"][1] for g in EP])),
             coh_board=float(np.mean([A[g]["coh"][2] / A[g]["coh"][3] for g in EP])),
             persist_ratio=float(np.mean([A[g]["per"][2] / A[g]["per"][1] for g in EP])), persist_corr=float(np.mean([A[g]["per"][3] / np.sqrt(A[g]["per"][4] * A[g]["per"][5]) for g in EP])),
             cm_psd=[[float(v) for v in row] for row in pc])
    dv = ov30[1:-1] / (0.5 * (ov30[:-2] + ov30[2:]))
    R["line625"] = float(ov48[12] / (0.5 * (ov48[11] + ov48[13])))
    R["line_dev30"] = [float(np.max(np.abs(dv[:-1] - 1))), float(fr30[1:-1][np.argmax(np.abs(dv[:-1] - 1))])]
    f_b = 1e3 / 4.0077
    R["alias"] = [(h, float(((h * f_b + FS / 2) % FS) - FS / 2)) for h in (1, 2, 3, 4)]
    R["rf_test"] = {g: json.loads(str(A[g]["rf_test"])) for g in G}
    pat_rms = float(np.sqrt(np.nanmean(np.concatenate([(A[g]["pattern"] ** 2).mean(1) for g in EP]))))
    exp_pat = float(np.sqrt(np.nanmean(np.concatenate([(A[g]["pooled"] ** 2) / np.maximum(A[g]["nwin"], 1) * (1 - 1 / NW) for g in EP]))))
    R.update(pat_rms=pat_rms, pat_exp=exp_pat)
    # ---- 3. correlation
    crow, CR = [], {}
    for g in G + ["lms"]:
        x = ALL[g]
        for pre, lab in (("", "raw"), ("r_", "alt-removed")):
            for tag in ("A", "B"):
                st = pair_stats(x, pre, tag)
                crow.append((nm(g), lab, "offset removed" if tag == "A" else "offset kept", *st))
                CR[(g, pre, tag)] = st
    tsv(OUT / "noise_correlation.tsv",
        ["sample", "waveform", "baseline", "within_board_mean", "within_board_std", "n_pairs", "same_crate_other_board_mean", "std", "n_pairs",
         "other_crate_mean", "std", "n_pairs"], crow)
    R["corr"] = {}
    for pre in ("", "r_"):
        for tag in ("A", "B"):
            R["corr"][pre + tag] = dict(ep=[float(np.nanmean([CR[(g, pre, tag)][i] for g in EP])) for i in (0, 1, 3, 6)],
                                        ta=[float(np.nanmean([CR[(g, pre, tag)][i] for g in TA])) for i in (0, 1, 3, 6)],
                                        lms=[float(CR[("lms", pre, tag)][i]) for i in (0, 1, 3, 6)],
                                        npairs=[int(np.mean([CR[(g, pre, tag)][i] for g in EP])) for i in (2, 5, 8)])
    g30 = 1 + 2 * sum((1 - k / NW) * rho[k] for k in range(1, NW))
    g30_raw = g30
    g30 = max(g30, 1e-3)               # DC power of the stationary acf is ~0 (gauge: min-norm); floor so the expectation stays defined
    drows = []
    for g in G + ["lms"]:
        x = ALL[g]
        lo = x["pooled"] / math.sqrt(NW)
        drows.append((nm(g), float(np.nanmedian(x["drift"][:, 0])), float(np.nanmedian(x["drift"][:, 1])),
                      float(np.nanmedian(lo)), float(np.nanmedian(lo) * math.sqrt(g30)), float(np.nanmedian(x["drift"][:, 0] / lo)),
                      float(np.nanmedian(x["drift"][:, 1] / (lo * math.sqrt(g30))))))
    tsv(OUT / "noise_drift.tsv", ["sample", "rms_between_events", "highpass_rms_consecutive_rows", "sigma/sqrt30", "sigma/sqrt30*sqrt(acf factor)",
                                  "between/white", "highpass/acf-corrected"], drows)
    R["drift"] = dict(ep=[float(np.mean([d[i] for d in drows[:len(EP)]])) for i in (1, 2, 3, 4, 5, 6)], lms=[float(v) for v in drows[-1][1:]], acf_factor=float(g30_raw), acf_factor_used=float(g30))
    # ---- 4. common mode
    cmrows, CMD = [], {}
    for g in G + ["lms"]:
        for pre, lab in (("", "raw"), ("r_", "alt-removed")):
            cr_ = common_fraction(ALL[g], pre)
            CMD[(g, pre)] = np.array([r[2:] for r in cr_], float)
            for r in cr_:
                cmrows.append((nm(g), lab, *r))
    tsv(OUT / "noise_commonmode_boards.tsv",
        ["sample", "waveform", "board", "n_channel_windows", "sigma_before", "sigma_after_LOO", "var_removed_frac(direct)", "R2_regression", "beta",
         "expected_after/before_if_independent", "common_variance_fraction_f(corrected for 1/(n-1))"], cmrows)
    R["cm"] = {}
    for pre in ("", "r_"):
        for grp, gs in (("ep", EP), ("ta", TA), ("lms", ["lms"])):
            c = np.concatenate([CMD[(g, pre)] for g in gs])
            R["cm"][pre + grp] = dict(s0=float(np.median(c[:, 0])), s1=float(np.median(c[:, 1])), gain=float(np.median(c[:, 2])), r2=float(np.median(c[:, 3])),
                                      beta=float(np.median(c[:, 4])), expect=float(np.median(c[:, 5])), f=float(np.median(c[:, 6])),
                                      f_p10=float(np.percentile(c[:, 6], 10)), f_p90=float(np.percentile(c[:, 6], 90)), n=len(c))
    # ---- 5. template numbers (e-p covariance)
    z = np.load(DUMP / f"{key(*SPLITS[1][1:])}.npz", allow_pickle=True)
    big = (z["pu_E"] >= 100) & (~z["pu_fb"].astype(bool)) & (z["pu_t"] > 100) & (z["pu_t"] < 200)
    a_adc = float(np.median(z["pu_A"][big]))
    sig_raw = float(np.nanmedian(np.concatenate([A[g]["pooled"] for g in EP])))
    sig_st = float(np.nanmedian(np.concatenate([per[g] for g in EP])))
    R["tpl"] = template_numbers(rho, sig_raw, sig_st, a_adc)
    R["tpl"].update(a_adc=a_adc, sig_raw=sig_raw, sig_st=sig_st)
    # ---- 6. fix test on real pulses
    F = {}
    for g, r, s in SPLITS:
        f = SCRATCH / "fix" / f"{key(r, s)}.npz"
        if f.exists():
            F[g] = dict(np.load(f, allow_pickle=True))
    trows, drows2 = [], []
    for g in F:
        c = F[g]
        t = c["tab_t"]
        okm = np.all(np.isfinite(t[:, [2, 4, 6]]), 1)
        trows.append((LAB[g], int(c["n_t"]), len(t), int(okm.sum()), *[float(np.median(t[okm, j])) for j in (2, 4, 6)],
                      *[float(np.median(t[:, j])) for j in (3, 5, 7)], *[float(c["pool_" + k][0]) for k in ("raw", "own", "common")],
                      *[float(c["pool_" + k][1]) for k in ("raw", "own", "common")], float(c["dt_rms_own"]), float(c["dt_rms_common"]), float(c["dt_p99_own"]), float(c["dt_p99_common"])))
        for tag in ("raw", "own", "common"):
            d = c["dchi_" + tag]
            drows2.append((LAB[g], tag, int(c["n_single"]), float(np.median(d)), float(np.percentile(d, 97)), float(np.percentile(d, 99)),
                           float(np.mean(d > c["thr"]) * 100), float(c["thr"])))
    tsv(OUT / "noise_fix_sigma_t.tsv",
        ["split", "n_rows_E>=100", "n_modules_N>=100", "n_modules_fit_ok", "sigma_t_fit_med_raw_ns", "sigma_t_fit_med_own_ns", "sigma_t_fit_med_common_ns",
         "sigma_t_MAD_med_raw_ns", "sigma_t_MAD_med_own_ns", "sigma_t_MAD_med_common_ns", "pooled_fit_raw_ns", "pooled_fit_own_ns", "pooled_fit_common_ns",
         "pooled_MAD_raw_ns", "pooled_MAD_own_ns", "pooled_MAD_common_ns", "dt_own_vs_raw_1.4826MAD_ns", "dt_common_vs_raw_1.4826MAD_ns", "dt_own_vs_raw_p99_ns", "dt_common_vs_raw_p99_ns"], trows)
    tsv(OUT / "noise_fix_dchi2.tsv", ["split", "waveform", "n_singles", "dchi2_median", "dchi2_p97", "dchi2_p99", "pct_above_prod_thr", "prod_thr"], drows2)
    R["fix_t"], R["fix_d"] = trows, drows2
    # ---- PDF
    colr = dict(zip(G, ["#1f77b4", "#d62728", "#2ca02c", "#9467bd", "#ff7f0e"]))
    colr["lms"] = "k"
    with PdfPages(OUT / "noise_diag.pdf") as pdf:
        npg = [0]

        def sv(fig):
            pdf.savefig(fig)
            if os.environ.get("NOISE_DIAG_PNG"):
                fig.savefig(f"{os.environ['NOISE_DIAG_PNG']}/p{npg[0]}.png", dpi=60)
            npg[0] += 1

        fig, ax = plt.subplots(2, 2, figsize=(11, 8))
        for g in G:
            ax[0, 0].hist(A[g]["pooled"], np.linspace(0.5, 1.5, 81), histtype="step", color=colr[g], label=LAB[g])
            ax[0, 1].hist(per[g], np.linspace(0.2, 0.7, 81), histtype="step", color=colr[g])
            ax[1, 0].hist(A[g]["pooled"] / per[g], np.linspace(0.8, 4.0, 81), histtype="step", color=colr[g])
        ax[0, 0].hist(L["pooled"], np.linspace(0.5, 1.5, 81), histtype="step", color="k", ls="--", label="LMS events")
        ax[0, 0].set_xlabel("measured raw baseline sigma per module (ADC)"); ax[0, 0].legend(fontsize=7)
        ax[0, 0].set_title("pre-pulse window 0..29, pooled std (ddof 1)")
        ax[0, 1].set_xlabel("stored sig_ped per module (ADC)"); ax[0, 1].set_title("stored table")
        ax[1, 0].set_xlabel("measured raw sigma / stored sig_ped"); ax[1, 0].set_title("ratio per module")
        labs = ["raw std\n(ddof 0)", "raw\n1.4826 MAD", "smoothed\nstd", "smoothed\n1.4826 MAD", "WA clip rms\n(replica)", "stored\nsig_ped"]
        for gi, g in enumerate(G):
            d = [r for r in dec if r[0] == LAB[g]][0]
            ax[1, 1].bar(np.arange(6) + 0.15 * (gi - 2), [d[i] for i in (2, 3, 4, 5, 6, 7)], 0.15, color=colr[g], label=LAB[g])
        ax[1, 1].set_xticks(range(6)); ax[1, 1].set_xticklabels(labs, fontsize=7)
        ax[1, 1].set_ylabel("median over modules (ADC)"); ax[1, 1].set_title("same windows, each definition")
        fig.suptitle("Baseline sigma: measured vs stored"); fig.tight_layout(); sv(fig); plt.close(fig)
        fig, ax = plt.subplots(1, 3, figsize=(13, 4))
        for k, lab in (("rstd0", "raw std"), ("rmad", "raw 1.4826 MAD"), ("sstd0", "smoothed std"), ("smad", "smoothed 1.4826 MAD"), ("war", "WA clip rms")):
            h = np.sum([A[g]["hist_" + k] for g in EP], 0)
            ax[0].step(HIST_BINS[:-1], h / h.sum() / (HIST_BINS[1] - HIST_BINS[0]), where="post", label=lab)
        ax[0].set_xlim(0, 2.5); ax[0].set_xlabel("per-window estimate (ADC)"); ax[0].set_ylabel("density"); ax[0].legend(fontsize=7)
        ax[0].set_title("per-window estimators (30 samples), e-p splits", fontsize=9)
        for g in G:
            ax[1].plot(range(KACF), A[g]["acf30"][:KACF], "o-", ms=3, color=colr[g], label=LAB[g] + " raw")
        ax[1].plot(range(KACF), rho[:KACF], "k-", lw=2, label="e-p debiased")
        ax[1].plot(range(KACF), rho_r[:KACF], "k--", label="e-p debiased, alt. part removed")
        ax[1].plot(range(KACF), rho_l[:KACF], "k:", label="LMS debiased")
        ax[1].axhline(0, color="gray", lw=.5); ax[1].set_xlabel("lag (samples)"); ax[1].set_ylabel("noise acf"); ax[1].legend(fontsize=6)
        ax[1].set_title("pooled baseline autocorrelation (mean removed)", fontsize=9)
        ax[2].plot(fr30[1:], ov30[1:], "o-", color="C1", ms=3, label="e-p, 30-sample")
        ax[2].plot(fr30[1:], ov30r[1:], "s--", color="C2", ms=3, label="e-p, alternating part removed")
        ax[2].plot(fr30[1:], ov30t[1:], "^-", color="C4", ms=3, label="Ta, 30-sample")
        ax[2].plot(fr48, ov48, "k-", lw=1, label="e-p late-only 48-sample")
        ax[2].axhline(1, color="gray", lw=.5); ax[2].set_xlabel("frequency (MHz)"); ax[2].set_ylabel("PSD / band mean (white = 1)")
        ax[2].set_yscale("log"); ax[2].legend(fontsize=6); ax[2].set_title("overall spectrum", fontsize=9)
        fig.tight_layout(); sv(fig); plt.close(fig)
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
        for g in G:
            ax[0].plot(fr30[1:], spec30[g][0].mean(0)[1:], "o-", ms=3, color=colr[g], label=LAB[g])
            ax[1].plot(fr48, spec48[g][0].mean(0), "o-", ms=3, color=colr[g], label=f"{LAB[g]} ({len(spec48[g][1])} ch)")
        ax[0].plot(fr30[1:], spec30["lms"][0].mean(0)[1:], "k--", label="LMS events")
        lo, hi = np.percentile(brd, [5, 95], axis=0)
        ax[0].fill_between(fr30[1:], lo[1:], hi[1:], color="gray", alpha=.25, label="boards 5-95 % (e-p)")
        for a_ in ax[:2]:
            a_.axhline(1, color="gray", lw=.6); a_.set_xlabel("frequency (MHz); Nyquist 125"); a_.set_ylabel("PSD / band mean"); a_.set_yscale("log")
        ax[0].legend(fontsize=6); ax[1].legend(fontsize=6)
        ax[0].set_title("pre-pulse window, samples 0..29 (rect)"); ax[1].set_title("late-only channels, samples 0..47 (Hann)")
        for q_, (lab, c_) in enumerate((("board mean of other channels (LOO)", "C3"), ("channel minus LOO mean", "C0"), ("channel", "k"))):
            ax[2].plot(fr30[1:], pc[q_][1:], "o-", ms=3, color=c_, label=lab)
        ax[2].set_yscale("log"); ax[2].set_xlabel("frequency (MHz)"); ax[2].set_ylabel("periodogram (ADC^2 per bin)"); ax[2].legend(fontsize=6)
        ax[2].set_title("common mode in boards with >= 5 clean channels (e-p)", fontsize=9)
        fig.suptitle("Baseline spectra per run type (RF alias 0.5-1 MHz is inside the first bin)")
        fig.tight_layout(); sv(fig); plt.close(fig)
        # correlation matrices: two boards, physics (2.24 H2) and LMS, raw / alt removed
        xp = A["2p24"] if "2p24" in A else A[G[0]]
        iu = np.triu_indices(NCH, 1)
        bd_ = xp["board"]
        score = {}
        for b in np.unique(bd_):
            ii = np.flatnonzero(bd_ == b)
            score[b] = int((xp["np_mat"][np.ix_(ii, ii)][np.triu_indices(len(ii), 1)] >= 300).sum()) if len(ii) >= 14 else -1
        pick = sorted(score, key=score.get)[-2:]
        fig, ax = plt.subplots(2, 4, figsize=(16, 8))
        for j, (b, x_, lab) in enumerate([(pick[0], xp, "2.24 GeV H2"), (pick[0], L, "LMS"), (pick[1], xp, "2.24 GeV H2"), (pick[1], L, "LMS")]):
            ii = np.flatnonzero(bd_ == b)
            ii = ii[np.argsort(xp["chan"][ii])]
            for i2, pre in enumerate(("", "r_")):
                m = x_[pre + "corrA"][np.ix_(ii, ii)].astype(f64)
                m[x_["np_mat"][np.ix_(ii, ii)] < 300] = np.nan
                im = ax[i2, j].imshow(m, vmin=-0.3, vmax=1, cmap="RdBu_r")
                ax[i2, j].set_title(f"board {b}, {lab}, {'raw' if pre == '' else 'alternating part removed'}\nmean off-diag {np.nanmean(m[np.triu_indices(len(ii), 1)]):+.3f}", fontsize=8)
        fig.colorbar(im, ax=ax, shrink=.6, label="corr coefficient (offset removed; blank: < 300 common samples)")
        sv(fig); plt.close(fig)
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
        for g in G + ["lms"]:
            x_ = ALL[g]
            c = x_["corrA"][iu].astype(f64)
            gd = x_["np_mat"][iu] >= 300
            sb = (x_["board"][:, None] == x_["board"][None, :])[iu]
            ax[0].hist(c[gd & sb], np.linspace(-.5, 1, 60), histtype="step", color=colr[g], label=nm(g), density=True)
            ax[1].hist(c[gd & ~sb], np.linspace(-.5, 1, 60), histtype="step", color=colr[g], density=True)
        ax[0].set_title("same board, pairs with >= 300 common samples", fontsize=9); ax[1].set_title("other boards", fontsize=9)
        for a_ in ax[:2]:
            a_.set_xlabel("pair corr (offset removed)"); a_.set_yscale("log")
        ax[0].legend(fontsize=6)
        dch_ = np.abs(xp["chan"][:, None] - xp["chan"][None, :])[iu]
        sbx = (xp["board"][:, None] == xp["board"][None, :])[iu]
        for g in EP + ["lms"]:
            x_ = ALL[g]
            gd = x_["np_mat"][iu] >= 300
            for pre, ls in (("", "-"), ("r_", "--")):
                c = x_[pre + "corrA"][iu].astype(f64)
                ax[2].plot(range(1, 16), [float(np.nanmean(c[gd & sbx & (dch_ == d)])) for d in range(1, 16)], ls, color=colr[g], lw=1)
        ax[2].set_xlabel("channel distance within board"); ax[2].set_ylabel("mean corr"); ax[2].set_title("within-board corr vs channel distance (solid raw, dashed alt removed)", fontsize=8)
        fig.tight_layout(); sv(fig); plt.close(fig)
        # common mode
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
        for g in EP + TA + ["lms"]:
            for pre, mk in (("", "o"), ("r_", "x")):
                c = CMD[(g, pre)]
                ax[0].scatter(c[:, 0], c[:, 1], s=8, marker=mk, color=colr[g], label=nm(g) if pre == "" else None)
        ax[0].plot([0.4, 4], [0.4, 4], "k-", lw=.6)
        ax[0].set_xlim(0.4, 4); ax[0].set_ylim(0.4, 4)
        ax[0].set_xlabel("sigma before (ADC)"); ax[0].set_ylabel("sigma after LOO board-mean subtraction (ADC)"); ax[0].legend(fontsize=6)
        ax[0].set_title("per board, samples 0..29 mean removed (o raw, x alternating part removed)", fontsize=8)
        for g in EP + ["lms"]:
            ax[1].hist(CMD[(g, "")][:, 6], np.linspace(0, 1, 40), histtype="step", color=colr[g], label=nm(g))
        ax[1].set_xlabel("common-variance fraction f per board (corrected)"); ax[1].legend(fontsize=6)
        for g in G + ["lms"]:
            ax[2].hist(ALL[g]["drift"][:, 0] / (ALL[g]["pooled"] / math.sqrt(NW)), np.linspace(0, 8, 80), histtype="step", color=colr[g])
        ax[2].axvline(1, color="k", lw=.6); ax[2].axvline(math.sqrt(g30), color="gray", ls=":")
        ax[2].set_xlabel("RMS of window mean across events / (sigma / sqrt30)"); ax[2].set_title("pedestal wander (dotted: acf-corrected white)", fontsize=9)
        fig.tight_layout(); sv(fig); plt.close(fig)
        # fixed-pattern / alternation / noise-limited onset error
        fig, ax = plt.subplots(1, 3, figsize=(14, 4.5))
        for g in EP:
            ax[0].plot(np.arange(NW), np.sqrt((A[g]["pattern"] ** 2).mean(0)), color=colr[g], label=LAB[g])
        ax[0].axhline(exp_pat, color="k", ls=":", label="expected for no pattern")
        ax[0].set_xlabel("sample"); ax[0].set_ylabel("rms over channels of event-averaged baseline (ADC)"); ax[0].legend(fontsize=6)
        ax[0].set_title("fixed pattern vs sample", fontsize=9)
        for g in EP:
            ax[1].hist(np.concatenate([[1]]) * 0 + (A[g]["oe_rms"] / (A[g]["pooled"] * math.sqrt(2 / 15))), np.linspace(0.5, 3.5, 80), histtype="step", color=colr[g])
        ax[1].set_xlabel("rms of (even-sample mean - odd-sample mean) / white expectation"); ax[1].set_title("alternating (125 MHz) component", fontsize=9)
        t_ = R["tpl"]
        ax[2].bar(range(3), [t_["st_white_stored"], t_["st_white_raw"], t_["st_ls_cov"]], color=["C0", "C1", "C3"])
        ax[2].set_xticks(range(3)); ax[2].set_xticklabels(["stored sigma,\nwhite", "raw sigma,\nwhite", "raw sigma,\nmeasured covariance"], fontsize=7)
        ax[2].set_ylabel(f"noise-only onset sigma at A = {a_adc:.0f} ADC (ns)")
        ax[2].set_title("onset error from noise (E >= 100 MeV, 2.24 GeV H2 bank)", fontsize=9)
        fig.tight_layout(); sv(fig); plt.close(fig)
        # linear fix on real pulses
        if F:
            fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
            gs = list(F)
            xs = np.arange(len(gs))
            for j, (k, c_) in enumerate((("raw", "C0"), ("own", "C1"), ("common", "C3"))):
                ax[0].bar(xs + 0.25 * (j - 1), [float(np.median(F[g]["dchi_" + k])) for g in gs], 0.25, color=c_, label=k)
                ax[1].bar(xs + 0.25 * (j - 1), [float(np.percentile(F[g]["dchi_" + k], 97)) for g in gs], 0.25, color=c_)
                ax[2].bar(xs + 0.25 * (j - 1), [float(F[g]["pool_" + k][1]) * 1e3 for g in gs], 0.25, color=c_)
            for a_, lab in zip(ax, ("single-pulse dchi2 median", "single-pulse dchi2 97th percentile", "pooled sigma_t vs RF (1.4826 MAD, ps)")):
                a_.set_xticks(xs); a_.set_xticklabels([LAB[g] for g in gs], fontsize=7); a_.set_ylabel(lab)
            ax[0].legend(fontsize=7, title="alternating part subtracted from", loc="lower right")
            ax[0].set_title("raw waveform / own pre-pulse fit / event-common estimate", fontsize=8)
            fig.tight_layout(); sv(fig); plt.close(fig)
    json.dump(R, open(OUT / "numbers.json", "w"), indent=1, default=float)
    log(json.dumps({k: v for k, v in R.items() if k not in ("dec", "psd30_ep", "psd30_ta", "psd30_res", "psd48_ep", "psd48_ta", "cm_psd")}, indent=1, default=float)[:9000])


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("stage", choices=["all", "decode", "analyse", "fix", "report", "decode-one"])
    ap.add_argument("--workers", type=int, default=6)
    ap.add_argument("--force", action="store_true")
    ap.add_argument("--run", type=int)
    ap.add_argument("--split", type=int)
    ap.add_argument("--lms-only", action="store_true")
    a = ap.parse_args(argv)
    if a.stage == "decode-one":
        return decode_one(a.run, a.split, a.lms_only)
    if a.stage in ("all", "decode"):
        stage_decode(a)
    if a.stage in ("all", "analyse"):
        stage_analyse(a)
    if a.stage in ("all", "fix"):
        stage_fix(a)
    if a.stage in ("all", "report"):
        stage_report(a)


if __name__ == "__main__":
    main()
