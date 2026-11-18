#!/usr/bin/env python3
"""exotic_p8_energy_error.py -- P8 of the exotic-events impact packages: the energy error of each exotic category, per channel and per cluster.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p8_energy_error.py [--jobs 12] [--out <dir>] [--cache <npz>] [--only <substr>]
                                                                                      [--examples 10] [--evio-dir <dir>]

Channels.  Every flagged PbWO4 channel (exclusive categories f4 > f5 > f1 > f2 > f3, uniform P2 thresholds, exotic_p4_module_rates.flags_of through
exotic_p7_cluster_level.row_categories' rule) that is a member of a PbWO4-seeded cluster of P7's join: the benchmark's gated clusters (control v3, CONTROL_GATE,
M1 phase factor, time offsets, control-v3 per-run gains), reclustered here for the events that hold a flagged row exactly as P7 pass 1 does, and matched to P7's
cluster cache (output/pileup_rates/exotic_p7_cache.npz) by (run, split, event number, seed module, E_g) for the populations (elastic window, Moller pairs) and
rings.  Data: H2 run types from the 2-file sample (output/pileup_rates/exotic_sample), Ta X17 from the 19 full-scan npz (output/pileup_rates/exotic_full); the
Ta sample (P7's cross-check row, 4 of its 6 splits inside the full set) is not refitted.  Run 24937 excluded.

Refit.  The channel's waveform is re-read from EVIO (decode only, the event's split) and fitted by pileup_rates.fit with the exotic scan's stored noise
reference of that split (default options; the stored k and chi2/dof must be reproduced: the 'repro' fraction), then by the alternative of its category:
  f1 tail-baseline : 'bs'  pileup_rates.fit with BASE_SLOPE (sloped baseline on the gated rows, pedestal > EX_DB sigma_ref above b_ref)
  f3 smooth-wide   : 'st'  pileup_rates.stretch_fit (one time-stretched template, s = 1.0-3.0) on the default fit
  f4 clipped-multi : 'cs'  pileup_rates.fit with CLIP_SINGLE (8 ns); rule "never drop the main onset": when that fit has no pulse within MATCH_NS of the default main
                           onset (k = 0 included), the fallback MF.fit_group on the overflow-masked waveform with the default main onset and the default secondaries
                           outside the clip stretch +- EX_CLIP_NS (the 'fallback' count)
  f2 wide          : 'k1'  one template pulse (MF.fit_group, k = 1, started at the default main onset), and 'st' the stretched single pulse (stretch_fit)
  f5 saturated     : 'k1'  one template pulse on the unclipped samples (overflow masked: the fit extrapolates the clipped top), and 'int' the measured integral:
                           raw sum def. 1 (pileup_rates.raw_sums, WaveAnalyzer window of the k1 model; clipped samples at the overflow value) / (c_rs I_m), I_m the
                           template integral, c_rs = the split's median raw sum / (A I_m) of unclipped k = 1 reference members (up to 200); a lower bound (the clipped area
                           is missing), the clipped-area fraction = 1 - E_int / E_k1
Energy.  A cluster member hit's energy is the clusterer's hit energy (benchmark h_e: calibration x I_m x phase factor x the benchmark matched-filter
amplitude h_a, times the clusterer's non-linearity factor where it is on) times the control-v3 gain, E_hit (as P7's E_g); g = E_hit / h_a (MeV per ADC amplitude) for matched-filter hits, the module's median g over the split for WaveAnalyzer hits (h_k = 0).
Each hit is matched to the default pileup_rates pulse nearest in time (h_t = onset + t50); the alternative amplitude of that pulse is: multi-pulse variants
(bs, cs) the alternative pulse matched one-to-one within MATCH_NS (0 if none: the pulse is gone); single-pulse variants (st, k1, int) the single pulse goes to
the module's in-window hit (any cluster) nearest its time (onset + s t50) when that time is inside the benchmark hit window (100, 200) ns, every other pulse of the
module gets 0 (a module split over two clusters is credited once; a hit outside every saved cluster takes the pulse out of the sums, a one-sided loss: JUDGE
page); a member whose single fit has no pulse enters the cluster sums with dE = 0.  Stretched pulses count A x s (the area scales with the stretch).  dE_hit = g (A_alt - A_def) for
matched-filter hits, g A_alt - E_hit for WaveAnalyzer hits; dE (channel) = the sum over the module's hits in the cluster; E_cluster = P7's E_g.  Alternative pulses
without a default partner (new pulses) are not added (counted).  First-order propagation: clustering, seed choice and position are not redone.

Outputs (--out): exotic_P8_energy_error_20261118.pdf, .tsv (long table: chan = per-channel dE stats, cluster = per-cluster dE_cluster / E_cluster stats, scale =
the energy-scale shift sum dE / sum E_g over the calibration samples, module = per-module numbers, diag = refit diagnostics), _examples.tsv.  Pass-1 records:
--cache (default output/pileup_rates/exotic_p8_cache.npz; delete to rebuild; needed by P9).
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import time
import argparse
import multiprocessing as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import exotic_p1_quantities as p1                       # noqa: E402
import exotic_p4_module_rates as p4                     # noqa: E402
import exotic_p7_cluster_level as p7                    # noqa: E402
from exotic_p1_quantities import RT_LABEL, write_tsv    # noqa: E402
from exotic_p3_sample import text_page                  # noqa: E402

STEM = "exotic_P8_energy_error_20261118"
CATS = p7.CATS
CNAME = p7.CNAME
VARS = {1: ("bs",), 3: ("st",), 4: ("cs",), 2: ("k1", "st"), 5: ("k1", "int")}
VLABEL = {"bs": "sloped baseline (BASE_SLOPE)", "st": "stretched single (stretch_fit)", "cs": "CLIP_SINGLE 8 ns, main kept",
          "k1": "single template, k = 1", "int": "measured integral (lower bound)"}
PRIMARY = {1: "bs", 3: "st", 4: "cs", 2: "k1", 5: "k1"}  # the variant used in the headline and the combined correction (judgement: see JUDGE)
ALLV = ("bs", "st", "cs", "k1", "int")
MATCH_NS = 8.0                                          # alternative <-> default pulse match (MIN_SEP_NS)
HIT_NS = 4.0                                            # benchmark hit <-> default pulse match (h_t vs onset + t50; one sample)
N_REF = 200                                             # unclipped k = 1 reference members per split (c_rs)
REF_AMIN = 100.0                                        # ADC, reference members' minimum amplitude
CHUNK = 1500
SCAN_ROOT = p7.SCAN_ROOT
POPS = ("all", "elastic", "moller")
SUBS = (("all", None), ("hit=main", 1.0), ("hit!=main", 0.0))   # P10: the cluster hit's pulse is the record's main (largest A) pulse, or another pulse
POPL = {"all": "all clusters", "elastic": "ep elastic window", "moller": "Moller pairs"}
TA_MODS = ("W564", "W559", "W493", "W663", "W593")
H2_MODS = ("W335", "W613", "W307", "W457", "W398")         # W398: 3.49 GeV high-gain module where real elastic showers clip (P10)
RADBINS = p7.RADBINS
NRING = len(RADBINS) - 1
F32 = np.float32


def erb_mod():
    import energy_resolution_benchmark as erb
    return erb


# ---------------------------------------------------------------------------------------------------------------- fits
def fit_opt(ctx, Y, m, sg, br, bs=False, cs=None):
    """pileup_rates.fit of the rows Y with BASE_SLOPE (bs) / CLIP_SINGLE at margin cs (None: off), chunked; the module globals are restored (as exotic_p2_flags.fit_opt)."""
    import pileup_rates as pr
    old = (pr.BASE_SLOPE, pr.CLIP_SINGLE, pr.EX_CLIP_NS)
    out = []
    try:
        pr.BASE_SLOPE, pr.CLIP_SINGLE = bool(bs), cs is not None
        if cs is not None:
            pr.EX_CLIP_NS = float(cs)
        for c0 in range(0, Y.shape[0], CHUNK):
            sl = slice(c0, c0 + CHUNK)
            f = pr.fit(ctx.mf, Y[sl], m[sl], sg[sl], br[sl])
            f["gate"] = np.asarray(ctx.mf.gate_rows).copy() if bs else np.zeros(Y[sl].shape[0], bool)
            out.append(f)
    finally:
        pr.BASE_SLOPE, pr.CLIP_SINGLE, pr.EX_CLIP_NS = old
    return {q: np.concatenate([f[q] for f in out]) for q in out[0]} if out else {}


def sub(f, rows):
    return {q: v[rows] for q, v in f.items()}


def start_baseline(Y, sg, br):
    """Start baseline of run_ref: the 30-sample pedestal, b_ref where that is farther than START_NSIGMA sigma_ref."""
    import pileup_rates as pr
    ped = pr.pedestal(Y.astype(F32))
    far = np.abs(ped - br) > pr.START_NSIGMA * sg
    return np.where(far, br, ped).astype(F32)


def group_fit(ctx, Y, m, sg, br, t_init, v_init):
    """MF.fit_group on Y (overflow samples masked) from onsets t_init [n, K] (valid v_init).  Returns t [n, K], a [n, K], valid, model [n, ns] (raw ADC)."""
    import pileup_rates as pr
    n, ns = Y.shape
    start = start_baseline(Y, sg, br)
    W = (Y < pr.OVERFLOW).astype(F32)
    y = (Y.astype(F32) - start[:, None])
    t, a, v, model = np.zeros(t_init.shape), np.zeros(t_init.shape), np.zeros(t_init.shape, bool), np.full((n, ns), np.nan)
    clip = (W < 1).any(1)
    for msk in (False, True):
        r = np.flatnonzero(clip == msk)
        if not r.size:
            continue
        w = W[r] if msk else None
        yy = y[r] * W[r] if msk else y[r]
        t0, a0, v0, mod, db = ctx.mf.fit_group(yy, w, m[r], sg[r].astype(F32), t_init[r].astype(F32), v_init[r], ns)
        t[r], a[r], v[r] = t0, np.where(v0, a0, 0.0), v0
        S = sum(np.where(v0[:, k, None], a0[:, k, None], 0.0) * ctx.mf.basis(m[r], t0[:, k].astype(F32), ns) for k in range(t_init.shape[1]))
        model[r] = start[r, None] + pr.base_curve(db, ns) + S
    return t, a, v, model


def main_of(f):
    """Default main pulse index (largest A) and its onset."""
    j = np.argmax(np.where(np.arange(f["a"].shape[1])[None, :] < f["k"][:, None], f["a"], -1.0), axis=1)
    return j, f["t"][np.arange(j.size), j]


def single(t, a, n, K):
    T, A, V = np.full((n, K), np.nan), np.zeros((n, K)), np.zeros((n, K), bool)
    T[:, 0], A[:, 0], V[:, 0] = t, np.where(np.isfinite(a), a, 0.0), np.isfinite(a) & np.isfinite(t)
    return T, A, V


def alt_fits(ctx, Y, mi, sg, br, cat, fd, crs=None):
    """Alternative fits of rows Y (module mi, category cat [n]) given the default fit fd.  Returns {variant: dict(t, A (area-equivalent amplitude), v,
    model [n, ns] raw ADC, rows (indices into Y), extras)}; crs: c_rs of the f5 integral (None: 1)."""
    import pileup_rates as pr
    n, ns = Y.shape
    K = pr.KMAX
    m = ctx.bank.row[mi]
    out = {}
    jm, tm = main_of(fd)
    # f1: sloped baseline
    r = np.flatnonzero(cat == 1)
    if r.size:
        f = fit_opt(ctx, Y[r], m[r], sg[r], br[r], bs=True)
        out["bs"] = dict(rows=r, t=f["t"], A=f["a"], v=np.arange(K)[None, :] < f["k"][:, None], k=f["k"], chi2=f["chi2"], gate=f["gate"], bslope=f["bslope"],
                         model=pr.fit_baseline(f, ns) + f["comp"].sum(1))
    # f4: CLIP_SINGLE, main kept
    r = np.flatnonzero(cat == 4)
    if r.size:
        f = fit_opt(ctx, Y[r], m[r], sg[r], br[r], cs=pr.EX_CLIP_NS)
        T, A, V = f["t"].copy(), f["a"].copy(), np.arange(K)[None, :] < f["k"][:, None]
        model = f["ped"][:, None] + f["comp"].sum(1)
        has = (V & (np.abs(np.nan_to_num(T, nan=-1e9) - tm[r, None]) <= MATCH_NS)).any(1)
        fb = np.flatnonzero(~has)
        if fb.size:
            rr = r[fb]
            Yr = Y[rr].astype(np.float64)
            ov = Yr >= pr.OVERFLOW
            first, last = ov.argmax(1), ns - 1 - ov[:, ::-1].argmax(1)
            lo, hi = pr.DT * first - pr.EX_CLIP_NS, pr.DT * last + pr.EX_CLIP_NS
            td, vd = fd["t"][rr], np.arange(K)[None, :] < fd["k"][rr][:, None]
            keep = vd & ((np.arange(K)[None, :] == jm[rr][:, None]) | (td < lo[:, None]) | (td > hi[:, None]))
            g_ = keep.any(1)                            # rows with no default pulse keep the CLIP_SINGLE result
            fb, rr, keep, td = fb[g_], rr[g_], keep[g_], td[g_]
        if fb.size:
            t1, a1, v1, mod1 = group_fit(ctx, Y[rr], m[rr], sg[rr], br[rr], np.nan_to_num(td, nan=0.0), keep)
            o = np.argsort(np.where(v1, t1, np.inf), axis=1)
            T[fb] = np.where(np.take_along_axis(v1, o, 1), np.take_along_axis(t1, o, 1), np.nan)
            A[fb] = np.take_along_axis(a1, o, 1)
            V[fb] = np.take_along_axis(v1, o, 1)
            model[fb] = mod1
        out["cs"] = dict(rows=r, t=T, A=A, v=V, k=V.sum(1), chi2=f["chi2"], fallback=~has, k_cs=f["k"], model=model)
    # f2 / f3: stretched single
    r = np.flatnonzero((cat == 3) | (cat == 2))
    if r.size:
        sf = pr.stretch_fit(ctx.mf, Y[r], m[r], sub(fd, r), sg[r])
        T, A, V = single(sf["t0"], sf["A"] * sf["s"], r.size, K)
        out["st"] = dict(rows=r, t=T, A=A, v=V, s=sf["s"], A_raw=sf["A"], chi2s=sf["chi2s"], chi2k=sf["chi2k"], nm=sf["nm"], model=sf["model"])
    # f2 / f5: one template pulse at the default main onset (overflow masked)
    r = np.flatnonzero((cat == 2) | (cat == 5))
    if r.size:
        ti = np.zeros((r.size, 1))
        ti[:, 0] = np.nan_to_num(tm[r], nan=0.0)
        t1, a1, v1, mod1 = group_fit(ctx, Y[r], m[r], sg[r], br[r], ti, np.isfinite(tm[r])[:, None])
        T, A, V = single(np.where(v1[:, 0], t1[:, 0], np.nan), np.where(v1[:, 0], a1[:, 0], np.nan), r.size, K)
        out["k1"] = dict(rows=r, t=T, A=A, v=V, model=mod1)
        # f5: measured integral of the clipped waveform (window of the k1 model)
        r5 = np.flatnonzero(cat[r] == 5)
        if r5.size:
            rr = r[r5]
            shp = np.where(V[r5, :1], A[r5, :1], 0.0) * ctx.mf.basis(m[rr], np.nan_to_num(T[r5, 0], nan=0.0).astype(F32), ns)
            y = Y[rr].astype(np.float64) - fd["ped"][rr].astype(np.float64)[:, None]
            raw1 = pr.raw_sums(y, shp.astype(np.float64), np.nan_to_num(T[r5, 0], nan=0.0))[0]
            ai = np.where(V[r5, 0], raw1 / ((1.0 if crs is None else crs) * ctx.I[mi[rr]]), np.nan)
            Ti, Ai, Vi = single(T[r5, 0], ai, rr.size, K)
            out["int"] = dict(rows=rr, t=Ti, A=Ai, v=Vi, raw1=raw1, model=None)
    return out


def assign(hits, t50, tdef, adef, vdef, alt, multi, s_eq=1.0, mod_hits=(), jmain=-1):
    """dE of one channel row: hits [(h_t raw, E_hit, h_a, h_k, g, hit row)] of this cluster record, default pulses (tdef, adef, vdef), alternative (t, A, v);
    mod_hits [(hit row, h_t raw)] = every in-window hit of the module in the event (any cluster), jmain = the default main (largest A) pulse.
    multi (bs, cs): each default pulse gets the alternative pulse matched one-to-one within MATCH_NS (0 if none).  Single-pulse variants (st, k1, int): the
    single pulse is the module's hit when its time onset + s_eq t50 (s_eq = the stretch) is inside the benchmark window (PRE_LO, PRE_HI); it replaces the
    default pulse of the module's hit nearest in time among ALL its in-window hits of the event (so a module split over two clusters is credited once: only the
    record holding that hit gets the pulse), every other hit's pulse gets 0 (outside the window: every hit gets 0).
    Returns (dE, n_unmatched_hits, n_new_pulses, hit_main); dE NaN when a single-pulse alternative has no pulse; hit_main = a hit of this record is the default
    main pulse."""
    erb = erb_mod()
    jd = np.flatnonzero(vdef)
    td = tdef[jd] + t50
    ta, aa, va = alt
    ja = np.flatnonzero(va)
    if not multi and not ja.size:
        return np.nan, 0, 0, None                       # hit_main left to another variant of the row
    # default pulse of each hit
    hi_, nun = [], 0
    for ht, e, ha, hk, g, _ in hits:
        if not jd.size:
            hi_.append(-1)
            nun += 1
            continue
        d = np.abs(td - ht)
        i = int(jd[np.argmin(d)])
        if d.min() > HIT_NS:
            nun += 1
            if hk != 0 and ha > 0:                      # nearest by amplitude instead
                i = int(jd[np.argmin(np.abs(adef[jd] - ha))])
        hi_.append(i)
    anew = np.zeros(adef.size)
    used = set()
    if multi:
        pairs = sorted((abs(tdef[i] - ta[j]), i, j) for i in jd for j in ja if abs(tdef[i] - ta[j]) <= MATCH_NS)
        di = set()
        for _, i, j in pairs:
            if i in di or j in used:
                continue
            di.add(i)
            used.add(j)
            anew[i] = aa[j]
        nnew = len(set(ja.tolist()) - used)
    else:
        nnew = 0
        teq = ta[ja[0]] + s_eq * t50
        cand = list(mod_hits) or [(h[5], h[0]) for h in hits]
        if erb.PRE_LO < teq < erb.PRE_HI and cand:
            row = min(cand, key=lambda c: abs(c[1] - teq))[0]
            for k, h in enumerate(hits):
                if h[5] == row and hi_[k] >= 0:
                    anew[hi_[k]] = aa[ja[0]]
    dE = 0.0
    for (ht, e, ha, hk, g, _), i in zip(hits, hi_):
        a0, a1 = (adef[i], anew[i]) if i >= 0 else (0.0, 0.0)
        dE += g * (a1 - a0) if (hk != 0 and ha > 0 and i >= 0) else g * a1 - e
    return dE, nun, nnew, bool(jmain >= 0 and jmain in hi_)


# ---------------------------------------------------------------------------------------------------------------- pass 1
def split_hits_a(b, fn_b, split, evsel, geo, rt):
    """p7.split_hits with h_a (benchmark matched-filter amplitude) carried along."""
    erb = erb_mod()
    h_ev = b["h_ev"]
    lo, hi = np.searchsorted(h_ev, [evsel[0], evsel[-1] + 1]) if evsel.size else (0, 0)
    h = {q: b[q][lo:hi] for q in ("h_mod", "h_t", "h_e", "h_k", "h_extra", "h_a")} | {"h_ev": h_ev[lo:hi]}
    keep = np.isin(h["h_ev"], evsel)
    h = {q: v[keep] for q, v in h.items()}
    if erb.PHASE_RESID_BY_RT.get(rt):
        erb.PHASE_RESID = str(erb.TABLE_DIR / erb.PHASE_RESID_BY_RT[rt])
        h = erb.phase_resid(h, b, geo, split, str(fn_b))
    return h


def flagged_rows(z):
    """npz row indices, keys (event number * 4096 + module) and exclusive categories of the flagged PbWO4 rows (p7.row_categories' rule); sorted keys of every row."""
    midx = z["ch_midx"].astype(np.int64)
    sel = z["mod_pwo4"][midx]
    k = z["ch_k"].astype(np.int64)
    a, t = z["ch_a"], z["ch_t"]
    main = np.argmax(np.where(np.isfinite(t), a, -1.0), axis=1)
    t0 = t[np.arange(k.size), main].astype(np.float64)
    with np.errstate(invalid="ignore"):
        trunc = (np.ceil(t0 / p4.DT) - p4.LEAD + p4.LW > p4.NS) & np.isfinite(t0)
    c = dict(db=z["ch_ex_db"].astype(np.float64), tail=z["ch_ex_tail"].astype(np.float64), w50=z["ch_w50"].astype(np.float64), k=k,
             nm=z["ch_ex_nm"].astype(np.int64), c2s=z["ch_ex_chi2s"].astype(np.float64), c2k=z["ch_ex_chi2k"].astype(np.float64),
             nclip=z["ch_nclip"].astype(np.int64), trunc=trunc, s4=z["ch_ex_f4"].astype(bool))
    F = p4.flags_of(c, p4.UNI)
    cat = np.where(sel, F["cat"], 0).astype(np.int8)
    key = z["ev_number"][z["ch_ev"].astype(np.int64)].astype(np.int64) * 4096 + midx
    o = np.argsort(key, kind="stable")
    return key, cat, o


def read_rows(ctx, f, want):
    """Decode-only EVIO re-read of the events in want {event number: [module]}: {(event number, module): samples}."""
    from prad2py import dec
    out, left = {}, set(want)
    ch = ctx.open(f)
    while left and ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            ch.select_event(i)
            en = int(ch.info().event_number)
            if en not in left:
                continue
            left.discard(en)
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            need = set(want[en])
            for mi, s in ctx.channels(d["event"]):
                if mi in need:
                    out[(en, mi)] = s
    ch.close()
    return out


CREC = ("c_run", "c_split", "c_evn", "c_seed", "c_eg", "c_cnt")
RREC = ("evn", "midx", "cat", "role", "ci", "e_def", "nh", "nwa", "k_def", "chi2_def", "repro", "a_main", "t_main", "nclip", "bslope", "s", "k_alt", "k_cs",
        "fallback", "gate", "raw1", "iint", "g", "nun", "nnew", "sg", "br", "hit_main") + tuple(f"de_{v}" for v in ALLV) + tuple(f"am_{v}" for v in ALLV)


def work(job):
    """One exotic-scan split: recluster the events with a flagged row (P7 pass 1), refit the flagged members (default and alternatives), dE per member."""
    src, fn, split, evio_dir = job
    import pileup_rates as pr
    erb = erb_mod()
    t0 = time.monotonic()
    z = np.load(fn)
    run, rt = int(z["run"]), str(z["run_type"])
    key, cat, korder = flagged_rows(z)
    skey = key[korder]
    fl = np.flatnonzero(cat > 0)
    fn_b = p7.bench_file(rt, run)
    b = np.load(fn_b)
    bn = b["ev_number"]
    evsel = np.flatnonzero((b["ev_split"] == split) & np.isin(bn, z["ev_number"]))
    gate = erb.CONTROL_GATE[rt]
    cls, geo = erb.clusterers(run, str(b["calib"]), (gate,))
    gains = p7.run_gains(rt, run, geo["name"])
    t50 = erb.template_columns(p1.PYDIR / str(b["bank"]), geo["name"])[0]
    h = split_hits_a(b, fn_b, split, evsel, geo, rt)
    # module median g (MeV per benchmark amplitude) over the split: WaveAnalyzer hits
    mf_ = (h["h_k"] != 0) & (h["h_a"] > 0)
    gm = gains[h["h_mod"][mf_]] * h["h_e"][mf_].astype(np.float64) / h["h_a"][mf_]
    nmod = geo["x"].size
    gmod = np.full(nmod, np.nan)
    for m_ in np.unique(h["h_mod"][mf_]).tolist():
        gmod[m_] = np.median(gm[h["h_mod"][mf_] == m_])
    gmod = np.where(np.isfinite(gmod), gmod, np.nanmedian(gmod) if np.isfinite(gmod).any() else 1.0)
    # restrict to events with a flagged row
    fev = set((key[fl] // 4096).tolist())
    hev = np.isin(bn[h["h_ev"]], list(fev))
    h = {q: v[hev] for q, v in h.items()}
    ev = h["h_ev"]
    lo = np.flatnonzero(np.r_[True, ev[1:] != ev[:-1]]) if ev.size else np.zeros(0, int)
    hi = np.r_[lo[1:], ev.size]
    toff = geo["toff"][h["h_mod"]] if geo.get("toff") is not None else 0.0
    tt = (h["h_t"].astype(np.float64) - toff).tolist()
    mod, e = h["h_mod"].tolist(), h["h_e"].tolist()
    cl = cls[f"g{gate:g}"]
    pw = geo["pwo4"]
    C = {q: [] for q in CREC}
    R = []                                              # (evn, midx, npz row, cat, role, ci, e_def, hits)
    REF = []                                            # (evn, midx, npz row)
    rng = np.random.default_rng(run * 1000 + split)
    zk, zclip, za0 = z["ch_k"], z["ch_nclip"], z["ch_a"][:, 0]           # loaded once (lazy npz access decompresses per call)
    for s, f_ in zip(lo.tolist(), hi.tolist()):
        cl.clear()
        hrow = {}
        for i in range(s, f_):
            cl.add_hit(mod[i], e[i], tt[i])
            hrow.setdefault(mod[i], []).append(i)
        cl.form_clusters()
        en_ = int(bn[ev[s]])
        for r in cl.reconstruct_matched():
            seed = r.cluster.center.index
            if not pw[seed]:
                continue
            em, hl = {}, {}
            for x in r.cluster.hits:                   # hit row: same module, nearest input energy (the clusterer's non-linearity rescales x.energy)
                cand = hrow.get(x.index, [])
                i = min(cand, key=lambda i: abs(e[i] - x.energy)) if cand else None
                em[x.index] = em.get(x.index, 0.0) + x.energy * gains[x.index]
                hl.setdefault(x.index, []).append((i, x.energy))
            mi = list(em)
            me = np.fromiter(em.values(), float, len(em))
            eg = float(me.sum())
            smax = mi[int(me.argmax())]
            frow = []
            for m_ in mi:
                kk = en_ * 4096 + m_
                p = np.searchsorted(skey, kk)
                row = int(korder[p]) if p < skey.size and skey[p] == kk else -1
                frow.append(row)
            if not any(row >= 0 and cat[row] > 0 for row in frow):
                continue
            ci = len(C["c_run"])
            for q, v in zip(CREC, (run, split, en_, seed, np.float32(eg), len(mi))):
                C[q].append(v)
            for m_, row in zip(mi, frow):
                hits = []
                for i, xe in hl[m_]:
                    if i is None:
                        continue
                    ha, hk = float(h["h_a"][i]), int(h["h_k"][i])
                    eh = xe * gains[m_]                 # member hit energy as P7's E_g counts it (gain, phase factor, clusterer non-linearity)
                    g = eh / ha if (hk != 0 and ha > 0) else gmod[m_] * xe / e[i]
                    hits.append((float(h["h_t"][i]), eh, ha, hk, g, i))
                if row >= 0 and cat[row] > 0:
                    R.append((en_, m_, row, int(cat[row]), 0 if m_ == smax else 1, ci, em[m_], hits,
                              [(j_, float(h["h_t"][j_])) for j_ in hrow.get(m_, [])]))
                elif row >= 0 and pw[m_] and zk[row] == 1 and zclip[row] == 0 and za0[row] >= REF_AMIN:
                    REF.append((en_, m_, row))
    REF = [REF[i] for i in sorted(rng.choice(len(REF), min(N_REF, len(REF)), replace=False))] if REF else []
    # EVIO re-read
    ctx = p1.get_ctx(run)
    want = {}
    for x in R:
        want.setdefault(x[0], []).append(x[1])
    for x in REF:
        want.setdefault(x[0], []).append(x[1])
    fe = Path(evio_dir) / f"prad_{run:06d}.evio.{split:05d}"
    samp = read_rows(ctx, str(fe), want) if want else {}
    j = list(z["splits"]).index(split)
    b_ref, sig_ref = z["mod_b_ref"][j].astype(np.float64), z["mod_sig_ref"][j].astype(np.float64)
    # reference members: raw sum / (A I_m) of unclipped k = 1 members (c_rs)
    ref = [x for x in REF if (x[0], x[1]) in samp]
    rs_ratio = np.zeros(0)
    if ref:
        Y = np.stack([samp[(x[0], x[1])] for x in ref])
        mi = np.array([x[1] for x in ref])
        fdr = fit_opt(ctx, Y, ctx.bank.row[mi], sig_ref[mi], b_ref[mi])
        okr = fdr["k"] == 1
        y = Y.astype(np.float64) - fdr["ped"].astype(np.float64)[:, None]
        jm, tm = main_of(fdr)
        comp = fdr["comp"][np.arange(mi.size), jm].astype(np.float64)
        raw1 = pr.raw_sums(y, comp, np.nan_to_num(tm, nan=0.0))[0]
        rs_ratio = (raw1 / (fdr["a"][np.arange(mi.size), jm] * ctx.I[mi]))[okr]
    c_rs = float(np.median(rs_ratio)) if rs_ratio.size >= 20 else np.nan
    ok = [x for x in R if (x[0], x[1]) in samp]
    nmiss = len(R) - len(ok)
    out = {q: [] for q in RREC}
    Ys = np.zeros((0, pr.NS_EXPECT), np.uint16)
    if ok:
        Y = np.stack([samp[(x[0], x[1])] for x in ok])
        Ys = Y
        mi = np.array([x[1] for x in ok])
        rows = np.array([x[2] for x in ok])
        catv = np.array([x[3] for x in ok])
        sg, br = sig_ref[mi], b_ref[mi]
        m = ctx.bank.row[mi]
        fd = fit_opt(ctx, Y, m, sg, br)
        repro = (fd["k"] == z["ch_k"][rows]) & np.isclose(fd["chi2"], z["ch_chi2"][rows], rtol=2e-3, atol=1e-3)
        alt = alt_fits(ctx, Y, mi, sg, br, catv, fd, c_rs if np.isfinite(c_rs) else 1.0)
        jm, tm = main_of(fd)
        n = len(ok)
        per = {v: (np.full(n, -1), d) for v, d in alt.items()}
        for v, (pos, d) in per.items():
            pos[d["rows"]] = np.arange(d["rows"].size)
        K = pr.KMAX
        vdef = np.arange(K)[None, :] < fd["k"][:, None]
        for q in RREC:
            out[q] = [np.nan] * n if q not in ("evn", "midx", "cat", "role", "ci") else [0] * n
        for i, x in enumerate(ok):
            de = {}
            nun = nnew = 0
            hit_main = None
            for v in ALLV:
                pos, d = per.get(v, (None, None))
                if d is None or pos[i] < 0:
                    de[v] = np.nan
                    continue
                p = pos[i]
                a_ = (d["t"][p], d["A"][p], d["v"][p])
                de[v], nu, nn, hm = assign(x[7], t50[x[1]], fd["t"][i], fd["a"][i], vdef[i], a_, v in ("bs", "cs"),
                                           float(d["s"][p]) if v == "st" and np.isfinite(d["s"][p]) else 1.0, x[8], int(jm[i]) if fd["k"][i] else -1)
                nun, nnew = max(nun, nu), max(nnew, nn)
                if hm is not None and (v in ("bs", "cs") or hit_main is None):
                    hit_main = hm
            for q, val in (("evn", x[0]), ("midx", x[1]), ("cat", x[3]), ("role", x[4]), ("ci", x[5]), ("e_def", x[6]), ("nh", len(x[7])),
                           ("nwa", sum(hh[3] == 0 for hh in x[7])), ("k_def", fd["k"][i]), ("chi2_def", fd["chi2"][i]), ("repro", repro[i]),
                           ("a_main", fd["a"][i, jm[i]]), ("t_main", tm[i]), ("nclip", int((Y[i] >= pr.OVERFLOW).sum())), ("iint", ctx.I[x[1]]),
                           ("g", np.median([hh[4] for hh in x[7]]) if x[7] else np.nan), ("nun", nun), ("nnew", nnew), ("sg", sg[i]), ("br", br[i]),
                           ("hit_main", np.nan if hit_main is None else float(hit_main))):
                out[q][i] = val
            for v in ALLV:
                out[f"de_{v}"][i] = de[v]
                pos, d = per.get(v, (None, None))
                out[f"am_{v}"][i] = (d["A"][pos[i]].max() if d["v"][pos[i]].any() else 0.0) if (d is not None and pos[i] >= 0) else np.nan
            if "bs" in per and per["bs"][0][i] >= 0:
                p = per["bs"][0][i]
                out["bslope"][i], out["k_alt"][i], out["gate"][i] = per["bs"][1]["bslope"][p], per["bs"][1]["k"][p], per["bs"][1]["gate"][p]
            if "cs" in per and per["cs"][0][i] >= 0:
                p = per["cs"][0][i]
                out["k_alt"][i], out["k_cs"][i], out["fallback"][i] = per["cs"][1]["k"][p], per["cs"][1]["k_cs"][p], per["cs"][1]["fallback"][p]
            if "st" in per and per["st"][0][i] >= 0:
                out["s"][i] = per["st"][1]["s"][per["st"][0][i]]
            if "int" in per and per["int"][0][i] >= 0:
                out["raw1"][i] = per["int"][1]["raw1"][per["int"][0][i]]
    d = {q: np.array(v) for q, v in out.items()}
    d["Y"] = Ys
    c = {q: np.array(v) for q, v in C.items()}
    meta = dict(src=src, run=run, rt=rt, split=split, fn=str(fn), n_rows=len(R), n_missing=nmiss, n_clusters=len(c["c_run"]), wall=time.monotonic() - t0,
                n_ref=int(rs_ratio.size), c_rs=c_rs)
    print(f"  [{src}] run {run} split {split}: {len(c['c_run'])} flagged clusters, {len(R)} flagged members ({nmiss} not re-read), "
          f"repro {np.mean(d['repro']) if len(R) else np.nan:.4f}, {meta['wall']:.0f} s", flush=True)
    return meta, d, c, rs_ratio


def jobs_of(evio_dir, only=""):
    jobs = [(s, f, sp, evio_dir) for s, f, sp in p7.jobs_of() if not (s == "sample" and "Ta" in Path(f).name)]
    out = []
    for s, f, sp, ed in jobs:
        with np.load(f) as z:
            if s == "sample" and "Ta" in str(z["run_type"]):
                continue
        if only and only not in Path(f).name:
            continue
        out.append((s, f, sp, ed))
    return out


def collect(jobs, nproc, cache):
    if cache and Path(cache).exists():
        c = np.load(cache, allow_pickle=True)
        print(f"[cache] {cache}", flush=True)
        R = {q[2:]: c[q] for q in c.files if q.startswith("r_")}
        C = {q[2:]: c[q] for q in c.files if q.startswith("c_") and q[2:] in CREC}
        return list(c["metas"]), R, C, {k: v for k, v in c["rs"].item().items()}
    t0 = time.monotonic()
    with mp.get_context("fork").Pool(nproc) as pool:
        res = pool.map(work, sorted(jobs, key=lambda j: -os.path.getsize(j[1])), chunksize=1)
    metas = [r[0] for r in res]
    off = np.cumsum([0] + [len(r[2]["c_run"]) for r in res])
    R = {}
    for q in RREC + ("Y",):
        parts = [r[1][q] for r in res if len(r[1]["evn"])]
        R[q] = np.concatenate(parts) if parts else np.zeros(0)
    R["ci"] = np.concatenate([np.asarray(r[1]["ci"], np.int64) + off[i] for i, r in enumerate(res) if len(r[1]["evn"])])
    R["job"] = np.concatenate([np.full(len(r[1]["evn"]), i, np.int32) for i, r in enumerate(res) if len(r[1]["evn"])])
    C = {q: np.concatenate([np.asarray(r[2][q]) for r in res]) for q in CREC}
    rs = {}
    for m, r in zip(metas, res):
        rs.setdefault(m["rt"], []).append(r[3])
    rs = {k: np.concatenate(v) for k, v in rs.items()}
    print(f"[collect] {len(jobs)} splits, {R['evn'].size} members, {time.monotonic() - t0:.0f} s", flush=True)
    for q in ("evn", "midx", "cat", "role", "nh", "nwa", "k_def", "nun", "nnew"):
        R[q] = R[q].astype(np.int64)
    R["Y"] = R["Y"].astype(np.uint16)
    if cache:
        Path(cache).parent.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(cache, metas=np.array(metas, dtype=object), rs=np.array(rs, dtype=object), **{"r_" + q: v for q, v in R.items()},
                            **{"c_" + q: v for q, v in C.items()})
    return metas, R, C, rs


# ---------------------------------------------------------------------------------------------------------------- pass 2
def qstats(x):
    x = x[np.isfinite(x)]
    if not x.size:
        return dict(n=0, med=np.nan, q16=np.nan, q84=np.nan, mean=np.nan, sem=np.nan)
    q = np.quantile(x, (0.5, 0.16, 0.84))
    return dict(n=int(x.size), med=float(q[0]), q16=float(q[1]), q84=float(q[2]), mean=float(x.mean()), sem=float(x.std(ddof=1) / np.sqrt(x.size)) if x.size > 1 else np.nan)


def join(metas, R, C, pm, D):
    """Cache index of every P8 cluster (-1 if not found) by (run, split, event number, seed, E_g); P7 rows of the matching source only (H2 sample, Ta full)."""
    J = D["job"]
    src = np.array([m["src"] for m in pm])[J]
    rt = np.array([m["rt"] for m in pm])[J]
    run = np.array([m["run"] for m in pm], np.int64)[J]
    use = ((D["scat"] > 0) | (D["nmask"] > 0)) & (((rt != "2p24GeV_Ta_X17") & (src == "sample")) | ((rt == "2p24GeV_Ta_X17") & (src == "full")))
    idx = np.flatnonzero(use)
    lut = {}
    for i, k in zip(idx.tolist(), zip(run[idx].tolist(), D["split"][idx].tolist(), D["evn"][idx].tolist(), D["seed"][idx].tolist(),
                                       D["eg"][idx].astype(np.float32).tolist())):
        lut.setdefault(k, i)
    keys = zip(C["c_run"].tolist(), C["c_split"].tolist(), C["c_evn"].tolist(), C["c_seed"].tolist(), C["c_eg"].astype(np.float32).tolist())
    return np.array([lut.get(k, -1) for k in keys], np.int64), int(use.sum())


def scale_rows(lab, rt, P, D, X, dEc, cl_idx, tag, ringsel=None):
    """Energy-scale shift S = sum dE_cluster / sum E_g over the clusters of population masks P (elastic, moller), overall and per ring; sigma = sqrt(sum dE^2) / sum E_g."""
    rows = []
    for pop in ("elastic", "moller", "all"):
        if pop not in P:
            continue
        for rb in [None] + list(range(NRING)):
            m = P[pop] if rb is None else P[pop] & (X["ring"] == rb)
            den = float(D["eg"][m].astype(np.float64).sum())
            n = int(m.sum())
            sel = m[cl_idx]
            num = float(dEc[sel].sum())
            var = float((dEc[sel] ** 2).sum())
            rows.append(dict(table="scale", dataset=lab, run_type=rt, population=pop, ring="all" if rb is None else f"{RADBINS[rb]:g}-{RADBINS[rb + 1]:g}",
                             category=tag[0], variant=tag[1], n_clusters=n, n_flag=int(sel.sum()), S=num / den if den else np.nan,
                             S_err=np.sqrt(var) / den if den else np.nan, sum_dE=num, sum_E=den))
    return rows


def report(metas, R, C, rs, out, n_ex, evio_dir, nproc):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    pm, D = p7.collect([], 1, str(SCAN_ROOT / "exotic_p7_cache.npz"))
    with np.load(next((SCAN_ROOT / "exotic_sample").glob("run_*.npz"))) as z:
        names = [str(x) for x in z["mod_name"]]
    X = p7.derive(pm, D, names)
    dss = [d for d in p7.datasets(pm) if d[0] != "Ta_sample"]
    mol = p7.moller_keys()
    PP = {}
    for ds in dss:
        P, _, _ = p7.populations(pm, D, X, ds, mol)
        PP[ds[2]] = P
    cidx, n_cache_flag = join(metas, R, C, pm, D)
    ci = cidx[R["ci"]]
    good = ci >= 0
    print(f"[join] P8 clusters {cidx.size}, matched {int((cidx >= 0).sum())}; P7 flagged clusters of the same sources {n_cache_flag}; members matched {int(good.sum())} / {good.size}", flush=True)
    rtm = np.where(good, X["rt"][np.maximum(ci, 0)], "")
    egm = np.where(good, D["eg"][np.maximum(ci, 0)].astype(np.float64), np.nan)
    ring = np.where(good, X["ring"][np.maximum(ci, 0)], -1)
    pop = {p: np.zeros(ci.size, bool) for p in POPS}
    for rt, P in PP.items():
        for p in POPS:
            if p in P:
                pop[p] |= good & P[p][np.maximum(ci, 0)] & (rtm == rt)
    rts = [ds[2] for ds in dss]
    rows = []
    # ---- per channel
    for c in CATS:
        for v in VARS[c]:
            de = R[f"de_{v}"].astype(np.float64)
            for rt in rts:
                for role, rm in (("any", None), ("seed", 0), ("neighbour", 1)):
                  for sub, sm in (SUBS if c == 1 else SUBS[:1]):
                   for p in POPS:
                    for rb in [None] + (list(range(NRING)) if role == "any" and sub == "all" else []):
                        m = good & (R["cat"] == c) & (rtm == rt) & pop[p] & (True if rm is None else R["role"] == rm) & (True if sm is None else R["hit_main"] == sm)
                        m &= True if rb is None else ring == rb
                        if sub != "all" and not m.any():
                            continue
                        s1, s2 = qstats(de[m]), qstats(de[m] / egm[m])
                        rows.append(dict(table="chan", run_type=rt, category=CNAME[c], variant=v, role=role, sub=sub, population=p,
                                         ring="all" if rb is None else f"{RADBINS[rb]:g}-{RADBINS[rb + 1]:g}", n=s1["n"],
                                         dE_med=s1["med"], dE_q16=s1["q16"], dE_q84=s1["q84"], dE_mean=s1["mean"], dE_sem=s1["sem"],
                                         f_med=s2["med"], f_q16=s2["q16"], f_q84=s2["q84"], f_mean=s2["mean"], f_sem=s2["sem"],
                                         rel_med=qstats(de[m] / R["e_def"][m])["med"] if m.any() else np.nan))
    # ---- per cluster (members of category c corrected) and scale shift
    dEc_all = {}
    for c in CATS:
        for v in VARS[c]:
            de = np.nan_to_num(R[f"de_{v}"].astype(np.float64))
            m = good & (R["cat"] == c)
            u, inv = np.unique(ci[m], return_inverse=True)
            dEc = np.bincount(inv, weights=de[m], minlength=u.size) if u.size else np.zeros(0)
            seedc = np.zeros(u.size, bool)
            if u.size:
                seedc[np.unique(inv[R["role"][m] == 0])] = True
            dEc_all[(c, v)] = (u, dEc, seedc)
            for rt in rts:
                P = PP[rt]
                inrt = X["rt"][u] == rt
                for p in POPS:
                    if p not in P:
                        continue
                    for role, rsel in (("any", np.ones(u.size, bool)), ("seed", seedc), ("neighbour", ~seedc)):
                        for rb in [None] + list(range(NRING)):
                            mm = inrt & P[p][u] & rsel & (True if rb is None else X["ring"][u] == rb)
                            s = qstats(dEc[mm] / D["eg"][u][mm])
                            rows.append(dict(table="cluster", run_type=rt, category=CNAME[c], variant=v, role=role, population=p,
                                             ring="all" if rb is None else f"{RADBINS[rb]:g}-{RADBINS[rb + 1]:g}", n=s["n"], f_med=s["med"], f_q16=s["q16"],
                                             f_q84=s["q84"], f_mean=s["mean"], f_sem=s["sem"]))
                rows += scale_rows(rt, rt, {p: P[p] & (X["rt"] == rt) for p in ("elastic", "moller") if p in P}, D, X, dEc, u, (CNAME[c], v))
            for sub, sm in (SUBS[1:] if c == 1 else ()):   # P10 (f1): hit pulse = the record's main (flagged) pulse vs another pulse of the record
                ms = m & (R["hit_main"] == sm)
                if not ms.any():
                    continue
                us, invs = np.unique(ci[ms], return_inverse=True)
                dEs = np.bincount(invs, weights=de[ms], minlength=us.size)
                for rt in rts:
                    P = PP[rt]
                    rows += scale_rows(rt, rt, {p: P[p] & (X["rt"] == rt) for p in ("elastic", "moller") if p in P}, D, X, dEs, us,
                                       (CNAME[c] + " [" + sub + "]", v))
    # combined: every category with its primary variant
    de = np.zeros(ci.size)
    for c in CATS:
        de = np.where(R["cat"] == c, np.nan_to_num(R[f"de_{PRIMARY[c]}"].astype(np.float64)), de)
    u, inv = np.unique(ci[good], return_inverse=True)
    dEc_comb = np.bincount(inv, weights=de[good], minlength=u.size)
    for rt in rts:
        P = PP[rt]
        rows += scale_rows(rt, rt, {p: P[p] & (X["rt"] == rt) for p in ("elastic", "moller") if p in P}, D, X, dEc_comb, u, ("all (primary)", "primary"))
    # ---- per module
    modrows = []
    nm = np.array(names)
    for rt in rts:
        P = PP[rt]
        mods = list(TA_MODS if "Ta" in rt else H2_MODS)
        cnt = np.bincount(R["midx"][good & (rtm == rt)], minlength=len(names))
        mods += [str(nm[i]) for i in np.argsort(-cnt)[:8] if cnt[i] > 0 and str(nm[i]) not in mods]
        for mn in mods:
            mi = names.index(mn)
            for c in (0,) + CATS:
                m = good & (rtm == rt) & (R["midx"] == mi) & ((R["cat"] > 0) if c == 0 else (R["cat"] == c))
                dv = de[m] if c == 0 else np.asarray(R[f"de_{PRIMARY[c]}"], np.float64)[m]
                s1, s2 = qstats(dv), qstats(dv / egm[m])
                el = P["elastic"] & (X["rt"] == rt) & (D["seed"] == mi)
                den = float(D["eg"][el].astype(np.float64).sum())
                # elastic clusters seeded by the module: shift from all flagged members (primary variants; c = 0) or category c only
                if c == 0:
                    uu, dd = u, dEc_comb
                else:
                    uu, dd = dEc_all[(c, PRIMARY[c])][:2]
                sel = el[uu]
                modrows.append(dict(table="module", run_type=rt, module=mn, category=CNAME[c], variant="primary" if c == 0 else PRIMARY[c], n=s1["n"],
                                    n_seed=int((m & (R["role"] == 0)).sum()), n_elastic=int((m & pop["elastic"]).sum()), dE_med=s1["med"], dE_q16=s1["q16"],
                                    dE_q84=s1["q84"], dE_mean=s1["mean"], f_mean=s2["mean"], f_med=s2["med"], n_el_seeded=int(el.sum()),
                                    S_el_seeded=float(dd[sel].sum()) / den if den else np.nan, S_el_err=float(np.sqrt((dd[sel] ** 2).sum())) / den if den else np.nan))
    rows += modrows
    # ---- diagnostics
    for rt in rts:
        crs = rs.get(rt, np.zeros(0))
        for c in CATS:
            m = good & (rtm == rt) & (R["cat"] == c)
            n = int(m.sum())
            if not n:
                continue
            d = dict(table="diag", run_type=rt, category=CNAME[c], n=n, repro=float(np.mean(R["repro"][m].astype(float))),
                     unmatched_hit=float(np.mean(R["nun"][m] > 0)), wa_hit=float(np.mean(R["nwa"][m] > 0)), multi_hit=float(np.mean(R["nh"][m] > 1)),
                     k_def_med=float(np.median(R["k_def"][m])), c_rs=float(np.median(crs)) if crs.size else np.nan, n_crs=int(crs.size))
            for v in VARS[c]:
                ok = np.isfinite(R[f"de_{v}"][m].astype(float))
                d[f"{v}_ok"] = float(ok.mean())
                r_ = R[f"am_{v}"][m].astype(float) / R["a_main"][m].astype(float)
                d[f"{v}_amain_ratio_med"] = float(np.nanmedian(r_)) if np.isfinite(r_).any() else np.nan
            if c == 1:
                d.update(gate=float(np.nanmean(R["gate"][m].astype(float))), k_alt0=float(np.mean(R["k_alt"][m] == 0)),
                         k_drop=float(np.mean(R["k_alt"][m] < R["k_def"][m])), bslope_med=float(np.nanmedian(R["bslope"][m].astype(float))))
            if c == 4:
                d.update(fallback=float(np.nanmean(R["fallback"][m].astype(float))), k_cs0=float(np.mean(R["k_cs"][m] == 0)), new_pulses=float(np.mean(R["nnew"][m] > 0)))
            if c in (2, 3):
                d.update(s_med=float(np.nanmedian(R["s"][m].astype(float))))
            if c == 5:
                fr = 1.0 - R["am_int"][m].astype(float) / R["am_k1"][m].astype(float)
                d.update(clip_area_frac_med=float(np.nanmedian(fr)) if np.isfinite(fr).any() else np.nan)
            rows.append(d)
    rows.append(dict(table="diag", run_type="join", category="", n=int(good.size), repro=float(good.mean())))
    # ---- examples
    exd = examples(R, good, rtm, pop, egm, metas, n_ex, nproc, names)
    cols = ("table", "dataset", "run_type", "population", "ring", "module", "category", "variant", "role", "sub", "n", "n_clusters", "n_flag", "n_seed", "n_elastic",
            "dE_med", "dE_q16", "dE_q84", "dE_mean", "dE_sem", "f_med", "f_q16", "f_q84", "f_mean", "f_sem", "rel_med", "S", "S_err", "sum_dE", "sum_E",
            "n_el_seeded", "S_el_seeded", "S_el_err", "repro", "unmatched_hit", "wa_hit", "multi_hit", "k_def_med", "c_rs", "n_crs", "gate", "k_alt0", "k_drop",
            "bslope_med", "fallback", "k_cs0", "new_pulses", "s_med", "clip_area_frac_med") + tuple(f"{v}_{q}" for v in ALLV for q in ("ok", "amain_ratio_med"))
    write_tsv(out / f"{STEM}.tsv", rows, cols)
    write_tsv(out / f"{STEM}_examples.tsv", [{q: e[q] for q in EXCOLS} for e in exd], EXCOLS)
    pages(out, plt, PdfPages, rows, R, good, rtm, pop, egm, ring, rts, exd, metas, rs)
    for L in (headline_lines(rows, rts), chan_lines(rows, rts)):
        print("\n".join(L))
    return rows


# ---------------------------------------------------------------------------------------------------------------- examples
EXCOLS = ("run_type", "run", "split", "evn", "module", "category", "role", "population", "variant", "e_def", "e_cluster", "dE", "dE_frac", "k_def", "k_alt",
          "am_stored", "am_refit")


def pick_ex(R, good, rtm, pop, c, n, rng):
    """n members of category c in elastic / Moller clusters, round robin over the run types (random within); filled from all clusters if needed."""
    out = []
    for pref in (pop["elastic"] | pop["moller"], good):
        pools = {rt: list(rng.permutation(np.flatnonzero(pref & (R["cat"] == c) & (rtm == rt) & np.isfinite(R[f"de_{PRIMARY[c]}"].astype(float)))))
                 for rt in p1.TYPES}
        while len(out) < n and any(pools.values()):
            for rt in p1.TYPES:
                while pools[rt] and pools[rt][-1] in out:
                    pools[rt].pop()
                if pools[rt] and len(out) < n:
                    out.append(int(pools[rt].pop()))
    return out


def example_job(a):
    """Refit of one stored member waveform (default and its category's alternatives) for a panel."""
    rec = a
    import pileup_rates as pr
    ctx = p1.get_ctx(rec["run"])
    Y = rec["Y"][None].astype(np.uint16)
    mi = np.array([rec["midx"]])
    sg, br = np.array([rec["sg"]]), np.array([rec["br"]])
    fd = fit_opt(ctx, Y, ctx.bank.row[mi], sg, br)
    alt = alt_fits(ctx, Y, mi, sg, br, np.array([rec["cat"]]), fd, rec["c_rs"])
    rec["def"] = dict(base=pr.fit_baseline(fd, Y.shape[1])[0], comp=fd["comp"][0][:fd["k"][0]], k=int(fd["k"][0]), a=fd["a"][0][:fd["k"][0]].tolist(),
                      t=fd["t"][0][:fd["k"][0]].tolist(), chi2=float(fd["chi2"][0]))
    rec["alt"] = {v: dict(model=None if d["model"] is None else d["model"][0], k=int(d["v"][0].sum()), A=d["A"][0][d["v"][0]].tolist(), t=d["t"][0][d["v"][0]].tolist(),
                          extra={q: float(np.asarray(d[q])[0]) for q in ("s", "bslope", "raw1") if q in d}) for v, d in alt.items()}
    return rec


def examples(R, good, rtm, pop, egm, metas, n_ex, nproc, names):
    if not n_ex:
        return []
    rng = np.random.default_rng(20261118)
    recs = []
    cmeta = {i: m for i, m in enumerate(metas)}
    for c in CATS:
        for i in pick_ex(R, good, rtm, pop, c, n_ex, rng):
            m = cmeta[int(R["job"][i])]
            v = PRIMARY[c]
            recs.append(dict(i=i, module=names[int(R["midx"][i])], run=m["run"], split=m["split"], run_type=m["rt"], evn=int(R["evn"][i]), midx=int(R["midx"][i]), cat=c, category=CNAME[c],
                             role="seed" if R["role"][i] == 0 else "neighbour", population="elastic" if pop["elastic"][i] else ("moller" if pop["moller"][i] else "other"),
                             variant=v, e_def=float(R["e_def"][i]), e_cluster=float(egm[i]), dE=float(R[f"de_{v}"][i]), dE_frac=float(R[f"de_{v}"][i] / egm[i]),
                             de_all={w: float(R[f"de_{w}"][i]) for w in VARS[c]}, k_def=int(R["k_def"][i]), k_alt=R["k_alt"][i], am_stored=float(R[f"am_{v}"][i]),
                             Y=R["Y"][i], sg=float(R["sg"][i]), br=float(R["br"][i]), c_rs=float(m.get("c_rs", 1.0)) if np.isfinite(m.get("c_rs", np.nan)) else 1.0))
    with mp.get_context("fork").Pool(nproc) as pool:
        exd = pool.map(example_job, recs, chunksize=1)
    for e in exd:
        al = e["alt"].get(e["variant"])
        e["am_refit"] = max(al["A"]) if al and al["A"] else 0.0
    return exd


# ---------------------------------------------------------------------------------------------------------------- text
def get(rows, **kw):
    for r in rows:
        if all(r.get(k) == v for k, v in kw.items()):
            return r
    return None


def pc(x, nd=4):
    return f"{100 * x:+.{nd}f}" if np.isfinite(x) else "-"


def sci(r):
    return f"{pc(r['S'])} +- {100 * r['S_err']:.{4}f}" if r and r["n_clusters"] else "-"


RINGL = ["all"] + [f"{RADBINS[i]:g}-{RADBINS[i + 1]:g}" for i in range(NRING)]


def headline_lines(rows, rts):
    L = ["Energy-scale shift S = sum dE_cluster / sum E_g (%) over every cluster of the calibration sample that correcting the category would produce;",
         "+- = sqrt(sum dE_cluster^2) / sum E_g (statistical, clusters independent).  Rings: seed radius (mm).  n_flag = clusters with a member of the category.", ""]
    for rt in rts:
        for p in ("elastic", "moller"):
            if not get(rows, table="scale", run_type=rt, population=p, ring="all"):
                continue
            n = get(rows, table="scale", run_type=rt, population=p, ring="all")["n_clusters"]
            L.append(f"{RT_LABEL[rt]}, {POPL[p]} (N = {n} clusters)")
            L.append(f"  {'category / variant':<44}{'n_flag':>8}" + "".join(f"{r_:>24}" for r_ in RINGL))
            for c in (4, 5, 1, "f1sub", 2, 3, "all"):
                if c == "f1sub":
                    for sub in ("hit=main", "hit!=main"):
                        cn = CNAME[1] + " [" + sub + "]"
                        rr = [get(rows, table="scale", run_type=rt, population=p, ring=r_, category=cn, variant="bs") for r_ in RINGL]
                        if rr[0] is not None:
                            L.append(f"  {'  ' + sub + ' (P10 split)':<44}{rr[0]['n_flag']:>8}" + "".join(f"{sci(r):>24}" for r in rr))
                    continue
                cn = "all (primary)" if c == "all" else CNAME[c]
                for v in (("primary",) if c == "all" else VARS[c]):
                    rr = [get(rows, table="scale", run_type=rt, population=p, ring=r_, category=cn, variant=v) for r_ in RINGL]
                    if rr[0] is None:
                        continue
                    tagv = "" if c == "all" else (" *" if v == PRIMARY[c] else "")
                    L.append(f"  {cn + ' / ' + v + tagv:<44}{rr[0]['n_flag']:>8}" + "".join(f"{sci(r):>24}" for r in rr))
            L.append("")
    L.append("* = primary variant (used in 'all (primary)').  f1 split (P10): hit=main = the clustered hit's pulse is the record's main (largest A) pulse,")
    L.append("hit!=main = the clustered hit is another pulse of the record (e.g. the flag sits on an earlier main pulse); dE is always the clustered pulse's own change.")
    L.append("Members whose single-pulse alternative has no pulse (k1_ok / st_ok < 1, Table 5) enter the cluster sums with dE = 0.")
    return L


def chan_lines(rows, rts, pops=("all", "elastic", "moller")):
    L = ["Per flagged member record (a channel that is a member of k clusters counts k times, with its share): dE = E_alt - E_default (MeV): median [16 %, 84 %],",
         "mean +- s.e.m.; dE / E_cluster (%): median, mean; dE / E_member median (%).  Per-ring rows (role any) and the hit=main / hit!=main split: TSV table chan.", ""]
    L.append(f"{'run type':<18}{'category / variant':<38}{'role':<10}{'pop.':<9}{'N':>7}  {'dE med [16, 84] MeV':<28}{'mean +- sem':<18}{'dE/E_cl med':>12}{'mean':>10}{'dE/E_m':>9}")
    for c in CATS:
        for v in VARS[c]:
            for rt in rts:
                for role, sub in [("any", "all"), ("seed", "all"), ("neighbour", "all")] + ([("any", "hit=main"), ("any", "hit!=main")] if c == 1 else []):
                    for p in pops:
                        r = get(rows, table="chan", run_type=rt, category=CNAME[c], variant=v, role=role, population=p, sub=sub, ring="all")
                        role_ = role if sub == "all" else sub
                        if r is None or not r["n"]:
                            continue
                        L.append(f"{RT_LABEL[rt]:<18}{CNAME[c] + ' / ' + v:<38}{role_:<10}{p:<9}{r['n']:>7}  "
                                 f"{r['dE_med']:>7.1f} [{r['dE_q16']:>7.1f}, {r['dE_q84']:>7.1f}]{'':<5}{r['dE_mean']:>7.2f} +- {r['dE_sem'] if np.isfinite(r['dE_sem']) else 0:<6.2f}"
                                 f"{pc(r['f_med'], 3):>12}{pc(r['f_mean'], 3):>10}{pc(r['rel_med'], 1):>9}")
        L.append("")
    return L


def cluster_lines(rows, rts):
    L = ["Clusters containing >= 1 member of the category (that category corrected, primary variant): dE_cluster / E_cluster (%): N, median [16 %, 84 %], mean;",
         "role seed = the category's member is the cluster's largest member; neighbour = only neighbours of the category.", ""]
    L.append(f"{'run type':<18}{'category':<22}{'role':<10}{'pop.':<9}{'N':>7}  {'med [16, 84] %':<34}{'mean %':>10}")
    for c in CATS:
        for rt in rts:
            for role in ("any", "seed", "neighbour"):
                for p in POPS:
                    r = get(rows, table="cluster", run_type=rt, category=CNAME[c], variant=PRIMARY[c], role=role, population=p, ring="all")
                    if r is None or not r["n"]:
                        continue
                    L.append(f"{RT_LABEL[rt]:<18}{CNAME[c]:<22}{role:<10}{p:<9}{r['n']:>7}  {pc(r['f_med'], 3)} [{pc(r['f_q16'], 3)}, {pc(r['f_q84'], 3)}]{'':<6}{pc(r['f_mean'], 3):>10}")
        L.append("")
    return L


def module_lines(rows, rts):
    L = ["Per module (flagged member channels of the module, primary variant): N (as seed, in elastic clusters), dE median [16, 84] MeV, mean dE / E_cluster;",
         "S_el(seeded) = sum dE_cluster / sum E_g over the elastic-window clusters seeded by the module (every flagged member corrected for 'any', the category's",
         "members only otherwise), +- statistical; n_el = elastic clusters seeded by the module.", ""]
    L.append(f"{'run type':<18}{'module':<8}{'category':<22}{'N':>7}{'seed':>6}{'el.':>6}  {'dE med [16, 84] MeV':<28}{'<dE/E_cl> %':>12}{'n_el':>8}  {'S_el(seeded) %':<24}")
    for rt in rts:
        for r in [x for x in rows if x["table"] == "module" and x["run_type"] == rt and x["n"] > 0]:
            cn = "any" if r["category"] == "any" else r["category"]
            L.append(f"{RT_LABEL[rt]:<18}{r['module']:<8}{cn:<22}{r['n']:>7}{r['n_seed']:>6}{r['n_elastic']:>6}  {r['dE_med']:>7.1f} [{r['dE_q16']:>7.1f}, {r['dE_q84']:>7.1f}]{'':<5}"
                     f"{pc(r['f_mean'], 3):>12}{r['n_el_seeded']:>8}  {pc(r['S_el_seeded'], 3)} +- {100 * r['S_el_err']:.3f}")
        L.append("")
    return L


def diag_lines(rows, rts, metas):
    L = ["Refit diagnostics per run type and category (fractions of the category's member channels):",
         "repro = default single-split refit reproduces the stored k and chi2/dof; unm. = >= 1 hit not within 4 ns of a default pulse (matched by amplitude);",
         "WA = a WaveAnalyzer hit (clipped, no matched-filter pulse); multi = > 1 hit of the module in the cluster; <v>_ok = alternative available;",
         "A ratio = largest alternative (area-equivalent) amplitude / default A_main, median.", ""]
    for rt in rts:
        for c in CATS:
            r = get(rows, table="diag", run_type=rt, category=CNAME[c])
            if r is None:
                continue
            s = (f"{RT_LABEL[rt]:<18}{CNAME[c]:<22}N {r['n']:>7}  repro {r['repro']:.4f}  unm. {r['unmatched_hit']:.4f}  WA {r['wa_hit']:.3f}  multi {r['multi_hit']:.3f}  "
                 f"k_def med {r['k_def_med']:.0f}  " + "  ".join(f"{v}_ok {r[f'{v}_ok']:.3f} A ratio {r[f'{v}_amain_ratio_med']:.3f}" for v in VARS[c]))
            if c == 1:
                s += f"  gated {r['gate']:.3f}  k_bs<k_def {r['k_drop']:.3f}  k_bs=0 {r['k_alt0']:.4f}  slope med {r['bslope_med']:.1f} ADC/100 ns"
            if c == 4:
                s += f"  CLIP_SINGLE lost main {r['fallback']:.3f} (k_cs=0 {r['k_cs0']:.3f})  new pulses {r['new_pulses']:.3f}"
            if c in (2, 3):
                s += f"  stretch s med {r['s_med']:.2f}"
            if c == 5:
                s += f"  clipped-area fraction 1 - E_int/E_k1 med {r['clip_area_frac_med']:.3f}"
            L.append(s)
        crs = [x for x in rows if x["table"] == "diag" and x["run_type"] == rt and np.isfinite(x.get("c_rs", np.nan))]
        if crs:
            L.append(f"{'':<18}c_rs (raw sum / (A I_m), unclipped k = 1 reference members) median {crs[0]['c_rs']:.4f} (n {crs[0]['n_crs']})")
        L.append("")
    j = get(rows, table="diag", run_type="join")
    L.append(f"Join: P8 member channels whose cluster is found in P7's cache: {100 * j['repro']:.3f} % of {j['n']}.")
    L.append("Members not re-read from EVIO per split: " + ", ".join(f"{m['run']}/{m['split']}: {m['n_missing']}" for m in metas if m["n_missing"]) if any(m["n_missing"] for m in metas)
             else "Every flagged member waveform was re-read from EVIO.")
    return L


JUDGE = [          # judgement from the 2026-11-18 run (numbers: Tables 1, 2, 5; examples)
    "Better estimate per category (one line each) and credibility of the alternative:",
    "",
    "f1 tail-baseline: 'bs' (sloped baseline) is the better estimate - it fits the falling baseline instead of small tail pulses and keeps the clustered pulse (Ta: k drops",
    "   in 60 % of the rows, A_main x 0.94 median, residuals no worse in the examples); credible.  dE = -8.6 % of the member energy (Ta, median), -1 to -3 % (H2).",
    "   P10 split: the clustered hit is the record's main pulse in 30 % of the Ta rows (dE -4.7 MeV median) and another pulse (flag on an earlier pulse) in 70 %",
    "   (dE -0.4 MeV median); dE is always the clustered pulse's own change; Ta elastic S -0.0076 % (hit = main) and -0.0068 % (hit != main).",
    "f3 smooth-wide: 'st' is credible only where the stretch is ~1 (Ta: s median 1.0, dE -0.3 % of E_cluster (median) in the elastic window, -0.003 % on the scale); in H2",
    "   (s median 1.4-1.9) the stretched single swallows real close pile-up (P2: type 3 is shape-indistinguishable) and adds +60-70 % to the member - NOT credible",
    "   there; the default is the better estimate for H2.",
    "f4 clipped-multi: 'cs' leaves the main amplitude unchanged in the median (A ratio 1.00, dE 0); where it rearranges a genuine pile-up train (CLIP_SINGLE lost the",
    "   main in 4-8 % of H2 rows; e.g. W299 3.49 GeV, +704 MeV with +2000 ADC residuals) it is worse than the default - the default is the better estimate, cs an upper bound.",
    "f2 wide: neither 'k1' nor 'st' is credible - H2 f2 are merged k 4-6 trains and Ta f2 baseline wobble with a real small in-window pulse; the single pulse",
    "   usually lands on another hit of the module (or outside the window), so the flagged member loses its energy: median dE -100 % of the member (3.49 GeV st -38 %),",
    "   while for the record that receives it the single pulse sums the whole train.",
    "   The default (or a veto) is the better estimate.",
    "f5 saturated single: 'k1' (extrapolation) is the better of the two for mild clipping (clipped-area fraction 2-4 % median, k1 and the integral agree to a few %;",
    "   W398 at 3.49 GeV, real elastic showers with 1 overflow sample, k = 1: dE_k1 = 0, dE_int = -23 to -34 MeV = -0.7 to -1.0 % of E_cluster, the integral's missing",
    "   clipped area); for deep saturation (tens of overflow samples, default A >> 4095 ADC, e.g. a 20.9 GeV default at 2.24 GeV) neither the extrapolation nor the",
    "   integral (a lower bound) is credible: the template does not describe the flat top.  These are P7's single-module fakes - veto rather than correct.",
    "",
    "Headline (Table 1): only Ta X17 moves measurably - elastic f1 -0.0144 % (inner ring -0.0153 %), f3 -0.0028 %, all credible -0.017 %; Moller f1 -0.0051 %",
    "(inner ring -0.0114 %).  Every credible H2 cell is below 0.01 % (largest: 3.49 GeV f5 k1 ring 320-450, +0.008 +- 0.008 %, one cluster).  No cell is above 0.1 %",
    "except the non-credible 2.24 GeV H2 f3 outer ring (450-700 mm, one cluster, +0.40 +- 0.42 %).  'all (primary)' adds every category's primary variant,",
    "including the non-credible H2 f3 / f2 / f5 rows: use the per-category rows with this page.",
    "",
    "Caveats: first-order propagation (clustering, seed, position and elastic / Moller membership not redone; a large dE can move a cluster out of the window);",
    "new pulses of the multi-pulse alternatives are not added (Table 5); a single pulse credited to a module hit outside every saved cluster is lost (one-sided:",
    "Ta elastic f3 st 196 of 1784 records, -0.0006 % of S; Ta elastic f2 k1 61 of 81 records, -0.0004 %; H2 <= 0.00001 %); H2 from the 2-file sample (f4 / f5 / f2",
    "counts small); Ta sample not refitted (subset of the full set).",
]


# ---------------------------------------------------------------------------------------------------------------- pages
CCOL = p7.CCOL
RTCOL = dict(zip(p1.TYPES, ("tab:blue", "tab:green", "tab:purple", "tab:orange")))


def pages(out, plt, PdfPages, rows, R, good, rtm, pop, egm, ring, rts, exd, metas, rs):
    head = [l for l in __doc__.split("\n")[6:33]]
    with PdfPages(out / f"{STEM}.pdf") as pdf:
        text_page(pdf, plt, "P8 energy error per exotic category: method", head, size=7.6)
        text_page(pdf, plt, "Table 1 (headline): energy-scale shift of the calibration samples per category, run type and ring", headline_lines(rows, rts), size=6.0)
        if JUDGE:
            text_page(pdf, plt, "Which alternative is the better estimate (per category)", JUDGE, size=8.4)
        cl = chan_lines(rows, rts)
        for k in range(0, len(cl), 70):
            text_page(pdf, plt, "Table 2: per-channel dE per category, run type, role, population" + (" (cont.)" if k else ""), cl[:3] + cl[k + 3 if k == 0 else k:k + 70], size=6.2)
        text_page(pdf, plt, "Table 3: dE_cluster / E_cluster of the clusters containing the category", cluster_lines(rows, rts), size=6.4)
        ml = module_lines(rows, rts)
        for k in range(0, len(ml), 80):
            text_page(pdf, plt, "Table 4: per module (Ta inner ring, H2 hot modules, top flagged modules)" + (" (cont.)" if k else ""), ml[k:k + 80], size=6.0)
        text_page(pdf, plt, "Table 5: refit diagnostics", diag_lines(rows, rts, metas), size=6.4)
        ring_page(pdf, plt, rows, rts)
        hist_pages(pdf, plt, R, good, rtm, pop, egm, rts)
        example_pages(pdf, plt, exd)


def ring_page(pdf, plt, rows, rts):
    fig, axs = plt.subplots(2, 2, figsize=(15, 10.5))
    x = np.arange(NRING)
    for ax, rt in zip(axs.ravel(), rts):
        for j, c in enumerate(CATS + ("all",)):
            cn = "all (primary)" if c == "all" else CNAME[c]
            v = "primary" if c == "all" else PRIMARY[c]
            col = "k" if c == "all" else CCOL[c]
            for p, mk, dx in (("elastic", "o", -0.15), ("moller", "s", 0.15)):
                rr = [get(rows, table="scale", run_type=rt, population=p, ring=r_, category=cn, variant=v) for r_ in RINGL[1:]]
                if rr[0] is None:
                    continue
                s = np.array([100 * r["S"] if r["n_clusters"] else np.nan for r in rr])
                e = np.array([100 * r["S_err"] if r["n_clusters"] else np.nan for r in rr])
                ax.errorbar(x + dx + 0.04 * (j - 2.5), s, yerr=e, fmt=mk, color=col, mfc=col if p == "elastic" else "white", ms=5, capsize=2, lw=1,
                            label=f"{cn} ({p})" if p == "elastic" or rt == "2p24GeV_Ta_X17" else None)
        ax.axhline(0, color="0.5", lw=0.8)
        for yv in (-0.1, 0.1):
            ax.axhline(yv, color="0.7", lw=0.6, ls=":")
        ax.set_xticks(x)
        ax.set_xticklabels(RINGL[1:])
        ax.set_xlabel("seed module radius (mm)")
        ax.set_ylabel("energy-scale shift S (%)")
        ax.set_title(RT_LABEL[rt])
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7, ncol=2)
    fig.suptitle("Energy-scale shift of the calibration samples if each category were corrected (primary variant), per ring\n"
                 "filled = ep elastic window, open = Moller pairs (Ta only); dotted = +-0.1 %; error = statistical", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    pdf.savefig(fig)
    plt.close(fig)


def hist_pages(pdf, plt, R, good, rtm, pop, egm, rts):
    for what in ("cl", "mem"):
        fig, axs = plt.subplots(2, 3, figsize=(16, 10))
        for ax, c in zip(axs.ravel(), CATS):
            for v, ls in zip(VARS[c], ("-", "--")):
                de = R[f"de_{v}"].astype(float)
                den = egm if what == "cl" else R["e_def"].astype(float)
                allx = 100 * de[good & (R["cat"] == c)] / den[good & (R["cat"] == c)]
                allx = allx[np.isfinite(allx)]
                if not allx.size:
                    continue
                lo, hi = np.quantile(allx, (0.01, 0.99))
                lo, hi = min(lo, -1e-3), max(hi, 1e-3)
                bins = np.linspace(lo, hi, 61)
                for rt in rts:
                    for p, lw in (("all", 1.0), ("elastic", 2.0)):
                        m = good & (R["cat"] == c) & (rtm == rt) & pop[p]
                        x = 100 * de[m] / den[m]
                        x = x[np.isfinite(x)]
                        if x.size:
                            ax.hist(np.clip(x, lo, hi), bins=bins, histtype="step", ls=ls, lw=lw, color=RTCOL[rt], density=True,
                                    label=f"{RT_LABEL[rt]} {p} {v} (n {x.size})")
            ax.set_yscale("log")
            ax.set_xlabel("dE / E_cluster (%)" if what == "cl" else "dE / E_member (%)")
            ax.set_title(CNAME[c] + "  (" + ", ".join(f"{v}: {'solid' if i == 0 else 'dashed'}" for i, v in enumerate(VARS[c])) + ")", fontsize=10)
            ax.legend(fontsize=6)
            ax.grid(alpha=0.3)
        axs.ravel()[-1].axis("off")
        fig.suptitle(("dE / E_cluster" if what == "cl" else "dE / E_member (relative change of the member's energy)") +
                     " per category: thin = all clusters, thick = ep elastic window; density, clipped to the 1-99 % range of the category", fontsize=12)
        fig.tight_layout(rect=(0, 0, 1, 0.95))
        pdf.savefig(fig)
        plt.close(fig)


def example_pages(pdf, plt, exd):
    TT = 4.0 * np.arange(100)
    ACOL = {"bs": "tab:blue", "cs": "tab:blue", "st": "tab:green", "k1": "tab:blue"}
    for k0 in range(0, len(exd), 4):
        fig = plt.figure(figsize=(16, 11))
        gs = fig.add_gridspec(4, 2, height_ratios=(3, 1, 3, 1), hspace=0.45)
        for j, e in enumerate(exd[k0:k0 + 4]):
            ax = fig.add_subplot(gs[2 * (j // 2), j % 2])
            axr = fig.add_subplot(gs[2 * (j // 2) + 1, j % 2], sharex=ax)
            Y = e["Y"].astype(float)
            ov = Y >= 4095
            d = e["def"]
            mdef = d["base"] + (np.sum(d["comp"], 0) if d["k"] else 0.0)
            ax.plot(TT[~ov], Y[~ov], "k.", ms=3.5, label="data")
            ax.plot(TT[ov], Y[ov], "rx", ms=5, label="overflow")
            ax.plot(TT, mdef, "r-", lw=1.4, label=f"default k={d['k']}: A " + ", ".join(f"{a:.0f}" for a in d["a"]))
            for cmp in d["comp"]:
                ax.plot(TT, d["base"] + cmp, "r:", lw=0.8)
            axr.plot(TT[~ov], (Y - mdef)[~ov], "r.", ms=3, label="data - default")
            for v, al in e["alt"].items():
                if al["model"] is None:
                    continue
                ex = al["extra"]
                lab = f"{v} ({VLABEL[v]}) k={al['k']}: A " + ", ".join(f"{a:.0f}" for a in al["A"]) + (f", s={ex['s']:.1f}" if "s" in ex else "") + (
                    f", slope {ex['bslope']:.0f}/100ns" if "bslope" in ex else "")
                col = ACOL[v] if not (v == "st" and e["cat"] == 2) else "tab:green"
                ax.plot(TT, al["model"], "-", color=col, lw=1.2, label=lab)
                axr.plot(TT[~ov], (Y - al["model"])[~ov], ".", color=col, ms=3, label=f"data - {v}")
            axr.axhline(0, color="0.6", lw=0.8)
            de = ", ".join(f"dE_{v} {e['de_all'][v]:+.1f} MeV ({100 * e['de_all'][v] / e['e_cluster']:+.2f} % E_cl)" for v in e["de_all"])
            ax.set_title(f"{RT_LABEL[e['run_type']]} run {e['run']} ev {e['evn']} {e['module']} {e['category']} ({e['role']}, {e['population']})\n"
                         f"E_member {e['e_def']:.0f} MeV, E_cluster {e['e_cluster']:.0f} MeV; {de}", fontsize=8)
            ax.legend(fontsize=6, loc="upper right")
            axr.legend(fontsize=6, loc="upper right")
            ax.set_ylabel("ADC")
            axr.set_ylabel("resid. (ADC)")
            axr.set_xlabel("time (ns)")
        fig.suptitle("Examples (elastic / Moller clusters first): default pileup_rates fit (red) vs the category's alternative(s); residuals below", fontsize=11)
        pdf.savefig(fig)
        plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--cache", default=str(SCAN_ROOT / "exotic_p8_cache.npz"), help="pass-1 member records (reused if present; delete to rebuild)")
    ap.add_argument("--jobs", type=int, default=12)
    ap.add_argument("--examples", type=int, default=10, help="example panels per category (0 = none)")
    ap.add_argument("--evio-dir", default=p1.EVIO_DIR)
    ap.add_argument("--only", default="", help="debug: substring of the exotic file names to process (no cache written)")
    a = ap.parse_args()
    jobs = jobs_of(a.evio_dir, a.only)
    if a.only:
        a.cache = ""
    metas, R, C, rs = collect(jobs, a.jobs, a.cache)
    report(metas, R, C, rs, a.out, a.examples, a.evio_dir, a.jobs)


if __name__ == "__main__":
    main()
