#!/usr/bin/env python3
"""exotic_p1_quantities.py -- P1 of the exotic-waveform study: candidate quantities Q1-Q14 for the user's 12 examples and for random PbWO4 channels.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p1_quantities.py [--sample-dir <dir>] [--full-scan-dir <dir>] [--evio-dir <dir>] [--out <dir>] [--cache <pkl>] [--jobs 3]
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p1_quantities.py --plot-only --cache <pkl>        # redo the TSV / PDF from a cached compute stage

Samples (seed SEED), PbWO4 kept channels (rows of the scan npz) of the 2-file rawsum_sample scans, per run type: POP_RAND random (any k), POP_K2 random with k >= 2,
POP_CLIP random with nclip > 0 (the last two drawn from the channels not already chosen).  The 12 examples are decoded from EVIO by event number.
Every channel is refitted from the raw samples with the default fit (pileup_rates switches off) and must reproduce the stored k and chi2/dof.
Quantities (sigma_ref = the module's noise reference of the EVIO split, b_ref its baseline reference; Y = raw samples, DT = 4 ns, 100 samples):
 Q1  dB      (pedestal30 - b_ref) / sigma_ref, pedestal30 = pileup_rates.pedestal (clipped mean of samples 0-29) [= stored ch_zb]
 Q2  slope   least-squares slope of the pulse-free samples, in sigma_ref per 100 ns; slope_sig = slope / its statistical error (white noise sigma_ref).
             Pulse-free = outside the union of the fit windows (32 samples, 3 before the onset) of every fitted pulse and not overflow; needs >= 10.
 Q3  noise   sqrt(mean (Y - line)^2) over those samples / sigma_ref
 Q4  tail    (mean of Y[90:100] - pedestal30) / sigma_ref; NaN unless >= 8 of those 10 samples are pulse-free
 Q5  w50     FWHM(Y - fit baseline) / FWHM(template at the main pulse's fitted onset), main = largest A [= stored ch_w50]
 Q6  rawsum  sum of (Y - fit baseline) over the template walk window (pileup_rates.raw_sums def. 1, on the smoothed main component) / (A_main I_m), I_m = iint_wa
             of the module's template [stored ch_rs_raw1 / (A I_m)]; for k >= 2 the sum still contains the other pulses
 Q7  span    max - min onset of the fitted pulses (ns), 0 for k = 1
 Q8  r       largest A_other / A_main (0 for k = 1)
 Q9  rho1    lag-1 autocorrelation sum z_i z_(i+1) / sum z_i^2 of the normalised residual z = (Y - model) / sqrt(V) of a single-pulse refit over its fit window
             (pairs of adjacent non-overflow samples); refit = MF.fit_group with one pulse started at the main onset (amplitude, onset, baseline free)
 Q10 chi1    chi2/dof of that single-pulse refit (same V, mask and dof = n - 3 as the stored chi2); chi1_ratio = chi1 / stored chi2/dof
 Q11 nclip   samples at 4095
 Q12 dclip   signed distance (ns) from the nearest non-main fitted onset to the clipped stretch [first, last overflow sample] x DT: onset left of it: first - t,
             right of it: t - last, inside: negative (-min(t - first, last - t)); NaN without overflow or without a second pulse
 Q13 E/Eb    A_main x gain / run beam energy
 Q14 Esec/Eb largest A_other x gain / run beam energy (0 for k = 1)
Outputs (--out): exotic_P1_quantities_20261117.pdf / .tsv / _examples.tsv.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import time
import pickle
import argparse
import tempfile
import multiprocessing as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402

PYDIR = HERE.parent.parent.parent
EVIO_DIR = "/Users/L00338853/work/PRad/data/evio"
STEM = "exotic_P1_quantities_20261117"
SEED = 20261117
POP_RAND, POP_K2, POP_CLIP = 500, 200, 100
TYPES = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]
EXCLUDE_RUNS = frozenset({24937})                    # user decision (2026-11-18): run 24937 is not a good run; excluded from the P4/P5/P6 inputs (its files stay on disk)
EXCL_NOTE = "run 24937 excluded (user: not a good run)"
RT_LABEL = {"0p73GeV_H2": "0.73 GeV H2", "2p24GeV_H2": "2.24 GeV H2", "3p49GeV_H2": "3.49 GeV H2", "2p24GeV_Ta_X17": "2.24 GeV Ta X17"}
# the user's examples: (id, exotic type, run, event, module, run type)
EXAMPLES = [
    ("1a", 1, 24603, 29267174, "W494", "0p73GeV_H2"),
    ("1b", 1, 24603, 29382036, "W634", "0p73GeV_H2"),
    ("1c", 1, 24911, 76390786, "W565", "3p49GeV_H2"),
    ("1d", 1, 26138, 20705487, "W494", "2p24GeV_Ta_X17"),
    ("2a", 2, 24590, 4084847, "W121", "0p73GeV_H2"),
    ("2b", 2, 24590, 6138838, "W818", "0p73GeV_H2"),
    ("2c", 2, 24590, 4060795, "W818", "0p73GeV_H2"),
    ("3a", 3, 24911, 45926822, "W775", "3p49GeV_H2"),
    ("3b", 3, 25299, 130824709, "W1118", "2p24GeV_H2"),
    ("3c", 3, 25228, 39246232, "W121", "2p24GeV_H2"),
    ("4a", 4, 25257, 32598460, "W877", "2p24GeV_H2"),
    ("4b", 4, 25228, 78491010, "W1000", "2p24GeV_H2"),
]
QCOLS = ["dB", "slope", "slope_sig", "noise", "tail", "w50", "rawsum", "span", "r", "rho1", "chi1", "chi1_ratio", "nclip", "dclip", "E_Eb", "Esec_Eb"]
QNAME = {"dB": "Q1 dB", "slope": "Q2 slope (sigma_ref / 100 ns)", "slope_sig": "Q2 slope significance", "noise": "Q3 noise ratio", "tail": "Q4 tail (sigma_ref)",
         "w50": "Q5 w50 ratio", "rawsum": "Q6 raw sum / (A I_m)", "span": "Q7 onset span (ns)", "r": "Q8 r = A_sec / A_main", "rho1": "Q9 lag-1 autocorr.",
         "chi1": "Q10 single-pulse chi2/dof", "chi1_ratio": "Q10 chi2_1 / chi2", "nclip": "Q11 nclip", "dclip": "Q12 dist to clip (ns)", "E_Eb": "Q13 E / E_beam",
         "Esec_Eb": "Q14 E_sec / E_beam"}
NS, DT = 100, 4.0
OVER = 4095


# ---------------------------------------------------------------------------------------------------------------- quantities
class Env:
    """Per-process cache of event_display Ctx and per-run-type iint."""
    ctxs: dict = {}
    iint: dict = {}


def get_ctx(run):
    import event_display as ed
    if run not in Env.ctxs:
        Env.ctxs[run] = ed.Ctx(run)
        ctx = Env.ctxs[run]
        z = np.load(ed.pr.BANK_ROOT / ctx.rt / f"templates_{ctx.rt}.npz", allow_pickle=False)
        d = {str(n): float(v) for n, v, ok in zip(z["names"], z["iint_wa"], z["valid"]) if ok > 0}
        med = float(np.median(list(d.values())))
        ctx.I = np.array([d.get(str(n), med) for n in ctx.names])
    return Env.ctxs[run]


def single_fit(ctx, Y, mi, sg, br, t_main):
    """Single-pulse refit started at t_main (MF.fit_group).  Returns chi2/dof, lag-1 autocorrelation of the normalised residual, residual array (ADC), window mask."""
    import pileup_rates as pr
    F32 = np.float32
    Yf = Y.astype(F32)[None]
    ped = pr.pedestal(Yf)
    far = (np.abs(ped - br) > pr.START_NSIGMA * sg) if np.isfinite(br) else np.zeros(1, bool)
    start = np.where(far, br, ped).astype(F32)
    y = Yf - start[:, None]
    W = (Yf < pr.OVERFLOW).astype(F32)
    clipped = bool((Yf >= pr.OVERFLOW).any())
    w = W if clipped else None
    if clipped:
        y = y * W
    m = ctx.bank.row[[mi]]
    sig = np.array([sg], F32)
    t0, a, v, model, db = ctx.mf.fit_group(y, w, m, sig, np.array([[t_main]], F32), np.array([[True]]), NS)
    if not v[0, 0]:
        return np.nan, np.nan, None, None
    mask = (pr.MF.wmask(t0[:, 0], NS) & (W > 0))[0]
    noise = np.sqrt(sig.astype(np.float64) ** 2 + 1.0 / 12.0)
    var = ctx.mf.variance(m, t0.astype(F32), a, v, noise, NS)[0]
    res = (y - model)[0].astype(np.float64)
    dof = int(mask.sum()) - 3
    chi = float((mask * res ** 2 / var).sum() / dof) if dof > 0 else np.nan
    z = res / np.sqrt(var)
    pair = mask[:-1] & mask[1:]
    rho = float((z[:-1] * z[1:])[pair].sum() / max((z[:-1] ** 2)[pair].sum(), 1e-12)) if pair.sum() >= 6 else np.nan
    return chi, rho, res, mask


def quantities(ctx, r, i, b_ref, sig_ref):
    """Q1-Q14 of channel row i of ed.fit_event result r; returns (dict, extras for plots)."""
    import pileup_rates as pr
    mi = int(r["mi"][i])
    Y = r["Y"][i].astype(np.float64)
    k = int(r["k"][i])
    a, t = r["a"][i][:k].astype(np.float64), r["t"][i][:k].astype(np.float64)
    comp = r["comp"][i][:k]
    br, sg = float(b_ref[mi]), float(sig_ref[mi])
    ped30 = float(pr.pedestal(Y[None])[0])
    ov = Y >= OVER
    q = {n: np.nan for n in QCOLS}
    q["dB"] = (ped30 - br) / sg
    t0 = np.zeros((1, pr.KMAX))
    v = np.zeros((1, pr.KMAX), bool)
    t0[0, :k], v[0, :k] = t, True
    win = pr.MF.windows(t0, v, NS)[0]
    free = ~win & ~ov
    idx = np.arange(NS)
    if free.sum() >= 10:
        tt = DT * idx[free]
        p = np.polyfit(tt, Y[free], 1)
        res = Y[free] - np.polyval(p, tt)
        q["slope"] = p[0] * 100.0 / sg
        q["slope_sig"] = q["slope"] / (100.0 / np.sqrt(((tt - tt.mean()) ** 2).sum()))
        q["noise"] = float(np.sqrt((res ** 2).mean()) / sg)
    last = slice(90, 100)
    if free[last].sum() >= 8:
        q["tail"] = float((Y[last][free[last]] - ped30).mean() / sg)
    main = int(np.argmax(a))
    am, tm = a[main], t[main]
    bank = ctx.bank
    f = dict(a=r["a"][i:i + 1], t=r["t"][i:i + 1], comp=r["comp"][i:i + 1])
    w50, _ = pr.shape_vars(Y[None], np.array([r["ped"][i]]), bank.row[[mi]], f, bank)
    q["w50"] = float(w50[0])
    raw1, _, _, _ = pr.raw_sums((Y - r["ped"][i])[None], comp[main][None].astype(np.float64), np.array([tm]))
    q["rawsum"] = float(raw1[0] / (am * ctx.I[mi]))
    oth = np.arange(k) != main
    q["span"] = float(t.max() - t.min())
    q["r"] = float(a[oth].max() / am) if k > 1 else 0.0
    chi1, rho, res1, mask1 = single_fit(ctx, Y, mi, sg, br, tm)
    q["rho1"], q["chi1"] = rho, chi1
    q["chi1_ratio"] = chi1 / float(r["chi2"][i]) if np.isfinite(chi1) and r["chi2"][i] > 0 else np.nan
    q["nclip"] = int(ov.sum())
    if ov.any() and k > 1:
        c0, c1 = DT * idx[ov].min(), DT * idx[ov].max()
        d = np.maximum(c0 - t[oth], t[oth] - c1)
        q["dclip"] = float(d[np.argmin(d)])
    q["E_Eb"] = float(am * ctx.gain[mi] / ctx.beam)
    q["Esec_Eb"] = float(a[oth].max() * ctx.gain[mi] / ctx.beam) if k > 1 else 0.0
    return q, dict(res1=res1, mask1=mask1)


# ---------------------------------------------------------------------------------------------------------------- worker
def work(task):
    """One task: ('rand', run, files [(path, split, b_ref, sig_ref)], picks {ev: {midx: stored}}) or ('ex', run, evs [(id, ev, midx)], full-scan npz path)."""
    import event_display as ed
    from prad2py import dec
    import pileup_rates as pr
    kind, run = task[0], task[1]
    t_start = time.time()
    ctx = get_ctx(run)
    out = []
    if kind == "rand":
        files, picks = task[2], task[3]
        left = set(picks)
        for f, split, b_ref, sig_ref in files:
            ch = ctx.open(f)
            while left and ch.read() == dec.Status.success:
                if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
                    continue
                for i in range(ch.get_n_events()):
                    ch.select_event(i)
                    en = int(ch.info().event_number)
                    if en not in picks or en not in left:
                        continue
                    left.discard(en)
                    d = ch.decode_event(i, with_ssp=False)
                    if not d["ok"]:
                        continue
                    r = ed.fit_event(ctx, d["event"], b_ref, sig_ref)
                    for midx, st in picks[en].items():
                        row = r["row"].get(midx)
                        rec = dict(run=run, rt=ctx.rt, ev=en, midx=midx, name=str(ctx.names[midx]), split=split, pop=st["pop"], ok=False, **{n: np.nan for n in QCOLS})
                        if row is not None:
                            rec.update(k=int(r["k"][row]), chi2=float(r["chi2"][row]), amain=float(r["a"][row].max()), st_k=st["k"], st_chi2=st["chi2"])
                            rec["ok"] = rec["k"] == st["k"] and bool(np.isclose(rec["chi2"], st["chi2"], rtol=2e-3, atol=1e-3))
                            q, _ = quantities(ctx, r, row, b_ref, sig_ref)
                            rec.update(q)
                        out.append(rec)
            ch.close()
        for en in left:
            for midx, st in picks[en].items():
                out.append(dict(run=run, rt=ctx.rt, ev=en, midx=midx, name=str(ctx.names[midx]), pop=st["pop"], ok=False, missing=True, **{n: np.nan for n in QCOLS}))
    else:
        evs, scan_path = task[2], task[3]
        z = np.load(scan_path)
        splits = [int(s) for s in z["splits"]]
        evnum = z["ev_number"]
        for ident, en, mname in evs:
            midx = int(ctx.index[mname])
            f, si, files = ed.locate(ctx, EVIO_DIR, en)
            spl = int(f.rsplit(".", 1)[1])
            if spl in splits:
                b_ref, sig_ref = z["mod_b_ref"][splits.index(spl)], z["mod_sig_ref"][splits.index(spl)]
                refsrc = "scan"
            else:
                _, b_ref, sig_ref = ed.noise_reference(ctx, f)
                refsrc = "computed"
            fe = ed.find_event(ctx, f, en)
            rec = dict(id=ident, run=run, rt=ctx.rt, ev=en, midx=midx, name=mname, split=spl, refsrc=refsrc, ok=False, **{n: np.nan for n in QCOLS})
            if fe is None:
                out.append(rec)
                continue
            r = ed.fit_event(ctx, fe, b_ref, sig_ref)
            row = r["row"].get(midx)
            if row is None:
                out.append(rec)
                continue
            rec.update(k=int(r["k"][row]), chi2=float(r["chi2"][row]), amain=float(r["a"][row].max()))
            ie = np.nonzero(evnum == en)[0]
            st = None
            if ie.size:
                sel = np.nonzero((z["ch_ev"] == ie[0]) & (z["ch_midx"] == midx))[0]
                if sel.size:
                    st = (int(z["ch_k"][sel[0]]), float(z["ch_chi2"][sel[0]]))
            rec["st_k"], rec["st_chi2"] = st if st else (None, None)
            rec["ok"] = st is not None and rec["k"] == st[0] and bool(np.isclose(rec["chi2"], st[1], rtol=2e-3, atol=1e-3))
            q, ex = quantities(ctx, r, row, b_ref, sig_ref)
            rec.update(q)
            kk = rec["k"]
            rec["plot"] = dict(Y=r["Y"][row].astype(float), ped=float(r["ped"][row]), a=r["a"][row][:kk].copy(), t=r["t"][row][:kk].copy(), comp=r["comp"][row][:kk].copy(),
                               res1=ex["res1"], mask1=ex["mask1"], gain=float(ctx.gain[midx]), beam=float(ctx.beam), b_ref=float(b_ref[midx]), sig_ref=float(sig_ref[midx]))
            out.append(rec)
    print(f"[work] {kind} run {run}: {len(out)} records, {time.time() - t_start:.0f} s", flush=True)
    return out


# ---------------------------------------------------------------------------------------------------------------- sampling
def build_tasks(sample_dir, full_dir):
    rng = np.random.default_rng(SEED)
    runs = {}
    for fn in sorted(Path(sample_dir).glob("run_*.npz")):
        z = np.load(fn)
        run = int(z["run"])
        pw = z["mod_pwo4"][z["ch_midx"].astype(np.int64)]
        runs[run] = dict(rt=str(z["run_type"]), pw=pw, k=z["ch_k"], nclip=z["ch_nclip"], midx=z["ch_midx"], ev=z["ch_ev"], evnum=z["ev_number"],
                         chi2=z["ch_chi2"], splits=[int(s) for s in z["splits"]], b=z["mod_b_ref"], s=z["mod_sig_ref"], path=fn)
    chosen = {run: {} for run in runs}          # run -> {(ev number, midx): stored}
    summary = {}
    for rt in TYPES:
        rr = [r for r in sorted(runs) if runs[r]["rt"] == rt]
        taken = set()
        for pop, n, cond in (("rand", POP_RAND, lambda d: d["pw"]), ("k2", POP_K2, lambda d: d["pw"] & (d["k"] >= 2)), ("clip", POP_CLIP, lambda d: d["pw"] & (d["nclip"] > 0))):
            pools = []
            for run in rr:
                d = runs[run]
                idx = np.nonzero(cond(d))[0]
                pools.append((run, idx))
            tot = sum(p[1].size for p in pools)
            summary[(rt, pop)] = tot
            cum = np.cumsum([0] + [p[1].size for p in pools])
            got = 0
            for g in rng.permutation(tot):
                if got >= n:
                    break
                j = int(np.searchsorted(cum, g, side="right") - 1)
                run, idx = pools[j]
                row = int(idx[g - cum[j]])
                d = runs[run]
                key = (run, row)
                if key in taken:
                    continue
                taken.add(key)
                chosen[run][(int(d["evnum"][d["ev"][row]]), int(d["midx"][row]))] = dict(pop=pop, k=int(d["k"][row]), chi2=float(d["chi2"][row]))
                got += 1
    tasks = []
    for run, c in chosen.items():
        d = runs[run]
        picks = {}
        for (en, mi), st in c.items():
            picks.setdefault(en, {})[mi] = st
        files = [(str(Path(EVIO_DIR) / f"prad_{run:06d}.evio.{s:05d}"), s, d["b"][j], d["s"][j]) for j, s in enumerate(d["splits"])]
        tasks.append(("rand", run, files, picks))
    ex_runs = {}
    for ident, typ, run, en, mod, rt in EXAMPLES:
        ex_runs.setdefault(run, []).append((ident, en, mod))
    ex = [("ex", run, evs, str(Path(full_dir) / f"run_{run}.npz")) for run, evs in ex_runs.items()]
    return tasks, ex, summary


def compute(a):
    t0 = time.time()
    tasks, ex, summary = build_tasks(a.sample_dir, a.full_scan_dir)
    print(f"[sample] pools {summary}", flush=True)
    allt = sorted(tasks, key=lambda x: -sum(len(v) for v in x[3].values())) + ex
    with mp.get_context("fork").Pool(a.jobs) as pool:
        res = pool.map(work, allt, chunksize=1)
    rows = [r for rr, t in zip(res, allt) if t[0] == "rand" for r in rr]
    exs = [r for rr, t in zip(res, allt) if t[0] == "ex" for r in rr]
    exs.sort(key=lambda r: [e[0] for e in EXAMPLES].index(r["id"]))
    return dict(rows=rows, exs=exs, summary=summary, elapsed=time.time() - t0)


# ---------------------------------------------------------------------------------------------------------------- output
def fmt(x):
    if isinstance(x, (bool, np.bool_)):
        return str(int(x))
    if isinstance(x, (float, np.floating)):
        return "nan" if not np.isfinite(x) else f"{x:.5g}"
    return str(x)


def write_tsv(path, rows, cols):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(fmt(r.get(c, "")) for c in cols) + "\n")


# ---------------------------------------------------------------------------------------------------------------- cuts
def _f(r, c):
    v = r.get(c, np.nan)
    return np.nan if v is None else v


# proposed cuts: key -> (description, function(record) -> bool; NaN never flags)
CUT = {
    1: ("dB > 8 and tail < -5 sigma_ref", lambda r: _f(r, "dB") > 8 and _f(r, "tail") < -5),
    2: ("w50 ratio > 3", lambda r: _f(r, "w50") > 3),
    3: ("rawsum > 2.2 and chi1_ratio > 50", lambda r: _f(r, "rawsum") > 2.2 and _f(r, "chi1_ratio") > 50),
    4: ("nclip >= 1 and E_sec/E_beam > 0.5", lambda r: _f(r, "nclip") >= 1 and _f(r, "Esec_Eb") > 0.5),
}
# per type: histogram quantities (name, log x, cut line or None), scatter (x, y, log x, log y), cut lines on the scatter (x, y)
PAGES = {
    1: dict(title="Type 1: noisy descending low-energy forms", hist=[("dB", False, 8), ("tail", False, -5), ("slope_sig", False, None)], sc=("dB", "tail", False, False), cl=(8, -5)),
    2: dict(title="Type 2: smooth / noisy wide pulse", hist=[("w50", False, 3), ("rawsum", False, None), ("dB", False, None)], sc=("w50", "dB", False, False), cl=(3, None)),
    3: dict(title="Type 3: smooth peak region, multiple close pulses (over-fit)", hist=[("rawsum", False, 2.2), ("chi1_ratio", True, 50), ("rho1", False, None)], sc=("rawsum", "chi1_ratio", False, True), cl=(2.2, 50)),
    4: dict(title="Type 4: clipped at 4095, described by several pulses", hist=[("nclip", False, 1), ("Esec_Eb", False, 0.5), ("E_Eb", False, None)], sc=("dclip", "Esec_Eb", False, False), cl=(None, 0.5)),
}
EX_TYPE = {e[0]: e[1] for e in EXAMPLES}
POPCOL = {"rand": "0.45", "k2": "tab:blue", "clip": "tab:orange"}
TCOL = dict(zip(TYPES, ["tab:blue", "tab:green", "tab:purple", "tab:orange"]))


def flag_all(recs):
    for r in recs:
        for ty, (_, fn) in CUT.items():
            r[f"c{ty}"] = bool(fn(r))


def full_sample(sample_dir):
    """Fractions over ALL kept PbWO4 channels of the sample scans (stored keys only) per run type, for the cut parts the stored keys can express."""
    import pileup_rates as pr
    acc = {rt: dict(n=0, dB8=0, w3=0, rs22_k5=0, clip=0, clip_k2=0, clip_c4=0, ret_clip=0, ret_c4=0) for rt in TYPES}
    for fn in sorted(Path(sample_dir).glob("run_*.npz")):
        z = np.load(fn)
        rt = str(z["run_type"])
        names = [str(n) for n in z["mod_name"]]
        cfg = pr.run_config(int(z["run"]))
        gain = pr.module_gain(cfg, rt, names)
        zz = np.load(pr.BANK_ROOT / rt / f"templates_{rt}.npz", allow_pickle=False)
        d = {str(n): float(v) for n, v, ok in zip(zz["names"], zz["iint_wa"], zz["valid"]) if ok > 0}
        med = float(np.median(list(d.values())))
        iint = np.array([d.get(n, med) for n in names])
        mi = z["ch_midx"].astype(np.int64)
        pw = z["mod_pwo4"][mi]
        k, a, nc = z["ch_k"], z["ch_a"].astype(np.float64), z["ch_nclip"]
        am = a.max(axis=1)
        sec = np.where(np.arange(6)[None, :] < k[:, None], a, 0.0)
        sec[np.arange(k.size), a.argmax(axis=1)] = 0.0
        esec = sec.max(axis=1) * gain[mi] / cfg["beam"]
        with np.errstate(invalid="ignore", divide="ignore"):
            rs = z["ch_rs_raw1"] / (am * iint[mi])
            c = acc[rt]
            c["n"] += int(pw.sum())
            c["dB8"] += int((pw & (z["ch_zb"] > 8)).sum())
            c["w3"] += int((pw & (z["ch_w50"] > 3)).sum())
            c["rs22_k5"] += int((pw & (rs > 2.2) & (k >= 5)).sum())
            c["clip"] += int((pw & (nc >= 1)).sum())
            c["clip_k2"] += int((pw & (nc >= 1) & (k >= 2)).sum())
            c["clip_c4"] += int((pw & (nc >= 1) & (esec > 0.5)).sum())
            c["ret_clip"] += int((pw & (nc >= 1) & (nc < pr.SAT_NCLIP)).sum())
            c["ret_c4"] += int((pw & (nc >= 1) & (nc < pr.SAT_NCLIP) & (esec > 0.5)).sum())
        print(f"[full] run {z['run']} {rt}", flush=True)
    return acc


# ---------------------------------------------------------------------------------------------------------------- plots
def ex_page(pdf, plt, exs):
    n = len(exs)
    fig = plt.figure(figsize=(44, 16))
    outer = fig.add_gridspec(2, 6, hspace=0.32, wspace=0.12, left=0.03, right=0.995, top=0.86, bottom=0.02)
    fig.text(0.03, 0.975, f"P1: the {n} examples, refitted from EVIO with the default fit (verified against the stored scan k and chi2/dof)", fontsize=22, fontweight="bold")
    fig.text(0.03, 0.945, "black dots: raw samples (red x: overflow); red: total fit; dashed colours: fitted pulses in time order; dotted green: single-pulse refit (Q9, Q10); grey: fit baseline; "
             "shaded: fit windows.  Lower strip: data - fit.  Text: Q1-Q14 (definitions in the script header and on the table page).", fontsize=12)
    for j, e in enumerate(exs):
        gs = outer[j // 6, j % 6].subgridspec(4, 1, height_ratios=[3, 1, 0.45, 1.5], hspace=0.25)
        ax, axr, axt = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[3])
        p = e["plot"]
        Y, ov = p["Y"], p["Y"] >= OVER
        tt = DT * np.arange(NS)
        ax.plot(tt[~ov], Y[~ov], "k.", ms=5)
        if ov.any():
            ax.plot(tt[ov], Y[ov], "x", color="tab:red", ms=7, mew=1.6)
        ax.axhline(p["ped"], color="0.6", lw=0.9)
        tot = p["ped"] + p["comp"].sum(0)
        ax.plot(tt, tot, "r-", lw=1.4)
        for c, comp in zip(("tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:brown", "tab:pink"), p["comp"]):
            ax.plot(tt, p["ped"] + comp, "--", color=c, lw=0.9)
        if p["res1"] is not None:
            ax.plot(tt[~ov], (Y - p["res1"])[~ov], ":", color="darkgreen", lw=1.5)
        for tp in p["t"]:
            ax.axvspan(tp - 3 * DT + 0.0, tp + 29 * DT, color="0.85", alpha=0.18, lw=0)
        yl = (min(Y.min(), tot.min()), Y.max())
        ax.set_ylim(yl[0] - 0.05 * (yl[1] - yl[0]), min(yl[1] * 1.08 + 5, 4500))
        ax.set_xlim(0, DT * NS)
        ax.tick_params(labelsize=9, labelbottom=False)
        ax.set_title(f"{e['id']}  (type {EX_TYPE[e['id']]})  run {e['run']}  ev {e['ev']}  {e['name']}  {e['rt']}\nk = {e['k']}  chi2/dof = {e['chi2']:.2f}  "
                     f"stored {e['st_k']} / {e['st_chi2']:.2f}  {'OK' if e['ok'] else 'NOT REPRODUCED'}   A_main {e['amain']:.0f}  E = {e['E_Eb'] * p['beam']:.0f} MeV",
                     fontsize=11, loc="left", pad=4)
        res = Y - tot
        axr.plot(tt[~ov], res[~ov], "k.-", lw=0.6, ms=3)
        axr.axhline(0, color="r", lw=0.8)
        axr.set_xlim(0, DT * NS)
        axr.tick_params(labelsize=9)
        axr.set_xlabel("t (ns)", fontsize=9)
        axt.axis("off")
        q = {c: e[c] for c in QCOLS}
        f = lambda c, s: f"{q[c]:{s}}" if np.isfinite(q[c]) else "nan"      # noqa: E731
        t1 = (f"Q1 dB {f('dB', '.1f')}   Q2 slope {f('slope', '.2f')} ({f('slope_sig', '.1f')} sig)   Q3 noise {f('noise', '.2f')}   Q4 tail {f('tail', '.1f')}\n"
              f"Q5 w50 {f('w50', '.2f')}   Q6 rawsum {f('rawsum', '.2f')}   Q7 span {f('span', '.0f')} ns   Q8 r {f('r', '.2f')}\n"
              f"Q9 rho1 {f('rho1', '.2f')}   Q10 chi1 {f('chi1', '.3g')} (x{f('chi1_ratio', '.3g')})   Q11 nclip {int(q['nclip'])}   Q12 dclip {f('dclip', '.1f')} ns\n"
              f"Q13 E/Eb {f('E_Eb', '.3f')}   Q14 Esec/Eb {f('Esec_Eb', '.3f')}   sigma_ref {p['sig_ref']:.2f}  b_ref {p['b_ref']:.1f}")
        axt.text(0.0, 1.0, t1, transform=axt.transAxes, fontsize=9.5, va="top", family="monospace", linespacing=1.5)
    pdf.savefig(fig)
    plt.close(fig)


def type_page(pdf, plt, ty, rows, exs, acc):
    cfg = PAGES[ty]
    R = [r for r in rows if r["pop"] == "rand"]
    K2 = [r for r in rows if r["pop"] == "k2"]
    CL = [r for r in rows if r["pop"] == "clip"]
    E = [e for e in exs if EX_TYPE[e["id"]] == ty]
    fig = plt.figure(figsize=(20, 14))
    gs = fig.add_gridspec(3, 2, width_ratios=[1, 1.15], height_ratios=[1, 1, 1], hspace=0.42, wspace=0.18, left=0.06, right=0.985, top=0.9, bottom=0.05)
    fig.text(0.02, 0.965, cfg["title"], fontsize=20, fontweight="bold")
    fig.text(0.02, 0.94, f"examples (red): {', '.join(e['id'] for e in E)}.   Histograms: the random PbWO4 sample ({POP_RAND} channels per run type, any k; step lines, run types as in the legend).   "
             f"Proposed cut: {CUT[ty][0]} (black dashed).", fontsize=11)
    for i, (qn, logx, cut) in enumerate(cfg["hist"]):
        ax = fig.add_subplot(gs[i, 0])
        allv = np.array([r[qn] for r in R] + [e[qn] for e in E], float)
        allv = allv[np.isfinite(allv)]
        if logx:
            allv = allv[allv > 0]
            lo, hi = np.log10(np.percentile(allv, 0.3)), np.log10(allv.max()) + 0.05
            bins = np.logspace(lo, hi, 41)
        else:
            lo, hi = np.percentile(allv, 0.3), np.percentile(allv, 99.7)
            lo, hi = min(lo, min(e[qn] for e in E if np.isfinite(e[qn]))), max(hi, max(e[qn] for e in E if np.isfinite(e[qn])))
            if qn in ("nclip",):
                bins = np.arange(-0.5, max(hi, 3) + 1.5, 1.0)
            else:
                pad = 0.03 * (hi - lo)
                bins = np.linspace(lo - pad, hi + pad, 61)
        for rt in TYPES:
            v = np.array([r[qn] for r in R if r["rt"] == rt], float)
            v = v[np.isfinite(v)]
            ax.hist(np.clip(v, bins[0], bins[-1]), bins=bins, histtype="step", lw=1.6, color=TCOL[rt], label=f"{RT_LABEL[rt]} (n={v.size})")
        ax.set_yscale("log")
        if logx:
            ax.set_xscale("log")
        ymax = ax.get_ylim()[1]
        ax.set_ylim(0.5, ymax * 2.5)
        for e in E:
            if np.isfinite(e[qn]):
                x = min(max(e[qn], bins[0]), bins[-1])
                ax.axvline(x, color="red", lw=1.4)
                ax.text(x, ymax * 1.4, e["id"], color="red", fontsize=11, ha="center", va="bottom", fontweight="bold")
        if cut is not None:
            ax.axvline(cut, color="k", ls="--", lw=1.6)
        ax.set_xlabel(QNAME[qn], fontsize=12)
        ax.set_ylabel("channels", fontsize=11)
        ax.tick_params(labelsize=10)
        if i == 0:
            ax.legend(fontsize=9, loc="upper left", frameon=False)
    # scatter
    xn, yn, lx, ly = cfg["sc"]
    ax = fig.add_subplot(gs[0:2, 1])
    for pop, lab, rr in (("rand", f"random (any k), {len(R)}", R), ("k2", f"random k >= 2, {len(K2)}", K2), ("clip", f"random nclip > 0, {len(CL)}", CL)):
        x = np.array([r[xn] for r in rr], float)
        y = np.array([r[yn] for r in rr], float)
        ok = np.isfinite(x) & np.isfinite(y)
        if lx:
            ok &= x > 0
        if ly:
            ok &= y > 0
        ax.plot(x[ok], y[ok], "o", ms=3.2, alpha=0.45, color=POPCOL[pop], mew=0, label=lab)
    for e in E:
        if np.isfinite(e[xn]) and np.isfinite(e[yn]):
            ax.plot(e[xn], e[yn], "o", ms=9, mfc="none", mec="red", mew=2)
            ax.annotate(e["id"], (e[xn], e[yn]), xytext=(6, 6), textcoords="offset points", color="red", fontsize=12, fontweight="bold")
    base = CL if ty == 4 else R               # axis range: the random population (the clipped one for type 4, which has no random dclip) and the examples; other points outside are cut off
    for axis_, name_, lg in (("x", xn, lx), ("y", yn, ly)):
        v = np.array([r[name_] for r in base] + [e[name_] for e in E], float)
        v = v[np.isfinite(v) & ((v > 0) if lg else True)]
        lo, hi = np.percentile(v, 0.3), np.percentile(v, 99.7)
        lo, hi = min(lo, np.nanmin([e[name_] for e in E])), max(hi, np.nanmax([e[name_] for e in E]))
        pad = (hi / lo) ** 0.08 if lg else None
        lim = (lo / pad, hi * pad) if lg else (lo - 0.08 * (hi - lo), hi + 0.08 * (hi - lo))
        (ax.set_xlim if axis_ == "x" else ax.set_ylim)(*lim)
    if lx:
        ax.set_xscale("log")
    if ly:
        ax.set_yscale("log")
    cx, cy = cfg["cl"]
    if cx is not None:
        ax.axvline(cx, color="k", ls="--", lw=1.5)
    if cy is not None:
        ax.axhline(cy, color="k", ls="--", lw=1.5)
    ax.set_xlabel(QNAME[xn], fontsize=12)
    ax.set_ylabel(QNAME[yn], fontsize=12)
    ax.legend(fontsize=9, loc="best")
    ax.tick_params(labelsize=10)
    # flag table
    axt = fig.add_subplot(gs[2, 1])
    axt.axis("off")
    key = f"c{ty}"
    lines = [f"cut {ty}: {CUT[ty][0]}", f"examples flagged: {sum(e[key] for e in E)} of {len(E)}  ({', '.join(e['id'] for e in E if e[key]) or '-'})", "",
             f"{'run type':18s}{'random':>16s}{'k>=2 pop':>16s}{'nclip>0 pop':>16s}"]
    for rt in TYPES:
        def cell(pop):
            rr = [r for r in rows if r["rt"] == rt and r["pop"] == pop]
            n = len(rr)
            if n < 5:
                return f"{'n=' + str(n):>16s}"
            m = sum(r[key] for r in rr)
            return f"{m:4d}/{n:<4d} {100 * m / n:5.1f}%".rjust(16)
        lines.append(f"{RT_LABEL[rt]:18s}{cell('rand')}{cell('k2')}{cell('clip')}")
    nr = sum(r[key] for r in R)
    lines += [f"{'all (random)':18s}{nr:>9d}/{len(R):<4d} {100 * nr / len(R):5.2f}%", "", "all kept PbWO4 channels of the 2-file sample, stored keys only:"]
    pr_ = {1: ("dB > 8 (zb)", "dB8"), 2: ("w50 > 3", "w3"), 3: ("k >= 5 and rawsum > 2.2 (proxy: no chi1 ratio)", "rs22_k5"), 4: ("nclip >= 1 and E_sec/E_b > 0.5", "clip_c4")}[ty]
    lines.append(f"  {pr_[0]}:")
    for rt in TYPES:
        c = acc[rt]
        lines.append(f"    {RT_LABEL[rt]:18s}{c[pr_[1]]:>9d} / {c['n']:<9d} = {100 * c[pr_[1]] / c['n']:.4f} %")
    if ty == 4:
        lines.append("  clipped channels: nclip>=1 share of all, share with k>=2,")
        lines.append("  and the chain-kept part (nclip < SAT_NCLIP) with its cut-flagged count:")
        for rt in TYPES:
            c = acc[rt]
            lines.append(f"    {RT_LABEL[rt]:18s}nclip>=1 {100 * c['clip'] / c['n']:.3f} %   k>=2 {100 * c['clip_k2'] / max(c['clip'], 1):.1f} % of them   chain-kept {c['ret_clip']}, flagged {c['ret_c4']}")
    axt.text(0.0, 1.0, "\n".join(lines), transform=axt.transAxes, fontsize=9.6, va="top", family="monospace", linespacing=1.4)
    pdf.savefig(fig)
    plt.close(fig)


def table_page(pdf, plt, exs, elapsed, nfail, summary):
    fig = plt.figure(figsize=(30, 14))
    fig.text(0.02, 0.96, "Q values of the examples and the proposed cuts", fontsize=20, fontweight="bold")
    hdr = ["id", "type", "run", "event", "module", "run type", "k", "chi2/dof"] + [QNAME[c].split(" ", 1)[0] + " " + c for c in QCOLS] + ["cut1", "cut2", "cut3", "cut4"]
    cells = []
    for e in exs:
        row = [e["id"], str(EX_TYPE[e["id"]]), str(e["run"]), str(e["ev"]), e["name"], e["rt"], str(e["k"]), f"{e['chi2']:.2f}"]
        for c in QCOLS:
            v = e[c]
            row.append("nan" if not np.isfinite(v) else (f"{v:.3g}" if abs(v) < 1e4 else f"{v:.2e}"))
        row += ["Y" if e[f"c{t}"] else "." for t in (1, 2, 3, 4)]
        cells.append(row)
    ax = fig.add_axes([0.01, 0.42, 0.98, 0.5])
    ax.axis("off")
    tb = ax.table(cellText=cells, colLabels=hdr, loc="upper center", cellLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(8.5)
    tb.scale(1, 1.9)
    txt = ("Proposed cuts (PbWO4 only; NaN never flags):\n" + "\n".join(f"  cut {t}: {CUT[t][0]}" for t in CUT) + "\n\n"
           "Q1 dB = (pedestal30 - b_ref) / sigma_ref;  Q2 slope of the pulse-free samples (sigma_ref per 100 ns) and its significance;  Q3 rms of the pulse-free samples about that line / sigma_ref;\n"
           "Q4 tail = mean of the last 10 samples - pedestal30, in sigma_ref (nan if the last 10 samples are inside a fit window);  Q5 w50 = FWHM(raw) / FWHM(template at the main pulse's fit);\n"
           "Q6 rawsum = walk-window raw sum / (A_main I_m);  Q7 onset span (ns);  Q8 r = A_sec,max / A_main;  Q9 lag-1 autocorrelation of the normalised residual of the single-pulse refit;  Q10 its chi2/dof and ratio to the stored chi2/dof;\n"
           "Q11 nclip;  Q12 signed distance of the nearest secondary onset to the clipped stretch (ns, negative inside);  Q13 E_main / E_beam;  Q14 largest secondary E / E_beam.\n"
           f"Samples: seed {SEED}; per run type {POP_RAND} random kept PbWO4 channels, {POP_K2} with k >= 2, {POP_CLIP} with nclip > 0 (pools: " +
           ", ".join(f"{RT_LABEL[rt]} nclip>0 {summary[(rt, 'clip')]}" for rt in TYPES) + ").  "
           f"Refits not reproducing the stored k and chi2/dof: {nfail}.  Compute time {elapsed:.0f} s.")
    fig.text(0.02, 0.38, txt, fontsize=10.5, va="top", family="monospace", linespacing=1.5)
    pdf.savefig(fig)
    plt.close(fig)


def report(data, stem, sample_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    rows, exs = data["rows"], data["exs"]
    flag_all(rows)
    flag_all(exs)
    good = [r for r in rows if r.get("ok")]
    nfail = len(rows) - len(good) + sum(not e["ok"] for e in exs)
    acc = full_sample(sample_dir)
    with PdfPages(str(stem) + ".pdf") as pdf:
        ex_page(pdf, plt, exs)
        for ty in (1, 2, 3, 4):
            type_page(pdf, plt, ty, good, exs, acc)
        table_page(pdf, plt, exs, data["elapsed"], nfail, data["summary"])
    base = ["rt", "run", "ev", "name", "midx", "split", "pop", "ok", "k", "st_k", "chi2", "st_chi2", "amain"]
    cols = base + QCOLS + ["c1", "c2", "c3", "c4"]
    write_tsv(str(stem) + ".tsv", rows, cols)
    write_tsv(str(stem) + "_examples.tsv", exs, ["id"] + base[:-1] + ["refsrc"] + ["amain"] + QCOLS + ["c1", "c2", "c3", "c4"])
    # flag summary for the report
    print("[flags] random population:")
    for ty in CUT:
        for rt in TYPES:
            rr = [r for r in good if r["rt"] == rt and r["pop"] == "rand"]
            print(f"  cut{ty} {rt}: {sum(r[f'c{ty}'] for r in rr)}/{len(rr)}")
    print(f"[done] {len(rows)} sample rows, {len(exs)} examples, {nfail} not reproduced")
    for t, c in acc.items():
        print("[full]", t, c)


def main():
    global EVIO_DIR
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir", default=str(PYDIR / "output" / "pileup_rates" / "rawsum_sample"))
    ap.add_argument("--full-scan-dir", default=str(PYDIR / "output" / "pileup_rates" / "scan"))
    ap.add_argument("--evio-dir", default=EVIO_DIR)
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--cache", default=str(Path(tempfile.gettempdir()) / "exotic_p1_cache.pkl"))
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    EVIO_DIR = a.evio_dir
    if a.plot_only:
        data = pickle.load(open(a.cache, "rb"))
    else:
        data = compute(a)
        pickle.dump(data, open(a.cache, "wb"))
    report(data, Path(a.out) / STEM, a.sample_dir)


if __name__ == "__main__":
    main()
