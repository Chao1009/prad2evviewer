#!/usr/bin/env python3
"""exotic_p2_flags.py -- P2 of the exotic-waveform study: the four exotic flags of pileup_rates.exotic_vars on the P1 channels, and the BASE_SLOPE / CLIP_SINGLE fit options on the examples.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p2_flags.py [--sample-dir <dir>] [--full-scan-dir <dir>] [--evio-dir <dir>] [--out <dir>] [--cache <pkl>] [--jobs 3]
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p2_flags.py --plot-only --cache <pkl>        # redo the TSV / PDF from a cached compute stage

Channels: those of exotic_p1_quantities (its 12 EXAMPLES, and its samples of the 2-file rawsum_sample scans: per run type POP_RAND random, POP_K2 random with k >= 2, POP_CLIP
random with nclip > 0, seed SEED, PbWO4 kept channels; build_tasks is imported).  Every channel is refitted from the raw samples with the default fit (fit_event) and must reproduce the
stored k and chi2/dof; the flags are computed on a single-row pileup_rates.fit of the same channel (same b_ref, sigma_ref), which must reproduce the fit_event k and chi2/dof.
Flags, as implemented in pileup_rates.exotic_vars (default thresholds EX_DB 8, EX_TAIL -5, EX_W50 3, EX_STRETCH_K 3, EX_STRETCH_W50 1.5, EX_CLIP_NS 8 ns, EX_ESEC 0.5):
 f1 tail-baseline  ex_db > 8 and ex_tail < -5; ex_db = (pedestal30 - b_ref) / sigma_ref, ex_tail = (mean of the pulse-free samples of Y[90:100] - pedestal30) / sigma_ref (>= 8 of 10 pulse-free)
 f2 wide           w50 > 3 (w50 = FWHM(Y - fit baseline) / FWHM(template at the main pulse), pileup_rates.shape_vars)
 f3 smooth-wide    k >= 2 and ex_chi2s <= ex_chi2k + 2 (k - 1) ln N_M.  ex_chi2s = chi2 of the best single time-stretched-template pulse (s on 1.0 ... 3.0, 1 ns onset grid; stretch_fit) on the
                   sample set M of the stored chi2 (fit windows of the fitted pulses, no overflow), ex_chi2k = the k-pulse fit's chi2 on M, N_M = ex_nm = |M| (BIC: the k-pulse fit has
                   2 (k - 1) more parameters).  The stretch fit runs only for k >= 3 or w50 > 1.5 (else ex_nm = 0, flag off)
 f4 clipped-multi  nclip >= 1 and (ex_dclip <= 8 ns or ex_esec > 0.5); ex_dclip = min over the secondary pulses of max(c0 - t, t - c1), [c0, c1] = 4 ns x [first, last overflow sample] (negative inside),
                   ex_esec = largest secondary A x gain / E_beam
The flags are always computed with the default EX_CLIP_NS = 8 ns (EX_CLIP_NS is also the CLIP_SINGLE margin, changed only inside the CLIP_SINGLE refits).
Examples: the flags on the default fit (hits against the example's exotic type), then the fit options pileup_rates.BASE_SLOPE for 1a-1d and CLIP_SINGLE (margin 8 and 20 ns)
for 4a, 4b (all examples carry the three option refits in the TSV).  BASE_SLOPE acts only on the gated rows (30-sample pedestal > EX_DB = 8 sigma_ref above b_ref; mf.gate_rows after the fit),
which get a joint sloped-baseline fit over every unclipped sample; every other row is fitted as by default.  Every sample channel is also refitted under BASE_SLOPE (bs_gate, bs_k, ... in the TSV):
gated and k-changed counts per run type and population on the BASE_SLOPE page.  Clip population (nclip > 0 sample): dclip distribution of the k >= 2 channels, fractions within 8 ... 40 ns, and the
CLIP_SINGLE refit at 8 and 20 ns (channels changing k, channels ending at k = 1).  Flag fractions per run type and population with the share where the stretch fit ran.
Outputs (--out): exotic_P2_flags_20261117.pdf / .tsv (one row per channel, examples first; consistency columns p1_c1, p1_c2 from exotic_P1_quantities_20261117*.tsv).
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import csv
import time
import pickle
import argparse
import tempfile
import textwrap
import multiprocessing as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import exotic_p1_quantities as p1                       # noqa: E402
from exotic_p1_quantities import EXAMPLES, TYPES, RT_LABEL, build_tasks, get_ctx, SEED, fmt, write_tsv   # noqa: E402

STEM = "exotic_P2_flags_20261117"
P1_STEM = "exotic_P1_quantities_20261117"
EX_TYPE = {e[0]: e[1] for e in EXAMPLES}
NS, DT, OVER = 100, 4.0, 4095
OPTS = (("bs", dict(bs=True, cs=None)), ("cs8", dict(bs=False, cs=8.0)), ("cs20", dict(bs=False, cs=20.0)))
MARGINS = (8, 12, 16, 20, 30, 40)
FLAGS = ("f1", "f2", "f3", "f4")
EXK = ["ex_db", "ex_tail", "w50", "ex_s", "ex_chi2s", "ex_chi2k", "ex_nm", "ex_dclip", "ex_esec"]
POPCOL = {"rand": "0.45", "k2": "tab:blue", "clip": "tab:orange"}
TCOL = dict(zip(TYPES, ["tab:blue", "tab:green", "tab:purple", "tab:orange"]))


# ---------------------------------------------------------------------------------------------------------------- per-channel computation
def fit_opt(ctx, Y, m, sg, br, bs=False, cs=None):
    """pileup_rates.fit of the single row Y [1, ns] with the options BASE_SLOPE (bs) and CLIP_SINGLE at margin cs (ns; None: off); the module globals are restored."""
    import pileup_rates as pr
    old = (pr.BASE_SLOPE, pr.CLIP_SINGLE, pr.EX_CLIP_NS)
    try:
        pr.BASE_SLOPE = bool(bs)
        pr.CLIP_SINGLE = cs is not None
        if cs is not None:
            pr.EX_CLIP_NS = float(cs)
        f = pr.fit(ctx.mf, Y, m, sg, br)
        f["gate"] = bool(np.asarray(ctx.mf.gate_rows)[0]) if bs else False       # BASE_SLOPE: row fitted with the sloped baseline (pedestal > EX_DB sigma_ref above b_ref)
        return f
    finally:
        pr.BASE_SLOPE, pr.CLIP_SINGLE, pr.EX_CLIP_NS = old


def flags_of(ctx, f, Y, m, sg, br, mi):
    """exotic_vars of the single row (default EX_CLIP_NS) for the fit result f; w50 on the fit baseline as quality_vars."""
    import pileup_rates as pr
    assert pr.EX_CLIP_NS == 8.0
    w50, _ = pr.shape_vars(Y, pr.fit_baseline(f, Y.shape[1]) if "bslope" in f else f["ped"], m, f, ctx.bank)
    ex = pr.exotic_vars(ctx.mf, Y, m, f, br, sg, ctx.gain[[mi]], ctx.beam, w50)
    out = {c: float(ex[c][0]) for c in ("ex_db", "ex_tail", "ex_s", "ex_chi2s", "ex_chi2k", "ex_dclip", "ex_esec")}
    out["ex_nm"] = int(ex["ex_nm"][0])
    out["w50"] = float(w50[0])
    for i in (1, 2, 3, 4):
        out[f"f{i}"] = bool(ex[f"ex_f{i}"][0])
    return out


def summary_of(ctx, f, mi):
    k = int(f["k"][0])
    a = f["a"][0][:k].astype(float)
    return dict(k=k, chi2=float(f["chi2"][0]), amain=float(a.max()) if k else np.nan, a=a, t=f["t"][0][:k].astype(float), E=a * float(ctx.gain[mi]),
                bslope=float(f["bslope"][0]) if "bslope" in f else np.nan, gate=bool(f.get("gate", False)))


def channel(ctx, r, row, b_ref, sig_ref, full_opts, clip_opts):
    """Record of channel row of the fit_event result r: default fit (single-row refit), flags; the option refits for examples (full_opts) or the clip population (clip_opts)."""
    import pileup_rates as pr
    mi = int(r["mi"][row])
    Y = r["Y"][row][None]
    m = ctx.bank.row[[mi]]
    sg, br = sig_ref[[mi]], b_ref[[mi]]
    f = fit_opt(ctx, Y, m, sg, br)
    rec = summary_of(ctx, f, mi)
    rec["ok1"] = rec["k"] == int(r["k"][row]) and bool(np.isclose(rec["chi2"], r["chi2"][row], rtol=1e-6, atol=1e-9))
    rec["nclip"] = int(pr.clipped(Y)[0])
    rec.update(flags_of(ctx, f, Y, m, sg, br, mi))
    opts = {}
    if True:
        for name, kw in OPTS:
            if not full_opts and kw["cs"] is not None and not clip_opts:
                continue
            fo = fit_opt(ctx, Y, m, sg, br, **kw)
            s = summary_of(ctx, fo, mi)
            if full_opts:
                s.update(flags_of(ctx, fo, Y, m, sg, br, mi))
                s["model"] = fo["comp"][0][:s["k"]].astype(float)
                s["base"] = pr.fit_baseline(fo, NS)[0]
                s["stretch"] = stretch_model(ctx, s, fo, Y, m, sg)
            opts[name] = s
    rec["opts"] = opts
    if full_opts:
        rec["model"] = f["comp"][0][:rec["k"]].astype(float)
        rec["base"] = pr.fit_baseline(f, NS)[0]
        rec["stretch"] = stretch_model(ctx, rec, f, Y, m, sg)
    return rec


def stretch_model(ctx, s, f, Y, m, sg):
    import pileup_rates as pr
    if s["k"] < 1 or not (s["k"] >= pr.EX_STRETCH_K or s.get("w50", 0) > pr.EX_STRETCH_W50):
        return None
    mod = pr.stretch_fit(ctx.mf, Y, m, f, sg, rows=[0])["model"][0]
    return mod if np.isfinite(mod).all() else None


# ---------------------------------------------------------------------------------------------------------------- worker
def work(task):
    """('rand', run, files, picks) or ('ex', run, evs, full-scan npz): records as in exotic_p1_quantities.work, with flags (and option refits for the examples and the clip population)."""
    import event_display as ed
    from prad2py import dec
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
                        rec = dict(run=run, rt=ctx.rt, ev=en, midx=midx, name=str(ctx.names[midx]), pop=st["pop"], id="", ok=False, ok1=False)
                        if row is not None:
                            k, chi2 = int(r["k"][row]), float(r["chi2"][row])
                            rec["ok"] = k == st["k"] and bool(np.isclose(chi2, st["chi2"], rtol=2e-3, atol=1e-3))
                            rec.update(st_k=st["k"], st_chi2=st["chi2"])
                            rec.update(channel(ctx, r, row, b_ref, sig_ref, False, st["pop"] == "clip"))
                        out.append(rec)
            ch.close()
        for en in left:
            for midx, st in picks[en].items():
                out.append(dict(run=run, rt=ctx.rt, ev=en, midx=midx, name=str(ctx.names[midx]), pop=st["pop"], id="", ok=False, ok1=False, missing=True))
    else:
        evs, scan_path = task[2], task[3]
        z = np.load(scan_path)
        splits = [int(s) for s in z["splits"]]
        evnum = z["ev_number"]
        for ident, en, mname in evs:
            midx = int(ctx.index[mname])
            f, si, files = ed.locate(ctx, p1.EVIO_DIR, en)
            spl = int(f.rsplit(".", 1)[1])
            if spl in splits:
                b_ref, sig_ref = z["mod_b_ref"][splits.index(spl)], z["mod_sig_ref"][splits.index(spl)]
            else:
                _, b_ref, sig_ref = ed.noise_reference(ctx, f)
            fe = ed.find_event(ctx, f, en)
            rec = dict(id=ident, run=run, rt=ctx.rt, ev=en, midx=midx, name=mname, pop="example", ok=False, ok1=False)
            r = ed.fit_event(ctx, fe, b_ref, sig_ref) if fe is not None else None
            row = r["row"].get(midx) if r is not None else None
            if row is None:
                out.append(rec)
                continue
            k, chi2 = int(r["k"][row]), float(r["chi2"][row])
            ie = np.nonzero(evnum == en)[0]
            st = None
            if ie.size:
                sel = np.nonzero((z["ch_ev"] == ie[0]) & (z["ch_midx"] == midx))[0]
                if sel.size:
                    st = (int(z["ch_k"][sel[0]]), float(z["ch_chi2"][sel[0]]))
            rec["st_k"], rec["st_chi2"] = st if st else (None, None)
            rec["ok"] = st is not None and k == st[0] and bool(np.isclose(chi2, st[1], rtol=2e-3, atol=1e-3))
            rec.update(channel(ctx, r, row, b_ref, sig_ref, True, False))
            rec["Y"] = r["Y"][row].astype(float)
            rec["gain"], rec["beam"] = float(ctx.gain[midx]), float(ctx.beam)
            out.append(rec)
    print(f"[work] {kind} run {run}: {len(out)} records, {time.time() - t_start:.0f} s", flush=True)
    return out


def compute(a):
    t0 = time.time()
    p1.EVIO_DIR = a.evio_dir
    tasks, ex, summary = build_tasks(a.sample_dir, a.full_scan_dir)
    print(f"[sample] pools {summary}", flush=True)
    allt = sorted(tasks, key=lambda x: -sum(len(v) for v in x[3].values())) + ex
    with mp.get_context("fork").Pool(a.jobs) as pool:
        res = pool.map(work, allt, chunksize=1)
    rows = [r for rr, t in zip(res, allt) if t[0] == "rand" for r in rr]
    exs = [r for rr, t in zip(res, allt) if t[0] == "ex" for r in rr]
    exs.sort(key=lambda r: [e[0] for e in EXAMPLES].index(r["id"]))
    return dict(rows=rows, exs=exs, summary=summary, elapsed=time.time() - t0)


# ---------------------------------------------------------------------------------------------------------------- analysis helpers
def has(r):
    return "k" in r and r.get("ok") and r.get("ok1")


def fl(r, i):
    return bool(r.get(f"f{i}", False))


def anyflag(r):
    return any(fl(r, i) for i in (1, 2, 3, 4))


def read_p1():
    """P1 cut1 / cut2 by (run, ev, midx) and by example id; None when the P1 files are missing."""
    try:
        by_key, by_id = {}, {}
        with open(HERE / f"{P1_STEM}.tsv") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                by_key[(int(r["run"]), int(r["ev"]), int(r["midx"]))] = (int(r["c1"]), int(r["c2"]))
        with open(HERE / f"{P1_STEM}_examples.tsv") as f:
            for r in csv.DictReader(f, delimiter="\t"):
                by_id[r["id"]] = (int(r["c1"]), int(r["c2"]))
        return by_key, by_id
    except OSError:
        return None


def consistency(rows, exs):
    p = read_p1()
    mism, n = [], 0
    for r in rows + exs:
        r["p1_c1"] = r["p1_c2"] = ""
    if p is None:
        return mism, n
    for r in rows + exs:
        v = p[1].get(r["id"]) if r["pop"] == "example" else p[0].get((r["run"], r["ev"], r["midx"]))
        if v is None or "k" not in r:
            continue
        r["p1_c1"], r["p1_c2"] = v
        n += 1
        if (int(fl(r, 1)), int(fl(r, 2))) != v:
            mism.append(r)
    return mism, n


def popstats(rows, rt, pop):
    rr = [r for r in rows if r["rt"] == rt and r["pop"] == pop and has(r)]
    n = len(rr)
    out = dict(n=n, stretch=sum(r["ex_nm"] > 0 for r in rr) / n if n else np.nan)
    for i in (1, 2, 3, 4):
        out[f"f{i}"] = sum(fl(r, i) for r in rr) / n if n else np.nan
    out["any"] = sum(anyflag(r) for r in rr) / n if n else np.nan
    return out


def clip_stats(rows):
    """Per run type: dclip distribution of the clip-population channels with nclip > 0 and a secondary pulse, and the CLIP_SINGLE refit counts."""
    out = {}
    for rt in TYPES:
        cl = [r for r in rows if r["rt"] == rt and r["pop"] == "clip" and has(r)]
        sec = [r for r in cl if r["k"] >= 2 and np.isfinite(r["ex_dclip"])]
        d = np.array([r["ex_dclip"] for r in sec])
        low = np.array([r["ex_esec"] <= 0.5 for r in sec], bool)
        s = dict(n=len(cl), nsec=len(sec), d=d, low=low, nlow=int(low.sum()))
        s["q"] = np.percentile(d, [10, 25, 50, 75, 90]) if d.size else np.full(5, np.nan)
        s["fr"] = {M: float((d <= M).mean()) if d.size else np.nan for M in MARGINS}
        s["fr_low"] = {M: float((d[low] <= M).mean()) if low.any() else np.nan for M in MARGINS}
        s["cnt_low"] = {M: int((d[low] <= M).sum()) for M in MARGINS}
        for name in ("cs8", "cs20"):
            ok = [r for r in cl if name in r["opts"]]
            s[name] = dict(chg=sum(r["opts"][name]["k"] != r["k"] for r in ok), k1=sum(r["opts"][name]["k"] == 1 for r in ok), n=len(ok),
                           k1_def=sum(r["k"] == 1 for r in ok), chg_sec=sum(r["opts"][name]["k"] != r["k"] for r in ok if r["k"] >= 2))
        out[rt] = s
    return out


def binom(p, n):
    return np.sqrt(max(p * (1 - p), 0) / n) if n else np.nan


# ---------------------------------------------------------------------------------------------------------------- plots
def draw(ax, axr, Y, base, comp, stretch, title, tcol="0.0"):
    tt = DT * np.arange(NS)
    ov = Y >= OVER
    ax.plot(tt[~ov], Y[~ov], "k.", ms=5)
    if ov.any():
        ax.plot(tt[ov], Y[ov], "x", color="tab:red", ms=7, mew=1.6)
    ax.plot(tt, base, color="0.6", lw=0.9)
    tot = base + comp.sum(0)
    ax.plot(tt, tot, "r-", lw=1.4)
    for c, cp in zip(("tab:blue", "tab:orange", "tab:green", "tab:purple", "tab:brown", "tab:pink"), comp):
        ax.plot(tt, base + cp, "--", color=c, lw=0.9)
    if stretch is not None:
        ax.plot(tt[~ov], stretch[~ov], ":", color="darkgreen", lw=1.8)
    yl = (min(Y.min(), tot.min()), Y.max())
    ax.set_ylim(yl[0] - 0.05 * (yl[1] - yl[0]), min(yl[1] * 1.08 + 5, 4500))
    ax.set_xlim(0, DT * NS)
    ax.tick_params(labelsize=9, labelbottom=False)
    ax.set_title(title, fontsize=10.5, loc="left", pad=4, color=tcol)
    if axr is not None:
        axr.plot(tt[~ov], (Y - tot)[~ov], "k.-", lw=0.6, ms=3)
        axr.axhline(0, color="r", lw=0.8)
        axr.set_xlim(0, DT * NS)
        axr.tick_params(labelsize=9)
        axr.set_xlabel("t (ns)", fontsize=9)


def flag_text(e):
    """Text lines of the four flags of record e (default fit)."""
    import pileup_rates as pr
    mk = lambda i: "Y" if fl(e, i) else "."                  # noqa: E731
    n = lambda v, s: f"{v:{s}}" if np.isfinite(v) else "nan"  # noqa: E731
    if e["ex_nm"] > 0:
        bic = e["ex_chi2k"] + 2 * (e["k"] - 1) * np.log(e["ex_nm"])
        f3 = f"chi2_s {e['ex_chi2s']:.1f} (s {e['ex_s']:.1f}) vs chi2_k + 2(k-1)ln N_M = {e['ex_chi2k']:.1f} + {2 * (e['k'] - 1) * np.log(e['ex_nm']):.1f} = {bic:.1f}"
    else:
        f3 = "stretch fit not run (k < 3 and w50 <= 1.5)"
    return (f"f1 {mk(1)}  db {n(e['ex_db'], '.1f')} (> {pr.EX_DB:g}), tail {n(e['ex_tail'], '.1f')} (< {pr.EX_TAIL:g})\n"
            f"f2 {mk(2)}  w50 {n(e['w50'], '.2f')} (> {pr.EX_W50:g})\n"
            f"f3 {mk(3)}  {f3}\n"
            f"f4 {mk(4)}  nclip {e['nclip']}, dclip {n(e['ex_dclip'], '.1f')} (<= {pr.EX_CLIP_NS:g}), esec {n(e['ex_esec'], '.3f')} (> {pr.EX_ESEC:g})")


def ex_page(pdf, plt, exs):
    fig = plt.figure(figsize=(44, 18))
    outer = fig.add_gridspec(2, 6, hspace=0.3, wspace=0.12, left=0.03, right=0.995, top=0.86, bottom=0.02)
    fig.text(0.03, 0.975, "P2: the 12 examples with the four exotic flags (default fit, pileup_rates.exotic_vars)", fontsize=22, fontweight="bold")
    fig.text(0.03, 0.945, "black dots: raw samples (red x: overflow); red: total fit; dashed colours: fitted pulses; dotted green: best single stretched-template pulse (where the stretch fit ran); grey: fit baseline.  "
             "Text: f1-f4 (Y = flag set) with their values; the title colour is green when the flag of the example's type is set, red when it is not.", fontsize=12)
    for j, e in enumerate(exs):
        gs = outer[j // 6, j % 6].subgridspec(3, 1, height_ratios=[3, 1, 1.6], hspace=0.25)
        ax, axr, axt = fig.add_subplot(gs[0]), fig.add_subplot(gs[1]), fig.add_subplot(gs[2])
        hit = fl(e, EX_TYPE[e["id"]])
        draw(ax, axr, e["Y"], e["base"], e["model"], e["stretch"],
             f"{e['id']} (type {EX_TYPE[e['id']]}) run {e['run']} ev {e['ev']} {e['name']} {e['rt']}\nk = {e['k']}  chi2/dof = {e['chi2']:.2f}  stored {e['st_k']} / {e['st_chi2']:.2f}  "
             f"{'OK' if has(e) else 'NOT REPRODUCED'}   E_main {e['amain'] * e['gain']:.0f} MeV   f{EX_TYPE[e['id']]} {'HIT' if hit else 'MISS'}", "darkgreen" if hit else "red")
        axt.axis("off")
        axt.text(0.0, 1.0, flag_text(e), transform=axt.transAxes, fontsize=9.3, va="top", family="monospace", linespacing=1.5)
    pdf.savefig(fig)
    plt.close(fig)


def table_page(pdf, plt, exs, mism, ncons, nfail, nfail1, elapsed):
    fig = plt.figure(figsize=(30, 14))
    fig.text(0.02, 0.96, "Flags of the 12 examples (default fit)", fontsize=20, fontweight="bold")
    hdr = ["id", "type", "run", "event", "module", "k", "chi2/dof", "f1", "f2", "f3", "f4", "db", "tail", "w50", "s", "chi2_s", "chi2_k", "N_M", "BIC pen.", "dclip", "esec", "hit f(type)", "extra flags"]
    cells = []
    for e in exs:
        ty = EX_TYPE[e["id"]]
        pen = 2 * (e["k"] - 1) * np.log(e["ex_nm"]) if e["ex_nm"] > 0 else np.nan
        g = lambda v, s=".3g": "nan" if not np.isfinite(v) else f"{v:{s}}"      # noqa: E731
        cells.append([e["id"], str(ty), str(e["run"]), str(e["ev"]), e["name"], str(e["k"]), f"{e['chi2']:.2f}"] + ["Y" if fl(e, i) else "." for i in (1, 2, 3, 4)] +
                     [g(e["ex_db"]), g(e["ex_tail"]), g(e["w50"]), g(e["ex_s"], ".1f"), g(e["ex_chi2s"]), g(e["ex_chi2k"]), str(e["ex_nm"]), g(pen), g(e["ex_dclip"]), g(e["ex_esec"]),
                      "HIT" if fl(e, ty) else "MISS", ",".join(f"f{i}" for i in (1, 2, 3, 4) if fl(e, i) and i != ty) or "-"])
    ax = fig.add_axes([0.01, 0.5, 0.98, 0.43])
    ax.axis("off")
    tb = ax.table(cellText=cells, colLabels=hdr, loc="upper center", cellLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(8.5)
    tb.scale(1, 1.8)
    hits = sum(fl(e, EX_TYPE[e["id"]]) for e in exs)
    miss = [e["id"] for e in exs if not fl(e, EX_TYPE[e["id"]])]
    txt = (f"Flags: f1 db > 8 and tail < -5;  f2 w50 > 3;  f3 k >= 2 and chi2_s <= chi2_k + 2(k-1) ln N_M;  f4 nclip >= 1 and (dclip <= 8 ns or esec > 0.5).\n"
           f"Hits on the example's own type: {hits} of {len(exs)}; misses: {', '.join(miss) or 'none'}.\n"
           f"f3 set on: {', '.join(e['id'] for e in exs if fl(e, 3)) or 'none'}; not set on: {', '.join(e['id'] for e in exs if not fl(e, 3)) or 'none'} "
           f"(type-3 examples with f3 set: {sum(fl(e, 3) for e in exs if EX_TYPE[e['id']] == 3)} of 3).\n"
           f"Consistency with P1: flag 1 / flag 2 against the P1 cut1 / cut2 on {ncons} channels: {len(mism)} mismatches.\n"
           f"Refits not reproducing the stored k and chi2/dof (fit_event): {nfail}; single-row pileup_rates.fit not reproducing the fit_event k and chi2/dof: {nfail1}.  Compute time {elapsed:.0f} s.")
    fig.text(0.02, 0.44, txt, fontsize=11.5, va="top", family="monospace", linespacing=1.6)
    pdf.savefig(fig)
    plt.close(fig)


def opt_table(ax, hdr, cells, size=9):
    ax.axis("off")
    tb = ax.table(cellText=cells, colLabels=hdr, loc="upper center", cellLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(size)
    tb.scale(1, 1.55)


def bs_stats(rows):
    """BASE_SLOPE on the sample rows: per (run type, pop) the channels, gated ones, and gated ones whose k changes."""
    out = {}
    for rt in TYPES:
        for pop in ("rand", "k2", "clip"):
            rr = [r for r in rows if r["rt"] == rt and r["pop"] == pop and has(r) and "bs" in r["opts"]]
            g = [r for r in rr if r["opts"]["bs"]["gate"]]
            ch = [r for r in rr if r["opts"]["bs"]["k"] != r["k"]]
            out[(rt, pop)] = dict(n=len(rr), gated=len(g), changed=len(ch), changed_gated=sum(r["opts"]["bs"]["k"] != r["k"] for r in g),
                                  to1=sum(r["opts"]["bs"]["k"] <= 1 for r in g), kdef1=sum(r["k"] <= 1 for r in g))
    return out


def base_slope_page(pdf, plt, exs, bss):
    E = [e for e in exs if EX_TYPE[e["id"]] == 1]
    fig = plt.figure(figsize=(28, 20))
    gs = fig.add_gridspec(4, len(E), height_ratios=[1, 1, 0.75, 1.0], hspace=0.25, wspace=0.12, left=0.04, right=0.99, top=0.88, bottom=0.02)
    fig.text(0.03, 0.965, "BASE_SLOPE before / after on the type-1 examples (top: default fit, middle: pileup_rates.BASE_SLOPE = True)", fontsize=20, fontweight="bold")
    fig.text(0.03, 0.935, "BASE_SLOPE acts only on gated rows (30-sample pedestal > 8 sigma_ref above b_ref): joint sloped-baseline fit over every unclipped sample, no BASE_CLIP, dof one less; other rows are fitted as by default.\n"
             "red: total fit; dashed: pulses; grey: fit baseline (offset + slope t / 100 ns); dotted green: stretched template.", fontsize=12)
    cells = []
    for j, e in enumerate(E):
        for i, (name, s) in enumerate((("default", e), ("BASE_SLOPE", e["opts"]["bs"]))):
            ax = fig.add_subplot(gs[i, j])
            draw(ax, None, e["Y"], s["base"], s["model"], s["stretch"], f"{e['id']} {name}: k = {s['k']}  chi2/dof = {s['chi2']:.2f}  bslope {s['bslope']:.2f}  gated {'yes' if s['gate'] else 'no'}" if name != "default" else
                 f"{e['id']} {name}: k = {s['k']}  chi2/dof = {s['chi2']:.2f}")
            ax.tick_params(labelbottom=(i == 1))
            if i == 1:
                ax.set_xlabel("t (ns)")
        for name, s in (("default", e), ("BASE_SLOPE", e["opts"]["bs"])):
            cells.append([e["id"], name, "yes" if (name != "default" and s["gate"]) else "-", str(s["k"]), f"{s['chi2']:.2f}", ", ".join(f"{x:.0f}" for x in s["a"]), ", ".join(f"{x:.0f}" for x in s["E"]),
                          "-" if name == "default" else f"{s['bslope']:.2f}"] + ["Y" if s[f"f{i}"] else "." for i in (1, 2, 3, 4)] + [f"{s['ex_tail']:.1f}", f"{s['w50']:.2f}"])
    opt_table(fig.add_subplot(gs[2, :]), ["id", "fit", "gated", "k", "chi2/dof", "amplitudes (ADC)", "E per pulse (MeV)", "bslope (ADC/100 ns)", "f1", "f2", "f3", "f4", "tail (sigma)", "w50"], cells, 10)
    lines = [f"BASE_SLOPE on the sample channels (kept PbWO4, default k and chi2 reproduced): gated rows and rows whose k changes (only gated rows can change)",
             f"{'run type':18s}{'pop':>6s}{'N':>6s}{'gated':>7s}{'gated %':>9s}{'k changed':>10s}{'gated k<=1 (default / BASE_SLOPE)':>38s}"]
    for rt in TYPES:
        for pop in ("rand", "k2", "clip"):
            c = bss[(rt, pop)]
            lines.append(f"{RT_LABEL[rt] if pop == 'rand' else '':18s}{pop:>6s}{c['n']:6d}{c['gated']:7d}{100 * c['gated'] / max(c['n'], 1):9.2f}{c['changed']:10d}{c['kdef1']:>20d} /{c['to1']:3d}")
    axt = fig.add_subplot(gs[3, :])
    axt.axis("off")
    axt.text(0.0, 1.0, "\n".join(lines), transform=axt.transAxes, fontsize=10, va="top", family="monospace", linespacing=1.35)
    pdf.savefig(fig)
    plt.close(fig)


def clip_single_page(pdf, plt, exs):
    E = [e for e in exs if EX_TYPE[e["id"]] == 4]
    fig = plt.figure(figsize=(28, 14))
    gs = fig.add_gridspec(2, 3, height_ratios=[2, 1.2], hspace=0.3, wspace=0.1, left=0.04, right=0.99, top=0.88, bottom=0.03)
    fig.text(0.03, 0.965, "CLIP_SINGLE before / after on the type-4 examples (default, margin 8 ns, margin 20 ns)", fontsize=20, fontweight="bold")
    fig.text(0.03, 0.935, "red: total fit; dashed: pulses; grey: baseline.  CLIP_SINGLE forbids secondary / shoulder onsets within the margin of the overflow stretch [first, last overflow sample] x 4 ns (red x: overflow samples).", fontsize=12)
    cells = []
    sub = gs[0, :].subgridspec(len(E), 3, hspace=0.35, wspace=0.12)
    for j, e in enumerate(E):
        for i, (name, s) in enumerate((("default", e), ("CLIP_SINGLE 8 ns", e["opts"]["cs8"]), ("CLIP_SINGLE 20 ns", e["opts"]["cs20"]))):
            ax = fig.add_subplot(sub[j, i])
            draw(ax, None, e["Y"], s["base"], s["model"], None, f"{e['id']} {name}: k = {s['k']}  chi2/dof = {s['chi2']:.2f}   A = " + ", ".join(f"{x:.0f}" for x in s["a"]))
            ax.tick_params(labelbottom=(j == len(E) - 1))
            cells.append([e["id"], name, str(s["k"]), f"{s['chi2']:.2f}", ", ".join(f"{x:.0f}" for x in s["a"]), ", ".join(f"{x:.0f}" for x in s["E"]),
                          ", ".join(f"{x:.0f}" for x in s["t"]), "Y" if s["f4"] else ".", f"{s['ex_dclip']:.1f}", f"{s['ex_esec']:.3f}"])
    opt_table(fig.add_subplot(gs[1, :]), ["id", "fit", "k", "chi2/dof", "A (ADC, time order)", "E (MeV)", "onsets (ns)", "f4", "dclip (ns)", "esec"], cells, 10)
    pdf.savefig(fig)
    plt.close(fig)


def clip_recommendation(cs):
    """One line: the recommended CLIP_SINGLE margin and the numbers behind it (2.24 and 3.49 GeV H2 clip populations)."""
    big = [rt for rt in ("2p24GeV_H2", "3p49GeV_H2")]
    chg = {nm: sum(cs[rt][nm]["chg"] for rt in big) for nm in ("cs8", "cs20")}
    k1 = {nm: sum(cs[rt][nm]["k1"] - cs[rt][nm]["k1_def"] for rt in big) for nm in ("cs8", "cs20")}
    n = sum(cs[rt]["n"] for rt in big)
    return (f"Recommended CLIP_SINGLE margin: keep 8 ns.  On the 2.24 + 3.49 GeV H2 clip populations ({n} channels) 8 ns changes k in {chg['cs8']} (+{k1['cs8']} new k = 1), 20 ns in {chg['cs20']} "
            f"(+{k1['cs20']} new k = 1); dclip has no gap ({', '.join(f'{100 * np.mean(np.concatenate([cs[rt]['d'] for rt in big]) <= M):.0f}% <= {M} ns' for M in (8, 20, 40))}), "
            f"so a wider margin removes more plausibly real secondaries, and neither margin lowers k on 4a / 4b.")


def clip_margin_page(pdf, plt, cs):
    fig = plt.figure(figsize=(22, 14))
    gs = fig.add_gridspec(2, 1, height_ratios=[1, 1.25], hspace=0.2, left=0.06, right=0.985, top=0.9, bottom=0.03)
    fig.text(0.02, 0.965, "CLIP_SINGLE margin: dclip of clipped channels with a secondary pulse, and the refit on the clip population", fontsize=20, fontweight="bold")
    fig.text(0.02, 0.94, "Population: the P1 random nclip > 0 channels (up to 100 per run type), kept PbWO4, default fit k >= 2.  dclip = signed distance (ns) of the nearest secondary onset to the clipped stretch (negative inside).", fontsize=11.5)
    ax = fig.add_subplot(gs[0, 0])
    for rt in TYPES:
        d = np.sort(cs[rt]["d"])
        if d.size:
            ax.step(d, np.arange(1, d.size + 1) / d.size, where="post", color=TCOL[rt], lw=1.8, label=f"{RT_LABEL[rt]} (n={d.size})")
    for M in MARGINS:
        ax.axvline(M, color="0.7", lw=0.8)
    ax.set_xlim(-60, 120)
    ax.set_xlabel("dclip (ns)")
    ax.set_ylabel("cumulative fraction")
    ax.legend(fontsize=10, loc="lower right")
    axt = fig.add_subplot(gs[1, 0])
    axt.axis("off")
    lines = [f"{'run type':18s}{'n clip':>7s}{'k>=2':>6s}  dclip quantiles 10/25/50/75/90 (ns)   " + "".join(f"{'<=' + str(M):>7s}" for M in MARGINS)]
    for rt in TYPES:
        s = cs[rt]
        lines.append(f"{RT_LABEL[rt]:18s}{s['n']:7d}{s['nsec']:6d}  " + " ".join(f"{x:7.1f}" for x in s["q"]) + "   " + "".join(f"{100 * s['fr'][M]:6.0f}%" for M in MARGINS))
    lines += ["", "of which esec <= 0.5 (not caught by the energy condition of f4): n, and fraction within the margin (count)"]
    for rt in TYPES:
        s = cs[rt]
        lines.append(f"{RT_LABEL[rt]:18s}{s['nlow']:7d}" + " " * 52 + "".join(f"{100 * s['fr_low'][M]:4.0f}%({s['cnt_low'][M]:2d})" for M in MARGINS))
    lines += ["", "CLIP_SINGLE refit of every channel of the clip population (all k, stored default k in the first column):",
              f"{'run type':18s}{'n':>5s}{'default k=1':>12s}{'k>=2':>6s}  {'margin':>6s}{'k changed':>10s}{'of k>=2':>8s}{'ends k=1':>9s}"]
    for rt in TYPES:
        s = cs[rt]
        for nm, M in (("cs8", 8), ("cs20", 20)):
            c = s[nm]
            lines.append(f"{RT_LABEL[rt] if M == 8 else '':18s}{c['n']:5d}{c['k1_def']:12d}{c['n'] - c['k1_def']:6d}  {M:4d}ns{c['chg']:10d}{c['chg_sec']:8d}{c['k1']:9d}")
    lines += [""] + textwrap.wrap(clip_recommendation(cs), 190)
    axt.text(0.0, 1.0, "\n".join(lines), transform=axt.transAxes, fontsize=10.5, va="top", family="monospace", linespacing=1.45)
    pdf.savefig(fig)
    plt.close(fig)


def fraction_page(pdf, plt, rows, exs, summary):
    fig = plt.figure(figsize=(22, 12))
    fig.text(0.02, 0.96, "Flag fractions per run type and population (default fit, PbWO4)", fontsize=20, fontweight="bold")
    lines = [f"{'run type':18s}{'pop':>6s}{'N':>6s}{'f1 %':>12s}{'f2 %':>12s}{'f3 %':>12s}{'f4 %':>12s}{'any %':>12s}{'stretch ran %':>15s}"]
    for rt in TYPES:
        for pop in ("rand", "k2", "clip"):
            s = popstats(rows, rt, pop)
            n = s["n"]
            def c(key):
                if not n:
                    return "-"
                p = s[key]
                if p == 0:
                    return f"<{300 / n:.2f}"                    # 95% CL upper limit for zero counts (3 / N)
                return f"{100 * p:.2f}+-{100 * binom(p, n):.2f}" if pop == "rand" else f"{100 * p:.1f}"
            lines.append(f"{RT_LABEL[rt] if pop == 'rand' else '':18s}{pop:>6s}{n:6d}" + "".join(f"{c(k):>12s}" for k in ("f1", "f2", "f3", "f4", "any")) + f"{100 * s['stretch']:15.1f}")
    R = [r for r in rows if r["pop"] == "rand" and has(r)]
    lines += ["", f"random population, all run types: N = {len(R)}; " + ", ".join(f"f{i} {sum(fl(r, i) for r in R)}" for i in (1, 2, 3, 4)) + f", any {sum(anyflag(r) for r in R)}",
              "rand: fraction +- binomial error (percent); zero counts: <3/N (95% CL upper limit, percent); k2 / clip: fraction only (selected populations, not rates).  stretch ran = share with ex_nm > 0 (k >= 3 or w50 > 1.5).",
              "", "random PbWO4 sample sizes (pools): " + ", ".join(f"{RT_LABEL[rt]} k>=2 {summary[(rt, 'k2')]}, nclip>0 {summary[(rt, 'clip')]}" for rt in TYPES)]
    fig.text(0.02, 0.9, "\n".join(lines), fontsize=11.5, va="top", family="monospace", linespacing=1.55)
    pdf.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------- output
def tsv_rows(rows, exs):
    out = []
    for r in exs + rows:
        d = dict(r)
        for name in ("bs", "cs8", "cs20"):
            o = r.get("opts", {}).get(name)
            if o:
                d.update({f"{name}_k": o["k"], f"{name}_chi2": o["chi2"], f"{name}_amain": o["amain"], f"{name}_Emain": o["amain"] * float(o["E"][0] / o["a"][0]) if o["k"] else np.nan})
                if name == "bs":
                    d["bs_bslope"], d["bs_gate"] = o["bslope"], o["gate"]
                for i in (1, 2, 3, 4):
                    if f"f{i}" in o:
                        d[f"{name}_f{i}"] = o[f"f{i}"]
        d["any"] = anyflag(r) if "k" in r else ""
        d["stretch_ran"] = (r["ex_nm"] > 0) if "k" in r else ""
        out.append(d)
    return out


def report(data, stem, sample_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    rows, exs = data["rows"], data["exs"]
    mism, ncons = consistency(rows, exs)
    nfail = sum(not r.get("ok") for r in rows + exs)
    nfail1 = sum(("k" in r) and not r.get("ok1") for r in rows + exs)
    cs = clip_stats(rows)
    bss = bs_stats(rows)
    with PdfPages(str(stem) + ".pdf") as pdf:
        ex_page(pdf, plt, exs)
        table_page(pdf, plt, exs, mism, ncons, nfail, nfail1, data["elapsed"])
        base_slope_page(pdf, plt, exs, bss)
        clip_single_page(pdf, plt, exs)
        clip_margin_page(pdf, plt, cs)
        fraction_page(pdf, plt, rows, exs, data["summary"])
    base = ["rt", "run", "ev", "name", "midx", "pop", "id", "ok", "ok1", "k", "st_k", "chi2", "st_chi2", "amain", "nclip"] + EXK + ["f1", "f2", "f3", "f4", "any", "stretch_ran", "p1_c1", "p1_c2"]
    opt = [f"{n}_{q}" for n in ("bs", "cs8", "cs20") for q in ("k", "chi2", "amain", "Emain")] + ["bs_bslope", "bs_gate"] + [f"{n}_f{i}" for n in ("bs", "cs8", "cs20") for i in (1, 2, 3, 4)]
    write_tsv(str(stem) + ".tsv", tsv_rows(rows, exs), base + opt)
    print(f"[flags] examples:")
    for e in exs:
        print(f"  {e['id']} k {e['k']} chi2 {e['chi2']:.2f} " + " ".join(f"f{i}={int(fl(e, i))}" for i in (1, 2, 3, 4)) +
              f" db {e['ex_db']:.1f} tail {e['ex_tail']:.1f} w50 {e['w50']:.2f} s {e['ex_s']:.1f} chi2s {e['ex_chi2s']:.1f} chi2k {e['ex_chi2k']:.1f} nm {e['ex_nm']} dclip {e['ex_dclip']:.1f} esec {e['ex_esec']:.3f}")
    print("[options] examples:")
    for e in exs:
        for n, o in e["opts"].items():
            print(f"  {e['id']} {n}: k {o['k']} chi2 {o['chi2']:.2f} A {np.round(o['a'], 0).tolist()} E {np.round(o['E'], 0).tolist()} bslope {o['bslope']:.2f} f1-4 " + "".join(str(int(o[f'f{i}'])) for i in (1, 2, 3, 4)))
    for rt in TYPES:
        for pop in ("rand", "k2", "clip"):
            c = bss[(rt, pop)]
            print(f"[bs] {rt} {pop}: N {c['n']} gated {c['gated']} k changed {c['changed']} (gated k<=1: default {c['kdef1']}, BASE_SLOPE {c['to1']})")
    print("[clip] " + clip_recommendation(cs))
    print(f"[consistency] compared {ncons}, mismatches {len(mism)}")
    for r in mism[:20]:
        print(f"  {r['pop']} {r.get('id', '')} run {r['run']} ev {r['ev']} {r['name']}: f1 {int(fl(r, 1))} f2 {int(fl(r, 2))} vs P1 {r['p1_c1']} {r['p1_c2']}")
    print(f"[done] {len(rows)} sample rows, {len(exs)} examples, {nfail} not reproduced (fit_event), {nfail1} single-row refits differ, tsv rows {len(rows) + len(exs)}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir", default=str(p1.PYDIR / "output" / "pileup_rates" / "rawsum_sample"))
    ap.add_argument("--full-scan-dir", default=str(p1.PYDIR / "output" / "pileup_rates" / "scan"))
    ap.add_argument("--evio-dir", default=p1.EVIO_DIR)
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--cache", default=str(Path(tempfile.gettempdir()) / "exotic_p2_cache.pkl"))
    ap.add_argument("--jobs", type=int, default=3)
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    p1.EVIO_DIR = a.evio_dir
    if a.plot_only:
        data = pickle.load(open(a.cache, "rb"))
    else:
        data = compute(a)
        pickle.dump(data, open(a.cache, "wb"))
    report(data, Path(a.out) / STEM, a.sample_dir)


if __name__ == "__main__":
    main()
