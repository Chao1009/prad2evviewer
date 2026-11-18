#!/usr/bin/env python3
"""exotic_p3_sample.py -- P3 of the exotic-waveform study: rates, overlaps and examples of the four exotic flags on the 2-file sample (run_exotic_sample.sh).

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p3_sample.py [--sample-dir <dir>] [--evio-dir <dir>] [--out <dir>] [--cache <pkl>] [--jobs 8] [--plot-only]

Population: PbWO4 kept channels (scan rows with A_main >= KEEP_ADC; trigger gate and LMS veto as pileup_rates.scan) of the 13 runs x 2 splits of
$OUT_DIR/exotic_sample/run_<run>.npz (default fit, PILEUP_EXOTIC=1, options off), per run type.  Flags ch_ex_f1 (tail-baseline), f2 (wide), f3 (smooth-wide, k >= 2),
f4 (clipped-multi) as defined in pileup_rates.exotic_vars; every number below is a fraction of the kept PbWO4 channels of the run type unless stated.
Tables: flag counts and fractions; overlaps (all 16 flag patterns); k of flagged vs unflagged channels (per flag); f2 / f3 per-module rates (top 20 by rate over all run
types, modules with >= MIN_N channels; module maps per run type); f4 energy distributions E_main / E_beam and E_sec / E_beam (largest secondary) and leakage into the elastic window
|E / E_beam - 1| <= WIN3 (3 sigma of the HyCal resolution, as single_peak_beam_energy_examples): channels with E_main, with E_sec, with either, and with any secondary in the window.
Examples: per flag and run type N_EX random (seed SEED) flagged channels, redrawn from the EVIO files (event number, module), refitted with the default fit and required to reproduce the
stored k and chi2/dof (SPARE extra draws replace failures).  Each panel: data (red x: overflow), default fit (red total, dashed components, grey baseline, dotted green: the stretched
single-pulse fit of flag 3 where it ran), residual below.  f1 examples add the PILEUP_BASE_SLOPE refit, f4 examples the PILEUP_CLIP_SINGLE refit (8 ns margin, the default EX_CLIP_NS)
as a second panel (refits by exotic_p2_flags.fit_opt on the single row; same b_ref, sigma_ref).
Outputs (--out): exotic_P3_sample_20261117.pdf (table pages, then examples), .tsv (rates, overlaps, k distributions, f4 leakage), _modules.tsv, _examples.tsv.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import json
import time
import pickle
import argparse
import itertools
import tempfile
import multiprocessing as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import exotic_p1_quantities as p1                       # noqa: E402
import exotic_p2_flags as p2                            # noqa: E402
from exotic_p1_quantities import TYPES, RT_LABEL, get_ctx, fmt, write_tsv   # noqa: E402

STEM = "exotic_P3_sample_20261117"
SEED = 20261118
N_EX, SPARE = 20, 4
MIN_N = 100                                             # modules with fewer kept channels (pooled over run types) are left out of the top-20 lists
FLAGS = (1, 2, 3, 4)
FNAME = {1: "f1 tail-baseline", 2: "f2 wide", 3: "f3 smooth-wide", 4: "f4 clipped-multi"}
# energy window |E/E_beam - 1| <= WIN3 per run type: 3 sigma of the HyCal resolution (as single_peak_beam_energy_examples.py)
WIN3 = {"0p73GeV_H2": 0.095, "2p24GeV_H2": 0.052, "3p49GeV_H2": 0.044, "2p24GeV_Ta_X17": 0.060}
KBINS = ["1", "2", "3", "4", ">=5"]


def kbin(k):
    return np.minimum(k, 5) - 1


# ---------------------------------------------------------------------------------------------------------------- population
def load_pop(sample_dir):
    """Per-run-type arrays of the PbWO4 kept channels, and the per-run records needed for the examples."""
    import pileup_rates as pr
    pop = {rt: dict(run=[], ev=[], evnum=[], midx=[], k=[], chi2=[], emain=[], esec=[], eany=[], f1=[], f2=[], f3=[], f4=[], split=[]) for rt in TYPES}
    runs, names, geo_names = {}, None, None
    for fn in sorted(Path(sample_dir).glob("run_*.npz")):
        z = np.load(fn)
        run, rt = int(z["run"]), str(z["run_type"])
        if rt not in pop:
            continue
        nm = [str(n) for n in z["mod_name"]]
        names = names or nm
        assert nm == names
        cfg = pr.run_config(run)
        gain = pr.module_gain(cfg, rt, nm)
        midx = z["ch_midx"].astype(np.int64)
        sel = z["mod_pwo4"][midx]
        k, a = z["ch_k"].astype(np.int64), z["ch_a"].astype(np.float64)
        valid = np.arange(a.shape[1])[None, :] < k[:, None]
        e = np.where(valid, a, 0.0) * gain[midx][:, None] / cfg["beam"]            # E / E_beam per pulse (nan gain: nan)
        main = np.argmax(np.where(valid, a, -1.0), axis=1)
        rows = np.arange(k.size)
        emain = e[rows, main]
        e2 = e.copy()
        e2[rows, main] = 0.0
        esec = e2.max(axis=1)
        inw = np.abs(e2 - 1) <= WIN3[rt]
        inw &= valid
        inw[rows, main] = False
        d = pop[rt]
        d["run"].append(np.full(k.size, run)[sel])
        d["ev"].append(z["ch_ev"].astype(np.int64)[sel])
        d["evnum"].append(z["ev_number"][z["ch_ev"].astype(np.int64)][sel])
        d["midx"].append(midx[sel])
        d["k"].append(k[sel]), d["chi2"].append(z["ch_chi2"].astype(np.float64)[sel])
        d["emain"].append(emain[sel]), d["esec"].append(esec[sel]), d["eany"].append(inw.any(axis=1)[sel])
        for i in FLAGS:
            d[f"f{i}"].append(z[f"ch_ex_f{i}"].astype(bool)[sel])
        sp = np.array([int(s) for s in z["splits"]])
        runs[run] = dict(rt=rt, splits=sp, b=z["mod_b_ref"], s=z["mod_sig_ref"], beam=float(cfg["beam"]), wall=float(z["wall_s"]), n_pass=int(z["n_pass"]), n_seen=int(z["n_seen"]))
    for rt in TYPES:
        pop[rt] = {k: (np.concatenate(v) if v else np.zeros(0, bool if k in ('eany', 'f1', 'f2', 'f3', 'f4') else float if k in ('chi2', 'emain', 'esec') else np.int64)) for k, v in pop[rt].items()}
    return pop, runs, names


def stats(pop, names):
    """All rates and tables as dicts of plain numbers."""
    S = dict(rate={}, pat={}, kdist={}, mod={}, f4={})
    for rt in TYPES:
        d = pop[rt]
        n = d["k"].size
        F = {i: d[f"f{i}"] for i in FLAGS}
        anyf = np.any([F[i] for i in FLAGS], axis=0)
        S["rate"][rt] = dict(n=n, any=int(anyf.sum()), **{f"f{i}": int(F[i].sum()) for i in FLAGS})
        code = sum(F[i].astype(int) << (i - 1) for i in FLAGS)
        S["pat"][rt] = {c: int((code == c).sum()) for c in range(16)}
        S["kdist"][rt] = {}
        for i in FLAGS:
            for lab, sel in (("flagged", F[i]), ("unflagged", ~F[i])):
                S["kdist"][rt][(i, lab)] = np.bincount(kbin(d["k"][sel]), minlength=5)[:5]
        f4 = F[4]
        em, es = d["emain"][f4], d["esec"][f4]
        win = WIN3[rt]
        inm, ins = np.abs(em - 1) <= win, np.abs(es - 1) <= win
        S["f4"][rt] = dict(n=int(f4.sum()), em=em, es=es, main=int(inm.sum()), sec=int(ins.sum()), either=int((inm | ins).sum()), anysec=int((inm | d["eany"][f4]).sum()),
                           anysec_only=int(d["eany"][f4].sum()), win=win)
        S["mod"][rt] = {i: np.bincount(d["midx"][sel], minlength=len(names)) for i, sel in (("n", np.ones(n, bool)), ("f2", F[2]), ("f3", F[3]))}
    return S


# ---------------------------------------------------------------------------------------------------------------- examples
def draw_picks(pop, runs, rng):
    """Per run type and flag: N_EX + SPARE random flagged channels (ordered); returns {(rt, flag): [(run, evnum, midx, k, chi2)]}."""
    out = {}
    for rt in TYPES:
        d = pop[rt]
        for i in FLAGS:
            idx = np.nonzero(d[f"f{i}"])[0]
            sel = rng.permutation(idx)[:N_EX + SPARE]
            out[(rt, i)] = [(int(d["run"][j]), int(d["evnum"][j]), int(d["midx"][j]), int(d["k"][j]), float(d["chi2"][j])) for j in sel]
    return out


def work(task):
    """('run', run, files [(path, split, b_ref, sig_ref)], picks {evnum: {midx: stored k, chi2}}): channel records with the default fit, flags and option refits."""
    import event_display as ed
    from prad2py import dec
    run, files, picks = task[1], task[2], task[3]
    t_start = time.time()
    ctx = get_ctx(run)
    out, left = [], set(picks)
    for f, split, b_ref, sig_ref in files:
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
                r = ed.fit_event(ctx, d["event"], b_ref, sig_ref)
                for midx, st in picks[en].items():
                    row = r["row"].get(midx)
                    rec = dict(run=run, rt=ctx.rt, ev=en, midx=midx, name=str(ctx.names[midx]), ok=False, ok1=False, split=split)
                    if row is not None:
                        k, chi2 = int(r["k"][row]), float(r["chi2"][row])
                        rec["ok"] = k == st["k"] and bool(np.isclose(chi2, st["chi2"], rtol=2e-3, atol=1e-3))
                        rec.update(st_k=st["k"], st_chi2=st["chi2"])
                        rec.update(p2.channel(ctx, r, row, b_ref, sig_ref, True, False))
                        rec["Y"] = r["Y"][row].astype(float)
                        rec["gain"], rec["beam"] = float(ctx.gain[midx]), float(ctx.beam)
                    out.append(rec)
        ch.close()
    for en in left:
        for midx, st in picks[en].items():
            out.append(dict(run=run, rt=ctx.rt, ev=en, midx=midx, name=str(ctx.names[midx]), ok=False, ok1=False, missing=True))
    print(f"[work] run {run}: {len(out)} records, {time.time() - t_start:.0f} s", flush=True)
    return out


def compute(a, pop, runs, picks):
    t0 = time.time()
    need = {}                                           # run -> {evnum: {midx: stored}}
    for lst in picks.values():
        for run, en, midx, k, chi2 in lst:
            need.setdefault(run, {}).setdefault(en, {})[midx] = dict(k=k, chi2=chi2)
    tasks = []
    for run, pk in need.items():
        r = runs[run]
        files = [(str(Path(a.evio_dir) / f"prad_{run:06d}.evio.{s:05d}"), int(s), r["b"][j], r["s"][j]) for j, s in enumerate(r["splits"])]
        tasks.append(("run", run, files, pk))
    tasks.sort(key=lambda t: -sum(len(v) for v in t[3].values()))
    with mp.get_context("fork").Pool(a.jobs) as pool:
        res = pool.map(work, tasks, chunksize=1)
    recs = {(r["run"], r["ev"], r["midx"]): r for rr in res for r in rr}
    return recs, time.time() - t0


def select_examples(picks, recs):
    """First N_EX reproducing channels per (run type, flag) in draw order; flags of the single-row refit must agree with the scan flag (counted otherwise)."""
    ex, fails = {}, dict(not_reproduced=0, missing=0, flag_mismatch=0)
    for (rt, i), lst in picks.items():
        got = []
        for run, en, midx, k, chi2 in lst:
            r = recs.get((run, en, midx))
            if r is None or r.get("missing"):
                fails["missing"] += 1
                continue
            if not (r["ok"] and r["ok1"]):
                fails["not_reproduced"] += 1
                continue
            if len(got) < N_EX:
                if not p2.fl(r, i):
                    fails["flag_mismatch"] += 1
                got.append(r)
        ex[(rt, i)] = got
    return ex, fails


# ---------------------------------------------------------------------------------------------------------------- output
def pct(x, n):
    return f"{100 * x / max(n, 1):.3f}"


def tsv_stats(S, names):
    rows = []
    for rt in TYPES:
        r = S["rate"][rt]
        for key in ("n", "f1", "f2", "f3", "f4", "any"):
            rows.append(dict(section="rate", run_type=rt, item=key, count=r[key], percent=pct(r[key], r["n"]) if key != "n" else ""))
        for c in range(16):
            lab = "+".join(f"f{i}" for i in FLAGS if c >> (i - 1) & 1) or "none"
            rows.append(dict(section="overlap", run_type=rt, item=lab, count=S["pat"][rt][c], percent=pct(S["pat"][rt][c], r["n"])))
        for i in FLAGS:
            for lab in ("flagged", "unflagged"):
                h = S["kdist"][rt][(i, lab)]
                for b, v in zip(KBINS, h):
                    rows.append(dict(section="k_dist", run_type=rt, item=f"f{i} {lab} k={b}", count=int(v), percent=pct(v, h.sum())))
        f = S["f4"][rt]
        for key in ("n", "main", "sec", "either", "anysec"):
            rows.append(dict(section="f4_leakage", run_type=rt, item=key, count=f[key], percent=pct(f[key], f["n"]) if key != "n" else "", window=f["win"]))
    return rows


def module_rows(S, names):
    rows = []
    for rt in TYPES:
        m = S["mod"][rt]
        for j, nme in enumerate(names):
            if m["n"][j]:
                rows.append(dict(run_type=rt, module=nme, midx=j, n=int(m["n"][j]), f2=int(m["f2"][j]), f3=int(m["f3"][j]),
                                 f2_percent=pct(m["f2"][j], m["n"][j]), f3_percent=pct(m["f3"][j], m["n"][j])))
    for j, nme in enumerate(names):
        n = sum(int(S["mod"][rt]["n"][j]) for rt in TYPES)
        if n:
            f2, f3 = (sum(int(S["mod"][rt][q][j]) for rt in TYPES) for q in ("f2", "f3"))
            rows.append(dict(run_type="all", module=nme, midx=j, n=n, f2=f2, f3=f3, f2_percent=pct(f2, n), f3_percent=pct(f3, n)))
    return rows


def top_modules(S, names, q, top=20):
    n = sum(S["mod"][rt]["n"] for rt in TYPES)
    c = sum(S["mod"][rt][q] for rt in TYPES)
    rate = np.where(n >= MIN_N, c / np.maximum(n, 1), -1.0)
    o = np.argsort(-rate, kind="stable")[:top]
    return [(names[j], int(n[j]), int(c[j]), 100 * c[j] / n[j]) for j in o if rate[j] >= 0]


def text_page(pdf, plt, title, lines, size=9.5):
    fig = plt.figure(figsize=(16, 11))
    fig.text(0.04, 0.965, title, fontsize=17, fontweight="bold", va="top")
    fig.text(0.04, 0.92, "\n".join(lines), fontsize=size, va="top", family="monospace", linespacing=1.45)
    pdf.savefig(fig)
    plt.close(fig)


def rate_lines(S):
    L = [f"{'run type':<18}{'N kept PbWO4':>14}" + "".join(f"{'f' + str(i) + ' n':>9}{'%':>8}" for i in FLAGS) + f"{'any n':>9}{'%':>8}"]
    for rt in TYPES:
        r = S["rate"][rt]
        L.append(f"{RT_LABEL[rt]:<18}{r['n']:>14d}" + "".join(f"{r[f'f{i}']:>9d}{100 * r[f'f{i}'] / max(r['n'], 1):>8.3f}" for i in FLAGS) + f"{r['any']:>9d}{100 * r['any'] / max(r['n'], 1):>8.3f}")
    n = sum(S["rate"][rt]["n"] for rt in TYPES)
    L.append(f"{'all':<18}{n:>14d}" + "".join(f"{sum(S['rate'][rt][f'f{i}'] for rt in TYPES):>9d}{100 * sum(S['rate'][rt][f'f{i}'] for rt in TYPES) / n:>8.3f}" for i in FLAGS)
             + f"{sum(S['rate'][rt]['any'] for rt in TYPES):>9d}{100 * sum(S['rate'][rt]['any'] for rt in TYPES) / n:>8.3f}")
    return L


def overlap_lines(S):
    L = [f"{'pattern':<18}" + "".join(f"{RT_LABEL[rt]:>22}" for rt in TYPES)]
    for c in range(16):
        lab = "+".join(f"f{i}" for i in FLAGS if c >> (i - 1) & 1) or "none"
        L.append(f"{lab:<18}" + "".join(f"{S['pat'][rt][c]:>12d}{100 * S['pat'][rt][c] / max(S['rate'][rt]['n'], 1):>10.4f}%" for rt in TYPES))
    L.append("")
    L.append("pairwise (channels with both flags / channels with the row flag, %), all run types pooled:")
    L.append(f"{'':<8}" + "".join(f"{'f' + str(j):>9}" for j in FLAGS))
    for i in FLAGS:
        tot = sum(S["rate"][rt][f"f{i}"] for rt in TYPES)
        both = [sum(sum(S["pat"][rt][c] for c in range(16) if (c >> (i - 1) & 1) and (c >> (j - 1) & 1)) for rt in TYPES) for j in FLAGS]
        L.append(f"{'f' + str(i) + ' n=' + str(tot):<14}" + "".join(f"{100 * b / max(tot, 1):>9.1f}" for b in both))
    return L


def k_page(pdf, plt, S):
    fig, axs = plt.subplots(len(TYPES), 4, figsize=(18, 13), sharey=False)
    fig.suptitle("k distribution of flagged vs unflagged kept PbWO4 channels (fraction within each group)", fontsize=16, fontweight="bold")
    x = np.arange(5)
    for r, rt in enumerate(TYPES):
        for c, i in enumerate(FLAGS):
            ax = axs[r, c]
            for off, lab, col in ((-0.2, "flagged", "tab:red"), (0.2, "unflagged", "0.5")):
                h = S["kdist"][rt][(i, lab)]
                ax.bar(x + off, h / max(h.sum(), 1), width=0.4, color=col, label=f"{lab} (n={int(h.sum())})")
            ax.set_xticks(x)
            ax.set_xticklabels(KBINS)
            ax.set_yscale("log")
            ax.set_title(f"{RT_LABEL[rt]}: {FNAME[i]}", fontsize=10)
            ax.legend(fontsize=8)
            if r == len(TYPES) - 1:
                ax.set_xlabel("k")
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def module_pages(pdf, plt, S, names):
    import pileup_rates as pr
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    geo = {m["n"]: m["geo"] for m in json.load(open(pr.HYCAL_MAP))}
    pw = [i for i, n in enumerate(names) if n.startswith("W")]
    rect = {i: Rectangle((geo[names[i]]["x"] - geo[names[i]]["sx"] / 2, geo[names[i]]["y"] - geo[names[i]]["sy"] / 2), geo[names[i]]["sx"], geo[names[i]]["sy"]) for i in pw}
    for q in ("f2", "f3"):
        top = top_modules(S, names, q)
        lines = [f"top 20 modules by {q} rate, run types pooled, dominated by Ta X17 (modules with >= {MIN_N} kept channels; per-type rates in _modules.tsv)", f"{'module':<10}{'n kept':>10}{'n ' + q:>10}{'rate %':>10}"]
        lines += [f"{a:<10}{b:>10d}{c:>10d}{d:>10.2f}" for a, b, c, d in top]
        fig = plt.figure(figsize=(18, 13))
        fig.suptitle(f"{FNAME[int(q[1])]}: per-module rate (kept PbWO4 channels)", fontsize=16, fontweight="bold")
        fig.text(0.02, 0.94, "\n".join(lines), fontsize=9, va="top", family="monospace", linespacing=1.35)
        vals = {}
        for rt in TYPES:
            n, c = S["mod"][rt]["n"], S["mod"][rt][q]
            vals[rt] = np.where(n >= 30, 100.0 * c / np.maximum(n, 1), np.nan)
        vmax = max(1.0, np.nanpercentile(np.concatenate([v[pw] for v in vals.values()]), 99))
        for t, rt in enumerate(TYPES):
            ax = fig.add_axes([0.40 + 0.30 * (t % 2), 0.52 if t < 2 else 0.06, 0.28, 0.36])
            has = [i for i in pw if np.isfinite(vals[rt][i])]
            ax.add_collection(PatchCollection([rect[i] for i in pw if i not in has], facecolor="0.85", edgecolor="0.6", lw=0.2))
            pc = PatchCollection([rect[i] for i in has], cmap="viridis", edgecolor="0.3", lw=0.2)
            pc.set_array(vals[rt][has])
            pc.set_clim(0, vmax)
            ax.add_collection(pc)
            ax.set_xlim(-600, 600), ax.set_ylim(-600, 600)
            ax.set_aspect("equal")
            ax.set_title(f"{RT_LABEL[rt]} ({q} %, modules with >= 30 channels; grey: fewer)", fontsize=9)
            fig.colorbar(pc, ax=ax, fraction=0.04)
        pdf.savefig(fig)
        plt.close(fig)


def f4_page(pdf, plt, S):
    fig, axs = plt.subplots(2, 4, figsize=(18, 9.5))
    fig.suptitle("f4 clipped-multi: E_main / E_beam (top) and E_sec / E_beam, largest secondary (bottom); shaded: elastic window |E/E_beam - 1| <= 3 sigma", fontsize=14, fontweight="bold")
    bins = np.linspace(0, 2.0, 81)
    for c, rt in enumerate(TYPES):
        f = S["f4"][rt]
        for r, (key, lab) in enumerate((("em", "E_main"), ("es", "E_sec"))):
            ax = axs[r, c]
            v = f[key][np.isfinite(f[key])]
            ax.hist(np.clip(v, 0, 2.0), bins=bins, color="tab:orange")
            ax.axvspan(1 - f["win"], 1 + f["win"], color="tab:blue", alpha=0.2)
            ax.set_xlabel(f"{lab} / E_beam (clipped at 2)")
            ax.set_title(f"{RT_LABEL[rt]} (n={f['n']})", fontsize=10)
            ax.set_yscale("log")
    fig.tight_layout(rect=(0, 0, 1, 0.95))
    pdf.savefig(fig)
    plt.close(fig)
    L = [f"{'run type':<18}{'window':>9}{'N f4':>8}{'E_main in':>11}{'%':>7}{'E_sec in':>10}{'%':>7}{'either':>8}{'%':>7}{'main or any sec':>17}{'%':>7}"]
    for rt in TYPES:
        f = S["f4"][rt]
        L.append(f"{RT_LABEL[rt]:<18}{'+-' + str(f['win']):>9}{f['n']:>8d}{f['main']:>11d}{100 * f['main'] / max(f['n'], 1):>7.1f}{f['sec']:>10d}{100 * f['sec'] / max(f['n'], 1):>7.1f}"
                 f"{f['either']:>8d}{100 * f['either'] / max(f['n'], 1):>7.1f}{f['anysec']:>17d}{100 * f['anysec'] / max(f['n'], 1):>7.1f}")
    return L


def panel_title(e, i, tag, o=None):
    n = lambda v, s: f"{v:{s}}" if np.isfinite(v) else "nan"          # noqa: E731
    head = f"{RT_LABEL[e['rt']]} run {e['run']} ev {e['ev']} {e['name']}"
    if o is None:
        base = f"k={e['k']} chi2={e['chi2']:.2f} E_main={e['amain'] * e['gain']:.0f} MeV ({100 * e['amain'] * e['gain'] / e['beam']:.0f}% Eb)"
        fv = [f"f{j}={'Y' if p2.fl(e, j) else '.'}" for j in FLAGS]
        v = {1: f"db={n(e['ex_db'], '.1f')} tail={n(e['ex_tail'], '.1f')}", 2: f"w50={n(e['w50'], '.2f')}",
             3: f"chi2s={n(e['ex_chi2s'], '.1f')} chi2k={n(e['ex_chi2k'], '.1f')} nm={e['ex_nm']}",
             4: f"nclip={e['nclip']} dclip={n(e['ex_dclip'], '.1f')} esec={n(e['ex_esec'], '.2f')}"}[i]
        return f"{head}\n{base}\n{' '.join(fv)}  {v}"
    ex = f"bslope={o['bslope']:.2f} gate={int(o['gate'])}" if tag == "bs" else ""
    return (f"{'PILEUP_BASE_SLOPE=1' if tag == 'bs' else 'PILEUP_CLIP_SINGLE=1 (8 ns)'} refit\nk={o['k']} chi2={o['chi2']:.2f} A={o['amain']:.0f} E={o['amain'] * e['gain']:.0f} MeV {ex}"
            f"\n(default: k={e['k']} chi2={e['chi2']:.2f} A={e['amain']:.0f})")


def example_pages(pdf, plt, ex):
    ncol, nrow = 4, 3
    for rt in TYPES:
        for i in FLAGS:
            lst = ex[(rt, i)]
            two = i in (1, 4)
            tag = "bs" if i == 1 else "cs8"
            per = (ncol // 2 if two else ncol) * nrow
            if not lst:
                fig = plt.figure(figsize=(12, 3))
                fig.text(0.05, 0.5, f"{RT_LABEL[rt]} {FNAME[i]}: no flagged channels in the sample", fontsize=14)
                pdf.savefig(fig)
                plt.close(fig)
                continue
            for p0 in range(0, len(lst), per):
                chunk = lst[p0:p0 + per]
                fig = plt.figure(figsize=(24, 17))
                outer = fig.add_gridspec(nrow, ncol, hspace=0.5, wspace=0.2, left=0.03, right=0.995, top=0.88, bottom=0.03)
                fig.text(0.03, 0.975, f"{RT_LABEL[rt]}: {FNAME[i]} examples {p0 + 1}-{p0 + len(chunk)} of {len(lst)} (random; flagged channels of the sample)", fontsize=18, fontweight="bold")
                fig.text(0.03, 0.95, "black dots: raw samples (red x: overflow); red: total default fit; dashed: fitted pulses; grey: fit baseline; dotted green: stretched single pulse (f3, where run); "
                         + ("second panel of each example: " + ("PILEUP_BASE_SLOPE=1 refit." if i == 1 else "PILEUP_CLIP_SINGLE=1 refit (8 ns margin).") if two else ""), fontsize=11)
                for j, e in enumerate(chunk):
                    for part in range(2 if two else 1):
                        pos = j * (2 if two else 1) + part
                        gs = outer[pos // ncol, pos % ncol].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.12)
                        ax, axr = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
                        if part == 0:
                            p2.draw(ax, axr, e["Y"], e["base"], e["model"], e["stretch"], panel_title(e, i, tag), "0.0")
                        else:
                            o = e["opts"][tag]
                            p2.draw(ax, axr, e["Y"], o["base"], o["model"], None, panel_title(e, i, tag, o), "darkblue")
                    ax.title.set_fontsize(7.5)
                for a in fig.axes:
                    a.title.set_fontsize(7.5)
                pdf.savefig(fig)
                plt.close(fig)


def example_rows(ex):
    rows = []
    for (rt, i), lst in ex.items():
        for j, e in enumerate(lst):
            r = dict(run_type=rt, flag=f"f{i}", n_in_flag=j + 1, run=e["run"], ev=e["ev"], module=e["name"], midx=e["midx"], k=e["k"], chi2=e["chi2"], amain=e["amain"],
                     Emain_MeV=e["amain"] * e["gain"], Emain_Eb=e["amain"] * e["gain"] / e["beam"], nclip=e["nclip"], w50=e["w50"], ex_db=e["ex_db"], ex_tail=e["ex_tail"],
                     ex_s=e["ex_s"], ex_chi2s=e["ex_chi2s"], ex_chi2k=e["ex_chi2k"], ex_nm=e["ex_nm"], ex_dclip=e["ex_dclip"], ex_esec=e["ex_esec"],
                     **{f"f{q}": int(p2.fl(e, q)) for q in FLAGS})
            for tag in ("bs", "cs8"):
                o = e["opts"][tag]
                r.update({f"{tag}_k": o["k"], f"{tag}_chi2": o["chi2"], f"{tag}_amain": o["amain"], f"{tag}_Emain_MeV": o["amain"] * e["gain"]})
            r["bs_bslope"], r["bs_gate"] = e["opts"]["bs"]["bslope"], int(e["opts"]["bs"]["gate"])
            rows.append(r)
    return rows


def report(data, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    S, names, ex, fails = data["S"], data["names"], data["ex"], data["fails"]
    stem = Path(out) / STEM
    notes = [f"example counts per run type and flag (target {N_EX}): " + ", ".join(f"{RT_LABEL[rt]} " + "/".join(str(len(ex[(rt, i)])) for i in FLAGS) for rt in TYPES) + "  (f1/f2/f3/f4)",
             f"example draws replaced or flagged: {fails}; scan wall time of the 13 runs (sum): {data['scan_wall']:.0f} s"]
    with PdfPages(str(stem) + ".pdf") as pdf:
        text_page(pdf, plt, "P3: exotic flags on the 2-file sample (default fit, PILEUP_EXOTIC=1, options off)",
                  ["population: kept PbWO4 channels (A_main >= KEEP_ADC) of the 13 runs x 2 splits of exotic_sample, per run type; fractions of that population",
                   "f1 tail-baseline (ex_db > 8, ex_tail < -5); f2 wide (w50 > 3); f3 smooth-wide (k >= 2, stretched single pulse within the BIC allowance); f4 clipped-multi (nclip >= 1, secondary within 8 ns of the clip or E_sec > 0.5 E_beam)", ""]
                  + rate_lines(S) + [""] + notes, size=10)
        text_page(pdf, plt, "Overlaps of the four flags (all 16 patterns; count and % of kept PbWO4 channels)", overlap_lines(S))
        k_page(pdf, plt, S)
        module_pages(pdf, plt, S, names)
        f4_lines = f4_page(pdf, plt, S)
        text_page(pdf, plt, "f4 clipped-multi: leakage of the flagged channels into the elastic window (+-3 sigma)",
                  ["E_main in / E_sec in: flagged channels whose main pulse / largest secondary has |E/E_beam - 1| <= window; either: E_main or E_sec; main or any sec: E_main or any secondary pulse", ""] + f4_lines)
        example_pages(pdf, plt, ex)
    write_tsv(str(stem) + ".tsv", tsv_stats(S, names), ["section", "run_type", "item", "count", "percent", "window"])
    write_tsv(str(stem) + "_modules.tsv", module_rows(S, names), ["run_type", "module", "midx", "n", "f2", "f3", "f2_percent", "f3_percent"])
    cols = ["run_type", "flag", "n_in_flag", "run", "ev", "module", "midx", "k", "chi2", "amain", "Emain_MeV", "Emain_Eb", "nclip", "w50", "ex_db", "ex_tail", "ex_s", "ex_chi2s", "ex_chi2k",
            "ex_nm", "ex_dclip", "ex_esec", "f1", "f2", "f3", "f4"] + [f"{t}_{q}" for t in ("bs", "cs8") for q in ("k", "chi2", "amain", "Emain_MeV")] + ["bs_bslope", "bs_gate"]
    write_tsv(str(stem) + "_examples.tsv", example_rows(ex), cols)
    print("\n".join(rate_lines(S) + [""] + overlap_lines(S) + [""] + f4_lines + [""] + notes))
    for q in ("f2", "f3"):
        print(f"[top {q}] " + "; ".join(f"{a} {d:.1f}% ({c}/{b})" for a, b, c, d in top_modules(S, names, q)))
    for rt in TYPES:
        for i in FLAGS:
            print(f"[examples] {rt} f{i}: {len(ex[(rt, i)])}")
    for (rt, i) in (k for k in ex if k[1] in (1, 4)):
        tag = "bs" if i == 1 else "cs8"
        L = ex[(rt, i)]
        if L:
            dk = np.array([e["opts"][tag]["k"] - e["k"] for e in L])
            ch = np.array([e["opts"][tag]["chi2"] / e["chi2"] for e in L])
            print(f"[opt] {rt} f{i} {tag}: k changed {int((dk != 0).sum())}/{len(L)}, k lower {int((dk < 0).sum())}, ended k=1 {sum(e['opts'][tag]['k'] == 1 for e in L)}, "
                  f"ended k=0 (no pulse, chi2 nan) {sum(e['opts'][tag]['k'] == 0 for e in L)}, median chi2 ratio (finite) {np.nanmedian(ch) if np.isfinite(ch).any() else float('nan'):.2f}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir", default=str(p1.PYDIR / "output" / "pileup_rates" / "exotic_sample"))
    ap.add_argument("--evio-dir", default=p1.EVIO_DIR)
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--cache", default=str(Path(tempfile.gettempdir()) / "exotic_p3_cache.pkl"))
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    p1.EVIO_DIR = a.evio_dir
    if a.plot_only:
        data = pickle.load(open(a.cache, "rb"))
    else:
        t0 = time.time()
        pop, runs, names = load_pop(a.sample_dir)
        S = stats(pop, names)
        picks = draw_picks(pop, runs, np.random.default_rng(SEED))
        recs, el = compute(a, pop, runs, picks)
        ex, fails = select_examples(picks, recs)
        data = dict(S=S, names=names, ex=ex, fails=fails, scan_wall=sum(r["wall"] for r in runs.values()), elapsed=time.time() - t0)
        pickle.dump(data, open(a.cache, "wb"))
        print(f"[compute] {data['elapsed']:.0f} s")
    report(data, a.out)


if __name__ == "__main__":
    main()
