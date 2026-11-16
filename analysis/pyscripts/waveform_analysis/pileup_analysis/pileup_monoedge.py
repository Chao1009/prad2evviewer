#!/usr/bin/env python3
"""pileup_monoedge.py -- bounded old-vs-new comparison of the monotonic-edge rule and the noisy flag (pileup_rates.py MONO_EDGE, NOISE_FLAG).

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/pileup_monoedge.py splits [--evio-dir <dir>]      # 2 random EVIO splits per run (seed fixed): "<run> <file> <file>"
    bash waveform_analysis/pileup_analysis/run_monoedge_sample.sh                               # the old and new scans of those splits
    python3 waveform_analysis/pileup_analysis/pileup_monoedge.py report [--sample-dir <dir>] [--out <dir>] [--stem <stem>]

report: tables (old vs new R_N, noisy fraction, top noisy modules), 20 examples per N = 2..6 refitted with the NEW fit, 20 noisy-flagged examples
-> <stem>.pdf, <stem>_table.tsv, <stem>_examples.tsv, <stem>_noisy_modules.tsv.  Definitions are printed on the first PDF page.
"""
from __future__ import annotations

import os
import sys
import glob
import argparse
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
PYDIR = HERE.parent.parent
SEED = 20261115
SPLITS_PER_RUN = 2
EVIO_DIR = "/Users/L00338853/work/PRad/data/evio"


def splits_main(a):
    """Print '<run> <file> <file>' per run of run_pileup_scan.sh: SPLITS_PER_RUN random splits (seed SEED + run), sorted."""
    import re
    import numpy as np
    text = open(HERE / "run_pileup_scan.sh").read()
    runs = []
    for m in re.finditer(r"^RUNS_\w+=\(([\d ]+)\)", text, re.M):
        runs += [int(x) for x in m.group(1).split()]
    for run in runs:
        files = sorted(glob.glob(str(Path(a.evio_dir) / f"prad_{run:06d}.evio.*")))
        pick = sorted(np.random.default_rng(SEED + run).choice(len(files), size=min(SPLITS_PER_RUN, len(files)), replace=False).tolist())
        print(run, *[files[i] for i in pick])


if __name__ == "__main__" and len(sys.argv) > 1 and sys.argv[1] == "splits":
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd")
    ap.add_argument("--evio-dir", default=EVIO_DIR)
    splits_main(ap.parse_args())
    sys.exit(0)

# ---------------------------------------------------------------------------------------------------------------- report
os.environ["PILEUP_MONO_EDGE"] = "1"                    # examples are refitted with the NEW fit (and the noisy flag); set before pileup_rates is imported
os.environ["PILEUP_NOISE_FLAG"] = "1"
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "exotic_events"))
import numpy as np                                      # noqa: E402
import matplotlib                                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402
import pileup_multiplicity as pm                        # noqa: E402  (panel format, rate formatting, example selection; imports event_display and pileup_rates)

pr, ed = pm.pr, pm.ed
F32 = np.float32
NS_LIST, N_EX, N_PER_RT = pm.NS_LIST, pm.N_EX, pm.N_PER_RT
POPS = pm.POPS
VARIANTS = ("old", "new")
SHADE_RISE, SHADE_FALL = "#f4a259", "#5b8def"


def load_sample(sample_dir):
    """{run: {"old": npz dict, "new": npz dict}} (arrays loaded) of the bounded sample."""
    out = {}
    for fn in sorted((Path(sample_dir) / "new").glob("run_*.npz")):
        run = int(fn.stem.split("_")[1])
        if not (Path(sample_dir) / "old" / fn.name).exists():
            sys.exit(f"[ERROR] run {run}: no old scan {Path(sample_dir) / 'old' / fn.name}")
        out[run] = {v: dict(np.load(Path(sample_dir) / v / fn.name)) for v in VARIANTS}
        if not np.array_equal(out[run]["old"]["splits"], out[run]["new"]["splits"]):
            sys.exit(f"[ERROR] run {run}: old and new scans cover different splits")
    return out


def tables(sample):
    """res[pop][rt][variant] = dict(n = non-noisy counts by k [KMAX + 1], kept, noisy, mod_kept/mod_noisy per module, runs); sync = old/new matching."""
    nm = len(next(iter(sample.values()))["new"]["mod_name"])
    res = {p: {rt: {v: dict(n=np.zeros(pr.KMAX + 1, np.int64), kept=0, noisy=0, mk=np.zeros(nm, np.int64), mn=np.zeros(nm, np.int64)) for v in VARIANTS}
               for rt in pr.RUN_TYPES} for p, _ in POPS}
    sync = dict(both=0, old_only=0, new_only=0, noisy_diff=0, runs={rt: [] for rt in pr.RUN_TYPES}, n_seen={rt: 0 for rt in pr.RUN_TYPES}, n_pass={rt: 0 for rt in pr.RUN_TYPES})
    for run, d in sorted(sample.items()):
        rt = str(d["new"]["run_type"])
        sync["runs"][rt].append(run)
        sync["n_seen"][rt] += int(d["new"]["n_seen"])
        sync["n_pass"][rt] += int(d["new"]["n_pass"])
        for v in VARIANTS:
            z = d[v]
            k, mi, nz = z["ch_k"].astype(np.int64), z["ch_midx"].astype(np.int64), z["ch_noisy"].astype(bool)
            pw = z["mod_pwo4"][mi]
            for p, sel in POPS:
                m0 = np.ones(k.size, bool) if sel is None else (pw == sel)
                r = res[p][rt][v]
                m = m0 & ~nz
                r["n"] += np.bincount(k[m], minlength=pr.KMAX + 1)[:pr.KMAX + 1]
                r["kept"] += int(m0.sum())
                r["noisy"] += int((m0 & nz).sum())
                if p == POPS[0][0]:
                    r["mk"] += np.bincount(mi, minlength=nm)
                    r["mn"] += np.bincount(mi[nz], minlength=nm)
        key = lambda z: z["ch_ev"].astype(np.int64) * 4096 + z["ch_midx"].astype(np.int64)    # noqa: E731
        ko, kn = key(d["old"]), key(d["new"])
        both, io, inn = np.intersect1d(ko, kn, return_indices=True)
        sync["both"] += both.size
        sync["old_only"] += ko.size - both.size
        sync["new_only"] += kn.size - both.size
        sync["noisy_diff"] += int((d["old"]["ch_noisy"][io] != d["new"]["ch_noisy"][inn]).sum())
    return res, sync


def rate_rows(res):
    """Rows for the TSV and the table pages: per population, run type, N in 1..6 and '>=4'."""
    rows = []
    for p, _ in POPS:
        for rt in pr.RUN_TYPES:
            r = {}
            for v in VARIANTS:
                n, tot = res[p][rt][v]["n"], int(res[p][rt][v]["n"][1:].sum())
                for N in list(range(1, pr.KMAX + 1)) + [">=4"]:
                    c = int(n[4:].sum()) if N == ">=4" else int(n[N])
                    f = c / tot if tot else np.nan
                    r[(v, N)] = (c, tot, f, np.sqrt(f * (1 - f) / tot) if tot else np.nan)
            for N in list(range(1, pr.KMAX + 1)) + [">=4"]:
                rows.append(dict(pop=p, rt=rt, N=N, old=r[("old", N)], new=r[("new", N)],
                                 noisy_old=res[p][rt]["old"]["noisy"], kept_old=res[p][rt]["old"]["kept"],
                                 noisy_new=res[p][rt]["new"]["noisy"], kept_new=res[p][rt]["new"]["kept"]))
    return rows


def write_table(path, rows):
    with open(path, "w") as f:
        f.write("# R_N = waveforms with exactly N fitted pulses / waveforms with >= 1 pulse, noisy-flagged waveforms excluded from numerator and denominator; old = rule off, new = rule on, same waveforms\n")
        f.write("population\trun_type\tN\told_count\told_denominator\told_rate\told_err_binom\tnew_count\tnew_denominator\tnew_rate\tnew_err_binom\t"
                "noisy_old\tkept_old\tnoisy_frac_old\tnoisy_new\tkept_new\tnoisy_frac_new\n")
        for r in rows:
            o, n = r["old"], r["new"]
            fo = r["noisy_old"] / r["kept_old"] if r["kept_old"] else float("nan")
            fn = r["noisy_new"] / r["kept_new"] if r["kept_new"] else float("nan")
            f.write(f"{r['pop']}\t{r['rt']}\t{r['N']}\t{o[0]}\t{o[1]}\t{o[2]:.6g}\t{o[3]:.3g}\t{n[0]}\t{n[1]}\t{n[2]:.6g}\t{n[3]:.3g}\t"
                    f"{r['noisy_old']}\t{r['kept_old']}\t{fo:.6g}\t{r['noisy_new']}\t{r['kept_new']}\t{fn:.6g}\n")


def top_noisy(res, names, rts, n=10):
    """Top n modules by noisy count (NEW scan; flags are identical in OLD) over the run types rts: [(name, noisy, kept, fraction)]."""
    mn = sum(res[POPS[0][0]][rt]["new"]["mn"] for rt in rts)
    mk = sum(res[POPS[0][0]][rt]["new"]["mk"] for rt in rts)
    o = np.argsort(-mn, kind="stable")[:n]
    return [(str(names[i]), int(mn[i]), int(mk[i]), mn[i] / mk[i] if mk[i] else float("nan")) for i in o]


def write_noisy_modules(path, res, names):
    with open(path, "w") as f:
        f.write("scope\trank\tmodule\tnoisy_count\tkept_waveforms\tnoisy_fraction\n")
        for scope, rts in [("all run types", pr.RUN_TYPES)] + [(rt, [rt]) for rt in pr.RUN_TYPES]:
            for i, (nm, c, k, fr) in enumerate(top_noisy(res, names, rts), 1):
                f.write(f"{scope}\t{i}\t{nm}\t{c}\t{k}\t{fr:.6g}\n")


# ------------------------------------------------------------------------------------------------------------------ pages
def title_page(pdf, res, sync, ex_counts, noisy_n, files):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.03, 0.95, "Pile-up multiplicity with the monotonic-edge rule and the noisy flag (bounded old-vs-new sample)", fontsize=20, fontweight="bold")
    txt = (
        "Monotonic-edge rule (pileup_rates.py MONO_EDGE; OLD = rule off, NEW = rule on; nothing else differs)\n"
        f"  Dominant peak P = largest sample of the waveform (clipped waveform: first overflow sample ends the rise window, last one starts the fall window).\n"
        f"  sigma = module noise reference sigma_ref; d_j = y[j+1] - y[j]; a difference touching an overflow sample (>= {pr.OVERFLOW}) counts as monotonic.\n"
        f"  Rising-monotonic : every d ending at the {pr.MONO_RISE_N} samples up to P is > -2 sqrt(2) sigma.      Falling-monotonic: every d ending at the {pr.MONO_FALL_N} samples after P (40 ns, cut at the\n"
        "  waveform end) is < +2 sqrt(2) sigma.  Main pulse = stage A: highest candidate; stage C: largest fitted pulse; its onset t_main.\n"
        "  Rising-monotonic: no additional pulse with onset in [t_main, P].   Falling-monotonic: none with onset in [P, P + 40 ns].   Both: [t_main, P + 40 ns].\n"
        "  Pulses outside the stretch (earlier, or later than P + 40 ns) stay allowed.  Candidates inside the stretch are rejected in stage A and stage C, before the fit.\n"
        "  The joint fit may still move an admitted pulse's fitted onset (a few ns) into the stretch; the rule acts on candidates, not on final onsets.\n\n"
        "Noisy flag (NOISE_FLAG; computed from the single-pulse fit of the stage A main candidate, so it is the same in OLD and NEW; the multi-pulse fit still runs)\n"
        f"  Fit one pulse (main only) + baseline; r_j = y_j - model_j; V_j = sigma_ref^2 + 1/12 + [A shape_sigma(t_j - t0) sqrt(a_band / max(A, {pr.A_MIN_SCALE:g}))]^2 (as the chi2/dof).\n"
        f"  Sample marked when |r_j| > {pr.NOISE_NSIG:g} sqrt(V_j) (never an overflow sample).  Noisy = some {pr.NOISE_WIN}-sample (40 ns) window holds >= {pr.NOISE_MARKS} marked samples with >= {pr.NOISE_FLIPS} sign\n"
        "  flips between consecutive marked samples (time order).\n\n"
        "Rates: R_N = waveforms with exactly N fitted pulses / waveforms with >= 1 fitted pulse (kept channels: A_main >= KEEP_ADC), both on NON-NOISY waveforms; binomial error.\n"
        "  Noisy fraction = noisy / all kept waveforms.  Selection of pileup_rates.py defaults (trigger gate, LMS veto, bad-waveform removal); saturated and chi2/dof > 1000 waveforms included.\n"
        f"Sample: N is capped at KMAX = {pr.KMAX}.  2 random EVIO splits per run ({sum(len(v) for v in sync['runs'].values())} runs, seed {SEED}), same splits for OLD and NEW.  Events seen / passing the trigger:\n  "
        + "; ".join(f"{rt} {sync['n_seen'][rt]:,}/{sync['n_pass'][rt]:,}" for rt in pr.RUN_TYPES) + "\n"
        f"  Channels in both scans {sync['both']:,}; only in OLD {sync['old_only']:,}; only in NEW {sync['new_only']:,} (kept-channel set differs only by the {pr.KEEP_ADC:g} ADC cut); noisy flag differs on {sync['noisy_diff']} common channels.\n"
        "Examples: NEW fit, non-noisy, not saturated, chi2/dof <= 1000, random (seed fixed), 5 per run type first.  Panels: black dots raw samples (red x overflow);\n"
        "  grey pedestal; red total;\n"
        "  dashed pulses; orange band = protected stretch of a rising-monotonic edge, blue band = falling-monotonic one (from the largest pulse of the NEW fit).\n"
        "  Examples found: " + ", ".join(f"N={N}: {len(ex_counts[N])}" for N in NS_LIST) + f"; noisy: {noisy_n} (requested {N_EX} each).\n")
    fig.text(0.03, 0.915, txt, fontsize=9.4, va="top", family="monospace", linespacing=1.32)
    pdf.savefig(fig)
    plt.close(fig)


def table_page(pdf, res, rows, pop):
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0.02, 0.08, 0.96, 0.76])
    ax.axis("off")
    cols = list(range(1, pr.KMAX + 1)) + [">=4"]
    hdr = ["run type", "n_all old / new", ""] + [f"N={N}" if N != ">=4" else "N>=4" for N in cols]
    cells, kinds = [], []
    for kind in ("new R_N (%)", "old R_N (%)", "new count", "old count"):
        for rt in pr.RUN_TYPES:
            sub = {r["N"]: r for r in rows if r["pop"] == pop and r["rt"] == rt}
            v = "new" if kind.startswith("new") else "old"
            body = [fmt(sub[N][v]) if "R_N" in kind else f"{sub[N][v][0]:,}" for N in cols]
            cells.append([rt, f"{sub[1]['old'][1]:,} / {sub[1]['new'][1]:,}", kind] + body)
            kinds.append(kind)
    t = ax.table(cellText=cells, colLabels=hdr, loc="upper center", cellLoc="center", colWidths=[0.1, 0.13, 0.1] + [0.0875] * 7)
    t.auto_set_font_size(False)
    t.set_fontsize(9)
    t.scale(1, 1.5)
    col = {"new R_N (%)": "#eef3fa", "old R_N (%)": "#e8f3e8", "new count": "#ffffff", "old count": "#f6f2e8"}
    for (i, j), c in t.get_celld().items():
        if i == 0:
            c.set_text_props(fontweight="bold")
        else:
            c.set_facecolor(col[kinds[i - 1]])
    fig.suptitle(f"R_N old vs new (non-noisy waveforms): {pop}", fontsize=18, fontweight="bold")
    fig.text(0.03, 0.03, "R_N = n_N / n_all, noisy-flagged waveforms excluded, binomial error (statistics only, 2 splits per run).  Blue: NEW (rule on); green: OLD (rule off); same waveforms.", fontsize=10)
    pdf.savefig(fig)
    plt.close(fig)


def fmt(x):
    return pm.fmt_rate(x[2], x[3])


def noisy_page(pdf, res, names):
    fig = plt.figure(figsize=(16, 9))
    fig.suptitle("Noisy fraction and the modules with the most noisy waveforms", fontsize=18, fontweight="bold")
    lines = ["Noisy fraction = noisy waveforms / kept waveforms (NEW scan; the flag is identical in OLD on the common channels)", ""]
    lines.append(f"{'run type':16s}" + "".join(f"{p:>34s}" for p, _ in POPS))
    for rt in pr.RUN_TYPES:
        cells = []
        for p, _ in POPS:
            r = res[p][rt]["new"]
            cells.append(f"{r['noisy']:,} / {r['kept']:,} = {100 * r['noisy'] / max(r['kept'], 1):.4f} %")
        lines.append(f"{rt:16s}" + "".join(f"{c:>34s}" for c in cells))
    fig.text(0.03, 0.92, "\n".join(lines), fontsize=10.5, va="top", family="monospace")
    scopes = [("all run types", pr.RUN_TYPES)] + [(rt, [rt]) for rt in pr.RUN_TYPES]
    for j, (scope, rts) in enumerate(scopes):
        txt = f"{scope}\n{'rank module':>12s}{'noisy':>8s}{'kept':>9s}{'frac %':>9s}\n" + "\n".join(
            f"{i:>3d} {nm:>8s}{c:8d}{k:9d}{100 * fr:9.3f}" for i, (nm, c, k, fr) in enumerate(top_noisy(res, names, rts), 1))
        fig.text(0.02 + 0.2 * j, 0.62, txt, fontsize=7.8, va="top", family="monospace")
    fig.text(0.03, 0.03, "Top 10 modules by noisy count, with their noisy fraction among the module's kept waveforms (bounded sample).", fontsize=10)
    pdf.savefig(fig)
    plt.close(fig)


# --------------------------------------------------------------------------------------------------------------- examples
def event_files(sample, evio_dir):
    """{run: [file per split of the scan, in scan order]}."""
    return {run: [str(Path(evio_dir) / f"prad_{run:06d}.evio.{int(s):05d}") for s in d["new"]["splits"]] for run, d in sample.items()}


def candidates(sample, files, ctxs, rng):
    """(N-candidates {N: [...]}, noisy candidates [...]) of the NEW scans: dict(run, rt, ev, midx, file, k)."""
    cand = {N: [] for N in NS_LIST}
    noisy = []
    for run, d in sorted(sample.items()):
        z = d["new"]
        rt = str(z["run_type"])
        ctx = ctxs.setdefault(run, ed.Ctx(run))
        fl = files[run]
        firsts = np.array([ed.first_event_number(ctx, f) for f in fl], dtype=np.int64)
        evn = z["ev_number"][z["ch_ev"]]
        fidx = np.clip(np.searchsorted(firsts, evn, side="right") - 1, 0, len(fl) - 1)
        k, nz, ncl = z["ch_k"], z["ch_noisy"].astype(bool), z["ch_nclip"]
        ok = (ncl < pr.SAT_NCLIP) & ~(z["ch_chi2"] > pr.CHI2_MAX) & ~nz
        for i in range(len(fl)):
            for N in NS_LIST:
                idx = np.nonzero(ok & (fidx == i) & (k == N))[0]
                for j in idx[rng.permutation(idx.size)[:200]]:
                    cand[N].append(dict(run=run, rt=rt, ev=int(evn[j]), midx=int(z["ch_midx"][j]), file=fl[i], k=N))
            idx = np.nonzero((ncl < pr.SAT_NCLIP) & nz & (fidx == i))[0]
            for j in idx[rng.permutation(idx.size)[:200]]:
                noisy.append(dict(run=run, rt=rt, ev=int(evn[j]), midx=int(z["ch_midx"][j]), file=fl[i], k=int(k[j])))
        print(f"[cand] run {run} {rt}", flush=True)
    return cand, noisy


def spread(pool, rng):
    """(first, rest): N_PER_RT per run type first, the others after (as pileup_multiplicity.choose)."""
    order = rng.permutation(len(pool))
    first = []
    for rt in pr.RUN_TYPES:
        first += [pool[i] for i in order if pool[i]["rt"] == rt][:N_PER_RT]
    used = {id(x) for x in first}
    return first, [pool[i] for i in order if id(pool[i]) not in used]


def extras(ctx, e, b_ref, sig_ref):
    """Protected stretch of the NEW fit's main pulse and the single-pulse noisy check of one example (the front end of pileup_rates.run_ref)."""
    mi = e["midx"]
    Yf = e["Y"].astype(F32)[None]
    ped = pr.pedestal(Yf)
    sig = np.array([sig_ref[mi]], F32)
    ped = np.where(np.abs(ped - b_ref[mi]) > pr.START_NSIGMA * sig, F32(b_ref[mi]), ped).astype(F32)
    y = Yf - ped[:, None]
    ov = Yf >= pr.OVERFLOW
    w = (~ov).astype(F32) if ov.any() else None
    ys = y * w if w is not None else y
    m = ctx.bank.row[[mi]]
    tm = float(e["t"][int(np.argmax(e["a"]))])
    prot = ctx.mf.protected(ys, w, sig, np.array([tm]), np.array([True]))
    t0, valid = ctx.mf.find_candidates(ys, w, m, sig)
    noisy, marked, model = ctx.mf.noise_check(ys, w, m, sig, t0[:, 0], valid[:, 0])
    return dict(prot=[float(x[0]) for x in prot], noisy=bool(noisy[0]), marked=marked[0], model=(model[0] + ped[0]).astype(float))


def draw_panel(ax, e, ex, noisy_view=False):
    pm.draw_panel(ax, e)
    lo1, hi1, lo2, hi2 = ex["prot"]
    note = []
    if hi1 >= lo1:
        ax.axvspan(lo1, hi1, color=SHADE_RISE, alpha=0.20, lw=0)
        note.append("rise")
    if hi2 >= lo2:
        ax.axvspan(lo2, hi2, color=SHADE_FALL, alpha=0.15, lw=0)
        note.append("fall")
    if noisy_view:
        ax.plot(ed.TT, ex["model"], ":", color="tab:green", lw=1.6)
        mk = ex["marked"]
        ax.plot(ed.TT[mk], e["Y"][mk], "o", mfc="none", mec="magenta", ms=9, mew=1.6)
        extra = f"NOISY: {int(mk.sum())} marked samples (single-pulse fit dotted green)"
    else:
        extra = "protected: " + ("+".join(note) if note else "none")
    t = ax.get_title(loc="left")
    ax.set_title(extra + "\n" + t, fontsize=7.5, loc="left", pad=3)


def example_pages(pdf, ex, extra_of, noisy_list, noisy_extra):
    for N in NS_LIST:
        lst = ex[N]
        fig, axs = plt.subplots(4, 5, figsize=(22, 17))
        for ax in axs.flat:
            ax.axis("off")
        for ax, e in zip(axs.flat, lst):
            ax.axis("on")
            draw_panel(ax, e, extra_of[id(e)])
        fig.suptitle(f"NEW fit, N = {N} fitted pulses: {len(lst)} random examples (requested {N_EX})" + ("  -- FEWER THAN REQUESTED" if len(lst) < N_EX else ""), fontsize=18, fontweight="bold")
        fig.text(0.01, 0.003, "black dots: raw samples (red x: overflow >= 4095); grey: fitted pedestal; red: total model; dashed: pulses in time order; orange band: protected stretch of a rising-monotonic "
                 "edge, blue band: falling-monotonic one. x: t (ns), y: ADC.", fontsize=11)
        fig.tight_layout(rect=(0, 0.015, 1, 0.96), h_pad=1.2)
        pdf.savefig(fig)
        plt.close(fig)
    fig, axs = plt.subplots(4, 5, figsize=(22, 17))
    for ax in axs.flat:
        ax.axis("off")
    for ax, e in zip(axs.flat, noisy_list):
        ax.axis("on")
        draw_panel(ax, e, noisy_extra[id(e)], noisy_view=True)
    fig.suptitle(f"Noisy-flagged waveforms (NEW multi-pulse fit shown): {len(noisy_list)} random examples (requested {N_EX})" + ("  -- FEWER THAN REQUESTED" if len(noisy_list) < N_EX else ""),
                 fontsize=18, fontweight="bold")
    fig.text(0.01, 0.003, "dotted green: single-pulse fit (main pulse + baseline) used by the flag; magenta rings: samples with |r| > 5 sqrt(V); other symbols as on the N pages.  "
             "The title shows the multi-pulse fit's N pulses (A ADC @ t0 ns).", fontsize=11)
    fig.tight_layout(rect=(0, 0.015, 1, 0.96), h_pad=1.2)
    pdf.savefig(fig)
    plt.close(fig)


def report_main(a):
    t_start = time.time()
    out = Path(a.out)
    sample = load_sample(a.sample_dir)
    res, sync = tables(sample)
    rows = rate_rows(res)
    names = next(iter(sample.values()))["new"]["mod_name"]
    write_table(out / f"{a.stem}_table.tsv", rows)
    write_noisy_modules(out / f"{a.stem}_noisy_modules.tsv", res, names)
    print(f"[table] {time.time() - t_start:.0f} s", flush=True)
    ex = {N: [] for N in NS_LIST}
    nlist, extra_of, noisy_extra = [], {}, {}
    if not a.no_examples:
        files = event_files(sample, a.evio_dir)
        rng = np.random.default_rng(SEED)
        ctxs = {}
        cand, ncand = candidates(sample, files, ctxs, rng)
        print({N: len(v) for N, v in cand.items()}, "noisy", len(ncand), flush=True)
        pools = {N: spread(cand[N], rng) for N in NS_LIST}
        npool = spread(ncand, rng)
        refs = {}
        for run, d in sample.items():
            for f, br, sr in zip(files[run], d["new"]["mod_b_ref"], d["new"]["mod_sig_ref"]):
                refs[f] = (br, sr)
        need, wanted = {}, {}
        for first, rest in list(pools.values()) + [npool]:
            for x in first + rest[:40]:
                need.setdefault(x["file"], set()).add(x["ev"])
                wanted.setdefault(x["file"], {}).setdefault(x["ev"], []).append(x)
        store = pm.decode_wanted(ctxs, refs, wanted, need)
        def build(x, N):
            r = store.get((x["file"], x["ev"]))
            e = pm.make_example(ctxs[x["run"]], r, x, N) if r is not None else None
            if e is None:
                return None
            return e, extras(ctxs[x["run"]], e, *refs[x["file"]])
        for N, (first, rest) in pools.items():
            for x in first + rest:
                if len(ex[N]) >= N_EX:
                    break
                b = build(x, N)
                if b is not None and not b[1]["noisy"]:
                    ex[N].append(b[0])
                    extra_of[id(b[0])] = b[1]
            print(f"[ex] N={N}: {len(ex[N])} of {N_EX}", flush=True)
        first, rest = npool
        for x in first + rest:
            if len(nlist) >= N_EX:
                break
            b = build(x, x["k"])
            if b is not None and b[1]["noisy"]:
                nlist.append(b[0])
                noisy_extra[id(b[0])] = b[1]
        print(f"[ex] noisy: {len(nlist)} of {N_EX}", flush=True)
        with open(out / f"{a.stem}_examples.tsv", "w") as f:
            f.write("kind\tN\trun\tevent\tmodule\trun_type\tchi2_dof\tpulses_A_t0\tprotected_stretch_ns\tnoisy_marked_samples\n")
            for kind, lst, xt in [(f"N={N}", ex[N], extra_of) for N in NS_LIST] + [("noisy", nlist, noisy_extra)]:
                for e in lst:
                    p = xt[id(e)]["prot"]
                    st = ",".join(f"[{p[i]:.0f},{p[i + 1]:.0f}]" for i in (0, 2) if p[i + 1] >= p[i]) or "none"
                    f.write(f"{kind}\t{e['N']}\t{e['run']}\t{e['ev']}\t{e['name']}\t{e['rt']}\t{e['chi2']:.3f}\t" + ";".join(f"{x:.1f}@{y:.1f}" for x, y in zip(e["a"], e["t"]))
                            + f"\t{st}\t{int(xt[id(e)]['marked'].sum())}\n")
    with PdfPages(out / f"{a.stem}.pdf") as pdf:
        title_page(pdf, res, sync, ex, len(nlist), None)
        for p, _ in POPS:
            table_page(pdf, res, rows, p)
        noisy_page(pdf, res, names)
        if not a.no_examples:
            example_pages(pdf, ex, extra_of, nlist, noisy_extra)
    print(f"[sync] both {sync['both']} old_only {sync['old_only']} new_only {sync['new_only']} noisy_diff {sync['noisy_diff']}")
    print(f"[done] {time.time() - t_start:.0f} s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["splits", "report"])
    ap.add_argument("--sample-dir", default=str(PYDIR / "output" / "pileup_rates" / "monoedge_sample"))
    ap.add_argument("--evio-dir", default=EVIO_DIR)
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--stem", default="pileup_multiplicity_20261115_monoedge")
    ap.add_argument("--no-examples", action="store_true")
    a = ap.parse_args()
    report_main(a)


if __name__ == "__main__":
    main()
