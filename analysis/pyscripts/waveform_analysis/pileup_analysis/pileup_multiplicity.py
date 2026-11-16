#!/usr/bin/env python3
"""pileup_multiplicity.py -- fraction of waveforms fitted with N = 1..6 pulses per run type, plus 20 refitted examples per N.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/pileup_multiplicity.py [--scan-dir <dir>] [--evio-dir <dir>] [--out <dir>]

Table: from the pileup_rates.py scan files (<scan-dir>/run_*.npz, ch_k = number of fitted pulses of each kept channel).
R_N = n_N / n_all, n_all = kept channels (>= 1 fitted pulse, A_main >= KEEP_ADC); pulse-weighted P_N = N n_N / sum_k k n_k.
Examples: random (seed SEED) kept channels with k == N, not saturated (nclip < SAT_NCLIP) and chi2/dof <= CHI2_MAX, from FILES_PER_RUN
random EVIO splits per run (bounds the noise-reference reads); the event is decoded again and refitted with pileup_rates.fit
(event_display.noise_reference / find-event logic, same per-file noise reference as the scan).  Spread: N_PER_RT per run type, rest filled from any.
"""
from __future__ import annotations

import os
import sys
import glob
import argparse
import time
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "exotic_events"))
import numpy as np                                      # noqa: E402
import matplotlib                                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402
import event_display as ed                              # noqa: E402  (imports pileup_rates as pr)
from prad2py import dec                                 # noqa: E402

pr = ed.pr
PYDIR = HERE.parent.parent
SEED = 20261115
NS_LIST = (2, 3, 4, 5, 6)
N_EX = 20
N_PER_RT = 5
FILES_PER_RUN = 2
POPS = (("all HyCal", None), ("PbWO4", True), ("lead glass", False))
COLORS = ed.PULSE_COLORS


def veto_mask(k, t, M, W):
    """True where among the sorted fitted onsets some t[i+M] - t[i] <= W (more than M pulses within W ns)."""
    slot = np.arange(t.shape[1])[None, :]
    ts = np.sort(np.where(slot < k[:, None], t, np.float32(np.inf)), axis=1)
    v = np.zeros(k.size, bool)
    with np.errstate(invalid="ignore"):
        for i in range(t.shape[1] - M):
            v |= (k >= i + M + 1) & (ts[:, i + M] - ts[:, i] <= W)
    return v


def tables(scan_dir, veto=None):
    """Per population and run type: counts n_1..n_6, sat, bfit (excluded-class counts inside the counts).
    veto = (M, W): waveforms with more than M onsets within W ns are dropped from every count; their numbers go in 'vn'."""
    res = {p: {rt: dict(n=np.zeros(pr.KMAX + 1, np.int64), vn=np.zeros(pr.KMAX + 1, np.int64), sat=0, bfit=0, runs=[]) for rt in pr.RUN_TYPES} for p, _ in POPS}
    for fn in sorted(Path(scan_dir).glob("run_*.npz")):
        z = np.load(fn)
        rt, run = str(z["run_type"]), int(z["run"])
        k, mi, chi, ncl = z["ch_k"], z["ch_midx"], z["ch_chi2"], z["ch_nclip"]
        pw = z["mod_pwo4"][mi]
        vz = veto_mask(k.astype(np.int64), z["ch_t"], *veto) if veto else np.zeros(k.size, bool)
        for p, sel in POPS:
            m0 = np.ones(k.size, bool) if sel is None else (pw == sel)
            m = m0 & ~vz
            r = res[p][rt]
            r["vn"] += np.bincount(k[m0 & vz], minlength=pr.KMAX + 1)[:pr.KMAX + 1]
            r["n"] += np.bincount(k[m], minlength=pr.KMAX + 1)[:pr.KMAX + 1]
            r["sat"] += int((m & (ncl >= pr.SAT_NCLIP)).sum())
            r["bfit"] += int((m & (chi > pr.CHI2_MAX)).sum())
            if p == POPS[0][0]:
                r["runs"].append(run)
    return res


def rate_rows(res):
    rows = []
    for p, _ in POPS:
        for rt in pr.RUN_TYPES:
            r = res[p][rt]
            n, tot = r["n"], int(r["n"][1:].sum())
            pulses = float((np.arange(pr.KMAX + 1) * n).sum())
            for N in range(1, pr.KMAX + 1):
                f = n[N] / tot if tot else np.nan
                rows.append(dict(pop=p, rt=rt, N=N, count=int(n[N]), denom=tot, rate=f, err=np.sqrt(f * (1 - f) / tot) if tot else np.nan,
                                 pw=N * n[N] / pulses if pulses else np.nan, pw_err=np.nan))
    return rows


def fmt_rate(f, e):
    if not np.isfinite(f):
        return "-"
    p, q = 100 * f, 100 * e
    return f"{p:.3f} +- {q:.3f}" if p >= 1 else f"{p:.4f} +- {q:.4f}"


def table_page(pdf, res, rows, pop, note):
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0.02, 0.2, 0.96, 0.62])
    ax.axis("off")
    hdr = ["run type", "runs", "n_all"] + [f"R_{N} (%)" for N in range(1, 7)]
    cells = []
    for rt in pr.RUN_TYPES:
        sub = {r["N"]: r for r in rows if r["pop"] == pop and r["rt"] == rt}
        cells.append([rt, " ".join(map(str, sorted(res[POPS[0][0]][rt]["runs"]))), f"{sub[1]['denom']:,}"] + [fmt_rate(sub[N]["rate"], sub[N]["err"]) for N in range(1, 7)])
    cells2 = []
    for rt in pr.RUN_TYPES:
        sub = {r["N"]: r for r in rows if r["pop"] == pop and r["rt"] == rt}
        cells2.append([rt, "counts", f"sat {res[pop][rt]['sat']:,}; chi2>1000 {res[pop][rt]['bfit']:,}"] + [f"{sub[N]['count']:,}" for N in range(1, 7)])
    cells3 = []
    for rt in pr.RUN_TYPES:
        sub = {r["N"]: r for r in rows if r["pop"] == pop and r["rt"] == rt}
        cells3.append([rt, "pulse-weighted", ""] + [f"{100 * sub[N]['pw']:.3g}" for N in range(1, 7)])
    allc = cells + cells2 + cells3
    t = ax.table(cellText=allc, colLabels=hdr, loc="upper center", cellLoc="center", colWidths=[0.11, 0.14, 0.15] + [0.1, 0.1, 0.1, 0.1, 0.1, 0.1])
    t.auto_set_font_size(False)
    t.set_fontsize(9.5)
    t.scale(1, 1.9)
    for (i, j), c in t.get_celld().items():
        if i == 0:
            c.set_text_props(fontweight="bold")
        if 1 <= i <= 4:
            c.set_facecolor("#eef3fa")
        elif 9 <= i:
            c.set_facecolor("#f6f2e8")
    fig.suptitle(f"Fraction of waveforms fitted with N pulses: {pop}", fontsize=18, fontweight="bold")
    fig.text(0.03, 0.03, textwrap.fill(note, 190), fontsize=10)
    fig.text(0.03, 0.855, "Blue rows: R_N = n_N / n_all (waveform-count rate, primary) +- binomial error, in %.  White rows: n_N.  Yellow rows: pulse-weighted fraction N n_N / sum_k k n_k (%).", fontsize=10)
    pdf.savefig(fig)
    plt.close(fig)


def veto_title_page(pdf, res, M, W, vex):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.04, 0.93, f"Pile-up multiplicity with waveform veto: no more than {M} pulses within {W:g} ns", fontsize=22, fontweight="bold")
    lines = [f"Veto: a waveform is removed (numerator and denominator) if any window of {W:g} ns holds more than {M} fitted onsets,",
             f"      i.e. among its sorted t0, t0[i+{M}] - t0[i] <= {W:g} ns for some i.  Everything else as in pileup_multiplicity_20261115.pdf.", "",
             "Vetoed waveforms (all HyCal), per run type and N (they are not in the tables' n_N or n_all):", ""]
    lines.append(f"{'run type':14s}" + "".join(f"{'N=' + str(N):>10s}" for N in range(1, pr.KMAX + 1)) + f"{'total':>11s}{'kept':>13s}{'% vetoed':>10s}")
    for rt in pr.RUN_TYPES:
        r = res[POPS[0][0]][rt]
        v, kept = int(r["vn"].sum()), int(r["n"].sum())
        lines.append(f"{rt:14s}" + "".join(f"{int(r['vn'][N]):10,d}" for N in range(1, pr.KMAX + 1)) + f"{v:11,d}{kept:13,d}{100 * v / max(v + kept, 1):10.3f}")
    lines += ["", "Tables: blue rows = new R_N (veto applied); green rows = old R_N without veto; white = new counts; orange = vetoed counts.",
              "Examples: the 100 refitted examples of the original PDF are reused unchanged (no replacements).  Vetoed ones are listed on the next page and",
              f"marked VETOED in red on the example pages.  Vetoed examples per N: " + ", ".join(f"N={N}: {sum(vex[N])}/{len(vex[N])}" for N in NS_LIST)]
    fig.text(0.04, 0.87, "\n".join(lines), fontsize=11.5, va="top", family="monospace", linespacing=1.4)
    pdf.savefig(fig)
    plt.close(fig)


def veto_table_page(pdf, res, rows, res0, rows0, pop, M, W):
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0.02, 0.1, 0.96, 0.72])
    ax.axis("off")
    hdr = ["run type", "n_all", ""] + [f"N={N}" for N in range(1, 7)]
    cells, kinds = [], []
    for kind, rs, rr in (("new R_N (%)", res, rows), ("old R_N (%)", res0, rows0), ("new count", res, rows), ("vetoed count", res, rows)):
        for rt in pr.RUN_TYPES:
            sub = {r["N"]: r for r in rr if r["pop"] == pop and r["rt"] == rt}
            if kind == "vetoed count":
                body = [f"{int(res[pop][rt]['vn'][N]):,}" for N in range(1, 7)]
            elif kind == "new count":
                body = [f"{sub[N]['count']:,}" for N in range(1, 7)]
            else:
                body = [fmt_rate(sub[N]["rate"], sub[N]["err"]) for N in range(1, 7)]
            cells.append([rt, f"{sub[1]['denom']:,}", kind] + body)
            kinds.append(kind)
    t = ax.table(cellText=cells, colLabels=hdr, loc="upper center", cellLoc="center", colWidths=[0.11, 0.11, 0.12] + [0.11] * 6)
    t.auto_set_font_size(False)
    t.set_fontsize(9.5)
    t.scale(1, 1.55)
    col = {"new R_N (%)": "#eef3fa", "old R_N (%)": "#e8f3e8", "new count": "#ffffff", "vetoed count": "#fbeee0"}
    for (i, j), c in t.get_celld().items():
        if i == 0:
            c.set_text_props(fontweight="bold")
        else:
            c.set_facecolor(col[kinds[i - 1]])
    fig.suptitle(f"Veto: no more than {M} pulses within {W:g} ns -- {pop}", fontsize=18, fontweight="bold")
    fig.text(0.03, 0.03, "New: veto applied, n_all excludes vetoed waveforms.  Old: no veto (original table).  Errors binomial, statistics only.", fontsize=10)
    pdf.savefig(fig)
    plt.close(fig)


def veto_examples(pdf, orig_pdf, ex_tsv, M, W):
    """List of vetoed examples + original example pages (copied, VETOED overlay).  Returns per-N veto flags."""
    from pypdf import PdfReader
    rows = [l.rstrip("\n").split("\t") for l in open(ex_tsv)][1:]
    flags = {N: [] for N in NS_LIST}
    for r in rows:
        tt = np.array([float(x.split("@")[1]) for x in r[6].split(";")], np.float32)[None, :]
        pad = np.full((1, pr.KMAX), np.inf, np.float32)
        pad[0, :tt.shape[1]] = tt[0]
        flags[int(r[0])].append(bool(veto_mask(np.array([tt.shape[1]]), pad, M, W)[0]))
    return flags, rows, PdfReader(orig_pdf)


def veto_list_page(pdf, flags, rows, M, W):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.04, 0.93, f"Examples removed by the veto ({M} pulses / {W:g} ns)", fontsize=20, fontweight="bold")
    lst, i = [], 0
    for r in rows:
        N = int(r[0])
        pos = sum(1 for q in rows[:i] if q[0] == r[0])
        if flags[N][pos]:
            lst.append(f"N={N} panel {pos + 1:2d} (row {pos // 5 + 1}, col {pos % 5 + 1}): run {r[1]} ev {r[2]} {r[3]} {r[4]}  {r[6]}")
        i += 1
    txt = "\n".join(textwrap.shorten(x, 150) for x in lst) if lst else "(none)"
    fig.text(0.04, 0.88, f"{len(lst)} of {len(rows)} examples vetoed.  Panel order = left to right, top to bottom on the example pages.\n\n" + txt, fontsize=9, va="top", family="monospace")
    pdf.savefig(fig)
    plt.close(fig)


def title_page(pdf, res, tot_files):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.04, 0.93, "Pile-up multiplicity: fraction of waveforms with N fitted pulses", fontsize=22, fontweight="bold")
    txt = (
        "Definitions\n"
        "  R_N = n_N / n_all : n_N = waveforms whose matched-filter fit finds exactly N pulses (N = 1..6); n_all = waveforms with >= 1 fitted pulse\n"
        "        (kept channels). PRIMARY. Error: binomial sqrt(R(1-R)/n_all).\n"
        "  P_N = N n_N / sum_k k n_k : fraction of all fitted pulses that sit in N-pulse waveforms (secondary, pulse-weighted).\n\n"
        "Selection (pileup_rates.py defaults, scan outputs of 2026-11-11, all EVIO splits of each run)\n"
        f"  HyCal channels with 100 samples, not flagged bad (sample == 0 or > {pr.UNDERSHOOT_ADC:g} ADC below baseline); fitted when raw max - baseline >= {pr.PRECUT_ADC:g} ADC;\n"
        f"  kept when the largest fitted amplitude >= {pr.KEEP_ADC:g} ADC. Trigger gate of the run type, LMS-bit veto; no lone-channel / array-wide-event veto.\n"
        f"  Saturated waveforms (nclip >= {pr.SAT_NCLIP}) and chi2/dof > {pr.CHI2_MAX:g} waveforms are INCLUDED in the counts (their numbers are in the table);\n"
        "  they are excluded from the example pages.\n"
        f"  KMAX = {pr.KMAX}: the fit never returns more than {pr.KMAX} pulses, so N = 6 includes waveforms that hit the cap.\n"
        "  Table 1: all HyCal (PbWO4 + lead glass); tables 2, 3: PbWO4 and lead glass separately.\n\n"
        "Run types and runs (run_pileup_scan.sh)\n  " + "\n  ".join(f"{rt}: " + " ".join(map(str, sorted(res[POPS[0][0]][rt]['runs']))) for rt in pr.RUN_TYPES) + "\n\n"
        "Examples (pages after the tables)\n"
        f"  {N_EX} random waveforms per N = 2..6 (seed {SEED}), {N_PER_RT} per run type where available, drawn from {FILES_PER_RUN} random EVIO splits per run (limits the\n"
        "  per-file noise-reference reads); each event is decoded again and refitted with the same fit and noise reference as the scan; panels with a\n"
        "  different k than the scan are replaced. Black dots: raw samples (red x: overflow >= 4095); grey: fitted pedestal; red: total model; dashed: pulses.\n"
        "  ")
    fig.text(0.04, 0.88, txt, fontsize=11.5, va="top", family="monospace", linespacing=1.35)
    pdf.savefig(fig)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------------------ examples
def split_ranges(ctx, files):
    firsts = [ed.first_event_number(ctx, f) for f in files]
    return firsts


def candidates(scan_dir, evio_dir, rng):
    """{N: list of dict(run, rt, ev, midx, file)} candidates (non-saturated, chi2 ok) in FILES_PER_RUN random splits per run."""
    cand = {N: [] for N in NS_LIST}
    chosen = {}
    for fn in sorted(Path(scan_dir).glob("run_*.npz")):
        z = np.load(fn)
        rt, run = str(z["run_type"]), int(z["run"])
        files = sorted(glob.glob(str(Path(evio_dir) / f"prad_{run:06d}.evio.*")))
        ctx = ed.Ctx(run)
        firsts = [ed.first_event_number(ctx, f) for f in files]
        evn = z["ev_number"][z["ch_ev"]]
        fidx = np.searchsorted(np.array([f if f is not None else 2**62 for f in firsts]), evn, side="right") - 1
        ok = (z["ch_nclip"] < pr.SAT_NCLIP) & ~(z["ch_chi2"] > pr.CHI2_MAX)
        k = z["ch_k"]
        present = [i for i in np.unique(fidx) if i >= 0 and (ok & (fidx == i) & (k >= 2)).sum() > 0]
        pick = sorted(rng.choice(present, size=min(FILES_PER_RUN, len(present)), replace=False).tolist())
        chosen[run] = [files[i] for i in pick]
        print(f"[cand] run {run} {rt}: files {[os.path.basename(files[i]) for i in pick]}", flush=True)
        for N in NS_LIST:
            for i in pick:
                idx = np.nonzero(ok & (fidx == i) & (k == N))[0]
                for j in idx[rng.permutation(idx.size)[:200]]:        # at most 200 random per (run, file, N)
                    cand[N].append(dict(run=run, rt=rt, ev=int(evn[j]), midx=int(z["ch_midx"][j]), file=files[i]))
    return cand


def choose(cand, rng):
    """N_EX per N: N_PER_RT per run type first, remainder from the others; returns {N: ordered pool (more than N_EX, for replacement)}."""
    out = {}
    for N in NS_LIST:
        pool = cand[N]
        order = rng.permutation(len(pool))
        by = {rt: [pool[i] for i in order if pool[i]["rt"] == rt] for rt in pr.RUN_TYPES}
        first = []
        for rt in pr.RUN_TYPES:
            first += by[rt][:N_PER_RT]
        used = {id(x) for x in first}
        rest = [pool[i] for i in order if id(pool[i]) not in used]
        out[N] = (first, rest)
    return out


def refit(chosen_pool, evio_dir):
    """Refit the events needed; returns {(N, idx): example dict}.  Replacement from the rest pool when k differs."""
    # one noise reference per file, one decode pass per file for all wanted event numbers
    wanted = {}
    for N, (first, rest) in chosen_pool.items():
        for x in first + rest[:40]:
            wanted.setdefault(x["file"], {}).setdefault(x["ev"], []).append(x)
    ctxs, refs, evs = {}, {}, {}
    for f, d in sorted(wanted.items()):
        run = next(iter(d.values()))[0]["run"]
        ctx = ctxs.setdefault(run, ed.Ctx(run))
        t0 = time.time()
        refs[f] = ed.noise_reference(ctx, f)[1:]
        print(f"[ref] {os.path.basename(f)} {time.time() - t0:.0f}s, {len(d)} candidate events", flush=True)
    return ctxs, refs, wanted


def decode_wanted(ctxs, refs, wanted, need):
    """Fit and store, for the events in need = {file: set(ev)}, the channel results (only the needed channels)."""
    store = {}
    for f, evset in need.items():
        run = next(iter(wanted[f].values()))[0]["run"]
        ctx = ctxs[run]
        b_ref, sig_ref = refs[f]
        left = set(evset)
        ch = ctx.open(f)
        while left and ch.read() == dec.Status.success:
            if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
                continue
            for i in range(ch.get_n_events()):
                ch.select_event(i)
                en = int(ch.info().event_number)
                if en in left:
                    d = ch.decode_event(i, with_ssp=False)
                    if d["ok"]:
                        store[(f, en)] = ed.fit_event(ctx, d["event"], b_ref, sig_ref)
                    left.discard(en)
        ch.close()
    return store


def make_example(ctx, r, x, N):
    i = r["row"].get(x["midx"])
    if i is None or int(r["k"][i]) != N:
        return None
    return dict(run=x["run"], rt=x["rt"], ev=x["ev"], midx=x["midx"], name=ctx.names[x["midx"]], pwo4=bool(ctx.pwo4[x["midx"]]), N=N,
                Y=r["Y"][i].astype(float), ped=float(r["ped"][i]), a=r["a"][i][:N].copy(), t=r["t"][i][:N].copy(), comp=r["comp"][i][:N].copy(), chi2=float(r["chi2"][i]),
                nclip=int(r["nclip"][i]))


def draw_panel(ax, e):
    y = e["Y"]
    ov = y >= 4095
    ax.plot(ed.TT[~ov], y[~ov], "k.", ms=4)
    if ov.any():
        ax.plot(ed.TT[ov], y[ov], "x", color="tab:red", ms=6, mew=1.5)
    ax.axhline(e["ped"], color="0.55", lw=1.0)
    ax.plot(ed.TT, e["ped"] + e["comp"].sum(0), "r-", lw=1.4)
    for j in range(e["N"]):
        ax.plot(ed.TT, e["ped"] + e["comp"][j], "--", color=COLORS[j], lw=1.2)
    pulses = "; ".join(f"{a:.0f}@{t:.0f}" for a, t in zip(e["a"], e["t"]))
    ttl = f"run {e['run']} ev {e['ev']} {e['name']} ({'PbWO4' if e['pwo4'] else 'LG'}) {e['rt']}\nchi2/dof {e['chi2']:.2f}   (A ADC @ t0 ns)\n" + textwrap.fill(pulses, 46)
    ax.set_title(ttl, fontsize=7.5, loc="left", pad=3)
    ax.tick_params(labelsize=7)
    ax.set_xlim(0, ed.DT * ed.NS)


def examples_pages(pdf, ex):
    for N in NS_LIST:
        lst = ex[N]
        fig, axs = plt.subplots(4, 5, figsize=(22, 17))
        for ax in axs.flat:
            ax.axis("off")
        for ax, e in zip(axs.flat, lst):
            ax.axis("on")
            draw_panel(ax, e)
        fig.suptitle(f"N = {N} fitted pulses: {len(lst)} random examples (requested {N_EX})", fontsize=18, fontweight="bold")
        fig.text(0.01, 0.003, "black dots: raw samples (red x: overflow >= 4095); grey line: fitted pedestal; red: total model; dashed: pulse components in time order "
                 "(blue, orange, green, purple, brown, pink). x axis: t (ns), y: ADC.", fontsize=11)
        fig.tight_layout(rect=(0, 0.015, 1, 0.96), h_pad=1.2)
        pdf.savefig(fig)
        plt.close(fig)


def veto_main(a):
    from pypdf import PdfReader, PdfWriter
    out, M, W = Path(a.out), a.veto_m, a.veto_w
    t0 = time.time()
    res0 = tables(a.scan_dir)
    rows0 = rate_rows(res0)
    res = tables(a.scan_dir, (M, W))
    rows = rate_rows(res)
    with open(out / f"{a.stem}_table.tsv", "w") as f:
        f.write(f"# veto: more than {M} fitted onsets within {W:g} ns removes the waveform from counts and denominator\n")
        f.write("population\trun_type\tN\tcount\tdenominator\trate_waveform\trate_err_binom\tpulse_weighted_fraction\tn_vetoed\told_count\told_denominator\told_rate_waveform\n")
        for r, r0 in zip(rows, rows0):
            f.write(f"{r['pop']}\t{r['rt']}\t{r['N']}\t{r['count']}\t{r['denom']}\t{r['rate']:.6g}\t{r['err']:.3g}\t{r['pw']:.6g}\t{int(res[r['pop']][r['rt']]['vn'][r['N']])}\t{r0['count']}\t{r0['denom']}\t{r0['rate']:.6g}\n")
    print(f"[table] {time.time() - t0:.0f}s", flush=True)
    orig_pdf, ex_tsv = out / f"{a.orig_stem}.pdf", out / f"{a.orig_stem}_examples.tsv"
    flags, rows_ex, rd = veto_examples(None, orig_pdf, ex_tsv, M, W)
    tmp = out / f".{a.stem}_front.pdf"
    with PdfPages(tmp) as pdf:
        veto_title_page(pdf, res, M, W, flags)
        for p, _ in POPS:
            veto_table_page(pdf, res, rows, res0, rows0, p, M, W)
        veto_list_page(pdf, flags, rows_ex, M, W)
    # overlay pages: VETOED stamps over the original example pages (panel grid 4 x 5, tight_layout axes boxes taken from a re-created layout)
    from matplotlib.transforms import Bbox
    w = PdfWriter()
    for pg in PdfReader(tmp).pages:
        w.add_page(pg)
    for N, pg in zip(NS_LIST, rd.pages[4:]):
        fl = flags[N]
        if any(fl):
            ov = out / f".{a.stem}_ov{N}.pdf"
            fig, axs = plt.subplots(4, 5, figsize=(22, 17))
            for ax in axs.flat:
                ax.axis("off")
            fig.tight_layout(rect=(0, 0.015, 1, 0.96), h_pad=1.2)
            fig.patch.set_alpha(0)
            for ax, v in zip(axs.flat, fl):
                if v:
                    bb = ax.get_position()
                    fig.text((bb.x0 + bb.x1) / 2, (bb.y0 + bb.y1) / 2 - 0.01, "VETOED", color="red", alpha=0.45, fontsize=34, fontweight="bold", ha="center", va="center", rotation=25)
            fig.savefig(ov, transparent=True)
            plt.close(fig)
            pg.merge_page(PdfReader(ov).pages[0])
            ov.unlink()
        w.add_page(pg)
    with open(out / f"{a.stem}.pdf", "wb") as f:
        w.write(f)
    tmp.unlink()
    print(f"[done] {time.time() - t0:.0f}s", flush=True)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-dir", default=str(PYDIR / "output" / "pileup_rates" / "scan"))
    ap.add_argument("--evio-dir", default="/Users/L00338853/work/PRad/data/evio")
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--stem", default="pileup_multiplicity_20261115")
    ap.add_argument("--no-examples", action="store_true")
    ap.add_argument("--veto-m", type=int, default=0, help="veto waveforms with more than M fitted onsets within --veto-w ns (default off)")
    ap.add_argument("--veto-w", type=float, default=0.0)
    ap.add_argument("--orig-stem", default="pileup_multiplicity_20261115", help="veto mode: stem of the original PDF/examples TSV whose example pages are reused")
    a = ap.parse_args()
    if a.veto_m:
        return veto_main(a)
    out = Path(a.out)
    t0 = time.time()
    res = tables(a.scan_dir)
    rows = rate_rows(res)
    with open(out / f"{a.stem}_table.tsv", "w") as f:
        f.write("population\trun_type\tN\tcount\tdenominator\trate_waveform\trate_err_binom\tpulse_weighted_fraction\tn_saturated_in_denominator\tn_chi2gt1000_in_denominator\n")
        for r in rows:
            q = res[r["pop"]][r["rt"]]
            f.write(f"{r['pop']}\t{r['rt']}\t{r['N']}\t{r['count']}\t{r['denom']}\t{r['rate']:.6g}\t{r['err']:.3g}\t{r['pw']:.6g}\t{q['sat']}\t{q['bfit']}\n")
    print(f"[table] {time.time() - t0:.0f}s", flush=True)
    ex = {N: [] for N in NS_LIST}
    if not a.no_examples:
        rng = np.random.default_rng(SEED)
        cand = candidates(a.scan_dir, a.evio_dir, rng)
        print({N: len(v) for N, v in cand.items()}, flush=True)
        pool = choose(cand, rng)
        ctxs, refs, wanted = refit(pool, a.evio_dir)
        need = {}
        for N, (first, rest) in pool.items():
            for x in first + rest[:40]:
                need.setdefault(x["file"], set()).add(x["ev"])
        store = decode_wanted(ctxs, refs, wanted, need)
        for N, (first, rest) in pool.items():
            for x in first + rest:
                if len(ex[N]) >= N_EX:
                    break
                r = store.get((x["file"], x["ev"]))
                e = make_example(ctxs[x["run"]], r, x, N) if r is not None else None
                if e is not None:
                    ex[N].append(e)
            print(f"[ex] N={N}: {len(ex[N])} of {N_EX}", flush=True)
        with open(out / f"{a.stem}_examples.tsv", "w") as f:
            f.write("N\trun\tevent\tmodule\trun_type\tchi2_dof\tpulses_A_t0\n")
            for N in NS_LIST:
                for e in ex[N]:
                    f.write(f"{N}\t{e['run']}\t{e['ev']}\t{e['name']}\t{e['rt']}\t{e['chi2']:.3f}\t" + ";".join(f"{x:.1f}@{y:.1f}" for x, y in zip(e["a"], e["t"])) + "\n")
    note = ("Counts from the full scans of the listed runs (all EVIO splits). Saturated (nclip >= 3) and chi2/dof > 1000 waveforms are included; N = 6 includes waveforms that "
            "hit the KMAX = 6 cap. Errors are binomial and cover statistics only.")
    with PdfPages(out / f"{a.stem}.pdf") as pdf:
        title_page(pdf, res, 0)
        for p, _ in POPS:
            table_page(pdf, res, rows, p, note)
        if not a.no_examples:
            examples_pages(pdf, ex)
    print(f"[done] {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
