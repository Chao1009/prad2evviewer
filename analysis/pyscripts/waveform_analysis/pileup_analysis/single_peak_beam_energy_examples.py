#!/usr/bin/env python3
"""single_peak_beam_energy_examples.py -- 20 single-pulse waveforms per beam type with E_main within beam energy +-10 % and chi2/dof > 20 (default MF fit).

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/single_peak_beam_energy_examples.py [--sample-dir <dir>] [--evio-dir <dir>] [--out <dir>]

Candidates come from the OLD (all switches off) bounded scans of run_monoedge_sample.sh: kept channels with k == 1 and E = A x gain
(module_gain of pileup_rates) within +-WIN of the run's beam energy and chi2/dof > CHI2_MIN.  The first N_EX per run type (run, event order) are decoded again,
refitted with pr.MONO_EDGE off and verified against the stored k and chi2/dof; a failed check is replaced by the next candidate.
Counts: kept channels = fitted channels with A_main >= KEEP_ADC (the scan's examined population), then k == 1, then in the energy window, then chi2/dof > CHI2_MIN.
-> --stem (default single_peak_beam_energy_chi2gt20_3sigma_20261116) .pdf / .tsv; --window 0.10 --stem single_peak_beam_energy_chi2gt20_20261116 reproduces the first version
"""
from __future__ import annotations

import os
import sys
import time
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
os.environ["PILEUP_SHAPE_SCALE"] = "1"
import numpy as np                                      # noqa: E402
import matplotlib                                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402
import pileup_monoedge as pmo                           # noqa: E402  (sets PILEUP_MONO_EDGE=1 at import; toggled off below)
from prad2py import dec                                 # noqa: E402

pm, pr, ed = pmo.pm, pmo.pr, pmo.ed
CHI2_MIN, N_EX, SPARE = 20.0, 20, 15
# energy window |E/E_beam - 1| <= WIN per run type: 3 sigma of the HyCal resolution (control-v3 fit 0.90 + 2.05/sqrt(E) + 1.35/E in quadrature, % with E in GeV; Ta_X17: measured 1.99 %)
WIN3 = {"0p73GeV_H2": 0.095, "2p24GeV_H2": 0.052, "3p49GeV_H2": 0.044, "2p24GeV_Ta_X17": 0.060}
FULL_FIRST = {"0p73GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"}   # run types taken from the all-split scans (as the first +-10 % version, where the 2-file sample fell short)
TYPES = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]
STEM = "single_peak_beam_energy_chi2gt20_3sigma_20261116"
KEYS = ("ch_ev", "ch_midx", "ch_k", "ch_a", "ch_t", "ch_chi2", "ch_nclip")


def scan_type(scan_dir, evio_dir, rt, half):
    """Counts and ordered candidate list of one run type."""
    cnt = dict(kept=0, k1=0, win=0, hi=0, runs=[], A=[])
    cands = []
    for fn in sorted(Path(scan_dir).glob("run_*.npz")):
        z = np.load(fn)
        if str(z["run_type"]) != rt:
            continue
        run = int(z["run"])
        cfg = pr.run_config(run)
        names = [str(n) for n in z["mod_name"]]
        gain = pr.module_gain(cfg, rt, names)
        k, a, chi2, midx = z["ch_k"], z["ch_a"][:, 0], z["ch_chi2"], z["ch_midx"]
        # amplitude of the single pulse; E from the module calibration
        e = a * gain[midx]
        k1 = k == 1
        win = k1 & (np.abs(e / cfg["beam"] - 1) <= half)
        hi = win & (chi2 > CHI2_MIN)
        cnt["kept"] += k.size
        cnt["k1"] += int(k1.sum())
        cnt["win"] += int(win.sum())
        cnt["hi"] += int(hi.sum())
        cnt["A"] += [float(x) for x in a[hi]]
        cnt["runs"].append(run)
        files = [str(Path(evio_dir) / f"prad_{run:06d}.evio.{int(s):05d}") for s in z["splits"]]
        evn = z["ev_number"]
        ctx = ed.Ctx(run)
        firsts = np.array([ed.first_event_number(ctx, f) for f in files], dtype=np.int64)
        idx = np.nonzero(hi)[0]
        idx = idx[np.lexsort((midx[idx], z["ch_ev"][idx]))]
        for j in idx[:N_EX + SPARE]:
            en = int(evn[z["ch_ev"][j]])
            fi = int(np.clip(np.searchsorted(firsts, en, side="right") - 1, 0, len(files) - 1))
            cands.append(dict(run=run, rt=rt, ev=en, midx=int(midx[j]), file=files[fi], split=fi, ref=(z["mod_b_ref"][fi], z["mod_sig_ref"][fi]), chi2=float(chi2[j]), k=1, beam=cfg["beam"]))
        print(f"[scan] run {run} {rt}: kept {k.size}, k=1 {int(k1.sum())}, window {int(win.sum())}, chi2>{CHI2_MIN:g} {int(hi.sum())}", flush=True)
    return cnt, cands


def refit(cands):
    """Decode each needed event once, default fit.  {(file, ev): result}, ctxs."""
    need, ref_of = {}, {}
    for c in cands:
        need.setdefault((c["run"], c["file"], c["split"]), set()).add(c["ev"])
        ref_of[(c["run"], c["file"])] = c["ref"]
    ctxs, store = {}, {}
    pr.MONO_EDGE = False
    for (run, f, si), evs in sorted(need.items()):
        ctx = ctxs.setdefault(run, ed.Ctx(run))
        b_ref, sig_ref = ref_of[(run, f)]
        left = set(evs)
        ch = ctx.open(f)
        t0 = time.time()
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
        print(f"[fit] run {run} {Path(f).name}: {len(evs)} events, {time.time() - t0:.0f} s", flush=True)
    return store, ctxs


def build(c, r, ctx):
    i = r["row"].get(c["midx"])
    if i is None or int(r["k"][i]) != 1 or not np.isclose(float(r["chi2"][i]), c["chi2"], rtol=2e-3, atol=1e-3):
        return None
    mi = c["midx"]
    a, t = float(r["a"][i][0]), float(r["t"][i][0])
    _, sig_ref = c["ref"]
    sig = np.array([sig_ref[mi]], np.float32)
    m = ctx.bank.row[[mi]]
    var = ctx.mf.variance(m, np.array([[t]], np.float32), np.array([[a]], np.float32), np.array([[True]]), sig, ed.NS)[0]
    ov = r["Y"][i] >= pr.OVERFLOW
    dof = int((pr.MF.wmask(np.array([t], np.float32), ed.NS)[0] & ~ov).sum()) - 3      # as pileup_rates.fit: window samples - 2 per pulse - 1
    return dict(c, name=str(ctx.names[mi]), pwo4=bool(ctx.pwo4[mi]), Y=r["Y"][i].astype(float), ped=float(r["ped"][i]), a=a, t=t,
                comp=r["comp"][i][0].copy(), E=a * float(ctx.gain[mi]), dof=dof, sd=np.sqrt(var), nclip=int(r["nclip"][i]))


def page(pdf, rt, chunk, p0, n, win):
    fig = plt.figure(figsize=(16, 12.5))
    fig.text(0.04, 0.97, f"{rt}: single-pulse fit, E within beam energy +-{100 * win:.1f} %, chi2/dof > {CHI2_MIN:g}", fontsize=15, fontweight="bold")
    fig.text(0.04, 0.948, f"examples {p0 + 1}-{p0 + len(chunk)} of {n}.  Black dots: raw samples (red x: overflow); red: MF fit; grey: pedestal; green dotted: onset t0.  Lower strip: data - fit with +-sqrt(V) band.", fontsize=9)
    outer = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.12, left=0.05, right=0.99, top=0.9, bottom=0.04)
    for n_, it in enumerate(chunk):
        gs = outer[n_ // 2, n_ % 2].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.08)
        ax, axr = fig.add_subplot(gs[0]), fig.add_subplot(gs[1])
        y = it["Y"]
        ov = y >= pr.OVERFLOW
        ax.plot(ed.TT[~ov], y[~ov], "k.", ms=4)
        if ov.any():
            ax.plot(ed.TT[ov], y[ov], "x", color="tab:red", ms=6, mew=1.5)
        ax.axhline(it["ped"], color="0.55", lw=1.0)
        model = it["ped"] + it["comp"]
        ax.plot(ed.TT, model, "r-", lw=1.4)
        ax.axvline(it["t"], color="tab:green", ls=":", lw=1.6)
        ax.set_title(f"run {it['run']} ev {it['ev']} {it['name']} ({'PbWO4' if it['pwo4'] else 'LG'})  A {it['a']:.0f} ADC  E {it['E']:.0f} MeV  t0 {it['t']:.1f} ns  "
                     f"chi2/dof {it['chi2']:.1f}" + (f"  nclip {it['nclip']}" if it["nclip"] else ""), fontsize=8.5, loc="left", pad=3)
        ax.set_xlim(0, ed.DT * ed.NS)
        ax.tick_params(labelsize=7, labelbottom=False)
        res = y - model
        axr.fill_between(ed.TT, -it["sd"], it["sd"], color="0.8", lw=0)
        axr.plot(ed.TT[~ov], res[~ov], "k.-", lw=0.6, ms=3)
        axr.axhline(0, color="r", lw=0.8)
        axr.axvline(it["t"], color="tab:green", ls=":", lw=1.2)
        axr.set_xlim(0, ed.DT * ed.NS)
        axr.tick_params(labelsize=7)
        axr.set_ylabel("data - fit", fontsize=7)
        axr.set_xlabel("t (ns)", fontsize=8)
    pdf.savefig(fig)
    plt.close(fig)


def title_page(pdf, cnts, got, elapsed):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.03, 0.94, "Single-pulse waveforms at beam energy with chi2/dof > 20 (default MF fit)", fontsize=20, fontweight="bold")
    fig.text(0.03, 0.88, "Population: kept channels of the bounded 2-file-per-run sample scans (run_monoedge_sample.sh, old = all switches off; A_main >= %g ADC).\n"
             "Selection: k == 1 fitted pulse, E = A x gain (module calibration, pileup_rates.module_gain) within the window of the table (|E/E_beam - 1| <=), chi2/dof > %g.\n"
             "Scan source per row: 2p24GeV_H2 = 2-file sample; the other types = all-split default scans of the same runs.\n"
             "The first %d per run type (run, event order) are refitted from the raw waveform and shown; the refit must reproduce the stored k and chi2/dof." % (pr.KEEP_ADC, CHI2_MIN, N_EX),
             fontsize=11, va="top", family="monospace", linespacing=1.5)
    hdr = ["beam type", "scan", "window", "runs", "kept channels", "k = 1", "k = 1, E in window", "chi2/dof > 20", "fraction", "A min/med/max (ADC)", "shown"]
    rows = []
    for rt in TYPES:
        c = cnts[rt]
        rows.append([rt, c["src"], f"+-{100 * c['win_half']:.1f} %", f"{len(c['runs'])} runs ({c['runs'][0]}-{c['runs'][-1]})", f"{c['kept']:,}", f"{c['k1']:,}", f"{c['win']:,}", f"{c['hi']:,}",
                     f"{100 * c['hi'] / max(c['win'], 1):.3f} %", f"{min(c['A']):.0f} / {np.median(c['A']):.0f} / {max(c['A']):.0f}" if c["A"] else "-", str(got[rt])])
    ax = fig.add_axes([0.03, 0.3, 0.94, 0.35])
    ax.axis("off")
    tb = ax.table(cellText=rows, colLabels=hdr, loc="center", cellLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(8.5)
    tb.scale(1, 2.2)
    fig.text(0.03, 0.28, f"fraction = (chi2/dof > 20) / (k = 1 with E in window).  Run time {elapsed:.0f} s.", fontsize=10)
    pdf.savefig(fig)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir", default=str(pmo.PYDIR / "output" / "pileup_rates" / "monoedge_sample"))
    ap.add_argument("--full-scan-dir", default=str(pmo.PYDIR / "output" / "pileup_rates" / "scan"), help="all-split default scans, used for a run type the 2-file sample cannot fill")
    ap.add_argument("--window", default=None, help="half-width of the energy window |E/E_beam - 1|: one number for all run types (0.10 = first version), or "
                    "type=value,... (default per type: 3 sigma, %s)" % ",".join(f"{k}={v}" for k, v in WIN3.items()))
    ap.add_argument("--stem", default="single_peak_beam_energy_chi2gt20_3sigma_20261116", help="output file stem")
    ap.add_argument("--evio-dir", default=pmo.EVIO_DIR)
    ap.add_argument("--out", default=str(HERE))
    a = ap.parse_args()
    global STEM
    STEM = a.stem
    win = dict(WIN3)
    if a.window:
        if "=" in a.window:
            win.update({kv.split("=")[0]: float(kv.split("=")[1]) for kv in a.window.split(",")})
        else:
            win = {rt: float(a.window) for rt in TYPES}
    t0 = time.time()
    cnts, allc = {}, []
    for rt in TYPES:
        full = rt in FULL_FIRST
        cnts[rt], c = scan_type(a.full_scan_dir if full else Path(a.sample_dir) / "old", a.evio_dir, rt, win[rt])
        cnts[rt]["src"] = "all splits" if full else "2-file sample"
        cnts[rt]["win_half"] = win[rt]
        if cnts[rt]["hi"] < N_EX and not full:
            print(f"[scan] {rt}: sample gives {cnts[rt]['hi']} < {N_EX}; using the all-split scans", flush=True)
            cnts[rt]["sample"] = {k: cnts[rt][k] for k in ("kept", "k1", "win", "hi")}
            cnts[rt], c = scan_type(a.full_scan_dir, a.evio_dir, rt, win[rt])
            cnts[rt]["src"] = "all splits (sample fell short)"
            cnts[rt]["win_half"] = win[rt]
        allc += c
    store, ctxs = refit(allc)
    final, nfail = {}, 0
    for rt in TYPES:
        items = []
        for c in [c for c in allc if c["rt"] == rt]:
            if len(items) >= N_EX:
                break
            r = store.get((c["file"], c["ev"]))
            it = build(c, r, ctxs[c["run"]]) if r is not None else None
            if it is None:
                nfail += 1
                print(f"[verify FAIL] {c['run']} ev {c['ev']} m {c['midx']}", flush=True)
                continue
            items.append(it)
        final[rt] = items
    stem = Path(a.out) / STEM
    with PdfPages(str(stem) + ".pdf") as pdf:
        title_page(pdf, cnts, {k: len(v) for k, v in final.items()}, time.time() - t0)
        for rt in TYPES:
            for p0 in range(0, len(final[rt]), 4):
                page(pdf, rt, final[rt][p0:p0 + 4], p0, len(final[rt]), cnts[rt]["win_half"])
    with open(str(stem) + ".tsv", "w") as f:
        f.write("run_type\trun\tevent\tmodule\tA\tE_MeV\tt0\tchi2\tdof\tnclip\n")
        for rt in TYPES:
            for it in final[rt]:
                f.write(f"{rt}\t{it['run']}\t{it['ev']}\t{it['name']}\t{it['a']:.1f}\t{it['E']:.1f}\t{it['t']:.2f}\t{it['chi2'] * it['dof']:.2f}\t{it['dof']}\t{it['nclip']}\n")
    print(f"[done] {sum(len(v) for v in final.values())} examples, {nfail} verify failures, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
