#!/usr/bin/env python3
"""t0_shape_err_sensitivity.py -- do the 50 waveforms of monoedge_fit_comparison keep their extra pulse when the template-shape error is inflated?

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/t0_shape_err_sensitivity.py [--sample-dir <dir>] [--evio-dir <dir>] [--out <dir>]

Parameter: sigma_s = bank shape_sigma (per module, floor SHAPE_SIGMA_MIN = 0.005 of the peak) in the sigma(t) variance
V_j = sig^2 + (SHAPE_SCALE * sum_k a_k sigma_s(t_j - t0_k) sqrt(aband / max(a_k, 20)))^2; knob pileup_rates.SHAPE_SCALE (env PILEUP_SHAPE_SCALE), default 1.0.
The same scale multiplies the shape term of the candidate-peak test.  Each waveform of monoedge_fit_comparison_20261116.tsv is refitted at scale 1, 2, 3
(MONO_EDGE off).  Pass: groups A (H2_small) and C (Ta_small) drop to one pulse, B (H2_large) and D (Ta_large) keep two.
-> T0_shape_err_sensitivity_20261116.pdf / .tsv
"""
from __future__ import annotations

import sys
import time
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import matplotlib                                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402
import monoedge_fit_comparison as mfc                   # noqa: E402
from prad2py import dec                                 # noqa: E402

pmo, pm, pr, ed = mfc.pmo, mfc.pm, mfc.pr, mfc.ed
SCALES = (1.0, 2.0, 3.0)
GROUPS = [("H2_small", "A", 1), ("H2_large", "B", 2), ("Ta_small", "C", 1), ("Ta_large", "D", 2)]   # key, label, expected pulses
GTITLE = {k: t for k, t, *_ in mfc.CATS}
PER = 3


def load_list(tsv):
    rows = []
    with open(tsv) as f:
        hdr = f.readline().split()
        for ln in f:
            d = dict(zip(hdr, ln.split()))
            rows.append(dict(run=int(d["run"]), ev=int(d["event"]), name=d["module"], rt=d["run_type"], cat=d["category"],
                             k_off=int(d["k_off"]), chi2_off=float(d["chi2_off"]), dt=float(d["removed_dt_ns"]), r=float(d["removed_r"])))
    return rows


def locate(sample_dir, evio_dir, rows):
    """Attach file/split/midx to each TSV row (candidate table of monoedge_fit_comparison)."""
    cands = mfc.collect(sample_dir, evio_dir)
    names = {}
    out = []
    for r in rows:
        names.setdefault(r["run"], ed.Ctx(r["run"]).names)
        hit = [c for c in cands if c["run"] == r["run"] and c["ev"] == r["ev"] and names[r["run"]][c["midx"]] == r["name"]]
        if not hit:
            raise SystemExit(f"cannot locate {r['run']} {r['ev']} {r['name']}")
        out.append(dict(r, file=hit[0]["file"], split=hit[0]["split"], midx=hit[0]["midx"]))
    return out


def fit_scales(sample_dir, items):
    refs = {}
    for fn in sorted((Path(sample_dir) / "new").glob("run_*.npz")):
        z = np.load(fn)
        for i in range(len(z["splits"])):
            refs[(int(z["run"]), i)] = (z["mod_b_ref"][i], z["mod_sig_ref"][i])
    need = {}
    for c in items:
        need.setdefault((c["run"], c["file"], c["split"]), set()).add(c["ev"])
    ctxs, store = {}, {}
    pr.MONO_EDGE = False
    try:
        for (run, f, si), evs in sorted(need.items()):
            ctx = ctxs.setdefault(run, ed.Ctx(run))
            b_ref, sig_ref = refs[(run, si)]
            left, ch, t0 = set(evs), ctx.open(f), time.time()
            while left and ch.read() == dec.Status.success:
                if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
                    continue
                for i in range(ch.get_n_events()):
                    ch.select_event(i)
                    en = int(ch.info().event_number)
                    if en in left:
                        d = ch.decode_event(i, with_ssp=False)
                        if d["ok"]:
                            out = {}
                            for s in SCALES:
                                pr.SHAPE_SCALE = s
                                out[s] = ed.fit_event(ctx, d["event"], b_ref, sig_ref)
                            store[(f, en)] = out
                        left.discard(en)
            ch.close()
            print(f"[fit] run {run} {Path(f).name}: {len(evs)} events, {time.time() - t0:.0f} s", flush=True)
    finally:
        pr.SHAPE_SCALE = 1.0
    return store, ctxs


def build(items, store, ctxs):
    res, missing = [], 0
    for c in items:
        r = store.get((c["file"], c["ev"]))
        ex = None if r is None else {s: mfc.make(ctxs[c["run"]], c, r[s], "") for s in SCALES}
        if ex is None or any(v is None for v in ex.values()):
            missing += 1
            print(f"[missing] {c['run']} {c['ev']} {c['name']}", flush=True)
            continue
        res.append(dict(c, ex=ex))
    return res, missing


def panel(ax, axr, e, s):
    pm.draw_panel(ax, e)
    for t in e["t"]:
        ax.axvline(t, color="m", ls=":", lw=1.2)
    ax.set_title(f"scale x{s:g}   k = {e['N']}   chi2/dof = {e['chi2']:.2f}", fontsize=9, loc="left", pad=3)
    ov = e["Y"] >= pr.OVERFLOW
    res = e["Y"] - (e["ped"] + e["comp"].sum(0))
    axr.plot(ed.TT[~ov], res[~ov], "k.-", lw=0.6, ms=3)
    axr.axhline(0, color="r", lw=0.8)
    axr.set_xlim(0, ed.DT * ed.NS)
    axr.tick_params(labelsize=7)
    ax.tick_params(labelsize=7, labelbottom=False)
    axr.set_ylabel("data - fit", fontsize=7)
    axr.set_xlabel("t (ns)", fontsize=8)


def summary(res):
    tab = {}
    for key, lab, want in GROUPS:
        g = [x for x in res if x["cat"] == key]
        ks = {s: np.array([x["ex"][s]["N"] for x in g]) for s in SCALES}
        chi = {s: float(np.median([x["ex"][s]["chi2"] for x in g])) if g else np.nan for s in SCALES}
        if lab in "AC":
            ok = {s: int((ks[s] == 1).sum()) for s in SCALES}     # dropped to one pulse
        else:
            ok = {s: int((ks[s] >= 2).sum()) for s in SCALES}
        tab[key] = dict(n=len(g), ks=ks, chi=chi, ok=ok, want=want, lab=lab)
    return tab


def title_page(pdf, tab, nmiss, nmatch, ntot, elapsed):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.03, 0.95, "Template-shape error inflation: refits of the 50 mono-edge waveforms", fontsize=20, fontweight="bold")
    fig.text(0.03, 0.905,
             "Parameter: sigma_s = bank shape_sigma (per module, floor SHAPE_SIGMA_MIN = 0.005); V_j = sig^2 + (SCALE * sum_k a_k sigma_s(t-t0_k) sqrt(aband/max(a_k,20)))^2.\n"
             "Knob pileup_rates.SHAPE_SCALE (env PILEUP_SHAPE_SCALE), default 1.0, here set to 1, 2, 3 at run time (mono-edge rule off).  The candidate-peak test uses the same scale.\n"
             "Pass criterion: A (H2 small) and C (Ta small) drop to one pulse; B (H2 large) and D (Ta large) keep two (k >= 2).", fontsize=10.5, va="top", family="monospace")
    hdr = ["group", "n", "k=1 @x1", "k=1 @x2", "k=1 @x3", "k>=2 @x1", "k>=2 @x2", "k>=2 @x3", "med chi2/dof x1", "x2", "x3"]
    rows = []
    for key, lab, want in GROUPS:
        t = tab[key]
        n1 = [int((t["ks"][s] == 1).sum()) for s in SCALES]
        n2 = [int((t["ks"][s] >= 2).sum()) for s in SCALES]
        rows.append([f"{lab} {key}", t["n"]] + n1 + n2 + [f"{t['chi'][s]:.2f}" for s in SCALES])
    ax = fig.add_axes([0.03, 0.52, 0.94, 0.28])
    ax.axis("off")
    tb = ax.table(cellText=rows, colLabels=hdr, loc="upper center", cellLoc="center")
    tb.auto_set_font_size(False)
    tb.set_fontsize(11)
    tb.scale(1, 1.8)
    lines = []
    for key, lab, want in GROUPS:
        t = tab[key]
        verdict = "PASS" if (t["n"] and all(t["ok"][s] == t["n"] for s in SCALES[1:])) else "no"
        lines.append(f"  {lab} ({key}): wanted k {'= 1' if lab in 'AC' else '>= 2'};  matching n at x1/x2/x3 = {t['ok'][1.0]}/{t['ok'][2.0]}/{t['ok'][3.0]} of {t['n']}   -> "
                     f"{'all match at x2 and x3' if verdict == 'PASS' else 'not all match at x2/x3'}")
    fig.text(0.03, 0.46, "Measured verdict per group (fraction meeting the criterion):\n" + "\n".join(lines) +
             f"\n\nSanity: x1 refit reproduces the stored rule-off scan (k_off, chi2_off within 0.2% / 1e-3): {nmatch} of {ntot}.   Waveforms not refitted: {nmiss}.   Run time {elapsed:.0f} s.",
             fontsize=11, va="top", family="monospace", linespacing=1.4)
    pdf.savefig(fig)
    plt.close(fig)


def pages(pdf, res):
    for key, lab, want in GROUPS:
        items = [x for x in res if x["cat"] == key]
        for p0 in range(0, len(items), PER):
            chunk = items[p0:p0 + PER]
            fig = plt.figure(figsize=(16, 12.5))
            outer = fig.add_gridspec(PER, 1, hspace=0.32, left=0.04, right=0.99, top=0.86, bottom=0.04)
            fig.text(0.04, 0.97, f"Group {lab}: {GTITLE[key]}", fontsize=15, fontweight="bold")
            fig.text(0.04, 0.948, f"waveforms {p0 + 1}-{p0 + len(chunk)} of {len(items)}.   Columns: shape-error scale x1, x2, x3 (mono-edge off); magenta dotted = fitted onsets.", fontsize=9)
            for i, it in enumerate(chunk):
                gs = outer[i].subgridspec(2, 3, height_ratios=[3, 1], hspace=0.08, wspace=0.12)
                axs = [(fig.add_subplot(gs[0, j]), fig.add_subplot(gs[1, j])) for j in range(3)]
                for (a, r), s in zip(axs, SCALES):
                    panel(a, r, it["ex"][s], s)
                lim = (min(a.get_ylim()[0] for a, _ in axs), max(a.get_ylim()[1] for a, _ in axs))
                rl = (min(r.get_ylim()[0] for _, r in axs), max(r.get_ylim()[1] for _, r in axs))
                for a, r in axs:
                    a.set_ylim(lim)
                    r.set_ylim(rl)
                axs[0][0].text(0.0, 1.17, f"run {it['run']}  ev {it['ev']}  {it['name']}  {it['rt']}   removed pulse (mono-edge study): dt = {it['dt']:.1f} ns, r = {it['r']:.3f}",
                               transform=axs[0][0].transAxes, fontsize=10.5, fontweight="bold")
            pdf.savefig(fig)
            plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--sample-dir", default=str(pmo.PYDIR / "output" / "pileup_rates" / "monoedge_sample"))
    ap.add_argument("--evio-dir", default=pmo.EVIO_DIR)
    ap.add_argument("--out", default=str(HERE))
    a = ap.parse_args()
    t0 = time.time()
    rows = load_list(HERE / "monoedge_fit_comparison_20261116.tsv")
    items = locate(a.sample_dir, a.evio_dir, rows)
    store, ctxs = fit_scales(a.sample_dir, items)
    res, nmiss = build(items, store, ctxs)
    nmatch = sum(x["ex"][1.0]["N"] == x["k_off"] and np.isclose(x["ex"][1.0]["chi2"], x["chi2_off"], rtol=2e-3, atol=1e-3) for x in res)
    tab = summary(res)
    stem = Path(a.out) / "T0_shape_err_sensitivity_20261116"
    with PdfPages(str(stem) + ".pdf") as pdf:
        title_page(pdf, tab, nmiss, nmatch, len(res), time.time() - t0)
        pages(pdf, res)
    with open(str(stem) + ".tsv", "w") as f:
        f.write("run\tevent\tmodule\tgroup\tk_x1\tk_x2\tk_x3\tchi2_x1\tchi2_x2\tchi2_x3\n")
        for x in res:
            e = x["ex"]
            f.write("\t".join([str(x["run"]), str(x["ev"]), x["name"], x["cat"]] + [str(e[s]["N"]) for s in SCALES] + [f"{e[s]['chi2']:.4g}" for s in SCALES]) + "\n")
    for key, lab, want in GROUPS:
        t = tab[key]
        print(f"[{lab} {key}] n={t['n']} k==1: {[int((t['ks'][s] == 1).sum()) for s in SCALES]} k>=2: {[int((t['ks'][s] >= 2).sum()) for s in SCALES]} med chi2 {[round(t['chi'][s], 3) for s in SCALES]}")
    print(f"[done] x1 matches k_off/chi2_off: {nmatch} of {len(res)}; missing {nmiss}; {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
