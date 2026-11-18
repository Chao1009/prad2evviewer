#!/usr/bin/env python3
"""pileup_noisybase.py -- effect of the noisy-baseline flag (pileup_rates.py PILEUP_NOISY_BASE=1) on the pile-up multiplicity, with 20 flagged examples.

    source setup_env.sh
    bash waveform_analysis/pileup_analysis/run_noisybase_sample.sh        # scans the 2-split sample (same splits as the monoedge sample) with the option on
    python3 waveform_analysis/pileup_analysis/pileup_noisybase.py report [--stem <stem>]

Basis: the full scans do not store the baseline std, so the 2-split sample is used: old = monoedge_sample/old (option off), new = noisybase_sample (option on),
the same splits.  Noisy baseline: s_b = std (ddof = 1) of samples 0-29 > 5 sigma_ref (module, file); the waveform is not fitted and leaves R_N numerator and denominator.
Output: <stem>.pdf (new definition/result pages, the original 20261115 pages reused, sample tables old vs new, the 20 noisy-baseline examples),
<stem>_table.tsv, <stem>_examples.tsv, <stem>_modules.tsv.
"""
from __future__ import annotations

import sys
import argparse
import textwrap
import time
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
import pileup_multiplicity as pm                        # noqa: E402

pr, ed = pm.pr, pm.ed
PYDIR = HERE.parent.parent
SEED = 20261116
POPS = pm.POPS
K = pr.KMAX


def load(old_dir, new_dir):
    """Per population / run type: n_old, n_new [KMAX+1]; flagged counts; per module counts; example pool."""
    R = {p: {rt: dict(old=np.zeros(K + 1, np.int64), new=np.zeros(K + 1, np.int64), dig=0, nb=0, pre=0, pulse=0, runs=[]) for rt in pr.RUN_TYPES} for p, _ in POPS}
    mod = {}
    pool = []
    for fn in sorted(Path(new_dir).glob("run_*.npz")):
        zn = np.load(fn)
        zo = np.load(Path(old_dir) / fn.name)
        rt, run = str(zn["run_type"]), int(zn["run"])
        names, pw = zn["mod_name"], zn["mod_pwo4"]
        assert zn["n_pass"] == zo["n_pass"], fn
        nd = zn["digit_count"] - zn["bad_count"]
        for p, sel in POPS:
            m = np.ones(names.size, bool) if sel is None else (pw == sel)
            r = R[p][rt]
            for tag, z in (("old", zo), ("new", zn)):
                k = z["ch_k"]
                r[tag] += np.bincount(k[m[z["ch_midx"]]], minlength=K + 1)[:K + 1]
            r["dig"] += int(nd[m].sum())
            r["nb"] += int(zn["mod_nb_count"][m].sum())
            r["pre"] += int(zn["mod_nb_pre"][m].sum())
            r["pulse"] += int(zn["mod_nb_pulse"][m].sum())
            if p == POPS[0][0]:
                r["runs"].append(run)
        for i in range(names.size):
            d = mod.setdefault(str(names[i]), dict(pw=bool(pw[i]), dig=0, nb=0, pre=0, pulse=0))
            d["dig"] += int(nd[i])
            d["nb"] += int(zn["mod_nb_count"][i])
            d["pre"] += int(zn["mod_nb_pre"][i])
            d["pulse"] += int(zn["mod_nb_pulse"][i])
        if "nbx_samples" in zn.files:
            for j in range(zn["nbx_midx"].size):
                pool.append(dict(Y=zn["nbx_samples"][j], run=run, ev=int(zn["nbx_evnum"][j]), name=str(names[zn["nbx_midx"][j]]), pwo4=bool(pw[zn["nbx_midx"][j]]),
                                 rt=rt, z=float(zn["nbx_z"][j]), pulse=bool(zn["nbx_pulse"][j]), key=float(zn["nbx_key"][j])))
    return R, mod, pool


def rn_rows(R):
    rows = []
    for p, _ in POPS:
        for rt in pr.RUN_TYPES:
            r = R[p][rt]
            for tag in ("old", "new"):
                n = r[tag]
                tot = int(n[1:].sum())
                for N in range(1, K + 1):
                    f = n[N] / tot if tot else np.nan
                    rows.append(dict(pop=p, rt=rt, basis=tag, N=N, count=int(n[N]), denom=tot, rate=f, err=np.sqrt(f * (1 - f) / tot) if tot else np.nan))
    return rows


def pct(a, b):
    return f"{100 * a / b:.3f}" if b else "-"


def table_page(pdf, R, rows, pop):
    fig = plt.figure(figsize=(16, 9))
    ax = fig.add_axes([0.02, 0.22, 0.96, 0.60])
    ax.axis("off")
    hdr = ["run type", "sample", "n_all"] + [f"R_{N} (%)" for N in range(1, 7)] + ["1 - R_1 (%)"]
    cells = []
    for rt in pr.RUN_TYPES:
        for tag in ("old", "new"):
            sub = {r["N"]: r for r in rows if r["pop"] == pop and r["rt"] == rt and r["basis"] == tag}
            tot = sub[1]["denom"]
            lab = "old (no flag)" if tag == "old" else "new (noisy base. excl.)"
            cells.append([rt if tag == "old" else "", lab, f"{tot:,}"] + [pm.fmt_rate(sub[N]["rate"], sub[N]["err"]) for N in range(1, 7)] + [f"{100 * (1 - sub[1]['rate']):.4f}"])
    t = ax.table(cellText=cells, colLabels=hdr, loc="upper center", cellLoc="center", colWidths=[0.12, 0.15, 0.09] + [0.085] * 6 + [0.08])
    t.auto_set_font_size(False)
    t.set_fontsize(9)
    t.scale(1, 1.9)
    for (i, j), c in t.get_celld().items():
        if i == 0:
            c.set_text_props(fontweight="bold")
        elif i % 2 == 0:
            c.set_facecolor("#eef3fa")
    fl = [f"{rt}: {R[pop][rt]['nb']:,} flagged of {R[pop][rt]['dig']:,} non-bad digitised ({pct(R[pop][rt]['nb'], R[pop][rt]['dig'])} %)" for rt in pr.RUN_TYPES]
    fig.suptitle(f"Noisy-baseline exclusion, 2-split sample: {pop}", fontsize=18, fontweight="bold")
    fig.text(0.03, 0.86, "R_N = n_N / n_all +- binomial error (%), same splits for both rows; old = option off, new = waveforms with s_b > 5 sigma_ref not fitted and removed from n_N and n_all.", fontsize=10)
    fig.text(0.03, 0.04, "Flagged waveforms (all non-bad digitised channels):\n  " + "\n  ".join(fl), fontsize=10, family="monospace")
    pdf.savefig(fig)
    plt.close(fig)


def summary_page(pdf, R, mod, ncheck):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.04, 0.94, "Noisy baseline: definition and result (2-split sample)", fontsize=22, fontweight="bold")
    df = ("Definition\n"
          f"  s_b = standard deviation (ddof = 1) of the {pr.PED_NSAMPLES} pedestal samples (0-{pr.PED_NSAMPLES - 1}) of the raw waveform; sigma_ref = the module's noise reference (per module, per file; first 20000 events).\n"
          f"  Noisy baseline: s_b > {pr.NB_NSIG:g} sigma_ref.  Flagged waveforms are NOT fitted (no pulses) and leave R_N numerator and denominator.  Applied after the bad-waveform removal.\n"
          "  pileup_rates.py option PILEUP_NOISY_BASE=1 (default off, output unchanged); MONO_EDGE and NOISE_FLAG off.  Basis: the full scans do not store the baseline std, so the\n"
          "  2-split sample (same splits as the monoedge sample) was re-scanned with the option on; 'old' is the monoedge-sample scan without it.\n"
          f"  Early-pulse caveat: a real pulse inside samples 0-{pr.PED_NSAMPLES - 1} inflates s_b.  'pulse-like' = max - median > {pr.NB_PULSE_EXC:g} sigma_ref and all samples beyond 3 sigma_ref of the\n"
          f"  median within <= {pr.NB_PULSE_SPAN} consecutive samples; counted among flagged waveforms with raw max - baseline >= {pr.PRECUT_ADC:g} ADC (those that would have been fitted).\n")
    fig.text(0.04, 0.89, df, fontsize=10.5, va="top", family="monospace", linespacing=1.3)
    hdr = ["run type", "population", "flagged / non-bad", "fraction (%)", "kept removed", "of old n_all (%)", "dR_1 (pp)", "pulse-like / flagged>=8 ADC"]
    cells = []
    for rt in pr.RUN_TYPES:
        for p, _ in POPS:
            r = R[p][rt]
            no, nn = int(r["old"][1:].sum()), int(r["new"][1:].sum())
            d1 = 100 * (r["new"][1] / nn - r["old"][1] / no) if no and nn else np.nan
            cells.append([rt if p == POPS[0][0] else "", p, f"{r['nb']:,} / {r['dig']:,}", pct(r["nb"], r["dig"]), f"{no - nn:,}", pct(no - nn, no), f"{d1:+.4f}", f"{r['pulse']:,} / {r['pre']:,} ({pct(r['pulse'], r['pre'])} %)"])
    ax = fig.add_axes([0.02, 0.05, 0.96, 0.45])
    ax.axis("off")
    t = ax.table(cellText=cells, colLabels=hdr, loc="upper center", cellLoc="center", colWidths=[0.12, 0.1, 0.17, 0.09, 0.1, 0.11, 0.09, 0.22])
    t.auto_set_font_size(False)
    t.set_fontsize(8.5)
    t.scale(1, 1.25)
    for (i, j), c in t.get_celld().items():
        if i == 0:
            c.set_text_props(fontweight="bold")
        elif ((i - 1) // 3) % 2 == 0:
            c.set_facecolor("#eef3fa")
    pdf.savefig(fig)
    plt.close(fig)
    top = sorted(((v["nb"] / v["dig"], n, v) for n, v in mod.items() if v["dig"] >= 5000), reverse=True)[:10]
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.04, 0.94, "Noisy baseline: top 10 modules and checks", fontsize=22, fontweight="bold")
    lines = ["module  type   flagged / non-bad digitised   fraction (%)   flagged>=8 ADC   pulse-like"]
    for f, n, v in top:
        lines.append(f"{n:6s}  {'PbWO4' if v['pw'] else 'LG   '}  {v['nb']:8,} / {v['dig']:11,}   {100 * f:10.3f}   {v['pre']:12,}   {v['pulse']:9,}")
    lines += ["", "(modules with >= 5000 non-bad digitised waveforms in the sample, ranked by flagged fraction)", "",
              f"check: {ncheck[0]:,} of {ncheck[1]:,} kept channels of the new scan have the same fit (k) as the old scan; the rest of the old kept channels are the removed ones."]
    fig.text(0.04, 0.88, "\n".join(lines), fontsize=11, va="top", family="monospace")
    pdf.savefig(fig)
    plt.close(fig)
    return top


def choose(pool):
    """20 examples, 5 per run type (fixed seed): 2 pulse-like + 3 broadband where available, the rest from any pool of that run type, then from anything."""
    rng = np.random.default_rng(SEED)
    chosen = []
    for rt in pr.RUN_TYPES:
        sub = [e for e in pool if e["rt"] == rt]
        sub.sort(key=lambda e: e["key"])
        pu = [e for e in sub if e["pulse"]]
        br = [e for e in sub if not e["pulse"]]
        pick = [pu[i] for i in rng.permutation(len(pu))[:2]] + [br[i] for i in rng.permutation(len(br))[:3]]
        rest = [e for e in sub if all(e is not q for q in pick)]
        pick += [rest[i] for i in rng.permutation(len(rest))[:5 - len(pick)]]
        chosen += pick
    if len(chosen) < 20:
        rest = [e for e in pool if all(e is not q for q in chosen)]
        chosen += [rest[i] for i in rng.permutation(len(rest))[:20 - len(chosen)]]
    return chosen


def examples_page(pdf, ex, nflag):
    fig, axs = plt.subplots(4, 5, figsize=(22, 17))
    for ax in axs.flat:
        ax.axis("off")
    for ax, e in zip(axs.flat, ex):
        ax.axis("on")
        y = e["Y"]
        ov = y >= 4095
        ax.axvspan(0, pr.PED_NSAMPLES * ed.DT, color="tab:orange", alpha=0.18)
        ax.plot(ed.TT[~ov], y[~ov], "k.-", ms=4, lw=0.5)
        if ov.any():
            ax.plot(ed.TT[ov], y[ov], "x", color="tab:red", ms=6, mew=1.5)
        sb = float(np.std(y[:pr.PED_NSAMPLES].astype(float), ddof=1))
        ax.set_title(f"run {e['run']} ev {e['ev']} {e['name']} ({'PbWO4' if e['pwo4'] else 'LG'}) {e['rt']}\ns_b = {sb:.1f} ADC = {e['z']:.1f} sigma_ref   "
                     f"[{'pulse-like (early pulse)' if e['pulse'] else 'broadband noise'}]", fontsize=8, loc="left", pad=3)
        ax.tick_params(labelsize=7)
        ax.set_xlim(0, ed.DT * ed.NS)
    fig.suptitle("Noisy-baseline waveforms (s_b > 5 sigma_ref): 20 random examples, not fitted", fontsize=18, fontweight="bold")
    fig.text(0.01, 0.003, f"black: raw samples (red x: overflow >= 4095); orange band: the {pr.PED_NSAMPLES}-sample pedestal window used for s_b.  No fit is drawn (none is done).  Seed {SEED}; "
             "5 per run type (2 pulse-like + 3 broadband where available) drawn from flagged waveforms with raw max - baseline >= 8 ADC in the sample.  x: t (ns), y: ADC.", fontsize=11)
    fig.tight_layout(rect=(0, 0.015, 1, 0.96), h_pad=1.2)
    pdf.savefig(fig)
    plt.close(fig)


def consistency(old_dir, new_dir):
    """kept channels of the new scan that match an old kept channel (event, module) with the same k."""
    same = tot = 0
    for fn in sorted(Path(new_dir).glob("run_*.npz"))[:4]:
        zn, zo = np.load(fn), np.load(Path(old_dir) / fn.name)
        kn = zn["ch_ev"].astype(np.int64) * 4096 + zn["ch_midx"]
        ko = zo["ch_ev"].astype(np.int64) * 4096 + zo["ch_midx"]
        o = np.argsort(ko)
        i = np.searchsorted(ko[o], kn)
        i = np.minimum(i, ko.size - 1)
        ok = ko[o][i] == kn
        same += int((ok & (zo["ch_k"][o][i] == zn["ch_k"])).sum())
        tot += kn.size
    return same, tot


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["report"])
    ap.add_argument("--old-dir", default=str(PYDIR / "output" / "pileup_rates" / "monoedge_sample" / "old"))
    ap.add_argument("--new-dir", default=str(PYDIR / "output" / "pileup_rates" / "noisybase_sample"))
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--stem", default="pileup_multiplicity_20261116_noisybase")
    ap.add_argument("--orig-stem", default="pileup_multiplicity_20261115")
    a = ap.parse_args()
    from pypdf import PdfReader, PdfWriter
    out, t0 = Path(a.out), time.time()
    R, mod, pool = load(a.old_dir, a.new_dir)
    ncheck = consistency(a.old_dir, a.new_dir)
    rows = rn_rows(R)
    print(f"[load] {time.time() - t0:.0f}s, pool {len(pool)}, consistency {ncheck}", flush=True)
    with open(out / f"{a.stem}_table.tsv", "w") as f:
        f.write("# 2-split sample; old = noisy-baseline option off, new = s_b > 5 sigma_ref waveforms not fitted and removed from numerator and denominator\n")
        f.write("population\trun_type\tN\told_count\told_denominator\told_rate\tnew_count\tnew_denominator\tnew_rate\tnew_rate_err_binom\tflagged\tnon_bad_digitised\tflagged_ge8ADC\tpulse_like\n")
        for p, _ in POPS:
            for rt in pr.RUN_TYPES:
                r = R[p][rt]
                for N in range(1, K + 1):
                    o = next(x for x in rows if x["pop"] == p and x["rt"] == rt and x["basis"] == "old" and x["N"] == N)
                    n = next(x for x in rows if x["pop"] == p and x["rt"] == rt and x["basis"] == "new" and x["N"] == N)
                    f.write(f"{p}\t{rt}\t{N}\t{o['count']}\t{o['denom']}\t{o['rate']:.6g}\t{n['count']}\t{n['denom']}\t{n['rate']:.6g}\t{n['err']:.3g}\t{r['nb']}\t{r['dig']}\t{r['pre']}\t{r['pulse']}\n")
    with open(out / f"{a.stem}_modules.tsv", "w") as f:
        f.write("module\tis_pbwo4\tnon_bad_digitised\tflagged\tfraction\tflagged_ge8ADC\tpulse_like\n")
        for n, v in sorted(mod.items(), key=lambda kv: -kv[1]["nb"] / max(kv[1]["dig"], 1)):
            f.write(f"{n}\t{int(v['pw'])}\t{v['dig']}\t{v['nb']}\t{v['nb'] / max(v['dig'], 1):.6g}\t{v['pre']}\t{v['pulse']}\n")
    ex = choose(pool)
    with open(out / f"{a.stem}_examples.tsv", "w") as f:
        f.write("panel\trun\tevent\tmodule\tpopulation\trun_type\ts_b_over_sigma_ref\tclass\n")
        for i, e in enumerate(ex):
            f.write(f"{i + 1}\t{e['run']}\t{e['ev']}\t{e['name']}\t{'PbWO4' if e['pwo4'] else 'LG'}\t{e['rt']}\t{e['z']:.2f}\t{'pulse-like' if e['pulse'] else 'broadband'}\n")
    tmp = out / f".{a.stem}_front.pdf"
    with PdfPages(tmp) as pdf:
        summary_page(pdf, R, mod, ncheck)
    tmp2 = out / f".{a.stem}_mid.pdf"
    with PdfPages(tmp2) as pdf:
        for p, _ in POPS:
            table_page(pdf, R, rows, p)
    tmp3 = out / f".{a.stem}_ex.pdf"
    with PdfPages(tmp3) as pdf:
        examples_page(pdf, ex, None)
    orig = PdfReader(out / f"{a.orig_stem}.pdf")
    w = PdfWriter()
    for pg in PdfReader(tmp).pages:
        w.add_page(pg)
    w.add_page(orig.pages[0])
    for pg in PdfReader(tmp2).pages:
        w.add_page(pg)
    for pg in orig.pages[1:]:
        w.add_page(pg)
    for pg in PdfReader(tmp3).pages:
        w.add_page(pg)
    with open(out / f"{a.stem}.pdf", "wb") as f:
        w.write(f)
    for t in (tmp, tmp2, tmp3):
        t.unlink()
    print(f"[done] {time.time() - t0:.0f}s", flush=True)


if __name__ == "__main__":
    main()
