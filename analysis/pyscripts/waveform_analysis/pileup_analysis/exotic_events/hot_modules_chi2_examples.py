#!/usr/bin/env python3
"""Hot-module (W121, W493) chi2/dof examples: 20 waveforms per module spanning the module's chi2/dof distribution.

Selection (per module, pooled over all runs of the full scans; every kept channel record counts):
  16 quantile picks at the (i+0.5)/16 quantiles of the module's chi2/dof distribution (3.1 % ... 96.9 %), plus 4 labelled tail picks:
  Q99.9, 'hi' (median of the chi2/dof > 1000 records, none if none exist -> Q99.5), 'max' (largest chi2/dof), 'sat' (largest chi2/dof among
  saturated records, nclip >= SAT_NCLIP; if none exist, the 99.95 % quantile).
  For each target value the candidates are the records within +-0.25 % of ranks (or the nearest ones); the pick prefers the run type used
  least so far, then a split already needed for another pick (fewer noise references to build), then a seeded random draw.
Fit: pileup_rates MF default options; PILEUP_NOISE_FLAG=1 only adds the noisy flag (the fit is unchanged), the noisy flag is shown in the panel title.
"""
from __future__ import annotations
import os, sys, glob, time, argparse, textwrap
from pathlib import Path
os.environ["PILEUP_NOISE_FLAG"] = "1"
HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import event_display as ed
import pileup_multiplicity as pm

pr = ed.pr
PYDIR = HERE.parent.parent.parent
MODS = ("W121", "W493")
NQ, SEED = 16, 20261116
COL = ed.PULSE_COLORS

_orig_fit = pr.fit
_last = {}


def _fit(*a, **k):
    f = _orig_fit(*a, **k)
    _last["noisy"] = f.get("noisy")
    return f


pr.fit = _fit
_orig_fe = ed.fit_event


def fit_event(ctx, fe, b_ref, sig_ref):
    _last["noisy"] = None
    r = _orig_fe(ctx, fe, b_ref, sig_ref)
    nz = np.zeros(len(r["mi"]), bool)
    if _last["noisy"] is not None and r["fitted"].any():
        nz[np.nonzero(r["fitted"])[0]] = _last["noisy"]
    r["noisy"] = nz
    return r


ed.fit_event = fit_event


def load(scan_dir, evio_dir):
    """Per run: records of the two modules, event numbers, file index, k histograms of everything."""
    recs = {m: [] for m in MODS}
    kh_all, kh_pwo, kh_lg, kh_mod = np.zeros(8), np.zeros(8), np.zeros(8), {m: np.zeros(8) for m in MODS}
    for fn in sorted(Path(scan_dir).glob("run_*.npz")):
        z = np.load(fn)
        run, rt = int(z["run"]), str(z["run_type"])
        names = list(z["mod_name"])
        files = sorted(glob.glob(str(Path(evio_dir) / f"prad_{run:06d}.evio.*")))
        ctx = ed.Ctx(run)
        firsts = np.array([f if f is not None else 2**62 for f in (ed.first_event_number(ctx, f) for f in files)])
        mi, k = z["ch_midx"], z["ch_k"].astype(int)
        pw = z["mod_pwo4"][mi]
        for a, s in ((kh_all, np.ones(mi.size, bool)), (kh_pwo, pw), (kh_lg, ~pw)):
            a += np.bincount(k[s], minlength=8)[:8]
        for m in MODS:
            i = names.index(m)
            s = np.nonzero(mi == i)[0]
            kh_mod[m] += np.bincount(k[s], minlength=8)[:8]
            evn = z["ev_number"][z["ch_ev"][s]]
            fi = np.searchsorted(firsts, evn, side="right") - 1
            recs[m].append(dict(run=np.full(s.size, run), rt=np.array([rt] * s.size), ev=evn, file=np.array([files[j] for j in fi]) if s.size else np.array([]),
                                chi2=z["ch_chi2"][s].astype(float), k=k[s], nclip=z["ch_nclip"][s].astype(int), midx=np.full(s.size, i), pwo4=bool(z["mod_pwo4"][i])))
        print(f"[load] run {run} {rt}", flush=True)
    out = {}
    for m in MODS:
        out[m] = {q: np.concatenate([r[q] for r in recs[m]]) for q in ("run", "rt", "ev", "file", "chi2", "k", "nclip", "midx")}
        out[m]["pwo4"] = recs[m][0]["pwo4"]
    return out, dict(all=kh_all, pwo=kh_pwo, lg=kh_lg, mod=kh_mod)


def select(d, rng, used_files):
    chi = d["chi2"]
    ok = np.nonzero(np.isfinite(chi))[0]
    order = ok[np.argsort(chi[ok], kind="stable")]
    n = order.size
    targets = [(f"Q{100 * (i + 0.5) / NQ:.1f}", (i + 0.5) / NQ) for i in range(NQ)]
    hi = np.nonzero(chi[order] > pr.CHI2_MAX)[0]
    sat = np.nonzero(d["nclip"][order] >= pr.SAT_NCLIP)[0]
    tails = [("Q99.9", ("q", 0.999))]
    tails.append(("chi2>1000 (median)", ("r", int(hi[hi.size // 2]))) if hi.size else ("Q99.5 (no chi2>1000)", ("q", 0.995)))
    tails.append(("max chi2", ("r", n - 1)))
    if sat.size:
        tails.append(("saturated (max chi2)", ("r", int(sat[np.argmax(chi[order][sat])]))))
    else:
        tails.append(("Q99.95 (no saturated)", ("q", 0.9995)))
    slots = [(lab, ("q", q)) for lab, q in targets] + tails
    rt_used, picks, taken = {}, [], set()
    for lab, (kind, v) in slots:
        r0 = int(round(v * (n - 1))) if kind == "q" else v
        w = max(int(0.0025 * n), 3) if kind == "q" else 0
        cand = [order[j] for j in range(max(r0 - w, 0), min(r0 + w + 1, n))]
        if kind == "r":
            cand = [order[r0]]
        cand = [c for c in cand if c not in taken]
        if not cand:
            continue
        cand = np.array(cand)
        if kind == "q" and lab.startswith("Q99.9"):
            pass
        sc = np.array([(rt_used.get(d["rt"][c], 0), d["file"][c] not in used_files, rng.random()) for c in cand], dtype=object)
        best = cand[min(range(len(cand)), key=lambda i: sc[i].tolist() if False else (sc[i][0], sc[i][1], sc[i][2]))]
        taken.add(best)
        rt_used[d["rt"][best]] = rt_used.get(d["rt"][best], 0) + 1
        used_files.add(d["file"][best])
        picks.append(dict(label=lab, i=int(best), pct=100.0 * (np.searchsorted(chi[ok][np.argsort(chi[ok])], chi[best]) / n)))
    return picks


def refit_all(picksets, data, evio_dir):
    need, wanted, ctxs, refs = {}, {}, {}, {}
    for m in MODS:
        for p in picksets[m]:
            i = p["i"]
            f, ev, run = str(data[m]["file"][i]), int(data[m]["ev"][i]), int(data[m]["run"][i])
            need.setdefault(f, set()).add(ev)
            wanted.setdefault(f, {}).setdefault(ev, []).append(dict(run=run))
            ctxs.setdefault(run, ed.Ctx(run))
    for f in sorted(need):
        run = next(iter(wanted[f].values()))[0]["run"]
        t0 = time.time()
        refs[f] = ed.noise_reference(ctxs[run], f)[1:]
        print(f"[ref] {os.path.basename(f)} {time.time() - t0:.0f}s", flush=True)
    store = pm.decode_wanted(ctxs, refs, wanted, need)
    return ctxs, store


def draw(ax, e):
    y = e["Y"]
    ov = y >= 4095
    ax.plot(ed.TT[~ov], y[~ov], "k.", ms=4)
    if ov.any():
        ax.plot(ed.TT[ov], y[ov], "x", color="tab:red", ms=6, mew=1.5)
    ax.axhline(e["ped"], color="0.55", lw=1.0)
    if e["N"]:
        ax.plot(ed.TT, e["ped"] + e["comp"].sum(0), "r-", lw=1.4)
    for j in range(e["N"]):
        ax.plot(ed.TT, e["ped"] + e["comp"][j], "--", color=COL[j], lw=1.2)
    pulses = "; ".join(f"{a:.0f}@{t:.0f}" for a, t in zip(e["a"], e["t"]))
    ttl = (f"[{e['label']}] run {e['run']} ev {e['ev']} {e['name']} {e['rt']}\nchi2/dof {e['chi2']:.4g}   N = {e['N']}   noisy: {'yes' if e['noisy'] else 'no'}"
           f"   clipped {e['nclip']}\n(A ADC @ t0 ns) " + textwrap.fill(pulses, 48))
    ax.set_title(ttl, fontsize=7.2, loc="left", pad=3)
    ax.tick_params(labelsize=7)
    ax.set_xlim(0, ed.DT * ed.NS)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-dir", default=str(PYDIR / "output" / "pileup_rates" / "scan"))
    ap.add_argument("--evio-dir", default="/Users/L00338853/work/PRad/data/evio")
    ap.add_argument("--stem", default=str(HERE / "hot_modules_W121_W493_chi2_examples"))
    a = ap.parse_args()
    t00 = time.time()
    data, kh = load(a.scan_dir, a.evio_dir)
    rng = np.random.default_rng(SEED)
    used = set()
    picks = {m: select(data[m], rng, used) for m in MODS}
    ctxs, store = refit_all(picks, data, a.evio_dir)
    ex = {m: [] for m in MODS}
    mism = 0
    for m in MODS:
        d = data[m]
        for p in picks[m]:
            i = p["i"]
            f, ev, run, mi = str(d["file"][i]), int(d["ev"][i]), int(d["run"][i]), int(d["midx"][i])
            r = store.get((f, ev))
            if r is None or mi not in r["row"]:
                print(f"[warn] {m} {p['label']}: not refit", flush=True)
                continue
            j = r["row"][mi]
            N = int(r["k"][j])
            e = dict(label=p["label"], run=run, rt=str(d["rt"][i]), ev=ev, name=m, N=N, Y=r["Y"][j].astype(float), ped=float(r["ped"][j]), a=r["a"][j][:N].copy(), t=r["t"][j][:N].copy(),
                     comp=r["comp"][j][:N].copy(), chi2=float(r["chi2"][j]), nclip=int(r["nclip"][j]), noisy=bool(r["noisy"][j]), scan_chi2=float(d["chi2"][i]), scan_k=int(d["k"][i]), pct=p["pct"])
            if not np.isclose(e["chi2"], e["scan_chi2"], rtol=1e-3) or N != e["scan_k"]:
                mism += 1
                print(f"[warn] refit differs: {m} {p['label']} run {run} ev {ev}: scan chi2 {e['scan_chi2']:.4g} k {e['scan_k']} refit {e['chi2']:.4g} N {N}", flush=True)
            ex[m].append(e)
    print(f"[refit] mismatches {mism}", flush=True)
    with open(a.stem + ".tsv", "w") as fo:
        fo.write("module\tlabel\trun\tevent\trun_type\tchi2_dof_refit\tchi2_dof_scan\tN_refit\tN_scan\tpercentile_in_module\tn_clipped\tnoisy\tpulses_A_t0\n")
        for m in MODS:
            for e in ex[m]:
                fo.write(f"{m}\t{e['label']}\t{e['run']}\t{e['ev']}\t{e['rt']}\t{e['chi2']:.5g}\t{e['scan_chi2']:.5g}\t{e['N']}\t{e['scan_k']}\t{e['pct']:.2f}\t{e['nclip']}\t{int(e['noisy'])}\t"
                         + ";".join(f"{x:.1f}@{y:.1f}" for x, y in zip(e["a"], e["t"])) + "\n")
    with PdfPages(a.stem + ".pdf") as pdf:
        for m in MODS:
            d = data[m]
            chi = d["chi2"][np.isfinite(d["chi2"])]
            fig, (a1, a2) = plt.subplots(1, 2, figsize=(17, 7.5), gridspec_kw=dict(width_ratios=[1.5, 1]))
            bins = np.logspace(-1, 5, 61)
            a1.hist(np.clip(chi, bins[0], bins[-1]), bins=bins, color="0.7")
            a1.set_xscale("log"), a1.set_yscale("log")
            for n_, e in enumerate(ex[m]):
                a1.axvline(np.clip(e["scan_chi2"], bins[0], bins[-1]), color="tab:red", lw=0.8, alpha=0.8)
                a1.text(np.clip(e["scan_chi2"], bins[0], bins[-1]), a1.get_ylim()[1] * (0.6 if n_ % 2 else 0.9), str(n_ + 1), color="tab:red", fontsize=8, ha="center")
            a1.axvline(pr.CHI2_MAX, color="tab:blue", ls=":", lw=1)
            a1.set_xlabel("chi2/dof (all runs pooled; values outside 0.1-1e5 in the end bins)"), a1.set_ylabel("records")
            a1.set_title(f"{m} ({'PbWO4' if d['pwo4'] else 'lead glass'}): chi2/dof of {chi.size} fitted records, median {np.median(chi):.2f}, Q99 {np.percentile(chi, 99):.1f}, "
                         f"{(chi > pr.CHI2_MAX).sum()} > 1000, {(d['nclip'] >= pr.SAT_NCLIP).sum()} saturated\nred lines: the 20 picks (numbers = panel order)", fontsize=9)
            Ns = np.arange(1, 7)
            for lab, h, c in ((m, kh["mod"][m], "tab:red"), ("HyCal average", kh["all"], "k"), ("PbWO4" if d["pwo4"] else "lead glass", kh["pwo"] if d["pwo4"] else kh["lg"], "tab:green")):
                a2.plot(Ns, h[1:7] / h[1:7].sum(), "o-", color=c, label=f"{lab} (N={int(h[1:7].sum())})")
            a2.set_yscale("log"), a2.set_xlabel("N fitted pulses"), a2.set_ylabel("fraction of records (N >= 1)")
            a2.set_title("pulse multiplicity, all runs pooled, all kept channels", fontsize=9)
            a2.legend(fontsize=8)
            fig.tight_layout()
            pdf.savefig(fig)
            plt.close(fig)
            fig, axs = plt.subplots(4, 5, figsize=(22, 17))
            for ax in axs.flat:
                ax.axis("off")
            for n_, (ax, e) in enumerate(zip(axs.flat, ex[m])):
                ax.axis("on")
                draw(ax, e)
                ax.text(0.99, 0.95, str(n_ + 1), transform=ax.transAxes, color="tab:red", fontsize=11, ha="right", va="top")
            fig.suptitle(f"{m}: {len(ex[m])} waveforms ordered by chi2/dof (quantile picks then tail picks)", fontsize=18, fontweight="bold")
            fig.text(0.01, 0.003, "black dots: raw samples (red x: overflow >= 4095); grey: fitted pedestal; red: total model; dashed: pulse components in time order; x: t (ns), y: ADC. "
                     "Quantile labels Qxx = percentile of the module's pooled chi2/dof.", fontsize=11)
            fig.tight_layout(rect=(0, 0.015, 1, 0.96), h_pad=1.2)
            pdf.savefig(fig)
            plt.close(fig)
    print(f"[done] {time.time() - t00:.0f}s", flush=True)


if __name__ == "__main__":
    main()
