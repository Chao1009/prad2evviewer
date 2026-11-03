#!/usr/bin/env python3
"""
waveform_good_bad.py - good/bad matched-filter fits of sub-50-ADC pulses, per module.

Companion of waveform_examples.py (same event/waveform selection, same MF
fit and chi2/dof; see that file's docstring).  Two steps:

  collect EVIO OUT.pkl --mf-templates NPZ
      Reads one split, MF-fits EVERY qualifying waveform (RawSum bit, 1 HyCal
      cluster, < 20 fired channels, clean single peak, 20 <= h < 50 ADC) with
      the given spline bank, and stores per module: counts of good
      (chi2/dof < 2), bad (> 3), in-between, no-MF-pulse, plus up to
      --keep good and --keep bad examples (with the fitted curve).
  rates EVIO OUT.pkl --mf-templates NPZ
      Same selection and fit as collect with no examples (adds per-module
      amplitude-bin x class counts and a coarse chi2 histogram); skips an
      existing OUT.pkl.  Aggregated by waveform_rates_aggregate.py.
  render OUT_PDF PKL [PKL ...]
      Pools the pickles in the given order (first found wins) and writes ONE
      PDF, a page per module (name order), 2x2: top row good, bottom row bad,
      empty slots labelled "none found"; plus a summary txt next to it.
"""

from __future__ import annotations

import argparse
import pickle
import sys
from pathlib import Path

import numpy as np

GOOD_LT, BAD_GT = 2.0, 3.0
# rates mode: coarse chi2/dof histogram (last bin = overflow) and amplitude bins
CHI2_EDGES = np.array([0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 4.0, 6.0, 10.0, 1e9])
AMP_EDGES = np.array([20.0, 30.0, 40.0, 50.0])


def _module_sort_key(name: str):
    prefix = name.rstrip("0123456789")
    num = name[len(prefix):]
    return ({"W": 0, "G": 1}.get(prefix, 2), prefix, int(num) if num else 0)


# ------------------------------------------------------------------ collect --
def collect(args) -> int:
    import _common as C
    import mf_pulse
    from fit_pulse_template import _resolve_trigger_bits
    from waveform_examples import _bits_file

    accept_mask = _resolve_trigger_bits("SSP0", _bits_file())
    p = C.setup_pipeline(evio_path=args.evio_path)
    clk_ns = float(p.cfg.wave_cfg.clk_ns)
    run = C.extract_run_number(args.evio_path)
    split = Path(args.evio_path).name.rsplit(".", 1)[-1]

    hy_mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    wcfg = p.wave_ana.cfg
    bank, binfo = mf_pulse.build_bank(
        [m_.name for m_ in hy_mods], [bool(m_.is_pwo4()) for m_ in hy_mods],
        args.mf_templates,
        dict(smooth_order=int(wcfg.smooth_order),
             int_tail_ratio=float(wcfg.int_tail_ratio),
             tail_break_n=int(wcfg.tail_break_n)))
    rec = mf_pulse.MFReconstructor(
        bank, nsigma=5.0, min_amp=5.0, kmax=4, r_model=0.3,
        ped_nsamples=int(wcfg.ped_nsamples))
    overflow = int(wcfg.overflow)
    fallback = set(binfo["fallback"])
    multi_pulse = p.hc_clusterer.get_config().seed_time_window > 0.0

    # waveform selection = waveform_examples.py defaults
    H_MIN, H_RMS, PRE, POST, H_LT, N_FIRED = 20.0, 10.0, 8, 40, 50.0, 20

    mods: dict[str, dict] = {}
    tot = dict(good=0, bad=0, mid=0, nofit=0)
    buf: list = []

    def entry(name):
        return mods.setdefault(name, dict(
            good=0, bad=0, mid=0, nofit=0, ex_good=[], ex_bad=[],
            ab=np.zeros((3, 4), np.int64),
            hist=np.zeros(CHI2_EDGES.size - 1, np.int64)))

    def flush():
        if not buf:
            return
        Y = np.stack([b["samples"] for b in buf])
        mid = np.array([b["midx"] for b in buf], dtype=np.int64)
        res = rec.run(Y, mid, overflow)
        starts = np.searchsorted(res.p_ch, np.arange(len(buf) + 1))
        for i, b in enumerate(buf):
            e = entry(b["name"])
            sl = slice(starts[i], starts[i + 1])
            ia = min(int((b["height"] - AMP_EDGES[0]) // 10.0), 2)
            if starts[i + 1] == starts[i]:
                e["nofit"] += 1
                e["ab"][ia, 3] += 1
                tot["nofit"] += 1
                continue
            chi2 = float(res.p_chi2[sl][0])
            kind = "good" if chi2 < GOOD_LT else ("bad" if chi2 > BAD_GT else "mid")
            e[kind] += 1
            e["ab"][ia, ("good", "bad", "mid").index(kind)] += 1
            e["hist"][min(int(np.searchsorted(CHI2_EDGES, chi2, "right")) - 1,
                          e["hist"].size - 1)] += 1
            tot[kind] += 1
            if kind == "mid" or len(e["ex_" + kind]) >= args.keep:
                continue
            ped, amp, t0 = float(res.ped[i]), res.p_amp[sl].copy(), res.p_t0[sl].copy()
            tt = np.arange(0.0, Y.shape[1] * clk_ns, 0.5)
            fit = np.full_like(tt, ped)
            for a_, t0_ in zip(amp, t0):
                fit += a_ * bank.shape(np.array([b["midx"]])[:, None],
                                       (tt - t0_)[None, :])[0]
            e["ex_" + kind].append(dict(
                run=run, split=split, evn=b["evn"], samples=b["samples"],
                ped_wa=b["ped"], height=b["height"], mf_ped=ped, amp=amp, t0=t0,
                chi2=chi2, npk=int(starts[i + 1] - starts[i]),
                fallback=b["name"] in fallback, tt=tt, fit=fit.astype(np.float32)))
        buf.clear()

    stats = C.LoopStats()
    cnt = dict(trigger=0, fired=0, cluster=0)
    for fadc_evt, _ in C.iter_physics_events(
            p, stats, with_ssp=False, max_events=args.max_events,
            accept=lambda b: bool(b & accept_mask), progress=None):
        cnt["trigger"] += 1
        chans, n_fired = [], 0
        for tag, s, c, cd in C.iter_fadc_channels(fadc_evt):
            samples = np.asarray(cd.samples, dtype=np.uint16)
            ped, rms, peaks = p.wave_ana.analyze(samples)
            chans.append((tag, s, c, samples, ped, rms, peaks))
            if peaks:
                pk = peaks[0]
                if pk.quality == 0 and not pk.overflow and pk.height >= H_MIN:
                    n_fired += 1
        if n_fired >= N_FIRED:
            continue
        cnt["fired"] += 1
        hy = []
        p.hc_clusterer.clear()
        for tag, s, c, samples, ped, rms, peaks in chans:
            crate = p.crate_map.get(tag)
            if crate is None:
                continue
            mod = p.hycal.module_by_daq(crate, s, c)
            if mod is None or not mod.is_hycal():
                continue
            C.add_hycal_peaks(p.hc_clusterer, mod, peaks, multi_pulse)
            hy.append((mod, samples, ped, rms, peaks))
        p.hc_clusterer.form_clusters()
        if len(p.hc_clusterer.reconstruct_hits()) != 1:
            continue
        cnt["cluster"] += 1
        evn = int(fadc_evt.info.event_number)
        for mod, samples, ped, rms, peaks in hy:
            if len(peaks) != 1:
                continue
            pk = peaks[0]
            if pk.quality != 0 or pk.overflow:
                continue
            if pk.pos - PRE < 0 or pk.pos + POST + 1 > samples.shape[0]:
                continue
            h = float(pk.height)
            if h < H_MIN or h < H_RMS * rms or h >= H_LT:
                continue
            buf.append(dict(name=mod.name, midx=int(mod.index), evn=evn,
                            samples=samples.copy(), ped=float(ped), height=h))
        if len(buf) >= 20000:
            flush()
    flush()

    out = dict(run=run, split=split, n_trigger=cnt["trigger"], template=args.mf_templates,
               sha1=binfo["sha1"], run_type=binfo["run_type"],
               n_fallback=len(fallback), cnt=cnt, tot=tot, n_phys=stats.n_phys,
               mods=mods)
    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    with open(args.out, "wb") as f:
        pickle.dump(out, f)
    print(f"[done] {args.out} run {run} split {split} phys={stats.n_phys} "
          f"{cnt} good/bad/mid/nofit={tot} modules={len(mods)} "
          f"{stats.elapsed:.0f}s", flush=True)
    return 0


# ------------------------------------------------------------------- render --
def render(args) -> int:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    data = []
    for pk in args.pkls:
        with open(pk, "rb") as f:
            data.append(pickle.load(f))
    names = sorted({n for d in data for n in d["mods"]}, key=_module_sort_key)
    sel = {}                       # name -> (good list, bad list)
    counts = {}                    # name -> [good, bad, mid, nofit]
    for n in names:
        g, b = [], []
        c = [0, 0, 0, 0]
        for d in data:
            e = d["mods"].get(n)
            if e is None:
                continue
            for k, key in enumerate(("good", "bad", "mid", "nofit")):
                c[k] += e[key]
            g += e["ex_good"]
            b += e["ex_bad"]
        sel[n] = (g[:2], b[:2])
        counts[n] = c

    cut_txt = ("trigger & SSP0 | 1 HyCal cluster | fired channels < 20 | "
               "1 clean peak, 20 <= h < 50 ADC, h >= 10 x ped rms\n"
               f"MF spline fit (mf_pulse, nsigma 5 / min_amp 5 / kmax 4 / r_model 0.3), "
               f"each run with its own R5 template bank; good chi2/dof < {GOOD_LT:g}, "
               f"bad > {BAD_GT:g}")
    out_pdf = Path(args.out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    with PdfPages(out_pdf) as pdf:
        for n in names:
            fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
            g, b = sel[n]
            for row, (lab, lst, col) in enumerate(
                    (("good", g, "C2"), ("bad", b, "C3"))):
                for j in range(2):
                    ax = axes[row, j]
                    if j >= len(lst):
                        ax.text(0.5, 0.5, f"none found ({lab})", ha="center",
                                va="center", transform=ax.transAxes, color="0.5")
                        ax.set_xticks([]); ax.set_yticks([])
                        continue
                    ex = lst[j]
                    y = ex["samples"].astype(float)
                    x = np.arange(y.size) * 4.0
                    ax.plot(x, y, "-", color="C0", lw=1.0, marker=".", ms=3)
                    ax.axhline(ex["mf_ped"], color="0.4", ls="--", lw=1.0)
                    ax.plot(ex["tt"], ex["fit"], "-", color=col, lw=1.6)
                    txt = "\n".join(f"A{k + 1} {a_:.1f} ADC  t0 {t_:.1f} ns"
                                    for k, (a_, t_) in enumerate(zip(ex["amp"], ex["t0"])))
                    txt += f"\nped {ex['mf_ped']:.1f}  chi2/dof {ex['chi2']:.2f}"
                    if ex["fallback"]:
                        txt += "\n[global template]"
                    ax.text(0.97, 0.55, txt, transform=ax.transAxes, fontsize=7.5,
                            ha="right", va="top", color=col,
                            bbox=dict(fc="w", ec="0.8", alpha=0.85))
                    ax.set_title(f"{lab.upper()}  run {ex['run']} split {ex['split']}  "
                                 f"event {ex['evn']}", fontsize=9)
                    ax.set_xlabel("time (ns)")
                    ax.set_ylabel("ADC")
                    ax.grid(alpha=0.3)
            c = counts[n]
            fig.suptitle(f"{n}  -  found good {c[0]}  bad {c[1]}  in-between {c[2]}  "
                         f"no-MF-pulse {c[3]}  (all pooled splits)", fontsize=12)
            fig.text(0.5, 0.005, cut_txt, ha="center", va="bottom", fontsize=7.5)
            fig.tight_layout(rect=(0, 0.06, 1, 0.96))
            pdf.savefig(fig)
            plt.close(fig)

    full = sum(1 for n in names if len(sel[n][0]) == 2 and len(sel[n][1]) == 2)
    none = sum(1 for n in names if not sel[n][0] and not sel[n][1])
    part = len(names) - full - none
    out_txt = out_pdf.with_suffix(".txt")
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write("# good/bad MF fits of sub-50-ADC pulses (waveform_good_bad.py)\n")
        f.write("# " + cut_txt.replace("\n", "\n# ") + "\n")
        f.write(f"inputs\t{len(data)} splits: "
                + " ".join(f"{d['run']}_{d['split']}" for d in data) + "\n")
        f.write(f"modules_with_any_qualifying_waveform\t{len(names)}\n")
        f.write(f"pdf_pages\t{len(names)}\n")
        f.write(f"modules_full_2good_2bad\t{full}\n")
        f.write(f"modules_partial\t{part}\n")
        f.write(f"modules_none\t{none}\n")
        T = np.sum([[d["tot"][k] for k in ("good", "bad", "mid", "nofit")]
                    for d in data], axis=0)
        fit_n = max(int(T[:3].sum()), 1)
        f.write(f"overall_good_bad_mid_nofit\t{T[0]}\t{T[1]}\t{T[2]}\t{T[3]}\n")
        f.write(f"overall_fraction_good_bad_mid(of fitted)\t{T[0]/fit_n:.3f}\t"
                f"{T[1]/fit_n:.3f}\t{T[2]/fit_n:.3f}\n")
        f.write("# per run/split: run split template good bad mid nofit "
                "frac_good frac_bad frac_mid\n")
        for d in data:
            t = d["tot"]
            n_ = max(t["good"] + t["bad"] + t["mid"], 1)
            f.write(f"run\t{d['run']}\t{d['split']}\t{Path(d['template']).parent.name}\t"
                    f"{t['good']}\t{t['bad']}\t{t['mid']}\t{t['nofit']}\t"
                    f"{t['good']/n_:.3f}\t{t['bad']/n_:.3f}\t{t['mid']/n_:.3f}\n")
        f.write("module\tgood_found\tbad_found\tmid\tnofit\tgood_shown\tbad_shown\n")
        for n in names:
            c = counts[n]
            f.write(f"{n}\t{c[0]}\t{c[1]}\t{c[2]}\t{c[3]}\t{len(sel[n][0])}\t"
                    f"{len(sel[n][1])}\n")
    print(f"[done] {out_pdf} ({len(names)} pages) full={full} partial={part} "
          f"none={none}", flush=True)
    return 0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    c = sub.add_parser("collect")
    c.add_argument("evio_path")
    c.add_argument("out")
    c.add_argument("--mf-templates", required=True)
    c.add_argument("--keep", type=int, default=2,
                   help="Good and bad examples kept per module per split.")
    c.add_argument("--max-events", type=int, default=0)
    q = sub.add_parser("rates", help="collect with --keep 0 (counts, chi2 hist, "
                       "amplitude bins only; no examples)")
    q.add_argument("evio_path")
    q.add_argument("out")
    q.add_argument("--mf-templates", required=True)
    q.add_argument("--max-events", type=int, default=0)
    r = sub.add_parser("render")
    r.add_argument("out_pdf")
    r.add_argument("pkls", nargs="+")
    args = ap.parse_args()
    if args.cmd == "rates":
        args.keep = 0
        if Path(args.out).exists():
            print(f"[skip] {args.out} exists", flush=True)
            return 0
        return collect(args)
    return collect(args) if args.cmd == "collect" else render(args)


if __name__ == "__main__":
    sys.exit(main())
