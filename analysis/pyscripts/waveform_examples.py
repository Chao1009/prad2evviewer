#!/usr/bin/env python3
"""
waveform_examples.py - per-module example waveforms for one EVIO split.

Writes one PDF page per HyCal module that has at least one qualifying
waveform; each page is a 2x2 grid of up to --n-examples waveforms (raw
samples vs ADC, pedestal line).  A summary .txt with the cut flow and the
per-module example counts is written next to the PDF.

Event selection (in this order):
  1. Trigger: the trigger word has any bit of --trigger-accept set
     (default SSP0 = bit 8, "Raw Sum", resolved from trigger_bits.json).
  2. Multiplicity: fewer than --n-fired-lt channels "fired".  A channel
     fired when its first WaveAnalyzer peak is clean (quality == 0, no
     overflow) and its height >= --module-count-height-min (default -1 =
     reuse --height-min).  Same definition as fit_pulse_template.py's
     --max-modules-per-event filter (all FADC channels counted).
  3. Exactly --n-clusters HyCal clusters (default 1), built with the
     pipeline-default clustering config and calibration, the same way as
     _common.reconstruct_hycal().

Waveform selection (per HyCal module, fit_pulse_template.py's "good pulse"
cuts at their defaults plus an upper height cut):
  exactly one peak, quality == 0, no overflow, [pos - pre, pos + post]
  inside the buffer, height >= --height-min, height >= --height-rms-mult x
  pedestal rms, and height < --height-lt (default 50 ADC).  "height" is the
  WaveAnalyzer peak height above pedestal (ADC counts).

Matched-filter fit (mf_pulse.MFReconstructor.run, the R5 --pulse-method mf
code, parameters at the benchmark defaults: nsigma 5, min_amp 5, kmax 4,
r_model 0.3): every example waveform is fitted with the per-module spline
template of --mf-templates (a spline_templates_v1 .npz; modules without their
own shape use the bank's average shape, marked "global template" on the panel).
Fitted per pulse: amplitude a (ADC, unit-peak template) and onset t0 (ns); the
pedestal is the MF's own clipped mean of samples 0-29 (not a fit parameter).
chi2/dof uses sigma^2 = ped_sigma^2 + (0.03 x model)^2 over the fit window.

Usage:
  python waveform_examples.py /path/prad_024386.evio.00010 out.pdf \
      --mf-templates output/templates_spline_best50_20261026/3p49GeV_H2/templates_3p49GeV_H2.npz
"""

from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path

import numpy as np

import _common as C
import mf_pulse
from fit_pulse_template import _resolve_trigger_bits


def _bits_file() -> Path:
    if os.environ.get("PRAD2_DATABASE_DIR"):
        return Path(os.environ["PRAD2_DATABASE_DIR"]) / "trigger_bits.json"
    return Path(__file__).resolve().parents[2] / "database" / "trigger_bits.json"


def _module_sort_key(name: str):
    prefix = name.rstrip("0123456789")
    num = name[len(prefix):]
    return ({"W": 0, "G": 1}.get(prefix, 2), prefix, int(num) if num else 0)


def main() -> int:
    ap = argparse.ArgumentParser(
        description="Per-module example waveforms (one PDF page per module).",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    ap.add_argument("evio_path", help="One EVIO split file.")
    ap.add_argument("out_pdf", help="Output PDF path (summary .txt alongside).")
    ap.add_argument("--mf-templates", required=True, metavar="NPZ",
                    help="Spline template bank (spline_templates_v1 .npz).")
    ap.add_argument("--mf-nsigma", type=float, default=5.0)
    ap.add_argument("--mf-min-amp", type=float, default=5.0)
    ap.add_argument("--mf-kmax", type=int, default=4)
    ap.add_argument("--mf-r-model", type=float, default=0.3)
    ap.add_argument("--max-events", type=int, default=0,
                    help="Stop after N physics events (0 = whole split).")
    ap.add_argument("--n-examples", type=int, default=4,
                    help="Waveforms per module page (2x2 grid).")
    ap.add_argument("--trigger-accept", default="SSP0",
                    help="Trigger names/bits (trigger_bits.json); event "
                         "passes if any is set.")
    ap.add_argument("--n-fired-lt", type=int, default=20,
                    help="Keep events with fewer than this many fired channels.")
    ap.add_argument("--n-clusters", type=int, default=1,
                    help="Required number of HyCal clusters.")
    ap.add_argument("--height-lt", type=float, default=50.0,
                    help="Waveform peak height above pedestal must be < this (ADC).")
    # fit_pulse_template.py defaults
    ap.add_argument("--height-min", type=float, default=20.0,
                    help="Minimum peak height (ADC, pedsub).")
    ap.add_argument("--height-rms-mult", type=float, default=10.0,
                    help="Also require peak height >= this x ped rms.")
    ap.add_argument("--pre-samples", type=int, default=8,
                    help="Peak must have this many samples before it.")
    ap.add_argument("--post-samples", type=int, default=40,
                    help="Peak must have this many samples after it.")
    ap.add_argument("--module-count-height-min", type=float, default=-1.0,
                    help="Height for a channel to count as fired "
                         "(-1 = reuse --height-min).")
    args = ap.parse_args()

    accept_mask = _resolve_trigger_bits(args.trigger_accept, _bits_file())
    fired_hmin = (args.module_count_height_min
                  if args.module_count_height_min >= 0 else args.height_min)

    p = C.setup_pipeline(evio_path=args.evio_path)
    clk_ns = float(p.cfg.wave_cfg.clk_ns)
    # Matched-filter bank/reconstructor, built as in the benchmark (setup block).
    hy_mods = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    wcfg = p.wave_ana.cfg
    bank, binfo = mf_pulse.build_bank(
        [m_.name for m_ in hy_mods], [bool(m_.is_pwo4()) for m_ in hy_mods],
        args.mf_templates,
        dict(smooth_order=int(wcfg.smooth_order),
             int_tail_ratio=float(wcfg.int_tail_ratio),
             tail_break_n=int(wcfg.tail_break_n)))
    rec = mf_pulse.MFReconstructor(
        bank, nsigma=args.mf_nsigma, min_amp=args.mf_min_amp,
        kmax=args.mf_kmax, r_model=args.mf_r_model,
        ped_nsamples=int(wcfg.ped_nsamples))
    overflow = int(wcfg.overflow)
    fallback = set(binfo["fallback"])
    print(f"[setup] MF templates {args.mf_templates} sha1 {binfo['sha1'][:10]} "
          f"run_type {binfo['run_type']} {binfo['n_module']} own, "
          f"{len(fallback)} fallback", flush=True)
    run = C.extract_run_number(args.evio_path)
    split = Path(args.evio_path).name.rsplit(".", 1)[-1]

    # module name -> list of example dicts; also count all qualifying pulses
    examples: dict[str, list] = {}
    n_qual: dict[str, int] = {}
    mtype: dict[str, str] = {}
    cnt = dict(trigger=0, fired=0, cluster=0, with_example=0)
    multi_pulse = p.hc_clusterer.get_config().seed_time_window > 0.0
    stats = C.LoopStats()

    def _progress(st):
        print(f"  [progress] phys={st.n_phys} trig={cnt['trigger']} "
              f"mult={cnt['fired']} 1cl={cnt['cluster']} "
              f"modules={len(examples)} {st.elapsed:.0f}s", flush=True)

    for fadc_evt, _ in C.iter_physics_events(
            p, stats, with_ssp=False, max_events=args.max_events,
            accept=lambda b: bool(b & accept_mask), progress=_progress,
            progress_every=20000):
        cnt["trigger"] += 1

        # One waveform pass per channel, reused by every cut below.
        chans = []
        n_fired = 0
        for tag, s, c, cd in C.iter_fadc_channels(fadc_evt):
            samples = np.asarray(cd.samples, dtype=np.uint16)
            ped, rms, peaks = p.wave_ana.analyze(samples)
            chans.append((tag, s, c, samples, ped, rms, peaks))
            if peaks:
                pk = peaks[0]
                if pk.quality == 0 and not pk.overflow and pk.height >= fired_hmin:
                    n_fired += 1
        if n_fired >= args.n_fired_lt:
            continue
        cnt["fired"] += 1

        # HyCal clustering, as in _common.reconstruct_hycal().
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
        if len(p.hc_clusterer.reconstruct_hits()) != args.n_clusters:
            continue
        cnt["cluster"] += 1

        evn = int(fadc_evt.info.event_number)
        took = False
        for mod, samples, ped, rms, peaks in hy:
            if len(peaks) != 1:
                continue
            pk = peaks[0]
            if pk.quality != 0 or pk.overflow:
                continue
            if pk.pos - args.pre_samples < 0 or \
                    pk.pos + args.post_samples + 1 > samples.shape[0]:
                continue
            h = float(pk.height)
            if h < args.height_min or h < args.height_rms_mult * rms \
                    or h >= args.height_lt:
                continue
            name = mod.name
            n_qual[name] = n_qual.get(name, 0) + 1
            lst = examples.setdefault(name, [])
            mtype[name] = mod.type.name
            if len(lst) < args.n_examples:
                lst.append(dict(evn=evn, midx=int(mod.index), samples=samples.copy(), ped=float(ped),
                                rms=float(rms), pos=int(pk.pos), height=h))
                took = True
        cnt["with_example"] += took

    # ---- matched-filter fit of every kept example ----
    flat = [(n_, e_) for n_ in examples for e_ in examples[n_]]
    if flat:
        Y = np.stack([e_["samples"] for _, e_ in flat])
        mid = np.array([e_["midx"] for _, e_ in flat], dtype=np.int64)
        res = rec.run(Y, mid, overflow)
        starts = np.searchsorted(res.p_ch, np.arange(len(flat) + 1))
        for i, (n_, e_) in enumerate(flat):
            sl = slice(starts[i], starts[i + 1])
            e_["mf"] = dict(ped=float(res.ped[i]), sigma=float(res.sigma[i]),
                            t0=res.p_t0[sl].copy(), amp=res.p_amp[sl].copy(),
                            chi2=(float(res.p_chi2[sl][0]) if starts[i + 1] > starts[i]
                                  else float("nan")),
                            npk=int(starts[i + 1] - starts[i]),
                            fallback=n_ in fallback)

    # ---- PDF ----
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages

    cut_txt = (f"trigger & {args.trigger_accept} (0x{accept_mask:x}) | "
               f"{args.n_clusters} HyCal cluster | fired channels < {args.n_fired_lt} "
               f"(clean 1st peak, h >= {fired_hmin:g} ADC)\n"
               f"waveform: 1 clean peak, {args.height_min:g} <= h < {args.height_lt:g} ADC "
               f"above pedestal, h >= {args.height_rms_mult:g} x ped rms, "
               f"{args.pre_samples}/{args.post_samples} samples pre/post peak\n"
               f"MF fit: {Path(args.mf_templates).name} (nsigma {args.mf_nsigma:g}, "
               f"min_amp {args.mf_min_amp:g}, kmax {args.mf_kmax}, r_model {args.mf_r_model:g})")
    out_pdf = Path(args.out_pdf)
    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    names = sorted(examples, key=_module_sort_key)
    with PdfPages(out_pdf) as pdf:
        for name in names:
            fig, axes = plt.subplots(2, 2, figsize=(11, 8.5))
            for ax, ex in zip(axes.flat, examples[name]):
                y = ex["samples"].astype(float)
                x = np.arange(y.size) * clk_ns
                ax.plot(x, y, "-", color="C0", lw=1.0, marker=".", ms=3)
                ax.axhline(ex["ped"], color="C3", ls="--", lw=1.0,
                           label=f"ped {ex['ped']:.1f} (rms {ex['rms']:.2f})")
                ax.plot([ex["pos"] * clk_ns], [ex["ped"] + ex["height"]], "v",
                        color="C1", ms=6, label=f"peak h {ex['height']:.1f} ADC")
                mf = ex["mf"]
                tt = np.arange(0.0, y.size * clk_ns, 0.5)
                fit = np.full_like(tt, mf["ped"])
                for a_, t0_ in zip(mf["amp"], mf["t0"]):
                    fit += a_ * bank.shape(np.array([ex["midx"]])[:, None],
                                           (tt - t0_)[None, :])[0]
                if mf["npk"]:
                    ax.plot(tt, fit, "-", color="C2", lw=1.6, label="MF spline fit")
                    ptxt = "\n".join(f"A{k + 1} {a_:.1f} ADC  t0 {t0_:.1f} ns"
                                     for k, (a_, t0_) in enumerate(zip(mf["amp"], mf["t0"])))
                    ftxt = (f"{ptxt}\nped {mf['ped']:.1f}  "
                            f"chi2/dof {mf['chi2']:.2f}  npk {mf['npk']}")
                else:
                    ftxt = f"MF: no pulse accepted  (ped {mf['ped']:.1f})"
                if mf["fallback"]:
                    ftxt += "\n[global template]"
                ax.text(0.97, 0.55, ftxt, transform=ax.transAxes, fontsize=7.5,
                        ha="right", va="top", color="C2",
                        bbox=dict(fc="w", ec="0.8", alpha=0.85))
                ax.set_title(f"{name}  event {ex['evn']}", fontsize=10)
                ax.set_xlabel("time (ns)")
                ax.set_ylabel("ADC")
                ax.legend(fontsize=8, loc="upper right")
                ax.grid(alpha=0.3)
            for ax in list(axes.flat)[len(examples[name]):]:
                ax.axis("off")
            fig.suptitle(f"Run {run}  split {split}  -  {name} ({mtype[name]}): "
                         f"{len(examples[name])} of {n_qual[name]} qualifying waveforms",
                         fontsize=12)
            fig.text(0.5, 0.005, cut_txt, ha="center", va="bottom", fontsize=7.5)
            fig.tight_layout(rect=(0, 0.07, 1, 0.96))
            pdf.savefig(fig)
            plt.close(fig)

    out_txt = out_pdf.with_suffix(".txt")
    with open(out_txt, "w", encoding="utf-8") as f:
        f.write(f"# run {run} split {split}  file {args.evio_path}\n")
        f.write("# cuts: " + cut_txt.replace("\n", "\n#       ") + "\n")
        f.write(f"physics_events\t{stats.n_phys}\n")
        f.write(f"trigger_pass\t{cnt['trigger']}\n")
        f.write(f"fired_lt_pass\t{cnt['fired']}\n")
        f.write(f"single_cluster_pass\t{cnt['cluster']}\n")
        f.write(f"events_with_qualifying_waveform_kept\t{cnt['with_example']}\n")
        f.write(f"qualifying_waveforms_total\t{sum(n_qual.values())}\n")
        f.write(f"modules_with_examples\t{len(names)}\n")
        f.write(f"pdf_pages\t{len(names)}\n")
        f.write(f"mf_templates\t{args.mf_templates}\n")
        f.write(f"mf_template_sha1\t{binfo['sha1']}\n")
        f.write(f"mf_template_run_type\t{binfo['run_type']}\n")
        f.write(f"mf_template_own_modules\t{binfo['n_module']}\n")
        f.write(f"mf_template_fallback_modules\t{len(fallback)}\n")
        if flat:
            allm = [e_["mf"] for _, e_ in flat]
            ch = np.array([m_["chi2"] for m_ in allm if m_["npk"]])
            npk = np.array([m_["npk"] for m_ in allm])
            nfb = sum(1 for m_ in allm if m_["fallback"])
            f.write("# MF fit quality over the example waveforms (not all qualifying)\n")
            f.write(f"mf_examples_fitted\t{len(allm)}\n")
            f.write(f"mf_npk_0_1_2_3plus\t{int((npk == 0).sum())}\t{int((npk == 1).sum())}"
                    f"\t{int((npk == 2).sum())}\t{int((npk >= 3).sum())}\n")
            f.write(f"mf_examples_on_global_template\t{nfb}\n")
            if ch.size:
                q = np.percentile(ch, [5, 25, 50, 75, 95])
                f.write("mf_chi2_dof_p5_p25_p50_p75_p95\t"
                        + "\t".join(f"{v:.2f}" for v in q) + "\n")
                f.write(f"mf_chi2_dof_gt3_fraction\t{float((ch > 3).mean()):.3f}\n")
        f.write("module\ttype\tn_examples\tn_qualifying\n")
        for name in names:
            f.write(f"{name}\t{mtype[name]}\t{len(examples[name])}\t{n_qual[name]}\n")
    print(f"[done] {out_pdf} ({len(names)} pages), {out_txt}", flush=True)
    print(f"[done] phys={stats.n_phys} trig={cnt['trigger']} mult={cnt['fired']} "
          f"1cl={cnt['cluster']} elapsed={stats.elapsed:.0f}s", flush=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
