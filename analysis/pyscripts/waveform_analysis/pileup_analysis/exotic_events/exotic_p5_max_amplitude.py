#!/usr/bin/env python3
"""exotic_p5_max_amplitude.py -- P5 of the exotic-waveform study: distribution of the maximum amplitude of the channels flagged f1..f5, per run type.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p5_max_amplitude.py [--jobs 8] [--plot-only] [--cache <pkl>] [--out <dir>]

Population: PbWO4 kept channels (A_main >= KEEP_ADC, gate and LMS veto as in pileup_rates.scan) with the uniform (P2) flags: f1, f3, f4 as stored (ch_ex_f1/3/4), f2 as stored
(ch_ex_f2) and not truncated (ch_trunc_main false; the t0 > 284 ns exclusion), not clipped (nclip == 0) and not f4 (f4 first, then f5 = nclip >= 1 and not f4 (saturated single; raw peak ~ 4095 - pedestal by construction, the fitted A is the extrapolation); P4 f2_eligible).  Reference row: all kept channels, as a random 1 % of the kept rows (--ref-frac).
Two definitions per channel:
  (a) raw peak = max over the 100 samples minus the pipeline pedestal (pileup_rates.pedestal, clipped mean of the first 30 samples), ADC.  Neither the full-scan pickles nor the npz
      store waveforms, so the flagged (event, module) rows are re-read from the EVIO files (decode only, no fit): rows whose event cannot be decoded have no raw peak.
  (b) fitted amplitude A_max = max over the k fitted pulses of ch_a (ADC).  The pipeline's main pulse is the largest-A pulse, so A_main == A_max by construction (checked).
E / E_beam = quantity x module_gain / beam energy (module_gain as pileup_rates: calibration factor x I_int); the raw peak is converted with the same ADC-to-MeV factor.
Data source per run type: H2 types from the 2-file sample npz (exotic_sample/, the only H2 files with the fitted amplitudes: the full-scan H2 pickles carry no amplitude);
Ta X17 from the 19 full-scan npz (exotic_full/).  Run 24937 is excluded (user: not a good run).
Categories are the exclusive ones of exotic_p4_module_rates.exclusive_category (priority f4 > f5 > f1 > f2 > f3).  1-sigma band: the 16 % and 84 % quantiles; median +- 1.4826 MAD is given alongside (cols mad_lo, mad_hi).
Outputs (--out): exotic_P5_max_amplitude_20261118.pdf (raw-peak page, fitted-amplitude page, tables) and .tsv.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import time
import pickle
import argparse
import tempfile
import multiprocessing as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import exotic_p1_quantities as p1                       # noqa: E402
from exotic_p1_quantities import TYPES, RT_LABEL, write_tsv, EXCLUDE_RUNS, EXCL_NOTE     # noqa: E402
from exotic_p4_module_rates import f2_eligible, exclusive_category, ORDER_TXT          # noqa: E402

STEM = "exotic_P5_max_amplitude_20261118"
FLAGS = (1, 2, 3, 4, 5)
ROWS = ("f1", "f2", "f3", "f4", "f5", "ref")
FNAME = {"f1": "f1 tail-baseline", "f2": "f2 wide", "f3": "f3 smooth-wide", "f4": "f4 clipped-multi", "f5": "f5 saturated single", "ref": "all kept (unflagged ref.)"}
SAT = 4095                                              # 12-bit digitiser full scale (max sample)
SEED = 20261118


def source_files(sample_dir, full_dir):
    """[(rt-agnostic npz path, kind)]: H2 types from the sample, Ta X17 from the full scan."""
    out = [(str(f), "sample") for f in sorted(Path(sample_dir).glob("run_*.npz"))]
    out += [(str(f), "full") for f in sorted(Path(full_dir).glob("run_*_s*.npz")) if ".part" not in f.name]
    return [(f, kind) for f, kind in out if int(Path(f).name.split("_")[1].split(".")[0]) not in EXCLUDE_RUNS]      # run 24937 excluded (user: not a good run)


def raw_peaks(z, run, midx, idx, evio_dir):
    """Raw peak (max sample - pipeline pedestal) and max sample of the rows idx of the scan npz z, by an EVIO decode-only re-read of their events (nan where the event is not decoded).  Shared with P6."""
    import event_display as ed
    from prad2py import dec
    import pileup_rates as pr
    evnum = z["ev_number"][z["ch_ev"].astype(np.int64)[idx]]
    raw = np.full(idx.size, np.nan)
    smax = np.full(idx.size, np.nan)
    want = {}                                            # evnum -> {midx: position in idx}
    for p, (en, m) in enumerate(zip(evnum.tolist(), midx[idx].tolist())):
        want.setdefault(en, {})[m] = p
    ctx = p1.get_ctx(run)
    left = set(want)
    for s in z["splits"]:
        if not left:
            break
        ch = ctx.open(str(Path(evio_dir) / f"prad_{run:06d}.evio.{int(s):05d}"))
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
                chs = ctx.channels(d["event"])
                mi = np.array([c[0] for c in chs])
                Y = np.stack([c[1] for c in chs])
                ped = pr.pedestal(Y)
                ym = Y.max(axis=1).astype(np.float64)
                where = {int(m): j for j, m in enumerate(mi)}
                for m, p in want[en].items():
                    j = where.get(m)
                    if j is not None:
                        raw[p] = ym[j] - ped[j]
                        smax[p] = ym[j]
        ch.close()
    return raw, smax


def work(task):
    """One npz: rows of the flagged channels and of the reference sample, with A_max, E/E_beam and the raw peak (re-read from EVIO)."""
    fn, kind, evio_dir, ref_frac = task
    import event_display as ed
    from prad2py import dec
    import pileup_rates as pr
    t_start = time.time()
    z = np.load(fn)
    run, rt = int(z["run"]), str(z["run_type"])
    if (kind == "sample") == ("Ta" in rt):               # H2 from the sample, Ta from the full scan only
        return None
    midx = z["ch_midx"].astype(np.int64)
    sel = z["mod_pwo4"][midx]
    k = z["ch_k"].astype(np.int64)
    a = z["ch_a"].astype(np.float64)
    valid = np.arange(a.shape[1])[None, :] < k[:, None]
    amax = np.where(valid, a, -np.inf).max(axis=1)
    t = z["ch_t"]
    main = np.argmax(np.where(np.isfinite(t), a, -1.0), axis=1)
    amain = a[np.arange(k.size), main]
    names = [str(n) for n in z["mod_name"]]
    cfg = pr.run_config(run)
    gain = pr.module_gain(cfg, rt, names)
    efac = gain[midx] / cfg["beam"]
    trunc = z["ch_trunc_main"].astype(bool)
    f = {i: z[f"ch_ex_f{i}"].astype(bool) & sel for i in (1, 2, 3, 4)}
    f[5] = (z["ch_nclip"] >= 1) & sel & ~f[4]              # saturated single (offline flag, P4 flags_of)
    f[2] &= f2_eligible(trunc, z["ch_nclip"]) & ~f[4]     # truncated and clipped records are not f2 (shared rule, P4)
    cat = exclusive_category(f)                          # exclusive categories f4 > f5 > f1 > f2 > f3 (shared with P4/P6); f below = the exclusive sets
    f = {i: (cat == i) & sel for i in FLAGS}
    rng = np.random.default_rng(SEED + run * 1000 + int(z["splits"][0]))
    ref = sel & (rng.random(k.size) < ref_frac)
    need = ref.copy()
    for i in FLAGS:
        need |= f[i]
    idx = np.nonzero(need)[0]
    raw, smax = raw_peaks(z, run, midx, idx, evio_dir)
    rec = dict(run=run, rt=rt, kind=kind, n_kept=int(sel.sum()), n_ref=int(ref.sum()), n_main_ne_max=int((np.abs(amain[sel] - amax[sel]) > 1e-3).sum()),
               n_noraw=int(np.isnan(raw).sum()), wall=time.time() - t_start,
               amax=amax[idx].astype(np.float32), e=(amax[idx] * efac[idx]).astype(np.float32), raw=raw.astype(np.float32), smax=smax.astype(np.float32),
               efac=efac[idx].astype(np.float32), ref=ref[idx], **{f"f{i}": f[i][idx] for i in FLAGS})
    print(f"[work] {os.path.basename(fn)} {rt}: {idx.size} rows, {rec['n_noraw']} without raw peak, {rec['wall']:.0f} s", flush=True)
    return rec


def gather(recs):
    """Per run type and row label: arrays (amax, e, raw, rawe, smax) and the list of source files."""
    D = {rt: {r: dict(amax=[], e=[], raw=[], rawe=[], smax=[]) for r in ROWS} for rt in TYPES}
    src = {rt: [] for rt in TYPES}
    chk = dict(n_main_ne_max=0, n_noraw=0, n_rows=0)
    for r in recs:
        if r is None:
            continue
        rt = r["rt"]
        src[rt].append((os.path.basename(r["kind"]), r["run"]))
        chk["n_main_ne_max"] += r["n_main_ne_max"]
        chk["n_noraw"] += r["n_noraw"]
        chk["n_rows"] += r["amax"].size
        for lab in ROWS:
            m = r["ref"] if lab == "ref" else r[lab]
            d = D[rt][lab]
            d["amax"].append(r["amax"][m]), d["e"].append(r["e"][m]), d["raw"].append(r["raw"][m]), d["smax"].append(r["smax"][m])
            d["rawe"].append(r["raw"][m] * r["efac"][m])
    for rt in TYPES:
        for lab in ROWS:
            D[rt][lab] = {q: (np.concatenate(v) if v else np.zeros(0, np.float32)) for q, v in D[rt][lab].items()}
    return D, src, chk


def quant(x):
    """N, median, q16, q84, MAD band (median -/+ 1.4826 MAD)."""
    x = x[np.isfinite(x)]
    if x.size == 0:
        return dict(N=0)
    med = float(np.median(x))
    s = 1.4826 * float(np.median(np.abs(x - med)))
    q16, q84 = np.percentile(x, [16, 84])
    return dict(N=int(x.size), median=med, q16=float(q16), q84=float(q84), mad_lo=med - s, mad_hi=med + s, sigma_mad=s)


DEFS = (("raw_adc", "raw", "raw peak [ADC]"), ("fit_adc", "amax", "fitted A_max [ADC]"), ("raw_efrac", "rawe", "raw peak E/E_beam"), ("fit_efrac", "e", "fitted A_max E/E_beam"))


def table(D):
    rows = []
    for rt in TYPES:
        for lab in ROWS:
            d = D[rt][lab]
            sat = d["smax"][np.isfinite(d["smax"])]
            for dn, q, _ in DEFS:
                st = quant(d[q])
                st.update(rt=rt, flag=lab, defn=dn, n_rows=int(d["amax"].size), frac_smax_ge_4095=(float((sat >= SAT).mean()) if sat.size else float("nan")),
                          frac_amax_gt_4095=(float((d["amax"] > SAT).mean()) if d["amax"].size else float("nan")),
                          source="sample npz (2 files/run)" if "H2" in rt else "full scan npz (19 files)")
                rows.append(st)
    return rows


def hist_page(pdf, plt, D, q, title, xlabel, edges, log=True):
    fig, ax = plt.subplots(5, 4, figsize=(16, 13), sharex=True)
    for j, rt in enumerate(TYPES):
        for i, fl in enumerate(("f1", "f2", "f3", "f4", "f5")):
            a = ax[i, j]
            x = D[rt][fl][q]
            x = x[np.isfinite(x) & (x > 0)]
            r = D[rt]["ref"][q]
            r = r[np.isfinite(r) & (r > 0)]
            if r.size:
                a.hist(r, bins=edges, weights=np.full(r.size, 1.0 / r.size), histtype="step", color="0.55", lw=1.0, label="all kept")
            if x.size:
                a.hist(x, bins=edges, weights=np.full(x.size, 1.0 / x.size), color="C%d" % i, alpha=0.75, label=fl)
                s = quant(x)
                a.axvline(s["median"], color="k", lw=1.4)
                for v in (s["q16"], s["q84"]):
                    a.axvline(v, color="k", lw=1.0, ls="--")
                a.text(0.97, 0.95, f"N={s['N']}\nmed {s['median']:.3g}\n16% {s['q16']:.3g}\n84% {s['q84']:.3g}", transform=a.transAxes, ha="right", va="top", fontsize=8.5,
                       family="monospace", bbox=dict(fc="w", ec="none", alpha=0.7))
            else:
                a.text(0.5, 0.5, "N=0", transform=a.transAxes, ha="center", fontsize=10)
            if log:
                a.set_xscale("log")
            if i == 0:
                a.set_title(RT_LABEL[rt], fontsize=11)
            if j == 0:
                a.set_ylabel(f"{FNAME[fl]}\nfraction / bin", fontsize=9)
            if i == 4:
                a.set_xlabel(xlabel, fontsize=9)
            a.tick_params(labelsize=8)
    fig.suptitle(title + "   (solid line: median, dashed: 16 % / 84 % quantiles; grey: all kept channels, same-normalisation)", fontsize=12)
    fig.tight_layout(rect=(0, 0, 1, 0.965))
    pdf.savefig(fig)
    plt.close(fig)


def tbl_lines(rows, dn, fmt):
    L = [f"{'':>26s}" + "".join(f"{RT_LABEL[rt]:>34s}" for rt in TYPES), f"{'':>26s}" + "".join(f"{'N   median [16%, 84%]':>34s}" for _ in TYPES)]
    for lab in ROWS:
        c = ""
        for rt in TYPES:
            r = next(x for x in rows if x["rt"] == rt and x["flag"] == lab and x["defn"] == dn)
            txt = f"{fmt(r['median'])} [{fmt(r['q16'])}, {fmt(r['q84'])}]" if r["N"] else "-"
            c += f"{r['N']:>8d}  {txt:>24s}"
        L.append(f"{FNAME[lab]:>26s}" + c)
    return L


def report(data, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    D, src, chk = data["D"], data["src"], data["chk"]
    rows = table(D)
    cols = ["rt", "flag", "defn", "N", "n_rows", "median", "q16", "q84", "mad_lo", "mad_hi", "sigma_mad", "frac_smax_ge_4095", "frac_amax_gt_4095", "source"]
    write_tsv(str(Path(out) / (STEM + ".tsv")), rows, cols)
    e_a = np.logspace(0.7, 4.3, 61)
    e_e = np.logspace(-3.5, 0.5, 61)
    with PdfPages(str(Path(out) / (STEM + ".pdf"))) as pdf:
        hist_page(pdf, plt, D, "raw", "Definition (a): raw peak = max sample - pedestal", "raw peak [ADC]", e_a)
        hist_page(pdf, plt, D, "amax", "Definition (b): largest fitted amplitude A_max (= A_main)", "A_max [ADC]", e_a)
        f3 = lambda v: f"{v:.0f}" if v >= 100 else f"{v:.1f}"          # noqa: E731
        f4 = lambda v: f"{v:.4f}" if v < 0.1 else f"{v:.3f}"           # noqa: E731
        L = ["Median [16 %, 84 % quantiles]: 1-sigma band = 16/84 % quantiles (median +- 1.4826 MAD is in the TSV).",
             "Flagged = PbWO4 kept channels, uniform P2 thresholds, EXCLUSIVE categories (priority " + ORDER_TXT + "; f2 not truncated and not clipped).", f"*** {EXCL_NOTE} ***", ""]
        L += ["(a) raw peak [ADC]"] + tbl_lines(rows, "raw_adc", f3) + ["", "(b) fitted A_max [ADC]  (A_main identical: the main pulse is the largest-A pulse)"] + tbl_lines(rows, "fit_adc", f3)
        L += ["", "(a) raw peak, E / E_beam"] + tbl_lines(rows, "raw_efrac", f4) + ["", "(b) fitted A_max, E / E_beam"] + tbl_lines(rows, "fit_efrac", f4)
        sat = [f"{FNAME[lab]:>26s}" + "".join(f"{next(x for x in rows if x['rt'] == rt and x['flag'] == lab and x['defn'] == 'raw_adc')['frac_smax_ge_4095']:>34.3f}" for rt in TYPES) for lab in ROWS]
        L += ["", "fraction of rows with max sample >= 4095 (saturated record)"] + sat
        L += ["", "Sources: H2 types, sample npz (exotic_sample, 2 files per run; the H2 pickles of exotic_full carry no amplitude); Ta X17, 19 full-scan npz; reference row = random "
              f"{data['ref_frac']:.0%} of kept rows.", f"Raw peak re-read from EVIO for the selected rows; {chk['n_noraw']} of {chk['n_rows']} rows without a raw peak (undecodable event).",
              f"A_main != A_max rows: {chk['n_main_ne_max']}.  E/E_beam: module_gain x amplitude / beam energy (raw peak with the same factor)."]
        fig = plt.figure(figsize=(16, 12))
        fig.text(0.03, 0.97, "P5: maximum amplitude of the flagged channels", fontsize=17, fontweight="bold", va="top")
        fig.text(0.03, 0.93, "\n".join(L), fontsize=8.2, va="top", family="monospace", linespacing=1.3)
        pdf.savefig(fig)
        plt.close(fig)
    print(f"[report] {Path(out) / (STEM + '.pdf')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir", default=str(p1.PYDIR / "output" / "pileup_rates" / "exotic_sample"))
    ap.add_argument("--full-dir", default=str(p1.PYDIR / "output" / "pileup_rates" / "exotic_full"))
    ap.add_argument("--evio-dir", default=p1.EVIO_DIR)
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--cache", default=str(Path(tempfile.gettempdir()) / "exotic_p5_cache.pkl"))
    ap.add_argument("--ref-frac", type=float, default=0.01)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    if a.plot_only:
        data = pickle.load(open(a.cache, "rb"))
    else:
        t0 = time.time()
        tasks = [(f, kind, a.evio_dir, a.ref_frac) for f, kind in source_files(a.sample_dir, a.full_dir)]
        with mp.get_context("fork").Pool(a.jobs) as pool:
            recs = pool.map(work, tasks, chunksize=1)
        D, src, chk = gather(recs)
        data = dict(D=D, src=src, chk=chk, ref_frac=a.ref_frac, elapsed=time.time() - t0)
        pickle.dump(data, open(a.cache, "wb"))
        print(f"[compute] {data['elapsed']:.0f} s")
    report(data, a.out)


if __name__ == "__main__":
    main()
