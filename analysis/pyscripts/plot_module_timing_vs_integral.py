#!/usr/bin/env python3
"""plot_module_timing_vs_integral.py - per-module t_module - t_veto vs fitted pulse integral.

Checks whether dt = t_module - t_veto depends on the module pulse integral (energy).  One page
per PbWO4 module, four panels (0.73 GeV H2, 2.24 GeV H2, 3.49 GeV H2, X17 Ta): 2D histogram of
dt (vertical) vs the fitted integral (horizontal, linear) with the median and 16-84 % band per
integral bin, a weighted linear fit of the medians, and one summary page of the slopes first.
X17 Ta has no veto readout: its panel shows dt = t_module - t_seed (plot_module_timing.seed_refs:
highest-m_E PbWO4 pulse of the event other than the module, m_E > --ref-min-e; for the seed module
the second-highest, i.e. t_seed - t_2nd).

DEFINITIONS (event selection, t_evt, D0, run exclusion and veto pick as plot_module_timing.py
--ref veto, whose functions are imported)
  Fit integral      Eu = A * I_int(phase) [* f(A)] (ADC*samples) of the module's dominant mf
                    pulse, recovered exactly as m_E / cal_factor (energize() is linear; factor
                    from the calibration file named in each npz).  Rows where the dump fell back
                    to WaveAnalyzer (m_fb, raw-sum integral) are excluded.
  Veto alignment    --veto-align counter|centroid|type as plot_module_timing.py (shared function
                    measure_counter_offsets; default counter = module-referenced shifts c_k,
                    v_t' = v_t - c_k).  Gate |v_t' - t_evt - D0| < --veto-gate, integral >
                    --veto-min-int, largest integral wins.  --veto-counter Vk: only pulses of
                    that counter.  Shifts and both measurements in the _veto_offsets JSON.
  Panel window      dt within +-W (--dt-half, 10 ns) of the module's peak (smoothed argmax of the
                    0.2 ns histogram over the plot_module_timing histogram range); N = entries in it.
  Integral range    per run type: from the 0.5 % percentile of the in-window integrals of the page
                    modules up to 1.2 x E_beam / (median PbWO4 cal_factor of the type), 150 linear
                    bins.  Profile: bins with >= --min-bin entries.
  Slope             weighted least squares of median dt vs mean integral of the bin, weights
                    1/err^2, err = 1.2533 * (p84 - p16) / 2 / sqrt(n); error scaled by
                    sqrt(chi2/ndf) when > 1; quoted in ns per 1000 ADC*samples.
  Pages             modules with N > --min-n in at least one run type, N counted exactly as the
                    plot_module_timing.py --ref veto histograms (all rows, same ranges; X17 included),
                    ordered by id, so both PDFs align page for page.

Outputs (--out-dir): module_timing_vs_integral{sfx}.pdf, .tsv, _veto_offsets{sfx}.json, where
sfx = "" | "_V<k>" (--veto-counter) | "_typealign" / "_centroidalign" (--veto-align).

Usage:
  MPLBACKEND=Agg python3 plot_module_timing_vs_integral.py [--veto-align counter|centroid|type]
      [--veto-counter any|V1|V2|V3|V4] [--dump-dir DIR] [--out-dir DIR] [--jobs 12]
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plot_module_timing as PMT  # noqa: E402
from plot_module_timing import (LABEL, NMAX_JOBS, RUN_NAME, TAGS, TYPE_OF_RUN, W_BIN, centres, f,  # noqa: E402
                                fit_core, pick_veto, seed_refs)

NIB = 150              # integral bins (display and profile)
E_BEAM = {"0p73": 730.0, "2p24": 2240.0, "3p49": 3490.0, "x17": 2240.0}   # MeV, integral axis upper edge


# ---- workers ----------------------------------------------------------------

def pass2(job):
    """Rows with a reference time: per-module counts in the plot_module_timing histogram range
    (all rows) and, for non-fallback rows in range, (module, dt, fit integral, module is seed)."""
    path, cfg, d0, off, seed, lo, nb = job
    d, n_ev, sel, m_ev, m_idx, t, E, t_evt, t_loo, v = PMT._load(path, cfg)
    npw = cfg["npw"]
    out = dict(run=int(d["run"]), npage=np.zeros(npw, np.int64), p=np.zeros(0, np.int16),
               dt=np.zeros(0, np.float32), eu=np.zeros(0, np.float32), seed=np.zeros(0, bool), n_fb=0)
    is_seed = np.zeros(len(m_ev), bool)
    if seed:
        t_ref, is_seed = seed_refs(m_ev, t, E, sel, cfg["ref_min_e"])
    elif np.isfinite(d0):
        vt = v["t"] - off[v["ch"]]
        keep = np.ones(len(vt), bool) if cfg["counter"] == 0 else (v["ch"] == cfg["counter"])
        t_ref = pick_veto(v["ev"][keep], v["ch"][keep], vt[keep], v["i"][keep], sel, t_evt, d0,
                          cfg["gate"], cfg["veto_min_int"], n_ev)[0][m_ev]
    else:
        return out
    pmap = cfg["pmap"]
    r = np.flatnonzero(sel[m_ev] & (pmap[m_idx] >= 0) & np.isfinite(t_ref))
    dt = t[r] - t_ref[r]
    b = np.floor((dt - lo) / W_BIN).astype(np.int64)
    w = (b >= 0) & (b < nb)
    r, dt = r[w], dt[w]
    out["npage"] = np.bincount(pmap[m_idx[r]], minlength=npw).astype(np.int64)
    cal = {m["name"]: float(m["factor"]) for m in json.load(open(str(d["calib_file"])))}
    cf = np.array([cal.get(str(n), 0.0) for n in d["module_names"]])
    fb = d["m_fb"][r]
    out["n_fb"] = int(fb.sum())
    k = ~fb & (cf[m_idx[r]] > 0)
    r, dt = r[k], dt[k]
    out["p"] = pmap[m_idx[r]].astype(np.int16)
    out["dt"] = dt.astype(np.float32)
    out["eu"] = (E[r] / cf[m_idx[r]]).astype(np.float32)
    out["seed"] = is_seed[r]
    return out


# ---- per-module analysis ----------------------------------------------------

def window_centre(dt, lo, nb):
    h = np.bincount(np.clip(((dt - lo) / W_BIN).astype(np.int64), 0, nb - 1), minlength=nb)
    k = int(round(1.0 / W_BIN)) | 1
    sm = np.convolve(h, np.ones(k) / k, mode="same")
    return lo + (int(np.argmax(sm)) + 0.5) * W_BIN


def profile(dt, eu, edges, min_bin):
    """Per integral bin: n, mean integral, median, p16, p84, median error."""
    b = np.digitize(eu, edges) - 1
    rows = []
    for k in range(len(edges) - 1):
        s = dt[b == k]
        if s.size < min_bin:
            continue
        p16, med, p84 = np.percentile(s, [16, 50, 84])
        rows.append((s.size, float(eu[b == k].mean()), med, p16, p84,
                     1.2533 * 0.5 * (p84 - p16) / np.sqrt(s.size)))
    return np.array(rows).reshape(-1, 6)


def linfit(x, y, e):
    """Weighted straight line; returns slope, slope_err (scaled by sqrt(chi2/ndf) if > 1), chi2/ndf."""
    if len(x) < 3:
        return np.nan, np.nan, np.nan
    w = 1.0 / np.maximum(e, 1e-6) ** 2
    A = np.vstack([np.ones_like(x), x]).T
    cov = np.linalg.inv(A.T @ (A * w[:, None]))
    p = cov @ (A.T @ (w * y))
    chi2 = float(np.sum(w * (y - A @ p) ** 2))
    ndf = len(x) - 2
    serr = np.sqrt(cov[1, 1]) * max(1.0, np.sqrt(chi2 / ndf))
    return float(p[1]), float(serr), chi2 / ndf


def analyse(dt, eu, sd, lo, nb, half, edges, min_bin):
    c = window_centre(dt, lo, nb)
    m = np.abs(dt - c) < half
    dt, eu = dt[m], eu[m]
    res = dict(N=int(m.sum()), c=c, dt=dt, eu=eu, prof=np.zeros((0, 6)), slope=np.nan,
               slope_err=np.nan, chi2=np.nan, med=float(np.median(dt)) if dt.size else np.nan,
               f_seed=float(sd[m].mean()) if dt.size else np.nan)
    pr = profile(dt, eu, edges, min_bin)
    res["prof"] = pr
    if len(pr):
        s, se, chi = linfit(pr[:, 1] / 1000.0, pr[:, 2], pr[:, 5])
        res.update(slope=s, slope_err=se, chi2=chi)
    return res


# ---- plotting ---------------------------------------------------------------

def module_page(pdf, plt, geo_row, R, edges, half, min_draw, kind):
    """kind[tag]: 'veto', 'seed' (t_m - t_seed, no veto) or None (no data)."""
    from matplotlib.colors import LogNorm
    nm, mid, x, y = geo_row
    fig, axs = plt.subplots(2, 2, figsize=(11, 8.5))
    fig.suptitle(f"{nm}   id {mid}   x = {x:.2f} mm, y = {y:.2f} mm", fontsize=13, weight="bold")
    for ax, tag in zip(axs.ravel(), TAGS):
        ax.set_title(LABEL[tag], fontsize=10)
        ax.tick_params(labelsize=8)
        r = R.get(tag)
        msg = None
        if kind[tag] is None:
            msg = "no veto readout" if tag == "x17" else "no data"
        elif r is None or r["N"] < min_draw:
            msg = f"N = {0 if r is None else r['N']} (< {min_draw}), not drawn"
        if msg:
            ax.set_xticks([])
            ax.set_yticks([])
            ax.text(0.5, 0.5, msg, transform=ax.transAxes, ha="center", va="center", fontsize=9,
                    color="0.35")
            continue
        ed = edges[tag]
        yb = r["c"] + np.linspace(-half, half, int(round(2 * half / W_BIN)) + 1)
        h = np.histogram2d(r["eu"], r["dt"], bins=(ed, yb))[0]
        h = np.ma.masked_less(h, 1)
        if h.count():
            ax.pcolormesh(ed, yb, h.T, cmap="viridis", norm=LogNorm(vmin=1, vmax=max(h.max(), 2)),
                          rasterized=True)
        pr = r["prof"]
        if len(pr):
            ax.plot(pr[:, 1], pr[:, 2], "o-", color="crimson", ms=2.5, lw=1.2, label="median")
            ax.plot(pr[:, 1], pr[:, 3], "--", color="crimson", lw=0.8, label="16-84 %")
            ax.plot(pr[:, 1], pr[:, 4], "--", color="crimson", lw=0.8)
        txt = [f"N = {r['N']}", f"median dt = {r['med']:.2f} ns"]
        if np.isfinite(r["slope"]):
            i0 = np.average(pr[:, 1], weights=1 / pr[:, 5] ** 2)
            m0 = np.average(pr[:, 2], weights=1 / pr[:, 5] ** 2)
            xs = np.array([pr[0, 1], pr[-1, 1]])
            ax.plot(xs, m0 + r["slope"] * (xs - i0) / 1000.0, "--", color="darkorange", lw=1.4, label="lin. fit")
            dd = pr[-1, 2] - pr[0, 2]
            de = np.hypot(pr[-1, 5], pr[0, 5])
            txt += [f"slope = {r['slope']:+.3f} +- {r['slope_err']:.3f} ns/1000",
                    f"dt(high) - dt(low) = {dd:+.2f} +- {de:.2f} ns",
                    f"  (bins at {pr[0, 1]:.0f} and {pr[-1, 1]:.0f})"]
        else:
            txt.append("< 3 bins with enough entries, no fit")
        if kind[tag] == "seed":
            txt.append(f"module = seed in {100 * r['f_seed']:.1f} %\n(those: t_seed - t_2nd)")
        ax.text(0.97, 0.97, "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
        ax.set_xlim(ed[0], ed[-1])
        ax.set_ylim(yb[0], yb[-1])
        ax.set_xlabel("fit integral A * I_int (ADC * samples)", fontsize=8)
        ax.set_ylabel("t_module - t_seed (ns), no veto" if kind[tag] == "seed" else "t_module - t_veto (ns)",
                      fontsize=8)
        ax.legend(fontsize=6.5, loc="lower right", framealpha=0.8)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


def robust_outliers(names, s, se):
    """Median, robust spread (1.4826 MAD) and modules with |s - median| > 3 x spread AND > 3 x own error."""
    ok = np.isfinite(s)
    if ok.sum() < 5:
        return np.nan, np.nan, []
    med = float(np.median(s[ok]))
    rs = 1.4826 * float(np.median(np.abs(s[ok] - med)))
    rows = [(names[i], s[i], se[i]) for i in np.flatnonzero(ok)
            if rs > 0 and abs(s[i] - med) > 3 * max(rs, se[i])]
    return med, rs, rows


def summary_page(pdf, plt, gp, page_mods, S, tags, kind, seed_lines, title, foot):
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(len(tags), 3, width_ratios=[3.0, 1.3, 1.6], hspace=0.55, wspace=0.3,
                          left=0.07, right=0.98, top=0.92, bottom=0.19)
    for j, tag in enumerate(tags):
        s, se = S[tag]["slope"][page_mods], S[tag]["slope_err"][page_mods]
        ok = np.isfinite(s)
        a0 = fig.add_subplot(gs[j, 0])
        a0.errorbar(gp["id"][page_mods][ok], s[ok], yerr=se[ok], fmt="o", ms=2, lw=0.6, color="C0")
        a0.axhline(0, color="0.5", lw=0.6)
        a0.set_ylabel("slope (ns / 1000)", fontsize=8)
        what = "t_module - t_seed (no veto)" if kind[tag] == "seed" else "t_module - t_veto"
        a0.set_title(f"{LABEL[tag]}: slope of median {what} vs fit integral ({ok.sum()} modules)", fontsize=8.5)
        a0.tick_params(labelsize=7)
        med, rs, orow = robust_outliers(gp["name"][page_mods], s, se)
        if ok.sum():
            lo, hi = np.percentile(s[ok], [1, 99])
            pad = 0.15 * (hi - lo if hi > lo else 1.0) + np.median(se[ok])
            a0.set_ylim(lo - pad, hi + pad)
        if j == len(tags) - 1:
            a0.set_xlabel("module id", fontsize=8)
        a1 = fig.add_subplot(gs[j, 1])
        if ok.sum():
            a1.hist(s[ok], bins=25, color="C0", alpha=0.8)
            if np.isfinite(med):
                a1.axvline(med, color="crimson", lw=1)
        a1.set_xlabel("slope (ns / 1000 ADC*samples)", fontsize=7)
        a1.tick_params(labelsize=7)
        a2 = fig.add_subplot(gs[j, 2])
        a2.axis("off")
        L = [f"median slope {f(med, 3)} ns/1000",
             f"robust spread (1.4826 MAD) {f(rs, 3)}",
             f"median stat. err {f(np.median(se[ok]) if ok.sum() else np.nan, 3)}",
             f"median dt(high)-dt(low) {f(np.nanmedian(S[tag]['ddt'][page_mods]) if ok.sum() else np.nan, 2)} ns",
             f"|slope - med| > 3 x max(spread, err): {len(orow)}"]
        nmax = 2 if kind[tag] == "seed" else 7
        L += [f"  {n:<6}{v:+8.3f} +- {e:.3f}" for n, v, e in orow[:nmax]]
        if len(orow) > nmax:
            L.append(f"  ... {len(orow) - nmax} more (see TSV)")
        if kind[tag] == "seed":
            L += seed_lines
        a2.text(0.0, 1.0, "\n".join(L), va="top", ha="left", fontsize=6.5, family="monospace",
                transform=a2.transAxes)
    fig.suptitle(title, fontsize=11)
    import textwrap
    fig.text(0.02, 0.01, "\n".join(w for ln in foot for w in textwrap.wrap(ln, 175, subsequent_indent="  ")),
             fontsize=6.3, family="monospace", va="bottom")
    pdf.savefig(fig)
    plt.close(fig)


# ---- main -------------------------------------------------------------------

def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--dump-dir", default=str(HERE / "output/module_timing_20261025/dump"))
    ap.add_argument("--out-dir", default=str(HERE / "output/module_timing_20261025"))
    ap.add_argument("--veto-align", choices=("counter", "centroid", "type"), default="counter")
    ap.add_argument("--veto-counter", choices=("any", "V1", "V2", "V3", "V4"), default="any")
    ap.add_argument("--trigger-mask", type=lambda s: int(s, 0), default=0xF00)
    ap.add_argument("--ref-min-e", type=float, default=10.0)
    ap.add_argument("--ref-min-n", type=int, default=3)
    ap.add_argument("--veto-min-int", type=float, default=2000.0)
    ap.add_argument("--veto-gate", type=float, default=10.0)
    ap.add_argument("--dt-half", type=float, default=10.0, help="dt half window around the module peak (ns)")
    ap.add_argument("--min-n", type=int, default=10000)
    ap.add_argument("--min-draw", type=int, default=500, help="draw a panel with at least this N")
    ap.add_argument("--min-bin", type=int, default=50, help="entries per integral bin for the profile")
    ap.add_argument("--geometry-run", type=int, default=25320)
    ap.add_argument("--jobs", type=int, default=NMAX_JOBS)
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    t_start = time.time()
    nproc = max(1, min(args.jobs, NMAX_JOBS))
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    counter = 0 if args.veto_counter == "any" else int(args.veto_counter[1])
    sfx = ("" if counter == 0 else f"_V{counter}") + ("" if args.veto_align == "counter" else f"_{args.veto_align}align")
    pdf_path = out_dir / f"module_timing_vs_integral{sfx}.pdf"
    tsv_path = out_dir / f"module_timing_vs_integral{sfx}.tsv"
    off_path = out_dir / f"module_timing_vs_integral_veto_offsets{sfx}.json"

    files, runs_of = [], {t: [] for t in TAGS}
    for fp in sorted(glob.glob(os.path.join(args.dump_dir, "*.npz"))):
        mt = re.search(RUN_NAME, fp)
        run = int(mt.group(1)) if mt else None
        if run not in TYPE_OF_RUN:
            continue
        files.append(fp)
        if run not in runs_of[TYPE_OF_RUN[run]]:
            runs_of[TYPE_OF_RUN[run]].append(run)
    if not files:
        sys.exit(f"no usable npz in {args.dump_dir}")
    for t in TAGS:
        runs_of[t].sort()

    geo = PMT.load_geometry(args.geometry_run)
    names = np.load(files[0])["module_names"]
    if [str(n) for n in names] != [str(n) for n in geo["name"]]:
        sys.exit("module_names of the dump differ from the hycal geometry")
    pidx = np.flatnonzero(geo["pwo4"])
    npw = len(pidx)
    pmap = np.full(len(names), -1, np.int64)
    pmap[pidx] = np.arange(npw)
    gp = {k: geo[k][pidx] for k in ("name", "id", "x", "y")}
    cfg = dict(mask=args.trigger_mask, ref_min_e=args.ref_min_e, ref_min_n=args.ref_min_n,
               veto_min_int=args.veto_min_int, gate=args.veto_gate, off=None, pmap=pmap, npw=npw,
               counter=counter, names_hash=hashlib.md5("\n".join(map(str, names)).encode()).hexdigest())

    # ---- pass 1: D0, run exclusion, counter offsets (plot_module_timing functions) ----
    t0 = time.time()
    p1 = {}
    for r in PMT._run_pool(PMT.pass1, [(fp, cfg) for fp in files], nproc, "pass1"):
        q = p1.setdefault(r["run"], dict(h_d=0, h_dc=0))
        q["h_d"] = q["h_d"] + r["h_d"]
        q["h_dc"] = q["h_dc"] + r["h_dc"]
    d0, d0fit, d0run, d0_stat, d0_hist = PMT.veto_d0({r: p1[r]["h_d"] for r in p1}, runs_of, args.veto_gate)
    off, cinfo = PMT.measure_counter_offsets(files, cfg, runs_of, TYPE_OF_RUN, d0, d0_hist,
                                             {r: p1[r]["h_dc"] for r in p1}, args.veto_align, nproc)
    seed_tags = [t for t in TAGS if runs_of[t] and not any(np.isfinite(d0_hist[r]) for r in runs_of[t])]
    vlo, vnb = PMT.veto_windows(runs_of, d0_hist, seed_tags)
    json.dump(dict(veto_align=args.veto_align, veto_counter=args.veto_counter, veto_min_int=args.veto_min_int,
                   definition="aligned v_t' = v_t - applied[k] (ns, V1..V4).  module_referenced: unaligned "
                              "dt = t_m - t_veto medians per (module, counter) fitted as a_m + b_k, shift = -b_k, "
                              "entry-weighted mean 0.  vt_minus_tevt_centroid: per-counter coincidence centroid "
                              "(used by --veto-align centroid as centroid - D0).",
                   run_types={LABEL[t]: cinfo[t] for t in TAGS if cinfo[t]["runs_used"]}),
              open(off_path, "w"), indent=1)
    print(f"pass 1 + counter alignment done in {time.time() - t0:.0f} s", flush=True)
    for t in TAGS:
        if cinfo[t]["runs_used"]:
            print(f"  {LABEL[t]}: D0 {f(d0[t], 2)} ns, runs {cinfo[t]['runs_used']}; shifts "
                  + " ".join(f"V{c} {v:+.2f}" for c, v in enumerate(cinfo[t]["applied"], 1))
                  + f"; residual rms {cinfo[t]['module_referenced'].get('residual_rms_ns')} ns", flush=True)

    # ---- pass 2: per-row dt and fit integral ----
    t0 = time.time()
    run_of = lambda fp: int(re.search(RUN_NAME, fp).group(1))
    jobs = [(fp, cfg, d0_hist[run_of(fp)], off[TYPE_OF_RUN[run_of(fp)]], TYPE_OF_RUN[run_of(fp)] in seed_tags,
             vlo[TYPE_OF_RUN[run_of(fp)]], vnb[TYPE_OF_RUN[run_of(fp)]]) for fp in files]
    acc = {t: dict(p=[], dt=[], eu=[], seed=[], n_fb=0, npage=np.zeros(npw, np.int64)) for t in TAGS}
    for r in PMT._run_pool(pass2, jobs, nproc, "pass2"):
        a = acc[TYPE_OF_RUN[r["run"]]]
        for k in ("p", "dt", "eu", "seed"):
            a[k].append(r[k])
        a["n_fb"] += r["n_fb"]
        a["npage"] += r["npage"]
    print(f"pass 2 done in {time.time() - t0:.0f} s", flush=True)
    kind = {t: ("seed" if t in seed_tags else "veto") if acc[t]["npage"].sum() else None for t in TAGS}
    tags = [t for t in TAGS if kind[t]]
    D = {}
    for t in tags:
        p = np.concatenate(acc[t]["p"]).astype(np.int64)
        o = np.argsort(p, kind="stable")
        D[t] = dict(dt=np.concatenate(acc[t]["dt"])[o], eu=np.concatenate(acc[t]["eu"])[o],
                    seed=np.concatenate(acc[t]["seed"])[o], start=np.searchsorted(p[o], np.arange(npw + 1)))
        acc[t]["p"] = acc[t]["dt"] = acc[t]["eu"] = acc[t]["seed"] = None
        print(f"  {LABEL[t]}: {len(p)} rows, {acc[t]['n_fb']} WaveAnalyzer-fallback rows excluded", flush=True)
    page_mods = np.flatnonzero(np.any([acc[t]["npage"] > args.min_n for t in tags], axis=0))
    page_mods = page_mods[np.argsort(gp["id"][page_mods])]
    print(f"modules with a page: {len(page_mods)} ("
          + ", ".join(f"{LABEL[t]} {int((acc[t]['npage'] > args.min_n).sum())}" for t in tags) + ")", flush=True)

    # ---- per-type integral range, per-module analysis ----
    edges, R, S, seed_lines = {}, {}, {}, []
    for t in tags:
        st = D[t]["start"]
        sl = lambda i: slice(st[i], st[i + 1])
        get = lambda i: (D[t]["dt"][sl(i)].astype(np.float64), D[t]["eu"][sl(i)].astype(np.float64),
                         D[t]["seed"][sl(i)])
        pool = []
        for i in page_mods:
            dt, eu, _ = get(i)
            if dt.size:
                c = window_centre(dt, vlo[t], vnb[t])
                pool.append(eu[np.abs(dt - c) < args.dt_half])
        pool = np.concatenate(pool) if pool else np.array([0.0, 1.0])
        f0 = np.load(next(fp for fp in files if TYPE_OF_RUN[run_of(fp)] == t))
        cal = {m["name"]: float(m["factor"]) for m in json.load(open(str(f0["calib_file"])))}
        cfs = np.array([cal.get(str(n), 0.0) for n in names[pidx]])
        lo, p995 = np.percentile(pool, [0.5, 99.5])
        hi = 1.2 * E_BEAM[t] / float(np.median(cfs[cfs > 0]))
        edges[t] = np.linspace(lo, hi, NIB + 1)
        print(f"  {LABEL[t]}: integral range [{lo:.0f}, {hi:.0f}] ADC*samples (99.5 % of the pool {p995:.0f}, "
              f"{100 * np.mean(pool > hi):.2f} % above)", flush=True)
        S[t] = {k: np.full(npw, np.nan) for k in ("slope", "slope_err", "chi2", "ddt", "ddt_err", "mu", "sigma")}
        for i in page_mods:
            dt, eu, sd = get(i)
            if dt.size == 0:
                continue
            r = analyse(dt, eu, sd, vlo[t], vnb[t], args.dt_half, edges[t], args.min_bin)
            R[(t, i)] = r
            pr = r["prof"]
            if r["N"] >= args.min_draw and np.isfinite(r["slope"]):
                S[t]["slope"][i], S[t]["slope_err"][i], S[t]["chi2"][i] = r["slope"], r["slope_err"], r["chi2"]
                S[t]["ddt"][i] = pr[-1, 2] - pr[0, 2]
                S[t]["ddt_err"][i] = np.hypot(pr[-1, 5], pr[0, 5])
            if kind[t] == "seed" and acc[t]["npage"][i] > args.min_n:
                # same histogram and fit as plot_module_timing.py (full range, 0.2 ns), so both PDFs agree
                hb = np.bincount(np.floor((D[t]["dt"][sl(i)] - vlo[t]) / W_BIN).astype(np.int64), minlength=vnb[t])
                fr = fit_core(centres(vlo[t], W_BIN, vnb[t]), hb[:vnb[t]])
                if fr["ok"]:
                    S[t]["mu"][i], S[t]["sigma"][i] = fr["mu"], fr["sigma"]
        if kind[t] == "seed":
            mu, sg = S[t]["mu"][page_mods], S[t]["sigma"][page_mods]
            ok = np.isfinite(mu)
            if ok.sum():
                mm = float(np.median(mu[ok]))
                far = [(gp["name"][page_mods][k], mu[k] - mm) for k in np.flatnonzero(ok) if abs(mu[k] - mm) > 3.0]
                seed_lines = ["", "core fit (as module_timing.pdf):",
                              f" {int(ok.sum())} mod., median sigma {np.median(sg[ok]):.2f} ns",
                              f" |mean - med| > 3 ns: {len(far)} (TSV)"]
                print("  X17 seed summary: " + " | ".join(seed_lines[2:4]) + "; " + ", ".join(f"{n} {d:+.2f}" for n, d in far),
                      flush=True)

    # ---- TSV ----
    cols = ["module", "id", "run_type", "reference", "N", "n_bins", "slope_ns_per_1000", "slope_err", "chi2_ndf",
            "int_low_bin", "dt_low", "int_high_bin", "dt_high", "dt_high_minus_low", "dt_hml_err",
            "median_dt", "mean_fit", "sigma_fit", "f_module_is_seed", "int_lo", "int_hi"]
    with open(tsv_path, "w") as fh:
        fh.write("\t".join(cols) + "\n")
        for i in page_mods:
            for t in tags:
                r = R.get((t, i))
                if r is None:
                    continue
                pr = r["prof"]
                hasp = len(pr) > 0
                row = [gp["name"][i], gp["id"][i], LABEL[t], "t_seed" if kind[t] == "seed" else "t_veto",
                       r["N"], len(pr),
                       f(S[t]["slope"][i], 4), f(S[t]["slope_err"][i], 4), f(S[t]["chi2"][i], 2),
                       f(pr[0, 1] if hasp else np.nan, 0), f(pr[0, 2] if hasp else np.nan, 3),
                       f(pr[-1, 1] if hasp else np.nan, 0), f(pr[-1, 2] if hasp else np.nan, 3),
                       f(S[t]["ddt"][i], 3), f(S[t]["ddt_err"][i], 3), f(r["med"], 3),
                       f(S[t]["mu"][i], 3), f(S[t]["sigma"][i], 3),
                       f(r["f_seed"], 4) if kind[t] == "seed" else "nan",
                       f(edges[t][0], 0), f(edges[t][-1], 0)]
                fh.write("\t".join(map(str, row)) + "\n")

    # ---- PDF ----
    t0 = time.time()
    plt = PMT._mpl()
    from matplotlib.backends.backend_pdf import PdfPages
    vtags = [t for t in tags if kind[t] == "veto"]
    vlab = "any of V1-V4" if counter == 0 else f"V{counter} only"
    alab = {"counter": "counters aligned (module-referenced)", "centroid": "counters aligned (centroid)",
            "type": "counters not aligned"}[args.veto_align]
    foot = [f"dt = t_module (mf 50 % leading edge, dominant pulse) - t_veto (in-time {vlab}, {alab}; integral > "
            f"{args.veto_min_int:g}, |v_t' - t_evt - D0| < {args.veto_gate:g} ns, largest integral).",
            "Veto shifts c_k (ns, v_t' = v_t - c_k): " + "; ".join(
                f"{LABEL[t]} " + " ".join(f"V{c} {v:+.2f}" for c, v in enumerate(cinfo[t]["applied"], 1))
                for t in vtags) + f"  ({off_path.name})",
            f"X17 Ta (no veto): dt = t_module - t_seed, seed = highest-m_E PbWO4 pulse of the event other than the module "
            f"(m_E > {args.ref_min_e:g} MeV, pulse level); seed module: t_seed - t_2nd.",
            f"Fit integral = m_E / cal_factor (ADC*samples), WaveAnalyzer-fallback rows excluded.  Window: module peak +- "
            f"{args.dt_half:g} ns.  Profile: median, 16-84 % per bin (>= {args.min_bin} entries), {NIB} bins.",
            "Slope: weighted linear fit of the medians, error scaled by sqrt(chi2/ndf) if > 1.  Outliers: |slope - median| > "
            "3 x max(1.4826 MAD, own error).  D0: " + ", ".join(f"{LABEL[t]} {d0[t]:.2f}" for t in vtags) + " ns; runs: "
            + "; ".join(f"{LABEL[t]} {cinfo[t]['runs_used']}" for t in vtags) + ".",
            "Integral ranges: " + ", ".join(f"{LABEL[t]} [{edges[t][0]:.0f}, {edges[t][-1]:.0f}]" for t in tags)
            + f".  Pages: {len(page_mods)} modules with N > {args.min_n} in a panel (as module_timing.pdf)."]
    with PdfPages(pdf_path) as pdf:
        summary_page(pdf, plt, gp, page_mods, S, tags, kind, seed_lines,
                     f"dt vs fit integral: per-module slopes  (veto {vlab}, {alab}; X17: t_module - t_seed)", foot)
        for i in page_mods:
            module_page(pdf, plt, (gp["name"][i], gp["id"][i], gp["x"][i], gp["y"][i]),
                        {t: R.get((t, i)) for t in tags}, edges, args.dt_half, args.min_draw, kind)
    print(f"pdf written in {time.time() - t0:.0f} s: {pdf_path}", flush=True)
    for t in tags:
        s = S[t]["slope"][page_mods]
        ok = np.isfinite(s)
        med, rs, orow = robust_outliers(gp["name"][page_mods], s, S[t]["slope_err"][page_mods])
        print(f"  {LABEL[t]}: {ok.sum()} slopes, median {f(med, 3)}, spread {f(rs, 3)} ns/1000, median "
              f"dt(high)-dt(low) {f(np.nanmedian(S[t]['ddt'][page_mods]) if ok.sum() else np.nan, 2)} ns; "
              f"{len(orow)} outliers: " + ", ".join(f"{n} {v:+.3f}+-{e:.3f}" for n, v, e in orow), flush=True)
    print(f"total {time.time() - t_start:.0f} s; outputs: {pdf_path}, {tsv_path}, {off_path}")


if __name__ == "__main__":
    main()
