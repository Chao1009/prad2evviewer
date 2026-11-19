#!/usr/bin/env python3
"""exotic_p4_module_rates.py -- P4 of the exotic-waveform study: rate of each exotic flag per PbWO4 module and run type on the full scan (run_exotic_full.sh).

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p4_module_rates.py [--scan-dir <dir>] [--out <dir>] [--cache <pkl>] [--jobs 8] [--plot-only]

Population: PbWO4 kept channels (scan rows with A_main >= KEEP_ADC; trigger gate and LMS veto as in pileup_rates.scan) of every $OUT_DIR/exotic_full/run_*_s*.npz (default fit,
PILEUP_EXOTIC=1, options off), per run type.  Rate of a flag in a module = flagged kept channels / kept channels of that module.

Run 24937 is excluded (exotic_p1_quantities.EXCLUDE_RUNS; user: not a good run).  Its files stay on disk.
Flags are rebuilt offline from the stored quantities (ch_ex_db, ch_ex_tail, ch_w50, ch_ex_chi2s/chi2k/nm, ch_nclip, ch_ex_dclip, ch_ex_esec, ch_k, ch_t):
  f1 tail-baseline  db > D1 and tail < T1                      (uniform D1 = 8, T1 = -5, in sigma_ref)
  f2 wide           w50 > W2, main pulse not truncated, record not clipped (nclip == 0), and not f4 (uniform W2 = 3)
  f3 smooth-wide    k >= 2, a stretch fit ran (nm > 0) and chi2s <= chi2k + LAM (k - 1) ln nm   (uniform LAM = 2, the P2 rule)
  f4 clipped-multi  nclip >= 1 and (dclip <= 8 ns or esec > 0.5)  (as stored; no per-type version: its thresholds are not distribution cuts)
  f5 saturated single  nclip >= 1 and not f4 (a clipped record without a qualifying secondary pulse; threshold-independent).  f4 + f5 = every clipped record (saturated pulses).
The raw flags overlap (f1/f3/f2 may coincide); each kept channel also gets ONE exclusive category by the priority f4 > f5 > f1 > f2 > f3 (flags_of returns both; maps, top-20 lists
and the exclusive rate table use the exclusive category).  f2 never takes a clipped or truncated row.
Truncation of f2 (user decision 3): the main (largest-A) fitted pulse is truncated when its fit window is cut by the end of the record, i.e. ceil(t0 / DT) - LEAD + LW > NS
(DT = 4 ns, LEAD = 3, LW = 32 samples, NS = 100 samples; t0 > 284 ns).  This is the pipeline's trunc_main criterion, recomputed from the stored t0 and checked against ch_trunc_main.
Truncated pulses are not eligible for f2 under either threshold set; the uniform set otherwise reproduces the stored ch_ex_f1..f4 (checked: f1, f3, f4 exactly, f2 before the exclusion).

Per-run-type thresholds (user decision 1): same-sigma relative to each run type's core.  Core = k = 1 kept PbWO4 channels with a non-truncated main pulse.  For each quantity
q in (db, tail, w50) the core median and sigma (1.4826 MAD) are taken per run type; the P2 uniform threshold is expressed in core sigmas of the reference (the three H2 run types pooled),
z_ref = (thr_uniform - median_ref) / sigma_ref, and the per-type threshold is median_rt + z_ref sigma_rt.  f3: the chi2 penalty scales with the core chi2/dof level of the type,
LAM_rt = 2 median_rt(chi2/dof, core) / median_ref(chi2/dof, core).  User decision 2: the 0.73 GeV db threshold is never raised above the uniform 8 (marginal f1 cases stay flagged).

Every rate and percentage carries its 68 % Wilson interval (PDF tables "p [lo, hi]"; TSV rate_lo68 / rate_hi68 or lo68_percent / hi68_percent; moved shares and f4 / clipped ratio binomial too;
f4_over_clipped rows of _summary.tsv: n_kept = clipped records (denominator), n_flag = f4).  Module maps carry no per-bin annotation (module intervals are in the TSV).
Outputs (--out): exotic_P4_module_rates_20261117.pdf (five main pages: module maps f1..f5, 2x2 run types, per-type thresholds, per-panel log colour scale; summary, threshold, top-20
tables; appendix with the six maps for the uniform thresholds), .tsv (module x run type x flag x threshold set), _thresholds.tsv, _summary.tsv.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import json
import time
import pickle
import argparse
import multiprocessing as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import exotic_p1_quantities as p1                       # noqa: E402
from exotic_p1_quantities import TYPES, RT_LABEL, write_tsv, EXCLUDE_RUNS, EXCL_NOTE     # noqa: E402
from exotic_p3_sample import text_page, pct_ci as ci_txt, ci_cols   # noqa: E402

STEM = "exotic_P4_module_rates_20261117"
FLAGS = (1, 2, 3, 4, 5)
PFLAGS = (1, 2, 3, 4, 5, 6)                           # maps, top lists, module table: f1..f5 and 6 = saturated (clipped record, f4 + f5)
MIN_COUNT = 5                                           # statistics criterion for a per-module rate to be quoted: n_kept >= MIN_TOP and n_flag >= MIN_COUNT (68 % relative half-width <~ 45 %)
RADBINS = (0, 120, 170, 230, 320, 450, 700)             # radial bins [mm] from the HyCal centre for the saturation-vs-radius table
CAT_ORDER = (4, 5, 1, 2, 3)                           # exclusive category priority: saturated first (unambiguous), f1 (a baseline tail is a physical cause), f2 wide, f3 last (residual over-fit class)
ORDER_TXT = "f4 > f5 > f1 > f2 > f3"
FNAME = {1: "f1 tail-baseline", 2: "f2 wide", 3: "f3 smooth-wide", 4: "f4 clipped-multi", 5: "f5 saturated single", 6: "saturated f4+f5 (nclip>=1)"}
UNI = dict(db=8.0, tail=-5.0, w50=3.0, lam=2.0)         # P2 uniform thresholds
NS, DT, LEAD, LW = 100, 4.0, 3, 32                      # record length (samples), sample period (ns), fit window (pileup_rates)
MIN_MAP = 30                                            # modules with fewer kept channels are grey on the maps
MIN_TOP = 100                                           # modules with fewer kept channels are left out of the top-20 lists
HB = dict(db=(-20.0, 20.0, 4000), tail=(-20.0, 20.0, 4000), w50=(0.0, 6.0, 1200), chi2=(0.0, 5.0, 1000))   # core histograms (lo, hi, nbins)
CAND = dict(db=5.0, tail=-3.0, w50=2.0)                 # rows kept for the offline flags: db > 5 and tail < -3, or w50 > 2, or a stretch fit ran with k >= 2, or nclip >= 1, or a stored flag


# ---------------------------------------------------------------------------------------------------------------- pass 1: per file
def load_file(fn):
    """Per-file partial: kept-channel counts per module, core histograms, and the candidate rows (all quantities the flags need)."""
    z = np.load(fn)
    run, rt = int(z["run"]), str(z["run_type"])
    midx = z["ch_midx"].astype(np.int64)
    sel = z["mod_pwo4"][midx]
    nmod = len(z["mod_name"])
    k = z["ch_k"].astype(np.int64)
    a, t = z["ch_a"], z["ch_t"]
    main = np.argmax(np.where(np.isfinite(t), a, -1.0), axis=1)
    t0 = t[np.arange(k.size), main].astype(np.float64)
    with np.errstate(invalid="ignore"):
        trunc = np.ceil(t0 / DT) - LEAD + LW > NS
    trunc &= np.isfinite(t0)
    same = bool((trunc[sel] == z["ch_trunc_main"][sel]).all())
    db, tail, w50 = (z[q][:].astype(np.float64) for q in ("ch_ex_db", "ch_ex_tail", "ch_w50"))
    c2s, c2k, nm = z["ch_ex_chi2s"].astype(np.float64), z["ch_ex_chi2k"].astype(np.float64), z["ch_ex_nm"].astype(np.int64)
    nclip, dclip, esec = z["ch_nclip"].astype(np.int64), z["ch_ex_dclip"].astype(np.float64), z["ch_ex_esec"].astype(np.float64)
    sf = [z[f"ch_ex_f{i}"].astype(bool) for i in (1, 2, 3, 4)]      # stored flags (f5 is offline-only)
    chi2 = z["ch_chi2"].astype(np.float64)
    out = dict(run=run, rt=rt, n=np.bincount(midx[sel], minlength=nmod), trunc_agree=same, nsel=int(sel.sum()), wall=float(z["wall_s"]), nmod=nmod, names=[str(x) for x in z["mod_name"]])
    core = sel & (k == 1) & ~trunc
    out["core_n"] = int(core.sum())
    for q, x in (("db", db), ("tail", tail), ("w50", w50), ("chi2", chi2)):
        lo, hi, nb = HB[q]
        v = x[core]
        v = v[np.isfinite(v)]
        out["h_" + q] = np.histogram(np.clip(v, lo, hi - 1e-9), bins=nb, range=(lo, hi))[0]
    with np.errstate(invalid="ignore"):
        cand = ((db > CAND["db"]) & (tail < CAND["tail"])) | (w50 > CAND["w50"]) | ((nm > 0) & (k >= 2)) | (nclip >= 1) | sf[0] | sf[1] | sf[2] | sf[3]
    cand &= sel
    d = dict(midx=midx, k=k, db=db, tail=tail, w50=w50, c2s=c2s, c2k=c2k, nm=nm, nclip=nclip, dclip=dclip, esec=esec, trunc=trunc,
             s1=sf[0], s2=sf[1], s3=sf[2], s4=sf[3])
    out["cand"] = {q: v[cand].astype(np.float32 if v.dtype == np.float64 else v.dtype) for q, v in d.items()}
    out["cand"]["tmain"] = t0[cand].astype(np.float32)
    # cross-check of the stored flags on all rows is by the candidate set (every stored flag is in it)
    return out


def load_part(fn):
    """The partial of one scan file: from <stem>.p4.pkl when run_exotic_full.sh LEAN=1 reduced it, else from the npz."""
    pk = fn[:-4] + ".p4.pkl"
    return pickle.load(open(pk, "rb")) if os.path.exists(pk) else load_file(fn)


# ---------------------------------------------------------------------------------------------------------------- thresholds
def hist_stats(h, lo, hi):
    """Median and 1.4826 MAD of a histogram (bin centres, linear interpolation inside the median bin)."""
    nb = h.size
    edges = np.linspace(lo, hi, nb + 1)
    c = np.cumsum(h)
    n = c[-1]
    i = np.searchsorted(c, n / 2.0)
    frac = (n / 2.0 - (c[i - 1] if i else 0)) / max(h[i], 1)
    med = edges[i] + frac * (edges[1] - edges[0])
    cen = 0.5 * (edges[:-1] + edges[1:])
    dev = np.abs(cen - med)
    o = np.argsort(dev, kind="stable")
    cc = np.cumsum(h[o])
    j = np.searchsorted(cc, n / 2.0)
    return float(med), float(1.4826 * dev[o][min(j, nb - 1)]), int(n)


def derive_thresholds(core):
    """core[rt][q] = histogram.  Returns per-type thresholds (dict rt -> dict) and the core statistics table rows."""
    stats = {rt: {q: hist_stats(core[rt][q], *HB[q][:2]) for q in HB} for rt in TYPES}
    ref = {q: hist_stats(sum(core[rt][q] for rt in TYPES if "H2" in rt), *HB[q][:2]) for q in HB}
    z = {q: (UNI[q] - ref[q][0]) / ref[q][1] for q in ("db", "tail", "w50")}
    thr = {}
    for rt in TYPES:
        s = stats[rt]
        t = dict(db=s["db"][0] + z["db"] * s["db"][1], tail=s["tail"][0] + z["tail"] * s["tail"][1], w50=s["w50"][0] + z["w50"] * s["w50"][1],
                 lam=UNI["lam"] * s["chi2"][0] / ref["chi2"][0])
        if rt == "0p73GeV_H2":
            t["db"] = min(t["db"], UNI["db"])           # decision 2: the marginal f1 cases at 0.73 GeV stay flagged
        thr[rt] = t
    return thr, stats, ref, z


# ---------------------------------------------------------------------------------------------------------------- flags
def f2_eligible(trunc, nclip):
    """Rows eligible for f2 (shared with P5): main pulse not window-truncated and record not clipped (nclip == 0, the f4 clipping definition).
    f4 is evaluated first, so a clipped record is f4 (when it meets the f4 conditions) or unflagged by f2, never both."""
    return ~np.asarray(trunc, bool) & (np.asarray(nclip) == 0)


def exclusive_category(F):
    """Exclusive category per row from the raw flags F[1..5] (bool arrays): one of 4, 5, 1, 2, 3 by the priority CAT_ORDER, 0 = none.  Shared with P5/P6."""
    cat = np.zeros(np.asarray(F[1]).shape, np.int8)
    for i in reversed(CAT_ORDER):
        cat[np.asarray(F[i], bool)] = i
    return cat


def flags_of(c, th, trunc_excl=True, clip_excl=True):
    """Raw (overlapping) flags f1..f5 of candidate rows c under thresholds th, plus F["cat"], the exclusive category (4, 5, 1, 2, 3 in that priority, 0 = none).
    f4 and f5 are mutually exclusive by definition (f5 = nclip >= 1 and not f4); f2 excludes truncated (trunc_excl) and clipped (clip_excl, nclip >= 1) rows; f1, f2, f3 may overlap."""
    with np.errstate(invalid="ignore", divide="ignore"):
        f1 = (c["db"] > th["db"]) & (c["tail"] < th["tail"])
        f2 = c["w50"] > th["w50"]
        if trunc_excl and clip_excl:
            f2 = f2 & f2_eligible(c["trunc"], c["nclip"])
        elif trunc_excl:                                # diagnostic split for the removal counts
            f2 = f2 & ~c["trunc"]
        f3 = (c["k"] >= 2) & (c["nm"] > 0) & (c["c2s"] <= c["c2k"] + th["lam"] * (c["k"] - 1) * np.log(np.maximum(c["nm"], 1)))
    f4 = c["s4"].astype(bool)
    f5 = (np.asarray(c["nclip"]) >= 1) & ~f4           # saturated single: clipped record that is not f4 (shared with P5/P6)
    f2 = f2 & ~f4
    if clip_excl:
        f2 = f2 & ~f5
    F = {1: f1, 2: f2, 3: f3, 4: f4, 5: f5}
    F["cat"] = exclusive_category(F)
    return F


def assemble(parts):
    """Concatenate the file partials per run type."""
    pop = {rt: dict(n=0, cand=[], core={q: 0 for q in HB}, nfile=0) for rt in TYPES}
    for p in parts:
        d = pop[p["rt"]]
        d["n"] = d["n"] + p["n"]
        d["nfile"] += 1
        d["cand"].append(p["cand"])
        for q in HB:
            d["core"][q] = d["core"][q] + p["h_" + q]
    for rt in TYPES:
        c = pop[rt]["cand"]
        pop[rt]["cand"] = {q: np.concatenate([x[q] for x in c]) for q in c[0]} if c else {}
    return pop


def rates(pop, thr, nmod):
    """R[set][rt] = dict(n [nmod], c{i} [nmod] flagged counts, plus totals); sets 'type' (per-type thresholds), 'uniform'; and f2 truncation numbers."""
    R = {"type": {}, "uniform": {}}
    trunc_info = {}
    for rt in TYPES:
        c = pop[rt]["cand"]
        mi = c["midx"].astype(np.int64)
        for name, th in (("type", thr[rt]), ("uniform", UNI)):
            F = flags_of(c, th)
            r = dict(n=pop[rt]["n"])
            for i in FLAGS:
                r[f"c{i}"] = np.bincount(mi[F[i]], minlength=nmod)
            r["c6"] = r["c4"] + r["c5"]
            anyf = np.any([F[i] for i in FLAGS], axis=0)
            r["cany"] = np.bincount(mi[anyf], minlength=nmod)
            for i in FLAGS:                             # exclusive categories (CAT_ORDER): every kept channel in exactly one of f4, f5, f1, f2, f3, none
                r[f"x{i}"] = np.bincount(mi[F["cat"] == i], minlength=nmod)
            r["x6"] = r["x4"] + r["x5"]
            r["xany"] = np.bincount(mi[F["cat"] > 0], minlength=nmod)
            assert int(r["xany"].sum()) == int(r["cany"].sum())
            r["moved"] = {(i, j): int((F[i] & (F["cat"] == j)).sum()) for i in FLAGS for j in FLAGS if i != j}      # raw flag i whose exclusive category is j (moved out of i into j)
            R[name][rt] = r
        Fu0 = flags_of(c, UNI, trunc_excl=False, clip_excl=False)
        Ft0 = flags_of(c, thr[rt], trunc_excl=False, clip_excl=False)
        Fu1 = flags_of(c, UNI, clip_excl=False)
        Ft1 = flags_of(c, thr[rt], clip_excl=False)
        Fu = flags_of(c, UNI)
        trunc_info[rt] = dict(uniform_before=int(Fu0[2].sum()), uniform_trunc_only=int(Fu1[2].sum()), type_trunc_only=int(Ft1[2].sum()),
                              uni_f1_and_f4=int((Fu[1] & Fu[4]).sum()), uni_f3_and_f4=int((Fu[3] & Fu[4]).sum()), uni_f2_and_f4=int((Fu[2] & Fu[4]).sum()),
                              uni_f5_and_f1=int((Fu[5] & Fu[1]).sum()), uni_f5_and_f2=int((Fu[5] & Fu[2]).sum()), uni_f5_and_f3=int((Fu[5] & Fu[3]).sum()), uni_f5_and_f4=int((Fu[5] & Fu[4]).sum()),
                              uniform_after=int(R["uniform"][rt]["c2"].sum()),
                              type_before=int(Ft0[2].sum()), type_after=int(R["type"][rt]["c2"].sum()),
                              stored_f2=int(c["s2"].sum()), stored_f1=int(c["s1"].sum()), stored_f3=int(c["s3"].sum()), stored_f4=int(c["s4"].sum()),
                              uni_f1=int(R["uniform"][rt]["c1"].sum()), uni_f3=int(R["uniform"][rt]["c3"].sum()), uni_f4=int(R["uniform"][rt]["c4"].sum()), uni_f5=int(R["uniform"][rt]["c5"].sum()),
                              f2_trunc_cands=int(((c["w50"] > UNI["w50"]) & c["trunc"]).sum()))
    return R, trunc_info


# ---------------------------------------------------------------------------------------------------------------- tables
def pctf(c, n):
    return 100.0 * c / max(n, 1)


def wilson(c, n, z=1.0):
    """Wilson score interval (z = 1: 68 %) of a binomial proportion c / n; arrays or scalars, returns (lo, hi) as fractions."""
    c, n = np.asarray(c, float), np.maximum(np.asarray(n, float), 1.0)
    p = c / n
    den = 1.0 + z * z / n
    mid = (p + z * z / (2 * n)) / den
    h = z * np.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    lo, hi = np.where(c <= 0, 0.0, np.clip(mid - h, 0.0, 1.0)), np.where(c >= n, 1.0, np.clip(mid + h, 0.0, 1.0))   # exact 0 at k = 0, exact 1 at k = n
    return (float(lo), float(hi)) if lo.ndim == 0 else (lo, hi)


def pct_ci(c, n):
    return ci_txt(c, n)


def ratio_ci(c, n):
    """Ratio c / n of two counts with c <= n (a share of n: binomial, Wilson 68 %), as 'r [lo, hi]' (fractions, 3 decimals)."""
    lo, hi = wilson(c, n)
    return f"{c / max(n, 1):.3f} [{float(lo):.3f}, {float(hi):.3f}]"


def summary_lines(R, ti):
    """Raw (overlapping) flag counts and rates per run type and threshold set, % [68 % Wilson interval]; two blocks of columns (f1-f3, f4-any) to fit the page."""
    W = 32
    L = []
    for blk in ((1, 2, 3), (4, 5, "any")):
        L += [f"{'run type':<17}{'N kept':>12}{'set':>9}" + "".join(f"{(('f' + str(i)) if i != 'any' else 'any') + ' n':>9}{'% [68 % Wilson]':>{W}}" for i in blk)]
        for rt in TYPES:
            for name in ("type", "uniform"):
                r = R[name][rt]
                n = int(r["n"].sum())
                L.append(f"{RT_LABEL[rt] if name == 'type' else '':<17}{n if name == 'type' else '':>12}{name:>9}"
                         + "".join(f"{int(r['c' + str(i)].sum()):>9d}{pct_ci(int(r['c' + str(i)].sum()), n):>{W}}" for i in blk))
        L.append("")
    return L


def excl_lines(R):
    """Exclusive-category rates (every kept channel in one category, order ORDER_TXT; raw overlapping rates are in the summary table) and the moved matrix, all with 68 % Wilson intervals."""
    W = 32
    L = [f"Exclusive categories, priority {ORDER_TXT} (each kept channel in exactly one category; none = unflagged), % [68 % Wilson interval] of the kept PbWO4 channels:"]
    for blk in ((4, 5, 1), (2, 3, "any")):
        L += [f"{'run type':<17}{'set':>9}{'N kept':>12}" + "".join(f"{(('f' + str(i)) if i != 'any' else 'any') + ' n':>9}{'excl %':>{W}}" for i in blk)]
        for rt in TYPES:
            for name in ("type", "uniform"):
                r = R[name][rt]
                n = int(r["n"].sum())
                L.append(f"{RT_LABEL[rt] if name == 'type' else '':<17}{name:>9}{n:>12d}" + "".join(f"{int(r['x%s' % i].sum()):>9d}{pct_ci(int(r['x%s' % i].sum()), n):>{W}}" for i in blk))
        L.append("")
    L.pop()
    L += ["", "Moved matrix: channels moved out of the raw flag (row) because a higher-priority flag took them (column = the category they were put in); count and share of the raw flag count, % [68 % Wilson interval].",
          "(f2 and f4/f5 never overlap: f2 excludes clipped records; f4 and f5 are exclusive by definition.)",
          f"{'run type':<17}{'set':>9}{'from':>6}" + "".join(f"{'-> f' + str(j):>35}" for j in CAT_ORDER)]
    for rt in TYPES:
        for name in ("uniform", "type"):
            for i in (1, 2, 3):                              # f4 and f5 are never moved
                tot = int(R[name][rt][f"c{i}"].sum())
                L.append(f"{RT_LABEL[rt] if (i == 1 and name == 'uniform') else '':<17}{name if i == 1 else '':>9}{'f' + str(i):>6}"
                         + "".join((f"{R[name][rt]['moved'][(i, j)]:>8d}{pct_ci(R[name][rt]['moved'][(i, j)], tot):>27}" if i != j else f"{'-':>35}") for j in CAT_ORDER))
    return L


def sat_lines(R, ti):
    """Saturated (clipped-record) ratio per run type with 68 % Wilson intervals; threshold-independent (uniform = per-type).  Full data (all scan files): f4 + f5 = nclip >= 1."""
    L = ["Saturated pulses (clipped record, nclip >= 1 = f4 + f5) as a percentage of the kept PbWO4 channels of the full scan; [68 % Wilson interval]:",
         f"{'run type':<17}{'N kept':>12}{'clipped n':>11}{'% [68 % Wilson]':>32}{'f4 n':>8}{'f4 % [68 % Wilson]':>34}{'f5 n':>8}{'f5 % [68 % Wilson]':>34}   {'f4 / clipped [68 % Wilson]':>28}"]
    for rt in TYPES:
        r = R["uniform"][rt]
        n = int(r["n"].sum())
        c4, c5 = int(r["c4"].sum()), int(r["c5"].sum())
        L.append(f"{RT_LABEL[rt]:<17}{n:>12d}{c4 + c5:>11d}{pct_ci(c4 + c5, n):>32}{c4:>8d}{pct_ci(c4, n):>34}{c5:>8d}{pct_ci(c5, n):>34}   {ratio_ci(c4, c4 + c5):>28}")
    L.append("f5 overlaps with the other flags (uniform thresholds; f5&f4 = f5&f2 = 0 by construction; counts):")
    L += [f"  {RT_LABEL[rt]:<17} f5&f1 {ti[rt]['uni_f5_and_f1']}, f5&f3 {ti[rt]['uni_f5_and_f3']}, f5&f2 {ti[rt]['uni_f5_and_f2']}, f5&f4 {ti[rt]['uni_f5_and_f4']}" for rt in TYPES]
    L.append(f"Per-module rates are quoted (saturated map and all top-20 lists) only for modules with n_kept >= {MIN_TOP} and n_flag >= {MIN_COUNT} (68 % relative half-width <~ 45 %); the others are grey (the f1..f5 maps colour modules with n_kept >= {MIN_MAP}).")
    return L


def radius_of(names):
    import pileup_rates as pr
    geo = {m["n"]: m["geo"] for m in json.load(open(pr.HYCAL_MAP))}
    return np.array([float(np.hypot(geo[n]["x"], geo[n]["y"])) if n in geo else np.nan for n in names])


def edge_adc(names):
    """{rt: [nmod]} elastic-edge amplitude beam energy / module gain (ADC) of a run of each type (first sample file; calibration and I_int as pileup_rates.module_gain); nan without a factor."""
    import pileup_rates as pr
    out = {}
    for f in sorted((p1.PYDIR / "output" / "pileup_rates" / "exotic_sample").glob("run_*.npz")):
        z = np.load(f)
        rt = str(z["run_type"])
        if int(z["run"]) in EXCLUDE_RUNS:
            continue
        if rt not in out:
            cfg = pr.run_config(int(z["run"]))
            out[rt] = cfg["beam"] / pr.module_gain(cfg, rt, names)
    return out


def radial_rows(R, names, key="c6"):
    """Rows (rt, r_lo, r_hi, n_mod, n_kept, n_flag): the flag count in PbWO4 radius bins from the HyCal centre."""
    rad = radius_of(names)
    pw = np.array([nm.startswith("W") for nm in names])
    rows = []
    for rt in TYPES:
        r = R["uniform"][rt]
        for lo, hi in zip(RADBINS[:-1], RADBINS[1:]):
            m = pw & (rad >= lo) & (rad < hi)
            rows.append((rt, lo, hi, int(m.sum()), int(r["n"][m].sum()), int(r[key][m].sum())))
    return rows


def radial_lines(R, names):
    L = ["Saturated fraction vs distance from the HyCal centre (PbWO4 modules), % [68 % Wilson interval] (count of clipped records / kept channels):",
         f"{'bin':<12}" + "".join(f"{RT_LABEL[rt]:>50}" for rt in TYPES)]
    rr = radial_rows(R, names)
    for lo, hi in zip(RADBINS[:-1], RADBINS[1:]):
        row = [x for x in rr if x[1] == lo and x[2] == hi]
        L.append(f"{f'{lo}-{hi} mm':<12}" + "".join(f"{(f'{pct_ci(c, n)} ({c}/{n})' if n else '-'):>50}" for rt in TYPES for _, _, _, _, n, c in [x for x in row if x[0] == rt]))
    rad = radius_of(names)
    pw = np.array([nm.startswith("W") for nm in names])
    ea = edge_adc(names)
    L.append("elastic-edge amplitude (beam energy / module gain, median over the modules of the bin, ADC; the 12-bit scale clips at ~4095 - pedestal):")
    for rt in TYPES:
        L.append(f"{RT_LABEL[rt]:<17}" + "".join(f"{np.nanmedian(ea[rt][pw & (rad >= lo) & (rad < hi)]):>14.0f}" for lo, hi in zip(RADBINS[:-1], RADBINS[1:])) + f"   (all PbWO4: {np.nanmedian(ea[rt][pw]):.0f})")
    return L


def thr_lines(thr, stats, ref, z, ti):
    L = ["core = k = 1 kept PbWO4 channels, main pulse not truncated; median / sigma (1.4826 MAD) of each quantity; per-type threshold = median + z_ref sigma (z_ref from the H2 pooled core, uniform value)",
         f"{'':<17}{'N core':>12}  {'db med':>8}{'sig':>7}  {'tail med':>9}{'sig':>7}  {'w50 med':>8}{'sig':>7}  {'chi2/dof med':>13}"]
    for rt in TYPES:
        s = stats[rt]
        L.append(f"{RT_LABEL[rt]:<17}{s['db'][2]:>12d}  {s['db'][0]:>8.3f}{s['db'][1]:>7.3f}  {s['tail'][0]:>9.3f}{s['tail'][1]:>7.3f}  {s['w50'][0]:>8.3f}{s['w50'][1]:>7.3f}  {s['chi2'][0]:>13.3f}")
    L.append(f"{'H2 pooled (ref)':<17}{ref['db'][2]:>12d}  {ref['db'][0]:>8.3f}{ref['db'][1]:>7.3f}  {ref['tail'][0]:>9.3f}{ref['tail'][1]:>7.3f}  {ref['w50'][0]:>8.3f}{ref['w50'][1]:>7.3f}  {ref['chi2'][0]:>13.3f}")
    L += ["", f"z_ref (sigma of the H2-pooled core between its median and the uniform threshold): db {z['db']:.2f}, tail {z['tail']:.2f}, w50 {z['w50']:.2f}", "",
          f"{'thresholds':<17}{'f1 db >':>10}{'f1 tail <':>11}{'f2 w50 >':>10}{'f3 LAM':>9}{'f4':>20}"]
    L.append(f"{'uniform (P2)':<17}{UNI['db']:>10.3f}{UNI['tail']:>11.3f}{UNI['w50']:>10.3f}{UNI['lam']:>9.3f}{'nclip>=1,8ns|0.5':>20}")
    for rt in TYPES:
        t = thr[rt]
        L.append(f"{RT_LABEL[rt]:<17}{t['db']:>10.3f}{t['tail']:>11.3f}{t['w50']:>10.3f}{t['lam']:>9.3f}{'(unchanged)':>20}")
    L += ["", "f3: flag if k >= 2 and chi2_stretch <= chi2_k + LAM (k - 1) ln(nm).  f2 excludes pulses whose fit window is cut by the end of the record (t0 > 284 ns) and clipped records (nclip >= 1).",
          "0.73 GeV H2 db threshold capped at the uniform 8 (decision 2).", "", "f2 exclusions (truncated main pulse, then clipped record nclip >= 1; f4 evaluated first), flags removed:",
          f"{'run type':<17}{'uniform before':>16}{'trunc':>9}{'clip':>9}{'after':>9}{'per-type before':>17}{'trunc':>9}{'clip':>9}{'after':>9}{'f1&f4':>8}{'f3&f4':>8}"]
    for rt in TYPES:
        a = ti[rt]
        L.append(f"{RT_LABEL[rt]:<17}{a['uniform_before']:>16d}{a['uniform_before'] - a['uniform_trunc_only']:>9d}{a['uniform_trunc_only'] - a['uniform_after']:>9d}{a['uniform_after']:>9d}{a['type_before']:>17d}{a['type_before'] - a['type_trunc_only']:>9d}{a['type_trunc_only'] - a['type_after']:>9d}{a['type_after']:>9d}{a['uni_f1_and_f4']:>8d}{a['uni_f3_and_f4']:>8d}")
    return L


def top_lists(R, names, name, top=20):
    """{(rt, i): [(module, n, c, rate %)]} top modules by rate with n >= MIN_TOP."""
    out = {}
    for rt in TYPES:
        n = R[name][rt]["n"]
        for i in PFLAGS:
            c = R[name][rt][f"x{i}"]
            rate = np.where((n >= MIN_TOP) & (c >= MIN_COUNT), c / np.maximum(n, 1), -1.0)
            o = np.argsort(-rate, kind="stable")[:top]
            lo, hi = wilson(c, n)
            out[(rt, i)] = [(names[j], int(n[j]), int(c[j]), 100.0 * c[j] / n[j], 100.0 * lo[j], 100.0 * hi[j]) for j in o if rate[j] > 0]
    return out


def top_pages(pdf, plt, R, names, name, label):
    T = top_lists(R, names, name)
    for i in PFLAGS:
        L = [f"top-20 modules by {FNAME[i]} EXCLUSIVE rate (order {ORDER_TXT}; {label}; modules with n_kept >= {MIN_TOP} and n_flag >= {MIN_COUNT}); module  n_kept  n_flag  rate % [68 % Wilson interval]", ""]
        rows = []
        for r in range(20):
            row = ""
            for rt in TYPES:
                lst = T[(rt, i)]
                row += (f"{lst[r][0]:<6}{lst[r][1]:>8d}{lst[r][2]:>6d}  {pct_ci(lst[r][2], lst[r][1])}".ljust(56) if r < len(lst) else " " * 56)
            rows.append(row)
        L.append("".join(f"{RT_LABEL[rt]:<56}" for rt in TYPES))
        L += rows
        text_page(pdf, plt, f"Top-20 modules: {FNAME[i]} (exclusive category)", L, size=8)


# ---------------------------------------------------------------------------------------------------------------- maps
def map_page(pdf, plt, R, names, name, i, label):
    import pileup_rates as pr
    from matplotlib.patches import Rectangle
    from matplotlib.collections import PatchCollection
    from matplotlib.colors import LogNorm
    geo = {m["n"]: m["geo"] for m in json.load(open(pr.HYCAL_MAP))}
    pw = [j for j, nme in enumerate(names) if nme.startswith("W")]
    rect = {j: Rectangle((geo[names[j]]["x"] - geo[names[j]]["sx"] / 2, geo[names[j]]["y"] - geo[names[j]]["sy"] / 2), geo[names[j]]["sx"], geo[names[j]]["sy"]) for j in pw}
    fig, axs = plt.subplots(2, 2, figsize=(17, 16))
    fig.suptitle(f"{FNAME[i]} (exclusive category, order {ORDER_TXT}): rate per PbWO4 module ({label});\nlog colour scale, own range per panel; the per-module 68 % Wilson intervals are in the TSV", fontsize=14, fontweight="bold")
    for ax, rt in zip(axs.flat, TYPES):
        n, c = R[name][rt]["n"], R[name][rt][f"x{i}"]
        rate = 100.0 * c / np.maximum(n, 1)
        if i == 6:                                      # saturation is rare: quote a module only with n_kept >= MIN_TOP and n_flag >= MIN_COUNT (zero counts need n_kept >= MIN_TOP)
            pos = [j for j in pw if n[j] >= MIN_TOP and c[j] >= MIN_COUNT]
            zero = [j for j in pw if n[j] >= MIN_TOP and c[j] == 0]
            few = [j for j in pw if j not in set(pos) and j not in set(zero)]
        else:
            ok = [j for j in pw if n[j] >= MIN_MAP]
            pos = [j for j in ok if c[j] > 0]
            zero = [j for j in ok if c[j] == 0]
            few = [j for j in pw if n[j] < MIN_MAP]
        ax.add_collection(PatchCollection([rect[j] for j in few], facecolor="0.8", edgecolor="0.6", lw=0.2))
        ax.add_collection(PatchCollection([rect[j] for j in zero], facecolor="white", edgecolor="0.7", lw=0.2))
        if pos:
            v = rate[pos]
            pc = PatchCollection([rect[j] for j in pos], cmap="viridis", edgecolor="0.3", lw=0.2, norm=LogNorm(vmin=max(v.min(), v.max() * 1e-3), vmax=v.max()))
            pc.set_array(v)
            ax.add_collection(pc)
            fig.colorbar(pc, ax=ax, fraction=0.045, label="rate (%)")
        tot = R[name][rt]
        ax.set_xlim(-600, 600), ax.set_ylim(-600, 600)
        ax.set_aspect("equal")
        ax.set_title(f"{RT_LABEL[rt]}: {int(c.sum())} / {int(n.sum())} = {pct_ci(int(c.sum()), int(n.sum()))} % (68 % Wilson); white: none, grey: " + (f"n < {MIN_TOP} or 0 < flag < {MIN_COUNT}" if i == 6 else f"< {MIN_MAP} ch.") + ")", fontsize=8.5)
    fig.tight_layout(rect=(0, 0, 1, 0.96))
    pdf.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------- output
def structure_line(R, names, name, rt, i):
    """Mean rate inside the innermost two module rings / outer edge / elsewhere is not defined here (no ring map); report the share of the flag in the top 10 modules."""
    c = R[name][rt][f"c{i}"]
    o = np.sort(c)[::-1]
    return f"{o[:10].sum() / max(o.sum(), 1):.2f}"


def tsv_modules(R, names, thr):
    rows = []
    for name in ("type", "uniform"):
        for rt in TYPES:
            r = R[name][rt]
            for j, nme in enumerate(names):
                if not nme.startswith("W") or not r["n"][j]:
                    continue
                for basis, pre in (("raw", "c"), ("excl", "x")):
                    for i in PFLAGS:
                        n_, c_ = int(r["n"][j]), int(r[f"{pre}{i}"][j])
                        lo, hi = wilson(c_, n_)
                        rows.append(dict(threshold_set=name, run_type=rt, module=nme, midx=j, basis=basis, flag=f"f{i}" if i < 6 else "sat", n_kept=n_, n_flag=c_, rate=c_ / n_, rate_lo68=float(lo), rate_hi68=float(hi),
                                         stat_ok=int(n_ >= MIN_TOP and c_ >= MIN_COUNT)))
    return rows


def report(data, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    R, names, thr, stats, ref, z, ti = (data[k] for k in ("R", "names", "thr", "stats", "ref", "z", "ti"))
    stem = Path(out) / STEM
    sl = summary_lines(R, ti)
    with PdfPages(str(stem) + ".pdf") as pdf:
        for i in PFLAGS:
            map_page(pdf, plt, R, names, "type", i, "per-run-type thresholds" if i != 6 else "threshold-independent; full scan")
        text_page(pdf, plt, "P4: exotic-flag rates per run type, full scan (kept PbWO4 channels; per-type vs uniform P2 thresholds)",
                  [f"population: {data['nfiles']} scan files of 12 runs, all EVIO splits; default fit, PILEUP_EXOTIC=1 only.  f2 excludes truncated main pulses and clipped records in both sets; f4 and f5 are the same in both.",
                   f"*** {EXCL_NOTE} ***", f"Categories are mutually exclusive in the maps and top-20 lists, assigned in the order {ORDER_TXT} (f4 saturated first, f5 saturated single, f1 baseline tail, f2 wide, f3 smooth-wide = residual over-fit class).", "Raw rows below are the overlapping flags; every rate is % [68 % Wilson interval].", ""]
                  + sl, size=8)
        text_page(pdf, plt, "P4: saturated pulses (clipped records) per run type and vs radius; rates with 68 % Wilson intervals", sat_lines(R, ti) + [""] + radial_lines(R, names), size=7.3)
        text_page(pdf, plt, f"P4: exclusive categories ({ORDER_TXT}) and the moved matrix; {EXCL_NOTE}", excl_lines(R), size=7.2)
        text_page(pdf, plt, "Thresholds: per-run-type rule, core statistics, f2 truncation exclusion", thr_lines(thr, stats, ref, z, ti), size=9)
        top_pages(pdf, plt, R, names, "type", "per-type thresholds")
        text_page(pdf, plt, "Appendix: the same maps with the uniform P2 thresholds", ["Six maps follow (f1..f5 and the saturated f4+f5); per-module rates of both threshold sets are in exotic_P4_module_rates_20261117.tsv (threshold_set = type | uniform)."])
        for i in PFLAGS:
            map_page(pdf, plt, R, names, "uniform", i, "uniform P2 thresholds" if i != 6 else "threshold-independent; full scan")
    write_tsv(str(stem) + ".tsv", tsv_modules(R, names, thr), ["threshold_set", "run_type", "module", "midx", "basis", "flag", "n_kept", "n_flag", "rate", "rate_lo68", "rate_hi68", "stat_ok"])
    rows = []
    for rt in TYPES:
        for q, key in (("f1 db >", "db"), ("f1 tail <", "tail"), ("f2 w50 >", "w50"), ("f3 LAM", "lam")):
            rows.append(dict(run_type=rt, quantity=q, uniform=UNI[key], per_type=thr[rt][key]))
        for q in ("db", "tail", "w50", "chi2"):
            rows.append(dict(run_type=rt, quantity=f"core {q} median", uniform="", per_type=stats[rt][q][0]))
            rows.append(dict(run_type=rt, quantity=f"core {q} sigma", uniform="", per_type=stats[rt][q][1]))
    write_tsv(str(stem) + "_thresholds.tsv", rows, ["run_type", "quantity", "uniform", "per_type"])
    rows = []
    for rt in TYPES:
        for name in ("type", "uniform"):
            r = R[name][rt]
            n = int(r["n"].sum())
            for key in ("1", "2", "3", "4", "5", "clipped", "any"):
                c = int(r["c4"].sum() + r["c5"].sum()) if key == "clipped" else int(r[f"c{key}"].sum())
                lo, hi = wilson(c, n)
                rows.append(dict(run_type=rt, threshold_set=name, flag=f"f{key}" if key.isdigit() else key, n_kept=n, n_flag=c, rate_percent=pctf(c, n), lo68_percent=100 * float(lo), hi68_percent=100 * float(hi)))
            for key in ("1", "2", "3", "4", "5", "any"):
                c = int(r[f"x{key}"].sum())
                lo, hi = wilson(c, n)
                rows.append(dict(run_type=rt, threshold_set=name, flag=f"excl_f{key}" if key.isdigit() else "excl_any", n_kept=n, n_flag=c, rate_percent=pctf(c, n), lo68_percent=100 * float(lo), hi68_percent=100 * float(hi)))
            c4, c5 = int(r["c4"].sum()), int(r["c5"].sum())
            rlo, rhi = wilson(c4, c4 + c5)
            rows.append(dict(run_type=rt, threshold_set=name, flag="f4_over_clipped", n_kept=c4 + c5, n_flag=c4, rate_percent=pctf(c4, c4 + c5), lo68_percent=100 * float(rlo), hi68_percent=100 * float(rhi)))
            for (i, j), v in r["moved"].items():
                lo, hi = wilson(v, n)
                tot = int(r[f"c{i}"].sum())
                olo, ohi = wilson(v, tot)
                rows.append(dict(run_type=rt, threshold_set=name, flag=f"moved_f{i}_to_f{j}", n_kept=n, n_flag=v, rate_percent=pctf(v, n), lo68_percent=100 * float(lo), hi68_percent=100 * float(hi),
                                 of_raw_percent=pctf(v, tot), of_raw_lo68_percent=100 * float(olo), of_raw_hi68_percent=100 * float(ohi)))
        a = ti[rt]
        for key in ("uniform_before", "uniform_trunc_only", "uniform_after", "type_trunc_only", "uni_f1_and_f4", "uni_f3_and_f4", "uni_f2_and_f4", "type_before", "type_after", "stored_f1", "stored_f2", "stored_f3", "stored_f4", "uni_f1", "uni_f3", "uni_f4", "uni_f5", "uni_f5_and_f1", "uni_f5_and_f2", "uni_f5_and_f3", "uni_f5_and_f4"):
            rows.append(dict(run_type=rt, threshold_set="check", flag=key, n_kept="", n_flag=a[key], rate_percent="", lo68_percent="", hi68_percent=""))
    for rt, lo_r, hi_r, nm_, n_, c_ in radial_rows(R, names):
        lo, hi = wilson(c_, n_)
        rows.append(dict(run_type=rt, threshold_set=f"radius_{lo_r}-{hi_r}mm", flag="sat", n_kept=n_, n_flag=c_, rate_percent=pctf(c_, n_), lo68_percent=100 * float(lo), hi68_percent=100 * float(hi)))
    write_tsv(str(stem) + "_summary.tsv", rows, ["run_type", "threshold_set", "flag", "n_kept", "n_flag", "rate_percent", "lo68_percent", "hi68_percent", "of_raw_percent", "of_raw_lo68_percent", "of_raw_hi68_percent"])
    print("\n".join(sl + [""] + excl_lines(R) + [""] + sat_lines(R, ti) + [""] + radial_lines(R, names) + [""] + thr_lines(thr, stats, ref, z, ti)))
    for name in ("type", "uniform"):
        T = top_lists(R, names, name, 5)
        for rt in TYPES:
            for i in PFLAGS:
                print(f"[top5 {name}] {rt} f{i}: " + "; ".join(f"{a} {pct_ci(c, b)} % ({c}/{b})" for a, b, c, d, lo, hi in T[(rt, i)]) + f"  [top-10 share {structure_line(R, names, name, rt, i)}]")


def reduce_file(fn):
    """--reduce <npz>: write <npz stem>.p4.pkl (the per-file partial of load_file: counts per module, core histograms, candidate rows); used by run_exotic_full.sh LEAN=1 to keep the disk use small."""
    out = str(fn)[:-4] + ".p4.pkl"
    pickle.dump(load_file(str(fn)), open(out + ".tmp", "wb"))
    os.replace(out + ".tmp", out)
    chk = pickle.load(open(out, "rb"))
    print(f"[reduce] {out}: {chk['nsel']} kept PbWO4 channels, {chk['cand']['midx'].size} candidate rows")


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "--reduce":
        return reduce_file(sys.argv[2])
    ap = argparse.ArgumentParser()
    ap.add_argument("--scan-dir", default=str(p1.PYDIR / "output" / "pileup_rates" / "exotic_full"))
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--cache", default=None)
    ap.add_argument("--jobs", type=int, default=8)
    ap.add_argument("--plot-only", action="store_true")
    a = ap.parse_args()
    cache = a.cache or str(Path(a.scan_dir) / "p4_cache.pkl")
    if a.plot_only:
        data = pickle.load(open(cache, "rb"))
    else:
        t0 = time.time()
        d = Path(a.scan_dir)
        files = sorted(set(str(f) for f in d.glob("run_*_s*.npz") if ".part" not in f.name) | set(str(f)[:-7] + ".npz" for f in d.glob("run_*_s*.p4.pkl")))
        nall = len(files)
        files = [f for f in files if int(Path(f).name.split("_")[1]) not in EXCLUDE_RUNS]
        print(f"[load] {EXCL_NOTE}: {nall - len(files)} of {nall} scan files dropped")
        if not files:
            sys.exit(f"[ERROR] no run_*_s*.npz in {a.scan_dir}")
        with mp.get_context("fork").Pool(a.jobs) as pool:
            parts = pool.map(load_part, files, chunksize=1)
        bad = [p["run"] for p in parts if not p["trunc_agree"]]
        print(f"[load] {len(files)} files, {time.time() - t0:.0f} s; truncation criterion differs from ch_trunc_main in runs: {bad or 'none'}")
        pop = assemble(parts)
        nmod = parts[0]["nmod"]
        import pileup_rates as pr
        names = parts[0]["names"]
        core = {rt: pop[rt]["core"] for rt in TYPES}
        thr, stats, ref, z = derive_thresholds(core)
        R, ti = rates(pop, thr, nmod)
        data = dict(R=R, names=names, thr=thr, stats=stats, ref=ref, z=z, ti=ti, nfiles=len(files), scan_wall=sum(p["wall"] for p in parts), trunc_bad=bad)
        pickle.dump(data, open(cache, "wb"))
        print(f"[compute] {time.time() - t0:.0f} s; sum of scan wall times {data['scan_wall']:.0f} s")
    report(data, a.out)


if __name__ == "__main__":
    main()
