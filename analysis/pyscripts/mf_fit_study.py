#!/usr/bin/env python3
"""mf_fit_study.py -- offline analysis of the free-tau pulse-fit study (plan_20261023).

Reads the benchmark outputs of the three arms (fixed = current mf template fit,
bounded / prior = per-pulse free two-tau refit, mf_fit_pulse.py) written under a root:

  <root>/dump/{bounded,prior}/<run>_<split>/pulses.npz      (--pulse-dump)
  <root>/res/{fixed,bounded,prior}/<run>_<split>/result.tsv

Subcommands (all take --root; MPLBACKEND=Agg):
  diag  plots in <root>/plots/, text <root>/diag_summary.txt, tables <root>/tables/*.tsv
  res   <root>/resolution_table.{md,tsv}: sigma_E/E per arm, paired delta, yields
  cost  <root>/cost.txt: microbenchmark rawsum / mf / free fit (needs setup_env.sh)
  intraw  integrated template vs rawsum integral plots/table in <root>/../int_vs_rawsum
  --arm sampled (diag, intraw): adds the sampled-template (spline bank) arm from the dump
  <root>/dump/sampled/<run>_<split>/pulses.npz (or <root>/../dump/sampled/..), a plain mf
  dump without --mf-shape: its mf rows (E_units, amp, chi2dof) are paired with the arms
  above by (split, event, midx), single rows on both sides, |dt| <= SAMPLED_DT_MAX.  Outputs
  get a "sampled" label (tables/chi2_sampled_<run>.tsv, plots/chi2_vs_amp_sampled_<run>.png,
  tables/eraw_sampled_<run>.tsv, plots/eraw_ratio_sampled_<run>.png,
  diag_summary_sampled.txt; int_vs_rawsum/*_sampled*); existing outputs are not rewritten.
  The existing arms are those of --root: use <study>/uniform (headline, uniform weights).
  --arm NAME (any name) with --arm-dir DIR reads the arm's dumps at DIR/<run>_<split>/pulses.npz
  (default <root>/dump/NAME, or <root>/../dump/NAME); outputs and labels use NAME.
  diag --out DIR writes plots/, tables/ and the summary under DIR (default: <root>).
  intraw --bins E0,E1,... (ADC edges): per run median Q_template/Q_raw and SE per bin, split by the
  fallback class of the arm's bank (--bank PATTERN with {runtype}), flatness verdict;
  <out>/intraw_bins_<arm>.{tsv,md} (arm = the --arm name, else fixed).

diag uses the mf rows with ft_status >= 0 (refit pulses); amplitude axis = fx_amp
(mf amplitude, ADC).  Status-0 means converged fits; failed rows keep the fixed value.
"""
from __future__ import annotations

import argparse
import os
import re
import sys
import time
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

DEFAULT_ROOT = HERE / "output" / "mf_fit_study_20261023"
RUNS = {"025320": ["00050", "00150", "00250"], "024503": ["00010", "00040", "00070"]}
DUMP_ARMS = ("bounded", "prior")
RES_ARMS = ("fixed", "bounded", "prior")
COL = {"fixed": "k", "bounded": "#1f77b4", "prior": "#d95f02", "sampled": "#2ca02c"}
AMP_EDGES = np.array([5, 10, 20, 50, 100, 200, 500, 1000, 2000, 4000, np.inf])
A_HI = 300.0              # "high amplitude" cut (ADC) for per-module and ratio studies
MIN_BIN = 20              # min pulses per amplitude bin for a statistic
MIN_MOD = 10              # min status-0 pulses at A > A_HI for a module median
SAMPLED_DT_MAX = 8.0      # ns, max |t_mf(sampled) - t_mf(existing arms)| for a pulse pair
NAN = float("nan")

MF_COLS = ("event", "midx", "time_ns", "amp", "chi2dof", "E_units", "fx_amp",
           "fx_chi2dof", "fx_E_units", "ft_status", "ft_tau_r", "ft_tau_f",
           "ft_sig_tau_r", "ft_sig_tau_f", "ft_sigw_tau_r", "ft_sigw_tau_f",
           "ft_bound")


# ------------------------------------------------------------------ helpers --
def bin_centers():
    return np.sqrt(AMP_EDGES[:-1] * np.where(np.isfinite(AMP_EDGES[1:]), AMP_EDGES[1:], 8000.0))


def bin_labels():
    e = AMP_EDGES
    return [f"{e[i]:g}-{e[i + 1]:g}" if np.isfinite(e[i + 1]) else f"{e[i]:g}+"
            for i in range(len(e) - 1)]


def bin_index(a):
    """Amplitude bin of each value, -1 outside (A < 5 or nan)."""
    a = np.asarray(a, float)
    ok = np.isfinite(a)
    return np.where(ok, np.digitize(np.where(ok, a, 0.0), AMP_EDGES[:-1]) - 1, -1)


def rsig(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    return 1.4826 * float(np.median(np.abs(x - np.median(x)))) if x.size else NAN


def bin_stat(bi, x, fn, minn=MIN_BIN):
    """fn(x[bi == k]) per amplitude bin (nan below minn finite values)."""
    out = np.full(len(AMP_EDGES) - 1, NAN)
    for k in range(len(out)):
        v = x[bi == k]
        v = v[np.isfinite(v)]
        if v.size >= minn:
            out[k] = fn(v)
    return out


def counts(bi, mask=None):
    m = np.ones(len(bi), bool) if mask is None else mask
    return np.bincount(bi[m & (bi >= 0)], minlength=len(AMP_EDGES) - 1)


def pct(p):
    return lambda v: float(np.percentile(v, p))


def fmt(x, nd=4):
    return "NA" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"


def write_tsv(path, header, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "w") as f:
        f.write("\t".join(header) + "\n")
        for r in rows:
            f.write("\t".join(fmt(v, 5) if isinstance(v, (float, np.floating)) else str(v)
                              for v in r) + "\n")


def savefig(fig, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=110)
    plt.close(fig)


def module_median(mi, x, nmod, minn):
    """Per-module median of x (nan when fewer than minn values)."""
    med = np.full(nmod, NAN)
    if mi.size == 0:
        return med
    o = np.argsort(mi, kind="stable")
    u, st = np.unique(mi[o], return_index=True)
    for m, part in zip(u, np.split(x[o], st[1:])):
        if part.size >= minn:
            med[m] = np.median(part)
    return med


# ------------------------------------------------------------------- loading --
class Run:
    """Pooled mf rows of one run, per arm; row i is the same pulse in every arm."""

    def __init__(self, run):
        self.run = run
        self.arms = {}
        self.splits = []
        self.tau0 = {}
        self.sig0 = {}
        self.names = None
        self.kappa = NAN


def _load_split(path):
    z = np.load(path)
    key = {k: z[k] for k in ("event", "midx", "method", "time_ns")}
    mf = z["method"] == 1
    d = {k: z[k][mf] for k in MF_COLS}
    rw = z["method"] == 0
    meta = dict(tau_r=z["tau_r"], tau_f=z["tau_f"], tau_sig_r=z["tau_sig_r"],
                tau_sig_f=z["tau_sig_f"], names=z["module_names"],
                kappa=float(z["kappa"]), mode=str(z["mf_shape"]),
                raw={k: z[k][rw] for k in ("event", "midx", "npk", "E_units")})
    return key, d, meta


def _key(sp, ev, mi):
    return (np.asarray(sp, np.int64) << 40) | (np.asarray(ev, np.int64) << 12) | np.asarray(mi, np.int64)


def pair_rows(k_ref, t_ref, k_s, t_s):
    """Index into the sampled rows for each reference row (-1: none).  A pair needs the key
    to be unique on both sides (single mf row of that channel/event) and |dt| <= SAMPLED_DT_MAX."""
    out = np.full(len(k_ref), -1, np.int64)
    if len(k_s) == 0 or len(k_ref) == 0:
        return out
    uk, ik, ck = np.unique(k_s, return_index=True, return_counts=True)
    _, inv_r, cnt_r = np.unique(k_ref, return_inverse=True, return_counts=True)
    j = np.clip(np.searchsorted(uk, k_ref), 0, uk.size - 1)
    hit = (uk[j] == k_ref) & (ck[j] == 1) & (cnt_r[inv_r] == 1)
    s = ik[j]
    dt = np.abs(t_s[s] - t_ref)
    hit &= np.isfinite(dt) & (dt <= SAMPLED_DT_MAX)
    out[hit] = s[hit]
    return out


def sampled_path(root, run, s, arm="sampled", arm_dir=None):
    bases = [Path(arm_dir)] if arm_dir else [root / "dump" / arm, root.parent / "dump" / arm]
    for base in bases:
        p = base / f"{run}_{s}" / "pulses.npz"
        if p.is_file():
            return p
    return None


def _load_sampled_split(path):
    z = np.load(path)
    mf, rw = z["method"] == 1, z["method"] == 0
    d = {k: z[k][mf] for k in ("event", "midx", "time_ns", "amp", "chi2dof", "E_units")}
    raw = {k: z[k][rw] for k in ("event", "midx", "npk", "E_units")}
    return d, raw, float(z["kappa"]), z["module_names"]


def load_sampled(root, run, R, arm="sampled", arm_dir=None):
    """Sampled-template mf rows paired to the rows of R (arrays of length R.n, nan where
    unpaired) plus ok mask.  The rawsum match uses the sampled dump's own rawsum rows."""
    parts, raws, kap = [], [], NAN
    empty = {k: np.zeros(0, np.int64) for k in ("event", "midx", "npk")}
    empty["E_units"] = np.zeros(0)
    for si, s in enumerate(R.splits):
        p = sampled_path(root, run, s, arm, arm_dir)
        if p is None:
            raws.append(empty)
            continue
        d, raw, kp, names = _load_sampled_split(p)
        if not np.array_equal(names, R.names):
            sys.exit(f"FATAL: {p}: module list differs from the existing arms")
        d["split"] = np.full(len(d["event"]), si, np.int8)
        parts.append(d)
        raws.append(raw)
        kap = kp
    if not parts:
        print(f"[diag] run {run}: no {arm} dumps for splits {R.splits}; skipped")
        return None
    S = {k: np.concatenate([d[k] for d in parts]) for k in parts[0]}
    S["eraw"] = match_raw(R, S, raws)
    base = R.arms[next(iter(R.arms))]
    j = pair_rows(_key(base["split"], base["event"], base["midx"]), base["time_ns"],
                  _key(S["split"], S["event"], S["midx"]), S["time_ns"])
    ok = j >= 0
    A = {"ok": ok, "kappa": kap, "splits": sorted({int(x) for x in S["split"]})}
    for k in ("E_units", "amp", "chi2dof", "eraw"):
        A[k] = np.where(ok, S[k][np.where(ok, j, 0)], NAN)
    print(f"[diag] run {run}: {arm} dumps in {len(A['splits'])}/{len(R.splits)} splits, "
          f"{len(S['event'])} mf rows, {int(ok.sum())} paired to the {R.n} rows "
          f"({int((ok & R.refit).sum())} of {int(R.refit.sum())} refit rows)")
    return A


def match_raw(R, b, raws):
    """E_raw (rawsum integral units) of the single-WaveAnalyzer-peak (npk == 1) rawsum row of
    the same (split, event, midx) for each mf row; nan when there is none."""
    rk = np.concatenate([_key(np.full(len(r["event"]), si), r["event"], r["midx"])
                         for si, r in enumerate(raws)])
    re_ = np.concatenate([r["E_units"] for r in raws])
    one = np.concatenate([r["npk"] for r in raws]) == 1
    rk, re_ = rk[one], re_[one]
    o = np.argsort(rk, kind="stable")
    rk, re_ = rk[o], re_[o]
    if rk.size > 1 and (np.diff(rk) == 0).any():
        print(f"[diag] WARNING run {R.run}: duplicate npk==1 rawsum rows for one (event, midx)")
    k = _key(b["split"], b["event"], b["midx"])
    i = np.clip(np.searchsorted(rk, k), 0, max(rk.size - 1, 0))
    hit = (rk[i] == k) if rk.size else np.zeros(len(k), bool)
    return np.where(hit, re_[i] if rk.size else NAN, NAN)


def load_run(root, run, splits):
    """Load the available arms of a run; verify identical pulse tables (fail loudly)."""
    avail = {a: [s for s in splits if (root / "dump" / a / f"{run}_{s}" / "pulses.npz").is_file()]
             for a in DUMP_ARMS}
    present = [a for a in DUMP_ARMS if avail[a]]
    if not present:
        print(f"[diag] run {run}: no dumps under {root / 'dump'}; skipped")
        return None
    use = [s for s in splits if all(s in avail[a] for a in present)]
    if not use:
        print(f"[diag] run {run}: arms {present} share no split; skipped")
        return None
    miss = [a for a in DUMP_ARMS if a not in present]
    if miss:
        print(f"[diag] run {run}: arm(s) {miss} missing; pair checks skipped")
    R = Run(run)
    R.splits = use
    parts = {a: [] for a in present}
    raws = []
    for si, s in enumerate(use):
        keys, metas = {}, {}
        for a in present:
            k, d, m = _load_split(root / "dump" / a / f"{run}_{s}" / "pulses.npz")
            if m["mode"] != a:
                sys.exit(f"FATAL: {run}_{s} dump in dump/{a} has mf_shape={m['mode']}")
            d["split"] = np.full(len(d["event"]), si, np.int8)
            keys[a], metas[a] = k, m
            parts[a].append(d)
        if len(present) == 2:
            a0, a1 = present
            for c in ("event", "midx", "method", "time_ns"):
                if not np.array_equal(keys[a0][c], keys[a1][c], equal_nan=True):
                    sys.exit(f"FATAL: {run}_{s}: pulse tables of {a0} and {a1} differ in "
                             f"'{c}' (rows/events/pulse identification not identical)")
            if not np.array_equal(parts[a0][-1]["ft_status"] >= 0,
                                  parts[a1][-1]["ft_status"] >= 0):
                sys.exit(f"FATAL: {run}_{s}: refit-row selection differs between arms")
            if not np.allclose(parts[a0][-1]["fx_amp"], parts[a1][-1]["fx_amp"],
                               rtol=1e-9, equal_nan=True):
                print(f"[diag] WARNING {run}_{s}: fx_amp differs between arms")
        m0 = metas[present[0]]
        raws.append(m0["raw"])
        if R.names is None:
            R.names, R.kappa = m0["names"], m0["kappa"]
            R.tau0 = {"r": m0["tau_r"], "f": m0["tau_f"]}
            R.sig0 = {"r": m0["tau_sig_r"], "f": m0["tau_sig_f"]}
        elif not (np.array_equal(R.names, m0["names"])
                  and np.allclose(R.tau0["r"], m0["tau_r"], equal_nan=True)):
            sys.exit(f"FATAL: {run}_{s}: module list / template taus differ from the first split")
    for a in present:
        R.arms[a] = {k: np.concatenate([p[k] for p in parts[a]]) for k in parts[a][0]}
    b = R.arms[present[0]]
    R.n = len(b["event"])
    R.refit = b["ft_status"] >= 0
    R.amp = b["fx_amp"]
    R.bi = bin_index(R.amp)
    R.eraw = match_raw(R, b, raws)
    print(f"[diag] run {run}: arms {present}, splits {use}, {R.n} mf rows, "
          f"{int(R.refit.sum())} refit")
    return R


# ---------------------------------------------------------------------- diag --
def shape_stats(R, a, w, excl_bound):
    """Per amplitude bin: spread, median sig / sigw, median shift, counts for tau_w.

    Pulses: status 0 (and, if excl_bound, no bound hit); reference = per-module median of
    the selected pulses at A > A_HI (modules with >= MIN_MOD of them only)."""
    D = R.arms[a]
    sel = (D["ft_status"] == 0) & (R.bi >= 0)
    if excl_bound:
        sel &= D["ft_bound"] == 0
    mi = D["midx"][sel].astype(np.int64)
    tau = D[f"ft_tau_{w}"][sel]
    hi = R.amp[sel] > A_HI
    med = module_median(mi[hi], tau[hi], len(R.names), MIN_MOD)
    ok = np.isfinite(med[mi])
    mi, tau, bi = mi[ok], tau[ok], R.bi[sel][ok]
    sg, sw = D[f"ft_sig_tau_{w}"][sel][ok], D[f"ft_sigw_tau_{w}"][sel][ok]
    spread = bin_stat(bi, tau - med[mi], rsig)
    sig = bin_stat(bi, sg, np.median)
    sigw = bin_stat(bi, sw, np.median)
    return dict(n=counts(bi), spread=spread, sig=sig, sigw=sigw,
                shift=bin_stat(bi, tau - R.tau0[w][mi], np.median),
                excess=np.sqrt(np.maximum(spread ** 2 - sig ** 2, 0.0)),
                excessw=np.sqrt(np.maximum(spread ** 2 - sigw ** 2, 0.0)),
                pull=spread / sig)


def diag_shape(R, plots, tabs, out):
    xc, labels = bin_centers(), bin_labels()
    variants = [(a, ex) for a in R.arms for ex in ((False, True) if a == "bounded" else (False,))]
    out.append(f"\n== (a) shape variation vs noise, run {R.run} ==")
    out.append(f"status 0 fits; spread = 1.4826*MAD of (tau_fit - module median), module median "
               f"from the selected fits at A>{A_HI:g} of modules with >= {MIN_MOD} of them; "
               f"sig = median per-fit noise-only sigma, sigw = incl. 3% model floor; "
               f"excess = sqrt(max(spread^2-sig^2,0)); shift = median(tau_fit - template). ns.")
    for w, wn in (("r", "tau_r"), ("f", "tau_f")):
        fig, ax = plt.subplots(1, 4, figsize=(17, 4))
        for a, ex in variants:
            st = shape_stats(R, a, w, ex)
            tag = a + ("_nobound" if ex else "")
            ls = "--" if ex else "-"
            ax[0].plot(xc, st["spread"], ls, marker="o", ms=3, color=COL[a], label=f"{tag} spread")
            ax[0].plot(xc, st["sig"], ":", color=COL[a], label=f"{tag} sig (noise)")
            ax[0].plot(xc, st["sigw"], "-.", color=COL[a], alpha=.6, label=f"{tag} sigw (as fit)")
            ax[1].plot(xc, st["excess"], ls, marker="o", ms=3, color=COL[a], label=tag)
            ax[2].plot(xc, st["pull"], ls, marker="o", ms=3, color=COL[a], label=tag)
            ax[3].plot(xc, st["shift"], ls, marker="o", ms=3, color=COL[a], label=tag)
            write_tsv(tabs / f"shape_{R.run}_{tag}_{wn}.tsv",
                      ["amp_bin", "n", "spread", "sig_noise", "sigw", "shift", "excess",
                       "excess_vs_sigw", "spread_over_sig"],
                      [(labels[k], int(st["n"][k]), st["spread"][k], st["sig"][k],
                        st["sigw"][k], st["shift"][k], st["excess"][k], st["excessw"][k],
                        st["pull"][k]) for k in range(len(xc))])
            out.append(f"-- {wn} {tag}")
            out.append("  bin          n   spread    sig     sigw   excess  spread/sig   shift")
            for k in range(len(xc)):
                out.append(f"  {labels[k]:<10} {int(st['n'][k]):6d} {fmt(st['spread'][k])} "
                           f"{fmt(st['sig'][k])} {fmt(st['sigw'][k])} {fmt(st['excess'][k])} "
                           f"{fmt(st['pull'][k], 2):>8} {fmt(st['shift'][k])}")
        ax[2].axhline(1, color="gray", lw=.8)
        ax[3].axhline(0, color="gray", lw=.8)
        for x, t in zip(ax, (f"{wn} spread / per-fit sigma (ns)", f"{wn} excess (ns)",
                             f"{wn} spread/sig", f"{wn} median(fit - template) (ns)")):
            x.set_xscale("log")
            x.set_xlabel("mf amplitude (ADC)")
            x.set_title(t, fontsize=9)
            x.grid(alpha=.3)
        ax[0].set_yscale("log")
        ax[0].legend(fontsize=6)
        ax[1].legend(fontsize=7)
        fig.suptitle(f"run {R.run}: shape variation vs noise ({wn})")
        savefig(fig, plots / f"shape_vs_noise_{R.run}_{wn}.png")
    # between-module vs within-module at A > A_HI
    out.append(f"\n-- between vs within modules, A>{A_HI:g}, status 0 (ns)")
    for a in R.arms:
        D = R.arms[a]
        for w in ("r", "f"):
            for ex in ((False, True) if a == "bounded" else (False,)):
                sel = (D["ft_status"] == 0) & (R.amp > A_HI)
                if ex:
                    sel &= D["ft_bound"] == 0
                mi = D["midx"][sel].astype(np.int64)
                tau = D[f"ft_tau_{w}"][sel]
                med = module_median(mi, tau, len(R.names), MIN_MOD)
                okm = np.isfinite(med[mi])
                mods = np.nonzero(np.isfinite(med))[0]
                msig = float(np.median(D[f"ft_sig_tau_{w}"][sel][okm])) if okm.any() else NAN
                tsig = float(np.median(R.sig0[w][mods])) if mods.size else NAN
                out.append(f"  {a}{'_nobound' if ex else ''} tau_{w}: n_mod={mods.size} "
                           f"within={fmt(rsig((tau - med[mi])[okm]))} "
                           f"between(median of modules)={fmt(rsig(med[mods]))} "
                           f"rms(med_fit - template)={fmt(rsig(med[mods] - R.tau0[w][mods]))} "
                           f"median per-fit sig={fmt(msig)} template sigma={fmt(tsig)}")


def diag_scatter(R, plots, nmod=6):
    base = R.arms[next(iter(R.arms))]
    ok = base["ft_status"] == 0
    cnt = np.bincount(base["midx"][ok].astype(np.int64), minlength=len(R.names))
    top = np.argsort(-cnt)[:nmod]
    xc = bin_centers()
    rng = np.random.default_rng(1)
    for a, D in R.arms.items():
        fig, ax = plt.subplots(2, len(top), figsize=(3.3 * len(top), 6), squeeze=False)
        for j, m in enumerate(top):
            sel = (D["ft_status"] == 0) & (D["midx"] == m) & (R.bi >= 0)
            idx = np.nonzero(sel)[0]
            if idx.size > 3000:
                idx = rng.choice(idx, 3000, replace=False)
            bi = R.bi[sel]
            for i, w in enumerate(("r", "f")):
                x = ax[i, j]
                t0, s0 = R.tau0[w][m], R.sig0[w][m]
                y = D[f"ft_tau_{w}"]
                x.scatter(R.amp[idx], y[idx], s=2, alpha=.3, color=COL[a], rasterized=True)
                x.axhline(t0, color="k", lw=1)
                x.axhspan(t0 - 3 * s0, t0 + 3 * s0, color="gray", alpha=.15)
                ms = bin_stat(bi, D[f"ft_sig_tau_{w}"][sel], np.median, 10)
                x.step(xc, t0 + ms, where="mid", color="r", lw=1)
                x.step(xc, t0 - ms, where="mid", color="r", lw=1)
                yy = y[sel]
                x.plot(xc, bin_stat(bi, yy, np.median, 10), "k--", lw=1)
                lo, hi = np.nanpercentile(yy, [0.5, 99.5]) if yy.size else (t0, t0)
                x.set_ylim(min(t0 - 4 * s0, lo), max(t0 + 4 * s0, hi))
                x.set_xscale("log")
                x.set_title(f"{R.names[m]} tau_{w} (n={int(cnt[m])})", fontsize=8)
                x.set_xlabel("mf amplitude (ADC)", fontsize=7)
        fig.suptitle(f"run {R.run} arm {a}: fitted tau (dots), template +-3 sigma (gray), "
                     f"median per-fit +-sig (red), median fit (dashed)", fontsize=9)
        savefig(fig, plots / f"tau_scatter_{R.run}_{a}.png")


def diag_chi2(R, plots, tabs, out):
    xc, labels = bin_centers(), bin_labels()
    per = {a: (D["ft_status"] == 0) & (R.bi >= 0) & np.isfinite(D["chi2dof"])
           & np.isfinite(D["fx_chi2dof"]) for a, D in R.arms.items()}
    com = np.logical_and.reduce(list(per.values()))
    fx = R.arms[next(iter(R.arms))]["fx_chi2dof"]
    out.append(f"\n== (b) chi2/dof vs amplitude, run {R.run} (status-0 fits only; failed rows keep "
               f"the fixed value and are excluded) ==")
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    rows = []
    sets = (("common subset (status 0 in all arms)", "common",
             [("fixed", fx, com, "k", "-")] + [(a, D["chi2dof"], com, COL[a], "-")
                                               for a, D in R.arms.items()]),
            ("per-arm status-0 subset", "per_arm",
             [(f"fixed|{a}", fx, per[a], COL[a], "--") for a in R.arms]
             + [(a, D["chi2dof"], per[a], COL[a], "-") for a, D in R.arms.items()]))
    for c, (tag, key, curves) in enumerate(sets):
        for name, v, m, col, ls in curves:
            med = bin_stat(R.bi[m], v[m], np.median)
            p84 = bin_stat(R.bi[m], v[m], pct(84))
            lab = name if "|" not in name else f"fixed (on {name[6:]} subset)"
            ax[0, c].plot(xc, med, ls, marker="o", ms=3, color=col, label=lab)
            ax[1, c].plot(xc, p84, ls, marker="o", ms=3, color=col, label=lab)
            n = counts(R.bi, m)
            for k in range(len(xc)):
                rows.append((key, name, labels[k], int(n[k]), med[k], p84[k]))
        ax[0, c].set_title(f"median chi2/dof, {tag}", fontsize=9)
        ax[1, c].set_title(f"84th percentile chi2/dof, {tag}", fontsize=9)
    for x in ax.ravel():
        x.set_xscale("log")
        x.set_yscale("log")
        x.set_xlabel("mf amplitude (ADC)")
        x.grid(alpha=.3)
        x.legend(fontsize=7)
    fig.suptitle(f"run {R.run}: chi2/dof (status-0 subsets only)")
    savefig(fig, plots / f"chi2_vs_amp_{R.run}.png")
    write_tsv(tabs / f"chi2_{R.run}.tsv", ["subset", "curve", "amp_bin", "n", "median", "p84"], rows)
    out.append("  subset     curve              bin          n    median      p84")
    for r in rows:
        if r[0] == "common":
            out.append(f"  {r[0]:<9} {r[1]:<18} {r[2]:<10} {r[3]:6d} {fmt(r[4], 3):>8} {fmt(r[5], 3):>8}")


def ratio_of(D):
    with np.errstate(all="ignore"):
        return D["E_units"] / D["fx_E_units"]


def diag_eratio(R, plots, tabs, out):
    xc, labels = bin_centers(), bin_labels()
    out.append(f"\n== (c) E_arm/E_fixed, run {R.run} (status 0) ==")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    rows = []
    for a, D in R.arms.items():
        m = (D["ft_status"] == 0) & (R.bi >= 0) & (D["fx_E_units"] > 0)
        r = ratio_of(D)
        q = [bin_stat(R.bi[m], r[m], pct(p)) for p in (16, 50, 84)]
        ax[0].plot(xc, q[1], marker="o", ms=3, color=COL[a], label=f"{a} median")
        ax[0].fill_between(xc, q[0], q[2], color=COL[a], alpha=.2)
        ax[1].plot(xc, q[2] - q[0], marker="o", ms=3, color=COL[a], label=a)
        n = counts(R.bi, m)
        out.append(f"-- {a}\n  bin          n     p16     p50     p84")
        for k in range(len(xc)):
            rows.append((a, labels[k], int(n[k]), q[0][k], q[1][k], q[2][k]))
            out.append(f"  {labels[k]:<10} {int(n[k]):6d} {fmt(q[0][k])} {fmt(q[1][k])} {fmt(q[2][k])}")
    ax[0].axhline(1, color="gray", lw=.8)
    ax[0].set_ylabel("E_arm / E_fixed (median, 16-84%)")
    ax[1].set_ylabel("p84 - p16")
    ax[1].set_yscale("log")
    for x in ax:
        x.set_xscale("log")
        x.set_xlabel("mf amplitude (ADC)")
        x.grid(alpha=.3)
        x.legend(fontsize=8)
    fig.suptitle(f"run {R.run}: energy ratio refit/fixed (status-0 fits)")
    savefig(fig, plots / f"eratio_vs_amp_{R.run}.png")
    write_tsv(tabs / f"eratio_{R.run}.tsv", ["arm", "amp_bin", "n", "p16", "p50", "p84"], rows)
    # vs tau deviation at A > A_HI
    fig, ax = plt.subplots(len(R.arms), 2, figsize=(10, 4 * len(R.arms)), squeeze=False)
    out.append(f"\n-- E ratio vs (tau_fit - template)/template sigma, A>{A_HI:g}, status 0")
    ed = np.linspace(-6, 6, 13)
    xm = 0.5 * (ed[1:] + ed[:-1])
    for i, (a, D) in enumerate(R.arms.items()):
        m = (D["ft_status"] == 0) & (R.amp > A_HI) & (D["fx_E_units"] > 0)
        r = ratio_of(D)[m]
        mi = D["midx"][m].astype(np.int64)
        for j, w in enumerate(("r", "f")):
            z = (D[f"ft_tau_{w}"][m] - R.tau0[w][mi]) / R.sig0[w][mi]
            x = ax[i, j]
            ok = np.isfinite(z) & np.isfinite(r)
            if ok.any():
                x.hist2d(np.clip(z[ok], -6, 6), np.clip(r[ok], .8, 1.2), bins=[48, 48],
                         range=[[-6, 6], [.8, 1.2]], cmap="viridis", norm=matplotlib.colors.LogNorm())
                bi = np.digitize(z, ed) - 1
                med = np.array([np.median(r[ok & (bi == k)]) if (ok & (bi == k)).sum() >= 20 else NAN
                                for k in range(len(xm))])
                x.plot(xm, med, "w-o", ms=3, lw=1)
                out.append(f"  {a} tau_{w}: n={int(ok.sum())} corr(z, ratio)="
                           f"{np.corrcoef(z[ok], r[ok])[0, 1]:.3f}; median ratio in z bins "
                           f"(centres {xm[0]:g}..{xm[-1]:g}): " + " ".join(fmt(v, 3) for v in med))
            x.set_title(f"{a}: tau_{w} deviation (units of template sigma)", fontsize=9)
            x.set_ylabel("E_arm/E_fixed")
    fig.suptitle(f"run {R.run}: E ratio vs tau deviation, A>{A_HI:g} (white: binned median)")
    savefig(fig, plots / f"eratio_vs_tau_{R.run}.png")


def diag_bound_fail(R, plots, tabs, out):
    xc, labels = bin_centers(), bin_labels()
    if "bounded" in R.arms:
        D = R.arms["bounded"]
        ok = (D["ft_status"] == 0) & (R.bi >= 0)
        b = D["ft_bound"].astype(np.int64)
        n = counts(R.bi, ok)
        fig, ax = plt.subplots(figsize=(6, 4))
        fr = {}
        for name, mask in (("lo_r", 1), ("hi_r", 2), ("lo_f", 4), ("hi_f", 8), ("any", 15)):
            f = counts(R.bi, ok & ((b & mask) > 0)) / np.maximum(n, 1)
            fr[name] = np.where(n >= MIN_BIN, f, NAN)
            ax.plot(xc, fr[name], marker="o", ms=3, lw=2 if name == "any" else 1, label=name)
        ax.set_xscale("log")
        ax.set_xlabel("mf amplitude (ADC)")
        ax.set_ylabel("fraction of status-0 fits at bound")
        ax.set_title(f"run {R.run}: bounded arm, bound-hit rate", fontsize=9)
        ax.grid(alpha=.3)
        ax.legend(fontsize=8)
        savefig(fig, plots / f"bound_hit_{R.run}.png")
        names = ("lo_r", "hi_r", "lo_f", "hi_f", "any")
        out.append(f"\n== (d) bound-hit rate, bounded arm, run {R.run} (of status-0 fits) ==")
        out.append("  bin          n   lo_r   hi_r   lo_f   hi_f    any")
        rows = []
        for k in range(len(xc)):
            rows.append((labels[k], int(n[k])) + tuple(fr[c][k] for c in names))
            out.append(f"  {labels[k]:<10} {int(n[k]):6d} "
                       + " ".join(f"{fmt(fr[c][k], 3):>6}" for c in names))
        write_tsv(tabs / f"bound_{R.run}.tsv", ["amp_bin", "n_status0", *names], rows)
    out.append(f"\n== (e) failure fraction, run {R.run} (of refit rows, status>=0) ==")
    fig, ax = plt.subplots(figsize=(6, 4))
    rows = []
    rf = R.refit & (R.bi >= 0)
    n = counts(R.bi, rf)
    for a, D in R.arms.items():
        out.append(f"-- {a}\n  bin          n  st1(noconv)  st2(nonfinite)  st3(box)   any")
        fr = {s: np.where(n >= MIN_BIN, counts(R.bi, rf & (D["ft_status"] == s)) / np.maximum(n, 1), NAN)
              for s in (1, 2, 3)}
        fr["any"] = np.where(n >= MIN_BIN, counts(R.bi, rf & (D["ft_status"] > 0)) / np.maximum(n, 1), NAN)
        ax.plot(xc, fr["any"], marker="o", ms=3, color=COL[a], label=f"{a} any")
        for s, ls in ((1, ":"), (2, "-."), (3, "--")):
            ax.plot(xc, fr[s], ls, color=COL[a], alpha=.7, label=f"{a} status {s}")
        for k in range(len(xc)):
            rows.append((a, labels[k], int(n[k]), fr[1][k], fr[2][k], fr[3][k], fr["any"][k]))
            out.append(f"  {labels[k]:<10} {int(n[k]):6d} {fmt(fr[1][k], 4):>10} {fmt(fr[2][k], 4):>12} "
                       f"{fmt(fr[3][k], 4):>9} {fmt(fr['any'][k], 4):>7}")
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("mf amplitude (ADC)")
    ax.set_ylabel("fraction of refit pulses")
    ax.set_title(f"run {R.run}: fit failures (1 not converged, 2 non-finite, 3 safety box)", fontsize=8)
    ax.grid(alpha=.3)
    ax.legend(fontsize=6)
    savefig(fig, plots / f"fail_frac_{R.run}.png")
    write_tsv(tabs / f"fail_{R.run}.tsv", ["arm", "amp_bin", "n_refit", "status1", "status2", "status3", "any"], rows)


def diag_raw(R, plots, tabs, out):
    """(f) mf integral energies vs the model-independent WaveAnalyzer integral E_raw."""
    xc, labels = bin_centers(), bin_labels()
    base = R.arms[next(iter(R.arms))]
    mraw = np.isfinite(R.eraw) & (R.eraw > 0) & (R.bi >= 0)
    out.append(f"\n== (f) E_mf/E_raw, run {R.run}: refit status-0 rows matched to the rawsum row of "
               f"the same (event, midx) with npk==1; {int((mraw & R.refit).sum())} of "
               f"{int(R.refit.sum())} refit rows matched; median and 1.4826*MAD ==")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    rows = []
    ranges = [("100-300", 100, 300), ("1000-4000", 1000, 4000)]
    for a, D in R.arms.items():
        m = mraw & (D["ft_status"] == 0) & (D["fx_E_units"] > 0)
        rf, ra = D["fx_E_units"] / R.eraw, D["E_units"] / R.eraw
        st = {n: (bin_stat(R.bi[m], v[m], np.median), bin_stat(R.bi[m], v[m], rsig))
              for n, v in (("fixed", rf), (a, ra))}
        n = counts(R.bi, m)
        for name, (med, rs) in st.items():
            ls = "--" if name == "fixed" else "-"
            col = COL["fixed"] if name == "fixed" else COL[a]
            ax[0].plot(xc, med, ls, marker="o", ms=3, color=col, alpha=.8 if name == "fixed" else 1,
                       label=f"fixed (on {a} subset)" if name == "fixed" else name)
            ax[1].plot(xc, rs, ls, marker="o", ms=3, color=col, label=f"{name} ({a} subset)")
        out.append(f"-- {a} (fixed evaluated on the same {a} status-0 rows)\n"
                   "  bin           n   fixed/raw  rms     |  arm/raw  rms")
        for k in range(len(xc)):
            rows.append((a, labels[k], int(n[k]), st["fixed"][0][k], st["fixed"][1][k],
                         st[a][0][k], st[a][1][k]))
            out.append(f"  {labels[k]:<10} {int(n[k]):6d}  {fmt(st['fixed'][0][k])} {fmt(st['fixed'][1][k])}"
                       f"  |  {fmt(st[a][0][k])} {fmt(st[a][1][k])}")
        for lab, lo, hi in ranges:
            s = m & (R.amp >= lo) & (R.amp < hi)
            if s.sum() >= MIN_BIN:
                out.append(f"  A {lab} (n={int(s.sum())}): fixed/raw {fmt(np.median(rf[s]))} "
                           f"(rms {fmt(rsig(rf[s]))}), {a}/raw {fmt(np.median(ra[s]))} "
                           f"(rms {fmt(rsig(ra[s]))})")
    ax[0].axhline(1, color="gray", lw=.8)
    ax[0].set_ylabel("median E_mf / E_raw")
    ax[1].set_ylabel("1.4826 MAD of E_mf / E_raw")
    ax[1].set_yscale("log")
    for x in ax:
        x.set_xscale("log")
        x.set_xlabel("mf amplitude (ADC)")
        x.grid(alpha=.3)
        x.legend(fontsize=7)
    fig.suptitle(f"run {R.run}: mf integral energy / rawsum (WaveAnalyzer) integral, npk==1")
    savefig(fig, plots / f"eraw_ratio_{R.run}.png")
    write_tsv(tabs / f"eraw_{R.run}.tsv",
              ["arm", "amp_bin", "n", "fixed_over_raw", "fixed_rms", "arm_over_raw", "arm_rms"], rows)


WAVE_TARGETS = (50, 200, 500, 1000, 2000, 3000)


def diag_waveforms(root, R, plots, out):
    """(g) six example pulses per run: data, fixed model, each arm's model, residuals."""
    from mf_pulse import DT, template, tmax_of
    s0 = R.splits[0]
    full = {a: np.load(root / "dump" / a / f"{R.run}_{s0}" / "pulses.npz") for a in R.arms}
    zb = full[next(iter(R.arms))]
    if "wave_row" not in zb.files or len(zb["wave_row"]) == 0:
        out.append(f"\n== (g) waveforms, run {R.run}: no wave rows in the dump ==")
        return
    wrow = zb["wave_row"]
    ok = (zb["method"][wrow] == 1) & np.isfinite(zb["fx_amp"][wrow])
    for z in full.values():
        ok &= z["ft_status"][wrow] == 0
    cand = np.nonzero(ok)[0]
    picked = []
    for t in WAVE_TARGETS:
        d = np.abs(np.log(zb["fx_amp"][wrow[cand]] / t))
        for c in cand[np.argsort(d, kind="stable")]:
            if int(c) not in picked:
                picked.append(int(c))
                break
    names = zb["module_names"]
    tr0, tf0 = zb["tau_r"], zb["tau_f"]
    tg = np.arange(0.0, 100 * DT, 0.5)
    ts = DT * np.arange(zb["wave"].shape[1])
    fig = plt.figure(figsize=(16, 9.5))
    gs = fig.add_gridspec(4, 3, height_ratios=[3, 1, 3, 1])
    out.append(f"\n== (g) example waveforms, run {R.run} (split {s0}; plots/waveforms_{R.run}.png) ==")
    for n_, c in enumerate(picked):
        r = int(wrow[c])
        m = int(zb["midx"][r])
        y = zb["wave"][c].astype(float) - float(zb["wave_ped"][c])
        sg = float(zb["wave_sigma"][c])
        t0f = float(zb["onset"][r]) * DT
        A_fx = float(zb["fx_amp"][r])
        mods = {"fixed": (A_fx, t0f, tr0[m], tf0[m])}
        info = {"fixed": f"chi2/dof {fmt(zb['fx_chi2dof'][r], 2)}"}
        for a, z in full.items():
            mods[a] = (float(z["amp"][r]), float(z["ft_t0"][r]), float(z["ft_tau_r"][r]),
                       float(z["ft_tau_f"][r]))
            info[a] = (f"chi2/dof {fmt(z['chi2dof'][r], 2)}  E/E_fx "
                       f"{fmt(z['E_units'][r] / z['fx_E_units'][r])}")
        i, j = divmod(n_, 3)
        axd = fig.add_subplot(gs[2 * i, j])
        axr = fig.add_subplot(gs[2 * i + 1, j], sharex=axd)
        axd.plot(ts, y, "ko", ms=3, label="data")
        axr.axhspan(-sg, sg, color="gray", alpha=.2)
        txt = []
        for name, (A, t0, tr, tf) in mods.items():
            col = COL[name]
            f = lambda t: A * template(t - t0, tr, tf, tmax_of(tr, tf))
            axd.plot(tg, f(tg), "-" if name != "fixed" else "--", color=col, lw=1.2, label=name)
            axr.plot(ts, y - f(ts), "o-", color=col, ms=2.5, lw=.8)
            txt.append(f"{name}: A={A:.0f} t0={t0:.1f} tau_r={tr:.2f} tau_f={tf:.2f}  {info[name]}")
        axd.text(0.98, 0.97, "\n".join(txt), transform=axd.transAxes, ha="right", va="top",
                 fontsize=6, family="monospace",
                 bbox=dict(fc="w", ec="0.8", alpha=.8))
        axd.set_title(f"{names[m]}  event {int(zb['event'][r])}  fx_amp {A_fx:.0f}", fontsize=8)
        axd.set_xlim(max(t0f - 30, 0), min(t0f + 140, ts[-1]))
        axd.legend(fontsize=6, loc="center right")
        axr.set_xlabel("t (ns)", fontsize=7)
        axr.set_ylabel("data - model", fontsize=7)
        axd.tick_params(labelsize=7)
        axr.tick_params(labelsize=7)
        out.append(f"  {names[m]} event {int(zb['event'][r])}: " + " | ".join(txt))
    fig.suptitle(f"run {R.run}: example pulses (pedestal-subtracted); gray band = +-ped sigma")
    savefig(fig, plots / f"waveforms_{R.run}.png")


def diag_sampled(R, S, plots, tabs, out, arm="sampled"):
    """Sampled-template arm beside fixed / bounded / prior: chi2/dof and E_mf/E_raw per
    amplitude bin (bins of the existing mf amplitude, so every arm sees the same pulses)."""
    xc, labels = bin_centers(), bin_labels()
    ok = S["ok"]
    col_arm = COL.get(arm, COL["sampled"])
    out.append(f"\n== {arm} arm, run {R.run}: {int(ok.sum())} paired rows ({int((ok & R.refit).sum())} "
               f"of {int(R.refit.sum())} refit rows); amp bins from the existing fixed-template "
               f"mf amplitude; {arm} kappa {fmt(S['kappa'], 4)} ==")
    sub = ok & (R.bi >= 0)
    fx = R.arms[next(iter(R.arms))]["fx_chi2dof"]
    per = [(D["ft_status"] == 0) & np.isfinite(D["chi2dof"]) for D in R.arms.values()]
    com = sub & np.logical_and.reduce(per) & np.isfinite(fx) & np.isfinite(S["chi2dof"])
    fvs = sub & R.refit & np.isfinite(fx) & np.isfinite(S["chi2dof"])
    sets = (("common (status 0 in all existing arms)", "common", com,
             [("fixed", fx, "k", "-")] + [(a, D["chi2dof"], COL[a], "-") for a, D in R.arms.items()]
             + [(arm, S["chi2dof"], col_arm, "-")]),
            (f"fixed vs {arm} (all refit-eligible rows)", f"fixed_vs_{arm}", fvs,
             [("fixed", fx, "k", "-"), (arm, S["chi2dof"], col_arm, "-")]))
    fig, ax = plt.subplots(2, 2, figsize=(11, 8))
    rows = []
    out.append(f"-- chi2/dof (median, p84), paired rows; ratio = fixed / {arm} median")
    for c, (tag, key, m, curves) in enumerate(sets):
        stats = {}
        for name, v, col, ls in curves:
            med = bin_stat(R.bi[m], v[m], np.median)
            p84 = bin_stat(R.bi[m], v[m], pct(84))
            stats[name] = med
            ax[0, c].plot(xc, med, ls, marker="o", ms=3, color=col, label=name)
            ax[1, c].plot(xc, p84, ls, marker="o", ms=3, color=col, label=name)
            n = counts(R.bi, m)
            for k in range(len(xc)):
                rows.append((key, name, labels[k], int(n[k]), med[k], p84[k]))
        ax[0, c].set_title(f"median chi2/dof, {tag}", fontsize=9)
        ax[1, c].set_title(f"84th percentile chi2/dof, {tag}", fontsize=9)
        out.append(f"  subset {key}: bin n " + " ".join(f"{nm}" for nm in stats) + f"  fixed/{arm}")
        n = counts(R.bi, m)
        for k in range(len(xc)):
            out.append(f"  {labels[k]:<10} {int(n[k]):6d} "
                       + " ".join(f"{fmt(stats[nm][k], 3):>8}" for nm in stats)
                       + f"  {fmt(stats['fixed'][k] / stats[arm][k], 2):>6}")
        hi = m & (R.amp > 500)
        if hi.sum() >= MIN_BIN:
            out.append(f"  A>500: n={int(hi.sum())} median chi2/dof fixed {fmt(np.median(fx[hi]), 3)} "
                       f"{arm} {fmt(np.median(S['chi2dof'][hi]), 3)} ratio "
                       f"{fmt(np.median(fx[hi]) / np.median(S['chi2dof'][hi]), 2)}")
    for x in ax.ravel():
        x.set_xscale("log")
        x.set_yscale("log")
        x.set_xlabel("mf amplitude (ADC)")
        x.grid(alpha=.3)
        x.legend(fontsize=7)
    fig.suptitle(f"run {R.run}: chi2/dof with the {arm}-template arm")
    savefig(fig, plots / f"chi2_vs_amp_{arm}_{R.run}.png")
    write_tsv(tabs / f"chi2_{arm}_{R.run}.tsv", ["subset", "curve", "amp_bin", "n", "median", "p84"], rows)
    # E_mf / E_raw (rawsum integral of the npk == 1 row of the same event and module)
    mraw = fvs & np.isfinite(R.eraw) & (R.eraw > 0) & np.isfinite(S["eraw"]) & (S["eraw"] > 0)
    mraw &= (R.arms[next(iter(R.arms))]["fx_E_units"] > 0) & (S["E_units"] > 0)
    fxE = R.arms[next(iter(R.arms))]["fx_E_units"]
    rf, rs = fxE / R.eraw, S["E_units"] / S["eraw"]
    with np.errstate(all="ignore"):
        ar = S["amp"] / R.amp
    st = {n: (bin_stat(R.bi[mraw], v[mraw], np.median), bin_stat(R.bi[mraw], v[mraw], rsig))
          for n, v in (("fixed", rf), (arm, rs))}
    amp_r = bin_stat(R.bi[mraw], ar[mraw], np.median)
    n = counts(R.bi, mraw)
    out.append(f"-- E_mf/E_raw, {int(mraw.sum())} rows with an npk==1 rawsum partner (E_units incl. kappa; "
               f"median and 1.4826*MAD); amp_ratio = median {arm} amp / fixed amp\n"
               f"  bin           n   fixed/raw  rms     |  {arm}/raw  rms  | amp_ratio")
    rows = []
    for k in range(len(xc)):
        rows.append((labels[k], int(n[k]), st["fixed"][0][k], st["fixed"][1][k],
                     st[arm][0][k], st[arm][1][k], amp_r[k]))
        out.append(f"  {labels[k]:<10} {int(n[k]):6d}  {fmt(st['fixed'][0][k])} {fmt(st['fixed'][1][k])}"
                   f"  |  {fmt(st[arm][0][k])} {fmt(st[arm][1][k])}  | {fmt(amp_r[k])}")
    for lab, lo, hi in (("100-300", 100, 300), ("300-1000", 300, 1000), ("1000-4000", 1000, 4000)):
        m = mraw & (R.amp >= lo) & (R.amp < hi)
        if m.sum() >= MIN_BIN:
            out.append(f"  A {lab} (n={int(m.sum())}): fixed/raw {fmt(np.median(rf[m]))} "
                       f"(rms {fmt(rsig(rf[m]))}), {arm}/raw {fmt(np.median(rs[m]))} "
                       f"(rms {fmt(rsig(rs[m]))})")
    fig, ax = plt.subplots(1, 2, figsize=(11, 4))
    for name, (med, rms) in st.items():
        col = COL.get(name, col_arm)
        ax[0].plot(xc, med, marker="o", ms=3, color=col, label=name)
        ax[1].plot(xc, rms, marker="o", ms=3, color=col, label=name)
    ax[0].axhline(1, color="gray", lw=.8)
    ax[0].set_ylabel("median E_mf / E_raw")
    ax[1].set_ylabel("1.4826 MAD of E_mf / E_raw")
    ax[1].set_yscale("log")
    for x in ax:
        x.set_xscale("log")
        x.set_xlabel("mf amplitude (ADC)")
        x.grid(alpha=.3)
        x.legend(fontsize=7)
    fig.suptitle(f"run {R.run}: mf integral energy / rawsum integral, fixed vs {arm}")
    savefig(fig, plots / f"eraw_ratio_{arm}_{R.run}.png")
    write_tsv(tabs / f"eraw_{arm}_{R.run}.tsv",
              ["amp_bin", "n", "fixed_over_raw", "fixed_rms", f"{arm}_over_raw", f"{arm}_rms",
               "amp_ratio"], rows)


def diag_splits(R, out):
    out.append(f"\n== per-split numbers, run {R.run} ==")
    base = R.arms[next(iter(R.arms))]
    for si, s in enumerate(R.splits):
        sp = base["split"] == si
        out.append(f"split {s}: mf rows {int(sp.sum())}, refit {int((sp & R.refit).sum())} "
                   f"(status -1 fraction {fmt(1 - (sp & R.refit).sum() / max(sp.sum(), 1))})")
        for a, D in R.arms.items():
            rf = sp & R.refit
            st = D["ft_status"]
            nrf = max(int(rf.sum()), 1)
            ok = rf & (st == 0)
            hi = ok & (R.amp > A_HI) & (D["fx_E_units"] > 0)
            r = ratio_of(D)[hi]
            out.append(f"  {a}: status1 {fmt((rf & (st == 1)).sum() / nrf)} status2 "
                       f"{fmt((rf & (st == 2)).sum() / nrf)} status3 {fmt((rf & (st == 3)).sum() / nrf)}"
                       f"  bound-any {fmt((ok & (D['ft_bound'] > 0)).sum() / max(int(ok.sum()), 1))}"
                       f"  A>{A_HI:g}: median E ratio {fmt(np.median(r) if r.size else NAN)}"
                       f" median chi2 fixed {fmt(np.nanmedian(base['fx_chi2dof'][hi]) if hi.any() else NAN, 3)}"
                       f" arm {fmt(np.nanmedian(D['chi2dof'][hi]) if hi.any() else NAN, 3)}")


def cmd_diag(args):
    root = Path(args.root)
    dest = Path(args.out) if args.out else root
    plots, tabs = dest / "plots", dest / "tables"
    out = [f"mf_fit_study diag, root {root}"]
    if args.arm:
        out[0] += f", arm {args.arm}"
        for run, splits in RUNS.items():
            R = load_run(root, run, splits)
            S = load_sampled(root, run, R, args.arm, args.arm_dir) if R is not None else None
            if S is not None:
                out.append(f"\n######## run {run}: splits {R.splits}, existing arms {list(R.arms)}")
                diag_sampled(R, S, plots, tabs, out, args.arm)
        dest.mkdir(parents=True, exist_ok=True)
        (dest / f"diag_summary_{args.arm}.txt").write_text("\n".join(out) + "\n")
        print(f"[diag] wrote {dest / f'diag_summary_{args.arm}.txt'}, plots in {plots}, tables in {tabs}")
        return
    for run, splits in RUNS.items():
        R = load_run(root, run, splits)
        if R is None:
            continue
        out.append(f"\n######## run {run}: splits {R.splits}, arms {list(R.arms)}, "
                   f"mf rows {R.n}, refit rows {int(R.refit.sum())}, fraction not refit "
                   f"(status -1) {fmt(1 - R.refit.mean())}")
        for a, D in R.arms.items():
            st = D["ft_status"][R.refit]
            out.append(f"  {a}: refit status counts 0/1/2/3 = "
                       + "/".join(str(int((st == s).sum())) for s in range(4)))
        diag_shape(R, plots, tabs, out)
        diag_scatter(R, plots)
        diag_chi2(R, plots, tabs, out)
        diag_eratio(R, plots, tabs, out)
        diag_bound_fail(R, plots, tabs, out)
        diag_raw(R, plots, tabs, out)
        diag_waveforms(root, R, plots, out)
        diag_splits(R, out)
    dest.mkdir(parents=True, exist_ok=True)
    (dest / "diag_summary.txt").write_text("\n".join(out) + "\n")
    print(f"[diag] wrote {dest / 'diag_summary.txt'}, plots in {plots}, tables in {tabs}")


# -------------------------------------------------------------------- headers --
def read_header(tsv):
    """n_phys line (n_phys, n_kept, elapsed_s) and the '# mf_shape' key=value line."""
    h = {"shape": None}
    for line in Path(tsv).read_text(errors="replace").splitlines():
        if not line.startswith("#"):
            if "n_phys" in h:
                break
            continue
        if line.startswith("# n_phys="):
            for k, v in re.findall(r"(\w+)=([\d.eE+-]+)", line):
                h[k] = float(v)
        elif line.startswith("# mf_shape "):
            h["shape"] = {k: (float(v) if re.fullmatch(r"[\d.eE+-]+", v) else v)
                          for k, v in re.findall(r"(\w+)=(\S+)", line)}
    return h


def res_files(root):
    """{(arm, run, split): path} of the existing result.tsv files."""
    d = {}
    for arm in RES_ARMS:
        for run, splits in RUNS.items():
            for s in splits:
                p = root / "res" / arm / f"{run}_{s}" / "result.tsv"
                if p.is_file():
                    d[(arm, run, s)] = p
    return d


# ------------------------------------------------------------------------ res --
def cmd_res(args):
    root = Path(args.root)
    import gate_scan_20261022 as G
    files = res_files(root)
    recs = {k: G.load_result(p, k[0], int(k[1]), int(k[2])) for k, p in files.items()}
    recs = {k: v for k, v in recs.items() if v is not None}
    md = ["# Resolution table: free-tau pulse fit arms (plan_20261023)", "",
          "sigma_E/E in %, post-calibration (elastic-window calibrated block, inner ring), "
          "path legacy = no gate, new = seed-time gated.  Delta = arm - fixed, paired per split "
          "(percentage points); mean +- SE over splits (SE = std(ddof=1)/sqrt(n), NA for n<2).", ""]
    tsv = []
    for run, splits in RUNS.items():
        have = [s for s in splits if all((a, run, s) in recs for a in RES_ARMS)]
        if not have:
            md += [f"## run {run}: no complete split (need fixed, bounded, prior)", ""]
            print(f"[res] run {run}: no complete split")
            continue
        r0 = recs[("fixed", run, have[0])]
        pref = r0["pref"]
        md += [f"## run {run} ({r0['rt']}), preferred fit: "
               f"{'Crystal Ball' if pref == 'cb' else 'Gaussian'}, splits {', '.join(have)}", ""]
        for kname, kd, title in (("pref", None, f"preferred ({pref})"),):
            for path in ("legacy", "new"):
                md += [f"### {title}, path {path}", "",
                       "| split | fixed | bounded | prior | d(bounded) | d(prior) |",
                       "|---|---|---|---|---|---|"]
                sv = {a: [G.sigma(recs[(a, run, s)], path, kd) * 100 for s in have] for a in RES_ARMS}
                for i, s in enumerate(have):
                    d = {a: sv[a][i] - sv["fixed"][i] for a in ("bounded", "prior")}
                    md.append(f"| {s} | {fmt(sv['fixed'][i], 3)} | {fmt(sv['bounded'][i], 3)} | "
                              f"{fmt(sv['prior'][i], 3)} | {fmt(d['bounded'], 3)} | {fmt(d['prior'], 3)} |")
                    for a in RES_ARMS:
                        tsv.append((run, kname, path, s, a, sv[a][i], sv[a][i] - sv["fixed"][i]))
                row = ["mean +- SE"]
                for a in RES_ARMS:
                    m, se, _ = G.mean_se(sv[a])
                    row.append(f"{fmt(m, 3)} +- {fmt(se, 3)}")
                    tsv.append((run, kname, path, "mean", a, m, NAN))
                    tsv.append((run, kname, path, "se", a, se, NAN))
                for a in ("bounded", "prior"):
                    m, se, _ = G.mean_se([x - y for x, y in zip(sv[a], sv["fixed"])])
                    row.append(f"{fmt(m, 3)} +- {fmt(se, 3)}")
                    tsv.append((run, kname, path, "mean_delta", a, m, NAN))
                    tsv.append((run, kname, path, "se_delta", a, se, NAN))
                md += ["| " + " | ".join(row) + " |", ""]
        # yields
        md += [f"### yields, run {run}", "",
               "n_kept (header), n_inner_events (calibrated block, legacy / new path) and n_in_fit "
               "(headline fit, legacy) per arm, order fixed / bounded / prior, arm/fixed in brackets.",
               "", "| split | n_kept | n_inner legacy | n_inner new | n_in_fit legacy |",
               "|---|---|---|---|---|"]
        qs = ("kept", "inner_legacy", "inner_new", "fit")
        tot = {(a, q): 0.0 for a in RES_ARMS for q in qs}
        for s in have:
            for a in RES_ARMS:
                r = recs[(a, run, s)]
                v = dict(kept=r["hdr"]["n_kept"],
                         inner_legacy=r["cal"]["legacy"].get("n_inner_events", NAN),
                         inner_new=r["cal"]["new"].get("n_inner_events", NAN),
                         fit=r["head"]["legacy"].get("n_in_fit", NAN))
                for q in qs:
                    tot[(a, q)] += v[q] if np.isfinite(v[q]) else 0.0
                    tsv.append((run, "yield_" + q, "-", s, a, v[q], NAN))
            md.append(f"| {s} | " + " | ".join(_trio(recs, run, s, q) for q in qs) + " |")
        trio = lambda q: (f"{fmt(tot[('fixed', q)], 0)} / " + " / ".join(
            f"{fmt(tot[(a, q)], 0)} ({fmt(tot[(a, q)] / tot[('fixed', q)] if tot[('fixed', q)] else NAN)})"
            for a in ("bounded", "prior")))
        md.append("| total | " + " | ".join(trio(q) for q in qs) + " |")
        for q in qs:
            for a in RES_ARMS:
                tsv.append((run, "yield_" + q, "-", "total", a, tot[(a, q)], NAN))
        md.append("")
    # job level / fit statistics
    md += ["## Job-level numbers and fit statistics from the result headers", "",
           "Summed over the splits of each run found; fail = n_fail/n_fit (refit pulses), st1-3 = "
           "n_status*/n_fit, bound fractions = n_at_*/(n_fit - n_fail) (bits are set on converged "
           "fits only).  t_mf_s = batched MFReconstructor.run, t_fit_s = FreeTauFitter.fit, elapsed_s "
           "= job wall time (machine may be loaded, see cost.txt).", "",
           "| run | arm | splits | n_phys | elapsed_s | n_chan | n_fit | fail | st1 | st2 | st3 | "
           "lo_r | hi_r | lo_f | hi_f | any bound | skip_multi | t_mf_s | t_fit_s |",
           "|" + "---|" * 19]
    for run in RUNS:
        for a in RES_ARMS:
            ks = [k for k in files if k[0] == a and k[1] == run]
            if not ks:
                continue
            hs = [read_header(files[k]) for k in ks]
            el = sum(h.get("elapsed_s", 0.0) for h in hs)
            nph = sum(h.get("n_phys", 0.0) for h in hs)
            sh = [h["shape"] for h in hs if h["shape"]]
            head = f"| {run} | {a} | {len(ks)} | {nph:.0f} | {el:.1f} |"
            if not sh:
                md.append(head + " NA |" * 14)
                continue
            S = lambda k: sum(x.get(k, 0.0) for x in sh)
            nf = max(S("n_fit"), 1.0)
            nc = max(nf - S("n_fail"), 1.0)
            md.append(head + f" {S('n_chan'):.0f} | {S('n_fit'):.0f} | " + " | ".join(
                fmt(v, 4) for v in (S("n_fail") / nf, S("n_status1") / nf, S("n_status2") / nf,
                                    S("n_status3") / nf, *(S(f"n_at_{b}") / nc for b in
                                                           ("lo_r", "hi_r", "lo_f", "hi_f")),
                                    S("n_any_bound") / nc))
                      + f" | {S('n_skip_multi'):.0f} | {S('t_mf_s'):.2f} | {S('t_fit_s'):.2f} |")
    md.append("")
    root.mkdir(parents=True, exist_ok=True)
    (root / "resolution_table.md").write_text("\n".join(md) + "\n")
    write_tsv(root / "resolution_table.tsv",
              ["run", "kind", "path", "split", "arm", "value", "delta_vs_fixed"], tsv)
    print("\n".join(md))
    print(f"[res] wrote {root / 'resolution_table.md'} and .tsv")


def _trio(recs, run, s, q):
    def v(a):
        r = recs[(a, run, s)]
        return {"kept": r["hdr"]["n_kept"],
                "inner_legacy": r["cal"]["legacy"].get("n_inner_events", NAN),
                "inner_new": r["cal"]["new"].get("n_inner_events", NAN),
                "fit": r["head"]["legacy"].get("n_in_fit", NAN)}[q]
    f = v("fixed")
    return f"{fmt(f, 0)} / " + " / ".join(
        f"{fmt(v(a), 0)} ({fmt(v(a) / f if f else NAN)})" for a in ("bounded", "prior"))


# ----------------------------------------------------------------------- cost --
def tmin(fn, rep):
    """(wall, cpu) seconds of the fastest of rep calls of fn()."""
    best = (float("inf"), 0.0)
    for _ in range(rep):
        w0, c0 = time.perf_counter(), time.process_time()
        fn()
        w, c = time.perf_counter() - w0, time.process_time() - c0
        if w < best[0]:
            best = (w, c)
    return best


def cmd_cost(args):
    root = Path(args.root)
    from prad2py import dec
    import mf_pulse
    import mf_fit_pulse
    dump = Path(args.dump) if args.dump else next(iter(sorted((root / "dump").glob("*/*/pulses.npz"))), None)
    if dump is None or not dump.is_file():
        sys.exit("cost: no dump (use --dump PATH)")
    z = np.load(dump)
    wr = z["wave_row"][:args.n]
    Y = z["wave"][:len(wr)]
    midx = z["midx"][wr].astype(np.int64)
    names = [str(x) for x in z["module_names"]]
    pwo4 = np.isfinite(z["tau_r"])
    cfg = dec.load_daq_config()
    wcfg = cfg.wave_cfg
    path = str(HERE / "output" / "pulse_templates_h50.dat")
    col = mf_pulse.resolve_column(path, "auto", prefer="25320")
    bank, _ = mf_pulse.build_bank(
        names, pwo4, path, dict(smooth_order=int(wcfg.smooth_order),
                                int_tail_ratio=float(wcfg.int_tail_ratio),
                                tail_break_n=int(wcfg.tail_break_n)), column=col)
    rec = mf_pulse.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4, r_model=0.3,
                                   ped_nsamples=int(wcfg.ped_nsamples))
    sr, sf, _ = mf_fit_pulse.load_tau_sigma(path, col, names, pwo4)
    n = len(Y)
    l0 = os.getloadavg()
    out = [f"mf_fit_study cost, dump {dump}",
           f"run {int(z['run'])}, {n} wave rows (single-peak unclipped pulses above the dump's "
           f"amplitude threshold), template column {col}",
           f"machine load average (1/5/15 min) at start: {l0[0]:.1f} {l0[1]:.1f} {l0[2]:.1f}, "
           f"cpu count {os.cpu_count()}"]
    wa = dec.WaveAnalyzer(wcfg)
    rows = [np.ascontiguousarray(Y[i]) for i in range(n)]
    w_wa, c_wa = tmin(lambda: [wa.analyze(r) for r in rows], args.rep)
    ov = int(wcfg.overflow)
    res = rec.run(Y, midx, ov)
    w_mf, c_mf = tmin(lambda: rec.run(Y, midx, ov), args.rep)
    pc = res.p_ch
    sel = np.nonzero(~res.wa_fb[pc] & (res.npk[pc] == 1) & ~res.clipped[pc])[0]
    row = pc[sel]
    y = Y[row].astype(np.float32) - res.ped[row][:, None]
    npul = len(sel)
    out += [f"mf found {len(pc)} pulses, {npul} single unclipped (fit input)", "",
            f"{'stage':<26}{'wall us/unit':>14}{'cpu us/unit':>14}   unit",
            f"{'rawsum (WaveAnalyzer)':<26}{w_wa / n * 1e6:14.2f}{c_wa / n * 1e6:14.2f}   channel",
            f"{'mf (MFReconstructor.run)':<26}{w_mf / n * 1e6:14.2f}{c_mf / n * 1e6:14.2f}   "
            f"channel (batch {n})"]
    fit = {}
    for mode in DUMP_ARMS:
        fi = mf_fit_pulse.FreeTauFitter(bank, sr, sf, mode=mode, nsig=3.0)
        args_fit = (y, midx[row], res.p_t0[sel], res.p_amp[sel], res.sigma[row])
        fr = fi.fit(*args_fit)
        w, c = tmin(lambda: fi.fit(*args_fit), args.rep)
        st = np.bincount(np.asarray(fr.status), minlength=4)
        fit[mode] = (w, c)
        out.append(f"{'free fit ' + mode:<26}{w / npul * 1e6:14.2f}{c / npul * 1e6:14.2f}   pulse "
                   f"(batch {npul}; status 0/1/2/3 = {'/'.join(map(str, st))}, "
                   f"mean niter {float(np.mean(fr.niter)):.1f})")
    out += ["", "ratios fit(per pulse) / other(per channel), wall (cpu in brackets):"]
    for mode in DUMP_ARMS:
        w, c = fit[mode][0] / npul, fit[mode][1] / npul
        out.append(f"  fit({mode})/mf = {w / (w_mf / n):.2f} ({c / (c_mf / n):.2f})   "
                   f"fit({mode})/rawsum = {w / (w_wa / n):.2f} ({c / (c_wa / n):.2f})")
    out.append(f"  mf/rawsum = {w_mf / w_wa:.2f} wall ({c_mf / c_wa:.2f} cpu)")
    l1 = os.getloadavg()
    out += ["note: every row is a pulse channel that is also refit; a job also runs mf/rawsum on "
            "pulse-free channels, so job-average per-channel costs are lower.",
            f"machine load average at end: {l1[0]:.1f} {l1[1]:.1f} {l1[2]:.1f}; fastest of {args.rep} "
            f"repeats, single process, numpy threads default"]
    files = res_files(root)
    out += ["", "job-level view from the result headers (all result.tsv found):",
            "  arm      run     split  elapsed_s  t_mf_s  n_chan  us/chan(mf)  t_fit_s  n_fit  us/pulse(fit)"]
    for (arm, r, s), p in sorted(files.items()):
        h = read_header(p)
        sh = h["shape"] or {}
        tm, nc = sh.get("t_mf_s", NAN), sh.get("n_chan", NAN)
        tf, nf = sh.get("t_fit_s", NAN), sh.get("n_fit", NAN)
        out.append(f"  {arm:<8} {r}  {s}  {fmt(h.get('elapsed_s', NAN), 1):>8} {fmt(tm, 2):>7} "
                   f"{fmt(nc, 0):>7} {fmt(tm / nc * 1e6, 2):>11} {fmt(tf, 2):>8} {fmt(nf, 0):>6} "
                   f"{fmt(tf / nf * 1e6, 2):>12}")
    root.mkdir(parents=True, exist_ok=True)
    (root / "cost.txt").write_text("\n".join(out) + "\n")
    print("\n".join(out))


# --------------------------------------------------------------------- intraw --
INTRAW_RUNS = dict(RUNS, **{"024590": ["00010"]})
DEFAULT_BANK = "output/templates_spline/{runtype}/templates_{runtype}.npz"
INT_BINS = [(5, 20), (20, 100), (100, 300), (300, 1000), (1000, 4000)]
SEED_SUF = " (seedmodel weights)"


def _arm_style(a):
    return COL.get(a.split()[0], COL["sampled"]), ("--" if "seedmodel" in a else "-")


def intraw_load(root, run, splits, arm=None, arm_dir=None):
    """Common-set Q_template (A*I_int, kappa removed) and Q_raw in ADC*ns for one run.
    arm: add that template arm (dumps of sampled_path), restrict every arm to the pulses paired to
    it, and skip the seed-model-weight variants."""
    from mf_pulse import LEAD, LW
    sampled = arm is not None
    out = []
    for si, s in enumerate(splits):
        pth = {a: root / "dump" / a / f"{run}_{s}" / "pulses.npz" for a in DUMP_ARMS}
        if not all(x.is_file() for x in pth.values()):
            continue
        ps_ = sampled_path(root, run, s, arm, arm_dir) if sampled else None
        if sampled and ps_ is None:
            continue
        z = {a: np.load(x) for a, x in pth.items()}
        zb = z["bounded"]
        ev, mi = zb["event"].astype(np.int64), zb["midx"].astype(np.int64)
        if not (np.array_equal(ev, z["prior"]["event"]) and np.array_equal(mi, z["prior"]["midx"])):
            sys.exit(f"FATAL: {run}_{s}: bounded/prior pulse tables differ")
        kap = float(zb["kappa"])
        ns = zb["wave"].shape[1] if zb["wave"].ndim == 2 else 100
        mf = (zb["method"] == 1) & (zb["ft_status"] >= 0)
        idx = np.nonzero(mf)[0]
        sel = (z["bounded"]["ft_status"][idx] == 0) & (z["prior"]["ft_status"][idx] == 0) \
            & ~zb["clipped"][idx]
        w0 = np.ceil(zb["onset"][idx]) - LEAD
        sel &= (w0 >= 0) & (w0 + LW <= ns)
        ridx = np.nonzero((zb["method"] == 0) & (zb["npk"] == 1))[0]
        rk = (ev[ridx] << 12) | mi[ridx]
        o = np.argsort(rk, kind="stable")
        rk, ridx = rk[o], ridx[o]
        k = (ev[idx] << 12) | mi[idx]
        i = np.clip(np.searchsorted(rk, k), 0, max(rk.size - 1, 0))
        r = ridx[i]
        pos = zb["onset"][r]
        sel &= (rk[i] == k) & ~zb["clipped"][r] & (pos >= 3) & (pos <= ns - 4) & (zb["E_units"][r] > 0)
        d = dict(A=zb["fx_amp"][idx][sel], midx=mi[idx][sel], split=np.full(int(sel.sum()), si),
                 Q={"fixed": zb["fx_E_units"][idx][sel] / kap * 4.0},
                 raw=zb["E_units"][r][sel] * 4.0, nmod=len(zb["module_names"]), names=zb["module_names"])
        for a in ("prior", "bounded"):
            d["Q"][a] = z[a]["E_units"][idx][sel] / kap * 4.0
        if sampled:
            ds, _raw, kap_s, names_s = _load_sampled_split(ps_)
            if not np.array_equal(names_s, zb["module_names"]):
                sys.exit(f"FATAL: {ps_}: module list differs from the existing arms")
            j = pair_rows((ev[idx][sel] << 12) | mi[idx][sel], zb["time_ns"][idx][sel],
                          (ds["event"].astype(np.int64) << 12) | ds["midx"].astype(np.int64),
                          ds["time_ns"])
            hit = j >= 0
            d["n_sel"] = int(sel.sum())
            for k in ("A", "midx", "split", "raw"):
                d[k] = d[k][hit]
            d["Q"] = {a: v[hit] for a, v in d["Q"].items()}
            d["Q"][arm] = ds["E_units"][j[hit]] / kap_s * 4.0
            print(f"[intraw] {run}_{s}: {arm} pairs {int(hit.sum())} of {d['n_sel']} common-set pulses")
        for a in () if sampled else DUMP_ARMS:       # optional seed-model-weight variant
            ps = root.parent / "dump" / a / f"{run}_{s}" / "pulses.npz"
            if ps.is_file():
                zs = np.load(ps)
                if np.array_equal(zs["event"], zb["event"]) and np.array_equal(zs["midx"], zb["midx"]):
                    ok = zs["ft_status"][idx][sel] == 0
                    d["Q"][a + SEED_SUF] = np.where(ok, zs["E_units"][idx][sel] / kap * 4.0, NAN)
        out.append(d)
    if not out:
        return None
    arms = [a for a in out[0]["Q"] if all(a in d["Q"] for d in out)]
    R = {k: np.concatenate([d[k] for d in out]) for k in ("A", "midx", "raw", "split")}
    R["Q"] = {a: np.concatenate([d["Q"][a] for d in out]) for a in arms}
    R["nmod"] = out[0]["nmod"]
    R["names"] = out[0]["names"]
    R["order"] = [a for a in ("fixed", "prior", "bounded", arm)
                  + tuple(x + SEED_SUF for x in ("prior", "bounded")) if a in arms]
    return R


def intraw_plots(run, R, out_dir, md, sfx=""):
    A, raw, Q, arms = R["A"], R["raw"], R["Q"], R["order"]
    ratio = {a: Q[a] / raw for a in arms}
    bi, xc, labels = bin_index(A), bin_centers(), bin_labels()
    kap_note = "kappa removed: Q = A*I_int, ADC*ns"
    # 1. 2-D density
    fig, ax = plt.subplots(1, len(arms), figsize=(4.2 * len(arms), 5.0), squeeze=False)
    ok0 = np.isfinite(Q["fixed"]) & (Q["fixed"] > 0)
    lo, hi = np.percentile(raw[ok0], [0.05, 99.99])
    be = np.logspace(np.log10(max(lo, 10)), np.log10(hi), 90)
    for x, a in zip(ax[0], arms):
        m = np.isfinite(Q[a]) & (Q[a] > 0)
        h = x.hist2d(raw[m], Q[a][m], bins=[be, be], cmap="viridis", norm=matplotlib.colors.LogNorm())
        x.plot(be, be, "r-", lw=.8)
        x.set(xscale="log", yscale="log", xlabel="Q_rawsum (ADC*ns)", title=f"{a}  n={int(m.sum())}")
        x.set_aspect("equal")
    ax[0][0].set_ylabel("Q_template (ADC*ns)")
    fig.suptitle(f"run {run}: template integral vs rawsum integral ({kap_note})", fontsize=9)
    savefig(fig, out_dir / f"q2d{sfx}_{run}.png")
    # 2. ratio and difference vs amplitude + low-amplitude loss marker
    fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
    med_fx = bin_stat(bi, ratio["fixed"], np.median)
    ref = np.median(ratio["fixed"][(A >= 300) & (A < 1000)])
    loss = [k for k in range(len(xc)) if AMP_EDGES[k + 1] <= 300 and np.isfinite(med_fx[k])
            and abs(med_fx[k] / ref - 1) > 0.01]
    for a in arms:
        col, ls = _arm_style(a)
        for j, v in enumerate((ratio[a], Q[a] - raw)):
            q = [bin_stat(bi, v, pct(p)) for p in (16, 50, 84)]
            ax[j].plot(xc, q[1], ls, marker="o", ms=3, color=col, label=a)
            ax[j].fill_between(xc, q[0], q[2], color=col, alpha=.12 if ls == "-" else .0)
    if loss:
        edge = AMP_EDGES[max(loss) + 1]
        for x in ax:
            x.axvline(edge, color="gray", ls=":", label=f"rawsum loss edge {edge:g} ADC" if x is ax[0] else None)
        rule = (f"rawsum loss edge = highest bin edge <= 300 ADC below which the fixed-arm median ratio "
                f"differs by >1% from its 300-1000 ADC median ({ref:.4f})")
    else:
        rule = "no rawsum loss edge found (fixed-arm median ratio within 1% of its 300-1000 ADC median)"
    fig.text(0.5, 0.005, rule, ha="center", fontsize=7)
    ax[0].axhline(1, color="gray", lw=.8)
    ax[0].set_ylabel("Q_template / Q_rawsum (median, 16-84%)")
    ax[1].axhline(0, color="gray", lw=.8)
    ax[1].set_ylabel("Q_template - Q_rawsum (ADC*ns)")
    ax[1].set_yscale("symlog", linthresh=100)
    for x in ax:
        x.set(xscale="log", xlabel="mf amplitude A (ADC)")
        x.grid(alpha=.3)
    ax[0].legend(fontsize=7, loc="upper right")
    fig.suptitle(f"run {run}: integrated template vs rawsum ({kap_note})", fontsize=9)
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(out_dir / f"ratio_vs_amp{sfx}_{run}.png", dpi=110)
    plt.close(fig)
    # 3. ratio histograms
    fig, ax = plt.subplots(1, len(INT_BINS), figsize=(3.6 * len(INT_BINS), 3.6))
    for x, (lo_, hi_) in zip(ax, INT_BINS):
        m = (A >= lo_) & (A < hi_)
        for a in arms:
            col, ls = _arm_style(a)
            v = ratio[a][m]
            x.hist(v[np.isfinite(v)], bins=100, range=(.6, 1.6), density=True, histtype="step",
                   color=col, ls=ls, label=a)
        x.set(title=f"A {lo_}-{hi_} ADC (n={int(m.sum())})", xlabel="Q_template / Q_rawsum")
        x.axvline(1, color="gray", lw=.6)
    ax[0].legend(fontsize=6)
    savefig(fig, out_dir / f"ratio_hist{sfx}_{run}.png")
    # 4. per-module median ratio at A > 300
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    hi_a = A > 300
    modspread = {}
    for a in arms:
        col, ls = _arm_style(a)
        v = ratio[a]
        ok = hi_a & np.isfinite(v)
        med = module_median(R["midx"][ok], v[ok], R["nmod"], 20)
        mm = med[np.isfinite(med)]
        rms = float(np.sqrt(np.mean((mm - mm.mean()) ** 2))) if mm.size else NAN
        modspread[a] = (mm.size, float(np.median(mm)) if mm.size else NAN, rms, rsig(mm))
        ax.hist(mm, bins=70, range=(.85, 1.15), histtype="step", color=col, ls=ls,
                label=f"{a}: n_mod={mm.size} RMS={rms:.4f} MAD={rsig(mm):.4f}")
    ax.set(xlabel="per-module median Q_template/Q_rawsum (A>300, >=20 pulses)", ylabel="modules")
    ax.legend(fontsize=6)
    ax.set_title(f"run {run}", fontsize=9)
    savefig(fig, out_dir / f"module_ratio{sfx}_{run}.png")
    # table
    md += [f"## run {run}", "", f"{int(len(A))} pulses in the common set, splits {sorted(set(R['split'].tolist()))}; "
           f"fixed-arm rawsum-loss edge: {AMP_EDGES[max(loss) + 1]:g} ADC" if loss else
           f"{int(len(A))} pulses in the common set; no rawsum-loss edge found", "",
           "| arm | A bin (ADC) | n | median ratio | p16 | p84 | median Q-Q_raw (ADC*ns) |", "|---|---|---|---|---|---|---|"]
    for a in arms:
        for lo_, hi_ in INT_BINS:
            m = (A >= lo_) & (A < hi_) & np.isfinite(ratio[a])
            if m.sum():
                q = np.percentile(ratio[a][m], [16, 50, 84])
                md.append(f"| {a} | {lo_}-{hi_} | {int(m.sum())} | {q[1]:.4f} | {q[0]:.4f} | {q[2]:.4f} | "
                          f"{np.median((Q[a] - raw)[m]):.1f} |")
            else:
                md.append(f"| {a} | {lo_}-{hi_} | 0 | NA | NA | NA | NA |")
    md += ["", "Module medians at A>300 (modules with >=20 pulses): spread across modules", "",
           "| arm | n_mod | median of modules | RMS | 1.4826 MAD |", "|---|---|---|---|---|"]
    md += [f"| {a} | {v[0]} | {fmt(v[1])} | {fmt(v[2])} | {fmt(v[3])} |" for a, v in modspread.items()]
    md.append("")


def parse_edges(s):
    e = [float(x) for x in s.split(",") if x.strip()]
    if len(e) < 2 or any(b <= a for a, b in zip(e, e[1:])):
        sys.exit(f"FATAL: --bins needs >= 2 increasing edges, got {s!r}")
    return e


def bank_classes(pattern, run, names):
    """Fallback class (0/1/2; -1 = module not in the bank) per module of the dump, from the
    `fallback` array (by module name) of the run type's bank; (None, path) if the bank is missing."""
    from gate_scan_20261022 import runtype_slug
    path = Path(pattern.replace("{runtype}", runtype_slug(int(run))))
    if not path.is_absolute() and not path.exists():
        path = HERE / path
    if not path.is_file():
        print(f"[intraw] run {run}: bank {path} not found; no class split")
        return None, path
    z = np.load(path)
    fb = dict(zip(z["names"].tolist(), z["fallback"].tolist()))
    return np.array([fb.get(str(n), -1) for n in names], np.int64), path


def bin_med_se(v):
    """(n, median, SE) with SE = 1.2533 * sigma / sqrt(n), sigma = 1.4826 MAD; nan below MIN_BIN."""
    v = v[np.isfinite(v)]
    if v.size < MIN_BIN:
        return int(v.size), NAN, NAN
    return int(v.size), float(np.median(v)), 1.2533 * rsig(v) / np.sqrt(v.size)


def flatness(stats, edges):
    """Over the bins with lower edge >= 100 ADC (and a statistic): (flat, spread %, worst pull)
    with flat = every pairwise |median_i - median_j| <= 2 hypot(SE_i, SE_j); spread = (max - min
    median) / min median; pull = max |difference| / hypot(SE_i, SE_j).  None for < 2 bins."""
    b = [s for s, lo in zip(stats, edges[:-1]) if lo >= 100 and np.isfinite(s[1])]
    if len(b) < 2:
        return None
    pull = max(abs(x[1] - y[1]) / np.hypot(x[2], y[2]) for i, x in enumerate(b) for y in b[i + 1:])
    med = [x[1] for x in b]
    return pull <= 2.0, (max(med) / min(med) - 1) * 100, pull


def intraw_bins(run, R, arm, edges, bank, tsv, md):
    """Rows for intraw_bins_<arm>.tsv and tables + flatness verdicts for the .md, one run."""
    A = R["A"]
    ratio = R["Q"][arm] / R["raw"]
    cls = None
    if bank:
        cm, path = bank_classes(bank, run, R["names"])
        if cm is not None:
            cls = cm[R["midx"]]
            md.append(f"Fallback classes from {path} (module counts 0/1/2: "
                      + "/".join(str(int((cm == c).sum())) for c in range(3)) + ").")
            md.append("")
    labels = [f"{lo:g}-{hi:g}" for lo, hi in zip(edges[:-1], edges[1:])]
    bi = np.digitize(A, edges) - 1
    sets = [("all", np.ones(len(A), bool))] + ([(str(c), cls == c) for c in range(3)] if cls is not None else [])
    stats = {}
    for name, m in sets:
        stats[name] = [bin_med_se(ratio[m & (bi == k)]) for k in range(len(labels))]
        for k, lab in enumerate(labels):
            tsv.append((run, name, lab, stats[name][k][0], stats[name][k][1], stats[name][k][2]))
    md.append("| A bin (ADC) | " + " | ".join(f"{n}: n | median +- SE" for n, _ in sets) + " |")
    md.append("|---|" + "---|---|" * len(sets))
    for k, lab in enumerate(labels):
        md.append(f"| {lab} | " + " | ".join(
            f"{stats[n][k][0]} | {fmt(stats[n][k][1])} +- {fmt(stats[n][k][2])}" for n, _ in sets) + " |")
    md.append("")
    for name, _ in sets:
        f = flatness(stats[name], edges)
        md.append(f"- flatness, class {name}, bins with A >= 100 (n >= {MIN_BIN}): "
                  + ("not enough bins" if f is None else
                     f"{'FLAT' if f[0] else 'NOT FLAT'}, max-min spread {f[1]:.3f} %, "
                     f"worst pairwise pull {f[2]:.2f} x hypot(SE) (limit 2)"))
    md.append("")


def cmd_intraw(args):
    root = Path(args.root).resolve()
    out_dir = Path(args.out) if args.out else root.parent / "int_vs_rawsum"
    out_dir.mkdir(parents=True, exist_ok=True)
    arm = args.arm
    smp = arm is not None
    sfx = f"_{arm}" if smp else ""
    md = ["# Integrated template vs rawsum integral" + (f", {arm}-template arm" if smp else ""), "",
          "Q_template = fx_E_units (fixed) or E_units (arm, status 0) divided by kappa (removed, so Q = "
          "A*I_int) times 4 ns; Q_raw = rawsum E_units * 4 ns of the npk==1 row of the same "
          "(split, event, midx).  Pulses: refit rows, status 0 in both arms, not clipped, mf window "
          "inside the waveform, rawsum peak position >= 3 samples from either end."
          + (f"  The {arm} arm is the E_units of a plain mf dump with the {arm}-template bank "
             f"(its own kappa removed; the bank's f(A) stays in), paired to the common-set pulses by "
             f"(event, midx) with a single mf row on both sides and |dt| <= {SAMPLED_DT_MAX:g} ns; "
             "every arm is evaluated on the pulses paired to it.  Seed-model-weight arms are left out."
             if smp else
             "  Arms marked 'seedmodel weights' come from the dump directory above the root (same pulses, "
             "fit weights from the seed model)."), ""]
    edges = parse_edges(args.bins) if args.bins else None
    barm = arm or "fixed"
    bank = args.bank or (DEFAULT_BANK if arm == "sampled" else None)
    btsv, bmd = [], [f"# Q_template / Q_raw per amplitude bin, {barm} arm", "",
                     f"Bins (ADC, [lo, hi)): {args.bins}.  median and SE = 1.2533 * sigma / sqrt(n), sigma = "
                     f"1.4826 MAD, per bin with n >= {MIN_BIN}.  Flatness: every pairwise difference of bin "
                     "medians (bins with A >= 100) <= 2 hypot(SE_i, SE_j); spread = (max median / min median "
                     "- 1) in %.  Q_template of the " + barm + " arm; Q = A*I_int (kappa removed).", ""]
    for run, splits in INTRAW_RUNS.items():
        R = intraw_load(root, run, splits, arm, args.arm_dir)
        if R is None:
            print(f"[intraw] run {run}: no complete split under {root / 'dump'}"
                  + (f" with a {arm} dump" if smp else "") + "; skipped")
            md += [f"## run {run}: no dumps found", ""]
            continue
        print(f"[intraw] run {run}: {len(R['A'])} pulses, arms {R['order']}")
        intraw_plots(run, R, out_dir, md, sfx)
        if edges:
            bmd += [f"## run {run}", "", f"{len(R['A'])} pulses in the common set.", ""]
            intraw_bins(run, R, barm, edges, bank, btsv, bmd)
    (out_dir / f"int_vs_rawsum{sfx}.md").write_text("\n".join(md) + "\n")
    print("\n".join(md))
    if edges:
        write_tsv(out_dir / f"intraw_bins_{barm}.tsv", ["run", "class", "bin", "n", "median", "se"], btsv)
        (out_dir / f"intraw_bins_{barm}.md").write_text("\n".join(bmd) + "\n")
        print("\n".join(bmd))
    print(f"[intraw] wrote {out_dir}")


# ----------------------------------------------------------------------- main --
def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    sub = ap.add_subparsers(dest="cmd", required=True)
    for name, fn in (("diag", cmd_diag), ("res", cmd_res), ("cost", cmd_cost), ("intraw", cmd_intraw)):
        p = sub.add_parser(name)
        p.add_argument("--root", default=str(DEFAULT_ROOT))
        p.set_defaults(fn=fn)
        if name in ("diag", "intraw"):
            p.add_argument("--arm", default=None,
                           help="add the template arm NAME (e.g. sampled) from <root>/dump/NAME "
                                "(or <root>/../dump/NAME); outputs are written with a "
                                "NAME label beside the existing ones")
            p.add_argument("--arm-dir", default=None,
                           help="dumps of --arm at DIR/<run>_<split>/pulses.npz")
        if name == "diag":
            p.add_argument("--out", default=None, help="default: <root> (plots/, tables/, summary)")
        if name == "intraw":
            p.add_argument("--out", default=None, help="default: <root>/../int_vs_rawsum")
            p.add_argument("--bins", default=None,
                           help="ADC bin edges E0,E1,...: per-bin median ratio, SE, fallback classes "
                                "and flatness -> <out>/intraw_bins_<arm>.{tsv,md}")
            p.add_argument("--bank", default=None,
                           help="bank npz pattern with {runtype} for the fallback classes "
                                f"(default for --arm sampled: {DEFAULT_BANK})")
        if name == "cost":
            p.add_argument("--dump", default=None, help="pulses.npz (default: first under root/dump)")
            p.add_argument("--n", type=int, default=4000, help="max wave rows")
            p.add_argument("--rep", type=int, default=3, help="repeats (fastest reported)")
    args = ap.parse_args(argv)
    args.fn(args)


if __name__ == "__main__":
    main()
