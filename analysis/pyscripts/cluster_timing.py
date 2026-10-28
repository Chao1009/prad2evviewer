#!/usr/bin/env python3
"""cluster_timing.py - RF-referenced timing of single-cluster events (seed time and energy-weighted time).

PURPOSE
    sigma_t per cluster: for every row of a benchmark result_clusters.tsv written with
    --write-cluster-tsv --cluster-times (16 columns), the cluster time relative to the accelerator RF,
    phi = fold(t - t_RF), after the per-run RF offsets of rf_time_calib.py.  Same estimator and same
    t_RF matching as the module-level analysis (plot_module_timing_rf.py, rf_time_calib.py), so the
    seed-definition sigma is directly comparable with the module sigma of module_timing_rf.tsv.
    Plan: plans/plan_best50_pipeline_20261026.md (F10-F12, section 3 "sigma_t per cluster", S6).

INPUT
    --scan-dir DIR/scan_<run6>/split_<n>/result_clusters.tsv   columns event_num config_name
        seed_module_id seed_x_mm seed_y_mm seed_row seed_col seed_energy_MeV cluster_E_raw_MeV
        cluster_E_cal_MeV cluster_x_mm cluster_y_mm evnum seed_time_ns n_members members
        (members = name:E_MeV:t_ns;...  hit times as fed to the clusterer: mf p_time, ns, no time
        calibration).  Rows of config_name == --path only.  Files without the 16 columns are skipped.
    --rf-dir  <run6>_<split5>.npz   rf_time_dump.py (keys ev_num trig ts rf_off rf_t rf_ch)
    --offsets-dir  <run5>.json (rf_time_calib.py: modules[{name, offset_ns}], usable, period_ns) and
        rf_calib_runs.tsv (P_ns, usable).  Runs with usable false (either file) or without offsets / RF
        files are skipped and logged.

DEFINITIONS (fold(x) = x - P*round(x/P), (-P/2, P/2], RF.fold; P = rf_calib_runs.tsv P_ns of the run)
  Offset sign       rf_time_calib: o_{m,r} is such that fold(t_m - t_RF - o_{m,r}) is centred at 0, so the
                    corrected time of a module is t - o.  o is defined modulo P (offset_ns is folded).
  t_RF              mirrors RF.prep: the event is found by evnum in the RF split (searchsorted on ev_num), selected
                    like the module rows: (trig & --trigger-mask) != 0, (trig & PMT.TRIG_EXCL) == 0, >= 1
                    channel-A hit (rf_ch == 0); t_RF = the channel-A hit of the event nearest in time to the
                    RAW seed time seed_time_ns (the keyed searchsorted over the (event, time)-sorted A hits of
                    RF.prep, same KEY_BIG).  Both definitions use this one t_RF (seed-anchored), as the module
                    analysis anchors on the pulse itself.
  seed              t = seed_time_ns - o_seed;               phi_seed = fold(t - t_RF).
  ewt               members i with an offset o_i and finite t_i and E_i > 0 (seed must be among them):
                    u_i = t_i - o_i, unwrapped to the seed's corrected time (u_i -> t_s + fold(u_i - t_s)),
                    because o_i is only known modulo P; t_ewt = sum E_i u_i / sum E_i;  phi_ewt = fold(t_ewt - t_RF).
                    (A plain mean of t_i - o_i, as written in the plan, would mix RF buckets of members.)
  Rows              a row enters both definitions iff: event in the RF split, trigger selected, >= 1 A hit, seed
                    time finite, seed has an offset, seed among the valid members, sum E_i > 0.  Row N is therefore
                    identical for seed and ewt.
  Tier              by the n_members column: all, 1, 2-4, 5+ (a row is in "all" and in its tier).
  Histogram         NB = 80 bins over one period (RF._bin).  Per run: sum of the per-split histograms; per type
                    (run = ALL): sum of the per-run histograms (= concatenating the offset-corrected phi),
                    P of the type = N-weighted mean of the run P.
  Fit               RF.fit_unwrapped on the histogram (unwrapped around the 5-bin smoothed argmax, fit_core +-2 sigma
                    Gaussian core, bg=False, smooth 0.25 ns); mu = fit mean, or the argmax centre on failure.  Histograms
                    with N < --min-n are not fitted (fit_ok 0, sigma empty).
  frac_members_with_offset   sum over the rows used of members with offset and finite t / sum of n_members.

OUTPUT (--out-dir)
  cluster_timing.tsv         run type path definition tier N P_ns mu_ns sigma_ns sigma_err_ns fit_ok
                             median_n_members frac_members_with_offset   (run = 5-digit run, or ALL per type)
  cluster_timing_splits.tsv  the same per split (columns run split type path ...), for pairing
  cluster_timing.pdf         one page per type: folded histograms with fits, rows seed / ewt, columns tiers
  run.log                    counts: rows read, missing RF event, no A hit, trigger, seed without offset, non-finite,
                             members dropped, ...

Usage:
  MPLBACKEND=Agg python3 cluster_timing.py --scan-dir DIR --rf-dir DIR --offsets-dir DIR --out-dir DIR
      [--path legacy|gated] [--jobs 12] [--min-n 300] [--trigger-mask 0xf00]
"""

from __future__ import annotations

import argparse
import csv
import glob
import json
import math
import multiprocessing as mp
import os
import re
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import plot_module_timing as PMT  # noqa: E402
import plot_module_timing_rf as RF  # noqa: E402
from plot_module_timing import NMAX_JOBS, TAGS, LABEL, TYPE_OF_RUN, TRIG_EXCL  # noqa: E402
from plot_module_timing_rf import NB, _bin, fit_unwrapped, fold  # noqa: E402

NCOL = 16
DEFS = ("seed", "ewt")
TIERS = ("all", "1", "2-4", "5+")
NMEM_MAX = 64                     # n_members histogram size (clipped)
COUNTERS = ("rows_read", "rows_other_path", "bad_line", "no_members", "n_members_mismatch",
            "rf_event_missing", "no_a_hit", "trig_rejected", "seed_time_nonfinite", "seed_no_offset",
            "seed_not_member", "ewt_nonfinite", "rows_used", "members_total", "members_used",
            "members_no_offset", "members_nonfinite")


def tier_of(n):
    return 1 if n == 1 else 2 if n <= 4 else 3


def new_acc():
    return dict(H=np.zeros((2, 4, NB), np.int64), nm=np.zeros((4, NMEM_MAX + 1), np.int64),
                mem_used=np.zeros(4, np.int64), mem_tot=np.zeros(4, np.int64),
                c={k: 0 for k in COUNTERS})


# ---- worker -----------------------------------------------------------------

def match_rf(rp, ev, t, mask):
    """Mirror of the RF part of RF.prep for rows (ev = TI event numbers, t = raw time used for the nearest-hit
    search).  Returns (code, t_rf): code 0 ok, 1 event not in the RF split, 2 trigger rejected, 3 no A hit."""
    r = np.load(rp)
    rf_ev, off = r["ev_num"], r["rf_off"].astype(np.int64)
    rf_t, rf_ch = r["rf_t"].astype(np.float64), r["rf_ch"]
    n_rf = len(rf_ev)
    code = np.zeros(len(ev), np.int8)
    t_rf = np.full(len(ev), np.nan)
    if n_rf == 0 or not len(ev):
        code[:] = 1
        return code, t_rf
    pos = np.minimum(np.searchsorted(rf_ev, ev), n_rf - 1)
    matched = rf_ev[pos] == ev
    ha = np.flatnonzero(rf_ch == 0)
    ev_a = np.repeat(np.arange(n_rf), np.diff(off))[ha]
    ta = rf_t[ha]
    o = np.lexsort((ta, ev_a))
    ev_a, ta = ev_a[o], ta[o]
    nA = np.bincount(ev_a, minlength=n_rf)
    trig = r["trig"][pos]
    trig_ok = ((trig & np.uint32(mask)) != 0) & ((trig & np.uint32(TRIG_EXCL)) == 0)
    code[~matched] = 1
    code[matched & ~trig_ok] = 2
    code[matched & trig_ok & (nA[pos] < 1)] = 3
    ok = code == 0
    if ok.any():
        K = ev_a * RF.KEY_BIG + ta
        ev_row, tt = pos[ok], t[ok]
        i = np.searchsorted(K, ev_row * RF.KEY_BIG + tt)
        il, ih = np.clip(i - 1, 0, len(K) - 1), np.clip(i, 0, len(K) - 1)
        dl = np.where(ev_a[il] == ev_row, np.abs(ta[il] - tt), np.inf)
        dh = np.where(ev_a[ih] == ev_row, np.abs(ta[ih] - tt), np.inf)
        t_rf[ok] = np.where(dh < dl, ta[ih], ta[il])
    return code, t_rf


def process_file(job):
    """One result_clusters.tsv -> per (definition, tier) phi histograms and counters."""
    path, rp, run, split, P, offs, want, mask = job
    acc = new_acc()
    c = acc["c"]
    out = dict(run=run, split=split, path=path, acc=acc, error=None)
    try:
        with open(path) as fh:
            head = fh.readline().rstrip("\n").split("\t")
            if len(head) != NCOL:
                out["error"] = f"{len(head)} columns in the header (need {NCOL})"
                return out
            evs, st, nmem, seed_name, mem = [], [], [], [], []
            for line in fh:
                f = line.rstrip("\n").split("\t")
                if len(f) != NCOL:
                    c["bad_line"] += 1
                    continue
                if f[1] != want:
                    c["rows_other_path"] += 1
                    continue
                c["rows_read"] += 1
                try:
                    evs.append(int(f[12]))
                    st.append(float(f[13]))
                    nmem.append(int(f[14]))
                except ValueError:
                    c["rows_read"] -= 1
                    c["bad_line"] += 1
                    continue
                s_ = f[2]   # _mod_name writes IDs > 2152 as bare numbers (W1153-W1156)
                seed_name.append(f"W{int(s_) - 1000}" if s_.isdigit() and int(s_) > 2152 else s_)
                mem.append(f[15])
    except OSError as e:
        out["error"] = str(e)
        return out
    ev = np.asarray(evs, np.int64)
    st = np.asarray(st, np.float64)
    fin = np.isfinite(st)
    code, t_rf = match_rf(rp, ev, np.where(fin, st, 0.0), mask)
    c["rf_event_missing"] += int((code == 1).sum())
    c["trig_rejected"] += int((code == 2).sum())
    c["no_a_hit"] += int((code == 3).sum())
    for k in np.flatnonzero(code == 0):
        if not fin[k]:
            c["seed_time_nonfinite"] += 1
            continue
        if not mem[k]:
            c["no_members"] += 1
            continue
        o_s = offs.get(seed_name[k])
        if o_s is None:
            c["seed_no_offset"] += 1
            continue
        items = mem[k].split(";")
        n_col = nmem[k]
        if len(items) != n_col:
            c["n_members_mismatch"] += 1
        ts = st[k] - o_s
        sw = su = 0.0
        n_used = 0
        seed_in = False
        n_nooff = n_nonfin = 0
        for it in items:
            nm, e, t = it.split(":")
            e, t = float(e), float(t)
            o = offs.get(nm)
            if o is None:
                n_nooff += 1
                continue
            if not (math.isfinite(t) and math.isfinite(e) and e > 0):
                n_nonfin += 1
                continue
            u = t - o
            u = ts + float(fold(u - ts, P))
            sw += e
            su += e * u
            n_used += 1
            if nm == seed_name[k]:
                seed_in = True
        c["members_total"] += len(items)
        c["members_no_offset"] += n_nooff
        c["members_nonfinite"] += n_nonfin
        if not seed_in:
            c["seed_not_member"] += 1
            continue
        if not sw > 0:
            c["ewt_nonfinite"] += 1
            continue
        c["rows_used"] += 1
        c["members_used"] += n_used
        phi_s = float(fold(ts - t_rf[k], P))
        phi_e = float(fold(su / sw - t_rf[k], P))
        bs, be = int(_bin(np.asarray(phi_s), P)), int(_bin(np.asarray(phi_e), P))
        for ti in (0, tier_of(n_col)):
            acc["H"][0, ti, bs] += 1
            acc["H"][1, ti, be] += 1
            acc["nm"][ti, min(n_col, NMEM_MAX)] += 1
            acc["mem_used"][ti] += n_used
            acc["mem_tot"][ti] += len(items)
    return out


# ---- aggregation ------------------------------------------------------------

def add_acc(a, b):
    a["H"] += b["H"]
    a["nm"] += b["nm"]
    a["mem_used"] += b["mem_used"]
    a["mem_tot"] += b["mem_tot"]
    for k in COUNTERS:
        a["c"][k] += b["c"][k]


def hist_median(h):
    n = h.sum()
    if not n:
        return float("nan")
    cs = np.cumsum(h)
    return float(np.searchsorted(cs, 0.5 * n))


def fit_all(items, min_n):
    """items: list of (key, hist, P) -> dict key -> fit result (None when N < min_n)."""
    out = {}
    for key, h, P in items:
        out[key] = fit_unwrapped((h.astype(float), P)) if h.sum() >= min_n else None
    return out


def stat_row(h, P, fit, nm, mem_used, mem_tot):
    N = int(h.sum())
    if fit is None:
        mu = RF.unwrap(h, P)[1] if N else float("nan")
        sg = se = float("nan")
        ok = 0
    else:
        mu = fit["mu_u"]
        ok = int(bool(fit["ok"]))
        sg = fit["sigma"] if ok else float("nan")
        se = fit["sigma_err"] if ok else float("nan")
    med = hist_median(nm)
    frac = mem_used / mem_tot if mem_tot else float("nan")
    return N, mu, sg, se, ok, med, frac


def fmt(x, nd=5):
    return "" if x is None or (isinstance(x, float) and not np.isfinite(x)) else (
        f"{x:.{nd}f}" if isinstance(x, float) else str(x))


HEAD = ["run", "type", "path", "definition", "tier", "N", "P_ns", "mu_ns", "sigma_ns", "sigma_err_ns", "fit_ok",
        "median_n_members", "frac_members_with_offset"]


def table_rows(prefix, typ, path, acc, P, fits):
    rows = []
    for di, dn in enumerate(DEFS):
        for ti, tn in enumerate(TIERS):
            N, mu, sg, se, ok, med, frac = stat_row(acc["H"][di, ti], P, fits[(di, ti)], acc["nm"][ti],
                                                    acc["mem_used"][ti], acc["mem_tot"][ti])
            rows.append([*prefix, typ, path, dn, tn, str(N), fmt(P, 6), fmt(mu), fmt(sg), fmt(se), str(ok),
                         fmt(med, 1), fmt(frac, 4)])
    return rows


# ---- PDF --------------------------------------------------------------------

def render_pdf(path, pages, min_n):
    from matplotlib.backends.backend_pdf import PdfPages
    plt = PMT._mpl()
    with PdfPages(path) as pdf:
        for tag, acc, P, fits, runs in pages:
            fig, axs = plt.subplots(2, 4, figsize=(14, 7.5), sharey="row")
            fig.suptitle(f"{LABEL[tag]}: cluster time - t_RF, offset-corrected, folded (P = {P:.5f} ns); runs "
                         + " ".join(str(r) for r in runs), fontsize=11, weight="bold")
            w = P / NB
            for di, dn in enumerate(DEFS):
                for ti, tn in enumerate(TIERS):
                    ax = axs[di, ti]
                    h = acc["H"][di, ti]
                    N = int(h.sum())
                    ax.set_title(f"{dn}, members {tn}", fontsize=10)
                    ax.tick_params(labelsize=8)
                    if N == 0:
                        ax.text(0.5, 0.5, "no data", transform=ax.transAxes, ha="center", va="center", fontsize=9)
                        continue
                    a, cen, x, cnt = RF.unwrap(h, P)
                    ax.stairs(cnt, np.r_[x - w / 2, x[-1] + w / 2], color="C0", lw=0.9)
                    fit = fits[(di, ti)]
                    txt = [f"N = {N}"]
                    if fit is not None and fit["ok"]:
                        xg = np.linspace(x[0], x[-1], 500)
                        g = lambda xx: fit["amp"] * np.exp(-0.5 * ((xx - fit["mu"]) / fit["sigma"]) ** 2)
                        ax.plot(xg, g(xg), color="crimson", lw=0.8, ls="--", alpha=0.6)
                        xw = np.linspace(fit["lo"], fit["hi"], 200)
                        ax.plot(xw, g(xw), color="crimson", lw=1.8)
                        txt += [f"mean = {fit['mu']:.3f} ns", f"sigma = {fit['sigma']:.3f} +- {fit['sigma_err']:.3f} ns"]
                    else:
                        txt += ["not fitted" if fit is None else "fit failed"]
                    med = hist_median(acc["nm"][ti])
                    txt += [f"median members {med:.0f}",
                            f"frac with offset {acc['mem_used'][ti] / max(acc['mem_tot'][ti], 1):.3f}"]
                    ax.text(0.97, 0.96, "\n".join(txt), transform=ax.transAxes, ha="right", va="top", fontsize=7.5,
                            bbox=dict(boxstyle="round", fc="white", ec="0.7", alpha=0.85))
                    ax.set_ylim(bottom=0)
                    ax.set_xlabel("t - t_RF (ns, folded)", fontsize=8)
                    if ti == 0:
                        ax.set_ylabel(f"entries / {w:.3f} ns", fontsize=8)
            fig.tight_layout(rect=(0, 0, 1, 0.95))
            pdf.savefig(fig)
            plt.close(fig)


# ---- main -------------------------------------------------------------------

def parse_args(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--scan-dir", required=True, help="DIR/scan_<run6>/split_<n>/result_clusters.tsv")
    ap.add_argument("--rf-dir", required=True, help="rf_time_dump.py npz: <run6>_<split5>.npz")
    ap.add_argument("--offsets-dir", required=True, help="rf_time_calib.py output: <run5>.json, rf_calib_runs.tsv")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--path", choices=("legacy", "gated"), default="legacy", help="cluster path (config_name)")
    ap.add_argument("--jobs", type=int, default=NMAX_JOBS)
    ap.add_argument("--min-n", type=int, default=300, help="histograms with fewer entries are not fitted")
    ap.add_argument("--trigger-mask", type=lambda s: int(s, 0), default=0xf00,
                    help="RF.prep event selection: (trig & mask) != 0 (default 0xf00)")
    return ap.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    logf = open(out_dir / "run.log", "w")

    def log(msg=""):
        print(msg, flush=True)
        logf.write(msg + "\n")
        logf.flush()

    t0 = time.time()
    log(f"cluster_timing: scan {args.scan_dir}  rf {args.rf_dir}  offsets {args.offsets_dir}  path {args.path}  "
        f"min-n {args.min_n}  trigger-mask {args.trigger_mask:#x}")
    # run table
    runs_tsv = Path(args.offsets_dir) / "rf_calib_runs.tsv"
    run_tab = {}
    with open(runs_tsv) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            run_tab[int(r["run"])] = (float(r["P_ns"]), int(float(r["usable"])) != 0, r.get("quality", ""))
    # files
    files = sorted(glob.glob(os.path.join(args.scan_dir, "scan_*", "split_*", "result_clusters.tsv")))
    log(f"{len(files)} result_clusters.tsv found")
    jobs, runinfo, skipped = [], {}, {}
    for fp in files:
        m = re.search(r"scan_0*(\d+)[/\\]split_0*(\d+)[/\\]", fp)
        if not m:
            log(f"skip (name not scan_<run>/split_<n>): {fp}")
            continue
        run, split = int(m.group(1)), int(m.group(2))
        if run not in runinfo and run not in skipped:
            jp = Path(args.offsets_dir) / f"{run}.json"
            if run not in run_tab or not jp.is_file():
                skipped[run] = "no rf_calib_runs.tsv row or offsets json"
            elif run not in TYPE_OF_RUN:
                skipped[run] = "run not in plot_module_timing.RUN_TYPES"
            else:
                j = json.load(open(jp))
                P, usable_tsv, qual = run_tab[run]
                if not (usable_tsv and j.get("usable", False)):
                    skipped[run] = f"RF calibration not usable (quality {qual or j.get('quality')})"
                else:
                    runinfo[run] = (P, {m_["name"]: float(m_["offset_ns"]) for m_ in j["modules"]})
        if run in skipped:
            continue
        rp = Path(args.rf_dir) / f"{run:06d}_{split:05d}.npz"
        if not rp.is_file():
            log(f"skip {fp}: no RF file {rp.name}")
            continue
        P, offs = runinfo[run]
        jobs.append((fp, str(rp), run, split, P, offs, args.path, args.trigger_mask))
    for r, why in sorted(skipped.items()):
        log(f"run {r} skipped: {why}")
    log(f"{len(jobs)} files in {len(runinfo)} runs to process, {args.jobs} jobs")
    if not jobs:
        log("nothing to do")
        return 1

    if args.jobs > 1:
        with mp.Pool(min(args.jobs, len(jobs))) as pool:
            res = []
            for k, r in enumerate(pool.imap_unordered(process_file, jobs), 1):
                res.append(r)
                if k % 20 == 0 or k == len(jobs):
                    log(f"  {k}/{len(jobs)} files, {time.time() - t0:.0f} s")
    else:
        res = [process_file(j) for j in jobs]
    res.sort(key=lambda r: (r["run"], r["split"]))
    bad = [r for r in res if r["error"]]
    for r in bad:
        log(f"skip {r['path']}: {r['error']}")
    res = [r for r in res if not r["error"]]

    # aggregate
    per_run, per_type, tot = {}, {}, new_acc()
    split_items = []
    for r in res:
        a = r["acc"]
        per_run.setdefault(r["run"], new_acc())
        add_acc(per_run[r["run"]], a)
        add_acc(tot, a)
        for di in range(2):
            for ti in range(4):
                split_items.append(((r["run"], r["split"], di, ti), a["H"][di, ti], runinfo[r["run"]][0]))
    type_runs = {}
    for run in sorted(per_run):
        type_runs.setdefault(TYPE_OF_RUN[run], []).append(run)
    for tag, rs in type_runs.items():
        acc = new_acc()
        for run in rs:
            add_acc(acc, per_run[run])
        nn = np.array([per_run[run]["H"][0, 0].sum() for run in rs], float)
        Pt = float(np.average([runinfo[run][0] for run in rs], weights=nn if nn.sum() else None))
        per_type[tag] = (acc, Pt, rs)

    # fits
    fits_split = fit_all(split_items, args.min_n)
    fits_run = {run: fit_all([((di, ti), per_run[run]["H"][di, ti], runinfo[run][0])
                              for di in range(2) for ti in range(4)], args.min_n) for run in per_run}
    fits_type = {tag: fit_all([((di, ti), per_type[tag][0]["H"][di, ti], per_type[tag][1])
                               for di in range(2) for ti in range(4)], args.min_n) for tag in per_type}

    rows = []
    for run in sorted(per_run):
        rows += table_rows([str(run)], TYPE_OF_RUN[run], args.path, per_run[run], runinfo[run][0], fits_run[run])
    for tag in TAGS:
        if tag in per_type:
            acc, Pt, rs = per_type[tag]
            rows += table_rows(["ALL"], tag, args.path, acc, Pt, fits_type[tag])
    with open(out_dir / "cluster_timing.tsv", "w") as fh:
        fh.write("\t".join(HEAD) + "\n")
        for r in rows:
            fh.write("\t".join(r) + "\n")
    with open(out_dir / "cluster_timing_splits.tsv", "w") as fh:
        fh.write("\t".join([HEAD[0], "split", *HEAD[1:]]) + "\n")
        for r in res:
            run, a = r["run"], r["acc"]
            fits = {(di, ti): fits_split[(run, r["split"], di, ti)] for di in range(2) for ti in range(4)}
            for row in table_rows([str(run)], TYPE_OF_RUN[run], args.path, a, runinfo[run][0], fits):
                fh.write("\t".join([row[0], str(r["split"]), *row[1:]]) + "\n")

    try:
        render_pdf(str(out_dir / "cluster_timing.pdf"),
                   [(tag, per_type[tag][0], per_type[tag][1], fits_type[tag], per_type[tag][2])
                    for tag in TAGS if tag in per_type], args.min_n)
    except Exception as e:  # the TSVs are the product; report a failed PDF
        log(f"WARNING: PDF not written: {e!r}")

    # log
    log("")
    log(f"{'run':>6} {'files':>5} " + " ".join(f"{k[:11]:>11}" for k in ("rows_read", "rf_event_missing", "no_a_hit",
                                                                         "trig_rejected", "seed_no_offset",
                                                                         "rows_used")))
    nf = {}
    for r in res:
        nf[r["run"]] = nf.get(r["run"], 0) + 1
    for run in sorted(per_run):
        c = per_run[run]["c"]
        log(f"{run:>6} {nf[run]:>5} " + " ".join(f"{c[k]:>11}" for k in ("rows_read", "rf_event_missing", "no_a_hit",
                                                                          "trig_rejected", "seed_no_offset",
                                                                          "rows_used")))
    log("")
    log("totals:")
    for k in COUNTERS:
        log(f"  {k:<22} {tot['c'][k]}")
    mt, mu_, mn = tot["c"]["members_total"], tot["c"]["members_used"], tot["c"]["members_no_offset"]
    log(f"  members dropped (no offset) {mn}, (non-finite t/E) {tot['c']['members_nonfinite']} of {mt} "
        f"(members of rows that reached the member loop); used {mu_}")
    log(f"  frac_members_with_offset (rows used) = {tot['mem_used'][0] / max(tot['mem_tot'][0], 1):.4f}")
    log("")
    log("pooled per type (tier all):")
    for tag in TAGS:
        if tag in per_type:
            for di, dn in enumerate(DEFS):
                N, mu, sg, se, ok, med, frac = stat_row(per_type[tag][0]["H"][di, 0], per_type[tag][1],
                                                        fits_type[tag][(di, 0)], per_type[tag][0]["nm"][0],
                                                        per_type[tag][0]["mem_used"][0], per_type[tag][0]["mem_tot"][0])
                log(f"  {tag:<5} {dn:<5} N {N:>8}  mu {mu:+.3f}  sigma {sg:.3f} +- {se:.3f}  fit_ok {ok}  "
                    f"median n_members {med:.0f}  frac {frac:.3f}")
    log(f"done in {time.time() - t0:.0f} s: {out_dir / 'cluster_timing.tsv'}")
    logf.close()
    return 0


if __name__ == "__main__":
    sys.exit(main())
