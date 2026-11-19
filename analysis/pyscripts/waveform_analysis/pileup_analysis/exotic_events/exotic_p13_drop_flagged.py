#!/usr/bin/env python3
"""exotic_p13_drop_flagged.py -- P13: the exotic-events impact studies P7, P8, P9 and P10 recomputed under a channel-level drop of f1-f5.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p13_drop_flagged.py [--stage hooks|run|analyse|report|all]    [--jobs 2]
                                                                                        [--root <dir>] [--out <dir>] [--nboot 200] [--nproc 12]

Purpose, inputs, outputs.  Policy under test (user decision 2026-11-18): a kept PbWO4 channel in any exclusive category f1-f5 (exotic_p4_module_rates.flags_of,
uniform P2 thresholds, f4 > f5 > f1 > f2 > f3, through exotic_p7_cluster_level.row_categories) is not used: every benchmark hit of that (run, EVIO event number,
module) is removed from the hit list before clustering and before the gain solve.  No cluster veto, no correction.  Data scope as P7-P10: H2 run types on the
2-file exotic sample (output/pileup_rates/exotic_sample), Ta X17 on the 19 full-scan splits (output/pileup_rates/exotic_full), run 24937 excluded, PbWO4 only,
0.73 GeV calibration sample at seed radius >= 107 mm.  Outputs (--out): exotic_P13_drop_flagged_20261118.pdf and .tsv; intermediates under --root
(default output/pileup_rates/exotic_p13: hooks, the arm reports, logs; arm A is a symlink to P9's A reports, read only) and the recluster cache
output/pileup_rates/exotic_p13_cache_v2.npz (exotic_p13_recluster.py).

Arms (names as P9):
  A   default control v3 (P9's A reports, reused)
  D   drop: the flagged channels' hits removed before clustering in both benchmark passes (gain solve and evaluation; ep and the Ta Moller calibration)
  D'  ("Dp" in files) drop in the first pass only (the gain-solve input; ep only), clustering of the evaluated clusters unchanged: separates the calibration effect
V3 (P9) is quoted for reference from exotic_P9_calibration_impact_20261118.tsv.

Procedure.  Stage hooks writes one EXOTIC_HOOK npz per run type and arm (splits, empty veto / corr, rmin, drop [n, 3] = (run, event number, module) of every
flagged PbWO4 row of the in-scope exotic-scan files, drop_calib_only for D').  Stage run calls energy_resolution_benchmark.py report (ep chain) per run type
and arm and report --calib moller for D (Ta X17), as P9.  Stage analyse runs exotic_p13_recluster.py (P7's join and clustering with and without the
dropped hits: cluster fates, energy removed, elastic-window yield) and then: the benchmark's own arm values, paired pooled differences D - A and
D' - A on the clusters present in both arms (same run, split, event number, seed; bootstrap errors, P9's fit code), the clusters of A that are lost or
re-seeded under D (not paired; count and energy), per-module gain changes against sigma_g of A, the Ta block H, and the Moller-pair yield (pairs of A
missing in D).  Results: the PDF sections are named after the study they recompute (P7, P8, P9, P10), then the gate page.
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import math
import argparse
import multiprocessing as mp
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import exotic_p7_cluster_level as p7                    # noqa: E402
import exotic_p9_calibration_impact as p9               # noqa: E402
from exotic_p1_quantities import RT_LABEL, EXCLUDE_RUNS # noqa: E402
from exotic_p3_sample import text_page                  # noqa: E402

STEM = "exotic_P13_drop_flagged_20261118"
RTS, TA = p9.RTS, p9.TA
EP_ARMS = ("A", "D", "Dp")
EP_CMP = ("D", "Dp")
M_ARMS = ("A", "D")
M_CMP = ("D",)
ARM_LABEL = {"A": "A default", "D": "D drop", "Dp": "D' drop calib", "V3": "V3 veto"}
SCAN_ROOT = p9.SCAN_ROOT
P9_ROOT = SCAN_ROOT / "exotic_p9"
P9_TSV = HERE / "exotic_P9_calibration_impact_20261118.tsv"
RADBINS = p9.RADBINS
Z_SIG = p9.Z_SIG
fnum, one, get, dz, pm, tsv, write_rows = p9.fnum, p9.one, p9.get, p9.dz, p9.pm, p9.tsv, p9.write_rows


def wilson(k, n, z=1.0):
    """(fraction, lo, hi) Wilson interval (68 % at z = 1)."""
    if n <= 0:
        return np.nan, np.nan, np.nan
    p = k / n
    d = 1.0 + z * z / n
    c = (p + z * z / (2 * n)) / d
    h = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return p, max(c - h, 0.0), min(c + h, 1.0)


# ---------------------------------------------------------------------------------------------------------------- stage hooks
def drop_keys():
    """{run type: [n, 3] (run, event number, module)} of every flagged PbWO4 row of the in-scope exotic-scan files, and per-category row counts."""
    files = [("sample", f) for f in sorted((SCAN_ROOT / "exotic_sample").glob("run_*.npz"))]
    files += [("full", f) for f in sorted((SCAN_ROOT / "exotic_full").glob("run_*_s*.npz"))]
    out = {rt: [] for rt in RTS}
    cnt = {rt: {c: 0 for c in (1, 2, 3, 4, 5)} for rt in RTS}
    for src, f in files:
        with np.load(f) as z:
            run, rt = int(z["run"]), str(z["run_type"])
            if run in EXCLUDE_RUNS or not p9.in_scope(rt, src):
                continue
            key, cat, _, _ = p7.row_categories(z)
        out[rt].append(np.stack([np.full(key.size, run, np.int64), key // 4096, key % 4096], 1))
        for c in cnt[rt]:
            cnt[rt][c] += int((cat == c).sum())
    keys = {rt: np.unique(np.concatenate(v), axis=0) if v else np.zeros((0, 3), np.int64) for rt, v in out.items()}
    return keys, cnt


def stage_hooks(root):
    sp = p9.splits_of()
    dk, cnt = drop_keys()
    (Path(root) / "hooks").mkdir(parents=True, exist_ok=True)
    rows = []
    for rt in RTS:
        for arm in ("D", "Dp"):
            np.savez(p9.hook_path(root, rt, arm), splits=np.array(sp[rt], np.int64).reshape(-1, 2), veto=np.zeros((0, 3), np.int64),
                     rmin=np.float64(p9.RMIN.get(rt, 0.0)), corr=np.zeros((0, 4), np.int64), corr_rho=np.zeros(0), drop=dk[rt],
                     drop_calib_only=np.bool_(arm == "Dp"))
        rows.append(dict(table="hook", section="P9", run_type=rt, quantity="splits", value=len(sp[rt])))
        rows.append(dict(table="hook", section="P9", run_type=rt, quantity="drop_channels", value=len(dk[rt])))
        rows += [dict(table="hook", section="P9", run_type=rt, quantity=f"drop_rows_f{c}", value=v) for c, v in cnt[rt].items()]
        print(f"[hooks] {rt}: {len(sp[rt])} splits, {len(dk[rt])} dropped channels {cnt[rt]}", flush=True)
    write_rows(Path(root) / "hooks" / "hooks.tsv", rows)


# ---------------------------------------------------------------------------------------------------------------- stage run
def link_a(root):
    """<root>/A -> P9's A reports (read only, never rerun here)."""
    a = Path(root) / "A"
    if not a.exists():
        a.symlink_to(P9_ROOT / "A", target_is_directory=True)


def run_cmds(root):
    py, bench = sys.executable, str(p9.WA / "energy_resolution_benchmark.py")
    cmds = []
    for arm in EP_CMP:
        for rt in RTS:
            hk, out = p9.hook_path(root, rt, arm), p9.arm_root(root, rt, arm)
            cmds.append((f"{rt}_{arm}", [py, bench, "report", "--set", f"EXOTIC_HOOK={str(hk)!r}", "--out", str(out), str(p9.BENCH_SCAN / rt)],
                         out / rt / "p9_clusters.npz", out))
    for arm in M_CMP:
        hk, out = p9.hook_path(root, TA, arm), Path(root) / arm / "moller"
        cmds.append((f"moller_{arm}", [py, bench, "report", "--calib", "moller", "--set", f"EXOTIC_HOOK={str(hk)!r}", "--out", str(out), "--cache",
                                       str(out / "cache"), str(p9.BENCH_SCAN / TA)], out / f"{TA}_moller" / "p9_moller.npz", out))
    return cmds


def stage_run(root, jobs):
    link_a(root)
    cmds = [c for c in run_cmds(root) if not c[2].exists()]
    cmds.sort(key=lambda c: (0 if (TA in c[0] or "moller" in c[0]) else 1))
    print(f"[run] {len(cmds)} reports, {jobs} at a time", flush=True)
    with mp.get_context("spawn").Pool(jobs) as pool:
        for name, rc, info in pool.imap_unordered(p9.run_one, cmds):
            print(f"[run] {name}: rc {rc} ({info})", flush=True)
            if rc:
                raise SystemExit(f"[ERROR] {name} failed, see {info}")


# ---------------------------------------------------------------------------------------------------------------- stage analyse (P9, P10 Moller)
class _Arms:
    """Temporarily point P9's arm constants at P13's arms (P9's helpers read them at call time)."""
    def __init__(self, arms, cmp):
        self.new = (tuple(arms), tuple(cmp))

    def __enter__(self):
        self.old = (p9.ARMS, p9.CMP)
        p9.ARMS, p9.CMP = self.new

    def __exit__(self, *a):
        p9.ARMS, p9.CMP = self.old


def bench_rows(root):
    """The benchmark's own arm values (mean +- SE over the splits, control-gate path) and per-split paired differences, as P9."""
    rows = []
    for rt in RTS:
        gp = p9.gate_path(rt)
        P = {}
        for arm in EP_ARMS:
            s = next(r for r in tsv(p9.ep_dir(root, rt, arm) / "summary.tsv") if r["path"] == gp)
            P[arm] = {(r["run"], r["split"]): r for r in tsv(p9.ep_dir(root, rt, arm) / "splits.tsv") if r["path"] == gp and r["variant"] == "all"}
            for q in ("res_cb", "mu_ep", "tail", "res_cb_oos", "mu_ep_oos"):
                rows.append(dict(table="bench", section="P9", run_type=rt, calib="ep", arm=arm, quantity=q, value=fnum(s[q]), se=fnum(s[f"{q}_se"]),
                                 n=fnum(s["n_used"]), n_splits=fnum(s["n_splits"])))
        for c in EP_CMP:
            for q in ("res_cb", "mu_ep", "tail", "res_cb_oos", "mu_ep_oos"):
                d = np.array([fnum(P[c][k][q]) - fnum(P["A"][k][q]) for k in P["A"] if k in P[c]])
                d = d[np.isfinite(d)]
                rows.append(dict(table="bench_paired", section="P9", run_type=rt, calib="ep", comparison=f"{c}-A", quantity=q,
                                 value=float(d.mean()) if d.size else np.nan, se=float(d.std(ddof=1) / math.sqrt(d.size)) if d.size > 1 else np.nan, n_splits=int(d.size)))
    for arm in M_ARMS:
        for r in tsv(Path(root) / arm / "moller" / f"{TA}_moller" / "summary.tsv"):
            if r["path"] == p9.gate_path(TA) and r["quantity"] in ("sum", "cluster"):
                for q in ("res", "res_pooled", "peak"):
                    rows.append(dict(table="bench", section="P9", run_type=TA, calib="moller", arm=arm, protocol=r["protocol"], quantity=f"{r['quantity']}:{q}",
                                     value=fnum(r[q]), se=fnum(r.get(f"{q}_se", "nan")), n=fnum(r["n"]), n_splits=fnum(r["n_splits"])))
    return rows


def boot_rows(root, nboot, nproc):
    """Paired pooled differences X - A on the shared clusters (ep: D, D'; Moller: D), bootstrap errors (P9's fit and bootstrap jobs)."""
    with _Arms(EP_ARMS, EP_CMP):
        p9._G["ep"] = p9.prep_ep(root)
    with _Arms(M_ARMS, M_CMP):
        p9._G["m"] = p9.prep_moller(root)
    ep_jobs = [(rt, c, prot, "overall", -1, nboot) for rt in RTS for c in EP_CMP for prot in ("in", "oos")]
    ep_jobs += [(rt, c, "oos", "ring", b, max(nboot // 2, 50)) for rt in RTS for c in EP_CMP for b in range(len(RADBINS) - 1)]
    m_jobs = [(c, prot, kind, None, nboot) for c in M_CMP for prot in ("in", "oos") for kind in ("sum", "cluster")]
    m_jobs += [(c, "oos", "cluster", b, max(nboot // 2, 50)) for c in M_CMP for b in range(len(RADBINS) - 1)]
    with mp.get_context("fork").Pool(nproc) as pool:
        ep_res = pool.map(p9.ep_boot_job, ep_jobs, chunksize=1)
        m_res = pool.map(p9.m_boot_job, m_jobs, chunksize=1)
    rows = []
    for (rt, c, prot, axis, b, nb), r in ep_res:
        S = p9._G["ep"][rt, c, prot]
        base = dict(table="paired", section="P9", run_type=rt, calib="ep", comparison=f"{c}-A", protocol=prot, axis=axis,
                    bin="all" if axis == "overall" else f"{RADBINS[b]:g}-{RADBINS[b + 1]:g}", n_A=S["nA"], n_X=S["nX"])
        if r is None:
            rows.append(base | dict(quantity="n_shared", value=int((S["ring"] == b).sum())))
            continue
        for q in p9.EPQ:
            dv = r["X"][q] - r["A"][q]
            rows.append(base | dict(quantity=q, n_shared=r["A"]["n"], value_A=r["A"][q], value_X=r["X"][q], diff=dv, se=r["sd"][q],
                                    z=dv / r["sd"][q] if r["sd"][q] > 0 else np.nan, nboot=r["nboot"]))
    for (c, prot, kind, b, nb), r in m_res:
        S = p9._G["m"][c]
        base = dict(table="paired", section="P9", run_type=TA, calib="moller", comparison=f"{c}-A", protocol=prot, axis="overall" if b is None else "ring",
                    bin=(kind if b is None else f"{RADBINS[b]:g}-{RADBINS[b + 1]:g}"), n_A=S["nA"], n_X=S["nX"])
        if r is None:
            rows.append(base | dict(quantity="n_shared", value=0))
            continue
        for q in p9.MQ:
            dv = r["X"][q] - r["A"][q]
            rows.append(base | dict(quantity=q, n_shared=r["A"]["n"], value_A=r["A"][q], value_X=r["X"][q], diff=dv, se=r["sd"][q],
                                    z=dv / r["sd"][q] if r["sd"][q] > 0 else np.nan, nboot=r["nboot"]))
    return rows


def boot_frac(num, den, nboot=400, seed=1):
    """sum(num) / sum(den) and its bootstrap sd over the clusters."""
    if den.size == 0 or den.sum() <= 0:
        return np.nan, np.nan
    rng = np.random.default_rng(seed)
    v = [num[k].sum() / den[k].sum() for k in (rng.integers(0, den.size, den.size) for _ in range(nboot))]
    return float(num.sum() / den.sum()), float(np.std(v, ddof=1))


def lost_rows(root):
    """Clusters of A (used in the pool, in and out of sample) absent from the D dump: re-seeded (another cluster of the same event is in D) or gone; and the
    clusters of D absent from A (new).  Not paired; counts with Wilson 68 %, energy share with a bootstrap error."""
    rows = []
    for rt in RTS:
        dA, dD = p9.load_dump(root, rt, "A"), p9.load_dump(root, rt, "D")
        kA, kD = p9.keys_of(dA), p9.keys_of(dD)
        setD = set(map(tuple, kD.tolist()))
        evD = set(map(tuple, kD[:, :3].tolist()))
        setA = set(map(tuple, kA.tolist()))
        for prot, eq, uq in (("in", "e", "used"), ("oos", "e_oos", "used_oos")):
            ua = np.flatnonzero(dA[uq])
            inD = np.array([tuple(k) in setD for k in kA[ua].tolist()], bool)
            ev = np.array([tuple(k[:3]) in evD for k in kA[ua].tolist()], bool)
            reseed, gone = ~inD & ev, ~inD & ~ev
            usedD = np.zeros(ua.size, bool)
            if inD.any():
                iD = {tuple(k): i for i, k in enumerate(kD.tolist())}
                usedD[inD] = [bool(dD[uq][iD[tuple(k)]]) for k in kA[ua][inD].tolist()]
            e = dA[eq][ua].astype(float)
            n = int(ua.size)
            uD = np.flatnonzero(dD[uq])
            new = np.array([tuple(k) not in setA for k in kD[uD].tolist()], bool)
            for what, m in (("lost", ~inD), ("reseeded", reseed), ("gone", gone), ("unused_in_D", inD & ~usedD)):
                f, lo, hi = wilson(int(m.sum()), n)
                ef, ese = boot_frac(e * m, e, seed=len(rows) + 1)
                rows.append(dict(table="lost", section="P9", run_type=rt, calib="ep", protocol=prot, quantity=what, k=int(m.sum()), n=n, frac=f, lo=lo, hi=hi,
                                 e_frac=ef, e_frac_se=ese, e_median=float(np.median(e[m])) if m.any() else np.nan))
            f, lo, hi = wilson(int(new.sum()), int(uD.size))
            rows.append(dict(table="lost", section="P9", run_type=rt, calib="ep", protocol=prot, quantity="new_in_D", k=int(new.sum()), n=int(uD.size), frac=f, lo=lo, hi=hi,
                             e_median=float(np.median(dD[eq][uD][new])) if new.any() else np.nan))
    return rows


def gain_rows(root, M):
    """Per (calibration, comparison, run, module): dg = g_X - g_A, sigma_g of A; free = prior fraction of A < 0.5; free_X = that of X."""
    idx = {n: i for i, n in enumerate(M["name"])}
    rows = []
    for rt, calib, arms in [(rt, "ep", EP_ARMS) for rt in RTS] + [(TA, "moller", M_ARMS)]:
        G = {arm: p9.gains_of(root, rt, arm, calib) for arm in arms}
        for c in arms[1:]:
            for k, (g, sg, ns, pf, tab) in G["A"].items():
                i = idx.get(k[1])
                if k not in G[c] or i is None or not M["pwo4"][i]:
                    continue
                gx, pfx = G[c][k][0], G[c][k][3]
                rows.append(dict(table="gain", section="P9", run_type=rt, calib=calib, comparison=f"{c}-A", run=k[0], module=k[1],
                                 ring=f"{RADBINS[M['ring'][i]]:g}-{RADBINS[M['ring'][i] + 1]:g}", region=p9.erb_mod().REGION_NAMES[M["reg"][i]] if M["reg"][i] >= 0 else "LG",
                                 pool=int(M["pool"][i]), solve=tab, g_A=g, sigma_g=sg, n_seed=ns, prior_frac=pf, prior_frac_X=pfx, free=int(pf < 0.5),
                                 free_X=int(pfx < 0.5), dg=gx - g, z=(gx - g) / sg if sg > 0 else np.nan))
    return rows


def moller_yield_rows(root, M):
    """Ta Moller pairs (selected, gate path) of A missing under D (lost) and new under D, overall and per ring of the inner member's seed radius; Wilson 68 %."""
    dA, dD = p9.load_moller(root, "A"), p9.load_moller(root, "D")
    pk = lambda d: np.stack([d["run"][0::2], d["split"][0::2], d["evn"][0::2]], 1).astype(np.int64)   # noqa: E731
    kA, kD = pk(dA), pk(dD)
    sD, sA = set(map(tuple, kD.tolist())), set(map(tuple, kA.tolist()))
    lost = np.array([tuple(k) not in sD for k in kA.tolist()], bool)
    new = np.array([tuple(k) not in sA for k in kD.tolist()], bool)
    rad = M["rad"]
    rin = np.minimum(rad[dA["seed"][0::2].astype(np.int64)], rad[dA["seed"][1::2].astype(np.int64)])
    ring = np.clip(np.digitize(rin, RADBINS) - 1, 0, len(RADBINS) - 2)
    rows = []
    for b in [-1] + list(range(len(RADBINS) - 1)):
        m = np.ones(lost.size, bool) if b < 0 else ring == b
        f, lo, hi = wilson(int((lost & m).sum()), int(m.sum()))
        rows.append(dict(table="yield_moller", section="P10", run_type=TA, axis="overall" if b < 0 else "ring",
                         bin="all" if b < 0 else f"{RADBINS[b]:g}-{RADBINS[b + 1]:g}", quantity="lost", k=int((lost & m).sum()), n=int(m.sum()), frac=f, lo=lo, hi=hi))
    f, lo, hi = wilson(int(new.sum()), int(new.size))
    rows.append(dict(table="yield_moller", section="P10", run_type=TA, axis="overall", bin="all", quantity="new_in_D", k=int(new.sum()), n=int(new.size), frac=f, lo=lo, hi=hi))
    return rows


def p9_reference_rows():
    """P9's V3 - A rows (paired, oos, overall; ep and the Moller sum) and its V3 - A gain extremes, for reference."""
    rows = []
    for r in tsv(P9_TSV):
        if r.get("table") == "paired" and r.get("comparison") == "V3-A" and r.get("axis") == "overall" and r.get("protocol") in ("oos", "in"):
            rows.append(dict(r) | dict(section="P9", source="P9"))
        if r.get("table") == "bench" and r.get("arm") == "V3":
            rows.append(dict(r) | dict(section="P9", source="P9"))
    return rows


def p10_reference():
    """{(dataset, population, veto): (k, n, frac, lo, hi)} of P10's overall veto table (V3, V4)."""
    out = {}
    for r in tsv(p9.P10_TSV):
        if r["table"] == "veto" and r["axis"] == "overall" and r["veto"] in ("V3", "V4"):
            out[r["dataset"], r["population"], r["veto"]] = (int(fnum(r["k"])), int(fnum(r["n"])), fnum(r["frac"]), fnum(r["lo"]), fnum(r["hi"]))
    return out


# ---------------------------------------------------------------------------------------------------------------- report: P9 section
def pct_ci(k, n, f, lo, hi, nd=3):
    """'f % [lo, hi] (k/n)' with f, lo, hi fractions."""
    f, lo, hi = fnum(f), fnum(lo), fnum(hi)
    if not math.isfinite(f):
        return "n/a"
    return f"{100 * f:.{nd}f} % [{100 * lo:.{nd}f}, {100 * hi:.{nd}f}] ({int(fnum(k))}/{int(fnum(n))})"


def table_lines(rows):
    L = ["Arm values: the benchmark's own numbers on these splits (control-gate path, mean +- SE over the splits; in = in sample, oos = gains / position fitted on the",
         "other half).  Differences: paired, pooled fit of the clusters present and used in both arms (same run, split, event number, seed module), bootstrap error;",
         "value +- error (z).  V3 - A from P9 for reference.  Full-data control v3 oos CB sigma/E: " + ", ".join(f"{RT_LABEL[rt]} {p9.BASELINE[rt]}" for rt in RTS) + " %.", ""]
    for rt in RTS:
        L.append(f"{RT_LABEL[rt]}  (splits {int(fnum(one(rows, table='hook', run_type=rt, quantity='splits').get('value')))}, dropped channels "
                 f"{int(fnum(one(rows, table='hook', run_type=rt, quantity='drop_channels').get('value')))}" + (f", seed radius >= {p9.RMIN[rt]:g} mm" if rt in p9.RMIN else "") + ")")
        L.append(f"  {'arm':<15}{'n used':>9}{'sigma/E in (%)':>20}{'sigma/E oos (%)':>20}{'mu in (%)':>20}{'mu oos (%)':>20}{'tail 3 sigma in (%)':>21}")
        for arm in EP_ARMS + ("V3",):
            b = lambda q: one(rows, table="bench", run_type=rt, calib="ep", arm=arm, quantity=q)   # noqa: E731
            if not b("res_cb"):
                continue
            L.append(f"  {ARM_LABEL[arm]:<15}{int(fnum(b('res_cb').get('n'))):>9}" + "".join(f"{pm(b(q).get('value'), b(q).get('se')):>20}"
                                                                                        for q in ("res_cb", "res_cb_oos", "mu_ep", "mu_ep_oos", "tail")))
        L.append(f"  {'paired diff':<15}{'n shared':>9}{'d sigma/E oos':>27}{'d mu oos':>27}{'d tail lo oos':>27}{'d tail hi oos':>27}{'d sigma/E in':>27}")
        for c in EP_CMP + ("V3",):
            g = lambda q, prot="oos": one(rows, table="paired", run_type=rt, calib="ep", comparison=f"{c}-A", protocol=prot, axis="overall", quantity=q)   # noqa: E731
            if not g("res_cb"):
                continue
            L.append(f"  {ARM_LABEL[c] + ' - A':<15}{int(fnum(g('res_cb').get('n_shared'))):>9}" + "".join(f"{dz(g(q)):>27}" for q in p9.EPQ)
                     + f"{dz(g('res_cb', 'in')):>27}")
        L.append("")
    L.append(f"{RT_LABEL[TA]}, Moller two-cluster calibration (pairs, gate path; sum = E1 + E2 vs Ebeam + m_e - losses, cluster = E / E_pred); D' not defined here")
    L.append(f"  {'arm':<15}{'pairs':>8}{'sum sigma/E oos (%)':>24}{'pooled':>22}{'peak oos (%)':>14}{'cluster sigma/E oos (%)':>26}")
    for arm in M_ARMS + ("V3",):
        b = lambda q: one(rows, table="bench", calib="moller", arm=arm, protocol="out-of-sample", quantity=q)   # noqa: E731
        if not b("sum:res"):
            continue
        L.append(f"  {ARM_LABEL[arm]:<15}{int(fnum(b('sum:res').get('n'))):>8}{pm(b('sum:res').get('value'), b('sum:res').get('se')):>24}"
                 f"{pm(b('sum:res_pooled').get('value'), b('sum:res_pooled').get('se')):>22}{fnum(b('sum:peak').get('value')):>14.3f}"
                 f"{pm(b('cluster:res').get('value'), b('cluster:res').get('se')):>26}")
    L.append(f"  {'paired diff':<15}{'pairs':>8}{'d sum sigma/E oos':>27}{'d sum peak oos':>27}{'d sum tail lo':>27}{'d sum tail hi':>27}{'d cluster sigma/E oos':>27}")
    for c in M_CMP + ("V3",):
        g = lambda q, kind="sum": one(rows, table="paired", calib="moller", comparison=f"{c}-A", protocol="oos", axis="overall", bin=kind, quantity=q)   # noqa: E731
        if not g("res"):
            continue
        L.append(f"  {ARM_LABEL[c] + ' - A':<15}{int(fnum(g('res').get('n_shared'))):>8}" + "".join(f"{dz(g(q)):>27}" for q in p9.MQ) + f"{dz(g('res', 'cluster')):>27}")
    L.append("  The Moller peak is quoted without a per-split error by the benchmark; its paired difference carries the bootstrap error.")
    return L


def lost_lines(rows):
    L = ["Clusters of A (used in the ep pool) that are not in the D dump (no cluster with the same run, split, event number and seed): 're-seeded' = the event still",
         "has a D pool cluster with another seed (its seed was dropped, a neighbour seeds it), 'gone' = no D pool cluster in that event.  'unused in D' = same key in",
         "D but outside D's used set.  These are not in the paired fits above.  Fractions Wilson 68 %; energy share = sum E of those / sum E of A used, bootstrap error.", ""]
    L.append(f"  {'run type':<18}{'protocol':<9}{'kind':<13}{'fraction of A used clusters':>40}{'energy share (%)':>22}{'median E (MeV)':>16}")
    for rt in RTS:
        for prot in ("in", "oos"):
            for q in ("lost", "reseeded", "gone", "unused_in_D", "new_in_D"):
                r = one(rows, table="lost", run_type=rt, protocol=prot, quantity=q)
                if not r:
                    continue
                es = pm(100 * fnum(r.get("e_frac")), 100 * fnum(r.get("e_frac_se")), 4) if q != "new_in_D" else "(of D used)"
                L.append(f"  {RT_LABEL[rt]:<18}{prot:<9}{q:<13}{pct_ci(r['k'], r['n'], r['frac'], r['lo'], r['hi'], 4):>40}{es:>22}{fnum(r.get('e_median')):>16.1f}")
        L.append("")
    return L


def ring_lines(rows):
    L = [f"Per P4 ring of the seed-module centre radius, out of sample, paired as above (rings with < {p9.MIN_RING} shared clusters omitted; ep pool reaches 200 mm;",
         "Moller: the cluster position radius).  Each entry: difference +- bootstrap error (z).", ""]
    for rt, calib, cmps in [(rt, "ep", EP_CMP) for rt in RTS] + [(TA, "moller", M_CMP)]:
        q1, q2 = ("res_cb", "mu_ep") if calib == "ep" else ("res", "peak")
        L.append(f"{RT_LABEL[rt]} ({calib})")
        L.append(f"  {'ring (mm)':<11}{'n':>8}{'sigma/E A (%)':>14}" + "".join(f"{'d sigma/E ' + c:>27}" for c in cmps) + f"{'mu A (%)':>10}" + "".join(f"{'d mu ' + c:>27}" for c in cmps))
        for b in range(len(RADBINS) - 1):
            lab = f"{RADBINS[b]:g}-{RADBINS[b + 1]:g}"
            g = lambda c, q: one(rows, table="paired", run_type=rt, calib=calib, comparison=f"{c}-A", protocol="oos", axis="ring", bin=lab, quantity=q)   # noqa: E731
            r0 = g("D", q1)
            if not r0 or "diff" not in r0:
                continue
            L.append(f"  {lab:<11}{int(fnum(r0.get('n_shared'))):>8}{fnum(r0.get('value_A')):>14.3f}" + "".join(f"{dz(g(c, q1)):>27}" for c in cmps)
                     + f"{fnum(g('D', q2).get('value_A')):>10.3f}" + "".join(f"{dz(g(c, q2)):>27}" for c in cmps))
        L.append("")
    return L


def gain_lines(rows, top=12):
    gs = p9.gain_summary(rows)
    L = ["Per-module gain change dg = g_X - g_A per run (ep: pool solve, ring-scoped solve where only it sets the module; Moller: per-run in-sample solve), in units of",
         "sigma_g of A.  Free = prior fraction of A < 0.5.  Counts of module-runs; 'beyond 2 sigma' is the gate criterion.", ""]
    L.append(f"  {'calibration':<26}{'comparison':<8}{'free module-runs':>18}{'|dg| > 1 sigma':>16}{'|dg| > 2 sigma':>16}   largest |dg| / sigma_g")
    for (rt, calib, c), o in sorted(gs.items(), key=lambda x: (RTS.index(x[0][0]), x[0][1], x[0][2])):
        w = o["worst"]
        L.append(f"  {RT_LABEL[rt] + ' ' + calib:<26}{c:<8}{o['n']:>18}{o['n1']:>16}{o['n2']:>16}   {w['module']} run {w['run']} dg {fnum(w['dg']):+.5f} +- {fnum(w['sigma_g']):.5f} "
                 f"({100 * fnum(w['dg']) / fnum(w['g_A']):+.3f} +- {100 * fnum(w['sigma_g']) / fnum(w['g_A']):.3f} %, {fnum(w['z']):+.2f} sigma)")
    L += ["", f"Largest |dg| / sigma_g (free modules, top {top}):"]
    fr = sorted([r for r in rows if r.get("table") == "gain" and int(fnum(r["free"]))], key=lambda r: -abs(fnum(r["z"])))[:top]
    L.append(f"  {'calibration':<26}{'cmp':<6}{'run':>7}{'module':>8}{'ring':>10}{'g_A':>10}{'sigma_g':>10}{'dg':>11}{'dg %':>18}{'z':>8}{'n_seed':>8}")
    for r in fr:
        L.append(f"  {RT_LABEL[r['run_type']] + ' ' + r['calib']:<26}{r['comparison']:<6}{r['run']:>7}{r['module']:>8}{r['ring']:>10}{fnum(r['g_A']):>10.5f}"
                 f"{fnum(r['sigma_g']):>10.5f}{fnum(r['dg']):>+11.5f}{pm(100 * fnum(r['dg']) / fnum(r['g_A']), 100 * fnum(r['sigma_g']) / fnum(r['g_A'])):>18}"
                 f"{fnum(r['z']):>+8.2f}{int(fnum(r['n_seed'])):>8}")
    L += ["", "Ta X17 central block H (the inner ring around the beam hole, where Ta f1 / f3 concentrate): free (prior fraction < 0.5) module-runs in A and in X, and dg / sigma_g"]
    for calib in ("ep", "moller"):
        for c in (EP_CMP if calib == "ep" else M_CMP):
            rr = get(rows, table="gain", run_type=TA, calib=calib, comparison=f"{c}-A", region="H")
            if not rr:
                L.append(f"  {calib:<7}{c}-A: no block H module in the solve")
                continue
            z = np.array([fnum(r["z"]) for r in rr])
            fa, fx = sum(int(fnum(r["free"])) for r in rr), sum(int(fnum(r["free_X"])) for r in rr)
            fz = np.array([fnum(r["z"]) for r in rr if int(fnum(r["free"]))])
            pfa = np.array([fnum(r["prior_frac"]) for r in rr])
            pfx = np.array([fnum(r["prior_frac_X"]) for r in rr])
            L.append(f"  {calib:<7}{c + '-A':<6} {len(rr)} module-runs; free in A {fa}, in X {fx}; prior fraction median A {np.median(pfa):.2f} -> X {np.median(pfx):.2f}; "
                     f"|z| max over all {np.nanmax(np.abs(z)):.2f}, over free {np.nanmax(np.abs(fz)) if fz.size else float('nan'):.2f}")
    return L


# ---------------------------------------------------------------------------------------------------------------- report: P10 Moller, gate, caveats
def moller_yield_lines(rows, p10):
    L = ["Ta X17 Moller pairs (moller_calib selection, gate path, the 19 splits) of A that are missing under D (lost), per ring of the inner member's seed radius;",
         "pairs selected under D but not under A (new).  Wilson 68 %.  P10's veto costs on the same pairs for comparison.", ""]
    L.append(f"  {'axis':<9}{'bin':<10}{'lost under D':>42}")
    for r in get(rows, table="yield_moller", quantity="lost"):
        if not int(fnum(r["n"])):
            continue
        L.append(f"  {r['axis']:<9}{r['bin']:<10}{pct_ci(r['k'], r['n'], r['frac'], r['lo'], r['hi'], 3):>42}")
    r = one(rows, table="yield_moller", quantity="new_in_D")
    L.append(f"  {'overall':<9}{'new':<10}{pct_ci(r.get('k'), r.get('n'), r.get('frac'), r.get('lo'), r.get('hi'), 3):>42}  (of D pairs)")
    L.append("")
    for v in ("V3", "V4"):
        c = p10.get((TA, "moller_pair", v))
        if c:
            L.append(f"  P10 {v} (veto): {pct_ci(c[0], c[1], c[2], c[3], c[4], 3)}  (P10 counts the full Moller selection of the P7 cache, not these 19-split arm runs)")
    return L


def gate_eval(rows):
    """{(rt, calib, cmp): flags and values}: does the arm move sigma/E, mu or a tail (oos, overall, |z| > Z_SIG) or any free gain beyond Z_SIG sigma_g?"""
    gs = p9.gain_summary(rows)
    out = {}
    for rt, calib, cmps in [(rt, "ep", EP_CMP) for rt in RTS] + [(TA, "moller", M_CMP)]:
        for c in cmps:
            if calib == "ep":
                g = lambda q: one(rows, table="paired", run_type=rt, calib="ep", comparison=f"{c}-A", protocol="oos", axis="overall", quantity=q)   # noqa: E731
                qs = ("res_cb", "mu_ep", "tail_lo3", "tail_hi3")
            else:
                g = lambda q: one(rows, table="paired", calib="moller", comparison=f"{c}-A", protocol="oos", axis="overall", bin="sum", quantity=q)   # noqa: E731
                qs = ("res", "peak", "tail_lo3", "tail_hi3")
            R = {k: g(q) for k, q in zip(("res", "mu", "tlo", "thi"), qs)}
            o = gs.get((rt, calib, f"{c}-A"), dict(n=0, n1=0, n2=0, worst=None))
            e = {k: dict(d=fnum(r.get("diff")), se=fnum(r.get("se")), z=fnum(r.get("z"))) for k, r in R.items()}
            out[rt, calib, c] = dict(e, flag=any(abs(v["z"]) > Z_SIG for v in e.values() if math.isfinite(v["z"])) or o["n2"] > 0, n2=o["n2"], n=o["n"], worst=o["worst"])
    return out


def gate_lines(rows, ep_cost):
    """ep_cost: {run type: (k, n, f, lo, hi)} elastic-window loss under D (P10 section, from exotic_p13_recluster)."""
    G = gate_eval(rows)
    L = [f"Does the arm move sigma/E, the peak or a 3-sigma tail (out of sample, paired, |z| > {Z_SIG:g}) or any free module gain by more than {Z_SIG:g} sigma_g?",
         "Entries: difference +- bootstrap error (percentage points for sigma/E and tails, % for the peak).  Cost = loss of the ep elastic window (P10 section, A's window,",
         "fixed gains; gross = A clusters lost, Wilson 68 %; net = lost minus clusters entering the window, Poisson error) or of the Moller pairs under D (Wilson 68 %).", ""]
    L.append(f"  {'calibration':<24}{'arm':<4}{'d sigma/E':>22}{'d mu':>22}{'d tail lo':>22}{'d tail hi':>22}{'gains > 2 sg':>14}  {'moves?':<7}{'yield cost (D): gross; net of entries':>56}")
    for (rt, calib, c), o in G.items():
        f = lambda k: f"{o[k]['d']:+.4f} +- {o[k]['se']:.4f}" if math.isfinite(o[k]["d"]) else "n/a"   # noqa: E731
        if calib == "moller":
            r = one(rows, table="yield_moller", quantity="lost", axis="overall")
            cost = pct_ci(r.get("k"), r.get("n"), r.get("frac"), r.get("lo"), r.get("hi"), 3)
        else:
            cc = ep_cost.get(rt)
            cost = pct_ci(*cc, nd=3) if cc else "n/a"
            r = one(rows, table="p10_elastic_loss", run_type=rt, axis="overall")
            if r and c == "D":
                kl, ke, n = fnum(r["n_lost"]), fnum(r["n_enter"]), fnum(r["n_A"])
                cost += f"; net {100 * (kl - ke) / n:.3f} +- {100 * math.sqrt(kl + ke) / n:.3f} %"
        L.append(f"  {RT_LABEL[rt] + ' ' + calib:<24}{c:<4}{f('res'):>22}{f('mu'):>22}{f('tlo'):>22}{f('thi'):>22}{o['n2']:>7} of {o['n']:<5}  {'yes' if o['flag'] else 'no':<7}"
                 f"{cost if c == 'D' else '(as D)':>56}")
    L += [""] + RECOMMEND
    return L


RECOMMEND = [   # judgement from the 2026-11-18 run (numbers: this page, the P7 / P8 / P9 / P10 sections; intervals Wilson 68 %)
    "Reading.  H2 2.24 / 3.49 GeV: D leaves sigma/E and the peak where they are (|z| <= 0.7); the only |z| > 2 entries are the 2.24 GeV low tail, which D lowers",
    "  (-0.017 +- 0.006 pp, z -2.9, the removed fakes), and D' at 2.24 GeV (peak +0.0014 +- 0.0006 %, z 2.3).  One of 2121 free module-runs moves by more than 2 sigma_g",
    "  (W522 run 25299 -0.53 +- 0.21 %).  D removes every fake that V3 removes (23 / 14 in the elastic window, P7 section) but also real clusters: elastic-window",
    "  loss 0.029 % [0.026, 0.031] at 2.24 and 0.030 % [0.026, 0.034] at 3.49 GeV, against V3 0.0044 % [0.0036, 0.0055] and 0.0067 % [0.0051, 0.0087] (P10).",
    "H2 0.73 GeV: no change of sigma/E or the peak (|z| <= 0.5); two module-runs move beyond 2 sigma_g (W702 run 24603 -1.47 +- 0.34 %, W556 -0.84 +- 0.32 %);",
    "  elastic-window loss 0.100 % [0.083, 0.120], against V4 0.036 % [0.026, 0.049] and V3 0 % [0, 0.0036] (P10).",
    "Ta X17: sigma/E and the peak do not move (ep d sigma/E +0.010 +- 0.012 pp, d mu +0.005 +- 0.010 %), but the low tail grows (+0.24 +- 0.06 pp, z 3.8): the",
    "  dropped f1 / f3 neighbours take 0.40 +- 0.01 % of the elastic energy (P8 section, S first order) and the gain refit only recentres the peak (D' alone",
    "  moves it by +0.052 +- 0.006 %).  Moller: pair-sum peak +0.028 +- 0.011 % (z 2.6), sigma/E -0.014 +- 0.012 pp; W659 run 26230 -1.08 +- 0.44 % (-2.5 sigma_g).",
    "  Cost at fixed gains: 1.31 % [1.28, 1.34] of the elastic window gross, 0.87 +- 0.04 % net of the 570 clusters that enter it (P10 section); Moller pairs",
    "  0.24 % [0.21, 0.29] lost, 0.34 % [0.29, 0.39] of the D pairs new.  For comparison V4 costs 0.286 % [0.272, 0.302] of the elastic window and",
    "  0.122 % [0.096, 0.154] of the Moller pairs (P10).  Block H stays prior-held in the Moller solve (0 of 64 module-runs free in A and D) and largely in ep",
    "  (4 free in A, 3 in D, |dg| <= 0.56 sigma_g): D does not make the inner ring measurable.",
    "",
    "Recommendation: keep the current gate (V3 on 2.24 / 3.49 GeV H2 only; no correction; flags stay event markers).  Channel-level drop buys no measurable",
    "  resolution or peak improvement anywhere.  On H2 it loses 0.029 % [0.026, 0.031] / 0.030 % [0.026, 0.034] of the elastic window where V3 loses",
    "  0.0044 % [0.0036, 0.0055] / 0.0067 % [0.0051, 0.0087] for the same fakes; at 0.73 GeV it loses 0.100 % [0.083, 0.120] and moves two gains beyond 2 sigma_g;",
    "  on Ta it widens the low tail (z 3.8), shifts the Moller peak (z 2.6) and loses 0.87 +- 0.04 % (net, fixed gains) of the elastic window and",
    "  0.24 % [0.21, 0.29] of the Moller pairs.",
]
CAVEATS = [
    "Scope: the P7-P10 scope (H2 on the 2-file sample: 4 / 10 / 4 splits; Ta X17 on 19 splits of 25800 and 26230).  The arms are refit on these splits only, so",
    "  their values sit above the full-data control v3 numbers; the comparison is between arms on identical inputs.  Every benchmark event of these splits is in",
    "  the exotic scan (coverage 100 %), so every flagged channel of the evaluated events is known.",
    "Drop rule: all hits of a flagged (run, event, module) are removed, including pulses other than the flagged one (a Ta f1 row is often flagged on an early",
    "  pulse, P10).  The benchmark's below-threshold fill (BELOW_THR) still sees the dropped hits as stored, so a dropped module gets no virtual member.",
    "Re-seeding: a dropped seed lets a former member seed a smaller cluster (median E_D / E_A per sample with its bootstrap error on the P7 page), or the cluster",
    "  vanishes, or its remaining members join a neighbouring cluster ('merged'; each D cluster is matched to at most one A cluster).  Such clusters change their",
    "  key and are not paired (P9 'lost' page).  Events whose fake second cluster vanishes can enter the single-cluster selection ('new in D').",
    "S2 (P8 section) sums each matched D cluster once; vanished, merged and unmatched clusters add 0, so S2 can exceed S1 in size (Ta all clusters).",
    "D' (ep only): the gain solve (first pass, in and out of sample) sees the dropped hits removed, the evaluated clusters do not; its gains equal D's, so the",
    "  D' - A gain rows repeat D - A.  The Moller calibration has no D' (one candidate set serves the solve and the evaluation).",
    "The P7 / P8 / P10 sections are first order in the gains (control-v3 gains, no refit), as P7-P10; the P9 section refits them.  The ep elastic-window",
    "  loss of the P10 section uses A's per-ring window (no refit); the benchmark-pool counts on the P9 'lost' page are the refit counterpart.",
    "Paired errors: bootstrap of the shared clusters with the gains fixed per arm; gain changes against A's sigma_g, with A and D sharing most of the data, so",
    "  a 2 sigma_g move is rarer than for independent solves.",
    "Moller cache: the hook line in moller_calib.pair_job changes its source hash, so the default candidate caches rebuild once (results unchanged).",
    "Member counts (P8 section) count cluster memberships (a module in two clusters counts twice); P7's row counts count rows.",
]


def fig_diffs(pdf, plt, rows):
    fig, axs = plt.subplots(1, 2, figsize=(16, 6.5))
    cols = {"D": "tab:green", "Dp": "tab:olive", "V3": "tab:red"}
    cal = [(rt, "ep") for rt in RTS] + [(TA, "moller")]
    for ax, (q_ep, q_m, title) in zip(axs, (("res_cb", "res", "d sigma/E (oos, percentage points)"), ("mu_ep", "peak", "d peak (oos, %)"))):
        for j, c in enumerate(("D", "Dp", "V3")):
            xs, ys, es = [], [], []
            for i, (rt, calib) in enumerate(cal):
                r = one(rows, table="paired", run_type=rt, calib=calib, comparison=f"{c}-A", protocol="oos", axis="overall",
                        quantity=q_ep if calib == "ep" else q_m, **({} if calib == "ep" else {"bin": "sum"}))
                if r and "diff" in r:
                    xs.append(i + (j - 1) * 0.2)
                    ys.append(fnum(r.get("diff")))
                    es.append(fnum(r.get("se")))
            ax.errorbar(xs, ys, es, fmt="o", color=cols[c], label=f"{ARM_LABEL[c]} - A", capsize=3)
        ax.axhline(0, color="k", lw=0.8)
        ax.set_xticks(range(len(cal)))
        ax.set_xticklabels([RT_LABEL[rt] + ("\nMoller sum" if calib == "moller" else "\nep") for rt, calib in cal], fontsize=9)
        ax.set_title(title)
        ax.legend(fontsize=9)
    fig.suptitle("P9 under drop: paired differences to the default arm (bootstrap errors; V3 from P9)", fontsize=14)
    fig.tight_layout()
    pdf.savefig(fig)
    plt.close(fig)


# ---------------------------------------------------------------------------------------------------------------- stages recluster / analyse / report, main
def recluster_rows(nproc):
    """P7, P8 and P10 (ep) under drop (exotic_p13_recluster.py, own cache)."""
    import exotic_p13_recluster as rc
    return [dict(r) for r in rc.tables(*rc.collect(nproc, str(rc.CACHE)))]


def ep_cost(rows):
    """{run type: (k, n, f, lo, hi)}: ep elastic-window clusters of A lost under D (P10 section rows of exotic_p13_recluster)."""
    import exotic_p13_recluster as rc
    return rc.ep_cost(rows)


def stage_analyse(root, out, nboot, nproc):
    link_a(root)
    M = p9.module_info(root)
    rows = [dict(r) for r in tsv(Path(root) / "hooks" / "hooks.tsv")]
    rows += bench_rows(root)
    rows += boot_rows(root, nboot, nproc)
    rows += lost_rows(root)
    rows += gain_rows(root, M)
    rows += moller_yield_rows(root, M)
    rows += p9_reference_rows()
    rows += recluster_rows(nproc)
    write_rows(out / f"{STEM}.tsv", rows)
    print(f"[tsv] {out / f'{STEM}.tsv'}: {len(rows)} rows", flush=True)
    report(rows, out)


def report(rows, out):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import exotic_p13_recluster as rc
    pdf_path = out / f"{STEM}.pdf"
    with PdfPages(pdf_path) as pdf:
        text_page(pdf, plt, "P13: all impact studies recomputed under channel-level drop of f1-f5",
                  [x.rstrip() for x in __doc__.split("Purpose, inputs, outputs.")[1].splitlines()] + [""] + hook_lines(rows), size=8.4)
        for sec in ("P7", "P8"):
            rc.pages(pdf, plt, rows, sec)
        text_page(pdf, plt, "P9 under drop: sigma/E, peak, tails (arms A, D, D'; V3 from P9)", table_lines(rows), size=7.2)
        fig_diffs(pdf, plt, rows)
        text_page(pdf, plt, "P9 under drop: per ring", ring_lines(rows), size=7.4)
        text_page(pdf, plt, "P9 under drop: clusters of A lost or re-seeded under D (not paired)", lost_lines(rows), size=7.6)
        text_page(pdf, plt, "P9 under drop: per-module gain changes, Ta block H", gain_lines(rows), size=7.4)
        rc.pages(pdf, plt, rows, "P10")
        text_page(pdf, plt, "P10 under drop: Ta Moller pairs", moller_yield_lines(rows, p10_reference()), size=8.6)
        text_page(pdf, plt, "Gate under drop", gate_lines(rows, ep_cost(rows)), size=7.6)
        text_page(pdf, plt, "Caveats and open questions", CAVEATS, size=9.0)
    print(f"[pdf] {pdf_path}", flush=True)


def hook_lines(rows):
    L = ["Dropped channels per run type (in-scope exotic-scan rows, exclusive category; one row per channel and event):", ""]
    L.append(f"  {'run type':<18}{'splits':>8}{'channels':>10}" + "".join(f"{'f' + str(c):>9}" for c in (4, 5, 1, 2, 3)))
    for rt in RTS:
        g = lambda q: int(fnum(one(rows, table="hook", run_type=rt, quantity=q).get("value")))   # noqa: E731
        L.append(f"  {RT_LABEL[rt]:<18}{g('splits'):>8}{g('drop_channels'):>10}" + "".join(f"{g(f'drop_rows_f{c}'):>9}" for c in (4, 5, 1, 2, 3)))
    return L


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--stage", default="all", choices=("hooks", "run", "analyse", "report", "all"))
    ap.add_argument("--jobs", type=int, default=2, help="benchmark reports at a time (each uses the benchmark's NPROC workers)")
    ap.add_argument("--root", default=str(SCAN_ROOT / "exotic_p13"))
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--nboot", type=int, default=p9.NBOOT)
    ap.add_argument("--nproc", type=int, default=12)
    a = ap.parse_args()
    if a.stage in ("hooks", "all"):
        stage_hooks(a.root)
    if a.stage in ("run", "all"):
        stage_run(a.root, a.jobs)
    if a.stage in ("analyse", "all"):
        stage_analyse(a.root, Path(a.out), a.nboot, a.nproc)
    if a.stage == "report":
        report([dict(r) for r in tsv(Path(a.out) / f"{STEM}.tsv")], Path(a.out))


if __name__ == "__main__":
    main()
