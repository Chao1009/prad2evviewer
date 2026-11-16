#!/usr/bin/env python3
"""
transfer.py - transfer to the production HV through LMS ratios and pool validation (plan P14.4: D5).

Called through the snake_scan.py dispatcher (convention add_parser/run):
  snake_scan.py transfer --rt RT [--gains PATH] [--lscan PATH] [--merged] [--tag TAG] [--no-prior]
Inputs (SS/out): edge_23701.tsv (module, s, s_err), lms_scan.tsv (L, err; 23701 alone, primary), lms_prod_<rt>.tsv (L, err), and the
benchmark's ring/gains.tsv (V3/v3_defaults/<rt>/ring/gains.tsv; pool and region rows with prior_frac < 0.5 are
the solved gains).  --merged replaces L_scan by lms_scan_merged.tsv (23701
merged with 23744 per module within 2 SE; the sensitivity of the 23744 merge, tag "merged").  Outputs: out/transfer_<rt>[_tag].tsv (module, region, s, s_err, L_scan,
L_prod, g_T, sigma_T, g_pool, r, + sigma_pool, group, exp_r, L_scan_err, L_prod_err, L_source) and, unless
--no-prior, out/prior_<rt>.tsv (module, g, sigma_g; O and T modules only).

transfer(edge, lscan, lprod, gains, regions, rt, ...) -> dict(rows, N, spread, ...):
  x_i  = g_pool,i s_i L_prod,i / L_scan,i  over the scanned modules with a solved gain (D5)
  N_rt = median x_i over the pool group;  g_T,i = N_rt (L_scan,i / L_prod,i) / s_i;  r_i = g_pool,i / g_T,i
  spread = 1.48 MAD of r_i;  expected = per-module hypot(sigma_g/g, s_err/s, L_scan err, L_prod err) (relative)
  gate: spread <= GATE (2 %), assessed only with >= MIN_POOL pool modules.
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

SS = Path(__file__).resolve().parent
OUT = SS / "out"
PYSCRIPTS = SS.parents[1]
V3 = PYSCRIPTS / "output" / "energy_resolution" / "control_phase2_v3_20261114" / "v3_defaults"
GATE = 0.02
MIN_POOL = 3
PRIOR_FRAC_MAX = 0.5
MAD_K = 1.4826
for p in (str(SS),):
    if p not in sys.path:
        sys.path.insert(0, p)


def read_tsv(path):
    """(rows as dicts of str, comment lines)."""
    lines = Path(path).read_text().splitlines()
    com = [ln for ln in lines if ln.startswith("#")]
    body = [ln for ln in lines if ln and not ln.startswith("#")]
    head = body[0].split("\t")
    return [dict(zip(head, ln.split("\t"))) for ln in body[1:]], com


def read_edge(path):
    return {r["module"]: (float(r["s"]), float(r["s_err"])) for r in read_tsv(path)[0]}


def read_L(path):
    return {r["module"]: (float(r["L"]), float(r["err"]), r.get("source", "")) for r in read_tsv(path)[0]}


def read_pool_gains(path):
    """module -> (g, sigma_g, n_runs): inverse-variance mean over the runs' rows with prior_frac < 0.5.
    Region modules solved once over all runs have identical rows; then sigma is that row's sigma (not / sqrt(n))."""
    by = {}
    for r in read_tsv(path)[0]:
        if float(r["prior_frac"]) < PRIOR_FRAC_MAX:
            by.setdefault(r["module"], []).append((float(r["g"]), float(r["sigma_g"])))
    out = {}
    for m, v in by.items():
        g = np.array([a for a, _ in v])
        sg = np.array([b for _, b in v])
        w = 1.0 / sg ** 2
        gm = float((w * g).sum() / w.sum())
        sm = float(sg.mean()) if np.ptp(g) < 1e-6 else float(1.0 / np.sqrt(w.sum()))
        out[m] = (gm, sm, len(v))
    return out


def robust_spread(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if len(x) < 2:
        return float("nan")
    return float(MAD_K * np.median(np.abs(x - np.median(x))))


def transfer(edge, lscan, lprod, gpool, regions, groups=("pool",), gate=GATE, min_pool=MIN_POOL):
    """edge {m: (s, s_err)}, lscan/lprod {m: (L, err, src)}, gpool {m: (g, sigma, n)}, regions {m: region}.
    Returns dict(rows, N, n_N, spread, expected, gate, verdict, n_pool); rows in edge-table order."""
    mods = [m for m in edge if m in lscan and m in lprod]
    x = {}
    for m in mods:
        if m in gpool and regions[m] in groups:
            s = edge[m][0]
            x[m] = gpool[m][0] * s * lprod[m][0] / lscan[m][0]
    N = float(np.median(list(x.values()))) if x else float("nan")
    rows = []
    for m in mods:
        s, se = edge[m]
        Ls, Lse = lscan[m][0], lscan[m][1]
        Lp, Lpe = lprod[m][0], lprod[m][1]
        gT = N * (Ls / Lp) / s
        rel_T = math.hypot(math.hypot(se / s, Lse / Ls), Lpe / Lp)
        row = dict(module=m, region=regions[m], s=s, s_err=se, L_scan=Ls, L_prod=Lp, g_T=gT, sigma_T=gT * rel_T,
                   g_pool=float("nan"), r=float("nan"), sigma_pool=float("nan"), group="-", exp_r=float("nan"),
                   L_scan_err=Lse, L_prod_err=Lpe, L_source=lscan[m][2])
        if m in gpool:
            g, sg, _ = gpool[m]
            row.update(g_pool=g, sigma_pool=sg, r=g / gT, group=regions[m] if regions[m] in ("pool", "H") else "-",
                       exp_r=math.hypot(sg / g, rel_T))
        rows.append(row)
    return dict(rows=rows, N=N, n_N=len(x), **summarize(rows, groups, gate, min_pool))


def summarize(rows, groups=("pool",), gate=GATE, min_pool=MIN_POOL):
    sel = [r for r in rows if r["region"] in groups and np.isfinite(r["r"])]
    n = len(sel)
    spread = robust_spread([r["r"] for r in sel])
    exp = float(math.sqrt(np.mean([r["exp_r"] ** 2 for r in sel]))) if sel else float("nan")
    if n < min_pool:
        verdict = "cannot be assessed"
    else:
        verdict = "pass" if spread <= gate else "fail"
    return dict(n_pool=n, spread=spread, expected=exp, sd=float(np.std([r["r"] for r in sel], ddof=1)) if n > 1 else float("nan"),
                verdict=verdict)


def attribution(rows, lscan_free=False):
    """Spread of x_i = g s L_prod/L_scan variants over the rows with a solved gain (pool group): all factors, without the
    L ratio (g s), without s (g L_prod/L_scan), and g alone; each normalised by its median."""
    out = {}
    sel = [r for r in rows if r["region"] == "pool" and np.isfinite(r["g_pool"])]
    if len(sel) < 2:
        return out
    g = np.array([r["g_pool"] for r in sel]); s = np.array([r["s"] for r in sel])
    q = np.array([r["L_prod"] / r["L_scan"] for r in sel])
    for k, v in (("g s q", g * s * q), ("g s", g * s), ("g q", g * q), ("g", g), ("s q", s * q)):
        out[k] = robust_spread(v / np.median(v))
    out["sd(g s q)"] = float(np.std(g * s * q / np.median(g * s * q), ddof=1))
    return out


FLAGGED = {"W594": "7-member 3x3 (two neighbours are beam-hole positions; containment differs)"}
FIELDS = ("module", "region", "s", "s_err", "L_scan", "L_prod", "g_T", "sigma_T", "g_pool", "r", "sigma_pool", "group",
          "exp_r", "L_scan_err", "L_prod_err", "L_source")


def _f(v):
    return v if isinstance(v, str) else ("nan" if not np.isfinite(v) else f"{v:.6g}")


def write_transfer(path, res, header):
    with open(path, "w") as fh:
        for h in header:
            fh.write(f"# {h}\n")
        fh.write("\t".join(FIELDS) + "\n")
        for r in res["rows"]:
            fh.write("\t".join(_f(r[k]) for k in FIELDS) + "\n")


def write_prior(path, res, header):
    """O and T dwell modules: g = g_T, sigma_g = max(pool spread, per-module sigma_T / g_T x g_T); relative terms."""
    sp = res["spread"] if np.isfinite(res["spread"]) else 0.0
    with open(path, "w") as fh:
        for h in header:
            fh.write(f"# {h}\n")
        fh.write("module\tg\tsigma_g\n")
        for r in res["rows"]:
            if r["region"] in ("O", "T"):
                sg = max(sp * r["g_T"], r["sigma_T"])
                fh.write(f"{r['module']}\t{r['g_T']:.6f}\t{sg:.6f}\n")


def region_of(map_path=None):
    import cosmic_mip as CM
    g = CM.region_groups(map_path) if map_path else CM.region_groups()
    return {m: k for k, v in g.items() for m in v}


def build(rt, gains=None, merged=False, lscan_path=None, out_dir=OUT, tag="", prior=True, quiet=False):
    gains = Path(gains) if gains else V3 / rt / "ring" / "gains.tsv"
    edge_p, lp = out_dir / "edge_23701.tsv", out_dir / f"lms_prod_{rt}.tsv"
    ls_p = Path(lscan_path) if lscan_path else out_dir / ("lms_scan_merged.tsv" if merged else "lms_scan.tsv")
    reg = region_of()
    ls = read_L(ls_p)
    res = transfer(read_edge(edge_p), ls, read_L(lp), read_pool_gains(gains), reg)
    hdr = [f"source: edge {edge_p.name}; L_scan {ls_p.name}; "
           f"L_prod {lp.name}; gains {gains}; map hycal_map_til_run023849.json (regions as energy_resolution_benchmark)",
           f"N_rt = {res['N']:.6f} (median of g s L_prod/L_scan over {res['n_N']} pool modules); pool spread (1.48 MAD of r) = "
           f"{res['spread']:.4f} over {res['n_pool']} modules, expected {res['expected']:.4f}; gate (<= {GATE}) {res['verdict']}",
           "g_T = N (L_scan/L_prod)/s; sigma_T = g_T hypot(s_err/s, L_scan err, L_prod err) (N error not included); r = g_pool/g_T; "
           "g_pool = inverse-variance mean over the run rows with prior_frac < 0.5 of ring/gains.tsv (nan: none)"]
    sfx = f"_{tag}" if tag else ""
    write_transfer(out_dir / f"transfer_{rt}{sfx}.tsv", res, hdr)
    if prior:
        write_prior(out_dir / f"prior_{rt}.tsv", res, hdr[:2] + ["g = g_T; sigma_g = max(pool spread x g_T, sigma_T); O and T dwell modules only"])
    return res


def mip_groups(path=OUT / "mip_23744.tsv"):
    return {r["group"]: r for r in read_tsv(path)[0]}


def study_rows(rts=("2p24GeV_H2", "3p49GeV_H2"), out_dir=OUT):
    """Rows of ER/improvement/P14_snake_scan.tsv: per run type, L_scan variant (primary 23701 alone | merged with 23744) and N group
    (pool | pool+H | pool+H+M): n, N, spread, sd, expected, verdict, attribution spreads, leave-one-out range (pool)."""
    reg = region_of()
    edge = read_edge(out_dir / "edge_23701.tsv")
    scans = {"primary": read_L(out_dir / "lms_scan.tsv"), "merged": read_L(out_dir / "lms_scan_merged.tsv")}
    rows = []
    for rt in rts:
        gp = read_pool_gains(V3 / rt / "ring" / "gains.tsv")
        lp = read_L(out_dir / f"lms_prod_{rt}.tsv")
        for sv, ls in scans.items():
            for gname, groups in (("pool", ("pool",)), ("pool+H", ("pool", "H")), ("pool+H-W594", ("pool", "H")),
                                  ("pool+H+M", ("pool", "H", "M"))):
                ed = {m: v for m, v in edge.items() if m not in FLAGGED} if gname.endswith("-W594") else edge
                res = transfer(ed, ls, lp, gp, reg, groups=groups)
                sel = [r for r in res["rows"] if r["region"] in groups and np.isfinite(r["r"])]
                at = attribution([dict(r, region="pool") for r in sel]) if len(sel) > 1 else {}
                loo = []
                if gname == "pool":
                    for drop in sel:
                        gg = {m: v for m, v in gp.items() if m != drop["module"]}
                        loo.append(transfer(edge, ls, lp, gg, reg, groups=groups)["spread"])
                rows.append(dict(rt=rt, L_scan=sv, N_group=gname, n=res["n_pool"], N=res["N"], spread=res["spread"],
                                 sd=res["sd"], expected=res["expected"], verdict=res["verdict"],
                                 sp_gsq=at.get("g s q", float("nan")), sp_gs=at.get("g s", float("nan")),
                                 sp_gq=at.get("g q", float("nan")), sp_g=at.get("g", float("nan")),
                                 sp_sq=at.get("s q", float("nan")),
                                 loo_min=min(loo) if loo else float("nan"), loo_max=max(loo) if loo else float("nan"),
                                 flags=" ".join(f"{m}:{FLAGGED[m]}" for m in FLAGGED if any(r["module"] == m for r in sel)) or "-",
                                 r_by_module=" ".join(f"{r['module']}={r['r']:.4f}" for r in sel)))
    return rows


def write_study(path, rows):
    keys = list(rows[0])
    with open(path, "w") as fh:
        fh.write("# P14.4 transfer study: pool/region spread of r = g_pool/g_T by run type, L_scan variant and N group; "
                 "sp_* = 1.48 MAD of x/median(x) for x = g s q, g s, g q, g, s q with q = L_prod/L_scan (group rows only); "
                 "loo = leave-one-out range of the pool spread; flags = modules in the group with a caveat (W594's 3x3 has 7 members); "
                 "L_scan primary = 23701 alone (verdict source), merged = 23701 + 23744 per-module merge (sensitivity); gains V3/v3_defaults/<rt>/ring/gains.tsv\n")
        fh.write("\t".join(keys) + "\n")
        for r in rows:
            fh.write("\t".join(_f(r[k]) for k in keys) + "\n")


def add_parser(sub):
    p = sub.add_parser("transfer", help="LMS-ratio transfer, pool validation, prior table (D5)")
    p.add_argument("--rt", required=True)
    p.add_argument("--gains", default=None, help="ring/gains.tsv (default V3/v3_defaults/<rt>/ring/gains.tsv)")
    p.add_argument("--lscan", default=None, help="L_scan table (default out/lms_scan.tsv, 23701 alone)")
    p.add_argument("--merged", action="store_true", help="L_scan from out/lms_scan_merged.tsv (23701 + 23744 merge; sensitivity)")
    p.add_argument("--tag", default="", help="suffix of the transfer table; use with sensitivities")
    p.add_argument("--no-prior", action="store_true")
    p.add_argument("--study", default=None, help="write the study table (both run types, all variants) to this TSV and exit")
    return p


def run(a):
    if a.study:
        rows = study_rows()
        write_study(a.study, rows)
        print(f"study: {len(rows)} rows -> {a.study}")
        return 0
    tag = a.tag or ("merged" if a.merged else "")
    res = build(a.rt, a.gains, a.merged, a.lscan, OUT, tag, prior=not (a.no_prior or tag))
    for r in res["rows"]:
        print(f"{r['module']:5s} {r['region']:5s} s={r['s']:.4f} Ls/Lp={r['L_scan']/r['L_prod']:.4f} g_T={r['g_T']:.4f} "
              f"g_pool={r['g_pool']:.4f} r={r['r']:.4f}")
    print(f"{a.rt}{' ' + tag if tag else ''}: N={res['N']:.5f} (n={res['n_N']}), pool spread={res['spread']:.4f} over "
          f"{res['n_pool']} modules, expected {res['expected']:.4f}, gate (<= {GATE}): {res['verdict']}")
    return 0
