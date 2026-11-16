#!/usr/bin/env python3
"""p2_tables.py: P2 tables from output/energy_resolution/variants/P2/<run type>/{regions,regions_splits}.tsv.
Writes ER/improvement/P2_regions.{md,tsv}; pool continuity against baseline_20261112/pre.  Usage: p2_tables.py"""
import csv, math, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[4]                 # analysis/pyscripts
OUT = ROOT / "output/energy_resolution"
P2, BASE = OUT / "variants/P2", OUT / "baseline_20261112/pre"
DST = Path(__file__).resolve().parents[1]
RTS = {"0p73GeV_H2": "g12", "2p24GeV_H2": "g12", "3p49GeV_H2": "g12", "3p49GeV_H2_old": "g12", "2p24GeV_Ta_X17": "g8"}
REGS = ["pool", "I1", "I2", "I3", "pool_edge", "pool_centre", "H", "M", "O", "T"]


def rd(p):
    with open(p) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def f(x, d=3):
    return "nan" if x is None or not math.isfinite(x) else f"{x:.{d}f}"


def main():
    md, tsv = [], []
    md.append("# P2: region-resolved resolution and separation\n")
    md.append("Run: `report --set EVAL_REGIONS=True --set MOLLER_PAIRS=True --set MOLLER_BEAM_OFFSET=<mean of the run type's runs in "
              "`output/nonlin_20261104/step3/beam_offsets.tsv`, (0,0) where absent>`, gates 8 12 (Ta 8). "
              "Values: mean ± SE over per-split fits, out-of-sample columns marked oos. sep_pi = 145 MeV / sigma_E.\n")
    for rt, gate in RTS.items():
        d = P2 / rt
        if not (d / "regions.tsv").exists():
            continue
        R = rd(d / "regions.tsv")
        get = lambda p, k, r: next((x for x in R if x["path"] == p and x["kind"] == k and x["region"] == r), None)
        md.append(f"\n## {rt}\n")
        for p in ("legacy", gate):
            md.append(f"\n### ep, path {p}\n")
            md.append("| region | status | n_range | modules | prior frac | res_cb (%) | res_cb oos | mu (%) | sigma_E (MeV) | sep_pi | tail_hi3 | leak_cut | eff_cut |")
            md.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
            for r in REGS + ["T_lg", "T_nolg"]:
                x = get(p, "ep", r)
                if x is None:
                    continue
                tsv.append([rt, p, "ep", r, x["status"], x["n_range"], x["res_cb"], x["res_cb_se"], x["res_cb_oos"], x["mu"], x["sep_pi"], x["tail_hi3"], x["leak_cut"], x["eff_cut"], x["prior_frac_mean"]])
                md.append(f"| {r} | {x['status']} | {x['n_range']} | {x['n_modules']} | {f(fl(x['prior_frac_mean']))} | "
                          f"{f(fl(x['res_cb']))} ± {f(fl(x['res_cb_se']))} | {f(fl(x['res_cb_oos']))} | {f(fl(x['mu']))} ± {f(fl(x['mu_se']))} | "
                          f"{f(fl(x['sigma_e']),1)} | {f(fl(x['sep_pi']),2)} | {f(fl(x['tail_hi3']))} | {f(fl(x['leak_cut']),1)} | {f(fl(x['eff_cut']),1)} |")
        md.append(f"\n### moller2 ((E1+E2)/E_beam) and moller (E/E_ee), path legacy\n")
        md.append("| kind | region | status | n | n_range | sigma/E (%) | mu (%) | mu(moller) - mu(ep) (pp) |")
        md.append("|---|---|---|---|---|---|---|---|")
        for k in ("moller2", "moller"):
            for r in (["all"] if k == "moller2" else REGS):
                x = get("legacy", k, r)
                if x is None or x["n"] in ("0",):
                    continue
                ep = get("legacy", "ep", r)
                dm = fl(x["mu"]) - fl(ep["mu"]) if ep and k == "moller" and abs(fl(x["mu"])) < 20 else float("nan")
                tsv.append([rt, "legacy", k, r, x["status"], x["n"], x["res_cb"], "", x["res_g2"], x["mu"], "", "", "", "", dm])
                md.append(f"| {k} | {r} | {x['status']} | {x['n']} | {x['n_range']} | g2 {f(fl(x['res_g2']))} | {f(fl(x['mu']))} ± {f(fl(x['mu_se']))}{' (unreliable fit)' if abs(fl(x['mu'])) >= 20 else ''} | {f(dm)} |")
        # T lead-glass split
        tl, tn, tt = get("legacy", "ep", "T_lg"), get("legacy", "ep", "T_nolg"), get("legacy", "ep", "T")
        if tl and tn:
            nl, nn = int(tl["n"]), int(tn["n"])
            md.append(f"\nT clusters (legacy): with lead-glass coincidence {nl} ({100*nl/max(nl+nn,1):.1f} %), median E/E' {f(fl(tl['med']),4)} "
                      f"(deficit {f(100*(1-fl(tl['med'])),2)} %); without {nn}, median {f(fl(tn['med']),4)} (deficit {f(100*(1-fl(tn['med'])),2)} %).")
        # drift: mu(pool) against split order per run
        S = rd(d / "regions_splits.tsv")
        md.append(f"\n### drift: mu of the pool, legacy, against split order per run\n")
        md.append("| run | splits | mu mean (%) | slope (pp / split) ± SE | mu first-last (pp) |")
        md.append("|---|---|---|---|---|")
        for run in sorted({x["run"] for x in S}):
            rows = sorted([x for x in S if x["run"] == run and x["path"] == "legacy" and x["kind"] == "ep" and x["region"] == "pool" and math.isfinite(fl(x["mu"]))], key=lambda x: int(x["split"]))
            if len(rows) < 3:
                md.append(f"| {run} | {len(rows)} | - | - | - |")
                continue
            xs, ys = np.array([int(x["split"]) for x in rows], float), np.array([fl(x["mu"]) for x in rows])
            A = np.polyfit(xs, ys, 1)
            res = ys - np.polyval(A, xs)
            se = math.sqrt(res.var(ddof=2) / ((xs - xs.mean()) ** 2).sum()) if len(xs) > 2 else float("nan")
            md.append(f"| {run} | {len(rows)} | {f(ys.mean())} | {f(A[0],4)} ± {f(se,4)} | {f(ys[-1]-ys[0])} |")
        # continuity with the baseline pool row
        b = rd(BASE / rt / "summary.tsv")
        md.append("")
        pool = get("legacy", "ep", "pool")
        bl = [x for x in b if x.get("path") == "legacy"] or b[:1]
        if bl and "res_cb" in bl[0] and pool:
            md.append(f"Continuity: pool legacy res_cb {f(fl(pool['res_cb']),4)} vs baseline pre {f(fl(bl[0]['res_cb']),4)}.")
    (DST / "P2_regions.md").write_text("\n".join(md) + "\n")
    with open(DST / "P2_regions.tsv", "w") as fh:
        fh.write("run_type\tpath\tkind\tregion\tstatus\tn_range\tres_cb\tres_cb_se\tres_cb_oos_or_g2\tmu\tsep_pi\ttail_hi3\tleak_cut\teff_cut\tprior_frac_or_dmu\n")
        for r in tsv:
            fh.write("\t".join(map(str, r)) + "\n")


if __name__ == "__main__":
    main()
