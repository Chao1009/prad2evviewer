#!/usr/bin/env python3
"""P16 D3 table: variant - control, paired per (run, split), control-gate path, pool from regions_splits.tsv and edge / centre refitted from the
two reports' thr_clusters.npz (THR_DUMP=True) with the same fit_path.

python p16_d3.py <control root> <variant root> [--cross <even-table root> <odd-table root>] [--rts ...] [--tsv out.tsv]

Columns (pp or % points; mean +- SE over pairs): res_cb in sample and out of sample; tail_hi_f, tail_lo_f, tail_hi3, leak_cut, eff_cut out of
sample; peak shift (mu in sample); --cross: each split from the report whose table was built on the other half of the splits (the
table out of sample too); yield change 1 - sum yield_n(variant) / sum yield_n(control) is printed as a loss (%; splits.tsv).
"""
import argparse, csv, math
import numpy as np
import p16_common as pc
bm = pc.bm

QP = (("res_cb", "res_cb"), ("res_cb_oos", "res_cb_oos"), ("tail_hi_f_oos", "tail_hi_f_oos"), ("tail_lo_f_oos", "tail_lo_f_oos"),
      ("tail_hi3_oos", "tail_hi3_oos"), ("leak_cut_oos", "leak_cut_oos"), ("eff_cut_oos", "eff_cut_oos"), ("mu", "mu"))
QD = ("res_cb", "res_cb_oos", "tail_hi_f_oos", "tail_lo_f_oos", "tail_hi3_oos", "leak_cut_oos", "eff_cut_oos", "mu")


def rows(path):
    return list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))


def gate(rt):
    return f"g{bm.CONTROL_GATE[rt]:g}"


def pool(root, rt):
    return {(r["run"], r["split"]): r for r in rows(f"{root}/{rt}/regions_splits.tsv")
            if r["path"] == gate(rt) and r["kind"] == "ep" and r["region"] == "pool"}


def yields(root, rt):
    return {(r["run"], r["split"]): float(r["yield_n"]) for r in rows(f"{root}/{rt}/splits.tsv") if r["variant"] == "all" and r["path"] == gate(rt)}


def parity(keys):
    """(run, split) -> ordinal of the split within its run modulo 2 (the report's even / odd halves)."""
    out = {}
    for run in sorted({k[0] for k in keys}):
        sp = sorted((k for k in keys if k[0] == run), key=lambda k: int(k[1]))
        out |= {k: i % 2 for i, k in enumerate(sp)}
    return out


def merge(a, b, par):
    """Per split: from a (table built on the even half) where the split is odd, else from b (odd-half table)."""
    return {k: (a if par[k] == 1 else b)[k] for k in par if k in a and k in b}


def dump_fits(root, rt, geo):
    """{region: {(run, split): {q: value}}} for edge / centre (and pool, as a check) from the dump."""
    z = pc.load(f"{root}/{rt}/thr_clusters.npz")
    u, v = bm.uv(z, geo)
    edge = (np.abs(u) > bm.EDGE_UV) | (np.abs(v) > bm.EDGE_UV)
    out = {}
    for k in sorted(set(zip(z["run"].tolist(), z["split"].tolist()))):
        sp = (z["run"] == k[0]) & (z["split"] == k[1])
        for rg, rm in (("pool", np.ones(edge.size, bool)), ("pool_edge", edge), ("pool_centre", ~edge)):
            fi, fo = pc.fit(z, sp & rm, geo, rt), pc.fit(z, sp & rm, geo, rt, "e_oos", True) if np.isfinite(z["e_oos"][sp]).any() else None
            if fi is None:
                continue
            d = dict(res_cb=fi["res_cb"], mu=fi["mu_ep"])
            if fo is not None:
                d |= {f"{q}_oos": fo[q] for q in ("res_cb", "tail_hi_f", "tail_lo_f", "tail_hi3", "leak_cut", "eff_cut")}
            out.setdefault(rg, {})[(str(k[0]), str(k[1]))] = d
    return out, float(np.mean(z["er"])), float(np.mean(z["nmem"]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("control")
    ap.add_argument("variant", help="variant root (ignored with --cross)")
    ap.add_argument("--cross", nargs=2, metavar=("EVEN_ROOT", "ODD_ROOT"), help="variant roots built with the even / odd half table")
    ap.add_argument("--rts", nargs="+", default=["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old"])
    ap.add_argument("--tsv")
    a = ap.parse_args()
    out = []
    for rt in a.rts:
        geo = pc.geometry(rt)
        pc_ = pool(a.control, rt)
        dc, *ac = dump_fits(a.control, rt, geo)
        if a.cross:                    # each split scored with the table built on the other half
            par = parity(list(pc_))
            pv = merge(pool(a.cross[0], rt), pool(a.cross[1], rt), par)
            yv = merge(yields(a.cross[0], rt), yields(a.cross[1], rt), par)
            de, *ae = dump_fits(a.cross[0], rt, geo)
            do, *ao = dump_fits(a.cross[1], rt, geo)
            dv, av = {rg: merge(de[rg], do[rg], par) for rg in de}, [0.5 * (x + y) for x, y in zip(ae, ao)]
        else:
            pv, yv = pool(a.variant, rt), yields(a.variant, rt)
            dv, *av = dump_fits(a.variant, rt, geo)
        res = {"pool": {q: pc.mse([float(pv[k][c]) - float(pc_[k][c]) for k in pc_ if k in pv]) for q, c in QP}}
        chk = max(abs(dc["pool"][k]["res_cb"] - float(pc_[k]["res_cb"])) for k in dc["pool"] if k in pc_)
        for rg in ("pool_edge", "pool_centre"):
            res[rg] = {q: pc.mse([dv[rg][k][q] - dc[rg][k][q] for k in dc[rg] if k in dv[rg] and q in dc[rg][k] and q in dv[rg][k]]) for q in QD}
        yc = yields(a.control, rt)
        yl = 100 * (1 - sum(yv[k] for k in yv if k in yc) / sum(yc[k] for k in yv if k in yc))
        ctrl = {q: pc.mse([float(pc_[k][c]) for k in pc_]) for q, c in QP}
        print(f"== {rt}: yield loss {yl:+.2f} %, mean er +{av[0] - ac[0]:.1f} MeV, members +{av[1] - ac[1]:.2f} per cluster; dump-vs-official pool res_cb max |diff| {chk:.2g}")
        print("   control pool: " + "  ".join(f"{q} {ctrl[q][0]:.3f}" for q in QD))
        for rg, d in res.items():
            print(f"   {rg:11s} " + "  ".join(f"{q} {m:+.4f}+-{s:.4f}" for q, (m, s, n) in d.items()))
            out += [dict(rt=rt, region=rg, q=q, mean=m, se=s, n=n) for q, (m, s, n) in d.items()]
        out.append(dict(rt=rt, region="pool", q="yield_loss_pct", mean=yl, se=float("nan"), n=0))
        out += [dict(rt=rt, region="pool", q=f"control_{q}", mean=m, se=s, n=n) for q, (m, s, n) in ctrl.items()]
        out.append(dict(rt=rt, region="pool", q="added_er_MeV", mean=av[0] - ac[0], se=float("nan"), n=0))
    if a.tsv:
        with open(a.tsv, "w") as f:
            f.write(f"# P16 D3 (p16_d3.py): variant {a.cross or a.variant} - control {a.control}, paired per (run, split), control-gate path; pool from regions_splits.tsv, edge / centre refitted from thr_clusters.npz\n")
            f.write("run_type\tregion\tquantity\tmean\tse\tn\n")
            for r in out:
                f.write(f"{r['rt']}\t{r['region']}\t{r['q']}\t{r['mean']:.6g}\t{r['se']:.6g}\t{r['n']}\n")


if __name__ == "__main__":
    main()
