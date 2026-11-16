"""P9 tables: paired variant - control Delta +- SE per (run type, path, region) from two reports' regions_splits.tsv.

python p9_table.py <control root> <variants root> <variant> [<variant> ...] [--rts 0p73GeV_H2 ...] [--regions pool pool_edge pool_centre]
        [--tsv out.tsv]
Both roots hold <run type>/regions_splits.tsv; the path is the control gate of the run type.  Values are out of sample (_oos) unless
the column name says in-sample.  SE = std(ddof 1) / sqrt(n) over paired (run, split) fits (D3).
"""
import argparse, csv, math
import numpy as np

GATE = {"0p73GeV_H2": "g20", "2p24GeV_H2": "g16", "3p49GeV_H2": "g16", "3p49GeV_H2_old": "g16", "2p24GeV_Ta_X17": "g8"}
QS = ("res_cb", "eff_cut", "tail_hi3", "tail_hi_f", "tail_lo_f", "leak_cut", "res_cb_oos", "eff_cut_oos", "tail_hi3_oos", "tail_hi_f_oos", "tail_lo_f_oos", "leak_cut_oos")


def rd(p):
    return {(r["run"], r["split"], r["path"], r["kind"], r["region"]): r
            for r in csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")}


def yloss(c, v, rt):
    """1 - sum yield_n(variant) / sum yield_n(control), %, control gate path (splits.tsv, variant all)."""
    t = lambda d: sum(float(r["yield_n"]) for r in csv.DictReader((l for l in open(f"{d}/{rt}/splits.tsv") if not l.startswith("#")), delimiter="\t")
                      if r["variant"] == "all" and r["path"] == GATE[rt])
    return 100 * (1 - t(v) / t(c))


def mse(v):
    v = np.array([x for x in v if math.isfinite(x)], float)
    return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size))) if v.size > 1 else (float("nan"),) * 2


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("control"); ap.add_argument("root"); ap.add_argument("variants", nargs="+")
    ap.add_argument("--rts", nargs="+", default=["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old"])
    ap.add_argument("--regions", nargs="+", default=["pool", "pool_edge", "pool_centre"])
    ap.add_argument("--tsv"); ap.add_argument("--md")
    a = ap.parse_args()
    out = []
    for rt in a.rts:
        ctl = rd(f"{a.control}/{rt}/regions_splits.tsv")
        for var in a.variants:
            try:
                v = rd(f"{a.root}/{var}/{rt}/regions_splits.tsv")
            except FileNotFoundError:
                continue
            for reg in a.regions:
                ks = [k for k in ctl if k[2:] == (GATE[rt], "ep", reg) and k in v]
                row = dict(run_type=rt, variant=var, path=GATE[rt], region=reg, n=len(ks))
                for q in QS:
                    c = mse([float(ctl[k][q]) for k in ks])
                    d = mse([float(v[k][q]) - float(ctl[k][q]) for k in ks])
                    row |= {f"{q}_ctl": c[0], f"d_{q}": d[0], f"d_{q}_se": d[1]}
                g = mse([float(v[k]["res_cb_oos"]) - float(v[k]["res_cb"]) - float(ctl[k]["res_cb_oos"]) + float(ctl[k]["res_cb"]) for k in ks])
                row |= {"d_gap": g[0], "d_gap_se": g[1], "yield_loss": yloss(a.control, f"{a.root}/{var}", rt)}
                out.append(row)
    for r in out:
        print(f"{r['run_type']:15s} {r['variant']:18s} {r['region']:11s} " + " ".join(
            f"{q.replace('_oos','~')}:{r['d_'+q]:+.3f}({r['d_'+q+'_se']:.3f})" for q in QS if q.endswith("_oos") or q == "res_cb"))
    if a.md:
        f = lambda r, q: f"{r['d_' + q]:+.3f} ± {r['d_' + q + '_se']:.3f}"
        with open(a.md, "w") as fh:
            fh.write("| run type | variant | region | control res_cb oos | Δ res_cb oos | Δ eff_cut oos | Δ tail_hi3 oos | Δ tail_hi_f oos | Δ tail_lo_f oos | Δ leak_cut oos | Δ res_cb in | Δ gap | yield loss % |\n" + "|---|" * 0
                     + "|" + "---|" * 13 + "\n")
            for r in out:
                fh.write(f"| {r['run_type']} | {r['variant']} | {r['region']} | {r['res_cb_oos_ctl']:.3f} | " + " | ".join(
                    f(r, q) for q in ("res_cb_oos", "eff_cut_oos", "tail_hi3_oos", "tail_hi_f_oos", "tail_lo_f_oos", "leak_cut_oos", "res_cb", "gap"))
                    + f" | {r['yield_loss']:.2f} |\n")
    if a.tsv:
        with open(a.tsv, "w") as f:
            f.write("\t".join(out[0]) + "\n")
            f.writelines("\t".join(f"{x:.6g}" if isinstance(x, float) else str(x) for x in r.values()) + "\n" for r in out)


if __name__ == "__main__":
    main()
