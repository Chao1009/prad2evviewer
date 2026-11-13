"""Paired control-vs-variant comparison of two benchmark reports (D3 of the energy-resolution plan).

Usage: python3 waveform_analysis/compare_reports.py <control dir> <variant dir> [--paths legacy g12 ...]
Each dir is one run type's report dir (regions_splits.tsv, splits.tsv); writes <variant dir>/delta.tsv:
per (path, kind, region) the mean +- SE over paired splits of variant - control, the U2 in/oos gap change,
the yield loss (%), and real/null verdicts (|d| >= 2 SE, D3; in and out of sample).
"""
import argparse, csv, math
import numpy as np

QS = ("res_cb", "res_cb_oos", "res_g2", "tail_hi3", "tail_hi3_oos", "tail_lo3", "leak_cut", "leak_cut_oos", "eff_cut", "mu")
VERDICTS = ("res_cb", "res_cb_oos", "tail_hi3", "tail_hi3_oos", "leak_cut", "leak_cut_oos")


def tsv_rows(path):
    """Rows of a tab-separated file (comment lines starting with # skipped) as dicts of strings."""
    return list(csv.DictReader((l for l in open(path) if not l.startswith("#")), delimiter="\t"))


def fmt(v):
    return v if isinstance(v, str) else str(v) if isinstance(v, (int, np.integer)) else f"{v:.6g}"


def mean_se(v):
    """(mean, standard error) of the finite values."""
    v = np.array([x for x in v if math.isfinite(x)], dtype=float)
    if v.size == 0:
        return float("nan"), float("nan")
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else float("nan")


def verdict(d, se):
    return "real" if math.isfinite(d) and math.isfinite(se) and abs(d) >= 2 * se and (d or se) else "null"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("control"); ap.add_argument("variant"); ap.add_argument("--paths", nargs="+")
    a = ap.parse_args()
    key = lambda r: (r["run"], r["split"], r["path"], r["kind"], r["region"])
    rd = lambda d: {key(r): r for r in tsv_rows(f"{d}/regions_splits.tsv")}
    ctl, var = rd(a.control), rd(a.variant)
    ovr = open(f"{a.variant}/regions_splits.tsv").readline().strip().removeprefix("# overrides:").strip()
    yl = lambda d: {(r["run"], r["split"], r["path"]): float(r["yield_n"]) for r in tsv_rows(f"{d}/splits.tsv") if r["variant"] == "all"}
    yc, yv = yl(a.control), yl(a.variant)
    groups = {}
    for k in ctl:
        if k in var and (not a.paths or k[2] in a.paths):
            groups.setdefault(k[2:], []).append(k)
    head = ["path", "kind", "region", "n"] + [c for q in QS for c in (f"d_{q}", f"d_{q}_se")] + ["gap", "gap_se", "yield_loss"] \
        + [f"verdict_{q}" for q in VERDICTS]
    rows = []
    for (path, kind, region), ks in groups.items():
        diff = lambda f: [f(var[k]) - f(ctl[k]) for k in ks]
        g = lambda q: (lambda r: float(r[q]))
        gap = lambda r: float(r["res_cb_oos"]) - float(r["res_cb"])
        st = {q: mean_se(diff(g(q))) for q in QS}
        st["gap"] = mean_se(diff(gap))
        nq = "mu" if kind.startswith("moller") else "res_cb"
        n = sum(math.isfinite(float(var[k][nq])) and math.isfinite(float(ctl[k][nq])) for k in ks)
        ps = [p for p in yc if p[2] == path and p in yv]
        yloss = 100 * (1 - sum(yv[p] for p in ps) / sum(yc[p] for p in ps)) if ps and sum(yc[p] for p in ps) else float("nan")
        v = [verdict(*st[q]) for q in VERDICTS]
        rows.append([path, kind, region, n] + [x for q in QS for x in st[q]] + list(st["gap"]) + [yloss] + v)
    with open(f"{a.variant}/delta.tsv", "w") as f:
        f.write(f"# control: {a.control}; variant: {a.variant}; overrides: {ovr}\n" + "\t".join(head) + "\n")
        f.writelines("\t".join(fmt(x) for x in r) + "\n" for r in rows)
    ix = {h: i for i, h in enumerate(head)}
    pm = lambda r, q: f"{fmt(r[ix['d_' + q]])}±{fmt(r[ix['d_' + q + '_se']])}"
    print("path kind region n d_res_cb d_tail_hi3 d_leak_cut yield_loss verdicts(" + ",".join(VERDICTS) + ")")
    for r in rows:
        print(*r[:4], pm(r, "res_cb"), pm(r, "tail_hi3"), pm(r, "leak_cut"), fmt(r[ix["yield_loss"]]), "/".join(r[-len(VERDICTS):]))


if __name__ == "__main__":
    main()
