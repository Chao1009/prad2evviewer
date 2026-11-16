#!/usr/bin/env python3
"""p7_table.py: P7 table (D12).  Reads output/energy_resolution/variants/P7/<variant>/<RT>/{regions_splits.tsv,gains_regions.tsv} against the
EVAL_REGIONS control variants/P7/ctrl/<RT> (same code, TRAIN_REGIONS = ()) and writes improvement/P7_beyond_pool.tsv.
Per (variant, run type, path, region): paired variant - control over the (run, split) fits (mean +- SE, SE = std(ddof 1) / sqrt(n)) of res_cb in sample and
out of sample (oos = even / odd cross fit), mu, tail_hi_f, tail_lo_f, leak_cut, eff_cut, sep_pi (all out of sample unless marked in), the U2 gap
(oos - in) and the control's gap, the control's values, and the gain statistics of the region (modules in the region, with a gain row, free = prior
fraction < 0.5, mean prior fraction) from gains_regions.tsv (mean over runs).  Verdict per D3 / U2 (res_cb oos: better / worse at >= 2 SE and >= 0.01 pp;
guards: eff_cut down, tail_lo_f up, tail_hi_f up, U2 gap).
Usage: p7_table.py [variant ...]   (default: every variant dir found; prints markdown)."""
import csv, math, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
V7 = ROOT / "output/energy_resolution/variants/P7"
ER = Path(__file__).resolve().parents[1]
GATE = {"0p73GeV_H2": "g20", "2p24GeV_H2": "g16", "3p49GeV_H2": "g16", "3p49GeV_H2_old": "g16", "2p24GeV_Ta_X17": "g8"}
REGIONS = ("pool", "I1", "I2", "H", "M", "O", "T")
Q = ("res_cb", "res_cb_oos", "mu_oos", "tail_hi_f_oos", "tail_lo_f_oos", "leak_cut_oos", "eff_cut_oos", "sep_pi_oos", "tail_hi_f", "tail_lo_f", "leak_cut", "eff_cut")


def rd(p):
    with open(p) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def mean_se(v):
    v = np.array([x for x in v if math.isfinite(x)], float)
    if v.size == 0:
        return float("nan"), float("nan")
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size)) if v.size > 1 else float("nan")


def modreg():
    return {r["module"]: r["region"] for r in rd(V7 / "module_regions.tsv")}


MR = None


def gain_stats(d, region):
    """Per run of the report's gains.tsv (the region solve's, d / ring, when GAIN_REGION_SCOPE ran), averaged over runs: modules of the region,
    with a gain row, free (prior fraction < 0.5), seeded (>= 30 clusters seeded in the module inside its window, MIN_PER_MODULE), the mean prior
    fraction over all the region's modules (no row = 1) and the median n_seed of the seeded modules."""
    global MR
    MR = MR or modreg()
    src = d / "ring" if (d / "ring" / "gains.tsv").exists() else d
    rows = [r for r in rd(src / "gains.tsv") if r["moller"] == "0" and r["nonlinear"] == "0"]
    mods = [m for m, r in MR.items() if r in region.split("+")]
    out = []
    for run in sorted({r["run"] for r in rows}):
        got = {r["module"]: r for r in rows if r["run"] == run}
        pf = np.array([fl(got[m]["prior_frac"]) if m in got else 1.0 for m in mods])
        ns = np.array([fl(got[m]["n_seed"]) if m in got else 0.0 for m in mods])
        out.append((len(mods), sum(m in got for m in mods), int((pf < 0.5).sum()), int((ns >= 30).sum()), float(pf.mean()),
                    float(np.median(ns[ns >= 30])) if (ns >= 30).any() else float("nan")))
    if not out:
        return {}
    m = np.nanmean(np.array(out, float), axis=0)
    return dict(n_modules=int(m[0]), n_rows=m[1], n_free=m[2], n_seeded30=m[3], prior_frac=m[4], nseed_med=m[5], n_runs=len(out))


def verdict(row):
    d, se = row["d_res_cb_oos"], row["d_res_cb_oos_se"]
    core = "better" if d <= -2 * se and d <= -0.01 else "worse" if d >= 2 * se and d >= 0.01 else "null"
    g = []
    if row["d_eff_cut_oos"] < -2 * row["d_eff_cut_oos_se"] and row["d_eff_cut_oos"] < -0.1:
        g.append("eff_cut down")
    if row["d_tail_lo_f_oos"] > 2 * row["d_tail_lo_f_oos_se"] and row["d_tail_lo_f_oos"] > 0.05:
        g.append("tail_lo_f up")
    if row["d_tail_hi_f_oos"] > 2 * row["d_tail_hi_f_oos_se"] and row["d_tail_hi_f_oos"] > 0.02:
        g.append("tail_hi_f up")
    if row["d_gap"] > max(2 * row["d_gap_se"], 0.02):
        g.append("U2 gap")
    return core + ("; guard: " + ", ".join(g) if g else "")


def rows_for(v, rt):
    d, c = V7 / v / rt, V7 / "ctrl" / rt
    if not (d / "regions_splits.tsv").exists() or not (c / "regions_splits.tsv").exists():
        return []
    key = lambda r: (r["run"], r["split"], r["path"], r["kind"], r["region"])   # noqa: E731
    cv = {key(r): r for r in rd(c / "regions_splits.tsv")}
    vv = {key(r): r for r in rd(d / "regions_splits.tsv")}
    ovr = open(d / "regions_splits.tsv").readline().strip().removeprefix("# overrides:").strip()
    out = []
    for path in ("legacy", GATE[rt]):
        for rg in REGIONS:
            ks = [k for k in cv if k in vv and k[2] == path and k[3] == "ep" and k[4] == rg]
            if not ks:
                continue
            r = dict(variant=v, run_type=rt, path=path, region=rg, overrides=ovr)
            for q in Q:
                dm, ds = mean_se([fl(vv[k][q]) - fl(cv[k][q]) for k in ks])
                r[f"d_{q}"], r[f"d_{q}_se"] = dm, ds
                r[f"ctrl_{q}"] = mean_se([fl(cv[k][q]) for k in ks])[0]
                r[f"var_{q}"] = mean_se([fl(vv[k][q]) for k in ks])[0]
            gap = lambda t, k: fl(t[k]["res_cb_oos"]) - fl(t[k]["res_cb"])   # noqa: E731
            r["ctrl_gap"] = mean_se([gap(cv, k) for k in ks])[0]
            r["var_gap"] = mean_se([gap(vv, k) for k in ks])[0]
            r["d_gap"], r["d_gap_se"] = mean_se([gap(vv, k) - gap(cv, k) for k in ks])
            r["n_pairs"] = sum(math.isfinite(fl(vv[k]["res_cb"])) and math.isfinite(fl(cv[k]["res_cb"])) for k in ks)
            r["n_range"] = mean_se([fl(vv[k]["n_range"]) for k in ks])[0] * len(ks)
            reg_ = "I1+I2+I3" if rg == "pool" else rg
            gs = gain_stats(d, reg_)
            r |= {f"gc_{k}": x for k, x in gain_stats(c, reg_).items()}
            r |= {f"g_{k}": x for k, x in gs.items()}
            r["verdict"] = verdict(r)
            out.append(r)
    return out


def main():
    vs = sys.argv[1:] or sorted(p.name for p in V7.iterdir() if p.is_dir() and p.name not in ("ctrl", "logs"))
    allr = []
    for v in vs:
        for rt in GATE:
            allr += rows_for(v, rt)
    if not allr:
        sys.exit("no rows")
    cols = list(allr[0].keys())
    for r in allr:
        for k in r:
            if k not in cols:
                cols.append(k)
    with open(ER / "P7_beyond_pool.tsv", "w") as f:
        f.write("# P7 beyond the pool; delta = variant - EVAL_REGIONS control (variants/P7/ctrl), paired over (run, split); res_cb etc. in % (oos = even/odd cross fit), deltas in pp; g_* from gains_regions.tsv (mean over runs; free = prior fraction < 0.5)\n")
        f.write("\t".join(cols) + "\n")
        for r in allr:
            f.write("\t".join("" if r.get(c) is None else (f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c])) for c in cols) + "\n")
    print("| variant | run type | path | region | ctrl oos | d oos | d in | d gap (ctrl gap) | d mu oos | d tail_hi_f | d tail_lo_f | d leak | d eff | d sep_pi | seeded/free/mod (prior) | verdict |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    pm = lambda r, q, p=3: f"{r['d_' + q]:+.{p}f} ± {r['d_' + q + '_se']:.{p}f}"   # noqa: E731
    for r in allr:
        gf = f"{r.get('g_n_seeded30', float('nan')):.0f}/{r.get('g_n_free', float('nan')):.0f}/{r.get('g_n_modules', 0)} ({r.get('g_prior_frac', float('nan')):.2f})"
        print(f"| {r['variant']} | {r['run_type']} | {r['path']} | {r['region']} | {r['ctrl_res_cb_oos']:.3f} | {pm(r, 'res_cb_oos')} | {pm(r, 'res_cb')} | "
              f"{pm(r, 'res_cb_oos')[:0]}{r['d_gap']:+.3f} ({r['ctrl_gap']:+.3f}) | {pm(r, 'mu_oos')} | {pm(r, 'tail_hi_f_oos', 2)} | {pm(r, 'tail_lo_f_oos', 2)} | "
              f"{pm(r, 'leak_cut_oos', 2)} | {pm(r, 'eff_cut_oos', 2)} | {pm(r, 'sep_pi_oos', 2)} | {gf} | {r['verdict']} |")


if __name__ == "__main__":
    main()
