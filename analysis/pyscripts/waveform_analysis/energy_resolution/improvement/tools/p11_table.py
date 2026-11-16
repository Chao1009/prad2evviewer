"""P11 table: transition-ring (T) clusters with and without a lead-glass coincidence, from variants/P11/<RT>/t_clusters.npz (benchmark
report --set EVAL_REGIONS=True --set T_LG_STUDY=True).  Run from analysis/pyscripts:  python3 <this file>   -> improvement/P11_boundary_lg.tsv
Per (run, split) T fits with the benchmark's fit_path (same arguments as region_rows), mean +- SE over splits.  Variants:
  a    all T clusters (current treatment)            b    clusters with a coincidence dropped
  c0   lg clusters E * (1 + k)                       (flag only, one factor)
  c1   lg clusters E * (1 + k n_lg)                  (n_lg = lead-glass modules with a hit within LG_COINC_MM; the scans stored no lead-glass energy)
  d1   lg clusters E * (1 + k n_lg (60 - d_lg) / 30) (d_lg = nearest such module, mm)
  f    seed-layer cut: clusters seeded in the outermost PbWO4 layer (layer 0) dropped, layer 1 kept (lg flag ignored)
k: half of the splits (even / odd rank of the split within its run) set k so that the median of the corrected lg E / m of that half equals
the median of its non-lg clusters in 0.5 m .. 1.3 m; the other half is corrected with it (held out)."""
import math, sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import energy_resolution_benchmark as B

O, IMP = Path("output/energy_resolution"), Path("waveform_analysis/energy_resolution/improvement")
RTS = B.RUN_TYPES
Q = ("res_cb", "mu", "tail_hi3", "tail_lo_f", "sep_pi", "eff_cut", "leak_cut", "sigma_e")
VARS = ("a", "b", "c0", "c1", "d1", "f")


def load(rt, path):
    z = np.load(O / "variants/P11" / rt / "t_clusters.npz")
    d = {k[len(path) + 1:]: z[k] for k in z.files if k.startswith(path + "_")}
    return d | dict(hz=float(z["hz"]), ebeam=float(z["ebeam"]), n_mod=int(z["n_mod"]))


def weight(v, d):
    n, dl = d["nlg"].astype(float), np.nan_to_num(d["dlg"].astype(float), nan=60.0)
    return {"c0": np.ones_like(n), "c1": n, "d1": n * (60.0 - dl) / 30.0}[v]


def fit_k(v, d, sel):
    """k of variant v from the clusters of d in the mask sel (a training half)."""
    m, e = d["m"], d["e"] / d["m"]
    lg, w = d["lg"].astype(bool), weight(v, d)
    win = (e > 0.5) & (e < 1.3)
    tgt = np.median(e[sel & ~lg & win])
    t = sel & lg
    f = lambda k: np.median((e * (1.0 + k * w))[t]) - tgt   # noqa: E731
    lo, hi = 0.0, 4.0
    if f(lo) >= 0:
        return 0.0
    for _ in range(60):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if f(mid) < 0 else (lo, mid)
    return 0.5 * (lo + hi)


def fit_split(d, sel, mod, rt, e):
    """Fit of the clusters in the mask sel (energies e); mod = sel, kept for the seed list."""
    kE, kS, m = float(d["kE"][sel][0]), float(d["kS"][sel][0]), float(d["m"][sel][0])
    geo = dict(hz=d["hz"], ebeam=d["ebeam"], x=np.zeros(d["n_mod"]))
    sids = np.unique(d["seed"][mod].astype(np.int64))
    part = {q: d[q][sel] for q in ("seed", "x", "y", "pu", "clip")} | {"e": e[sel]}
    r = B.fit_path(part, dict.fromkeys(sids.tolist(), m), m, m, geo, kE, kS, False, B.SEP_CUT_MEV[rt])[0]
    return dict(r, mu=r["mu_ep"]) if r["n_range"] >= B.REGION_MIN_FIT else None


def run(rt, path):
    d = load(rt, path)
    z = np.load(next((O / "scan" / rt).glob("*.npz")))
    geo = B.clusterers(int(z["run"]), str(z["calib"]), ())[1]
    pw, rw, cl = geo["pwo4"], geo["row"], geo["col"]
    layer = np.minimum.reduce([rw - rw[pw].min(), rw[pw].max() - rw, cl - cl[pw].min(), cl[pw].max() - cl])
    outer = layer[d["seed"].astype(int)] == 0
    key = d["run"] * 10_000 + d["split"]
    ks = {}
    for r in np.unique(d["run"]):
        sp = sorted(np.unique(d["split"][d["run"] == r]))
        ks.update({int(r) * 10_000 + int(s): i % 2 for i, s in enumerate(sp)})
    half = np.array([ks[int(x)] for x in key])
    lg = d["lg"].astype(bool)
    kfit = {v: [fit_k(v, d, half == h) for h in (0, 1)] for v in ("c0", "c1", "d1")}
    res = {v: [] for v in VARS}
    for kk in np.unique(key):
        s = key == kk
        h = int(half[s][0])
        for v in VARS:
            if v == "a":
                e, mod = d["e"], s
            elif v == "b":
                e, mod = d["e"], s & ~lg
            elif v == "f":
                e, mod = d["e"], s & ~outer
            else:
                e, mod = d["e"] * np.where(lg, 1.0 + kfit[v][1 - h] * weight(v, d), 1.0), s   # k of the other half
            r = fit_split(d, mod, mod, rt, e)
            res[v].append(None if r is None else dict(r, key=int(kk), n_all=int(s.sum()), n_lg=int((s & lg).sum())))
    return d, res, kfit


def msd(v):
    return B.mean_se(v)


def main():
    path_of = {rt: f"g{B.CONTROL_GATE[rt]:g}" for rt in RTS}
    rows, out = [], []
    for rt in RTS:
        for path in ("legacy", path_of[rt]):
            if not (O / "variants/P11" / rt / "t_clusters.npz").exists():
                continue
            d, res, kfit = run(rt, path)
            ref = {x["key"]: x for x in res["a"] if x}
            for v in VARS:
                ok = [x for x in res[v] if x]
                row = dict(run_type=rt, path=path, variant=v, n_splits=len(ok), n_clusters=int(d["e"].size), lg_pct=100.0 * float(d["lg"].mean()),
                           k_even=kfit[v][0] if v in kfit else float("nan"), k_odd=kfit[v][1] if v in kfit else float("nan"))
                row["n_range"] = sum(x["n_range"] for x in ok)
                row["range_a"] = sum(x["n_range"] for x in ref.values())
                row["kept_pct"] = 100.0 * row["n_range"] / row["range_a"] if row["range_a"] else float("nan")
                for q in Q:
                    m, se = msd([x[q] for x in ok])
                    pr = [x[q] - ref[x["key"]][q] for x in ok if x["key"] in ref]
                    dm, dse = msd(pr)
                    row |= {q: m, q + "_se": se, "d_" + q: dm, "d_" + q + "_se": dse}
                eff_a = row["eff_cut"] if v == "a" else eff_a
                row["net_pct"] = row["kept_pct"] * row["eff_cut"] / (eff_a if v != "a" else row["eff_cut"])   # clusters passing the E cut, % of a
                rows.append(row)
    cols = list(rows[0])
    with open(IMP / "P11_boundary_lg.tsv", "w") as f:
        f.write("# transition ring (T), tools/p11_table.py; d_* = variant - a per split (mean, SE); k_* = scale factors fitted on the even / odd half\n")
        f.write("\t".join(cols) + "\n")
        f.writelines("\t".join(f"{r[c]:.6g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n" for r in rows)
    for r in rows:
        if r["path"] != "legacy":
            print(f"{r['run_type']:16s} {r['path']:4s} {r['variant']:3s} lg {r['lg_pct']:5.1f}% kept {r['kept_pct']:5.1f}% res {r['res_cb']:.3f}+-{r['res_cb_se']:.3f} "
                  f"net {r['net_pct']:5.1f}% leak {r['leak_cut']:.1f} d {r['d_res_cb']:+.3f}+-{r['d_res_cb_se']:.3f} tail_hi3 {r['tail_hi3']:.2f} sep {r['sep_pi']:.2f} eff {r['eff_cut']:.1f} k {r['k_even']:.2f}/{r['k_odd']:.2f}")


if __name__ == "__main__":
    main()
