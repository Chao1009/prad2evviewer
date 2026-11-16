"""Paired D3 table of two reports at their control gate paths.  Run from analysis/pyscripts:
  python3 waveform_analysis/energy_resolution/improvement/tools/p2v2_table.py <control dir> <control path> <variant dir> <variant path>   (one run type's dirs; prints one line)
  python3 waveform_analysis/energy_resolution/improvement/tools/p2v2_table.py control   (v2 control table -> improvement/control_phase2_v2.tsv)
  python3 waveform_analysis/energy_resolution/improvement/tools/p2v2_table.py v2      (v2 vs v1 for the five run types -> improvement/control_phase2_v2_delta.tsv)
Pairs per (run, split) the pool ep rows; d = variant - control, mean +- SE, real if |d| >= 2 SE (D3)."""
import csv, math, sys
from pathlib import Path
import numpy as np

O = Path("output/energy_resolution")
IMP = Path("waveform_analysis/energy_resolution/improvement")
V1, V2 = O / "control_phase2_20261113", O / "control_phase2_v2_20261113"
GATE1 = {"0p73GeV_H2": "g20", "2p24GeV_H2": "g16", "3p49GeV_H2": "g16", "3p49GeV_H2_old": "g16", "2p24GeV_Ta_X17": "g8"}
GATE2 = GATE1 | {"2p24GeV_Ta_X17": "g4"}
Q = ("res_cb", "res_cb_oos", "tail_hi_f", "tail_lo_f", "tail_hi3", "leak_cut", "eff_cut", "mu")
rows = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def mse(v):
    v = np.array([x for x in v if math.isfinite(x)])
    return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size)), int(v.size)) if v.size > 1 else (float("nan"), float("nan"), int(v.size))


def pool(d, path):
    R = {(r["run"], r["split"]): r for r in rows(d / "regions_splits.tsv") if r["path"] == path and r["kind"] == "ep" and r["region"] == "pool"}
    S = {(r["run"], r["split"]): r for r in rows(d / "splits.tsv") if r["path"] == path and r["variant"] == "all"}
    return R, S


def pair(cd, cp, vd, vp):
    """{quantity: (control mean, variant mean, d, SE, n)}, plus pu (pile-up flag) and yield (variant yield_n / control, % change)."""
    Rc, Sc = pool(cd, cp)
    Rv, Sv = pool(vd, vp)
    ks = sorted(set(Rc) & set(Rv))
    out = {}
    for q in Q:
        c = np.array([float(Rc[k][q]) for k in ks])
        v = np.array([float(Rv[k][q]) for k in ks])
        ok = np.isfinite(c) & np.isfinite(v)
        out[q] = (mse(c[ok])[0], mse(v[ok])[0], *mse(v[ok] - c[ok]))
    pu = [(float(Sc[k]["pu"]), float(Sv[k]["pu"])) for k in ks if k in Sc and k in Sv]
    out["pu"] = (mse([a for a, b in pu])[0], mse([b for a, b in pu])[0], *mse([b - a for a, b in pu]))
    yl = [(float(Sc[k]["yield_n"]), float(Sv[k]["yield_n"])) for k in ks if k in Sc and k in Sv]
    out["yield"] = (100.0 * sum(b for a, b in yl) / sum(a for a, b in yl) - 100.0,)
    return out


def verdict(d, se):
    return "real" if math.isfinite(d) and math.isfinite(se) and d != 0 and abs(d) >= 2 * se else "null"


def line(m):
    return "  ".join(f"{q} {m[q][2]:+.4f}+-{m[q][3]:.4f}" for q in (*Q, "pu")) + f"  yield {m['yield'][0]:+.2f}%  n={m['res_cb'][4]}"


def control_table():
    """control_phase2_v2.tsv: one row per run type at the v2 control gate (columns as control_phase2.tsv, tools/p2_control.py), plus tail_hi_f."""
    out = []
    for rt in GATE2:
        S = {r["path"]: r for r in rows(V2 / rt / "summary.tsv")}
        R = [r for r in rows(V2 / rt / "regions_splits.tsv") if r["path"] == GATE2[rt] and r["kind"] == "ep" and r["region"] == "pool"]
        s = S[GATE2[rt]]
        row = dict(run_type=rt, gate=GATE2[rt], n_splits=s["n_splits"], n_used=s["n_used"])
        for q in ("res_cb", "res_cb_oos"):
            row[q], row[q + "_se"] = float(s[q]), float(s[q + "_se"])
        for q in ("tail_hi3", "tail_hi_f", "tail_lo_f", "leak_cut", "eff_cut"):
            row[q], row[q + "_se"], _ = mse([float(r[q]) for r in R])
        row["pu_pct"], row["pu_se"] = float(s["pu"]), float(s["pu_se"])
        row["yield_pct"] = float(s["yield_pct"])
        out.append(row)
        print(rt, " ".join(f"{p}={float(r['res_cb_oos']):.4f}" for p, r in S.items()))
    with open(IMP / "control_phase2_v2.tsv", "w") as f:
        f.write("\t".join(out[0]) + "\n")
        f.writelines("\t".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in r.values()) + "\n" for r in out)


if __name__ == "__main__":
    if sys.argv[1:2] == ["control"]:
        control_table()
    elif sys.argv[1:2] == ["v2"]:
        out = []
        for rt in GATE1:
            m = pair(V1 / rt, GATE1[rt], V2 / rt, GATE2[rt])
            row = dict(run_type=rt, gate_v1=GATE1[rt], gate_v2=GATE2[rt], n_pairs=m["res_cb"][4], yield_pct=m["yield"][0])
            for q in (*Q, "pu"):
                row |= {f"{q}_v1": m[q][0], f"{q}_v2": m[q][1], f"d_{q}": m[q][2], f"d_{q}_se": m[q][3], f"verdict_{q}": verdict(m[q][2], m[q][3])}
            out.append(row)
            print(rt, line(m))
        with open(IMP / "control_phase2_v2_delta.tsv", "w") as f:
            f.write("\t".join(out[0]) + "\n")
            f.writelines("\t".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in r.values()) + "\n" for r in out)
    else:
        cd, cp, vd, vp = sys.argv[1:5]
        print(line(pair(Path(cd), cp, Path(vd), vp)))
