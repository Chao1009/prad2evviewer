"""P10 D3 table: variant reports vs control v2, per run and pooled (pool ep, control gate path).  Run from analysis/pyscripts:
  python3 waveform_analysis/energy_resolution/improvement/tools/p10_table.py <label> <RT> <variant dir> [<control dir>] [--path g16] [--runs 25228 ...]
Prints one line per run and a pooled line; appends rows to improvement/P10_phase_corr.tsv unless --no-tsv.  d = variant - control (paired per run and split),
mean +- SE, real if |d| >= 2 SE (D3); guards tail_hi_f, tail_lo_f, eff_cut, yield."""
import argparse, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).parent))
import numpy as np
import p2v2_table as t

ap = argparse.ArgumentParser()
ap.add_argument("label"); ap.add_argument("rt"); ap.add_argument("var")
ap.add_argument("ctl", nargs="?"); ap.add_argument("--path", default="g16"); ap.add_argument("--runs", nargs="*"); ap.add_argument("--no-tsv", action="store_true")
a = ap.parse_args()
ctl = Path(a.ctl) if a.ctl else t.V2 / a.rt
var = Path(a.var)
Rc, Sc = t.pool(ctl, a.path)
Rv, Sv = t.pool(var, a.path)
ks = sorted(set(Rc) & set(Rv))
runs = sorted({k[0] for k in ks})
if a.runs:
    runs = [r for r in runs if r in a.runs]
ovr = open(var / "regions_splits.tsv").readline().strip().removeprefix("# overrides:").strip()
bank = sorted({str(np.load(f, mmap_mode="r")["bank_sha1"])[:8] for f in Path("output/energy_resolution/scan").glob(f"{a.rt}/*.npz")})
Q = t.Q
cols = ["label", "run_type", "run", "n", "control_res_cb"] + [c for q in Q for c in (f"d_{q}", f"d_{q}_se")] + ["d_pu", "d_pu_se", "yield_change_pct", "verdict", "overrides"]
out = []
for sel, name in [([k for k in ks if k[0] == r], r) for r in runs] + [([k for k in ks if k[0] in runs], "pooled")]:
    rec = [a.label, a.rt, name, len(sel)]
    c0 = [float(Rc[k]["res_cb"]) for k in sel]
    rec.append(np.nanmean(c0))
    ds = {}
    for q in Q:
        d = np.array([float(Rv[k][q]) - float(Rc[k][q]) for k in sel])
        m, se, _ = t.mse(d[np.isfinite(d)])
        ds[q] = (m, se)
        rec += [m, se]
    pu = np.array([float(Sv[k]["pu"]) - float(Sc[k]["pu"]) for k in sel if k in Sc and k in Sv])
    m, se, _ = t.mse(pu)
    rec += [m, se]
    yc = sum(float(Sc[k]["yield_n"]) for k in sel if k in Sc and k in Sv)
    yv = sum(float(Sv[k]["yield_n"]) for k in sel if k in Sc and k in Sv)
    rec.append(100.0 * (yv / yc - 1.0))
    v = [f"{q} {t.verdict(*ds[q])} {ds[q][0]:+.3f}" for q in ("res_cb", "res_cb_oos", "tail_hi_f", "leak_cut") if t.verdict(*ds[q]) == "real"]
    rec += ["; ".join(v) or "null", ovr]
    out.append(rec)
    print(f"{a.label:10s} {a.rt} {name:>7} n={len(sel):2d}  " + "  ".join(f"{q} {ds[q][0]:+.4f}+-{ds[q][1]:.4f}" for q in ("res_cb", "res_cb_oos", "tail_hi_f", "tail_lo_f", "leak_cut", "eff_cut", "mu"))
          + f"  yield {rec[-3]:+.2f}%  pu {rec[-5]:+.2f}  [{rec[-2]}]")
if not a.no_tsv:
    p = t.IMP / "P10_phase_corr.tsv"
    new = not p.exists()
    with open(p, "a") as f:
        if new:
            f.write("# P10: residual onset phase and PHASE_CORR D3 table; pool ep, control gate path g16 (see P10_phase_corr.md); scan bank_sha1 of the run type: see md\n")
            f.write("\t".join(cols) + "\n")
        for r in out:
            f.write("\t".join(f"{x:.6g}" if isinstance(x, (float, np.floating)) else str(x) for x in r) + "\n")
