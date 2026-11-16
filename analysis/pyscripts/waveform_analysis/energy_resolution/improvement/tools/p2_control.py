"""Phase-2 control table. Run from analysis/pyscripts: python3 waveform_analysis/energy_resolution/improvement/tools/p2_control.py
Reads output/energy_resolution/control_phase2_20261113/<RT>/{summary.tsv,regions_splits.tsv} and the P5 table; writes
improvement/control_phase2.tsv (one row per run type at the control gate; res_cb in/oos from summary.tsv, shares from the pool rows of
regions_splits.tsv, mean +- SE over splits) and prints the P5 comparison and the gate scan (res_cb_oos per path)."""
import csv, math
from pathlib import Path
import numpy as np

O = Path("output/energy_resolution/control_phase2_20261113")
IMP = Path("waveform_analysis/energy_resolution/improvement")
RTS = ("0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old", "2p24GeV_Ta_X17")
GATE = {"0p73GeV_H2": "g20", "2p24GeV_H2": "g16", "3p49GeV_H2": "g16", "3p49GeV_H2_old": "g16", "2p24GeV_Ta_X17": "g8"}   # CONTROL_GATE
P5 = {"0p73GeV_H2": (3.365, 3.397), "2p24GeV_H2": (1.847, 1.858), "3p49GeV_H2": (1.538, 1.561), "3p49GeV_H2_old": (1.509, 1.532),
      "2p24GeV_Ta_X17": (2.004, 2.019)}   # P5_budget.md §1: best gate + log weight (in, oos)
rows = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def mse(v):
    v = np.array([x for x in v if math.isfinite(x)])
    return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size))) if v.size > 1 else (float("nan"),) * 2


out = []
for rt in RTS:
    S = {r["path"]: r for r in rows(O / rt / "summary.tsv")}
    R = [r for r in rows(O / rt / "regions_splits.tsv") if r["path"] == GATE[rt] and r["kind"] == "ep" and r["region"] == "pool"]
    s = S[GATE[rt]]
    row = dict(run_type=rt, gate=GATE[rt], n_splits=s["n_splits"], n_used=s["n_used"])
    for q in ("res_cb", "res_cb_oos"):
        row[q], row[q + "_se"] = float(s[q]), float(s[q + "_se"])
    for q in ("tail_hi3", "tail_lo_f", "leak_cut", "eff_cut"):
        row[q], row[q + "_se"] = mse([float(r[q]) for r in R])
    row["pu_pct"], row["pu_se"] = float(s["pu"]), float(s["pu_se"])
    row["yield_pct"] = float(s["yield_pct"])
    out.append(row)
    print(f"{rt}: scan res_cb_oos " + " ".join(f"{p}={float(r['res_cb_oos']):.4f}" for p, r in S.items()))
    print(f"   control {GATE[rt]} {row['res_cb']:.4f}+-{row['res_cb_se']:.4f} / oos {row['res_cb_oos']:.4f}; P5 {P5[rt]}; n_pool rows {len(R)}")
with open(IMP / "control_phase2.tsv", "w") as f:
    f.write("\t".join(out[0]) + "\n")
    f.writelines("\t".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in r.values()) + "\n" for r in out)
