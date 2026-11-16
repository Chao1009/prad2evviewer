"""P8 table. Run from analysis/pyscripts: python3 waveform_analysis/energy_resolution/improvement/tools/p8_table.py
Paired deltas (variant - phase-2 control, per (run, split), mean +- SE) of the pool ep rows of regions_splits.tsv for every variant
under output/energy_resolution/variants/P8/<variant>/<RT>, at the legacy path and the control gate; yield loss from splits.tsv; Ta also
for run 26138 alone.  Writes improvement/P8_secondaries_acceptance.tsv and prints one line per row."""
import csv, math, sys
from pathlib import Path
import numpy as np

C = Path("output/energy_resolution/control_phase2_20261113")
V = Path("output/energy_resolution/variants/P8")
SCAN = Path("output/energy_resolution/scan")
IMP = Path("waveform_analysis/energy_resolution/improvement")
GATE = {"0p73GeV_H2": "g20", "2p24GeV_H2": "g16", "3p49GeV_H2": "g16", "3p49GeV_H2_old": "g16", "2p24GeV_Ta_X17": "g8"}   # CONTROL_GATE
Q = ("res_cb", "res_cb_oos", "sep_pi", "sep_pi_oos", "tail_hi3", "leak_cut", "leak_cut_oos", "tail_lo_f", "eff_cut", "mu")
rows = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def mse(v):
    v = np.array([x for x in v if math.isfinite(x)])
    return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size))) if v.size > 1 else (float("nan"),) * 2


def banks(rt):
    return ",".join(sorted({str(np.load(f)["bank_sha1"])[:8] for f in sorted((SCAN / rt).glob("run_*.npz")) if "24937" not in f.name}))


out = []
for vdir in sorted(p for p in V.iterdir() if p.is_dir() and p.name not in ("logs",)):
    for rt in GATE:
        d = vdir / rt
        if not (d / "regions_splits.tsv").exists():
            continue
        ovr = open(d / "regions_splits.tsv").readline().strip().removeprefix("# overrides:").strip()
        sel = lambda p, path: {(r["run"], r["split"]): r for r in rows(p / "regions_splits.tsv")
                               if r["path"] == path and r["kind"] == "ep" and r["region"] == "pool"}
        yl = lambda p, path: {(r["run"], r["split"]): (float(r["yield_n"]), float(r["pu"])) for r in rows(p / "splits.tsv")
                              if r["path"] == path and r["variant"] == "all"}
        for path in ("legacy", GATE[rt]):
            cv, vv = sel(C / rt, path), sel(d, path)
            yc, yv = yl(C / rt, path), yl(d, path)
            for runs in ((None, "26138") if rt.endswith("Ta_X17") else (None,)):
                ks = [k for k in cv if k in vv and (runs is None or k[0] == runs)]
                row = dict(variant=vdir.name, run_type=rt, runs=runs or "all", path=path, n_pairs=len(ks))
                for q in Q:
                    row[f"d_{q}"], row[f"d_{q}_se"] = mse([float(vv[k][q]) - float(cv[k][q]) for k in ks])
                    row[f"v_{q}"] = "real" if abs(row[f"d_{q}"]) >= 2 * row[f"d_{q}_se"] else "null"
                ky = [k for k in ks if k in yc and k in yv]
                row["yield_loss_pct"] = 100 * (1 - sum(yv[k][0] for k in ky) / sum(yc[k][0] for k in ky)) if ky else float("nan")
                row["pu_ctl"], row["pu_var"] = (float(np.mean([y[k][1] for k in ky])) for y in (yc, yv))
                row["bank_sha1"], row["overrides"] = banks(rt), ovr
                out.append(row)
                f2 = lambda q: f"{row['d_' + q]:+.4f}±{row['d_' + q + '_se']:.4f}"
                print(f"{vdir.name:22s} {rt:15s} {row['runs']:5s} {path:6s} n{len(ks):3d} res_cb {f2('res_cb')} oos {f2('res_cb_oos')} "
                      f"sep {f2('sep_pi')} hi3 {f2('tail_hi3')} leak {f2('leak_cut')} lo_f {f2('tail_lo_f')} eff {f2('eff_cut')} "
                      f"yloss {row['yield_loss_pct']:+.2f} pu {row['pu_ctl']:.1f}->{row['pu_var']:.1f}")
if out:
    with open(IMP / "P8_secondaries_acceptance.tsv", "w") as f:
        f.write("\t".join(out[0]) + "\n")
        f.writelines("\t".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in r.values()) + "\n" for r in out)
