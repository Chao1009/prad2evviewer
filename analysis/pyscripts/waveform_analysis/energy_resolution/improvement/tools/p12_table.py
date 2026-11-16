#!/usr/bin/env python3
"""P12 loss table: reads variants/P12/{raw,m1,m8,m12,m16}; writes variants/P12/p12_loss_vs_cut.tsv and p12_loss_fit.tsv and prints them.
Model: effective threshold t = max(cut, T_HW) with T_HW the median readout edge (MeV); res^2(t) = r0^2 + b t^2 and
mu_none(t) = mu0 - s t, fitted on the points t = T_HW (cuts 1 and 8 pooled), 12, 16; the loss at the readout edge is b T_HW^2 and s T_HW."""
import csv
import sys
from pathlib import Path
import numpy as np

V = Path(__file__).resolve().parents[4] / "output/energy_resolution/variants/P12"
ENERGY = {"0p73GeV_H2": ("0p73", 728.0), "2p24GeV_H2": ("2p24", 2240.0), "3p49GeV_H2": ("3p49", 3490.0)}


def rd(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def pick(rows, **kw):
    return next(r for r in rows if all(r[k] == v for k, v in kw.items()))


def main():
    out = [("run_type", "cut_mev", "res_cb", "res_cb_edge", "res_cb_centre", "tail_lo3", "tail_lo3_edge", "tail_lo3_centre", "mu_none_pct", "mu_none_se")]
    fits = [("run_type", "T_hw_mev", "r0_pct", "b_pct2_per_mev2", "res_loss_at_T_pct", "loss_frac_of_res_quad", "mu_slope_pct_per_mev", "mean_loss_at_T_pct",
             "mean_loss_at_T_mev", "res_cb_ctrl", "loss_frac_if_linear_in_t", "loss_frac_if_t3")]
    for rt, (tag, E) in ENERGY.items():
        thw = float(np.median(10 * np.load(V / f"raw/edge_{rt}.npz")["k_mev_per_adc"][np.isfinite(np.load(V / f"raw/edge_{rt}.npz")["k_mev_per_adc"])]))
        pts = {}
        for c in (1, 8, 12, 16):
            rg, sp = rd(V / f"m{c}/{rt}/regions.tsv"), [r for r in rd(V / f"m{c}/{rt}/splits.tsv") if r["path"] == "legacy" and r["variant"] == "all" and r["mu_ep_none"] != "nan"]  # none arm: legacy only
            g = lambda reg, q: float(pick(rg, path="g12", kind="ep", region=reg)[q])   # noqa: E731
            m = np.array([float(r["mu_ep_none"]) for r in sp])
            pts[c] = dict(res=g("pool", "res_cb"), t=g("pool", "tail_lo3"), mu=m.mean(), mse=m.std(ddof=1) / np.sqrt(len(m)))
            out.append((rt, c, pts[c]["res"], g("pool_edge", "res_cb"), g("pool_centre", "res_cb"), pts[c]["t"], g("pool_edge", "tail_lo3"),
                        g("pool_centre", "tail_lo3"), pts[c]["mu"], pts[c]["mse"]))
        x = np.array([thw, 12.0, 16.0])
        r2 = np.array([np.mean([pts[1]["res"] ** 2, pts[8]["res"] ** 2]), pts[12]["res"] ** 2, pts[16]["res"] ** 2])
        mu = np.array([np.mean([pts[1]["mu"], pts[8]["mu"]]), pts[12]["mu"], pts[16]["mu"]])
        b, a = np.polyfit(x ** 2, r2, 1)
        s, mu0 = np.polyfit(x, mu, 1)
        loss = np.sqrt(max(b, 0.0) * thw ** 2)
        alt = []
        for p in (1, 3):                      # sensitivity of the extrapolation to the power of t in the variance
            bp = np.polyfit(x ** p, r2, 1)[0]
            alt.append(np.sqrt(max(bp, 0.0) * thw ** p) / np.sqrt(r2[0]))
        fits.append((rt, thw, np.sqrt(max(a, 0.0)), b, loss, loss / np.sqrt(np.mean(r2[:1])), -s, -s * thw, -s * thw * E / 100.0, np.sqrt(r2[0]), *alt))
    for name, tab in (("p12_loss_vs_cut.tsv", out), ("p12_loss_fit.tsv", fits)):
        with open(V / name, "w") as f:
            for r in tab:
                f.write("\t".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in r) + "\n")
        print(name)
        for r in tab:
            print("\t".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in r))


if __name__ == "__main__":
    main()
