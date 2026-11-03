#!/usr/bin/env python3
"""accidental_comb_period.py - bunch period from the no-gate accidental comb (recreated 10-30).

Input  module_timing_nogate_pooled.npz written by plot_module_timing.py --veto-gate none
       (per-module all-pairs t_module - t_veto histograms per run type, 0.5 ns bins; modules are
       aligned on their own in-time peak before pooling).
Method teeth k = +-2 ... +-24 (k = 0, +-1 hold the in-time signal) at c + k P around the main
       peak; per tooth a Gaussian (sigma fixed 0.55 ns) + constant fit in +-1.5 ns, centroid = fitted mean; then a
       weighted straight line centroid = k P + c (tooth errors from the fit; sigma_P scaled by
       sqrt(chi2/ndf) when > 1).  P_RF = 128.247 / 32 ns.
Usage  python3 accidental_comb_period.py NPZ --out-prefix PATH   (writes PATH.png and PATH.txt)
"""
import argparse
import numpy as np
from scipy.optimize import curve_fit

P_RF = 128.247 / 32
HW = 1.5


SIG = 0.55   # tooth sigma fixed (module + veto + bunch length, 10-25 note)


def g(x, a, m, b):
    return a * np.exp(-0.5 * ((x - m) / SIG) ** 2) + b


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("npz")
    ap.add_argument("--out-prefix", required=True)
    a = ap.parse_args()
    d = np.load(a.npz)
    x = d["x"]
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    tags = [str(t) for t in d["tags"]]
    fig, axs = plt.subplots(len(tags), 2, figsize=(12, 3.2 * len(tags)), squeeze=False)
    out = ["Bunch period from the no-gate accidental comb (t_module - t_veto, all pairs)",
           f"teeth k = +-2..+-24, Gaussian(sigma 0.55)+const fit +-{HW} ns per tooth, weighted fit centroid = kP + c; "
           f"P_RF = 128.247/32 = {P_RF:.5f} ns", ""]
    for r, t in enumerate(tags):
        d0 = float(d["d0"][r])
        # align every module on its own in-time peak (module time offsets smear the teeth when pooled)
        hm = d[f"hm_{t}"].astype(float)
        h = np.zeros(hm.shape[1])
        cen = len(x) // 2
        for row in hm:
            if row.sum() < 200:
                continue
            smr = np.convolve(row, np.ones(3) / 3, "same")
            mk = (x > -d0 - 15) & (x < -d0 + 15)
            pk = x[mk][np.argmax(smr[mk])]
            sh = int(round((pk - (-d0)) / (x[1] - x[0])))
            h += np.roll(row, -sh)
        m = (x > -d0 - 15) & (x < -d0 + 15)
        sm = np.convolve(h, np.ones(3) / 3, "same")
        c0 = x[m][np.argmax(sm[m])]
        P0 = P_RF
        ks, mu, er = [], [], []
        for k in list(range(-24, -1)) + list(range(2, 25)):
            e = c0 + k * P0
            w = np.abs(x - e) < HW + 0.01
            if w.sum() < 5 or h[w].sum() < 50:
                continue
            try:
                p, cv = curve_fit(g, x[w], h[w], p0=[h[w].max() - h[w].min(), e, h[w].min()],
                                  sigma=np.sqrt(np.maximum(h[w], 1)), absolute_sigma=True, maxfev=4000)
            except Exception:
                continue
            if p[0] <= 0 or abs(p[1] - e) > HW or not np.isfinite(cv[1, 1]):
                continue
            ks.append(k); mu.append(p[1]); er.append(np.sqrt(cv[1, 1]))
        ks, mu, er = map(np.array, (ks, mu, er))
        wt = 1 / er ** 2
        A = np.vstack([ks, np.ones_like(ks)]).T
        cov = np.linalg.inv(A.T @ (A * wt[:, None]))
        P, c = cov @ (A.T @ (wt * mu))
        chi2 = float(np.sum(wt * (mu - A @ [P, c]) ** 2))
        ndf = len(ks) - 2
        sP = np.sqrt(cov[0, 0]) * max(1.0, np.sqrt(chi2 / ndf))
        res = mu - (P * ks + c)
        out.append(f"{t}: teeth {len(ks)}  P = {P:.5f} +- {sP:.5f} ns  (P-P_RF)/sigma = {(P-P_RF)/sP:+.1f}  "
                   f"(P-4.000)/sigma = {(P-4.0)/sP:+.1f}  chi2/ndf = {chi2:.1f}/{ndf}  c = {c:.3f}  "
                   f"centroid rms resid = {res.std():.3f} ns")
        ax = axs[r, 0]
        mm = np.abs(x - c0) < 100
        ax.step(x[mm], h[mm], where="mid", lw=0.7)
        ax.set_yscale("log"); ax.set_title(f"{t}: all-pairs t_module - t_veto", fontsize=9)
        ax.set_xlabel("ns")
        ax = axs[r, 1]
        ax.errorbar(ks, res * 1e3, er * 1e3, fmt="o", ms=3)
        ax.axhline(0, c="k", lw=0.5)
        ax.set_title(f"P = {P:.5f} +- {sP:.5f} ns ({len(ks)} teeth); RF {P_RF:.5f}", fontsize=9)
        ax.set_xlabel("tooth k"); ax.set_ylabel("centroid - (kP + c)  [ps]")
    fig.tight_layout()
    fig.savefig(a.out_prefix + ".png", dpi=110)
    open(a.out_prefix + ".txt", "w").write("\n".join(out) + "\n")
    print("\n".join(out))


if __name__ == "__main__":
    main()
