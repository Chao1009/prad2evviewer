#!/usr/bin/env python3
"""Fit R(E) = sigma_E/E = sqrt(a^2 + b^2/E + c^2/E^2) (a, b in %, c in % GeV) to the best-gate out-of-sample Crystal Ball
resolution of three H2 run types; the tantalum X17 point is drawn only.  Reads energy_resolution/<run type>/summary.tsv, writes
energy_resolution/resolution_fit.{pdf,png,tsv}.  Also fits the two-term R(E) = sqrt(a^2 + b^2/E) to the same points.  Usage: energy_resolution_fit.py [directory of energy_resolution]"""
import csv
import json
import os
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit

H2 = {"0p73GeV_H2": 24590, "2p24GeV_H2": 25228, "3p49GeV_H2": 24698}      # run type: one of its runs (beam energy from runinfo)
X17 = ("2p24GeV_Ta_X17", 25800)
YIELD_MIN = 99.0                                                            # percent of legacy yield, as in energy_resolution_benchmark.py
# x is the beam energy; the measured E'/E_beam = 0.997-0.999, so using E' instead would not matter.
plt.rcParams.update({"axes.labelsize": 16, "axes.labelweight": "bold", "xtick.labelsize": 14, "ytick.labelsize": 14,
                     "legend.fontsize": 14, "lines.linewidth": 2.0})


def beam_gev(run):
    j = json.load(open(Path(os.environ["PRAD2_DATABASE_DIR"]) / "runinfo" / "general.json"))
    cs = sorted((c for c in j["configurations"] if "from_run" in c and c["from_run"] <= run), key=lambda c: c["from_run"])
    return [c["beam_energy"] for c in cs if "beam_energy" in c][-1] / 1000.0


def point(root, rt, run):
    """(beam energy GeV, best-gate res_cb_oos %, standard error, gate name)."""
    with open(root / rt / "summary.tsv") as f:
        S = {r["path"]: r for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t")}
    v = lambda p, k: float(S[p][k])                                                      # noqa: E731
    ok = sorted((float(p[1:]), p) for p in S if p != "legacy" and v(p, "yield_pct") >= YIELD_MIN
                and np.isfinite(v(p, "res_cb_oos")) and np.isfinite(v(p, "res_cb_oos_se")))
    lo = min(ok, key=lambda x: v(x[1], "res_cb_oos"))[1]
    p = next(p for g, p in ok if v(p, "res_cb_oos") <= v(lo, "res_cb_oos") + v(lo, "res_cb_oos_se"))
    return beam_gev(run), v(p, "res_cb_oos"), v(p, "res_cb_oos_se"), p


def model(E, a, b, c):
    return np.sqrt(a**2 + b**2 / E + c**2 / E**2)


def model2(E, a, b):
    return np.sqrt(a**2 + b**2 / E)


def main():
    root = Path(sys.argv[1]) if len(sys.argv) > 1 else Path(__file__).resolve().parent / "energy_resolution"
    E, R, s, gate = zip(*(point(root, rt, run) for rt, run in H2.items()))
    q, cov = curve_fit(model, E, R, p0=[1, 1, 1], sigma=s, absolute_sigma=True, bounds=(0, np.inf), maxfev=20000)
    err = np.sqrt(np.diag(cov))
    q2, cov2 = curve_fit(model2, E, R, p0=[1, 1], sigma=s, absolute_sigma=True, bounds=(0, np.inf), maxfev=20000)
    err2 = np.sqrt(np.diag(cov2))
    chi2 = np.sum(((np.asarray(R) - model2(np.asarray(E), *q2)) / np.asarray(s)) ** 2)
    dof2 = len(E) - len(q2)
    ex, rx, sx, gx = point(root, *X17)
    print("gates:", gate, gx, "E:", E, ex, "se %:", s, sx)
    with open(root / "resolution_fit.tsv", "w") as f:
        f.write("a\ta_err\tb\tb_err\tc\tc_err\n" + "\t".join(f"{x:.4g}" for pair in zip(q, err) for x in pair) + "\n")
        f.write("\t".join([f"{x:.4g}" for pair in zip(q2, err2) for x in pair] + ["", ""]) + "\n")
    fig, ax = plt.subplots(figsize=(9, 7))
    Eg = np.linspace(0.5, 3.8, 300)
    ax.plot(Eg, model(Eg, *q), color="C3", label=r"$a\oplus b/\sqrt{E}\oplus c/E$")
    ax.plot(Eg, model2(Eg, *q2), color="C2", ls="--", label=r"$a\oplus b/\sqrt{E}$")
    ax.errorbar(E, R, s, color="C0", marker="o", ms=6, ls="", capsize=4, elinewidth=1.5, label="hydrogen, best gate")
    ax.errorbar(ex, rx, sx, color="C1", marker="*", ms=12, mfc="none", ls="", capsize=4, elinewidth=1.5, label="2.24 GeV X17, not fitted")
    ax.text(0.97, 0.97, f"three-term, zero degrees of freedom:\na = ({q[0]:.2f} ± {err[0]:.2f}) %\n"
            f"b = ({q[1]:.2f} ± {err[1]:.2f}) % $\\sqrt{{\\mathrm{{GeV}}}}$\nc = ({q[2]:.2f} ± {err[2]:.2f}) % GeV\n"
            f"two-term, one degree of freedom:\n$\\chi^2$/dof = {chi2:.2f}/{dof2}\n"
            f"a = ({q2[0]:.2f} ± {err2[0]:.2f}) %\nb = ({q2[1]:.2f} ± {err2[1]:.2f}) % $\\sqrt{{\\mathrm{{GeV}}}}$",
            transform=ax.transAxes, ha="right", va="top", fontsize=14, bbox=dict(fc="white", ec="gray"))
    ax.set_xlabel("beam energy (GeV)")
    ax.set_ylabel(r"$\sigma_E/E$ (%)")
    ax.set_ylim(0, 4.5)
    ax.legend(loc="lower left")
    ax.text(0.97, 0.03, f"error bars: standard error over splits ({min(s + (sx,)):.3f}–{max(s + (sx,)):.3f} %)",
            transform=ax.transAxes, ha="right", va="bottom", fontsize=12, color="dimgray")
    fig.tight_layout()
    fig.savefig(root / "resolution_fit.pdf")
    fig.savefig(root / "resolution_fit.png", dpi=100)


if __name__ == "__main__":
    main()
