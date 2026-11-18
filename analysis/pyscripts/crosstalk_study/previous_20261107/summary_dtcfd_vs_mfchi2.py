#!/usr/bin/env python3
"""summary_dtcfd_vs_mfchi2.py - one figure: calibrated CFD log-normal-fit dt = T_module - T_seed (x) vs MF chi2/dof of the non-seed module (y).
Rows: e-p (cache/<slug>) and single-cluster 100 MeV (cache/<slug>_e100); columns 0.73, 2.24, 3.49 GeV H2.  From the existing caches only.
Writes figures/summary/dtcfd_vs_mfchi2.{png,pdf} and tables/dtcfd_vs_mfchi2.md.   Run from analysis/pyscripts after setup_env.sh."""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import crosstalk_study as C                      # noqa: E402

TYPES = [("0p73GeV_H2", "0.73 GeV H2"), ("2p24GeV_H2", "2.24 GeV H2"), ("3p49GeV_H2", "3.49 GeV H2")]
CT = (7.0, 15.0)


def load(slug, cache_type, sel, nmod=None):
    import glob
    fs = sorted(glob.glob(str(HERE / "cache" / cache_type / "*_legacy.npz")))
    P = [dict(np.load(f)) for f in fs]
    names = P[0]["module_names"]
    d = {k: np.concatenate([p[k] for p in P]) for k in P[0] if P[0][k].ndim and k != "module_names"}
    nev = sum(int(p["n_clusters"]) for p in P)
    if nmod:                                         # events with fewer than `nmod` fired modules (cache column nmod_wa20, see crosstalk_study.add_nmod)
        hk, ck = d["nmod_wa20"] < nmod, d["cl_nmod_wa20"] < nmod
        d = {k: (v[hk] if len(v) == len(hk) else v[ck] if len(v) == len(ck) else v) for k, v in d.items()}
        nev = int(ck.sum())
    if sel == "ep":
        sig = sum(float(p["sigma_cb"]) * int(p["n_inner"]) for p in P) / sum(int(p["n_inner"]) for p in P)
        hm = (d["moller"] == 0) & (np.abs(d["r"] - 1) < 3 * sig)
        nev = int(((d["cl_moller"] == 0) & (np.abs(d["cl_r"] - 1) < 3 * sig)).sum())
        d = {k: (v[hm] if len(v) == len(hm) else v) for k, v in d.items()}
    ok = (d["algo"] == 2) & (d["seed_algo"] == 2) & np.isfinite(d["t_wa"]) & np.isfinite(d["chi2"])
    off, _ = C.load_offsets("auto", names)
    dt = (d["t_wa"] - off[d["mod_i"]]) - (d["seed_t_wa"] - off[d["seed_mod_i"]])
    return dt[ok], d["chi2"][ok], nev


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    nm = int(sys.argv[sys.argv.index("--nmod") + 1]) if "--nmod" in sys.argv else None     # --nmod 20: events with n_mod < 20, both rows from cache/<slug>_e100
    sfx = f"_nmod{nm}" if nm else ""
    data = {}
    for slug, _ in TYPES:
        data[(slug, "ep")] = load(slug, slug + "_e100", "ep", nm) if nm else load(slug, slug, "ep")
        data[(slug, "single")] = load(slug, slug + "_e100", "single", nm)
    zb = np.arange(-5, 25.01, 0.25)
    allc = np.concatenate([c[(dt >= zb[0]) & (dt <= zb[-1])] for dt, c, _ in data.values()])
    ymax = float(np.ceil(np.percentile(allc, 99)))
    yb = np.linspace(0, ymax, 60)
    fig, axs = plt.subplots(2, 3, figsize=(15, 8.5), sharex=True, sharey=True)
    rows = []
    for ri, sel in enumerate(("ep", "single")):
        for ci, (slug, lab) in enumerate(TYPES):
            dt, chi, nev = data[(slug, sel)]
            ax = axs[ri, ci]
            h = ax.hist2d(dt, np.clip(chi, 0, ymax), bins=[zb, yb], norm=LogNorm(vmin=1), cmap="viridis")
            ax.axvspan(CT[0], CT[1], color="orange", alpha=0.15)
            ax.set_xlim(zb[0], zb[-1]); ax.set_ylim(0, ymax)
            ax.set_title(f"{lab}, {'e-p' if sel == 'ep' else 'single > 100 MeV'}\nN = {len(dt):,}", fontsize=10)
            if ri == 1:
                ax.set_xlabel("dt = T_module - T_seed, CFD fit  [ns]")
            if ci == 0:
                ax.set_ylabel(f"MF chi2/dof (linear, clipped at {ymax:.0f})")
            w = (dt >= CT[0]) & (dt <= CT[1]); c2 = np.abs(dt) < 2
            rows.append((lab, "e-p" if sel == "ep" else "single 100 MeV", nev, len(dt), int(w.sum()), np.median(chi[w]), np.median(chi[c2]),
                         tuple(np.mean(chi[w] > t) for t in (2, 3, 5, 10)) + tuple(np.mean(chi[c2] > t) for t in (2, 3, 5, 10)) + (int(c2.sum()),)))
    fig.colorbar(h[3], ax=axs, label="hits", shrink=0.8)
    out = HERE / "figures" / "summary"; out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"dtcfd_vs_mfchi2{sfx}.png", dpi=150, bbox_inches="tight"); fig.savefig(out / f"dtcfd_vs_mfchi2{sfx}.pdf", bbox_inches="tight")
    L = [f"# dt (calibrated CFD log-normal fit) vs MF chi2/dof; common y range 0..{ymax:.0f} (99 % quantile of all six panels in -5..25 ns)", "",
         "Window +7..+15 ns, core |dt| < 2 ns; hits = non-seed hits with time_algo 2 at seed and module; N_ev = events of the selection.", "",
         "Core fractions = false-rejection rate if chi2/dof were used as a crosstalk cut.", "",
         "| sample | selection | N_ev | N hits | window hits | core hits | median window | median core | window > 2 | > 3 | > 5 | > 10 | core > 2 | > 3 | > 5 | > 10 |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for r in rows:
        f = r[7]
        L.append(f"| {r[0]} | {r[1]} | {r[2]:,} | {r[3]:,} | {r[4]:,} | {f[8]:,} | {r[5]:.3f} | {r[6]:.3f} | "
                 + " | ".join(f"{100 * x:.1f} %" for x in f[:8]) + " |")
    (HERE / "tables").mkdir(exist_ok=True)
    (HERE / "tables" / f"dtcfd_vs_mfchi2{sfx}.md").write_text("\n".join(L) + "\n")
    print("\n".join(L))


if __name__ == "__main__":
    main()
