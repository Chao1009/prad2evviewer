#!/usr/bin/env python3
"""profile_mfchi2_vs_dtcfd.py - median MF chi2/dof (curve) with the 16-84 % interval (band) vs calibrated CFD-fit dt, for 0.73/2.24/3.49 GeV H2,
single-cluster > 100 MeV (cache/<slug>_e100).  Writes figures/summary/mfchi2_profile_vs_dtcfd.{png,pdf} and tables/mfchi2_profile_vs_dtcfd.tsv.
Run from analysis/pyscripts after setup_env.sh."""
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import summary_dtcfd_vs_mfchi2 as S                    # noqa: E402

CT = (7.0, 15.0)
NM = int(sys.argv[sys.argv.index("--nmod") + 1]) if "--nmod" in sys.argv else None     # --nmod 20: only events with n_mod < 20
SFX = f"_nmod{NM}" if NM else ""
COL = {"0p73GeV_H2": "C0", "2p24GeV_H2": "C1", "3p49GeV_H2": "C3"}


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    edges = np.arange(-5, 25.01, 0.5)
    ctr = 0.5 * (edges[1:] + edges[:-1])
    fig, axs = plt.subplots(1, 2, figsize=(13, 5.2), sharex=True)
    tsv = ["energy\tdt_centre_ns\tN\tmedian\tp16\tp84"]
    for slug, lab in S.TYPES:
        dt, chi, _ = S.load(slug, slug + "_e100", "single", NM)
        k = np.digitize(dt, edges) - 1
        o = np.argsort(k, kind="stable"); ks, cs = k[o], chi[o]
        st = np.searchsorted(ks, np.arange(len(ctr))); en = np.searchsorted(ks, np.arange(len(ctr)), "right")
        x, m, lo, hi = [], [], [], []
        for i in range(len(ctr)):
            n = en[i] - st[i]
            if n < 100:
                continue
            q = np.percentile(cs[st[i]:en[i]], [16, 50, 84])
            x.append(ctr[i]); lo.append(q[0]); m.append(q[1]); hi.append(q[2])
            tsv.append(f"{lab}\t{ctr[i]:.2f}\t{n}\t{q[1]:.3f}\t{q[0]:.3f}\t{q[2]:.3f}")
        w = (dt >= CT[0]) & (dt <= CT[1]); c2 = np.abs(dt) < 2
        qw = np.percentile(chi[w], [16, 50, 84]); qc = np.percentile(chi[c2], [16, 50, 84])
        print(f"{lab}: N = {len(dt):,}; core |dt|<2: median {qc[1]:.3f}, p16 {qc[0]:.3f}, p84 {qc[2]:.3f} (N={int(c2.sum()):,}); "
              f"window 7..15: median {qw[1]:.3f}, p16 {qw[0]:.3f}, p84 {qw[2]:.3f} (N={int(w.sum()):,})")
        for ax in axs:
            ax.fill_between(x, lo, hi, color=COL[slug], alpha=0.25, lw=0)
            ax.plot(x, m, color=COL[slug], lw=1.8, label=f"{lab}, N = {len(dt):,}")
    for ax, lg in zip(axs, (False, True)):
        ax.axvspan(CT[0], CT[1], color="orange", alpha=0.12)
        ax.set_xlim(-5, 25); ax.set_xlabel("dt = T_module - T_seed, CFD log-normal fit  [ns]")
        ax.set_ylabel("MF chi2/dof: median and 16-84 % band" + (" (log)" if lg else ""))
        if lg:
            ax.set_yscale("log")
        else:
            ax.set_ylim(0, 20)
    axs[0].legend(loc="upper left", fontsize=9)
    fig.suptitle(f"single-cluster events > 100 MeV{' with n_mod < ' + str(NM) if NM else ''}, calibrated CFD fit-only dt; curve = median, band = 16-84 % interval (not a Gaussian sigma); bins with < 100 hits skipped", fontsize=9)
    fig.tight_layout()
    out = HERE / "figures" / "summary"; out.mkdir(parents=True, exist_ok=True)
    fig.savefig(out / f"mfchi2_profile_vs_dtcfd{SFX}.png", dpi=150); fig.savefig(out / f"mfchi2_profile_vs_dtcfd{SFX}.pdf")
    (HERE / "tables").mkdir(exist_ok=True)
    (HERE / "tables" / f"mfchi2_profile_vs_dtcfd{SFX}.tsv").write_text("\n".join(tsv) + "\n")


if __name__ == "__main__":
    main()
