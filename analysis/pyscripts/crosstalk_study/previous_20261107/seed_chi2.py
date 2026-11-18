#!/usr/bin/env python3
"""seed_chi2.py - MF chi2/dof of the seed module (its own method-1 pulse of the benchmark pulse dump, the pulse that gives seed_t / seed_amp) added as hit-level
cache column `seed_chi2` of cache/<slug>_e100, and figures fig5_seedchi2_vs_amp / fig5b (log y) / fig5c_seedchi2_vs_dt for the e-p samples (single-cluster > 100 MeV, no Moller tag,
|r-1| < 3 sigma, calibrated CFD fit-only dt).  Writes tables/seed_chi2.md.  Run from analysis/pyscripts after setup_env.sh.  Usage: seed_chi2.py [slug ...]"""
import glob
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
SCR = Path("/private/tmp/claude-505/-Users-L00338853-work-PRad-prad2evviewer-analysis-pyscripts/db0289d2-4ae1-47bc-9d5f-4788ee635668/scratchpad")
sys.path.insert(0, str(HERE))
import crosstalk_study as C                          # noqa: E402

SAMPLES = [("0p73GeV_H2", "0.73 GeV H2"), ("2p24GeV_H2", "2.24 GeV H2"), ("3p49GeV_H2", "3.49 GeV H2"),
           ("2p24GeV_Ta_X17_25800", "Ta X17 25800"), ("2p24GeV_Ta_X17_26138", "Ta X17 26138"), ("2p24GeV_Ta_X17_26230", "Ta X17 26230")]
CT = (7.0, 15.0)


def add_seed_chi2(cachef, pulsef):
    z = dict(np.load(cachef))
    pz = np.load(pulsef, allow_pickle=True)
    k1 = pz["method"] == 1
    ev0 = int(z["ev"].min())
    key = ((pz["event"][k1].astype(np.int64) - ev0) * 4096 + pz["midx"][k1]).astype(np.float64) * 1000.0 + pz["time_ns"][k1]
    o = np.argsort(key); kk, ch = key[o], pz["chi2dof"][k1][o]
    q = ((z["ev"].astype(np.int64) - ev0) * 4096 + z["seed_mod_i"]).astype(np.float64) * 1000.0 + z["seed_t"]
    j = np.clip(np.searchsorted(kk, q), 1, len(kk) - 1)
    b = np.where(np.abs(kk[j - 1] - q) <= np.abs(kk[j] - q), j - 1, j)
    ok = np.abs(kk[b] - q) < 0.01
    z["seed_chi2"] = np.where(ok, ch[b], np.nan).astype(np.float32)
    np.savez_compressed(cachef, **z)
    return z


def load(slug):
    cdir = HERE / "cache" / f"{slug}_e100"
    fs = sorted(glob.glob(str(cdir / "*_legacy.npz")))
    P, nmiss, ntot = [], 0, 0
    for f in fs:
        z = dict(np.load(f))
        if "seed_chi2" not in z:
            stem = Path(f).name.split("_legacy")[0]
            z = add_seed_chi2(f, SCR / f"ct_{slug}_e100" / stem / "pulse.npz")
        nmiss += int(np.isnan(z["seed_chi2"]).sum()); ntot += len(z["seed_chi2"])
        P.append(z)
    names = P[0]["module_names"]
    d = {k: np.concatenate([p[k] for p in P]) for k in P[0] if P[0][k].ndim and k != "module_names"}
    d["part"] = np.repeat(np.arange(len(P), dtype=np.int32), [len(p["t"]) for p in P])
    sig = sum(float(p["sigma_cb"]) * int(p["n_inner"]) for p in P) / sum(int(p["n_inner"]) for p in P)
    hm = (d["moller"] == 0) & (np.abs(d["r"] - 1) < 3 * sig)
    nev = int(((d["cl_moller"] == 0) & (np.abs(d["cl_r"] - 1) < 3 * sig)).sum())
    d = {k: (v[hm] if len(v) == len(hm) else v) for k, v in d.items()}
    ok = (d["algo"] == 2) & (d["seed_algo"] == 2) & np.isfinite(d["t_wa"])
    off, _ = C.load_offsets("auto", names)
    dt = (d["t_wa"] - off[d["mod_i"]]) - (d["seed_t_wa"] - off[d["seed_mod_i"]])
    ekey = d["part"].astype(np.int64) * (1 << 32) + d["ev"]
    return dict(dt=dt[ok], amp=d["amp"][ok], sc=d["seed_chi2"][ok], ekey=ekey[ok], nev=nev, nmiss=nmiss, ntot=ntot)


def pct(x):
    q = np.nanpercentile(x, [16, 50, 84]); return f"{q[1]:.3f} [{q[0]:.3f}, {q[2]:.3f}]"


def main():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    from scipy.stats import spearmanr
    want = sys.argv[1:] or [s for s, _ in SAMPLES]
    L = ["# Seed-module MF chi2/dof (e-p, single-cluster > 100 MeV, calibrated CFD fit-only hits)", "",
         "Seed chi2/dof = chi2/dof of the seed's own MF pulse (benchmark pulse dump). Event-level values: one per event with at least one selected non-seed hit. Hits: one per (event, non-seed member).", "",
         "| sample | events | hits | seeds without MF chi2 | seed chi2 events: median [p16, p84] | window hits +7..+15 ns | core hits abs(dt)<2 ns | events > 2 | > 3 | > 5 | Spearman(seed chi2, amp) | Spearman(seed chi2, dt) |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    ab = np.logspace(np.log10(5), np.log10(5000), 60)
    zb = np.arange(-5, 25.01, 0.25)
    for slug, lab in SAMPLES:
        if slug not in want:
            continue
        r = load(slug)
        dt, amp, sc, ek = r["dt"], r["amp"], r["sc"], r["ekey"]
        fin = np.isfinite(sc)
        dt, amp, sc, ek = dt[fin], amp[fin], sc[fin], ek[fin]
        _, first = np.unique(ek, return_index=True)
        ev_sc = sc[first]
        w = (dt >= CT[0]) & (dt <= CT[1]); c2 = np.abs(dt) < 2
        ymax = float(np.ceil(np.percentile(sc, 99)))
        out = HERE / "figures" / slug
        out.mkdir(parents=True, exist_ok=True)
        cap = f"{lab} e-p, N_ev = {r['nev']:,}, hits = {len(sc):,}\n(one entry per non-seed member: an event with 10 members enters 10 times at the same y)"
        for fn, ylog in (("fig5_seedchi2_vs_amp.png", False), ("fig5b_seedchi2_vs_amp_logy.png", True)):
            fig, ax = plt.subplots(figsize=(7.5, 5.5))
            yb = np.logspace(np.log10(max(sc[sc > 0].min(), 0.05)), np.log10(sc.max()), 60) if ylog else np.linspace(0, ymax, 60)
            h = ax.hist2d(np.clip(amp, ab[0], ab[-1]), np.clip(sc, yb[0], yb[-1]), bins=[ab, yb], norm=LogNorm(), cmap="viridis")
            fig.colorbar(h[3], ax=ax, label="hits")
            ax.set_xscale("log"); ax.set_xlim(ab[0], ab[-1])
            if ylog:
                ax.set_yscale("log")
            else:
                ax.set_ylim(0, ymax)
            ax.set_xlabel("non-seed module MF amplitude  [ADC]")
            ax.set_ylabel("seed MF chi2/dof" + (" (log)" if ylog else f" (linear, clipped at {ymax:.0f})"))
            ax.set_title(cap, fontsize=7.5); fig.tight_layout(); fig.savefig(out / fn, dpi=150); plt.close(fig)
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        h = ax.hist2d(dt, np.clip(sc, 0, ymax), bins=[zb, np.linspace(0, ymax, 60)], norm=LogNorm(), cmap="viridis")
        fig.colorbar(h[3], ax=ax, label="hits"); ax.axvspan(CT[0], CT[1], color="orange", alpha=0.15)
        ax.set_xlim(zb[0], zb[-1]); ax.set_ylim(0, ymax)
        ax.set_xlabel("dt = T_module - T_seed, CFD log-normal fit  [ns]"); ax.set_ylabel(f"seed MF chi2/dof (linear, clipped at {ymax:.0f})")
        ax.set_title(cap, fontsize=7.5); fig.tight_layout(); fig.savefig(out / "fig5c_seedchi2_vs_dt.png", dpi=150); plt.close(fig)
        rng = np.random.default_rng(1); ix = rng.choice(len(sc), min(len(sc), 1_000_000), replace=False)
        rho_a = spearmanr(sc[ix], amp[ix])[0]; rho_d = spearmanr(sc[ix], dt[ix])[0]
        L.append(f"| {lab} | {len(ev_sc):,} | {len(sc):,} | {r['nmiss']:,} of {r['ntot']:,} | {pct(ev_sc)} | {pct(sc[w])} (N={int(w.sum()):,}) | {pct(sc[c2])} (N={int(c2.sum()):,}) | "
                 f"{100 * np.mean(ev_sc > 2):.1f} % | {100 * np.mean(ev_sc > 3):.1f} % | {100 * np.mean(ev_sc > 5):.1f} % | {rho_a:+.3f} | {rho_d:+.3f} |")
        print(L[-1], flush=True)
    if len(want) == len(SAMPLES):
        (HERE / "tables").mkdir(exist_ok=True)
        (HERE / "tables" / "seed_chi2.md").write_text("\n".join(L) + "\n")
    else:
        print("\n".join(L))


if __name__ == "__main__":
    main()
