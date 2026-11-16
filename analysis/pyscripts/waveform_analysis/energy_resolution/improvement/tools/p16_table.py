#!/usr/bin/env python3
"""P16 slice 2: below-threshold (KTeV-style) table per run type from a control report's thr_clusters.npz (THR_DUMP=True).

python p16_table.py <control root> <out dir> [--rts ...] [--half even|odd]

Entries: the usable PbWO4 window modules (5 x 5 around the seed, seed excluded; threshold known) of the peak clusters (|E / E' - 1| < 2 s,
s the iterated core width), binned in d (distance of the module centre from the cluster position, module sizes, DB) and in the module threshold
T_m = TET x k_m (raw MeV, classes split at T_EDGES).  Units x = e / T_m (raw member energy over raw threshold: gain independent).
Per bin: P_a = share without a stored hit; the true deposits below X_STAR (full efficiency above) are taken as a truncated exponential on
[0, X_STAR] whose mass is P_a + the share of read-out members below X_STAR and whose edge density is the members' density at X_STAR
(exponential fit on [X_STAR, X_STAR + 0.6]); expected energy of an absent module c = (mass x truncated mean - read-out energy below
X_STAR) / P_a, clipped to [0, X_STAR]; 0 where P_a < PA_MIN or > PA_MAX.  The report adds c x T_m per absent module (BELOW_THR).  --half: only the even / odd split ordinal
of each run (overfitting check, U2).  Prints per bin P_a and c, the mean added energy per cluster and sum P_a (1 - P_a) c^2 (the variance a
presence-based correction can remove, independent modules)."""
import argparse, math
from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import p16_common as pc
bm = pc.bm

DB = np.arange(0.5, 3.51, 0.125)
T_EDGES = np.array([9.0, 9.8])
X_STAR, X_FIT = 1.2, 0.6
MIN_N = 200
PA_MIN = 0.005                    # bins where almost every module is read out carry none either (c = L / P_a unstable)
PA_MAX = 0.95                     # bins where almost every module is absent carry no correction (the model's mean there is not constrained)


def trunc_mean(R):
    """Mean (units of X_STAR) of a density A exp(-t x / X_STAR) on [0, X_STAR] with edge density x X_STAR / mass = R."""
    h = lambda t: (t / math.expm1(t) if abs(t) > 1e-9 else 1.0) - R   # noqa: E731  decreasing in t, h(0) = 1 - R
    t = brentq(h, -60.0, 60.0) if h(60.0) < 0 < h(-60.0) else (60.0 if R <= 1 else -60.0)
    return 0.5 if abs(t) < 1e-6 else 1.0 / t - 1.0 / math.expm1(t)


def core(r, k=2.5, it=8):
    """Median and width (1.4826 MAD) of r iterated inside +- k x width from the mode and 3 %: the elastic core without the radiative and
    Moller tails."""
    h, e = np.histogram(r, np.arange(-0.3, 0.3001, 0.004))
    med, s = 0.5 * (e[h.argmax()] + e[h.argmax() + 1]), 0.03            # start at the mode with 3 %
    for _ in range(it):
        w = r[np.abs(r - med) < k * s]
        med, s = np.median(w), 1.4826 * np.median(np.abs(w - np.median(w)))
    return med, s


def edge_density(x, N):
    """Density per entry per unit x at X_STAR from the read-out x in [X_STAR, X_STAR + X_FIT] (log-linear fit of 6 bins)."""
    e = np.linspace(X_STAR, X_STAR + X_FIT, 7)
    h = np.histogram(x, e)[0].astype(float)
    if (h > 0).sum() < 3:
        return h.sum() / (X_FIT * N)
    c = 0.5 * (e[1:] + e[:-1])
    k = h > 0
    b = np.polyfit(c[k] - X_STAR, np.log(h[k] / (N * (e[1] - e[0]))), 1, w=np.sqrt(h[k]))
    return float(np.exp(b[1]))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("root")
    ap.add_argument("out")
    ap.add_argument("--rts", nargs="+", default=["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old"])
    ap.add_argument("--half", choices=("even", "odd"))
    a = ap.parse_args()
    Path(a.out).mkdir(parents=True, exist_ok=True)
    for rt in a.rts:
        z, geo = pc.load(f"{a.root}/{rt}/thr_clusters.npz"), pc.geometry(rt)
        T = pc.thresholds(rt)
        ep = bm.e_elastic(bm.theta_deg(geo["x"], geo["y"], geo["hz"]), geo["ebeam"], bm.pool_mass(rt))[z["seed"].astype(np.int64)]   # target mass of the run type (Ta: v3)
        r = z["e"] / ep - 1
        med, s = core(r[z["used"]])
        pk = z["used"] & (np.abs(r - med) < 2 * s)
        if a.half:
            od = {}
            for run in np.unique(z["run"]):
                sp = sorted(set(z["split"][z["run"] == run].tolist()))
                od |= {(int(run), x): i % 2 for i, x in enumerate(sp)}
            par = np.array([od[(int(u), int(v))] for u, v in zip(z["run"], z["split"])])
            pk &= par == (0 if a.half == "even" else 1)
        m = z["win_m"].astype(np.int64)
        mm = np.where(m >= 0, m, 0)
        usable = np.isfinite(T)
        Tm = np.where(m >= 0, T[mm], np.nan)
        d, size = pc.window_d(z, geo)
        ok = (m >= 0) & usable[mm] & pk[:, None]
        me = np.zeros(m.shape)                                   # raw member energy of each window module (0 when not a member)
        for c in np.flatnonzero(pk):
            hit = np.nonzero(m[c][:, None] == z["mi"][c][None, :])
            me[c, hit[0]] = z["me"][c][hit[1]]
        x = me / np.where(np.isfinite(Tm), Tm, 1.0)
        ib, it = np.digitize(d, DB) - 1, np.digitize(Tm, T_EDGES)
        nd, nt = DB.size - 1, T_EDGES.size + 1
        C, PA, NN = np.zeros((nd, nt)), np.full((nd, nt), np.nan), np.zeros((nd, nt), np.int64)
        pres = z["win_p"]
        for b in range(nd):
            for t in range(nt):
                k = ok & (ib == b) & (it == t)
                N = int(k.sum())
                NN[b, t] = N
                if N < MIN_N:
                    continue
                pa = float((~pres[k]).mean())
                PA[b, t] = pa
                if pa <= 0:
                    continue
                xp = x[k & pres & (me > 0)]
                lowm = pa + (xp < X_STAR).sum() / N
                rho = edge_density(xp, N)
                xbar = X_STAR * trunc_mean(rho * X_STAR / lowm)
                L = lowm * xbar - xp[xp < X_STAR].sum() / N
                C[b, t] = min(max(L / pa, 0.0), X_STAR) if PA_MIN <= pa <= PA_MAX else 0.0
        # per cluster added energy and removable variance (peak clusters)
        add = np.where(ok & ~pres & (ib >= 0) & (ib < nd), C[np.clip(ib, 0, nd - 1), np.clip(it, 0, nt - 1)] * np.nan_to_num(Tm), 0.0)
        pa_e = np.where(ok & (ib >= 0) & (ib < nd), np.nan_to_num(PA[np.clip(ib, 0, nd - 1), np.clip(it, 0, nt - 1)]), 0.0)
        c_e = np.where(ok & (ib >= 0) & (ib < nd), C[np.clip(ib, 0, nd - 1), np.clip(it, 0, nt - 1)] * np.nan_to_num(Tm), 0.0)
        var = (pa_e * (1 - pa_e) * c_e ** 2).sum(1)[pk]
        E = float(np.median(z["er"][pk]))
        sig = s * E
        print(f"== {rt}{' ' + a.half if a.half else ''}: {pk.sum()} peak clusters, E {E:.0f} MeV, sigma_E ~ {sig:.1f} MeV; added per cluster "
              f"mean {add.sum(1)[pk].mean():.1f} MeV (rms {add.sum(1)[pk].std():.1f}); removable variance <sum P(1-P)c^2> {var.mean():.1f} MeV^2 "
              f"-> sigma {sig:.1f} -> {math.sqrt(max(sig ** 2 - var.mean(), 0)):.1f} MeV (upper bound, independent modules)", flush=True)
        for b in range(nd):
            if np.isfinite(PA[b]).any() and np.nanmax(PA[b]) > 0.005 and np.nanmin(PA[b]) < 0.995:
                print(f"   d {DB[b]:.3f}  N {NN[b].tolist()}  P_a {np.round(PA[b], 3).tolist()}  c(x T) {np.round(C[b], 3).tolist()}")
        np.savez(Path(a.out) / f"{rt}.npz", c=C, p_absent=PA, n=NN, d_edges=DB, t_edges=T_EDGES, T=np.nan_to_num(T, nan=0.0), usable=usable,
                 half=2, size=size, x_star=X_STAR, run_type=rt, source=str(Path(a.root).resolve()), split_half=a.half or "all")


if __name__ == "__main__":
    main()
