"""P3 (A): width and tail excess per (dt, A2/A1) cell of the secondary pulse, against pool clusters without one.

Usage: python3 p3_cells.py <P3 variants dir> <out tsv> [run types ...]     (reads <dir>/<rt>/secondaries.npz, SECONDARY_DUMP=True)
Per cluster, E/E' is the scan-level ratio (no gains, no position correction).  Reference: n_extra == 0.  mu, sigma: iterated
Gaussian window (mu +- 2 sigma) of the reference peak.  Per cell: share, peak sigma (cell, window of its own) minus reference,
tail_lo3 / tail_hi3 = % of in-range (0.6-1.4) clusters below mu - 3 sigma_ref / above mu + 3 sigma_ref, minus the reference; d_med_pct: shift of the median within mu +- 0.1 (pp); e_extra: median secondary energy.
SE: binomial for tails, 30 bootstraps for the width.  flag = |excess| >= 2 SE with n >= 500.
"""
import sys
import numpy as np

DT = (-np.inf, 0, 12, 16, 20, 30, np.inf)
AR = (0, 0.04, 0.07, 0.12, 0.25, 1.0, np.inf)


def peak(x, mu=1.0, sg=0.04):
    for _ in range(8):
        w = x[np.abs(x - mu) < 2 * sg]
        if w.size < 50:
            return np.nan, np.nan
        mu, sg = w.mean(), w.std() * 1.15      # 2-sigma truncation shrinks std by ~0.88
    return mu, sg


def main():
    d, out, rts = sys.argv[1], sys.argv[2], sys.argv[3:] or ["2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]
    rng = np.random.default_rng(1)
    beam = {"2p24GeV_H2": 2239.5, "3p49GeV_H2": 3487.0, "2p24GeV_Ta_X17": 2239.5}
    rows = ["run_type\tdt_lo\tdt_hi\tar_lo\tar_hi\tn\tshare\td_sigma_pct\td_sigma_se\ttail_lo3\td_tail_lo3\td_tail_lo3_se\ttail_hi3\td_tail_hi3\t"
            "d_tail_hi3_se\tflag\td_med_pct\te_extra_med_mev\te_extra_pct_beam"]
    for rt in rts:
        z = np.load(f"{d}/{rt}/secondaries.npz")
        r = z["e_over_ep"]
        ok = (r > 0.6) & (r < 1.4)
        ref = ok & (z["n_extra"] == 0)
        mu, sg = peak(r[ref])
        sg0 = sg
        tl = lambda m: (np.mean(r[m] < mu - 3 * sg0) * 100, np.mean(r[m] > mu + 3 * sg0) * 100)
        rl, rh = tl(ref)
        print(f"{rt}: ref n={ref.sum()} mu={mu:.4f} sigma={100 * sg:.3f}% tail_lo3={rl:.2f} tail_hi3={rh:.2f}; "
              f"with secondary {np.mean(z['n_extra'][ok] > 0) * 100:.1f}% of in-range")
        nref = ref.sum()
        for i in range(len(DT) - 1):
            for j in range(len(AR) - 1):
                m = ok & (z["n_extra"] > 0) & (z["dt"] >= DT[i]) & (z["dt"] < DT[i + 1]) & (z["ar"] >= AR[j]) & (z["ar"] < AR[j + 1])
                n = int(m.sum())
                if n < 100:
                    continue
                cl, ch = tl(m)
                sc = []
                x = r[m]
                if n >= 500:
                    for _ in range(30):
                        sc.append(peak(rng.choice(x, n))[1])
                dsg = (peak(x)[1] - sg0) * 100
                se = lambda p, nn, p0: 100 * np.sqrt(p * (1 - p) / nn + p0 * (1 - p0) / nref)
                sl, sh = se(cl / 100, n, rl / 100), se(ch / 100, n, rh / 100)
                ssg = np.nanstd(sc) * 100 if sc else np.nan
                flag = ""
                if n >= 500:
                    flag = ",".join(k for k, v, s in (("sigma", dsg, ssg), ("lo", cl - rl, sl), ("hi", ch - rh, sh)) if abs(v) >= 2 * s)
                rows.append("\t".join(str(v) for v in (rt, DT[i], DT[i + 1], AR[j], AR[j + 1], n, f"{100 * n / ok.sum():.2f}", f"{dsg:.4f}",
                                                       f"{ssg:.4f}", f"{cl:.3f}", f"{cl - rl:.3f}", f"{sl:.3f}", f"{ch:.3f}", f"{ch - rh:.3f}",
                                                       f"{sh:.3f}", flag, f"{100 * (np.median(x[(x > mu - 0.1) & (x < mu + 0.1)]) - np.median(r[ref & (r > mu - 0.1) & (r < mu + 0.1)])):.3f}",
                                                       f"{np.median(z['e_extra'][m]):.1f}", f"{100 * np.median(z['e_extra'][m]) / beam[rt]:.2f}")))
    open(out, "w").write("\n".join(rows) + "\n")


if __name__ == "__main__":
    main()
