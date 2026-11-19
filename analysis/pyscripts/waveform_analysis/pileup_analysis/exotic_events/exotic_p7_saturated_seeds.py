#!/usr/bin/env python3
"""exotic_p7_saturated_seeds.py -- P7 follow-up: the saturated (f4/f5) seeds of the elastic-window clusters of the 2.24 and 3.49 GeV H2 sample.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p7_saturated_seeds.py

Reads output/pileup_rates/exotic_p7_cache.npz (read only; elastic-window selection exactly as exotic_p7_cluster_level.populations) and joins each cluster's
seed row to the sample npz (exotic_sample/run_*.npz).  Writes exotic_P7_saturated_seeds_20261118.tsv (one row per cluster) next to this script.
Population (a) = n_members <= 2, (b) = n_members >= 3.  A_main = largest fitted amplitude (ADC); full scale = 4095 - ped.
Appended columns (repeated on every row of the run type): n_sat / n_elastic = saturated seeds / elastic-window clusters, sat_frac with its Wilson 68 % interval
(sat_lo68, sat_hi68; fractions, not %); pop_n = rows of this row's population (a / b) in the run type, pop_frac = pop_n / n_sat with Wilson 68 % (pop_lo68, pop_hi68).
"""
import sys
from pathlib import Path
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
import numpy as np
import exotic_p4_module_rates as p4
import exotic_p7_cluster_level as p7
from exotic_p1_quantities import write_tsv

OUT = HERE / "exotic_P7_saturated_seeds_20261118.tsv"
COLS = ("dataset", "run", "evn", "module", "scat", "n_members", "seed_frac", "k", "a_main", "a_over_fs", "chi2_dof", "nclip", "w50", "pop",
        "n_sat", "n_elastic", "sat_frac", "sat_lo68", "sat_hi68", "pop_n", "pop_frac", "pop_lo68", "pop_hi68")


def main():
    metas, D = p7.collect([], 1, str(p7.SCAN_ROOT / "exotic_p7_cache.npz"))
    with np.load(next((p7.SCAN_ROOT / "exotic_sample").glob("run_*.npz"))) as z:
        names = [str(x) for x in z["mod_name"]]
    X = p7.derive(metas, D, names)
    rows = []
    for ds in p7.datasets(metas):
        if ds[0] not in ("2p24GeV_H2", "3p49GeV_H2"):
            continue
        P, _, _ = p7.populations(metas, D, X, ds, (None, {}, ""))
        sel = np.flatnonzero(P["elastic"] & np.isin(D["scat"], (4, 5)))
        cache = {}
        for j in sel:
            m = metas[int(D["job"][j])]
            fn = m["fn"]
            if fn not in cache:
                z = np.load(fn)
                cache[fn] = {q: z[q] for q in ("ev_number", "ch_midx", "ch_ev", "ch_k", "ch_a", "ch_t", "ch_w50", "ch_nclip", "ch_chi2", "ch_ped")}
            z = cache[fn]
            ei = np.flatnonzero(z["ev_number"] == D["evn"][j])
            assert ei.size == 1
            r = np.flatnonzero((z["ch_ev"] == ei[0]) & (z["ch_midx"] == D["seed"][j]) & (z["ch_nclip"] > 0))
            assert r.size == 1, (ds[0], m["run"], D["evn"][j], r.size)
            i = r[0]
            a = float(np.nanmax(np.where(np.isfinite(z["ch_t"][i]), z["ch_a"][i], -1.0)))
            ped = float(z["ch_ped"][i])
            rows.append(dict(dataset=ds[0], run=m["run"], evn=int(D["evn"][j]), module=names[int(D["seed"][j])], scat=int(D["scat"][j]), n_members=int(D["cnt"][j]),
                             seed_frac=float(X["f_seed"][j]), k=int(z["ch_k"][i]), a_main=a, a_over_fs=a / (4095.0 - ped), chi2_dof=float(z["ch_chi2"][i]),
                             nclip=int(z["ch_nclip"][i]), w50=float(z["ch_w50"][i]), pop="a" if D["cnt"][j] <= 2 else "b"))
        n_el, n_sat = int(P["elastic"].sum()), int(sel.size)
        sl, sh = p4.wilson(n_sat, n_el)
        pn = {q: sum(r["pop"] == q for r in rows if r["dataset"] == ds[0]) for q in "ab"}
        pw = {q: p4.wilson(pn[q], n_sat) for q in "ab"}
        for r in rows:
            if r["dataset"] == ds[0]:
                r.update(n_sat=n_sat, n_elastic=n_el, sat_frac=n_sat / n_el, sat_lo68=float(sl), sat_hi68=float(sh), pop_n=pn[r["pop"]], pop_frac=pn[r["pop"]] / n_sat,
                         pop_lo68=float(pw[r["pop"]][0]), pop_hi68=float(pw[r["pop"]][1]))
        print(f"{ds[0]} elastic {n_el} sat {n_sat} = {100 * n_sat / n_el:.4f} % [{100 * sl:.4f}, {100 * sh:.4f}] (Wilson 68 %); "
              + ", ".join(f"pop {q}: {pn[q]} = {100 * pn[q] / n_sat:.1f} % [{100 * pw[q][0]:.1f}, {100 * pw[q][1]:.1f}]" for q in "ab"))
    write_tsv(OUT, rows, COLS)
    print("wrote", OUT, len(rows))


if __name__ == "__main__":
    main()
