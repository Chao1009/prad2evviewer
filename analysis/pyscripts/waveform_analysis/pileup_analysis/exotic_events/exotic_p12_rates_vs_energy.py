#!/usr/bin/env python3
"""exotic_p12_rates_vs_energy.py -- P12 of the exotic-waveform study: exotic-category rates versus energy slice, per channel and per cluster, and the impact per slice.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/exotic_events/exotic_p12_rates_vs_energy.py [--jobs 4] [--out <dir>] [--cache <npz>] [--plot-only]

Slices (MeV, as given; empty slices are dropped per run type): 10-50, 50-100, 100-300, 300-500, 500-1000, 1000-1500, 1500-2500, 2500-3500, 3500-4000, >4000.
Categories: the exclusive f4 > f5 > f1 > f2 > f3 of exotic_p4_module_rates.flags_of with the uniform (P2) thresholds, via exotic_p7_cluster_level.row_categories; run 24937 excluded; PbWO4 only.

1. Per channel.  Population: kept PbWO4 channels (A_main >= 10 ADC).  Energy = A_main x module gain (MeV; pileup_rates.module_gain, as P5), the default-fit main pulse (largest A).
   Rate = channels of the category / kept channels in the slice, Wilson 68 % interval.  Data: H2 run types from the sample npz (output/pileup_rates/exotic_sample; the H2 full-scan
   pickles carry no amplitude), Ta X17 from the 19 full-scan npz (output/pileup_rates/exotic_full).  Channels below 10 MeV (high-gain modules at the 10 ADC cut) are not in a slice.
2. Per cluster.  P7 cache (output/pileup_rates/exotic_p7_cache.npz; decoder exotic_p7_cluster_level): all PbWO4-seeded clusters, E_cluster = E_g (MeV, gain-corrected member sum).  Fraction with a
   flagged seed (fake-cluster candidates) and with >= 1 flagged member (seed or neighbour), per category and any; the cluster is counted in every category it holds, so the category columns can exceed 'any'.
3. Impact per slice.  Policy (user, 2026-11-18): flagged channels are dropped from clustering and calibration, no correction.  The impact is the energy the drop removes:
   -<fE>, minus the mean over all clusters of the slice of the flagged-energy fraction (sum of flagged member energies / E_g; P7 f_all), per category and any.  No P8 join.

Outputs (--out): exotic_P12_rates_vs_energy_20261118.pdf (per-channel rates + table, per-cluster fractions + table, impact + table) and .tsv (long: table = chan | cluster | impact).
The per-channel records are an intermediate cache (--cache, default output/pileup_rates/exotic_p12_chan_cache.npz; delete to rebuild).
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
import sys
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import exotic_p1_quantities as p1                       # noqa: E402
import exotic_p4_module_rates as p4                     # noqa: E402
import exotic_p7_cluster_level as p7                    # noqa: E402
from exotic_p1_quantities import TYPES, RT_LABEL, write_tsv, EXCLUDE_RUNS     # noqa: E402
from exotic_p3_sample import text_page                  # noqa: E402

STEM = "exotic_P12_rates_vs_energy_20261118"
EDGES = (10.0, 50.0, 100.0, 300.0, 500.0, 1000.0, 1500.0, 2500.0, 3500.0, 4000.0, np.inf)
CATS = p7.CATS                                          # 4, 5, 1, 2, 3
CNAME = p7.CNAME
CCOL = p7.CCOL
SCAN_ROOT = p1.PYDIR / "output" / "pileup_rates"
SOURCE = {"0p73GeV_H2": "H2 sample", "2p24GeV_H2": "H2 sample", "3p49GeV_H2": "H2 sample", "2p24GeV_Ta_X17": "Ta X17 full scan (19 npz)"}


def slab(i):
    lo, hi = EDGES[i], EDGES[i + 1]
    return f"{lo:g}-{hi:g}" if np.isfinite(hi) else f">{lo:g}"


def sl_index(e):
    return np.digitize(e, EDGES) - 1                    # -1: below the first edge


# ------------------------------------------------------------------------------------------------------------ per channel
def chan_file(fn):
    """(run type, E MeV, exclusive category) of the kept PbWO4 channels of one npz."""
    import pileup_rates as pr
    z = np.load(fn)
    rt, run = str(z["run_type"]), int(z["run"])
    midx = z["ch_midx"].astype(np.int64)
    sel = z["mod_pwo4"][midx]
    k = z["ch_k"].astype(np.int64)
    a, t = z["ch_a"].astype(np.float64), z["ch_t"]
    main = np.argmax(np.where(np.isfinite(t), a, -1.0), axis=1)
    amain = a[np.arange(k.size), main]
    t0 = t[np.arange(k.size), main].astype(np.float64)
    with np.errstate(invalid="ignore"):
        trunc = (np.ceil(t0 / p4.DT) - p4.LEAD + p4.LW > p4.NS) & np.isfinite(t0)
    c = dict(db=z["ch_ex_db"].astype(np.float64), tail=z["ch_ex_tail"].astype(np.float64), w50=z["ch_w50"].astype(np.float64), k=k,
             nm=z["ch_ex_nm"].astype(np.int64), c2s=z["ch_ex_chi2s"].astype(np.float64), c2k=z["ch_ex_chi2k"].astype(np.float64),
             nclip=z["ch_nclip"].astype(np.int64), trunc=trunc, s4=z["ch_ex_f4"].astype(bool))
    cat = p4.flags_of(c, p4.UNI)["cat"]
    gain = pr.module_gain(pr.run_config(run), rt, [str(n) for n in z["mod_name"]])
    e = amain * gain[midx]
    ok = sel & np.isfinite(e)
    return rt, run, e[ok].astype(np.float32), cat[ok].astype(np.int8), int(sel.sum() - ok.sum())


def chan_collect(cache):
    if cache and Path(cache).exists():
        c = np.load(cache, allow_pickle=True)
        print(f"[cache] {cache}", flush=True)
        return {rt: (c[f"e_{rt}"], c[f"c_{rt}"]) for rt in TYPES}, dict(c["meta"].item())
    files = [(f, "sample") for f in sorted((SCAN_ROOT / "exotic_sample").glob("run_*.npz"))]
    files += [(f, "full") for f in sorted((SCAN_ROOT / "exotic_full").glob("run_*_s*.npz")) if ".part" not in f.name]
    acc = {rt: ([], []) for rt in TYPES}
    meta = dict(nfile={rt: 0 for rt in TYPES}, runs={rt: set() for rt in TYPES}, nan={rt: 0 for rt in TYPES})
    for f, kind in files:
        with np.load(f) as z:
            run, rt = int(z["run"]), str(z["run_type"])
        if run in EXCLUDE_RUNS or (kind == "sample") == ("Ta" in rt):      # H2 from the sample, Ta X17 from the full scan only
            continue
        rt, run, e, cat, nnan = chan_file(f)
        acc[rt][0].append(e)
        acc[rt][1].append(cat)
        meta["nfile"][rt] += 1
        meta["runs"][rt].add(run)
        meta["nan"][rt] += nnan
        print(f"  {f.name}: {rt} {e.size} channels", flush=True)
    meta['below'] = {rt: int(sum((x < EDGES[0]).sum() for x in acc[rt][0])) for rt in TYPES}
    out = {rt: (np.concatenate(acc[rt][0]), np.concatenate(acc[rt][1])) for rt in TYPES}
    if cache:
        np.savez_compressed(cache, meta=np.array(meta, dtype=object), **{f"e_{rt}": out[rt][0] for rt in TYPES}, **{f"c_{rt}": out[rt][1] for rt in TYPES})
    return out, meta


def chan_rows(CH, meta):
    rows = []
    for rt in TYPES:
        e, cat = CH[rt]
        si = sl_index(e)
        for i in range(len(EDGES) - 1):
            m = si == i
            n = int(m.sum())
            if not n:
                continue
            for c in (0,) + CATS:
                k = int((cat[m] > 0).sum()) if c == 0 else int((cat[m] == c).sum())
                lo, hi = p4.wilson(k, n)
                rows.append(dict(table="chan", run_type=rt, source=SOURCE[rt], slice=slab(i), e_lo=EDGES[i], e_hi=EDGES[i + 1], category=CNAME[c], n=n, k=k,
                                 rate=k / n, rate_lo=float(lo), rate_hi=float(hi)))
    return rows


def fill_empty_chan(rows):
    """Channel rows (N = 0, k = 0, rate and interval NaN) for every (run type, slice) that has clusters but no channels, so the channel and cluster tables list identical slices."""
    have = {(r["run_type"], r["slice"]) for r in rows if r["table"] == "chan"}
    add, seen = [], set()
    for r in rows:
        key = (r["run_type"], r["slice"])
        if r["table"] == "cluster" and key not in have and key not in seen:
            seen.add(key)
            for c in (0,) + CATS:
                add.append(dict(table="chan", run_type=r["run_type"], source=r["source"], slice=r["slice"], e_lo=r["e_lo"], e_hi=r["e_hi"], category=CNAME[c], n=0, k=0, rate=np.nan, rate_lo=np.nan, rate_hi=np.nan))
    ch = sorted([r for r in rows if r["table"] == "chan"] + add, key=lambda r: (TYPES.index(r["run_type"]), r["e_lo"]))      # stable: category order kept
    return ch + [r for r in rows if r["table"] != "chan"]


# ------------------------------------------------------------------------------------------------------------ per cluster
def cluster_load():
    pm, D = p7.collect([], 1, str(SCAN_ROOT / "exotic_p7_cache.npz"))
    with np.load(next((SCAN_ROOT / "exotic_sample").glob("run_*.npz"))) as z:
        names = [str(x) for x in z["mod_name"]]
    X = p7.derive(pm, D, names)
    ds = {}
    for rt in TYPES:
        src = "full" if "Ta" in rt else "sample"
        ds[rt] = (X["rt"] == rt) & (X["src"] == src)
    return pm, D, X, ds


def cluster_rows(pm, D, X, ds):
    rows = []
    scat, nmask = D["scat"].astype(np.int64), D["nmask"].astype(np.int64)
    eg = D["eg"].astype(np.float64)
    si = sl_index(eg)
    fE = {0: (X["f_all"], scat > 0, (scat > 0) | (nmask > 0))}
    for c in CATS:
        fE[c] = (np.where(scat == c, X["f_seed"], 0.0) + X["f_nb"][:, c - 1], scat == c, (scat == c) | ((nmask >> (c - 1) & 1) > 0))
    for rt in TYPES:
        for i in range(len(EDGES) - 1):
            m = ds[rt] & (si == i)
            n = int(m.sum())
            if not n:
                continue
            sume = float(eg[m].sum())
            base = dict(run_type=rt, source=SOURCE[rt], slice=slab(i), e_lo=EDGES[i], e_hi=EDGES[i + 1], n=n)
            for c in (0,) + CATS:
                f, seed, anym = fE[c]
                for role, flag in (("seed", seed), ("any", anym)):
                    k = int((m & flag).sum())
                    lo, hi = p4.wilson(k, n)
                    rows.append(dict(base, table="cluster", category=CNAME[c], role=role, k=k, rate=k / n, rate_lo=float(lo), rate_hi=float(hi)))
                fm = f[m]
                rows.append(dict(base, table="impact", category=CNAME[c], role="fE", k=int((m & anym).sum()), fE_mean=-float(fm.mean()), fE_sem=float(fm.std(ddof=1) / np.sqrt(n)) if n > 1 else np.nan))
    return rows


# ------------------------------------------------------------------------------------------------------------ output
def get(rows, **kw):
    return next((r for r in rows if all(r.get(k) == v for k, v in kw.items())), None)


def pct(x, d=4):
    return f"{100 * x:.{d}f}" if x == x else "nan"


def slices_of(rows, rt, table):
    seen = []
    for r in rows:
        if r["table"] == table and r["run_type"] == rt and r["slice"] not in seen:
            seen.append(r["slice"])
    return seen


SMALL_N = 100                                           # slices with fewer entries are marked '*' (interval wide, rate poorly known)
FLAG_ORDER = (1, 2, 3, 4, 5)                            # table columns f1..f5, then any


def fmt_ci(r):
    """'rate % [lo, hi]' (Wilson 68 %) with ~3 significant figures; zero-count cells keep their upper limit."""
    if r["n"] == 0:
        return "-"
    def f(x):
        v = 100 * x
        if v <= 0:
            return "0"
        return f"{v:.{min(6, max(0, 2 - int(np.floor(np.log10(v)))))}f}"
    return f"{f(r['rate'])} [{f(r['rate_lo'])}, {f(r['rate_hi'])}]"


def flag_lines(rows, table, role, rt):
    """Lines of the per-flag x per-slice table of one run type: rows = slices, columns = N, f1..f5, any."""
    cols = [CNAME[c] for c in FLAG_ORDER] + ["any"]
    w = 32
    L = [f"{'slice (MeV)':<12}{'N':>11}  " + "".join(f"{c.split()[0] + ' % [lo, hi]':>{w}}" for c in cols)]
    for sl in slices_of(rows, rt, table):
        rr = [get(rows, table=table, run_type=rt, slice=sl, category=c, **({"role": role} if role else {})) for c in cols]
        n = rr[0]["n"]
        L.append(f"{sl:<12}{n:>10d}{'*' if n < SMALL_N else ' '} " + "".join(f"{fmt_ci(r):>{w}}" for r in rr))
    return L


def table_pages(pdf, plt, rows):
    """One page per run type at channel level, then one per run type at cluster level (user, 2026-11-18)."""
    for table, role, what, pop in (("chan", None, "channel", "kept PbWO4 channels (A_main >= 10 ADC) of the slice, energy = A_main x gain"),
                                   ("cluster", "any", "cluster", "PbWO4-seeded clusters of the slice (P7 cache), E_cluster = E_g; a cluster counts if it holds >= 1 flagged member (seed or neighbour)")):
        for rt in TYPES:
            L = [f"Rate (%) of {pop}.", f"{RT_LABEL[rt]}  [{SOURCE[rt]}];  run 24937 excluded.  Cell = rate % [lo, hi], Wilson 68 % interval; exclusive categories f4 > f5 > f1 > f2 > f3 "
                 + ("(a cluster is counted in every category it holds, so the sum of f1-f5 can exceed 'any')." if table == "cluster" else "(each channel in one category; any = sum)."), f"'*' after N: fewer than {SMALL_N} entries in the slice (small-N, rate poorly determined).  Zero-count cells show the upper limit.", ""]
            L += flag_lines(rows, table, role, rt)
            text_page(pdf, plt, f"P12 {what} level, {RT_LABEL[rt]}: rate per flag per energy slice", L, size=7.5)


def impact_table(rows):
    L = ["Impact per E_cluster slice of dropping the flagged channels (no correction): energy removed = -<fE>, the mean flagged-energy fraction of all clusters (%, +- s.e.m.), per category and any.",
         "A cluster holding two categories adds to both.  fE = sum of flagged member energies / E_g.", ""]
    for rt in TYPES:
        L += [f"{RT_LABEL[rt]}  [{SOURCE[rt]}]", f"  {'slice (MeV)':<12}{'N clus':>10}{'-<fE> any %':>13}{'(sem)':>10}{'dominant':>10}" + "".join(f"{CNAME[c].split()[0] + ' (+-sem)':>21}" for c in CATS)]
        for s in slices_of(rows, rt, "impact"):
            a = get(rows, table="impact", run_type=rt, slice=s, category="any", role="fE")
            cs = {c: get(rows, table="impact", run_type=rt, slice=s, category=CNAME[c], role="fE") for c in CATS}
            dom = min(CATS, key=lambda c: cs[c]["fE_mean"])
            L.append(f"  {s:<12}{a['n']:>10d}{pct(a['fE_mean'], 5):>13}{pct(a['fE_sem'], 5):>10}{(CNAME[dom].split()[0] if a['fE_mean'] < 0 else '-'):>10}"
                     + "".join(f"{pct(cs[c]['fE_mean'], 5) + '+-' + pct(cs[c]['fE_sem'], 5):>21}" for c in CATS))
        L.append("")
    return L


def draw_rates(axs, rows, table, role, ycol_any=True):
    """2x2: category rate vs slice index (log y); table 'chan' or 'cluster'."""
    for ax, rt in zip(np.ravel(axs), TYPES):
        sl = slices_of(rows, rt, table)
        x = np.arange(len(sl))
        for c in CATS:
            rr = [get(rows, table=table, run_type=rt, slice=s, category=CNAME[c], **({"role": role} if role else {})) for s in sl]
            y = np.array([100 * r["rate"] for r in rr])
            ok = np.array([r["k"] > 0 for r in rr])
            lo = np.array([100 * (r["rate"] - r["rate_lo"]) for r in rr])
            hi = np.array([100 * (r["rate_hi"] - r["rate"]) for r in rr])
            ax.errorbar(x[ok], y[ok], yerr=[lo[ok], hi[ok]], fmt="-o", ms=4, lw=0.9, capsize=1.5, color=CCOL[c], label=CNAME[c])
        ra = [get(rows, table=table, run_type=rt, slice=s, category="any", **({"role": role} if role else {})) for s in sl]
        ya = np.array([100 * r["rate"] for r in ra])
        oka = np.array([r["k"] > 0 for r in ra])
        ax.errorbar(x[oka], ya[oka], yerr=[[max(0.0, 100 * (r["rate"] - r["rate_lo"])) for r in ra if r["k"]], [max(0.0, 100 * (r["rate_hi"] - r["rate"])) for r in ra if r["k"]]], fmt="k--", lw=1.0, capsize=1.5, label="any")
        ax.set_yscale("log")
        ax.set_xticks(x)
        ax.set_xticklabels(sl, rotation=45, ha="right", fontsize=7)
        ax.set_title(f"{RT_LABEL[rt]}  [{SOURCE[rt]}]", fontsize=9)
        ax.grid(alpha=0.25, lw=0.4)


def report(rows, meta, out, extra):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    out = Path(out)
    out.mkdir(parents=True, exist_ok=True)
    with PdfPages(out / (STEM + ".pdf")) as pdf:
        head = ["P12: exotic categories versus energy slice.  Slices (MeV): " + ", ".join(slab(i) for i in range(len(EDGES) - 1)),
                "Exclusive categories f4 > f5 > f1 > f2 > f3 (uniform thresholds); PbWO4 only; run 24937 excluded.",
                "Per channel: energy = default-fit A_main x module gain; kept = A_main >= 10 ADC.  H2 from the sample (H2 full-scan pickles carry no amplitude); Ta X17 from the 19 full-scan npz.",
                "Per cluster: P7 cache, E_cluster = E_g; impact = energy removed by dropping flagged channels, -<fE>, from the P7 cache.", ""]
        for rt in TYPES:
            head.append(f"  {RT_LABEL[rt]:<18}{SOURCE[rt]:<28}files {meta['nfile'][rt]:>3}, runs {len(meta['runs'][rt]):>3}, kept channels {int(sum(r['n'] for r in rows if r['table'] == 'chan' and r['run_type'] == rt and r['category'] == 'any')):>12d} (channels below 10 MeV, not in a slice: {meta['below'][rt]}; non-finite gain: {meta['nan'][rt]})")
        head += ["", "Slices above E_beam (e.g. > 4000 MeV at 2.24 GeV) are populated mainly by f4/f5 clipped records whose fitted A_main is an extrapolation: their rate near 100 % is by construction, not a physical rise.", extra]
        text_page(pdf, plt, "P12 rates vs energy slice: data and definitions", head, size=9)
        fig, axs = plt.subplots(2, 2, figsize=(16, 11), sharey=True, layout="constrained")
        draw_rates(axs, rows, "chan", None)
        axs[0, 0].legend(fontsize=8)
        for ax in axs[1]:
            ax.set_xlabel("channel energy A_main x gain (MeV)")
        for ax in axs[:, 0]:
            ax.set_ylabel("rate per kept channel (%)")
        fig.suptitle("Per-channel exclusive category rate vs energy slice (Wilson 68 %; markers only where the category has entries); slices above E_beam: f4/f5 extrapolated amplitudes", fontsize=11)
        pdf.savefig(fig)
        plt.close(fig)
        fig, axs = plt.subplots(2, 2, figsize=(16, 11), sharey=True, layout="constrained")
        draw_rates(axs, rows, "cluster", "any")
        axs[0, 0].legend(fontsize=8, title="clusters with >= 1 flagged member")
        for ax in axs[1]:
            ax.set_xlabel("cluster energy E_cluster (MeV)")
        for ax in axs[:, 0]:
            ax.set_ylabel("fraction of clusters (%)")
        fig.suptitle("Per-cluster fraction with >= 1 flagged member (seed or neighbour) by category; dashed = any", fontsize=12)
        pdf.savefig(fig)
        plt.close(fig)
        fig, axs = plt.subplots(2, 2, figsize=(16, 11), layout="constrained")
        for ax, rt in zip(axs.ravel(), TYPES):
            sl = slices_of(rows, rt, "impact")
            x = np.arange(len(sl))
            for c in (0,) + CATS:
                rr = [get(rows, table="impact", run_type=rt, slice=s, category=CNAME[c], role="fE") for s in sl]
                y = -100 * np.array([r["fE_mean"] for r in rr])
                ye = 100 * np.array([r["fE_sem"] for r in rr])
                ok = y > 0
                ax.errorbar(x[ok], y[ok], yerr=np.minimum(ye[ok], 0.999 * y[ok]), fmt="--" if c == 0 else "-o", ms=4, lw=0.9, capsize=1.5, color="k" if c == 0 else CCOL[c], label=CNAME[c])
            ax.set_yscale("log")
            ax.set_xticks(x)
            ax.set_xticklabels(sl, rotation=45, ha="right", fontsize=7)
            ax.set_title(f"{RT_LABEL[rt]}: energy removed by the drop, -<fE> (%, plotted as magnitude)", fontsize=9)
            ax.set_xlabel("cluster energy E_cluster (MeV)")
            ax.grid(alpha=0.25, lw=0.4)
        axs[0, 0].legend(fontsize=7, loc="upper left")
        fig.suptitle("Impact per E_cluster slice: energy removed from the cluster by dropping the flagged channels, -<fE> (magnitude plotted; error bars = s.e.m.)", fontsize=12)
        pdf.savefig(fig)
        plt.close(fig)
        text_page(pdf, plt, "Impact per slice (table)", impact_table(rows), size=6.0)
        table_pages(pdf, plt, rows)
    for r in rows:
        r["small_n"] = int(r["n"] < SMALL_N)
    cols = ["table", "run_type", "source", "slice", "e_lo", "e_hi", "category", "role", "n", "k", "rate", "rate_lo", "rate_hi", "fE_mean", "fE_sem", "small_n"]
    write_tsv(out / (STEM + ".tsv"), rows, cols)
    print("wrote", out / (STEM + ".pdf"), out / (STEM + ".tsv"))


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--out", default=str(HERE))
    ap.add_argument("--cache", default=str(SCAN_ROOT / "exotic_p12_chan_cache.npz"))
    ap.add_argument("--plot-only", action="store_true", help="accepted for symmetry; the cache is reused whenever it exists")
    a = ap.parse_args()
    CH, meta = chan_collect(a.cache)
    rows = chan_rows(CH, meta)
    pm, D, X, ds = cluster_load()
    rows += cluster_rows(pm, D, X, ds)
    rows = fill_empty_chan(rows)
    extra = "Impact: -<fE> = energy removed from each cluster by dropping its flagged channels (user policy), P7 cache; the P8 credible-correction join is not used."
    report(rows, meta, a.out, extra)


if __name__ == "__main__":
    main()
