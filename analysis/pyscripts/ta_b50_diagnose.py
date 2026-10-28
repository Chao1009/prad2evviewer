#!/usr/bin/env python3
"""ta_b50_diagnose.py -- why the best-50 (B50) bank + gains (R4) is worse than R3 on the Ta (X17) runs.

Ta runs are always reported per run (25800, 26138, 26230), never pooled.  Subcommands (each writes
into output/ta_b50_diag_20261028/):

  modules    (i)  per-run, per-seed-module paired R4 - R3 comparison.  Events of R3 and R4 are
                  paired (the two scans keep slightly different event lists, so the per-event
                  table is aligned by seed module + cluster position + raw energy); per module
                  the post-calibration CB sigma/mu of both is fitted (benchmark fit_peak_cb, as
                  plot_module_resolution.py), delta = R4 - R3 with a paired bootstrap SE.
                  Also per-radius groups and the per-split delta from compare/r4_vs_r3_splits.tsv.
  templates  (ii) per-module R3 production vs B50 Ta template: max|df|, d rise, d t_pk, and the
                  Spearman correlation with the module delta of `modules`.
  gains      (iii) B50 / percrystal_spl gain ratio per module for the three Ta runs, outliers,
                  Spearman with the module delta; H2 ratio spread for reference.
  select     (iv) cut counts per Ta run from the shipped select logs, and the amplitude composition of the
                  build pool by radius (Ta vs 2.24 GeV H2).
  selvar     (iv) per-split comparison of relaxed-cut select variants (output/ta_b50_diag_20261028/
                  select_variants/<name>_<run>_<split>.npz, made with spline_template_extract.py select and
                  the B50 options minus one cut) against the shipped B50 select.
  runs       per-run (never pooled) paired table R3 / R4 / fix, legacy and gated CB, with the
                  |Delta| <= 1 SE promotion check (reads the compare_scan_dirs *_splits.tsv).

`modules --cand DIR --tag T` repeats (i) for another candidate scan (e.g. a fix scan) against R3.

Usage: python3 ta_b50_diagnose.py modules [--boot 100] [--min-n 150]
"""
from __future__ import annotations

import argparse
import contextlib
import io
import re
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
OUT = HERE / "output/ta_b50_diag_20261028"
R3 = HERE / "output/resolution_fit_20261026/r3_splits"
R4 = HERE / "output/resolution_fit_20261027_best50/r4_splits"
TA_RUNS = (25800, 26138, 26230)


def _pmr():
    import plot_module_resolution as P
    return P


# ---------------------------------------------------------------------------------------------
# event pairing
# ---------------------------------------------------------------------------------------------
def align(x: pd.DataFrame, y: pd.DataFrame, look: int = 6) -> np.ndarray:
    """Pair the rows of two per-event cluster tables of the same config (R3 and R4 keep slightly
    different event lists).  Match: same seed module, cluster x and y within 8 mm, raw energy
    within 20 %, sequential with a +-look skip window.  Returns an (n, 2) index array."""
    sx, sy = x.seed_module_id.values, y.seed_module_id.values
    ex, ey = x.cluster_E_raw_MeV.values, y.cluster_E_raw_MeV.values
    px, py = x.cluster_x_mm.values, y.cluster_x_mm.values
    qx, qy = x.cluster_y_mm.values, y.cluster_y_mm.values

    def ok(i, j):
        return (sx[i] == sy[j] and abs(px[i] - py[j]) < 8 and abs(qx[i] - qy[j]) < 8
                and abs(ex[i] - ey[j]) < 0.2 * ex[i])
    i = j = 0
    out = []
    while i < len(x) and j < len(y):
        if ok(i, j):
            out.append((i, j)); i += 1; j += 1; continue
        for d in range(1, look + 1):
            if i + d < len(x) and ok(i + d, j):
                i += d; break
            if j + d < len(y) and ok(i, j + d):
                j += d; break
        else:
            i += 1; j += 1
    return np.array(out)


def load_paired(run: int, path: str = "legacy", r3: Path = R3, r4: Path = R4, inner=200.0,
                rowcol=(15, 20, 15, 20)):
    """Paired post-calibration events of one run (all complete splits) in the benchmark's
    calibration pool (seed radius, min seed energy, outside the row/col rectangle, finite E_cal in
    both scans).  Returns (DataFrame[mod, e3, e4, raw3, raw4, split], hdr)."""
    P = _pmr()
    cfg = "legacy" if path == "legacy" else "gated"
    r0, r1, c0, c1 = rowcol
    rows, hdr = [], None
    for d in sorted((r3 / f"scan_{run:06d}").glob("split_*")):
        d4 = r4 / f"scan_{run:06d}" / d.name
        if not (all((d / f).exists() for f in P.NEEDED) and all((d4 / f).exists() for f in P.NEEDED)):
            continue
        h = P.parse_header(d / "result.tsv")
        hdr = hdr or h
        a = pd.read_csv(d / "result_clusters.tsv", sep="\t")
        b = pd.read_csv(d4 / "result_clusters.tsv", sep="\t")
        a, b = (t[t.config_name == cfg].reset_index(drop=True) for t in (a, b))
        m = align(a, b)
        xa, xb = a.iloc[m[:, 0]].reset_index(drop=True), b.iloc[m[:, 1]].reset_index(drop=True)
        inrect = (xa.seed_row + 1).between(r0, r1) & (xa.seed_col + 1).between(c0, c1)
        keep = ((xa.seed_x_mm ** 2 + xa.seed_y_mm ** 2 < inner ** 2) & (xa.seed_energy_MeV >= h["min_seed"])
                & ~inrect & np.isfinite(xa.cluster_E_cal_MeV) & np.isfinite(xb.cluster_E_cal_MeV))
        rows.append(pd.DataFrame({"mod": xa.seed_module_id[keep], "e3": xa.cluster_E_cal_MeV[keep],
                                  "e4": xb.cluster_E_cal_MeV[keep], "raw3": xa.cluster_E_raw_MeV[keep],
                                  "raw4": xb.cluster_E_raw_MeV[keep],
                                  "split": int(d.name.split("_")[1]), "n_aligned": len(m),
                                  "n3": len(a), "n4": len(b)}))
    return pd.concat(rows, ignore_index=True), hdr


def _iqr_sd(x) -> float:
    x = np.asarray(x)
    return float((np.quantile(x, .75) - np.quantile(x, .25)) / 1.349) if len(x) > 3 else float("nan")


def cb_sigma(vals, hdr) -> tuple[float, float]:
    """(mu, sigma/mu in %) of the benchmark CB fit; nan if the fit fails."""
    P = _pmr()
    try:
        f = P.fit_one(np.asarray(vals), hdr)
    except Exception:
        return float("nan"), float("nan")
    return f["cb_mu"], 100 * P.ratio(f["cb_sigma"], f["cb_mu"])


def paired_delta(e3, e4, hdr, boot, rng):
    """(s3, s4, delta, SE) in %, CB sigma/mu of both arrays; SE of delta by paired event bootstrap."""
    s3 = cb_sigma(e3, hdr)[1]
    s4 = cb_sigma(e4, hdr)[1]
    d = []
    n = len(e3)
    for _ in range(boot):
        k = rng.integers(0, n, n)
        d.append(cb_sigma(e4[k], hdr)[1] - cb_sigma(e3[k], hdr)[1])
    d = np.array(d)
    return s3, s4, s4 - s3, float(np.nanstd(d))


# ---------------------------------------------------------------------------------------------
# (i) modules
# ---------------------------------------------------------------------------------------------
def cmd_modules(args):
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(1)
    P = _pmr()
    geo = P.load_geometry(25800)      # indexed by module name (W494 ...)
    allrows, grp_rows = [], []
    for run in TA_RUNS:
        df, hdr = load_paired(run, args.path, r4=Path(args.cand))
        n_al = df.groupby("split").first()
        print(f"run {run}: {len(df)} paired pool events, {df.split.nunique()} splits; "
              f"aligned/kept per split {int(n_al.n_aligned.sum())}/{int(n_al.n3.sum())} (R3), {int(n_al.n4.sum())} (R4)")
        # robust per-event shift (does the bank change single energies?)
        rel = (df.e4 - df.e3) / df.e3
        print(f"   per-event (E4-E3)/E3: median {100*rel.median():+.3f} %  IQR/1.35 {100*(rel.quantile(.75)-rel.quantile(.25))/1.35:.3f} %"
              f"   raw: median {100*((df.raw4-df.raw3)/df.raw3).median():+.3f} %")
        tot = paired_delta(df.e3.values, df.e4.values, hdr, args.boot, rng)
        print(f"   all modules: CB s/mu R3 {tot[0]:.3f}  R4 {tot[1]:.3f}  delta {tot[2]:+.3f} +- {tot[3]:.3f} pp")
        grp_rows.append(dict(run=run, group="all", n=len(df), s3=tot[0], s4=tot[1], delta=tot[2], se=tot[3]))
        mods = df["mod"].unique()
        for m in mods:
            s = df[df["mod"] == m]
            r = float(geo.loc[m].r) if m in geo.index else np.nan
            row = dict(run=run, module=m, r=r, n=len(s), nsplit=s.split.nunique(), med_raw3=s.raw3.median(),
                       frac_shift=100 * ((s.e4 - s.e3) / s.e3).median(),
                       shift_sd=100 * _iqr_sd((s.e4 - s.e3) / s.e3),
                       raw_shift_sd=100 * _iqr_sd((s.raw4 - s.raw3) / s.raw3),
                       s3=np.nan, s4=np.nan, delta=np.nan, se=np.nan)
            if len(s) >= args.min_n:
                row["s3"], row["s4"], row["delta"], row["se"] = paired_delta(
                    s.e3.values, s.e4.values, hdr, args.boot, rng)
            allrows.append(row)
        # radius groups (the 8 modules with module medians are the r<80 ones)
        rr = df["mod"].map(lambda m: float(geo.loc[m].r) if m in geo.index else np.nan)
        for name, mask in (("r<80 mm", rr < 80), ("80-100 mm", (rr >= 80) & (rr < 100)),
                           ("100-200 mm", rr >= 100)):
            s = df[mask]
            if len(s) >= args.min_n:
                t = paired_delta(s.e3.values, s.e4.values, hdr, args.boot, rng)
                grp_rows.append(dict(run=run, group=name, n=len(s), s3=t[0], s4=t[1], delta=t[2], se=t[3]))
    mt = pd.DataFrame(allrows)
    gt = pd.DataFrame(grp_rows)
    tag = f"_{args.tag}" if args.tag else ""
    mt.to_csv(OUT / f"modules{tag}_{args.path}.tsv", sep="\t", index=False, float_format="%.4f")
    gt.to_csv(OUT / f"groups{tag}_{args.path}.tsv", sep="\t", index=False, float_format="%.4f")
    print("\nper radius group (CB s/mu %, paired delta pp):")
    print(gt.round(3).to_string(index=False))
    for run in TA_RUNS:
        t = mt[(mt.run == run) & mt.delta.notna()].sort_values("delta", ascending=False)
        print(f"\nrun {run}: {len(t)} fitted modules (n >= {args.min_n}); delta mean {t.delta.mean():+.3f} "
              f"median {t.delta.median():+.3f}; modules with delta > +2 SE: {(t.delta > 2 * t.se).sum()}, "
              f"< -2 SE: {(t.delta < -2 * t.se).sum()}; weighted mean "
              f"{np.average(t.delta, weights=1 / t.se ** 2):+.3f} +- {1 / np.sqrt((1 / t.se ** 2).sum()):.3f}")
        print(t[["module", "r", "n", "s3", "s4", "delta", "se"]].round(3).head(8).to_string(index=False))


# ---------------------------------------------------------------------------------------------
# (ii) templates
# ---------------------------------------------------------------------------------------------
TPL3 = HERE / "output/templates_spline/2p24GeV_Ta_X17/templates_2p24GeV_Ta_X17.npz"
TPL4 = HERE / "output/templates_spline_best50_20261026/2p24GeV_Ta_X17/templates_2p24GeV_Ta_X17.npz"


def template_table(p3=TPL3, p4=TPL4, only=None) -> pd.DataFrame:
    """Per-module R3 (production) vs B50 template differences (normalised shapes, onset frame):
    max|df| in -5..60 ns (% of peak), d rise 10-90 / w99 / FWHM / t_pk (ns), relative d iint_wa."""
    import template_doublet_audit as T
    z3, z4 = np.load(p3, allow_pickle=True), np.load(p4, allow_pickle=True)
    names = [str(n) for n in z3["names"]]
    assert names == [str(n) for n in z4["names"]]
    rows = []
    for i, n in enumerate(names):
        ok = bool(z3["valid"][i]) and bool(z4["valid"][i]) and (only is None or n in only)
        r = dict(module=n, valid3=int(z3["valid"][i]), valid4=int(z4["valid"][i]),
                 fb3=int(z3["fallback"][i]), fb4=int(z4["fallback"][i]),
                 n_pool3=int(z3["n_pool"][i]), n_pool4=int(z4["n_pool"][i]))
        if ok:
            m3, m4 = T.shape_metrics(z3["shape"][i]), T.shape_metrics(z4["shape"][i])
            r.update(maxdf=T.maxdiff(z3["shape"][i], z4["shape"][i]), d_rise=m4[0] - m3[0],
                     d_w99=m4[1] - m3[1], d_fwhm=m4[2] - m3[2], d_tpk=m4[3] - m3[3],
                     d_tpk_npz=float(z4["t_pk"][i] - z3["t_pk"][i]),
                     d_iint=float(z4["iint_wa"][i] / z3["iint_wa"][i] - 1),
                     lam4=float(z4["lam_module"][i]) if "lam_module" in z4.files else np.nan)
        rows.append(r)
    return pd.DataFrame(rows)


def cmd_templates(args):
    from scipy.stats import spearmanr
    OUT.mkdir(parents=True, exist_ok=True)
    geo = _pmr().load_geometry(25800)
    slug = args.run_type
    t = template_table(HERE / f"output/templates_spline/{slug}/templates_{slug}.npz",
                       HERE / f"output/templates_spline_best50_20261026/{slug}/templates_{slug}.npz",
                       only=set(geo.index[geo.r < 200])).set_index("module")   # inner ring: where Ta clusters live
    t.to_csv(OUT / f"template_diff_{slug}.tsv", sep="\t", float_format="%.5f")
    ok = t.dropna(subset=["maxdf"])
    fbmask = (ok.fb3 == 0) & (ok.fb4 == 0)
    print(f"{len(t)} modules, {len(ok)} valid in both, {int(fbmask.sum())} non-fallback in both")
    for c in ("maxdf", "d_rise", "d_w99", "d_fwhm", "d_tpk", "d_tpk_npz", "d_iint"):
        x = ok[c]
        print(f"  {c:10s} median {x.median():+.4f}  p16/p84 {x.quantile(.16):+.4f}/{x.quantile(.84):+.4f}  |x| median {x.abs().median():.4f}  max|x| {x.abs().max():.4f}")
    if slug != "2p24GeV_Ta_X17":
        return
    mt = pd.read_csv(OUT / f"modules_{args.path}.tsv", sep="\t")
    out = []
    for run in TA_RUNS:
        m = mt[mt.run == run].set_index("module").join(ok, how="left")
        for qty, col in (("module delta (CB, pp)", "delta"), ("raw-energy shift SD (%)", "raw_shift_sd"),
                         ("post-energy shift SD (%)", "shift_sd")):
            s = m.dropna(subset=[col])
            for tc in ("maxdf", "d_rise", "d_w99", "d_tpk", "d_iint"):
                s2 = s.dropna(subset=[tc])
                x = s2[tc].abs() if tc != "d_iint" else s2[tc]
                x = s2[tc].abs()
                rho, p = spearmanr(x, s2[col]) if len(s2) > 4 else (np.nan, np.nan)
                out.append(dict(run=run, y=qty, x="|" + tc + "|", n=len(s2), rho=rho, p=p))
    o = pd.DataFrame(out)
    o.to_csv(OUT / "template_spearman.tsv", sep="\t", index=False, float_format="%.4f")
    print(o.round(3).to_string(index=False))


# ---------------------------------------------------------------------------------------------
# (iii) gains
# ---------------------------------------------------------------------------------------------
G_B50 = HERE / "output/percrystal_best50_20261026"
G_SPL = HERE / "output/percrystal_spl"


def load_gains(d: Path, run: int, path="legacy") -> pd.DataFrame:
    t = pd.read_csv(d / f"{run:06d}" / f"gains_{run:06d}_all.tsv", sep="\t", comment="#")
    return t[t.path == path].set_index("module")


def cmd_gains(args):
    from scipy.stats import spearmanr
    OUT.mkdir(parents=True, exist_ok=True)
    mt = pd.read_csv(OUT / f"modules_{args.path}.tsv", sep="\t")
    rows = []
    for run in (*TA_RUNS, 25228, 25320, 24386):    # Ta + two 2.24 GeV H2 + one 3.49 GeV H2 for reference
        a, b = load_gains(G_SPL, run), load_gains(G_B50, run)
        j = a.join(b, lsuffix="_spl", rsuffix="_b50", how="inner")
        free = j[(j.prior_frac_spl < 0.5) & (j.prior_frac_b50 < 0.5)].copy()
        free["ratio"] = free.g_b50 / free.g_spl
        med = free.ratio.median()
        free["dev"] = free.ratio / med - 1
        mad = 1.4826 * (free.dev - free.dev.median()).abs().median()
        out = free[(free.dev - free.dev.median()).abs() > 4 * mad]
        r = dict(run=run, n_free=len(free), median_ratio=med, rms_pct=100 * free.dev.std(), mad_pct=100 * mad,
                 n_out4=len(out), max_abs_dev_pct=100 * free.dev.abs().max())
        if run in TA_RUNS:
            m = mt[mt.run == run].set_index("module")
            # 3x3 neighbourhood rms of the ratio deviation around each seed (members enter the cluster sum)
            s3 = {}
            for mod in m.index:
                if mod in j.index:
                    rr, cc = j.loc[mod, ["row_spl", "col_spl"]]
                    nb = free[(free.row_spl - rr).abs().le(1) & (free.col_spl - cc).abs().le(1)]
                    s3[mod] = 100 * np.sqrt((nb.dev ** 2).mean()) if len(nb) else np.nan
            m = m.join(free[["dev"]]).assign(seed_dev_abs=lambda x: 100 * x.dev.abs(), nb_rms=pd.Series(s3))
            for col in ("seed_dev_abs", "nb_rms"):
                s2 = m.dropna(subset=[col, "delta"])
                rho, p = spearmanr(s2[col], s2["delta"]) if len(s2) > 4 else (np.nan, np.nan)
                r[f"spearman_{col}_vs_delta"] = f"{rho:+.2f} (p {p:.2f}, n {len(s2)})"
            r["seed_dev_pct_central8"] = ", ".join(
                f"{mm} {100*m.loc[mm, 'dev']:+.2f}" for mm in m.index[m.r < 80] if np.isfinite(m.loc[mm, "dev"]))
            m.reset_index().to_csv(OUT / f"gains_vs_delta_{run}.tsv", sep="\t", index=False, float_format="%.4f")
        rows.append(r)
    o = pd.DataFrame(rows)
    o.to_csv(OUT / "gains_summary.tsv", sep="\t", index=False, float_format="%.4f")
    for r in rows:
        print({k: (round(v, 4) if isinstance(v, float) else v) for k, v in r.items()})


# ---------------------------------------------------------------------------------------------
# (iv) selection
# ---------------------------------------------------------------------------------------------
B50 = HERE / "output/templates_spline_best50_20261026"
AMP_BINS = (0, 600, 1200, 1e9)
AMP_LAB = ("<600", "600-1200", ">1200")


def _first_int(rx, txt, default=0):
    m = re.search(rx, txt)
    return int(m.group(1).replace(",", "")) if m else default


def parse_select_log(p: Path) -> dict:
    t = p.read_text(errors="replace")
    return dict(
        c3b=_first_int(r"c3b_mf_onset\s+(\d+)", t),
        strict_rej=_first_int(r"mf-strict: (\d+) c3b survivors rejected", t),
        mult_ev_rej=_first_int(r"max-modules-per-event \d+ \(A >= 50 ADC\): .*?; (\d+) of \d+ events rejected", t),
        mult_events=_first_int(r"max-modules-per-event \d+ \(A >= 50 ADC\): .*?; \d+ of (\d+) events rejected", t),
        mult_rows=_first_int(r"(\d+) channel rows dropped", t),
        emax_rej=_first_int(r"(\d+) rejected above E_max", t),
        bit24=_first_int(r"exclude-trigger-bits \S+: (\d+) of", t),
        stored=_first_int(r"best-K \d+: (\d+) stored", t),
        full=_first_int(r"stored, (\d+) of 1152 modules full", t),
        occ_mean=float(re.search(r"per-event count mean ([\d.]+)", t).group(1)) if "per-event count mean" in t else np.nan)


def pool_table(slug: str, per_run=False, n=50) -> pd.DataFrame:
    """The build pool: top-n by mf_energy of the module's tier >= 1 band-b0 candidates (b0 = 300-2500 ADC),
    over all splits of the run type (per_run: of each source run on its own); `run` = source run; with r."""
    import glob
    P = _pmr()
    rows = []
    for f in sorted(glob.glob(str(B50 / slug / "pulses" / "*.npz"))):
        z = np.load(f, allow_pickle=True)
        k = z["band_b0"] & (z["tier"] >= 1)
        rows.append(pd.DataFrame({"run": int(z["run"]), "name": z["name"][k], "amp": z["mf_amp"][k],
                                  "E": z["mf_energy"][k]}))
    d = pd.concat(rows)
    geo = P.load_geometry(25800)
    top = d.sort_values("E", ascending=False).groupby(["run", "name"] if per_run else ["name"]).head(n).copy()
    top["r"] = top.name.map(geo.r)
    top["bin"] = pd.cut(top.amp, AMP_BINS, labels=AMP_LAB, right=False)
    return top


def cmd_select(args):
    OUT.mkdir(parents=True, exist_ok=True)
    # --- cut counts per Ta split / run, from the shipped select logs
    rows = []
    for f in sorted((B50 / "logs").glob("select_02[568]*_*.log")):
        m = re.match(r"select_0*(\d+)_(\d+)\.log", f.name)
        run = int(m.group(1))
        if run not in TA_RUNS:
            continue
        rows.append(dict(run=run, split=int(m.group(2)), **parse_select_log(f)))
    L = pd.DataFrame(rows)
    L.to_csv(OUT / "select_cuts_per_split.tsv", sep="\t", index=False)
    g = L.groupby("run").agg(splits=("split", "size"), c3b=("c3b", "sum"), strict_rej=("strict_rej", "sum"),
                             mult_ev_rej=("mult_ev_rej", "sum"), mult_events=("mult_events", "sum"),
                             emax_rej=("emax_rej", "sum"), bit24=("bit24", "sum"), stored=("stored", "sum"),
                             full=("full", "mean"), occ=("occ_mean", "mean"))
    g["strict_pct_of_c3b"] = 100 * g.strict_rej / g.c3b
    g["mult_pct_events"] = 100 * g.mult_ev_rej / g.mult_events
    print("cut counts per Ta run (sum over splits):")
    print(g.round(2).to_string())
    g.to_csv(OUT / "select_cuts_per_run.tsv", sep="\t", float_format="%.3f")
    # --- pool amplitude composition by radius
    out = []
    for slug, tag in (("2p24GeV_Ta_X17", "Ta"), ("2p24GeV_H2", "H2 2.24")):
        top = pool_table(slug, per_run=args.per_run_pool)
        top["reg"] = np.where(top.r < 200, "r<200 mm", "r>=200 mm")
        for (run, reg), s in top.groupby(["run", "reg"]):
            fr = s.bin.value_counts(normalize=True)
            out.append(dict(type=tag, run=run, region=reg, n=len(s), n_mod=s.name.nunique(),
                            **{f"frac{l}": 100 * fr.get(l, 0) for l in AMP_LAB}))
        for run, s in top.groupby("run"):
            fr = s.bin.value_counts(normalize=True)
            out.append(dict(type=tag, run=run, region="all", n=len(s), n_mod=s.name.nunique(),
                            **{f"frac{l}": 100 * fr.get(l, 0) for l in AMP_LAB}))
    o = pd.DataFrame(out)
    o.to_csv(OUT / "select_pool_amp.tsv", sep="\t", index=False, float_format="%.2f")
    print("\ntop-50 pool amplitude fractions (%) per source run and region (pool as built, all runs of the type unless --per-run-pool):")
    print(o.round(1).to_string(index=False))


def cmd_selvar(args):
    """Compare the per-split select output of relaxed-cut variants (select_variants/<name>_<run>_<split>.npz)
    with the shipped B50 select of the same split: stored candidates by amplitude bin and region, and the
    candidates new to the variant (key event_number, midx)."""
    P = _pmr()
    geo = P.load_geometry(25800)
    V = OUT / "select_variants"
    rows = []

    def load(p):
        z = np.load(p, allow_pickle=True)
        d = pd.DataFrame({"ev": z["event_number"], "mi": z["midx"], "name": z["name"], "amp": z["mf_amp"],
                          "E": z["mf_energy"], "b0": z["band_b0"], "tier": z["tier"], "npk": z["mf_npk"]})
        d["r"] = d.name.map(geo.r)
        d["reg"] = np.where(d.r < 200, "r<200", "r>=200")
        d["bin"] = pd.cut(d.amp, AMP_BINS, labels=AMP_LAB, right=False)
        d["pool"] = d.b0 & (d.tier >= 1)
        return d
    for p in sorted(V.glob("*_*_*.npz")):
        name, run, split = p.stem.split("_")
        base = load(B50 / "2p24GeV_Ta_X17/pulses" / f"0{run}_{split}.npz")
        var = load(p)
        key = lambda d: set(zip(d.ev, d.mi))
        new = ~pd.Series(list(zip(var.ev, var.mi))).isin(key(base)).values
        gone = ~pd.Series(list(zip(base.ev, base.mi))).isin(key(var)).values
        for tag, d in (("base", base), ("variant", var)):
            for reg, s in d[d.pool].groupby("reg"):
                fr = s.bin.value_counts(normalize=True)
                rows.append(dict(variant=name, run=int(run), what=tag, region=reg, n_pool_cand=len(s),
                                 **{f"frac{l}": 100 * fr.get(l, 0) for l in AMP_LAB}))
        nn = var[new]
        rows.append(dict(variant=name, run=int(run), what="new_in_variant", region="all", n_pool_cand=len(nn),
                         **{f"frac{l}": 100 * (nn.bin == l).mean() if len(nn) else np.nan for l in AMP_LAB}))
        rows.append(dict(variant=name, run=int(run), what="lost_vs_base", region="all", n_pool_cand=int(gone.sum()),
                         **{f"frac{l}": np.nan for l in AMP_LAB}))
        print(f"{name} {run}.{split}: stored base {len(base)} variant {len(var)}; new {int(new.sum())} "
              f"(amp <600 {int((nn.bin == '<600').sum())}, 600-1200 {int((nn.bin == '600-1200').sum())}, "
              f">1200 {int((nn.bin == '>1200').sum())}; r<200 {int((nn.r < 200).sum())}); lost {int(gone.sum())}")
    o = pd.DataFrame(rows)
    o.to_csv(OUT / "select_variants_summary.tsv", sep="\t", index=False, float_format="%.2f")
    print(o.round(1).to_string(index=False))


def cmd_runs(args):
    """Per-run (never pooled) paired table from compare_scan_dirs *_splits.tsv files: sigma_E/E (CB post, %),
    Delta vs R3 = mean +- SE over the splits of the run (SD/sqrt(n)), legacy and gated; the promotion rule
    |Delta| <= 1 SE (legacy CB, per run) is evaluated for each candidate."""
    cmp_dir = HERE / "output/resolution_fit_20261028_ta/compare"
    rows = []
    for v in args.variants:
        t = pd.read_csv(cmp_dir / f"{v}_vs_r3_splits.tsv", sep="\t")
        t = t[t.type.str.contains("Ta")]
        for run in TA_RUNS:
            for path in ("legacy", "gated"):
                s = t[(t.run == run) & (t.path == path)]
                n = len(s)
                r3 = s["R3_CB_post"].values
                row = dict(variant=v, run=run, path=path, n=n, R3=r3.mean())
                for lab in ("R4", f"Fix_{v}"):
                    x = s[f"{lab}_CB_post"].values
                    d = x - r3
                    se = d.std(ddof=1) / np.sqrt(n)
                    row[lab] = x.mean()
                    row[f"d_{lab}"] = d.mean()
                    row[f"se_{lab}"] = se
                    row[f"within1SE_{lab}"] = "yes" if abs(d.mean()) <= se else "no"
                dd = s[f"Fix_{v}_CB_post"].values - s["R4_CB_post"].values
                row["d_fix_minus_R4"] = dd.mean()
                row["se_fix_minus_R4"] = dd.std(ddof=1) / np.sqrt(n)
                rows.append(row)
    o = pd.DataFrame(rows)
    o.to_csv(OUT / "runs_table.tsv", sep="\t", index=False, float_format="%.4f")
    for v in args.variants:
        print(f"\n== {v}: sigma_E/E CB post (%), delta vs R3 = mean +- SE (pp), per run")
        print("run  path    n   R3     R4     d(R4-R3)          within1SE  Fix    d(Fix-R3)         within1SE  d(Fix-R4)")
        for _, r in o[o.variant == v].iterrows():
            print(f"{r.run} {r.path:6s} {r.n:3d}  {r.R3:.3f}  {r.R4:.3f}  {r.d_R4:+.3f} +- {r.se_R4:.3f}  {r['within1SE_R4']:>3s}      "
                  f"{r[f'Fix_{v}']:.3f}  {r[f'd_Fix_{v}']:+.3f} +- {r[f'se_Fix_{v}']:.3f}  {r[f'within1SE_Fix_{v}']:>3s}      "
                  f"{r.d_fix_minus_R4:+.3f} +- {r.se_fix_minus_R4:.3f}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("modules", help="(i) per-run, per-module paired R4 - R3")
    p.add_argument("--boot", type=int, default=100)
    p.add_argument("--min-n", type=int, default=150)
    p.add_argument("--path", choices=("legacy", "gated"), default="legacy")
    p.add_argument("--cand", default=str(R4), help="candidate scan dir compared with R3 (default: R4 = B50 scan)")
    p.add_argument("--tag", default="", help="suffix of the output tables (default none: the R4 tables used by templates/gains)")
    p.set_defaults(fn=cmd_modules)
    p = sub.add_parser("gains", help="(iii) B50 / spl gain ratios and their relation to the module delta")
    p.add_argument("--path", choices=("legacy", "gated"), default="legacy")
    p.set_defaults(fn=cmd_gains)
    p = sub.add_parser("select", help="(iv) Ta pool amplitude composition, cut counts per run")
    p.add_argument("--per-run-pool", action="store_true", help="pool of each source run on its own (default: as built, all runs of the type; run = source run of the pulse)")
    p.set_defaults(fn=cmd_select)
    p = sub.add_parser("selvar", help="(iv) relaxed-cut select variants vs the shipped B50 select (one split each)")
    p.set_defaults(fn=cmd_selvar)
    p = sub.add_parser("runs", help="per-run paired table R3 / R4 / fix from compare_scan_dirs splits tsv")
    p.add_argument("--variants", nargs="+", default=["phase50", "n200"])
    p.set_defaults(fn=cmd_runs)
    p = sub.add_parser("templates", help="(ii) R3 vs B50 Ta template differences; Spearman with (i)")
    p.add_argument("--run-type", default="2p24GeV_Ta_X17",
                   help="bank slug (other than Ta: differences only, no Spearman; baseline for the Ta numbers)")
    p.add_argument("--path", choices=("legacy", "gated"), default="legacy")
    p.set_defaults(fn=cmd_templates)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
