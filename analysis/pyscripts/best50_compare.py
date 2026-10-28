#!/usr/bin/env python3
"""Comparison tables for the best-50 bank pipeline (plans/plan_best50_pipeline_20261026.md).

Subcommands (each writes <out>.md and <out>.tsv; numpy/pandas only):
  pileup  --ref DIR --new DIR [--arms L=DIR ...] --out PREFIX
          joins <slug>/eff_<slug>.json grids (benchmark_pileup_mf.py) per block;
          also writes PREFIX_cells.tsv (cell level)
  timing  --new TSV --ref TSV --out PREFIX
          module_timing_rf.tsv (long: module id type ...) or module_timing_veto.tsv (wide: sigma_<t>)
  rfcalib --new TSV --ref TSV --out PREFIX
          rf_calib_modules.tsv joined on (run, module); P_ns from rf_calib_runs.tsv next to the file
  cluster --new DIR --ref DIR --out PREFIX
          DIR/cluster_timing.tsv joined on (run, type, path, definition, tier)

Deltas are always new - ref.
"""
import argparse
import json
import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

SLUGS = ["0p73GeV_H2", "2p24GeV_C", "2p24GeV_H2", "2p24GeV_Ta_X17", "3p49GeV_H2"]
# (key in json, sub-key, window?, label)
BLOCKS = [
    ("headline", False, "r_model 0.3 (headline)"),
    ("headline", True, "r_model 0.3, window rule 40:0.16"),
    ("second", False, "r_model 0.1"),
    ("second", True, "r_model 0.1, window rule 40:0.16"),
]
P_DEFAULT = 4.0077


# ---------------------------------------------------------------- helpers
def write_out(prefix, md_lines, df):
    prefix = Path(prefix)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    Path(str(prefix) + ".md").write_text("\n".join(md_lines) + "\n")
    df.to_csv(str(prefix) + ".tsv", sep="\t", index=False, float_format="%.6g")
    print(f"[done] {prefix}.md  {prefix}.tsv  ({len(df)} rows)")


def md_table(headers, rows):
    out = ["| " + " | ".join(headers) + " |",
           "|" + "|".join("---:" if i else "---" for i in range(len(headers))) + "|"]
    out += ["| " + " | ".join(str(c) for c in r) + " |" for r in rows]
    return out


def fnum(x, nd=3, sign=False):
    if x is None or (isinstance(x, float) and not np.isfinite(x)):
        return "nan"
    return f"{x:+.{nd}f}" if sign else f"{x:.{nd}f}"


def fold(x, P):
    return (np.asarray(x, float) + P / 2) % P - P / 2


def wmedian(v, w):
    v = np.asarray(v, float)
    w = np.asarray(w, float)
    m = np.isfinite(v) & np.isfinite(w) & (w > 0)
    if not m.any():
        return np.nan
    v, w = v[m], w[m]
    o = np.argsort(v)
    v, w = v[o], w[o]
    c = np.cumsum(w)
    return float(v[np.searchsorted(c, 0.5 * c[-1])])


def med_pct(x):
    x = np.asarray(x, float)
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan, np.nan
    return tuple(float(v) for v in np.percentile(x, [50, 16, 84]))


# ---------------------------------------------------------------- pileup
def load_eff(d, slug):
    p = Path(d) / slug / f"eff_{slug}.json"
    if not p.exists():
        return None
    j = json.loads(p.read_text())
    md = (Path(d) / slug / f"eff_{slug}.md")
    min_sep = None
    if md.exists():
        m = re.search(r"MIN_SEP\s+([0-9.]+)\s*ns", md.read_text())
        min_sep = float(m.group(1)) if m else None
    return j, min_sep


def grid(j, key, win):
    """dict metric -> (nd, nr) array in % for eff/n2, per event for fake."""
    res = j["results_window" if win else "results"][key]["cells"]
    dts, rs = j["dt_ns"], j["ratios"]
    out = {k: np.full((len(dts), len(rs)), np.nan) for k in ("eff", "fake", "n2mf", "n2wa")}
    for c in res:
        i, k = dts.index(c["dt_ns"]), rs.index(c["ratio"])
        out["eff"][i, k] = 100 * c["eff_match"]
        out["fake"][i, k] = c["fake_rate"]
        out["n2mf"][i, k] = 100 * c["n2_mf"]
        out["n2wa"][i, k] = 100 * c["n2_wa"]
    return out


def pooled(new, ref):
    d = new - ref
    sel = ref < 100 - 1e-9
    return (float(np.nanmean(d)), int(d.size),
            float(np.nanmean(d[sel])) if sel.any() else np.nan, int(sel.sum()))


def grid_md(dts, rs, fmt_cell):
    rows = [[f"{dt:g}"] + [fmt_cell(i, k) for k in range(len(rs))] for i, dt in enumerate(dts)]
    return md_table(["dt (ns)"] + [f"{r:g}" for r in rs], rows)


def cmd_pileup(a):
    arms = {}
    for s in a.arms or []:
        lab, _, d = s.partition("=")
        arms[lab] = d
    new_lab = "B"
    L = ["# Pile-up detection efficiency: new bank vs reference", "",
         f"- ref: `{a.ref}`", f"- new (arm {new_lab}): `{a.new}`"]
    for lab, d in arms.items():
        L.append(f"- arm {lab}: `{d}`")
    L += ["- eff_match, n2_mf, n2_wa in % of events (50 per cell); fake_rate in unmatched mf pulses per event.",
          "- eff cells are `new (delta)` with delta = new - ref in points; pooled values are the mean of the cell deltas.",
          "- pooled `ref<100` = only the cells whose ref eff_match is below 100 %.", ""]
    cells, pool = [], []
    data = {}
    for slug in SLUGS:
        r = load_eff(a.ref, slug)
        n = load_eff(a.new, slug)
        if r is None or n is None:
            L.append(f"(slug {slug} missing in ref or new, skipped)")
            continue
        data[slug] = (r, n, {lab: load_eff(d, slug) for lab, d in arms.items()})
    # ---- top summary
    summ_rows = []
    for slug, ((jr, _), (jn, _), ar) in data.items():
        gr, gn = grid(jr, "headline", False), grid(jn, "headline", False)
        pb = pooled(gn["eff"], gr["eff"])
        d = (gn["eff"] - gr["eff"])
        idx = np.argsort(d, axis=None)[:3]
        worst = []
        for ix in idx:
            i, k = np.unravel_index(ix, d.shape)
            if not d[i, k] < 0:
                continue
            worst.append(f"dt {jr['dt_ns'][i]:g} r {jr['ratios'][k]:g}: {gr['eff'][i, k]:.0f}->{gn['eff'][i, k]:.0f} ({d[i, k]:+.0f})")
        row = [slug, f"{pb[0]:+.2f}", f"{pb[2]:+.2f} (n={pb[3]})"]
        for lab in arms:
            x = ar[lab]
            if x is None:
                row += ["n/a", "n/a"]
            else:
                p = pooled(grid(x[0], "headline", False)["eff"], gr["eff"])
                row += [f"{p[0]:+.2f}", f"{p[2]:+.2f}"]
        row.append("; ".join(worst) if worst else "none negative")
        summ_rows.append(row)
    hdr = ["slug", "B pooled all", "B pooled ref<100"]
    for lab in arms:
        hdr += [f"{lab} pooled all", f"{lab} pooled ref<100"]
    hdr.append("B worst 3 cells (ref->new, delta)")
    L += ["## Summary: pooled delta eff_match (points), headline block r_model 0.3", ""]
    L += md_table(hdr, summ_rows) + [""]
    # ---- per slug
    for slug, ((jr, ms_r), (jn, ms_n), ar) in data.items():
        sr, sn = jr["settings"], jn["settings"]
        L += [f"## {slug}", "",
              f"- template sha1 ref `{sr['template_sha1'][:12]}`, new `{sn['template_sha1'][:12]}`"
              f" (modules from file {sr['n_module_from_file']} / {sn['n_module_from_file']},"
              f" fallback {sr['n_fallback_modules']} / {sn['n_fallback_modules']})",
              f"- MIN_SEP_NS ref {ms_r:g}, new {ms_n:g}" if ms_r and ms_n else "- MIN_SEP_NS not parsed",
              f"- r_model headline {sn['r_models']['headline']:g}, second {sn['r_models']['second']:g};"
              f" tolerance {sn['tolerance_ns']:g} ns (ref {sr['tolerance_ns']:g});"
              f" nsigma {sn['nsigma']:g}, min_amp {sn['min_amp']:g}, kmax {sn['kmax']}",
              f"- synth: ref `{sr['synth_dir']}`, new `{sn['synth_dir']}`", ""]
        dts, rs = jr["dt_ns"], jr["ratios"]
        for key, win, label in BLOCKS:
            gr, gn = grid(jr, key, win), grid(jn, key, win)
            eq = bool(np.array_equal(np.nan_to_num(gr["n2wa"], nan=-1), np.nan_to_num(gn["n2wa"], nan=-1)))
            p = pooled(gn["eff"], gr["eff"])
            L += [f"### {slug}: {label}", "",
                  "eff_match new (delta new - ref), % :", ""]
            L += grid_md(dts, rs, lambda i, k: f"{gn['eff'][i, k]:.0f} ({gn['eff'][i, k] - gr['eff'][i, k]:+.0f})"
                         if gn['eff'][i, k] != gr['eff'][i, k] else f"{gn['eff'][i, k]:.0f} (0)") + [""]
            L += ["fake_rate new (ref), per event:", ""]
            L += grid_md(dts, rs, lambda i, k: f"{gn['fake'][i, k]:.2f} ({gr['fake'][i, k]:.2f})") + [""]
            L += ["n2_mf new (ref), %:", ""]
            L += grid_md(dts, rs, lambda i, k: f"{gn['n2mf'][i, k]:.0f} ({gr['n2mf'][i, k]:.0f})") + [""]
            L.append(f"- pooled delta eff_match: all {p[1]} cells {p[0]:+.3f} points; "
                     f"cells with ref < 100 % ({p[3]} cells) {p[2]:+.3f} points")
            L.append(f"- n2_wa new == ref cell by cell: {'YES' if eq else 'NO'}")
            pool.append(dict(slug=slug, block=label, arm=new_lab, n_cells=p[1], pooled_d_all=p[0],
                             n_ref_lt100=p[3], pooled_d_ref_lt100=p[2], n2wa_equal=int(eq)))
            for lab in arms:
                x = ar[lab]
                if x is None:
                    L.append(f"- arm {lab}: missing")
                    continue
                ga = grid(x[0], key, win)
                pa = pooled(ga["eff"], gr["eff"])
                eqa = bool(np.array_equal(np.nan_to_num(gr["n2wa"], nan=-1), np.nan_to_num(ga["n2wa"], nan=-1)))
                L.append(f"- arm {lab} (bank sha1 {x[0]['settings']['template_sha1'][:12]}, synth sha1 "
                         f"{x[0]['settings']['synth_manifest']['template_sha1'][:12]}): pooled delta vs ref all {pa[0]:+.3f}, "
                         f"ref<100 {pa[2]:+.3f} points; n2_wa == ref: {'YES' if eqa else 'NO (injected waveforms differ in the new synth)'}")
                pool.append(dict(slug=slug, block=label, arm=lab, n_cells=pa[1], pooled_d_all=pa[0],
                                 n_ref_lt100=pa[3], pooled_d_ref_lt100=pa[2], n2wa_equal=int(eqa)))
            L.append("")
            for i, dt in enumerate(dts):
                for k, r in enumerate(rs):
                    cells.append(dict(slug=slug, block=label, arm=new_lab, dt_ns=dt, ratio=r,
                                      eff_ref=gr["eff"][i, k], eff_new=gn["eff"][i, k],
                                      d_eff=gn["eff"][i, k] - gr["eff"][i, k],
                                      fake_ref=gr["fake"][i, k], fake_new=gn["fake"][i, k],
                                      n2mf_ref=gr["n2mf"][i, k], n2mf_new=gn["n2mf"][i, k],
                                      n2wa_ref=gr["n2wa"][i, k], n2wa_new=gn["n2wa"][i, k]))
                    for lab in arms:
                        x = ar[lab]
                        if x is None:
                            continue
                        ga = grid(x[0], key, win)
                        cells.append(dict(slug=slug, block=label, arm=lab, dt_ns=dt, ratio=r,
                                          eff_ref=gr["eff"][i, k], eff_new=ga["eff"][i, k],
                                          d_eff=ga["eff"][i, k] - gr["eff"][i, k],
                                          fake_ref=gr["fake"][i, k], fake_new=ga["fake"][i, k],
                                          n2mf_ref=gr["n2mf"][i, k], n2mf_new=ga["n2mf"][i, k],
                                          n2wa_ref=gr["n2wa"][i, k], n2wa_new=ga["n2wa"][i, k]))
    allwa = all(p["n2wa_equal"] for p in pool if p["arm"] == new_lab)
    L.insert(L.index("## Summary: pooled delta eff_match (points), headline block r_model 0.3"),
             f"- n2_wa of {new_lab} equals ref cell by cell in every slug and block: {'YES' if allwa else 'NO'}\n")
    write_out(a.out, L, pd.DataFrame(pool))
    pd.DataFrame(cells).to_csv(str(a.out) + "_cells.tsv", sep="\t", index=False, float_format="%.6g")


# ---------------------------------------------------------------- timing
def _stats_row(label, new_s, ref_s, new_n, ref_n, dsig, extra=None):
    m, lo, hi = med_pct(dsig)
    return dict(group=label, n_fit_both=int(np.isfinite(dsig).sum()),
                med_sigma_new=float(np.nanmedian(new_s)) if len(new_s) else np.nan,
                med_sigma_ref=float(np.nanmedian(ref_s)) if len(ref_s) else np.nan,
                med_dsigma=m, p16_dsigma=lo, p84_dsigma=hi,
                wmed_sigma_new=wmedian(new_s, new_n), wmed_sigma_ref=wmedian(ref_s, ref_n),
                **(extra or {}))


def cmd_timing(a):
    new = pd.read_csv(a.new, sep="\t")
    ref = pd.read_csv(a.ref, sep="\t")
    rows, L = [], ["# Module timing: new vs ref", "", f"- new: `{a.new}`", f"- ref: `{a.ref}`",
                   "- sigma in ns; Delta = new - ref over modules with fit_ok in both files; p16/p84 are percentiles of the per-module delta.",
                   "- N-weighted median: weights = N of the respective file."]
    if "type" in new.columns:
        kind = "rf"
        df = new.merge(ref, on=["module", "type"], suffixes=("_n", "_r"))
        L += [f"- kind: RF-referenced (long). Offset shift = median of (new - ref) offset folded into (-P/2, P/2], P = {a.period:g} ns.",
              f"- joined module-type rows: {len(df)} (new {len(new)}, ref {len(ref)})", ""]
        groups = [(t, g) for t, g in df.groupby("type", sort=True)]
        sig, ok, nn = "sigma", "fit_ok", "N"
    else:
        kind = "veto"
        types = [c[len("sigma_"):] for c in new.columns if c.startswith("sigma_")]
        df = new.merge(ref, on=["module", "id"], suffixes=("_n", "_r"))
        L += [f"- kind: veto-referenced (wide), types {types}. Mean shift = median of (new - ref) mean_<t>.",
              f"- joined modules: {len(df)} (new {len(new)}, ref {len(ref)})", ""]
        groups = []
        for t in types:
            g = pd.DataFrame({"module": df["module"],
                              "sigma_n": df[f"sigma_{t}_n"], "sigma_r": df[f"sigma_{t}_r"],
                              "fit_ok_n": df[f"fit_ok_{t}_n"], "fit_ok_r": df[f"fit_ok_{t}_r"],
                              "N_n": df[f"N_{t}_n"], "N_r": df[f"N_{t}_r"],
                              "mean_n": df[f"mean_{t}_n"], "mean_r": df[f"mean_{t}_r"]})
            groups.append((t, g))
    for t, g in groups:
        sn_, sr_ = g["sigma_n"].astype(float), g["sigma_r"].astype(float)
        both = (g["fit_ok_n"] == 1) & (g["fit_ok_r"] == 1) & np.isfinite(sn_) & np.isfinite(sr_)
        gb = g[both]
        d = (gb["sigma_n"] - gb["sigma_r"]).to_numpy(float)
        extra = dict(n_modules=len(g), n_fit_new=int(((g["fit_ok_n"] == 1)).sum()),
                     n_fit_ref=int(((g["fit_ok_r"] == 1)).sum()))
        if kind == "rf":
            for o in ("offset_raw", "offset_tcal"):
                extra[f"med_d_{o}"] = float(np.nanmedian(fold(gb[f"{o}_n"] - gb[f"{o}_r"], a.period))) if len(gb) else np.nan
        else:
            extra["med_d_mean"] = float(np.nanmedian(gb["mean_n"] - gb["mean_r"])) if len(gb) else np.nan
        r = _stats_row(t, gb["sigma_n"].to_numpy(float), gb["sigma_r"].to_numpy(float),
                       gb["N_n"].to_numpy(float), gb["N_r"].to_numpy(float), d, extra)
        rows.append(r)
    out = pd.DataFrame(rows)
    out.insert(0, "kind", kind)
    ex = [c for c in out.columns if c.startswith("med_d_o") or c == "med_d_mean"]
    hdr = ["type", "n_fit both (new/ref fit)", "med sigma new", "med sigma ref", "med dsigma [p16, p84]",
           "Nw-med sigma new", "Nw-med sigma ref"] + ex
    tab = []
    for _, r in out.iterrows():
        tab.append([r["group"], f"{r['n_fit_both']} ({r['n_fit_new']}/{r['n_fit_ref']} of {r['n_modules']})",
                    fnum(r["med_sigma_new"]), fnum(r["med_sigma_ref"]),
                    f"{fnum(r['med_dsigma'], 4, True)} [{fnum(r['p16_dsigma'], 4, True)}, {fnum(r['p84_dsigma'], 4, True)}]",
                    fnum(r["wmed_sigma_new"]), fnum(r["wmed_sigma_ref"])] + [fnum(r[c], 4, True) for c in ex])
    L += md_table(hdr, tab)
    write_out(a.out, L, out)


# ---------------------------------------------------------------- rfcalib
def cmd_rfcalib(a):
    new = pd.read_csv(a.new, sep="\t", dtype={"run": str})
    ref = pd.read_csv(a.ref, sep="\t", dtype={"run": str})
    P = {}
    for src in (a.ref, a.new):  # new overrides ref
        p = Path(src).parent / "rf_calib_runs.tsv"
        if p.exists():
            r = pd.read_csv(p, sep="\t", dtype={"run": str})
            P.update(dict(zip(r["run"], r["P_ns"])))
    df = new.merge(ref, on=["run", "module"], suffixes=("_n", "_r"))
    meas = {"e100", "allE"}
    df["P"] = df["run"].map(P).fillna(P_DEFAULT)
    df["both_meas"] = df["tier_n"].isin(meas) & df["tier_r"].isin(meas)
    df["dsig"] = np.where(df["both_meas"], df["sigma_ns_n"] - df["sigma_ns_r"], np.nan)
    df["doff"] = np.where(df["both_meas"], fold(df["offset_ns_n"] - df["offset_ns_r"], df["P"].to_numpy()), np.nan)
    df["tier_changed"] = df["tier_n"] != df["tier_r"]
    L = ["# RF calibration per (run, module): new vs ref", "", f"- new: `{a.new}`", f"- ref: `{a.ref}`",
         f"- joined (run, module) pairs: {len(df)} (new {len(new)}, ref {len(ref)})",
         "- measured pair = tier in {e100, allE} in both files; dsigma and doffset are over measured pairs with finite values;",
         "  doffset = new - ref folded into (-P/2, P/2] with P_ns of the run from rf_calib_runs.tsv.",
         "- tier changed = fraction of all joined pairs of the group whose tier differs.",
         "- rows: per run, per type over all runs, and per type restricted to pairs with tier e100 in both (type_e100_both).", ""]

    def summarize(g, label):
        mk = g["both_meas"]
        gs = g[mk & np.isfinite(g["dsig"])]
        return dict(group=label, n_pairs=len(g), n_measured=int(mk.sum()),
                    med_sigma_new=float(np.nanmedian(gs["sigma_ns_n"])) if len(gs) else np.nan,
                    med_sigma_ref=float(np.nanmedian(gs["sigma_ns_r"])) if len(gs) else np.nan,
                    med_dsigma=float(np.nanmedian(gs["dsig"])) if len(gs) else np.nan,
                    med_doffset=float(np.nanmedian(g.loc[mk, "doff"])) if mk.any() else np.nan,
                    frac_tier_changed=float(g["tier_changed"].mean()) if len(g) else np.nan)

    rows = []
    for (run, typ), g in df.groupby(["run", "type_n"], sort=True):
        rows.append(dict(level="run", run=run, type=typ, **summarize(g, f"{run} {typ}")))
    for typ, g in df.groupby("type_n", sort=True):
        rows.append(dict(level="type_all_runs", run="ALL", type=typ, **summarize(g, typ)))
        ge = g[(g["tier_n"] == "e100") & (g["tier_r"] == "e100")]
        rows.append(dict(level="type_e100_both", run="ALL", type=typ, **summarize(ge, typ + " e100")))
    out = pd.DataFrame(rows).drop(columns=["group"])
    hdr = ["level", "run", "type", "pairs", "measured", "med sigma new", "med sigma ref", "med dsigma",
           "med doffset", "frac tier changed"]
    tab = [[r["level"], r["run"], r["type"], r["n_pairs"], r["n_measured"], fnum(r["med_sigma_new"]),
            fnum(r["med_sigma_ref"]), fnum(r["med_dsigma"], 4, True), fnum(r["med_doffset"], 4, True),
            fnum(r["frac_tier_changed"], 4)] for _, r in out.iterrows()]
    L += md_table(hdr, tab)
    write_out(a.out, L, out)


# ---------------------------------------------------------------- cluster
def cmd_cluster(a):
    key = ["run", "type", "path", "definition", "tier"]
    rd = lambda d: pd.read_csv(Path(d) / "cluster_timing.tsv", sep="\t", dtype={"run": str})
    new, ref = rd(a.new), rd(a.ref)
    df = new.merge(ref, on=key, how="outer", suffixes=("_n", "_r"))
    both = (df["fit_ok_n"] == 1) & (df["fit_ok_r"] == 1)
    df["dsigma_ns"] = np.where(both, df["sigma_ns_n"] - df["sigma_ns_r"], np.nan)
    tiers = ["all", "1", "2-4", "5+"]
    df["tier"] = df["tier"].astype(str)
    L = ["# sigma_t per cluster: new vs ref", "", f"- new: `{a.new}/cluster_timing.tsv`",
         f"- ref: `{a.ref}/cluster_timing.tsv`",
         "- sigma in ns (RF-referenced fit_unwrapped); delta = new - ref where fit_ok in both. The ref may cover fewer splits than new.",
         "- frac_off = fraction of cluster members with an RF offset (new / ref).", ""]

    def cell(r, s):
        return f"{r['sigma_ns_' + s]:.3f} ± {r['sigma_err_ns_' + s]:.3f}" if np.isfinite(r["sigma_ns_" + s]) else "nan"

    def fo(r):
        return f"{fnum(r['frac_members_with_offset_n'], 3)} / {fnum(r['frac_members_with_offset_r'], 3)}"

    for typ in sorted(df["type"].dropna().unique()):
        for path in sorted(df.loc[df["type"] == typ, "path"].dropna().unique()):
            g = df[(df["type"] == typ) & (df["path"] == path) & (df["run"] == "ALL")]
            tab = []
            for dfn in ["seed", "ewt"]:
                for t in tiers:
                    r = g[(g["definition"] == dfn) & (g["tier"] == t)]
                    if r.empty:
                        continue
                    r = r.iloc[0]
                    tab.append([dfn, t, fnum(r["N_n"], 0), fnum(r["N_r"], 0), cell(r, "n"), cell(r, "r"),
                                fnum(r["dsigma_ns"], 3, True), fo(r)])
            L += [f"## {typ}, path {path}: pooled over runs", ""]
            L += md_table(["definition", "tier", "N new", "N ref", "sigma new", "sigma ref", "delta", "frac_off new / ref"], tab) + [""]
    for typ in sorted(df["type"].dropna().unique()):
        for path in sorted(df.loc[df["type"] == typ, "path"].dropna().unique()):
            g = df[(df["type"] == typ) & (df["path"] == path) & (df["run"] != "ALL") & (df["tier"] == "all")]
            tab = []
            for _, r in g.sort_values(["definition", "run"]).iterrows():
                tab.append([r["definition"], r["run"], fnum(r["N_n"], 0), fnum(r["N_r"], 0), cell(r, "n"), cell(r, "r"),
                            fnum(r["dsigma_ns"], 3, True), fo(r)])
            L += [f"## {typ}, path {path}: per run, tier all", ""]
            L += md_table(["definition", "run", "N new", "N ref", "sigma new", "sigma ref", "delta", "frac_off new / ref"], tab) + [""]
    cols = key + [c for c in df.columns if c not in key]
    write_out(a.out, L, df[cols].sort_values(key))


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pileup")
    p.add_argument("--ref", required=True)
    p.add_argument("--new", required=True)
    p.add_argument("--arms", nargs="*", default=[], metavar="LABEL=DIR")
    p.add_argument("--out", required=True)
    p.set_defaults(f=cmd_pileup)
    p = sub.add_parser("timing")
    p.add_argument("--new", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--period", type=float, default=P_DEFAULT, help="RF period (ns) to fold offset shifts")
    p.set_defaults(f=cmd_timing)
    p = sub.add_parser("rfcalib")
    p.add_argument("--new", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(f=cmd_rfcalib)
    p = sub.add_parser("cluster")
    p.add_argument("--new", required=True)
    p.add_argument("--ref", required=True)
    p.add_argument("--out", required=True)
    p.set_defaults(f=cmd_cluster)
    a = ap.parse_args()
    a.f(a)


if __name__ == "__main__":
    sys.exit(main())
