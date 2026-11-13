"""cmp_templates.py -- 10-90% rise, 90-10% fall, sigma(t) and pulse jitter of the extracted spline templates.

    python3 cmp_templates.py

Reads BANK_ROOT/<run type>/templates_<run type>.npz (spline_templates_v3; `shape` [module, t_grid],
unit peak, t_grid in ns after onset) for each RUN_TYPES entry; a missing file is reported and skipped.
Counted modules: fallback == 0 (own template) and not dead; single-pulse modules are counted (their rise
and fall come from their own shape).  Dead modules are hatched on the maps; modules on the run-type
average template (fallback != 0) are not counted, since their shape is not their own (grey on the maps).
Per counted module: baseline = shape at t_grid[0], normalise to the peak,
rise = t(90%) - t(10%) on the leading edge and fall = t(10%) - t(90%) after the peak, by linear
interpolation; a module without a crossing is skipped.
Rule: tables and distribution plots go to pulse_templates/timing/, maps to waveform_analysis/figures/.
timing/: timing_summary.tsv (per run type: counts, rise/fall median and sigmaMAD), timing_modules.tsv,
timing_dist.png/.pdf, timing_sigma.tsv (median shape_sigma at t_pk and t50, % of peak, over counted
modules with their own sigma) and timing_sigma.png/.pdf (array-median sigma(t) per run type),
timing_maxdiff.tsv (module, x, y, n_runtypes, rise/fall max - min over run types, rise/fall sample standard
deviation (ddof = 1) over run types, needs >= 3, else NaN; timing_summary.tsv row across_run_types holds its median and p90).
figures/timing_std_map: heat maps of that standard deviation (grey = fewer than 3 run types, hatched = dead in any).
figures/: timing_maxdiff_map (HyCal heat maps of that max - min; >= 2 run types needed, else grey) and
timing_map_<run type> (PbWO4 rise left, fall right; colour limits = 1st-99th percentile pooled over the
run types read), timing_map_all (all run types, one row each, one colour bar per column), sigma_map_<run type> (PbWO4 sigma(t) at the pulse peak left, at the 50% rising-edge point
right; one colour scale per quantity shared by all run types, 0 to the 99th percentile; dead hatched
"////", modules still on the borrowed run-type median sigma (fewer than 50 band pulses) hatched "\\\\") and sigma_map_all (all run types, one row each);
timing/sigma_modules.tsv (module, x, y, per run type sigma at peak, at t50 in % of peak, borrowed flag).
jitter_map_<run type> (three panels: rising-edge time jitter sigma(t) / slope at 10, 50, 90% of the peak,
shared scales, same hatching) and jitter_map_all; timing/jitter_modules.tsv; timing_summary.tsv also carries
the median and 90th percentile of the 50% jitter over own-sigma modules, and of rise_jitter_ns / fall_jitter_ns
over the live modules that have them.  sigma(t) is the template's own continuous shape_sigma (band-pool spread
of every module); only dead and still-borrowed modules are hatched.
rise_jitter_map_<run type> / fall_jitter_map_<run type> (+ _all.png): PbWO4 maps of the per-pulse pulse-to-pulse
jitter (1.4826 MAD, ns) of the 10-90 % rise / 90-10 % fall time of the raw band-pool pulses, one scale per
quantity shared by all run types (0 to the 99th percentile); grey = fewer than 50 band pulses, hatched = dead;
jitter_modules.tsv also lists rise_jitter_ns, fall_jitter_ns, rise_med_ns, fall_med_ns.
Geometry from database/hycal_map.json.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.rcParams.update({
    "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
    "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
    "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})

HERE = Path(__file__).resolve().parent
FIGDIR = HERE / "figures"
BANK_ROOT = HERE / "pulse_templates"
OUT = BANK_ROOT / "timing"
HYCAL_MAP = HERE.parent.parent.parent / "database" / "hycal_map.json"
LO, HI = 0.10, 0.90
RUN_TYPES = ("0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17")
MAP_LIM = 350.0    # heat maps: fixed x and y range is [-MAP_LIM, MAP_LIM] mm


def _cross(t, y, i0, i1, level):
    """First upward (i0 < i1) or downward (i0 > i1) crossing of `level` between
    sample indices i0 and i1; linear interpolation.  NaN if none."""
    idx = np.arange(i0, i1 + (1 if i1 > i0 else -1), 1 if i1 > i0 else -1)
    for a, b in zip(idx[:-1], idx[1:]):
        ya, yb = y[a], y[b]
        if (ya - level) * (yb - level) <= 0 and ya != yb:
            return t[a] + (level - ya) / (yb - ya) * (t[b] - t[a])
    return np.nan


def rise_fall(t, y):
    y = np.asarray(y, dtype=float)
    y = y - y[0]
    pk = int(np.argmax(y))
    if y[pk] <= 0:
        return np.nan, np.nan
    y = y / y[pk]
    # leading edge: walk back from the peak so a noisy pre-onset wiggle cannot
    # fake an early crossing
    r = _cross(t, y, pk, 0, HI)
    r_lo = _cross(t, y, pk, 0, LO)
    rise = r - r_lo
    f_hi = _cross(t, y, pk, len(t) - 1, HI)
    f_lo = _cross(t, y, pk, len(t) - 1, LO)
    return rise, f_lo - f_hi


def jitter(t, y, sg):
    """Rising-edge time jitter sigma(t_f) / |dS/dt(t_f)| in ns at f = 0.1, 0.5, 0.9 of the peak (crossings by
    walking back from the peak, as rise_fall; slope = centred difference on the grid), and the approximation
    sigma(t_50) * rise time / 0.8.  Returns four values, NaN where there is no crossing."""
    y = np.asarray(y, dtype=float) - y[0]
    pk = int(np.argmax(y))
    if y[pk] <= 0:
        return (np.nan,) * 4
    y = y / y[pk]
    slope = np.gradient(y, t)
    tf = [_cross(t, y, pk, 0, f) for f in (LO, 0.5, HI)]
    if not np.all(np.isfinite(tf)):
        return (np.nan,) * 4
    out = [np.interp(x, t, sg) / abs(np.interp(x, t, slope)) for x in tf]
    return (*out, np.interp(tf[1], t, sg) * (tf[2] - tf[0]) / (HI - LO))


def stats(x):
    x = x[np.isfinite(x)]
    if x.size == 0:
        return np.nan, np.nan
    med = np.median(x)
    return med, 1.4826 * np.median(np.abs(x - med))


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    res = {}          # run type -> (rise[], fall[])
    sig = {}          # run type -> (t_grid, array-median sigma, sigma at t_pk, sigma at t50)
    rows, summ = [], []
    permod = {}       # module name -> {run type: (rise, fall)}, finite pairs only
    excl = {}         # run type -> set of dead module names
    jit = {}          # run type -> module name -> (jitter at 10, 50, 90 %, approximation at 50 %, borrowed, dead)
    pj = {}           # run type -> module name -> (rise jitter, fall jitter, rise median, fall median, dead), ns
    sigmod = {}       # run type -> module name -> (sigma at peak %, sigma at t50 %, borrowed, dead)
    for rt in RUN_TYPES:
        fp = BANK_ROOT / rt / f"templates_{rt}.npz"
        if not fp.exists():
            print(f"[{rt}] {fp} not found, run type skipped")
            continue
        d = np.load(fp, allow_pickle=False)
        t, shape, names = d["t_grid"], d["shape"], d["names"]   # load once (npz is lazy)
        bad = d["dead"] == 1
        use = (d["fallback"] == 0) & ~bad
        excl[rt] = {str(n) for n in names[bad]}
        rs, fs, skipped = [], [], []
        for i in np.where(use)[0]:
            r, f = rise_fall(t, shape[i])
            rows.append((rt, str(names[i]), r, f))
            if np.isfinite(r) and np.isfinite(f):
                rs.append(r)
                fs.append(f)
                permod.setdefault(str(names[i]), {})[rt] = (r, f)
            else:
                skipped.append(str(names[i]))
        rs, fs = np.array(rs), np.array(fs)
        res[rt] = (rs, fs)
        mr, sr = stats(rs)
        mf, sf = stats(fs)
        nd, ns = int((d["dead"] == 1).sum()), int((use & (d["single"] == 1)).sum())
        summ.append((rt, rs.size, len(skipped), nd, ns, mr, sr, mf, sf))
        if skipped:
            print(f"[{rt}] skipped (no valid crossing): {' '.join(skipped)}")
        print(f"[{rt}] {int(use.sum())} modules counted ({ns} single-pulse included); excluded: {nd} dead; "
              f"{int(((d['fallback'] != 0) & ~bad).sum())} run-type average not counted")
        own = use & (d["sigma_borrowed"] == 0)
        sg = d["shape_sigma"][own].astype(float)
        at = lambda k: np.array([np.interp(x, t, s_) for x, s_ in zip(d[k][own], sg)])   # noqa: E731
        sig[rt] = (t, np.nanmedian(sg, axis=0), np.nanmedian(at("t_pk")), np.nanmedian(at("t50")), int(own.sum()))
        allsg = d["shape_sigma"].astype(float)
        jit[rt] = {str(n): (*jitter(t, sh, s_), int(br), int(dd))
                   for n, sh, s_, br, dd in zip(names, shape, allsg, d["sigma_borrowed"], d["dead"])}
        pj[rt] = {str(n): (*(float(d[k][i]) for k in ("rise_jitter_ns", "fall_jitter_ns", "rise_med_ns", "fall_med_ns")), int(dd))
                  for i, (n, dd) in enumerate(zip(names, d["dead"]))}
        sigmod[rt] = {str(n): (100 * np.interp(a, t, s_), 100 * np.interp(b, t, s_), int(br), int(dd))
                      for n, a, b, s_, br, dd in zip(names, d["t_pk"], d["t50"], allsg, d["sigma_borrowed"], d["dead"])}
    if not res:
        raise SystemExit(f"no templates_<run type>.npz under {BANK_ROOT}")

    for k, s_ in enumerate(summ):
        j50 = np.array([v[1] for v in jit[s_[0]].values() if not v[4] and np.isfinite(v[1])])
        pv = np.array([[v[0], v[1]] for v in pj[s_[0]].values() if not v[4]])
        pq = [f(pv[np.isfinite(pv[:, c]), c], q) for c in (0, 1) for f, q in ((np.median, None), (np.percentile, 90))]
        summ[k] = (*s_, np.median(j50), np.percentile(j50, 90), pq[0], pq[1], pq[2], pq[3])
    hdr = ["run_type", "n_used", "n_skipped", "n_dead", "n_single", "rise_median_ns", "rise_sigmaMAD_ns",
           "fall_median_ns", "fall_sigmaMAD_ns", "jitter_50_median_ns", "jitter_50_p90_ns",
           "rise_jitter_median_ns", "rise_jitter_p90_ns", "fall_jitter_median_ns", "fall_jitter_p90_ns",
           "rise_std_median_ns", "rise_std_p90_ns", "fall_std_median_ns", "fall_std_p90_ns"]
    with open(OUT / "timing_summary.tsv", "w") as fh:
        fh.write("\t".join(hdr) + "\n")
        for s in summ:
            fh.write("\t".join([s[0]] + [str(v) for v in s[1:5]] + [f"{v:.3f}" for v in s[5:]] + ["nan"] * 4) + "\n")
    with open(OUT / "timing_modules.tsv", "w") as fh:
        fh.write("run_type\tmodule\trise_ns\tfall_ns\n")
        for rt, nm, r, f in rows:
            fh.write(f"{rt}\t{nm}\t{r:.3f}\t{f:.3f}\n")
    with open(OUT / "timing_sigma.tsv", "w") as fh:
        fh.write("run_type\tn_modules\tsigma_at_peak_pct\tsigma_at_t50_pct\n")
        for rt, v in sig.items():
            fh.write(f"{rt}\t{v[4]}\t{100 * v[2]:.3f}\t{100 * v[3]:.3f}\n")

    print("\n" + "  ".join(f"{h:>16}" for h in hdr))
    for s in summ:
        print("  ".join([f"{s[0]:>16}"] + [f"{v:>16d}" for v in s[1:5]] + [f"{v:>16.3f}" for v in s[5:]]))
    print("\nmedian shape_sigma (% of peak): run_type n_modules at_peak at_t50")
    for rt, v in sig.items():
        print(f"  {rt:>16} {v[4]:>4d} {100 * v[2]:>7.3f} {100 * v[3]:>7.3f}")

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(len(res), 2, figsize=(13, 3.6 * len(res)), squeeze=False)
    for j, (name, k) in enumerate((("10-90% rise time (ns)", 0), ("90-10% fall time (ns)", 1))):
        allv = np.concatenate([v[k] for v in res.values()])
        lo, hi = np.percentile(allv, [0.5, 99.5])
        pad = 0.1 * (hi - lo)
        bins = np.linspace(lo - pad, hi + pad, 41)
        for i, (rt, v) in enumerate(res.items()):
            ax = axs[i, j]
            ax.hist(v[k], bins=bins, color=f"C{i}", alpha=0.8)
            med, sg = stats(v[k])
            ax.axvline(med, color="k", lw=1)
            ax.text(0.97, 0.9, f"{rt}  N={v[k].size}\nmed={med:.2f}  $\\sigma_{{MAD}}$={sg:.2f}",
                    transform=ax.transAxes, ha="right", va="top", fontsize=14)
            ax.set_ylabel("modules")
            if i == len(res) - 1:
                ax.set_xlabel(name)
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"timing_dist.{ext}", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(9, 5.5))
    for i, (rt, v) in enumerate(sig.items()):
        ax.plot(v[0], 100 * v[1], color=f"C{i}", label=f"{rt} (N={v[4]})")
    ax.set_xlabel("t after onset (ns)")
    ax.set_ylabel("array-median shape_sigma (% of peak)")
    ax.set_xlim(-8, 150)
    ax.legend()
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(OUT / f"timing_sigma.{ext}", dpi=150)
    plt.close(fig)
    print(f"\nwrote {OUT}/ timing_summary.tsv timing_modules.tsv timing_sigma.tsv "
          "timing_dist.png/.pdf timing_sigma.png/.pdf")
    maxdiff_map(permod, set().union(*excl.values()))
    run_type_maps(permod, excl)
    sigma_maps(sigmod)
    jitter_maps(jit, pj)
    pulse_jitter_maps(pj)


def module_collection(rec, vals, cmap):
    """PatchCollection of module rectangles; rec rows (name, x, y, sx, sy, ...); NaN in vals -> cmap bad colour."""
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle

    pc = PatchCollection([Rectangle((r[1] - r[3] / 2, r[2] - r[4] / 2), r[3], r[4]) for r in rec],
                         cmap=cmap, edgecolor="none")
    pc.set_array(np.ma.masked_invalid(vals))
    return pc


def colour_bar(ax, pc, label):
    cb = ax.figure.colorbar(pc, ax=ax, fraction=0.046, pad=0.02)
    cb.set_label(label, fontsize=16, fontweight="bold")
    cb.ax.tick_params(labelsize=14)
    return cb


def hatch_excluded(ax, rec, names):
    """Grey hatched overlay on the modules of rec whose name is in names (dead)."""
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Rectangle

    ex = [r for r in rec if r[0] in names]
    if ex:
        ax.add_collection(PatchCollection(
            [Rectangle((r[1] - r[3] / 2, r[2] - r[4] / 2), r[3], r[4]) for r in ex],
            facecolor="0.82", edgecolor="0.4", hatch="////", linewidth=0))


def maxdiff_map(permod, excluded):
    """Per-module max-min of rise and fall over run types; heat maps at true positions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    mods = [m for m in json.load(open(HYCAL_MAP)) if m["t"] in ("PbWO4", "PbGlass")]
    FIGDIR.mkdir(parents=True, exist_ok=True)
    rec = []          # (name, x, y, sx, sy, n, rise_maxdiff, fall_maxdiff)
    for m in mods:
        g = m["geo"]
        v = np.array(list(permod.get(m["n"], {}).values())).reshape(-1, 2)
        n = len(v)
        d = (np.ptp(v[:, 0]), np.ptp(v[:, 1])) if n >= 2 else (np.nan, np.nan)
        sd = (np.std(v[:, 0], ddof=1), np.std(v[:, 1], ddof=1)) if n >= 3 else (np.nan, np.nan)
        rec.append((m["n"], g["x"], g["y"], g["sx"], g["sy"], n, *d, *sd))
    with open(OUT / "timing_maxdiff.tsv", "w") as fh:
        fh.write("module\tx\ty\tn_runtypes\trise_maxdiff_ns\tfall_maxdiff_ns\trise_std_ns\tfall_std_ns\n")
        for r in rec:
            fh.write(f"{r[0]}\t{r[1]:.3f}\t{r[2]:.3f}\t{r[5]}\t{r[6]:.3f}\t{r[7]:.3f}\t{r[8]:.3f}\t{r[9]:.3f}\n")

    fig, axs = plt.subplots(1, 2, figsize=(19, 9.5))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.82")
    for ax, (name, k) in zip(axs, (("10-90% rise time", 6), ("90-10% fall time", 7))):
        vals = np.array([r[k] for r in rec])
        pc = module_collection(rec, vals, cmap)
        ok = vals[np.isfinite(vals)]
        pc.set_clim(0, np.percentile(ok, 99) if ok.size else 1)
        ax.add_collection(pc)
        hatch_excluded(ax, rec, excluded)
        ax.set_xlim(-MAP_LIM, MAP_LIM)
        ax.set_ylim(-MAP_LIM, MAP_LIM)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(f"{name}: largest difference between run types\n"
                     f"(grey = fewer than 2 own templates, hatched = dead in any run type)\n"
                     f"run-type average templates not counted; N={ok.size}")
        colour_bar(ax, pc, "max - min (ns)")
    fig.tight_layout()
    for ext in ("png", "pdf"):
        fig.savefig(FIGDIR / f"timing_maxdiff_map.{ext}", dpi=150)
    plt.close(fig)
    n2 = sum(1 for r in rec if r[5] >= 2)
    parts = []
    for nm, k in (("rise", 6), ("fall", 7)):
        v = np.array([r[k] for r in rec if r[5] >= 2])
        if v.size:
            parts.append(f"{nm} maxdiff median {np.median(v):.3f} ns, largest {v.max():.3f} ns, "
                         f"p90 {np.percentile(v, 90):.3f} ns")
    print(f"max-diff over {n2} modules with >=2 run types: " + "; ".join(parts))
    print(f"wrote {FIGDIR}/ timing_maxdiff_map.png/.pdf; {OUT}/timing_maxdiff.tsv")
    std_map(rec, excluded)


def std_map(rec, excluded):
    """Sample standard deviation (ddof = 1) of rise and fall over run types (>= 3 own templates); heat maps."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    n3 = sum(1 for r in rec if r[5] >= 3)
    q = []
    for k in (8, 9):
        v = np.array([r[k] for r in rec if r[5] >= 3])
        q += [np.median(v), np.percentile(v, 90)] if v.size else [np.nan] * 2
    with open(OUT / "timing_summary.tsv", "a") as fh:
        fh.write("\t".join(["across_run_types"] + ["nan"] * 14 + [f"{x:.3f}" for x in q]) + "\n")
    print(f"std over run types, {n3} modules with >=3 run types: rise median {q[0]:.3f} ns, p90 {q[1]:.3f} ns; "
          f"fall median {q[2]:.3f} ns, p90 {q[3]:.3f} ns")

    fig, axs = plt.subplots(1, 2, figsize=(19, 9.5))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.82")
    for ax, (name, k) in zip(axs, (("10-90% rise time", 8), ("90-10% fall time", 9))):
        vals = np.array([r[k] for r in rec])
        pc = module_collection(rec, vals, cmap)
        ok = vals[np.isfinite(vals)]
        pc.set_clim(0, np.percentile(ok, 99) if ok.size else 1)
        ax.add_collection(pc)
        hatch_excluded(ax, rec, excluded)
        ax.set_xlim(-MAP_LIM, MAP_LIM)
        ax.set_ylim(-MAP_LIM, MAP_LIM)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(f"{name}; N={ok.size} modules with at least 3 run types")
        colour_bar(ax, pc, "standard deviation (ns)")
    fig.suptitle("standard deviation of the rise / fall time across run types")
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.legend(handles=[Patch(facecolor="0.82", edgecolor="0.4", label="fewer than 3 run types"),
                        Patch(facecolor="0.82", edgecolor="0.4", hatch="////", label="dead module")],
               loc="lower center", ncol=2, fontsize=14)
    for ext in ("png", "pdf"):
        fig.savefig(FIGDIR / f"timing_std_map.{ext}", dpi=150)
    plt.close(fig)
    print(f"wrote {FIGDIR}/ timing_std_map.png/.pdf")


def run_type_maps(permod, excl):
    """Per run type read: PbWO4 rise (left) and fall (right) at module positions.  Colour limits are the
    1st-99th percentile pooled over those run types, so the maps compare; grey = run-type average template
    (not counted), hatched = dead.  timing_map_all.png/.pdf stacks the run types, one row each, one colour bar
    per column."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    mods = [m for m in json.load(open(HYCAL_MAP)) if m["t"] == "PbWO4"]
    FIGDIR.mkdir(parents=True, exist_ok=True)
    rec = [(m["n"], m["geo"]["x"], m["geo"]["y"], m["geo"]["sx"], m["geo"]["sy"]) for m in mods]
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.82")
    panels = (("10-90% rise time", 0, "rise time (ns)"), ("90-10% fall time", 1, "fall time (ns)"))
    lims = {}
    for k in (0, 1):
        pooled = np.array([pair[k] for mod in permod.values() for pair in mod.values()])
        lims[k] = np.percentile(pooled, [1, 99])
    print("colour limits (1st-99th pooled): " + ", ".join(
        f"{nm} {lims[k][0]:.2f} to {lims[k][1]:.2f} ns" for nm, k, _ in panels))
    def draw(ax, rt, k):
        vals = np.array([permod.get(r[0], {}).get(rt, (np.nan, np.nan))[k] for r in rec], dtype=float)
        pc = module_collection(rec, vals, cmap)
        pc.set_clim(*lims[k])
        ax.add_collection(pc)
        hatch_excluded(ax, rec, excl[rt])
        ax.set_xlim(-MAP_LIM, MAP_LIM)
        ax.set_ylim(-MAP_LIM, MAP_LIM)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        return pc

    for rt in excl:
        n = sum(rt in permod.get(r[0], {}) for r in rec)
        fig, axs = plt.subplots(1, 2, figsize=(19, 9.5))
        for ax, (name, k, label) in zip(axs, panels):
            pc = draw(ax, rt, k)
            ax.set_title(f"{name}  (N={n} of {len(rec)} PbWO4)\n"
                         f"grey = run-type average template, not counted; hatched = dead")
            colour_bar(ax, pc, label)
        fig.suptitle(rt)
        fig.tight_layout()
        for ext in ("png", "pdf"):
            fig.savefig(FIGDIR / f"timing_map_{rt}.{ext}", dpi=150)
        plt.close(fig)
        print(f"[{rt}] {n} of {len(rec)} PbWO4 modules coloured")

    # overview: one row per run type read, rise left and fall right; one colour bar per column
    fig, axs = plt.subplots(len(excl), 2, figsize=(19, 9.5 * len(excl)), squeeze=False, layout="constrained")
    pcs = {}
    for i, rt in enumerate(excl):
        n = sum(rt in permod.get(r[0], {}) for r in rec)
        for ax, (name, k, label) in zip(axs[i], panels):
            pcs[k] = draw(ax, rt, k)
            ax.set_title(f"{rt}: {name}  (N={n} of {len(rec)} PbWO4)")
    for j, (_, k, label) in enumerate(panels):
        cb = fig.colorbar(pcs[k], ax=list(axs[:, j]), fraction=0.046, pad=0.02)
        cb.set_label(label, fontsize=16, fontweight="bold")
        cb.ax.tick_params(labelsize=14)
    legend = [Patch(facecolor="0.82", edgecolor="0.4", hatch="////", label="dead module"),
              Patch(facecolor="0.82", edgecolor="0.4", label="run-type average template, not counted")]
    fig.legend(handles=legend, loc="outside lower center", ncol=2, fontsize=14)
    for ext in ("png", "pdf"):
        fig.savefig(FIGDIR / f"timing_map_all.{ext}", dpi=110)
    plt.close(fig)
    print(f"wrote {FIGDIR}/ timing_map_<run type>.png/.pdf timing_map_all.png/.pdf")


def sigma_maps(sigmod):
    """PbWO4 heat maps of the pulse-to-pulse fluctuation sigma(t) in percent of peak at the pulse peak and at
    the 50% rising-edge point.  Colour scale per quantity: 0 to the 99th percentile over all run types
    (modules with their own sigma); dead modules hatched "////", borrowed sigma (incl. single-pulse) "\\\\"."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Patch, Rectangle

    mods = [m for m in json.load(open(HYCAL_MAP)) if m["t"] == "PbWO4"]
    rec = [(m["n"], m["geo"]["x"], m["geo"]["y"], m["geo"]["sx"], m["geo"]["sy"]) for m in mods]
    nan = (np.nan, np.nan, 0, 0)
    panels = (("at the pulse peak", 0), ("at the 50% rising-edge point", 1))
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.82")
    vmax = {}
    for k in (0, 1):
        own = np.array([v[k] for sm in sigmod.values() for v in sm.values() if not v[2] and not v[3]])
        vmax[k] = np.nanpercentile(own, 99)
    print("shared colour scale vmax (99th percentile, % of peak): " + ", ".join(
        f"{nm} {vmax[k]:.3f}" for nm, k in panels))
    with open(OUT / "sigma_modules.tsv", "w") as fh:
        fh.write("module\tx\ty\t" + "\t".join(f"{rt}_{c}" for rt in sigmod
                 for c in ("sigma_peak_pct", "sigma_t50_pct", "borrowed")) + "\n")
        for r in rec:
            cells = []
            for sm in sigmod.values():
                v = sm.get(r[0], nan)
                cells += [f"{v[0]:.3f}", f"{v[1]:.3f}", str(v[2])]
            fh.write(f"{r[0]}\t{r[1]:.3f}\t{r[2]:.3f}\t" + "\t".join(cells) + "\n")

    def draw(ax, rt, k, label):
        sm = sigmod[rt]
        vals = np.array([sm.get(r[0], nan)[k] if not sm.get(r[0], nan)[3] else np.nan for r in rec])
        pc = module_collection(rec, vals, cmap)
        pc.set_clim(0, vmax[k])
        ax.add_collection(pc)
        for hatch, sel in (("////", lambda v: v[3]), ("\\\\\\\\", lambda v: v[2] and not v[3])):
            ex = [r for r in rec if r[0] in sm and sel(sm[r[0]])]
            if ex:
                ax.add_collection(PatchCollection(
                    [Rectangle((r[1] - r[3] / 2, r[2] - r[4] / 2), r[3], r[4]) for r in ex],
                    facecolor="0.82", edgecolor="0.4", hatch=hatch, linewidth=0))
        ax.set_xlim(-MAP_LIM, MAP_LIM)
        ax.set_ylim(-MAP_LIM, MAP_LIM)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(f"pulse-to-pulse fluctuation sigma(t) {label}")
        colour_bar(ax, pc, f"sigma(t) (% of peak), maximum {vmax[k]:.2f} = 99th percentile")

    legend = [Patch(facecolor="0.82", edgecolor="0.4", hatch="////", label="dead module"),
              Patch(facecolor="0.82", edgecolor="0.4", hatch="\\\\\\\\",
                    label="borrowed sigma(t) (fewer than 50 band pulses)")]
    for rt in sigmod:
        fig, axs = plt.subplots(1, 2, figsize=(19, 9.5))
        for ax, (label, k) in zip(axs, panels):
            draw(ax, rt, k, label)
        fig.suptitle(f"{rt}: pulse-to-pulse fluctuation sigma(t)")
        fig.legend(handles=legend, loc="lower center", ncol=2, fontsize=14)
        fig.tight_layout(rect=(0, 0.04, 1, 1))
        for ext in ("png", "pdf"):
            fig.savefig(FIGDIR / f"sigma_map_{rt}.{ext}", dpi=150)
        plt.close(fig)
    fig, axs = plt.subplots(len(sigmod), 2, figsize=(19, 9 * len(sigmod)), squeeze=False)
    for i, rt in enumerate(sigmod):
        for ax, (label, k) in zip(axs[i], panels):
            draw(ax, rt, k, label)
            ax.set_title(f"{rt}: pulse-to-pulse fluctuation sigma(t) {label}")
    fig.legend(handles=legend, loc="lower center", ncol=2, fontsize=14)
    fig.tight_layout(rect=(0, 0.01, 1, 1))
    fig.savefig(FIGDIR / "sigma_map_all.png", dpi=110)
    plt.close(fig)
    for rt, sm in sigmod.items():
        own = [v for v in sm.values() if not v[2] and not v[3]]
        q = lambda k: (np.nanmedian([v[k] for v in own]), np.nanpercentile([v[k] for v in own], 90))   # noqa: E731
        nb = sum(1 for r in rec if r[0] in sm and sm[r[0]][2] and not sm[r[0]][3])
        nd = sum(1 for r in rec if r[0] in sm and sm[r[0]][3])
        print(f"[{rt}] PbWO4 own-sigma {sum(1 for r in rec if r[0] in sm and not sm[r[0]][2] and not sm[r[0]][3])}: "
              f"peak median {q(0)[0]:.3f} p90 {q(0)[1]:.3f}; t50 median {q(1)[0]:.3f} p90 {q(1)[1]:.3f}; "
              f"hatched borrowed {nb}, dead {nd}")
    print(f"wrote {FIGDIR}/ sigma_map_<run type>.png/.pdf sigma_map_all.png; {OUT}/sigma_modules.tsv")


def jitter_maps(jit, pj):
    """PbWO4 heat maps of the rising-edge time jitter sigma(t) / slope at 10, 50 and 90% of the peak (ns), one
    colour scale per panel shared by all run types (0 to the 99th percentile over own-sigma modules); dead
    hatched "////", borrowed sigma hatched "\\\\"; also jitter_modules.tsv and the median / 90th percentile."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.collections import PatchCollection
    from matplotlib.patches import Patch, Rectangle

    mods = [m for m in json.load(open(HYCAL_MAP)) if m["t"] == "PbWO4"]
    rec = [(m["n"], m["geo"]["x"], m["geo"]["y"], m["geo"]["sx"], m["geo"]["sy"]) for m in mods]
    nan = (np.nan,) * 4 + (0, 0)
    fracs = (10, 50, 90)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.82")
    vmax = [np.nanpercentile([v[k] for jm in jit.values() for n, v in jm.items()
                              if n in {r[0] for r in rec} and not v[4] and not v[5]], 99) for k in range(3)]
    print("jitter shared colour scale vmax (99th percentile, ns): " + ", ".join(
        f"{f}% {vmax[k]:.3f}" for k, f in enumerate(fracs)))
    cols = ("jitter_10_ns", "jitter_50_ns", "jitter_90_ns", "approx_50_ns", "borrowed",
            "rise_jitter_ns", "fall_jitter_ns", "rise_med_ns", "fall_med_ns")
    with open(OUT / "jitter_modules.tsv", "w") as fh:
        fh.write("module\tx\ty\t" + "\t".join(f"{rt}_{c}" for rt in jit for c in cols) + "\n")
        for r in rec:
            cells = []
            for rt, jm in jit.items():
                v = jm.get(r[0], nan)
                cells += [f"{x:.4f}" for x in v[:4]] + [str(v[4])] + [f"{x:.4f}" for x in pj[rt].get(r[0], (np.nan,) * 5)[:4]]
            fh.write(f"{r[0]}\t{r[1]:.3f}\t{r[2]:.3f}\t" + "\t".join(cells) + "\n")

    def draw(ax, rt, k):
        jm = jit[rt]
        vals = np.array([jm.get(r[0], nan)[k] if not jm.get(r[0], nan)[5] else np.nan for r in rec])
        pc = module_collection(rec, vals, cmap)
        pc.set_clim(0, vmax[k])
        ax.add_collection(pc)
        for hatch, sel in (("////", lambda v: v[5]), ("\\\\\\\\", lambda v: v[4] and not v[5])):
            ex = [r for r in rec if r[0] in jm and sel(jm[r[0]])]
            if ex:
                ax.add_collection(PatchCollection(
                    [Rectangle((r[1] - r[3] / 2, r[2] - r[4] / 2), r[3], r[4]) for r in ex],
                    facecolor="0.82", edgecolor="0.4", hatch=hatch, linewidth=0))
        ax.set_xlim(-MAP_LIM, MAP_LIM)
        ax.set_ylim(-MAP_LIM, MAP_LIM)
        ax.set_aspect("equal")
        ax.set_xlabel("x (mm)")
        ax.set_ylabel("y (mm)")
        ax.set_title(f"{rt}: time jitter of the rising edge, sigma(t) / slope\nat {fracs[k]}% of the peak")
        colour_bar(ax, pc, f"time jitter (ns), maximum {vmax[k]:.3f} = 99th percentile")

    legend = [Patch(facecolor="0.82", edgecolor="0.4", hatch="////", label="dead module"),
              Patch(facecolor="0.82", edgecolor="0.4", hatch="\\\\\\\\",
                    label="borrowed sigma(t) (fewer than 50 band pulses)")]
    for rt in jit:
        fig, axs = plt.subplots(1, 3, figsize=(28, 9.5))
        for k, ax in enumerate(axs):
            draw(ax, rt, k)
        fig.suptitle(f"{rt}: time jitter of the rising edge, sigma(t) / slope")
        fig.legend(handles=legend, loc="lower center", ncol=2, fontsize=14)
        fig.tight_layout(rect=(0, 0.04, 1, 0.97))
        for ext in ("png", "pdf"):
            fig.savefig(FIGDIR / f"jitter_map_{rt}.{ext}", dpi=150)
        plt.close(fig)
    fig, axs = plt.subplots(len(jit), 3, figsize=(28, 9 * len(jit)), squeeze=False)
    for i, rt in enumerate(jit):
        for k, ax in enumerate(axs[i]):
            draw(ax, rt, k)
    fig.legend(handles=legend, loc="lower center", ncol=2, fontsize=14)
    fig.tight_layout(rect=(0, 0.01, 1, 1))
    fig.savefig(FIGDIR / "jitter_map_all.png", dpi=110)
    plt.close(fig)
    for rt, jm in jit.items():
        own = np.array([v[:4] for v in jm.values() if not v[4] and not v[5]])
        med, p90 = np.nanmedian(own, axis=0), np.nanpercentile(own, 90, axis=0)
        ratio = np.nanmedian(own[:, 1] / own[:, 3])
        print(f"[{rt}] jitter own-sigma {len(own)}: " + "; ".join(
            f"{f}% median {med[k]:.3f} p90 {p90[k]:.3f}" for k, f in enumerate(fracs))
            + f"; median exact/approx at 50% {ratio:.3f}")
    print(f"wrote {FIGDIR}/ jitter_map_<run type>.png/.pdf jitter_map_all.png; {OUT}/jitter_modules.tsv")


def pulse_jitter_maps(pj):
    """PbWO4 maps of rise_jitter_ns (10-90 % rise time) and fall_jitter_ns (90-10 % fall time) per run type
    (rise_jitter_map_<run type>, fall_jitter_map_<run type>, + _all.png); one colour scale per quantity, 0 to the
    99th percentile over all run types; grey = fewer than 50 band pulses, hatched "////" = dead."""
    import matplotlib.pyplot as plt
    from matplotlib.patches import Patch

    mods = [m for m in json.load(open(HYCAL_MAP)) if m["t"] == "PbWO4"]
    rec = [(m["n"], m["geo"]["x"], m["geo"]["y"], m["geo"]["sx"], m["geo"]["sy"]) for m in mods]
    nan = (np.nan,) * 4 + (0,)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("0.82")
    legend = [Patch(facecolor="0.82", edgecolor="0.4", hatch="////", label="dead module"),
              Patch(facecolor="0.82", edgecolor="0.4", label="fewer than 50 band pulses")]
    for key, k, what in (("rise", 0, "10\u201390 % rise time"), ("fall", 1, "90\u201310 % fall time")):
        vmax = np.nanpercentile([v[k] for pm in pj.values() for v in pm.values() if not v[4]], 99)

        def draw(ax, rt):
            pm = pj[rt]
            vals = np.array([pm.get(r[0], nan)[k] if not pm.get(r[0], nan)[4] else np.nan for r in rec])
            pc = module_collection(rec, vals, cmap)
            pc.set_clim(0, vmax)
            ax.add_collection(pc)
            hatch_excluded(ax, rec, {n for n, v in pm.items() if v[4]})
            ax.set_xlim(-MAP_LIM, MAP_LIM)
            ax.set_ylim(-MAP_LIM, MAP_LIM)
            ax.set_aspect("equal")
            ax.set_xlabel("x (mm)")
            ax.set_ylabel("y (mm)")
            ax.set_title(f"{rt}: pulse-to-pulse jitter of the {what}", fontsize=16)
            colour_bar(ax, pc, f"jitter (ns), maximum {vmax:.2f} = 99th percentile")

        for rt in pj:
            fig, ax = plt.subplots(figsize=(11, 9.5))
            draw(ax, rt)
            fig.legend(handles=legend, loc="lower center", ncol=2, fontsize=14)
            fig.tight_layout(rect=(0, 0.05, 1, 1))
            for ext in ("png", "pdf"):
                fig.savefig(FIGDIR / f"{key}_jitter_map_{rt}.{ext}", dpi=150)
            plt.close(fig)
        fig, axs = plt.subplots(len(pj), 1, figsize=(11, 9.5 * len(pj)), squeeze=False)
        for i, rt in enumerate(pj):
            draw(axs[i, 0], rt)
        fig.legend(handles=legend, loc="lower center", ncol=2, fontsize=14)
        fig.tight_layout(rect=(0, 0.01, 1, 1))
        fig.savefig(FIGDIR / f"{key}_jitter_map_all.png", dpi=110)
        plt.close(fig)
        for rt, pm in pj.items():
            x = np.array([v[k] for v in pm.values() if not v[4]])
            x = x[np.isfinite(x)]
            print(f"[{rt}] {key}_jitter: {x.size} modules, median {np.median(x):.3f} p90 {np.percentile(x, 90):.3f} ns")
    print(f"wrote {FIGDIR}/ rise_jitter_map_<run type> fall_jitter_map_<run type> .png/.pdf and _all.png")


if __name__ == "__main__":
    main()
