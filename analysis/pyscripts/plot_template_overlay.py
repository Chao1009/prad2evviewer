#!/usr/bin/env python3
"""plot_template_overlay.py -- per-module overlay of the spline template banks.

Reads output/templates_spline/<rt>/templates_<rt>.npz (spline_templates_v1) for
the run types 0p73GeV_H2, 2p24GeV_H2, 3p49GeV_H2 and 2p24GeV_Ta_X17 (read-only)
and writes one PDF:

  page 1      cover: tier fractions per run type, module-to-module spread of
              I_int (iint_wa) and of the template FWHM, median template per run
              type with its difference to the 2.24 GeV H2 median, summary table
  pages 2...  one panel per PbWO4 module (4 per page, ordered by module id):
              the unit-peak templates vs t - onset (ns), legend with N (n_pool),
              tier (fb0/fb1/fb2), t_pk, t50, I_int; fb2 (run-type average) rows
              dashed; lower panel = template - reference template in % of peak

Templates are evaluated as mf_pulse.SampledTemplateBank.from_npz does: the
row's cubic B-spline held at its last value beyond the knot range, 0 for
t <= 0; a row with valid == 0 takes the average spline and counts as fb2.
Reference per module: 2.24 GeV H2, or the best-statistics (largest n_pool)
non-fb2 template when the 2.24 GeV H2 row is fb2.

Region summary (printed and on the cover): median over modules of the median
|delta| (% of peak) in rise 0 < t < t_pk-3, peak |t - t_pk| <= 3 and tail
t_pk+20 <= t <= 200 ns (t_pk of the reference), pairs with an fb2 side excluded.

    MPLBACKEND=Agg python plot_template_overlay.py [--limit 8]
"""

import argparse
import os
import sys

import numpy as np

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages  # noqa: E402
from scipy.interpolate import BSpline  # noqa: E402

# (slug, label, colour) -- fixed order and colour per run type
RUNS = [("0p73GeV_H2", "0.73 GeV H2", "#2a78d6"),
        ("2p24GeV_H2", "2.24 GeV H2", "#eb6834"),
        ("3p49GeV_H2", "3.49 GeV H2", "#1baf7a"),
        ("2p24GeV_Ta_X17", "X17 Ta pooled", "#eda100")]
REF = "2p24GeV_H2"
T_PLOT = np.arange(-10.0, 200.0 + 1e-9, 0.25)
T_FINE = np.arange(0.0, 200.0 + 1e-9, 0.05)
TIER_COL = ["#184f95", "#5598e7", "#b7d3f6"]       # fb0 / fb1 / fb2


def eval_rows(z, t):
    """[n_rows, len(t)] unit-peak templates (as from_npz builds them) and tiers."""
    k = int(z["spline_k"])

    def ev(kt, c):
        s = BSpline(np.asarray(kt, np.float64), np.asarray(c, np.float64), k)
        tmax = float(s.t[-k - 1])
        return np.where(t > 0, s(np.minimum(t, tmax)), 0.0)

    at, ac = z["avg_spline_t"], z["avg_spline_c"]
    avg = ev(at[np.isfinite(at)], ac[np.isfinite(ac)])
    nt, nc, valid = z["spline_nt"], z["spline_nc"], z["valid"]
    out = np.empty((len(nt), len(t)))
    tier = np.asarray(z["fallback"], dtype=int).copy()
    for r in range(len(nt)):
        if valid[r] > 0:
            out[r] = ev(z["spline_t"][r, :nt[r]], z["spline_c"][r, :nc[r]])
        else:
            out[r], tier[r] = avg, 2
    return out, tier


def fwhm(y, t):
    """FWHM (ns) of each row of y on grid t (linear interpolation at 0.5)."""
    res = np.full(len(y), np.nan)
    for i, v in enumerate(y):
        ip = int(np.argmax(v))
        if v[ip] <= 0.5:
            continue
        up = np.nonzero(v[:ip + 1] < 0.5)[0]
        dn = np.nonzero(v[ip:] < 0.5)[0]
        if len(up) == 0 or len(dn) == 0:
            continue
        a, b = up[-1], ip + dn[0]
        t1 = t[a] + (0.5 - v[a]) / (v[a + 1] - v[a]) * (t[a + 1] - t[a])
        t2 = t[b - 1] + (v[b - 1] - 0.5) / (v[b - 1] - v[b]) * (t[b] - t[b - 1])
        res[i] = t2 - t1
    return res


def load(base):
    banks = {}
    for slug, _, _ in RUNS:
        path = os.path.join(base, slug, f"templates_{slug}.npz")
        z = np.load(path, allow_pickle=False)
        names = [str(x) for x in z["names"]]
        y, tier = eval_rows(z, T_PLOT)
        yf, _ = eval_rows(z, T_FINE)
        banks[slug] = dict(
            row={nm: i for i, nm in enumerate(names)}, names=names,
            mid=np.asarray(z["module_id"]), y=y, tier=tier,
            n=np.asarray(z["n_pool"]), t_pk=np.asarray(z["t_pk"]),
            t50=np.asarray(z["t50"]), iint=np.asarray(z["iint_wa"]),
            fwhm=fwhm(yf, T_FINE))
        print(f"{slug}: {len(names)} rows, tiers "
              + " ".join(f"fb{k}={int((tier == k).sum())}" for k in range(3)))
    return banks


def module_list(banks):
    ids = {}
    for b in banks.values():
        for nm, m in zip(b["names"], b["mid"]):
            if nm.startswith("W"):                     # PbWO4 only
                ids[nm] = int(m)
    return sorted(ids, key=lambda nm: ids[nm])


def pick_ref(banks, nm):
    b = banks[REF]
    if nm in b["row"] and b["tier"][b["row"][nm]] < 2:
        return REF
    best, nbest = REF, -1
    for slug, _, _ in RUNS:
        bb = banks[slug]
        r = bb["row"].get(nm)
        if r is not None and bb["tier"][r] < 2 and bb["n"][r] > nbest:
            best, nbest = slug, bb["n"][r]
    return best


def region_masks(t_pk):
    return {"rise": (T_PLOT > 0) & (T_PLOT < t_pk - 3),
            "peak": np.abs(T_PLOT - t_pk) <= 3,
            "tail": (T_PLOT >= t_pk + 20) & (T_PLOT <= 200)}


def module_stats(banks, mods):
    """Per module: reference, {slug: {region: med|d|, 'max': max|d|}} (non-fb2 pairs)."""
    stats = {}
    for nm in mods:
        ref = pick_ref(banks, nm)
        br = banks[ref]
        rr = br["row"].get(nm)
        st = {"ref": ref, "d": {}}
        if rr is None or br["tier"][rr] == 2:
            stats[nm] = st
            continue
        msk = region_masks(br["t_pk"][rr])
        for slug, _, _ in RUNS:
            b = banks[slug]
            r = b["row"].get(nm)
            if slug == ref or r is None or b["tier"][r] == 2:
                continue
            d = 100.0 * np.abs(b["y"][r] - br["y"][rr])
            st["d"][slug] = {k: float(np.median(d[m])) for k, m in msk.items()}
            st["d"][slug]["max"] = float(d[(T_PLOT > 0)].max())
        stats[nm] = st
    return stats


def region_summary(stats):
    rows = {}
    for slug, _, _ in RUNS + [("all", "", "")]:
        v = [s["d"][sl] for s in stats.values() for sl in s["d"]
             if slug in ("all", sl)]
        if v:
            rows[slug] = {k: float(np.median([x[k] for x in v]))
                          for k in ("rise", "peak", "tail", "max")} | {"n": len(v)}
    return rows


def cover(pdf, banks, mods, summ):
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(3, 3, height_ratios=[1, 1, 0.9], hspace=0.55, wspace=0.3)
    fig.suptitle("HyCal PbWO4 spline templates: run-type comparison "
                 f"({len(mods)} modules)", fontsize=12)
    labels = [lab for _, lab, _ in RUNS]

    ax = fig.add_subplot(gs[0, 0])
    left = np.zeros(len(RUNS))
    for k in range(3):
        fr = np.array([np.mean(banks[s]["tier"][[banks[s]["row"][m] for m in mods
                                                   if m in banks[s]["row"]]] == k)
                       for s, _, _ in RUNS])
        ax.barh(labels, fr, left=left, color=TIER_COL[k], edgecolor="white",
                linewidth=1, label=f"fb{k}")
        for i, (l0, f) in enumerate(zip(left, fr)):
            if f > 0.06:
                ax.text(l0 + f / 2, i, f"{100 * f:.0f}%", ha="center", va="center",
                        fontsize=7, color="white" if k < 2 else "black")
        left += fr
    ax.invert_yaxis()
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of modules")
    ax.set_title("fallback tier", fontsize=9)
    ax.legend(fontsize=7, ncol=3, loc="lower center", bbox_to_anchor=(0.5, 1.12),
              frameon=False)
    ax.tick_params(labelsize=7)

    for col, key, xl in ((1, "iint", "I_int (iint_wa)"), (2, "fwhm", "FWHM (ns)")):
        ax = fig.add_subplot(gs[0, col])
        allv = np.concatenate([banks[s][key][banks[s]["tier"] < 2] for s, _, _ in RUNS])
        lo, hi = np.nanpercentile(allv, [0.5, 99.5])
        bins = np.linspace(lo, hi, 50)
        for s, lab, c in RUNS:
            v = banks[s][key][banks[s]["tier"] < 2]
            med = np.nanmedian(v)
            rs = 100 * 0.7413 * np.subtract(*np.nanpercentile(v, [75, 25])) / med
            ax.hist(v, bins=bins, histtype="step", color=c, lw=1.5,
                    label=f"{lab}: med {med:.3g}, spread {rs:.1f}%")
        ax.set_xlabel(xl)
        ax.set_ylabel("modules (fb0+fb1)", fontsize=8)
        ax.set_ylim(0, ax.get_ylim()[1] * 1.6)
        ax.set_title(f"module-to-module spread of {xl.split()[0]}", fontsize=9)
        ax.legend(fontsize=5.5, frameon=False, loc="upper left")
        ax.tick_params(labelsize=7)

    sub = gs[1, :].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
    ax, axd = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
    meds = {}
    for s, lab, c in RUNS:
        b = banks[s]
        meds[s] = np.median(b["y"][b["tier"] < 2], axis=0)
        ax.plot(T_PLOT, meds[s], color=c, lw=1.5, label=lab)
    for s, lab, c in RUNS:
        if s != REF:
            axd.plot(T_PLOT, 100 * (meds[s] - meds[REF]), color=c, lw=1.2)
    axd.axhline(0, color="0.5", lw=0.8)
    ax.set_title("median template per run type (fb0+fb1 modules)", fontsize=9)
    ax.set_ylabel("unit-peak shape", fontsize=7)
    ax.legend(fontsize=7, frameon=False)
    ax.tick_params(labelbottom=False, labelsize=7)
    axd.set_ylabel("Δ vs 2.24 H2\n(% of peak)", fontsize=7)
    axd.set_xlabel("t - onset (ns)")
    axd.tick_params(labelsize=7)
    for a in (ax, axd):
        a.set_xlim(T_PLOT[0], T_PLOT[-1])
        a.grid(alpha=0.3)

    ax = fig.add_subplot(gs[2, :])
    ax.axis("off")
    hdr = (f"{'run type':<16}{'fb0':>6}{'fb1':>6}{'fb2':>6}{'med N':>8}"
           f"{'I_int med':>11}{'FWHM med':>10}   |{'rise':>7}{'peak':>7}{'tail':>7}"
           f"{'max':>7}{'pairs':>7}")
    lines = [hdr, "-" * len(hdr)]
    for s, lab, _ in RUNS:
        b = banks[s]
        sel = [b["row"][m] for m in mods if m in b["row"]]
        tr = b["tier"][sel]
        g = b["tier"] < 2
        sm = summ.get(s)
        dd = (f"{sm['rise']:>7.2f}{sm['peak']:>7.2f}{sm['tail']:>7.2f}{sm['max']:>7.2f}"
              f"{sm['n']:>7d}" if sm else f"{'(reference)':>35}")
        lines.append(f"{lab:<16}" + "".join(f"{np.mean(tr == k):>6.1%}" for k in range(3))
                     + f"{np.median(b['n'][sel]):>8.0f}{np.nanmedian(b['iint'][g]):>11.4f}"
                     f"{np.nanmedian(b['fwhm'][g]):>10.2f}   |" + dd)
    sm = summ.get("all")
    if sm:
        lines.append(f"{'all pairs':<16}{'':>47}   |{sm['rise']:>7.2f}{sm['peak']:>7.2f}"
                     f"{sm['tail']:>7.2f}{sm['max']:>7.2f}{sm['n']:>7d}")
    lines.append("")
    lines.append("right block: median over modules of median |template - reference| "
                 "(% of peak) per region; rise 0<t<t_pk-3, peak |t-t_pk|<=3,")
    lines.append("tail t_pk+20..200 ns; max = median of per-module max |delta|; "
                 "reference 2.24 GeV H2 (best-N non-fb2 if that is fb2); fb2 pairs excluded.")
    lines.append("I_int / FWHM spread = 0.7413 * IQR / median over fb0+fb1 modules.")
    ax.text(0, 1, "\n".join(lines), family="monospace", fontsize=7, va="top")
    pdf.savefig(fig)
    plt.close(fig)


def module_page(pdf, banks, mods, stats, flag_pct):
    fig = plt.figure(figsize=(11, 8.5))
    gs = fig.add_gridspec(2, 2, hspace=0.3, wspace=0.2, left=0.06, right=0.98,
                          top=0.95, bottom=0.06)
    for j, nm in enumerate(mods):
        sub = gs[j // 2, j % 2].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
        ax, axd = fig.add_subplot(sub[0]), fig.add_subplot(sub[1])
        st = stats[nm]
        ref = st["ref"]
        br = banks[ref]
        rr = br["row"].get(nm)
        mid, dmax = None, 0.0
        for s, lab, c in RUNS:
            b = banks[s]
            r = b["row"].get(nm)
            if r is None:
                continue
            mid = int(b["mid"][r])
            fb = int(b["tier"][r])
            ls = "--" if fb == 2 else "-"
            ax.plot(T_PLOT, b["y"][r], color=c, lw=1.2, ls=ls,
                    label=(f"{lab:<14} N={b['n'][r]:<5d} fb{fb}  t_pk {b['t_pk'][r]:5.2f}"
                           f"  t50 {b['t50'][r]:5.2f}  I {b['iint'][r]:.4f}"))
            if s != ref and rr is not None:
                d = 100 * (b["y"][r] - br["y"][rr])
                dmax = max(dmax, float(np.abs(d).max()))
                axd.plot(T_PLOT, d, color=c, lw=1.0, ls=ls)
        mx = max((d["max"] for d in st["d"].values()), default=None)
        flag = mx is not None and mx > flag_pct
        mtxt = "n/a (no non-fb2 pair)" if mx is None else f"{mx:.1f}%"
        ax.set_title(f"{nm} (id {mid})   max|Δ| {mtxt}"
                     + ("  ** run types disagree **" if flag else ""),
                     fontsize=9, color="#b00020" if flag else "black")
        ax.legend(fontsize=5.5, frameon=False, loc="upper right",
                  prop={"family": "monospace", "size": 5.5})
        ax.set_ylim(-0.05, 1.08)
        ax.tick_params(labelbottom=False, labelsize=7)
        ax.set_ylabel("unit-peak shape", fontsize=7)
        axd.axhline(0, color="0.5", lw=0.8)
        lim = min(max(5.0, 1.15 * dmax), 40.0)
        axd.set_ylim(-lim, lim)
        axd.set_ylabel(f"Δ vs {dict((s, l) for s, l, _ in RUNS)[ref]}\n(% of peak)",
                       fontsize=6)
        axd.set_xlabel("t - onset (ns)", fontsize=7)
        axd.tick_params(labelsize=7)
        for a in (ax, axd):
            a.set_xlim(T_PLOT[0], T_PLOT[-1])
            a.grid(alpha=0.3)
    pdf.savefig(fig)
    plt.close(fig)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--base", default="output/templates_spline")
    ap.add_argument("--out", default="output/templates_spline/overlay/template_overlay.pdf")
    ap.add_argument("--limit", type=int, default=0, help="first N modules only (test)")
    ap.add_argument("--flag-pct", type=float, default=10.0,
                    help="flag modules whose max |delta| (non-fb2 pairs) exceeds this %%")
    a = ap.parse_args()

    banks = load(a.base)
    mods = module_list(banks)
    if a.limit:
        mods = mods[:a.limit]
    stats = module_stats(banks, mods)
    summ = region_summary(stats)

    print("median |delta| (% of peak) vs reference, fb2 pairs excluded:")
    for s, v in summ.items():
        print(f"  {s:<16} rise {v['rise']:.2f}  peak {v['peak']:.2f}  tail {v['tail']:.2f}"
              f"  max {v['max']:.2f}  (pairs {v['n']})")
    flagged = sorted(((max((d["max"] for d in s["d"].values()), default=0.0), nm)
                      for nm, s in stats.items()), reverse=True)
    bad = [(m, nm) for m, nm in flagged if m > a.flag_pct]
    print(f"modules with max |delta| > {a.flag_pct}%: {len(bad)}")
    for m, nm in bad[:30]:
        worst = max(stats[nm]["d"].items(), key=lambda kv: kv[1]["max"])
        print(f"  {nm:<6} max {m:6.1f}%  worst {worst[0]} vs {stats[nm]['ref']}  "
              f"(rise {worst[1]['rise']:.1f} peak {worst[1]['peak']:.1f} "
              f"tail {worst[1]['tail']:.1f})")

    os.makedirs(os.path.dirname(os.path.abspath(a.out)), exist_ok=True)
    with PdfPages(a.out) as pdf:
        cover(pdf, banks, mods, summ)
        for i in range(0, len(mods), 4):
            module_page(pdf, banks, mods[i:i + 4], stats, a.flag_pct)
    npage = 1 + -(-len(mods) // 4)
    print(f"wrote {a.out}: {len(mods)} modules, {npage} pages")
    return 0


if __name__ == "__main__":
    sys.exit(main())
