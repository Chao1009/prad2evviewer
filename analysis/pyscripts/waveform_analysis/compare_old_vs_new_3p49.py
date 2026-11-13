#!/usr/bin/env python3
"""Compare old (24386/24503/24512) and new (24698/24911/24937) 3p49GeV_H2 select files and template banks.
Writes pulse_templates/3p49GeV_H2_old/compare_old_vs_new.pdf and compare.txt."""
import glob, json, sys, os
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
# plot style: large bold axis labels, large legends
plt.rcParams.update({
    "axes.labelsize": 16, "axes.labelweight": "bold", "axes.titlesize": 16, "axes.titleweight": "bold",
    "xtick.labelsize": 14, "ytick.labelsize": 14, "legend.fontsize": 14, "axes.linewidth": 1.8,
    "xtick.major.width": 1.5, "ytick.major.width": 1.5, "xtick.major.size": 6, "ytick.major.size": 6,
    "lines.linewidth": 2.0, "figure.titlesize": 18, "figure.titleweight": "bold"})
from matplotlib.collections import PatchCollection
from matplotlib.patches import Rectangle
from matplotlib.backends.backend_pdf import PdfPages
from multiprocessing import Pool

HERE = Path(__file__).resolve().parent
PT = HERE / "pulse_templates"
NEW, OLD = PT / "3p49GeV_H2", PT / "3p49GeV_H2_old"
DB = Path(os.environ.get("PRAD2_DATABASE_DIR") or HERE.parents[2] / "database")
GEO = {m["n"]: m["geo"] for m in json.load(open(DB / "hycal_map.json"))}
TG = None
OFF = 576   # select midx of PbWO4 module i is i + 576 (checked against names)
OUT = []


def say(s=""):
    print(s, flush=True); OUT.append(s)


def per_file(f):
    d = np.load(f, allow_pickle=True)
    meta = json.loads(str(d["meta"]))
    midx, amp, er, rms = d["midx"], d["amp"], d["e_rough"], d["rms"]
    run = int(d["run"])
    nm = d["cutflow"].shape[0]; midx = midx - OFF
    keep = midx >= 0
    midx, amp, er, rms = midx[keep], amp[keep], er[keep], rms[keep]
    n = np.bincount(midx, minlength=1152)
    return dict(run=run, meta=meta, midx=midx, amp=amp, er=er, rms=rms, n=n, cutflow=d["cutflow"][OFF:OFF + 1152],
                names=d["cutflow_names"], pname=d["name"][:1])


def load_set(files):
    with Pool(min(8, len(files))) as p:
        parts = p.map(per_file, sorted(files))
    S = dict(runs=[x["run"] for x in parts], meta=[x["meta"] for x in parts])
    for k in ("midx", "amp", "er", "rms"):
        S[k] = np.concatenate([x[k] for x in parts])
    S["n"] = sum(x["n"] for x in parts)
    S["cutflow"] = sum(x["cutflow"] for x in parts)
    S["cf_names"] = [str(s) for s in parts[0]["names"]]
    return S


def peak_pos(er_rel, midx, nm, lo=0.5, hi=1.2, bw=0.01, nmin=60):
    edges = np.arange(lo, hi + bw / 2, bw)
    out = np.full(nm, np.nan); cnt = np.zeros(nm, int)
    sel = er_rel >= lo
    order = np.argsort(midx[sel], kind="stable")
    m, e = midx[sel][order], er_rel[sel][order]
    bnd = np.searchsorted(m, np.arange(nm + 1))
    for i in range(nm):
        x = e[bnd[i]:bnd[i + 1]]; cnt[i] = x.size
        if x.size < nmin: continue
        h, _ = np.histogram(x, edges)
        hs = np.convolve(h, np.ones(3) / 3, "same")
        j = int(np.argmax(hs))
        c = 0.5 * (edges[j] + edges[j + 1])
        w = x[np.abs(x - c) < 0.06]          # refine: mean within +-0.06 E_beam of the mode
        out[i] = w.mean() if w.size >= 20 else c
    return out, cnt


def rise_fall(f):
    f = np.nan_to_num(np.asarray(f, float)); ip = int(np.argmax(f))
    def cross(level, step):
        j = ip
        while 0 <= j + step < f.size and f[j] >= level: j += step
        if f[j] >= level: return np.nan
        k = j - step
        return float(TG[j] + (level - f[j]) / (f[k] - f[j]) * (TG[k] - TG[j]))
    pk = f[ip]
    return cross(.9 * pk, -1) - cross(.1 * pk, -1), cross(.1 * pk, 1) - cross(.9 * pk, 1)


def mapfig(ax, names, vals, title, cmap="viridis", vmin=None, vmax=None, label="", lim=600):
    rect = [Rectangle((GEO[n]["x"] - GEO[n]["sx"] / 2, GEO[n]["y"] - GEO[n]["sy"] / 2), GEO[n]["sx"], GEO[n]["sy"]) for n in names]
    pc = PatchCollection(rect, cmap=cmap, edgecolor="0.7", lw=0.1)
    v = np.ma.masked_invalid(vals); pc.set_array(v); pc.set_clim(vmin, vmax)
    ax.add_collection(pc)
    ax.set(xlim=(-lim, lim), ylim=(-lim, lim), aspect="equal", title=title, xlabel="x (mm)", ylabel="y (mm)")
    cb = plt.colorbar(pc, ax=ax, fraction=0.046)
    cb.set_label(label, fontsize=16, fontweight="bold"); cb.ax.tick_params(labelsize=14)


def main():
    global TG
    tn = np.load(NEW / "templates_3p49GeV_H2.npz", allow_pickle=True)
    to = np.load(OLD / "templates_3p49GeV_H2.npz", allow_pickle=True)
    names = np.array([str(x) for x in tn["names"]])
    assert (names == np.array([str(x) for x in to["names"]])).all()
    TG = tn["t_grid"]; nm = names.size
    pwo = np.array([n.startswith("W") for n in names])
    New = load_set(glob.glob(str(NEW / "pulses/run_*.npz")))
    Old = load_set(glob.glob(str(OLD / "pulses/run_*.npz")))
    # consistent module ordering check
    assert New["n"].size == nm == Old["n"].size
    mn, mo = New["meta"][0], Old["meta"][0]
    say("# 3p49GeV_H2 old vs new")
    for tag, S in (("old", Old), ("new", New)):
        for r, m in zip(S["runs"], S["meta"]):
            say(f"{tag} run {r}: E_beam {m['e_beam_mev']:.1f} MeV, calib {Path(m['calib_file']).name}, "
                f"trigger {m['trigger']}, dead {','.join(m['dead_modules'])}")
    Eb = {"old": Old["meta"][0]["e_beam_mev"], "new": New["meta"][0]["e_beam_mev"]}
    # per-run E_beam / calib consistency
    for tag, S in (("old", Old), ("new", New)):
        say(f"{tag}: distinct E_beam {sorted(set(round(m['e_beam_mev'],1) for m in S['meta']))}, calib {sorted(set(Path(m['calib_file']).name for m in S['meta']))}")
    nO, nN = Old["n"].astype(float), New["n"].astype(float)
    fO, fN = nO / nO[pwo].sum(), nN / nN[pwo].sum()
    say(f"pulses stored (PbWO4): old {int(nO[pwo].sum())}, new {int(nN[pwo].sum())}")
    ok = pwo & (nO > 20) & (nN > 20)
    ratio = np.where(ok, fO / np.where(fN > 0, fN, np.nan), np.nan)
    rel = {t: S["er"] / Eb[t] for t, S in (("old", Old), ("new", New))}
    pO, cO = peak_pos(rel["old"], Old["midx"], nm)
    pN, cN = peak_pos(rel["new"], New["midx"], nm)
    okp = np.isfinite(pO) & np.isfinite(pN) & pwo
    pr = np.where(okp, pO / pN, np.nan)
    # inner modules: largest e-p pool in both
    score = np.where(okp, np.minimum(cO, cN), -1)
    inner = np.argsort(-score)[:3]
    say(f"3 'inner' modules (largest min(old,new) count above 0.5 E_beam): {', '.join(names[inner])}")
    rr = np.hypot(*(np.array([[GEO[n]["x"], GEO[n]["y"]] for n in names]).T))
    say(f"e-p peak (e_rough/E_beam mode): modules with both {int(okp.sum())}; ratio old/new median {np.nanmedian(pr):.4f}, "
        f"mean {np.nanmean(pr):.4f}, 16-84% {np.nanpercentile(pr,16):.4f}-{np.nanpercentile(pr,84):.4f}")
    for lo, hi in ((0, 150), (150, 250), (250, 400), (400, 1e9)):
        s = okp & (rr >= lo) & (rr < hi)
        if s.sum(): say(f"  r {lo}-{hi if hi < 1e8 else 'inf'} mm ({int(s.sum())} mod): peak ratio median {np.median(pr[s]):.4f}, "
                        f"pulse-fraction ratio median {np.nanmedian(ratio[s]):.3f}")
    hs = okp & (np.minimum(cO, cN) >= 2000)
    say(f"  high-stat modules (>=2000 pulses above 0.5 E_beam in both, {int(hs.sum())}): peak ratio median {np.median(pr[hs]):.4f}, "
        f"16-84% {np.percentile(pr[hs],16):.4f}-{np.percentile(pr[hs],84):.4f}; " + ", ".join(f"{names[i]} {pO[i]:.4f}/{pN[i]:.4f}" for i in inner))
    say(f"pulse-fraction ratio old/new per module: median {np.nanmedian(ratio):.3f}, 16-84% {np.nanpercentile(ratio,16):.3f}-{np.nanpercentile(ratio,84):.3f}")
    # global peak, pooled
    def pooled_peak(S, t):
        x = (S["er"] / Eb[t]); x = x[(x > .5) & (x < 1.2)]
        h, e = np.histogram(x, np.arange(.5, 1.2001, .01)); h = np.convolve(h, np.ones(3) / 3, "same")
        return 0.5 * (e[np.argmax(h)] + e[np.argmax(h) + 1])
    say(f"pooled e-p peak e_rough/E_beam: old {pooled_peak(Old,'old'):.3f}, new {pooled_peak(New,'new'):.3f}")
    # pedestal rms
    for t, S in (("old", Old), ("new", New)):
        r = S["rms"]
        say(f"ped rms {t}: median {np.median(r):.3f}, 16/84% {np.percentile(r,16):.3f}/{np.percentile(r,84):.3f}")
    # cut-flow
    cfO, cfN = Old["cutflow"][pwo].sum(0), New["cutflow"][pwo].sum(0)
    say("\ncut flow (PbWO4 modules), count and fraction of 'nseen' (old | new) and of previous cut:")
    for i, nme in enumerate(Old["cf_names"]):
        pv = max(i - 1, 0)
        say(f"  {nme:16s} {cfO[i]:>12d} {cfO[i]/cfO[0]:.5f} | {cfN[i]:>12d} {cfN[i]/cfN[0]:.5f}   step {cfO[i]/cfO[pv]:.4f} | {cfN[i]/cfN[pv]:.4f}")

    # --- template level
    def tmpl(t):
        valid = (t["valid"] == 1) & (t["dead"] == 0) & pwo
        own = valid & (t["fallback"] != 2) & (t["single"] == 0)
        rf = np.array([rise_fall(t["shape"][i]) if own[i] else (np.nan, np.nan) for i in range(nm)])
        return own, rf
    ownO, rfO = tmpl(to); ownN, rfN = tmpl(tn)
    both = ownO & ownN
    dr, df = rfO[:, 0] - rfN[:, 0], rfO[:, 1] - rfN[:, 1]
    dr[~both] = np.nan; df[~both] = np.nan
    say(f"\ntemplates: own (non-fallback, non-single) modules old {int(ownO.sum())}, new {int(ownN.sum())}, both {int(both.sum())}")
    for t, n_ in (("old", to), ("new", tn)):
        say(f"  {t}: pool sum {int(n_['n_pool'][pwo].sum())}, median pool {int(np.median(n_['n_pool'][pwo]))}, "
            f"fallback==2 {int(((n_['fallback']==2)&pwo&(n_['dead']==0)).sum())}, single {int(((n_['single']==1)&pwo&(n_['dead']==0)).sum())}, "
            f"sigma_borrowed {int((n_['sigma_borrowed'][pwo]==1).sum())}, runs {json.loads(str(n_['meta']))['runs']}")
    say(f"  rise old/new median {np.nanmedian(rfO[both,0]):.3f}/{np.nanmedian(rfN[both,0]):.3f} ns, d_rise median {np.nanmedian(dr):+.3f} rms {np.nanstd(dr):.3f} ns")
    say(f"  fall old/new median {np.nanmedian(rfO[both,1]):.3f}/{np.nanmedian(rfN[both,1]):.3f} ns, d_fall median {np.nanmedian(df):+.3f} rms {np.nanstd(df):.3f} ns")
    sh_o, sh_n = to["shape"].astype(float), tn["shape"].astype(float)
    # shapes assumed peak-normalised; renormalise to be safe
    sh_o /= np.maximum(sh_o.max(1, keepdims=True), 1e-9); sh_n /= np.maximum(sh_n.max(1, keepdims=True), 1e-9)
    dsh = np.where(both, np.abs(sh_o - sh_n).max(1), np.nan)
    say(f"  max|dshape| (peak-normalised): median {np.nanmedian(dsh):.4f}, 90% {np.nanpercentile(dsh,90):.4f}, max {np.nanmax(dsh):.3f}")
    for lo, hi in ((0, 150), (150, 250), (250, 400), (400, 1e9)):
        s = both & (rr >= lo) & (rr < hi)
        if s.sum(): say(f"  r {lo}-{hi if hi<1e8 else 'inf'} mm ({int(s.sum())}): d_rise med {np.nanmedian(dr[s]):+.3f}, d_fall med {np.nanmedian(df[s]):+.3f}, max|dshape| med {np.nanmedian(dsh[s]):.4f}")
    # onset offsets / t_pk
    dtpk = np.where(both, to["t_pk"] - tn["t_pk"], np.nan)
    say(f"  t_pk old-new: median {np.nanmedian(dtpk):+.3f} ns, rms {np.nanstd(dtpk):.3f}")
    say(f"  align_rms_ns median old {np.nanmedian(to['align_rms_ns'][both]):.3f}, new {np.nanmedian(tn['align_rms_ns'][both]):.3f}")
    sgO = np.median(to["shape_sigma"][both], 0); sgN = np.median(tn["shape_sigma"][both], 0)

    # --- figures
    with PdfPages(OUT_PDF) as pdf:
        fig, ax = plt.subplots(2, 2, figsize=(14, 12))
        mapfig(ax[0, 0], names, ratio, "pulses per module: fraction old / fraction new", "RdBu_r", 0.5, 1.5, "ratio")
        ax[0, 1].hist(ratio[np.isfinite(ratio)], 60, range=(0, 3)); ax[0, 1].set(xlabel="fraction old / fraction new", title="pulse-fraction ratio per module")
        ax[1, 0].hist(rel["old"], 140, range=(0, 1.4), histtype="step", density=True, label="old"); ax[1, 0].hist(rel["new"], 140, range=(0, 1.4), histtype="step", density=True, label="new")
        ax[1, 0].set(yscale="log", xlabel="e_rough / E_beam", title="pooled e_rough"); ax[1, 0].legend()
        ax[1, 1].hist(Old["amp"], 200, range=(0, 4095), histtype="step", density=True, label="old"); ax[1, 1].hist(New["amp"], 200, range=(0, 4095), histtype="step", density=True, label="new")
        ax[1, 1].set(yscale="log", xlabel="amplitude (ADC)", title="pooled amplitude"); ax[1, 1].legend()
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        fig, ax = plt.subplots(3, 2, figsize=(14, 16))
        for k, i in enumerate(inner):
            a = ax[k, 0]
            for t, S, c in (("old", Old, "C0"), ("new", New, "C1")):
                s = S["midx"] == i
                a.hist(S["er"][s] / Eb[t], 100, range=(0, 1.3), histtype="step", density=True, label=f"{t} (E_beam {Eb[t]:.0f})", color=c)
            a.set(yscale="log", title=f"{names[i]}: e_rough/E_beam", xlabel="e_rough / E_beam"); a.legend()
            a = ax[k, 1]
            for t, S, c in (("old", Old, "C0"), ("new", New, "C1")):
                s = S["midx"] == i
                a.hist(S["amp"][s], 100, range=(0, 4095), histtype="step", density=True, label=t, color=c)
            a.set(yscale="log", title=f"{names[i]}: amplitude", xlabel="ADC"); a.legend()
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        fig, ax = plt.subplots(2, 2, figsize=(14, 12))
        ax[0, 0].plot(pN[okp], pO[okp], ".", ms=3); ax[0, 0].plot([.8, 1.1], [.8, 1.1], "k-", lw=.5)
        ax[0, 0].set(xlabel="new peak e_rough/E_beam", ylabel="old", title="e-p peak per module (mode above 0.5 E_beam)")
        mapfig(ax[0, 1], names, pr, "e-p peak ratio old/new", "RdBu_r", 0.97, 1.03, "ratio")
        ax[1, 0].hist(pr[okp], 60); ax[1, 0].set(xlabel="peak ratio old/new", title=f"median {np.nanmedian(pr):.4f}")
        ax[1, 1].hist(Old["rms"], 100, range=(0, 1.5), histtype="step", density=True, label="old"); ax[1, 1].hist(New["rms"], 100, range=(0, 1.5), histtype="step", density=True, label="new")
        ax[1, 1].set(xlabel="pedestal rms (ADC)", title="pedestal rms", yscale="log"); ax[1, 1].legend()
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        fig, ax = plt.subplots(figsize=(13, 7)); ax.axis("off")
        txt = "cut flow (PbWO4, fraction of nseen): old | new   [step ratio old|new]\n\n"
        for i, nme in enumerate(Old["cf_names"]):
            pv = max(i - 1, 0)
            txt += f"{nme:16s} {cfO[i]/cfO[0]:.5f} | {cfN[i]/cfN[0]:.5f}    [{cfO[i]/cfO[pv]:.4f} | {cfN[i]/cfN[pv]:.4f}]\n"
        ax.text(0, 1, txt, family="monospace", fontsize=14, va="top"); pdf.savefig(fig); plt.close(fig)

        fig, ax = plt.subplots(2, 3, figsize=(21, 12))
        ax[0, 0].plot(rfN[both, 0], rfO[both, 0], ".", ms=3); ax[0, 0].set(xlabel="new rise (ns)", ylabel="old", title="10-90% rise")
        ax[0, 1].plot(rfN[both, 1], rfO[both, 1], ".", ms=3); ax[0, 1].set(xlabel="new fall (ns)", ylabel="old", title="90-10% fall")
        for a in ax[0, :2]:
            lim = a.get_xlim(); a.plot(lim, lim, "k-", lw=.5)
        ax[0, 2].hist(dsh[both], 60); ax[0, 2].set(xlabel="max |dshape|", title="old-new max |shape| difference")
        mapfig(ax[1, 0], names, dr, "d_rise old-new (ns)", "RdBu_r", -1, 1, "ns")
        mapfig(ax[1, 1], names, df, "d_fall old-new (ns)", "RdBu_r", -3, 3, "ns")
        mapfig(ax[1, 2], names, dsh, "max |dshape|", "viridis", 0, np.nanpercentile(dsh, 98), "")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)

        fig, ax = plt.subplots(2, 2, figsize=(14, 12))
        ax[0, 0].plot(TG, sgO, label="old"); ax[0, 0].plot(TG, sgN, label="new"); ax[0, 0].set(xlim=(-10, 150), xlabel="t (ns)", title="array-median shape_sigma(t)"); ax[0, 0].legend()
        ax[0, 1].plot(TG, sgO / np.maximum(sgN, 1e-9)); ax[0, 1].set(xlim=(-10, 150), ylim=(0, 3), title="sigma old/new", xlabel="t (ns)")
        mapfig(ax[1, 0], names, np.where(pwo & (to["dead"] == 0), to["n_pool"], np.nan), "old: pool size", "viridis", 0, None, "", )
        mapfig(ax[1, 1], names, np.where(pwo & (tn["dead"] == 0), tn["n_pool"], np.nan), "new: pool size", "viridis", 0, None, "")
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
    open(OUT_TXT, "w").write("\n".join(OUT) + "\n")


OUT_PDF = OLD / "compare_old_vs_new.pdf"; OUT_TXT = OLD / "compare.txt"
if __name__ == "__main__":
    main()
