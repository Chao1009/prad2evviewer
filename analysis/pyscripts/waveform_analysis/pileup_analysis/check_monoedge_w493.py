#!/usr/bin/env python3
"""check_monoedge_w493.py -- read-only checks on existing pileup scans (no refit).

Check 1: pulses removed by the monotonic-edge rule (old vs new bounded sample): dt = t_removed - t_main, r = A_removed/A_main.
Check 2: W493 (+3x3 neighbours), W121 and PbWO4 average onset-time distributions in the full scans.

    python3 pileup_analysis/check_monoedge_w493.py
-> checks_monoedge_removed_W493_timing_20261116.pdf, .tsv
Window: 100 samples x 4 ns = 400 ns; t0+t_pk>372 ns excluded, t0>284 ns flagged truncated in the fitter.
"""
import sys, json, time
from pathlib import Path
from collections import defaultdict
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

sys.dont_write_bytecode = True
HERE = Path(__file__).resolve().parent
OUT = HERE.parent.parent / "output" / "pileup_rates"
MAP = HERE.parent.parent.parent.parent / "database" / "hycal_map.json"
STEM = HERE / "checks_monoedge_removed_W493_timing_20261116"
TOL = 2.0
T0 = time.time()
rows = []   # tsv rows: section, group, quantity, value


def tsv(*r):
    rows.append(r)


# ---------------------------------------------------------------------------------------------------- check 1
def pulses(z, sel):
    t, a, k = z["ch_t"][sel], z["ch_a"][sel], z["ch_k"][sel]
    return t, a, k


def check1():
    rem, sur = defaultdict(list), defaultdict(list)   # run_type -> list of (dt, r, Amain)
    nmatch = nchg = 0
    for fn in sorted((OUT / "monoedge_sample" / "new").glob("run_*.npz")):
        zn = np.load(fn)
        zo = np.load(OUT / "monoedge_sample" / "old" / fn.name)
        rt = str(zn["run_type"])
        kn_all = {q: zn[q] for q in ("ch_ev", "ch_midx", "ch_k")}
        ko_all = {q: zo[q] for q in ("ch_ev", "ch_midx", "ch_k")}
        if (kn_all["ch_ev"].shape == ko_all["ch_ev"].shape and np.array_equal(kn_all["ch_ev"], ko_all["ch_ev"])
                and np.array_equal(kn_all["ch_midx"], ko_all["ch_midx"])):
            io = in_ = np.arange(ko_all["ch_ev"].size)
        else:   # match by (ev, midx) key
            ko = ko_all["ch_ev"].astype(np.int64) * 4096 + ko_all["ch_midx"]
            kn = kn_all["ch_ev"].astype(np.int64) * 4096 + kn_all["ch_midx"]
            _, io, in_ = np.intersect1d(ko, kn, return_indices=True)
        nmatch += io.size
        ko_k, kn_k = zo["ch_k"][io], zn["ch_k"][in_]
        nchg += int((ko_k > kn_k).sum())
        s = ko_k >= 2
        to, ao = zo["ch_t"][io][s].astype(np.float64), zo["ch_a"][io][s].astype(np.float64)
        tn, kn_s = zn["ch_t"][in_][s].astype(np.float64), kn_k[s]
        ko_s = ko_k[s]
        col = np.arange(6)[None, :]
        vo = col < ko_s[:, None]
        vn = col < kn_s[:, None]
        ao = np.where(vo, ao, -1.0)
        im = ao.argmax(1)
        r_ = np.arange(ao.shape[0])
        tm, am = to[r_, im], ao[r_, im]
        d = np.abs(to[:, :, None] - np.where(vn, tn, np.inf)[:, None, :])
        surv = (d <= TOL).any(2)
        sec = vo & (col != im[:, None])
        for name, m, store in (("rem", sec & ~surv, rem), ("sur", sec & surv, sur)):
            rr, cc = np.nonzero(m)
            store[rt].append(np.c_[to[rr, cc] - tm[rr], ao[rr, cc] / am[rr], am[rr]])
        print(f"  check1 run {int(zn['run'])} {rt}: matched {io.size}, k decreased {(ko_k > kn_k).sum()}", flush=True)
    return {n: {rt: np.vstack(v) for rt, v in s.items()} for n, s in (("rem", rem), ("sur", sur))}, nmatch, nchg


def plot1(pdf, res, nmatch, nchg):
    rem, sur = res["rem"], res["sur"]
    groups = {"all": None, "Ta_X17": lambda rt: "Ta" in rt, "non-Ta": lambda rt: "Ta" not in rt}
    # tsv + title page
    lines = [f"Check 1: pulses removed by monotonic-edge rule (sample old vs new)",
             f"matched channels {nmatch}, with k_old>k_new {nchg}; pulse match |dt|<={TOL} ns; main = largest A of old fit",
             "dt = t_removed - t_main (onset, ns); r = A_removed/A_main",
             "protected fall window = 0 < dt <= 40 ns (onset-relative).",
             "Criterion: mismatch splits -> small r, dt clustered just after peak; real pile-up -> broad dt, larger r."]
    fig = plt.figure(figsize=(11, 8.5)); fig.text(0.05, 0.95, "\n".join(lines), va="top", family="monospace", fontsize=9)
    ytxt = []
    summ = {}
    for gname, sel in groups.items():
        rt_list = [rt for rt in sorted(set(rem) | set(sur)) if sel is None or sel(rt)]
        R = np.vstack([rem[rt] for rt in rt_list if rt in rem]) if any(rt in rem for rt in rt_list) else np.zeros((0, 3))
        S = np.vstack([sur[rt] for rt in rt_list if rt in sur]) if any(rt in sur for rt in rt_list) else np.zeros((0, 3))
        summ[gname] = (R, S)
    for rt in sorted(set(rem) | set(sur)):
        summ[rt] = (rem.get(rt, np.zeros((0, 3))), sur.get(rt, np.zeros((0, 3))))
    hdr = f"{'group':14s} {'n_rem':>8s} {'f(0<dt<=40)':>12s} {'r<0.05':>8s} {'r<0.1':>8s} {'r<0.2':>8s} {'n r>0.2':>8s} | {'n_surv':>9s} {'f(0-40)':>8s} {'r<0.2':>7s} {'med dt rem':>10s} {'med dt sur':>10s}"
    ytxt.append(hdr)
    for g, (R, S) in summ.items():
        def fr(X, c):
            return (c(X).mean() if len(X) else float("nan"))
        w = lambda X: (X[:, 0] > 0) & (X[:, 0] <= 40)
        vals = dict(n_removed=len(R), f_dt_0_40=fr(R, w), f_r_lt_0p05=fr(R, lambda X: X[:, 1] < .05), f_r_lt_0p1=fr(R, lambda X: X[:, 1] < .1),
                    f_r_lt_0p2=fr(R, lambda X: X[:, 1] < .2), n_r_gt_0p2=int((R[:, 1] >= .2).sum()), n_survive=len(S),
                    surv_f_dt_0_40=fr(S, w), surv_f_r_lt_0p2=fr(S, lambda X: X[:, 1] < .2),
                    med_dt_removed=float(np.median(R[:, 0])) if len(R) else float("nan"), med_dt_survive=float(np.median(S[:, 0])) if len(S) else float("nan"),
                    med_r_removed=float(np.median(R[:, 1])) if len(R) else float("nan"), med_r_survive=float(np.median(S[:, 1])) if len(S) else float("nan"))
        for q, v in vals.items():
            tsv("check1", g, q, v)
        ytxt.append(f"{g:14s} {vals['n_removed']:8d} {vals['f_dt_0_40']:12.3f} {vals['f_r_lt_0p05']:8.3f} {vals['f_r_lt_0p1']:8.3f} {vals['f_r_lt_0p2']:8.3f} {vals['n_r_gt_0p2']:8d} | {vals['n_survive']:9d} {vals['surv_f_dt_0_40']:8.3f} {vals['surv_f_r_lt_0p2']:7.3f} {vals['med_dt_removed']:10.1f} {vals['med_dt_survive']:10.1f}")
    fig.text(0.03, 0.55, "\n".join(ytxt), va="top", family="monospace", fontsize=6.3)
    pdf.savefig(fig); plt.close(fig)
    print("\n".join(ytxt))
    # per run type plots
    for g, (R, S) in summ.items():
        if g in ("all", "non-Ta") and len(summ) > 3 and False:
            continue
        fig, ax = plt.subplots(2, 2, figsize=(11, 8.5))
        fig.suptitle(f"Check 1: {g}   removed n={len(R)}, surviving n={len(S)}")
        a = ax[0, 0]
        lr = np.log10(np.clip(R[:, 1], 1e-4, None)) if len(R) else R[:, 1]
        if len(R):
            h = a.hist2d(R[:, 0], lr, bins=[np.arange(-100, 300, 4), np.linspace(-4, 0.5, 60)], cmin=1, norm=matplotlib.colors.LogNorm())
            plt.colorbar(h[3], ax=a)
        a.axvspan(0, 40, color="#f4a259", alpha=.15); a.axhline(np.log10(.2), c="r", lw=.7)
        a.set_xlabel("dt = t_removed - t_main (ns)"); a.set_ylabel("log10 r"); a.set_title("removed: dt vs r")
        a = ax[0, 1]
        if len(S):
            h = a.hist2d(S[:, 0], np.log10(np.clip(S[:, 1], 1e-4, None)), bins=[np.arange(-100, 300, 4), np.linspace(-4, 0.5, 60)], cmin=1, norm=matplotlib.colors.LogNorm())
            plt.colorbar(h[3], ax=a)
        a.axvspan(0, 40, color="#f4a259", alpha=.15); a.axhline(np.log10(.2), c="r", lw=.7)
        a.set_xlabel("dt (ns)"); a.set_ylabel("log10 r"); a.set_title("surviving: dt vs r")
        a = ax[1, 0]
        b = np.arange(-100, 300, 4)
        for X, lab, c in ((R, "removed", "C3"), (S, "surviving", "C0")):
            if len(X):
                a.hist(X[:, 0], bins=b, histtype="step", density=True, label=lab, color=c)
        a.axvspan(0, 40, color="#f4a259", alpha=.15); a.set_yscale("log"); a.set_xlabel("dt (ns)"); a.set_title("dt (normalized)"); a.legend()
        a = ax[1, 1]
        b = np.logspace(-4, 0.5, 50)
        for X, lab, c in ((R, "removed", "C3"), (S, "surviving", "C0")):
            if len(X):
                a.hist(X[:, 1], bins=b, histtype="step", density=True, label=lab, color=c)
        a.axvline(.2, c="r", lw=.7); a.set_xscale("log"); a.set_yscale("log"); a.set_xlabel("r = A_rem/A_main"); a.set_title("r (normalized)"); a.legend()
        fig.tight_layout(rect=(0, 0, 1, .96)); pdf.savefig(fig); plt.close(fig)
    # removed by A_main
    fig, ax = plt.subplots(figsize=(11, 8.5))
    for g in [g for g in summ if g not in ("all", "non-Ta", "Ta_X17")]:
        R = summ[g][0]
        if len(R):
            ax.hist(np.log10(np.clip(R[:, 2], 1, None)), bins=40, histtype="step", label=f"{g} (n={len(R)})")
    ax.set_xlabel("log10 A_main (ADC)"); ax.set_yscale("log"); ax.set_title("removed pulses: main amplitude"); ax.legend(fontsize=7)
    pdf.savefig(fig); plt.close(fig)


# ---------------------------------------------------------------------------------------------------- check 2
BINS = np.arange(0, 404, 4.0)


def check2():
    mp = json.load(open(MAP))
    byname = {m["n"]: m for m in mp}
    c = byname["W493"]["geo"]
    pitch = byname["W493"]["geo"]["sx"]

    def nbrs(n):
        g = byname[n]["geo"]
        return [m["n"] for m in mp if m["t"] == "PbWO4" and abs(m["geo"]["x"] - g["x"]) < 1.5 * pitch and abs(m["geo"]["y"] - g["y"]) < 1.5 * pitch]
    n493 = nbrs("W493")
    print("W493 3x3:", n493)
    mods = sorted(set(n493 + ["W121"]))
    H = defaultdict(lambda: np.zeros(BINS.size - 1))
    C = defaultdict(lambda: np.zeros(4))   # n, >350, 140-156, >284
    for fn in sorted((OUT / "scan").glob("run_*.npz")):
        z = np.load(fn)
        rt = str(z["run_type"]); names = z["mod_name"]; pwo = z["mod_pwo4"]
        midx, k, t = z["ch_midx"], z["ch_k"], z["ch_t"]
        col = np.arange(6)[None, :]
        valid = col < k[:, None]
        tv = np.where(valid, t, np.nan)
        for g, msk in [("PbWO4", pwo[midx])] + [(m, midx == int(np.nonzero(names == m)[0][0])) for m in mods if (names == m).any()]:
            x = tv[msk].ravel(); x = x[np.isfinite(x)]
            H[(rt, g)] += np.histogram(np.clip(x, 0, 399.9), BINS)[0]
            C[(rt, g)] += [x.size, (x > 350).sum(), ((x >= 140) & (x <= 156)).sum(), (x > 284).sum()]
        print(f"  check2 {int(z['run'])} {rt} {time.time() - T0:.0f}s", flush=True)
    return H, C, n493, mods


def plot2(pdf, H, C, n493):
    rts = sorted({k[0] for k in H})
    for g in sorted({k[1] for k in H}):
        for rt in rts:
            if (rt, g) in C:
                n, a, b, c = C[(rt, g)]
                tsv("check2", f"{rt}|{g}", "n_pulses", int(n)); tsv("check2", f"{rt}|{g}", "f_t0_gt350", a / max(n, 1))
                tsv("check2", f"{rt}|{g}", "f_t0_140_156", b / max(n, 1)); tsv("check2", f"{rt}|{g}", "f_t0_gt284", c / max(n, 1))
    for rt in rts:
        fig, ax = plt.subplots(1, 2, figsize=(11, 8.5), gridspec_kw=dict(width_ratios=[2, 1.2]))
        a = ax[0]
        x = BINS[:-1] + 2
        for g in sorted({k[1] for k in H if k[0] == rt}):
            h = H[(rt, g)]
            y = h / max(h.sum(), 1)
            if g == "W493": a.step(x, y, c="r", lw=2, label="W493")
            elif g == "W121": a.step(x, y, c="b", lw=1.6, label="W121")
            elif g == "PbWO4": a.step(x, y, c="k", lw=2, label="PbWO4 avg")
            else: a.step(x, y, c="gray", lw=.7, label="W493 nbrs" if g == sorted([m for m in n493 if m != "W493"])[0] else None)
        for v, ls in ((140, ":"), (156, ":"), (284, "--"), (350, "-."), (372, "-")):
            a.axvline(v, c="g", ls=ls, lw=.8)
        a.set_yscale("log"); a.set_xlabel("pulse onset t0 (ns); last bin includes overflow"); a.set_ylabel("fraction of pulses / 4 ns")
        a.set_title(f"{rt}: green 140/156 window, 284 trunc flag, 350, 372 excl."); a.legend(fontsize=8)
        a = ax[1]; a.axis("off")
        L = [f"{'module':10s} {'n':>9s} {'>350':>8s} {'140-156':>8s} {'>284':>8s}"]
        for g in ["W493"] + [m for m in n493 if m != "W493"] + ["W121", "PbWO4"]:
            if (rt, g) in C:
                n, aa, b, c = C[(rt, g)]
                L.append(f"{g:10s} {int(n):9d} {aa / max(n, 1):8.4f} {b / max(n, 1):8.4f} {c / max(n, 1):8.4f}")
        a.text(0, 1, "\n".join(L), va="top", family="monospace", fontsize=7)
        fig.tight_layout(); pdf.savefig(fig); plt.close(fig)
        print(rt, "\n".join(L))


def main():
    res, nmatch, nchg = check1()
    H, C, n493, mods = check2()
    with PdfPages(str(STEM) + ".pdf") as pdf:
        plot1(pdf, res, nmatch, nchg)
        plot2(pdf, H, C, n493)
    with open(str(STEM) + ".tsv", "w") as f:
        f.write("section\tgroup\tquantity\tvalue\n")
        for r in rows:
            f.write("\t".join(str(x) if not isinstance(x, float) else f"{x:.5g}" for x in r) + "\n")
    print(f"done in {time.time() - T0:.0f}s")


if __name__ == "__main__":
    main()
