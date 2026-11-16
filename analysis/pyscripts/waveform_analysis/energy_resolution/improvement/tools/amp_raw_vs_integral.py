"""Template-free amplitude estimators (raw peak, cfd7, sum5) vs the benchmark integral, 2 x 3 panels (run 25320 P1b scan).

Run (from analysis/pyscripts, any python with numpy + matplotlib; no prad2py needed):
    python3 waveform_analysis/energy_resolution/improvement/tools/amp_raw_vs_integral.py [scan.npz]
Input: output/energy_resolution/scan_variants/amp_p1b/2p24GeV_H2/run_25320_s50.npz (h_apk_raw raw peak sample, h_c7 = cfd7,
h_s5 = sum5, h_e MeV, h_a template amplitude, h_k pulses in the channel, h_t time, h_mod, mod_pwo4, calib, bank).
Output: waveform_analysis/figures/amp_raw_vs_integral.{png,pdf}; numbers to stdout.

Integral = h_e / factor (ADC x sample), as in amp_vs_integral.py.  Hits: PbWO4, matched filter (h_k > 0), 100 < h_t < 200 ns.
Multi-pulse = h_k >= 2.  Phase x = ((h_t - t50) mod 4 ns) / 4 ns (t50 from the template bank, as p1b_corr_table.py).
Panel (f): per-module median of raw/integral for single-pulse hits with integral in [F_LO, F_HI], modules with >= NMOD hits.
Inset in (d): raw/h_a / module median vs phase (single-pulse hits, h_a > 200).
"""
import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.abspath(os.path.join(HERE, *[".."] * 4))
WA = os.path.join(PYDIR, "waveform_analysis")
FILE = sys.argv[1] if len(sys.argv) > 1 else f"{PYDIR}/output/energy_resolution/scan_variants/amp_p1b/2p24GeV_H2/run_25320_s50.npz"
EST = [("h_apk_raw", "raw peak", "r"), ("h_c7", "cfd7", "tab:blue"), ("h_s5", "sum5", "tab:green")]
IE = np.geomspace(20, 4e4, 121)
F_LO, F_HI, NMOD, NPH = 300.0, 1500.0, 50, 16


def rwidth(r):
    return 1.4826 * np.median(np.abs(r - np.median(r)))


def binned(x, y, edges, nmin=50):
    c = np.sqrt(edges[1:] * edges[:-1]); i = np.digitize(x, edges) - 1
    med, w = np.full(c.size, np.nan), np.full(c.size, np.nan)
    for b in range(c.size):
        s = y[i == b]
        if s.size >= nmin:
            med[b], w[b] = np.median(s), rwidth(s)
    return c, med, w


z = np.load(FILE)
fac = {e["name"]: e["factor"] for e in json.load(open(z["calib"].item()))}
names = np.array([str(n) for n in z["mod_name"]])
g = np.array([fac.get(n, 0.0) for n in names])
m = z["h_mod"].astype(int)
t = z["h_t"].astype(np.float64)
ok = (z["h_k"] > 0) & z["mod_pwo4"][m] & (g[m] > 0) & (z["h_a"] > 0) & (z["h_e"] > 0) & (t > 100) & (t < 200)
for c_, _, _ in EST:
    ok &= np.isfinite(z[c_]) & (z[c_] > 0)
I = z["h_e"][ok] / g[m[ok]]
V = {c_: z[c_][ok].astype(np.float64) for c_, _, _ in EST}
K, A, M, T = z["h_k"][ok], z["h_a"][ok].astype(np.float64), m[ok], t[ok]
mp = K >= 2
print(f"{FILE}\n{ok.sum()} hits, {mp.mean()*100:.1f} % multi-pulse (h_k>=2)")

fig, ax = plt.subplots(2, 3, figsize=(16, 9.5), constrained_layout=True)
rows = []
for a, (col, lab, _) in zip(ax[0], EST):
    v = V[col]
    ae = np.geomspace(3, 4.2e3, 121)
    H, _, _ = np.histogram2d(I, v, bins=[IE, ae])
    pm = a.pcolormesh(IE, ae, H.T, norm=LogNorm(vmin=1), cmap="viridis", rasterized=True)
    Hm, _, _ = np.histogram2d(I[mp], v[mp], bins=[IE, ae])
    a.contour(np.sqrt(IE[1:] * IE[:-1]), np.sqrt(ae[1:] * ae[:-1]), Hm.T, levels=[3, 30, 300], colors="w",
              linewidths=[0.6, 0.9, 1.2])
    c, med, _ = binned(I, v, IE)
    a.plot(c, med, "r", lw=1.6, label="median per bin")
    a.plot([], [], "w", lw=1, label=f"multi-pulse (h_k>=2) contours, {mp.mean()*100:.1f} %")
    a.set_xscale("log"); a.set_yscale("log"); a.set_xlim(IE[0], IE[-1]); a.set_ylim(ae[0], ae[-1])
    a.set_xlabel("integral (ADC x sample)"); a.set_ylabel(f"{lab} (ADC)")
    a.set_title(f"({'abc'[len(rows)]}) {lab} vs integral; median ratio {np.median(v/I):.4f}", fontsize=10)
    a.legend(loc="upper left", fontsize=8, facecolor="0.25", labelcolor="w", framealpha=0.8)
    fig.colorbar(pm, ax=a, pad=0.01, label="hits / bin")
    rows.append((lab, np.median(v / I), rwidth(v / I)))

# (d) raw/integral vs integral
a = ax[1, 0]
r = V["h_apk_raw"] / I
c, med, w = binned(I, r, IE)
a.plot(c, med, "r", lw=1.6, label=f"raw peak: med {np.median(r):.4f}, width {rwidth(r):.4f} ({rwidth(r)/np.median(r)*100:.1f} %)")
a.fill_between(c, med - w, med + w, color="r", alpha=0.2)
a.set_xscale("log"); a.set_xlabel("integral (ADC x sample)"); a.set_ylabel("raw peak / integral (1/sample)")
a.set_ylim(0.10, 0.20)
a.set_title("(d) raw/integral vs integral, median +- robust width", fontsize=10)
a.legend(loc="lower left", fontsize=8); a.grid(alpha=0.3)

# (e) all three
a = ax[1, 1]
for col, lab, cl in EST:
    r = V[col] / I
    c, med, w = binned(I, r, IE)
    a.plot(c, med, color=cl, lw=1.6, label=f"{lab}: med {np.median(r):.4f}, width {rwidth(r):.4f} ({rwidth(r)/np.median(r)*100:.1f} %)")
    a.fill_between(c, med - w, med + w, color=cl, alpha=0.15)
a.set_xscale("log"); a.set_yscale("log"); a.set_ylim(0.1, 0.9); a.set_xlabel("integral (ADC x sample)"); a.set_ylabel("estimator / integral (1/sample)")
a.set_title("(e) ratio vs integral, three estimators", fontsize=10)
a.legend(loc="lower right", fontsize=8); a.grid(alpha=0.3, which="both")

# (f) per-module median at fixed integral
a = ax[1, 2]
sel = (K == 1) & (I > F_LO) & (I < F_HI)
Ms, rs = M[sel], (V["h_apk_raw"] / I)[sel]
o = np.argsort(Ms, kind="stable"); Ms, rs = Ms[o], rs[o]
u, st, cnt = np.unique(Ms, return_index=True, return_counts=True)
mr = np.array([np.median(rs[s:s + n]) for s, n in zip(st, cnt) if n >= NMOD])
a.hist(mr, bins=40, color="tab:orange")
mmed, msd = np.median(mr), rwidth(mr)
a.set_xlabel(f"module median raw/integral, integral {F_LO:.0f}-{F_HI:.0f}, h_k=1"); a.set_ylabel("modules")
a.set_title(f"(f) {mr.size} modules: median {mmed:.4f}, robust sd {msd:.4f} ({msd/mmed*100:.1f} %)", fontsize=10)

# inset in (d): phase saw-tooth
t50 = None
try:
    b = np.load(os.path.join(PYDIR, str(z["bank"])), allow_pickle=False)
    okb = b["valid"] > 0
    row = {str(s): i for i, s in enumerate(b["names"]) if okb[i]}
    med50 = float(np.median(b["t50"][okb]))
    t50 = np.array([b["t50"][row[n]] if n in row else med50 for n in names], dtype=np.float64)
except Exception as e:
    print("phase: bank unreadable:", e)
amp_ph = np.nan
if t50 is not None:
    s1 = (K == 1) & (A > 200)
    x = ((T[s1] - t50[M[s1]]) % 4.0) / 4.0
    q = V["h_apk_raw"][s1] / A[s1]
    mm = M[s1]
    scale = np.full(len(names), np.nan)
    for j in np.unique(mm):
        k = mm == j
        if k.sum() >= 50:
            scale[j] = np.median(q[k])
    qn = q / scale[mm]; good = np.isfinite(qn)
    b_ = np.minimum((x[good] * NPH).astype(int), NPH - 1)
    curve = np.array([np.median(qn[good][b_ == i]) for i in range(NPH)])
    amp_ph = curve.max() - curve.min()
    ins = ax[1, 0].inset_axes([0.62, 0.52, 0.34, 0.25])
    ins.plot((np.arange(NPH) + 0.5) / NPH, curve, "ko-", ms=3)
    ins.set_xlabel("onset phase", fontsize=7); ins.set_ylabel("raw/h_a (norm.)", fontsize=7)
    ins.tick_params(labelsize=7); ins.grid(alpha=0.3)
    print("phase curve", np.round(curve, 4))
fig.suptitle("Template-free amplitude estimators vs integral: run 25320 s50, PbWO4 MF hits, 100 < t < 200 ns")
for ext in ("png", "pdf"):
    fig.savefig(f"{WA}/figures/amp_raw_vs_integral.{ext}", dpi=150)
for lab, md, w in rows:
    print(f"{lab:9s} median ratio to integral {md:.5f}  robust width {w:.5f} ({w/md*100:.2f} %)")
print(f"per-module raw/integral at integral {F_LO:.0f}-{F_HI:.0f} (h_k=1, >={NMOD} hits): {mr.size} modules, median {mmed:.5f}, "
      f"robust sd {msd:.5f} ({msd/mmed*100:.2f} %), range {mr.min():.5f}..{mr.max():.5f}")
print(f"phase saw-tooth raw/h_a peak-to-peak {amp_ph:.4f} (semi-amplitude {amp_ph/2:.4f})")
