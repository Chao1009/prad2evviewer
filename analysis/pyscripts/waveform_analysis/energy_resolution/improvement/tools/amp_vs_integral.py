"""Fitted pulse amplitude vs pulse integral of the matched-filter pulses (benchmark scan banks), 2 x 3 panels.

Run (from analysis/pyscripts, any python with numpy + matplotlib; no prad2py needed):
    python3 waveform_analysis/energy_resolution/improvement/tools/amp_vs_integral.py
Inputs: output/energy_resolution/scan/<run type>/*.npz (keys h_a amplitude [ADC], h_e energy [MeV], h_k pulses in the
channel, h_mod, mod_name, mod_pwo4, calib), and the calibration JSON named in each file (factor = MeV per integral unit).
Output: waveform_analysis/figures/amp_vs_integral.{png,pdf}; the numbers go to stdout.

Integral = h_e / factor = the benchmark's  iint_wa x I(phi)/mean(I) x amplitude, in ADC x sample (1 sample = 4 ns); so
amplitude / integral = 1 / (module template integral x onset-phase factor): a template property, not an independent
measurement.  Only matched-filter PbWO4 hits (h_k > 0; WaveAnalyzer fallback hits have no amplitude), a random sample of
at most N_PER_FILE per scan file.  Multi-pulse = h_k > 1 (more than one pulse fitted in the channel waveform).
"""
import glob
import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

HERE = os.path.dirname(os.path.abspath(__file__))
PYDIR = os.path.abspath(os.path.join(HERE, *[".."] * 4))
WA = os.path.join(PYDIR, "waveform_analysis")
RTS = [("0p73GeV_H2", "0.73 GeV"), ("2p24GeV_H2", "2.24 GeV"), ("3p49GeV_H2", "3.49 GeV"),
       ("3p49GeV_H2_old", "3.49 GeV old"), ("2p24GeV_Ta_X17", "Ta X17 (2.24 GeV)")]
N_PER_FILE, SEED = 400_000, 1
IE = np.geomspace(20, 4e4, 121)     # integral bins
AE = np.geomspace(3, 4.2e3, 121)  # amplitude bins


def load(rt):
    I, A, K, M = [], [], [], []
    rng = np.random.default_rng(SEED)
    for f in sorted(glob.glob(f"{PYDIR}/output/energy_resolution/scan/{rt}/*.npz")):
        z = np.load(f)
        fac = {e["name"]: e["factor"] for e in json.load(open(z["calib"].item()))}
        names = z["mod_name"]
        g = np.array([fac.get(str(n), 0.0) for n in names])
        m = z["h_mod"].astype(int)
        ok = (z["h_k"] > 0) & z["mod_pwo4"][m] & (g[m] > 0) & (z["h_a"] > 0) & (z["h_e"] > 0)
        idx = np.nonzero(ok)[0]
        idx = np.sort(rng.choice(idx, N_PER_FILE, replace=False)) if idx.size > N_PER_FILE else idx
        mm = m[idx]
        I.append(z["h_e"][idx] / g[mm]); A.append(z["h_a"][idx]); K.append(z["h_k"][idx])
        M.append(np.array([str(n) for n in names])[mm])
    return np.concatenate(I), np.concatenate(A), np.concatenate(K), np.concatenate(M)


def binned_median(x, y, edges, nmin=50):
    c = np.sqrt(edges[1:] * edges[:-1])
    i = np.digitize(x, edges) - 1
    med = np.full(c.size, np.nan)
    for b in range(c.size):
        s = y[i == b]
        if s.size >= nmin:
            med[b] = np.median(s)
    return c, med


def rwidth(r):
    return 1.4826 * np.median(np.abs(r - np.median(r)))


fig, ax = plt.subplots(2, 3, figsize=(16, 9.5), constrained_layout=True)
cols = plt.cm.tab10(np.arange(5))
rows, out = [], []
data = {}
for k, (rt, lab) in enumerate(RTS):
    I, A, K, M = load(rt)
    data[rt] = (I, A, K, M)
    a = ax.flat[k]
    H, _, _ = np.histogram2d(I, A, bins=[IE, AE])
    pm = a.pcolormesh(IE, AE, H.T, norm=LogNorm(vmin=1), cmap="viridis", rasterized=True)
    mp = K > 1
    if mp.sum() > 500:
        Hm, _, _ = np.histogram2d(I[mp], A[mp], bins=[IE, AE])
        cI, cA = np.sqrt(IE[1:] * IE[:-1]), np.sqrt(AE[1:] * AE[:-1])
        a.contour(cI, cA, Hm.T, levels=[3, 30, 300], colors="w", linewidths=[0.6, 0.9, 1.2], linestyles="-")
    c, med = binned_median(I, A, IE)
    a.plot(c, med, color="r", lw=1.6, label="median amplitude per bin")
    r = A / I
    a.plot([], [], color="w", lw=1.0, label=f"multi-pulse (k>1) contours, {mp.mean()*100:.1f} %")
    a.set_xscale("log"); a.set_yscale("log")
    a.set_xlim(IE[0], IE[-1]); a.set_ylim(AE[0], AE[-1])
    a.set_xlabel("pulse integral (ADC x sample, 4 ns)"); a.set_ylabel("fitted amplitude (ADC)")
    a.set_title(f"{lab}: {len(I)/1e6:.2f} M pulses; A/I = {np.median(r):.4f}", fontsize=10)
    leg = a.legend(loc="upper left", fontsize=8, facecolor="0.25", labelcolor="w", framealpha=0.8)
    fig.colorbar(pm, ax=a, pad=0.01, label="pulses / bin")
    rows.append((lab, len(I), np.median(r), rwidth(r), mp.mean()))

a = ax.flat[5]
for k, (rt, lab) in enumerate(RTS):
    I, A, K, M = data[rt]
    r = A / I
    c, med = binned_median(I, r, IE)
    a.plot(c, med, color=cols[k], lw=1.6,
           label=f"{lab}: med {np.median(r):.4f}, 1.4826 MAD {rwidth(r):.4f} ({rwidth(r)/np.median(r)*100:.1f} %)")
a.set_xscale("log"); a.set_xlim(IE[0], IE[-1]); a.set_ylim(0.120, 0.160)
a.set_xlabel("pulse integral (ADC x sample, 4 ns)"); a.set_ylabel("amplitude / integral (1/sample)")
a.set_title("ratio vs integral, median per integral bin", fontsize=10)
a.legend(loc="upper left", fontsize=7)
a.grid(alpha=0.3)

# per-module median ratio at 2.24 GeV H2
I, A, K, M = data["2p24GeV_H2"]
mods, inv = np.unique(M, return_inverse=True)
r = A / I
mr = np.array([np.median(r[inv == j]) for j in range(len(mods)) if (inv == j).sum() >= 200])
ins = a.inset_axes([0.60, 0.10, 0.37, 0.30])
ins.hist(mr, bins=40, color=cols[1])
ins.set_xlabel("per-module median A/I, 2.24 GeV", fontsize=7); ins.set_ylabel("modules", fontsize=7)
ins.tick_params(labelsize=7)
mod_med, mod_sd = np.median(mr), rwidth(mr)
mod_rng = (mr.min(), mr.max())
fig.suptitle("Matched-filter pulses: fitted amplitude vs integral (sample of the scan banks, PbWO4, 100 < t < 200 ns)")
for ext in ("png", "pdf"):
    fig.savefig(f"{WA}/figures/amp_vs_integral.{ext}", dpi=150)
print(f"{'run type':20s} {'N':>9s} {'median A/I':>11s} {'1.4826MAD':>10s} {'rel %':>6s} {'multi %':>8s}")
for lab, n, m, w, f in rows:
    print(f"{lab:20s} {n:9d} {m:11.5f} {w:10.5f} {w/m*100:6.2f} {f*100:8.2f}")
print(f"2.24 GeV per-module median A/I ({mr.size} modules, >=200 pulses): median {mod_med:.5f}, robust sd {mod_sd:.5f} "
      f"({mod_sd/mod_med*100:.2f} %), range {mod_rng[0]:.5f}..{mod_rng[1]:.5f}")
