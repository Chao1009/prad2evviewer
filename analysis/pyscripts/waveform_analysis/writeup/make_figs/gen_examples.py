#!/usr/bin/env python3
"""gen_examples.py -- example stage-1 pile-up waveforms with their truth components (read-only).

    python3 waveform_analysis/writeup/make_figs/gen_examples.py [stage1.npz] [bank.npz]

Reads the first rows of output/pileup_stage1/pileup_stage1.npz (written by
pileup_analysis/pileup_generator_stage1.py) by streaming the zip members, so the
1.9 GB file is never fully decompressed.  Nothing is generated or written except
writeup/figs/gen_examples.pdf.

Events cycle through the generator's (module, A1, (r, dt)) combinations with (r, dt)
fastest; the first 141 events are W493 at A1 = 300 ADC, the next 141 W493 at A1 = 1000 ADC
(control first, then r x dt).  The plotted events are picked from the stored truth fields.

Truth components are drawn as ped + A1 s(t - t1) and ped + A2 s(t - t2) with s the module's
unit-peak spline, evaluated as the generator's Bank class does (B-spline on a 0.05 ns grid,
end value held, linear interpolation, s = 0 outside (0, T_end)).  Two limits: (1) the
per-event template smearing delta(t) is not stored, so the curves are the unsmeared
components; (2) the bank that generated the file (meta bank_sha1) is no longer on disk,
so s comes from the bank given on the command line (default: the current 2.24 GeV H2
bank).  The script prints both sha1 values and checks the bank's a_band of each plotted
module against the a_band stored in meta (a changed bank build would change it); per
event it prints the rms of (waveform - unsmeared sum) / noise_sigma over the pre-onset
samples (about 1 when the pedestal and noise are as recorded) and the mean residual over
the pulse region in units of A1 (about 0 when the template shape is the generator's).
"""
from __future__ import annotations

import hashlib
import json
import sys
import zipfile
from pathlib import Path

import numpy as np

PYDIR = Path(__file__).resolve().parents[3]                    # analysis/pyscripts
NPZ = PYDIR / "output" / "pileup_stage1" / "pileup_stage1.npz"
BANK = PYDIR / "waveform_analysis" / "pulse_templates" / "2p24GeV_H2" / "templates_2p24GeV_H2.npz"
OUT = PYDIR / "waveform_analysis" / "writeup" / "figs" / "gen_examples.pdf"
NROWS = 290                    # rows read from the start of the file (covers W493, A1 = 300 and 1000)
LUT_DT = 0.05                  # ns, as in the generator
CLK_NS = 4.0
# (module, A1, r, dt): the six panels
PICK = [("W493", 1000.0, 0.0, 0.0), ("W493", 1000.0, 1.0, 16.0), ("W493", 1000.0, 0.3, 48.0),
        ("W493", 1000.0, 0.1, 8.0), ("W493", 1000.0, 0.5, -32.0), ("W493", 1000.0, 2.0, -64.0)]
BLUE, ORANGE, VERM = "#0072B2", "#E69F00", "#D55E00"          # Okabe-Ito


def head_rows(zf, name, n):
    """First n rows of a C-order .npy member, read from the zip stream."""
    with zf.open(name + ".npy") as f:
        ver = np.lib.format.read_magic(f)
        rd = np.lib.format.read_array_header_1_0 if ver == (1, 0) else np.lib.format.read_array_header_2_0
        shape, fortran, dtype = rd(f)
        if fortran:
            raise SystemExit(f"{name}: Fortran order not supported")
        if shape == ():
            return np.frombuffer(f.read(dtype.itemsize), dtype=dtype)[0]
        per = int(np.prod(shape[1:])) if len(shape) > 1 else 1
        m = min(n, shape[0])
        buf = f.read(m * per * dtype.itemsize)
        return np.frombuffer(buf, dtype=dtype).reshape((m,) + tuple(shape[1:]))


def bspline(t, c, k, x):
    """De Boor evaluation of the B-spline (t, c, k) at x inside [t[k], t[-k-1]] (numpy only;
    same values as scipy.interpolate.BSpline there)."""
    t, c, x = (np.asarray(v, dtype=np.float64) for v in (t, c, x))
    n = c.size
    i = np.clip(np.searchsorted(t, x, side="right") - 1, k, n - 1)
    d = np.stack([c[i - k + j] for j in range(k + 1)], axis=0)
    for r in range(1, k + 1):
        for j in range(k, r - 1, -1):
            lo, hi = t[i - k + j], t[i + 1 + j - r]
            a = np.where(hi > lo, (x - lo) / np.where(hi > lo, hi - lo, 1.0), 0.0)
            d[j] = (1 - a) * d[j - 1] + a * d[j]
    return d[k]


def shape_fn(bank_path, name):
    z = np.load(bank_path, allow_pickle=False)
    k, t_end = int(z["spline_k"]), float(z["T_end"])
    r = [str(s) for s in z["names"]].index(name)
    t = z["spline_t"][r, :z["spline_nt"][r]]
    c = z["spline_c"][r, :z["spline_nc"][r]]
    n_lut = int(round(t_end / LUT_DT)) + 1
    lut = bspline(t, c, k, np.minimum(LUT_DT * np.arange(n_lut), t[-k - 1]))
    lut = lut.astype(np.float32).astype(np.float64)

    def f(tt):
        x = np.asarray(tt, dtype=np.float64) / LUT_DT
        i = np.clip(np.floor(x), 0, n_lut - 2).astype(np.int64)
        fr = x - i
        v = lut[i] * (1 - fr) + lut[i + 1] * fr
        return np.where((tt > 0) & (tt < t_end), v, 0.0)
    return f


def main():
    npz = Path(sys.argv[1]) if len(sys.argv) > 1 else NPZ
    bank = Path(sys.argv[2]) if len(sys.argv) > 2 else BANK
    zf = zipfile.ZipFile(npz)
    meta = json.loads(str(head_rows(zf, "meta", 1)))
    d = {k: head_rows(zf, k, NROWS) for k in
         ("waveforms", "module", "a1", "a2", "ratio", "dt_ns", "t1_ns", "t2_ns", "pedestal", "noise_sigma", "smear_rms")}
    sha = hashlib.sha1(bank.read_bytes()).hexdigest()
    print(f"[bank] generation sha1 {meta['bank_sha1']}; plotted with {bank} sha1 {sha}"
          + ("" if sha == meta["bank_sha1"] else "  (different file)"))
    s = {nm: shape_fn(bank, nm) for nm in sorted({p[0] for p in PICK})}
    zb = np.load(bank, allow_pickle=False)
    for nm in s:
        ab = float(zb["a_band"][[str(x) for x in zb["names"]].index(nm)])
        print(f"[bank] {nm} a_band bank {ab!r} meta {meta['smear']['a_band'][nm]!r} "
              + ("equal" if ab == meta["smear"]["a_band"][nm] else "DIFFERENT"))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({"font.size": 10, "axes.labelsize": 11, "legend.fontsize": 9,
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axs = plt.subplots(2, 3, figsize=(10.5, 6.2), sharex=True)
    tt = np.arange(0.0, 400.0, 0.1)
    ts = CLK_NS * np.arange(d["waveforms"].shape[1])
    for ax, (nm, a1, r, dt) in zip(axs.flat, PICK):
        hit = np.nonzero((d["module"] == nm) & (d["a1"] == np.float32(a1)) & (d["ratio"] == np.float32(r))
                         & (d["dt_ns"] == np.float32(dt)))[0]
        if not hit.size:
            raise SystemExit(f"no event for {(nm, a1, r, dt)} in the first {NROWS} rows")
        e = int(hit[0])
        ped, t1, t2, A1, A2 = (float(d[k][e]) for k in ("pedestal", "t1_ns", "t2_ns", "a1", "a2"))
        prim = A1 * s[nm](tt - t1)
        sec = A2 * s[nm](tt - t2) if r > 0 else np.zeros_like(tt)
        y = d["waveforms"][e].astype(np.float64)
        model = ped + A1 * s[nm](ts - t1) + (A2 * s[nm](ts - t2) if r > 0 else 0.0)
        first = min(t1, t2) if r > 0 else t1
        pre = ts <= first
        reg = (ts > first) & (ts < max(t1, t2 if r > 0 else t1) + 120.0)
        pull = np.sqrt(np.mean(((y[pre] - model[pre]) / float(d["noise_sigma"][e])) ** 2))
        mres = float(np.mean(y[reg] - model[reg])) / A1
        print(f"[event {e}] {nm} A1={A1:g} r={r:g} dt={dt:g} t1={t1:.2f} t2={t2:.2f} ped={ped:.2f} "
              f"smear_rms={float(d['smear_rms'][e]):.2e} pre-onset rms((y - sum)/sigma)={pull:.2f} "
              f"pulse-region mean(y - sum)/A1={mres:+.4f} max={int(y.max())}")
        ax.plot(tt, ped + prim, color=BLUE, lw=1.4, ls="--", label=r"primary $p + A_1 s(t-t_1)$")
        if r > 0:
            ax.plot(tt, ped + sec, color=ORANGE, lw=1.4, ls="-.", label=r"secondary $p + A_2 s(t-t_2)$")
        ax.plot(tt, ped + prim + sec, color=VERM, lw=1.0, label="sum (no noise, no smearing)")
        ax.plot(ts, y, "o", color="black", ms=2.6, label="generated samples")
        title = f"{nm}, $A_1$ = {A1:g} ADC, " + ("control ($r$ = 0)" if r == 0 else f"$r$ = {r:g}, $\\Delta t$ = {dt:+g} ns")
        ax.set_title(title, fontsize=10)
        ax.set_xlim(0, 400)
        ax.grid(alpha=0.25, lw=0.5)
    for ax in axs[-1]:
        ax.set_xlabel("time from sample 0 (ns)")
    for ax in axs[:, 0]:
        ax.set_ylabel("ADC counts")
    h, lab = axs.flat[1].get_legend_handles_labels()
    fig.legend(h, lab, loc="lower center", ncol=4, frameon=False, bbox_to_anchor=(0.5, 0.0))
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    OUT.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT)
    print(f"[done] {OUT}")


if __name__ == "__main__":
    main()
