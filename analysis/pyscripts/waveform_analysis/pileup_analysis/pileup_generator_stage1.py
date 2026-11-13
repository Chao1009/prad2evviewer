#!/usr/bin/env python3
"""pileup_generator_stage1.py -- fully synthetic two-pulse PbWO4 waveforms (stage 1).

    python3 waveform_analysis/pileup_analysis/pileup_generator_stage1.py [--n-events N] [out.npz]

Each event is one FADC channel, 100 samples of 4 ns (the real HyCal readout window):

    y(t) = ped + A1 * s_m(t - t1) + A2 * s_m(t - t2) + noise,  rounded; never above 4095 (see below)

s_m is module m's own unit-peak spline template from the 2.24 GeV H2 template bank
(read by Bank below, the same shapes the matched filter uses; s = 0 for
t <= 0, i.e. t1 and t2 are template onsets in the bank convention, which is what
MFReconstructor reports as p_t0).  t1 = T1_NS + U[0, 4) ns (random sub-sample phase),
t2 = t1 + dt (dt < 0: the secondary arrives first; the main pulse is always the A1 one at t1);
with one template for both pulses dt is both the onset-to-onset and the peak-to-peak gap.  A2 = r * A1.  The r = 0 control has no secondary (A2 = 0,
t2 = NaN).  Noise: white Gaussian, sigma chosen so that after rounding the sample
std matches the module's measured pre-pulse std (~0.9 ADC); the real noise's lag-1
anticorrelation (~-0.3) is not modelled.

Template smearing (SMEAR_K = 1): the bank's shape_sigma (per-module pulse-to-pulse scatter
of real pulses around the unit-peak template, as a fraction of the pulse peak, measured at the
band amplitude a_band) scaled by sqrt(a_band / A1), the measured 1 / sqrt(A) law of the
fluctuation, is the template's uncertainty.  Per event one smooth random perturbation delta(t) is drawn and BOTH
pulses of that event use s'(t) = s(t) + SMEAR_K * delta(t) (same realisation within an
event, independent across events; 0 for t <= 0, so the onset stays put; the peak is not
renormalised).  delta is the module's own cubic B-spline basis with jittered coefficients,
delta = sum_i B_i(t) eps_i, eps_i ~ N(0, (a * shape_sigma(g_i))^2) with g_i the Greville
abscissa of basis function i, multiplied by sqrt(a_band / A1) of the event, and one scale a per
module chosen so that the rms over two points, the peak (t_pk) and the 50 % leading edge (t50), of
std_events[delta(t)] / (shape_sigma(t) sqrt(a_band / A1)) is exactly 1 for every A1 (the ratio at each
point, analytic and from 2000 random draws, and the achieved and target rms per A1 are printed and
stored in meta).  The bank carries no covariance, so the correlation of delta
is the spline's (knot spacing 1 ns on the rise, up to 16 ns in the tail), not a measured
one.  For speed delta is tabulated per event on a 0.2 ns grid (0-250 ns after onset) and
linearly interpolated at the sample times.  The eps draw follows the noise draw of each
chunk.

Unlike pileup_generator.py (which injects one template pulse on a real EVIO base
pulse), nothing here reads data: both pulses are templates, so all truth is exact.

No clipping: combinations whose noiseless unsmeared waveform peaks above 4095 - 4 noise
sigma are dropped from the grid (printed at startup, stored in meta dropped_combos), and any
event with a sample >= 4095 (rare: noise, smearing) is redrawn (noise and smear, same phase).

Events cycle through the (module, A1, (r, dt)) combinations, so every combination
gets N // n_comb or N // n_comb + 1 events ((r, dt) varies fastest, then A1, then
module, so N < 1065 covers only the first modules).  Fixed seed; phase and noise per event.

Output npz (np.load(path, allow_pickle=False)):
  waveforms  uint16  [N, 100]   raw ADC (pedestal included), same dtype/length as
                                pileup_generator.py and the real waveforms
  module     <U4     [N]        module name (row in MODULES = mf bank row)
  a1, a2     float32 [N]        peak heights above pedestal (ADC)
  ratio      float32 [N]        a2 / a1 (0 = control)
  dt_ns      float32 [N]        t2 - t1 (0 for the control)
  t1_ns      float32 [N]        main-pulse onset (ns from sample 0)
  t2_ns      float32 [N]        secondary onset (NaN for the control)
  pedestal   float32 [N]        injected pedestal mean (ADC)
  noise_sigma float32 [N]       injected noise sigma before rounding (ADC)
  n_clipped  int16   [N]        samples at 4095 (always 0: kept for readers)
  smear_k    float32 ()         SMEAR_K used
  smear_rms  float32 [N]        rms over the pulse region of SMEAR_K * delta (unit-peak
                                units), the event's template error
  meta       <U..    ()         JSON: bank path (relative to analysis/pyscripts), sha1, grids, modules, seed, units
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent.parent                             # analysis/pyscripts

BANK = PYDIR / "waveform_analysis" / "pulse_templates" / "2p24GeV_H2" / "templates_2p24GeV_H2.npz"
DEF_OUT = PYDIR / "output" / "pileup_stage1" / "pileup_stage1.npz"
LUT_DT = 0.05                 # ns, template lookup step

# Second ring around the 2x2 beam hole (max(|x|,|y|) = 2.5 pitches, 20 modules; the
# 12-module ring touching the hole is excluded), all with their own template
# (valid 1, fallback 0).  The 8 nearest to the hole (r = 52.9 mm) are two per side;
# take one per side (top, right, bottom, left) plus the second top one.
# (name, x mm, y mm, pedestal ADC, noise sigma ADC), from samples 0-29 of the module's
# pulses of the bank read (pulse_templates/2p24GeV_H2: pools of e_rough >= 0.5 E_beam pulses, <= 1000 per module):
# pedestal = mean; sigma = sqrt(std^2 - 1/12), std = mean per-pulse sample std (0.86-0.97).
MODULES = [
    ("W493", -10.38, 51.88, 170.18, 0.84),
    ("W564", 51.92, 10.38, 170.21, 0.92),
    ("W664", 10.38, -51.88, 183.40, 0.89),
    ("W593", -51.92, -10.38, 145.14, 0.81),
    ("W494", 10.38, 51.88, 173.81, 0.88),
]
A1_ADC = [300.0, 1000.0, 3000.0]
# pileup_generator.py's grid (ratio 0.1-2.0, |dt| 8-128 ns) plus ratio 0.05; dt < 0: the secondary arrives first
# (for dt <= -32 it starts inside the 30-sample pedestal window: intended, the pedestal estimate is then biased)
RATIOS = [0.05, 0.1, 0.2, 0.3, 0.5, 1.0, 2.0]
DT_NS = [-128.0, -96.0, -64.0, -48.0, -32.0, -24.0, -20.0, -16.0, -12.0, -8.0,
         8.0, 12.0, 16.0, 20.0, 24.0, 32.0, 48.0, 64.0, 96.0, 128.0]
NSAMP = 100                   # samples per waveform
CLK_NS = 4.0                  # ns per sample
T1_NS = 152.0                 # main onset before the phase (real median mf t0 ~153 ns)
OVERFLOW = 4095
SEED = 12345
DEF_N = 9900                  # 10 events per combination
SMEAR_K = 1.0                 # template smearing strength in units of shape_sigma
SMEAR_REGION = 0.01           # pulse region for smear_rms: s >= this x peak
SMEAR_NCHECK = 2000           # random draws for the numerical std-ratio check
SMEAR_BLOCK = 4096            # events per delta-grid block (block x 1251 float32 = 20 MB)
SMEAR_DT = 0.2                # ns, grid on which delta is tabulated per event, then linearly
#                               interpolated (basis knots are >= 1 ns apart; 4x cheaper than 0.05)
SMEAR_TMAX = 250.0            # ns after onset; window end 400 ns - earliest onset 152 ns = 248
CHUNK = 100000              # events per block in generate(); float64 work is CHUNK x 100 x 8 B


def combinations(bank):
    """(kept, dropped) lists of (module row, a1, ratio, dt): per (module, a1), the r = 0 control then
    r x dt.  Dropped: noiseless unsmeared peak ped + a1 s + r a1 s(shifted) above OVERFLOW - 4 sigma."""
    pairs = [(0.0, 0.0)] + [(r, d) for r in RATIOS for d in DT_NS]
    t = np.arange(0.0, NSAMP * CLK_NS, 0.1)
    kept, dropped = [], []
    for m in range(len(MODULES)):
        s1 = bank.shape(m, t - T1_NS, dtype=np.float64)
        for a1 in A1_ADC:
            for r, d in pairs:
                y = a1 * s1 + (r * a1 * bank.shape(m, t - T1_NS - d, dtype=np.float64) if r > 0 else 0.0)
                top = MODULES[m][3] + np.max(y)
                (dropped if top > OVERFLOW - 4 * MODULES[m][4] else kept).append((m, a1, r, d))
    return kept, dropped


class Bank:
    """Unit-peak template of each MODULES entry: B-spline tabulated on the LUT_DT grid (float32),
    linearly interpolated by shape(); t_pk = peak time, t50 = 50 % leading-edge time (ns)."""

    def __init__(self):
        from scipy.interpolate import BSpline
        z = np.load(BANK, allow_pickle=False)
        k, self.T_end = int(z["spline_k"]), float(z["T_end"])
        self.n_lut = int(round(self.T_end / LUT_DT)) + 1
        tg = LUT_DT * np.arange(self.n_lut)
        row = {str(s): i for i, s in enumerate(z["names"])}
        self.names = [m[0] for m in MODULES]
        lut = []
        for nm in self.names:
            r = row.get(nm)
            if r is None or z["valid"][r] <= 0 or z["fallback"][r] != 0:
                raise SystemExit(f"[ERROR] {nm} has no own template in {BANK}")
            t = z["spline_t"][r, :z["spline_nt"][r]]
            s = BSpline(t, z["spline_c"][r, :z["spline_nc"][r]], k)
            lut.append(s(np.minimum(tg, t[-k - 1])))         # hold the end value beyond the base interval
        self.lut = np.asarray(lut, dtype=np.float64).astype(np.float32)
        self.t_pk = np.argmax(self.lut, axis=1) * LUT_DT
        self._flat = self.lut.reshape(-1)
        mi = np.arange(len(self.names))
        g = np.linspace(0.0, 1.0, 4001)[None, :] * self.t_pk[:, None]
        sg = self.shape(mi[:, None], g, np.float64)
        i1 = np.argmax(sg >= 0.5, axis=1)
        i0 = np.maximum(i1 - 1, 0)
        s0, s1 = sg[mi, i0], sg[mi, i1]
        fr = np.where(s1 > s0, (0.5 - s0) / np.maximum(s1 - s0, 1e-300), 0.0)
        self.t50 = g[mi, i0] + fr * (g[mi, i1] - g[mi, i0])

    def shape(self, m, t, dtype=np.float32):
        """Unit-peak shape at module m, t ns after onset (m, t broadcast); 0 outside (0, T_end)."""
        t = np.asarray(t, dtype=dtype)
        x = t * dtype(1.0 / LUT_DT)
        i = np.clip(np.floor(x), 0, self.n_lut - 2)
        f = x - i
        b = np.asarray(m, dtype=np.int64) * self.n_lut + i.astype(np.int64)
        v = self._flat[b].astype(dtype) * (1 - f) + self._flat[b + 1].astype(dtype) * f
        return np.where((t > 0) & (t < self.T_end), v, dtype(0))


def load_smear(bank):
    """Per module: basis B [nco, nc] (float32, SMEAR_DT grid, t = 0..SMEAR_TMAX), coefficient
    scales s [nc] (calibrated), region index array, calibration ratios at peak and t50 (analytic, numeric),
    a_band, and per A1 the achieved and target rms over (peak, t50) of delta (unit-peak units, numeric draws)."""
    from scipy.interpolate import BSpline
    z = np.load(BANK, allow_pickle=False)
    row = {str(nm): r for r, nm in enumerate(z["names"])}
    k = int(z["spline_k"])
    n = bank.n_lut
    tg = LUT_DT * np.arange(n)
    tgrid = np.asarray(z["t_grid"], dtype=np.float64)       # 0.1 ns step from -8 ns (bank grid is 0.05)
    if "a_band" not in z.files:
        raise SystemExit(f"[ERROR] {BANK} has no a_band: rebuild the bank with the current template_extract.py")
    out = []
    for name in bank.names:
        r = row[name]
        t = z["spline_t"][r, :z["spline_nt"][r]]
        nc = int(z["spline_nc"][r])
        B = np.empty((n, nc))
        for i in range(nc):
            c = np.zeros(nc)
            c[i] = 1.0
            B[:, i] = np.nan_to_num(BSpline(t, c, k, extrapolate=False)(tg))
        se_raw = np.asarray(z["shape_sigma"][r], dtype=np.float64)
        good = np.isfinite(se_raw)
        se = np.interp(tg, tgrid[good], se_raw[good], left=np.nan, right=np.nan)  # onto the bank grid
        sh = np.asarray(bank.lut[bank.names.index(name)], dtype=np.float64)
        ok = np.isfinite(se)
        reg = np.flatnonzero(ok & (sh >= SMEAR_REGION * sh.max()) & (tg > 0))
        grev = np.array([t[i + 1:i + k + 1].mean() for i in range(nc)])
        s_i = np.interp(grev, tg[ok], se[ok])               # local shape_sigma at each coefficient
        std1 = np.sqrt((B ** 2 * s_i ** 2).sum(axis=1))      # pointwise std of delta for a = 1
        mi = bank.names.index(name)
        cal = np.round(np.array([bank.t_pk[mi], bank.t50[mi]]) / LUT_DT).astype(int)   # peak, t50 on the LUT grid
        a = 1.0 / np.sqrt(np.mean((std1[cal] / se[cal]) ** 2))   # rms of the ratio over (peak, t50) = 1
        s_i = a * s_i
        eps = np.random.default_rng(SEED + 1).standard_normal((SMEAR_NCHECK, nc)) * s_i
        std_num = (eps @ B.T).std(axis=0, ddof=1)            # separate rng: main stream untouched
        a_band = float(z["a_band"][r])
        fac = {a1: np.sqrt(a_band / a1) for a1 in A1_ADC}   # the 1 / sqrt(A) law: delta scales with sqrt(a_band / A1)
        achieved = {a1: float(np.sqrt(np.mean(((eps @ B.T).std(axis=0, ddof=1)[cal] * f) ** 2))) for a1, f in fac.items()}
        target = {a1: float(np.sqrt(np.mean((se[cal] * f) ** 2))) for a1, f in fac.items()}
        stp = int(round(SMEAR_DT / LUT_DT))
        nco = int(SMEAR_TMAX / SMEAR_DT) + 1
        out.append(dict(B=B[::stp][:nco].astype(np.float32), s=s_i,
                        reg=np.unique(reg[reg < (nco - 1) * stp] // stp), nc=nc, scale=float(a), a_band=a_band, achieved=achieved, target=target,
                        ratio_num=std_num[cal] / se[cal], ratio_ana=a * std1[cal] / se[cal]))
    return out


def delta_at(D, t):
    """delta on the SMEAR_DT grid (D [n, nco], t = 0..SMEAR_TMAX) linearly interpolated at
    t [n, ns] (ns after onset); 0 outside (0, SMEAR_TMAX), like bank.shape outside the template."""
    x = t * (1.0 / SMEAR_DT)
    i = np.clip(np.floor(x), 0, D.shape[1] - 2)
    f = (x - i).astype(np.float32)
    i = i.astype(np.int64)
    v = (np.take_along_axis(D, i, 1) * (1 - f) + np.take_along_axis(D, i + 1, 1) * f)
    return np.where((t > 0) & (t < SMEAR_DT * (D.shape[1] - 1)), v, 0.0)


def generate(n_events, bank, smear, comb):
    idx = np.arange(n_events) % len(comb)
    m, a1, r, dt = (np.array(c)[idx] for c in zip(*comb))
    m = m.astype(np.int64)
    rng = np.random.default_rng(SEED)
    t1 = T1_NS + rng.uniform(0.0, CLK_NS, n_events)
    t2 = np.where(r > 0, t1 + dt, np.nan)
    a2 = r * a1
    ped = np.array([x[3] for x in MODULES])[m]
    sig = np.array([x[4] for x in MODULES])[m]
    t = np.arange(NSAMP) * CLK_NS
    wf = np.empty((n_events, NSAMP), dtype=np.uint16)     # preallocated, filled per chunk
    sm_rms = np.zeros(n_events, dtype=np.float32)
    ncmax = max(q["nc"] for q in smear)
    fac = np.sqrt(np.array([q["a_band"] for q in smear])[m] / a1)     # sqrt(a_band / A1): the 1 / sqrt(A) law

    def draw(ev):
        """Noise + smear realisation of events ev: rounded waveforms [n, NSAMP], smear_rms [n]."""
        m_c, a1_c, a2_c, r_c, t1_c, t2_c = m[ev], a1[ev], a2[ev], r[ev], t1[ev], t2[ev]
        y = ped[ev, None] + sig[ev, None] * rng.standard_normal((ev.size, NSAMP))
        z_c = rng.standard_normal((ev.size, ncmax)) * fac[ev, None]    # eps = s * z sqrt(a_band / A1), drawn after the noise
        sm = np.zeros(ev.size, dtype=np.float32)
        for k in range(len(MODULES)):
            s = m_c == k
            y[s] += a1_c[s, None] * bank.shape(k, t[None, :] - t1_c[s, None], dtype=np.float64)
            s2 = s & (r_c > 0)
            y[s2] += a2_c[s2, None] * bank.shape(k, t[None, :] - t2_c[s2, None], dtype=np.float64)
            q = smear[k]
            idx = np.flatnonzero(s)
            for b in range(0, idx.size, SMEAR_BLOCK):
                e = idx[b:b + SMEAR_BLOCK]
                D = SMEAR_K * ((z_c[e, :q["nc"]] * q["s"]).astype(np.float32) @ q["B"].T)
                sm[e] = np.sqrt((D[:, q["reg"]].astype(np.float64) ** 2).mean(axis=1))
                y[e] += a1_c[e, None] * delta_at(D, t[None, :] - t1_c[e, None])
                two = r_c[e] > 0
                y[e[two]] += a2_c[e[two], None] * delta_at(D[two], t[None, :] - t2_c[e[two], None])
        return np.rint(y), sm

    n_redraw = 0
    for lo in range(0, n_events, CHUNK):
        ev = np.arange(lo, min(lo + CHUNK, n_events))
        y, sm = draw(ev)
        bad = (y >= OVERFLOW).any(axis=1)                     # redraw until no sample reaches 4095
        while bad.any():
            n_redraw += int(bad.sum())
            y[bad], sm[bad] = draw(ev[bad])
            bad = (y >= OVERFLOW).any(axis=1)
        wf[ev], sm_rms[ev] = y.astype(np.uint16), sm
    print(f"[clip] {n_redraw} event draws rejected for a sample >= {OVERFLOW}")
    f32 = lambda v: np.asarray(v, dtype=np.float32)     # noqa: E731
    return dict(waveforms=wf, module=np.array([MODULES[k][0] for k in m]),
                a1=f32(a1), a2=f32(a2), ratio=f32(r), dt_ns=f32(dt), t1_ns=f32(t1),
                t2_ns=f32(t2), pedestal=f32(ped), noise_sigma=f32(sig),
                n_clipped=np.zeros(n_events, dtype=np.int16), smear_k=np.float32(SMEAR_K), smear_rms=sm_rms)


def report(d, nc):
    n = len(d["a1"])
    print(f"\n[grid] modules {[x[0] for x in MODULES]}")
    print(f"[grid] A1 (ADC) {A1_ADC}")
    print(f"[grid] ratio r = A2/A1 {RATIOS}  (+ r = 0 control, no secondary)")
    print(f"[grid] dt (ns) {DT_NS}")
    print(f"[grid] {nc} combinations = {len(MODULES)} modules x {len(A1_ADC)} A1 x "
          f"({len(RATIOS)} x {len(DT_NS)} + 1), minus dropped; {n} events")
    key = np.char.add(np.char.add(d["module"], "|"), d["a1"].astype(str))
    _, per_ma = np.unique(key, return_counts=True)
    print(f"[counts] per (module, A1): min {per_ma.min()} max {per_ma.max()}")
    print("[counts] events per (r, dt), summed over the "
          f"{len(MODULES) * len(A1_ADC)} (module, A1); per combination = this / "
          f"{len(MODULES) * len(A1_ADC)}")
    print("   r \\ dt " + "".join(f"{x:>6g}" for x in DT_NS))
    print(f"  {0.0:>6g}  " + f"{int((d['ratio'] == 0).sum()):>6d}  (control)")
    for r in RATIOS:
        row = [int(((d["ratio"] == np.float32(r)) & (d["dt_ns"] == np.float32(x))).sum())
               for x in DT_NS]
        print(f"  {r:>6g}  " + "".join(f"{c:>6d}" for c in row))
    print(f"[clip] max sample {int(d['waveforms'].max())} (limit {OVERFLOW})")


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("out", nargs="?", default=str(DEF_OUT), help=f"output npz (default {DEF_OUT})")
    ap.add_argument("--n-events", type=int, default=DEF_N,
                    help=f"total number of events generated (default {DEF_N})")
    a = ap.parse_args()
    if a.n_events <= 0:
        ap.error("--n-events must be > 0")
    bank = Bank()
    smear = load_smear(bank)
    comb, dropped = combinations(bank)
    print(f"[grid] dropped {len(dropped)} combinations (peak > {OVERFLOW} - 4 sigma): "
          "(module, A1, r, dt) " + ", ".join(f"({MODULES[m][0]}, {a:g}, {r:g}, {d:g})" for m, a, r, d in dropped))
    print(f"[grid] {len(comb)} combinations left; {a.n_events // len(comb)} events each "
          f"(--n-events {10000 * len(comb)} gives 10000 per combination)")
    d = generate(a.n_events, bank, smear, comb)
    meta = dict(
        generator="pileup_generator_stage1.py", run_type="2p24GeV_H2", bank=str(BANK.relative_to(PYDIR)),
        bank_sha1=hashlib.sha1(BANK.read_bytes()).hexdigest(),
        modules=[dict(name=x[0], x_mm=x[1], y_mm=x[2], ped=x[3], noise_sigma=x[4]) for x in MODULES],
        a1_adc=A1_ADC, ratios=RATIOS, dt_ns=DT_NS, n_combinations=len(comb),
        max_adc=OVERFLOW,
        dropped_combos=[[MODULES[m][0], a, r, d] for m, a, r, d in dropped],
        n_events=a.n_events, nsamples=NSAMP, clk_ns=CLK_NS, t1_nominal_ns=T1_NS,
        t1_phase="uniform [0, 4) ns", seed=SEED, overflow=OVERFLOW,
        noise="white Gaussian, sigma = noise_sigma, then rint; events with a sample >= 4095 redrawn",
        onset="t1/t2 = template onsets (bank convention, s = 0 for t <= 0) = mf p_t0; "
              "peak at t + t_pk[module], 50% leading edge at t + t50[module]",
        smear=dict(
            k=SMEAR_K, model="s'(t) = s(t) + K*delta(t) for both pulses of an event; delta = "
            "B-spline basis of the module's template with coefficients eps_i ~ N(0, s_i^2), "
            "s_i = a * shape_sigma at the Greville abscissa; delta = 0 for t <= 0; peak not "
            "renormalised",
            source="shape_sigma * sqrt(a_band / A1)",
            a_band={x[0]: q["a_band"] for x, q in zip(MODULES, smear)},
            achieved_rms_by_a1={x[0]: {f"{a1:g}": q["achieved"][a1] for a1 in A1_ADC} for x, q in zip(MODULES, smear)},
            target_rms_by_a1={x[0]: {f"{a1:g}": q["target"][a1] for a1 in A1_ADC} for x, q in zip(MODULES, smear)},
            region=f"s >= {SMEAR_REGION} x peak with finite shape_sigma, t > 0",
            calibration="a per module: rms over (peak, t50) of std_events[delta]/(shape_sigma sqrt(a_band/A1)) = 1",
            scale_a={x[0]: q["scale"] for x, q in zip(MODULES, smear)},
            ratio_peak_analytic={x[0]: float(q["ratio_ana"][0]) for x, q in zip(MODULES, smear)},
            ratio_peak_numeric={x[0]: float(q["ratio_num"][0]) for x, q in zip(MODULES, smear)},
            ratio_t50_analytic={x[0]: float(q["ratio_ana"][1]) for x, q in zip(MODULES, smear)},
            ratio_t50_numeric={x[0]: float(q["ratio_num"][1]) for x, q in zip(MODULES, smear)},
            n_check=SMEAR_NCHECK,
            seed_usage="main rng: per chunk noise draw, then eps z [chunk, nc]; "
                       "calibration check uses its own rng(seed + 1)"),
        t_pk_ns={x[0]: float(bank.t_pk[k]) for k, x in enumerate(MODULES)},
        t50_ns={x[0]: float(bank.t50[k]) for k, x in enumerate(MODULES)})
    out = Path(a.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(out, meta=np.array(json.dumps(meta)), **d)
    report(d, len(comb))
    for x, q in zip(MODULES, smear):
        nz = d["module"] == x[0]
        print(f"[smear] {x[0]} K={SMEAR_K:g} a={q['scale']:.3f} std-ratio peak analytic "
              f"{q['ratio_ana'][0]:.4f} numeric({SMEAR_NCHECK}) {q['ratio_num'][0]:.4f}, t50 analytic "
              f"{q['ratio_ana'][1]:.4f} numeric {q['ratio_num'][1]:.4f}; "
              f"median smear_rms {np.median(d['smear_rms'][nz]):.2e}")
        print("        achieved/target rms over (peak, t50) per A1 (a_band {:.0f} ADC): ".format(q["a_band"]) + "; ".join(
            f"{a1:g}: {100 * q['achieved'][a1]:.3f} % / {100 * q['target'][a1]:.3f} % = {q['achieved'][a1] / q['target'][a1]:.3f}"
            for a1 in A1_ADC))
    print(f"[done] wrote {out} ({a.n_events} events, waveforms {d['waveforms'].shape} uint16)")


if __name__ == "__main__":
    main()
