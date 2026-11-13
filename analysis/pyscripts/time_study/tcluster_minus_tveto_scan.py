#!/usr/bin/env python3
"""Method scan for t_cluster - t_veto (one split per sample; results never pooled across runs).

Loops over  sample x veto x cluster-time variant x veto-time variant x event selection  and, for each combination, fits the
histogram of dt = t_cluster - t_veto:  0.25 ns bins, Gauss + constant within +-5 ns of the (smoothed) mode, sigma_fit +- error,
robust sigma = half the 16-84 % range of the hits within +-5 ns of the mode, sigma_core = iterative +-2 sigma Gaussian without
constant on 0.2 ns bins (the earlier module-timing method, plot_module_timing.fit_core), and background/peak = c/a.

Inputs (read only): benchmark outputs <SCRATCH>/<dir>/<run_split>/{result_clusters.tsv,pulse.npz}; veto dumps
output/module_timing_20261030/dump/<run_split>.npz; database/hycal_time_offsets/time_calib.json; optionally the veto
log-normal-fit pass (--veto-fit-pass, raw EVIO via prad2py, cached in <SCRATCH>/vetofit_<run_split>.npz).

  python3 tcluster_minus_tveto_scan.py --scratch SCRATCH [--sample SLUG ...]            # table + figures
  python3 tcluster_minus_tveto_scan.py --scratch SCRATCH --veto-fit-pass SLUG           # needs prad2py (source setup_env.sh)

NOTE seed_time_ns of the e100 benchmark outputs is RAW (the benchmark ran with time calibration off); the calibrated
variants subtract the per-module offset of time_calib.json (t - offset, like the benchmark's --time-calib).
"""
import argparse
import sys
from pathlib import Path

import numpy as np
from scipy.optimize import curve_fit

HERE = Path(__file__).resolve().parent
PY = HERE.parent
DUMP = PY / "output/module_timing_20261030/dump"
TCAL = PY.parents[1] / "database/hycal_time_offsets/time_calib.json"

# slug -> (benchmark dir, run_split, note)
SAMPLES = {
    "0p73GeV_H2": ("ct_0p73GeV_H2_e100", "024590_00050", ""),
    "2p24GeV_H2": ("ct_2p24GeV_H2_25320_e100", "025320_00050", ""),
    "2p24GeV_H2_25228": ("ct_2p24GeV_H2_e100", "025228_00100", "no coincidence peak"),
    "3p49GeV_H2": ("ct_3p49GeV_H2_e100", "024503_00050", ""),
}
WIN = 5.0          # ns, fit / robust window around the mode
GATE = 10.0        # ns, in-time window of the gated veto variants
HMIN_SEED = 500.0  # ADC, seed amplitude cut
HMEM = 50.0        # ADC, member amplitude for the weighted-mean cluster time


# ---------------------------------------------------------------- loading
def load_offsets(names):
    import json
    tab = {m["name"]: float(m["offset_ns"]) for m in json.load(open(TCAL))["modules"] if m.get("resolved")}
    return np.array([tab.get(str(n), 0.0) for n in names])


class Lookup:
    """Nearest-in-time row of a (event, module) table: key = ((event-ev0)*4096+module)*1000 + t + 100."""
    def __init__(self, ev, mi, t, ev0):
        self.ev0 = ev0
        k = ((ev.astype(np.int64) - ev0) * 4096 + mi) * 1000.0 + t + 100.0
        self.o = np.argsort(k)
        self.k = k[self.o]

    def __call__(self, ev, mi, t, tol):
        q = ((ev.astype(np.int64) - self.ev0) * 4096 + mi) * 1000.0 + t + 100.0
        j = np.clip(np.searchsorted(self.k, q), 1, len(self.k) - 1)
        b = np.where(np.abs(self.k[j - 1] - q) <= np.abs(self.k[j] - q), j - 1, j)
        ok = np.abs(self.k[b] - q) < tol
        return np.where(ok, self.o[b], -1)


def load_sample(scratch, slug):
    """Everything per cluster (legacy rows): times, amplitudes, selection flags; plus the veto pulse table."""
    d, rs, _ = SAMPLES[slug]
    base = Path(scratch) / d / rs
    c = np.genfromtxt(base / "result_clusters.tsv", delimiter="\t", names=True, dtype=None, encoding=None)
    c = c[c["config_name"] == "legacy"]
    p = np.load(base / "pulse.npz")
    names = [str(x) for x in p["module_names"]]
    idx = {n: i for i, n in enumerate(names)}
    off = load_offsets(names)
    ev0 = int(p["event"].min())
    m0, m1 = p["method"] == 0, p["method"] == 1
    L0 = Lookup(p["event"][m0], p["midx"][m0], p["time_ns"][m0], ev0)
    L1 = Lookup(p["event"][m1], p["midx"][m1], p["time_ns"][m1], ev0)
    t0, a0 = p["time_ns"][m0], p["amp"][m0]
    a1 = p["amp"][m1]
    n = len(c)
    S = dict(ev=c["evnum"].astype(np.int64), seed=np.array([m.split(":")[0] for m in c["members"]]), moller=(c["moller_pair"] | c["moller_kin"] | c["moller_rej"]) != 0,
             seed_E=c["seed_energy_MeV"].astype(float))
    smi = np.array([idx[s] for s in S["seed"]])
    S["seed_mi"] = smi
    S["t_mf_raw"] = c["seed_time_ns"].astype(float)
    S["off"] = off[smi]
    # seed WaveAnalyzer CFD time (method-0 row nearest to the MF time, within 10 ns) and MF amplitude
    r0 = L0(S["ev"], smi, S["t_mf_raw"], 10.0)
    S["t_cfd_raw"] = np.where(r0 >= 0, t0[np.maximum(r0, 0)], np.nan)
    r1 = L1(S["ev"], smi, S["t_mf_raw"], 0.01)
    S["seed_amp"] = np.where(r1 >= 0, a1[np.maximum(r1, 0)], np.nan)
    # members: energy-weighted mean of calibrated / raw CFD times of members with MF amplitude >= HMEM
    e_, mi_, tm_, row = [], [], [], []
    for i, s in enumerate(c["members"]):
        for tok in s.split(";"):
            nm, e, t = tok.split(":")
            e_.append(float(e)); mi_.append(idx[nm]); tm_.append(float(t)); row.append(i)
    e_, mi_, tm_, row = np.array(e_), np.array(mi_), np.array(tm_), np.array(row)
    evm = S["ev"][row]
    ra = L1(evm, mi_, tm_, 0.01)
    amp = np.where(ra >= 0, a1[np.maximum(ra, 0)], 0.0)
    rw = L0(evm, mi_, tm_, 10.0)
    tw = np.where(rw >= 0, t0[np.maximum(rw, 0)], np.nan)
    use = (amp >= HMEM) & np.isfinite(tw)
    for lab, tt in (("raw", tw), ("cal", tw - off[mi_]), ("mfraw", tm_), ("mfcal", tm_ - off[mi_])):
        u = use if lab[:2] != "mf" else (amp >= HMEM)
        num = np.bincount(row[u], (e_ * tt)[u], n)
        den = np.bincount(row[u], e_[u], n)
        S["wm_" + lab] = np.where(den > 0, num / np.maximum(den, 1e-9), np.nan)
    # veto pulses
    v = np.load(DUMP / f"{rs}.npz")
    S["veto"] = dict(ev=v["ev_num"][v["v_ev"]].astype(np.int64), ch=v["v_ch"].astype(int), t=v["v_t"].astype(float),
                     h=v["v_h"].astype(float), i=v["v_int"].astype(float))
    tr = dict(zip(v["ev_num"].tolist(), v["trig"].tolist()))
    trg = np.array([tr.get(int(e), 0) for e in S["ev"]], dtype=np.int64)
    S["trig_ok"] = ((trg & 0xf00) != 0) & ((trg & 0x3000000) == 0)
    vf = Path(scratch) / f"vetofit_{rs}.npz"
    if vf.exists():
        z = np.load(vf)
        S["vetofit"] = dict(ev=z["ev"].astype(np.int64), ch=z["ch"].astype(int), t=z["t"], h=z["h"], i=z["i"], algo=z["algo"])
    return S


# ---------------------------------------------------------------- cluster-time variants (one function each)
def cl_mf_cal(S): return S["t_mf_raw"] - S["off"]
def cl_mf_raw(S): return S["t_mf_raw"].copy()
def cl_cfd_cal(S): return S["t_cfd_raw"] - S["off"]
def cl_cfd_raw(S): return S["t_cfd_raw"].copy()
def cl_wmcfd_cal(S): return S["wm_cal"]
def cl_wmcfd_raw(S): return S["wm_raw"]
def cl_wmmf_cal(S): return S["wm_mfcal"]


CL_VARIANTS = {"MF_cal": cl_mf_cal, "MF_raw": cl_mf_raw, "CFD_cal": cl_cfd_cal, "CFD_raw": cl_cfd_raw,
               "wmCFD_cal": cl_wmcfd_cal, "wmCFD_raw": cl_wmcfd_raw, "wmMF_cal": cl_wmmf_cal}


# ---------------------------------------------------------------- veto-time selection
def pick_veto(tab, ev, tcl, chs, hmin=0.0, imin=0.0, gate=None, shift=None, algo2=False):
    """Per cluster row the veto pulse (largest integral among counters `chs` that pass the cuts) -> (t_v, h_v), nan if none.
    gate = centre of the accepted t_cluster - t_veto (+-GATE ns); shift[k] is subtracted from the pulse time (counter alignment)."""
    o = np.argsort(ev)
    r = np.searchsorted(ev[o], tab["ev"])
    r[r == len(ev)] = 0
    ok = (ev[o][r] == tab["ev"]) & np.isin(tab["ch"], chs) & (tab["h"] >= hmin) & (tab["i"] >= imin)
    if algo2:
        ok &= tab["algo"] == 2
    row = o[r]
    t = tab["t"] - (0.0 if shift is None else np.asarray(shift)[tab["ch"] - 1])
    ok &= np.isfinite(tcl[row])
    if gate is not None:
        ok &= np.abs(tcl[row] - t - gate) <= GATE
    ii = np.flatnonzero(ok)
    ii = ii[np.lexsort((tab["i"][ii], row[ii]))]
    last = np.r_[row[ii][1:] != row[ii][:-1], True] if len(ii) else np.zeros(0, bool)
    ii = ii[last]
    tv, hv = np.full(len(ev), np.nan), np.full(len(ev), np.nan)
    tv[row[ii]], hv[row[ii]] = t[ii], tab["h"][ii]
    return tv, hv


def walk_correct(dt, h, centre):
    """In-sample time-walk correction: quadratic in ln(h) fitted to the hits within +-4 ns of `centre`, subtracted (centre kept)."""
    m = np.isfinite(dt) & np.isfinite(h) & (h > 0) & (np.abs(dt - centre) < 4)
    if m.sum() < 50:
        return dt
    for _ in range(3):
        co = np.polyfit(np.log(h[m]), dt[m], 2)
        res = dt - np.polyval(co, np.log(np.where(h > 0, h, 1)))
        m = np.isfinite(res) & (h > 0) & (np.abs(res) < 3 * max(np.std(res[m]), 0.2)) & (np.abs(dt - centre) < 4)
    return res + centre


# ---------------------------------------------------------------- histogram fits
def gauss_c(x, a, mu, s, c): return a * np.exp(-0.5 * ((x - mu) / s) ** 2) + c
def gauss(x, a, mu, s): return a * np.exp(-0.5 * ((x - mu) / s) ** 2)


def find_mode(dt, lo=-300.0, hi=100.0):
    """Mode of dt on 0.25 ns bins smoothed by a 1 ns box (robust against single-bin spikes)."""
    h, e = np.histogram(dt, bins=np.arange(lo, hi + 0.25, 0.25))
    hs = np.convolve(h, np.ones(4), "same")
    return 0.5 * (e[:-1] + e[1:])[np.argmax(hs)]


def fit_row(dt, mode=None):
    """dict(N, Npk, mu, sig, err, rob, core, bp) for the array dt (nan ignored)."""
    dt = dt[np.isfinite(dt)]
    R = dict(N=len(dt), Npk=0, mu=np.nan, sig=np.nan, err=np.nan, rob=np.nan, core=np.nan, bp=np.nan)
    if len(dt) < 30:
        return R
    mode = find_mode(dt) if mode is None else mode
    w = dt[np.abs(dt - mode) <= WIN]
    R["Npk"] = len(w)
    if len(w) >= 10:
        R["rob"] = 0.5 * (np.percentile(w, 84) - np.percentile(w, 16))
    h, e = np.histogram(dt, bins=np.arange(mode - WIN, mode + WIN + 1e-9, 0.25))
    x = 0.5 * (e[:-1] + e[1:])
    try:
        p, cv = curve_fit(gauss_c, x, h, p0=(max(h.max() - np.median(h), 1), mode, 1.0, np.median(h)), sigma=np.sqrt(np.maximum(h, 1.0)),
                          bounds=([0, mode - WIN, 0.05, 0], [np.inf, mode + WIN, WIN, np.inf]), maxfev=2000)
        R.update(mu=p[1], sig=abs(p[2]), err=np.sqrt(cv[2, 2]), bp=p[3] / p[0] if p[0] > 0 else np.nan)
    except Exception:
        pass
    R["core"] = fit_core(dt, mode)
    return R


def fit_core(dt, mode):
    """Iterative +-2 sigma Gaussian (no constant) on 0.2 ns bins, <= 5 refits, start sigma 1 ns."""
    mu, s = mode, 1.0
    for _ in range(5):
        lo, hi = mu - max(2 * s, 0.6), mu + max(2 * s, 0.6)
        h, e = np.histogram(dt, bins=np.arange(np.floor(lo / 0.2) * 0.2, hi + 0.2, 0.2))
        x = 0.5 * (e[:-1] + e[1:])
        if h.sum() < 20 or len(x) < 4:
            return np.nan
        try:
            p, _ = curve_fit(gauss, x, h, p0=(h.max(), mu, s), sigma=np.sqrt(np.maximum(h, 1.0)), maxfev=2000)
        except Exception:
            return np.nan
        if abs(p[1] - mu) < 1e-3 and abs(abs(p[2]) - s) < 1e-3:
            mu, s = p[1], abs(p[2]); break
        mu, s = p[1], max(abs(p[2]), 0.05)
    return s


def old_fit(dt, bin_ns=2.0):
    """The simple script's fit: 2 ns bins, Gauss + constant within +-10 ns of the mode -> sigma."""
    dt = dt[np.isfinite(dt)]
    if len(dt) < 30:
        return np.nan
    h, e = np.histogram(dt, bins=np.arange(-300, 100, bin_ns))
    x = 0.5 * (e[1:] + e[:-1])
    mode = x[np.argmax(h)]
    w = abs(x - mode) <= 10
    bg = np.median(h)
    best = (np.inf, np.nan)
    for s0 in (4.0, 1.5, 0.7):                      # the 2 ns-bin fit is fragile: keep the start with the smallest residual
        try:
            p, _ = curve_fit(gauss_c, x[w], h[w], p0=(h.max() - bg, mode, s0, bg), bounds=([0, -np.inf, 0.2, 0], np.inf))
        except Exception:
            continue
        rss = np.sum((gauss_c(x[w], *p) - h[w]) ** 2)
        if rss < best[0]:
            best = (rss, abs(p[2]))
    return best[1]


def per_module_sigma(dt, mods, nmin=30):
    """Median over seed modules (>= nmin hits within +-WIN of the global mode) of the robust sigma and of the core-fit sigma."""
    ok = np.isfinite(dt)
    mode = find_mode(dt[ok])
    rob, core = [], []
    for m in np.unique(mods[ok]):
        x = dt[ok & (mods == m)]
        loc = x[np.abs(x - mode) < 10]
        if len(loc) < nmin:
            continue
        m0 = np.median(loc)
        w = x[np.abs(x - m0) <= WIN]
        if len(w) < nmin:
            continue
        rob.append(0.5 * (np.percentile(w, 84) - np.percentile(w, 16)))
        c = fit_core(x, m0)
        if np.isfinite(c):
            core.append(c)
    return (np.median(rob) if rob else np.nan, np.median(core) if core else np.nan, len(rob))


# ---------------------------------------------------------------- scan
SELECTIONS = ["e100", "e100+A500", "e100+A500+noMoller", "e100+trig"]


def sel_mask(S, name):
    m = np.ones(len(S["ev"]), bool)
    if "A500" in name:
        m &= S["seed_amp"] >= HMIN_SEED
    if "noMoller" in name:
        m &= ~S["moller"]
    if "trig" in name:
        m &= S["trig_ok"]
    return m


def veto_variants(S, tcl, sel, gate, chs, shift):
    """name -> (t_v, h_v) for one cluster-time variant, selection and counter set (gate = first-pass mode)."""
    T, ev = S["veto"], S["ev"]
    out = {"LI": pick_veto(T, ev, tcl, chs, shift=shift),
           "gate": pick_veto(T, ev, tcl, chs, gate=gate, shift=shift),
           "h50": pick_veto(T, ev, tcl, chs, hmin=50, shift=shift),
           "h100": pick_veto(T, ev, tcl, chs, hmin=100, shift=shift),
           "gate+h50": pick_veto(T, ev, tcl, chs, hmin=50, gate=gate, shift=shift),
           "gate+h100": pick_veto(T, ev, tcl, chs, hmin=100, gate=gate, shift=shift),
           "gate+int2000": pick_veto(T, ev, tcl, chs, imin=2000, gate=gate, shift=shift)}
    if "vetofit" in S:
        out["gate+lnfit"] = pick_veto(S["vetofit"], ev, tcl, chs, gate=gate, shift=shift, algo2=True)
    return out


def scan_sample(S, slug):
    rows, cache = [], {}
    n = len(S["ev"])
    for sel in SELECTIONS:
        msel = sel_mask(S, sel)
        for cl, fcl in CL_VARIANTS.items():
            tcl = np.where(msel, fcl(S), np.nan)
            mus = {}
            for v in ("V1", "V2", "V3", "V4", "Vall_noalign", "Vall"):
                if v in ("V1", "V2", "V3", "V4"):
                    chs, shift = [int(v[1])], None
                else:
                    chs = [1, 2, 3, 4]
                    shift = None if v == "Vall_noalign" else np.array([mus[k] - np.mean(list(mus.values())) for k in (1, 2, 3, 4)])
                tv0, _ = pick_veto(S["veto"], S["ev"], tcl, chs, shift=shift)
                dt0 = tcl - tv0
                gate = find_mode(dt0[np.isfinite(dt0)]) if np.isfinite(dt0).sum() > 30 else 0.0
                for vt, (tv, hv) in veto_variants(S, tcl, sel, gate, chs, shift).items():
                    dt = tcl - tv
                    R = fit_row(dt)
                    rows.append(dict(sample=slug, veto=v, cl=cl, vt=vt, sel=sel, **R))
                    cache[(sel, cl, v, vt)] = dt
                    if vt == "gate" and v in ("V1", "V2", "V3", "V4") and np.isfinite(R["mu"]):
                        mus[int(v[1])] = R["mu"]
                    if vt in ("gate", "gate+int2000", "gate+lnfit") and R["Npk"] > 100:
                        dtw = walk_correct(dt, hv, R["mu"])
                        rows.append(dict(sample=slug, veto=v, cl=cl, vt=vt + "+walk(in-sample)", sel=sel, **fit_row(dtw)))
                if v in ("V1", "V2", "V3", "V4") and int(v[1]) not in mus:
                    mus[int(v[1])] = gate
    return rows, cache


def mark_best(rows):
    best = {}
    for r in rows:
        if r["Npk"] >= 100 and np.isfinite(r["sig"]) and r["err"] < 0.3 * r["sig"] and not r["veto"].endswith("noalign"):
            k = (r["sample"], r["veto"])
            if k not in best or r["sig"] < best[k]["sig"]:
                best[k] = r
    nw = {}
    for r in rows:
        if r["Npk"] >= 100 and np.isfinite(r["sig"]) and r["err"] < 0.3 * r["sig"] and not r["veto"].endswith("noalign") and "walk" not in r["vt"]:
            k = (r["sample"], r["veto"])
            if k not in nw or r["sig"] < nw[k]["sig"]:
                nw[k] = r
    for r in rows:
        k = (r["sample"], r["veto"])
        r["best"] = "**BEST**" if best.get(k) is r else ("best (no walk)" if nw.get(k) is r else "")
    return best, nw


def write_table(rows, path, header):
    f = lambda x, d=2: "" if not np.isfinite(x) else f"{x:.{d}f}"
    L = header + ["", "| sample | veto | cluster time | veto time | selection | N | N(+-5ns) | mu (ns) | sigma_fit (ns) | err | robust sigma | core sigma (0.2 ns) | bg/peak | best |",
                  "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    srt = sorted(rows, key=lambda r: (r["sample"], r["veto"], not np.isfinite(r["sig"]), r["sig"] if np.isfinite(r["sig"]) else 0))
    for r in srt:
        L.append(f"| {r['sample']} | {r['veto']} | {r['cl']} | {r['vt']} | {r['sel']} | {r['N']} | {r['Npk']} | {f(r['mu'])} | {f(r['sig'], 3)} | {f(r['err'], 3)} | "
                 f"{f(r['rob'], 3)} | {f(r['core'], 3)} | {f(r['bp'], 3)} | {r['best']} |")
    Path(path).write_text("\n".join(L) + "\n")


# ---------------------------------------------------------------- ingredient ladder (worked example) and per-module median
def ladder(S, slug, veto="V2"):
    """sigma after each ingredient, in the order applied; chs per `veto`."""
    chs = [int(veto[1])] if veto[1:].isdigit() else [1, 2, 3, 4]
    out = []
    msel = sel_mask(S, "e100")
    raw, cal = np.where(msel, cl_mf_raw(S), np.nan), np.where(msel, cl_mf_cal(S), np.nan)

    def dtof(tcl, **kw):
        tv, hv = pick_veto(S["veto"], S["ev"], tcl, chs, **kw)
        return tcl - tv, hv
    d, _ = dtof(raw); out.append(("today: MF raw, largest integral, 2 ns bins, Gauss+const +-10 ns", old_fit(d), len(d[np.isfinite(d)])))
    d, _ = dtof(cal); out.append(("+ time calibration (seed-module offset)", old_fit(d), len(d[np.isfinite(d)])))
    g = find_mode(d[np.isfinite(d)])
    d1, _ = dtof(cal, gate=g); out.append(("+ veto in-time gate +-10 ns (largest integral in the gate)", old_fit(d1), len(d1[np.isfinite(d1)])))
    d, hv = dtof(cal, imin=2000, gate=g); out.append(("+ veto integral > 2000", old_fit(d), len(d[np.isfinite(d)])))
    if len(chs) > 1:
        mus = {}
        for k in (1, 2, 3, 4):
            dk, _ = pick_veto(S["veto"], S["ev"], cal, [k], imin=2000, gate=g)[0], None
            dk = cal - dk
            mus[k] = fit_row(dk)["mu"]
        sh = np.array([mus[k] - np.nanmean(list(mus.values())) for k in (1, 2, 3, 4)])
        d, hv = dtof(cal, imin=2000, gate=g, shift=sh); out.append(("+ per-counter alignment (pooled V1-4)", old_fit(d), len(d[np.isfinite(d)])))
    else:
        out.append(("+ per-counter alignment: no effect on a single counter", old_fit(d), len(d[np.isfinite(d)])))
    R = fit_row(d); out.append(("+ 0.25 ns bins, Gauss+const +-5 ns", R["sig"], R["Npk"]))
    out.append(("  same, core fit 0.2 ns bins (+-2 sigma Gauss, no const)", R["core"], R["Npk"]))
    out.append(("  same, robust sigma (+-5 ns)", R["rob"], R["Npk"]))
    return out


def permodule(S, veto="V2"):
    chs = [int(veto[1])] if veto[1:].isdigit() else [1, 2, 3, 4]
    res = []
    for lab, fcl in (("MF_raw", cl_mf_raw), ("MF_cal", cl_mf_cal)):
        tcl = np.where(sel_mask(S, "e100"), fcl(S), np.nan)
        tv, _ = pick_veto(S["veto"], S["ev"], tcl, chs)
        g = find_mode((tcl - tv)[np.isfinite(tcl - tv)])
        tv, _ = pick_veto(S["veto"], S["ev"], tcl, chs, imin=2000, gate=g)
        res.append((lab, per_module_sigma(tcl - tv, S["seed_mi"])))
    return res


# ---------------------------------------------------------------- figure
def figure(S, slug, combo, out, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    sel, cl, vt = combo
    msel = sel_mask(S, sel)
    tcl = np.where(msel, CL_VARIANTS[cl](S), np.nan)
    fig, axs = plt.subplots(2, 2, figsize=(10, 7))
    for k, ax in zip((1, 2, 3, 4), axs.flat):
        tv0, _ = pick_veto(S["veto"], S["ev"], tcl, [k])
        d0 = tcl - tv0
        g = find_mode(d0[np.isfinite(d0)])
        tv, hv = veto_variants(S, tcl, sel, g, [k], None)[vt.replace("+walk(in-sample)", "")]
        dt = tcl - tv
        if vt.endswith("+walk(in-sample)"):
            dt = walk_correct(dt, hv, fit_row(dt)["mu"])
        R = fit_row(dt)
        dtf = dt[np.isfinite(dt)]
        mode = find_mode(dtf)
        h, e = np.histogram(dtf, bins=np.arange(mode - 15, mode + 15, 0.5))
        ax.stairs(h, e, fill=True, color="C0")
        if np.isfinite(R["mu"]):
            xs = np.linspace(R["mu"] - WIN, R["mu"] + WIN, 200)
            hh, ee = np.histogram(dtf, bins=np.arange(mode - WIN, mode + WIN + 1e-9, 0.25))
            c = np.median(hh) if R["bp"] != R["bp"] else R["bp"]
            a = hh.max()  # amplitude refit for drawing on 0.5 ns bins
            p = curve_fit(gauss_c, 0.5 * (ee[1:] + ee[:-1]), hh, p0=(a, R["mu"], R["sig"], max(np.median(hh), 0)),
                          sigma=np.sqrt(np.maximum(hh, 1.0)), maxfev=2000)[0]
            ax.plot(xs, 2 * gauss_c(xs, *p), "r-")   # 0.5 ns bins hold twice the counts of 0.25 ns bins
        ax.set_title(f"V{k}")
        ax.set_xlabel("t_cluster - t_veto (ns)"); ax.set_ylabel("events / 0.5 ns")
        ax.text(0.03, 0.95, f"N = {R['N']}\n$\\mu$ = {R['mu']:.2f} ns\n$\\sigma$ = {R['sig']:.2f} +- {R['err']:.2f} ns\nrobust {R['rob']:.2f} ns",
                transform=ax.transAxes, va="top")
    fig.suptitle(f"{label}: {cl}, {vt}, {sel}  (0.5 ns bins; fit 0.25 ns, Gauss+const +-5 ns)")
    fig.tight_layout(); fig.savefig(out, dpi=130); plt.close(fig)


# ---------------------------------------------------------------- veto log-normal-fit pass (raw EVIO)
def veto_fit_pass(scratch, slug):
    """WaveAnalyzer.analyze_fit on the veto channels (crate 3 slot 3 ch 0-3) of the events of the sample's cluster file."""
    sys.path.insert(0, str(PY)); sys.path.insert(0, str(PY / "crosstalk_study"))
    import _common as C
    import crosstalk_study as X
    from prad2py import dec
    d, rs, _ = SAMPLES[slug]
    run, split = int(rs[:6]), int(rs[7:])
    evs = np.genfromtxt(Path(scratch) / d / rs / "result_clusters.tsv", delimiter="\t", names=True, dtype=None, encoding=None,
                        usecols=("evnum",))["evnum"]
    want = set(int(e) for e in evs)
    evio, calib = X.evio_of(run, split)
    p = C.setup_pipeline(evio_path=evio, hc_calib_file=calib)
    ch = dec.EvChannel(); ch.set_config(p.cfg)
    if ch.open_auto(evio) != dec.Status.success:
        sys.exit(f"cannot open {evio}")
    rows = []
    while ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            dd = ch.decode_event(i, with_ssp=False)
            if not dd["ok"]:
                continue
            e = dd["event"]
            evn = int(e.info.event_number)
            if evn not in want:
                continue
            for ri in range(e.nrocs):
                roc = e.roc(ri)
                if not roc.present or p.crate_map.get(roc.tag) != 3:
                    continue
                for c in range(4):
                    cd = roc.slot(3).channel(c)
                    if cd.nsamples <= 0:
                        continue
                    _, _, pks, fits = p.wave_ana.analyze_fit(cd.samples)
                    for pk, ft in zip(pks, fits):
                        rows.append((evn, c + 1, pk.time, pk.time_algo, pk.height, pk.integral,
                                     ft.chi2_per_dof if ft.ok else np.nan))
    A = np.array(rows, dtype=float).reshape(-1, 7)
    np.savez(Path(scratch) / f"vetofit_{rs}.npz", ev=A[:, 0].astype(np.int64), ch=A[:, 1].astype(int), t=A[:, 2], algo=A[:, 3].astype(int),
             h=A[:, 4], i=A[:, 5], chi2=A[:, 6])
    print(f"{slug}: {len(A)} veto pulses, time_algo counts {np.bincount(A[:, 3].astype(int), minlength=3)}")


# ---------------------------------------------------------------- main
def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--scratch", required=True)
    ap.add_argument("--sample", nargs="*", default=list(SAMPLES))
    ap.add_argument("--veto-fit-pass", metavar="SLUG")
    ap.add_argument("--table", default=str(HERE / "tables/method_scan.md"))
    ap.add_argument("--figdir", default=str(HERE / "figures"))
    a = ap.parse_args()
    if a.veto_fit_pass:
        return veto_fit_pass(a.scratch, a.veto_fit_pass)
    allrows, extra, bestcombo = [], [], {}
    for slug in a.sample:
        S = load_sample(a.scratch, slug)
        rows, _ = scan_sample(S, slug)
        allrows += rows
        best, nw = mark_best(rows)
        print(f"== {slug}: {len(rows)} combinations, {len(S['ev'])} clusters")
        for k, r in sorted(best.items()):
            print(f"   {k[1]}: sigma_fit {r['sig']:.3f}+-{r['err']:.3f} robust {r['rob']:.3f}  [{r['cl']} | {r['vt']} | {r['sel']}] Npk={r['Npk']}")
            r = nw[k]
            print(f"      no walk: sigma_fit {r['sig']:.3f}+-{r['err']:.3f} robust {r['rob']:.3f}  [{r['cl']} | {r['vt']} | {r['sel']}] Npk={r['Npk']}")
        for v in ("V1", "V2", "V3", "V4", "Vall"):
            extra.append((slug, v, ladder(S, slug, v)))
        extra.append((slug, "permodule", permodule(S, "Vall")))
        # one combination per sample for the figure: smallest mean sigma over V1-V4
        sc = {}
        for r in rows:
            if r["veto"] in ("V1", "V2", "V3", "V4") and r["Npk"] >= 60 and np.isfinite(r["sig"]) and r["err"] < 0.4 * r["sig"] and "walk" not in r["vt"]:
                sc.setdefault((r["sel"], r["cl"], r["vt"]), []).append(r["sig"])
        sc = {k: np.mean(v) for k, v in sc.items() if len(v) >= 3}
        bestcombo[slug] = min(sc, key=sc.get)
        print("   figure combination:", bestcombo[slug], f"mean sigma {sc[bestcombo[slug]]:.3f}")
        figure(S, slug, bestcombo[slug], Path(a.figdir) / f"tcluster_minus_tveto_{slug}_best.png", slug)
    hdr = ["# t_cluster - t_veto method scan", "",
           "Generated by `tcluster_minus_tveto_scan.py`. One split per sample; fit: 0.25 ns bins, Gauss + constant within +-5 ns of the mode; "
           "robust sigma = half the 16-84 % range of the hits within +-5 ns of the mode (truncated, so it understates widths above ~2 ns); "
           "core sigma = iterative +-2 sigma Gaussian, no constant, 0.2 ns bins. BEST = smallest sigma_fit per (sample, veto) among rows with "
           ">= 100 hits in +-5 ns, err < 30 % of sigma, not the unaligned pooled row. `walk(in-sample)` variants fit the correction on the same "
           "hits, so they are optimistic. Cluster time: MF/CFD = seed matched-filter / WaveAnalyzer CFD time, cal/raw = with/without the "
           "time_calib.json offset, wm = energy-weighted mean over members with amplitude >= 50 ADC. Veto time: LI = largest-integral pulse; "
           "gate = within +-10 ns of the first-pass mode; h50/h100 = pulse height cut; int2000 = integral > 2000; lnfit = log-normal-fit CFD (time_algo 2) "
           "(only if the veto-fit pass was run). Vall = largest integral over V1-V4, with (Vall) or without (Vall_noalign) per-counter alignment."]
    write_table(allrows, a.table, hdr)
    import pickle
    pickle.dump(extra, open(Path(a.scratch) / "scan_extra.pkl", "wb"))
    for slug, v, x in extra:
        print("##", slug, v)
        if v == "permodule":
            for lab, (rob, core, nm) in x:
                print(f"   per-module median ({lab}, gate+int2000, Vall): robust {rob:.3f} core {core:.3f} over {nm} modules")
        else:
            for lab, s, n in x:
                print(f"   {s:.3f}  N={n}  {lab}")


if __name__ == "__main__":
    main()
