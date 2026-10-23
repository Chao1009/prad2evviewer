#!/usr/bin/env python3
"""mf_validate.py -- pulse-level validation of the matched-filter (mf) pulse
reconstruction against the WaveAnalyzer rawsum path (plan_20261022.md, Phase 2).

  python3 mf_validate.py pulses DUMP.npz [DUMP2.npz ...] --out DIR
                         [--time-calib PATH] [--hycal-map PATH]
  python3 mf_validate.py timecal DUMP.npz [...] --method {0,1} --out OUT.json [--ref W495] [--min-entries 100]
                         [--height-min 100] [--ebeam-frac 0.8] [--ebeam MeV] [--hycal-map PATH]
                         (per-module time offsets, port of hycal_module_time_calib.cpp; `pulses` takes
                          --time-calib-mf PATH to apply such a table to mf times in V7)
  python3 mf_validate.py pileup SYNTH_DIR --out DIR [--templates PATH]

`pulses` reads the dump NPZ written by the benchmark with --pulse-dump (one
file per run; both methods on the same events) and computes
  V5  template (h50) vs background shape on isolated waveforms
  V6  E_units(mf)/E_units(rawsum) linearity + integration-range effect
  V7  dt vs neighbour energy (seed/neighbour), per-module t_mf - t_rawsum
  V8  rawsum finding fraction, extra mf pulses, fake rate, saturation
Outputs per run in DIR: pulses_<run>.json/.md, PNGs, optionally
mf_time_offsets_<run>.json; plus pulses_summary.md over all dumps.
"""
import argparse
import json
import math
import os
import sys

import numpy as np
import matplotlib
matplotlib.use(os.environ.get("MPLBACKEND", "Agg"))
import matplotlib.pyplot as plt

DT_SAMPLE = 4.0                      # ns per sample
MATCH_NS = 8.0                       # V8 matching half-window
BG_TAU = (7.38, 38.55)               # background shape (pulse_templates_024177)
E_BINS_V8 = [5, 20, 100, 500, np.inf]
E_BINS_V7 = [0, 20, 100, np.inf]
AMP_BINS = [5, 20, 50, 100, 200, 500, 1000, 2000, np.inf]
SNR_BINS = [0, 5, 10, 20, 50, 100, 200, 500, 1000, np.inf]
V6_EDGES = np.geomspace(5.0, 4000.0, 13)
PARENT_AMP_BINS = [0, 20, 50, 100, 500, 2000, np.inf]


# ----------------------------------------------------------------- helpers
def mad(x):
    return float(np.median(np.abs(x - np.median(x)))) if len(x) else float("nan")


def sig_mad(x):
    return 1.4826 * mad(np.asarray(x))


def med(x):
    return float(np.median(x)) if len(x) else float("nan")


def label_bins(edges):
    out = []
    for lo, hi in zip(edges[:-1], edges[1:]):
        out.append(f">{lo:g}" if np.isinf(hi) else f"{lo:g}-{hi:g}")
    return out


def group_apply(g, v, fn, minn=1):
    """{group id: fn(values)} for groups with >= minn members."""
    if len(g) == 0:
        return {}
    o = np.argsort(g, kind="stable")
    u, st = np.unique(g[o], return_index=True)
    return {int(k): fn(p) for k, p in zip(u, np.split(v[o], st[1:])) if len(p) >= minn}


def last_per_group(key, order_val):
    """Index (into key) of the row with the largest order_val in each key group."""
    o = np.lexsort((order_val, key))
    ks = key[o]
    last = np.r_[ks[1:] != ks[:-1], True]
    return o[last]


def sanitize(o):
    if isinstance(o, dict):
        return {str(k): sanitize(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [sanitize(v) for v in o]
    if isinstance(o, np.ndarray):
        return sanitize(o.tolist())
    if isinstance(o, (np.floating, float)):
        return None if not math.isfinite(float(o)) else float(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.bool_):
        return bool(o)
    return o


def two_tau(t, tr, tf):
    """p(t) = (1-exp(-t/tr)) exp(-t/tf) / T_max for t>0 (unit peak), else 0."""
    u = tr / (tr + tf)
    tmax = (1.0 - u) * u ** (tr / tf)
    tt = np.maximum(t, 0.0)
    return np.where(t > 0, (1.0 - np.exp(-tt / tr)) * np.exp(-tt / tf) / tmax, 0.0)


def db_dir():
    d = os.environ.get("PRAD2_DATABASE_DIR")
    if d:
        return d
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "database")


# ------------------------------------------------------------- input / aux
def load_dump(path):
    z = np.load(path, allow_pickle=False)
    d = {k: z[k] for k in z.files}
    n = len(d["event"])
    nan, f0 = np.full(n, np.nan), np.zeros(n, bool)
    for k, dflt in (("chi2dof", nan), ("sigma_raw", nan), ("ped_rms", nan),
                    ("onset", nan), ("clipped", f0), ("wa_fallback", f0),
                    ("npk", np.ones(n, np.int32))):
        d.setdefault(k, dflt)
    d["event"] = d["event"].astype(np.int64)
    d["midx"] = d["midx"].astype(np.int64)
    d["method"] = d["method"].astype(np.int8)
    d["n"] = n
    run = str(np.asarray(d["run"]).item()) if "run" in d else os.path.splitext(os.path.basename(path))[0]
    d["run"] = "".join(ch for ch in run if ch.isalnum() or ch in "_-") or "run"
    if "module_names" not in d:
        d["module_names"] = np.array([f"M{i}" for i in range(int(d["midx"].max()) + 1)])
    return d


def load_hycal_map(path, names):
    """-> (xy[nmod,2] NaN if unknown, pitch). hycal_map.json = list of
    {'n': name, 't': type, 'geo': {'x','y','sx','sy',...}}; mapped by name."""
    if path is None or not os.path.isfile(path):
        return None, None
    ent = {e["n"]: e for e in json.load(open(path))}
    xy = np.full((len(names), 2), np.nan)
    sx = []
    for i, nm in enumerate(names):
        e = ent.get(str(nm))
        if e is not None and e.get("t") == "PbWO4":
            xy[i] = (e["geo"]["x"], e["geo"]["y"])
            sx.append(e["geo"]["sx"])
    return xy, (float(np.median(sx)) if sx else None)


def load_time_calib(path, names):
    """time_calib.json = {'modules': [{'name','module_id','offset_ns','resolved',..}]}.
    -> (off[nmod] ns, 0 for unresolved/unlisted, raw json). Benchmark convention:
    t_corrected = t_raw - offset (load_time_offsets / B:1188-1194)."""
    if path is None or not os.path.isfile(path):
        return None, None
    js = json.load(open(path))
    by = {e["name"]: e for e in js["modules"]}
    off = np.zeros(len(names))
    for i, nm in enumerate(names):
        e = by.get(str(nm))
        if e is not None and e.get("resolved"):
            off[i] = float(e["offset_ns"])
    return off, js


# ---------------------------------------------------------------- matching
def match_pulses(d):
    """V8 rule: per (event, midx) greedy one-to-one rawsum<->mf by |dt|<=8 ns.
    Returns (ri, fi) = matched row indices of rawsum and mf pulses."""
    key = d["event"] * 4096 + d["midx"]
    t = d["time_ns"]
    r = np.flatnonzero(d["method"] == 0)
    f = np.flatnonzero(d["method"] == 1)
    r = r[np.argsort(key[r], kind="stable")]
    f = f[np.argsort(key[f], kind="stable")]
    kr, sr, cr = np.unique(key[r], return_index=True, return_counts=True)
    kf, sf, cf = np.unique(key[f], return_index=True, return_counts=True)
    _, a, b = np.intersect1d(kr, kf, return_indices=True)
    ri, fi = [], []
    one = (cr[a] == 1) & (cf[b] == 1)
    r1, f1 = r[sr[a[one]]], f[sf[b[one]]]
    ok = np.abs(t[r1] - t[f1]) <= MATCH_NS
    ri.append(r1[ok])
    fi.append(f1[ok])
    lr, lf = [], []
    for ia, ib in zip(a[~one], b[~one]):
        rr = r[sr[ia]:sr[ia] + cr[ia]]
        ff = f[sf[ib]:sf[ib] + cf[ib]]
        dt = np.abs(t[rr][:, None] - t[ff][None, :])
        usedr, usedf = set(), set()
        for k in np.argsort(dt, axis=None):
            i, j = divmod(int(k), len(ff))
            if dt[i, j] > MATCH_NS:
                break
            if i in usedr or j in usedf:
                continue
            usedr.add(i)
            usedf.add(j)
            lr.append(rr[i])
            lf.append(ff[j])
    ri.append(np.array(lr, np.int64))
    fi.append(np.array(lf, np.int64))
    return np.concatenate(ri), np.concatenate(fi)


# ------------------------------------------------------------------- V5
def fit_waveforms(wave, ped, sig, onset, tr, tf, ext_rms=None, a_ref=None):
    """(a, t0) free fit per waveform: t0 grid 0.25 ns over +-8 ns around the
    onset, linear weighted-LS amplitude, window onset-2..onset+20 samples,
    variance sigma^2 + (0.03 a s)^2 ("per-sample floor").
    With ext_rms (+ a_ref, the amplitude of a first pass) the template
    extraction convention is used instead (WaveAnalyzer.cpp:1265-1266):
    constant sigma = max(max(ped_rms, 1), 0.03 A) ADC, chi2 reported with the
    sigma of the final A. Returns dict of arrays."""
    n = len(wave)
    W = 23
    i0 = np.clip(np.round(onset).astype(int) - 2, 0, wave.shape[1] - W)
    idx = i0[:, None] + np.arange(W)[None, :]
    y = np.take_along_axis(wave.astype(float), idx, 1) - ped[:, None]
    ext = ext_rms is not None
    sg0 = np.maximum(ext_rms, 1.0) if ext else np.maximum(sig, 0.5)
    if ext:
        sg0 = np.maximum(sg0, 0.03 * a_ref)
    sg = np.broadcast_to(sg0[:, None], idx.shape)
    ts = idx * DT_SAMPLE
    grid = np.arange(-8.0, 8.0 + 1e-9, 0.25)
    out = {k: np.zeros(n) for k in ("a", "t0", "chi2")}
    res = np.zeros((n, W))
    for c0 in range(0, n, 1500):
        sl = slice(c0, min(n, c0 + 1500))
        t0g = onset[sl, None] * DT_SAMPLE + grid[None, :]
        s = two_tau(ts[sl, None, :] - t0g[:, :, None], tr[sl, None, None], tf[sl, None, None])
        w = 1.0 / sg[sl, None, :] ** 2
        den = (s * s * w).sum(2)
        a = np.where(den > 0, (s * y[sl, None, :] * w).sum(2) / np.where(den > 0, den, 1), 0.0)
        a = np.maximum(a, 0.0)
        mod = a[:, :, None] * s
        chi = ((y[sl, None, :] - mod) ** 2 / (sg[sl, None, :] ** 2 + (0 if ext else (0.03 * mod) ** 2))).sum(2)
        k = chi.argmin(1)
        r = np.arange(len(k))
        out["a"][sl], out["t0"][sl] = a[r, k], t0g[r, k]
        out["chi2"][sl] = chi[r, k] / (W - 2)
        res[sl] = y[sl] - mod[r, k]
    if ext:
        sgf = np.maximum(np.maximum(ext_rms, 1.0), 0.03 * out["a"])
        out["chi2"] = (res ** 2).sum(1) / sgf ** 2 / (W - 2)
    out["res_norm"] = np.where(out["a"][:, None] > 0, res / np.maximum(out["a"][:, None], 1e-9), np.nan)
    out["trel"] = ts - out["t0"][:, None]
    return out


def v5_template(d, names, out_png):
    if "wave" not in d or len(d["wave"]) == 0:
        return {"note": "no waveforms in dump"}
    row = d["wave_row"].astype(int)
    mi = d["midx"][row]
    tr_m, tf_m = d.get("tau_r"), d.get("tau_f")
    if tr_m is None or tf_m is None:
        return {"note": "no tau_r/tau_f in dump"}
    onset = d["onset"][row]
    wave = d["wave"]
    ped = d["wave_ped"].astype(float)
    sg = d["wave_sigma"].astype(float)
    # selection: A >= 50, valid taus and onset, no clipped sample, onset in range
    sel = ((d["amp"][row] >= 50) & np.isfinite(tr_m[mi]) & np.isfinite(tf_m[mi])
           & np.isfinite(onset) & (wave.max(1) < 4095) & (onset > -2) & (onset < wave.shape[1] - 5))
    ii = np.flatnonzero(sel)
    if len(ii) == 0:
        return {"note": "no waveforms pass the selection"}
    args = (wave[ii], ped[ii], sg[ii], onset[ii])
    fh = fit_waveforms(*args, tr_m[mi[ii]].astype(float), tf_m[mi[ii]].astype(float))
    fb = fit_waveforms(*args, np.full(len(ii), BG_TAU[0]), np.full(len(ii), BG_TAU[1]))
    m = mi[ii]
    prms = d["ped_rms"][row][ii].astype(float)
    have_ext = bool(np.isfinite(prms).all())
    if have_ext:
        eh = fit_waveforms(*args, tr_m[mi[ii]].astype(float), tf_m[mi[ii]].astype(float), prms, fh["a"])
        eb = fit_waveforms(*args, np.full(len(ii), BG_TAU[0]), np.full(len(ii), BG_TAU[1]), prms, fb["a"])
    cnt = group_apply(m, fh["chi2"], len)
    for thr in (50, 20):
        mods = [k for k, v in cnt.items() if v >= thr]
        if len(mods) >= 10 or thr == 20:
            break
    mh = group_apply(m, fh["chi2"], np.median)
    mb = group_apply(m, fb["chi2"], np.median)
    mods = [k for k in mods if k in mh]
    h = np.array([mh[k] for k in mods])
    b = np.array([mb[k] for k in mods])
    res = {"convention": "per-sample floor (top-level keys); see extraction_convention",
           "n_waveforms": int(len(ii)), "min_waveforms_per_module": thr,
           "n_modules": len(mods),
           "frac_h50_chi2_le_1.5": float(np.mean(h <= 1.5)) if len(h) else None,
           "frac_h50_lt_bg": float(np.mean(h < b)) if len(h) else None,
           "median_of_module_chi2_h50": med(h), "median_of_module_chi2_bg": med(b),
           "pass_chi2": bool(len(h) and np.mean(h <= 1.5) >= 0.95),
           "pass_h50_lt_bg": bool(len(h) and np.mean(h < b) >= 0.99),
           "per_module": {str(names[k]): [cnt[k], mh[k], mb[k]] for k in mods}}
    if have_ext:
        xh = group_apply(m, eh["chi2"], np.median)
        xb = group_apply(m, eb["chi2"], np.median)
        h2 = np.array([xh[k] for k in mods])
        b2 = np.array([xb[k] for k in mods])
        res["extraction_convention"] = {
            "frac_h50_chi2_le_1.5": float(np.mean(h2 <= 1.5)) if len(h2) else None,
            "frac_h50_lt_bg": float(np.mean(h2 < b2)) if len(h2) else None,
            "median_of_module_chi2_h50": med(h2), "median_of_module_chi2_bg": med(b2),
            "pass_chi2": bool(len(h2) and np.mean(h2 <= 1.5) >= 0.95),
            "pass_h50_lt_bg": bool(len(h2) and np.mean(h2 < b2) >= 0.99),
            "per_module": {str(names[k]): [cnt[k], xh[k], xb[k]] for k in mods}}
    edges = np.arange(-12, 82, 2.0)
    prof = {}
    for lab, f in (("h50", fh), ("bg", fb)):
        x, r = f["trel"].ravel(), f["res_norm"].ravel()
        ok = np.isfinite(r)
        bi = np.digitize(x[ok], edges) - 1
        prof[lab] = [float(np.mean(r[ok][bi == k])) if np.sum(bi == k) > 20 else float("nan")
                     for k in range(len(edges) - 1)]
    res["residual_profile"] = {"t_rel_ns_left_edges": edges[:-1], **prof}
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    c = edges[:-1] + 1
    ax[0].plot(c, prof["h50"], "o-", label="h50 template")
    ax[0].plot(c, prof["bg"], "s--", label="background shape")
    ax[0].axhline(0, color="k", lw=0.5)
    ax[0].set(xlabel="t - t0 [ns]", ylabel="mean (y - a s)/a", title="V5 normalised residual")
    ax[0].legend()
    ax[1].scatter(b, h, s=6)
    lim = [0, max(2.0, float(np.nanmax(np.r_[h, b, 0.0])) * 1.05)]
    ax[1].plot(lim, lim, "k--", lw=0.6)
    ax[1].axhline(1.5, color="r", lw=0.6)
    ax[1].set(xlabel="module median chi2/dof (background shape)",
              ylabel="module median chi2/dof (h50)", title="V5 per module", xlim=lim, ylim=lim)
    ax[2].hist(np.clip(h, 0, 5), bins=40)
    ax[2].axvline(1.5, color="r")
    ax[2].set(xlabel="module median chi2/dof (h50)", ylabel="modules")
    fig.tight_layout()
    fig.savefig(out_png, dpi=110)
    plt.close(fig)
    return res


# ------------------------------------------------------------------- V6
def v6_linearity(d, ri, fi, names, pngs):
    Er = d["E_MeV"][ri]
    ur, uf = d["E_units"][ri], d["E_units"][fi]
    ok = (ur > 0) & (Er > 0)
    ri, fi, Er, ur, uf = ri[ok], fi[ok], Er[ok], ur[ok], uf[ok]
    ratio = uf / ur
    m = d["midx"][ri]
    res = {"n_pairs": int(len(ri)), "kappa_in_dump": float(d["kappa"]) if "kappa" in d else None}
    nb = len(V6_EDGES) - 1
    bi = np.digitize(Er, V6_EDGES) - 1
    rows = []
    for k in range(nb):
        s = bi == k
        mm = group_apply(m[s], ratio[s], np.median, 5)
        v = np.array(list(mm.values()))
        rows.append({"E_lo": V6_EDGES[k], "E_hi": V6_EDGES[k + 1], "n": int(s.sum()),
                     "pooled_median": med(ratio[s]),
                     "p16": float(np.percentile(ratio[s], 16)) if s.any() else float("nan"),
                     "p84": float(np.percentile(ratio[s], 84)) if s.any() else float("nan"),
                     "median_of_module_medians": med(v), "n_modules": len(v)})
    res["bins"] = rows
    s = (Er >= 500) & (Er <= 2000)
    pooled = med(ratio[s])
    mm = group_apply(m[s], ratio[s], np.median, 10)
    v = np.array(list(mm.values()))
    res.update({"n_pairs_0p5_2GeV": int(s.sum()), "pooled_median_0p5_2GeV": pooled,
                "pass_linearity": bool(s.any() and abs(pooled - 1) <= 0.01),
                "kappa_suggested": (1.0 / pooled) if s.any() and abs(pooled - 1) > 0.01 else None,
                "n_modules_spread": len(v), "module_median_MAD": mad(v),
                "module_median_sigma_MAD": sig_mad(v),
                "pass_spread": bool(len(v) and sig_mad(v) <= 0.02),
                "per_module_median_0p5_2GeV": {str(names[k]): x for k, x in mm.items()}})
    # integration-range effect: E_rawsum/E_mf vs mf amplitude and amp/ped_rms
    inv = ur / uf
    amp = d["amp"][fi]
    snr = amp / np.maximum(d["ped_rms"][fi], 1e-9)
    rng = {}
    for name, x, edges in (("by_amp_ADC", amp, AMP_BINS), ("by_amp_over_ped_rms", snr, SNR_BINS)):
        b = np.digitize(x, edges) - 1
        rng[name] = [{"lo": edges[k], "hi": edges[k + 1], "n": int((b == k).sum()),
                      "median_raw_over_mf": med(inv[b == k]),
                      "p16": float(np.percentile(inv[b == k], 16)) if (b == k).any() else float("nan"),
                      "p84": float(np.percentile(inv[b == k], 84)) if (b == k).any() else float("nan")}
                     for k in range(len(edges) - 1)]
    res["integration_range"] = rng
    # plots
    fig, ax = plt.subplots(figsize=(6.5, 4.5))
    c = np.sqrt(V6_EDGES[:-1] * V6_EDGES[1:])
    ax.fill_between(c, [r["p16"] for r in rows], [r["p84"] for r in rows], alpha=0.25, label="pooled 16-84%")
    ax.plot(c, [r["pooled_median"] for r in rows], "o-", label="pooled median")
    ax.plot(c, [r["median_of_module_medians"] for r in rows], "s--", label="median of module medians")
    ax.axhspan(0.99, 1.01, color="g", alpha=0.15)
    ax.set_xscale("log")
    ax.set(xlabel="E_rawsum [MeV]", ylabel="E_units(mf) / E_units(rawsum)", title="V6 linearity")
    ax.legend()
    fig.tight_layout()
    fig.savefig(pngs[0], dpi=110)
    plt.close(fig)
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    for a, (name, edges, xl) in zip(ax, (("by_amp_ADC", AMP_BINS, "mf amplitude [ADC]"),
                                          ("by_amp_over_ped_rms", SNR_BINS, "mf amplitude / ped_rms"))):
        r = rng[name]
        x = np.arange(len(r))
        a.errorbar(x, [q["median_raw_over_mf"] for q in r],
                   yerr=[[q["median_raw_over_mf"] - q["p16"] for q in r],
                         [q["p84"] - q["median_raw_over_mf"] for q in r]], fmt="o-", capsize=3)
        a.set_xticks(x)
        a.set_xticklabels(label_bins(edges), rotation=45)
        a.axhline(1, color="k", lw=0.5)
        a.set(xlabel=xl, ylabel="E_units(rawsum) / E_units(mf)", title="integration-range effect (median, 16-84%)")
    fig.tight_layout()
    fig.savefig(pngs[1], dpi=110)
    plt.close(fig)
    return res


# ------------------------------------------------------------------- V7
def v7_dt(d, ri, fi, names, xy, pitch, off, run, out_dir, png, off_mf=None):
    res = {}
    if xy is not None and pitch:
        pos = np.nan_to_num(xy, nan=1e9)
        dist = np.hypot(pos[:, None, 0] - pos[None, :, 0], pos[:, None, 1] - pos[None, :, 1])
        NB = (dist <= 1.5 * pitch) & ~np.eye(len(pos), dtype=bool)
        data = {}
        for meth, lab in ((0, "rawsum"), (1, "mf")):
            idx = np.flatnonzero(d["method"] == meth)
            key = d["event"][idx] * 4096 + d["midx"][idx]
            idx = idx[last_per_group(key, d["E_MeV"][idx])]
            ev, mi, E, t = d["event"][idx], d["midx"][idx], d["E_MeV"][idx], d["time_ns"][idx]
            o = np.lexsort((E, ev))
            ev, mi, E, t = ev[o], mi[o], E[o], t[o]
            _, st = np.unique(ev, return_index=True)
            en = np.r_[st[1:], len(ev)] - 1
            eid = np.repeat(np.arange(len(st)), en + 1 - st)
            sp = en[eid]
            keep = (E[sp] > 100) & (np.arange(len(ev)) != sp) & NB[mi[sp], mi]
            tc = t - (off[mi] if off is not None else 0.0)
            o2 = off_mf if (lab == "mf" and off_mf is not None) else off
            tm = t - (o2[mi] if o2 is not None else 0.0)
            data[lab] = (E[keep], t[keep] - t[sp][keep], tc[keep] - tc[sp][keep], int((E[en] > 100).sum()),
                         tm[keep] - tm[sp][keep])
        tab = {}
        for lab, (E, dr, dc, nseed, dm) in data.items():
            rows = []
            for lo, hi in zip(E_BINS_V7[:-1], E_BINS_V7[1:]):
                s = (E >= lo) & (E < hi)
                rows.append({"E_lo": lo, "E_hi": hi, "n": int(s.sum()), "sigma_MAD_raw": sig_mad(dr[s]),
                             "sigma_MAD_corr": sig_mad(dc[s]) if off is not None else None,
                             "sigma_MAD_own_offsets": sig_mad(dm[s]) if off_mf is not None else None})
            rows.append({"E_lo": "all", "E_hi": "", "n": int(len(E)), "sigma_MAD_raw": sig_mad(dr),
                         "sigma_MAD_corr": sig_mad(dc) if off is not None else None,
                         "sigma_MAD_own_offsets": sig_mad(dm) if off_mf is not None else None})
            tab[lab] = {"n_seeds": data[lab][3], "bins": rows}
        res["sigma_MAD"] = tab
        key_ = "sigma_MAD_corr" if off is not None else "sigma_MAD_raw"
        res["mf_le_rawsum_every_bin"] = bool(all(
            (a[key_] is None or b[key_] is None or not np.isfinite(a[key_]) or not np.isfinite(b[key_])
             or a[key_] <= b[key_]) for a, b in zip(tab["mf"]["bins"], tab["rawsum"]["bins"])))
        res["time_calib_used"] = off is not None
        res["mf_offsets_used"] = off_mf is not None
        if off_mf is not None:
            res["mf_own_le_rawsum_every_bin"] = bool(all(
                (np.isfinite(a["sigma_MAD_own_offsets"]) and np.isfinite(b["sigma_MAD_own_offsets"]) and
                 b["sigma_MAD_own_offsets"] < a["sigma_MAD_own_offsets"])
                for a, b in zip(tab["rawsum"]["bins"], tab["mf"]["bins"])))
        fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
        lab_bins = label_bins(E_BINS_V7) + ["all"]
        x = np.arange(len(lab_bins))
        for j, (ky, ttl) in enumerate((("sigma_MAD_raw", "raw times"), ("sigma_MAD_corr", "time_calib offsets"))):
            for lab, sh in (("rawsum", -0.2), ("mf", 0.2)):
                v = [r[ky] if r[ky] is not None else np.nan for r in tab[lab]["bins"]]
                ax[j].bar(x + sh, v, 0.4, label=lab)
            ax[j].set_xticks(x)
            ax[j].set_xticklabels(lab_bins)
            ax[j].set(xlabel="neighbour energy [MeV]", ylabel="sigma_MAD(t_nb - t_seed) [ns]", title=ttl)
            ax[j].legend()
        fig.tight_layout()
        fig.savefig(png[0], dpi=110)
        plt.close(fig)
    else:
        res["note"] = "no hycal_map.json: neighbour analysis skipped"
    # per-module median t_mf - t_rawsum (matched, rawsum amp > 10 ped_rms)
    good = d["amp"][ri] > 10 * d["ped_rms"][ri]
    dt = d["time_ns"][fi] - d["time_ns"][ri]
    mm = group_apply(d["midx"][ri][good], dt[good], np.median, 10)
    v = np.array(list(mm.values()))
    rms = float(np.sqrt(np.mean(v ** 2))) if len(v) else float("nan")
    res.update({"offset_n_modules": len(v), "offset_median_all": med(v), "offset_rms_over_modules": rms,
                "per_module_t_mf_minus_t_rawsum": {str(names[k]): x for k, x in mm.items()}})
    fig, ax = plt.subplots(figsize=(6, 4))
    if len(v):
        ax.hist(v, bins=40)
    ax.set(xlabel="module median t_mf - t_rawsum [ns]", ylabel="modules",
           title=f"V7 per-module offset, RMS {rms:.3f} ns")
    fig.tight_layout()
    fig.savefig(png[1], dpi=110)
    plt.close(fig)
    res["offset_table_written"] = None
    if len(v) and rms > 0.3:
        path = os.path.join(out_dir, f"mf_time_offsets_{run}.json")
        res["offset_table_written"] = path
    return res, {int(k): x for k, x in mm.items()} if (len(v) and rms > 0.3) else None


def write_offset_table(path, names, off_js, deltas):
    """Same schema as time_calib.json; offset_ns = old offset + median(t_mf - t_rawsum)
    so that  t_mf - offset_new  lands on the corrected rawsum time."""
    by_delta = {str(names[k]): x for k, x in deltas.items()}
    if off_js is None:
        js = {"offset_units": "ns", "default": 0.0, "modules": [
            {"name": str(n), "offset_ns": 0.0, "resolved": False} for n in names]}
    else:
        js = json.loads(json.dumps(off_js))
    for e in js["modules"]:
        if e["name"] in by_delta:
            e["offset_ns"] = (float(e.get("offset_ns", 0.0)) if e.get("resolved") else 0.0) + by_delta[e["name"]]
            e["resolved"] = True
    js["note"] = "mf offsets = time_calib offset + per-module median(t_mf - t_rawsum); made by mf_validate.py"
    with open(path, "w") as fh:
        json.dump(js, fh, indent=1)


# ------------------------------------------------------------------- V8
def v8_finding(d, ri, fi, names, xy, pngs):
    res = {}
    r_all = np.flatnonzero(d["method"] == 0)
    f_all = np.flatnonzero(d["method"] == 1)
    matched_r = np.zeros(d["n"], bool)
    matched_f = np.zeros(d["n"], bool)
    matched_r[ri] = True
    matched_f[fi] = True
    Er = d["E_MeV"][r_all]
    ff = matched_r[r_all]
    rows = []
    for lo, hi in zip(E_BINS_V8[:-1], E_BINS_V8[1:]):
        s = (Er >= lo) & (Er < hi)
        rows.append({"E_lo": lo, "E_hi": hi, "n_rawsum": int(s.sum()),
                     "f_found": float(ff[s].mean()) if s.any() else float("nan")})
    res["f_found_pooled"] = rows
    m = d["midx"][r_all]
    per_mod = {}
    for k, (lo, hi) in enumerate(zip(E_BINS_V8[:-1], E_BINS_V8[1:])):
        s = (Er >= lo) & (Er < hi)
        per_mod[label_bins(E_BINS_V8)[k]] = {str(names[a]): b for a, b in
                                            group_apply(m[s], ff[s].astype(float), np.mean, 1).items()}
    res["f_found_per_module"] = per_mod
    hi20 = Er >= 20
    cnt = group_apply(m[hi20], ff[hi20], len)
    fm = group_apply(m[hi20], ff[hi20].astype(float), np.mean)
    for thr in (10, 1):
        mods = [k for k in fm if cnt[k] >= thr]
        if mods:
            break
    f_arr = np.array([fm[k] for k in mods])
    res["f_found_above_20MeV_min_pulses"] = thr
    res["n_modules_above_20MeV"] = len(mods)
    res["pooled_f_found_above_20MeV"] = float(ff[hi20].mean()) if hi20.any() else float("nan")
    if xy is not None:
        inner = [k for k in mods if np.isfinite(xy[k, 0]) and np.hypot(*xy[k]) < 200.0]
        fi_arr = np.array([fm[k] for k in inner])
        res.update({"n_inner_modules": len(inner),
                    "frac_inner_modules_f_found_ge_0p99": float(np.mean(fi_arr >= 0.99)) if len(inner) else None,
                    "pass_inner": bool(len(inner) and np.mean(fi_arr >= 0.99) >= 0.95)})
    else:
        res["note_inner"] = "no hycal_map: inner-ring test skipped"
    # extra (unmatched) mf pulses
    key = d["event"] * 4096 + d["midx"]
    mkey = key[f_all]
    parent = f_all[last_per_group(mkey, d["amp"][f_all])]
    uk, inv = np.unique(mkey, return_inverse=True)
    pk_of_key = np.empty(len(uk), np.int64)
    pk_of_key[np.searchsorted(uk, key[parent])] = parent
    ex = f_all[~matched_f[f_all]]
    par = pk_of_key[np.searchsorted(uk, key[ex])]
    orphan = par == ex
    dtp = d["time_ns"][ex] - d["time_ns"][par]
    Ef = d["E_MeV"][ex]
    pa = d["amp"][par]
    eb = [0] + E_BINS_V8
    ext = {"n_mf": int(len(f_all)), "n_extra": int(len(ex)), "n_extra_orphan_no_rawsum_in_channel": int(orphan.sum()),
           "by_E_MeV": [], "by_parent_amp": []}
    for lo, hi in zip(eb[:-1], eb[1:]):
        s = (Ef >= lo) & (Ef < hi)
        nmf = int(((d["E_MeV"][f_all] >= lo) & (d["E_MeV"][f_all] < hi)).sum())
        ext["by_E_MeV"].append({"E_lo": lo, "E_hi": hi, "n_extra": int(s.sum()), "n_mf": nmf})
    for lo, hi in zip(PARENT_AMP_BINS[:-1], PARENT_AMP_BINS[1:]):
        s = (pa >= lo) & (pa < hi)
        ext["by_parent_amp"].append({"amp_lo": lo, "amp_hi": hi, "n_extra": int(s.sum())})
    nonorph = ~orphan
    band = nonorph & (dtp >= 6) & (dtp <= 18)
    ext["n_with_parent"] = int(nonorph.sum())
    ext["afterpulse_band_6_18ns_n"] = int(band.sum())
    ext["afterpulse_band_6_18ns_frac"] = float(band.sum() / nonorph.sum()) if nonorph.any() else None
    pm = group_apply(d["midx"][ex][band], np.ones(band.sum()), len)
    ext["afterpulse_band_per_module"] = {str(names[k]): v for k, v in pm.items()}
    res["extras"] = ext
    # fake rate: mf onset in samples 0-20 minus rawsum pos in 0-20, per channel-event
    ak, ainv = np.unique(key, return_inverse=True)
    on = d["onset"]
    inwin = np.isfinite(on) & (on >= 0) & (on <= 20)
    nmf = np.bincount(ainv[(d["method"] == 1) & inwin], minlength=len(ak))
    nrw = np.bincount(ainv[(d["method"] == 0) & inwin], minlength=len(ak))
    exc = np.maximum(nmf - nrw, 0).sum()
    res["fake_rate"] = {"excess_mf_pulses": int(exc), "n_channel_events": int(len(ak)),
                        "rate": float(exc / len(ak)) if len(ak) else None,
                        "pass": bool(len(ak) and exc / len(ak) <= 1e-3)}
    cl = d["clipped"][r_all]
    res["saturation"] = {"n_rawsum": int(len(r_all)), "frac_clipped": float(cl.mean()) if len(cl) else None,
                         "n_rawsum_gt_1GeV": int((Er > 1000).sum()),
                         "frac_clipped_gt_1GeV": float(cl[Er > 1000].mean()) if (Er > 1000).any() else None,
                         "mf_wa_fallback_n": int(d["wa_fallback"][f_all].sum()),
                         "mf_wa_fallback_frac": float(d["wa_fallback"][f_all].mean()) if len(f_all) else None}
    # plots
    fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
    x = np.arange(len(rows))
    ax[0].bar(x, [r["f_found"] for r in rows])
    ax[0].set_xticks(x)
    ax[0].set_xticklabels(label_bins(E_BINS_V8))
    ax[0].axhline(0.99, color="r", lw=0.6)
    ax[0].set(xlabel="E_rawsum [MeV]", ylabel="f_found (pooled)", title="V8 finding fraction", ylim=(0.8, 1.01))
    if len(f_arr):
        ax[1].hist(np.clip(f_arr, 0.8, 1.0), bins=40)
    ax[1].axvline(0.99, color="r")
    ax[1].set(xlabel=f"module f_found above 20 MeV (n>={thr})", ylabel="modules")
    fig.tight_layout()
    fig.savefig(pngs[0], dpi=110)
    plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(15, 4))
    h = dtp[nonorph]
    ax[0].hist(h[(h >= -40) & (h <= 80)], bins=np.arange(-40, 82, 2.0))
    ax[0].axvspan(6, 18, color="r", alpha=0.15)
    ax[0].set(xlabel="t_extra - t_parent [ns]", ylabel="extra mf pulses", title="V8 extras vs parent")
    for a, k, lbls, xl in ((ax[1], "by_E_MeV", label_bins(eb), "E_mf [MeV]"),
                           (ax[2], "by_parent_amp", label_bins(PARENT_AMP_BINS), "parent amplitude [ADC]")):
        a.bar(np.arange(len(ext[k])), [r["n_extra"] for r in ext[k]])
        a.set_xticks(np.arange(len(ext[k])))
        a.set_xticklabels(lbls, rotation=30)
        a.set(xlabel=xl, ylabel="extra mf pulses")
    fig.tight_layout()
    fig.savefig(pngs[1], dpi=110)
    plt.close(fig)
    return res


# ---------------------------------------------------------------- driver
def fmt(x, p=3):
    return "n/a" if x is None or (isinstance(x, float) and not math.isfinite(x)) else (
        f"{x:.{p}g}" if isinstance(x, float) else str(x))


def checks_of(R):
    c = []
    v5 = R.get("V5", {})
    for lab, v in (("per-sample floor", v5), ("extraction convention", v5.get("extraction_convention"))):
        if v and "frac_h50_chi2_le_1.5" in v:
            c.append((f"V5 [{lab}] modules with median chi2/dof(h50) <= 1.5 (>=95%)",
                      f"frac {fmt(v['frac_h50_chi2_le_1.5'])}, {v5['n_modules']} modules (>= {v5['min_waveforms_per_module']} wf), "
                      f"median chi2 h50 {fmt(v['median_of_module_chi2_h50'])} / bg {fmt(v['median_of_module_chi2_bg'])}", v["pass_chi2"]))
            c.append((f"V5 [{lab}] modules with h50 < background (>=99%)", fmt(v["frac_h50_lt_bg"]), v["pass_h50_lt_bg"]))
    v6 = R.get("V6", {})
    if "pooled_median_0p5_2GeV" in v6:
        c.append(("V6 pooled E_mf/E_rawsum, 0.5-2 GeV (1 +- 0.01)", fmt(v6["pooled_median_0p5_2GeV"], 4)
                  + (f"; kappa={fmt(v6['kappa_suggested'], 4)}" if v6["kappa_suggested"] else ""), v6["pass_linearity"]))
        c.append(("V6 module spread sigma_MAD (<= 2%)", f"{fmt(v6['module_median_sigma_MAD'])} ({v6['n_modules_spread']} modules)", v6["pass_spread"]))
    v7 = R.get("V7", {})
    if "sigma_MAD" in v7:
        sm = v7["sigma_MAD"]
        rb, mb = sm["rawsum"]["bins"], sm["mf"]["bins"]
        bl = ["<20", "20-100", ">100"]
        if v7["time_calib_used"]:
            c.append(("V7 overall sigma_MAD, rawsum time_calib offsets applied to both [ns] rawsum / mf",
                      f"{fmt(rb[-1]['sigma_MAD_corr'])} / {fmt(mb[-1]['sigma_MAD_corr'])}; per bin " + ", ".join(
                          f"{l} MeV {fmt(x['sigma_MAD_corr'])} / {fmt(y['sigma_MAD_corr'])}" for l, x, y in zip(bl, rb, mb)),
                      v7["mf_le_rawsum_every_bin"]))
        if v7.get("mf_offsets_used"):
            c.append(("V7 overall sigma_MAD, mf with mf offsets vs rawsum with rawsum offsets [ns] rawsum / mf",
                      f"{fmt(rb[-1]['sigma_MAD_own_offsets'])} / {fmt(mb[-1]['sigma_MAD_own_offsets'])}; per bin " + ", ".join(
                          f"{l} MeV {fmt(x['sigma_MAD_own_offsets'])} / {fmt(y['sigma_MAD_own_offsets'])}" for l, x, y in zip(bl, rb, mb)),
                      v7["mf_own_le_rawsum_every_bin"]))
        ok = all(not (np.isfinite(x["sigma_MAD_raw"]) and np.isfinite(y["sigma_MAD_raw"])) or y["sigma_MAD_raw"] <= x["sigma_MAD_raw"]
                 for x, y in zip(rb, mb))
        c.append(("V7 overall sigma_MAD, raw times [ns] rawsum / mf",
                  f"{fmt(rb[-1]['sigma_MAD_raw'])} / {fmt(mb[-1]['sigma_MAD_raw'])}; per bin " + ", ".join(
                      f"{l} MeV {fmt(x['sigma_MAD_raw'])} / {fmt(y['sigma_MAD_raw'])}" for l, x, y in zip(bl, rb, mb)), bool(ok)))
    c.append(("V7 RMS of module median t_mf - t_rawsum (<= 0.3 ns)", fmt(v7.get("offset_rms_over_modules")),
              None if v7.get("offset_rms_over_modules") is None else bool(v7["offset_rms_over_modules"] <= 0.3)))
    v8 = R.get("V8", {})
    if "pass_inner" in v8:
        c.append(("V8 inner modules with f_found(>20 MeV) >= 0.99 (>=95%)",
                  f"frac {fmt(v8['frac_inner_modules_f_found_ge_0p99'])}, {v8['n_inner_modules']} modules", v8["pass_inner"]))
    if "fake_rate" in v8:
        c.append(("V8 fake rate (<= 1e-3)", fmt(v8["fake_rate"]["rate"]), v8["fake_rate"]["pass"]))
        c.append(("V8 afterpulse band +6..+18 ns, fraction of extras", fmt(v8["extras"]["afterpulse_band_6_18ns_frac"]), None))
        c.append(("Saturation: clipped fraction all / >1 GeV; mf WA fallbacks",
                  f"{fmt(v8['saturation']['frac_clipped'])} / {fmt(v8['saturation']['frac_clipped_gt_1GeV'])}; "
                  f"{v8['saturation']['mf_wa_fallback_n']}", None))
    return c


def run_pulses(args):
    os.makedirs(args.out, exist_ok=True)
    db = db_dir()
    map_path = args.hycal_map or os.path.join(db, "hycal_map.json")
    tc_path = args.time_calib or os.path.join(db, "hycal_time_offsets", "time_calib.json")
    summary = ["# mf_validate pulses summary", ""]
    allc = {}
    for path in args.dumps:
        d = load_dump(path)
        run, names = d["run"], d["module_names"]
        P = lambda s: os.path.join(args.out, f"pulses_{run}_{s}.png")
        xy, pitch = load_hycal_map(map_path, names)
        off, off_js = load_time_calib(tc_path, names)
        off_mf, _ = load_time_calib(args.time_calib_mf, names)
        print(f"[{run}] {d['n']} pulses, map={'ok' if xy is not None else 'missing'}, "
              f"time_calib={'ok' if off is not None else 'missing'}", flush=True)
        ri, fi = match_pulses(d)
        R = {"run": run, "dump": path, "n_pulses": d["n"], "n_rawsum": int((d["method"] == 0).sum()),
             "n_mf": int((d["method"] == 1).sum()), "n_matched": int(len(ri))}
        R["V5"] = v5_template(d, names, P("v5"))
        R["V6"] = v6_linearity(d, ri, fi, names, [P("v6_linearity"), P("v6_intrange")])
        R["V7"], deltas = v7_dt(d, ri, fi, names, xy, pitch, off, run, args.out, [P("v7_dt"), P("v7_offsets")], off_mf)
        if deltas:
            write_offset_table(R["V7"]["offset_table_written"], names, off_js, deltas)
        R["V8"] = v8_finding(d, ri, fi, names, xy, [P("v8_found"), P("v8_extras")])
        with open(os.path.join(args.out, f"pulses_{run}.json"), "w") as fh:
            json.dump(sanitize(R), fh, indent=1)
        ch = checks_of(R)
        allc[run] = ch
        lines = [f"## Run {run}", "", f"pulses {R['n_pulses']} (rawsum {R['n_rawsum']}, mf {R['n_mf']}), "
                 f"matched pairs {R['n_matched']}", "", "| check | value | pass |", "|---|---|---|"]
        lines += [f"| {a} | {b} | {'PASS' if p else 'FAIL' if p is not None else '-'} |" for a, b, p in ch]
        lines += ["", "V5 note: `sigma_raw` = 1.4826*MAD on integer ADC samples is quantised to {0.3 floor, 0.74, 1.48, 2.22}, "
                  "which inflates the per-sample-floor chi2 on quiet channels (the extraction convention uses max(ped_rms, 1) instead)."]
        if R["V7"].get("offset_table_written"):
            lines += ["", f"mf time-offset table written: `{R['V7']['offset_table_written']}`"]
        lines += ["", "V6 bins (E_rawsum MeV: pooled median E_mf/E_raw, n): " + "; ".join(
            f"{r['E_lo']:.3g}-{r['E_hi']:.3g}: {fmt(r['pooled_median'], 4)} ({r['n']})" for r in R["V6"]["bins"]), ""]
        lines += ["Integration range (median E_raw/E_mf by mf amplitude): " + "; ".join(
            f"{a}: {fmt(r['median_raw_over_mf'], 4)}" for a, r in zip(label_bins(AMP_BINS), R["V6"]["integration_range"]["by_amp_ADC"])), ""]
        with open(os.path.join(args.out, f"pulses_{run}.md"), "w") as fh:
            fh.write("\n".join(lines) + "\n")
        summary += lines
    summary += ["", "## Combined", "", "| check | " + " | ".join(allc) + " |", "|---|" + "---|" * len(allc)]
    names_all = []
    for ch in allc.values():
        for a, _, _ in ch:
            if a not in names_all:
                names_all.append(a)
    for a in names_all:
        row = []
        for ch in allc.values():
            m = [(b, p) for x, b, p in ch if x == a]
            row.append("-" if not m else f"{m[0][0]} {'PASS' if m[0][1] else 'FAIL' if m[0][1] is not None else ''}".strip())
        summary.append(f"| {a} | " + " | ".join(row) + " |")
    with open(os.path.join(args.out, "pulses_summary.md"), "w") as fh:
        fh.write("\n".join(summary) + "\n")
    print(f"wrote {args.out}/pulses_summary.md")


# ------------------------------------------------------------------- V9
V9_MATCH_NS = 8.0


def _sig68(x):
    """Half-width of the central 68% interval."""
    x = np.asarray(x, dtype=np.float64)
    if len(x) < 2:
        return float("nan")
    lo, hi = np.percentile(x, [15.865, 84.135])
    return float(0.5 * (hi - lo))


def load_wave_cfg(path=None):
    path = path or os.path.join(db_dir(), "daq_config.json")
    with open(path) as fh:
        a = json.load(fh)["fadc250_waveform"]["analyzer"]
    return dict(smooth_order=int(a["smooth_order"]), int_tail_ratio=float(a["int_tail_ratio"]),
                tail_break_n=int(a["tail_break_n"]), ped_nsamples=int(a["ped_nsamples"]),
                overflow=int(a["overflow"]))


def _t50_of(tr, tf):
    """Leading-edge 50% time after onset of the unit model pulse (as TemplateBank.t50)."""
    import mf_pulse
    tm = mf_pulse.tmax_of(tr, tf)
    tpk = tr * math.log1p(tf / tr)
    g = np.linspace(0.0, 1.0, 4001) * tpk
    s = mf_pulse.template(g, tr, tf, tm)
    i1 = int(np.argmax(s >= 0.5))
    i0 = max(i1 - 1, 0)
    fr = (0.5 - s[i0]) / max(s[i1] - s[i0], 1e-300) if s[i1] > s[i0] else 0.0
    return float(g[i0] + fr * (g[i1] - g[i0]))


def _match_two(tm, tt):
    """Greedy one-to-one match of mf times tm to the 2 truth times tt within V9_MATCH_NS.
    Returns the mf pulse index for each truth (-1 when unmatched)."""
    out = [-1, -1]
    if len(tm) == 0:
        return out
    dtm = np.abs(np.asarray(tm)[:, None] - np.asarray(tt)[None, :])
    used = set()
    for k in np.argsort(dtm, axis=None, kind="stable"):
        i, j = divmod(int(k), 2)
        if dtm[i, j] > V9_MATCH_NS:
            break
        if i in used or out[j] >= 0:
            continue
        used.add(i)
        out[j] = i
    return out


def v9_cell(R):
    """R: per-event dicts of one (dT, ratio) cell -> per-cell statistics."""
    n = len(R)
    fi = fb = both = cnt2 = clean = nfp = evfp = 0
    ea, ea_b, et50, eon = [], [], [], []
    for r in R:
        mi, mb, npk = r["m_inj"], r["m_base"], r["npk"]
        fi += mi >= 0
        fb += mb >= 0
        both += (mi >= 0) and (mb >= 0)
        cnt2 += npk >= 2
        fp = npk - ((mi >= 0) + (mb >= 0))
        nfp += fp
        evfp += fp > 0
        clean += (mi >= 0) and (mb >= 0) and fp == 0
        if mi >= 0:
            ea.append(100.0 * (r["amp"][mi] - r["h_inj"]) / r["h_inj"])
            et50.append(r["t50"][mi] - r["t50_inj"])
            eon.append(r["t0"][mi] - r["on_inj"])
        if mb >= 0:
            ea_b.append(100.0 * (r["amp"][mb] - r["h_base"]) / r["h_base"])
    pc = lambda k: 100.0 * k / n if n else float("nan")
    return dict(
        n_events=n, found_inj_pct=pc(fi), found_base_pct=pc(fb), both_found_pct=pc(both),
        count_ge2_pct=pc(cnt2), both_clean_pct=pc(clean),
        fp_per_event=nfp / n if n else float("nan"), fp_event_pct=pc(evfp),
        n_inj_matched=len(ea),
        amp_inj_bias_pct=med(ea), amp_inj_mean_pct=float(np.mean(ea)) if ea else float("nan"),
        amp_inj_w68_pct=_sig68(ea),
        amp_base_bias_pct=med(ea_b), amp_base_w68_pct=_sig68(ea_b),
        t50_inj_bias_ns=med(et50), t50_inj_w68_ns=_sig68(et50),
        onset_inj_bias_ns=med(eon), onset_inj_w68_ns=_sig68(eon))


def _grid(cells, key, dts, ratios):
    g = np.full((len(ratios), len(dts)), np.nan)
    for c in cells:
        g[ratios.index(c["ratio"]), dts.index(c["dt_ns"])] = c.get(key, np.nan)
    return g


def _heat(ax, g, dts, ratios, title, cmap, vmin=None, vmax=None, fmt_s="{:.0f}"):
    im = ax.imshow(g, origin="lower", aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax)
    ax.set_xticks(range(len(dts)))
    ax.set_xticklabels([f"{d:g}" for d in dts], fontsize=7)
    ax.set_yticks(range(len(ratios)))
    ax.set_yticklabels([f"{r:g}" for r in ratios], fontsize=8)
    ax.set_xlabel("dT peak-to-peak (ns)")
    ax.set_title(title, fontsize=9)
    for i in range(len(ratios)):
        for j in range(len(dts)):
            if np.isfinite(g[i, j]):
                ax.text(j, i, fmt_s.format(g[i, j]), ha="center", va="center", fontsize=6)
    return im


def _v9_cells(rec, Z, Y, mid, overflow):
    """Reconstruct all waveforms and return (cells, mf stats). Cells follow the NPZ order."""
    res = rec.run(Y, mid, overflow)
    bounds = np.searchsorted(res.p_ch, np.arange(len(Y) + 1))   # run() sorts by (channel, t0)
    p_t0, p_t50, p_amp = res.p_t0, res.p_time, res.p_amp
    cells, row = [], 0
    for z in Z:
        # truth t50 = truth onset + t50 of the injected (global) template
        t50g = _t50_of(float(z["meta_tau_r_ns"]), float(z["meta_tau_f_ns"]))
        R = []
        for i in range(len(z["waveforms"])):
            a, b = bounds[row], bounds[row + 1]
            ton = z["truth_onset_t0_ns"][i].astype(np.float64)
            t50t = ton + t50g
            mi, mb = _match_two(p_t50[a:b], t50t)
            R.append(dict(npk=int(b - a), m_inj=mi, m_base=mb, amp=p_amp[a:b], t50=p_t50[a:b],
                          t0=p_t0[a:b], h_inj=float(z["truth_heights_adc"][i, 0]),
                          h_base=float(z["truth_heights_adc"][i, 1]),
                          t50_inj=t50t[0], on_inj=ton[0]))
            row += 1
        c = v9_cell(R)
        c.update(dt_ns=float(z["meta_dt_ns"]), ratio=round(float(z["meta_ratio"]), 4),
                 truth_t50_offset_ns=t50g)
        cells.append(c)
    return cells, res.stats


def _v9_compare(cells, sd):
    """Attach the WA / prototype-MF / deconv numbers of the same cells from the existing JSONs."""
    cmp_ = {}
    for fn, tag in (("mf_benchmark.json", "mf"), ("deconv_benchmark.json", "dc")):
        p = os.path.join(sd, fn)
        if os.path.exists(p):
            with open(p) as fh:
                cmp_[tag] = {(c["dt_ns"], round(c["ratio"], 4)): c for c in json.load(fh)["cells"]}
    for c in cells:
        k = (c["dt_ns"], c["ratio"])
        mb, dc = cmp_.get("mf", {}).get(k), cmp_.get("dc", {}).get(k)
        if mb:
            c["wa_pct"] = mb["wa_eff_pct"]
            c["proto_mf_pct"] = mb["mf_eff_pct"]
        if dc:
            c["deconv_ok_pct"] = 100.0 * dc["deconv_ok_eff"]
            c["deconv_inj_bias_pct"] = dc["median_inj_err_pct"]
            c["deconv_inj_absmed_pct"] = dc["median_abs_inj_err_pct"]
    cells.sort(key=lambda c: (c["dt_ns"], c["ratio"]))


def _v9_targets(cells):
    def tgt(name, sel, key, ok):
        sub = [c for c in cells if sel(c)]
        bad = [c for c in sub if not (np.isfinite(c[key]) and ok(c[key]))]
        return dict(name=name, n_cells=len(sub), passed=not bad,
                    failing=[dict(dt_ns=c["dt_ns"], ratio=c["ratio"], value=c[key]) for c in bad])
    T = [
        tgt("both found >= 95% at dT>=24 ns, ratio>=1", lambda c: c["dt_ns"] >= 24 and c["ratio"] >= 1,
            "both_found_pct", lambda v: v >= 95.0),
        tgt("both found >= 90% at dT>=48 ns, ratio>=0.3", lambda c: c["dt_ns"] >= 48 and c["ratio"] >= 0.3,
            "both_found_pct", lambda v: v >= 90.0),
        tgt("|injected amp bias (median)| <= 3% at dT>=48 ns (all ratios)", lambda c: c["dt_ns"] >= 48,
            "amp_inj_bias_pct", lambda v: abs(v) <= 3.0),
        tgt("|injected amp bias (median)| <= 3% at dT>=48 ns, ratio>=0.3 (supplementary)",
            lambda c: c["dt_ns"] >= 48 and c["ratio"] >= 0.3, "amp_inj_bias_pct", lambda v: abs(v) <= 3.0),
    ]
    wa_ref = [dict(dt_ns=c["dt_ns"], ratio=c["ratio"], mf=c["both_found_pct"], wa=c["wa_pct"],
                   mf_minus_wa=c["both_found_pct"] - c["wa_pct"])
              for c in cells if c["dt_ns"] >= 24 and c["ratio"] >= 1 and "wa_pct" in c]
    T[0]["no_loss_vs_wa"] = all(w["mf_minus_wa"] >= 0 for w in wa_ref)
    T[0]["cells_worse_than_wa"] = [w for w in wa_ref if w["mf_minus_wa"] < 0]
    return T


def _v9_write(out_dir, tag, title, sd, ns, conv, stats, cells, T):
    dts = sorted({c["dt_ns"] for c in cells})
    ratios = sorted({c["ratio"] for c in cells})
    out = dict(synth_dir=sd, variant=title, conventions=conv, mf_stats=stats, targets=T, cells=cells)
    with open(os.path.join(out_dir, f"pileup_v9{tag}.json"), "w") as fh:
        json.dump(sanitize(out), fh, indent=1)

    L = [f"# V9 pile-up truth test{' - ' + title if title else ''}", "",
         f"`{sd}`, {sum(c['n_events'] for c in cells)} events in {len(cells)} cells, {ns} samples per waveform.",
         "", "## Conventions", ""]
    L += [f"- **{k}**: {v}" for k, v in conv.items()]
    L += ["", "## Targets", "", "| target | cells | result | failing cells (dT, ratio: value) |", "|---|---|---|---|"]
    for t in T:
        fl = "; ".join(f"{f['dt_ns']:g}, {f['ratio']:g}: {fmt(f['value'], 3)}" for f in t["failing"]) or "-"
        L.append(f"| {t['name']} | {t['n_cells']} | {'PASS' if t['passed'] else 'FAIL'} | {fl} |")
    L += ["", "No loss vs WA (dT>=24, ratio>=1): " + ("yes" if T[0]["no_loss_vs_wa"] else "NO, worse in " + "; ".join(
        f"{w['dt_ns']:g},{w['ratio']:g} (mf {w['mf']:.0f} vs WA {w['wa']:.0f})" for w in T[0]["cells_worse_than_wa"])), ""]
    L += ["## Cells", "",
          "both = both truth pulses matched (%); n>=2 = at least 2 mf pulses (%); clean = both matched and no extra "
          "pulse (%); fp/ev = unmatched mf pulses per event; n_inj = events with the injected pulse matched; "
          "amp bias/w68 = (A_mf-A_true)/A_true in % (median, 68% half-width); t50/onset in ns; WA = WaveAnalyzer "
          ">=2 peaks (%); deconv = deconv ok (%) and its median injected-amplitude error (%).", "",
          "| dT | ratio | both | n>=2 | clean | fp/ev | n_inj | amp bias | amp w68 | t50 bias | t50 w68 | onset bias | WA | deconv | deconv amp |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for c in cells:
        pc0 = lambda k: "n/a" if k not in c else f"{c[k]:.0f}"
        L.append(f"| {c['dt_ns']:g} | {c['ratio']:g} | {c['both_found_pct']:.0f} | {c['count_ge2_pct']:.0f} | "
                 f"{c['both_clean_pct']:.0f} | {c['fp_per_event']:.2f} | {c['n_inj_matched']} | "
                 f"{fmt(c['amp_inj_bias_pct'], 3)} | {fmt(c['amp_inj_w68_pct'], 3)} | {fmt(c['t50_inj_bias_ns'], 2)} | "
                 f"{fmt(c['t50_inj_w68_ns'], 2)} | {fmt(c['onset_inj_bias_ns'], 2)} | {pc0('wa_pct')} | "
                 f"{pc0('deconv_ok_pct')} | {fmt(c.get('deconv_inj_bias_pct'), 3)} |")
    with open(os.path.join(out_dir, f"pileup_v9{tag}.md"), "w") as fh:
        fh.write("\n".join(L) + "\n")

    fig, ax = plt.subplots(1, 3, figsize=(15, 3.8), constrained_layout=True)
    for a, key, ttl in zip(ax, ("both_found_pct", "wa_pct", "deconv_ok_pct"),
                           ("mf: both truth pulses found (%)", "WaveAnalyzer: >=2 peaks (%)", "deconv ok (%)")):
        im = _heat(a, _grid(cells, key, dts, ratios), dts, ratios, ttl, "viridis", 0, 100)
    ax[0].set_ylabel("ratio (inj / base)")
    fig.colorbar(im, ax=ax, shrink=0.8)
    fig.suptitle(title or "primary", fontsize=9)
    fig.savefig(os.path.join(out_dir, f"pileup_v9{tag}_found.png"), dpi=130)
    plt.close(fig)
    fig, ax = plt.subplots(1, 3, figsize=(15, 3.8), constrained_layout=True)
    for a, key, ttl, vm in zip(ax, ("amp_inj_bias_pct", "amp_inj_w68_pct", "t50_inj_bias_ns"),
                               ("mf injected amp bias, median (%)", "mf injected amp 68% half-width (%)",
                                "mf injected t50 bias (ns)"), (10, None, 4)):
        g = _grid(cells, key, dts, ratios)
        if vm:
            im = _heat(a, g, dts, ratios, ttl, "RdBu_r", -vm, vm, "{:.1f}")
        else:
            im = _heat(a, g, dts, ratios, ttl, "magma_r", 0, 30, "{:.1f}")
        fig.colorbar(im, ax=a, shrink=0.8)
    ax[0].set_ylabel("ratio (inj / base)")
    fig.suptitle(title or "primary", fontsize=9)
    fig.savefig(os.path.join(out_dir, f"pileup_v9{tag}_bias.png"), dpi=130)
    plt.close(fig)
    for t in T:
        print(f"[V9{tag}] {'PASS' if t['passed'] else 'FAIL'}  {t['name']}  ({len(t['failing'])}/{t['n_cells']} failing)")
    print(f"wrote {out_dir}/pileup_v9{tag}.json, .md, _found.png, _bias.png")


def run_pileup(args):
    import glob
    import mf_pulse
    os.makedirs(args.out, exist_ok=True)
    sd = args.synth_dir
    with open(os.path.join(sd, "manifest.json")) as fh:
        man = json.load(fh)
    files = sorted(glob.glob(os.path.join(sd, "pileup_*.npz")))
    if not files:
        sys.exit(f"no pileup_*.npz in {sd}")
    h50 = args.templates or os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                         "output", "pulse_templates_h50.dat")
    wcfg = load_wave_cfg(args.daq_config)

    Z = []
    for f in files:
        z = np.load(f, allow_pickle=False)
        Z.append({k: z[k] for k in z.files})
    ns = Z[0]["waveforms"].shape[1]               # the manifest n_samples is not trusted
    names = sorted({str(x) for z in Z for x in z["base_channel_names"]})
    pwo4 = [nm.startswith("W") for nm in names]
    bank, info = mf_pulse.build_bank(names, pwo4, h50, wcfg, column="25320",
                                     int_tail_ratio=args.int_tail_ratio)
    idx = {nm: i for i, nm in enumerate(names)}
    Y = np.concatenate([z["waveforms"] for z in Z])
    mid = np.array([idx[str(x)] for z in Z for x in z["base_channel_names"]], dtype=np.int64)
    onsets = np.concatenate([z["truth_onset_t0_ns"][:, 1] for z in Z])
    print(f"[pileup] {len(files)} cells, {len(Y)} waveforms x {ns} samples (manifest says {man.get('n_samples')}); "
          f"{len(names)} base modules, {info['n_module']} h50 templates, "
          f"aggregate fallback: {','.join(info['fallback']) or 'none'}", flush=True)

    conv = dict(
        samples_per_waveform=ns, manifest_n_samples=man.get("n_samples"),
        template="per-module h50 (25320) of the base channel, aggregate fallback for "
                 + (",".join(info["fallback"]) or "none")
                 + f" ({info['n_module']} of {len(names)} base modules from file); the injected pulse itself "
                 f"uses the global PbWO4 template (tau_r {float(Z[0]['meta_tau_r_ns']):.3f}, "
                 f"tau_f {float(Z[0]['meta_tau_f_ns']):.3f} ns)",
        amplitude="mf A (unit-peak model) vs truth_heights_adc (peak-normalised, h_inj = ratio * base_height); "
                  "bias = (A_mf - A_true)/A_true, median and 68% half-width",
        time="primary: mf t50 (onset + leading-edge 50% time of the module template) vs truth onset + t50 of "
             "the global template; also onset (mf t0 vs truth_onset_t0_ns[:,0]); truth peak times are not used",
        matching=f"one-to-one greedy on t50 within +-{V9_MATCH_NS:g} ns to the injected and the base-ref truth pulse; "
                 "both_found = both truth pulses matched; both_clean = both matched and no unmatched mf pulse; "
                 "fp = mf pulses not matched to either truth; n>=2 = at least two mf pulses (the count-based "
                 "definition behind the WA and deconv numbers)",
        base_truth_note="the base pulse truth (column 1) is WaveAnalyzer/fit_pulse_shape output, not exact",
        base_onset_range_ns=[float(onsets.min()), float(onsets.max())])
    variants = [("", "", wcfg["ped_nsamples"])]
    if args.ped_diag > 0:
        variants.append((f"_ped{args.ped_diag}", f"diagnostic: pedestal window {args.ped_diag} samples", args.ped_diag))
    summary = {}
    for tag, title, ped_n in variants:
        rec = mf_pulse.MFReconstructor(bank, nsigma=args.nsigma, min_amp=args.min_amp, kmax=args.kmax,
                                       r_model=args.r_model, ped_nsamples=ped_n)
        cells, stats = _v9_cells(rec, Z, Y, mid, wcfg["overflow"])
        _v9_compare(cells, sd)
        T = _v9_targets(cells)
        vconv = dict(conv, mf_params=dict(
            nsigma=args.nsigma, min_amp=args.min_amp, kmax=args.kmax, r_model=args.r_model,
            int_tail_ratio=info["int_tail_ratio"], tail_break_n=info["tail_break_n"], ped_nsamples=ped_n,
            h50=h50, h50_sha1=info["sha1"]))
        if ped_n >= wcfg["ped_nsamples"] // 2 and onsets.min() < DT_SAMPLE * ped_n:
            vconv["pedestal_warning"] = (
                f"the base pulse onset ({onsets.min():.0f}-{onsets.max():.0f} ns) lies inside the first "
                f"{ped_n} samples ({DT_SAMPLE * ped_n:.0f} ns) that mf uses for the pedestal; these 100-sample "
                "synthetic waveforms start the pulse at sample 6-26, whereas data pulses sit later in the window")
        print(f"[pileup{tag}] mf: {stats}", flush=True)
        _v9_write(args.out, tag, title, sd, ns, vconv, stats, cells, T)
        summary[tag or "primary"] = [dict(name=t["name"], passed=t["passed"], n_failing=len(t["failing"])) for t in T]
    with open(os.path.join(args.out, "pileup_v9_summary.json"), "w") as fh:
        json.dump(summary, fh, indent=1)


# --------------------------------------------------------------- timecal
# Python port of analysis/tools/hycal_module_time_calib.cpp (:196-198 histogram
# set, :237-316 event loop, :327-380 fits, :403-560 network solve + JSON).
TC_NOTE = [
    "Port of hycal_module_time_calib.cpp working on benchmark pulse dumps, not ROOT replay trees.",
    "Event gate approximated: sum of E_MeV over the largest in-window (time-win) PbWO4 pulse per module of the chosen "
    "method > ebeam_frac*Ebeam and >= 4 PbWO4 modules with a pulse; the C++ requires exactly one cluster with nblocks>3 "
    "and E>0.8*Ebeam, TBIT_sum and nch<=70 (trigger bit and nch are not in the dump, nch<=70 is not applied).",
    "Pulses: PbWO4 channels with npk==1 for the method and amp >= height-min (rawsum pk.height, mf fitted A); raw time_ns.",
    "Pairs: 8-neighbours by row/col in the PbWO4 sector (is_neighbor same-sector rule) from hycal_map.json; dt = t(lower "
    "module index) - t(higher index); 200 bins over [-10,10] ns; entries (incl. overflow) must exceed min-entries.",
    "Fit: chi2 Gaussian (sqrt(n) errors, empty bins skipped, scipy curve_fit instead of ROOT Minuit) around the peak bin with "
    "half-width max(1, RMS), fallback full range; a fit whose mean lies outside its fit range is treated as failed (Python-only; the C++ has no such check); weight error = GetMeanError (RMS/sqrt(N in range)) floored at 0.01 ns.",
    "Solve: weighted least squares, reference offset 0, BFS connectivity, SVD pseudo-inverse; offset_error = sqrt(cov_ii). "
    "t_corrected = t_raw - offset_ns."]


def _gaus(x, a, mu, sg):
    return a * np.exp(-0.5 * ((x - mu) / sg) ** 2)


def fit_pair_hist(vals, n_entries, nb=200, lo=-10.0, hi=10.0, err_floor=0.01):
    """C++ fit block for one pair. vals = in-range dt values; returns (mu, sigma, mean_err, chi2/ndf) or None."""
    from scipy.optimize import curve_fit
    cnt, edges = np.histogram(vals, bins=nb, range=(lo, hi))
    ctr = 0.5 * (edges[1:] + edges[:-1])
    mean, rms = float(vals.mean()), float(vals.std())
    if cnt.max() <= 0 or not np.isfinite(mean):
        return None
    pk = ctr[int(np.argmax(cnt))]
    win = max(1.0, max(0.5, rms))
    mean_err = max(rms / math.sqrt(len(vals)), err_floor) if np.isfinite(rms) else err_floor
    for flo, fhi in ((max(lo, pk - win), min(hi, pk + win)), (lo, hi)):
        sel = (ctr >= flo) & (ctr <= fhi) & (cnt > 0)
        if sel.sum() < 3:
            continue
        x, y = ctr[sel], cnt[sel].astype(float)
        inr = (ctr >= flo) & (ctr <= fhi)
        w = cnt[inr] / max(cnt[inr].sum(), 1)
        s0 = float(np.sqrt(max(np.sum(w * (ctr[inr] - np.sum(w * ctr[inr])) ** 2), 0.01)))
        try:
            p, _ = curve_fit(_gaus, x, y, p0=(y.max(), x[np.argmax(y)], s0), sigma=np.sqrt(y), absolute_sigma=True, maxfev=2000)
        except Exception:
            continue
        mu, sg = float(p[1]), abs(float(p[2]))
        if not (np.isfinite(mu) and np.isfinite(sg) and sg > 0):
            continue
        if not (flo <= mu <= fhi):           # fit mean outside its window: treat as failed
            continue
        ndf = len(x) - 3
        chi2 = float(np.sum((y - _gaus(x, *p)) ** 2 / y))
        c2n = chi2 / ndf if ndf > 0 else float("inf")
        return mu, sg, mean_err, c2n
    return None


def timecal_pairs(d, method, rc, args, ebeam):
    """dt values and pair keys (lo*nmod+hi by module index) of one dump, plus n gated events."""
    nmod = len(rc)
    is_w = rc[:, 0] >= 0
    m = (d["method"] == method) & is_w[d["midx"]]
    t, key = d["time_ns"], d["event"] * 4096 + d["midx"]
    # event gate: largest in-window pulse per (event, module), summed per event
    w = np.flatnonzero(m & (t > args.time_win[0]) & (t < args.time_win[1]))
    w = w[last_per_group(key[w], d["E_MeV"][w])]
    ue, inv = np.unique(d["event"][w], return_inverse=True)
    esum, ncnt = np.bincount(inv, weights=d["E_MeV"][w]), np.bincount(inv)
    good_ev = ue[(esum > args.ebeam_frac * ebeam) & (ncnt >= 4)]
    sel = np.flatnonzero(m & (d["npk"] == 1) & (d["amp"] >= args.height_min) & np.isin(d["event"], good_ev))
    sel = sel[np.lexsort((d["midx"][sel], d["event"][sel]))]
    ev, mi, tt = d["event"][sel], d["midx"][sel], t[sel]
    dts, keys = [], []
    for k in range(1, 80):
        if len(ev) <= k:
            break
        same = ev[:-k] == ev[k:]
        if not same.any():
            break
        i = np.flatnonzero(same)
        a, b = mi[i], mi[i + k]
        nb = (np.abs(rc[a, 0] - rc[b, 0]) <= 1) & (np.abs(rc[a, 1] - rc[b, 1]) <= 1) & (a != b)
        i, a, b = i[nb], a[nb], b[nb]
        dts.append(tt[i] - tt[i + k])
        keys.append(a.astype(np.int64) * nmod + b)
    z = np.zeros(0)
    return (np.concatenate(dts) if dts else z), (np.concatenate(keys) if keys else z.astype(np.int64)), len(good_ev)


def run_timecal(args):
    db = db_dir()
    map_path = args.hycal_map or os.path.join(db, "hycal_map.json")
    ent = json.load(open(map_path))
    ds = [load_dump(p) for p in args.dumps]
    names = ds[0]["module_names"]
    nmod = len(names)
    by = {e["n"]: e for e in ent}
    rc = np.full((nmod, 2), -1, np.int64)       # (row, col) of PbWO4 modules, -1 otherwise
    for i, nm in enumerate(names):
        e = by.get(str(nm))
        if e is not None and e.get("t") == "PbWO4":
            rc[i] = (e["geo"]["row"], e["geo"]["col"])
    wmods = np.flatnonzero(rc[:, 0] >= 0)
    pairs_dt, pairs_key, n_ev = [], [], 0
    for d in ds:
        eb = args.ebeam or (float(d["ebeam"]) if "ebeam" in d else (2239.51 if d["run"].startswith("25") else None))
        if eb is None:
            sys.exit(f"timecal: no beam energy for run {d['run']}; pass --ebeam")
        dt, ky, ne = timecal_pairs(d, args.method, rc, args, eb)
        pairs_dt.append(dt)
        pairs_key.append(ky)
        n_ev += ne
        print(f"[timecal {d['run']}] Ebeam {eb:g}: {ne} gated events, {len(dt)} pair fills", flush=True)
    dt = np.concatenate(pairs_dt)
    ky = np.concatenate(pairs_key)
    o = np.argsort(ky, kind="stable")
    uk, st = np.unique(ky[o], return_index=True)
    parts = np.split(dt[o], st[1:])
    rows = []
    for k, v in zip(uk, parts):
        if len(v) <= args.min_entries:
            continue
        vin = v[(v >= -10) & (v < 10)]
        if len(vin) < 3:
            continue
        r = fit_pair_hist(vin, len(v))
        if r is not None and np.isfinite(r[3]):
            rows.append((int(k // nmod), int(k % nmod), *r))
    print(f"[timecal] {len(uk)} pairs with entries, {len(rows)} valid fits (min-entries {args.min_entries})", flush=True)
    ref = int(np.flatnonzero(names == args.ref)[0]) if args.ref in set(map(str, names)) else None
    if ref is None:
        sys.exit(f"timecal: reference {args.ref} not in module_names")
    adj = {}
    for i, j, *_ in rows:
        adj.setdefault(i, []).append(j)
        adj.setdefault(j, []).append(i)
    conn, queue = {ref}, [ref]
    while queue:
        c = queue.pop(0)
        for nx in adj.get(c, []):
            if nx not in conn:
                conn.add(nx)
                queue.append(nx)
    off, err = {ref: 0.0}, {ref: 0.0}
    unk = sorted(conn - {ref})
    if unk:
        ix = {m: n for n, m in enumerate(unk)}
        N = np.zeros((len(unk), len(unk)))
        rhs = np.zeros(len(unk))
        for i, j, mu, sg, me, c2 in rows:
            if i not in conn or j not in conn:
                continue
            wt = 1.0 / max(me, 0.01) ** 2
            if i != ref:
                N[ix[i], ix[i]] += wt
                rhs[ix[i]] += wt * mu
            if j != ref:
                N[ix[j], ix[j]] += wt
                rhs[ix[j]] -= wt * mu
            if i != ref and j != ref:
                N[ix[i], ix[j]] -= wt
                N[ix[j], ix[i]] -= wt
        U, sv, Vt = np.linalg.svd(N)
        inv_s = np.where(sv > sv.max() * 1e-14, 1.0 / sv, 0.0)
        cov = (Vt.T * inv_s) @ U.T
        sol = cov @ rhs
        for m, n in ix.items():
            off[m], err[m] = float(sol[n]), float(math.sqrt(max(cov[n, n], 0.0)))
    mods = []
    for m in wmods:
        nm = str(names[m])
        mods.append({"name": nm, "module_id": 1000 + int(nm[1:]), "offset_ns": off.get(m, 0.0),
                     "offset_error_ns": err.get(m, 0.0), "connected": m in conn, "resolved": m in off})
    js = {"reference_module": 1000 + int(args.ref[1:]), "reference_name": args.ref, "offset_units": "ns", "default": 0.0,
          "modules": mods, "note": TC_NOTE,
          "settings": {"method": args.method, "min_entries": args.min_entries, "height_min": args.height_min,
                       "ebeam_frac": args.ebeam_frac, "time_win_ns": list(args.time_win), "n_gated_events": n_ev,
                       "n_valid_pairs": len(rows), "n_connected": len(conn), "inputs": list(args.dumps)}}
    # control comparison against the existing table
    cmp_path = os.path.join(db, "hycal_time_offsets", "time_calib.json")
    if os.path.isfile(cmp_path):
        ref_by = {e["name"]: e for e in json.load(open(cmp_path))["modules"] if e.get("resolved")}
        dd = np.array([e["offset_ns"] - ref_by[e["name"]]["offset_ns"] for e in mods if e["resolved"] and e["name"] in ref_by])
        js["comparison_to_database_time_calib"] = {
            "n_modules": len(dd), "median_diff_ns": med(dd), "rms_diff_ns": float(np.sqrt(np.mean(dd ** 2))) if len(dd) else None,
            "frac_within_0p2ns": float(np.mean(np.abs(dd) <= 0.2)) if len(dd) else None,
            "note": "diff = new - database offset, modules resolved in both"}
        print("[timecal] vs database time_calib.json:", js["comparison_to_database_time_calib"])
    os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
    with open(args.out, "w") as fh:
        json.dump(sanitize(js), fh, indent=1)
    print(f"[timecal] connected {len(conn)}/{len(wmods)}; wrote {args.out}")


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("pulses", help="V5-V8 on pulse dump NPZ files")
    p.add_argument("dumps", nargs="+", help="dump NPZ files (one per run)")
    p.add_argument("--out", required=True, help="output directory")
    p.add_argument("--time-calib", default=None, help="time_calib.json (default: $PRAD2_DATABASE_DIR/hycal_time_offsets/time_calib.json)")
    p.add_argument("--hycal-map", default=None, help="hycal_map.json (default: $PRAD2_DATABASE_DIR/hycal_map.json)")
    p.add_argument("--time-calib-mf", default=None, help="time_calib-schema table applied to mf times only in the V7 "
                   "'mf with mf offsets vs rawsum with rawsum offsets' row (e.g. the timecal output)")
    p.set_defaults(fn=run_pulses)
    t = sub.add_parser("timecal", help="per-module time offsets from dumps (port of hycal_module_time_calib.cpp)")
    t.add_argument("dumps", nargs="+", help="dump NPZ files")
    t.add_argument("--method", type=int, choices=(0, 1), required=True, help="0 = rawsum, 1 = mf")
    t.add_argument("--out", required=True, help="output JSON (time_calib.json schema)")
    t.add_argument("--ref", default="W495", help="reference module name (offset 0)")
    t.add_argument("--min-entries", type=int, default=100, help="pair histogram needs > this many entries (C++: 100)")
    t.add_argument("--height-min", type=float, default=100.0, help="amplitude cut in ADC (C++: peak_height >= 100)")
    t.add_argument("--ebeam-frac", type=float, default=0.8)
    t.add_argument("--ebeam", type=float, default=None, help="beam energy MeV (default: dump key 'ebeam', else 2239.51 for runs 25xxx)")
    t.add_argument("--time-win", type=float, nargs=2, default=(130.0, 200.0), metavar=("LO", "HI"),
                   help="HyCal time window [ns] for the cluster-energy gate (runinfo hc_time_window)")
    t.add_argument("--hycal-map", default=None, help="hycal_map.json (default: $PRAD2_DATABASE_DIR/hycal_map.json)")
    t.set_defaults(fn=run_timecal)
    q = sub.add_parser("pileup", help="V9 pile-up truth test on a pileup_generator output directory")
    q.add_argument("synth_dir", help="pileup_synth_* directory (manifest.json + pileup_*.npz)")
    q.add_argument("--out", required=True, help="output directory")
    q.add_argument("--templates", default=None, help="h50 table (default: output/pulse_templates_h50.dat)")
    q.add_argument("--daq-config", default=None, help="daq_config.json (default: $PRAD2_DATABASE_DIR/daq_config.json)")
    q.add_argument("--nsigma", type=float, default=5.0)
    q.add_argument("--min-amp", type=float, default=5.0)
    q.add_argument("--kmax", type=int, default=4)
    q.add_argument("--r-model", type=float, default=0.05)
    q.add_argument("--ped-diag", type=int, default=5, help="also run with this pedestal window in samples (0 = off)")
    q.add_argument("--int-tail-ratio", type=float, default=None, help="override wave_cfg int_tail_ratio")
    q.set_defaults(fn=run_pileup)
    args = ap.parse_args()
    args.fn(args)


if __name__ == "__main__":
    main()
