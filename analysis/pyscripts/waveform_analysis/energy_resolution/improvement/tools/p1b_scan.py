#!/usr/bin/env python3
"""p1b_scan.py <out.npz> <evio files...>: AMP P1b rescan (diagnostic only; the benchmark is not edited, no default changes).

Runs `energy_resolution_benchmark.py scan --set PEAK_COLUMNS=True` with two wrappers installed in the module: pulse_peaks (extra per-hit columns
computed from the same subtracted waveform and peak sample as h_apk_raw) and numpy's savez_compressed (adds the columns to the scan file).
h_apk_par0 is replaced by h_id (hit id, float32, exact below 2^24); every other column is the benchmark's own and must be identical to a plain scan.
Extra columns (ADC, baseline-subtracted, other fitted pulses subtracted; NaN outside the waveform / clipped / WaveAnalyzer hits):
  h_s3 h_s4 h_s5 h_s7   fixed sums of the samples [p-1,p+1], [p-1,p+2], [p-1,p+3], [p-2,p+4] around the peak sample p
  h_c5 h_c7             fixed sums [jc, jc+4], [jc-1, jc+5] anchored at the CFD sample jc (last sample below half of the peak before p; template free)
  h_psd h_c1 h_sref     std of the clipped 30 pedestal samples of the waveform, their lag-1 correlation, the benchmark's per-module noise reference
  h_sq h_st             (every 20th hit and the stored windows only, NaN elsewhere) sum of the unit-peak template squared / unrenormalised over the fit window (white-noise MF: sigma_a = sigma / sqrt(h_sq))
  h_dp                  p - (t0 + t_pk) / DT, the peak sample against the template peak (samples)
  sub_id sub_win sub_mwin   for every 4th hit (by id) with h_a >= SUB_AMIN: samples [p-5, p+10] and the fitted model a * template at them
"""
import os
import pickle
import sys
from pathlib import Path
import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[2]))
import energy_resolution_benchmark as bm        # noqa: E402

SUB_AMIN, W0, W1 = 40.0, 5, 11                  # stored window: samples p - W0 .. p + W1 - 1
_orig_pp, _orig_rr, _orig_save = bm.pulse_peaks, bm.run_ref, np.savez_compressed
ST = {"sig": None, "n": 0}
COLS = ("s3", "s4", "s5", "s7", "c5", "c7", "psd", "c1", "sref", "sq", "st", "dp")
EX = {c: [] for c in COLS}
SUB = {"id": [], "win": [], "mwin": []}


def rr(mf, Y, midx, sig, b_ref):
    ST["sig"] = np.asarray(sig, dtype=np.float64)
    return _orig_rr(mf, Y, midx, sig, b_ref)


def wsum(ys, lo, hi):
    """Sum of ys[h, lo..hi] (inclusive, int arrays); NaN when the window leaves the waveform."""
    ns = ys.shape[1]
    c = np.concatenate([np.zeros((ys.shape[0], 1)), np.cumsum(ys, axis=1)], axis=1)
    h = np.arange(ys.shape[0])
    ok = (lo >= 0) & (hi <= ns - 1)
    lo, hi = np.clip(lo, 0, ns - 1), np.clip(hi, 0, ns - 1)
    return np.where(ok, c[h, hi + 1] - c[h, lo], np.nan)


def pp(bank, mf, Y, base, m, t0, a, v, r, j):
    raw, par, par0, tc = _orig_pp(bank, mf, Y, base, m, t0, a, v, r, j)
    DT, F32 = bm.DT, bm.F32
    n, ns = Y.shape
    y = Y.astype(np.float64) - base.astype(np.float64)[:, None]
    ii = np.arange(ns)[None, :]
    oth = np.zeros((r.size, ns))
    multi = v.sum(axis=1) >= 2
    rows2 = np.flatnonzero(multi)
    if rows2.size:
        ridx = np.full(n, -1)
        ridx[rows2] = np.arange(rows2.size)
        contrib = np.stack([np.where(v[rows2, k], a[rows2, k], 0.0)[:, None] * mf.basis(m[rows2], t0[rows2, k].astype(F32), ns)
                            for k in range(bm.KMAX)], axis=1)
        hm = np.flatnonzero(multi[r])
        oth[hm] = contrib[ridx[r[hm]]].sum(axis=1) - contrib[ridx[r[hm]], j[hm]]
    ys = y[r] - oth
    h = np.arange(r.size)
    near = np.abs(DT * ii - (t0[r, j] + bank.t_pk[m[r]])[:, None]) <= bm.PEAK_SEARCH_NS
    p = np.where(near, ys, -np.inf).argmax(axis=1)
    clipped = (Y >= bm.OVERFLOW).any(axis=1)[r]
    ok = ~clipped
    assert np.allclose(ys[h, p][ok], raw[ok]), "peak sample mismatch with pulse_peaks"
    half = bm.PEAK_CFD_FRAC * ys[h, p]
    jj = np.where((ys < half[:, None]) & (ii < p[:, None]), ii, -1).max(axis=1)
    nanc = lambda x: np.where(clipped, np.nan, x)             # noqa: E731
    s3, s4, s5, s7 = (nanc(wsum(ys, p + lo, p + hi)) for lo, hi in ((-1, 1), (-1, 2), (-1, 3), (-2, 4)))
    c5 = nanc(np.where(jj >= 0, wsum(ys, jj, jj + 4), np.nan))
    c7 = nanc(np.where(jj >= 0, wsum(ys, jj - 1, jj + 5), np.nan))
    x = Y[:, :bm.PED_NSAMPLES].astype(np.float64)               # pedestal statistics per waveform
    med = np.median(x, axis=1)
    sg = np.maximum(1.4826 * np.median(np.abs(x - med[:, None]), axis=1), bm.SIGMA_MIN)
    good = np.abs(x - med[:, None]) <= np.maximum(2.0 * sg, 1.0)[:, None]
    ped = (x * good).sum(1) / np.maximum(good.sum(1), 1)
    res = (x - ped[:, None]) * good
    psd = np.sqrt((res ** 2).sum(1) / np.maximum(good.sum(1) - 1, 1))
    pair = good[:, :-1] & good[:, 1:]
    c1 = (res[:, :-1] * res[:, 1:] * pair).sum(1) / np.maximum((res ** 2).sum(1), 1e-9) * (good.sum(1) / np.maximum(pair.sum(1), 1))
    base_id = ST["n"]
    ids = base_id + np.arange(r.size)
    sub = np.flatnonzero((ids % 20 == 0) | (ok & (a[r, j] >= SUB_AMIN) & (ids % 4 == 0)))   # template quantities on a subsample only (the basis spline is the slow part)
    bas = mf.basis(m[r][sub], t0[r, j][sub].astype(F32), ns).astype(np.float64)
    w0 = (np.ceil(t0[r, j][sub] / DT).astype(np.int64) - bm.LEAD)[:, None]
    wm = (ii >= w0) & (ii < w0 + bm.LW)
    sq, st = np.full(r.size, np.nan), np.full(r.size, np.nan)
    sq[sub], st[sub] = (bas ** 2 * wm).sum(1), (bas * wm).sum(1)
    dp = p - (t0[r, j] + bank.t_pk[m[r]]) / DT
    sref = ST["sig"][r]
    ST["n"] += r.size
    for c, val in zip(COLS, (s3, s4, s5, s7, c5, c7, psd[r], c1[r], sref, sq, st, dp)):
        EX[c].append(np.asarray(val, dtype=np.float32))
    pos = np.full(r.size, -1)
    pos[sub] = np.arange(sub.size)
    sel = np.flatnonzero(ok & (a[r, j] >= SUB_AMIN) & (ids % 4 == 0) & (p >= W0) & (p + W1 <= ns))
    if sel.size:
        k = p[sel][:, None] + np.arange(-W0, W1)[None, :]
        SUB["id"].append(ids[sel])
        SUB["win"].append(ys[sel[:, None], k].astype(np.float32))
        SUB["mwin"].append((a[r, j][sel][:, None] * bas[pos[sel][:, None], k]).astype(np.float32))
    return raw, par, ids.astype(np.float64), tc


def save(out, **kw):
    try:
        ids = kw.pop("h_apk_par0").astype(np.float64)
        assert np.nanmax(ids) < 2 ** 24, "hit id beyond float32 range"
        ex = {c: np.concatenate(EX[c]) for c in COLS}
        fin = np.isfinite(ids)
        ii = np.where(fin, ids, 0).astype(np.int64)
        for c in COLS:
            kw["h_" + c] = np.where(fin, ex[c][ii], np.nan).astype(np.float32)
        kw["h_id"] = ids.astype(np.float32)
        for q in SUB:
            kw["sub_" + q] = np.concatenate(SUB[q]) if SUB[q] else np.zeros(0)
        print(f"[p1b] extra columns merged: {len(ids)} hits, {kw['sub_id'].size} stored windows", flush=True)
    except Exception as e:                         # keep the work: raw extras to a pickle next to the output
        print(f"[p1b] MERGE FAILED {e!r}; extras pickled", flush=True)
        pickle.dump((EX, SUB), open(str(out) + ".extras.pkl", "wb"))
    return _orig_save(out, **kw)


if __name__ == "__main__":
    bm.pulse_peaks, bm.run_ref, np.savez_compressed = pp, rr, save
    sys.argv = [sys.argv[0], "scan", "--set", "PEAK_COLUMNS=True"] + sys.argv[1:]
    bm.main()
