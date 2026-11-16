#!/usr/bin/env python3
"""
cosmic_mip.py - vertical-track MIP peaks per region group and the LMS events of run 23744 (plan P14.5, D9, D4).

Usage (from analysis/pyscripts, `source setup_env.sh`, MPLBACKEND=Agg):
    python SS/cosmic_mip.py EVIO/prad_023744.evio.00000 [more files] --map SS/hycal_map_til_run023849.json
           [--out SS/out/mip_23744.tsv] [--lms-out SS/out/lms_23744.npz] [--max-events N] [--cache DIR]
One pass per file extracts (a) the vertical-track hits and (b) the LMS / alpha events.  Per-file hit tables are kept
in --cache (default: none) so that the group fit can be rerun without decoding.

Track selection (D9), per physics event without trigger bit 24 / 25 and with <= 200 channels:
  hit = PbWO4 channel whose peak sample max(x) - median(x[:8]) exceeds 10 ADC (lead glass is not counted);
  vertical track = a column holding a run of >= 3 contiguous rows of hits, and <= 2 hits elsewhere (other columns, or
  the same column outside the run).  Interior hits = the run without its two end modules; a hit enters the MIP
  statistics when WaveAnalyzer finds exactly one peak with integral > 0 (the event's integral is that peak's).

API:
  tracks(evio, map, max_events=0, pipe=None) -> hits    structured array, one row per interior hit:
        module (PbWO4 name), peak_adc, integral_adc (WaveAnalyzer), column_len (length of the contiguous run)
  scan_file(evio, map, max_events=0) -> dict(hits, n_events, n_tracks, n_multi, lms (dict), ...)
  region_groups(map) -> {group: set(module names)}   groups T, O, M, H, pool
  mip_peaks(hits, groups) -> list[(group, mpv, err, n)]   Gaussian fit in [0.7, 1.5] x the mode (centre of the
        maximum bin of 10 ADC bins, as analysis/tools/cosmic_test.cpp:98-110) of the integral spectrum

Region groups: the energy_resolution_benchmark.py definitions (module_regions, read-only, constants replicated):
  H = rows 15-20 and columns 15-20 (1-based; EXCLUDE_ROWCOL), T = outermost 2 PbWO4 rows/columns
  (TRANSITION_LAYERS) not in H, then by the module-centre radius: pool = I1 + I2 + I3 (r < 200 mm; REGION_EDGES
  120 | 160 | 200), M = 200-300 mm, O = >= 300 mm.  The groups are disjoint; the table also has an `all` row.

LMS npz layout (--lms-out; PbWO4 modules in hycal_map module-index order of the derived map, names in mod_names):
  mod_names[nW], mod_x[nW], mod_y[nW]            module names and centres
  ts[n], ev[n], lms_file[n]                       TI timestamp, event number, file index of the LMS events
  lms_I[n, nW] float32                            WaveAnalyzer integral (exactly 1 peak, > 0; else NaN)
  lms_ref[n, 3] float32                           same for LMS1, LMS2, LMS3 (by name through the map)
  alpha_ts[m], alpha_ev[m], alpha_file[m], alpha_ref[m, 3]   alpha-source events (bit 25, not bit 24), reference PMTs
  files[nf] (basenames), n_all[nf], n_lms[nf], n_flag_phys[nf]  per file: physics events, LMS events kept,
        excluded nch > 900 events with neither bit 24 nor bit 25 (in 23744 the trigger bits are 0x100 only: about
        5.8 k per file, 3 % of the events, not the 910 of the brief; alpha events with nch > 900 stay alpha events)
  trig_hi_keys[k], trig_hi_counts[k]              trigger bits of all nch > 900 events
LMS event = bit 24 set and nch > 900 (D4).  Events with nch > 900 and neither bit 24 nor 25 are counted
(n_flag_phys) and excluded.  Bit-24 events with nch <= 900 are counted in n_lms_low and excluded.
"""
from __future__ import annotations

import argparse
import json
import os
import sys
import time
from collections import Counter
from pathlib import Path

import numpy as np

SS = Path(__file__).resolve().parent
PYSCRIPTS = SS.parents[1]
if str(PYSCRIPTS) not in sys.path:
    sys.path.insert(0, str(PYSCRIPTS))
if str(SS) not in sys.path:
    sys.path.insert(0, str(SS))
import snake_scan as SSN  # noqa: E402

TBIT_LMS, TBIT_ALPHA = 1 << 24, 1 << 25
HIT_ADC = 10.0            # D9 hit threshold on the peak sample
MIN_RUN = 3               # contiguous modules of a vertical track
MAX_OFF = 2               # hits outside the run
MAX_NCH_TRACK = 200       # events with more channels are not track candidates
NCH_LMS = 900             # LMS event: bit 24 and more channels than this
BIN_ADC = 10.0
HIST_MAX = 800.0
GROUPS = ("T", "O", "M", "H", "pool")
MEV_PER_ADC = 0.114       # D3: half-plateau crossing of the endpoint, 2108.85 MeV at 18.5 k ADC (P14.2 measured 15.4 k ADC: 0.137)
HIT_DTYPE = [("module", "U8"), ("peak_adc", "f4"), ("integral_adc", "f4"), ("column_len", "i2")]


# ---------------------------------------------------------------------------------------------------------------
# region groups (energy_resolution_benchmark.module_regions, constants replicated)
# ---------------------------------------------------------------------------------------------------------------

def region_groups(map_path=SSN.DERIVED_MAP) -> dict:
    ents = [e for e in json.loads(Path(map_path).read_text()) if e["t"] == "PbWO4"]
    x = np.array([e["geo"]["x"] for e in ents])
    y = np.array([e["geo"]["y"] for e in ents])
    row = np.array([e["geo"]["row"] for e in ents])      # 1-based, as BM's row + 1
    col = np.array([e["geo"]["col"] for e in ents])
    code = np.digitize(np.hypot(x, y), (120.0, 160.0, 200.0, 300.0))
    h = (15 <= row) & (row <= 20) & (15 <= col) & (col <= 20)
    n = 2
    t = (row < row.min() + n) | (row > row.max() - n) | (col < col.min() + n) | (col > col.max() - n)
    reg = np.where(h, "H", np.where(t, "T", np.where(code <= 2, "pool", np.where(code == 3, "M", "O"))))
    return {g: {e["n"] for e, r in zip(ents, reg) if r == g} for g in GROUPS}


# ---------------------------------------------------------------------------------------------------------------
# one pass over a file
# ---------------------------------------------------------------------------------------------------------------

def _find_track(hit_mods, rc):
    """hit_mods: list of module indices with a hit; rc: index -> (row, col).  Returns the run (list of module
    indices ordered by row) or None."""
    bycol = {}
    for m in hit_mods:
        r, c = rc[m]
        bycol.setdefault(c, []).append((r, m))
    for c, lst in bycol.items():
        if len(lst) < MIN_RUN:
            continue
        lst.sort()
        best, cur = [], [lst[0]]
        for a, b in zip(lst, lst[1:]):
            if b[0] == a[0] + 1:
                cur.append(b)
            else:
                best = cur if len(cur) > len(best) else best
                cur = [b]
        best = cur if len(cur) > len(best) else best
        if len(best) >= MIN_RUN and len(hit_mods) - len(best) <= MAX_OFF:
            return [m for _, m in best]
    return None


def scan_file(evio, map_path=SSN.DERIVED_MAP, max_events: int = 0, pipe=None) -> dict:
    from prad2py import dec
    evio = str(evio)
    map_path = Path(map_path).resolve()          # a relative path would be resolved by the C++ side under database/
    if not map_path.exists():
        raise FileNotFoundError(map_path)
    pipe = pipe or SSN.open_pipeline(str(map_path), evio)
    hy = pipe.hycal
    mods = [hy.module(i) for i in range(hy.module_count())]
    pw_idx = [i for i, m in enumerate(mods) if m.is_pwo4()]
    col_of = {mi: k for k, mi in enumerate(pw_idx)}
    ref_idx = {i: int(m.name[3]) - 1 for i, m in enumerate(mods)
               if len(m.name) == 4 and m.name.startswith("LMS") and m.name[3] in "123"}
    rc = {i: (int(mods[i].row), int(mods[i].column)) for i in pw_idx}
    nW = len(pw_idx)
    names = [mods[i].name for i in pw_idx]
    ana = pipe.wave_ana
    cache = {}

    def mod_of(crate, slot, c):
        k = (crate, slot, c)
        if k not in cache:
            m = hy.module_by_daq(crate, slot, c)
            cache[k] = None if m is None else m.index
        return cache[k]

    ch = SSN.open_evio(evio, pipe)
    hits, lms = [], dict(ts=[], ev=[], I=[], ref=[])
    alp = dict(ts=[], ev=[], ref=[])
    st = Counter()
    trig_hi = Counter()
    peaks_all = []
    while ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            ch.select_event(i)
            inf = ch.info()
            tb = int(inf.trigger_bits)
            e = ch.fadc()
            chans = []           # (module index, samples)
            for ri in range(e.nrocs):
                roc = e.roc(ri)
                if not roc.present:
                    continue
                crate = pipe.crate_map.get(roc.tag)
                if crate is None:
                    continue
                for s in roc.present_slots():
                    slot = roc.slot(s)
                    for c in slot.present_channels():
                        mi = mod_of(crate, s, c)
                        if mi is not None:
                            chans.append((mi, slot.channel(c)))
            nch = len(chans)
            st["n_all"] += 1
            is_l, is_a = bool(tb & TBIT_LMS), bool(tb & TBIT_ALPHA)
            if nch > NCH_LMS:
                trig_hi[tb] += 1
                if is_l:
                    row = np.full(nW, np.nan, np.float32)
                    ref = np.full(3, np.nan, np.float32)
                    for mi, cd in chans:
                        if cd.nsamples <= 0 or not (mi in ref_idx or mi in col_of):
                            continue
                        _, _, pks = ana.analyze(cd.samples)
                        if len(pks) == 1 and pks[0].integral > 0:
                            if mi in ref_idx:
                                ref[ref_idx[mi]] = pks[0].integral
                            else:
                                row[col_of[mi]] = pks[0].integral
                    lms["ts"].append(int(inf.timestamp)); lms["ev"].append(int(inf.event_number))
                    lms["I"].append(row); lms["ref"].append(ref)
                elif is_a:
                    ref = np.full(3, np.nan, np.float32)
                    for mi, cd in chans:
                        if mi in ref_idx and cd.nsamples > 0:
                            _, _, pks = ana.analyze(cd.samples)
                            if len(pks) == 1 and pks[0].integral > 0:
                                ref[ref_idx[mi]] = pks[0].integral
                    alp["ts"].append(int(inf.timestamp)); alp["ev"].append(int(inf.event_number)); alp["ref"].append(ref)
                else:
                    st["n_flag_phys"] += 1
                continue
            if is_l:
                st["n_lms_low"] += 1
                continue
            if is_a:
                ref = np.full(3, np.nan, np.float32)
                for mi, cd in chans:
                    if mi in ref_idx and cd.nsamples > 0:
                        _, _, pks = ana.analyze(cd.samples)
                        if len(pks) == 1 and pks[0].integral > 0:
                            ref[ref_idx[mi]] = pks[0].integral
                alp["ts"].append(int(inf.timestamp)); alp["ev"].append(int(inf.event_number)); alp["ref"].append(ref)
                continue
            if nch > MAX_NCH_TRACK or nch < MIN_RUN:
                continue
            hv = {}
            for mi, cd in chans:
                if mi in col_of and cd.nsamples >= 9:
                    x = np.asarray(cd.samples, np.float64)
                    v = x.max() - np.median(x[:8])
                    if v > HIT_ADC:
                        hv[mi] = (v, cd)
            if len(hv) < MIN_RUN:
                continue
            run = _find_track(list(hv), rc)
            if run is None:
                continue
            st["n_tracks"] += 1
            st["track_len_sum"] += len(run)
            for mi in run[1:-1]:
                v, cd = hv[mi]
                _, _, pks = ana.analyze(cd.samples)
                if len(pks) == 1 and pks[0].integral > 0:
                    hits.append((mods[mi].name, v, pks[0].integral, len(run)))
                else:
                    st["n_multi"] += 1
            if max_events and st["n_all"] >= max_events:
                break
        if max_events and st["n_all"] >= max_events:
            break
    return dict(hits=np.array(hits, HIT_DTYPE), counts=dict(st), trig_hi=dict(trig_hi), lms=lms, alpha=alp,
                names=names,
                mod_x=np.array([mods[i].x for i in pw_idx], np.float32),
                mod_y=np.array([mods[i].y for i in pw_idx], np.float32),
                W=np.array([mods[i].cal_factor for i in pw_idx], np.float32))


def tracks(evio, map_path=SSN.DERIVED_MAP, max_events: int = 0, pipe=None):
    """hits(module, peak_adc, integral_adc, column_len) of the interior modules of the vertical tracks."""
    return scan_file(evio, map_path, max_events, pipe)["hits"]


# ---------------------------------------------------------------------------------------------------------------
# MIP peak fit
# ---------------------------------------------------------------------------------------------------------------

def fit_peak(vals, bin_w=BIN_ADC, hi=HIST_MAX):
    """Gaussian fit in [0.7, 1.5] x mode; returns (mean, err, n_in_window, mode).  NaN when the fit fails."""
    from scipy.optimize import curve_fit
    vals = np.asarray(vals, float)
    nb = int(round(hi / bin_w))
    cnt, edges = np.histogram(vals, bins=nb, range=(0, hi))
    if cnt.sum() < 10:
        return np.nan, np.nan, 0, np.nan
    ctr = 0.5 * (edges[1:] + edges[:-1])
    mode = ctr[np.argmax(cnt)]
    sel = (ctr >= 0.7 * mode) & (ctr <= 1.5 * mode) & (cnt > 0)
    if sel.sum() < 4:
        return np.nan, np.nan, int(sel.sum()), mode
    f = lambda x, a, m, s: a * np.exp(-0.5 * ((x - m) / s) ** 2)   # noqa: E731
    try:
        p, cov = curve_fit(f, ctr[sel], cnt[sel], p0=(cnt.max(), mode, 0.2 * mode),
                           sigma=np.sqrt(cnt[sel]), absolute_sigma=True, maxfev=5000)
        return float(p[1]), float(np.sqrt(cov[1, 1])), int(cnt[(ctr >= 0.7 * mode) & (ctr <= 1.5 * mode)].sum()), mode
    except Exception:
        return np.nan, np.nan, int(sel.sum()), mode


def mip_peaks(hits, groups):
    """[(group, mpv, err, n)] of the integral spectrum per group (n = hits in the group)."""
    out = []
    for g, mem in groups.items():
        sel = np.isin(hits["module"], list(mem))
        m, e, _, _ = fit_peak(hits["integral_adc"][sel])
        out.append((g, m, e, int(sel.sum())))
    return out


# ---------------------------------------------------------------------------------------------------------------
# outputs
# ---------------------------------------------------------------------------------------------------------------

def group_table(per_file, names, W, groups, files):
    """rows: dict per group (+ all) with pooled and per-file fits; per_file = list of hits arrays."""
    wmap = dict(zip(names, W))
    allh = np.concatenate(per_file) if per_file else np.array([], HIT_DTYPE)
    rows = []
    for g in list(groups) + ["all"]:
        mem = groups.get(g)

        def pick(h):
            return h if mem is None else h[np.isin(h["module"], list(mem))]
        h = pick(allh)
        mpv, err, _, _ = fit_peak(h["integral_adc"])
        w = np.array([wmap[m] for m in h["module"]], float)
        mev, mev_e, _, _ = fit_peak(h["integral_adc"] * w, bin_w=BIN_ADC * 0.13, hi=HIST_MAX * 0.13)
        pf = []
        for hf in per_file:
            hh = pick(hf)
            m, e, _, _ = fit_peak(hh["integral_adc"])
            pf.append((m, e, len(hh)))
        rows.append(dict(group=g, mpv=mpv, err=err, n_hits=len(h), n_modules=len(set(h["module"].tolist())),
                         peak_median=float(np.median(h["peak_adc"])) if len(h) else np.nan,
                         mev_nominal=mpv * MEV_PER_ADC, mev_nominal_err=err * MEV_PER_ADC, mev_w=mev, mev_w_err=mev_e,
                         w_median=float(np.median(w)) if len(w) else np.nan, per_file=pf,
                         n_members=len(mem) if mem is not None else len(names)))
    return rows


def write_mip_tsv(path, rows, files, map_path, stats):
    fn = [Path(f).name for f in files]
    cols = ["group", "mpv", "err", "n_hits", "n_modules", "n_members", "peak_median", "mev_nominal", "mev_nominal_err",
            "mev_w", "mev_w_err", "w_median"]
    for f in fn:
        s = f.split(".")[-1]
        cols += [f"mpv_{s}", f"err_{s}", f"n_{s}"]
    with open(path, "w") as fh:
        fh.write(f"# source: {' '.join(fn)} ; map {Path(map_path).name} ; evio dir {Path(files[0]).resolve().parent}\n")
        fh.write("# " + "; ".join(f"{Path(f).name}: {s}" for f, s in zip(files, stats)) + "\n")
        fh.write("# mpv: Gaussian core fit in [0.7,1.5] x mode of the WaveAnalyzer integral (ADC), interior hits of "
                 "vertical tracks (D9); mev_nominal = mpv x 0.114 MeV/ADC; mev_w = fit of integral x per-module "
                 "production W (hycal calibration of the run); T,O,M,H,pool as energy_resolution_benchmark regions"
                 "; note: MEV_PER_ADC = 0.114 is the plan's value from an 18.5 k ADC edge, P14.2 measured the edge at 15.4 k ADC (about 0.137 MeV/ADC), so the mev_nominal columns are about 17 % low (group ratios unaffected)\n")
        fh.write("\t".join(cols) + "\n")
        for r in rows:
            v = [r["group"]] + [f"{r[k]:.4g}" if isinstance(r[k], float) else str(r[k]) for k in
                                ("mpv", "err", "n_hits", "n_modules", "n_members", "peak_median", "mev_nominal",
                                 "mev_nominal_err", "mev_w", "mev_w_err", "w_median")]
            for m, e, n in r["per_file"]:
                v += [f"{m:.4g}", f"{e:.3g}", str(n)]
            fh.write("\t".join(v) + "\n")


def write_lms_npz(path, results, files):
    nW = len(results[0]["names"])
    def cat(key, sub, shape, dtype):
        a = [np.asarray(r[key][sub], dtype).reshape((len(r[key][sub]),) + shape[1:]) for r in results]
        return np.concatenate(a)
    lms_file = np.concatenate([np.full(len(r["lms"]["ts"]), k, np.int16) for k, r in enumerate(results)])
    alp_file = np.concatenate([np.full(len(r["alpha"]["ts"]), k, np.int16) for k, r in enumerate(results)])
    th = Counter()
    for r in results:
        th.update(r["trig_hi"])
    keys = np.array(sorted(th), np.int64)
    np.savez_compressed(
        path, mod_names=np.array(results[0]["names"]), mod_x=results[0]["mod_x"], mod_y=results[0]["mod_y"],
        ts=cat("lms", "ts", (-1,), np.int64), ev=cat("lms", "ev", (-1,), np.int64), lms_file=lms_file,
        lms_I=cat("lms", "I", (-1, nW), np.float32), lms_ref=cat("lms", "ref", (-1, 3), np.float32),
        alpha_ts=cat("alpha", "ts", (-1,), np.int64), alpha_ev=cat("alpha", "ev", (-1,), np.int64), alpha_file=alp_file,
        alpha_ref=cat("alpha", "ref", (-1, 3), np.float32),
        files=np.array([Path(f).name for f in files]),
        n_all=np.array([r["counts"].get("n_all", 0) for r in results], np.int64),
        n_lms=np.array([len(r["lms"]["ts"]) for r in results], np.int64),
        n_flag_phys=np.array([r["counts"].get("n_flag_phys", 0) for r in results], np.int64),
        n_lms_low=np.array([r["counts"].get("n_lms_low", 0) for r in results], np.int64),
        trig_hi_keys=keys, trig_hi_counts=np.array([th[k] for k in keys], np.int64))


def _scan_cached(f, map_path, max_events, cache):
    if cache:
        p = Path(cache) / (Path(f).name + ".pkl")
        if p.exists():
            import pickle
            return pickle.load(open(p, "rb"))
    t0 = time.time()
    r = scan_file(f, map_path, max_events)
    r["seconds"] = time.time() - t0
    if cache:
        import pickle
        Path(cache).mkdir(parents=True, exist_ok=True)
        pickle.dump(r, open(p, "wb"))
    return r


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[1])
    ap.add_argument("evio", nargs="+")
    ap.add_argument("--map", default=str(SSN.DERIVED_MAP))
    ap.add_argument("--out", default="")
    ap.add_argument("--lms-out", default="")
    ap.add_argument("--max-events", type=int, default=0)
    ap.add_argument("--cache", default="")
    a = ap.parse_args(argv)
    run = Path(a.evio[0]).name.split(".")[0].split("_")[-1].lstrip("0")
    out = Path(a.out) if a.out else SS / "out" / f"mip_{run}.tsv"
    res = [_scan_cached(f, a.map, a.max_events, a.cache) for f in a.evio]
    groups = region_groups(a.map)
    rows = group_table([r["hits"] for r in res], res[0]["names"], res[0]["W"], groups, a.evio)
    stats = []
    for f, r in zip(a.evio, res):
        c = r["counts"]
        s = (f"{c.get('n_all', 0)} events, {c.get('n_tracks', 0)} tracks, {len(r['hits'])} interior hits, "
             f"{c.get('n_multi', 0)} multi-peak, {len(r['lms']['ts'])} LMS, {c.get('n_flag_phys', 0)} flagged nch>900 "
             f"without bit 24/25, {len(r['alpha']['ts'])} alpha")
        stats.append(s)
        print(f"{Path(f).name}: {s} ({r.get('seconds', 0):.0f} s)", flush=True)
    write_mip_tsv(out, rows, a.evio, a.map, stats)
    if a.lms_out:
        write_lms_npz(a.lms_out, res, a.evio)
    for r in rows:
        print(f"{r['group']:5s} mpv {r['mpv']:.1f} +- {r['err']:.1f} ADC  n {r['n_hits']}  modules {r['n_modules']}  "
              f"peak median {r['peak_median']:.1f}  MeV {r['mev_nominal']:.1f} / {r['mev_w']:.1f}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
