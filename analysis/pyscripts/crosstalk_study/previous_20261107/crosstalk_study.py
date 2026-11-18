#!/usr/bin/env python3
"""crosstalk_study.py - dt = T_module - T_seed, amplitude and chi2/dof of the non-seed cluster members (one run type).

INPUTS
    output/r5ri_20261104/jobs/final.txt   the chain's final-stage job line of each split (config source; gains and poscorr map of the
                                          chain are reused, nothing in it is modified); database/ and runinfo/ are read only.
    EVIO splits under /Users/L00338853/work/PRad/data/evio (read by the benchmark and by the WaveAnalyzer pass).
    The frozen benchmark_hycal_timing_test_20261021.py (sha 3293ab96...) is re-run once per split with the chain's final-stage arguments
    plus --cluster-times --pulse-dump (default-off diagnostics); an existing output directory (<work>/<run>_<split>/ with rc and
    pulse.npz of the same job line) is ingested instead of re-running.

PER-HIT QUANTITIES (cache/<run-type>/<run>_<split>_<path>.npz, see README for the column list)
    cluster members, MF t / amp / chi2dof, WaveAnalyzer CFD time t_wa (method-0 peak nearest to the MF time within 10 ns),
    time_algo and log-normal fit chi2/dof of that peak (per-channel dec.WaveAnalyzer.analyze_fit pass over the EVIO),
    cluster energy (E_raw, includes gain file and pos-corr), E'_ep at the seed centre, r = E/E', Moller tag bits, module indices.

OPTIONS
    --selection ep|single   ep (default): single-cluster, not Moller-tagged (pair|kin|rej), |r-1| < 3 sigma_CB;  single: all
    --time-source cfd|mf    cfd (default): WaveAnalyzer CFD time;  mf: matched-filter time
    --cfd-algo fit|all      cfd: only time_algo == 2 at seed and module (default) | any algorithm
    --tcal off|auto|<json>  per-module time offsets applied to seed and module (auto = database/hycal_time_offsets/time_calib.json);
                            with offsets, the raw figures go to <out>_raw
    --selfref               also derive self-referenced offsets (cross-check) -> crosstalk_study/tcal/
    --path legacy|gated     --reuse (use the cache)  --jobs N  --work DIR (keep the benchmark outputs)  --splits RUN/SPLIT ...

USAGE: see README.md (from analysis/pyscripts after `source setup_env.sh`).
"""
import argparse
import re
import shutil
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

PY = Path(__file__).resolve().parents[1]
FINAL = PY / "output/r5ri_20261104/jobs/final.txt"
HERE = Path(__file__).resolve().parent
DEFAULT_SPLITS = ["25228/50"]
CACHE_KEYS = ("r", "moller", "cl_r", "cl_moller", "sigma_cb", "n_inner", "t_wa", "seed_t_wa", "mod_i", "seed_mod_i", "module_names")   # columns every analysis needs
WA_MATCH = 10.0        # ns: max |t_WaveAnalyzer - t_MF| for the WaveAnalyzer peak matched to an MF pulse
TYPE_SPLITS = {"0p73GeV_H2": [f"24590/{i}" for i in range(10, 101, 10)] + [f"24603/{i}" for i in range(50, 301, 50)],
               "2p24GeV_H2": [f"25228/{i}" for i in range(50, 551, 50)],
               "3p49GeV_H2": [f"24386/{i}" for i in (10, 20, 30, 40, 50, 70, 80, 90, 100, 110)]
               + [f"24503/{i}" for i in range(10, 101, 10)] + [f"24512/{i}" for i in range(10, 71, 10)]}


def run_splits(rtype):
    """'<run-type>_<run>' (e.g. 2p24GeV_Ta_X17_25800): all splits of that single run in the chain's final.txt (Ta is analysed per run)."""
    m = re.match(r".*_(\d{5})$", rtype)
    if not m:
        return None
    return [f"{int(m.group(1))}/{int(x)}" for x in re.findall(r"scan_%06d/split_(\d+) " % int(m.group(1)), FINAL.read_text())]


def final_cmd(run, split):
    stem = f"output/r5ri_20261104/scan/splits/scan_{run:06d}/split_{split:05d}"
    for ln in FINAL.read_text().splitlines():
        if ln.startswith(stem + " "):
            return ln.split(" ", 1)[1], stem
    sys.exit(f"[ERROR] {run}/{split} not in {FINAL}")


def run_benchmark(run, split, work, thr=None):
    cmd, stem = final_cmd(run, split)
    out = work / f"{run:06d}_{split:05d}"
    out.mkdir(parents=True, exist_ok=True)
    cmd = cmd.replace(stem + "/result", str(out / "result"))
    cmd = cmd.split(" ", 1)[1] if cmd.startswith("/") else cmd          # drop the chain's interpreter path
    full = f"{sys.executable} {cmd} --cluster-times --pulse-dump {out / 'pulse.npz'}" + (f" --signal-min {thr:g}" if thr else "")
    t0 = time.time()
    if (out / "rc").exists() and (out / "rc").read_text().split()[0] in ("0", "rc=0") and (out / "pulse.npz").exists():
        return out, 0.0                          # existing benchmark output of the same job line (kept with --work)
    with open(out / "benchmark.log", "w") as lg:
        rc = subprocess.run(full, shell=True, cwd=PY, stdout=lg, stderr=subprocess.STDOUT).returncode
    (out / "rc").write_text(f"{rc}\n")
    if rc:
        sys.exit(f"[ERROR] benchmark rc={rc}, see {out / 'benchmark.log'}")
    return out, time.time() - t0


def chain_sigma(stem, path):
    """(sigma_over_E_cb, n_inner_events) of the chain's result.tsv for this path (CB core sigma / mu of the inner ring)."""
    path = {"legacy": "legacy", "gated": "new"}[path]            # the benchmark's tables call the gated path "new"
    f = PY / stem / "result.tsv"
    if not f.exists():
        return float("nan"), 0
    lines = f.read_text().splitlines()
    for i, ln in enumerate(lines):
        if ln.startswith("path\tn_modules\tn_inner_events"):
            h = ln.split("\t")
            for ln2 in lines[i + 1:]:
                p = ln2.split("\t")
                if p[0] == path:
                    return float(p[h.index("sigma_over_E_cb")]), int(p[h.index("n_inner_events")])
            break
    return float("nan"), 0


def reduce_split(out, path, stem):
    """-> dict of per-hit arrays (one entry per non-seed cluster member) plus per-cluster arrays (cl_*)."""
    sys.path.insert(0, str(PY))
    import benchmark_hycal_timing_test_20261021 as B          # frozen; only expected_energy_elastic is used
    txt = (out / "result.tsv").read_text()[:20000]
    ebeam = float(re.search(r"resolved Ebeam=([0-9.]+)", txt).group(1))
    mtar = float(re.search(r"mass_MeV=([0-9.]+)", txt).group(1))
    import pandas as pd
    wt = pd.read_csv(out / "result_window.tsv", sep="\t", comment="#")
    wt = wt[wt["path"] == {"legacy": "legacy", "gated": "new"}[path]]
    hz = float(np.median(wt["r_mm"].to_numpy() / np.tan(np.radians(wt["theta_deg"].to_numpy()))))   # HyCal z (mm) from the benchmark's own table
    z = np.load(out / "pulse.npz", allow_pickle=True)
    names = z["module_names"]
    ev, mi, tt, aa, cc, me = (z[k] for k in ("event", "midx", "time_ns", "amp", "chi2dof", "method"))
    k = me == 1
    k0 = me == 0                                    # WaveAnalyzer peaks (time_ns = CFD pickoff, algo not recorded)
    key0 = z["event"][k0].astype(np.int64) * 4096 + z["midx"][k0]
    o0 = np.argsort(key0, kind="stable")
    key0, tt0 = key0[o0], z["time_ns"][k0][o0]
    ev, mi, tt, aa, cc = ev[k], mi[k], tt[k], aa[k], cc[k]
    key = ev.astype(np.int64) * 4096 + mi
    o = np.argsort(key, kind="stable")
    key, tt, aa, cc = key[o], tt[o], aa[o], cc[o]
    idx = {n: i for i, n in enumerate(names)}

    def wa(evn, name, t):
        """WaveAnalyzer peak time nearest to the MF time t within WA_MATCH ns (nan if none)."""
        kk = int(evn) * 4096 + idx[name]
        a, b = np.searchsorted(key0, kk), np.searchsorted(key0, kk, "right")
        if a == b:
            return np.nan
        j = a + int(np.argmin(np.abs(tt0[a:b] - t)))
        return tt0[j] if abs(tt0[j] - t) < WA_MATCH else np.nan

    def lookup(evn, name, t):
        kk = int(evn) * 4096 + idx[name]
        a, b = np.searchsorted(key, kk), np.searchsorted(key, kk, "right")
        if a == b:
            return None
        j = a + int(np.argmin(np.abs(tt[a:b] - t)))
        return (tt[j], aa[j], cc[j]) if abs(tt[j] - t) < 0.01 else None

    cols = {c: [] for c in ("ev", "seed_t", "seed_amp", "seed_E", "t", "amp", "chi2", "E", "mod_i", "seed_mod_i",
                            "cl_E", "eprime", "r", "moller", "seed_t_wa", "t_wa")}
    ccols = {c: [] for c in ("cl_ev", "cl_r", "cl_moller")}
    n_cl = n_nomatch = 0
    with open(out / "result_clusters.tsv") as f:
        hdr = next(f).rstrip("\n").split("\t")
        ic, ie, im = hdr.index("config_name"), hdr.index("evnum"), hdr.index("members")
        iE, isx, isy = hdr.index("cluster_E_raw_MeV"), hdr.index("seed_x_mm"), hdr.index("seed_y_mm")
        ipr, ikn, irj = (hdr.index(c) for c in ("moller_pair", "moller_kin", "moller_rej"))
        for ln in f:
            p = ln.rstrip("\n").split("\t")
            if p[ic] != path or not p[im]:
                continue
            mem = [m.split(":") for m in p[im].split(";")]
            ev_ = int(p[ie])
            s = lookup(ev_, mem[0][0], float(mem[0][2]))
            if s is None:
                continue
            swa = wa(ev_, mem[0][0], s[0])
            n_cl += 1
            # cluster_E_raw_MeV already includes --gain-file and --pos-corr (cal/raw = 1 with --gain-refine off) and is the
            # energy the benchmark's moller_kin tag uses; E' = elastic kinematic energy at the seed-module centre.
            ecl = float(p[iE])
            th = np.degrees(np.arctan(np.hypot(float(p[isx]), float(p[isy])) / hz))
            ep = B.expected_energy_elastic(float(th), ebeam, mtar)
            tag = int(p[ipr]) + 2 * int(p[ikn]) + 4 * int(p[irj])      # bit0 pair, bit1 kin, bit2 rej
            ccols["cl_ev"].append(ev_); ccols["cl_r"].append(ecl / ep); ccols["cl_moller"].append(tag)
            for n, e, t in mem[1:]:
                r = lookup(ev_, n, float(t))
                if r is None:
                    n_nomatch += 1
                    continue
                for c, v in zip(cols, (ev_, s[0], s[1], float(mem[0][1]), r[0], r[1], r[2], float(e), idx[n], idx[mem[0][0]],
                                       ecl, ep, ecl / ep, tag, swa, wa(ev_, n, r[0]))):
                    cols[c].append(v)
    dts = {"ev": np.int32, "moller": np.int8, "mod_i": np.int16, "seed_mod_i": np.int16}
    res = {c: np.asarray(v, dtype=dts.get(c, np.float32)) for c, v in cols.items() }
    res.update({c: np.asarray(v, dtype={"cl_ev": np.int32, "cl_moller": np.int8}.get(c, np.float32)) for c, v in ccols.items()})
    res["sigma_cb"], res["n_inner"] = (np.array(x) for x in chain_sigma(stem, path))
    res["module_names"] = np.asarray(names).astype(str)
    res["ebeam"], res["eprime_hz"] = np.array(ebeam), np.array(hz)
    res["n_clusters"] = np.array(n_cl)
    res["n_nomatch"] = np.array(n_nomatch)
    res["n_events_dump"] = np.array(len(np.unique(z["event"])))
    return res


def evio_of(run, split):
    cmd, _ = final_cmd(run, split)
    m = re.search(r"\S*prad_%06d\.evio\.%05d" % (run, split), cmd)
    return m.group(0), re.search(r"--hc-calib-file (\S+)", cmd).group(1)


def wa_fit_split(job):
    """Per-channel WaveAnalyzer pass (prad2py dec.WaveAnalyzer.analyze_fit) over the EVIO split for the (event, module) pairs of the
    cached hits and their seeds; adds time_algo (0 peak / 1 linear CFD / 2 log-normal fit CFD) and the log-normal fit chi2/dof of
    the WaveAnalyzer peak matched to t_wa (|dt| < 0.01 ns).  Updates the cache npz in place (columns algo, lnchi2, seed_algo, seed_lnchi2)."""
    cachef, run, split = job
    sys.path.insert(0, str(PY))
    import _common as C
    from prad2py import dec
    z = dict(np.load(cachef))
    mnames = [str(x) for x in z["module_names"]]
    evio, calib = evio_of(run, split)
    p = C.setup_pipeline(evio_path=evio, hc_calib_file=calib)
    names = {m: i for i, m in enumerate(mnames)}
    daq = {}
    for n, i in names.items():
        m = p.hycal.module_by_name(n)
        if m is not None and m.is_hycal():
            daq[i] = (m.daq.crate, m.daq.slot, m.daq.channel)
    pairs = np.unique(np.concatenate([z["ev"].astype(np.int64) * 4096 + z["mod_i"], z["ev"].astype(np.int64) * 4096 + z["seed_mod_i"]]))
    need = {}
    for k in pairs:
        need.setdefault(int(k) // 4096, []).append(int(k) % 4096)
    ev0 = int(z["ev"].min())
    rows = []
    ch = dec.EvChannel()
    ch.set_config(p.cfg)
    if ch.open_auto(evio) != dec.Status.success:
        sys.exit(f"[ERROR] cannot open {evio}")
    ana = p.wave_ana
    while ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            e = d["event"]
            evn = int(e.info.event_number)
            want = need.get(evn)
            if not want:
                continue
            rocs = {}
            for ri in range(e.nrocs):
                roc = e.roc(ri)
                if roc.present and p.crate_map.get(roc.tag) is not None:
                    rocs[p.crate_map[roc.tag]] = roc
            for mi in want:
                a = daq.get(mi)
                roc = rocs.get(a[0]) if a else None
                if roc is None:
                    continue
                cd = roc.slot(a[1]).channel(a[2])
                if cd.nsamples <= 0:
                    continue
                _, _, pks, fits = ana.analyze_fit(cd.samples)
                for pk, ft in zip(pks, fits):
                    rows.append(((evn - ev0) * 4096 + mi, pk.time, pk.time_algo, ft.chi2_per_dof if ft.ok else np.nan))
    T = np.array(rows, dtype=np.float64).reshape(-1, 4)
    T = T[np.lexsort((T[:, 1], T[:, 0]))]
    kk = T[:, 0] * 1000.0 + T[:, 1]

    def join(ev, mi, t):
        q = ((ev.astype(np.int64) - ev0) * 4096 + mi).astype(np.float64) * 1000.0 + t
        j = np.clip(np.searchsorted(kk, q), 1, max(len(kk) - 1, 1))
        best = np.where(np.abs(kk[j - 1] - q) <= np.abs(kk[j] - q), j - 1, j) if len(kk) > 1 else np.zeros(len(q), int)
        ok = (np.abs(kk[best] - q) < 0.01) & ~np.isnan(t) if len(kk) else np.zeros(len(q), bool)
        return np.where(ok, T[best, 2], -1).astype(np.int8), np.where(ok, T[best, 3], np.nan).astype(np.float32)
    z["algo"], z["lnchi2"] = join(z["ev"], z["mod_i"], z["t_wa"])
    z["seed_algo"], z["seed_lnchi2"] = join(z["ev"], z["seed_mod_i"], z["seed_t_wa"])
    np.savez_compressed(cachef, **z)
    return cachef, len(T)


def add_nmod(cachef, pulsef):
    """Adds per-event numbers of fired PbWO4 modules from the benchmark pulse dump (all fired channels of every kept event are in the dump):
    nmod_wa20 = distinct modules with a WaveAnalyzer peak (method 0) of height >= 20 ADC (the N_FIRED criterion of waveform_good_bad.py / H_MIN = 20),
    nmod_mf5 = distinct modules with an MF pulse (method 1, amplitude >= min_amp 5 ADC).  Hit-level columns nmod_*, cluster-level cl_nmod_*."""
    z = dict(np.load(cachef))
    pz = np.load(pulsef, allow_pickle=True)
    ev, mi, me, amp = pz["event"], pz["midx"], pz["method"], pz["amp"]
    for tag, sel in (("wa20", (me == 0) & (amp >= 20)), ("mf5", me == 1)):
        key = np.unique(ev[sel].astype(np.int64) * 4096 + mi[sel])
        u, c = np.unique(key // 4096, return_counts=True)
        def look(e):
            j = np.clip(np.searchsorted(u, e), 0, max(len(u) - 1, 0))
            return np.where(u[j] == e, c[j], 0).astype(np.int16) if len(u) else np.zeros(len(e), np.int16)
        z[f"nmod_{tag}"] = look(z["ev"].astype(np.int64)); z[f"cl_nmod_{tag}"] = look(z["cl_ev"].astype(np.int64))
    np.savez_compressed(cachef, **z)
    return z


def load_or_make(run, split, rtype, path, reuse, work, need_ep=True, thr=None):
    cache = HERE / "cache" / rtype
    cache.mkdir(parents=True, exist_ok=True)
    f = cache / f"{run:06d}_{split:05d}_{path}.npz"
    if reuse and f.exists():
        d = dict(np.load(f))
        if all(k in d for k in CACHE_KEYS):
            return d, None
        print(f"[cache] {f.name} lacks columns {[k for k in CACHE_KEYS if k not in d]}: re-reducing "
              f"(benchmark output of {work} is reused when present, else the benchmark is rerun)", flush=True)
    out, el = run_benchmark(run, split, work, thr)
    res = reduce_split(out, path, final_cmd(run, split)[1])
    np.savez_compressed(f, **res)
    print(f"[split] {run}/{split}: {int(res['n_clusters'])} clusters, {len(res['t'])} non-seed hits "
          f"({int(res['n_nomatch'])} without MF row), benchmark {el:.0f} s", flush=True)
    return res, el


def core_widths(dt):
    """(robust sigma = half of the 16-84 % quantile range, sigma of a Gaussian fit to median +- 3 ns with 0.1 ns bins)."""
    from scipy.optimize import curve_fit
    q16, q50, q84 = np.percentile(dt, [16, 50, 84])
    h, e = np.histogram(dt, bins=np.arange(q50 - 3, q50 + 3.001, 0.1))
    c = 0.5 * (e[1:] + e[:-1])
    try:
        p, _ = curve_fit(lambda x, A, m, sg: A * np.exp(-0.5 * ((x - m) / sg) ** 2), c, h, p0=[h.max(), q50, 0.5 * (q84 - q16)])
        g = abs(float(p[2]))
    except Exception:
        g = float("nan")
    return 0.5 * float(q84 - q16), g


def load_offsets(mode, mnames):
    """--tcal off|auto|<path> -> (per-module offset array ns indexed like mnames, label) or (None, 'off (raw)').
    auto = database/hycal_time_offsets/time_calib.json, applied as t = t - offset like the benchmark (resolved entries only)."""
    if mode == "off":
        return None, "off (raw)"
    import hashlib, json
    f = Path(mode) if mode != "auto" else PY.parents[1] / "database/hycal_time_offsets/time_calib.json"
    raw = f.read_bytes()
    tab = {m["name"]: float(m["offset_ns"]) for m in json.loads(raw)["modules"] if m.get("resolved")}
    off = np.array([tab.get(str(n), 0.0) for n in mnames])
    return off, f"{f.name} (sha1 {hashlib.sha1(raw).hexdigest()[:10]}, {len(tab)} resolved)"


def selfref_offsets(tm, ts, mi, si, nmod, official=None, label="", passes=0, win=20.0, nmin=50):
    """Self-referenced per-module offsets (cross-check only).  Per (module, seed-module) pair the median of dt = T_module - T_seed over
    |dt| < win ns (then 5, 3, 3 ns around the current solution) estimates o_module - o_seed; the offsets solve the weighted (sqrt(n)) least-squares network of these pair medians
    (pairs with >= 20 hits), gauge: median offset of the measured modules (>= nmin hits) = 0.  Same idea as the neighbour-pair network
    of hycal_module_time_calib.cpp, on the cached e-p sample; the CFD/MF time basis is whatever is passed in."""
    from scipy.sparse import coo_matrix
    from scipy.sparse.linalg import lsqr
    n = np.bincount(mi, minlength=nmod)
    off = np.zeros(nmod)
    its = 0
    for w_ in (win, 5.0, 3.0, 3.0):                  # re-centre the pair windows on the current solution (tail-robust)
        dt = (tm - off[mi]) - (ts - off[si])
        k = np.abs(dt) < w_
        key = mi[k] * 4096 + si[k]
        o = np.lexsort((dt[k], key))
        key, dk = key[o], dt[k][o]
        u, st, cnt = np.unique(key, return_index=True, return_counts=True)
        good = cnt >= 20
        u, st, cnt = u[good], st[good], cnt[good]
        med = dk[st + cnt // 2]
        pm, ps = u // 4096, u % 4096
        w = np.sqrt(cnt)
        r = np.arange(len(u))
        A = coo_matrix((np.concatenate([w, -w, [1e3] * nmod]),
                        (np.concatenate([r, r, np.full(nmod, len(u))]), np.concatenate([pm, ps, np.arange(nmod)]))),
                       shape=(len(u) + 1, nmod)).tocsr()
        sol = lsqr(A, np.concatenate([w * med, [0.0]]), atol=1e-10, btol=1e-10, iter_lim=5000)
        off = off + sol[0]
        its += int(sol[2])
    it = its
    meas = n >= nmin
    off -= np.median(off[meas])
    dtc = (tm - off[mi]) - (ts - off[si]); k = np.abs(dtc) < win
    o = np.lexsort((dtc[k], mi[k])); mk, dk = mi[k][o], dtc[k][o]
    st = np.searchsorted(mk, np.arange(nmod)); en = np.searchsorted(mk, np.arange(nmod), "right")
    se = np.full(nmod, np.nan)
    for m in np.nonzero(meas)[0]:
        se[m] = 1.2533 * np.std(dk[st[m]:en[m]]) / np.sqrt(en[m] - st[m])
    res = dict(off=off, n=n, se=se, passes=int(it))
    if official is not None:
        both = meas & (official != 0)
        dd = off[both] - official[both]
        dd -= np.median(dd)
        res["rms_vs_official"], res["n_cmp"] = float(np.sqrt(np.mean(dd ** 2))), int(both.sum())
        gp = (official[pm] != 0) & (official[ps] != 0)
        res["rms_pair"] = float(np.sqrt(np.mean(((off[pm] - off[ps]) - (official[pm] - official[ps]))[gp] ** 2)))   # local (neighbour-pair) agreement
    return res


def write_ratio(row, label):
    """tables/ratio_<slug>_<path>.json for this run, and tables/crosstalk_ratio.md collecting all of them (main sample only)."""
    import json
    td = HERE / "tables"; td.mkdir(exist_ok=True)
    if "off" not in label[:3]:
        tag = (("_single" if row["selection"] != "ep" else "") + (f"_e{row['cluster_thr']:g}" if row.get("cluster_thr") else "")) + (f"_nmod{row['nmod_max']}{row['nmod_def']}" if row.get("nmod_max") else "") + ("" if (row["time_source"] == "cfd" and row["cfd_algo"] == "fit")
                                                                else f"_{row['time_source']}" + ("" if row["time_source"] == "mf" else "all"))
        (td / f"ratio_{row['run_type']}_{row['path']}{tag}.json").write_text(json.dumps(row))
    L = ["# Crosstalk and core ratios (selection column: ep = single-cluster, no Moller tag, |r-1|<3 sigma; single = all single-cluster events; calibrated, WaveAnalyzer CFD log-normal fit only; dt = T_module - T_seed)", "",
         "(a) events with >= 1 non-seed hit in the region / total events of the selection; (b) hits in the region / total non-seed hits; h/ev = region hits / total e-p events.", "",
         "| run type | selection | cluster_thr_MeV | n_mod | path | time | total events | total non-seed hits | crosstalk window [ns] | (a) | (b) | h/ev | core window | (a') | (b') | h/ev' |",
         "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for f in sorted(td.glob("ratio_*.json")):
        r = json.loads(f.read_text())
        tl = "MF" if r["time_source"] == "mf" else ("CFD fit" if r["cfd_algo"] == "fit" else "CFD any")
        L.append(f"| {r['run_type']} | {'e-p' if r['selection'] == 'ep' else 'single'} | {('%g' % r['cluster_thr']) if r.get('cluster_thr') else '500 x kE'} | {('<%d (%s)' % (r['nmod_max'], r['nmod_def'])) if r.get('nmod_max') else 'all'} | {r['path']} | {tl} | {r['n_ev']:,} | {r['n_hits']:,} | {r['ct_lo']:g}..{r['ct_hi']:g} | "
                 f"{100 * r['n_ev_window'] / r['n_ev']:.2f} % | {100 * r['n_hits_window'] / r['n_hits']:.2f} % | {r['n_hits_window'] / r['n_ev']:.3f} | "
                 f"abs(dt) < {r['core']:g} | {100 * r['n_ev_core'] / r['n_ev']:.2f} % | {100 * r['n_hits_core'] / r['n_hits']:.2f} % | "
                 f"{r['n_hits_core'] / r['n_ev']:.3f} |")
    (td / "crosstalk_ratio.md").write_text("\n".join(L) + "\n")


def plot_all(d, outdir, title, extra="", a_algo_fit=False, ct=(7.0, 15.0), short=""):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.colors import LogNorm
    outdir.mkdir(parents=True, exist_ok=True)
    dt = d["t"] - d["seed_t"]                      # dt = T_module - T_seed, ns (positive = module later than the seed)
    amp, sa, chi = d["amp"], d["seed_amp"], d["chi2"]
    n = len(dt)
    summ = {}
    # fig 1a
    out10 = np.mean(np.abs(dt) > 10)
    summ_w = core_widths(dt)
    wtxt = (f"core: robust sigma (16-84 %) {summ_w[0]:.3f} ns, Gauss(+-3 ns) sigma {summ_w[1]:.3f} ns"
            + (f"\n{extra}" if extra else ""))
    summ["n"], summ["frac_out10"] = n, out10
    summ["w_robust"], summ["w_gauss"] = summ_w
    summ["n_late10"], summ["n_early10"] = int(np.sum(dt > 10)), int(np.sum(dt < -10))   # module later / earlier than seed
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.hist(np.clip(dt, -60, 60), bins=np.arange(-60, 60.01, 0.5), color="C0")
    ax.set_yscale("log")
    ax.axvline(-10, c="k", ls=":"); ax.axvline(10, c="k", ls=":")
    ax.axvspan(7, 15, color="orange", alpha=0.15, label="crosstalk region (10-30 note)"); ax.legend(loc="upper left", fontsize=8)  # fig1a only
    ax.set_xlabel("dt = T_module - T_seed  [ns]  (clipped to +-60)"); ax.set_ylabel("non-seed hits / 0.5 ns")
    ax.text(0.97, 0.95, f"N = {n}\n|dt| > 10 ns: {100 * out10:.2f} %\n{wtxt}", transform=ax.transAxes,
            ha="right", va="top", fontsize=9, bbox=dict(fc="w", ec="0.7"))
    ax.set_title(title, fontsize=9); fig.tight_layout(); fig.savefig(outdir / "fig1a_dt_hist.png", dpi=170); plt.close(fig)

    # fig 1c: zoom of fig 1a, -5..25 ns, 0.25 ns bins, linear and log y
    zb = np.arange(-5, 25.01, 0.25)
    nwin = int(np.sum((dt >= zb[0]) & (dt <= zb[-1])))
    fig, axs = plt.subplots(1, 2, figsize=(11, 5), sharex=True)
    for ax, ylog in zip(axs, (False, True)):
        ax.hist(dt, bins=zb, color="C0")
        ax.axvspan(ct[0], ct[1], color="orange", alpha=0.15)
        if ylog:
            ax.set_yscale("log")
        ax.set_xlabel("dt = T_module - T_seed  [ns]"); ax.set_ylabel("non-seed hits / 0.25 ns" + (" (log)" if ylog else ""))
    fig.suptitle(short or title.split("\n")[0], fontsize=10); fig.tight_layout()
    fig.savefig(outdir / "fig1c_dt_hist_zoom.png", dpi=170); plt.close(fig)

    # fig 1b: dt = T_module - T_seed on -5..25 ns (0.25 ns bins) vs non-seed module amplitude (log y), median amp per dt bin
    ab = np.logspace(np.log10(max(amp[amp > 0].min(), 5)), np.log10(amp.max()), 60)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    h = ax.hist2d(dt, np.clip(amp, ab[0], ab[-1]), bins=[zb, ab], norm=LogNorm(), cmap="viridis")
    fig.colorbar(h[3], ax=ax, label="hits")
    ax.set_yscale("log"); ax.set_xlim(zb[0], zb[-1]); ax.set_ylim(ab[0], ab[-1])
    ax.set_xlabel("dt = T_module - T_seed  [ns]  (window -5 to 25 ns, 0.25 ns bins)")
    ax.set_ylabel("non-seed module MF amplitude  [ADC]")
    ax.set_title(title, fontsize=9); fig.tight_layout(); fig.savefig(outdir / "fig1b_dt_vs_amp.png", dpi=170); plt.close(fig)

    # fig 2: dt = T_module - T_seed on -5..25 ns (0.25 ns bins) vs MF chi2/dof of the non-seed module (linear y, log colour)
    ok = np.isfinite(chi)
    win = ok & (dt >= zb[0]) & (dt <= zb[-1])
    ymax = float(np.ceil(np.percentile(chi[win], 99)))          # upper limit covering ~99 % of in-window hits
    frac_above = float(np.mean(chi[win] > ymax))
    yb = np.linspace(0, ymax, 60)
    fig, ax = plt.subplots(figsize=(7.5, 5.5))
    h = ax.hist2d(dt[ok], np.clip(chi[ok], 0, ymax), bins=[zb, yb], norm=LogNorm(), cmap="viridis")
    fig.colorbar(h[3], ax=ax, label="hits")
    inn, outt = ok & (np.abs(dt) < 5), ok & (np.abs(dt) > 20)
    summ["chi_in"], summ["chi_out"] = np.median(chi[inn]), (np.median(chi[outt]) if outt.any() else np.nan)
    summ["n_in"], summ["n_out"] = int(inn.sum()), int(outt.sum())
    summ["chi_ymax"], summ["frac_chi_above"] = ymax, frac_above
    ax.text(0.97, 0.55, f"median chi2/dof: |dT|<5 ns {summ['chi_in']:.2f} (N={summ['n_in']})\n"
            f"|dT|>20 ns {summ['chi_out']:.2f} (N={summ['n_out']})", transform=ax.transAxes, ha="right", va="center",
            fontsize=8, bbox=dict(fc="w", ec="0.7"))
    ax.text(0.03, 0.97, f"chi2/dof > {ymax:.0f}: {100 * frac_above:.2f} % of in-window hits\n"
            f"(clipped into top row)", transform=ax.transAxes, ha="left", va="top", fontsize=8,
            bbox=dict(fc="w", ec="0.7"))
    ax.set_xlim(zb[0], zb[-1]); ax.set_ylim(0, ymax)
    ax.set_xlabel("dt = T_module - T_seed  [ns]  (window -5 to 25 ns, 0.25 ns bins)")
    ax.set_ylabel(f"MF chi2/dof of the non-seed module fit  (linear, clipped at {ymax:.0f})")
    ax.set_title(title, fontsize=9)
    fig.tight_layout(); fig.savefig(outdir / "fig2_dt_vs_chi2.png", dpi=170); plt.close(fig)

    # fig 2b: dt vs log-normal-fit chi2/dof of the non-seed module (CFD fit-only sets)
    if "lnchi2" in d and a_algo_fit:
        lc = d["lnchi2"]; okc = np.isfinite(lc)
        winc = okc & (dt >= zb[0]) & (dt <= zb[-1])
        ym = float(np.ceil(np.percentile(lc[winc], 99)))
        fig, ax = plt.subplots(figsize=(7.5, 5.5))
        h = ax.hist2d(dt[okc], np.clip(lc[okc], 0, ym), bins=[zb, np.linspace(0, ym, 60)], norm=LogNorm(), cmap="viridis")
        fig.colorbar(h[3], ax=ax, label="hits")
        ax.set_xlim(zb[0], zb[-1]); ax.set_ylim(0, ym)
        ax.set_xlabel("dt = T_module - T_seed  [ns]  (window -5 to 25 ns, 0.25 ns bins)")
        ax.set_ylabel(f"log-normal fit chi2/dof, non-seed module\n(linear, clipped at {ym:.0f})")
        ax.set_title(title, fontsize=9); fig.tight_layout(); fig.savefig(outdir / "fig2b_dt_vs_lnchi2.png", dpi=170); plt.close(fig)
        summ["lnchi2_median"] = float(np.median(lc[okc]))

    # fig 3
    lo, hi = 5, 5000.0                              # fixed axes: identical for all run types and for fig3 / fig4 / fig4b
    b = np.logspace(np.log10(lo), np.log10(hi), 90)
    fig, ax = plt.subplots(figsize=(7.5, 6.5))
    h = ax.hist2d(np.clip(sa, lo, hi), np.clip(amp, lo, hi), bins=[b, b], norm=LogNorm(), cmap="viridis")
    fig.colorbar(h[3], ax=ax, label="hits")
    x = np.array([lo, hi])
    ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
    ax.set_xlabel("seed module MF amplitude  [ADC]"); ax.set_ylabel("non-seed module MF amplitude  [ADC]")
    ax.set_title(title, fontsize=9)
    fig.tight_layout(); fig.savefig(outdir / "fig3_amp_vs_amp.png", dpi=170); plt.close(fig)

    # fig 4 / 4b: seed vs non-seed amplitude for the crosstalk window and for the core, binning and axes of fig 3
    wsel = (dt >= ct[0]) & (dt <= ct[1])
    csel = np.abs(dt) < 5
    for sel, fn, lab in ((wsel, "fig4_amp2d_crosstalk.png", f"crosstalk window {ct[0]:g} <= dt <= {ct[1]:g} ns"),
                         (csel, "fig4b_amp2d_core.png", "core |dt| < 5 ns")):
        fig, ax = plt.subplots(figsize=(7.5, 6.5))
        h = ax.hist2d(np.clip(sa[sel], lo, hi), np.clip(amp[sel], lo, hi), bins=[b, b], norm=LogNorm(), cmap="viridis")
        fig.colorbar(h[3], ax=ax, label="hits")
        ax.set_xscale("log"); ax.set_yscale("log"); ax.set_xlim(lo, hi); ax.set_ylim(lo, hi)
        ax.set_xlabel("seed module MF amplitude  [ADC]"); ax.set_ylabel("non-seed module MF amplitude  [ADC]")
        ax.text(0.03, 0.95, f"N = {int(sel.sum()):,}", transform=ax.transAxes, ha="left", va="top", fontsize=9,
                bbox=dict(fc="w", ec="0.7"))
        ax.set_title(title + "\n" + lab, fontsize=9)
        fig.tight_layout(); fig.savefig(outdir / fn, dpi=170); plt.close(fig)
    summ["ct_n"] = int(wsel.sum())
    summ["ct_med_amp"] = float(np.median(amp[wsel])) if wsel.any() else float("nan")
    summ["ct_med_ratio"] = float(np.median(amp[wsel] / sa[wsel])) if wsel.any() else float("nan")
    summ["ct_frac10"] = float(np.mean(amp[wsel] < 10)) if wsel.any() else float("nan")
    summ["core_n"] = int(csel.sum())
    summ["core_med_amp"] = float(np.median(amp[csel])); summ["core_med_ratio"] = float(np.median(amp[csel] / sa[csel]))
    summ["core_frac10"] = float(np.mean(amp[csel] < 10))
    return summ


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--run-type", default="2p24GeV_H2")
    ap.add_argument("--out", required=True)
    ap.add_argument("--splits", nargs="+", default=None, metavar="RUN/SPLIT",
                    help="default: all splits of the run type in the reference chain (TYPE_SPLITS), else 25228/50")
    ap.add_argument("--jobs", type=int, default=1, help="parallel benchmark jobs")
    ap.add_argument("--path", choices=["legacy", "gated"], default="legacy")
    ap.add_argument("--reuse", action="store_true")
    ap.add_argument("--time-source", choices=["cfd", "mf"], default="cfd",
                    help="cfd: WaveAnalyzer CFD time (method-0 peak of the pulse dump nearest to the MF pulse); mf: matched-filter time")
    ap.add_argument("--cfd-algo", choices=["fit", "all"], default="fit",
                    help="cfd time source: fit = only hits with the log-normal fit CFD (time_algo 2) at seed and module (default); "
                         "all = any time_algo (0 peak / 1 linear CFD / 2 fit)")
    ap.add_argument("--cluster-thr", type=float, default=None, metavar="MEV",
                    help="single-cluster threshold: passes --signal-min MEV to the benchmark (default: the benchmark's 500*kE MeV); "
                         "cache goes to cache/<run-type>_e<MEV>")
    ap.add_argument("--selection", choices=["ep", "single"], default="ep",
                    help="ep: single-cluster, no Moller tag (pair|kin|rej), |r-1| < 3 sigma (default); single: all single-cluster events")
    ap.add_argument("--tcal", default="auto", help="off | auto (database/hycal_time_offsets/time_calib.json) | path of a time_calib json")
    ap.add_argument("--ct-window", nargs=2, type=float, default=[7.0, 15.0], metavar=("LO", "HI"),
                    help="crosstalk window in dt = T_module - T_seed [ns] for fig4 (default 7 15)")
    ap.add_argument("--core-window", type=float, default=2.0, metavar="W", help="core |dt| < W ns for the core ratios (default 2)")
    ap.add_argument("--nmod-max", type=int, default=None, metavar="N",
                    help="keep only events with fewer than N fired modules (see --nmod-def); needs the benchmark outputs in --work the first time")
    ap.add_argument("--nmod-def", choices=["wa20", "mf5"], default="wa20",
                    help="wa20: modules with a WaveAnalyzer peak >= 20 ADC (N_FIRED of waveform_good_bad.py); mf5: modules with an MF pulse (>= 5 ADC)")
    ap.add_argument("--no-raw", action="store_true", help="do not write the <out>_raw sibling")
    ap.add_argument("--selfref", action="store_true", help="also derive self-referenced offsets (cross-check) -> crosstalk_study/tcal/")
    ap.add_argument("--work", default=None)
    a = ap.parse_args()
    work = Path(a.work) if a.work else Path(tempfile.mkdtemp(prefix="sup-ct-"))
    a.cache_type = a.run_type + (f"_e{a.cluster_thr:g}" if a.cluster_thr else "")
    a.splits = a.splits or TYPE_SPLITS.get(a.run_type) or run_splits(a.run_type) or DEFAULT_SPLITS
    parts, tot = [], 0.0
    ep_sel = a.selection == "ep"
    from concurrent.futures import ThreadPoolExecutor
    def one(s):
        run, sp = (int(x) for x in s.split("/"))
        return load_or_make(run, sp, a.cache_type, a.path, a.reuse, work, ep_sel, a.cluster_thr)
    with ThreadPoolExecutor(a.jobs) as ex:
        for r, el in ex.map(one, a.splits):
            parts.append(r)
            tot += el or 0.0
    if a.time_source == "cfd" and a.cfd_algo == "fit":
        todo = [(HERE / "cache" / a.cache_type / f"{int(sp.split('/')[0]):06d}_{int(sp.split('/')[1]):05d}_{a.path}.npz",
                 int(sp.split("/")[0]), int(sp.split("/")[1])) for sp, r in zip(a.splits, parts) if "algo" not in r]
        if todo:
            from prad2py import dec as _dec
            if not hasattr(_dec.WaveAnalyzer, "analyze_fit"):
                sys.exit("[ERROR] prad2py lacks dec.WaveAnalyzer.analyze_fit (see crosstalk_study/README_prad2py_fit.md); "
                         "use --cfd-algo all for the any-algorithm CFD times")
            import multiprocessing as mp
            t0 = time.time()
            with mp.get_context("spawn").Pool(max(1, min(a.jobs, len(todo)))) as pool:
                for cf, nrow in pool.imap_unordered(wa_fit_split, todo):
                    print(f"[wa-fit] {cf.name}: {nrow:,} WaveAnalyzer peaks matched", flush=True)
            print(f"[wa-fit] {len(todo)} splits, {time.time() - t0:.0f} s", flush=True)
            parts = [dict(np.load(HERE / "cache" / a.cache_type / f"{int(sp.split('/')[0]):06d}_{int(sp.split('/')[1]):05d}_{a.path}.npz"))
                     for sp in a.splits]
    mnames = parts[0].get("module_names")
    assert all(np.array_equal(p["module_names"], mnames) for p in parts), "module name tables differ between splits"
    d = {k: np.concatenate([p[k] for p in parts]) for k in parts[0] if parts[0][k].ndim and k != "module_names"}
    d["part"] = np.repeat(np.arange(len(parts), dtype=np.int32), [len(p["t"]) for p in parts])      # split index (ev ids repeat across splits)
    thr_txt = f"{a.cluster_thr:g} MeV" if a.cluster_thr else "500 MeV x kE"
    ncl = sum(int(p["n_clusters"]) for p in parts)
    nm_info = None
    if a.nmod_max:
        for si, (sp, p_) in enumerate(zip(a.splits, parts)):
            if f"nmod_{a.nmod_def}" not in p_:
                run_, spl_ = (int(x) for x in sp.split("/"))
                cf = HERE / "cache" / a.cache_type / f"{run_:06d}_{spl_:05d}_{a.path}.npz"
                parts[si] = add_nmod(cf, work / f"{run_:06d}_{spl_:05d}" / "pulse.npz")
        d = {k: np.concatenate([p[k] for p in parts]) for k in parts[0] if parts[0][k].ndim and k != "module_names"}
        d["part"] = np.repeat(np.arange(len(parts), dtype=np.int32), [len(p["t"]) for p in parts])
        cn = d[f"cl_nmod_{a.nmod_def}"]
        nm_info = dict(n_before=int(len(cn)), median=float(np.median(cn)), p90=float(np.percentile(cn, 90)), frac_pass=float(np.mean(cn < a.nmod_max)))
        print(f"[nmod] definition {a.nmod_def}; single-cluster events before: {len(cn):,}; n_mod median {nm_info['median']:.0f}, p90 {nm_info['p90']:.0f}; "
              f"pass n_mod < {a.nmod_max}: {int((cn < a.nmod_max).sum()):,} ({100 * nm_info['frac_pass']:.1f} %)")
        hk, ck = d[f"nmod_{a.nmod_def}"] < a.nmod_max, cn < a.nmod_max
        d = {k: (v[hk] if len(v) == len(hk) else v[ck] if len(v) == len(ck) else v) for k, v in d.items()}
        ncl = int(ck.sum())
    sl = a.splits[0] if len(a.splits) == 1 else f"{a.splits[0]}..{a.splits[-1]}"
    head = f"{a.run_type}, {len(a.splits)} split(s) {sl}, {a.path} clusters\n"
    if ep_sel:
        sig = float(np.sum([float(p["sigma_cb"]) * int(p["n_inner"]) for p in parts]) / sum(int(p["n_inner"]) for p in parts))
        rw = 3 * sig
        cm = d["cl_moller"] == 0
        ce = cm & (np.abs(d["cl_r"] - 1) < rw)
        hm = (d["moller"] == 0) & (np.abs(d["r"] - 1) < rw)
        print(f"[cutflow] single-cluster events {ncl:,}")
        for lab, b in (("pair", 1), ("kin", 2), ("rej", 4)):
            print(f"[cutflow]   tag {lab}: {int(np.sum((d['cl_moller'] & b) > 0)):,}")
        print(f"[cutflow] not Moller-tagged (pair|kin|rej == 0): {int(cm.sum()):,}")
        print(f"[cutflow] sigma_CB/E = {sig:.5f} (chain result.tsv, n_inner-weighted), window |r-1| < {rw:.4f} "
              f"-> r in [{1 - rw:.4f}, {1 + rw:.4f}]: {int(ce.sum()):,}")
        nev = int(ce.sum())
        d = {k: (v[hm] if len(v) == len(hm) else v) for k, v in d.items()}
        ttl = (f"{head}e-p single-cluster events (no Møller tag, |r−1| < 3σ), N_ev = {nev:,}, {len(d['t']):,} non-seed hits")
    else:
        ttl = f"{head}single-cluster events (1 cluster > {thr_txt}), N_ev = {ncl:,}, {len(d['t']):,} non-seed hits"
    nmod = len(mnames)
    dsel = d
    off, off_label = load_offsets(a.tcal, mnames)

    def build(src):
        """-> (dict with t/seed_t of the time source, title suffix)."""
        x = dict(dsel)
        if src == "cfd":
            nh = len(x["t"])
            bad_s, bad_m = np.isnan(x["seed_t_wa"]), np.isnan(x["t_wa"])
            ok = ~(bad_s | bad_m)
            print(f"[cfd] hits {nh:,}: seed without WA match {int(bad_s.sum()):,} ({100 * bad_s.mean():.3f} %), module without WA "
                  f"match {int(bad_m.sum()):,} ({100 * bad_m.mean():.3f} %), kept {int(ok.sum()):,}; WA-MF time difference: median "
                  f"{np.nanmedian(x['t_wa'] - x['t']):.3f} ns, |diff| > 1 ns {100 * np.nanmean(np.abs(x['t_wa'] - x['t']) > 1):.2f} % "
                  f"(match window {WA_MATCH:g} ns)")
            if "algo" in x:
                for lab, al in (("seed", x["seed_algo"]), ("non-seed hit", x["algo"])):
                    c = [int(np.sum(al == v)) for v in (0, 1, 2)]
                    print(f"[algo] {lab}: time_algo 0 (peak) {100 * c[0] / len(al):.2f} %, 1 (linear CFD) {100 * c[1] / len(al):.2f} %, "
                          f"2 (log-normal fit) {100 * c[2] / len(al):.2f} %, unmatched {100 * np.mean(al < 0):.3f} % (N={len(al):,})")
                if a.cfd_algo == "fit":
                    ok = ok & (x["algo"] == 2) & (x["seed_algo"] == 2)
                    print(f"[algo] fit-only (time_algo 2 at seed and module): {int(ok.sum()):,} of {nh:,} hits ({100 * ok.mean():.2f} %)")
            x = {k: (v[ok] if len(v) == len(ok) else v) for k, v in x.items()}
            x["seed_t"], x["t"] = x["seed_t_wa"], x["t_wa"]
        return x

    def report(tag, s):
        print(f"[summary:{tag}] fraction |dt| > 10 ns = {100 * s['frac_out10']:.3f} %; dt > +10 ns (module later) {s['n_late10']:,}, "
              f"dt < -10 ns (module earlier) {s['n_early10']:,}; core robust sigma {s['w_robust']:.3f} ns, Gauss {s['w_gauss']:.3f} ns")
        print(f"[summary:{tag}] median chi2/dof |dt|<5: {s['chi_in']:.3f} (N={s['n_in']}); |dt|>20: {s['chi_out']:.3f} (N={s['n_out']}); "
              f"fig2 limit {s['chi_ymax']:.0f}, fraction above {100 * s['frac_chi_above']:.3f} %")

    if a.selfref:
        outd = HERE / "tcal"; outd.mkdir(exist_ok=True)
        for src in ("mf", "cfd"):
            x = build(src)
            r_ = selfref_offsets(x["t"], x["seed_t"], x["mod_i"].astype(np.int64), x["seed_mod_i"].astype(np.int64), nmod,
                                 official=load_offsets("auto", mnames)[0])
            f = outd / f"{a.run_type}_{src}_selfref_offsets.tsv"
            with open(f, "w") as fo:
                fo.write("module\tn\toffset_ns\tse\n")
                for m in np.nonzero(r_["n"] >= 50)[0]:
                    fo.write(f"{mnames[m]}\t{r_['n'][m]}\t{r_['off'][m]:.4f}\t{r_['se'][m]:.4f}\n")
            print(f"[selfref] {a.run_type} {src} ({a.selection}, {a.path}): {int((r_['n'] >= 50).sum())} modules, pair-network LSQ ({r_['passes']} lsqr iterations); RMS(selfref - time_calib.json, median removed) = "
                  f"{r_['rms_vs_official']:.3f} ns over {r_['n_cmp']} modules; RMS over seed-neighbour pairs of "
                  f"(o_m - o_s)(selfref) - (o_m - o_s)(official) = {r_['rms_pair']:.3f} ns -> {f}")
    x = build(a.time_source)
    srcl = ("matched filter" if a.time_source == "mf" else
            "WaveAnalyzer CFD, log-normal fit only (time_algo 2 at seed and module)" if a.cfd_algo == "fit" else
            "WaveAnalyzer CFD, any time_algo (0 peak / 1 linear / 2 log-normal fit)")
    ttl = ttl.replace(f"{len(dsel['t']):,} non-seed hits", f"{len(x['t']):,} non-seed hits")
    ttl += f"\ntime: {srcl}; amplitude, chi2/dof: MF"
    fitmode = a.time_source == "cfd" and a.cfd_algo == "fit"
    t_raw, ts_raw = x["t"], x["seed_t"]
    wr = None
    ct = tuple(a.ct_window)
    if off is not None:
        x["t"], x["seed_t"] = t_raw, ts_raw
        rawdir = Path(tempfile.mkdtemp(prefix="sup-ct-raw-")) if a.no_raw else Path(str(a.out) + "_raw")
        s0 = plot_all(x, rawdir, ttl + "\ntcal: off (raw)", a_algo_fit=fitmode, ct=ct)
        if a.no_raw:
            shutil.rmtree(rawdir, ignore_errors=True)
        report("raw", s0)
        wr = f"raw: robust sigma {s0['w_robust']:.3f} ns, Gauss {s0['w_gauss']:.3f} ns"
        x["t"], x["seed_t"] = t_raw - off[x["mod_i"]], ts_raw - off[x["seed_mod_i"]]
    nev_all = nev if ep_sel else ncl
    s = plot_all(x, Path(a.out), ttl + f"\ntcal: {off_label}", extra=wr or "", a_algo_fit=fitmode, ct=ct,
                 short=f"{a.run_type} e-p single-cluster events, N_ev = {nev_all:,}" if ep_sel else f"{a.run_type} single-cluster events, N_ev = {nev_all:,}")
    dtx = x["t"] - x["seed_t"]
    inw = (dtx >= ct[0]) & (dtx <= ct[1])
    ekey = x["part"].astype(np.int64) * (1 << 32) + x["ev"]
    inc = np.abs(dtx) < a.core_window
    n_ev_core = int(len(np.unique(ekey[inc])))
    n_ev_win = int(len(np.unique(ekey[inw])))
    n_ev_hit = int(len(np.unique(ekey)))
    row = dict(run_type=a.run_type, path=a.path, selection=a.selection, time_source=a.time_source, cfd_algo=a.cfd_algo, ct_lo=ct[0], ct_hi=ct[1],
               cluster_thr=a.cluster_thr or 0, nmod_max=a.nmod_max or 0, nmod_def=a.nmod_def if a.nmod_max else "",
               nmod_median=(nm_info or {}).get("median", 0), nmod_p90=(nm_info or {}).get("p90", 0), n_ev_before_nmod=(nm_info or {}).get("n_before", 0), n_ev=nev_all, n_ev_with_hits=n_ev_hit, n_ev_window=n_ev_win, n_hits=len(dtx), n_hits_window=int(inw.sum()),
               core=a.core_window, n_ev_core=n_ev_core, n_hits_core=int(inc.sum()))
    print(f"[ratio] {a.run_type}: events with >= 1 non-seed hit in {ct[0]:g}..{ct[1]:g} ns: {n_ev_win:,} / {nev_all:,} = {100 * n_ev_win / nev_all:.2f} %; "
          f"window hits {int(inw.sum()):,} / {len(dtx):,} = {100 * inw.mean():.2f} %; window hits per event {inw.sum() / nev_all:.3f} "
          f"({inw.sum() / max(n_ev_win, 1):.2f} per event with a window hit)")
    print(f"[ratio] {a.run_type}: |dt| < {a.core_window:g} ns: events with >= 1 hit {n_ev_core:,} / {nev_all:,} = {100 * n_ev_core / nev_all:.2f} %; "
          f"hits {int(inc.sum()):,} / {len(dtx):,} = {100 * inc.mean():.2f} %; hits per event {inc.sum() / nev_all:.3f}")
    write_ratio(row, off_label)
    print(f"[summary] {ttl}\n[summary] tcal: {off_label}")
    report("cal" if off is not None else "raw", s)
    print(f"[fig4] crosstalk window {ct[0]:g}..{ct[1]:g} ns: N = {s['ct_n']:,}; median non-seed amp {s['ct_med_amp']:.2f} ADC; "
          f"median non-seed/seed ratio {s['ct_med_ratio']:.4f}; fraction with non-seed amp < 10 ADC {100 * s['ct_frac10']:.2f} %")
    print(f"[fig4b] core |dt| < 5 ns: N = {s['core_n']:,}; median non-seed amp {s['core_med_amp']:.2f} ADC; median ratio "
          f"{s['core_med_ratio']:.4f}; fraction < 10 ADC {100 * s['core_frac10']:.2f} %")
    print(f"[summary] median non-seed/seed amplitude = {np.median(x['amp'] / x['seed_amp']):.4f}")
    print(f"[summary] benchmark wall time total {tot:.0f} s")
    if not a.work:
        shutil.rmtree(work, ignore_errors=True)


if __name__ == "__main__":
    main()
