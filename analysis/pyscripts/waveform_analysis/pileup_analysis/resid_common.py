"""resid_common.py -- shared pass of t1_residual_profile.py and t2_flat_top.py: the single-pulse waveform set S1, cached.

Pass (first N physics events of each EVIO split of the monoedge sample, 13 runs x 2 splits; splits and per-file noise references come from
output/pileup_rates/monoedge_sample/new/run_*.npz): trigger gate and LMS veto and the array-wide veto (>= VETO_NKEPT kept channels) as pileup_rates.scan,
fit with the DEFAULT pileup_rates fit (MONO_EDGE and NOISE_FLAG off, SHAPE_SCALE 1; forced below, whatever the environment says), every waveform of
the sample.  S1 = fit k == 1, no overflow sample (4095), fitted amplitude A in [A_LO, A_HI] ADC.  "Pruned" = waveform without a fitted pulse
(k == 0) although its raw peak (max - 30-sample pedestal) is in [A_LO, A_HI] ADC and has no overflow sample (reference count only).
Per S1 waveform the cache keeps: run, event, module, samples Y, pedestal, sigma_ref, A, t0, chi2/dof, residual Y - ped - fit, sigma(t)^2 variance.
Cache: <cache-dir>/resid_cache_N<N>_P<Npwo>.npz (default output/resid_study/ next to the other pileup outputs; not for git).
"""
from __future__ import annotations

import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")
for _v, _d in (("PILEUP_MONO_EDGE", "0"), ("PILEUP_NOISE_FLAG", "0"), ("PILEUP_SHAPE_SCALE", "1")):
    os.environ[_v] = _d                                  # default fit behaviour only

import sys
import time
import argparse
from pathlib import Path

HERE = Path(__file__).resolve().parent
PYDIR = HERE.parent.parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "exotic_events"))
import numpy as np                                      # noqa: E402

EVIO_DIR = "/Users/L00338853/work/PRad/data/evio"
SAMPLE_DIR = PYDIR / "output" / "pileup_rates" / "monoedge_sample"
CACHE_DIR = PYDIR / "output" / "resid_study"
A_LO, A_HI = 50.0, 3000.0
PREFILTER = 40.0                                         # ADC: raw peak (max - median of the 30 pedestal samples) needed to fit a waveform in this pass
AMP_EDGES = np.array([50.0, 200.0, 800.0, 3000.0])
AMP_LABELS = ["50-200", "200-800", "800-3000"]
CLASSES = ["PbWO4", "LG"]
NS, DT = 100, 4.0
TT = DT * np.arange(NS)
N_DEFAULT, N_PWO_DEFAULT = 60000, 1500
STAMP = "20261116"


def add_common_args(ap):
    ap.add_argument("--n-events", type=int, default=N_DEFAULT, help="first N decoded events of each EVIO split (counted before the trigger gate), used for lead glass (default %(default)s); names the cache")
    ap.add_argument("--n-events-pwo", type=int, default=N_PWO_DEFAULT, help="first N decoded events of each split (before the trigger gate) used for PbWO4 (plentiful; default %(default)s); names the cache")
    ap.add_argument("--cache-dir", default=str(CACHE_DIR), help="directory of the cache npz (default %(default)s)")
    ap.add_argument("--rebuild", action="store_true", help="redo the pass even if the cache exists")
    ap.add_argument("--nproc", type=int, default=16, help="worker processes of the pass (default %(default)s)")
    ap.add_argument("--sample-dir", default=str(SAMPLE_DIR), help="monoedge sample scans (splits, noise references) (default %(default)s)")
    ap.add_argument("--evio-dir", default=EVIO_DIR)
    ap.add_argument("--out", default=str(HERE), help="output directory of PDF and TSV (default %(default)s)")


def _task(job):
    """One EVIO split: returns dict of S1 arrays + counts.  Imports of the fit code are here so that cache users need no prad2py."""
    run, si, path, b_ref, sig_ref, nmax, n_pwo = job
    import event_display as ed
    from prad2py import dec
    pr = ed.pr
    assert not pr.MONO_EDGE and not pr.NOISE_FLAG and pr.SHAPE_SCALE == 1.0
    t_start = time.time()
    ctx = ed.Ctx(run)
    rt = ctx.rt
    gate = (lambda b: (b & pr.PHYSICS_MASK_ANY) != 0) if rt in pr.PHYSICS_ANY else (lambda b: b == pr.PHYSICS_BITS_EXACT)
    nm = len(ctx.names)
    S = {k: [] for k in ("ev", "mi", "Y", "ped", "sig", "a", "t0", "chi2", "resid", "var")}
    pruned = np.zeros(nm, np.int64)
    cnt = dict(n_seen=0, n_pass=0, n_dig=0, n_fit=0, n_veto=0)
    rows, evs, mis, evnum = [], [], [], []

    def flush():
        if not rows:
            return
        Y, ev, mi = np.stack(rows), np.array(evs), np.array(mis)
        rows.clear(), evs.clear(), mis.clear()
        x = np.median(Y[:, :pr.PED_NSAMPLES], axis=1)
        raw0 = Y.max(axis=1) - x
        cnt["n_dig"] += len(mi)
        nk = np.bincount(ev[raw0 >= pr.KEEP_ADC], minlength=ev.max() + 1)[ev]           # proxy of the kept-channel count of the array-wide veto
        ok = nk < pr.VETO_NKEPT
        cnt["n_veto"] += int((~ok & (raw0 >= PREFILTER)).sum())
        # only waveforms that can be S1 or pruned (raw peak >= PREFILTER ADC) are fitted; PbWO4 only in the first n_pwo events (the fit of a waveform does not depend on the others)
        pre = (raw0 >= PREFILTER) & ok & (~ctx.pwo4[mi] | (ev < n_pwo))
        s = np.nonzero(pre)[0]
        s = s[~pr.bad_waveform(Y[s])]
        if s.size == 0:
            return
        m = ctx.bank.row[mi[s]]
        f = pr.fit(ctx.mf, Y[s], m, sig_ref[mi[s]], b_ref[mi[s]])
        f.pop("far")
        cnt["n_fit"] += s.size
        ovf = (Y[s] >= pr.OVERFLOW).any(axis=1)
        raw = Y[s].max(axis=1) - f["ped"]
        a1 = f["a"][:, 0]                                    # pulses are in time order; k == 1: the only one
        s1 = (f["k"] == 1) & ~ovf & (a1 >= A_LO) & (a1 <= A_HI)
        pr_ = (f["k"] == 0) & ~ovf & (raw >= A_LO) & (raw <= A_HI)
        pruned[:] += np.bincount(mi[s][pr_], minlength=nm)
        idx = np.nonzero(s1)[0]
        if idx.size == 0:
            return
        k = f["k"][idx]
        valid = np.arange(pr.KMAX)[None, :] < k[:, None]
        t0 = np.where(valid, np.nan_to_num(f["t"][idx]), 0.0)
        a = f["a"][idx]
        noise = np.sqrt(f["sig"][idx].astype(np.float64) ** 2 + 1.0 / 12.0)
        var = ctx.mf.variance(m[idx], t0.astype(np.float32), a, valid, noise, NS)
        Ys = Y[s[idx]]
        res = Ys - f["ped"][idx][:, None] - f["comp"][idx].sum(axis=1)
        for key, v in (("ev", np.array(evnum)[ev[s[idx]]]), ("mi", mi[s[idx]]), ("Y", Ys), ("ped", f["ped"][idx]), ("sig", f["sig"][idx]), ("a", a[:, 0]),
                       ("t0", t0[:, 0]), ("chi2", f["chi2"][idx]), ("resid", res), ("var", var)):
            S[key].append(np.asarray(v))

    ch = ctx.open(path)
    done = False
    while not done and ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            cnt["n_seen"] += 1
            done = cnt["n_seen"] >= nmax
            fe = d["event"]
            bits = int(fe.info.trigger_bits)
            if gate(bits) and not (bits & pr.LMS_BIT):
                cnt["n_pass"] += 1
                ev = len(evnum)
                evnum.append(int(fe.info.event_number))
                for mi, smp in ctx.channels(fe):
                    rows.append(smp)
                    evs.append(ev)
                    mis.append(mi)
                if len(evnum) % pr.BLOCK == 0:
                    flush()
            if done:
                break
    flush()
    ch.close()
    out = {k: (np.concatenate(v) if v else np.zeros((0, NS) if k in ("Y", "resid", "var") else 0)) for k, v in S.items()}
    out["run"] = np.full(out["mi"].size, run, np.int32)
    out["split"] = np.full(out["mi"].size, si, np.int8)
    out["rt"] = rt
    out["names"] = ctx.names
    out["pwo4"] = ctx.pwo4
    out["pruned"] = pruned
    out["cnt"] = cnt
    out["wall"] = time.time() - t_start
    return out


def build(a):
    """Run the pass, save and return the cache dict."""
    import multiprocessing as mp
    jobs = []
    for fn in sorted((Path(a.sample_dir) / "new").glob("run_*.npz")):
        z = np.load(fn)
        run = int(z["run"])
        for i, sp in enumerate(z["splits"]):
            jobs.append((run, i, str(Path(a.evio_dir) / f"prad_{run:06d}.evio.{int(sp):05d}"), z["mod_b_ref"][i], z["mod_sig_ref"][i], a.n_events, a.n_events_pwo))
    jobs.sort(key=lambda j: -j[0] % 7)                   # mix run types across the pool
    t0 = time.time()
    print(f"[pass] {len(jobs)} splits, first {a.n_events} events each, {a.nproc} processes", flush=True)
    res = []
    with mp.get_context("spawn").Pool(min(a.nproc, len(jobs))) as pool:
        for r in pool.imap_unordered(_task, jobs):
            res.append(r)
            print(f"[pass] run {r['run'][0] if r['run'].size else '?'} {r['rt']}: {r['cnt']} S1 {r['mi'].size}  {r['wall']:.0f} s  ({len(res)}/{len(jobs)})", flush=True)
    names = res[0]["names"]
    assert all(np.array_equal(r["names"], names) for r in res)
    cat = lambda k: np.concatenate([r[k] for r in res])   # noqa: E731
    out = {k: cat(k).astype(dt) for k, dt in (("run", np.int32), ("split", np.int8), ("ev", np.int64), ("mi", np.int16), ("Y", np.uint16), ("ped", np.float32), ("sig", np.float32),
                                                ("a", np.float32), ("t0", np.float32), ("chi2", np.float32), ("resid", np.float32), ("var", np.float32))}
    out["rt"] = np.concatenate([[r["rt"]] * r["mi"].size for r in res]).astype(str)
    out["names"], out["pwo4"] = names, res[0]["pwo4"]
    out["pruned"] = sum(r["pruned"] for r in res)
    for k in ("n_seen", "n_pass", "n_dig", "n_fit", "n_veto"):
        out[k] = np.int64(sum(r["cnt"][k] for r in res))
    out["n_events"], out["n_events_pwo"], out["wall"] = np.int64(a.n_events), np.int64(a.n_events_pwo), np.float64(time.time() - t0)
    p = cache_path(a)
    os.makedirs(p.parent, exist_ok=True)
    np.savez_compressed(p, **out)
    print(f"[pass] saved {p}  {out['mi'].size} S1 waveforms, {time.time() - t0:.0f} s wall", flush=True)
    return out


def cache_path(a):
    return Path(a.cache_dir) / f"resid_cache_N{a.n_events}_P{a.n_events_pwo}.npz"


def load_s1(a):
    """S1 cache dict (arrays), built first if needed.  Adds cls (0 PbWO4, 1 LG), name per waveform, ab (amplitude bin 0..2)."""
    p = cache_path(a)
    if a.rebuild or not p.exists():
        d = build(a)
    else:
        d = dict(np.load(p))
        print(f"[cache] {p}: {d['mi'].size} S1 waveforms", flush=True)
    # S1 cut at load time: the fit window (MF.wmask: 32 samples from ceil(t0/4) - 3) must lie inside samples 0..NS-1
    w0 = np.ceil(d["t0"] / DT).astype(np.int64) - 3
    inside = (w0 >= 0) & (w0 + 32 <= NS)
    assert np.array_equal(inside, window_mask(d["t0"]).sum(axis=1) == 32)
    cls0 = np.where(d["pwo4"][d["mi"]], 0, 1)
    d["n_trunc"] = np.array([int(((cls0 == c) & ~inside).sum()) for c in range(2)])
    n0 = inside.size
    for k, v in list(d.items()):
        if isinstance(v, np.ndarray) and v.ndim >= 1 and v.shape[0] == n0 and k not in ("names", "pwo4", "pruned"):
            d[k] = v[inside]
    print(f"[S1] removed truncated-window waveforms: PbWO4 {d['n_trunc'][0]}, LG {d['n_trunc'][1]}", flush=True)
    d["cls"] = np.where(d["pwo4"][d["mi"]], 0, 1).astype(np.int8)
    d["ab"] = (np.digitize(d["a"], AMP_EDGES[1:-1])).astype(np.int8)
    d["modname"] = d["names"][d["mi"]]
    return d


def counts_text(d):
    """S1 counts per class and the pruned-candidate counts."""
    lines = []
    for c, cn in enumerate(CLASSES):
        sel = d["cls"] == c
        pru = int(d["pruned"][(d["pwo4"] == (c == 0))].sum())
        lines.append(f"{cn}: S1 {int(sel.sum())} after removing {int(d['n_trunc'][c])} with the fit window outside the waveform (A bins {AMP_LABELS[0]}/{AMP_LABELS[1]}/{AMP_LABELS[2]}: "
                     + "/".join(str(int((sel & (d['ab'] == b)).sum())) for b in range(3)) + f"), pruned-candidate waveforms {pru}")
    return lines


def window_mask(t0):
    """Fit window of a pulse with onset t0 (pileup_rates.MF.wmask: 32 samples starting 3 before the onset) [n, NS]."""
    w0 = (np.ceil(t0 / DT).astype(np.int64) - 3)[:, None]
    ii = np.arange(NS)[None, :]
    return (ii >= w0) & (ii < w0 + 32)


def run_runs_note(d):
    runs = sorted(set(int(x) for x in d["run"]))
    return f"N = {int(d['n_events'])} first events per EVIO split for lead glass ({int(d['n_events_pwo'])} for PbWO4), {len(runs)} runs x 2 splits (monoedge sample); {int(d['n_pass'])} events after trigger/LMS gate, {int(d['n_veto'])} fitted-candidate rows lost to the array-wide veto"


def panel(ax, Y, ped, model, t0=None, title=None):
    """Data (black dots), pedestal (grey), total fit (red) of one waveform."""
    ax.plot(TT, Y, "k.", ms=3.5)
    ax.axhline(ped, color="0.6", lw=0.8)
    ax.plot(TT, ped + model, "r-", lw=1.0)
    if t0 is not None:
        ax.axvline(t0, color="tab:blue", ls=":", lw=0.8)
    if title:
        ax.set_title(title, fontsize=8, loc="left", pad=2)
    ax.tick_params(labelsize=7)
    ax.set_xlim(0, TT[-1])


def main():
    ap = argparse.ArgumentParser(description="Build (or check) the S1 cache shared by t1_residual_profile.py and t2_flat_top.py.")
    add_common_args(ap)
    a = ap.parse_args()
    d = load_s1(a)
    print(run_runs_note(d))
    print("\n".join(counts_text(d)))


if __name__ == "__main__":
    main()


# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# T1x extension (high-amplitude pass).  The S1 cache above stops at A = 3000 ADC and has no overflow waveforms; this pass selects PbWO4 single-pulse rows with
# A >= HI_LO from the stored scans of the 2-file sample (sample-dir/old/run_*.npz: ch_k == 1, ch_a[:, 0], ch_ev -> ev_number, ch_midx), draws at most
# --hi-cap rows per (run type, amplitude bin 2000-3000 / 3000-4095), decodes only those events and refits their waveforms with the same default fit.
# Waveforms with overflow samples (nclip > 0) are kept and flagged ('ovf'); the analysis masks their overflow samples.
# Cache: <cache-dir>/resid_hi_cache_C<cap>.npz.  The functions above are not changed.
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
HI_LO, HI_MID, HI_TOP = 2000.0, 3000.0, 4096.0
HI_CAP_DEFAULT = 2500
HI_SEED = 20261116


def add_hi_args(ap):
    ap.add_argument("--hi-cap", type=int, default=HI_CAP_DEFAULT, help="max rows per (run type, amplitude bin 2000-3000 / 3000-4095) drawn from the sample scans for the high-amplitude pass (default %(default)s)")
    ap.add_argument("--rebuild-hi", action="store_true", help="redo the high-amplitude pass even if its cache exists")


def hi_cache_path(a):
    return Path(a.cache_dir) / f"resid_hi_cache_C{a.hi_cap}.npz"


def _task_hi(job):
    """One EVIO split: decode the wanted events (event number -> module indices), fit, return the rows with fitted k == 1, A >= HI_LO."""
    run, si, path, b_ref, sig_ref, want = job
    import event_display as ed
    from prad2py import dec
    pr = ed.pr
    assert not pr.MONO_EDGE and not pr.NOISE_FLAG and pr.SHAPE_SCALE == 1.0
    t_start = time.time()
    ctx = ed.Ctx(run)
    ch = ctx.open(path)
    rows, mis, evnum = [], [], []
    left = len(want)
    while left > 0 and ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            en = int(ch.info().event_number) if ch.get_n_events() == 1 else None
            if en is None:
                d = ch.decode_event(i, with_ssp=False)
                if not d["ok"]:
                    continue
                en = int(d["event"].info.event_number)
            if en not in want:
                continue
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            need = set(want[en])
            for mi, smp in ctx.channels(d["event"]):
                if mi in need:
                    rows.append(smp)
                    mis.append(mi)
                    evnum.append(en)
            left -= 1
    ch.close()
    out = {k: [] for k in ("ev", "mi", "Y", "ped", "sig", "a", "t0", "chi2", "resid", "var", "ovf")}
    n_fit = 0
    if rows:
        Y, mi, ev = np.stack(rows), np.array(mis), np.array(evnum, np.int64)
        s = np.nonzero(~pr.bad_waveform(Y))[0]
        if s.size:
            m = ctx.bank.row[mi[s]]
            f = pr.fit(ctx.mf, Y[s], m, sig_ref[mi[s]], b_ref[mi[s]])
            f.pop("far")
            n_fit = s.size
            a1 = f["a"][:, 0]
            idx = np.nonzero((f["k"] == 1) & (a1 >= HI_LO) & (a1 < HI_TOP))[0]
            if idx.size:
                valid = np.arange(pr.KMAX)[None, :] < f["k"][idx][:, None]
                t0 = np.where(valid, np.nan_to_num(f["t"][idx]), 0.0)
                a = f["a"][idx]
                noise = np.sqrt(f["sig"][idx].astype(np.float64) ** 2 + 1.0 / 12.0)
                var = ctx.mf.variance(m[idx], t0.astype(np.float32), a, valid, noise, NS)
                Ys = Y[s[idx]]
                res = Ys - f["ped"][idx][:, None] - f["comp"][idx].sum(axis=1)
                for key, v in (("ev", ev[s[idx]]), ("mi", mi[s[idx]]), ("Y", Ys), ("ped", f["ped"][idx]), ("sig", f["sig"][idx]), ("a", a[:, 0]), ("t0", t0[:, 0]),
                               ("chi2", f["chi2"][idx]), ("resid", res), ("var", var), ("ovf", (Ys >= pr.OVERFLOW).any(axis=1))):
                    out[key].append(np.asarray(v))
    out = {k: (np.concatenate(v) if v else np.zeros((0, NS) if k in ("Y", "resid", "var") else 0)) for k, v in out.items()}
    out["run"] = np.full(out["mi"].size, run, np.int32)
    out["rt"] = ctx.rt
    out["n_want"], out["n_found_ev"], out["n_rows"], out["n_fit"] = len(want), len(want) - left, len(rows), n_fit
    out["wall"] = time.time() - t_start
    return out


def build_hi(a):
    """Select rows from the sample scans, run the pass, save and return the cache dict."""
    import multiprocessing as mp
    rng = np.random.default_rng(HI_SEED)
    sc = {}
    pool = {}                                           # (rt, bin) -> list of (run, ev_number, mi)
    for fn in sorted((Path(a.sample_dir) / "old").glob("run_*.npz")):
        z = np.load(fn)
        run, rt = int(z["run"]), str(z["run_type"])
        sc[run] = (z["splits"], z["mod_b_ref"], z["mod_sig_ref"])
        A = z["ch_a"][:, 0]
        sel = np.nonzero((z["ch_k"] == 1) & (A >= HI_LO) & (A < HI_TOP) & z["mod_pwo4"][z["ch_midx"]])[0]
        for b, (lo, hi) in enumerate(((HI_LO, HI_MID), (HI_MID, HI_TOP))):
            s = sel[(A[sel] >= lo) & (A[sel] < hi)]
            pool.setdefault((rt, b), []).append(np.stack([np.full(s.size, run), z["ev_number"][z["ch_ev"][s]], z["ch_midx"][s].astype(np.int64)], axis=1))
    want = {}                                           # run -> {ev_number: [mi]}
    npool, npick = {}, {}
    for key, v in sorted(pool.items()):
        v = np.concatenate(v)
        npool[key] = len(v)
        if len(v) > a.hi_cap:
            v = v[rng.choice(len(v), a.hi_cap, replace=False)]
        npick[key] = len(v)
        for run, en, mi in v:
            want.setdefault(int(run), {}).setdefault(int(en), []).append(int(mi))
    print("[hi] pool / picked rows per (run type, bin):", {f"{k[0]}:{('2000-3000', '3000-4095')[k[1]]}": (npool[k], npick[k]) for k in sorted(npool)}, flush=True)
    jobs = []
    for run, (splits, b_ref, sig_ref) in sc.items():
        if run not in want:
            continue
        for i, sp in enumerate(splits):
            jobs.append((run, i, str(Path(a.evio_dir) / f"prad_{run:06d}.evio.{int(sp):05d}"), b_ref[i], sig_ref[i], want[run]))
    jobs.sort(key=lambda j: -len(j[5]))
    t0 = time.time()
    print(f"[hi] {len(jobs)} splits, {sum(len(j[5]) for j in jobs)} (event, split) lookups, {a.nproc} processes", flush=True)
    res = []
    with mp.get_context("spawn").Pool(min(a.nproc, len(jobs))) as pool_:
        for r in pool_.imap_unordered(_task_hi, jobs):
            res.append(r)
            print(f"[hi] run {r['run'][0] if r['run'].size else '?'} {r['rt']}: events found {r['n_found_ev']}/{r['n_want']}, rows {r['n_rows']}, S1 {r['mi'].size}  {r['wall']:.0f} s  ({len(res)}/{len(jobs)})", flush=True)
    cat = lambda k: np.concatenate([r[k] for r in res])   # noqa: E731
    out = {k: cat(k).astype(dt) for k, dt in (("run", np.int32), ("ev", np.int64), ("mi", np.int16), ("Y", np.uint16), ("ped", np.float32), ("sig", np.float32), ("a", np.float32),
                                                ("t0", np.float32), ("chi2", np.float32), ("resid", np.float32), ("var", np.float32), ("ovf", bool))}
    out["rt"] = np.concatenate([[r["rt"]] * r["mi"].size for r in res]).astype(str)
    out["pool_keys"] = np.array([f"{k[0]}:{k[1]}" for k in sorted(npool)])
    out["pool_n"] = np.array([npool[k] for k in sorted(npool)])
    out["pick_n"] = np.array([npick[k] for k in sorted(npool)])
    out["hi_cap"], out["wall"] = np.int64(a.hi_cap), np.float64(time.time() - t0)
    p = hi_cache_path(a)
    os.makedirs(p.parent, exist_ok=True)
    np.savez_compressed(p, **out)
    print(f"[hi] saved {p}  {out['mi'].size} waveforms, {time.time() - t0:.0f} s wall", flush=True)
    return out


def load_hi(a):
    """High-amplitude cache dict, built first if needed; truncated fit windows removed as in load_s1."""
    p = hi_cache_path(a)
    if a.rebuild_hi or not p.exists():
        d = build_hi(a)
    else:
        d = dict(np.load(p))
        print(f"[cache] {p}: {d['mi'].size} high-amplitude waveforms", flush=True)
    w0 = np.ceil(d["t0"] / DT).astype(np.int64) - 3
    inside = (w0 >= 0) & (w0 + 32 <= NS)
    n0 = inside.size
    d["n_trunc"] = int((~inside).sum())
    for k, v in list(d.items()):
        if isinstance(v, np.ndarray) and v.ndim >= 1 and v.shape[0] == n0 and k not in ("names", "pwo4", "pool_keys", "pool_n", "pick_n"):
            d[k] = v[inside]
    return d
