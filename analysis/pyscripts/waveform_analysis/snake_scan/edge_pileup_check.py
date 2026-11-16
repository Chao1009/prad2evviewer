#!/usr/bin/env python3
"""
edge_pileup_check.py - P14.2b: all-event cross-check of the 3x3 endpoint edge against the pile-up rejection of edge_pass.

  collect  <evio> --dwell DWELL --cache DIR     one decode per file; per dwell module and event keep, for each of the 9
                                                 members, the number of peaks and up to NP=3 peaks (time, height, integral)
  analyze  --cache DIR --out TSV                refit the edge for several rejection settings from the caches

Settings (SETTINGS below): `none` (no rejection, member value = sum of all peak integrals as edge_pass does for
>= 2 peaks), `w1` (edge_pass: reject if any fired member has >= 2 peaks; member value = its integral),
and looser ones that reject only if a secondary peak is above a height threshold and/or inside a time window.
Member value of a kept event = integral of the highest peak (secondary peaks below the threshold are noise, not added).
Imports fit_edge, make_hist, constants from edge_pass; edge_pass.py itself is not modified.
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np

SS = Path(__file__).resolve().parent
if str(SS) not in sys.path:
    sys.path.insert(0, str(SS))
import snake_scan as SSN  # noqa: E402
import edge_pass as EP  # noqa: E402

NP = 3


def collect(evio, rows, map_path=SSN.DERIVED_MAP, max_events=0):
    from prad2py import dec
    fname = Path(evio).name
    ranges = sorted((r["ev_first"], r["ev_last"], r["module"]) for r in rows if r["kind"] == "dwell" and r["file"] == fname)
    nb = EP.neighbours(map_path)
    pipe = SSN.open_pipeline(map_path, str(evio))
    hy, wa = pipe.hycal, pipe.wave_ana
    ch = SSN.open_evio(evio, pipe)
    cache = {}

    def modname(crate, slot, c):
        k = (crate, slot, c)
        if k not in cache:
            m = hy.module_by_daq(crate, slot, c)
            cache[k] = m.name if (m is not None and m.is_pwo4()) else None
        return cache[k]

    acc = {m: dict(members=nb[m], npk=[], pk=[], npile_noc=0) for _, _, m in ranges}
    ri_cur, nsc, done = 0, 0, False
    while not done and ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            ch.select_event(i)
            inf = ch.info()
            evn = int(inf.event_number)
            while ri_cur < len(ranges) and evn > ranges[ri_cur][1]:
                ri_cur += 1
            if ri_cur >= len(ranges):
                done = True
                break
            if evn < ranges[ri_cur][0] or int(inf.trigger_bits) & (SSN.TBIT_LMS | SSN.TBIT_ALPHA):
                continue
            nsc += 1
            mod = ranges[ri_cur][2]
            members = acc[mod]["members"]
            pos = {n: k for k, n in enumerate(members)}
            npk = np.zeros(len(members), np.uint8)
            pk = np.zeros((len(members), NP, 3), np.float32)      # time, height, integral
            central_ok = False
            e = ch.fadc()
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
                        nm = modname(crate, s, c)
                        k = pos.get(nm)
                        if k is None:
                            continue
                        cd = slot.channel(c)
                        if cd.nsamples < 9:
                            continue
                        _, _, peaks = wa.analyze(np.asarray(cd.samples, np.uint16))
                        npk[k] = min(len(peaks), 255)
                        for j, p in enumerate(peaks[:NP]):
                            pk[k, j] = (p.time, p.height, p.integral)
                        if len(peaks) > NP:     # keep the NP highest
                            top = sorted(peaks, key=lambda p: -p.height)[:NP]
                            for j, p in enumerate(top):
                                pk[k, j] = (p.time, p.height, p.integral)
                        if peaks and nm == mod:
                            central_ok = True
            if not central_ok:
                continue
            acc[mod]["npk"].append(npk)
            acc[mod]["pk"].append(pk)
            if max_events and nsc >= max_events:
                done = True
                break
    out = {}
    for m, a in acc.items():
        out[m] = dict(members=np.array(a["members"]), npk=np.array(a["npk"], np.uint8).reshape(-1, len(a["members"])),
                      pk=np.array(a["pk"], np.float32).reshape(-1, len(a["members"]), NP, 3))
    return out


def _collect_one(a):
    evio, rows, map_path, cache, mx = a
    t0 = time.process_time()
    res = collect(evio, rows, map_path, mx)
    for m, d in res.items():
        np.savez_compressed(Path(cache) / f"{Path(evio).name}_{m}.npz", **d)
    return Path(evio).name, {m: len(d["npk"]) for m, d in res.items()}, time.process_time() - t0


def cmd_collect(args):
    rows = SSN.read_dwell_tsv(args.dwell)
    Path(args.cache).mkdir(parents=True, exist_ok=True)
    jobs = [(f, rows, args.map, args.cache, args.max_events) for f in args.evio
            if any(r["file"] == Path(f).name and r["kind"] == "dwell" for r in rows)]
    if args.jobs > 1:
        from multiprocessing import Pool
        with Pool(min(args.jobs, 2)) as p:
            res = p.map(_collect_one, jobs)
    else:
        res = [_collect_one(j) for j in jobs]
    for r in res:
        print(f"[collect] {r[0]}: {r[1]} cpu {r[2]:.0f} s", flush=True)


# ---------------------------------------------------------------------------------------------------------------
# settings: reject(npk, pk, c) -> bool mask (True = pile-up).  c = index of the central module among the members.
# secondary peak of a member = every stored peak except the highest one.
# ---------------------------------------------------------------------------------------------------------------

T0 = [154.5]      # in-time peak position (ns), set from the data (median time of the central highest peak)


def _main_idx(pk):
    """index of the in-time (main) peak of each member: the stored peak closest to T0 (height > 0)."""
    d = np.where(pk[..., 1] > 0, np.abs(pk[..., 0] - T0[0]), np.inf)
    return np.argmin(d, axis=2)


def _secondary(npk, pk):
    """(sec_h, sec_dt): height and |time - main time| of every stored peak except the in-time one (0 / inf where absent)."""
    im = _main_idx(pk)
    tm = np.take_along_axis(pk[..., 0], im[..., None], 2)
    sec = (pk[..., 1] > 0) & (np.arange(NP)[None, None, :] != im[..., None])
    return np.where(sec, pk[..., 1], 0.0), np.where(sec, np.abs(pk[..., 0] - tm), np.inf)


def make_reject(kind, thr=0.0, tw=np.inf, scope="any"):
    """kind none: nothing rejected.  w1: npk >= 2 in any fired member (edge_pass).  sec: a secondary peak with
    height >= thr (ADC) and |dt to the in-time peak| <= tw (ns), in the central module only (scope 'central') or any member."""
    def f(npk, pk, c):
        if kind == "none":
            return np.zeros(len(npk), bool)
        if kind == "w1":
            return (npk >= 2).any(1)
        sh, sdt = _secondary(npk, pk)
        bad = (sh >= thr) & (sh > 0) & (sdt <= tw)
        if scope == "central":
            bad = bad[:, c]
        return bad.reshape(len(bad), -1).any(1)
    return f


def member_values(npk, pk, mode):
    """mode 'main': integral of the in-time peak; 'sum': sum of all stored peak integrals (edge_pass for >= 2 peaks)."""
    if mode == "sum":
        return pk[..., 2].sum(2)
    return np.take_along_axis(pk[..., 2], _main_idx(pk)[..., None], 2)[..., 0]


def load_events(cache):
    ev = {}
    for f in sorted(Path(cache).glob("*.npz")):
        m = f.name.split("_")[-1][:-4]
        d = np.load(f)
        a = ev.setdefault(m, dict(members=list(d["members"]), npk=[], pk=[]))
        a["npk"].append(d["npk"])
        a["pk"].append(d["pk"])
    for a in ev.values():
        a["npk"] = np.concatenate(a["npk"])
        a["pk"] = np.concatenate(a["pk"])
    return ev


def hist_for(a, m, W, rej, mode, extra_scale=None):
    c = a["members"].index(m)
    mask = rej(a["npk"], a["pk"], c)
    I = member_values(a["npk"], a["pk"], mode)[~mask]
    return EP.make_hist(dict(module=m, members=a["members"], I=I, n_pileup=int(mask.sum())), W, extra_scale)


def fit_setting(ev, W, rej, mode):
    """Two-pass fit as edge_pass.run: pass 1 with W, pass 2 with neighbours' W_j / s_j."""
    h1 = {m: hist_for(a, m, W, rej, mode) for m, a in ev.items()}
    f1 = {m: EP.fit_edge(h) for m, h in h1.items()}
    s1 = {m: v[0] / EP.E_BEAM for m, v in f1.items()}
    out = {}
    for m, a in ev.items():
        scale = {n: s1[n] for n in a["members"] if n != m and n in s1}
        h = hist_for(a, m, W, rej, mode, scale)
        E, sig, pl, sl, err, spread = EP.fit_edge(h)
        s = E / EP.E_BEAM
        out[m] = dict(E=E, s=s, s_err=float(np.hypot(err / EP.E_BEAM, 0.5 * spread * s)), spread=spread, sigma=sig,
                      n_keep=h["n_single"], n_rej=h["n_pileup"], s_pass1=s1[m])
    return out


# (name, reject, member-value mode, description)
SETTINGS = [
    ("none_sum", make_reject("none"), "sum", "no rejection, member value = sum of all peak integrals"),
    ("none_main", make_reject("none"), "main", "no rejection, member value = in-time peak integral"),
    ("w1", make_reject("w1"), "main", "wave 1: reject if any fired member has >= 2 peaks"),
    ("loose1", make_reject("sec", 50.0, 100.0, "any"), "main", "reject if a secondary peak >= 50 ADC lies within 100 ns of the in-time peak (any member)"),
    ("loose2", make_reject("sec", 200.0, 60.0, "any"), "main", "reject if a secondary peak >= 200 ADC lies within 60 ns of the in-time peak (any member)"),
]


def cmd_analyze(args):
    import math
    import transfer as TR
    ev = load_events(args.cache)
    c0 = []
    for m, a in ev.items():
        c = a["members"].index(m)
        c0.append(np.median(a["pk"][np.arange(len(a["pk"])), c, np.argmax(a["pk"][:, c, :, 1], 1), 0]))
    T0[0] = float(np.median(c0))
    print(f"[analyze] T0 = {T0[0]:.1f} ns, modules {sorted(ev)}")
    W = EP.load_calib(EP.resolve_calib(23701))
    res = {}
    for name, rej, mode, desc in SETTINGS:
        t0 = time.process_time()
        res[name] = fit_setting(ev, W, rej, mode)
        print(f"[analyze] {name}: " + " ".join(f"{m}:{v['s']:.4f}" for m, v in sorted(res[name].items())) +
              f" ({time.process_time() - t0:.0f} s)", flush=True)
    mods = sorted(ev, key=lambda n: int(n[1:]))
    names = [n for n, *_ in SETTINGS]
    cols = ["module", "setting", "E_edge_MeV", "s", "s_err", "window_spread", "sigma_MeV", "n_keep", "n_rej", "frac_rej",
            "dE_vs_w1_MeV", "ds_vs_w1_pct", "sys_spread_pct", "sys_shift_none_pct", "s_err_w1_pct"]
    rows = []
    for m in mods:
        sw = np.array([res[n][m]["s"] for n in names])
        w1 = res["w1"][m]
        spread = 100 * (sw.max() - sw.min()) / w1["s"]
        shift = 100 * (res["none_sum"][m]["s"] / w1["s"] - 1)
        for n in names:
            r = res[n][m]
            rows.append([m, n, f"{r['E']:.1f}", f"{r['s']:.5f}", f"{r['s_err']:.5f}", f"{r['spread']:.5f}", f"{r['sigma']:.1f}",
                         r["n_keep"], r["n_rej"], f"{r['n_rej'] / (r['n_keep'] + r['n_rej']):.4f}",
                         f"{r['E'] - w1['E']:+.1f}", f"{100 * (r['s'] / w1['s'] - 1):+.3f}", f"{spread:.3f}", f"{shift:+.3f}",
                         f"{100 * w1['s_err'] / w1['s']:.3f}"])
    # transfer at 2.24 and 3.49 GeV with each setting's s (L_scan primary)
    reg = TR.region_of()
    ls = TR.read_L(TR.OUT / "lms_scan.tsv")
    trf = []
    for rt in ("2p24GeV_H2", "3p49GeV_H2"):
        lp = TR.read_L(TR.OUT / f"lms_prod_{rt}.tsv")
        gp = TR.read_pool_gains(TR.V3 / rt / "ring" / "gains.tsv")
        for n in names:
            edge = {m: (res[n][m]["s"], res[n][m]["s_err"]) for m in mods}
            t = TR.transfer(edge, ls, lp, gp, reg)
            ts = TR.transfer(edge, ls, lp, gp, reg, groups=("pool", "H"))
            trf.append((rt, n, t["spread"], t["sd"], t["expected"], ts["spread"]))
            print(f"[transfer] {rt} {n}: pool spread {100 * t['spread']:.2f} % sd {100 * t['sd']:.2f} % pool+H {100 * ts['spread']:.2f} %")
    pool = [m for m in mods if reg.get(m) == "pool"]
    with open(args.out, "w") as f:
        f.write(f"# P14.2b edge pile-up cross-check, run 23701 (9 local files), T0 = {T0[0]:.1f} ns; settings: " +
                "; ".join(f"{n} = {d}" for n, _, _, d in SETTINGS) + "\n")
        f.write("# E_edge: 3x3 sum in MeV of the production W, two-pass fit as edge_pass (s = E_edge/2108.85); sys_spread_pct = (max-min of s over "
                "the settings)/s_w1; sys_shift_none_pct = s(none_sum)/s(w1) - 1; s_err_w1_pct = wave-1 s_err/s\n")
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
        f.write("# transfer: pool spread (1.48 MAD of r), sd and expected, and pool+H spread for each setting (L_scan = 23701 alone; N_rt refit)\n")
        f.write("# rt\tsetting\tpool_spread_pct\tpool_sd_pct\tpool_expected_pct\tpoolH_spread_pct\n")
        for rt, n, sp, sd, ex, sph in trf:
            f.write(f"# {rt}\t{n}\t{100 * sp:.3f}\t{100 * sd:.3f}\t{100 * ex:.3f}\t{100 * sph:.3f}\n")
    print(f"[analyze] wrote {args.out}; pool modules {pool}")


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("collect")
    p.add_argument("evio", nargs="+")
    p.add_argument("--dwell", required=True)
    p.add_argument("--map", default=str(SSN.DERIVED_MAP))
    p.add_argument("--cache", required=True)
    p.add_argument("--max-events", type=int, default=0)
    p.add_argument("-j", "--jobs", type=int, default=1)
    p.set_defaults(fn=cmd_collect)
    q = sub.add_parser("analyze")
    q.add_argument("--cache", required=True)
    q.add_argument("--out", required=True)
    q.set_defaults(fn=cmd_analyze)
    args = ap.parse_args()
    return args.fn(args)


if __name__ == "__main__":
    sys.exit(main() or 0)
