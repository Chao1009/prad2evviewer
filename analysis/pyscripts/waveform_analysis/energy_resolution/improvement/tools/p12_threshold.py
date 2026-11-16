#!/usr/bin/env python3
"""P12 readout threshold.  Usage (source setup_env.sh first):
  p12_threshold.py raw  <run> <evio split> <out.npz> [max physics events]   decode: PbWO4 channels read out per event, height above pedestal
  p12_threshold.py edge <scan npz> [<scan npz> ...] <out.npz>              per module: lower edge of h_a (ADC) and h_e (MeV), MeV per ADC
Reads the BM helpers by import; writes nothing else."""
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[3]))
import energy_resolution_benchmark as bm
from prad2py import dec

NPED = 8          # leading samples whose median is the pedestal (the DAQ pedestal itself is not in the decoded data)


def raw(run, path, out, nmax=40000):
    cfg = bm.run_config(run)
    core = bm.build(run, str(bm.DB_DIR / cfg["calib"]))
    hy, cm = core.hycal, {int(r.tag): int(r.crate) for r in core.daq_cfg.roc_tags}
    mods = [hy.module(i) for i in range(hy.module_count())]
    ch = dec.EvChannel()
    ch.set_config(core.daq_cfg)
    ch.open_auto(path)
    gate = (lambda b: b & bm.PHYSICS_MASK_ANY != 0 and not b & bm.LMS_BIT) if bm.RUN_TYPE_OF[run] in bm.PHYSICS_ANY else (lambda b: b == bm.PHYSICS_BITS)
    mi_l, h_l, ev_l, ns_l, nsat_l, k_l, n = [], [], [], [], [], [], 0
    while ch.read() == dec.Status.success and n < nmax:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"] or not gate(int(d["event"].info.trigger_bits)):
                continue
            fe, k = d["event"], 0
            for ri in range(fe.nrocs):
                roc = fe.roc(ri)
                if not roc.present or roc.tag not in cm:
                    continue
                for sl in roc.present_slots():
                    s = roc.slot(sl)
                    for c in s.present_channels():
                        m = hy.module_by_daq(cm[roc.tag], sl, c)
                        if m is None or not m.is_pwo4():
                            continue
                        cd, k = s.channel(c), k + 1
                        y = np.asarray(cd.samples, float)
                        mi_l.append(m.index); ns_l.append(cd.nsamples)
                        h_l.append(y.max() - np.median(y[:NPED]) if y.size >= NPED else np.nan)
                        # consecutive samples above TET (firmware NSAT) in the decoded window
                        nsat_l.append(int(((y - np.median(y[:NPED])) > 10.0).sum()) if y.size else 0)
            k_l.append(k); n += 1
    np.savez_compressed(out, run=run, n_events=n, n_chan=np.array(k_l), mod=np.array(mi_l), nsamples=np.array(ns_l), height=np.array(h_l),
                        n_above_tet=np.array(nsat_l), names=np.array([m.name for m in mods]), pwo4=np.array([m.is_pwo4() for m in mods]))
    print(f"[done] {out}: {n} events, {np.mean(k_l):.2f} PbWO4 channels per event")


def edge(files, out):
    A, E, M = [], [], []
    for f in files:
        z = np.load(f, allow_pickle=False)
        pw = z["mod_pwo4"][z["h_mod"]] & np.isfinite(z["h_a"]) & (z["h_a"] > 0) & ~z["h_extra"]
        A.append(z["h_a"][pw]); E.append(z["h_e"][pw]); M.append(z["h_mod"][pw].astype(int)); names = z["mod_name"]
    A, E, M = map(np.concatenate, (A, E, M))
    n, q = len(names), (0.5, 2, 5, 50)
    cnt = np.bincount(M, minlength=n)
    res = {"count": cnt, "k_mev_per_adc": np.full(n, np.nan)}
    for nm, v in (("a", A), ("e", E)):
        for p in q:
            res[f"{nm}_p{p:g}"] = np.full(n, np.nan)
    for m in np.nonzero(cnt >= 50)[0]:
        s = M == m
        res["k_mev_per_adc"][m] = np.median(E[s] / A[s])
        for nm, v in (("a", A), ("e", E)):
            for p, x in zip(q, np.percentile(v[s], q)):
                res[f"{nm}_p{p:g}"][m] = x
    np.savez_compressed(out, names=names, **res)
    print(f"[done] {out}: {int((cnt >= 50).sum())} modules with >= 50 hits")


if __name__ == "__main__":
    a = sys.argv[1:]
    if a[:1] == ["raw"]:
        raw(int(a[1]), a[2], a[3], int(a[4]) if len(a) > 4 else 40000)
    elif a[:1] == ["edge"]:
        edge(a[1:-1], a[-1])
    else:
        sys.exit(__doc__)
