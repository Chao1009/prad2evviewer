#!/usr/bin/env python3
"""examples_mfchi2_gt5.py - 50 random examples of non-seed hits with |dt(CFD fit, calibrated)| < 2 ns and MF chi2/dof > 5: raw waveform, MF fit,
residual, markers.  Sample: 2.24 GeV H2, single-cluster > 100 MeV (cache/2p24GeV_H2_e100).  Re-reads the EVIO for the selected (split, event, channel)
and rebuilds the benchmark's MF reconstruction (same template bank, nsigma 5, min_amp 5, kmax 4, r_model 0.3, mad noise, same pedestal).
Output: figures/examples/mfchi2_gt5_core_2p24GeV_H2.pdf (+ PNG pages in a temporary folder).  Run from analysis/pyscripts after setup_env.sh."""
import re
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
PY = HERE.parent
sys.path.insert(0, str(HERE)); sys.path.insert(0, str(PY))
import crosstalk_study as C                          # noqa: E402

SLUG, CACHE = "2p24GeV_H2", "2p24GeV_H2_e100"
SEED = 20261107
N_EX = 50


def population():
    off = None
    rows = []
    for si, s in enumerate(C.TYPE_SPLITS[SLUG]):
        run, sp = map(int, s.split("/"))
        z = dict(np.load(HERE / "cache" / CACHE / f"{run:06d}_{sp:05d}_legacy.npz"))
        names = z["module_names"]
        if off is None:
            off, _ = C.load_offsets("auto", names)
        ok = (z["algo"] == 2) & (z["seed_algo"] == 2) & np.isfinite(z["t_wa"]) & np.isfinite(z["chi2"])
        dt = (z["t_wa"] - off[z["mod_i"]]) - (z["seed_t_wa"] - off[z["seed_mod_i"]])
        sel = ok & (np.abs(dt) < 2) & (z["chi2"] > 5)
        idx = np.nonzero(sel)[0]
        rows.append(dict(run=run, sp=sp, z=z, idx=idx, dt=dt, names=names, off=off))
    return rows


def pick(rows):
    rng = np.random.default_rng(SEED)
    quota = [N_EX // len(rows) + (1 if i < N_EX % len(rows) else 0) for i in range(len(rows))]
    out = []
    for r, q in zip(rows, quota):
        order = rng.permutation(r["idx"]); evs, mods = set(), set(); n = 0
        for i in order:
            e, m = int(r["z"]["ev"][i]), int(r["z"]["mod_i"][i])
            if e in evs or m in mods:
                continue
            evs.add(e); mods.add(m); out.append((r, int(i))); n += 1
            if n == q:
                break
    return out


def grab(job):
    """EVIO samples of the selected channels: job = (run, split, [(ev, mod_i), ...], module_names) -> {(ev, mod_i): samples}."""
    run, sp, want, mnames = job
    import _common as C_
    from prad2py import dec
    evio, calib = C.evio_of(run, sp)
    p = C_.setup_pipeline(evio_path=evio, hc_calib_file=calib)
    daq = {}
    for i, n in enumerate(mnames):
        m = p.hycal.module_by_name(str(n))
        if m is not None and m.is_hycal():
            daq[i] = (m.daq.crate, m.daq.slot, m.daq.channel)
    need = {}
    for e, m in want:
        need.setdefault(e, []).append(m)
    got = {}
    ch = dec.EvChannel(); ch.set_config(p.cfg)
    ch.open_auto(evio)
    while ch.read() == dec.Status.success and len(got) < len(want):
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            d = ch.decode_event(i, with_ssp=False)
            if not d["ok"]:
                continue
            e = d["event"]; evn = int(e.info.event_number)
            if evn not in need:
                continue
            rocs = {p.crate_map[e.roc(ri).tag]: e.roc(ri) for ri in range(e.nrocs) if e.roc(ri).present and p.crate_map.get(e.roc(ri).tag) is not None}
            for mi in need[evn]:
                a = daq.get(mi); roc = rocs.get(a[0]) if a else None
                if roc is None:
                    continue
                cd = roc.slot(a[1]).channel(a[2])
                if cd.nsamples > 0:
                    got[(evn, mi)] = np.asarray(cd.samples, dtype=np.uint16).copy()
    return got


def main():
    import multiprocessing as mp
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.backends.backend_pdf import PdfPages
    import mf_pulse
    import _common as C_

    rows = population()
    npop = sum(len(r["idx"]) for r in rows)
    sel = pick(rows)
    print(f"[population] {npop:,} hits with |dt|<2 ns and MF chi2/dof > 5 (2.24 GeV H2, single > 100 MeV, fit-only calibrated); drawn {len(sel)} (seed {SEED})", flush=True)
    jobs = {}
    for r, i in sel:
        jobs.setdefault((r["run"], r["sp"]), (r, []))[1].append((int(r["z"]["ev"][i]), int(r["z"]["mod_i"][i])))
    with mp.get_context("spawn").Pool(min(6, len(jobs))) as pool:
        res_g = pool.map(grab, [(k[0], k[1], v[1], v[0]["names"]) for k, v in jobs.items()])
    samples = {}
    for k, g in zip(jobs, res_g):
        samples[k] = g
    # MF reconstruction exactly as the benchmark
    r0 = sel[0][0]
    cmd, _ = C.final_cmd(r0["run"], r0["sp"])
    mf_path = re.search(r"--mf-templates (\S+)", cmd).group(1)
    evio, calib = C.evio_of(r0["run"], r0["sp"])
    p = C_.setup_pipeline(evio_path=evio, hc_calib_file=calib)
    hy = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    wcfg = p.wave_ana.cfg
    col = mf_pulse.resolve_column(mf_path, "auto", prefer="25320")
    bank, binfo = mf_pulse.build_bank([m.name for m in hy], [bool(m.is_pwo4()) for m in hy], str(PY / mf_path) if not Path(mf_path).is_absolute() else mf_path,
                                      dict(smooth_order=int(wcfg.smooth_order), int_tail_ratio=float(wcfg.int_tail_ratio), tail_break_n=int(wcfg.tail_break_n)),
                                      column=col, int_tail_ratio=None)
    rec = mf_pulse.MFReconstructor(bank, nsigma=5.0, min_amp=5.0, kmax=4, r_model=0.3, ped_nsamples=int(wcfg.ped_nsamples))
    ex = []
    for r, i in sel:
        z = r["z"]; ev, mi = int(z["ev"][i]), int(z["mod_i"][i])
        smp = samples[(r["run"], r["sp"])].get((ev, mi))
        if smp is None:
            print("[warn] missing samples", r["run"], r["sp"], ev, mi); continue
        ex.append(dict(r=r, i=i, ev=ev, mi=mi, y=smp))
    Y = np.stack([e["y"] for e in ex]); mid = np.array([e["mi"] for e in ex])
    res = rec.run(Y, mid, int(wcfg.overflow))
    starts = np.searchsorted(res.p_ch, np.arange(len(ex) + 1))
    dchi = []
    ns = Y.shape[1]
    tt = np.arange(0.0, ns * 4.0, 0.5)
    for k, e in enumerate(ex):
        z, i = e["r"]["z"], e["i"]
        sl = slice(starts[k], starts[k + 1])
        e["t0"], e["amp"], e["chi2_re"], e["ped"] = res.p_t0[sl], res.p_amp[sl], (res.p_chi2[sl][0] if starts[k + 1] > starts[k] else np.nan), float(res.ped[k])
        e["K"] = int(starts[k + 1] - starts[k]); e["tmf"] = res.p_time[sl]
        fit = np.full_like(tt, e["ped"])
        for a_, t0_ in zip(e["amp"], e["t0"]):
            fit += a_ * bank.shape(np.array([e["mi"]])[:, None], (tt - t0_)[None, :])[0]
        e["fit"] = fit
        fs = np.full(ns, e["ped"], dtype=float)
        tsm = 4.0 * np.arange(ns)
        for a_, t0_ in zip(e["amp"], e["t0"]):
            fs += a_ * bank.shape(np.array([e["mi"]])[:, None], (tsm - t0_)[None, :])[0]
        e["res"] = e["y"].astype(float) - fs
        e["win"] = [(int(w0), int(w0) + mf_pulse.LW) for w0 in np.unique(rec._window_start(e["t0"]))] if e["K"] else []
        e["chi2_cache"] = float(z["chi2"][i]); dchi.append(abs(e["chi2_re"] - e["chi2_cache"]))
    dchi = np.array(dchi)
    print(f"[chi2] recomputed vs cached chi2/dof over {len(ex)} examples: max |diff| = {np.nanmax(dchi):.3g}, median {np.nanmedian(dchi):.3g}; "
          f"K != 1: {sum(e['K'] != 1 for e in ex)}", flush=True)
    # figure
    out = HERE / "figures" / "examples"; out.mkdir(parents=True, exist_ok=True)
    pdf = out / "mfchi2_gt5_core_2p24GeV_H2.pdf"
    page_png = Path("/private/tmp/claude-505/-Users-L00338853-work-PRad-prad2evviewer-analysis-pyscripts/db0289d2-4ae1-47bc-9d5f-4788ee635668/scratchpad/ex_pages")
    page_png.mkdir(exist_ok=True)
    with PdfPages(pdf) as pp:
        fig = plt.figure(figsize=(8.27, 11.69)); fig.add_subplot(111).axis("off")
        cap = (f"MF chi2/dof > 5 while the CFD dt is within +-2 ns: 50 examples, 2.24 GeV H2\n\n"
               f"Selection: run 25228, all 11 splits (50..550), single-cluster events with exactly one cluster > 100 MeV (benchmark --signal-min 100,\n"
               f"legacy path); non-seed cluster members with WaveAnalyzer log-normal-fit CFD time (time_algo 2) at module and seed;\n"
               f"dt = T_module - T_seed after time_calib.json offsets, |dt| < 2 ns; MF chi2/dof > 5.\n"
               f"Population: {npop:,} hits. Drawn at random (numpy default_rng seed {SEED}), 4-5 per split, at most one per event and per module.\n\n"
               f"Fit: matched filter (mf_pulse), per-module spline template bank {mf_path}, nsigma 5, min amplitude 5 ADC, kmax 4,\n"
               f"r_model 0.3, noise model mad (sigma^2 = sigma_ped^2 + (0.03 model)^2), pedestal = clipped mean of the first {int(wcfg.ped_nsamples)} samples.\n"
               f"Fit window: 32 samples (128 ns) starting 3 samples before the onset ceil(t0/4); amplitude A and onset t0 free;\n"
               f"dof = valid samples in the window - 2 (single pulse). Blue: raw samples (ADC, 4 ns/sample); red: MF fit; grey band: fit window;\n"
               f"lower panel: data - fit. Vertical lines: t_MF = onset + t50 (red), t_CFD module (green), t_CFD seed (black dashed);\n"
               f"times are raw (uncalibrated) ns on the readout axis; dt in the text is calibrated.\n"
               f"chi2/dof recomputed here vs cached: max |diff| = {np.nanmax(dchi):.2g} over {len(ex)} examples.")
        fig.text(0.06, 0.95, cap, va="top", fontsize=6.8)
        pp.savefig(fig); fig.savefig(page_png / "p0.png", dpi=105); plt.close(fig)
        for pg in range(0, len(ex), 6):
            fig = plt.figure(figsize=(8.27, 11.69))
            outer = fig.add_gridspec(3, 2, hspace=0.28, wspace=0.22, left=0.07, right=0.98, top=0.97, bottom=0.04)
            for j, e in enumerate(ex[pg:pg + 6]):
                inner = outer[j // 2, j % 2].subgridspec(2, 1, height_ratios=[3, 1], hspace=0.05)
                ax = fig.add_subplot(inner[0]); ar = fig.add_subplot(inner[1], sharex=ax)
                z, i = e["r"]["z"], e["i"]
                ax.plot(4.0 * np.arange(ns), e["y"], ".", ms=3, color="C0"); ax.plot(tt, e["fit"], "-", lw=1, color="C3")
                for w0, w1 in e["win"]:
                    ax.axvspan(4.0 * w0, 4.0 * w1, color="0.85", zorder=0); ar.axvspan(4.0 * w0, 4.0 * w1, color="0.85", zorder=0)
                for t_ in e["tmf"]:
                    ax.axvline(t_, color="C3", lw=0.8)
                ax.axvline(float(z["t_wa"][i]), color="C2", lw=0.8); ax.axvline(float(z["seed_t_wa"][i]), color="k", ls="--", lw=0.8)
                ar.plot(4.0 * np.arange(ns), e["res"], ".-", ms=2, lw=0.5, color="0.3"); ar.axhline(0, color="C3", lw=0.5)
                dtc = float(e["r"]["dt"][i])
                lab = (f"{e['r']['run']}/{e['r']['sp']} ev {e['ev']} {e['r']['names'][e['mi']]}\n"
                       f"A={float(z['amp'][i]):.1f} chi2/dof={e['chi2_cache']:.1f} K={e['K']}\n"
                       f"dt_CFD={dtc:+.2f} tCFD-tMF={float(z['t_wa'][i]) - float(z['t'][i]):+.2f} ns\n"
                       f"seed A={float(z['seed_amp'][i]):.0f} algo={int(z['algo'][i])} lnchi2={float(z['lnchi2'][i]):.1f}")
                ax.text(0.99, 0.97, lab, transform=ax.transAxes, ha="right", va="top", fontsize=6, bbox=dict(fc="w", ec="0.7", alpha=0.85))
                ax.tick_params(labelbottom=False, labelsize=6); ar.tick_params(labelsize=6)
                ar.set_xlabel("ns", fontsize=6); ax.set_ylabel("ADC", fontsize=6); ar.set_ylabel("data-fit", fontsize=6)
            pp.savefig(fig); fig.savefig(page_png / f"p{pg // 6 + 1}.png", dpi=105); plt.close(fig)
    print(f"[done] {pdf}  ({len(ex)} examples, {len(ex) // 6 + (len(ex) % 6 > 0)} pages + caption)")
    np.savez_compressed(page_png / "examples_meta.npz", chi2_re=[e["chi2_re"] for e in ex], chi2_cache=[e["chi2_cache"] for e in ex])


if __name__ == "__main__":
    main()
