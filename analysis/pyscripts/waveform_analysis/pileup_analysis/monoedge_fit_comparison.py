#!/usr/bin/env python3
"""monoedge_fit_comparison.py -- 50 waveforms where the monotonic-edge rule changed the fit, each refitted from the raw waveform with the rule OFF and ON.

    source setup_env.sh
    python3 waveform_analysis/pileup_analysis/monoedge_fit_comparison.py [--sample-dir <dir>] [--evio-dir <dir>] [--out <dir>]

Candidates: channels matched row-for-row (run, ev, midx) in the old (rule off) and new (rule on) bounded scans with k_old > k_new, not saturated
(nclip < SAT_NCLIP), chi2/dof <= CHI2_MAX in both.  Removed pulse = old pulse without a new pulse within TOL ns (largest r if several);
dt = t_removed - t_main(old, largest amplitude), r = A_removed / A_main.  Four categories (fixed seed SEED), quotas below.
Each selected waveform is decoded again and refitted twice (pr.MONO_EDGE False / True); a refit must reproduce the stored k and chi2/dof, else the
next candidate of the category replaces it.  -> monoedge_fit_comparison_20261116.pdf / .tsv
"""
from __future__ import annotations

import sys
import time
import argparse
import textwrap
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.dont_write_bytecode = True
sys.path.insert(0, str(HERE))
import numpy as np                                      # noqa: E402
import matplotlib                                       # noqa: E402
matplotlib.use("Agg")
import matplotlib.pyplot as plt                         # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages    # noqa: E402
import pileup_monoedge as pmo                           # noqa: E402  (decode helpers, protected stretch; sets env MONO_EDGE=1 before pileup_rates import, toggled below)
from prad2py import dec                                 # noqa: E402

pm, pr, ed = pmo.pm, pmo.pr, pmo.ed
SEED = 20261116
TOL = 2.0
SPARE = 12
H2 = ("0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2")
TA = ("2p24GeV_Ta_X17",)
# (key, title, run types, quota per run type, r range, dt range)
CATS = [
    ("H2_small", "A. H2, small removed pulse (r < 0.1, dt 0-60 ns): template-mismatch case", H2, (7, 7, 6), (0, 0.1), (0, 60)),
    ("H2_large", "B. H2, large removed pulse (r >= 0.2)", H2, (3, 3, 4), (0.2, 1e9), (-1e9, 1e9)),
    ("Ta_small", "C. Ta_X17, small removed pulse (r < 0.1, dt 0-60 ns)", TA, (10,), (0, 0.1), (0, 60)),
    ("Ta_large", "D. Ta_X17, large removed pulse (r >= 0.2): suspected real pile-up", TA, (10,), (0.2, 1e9), (-1e9, 1e9)),
]
SHADE_RISE, SHADE_FALL = pmo.SHADE_RISE, pmo.SHADE_FALL
KEYS = ("ch_ev", "ch_midx", "ch_k", "ch_a", "ch_t", "ch_chi2", "ch_nclip")


def collect(sample_dir, evio_dir):
    """All candidates: list of dict(run, rt, ev, midx, file, cat-features)."""
    out = []
    for fn in sorted((Path(sample_dir) / "new").glob("run_*.npz")):
        zn_f, zo_f = np.load(fn), np.load(Path(sample_dir) / "old" / fn.name)
        run, rt = int(zn_f["run"]), str(zn_f["run_type"])
        zn, zo = {k: zn_f[k] for k in KEYS}, {k: zo_f[k] for k in KEYS}
        evn = zn_f["ev_number"]
        assert np.array_equal(zn_f["splits"], zo_f["splits"])
        files = [str(Path(evio_dir) / f"prad_{run:06d}.evio.{int(s):05d}") for s in zn_f["splits"]]
        ko_ = zo["ch_ev"].astype(np.int64) * 4096 + zo["ch_midx"]
        kn_ = zn["ch_ev"].astype(np.int64) * 4096 + zn["ch_midx"]
        _, io, inn = np.intersect1d(ko_, kn_, return_indices=True)
        ko, kn = zo["ch_k"][io], zn["ch_k"][inn]
        ok = (ko > kn) & (zo["ch_nclip"][io] < pr.SAT_NCLIP) & (zn["ch_nclip"][inn] < pr.SAT_NCLIP) & (zo["ch_chi2"][io] <= pr.CHI2_MAX) & (zn["ch_chi2"][inn] <= pr.CHI2_MAX)
        sel = np.nonzero(ok)[0]
        ctx = ed.Ctx(run)
        firsts = np.array([ed.first_event_number(ctx, f) for f in files], dtype=np.int64)
        for s in sel:
            a_o, t_o, a_n, t_n = zo["ch_a"][io[s]], zo["ch_t"][io[s]], zn["ch_a"][inn[s]], zn["ch_t"][inn[s]]
            vo, vn = np.arange(6) < ko[s], np.arange(6) < kn[s]
            im = int(np.argmax(np.where(vo, a_o, -1)))
            surv = (np.abs(t_o[:, None] - np.where(vn, t_n, np.inf)[None, :]) <= TOL).any(1)
            rem = vo & ~surv
            rem[im] = False
            if not rem.any():
                continue
            rr = np.where(rem, a_o / a_o[im], -1)
            j = int(np.argmax(rr))
            e = int(evn[zo["ch_ev"][io[s]]])
            fi = int(np.clip(np.searchsorted(firsts, e, side="right") - 1, 0, len(files) - 1))
            out.append(dict(run=run, rt=rt, ev=e, midx=int(zo["ch_midx"][io[s]]), file=files[fi], split=fi, k_off=int(ko[s]), k_on=int(kn[s]),
                            chi2_off=float(zo["ch_chi2"][io[s]]), chi2_on=float(zn["ch_chi2"][inn[s]]), dt=float(t_o[j] - t_o[im]), r=float(rr[j])))
        print(f"[cand] run {run} {rt}: {len(sel)} k decreased, {sum(c['run'] == run for c in out)} with removed pulse", flush=True)
    return out


def pick(cands, rng):
    """{cat key: ordered candidate list (quota + SPARE per run type), quota map}."""
    res = {}
    for key, title, rts, quota, (rlo, rhi), (dlo, dhi) in CATS:
        lst = []
        for rt, q in zip(rts, quota):
            pool = [c for c in cands if c["rt"] == rt and rlo <= c["r"] < rhi and dlo <= c["dt"] <= dhi]
            print(f"[pool] {key} {rt}: {len(pool)}", flush=True)
            order = rng.permutation(len(pool))
            lst += [dict(pool[i], cat=key, quota_rt=rt, q=q, rank=n) for n, i in enumerate(order[:q + SPARE])]
        res[key] = lst
    return res


def fit_both(sample_dir, chosen):
    """Decode each needed event once; fit with the rule off and on.  Returns {(file, ev): {"off": r, "on": r}}."""
    refs = {}
    for fn in sorted((Path(sample_dir) / "new").glob("run_*.npz")):
        z = np.load(fn)
        run = int(z["run"])
        for i, s in enumerate(z["splits"]):
            refs[(run, i)] = (z["mod_b_ref"][i], z["mod_sig_ref"][i])
    need = {}
    for lst in chosen.values():
        for c in lst:
            need.setdefault((c["run"], c["file"], c["split"]), set()).add(c["ev"])
    ctxs, store = {}, {}
    for (run, f, si), evs in sorted(need.items()):
        ctx = ctxs.setdefault(run, ed.Ctx(run))
        b_ref, sig_ref = refs[(run, si)]
        left = set(evs)
        ch = ctx.open(f)
        t0 = time.time()
        while left and ch.read() == dec.Status.success:
            if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
                continue
            for i in range(ch.get_n_events()):
                ch.select_event(i)
                en = int(ch.info().event_number)
                if en in left:
                    d = ch.decode_event(i, with_ssp=False)
                    if d["ok"]:
                        fe = d["event"]
                        out = {}
                        for name, flag in (("off", False), ("on", True)):
                            pr.MONO_EDGE = flag
                            out[name] = ed.fit_event(ctx, fe, b_ref, sig_ref)
                        store[(f, en)] = out
                    left.discard(en)
        ch.close()
        print(f"[fit] run {run} {Path(f).name}: {len(evs)} events, {time.time() - t0:.0f} s", flush=True)
    pr.MONO_EDGE = False
    return store, ctxs, refs


def make(ctx, c, r, tag):
    i = r["row"].get(c["midx"])
    if i is None:
        return None
    N = int(r["k"][i])
    return dict(run=c["run"], rt=c["rt"], ev=c["ev"], midx=c["midx"], name=ctx.names[c["midx"]], pwo4=bool(ctx.pwo4[c["midx"]]), N=N, Y=r["Y"][i].astype(float),
                ped=float(r["ped"][i]), a=r["a"][i][:N].copy(), t=r["t"][i][:N].copy(), comp=r["comp"][i][:N].copy(), chi2=float(r["chi2"][i]), nclip=int(r["nclip"][i]))


def assemble(chosen, store, ctxs, refs):
    """Verified items per category (quota order, spares fill failed verifications); each item has the OFF and ON examples."""
    final, nfail = {}, 0
    for key, title, rts, quota, _, _ in CATS:
        items = []
        for rt, q in zip(rts, quota):
            got = 0
            for c in [c for c in chosen[key] if c["quota_rt"] == rt]:
                if got >= q:
                    break
                r = store.get((c["file"], c["ev"]))
                if r is None:
                    nfail += 1
                    continue
                ctx = ctxs[c["run"]]
                eo, en = make(ctx, c, r["off"], "off"), make(ctx, c, r["on"], "on")
                if eo is None or en is None or eo["N"] != c["k_off"] or en["N"] != c["k_on"] or not np.isclose(eo["chi2"], c["chi2_off"], rtol=2e-3, atol=1e-3) \
                        or not np.isclose(en["chi2"], c["chi2_on"], rtol=2e-3, atol=1e-3):
                    nfail += 1
                    print(f"[verify FAIL] {c['run']} ev {c['ev']} m {c['midx']}: stored k {c['k_off']}/{c['k_on']} chi2 {c['chi2_off']:.3f}/{c['chi2_on']:.3f}", flush=True)
                    continue
                refb, refs_ = refs[(c["run"], c["split"])]
                en_x = pmo.extras(ctx, en, refb, refs_)
                rm = [j for j in range(eo["N"]) if not (np.abs(eo["t"][j] - en["t"]) <= TOL).any()]
                items.append(dict(c, off=eo, on=en, prot=en_x["prot"], removed=eo["t"][rm]))
                got += 1
            print(f"[assemble] {key} {rt}: {got} of {q}", flush=True)
        final[key] = items
    return final, nfail


def panel(ax, axr, e, title, prot=None, removed=None):
    pm.draw_panel(ax, e)
    pulses = "\n".join(f"{j + 1}: A {a:.0f} @ t0 {t:.0f}" for j, (a, t) in enumerate(zip(e["a"], e["t"])))
    ax.set_title(f"{title}   N = {e['N']}   chi2/dof = {e['chi2']:.2f}", fontsize=9, loc="left", pad=3)
    ax.text(0.99, 0.97, pulses, transform=ax.transAxes, fontsize=7, va="top", ha="right", family="monospace", bbox=dict(fc="white", alpha=0.75, lw=0))
    if prot is not None:
        lo1, hi1, lo2, hi2 = prot
        if hi1 >= lo1:
            ax.axvspan(lo1, hi1, color=SHADE_RISE, alpha=0.22, lw=0)
        if hi2 >= lo2:
            ax.axvspan(lo2, hi2, color=SHADE_FALL, alpha=0.18, lw=0)
        for t in removed:
            ax.axvline(t, color="m", ls=":", lw=1.8)
    ov = e["Y"] >= pr.OVERFLOW
    res = e["Y"] - (e["ped"] + e["comp"].sum(0))
    axr.plot(ed.TT[~ov], res[~ov], "k.-", lw=0.6, ms=3)
    axr.axhline(0, color="r", lw=0.8)
    axr.set_xlim(0, ed.DT * ed.NS)
    axr.tick_params(labelsize=7)
    ax.tick_params(labelsize=7, labelbottom=False)
    axr.set_ylabel("data - fit", fontsize=7)
    axr.set_xlabel("t (ns)", fontsize=8)


def title_page(pdf, final, nfail, elapsed):
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.03, 0.95, "Monotonic-edge rule: 50 waveform fits with the rule OFF and ON", fontsize=20, fontweight="bold")
    txt = (
        "The rule (pileup_rates.py MONO_EDGE, PILEUP_MONO_EDGE=1).  Dominant peak P = largest sample.  If the 3 differences up to P do not fall by more than 2 sqrt(2) sigma_ref\n"
        "  (rising-monotonic), no additional pulse may start in [t_main, P]; if the 10 differences after P do not rise by more than 2 sqrt(2) sigma_ref (falling-monotonic), none may\n"
        "  start in [P, P + 40 ns].  The main pulse is exempt, pulses outside the stretch stay allowed.  Rule off = default fit; nothing else differs.\n\n"
        "Panels (one waveform per row): LEFT = rule OFF refit, RIGHT = rule ON refit, both from the raw waveform (verified to reproduce the stored scan k and chi2/dof).\n"
        "  Black dots: raw samples (red x: overflow).  Grey line: pedestal.  Dashed colours: pulse components (blue, orange, green, purple, brown, pink; time order).  Red: total fit.\n"
        "  Small panel below: data - fit (overflow samples omitted).  Text box: (A ADC, t0 ns) of each pulse.  Title: N pulses and chi2/dof.\n"
        "  RIGHT panel only: orange band = protected stretch of a rising-monotonic edge, blue band = falling-monotonic one (from the largest ON pulse);\n"
        "  magenta dotted lines = onsets of the pulses of the OFF fit that the rule removed.\n"
        "Row title: run, event, module, run type, removed pulse dt = t_removed - t_main (ns) and r = A_removed / A_main, from the OFF fit.\n\n"
        "Selection (seed %d): channels of the bounded old-vs-new sample (13 runs x 2 splits) with k_old > k_new, not saturated (nclip < %d), chi2/dof <= %g in both fits;\n"
        "  random within each category; the removed pulse is the old-fit pulse with no new-fit pulse within %.0f ns (largest r if several).\n" % (SEED, pr.SAT_NCLIP, pr.CHI2_MAX, TOL)
        + "Categories (n shown):\n" + "\n".join(f"  {t}   n = {len(final[k])}" for k, t, *_ in CATS)
        + f"\nCandidates that failed the k / chi2 reproduction check and were replaced: {nfail}.   Run time {elapsed:.0f} s.\n")
    fig.text(0.03, 0.9, txt, fontsize=9.6, va="top", family="monospace", linespacing=1.4)
    pdf.savefig(fig)
    plt.close(fig)


def pages(pdf, final):
    rows = []
    PER = 3
    for key, title, *_ in CATS:
        items = final[key]
        for p0 in range(0, len(items), PER):
            chunk = items[p0:p0 + PER]
            fig = plt.figure(figsize=(16, 12.5))
            outer = fig.add_gridspec(PER, 1, hspace=0.32, left=0.04, right=0.99, top=0.86, bottom=0.04)
            fig.text(0.04, 0.97, title, fontsize=15, fontweight="bold")
            fig.text(0.04, 0.948, f"waveforms {p0 + 1}-{p0 + len(chunk)} of {len(items)}.   LEFT rule OFF, RIGHT rule ON; shaded = protected stretch, magenta dotted = pulse removed by the rule.", fontsize=9)
            for i, it in enumerate(chunk):
                gs = outer[i].subgridspec(2, 2, height_ratios=[3, 1], hspace=0.08, wspace=0.12)
                a_off, r_off = fig.add_subplot(gs[0, 0]), fig.add_subplot(gs[1, 0])
                a_on, r_on = fig.add_subplot(gs[0, 1]), fig.add_subplot(gs[1, 1])
                panel(a_off, r_off, it["off"], "rule OFF")
                panel(a_on, r_on, it["on"], "rule ON", it["prot"], it["removed"])
                ymax = max(it["off"]["Y"].max(), it["on"]["Y"].max())
                lim = (min(a_off.get_ylim()[0], a_on.get_ylim()[0]), max(a_off.get_ylim()[1], a_on.get_ylim()[1]))
                a_off.set_ylim(lim)
                a_on.set_ylim(lim)
                rl = (min(r_off.get_ylim()[0], r_on.get_ylim()[0]), max(r_off.get_ylim()[1], r_on.get_ylim()[1]))
                r_off.set_ylim(rl)
                r_on.set_ylim(rl)
                a_off.text(0.0, 1.17, f"run {it['run']}  ev {it['ev']}  {it['off']['name']} ({'PbWO4' if it['off']['pwo4'] else 'LG'})  {it['rt']}   removed pulse: dt = {it['dt']:.1f} ns, r = {it['r']:.3f}",
                           transform=a_off.transAxes, fontsize=10.5, fontweight="bold")
                rows.append((it["run"], it["ev"], it["off"]["name"], it["rt"], key, it["k_off"], it["k_on"], it["off"]["chi2"], it["on"]["chi2"], it["dt"], it["r"]))
            pdf.savefig(fig)
            plt.close(fig)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--sample-dir", default=str(pmo.PYDIR / "output" / "pileup_rates" / "monoedge_sample"))
    ap.add_argument("--evio-dir", default=pmo.EVIO_DIR)
    ap.add_argument("--out", default=str(HERE))
    a = ap.parse_args()
    t0 = time.time()
    rng = np.random.default_rng(SEED)
    cands = collect(a.sample_dir, a.evio_dir)
    chosen = pick(cands, rng)
    store, ctxs, refs = fit_both(a.sample_dir, chosen)
    final, nfail = assemble(chosen, store, ctxs, refs)
    stem = Path(a.out) / "monoedge_fit_comparison_20261116"
    with PdfPages(str(stem) + ".pdf") as pdf:
        title_page(pdf, final, nfail, time.time() - t0)
        rows = pages(pdf, final)
    with open(str(stem) + ".tsv", "w") as f:
        f.write("run\tevent\tmodule\trun_type\tcategory\tk_off\tk_on\tchi2_off\tchi2_on\tremoved_dt_ns\tremoved_r\n")
        for r in rows:
            f.write("\t".join(f"{x:.4g}" if isinstance(x, float) else str(x) for x in r) + "\n")
    print(f"[done] {len(rows)} comparisons, {time.time() - t0:.0f} s", flush=True)


if __name__ == "__main__":
    main()
