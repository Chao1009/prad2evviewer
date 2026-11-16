#!/usr/bin/env python3
"""p6_table.py: P6 gain-arm table (D12).  Reads output/energy_resolution/variants/P6/<variant>/<RT>/{delta.tsv,summary.tsv,gains.tsv,regions_splits.tsv}
against the phase-2 control (or the gains-only control g0 for the gains-only variants) and writes improvement/P6_gain_arm.tsv; prints a markdown table.
Usage: p6_table.py [variant ...]  (default: every variant dir found).  Fitted factors / events per module come from gains.tsv (free = prior_frac < 0.5)."""
import csv, math, sys
from pathlib import Path
import numpy as np

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "output/energy_resolution"
V6 = OUT / "variants/P6"
CTRL = OUT / "control_phase2_20261113"
ER = Path(__file__).resolve().parents[1]
GATE = {"0p73GeV_H2": "g20", "2p24GeV_H2": "g16", "3p49GeV_H2": "g16", "3p49GeV_H2_old": "g16", "2p24GeV_Ta_X17": "g8"}
GONLY = ("g0_ctrl_gonly", "g1_seed_gonly")


def rd(p):
    with open(p) as f:
        return list(csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t"))


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def one(rows, **kw):
    return next((r for r in rows if all(r[k] == v for k, v in kw.items())), None)


def gains_stats(d):
    rows = rd(d / "gains.tsv")
    rows = [r for r in rows if r["moller"] == "0" and r["nonlinear"] == "0"]
    runs = sorted({r["run"] for r in rows})
    free = [r for r in rows if fl(r["prior_frac"]) < 0.5]
    nm = np.array([fl(r["n_member"]) for r in free])
    return dict(n_fit=int(round(np.mean([sum(r["run"] == k for r in rows) for k in runs]))),
                n_free=int(round(np.mean([sum(r["run"] == k for r in free) for k in runs]))),
                nmem_med=float(np.median(nm)) if nm.size else float("nan"), nmem_min=float(nm.min()) if nm.size else float("nan"))


def row(v, rt, base=None):
    d = V6 / v / rt
    if not (d / "delta.tsv").exists():
        return []
    ctl = base or CTRL / rt
    out = []
    ovr = open(d / "regions_splits.tsv").readline().strip().removeprefix("# overrides:").strip()
    gs = gains_stats(d)
    for path in ("legacy", GATE[rt]):
        dl = one(rd(d / "delta.tsv"), path=path, kind="ep", region="pool")
        sv, sc = one(rd(d / "summary.tsv"), path=path), one(rd(ctl / "summary.tsv"), path=path)
        if not (dl and sv and sc):
            continue
        g = lambda q: fl(dl[q])
        out.append(dict(variant=v, run_type=rt, path=path, overrides=ovr,
                        ctrl_in=fl(sc["res_cb"]), ctrl_oos=fl(sc["res_cb_oos"]), var_in=fl(sv["res_cb"]), var_oos=fl(sv["res_cb_oos"]),
                        d_in=g("d_res_cb"), d_in_se=g("d_res_cb_se"), d_oos=g("d_res_cb_oos"), d_oos_se=g("d_res_cb_oos_se"),
                        d_tail_hi3_oos=g("d_tail_hi3_oos"), d_tail_hi3_oos_se=g("d_tail_hi3_oos_se"), d_leak_oos=g("d_leak_cut_oos"),
                        d_leak_oos_se=g("d_leak_cut_oos_se"), d_eff_cut=g("d_eff_cut"), d_eff_cut_se=g("d_eff_cut_se"),
                        d_tail_lo_f=g("d_tail_lo_f"), d_tail_lo_f_se=g("d_tail_lo_f_se"), ctrl_gap=fl(sc["res_cb_oos"]) - fl(sc["res_cb"]),
                        var_gap=fl(sv["res_cb_oos"]) - fl(sv["res_cb"]), d_gap=g("gap"), d_gap_se=g("gap_se"), yield_loss=g("yield_loss"),
                        n_pairs=int(g("n")), **gs, verdict=judge(dl)))
    return out


def judge(dl):
    """D3 / U2 verdict on the out-of-sample res_cb and the guards: gain / null / worse, then the U2 gap rule and the guards."""
    d, se = fl(dl["d_res_cb_oos"]), fl(dl["d_res_cb_oos_se"])
    core = "better" if d <= -2 * se and d <= -0.01 else "worse" if d >= 2 * se and d >= 0.01 else "null"
    gap = fl(dl["gap"]) > max(2 * fl(dl["gap_se"]), 0.02)
    guards = []
    if fl(dl["d_eff_cut"]) < -2 * fl(dl["d_eff_cut_se"]) and fl(dl["d_eff_cut"]) < -0.1:
        guards.append("eff_cut down")
    if fl(dl["d_tail_lo_f"]) > 2 * fl(dl["d_tail_lo_f_se"]) and fl(dl["d_tail_lo_f"]) > 0.05:
        guards.append("tail_lo_f up")
    if fl(dl["yield_loss"]) > 1.0:
        guards.append("yield loss > 1 %")
    return core + ("; U2 gap grows" if gap else "") + ("; guard: " + ", ".join(guards) if guards else "")


def main():
    vs = sys.argv[1:] or sorted(p.name for p in V6.iterdir() if p.is_dir() and p.name != "logs")
    rows = []
    for v in vs:
        for rt in GATE:
            rows += row(v, rt, V6 / "g0_ctrl_gonly" / rt if v == "g1_seed_gonly" else None)
    if not rows:
        sys.exit("no rows")
    cols = list(rows[0])
    with open(ER / "P6_gain_arm.tsv", "w") as f:
        f.write("# P6 gain-solve arms; delta = variant - control (phase-2 control; gains-only variant g1 against g0); path legacy and the control gate; res_cb in %, deltas in pp, SE paired over splits\n")
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[c]:.5g}" if isinstance(r[c], float) else str(r[c]) for c in cols) + "\n")
    p = lambda a, b, n=3: f"{a:+.{n}f} ± {b:.{n}f}"
    print("| variant | run type | path | ctrl in / oos | var in / oos | Δ in | Δ oos | Δ tail_hi3 oos | Δ leak_cut oos | Δ eff_cut | Δ tail_lo_f | Δ gap | yield loss | free / fit | events per free module med / min | verdict |")
    print("|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for r in rows:
        print(f"| {r['variant']} | {r['run_type']} | {r['path']} | {r['ctrl_in']:.3f} / {r['ctrl_oos']:.3f} | {r['var_in']:.3f} / {r['var_oos']:.3f} | {p(r['d_in'], r['d_in_se'])} | "
              f"{p(r['d_oos'], r['d_oos_se'])} | {p(r['d_tail_hi3_oos'], r['d_tail_hi3_oos_se'])} | {p(r['d_leak_oos'], r['d_leak_oos_se'], 2)} | {p(r['d_eff_cut'], r['d_eff_cut_se'], 2)} | "
              f"{p(r['d_tail_lo_f'], r['d_tail_lo_f_se'], 2)} | {p(r['d_gap'], r['d_gap_se'])} | {r['yield_loss']:+.2f} | {r['n_free']} / {r['n_fit']} | {r['nmem_med']:.0f} / {r['nmem_min']:.0f} | {r['verdict']} |")


main()
