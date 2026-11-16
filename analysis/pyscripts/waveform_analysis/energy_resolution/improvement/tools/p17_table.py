"""P17 table. Run from analysis/pyscripts: python3 waveform_analysis/energy_resolution/improvement/tools/p17_table.py [<run type> ...]
Pairs, per (run, split), the pool ep rows of every asymmetric-window path of output/energy_resolution/variants/P17/<RT> (gate tag =
late + early/100, tags in ASYM below) with the control gate path of control_phase2_20261113/<RT> (D3: variant - control, mean +- SE,
real if |d| >= 2 SE).  Also the symmetric Ta gates of variants/P17/sym_extra.  Writes improvement/P17_asym_window.tsv (all run types
given) and prints the tables."""
import csv, math, sys
from pathlib import Path
import numpy as np

O = Path("output/energy_resolution")
CTL = O / "control_phase2_20261113"
IMP = Path("waveform_analysis/energy_resolution/improvement")
GATE = {"0p73GeV_H2": "g20", "2p24GeV_H2": "g16", "3p49GeV_H2": "g16", "3p49GeV_H2_old": "g16", "2p24GeV_Ta_X17": "g8"}
ASYM = eval((O / "variants/P17/tags.txt").read_text())          # tag -> (early, late)
Q = ("res_cb", "res_cb_oos", "tail_hi3", "tail_hi3_oos", "tail_lo_f", "leak_cut", "eff_cut")
rows = lambda p: list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def mse(v):
    v = np.array([x for x in v if math.isfinite(x)])
    return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size)), int(v.size)) if v.size > 1 else (float("nan"), float("nan"), 0)


def pool(d, path):
    R = {(r["run"], r["split"]): r for r in rows(d / "regions_splits.tsv") if r["path"] == path and r["kind"] == "ep" and r["region"] == "pool"}
    S = {(r["run"], r["split"]): r for r in rows(d / "splits.tsv") if r["path"] == path and r["variant"] == "all"}
    return R, S


def pair(cd, cp, vd, vp):
    """dict of (control mean, variant mean, d, SE, n) per quantity, plus pile-up flag and yield."""
    Rc, Sc = pool(cd, cp)
    Rv, Sv = pool(vd, vp)
    ks = sorted(set(Rc) & set(Rv))
    out = {}
    for q in Q:
        c = [float(Rc[k][q]) for k in ks]
        v = [float(Rv[k][q]) for k in ks]
        ok = [i for i in range(len(ks)) if math.isfinite(c[i]) and math.isfinite(v[i])]
        out[q] = (mse([c[i] for i in ok])[0], mse([v[i] for i in ok])[0], *mse([v[i] - c[i] for i in ok]))
    pu = [(float(Sc[k]["pu"]), float(Sv[k]["pu"])) for k in ks if k in Sc and k in Sv]
    out["pu"] = (mse([a for a, b in pu])[0], mse([b for a, b in pu])[0], *mse([b - a for a, b in pu]))
    yl = [(float(Sc[k]["yield_n"]), float(Sv[k]["yield_n"])) for k in ks if k in Sc and k in Sv]
    out["yield"] = (100.0 * sum(b for a, b in yl) / sum(a for a, b in yl) - 100.0,)
    return out


def verdict(d, se):
    return "real" if math.isfinite(d) and math.isfinite(se) and abs(d) >= 2 * se else "null"


def main():
    rts = sys.argv[1:] or list(GATE)
    out = []
    for rt in rts:
        vd = O / "variants/P17" / rt
        if not (vd / "regions_splits.tsv").exists():
            continue
        paths = {r["path"] for r in rows(vd / "regions_splits.tsv")}
        S = {r["path"]: r for r in rows(vd / "summary.tsv")}
        tags = sorted((t, ASYM[t]) for t in ASYM if f"g{t:g}" in paths)
        for t, (e, l) in tags:
            p = f"g{t:g}"
            m = pair(CTL / rt, GATE[rt], vd, p)
            if rt == "2p24GeV_Ta_X17" and (O / "variants/P17/sym_extra/2p24GeV_Ta_X17/regions_splits.tsv").exists():     # also against the symmetric 4 ns gate
                m4 = pair(O / "variants/P17/sym_extra/2p24GeV_Ta_X17", "g4", vd, p)
                s4 = {f"s4_d_{q}": m4[q][2] for q in ("res_cb_oos", "tail_hi3", "leak_cut", "pu")} | {f"s4_d_{q}_se": m4[q][3] for q in ("res_cb_oos", "tail_hi3", "leak_cut", "pu")}
            else:
                s4 = {}
            out.append(dict(run_type=rt, early=e, late=l, control=GATE[rt], n=m["res_cb_oos"][4], **s4, **{f"{q}_ctl": m[q][0] for q in Q},
                            **{f"{q}": m[q][1] for q in Q}, **{f"d_{q}": m[q][2] for q in Q}, **{f"d_{q}_se": m[q][3] for q in Q},
                            pu_ctl=m["pu"][0], pu=m["pu"][1], d_pu=m["pu"][2], d_pu_se=m["pu"][3], yield_vs_ctl=m["yield"][0],
                            yield_vs_legacy=float(S[p]["yield_pct"])))
    sx = O / "variants/P17/sym_extra/2p24GeV_Ta_X17"
    if (sx / "regions_splits.tsv").exists() and (not sys.argv[1:] or "2p24GeV_Ta_X17" in rts):
        for g in (3, 4):
            m = pair(CTL / "2p24GeV_Ta_X17", "g8", sx, f"g{g}")
            out.append(dict(run_type="2p24GeV_Ta_X17", early=float(g), late=float(g), control="g8", n=m["res_cb_oos"][4], symmetric=1,
                            **{f"{q}_ctl": m[q][0] for q in Q}, **{f"{q}": m[q][1] for q in Q}, **{f"d_{q}": m[q][2] for q in Q},
                            **{f"d_{q}_se": m[q][3] for q in Q}, pu_ctl=m["pu"][0], pu=m["pu"][1], d_pu=m["pu"][2], d_pu_se=m["pu"][3],
                            yield_vs_ctl=m["yield"][0], yield_vs_legacy=float("nan")))
    if not out:
        sys.exit("no P17 reports yet")
    head = list(dict.fromkeys(k for r in out for k in r))
    with open(IMP / "P17_asym_window.tsv", "w") as f:
        f.write("\t".join(head) + "\n")
        for r in out:
            f.write("\t".join(f"{r.get(h, ''):.6g}" if isinstance(r.get(h), float) else str(r.get(h, "")) for h in head) + "\n")
    for rt in dict.fromkeys(r["run_type"] for r in out):
        print(f"\n## {rt} (control {GATE[rt]}; d = variant - control; oos res_cb in %, SE; tails/leak/eff in %; pu = pile-up flag)")
        print("| early | late | res_cb oos | d oos | d in | d tail_hi3 | d tail_lo_f | d leak | d eff | d pu | yield vs ctl % |\n|---|---|---|---|---|---|---|---|---|---|---|")
        for r in sorted((r for r in out if r["run_type"] == rt), key=lambda r: (r["late"], r["early"])):
            f = lambda q, n=3: f"{r['d_' + q]:+.{n}f} ± {r['d_' + q + '_se']:.{n}f}"
            print(f"| {r['early']:g}{'s' if r.get('symmetric') else ''} | {r['late']:g} | {r['res_cb_oos']:.3f} | {f('res_cb_oos')} | {f('res_cb')} | {f('tail_hi3')} | {f('tail_lo_f', 2)} "
                  f"| {f('leak_cut', 2)} | {f('eff_cut', 2)} | {f('pu', 2)} | {r['yield_vs_ctl']:+.2f} |")


main()
