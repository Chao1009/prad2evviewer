"""Tables for the phase-2 control v3 (v3 = v2 + P16 threshold correction + M1 + P7).  Run from analysis/pyscripts:
  python3 waveform_analysis/energy_resolution/improvement/tools/p2v3_table.py pair <cdir> <cpath> <vdir> <vpath> [region]   (one run type's dirs; prints one line)
  python3 waveform_analysis/energy_resolution/improvement/tools/p2v3_table.py control  (v3 control table -> improvement/control_phase2_v3.tsv)
  python3 waveform_analysis/energy_resolution/improvement/tools/p2v3_table.py delta    (v3 vs v2, pool -> improvement/control_phase2_v3_delta.tsv)
  python3 waveform_analysis/energy_resolution/improvement/tools/p2v3_table.py steps    (step attribution v2 -> +P16 -> +M1 -> +P7 -> improvement/control_phase2_v3_attr.tsv)
  python3 waveform_analysis/energy_resolution/improvement/tools/p2v3_table.py gates    (every v3 path vs the v3 gate path -> improvement/control_phase2_v3_gates.tsv)
  python3 waveform_analysis/energy_resolution/improvement/tools/p2v3_table.py final    (legacy / v1 / v2 / v3 in-sample and oos -> improvement/control_phase2_v3_final.tsv)
Any subcommand takes `--r3 <dir>` to override the v3 root.  Missing dirs are skipped with a [skip] note.
Pairs per (run, split) the ep rows of a region; d = variant - control, mean +- SE, real if |d| >= 2 SE (D3)."""
import csv, math, sys
from functools import lru_cache
from pathlib import Path
import numpy as np

O = Path("output/energy_resolution")
IMP = Path("waveform_analysis/energy_resolution/improvement")
BASE, V1, V2 = O / "baseline_20261112" / "pre", O / "control_phase2_20261113", O / "control_phase2_v2_20261113"
R3 = O / "control_phase2_v3_20261114"
if "--r3" in sys.argv:
    i = sys.argv.index("--r3")
    R3 = Path(sys.argv[i + 1])
    del sys.argv[i:i + 2]
GATE_V1 = {"0p73GeV_H2": "g20", "2p24GeV_H2": "g16", "3p49GeV_H2": "g16", "3p49GeV_H2_old": "g16", "2p24GeV_Ta_X17": "g8"}
GATE_V2 = GATE_V1 | {"2p24GeV_Ta_X17": "g4"}
GATE_V3 = dict(GATE_V2)
RTS = tuple(GATE_V1)
Q = ("res_cb", "res_cb_oos", "tail_hi_f_oos", "tail_lo_f_oos", "tail_hi3_oos", "leak_cut_oos", "eff_cut_oos", "mu_oos", "sep_pi_oos",
     "tail_hi_f", "tail_lo_f", "leak_cut", "eff_cut")
QP = (*Q, "pu")
P7CTRL = O / "variants" / "P7" / "ctrl"          # control v2 with EVAL_REGIONS (ring rows of v2; pool rows equal v2)
STEPS = {rt: [("v2", V2, GATE_V2[rt]), ("+P16", R3 / "s1_p16", GATE_V3[rt]), ("+M1", R3 / "s2_m1", GATE_V3[rt]), ("+P7", R3 / "v3", GATE_V3[rt])] for rt in RTS}
# Ta: its own P16 table fails D3 (eff_cut, tail_lo_f), so its chain is v2 -> 2.24 GeV correction under P16 (no Ta table) -> + Ta M1 -> v3 (no P7 ring)
STEPS["2p24GeV_Ta_X17"] = [("v2", V2, "g4"), ("+P16 (2.24 map only)", R3 / "s1_ta_nop16", "g4"), ("+M1", R3 / "s2_ta_nop16", "g4"), ("+P7 (none)", R3 / "v3", "g4")]
NAN = float("nan")


@lru_cache(None)
def rows(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t"))


def fl(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return NAN


def mse(v):
    v = np.array([x for x in v if math.isfinite(x)])
    return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size)), int(v.size)) if v.size > 1 else (NAN, NAN, int(v.size))


def verdict(d, se):
    return "real" if math.isfinite(d) and math.isfinite(se) and d != 0 and abs(d) >= 2 * se else "null"


def ep(d, path, region):
    return {(r["run"], r["split"]): r for r in rows(d / "regions_splits.tsv") if r["path"] == path and r["kind"] == "ep" and r["region"] == region}


def spl(d, path):
    return {(r["run"], r["split"]): r for r in rows(d / "splits.tsv") if r["path"] == path and r["variant"] == "all"}


def pair(cd, cp, vd, vp, region="pool"):
    """{quantity: (control mean, variant mean, d, SE, n)} for Q; plus pu (pool only), yield (pool only, variant/control yield_n % change, 1-tuple)
    and gap ((variant oos-in) - (control oos-in) of res_cb: (d, SE, n)).  None when either side has no rows of the region."""
    Rc, Rv = ep(cd, cp, region), ep(vd, vp, region)
    ks = sorted(set(Rc) & set(Rv))
    if not ks:
        return None
    out = {}
    for q in Q:
        c = np.array([fl(Rc[k][q]) for k in ks])
        v = np.array([fl(Rv[k][q]) for k in ks])
        ok = np.isfinite(c) & np.isfinite(v)
        out[q] = (mse(c[ok])[0], mse(v[ok])[0], *mse(v[ok] - c[ok]))
    g = [(fl(Rv[k]["res_cb_oos"]) - fl(Rv[k]["res_cb"])) - (fl(Rc[k]["res_cb_oos"]) - fl(Rc[k]["res_cb"])) for k in ks]
    out["gap"] = mse(g)
    if region == "pool":
        Sc, Sv = spl(cd, cp), spl(vd, vp)
        pu = [(fl(Sc[k]["pu"]), fl(Sv[k]["pu"])) for k in ks if k in Sc and k in Sv]
        pu = [(a, b) for a, b in pu if math.isfinite(a) and math.isfinite(b)]
        out["pu"] = (mse([a for a, b in pu])[0], mse([b for a, b in pu])[0], *mse([b - a for a, b in pu]))
        yl = [(fl(Sc[k]["yield_n"]), fl(Sv[k]["yield_n"])) for k in ks if k in Sc and k in Sv]
        out["yield"] = (100.0 * sum(b for a, b in yl) / sum(a for a, b in yl) - 100.0,) if yl and sum(a for a, b in yl) else (NAN,)
    else:
        out["pu"] = (NAN,) * 4 + (0,)
        out["yield"] = (NAN,)
    return out


def line(m):
    return ("  ".join(f"{q} {m[q][2]:+.4f}+-{m[q][3]:.4f}" for q in QP) + f"  yield {m['yield'][0]:+.2f}%"
            + f"  gap {m['gap'][0]:+.4f}+-{m['gap'][1]:.4f}  n={m['res_cb'][4]}")


def dcols(m):
    """d_q, d_q_se, verdict_q for Q + pu, yield_pct, d_gap, d_gap_se."""
    row = {}
    for q in QP:
        row |= {f"d_{q}": m[q][2], f"d_{q}_se": m[q][3], f"verdict_{q}": verdict(m[q][2], m[q][3])}
    row |= {"yield_pct": m["yield"][0], "d_gap": m["gap"][0], "d_gap_se": m["gap"][1]}
    return row


def write(name, out):
    if not out:
        print(f"[skip] nothing to write for {name}")
        return
    with open(IMP / name, "w") as f:
        f.write("\t".join(out[0]) + "\n")
        f.writelines("\t".join(f"{v:.6g}" if isinstance(v, float) else str(v) for v in r.values()) + "\n" for r in out)
    print(f"wrote {IMP / name}")


def have(d, *names):
    ok = all((d / n).exists() for n in names)
    if not ok:
        print(f"[skip] {d} missing")
    return ok


def summ(d, path):
    return {r["path"]: r for r in rows(d / "summary.tsv")}.get(path)


def control_table():
    out = []
    for rt in RTS:
        d, g = R3 / "v3" / rt, GATE_V3[rt]
        if not have(d, "summary.tsv", "regions_splits.tsv"):
            continue
        s = summ(d, g)
        if s is None:
            print(f"[skip] {d}: no path {g}")
            continue
        R = list(ep(d, g, "pool").values())
        row = dict(run_type=rt, gate=g, n_splits=s["n_splits"], n_used=s["n_used"])
        for q in ("res_cb", "res_cb_oos"):
            row[q], row[q + "_se"] = float(s[q]), float(s[q + "_se"])
        for sfx in ("", "_oos"):
            for q in ("tail_hi_f", "tail_lo_f", "tail_hi3", "leak_cut", "eff_cut", "sep_pi"):
                row[q + sfx], row[q + sfx + "_se"], _ = mse([fl(r[q + sfx]) for r in R])
        row["pu_pct"], row["pu_se"], row["yield_pct"] = float(s["pu"]), float(s["pu_se"]), float(s["yield_pct"])
        out.append(row)
        print(f"| {rt} | {g} | {row['res_cb']:.4f} +- {row['res_cb_se']:.4f} | {row['res_cb_oos']:.4f} +- {row['res_cb_oos_se']:.4f} | "
              f"{row['tail_hi_f_oos']:.4f} | {row['tail_lo_f_oos']:.2f} | {row['leak_cut_oos']:.3f} | {row['eff_cut_oos']:.2f} | {row['pu_pct']:.2f} | {row['yield_pct']:.2f} |")
    write("control_phase2_v3.tsv", out)


def delta_table():
    out = []
    for rt in RTS:
        cd, vd = V2 / rt, R3 / "v3" / rt
        if not (have(cd, "regions_splits.tsv", "splits.tsv") and have(vd, "regions_splits.tsv", "splits.tsv")):
            continue
        m = pair(cd, GATE_V2[rt], vd, GATE_V3[rt])
        if m is None:
            print(f"[skip] {rt}: no paired rows")
            continue
        row = dict(run_type=rt, gate_v2=GATE_V2[rt], gate_v3=GATE_V3[rt], n_pairs=m["res_cb"][4], yield_pct=m["yield"][0])
        for q in QP:
            row |= {f"{q}_v2": m[q][0], f"{q}_v3": m[q][1], f"d_{q}": m[q][2], f"d_{q}_se": m[q][3], f"verdict_{q}": verdict(m[q][2], m[q][3])}
        row |= {"d_gap": m["gap"][0], "d_gap_se": m["gap"][1]}
        out.append(row)
        print(rt, line(m))
    write("control_phase2_v3_delta.tsv", out)


def steps_table():
    out = []

    def add(rt, region, a, b, m):
        out.append(dict(rt=rt, region=region, step_from=a, step_to=b, n_pairs=m["res_cb"][4], **dcols(m)))
        return m

    for rt in RTS:
        st = STEPS[rt]
        ok = [s for s in st if (s[1] / rt / "regions_splits.tsv").exists() and (s[1] / rt / "splits.tsv").exists()]
        for s in st:
            if s not in ok:
                print(f"[skip] {s[1] / rt} missing")
        for region in ("pool", "H", "M"):
            sq = ok if region == "pool" else [("v2 (P7 ctrl)", P7CTRL, s[2]) if s[1] == V2 else s for s in ok]   # v2 has no ring rows: v2 + EVAL_REGIONS (pool = v2)
            seq = [s for s in sq if ep(s[1] / rt, s[2], region)]
            if len(seq) < 2:
                continue
            ms = []
            for a, b in zip(seq, seq[1:]):
                m = pair(a[1] / rt, a[2], b[1] / rt, b[2], region)
                if m is not None:
                    ms.append(add(rt, region, a[0], b[0], m))
            tot = pair(seq[0][1] / rt, seq[0][2], seq[-1][1] / rt, seq[-1][2], region)
            if tot is not None:
                add(rt, region, seq[0][0], "total", tot)
            if region == "pool" and tot is not None:
                print(f"{rt} check: sum consecutive d_res_cb_oos {sum(m['res_cb_oos'][2] for m in ms):+.4f} (n={[m['res_cb_oos'][4] for m in ms]})"
                      f" vs total {tot['res_cb_oos'][2]:+.4f} (n={tot['res_cb_oos'][4]})")
                for r in out[-len(ms) - 1:]:
                    print(f"  {r['step_from']} -> {r['step_to']}  d_res_cb {r['d_res_cb']:+.4f}+-{r['d_res_cb_se']:.4f}  d_res_cb_oos {r['d_res_cb_oos']:+.4f}+-{r['d_res_cb_oos_se']:.4f}")
    write("control_phase2_v3_attr.tsv", out)


def gates_table():
    out = []
    for rt in RTS:
        d = R3 / "v3" / rt
        if not have(d, "summary.tsv", "regions_splits.tsv"):
            continue
        S = {r["path"]: r for r in rows(d / "summary.tsv")}
        g = GATE_V3[rt]
        if g not in S:
            print(f"[skip] {d}: no gate path {g}")
            continue
        Rg = ep(d, g, "pool")
        for p, s in S.items():
            if p == "legacy":
                continue
            Rp = ep(d, p, "pool")
            ks = sorted(set(Rg) & set(Rp))
            dd, se, n = mse([fl(Rp[k]["res_cb_oos"]) - fl(Rg[k]["res_cb_oos"]) for k in ks])
            out.append(dict(run_type=rt, path=p, gate=g, res_cb_oos=float(s["res_cb_oos"]), res_cb_oos_se=float(s["res_cb_oos_se"]),
                            d_vs_gate=dd, d_vs_gate_se=se, n_pairs=n, verdict=verdict(dd, se)))
            print(f"{rt} {p:>6} oos {float(s['res_cb_oos']):.4f}+-{float(s['res_cb_oos_se']):.4f}  vs {g} {dd:+.4f}+-{se:.4f} {verdict(dd, se)} n={n}")
    write("control_phase2_v3_gates.tsv", out)


def final_table():
    stages = [("legacy", lambda rt: (BASE / rt, "legacy")), ("v1", lambda rt: (V1 / rt, GATE_V1[rt])),
              ("v2", lambda rt: (V2 / rt, GATE_V2[rt])), ("v3", lambda rt: (R3 / "v3" / rt, GATE_V3[rt]))]
    out, cell = [], {}
    for rt in RTS:
        for name, f in stages:
            d, g = f(rt)
            if not have(d, "summary.tsv"):
                continue
            s = summ(d, g)
            if s is None:
                print(f"[skip] {d}: no path {g}")
                continue
            r = dict(rt=rt, stage=name, gate=g, n_splits=s["n_splits"])
            for q in ("res_cb", "res_cb_se", "res_cb_oos", "res_cb_oos_se"):
                r[q] = float(s[q])
            out.append(r)
            cell[rt, name] = f"{r['res_cb']:.3f}+-{r['res_cb_se']:.3f} / {r['res_cb_oos']:.3f}+-{r['res_cb_oos_se']:.3f}"
    print("| run type | " + " | ".join(n for n, _ in stages) + " |\n|---|" + "---|" * len(stages))
    for rt in RTS:
        print(f"| {rt} | " + " | ".join(cell.get((rt, n), "-") for n, _ in stages) + " |")
    write("control_phase2_v3_final.tsv", out)


if __name__ == "__main__":
    cmd = sys.argv[1:2]
    if cmd == ["control"]:
        control_table()
    elif cmd == ["delta"]:
        delta_table()
    elif cmd == ["steps"]:
        steps_table()
    elif cmd == ["gates"]:
        gates_table()
    elif cmd == ["final"]:
        final_table()
    elif cmd == ["pair"] and len(sys.argv) >= 6:
        m = pair(Path(sys.argv[2]), sys.argv[3], Path(sys.argv[4]), sys.argv[5], *sys.argv[6:7])
        print(line(m) if m else "[skip] no paired rows")
    else:
        print(__doc__)
