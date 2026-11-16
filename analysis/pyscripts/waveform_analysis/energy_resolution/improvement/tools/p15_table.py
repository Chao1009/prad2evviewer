#!/usr/bin/env python3
"""p15_table.py: P15 tables (calibrated time-gate scan) from output/energy_resolution/variants/P15/<run type>/summary.tsv.
Writes energy_resolution/improvement/P15_tgate_calibrated.tsv and the tables printed to stdout (pasted into the .md).  Usage: p15_table.py"""
import csv, math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[4]
OUT = ROOT / "output/energy_resolution"
ER = Path(__file__).resolve().parents[1]
RTS = ["0p73GeV_H2", "2p24GeV_H2", "2p24GeV_Ta_X17", "3p49GeV_H2", "3p49GeV_H2_old"]
YIELD_MIN = 99.0
GATES = [1, 2, 3, 4, 5, 6, 8, 10, 12, 16, 20, 24]


def rd(p):
    with open(p) as f:
        return {r["path"]: r for r in csv.DictReader((l for l in f if not l.startswith("#")), delimiter="\t")}


def fl(r, k):
    try:
        return float(r[k])
    except (TypeError, ValueError, KeyError):
        return float("nan")


def best(S):
    """ERB summary() rule: smallest gate within 1 SE of the out-of-sample minimum (yield >= YIELD_MIN); and the plain minimum."""
    gs = sorted((float(p[1:]), p) for p in S if p != "legacy")
    ok = [(g, p) for g, p in gs if fl(S[p], "yield_pct") >= YIELD_MIN and math.isfinite(fl(S[p], "res_cb_oos")) and math.isfinite(fl(S[p], "res_cb_oos_se"))]
    lo = min(ok, key=lambda x: fl(S[x[1]], "res_cb_oos"))
    rule = next(x for x in ok if fl(S[x[1]], "res_cb_oos") <= fl(S[lo[1]], "res_cb_oos") + fl(S[lo[1]], "res_cb_oos_se"))
    return rule, lo


def main():
    tsv, md, bt = [], [], []
    for rt in RTS:
        C = rd(OUT / "variants/P15" / rt / "summary.tsv")
        R = rd(ER / ".." / rt / "summary.tsv")
        B = rd(OUT / "baseline_20261112/pre" / rt / "summary.tsv")
        cont = {k: (fl(C["legacy"], k), fl(B["legacy"], k)) for k in ("res_cb", "res_cb_oos", "n_used")}
        L, Lo = fl(C["legacy"], "res_cb"), fl(C["legacy"], "res_cb_oos")
        md.append(f"\n### {rt}\n\nLegacy continuity (P15 vs baseline pre/): " + ", ".join(f"{k} {a:.5g} vs {b:.5g}" for k, (a, b) in cont.items()) + "\n")
        md.append("| gate (ns) | res_cb in | SE | res_cb oos | SE | tail | pile-up | yield % | d in vs legacy | d oos vs legacy | raw oos (12 Nov) | cal - raw oos |")
        md.append("|---|---|---|---|---|---|---|---|---|---|---|---|")
        gates = [g for g in GATES if f"g{g:g}" in C]
        for p in ["legacy"] + [f"g{g:g}" for g in gates]:
            r = C[p]
            rr = R.get(p) or (B.get(p) if p in B else None)
            raw = fl(rr, "res_cb_oos") if rr else float("nan")
            o = fl(r, "res_cb_oos")
            row = dict(run_type=rt, path=p, res_cb=fl(r, "res_cb"), res_cb_se=fl(r, "res_cb_se"), res_cb_oos=o, res_cb_oos_se=fl(r, "res_cb_oos_se"),
                       tail=fl(r, "tail"), pu=fl(r, "pu"), yield_pct=fl(r, "yield_pct"), d_in=fl(r, "res_cb") - L, d_oos=o - Lo, raw_oos=raw, cal_minus_raw_oos=o - raw)
            tsv.append(row)
            md.append(f"| {p[1:] if p != 'legacy' else 'legacy'} | {row['res_cb']:.3f} | {row['res_cb_se']:.3f} | {o:.3f} | {row['res_cb_oos_se']:.3f} | {row['tail']:.2f} | {row['pu']:.2f} | {row['yield_pct']:.2f} | {row['d_in']:+.3f} | {row['d_oos']:+.3f} | {raw:.3f} | {o - raw:+.3f} |")
        (cg, cp), (mg, mp) = best({p: C[p] for p in C if p == "legacy" or float(p[1:]) in gates})
        (rg, rp), (rmg, rmp) = best(R)
        co, ro = fl(C[cp], "res_cb_oos"), fl(R[rp], "res_cb_oos")
        bt.append((rt, cg, co, fl(C[cp], "res_cb_oos_se"), mg, fl(C[mp], "res_cb_oos"), rg, ro, rmg, fl(R[rmp], "res_cb_oos"), Lo, co - Lo, co - ro))
    md.append("\n### Best gate per run type (out-of-sample res_cb, %)\n")
    md.append("| run type | cal rule-best (ns) | res_cb oos | SE | cal scan min (ns) | res_cb oos | raw rule-best (ns) | res_cb oos | raw scan min (ns) | res_cb oos | legacy oos | cal best - legacy | cal best - raw best |")
    md.append("|---|---|---|---|---|---|---|---|---|---|---|---|---|")
    for b in bt:
        md.append(f"| {b[0]} | {b[1]:g} | {b[2]:.3f} | {b[3]:.3f} | {b[4]:g} | {b[5]:.3f} | {b[6]:g} | {b[7]:.3f} | {b[8]:g} | {b[9]:.3f} | {b[10]:.3f} | {b[11]:+.3f} | {b[12]:+.3f} |")
    print("\n".join(md))
    with open(ER / "P15_tgate_calibrated.tsv", "w") as f:
        w = csv.DictWriter(f, fieldnames=list(tsv[0]), delimiter="\t")
        w.writeheader()
        for r in tsv:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})


if __name__ == "__main__":
    main()
