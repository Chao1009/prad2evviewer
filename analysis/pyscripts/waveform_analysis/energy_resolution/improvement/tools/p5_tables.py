"""P5 tables from the phase-1 outputs (no fits of its own). Run from analysis/pyscripts:
python3 waveform_analysis/energy_resolution/improvement/tools/p5_tables.py
Writes improvement/P5_latest.tsv (legacy vs best calibrated gate, P15), P5_budget.tsv (one row per term and run type:
paired d res_cb +- SE, quadrature size q = sqrt(R^2 - (R + d)^2)
(applied corrections: d = cost of switching off, q = sqrt((R + d)^2 - R^2)) +- SE, status, source) and P5_separation.tsv (U7 columns per
region and path, P2 regions plus the P15 best-gate pool row). Quadrature sizes are taken against the control's pool res_cb."""
import csv, math, sys
from pathlib import Path
import numpy as np

O = Path("output/energy_resolution")
IMP = Path("waveform_analysis/energy_resolution/improvement")
RTS = ("0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old", "2p24GeV_Ta_X17")
SEP = ("n_range", "res_cb", "res_cb_oos", "sigma_e", "sep_pi", "tail_lo3", "tail_hi3", "tail_lo_f", "tail_hi_f", "eff_cut", "leak_cut")
REGIONS = ("pool", "pool_centre", "pool_edge", "I1", "I2", "M", "O", "T", "H")


def rows(p):
    return list(csv.DictReader((l for l in open(p) if not l.startswith("#")), delimiter="\t")) if Path(p).exists() else []


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def mse(v):
    v = np.array([x for x in v if math.isfinite(x)])
    return (float(v.mean()), float(v.std(ddof=1) / math.sqrt(v.size))) if v.size > 1 else (float("nan"),) * 2


def best_gate(rt):
    """P15 best gate: the path named in P15_tgate_calibrated.tsv (column best_path) if present, else the lowest res_cb_oos."""
    for r in rows(IMP / "P15_tgate_calibrated.tsv"):
        if r.get("run_type", r.get("RT")) == rt and r.get("best_path"):
            return r["best_path"]
    s = [r for r in rows(O / f"variants/P15/{rt}/summary.tsv") if r["path"].startswith("g")]
    return min(s, key=lambda r: f(r["res_cb_oos"]))["path"] if s else None


def paired(ctrl, var, q="res_cb", pc="legacy", pv="legacy"):
    """Paired variant - control over (run, split) of two splits.tsv files (variant 'all')."""
    c = {(r["run"], r["split"]): f(r[q]) for r in rows(ctrl) if r["path"] == pc and r["variant"] == "all"}
    v = {(r["run"], r["split"]): f(r[q]) for r in rows(var) if r["path"] == pv and r["variant"] == "all"}
    d = [v[k] - c[k] for k in c if k in v and math.isfinite(v[k] - c[k])]
    return mse(d), len(d)


def summ(d, path="legacy"):
    return next((r for r in rows(Path(d) / "summary.tsv") if r["path"] == path), {})


def latest():
    out = []
    for rt in RTS:
        p15, g = O / f"variants/P15/{rt}", best_gate(rt)
        if g is None:
            continue
        lg, bg = summ(p15), summ(p15, g)
        (d, se), n = paired(p15 / "splits.tsv", p15 / "splits.tsv", pv=g)
        (do, seo), _ = paired(p15 / "splits.tsv", p15 / "splits.tsv", "res_cb_oos", pv=g)
        out.append(dict(run_type=rt, legacy=lg["res_cb"], legacy_se=lg["res_cb_se"], legacy_oos=lg["res_cb_oos"], best_path=g,
                        best=bg["res_cb"], best_se=bg["res_cb_se"], best_oos=bg["res_cb_oos"], d=d, d_se=se, d_oos=do,
                        d_oos_se=seo, n_pairs=n, yield_pct=bg["yield_pct"]))
    return out


# (term, run type, control report dir, variant report dir, control path, variant path, status note, source)
def terms():
    t = []
    pre, P = O / "baseline_20261112/pre", O / "variants"
    for rt in RTS:
        t.append(("log_weight_thres R7ri (P4b)", rt, pre / rt, P / f"P4/b_logw_r7ri/{rt}", "g12", "g12", "lever", "P4_switches.tsv"))
        t.append(("Moller rejection (P1a)", rt, pre / rt, P / f"P1/a_moller_reject/{rt}", "legacy", "legacy", "lever", "P1_reconcile.tsv"))
        t.append(("non-linearity (P1g)", rt, pre / rt, P / f"P1/g_nonlinear/{rt}", "legacy", "legacy", "lever", "P1_reconcile.tsv"))
        g = best_gate(rt)
        if g:
            t.append((f"time gate {g} calibrated (P15)", rt, P / f"P15/{rt}", P / f"P15/{rt}", "legacy", g, "lever", "P15"))
    sv = O / "scan_variants"
    c0 = O / "baseline_20261112/c0/2p24GeV_H2"
    t += [("onset-phase correction (P4a, applied)", "2p24GeV_H2", c0, P / "P4/a_phase_off/2p24GeV_H2", "legacy", "legacy", "applied", "P4 a, run 25320"),
          ("MIN_SEP_NS 12 (P3B)", "2p24GeV_H2", c0, sv / "minsep12/report/2p24GeV_H2", "legacy", "legacy", "lever", "P3 B, run 25320"),
          ("R_MODEL 0.3 (P3B)", "2p24GeV_H2", c0, sv / "rmodel0p3/report/2p24GeV_H2", "legacy", "legacy", "lever", "P3 B, run 25320")]
    for rt in ("2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"):
        for path in (("legacy", "g8") if rt.endswith("X17") else ("legacy", "g12")):
            t.append((f"R_MODEL 0.3, D5 extension", rt, P / f"P5/ext_ctrl/{rt}", P / f"P5/ext_rmodel/{rt}", path, path, "lever", "P5 extension"))
    for rt in ("2p24GeV_H2", "3p49GeV_H2"):
        t.append(("R_MODEL 0.3 on calibrated g16 + log weight", rt, P / f"P5/toff_logw_2run/{rt}" if rt.startswith("3p49") else P / f"P5/toff_logw/{rt}",
                  P / f"P5/toff_logw_rmodel/{rt}", "g16", "g16", "lever", "P5 combination"))
    for rt in ("2p24GeV_H2", "3p49GeV_H2", "3p49GeV_H2_old"):
        t.append(("log weight on calibrated g16", rt, P / f"P15/{rt}", P / f"P5/toff_logw/{rt}", "g16", "g16", "lever", "P5 combination"))
    return t


def budget():
    out = []
    for term, rt, cd, vd, cp, vp, status, src in terms():
        if not (Path(cd) / "splits.tsv").exists() or not (Path(vd) / "splits.tsv").exists():
            continue
        (d, se), n = paired(Path(cd) / "splits.tsv", Path(vd) / "splits.tsv", pc=cp, pv=vp)
        R = f(summ(cd, cp).get("res_cb"))
        q2 = (R + d) ** 2 - R * R if status == "applied" else R * R - (R + d) ** 2   # applied: d = cost of switching it off
        q = math.copysign(math.sqrt(abs(q2)), q2)
        q_se = abs(R + d) * se / abs(q) if q else float("nan")
        verdict = "real" if se > 0 and abs(d) >= 2 * se else "null"
        out.append(dict(term=term, run_type=rt, path=vp, n_pairs=n, R_ctrl=R, d=d, d_se=se, q=q, q_se=q_se, verdict=verdict,
                        status=status, source=src))
    return out


def sep_rows(rt, d, src, paths, regions=REGIONS):
    """U7 columns per (path, region): mean +- SE over splits of regions_splits.tsv, status from regions.tsv."""
    st = {(r["path"], r["region"]): r["status"] for r in rows(Path(d) / "regions.tsv") if r["kind"] == "ep"}
    by = {}
    for r in rows(Path(d) / "regions_splits.tsv"):
        if r["kind"] == "ep" and r["path"] in paths and r["region"] in regions:
            by.setdefault((r["path"], r["region"]), []).append(r)
    out = []
    for (p, reg), rs in sorted(by.items(), key=lambda kv: (kv[0][0], regions.index(kv[0][1]))):
        row = dict(run_type=rt, source=src, path=p, region=reg, status=st.get((p, reg), ""), n_range=sum(f(r["n_range"]) for r in rs))
        for k in SEP[1:]:
            if k.endswith("_se"):
                continue
            m, se = mse([f(r[k]) for r in rs])
            row[k], row[k + "_se"] = m, se
        out.append(row)
    return out


def separation():
    out, P = [], O / "variants"
    for rt in RTS:
        out += sep_rows(rt, P / f"P2/{rt}", "P2 raw times", ("legacy", "g12", "g8"))
        out += sep_rows(rt, P / f"P15/{rt}", "P15 calibrated times", (best_gate(rt),), ("pool",))
        for tag, src in (("toff_logw", "P5 calibrated g16 + log weight"), ("toff_logw_rmodel", "P5 calibrated g16 + log weight + R_MODEL 0.3")):
            out += sep_rows(rt, P / f"P5/{tag}/{rt}", src, ("g16",), ("pool",))
    return out


def write(name, data):
    if not data:
        print(f"[skip] {name}: no rows", file=sys.stderr)
        return
    with open(IMP / name, "w") as fh:
        w = csv.DictWriter(fh, fieldnames=list(data[0]), delimiter="\t", lineterminator="\n")
        w.writeheader()
        for r in data:
            w.writerow({k: (f"{v:.5g}" if isinstance(v, float) else v) for k, v in r.items()})
    print(f"[done] {IMP / name}: {len(data)} rows")


if __name__ == "__main__":
    write("P5_latest.tsv", latest())
    write("P5_budget.tsv", budget())
    write("P5_separation.tsv", separation())
