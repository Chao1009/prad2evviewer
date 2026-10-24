#!/usr/bin/env python3
"""compare_scan_dirs.py - pair driver-layout scans per (run, split) and tabulate sigma_E/E.

    python compare_scan_dirs.py --scans A=output/scan_20261023_pc B=output/scan_20261024_spl \
        --out output/scan_20261024_spl/compare/arms

Inputs: scan directories written by benchmark_hycal_timing_scan.py
(<dir>/scan_<run:06d>/split_<n:05d>/result.tsv), read through
gate_scan_20261022.load_result.  Only (run, split) pairs present in EVERY scan are used
(the rest are listed in the report); runs in --exclude-runs (default 24578) are skipped.

Quantities (calibrated-block sigma/E, fractions in the TSVs are printed in percent):
  post  = per-crystal gains with the seed-indexed refinement   (sigma_post, sigma_cb)
  pre   = pooled gains only, no seed-indexed step              (sigma_pre,  sigma_cb_pre)
  fit   = the run type's preferred fit (gate_scan_20261022.R_PREF): CB for 3.49 GeV H2,
          2.24 GeV C and Ta; Gauss for 2.24 and 0.73 GeV H2
  path  = legacy and gated ("new" in result.tsv), reported side by side.
Per run type the mean is over splits; delta = cmp - ref is paired per split with
SE = SD/sqrt(n) over splits (SD = the split-to-split SD of delta).

Arms (plan_20261024 §8.3), by position: scan 1 post = A, scan 1 pre = A', scan 2 post = B,
scan 2 pre = C (further scans: <label> and <label>'). Rules, with the second scan as the
candidate:
  (1) choice B vs C per run type (legacy path): C iff sigma_B - sigma_C > 2 SE (C lower by
      more than 2 SE), else B; the same arm is then reported on the gated path;
  (2) template verdict, chosen arm vs A: 2.24 / 3.49 GeV H2: Delta <= -0.02 pp and
      |Delta| > 2 SE = success; |Delta| <= 2 SE = neutral; Delta > 0 by > 2 SE = investigate
      (also for 2.24 GeV C); Ta is reported only; other significant changes are listed as
      "better/worse (not graded)";
  (3) eps_yield, eps_peak (harness definitions, G.stats_of): eps_yield = sum of gated
      n_events of the scan / sum of legacy n_events of scan 1; eps_peak = sum of the
      peak_yield "new" / "legacy" counts within the scan.  Floor 0.99.

Outputs: PREFIX.md, PREFIX.tsv (per run type x path x comparison), PREFIX_splits.tsv
(one row per split and path, all four sigma kinds of every scan).
"""
from __future__ import annotations

import argparse
import math
import re
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import gate_scan_20261022 as G                       # noqa: E402
from benchmark_hycal_timing_scan import _run_type    # noqa: E402

PATHS = (("legacy", "legacy"), ("new", "gated"))     # (result.tsv path, report name)
KINDS = ("G_pre", "G_post", "CB_pre", "CB_post")
TYPE_ORDER = ["3.49 GeV H2", "2.24 GeV H2", "0.73 GeV H2", "2.24 GeV C", "2.24 GeV Ta (X17)"]
GRADED = ("3.49 GeV H2", "2.24 GeV H2")              # success rule (2)
WORSE_CHECK = ("3.49 GeV H2", "2.24 GeV H2", "0.73 GeV H2", "2.24 GeV C")   # "any H2 or C"
MIN_SUCCESS_PP = 0.02


# ---------------------------------------------------------------------------
# Loading
# ---------------------------------------------------------------------------
def load_scan(d: Path, exclude: set[int]) -> dict:
    out = {}
    for t in sorted(Path(d).glob("scan_*/split_*/result.tsv")):
        m1 = re.match(r"scan_0*(\d+)$", t.parent.parent.name)
        m2 = re.match(r"split_(\d+)$", t.parent.name)
        if not (m1 and m2):
            continue
        run, sp = int(m1.group(1)), int(m2.group(1))
        if run in exclude or _run_type(run) not in G.R_PREF:
            continue
        rec = G.load_result(t, str(d), run, sp)
        if rec is not None:
            out[(run, sp)] = rec
    return out


def pref_fit(run: int) -> str:
    return "cb" if G.R_PREF.get(_run_type(run), "gauss") == "cb" else "gauss"


def sig(rec, path, run, kind) -> float:
    """sigma/E in percent; kind in KINDS or 'pref_post' / 'pref_pre'."""
    if kind.startswith("pref_"):
        kind = ("CB_" if pref_fit(run) == "cb" else "G_") + kind[5:]
    return G.sigma(rec, path, kind) * 100


def calfield(rec, path, name) -> float:
    c = rec["cal"].get(path) or {}
    v = c.get(name, float("nan"))
    return float(v) if v is not None else float("nan")


def mean_sd_se(x):
    a = np.asarray([v for v in x if np.isfinite(v)], dtype=float)
    n = len(a)
    if n == 0:
        return float("nan"), float("nan"), float("nan"), 0
    if n == 1:
        return float(a[0]), float("nan"), float("nan"), 1
    sd = float(a.std(ddof=1))
    return float(a.mean()), sd, sd / math.sqrt(n), n


def fnum(x, nd=3, signed=False) -> str:
    if x is None or not np.isfinite(x):
        return "NA"
    s = f"{x:+.{nd}f}" if signed else f"{x:.{nd}f}"
    return s.replace("-", "−")


def dse(m, se, nd=3) -> str:
    return fnum(m, nd, True) + (f" ± {se:.{nd}f}" if np.isfinite(se) else "")


# ---------------------------------------------------------------------------
# Arms
# ---------------------------------------------------------------------------
def arm_defs(labels: list[str]) -> list[dict]:
    """Arms in report order: [{name, scan(index), kind('post'|'pre')}]."""
    names = ["A", "A'", "B", "C"]
    arms = []
    for i, lab in enumerate(labels):
        pn, qn = (names[2 * i], names[2 * i + 1]) if i < 2 else (lab, lab + "'")
        arms += [{"name": pn, "scan": i, "kind": "post"}, {"name": qn, "scan": i, "kind": "pre"}]
    return arms


def comparisons(arms: list[dict]) -> list[tuple[str, str]]:
    """(ref, cmp) pairs, delta = cmp - ref."""
    by = {a["name"]: a for a in arms}
    names = [a["name"] for a in arms]
    pairs = [("A", n) for n in names[1:]]                      # everything against A
    for a in arms:                                             # pre vs post within each scan
        if a["kind"] == "pre" and a["scan"] > 0:
            post = next(x["name"] for x in arms if x["scan"] == a["scan"] and x["kind"] == "post")
            pairs.append((a["name"], post))                    # delta = post - pre (B - C)
    return [p for p in pairs if p[0] in by and p[1] in by]


# ---------------------------------------------------------------------------
# Tables
# ---------------------------------------------------------------------------
def type_keys(keys, rt):
    return [k for k in keys if _run_type(k[0]) == rt]


def arm_sig(scans, labels, arm, key, path) -> float:
    rec = scans[labels[arm["scan"]]][key]
    return sig(rec, path, key[0], "pref_" + arm["kind"])


def eps_of(scans, labels, ks, i) -> tuple[float, float]:
    """eps_yield, eps_peak of scan i over splits ks (harness definitions)."""
    ref, cur = scans[labels[0]], scans[labels[i]]
    nn = sum(cur[k]["head"]["new"].get("n_events", float("nan")) for k in ks)
    nr = sum(ref[k]["head"]["legacy"].get("n_events", float("nan")) for k in ks)
    pn = sum(cur[k]["hdr"]["peak"]["new"] for k in ks)
    pl = sum(cur[k]["hdr"]["peak"]["legacy"] for k in ks)
    ey = nn / nr if (ks and nr > 0) else float("nan")
    ep = pn / pl if (ks and np.isfinite(pl) and pl > 0) else float("nan")
    return ey, ep


def type_table(scans, labels, arms, keys):
    """List of dict rows: per run type x path x comparison (ref, cmp)."""
    rows, armstat = [], {}
    for rt in TYPE_ORDER:
        ks = type_keys(keys, rt)
        if not ks:
            continue
        fit = "CB" if G.R_PREF[rt] == "cb" else "Gauss"
        slug = G.runtype_slug(ks[0][0])
        for p, pname in PATHS:
            for a in arms:
                v = [arm_sig(scans, labels, a, k, p) for k in ks]
                m, sd, se, n = mean_sd_se(v)
                lab = labels[a["scan"]]
                nmod = np.nanmean([calfield(scans[lab][k], p, "n_modules") for k in ks])
                nev = np.nanmean([calfield(scans[lab][k], p, "n_inner_events") for k in ks])
                ey, ep = eps_of(scans, labels, ks, a["scan"])
                armstat[(rt, pname, a["name"])] = dict(mean=m, se_mean=se, n=n, nmod=nmod,
                                                       nev=nev, eps_yield=ey, eps_peak=ep)
            for ref, cmp_ in comparisons(arms):
                ar = next(x for x in arms if x["name"] == ref)
                ac = next(x for x in arms if x["name"] == cmp_)
                d = [arm_sig(scans, labels, ac, k, p) - arm_sig(scans, labels, ar, k, p)
                     for k in ks]
                dm, dsd, dse_, n = mean_sd_se(d)
                sr, sc = armstat[(rt, pname, ref)], armstat[(rt, pname, cmp_)]
                rows.append(dict(type=rt, slug=slug, path=pname, fit=fit, n_splits=n,
                                 ref=ref, cmp=cmp_, ref_mean=sr["mean"], cmp_mean=sc["mean"],
                                 delta=dm, delta_sd=dsd, delta_se=dse_,
                                 nmod_ref=sr["nmod"], nmod_cmp=sc["nmod"],
                                 nev_ref=sr["nev"], nev_cmp=sc["nev"],
                                 eps_yield=sc["eps_yield"], eps_peak=sc["eps_peak"]))
    return rows, armstat


def decide(rows, arms):
    """Rules (1) and (2) per run type; needs >= 2 scans (arms B and C)."""
    names = [a["name"] for a in arms]
    if not all(n in names for n in ("B", "C")):
        return {}
    get = {(r["type"], r["path"], r["ref"], r["cmp"]): r for r in rows}
    out = {}
    for rt in TYPE_ORDER:
        bc = get.get((rt, "legacy", "C", "B"))
        if bc is None:
            continue
        d, se = bc["delta"], bc["delta_se"]
        pick_c = bool(np.isfinite(d) and np.isfinite(se) and d > 2 * se)   # B - C > 2 SE
        chosen = "C" if pick_c else "B"
        res = {"chosen": chosen, "d_bc": d, "se_bc": se}
        for pname in ("legacy", "gated"):
            r = get.get((rt, pname, "A", chosen))
            if r is None:
                continue
            res[pname] = r
        r = res.get("legacy")
        if r is not None:
            d, se = r["delta"], r["delta_se"]
            sig2 = np.isfinite(se) and abs(d) > 2 * se
            if rt == "2.24 GeV Ta (X17)":
                v = "report only (Ta, <= 34 modules)"
            elif not sig2:
                v = "neutral"
            elif d > 0:
                v = "investigate (worse by > 2 SE)" if rt in WORSE_CHECK else "worse (not graded)"
            elif rt in GRADED and d <= -MIN_SUCCESS_PP:
                v = "success"
            else:
                v = "better (not graded)" if rt not in GRADED else "better, |Delta| < 0.02 pp"
            res["verdict"] = v
        out[rt] = res
    return out


# ---------------------------------------------------------------------------
# Output
# ---------------------------------------------------------------------------
TSV_COLS = ["type", "slug", "path", "fit", "n_splits", "ref", "cmp", "ref_mean", "cmp_mean",
            "delta", "delta_sd", "delta_se", "nmod_ref", "nmod_cmp", "nev_ref", "nev_cmp",
            "eps_yield", "eps_peak"]


def write_tsv(path: Path, cols, rows, nd=5):
    with open(path, "w") as f:
        f.write("\t".join(cols) + "\n")
        for r in rows:
            f.write("\t".join(f"{r[c]:.{nd}f}" if isinstance(r[c], (float, np.floating))
                              else str(r[c]) for c in cols) + "\n")


def split_rows(scans, labels, keys):
    rows = []
    for k in keys:
        for p, pname in PATHS:
            r = dict(type=_run_type(k[0]), run=k[0], split=k[1], path=pname,
                     fit="CB" if pref_fit(k[0]) == "cb" else "Gauss")
            for lab in labels:
                rec = scans[lab][k]
                for kind in KINDS:
                    r[f"{lab}_{kind}"] = sig(rec, p, k[0], kind)
                r[f"{lab}_nmod"] = calfield(rec, p, "n_modules")
                r[f"{lab}_nev"] = calfield(rec, p, "n_inner_events")
                r[f"{lab}_n_events"] = rec["head"][p].get("n_events", float("nan"))
                r[f"{lab}_peak"] = rec["hdr"]["peak"][p]
            rows.append(r)
    cols = ["type", "run", "split", "path", "fit"] + [
        f"{lab}_{x}" for lab in labels for x in (*KINDS, "nmod", "nev", "n_events", "peak")]
    return rows, cols


def report(scans, labels, dirs, arms, keys, rows, armstat, dec, skipped, args) -> str:
    L = ["# Scan comparison", ""]
    L += ["| label | directory | splits loaded | paired |", "|---|---|---|---|"]
    for lab in labels:
        L.append(f"| {lab} | `{dirs[lab]}` | {len(scans[lab])} | {len(keys)} |")
    L.append("")
    if skipped:
        L += [f"Not paired (missing in at least one scan): {len(skipped)}: "
              + ", ".join(f"{r}.{s:05d}" for r, s in skipped[:30])
              + (" ..." if len(skipped) > 30 else ""), ""]
    L += ["Arms: " + "; ".join(
        f"**{a['name']}** = {labels[a['scan']]} {a['kind']}" for a in arms)
        + ". post = with the seed-indexed step, pre = pooled gains only.", "",
        "Preferred fit per run type: CB for 3.49 GeV H2, 2.24 GeV C, Ta; Gauss for 2.24 and "
        "0.73 GeV H2. sigma/E in percent, mean over splits; Δ = cmp − ref, paired per "
        "split, mean ± SE (SD/√n over splits) in pp, the split-to-split SD in brackets. "
        "Modules and events are the calibrated-block values per split.", ""]

    get = {(r["type"], r["path"], r["ref"], r["cmp"]): r for r in rows}
    for pname in ("legacy", "gated"):
        L += [f"## sigma_E/E by arm, {pname} path (%)", ""]
        head = ["Run type", "fit", "splits"] + [f"{a['name']}" for a in arms] + ["modules (A, B)",
                                                                              "events (A, B)"]
        L += ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
        for rt in TYPE_ORDER:
            if (rt, pname, arms[0]["name"]) not in armstat:
                continue
            s0 = armstat[(rt, pname, "A")]
            row = [rt, "CB" if G.R_PREF[rt] == "cb" else "Gauss", str(s0["n"])]
            row += [fnum(armstat[(rt, pname, a["name"])]["mean"]) for a in arms]
            sb = armstat.get((rt, pname, "B"), s0)
            row += [f"{s0['nmod']:.1f}, {sb['nmod']:.1f}",
                    f"{s0['nev']/1e3:.1f}k, {sb['nev']/1e3:.1f}k"]
            L.append("| " + " | ".join(row) + " |")
        L.append("")
        L += [f"### Paired Δ (pp), {pname} path", ""]
        head = ["Run type"] + [f"{c} − {r}" for r, c in comparisons(arms)]
        L += ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
        for rt in TYPE_ORDER:
            cells = []
            for r, c in comparisons(arms):
                x = get.get((rt, pname, r, c))
                if x is not None:
                    cells.append(f"{dse(x['delta'], x['delta_se'])} [{fnum(x['delta_sd'])}]")
            if cells:
                L.append("| " + " | ".join([rt] + cells) + " |")
        L.append("")

    if dec:
        L += ["## Choice and verdict (plan §8.3)", "",
              "Rule 1: B or C, decided on the legacy path (C iff B − C > 2 SE); the same arm "
              "is shown on the gated path. Rule 2: chosen arm vs A.", ""]
        head = ["Run type", "B − C (pp)", "chosen", "Δ chosen − A legacy (pp)",
                "Δ gated (pp)", "verdict", "eps_yield", "eps_peak"]
        L += ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
        for rt in TYPE_ORDER:
            if rt not in dec:
                continue
            d = dec[rt]
            lg, gt = d.get("legacy"), d.get("gated")
            ok = lambda v: ("" if not np.isfinite(v) or v >= G.EPS_MIN else " (< 0.99)")   # noqa: E731
            ey = armstat[(rt, "gated", d["chosen"])]["eps_yield"]
            ep = armstat[(rt, "gated", d["chosen"])]["eps_peak"]
            L.append("| " + " | ".join([
                rt, dse(d["d_bc"], d["se_bc"]), d["chosen"],
                dse(lg["delta"], lg["delta_se"]) if lg else "NA",
                dse(gt["delta"], gt["delta_se"]) if gt else "NA",
                d.get("verdict", "NA"), fnum(ey, 4) + ok(ey), fnum(ep, 4) + ok(ep)]) + " |")
        L += ["",
              "eps_yield = sum of gated n_events of the chosen arm's scan / sum of legacy "
              "n_events of scan 1; eps_peak = peak_yield new / legacy counts within the chosen "
              "arm's scan (harness definitions; pre and post arms of one scan share both). "
              "Floor 0.99 (rule 3). Success is graded at 2.24 and 3.49 GeV H2; 2.24 GeV C is "
              "checked only for a worse result; Ta is reported, not interpreted.", ""]

    two = len(scans) > 1
    head = ["run", "type", "splits"] + [a["name"] for a in arms] + (["\u0394 (B \u2212 A) pp"] if two else [])
    L += ["## Per run, preferred fit, legacy path (%)", "",
          "| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for run in sorted({k[0] for k in keys}):
        ks = [k for k in keys if k[0] == run]
        cells = [fnum(mean_sd_se([arm_sig(scans, labels, a, k, "legacy") for k in ks])[0])
                 for a in arms]
        row = [str(run), _run_type(run), str(len(ks))] + cells
        if two:
            ar, ac = arms[0], arms[2]
            m, _sd, se, _n = mean_sd_se([arm_sig(scans, labels, ac, k, "legacy")
                                         - arm_sig(scans, labels, ar, k, "legacy") for k in ks])
            row.append(dse(m, se))
        L.append("| " + " | ".join(row) + " |")
    L.append("")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scans", nargs="+", required=True, metavar="LABEL=DIR",
                    help="two or more driver-layout scan directories; the first is the "
                         "reference (arm A), the second the candidate (arms B, C)")
    ap.add_argument("--out", required=True, metavar="PREFIX",
                    help="writes PREFIX.md, PREFIX.tsv and PREFIX_splits.tsv")
    ap.add_argument("--exclude-runs", nargs="*", type=int, default=[24578],
                    help="runs to skip (default: 24578)")
    args = ap.parse_args()

    labels, dirs = [], {}
    for s in args.scans:
        if "=" not in s:
            ap.error(f"--scans entry {s!r} is not LABEL=DIR")
        lab, d = s.split("=", 1)
        if lab in dirs:
            ap.error(f"duplicate label {lab!r}")
        labels.append(lab)
        dirs[lab] = d
    scans = {lab: load_scan(Path(dirs[lab]), set(args.exclude_runs)) for lab in labels}
    print("splits loaded:", {k: len(v) for k, v in scans.items()})
    if any(len(v) == 0 for v in scans.values()):
        sys.exit("ERROR: a scan directory has no usable result.tsv")
    keys = sorted(set.intersection(*[set(s) for s in scans.values()]))
    skipped = sorted(set.union(*[set(s) for s in scans.values()]) - set(keys))
    if skipped:
        print(f"WARNING: {len(skipped)} (run, split) not present in every scan", file=sys.stderr)

    arms = arm_defs(labels)
    rows, armstat = type_table(scans, labels, arms, keys)
    dec = decide(rows, arms)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    write_tsv(Path(f"{out}.tsv"), TSV_COLS, rows)
    srows, scols = split_rows(scans, labels, keys)
    write_tsv(Path(f"{out}_splits.tsv"), scols, srows)
    md = report(scans, labels, dirs, arms, keys, rows, armstat, dec, skipped, args)
    Path(f"{out}.md").write_text(md)
    print(md)
    print(f"wrote {out}.md, {out}.tsv, {out}_splits.tsv")


if __name__ == "__main__":
    main()
