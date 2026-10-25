#!/usr/bin/env python3
"""regauss_tables.py - old (fixed +-150*kS MeV) vs new (iterative +-2 sigma) Gaussian from a scan.

    python3 regauss_tables.py --scan-dir output/scan_20261024_spline --out output/regauss_20261025/spline

For every (run, split) of a driver-layout scan directory (<dir>/scan_<run:06d>/split_<n:05d>/
result.tsv, read through compare_scan_dirs.load_scan) the Gaussian of the four inner-ring
spectra (legacy / gated x pre / post recalibration) is refit from the histogram stored in
result.tsv ("inner-ring histogram bin centres + counts") with the benchmark's
fit_gauss_reported in both definitions:
  fixed150 = the old window, +-150*kS MeV around the histogram argmax (reproduces the stored
             Gaussian to the printed precision);
  2sigma   = iterative refit in [mu-2 sigma, mu+2 sigma] until mu and sigma change by < 1 %
             (max 5 refits), started from the fixed fit; the --gauss-window default.
Per run type the mean over splits is tabulated: sigma/E (%), mu (MeV), the window half-width
(fixed150: the nominal 150*kS; 2sigma: mean of the final (hi-lo)/2) and the mean number of
refits.  The run type's preferred fit (gate_scan_20261022.R_PREF: Gaussian for 2.24 and
0.73 GeV H2, CB elsewhere) is flagged; the stored CB sigma/E is listed for reference (not
refitted).  n_fb = splits whose 2sigma refit did not run (fixed-window value kept / NaN); n_cap = splits
that did not converge to 1 % within 5 refits (the fit's converged flag).

Outputs: OUT.md, OUT.tsv (one row per run type x path x pre/post).
"""
from __future__ import annotations

import argparse
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import compare_scan_dirs as C                      # noqa: E402
import gate_scan_20261022 as G                     # noqa: E402

PATHS = (("legacy", "legacy"), ("new", "gated"))
WHENS = ("post", "pre")
COLS = ["type", "slug", "pref_fit", "path", "when", "n_splits", "sigmaE_fixed150_pct",
        "sigmaE_2sigma_pct", "delta_pp", "mu_fixed150", "mu_2sigma", "hw_fixed150_MeV",
        "hw_2sigma_MeV", "n_iter_mean", "n_fallback", "n_at_cap", "sigmaE_CB_stored_pct", "pref_is_gauss"]


def mean(x):
    a = np.asarray([v for v in x if np.isfinite(v)], dtype=float)
    return float(a.mean()) if len(a) else float("nan")


def collect(scan, exclude):
    """{(run, split): {'fixed': fits, 'new': fits, 'kS': kS, 'rec': rec}} (records without a
    stored histogram are skipped and counted)."""
    recs = C.load_scan(Path(scan), exclude)
    out, miss = {}, 0
    for key, rec in recs.items():
        fx = C.refit_gauss(rec, "fixed150")
        if fx is None:
            miss += 1
            continue
        nw = C.refit_gauss(rec, "2sigma")
        hdr = "\n".join(l for l in Path(rec["path"]).read_text(errors="replace").splitlines()
                        if l.startswith("#"))
        kS = float(C._RX_KS.search(hdr).group(2))
        out[key] = dict(fixed=fx, new=nw, kS=kS, rec=rec)
    if miss:
        print(f"WARNING: {miss} split(s) without stored histogram skipped", file=sys.stderr)
    return out


def rows_of(data):
    rows = []
    for rt in C.TYPE_ORDER:
        ks = [k for k in sorted(data) if C._run_type(k[0]) == rt]
        if not ks:
            continue
        pref = "CB" if G.R_PREF[rt] == "cb" else "Gauss"
        for p, pname in PATHS:
            for w in WHENS:
                def per(which, i):
                    return [data[k][which][(p, w)][i] for k in ks]
                s_fx = [100 * a / b if b > 0 else float("nan")
                        for a, b in zip(per("fixed", 1), per("fixed", 0))]
                s_nw = [100 * a / b if b > 0 else float("nan")
                        for a, b in zip(per("new", 1), per("new", 0))]
                hw = [0.5 * (data[k]["new"][(p, w)][5] - data[k]["new"][(p, w)][4]) for k in ks]
                nit = per("new", 6)
                cb = [100 * G.sigma(data[k]["rec"], p, f"CB_{w}") for k in ks]
                rows.append(dict(
                    type=rt, slug=G.runtype_slug(ks[0][0]), pref_fit=pref, path=pname, when=w,
                    n_splits=len(ks), sigmaE_fixed150_pct=mean(s_fx), sigmaE_2sigma_pct=mean(s_nw),
                    delta_pp=mean([b - a for a, b in zip(s_fx, s_nw)]),
                    mu_fixed150=mean(per("fixed", 0)), mu_2sigma=mean(per("new", 0)),
                    hw_fixed150_MeV=mean([150.0 * data[k]["kS"] for k in ks]),
                    hw_2sigma_MeV=mean(hw), n_iter_mean=mean(nit),
                    n_fallback=int(sum(1 for v in nit if v == 0)),
                    n_at_cap=int(sum(1 for k in ks if not data[k]["new"][(p, w)][7])),
                    sigmaE_CB_stored_pct=mean(cb), pref_is_gauss=int(pref == "Gauss")))
    return rows


def fmt(x, nd=2):
    return "NA" if x is None or not np.isfinite(x) else f"{x:.{nd}f}"


def markdown(rows, scan):
    L = [f"# Gaussian definition: fixed +-150*kS MeV vs iterative +-2 sigma", "",
         f"Scan `{scan}`; histograms from result.tsv, fits by benchmark fit_gauss_reported. "
         "Mean over splits. sigma/E in %, delta = 2sigma - fixed150 in pp; window half-width "
         "in MeV (fixed150: nominal 150*kS; 2sigma: mean final (hi-lo)/2); iter = mean number "
         "of refits; n_fb = splits without a 2sigma refit; n_cap = splits that did not converge "
         "within 5 refits (converged flag). CB = stored CB sigma/E (unchanged). "
         "**Bold** = the run type's preferred fit is the Gaussian.", ""]
    for rt in C.TYPE_ORDER:
        rr = [r for r in rows if r["type"] == rt]
        if not rr:
            continue
        pref = rr[0]["pref_fit"]
        L += [f"## {rt}  (preferred fit: {pref}, {rr[0]['n_splits']} splits)", "",
              "| path | spectrum | sigma/E fixed150 | sigma/E 2sigma | delta pp | mu fixed150 | "
              "mu 2sigma | hw fixed150 | hw 2sigma | iter | n_fb | n_cap | CB sigma/E |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
        for r in rr:
            b = (lambda s: f"**{s}**") if r["pref_is_gauss"] and r["when"] == "post" else (lambda s: s)
            L.append("| " + " | ".join([
                r["path"], r["when"] + "-cal", b(fmt(r["sigmaE_fixed150_pct"])),
                b(fmt(r["sigmaE_2sigma_pct"])), fmt(r["delta_pp"], 3), fmt(r["mu_fixed150"], 1),
                fmt(r["mu_2sigma"], 1), fmt(r["hw_fixed150_MeV"], 0), fmt(r["hw_2sigma_MeV"], 0),
                fmt(r["n_iter_mean"], 1), str(r["n_fallback"]), str(r["n_at_cap"]), fmt(r["sigmaE_CB_stored_pct"])]) + " |")
        L.append("")
    return "\n".join(L) + "\n"


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0],
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--scan-dir", required=True, help="driver-layout scan directory")
    ap.add_argument("--out", required=True, metavar="PREFIX", help="writes PREFIX.md and PREFIX.tsv")
    ap.add_argument("--exclude-runs", nargs="*", type=int, default=[24578],
                    help="runs to skip (default: 24578)")
    args = ap.parse_args()
    data = collect(args.scan_dir, set(args.exclude_runs))
    if not data:
        sys.exit("ERROR: no usable result.tsv with a stored histogram")
    rows = rows_of(data)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)
    C.write_tsv(Path(f"{out}.tsv"), COLS, rows, nd=5)
    md = markdown(rows, args.scan_dir)
    Path(f"{out}.md").write_text(md)
    print(md)
    print(f"wrote {out}.md, {out}.tsv")


if __name__ == "__main__":
    main()
