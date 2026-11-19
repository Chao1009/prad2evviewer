#!/usr/bin/env python3
"""exotic_summary.py -- P6c of the exotic-waveform study: one PDF with the latest numbers and every plot of P1-P10 (P8, P10, P9 appended in the analysis order).

    python3 exotic_summary.py [--out exotic_summary_20261118.pdf]

Front matter (new pages, numbers read from the current TSVs, never retyped): definitions, rates, saturation, moved matrix, top modules, radial
structure, amplitudes, f2 diagnosis, decisions, P7 cluster level, saturated seeds in the elastic window (read from exotic_P7_saturated_seeds_*.tsv), P8 energy error, P10 yield impact, P9 calibration impact and the gate (read from the P8, P10, P9 TSVs), P12 rates vs energy slice (P12 TSV), P13 drop policy (P13 TSV).
Every percent / rate / fraction carries its statistical uncertainty (Wilson 68 %, z = 1, `x % [lo, hi]`, or +-) or a stated count (standing user rule).  Then, for each step P1..P7, P8, P10, P9, P12, P13, a section page and every page of that step's PDF (merged unchanged;
P1-P3 are NOT regenerated).  Reads only files in this directory plus the HyCal module map; writes only the output PDF.  No rescans, no git.
"""
from __future__ import annotations

import argparse
import io
import json
import os
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                                         # noqa: E402
from matplotlib.backends.backend_pdf import PdfPages                    # noqa: E402
from pypdf import PdfReader, PdfWriter                                  # noqa: E402

HERE = Path(__file__).resolve().parent
RTS = ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]
LAB = {"0p73GeV_H2": "0.73 GeV H2", "2p24GeV_H2": "2.24 GeV H2", "3p49GeV_H2": "3.49 GeV H2", "2p24GeV_Ta_X17": "2.24 GeV Ta X17"}
CATS = ["f1", "f2", "f3", "f4", "f5"]
CNAME = {"f1": "f1 tail-baseline", "f2": "f2 wide", "f3": "f3 smooth-wide", "f4": "f4 clipped+sec.", "f5": "f5 clipped alone"}
RADBINS = (0, 120, 170, 230, 320, 450, 700)

P4S = "exotic_P4_module_rates_20261117_summary.tsv"
P4M = "exotic_P4_module_rates_20261117.tsv"
P4T = "exotic_P4_module_rates_20261117_thresholds.tsv"
P5T = "exotic_P5_max_amplitude_20261118.tsv"
P6T = "exotic_P6_amp_vs_w50_20261118.tsv"
P7T = "exotic_P7_cluster_level_20261118.tsv"
P7S = "exotic_P7_saturated_seeds_20261118.tsv"

STEPS = [  # (tag, pdf, date, script, outcome)
    ("P1", "exotic_P1_quantities_20261117.pdf", "2026-11-17", "exotic_p1_quantities.py",
     "Quantities Q1-Q14 on the 12 user-identified examples (types 1-4) plus random samples; first cuts f1 db > 8 sigma and tail < -5 sigma, f2 w50 ratio > 3, f4 clipped. "
     "Made with the pre-exclusive definitions; includes run 24937."),
    ("P2", "exotic_P2_flags_20261117.pdf", "2026-11-17", "exotic_p2_flags.py",
     "Flag definitions tested on the examples; flags behind PILEUP_EXOTIC=1 (BASE_SLOPE and CLIP_SINGLE off); f3 stretched-template flag catches types 1 and 2, not type 3. "
     "Pre-exclusive definitions; includes run 24937."),
    ("P3", "exotic_P3_sample_20261117.pdf", "2026-11-17", "exotic_p3_sample.py",
     "2-file sample: rates, overlaps, 20 examples per flag, BASE_SLOPE/CLIP_SINGLE panels; any-flag 0.08-0.15 % in H2, 1.08 % in Ta X17. "
     "Pre-exclusive definitions; includes run 24937."),
    ("P4", "exotic_P4_module_rates_20261117.pdf", "2026-11-17 (regenerated in place 2026-11-18)", "exotic_p4_module_rates.py",
     "Full scan (12 runs, 99 files after dropping 24937 from 103): exclusive categories f4 > f5 > f1 > f2 > f3, per-module maps and top lists, saturation per module and radius, thresholds; run 24937 excluded."),
    ("P5", "exotic_P5_max_amplitude_20261118.pdf", "2026-11-18", "exotic_p5_max_amplitude.py",
     "Maximum amplitude per category and run type (raw peak and fitted A, median and 16/84 %): f4/f5 saturated at about 3930 ADC, f1 about 60-80 ADC, f2 about 10 raw ADC in Ta; "
     "H2 from the 2-file sample, Ta from the full scan."),
    ("P6", "exotic_P6_amp_vs_w50_20261118.pdf", "2026-11-18", "exotic_p6_amp_vs_w50.py",
     "Amplitude vs w50, saturation with E > 0.5 E_beam ratio, f2 high-amplitude diagnosis, 40 + 80 examples; the high-A f2 mode is merged k >= 4 pile-up (2.24 H2) or W335 (0.73 H2)."),
    ("P7", "exotic_P7_cluster_level_20261118.pdf", "2026-11-18", "exotic_p7_cluster_level.py",
     "Flags joined to the benchmark's own reconstructed clusters (control v3): fraction of clusters with a flagged module by role, category, ring and energy bin, for all clusters, the ep elastic window and Moller pairs; "
     "flagged energy fraction; fake-cluster candidates (flagged seed); no H2 calibration cell reaches 0.1 % of the cluster energy, Ta elastic f1-neighbour and f3-seed cells do. "
     "Saturated-seed table: pages 10-11 of this summary (exotic_p7_saturated_seeds.py)."),
    ("P8", "exotic_P8_energy_error_20261118.pdf", "2026-11-18", "exotic_p8_energy_error.py",
     "Energy error per category on 391,292 refitted flagged cluster members: scale shift Ta elastic f1 sloped baseline -0.0144 +- 0.0003 %, f3 stretched -0.0028 %, Moller f1 -0.0051 %, Ta inner-ring modules -0.018 to -0.027 %; "
     "every credible H2 cell below 0.01 %; credible alternatives per category (f1 sloped baseline, f3 stretched in Ta only, f5 k = 1 for mild clipping, default for f2, f4 and H2 f3)."),
    ("P10", "exotic_P10_yield_impact_20261118.pdf", "2026-11-18", "exotic_p10_yield_impact.py",
     "Veto losses V1-V4 per run type and ring, Ta Moller pairs, contamination removed by V3 (fakes), Ta inner-ring f1 rate vs occupancy and current, out-of-cluster Ta f1 rows. "
     "V3 costs 0 / 0.0044 / 0.0067 / 0 % of the elastic window (0.73 / 2.24 / 3.49 H2 / Ta); V4 0.286 % on Ta; V2 not viable."),
    ("P9", "exotic_P9_calibration_impact_20261118.pdf", "2026-11-18", "exotic_p9_calibration_impact.py",
     "Benchmark refit three ways (V3 / V4 veto, C corrected) against the default: no arm moves oos sigma/E, the peak or the tails beyond 2 sigma of the paired error, no free module gain moves by 2 sigma_g. "
     "Gate: V3 on 2.24 / 3.49 GeV H2, no V4, no correction; 0.73 GeV and Ta X17 neither."),
    ("P12", "exotic_P12_rates_vs_energy_20261118.pdf", "2026-11-18", "exotic_p12_rates_vs_energy.py",
     "Exotic-category rates per energy slice (10-50 ... > 4000 MeV), per run type and flag (f1-f5, any), per channel (kept PbWO4 channels, A_main x gain) and per cluster (flagged seed / >= 1 flagged member), "
     "Wilson 68 % intervals; impact per slice = flagged energy fraction removed by dropping flagged channels.  Tables on page 15 of this summary."),
    ("P13", "exotic_P13_drop_flagged_20261118.pdf", "2026-11-18", "exotic_p13_drop_flagged.py (+ exotic_p13_recluster.py)",
     "P7-P10 recomputed under channel-level drop of f1-f5 (arm D; D' = drop in the gain solve only), against the default control v3 (arm A): outcome text read from the P13 TSV in dyn_outcomes."),
]


# ----------------------------------------------------------------------------------------------------------------------------------------------
def read(name):
    return pd.read_csv(HERE / name, sep="\t")


def dyn_outcomes():
    """Section-page outcome texts whose numbers come from TSVs (P3, P8, P10): every percentage with its interval / error."""
    out = {}
    p3 = read("exotic_P3_sample_20261117.tsv")
    r = p3[(p3.section == "rate") & (p3["item"] == "any")].set_index("run_type")
    h = [f"{g3(r.loc[rt].percent)} % [{g3(r.loc[rt].percent_lo68)}, {g3(r.loc[rt].percent_hi68)}]" for rt in RTS[:3]]
    ta = r.loc[RTS[3]]
    out["P3"] = ("2-file sample: rates, overlaps, 20 examples per flag, BASE_SLOPE/CLIP_SINGLE panels; any-flag (68 % Wilson) " + ", ".join(h) + f" in H2 (0.73 / 2.24 / 3.49 GeV), {ta.percent:.3f} % [{ta.percent_lo68:.3f}, {ta.percent_hi68:.3f}] in Ta X17. "
                 "Pre-exclusive definitions; includes run 24937.")
    d = read("exotic_P8_energy_error_20261118.tsv")
    sc = d[d.table == "scale"]
    mo = d[d.table == "module"]

    def S(rt, pop, cat, var, ring="all"):
        x = sc[(sc.run_type == rt) & (sc.population == pop) & (sc.category == cat) & (sc.variant == var) & (sc.ring == ring)]
        return f"{100 * x.S.iloc[0]:+.4f} +- {100 * x.S_err.iloc[0]:.4f} %"
    ta_ = "2p24GeV_Ta_X17"
    inner = [mo[(mo.run_type == ta_) & (mo.module == m) & (mo.category == "any")].iloc[0] for m in ["W493", "W494", "W559", "W564", "W593", "W663"]]
    lo_, hi_ = min(inner, key=lambda r: r.S_el_seeded), max(inner, key=lambda r: r.S_el_seeded)
    n_join = int(d[(d.table == "diag") & (d.run_type == "join")].n.iloc[0])
    out["P8"] = (f"Energy error per category on {n_join:,} refitted flagged cluster members: scale shift Ta elastic f1 sloped baseline {S(ta_, 'elastic', 'f1 tail-baseline', 'bs')}, f3 stretched {S(ta_, 'elastic', 'f3 smooth-wide', 'st')}, "
                 f"Moller f1 {S(ta_, 'moller', 'f1 tail-baseline', 'bs')}, Ta inner-ring modules {100 * hi_.S_el_seeded:+.3f} +- {100 * hi_.S_el_err:.3f} to {100 * lo_.S_el_seeded:+.3f} +- {100 * lo_.S_el_err:.3f} %; "
                 "every credible H2 cell below 0.01 % (page 12 of this summary); credible alternatives per category (f1 sloped baseline, f3 stretched in Ta only, f5 k = 1 for mild clipping, default for f2, f4 and H2 f3).")
    p = read("exotic_P10_yield_impact_20261118.tsv")
    v = p[p.table == "veto"]

    def V(ds, pop, ve):
        x = v[(v.dataset == ds) & (v.population == pop) & (v.veto == ve) & (v.axis == "overall") & (v["bin"] == "all")].iloc[0]
        return f"{100 * x.frac:.4f} % [{100 * x.lo:.4f}, {100 * x.hi:.4f}]"
    out["P10"] = ("Veto losses V1-V4 per run type and ring, Ta Moller pairs, contamination removed by V3 (fakes), Ta inner-ring f1 rate vs occupancy and current, out-of-cluster Ta f1 rows. V3 costs "
                  + "; ".join(f"{V(ds, 'elastic', 'V3')} ({lab})" for ds, lab in zip(RTS, ["0.73 H2", "2.24 H2", "3.49 H2", "Ta"])) + f" of the elastic window (Wilson 68 %); V4 {V('2p24GeV_Ta_X17', 'elastic', 'V4')} on Ta; V2 not viable.")
    try:
        R = p13_numbers(p13_load())
        h1, h2, t = R["2p24GeV_H2"], R["3p49GeV_H2"], R["2p24GeV_Ta_X17"]
        out["P13"] = ("P7-P10 recomputed under channel-level drop of f1-f5 (arm D; D' = drop in the gain solve only) against the default control v3 (arm A).  Gate: drop policy NOT adopted; keep V3 on 2.24 / 3.49 GeV H2, no correction.  "
                      f"Elastic-window loss (Wilson 68 %) {ci1(h1['gross'][0], h1['gross'][1], 4)} / {ci1(h2['gross'][0], h2['gross'][1], 4)} at 2.24 / 3.49 GeV H2 (V3: {_v3('2p24GeV_H2')} / {_v3('3p49GeV_H2')}); "
                      f"Ta X17 low tail {t['tail']['diff']:+.2f} +- {t['tail']['se']:.2f} pp, window loss net {t['net'][0]:.2f} +- {t['net'][1]:.2f} %.  Summary page: page 16.")
    except Exception as ex:
        print("P13 outcome text not built:", ex)
    return out


def pct(x, nd=4):
    return f"{x:.{nd}f}" if x >= 1e-3 else f"{x:.2e}"


def wil(k, n):
    """Wilson 68 % interval (z = 1) of k / n in PERCENT (copy of exotic_p4_module_rates.wilson)."""
    k, n = float(k), max(float(n), 1.0)
    p = k / n
    den = 1.0 + 1.0 / n
    mid = (p + 1.0 / (2 * n)) / den
    h = np.sqrt(p * (1 - p) / n + 1.0 / (4 * n * n)) / den
    lo, hi = (0.0 if k == 0 else max(mid - h, 0.0)), (1.0 if k >= n else min(mid + h, 1.0))
    lo, hi = (0.0 if abs(lo) < 1e-12 else lo), hi
    return 100 * lo, 100 * hi


def g3(x):
    return f"{x:.3g}"


def ci1(k, n, nd=3):
    """'x % [lo, hi]' of k / n (Wilson 68 %)."""
    lo, hi = wil(k, n)
    f = (lambda v: f"{v:.{nd}f}") if nd < 3 else g3             # nd < 3: decimals; nd >= 3: three significant digits
    return f"{f(100 * k / max(n, 1))} % [{f(lo)}, {f(hi)}]"


def cik(k, n, nd=3):
    return ci1(k, n, nd) + f" ({int(k)} of {int(n):,})"


class Pages:
    def __init__(self, pdf):
        self.pdf = pdf

    def text(self, title, lines, size=8.2, sub=None):
        fig = plt.figure(figsize=(16, 11))
        fig.text(0.04, 0.965, title, fontsize=17, fontweight="bold", va="top")
        y = 0.925
        if sub:
            fig.text(0.04, y, sub, fontsize=9.5, va="top", style="italic")
            y -= 0.03
        import textwrap
        lim = int(1060 / (0.602 * size))
        lines = [w for ln in lines for w in (textwrap.wrap(ln, lim, subsequent_indent="    ") if len(ln) > lim else [ln])]
        fig.text(0.04, y, "\n".join(lines), fontsize=size, va="top", family="monospace", linespacing=1.42)
        self.pdf.savefig(fig)
        plt.close(fig)


def wrap(s, w=190, ind="  "):
    import textwrap
    return textwrap.wrap(s, w, initial_indent=ind, subsequent_indent=ind + "  ")


# ----------------------------------------------------------------------------------------------------------------------------------------------
def build_numbers():
    S = read(P4S)
    M = read(P4M)
    T = read(P4T)
    P5 = read(P5T)
    P6 = read(P6T)
    N = {}
    info = {"checks": []}
    for rt in RTS:
        N[rt] = int(S[(S.run_type == rt) & (S.threshold_set == "uniform") & (S.flag == "f1")].n_kept.iloc[0])

    def sv(rt, ts, flag, col="n_flag"):
        return S[(S.run_type == rt) & (S.threshold_set == ts) & (S.flag == flag)][col].iloc[0]

    # internal consistency: raw - moved out = exclusive
    for rt in RTS:
        for c in CATS:
            moved_out = sum(sv(rt, "uniform", f"moved_{c}_to_{d}") for d in CATS if d != c)
            if int(sv(rt, "uniform", c) - moved_out) != int(sv(rt, "uniform", f"excl_{c}")):
                info["checks"].append(f"raw - moved != exclusive for {rt} {c}")
        if int(sum(sv(rt, "uniform", f"excl_{c}") for c in CATS)) != int(sv(rt, "uniform", "excl_any")):
            info["checks"].append(f"sum of exclusive != excl_any for {rt}")
    return S, M, T, P5, P6, N, sv, info


def page_definitions(P):
    L = [
        "Study: exotic waveforms in HyCal, PbWO4 modules only, run 24937 (3.49 GeV H2, 'not a good run') excluded from P4-P6.  Flags are markers (BASE_SLOPE and CLIP_SINGLE off, PILEUP_EXOTIC=1 only).",
        "Template FWHM 23.5 ns (PbWO4, all run types), so w50 > 3 means data FWHM > 70 ns.  Uniform P2 thresholds are the reference; per-type thresholds are in the appendix table (page 3).",
        "",
        "The four user-identified types and their flags",
        "  type 1  tail below the baseline after the pulse   f1  db > 8 sigma  AND  tail < -5 sigma  (db = baseline step, tail = post-pulse level, in sigma of the k = 1 core)",
        "  type 2  anomalously wide smooth pulse             f2  w50 > 3 (data FWHM > 70 ns), NOT clipped (nclip = 0), NOT window-truncated (t0 > 284 ns excluded)",
        "  type 3  smooth-wide, mimics close pile-up         f3  k >= 2  AND  the stretched single template is preferred by BIC (chi2_stretch <= chi2_k + LAM (k-1) ln nm, LAM = 2)",
        "  type 4  clipped / saturated pulse                 f4  clipped AND (secondary pulse within 8 ns of the clip  OR  E_sec > 0.5 E_beam)",
        "          (added 2026-11-18)                        f5  clipped AND NOT f4  (a lone saturated pulse; kept separately to quantify saturation)",
        "",
        "Exclusive categories (user decision 2026-11-18): each channel counts once, priority  f4 > f5 > f1 > f2 > f3.  'raw' below = the flag before the priority is applied.",
        "  f2 additionally excludes clipped records and window-truncated pulses (page 3 gives the removed fractions and counts per run type, from the summary TSV check rows).",
        "",
        "Thresholds and exclusions at a glance",
        "  uniform (P2):  f1 db > 8, tail < -5;  f2 w50 > 3;  f3 LAM = 2;  f4 8 ns / 0.5 E_beam.   per-type = median + z_ref sigma of the k = 1 core of each run type (page 3).",
        "  per-module rates are quoted only where n_kept >= 100 and n_flag >= 5 (stat_ok in the module TSV); maps colour modules with n_kept >= 30.",
        "  every percentage / rate / fraction below carries its statistical uncertainty: Wilson 68 % (z = 1) written `x % [lo, hi]` (summary TSV columns lo68_percent, hi68_percent), or +-, or a count `k of N`.  Thresholds such as 0.5 E_beam, 0.1 % criteria and E/E_beam are definitions, not rates.",
        "  out of scope here: other items (Moller pairs appear in P7-P10); BASE_SLOPE and CLIP_SINGLE (off); Ta X17 per-type 'f1 db > 9.76' is a reference only.",
        "",
        "Page map: 1 definitions | 2 exclusive and raw rates, saturation, moved matrix | 3 thresholds, per-type rates | 4 top modules per category | 5 radial structure |",
        "          6 amplitudes (P5) | 7 f2 diagnosis (P6) | 8 decisions, open items and the gate | 9 P7 cluster level | 10-11 saturated seeds in the elastic window (P7) |",
        "          12 P8 energy error | 13 P10 yield impact | 14 P9 calibration impact and gate decision | 15 P12 rates vs energy slice | 16 P13 drop policy (P7-P10 recomputed) |",
        "          then P1, P2, P3, P4, P5, P6, P7, P8, P10, P9, P12, P13 (section page + every page of the step PDF; P10 precedes P9, the analysis order; P12, P13 last).",
        "",
        "P1-P3 pages below were made under the pre-exclusive definitions (f1-f4 overlapping, first cut values) and include run 24937.  They were NOT regenerated.",
    ]
    P.text("Exotic waveforms in HyCal: summary of P1-P13 (2026-11-18)", L, size=9.2,
           sub="Definitions as of 2026-11-18; numbers on the following pages are read from the current TSVs.")


def page_rates(P, S, N, sv):
    L = []
    cw = 48

    def cell(rt, key):
        n = int(sv(rt, "uniform", key))
        return f"{n:>8d} = {pct(sv(rt, 'uniform', key, 'rate_percent'), 5)} % [{pct(sv(rt, 'uniform', key, 'lo68_percent'), 5)}, {pct(sv(rt, 'uniform', key, 'hi68_percent'), 5)}]".rjust(cw)
    L.append("A. EXCLUSIVE category rates, uniform thresholds, full scan (P4).  k = n flagged channels, rate = k / N kept PbWO4 channels, % [Wilson 68 %]")
    L.append(f"{'category':<10}" + "".join(f"{LAB[rt] + f' (N = {N[rt]:,})':>{cw}}" for rt in RTS))
    for c in CATS + ["any"]:
        L.append(f"{c:<10}" + "".join(cell(rt, f"excl_{c}") for rt in RTS))
    L.append("")
    L.append("B. RAW (pre-priority) flag rates, same thresholds; k, % [Wilson 68 %]")
    L.append(f"{'category':<10}" + "".join(f"{LAB[rt] + f' (N = {N[rt]:,})':>{cw}}" for rt in RTS))
    for c in CATS + ["any"]:
        L.append(f"{c if c != 'any' else 'any (union)':<10}" + "".join(cell(rt, c) for rt in RTS))
    L.append("")
    L.append("C. Moved matrix (uniform): channels flagged raw in category a that the priority assigns to category b: k = n moved, % [Wilson 68 %] of the k_raw channels raw-flagged in a (summary TSV of_raw_percent, of_raw_lo68_percent, of_raw_hi68_percent)")
    for rt in RTS:
        mv = [(a, b, int(sv(rt, "uniform", f"moved_{a}_to_{b}"))) for a in CATS for b in CATS if a != b]
        nz = [f"{a}->{b}: {n} of {int(sv(rt, 'uniform', a))} = {pct(sv(rt, 'uniform', f'moved_{a}_to_{b}', 'of_raw_percent'), 3)} % [{pct(sv(rt, 'uniform', f'moved_{a}_to_{b}', 'of_raw_lo68_percent'), 3)}, {pct(sv(rt, 'uniform', f'moved_{a}_to_{b}', 'of_raw_hi68_percent'), 3)}]" for a, b, n in mv if n]
        zero = [sv(rt, "uniform", f"moved_{a}_to_{b}", "of_raw_hi68_percent") for a, b, n in mv if not n and int(sv(rt, "uniform", a)) > 0]
        for i in range(0, max(len(nz), 1), 3):
            L.append(f"  {LAB[rt] if i == 0 else '':<17}" + ("   |   ".join(nz[i:i + 3]) if nz else "none"))
        L.append(f"  {'':<17}(all other entries are 0 moved; their upper 68 % bound of the raw share is at most {pct(max(zero), 3)} %)" if zero else "")
    chk = {rt: {k: int(sv(rt, "check", k)) for k in ("uni_f1_and_f4", "uni_f2_and_f4", "uni_f3_and_f4", "uni_f5_and_f3")} for rt in RTS}
    L.append("   raw overlaps (check rows, counts): f1&f4 = " + "/".join(str(chk[rt]["uni_f1_and_f4"]) for rt in RTS) + ", f2&f4 = " + "/".join(str(chk[rt]["uni_f2_and_f4"]) for rt in RTS)
             + ", f3&f4 = " + "/".join(str(chk[rt]["uni_f3_and_f4"]) for rt in RTS) + ", f3&f5 = " + "/".join(str(chk[rt]["uni_f5_and_f3"]) for rt in RTS) + "  (0.73 / 2.24 / 3.49 H2 / Ta X17)")
    L.append("")
    L.append("D. SATURATION per run type: clipped records (f4 + f5), full scan P4 'clipped' rows, % of kept [Wilson 68 %]; P6 sample/full 'saturation' rows with the E > 0.5 E_beam ratio, k of N, % [68 %]")
    L.append(f"{'run type':<17}{'N kept':>12}{'n clip':>8}{'f4':>6}{'f5':>6}{'% of kept [68 %]':>34}   | P6 source, E>0.5 E_b pulses: n_clip / n  = % [68 %]")
    return L

def page_sat_rows(L, S, N, sv, P6):
    sat = P6[P6.section == "saturation"]
    for rt in RTS:
        nc = int(sv(rt, "uniform", "clipped"))
        lo, hi = sv(rt, "uniform", "clipped", "lo68_percent"), sv(rt, "uniform", "clipped", "hi68_percent")
        r0 = f"{LAB[rt]:<17}{N[rt]:>12d}{nc:>8d}{int(sv(rt, 'uniform', 'f4')):>6d}{int(sv(rt, 'uniform', 'f5')):>6d}{pct(sv(rt, 'uniform', 'clipped', 'rate_percent'), 5):>12} % [{pct(lo, 5)}, {pct(hi, 5)}]   | "
        rows = sat[sat.run_type == rt]
        first = True
        for _, r in rows.iterrows():
            src = "sample (2 files/run)" if str(r.source).startswith("sample") else "19 Ta npz of the full scan"
            s = f"{src}: {int(r.n_clip_E_gt_0p5)} / {int(r.n_E_gt_0p5)} = {pct(r.pct_clip_E_gt_0p5, 4)} % [{pct(r.lo68_E, 4)}, {pct(r.hi68_E, 4)}]  (n_kept {int(r.n_kept)}, clipped {pct(r.pct_clip, 5)} % [{pct(r.lo68, 5)}, {pct(r.hi68, 5)}])"
            L.append((r0 if first else " " * len(r0)) + s)
            first = False
    last = {rt: sat[(sat.run_type == rt)].iloc[-1] for rt in RTS}
    ci = lambda rt: f"{pct(last[rt].pct_clip_E_gt_0p5, 4)} % [{pct(last[rt].lo68_E, 4)}, {pct(last[rt].hi68_E, 4)}] ({int(last[rt].n_clip_E_gt_0p5)} of {int(last[rt].n_E_gt_0p5):,})"
    L.append(f"   Saturated share of E > 0.5 E_beam pulses (Wilson 68 %): {ci('2p24GeV_H2')} at 2.24 GeV H2, {ci('3p49GeV_H2')} at 3.49, {ci('0p73GeV_H2')} at 0.73, {ci('2p24GeV_Ta_X17')} in Ta (H2 from the 2-file sample).")
    L.append(f"   Concluded 2026-11-18: saturation up to 0.1 % of full-energy pulses for now (the 0.1 % is the working bound; measured at 2.24 GeV H2 {pct(last['2p24GeV_H2'].pct_clip_E_gt_0p5, 3)} % [{pct(last['2p24GeV_H2'].lo68_E, 3)}, {pct(last['2p24GeV_H2'].hi68_E, 3)}]); sample statistics are too small, keep f5 for the ifarm full scan.")
    return L

def page_thresholds(P, S, T, N, sv):
    L = ["Thresholds (thresholds TSV) and per-type EXCLUSIVE counts (summary TSV, threshold_set = type), vs uniform in brackets (counts; the rates with Wilson 68 % intervals are in the summary TSV and on page 2 for uniform)", ""]
    L.append(f"{'run type':<17}{'f1 db >':>9}{'f1 tail <':>11}{'f2 w50 >':>10}{'f3 LAM':>8}   | " + "".join(f"{c + ' type [uniform]':>29}" for c in CATS))
    for rt in RTS:
        t = T[T.run_type == rt].set_index("quantity")
        row = f"{LAB[rt]:<17}{t.loc['f1 db >', 'per_type']:>9.3f}{t.loc['f1 tail <', 'per_type']:>11.3f}{t.loc['f2 w50 >', 'per_type']:>10.3f}{t.loc['f3 LAM', 'per_type']:>8.3f}   | "
        for c in CATS:
            row += f"{int(sv(rt, 'type', 'excl_' + c)):>9d} [{int(sv(rt, 'uniform', 'excl_' + c)):>8d}]".rjust(29)
        L.append(row)
    L.append(f"{'uniform (P2)':<17}{8:>9.3f}{-5:>11.3f}{3:>10.3f}{2:>8.3f}")
    L.append("")
    L.append("k = 1 core of each run type (median / sigma, thresholds TSV): db med, db sigma, tail med, tail sigma, w50 med, w50 sigma, chi2/dof med")
    for rt in RTS:
        t = T[T.run_type == rt].set_index("quantity").per_type
        L.append(f"  {LAB[rt]:<17}" + "  ".join(f"{k.replace('core ', '')} {t[k]:.3f}" for k in ["core db median", "core db sigma", "core tail median", "core tail sigma", "core w50 median", "core w50 sigma", "core chi2 median"]))
    L.append("")
    L.append("f2 exclusions (summary TSV check rows, uniform): flags before the window-truncation and clip rules, then after")
    for rt in RTS:
        b, a, tr = int(sv(rt, "check", "uniform_before"), ), int(sv(rt, "check", "uniform_after")), int(sv(rt, "check", "uniform_trunc_only"))
        L.append(f"  {LAB[rt]:<17} before {b:>7d}  after truncation {tr:>7d}  after clip {a:>7d}   (truncation removes {b - tr} of {b} = {ci1(b - tr, b, 1)}, clip removes {tr - a} of {tr})")
    L.append("")
    L.append("Per-type thresholds move little: only Ta X17 f1 db 8 -> 9.76; the Ta excess persists under per-type thresholds.  0.73 GeV db threshold is capped at the uniform 8.")
    return L


def top_modules(M, rt, cat, k=3):
    x = M[(M.threshold_set == "uniform") & (M.basis == "excl") & (M.run_type == rt) & (M.flag == cat) & (M.stat_ok == 1)]
    x = x.sort_values("rate", ascending=False).head(k)
    return x


def page_top(P, M):
    L = ["Top modules per EXCLUSIVE category and run type (module TSV: threshold_set uniform, basis excl, stat_ok = 1, i.e. n_kept >= 100 and n_flag >= 5), top 3 by rate",
         "entries: module  n_flag / n_kept = rate % [Wilson 68 %]", ""]
    for c in CATS:
        L.append(f"{CNAME[c]}")
        for rt in RTS:
            x = top_modules(M, rt, c)
            s = "   ".join(f"{r.module} {int(r.n_flag)}/{int(r.n_kept)} = {pct(100 * r.rate, 3)} % [{pct(100 * r.rate_lo68, 3)}, {pct(100 * r.rate_hi68, 3)}]" for r in x.itertuples()) or "no module passes the statistics criterion"
            nok = int(((M.threshold_set == 'uniform') & (M.basis == 'excl') & (M.run_type == rt) & (M.flag == c) & (M.stat_ok == 1)).sum())
            L.append(f"   {LAB[rt]:<17}({nok:>4d} stat-ok modules)  {s}")
        L.append("")
    return L


def radial(M, rt, cat, basis="excl"):
    import json as _j
    mp = _j.load(open(Path(os.environ.get("PRAD2_DATABASE_DIR", HERE.parents[4] / "database")) / "hycal_map.json"))
    geo = {m["n"]: m["geo"] for m in mp}
    x = M[(M.threshold_set == "uniform") & (M.basis == basis) & (M.run_type == rt) & (M.flag == cat)].copy()
    x["r"] = [float(np.hypot(geo[n]["x"], geo[n]["y"])) if n in geo else np.nan for n in x.module]
    x = x[x.module.str.startswith("W") & x.r.notna()]
    out = []
    for lo, hi in zip(RADBINS[:-1], RADBINS[1:]):
        y = x[(x.r >= lo) & (x.r < hi)]
        out.append((int(y.n_kept.sum()), int(y.n_flag.sum())))
    return out


def page_radial(P, M):
    L = ["Radial structure: EXCLUSIVE category rate % [Wilson 68 %] in PbWO4 radius bins from the HyCal centre (module TSV summed over the modules of each bin, uniform thresholds, all modules); last column: ratio of the 320-450 to the 0-120 mm rate +/- propagated error", ""]
    L.append(f"{'category':<17}{'run type':<17}" + "".join(f"{f'{lo}-{hi} mm':>31}" for lo, hi in zip(RADBINS[:-1], RADBINS[1:])) + f"{'320-450 / 0-120':>20}")
    tr = {}
    for c in CATS:
        for rt in RTS:
            rr = radial(M, rt, c)
            rates = [100 * f / n if n else np.nan for n, f in rr]
            ivs = [wil(f, n) if n else (np.nan, np.nan) for n, f in rr]
            if rates[0] and rates[0] > 0 and rr[4][0]:
                ratio = rates[4] / rates[0]
                s0, s4 = (ivs[0][1] - ivs[0][0]) / 2, (ivs[4][1] - ivs[4][0]) / 2
                err = ratio * np.sqrt((s0 / rates[0]) ** 2 + ((s4 / rates[4]) ** 2 if rates[4] > 0 else 0.0))
                rtxt = f"{ratio:.3g} +/- {err:.2g}" if rates[4] > 0 else f"0 (< {ivs[4][1] / rates[0]:.2g})"
            else:
                ratio, rtxt = np.nan, "-"
            tr[(c, rt)] = (rates, ratio, rr)
            L.append(f"{CNAME[c]:<17}{LAB[rt]:<17}" + "".join(f"{(f'{g3(r_)} [{g3(iv[0])}, {g3(iv[1])}]' if n else '-'):>31}" for r_, iv, (n, f) in zip(rates, ivs, rr)) + f"{rtxt:>20}")
        L.append("")
    return L, tr

def page_amp(P, P5):
    L = ["Maximum amplitude per exclusive category (P5 TSV): median [16 %, 84 %], pedestal-subtracted ADC and fraction of the beam energy (E/E_beam).  raw = peak sample; fit = largest fitted A (= A_main)",
         "Source: H2 rows = 2-file sample npz (the H2 pickles carry no amplitudes); Ta X17 rows = full scan npz (19 files).  ref = 1 % random sample of kept channels.", ""]
    L.append(f"{'run type':<17}{'cat':<5}{'N':>8}{'raw ADC':>26}{'raw E/E_beam':>26}{'fit ADC':>26}{'fit E/E_beam':>26}  source")
    for rt in RTS:
        for c in CATS + ["ref"]:
            g = P5[(P5.rt == rt) & (P5.flag == c)].set_index("defn")
            if "raw_adc" not in g.index:
                continue
            n = int(g.loc["raw_adc", "N"])
            src = "19 Ta npz (of 99 files)" if str(g.loc["raw_adc", "source"]).startswith("full") else "sample (2 files/run)"
            if n == 0:
                L.append(f"{LAB[rt]:<17}{c:<5}{0:>8d}{'-':>26}{'-':>26}{'-':>26}{'-':>26}  {src}")
                continue

            def f(d, nd):
                r = g.loc[d]
                return f"{r['median']:.{nd}f} [{r.q16:.{nd}f}, {r.q84:.{nd}f}]".rjust(26)
            L.append(f"{LAB[rt]:<17}{c:<5}{n:>8d}{f('raw_adc', 1)}{f('raw_efrac', 3)}{f('fit_adc', 1)}{f('fit_efrac', 3)}  {src}")
        L.append("")
    g1 = P5[(P5.rt == "0p73GeV_H2") & (P5.flag == "f1")].set_index("defn")
    L.append("")
    L.append("Fraction of rows with raw peak >= 4095 (saturated), k of N, % [Wilson 68 %] (P5 TSV n_smax_ge_4095 / n_smax; f4, f5 and ref, plus any category with k > 0)")
    for rt in RTS:
        parts = []
        for c in CATS + ["ref"]:
            g = P5[(P5.rt == rt) & (P5.flag == c) & (P5.defn == "raw_adc")]
            if g.empty or int(g.n_smax.iloc[0]) == 0:
                continue
            r = g.iloc[0]
            if c in ("f4", "f5", "ref") or int(r.n_smax_ge_4095) > 0:
                parts.append(f"{c}: {int(r.n_smax_ge_4095)} of {int(r.n_smax):,} = {g3(100 * r.frac_smax_ge_4095)} % [{g3(100 * r.frac_smax_ge_4095_lo68)}, {g3(100 * r.frac_smax_ge_4095_hi68)}]")
        L.append(f"  {LAB[rt]:<17}" + "   ".join(parts))
    L.append("")
    L.append(f"Reading: f4/f5 raw = saturated (about 3930 ADC); f4 fitted 1-2 E_beam (14 at 0.73 GeV with N = 2).  f1 is a few percent of E_beam (median raw E/E_beam {100 * g1.loc['raw_efrac', 'median']:.0f} % [16-84 % quantile band {100 * g1.loc['raw_efrac', 'q16']:.0f}, {100 * g1.loc['raw_efrac', 'q84']:.0f}] at 0.73 GeV, N = {int(g1.loc['raw_efrac', 'N'])}; a quantile band of the amplitude distribution, not a statistical error).  f2 is bimodal at 0.73 and 2.24 GeV H2 (q84 about 1000 ADC, unclipped rows), 20-40 ADC at 3.49 GeV, about 10 raw / 31 fitted ADC in Ta; f3 is broad in H2 (q16 20-45, q84 450-560 ADC).")
    L.append("Small N: 0.73 GeV f4 (N = 2) and f5 (N = 0); a larger H2 N needs an amplitude-carrying rescan (pickles lack event ids).")
    return L


def page_f2(P, P6):
    g = P6[P6.section == "f2_high_group"]
    L = ["f2 diagnosis (P6 TSV, section f2_high_group): high-amplitude f2 = fitted A > 300 ADC.  'hi' = A > 300, 'lo' = A <= 300; k = number of fitted pulses; fractions in % [Wilson 68 %]", ""]
    L.append(f"{'run type':<17}{'source':<11}{'n_f2':>6}{'n hi':>6}{'hi frac % [68 %]':>26}{'A med(all f2)':>14}{'hi: k1/k2/k3/k>=4':>19}{'lo: k1/k2/k3/k>=4':>22}{'hi modules':>11}{'top1 share % [68 %]':>24}  w50 med hi/lo")
    pc = lambda x, lo, hi: f"{g3(100 * x)} [{g3(100 * lo)}, {g3(100 * hi)}]"
    for r in g.itertuples():
        src = "sample" if str(r.source).startswith("sample") else "full scan"
        L.append(f"{LAB[r.run_type]:<17}{src:<11}{int(r.n_f2):>6d}{int(r.n_f2_A_gt300):>6d}{pc(r.frac_A_gt300, r.frac_lo68, r.frac_hi68):>26}{r.A_median:>14.1f}"
                 f"{f'{int(r.k1_hi)}/{int(r.k2_hi)}/{int(r.k3_hi)}/{int(r.k4p_hi)}':>19}{f'{int(r.k1_lo)}/{int(r.k2_lo)}/{int(r.k3_lo)}/{int(r.k4p_lo)}':>22}"
                 f"{(int(r.n_modules_hi) if r.n_modules_hi == r.n_modules_hi else 0):>11d}{(pc(r.top1_share_hi, r.top1_share_hi_lo68, r.top1_share_hi_hi68) if r.top1_share_hi == r.top1_share_hi else '-'):>24}  {(f'{r.w50_median_hi:.2f}' if r.w50_median_hi == r.w50_median_hi else '-')}/{r.w50_median_lo:.2f}")
    L.append("")
    L.append("k >= 4 fraction within each group, k of N, % [Wilson 68 %] (hi: of the n hi channels, lo: of the n_f2 - n hi channels):")
    for r in g.itertuples():
        nh, nl = int(r.n_f2_A_gt300), int(r.n_f2 - r.n_f2_A_gt300)
        hi_t = f"hi {int(r.k4p_hi)} of {nh} = {pc(r.frac_k4p_hi, r.frac_k4p_hi_lo68, r.frac_k4p_hi_hi68)} %" if nh else "hi: none"
        lo_t = f"lo {int(r.k4p_lo)} of {nl} = {pc(r.frac_k4p_lo, r.frac_k4p_lo_lo68, r.frac_k4p_lo_hi68)} %" if nl else "lo: none"
        L.append(f"  {LAB[r.run_type]:<17}{'sample' if str(r.source).startswith('sample') else 'full scan':<11}{hi_t:<56}{lo_t}")
    L.append("")
    L.append("top modules of the high-A group (hi) and low-A group (lo):")
    for r in g.itertuples():
        src = "sample" if str(r.source).startswith("sample") else "full scan"
        L.append(f"  {LAB[r.run_type]:<17}{src:<10} hi: {r.top_modules_hi if isinstance(r.top_modules_hi, str) else '-'}")
        L.append(f"  {'':<27} lo: {r.top_modules_lo}")
    L.append("")
    L.append("w50 profile (P6 TSV, section w50_profile; 2-files/run sample): kept channels per w50 bin (N kept) and the fraction with k >= 2, % [Wilson 68 %]")
    wp = P6[P6.section == "w50_profile"]
    L.append(f"{'run type':<17}" + "".join(f"{f'w50 {a:g}-{b:g}':>37}" for a, b in [(0, 1.5), (1.5, 2), (2, 3), (3, 4), (4, 99)]))
    for rt in RTS:
        x = wp[wp.run_type == rt].sort_values("w50_lo")
        L.append(f"{LAB[rt]:<17}" + "".join(f"{f'N = {int(r.n_kept)}: {pc(r.frac_k_ge2, r.frac_k_ge2_lo68, r.frac_k_ge2_hi68)} %':>37}" for r in x.itertuples()))
    L.append("")
    L.append("Diagnosis")
    return L

def page_decisions(P):
    p3 = read("exotic_P3_sample_20261117.tsv")
    r3 = p3[(p3.section == "rate") & (p3["item"] == "any")].set_index("run_type")
    lk = p3[(p3.section == "f4_leakage") & (p3["item"] == "either")].set_index("run_type")
    nk = p3[(p3.section == "f4_leakage") & (p3["item"] == "n")].set_index("run_type")
    c3 = lambda rt: f"{g3(r3.loc[rt].percent)} % [{g3(r3.loc[rt].percent_lo68)}, {g3(r3.loc[rt].percent_hi68)}]"
    lkc = lambda rt: f"{lk.loc[rt].percent:.1f} % [{lk.loc[rt].percent_lo68:.1f}, {lk.loc[rt].percent_hi68:.1f}] ({int(lk.loc[rt]['count'])} of {int(nk.loc[rt]['count'])})"
    v = read(P10T)
    v = v[(v.table == "veto") & (v.population == "elastic") & (v.veto == "V3") & (v.axis == "overall")].set_index("dataset")
    v3 = lambda ds: f"{100 * v.loc[ds].frac:.4f} % [{100 * v.loc[ds].lo:.4f}, {100 * v.loc[ds].hi:.4f}]"
    sat = read(P6T)
    sat = sat[(sat.section == "saturation") & (sat.run_type == "2p24GeV_H2")].iloc[-1]
    gn = pd.read_csv(HERE / P9T, sep="\t", low_memory=False)
    gn = gn[(gn.table == "gain") & (gn.run_type == "2p24GeV_Ta_X17") & gn.module.isin(["W493", "W494", "W559", "W564", "W593", "W663"])]
    gm = gn.loc[gn.dg.abs().idxmax()]
    L = [
        "Decisions and open items (from the plan Log, notes/plan_pileup_fit_validation_20261116.md)", "",
        "Decisions",
        "  2026-11-17  user examples define four types; flags f1 (db > 8 sigma and tail < -5 sigma), f2 (w50), f3 (stretched single, BIC), f4 (clipped with a secondary) behind PILEUP_EXOTIC=1;",
        "              BASE_SLOPE / CLIP_SINGLE gated and off.",
        f"  2026-11-17  P3 on the 2-file sample: any-flag (Wilson 68 %) {c3('0p73GeV_H2')} / {c3('2p24GeV_H2')} / {c3('3p49GeV_H2')} in H2 (0.73 / 2.24 / 3.49 GeV), {c3('2p24GeV_Ta_X17')} in Ta X17;",
        f"              f4 leakage into the +-3 sigma elastic window (main or secondary in the window): {lkc('2p24GeV_H2')} at 2.24 GeV, {lkc('3p49GeV_H2')} at 3.49 GeV of the f4 channels (H2).",
        "  2026-11-18  P4 full scan (13 runs, 103 files, 6.0 h at 16 procs); per-type thresholds move little (only Ta db 8 -> 9.76); f2 truncation exclusion (t0 > 284 ns).",
        "  2026-11-18  f2 excludes clipped records and is evaluated after f4 (shared f2_eligible); P4 and P5 regenerated in place.",
        "  2026-11-18  f5 added (clipped, not f4); saturation by run type, module and radius; P6 added.",
        "  2026-11-18  USER: flags become mutually exclusive categories, priority f4 > f5 > f1 > f2 > f3 (supervisor's order, pending the user's objection); no separate f4-before-f3 rule.",
        "  2026-11-18  USER: template FWHM 23.5 ns (PbWO4, all run types), w50 > 3 = data FWHM > 70 ns; Ta w50 cut kept at 3 (tail fraction, not a distinct population); w50 > 3 examples added to P6.",
        "  2026-11-18  USER: run 24937 (3.49 GeV H2) is 'not a good run'; dropped from the exotic inputs (files stay on disk).  Saturation: up to 0.1 % of full-energy pulses for now (working bound; measured at 2.24 GeV H2 " + f"{pct(sat.pct_clip_E_gt_0p5, 3)} % [{pct(sat.lo68_E, 3)}, {pct(sat.hi68_E, 3)}], {int(sat.n_clip_E_gt_0p5)} of {int(sat.n_E_gt_0p5):,}" + "), keep f5 for the ifarm full scan.",
        "  2026-11-18  USER: this consolidated PDF after P6; then P7-P10 (cluster-level impact packages).",
        "  2026-11-18  USER (P7 feedback): single-module saturated clusters at about E_beam are fakes (high-gain module + non-shower single-channel pulse), to veto in P10;",
        "              0.73 GeV calibration sample restricted to R >~ 107 mm for P8/P9; Ta f1 rows outside any cluster left to P10.",
        "  2026-11-18  USER: gain check: W457 / W564 gains normal (ratio to median 1.00-1.03); W613 / W307 / W335 outside the r < 200 mm calibration pool;",
        "              saturated seeds in the elastic window are narrow non-template spikes, not high-gain modules.  Exotic study = one chapter of the write-up tex (figures in a new folder under writeup/figs), after P10.",
        "  2026-11-18  USER: saturated seeds in the elastic window are two populations (pages 10-11): (a) fakes, n_members <= 2; (b) real showers with a clipped high-gain seed, n_members >= 3 (W398 at 3.49 GeV).",
        "              P10 veto V3 = saturated seed AND <= 2 members: (a) only.  (b) is not vetoed: an energy-estimate issue for P8 / P9, reported per module as high-gain candidates.",
        "  2026-11-18  USER (P10 feedback): V3 (saturated seed with seed fraction > 0.95) adopted as policy on 2.24 / 3.49 GeV H2; V4 decided after P9; run 25800 beam current taken as 3.0 nA (elog, 'looks like yes');",
        "              W398 listed as the high-gain module for P8 / P9; the P11 write-up chapter drafted in parallel (P8/P9 subsection first a placeholder).",
        "  2026-11-18  P8 / P9 done (see pages 12 and 14); P10 done (page 13); P12 rates vs energy slice done (page 15).",
        "  2026-11-18  USER (standing rule, convention entry): every percentage, rate or fraction is quoted with its statistical uncertainty (Wilson 68 % or +-) or its count; applied to all pages here.",
        "  2026-11-18  USER: rates per flag (f1-f5, any) per energy slice, per channel and per cluster (P12, page 15).",
        "  2026-11-18  USER: policy request: do not use a kept channel in any flag f1-f5 (drop the channel's hit before clustering and the gain solve), and redo all calculations (P7-P10) under it: P13 (page 16).",
        *P13_DECISION(),
        f"  2026-11-18  GATE (P9/P10): 2.24 and 3.49 GeV H2 keep V3 (fakes only; loss {v3('2p24GeV_H2')} / {v3('3p49GeV_H2')} of the elastic window, Wilson 68 %), no V4, no in-place correction; 0.73 GeV H2 and Ta X17: no veto, no correction;",
        "              the flags stay event markers; f5 is kept for the farm (ifarm) scan.",
        "",
        "Open items",
        "  - P1-P3 pages below were made under the pre-exclusive definitions and include run 24937; not regenerated.",
        "  - f3 vs f4 priority was settled by the exclusive order; the clipped channels' role in the per-type threshold core is deferred.",
        "  - Statistics: scan npz kept only for 19 Ta X17 files; the other H2/Ta files are candidate-row pickles (counts per the pre-exclusion Log), so anything beyond the stored rows needs a rescan.",
        "  - Raw peak for the H2 full scan is missing (pickles lack event ids); H2 f2/f4 amplitude N is small (0.73 GeV f4 N = 2).  An amplitude-carrying rescan is needed.",
        "  - Saturation: the E > 0.5 E_beam ratio for H2 comes from the 2-file sample only (0.73 GeV: 2 clipped pulses, page 2 gives the interval); too small, keep f5 for the ifarm full scan.",
        "  - Count criterion on the f1-f5 colour maps (deferred, P6).",
        "  - Run-type assignment of 24937 looked Ta-like (check recorded 2026-11-18); excluded, not reassigned.",
        "  - 3.49 GeV H2 sample: runs 24698 and 24911 (the earlier one-run note is superseded, user 2026-11-18); thin, the 3.49 ratios and P6 rows rest on it.",
        "  - P13 done (page 16): drop policy not adopted; P7-P10 recomputed under it.",
        "  - P7 done (page 9-11); P8 (page 12), P10 (page 13) and P9 (page 14) done; the write-up chapter (P11) is being drafted separately.",
        "  - DEFERRED (P8/P9 Log): P8 rho = -1 credibility for the 81k Ta members; f5 mild-clipping threshold (nclip <= 3) to confirm; Ta inner-ring (block H) gains need a calibration that frees block H",
        f"    (prior-held in the Moller calibration, about half prior in ep: 6 of 36 entries free, max |dg| {100 * gm.dg:+.3f} +- {100 * gm.sigma_g:.3f} % (sigma_g), page 14); ifarm full scan (f5 kept for it); H2 small-N cells in P8.",
    ]
    P.text("Decisions and open items", L, size=8.8, sub="Dated from the plan Log; numbers elsewhere in this summary come from the TSVs.")


def p7_cells(F, ds):
    """judged cells (population elastic / moller, overall and E bins, category x seed / neighbour, n_flag >= 5): all, those with <fE> >= 0.1 %, and the row with the largest <fE> below 0.1 %."""
    c = F[(F.dataset == ds) & F.population.isin(["elastic", "moller"]) & F.axis.isin(["overall", "ebin"]) & (F.category != "any") & (F.role != "any") & (F.n_flag >= 5)]
    rest = c[c.ef_mean < 1e-3]
    return c, c[c.ef_mean >= 1e-3], (rest.loc[rest.ef_mean.idxmax()] if len(rest) else None)

def page_p7(P):
    t = read(P7T)
    F = t[t.table == "frac"]
    dss = [("0p73GeV_H2", "0.73 GeV H2 (sample)"), ("2p24GeV_H2", "2.24 GeV H2 (sample)"), ("3p49GeV_H2", "3.49 GeV H2 (sample)"), ("2p24GeV_Ta_X17", "2.24 GeV Ta X17 (full)")]
    g = lambda ds, pop, role, cat="any": F[(F.dataset == ds) & (F.population == pop) & (F.axis == "overall") & (F.category == cat) & (F.role == role)].iloc[0]
    L = ["Table 1a: clusters with at least one flagged module (exclusive categories f1-f5, any category), PbWO4-seeded clusters of the benchmark (control v3).  n flagged / n clusters = % [Wilson 68 %].",
         "Role: seed = largest member module, neighbour = other members.  Moller pairs exist for Ta X17 only.", ""]
    L.append(f"{'dataset':<24}{'population':<12}{'any role':>46}   {'seed':>46}   {'neighbour':>46}")
    for ds, lab in dss:
        for pop, pl in (("all", "all"), ("elastic", "elastic"), ("moller", "Moller")):
            if F[(F.dataset == ds) & (F.population == pop)].empty:
                continue
            row = f"{lab:<24}{pl:<12}"
            for role in ("any", "seed", "neighbour"):
                r = g(ds, pop, role)
                row += f"{f'{int(r.n_flag)} / {int(r.n_clusters)} = {pct(100 * r.frac, 3)} % [{pct(100 * r.frac_lo, 3)}, {pct(100 * r.frac_hi, 3)}]':>46}   "
            L.append(row)
    L += ["", "Table 1b: <fE> = mean over all clusters of the population of the energy fraction in flagged modules, % +/- standard error of the mean (n clusters as in 1a)", f"{'dataset':<24}{'population':<12}{'any role':>22}{'seed':>22}{'neighbour':>22}"]
    for ds, lab in dss:
        for pop, pl in (("all", "all"), ("elastic", "elastic"), ("moller", "Moller")):
            if F[(F.dataset == ds) & (F.population == pop)].empty:
                continue
            row = f"{lab:<24}{pl:<12}"
            for role in ("any", "seed", "neighbour"):
                r = g(ds, pop, role)
                row += f"{f'{pct(100 * r.ef_mean, 3)} +/- {pct(100 * r.ef_sem, 4)} %':>22}"
            L.append(row)
    cats = ["any", "f4 clipped-multi", "f5 saturated single", "f1 tail-baseline", "f2 wide", "f3 smooth-wide"]
    L += ["", "Fake-cluster candidates (seed flagged): rate per 1000 events [Poisson / Wilson 68 %], all clusters", f"{'dataset':<24}{'events':>10}" + "".join(f"{c:>26}" for c in ["any", "f4 clipped-multi", "f5 saturated single", "f1 tail-baseline", "f2 wide", "f3 smooth-wide"])]
    for ds, lab in dss:
        r0 = g(ds, "all", "seed")
        L.append(f"{lab:<24}{int(r0.n_events):>10d}" + "".join(f"{f'{g(ds, 'all', 'seed', c).per_kevent:.3f} [{g(ds, 'all', 'seed', c).rate_lo68:.3f}, {g(ds, 'all', 'seed', c).rate_hi68:.3f}]':>26}" for c in cats))
    sat = F_sat = t[t.table == "sat"]
    L += ["", "Saturated seeds (f4 or f5) in the ep elastic window: n / n elastic-window clusters = % [Wilson 68 %] (table sat; list on pages 10-11, recomputed from the P7 cache)"]
    for ds, lab in dss:
        r = sat[(sat.dataset == ds) & (sat.population == "elastic") & (sat.axis == "overall")].iloc[0]
        L.append(f"  {lab:<24}{int(r.n_flag):>4d} / {int(r.n_clusters):<8d} = {pct(100 * r.frac, 4)} % [{pct(100 * r.frac_lo, 4)}, {pct(100 * r.frac_hi, 4)}]")
    L += ["", "Verdict for P8 / P9 (calibration samples: elastic and Moller, overall and per E_g / E_beam bin; cells with >= 5 flagged clusters; criterion <fE> >= 0.1 %)"]
    for ds, lab in dss:
        c, hit, rest = p7_cells(F, ds)
        rtxt = f"{pct(100 * rest.ef_mean, 3)} +/- {pct(100 * rest.ef_sem, 4)} %" if rest is not None else "-"
        if hit.empty:
            L.append(f"  {lab:<24}no cell reaches 0.1 % (largest <fE> {rtxt})")
        else:
            ov = hit[hit.axis == "overall"]
            L.append(f"  {lab:<24}overall cells >= 0.1 %: " + "; ".join(f"{r.population} {r.category} {r.role} {pct(100 * r.ef_mean, 3)} +/- {pct(100 * r.ef_sem, 4)} % ({int(r.n_flag)} of {int(r.n_clusters)})" for r in ov.itertuples())
                     + f"; E-bin cells >= 0.1 %: {int((hit.axis == 'ebin').sum())}; every other cell < 0.1 % (largest {rtxt})")
    L += ["  => P8 / P9 priority: the Ta X17 elastic f1-neighbour and f3-seed cells; no H2 cell needs a correction at the 0.1 % level (the <fE> bound is a population average; P8 / P9 judge per module).",
          "  0.73 GeV elastic window restricted to angles where ep and Moller separate (R >~ 107 mm, the midpoint of E_moller / E' and 1 below 0.90).",
          "", "Source: " + P7T + " (table frac: overall rows, columns n_flag, n_clusters, frac, frac_lo, frac_hi, ef_mean, ef_sem, per_kevent, rate_lo68, rate_hi68; table sat: saturated-seed rows); saturated-seed list: " + P7S + "."]
    P.text("P7 cluster level: flagged modules in reconstructed clusters", L, size=7.0, sub="Join of the exclusive exotic categories (P4) to the benchmark's gated clusters; 2-file H2 sample, Ta X17 full scan; run 24937 excluded.")

def page_sat_table(P, S):
    L = ["Elastic-window clusters (ep, +-3 sigma) with a saturated (f4 or f5) seed, 2.24 and 3.49 GeV H2 sample.  n_members = member modules of the cluster; seed fraction = seed energy / cluster energy (gain-corrected);",
         "k = fitted pulses; A_main = largest fitted amplitude (ADC), A/FS = A_main / (4095 - ped); chi2/dof of the fit.  pop: (a) n_members <= 2, (b) n_members >= 3.", ""]
    L.append(f"{'run type':<12}{'run':>6}{'event':>11}{'module':>8}{'cat':>4}{'n_mem':>6}{'seed frac':>10}{'k':>3}{'A_main':>9}{'A/FS':>7}{'chi2/dof':>10}{'nclip':>6}{'w50':>6}  pop")
    for r in S.sort_values(["dataset", "pop", "module", "run", "evn"]).itertuples():
        L.append(f"{LAB[r.dataset]:<12}{r.run:>6d}{r.evn:>11d}{r.module:>8}{'f' + str(r.scat):>4}{r.n_members:>6d}{r.seed_frac:>10.4f}{r.k:>3d}{r.a_main:>9.0f}{r.a_over_fs:>7.2f}{r.chi2_dof:>10.1f}{r.nclip:>6d}{r.w50:>6.2f}  {r.pop}")
    L += ["", "Source: " + P7S + " (cache exotic_p7_cache.npz joined to the exotic_sample npz; script exotic_p7_saturated_seeds.py)."]
    P.text("Saturated seeds in the elastic window (P7): the clusters", L, size=7.4)


def page_sat_pops(P, S):
    L = ["Counts per population and run type (saturated seed in the elastic window); all: k of N elastic-window clusters, % [Wilson 68 %]; (a), (b): share of the saturated seeds, k of N_sat, % [Wilson 68 %]"]
    for ds, lab in (("2p24GeV_H2", "2.24 GeV H2"), ("3p49GeV_H2", "3.49 GeV H2")):
        x = S[S.dataset == ds]
        r0 = x.iloc[0]
        ka, kb = int((x['pop'] == 'a').sum()), int((x['pop'] == 'b').sum())
        L.append(f"{lab:<14}all: {len(x)} of {int(r0.n_elastic):,} = {pct(100 * r0.sat_frac, 4)} % [{pct(100 * r0.sat_lo68, 4)}, {pct(100 * r0.sat_hi68, 4)}]   (a) n_mem <= 2: {ka} of {len(x)} = {ci1(ka, len(x), 1)}   (b) n_mem >= 3: {kb} of {len(x)} = {ci1(kb, len(x), 1)}")
    L += ["", "Counts per module (a) | (b)"]
    for pop, nm in (("a", "(a)"), ("b", "(b)")):
        for ds, lab in (("2p24GeV_H2", "2.24 GeV H2"), ("3p49GeV_H2", "3.49 GeV H2")):
            c = S[(S.dataset == ds) & (S["pop"] == pop)].module.value_counts()
            L.append(f"  {nm} {lab:<12}" + ("  ".join(f"{m} {n}" for m, n in c.items()) if len(c) else "none"))
    a, b = S[S["pop"] == "a"], S[S["pop"] == "b"]
    w = b[b.module == "W398"]
    o = b[b.module != "W398"]
    L += ["", "Populations"]
    L += wrap(f"(a) fakes: n_members <= 2, i.e. one or two modules carry the whole cluster (seed fraction >= {a.seed_frac.min():.3f}); non-shower single-channel pulses with k {int(a.k.min())}-{int(a.k.max())}, chi2/dof {a.chi2_dof.min():.0f}-{a.chi2_dof.max():.0f}, "
              f"fitted A / full scale {a.a_over_fs.min():.2f}-{a.a_over_fs.max():.2f} ({int((a.a_over_fs <= 1.06).sum())} of {len(a)} rows at {a.a_over_fs[a.a_over_fs <= 1.06].min():.2f}-{a.a_over_fs[a.a_over_fs <= 1.06].max():.2f}, the rest above full scale).  {int((a.module == 'W457').sum())} of {len(a)} are on W457.  P10's V3 veto applies to population (a) only.")
    L += wrap(f"(b) real showers with a clipped high-gain seed: n_members >= 3.  The W398 rows at 3.49 GeV ({len(w)} clusters, {w.n_members.min()}-{w.n_members.max()} members) have seed fraction {w.seed_frac.min():.2f}-{w.seed_frac.max():.2f}, k = {int(w.k.min())}{"" if w.k.min() == w.k.max() else "-" + str(int(w.k.max()))}, "
              f"chi2/dof {w.chi2_dof.min():.1f}-{w.chi2_dof.max():.1f} (template-like, kept).")
    L += wrap(f"Not covered by that description: the other {len(o)} (b) rows ({', '.join(o.module)}; {o.n_members.min()}-{o.n_members.max()} members, seed fraction {o.seed_frac.min():.3f}-{o.seed_frac.max():.3f}, k {int(o.k.min())}-{int(o.k.max())}, chi2/dof {o.chi2_dof.min():.0f}-{o.chi2_dof.max():.0f}) "
              "look like (a) fakes with one or two low-energy neighbours; they are in (b) by the n_members >= 3 rule and are not vetoed unless the P10 veto definition says otherwise.")
    rec = b.module.value_counts()
    rec = rec[rec >= 2]
    L += wrap("Modules recurring in (b), high-gain candidates: " + (", ".join(f"{m} ({n} clusters, {', '.join(sorted(LAB[d] for d in set(b[b.module == m].dataset)))})" for m, n in rec.items()) or "none") +
              ".  Gain check (plan Log 2026-11-18): W457 / W564 gains normal; W613 / W307 / W335 outside the r < 200 mm calibration pool.")
    n24, n349 = S[S.dataset == "2p24GeV_H2"], S[S.dataset == "3p49GeV_H2"]
    L += ["", f"P7 report: 23 (2.24 GeV, all <= 2 modules) and 18 (3.49 GeV, 9 of 18 with <= 2 modules); here {len(n24)} ({int((n24.n_members <= 2).sum())} with <= 2 modules) and {len(n349)} ({int((n349.n_members <= 2).sum())} with <= 2 modules).", "Source: " + P7S + "."]
    P.text("Saturated seeds in the elastic window (P7): populations (a) and (b)", L, size=8.6)


# ----------------------------------------------------------------------------------------------------------------------------------------------
# P8 / P10 / P9 front-matter pages (2026-11-18): numbers read from the three TSVs
P8T = "exotic_P8_energy_error_20261118.tsv"
P9T = "exotic_P9_calibration_impact_20261118.tsv"
P10T = "exotic_P10_yield_impact_20261118.tsv"
RTL = {"0p73GeV_H2": "0.73 H2", "2p24GeV_H2": "2.24 H2", "3p49GeV_H2": "3.49 H2", "2p24GeV_Ta_X17": "2.24 Ta X17", "Ta_sample": "Ta X17 (sample)"}


def _one(df, **kw):
    for k, v in kw.items():
        df = df[df[k] == v]
    return df.iloc[0] if len(df) else None


def page_p8(P):
    d = read(P8T)
    sc = d[d.table == "scale"]
    mo = d[d.table == "module"]
    dg = d[d.table == "diag"]

    def S(rt, pop, cat, var, ring):
        r = _one(sc, run_type=rt, population=pop, category=cat, variant=var, ring=ring)
        return r

    def sfmt(r):
        return "-" if r is None or r.n_flag == 0 else f"{100 * r.S:+.4f} +- {100 * r.S_err:.4f}"

    L = ["Energy-scale shift S = sum(dE) / sum(E_cluster) over the elastic (ep, +-3 sigma) or Moller-pair clusters, dE = credible alternative minus default on the clustered hit's own pulse; in %, +- = stat. error.",
         "Source: " + P8T + " (table scale; module; diag).  391,292 flagged cluster-member records refitted (diag join row: default reproduced for all 391,292 records, repro = 1); alternatives: bs = sloped baseline, st = stretched single, k1 = single-pulse extrapolation, int = integral, cs = clip-aware.", "",
         "1. Scale table (table scale)",
         f"  {'run type':<12}{'pop.':<9}{'category [variant]':<34}{'n_flag':>8}{'S all rings (%)':>24}{'S ring 0-120 mm (%)':>26}"]
    rows = [("2p24GeV_Ta_X17", "elastic", "f1 tail-baseline", "bs"), ("2p24GeV_Ta_X17", "elastic", "f1 tail-baseline [hit=main]", "bs"), ("2p24GeV_Ta_X17", "elastic", "f1 tail-baseline [hit!=main]", "bs"),
            ("2p24GeV_Ta_X17", "elastic", "f3 smooth-wide", "st"), ("2p24GeV_Ta_X17", "elastic", "f2 wide", "k1"),
            ("2p24GeV_Ta_X17", "moller", "f1 tail-baseline", "bs"), ("2p24GeV_Ta_X17", "moller", "f3 smooth-wide", "st"), ("2p24GeV_Ta_X17", "moller", "f2 wide", "k1"),
            ("2p24GeV_H2", "elastic", "f1 tail-baseline", "bs"), ("3p49GeV_H2", "elastic", "f1 tail-baseline", "bs"),
            ("0p73GeV_H2", "elastic", "f3 smooth-wide", "st"), ("2p24GeV_H2", "elastic", "f3 smooth-wide", "st"), ("3p49GeV_H2", "elastic", "f3 smooth-wide", "st"),
            ("3p49GeV_H2", "elastic", "f5 saturated single", "k1"), ("3p49GeV_H2", "elastic", "f5 saturated single", "int"), ("2p24GeV_H2", "elastic", "f4 clipped-multi", "cs"), ("3p49GeV_H2", "elastic", "f4 clipped-multi", "cs")]
    for rt, pop, cat, var in rows:
        a, b = S(rt, pop, cat, var, "all"), S(rt, pop, cat, var, "0-120")
        nf = "-" if a is None else f"{int(a.n_flag)}"
        tag = "  (H2 f3: not credible)" if rt != "2p24GeV_Ta_X17" and cat.startswith("f3") else ""
        L.append(f"  {RTL[rt]:<12}{pop:<9}{(cat + ' [' + var + ']'):<34}{nf:>8}{sfmt(a):>24}{sfmt(b):>26}{tag}")
    # credible H2 cells (f1 bs, f5 k1, f5 int) and the non-credible small-N cell
    h2 = sc[(sc.run_type != "2p24GeV_Ta_X17") & (sc.n_flag > 0) & (((sc.variant == "bs") & (sc.category == "f1 tail-baseline")) | sc.variant.isin(["k1", "int"]) & sc.category.str.startswith("f5"))]
    hm = h2.loc[h2.S.abs().idxmax()]
    big = sc[(sc.run_type != "2p24GeV_Ta_X17") & (sc.n_flag > 0) & (sc.S.abs() >= 1e-3)]
    h10 = h2[h2.n_flag >= 10]
    h10m = h10.loc[h10.S.abs().idxmax()]
    L += ["", f"  Credible H2 cells (f1 sloped baseline, f5 k1 / integral; {len(h2)} cells with flagged clusters): with >= 10 flagged clusters ({len(h10)} cells) the largest |S| = {100 * abs(h10m.S):.4f} +- {100 * h10m.S_err:.4f} % "
          f"({RTL[h10m.run_type]} {h10m.category} [{h10m.variant}], ring {h10m.ring}, {int(h10m.n_flag)} flagged), below 0.01 %; over all {len(h2)} cells the largest is {100 * abs(hm.S):.4f} +- {100 * hm.S_err:.4f} % "
          f"({RTL[hm.run_type]} {hm.category} [{hm.variant}], ring {hm.ring}, {int(hm.n_flag)} flagged, S = {100 * hm.S:+.4f} +- {100 * hm.S_err:.4f} %, consistent with 0 at about {abs(hm.S) / hm.S_err:.1f} sigma); none reaches 0.1 %.",
          "  Non-credible H2 cells at or above 0.1 %: " + ("; ".join(f"{RTL[r.run_type]} ring {r.ring} {r.category} [{r.variant}]: {int(r.n_flag)} flagged, S = {100 * r.S:+.2f} +- {100 * r.S_err:.2f} %" for r in big.drop_duplicates(["run_type", "ring", "category"]).itertuples()
                                                                  if r.category != "all (primary)") or "none") + "  (statistics-limited, deferred).",
          "", "2. Module table: Ta elastic clusters seeded by the inner-ring (block H) modules, all categories, primary variants (table module, columns S_el_seeded +- S_el_err, n_el_seeded)",
          f"  {'module':<8}{'n flagged rows':>16}{'n elastic seeded':>18}{'S (%)':>22}"]
    inner = ["W493", "W494", "W559", "W564", "W593", "W663"]
    vals = []
    for m in inner:
        r = _one(mo, run_type="2p24GeV_Ta_X17", module=m, category="any")
        vals.append(100 * r.S_el_seeded)
        L.append(f"  {m:<8}{int(r.n):>16}{int(r.n_el_seeded):>18}{100 * r.S_el_seeded:>+14.4f} +- {100 * r.S_el_err:.4f}")
    errs = [100 * _one(mo, run_type="2p24GeV_Ta_X17", module=m, category="any").S_el_err for m in inner]
    L.append(f"  range over the six inner-ring modules: {min(vals):+.4f} to {max(vals):+.4f} % (each +- {min(errs):.4f} to {max(errs):.4f} %, see the table)  (about {abs(max(vals)):.3f} to {abs(min(vals)):.3f} % low).  Outer Ta modules (W527, W560, W594, W597, W629, W630) have 24-33 seeded elastic clusters each: statistics-limited.")
    w = _one(mo, run_type="3p49GeV_H2", module="W398", category="f5 saturated single")
    L.append(f"  High-gain candidate W398 (3.49 GeV): {int(w.n)} f5 records, all seeds in the elastic window, dE(k1) = {w.dE_med:+.1f} MeV (k = 1 extrapolation adds nothing); the integral estimate dE = -23 to -34 MeV "
             "(-0.7 to -1.0 % of E_cluster) is a lower bound (P8 PDF p.3 text, not in the TSV tables).")
    L += ["", "3. Credible alternative per category (table diag, median over the run type; A ratio = alternative A_main / default A_main; ok = fraction of records where the alternative converged)",
          f"  {'category':<22}{'credible alternative':<44}{'evidence (diag)':<100}"]
    def dr(rt, cat):
        return _one(dg, run_type=rt, category=cat)
    ta1, h1 = dr("2p24GeV_Ta_X17", "f1 tail-baseline"), dr("2p24GeV_H2", "f1 tail-baseline")
    ta3, h3 = dr("2p24GeV_Ta_X17", "f3 smooth-wide"), dr("3p49GeV_H2", "f3 smooth-wide")
    f5h = dr("3p49GeV_H2", "f5 saturated single")
    ta2, h2f = dr("2p24GeV_Ta_X17", "f2 wide"), dr("2p24GeV_H2", "f2 wide")
    f4h = dr("2p24GeV_H2", "f4 clipped-multi")
    L += [f"  {'f1 tail-baseline':<22}{'sloped baseline (bs): credible everywhere':<44}bs_ok = {ta1.bs_ok:.2f} (Ta, N = {int(ta1.n):,}), {h1.bs_ok:.2f} (2.24 H2, N = {int(h1.n):,}); A_main ratio {ta1.bs_amain_ratio_med:.3f} (Ta), {h1.bs_amain_ratio_med:.3f} (2.24 H2); {ci1(round(ta1.k_drop * ta1.n), ta1.n, 3)} of Ta rows lose a pulse (N = {int(ta1.n):,})",
          f"  {'f3 smooth-wide':<22}{'stretched single (st) in Ta only':<44}Ta: stretch s median {ta3.s_med:.1f}, A ratio {ta3.st_amain_ratio_med:.2f};  H2 (3.49): s median {h3.s_med:.1f}, A ratio {h3.st_amain_ratio_med:.2f} (swallows real close pile-up): default in H2",
          f"  {'f2 wide':<22}{'none (default kept)':<44}k1 / st alternatives not credible: H2 f2 = merged k 4-6 trains (s median {h2f.s_med:.0f}, st A ratio {h2f.st_amain_ratio_med:.1f}); Ta A ratio {ta2.k1_amain_ratio_med:.2f}",
          f"  {'f5 saturated single':<22}{'k = 1 extrapolation for mild clipping':<44}clipped-area fraction {100 * f5h.clip_area_frac_med:.1f} % (median over N = {int(f5h.n)} records, 3.49 H2), k1 A ratio {f5h.k1_amain_ratio_med:.3f}, integral {f5h.int_amain_ratio_med:.3f}; deep saturation not correctable: veto (V3)",
          f"  {'f4 clipped + sec.':<22}{'none (default kept)':<44}clip-aware refit reproduces the default (cs A ratio {f4h.cs_amain_ratio_med:.2f}, ok {f4h.cs_ok:.2f}); too few records in H2 (0.73: 4, 2.24: {int(f4h.n)}, 3.49: {int(dr('3p49GeV_H2', 'f4 clipped-multi').n)})", "",
          f"Reading: only Ta X17 moves measurably: elastic f1 {sfmt(S('2p24GeV_Ta_X17', 'elastic', 'f1 tail-baseline', 'bs', 'all'))} %, f3 {sfmt(S('2p24GeV_Ta_X17', 'elastic', 'f3 smooth-wide', 'st', 'all'))} %, Moller f1 {sfmt(S('2p24GeV_Ta_X17', 'moller', 'f1 tail-baseline', 'bs', 'all'))} % (inner ring {sfmt(S('2p24GeV_Ta_X17', 'moller', 'f1 tail-baseline', 'bs', '0-120'))} %), inner-ring modules {min(vals):+.4f} to {max(vals):+.4f} % (+- {min(errs):.4f} to {max(errs):.4f} %).  Credible alternatives: f1 sloped baseline; f3 stretched (Ta only);",
          "default for f2, f4 and H2 f3; f5 k1 only for mild clipping.  Deferred (plan Log): second-order propagation (done in P9), small-N H2 cells, rho = -1 credibility of the 81k Ta members."]
    P.text("P8 energy error: scale shifts and credible alternatives", L, size=7.0, sub="Per category, Ta X17 first; numbers read from " + P8T)


def _w(r, nd=3):
    return f"{100 * r.frac:.{nd}f} [{100 * r.lo:.{nd}f},{100 * r.hi:.{nd}f}] ({int(r.k)}/{int(r.n)})"


def page_p10(P):
    d = read(P10T)
    v = d[d.table == "veto"]

    def V(ds, pop, ve, axis="overall", bn="all"):
        return _one(v, dataset=ds, population=pop, veto=ve, axis=axis, bin=bn)

    L = ["Source: " + P10T + " (tables veto, ooc_class, ooc_widen, ooc_gate, occ_fit, occ_run, satmod).  Elastic = P7 +-3 sigma window; Moller = selected pairs (lost when either member is vetoed); hiE = clusters with E_g > 0.5 E_beam.  Wilson 68 %.", "",
         "Veto definitions (per cluster; seed = largest member; categories = exclusive f4 > f5 > f1 > f2 > f3)",
         "  V1  seed flagged (any category)        V2  any member flagged        V3  saturated seed (f4 or f5) AND seed energy fraction > 0.95 (the fakes)        V3s  = V3 with 1 member (sub-row)",
         "  V4  = V3 + seed in f1 or f3             V3b  saturated seed in a normal cluster (seed fraction <= 0.95): real showers (W398), NOT vetoed, listed only", "",
         "1. Loss in the ep elastic window, % (68 % Wilson interval) (k/n)",
         f"  {'run type':<16}{'V2 any member':<34}{'V4 = V3 + f1/f3 seed':<34}{'V3 saturated seed':<34}{'V3s':>6}{'V3b':>6}"]
    for ds in ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]:
        L.append(f"  {RTL[ds]:<16}{_w(V(ds, 'elastic', 'V2')):<34}{_w(V(ds, 'elastic', 'V4')):<34}{_w(V(ds, 'elastic', 'V3')):<34}{int(V(ds, 'elastic', 'V3s').k):>6}{int(V(ds, 'elastic', 'V3b').k):>6}")
    dif = [ds for ds in ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"] if V(ds, 'elastic', 'V1').k != V(ds, 'elastic', 'V4').k]
    L.append("  V1 (seed flagged) equals V4 in the elastic window except at " + ", ".join(f"{RTL[ds]}: V1 {_w(V(ds, 'elastic', 'V1'))}" for ds in dif) + " (the difference is the V3b clipped-seed showers, seed fraction <= 0.95, which V4 does not veto).")
    L += ["", "2. Ta X17 Moller pairs, % (Wilson 68 %) (k/n)", f"  {'dataset':<18}{'V2 any member':<34}{'V4':<34}{'V3':<34}"]
    for ds in ["2p24GeV_Ta_X17", "Ta_sample"]:
        L.append(f"  {RTL[ds]:<18}{_w(V(ds, 'moller_pair', 'V2')):<34}{_w(V(ds, 'moller_pair', 'V4')):<34}{_w(V(ds, 'moller_pair', 'V3')):<34}")
    r0 = V("2p24GeV_Ta_X17", "elastic", "V4", "ring", "0-120")
    rest = max((V("2p24GeV_Ta_X17", "elastic", "V4", "ring", b) for b in ["120-170", "170-230", "230-320", "320-450"]), key=lambda x: x.frac)
    L += ["", f"  Ta elastic V4 inner ring (0-120 mm): {_w(r0)} %; other rings: largest {_w(rest)} % (ring {rest['bin']} mm).  Ta Moller V4 inner-member ring 0-120: {_w(V('2p24GeV_Ta_X17', 'moller_pair', 'V4', 'ring', '0-120'))} %.",
          "", "3. Contamination removed by V3 (fakes: vetoed saturated seeds with seed fraction > 0.95; for V3 the 'loss' is not real ep yield)"]
    for ds in ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]:
        L.append(f"  {RTL[ds]:<14} elastic window: {int(V(ds, 'elastic', 'V3').k):>5} fakes removed (V3b kept: {int(V(ds, 'elastic', 'V3b').k)});   all clusters with E_g > 0.5 E_beam: {int(V(ds, 'hiE', 'V3').k):>5} removed of {int(V(ds, 'hiE', 'V3').n)} "
                 f"(= {100 * V(ds, 'hiE', 'V3').frac:.4f} % [{100 * V(ds, 'hiE', 'V3').lo:.4f}, {100 * V(ds, 'hiE', 'V3').hi:.4f}]), V3b kept {int(V(ds, 'hiE', 'V3b').k)}")
    sm = d[d.table == "satmod"]
    w398 = sm[(sm.dataset == "3p49GeV_H2") & (sm.bin == "W398")].iloc[0]
    w457 = sm[(sm.dataset == "2p24GeV_H2") & (sm.bin == "W457")].iloc[0]
    L.append(f"  satmod table: 2.24 GeV fakes dominated by W457 ({int(w457.n)} of {int(sm[sm.dataset == '2p24GeV_H2'].n.sum())} saturated-seed clusters); W398 (3.49 GeV): {int(w398.n)} clusters, {int(w398.k_in_window)} in the elastic window, the V3b showers (high-gain module).")
    of = d[d.table == "occ_fit"]
    f = _one(of, dataset="Ta split level, x = mean kept PbWO4 channels per event")
    fh = _one(of, dataset="H2 split level (all three energies), x = mean kept PbWO4 channels per event")
    orr = d[d.table == "occ_run"].set_index("dataset")
    r58, r62, r61 = orr.loc["25800"], orr.loc["26230"], orr.loc["26138"]
    L += ["", "4. Ta inner-ring (0-120 mm) f1 rate vs occupancy (table occ_fit, occ_run; y = a + b x, weighted LSQ, x = mean kept PbWO4 channels per event)",
          f"  Ta: slope b = {100 * f.b:.4f} +- {100 * f.sb:.4f} %/kept channel ({int(f.n)} splits of runs 25800, 26230, 26138), elasticity {f.elasticity:.2f} +/- {f.elasticity_err:.2f}, chi2/dof {f.chi2 / f.dof:.1f} (run-level scatter above the binomial error; 3 runs only).   H2 (18 splits): {100 * fh.b:.4f} +- {100 * fh.sb:.4f} %/channel.",
          "  Run rates (inner-ring f1, % [Wilson 68 %], k of N):",
          f"    25800 ({r58.a:g} nA, elog; RCDB blank) {_w(r58)} = {100 * r58.frac:.3f} %, {r58.b:.1f} kept ch/ev",
          f"    26230 ({r62.a:g} nA) {_w(r62)} = {100 * r62.frac:.3f} %, {r62.b:.1f} kept ch/ev",
          f"    26138 ({r61.a:g} nA) {_w(r61)} = {100 * r61.frac:.3f} %, {r61.b:.1f} kept ch/ev",
          f"  26230 / 25800: rate x{r62.frac / r58.frac:.2f} +/- {(r62.frac / r58.frac) * np.hypot((r62.hi - r62.lo) / 2 / r62.frac, (r58.hi - r58.lo) / 2 / r58.frac):.2f}, current x{r62.a / r58.a:.2f}, kept channels x{r62.b / r58.b:.2f}.  26138 has the same current as 26230 but fewer kept channels and a lower rate: current, occupancy, trigger version and inner-layer HV are confounded.",
          "", "5. Out-of-cluster Ta f1 rows (table ooc_class; ooc_widen; ooc_gate)"]
    oc = d[d.table == "ooc_class"]
    tot = int(oc.n.sum())
    L.append(f"  {tot} Ta f1 rows lie outside any cluster (sum of the six classes):")
    for r in oc.itertuples():
        L.append(f"    {r.bin:<92}{int(r.n):>8}  {100 * r.frac:>6.3f} % [{100 * r.lo:.3f}, {100 * r.hi:.3f}] (of N = {tot:,})")
    ow = d[d.table == "ooc_widen"].set_index("bin")
    L.append(f"  A trigger window wider by 32 / 64 ns would let {int(ow.loc['w=32'].k_in_window)} / {int(ow.loc['w=64'].k_in_window)} of these rows seed 30+ MeV clusters in the window (table ooc_widen).  Not in the TSV (P10 PDF p.11 text): 557,506 f1 rows in total,")
    L.append("  191,223 in a cluster, 114,637 of them with the flagged pulse > 5 ns outside the hit window (input to P8); only 120 outside rows reach 500 MeV.  No action for the yields.")
    P.text("P10 yield impact: veto definitions and losses", L, size=7.0, sub="Veto losses per run type, Ta Moller pairs, contamination removed, Ta inner-ring occupancy slope, out-of-cluster f1 rows; numbers read from " + P10T)


def _zfmt(r):
    if r is None:
        return "-"
    if r["se"] == 0 or pd.isna(r["z"]):
        return "0 (identical)"
    return f"{r['diff']:+.4f} +- {r['se']:.4f} (z {r['z']:+.2f})"


def page_p9(P):
    d = read(P9T) if False else pd.read_csv(HERE / P9T, sep="\t", low_memory=False)
    be = d[d.table == "bench"]
    pa = d[d.table == "paired"]
    bp = d[d.table == "bench_paired"].copy()
    bp["z"] = bp.value / bp.se
    L = ["Arms: A default control v3; V3 veto (saturated seed, seed fraction > 0.95; the policy on 2.24 / 3.49 GeV H2); V4 = V3 + f1/f3 seeds; C corrected (first-order P8 credible dE: f1 sloped baseline, f3 stretched in Ta only, f5 k1 mild clipping).",
         "Units: sigma/E = Crystal Ball sigma/E (%); mu = peak/E' - 1 (%); differences in percentage points.  Out-of-sample (oos) = the benchmark's own protocol.  0.73 GeV: calibration sample at R >= 107 mm (user decision).",
         "Source: " + P9T + " (tables bench, bench_paired, paired, gain).", "",
         "1. Arm A out-of-sample and the paired differences X - A (table paired: same clusters refitted pooled, bootstrap error, z = diff / se; ep pool, overall; Ta X17 Moller: summed-pair fit)",
         f"  {'run type':<22}{'A oos sigma/E (%)':>20}{'d sigma/E V3':>34}{'d sigma/E V4':>34}{'d sigma/E C':>34}",
         ]
    Ls2 = []
    rows = []
    for rt in ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]:
        a = _one(be, run_type=rt, arm="A", quantity="res_cb_oos")
        row_s = f"  {RTL[rt] + ' ep':<22}{f'{a.value:.3f} +- {a.se:.3f}':>20}"
        row_m = f"  {RTL[rt] + ' ep':<22}{f'{_one(be, run_type=rt, arm='A', quantity='mu_ep_oos').value:+.3f} +- {_one(be, run_type=rt, arm='A', quantity='mu_ep_oos').se:.3f}':>20}"
        for c in ["V3", "V4", "C"]:
            rs = _one(pa, run_type=rt, calib="ep", comparison=f"{c}-A", protocol="oos", axis="overall", bin="all", quantity="res_cb")
            rm = _one(pa, run_type=rt, calib="ep", comparison=f"{c}-A", protocol="oos", axis="overall", bin="all", quantity="mu_ep")
            row_s += f"{_zfmt(rs):>34}"
            row_m += f"{_zfmt(rm):>34}"
        rows.append(row_s)
        Ls2.append(row_m)
    a = _one(be, run_type="2p24GeV_Ta_X17", arm="A", quantity="sum:res_pooled", protocol="out-of-sample")
    am = _one(be, run_type="2p24GeV_Ta_X17", arm="A", quantity="sum:peak", protocol="out-of-sample")
    rs, rm = f"  {'2.24 Ta X17 Moller sum':<22}{f'{a.value:.3f} +- {a.se:.3f}':>20}", f"  {'2.24 Ta X17 Moller sum':<22}{f'{am.value:+.3f} +- {am.se:.3f}':>20}"
    for c in ["V3", "V4", "C"]:
        rs += f"{_zfmt(_one(pa, run_type='2p24GeV_Ta_X17', calib='moller', comparison=f'{c}-A', protocol='oos', axis='overall', bin='sum', quantity='res')):>34}"
        rm += f"{_zfmt(_one(pa, run_type='2p24GeV_Ta_X17', calib='moller', comparison=f'{c}-A', protocol='oos', axis='overall', bin='sum', quantity='peak')):>34}"
    rows.append(rs)
    Ls2.append(rm)
    L += rows + ["", f"  {'run type':<22}{'A oos mu (%)':>20}{'d mu V3':>34}{'d mu V4':>34}{'d mu C':>34}"] + Ls2
    cov = d[d.table == "coverage"]
    L.append("  (A: the benchmark's split mean +- standard error (table bench); differences: the pooled paired refit (table paired; Moller sum A = pooled fit).  'identical': the arm changes no cluster in the fit (table coverage in_pool: 0.73 GeV V3 "
             f"{int(_one(cov, run_type='0p73GeV_H2', arm='V3', quantity='in_pool').value)}, 0.73 GeV C {int(_one(cov, run_type='0p73GeV_H2', arm='C', quantity='in_pool').value)}, Ta V3 {int(_one(cov, run_type='2p24GeV_Ta_X17', arm='V3', quantity='in_pool').value)}).)")
    # bench_paired cross-check
    L += ["", "2. Cross-check on the benchmark protocol (table bench_paired: split-mean oos difference +- se, z = value / se)",
          f"  {'run type':<22}{'d sigma/E V3':>36}{'d sigma/E V4':>36}{'d sigma/E C':>36}{'d mu V3':>36}{'d mu V4':>36}{'d mu C':>36}"]
    for rt in ["0p73GeV_H2", "2p24GeV_H2", "3p49GeV_H2", "2p24GeV_Ta_X17"]:
        row = f"  {RTL[rt] + ' ep':<22}"
        for q in ["res_cb_oos", "mu_ep_oos"]:
            for c in ["V3", "V4", "C"]:
                r = _one(bp, run_type=rt, quantity=q, comparison=f"{c}-A")
                row += f"{('0 (identical)' if r.se == 0 else f'{r.value:+.4f} +- {r.se:.4f} (z {r.z:+.2f})'):>36}"
        L.append(row)
    zs = pa[(pa.protocol == "oos") & pa.bin.isin(["all", "sum"]) & pa.quantity.isin(["res_cb", "mu_ep", "tail_lo3", "tail_hi3", "res", "peak", "tail_lo3"])].z.abs()
    zi = pa[(pa.protocol == "in") & pa.bin.isin(["all", "sum"])].z.abs()
    zr = pa[(pa.axis == "ring") & pa.z.notna()].z.abs()
    iz = pa[(pa.protocol == "oos") & pa.bin.isin(["all", "sum"])].loc[zs.idxmax()]
    L.append(f"  Largest |z| in the overall out-of-sample paired rows (sigma/E, mu, both tails; ep and Moller sum): {zs.max():.2f} ({RTL[iz.run_type]} {iz.comparison} {iz.quantity}); in-sample overall (table paired) {zi.max():.2f}; per-ring rows {zr.max():.2f}.  Nothing in table paired reaches 2; the verdict rests on the oos rows (bench_paired in-sample rows can exceed 2, e.g. 3.49 C-A mu).")
    # gain
    g = d[d.table == "gain"].copy()
    g["zg"] = g.dg / g.sigma_g
    f = g[g.free == 1]
    L += ["", "3. Per-module gain change dg = g_X - g_A (table gain; free modules only; zg = dg / sigma_g of arm A; solve: ep pool / ring-scoped, Moller per-run)",
          f"  {'run type, calibration':<30}{'arm':<6}{'free gains':>11}{'max|dg| +- sg (%)':>18}{'max |zg|':>10}  module of max |zg| (run), dg/g_A"]
    for rt, cal in [("0p73GeV_H2", "ep"), ("2p24GeV_H2", "ep"), ("3p49GeV_H2", "ep"), ("2p24GeV_Ta_X17", "ep"), ("2p24GeV_Ta_X17", "moller")]:
        for c in ["V3", "V4", "C"]:
            x = f[(f.run_type == rt) & (f.calib == cal) & (f.comparison == f"{c}-A")]
            if x.dg.abs().max() == 0:
                L.append(f"  {RTL[rt] + ' ' + cal:<30}{c:<6}{len(x):>11}{'0':>12}{'0':>10}")
                continue
            m = x.loc[x.zg.abs().idxmax()]
            L.append(f"  {RTL[rt] + ' ' + cal:<30}{c:<6}{len(x):>11}{f"{100 * x.dg.abs().max():.2f} +- {100 * x.loc[x.dg.abs().idxmax()].sigma_g:.2f}":>17}%{abs(m.zg):>10.2f}  {m.module} ({int(m.run)}), dg = {100 * m.dg / m.g_A:+.2f} +- {100 * m.sigma_g / m.g_A:.2f} %, zg {m.zg:+.2f}")
    mx = f.loc[f.zg.abs().idxmax()]
    L.append(f"  Largest |zg| overall: {mx.module} ({RTL[mx.run_type]} {mx.calib}, {mx.comparison}, run {int(mx.run)}): dg/g = {100 * mx.dg / mx.g_A:+.2f} +- {100 * mx.sigma_g / mx.g_A:.2f} %, zg {mx.zg:+.2f}.  No free gain moves by 2 sigma_g.")
    hb = g[(g.run_type == "2p24GeV_Ta_X17") & g.module.isin(["W493", "W494", "W559", "W564", "W593", "W663"])]
    hf = hb.groupby("calib").agg(n=("free", "size"), nfree=("free", "sum"), maxdg=("dg", lambda s: s.abs().max()))
    L.append(f"  Ta inner-ring (block H) modules: free in the ep calibration for {int(hf.loc['ep'].nfree)} of {int(hf.loc['ep'].n)} (module, run, arm) entries and in the Moller calibration for {int(hf.loc['moller'].nfree)} of {int(hf.loc['moller'].n)}; "
             f"prior fraction {hb[hb.calib == 'moller'].prior_frac.min():.3f}-{hb[hb.calib == 'moller'].prior_frac.max():.3f} (Moller: prior-held), {hb[hb.calib == 'ep'].prior_frac.min():.2f}-{hb[hb.calib == 'ep'].prior_frac.max():.2f} (ep: about half prior); max |dg| over them {100 * hb.dg.abs().max():.3f} +- {100 * hb.loc[hb.dg.abs().idxmax()].sigma_g:.3f} % (sigma_g of that entry), so their dg is not measurable here.")
    vt = pd.read_csv(HERE / P10T, sep="\t")
    vt = vt[(vt.table == "veto") & (vt.population == "elastic") & (vt.veto == "V3") & (vt.axis == "overall")].set_index("dataset")
    tc = _one(pa, run_type="2p24GeV_Ta_X17", calib="ep", comparison="C-A", protocol="oos", axis="overall", bin="all", quantity="mu_ep")
    L += ["", "GATE DECISION (plan Log 2026-11-18 gate (P9/P10); the P9 tables above show no arm moves sigma/E, the peak or the tails beyond 2 sigma of the paired error)",
          "  2.24 and 3.49 GeV H2   keep V3 (user policy; fakes only; loss " + f"{100 * vt.loc['2p24GeV_H2'].frac:.4f} % [{100 * vt.loc['2p24GeV_H2'].lo:.4f}, {100 * vt.loc['2p24GeV_H2'].hi:.4f}] / {100 * vt.loc['3p49GeV_H2'].frac:.4f}" + " % [" + f"{100 * vt.loc['3p49GeV_H2'].lo:.4f}, {100 * vt.loc['3p49GeV_H2'].hi:.4f}" + "] of the elastic window (Wilson 68 %, table P10 veto), below 0.01 %, no efficiency correction); no V4; no in-place correction C.",
          "  0.73 GeV H2 and Ta X17 no veto, no correction (V3 is a no-op there; V4 and C do not move any measured quantity).",
          "  Flags stay event markers.  f5 is kept for the ifarm full scan.  The C scale shift is absorbed by the gains (Ta ep peak dmu(C - A) = " + f"{tc['diff']:+.4f} +- {tc.se:.4f} %" + "), so no correction is applied."]
    P.text("P9 calibration impact and gate decision", L, size=7.0, sub="Three arms (veto V3, veto V4, corrected C) against the default A; numbers read from " + P9T)


P12T = "exotic_P12_rates_vs_energy_20261118.tsv"


def page_p12(P):
    d = read(P12T)
    sl = ["100-300", "500-1000", ">4000"]
    cats = ["any", "f1 tail-baseline", "f2 wide", "f3 smooth-wide", "f4 clipped-multi", "f5 saturated single"]
    short = ["any", "f1", "f2", "f3", "f4", "f5"]
    L = ["Rate of each exclusive category (priority f4 > f5 > f1 > f2 > f3, uniform thresholds, run 24937 excluded) per energy slice, per run type; cell = rate % [Wilson 68 %] (k); N in the row label; '*' = small-N slice (TSV small_n = 1).",
         "Selected slices 100-300, 500-1000 and > 4000 MeV (all slices 10-50 ... > 4000 MeV, with the same columns, in the P12 TSV and in the PDF table pages after the section page).", ""]
    for tb, title in (("chan", "1. Per channel: kept PbWO4 channels (A_main >= 10 ADC), energy = A_main x module gain (MeV)"), ("cluster", "2. Per cluster: clusters with >= 1 flagged member (role any), slice = E_cluster (MeV)")):
        L.append(title)
        L.append(f"  {'run type, slice (MeV)':<40}" + "".join(f"{s:>38}" for s in short))
        for rt in RTS:
            for s in sl:
                x = d[(d.table == tb) & (d.run_type == rt) & (d.slice == s) & ((d.role == "any") | d.role.isna())]
                if x.empty:
                    continue
                n = int(x.n.iloc[0])
                cells = []
                for c in cats:
                    r = x[x.category == c].iloc[0]
                    cells.append(f"{100 * r.rate:.3g} % [{(0 if (int(r.k) == 0 or abs(r.rate_lo) < 1e-12) else 100 * r.rate_lo):.3g}, {100 * r.rate_hi:.3g}] ({int(r.k)})")
                star = "*" if int(x.small_n.max()) else " "
                L.append(f"  {LAB[rt] + ', ' + s + star + f'  N = {n:,}':<40}" + "".join(f"{c:>38}" for c in cells))
        L.append("")
    L += ["Caveats: the > 4000 MeV slice is populated by clipped pulses (f4, f5) whose fitted amplitude is an extrapolation above full scale (A > 4095 ADC), so its energy is not a measured energy; the slice rates for f4 / f5 there are rates of saturated channels,",
          "not of an energy range.  Slices with N < the small-N limit are marked '*' and their intervals are wide.  H2 per-channel rows come from the 2-file sample, Ta X17 from the 19-file full scan (as P5).",
          "The impact per slice (energy removed by dropping flagged channels, <fE> +/- SEM) is in the P12 PDF (impact page) and the TSV (table impact).", "",
          "Source: " + P12T + " (tables chan, cluster; columns run_type, slice, category, role, n, k, rate, rate_lo, rate_hi, small_n)."]
    P.text("P12 rates vs energy slice: per channel and per cluster", L, size=6.0, sub="Per flag, per energy slice, per run type; Wilson 68 % intervals; numbers read from " + P12T)


P13T = "exotic_P13_drop_flagged_20261118.tsv"


def p13_load():
    return pd.read_csv(HERE / P13T, sep="\t", low_memory=False)


def p13_numbers(d):
    """per run type: paired D-A oos diffs, free gains > 2 sigma_g, elastic-window cost (gross, net), Moller cost, S, fakes, lost / new clusters; all from the P13 TSV."""
    R = {}
    pr = d[(d.table == "paired") & (d.comparison == "D-A") & (d.protocol == "oos")]
    gn = d[(d.table == "gain") & (d.comparison == "D-A")]
    el = d[(d.table == "p10_elastic_loss") & (d.axis == "overall") & (d["sample"] == "elastic")].set_index("run_type")
    sc = d[d.table == "p8_scale"]
    p7 = d[(d.table == "p7_cluster") & (d.axis == "overall")]
    ls = d[d.table == "lost"]
    for rt in RTS:
        r = {}
        for q, key in (("sig", "res_cb"), ("mu", "mu_ep"), ("tail", "tail_lo3")):
            x = pr[(pr.run_type == rt) & (pr.calib == "ep") & (pr.quantity == key) & (pr["bin"] == "all")].iloc[0]
            r[q] = x
        g = gn[(gn.run_type == rt) & (gn.calib == "ep")]
        g = g[(g.free > 0) | (g.free_X > 0)]
        _b = d[(d.table == "bench") & (d.arm == "A") & (d.run_type == rt) & (d.quantity == "res_cb_oos")].iloc[0]
        r["A_val"], r["A_se"] = _b.value, _b.se
        r["gain_ep"] = (int((g.dg.abs() > 2 * g.sigma_g).sum()), len(g))
        e = el.loc[rt]
        r["gross"] = (int(e.n_lost), int(e.n_A), 100 * e.lost_lo, 100 * e.lost_hi)
        ne = int(e.n_enter)
        r["net"] = (100 * (e.n_lost - ne) / e.n_A, 100 * np.sqrt(e.n_lost + ne) / e.n_A, ne)
        r["S"] = sc[sc.run_type == rt].set_index("sample")
        r["p7"] = p7[p7.run_type == rt].set_index("sample")
        l = ls[(ls.run_type == rt) & (ls.protocol == "oos")].set_index("quantity")
        r["lost"] = l
        R[rt] = r
    ta = "2p24GeV_Ta_X17"
    R["moller"] = {q: pr[(pr.run_type == ta) & (pr.calib == "moller") & (pr["bin"] == "sum") & (pr.quantity == k)].iloc[0] for q, k in (("sig", "res"), ("mu", "peak"), ("tail", "tail_lo3"))}
    g = gn[(gn.run_type == ta) & (gn.calib == "moller")]
    g = g[(g.free > 0) | (g.free_X > 0)]
    R["moller_A_se"] = d[(d.table == "bench") & (d.arm == "A") & (d.run_type == ta) & (d.quantity == "sum:res_pooled")].se.iloc[-1]
    R["moller_gain"] = (int((g.dg.abs() > 2 * g.sigma_g).sum()), len(g))
    ym = d[(d.table == "yield_moller") & (d.axis == "overall")].set_index("quantity")
    R["moller_cost"] = ym
    return R


def P13_DECISION():
    try:
        R = p13_numbers(p13_load())
    except Exception as ex:                                              # keep the summary buildable if the P13 TSV is absent
        return [f"  2026-11-18  P13 outcome: P13 TSV not readable ({ex})."]
    h1, h2 = R["2p24GeV_H2"], R["3p49GeV_H2"]
    return [
        "  2026-11-18  P13 OUTCOME (arm D, drop of f1-f5 channels, vs the default A; ep elastic window; page 16): DROP POLICY NOT ADOPTED, V3 stays on 2.24 / 3.49 GeV H2 only, flags stay event markers.",
        f"              No measurable gain in sigma/E or the peak anywhere; window cost D vs V3 at 2.24 / 3.49 GeV H2: {ci1(h1['gross'][0], h1['gross'][1], 4)} / {ci1(h2['gross'][0], h2['gross'][1], 4)} vs {_v3('2p24GeV_H2')} / {_v3('3p49GeV_H2')};",
        f"              Ta X17: low tail widens ({R['2p24GeV_Ta_X17']['tail']['diff']:+.2f} +- {R['2p24GeV_Ta_X17']['tail']['se']:.2f} pp, z {R['2p24GeV_Ta_X17']['tail']['z']:+.1f}), elastic-window loss net {R['2p24GeV_Ta_X17']['net'][0]:.2f} +- {R['2p24GeV_Ta_X17']['net'][1]:.2f} % (Poisson), Moller pairs {ci1(int(R['moller_cost'].loc['lost'].k), int(R['moller_cost'].loc['lost'].n), 3)} lost.",
    ]


def _v3(ds):
    v = read(P10T)
    x = v[(v.table == "veto") & (v.population == "elastic") & (v.veto == "V3") & (v.axis == "overall") & (v.dataset == ds)].iloc[0]
    return f"{100 * x.frac:.4f} % [{100 * x.lo:.4f}, {100 * x.hi:.4f}]"


def page_p13(P):
    d = p13_load()
    R = p13_numbers(d)
    ta = "2p24GeV_Ta_X17"
    pm = lambda x, nd=3: f"{x['diff']:+.{nd}f} +- {x['se']:.{nd}f} (z {x['z']:+.1f})"
    L = ["Policy (user, 2026-11-18): a kept PbWO4 channel in any exclusive category f1-f5 is not used: its hit is removed before clustering and the gain solve (arm D); D' = drop in the gain solve only.  A = default control v3.",
         "Scope as P7-P10: H2 on the 2-file sample (4 / 10 / 4 splits), Ta X17 on the 19 full-scan splits, run 24937 excluded, 0.73 GeV at seed radius >= 107 mm.  Paired = clusters present in both arms (same run, split, event, seed), bootstrap error.", "",
         "1. D - A, out of sample, paired: sigma/E (CB, percentage points), peak mu (%), low 3-sigma tail (percentage points); free module gains moved by more than 2 sigma_g (of the free module-runs in A or D)",
         f"  {'calibration':<22}{'A sigma/E oos (%)':>19}{'d sigma/E':>31}{'d mu':>31}{'d low tail':>31}{'gains > 2 sigma_g':>20}"]
    for rt in RTS:
        r = R[rt]
        L.append(f"  {LAB[rt] + ' ep':<22}{f"{r['A_val']:.3f} +- {r['A_se']:.3f}":>19}{pm(r['sig'], 4):>31}{pm(r['mu'], 4):>31}{pm(r['tail'], 3):>31}{f'{r['gain_ep'][0]} of {r['gain_ep'][1]}':>20}")
    m = R["moller"]
    L.append(f"  {'2.24 GeV Ta X17 Moller':<22}{f"{m['sig'].value_A:.3f} +- {R['moller_A_se']:.3f}":>19}{pm(m['sig'], 4):>31}{pm(m['mu'], 4):>31}{pm(m['tail'], 3):>31}{f'{R['moller_gain'][0]} of {R['moller_gain'][1]}':>20}")
    L += ["  (A sigma/E oos: the benchmark's own value +- SE over the splits, arm A; Moller row: pair sum, pooled fit.)", "",
          "2. Cost: ep elastic window of A (fixed gains) lost under D, Wilson 68 %; net of the clusters entering the window (Poisson error); Ta Moller pairs lost / new (Wilson 68 %); V3 (P10) for comparison; V4 (P10) in the Reading below",
          f"  {'run type':<22}{'lost, gross % [lo, hi] (k of N)':<50}{'net %':<22}{'V3 (P10)':<34}"]
    for rt in RTS:
        r = R[rt]
        k, n = r["gross"][:2]
        L.append(f"  {LAB[rt]:<22}{ci1(k, n, 4) + f' ({k} of {n:,})':<50}{f'{r['net'][0]:.4f} +- {r['net'][1]:.4f}':<22}{_v3(rt):<34}")
    mc = R["moller_cost"]
    kl, nl = int(mc.loc["lost"].k), int(mc.loc["lost"].n)
    kn, nn = int(mc.loc["new_in_D"].k), int(mc.loc["new_in_D"].n)
    L.append(f"  {'Ta X17 Moller pairs':<22}{ci1(kl, nl, 3) + f' ({kl} of {nl:,}) lost':<50}{ci1(kn, nn, 3) + f' ({kn} of {nn:,}) of D pairs new':<56}")
    L += ["", "3. Relative energy change S = dE/E of the clusters under D (first order, flagged members' energy removed; negative = energy removed), per sample, per mille +- error; right: clusters of the P9 ep calibration sample (oos, not paired) of A lost under D and new in D, fraction [Wilson 68 %] (k of N)",
          f"  {'run type':<18}{'S all (per mille)':>20}{'S hiE':>18}{'S elastic':>18}{'  lost clusters':<44}{'new clusters in D':<44}"]
    for rt in RTS:
        r = R[rt]
        S = r["S"]
        cells = [f"{1000 * S.loc[s].S1:+.3f} +- {1000 * S.loc[s].S1_err:.3f}" for s in ("all", "hiE", "elastic")]
        lo = r["lost"].loc["lost"]
        nw = r["lost"].loc["new_in_D"]
        L.append(f"  {LAB[rt]:<18}{cells[0]:>20}{cells[1]:>18}{cells[2]:>18}  {ci1(int(lo.k), int(lo.n), 3) + f' ({int(lo.k)} of {int(lo.n):,})':<44}{ci1(int(nw.k), int(nw.n), 3) + f' ({int(nw.k)} of {int(nw.n):,})':<44}")
    S = R[ta]["S"]
    L.append(f"  Ta X17 Moller sample: S = {1000 * S.loc['moller'].S1:+.3f} +- {1000 * S.loc['moller'].S1_err:.3f} per mille.")
    L += ["", "4. V3 fakes (flagged saturated seed, P7 definition) among the clusters of A removed by D: fakes removed, k of N fakes [Wilson 68 %], whole sample and elastic window"]
    for rt in RTS:
        p = R[rt]["p7"]
        a, e = p.loc["all"], p.loc["elastic"]
        ka, na, ke = int(a.n_fake_removed), int(a.n_fake), int(e.n_fake_removed)
        L.append(f"  {LAB[rt]:<18}all: {ka} of {na} fakes removed, " + (ci1(ka, na, 1) if na else "no fakes") + f"; elastic window: {ke} of {int(e.n_fake)} fakes removed, {int(e.n):,} clusters, " + ci1(ke, int(e.n), 4) + " of the elastic clusters")
    h1, h2, t = R["2p24GeV_H2"], R["3p49GeV_H2"], R[ta]
    L += ["", "Recommendation (from the P13 gate page): keep the current gate, V3 on 2.24 / 3.49 GeV H2 only, no correction, flags stay event markers.  Channel-level DROP POLICY NOT ADOPTED:",
          f"  no measurable resolution or peak improvement anywhere; on H2 it removes the fakes V3 removes but also real clusters (window loss {ci1(h1['gross'][0], h1['gross'][1], 4)} at 2.24, {ci1(h2['gross'][0], h2['gross'][1], 4)} at 3.49 GeV, V3 {_v3('2p24GeV_H2')} / {_v3('3p49GeV_H2')});",
          f"  at 0.73 GeV window loss {ci1(R['0p73GeV_H2']['gross'][0], R['0p73GeV_H2']['gross'][1], 4)} and {R['0p73GeV_H2']['gain_ep'][0]} gains beyond 2 sigma_g; on Ta the low tail widens ({pm(t['tail'], 3)} pp), the Moller peak shifts ({pm(m['mu'], 3)} %),",
          f"  the elastic window loses {t['net'][0]:.3f} +- {t['net'][1]:.3f} % net and the Moller pairs {ci1(kl, nl, 3)}.", "",
          "Source: " + P13T + " (tables paired, gain, p10_elastic_loss, yield_moller, p8_scale, p7_cluster, lost); the full P13 report follows this page's section page."]
    P.text("P13 drop policy (P7-P10 recomputed)", L, size=6.6, sub="Arm D (flagged channels dropped) vs A (default); out of sample, paired; Wilson 68 % intervals or +- errors; numbers read from " + P13T)


def section_page(pdf, tag, step, date, script, outcome, fname):
    import textwrap
    fig = plt.figure(figsize=(16, 9))
    fig.text(0.06, 0.78, f"{tag}", fontsize=60, fontweight="bold", va="top")
    fig.text(0.06, 0.58, fname, fontsize=16, va="top", family="monospace")
    body = [f"Date: {date}", f"Script: {script}", "", "Outcome:"] + textwrap.wrap(outcome, 110)
    fig.text(0.06, 0.52, "\n".join(body), fontsize=14, va="top", linespacing=1.5)
    pdf.savefig(fig)
    plt.close(fig)


# ----------------------------------------------------------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default=str(HERE / "exotic_summary_20261118.pdf"))
    a = ap.parse_args()
    S, M, T, P5, P6, N, sv, info = build_numbers()

    buf = io.BytesIO()
    nfront = 0
    with PdfPages(buf) as pdf:
        P = Pages(pdf)
        page_definitions(P)
        # page 2
        L = page_rates(P, S, N, sv)
        L = page_sat_rows(L, S, N, sv, P6)
        L += ["", "Sources: A, B, C, D (full scan): " + P4S + " (threshold_set uniform: excl_*, raw flag, moved_*, clipped); D E>0.5 E_beam ratio: " + P6T + " section saturation."]
        P.text("Latest numbers 1: rates, saturation, moved matrix", L, size=7.6, sub="Run types: 0.73 / 2.24 / 3.49 GeV H2 and 2.24 GeV Ta X17; PbWO4 only; run 24937 excluded; full-scan N kept in the table.")
        L = page_thresholds(P, S, T, N, sv)
        L += ["", "Sources: " + P4T + " (thresholds, core stats); " + P4S + " (threshold_set type / uniform / check)."]
        P.text("Latest numbers 2: thresholds, per-type rates, f2 exclusions", L, size=8.0)
        L = page_top(P, M)
        L += ["Source: " + P4M + " (basis column: raw / excl; stat_ok column)."]
        P.text("Latest numbers 3: top modules per category", L, size=7.4)
        L, tr = page_radial(P, M)
        L += ["Sources: " + P4M + " (summed per radius bin with the HyCal module map); saturation rows by radius also in " + P4S + " (threshold_set radius_*).", ""]
        L.append("Sentences (checked against the table above):")
        L += RADIAL_SENT(tr)
        P.text("Latest numbers 4: radial structure", L, size=7.0)
        L = page_amp(P, P5)
        L += ["", "Source: " + P5T + " (rt, flag, defn, N, median, q16, q84, source columns)."]
        P.text("Latest numbers 5: maximum amplitude per category (P5)", L, size=7.6)
        L = page_f2(P, P6)
        L += F2_SENT(P6)
        L += ["", "Source: " + P6T + " (sections f2_high_group, w50_profile)."]
        P.text("Latest numbers 6: f2 diagnosis (P6)", L, size=7.2)
        page_decisions(P)
        page_p7(P)
        Ssat = pd.read_csv(HERE / P7S, sep="\t")
        page_sat_table(P, Ssat)
        page_sat_pops(P, Ssat)
        page_p8(P)
        page_p10(P)
        page_p9(P)
        page_p12(P)
        page_p13(P)
        nfront = 16
        sec = io.BytesIO()
    front = PdfReader(io.BytesIO(buf.getvalue()))
    W = PdfWriter()
    for pg in front.pages:
        W.add_page(pg)
    W.add_outline_item("Front matter: definitions and latest numbers", 0)
    ORDER = ["P1", "P2", "P3", "P4", "P5", "P6", "P7", "P8", "P10", "P9", "P12", "P13"]
    oc = dyn_outcomes()
    for tag, fn, date, script, outcome in sorted([(t, f, d, s, oc.get(t, o)) for t, f, d, s, o in STEPS], key=lambda x: ORDER.index(x[0])):
        sb = io.BytesIO()
        with PdfPages(sb) as pdf:
            section_page(pdf, tag, tag, date, script, outcome, fn)
        start = len(W.pages)
        W.add_page(PdfReader(io.BytesIO(sb.getvalue())).pages[0])
        for pg in PdfReader(str(HERE / fn)).pages:
            W.add_page(pg)
        W.add_outline_item(f"{tag}  {fn}", start)
        print(f"{tag}: section page + {len(PdfReader(str(HERE / fn)).pages)} pages from {fn}")
    with open(a.out, "wb") as fh:
        W.write(fh)
    print(f"front matter pages: {len(front.pages)}; total pages: {len(W.pages)}; wrote {a.out}")
    for c in info["checks"]:
        print("INCONSISTENT:", c)
    if not info["checks"]:
        print("internal checks: raw - moved = exclusive, sum of exclusive = any, for all run types: OK")


def RADIAL_SENT(tr):
    def kn(c, rt, i):
        n, f = tr[(c, rt)][2][i]
        return f, n

    def r(c, rt, i):
        return ci1(*kn(c, rt, i))

    def rs(cs, rt, i):                       # summed categories
        return ci1(sum(kn(c, rt, i)[0] for c in cs), kn(cs[0], rt, i)[1])
    ta, h = "2p24GeV_Ta_X17", "2p24GeV_H2"
    L = []
    imin = min(range(1, 5), key=lambda i: tr[("f3", ta)][0][i])
    L += wrap(f"Ta X17: f1 and f2 concentrate in the inner ring around the beam hole and fall with radius (f1 {r('f1', ta, 0)} at 0-120 mm to {r('f1', ta, 4)} at 320-450 mm; "
              f"f2 {r('f2', ta, 0)} to {r('f2', ta, 4)}); f3 is U-shaped ({r('f3', ta, 0)} inner, minimum {r('f3', ta, imin)} at {RADBINS[imin]}-{RADBINS[imin + 1]} mm, {r('f3', ta, 4)} at 320-450 mm).")
    f1k, f1n = kn("f1", "3p49GeV_H2", 4)
    L += wrap(f"H2: f3 rises toward the outer edge at all three energies, and f2, f4, f5 at 2.24 and 3.49 GeV (0.73 GeV: f4/f5 peak at 230-320 mm on few counts, f2 at 230-320 and 450-700 mm) (2.24 GeV: f3 {r('f3', h, 0)} to {r('f3', h, 4)}, f4 {r('f4', h, 0)} to {r('f4', h, 4)}, "
              f"f5 {r('f5', h, 0)} to {r('f5', h, 4)}, 0-120 to 320-450 mm); f1 falls with radius at 0.73 and 2.24 GeV but rises in the outer bins at 3.49 GeV ({ci1(f1k, f1n)} at 320-450 mm, {f1k} of {f1n:,} flags).")
    top = max(range(5), key=lambda i: tr[("f4", ta)][2][i][1] and (kn("f4", ta, i)[0] + kn("f5", ta, i)[0]) / tr[("f4", ta)][2][i][0] if tr[("f4", ta)][2][i][0] else -1)
    L += wrap(f"Saturation (f4 + f5) rises with radius and is not an inner-ring effect: at 2.24 GeV H2 it is {rs(['f4', 'f5'], h, 0)} at 0-120 mm and {rs(['f4', 'f5'], h, 4)} at 320-450 mm; the Ta X17 saturated fraction is at most {rs(['f4', 'f5'], ta, top)} (bin {RADBINS[top]}-{RADBINS[top + 1]} mm, the largest) in any bin. "
              "Elastic-edge amplitudes per bin are on the P4 saturation page.")
    return L

def F2_SENT(P6):
    g = P6[P6.section == "f2_high_group"].set_index(["run_type", "source"])
    s = lambda rt, src: g.loc[(rt, src)]
    a, b, ta = s("0p73GeV_H2", "sample (2 files/run)"), s("2p24GeV_H2", "sample (2 files/run)"), s("2p24GeV_Ta_X17", "full scan npz (19 files)")
    pc = lambda x, lo, hi, nd=0: f"{g3(100 * x)} % [{g3(100 * lo)}, {g3(100 * hi)}]"
    L = wrap(f"2.24 GeV H2: the high-A f2 mode ({int(b.n_f2_A_gt300)} of {int(b.n_f2)}, {pc(b.frac_A_gt300, b.frac_lo68, b.frac_hi68)}) is merged k >= 4 pile-up ({int(b.k4p_hi)} of {int(b.n_f2_A_gt300)} with k >= 4, {pc(b.frac_k4p_hi, b.frac_k4p_hi_lo68, b.frac_k4p_hi_hi68)}) spread over {int(b.n_modules_hi)} modules "
             f"(top-module share {pc(b.top1_share_hi, b.top1_share_hi_lo68, b.top1_share_hi_hi68)}); it is not a single bad module.")
    c3 = s('3p49GeV_H2', 'sample (2 files/run)')
    L += wrap(f"0.73 GeV H2: the high-A mode ({int(a.n_f2_A_gt300)} of {int(a.n_f2)}, {pc(a.frac_A_gt300, a.frac_lo68, a.frac_hi68)}) is mostly module W335 ({a.top_modules_hi.split(';')[0]} of the high-A channels, share {pc(a.top1_share_hi, a.top1_share_hi_lo68, a.top1_share_hi_hi68)}), also k >= 4 ({int(a.k4p_hi)} of {int(a.n_f2_A_gt300)}, {pc(a.frac_k4p_hi, a.frac_k4p_hi_lo68, a.frac_k4p_hi_hi68)}); "
              f"3.49 GeV H2 has {int(c3.n_f2_A_gt300)} of {int(c3.n_f2)} ({pc(c3.frac_A_gt300, c3.frac_lo68, c3.frac_hi68)}, small).")
    wp = P6[P6.section == "w50_profile"]
    wt = wp[wp.run_type == "2p24GeV_Ta_X17"].sort_values("w50_lo").reset_index(drop=True)
    kb = lambda i: f"{g3(100 * wt.frac_k_ge2[i])} % [{g3(100 * wt.frac_k_ge2_lo68[i])}, {g3(100 * wt.frac_k_ge2_hi68[i])}] ({int(wt.n_k_ge2[i])} of {int(wt.n_kept[i])})"
    L += wrap(f"Ta X17: f2 is not a distinct population: only {int(ta.n_f2_A_gt300)} of {int(ta.n_f2)} are high-A ({pc(ta.frac_A_gt300, ta.frac_lo68, ta.frac_hi68, 1)}), the median A is {ta.A_median:.1f} ADC, and {int(ta.k4p_lo)} of {int(ta.n_f2 - ta.n_f2_A_gt300)} low-A flags have k >= 4 "
              f"({pc(ta.frac_k4p_lo, ta.frac_k4p_lo_lo68, ta.frac_k4p_lo_hi68)}; baseline-wobble over-fit); the w50 profile shows a continuum with the same k >= 2 fraction below and above w50 = 3 ({kb(2)} at w50 2-3, {kb(3)} at w50 3-4, {kb(4)} above 4), so w50 = 3 is a position in the tail of that continuum, not a population boundary.")
    fac, ferr = {}, {}
    for rt in RTS:
        n = wp[wp.run_type == rt].sort_values("w50_lo").n_kept.tolist()
        fac[rt] = n[2] / n[3]
        ferr[rt] = fac[rt] * np.sqrt(1 / n[2] + 1 / n[3])
    L += wrap(f"Counts per w50 bin fall from the 2-3 bin to the 3-4 bin by a factor " + " / ".join(f"{fac[rt]:.3g} +/- {ferr[rt]:.2g}" for rt in RTS[:3]) + f" in H2 (0.73 / 2.24 / 3.49 GeV) but only {fac['2p24GeV_Ta_X17']:.3g} +/- {ferr['2p24GeV_Ta_X17']:.2g} in Ta X17 (Poisson errors; table above): "
              "w50 = 3 sits in a density gap for H2 and is an arbitrary cut in a continuum for Ta.")
    return L

if __name__ == "__main__":
    main()
