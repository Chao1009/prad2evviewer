#!/usr/bin/env python3
"""benchmark_pileup_mf.py -- matched-filter detection efficiency on synthetic pile-up cells.

Reads the cells written by pileup_generator.py (manifest.json + pileup_*_dt###_
ratio####.npz; col 0 = injected pulse, col 1 = base pulse), runs the matched
filter (mf_pulse.MFReconstructor on a bank built with mf_pulse.build_bank, a
two-tau h50 table or a sampled .npz) on every synthetic waveform and reports,
per (dt, ratio) cell and per secondary-acceptance r_model:

  eff_match  fraction of events with an mf pulse within +-tol of the injected
             onset AND a different mf pulse within +-tol of the base onset
  fake_rate  mean number of mf pulses per event within +-tol of neither onset
  n2_mf      fraction with >= 2 mf pulses (Table 4's definition on the mf path)
  n2_wa      fraction with >= 2 WaveAnalyzer peaks (len(wave_ana.analyze(wf)[2])
             on the raw uint16 waveform; Table 4's exact column, independent of
             the template and of r_model)
  eff_match_clean, fake_rate_clean
             the same two on the events whose base onset is >= ped_nsamples/2
             samples (60 ns for 30), i.e. whose pulse cannot bias the mf
             pedestal (median of the first ped_nsamples samples; the same rule
             plot_pileup_examples.py applies).  The generator keeps bases from
             t0 >= 22 ns, so in run 025308 about 40 % of the bases start inside
             that window; the mf then mis-subtracts the baseline and loses the
             injected pulse.  The four main columns use ALL events (Table 4
             compatible); n_clean is stored per cell.

Onset convention.  mf p_t0 is the template onset in ns (the instant the unit
template leaves zero; the template peak is at bank.t_pk[m] later).  The
generator's truth_onset_t0_ns is the same quantity: col 0 = peak_inj - t_peak_
offset of the template it injected, col 1 = onset of the C++ single-pulse fit of
the base pulse.  Both are compared to p_t0 directly.  The generator's
truth_peak_times_ns is NOT used for matching (col 1 holds WaveAnalyzer pk.time,
which sits within ~1 ns of the fit onset, not one t_pk after it).  The residual
p_t0 - truth onset is measured on the reference cell (ratio 0.1, dt 128, else
the sparsest cell) and written to the JSON and the .md.

Zero-point.  A sampled bank on cells that carry truth_onset_sampled_ns is
compared with that column (same onset definition as the bank); a two-tau bank, or
a sampled bank on two-tau cells, uses truth_onset_t0_ns, where --onset-shift
auto / A,B removes the constant offset measured on the reference cell (default
none; it warns above half the tolerance).

Window rule (--window-rule DT_MAX:RATIO, default off; production 40:0.16, the
rule of plot_pileup_examples.py).  Applied post hoc to the accepted pulses of
every event: a pulse that is not the event's largest (E = pulse_energy_units,
kappa 1) with 0 <= t_sec - t_dom < DT_MAX ns (t = p_time) and E_sec / E_dom <
RATIO is dropped before matching and counting.  The grids of the plain run stay
as they are; the grids with the rule are an additional block per r_model
("results_window" in the JSON, "with window rule" tables in the .md, n2_wa is
unaffected), plus a check of the eff_match change in the cells with injected
ratio >= 0.2 and dt >= 24 ns.

    MPLBACKEND=Agg python3 benchmark_pileup_mf.py output/pileup_synth_025308 \\
        --templates output/pulse_templates_h50.dat --out <dir> --table4-compare
"""

from __future__ import annotations

import argparse
import glob
import hashlib
import json
import math
import os
import re
import sys
import time
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import matplotlib
matplotlib.use(os.environ.get("MPLBACKEND", "Agg"))
import matplotlib.pyplot as plt                       # noqa: E402

import _common as C                                   # noqa: E402
import mf_pulse                                       # noqa: E402

HERE = Path(__file__).resolve().parent

# Table 4, "WaveAnalyzer detection efficiency" (% of events with >= 2 findPeaks
# peaks; run 025308, 50 events per cell).  Source:
# docs/technical_notes/waveform_analysis/fit_pulse_template.md lines 495-509
# (section 6, "Pile-Up Deconvolution Benchmark").  The doc prints only the
# ratios 0.1, 0.3, 0.5, 1.0, 2.0 (no 0.2 column) and the dt rows 8, 16, 20, 24,
# 32, 48, 64, 96, 128 (no 12 ns row): 45 of the 60 cells.
TABLE4_RATIOS = (0.1, 0.3, 0.5, 1.0, 2.0)
TABLE4 = {
    8:   (0, 0, 0, 0, 0),
    16:  (0, 0, 0, 0, 0),
    20:  (0, 6, 2, 16, 14),
    24:  (0, 4, 8, 100, 100),
    32:  (0, 8, 8, 100, 100),
    48:  (8, 8, 100, 100, 100),
    64:  (8, 6, 100, 100, 100),
    96:  (6, 8, 100, 100, 100),
    128: (16, 18, 100, 100, 100),
}

METRICS = ("eff_match", "fake_rate", "n2_mf", "n2_wa",
           "eff_match_clean", "fake_rate_clean")
TITLES = {"eff_match": "eff_match (% events, both pulses matched)",
          "fake_rate": "fake_rate (unmatched mf pulses / event)",
          "n2_mf": "n2_mf (% events with >= 2 mf pulses)",
          "n2_wa": "n2_wa (% events with >= 2 WaveAnalyzer peaks)",
          "eff_match_clean": "eff_match_clean (% events, base onset outside the mf pedestal window)",
          "fake_rate_clean": "fake_rate_clean (unmatched mf pulses / event, same subset)"}


def _wilson(k, n, z=1.96):
    if n <= 0:
        return (0.0, 1.0)
    ph = k / n
    den = 1.0 + z * z / n
    ctr = (ph + z * z / (2 * n)) / den
    hw = z * math.sqrt(ph * (1 - ph) / n + z * z / (4 * n * n)) / den
    return (max(0.0, ctr - hw), min(1.0, ctr + hw))


def _file_sha1(path):
    with open(path, "rb") as fh:
        return hashlib.sha1(fh.read()).hexdigest()


# ------------------------------------------------------------------ cells --
def find_cells(synth_dir):
    """[(dt_ns, ratio, npz_path)] from manifest.json, else a glob of the dir."""
    synth_dir = Path(synth_dir)
    man = {}
    mp = synth_dir / "manifest.json"
    if mp.exists():
        man = json.loads(mp.read_text())
    cells = []
    for f in man.get("files", []):
        pth = synth_dir / f["path"]
        if pth.exists():
            cells.append((float(f["dt_ns"]), float(f["ratio"]), pth))
    if not cells:
        rx = re.compile(r"_dt(\d+)_ratio(\d+)\.npz$")
        for pth in sorted(glob.glob(str(synth_dir / "pileup_*_dt*_ratio*.npz"))):
            m = rx.search(pth)
            if m:
                cells.append((float(m.group(1)), int(m.group(2)) / 1000.0,
                              Path(pth)))
    return sorted(cells, key=lambda c: (c[0], c[1])), man


def _manifest_cli_value(man, flag):
    m = re.search(re.escape(flag) + r"\s+(\S+)", str(man.get("cli", "")))
    return m.group(1) if m and os.path.exists(m.group(1)) else ""


# --------------------------------------------------------------- matching --
def match_events(res, n_ev, t_inj, t_base, tol):
    """Per-event match flags from an MFResult (pulses sorted by channel)."""
    pc = res.p_ch
    pt = res.p_t0
    near_i = np.abs(pt - t_inj[pc]) <= tol
    near_b = np.abs(pt - t_base[pc]) <= tol
    n_a = np.bincount(pc, weights=near_i, minlength=n_ev)
    n_b = np.bincount(pc, weights=near_b, minlength=n_ev)
    n_u = np.bincount(pc, weights=(near_i | near_b), minlength=n_ev)
    # a distinct pair (one pulse near each truth) exists iff both sets are
    # non-empty and their union holds >= 2 pulses
    match = (n_a >= 1) & (n_b >= 1) & (n_u >= 2)
    n_fake = np.bincount(pc, weights=~(near_i | near_b), minlength=n_ev)
    return match, n_fake


def parse_window_rule(txt):
    """'DT_MAX:RATIO' -> (dt_max_ns, ratio) or None."""
    if not txt:
        return None
    try:
        d, r = txt.split(":")
        return float(d), float(r)
    except ValueError:
        sys.exit(f"[ERROR] --window-rule expects DT_MAX:RATIO (e.g. 40:0.17), got {txt!r}")


def window_keep(res, bank, midx, n_ev, rule):
    """Per-pulse keep mask: False for a pulse that is not its event's largest
    (pulse_energy_units, kappa 1) with 0 <= t - t_dom < DT_MAX and E/E_dom < RATIO."""
    pc = res.p_ch
    if pc.size == 0:
        return np.ones(0, dtype=bool)
    dmax, rmax = rule
    E = mf_pulse.pulse_energy_units(bank, midx[pc], res.p_t0, res.p_amp, 1.0)
    emax = np.zeros(n_ev)
    np.maximum.at(emax, pc, E)
    # first pulse reaching the event maximum is the dominant one
    idx = np.arange(pc.size)
    isd = E >= emax[pc]
    dom = np.full(n_ev, pc.size, dtype=np.int64)
    np.minimum.at(dom, pc[isd], idx[isd])
    t_dom = res.p_time[dom[pc]]
    dt = res.p_time - t_dom
    with np.errstate(divide="ignore", invalid="ignore"):
        ratio = E / emax[pc]
    drop = (idx != dom[pc]) & (dt >= 0.0) & (dt < dmax) & (ratio < rmax)
    return ~drop


def cell_metrics(pc, pt, n_ev, t_inj, t_base, tol, clean, n_wa):
    """Metrics of one cell from the accepted pulses (channel, onset)."""
    res = SimpleNamespace(p_ch=pc, p_t0=pt)
    match, n_fake = match_events(res, n_ev, t_inj, t_base, tol)
    n_mf = np.bincount(pc, minlength=n_ev)
    nc = int(clean.sum())
    return dict(
        n_events=int(n_ev), n_match=int(match.sum()),
        eff_match=float(match.mean()), fake_rate=float(n_fake.mean()),
        n_fake_events=int((n_fake > 0).sum()),
        n2_mf=float((n_mf >= 2).mean()), n_mf2=int((n_mf >= 2).sum()),
        n2_wa=n_wa / n_ev, n_wa2=n_wa,
        mean_n_mf=float(n_mf.mean()), n_clean=nc,
        eff_match_clean=(float(match[clean].mean()) if nc else float("nan")),
        fake_rate_clean=(float(n_fake[clean].mean()) if nc else float("nan")))


def onset_residuals(res, m_of_ev, bank, d, okey, win=10.0, shift=(0.0, 0.0)):
    """p_t0 - truth onset (nearest pulse within `win` ns) for both columns,
    and the same for the base peak convention (p_t0 + t_pk - truth peak)."""
    to = d[okey].astype(np.float64) + np.asarray(shift)[None, :]
    tp = d["truth_peak_times_ns"].astype(np.float64)
    out = {}
    for col, tag in ((0, "injected"), (1, "base")):
        r_on, r_pk = [], []
        for i in range(len(m_of_ev)):
            sel = res.p_ch == i
            if not sel.any():
                continue
            t = res.p_t0[sel]
            j = int(np.argmin(np.abs(t - to[i, col])))
            if abs(t[j] - to[i, col]) <= win:
                r_on.append(t[j] - to[i, col])
                r_pk.append(t[j] + bank.t_pk[m_of_ev[i]] - tp[i, col])
        def st(x):
            x = np.asarray(x)
            if x.size == 0:
                return dict(n=0)
            mad = 1.4826 * np.median(np.abs(x - np.median(x)))
            return dict(n=int(x.size), median=float(np.median(x)),
                        sigma_mad=float(mad), mean=float(x.mean()),
                        rms=float(np.sqrt((x ** 2).mean())))
        out[tag] = dict(onset=st(r_on), peak_convention=st(r_pk))
    return out


# ------------------------------------------------------------------- main --
def build_parser():
    ap = argparse.ArgumentParser(
        description="Matched-filter efficiency / fake rate / '>= 2 peaks' "
                    "on pileup_generator.py cells.")
    ap.add_argument("synth_dir", help="pileup_generator output dir (manifest.json + cells)")
    ap.add_argument("--templates", required=True,
                    help="mf template table: spline_templates_v1 .npz or h50 .dat")
    ap.add_argument("--out", required=True, help="output directory")
    ap.add_argument("--slug", default="",
                    help="output name suffix (default: template file stem)")
    ap.add_argument("--tolerance-ns", type=float, default=4.0,
                    help="onset match tolerance (ns) [4]")
    ap.add_argument("--r-model", type=float, default=0.3,
                    help="headline mf secondary acceptance [0.3]")
    ap.add_argument("--r-model-2", type=float, default=0.10,
                    help="second table's r_model [0.10]; <= 0 disables")
    ap.add_argument("--onset-shift", default="none",
                    help="ns added to the truth onsets (injected,base) before "
                         "matching: 'none' [default], 'auto' (median reference-"
                         "cell residual per column) or 'A,B'")
    ap.add_argument("--nsigma", type=float, default=5.0)
    ap.add_argument("--min-amp", type=float, default=5.0)
    ap.add_argument("--kmax", type=int, default=4)
    ap.add_argument("--ped-nsamples", type=int, default=0,
                    help="mf pedestal window in samples [WaveConfig.ped_nsamples]")
    ap.add_argument("--column", default="auto",
                    help="h50 .dat column (ignored for .npz) [auto, prefer 25320]")
    ap.add_argument("--run-num", type=int, default=-1,
                    help="run for the pipeline setup [base_run_numbers of the first cell]")
    ap.add_argument("--daq-config", default="",
                    help="DAQ config [the one in the manifest cli, else database]")
    ap.add_argument("--hc-map-file", default="")
    ap.add_argument("--window-rule", default=None, metavar="DT_MAX:RATIO",
                    help="drop accepted secondaries with 0 <= dt < DT_MAX ns and "
                         "E/E_dom < RATIO before matching, as an additional "
                         "block of grids (default off; production 40:0.16)")
    ap.add_argument("--table4-compare", action="store_true",
                    help="gate G2-2: compare n2_wa with Table 4 of fit_pulse_template.md "
                         "and check n2_mf <= eff_match + fake_rate")
    return ap


def grid(cells_res, dts, ratios, key):
    g = np.full((len(dts), len(ratios)), np.nan)
    for (dt, r), v in cells_res.items():
        g[dts.index(dt), ratios.index(r)] = v[key]
    return g


def main(argv=None):
    a = build_parser().parse_args(argv)
    wrule = parse_window_rule(a.window_rule)
    t_start = time.monotonic()
    out = Path(a.out)
    out.mkdir(parents=True, exist_ok=True)

    cells, man = find_cells(a.synth_dir)
    if not cells:
        sys.exit(f"[ERROR] no cells found in {a.synth_dir}")
    dts = sorted({c[0] for c in cells})
    ratios = sorted({c[1] for c in cells})
    first = np.load(cells[0][2], allow_pickle=False)
    run_num = a.run_num if a.run_num > 0 else int(first["base_run_numbers"][0])
    daq = a.daq_config or _manifest_cli_value(man, "--daq-config")
    hcmap = a.hc_map_file or _manifest_cli_value(man, "--hc-map-file")

    p = C.setup_pipeline(evio_path="", run_num=run_num, daq_config=daq,
                         hc_map_file=hcmap)
    hy = [p.hycal.module(i) for i in range(p.hycal.module_count())]
    names = [m.name for m in hy]
    pwo4 = [bool(m.is_pwo4()) for m in hy]
    ix = {n: i for i, n in enumerate(names)}
    wcfg = p.wave_ana.cfg
    wave_cfg = dict(smooth_order=int(wcfg.smooth_order),
                    int_tail_ratio=float(wcfg.int_tail_ratio),
                    tail_break_n=int(wcfg.tail_break_n))
    tpath = str(a.templates)
    is_npz = tpath.endswith(".npz")
    column = (a.column if is_npz or a.column != "auto" else
              mf_pulse.resolve_column(tpath, "auto", prefer="25320"))
    bank, info = mf_pulse.build_bank(names, pwo4, tpath, wave_cfg, column=column)
    kind = info.get("kind", "?")
    slug = a.slug or re.sub(r"^templates_", "", Path(tpath).stem)
    print(f"[setup] synth {a.synth_dir}: {len(cells)} cells, dt {dts}, ratio "
          f"{ratios}; run {run_num}; templates {tpath} (kind {kind}, sha1 "
          f"{info['sha1'][:10]}, {info['n_module']} modules from file, "
          f"{len(info['fallback'])} fallback)", flush=True)
    man_tk = man.get("template_kind", man.get("meta_template_kind"))
    man_sha = man.get("template_sha1", man.get("meta_template_sha1"))
    if man_sha and man_sha != info["sha1"]:
        print(f"[WARN] manifest template sha1 {str(man_sha)[:10]} != --templates "
              f"{info['sha1'][:10]}: the synthetic pulses were injected with a "
              f"different template than the mf bank", flush=True)

    rmods = [("headline", float(a.r_model))]
    if a.r_model_2 and a.r_model_2 > 0:
        rmods.append(("second", float(a.r_model_2)))
    ped_n = a.ped_nsamples if a.ped_nsamples > 0 else int(wcfg.ped_nsamples)
    recs = {tag: mf_pulse.MFReconstructor(
        bank, nsigma=a.nsigma, min_amp=a.min_amp, kmax=a.kmax, r_model=rm,
        ped_nsamples=ped_n) for tag, rm in rmods}

    ref_cell = min(cells, key=lambda c: (abs(c[1] - 0.1) + 1e-3 * abs(c[0] - 128),
                                         -c[0]))
    # ---- reference cell: onset zero-point between the bank and the truth ----
    dr = np.load(ref_cell[2], allow_pickle=False)
    midx_r = np.array([ix[n] for n in dr["base_channel_names"]], dtype=np.int64)
    okey = ("truth_onset_sampled_ns"
            if kind == "sampled" and "truth_onset_sampled_ns" in dr.files
            else "truth_onset_t0_ns")
    res_r = recs["headline"].run(dr["waveforms"], midx_r)
    conv_raw = onset_residuals(res_r, midx_r, bank, dr, okey)
    med = [conv_raw[t]["onset"].get("median", 0.0) for t in ("injected", "base")]
    if a.onset_shift == "auto":
        shift = (float(med[0]), float(med[1]))
    elif a.onset_shift == "none":
        shift = (0.0, 0.0)
    else:
        shift = tuple(float(x) for x in a.onset_shift.split(","))
        if len(shift) != 2:
            sys.exit("[ERROR] --onset-shift: none, auto or 'A,B'")
    conv = dict(raw=conv_raw, shift_ns=dict(injected=shift[0], base=shift[1]),
                applied=onset_residuals(res_r, midx_r, bank, dr, okey, shift=shift),
                cell=dict(dt_ns=ref_cell[0], ratio=ref_cell[1]))
    print(f"[convention] truth onset key {okey}; reference cell dt {ref_cell[0]:g} ratio {ref_cell[1]:g}: "
          f"median p_t0 - truth onset injected {med[0]:+.2f} ns, base "
          f"{med[1]:+.2f} ns; shift applied {shift[0]:+.2f}, {shift[1]:+.2f}",
          flush=True)
    if a.onset_shift == "none" and max(abs(med[0]), abs(med[1])) > 0.5 * a.tolerance_ns:
        print(f"[WARN] onset zero-point differs from the truth by more than half "
              f"the tolerance; rerun with --onset-shift auto if the bank and "
              f"the injected/base shapes use different onset definitions",
              flush=True)
    results = {tag: {} for tag, _ in rmods}
    results_w = {tag: {} for tag, _ in rmods}
    n2wa = {}
    fb_frac = {}
    for dt, ratio, pth in cells:
        d = np.load(pth, allow_pickle=False)
        wf = d["waveforms"]
        n_ev = wf.shape[0]
        miss = [n for n in set(d["base_channel_names"].tolist()) if n not in ix]
        if miss:
            sys.exit(f"[ERROR] {pth.name}: channels not in the geometry: {miss}")
        midx = np.array([ix[n] for n in d["base_channel_names"]], dtype=np.int64)
        t_inj = d[okey][:, 0].astype(np.float64) + shift[0]
        t_base = d[okey][:, 1].astype(np.float64) + shift[1]
        n_wa = int(sum(len(p.wave_ana.analyze(wf[i])[2]) >= 2
                       for i in range(n_ev)))
        n2wa[(dt, ratio)] = n_wa
        if "inject_template_fallback" in d.files:
            fb_frac[(dt, ratio)] = float(np.mean(d["inject_template_fallback"]))
        for tag, rm in rmods:
            res = recs[tag].run(wf, midx)
            clean = t_base >= 0.5 * ped_n * mf_pulse.DT
            results[tag][(dt, ratio)] = cell_metrics(
                res.p_ch, res.p_t0, n_ev, t_inj, t_base, a.tolerance_ns, clean, n_wa)
            if wrule is not None:
                kp = window_keep(res, bank, midx, n_ev, wrule)
                cw = cell_metrics(res.p_ch[kp], res.p_t0[kp], n_ev, t_inj, t_base,
                                  a.tolerance_ns, clean, n_wa)
                cw["n_pulses_dropped"] = int((~kp).sum())
                cw["n_events_with_drop"] = int(np.unique(res.p_ch[~kp]).size)
                results_w[tag][(dt, ratio)] = cw
        print(f"[cell] dt {dt:5.0f} ratio {ratio:4.2f} n {n_ev}: " + "  ".join(
            f"{tag}: eff {results[tag][(dt, ratio)]['eff_match']:.2f} fake "
            f"{results[tag][(dt, ratio)]['fake_rate']:.2f} n2mf "
            f"{results[tag][(dt, ratio)]['n2_mf']:.2f}" for tag, _ in rmods)
            + f"  n2wa {n_wa / n_ev:.2f}"
            + ("  | window rule: " + "  ".join(
                f"{tag}: eff {results_w[tag][(dt, ratio)]['eff_match']:.2f} fake "
                f"{results_w[tag][(dt, ratio)]['fake_rate']:.2f}" for tag, _ in rmods)
               if wrule is not None else ""), flush=True)

    # ----------------------------------------------------------- G2-2 ----
    g22 = None
    if a.table4_compare:
        head = results["headline"]
        rows, n_cmp, n_in2, n_inci = [], 0, 0, 0
        for dt, vals in TABLE4.items():
            for r, tv in zip(TABLE4_RATIOS, vals):
                key = (float(dt), float(r))
                if key not in head:
                    continue
                c = head[key]
                obs = 100.0 * c["n2_wa"]
                lo, hi = _wilson(c["n_wa2"], c["n_events"])
                ok2 = abs(obs - tv) <= 2.0 + 1e-9
                okci = lo * 100 - 1e-9 <= tv <= hi * 100 + 1e-9
                n_cmp += 1
                n_in2 += ok2
                n_inci += okci
                rows.append(dict(dt_ns=dt, ratio=r, table4=tv, n2_wa=obs,
                                 within_2pct=bool(ok2), within_ci95=bool(okci),
                                 ci95=[100 * lo, 100 * hi]))
        viol = []
        for tag, _ in rmods:
            for key, c in results[tag].items():
                if c["n2_mf"] > c["eff_match"] + c["fake_rate"] + 1e-12:
                    viol.append(dict(r_model=dict(rmods)[tag], dt_ns=key[0],
                                     ratio=key[1], n2_mf=c["n2_mf"],
                                     eff_match=c["eff_match"],
                                     fake_rate=c["fake_rate"]))
        g22 = dict(n_cells=n_cmp, n_within_2pct=int(n_in2),
                   frac_within_2pct=(n_in2 / n_cmp if n_cmp else None),
                   n_within_ci95=int(n_inci),
                   frac_within_ci95=(n_inci / n_cmp if n_cmp else None),
                   cells=rows, inequality_violations=viol,
                   inequality_note=(
                       "holds by construction while 2*tol < MIN_SEP 12 ns; "
                       "a bookkeeping check only"))

    # ------------------------------------------------- window-rule check --
    wchk = None
    if wrule is not None:
        wchk = {}
        for tag, _ in rmods:
            rows = []
            for key in sorted(results[tag]):
                dt, ratio = key
                base_c, win_c = results[tag][key], results_w[tag][key]
                rows.append(dict(
                    dt_ns=dt, ratio=ratio, eff_match=base_c["eff_match"],
                    eff_match_window=win_c["eff_match"],
                    d_eff_points=100.0 * (win_c["eff_match"] - base_c["eff_match"]),
                    n_match=base_c["n_match"], n_match_window=win_c["n_match"],
                    d_fake=win_c["fake_rate"] - base_c["fake_rate"],
                    n_pulses_dropped=win_c["n_pulses_dropped"]))
            sel = [r for r in rows if r["ratio"] >= 0.2 - 1e-9 and r["dt_ns"] >= 24 - 1e-9]
            tot_b = sum(r["n_match"] for r in sel)
            tot_w = sum(r["n_match_window"] for r in sel)
            wchk[tag] = dict(
                n_cells=len(sel),
                n_events=sum(results[tag][(r["dt_ns"], r["ratio"])]["n_events"] for r in sel),
                n_match=tot_b, n_match_window=tot_w,
                d_eff_points_pooled=(100.0 * (tot_w - tot_b)
                                     / max(sum(results[tag][(r["dt_ns"], r["ratio"])]["n_events"]
                                               for r in sel), 1)),
                n_cells_changed=sum(1 for r in sel if r["n_match"] != r["n_match_window"]),
                worst_cell_d_eff_points=(min(r["d_eff_points"] for r in sel) if sel else None),
                changed_cells=[r for r in sel if r["n_match"] != r["n_match_window"]],
                all_cells=rows)
            print(f"[window rule {wrule[0]:g}:{wrule[1]:g}] r_model {dict(rmods)[tag]:g}: "
                  f"cells ratio >= 0.2 & dt >= 24: {wchk[tag]['n_cells']} cells, "
                  f"eff_match change pooled {wchk[tag]['d_eff_points_pooled']:+.3f} points, "
                  f"{wchk[tag]['n_cells_changed']} cells changed, worst "
                  f"{wchk[tag]['worst_cell_d_eff_points']:+.2f} points", flush=True)

    # ---------------------------------------------------------- outputs --
    meta = dict(
        synth_dir=str(a.synth_dir), run_num=run_num, n_cells=len(cells),
        templates=tpath, template_sha1=info["sha1"], template_kind=kind,
        template_column=column, n_module_from_file=info["n_module"],
        n_fallback_modules=len(info["fallback"]),
        tolerance_ns=a.tolerance_ns, nsigma=a.nsigma, min_amp=a.min_amp,
        kmax=a.kmax, r_models=dict(rmods), ped_nsamples=ped_n,
        daq_config=daq or "(database default)",
        wave_analyzer=dict(smooth_order=wave_cfg["smooth_order"],
                           int_tail_ratio=wave_cfg["int_tail_ratio"],
                           tail_break_n=wave_cfg["tail_break_n"]),
        synth_manifest={k: v for k, v in man.items()
                        if k not in ("files", "cli")},
        synth_cli=man.get("cli"),
        synth_inject_fallback_frac_mean=(float(np.mean(list(fb_frac.values())))
                                         if fb_frac else None),
        onset_convention=(f"mf p_t0 = template onset (ns); truth = {okey} "
                          "(col 0 injected, col 1 base); matching on onsets, "
                          "+-tolerance"),
        truth_onset_key=okey,
        onset_shift_arg=a.onset_shift, convention_check=conv,
        window_rule=(dict(dt_max_ns=wrule[0], ratio=wrule[1]) if wrule else None),
        wall_s=None)
    meta["wall_s"] = time.monotonic() - t_start
    js = dict(settings=meta, dt_ns=dts, ratios=ratios, g2_2=g22, results={})
    for tag, rm in rmods:
        js["results"][tag] = dict(r_model=rm, cells=[
            dict(dt_ns=k[0], ratio=k[1], **v)
            for k, v in sorted(results[tag].items())])
    if wrule is not None:
        js["results_window"] = {tag: dict(r_model=rm, cells=[
            dict(dt_ns=k[0], ratio=k[1], **v)
            for k, v in sorted(results_w[tag].items())]) for tag, rm in rmods}
        js["window_rule_check"] = wchk
    (out / f"eff_{slug}.json").write_text(json.dumps(js, indent=1))

    L = [f"# Matched-filter pile-up efficiency: {slug}", "",
         f"- synthetic cells: `{a.synth_dir}` (run {run_num}, {len(cells)} cells, "
         f"{results['headline'][next(iter(results['headline']))]['n_events']} "
         f"events per cell)",
         f"- mf templates: `{tpath}` (kind {kind}, sha1 {info['sha1'][:12]}, "
         f"{info['n_module']} modules from file, {len(info['fallback'])} fallback)",
         f"- match tolerance +-{a.tolerance_ns:g} ns on the onset; nsigma "
         f"{a.nsigma:g}, min_amp {a.min_amp:g}, kmax {a.kmax}; MIN_SEP "
         f"{mf_pulse.MIN_SEP_NS:g} ns",
         "- eff_match, n2_mf, n2_wa in % of events; fake_rate in unmatched "
         "mf pulses per event; rows dt (ns, generator dt: injected peak - base CFD time), "
         "columns amplitude ratio injected / base",
         "- n2_wa is template- and r_model-independent (Table 4's column)",
         f"- *_clean: events with base onset >= {0.5 * ped_n * mf_pulse.DT:g} ns "
         f"(half the {ped_n}-sample mf pedestal window); n_clean per cell "
         f"{min(c['n_clean'] for c in results['headline'].values())} to "
         f"{max(c['n_clean'] for c in results['headline'].values())}; the four "
         "main columns use all events", ""]
    if wrule is not None:
        L += [f"- window rule {wrule[0]:g}:{wrule[1]:g} (post hoc, accepted pulses "
              "that are not the event's largest: 0 <= t - t_dom < DT_MAX ns and "
              "E/E_dom < RATIO dropped); the first table block of each r_model is "
              "without the rule, the second with it", ""]
    L += ["## Onset convention check", "",
          f"cell dt {conv['cell']['dt_ns']:g} ns, ratio {conv['cell']['ratio']:g}, "
          "headline r_model; residual = mf p_t0 - truth onset (ns), nearest "
          "pulse within 10 ns.  Shift applied to the truth onsets: injected "
          f"{conv['shift_ns']['injected']:+.2f}, base {conv['shift_ns']['base']:+.2f} ns "
          f"(--onset-shift {a.onset_shift})", "",
          "| truth | shift | n | median | sigma_MAD | rms | median if peak convention (p_t0 + t_pk - truth peak) |",
          "|---|---|---:|---:|---:|---:|---:|"]
    for which in ("raw", "applied"):
        for tag in ("injected", "base"):
            o, pk = conv[which][tag]["onset"], conv[which][tag]["peak_convention"]
            if o.get("n"):
                L.append(f"| {tag} | {which} | {o['n']} | {o['median']:+.2f} | "
                         f"{o['sigma_mad']:.2f} | {o['rms']:.2f} | "
                         f"{pk['median']:+.2f} |")
    L.append("")
    for tag, rm in rmods:
        L += [f"## r_model = {rm:g} ({tag})", ""]
        for key in METRICS:
            g = grid(results[tag], dts, ratios, key)
            scale = 1.0 if key.startswith("fake_rate") else 100.0
            fmt = "{:.2f}" if key.startswith("fake_rate") else "{:.0f}"
            L += [f"### {TITLES[key]}", "",
                  "| dt (ns) | " + " | ".join(f"{r:g}" for r in ratios) + " |",
                  "|---:|" + "---:|" * len(ratios)]
            for i, dt in enumerate(dts):
                L.append(f"| {dt:g} | " + " | ".join(
                    fmt.format(scale * g[i, j]) for j in range(len(ratios))) + " |")
            L.append("")
        if wrule is not None:
            L += [f"## r_model = {rm:g} ({tag}), with window rule "
                  f"(0 <= dt < {wrule[0]:g} ns and E/E_dom < {wrule[1]:g} dropped)", ""]
            for key in METRICS:
                if key == "n2_wa":
                    continue
                g = grid(results_w[tag], dts, ratios, key)
                scale = 1.0 if key.startswith("fake_rate") else 100.0
                fmt = "{:.2f}" if key.startswith("fake_rate") else "{:.0f}"
                L += [f"### {TITLES[key]} with window rule", "",
                      "| dt (ns) | " + " | ".join(f"{r:g}" for r in ratios) + " |",
                      "|---:|" + "---:|" * len(ratios)]
                for i, dt in enumerate(dts):
                    L.append(f"| {dt:g} | " + " | ".join(
                        fmt.format(scale * g[i, j]) for j in range(len(ratios))) + " |")
                L.append("")
            g0 = grid(results[tag], dts, ratios, "eff_match")
            g1 = grid(results_w[tag], dts, ratios, "eff_match")
            L += ["### eff_match change with the window rule (points)", "",
                  "| dt (ns) | " + " | ".join(f"{r:g}" for r in ratios) + " |",
                  "|---:|" + "---:|" * len(ratios)]
            for i, dt in enumerate(dts):
                L.append(f"| {dt:g} | " + " | ".join(
                    f"{100 * (g1[i, j] - g0[i, j]):+.1f}" for j in range(len(ratios))) + " |")
            w = wchk[tag]
            L += ["", f"cells with ratio >= 0.2 and dt >= 24 ns: {w['n_cells']} cells, "
                      f"{w['n_events']} events, matched {w['n_match']} -> "
                      f"{w['n_match_window']} (pooled {w['d_eff_points_pooled']:+.3f} "
                      f"points; {w['n_cells_changed']} cells changed; worst cell "
                      f"{w['worst_cell_d_eff_points']:+.2f} points)", ""]
    if g22:
        L += ["## G2-2: n2_wa vs Table 4 (fit_pulse_template.md section 6)", "",
              f"{g22['n_within_2pct']}/{g22['n_cells']} cells within 2 % "
              f"({100 * g22['frac_within_2pct']:.0f} %); "
              f"{g22['n_within_ci95']}/{g22['n_cells']} with the Table 4 value "
              f"inside the Wilson 95 % interval of the observed fraction "
              f"({100 * g22['frac_within_ci95']:.0f} %). Table 4 prints "
              "ratios 0.1, 0.3, 0.5, 1, 2 and dt 8, 16, 20, 24, 32, 48, 64, 96, "
              "128 only.", "",
              "| dt (ns) | ratio | Table 4 % | n2_wa % | within 2 % | in 95 % CI |",
              "|---:|---:|---:|---:|:-:|:-:|"]
        for r in g22["cells"]:
            L.append(f"| {r['dt_ns']:g} | {r['ratio']:g} | {r['table4']:g} | "
                     f"{r['n2_wa']:.0f} | {'y' if r['within_2pct'] else 'N'} | "
                     f"{'y' if r['within_ci95'] else 'N'} |")
        L += ["", f"n2_mf <= eff_match + fake_rate violations: "
                  f"{len(g22['inequality_violations'])}", "",
              g22["inequality_note"], ""]
    (out / f"eff_{slug}.md").write_text("\n".join(L) + "\n")

    # heat maps
    nr = len(rmods)
    fig, axs = plt.subplots(nr, len(METRICS), figsize=(4.0 * len(METRICS), 4.2 * nr), squeeze=False,
                            constrained_layout=True)
    for ri, (tag, rm) in enumerate(rmods):
        for ci, key in enumerate(METRICS):
            ax = axs[ri, ci]
            g = grid(results[tag], dts, ratios, key)
            scale = 1.0 if key.startswith("fake_rate") else 100.0
            vmax = (max(0.5, np.nanmax(g)) if key.startswith("fake_rate")
                    else 1.0)
            im = ax.imshow(scale * g, origin="lower", aspect="auto",
                           cmap="viridis", vmin=0, vmax=scale * vmax)
            ax.set_xticks(range(len(ratios)), [f"{r:g}" for r in ratios])
            ax.set_yticks(range(len(dts)), [f"{x:g}" for x in dts])
            ax.set_xlabel("amplitude ratio")
            if ci == 0:
                ax.set_ylabel(f"r_model {rm:g}\ndt (ns)")
            ax.set_title(TITLES[key].split(" (")[0], fontsize=10)
            fmt = "{:.2f}" if key.startswith("fake_rate") else "{:.0f}"
            for i in range(len(dts)):
                for j in range(len(ratios)):
                    v = scale * g[i, j]
                    if np.isnan(v):
                        continue
                    ax.text(j, i, fmt.format(v), ha="center", va="center",
                            fontsize=7,
                            color="white" if v < 0.6 * scale * vmax else "black")
            fig.colorbar(im, ax=ax, fraction=0.046, pad=0.02)
    fig.suptitle(f"{slug}: +-{a.tolerance_ns:g} ns onset match, run {run_num}, "
                 f"{kind} templates", fontsize=11)
    fig.savefig(out / f"eff_{slug}.png", dpi=130)
    plt.close(fig)

    if g22:
        print("\n[G2-2] n2_wa vs Table 4 (% events with >= 2 findPeaks peaks)")
        print(f"  within 2 %: {g22['n_within_2pct']}/{g22['n_cells']} "
              f"({100 * g22['frac_within_2pct']:.0f} %)  (* = outside)")
        print(f"  Table 4 inside Wilson 95 % CI of ours: {g22['n_within_ci95']}/"
              f"{g22['n_cells']} ({100 * g22['frac_within_ci95']:.0f} %)")
        print(f"  n2_mf <= eff_match + fake_rate violations: "
              f"{len(g22['inequality_violations'])}")
        for v in g22["inequality_violations"]:
            print(f"    r_model {v['r_model']:g} dt {v['dt_ns']:g} ratio "
                  f"{v['ratio']:g}: n2_mf {v['n2_mf']:.3f} > eff "
                  f"{v['eff_match']:.3f} + fake {v['fake_rate']:.3f}")
    print(f"[done] {out}/eff_{slug}.{{json,md,png}}  wall "
          f"{time.monotonic() - t_start:.1f} s")


if __name__ == "__main__":
    main()
