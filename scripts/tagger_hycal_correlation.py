#!/usr/bin/env python3
"""
tagger_hycal_correlation.py
===========================

Two-step tagger → HyCal coincidence analysis.

Step 1  For each pair (T10R, Eₓ) with x ∈ {49, 50, 51, 52, 53}, fill the
        event-wise TDC-difference histogram ``ΔT = tdc(T10R) − tdc(Eₓ)``
        and fit a Gaussian around the coincidence peak.

Step 2  For each pair, apply a ±Nσ cut around the peak (default 3σ) and
        fill the peak-height / peak-integral histograms of HyCal module
        **W1156** (ROC 0x8C slot 7 ch 3) for the events that survive the
        cut.

Input is a raw evio file decoded in-process through ``prad2py.dec`` — no
replay ROOT file needed.  Output is a single ROOT file containing all
ΔT histograms with fit functions attached plus the selected-event
W1156 histograms.

Usage
-----
    python scripts/tagger_hycal_correlation.py \
        /data/stage6/prad_023671/prad_023671.evio.00000 \
        -o tagger_w1156_corr.root \
        -n 500000

The script depends on the ``prad2py`` pybind11 module (build with
``-DBUILD_PYTHON=ON``) and on PyROOT.  Simple pedestal-and-max peak
finding is done inline — good enough for monitoring-style plots;
use the full WaveAnalyzer / Replay machinery if you need calibrated
quantities.

Hard-coded channel layout (update if the DAQ map changes):

    TAGGER crate (ROC 0x008E), slot 18 in the V1190:
      T10R  = ch  0       (reference)
      E49   = ch 11
      E50   = ch 12
      E51   = ch 13
      E52   = ch 14
      E53   = ch 15

    HyCal module W1156: ROC 0x008C, slot 7, channel 3
"""
from __future__ import annotations

import argparse
import os
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Make the freshly-built prad2py module importable without the user
# having to set PYTHONPATH.  Mirrors the logic in scripts/tdc_viewer.py.
_HERE = Path(__file__).resolve().parent
for _cand in (
    _HERE.parent / "build" / "python",
    _HERE.parent / "build-release" / "python",
):
    if _cand.is_dir() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))

try:
    import prad2py                        # noqa: F401
    from prad2py import dec
except Exception as exc:                  # noqa: BLE001
    sys.exit(f"cannot import prad2py: {exc}\n"
             "build with -DBUILD_PYTHON=ON and retry.")

try:
    import ROOT
except ImportError:
    sys.exit("this script needs PyROOT (the 'ROOT' python module).")


# ---------------------------------------------------------------------------
# Channel layout — adjust if the DAQ map changes
# ---------------------------------------------------------------------------

TAGGER_SLOT   = 18
T10R_CH       = 0
E_CHANNELS: List[Tuple[str, int]] = [
    ("E49", 11), ("E50", 12), ("E51", 13), ("E52", 14), ("E53", 15),
]

W1156_ROC     = 0x8C
W1156_SLOT    = 7
W1156_CHANNEL = 3

# FADC peak finder — pedestal window and integration half-width (in samples).
PED_WINDOW    = 10
INT_HALFWIDTH = 8

# Histogram axis defaults.
W1156_H_BINS  = (200, 0.0, 4000.0)      # peak height in ADC counts
W1156_I_BINS  = (200, 0.0, 40000.0)     # peak integral in ADC·sample


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def first_tdc(tdc_evt, slot: int, channel: int) -> Optional[int]:
    """Earliest (smallest-TDC) hit on (slot, channel) for this event, or None.

    V1190 can fire multiple times per channel per event; we keep the first
    one — the usual convention for coincidence timing.
    """
    best: Optional[int] = None
    for i in range(tdc_evt.n_hits):
        h = tdc_evt.hit(i)
        if h.slot == slot and h.channel == channel:
            if best is None or h.value < best:
                best = int(h.value)
    return best


def hycal_peak(samples: np.ndarray) -> Optional[Tuple[float, float, int]]:
    """Return (height, integral, peak-bin) from a raw FADC250 sample array.

    Pedestal = mean of the first PED_WINDOW samples.  Height = maximum of
    (samples − ped).  Integral = sum of (samples − ped) in a window of
    ±INT_HALFWIDTH around the maximum.  Returns None on empty input.
    """
    if samples.size < PED_WINDOW + 1:
        return None
    ped = float(np.mean(samples[:PED_WINDOW]))
    sub = samples.astype(np.float64) - ped
    tmax = int(np.argmax(sub))
    height = float(sub[tmax])
    lo = max(0, tmax - INT_HALFWIDTH)
    hi = min(samples.size, tmax + INT_HALFWIDTH + 1)
    integral = float(np.sum(sub[lo:hi]))
    return height, integral, tmax


def w1156_peak(event) -> Optional[Tuple[float, float]]:
    """Extract (height, integral) for W1156 from a prad2py EventData."""
    roc = event.find_roc(W1156_ROC)
    if roc is None:
        return None
    slot = roc.slot(W1156_SLOT)
    if not slot.present:
        return None
    ch = slot.channel(W1156_CHANNEL)
    if ch.nsamples <= 0:
        return None
    pk = hycal_peak(ch.samples)
    if pk is None:
        return None
    return pk[0], pk[1]


# ---------------------------------------------------------------------------
# First-pass event collection
# ---------------------------------------------------------------------------

def collect(ch: "dec.EvChannel", max_events: int) -> Dict[str, np.ndarray]:
    """Loop the evio file once and return parallel numpy arrays, one entry
    per accepted event.  An event is accepted when it has *at least* one
    T10R TDC hit and a W1156 sample present — channels Eₓ may be missing
    (encoded as -1) and the per-pair fits will mask them out later.
    """
    t10r_list: List[int] = []
    e_lists: Dict[str, List[int]] = {name: [] for name, _ in E_CHANNELS}
    evnum_list: List[int] = []
    w_height: List[float] = []
    w_integral: List[float] = []

    n = 0
    while ch.read() == dec.Status.success:
        if not ch.scan() or ch.get_event_type() != dec.EventType.Physics:
            continue
        for i in range(ch.get_n_events()):
            e = ch.decode_event(i, with_tdc=True)
            if not e["ok"]:
                continue

            t0 = first_tdc(e["tdc"], TAGGER_SLOT, T10R_CH)
            if t0 is None:
                continue

            w = w1156_peak(e["event"])
            if w is None:
                continue

            # At least one Eₓ must fire to be worth keeping.
            e_vals = {name: first_tdc(e["tdc"], TAGGER_SLOT, cch)
                      for name, cch in E_CHANNELS}
            if all(v is None for v in e_vals.values()):
                continue

            t10r_list.append(t0)
            for name in e_lists:
                v = e_vals[name]
                e_lists[name].append(v if v is not None else -1)
            evnum_list.append(int(e["event"].info.event_number))
            w_height.append(w[0])
            w_integral.append(w[1])

            n += 1
            if n % 100_000 == 0:
                print(f"  pass 1: {n:>10,} events collected", flush=True)
            if max_events and n >= max_events:
                return {
                    "t10r":     np.asarray(t10r_list,   dtype=np.int64),
                    "evnum":    np.asarray(evnum_list,  dtype=np.int64),
                    "heights":  np.asarray(w_height,    dtype=np.float32),
                    "integrals":np.asarray(w_integral,  dtype=np.float32),
                    **{name: np.asarray(v, dtype=np.int64) for name, v in e_lists.items()},
                }
    return {
        "t10r":     np.asarray(t10r_list,   dtype=np.int64),
        "evnum":    np.asarray(evnum_list,  dtype=np.int64),
        "heights":  np.asarray(w_height,    dtype=np.float32),
        "integrals":np.asarray(w_integral,  dtype=np.float32),
        **{name: np.asarray(v, dtype=np.int64) for name, v in e_lists.items()},
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main(argv=None) -> int:
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("evio", help="Input evio file")
    ap.add_argument("-o", "--out", default="tagger_w1156_corr.root",
                    help="Output ROOT file (default: tagger_w1156_corr.root)")
    ap.add_argument("-n", "--max-events", type=int, default=0,
                    help="Stop after N physics events (0 = all)")
    ap.add_argument("--daq-config", default="",
                    help="Override daq_config.json path (empty = default)")
    ap.add_argument("--dt-bins",  type=int, default=400)
    ap.add_argument("--dt-range", type=float, default=200.0,
                    help="Half-range for ΔT histograms in TDC-LSB units")
    ap.add_argument("--nsigma",   type=float, default=3.0,
                    help="Timing cut half-width, in σ of the Gaussian peak")
    args = ap.parse_args(argv)

    # --- open evio -------------------------------------------------------
    cfg = dec.load_daq_config(args.daq_config)
    ch  = dec.EvChannel()
    ch.set_config(cfg)
    st = ch.open(args.evio)
    if st != dec.Status.success:
        sys.exit(f"cannot open {args.evio}: {st}")

    print(f"reading {args.evio}")
    data = collect(ch, args.max_events)
    ch.close()

    n = data["t10r"].size
    print(f"collected {n:,} events "
          f"(T10R + W1156 present, ≥1 Eₓ)")
    if n == 0:
        sys.exit("no events survived the initial filter — nothing to plot")

    # --- build histograms, fit, cut, fill W1156 --------------------------
    out = ROOT.TFile(args.out, "RECREATE")
    ROOT.gStyle.SetOptFit(1)
    ROOT.gStyle.SetOptStat(1110)

    results: Dict[str, Dict[str, float]] = {}
    canvas = ROOT.TCanvas("summary", "tagger–W1156 correlations", 1400, 900)
    canvas.Divide(len(E_CHANNELS), 3)

    for idx, (ename, _ech) in enumerate(E_CHANNELS, start=1):
        te  = data[ename]
        mask = te >= 0
        dt   = (data["t10r"][mask] - te[mask]).astype(np.float64)

        # --- ΔT histogram + Gaussian fit around the dominant peak -------
        hdt = ROOT.TH1D(
            f"dt_T10R_{ename}",
            f"#DeltaT = T10R - {ename};tdc(T10R) - tdc({ename}) [LSB];events",
            args.dt_bins, -args.dt_range, args.dt_range,
        )
        for v in dt:
            hdt.Fill(float(v))

        peak_bin = hdt.GetMaximumBin()
        peak_x   = hdt.GetXaxis().GetBinCenter(peak_bin)
        bw       = hdt.GetXaxis().GetBinWidth(1)

        fit = ROOT.TF1(f"gfit_{ename}", "gaus",
                       peak_x - 20 * bw, peak_x + 20 * bw)
        hdt.Fit(fit, "RQ", "", peak_x - 20 * bw, peak_x + 20 * bw)
        mu    = float(fit.GetParameter(1))
        sigma = max(abs(float(fit.GetParameter(2))), bw)  # guard against zero σ

        hdt.Write()
        fit.Write()

        # --- apply timing cut ------------------------------------------
        passes_dt = np.abs(dt - mu) < args.nsigma * sigma
        sel = np.where(mask)[0][passes_dt]
        n_sel = int(sel.size)

        h_height = ROOT.TH1F(
            f"W1156_height_{ename}",
            f"W1156 peak height, |#DeltaT - {mu:.1f}| < {args.nsigma:.1f}#sigma "
            f"(T10R-{ename});height [ADC];events",
            *W1156_H_BINS,
        )
        h_integ  = ROOT.TH1F(
            f"W1156_integral_{ename}",
            f"W1156 peak integral, |#DeltaT - {mu:.1f}| < {args.nsigma:.1f}#sigma "
            f"(T10R-{ename});integral [ADC#upoint sample];events",
            *W1156_I_BINS,
        )
        if n_sel > 0:
            for k in sel:
                h_height.Fill(float(data["heights"][k]))
                h_integ.Fill(float(data["integrals"][k]))
        h_height.Write()
        h_integ.Write()

        # --- summary canvas --------------------------------------------
        canvas.cd(idx)
        hdt.Draw()
        canvas.cd(idx + len(E_CHANNELS))
        h_height.Draw()
        canvas.cd(idx + 2 * len(E_CHANNELS))
        h_integ.Draw()

        results[ename] = {"mu": mu, "sigma": sigma, "n_sel": n_sel,
                          "n_total": int(mask.sum())}

    canvas.Write()
    out.Close()

    # --- terminal summary -----------------------------------------------
    print()
    print("=== Summary ===")
    print(f"{'pair':>8}  {'μ[LSB]':>10}  {'σ[LSB]':>9}  "
          f"{'n_total':>10}  {'n_selected':>12}  {'keep':>6}")
    for ename, r in results.items():
        frac = (100.0 * r["n_sel"] / r["n_total"]) if r["n_total"] else 0.0
        print(f"T10R-{ename:<3s}  {r['mu']:10.2f}  {r['sigma']:9.2f}  "
              f"{r['n_total']:10,}  {r['n_sel']:12,}  {frac:5.1f}%")
    print()
    print(f"histograms written to {args.out}")
    print("hints: open in ROOT with")
    print(f"  root -l {args.out}")
    print("  root [0] summary->Draw()")
    print("  root [1] dt_T10R_E51->Draw()")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
