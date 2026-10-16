#!/usr/bin/env python3
"""
npeaks_distribution.py — peaks-per-channel distribution across EVIO runs.

For each run, scan physics events, run WaveAnalyzer.analyze_result() on
every channel with waveform data, and count how many channels have exactly
k peaks (k = 0, 1, 2, ..., max_npeaks).  Characterises the pile-up
environment and reports two complementary fractions:

  Denominator 1 — channel-weighted:
    fraction(k) = N(channels with k peaks) / N(all read-out channels)
    Answers: what fraction of read-out channels have k peaks?
    Sums to 100 % over all k including k = 0.

  Denominator 3 — pulse-weighted:
    fraction(k) = (k × N(channels with k peaks)) / (total pulses)
    Answers: what fraction of all reconstructed pulses live in a k-peak channel?
    k = 0 contributes zero.  Sums to 100 % over k ≥ 1.

Usage
-----
    python npeaks_distribution.py <evio_path> [<evio_path> ...] \\
        [--max-events N] [--max-npeaks N] [--daq-config PATH] \\
        [--out out.json] [--material TYPE] [--hc-map-file PATH]

Each positional argument is one run column.  It may be a single split, a
glob, or a directory; all matching splits are concatenated into one run.
The run label is sniffed from the filename (prad_NNNNNN.evio.*).
"""

from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time
from pathlib import Path
from typing import Dict, List, Optional, Tuple


# ---------------------------------------------------------------------------
# Run-label helper
# ---------------------------------------------------------------------------

_RUN_PAT = re.compile(r"prad[_-](\d+)\.evio", re.IGNORECASE)


def _sniff_run_label(path: str) -> str:
    """Return zero-padded run number from 'prad_NNNNNN.evio.*', else basename."""
    m = _RUN_PAT.search(Path(path).name)
    if m:
        return m.group(1).zfill(6)
    return Path(path).name


def _resolve_evio_paths(raw: str) -> List[str]:
    """Expand a single CLI argument to an ordered, deduped list of EVIO paths.

    Three modes:
      * glob  ('*' in the string) — Python glob expansion, sorted.
      * directory                 — enumerate prad_NNNNNN.evio.<digits> entries.
      * anything else             — single-file pass-through.
    """
    if not raw:
        return []
    p = Path(raw)
    if "*" in raw:
        import glob as _glob
        matches = sorted(_glob.glob(raw))
        return matches if matches else [raw]
    if p.is_dir():
        import glob as _glob
        matches = sorted(_glob.glob(str(p / "prad_*.evio.*")))
        return matches if matches else [raw]
    return [raw]


# ---------------------------------------------------------------------------
# Per-run scanning
# ---------------------------------------------------------------------------

def scan_run(
    label: str,
    evio_files: List[str],
    dec,                         # prad2py.dec module
    wa,                          # dec.WaveAnalyzer instance
    cfg,                         # dec.DaqConfig instance
    max_events: int,
    max_npeaks: int,
    material_filter: Optional[str],
    hycal,                       # det.HyCalSystem or None
    crate_map: Dict[int, int],   # roc_tag -> crate
) -> Optional[Dict]:
    """Scan one run's EVIO files and return a result dict, or None on hard failure."""

    # counts[k] = number of (channel, event) pairs with exactly k peaks.
    # Index 0 .. max_npeaks; index max_npeaks means "≥ max_npeaks" (clamped).
    counts: List[int] = [0] * (max_npeaks + 1)
    n_events = 0

    # Event-level pile-up accumulators.
    n_events_with_pileup: int = 0
    # event_max_npeaks_hist[k] = events whose max channel npeaks == k (clamped).
    event_max_npeaks_hist: List[int] = [0] * (max_npeaks + 1)
    # event_n_piled_channels_hist[b] = events with b channels having ≥2 peaks.
    # Bins 0..4 are exact; bin 5 means "5+" (clamped).
    event_n_piled_channels_hist: List[int] = [0] * 6

    ch = dec.EvChannel()
    ch.set_config(cfg)

    # Module-type cache: (roc_tag, slot, chan) -> accepted bool
    _accept_cache: Dict[Tuple[int, int, int], bool] = {}

    def _accepted(roc_tag: int, crate: Optional[int], s: int, c: int) -> bool:
        """Return True if this channel passes the material filter (or no filter)."""
        if material_filter is None:
            return True
        key = (roc_tag, s, c)
        cached = _accept_cache.get(key)
        if cached is not None:
            return cached
        ok = False
        if crate is not None and hycal is not None:
            mod = hycal.module_by_daq(crate, s, c)
            if mod is not None:
                ok = (mod.type.name == material_filter)
        _accept_cache[key] = ok
        return ok

    n_files_ok = 0
    done = False
    t0_wall = time.monotonic()
    next_progress = 2000

    for fpath in evio_files:
        if ch.open_auto(fpath) != dec.Status.success:
            print(f"[ERROR] cannot open {fpath}", flush=True)
            continue
        n_files_ok += 1
        print(f"  [run {label}] opened {fpath}", flush=True)

        while ch.read() == dec.Status.success:
            if not ch.scan():
                continue
            if ch.get_event_type() != dec.EventType.Physics:
                continue

            for i in range(ch.get_n_events()):
                decoded = ch.decode_event(i, with_ssp=False)
                if not decoded["ok"]:
                    continue
                n_events += 1

                fadc_evt = decoded["event"]

                # Per-event pile-up tracking.
                event_max_np: int = 0
                event_n_piled: int = 0

                for ri in range(fadc_evt.nrocs):
                    roc = fadc_evt.roc(ri)
                    if not roc.present:
                        continue
                    crate = crate_map.get(roc.tag)
                    for s in roc.present_slots():
                        slot = roc.slot(s)
                        for c in slot.present_channels():
                            cd = slot.channel(c)
                            if cd.nsamples <= 0:
                                continue
                            if not _accepted(roc.tag, crate, s, c):
                                continue

                            import numpy as np
                            samples = np.asarray(cd.samples, dtype=np.uint16)
                            wres = wa.analyze_result(samples)
                            np_found = min(len(wres.peaks), max_npeaks)
                            counts[np_found] += 1
                            if np_found > event_max_np:
                                event_max_np = np_found
                            if np_found >= 2:
                                event_n_piled += 1

                # After scanning all channels in this event, update event-level hists.
                event_max_npeaks_hist[event_max_np] += 1
                piled_bin = min(event_n_piled, 5)
                event_n_piled_channels_hist[piled_bin] += 1
                if event_n_piled >= 1:
                    n_events_with_pileup += 1

                if n_events >= next_progress:
                    elapsed = time.monotonic() - t0_wall
                    rate = n_events / elapsed if elapsed > 0 else 0.0
                    print(
                        f"  [run {label}] events={n_events}  "
                        f"rate={rate:.0f} ev/s  elapsed={elapsed:.1f}s",
                        flush=True,
                    )
                    while next_progress <= n_events:
                        next_progress += 2000

                if max_events and n_events >= max_events:
                    done = True
                    break
            if done:
                break

        ch.close()
        if done:
            break

    if n_files_ok == 0:
        # Every file in this run failed to open.
        return None

    elapsed = time.monotonic() - t0_wall
    print(
        f"  [run {label}] done — {n_events} events  "
        f"{sum(counts)} channel-hits  {elapsed:.1f}s",
        flush=True,
    )
    return {
        "counts": counts,
        "n_events": n_events,
        "n_events_with_pileup": n_events_with_pileup,
        "event_max_npeaks_hist": event_max_npeaks_hist,
        "event_n_piled_channels_hist": event_n_piled_channels_hist,
    }


# ---------------------------------------------------------------------------
# Formatting helpers
# ---------------------------------------------------------------------------

def _fmt_int(n: int) -> str:
    return f"{n:,}"


def _fmt_pct(v: float) -> str:
    return f"{v:.1f}%"


def _col_widths(
    run_labels: List[str],
    col_count_w: int = 14,
    col_pct_w: int = 12,
) -> Tuple[int, int]:
    """Return (count_col_width, pct_col_width) wide enough for the header."""
    for lab in run_labels:
        col_count_w = max(col_count_w, len(f"{lab} (count)") + 2)
        col_pct_w   = max(col_pct_w,   len(f"{lab} (%)") + 2)
    return col_count_w, col_pct_w


# ---------------------------------------------------------------------------
# Table printing — Denominator 1 (channel-weighted)
# ---------------------------------------------------------------------------

def print_table_denom1(
    run_labels: List[str],
    results: Dict[str, Optional[Dict]],
    max_npeaks: int,
) -> None:
    print()
    print("=== Denominator 1: fraction of read-out CHANNELS with k peaks ===")

    cw, pw = _col_widths(run_labels)
    npeaks_w = 8

    # Header
    hdr = f"{'npeaks':>{npeaks_w}}"
    for lab in run_labels:
        hdr += f"  {(lab + ' (count)'):>{cw}}"
        hdr += f"  {(lab + ' (%)'):>{pw}}"
    print(hdr)
    print("-" * len(hdr))

    # Rows k = 0 .. max_npeaks
    for k in range(max_npeaks + 1):
        label_k = f"{k}+" if k == max_npeaks else str(k)
        row = f"{label_k:>{npeaks_w}}"
        for lab in run_labels:
            res = results.get(lab)
            if res is None:
                row += f"  {'—':>{cw}}  {'—':>{pw}}"
            else:
                n_total = sum(res["counts"])
                cnt = res["counts"][k]
                pct = 100.0 * cnt / n_total if n_total else 0.0
                row += f"  {_fmt_int(cnt):>{cw}}  {_fmt_pct(pct):>{pw}}"
        print(row)

    # Total row
    row = f"{'total':>{npeaks_w}}"
    for lab in run_labels:
        res = results.get(lab)
        if res is None:
            row += f"  {'—':>{cw}}  {'—':>{pw}}"
        else:
            n_total = sum(res["counts"])
            row += f"  {_fmt_int(n_total):>{cw}}  {'100.0%':>{pw}}"
    print(row)


# ---------------------------------------------------------------------------
# Table printing — Denominator 3 (pulse-weighted)
# ---------------------------------------------------------------------------

def print_table_denom3(
    run_labels: List[str],
    results: Dict[str, Optional[Dict]],
    max_npeaks: int,
) -> None:
    print()
    print("=== Denominator 3: fraction of PULSES from k-peak channels ===")

    # Column header labels differ from denom1
    col_count_w = max(14, max((len(f"{lab} (pulses)") + 2) for lab in run_labels))
    col_pct_w   = max(12, max((len(f"{lab} (%)") + 2) for lab in run_labels))
    npeaks_w = 8

    hdr = f"{'npeaks':>{npeaks_w}}"
    for lab in run_labels:
        hdr += f"  {(lab + ' (pulses)'):>{col_count_w}}"
        hdr += f"  {(lab + ' (%)'):>{col_pct_w}}"
    print(hdr)
    print("-" * len(hdr))

    # k = 1 .. max_npeaks only (k=0 contributes zero pulses)
    for k in range(1, max_npeaks + 1):
        label_k = f"{k}+" if k == max_npeaks else str(k)
        row = f"{label_k:>{npeaks_w}}"
        for lab in run_labels:
            res = results.get(lab)
            if res is None:
                row += f"  {'—':>{col_count_w}}  {'—':>{col_pct_w}}"
            else:
                counts = res["counts"]
                total_pulses = sum(kk * counts[kk] for kk in range(1, max_npeaks + 1))
                pulses_k = k * counts[k]
                pct = 100.0 * pulses_k / total_pulses if total_pulses else 0.0
                row += f"  {_fmt_int(pulses_k):>{col_count_w}}  {_fmt_pct(pct):>{col_pct_w}}"
        print(row)

    # Total row
    row = f"{'total':>{npeaks_w}}"
    for lab in run_labels:
        res = results.get(lab)
        if res is None:
            row += f"  {'—':>{col_count_w}}  {'—':>{col_pct_w}}"
        else:
            counts = res["counts"]
            total_pulses = sum(kk * counts[kk] for kk in range(1, max_npeaks + 1))
            row += f"  {_fmt_int(total_pulses):>{col_count_w}}  {'100.0%':>{col_pct_w}}"
    print(row)


# ---------------------------------------------------------------------------
# Pile-up summary
# ---------------------------------------------------------------------------

def print_summary(
    run_labels: List[str],
    results: Dict[str, Optional[Dict]],
    max_npeaks: int,
) -> None:
    print()
    print("=== Pile-up fraction summary ===")
    for lab in run_labels:
        res = results.get(lab)
        if res is None:
            print(f"  {lab}: (no data — run failed to open)")
            continue
        counts = res["counts"]
        n_channels = sum(counts)
        n_pileup_ch = sum(counts[k] for k in range(2, max_npeaks + 1))
        total_pulses = sum(k * counts[k] for k in range(1, max_npeaks + 1))
        pileup_pulses = sum(k * counts[k] for k in range(2, max_npeaks + 1))
        pw_ch  = 100.0 * n_pileup_ch   / n_channels     if n_channels     else 0.0
        pw_pls = 100.0 * pileup_pulses / total_pulses    if total_pulses   else 0.0
        print(
            f"  {lab}: channel-weighted = {pw_ch:.1f}%, "
            f"pulse-weighted = {pw_pls:.1f}%  "
            f"({_fmt_int(n_channels)} channels, {_fmt_int(total_pulses)} pulses)"
        )


# ---------------------------------------------------------------------------
# Table printing — event-level pile-up (three new tables)
# ---------------------------------------------------------------------------

def print_table_event_pileup_summary(
    run_labels: List[str],
    results: Dict[str, Optional[Dict]],
) -> None:
    """Table 3a — how many events (and what fraction) contain ≥1 piled channel."""
    print()
    print("=== Event-level pile-up ===")

    run_w = max(6, max(len(lab) for lab in run_labels)) + 2
    ev_w  = max(8,  len("events"))   + 2
    pu_w  = max(14, len("with pile-up")) + 2
    fr_w  = max(10, len("fraction"))  + 2

    hdr = (
        f"  {'Run':>{run_w}}"
        f"  {'events':>{ev_w}}"
        f"  {'with pile-up':>{pu_w}}"
        f"  {'fraction':>{fr_w}}"
    )
    print(hdr)
    print("  " + "-" * (len(hdr) - 2))

    for lab in run_labels:
        res = results.get(lab)
        if res is None:
            print(
                f"  {lab:>{run_w}}  {'—':>{ev_w}}  {'—':>{pu_w}}  {'—':>{fr_w}}"
            )
            continue
        n_ev = res["n_events"]
        n_pu = res["n_events_with_pileup"]
        frac = 100.0 * n_pu / n_ev if n_ev else 0.0
        print(
            f"  {lab:>{run_w}}"
            f"  {_fmt_int(n_ev):>{ev_w}}"
            f"  {_fmt_int(n_pu):>{pu_w}}"
            f"  {_fmt_pct(frac):>{fr_w}}"
        )


def print_table_event_max_npeaks(
    run_labels: List[str],
    results: Dict[str, Optional[Dict]],
    max_npeaks: int,
) -> None:
    """Table 3b — distribution of the max npeaks per event."""
    print()
    print("=== Max npeaks per event ===")

    col_w = max(12, max(len(f"{lab} (%)") + 2 for lab in run_labels))
    key_w = 6

    hdr = f"{'max_np':>{key_w}}"
    for lab in run_labels:
        hdr += f"  {(lab + ' (%)'):>{col_w}}"
    print(hdr)
    print("-" * len(hdr))

    for k in range(max_npeaks + 1):
        label_k = f"{k}+" if k == max_npeaks else str(k)
        row = f"{label_k:>{key_w}}"
        for lab in run_labels:
            res = results.get(lab)
            if res is None:
                row += f"  {'—':>{col_w}}"
            else:
                n_ev = res["n_events"]
                cnt  = res["event_max_npeaks_hist"][k]
                pct  = 100.0 * cnt / n_ev if n_ev else 0.0
                row += f"  {_fmt_pct(pct):>{col_w}}"
        print(row)


def print_table_event_n_piled_channels(
    run_labels: List[str],
    results: Dict[str, Optional[Dict]],
) -> None:
    """Table 3c — distribution of how many channels per event have ≥2 peaks."""
    print()
    print("=== Piled channels per event (channels with >=2 peaks) ===")

    col_w = max(12, max(len(f"{lab} (%)") + 2 for lab in run_labels))
    key_w = 8

    hdr = f"{'n_piled':>{key_w}}"
    for lab in run_labels:
        hdr += f"  {(lab + ' (%)'):>{col_w}}"
    print(hdr)
    print("-" * len(hdr))

    # Bins 0..4 are exact counts; bin 5 is "5+".
    bin_labels = ["0", "1", "2", "3", "4", "5+"]
    for b, blabel in enumerate(bin_labels):
        row = f"{blabel:>{key_w}}"
        for lab in run_labels:
            res = results.get(lab)
            if res is None:
                row += f"  {'—':>{col_w}}"
            else:
                n_ev = res["n_events"]
                cnt  = res["event_n_piled_channels_hist"][b]
                pct  = 100.0 * cnt / n_ev if n_ev else 0.0
                row += f"  {_fmt_pct(pct):>{col_w}}"
        print(row)


# ---------------------------------------------------------------------------
# Markdown output
# ---------------------------------------------------------------------------

def _md_row(cells: List[str]) -> str:
    return "| " + " | ".join(cells) + " |"


def _md_sep(n_cols: int) -> str:
    return "|" + "|".join(["---"] * n_cols) + "|"


def print_markdown_denom1(
    run_labels: List[str],
    results: Dict[str, Optional[Dict]],
    max_npeaks: int,
) -> None:
    print()
    print("### Denominator 1 markdown")
    print()

    header_cells = ["npeaks"]
    for lab in run_labels:
        header_cells += [f"{lab} (count)", f"{lab} (%)"]
    print(_md_row(header_cells))
    print(_md_sep(len(header_cells)))

    for k in range(max_npeaks + 1):
        label_k = f"{k}+" if k == max_npeaks else str(k)
        row_cells = [label_k]
        for lab in run_labels:
            res = results.get(lab)
            if res is None:
                row_cells += ["—", "—"]
            else:
                n_total = sum(res["counts"])
                cnt = res["counts"][k]
                pct = 100.0 * cnt / n_total if n_total else 0.0
                row_cells += [f"{cnt:,}", f"{pct:.1f}%"]
        print(_md_row(row_cells))

    # total row
    row_cells = ["**total**"]
    for lab in run_labels:
        res = results.get(lab)
        if res is None:
            row_cells += ["—", "—"]
        else:
            n_total = sum(res["counts"])
            row_cells += [f"**{n_total:,}**", "**100.0%**"]
    print(_md_row(row_cells))


def print_markdown_denom3(
    run_labels: List[str],
    results: Dict[str, Optional[Dict]],
    max_npeaks: int,
) -> None:
    print()
    print("### Denominator 3 markdown")
    print()

    header_cells = ["npeaks"]
    for lab in run_labels:
        header_cells += [f"{lab} (pulses)", f"{lab} (%)"]
    print(_md_row(header_cells))
    print(_md_sep(len(header_cells)))

    for k in range(1, max_npeaks + 1):
        label_k = f"{k}+" if k == max_npeaks else str(k)
        row_cells = [label_k]
        for lab in run_labels:
            res = results.get(lab)
            if res is None:
                row_cells += ["—", "—"]
            else:
                counts = res["counts"]
                total_pulses = sum(kk * counts[kk] for kk in range(1, max_npeaks + 1))
                pulses_k = k * counts[k]
                pct = 100.0 * pulses_k / total_pulses if total_pulses else 0.0
                row_cells += [f"{pulses_k:,}", f"{pct:.1f}%"]
        print(_md_row(row_cells))

    row_cells = ["**total**"]
    for lab in run_labels:
        res = results.get(lab)
        if res is None:
            row_cells += ["—", "—"]
        else:
            counts = res["counts"]
            total_pulses = sum(kk * counts[kk] for kk in range(1, max_npeaks + 1))
            row_cells += [f"**{total_pulses:,}**", "**100.0%**"]
    print(_md_row(row_cells))


def print_markdown_event_pileup_summary(
    run_labels: List[str],
    results: Dict[str, Optional[Dict]],
) -> None:
    """Markdown table — event-level pile-up fraction (ready to paste into tech note)."""
    print()
    print("### Event-level pile-up markdown")
    print()

    header_cells = ["Run", "events", "with pile-up", "fraction"]
    print(_md_row(header_cells))
    print(_md_sep(len(header_cells)))

    for lab in run_labels:
        res = results.get(lab)
        if res is None:
            print(_md_row([lab, "—", "—", "—"]))
            continue
        n_ev = res["n_events"]
        n_pu = res["n_events_with_pileup"]
        frac = 100.0 * n_pu / n_ev if n_ev else 0.0
        print(_md_row([lab, f"{n_ev:,}", f"{n_pu:,}", f"{frac:.1f}%"]))


# ---------------------------------------------------------------------------
# JSON output
# ---------------------------------------------------------------------------

def build_json_output(
    run_labels: List[str],
    results: Dict[str, Optional[Dict]],
    max_npeaks: int,
    max_events: int,
    material: str,
) -> Dict:
    runs_out: Dict = {}
    for lab in run_labels:
        res = results.get(lab)
        if res is None:
            runs_out[lab] = None
            continue
        counts = res["counts"]
        n_channels = sum(counts)
        total_pulses = sum(k * counts[k] for k in range(1, max_npeaks + 1))
        n_pileup_ch = sum(counts[k] for k in range(2, max_npeaks + 1))
        pileup_pulses = sum(k * counts[k] for k in range(2, max_npeaks + 1))

        cw: Dict[str, float] = {}
        pw: Dict[str, float] = {}
        for k in range(max_npeaks + 1):
            str_k = f"{k}+" if k == max_npeaks else str(k)
            cw[str_k] = round(100.0 * counts[k] / n_channels, 1) if n_channels else 0.0
        for k in range(1, max_npeaks + 1):
            str_k = f"{k}+" if k == max_npeaks else str(k)
            pulses_k = k * counts[k]
            pw[str_k] = round(100.0 * pulses_k / total_pulses, 1) if total_pulses else 0.0

        n_ev = res["n_events"]
        n_events_with_pileup = res["n_events_with_pileup"]
        event_max_npeaks_hist = res["event_max_npeaks_hist"]
        event_n_piled_channels_hist = res["event_n_piled_channels_hist"]

        runs_out[lab] = {
            "n_events": n_ev,
            "n_channels": n_channels,
            "n_pulses": total_pulses,
            "counts": {
                (f"{k}+" if k == max_npeaks else str(k)): counts[k]
                for k in range(max_npeaks + 1)
            },
            "channel_weighted": cw,
            "pulse_weighted": pw,
            "pileup_channel_weighted_pct": round(
                100.0 * n_pileup_ch / n_channels, 1) if n_channels else 0.0,
            "pileup_pulse_weighted_pct": round(
                100.0 * pileup_pulses / total_pulses, 1) if total_pulses else 0.0,
            "n_events_with_pileup": n_events_with_pileup,
            "event_pileup_fraction_pct": round(
                100.0 * n_events_with_pileup / n_ev, 1) if n_ev else 0.0,
            "event_max_npeaks_hist": {
                (f"{k}+" if k == max_npeaks else str(k)): event_max_npeaks_hist[k]
                for k in range(max_npeaks + 1)
            },
            "event_n_piled_channels_hist": {
                (str(b) if b < 5 else "5+"): event_n_piled_channels_hist[b]
                for b in range(6)
            },
        }

    return {
        "meta": {
            "max_events": max_events,
            "max_npeaks": max_npeaks,
            "material": material if material else None,
        },
        "runs": runs_out,
    }


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    ap = argparse.ArgumentParser(
        description="Peaks-per-channel distribution across EVIO runs.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    ap.add_argument(
        "evio_paths", nargs="+",
        help=(
            "One or more EVIO inputs.  Each may be a single split, a glob, "
            "or a directory.  If multiple, each is treated as a separate 'run' "
            "column labeled by its run number (sniffed from the filename "
            "prad_NNNNNN.evio.*)."
        ),
    )
    ap.add_argument(
        "--max-events", type=int, default=10000,
        help="Physics events to scan per run (0 = all).",
    )
    ap.add_argument(
        "--max-npeaks", type=int, default=8,
        help="Peaks are clamped to this max (last row is 'N+').",
    )
    ap.add_argument(
        "--daq-config", default="",
        help="DAQ config path.",
    )
    ap.add_argument(
        "--out", default="",
        help="Optional: write results to this JSON file.",
    )
    ap.add_argument(
        "--material", default="",
        help=(
            "Optional: restrict to one module type (PbWO4, PbGlass, LMS, Veto). "
            "Empty = all channels.  Requires --hc-map-file for module type lookup."
        ),
    )
    ap.add_argument(
        "--hc-map-file", default="",
        help="HyCal map (needed only if --material is set).",
    )
    args = ap.parse_args()

    # Deferred heavy import so --help works without prad2py bindings.
    try:
        from prad2py import dec, det
    except ImportError as exc:
        raise SystemExit(
            f"[ERROR] cannot import prad2py: {exc}\n"
            "        Build the python bindings (cmake -DBUILD_PYTHON=ON) and "
            "ensure the install directory is on PYTHONPATH."
        ) from exc

    material_filter: Optional[str] = args.material.strip() or None

    # Build one Pipeline to get a shared cfg + hycal (needed for material
    # filtering and for the WaveAnalyzer).  We use the first resolvable file.
    first_files: List[str] = []
    for raw in args.evio_paths:
        first_files = _resolve_evio_paths(raw)
        if first_files:
            break
    if not first_files:
        raise SystemExit("[ERROR] no EVIO files found for any of the given paths.")

    # We need det.PipelineBuilder to wire up cfg + hycal.
    b = det.PipelineBuilder()
    if args.daq_config:
        b.set_daq_config(args.daq_config)
    if args.hc_map_file:
        b.set_hycal_map(args.hc_map_file)
    b.set_run_number_from_evio(first_files[0])
    core = b.build()

    cfg      = core.daq_cfg
    hycal    = core.hycal if material_filter else None
    crate_map: Dict[int, int] = {
        int(r.tag): int(r.crate) for r in core.daq_cfg.roc_tags
    }

    # One shared WaveAnalyzer instance (stateless per call).
    import numpy as np  # noqa: F401 — imported here to keep top-level imports lean
    wa = dec.WaveAnalyzer(core.daq_cfg.wave_cfg)

    # -----------------------------------------------------------------------
    # Process each run
    # -----------------------------------------------------------------------

    # Maps run_label -> result dict (or None if all files failed).
    results: Dict[str, Optional[Dict]] = {}
    # Preserve insertion order for column ordering.
    run_labels: List[str] = []

    for raw in args.evio_paths:
        evio_files = _resolve_evio_paths(raw)
        if not evio_files:
            print(f"[WARN] no files found for input: {raw!r} — skipping", flush=True)
            continue

        label = _sniff_run_label(evio_files[0])
        # Deduplicate labels if the same run is passed twice.
        if label in results:
            suffix = 2
            base = label
            while f"{base}_{suffix}" in results:
                suffix += 1
            label = f"{base}_{suffix}"

        run_labels.append(label)
        print(f"\n[run {label}] {len(evio_files)} file(s)", flush=True)

        res = scan_run(
            label=label,
            evio_files=evio_files,
            dec=dec,
            wa=wa,
            cfg=cfg,
            max_events=args.max_events,
            max_npeaks=args.max_npeaks,
            material_filter=material_filter,
            hycal=hycal,
            crate_map=crate_map,
        )
        results[label] = res
        if res is None:
            print(f"[WARN] run {label}: all files failed to open — "
                  "column will show '—'", flush=True)

    if not run_labels:
        raise SystemExit("[ERROR] no runs could be processed.")

    # -----------------------------------------------------------------------
    # Output
    # -----------------------------------------------------------------------

    print_table_denom1(run_labels, results, args.max_npeaks)
    print_table_denom3(run_labels, results, args.max_npeaks)
    print_summary(run_labels, results, args.max_npeaks)
    print_table_event_pileup_summary(run_labels, results)
    print_table_event_max_npeaks(run_labels, results, args.max_npeaks)
    print_table_event_n_piled_channels(run_labels, results)

    print()
    print("--- Markdown ---")
    print_markdown_denom1(run_labels, results, args.max_npeaks)
    print_markdown_denom3(run_labels, results, args.max_npeaks)
    print_markdown_event_pileup_summary(run_labels, results)

    if args.out:
        payload = build_json_output(
            run_labels=run_labels,
            results=results,
            max_npeaks=args.max_npeaks,
            max_events=args.max_events,
            material=args.material,
        )
        out_path = Path(args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        with open(out_path, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)
        print(f"\n[output] JSON written to {out_path}", flush=True)

    return 0


if __name__ == "__main__":
    sys.exit(main())
