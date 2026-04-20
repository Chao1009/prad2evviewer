#!/usr/bin/env python3
"""
Waveform Histograms
===================
Parse an evio file with ``prad2py`` and accumulate per-channel peak
histograms identical to what the Waveform tab of the prad2 event
monitor shows:

    - integral histogram
    - height   histogram
    - position histogram  (peak time)

Peak finding is a faithful port of ``prad2dec/src/WaveAnalyzer.cpp``
(triangular-kernel smoothing, iterative pedestal with outlier
rejection, local-maxima search with baseline-subtracted threshold and
tail-cut integration).

Time cut (-t)
-------------
Without ``-t`` every peak above the height threshold fills every
histogram (height + integral + position).

With ``-t MIN[,MAX]`` only peaks inside the window fill the position
histogram, and the height & integral histograms are filled once per
event using the *best* peak in that window (highest integral) —
matching the server's Waveform-tab logic.

Output is a single JSON keyed by ``"{roc}_{slot}_{ch}"`` with the
binning parameters copied from ``database/config.json``.

Usage
-----
    python scripts/waveform_histograms.py RUN.evio.00000
    python scripts/waveform_histograms.py RUN.evio.00000 -t 170,190
    python scripts/waveform_histograms.py RUN.evio.00000 -t 150 -n 50000
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np


# ===========================================================================
#  prad2py discovery (mirrors tagger_viewer.py)
# ===========================================================================

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_DIR   = _SCRIPT_DIR.parent

for _cand in (
    _REPO_DIR / "build" / "python",
    _REPO_DIR / "build-release" / "python",
    _REPO_DIR / "build" / "Release" / "python",
):
    if _cand.is_dir() and str(_cand) not in sys.path:
        sys.path.insert(0, str(_cand))

try:
    import prad2py                                # type: ignore
    _HAVE_PRAD2PY = True
    _PRAD2PY_ERR  = ""
except Exception as _exc:
    prad2py = None                                # type: ignore
    _HAVE_PRAD2PY = False
    _PRAD2PY_ERR  = f"{type(_exc).__name__}: {_exc}"


# ===========================================================================
#  WaveAnalyzer — Python port of prad2dec/src/WaveAnalyzer.cpp
# ===========================================================================

@dataclass
class WaveConfig:
    resolution:      int   = 2
    threshold:       float = 5.0      # in pedestal-RMS units
    min_threshold:   float = 3.0      # absolute ADC floor above pedestal
    min_peak_ratio:  float = 0.3
    int_tail_ratio:  float = 0.1
    ped_nsamples:    int   = 30
    ped_flatness:    float = 1.0
    ped_max_iter:    int   = 3
    overflow:        int   = 4095
    clk_mhz:         float = 250.0
    max_peaks:       int   = 8


@dataclass
class Peak:
    pos:      int
    left:     int
    right:    int
    height:   float
    integral: float
    time:     float
    overflow: bool


def _smooth(raw: np.ndarray, res: int) -> np.ndarray:
    """Triangular-kernel smoothing with local edge renormalisation.
    Matches WaveAnalyzer::smooth() byte-for-byte."""
    n = raw.size
    if res <= 1:
        return raw.astype(np.float32)
    out = np.empty(n, dtype=np.float32)
    rf = float(res)
    for i in range(n):
        val  = float(raw[i])
        wsum = 1.0
        for j in range(1, res):
            if j > i or i + j >= n:
                continue
            w = 1.0 - j / (rf + 1.0)
            val  += w * (float(raw[i - j]) + float(raw[i + j]))
            wsum += 2.0 * w
        out[i] = val / wsum
    return out


def _pedestal(buf: np.ndarray, cfg: WaveConfig) -> Tuple[float, float]:
    n = min(cfg.ped_nsamples, buf.size)
    if n <= 0:
        return 0.0, 0.0
    s = buf[:n].astype(np.float64).copy()
    mean = float(s.mean())
    var  = float(s.var())
    rms  = np.sqrt(var) if var > 0 else 0.0
    for _ in range(cfg.ped_max_iter):
        cut = max(rms, cfg.ped_flatness)
        keep = np.abs(s - mean) < cut
        count = int(keep.sum())
        if count == s.size or count < 5:
            break
        s = s[keep]
        mean = float(s.mean())
        var  = float(s.var())
        rms  = np.sqrt(var) if var > 0 else 0.0
    return mean, rms


def _trend_sign(d: float) -> int:
    if abs(d) < 0.1:
        return 0
    return 1 if d > 0 else -1


def _find_peaks(raw: np.ndarray, buf: np.ndarray,
                ped_mean: float, ped_rms: float, thr: float,
                cfg: WaveConfig) -> List[Peak]:
    n = buf.size
    if n < 3:
        return []
    peaks: List[Peak] = []
    pk_range: List[Tuple[int, int]] = []

    i = 1
    while i < n - 1 and len(peaks) < cfg.max_peaks:
        tr1 = _trend_sign(float(buf[i]) - float(buf[i - 1]))
        tr2 = _trend_sign(float(buf[i]) - float(buf[i + 1]))
        if tr1 * tr2 < 0 or (tr1 == 0 and tr2 == 0):
            i += 1
            continue

        # walk through flat plateau on the right
        flat_end = i
        if tr2 == 0:
            while (flat_end < n - 1 and
                   _trend_sign(float(buf[flat_end]) - float(buf[flat_end + 1])) == 0):
                flat_end += 1
            if (flat_end >= n - 1 or
                _trend_sign(float(buf[flat_end]) - float(buf[flat_end + 1])) <= 0):
                i += 1
                continue
        peak_pos = (i + flat_end) // 2

        # expand rising/falling range
        left, right = i, flat_end
        while left > 0 and _trend_sign(float(buf[left]) - float(buf[left - 1])) > 0:
            left -= 1
        while (right < n - 1 and
               _trend_sign(float(buf[right]) - float(buf[right + 1])) >= 0):
            right += 1

        span = right - left
        if span <= 0:
            i += 1
            continue

        base = (float(buf[left])  * (right - peak_pos) +
                float(buf[right]) * (peak_pos - left)) / span
        smooth_height = float(buf[peak_pos]) - base
        if smooth_height < thr:
            i = right
            continue

        ped_height = float(buf[peak_pos]) - ped_mean
        if ped_height < thr or ped_height < 3.0 * ped_rms:
            i = right
            continue

        # integrate outward with tail cutoff
        integral = float(buf[peak_pos]) - ped_mean
        tail_cut = ped_height * cfg.int_tail_ratio
        int_left, int_right = peak_pos, peak_pos
        for j in range(peak_pos - 1, left - 1, -1):
            v = float(buf[j]) - ped_mean
            if v < tail_cut or v < ped_rms or v * ped_height < 0:
                int_left = j
                break
            integral += v
            int_left = j
        for j in range(peak_pos + 1, right + 1):
            v = float(buf[j]) - ped_mean
            if v < tail_cut or v < ped_rms or v * ped_height < 0:
                int_right = j
                break
            integral += v
            int_right = j

        # refine peak position on raw samples
        raw_pos = peak_pos
        raw_height = float(raw[peak_pos]) - ped_mean
        search = max(1, cfg.resolution) + (flat_end - i) // 2
        for j in range(1, search + 1):
            if peak_pos - j >= 0:
                h = float(raw[peak_pos - j]) - ped_mean
                if h > raw_height:
                    raw_height = h
                    raw_pos = peak_pos - j
            if peak_pos + j < n:
                h = float(raw[peak_pos + j]) - ped_mean
                if h > raw_height:
                    raw_height = h
                    raw_pos = peak_pos + j

        # secondary-peak overlap rejection
        rejected = False
        for k, (lk, rk) in enumerate(pk_range):
            if left <= rk and right >= lk:
                if smooth_height < peaks[k].height * cfg.min_peak_ratio:
                    rejected = True
                    break
        if rejected:
            i = right
            continue

        peaks.append(Peak(
            pos      = raw_pos,
            left     = int_left,
            right    = int_right,
            height   = raw_height,
            integral = integral,
            time     = raw_pos * 1e3 / cfg.clk_mhz,
            overflow = raw[raw_pos] >= cfg.overflow,
        ))
        pk_range.append((left, right))
        i = right

    return peaks


def analyze(samples: np.ndarray,
            cfg: WaveConfig) -> Tuple[float, float, List[Peak]]:
    """Return (ped_mean, ped_rms, [Peak, ...]). Mirrors WaveAnalyzer::Analyze."""
    n = samples.size
    if n <= 0:
        return 0.0, 0.0, []
    buf = _smooth(samples, cfg.resolution)
    ped_mean, ped_rms = _pedestal(buf, cfg)
    thr = max(cfg.threshold * ped_rms, cfg.min_threshold)
    peaks = _find_peaks(samples, buf, ped_mean, ped_rms, thr, cfg)
    return ped_mean, ped_rms, peaks


# ===========================================================================
#  Histogram container
# ===========================================================================

@dataclass
class Hist1D:
    nbins:   int
    bins:    np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    under:   int        = 0
    over:    int        = 0

    def __post_init__(self):
        if self.bins.size == 0:
            self.bins = np.zeros(self.nbins, dtype=np.int64)

    def fill(self, v: float, bmin: float, bstep: float):
        if v < bmin:
            self.under += 1
            return
        b = int((v - bmin) / bstep)
        if b >= self.nbins:
            self.over += 1
            return
        self.bins[b] += 1

    def to_json(self) -> Dict:
        return {
            "bins":      self.bins.tolist(),
            "underflow": int(self.under),
            "overflow":  int(self.over),
        }


@dataclass
class ChannelHists:
    roc:          int
    slot:         int
    channel:      int
    events:       int       = 0
    peak_events:  int       = 0
    tcut_events:  int       = 0
    height:       Hist1D    = None
    integral:     Hist1D    = None
    position:     Hist1D    = None


# ===========================================================================
#  Mapping helpers
# ===========================================================================

def load_daq_map(path: Path) -> Dict[Tuple[int, int, int], str]:
    """Returns (crate_idx, slot, ch) -> module name."""
    with open(path, encoding="utf-8") as f:
        entries = json.load(f)
    out = {}
    for e in entries:
        out[(int(e["crate"]), int(e["slot"]), int(e["channel"]))] = e["name"]
    return out


def load_roc_tag_map(path: Path) -> Dict[int, int]:
    """Returns roc_tag -> crate_idx (from daq_config.roc_tags, FADC ROCs only)."""
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    out = {}
    for r in cfg.get("roc_tags", []):
        tag = int(r["tag"], 16) if isinstance(r["tag"], str) else int(r["tag"])
        if r.get("type") == "roc":
            out[tag] = int(r["crate"])
    return out


def load_hist_config(path: Path) -> Dict:
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    return cfg.get("waveform", {})


def load_trigger_bit_map(path: Path) -> Dict[str, int]:
    """Returns name -> bit index."""
    if not path.is_file():
        return {}
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    return {e["name"]: int(e["bit"]) for e in cfg.get("trigger_bits", [])}


# ===========================================================================
#  Main driver
# ===========================================================================

def _mask_from_names(names: List[str], bitmap: Dict[str, int]) -> int:
    m = 0
    for n in names:
        if n in bitmap:
            m |= (1 << bitmap[n])
        else:
            print(f"  warning: trigger bit {n!r} not in trigger_bits.json", file=sys.stderr)
    return m


def run(evio_path: str,
        *,
        output_path: Path,
        hist_config: Dict,
        daq_map: Dict[Tuple[int, int, int], str],
        roc_to_crate: Dict[int, int],
        accept_mask: int,
        reject_mask: int,
        max_events: int,
        wcfg: WaveConfig,
        time_cut: Optional[Tuple[float, float]],
        daq_config_path: str = "",
        progress_every: int = 1000) -> None:
    if not _HAVE_PRAD2PY:
        print(f"ERROR: prad2py not importable ({_PRAD2PY_ERR})", file=sys.stderr)
        print("Build with: cmake -DBUILD_PYTHON=ON -S . -B build && cmake --build build",
              file=sys.stderr)
        sys.exit(1)

    # derive binning
    h_cfg  = hist_config.get("height_hist",   {"min": 0, "max": 4000,  "step": 10})
    i_cfg  = hist_config.get("integral_hist", {"min": 0, "max": 20000, "step": 100})
    p_cfg  = hist_config.get("time_hist",     {"min": 0, "max": 400,   "step": 4})
    thr_cfg = hist_config.get("thresholds",   {})
    hist_threshold      = float(thr_cfg.get("min_peak_height", 10.0))
    wcfg.min_peak_ratio = float(thr_cfg.get("min_secondary_peak_ratio",
                                             wcfg.min_peak_ratio))
    if time_cut is not None:
        t_min, t_max = time_cut
    else:
        t_min, t_max = None, None

    def nbins(cfg):
        span = cfg["max"] - cfg["min"]
        return max(1, int(np.ceil(span / cfg["step"])))
    h_nbins, i_nbins, p_nbins = nbins(h_cfg), nbins(i_cfg), nbins(p_cfg)

    channels: Dict[Tuple[int, int, int], ChannelHists] = {}
    total_events = 0
    t0 = time.time()

    dec = prad2py.dec
    cfg = dec.load_daq_config(daq_config_path) if daq_config_path else dec.load_daq_config()
    ch  = dec.EvChannel()
    ch.set_config(cfg)
    st = ch.open(evio_path)
    if st != dec.Status.success:
        print(f"ERROR: cannot open {evio_path}: {st}", file=sys.stderr)
        sys.exit(1)

    print(f"Reading {evio_path} …")

    stop = False
    while ch.read() == dec.Status.success:
        if not ch.scan():
            continue
        if ch.get_event_type() != dec.EventType.Physics:
            continue

        for i in range(ch.get_n_events()):
            ch.select_event(i)
            info = ch.info()
            tb = int(info.trigger_bits)
            if accept_mask and (tb & accept_mask) == 0:
                continue
            if reject_mask and (tb & reject_mask):
                continue
            fadc_evt = ch.fadc()

            for r in range(fadc_evt.nrocs):
                roc = fadc_evt.roc(r)
                roc_tag = int(roc.tag)
                if roc_tag not in roc_to_crate:
                    continue  # non-FADC ROC (tagger, GEM, TI-slave)
                crate = roc_to_crate[roc_tag]
                for s in roc.present_slots():
                    slot = roc.slot(s)
                    for c in slot.present_channels():
                        samples = slot.channel(c).samples   # uint16 numpy copy
                        if samples.size < 10:
                            continue                         # ADC1881M single-sample
                        key = (roc_tag, s, c)
                        hits = channels.get(key)
                        if hits is None:
                            hits = ChannelHists(
                                roc=roc_tag, slot=s, channel=c,
                                height  =Hist1D(h_nbins),
                                integral=Hist1D(i_nbins),
                                position=Hist1D(p_nbins),
                            )
                            channels[key] = hits
                        hits.events += 1

                        _, _, peaks = analyze(samples, wcfg)

                        if time_cut is None:
                            # no window — every qualifying peak fills every hist
                            any_peak = False
                            for p in peaks:
                                if p.height < hist_threshold:
                                    continue
                                any_peak = True
                                hits.position.fill(p.time,     p_cfg["min"], p_cfg["step"])
                                hits.integral.fill(p.integral, i_cfg["min"], i_cfg["step"])
                                hits.height.fill  (p.height,   h_cfg["min"], h_cfg["step"])
                            if any_peak:
                                hits.peak_events += 1
                        else:
                            # best peak in window for height/integral;
                            # position hist only gets peaks inside the window
                            best_int = -1.0
                            best_hgt = -1.0
                            any_peak = False
                            for p in peaks:
                                if p.height < hist_threshold:
                                    continue
                                any_peak = True
                                if t_min <= p.time <= t_max:
                                    hits.position.fill(p.time, p_cfg["min"], p_cfg["step"])
                                    if p.integral > best_int:
                                        best_int = p.integral
                                        best_hgt = p.height
                            if any_peak:
                                hits.peak_events += 1
                            if best_int >= 0:
                                hits.tcut_events += 1
                                hits.integral.fill(best_int, i_cfg["min"], i_cfg["step"])
                                hits.height.fill  (best_hgt, h_cfg["min"], h_cfg["step"])

            total_events += 1
            if max_events and total_events >= max_events:
                stop = True
                break

            if total_events % progress_every == 0:
                dt = time.time() - t0
                rate = total_events / dt if dt > 0 else 0.0
                print(f"  {total_events:>8d} events  "
                      f"{len(channels):>4d} channels  "
                      f"{rate:6.1f} ev/s", flush=True)
        if stop:
            break

    ch.close()
    dt = time.time() - t0
    rate = total_events / dt if dt > 0 else 0.0
    print(f"Done: {total_events} events  in  {dt:.1f} s  ({rate:.1f} ev/s)")

    # ── dump ────────────────────────────────────────────────────────────
    out = {
        "source_file":  str(evio_path),
        "total_events": total_events,
        "time_cut":     (None if time_cut is None
                         else {"min": t_min, "max": t_max}),
        "height_hist":  {"min": h_cfg["min"], "max": h_cfg["max"], "step": h_cfg["step"], "nbins": h_nbins},
        "integral_hist":{"min": i_cfg["min"], "max": i_cfg["max"], "step": i_cfg["step"], "nbins": i_nbins},
        "position_hist":{"min": p_cfg["min"], "max": p_cfg["max"], "step": p_cfg["step"], "nbins": p_nbins},
        "threshold":    hist_threshold,
        "wave_config":  wcfg.__dict__,
        "channels":     {},
    }
    for (roc_tag, s, c), hh in sorted(channels.items()):
        key = f"{roc_tag}_{s}_{c}"
        crate = roc_to_crate.get(roc_tag, -1)
        out["channels"][key] = {
            "module":        daq_map.get((crate, s, c)),
            "roc":           roc_tag,
            "slot":          s,
            "channel":       c,
            "events":        hh.events,
            "peak_events":   hh.peak_events,
            "tcut_events":   hh.tcut_events,
            "height_hist":   hh.height.to_json(),
            "integral_hist": hh.integral.to_json(),
            "position_hist": hh.position.to_json(),
        }
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w") as f:
        json.dump(out, f)
    print(f"Wrote {output_path}  ({len(channels)} channels)")


def main():
    ap = argparse.ArgumentParser(
        description="Parse an evio file and accumulate waveform-tab peak histograms.")
    ap.add_argument("evio_file",
                    help="Path to prad_*.evio.* file")
    ap.add_argument("-o", "--output", type=Path, default=None,
                    help="Output JSON path (default: <evio_basename>_hist.json)")
    ap.add_argument("-n", "--max-events", type=int, default=0,
                    help="Stop after this many physics events (0 = all)")
    ap.add_argument("-t", "--time-cut", type=str, default=None, metavar="MIN[,MAX]",
                    help="Time-cut window (ns). Without -t every peak fills every "
                         "histogram. With -t MIN,MAX height & integral hists take "
                         "the best peak (by integral) inside the window, and the "
                         "position hist keeps only peaks inside it. -t MIN is "
                         "treated as MIN,+infinity.")
    ap.add_argument("--config", type=Path,
                    default=_REPO_DIR / "database" / "config.json",
                    help="Main config.json (for waveform binning)")
    ap.add_argument("--daq-config", type=Path,
                    default=_REPO_DIR / "database" / "daq_config.json",
                    help="daq_config.json (for ROC-tag → crate mapping)")
    ap.add_argument("--daq-map", type=Path,
                    default=_REPO_DIR / "database" / "daq_map.json",
                    help="daq_map.json (for module-name lookup)")
    ap.add_argument("--trigger-bits", type=Path,
                    default=_REPO_DIR / "database" / "trigger_bits.json",
                    help="trigger_bits.json (for --accept/--reject-trigger names)")
    ap.add_argument("--accept-trigger", action="append", default=[],
                    metavar="NAME",
                    help="Require at least one of these trigger bits "
                         "(repeatable). Default: none (accept any).")
    ap.add_argument("--reject-trigger", action="append", default=None,
                    metavar="NAME",
                    help="Drop events with any of these trigger bits "
                         "(repeatable). Default: uses config.json setting.")
    ap.add_argument("--progress-every", type=int, default=1000,
                    help="Print progress every N events (default 1000)")
    args = ap.parse_args()

    evio_path = args.evio_file
    out_path = args.output or Path(Path(evio_path).name + "_hist.json")

    time_cut: Optional[Tuple[float, float]] = None
    if args.time_cut is not None:
        parts = [p.strip() for p in args.time_cut.split(",") if p.strip()]
        try:
            if len(parts) == 1:
                time_cut = (float(parts[0]), float("inf"))
            elif len(parts) == 2:
                time_cut = (float(parts[0]), float(parts[1]))
            else:
                raise ValueError
        except ValueError:
            ap.error(f"--time-cut expects MIN or MIN,MAX (got {args.time_cut!r})")
        if time_cut[0] >= time_cut[1]:
            ap.error(f"--time-cut MIN must be < MAX (got {time_cut})")

    hist_cfg     = load_hist_config(args.config)
    roc_to_crate = load_roc_tag_map(args.daq_config)
    daq_map      = load_daq_map(args.daq_map) if args.daq_map.is_file() else {}
    bit_map      = load_trigger_bit_map(args.trigger_bits)

    accept_names = args.accept_trigger or hist_cfg.get("accept_trigger_bits", []) or []
    if args.reject_trigger is None:
        reject_names = hist_cfg.get("reject_trigger_bits", []) or []
    else:
        reject_names = args.reject_trigger
    accept_mask = _mask_from_names(accept_names, bit_map) if accept_names else 0
    reject_mask = _mask_from_names(reject_names, bit_map) if reject_names else 0

    print(f"Accept trigger bits: {accept_names or '(any)'}  mask=0x{accept_mask:08x}")
    print(f"Reject trigger bits: {reject_names or '(none)'} mask=0x{reject_mask:08x}")
    if time_cut is None:
        print("Time cut: (none — every peak fills every histogram)")
    else:
        lo, hi = time_cut
        hi_s = "+inf" if hi == float("inf") else f"{hi:g}"
        print(f"Time cut: [{lo:g}, {hi_s}] ns")

    run(evio_path,
        output_path  = out_path,
        hist_config  = hist_cfg,
        daq_map      = daq_map,
        roc_to_crate = roc_to_crate,
        accept_mask  = accept_mask,
        reject_mask  = reject_mask,
        max_events   = args.max_events,
        wcfg         = WaveConfig(),
        time_cut        = time_cut,
        daq_config_path = str(args.daq_config) if args.daq_config.is_file() else "",
        progress_every  = args.progress_every)


if __name__ == "__main__":
    main()
