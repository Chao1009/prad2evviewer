#!/usr/bin/env python3
"""
HyCal Event Viewer
==================
Browses an evio file event-by-event.  Opens the file in evio's
random-access mode (no event processing at open time), indexes physics
sub-events, and lets the user step through them.  Two tabs:

* **Waveform** — per-module FADC display, stacked waveform, four
  accumulating histograms (peak height, integral, time, n-peaks).
  "Process next 10k" fills histograms in a background pass.
* **Cluster** — HyCal heatmap of per-module energy with cluster
  overlays (crosshair + energy label), cluster table, selector.
  Clustering uses ``prad2py.det.HyCalCluster`` on live ADC data;
  the per-run calibration is loaded by ``prad2py.det.PipelineBuilder``.

Usage
-----
    python scripts/hycal_event_viewer.py RUN.evio.00000
    python scripts/hycal_event_viewer.py             # File → Open…
"""
from __future__ import annotations

import argparse
import json
import re
import sys
import traceback
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, NamedTuple, Optional, Tuple

import numpy as np

from PyQt6.QtWidgets import (
    QAbstractItemView, QApplication, QMainWindow, QWidget,
    QVBoxLayout, QHBoxLayout, QFormLayout, QGridLayout,
    QLabel, QComboBox, QCheckBox, QCompleter, QFileDialog, QMessageBox,
    QProgressDialog, QSizePolicy, QStatusBar, QToolTip, QPushButton,
    QSpinBox, QSplitter, QTabWidget, QTableWidget,
    QTableWidgetItem, QHeaderView, QDockWidget, QGroupBox,
    QDialog, QDialogButtonBox, QLineEdit,
)
from PyQt6.QtCore import (
    Qt, QObject, QPointF, QRectF, QThread, pyqtSignal, QTimer,
)
from PyQt6.QtGui import (
    QAction, QBrush, QKeySequence, QPainter, QColor, QPen, QFont, QPolygonF,
    QShortcut, QDoubleValidator,
)

from hycal_geoview import (
    load_modules as load_geo_modules, load_daq_map, load_roc_tag_map,
    HyCalMapWidget, ColorRangeController, cmap_qcolor, series_qcolor,
    draw_wave_axes, AUX_TYPES, OVERLAY_BUTTON_QSS,
    apply_theme_palette, set_theme,
    available_themes, THEME, themed,
    setup_tuning_dock, add_config_rows, editor_value, set_editor_value,
    config_to_editors, editors_to_config, start_worker_thread,
)
from evio_io import EvioCursor, iter_physics_records, open_evio
from prad2_env import import_prad2py


# ---- prad2py discovery ----

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_DIR   = _SCRIPT_DIR.parent

prad2py, _PRAD2PY_ERR = import_prad2py()
_HAVE_PRAD2PY = prad2py is not None


def _check_evchannel_support() -> Optional[str]:
    """Return None if the loaded prad2py has the expected EvChannel API,
    else an error string suitable for showing to the user."""
    if not _HAVE_PRAD2PY:
        return _PRAD2PY_ERR
    try:
        ch = prad2py.dec.EvChannel()
        if not hasattr(ch, "open_auto"):
            return ("prad2py is missing open_auto — rebuild prad2py after "
                    "the EvChannel.cpp changes:\n"
                    "  cmake -DBUILD_PYTHON=ON -S . -B build && "
                    "cmake --build build --target prad2py")
    except Exception as e:
        return f"{type(e).__name__}: {e}"
    return None


# ---- WaveAnalyzer ----
# Direct exports of the C++ implementation in prad2py.dec (the server runs
# the same code).

WaveConfig = prad2py.dec.WaveConfig if _HAVE_PRAD2PY else None
Peak       = prad2py.dec.Peak       if _HAVE_PRAD2PY else None


def analyze(samples, cfg):
    """Run the C++ WaveAnalyzer on one channel's samples.
    Returns ``(ped_mean, ped_rms, peaks_list)``."""
    return prad2py.dec.WaveAnalyzer(cfg).analyze(samples)


# ---- Histogram accumulator ----

@dataclass
class Hist1D:
    nbins:   int
    bmin:    float = 0.0
    bstep:   float = 1.0
    bins:    np.ndarray = field(default_factory=lambda: np.zeros(0, np.int64))
    under:   int = 0
    over:    int = 0

    def __post_init__(self):
        if self.bins.size == 0:
            self.bins = np.zeros(self.nbins, dtype=np.int64)

    def fill(self, v: float):
        if v < self.bmin:
            self.under += 1
            return
        b = int((v - self.bmin) / self.bstep)
        if b >= self.nbins:
            self.over += 1
            return
        self.bins[b] += 1

    def reset(self):
        self.bins[:] = 0
        self.under = 0
        self.over = 0

    def to_json(self) -> Dict:
        return {"bins": self.bins.tolist(),
                "underflow": int(self.under), "overflow": int(self.over)}


@dataclass
class ChannelHists:
    roc:         int
    slot:        int
    channel:     int
    module:      Optional[str]
    events:      int = 0
    peak_events: int = 0
    height:      Optional[Hist1D] = None
    integral:    Optional[Hist1D] = None
    position:    Optional[Hist1D] = None
    npeaks:      Optional[Hist1D] = None

    def fill_peaks(self, peaks, flt: "WaveformFilter") -> int:
        """Fold one event's peaks that pass ``flt`` into the histograms and
        return how many passed."""
        kept = 0
        for p in peaks:
            if not flt.passes(p):
                continue
            self.height.fill(p.height)
            self.integral.fill(p.integral)
            self.position.fill(p.time)
            kept += 1
        self.npeaks.fill(kept)
        self.events += 1
        if kept > 0:
            self.peak_events += 1
        return kept


class _HistSpec(NamedTuple):
    attr: str                 # ChannelHists field
    cfg_key: Optional[str]    # binning key in monitor_config.json's "waveform"
    json_key: str             # key in the saved histogram JSON
    label: str                # placeholder title
    title: str                # plot title
    color: str
    default: Dict             # binning when cfg_key is absent


# The n-peaks binning is fixed; its left edge at -0.5 puts integer counts
# (0, 1, 2 …) on bin centres rather than at the left edge of each bar.
_HIST_SPECS: Tuple[_HistSpec, ...] = (
    _HistSpec("height", "height_hist", "height_hist", "Peak Height",
              "Peak Height [ADC]", "#e599f7",
              {"min": 0, "max": 4000, "step": 10}),
    _HistSpec("integral", "integral_hist", "integral_hist", "Peak Integral",
              "Peak Integral [ADC·sample]", "#00b4d8",
              {"min": 0, "max": 20000, "step": 100}),
    _HistSpec("position", "time_hist", "position_hist", "Peak Time",
              "Peak Time [ns]", "#51cf66",
              {"min": 0, "max": 400, "step": 4}),
    _HistSpec("npeaks", None, "npeaks_hist", "Peaks / Event",
              "Peaks / Event", "#ffa657",
              {"min": -0.5, "max": 10.5, "step": 1}),
)


# ---- Waveform peak filter ----

# Mirrors AppState::peak_quality_bits_def in src/app_state_init.cpp; the
# masks are read from prad2py.dec.Q_PEAK_* at runtime so they always agree
# with the C++ side.
def _resolve_quality_bits() -> List[Dict[str, object]]:
    if not _HAVE_PRAD2PY:
        return []
    out: List[Dict[str, object]] = []
    for name, label in (("PILED", "Pile-up"), ("DECONVOLVED", "Deconvolved")):
        bit = getattr(prad2py.dec, f"Q_PEAK_{name}", None)
        if bit is None:
            continue
        # Q_PEAK_* are masks (1 << bit_index).  Store the mask directly —
        # the filter compares peak.quality & mask, no shift needed.
        out.append({"mask": int(bit), "name": name, "label": label})
    return out


PEAK_QUALITY_BITS: List[Dict[str, object]] = _resolve_quality_bits()


@dataclass
class WaveformFilter:
    """Per-peak filter for the Waveform tab — replicates ``PeakFilter`` from
    ``src/app_state.h`` so monitor_config.json's ``waveform.filter`` JSON
    drives both the Python viewer and the live web monitor identically.

    A ``None`` bound means "no constraint" (matches the web monitor's empty
    input fields).  ``enable=False`` short-circuits ``passes()`` so the user
    can toggle filtering off without losing their tuned values.
    """
    enable:       bool = True
    time_min:     Optional[float] = None
    time_max:     Optional[float] = None
    integral_min: Optional[float] = None
    integral_max: Optional[float] = None
    height_min:   Optional[float] = None
    height_max:   Optional[float] = None
    q_accept:     int = 0   # bitmask, 0 = accept any
    q_reject:     int = 0   # bitmask, 0 = reject none

    def passes(self, pk) -> bool:
        if not self.enable:
            return True
        if self.time_min     is not None and pk.time     < self.time_min:     return False
        if self.time_max     is not None and pk.time     > self.time_max:     return False
        if self.integral_min is not None and pk.integral < self.integral_min: return False
        if self.integral_max is not None and pk.integral > self.integral_max: return False
        if self.height_min   is not None and pk.height   < self.height_min:   return False
        if self.height_max   is not None and pk.height   > self.height_max:   return False
        q = int(getattr(pk, "quality", 0))
        if self.q_reject and (q & self.q_reject):    return False
        if self.q_accept and not (q & self.q_accept): return False
        return True

    def copy(self) -> "WaveformFilter":
        return WaveformFilter(**self.__dict__)

    @classmethod
    def from_json(cls, j: Optional[Dict]) -> "WaveformFilter":
        """Parse a ``waveform.filter`` JSON object into a filter.  Mirrors
        ``PeakFilter::parse``; unknown keys are ignored so future schema
        additions don't break loading."""
        f = cls()
        j = j or {}

        def _axis(key: str):
            ax = j.get(key) or {}
            lo = ax.get("min")
            hi = ax.get("max")
            return (float(lo) if lo is not None else None,
                    float(hi) if hi is not None else None)

        f.time_min,     f.time_max     = _axis("time")
        f.integral_min, f.integral_max = _axis("integral")
        f.height_min,   f.height_max   = _axis("height")

        qb = j.get("quality_bits") or {}
        f.q_accept = _names_to_mask(qb.get("accept") or [])
        f.q_reject = _names_to_mask(qb.get("reject") or [])
        return f

    def to_json(self) -> Dict:
        out: Dict = {}

        def _axis(lo, hi):
            ax: Dict = {}
            if lo is not None: ax["min"] = lo
            if hi is not None: ax["max"] = hi
            return ax

        if (a := _axis(self.time_min,     self.time_max))     : out["time"]     = a
        if (a := _axis(self.integral_min, self.integral_max)) : out["integral"] = a
        if (a := _axis(self.height_min,   self.height_max))   : out["height"]   = a
        if self.q_accept or self.q_reject:
            out["quality_bits"] = {
                "accept": _mask_to_names(self.q_accept),
                "reject": _mask_to_names(self.q_reject),
            }
        return out


def _names_to_mask(names) -> int:
    m = 0
    for n in names or ():
        for d in PEAK_QUALITY_BITS:
            if d["name"] == n:
                m |= int(d["mask"])
                break
    return m


def _mask_to_names(mask: int) -> List[str]:
    if not mask:
        return []
    return [d["name"] for d in PEAK_QUALITY_BITS if mask & int(d["mask"])]


# ---- Config loaders ----

def load_hist_config(path: Path) -> Dict:
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    return cfg.get("waveform", {})


def load_trigger_bit_map(path: Path) -> Dict[str, int]:
    if not path.is_file():
        return {}
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    return {e["name"]: int(e["bit"]) for e in cfg.get("trigger_bits", [])}


def _mask_from_names(names: List[str], bitmap: Dict[str, int]) -> int:
    m = 0
    for n in names:
        if n in bitmap:
            m |= (1 << bitmap[n])
        else:
            print(f"  warning: trigger bit {n!r} not in trigger_bits.json",
                  file=sys.stderr)
    return m


def trigger_ok(tb: int, accept: int, reject: int) -> bool:
    """True if ``tb`` has one of the ``accept`` bits (any, when 0) and none
    of the ``reject`` bits."""
    return not ((accept and (tb & accept) == 0) or (reject and (tb & reject)))


# ---- Histogram filling ----

def _nbins(c: Dict) -> int:
    span = c["max"] - c["min"]
    return max(1, int(np.ceil(span / c["step"])))


def _make_hists(cfgs: Dict[str, Dict], roc: int, slot: int, channel: int,
                module: Optional[str]) -> ChannelHists:
    """ChannelHists with one Hist1D per ``cfgs`` entry (field -> binning)."""
    return ChannelHists(
        roc=roc, slot=slot, channel=channel, module=module,
        **{a: Hist1D(_nbins(c), c["min"], c["step"]) for a, c in cfgs.items()})


def _iter_channel_peaks(fadc_evt,
                        channels: Dict[Tuple[int, int, int], ChannelHists],
                        wcfg):
    """Yield ``(key, hits, peaks)`` for every channel of the event that has
    histograms in ``channels`` and at least 10 samples."""
    for r in range(fadc_evt.nrocs):
        roc = fadc_evt.roc(r)
        roc_tag = int(roc.tag)
        for s in roc.present_slots():
            slot = roc.slot(s)
            for c in slot.present_channels():
                key = (roc_tag, s, c)
                hits = channels.get(key)
                if hits is None:
                    continue
                samples = slot.channel(c).samples
                if samples.size < 10:
                    continue
                _, _, peaks = analyze(samples, wcfg)
                yield key, hits, peaks


# One flag per physics sub-event, set once its peaks are in the histograms
# so re-visiting an event does not count it twice.
def _is_folded(acc: Optional[np.ndarray], idx: int) -> bool:
    return acc is not None and 0 <= idx < acc.size and bool(acc[idx])


def _mark_folded(acc: Optional[np.ndarray], idx: int) -> None:
    if acc is not None and 0 <= idx < acc.size:
        acc[idx] = True


# ---- Indexer — background pass to locate all physics sub-events ----

class IndexerWorker(QObject):
    """Scans the file once to record (evio_idx, sub_idx) per physics
    sub-event.  No waveform decoding — Scan() only."""

    progressed = pyqtSignal(int, int)   # (evio_events_scanned, total_evio_events)
    finished   = pyqtSignal(object)     # {"index": list, "total_evio": int, "cancelled": bool}
    failed     = pyqtSignal(str)

    def __init__(self, evio_path: str, daq_config_path: str):
        super().__init__()
        self._path = evio_path
        self._daq_cfg_path = daq_config_path
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        try:
            self.finished.emit(self._run())
        except Exception as e:
            self.failed.emit(f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    def _run(self) -> Dict:
        ch, is_ra = open_evio(self._path, self._daq_cfg_path)

        # In RA mode we know the total upfront; sequential mode walks to EOF
        # so we just report a rolling count.
        total_evio = ch.get_random_access_event_count() if is_ra else 0
        index: List[Tuple[int, int]] = []

        progress_every = max(1, total_evio // 200) if total_evio else 500
        n_rec = 0

        def _on_record(idx: int):
            nonlocal n_rec
            n_rec = idx + 1
            if (n_rec % progress_every) == 0:
                self.progressed.emit(n_rec, total_evio or n_rec)

        try:
            for ei in iter_physics_records(ch, is_ra, lambda: self._cancel,
                                           _on_record):
                index.extend((ei, si) for si in range(ch.get_n_events()))
        finally:
            ch.close()
        self.progressed.emit(n_rec, total_evio or n_rec)
        return {"index": index,
                "total_evio": n_rec, "cancelled": self._cancel}


# ---- Batch processor — fills all-module hists for the next N events ----

class BatchWorker(QObject):
    """Reads events start_idx .. start_idx + n - 1 (no display updates) and
    fills histograms for every channel present.  Runs in its own thread with
    its own EvChannel handle (separate from the UI's)."""

    progressed = pyqtSignal(int, int, int)  # (done, target, peaks_found)
    finished   = pyqtSignal(int)            # events_processed
    failed     = pyqtSignal(str)

    def __init__(self, evio_path: str, daq_config_path: str,
                 index: List[Tuple[int, int]],
                 start_idx: int, count: int,
                 channels: Dict[Tuple[int, int, int], ChannelHists],
                 wcfg: WaveConfig, peak_filter: "WaveformFilter",
                 accept_mask: int, reject_mask: int,
                 accumulated: Optional[np.ndarray] = None):
        super().__init__()
        self._path = evio_path
        self._daq_cfg_path = daq_config_path
        self._index = index
        self._start = start_idx
        self._count = count
        self._channels = channels
        self._wcfg = wcfg
        # Snapshot the filter so user edits during the batch don't change
        # the cuts mid-run (each batch should fill against a stable cut).
        self._filter = peak_filter.copy()
        self._accept = accept_mask
        self._reject = reject_mask
        self._accumulated = accumulated     # shared bool array, mutated in place
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        try:
            self.finished.emit(self._run())
        except Exception as e:
            self.failed.emit(f"{type(e).__name__}: {e}\n{traceback.format_exc()}")

    def _run(self) -> int:
        cur = EvioCursor(self._path, self._daq_cfg_path)
        n_done = 0
        peaks_found = 0
        progress_every = max(1, self._count // 100)
        wcfg = self._wcfg
        flt = self._filter
        channels = self._channels

        def _fold_event(phys_idx: int, sub_idx: int) -> int:
            """Scan + select the loaded event, accumulate hists for every
            channel, return peaks found this event (or -1 if the event
            is rejected by trigger mask or dedup)."""
            nonlocal n_done
            if _is_folded(self._accumulated, phys_idx):
                n_done += 1
                return -1
            ch = cur.ch
            if not ch.scan():
                return -1
            ch.select_event(sub_idx)
            if not trigger_ok(int(ch.info().trigger_bits),
                              self._accept, self._reject):
                n_done += 1
                return -1
            pfound = 0
            for _, hits, peaks in _iter_channel_peaks(ch.fadc(), channels, wcfg):
                pfound += hits.fill_peaks(peaks, flt)
            _mark_folded(self._accumulated, phys_idx)
            n_done += 1
            return pfound

        # self._index is in file order, so a sequential cursor only ever
        # walks forward.
        try:
            for i in range(self._count):
                if self._cancel: break
                phys_idx = self._start + i
                if phys_idx >= len(self._index): break
                ev_idx, sub_idx = self._index[phys_idx]
                try:
                    cur.seek(ev_idx)
                except RuntimeError:
                    if cur.is_ra:
                        continue
                    raise
                pf = _fold_event(phys_idx, sub_idx)
                if pf > 0: peaks_found += pf
                if (i % progress_every) == 0:
                    self.progressed.emit(n_done, self._count, peaks_found)
        finally:
            cur.close()

        self.progressed.emit(n_done, self._count, peaks_found)
        return n_done


def _find_channel_samples(fadc_evt, roc_tag: int, slot: int, channel: int):
    """Return samples array for (roc, slot, ch), or None if not present."""
    for r in range(fadc_evt.nrocs):
        roc = fadc_evt.roc(r)
        if int(roc.tag) != roc_tag:
            continue
        if slot not in roc.present_slots():
            continue
        slot_data = roc.slot(slot)
        if channel not in slot_data.present_channels():
            continue
        return slot_data.channel(channel).samples
    return None


# ---- Plot widget helpers (overlay controls, framed canvas base) ----

def _overlay_checkbox_qss() -> str:
    """QSS for a compact checkbox drawn on top of a plot canvas."""
    return (
        f"QCheckBox{{color:{THEME.TEXT_DIM};background:{THEME.PANEL};"
        f"padding:2px 6px;border:1px solid {THEME.BORDER};border-radius:6px;}}"
        f"QCheckBox:hover{{color:{THEME.TEXT};"
        f"border:1px solid {THEME.ACCENT};}}"
        f"QCheckBox:checked{{color:{THEME.TEXT};}}"
        f"QCheckBox::indicator{{width:12px;height:12px;"
        f"border:1px solid {THEME.BORDER};border-radius:3px;"
        f"background:{THEME.BG};}}"
        f"QCheckBox::indicator:hover{{border:1px solid {THEME.ACCENT};}}"
        f"QCheckBox::indicator:checked{{background:{THEME.ACCENT};"
        f"border:1px solid {THEME.ACCENT};}}"
    )


def _overlay_button_qss() -> str:
    """QSS for a compact pushbutton drawn on top of a plot canvas."""
    return (
        f"QPushButton{{color:{THEME.TEXT_DIM};background:{THEME.PANEL};"
        f"padding:2px 8px;border:1px solid {THEME.BORDER};border-radius:6px;"
        f"font:bold 9pt Monospace;}}"
        f"QPushButton:hover{{color:{THEME.TEXT};"
        f"border:1px solid {THEME.ACCENT};}}"
        f"QPushButton:disabled{{color:{THEME.TEXT_MUTED};}}"
    )


class _PlotCanvas(QWidget):
    """Canvas with a framed plot rect inset by the subclass's PAD_L/R/T/B
    margins and a title above it."""

    def _plot_rect(self) -> QRectF:
        w, h = self.width(), self.height()
        return QRectF(self.PAD_L, self.PAD_T,
                      max(1.0, w - self.PAD_L - self.PAD_R),
                      max(1.0, h - self.PAD_T - self.PAD_B))

    def _paint_frame(self, p: QPainter, title: str, title_dy: int) -> QRectF:
        """Fill the background, outline the plot rect and draw ``title``
        ``title_dy`` px above it; returns the plot rect."""
        p.fillRect(self.rect(), QColor(THEME.BG))
        r = self._plot_rect()
        p.setPen(QColor(THEME.BORDER))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(r)
        if title:
            f = QFont("Monospace", 10); f.setBold(True)
            p.setFont(f)
            p.setPen(QColor(THEME.TEXT))
            p.drawText(int(r.left()), int(r.top() - title_dy), title)
        return r


# ---- Hist1DWidget — QPainter bar chart with optional log Y ----

class Hist1DWidget(_PlotCanvas):
    PAD_L, PAD_R, PAD_T, PAD_B = 58, 14, 20, 20

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumHeight(140)
        self._bins: np.ndarray = np.zeros(0, dtype=np.int64)
        self._bmin: float = 0.0
        self._bstep: float = 1.0
        self._under: int = 0
        self._over: int = 0
        self._title: str = ""
        self._xlabel: str = ""
        self._color: QColor = QColor(THEME.ACCENT)
        self._log_y: bool = False
        self._hover_idx: int = -1

        self._logy_cb = QCheckBox("log Y", self)
        self._logy_cb.setFont(QFont("Monospace", 9, QFont.Weight.Bold))
        self._logy_cb.setStyleSheet(_overlay_checkbox_qss())
        self._logy_cb.toggled.connect(self.set_log_y)
        self._logy_cb.adjustSize()
        self._logy_cb.raise_()

    def set_data(self, bins, bmin: float, bstep: float,
                 under: int = 0, over: int = 0,
                 title: str = "", xlabel: str = "",
                 color: Optional[str] = None):
        self._bins = np.asarray(bins, dtype=np.int64)
        self._bmin = float(bmin)
        self._bstep = float(bstep)
        self._under = int(under)
        self._over  = int(over)
        self._title = title
        self._xlabel = xlabel
        if color:
            self._color = QColor(color)
        self._hover_idx = -1
        self.update()

    def set_log_y(self, on: bool):
        if on != self._log_y:
            self._log_y = on
            self.update()

    def clear(self, title: str = ""):
        self._bins = np.zeros(0, dtype=np.int64)
        self._title = title
        self._under = self._over = 0
        self._hover_idx = -1
        self.update()

    def resizeEvent(self, ev):
        cb = self._logy_cb
        cb.adjustSize()
        cb.move(self.width() - cb.width() - 6, 4)
        super().resizeEvent(ev)

    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        r = self._paint_frame(p, self._title, 8)

        n = self._bins.size
        if n == 0 or self._bins.sum() == 0:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, "(no data)")
            return

        if self._log_y:
            vals = np.where(self._bins > 0,
                            np.log10(self._bins.astype(np.float64)), 0.0)
            ymin, ymax = 0.0, float(vals.max())
        else:
            vals = self._bins.astype(np.float64)
            ymin, ymax = 0.0, float(vals.max())
        if ymax <= 0:
            ymax = 1.0

        bar_w = r.width() / n
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(self._color)
        for i in range(n):
            v = vals[i]
            if v <= 0:
                continue
            h = (v - ymin) / (ymax - ymin) * r.height()
            if h < 1.0:
                continue
            x0 = r.left() + i * bar_w
            y0 = r.bottom() - h
            p.fillRect(QRectF(x0, y0, max(bar_w, 1.0), h), self._color)

        if 0 <= self._hover_idx < n:
            x0 = r.left() + self._hover_idx * bar_w
            p.setPen(QPen(QColor(THEME.SELECT_BORDER), 1.2))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(QRectF(x0, r.top(), max(bar_w, 1.0), r.height()))

        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 8))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = r.bottom() - frac * r.height()
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            if self._log_y:
                val = 10 ** (ymin + frac * (ymax - ymin))
            else:
                val = ymin + frac * (ymax - ymin)
            p.drawText(int(r.left() - self.PAD_L + 2), int(y + 4),
                       _fmt_count(val))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            x = r.left() + frac * r.width()
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            val = self._bmin + frac * n * self._bstep
            p.drawText(int(x - 24), int(r.bottom() + 14), f"{val:g}")

        entries = int(self._bins.sum())
        info = f"N={entries:,}"
        if self._under:
            info += f"  under={self._under:,}"
        if self._over:
            info += f"  over={self._over:,}"
        p.setFont(QFont("Monospace", 9))
        p.setPen(QColor(THEME.TEXT_DIM))
        info_rect = QRectF(r.left(), r.top() - 20,
                           max(1.0, self.width() - self._logy_cb.width()
                               - 20 - r.left()),
                           14)
        p.drawText(info_rect,
                   Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                   info)

    def mouseMoveEvent(self, ev):
        r = self._plot_rect()
        n = self._bins.size
        if n == 0 or not r.contains(ev.position()):
            if self._hover_idx != -1:
                self._hover_idx = -1
                QToolTip.hideText()
                self.update()
            return
        idx = int((ev.position().x() - r.left()) / r.width() * n)
        idx = max(0, min(n - 1, idx))
        if idx != self._hover_idx:
            self._hover_idx = idx
            self.update()
        lo = self._bmin + idx * self._bstep
        hi = lo + self._bstep
        QToolTip.showText(ev.globalPosition().toPoint(),
                          f"[{lo:g}, {hi:g})  count={int(self._bins[idx]):,}",
                          self)

    def leaveEvent(self, _ev):
        if self._hover_idx != -1:
            self._hover_idx = -1
            QToolTip.hideText()
            self.update()


def _fmt_count(v: float) -> str:
    av = abs(v)
    if av == 0:
        return "0"
    if av >= 1e6:
        return f"{v/1e6:.2f}M"
    if av >= 1e3:
        return f"{v/1e3:.1f}k"
    if av >= 10:
        return f"{v:.0f}"
    return f"{v:.2g}"


# ---- WaveformPlotWidget — draws the current event's raw FADC samples ----

class WaveformPlotWidget(_PlotCanvas):
    PAD_L, PAD_R, PAD_T, PAD_B = 52, 14, 22, 30
    MAX_STACK = 200

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumHeight(150)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self._samples: np.ndarray = np.zeros(0, dtype=np.float32)
        self._peaks: List[Peak] = []
        self._ped_mean: float = 0.0
        self._ped_rms: float = 0.0
        self._title: str = ""
        self._clk_mhz: float = 250.0
        # Peaks rejected by the active filter are drawn dimmer and the cut
        # regions get a faint overlay (as in resources/waveform.js).
        # ``_filter_show`` is the "show" toggle: when False the overlays are
        # hidden but rejected peaks are still dimmed.
        self._filter: Optional[WaveformFilter] = None
        self._filter_show: bool = True

        # --- stack mode state ---
        self._stack_enabled: bool = False
        self._stack_traces: List[np.ndarray] = []   # bounded by MAX_STACK
        self._stack_key: str = ""                   # channel key; a new key resets the stack

        # --- overlay controls (top-right) ---
        self._stack_cb = QCheckBox("Stack", self)
        self._stack_cb.setFont(QFont("Monospace", 9, QFont.Weight.Bold))
        self._stack_cb.setStyleSheet(_overlay_checkbox_qss())
        self._stack_cb.setToolTip(
            f"Overlay waveforms across events (up to {self.MAX_STACK}). "
            "Peaks and integral shading hidden in stack mode.")
        self._stack_cb.toggled.connect(self._on_stack_toggled)
        self._stack_cb.adjustSize()

        self._stack_clear_btn = QPushButton("Clear", self)
        self._stack_clear_btn.setFont(QFont("Monospace", 9, QFont.Weight.Bold))
        self._stack_clear_btn.setStyleSheet(_overlay_button_qss())
        self._stack_clear_btn.setToolTip("Drop all stacked waveforms")
        self._stack_clear_btn.clicked.connect(self.clear_stack)
        self._stack_clear_btn.setVisible(False)
        self._stack_clear_btn.adjustSize()

        self._stack_count_lbl = QLabel("", self)
        self._stack_count_lbl.setFont(QFont("Monospace", 9))
        self._stack_count_lbl.setStyleSheet(
            f"color:{THEME.TEXT_DIM};background:transparent;")
        self._stack_count_lbl.setVisible(False)
        self._stack_count_lbl.adjustSize()

    # ---- Public API ----

    def set_data(self, samples: np.ndarray, peaks: List[Peak],
                 ped_mean: float, ped_rms: float,
                 title: str, clk_mhz: float = 250.0,
                 stack_key: Optional[str] = None,
                 filter: Optional["WaveformFilter"] = None,
                 filter_show: bool = True):
        samples = np.asarray(samples, dtype=np.float32)
        self._samples = samples
        self._peaks = peaks
        self._ped_mean = ped_mean
        self._ped_rms  = ped_rms
        self._title = title
        self._clk_mhz = clk_mhz
        self._filter = filter
        self._filter_show = bool(filter_show)

        if self._stack_enabled:
            key = stack_key if stack_key is not None else title
            if key != self._stack_key:
                self._stack_traces = []
                self._stack_key = key
            if samples.size >= 2:
                self._stack_traces.append(samples.copy())
                if len(self._stack_traces) > self.MAX_STACK:
                    self._stack_traces = self._stack_traces[-self.MAX_STACK:]
            self._update_stack_counter()
        self.update()

    def clear(self, title: str = ""):
        self._samples = np.zeros(0, dtype=np.float32)
        self._peaks = []
        self._title = title
        self.update()

    def clear_stack(self):
        """Drop every accumulated trace but keep the current waveform."""
        self._stack_traces = []
        self._stack_key = ""
        self._update_stack_counter()
        self.update()

    def reset_stack_if_new_key(self, key: str):
        """Reset traces when the caller switches to a different channel.

        Lets _display_waveform report a module change even when the new
        module has no samples in the current event — otherwise the empty
        early-return path would leave the previous module's stacks behind.
        """
        if self._stack_enabled and key != self._stack_key:
            self._stack_traces = []
            self._stack_key = key
            self._update_stack_counter()
            self.update()

    def is_stacking(self) -> bool:
        return self._stack_enabled

    # ---- Internals ----

    def _on_stack_toggled(self, on: bool):
        self._stack_enabled = on
        self._stack_clear_btn.setVisible(on)
        self._stack_count_lbl.setVisible(on)
        if not on:
            self._stack_traces = []
            self._stack_key = ""
        self._update_stack_counter()
        self._layout_overlays()
        self.update()

    def _update_stack_counter(self):
        self._stack_count_lbl.setText(
            f"{len(self._stack_traces)}/{self.MAX_STACK}")
        self._stack_count_lbl.adjustSize()
        self._layout_overlays()

    def _layout_overlays(self):
        # top-right: [count]  [Clear]  [Stack]
        margin = 6
        x = self.width() - margin
        y = 4
        x -= self._stack_cb.width()
        self._stack_cb.move(x, y)
        if self._stack_clear_btn.isVisible():
            x -= self._stack_clear_btn.width() + 4
            self._stack_clear_btn.move(x, y)
        if self._stack_count_lbl.isVisible():
            x -= self._stack_count_lbl.width() + 6
            self._stack_count_lbl.move(x, y + 2)

    def resizeEvent(self, ev):
        self._stack_cb.adjustSize()
        self._stack_clear_btn.adjustSize()
        self._stack_count_lbl.adjustSize()
        self._layout_overlays()
        super().resizeEvent(ev)

    # ---- Painting ----

    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        title = self._title
        if title and self._stack_enabled:
            title += f" — Stacked ({len(self._stack_traces)})"
        r = self._paint_frame(p, title, 6)

        if self._stack_enabled:
            self._paint_stacked(p, r)
        else:
            self._paint_single(p, r)

    # --- single-event (default) view ---------------------------------

    def _peak_passes(self, pk) -> bool:
        """Delegate to the active WaveformFilter.passes() so the dimming
        logic stays in lockstep with the histogram / geo filtering."""
        if self._filter is None:
            return True
        return self._filter.passes(pk)

    def _paint_single(self, p: QPainter, r: QRectF):
        n = self._samples.size
        if n < 2:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter,
                       "(no waveform — click Next to load an event)")
            return

        ymin = float(self._samples.min())
        ymax = float(self._samples.max())
        if ymax - ymin < 5.0:
            ymax = ymin + 5.0
        pad_y = (ymax - ymin) * 0.05
        ymin -= pad_y; ymax += pad_y

        def to_sx(i: float) -> float:
            return r.left() + (i / (n - 1)) * r.width()

        def to_sy(v: float) -> float:
            return r.bottom() - (v - ymin) / (ymax - ymin) * r.height()

        # ---- cut-region overlays (drawn first, low layer) ---------------
        # Mirrors xRangeShapes / yRangeShapes in resources/waveform.js.
        ns_total = (n - 1) * 1000.0 / self._clk_mhz if self._clk_mhz > 0 else 0.0
        cut_fill = QColor(THEME.TEXT_MUTED); cut_fill.setAlphaF(0.18)
        cut_edge = QColor(THEME.HIGHLIGHT)

        def shade(band: QRectF, x1: float, y1: float, x2: float, y2: float):
            """Shade a cut band and draw its dashed edge (x1,y1)-(x2,y2)."""
            p.setPen(Qt.PenStyle.NoPen)
            p.fillRect(band, cut_fill)
            p.setPen(QPen(cut_edge, 1, Qt.PenStyle.DashLine))
            p.drawLine(int(x1), int(y1), int(x2), int(y2))

        f = self._filter
        if f is not None and self._filter_show:
            # Time cut (left/right shaded bands outside [time_min, time_max]).
            if ns_total > 0 and f.time_min is not None and f.time_min > 0:
                x1 = to_sx((f.time_min * self._clk_mhz / 1000.0))
                x1 = min(max(x1, r.left()), r.right())
                if x1 > r.left():
                    shade(QRectF(r.left(), r.top(), x1 - r.left(), r.height()),
                          x1, r.top(), x1, r.bottom())
            if ns_total > 0 and f.time_max is not None and f.time_max < ns_total:
                x2 = to_sx((f.time_max * self._clk_mhz / 1000.0))
                x2 = min(max(x2, r.left()), r.right())
                if x2 < r.right():
                    shade(QRectF(x2, r.top(), r.right() - x2, r.height()),
                          x2, r.top(), x2, r.bottom())

        # pedestal baseline
        y_ped = to_sy(self._ped_mean) if self._ped_mean != 0 else None
        if y_ped is not None:
            p.setPen(QPen(QColor(THEME.TEXT_DIM), 1, Qt.PenStyle.DashLine))
            p.drawLine(int(r.left()), int(y_ped), int(r.right()), int(y_ped))
            # threshold line (same formula as waveform.js: pm + max(5*pr, 3))
            thr_v = self._ped_mean + max(5.0 * self._ped_rms, 3.0)
            y_thr = to_sy(thr_v)
            p.setPen(QPen(QColor(THEME.TEXT_MUTED), 1, Qt.PenStyle.DotLine))
            p.drawLine(int(r.left()), int(y_thr), int(r.right()), int(y_thr))
            # Height-cut shading — top/bottom bands for samples above
            # height_max + ped or below height_min + ped (filter is on
            # sample-pedestal units; we paint in raw ADC).
            if f is not None and self._filter_show:
                if f.height_min is not None:
                    hcut_v = self._ped_mean + f.height_min
                    if hcut_v > thr_v:
                        y_hcut = to_sy(hcut_v)
                        y_hcut = min(max(y_hcut, r.top()), r.bottom())
                        if y_hcut > r.top():
                            shade(QRectF(r.left(), y_hcut,
                                         r.width(), r.bottom() - y_hcut),
                                  r.left(), y_hcut, r.right(), y_hcut)
                if f.height_max is not None:
                    hcut_v = self._ped_mean + f.height_max
                    y_hcut = to_sy(hcut_v)
                    y_hcut = min(max(y_hcut, r.top()), r.bottom())
                    if y_hcut < r.bottom():
                        shade(QRectF(r.left(), r.top(),
                                     r.width(), y_hcut - r.top()),
                              r.left(), y_hcut, r.right(), y_hcut)

        # Fill the integral area (between pedestal and waveform) per peak,
        # colour-coded with series_qcolor. Mirrors resources/waveform.js.
        # Peaks rejected by the active filter are drawn with reduced alpha
        # so the user can tell at a glance which peaks the geo / hists use.
        if self._peaks and y_ped is not None:
            for i, pk in enumerate(self._peaks):
                passes = self._peak_passes(pk)
                base = series_qcolor(i)
                fill = QColor(base)
                fill.setAlphaF(0.18 if passes else 0.06)
                poly = QPolygonF()
                j = max(0, int(pk.left))
                j_end = min(n - 1, int(pk.right))
                for k in range(j, j_end + 1):
                    poly.append(QPointF(to_sx(k),
                                        to_sy(float(self._samples[k]))))
                # close along the pedestal baseline
                poly.append(QPointF(to_sx(j_end), y_ped))
                poly.append(QPointF(to_sx(j), y_ped))
                p.setPen(Qt.PenStyle.NoPen)
                p.setBrush(fill)
                p.drawPolygon(poly)
                # outline the peak section with the solid palette colour
                outline = QColor(base)
                if not passes:
                    outline.setAlphaF(0.45)
                p.setPen(QPen(outline, 2))
                p.setBrush(Qt.BrushStyle.NoBrush)
                for k in range(j, j_end):
                    p.drawLine(QPointF(to_sx(k),
                                       to_sy(float(self._samples[k]))),
                               QPointF(to_sx(k + 1),
                                       to_sy(float(self._samples[k + 1]))))

        # waveform line (default accent, drawn under peak outlines)
        p.setPen(QPen(QColor(THEME.ACCENT), 1.4))
        for i in range(n - 1):
            p.drawLine(int(to_sx(i)),     int(to_sy(float(self._samples[i]))),
                       int(to_sx(i + 1)), int(to_sy(float(self._samples[i + 1]))))

        # peak markers (diamonds, coloured per peak; rejected peaks get a
        # hollow diamond to keep them readable but visually distinct).
        if self._peaks:
            for i, pk in enumerate(self._peaks):
                if pk.pos < 0 or pk.pos >= n:
                    continue
                passes = self._peak_passes(pk)
                col = series_qcolor(i)
                if not passes:
                    col.setAlphaF(0.55)
                p.setPen(QPen(col, 1.2))
                p.setBrush(col if passes else Qt.BrushStyle.NoBrush)
                cx = to_sx(pk.pos)
                cy = to_sy(float(self._samples[pk.pos]))
                diamond = QPolygonF([
                    QPointF(cx,     cy - 4),
                    QPointF(cx + 4, cy),
                    QPointF(cx,     cy + 4),
                    QPointF(cx - 4, cy),
                ])
                p.drawPolygon(diamond)

        draw_wave_axes(p, r, ymin, ymax, n, self._clk_mhz, self.PAD_L)

        # ped/rms/peak-count readout — drawn inside the plot at top-right to
        # stay clear of the Stack checkbox / Clear button in the widget's
        # top-right margin.
        info = (f"ped={self._ped_mean:.1f}  rms={self._ped_rms:.2f}  "
                f"peaks={len(self._peaks)}")
        p.setFont(QFont("Monospace", 9))
        fm = p.fontMetrics()
        tw = fm.horizontalAdvance(info)
        th = fm.height()
        pad = 4
        box = QRectF(r.right() - tw - 2 * pad - 2, r.top() + 4,
                     tw + 2 * pad, th + 2)
        bg = QColor(THEME.BG); bg.setAlphaF(0.70)
        p.fillRect(box, bg)
        p.setPen(QColor(THEME.TEXT_DIM))
        p.drawText(box,
                   Qt.AlignmentFlag.AlignCenter, info)

    # --- stacked overlay view -----------------------------------------

    def _paint_stacked(self, p: QPainter, r: QRectF):
        traces = self._stack_traces
        if not traces:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter,
                       "(stack is empty — step through events to accumulate)")
            return

        ymin = min(float(w.min()) for w in traces)
        ymax = max(float(w.max()) for w in traces)
        if ymax - ymin < 5.0:
            ymax = ymin + 5.0
        pad_y = (ymax - ymin) * 0.05
        ymin -= pad_y; ymax += pad_y

        # Width uses the max length so shorter traces still fit left-aligned.
        n_max = max(w.size for w in traces)

        def to_sx(i: float, n: int) -> float:
            return r.left() + (i / max(1, n - 1)) * r.width()

        def to_sy(v: float) -> float:
            return r.bottom() - (v - ymin) / (ymax - ymin) * r.height()

        # Dimmed stacked traces.
        dim = QColor(THEME.ACCENT); dim.setAlphaF(0.18)
        p.setPen(QPen(dim, 1))
        for w in traces[:-1]:
            n = w.size
            for i in range(n - 1):
                p.drawLine(int(to_sx(i, n)),     int(to_sy(float(w[i]))),
                           int(to_sx(i + 1, n)), int(to_sy(float(w[i + 1]))))

        # Latest trace drawn on top at full colour.
        latest = traces[-1]
        n = latest.size
        p.setPen(QPen(QColor(THEME.ACCENT), 1.4))
        for i in range(n - 1):
            p.drawLine(int(to_sx(i, n)),     int(to_sy(float(latest[i]))),
                       int(to_sx(i + 1, n)), int(to_sy(float(latest[i + 1]))))

        draw_wave_axes(p, r, ymin, ymax, n_max, self._clk_mhz, self.PAD_L)

        p.setPen(QColor(THEME.TEXT_DIM))
        p.drawText(QRectF(r.left(), r.top() - 20,
                          max(1.0, r.width() - 8), 14),
                   Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                   f"stack={len(traces)}/{self.MAX_STACK}")


# ---- WaveformGeoView — small HyCal overview for module selection ----

class WaveformGeoView(HyCalMapWidget):
    """Compact HyCal geo view with two colour-coding modes.

    * ``current``  — module colour = max peak integral in the current event.
    * ``overall``  — module colour = occupancy (events-with-peak / accumulated
      events) across all events the user has browsed / batched.

    Modules that have never been seen in an event are drawn in a flat grey
    so "no-data" stays visually distinct from "data, low value".  Clicking
    any module emits moduleClicked with its name.
    """

    MODE_CURRENT = "current"
    MODE_OVERALL = "overall"
    SELECT_PEN_WIDTH = 2.0

    # Resolved at paint time so the active theme wins; see :class:`THEME`.
    @property
    def UNAVAIL_COLOR(self) -> QColor:
        return QColor(THEME.BORDER)

    def __init__(self, parent=None):
        # margin_bottom must exceed the base's colour-bar anchor (cb_y =
        # h - 40) so the module rects clear the bar — leave ~16 px gap.
        # The tiny LMS / V blocks off to the left of HyCal get name labels.
        super().__init__(parent, show_colorbar=True, include_lms=True,
                         label_types=AUX_TYPES,
                         margin_top=4, margin_bottom=56,
                         min_size=(220, 280), shrink=0.90)
        self._available: set = set()
        self._mode = self.MODE_CURRENT
        self._current_vals: Dict[str, float] = {}
        self._overall_vals: Dict[str, float] = {}
        # Modules whose peaks in the current event were all rejected by the
        # peak filter: shaded on the waveform plot but absent from the geo's
        # max-integral.  Used by _tooltip_text.
        self._rejected_current: set = set()
        # Headless range controller: handles auto-fit logic.  Both vmin and
        # vmax are inline-editable on the colorbar so the user can pin the
        # palette to a fixed range when comparing events.
        self._range_ctrl = ColorRangeController(
            self, auto_fit="minmax", parent=self)

        # Top-left mode toggle.  Default label matches MODE_CURRENT.
        self._mode_btn = QPushButton("Current", self)
        self._mode_btn.setFixedSize(74, 22)
        _f = QFont("Consolas", 9); _f.setBold(True)
        self._mode_btn.setFont(_f)
        self._mode_btn.setToolTip(
            "Colour coding:\n"
            "  Current — max peak integral in the currently viewed event\n"
            "  Overall — occupancy (events-with-peak / accumulated events)")
        self._mode_btn.setStyleSheet(themed(OVERLAY_BUTTON_QSS))
        self._mode_btn.clicked.connect(self._toggle_mode)

    def set_available(self, names):
        self._available = set(names)
        self.update()

    def set_current_values(self, vals: Dict[str, float]):
        self._current_vals = vals
        if self._mode == self.MODE_CURRENT:
            self._apply_mode_values()

    def set_rejected_current(self, names: set):
        """Modules where the current event had peaks rejected by the cut.
        Pass an empty set to clear."""
        self._rejected_current = set(names) if names else set()

    def set_overall_values(self, vals: Dict[str, float]):
        self._overall_vals = vals
        if self._mode == self.MODE_OVERALL:
            self._apply_mode_values()

    def _toggle_mode(self):
        self._mode = (self.MODE_OVERALL if self._mode == self.MODE_CURRENT
                      else self.MODE_CURRENT)
        self._mode_btn.setText(
            "Overall" if self._mode == self.MODE_OVERALL else "Current")
        self._apply_mode_values()

    def _apply_mode_values(self):
        if self._mode == self.MODE_OVERALL:
            self.set_values(self._overall_vals)
            self.set_range(0.0, 1.0)          # occupancy fraction
        else:
            self.set_values(self._current_vals)
            # Re-fit per event; user can override via inline colorbar edit.
            self._range_ctrl.auto_fit(self._current_vals)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        self._mode_btn.move(6, 6)

    def _paint_modules(self, p):
        avail = self._available
        u_col = self.UNAVAIL_COLOR
        no_data = self.NO_DATA_COLOR
        stops = self.palette_stops()
        vals = self._values
        for name, rect in self._rects.items():
            if name not in avail:
                p.fillRect(rect, u_col)
                continue
            v = vals.get(name)
            if v is None:
                p.fillRect(rect, no_data)
                continue
            p.fillRect(rect, cmap_qcolor(self.value_to_t(v), stops))

    def _tooltip_text(self, name: str) -> str:
        if name not in self._available:
            return f"{name}  (not seen yet)"
        v = self._values.get(name)
        unit = ("occupancy" if self._mode == self.MODE_OVERALL
                else "max integral")
        if v is None:
            if (self._mode == self.MODE_CURRENT
                    and name in self._rejected_current):
                return f"{name}  ({unit}: — peak outside cut)"
            return f"{name}  ({unit}: —)"
        return f"{name}  {unit}={v:.3g}"


# ---- Cluster display wrappers ----


class _DisplayCluster:
    """Display-friendly wrapper around a bound ClusterHit.  We can't set
    arbitrary attributes on pybind11 objects, so copy the displayed hit
    fields and add the centre name + member-name list."""
    __slots__ = ("energy", "x", "y", "nblocks", "center_name", "members")

    def __init__(self, hit, center_name: str, members: List[str]):
        self.energy      = float(hit.energy)
        self.x           = float(hit.x)
        self.y           = float(hit.y)
        self.nblocks     = int(hit.nblocks)
        self.center_name = center_name or ""
        self.members     = list(members)


# ---- Cluster map widget — HyCal heatmap + cluster overlays ----


class HyCalClusterMap(HyCalMapWidget):
    """HyCal map coloured by per-module energy (MeV) with cluster overlays
    (crosshair + small circle + energy label) mirroring the web monitor."""

    def __init__(self, parent=None):
        super().__init__(parent, include_lms=True, enable_zoom_pan=True,
                         show_colorbar=True)
        self._clusters: List = []               # List[_DisplayCluster]
        self._selected_cluster: Optional[int] = None
        self._member_modules: set = set()       # module names of selected cluster

    # -- public API ------------------------------------------------------

    def set_clusters(self, clusters):
        self._clusters = list(clusters) if clusters else []
        self._recompute_membership()
        self.update()

    def set_selected_cluster(self, idx: Optional[int]):
        """Highlight one cluster; pass None for 'show all'."""
        if idx is not None and not (0 <= idx < len(self._clusters)):
            idx = None
        self._selected_cluster = idx
        self._recompute_membership()
        self.update()

    # -- internals -------------------------------------------------------

    def _recompute_membership(self):
        """Build set of module names that belong to the selected cluster:
        its ``members`` list, or just ``center_name`` when that is empty."""
        self._member_modules.clear()
        if self._selected_cluster is None:
            return
        cl = self._clusters[self._selected_cluster]
        members = getattr(cl, "members", None)
        if members:
            self._member_modules.update(members)
        elif getattr(cl, "center_name", ""):
            self._member_modules.add(cl.center_name)

    def _paint_modules(self, p):
        """Default colormap paint, but dim non-members (low alpha) when a
        cluster is selected."""
        dim = self._selected_cluster is not None and self._member_modules
        if not dim:
            super()._paint_modules(p)
            return
        stops = self.palette_stops()
        no_data = self.NO_DATA_COLOR
        for name, rect in self._rects.items():
            v = self._values.get(name)
            if v is None:
                col = QColor(no_data)
            else:
                col = cmap_qcolor(self.value_to_t(v), stops)
            if name not in self._member_modules:
                col = QColor(col.red(), col.green(), col.blue(), 60)
            p.fillRect(rect, col)

    def _paint_cluster_frames(self, p):
        """Draw a coloured border around every member module of every
        cluster (or just the selected cluster).  Mirrors the per-module
        border colouring in resources/cluster.js."""
        if not self._clusters:
            return
        sel = self._selected_cluster
        p.save()
        p.setBrush(Qt.BrushStyle.NoBrush)
        for i, cl in enumerate(self._clusters):
            if sel is not None and i != sel:
                continue
            members = getattr(cl, "members", None) or ()
            if not members:
                continue
            width = 2.5 if (sel is not None and i == sel) else 1.5
            p.setPen(QPen(series_qcolor(i), width))
            for name in members:
                rect = self._rects.get(name)
                if rect is not None:
                    p.drawRect(rect)
        p.restore()

    def _paint_overlays(self, p, w, h):
        # Per-cluster coloured frames first, so the hover border (drawn by
        # super) and the cluster crosshairs/labels stay visible on top.
        self._paint_cluster_frames(p)
        super()._paint_overlays(p, w, h)   # hover border
        if not self._clusters:
            return

        p.save()
        cross_pen = QPen(QColor("#ffd166"), 1.6)
        circle_pen = QPen(QColor("#ffd166"), 1.4)
        font = QFont("Monospace", 9, QFont.Weight.Bold)
        p.setFont(font)
        fm = p.fontMetrics()

        for i, cl in enumerate(self._clusters):
            if self._selected_cluster is not None and i != self._selected_cluster:
                continue
            pt = self.geo_to_canvas(cl.x, cl.y)
            # crosshair
            p.setPen(cross_pen)
            L = 9
            p.drawLine(QPointF(pt.x() - L, pt.y()),
                       QPointF(pt.x() + L, pt.y()))
            p.drawLine(QPointF(pt.x(), pt.y() - L),
                       QPointF(pt.x(), pt.y() + L))
            # circle scaled by log-ish energy
            r = max(6.0, min(20.0, 4.0 + 2.0 * (cl.energy ** 0.33)))
            p.setPen(circle_pen)
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawEllipse(pt, r, r)
            text = f"C{i}  {cl.energy:.0f}"
            tw = fm.horizontalAdvance(text)
            tx = pt.x() + r + 3
            ty = pt.y() - 3
            # background halo for readability
            halo = QColor(0, 0, 0, 140)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(halo)
            p.drawRect(QRectF(tx - 2, ty - fm.ascent(),
                              tw + 4, fm.height()))
            p.setPen(QColor("#ffd166"))
            p.drawText(QPointF(tx, ty), text)
        p.restore()


# ---- Cluster panel — selector + table + footer stats ----


class ClusterPanel(QWidget):
    """Right-side pane for the Cluster tab: combo + table + summary line."""

    clusterSelected = pyqtSignal(object)    # int | None  (None = show all)

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(6, 6, 6, 6)
        lay.setSpacing(4)

        self._selector = QComboBox()
        self._selector.addItem("All clusters", None)
        self._selector.currentIndexChanged.connect(self._on_selector_changed)
        lay.addWidget(self._selector)

        self._table = QTableWidget(0, 6)
        self._table.setHorizontalHeaderLabels(
            ["#", "center", "E [MeV]", "x [mm]", "y [mm]", "nblocks"])
        self._table.verticalHeader().setVisible(False)
        self._table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection)
        self._table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        hh = self._table.horizontalHeader()
        hh.setSectionResizeMode(QHeaderView.ResizeMode.ResizeToContents)
        hh.setStretchLastSection(True)
        self._table.itemSelectionChanged.connect(self._on_table_selection)
        lay.addWidget(self._table, stretch=1)

        self._summary = QLabel("no event")
        self._summary.setFont(QFont("Monospace", 9))
        lay.addWidget(self._summary)

        self._clusters: List = []

    def set_clusters(self, clusters):
        self._clusters = list(clusters) if clusters else []
        # Repopulate selector
        self._selector.blockSignals(True)
        self._selector.clear()
        self._selector.addItem("All clusters", None)
        for i, cl in enumerate(self._clusters):
            label = f"C{i}  E={cl.energy:.0f} MeV  {getattr(cl, 'center_name', '?')}"
            self._selector.addItem(label, i)
        self._selector.setCurrentIndex(0)
        self._selector.blockSignals(False)

        # Populate table
        self._table.blockSignals(True)
        self._table.setRowCount(len(self._clusters))
        chip_font = QFont("Monospace", 9, QFont.Weight.Bold)
        text_black = QBrush(QColor("#000000"))
        for i, cl in enumerate(self._clusters):
            row = [
                f"{i}",
                getattr(cl, "center_name", "?"),
                f"{cl.energy:.1f}",
                f"{cl.x:.1f}",
                f"{cl.y:.1f}",
                f"{cl.nblocks}",
            ]
            for c, v in enumerate(row):
                item = QTableWidgetItem(v)
                if c == 0:
                    # Colour chip linking the row to its cluster colour.
                    item.setBackground(QBrush(series_qcolor(i)))
                    item.setForeground(text_black)
                    item.setFont(chip_font)
                    item.setTextAlignment(Qt.AlignmentFlag.AlignCenter)
                elif c in (2, 3, 4, 5):
                    item.setTextAlignment(Qt.AlignmentFlag.AlignRight
                                          | Qt.AlignmentFlag.AlignVCenter)
                self._table.setItem(i, c, item)
        self._table.blockSignals(False)

        n = len(self._clusters)
        tot = sum(c.energy for c in self._clusters)
        self._summary.setText(f"{n} clusters   ΣE = {tot:.0f} MeV")

    def _on_selector_changed(self, _idx: int):
        data = self._selector.currentData()
        # Sync table selection without echoing
        self._table.blockSignals(True)
        if data is None:
            self._table.clearSelection()
        else:
            self._table.selectRow(int(data))
        self._table.blockSignals(False)
        self.clusterSelected.emit(data)

    def _on_table_selection(self):
        rows = self._table.selectionModel().selectedRows()
        if not rows:
            return
        i = rows[0].row()
        # Sync combo
        self._selector.blockSignals(True)
        self._selector.setCurrentIndex(i + 1)   # +1 for leading "All"
        self._selector.blockSignals(False)
        self.clusterSelected.emit(i)


# ---- Cut Settings dialog — modal editor for the active WaveformFilter ----


class _RangeRow(QWidget):
    """Two QLineEdits (min/max) with a validator.  Empty value = no
    constraint, matching the web monitor's <input type='number'> + null
    parsing in cut_dialog.js."""

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QHBoxLayout(self)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(8)
        self.min_edit = QLineEdit()
        self.max_edit = QLineEdit()
        for w, ph in ((self.min_edit, "min"), (self.max_edit, "max")):
            w.setPlaceholderText(ph)
            w.setValidator(QDoubleValidator(self))
            w.setMinimumWidth(80)
        lay.addWidget(QLabel("min"))
        lay.addWidget(self.min_edit, 1)
        lay.addSpacing(6)
        lay.addWidget(QLabel("max"))
        lay.addWidget(self.max_edit, 1)

    def set_values(self, lo: Optional[float], hi: Optional[float]):
        self.min_edit.setText("" if lo is None else _fmt_filter_num(lo))
        self.max_edit.setText("" if hi is None else _fmt_filter_num(hi))

    def values(self) -> Tuple[Optional[float], Optional[float]]:
        def _parse(s: str) -> Optional[float]:
            s = s.strip()
            if not s:
                return None
            try:
                return float(s)
            except ValueError:
                return None
        return _parse(self.min_edit.text()), _parse(self.max_edit.text())


def _fmt_filter_num(v: float) -> str:
    """Compact representation that round-trips through float() — mirrors
    JSON.stringify(num) on the web side."""
    if v == int(v):
        return f"{int(v)}"
    return f"{v:g}"


class CutSettingsDialog(QDialog):
    """Modal "Cut Settings" editor — replicates resources/cut_dialog.js
    (and viewer.html's #cut-dialog markup) field for field.

    The user edits the active filter's per-axis ranges and quality bits,
    then clicks Save to commit.  Cancel discards changes.  Reset reverts
    the form (without committing) to the defaults snapshotted at startup
    from monitor_config.json's ``waveform.filter`` block.
    """

    def __init__(self, parent: QWidget,
                 current: WaveformFilter,
                 default: WaveformFilter):
        super().__init__(parent)
        self.setWindowTitle("Cut Settings")
        self.setModal(True)
        self._default = default

        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)
        root.setSpacing(8)

        self._rows: Dict[str, _RangeRow] = {}
        for axis, label, suffix in (("time",     "Time",     " (ns)"),
                                    ("integral", "Integral", ""),
                                    ("height",   "Height",   "")):
            gb = QGroupBox(label + suffix)
            gl = QVBoxLayout(gb)
            gl.setContentsMargins(8, 4, 8, 6)
            row = _RangeRow()
            self._rows[axis] = row
            gl.addWidget(row)
            root.addWidget(gb)

        # Quality bits — two columns of checkboxes (accept / reject).
        qg = QGroupBox("Quality bits")
        qg_lay = QGridLayout(qg)
        qg_lay.setContentsMargins(8, 4, 8, 6)
        qg_lay.setHorizontalSpacing(20)
        accept_lbl = QLabel("Accept")
        reject_lbl = QLabel("Reject")
        for lbl in (accept_lbl, reject_lbl):
            f = QFont(); f.setBold(True)
            lbl.setFont(f)
        qg_lay.addWidget(accept_lbl, 0, 0)
        qg_lay.addWidget(reject_lbl, 0, 1)

        self._accept_checks: Dict[str, QCheckBox] = {}
        self._reject_checks: Dict[str, QCheckBox] = {}
        if not PEAK_QUALITY_BITS:
            note = QLabel("(no quality bits exposed by prad2py)")
            note.setStyleSheet(f"color:{THEME.TEXT_DIM}; font-style:italic;")
            qg_lay.addWidget(note, 1, 0, 1, 2)
        else:
            for r, b in enumerate(PEAK_QUALITY_BITS, start=1):
                acc = QCheckBox(str(b["label"]))
                rej = QCheckBox(str(b["label"]))
                self._accept_checks[str(b["name"])] = acc
                self._reject_checks[str(b["name"])] = rej
                acc.toggled.connect(self._sync_bit_mutex)
                rej.toggled.connect(self._sync_bit_mutex)
                qg_lay.addWidget(acc, r, 0)
                qg_lay.addWidget(rej, r, 1)
        hint = QLabel(
            "Empty = no constraint.  Accept: peak's set bits must overlap "
            "the accepted flags.  Reject: peak fails if any rejected flag is set.")
        hint.setStyleSheet(f"color:{THEME.TEXT_DIM}; font-size:10px;")
        hint.setWordWrap(True)
        qg_lay.addWidget(hint, qg_lay.rowCount(), 0, 1, 2)
        root.addWidget(qg)

        bb = QDialogButtonBox(
            QDialogButtonBox.StandardButton.Save
            | QDialogButtonBox.StandardButton.Cancel
            | QDialogButtonBox.StandardButton.Reset)
        bb.button(QDialogButtonBox.StandardButton.Save).setDefault(True)
        bb.button(QDialogButtonBox.StandardButton.Reset).setToolTip(
            "Restore the form to monitor_config.json's filter values "
            "(does not commit until you click Save).")
        bb.accepted.connect(self.accept)
        bb.rejected.connect(self.reject)
        bb.button(QDialogButtonBox.StandardButton.Reset).clicked.connect(
            self._on_reset)
        root.addWidget(bb)

        self._populate(current)

    def _populate(self, f: WaveformFilter):
        self._rows["time"]    .set_values(f.time_min,     f.time_max)
        self._rows["integral"].set_values(f.integral_min, f.integral_max)
        self._rows["height"]  .set_values(f.height_min,   f.height_max)
        for d in PEAK_QUALITY_BITS:
            name = str(d["name"]); mask = int(d["mask"])
            acc = self._accept_checks.get(name)
            rej = self._reject_checks.get(name)
            if acc is not None: acc.setChecked(bool(f.q_accept & mask))
            if rej is not None: rej.setChecked(bool(f.q_reject & mask))
        self._sync_bit_mutex()

    def _on_reset(self):
        self._populate(self._default)

    def _sync_bit_mutex(self):
        """A bit can be in Accept OR Reject but not both — disable the
        twin checkbox when its sibling is checked."""
        for name, acc in self._accept_checks.items():
            rej = self._reject_checks.get(name)
            if rej is None:
                continue
            acc.setEnabled(not rej.isChecked())
            rej.setEnabled(not acc.isChecked())

    def result_filter(self) -> WaveformFilter:
        """Build a WaveformFilter from the current form contents.  ``enable``
        is left at its default; the caller restores it from the toolbar's
        "apply" toggle."""
        f = WaveformFilter()
        f.time_min,     f.time_max     = self._rows["time"].values()
        f.integral_min, f.integral_max = self._rows["integral"].values()
        f.height_min,   f.height_max   = self._rows["height"].values()
        f.q_accept = _names_to_mask([n for n, cb in self._accept_checks.items()
                                     if cb.isChecked()])
        f.q_reject = _names_to_mask([n for n, cb in self._reject_checks.items()
                                     if cb.isChecked()])
        return f


# ---- Main window ----

_NATKEY_RE = re.compile(r"(\d+)")
def _natural_sort_key(s: str):
    return [int(p) if p.isdigit() else p.lower()
            for p in _NATKEY_RE.split(s or "")]


# Advanced-dock rows: (config field, min, max, step, tooltip); step None
# makes an integer field.
_WAVE_FIELDS = (
    ("peak_nsigma",      0.0,    50.0, 0.5,  "peak detection threshold (× pedestal RMS)"),
    ("min_peak_height",  0.0,  4096.0, 1.0,  "absolute floor on detected peak height (ADC)"),
    ("min_peak_ratio",   0.0,     1.0, 0.01, "secondary/primary peak ratio"),
    ("int_tail_ratio",   0.0,     1.0, 0.01, "tail integration cut"),
    ("ped_flatness",     0.0,  1000.0, 0.5,  "pedestal RMS ceiling"),
    ("clk_mhz",          1.0,  1000.0, 1.0,  "FADC clock (MHz)"),
    ("smooth_order",     1,      16,   None, "kernel order (1 = identity, N gives 2N-1 taps)"),
    ("ped_nsamples",     1,      64,   None, "samples to use for pedestal"),
    ("ped_max_iter",     1,     100,   None, "pedestal iteration cap"),
    ("overflow",         0,   65535,   None, "overflow cutoff (ADC)"),
)

_CLUSTER_FIELDS = (
    ("min_module_energy",  0.0, 1e4, 0.1, "single-module threshold (MeV)"),
    ("min_center_energy",  0.0, 1e4, 0.1, "seed threshold (MeV)"),
    ("min_cluster_energy", 0.0, 1e5, 0.1, "total cluster threshold (MeV)"),
    ("log_weight_thres",   0.0, 20.0, 0.1, "log-weight offset"),
    ("least_split",        0.0, 1.0, 0.01, "min fraction to keep a split hit"),
    ("seed_time_window",   -1.0, 200.0, 0.5,
        "Multi-pulse seed-time gate (ns).  ≤0 disables timing "
        "gating (legacy single-pulse-per-module mode).  >0 lets "
        "AddHit() be called once per pulse; FormClusters then "
        "groups neighbours within ±this window of the seed pulse.\n"
        "Persistent default: 'seed_time_window' under the 'hycal' "
        "block in database/reconstruction_config.json."),
    ("min_cluster_size", 1, 100, None, "min modules in cluster"),
    ("split_iter",       0, 100, None, "island-split iteration cap"),
)


class HyCalEventViewer(QMainWindow):

    def __init__(self,
                 *,
                 hist_config: Dict,
                 daq_map: Dict,
                 roc_to_crate: Dict,
                 accept_mask: int,
                 reject_mask: int,
                 daq_config_path: str,
                 hycal_modules: Optional[List] = None,
                 hycal_map_path: Optional[str] = None,
                 recon_config_path: Optional[str] = None):
        super().__init__()
        self._daq_map       = daq_map
        self._roc_to_crate  = roc_to_crate
        self._accept_mask   = accept_mask
        self._reject_mask   = reject_mask
        self._daq_cfg_path  = daq_config_path
        self._hycal_modules = hycal_modules or []
        self._hycal_map_path = hycal_map_path
        self._recon_cfg_path = recon_config_path

        # Bin configs by ChannelHists field — user config over defaults.
        self._hist_cfg: Dict[str, Dict] = {
            spec.attr: (hist_config.get(spec.cfg_key, spec.default)
                        if spec.cfg_key else spec.default)
            for spec in _HIST_SPECS}

        # Seed the analyzer config from daq_config.json's
        # `fadc250_waveform.analyzer` block.  This is the single source of
        # truth for peak-detection knobs (peak_nsigma, min_peak_height,
        # min_peak_ratio).  Falls back to plain defaults if the daq_config
        # can't be loaded — opening files later will fail loudly anyway.
        try:
            _dc = prad2py.dec.load_daq_config(self._daq_cfg_path or "")
            self._wcfg = WaveConfig(_dc.wave_cfg)
        except Exception:
            self._wcfg = WaveConfig()

        # Peak filter from monitor_config.json's `waveform.filter`; the
        # `_filter_default` snapshot backs the Cut-Settings "Reset" button.
        # `_filter_show` is the overlay toggle on the waveform plot.
        _flt_json = hist_config.get("filter") or {}
        self._filter         = WaveformFilter.from_json(_flt_json)
        self._filter_default = WaveformFilter.from_json(_flt_json)
        self._filter_show    = True

        # Debounce timer for Advanced-dock re-runs: coalesce rapid
        # slider drags into a single re-read + re-analyse.  Must exist
        # before the dock widgets are built (they connect to it via
        # _on_advanced_changed).
        self._adv_debounce_ms = 150
        self._adv_redraw_timer = QTimer(self)
        self._adv_redraw_timer.setSingleShot(True)
        self._adv_redraw_timer.timeout.connect(self._rerun_current_event)

        # File state
        self._evio_path: Optional[Path] = None
        self._index: List[Tuple[int, int]] = []
        self._current_idx: int = -1
        # Folded-event flags (see _is_folded); 1 byte/event, so 10 M
        # events ≈ 10 MB.
        self._accumulated: Optional[np.ndarray] = None

        # Per-channel accumulated hists, keyed by (roc, slot, ch)
        self._channels: Dict[Tuple[int, int, int], ChannelHists] = {}
        self._selected_key: Optional[Tuple[int, int, int]] = None

        # Browse handle, kept open across navigation
        self._reader: Optional[EvioCursor] = None

        # Worker threads
        self._idx_worker: Optional[IndexerWorker] = None
        self._idx_thread: Optional[QThread] = None
        self._batch_worker: Optional[BatchWorker] = None
        self._batch_thread: Optional[QThread] = None

        # HyCal clustering (see _build_hycal_pipeline).  `_pipeline` keeps
        # the C++ Pipeline object alive; `_hcsys` borrows its `hycal`.
        # DAQ → module lookups are cached by (roc_tag, slot, ch) so the
        # per-event hot loop doesn't go through pybind11 every hit.
        self._pipeline = None
        self._hcsys = None
        self._hccl  = None
        self._hc_cache: Dict[Tuple[int, int, int], object] = {}
        # True once a calibration file has actually been loaded (either via
        # PipelineBuilder's runinfo lookup or the manual menu override).
        # Gates the "no calibration" warning banner on the cluster tab.
        self._hycal_calib_loaded = False
        if not self._build_hycal_pipeline():
            self._hycal_init_fallback()

        apply_theme_palette(self)
        self._build_ui()
        self._make_menu()

    # -- UI --

    def _build_ui(self):
        self.setWindowTitle("HyCal Event Viewer")
        self.resize(1500, 1000)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(3)

        self._file_lbl = QLabel("(no file loaded)")
        self._file_lbl.setFont(QFont("Monospace", 10))
        self._file_lbl.setStyleSheet(themed("color:#8b949e;"))
        root.addWidget(self._file_lbl)

        # -- top control bar: navigation on the left, module picker on the right --
        top = QHBoxLayout()

        self._prev_btn = self._small_btn("◀ Prev", self._on_prev)
        self._next_btn = self._small_btn("Next ▶", self._on_next)
        self._prev_btn.setEnabled(False)
        self._next_btn.setEnabled(False)
        top.addWidget(self._prev_btn)
        top.addWidget(self._next_btn)

        top.addSpacing(12)
        top.addWidget(self._mk_label("Event:"))
        self._event_spin = QSpinBox()
        self._event_spin.setFont(QFont("Monospace", 10))
        self._event_spin.setMinimum(0)
        self._event_spin.setMaximum(0)
        self._event_spin.setStyleSheet(themed(
            "QSpinBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:6px;padding:2px 6px;}"))
        self._event_spin.editingFinished.connect(self._on_spin_jump)
        self._event_spin.setEnabled(False)
        top.addWidget(self._event_spin)
        self._total_lbl = QLabel(" / 0")
        self._total_lbl.setFont(QFont("Monospace", 10))
        self._total_lbl.setStyleSheet(themed("color:#8b949e;"))
        top.addWidget(self._total_lbl)

        top.addSpacing(18)
        self._batch_btn = self._small_btn("Process next 10k",
                                          self._on_batch_10k, primary=True)
        self._batch_btn.setEnabled(False)
        top.addWidget(self._batch_btn)

        self._batch_status = QLabel("")
        self._batch_status.setFont(QFont("Monospace", 10))
        self._batch_status.setStyleSheet(themed("color:#8b949e;"))
        top.addSpacing(8)
        top.addWidget(self._batch_status)

        top.addStretch(1)

        # Right cluster: module dropdown + reset hist.
        mod_lbl = QLabel("Module:")
        mod_lbl.setFont(QFont("Monospace", 11, QFont.Weight.Bold))
        mod_lbl.setStyleSheet(themed("color:#c9d1d9;"))
        top.addWidget(mod_lbl)
        self._combo = QComboBox()
        self._combo.setEditable(True)
        self._combo.setInsertPolicy(QComboBox.InsertPolicy.NoInsert)
        self._combo.setFont(QFont("Monospace", 11))
        self._combo.setMinimumContentsLength(32)
        self._combo.setStyleSheet(themed(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:6px;padding:2px 6px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
            "selection-background-color:#1f6feb;}"))
        comp = self._combo.completer()
        if comp is not None:
            comp.setFilterMode(Qt.MatchFlag.MatchContains)
            comp.setCaseSensitivity(Qt.CaseSensitivity.CaseInsensitive)
            comp.setCompletionMode(QCompleter.CompletionMode.PopupCompletion)
        self._combo.currentIndexChanged.connect(self._on_combo_changed)
        top.addWidget(self._combo)
        self._reset_btn = self._small_btn("Reset hist", self._reset_current_hists)
        self._reset_btn.setEnabled(False)
        top.addWidget(self._reset_btn)

        root.addLayout(top)

        self._info = QLabel("")
        self._info.setFont(QFont("Monospace", 10))
        self._info.setStyleSheet(themed("color:#8b949e;"))
        root.addWidget(self._info)

        # -- tabbed central area: Waveform + Cluster -----------------------
        self._tabs = QTabWidget()
        self._tabs.addTab(self._build_waveform_tab(), "Waveform")
        self._tabs.addTab(self._build_cluster_tab(),  "Cluster")
        root.addWidget(self._tabs, stretch=1)

        # -- advanced dock (hidden by default) -----------------------------
        self._adv_dock = self._build_advanced_dock()

        self.setStatusBar(QStatusBar())
        self._clear_plots()

        # Keyboard: ← / → to navigate prev / next.
        QShortcut(QKeySequence(Qt.Key.Key_Left),  self, activated=self._on_prev)
        QShortcut(QKeySequence(Qt.Key.Key_Right), self, activated=self._on_next)

    # ---- Tabs -----------------------------------------------------------

    def _build_waveform_tab(self) -> QWidget:
        """Geo + waveform on the left, four histograms stacked on the
        right."""
        tab = QWidget()
        lay = QVBoxLayout(tab)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(3)

        # Cut-Settings toolbar (the web monitor's .tcut-bar): dialog button
        # plus the "apply" / "show" toggles.
        cut_bar = QHBoxLayout()
        cut_bar.setContentsMargins(2, 0, 2, 0)
        cut_bar.setSpacing(8)
        self._cut_settings_btn = QPushButton("Cut Settings…")
        self._cut_settings_btn.setToolTip("Edit waveform peak filter ranges")
        self._cut_settings_btn.setStyleSheet(themed(
            "QPushButton{background:#21262d;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;"
            "padding:3px 10px;font:bold 9pt Monospace;}"
            "QPushButton:hover{background:#30363d;color:#e6edf3;}"))
        self._cut_settings_btn.clicked.connect(self._open_cut_settings_dialog)
        cut_bar.addWidget(self._cut_settings_btn)
        self._cut_apply_cb = QCheckBox("apply")
        self._cut_apply_cb.setChecked(self._filter.enable)
        self._cut_apply_cb.setToolTip(
            "Apply the peak filter to histograms / geo coloring / clustering. "
            "When off, every analyzed peak counts.")
        self._cut_apply_cb.toggled.connect(self._on_cut_apply_toggled)
        cut_bar.addWidget(self._cut_apply_cb)
        self._cut_show_cb = QCheckBox("show")
        self._cut_show_cb.setChecked(self._filter_show)
        self._cut_show_cb.setToolTip(
            "Show cut-range overlays on the waveform plot.  Independent of "
            "apply — overlays can be hidden while the filter is active.")
        self._cut_show_cb.toggled.connect(self._on_cut_show_toggled)
        cut_bar.addWidget(self._cut_show_cb)
        cut_bar.addStretch(1)
        lay.addLayout(cut_bar)

        split = QSplitter(Qt.Orientation.Horizontal)

        # Left: geo view (square, top) + waveform plot (bottom)
        left = QWidget()
        left_lay = QVBoxLayout(left)
        left_lay.setContentsMargins(0, 0, 0, 0)
        left_lay.setSpacing(4)
        self._geo = WaveformGeoView()
        if self._hycal_modules:
            self._geo.set_modules(self._hycal_modules)
        self._geo.moduleClicked.connect(self._on_geo_clicked)
        self._geo.setSizePolicy(QSizePolicy.Policy.Expanding,
                                QSizePolicy.Policy.Expanding)
        left_lay.addWidget(self._geo, stretch=3)
        self._wave = WaveformPlotWidget()
        left_lay.addWidget(self._wave, stretch=1)
        split.addWidget(left)

        # Right: four histograms stacked vertically (each wide, long in x).
        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(0, 0, 0, 0)
        right_lay.setSpacing(0)
        self._hist_w = {spec.attr: Hist1DWidget() for spec in _HIST_SPECS}
        for hist in self._hist_w.values():
            right_lay.addWidget(hist, stretch=1)
        split.addWidget(right)

        # Even 50/50 split between the geo+waveform column and the hist stack.
        split.setStretchFactor(0, 1)
        split.setStretchFactor(1, 1)
        split.setSizes([750, 750])
        lay.addWidget(split, stretch=1)
        return tab

    def _build_cluster_tab(self) -> QWidget:
        """HyCal heatmap + cluster panel.  Populated each event by
        ``_display_clusters``."""
        tab = QWidget()
        root = QVBoxLayout(tab)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(2)

        # Calibration warning banner — visible when no gain file has been
        # loaded (cal_factor==0 → every energize() returns 0 → no clusters).
        self._calib_warn_lbl = QLabel(
            "⚠ No HyCal calibration loaded — cluster energies will be 0 "
            "and no clusters will form.  Use File → Load HyCal calibration…")
        self._calib_warn_lbl.setStyleSheet(
            f"background:{THEME.DANGER}; color:#ffffff; "
            f"padding:4px 8px; font-weight: bold;")
        self._calib_warn_lbl.setWordWrap(True)
        self._calib_warn_lbl.setVisible(not self._hycal_calib_loaded)
        root.addWidget(self._calib_warn_lbl)

        body = QWidget()
        lay = QHBoxLayout(body)
        lay.setContentsMargins(0, 0, 0, 0)
        lay.setSpacing(4)

        self._cluster_map = HyCalClusterMap()
        if self._hycal_modules:
            self._cluster_map.set_modules(self._hycal_modules)
        self._cluster_map.setSizePolicy(QSizePolicy.Policy.Expanding,
                                        QSizePolicy.Policy.Expanding)
        lay.addWidget(self._cluster_map, stretch=3)

        self._cluster_panel = ClusterPanel()
        self._cluster_panel.clusterSelected.connect(
            self._cluster_map.set_selected_cluster)
        lay.addWidget(self._cluster_panel, stretch=1)

        root.addWidget(body, stretch=1)
        return tab

    # ---- Advanced tuning dock ------------------------------------------

    def _build_advanced_dock(self):
        """Right-side collapsible dock exposing WaveConfig + HyCalClusterConfig.

        Both configs live in prad2py bindings; changes trigger a re-run of
        the current event (``_rerun_current_event``) so the effect is
        immediate.
        """
        dock = QDockWidget("Advanced tuning", self)
        rlay = setup_tuning_dock(dock)

        # ---- WaveConfig ------------------------------------------------
        wg = QGroupBox("Waveform analyser")
        wf = QFormLayout(wg)
        self._adv_wave = add_config_rows(wf, self._wcfg, _WAVE_FIELDS,
                                         self._on_advanced_changed)

        # The peak filter is edited only in the Cut Settings dialog, so it
        # is deliberately absent from this dock.

        rlay.addWidget(wg)

        # ---- HyCalClusterConfig ---------------------------------------
        cg = QGroupBox("HyCal clustering")
        cf = QFormLayout(cg)
        self._adv_cluster: Dict[str, QWidget] = {}
        if self._hccl is not None:
            ccfg = self._hccl.get_config()
            self._adv_cluster = add_config_rows(cf, ccfg, _CLUSTER_FIELDS,
                                                self._on_advanced_changed)
            cbx = QCheckBox("corner_conn (include diagonal neighbors)")
            set_editor_value(cbx, getattr(ccfg, "corner_conn", False))
            cbx.toggled.connect(self._on_advanced_changed)
            cf.addRow(cbx)
            self._adv_cluster["corner_conn"] = cbx
        else:
            cf.addRow(QLabel("(HyCalSystem not initialized)"))
        rlay.addWidget(cg)

        rlay.addStretch(1)

        # Snapshot initial widget values (from WaveConfig + HyCalClusterConfig)
        # so "Reset to defaults" can restore them after arbitrary tuning.
        self._adv_wave_defaults = {n: editor_value(e)
                                   for n, e in self._adv_wave.items()}
        self._adv_cluster_defaults = {n: editor_value(e)
                                      for n, e in self._adv_cluster.items()}

        reset_btn = QPushButton("Reset to defaults")
        reset_btn.clicked.connect(self._reset_advanced_defaults)
        rlay.addWidget(reset_btn)

        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, dock)
        dock.hide()
        return dock

    def _reset_advanced_defaults(self):
        """Restore every Advanced-dock widget to its initial (post-load)
        value and re-run the current event once.  The peak filter lives in
        the Cut Settings dialog and has its own Reset button."""
        for editors, defaults in ((self._adv_wave, self._adv_wave_defaults),
                                  (self._adv_cluster, self._adv_cluster_defaults)):
            for name, ed in editors.items():
                set_editor_value(ed, defaults[name])
        # Single re-run after all widgets are restored.
        self._on_advanced_changed()

    def _on_advanced_changed(self, *_):
        """Push every dock value back into WaveConfig + HyCalClusterConfig,
        then schedule a debounced re-run of the current event."""
        editors_to_config(self._adv_wave, self._wcfg)
        if self._hccl is not None and self._adv_cluster:
            ccfg = self._hccl.get_config()
            editors_to_config(self._adv_cluster, ccfg)
            self._hccl.set_config(ccfg)

        self._adv_redraw_timer.start(self._adv_debounce_ms)

    def _rerun_current_event(self):
        """Re-read + re-analyse the current event with the latest config."""
        if self._current_idx >= 0:
            self._goto(self._current_idx)

    # ---- Cut Settings dialog --------------------------------------------

    def _open_cut_settings_dialog(self):
        """On Save, the new filter values replace ``self._filter``
        (preserving ``enable``, which stays under the toolbar's "apply"
        toggle) and the current event is re-analysed."""
        dlg = CutSettingsDialog(self, self._filter, self._filter_default)
        if dlg.exec() != QDialog.DialogCode.Accepted:
            return
        new_filter = dlg.result_filter()
        new_filter.enable = self._filter.enable
        self._filter = new_filter
        self._rerun_current_event()

    def _on_cut_apply_toggled(self, on: bool):
        """Toolbar "apply" checkbox → WaveformFilter.enable; re-runs the
        current event to refresh hists / geo / clusters."""
        self._filter.enable = bool(on)
        self._rerun_current_event()

    def _on_cut_show_toggled(self, on: bool):
        """Toolbar "show" checkbox: toggles cut overlays on the waveform
        plot only — a repaint, no re-analysis."""
        self._filter_show = bool(on)
        if self._current_idx >= 0 and self._reader is not None:
            self._display_waveform(self._reader.ch.fadc())
        else:
            self._wave.update()

    def _small_btn(self, text: str, slot, primary: bool = False) -> QPushButton:
        btn = QPushButton(text)
        bg = "#1f6feb" if primary else "#21262d"
        fg = "#ffffff" if primary else "#c9d1d9"
        btn.setStyleSheet(themed(
            f"QPushButton{{background:{bg};color:{fg};"
            f"border:1px solid #30363d;padding:5px 14px;"
            f"font:bold 10pt Monospace;border-radius:3px;}}"
            f"QPushButton:hover{{background:#30363d;}}"
            f"QPushButton:disabled{{background:#161b22;color:#484f58;}}"))
        btn.clicked.connect(slot)
        return btn

    def _mk_label(self, text: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setFont(QFont("Monospace", 10))
        lbl.setStyleSheet(themed("color:#c9d1d9;"))
        return lbl

    def _make_menu(self):
        mb = self.menuBar()
        mf = mb.addMenu("&File")

        a_open = QAction("Open &evio…", self)
        a_open.setShortcut("Ctrl+O")
        a_open.triggered.connect(self._open_evio_dialog)
        mf.addAction(a_open)

        a_calib = QAction("Load HyCal &calibration…", self)
        a_calib.triggered.connect(self._load_hycal_calib_dialog)
        mf.addAction(a_calib)

        self._a_save = QAction("&Save histograms as JSON…", self)
        self._a_save.setShortcut("Ctrl+S")
        self._a_save.triggered.connect(self._save_json_dialog)
        self._a_save.setEnabled(False)
        mf.addAction(self._a_save)

        mf.addSeparator()
        a_quit = QAction("&Quit", self)
        a_quit.setShortcut("Ctrl+Q")
        a_quit.triggered.connect(self.close)
        mf.addAction(a_quit)

        mv = mb.addMenu("&View")
        a_adv = self._adv_dock.toggleViewAction()
        a_adv.setText("&Advanced tuning")
        a_adv.setShortcut("Ctrl+T")
        mv.addAction(a_adv)

    # -- HyCal pipeline (PipelineBuilder) ------------------------------

    def _guess_database_dir(self) -> Optional[str]:
        """Pick a database dir for resolving recon-config-internal paths
        (e.g., 'runinfo': 'runinfo/general.json').  We prefer the
        parent of whichever full-path config the user actually pointed at
        — recon, daq, hycal-map, in that order — so the resolved paths
        stay self-consistent with the rest of their tree."""
        for p in (self._recon_cfg_path, self._daq_cfg_path,
                  self._hycal_map_path):
            if p:
                parent = Path(p).expanduser().resolve().parent
                if parent.is_dir():
                    return str(parent)
        return None

    def _build_hycal_pipeline(self, evio_path: Optional[str] = None) -> bool:
        """Wire up HyCal via prad2py.det.PipelineBuilder.

        Pulls daq / recon / runinfo / map / per-run calibration from the
        configured paths so cluster energies are populated without a
        manual ``Load HyCal calibration…`` step.  ``evio_path`` lets the
        builder pick the right run number from the file name (otherwise
        it uses the latest entry in runinfo).

        Returns True on success, False if the build raised — callers
        should then drop back to ``_hycal_init_fallback`` so the rest
        of the viewer still functions for waveform-only browsing.
        """
        if not _HAVE_PRAD2PY:
            return False
        try:
            b = prad2py.det.PipelineBuilder()
            # Anchor relative paths inside recon_config (e.g. "runinfo": …)
            # to whichever directory the configs we were given live in,
            # otherwise the builder resolves them against CWD and silently
            # drops calibration when run from outside the repo.
            db_dir = self._guess_database_dir()
            if db_dir:
                b.set_database_dir(db_dir)
            if self._daq_cfg_path:
                b.set_daq_config(self._daq_cfg_path)
            if self._recon_cfg_path:
                b.set_recon_config(self._recon_cfg_path)
            if self._hycal_map_path:
                b.set_hycal_map(self._hycal_map_path)
            if evio_path:
                b.set_run_number_from_evio(str(evio_path))
            pipeline = b.build()
        except Exception as exc:                      # noqa: BLE001
            print(f"[hycal] PipelineBuilder.build failed: "
                  f"{type(exc).__name__}: {exc}", file=sys.stderr)
            return False

        self._pipeline = pipeline
        self._hcsys = pipeline.hycal
        self._hccl = prad2py.det.HyCalCluster(pipeline.hycal)
        try:
            self._hccl.set_config(pipeline.hycal_cluster_cfg)
        except Exception as exc:                      # noqa: BLE001
            print(f"[hycal] HyCalCluster.set_config failed: "
                  f"{type(exc).__name__}: {exc}", file=sys.stderr)
        self._hc_cache.clear()
        # PipelineBuilder leaves hycal_calib_path empty when no calibration
        # file resolved (no runinfo / no energy_calib_file).  Truthy =
        # calibration actually loaded → suppress the cluster-tab warning.
        self._hycal_calib_loaded = bool(getattr(pipeline, "hycal_calib_path", "")
                                        or "")
        # Sync the Advanced dock's cluster editors to the new config so
        # the user sees the values that actually drive reconstruction.
        self._sync_advanced_dock_cluster_cfg()
        if hasattr(self, "_calib_warn_lbl"):
            self._calib_warn_lbl.setVisible(not self._hycal_calib_loaded)
        return True

    def _hycal_init_fallback(self) -> bool:
        """Fallback to a plain HyCalSystem.init when PipelineBuilder fails
        (e.g., daq_config missing).  Provides geometry-only operation: no
        calibration, so cluster energies stay at 0 until the user picks a
        file via File → Load HyCal calibration…"""
        self._pipeline = None
        if not (_HAVE_PRAD2PY and self._hycal_map_path):
            self._hcsys = None
            self._hccl = None
            return False
        try:
            sys_obj = prad2py.det.HyCalSystem()
            if not sys_obj.init(str(self._hycal_map_path)):
                print("[hycal] HyCalSystem.init returned False — "
                      "cluster tab will be empty", file=sys.stderr)
                self._hcsys = None
                self._hccl = None
                return False
            self._hcsys = sys_obj
            self._hccl = prad2py.det.HyCalCluster(sys_obj)
        except Exception as exc:                      # noqa: BLE001
            print(f"[hycal] init failed: {type(exc).__name__}: {exc}",
                  file=sys.stderr)
            self._hcsys = None
            self._hccl = None
            return False
        self._hc_cache.clear()
        self._hycal_calib_loaded = False
        return True

    def _sync_advanced_dock_cluster_cfg(self):
        """Refresh the Advanced-dock cluster editors, and the values their
        Reset restores, from _hccl's config.  Called after a pipeline
        rebuild so the dock matches the new run's recon-config defaults
        instead of stale values; a no-op before the dock exists (the
        pipeline is first built from __init__)."""
        if self._hccl is None or not getattr(self, "_adv_cluster", None):
            return
        try:
            ccfg = self._hccl.get_config()
        except Exception:
            return
        config_to_editors(ccfg, self._adv_cluster)
        self._adv_cluster_defaults = {n: editor_value(e)
                                      for n, e in self._adv_cluster.items()}

    # -- file open --

    def _open_evio_dialog(self):
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Open evio file", str(Path.cwd()),
            "evio files (*.evio *.evio.*);;All files (*)")
        if path_str:
            self.open_path(Path(path_str))

    def _load_hycal_calib_dialog(self):
        """Load HyCal per-module calibration constants (cal_factor,
        cal_base_energy, cal_non_linear) from a JSON file via
        HyCalSystem.LoadCalibration.  Success hides the "no calibration"
        warning banner on the cluster tab and re-runs the current event
        so cluster energies reflect the new calibration."""
        if self._hcsys is None:
            QMessageBox.warning(self, "HyCal system not ready",
                "Cannot load calibration — HyCalSystem.init() did not "
                "complete.  Check that hycal_map.json is present.")
            return
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Load HyCal calibration", str(Path.cwd()),
            "JSON files (*.json);;All files (*)")
        if not path_str:
            return
        try:
            nmatched = self._hcsys.load_calibration(path_str)
        except Exception as exc:             # noqa: BLE001
            QMessageBox.critical(self, "Calibration load failed",
                                 f"{type(exc).__name__}: {exc}")
            return
        if nmatched is None or int(nmatched) <= 0:
            QMessageBox.warning(self, "Calibration load failed",
                f"load_calibration returned {nmatched!r} — 0 modules "
                f"matched.  Check that the file format matches "
                f"HyCalSystem::LoadCalibration's expected schema.")
            return
        self._hycal_calib_loaded = True
        if hasattr(self, "_calib_warn_lbl"):
            self._calib_warn_lbl.setVisible(False)
        self.statusBar().showMessage(
            f"Calibration loaded: {Path(path_str).name} "
            f"({int(nmatched)} modules matched)", 5000)
        # Re-run the current event so cluster energies reflect the new
        # cal_factors.  No-op if no event has been loaded yet.
        if self._current_idx >= 0:
            self._goto(self._current_idx)

    def open_path(self, path: Path):
        err = _check_evchannel_support()
        if err:
            QMessageBox.critical(self, "prad2py issue", err)
            return
        if self._idx_thread is not None:
            QMessageBox.information(self, "Busy", "Already indexing.")
            return

        # Tear down any previous reader / hists
        self._close_reader()
        self._channels.clear()
        self._selected_key = None
        self._combo.blockSignals(True); self._combo.clear(); self._combo.blockSignals(False)
        self._index = []
        self._current_idx = -1
        self._accumulated = None
        self._clear_plots()

        # Rebuild the HyCal pipeline against this evio's run number so the
        # right per-run calibration is loaded.  Failures fall back to
        # whatever pipeline / map-only state was already in place — the
        # waveform side keeps working regardless.
        self._build_hycal_pipeline(evio_path=str(path))

        # Start indexer
        dlg = QProgressDialog(f"Indexing {path.name} …", "Cancel", 0, 100, self)
        dlg.setWindowTitle("Indexing")
        dlg.setWindowModality(Qt.WindowModality.WindowModal)
        dlg.setMinimumDuration(0)
        dlg.setAutoClose(True)
        dlg.setValue(0)
        dlg.show()
        QApplication.processEvents()

        worker = IndexerWorker(str(path), self._daq_cfg_path)

        def _on_progress(done: int, total: int):
            if total > 0:
                dlg.setMaximum(total)
                dlg.setValue(done)
            dlg.setLabelText(f"Indexing {path.name}\n"
                             f"evio events: {done:,} / {total:,}")

        def _on_thread_finished():
            dlg.close()
            self._idx_thread = None
            self._idx_worker = None

        worker.progressed.connect(_on_progress)
        self._idx_worker = worker
        self._idx_thread = start_worker_thread(
            self, worker,
            lambda res: self._on_index_done(path, res),
            lambda msg: self._on_index_failed(path, msg),
            dialog=dlg, on_thread_finished=_on_thread_finished)

    def _on_index_done(self, path: Path, res: Dict):
        self._evio_path = path
        self._index = res["index"]
        n_phys = len(self._index)
        cancelled = bool(res.get("cancelled"))
        self._accumulated = np.zeros(n_phys, dtype=bool)

        # Open reader handle for browse use
        ok, err = self._open_reader(str(path))
        if not ok:
            QMessageBox.critical(self, "Open failed",
                                 f"Indexing finished but reader open failed:\n{err}")
            return

        mode_note = ("" if self._reader.is_ra
                     else "   [sequential mode — Prev is slow]")
        self._file_lbl.setText(
            f"{path.name}   physics events: {n_phys:,}   "
            f"(evio blocks: {res['total_evio']:,})"
            + mode_note
            + ("   [indexing cancelled]" if cancelled else ""))
        self._info.setText("Select a module, then click Next to start browsing.")

        self._event_spin.setMaximum(max(0, n_phys - 1))
        self._event_spin.setValue(0)
        if n_phys > 0:
            self._prev_btn.setEnabled(True)
            self._next_btn.setEnabled(True)
            self._event_spin.setEnabled(True)
            self._batch_btn.setEnabled(True)
            self._a_save.setEnabled(True)
            self._reset_btn.setEnabled(True)
        self._total_lbl.setText(f" / {max(0, n_phys - 1):,}")

        self.statusBar().showMessage(
            f"Indexed {n_phys:,} physics events from {path.name}")

    def _on_index_failed(self, path: Path, msg: str):
        QMessageBox.critical(self, "Indexing failed", f"{path}\n\n{msg}")
        self.statusBar().showMessage(f"Failed to index {path.name}")

    # -- reader (browse handle) --

    def _open_reader(self, path: str) -> Tuple[bool, str]:
        try:
            self._reader = EvioCursor(path, self._daq_cfg_path)
            return True, ""
        except Exception as e:
            return False, f"{type(e).__name__}: {e}"

    def _close_reader(self):
        if self._reader is not None:
            try:
                self._reader.close()
            except Exception:
                pass
        self._reader = None

    # -- navigation --

    def _on_prev(self):
        if self._current_idx > 0:
            self._goto(self._current_idx - 1)
        elif self._current_idx == -1 and self._index:
            self._goto(0)

    def _on_next(self):
        if self._current_idx < len(self._index) - 1:
            self._goto(self._current_idx + 1)

    def _on_spin_jump(self):
        v = self._event_spin.value()
        if 0 <= v < len(self._index) and v != self._current_idx:
            self._goto(v)

    def _goto(self, phys_idx: int):
        if self._reader is None or not (0 <= phys_idx < len(self._index)):
            return
        ev_idx, sub_idx = self._index[phys_idx]
        # RA: jump in O(1).  Sequential: reopen on backward jumps, walk
        # forward to the target.
        try:
            self._reader.seek(ev_idx)
        except RuntimeError as e:
            self.statusBar().showMessage(str(e))
            return
        ch = self._reader.ch
        if not ch.scan():
            self.statusBar().showMessage(f"scan() failed at physics #{phys_idx}")
            return
        ch.select_event(sub_idx)
        info = ch.info()

        self._current_idx = phys_idx
        self._event_spin.blockSignals(True)
        self._event_spin.setValue(phys_idx)
        self._event_spin.blockSignals(False)

        # Apply trigger filter: skip updating hists but still show waveform
        trig_ok = trigger_ok(int(info.trigger_bits),
                             self._accept_mask, self._reject_mask)

        fadc_evt = ch.fadc()
        self._update_channel_list_from_event(fadc_evt)
        self._accumulate_and_display(fadc_evt, info, trig_ok)

    def _update_channel_list_from_event(self, fadc_evt):
        """Add any new (roc, slot, ch) seen in this event to the combo."""
        added = False
        for r in range(fadc_evt.nrocs):
            roc = fadc_evt.roc(r)
            roc_tag = int(roc.tag)
            if roc_tag not in self._roc_to_crate:
                continue
            crate = self._roc_to_crate[roc_tag]
            for s in roc.present_slots():
                slot = roc.slot(s)
                for c in slot.present_channels():
                    key = (roc_tag, s, c)
                    if key not in self._channels:
                        module = self._daq_map.get((crate, s, c))
                        self._channels[key] = _make_hists(
                            self._hist_cfg, roc_tag, s, c, module)
                        added = True
        if added:
            self._refresh_combo()
            self._geo.set_available({c.module for c in self._channels.values()
                                     if c.module})

    def _refresh_combo(self):
        """Re-populate combo from discovered channels, preserving selection."""
        prev_key = self._selected_key
        items: List[Tuple[str, Tuple[int, int, int], str]] = []  # (sort_key, key, label)
        for key, ch in self._channels.items():
            mod = ch.module or "(unmapped)"
            label = (f"{mod:<8}  roc=0x{ch.roc:02X}  s={ch.slot:>2}  "
                     f"ch={ch.channel:>2}")
            sort_key = (0 if ch.module else 1, _natural_sort_key(mod))
            items.append((sort_key, key, label))
        items.sort(key=lambda x: x[0])

        self._combo.blockSignals(True)
        self._combo.clear()
        self._combo_keys: List[Tuple[int, int, int]] = []
        sel_idx = 0
        for i, (_, key, label) in enumerate(items):
            self._combo.addItem(label)
            self._combo_keys.append(key)
            if key == prev_key:
                sel_idx = i
        if self._combo_keys:
            self._combo.setCurrentIndex(sel_idx)
            self._selected_key = self._combo_keys[sel_idx]
        self._combo.blockSignals(False)

    def _on_combo_changed(self, idx: int):
        if not hasattr(self, "_combo_keys"):
            return
        if 0 <= idx < len(self._combo_keys):
            self._selected_key = self._combo_keys[idx]
            hits = self._channels.get(self._selected_key)
            self._geo.set_selected(hits.module if hits else None)
            # Re-render: hists from cache, waveform from current event
            self._display_hists_for_selected()
            if self._current_idx >= 0 and self._reader is not None:
                fadc_evt = self._reader.ch.fadc()
                self._display_waveform(fadc_evt)

    def _on_geo_clicked(self, name: str):
        """Geo-view click: switch combo to the (first) channel for this module."""
        if not name:
            return
        for i, key in enumerate(getattr(self, "_combo_keys", [])):
            hits = self._channels.get(key)
            if hits and hits.module == name:
                self._combo.setCurrentIndex(i)   # triggers _on_combo_changed
                return
        # Module exists in geometry but hasn't been seen yet — ignore the click
        self.statusBar().showMessage(
            f"{name}: no events seen yet for this module", 2000)

    # -- accumulate + display --

    def _accumulate_and_display(self, fadc_evt, info, trig_ok: bool):
        # Analyse every channel present in the event.  Histogram fills are
        # dedup'd via self._accumulated: an event already folded in still
        # gets re-analysed for display, but isn't counted again.
        sel_peaks: List[Peak] = []
        wcfg = self._wcfg
        flt = self._filter
        sel_key = self._selected_key
        idx = self._current_idx
        do_fill = trig_ok and not _is_folded(self._accumulated, idx)

        current_vals: Dict[str, float] = {}   # module_name -> max peak integral
        module_energies: Dict[str, float] = {}  # name -> MeV (cluster tab)
        rejected_current: set = set()

        # Reset clustering state for this event.
        if self._hccl is not None:
            self._hccl.clear()

        for key, hits, peaks in _iter_channel_peaks(fadc_evt, self._channels,
                                                    wcfg):
            if key == sel_key:
                sel_peaks = peaks
            # Pick the best (highest-integral) peak that passes the
            # active filter; matches viewer_utils.h::bestPeakInWindow
            # plus the integral/height/quality cuts in PeakFilter.
            # ``best_time`` is needed for HyCalCluster.add_hit's
            # multi-pulse seed-time gating.
            best_int = 0.0
            best_time = 0.0
            any_peak = len(peaks) > 0
            for p in peaks:
                if not flt.passes(p):
                    continue
                if p.integral > best_int:
                    best_int = p.integral
                    best_time = float(p.time)
            max_int = best_int
            if max_int > 0 and hits.module:
                current_vals[hits.module] = max_int
            elif any_peak and hits.module:
                rejected_current.add(hits.module)

            # Feed the cluster tab: resolve channel → HyCal module
            # and push (module_idx, energy_MeV, time_ns).
            roc_tag, s, c = key
            crate = self._roc_to_crate.get(roc_tag)
            if self._hccl is not None and crate is not None and max_int > 0:
                mod = self._resolve_hycal_module(crate, s, c, key)
                if mod is not None:
                    energy = mod.energize(max_int)
                    if energy > 0:
                        self._hccl.add_hit(mod.index, energy, best_time)
                        module_energies[mod.name] = energy

            if do_fill:
                hits.fill_peaks(peaks, flt)

        if do_fill:
            _mark_folded(self._accumulated, idx)

        self._geo.set_rejected_current(rejected_current)
        self._geo.set_current_values(current_vals)
        # Overall occupancy only needs refreshing when hists actually changed.
        if do_fill:
            self._geo.set_overall_values(self._compute_overall_occupancy())

        if sel_key is None:
            self._set_info_line(info, peaks=None)
            self._wave.clear("(select a module to view its waveform)")
        else:
            self._set_info_line(info, peaks=sel_peaks)
        self._display_hists_for_selected()
        self._display_waveform(fadc_evt)
        self._display_clusters(module_energies)

    def _resolve_hycal_module(self, crate: int, slot: int, ch: int,
                              cache_key: Tuple[int, int, int]):
        """Look up a HyCal Module by DAQ address, caching on the ROC-tag key.
        Returns the Module or None if this channel isn't a HyCal module
        (LMS / Veto / scaler channels return None and are silently ignored)."""
        if self._hcsys is None:
            return None
        if cache_key in self._hc_cache:
            return self._hc_cache[cache_key]
        m = self._hcsys.module_by_daq(int(crate), int(slot), int(ch))
        # Filter: only HyCal PWO4/PbGlass modules contribute to clustering.
        if m is not None and not m.is_hycal():
            m = None
        self._hc_cache[cache_key] = m
        return m

    def _display_clusters(self, module_energies: Dict[str, float]):
        """Push per-event energies + reconstructed clusters to the Cluster tab."""
        if self._hccl is None:
            return
        # Push energies to the heatmap
        if module_energies:
            vmax = max(module_energies.values())
        else:
            vmax = 1.0
        self._cluster_map.set_values(module_energies)
        self._cluster_map.set_range(0.0, vmax if vmax > 0 else 1.0)

        # Reconstruct + wrap for display
        self._hccl.form_clusters()
        matched = self._hccl.reconstruct_matched()
        display: List["_DisplayCluster"] = []
        for rr in matched:
            mc = rr.cluster
            hit = rr.hit
            centre = mc.center
            cname = ""
            members = []
            try:
                cmod = self._hcsys.module(centre.index)
                cname = cmod.name
            except Exception:
                pass
            for h in mc.hits:
                try:
                    members.append(self._hcsys.module(h.index).name)
                except Exception:
                    pass
            display.append(_DisplayCluster(hit, cname, members))

        self._cluster_map.set_clusters(display)
        self._cluster_panel.set_clusters(display)

    def _compute_overall_occupancy(self) -> Dict[str, float]:
        """module_name -> events_with_peak / events_accumulated (skip empty)."""
        out: Dict[str, float] = {}
        for hits in self._channels.values():
            if hits.module and hits.events > 0:
                out[hits.module] = hits.peak_events / hits.events
        return out

    def _set_info_line(self, info, peaks: Optional[List[Peak]]):
        pieces = [
            f"event #{self._current_idx:,}",
            f"tb=0x{int(info.trigger_bits):X}",
            f"evnum={int(info.event_number)}",
        ]
        if peaks is not None:
            pieces.append(f"peaks(this view)={len(peaks)}")
        key = self._selected_key
        if key:
            hits = self._channels.get(key)
            if hits:
                pieces.append(f"accum={hits.events:,}")
                pieces.append(f"w/peak={hits.peak_events:,}")
        self._info.setText("   ".join(pieces))

    def _display_hists_for_selected(self):
        key = self._selected_key
        hits = self._channels.get(key) if key else None
        if hits is None:
            self._clear_hists()
            return
        mod = hits.module or "(unmapped)"
        for spec in _HIST_SPECS:
            h = getattr(hits, spec.attr)
            self._hist_w[spec.attr].set_data(
                h.bins, h.bmin, h.bstep, under=h.under, over=h.over,
                title=f"{mod}  —  {spec.title}", color=spec.color)

    def _display_waveform(self, fadc_evt):
        key = self._selected_key
        if not key:
            self._wave.clear("(select a module)")
            return
        roc_tag, slot, ch = key
        stack_key = f"{roc_tag:02X}_{slot}_{ch}"
        samples = _find_channel_samples(fadc_evt, roc_tag, slot, ch)
        if samples is None or samples.size == 0:
            # Keep the stacker intact on empty events (as in
            # resources/waveform.js) but still clear it when the user has
            # switched to a different module.
            if self._wave.is_stacking():
                self._wave.reset_stack_if_new_key(stack_key)
                return
            hits = self._channels.get(key)
            mod = hits.module if hits and hits.module else "(unmapped)"
            self._wave.clear(f"{mod} not present in event #{self._current_idx}")
            return
        ped_mean, ped_rms, peaks = analyze(samples, self._wcfg)
        hits = self._channels.get(key)
        mod = hits.module if hits and hits.module else "(unmapped)"
        self._wave.set_data(samples, peaks, ped_mean, ped_rms,
                            title=(f"{mod}   roc=0x{roc_tag:02X}  "
                                   f"slot={slot}  ch={ch}"),
                            clk_mhz=self._wcfg.clk_mhz,
                            stack_key=stack_key,
                            filter=self._filter,
                            filter_show=self._filter_show)

    def _clear_hists(self):
        for spec in _HIST_SPECS:
            self._hist_w[spec.attr].clear(spec.label)

    def _clear_plots(self):
        self._clear_hists()
        self._wave.clear("(open an evio file and click Next)")

    def _reset_current_hists(self):
        key = self._selected_key
        if not key: return
        hits = self._channels.get(key)
        if not hits: return
        for spec in _HIST_SPECS:
            getattr(hits, spec.attr).reset()
        hits.events = 0
        hits.peak_events = 0
        self._display_hists_for_selected()
        self._geo.set_overall_values(self._compute_overall_occupancy())
        self.statusBar().showMessage(
            f"Reset histograms for {hits.module or '(unmapped)'}")

    # -- batch 10k --

    def _on_batch_10k(self):
        if self._batch_thread is not None:
            QMessageBox.information(self, "Busy", "Batch already running.")
            return
        start_idx = max(0, self._current_idx + 1)
        if start_idx >= len(self._index):
            QMessageBox.information(self, "End of file",
                                    "Already at the last physics event.")
            return
        remaining = len(self._index) - start_idx
        count = min(10_000, remaining)

        worker = BatchWorker(
            evio_path=str(self._evio_path),
            daq_config_path=self._daq_cfg_path,
            index=self._index, start_idx=start_idx, count=count,
            channels=self._channels, wcfg=self._wcfg,
            peak_filter=self._filter,
            accept_mask=self._accept_mask, reject_mask=self._reject_mask,
            accumulated=self._accumulated,
        )

        # Modal progress dialog — blocks input to the main window until the
        # batch finishes (or the user cancels), so they can't switch modules
        # / reload / Prev / Next while hists are being filled underneath.
        dlg = QProgressDialog(
            f"Processing {count:,} events…", "Cancel", 0, count, self)
        dlg.setWindowTitle("Accumulating")
        dlg.setWindowModality(Qt.WindowModality.WindowModal)
        dlg.setMinimumDuration(0)
        dlg.setAutoClose(True)
        dlg.setAutoReset(False)
        dlg.setValue(0)
        dlg.show()
        QApplication.processEvents()

        def _on_progress(done: int, target: int, peaks: int):
            dlg.setValue(done)
            dlg.setLabelText(
                f"Processing {target:,} events\n"
                f"done: {done:,} / {target:,}   peaks found: {peaks:,}")
            self._batch_status.setText(
                f"batch: {done:,}/{target:,}  peaks={peaks:,}")
            # refresh hists + geo overall map incrementally
            self._display_hists_for_selected()
            self._geo.set_overall_values(self._compute_overall_occupancy())

        def _on_finished(n: int):
            self._current_idx = start_idx + n - 1
            self._event_spin.blockSignals(True)
            self._event_spin.setValue(max(0, self._current_idx))
            self._event_spin.blockSignals(False)
            self._batch_status.setText(
                f"batch done: {n:,} events processed")
            self._display_hists_for_selected()
            self._geo.set_overall_values(self._compute_overall_occupancy())
            # advance one more to show the next waveform
            if self._current_idx + 1 < len(self._index):
                self._goto(self._current_idx + 1)

        def _on_failed(msg: str):
            QMessageBox.critical(self, "Batch failed", msg)
            self._batch_status.setText("batch failed")

        def _cleanup():
            dlg.close()
            self._batch_btn.setEnabled(True)
            self._batch_worker = None
            self._batch_thread = None

        worker.progressed.connect(_on_progress)
        self._batch_worker = worker
        self._batch_btn.setEnabled(False)
        self._batch_status.setText(f"batch: 0/{count:,}")
        self._batch_thread = start_worker_thread(
            self, worker, _on_finished, _on_failed,
            dialog=dlg, on_thread_finished=_cleanup)

    # -- JSON save --

    def _save_json_dialog(self):
        if not self._channels:
            return
        default = (self._evio_path.name + ".waveform.json"
                   if self._evio_path else "waveform_hist.json")
        path_str, _ = QFileDialog.getSaveFileName(
            self, "Save histograms as JSON", str(Path.cwd() / default),
            "JSON files (*.json)")
        if not path_str:
            return
        try:
            out = {
                "source_file": str(self._evio_path) if self._evio_path else "",
                **{spec.json_key: {k: self._hist_cfg[spec.attr][k]
                                   for k in ("min", "max", "step")}
                   for spec in _HIST_SPECS},
                "filter":        self._filter.to_json(),
                "filter_active": self._filter.enable,
                "wave_config":   self._wcfg.__dict__.copy(),
                "channels":      {},
            }
            for (roc, slot, ch), hits in sorted(self._channels.items()):
                out["channels"][f"{roc}_{slot}_{ch}"] = {
                    "module":      hits.module,
                    "roc":         hits.roc,
                    "slot":        hits.slot,
                    "channel":     hits.channel,
                    "events":      hits.events,
                    "peak_events": hits.peak_events,
                    **{spec.json_key: getattr(hits, spec.attr).to_json()
                       for spec in _HIST_SPECS},
                }
            Path(path_str).parent.mkdir(parents=True, exist_ok=True)
            with open(path_str, "w", encoding="utf-8") as f:
                json.dump(out, f)
            self.statusBar().showMessage(f"Saved {path_str}")
        except Exception as ex:
            QMessageBox.warning(self, "Save failed", f"{path_str}\n\n{ex}")

    # -- close --

    def closeEvent(self, ev):
        if self._idx_worker is not None:
            self._idx_worker.request_cancel()
        if self._batch_worker is not None:
            self._batch_worker.request_cancel()
        for thr in (self._idx_thread, self._batch_thread):
            if thr is not None and thr.isRunning():
                thr.quit()
                thr.wait(3000)
        self._close_reader()
        super().closeEvent(ev)


# ---- Main ----

def main():
    ap = argparse.ArgumentParser(
        description="HyCal Event Viewer — browse an evio file event-by-event "
                    "with Waveform and Cluster tabs.")
    ap.add_argument("path", nargs="?", type=Path,
                    help="evio file to open (otherwise use File → Open…).")
    ap.add_argument("--config", type=Path,
                    default=_REPO_DIR / "database" / "monitor_config.json",
                    help="monitor_config.json (waveform binning).")
    ap.add_argument("--daq-config", type=Path,
                    default=_REPO_DIR / "database" / "daq_config.json",
                    help="daq_config.json (ROC-tag → crate mapping).")
    ap.add_argument("--hycal-map", type=Path,
                    default=_REPO_DIR / "database" / "hycal_map.json",
                    help="hycal_map.json (module geometry + DAQ map).")
    ap.add_argument("--recon-config", type=Path,
                    default=_REPO_DIR / "database" / "reconstruction_config.json",
                    help="reconstruction_config.json (runinfo pointer + "
                         "HyCal cluster config).  Resolved by PipelineBuilder "
                         "to load the per-run HyCal calibration automatically.")
    ap.add_argument("--trigger-bits", type=Path,
                    default=_REPO_DIR / "database" / "trigger_bits.json",
                    help="trigger_bits.json (for --accept/--reject-trigger).")
    ap.add_argument("--accept-trigger", action="append", default=[],
                    metavar="NAME",
                    help="Require at least one of these trigger bits (repeatable).")
    ap.add_argument("--reject-trigger", action="append", default=None,
                    metavar="NAME",
                    help="Drop events with any of these trigger bits (repeatable). "
                         "Default: uses monitor_config.json setting.")
    ap.add_argument("--theme", choices=available_themes(), default="dark",
                    help="Colour theme (default: dark)")
    args = ap.parse_args()

    set_theme(args.theme)

    hist_cfg      = load_hist_config(args.config)      if args.config.is_file()      else {}
    roc_to_crate  = load_roc_tag_map(args.daq_config)  if args.daq_config.is_file()  else {}
    daq_map       = load_daq_map(args.hycal_map)       if args.hycal_map.is_file()   else {}
    bit_map       = load_trigger_bit_map(args.trigger_bits)
    hycal_modules = (load_geo_modules(args.hycal_map)
                     if args.hycal_map.is_file() else [])

    accept_names = args.accept_trigger or hist_cfg.get("accept_trigger_bits", []) or []
    reject_names = (args.reject_trigger if args.reject_trigger is not None
                    else hist_cfg.get("reject_trigger_bits", []) or [])
    accept_mask = _mask_from_names(accept_names, bit_map) if accept_names else 0
    reject_mask = _mask_from_names(reject_names, bit_map) if reject_names else 0

    app = QApplication(sys.argv)
    win = HyCalEventViewer(
        hist_config       = hist_cfg,
        daq_map           = daq_map,
        roc_to_crate      = roc_to_crate,
        accept_mask       = accept_mask,
        reject_mask       = reject_mask,
        daq_config_path   = str(args.daq_config) if args.daq_config.is_file() else "",
        hycal_modules     = hycal_modules,
        hycal_map_path    = (str(args.hycal_map)
                             if args.hycal_map.is_file() else None),
        recon_config_path = (str(args.recon_config)
                             if args.recon_config.is_file() else None),
    )
    win.show()
    if args.path is not None:
        QTimer.singleShot(0, lambda: win.open_path(args.path))
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
