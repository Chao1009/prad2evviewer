#!/usr/bin/env python3
"""
GEM Event Viewer (PyQt6) — step through an EVIO file event-by-event,
running full GEM reconstruction (pedestal + CM + ZS + clustering) in
process via ``prad2py``.

Features:
  * Pre-scans the file on open to build an event index (progress dialog).
  * Prev / Next / Goto event# + slider for navigation.
  * Collapsible "Advanced tuning" dock with every GemSystem / GemCluster
    knob (ZS σ, CM threshold, min-cluster-hits, clustering, XY matching) —
    values are live; each change re-runs reconstruction on cached SSP
    data (no EVIO I/O).

Usage:
    python gem/gem_event_viewer.py [file.evio.00000] [-D daq_config.json]
        [-G gem_map.json] [-P gem_ped.txt] [--theme THEME]

If an EVIO path is given on the command line the viewer starts scanning
it immediately; otherwise use File → Open EVIO.

Export mode (render PNGs and exit, no GUI):
    python gem/gem_event_viewer.py file.evio --event N | --events SPEC
    python gem/gem_event_viewer.py --layout [--show-every K]
    python gem/gem_event_viewer.py --json FILE|DIR|GLOB ...
  with -o OUT, --det N, --width W, --height H, --verbose.
"""

from __future__ import annotations

import argparse
import os
import sys
import tempfile
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_SCRIPT_DIR = Path(__file__).resolve().parent
# Shared helpers live in scripts/, gem/'s sibling in the source and install
# trees.  A missing import here means the install is broken.
_SCRIPTS_DIR = _SCRIPT_DIR.parent / "scripts"
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))
from prad2_env import find_database_file, import_prad2py  # noqa: E402

prad2py, PRAD2PY_ERROR = import_prad2py()
HAVE_PRAD2PY = prad2py is not None
det = prad2py.det if prad2py else None


from PyQt6.QtCore import (  # noqa: E402
    QObject, QPointF, QRectF, QSize, Qt, QThread, QTimer, pyqtSignal,
)
from PyQt6.QtGui import (  # noqa: E402
    QAction, QColor, QFont, QImage, QKeySequence,
    QPainter, QPen,
)
from PyQt6.QtWidgets import (  # noqa: E402
    QApplication,
    QCheckBox,
    QComboBox,
    QDockWidget,
    QDoubleSpinBox,
    QFileDialog,
    QFormLayout,
    QFrame,
    QGridLayout,
    QGroupBox,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QMessageBox,
    QProgressDialog,
    QPushButton,
    QScrollArea,
    QSizePolicy,
    QSlider,
    QSpinBox,
    QStatusBar,
    QTabWidget,
    QVBoxLayout,
    QWidget,
)

import numpy as np  # noqa: E402

# Sibling imports — this file lives in gem/ alongside the helpers.
# gem_view imports gem_strip_map which requires prad2py.det, so wrap so
# the GUI can still *start* and show an error dialog when missing.
try:
    from gem_view import (  # noqa: E402
        build_apv_map,
        build_det_list_from_gemsys,
        build_strip_layout,
        build_zs_apvs_from_gemsys,
        draw_event_panels,
        draw_layout,
        load_gem_map,
        process_zs_hits,
    )
except Exception as _sib_exc:  # noqa: BLE001
    build_strip_layout = load_gem_map = None  # type: ignore
    build_apv_map = build_det_list_from_gemsys = build_zs_apvs_from_gemsys = None  # type: ignore
    draw_event_panels = draw_layout = process_zs_hits = None  # type: ignore
    if HAVE_PRAD2PY:
        PRAD2PY_ERROR = (PRAD2PY_ERROR + "\n" if PRAD2PY_ERROR else "") + \
                        f"sibling import: {type(_sib_exc).__name__}: {_sib_exc}"
    HAVE_PRAD2PY = False

from evio_io import EvioCursor, iter_physics_records, open_evio  # noqa: E402
from hycal_geoview import (  # noqa: E402
    THEME, add_config_rows, apply_theme_palette, available_themes,
    editors_to_config, set_editor_value, set_theme, setup_tuning_dock,
    start_worker_thread, themed,
)


# ---- Event index ----


@dataclass
class EventMeta:
    """Metadata for one physics sub-event in the source EVIO file.

    ``record_idx`` / ``subevt_idx`` locate the event in the file's
    record structure — needed by Stepper to re-read it.  Identity
    fields (event_number / trigger_*) come from the EventInfo bank and
    are shown in the UI.
    """
    record_idx: int
    subevt_idx: int
    event_number: int
    trigger_number: int
    trigger_bits: int


def _build_event_index(path: str, daq_config_path: str, progress=None,
                       cancel=None) -> Tuple[List[EventMeta], int]:
    """EventMeta for every Physics sub-event in ``path``, plus the number
    of records in the file (records walked, in sequential mode).

    ``progress(n_events, n_records)`` is called every
    ``ScanWorker.PROGRESS_EVERY`` events; ``cancel()`` stops the walk.
    Raises RuntimeError if the file cannot be opened.
    """
    ch, is_ra = open_evio(path, daq_config_path)
    events: List[EventMeta] = []
    n_walked = 0

    def _on_record(idx: int):
        nonlocal n_walked
        n_walked = idx + 1

    try:
        n_evio = ch.get_random_access_event_count() if is_ra else 0
        for rec in iter_physics_records(ch, is_ra, cancel, _on_record):
            for i in range(ch.get_n_events()):
                ch.select_event(i)
                info = ch.info()
                events.append(EventMeta(
                    record_idx=rec,
                    subevt_idx=i,
                    event_number=int(info.event_number),
                    trigger_number=int(info.trigger_number),
                    trigger_bits=int(info.trigger_bits),
                ))
                if progress is not None and \
                        len(events) % ScanWorker.PROGRESS_EVERY == 0:
                    progress(len(events), rec + 1)
    finally:
        ch.close()
    return events, (n_evio if is_ra else n_walked)


class ScanWorker(QObject):
    """Builds an EventMeta list for every Physics sub-event in a file.

    Runs on a QThread so the UI stays responsive.  Emits ``progress``
    every ``PROGRESS_EVERY`` physics events seen and ``finished`` once
    done.  Cancel via ``request_cancel()``.
    """

    PROGRESS_EVERY = 2000

    progress = pyqtSignal(int, int)           # (physics_seen, records_seen)
    finished = pyqtSignal(object, float)      # (List[EventMeta], elapsed_seconds)
    failed   = pyqtSignal(str)

    def __init__(self, path: str, daq_config_path: str):
        super().__init__()
        self._path = path
        self._daq = daq_config_path
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        start = time.monotonic()
        try:
            events, n_rec = _build_event_index(
                self._path, self._daq, self.progress.emit,
                lambda: self._cancel)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(f"{type(exc).__name__}: {exc}")
            return
        self.progress.emit(len(events), n_rec)
        self.finished.emit(events, time.monotonic() - start)


# ---- Pedestal generation (det.GemPedestal, shared with gem_dump -m ped) ----


class PedestalWorker(QObject):
    """Builds per-strip pedestals from up to ``max_events`` events of
    ``path`` that carry full-readout APVs (online-ZS APVs are skipped).
    Writes the APV-block pedestal text file (det.GemPedestal.write) to
    ``output_path``.
    """

    PROGRESS_EVERY = 50

    progress = pyqtSignal(int, int)           # (done, target)
    finished = pyqtSignal(str, int, int)      # (output_path, napvs, events_used)
    failed   = pyqtSignal(str)

    def __init__(self, path: str, daq_config_path: str,
                 output_path: str, max_events: int = 1000):
        super().__init__()
        self._path = path
        self._daq = daq_config_path
        self._out = output_path
        self._max = max_events
        self._cancel = False

    def request_cancel(self):
        self._cancel = True

    def run(self):
        try:
            ch, is_ra = open_evio(self._path, self._daq)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(f"{type(exc).__name__}: {exc}")
            return

        ped    = det.GemPedestal()
        n_used = 0
        target = self._max

        try:
            for _ in iter_physics_records(
                    ch, is_ra, lambda: self._cancel or n_used >= target):
                for i in range(ch.get_n_events()):
                    if n_used >= target:
                        break
                    ch.select_event(i)
                    if ped.accumulate(ch.gem()) == 0:
                        continue
                    n_used += 1
                    if n_used % self.PROGRESS_EVERY == 0:
                        self.progress.emit(n_used, target)
        except Exception as exc:  # noqa: BLE001
            self.failed.emit(f"{type(exc).__name__}: {exc}")
            return
        finally:
            ch.close()
        self.progress.emit(n_used, target)

        if n_used == 0:
            self.failed.emit("no full-readout events found in file")
            return
        napvs = ped.write(self._out)
        if napvs < 0:
            self.failed.emit(f"write failed: {self._out}")
            return
        self.finished.emit(self._out, napvs, n_used)


# ---- Stepper — reads one event by index from the EVIO file ----


class Stepper:
    """Positioned EVIO reader.  Fetches SspEventData for any event by its
    evio record index.

    Uses random-access mode when the file supports it — Prev/Next/Jump
    are then O(1).  For evio-4 files without a RA-friendly block index
    falls back to sequential mode, where backward jumps cost a
    close/reopen + walk.  A small LRU cache keeps repeat visits instant
    in both modes.
    """

    CACHE_SIZE = 16

    def __init__(self, path: str, daq_config_path: str):
        self._path = path
        self._daq = daq_config_path
        self._cur: Optional[EvioCursor] = None
        self._cache: Dict[int, object] = {}  # event_idx -> SspEventData
        self._cache_order: List[int] = []

    # --- lifecycle -------------------------------------------------------

    def open(self):
        self._cur = EvioCursor(self._path, self._daq)

    def close(self):
        if self._cur is not None:
            self._cur.close()
            self._cur = None

    def is_random_access(self) -> bool:
        return self._cur is not None and self._cur.is_ra

    # --- fetch -----------------------------------------------------------

    def get_ssp(self, evmeta: EventMeta, event_idx: int):
        """Return the SspEventData for ``evmeta`` (event_idx is the ordinal
        into the EventMeta list, used as cache key)."""
        if event_idx in self._cache:
            self._cache_order.remove(event_idx)
            self._cache_order.append(event_idx)
            return self._cache[event_idx]

        if self._cur is None:
            self.open()
        self._cur.seek(evmeta.record_idx)
        ch = self._cur.ch
        if not ch.scan():
            raise RuntimeError(f"scan failed on record {evmeta.record_idx}")

        ch.select_event(evmeta.subevt_idx)
        ssp = ch.gem()

        # LRU cache insert (keep CACHE_SIZE most recent SSP payloads).
        self._cache[event_idx] = ssp
        self._cache_order.append(event_idx)
        while len(self._cache_order) > self.CACHE_SIZE:
            old = self._cache_order.pop(0)
            self._cache.pop(old, None)

        return ssp


# ---- GEM event canvas — native QPainter, no matplotlib ----


def _render_png(path: str, width: int, height: int, draw_fn,
                *args, **kwargs) -> bool:
    """Draw ``draw_fn(painter, rect, *args, **kwargs)`` into a fresh
    ``width`` x ``height`` image and save it as PNG.  Always dark on
    white for printed output, regardless of the GUI theme."""
    image = QImage(width, height, QImage.Format.Format_ARGB32)
    image.fill(QColor("white"))
    p = QPainter(image)
    try:
        draw_fn(p, QRectF(image.rect()), *args,
                bg=QColor("white"), fg=QColor("#222"), **kwargs)
    finally:
        p.end()
    return image.save(path, "PNG")


class GemEventCanvas(QWidget):
    """Custom widget that renders multi-detector event views via
    ``gem_view.draw_event_panels``.  Stores the last rendered payload so
    it can re-paint on resize and export to PNG on demand."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(600, 300)
        self.setAutoFillBackground(False)
        self._payload: Optional[Tuple[dict, list, dict, Optional[dict], Optional[str], int]] = None
        self._bg = QColor(getattr(THEME, "BG", "white"))
        self._fg = QColor(getattr(THEME, "TEXT", "#222"))

    def set_event(self, detectors, det_list, det_hits, hole,
                  *, title=None, det_filter=-1):
        self._payload = (detectors, det_list, det_hits, hole, title, det_filter)
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        try:
            self._render(p, QRectF(self.rect()))
        finally:
            p.end()

    def _render(self, p: QPainter, rect: QRectF):
        if self._payload is None or draw_event_panels is None:
            p.fillRect(rect, self._bg)
            p.setPen(self._fg)
            p.setFont(self.font())
            p.drawText(rect, Qt.AlignmentFlag.AlignCenter,
                       "(open an EVIO file — File → Open EVIO…)")
            return
        detectors, det_list, det_hits, hole, title, det_filter = self._payload
        draw_event_panels(p, rect, detectors, det_list, det_hits, hole,
                          title=title, det_filter=det_filter,
                          bg=self._bg, fg=self._fg)

    def save_png(self, path: str, *, width: int = 2400, height: int = 900) -> bool:
        """Render the current event into a PNG file."""
        if self._payload is None:
            return False
        detectors, det_list, det_hits, hole, title, det_filter = self._payload
        return _render_png(path, width, height, draw_event_panels,
                           detectors, det_list, det_hits, hole,
                           title=title, det_filter=det_filter)


# ---- Raw APV view ----


class ApvPanel(QWidget):
    """Mini per-APV panel — 128 channels × 6 time samples drawn as 6 line
    traces (blue → red by time sample), ZS-survivor channels marked at the
    bottom, diagnostic badge in the title bar.  Data lives in a single
    (128, 6) float32 numpy array — caller sets it via ``set_frame``."""

    MIN_W = 180
    MIN_H = 110
    HINT_W = 260
    HINT_H = 200
    TITLE_H = 16
    HIT_ROW_H = 6

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMinimumSize(self.MIN_W, self.MIN_H)
        # Horizontal stretch to fill the grid cell (capped by RawApvTab
        # to ≤ viewport/COLS); height is locked via setFixedHeight so it
        # stays constant across filter toggles regardless of the layout.
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Fixed)
        self.setFixedHeight(self.HINT_H)
        self._title = ""
        self._badge = ""                # e.g. "no hits" for full-readout APVs with no ZS survivors
        self._frame: Optional[np.ndarray] = None    # (128, 6) float or int16
        self._hits:  Optional[np.ndarray] = None    # (128,) bool
        self._y_lo = 0.0
        self._y_hi = 1.0
        # Display options set by RawApvTab on every refresh.
        self._sample_mask: Tuple[bool, ...] = (True,) * 6
        self._y_fixed: Optional[Tuple[float, float]] = None
        self._thr_trace:  Optional[np.ndarray] = None    # (128,) float
        self._cm_trace:   Optional[np.ndarray] = None    # (6,)  int16
        self._signal_flag = False       # True → accent border (ZS survivors present)

    def sizeHint(self):
        return QSize(self.HINT_W, self.HINT_H)

    def set_frame(self, title: str, frame: np.ndarray,
                  hits: Optional[np.ndarray] = None,
                  badge: str = "",
                  *,
                  sample_mask: Tuple[bool, ...] = (True,) * 6,
                  y_fixed: Optional[Tuple[float, float]] = None,
                  thr_trace: Optional[np.ndarray] = None,
                  cm_trace: Optional[np.ndarray] = None,
                  signal_flag: bool = False):
        self._title = title
        self._badge = badge
        self._frame = frame
        self._hits  = hits
        self._sample_mask = sample_mask
        self._y_fixed = y_fixed
        self._thr_trace = thr_trace
        self._cm_trace  = cm_trace
        self._signal_flag = signal_flag

        self._y_lo, self._y_hi = (
            y_fixed if y_fixed is not None
            else self.compute_fixed_range({0: frame}, sample_mask))
        self.update()

    @staticmethod
    def compute_fixed_range(frames: Dict[int, np.ndarray],
                            sample_mask: Tuple[bool, ...]) -> Tuple[float, float]:
        """Padded span over every enabled (strip, ts) value in ``frames``
        (one frame for a panel's own auto-scale, all of them for a shared
        Y scale)."""
        lo, hi = float("inf"), float("-inf")
        use_idx = [i for i, on in enumerate(sample_mask) if on]
        if not use_idx:
            return 0.0, 1.0
        for f in frames.values():
            if f is None or f.size == 0:
                continue
            v = f[:, use_idx]
            if v.size == 0: continue
            lo = min(lo, float(np.min(v)))
            hi = max(hi, float(np.max(v)))
        if not np.isfinite(lo) or not np.isfinite(hi):
            return 0.0, 1.0
        if hi - lo < 8.0:
            mid = 0.5 * (lo + hi)
            lo, hi = mid - 4.0, mid + 4.0
        pad = 0.08 * (hi - lo)
        return lo - pad, hi + pad

    def paintEvent(self, _ev):
        p = QPainter(self)
        try:
            self._paint(p)
        finally:
            p.end()

    def _paint(self, p: QPainter):
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        w, h = self.width(), self.height()

        # Canvas: slightly softer than THEME.BG so panels read as
        # inset plot tiles rather than sitting flush with the window.
        bg = QColor(THEME.BG_SUBTLE)
        fg = QColor(getattr(THEME, "TEXT", "#c9d1d9"))
        dim = QColor(getattr(THEME, "TEXT_DIM", "#8b949e"))
        p.fillRect(0, 0, w, h, bg)

        # Frame border priority: badge (red) > signal_flag (accent) > default.
        # Badge warns about "no hits" full-readout APVs; signal_flag
        # highlights APVs with surviving ZS hits so the eye can spot
        # them at a glance when the Signal Only filter is off.
        border = QColor(getattr(THEME, "BORDER", "#30363d"))
        badge_col = QColor(getattr(THEME, "DANGER", "#ff6b6b"))
        accent_col = QColor(getattr(THEME, "ACCENT", "#ffd166"))
        border_w = 1
        if self._badge:
            border = badge_col
            border_w = 1
        elif self._signal_flag:
            border = accent_col
            border_w = 2
        p.setPen(QPen(border, border_w))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(0, 0, w - 1, h - 1)

        # Title
        p.setPen(fg)
        p.setFont(QFont("Monospace", 8, QFont.Weight.Bold))
        title_rect = QRectF(4, 2, w - 8, self.TITLE_H - 2)
        p.drawText(title_rect, Qt.AlignmentFlag.AlignLeft
                   | Qt.AlignmentFlag.AlignVCenter, self._title)
        if self._badge:
            p.setPen(badge_col)
            p.drawText(title_rect, Qt.AlignmentFlag.AlignRight
                       | Qt.AlignmentFlag.AlignVCenter, self._badge)

        # Plot area
        plot = QRectF(4, self.TITLE_H + 2,
                      w - 8, h - self.TITLE_H - self.HIT_ROW_H - 4)
        if self._frame is None or self._frame.size == 0:
            p.setPen(dim)
            p.setFont(QFont("Monospace", 8))
            p.drawText(plot, Qt.AlignmentFlag.AlignCenter, "(no data)")
            return

        # Axes + zero line
        p.setPen(QPen(dim, 0, Qt.PenStyle.DotLine))
        if self._y_lo < 0 < self._y_hi:
            zy = plot.bottom() - (0 - self._y_lo) / (self._y_hi - self._y_lo) * plot.height()
            p.drawLine(QPointF(plot.left(), zy), QPointF(plot.right(), zy))

        n_strips = self._frame.shape[0]
        n_ts     = self._frame.shape[1]
        span_y = max(self._y_hi - self._y_lo, 1e-6)
        step_x = plot.width() / max(n_strips - 1, 1)

        def to_y(v: float) -> float:
            return plot.bottom() - (v - self._y_lo) / span_y * plot.height()

        # Threshold curve: per-channel ZS cut (ped.noise × zero_sup_thres).
        # Drawn as a dashed grey line; also mirrored at -threshold so the
        # reader can see both ±nσ bands that bracket the zero line.
        if self._thr_trace is not None and self._thr_trace.size == n_strips:
            pen = QPen(QColor(getattr(THEME, "TEXT_DIM", "#8b949e")), 0.8)
            pen.setStyle(Qt.PenStyle.DashLine)
            p.setPen(pen)
            prev_hi = prev_lo = None
            for ch in range(n_strips):
                x = plot.left() + ch * step_x
                t = float(self._thr_trace[ch])
                yhi = to_y(+t)
                ylo = to_y(-t)
                if prev_hi is not None:
                    p.drawLine(prev_hi, QPointF(x, yhi))
                    p.drawLine(prev_lo, QPointF(x, ylo))
                prev_hi = QPointF(x, yhi)
                prev_lo = QPointF(x, ylo)

        # Colored time-sample traces — blue (t=0) → red (t=5).  Time
        # samples hidden by the sample-mask checkboxes are skipped.
        for ts in range(n_ts):
            if not self._sample_mask[ts]:
                continue
            frac = ts / max(n_ts - 1, 1)
            col = QColor.fromHsvF(0.66 * (1.0 - frac), 0.85, 0.95)
            pen = QPen(col, 0.9)
            p.setPen(pen)
            prev = None
            for ch in range(n_strips):
                x = plot.left() + ch * step_x
                y = to_y(float(self._frame[ch, ts]))
                if prev is not None:
                    p.drawLine(prev, QPointF(x, y))
                prev = QPointF(x, y)

        # CM overlay — drawn AFTER data traces so it sits on top.  One
        # bold dashed line per enabled time sample, colour-matched with
        # the corresponding data trace so the user can pair firmware CM
        # with the same-colour strip waveform.  Extends across the full
        # plot width because CM is a single value for all 128 strips.
        if self._cm_trace is not None and self._cm_trace.size == n_ts:
            for ts in range(n_ts):
                if not self._sample_mask[ts]:
                    continue
                frac = ts / max(n_ts - 1, 1)
                col = QColor.fromHsvF(0.66 * (1.0 - frac), 0.6, 1.0)
                pen = QPen(col, 1.4)
                pen.setStyle(Qt.PenStyle.DashLine)
                p.setPen(pen)
                y = to_y(float(self._cm_trace[ts]))
                p.drawLine(QPointF(plot.left(), y), QPointF(plot.right(), y))

        # ZS survivor tick row (directly below the plot)
        if self._hits is not None and self._hits.any():
            row_y = h - self.HIT_ROW_H - 2
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(getattr(THEME, "ACCENT", "#ffd166")))
            for ch in range(n_strips):
                if self._hits[ch]:
                    x = plot.left() + ch * step_x
                    p.drawRect(QRectF(x - 0.8, row_y, 1.6, self.HIT_ROW_H))

        # Y-range readouts: two tiny labels tucked into the top-left and
        # bottom-left of the plot area.  Useful both when Shared Y is on
        # (see the common scale) and off (see each panel's auto scale).
        p.setPen(dim)
        p.setFont(QFont("Monospace", 7))
        fm = p.fontMetrics()
        hi_txt = self._fmt_compact(self._y_hi)
        lo_txt = self._fmt_compact(self._y_lo)
        tx = plot.left() + 2
        p.drawText(QPointF(tx, plot.top() + fm.ascent() - 1), hi_txt)
        p.drawText(QPointF(tx, plot.bottom() - 2), lo_txt)

    @staticmethod
    def _fmt_compact(v: float) -> str:
        """Short Y-axis label: integer when |v| < 1000, else 1-decimal ke."""
        if abs(v) < 1000:
            return f"{v:.0f}"
        return f"{v/1000:.1f}k"


# Per-GEM faint background tints used in the "All" sub-tab so adjacent
# detector sections read as different rows even when their headers scroll
# off-screen.  Alpha is intentionally low (~0.13) — the tint should be
# noticeable in the gaps between panels without competing with trace data.
GEM_SECTION_TINTS: List[str] = [
    "rgba(0, 180, 216, 0.13)",   # cyan
    "rgba(81, 207, 102, 0.13)",  # green
    "rgba(255, 146, 43, 0.13)",  # orange
    "rgba(204, 93, 232, 0.13)",  # purple
]


def _gem_section_tint(det_id: int) -> str:
    return GEM_SECTION_TINTS[det_id % len(GEM_SECTION_TINTS)]


class RawApvTab(QWidget):
    """Sub-tabbed APV viewer — an "All" overview tab plus one tab per GEM
    detector, each a grid of ApvPanel.  Per-event data is cached in dicts
    keyed by GemSystem APV index, filled once per event; tab switches just
    repaint."""

    COLS = 4

    def __init__(self, parent=None):
        super().__init__(parent)
        lay = QVBoxLayout(self)
        lay.setContentsMargins(4, 4, 4, 4)
        lay.setSpacing(4)

        # -- toolbar ---------------------------------------------------
        bar = QHBoxLayout()
        bar.setSpacing(8)

        # Checked: show processed (pedestal + CM + software-ZS applied).
        # Unchecked: show raw firmware samples straight from SspEventData.
        self.process_cb = QCheckBox("Process")
        self.process_cb.setChecked(True)
        self.process_cb.toggled.connect(self._on_control_changed)
        bar.addWidget(self.process_cb)

        self.zs_only_cb = QCheckBox("Signal Only")
        self.zs_only_cb.setChecked(False)
        self.zs_only_cb.toggled.connect(self._on_control_changed)
        bar.addWidget(self.zs_only_cb)

        self.fixed_y_cb = QCheckBox("Shared Y")
        self.fixed_y_cb.setChecked(True)
        self.fixed_y_cb.setToolTip(
            "Share one Y-axis range across every visible APV so traces "
            "can be compared directly.  Uncheck for per-panel auto-scale.")
        self.fixed_y_cb.toggled.connect(self._on_control_changed)
        bar.addWidget(self.fixed_y_cb)

        self.thr_line_cb = QCheckBox("Threshold")
        self.thr_line_cb.setToolTip(
            "Draw the ±(ped.noise × ZS σ) cutoff curve as a dashed grey line.")
        self.thr_line_cb.toggled.connect(self._on_control_changed)
        bar.addWidget(self.thr_line_cb)

        self.cm_overlay_cb = QCheckBox("CM overlay")
        self.cm_overlay_cb.setChecked(False)
        self.cm_overlay_cb.setToolTip(
            "Overlay the firmware-reported online_cm[6] values as dashed "
            "lines colour-matched to each time sample — cross-check "
            "against software common-mode.")
        self.cm_overlay_cb.toggled.connect(self._on_control_changed)
        bar.addWidget(self.cm_overlay_cb)

        # Time-sample mask: one checkbox per sample.  Default all on.
        bar.addWidget(QLabel("Samples:"))
        self.sample_cbs: List[QCheckBox] = []
        for t in range(6):
            cb = QCheckBox(f"t{t}")
            cb.setChecked(True)
            cb.toggled.connect(self._on_control_changed)
            bar.addWidget(cb)
            self.sample_cbs.append(cb)

        self._status = QLabel("")
        self._status.setStyleSheet(f"color:{THEME.TEXT_DIM};")
        bar.addWidget(self._status)

        bar.addStretch(1)
        lay.addLayout(bar)

        # -- sub-tabs per detector -------------------------------------
        self._tabs = QTabWidget()
        lay.addWidget(self._tabs, stretch=1)

        # Data cache populated per event.  Keys are GemSystem APV indices.
        self._apv_meta: Dict[int, Dict] = {}
        self._processed: Dict[int, np.ndarray] = {}    # apv_idx → (128,6) float32
        self._raw:       Dict[int, np.ndarray] = {}    # apv_idx → (128,6) int16
        self._hits:      Dict[int, np.ndarray] = {}    # apv_idx → (128,) bool
        self._thr:       Dict[int, np.ndarray] = {}    # apv_idx → (128,) float32
        self._cm:        Dict[int, np.ndarray] = {}    # apv_idx → (6,)  int16
        self._no_hit_apvs: set = set()
        # det_name → list of apv_index in order
        self._grouped: Dict[str, List[int]] = {}
        # det_name → {apv_idx: ApvPanel}
        self._panels: Dict[str, Dict[int, ApvPanel]] = {}
        # det_name → list of apv_index in the packed grid order
        self._sorted_idx: Dict[str, List[int]] = {}
        # det_name → QGridLayout holding the panels (for repack)
        self._grids: Dict[str, QGridLayout] = {}
        # det_name → tab index (for grey-out)
        self._tab_index_of: Dict[str, int] = {}

        # ---- "All" sub-tab parallel state ----
        # ApvPanel can have only one Qt parent, so the All tab needs its
        # own copies of every panel.  They share the same per-event data
        # dicts above and are refreshed in lock-step from
        # _refresh_all_panels.
        self._all_panels:    Dict[str, Dict[int, ApvPanel]] = {}
        self._all_grids:     Dict[str, QGridLayout] = {}

    def reset_all(self):
        self._apv_meta.clear()
        self._processed.clear()
        self._raw.clear()
        self._hits.clear()
        self._thr.clear()
        self._cm.clear()
        self._no_hit_apvs.clear()
        self._grouped.clear()
        self._panels.clear()
        self._sorted_idx.clear()
        self._grids.clear()
        self._tab_index_of.clear()
        self._all_panels.clear()
        self._all_grids.clear()
        # QTabWidget.clear() keeps the pages alive; delete them (and their
        # ApvPanels) so a rebuild does not leave the old set behind.
        while self._tabs.count():
            page = self._tabs.widget(0)
            self._tabs.removeTab(0)
            page.deleteLater()
        self._status.setText("")

    def set_apv_metadata(self, apv_meta: List[Dict]):
        """Call once per run (after geometry load): fixes the per-detector
        groupings and builds the sub-tab skeletons.  ``apv_meta`` is a list
        of dicts keyed by ``apv_index`` (GemSystem index)."""
        self.reset_all()
        self._apv_meta = {int(m["apv_index"]): m for m in apv_meta}

        for idx, m in self._apv_meta.items():
            self._grouped.setdefault(m["det_name"], []).append(idx)

        # Build the "All" overview tab first so it lands at index 0.
        self._build_all_tab()

        for det_name in sorted(self._grouped.keys()):
            page = QScrollArea()
            page.setWidgetResizable(True)
            content = QWidget()
            grid = QGridLayout(content)
            panels, sorted_apvs = self._fill_apv_grid(grid, det_name)
            grid.setRowStretch(grid.rowCount(), 1)
            page.setWidget(content)
            self._panels[det_name] = panels
            self._sorted_idx[det_name] = sorted_apvs
            self._grids[det_name] = grid
            tab_i = self._tabs.addTab(page, det_name)
            self._tab_index_of[det_name] = tab_i
        self._apply_panel_max_width()

    def _fill_apv_grid(self, grid: QGridLayout, det_name: str
                       ) -> Tuple[Dict[int, ApvPanel], List[int]]:
        """Add an ApvPanel per APV of ``det_name`` to ``grid``, sorted by
        (plane, crate, mpd, adc_ch) so hardware-adjacent APVs land next to
        each other while X/Y planes stay grouped.  Returns the panels by
        APV index and the sorted APV indices."""
        grid.setHorizontalSpacing(4)
        grid.setVerticalSpacing(4)
        panels: Dict[int, ApvPanel] = {}
        sorted_apvs = sorted(
            self._grouped[det_name],
            key=lambda i: (self._apv_meta[i]["plane_type"],
                           self._apv_meta[i]["crate_id"],
                           self._apv_meta[i]["mpd_id"],
                           self._apv_meta[i]["adc_ch"]))
        for n, idx in enumerate(sorted_apvs):
            r, c = divmod(n, self.COLS)
            panel = ApvPanel()
            m = self._apv_meta[idx]
            panel.setToolTip(
                f"crate {m['crate_id']} mpd {m['mpd_id']} adc {m['adc_ch']}  "
                f"{m['det_name']} {m['plane_type']} pos={m['det_pos']}  "
                f"(GemSystem idx {idx})")
            grid.addWidget(panel, r, c)
            panels[idx] = panel
        # Equal stretch across the COLS data columns; each panel caps
        # at viewport/COLS via the maxWidth set in resizeEvent below.
        for c in range(self.COLS):
            grid.setColumnStretch(c, 1)
        return panels, sorted_apvs

    def _build_all_tab(self):
        """Construct the "All" overview sub-tab — vertical stack of GEM
        sections, each with a faint per-detector tint and separated from
        its neighbours by a thin horizontal line.  An empty section
        (no APVs visible under Signal Only) still renders as a tinted
        empty row so the user can see the GEM is present but quiet."""
        if not self._grouped:
            return

        # Sort detectors by det_id (numeric, stable across runs).  Falls
        # back to det_name when det_id is missing.
        det_id_of = {n: int(self._apv_meta[idxs[0]].get("det_id", -1))
                     for n, idxs in self._grouped.items()}
        det_names_ordered = sorted(self._grouped.keys(),
                                   key=lambda n: (det_id_of[n], n))

        page = QScrollArea()
        page.setWidgetResizable(True)
        content = QWidget()
        outer = QVBoxLayout(content)
        outer.setContentsMargins(2, 2, 2, 2)
        outer.setSpacing(0)

        for k, det_name in enumerate(det_names_ordered):
            if k > 0:
                sep = QFrame()
                sep.setFrameShape(QFrame.Shape.HLine)
                sep.setFixedHeight(1)
                sep.setStyleSheet(
                    f"background-color:{THEME.BORDER};"
                    f"color:{THEME.BORDER};border:none;")
                outer.addWidget(sep)

            det_id = det_id_of[det_name]
            section = QFrame()
            section.setObjectName(f"_all_sec_{det_id}_{k}")
            # Object-name selector keeps the rgba tint scoped to *this*
            # section so it doesn't leak into child widgets (the panels
            # paint their own background in ApvPanel._paint).
            section.setStyleSheet(
                f"QFrame#{section.objectName()} {{"
                f" background-color:{_gem_section_tint(det_id)}; }}")
            sec_lay = QVBoxLayout(section)
            sec_lay.setContentsMargins(8, 6, 8, 8)
            sec_lay.setSpacing(4)

            n_apvs = len(self._grouped[det_name])
            tag = f"GEM {det_id}" if det_id >= 0 else det_name
            header = QLabel(f"{tag} — {det_name}   ({n_apvs} APVs)")
            header.setStyleSheet(
                f"color:{THEME.TEXT};font-weight:600;padding:2px 0;")
            sec_lay.addWidget(header)

            grid = QGridLayout()
            grid.setContentsMargins(0, 0, 0, 0)
            sec_lay.addLayout(grid)
            panels, _ = self._fill_apv_grid(grid, det_name)

            # Keep the section visible even when every panel is hidden by
            # Signal Only — the tint + header alone signal "GEM N is
            # present but quiet this event".
            section.setMinimumHeight(70)

            outer.addWidget(section)
            self._all_panels[det_name]    = panels
            self._all_grids[det_name]     = grid

        outer.addStretch(1)
        page.setWidget(content)
        self._tabs.insertTab(0, page, "All")
        self._tabs.setCurrentIndex(0)

    def resizeEvent(self, ev):
        super().resizeEvent(ev)
        self._apply_panel_max_width()

    def _apply_panel_max_width(self):
        """Cap each panel at ``viewport_width / COLS`` so filtering (hide
        via setVisible) can't make survivors balloon past 1/COLS."""
        if not self._panels and not self._all_panels:
            return
        # Use the tab widget's content area width — scroll bar reserved.
        vp_w = self._tabs.width() if self._tabs.count() else self.width()
        max_w = max(ApvPanel.MIN_W, (vp_w - 24) // self.COLS)
        for panels in self._panels.values():
            for p in panels.values():
                p.setMaximumWidth(int(max_w))
        for panels in self._all_panels.values():
            for p in panels.values():
                p.setMaximumWidth(int(max_w))

    def set_event_data(self,
                       processed: Dict[int, np.ndarray],
                       raw:       Dict[int, np.ndarray],
                       hits:      Dict[int, np.ndarray],
                       no_hit_apvs: Optional[set] = None,
                       thresholds: Optional[Dict[int, np.ndarray]] = None,
                       cm_traces:  Optional[Dict[int, np.ndarray]] = None):
        """Push per-event data; call after process_event() returns.

        ``no_hit_apvs`` — GemSystem indices to highlight as suspicious
        (red frame + 'no hits' badge).
        ``thresholds``  — per-APV (128,) float32 array of ``ped.noise ×
        zero_sup_threshold`` values, drawn as a dashed grey curve when
        the user enables the Threshold toggle.
        ``cm_traces``   — per-APV (6,) int16 array of firmware online_cm
        values, drawn as dashed lines when CM overlay is enabled."""
        self._processed = processed
        self._raw       = raw
        self._hits      = hits
        self._no_hit_apvs = no_hit_apvs or set()
        self._thr = thresholds or {}
        self._cm  = cm_traces  or {}
        self._refresh_all_panels()

    def _on_control_changed(self, *_):
        self._refresh_all_panels()

    def _refresh_section(self, panels: Dict[int, ApvPanel],
                         sorted_ids: List[int],
                         grid: Optional[QGridLayout],
                         *,
                         source: Dict[int, np.ndarray],
                         sample_mask: Tuple[bool, ...],
                         signal_only: bool,
                         shared_range: Optional[Tuple[float, float]],
                         show_thr: bool,
                         show_cm: bool) -> Tuple[int, int]:
        """Push frames into one grid of panels and repack visible ones to
        the front.  Returns (total, shown).  Used for both per-detector
        tabs and the "All" tab's per-GEM sections."""
        visible_ordered: List[int] = []
        total = 0
        for idx in sorted_ids:
            panel = panels[idx]
            m = self._apv_meta[idx]
            total += 1
            has_zs = self._hits.get(idx)
            has_any_zs = bool(has_zs is not None and has_zs.any())
            if signal_only and not has_any_zs:
                panel.setVisible(False)
                continue
            panel.setVisible(True)
            visible_ordered.append(idx)

            frame = source.get(idx)
            title = (f"c{m['crate_id']} m{m['mpd_id']} a{m['adc_ch']}  "
                     f"{m['det_name']} {m['plane_type']} p{m['det_pos']}")
            badge = "no hits" if idx in self._no_hit_apvs else ""
            # Highlight signal panels only when showing all APVs —
            # under Signal Only every visible panel has hits, so
            # highlighting would be redundant.
            signal_flag = has_any_zs and not signal_only
            panel.set_frame(
                title, frame, has_zs, badge,
                sample_mask=sample_mask,
                y_fixed=shared_range,
                thr_trace=self._thr.get(idx) if show_thr else None,
                cm_trace=self._cm.get(idx) if show_cm else None,
                signal_flag=signal_flag,
            )

        # Repack: put visible panels into the front slots in sorted
        # order so filtering doesn't leave gaps.  Hidden panels are
        # removed from the layout entirely (they'll rejoin when
        # visible again).  The columns keep equal stretch so panels
        # still scale with the window.
        if grid is not None:
            for idx, panel in panels.items():
                grid.removeWidget(panel)
            for n, idx in enumerate(visible_ordered):
                r, c = divmod(n, self.COLS)
                grid.addWidget(panels[idx], r, c)

        return total, len(visible_ordered)

    def _refresh_all_panels(self):
        processed_view = self.process_cb.isChecked()
        signal_only    = self.zs_only_cb.isChecked()
        fixed_y        = self.fixed_y_cb.isChecked()
        show_thr       = self.thr_line_cb.isChecked() and processed_view
        show_cm        = self.cm_overlay_cb.isChecked()
        sample_mask    = tuple(cb.isChecked() for cb in self.sample_cbs)

        source = self._processed if processed_view else self._raw

        # If "Shared Y" is on, compute a single (lo, hi) across every
        # visible APV in the active view and share it.
        shared_range: Optional[Tuple[float, float]] = None
        if fixed_y:
            if signal_only:
                visible = {i: f for i, f in source.items()
                           if self._hits.get(i) is not None
                           and bool(self._hits[i].any())}
            else:
                visible = source
            shared_range = ApvPanel.compute_fixed_range(visible, sample_mask)

        shown = 0
        total = 0
        tab_bar = self._tabs.tabBar()
        active_col = tab_bar.palette().color(tab_bar.foregroundRole())
        dim_col = QColor(active_col)
        dim_col.setAlpha(100)

        # Per-detector tabs — canonical panels, contribute to status counts.
        for det_name, panels in self._panels.items():
            sorted_ids = self._sorted_idx.get(det_name, list(panels.keys()))
            t, v = self._refresh_section(
                panels, sorted_ids, self._grids.get(det_name),
                source=source, sample_mask=sample_mask,
                signal_only=signal_only, shared_range=shared_range,
                show_thr=show_thr, show_cm=show_cm)
            total += t
            shown += v
            tab_i = self._tab_index_of.get(det_name)
            if tab_i is not None:
                tab_bar.setTabTextColor(
                    tab_i, dim_col if v == 0 else active_col)

        # "All" tab — parallel panels for the overview view.  Same data
        # source and panel order, separate widgets (Qt parenting is
        # single-owner).  These don't add to the status count to avoid
        # double-reporting.
        for det_name, panels in self._all_panels.items():
            sorted_ids = self._sorted_idx.get(det_name, list(panels.keys()))
            self._refresh_section(
                panels, sorted_ids, self._all_grids.get(det_name),
                source=source, sample_mask=sample_mask,
                signal_only=signal_only, shared_range=shared_range,
                show_thr=show_thr, show_cm=show_cm)

        mode = "processed" if processed_view else "raw"
        self._status.setText(f"{shown}/{total} APVs  [{mode}]")


# ---- Advanced tuning dock ----


class AdvancedDock(QDockWidget):
    """Collapsible dock with every clustering / XY-match knob exposed as
    spinboxes.  Emits ``changed`` whenever any value changes and
    ``resetRequested`` when the user clicks "Reset to defaults"."""

    changed         = pyqtSignal()
    resetRequested  = pyqtSignal()

    def __init__(self, parent=None):
        super().__init__("Advanced tuning", parent)
        layout = setup_tuning_dock(self)
        # Library defaults until a GemSystem supplies its configs; the GUI
        # still starts without prad2py and reports that itself.
        cfg = det.ClusterConfig() if det is not None else None

        def emit(*_):
            self.changed.emit()

        # --- Thresholds (live during reconstruction) ----------------
        tg = QGroupBox("Thresholds")
        tf = QFormLayout(tg)
        self.zs_slider, self.zs_spin = self._mkfloat_slider(2.0, 15.0, 5.0, 0.1)
        tf.addRow("ZS σ", self._slider_row(self.zs_slider, self.zs_spin))
        self.cm_slider, self.cm_spin = self._mkfloat_slider(5.0, 50.0, 20.0, 0.5)
        tf.addRow("CM thr", self._slider_row(self.cm_slider, self.cm_spin))
        # ClusterConfig field name -> editor
        self.cluster_editors = add_config_rows(tf, cfg, [
            ("min_cluster_hits", 1, 10, None, "", "min cluster hits"),
        ], emit)
        layout.addWidget(tg)

        # --- Clustering ----------------------------------------------
        cg = QGroupBox("Clustering")
        cf = QFormLayout(cg)
        self.cluster_editors.update(add_config_rows(cf, cfg, [
            ("max_cluster_hits", 1, 100, None, ""),
            ("consecutive_thres", 0, 10, None, ""),
            ("split_thres", 0.0, 100.0, 0.1, "", "split_thres (ADC)"),
            ("cross_talk_width", 0.0, 64.0, 0.1, "", "cross_talk_width (mm)"),
        ], emit))
        layout.addWidget(cg)

        # --- XY matching ---------------------------------------------
        xg = QGroupBox("XY matching")
        xf = QFormLayout(xg)
        match_mode = QComboBox()
        match_mode.addItems(["0 — sorted pairing", "1 — cartesian + cuts"])
        set_editor_value(match_mode, cfg.match_mode if cfg is not None else 1)
        match_mode.currentIndexChanged.connect(emit)
        xf.addRow("match_mode", match_mode)
        self.cluster_editors["match_mode"] = match_mode
        self.cluster_editors.update(add_config_rows(xf, cfg, [
            ("match_adc_asymmetry", 0.0, 1.0, 0.05, ""),
            ("match_time_diff", 0.0, 200.0, 1.0, "", "match_time_diff (ns)"),
            ("ts_period", 1.0, 100.0, 0.5, "", "ts_period (ns)"),
        ], emit))
        layout.addWidget(xg)

        layout.addStretch(1)

        # Reset button — emits resetRequested; the main window is
        # responsible for reverting widget values to the initial config
        # (gem_map.json thresholds, reconstruction_config.json clustering).
        self._reset_btn = QPushButton("Reset to defaults")
        self._reset_btn.clicked.connect(lambda: self.resetRequested.emit())
        layout.addWidget(self._reset_btn)

        self.zs_spin.valueChanged.connect(emit)
        self.cm_spin.valueChanged.connect(emit)

    @staticmethod
    def _mkfloat_slider(lo: float, hi: float, val: float, step: float
                        ) -> Tuple[QSlider, QDoubleSpinBox]:
        """Build a (slider, spinbox) pair for a float parameter — slider
        holds integer ticks of ``step`` resolution, spinbox shows the
        float value, edits on either widget reach both."""
        n_ticks = int(round((hi - lo) / step))
        slider = QSlider(Qt.Orientation.Horizontal)
        slider.setRange(0, n_ticks)
        slider.setValue(int(round((val - lo) / step)))
        spin = QDoubleSpinBox()
        spin.setRange(lo, hi); spin.setSingleStep(step); spin.setDecimals(2)
        spin.setValue(val)
        spin.setMaximumWidth(90)

        def s2b(tick: int):
            spin.blockSignals(True)
            spin.setValue(lo + tick * step)
            spin.blockSignals(False)

        def b2s(v: float):
            slider.blockSignals(True)
            slider.setValue(int(round((v - lo) / step)))
            slider.blockSignals(False)

        slider.valueChanged.connect(s2b)
        spin.valueChanged.connect(b2s)
        return slider, spin

    @staticmethod
    def _slider_row(slider: QSlider, spin: QDoubleSpinBox) -> QWidget:
        w = QWidget()
        h = QHBoxLayout(w)
        h.setContentsMargins(0, 0, 0, 0)
        h.addWidget(slider, 1)
        h.addWidget(spin)
        return w

    def load_values(self, zs: float, cm: float, cluster: Dict[str, object]):
        """Show ZS / CM thresholds and ClusterConfig field values (by field
        name) without emitting ``changed``.  Only the dock's own signals
        are blocked, so each slider still follows its spin box."""
        blocked = self.blockSignals(True)
        try:
            self.zs_spin.setValue(float(zs))
            self.cm_spin.setValue(float(cm))
            for name, ed in self.cluster_editors.items():
                set_editor_value(ed, cluster[name])
        finally:
            self.blockSignals(blocked)

    def apply_to_system(self, gsys: "det.GemSystem"):
        """Push live threshold values into a GemSystem."""
        gsys.zero_sup_threshold = float(self.zs_spin.value())
        gsys.common_mode_threshold = float(self.cm_spin.value())

    def apply_to_recon(self, gsys: "det.GemSystem"):
        """Write the clustering / XY-match values into every per-detector
        ClusterConfig of a GemSystem; reconstruct() clusters with those,
        not with the GemCluster's own config."""
        cfgs = list(gsys.get_recon_configs())
        for cfg in cfgs:
            editors_to_config(self.cluster_editors, cfg)
        gsys.set_recon_configs(cfgs)


# ---- Config discovery helpers ----


def default_daq_config() -> Optional[Path]:
    return find_database_file("daq_config.json")


def default_gem_map() -> Optional[Path]:
    return find_database_file("gem_map.json")


def _build_gem_pipeline(daq_config: str, gem_map: str,
                        gem_ped: Optional[str] = None,
                        evio: Optional[str] = None):
    """Detector pipeline wired exactly as the server and replay wire it:
    reconstruction_config strip cuts and per-detector cluster configs,
    runinfo common-mode ranges, and the runinfo pedestal file unless
    *gem_ped* is given.  The runinfo entry follows the run number sniffed
    from *evio* (latest entry when unknown).

    The daq_config directory serves as the database directory.  Paths are
    made absolute first because the builder resolves relative ones against
    that directory, not the CWD.
    """
    daq_config = os.path.abspath(daq_config)
    b = det.PipelineBuilder()
    b.set_database_dir(os.path.dirname(daq_config))
    b.set_daq_config(daq_config)
    b.set_gem_map(os.path.abspath(gem_map))
    if gem_ped:
        b.set_gem_pedestal(os.path.abspath(gem_ped))
    if evio:
        b.set_run_number_from_evio(evio)
    return b.build()


def _gem_geometry(gem_map: str, gsys):
    """(detectors, apv_map, hole, raw) for *gem_map*.  *gsys* is an
    initialized GemSystem that supplies the active extent and hole offset."""
    layers, apvs, hole, raw = load_gem_map(gem_map)
    return (build_strip_layout(layers, apvs, hole, raw, gem_sys=gsys),
            build_apv_map(apvs), hole, raw)


# ---- Main window ----


class GemEventViewer(QMainWindow):
    REDRAW_DEBOUNCE_MS = 120

    def __init__(self,
                 initial_evio: Optional[str] = None,
                 daq_config_path: Optional[str] = None,
                 gem_map_path: Optional[str] = None,
                 gem_ped_path: Optional[str] = None):
        super().__init__()
        self.setWindowTitle("GEM Event Viewer")
        self.resize(1400, 820)
        apply_theme_palette(self)

        self._daq_config_path = str(daq_config_path or default_daq_config() or "")
        self._gem_map_path    = str(gem_map_path or default_gem_map() or "")
        # Pedestal override (--gem-ped / File → Choose GEM pedestal file);
        # empty means the runinfo per-run file.  _ped_in_use is the file
        # actually loaded ("" when none).
        self._gem_ped_path    = str(gem_ped_path or "")
        self._ped_in_use      = ""

        # Latched once the full-readout/no-pedestal warning dialog has been
        # shown, so it is not re-shown on every event.
        self._ped_warning_shown = False

        # Geometry (reloaded with the pipeline)
        self._detectors: Dict[int, dict] = {}
        self._hole: Optional[dict] = None
        self._gem_raw: dict = {}
        self._apv_map: dict = {}

        # GEM reconstruction objects (live); _gsys belongs to _pipeline.
        self._pipeline = None
        self._gsys: Optional[det.GemSystem] = None
        self._gcl: Optional[det.GemCluster] = None

        # Dock values as configured, captured after each pipeline build so
        # "Reset defaults" restores them rather than wholly-untuned ones.
        self._default_zs = 5.0
        self._default_cm = 20.0
        self._default_cluster: Optional[Dict[str, object]] = None

        # EVIO file state
        self._evio_path = ""
        self._events: List[EventMeta] = []
        self._current = -1
        self._stepper: Optional[Stepper] = None
        self._last_ssp = None  # cached SspEventData for threshold re-runs

        # Worker thread for scanning
        self._scan_thread: Optional[QThread] = None
        self._scan_worker: Optional[ScanWorker] = None
        self._progress: Optional[QProgressDialog] = None

        # Coalesces tuning-dock changes into one reconstruction pass.
        self._redraw_timer = QTimer(self)
        self._redraw_timer.setSingleShot(True)
        self._redraw_timer.timeout.connect(self._re_reconstruct_current)

        self._build_ui()
        self._load_geometry_and_gemsys()

        if initial_evio:
            QTimer.singleShot(50, lambda: self._open_evio(initial_evio))

    # ---- UI construction ----

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(4, 4, 4, 4)
        root.setSpacing(4)

        # --- Top row: navigation (file info lives in the status bar) ---
        top = QHBoxLayout()

        self.btn_prev = QPushButton("◀ Prev"); self.btn_prev.setShortcut("Left")
        self.btn_next = QPushButton("Next ▶"); self.btn_next.setShortcut("Right")
        self.btn_prev.clicked.connect(lambda: self._step(-1))
        self.btn_next.clicked.connect(lambda: self._step(+1))
        self.btn_prev.setEnabled(False); self.btn_next.setEnabled(False)
        top.addWidget(self.btn_prev)
        top.addWidget(self.btn_next)

        top.addWidget(QLabel("Goto #"))
        self.goto_spin = QSpinBox()
        self.goto_spin.setRange(0, 0)
        self.goto_spin.setEnabled(False)
        self.goto_spin.editingFinished.connect(self._on_goto)
        top.addWidget(self.goto_spin)

        self.slider = QSlider(Qt.Orientation.Horizontal)
        self.slider.setRange(0, 0)
        self.slider.setEnabled(False)
        self.slider.valueChanged.connect(self._on_slider)
        top.addWidget(self.slider, 2)

        # Save as PNG — also bound as a window-level Ctrl+S shortcut via
        # the QAction below so the keybinding works regardless of focus.
        self.act_save_png = QAction("Save as PNG…", self)
        self.act_save_png.setShortcut(QKeySequence("Ctrl+S"))
        self.act_save_png.triggered.connect(self._save_canvas_png)
        self.addAction(self.act_save_png)
        self.btn_save_png = QPushButton("Save as PNG…")
        self.btn_save_png.clicked.connect(self.act_save_png.trigger)
        top.addWidget(self.btn_save_png)

        self.btn_advanced = QPushButton("Advanced…")
        self.btn_advanced.setCheckable(True)
        self.btn_advanced.toggled.connect(self._toggle_advanced)
        top.addWidget(self.btn_advanced)

        root.addLayout(top)

        # Canvas + Raw APV tabs share the main area.  Event data flows
        # into both after each process_event() call.
        self.tabs = QTabWidget()
        self.raw_apv_tab = RawApvTab()
        self.tabs.addTab(self.raw_apv_tab, "Raw APV")
        self.canvas = GemEventCanvas()
        self.tabs.addTab(self.canvas, "Clustering")
        root.addWidget(self.tabs, 1)

        # --- Status bar: left = per-event info, right = file info ---
        self.setStatusBar(QStatusBar(self))
        self._status = QLabel("")
        self.statusBar().addPermanentWidget(self._status, 1)
        # Red warning badge — only visible when full-readout data is loaded
        # without a pedestal file.
        self._ped_badge = QLabel("")
        self._ped_badge.setStyleSheet(
            f"color:{THEME.DANGER}; font-weight: bold; padding: 0 8px;")
        self._ped_badge.hide()
        self.statusBar().addPermanentWidget(self._ped_badge)
        # File info (name, event count, mode) — right-aligned, dim text.
        self.file_label = QLabel("(no file loaded)")
        self.file_label.setStyleSheet(themed("color:#8b949e; padding: 0 8px;"))
        self.statusBar().addPermanentWidget(self.file_label)
        self._set_status("Open an EVIO file via File → Open EVIO… (Ctrl+O).")

        # --- Advanced dock (hidden by default) ---
        self.adv_dock = AdvancedDock(self)
        self.adv_dock.hide()
        self.addDockWidget(Qt.DockWidgetArea.RightDockWidgetArea, self.adv_dock)
        self.adv_dock.visibilityChanged.connect(
            lambda vis: self.btn_advanced.setChecked(bool(vis)))
        self.adv_dock.changed.connect(self._on_threshold_change)
        self.adv_dock.resetRequested.connect(self._reset_defaults)

        self._build_menu()

    def _build_menu(self):
        mb = self.menuBar()
        m_file = mb.addMenu("&File")

        act_open = QAction("&Open EVIO…", self)
        act_open.setShortcut(QKeySequence("Ctrl+O"))
        act_open.triggered.connect(self._pick_evio)
        m_file.addAction(act_open)

        act_map = QAction("Choose &gem_map.json…", self)
        act_map.triggered.connect(self._pick_gem_map)
        m_file.addAction(act_map)

        act_ped = QAction("Choose gem_&ped file…", self)
        act_ped.triggered.connect(self._pick_gem_ped)
        m_file.addAction(act_ped)

        m_file.addSeparator()
        act_quit = QAction("&Quit", self)
        act_quit.setShortcut(QKeySequence("Ctrl+Q"))
        act_quit.triggered.connect(self.close)
        m_file.addAction(act_quit)

        m_view = mb.addMenu("&View")
        self.act_adv = QAction("Show &Advanced tuning", self, checkable=True)
        self.act_adv.triggered.connect(self._toggle_advanced)
        m_view.addAction(self.act_adv)

    # ---- Geometry + GemSystem init ----

    def _load_geometry_and_gemsys(self, keep_dock: bool = False):
        """(Re)build the pipeline and geometry.  With *keep_dock* the
        Advanced dock keeps the user's values (they are re-applied on the
        next draw) instead of being reset to the configured defaults."""
        if not HAVE_PRAD2PY:
            self._fatal_prad2py_missing()
            return
        if not self._gem_map_path or not os.path.isfile(self._gem_map_path):
            QMessageBox.warning(
                self, "gem_map.json not found",
                "Could not locate gem_map.json — use File → Choose gem_map.json to set it.")
            return

        try:
            ped = (self._gem_ped_path
                   if os.path.isfile(self._gem_ped_path) else None)
            self._pipeline = _build_gem_pipeline(
                self._daq_config_path, self._gem_map_path, ped,
                self._evio_path or None)
            self._gsys = self._pipeline.gem
            self._gcl = det.GemCluster()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "GemSystem init failed", str(exc))
            self._pipeline = self._gsys = self._gcl = None
            return
        self._ped_in_use = self._pipeline.gem_pedestal_path

        try:
            (self._detectors, self._apv_map, self._hole,
             self._gem_raw) = _gem_geometry(self._gem_map_path, self._gsys)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Bad gem_map.json", str(exc))
            return

        # Keep the configured values for "Reset defaults" and show them:
        # thresholds from gem_map.json, clustering from
        # reconstruction_config.json.  Field values are copied because
        # get_recon_configs() returns views that set_recon_configs()
        # invalidates.
        had_dock = self._default_cluster is not None
        self._default_zs = float(self._gsys.zero_sup_threshold)
        self._default_cm = float(self._gsys.common_mode_threshold)
        cfg = self._gsys.get_recon_configs()[0]
        self._default_cluster = {n: getattr(cfg, n)
                                 for n in self.adv_dock.cluster_editors}
        if not (keep_dock and had_dock):
            self.adv_dock.load_values(self._default_zs, self._default_cm,
                                      self._default_cluster)

        have_ped = bool(self._ped_in_use and os.path.isfile(self._ped_in_use))
        ped_status = ("ped: " + os.path.basename(self._ped_in_use) if have_ped
                      else "no pedestals loaded (required for full-readout data)")
        self._set_status(
            f"GEM system ready — {self._gsys.get_n_detectors()} detectors, "
            f"{self._gsys.get_n_apvs()} APVs, {ped_status}.")

        # Ped-loaded state changed — reset the "already warned" latch so the
        # next full-readout event can warn again if still missing peds.
        self._ped_warning_shown = have_ped

        # Feed the Raw APV tab its static per-APV metadata.  Skip APVs
        # without a DAQ assignment (crate/mpd/adc all -1) — those are
        # placeholder slots in gem_map.json and carry no data.
        try:
            det_names = {d.id: d.name for d in self._gsys.get_detectors()}
        except Exception:
            det_names = {}
        apv_meta = []
        for i in range(self._gsys.get_n_apvs()):
            cfg = self._gsys.get_apv_config(i)
            if int(cfg.crate_id) < 0 or int(cfg.mpd_id) < 0 or int(cfg.adc_ch) < 0:
                continue
            apv_meta.append({
                "apv_index":  i,
                "crate_id":   int(cfg.crate_id),
                "mpd_id":     int(cfg.mpd_id),
                "adc_ch":     int(cfg.adc_ch),
                "det_id":     int(cfg.det_id),
                "plane_type": str(cfg.plane_type),
                "det_pos":    int(cfg.det_pos),
                "det_name":   det_names.get(int(cfg.det_id),
                                            f"det{cfg.det_id}"),
            })
        self.raw_apv_tab.set_apv_metadata(apv_meta)

    # ---- File picker / pre-scan ----

    def _pick_evio(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open EVIO file", "",
            "EVIO files (*.evio *.evio.*);;All files (*)")
        if path:
            self._open_evio(path)

    def _pick_gem_map(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose gem_map.json", self._gem_map_path or "",
            "JSON (*.json);;All files (*)")
        if path:
            self._gem_map_path = path
            self._load_geometry_and_gemsys()
            if self._events and self._current >= 0:
                self._show_event(self._current)

    def _pick_gem_ped(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Choose GEM pedestal file", self._gem_ped_path or "",
            "Pedestal text (*.txt *.dat);;All files (*)")
        if path:
            self._gem_ped_path = path
            self._load_geometry_and_gemsys()
            if self._events and self._current >= 0:
                self._show_event(self._current)

    def _open_evio(self, path: str):
        if not HAVE_PRAD2PY:
            self._fatal_prad2py_missing()
            return
        if self._gsys is None:
            QMessageBox.warning(self, "GEM system not ready",
                                "Load a gem_map.json before opening an EVIO file.")
            return
        if not os.path.isfile(path):
            QMessageBox.warning(self, "File not found", path)
            return

        self._close_current_run()

        self._evio_path = path
        # Rebuild so the runinfo pedestal / common-mode files follow the run.
        self._load_geometry_and_gemsys(keep_dock=True)
        if self._gsys is None:
            return
        self.file_label.setText(f"Scanning: {os.path.basename(path)}")
        self._set_status("Pre-scanning EVIO file for event index…")

        self._progress = QProgressDialog(
            f"Scanning {os.path.basename(path)}…", "Cancel", 0, 0, self)
        self._progress.setWindowTitle("Building event index")
        self._progress.setWindowModality(Qt.WindowModality.ApplicationModal)
        self._progress.setMinimumDuration(0)
        self._progress.setAutoClose(False)
        self._progress.setAutoReset(False)
        self._progress.canceled.connect(self._cancel_scan)
        self._progress.show()

        self._scan_thread = QThread(self)
        self._scan_worker = ScanWorker(path, self._daq_config_path)
        self._scan_worker.moveToThread(self._scan_thread)
        self._scan_thread.started.connect(self._scan_worker.run)
        self._scan_worker.progress.connect(self._on_scan_progress)
        self._scan_worker.finished.connect(self._on_scan_done)
        self._scan_worker.failed.connect(self._on_scan_failed)
        self._scan_thread.start()

    def _cancel_scan(self):
        if self._scan_worker is not None:
            self._scan_worker.request_cancel()

    def _on_scan_progress(self, phys_seen: int, records_seen: int):
        if self._progress is not None:
            self._progress.setLabelText(
                f"Scanning {os.path.basename(self._evio_path)}…\n"
                f"{phys_seen:,} physics events found in {records_seen:,} records")

    def _on_scan_failed(self, msg: str):
        if self._progress is not None:
            self._progress.close(); self._progress = None
        self._tear_down_worker()
        QMessageBox.critical(self, "Scan failed", msg)
        self.file_label.setText("(no file loaded)")
        self._set_status("Scan failed — see dialog.")

    def _on_scan_done(self, events: List[EventMeta], elapsed: float):
        if self._progress is not None:
            self._progress.close(); self._progress = None
        self._tear_down_worker()

        self._events = events
        if not events:
            self.file_label.setText(
                f"{os.path.basename(self._evio_path)} — no physics events")
            self._set_status("No physics events in file.")
            return

        self._stepper = Stepper(self._evio_path, self._daq_config_path)
        try:
            self._stepper.open()
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Reopen failed", str(exc))
            return

        # Configure navigation widgets
        last = len(events) - 1
        ev_nums = [e.event_number for e in events]
        lo_ev, hi_ev = min(ev_nums), max(ev_nums)
        self.slider.blockSignals(True); self.slider.setRange(0, last); self.slider.setValue(0); self.slider.setEnabled(True); self.slider.blockSignals(False)
        self.goto_spin.blockSignals(True); self.goto_spin.setRange(lo_ev, hi_ev); self.goto_spin.setValue(events[0].event_number); self.goto_spin.setEnabled(True); self.goto_spin.blockSignals(False)
        self.btn_prev.setEnabled(True); self.btn_next.setEnabled(True)

        mode_note = ("" if self._stepper.is_random_access()
                     else "  [sequential mode — Prev is slow]")
        self.file_label.setText(
            f"{os.path.basename(self._evio_path)} — {len(events):,} physics events "
            f"(scanned in {elapsed:.1f} s){mode_note}")
        self._show_event(0)

    def _tear_down_worker(self):
        if self._scan_thread is not None:
            self._scan_thread.quit()
            self._scan_thread.wait(2000)
            self._scan_thread = None
        self._scan_worker = None

    def _close_current_run(self):
        if self._stepper is not None:
            self._stepper.close()
            self._stepper = None
        self._events = []
        self._current = -1
        self._last_ssp = None
        self.slider.setEnabled(False); self.goto_spin.setEnabled(False)
        self.btn_prev.setEnabled(False); self.btn_next.setEnabled(False)

    # ---- Event navigation ----

    def _step(self, delta: int):
        if not self._events:
            return
        tgt = max(0, min(len(self._events) - 1, self._current + delta))
        if tgt != self._current:
            self._show_event(tgt)

    def _on_goto(self):
        """User typed an event number in the Goto spinbox.  Find the closest
        event (numbers may not be contiguous) and jump to it."""
        if not self._events:
            return
        want = int(self.goto_spin.value())
        # linear scan, stops once past ``want`` (events sorted by event_number)
        best = 0; best_d = abs(self._events[0].event_number - want)
        for i, e in enumerate(self._events):
            d = abs(e.event_number - want)
            if d < best_d:
                best_d = d; best = i
            if e.event_number >= want:
                break
        if best != self._current:
            self._show_event(best)

    def _on_slider(self, value: int):
        if not self._events:
            return
        if value != self._current:
            self._show_event(value)

    def _show_event(self, event_idx: int):
        if self._stepper is None or not self._events:
            return
        if not (0 <= event_idx < len(self._events)):
            return
        evmeta = self._events[event_idx]
        self._set_status(f"Fetching event #{evmeta.event_number}…")
        QApplication.processEvents()
        try:
            ssp = self._stepper.get_ssp(evmeta, event_idx)
        except Exception as exc:  # noqa: BLE001
            QMessageBox.critical(self, "Read failed", str(exc))
            return
        self._current = event_idx
        self._last_ssp = ssp

        # Reflect in nav widgets without looping back.
        self.slider.blockSignals(True); self.slider.setValue(event_idx); self.slider.blockSignals(False)
        self.goto_spin.blockSignals(True); self.goto_spin.setValue(evmeta.event_number); self.goto_spin.blockSignals(False)

        self._check_pedestal_requirement(ssp)
        self._re_reconstruct_current()

    def _start_auto_pedestals(self) -> None:
        """Generate peds from the currently-loaded EVIO file and apply them."""
        if not self._evio_path:
            return
        # Temp pedestal file in the system temp dir; kept alive until the
        # window closes.
        fd, tmp_path = tempfile.mkstemp(prefix="gem_ped_", suffix=".txt")
        os.close(fd)
        self._auto_ped_tmp = tmp_path

        dlg = QProgressDialog("Accumulating pedestals…", "Cancel", 0, 1000, self)
        dlg.setWindowTitle("Auto-generating pedestals")
        dlg.setWindowModality(Qt.WindowModality.WindowModal)
        dlg.setMinimumDuration(0)
        dlg.setAutoClose(True)
        dlg.setValue(0)
        dlg.show()
        QApplication.processEvents()

        worker = PedestalWorker(self._evio_path, self._daq_config_path,
                                tmp_path, max_events=1000)

        def _on_progress(done: int, target: int):
            dlg.setMaximum(target)
            dlg.setValue(done)
            dlg.setLabelText(
                f"Accumulating pedestals…\n{done:,} / {target:,} events")

        def _on_finished(out_path: str, napvs: int, n_used: int):
            dlg.close()
            try:
                self._gsys.load_pedestals(out_path,
                                          self._pipeline.gem_crate_remap)
            except Exception as exc:  # noqa: BLE001
                QMessageBox.critical(
                    self, "Pedestal load failed",
                    f"Generated {out_path} but load_pedestals raised:\n"
                    f"{type(exc).__name__}: {exc}")
                return
            self._gem_ped_path = out_path
            self._ped_in_use = out_path
            self._ped_badge.hide()
            self._set_status(
                f"Auto-pedestals applied: {napvs} APVs from {n_used:,} "
                f"events → {out_path}")
            self._re_reconstruct_current()

        def _on_failed(msg: str):
            dlg.close()
            QMessageBox.critical(self, "Pedestal generation failed", msg)

        worker.progress.connect(_on_progress)
        self._auto_ped_worker = worker
        self._auto_ped_thread = start_worker_thread(
            self, worker, _on_finished, _on_failed, dialog=dlg)

    def _check_pedestal_requirement(self, ssp) -> None:
        """If this event is full-readout (firmware did not run online ZS)
        and we have no pedestal file, warn the user — loudly once, and via
        a persistent status-bar indicator afterwards.

        Full-readout data without pedestals cannot produce meaningful hits:
        the default noise is 5000 ADC, nothing clears the ZS threshold, and
        every event looks empty.  Better to say so explicitly than leave
        the user staring at a silent canvas.
        """
        have_ped = bool(self._ped_in_use and os.path.isfile(self._ped_in_use))
        if have_ped or not ssp.has_full_readout():
            # Pedestals loaded, or online-ZS data that needs none.
            self._ped_badge.hide()
            return

        self._ped_badge.setText("⚠ NO PEDESTAL FILE — full-readout data will reconstruct empty")
        self._ped_badge.show()
        if not self._ped_warning_shown:
            self._ped_warning_shown = True
            reply = QMessageBox.question(
                self, "Pedestal file required",
                "This file contains <b>full-readout</b> GEM data "
                "(no online zero-suppression), but no pedestal file is "
                "loaded.<br><br>"
                "Without pedestals, zero-suppression uses a default noise "
                "value → every event will look empty.<br><br>"
                "Auto-generate pedestals from this file now?<br>"
                "(reads up to 1000 full-readout events, takes a few seconds)",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.Yes)
            if reply == QMessageBox.StandardButton.Yes:
                self._start_auto_pedestals()

    # ---- Reconstruction + drawing ----

    def _on_threshold_change(self):
        # Debounce: coalesce multiple rapid changes (e.g. slider drag) into
        # one reconstruction pass.
        self._redraw_timer.start(self.REDRAW_DEBOUNCE_MS)

    def _reset_defaults(self):
        if self._gsys is None or self._default_cluster is None:
            return
        self.adv_dock.load_values(self._default_zs, self._default_cm,
                                  self._default_cluster)
        self._re_reconstruct_current()

    def _refill_raw_apv_cache(self):
        """Pull per-APV processed + raw + hit-mask frames for the current
        event and hand them to the Raw APV tab.  Runs after ProcessEvent so
        ``get_apv_frame`` / ``get_apv_hit_mask`` are valid; raw data comes
        straight from the cached SspEventData (what firmware shipped)."""
        if self._gsys is None or self._last_ssp is None:
            return

        # Walk the SSP structure once → (crate, mpd, adc) → raw (128, 6)
        # ndarray.
        MAX_APVS_PER_MPD = 16
        raw_by_addr: Dict[Tuple[int, int, int], np.ndarray] = {}
        cm_by_addr:  Dict[Tuple[int, int, int], np.ndarray] = {}
        full_readout_addrs: set = set()
        ssp = self._last_ssp
        for m in range(ssp.nmpds):
            mpd = ssp.mpd(m)
            if not mpd.present:
                continue
            for a in range(MAX_APVS_PER_MPD):
                apv = mpd.apv(a)
                if not apv.present:
                    continue
                key = (int(mpd.crate_id), int(mpd.mpd_id),
                       int(apv.addr.adc_ch))
                # Copy so the array outlives the SSP object if cached.
                raw_by_addr[key] = np.asarray(apv.strips).copy()
                if apv.has_online_cm:
                    cm_by_addr[key] = np.asarray(apv.online_cm).copy()
                if apv.full_readout:
                    full_readout_addrs.add(key)

        zs_sigma = float(self._gsys.zero_sup_threshold)
        processed:  Dict[int, np.ndarray] = {}
        raw:        Dict[int, np.ndarray] = {}
        hits:       Dict[int, np.ndarray] = {}
        thresholds: Dict[int, np.ndarray] = {}
        cm_traces:  Dict[int, np.ndarray] = {}
        no_hit_fr:  set = set()
        for i in range(self._gsys.get_n_apvs()):
            try:
                processed[i] = self._gsys.get_apv_frame(i)
                hits[i]      = self._gsys.get_apv_hit_mask(i)
            except Exception:
                continue
            # Threshold curve: noise × ZS σ.
            thresholds[i] = self._gsys.get_apv_ped_noise(i) * zs_sigma
            cfg = self._gsys.get_apv_config(i)
            key = (int(cfg.crate_id), int(cfg.mpd_id), int(cfg.adc_ch))
            if key in raw_by_addr:
                raw[i] = raw_by_addr[key]
                if key in cm_by_addr:
                    cm_traces[i] = cm_by_addr[key]
                # Highlight "suspicious" APVs: full-readout (firmware
                # shipped all 128 strips, no online ZS) AND no channel
                # survived software ZS.
                if key in full_readout_addrs and not bool(hits[i].any()):
                    no_hit_fr.add(i)
        self.raw_apv_tab.set_event_data(processed, raw, hits, no_hit_fr,
                                         thresholds=thresholds,
                                         cm_traces=cm_traces)

    def _re_reconstruct_current(self):
        if self._gsys is None or self._gcl is None:
            return
        if self._last_ssp is None:
            return

        # Push dock values into the GemSystem and its cluster configs.
        self.adv_dock.apply_to_system(self._gsys)
        self.adv_dock.apply_to_recon(self._gsys)

        t0 = time.monotonic()
        self._gsys.clear()
        self._gsys.process_event(self._last_ssp)
        self._gsys.reconstruct(self._gcl)
        elapsed_ms = (time.monotonic() - t0) * 1000.0

        # Build the structures gem_view expects.
        det_list = build_det_list_from_gemsys(self._gsys)
        zs_apvs = build_zs_apvs_from_gemsys(self._gsys)
        det_hits = process_zs_hits(zs_apvs, self._apv_map,
                                   self._detectors, self._hole, self._gem_raw)

        evmeta = self._events[self._current]
        title = (f"GEM Event Viewer — "
                 f"ev #{evmeta.event_number}  "
                 f"trig #{evmeta.trigger_number}  "
                 f"bits 0x{evmeta.trigger_bits:X}")
        self.canvas.set_event(self._detectors, det_list, det_hits,
                              self._hole, title=title)

        # Feed the Raw APV tab — bulk numpy bindings keep this cheap.
        self._refill_raw_apv_cache()

        n_2d = sum(len(d.get("hits_2d", [])) for d in det_list)
        self._set_status(
            f"#{self._current + 1}/{len(self._events)}  "
            f"ev={evmeta.event_number}  trig={evmeta.trigger_number}  "
            f"bits=0x{evmeta.trigger_bits:X}  "
            f"2D hits: {n_2d}   reco: {elapsed_ms:.1f} ms")

    # ---- Misc ----

    def _toggle_advanced(self, checked: bool):
        self.adv_dock.setVisible(checked)
        self.act_adv.setChecked(checked)
        self.btn_advanced.setChecked(checked)

    def _save_canvas_png(self):
        if self.canvas is None or self._current < 0:
            QMessageBox.information(self, "Nothing to save",
                                    "Load an EVIO file and step to an event first.")
            return
        default = f"gem_event_{self._events[self._current].event_number}.png" \
                  if self._events else "gem_event.png"
        path, _ = QFileDialog.getSaveFileName(
            self, "Save canvas as PNG", default,
            "PNG image (*.png);;All files (*)")
        if not path:
            return
        if not self.canvas.save_png(path):
            QMessageBox.warning(self, "Save failed",
                                f"Could not write {path}")
        else:
            self._set_status(f"Saved canvas to {path}")

    def _set_status(self, text: str):
        self._status.setText(text)

    def _fatal_prad2py_missing(self):
        QMessageBox.critical(
            self, "prad2py not available",
            "The prad2py pybind11 module could not be imported.\n\n"
            "Build it with:\n"
            "    cmake -DBUILD_PYTHON=ON -S . -B build && cmake --build build\n\n"
            f"Details:\n{PRAD2PY_ERROR}")

    def closeEvent(self, ev):
        # Cancel any in-flight pre-scan before we tear down the run.
        if self._scan_worker is not None:
            self._scan_worker.request_cancel()
            self._tear_down_worker()
        # Clean up the auto-generated pedestal file if we made one.
        tmp = getattr(self, "_auto_ped_tmp", None)
        if tmp:
            try:
                os.remove(tmp)
            except OSError:
                pass
        self._close_current_run()
        super().closeEvent(ev)


# ---- Batch mode — render events / layout to PNG without a GUI ----


def _parse_event_spec(spec: str) -> List[int]:
    """Parse 'N' / 'N-M' / 'N,M-K,...' into a sorted unique list of indices."""
    out: List[int] = []
    for chunk in spec.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if "-" in chunk:
            a, b = chunk.split("-", 1)
            lo, hi = int(a), int(b)
            if lo > hi:
                lo, hi = hi, lo
            out.extend(range(lo, hi + 1))
        else:
            out.append(int(chunk))
    return sorted(set(out))


def _load_json_event(path: str) -> dict:
    import json as _json
    raw = open(path, "rb").read()
    if raw[:2] in (b"\xff\xfe", b"\xfe\xff"):
        text = raw.decode("utf-16")
    elif raw[:3] == b"\xef\xbb\xbf":
        text = raw.decode("utf-8-sig")
    else:
        text = raw.decode("utf-8")
    return _json.loads(text)


def _print_event_summary(det_list, det_hits):
    """Per-detector hit / cluster table printed in --verbose mode."""
    for dd in det_list:
        did = dd["id"]
        hits = det_hits.get(did, {"x": [], "y": []})
        xcl = dd.get("x_clusters", [])
        ycl = dd.get("y_clusters", [])
        print(f"\n  {dd['name']}: {len(hits['x'])} X hits, "
              f"{len(hits['y'])} Y hits, "
              f"{len(xcl)}+{len(ycl)} clusters, "
              f"{len(dd.get('hits_2d', []))} 2D hits")
        if xcl or ycl:
            print(f"  {'plane':>5} {'pos(mm)':>8} {'peak':>8} {'total':>8} "
                  f"{'size':>4} {'tbin':>4} {'xtalk':>5}  strips")
            for plane, cls in [("X", xcl), ("Y", ycl)]:
                for cl in cls:
                    strips = cl.get("hit_strips", [])
                    srange = f"{min(strips)}-{max(strips)}" if strips else ""
                    print(f"  {plane:>5} {cl['position']:>8.2f} "
                          f"{cl['peak_charge']:>8.1f} "
                          f"{cl['total_charge']:>8.1f} {cl['size']:>4} "
                          f"{cl['max_timebin']:>4} "
                          f"{'y' if cl.get('cross_talk') else '':>5}  {srange}")


def _output_path(name: str, out: Optional[str], n_items: int,
                 fallback: str) -> str:
    """PNG path for one of ``n_items`` rendered items: ``fallback`` without
    ``-o``, ``out`` itself for a single item unless it ends in a path
    separator, otherwise ``name`` inside the directory ``out``."""
    if not out:
        return fallback
    if n_items == 1 and not out.endswith(("/", "\\")):
        return out
    os.makedirs(out, exist_ok=True)
    return os.path.join(out, name)


def _resolve_path(user_value, finder):
    if user_value:
        return user_value
    try:
        p = finder()
        return str(p) if p else None
    except Exception:
        return None


def _run_batch_layout(args) -> int:
    gem_map = _resolve_path(args.gem_map, default_gem_map)
    if not gem_map or not os.path.isfile(gem_map):
        print("error: gem_map.json not found (pass -G <path>)", file=sys.stderr)
        return 2
    # Geometry only: a bare GemSystem supplies x_size and the hole offset.
    gsys = det.GemSystem()
    gsys.init(gem_map)
    detectors, _, hole, _ = _gem_geometry(gem_map, gsys)
    det_layout = detectors[min(detectors.keys())]

    out = args.output or "gem_layout.png"
    if not _render_png(out, args.width, args.height, draw_layout,
                       det_layout, hole, show_every=args.show_every,
                       title=f"PRad-II GEM Strip Layout ({det_layout['name']})"):
        print(f"error: failed to save {out}", file=sys.stderr)
        return 1
    print(f"wrote {out}")
    return 0


def _run_batch_evio(args) -> int:
    if not args.evio or not os.path.isfile(args.evio):
        print("error: EVIO file required (first positional arg)", file=sys.stderr)
        return 2
    daq_cfg = _resolve_path(args.daq_config, default_daq_config)
    gem_map = _resolve_path(args.gem_map, default_gem_map)
    if not daq_cfg or not os.path.isfile(daq_cfg):
        print("error: daq_config.json not found (pass -D)", file=sys.stderr); return 2
    if not gem_map or not os.path.isfile(gem_map):
        print("error: gem_map.json not found (pass -G)", file=sys.stderr); return 2
    if args.gem_ped and not os.path.isfile(args.gem_ped):
        print(f"error: pedestal file not found: {args.gem_ped}",
              file=sys.stderr)
        return 2

    indices: List[int] = []
    if args.event is not None:
        indices = [int(args.event)]
    elif args.events:
        indices = _parse_event_spec(args.events)
    if not indices:
        print("error: provide --event N or --events SPEC", file=sys.stderr); return 2

    pipe = _build_gem_pipeline(daq_cfg, gem_map, args.gem_ped, args.evio)
    gsys = pipe.gem
    gcl = det.GemCluster()
    detectors, apv_map, hole, raw = _gem_geometry(gem_map, gsys)

    print(f"Scanning {args.evio} …")
    events, _ = _build_event_index(args.evio, daq_cfg)
    if not events:
        print("error: no physics events found", file=sys.stderr); return 1
    print(f"  {len(events):,} physics events")

    stepper = Stepper(args.evio, daq_cfg)
    stepper.open()
    try:
        rendered = 0
        for idx in indices:
            if not (0 <= idx < len(events)):
                print(f"  skipping index {idx}: out of range", file=sys.stderr)
                continue
            evmeta = events[idx]
            try:
                ssp = stepper.get_ssp(evmeta, idx)
            except Exception as exc:  # noqa: BLE001
                print(f"  skipping index {idx}: {exc}", file=sys.stderr)
                continue
            gsys.clear()
            gsys.process_event(ssp)
            gsys.reconstruct(gcl)
            det_list = build_det_list_from_gemsys(gsys)
            zs_apvs = build_zs_apvs_from_gemsys(gsys)
            det_hits = process_zs_hits(zs_apvs, apv_map, detectors, hole, raw)

            if args.det >= 0:
                det_list = [d for d in det_list if d["id"] == args.det]

            title = (f"GEM Event #{evmeta.event_number}  "
                     f"trig #{evmeta.trigger_number}  "
                     f"bits 0x{evmeta.trigger_bits:X}")
            name = f"gem_event_{evmeta.event_number:06d}.png"
            out_path = _output_path(name, args.output, len(indices), name)
            ok = _render_png(out_path, args.width, args.height,
                             draw_event_panels, detectors, det_list,
                             det_hits, hole, title=title)
            if ok:
                print(f"  wrote {out_path}")
                rendered += 1
                if args.verbose:
                    _print_event_summary(det_list, det_hits)
            else:
                print(f"  error: failed to save {out_path}", file=sys.stderr)
        print(f"Done: {rendered}/{len(indices)} rendered.")
        return 0 if rendered > 0 else 1
    finally:
        stepper.close()


def _run_batch_json(args) -> int:
    import glob as globmod
    gem_map = _resolve_path(args.gem_map, default_gem_map)
    if not gem_map or not os.path.isfile(gem_map):
        print("error: gem_map.json not found (pass -G)", file=sys.stderr); return 2

    files: List[str] = []
    for arg in args.json:
        if os.path.isdir(arg):
            files += sorted(globmod.glob(os.path.join(arg, "gem_event*.json")))
        elif "*" in arg or "?" in arg:
            files += sorted(globmod.glob(arg))
        else:
            files.append(arg)
    files = [f for f in files if f.lower().endswith(".json")]
    if not files:
        print("error: no JSON files found", file=sys.stderr); return 2

    gsys = det.GemSystem()
    gsys.init(gem_map)
    detectors, apv_map, hole, raw = _gem_geometry(gem_map, gsys)

    rendered = 0
    for i, fpath in enumerate(files):
        try:
            event = _load_json_event(fpath)
        except Exception as exc:  # noqa: BLE001
            print(f"  [{i+1}/{len(files)}] {os.path.basename(fpath)} — "
                  f"parse error: {exc}", file=sys.stderr)
            continue
        if not isinstance(event, dict) or "detectors" not in event:
            print(f"  [{i+1}/{len(files)}] {os.path.basename(fpath)} — "
                  f"skipped (not an event file)")
            continue
        det_list = event.get("detectors", [])
        if args.det >= 0:
            det_list = [d for d in det_list if d["id"] == args.det]
        det_hits = process_zs_hits(event.get("zs_apvs", []), apv_map,
                                   detectors, hole, raw)
        ev_num = int(event.get("event_number", i))
        title = f"GEM Cluster View — Event #{ev_num}"

        stem = os.path.splitext(fpath)[0]
        out_path = _output_path(os.path.basename(stem) + ".png", args.output,
                                len(files), stem + ".png")
        ok = _render_png(out_path, args.width, args.height,
                         draw_event_panels, detectors, det_list, det_hits,
                         hole, title=title)
        print(f"  [{i+1}/{len(files)}] {os.path.basename(fpath)} -> "
              f"{out_path}" + (" (failed)" if not ok else ""))
        if ok:
            rendered += 1
            if args.verbose:
                _print_event_summary(det_list, det_hits)
    print(f"Done: {rendered}/{len(files)} rendered.")
    return 0 if rendered > 0 else 1


# ---- CLI entry point ----


def main():
    parser = argparse.ArgumentParser(
        description="PyQt6 GEM event viewer.  Launches an interactive GUI by "
                    "default; export-mode flags (--layout / --event / --events "
                    "/ --json) render PNGs and exit without showing a window.")
    parser.add_argument("inputs", nargs="*", metavar="FILE",
                        help="EVIO file to open on start, or with --json the "
                             "JSON files / directories / globs to render.")
    parser.add_argument("-D", "--daq-config", default=None,
                        help="Override daq_config.json path.")
    parser.add_argument("-G", "--gem-map", default=None,
                        help="Override gem_map.json path.")
    parser.add_argument("-P", "--gem-ped", default=None,
                        help="Pedestal file override (default: runinfo "
                             "per-run file; matters only for full-readout "
                             "data).")
    parser.add_argument("--theme", choices=available_themes(), default="dark",
                        help="Colour theme (GUI only, default: dark).")

    exp = parser.add_argument_group("export mode (no GUI)")
    exp.add_argument("--layout", action="store_true",
                     help="Render strip-layout PNG and exit.")
    exp.add_argument("--event", type=int, default=None,
                     help="Render a single event (index into the EVIO file).")
    exp.add_argument("--events", default=None,
                     help="Event spec for multi-event export: 'N', 'N-M', or "
                          "comma-separated mix (e.g. '10-20,30,45-50').")
    exp.add_argument("--json", nargs="*", default=None,
                     help="Render from gem_dump JSON files / directory / glob "
                          "(here or as FILE args) instead of EVIO.")
    exp.add_argument("-o", "--output", default=None,
                     help="Output PNG (single) or directory (multi).")
    exp.add_argument("--det", type=int, default=-1,
                     help="Export only detector N (default: all).")
    exp.add_argument("--width", type=int, default=None,
                     help="PNG width in pixels (default: 2400, layout: 1500).")
    exp.add_argument("--height", type=int, default=None,
                     help="PNG height in pixels (default: 900, layout: 1100).")
    exp.add_argument("--show-every", type=int, default=8,
                     help="Strip decimation for --layout (default: 8).")
    exp.add_argument("--verbose", action="store_true",
                     help="Print per-event cluster summary to stdout.")
    # Intermixed so the bin/gem_cluster_view alias (`--json "$@"`) accepts
    # flags and file paths in any order.
    args = parser.parse_intermixed_args()
    if args.json is not None:
        args.json += args.inputs
    elif len(args.inputs) > 1:
        parser.error("only one EVIO file can be opened")
    args.evio = args.inputs[0] if args.inputs else None

    batch = args.layout or args.event is not None or \
            args.events is not None or args.json is not None
    if batch:
        # Headless-safe Qt: offscreen platform plugin, no visible windows.
        os.environ.setdefault("QT_QPA_PLATFORM", "offscreen")
        app = QApplication.instance() or QApplication(sys.argv)
        if not HAVE_PRAD2PY:
            print(f"error: prad2py not available\n{PRAD2PY_ERROR}",
                  file=sys.stderr)
            return 2
        if args.layout:
            # Layout defaults: taller PNG since it's a single panel.
            if args.width is None: args.width = 1500
            if args.height is None: args.height = 1100
            return _run_batch_layout(args)
        if args.width is None: args.width = 2400
        if args.height is None: args.height = 900
        if args.json is not None:
            return _run_batch_json(args)
        return _run_batch_evio(args)

    set_theme(args.theme)
    app = QApplication(sys.argv)
    win = GemEventViewer(
        initial_evio=args.evio,
        daq_config_path=args.daq_config,
        gem_map_path=args.gem_map,
        gem_ped_path=args.gem_ped,
    )
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    rc = main()
    if rc is not None:
        sys.exit(rc)
