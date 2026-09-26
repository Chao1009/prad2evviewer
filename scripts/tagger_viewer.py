#!/usr/bin/env python3
"""
Tagger TDC Viewer (PyQt6)
=========================

Interactive viewer for the V1190 TDC banks (0xE107) produced by the tagger
crate (ROC 0x008E).  Two data sources:

  1. An evio file (``*.evio``, ``*.evio.*``) — decoded in-process via the
     ``prad2py`` pybind11 module.  Build with ``-DBUILD_PYTHON=ON``; the
     checkout's ``build/python/`` is found automatically (an installed
     tree relies on ``PYTHONPATH``).
  2. Live ET stream — subscribe to a running ``prad2_server`` WebSocket.

Displays:

  * A bar chart of hits-per-channel for the selected slot (auto-sized
    to 16 / 32 / 64 / 128 based on the highest channel actually hit —
    click a bar to pick that channel).
  * A TDC value histogram for the selected (slot, channel).
  * Event-wise correlation tabs: Δt = A − B, and 2-D tdc(A) vs tdc(B).
  * A tree on the left with per-slot / per-channel hit counts and the
    human-readable counter name loaded from database/tagger_map.json.

Usage
-----

    # Offline (evio file via prad2py)
    python scripts/tagger_viewer.py /data/.../prad_023667.evio.00000

    # Live (online ET via prad2_server)
    python scripts/tagger_viewer.py --live ws://clondaq6:5051

Other options: -n/--max-events N, -D/--daq-config PATH, --roc TAG,
--no-smoke-test, --theme dark|light (see --help).

Only PyQt6 and numpy are required.  Plots are drawn with QPainter, so
matplotlib / pyqtgraph are NOT needed.
"""

from __future__ import annotations

import argparse
import json as _json
import os
import sys
import time as _time
from pathlib import Path
from typing import Dict, Optional, Tuple

import numpy as np

from PyQt6.QtCore import (
    QObject, Qt, QRectF, QTimer, QUrl, pyqtSignal,
)
from PyQt6.QtGui import QAction, QColor, QFont, QImage, QPainter, QPen

from evio_io import open_evio
from hycal_geoview import (
    apply_theme_palette, set_theme, available_themes, THEME,
    start_worker_thread,
)
from prad2_env import PRAD2PY_HINT, find_database_file, import_prad2py
from PyQt6.QtWidgets import (
    QApplication,
    QComboBox,
    QDialog,
    QDialogButtonBox,
    QFileDialog,
    QFrame,
    QHBoxLayout,
    QInputDialog,
    QLabel,
    QMainWindow,
    QMessageBox,
    QPlainTextEdit,
    QProgressDialog,
    QPushButton,
    QSpinBox,
    QSplitter,
    QStatusBar,
    QTabWidget,
    QTreeWidget,
    QTreeWidgetItem,
    QVBoxLayout,
    QWidget,
)
from PyQt6.QtWebSockets import QWebSocket


# Hit dtypes (shared by the live stream and the evio loader)

# 16-byte packed record, matches both the prad2py numpy output and the
# per-hit payload carried by prad2_server's TDC WebSocket frames.
RAW_DTYPE = np.dtype(
    [
        ("event_num", "<u4"),
        ("trigger_bits", "<u4"),
        ("roc_tag", "<u2"),
        ("slot", "u1"),
        ("channel_edge", "u1"),  # bit 7 = edge, bits 6:0 = channel
        ("tdc", "<u4"),
    ]
)
assert RAW_DTYPE.itemsize == 16, "raw record must be 16 bytes"

# Internal representation used everywhere downstream (splits channel/edge
# so the rest of the code doesn't have to mask bits).
RECORD_DTYPE = np.dtype(
    [
        ("event_num", "<u4"),
        ("trigger_bits", "<u4"),
        ("roc_tag", "<u2"),
        ("slot", "u1"),
        ("channel", "u1"),
        ("edge", "u1"),
        ("tdc", "<u4"),
    ]
)


# Live stream frame parser (prad2_server tagger broadcast)

# Header is little-endian:
#   char magic[4] ("TGR1")
#   u32 flags, u32 n_hits, u32 first_seq, u32 last_seq, u32 dropped
STREAM_HEADER_DTYPE = np.dtype(
    [
        ("magic",     "S4"),
        ("flags",     "<u4"),
        ("n_hits",    "<u4"),
        ("first_seq", "<u4"),
        ("last_seq",  "<u4"),
        ("dropped",   "<u4"),
    ]
)
assert STREAM_HEADER_DTYPE.itemsize == 24, "tagger stream header must be 24 bytes"
STREAM_MAGIC = b"TGR1"


def parse_stream_frame(buf) -> Tuple[dict, np.ndarray]:
    """Parse one server-broadcast tagger frame. Returns (header_dict, raw_hits)
    where raw_hits has RAW_DTYPE (packed 16-byte records) ready to be
    translated into RECORD_DTYPE via raw_to_record()."""
    mv = memoryview(buf).tobytes() if not isinstance(buf, (bytes, bytearray)) else bytes(buf)
    if len(mv) < STREAM_HEADER_DTYPE.itemsize:
        raise ValueError(f"frame too short ({len(mv)} bytes)")
    hdr = np.frombuffer(mv, dtype=STREAM_HEADER_DTYPE, count=1)[0]
    if bytes(hdr["magic"]) != STREAM_MAGIC:
        raise ValueError(f"bad magic: {bytes(hdr['magic'])!r}")
    n = int(hdr["n_hits"])
    expected = STREAM_HEADER_DTYPE.itemsize + n * RAW_DTYPE.itemsize
    if len(mv) < expected:
        raise ValueError(f"frame truncated: {len(mv)} < {expected}")
    raw = np.frombuffer(mv, dtype=RAW_DTYPE, count=n,
                        offset=STREAM_HEADER_DTYPE.itemsize)
    return (
        {
            "flags":     int(hdr["flags"]),
            "n_hits":    n,
            "first_seq": int(hdr["first_seq"]),
            "last_seq":  int(hdr["last_seq"]),
            "dropped":   int(hdr["dropped"]),
        },
        raw,
    )


def raw_to_record(raw: np.ndarray) -> np.ndarray:
    """Unpack RAW_DTYPE rows into RECORD_DTYPE rows (splits channel_edge)."""
    hits = np.empty(raw.size, dtype=RECORD_DTYPE)
    hits["event_num"]    = raw["event_num"]
    hits["trigger_bits"] = raw["trigger_bits"]
    hits["roc_tag"]      = raw["roc_tag"]
    hits["slot"]         = raw["slot"]
    hits["channel"]      = raw["channel_edge"] & 0x7F
    hits["edge"]         = (raw["channel_edge"] >> 7) & 0x1
    hits["tdc"]          = raw["tdc"]
    return hits


# In-process evio loader via prad2py

prad2py, PRAD2PY_ERROR = import_prad2py()
HAVE_PRAD2PY = prad2py is not None


def load_hits_from_evio(
    path: str,
    *,
    max_events: int = 0,
    daq_config: str = "",
    roc_filter: int = -1,
    progress=None,           # callable(events_seen, hits_seen) -> bool
    progress_every: int = 10000,
) -> np.ndarray:
    """Event-wise loop that extracts only the TDC bank from an evio file.

    Drives the per-event decode from Python via
    ``EvChannel.select_event()`` + ``info()``/``tdc()``.  No FADC / SSP /
    VTP decoding happens — typically 5-10× faster than the full-event path.

    ``progress(events_seen, hits_seen)`` is invoked every
    ``progress_every`` physics events.  If it returns False, the loop
    stops early and the collected hits are returned.
    """
    if not HAVE_PRAD2PY:
        raise RuntimeError(
            "prad2py module not available "
            f"({PRAD2PY_ERROR or 'not importable'}).\n" + PRAD2PY_HINT
        )
    dec = prad2py.dec
    ch, _ = open_evio(path, daq_config)

    # Accumulate per-event batches as (event_num, trigger_bits, hits_array).
    # We convert the collected lists into one structured numpy array at the
    # end — cheap relative to the decoding itself.
    ev_nums_chunks = []
    trig_chunks    = []
    hits_chunks    = []

    n_physics = 0
    stop_requested = False

    while ch.read() == dec.Status.success:
        if not ch.scan():
            continue
        if ch.get_event_type() != dec.EventType.Physics:
            continue

        for i in range(ch.get_n_events()):
            ch.select_event(i)
            info    = ch.info()
            tdc_evt = ch.tdc()
            n = tdc_evt.n_hits
            if n > 0:
                hits = tdc_evt.hits_numpy   # structured array, length n
                if roc_filter >= 0:
                    hits = hits[hits["roc_tag"] == roc_filter]
                    n = hits.size
                if n > 0:
                    hits_chunks.append(hits)
                    ev_nums_chunks.append(
                        np.full(n, info.event_number, dtype=np.uint32))
                    trig_chunks.append(
                        np.full(n, info.trigger_bits, dtype=np.uint32))

            n_physics += 1
            if max_events and n_physics >= max_events:
                stop_requested = True
                break
            if progress is not None and (n_physics % progress_every) == 0:
                total_hits = sum(h.size for h in hits_chunks)
                keep = progress(n_physics, total_hits)
                if keep is False:
                    stop_requested = True
                    break

        if stop_requested:
            break

    ch.close()

    # Final "we're done" progress tick.
    if progress is not None:
        total_hits = sum(h.size for h in hits_chunks)
        progress(n_physics, total_hits)

    if not hits_chunks:
        return np.zeros(0, dtype=RECORD_DTYPE)

    ev_nums = np.concatenate(ev_nums_chunks)
    trigs   = np.concatenate(trig_chunks)
    hits    = np.concatenate(hits_chunks)
    out = np.empty(hits.size, dtype=RECORD_DTYPE)
    out["event_num"]    = ev_nums
    out["trigger_bits"] = trigs
    out["roc_tag"]      = hits["roc_tag"]
    out["slot"]         = hits["slot"]
    out["channel"]      = hits["channel"]
    out["edge"]         = hits["edge"]
    out["tdc"]          = hits["value"]
    return out


def round_up_channels(max_ch: int) -> int:
    """Round the observed max channel number up to a V1190-friendly power
    of two (16, 32, 64, 128). Channels are 0-indexed in the data, so we
    compare against (max_ch + 1)."""
    need = max(int(max_ch) + 1, 1)
    for n in (16, 32, 64, 128):
        if need <= n:
            return n
    return 128


# Channel-name map (database/tagger_map.json)

def load_channel_map(path: Optional[Path] = None) -> Dict[Tuple[int, int], str]:
    """Load slot+channel -> name mapping. Silently returns {} on any failure
    (missing file, bad JSON, etc.) — names are purely decorative."""
    p = (Path(path) if path is not None
         else find_database_file("tagger_map.json", use_env=False))
    if p is None or not p.is_file():
        return {}
    try:
        with p.open("r", encoding="utf-8") as f:
            data = _json.load(f)
        out: Dict[Tuple[int, int], str] = {}
        for e in data.get("channels", []):
            out[(int(e["slot"]), int(e["channel"]))] = str(e["name"])
        return out
    except Exception:
        return {}


# Plot widgets


class _TaggerPlot(QWidget):
    """QPainter plot with a title above a plot rectangle inset by
    ``MARGINS`` (left, top, right, bottom) pixels."""

    MARGINS = (0, 0, 0, 0)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._title = ""

    def setTitle(self, title: str):
        self._title = title
        self.update()

    def _plotRect(self) -> QRectF:
        left, top, right, bottom = self.MARGINS
        return QRectF(left, top, self.width() - left - right,
                      self.height() - top - bottom)

    def _begin_paint(self, frame: bool = True) -> Tuple[QPainter, QRectF]:
        """(painter, plot rect) with the canvas cleared and the plot frame
        (when ``frame``) and the title drawn, pen left at THEME.TEXT."""
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        p.fillRect(self.rect(), QColor(THEME.CANVAS))
        r = self._plotRect()
        p.setPen(QPen(QColor(THEME.TEXT)))
        if frame:
            p.drawRect(r)
        if self._title:
            p.setFont(QFont("Monospace", 10, QFont.Weight.Bold))
            p.drawText(int(r.left()), int(r.top() - 6), self._title)
        return p, r


class BarChart(_TaggerPlot):
    """
    Horizontal index → count bar chart painted with QPainter.
    Emits ``barClicked(index)`` when a bar is clicked.
    """

    MARGINS = (50, 18, 10, 30)

    barClicked = pyqtSignal(int)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._counts: np.ndarray = np.zeros(0, dtype=np.int64)
        self._highlight: Optional[int] = None
        self.setMinimumHeight(180)
        self.setMouseTracking(True)

    def setData(self, counts: np.ndarray):
        self._counts = np.asarray(counts, dtype=np.int64)
        self._highlight = None
        self.update()

    def setHighlight(self, idx: Optional[int]):
        self._highlight = idx
        self.update()

    def _indexAtX(self, x: float) -> Optional[int]:
        r = self._plotRect()
        if not r.contains(x, r.center().y()):
            return None
        n = self._counts.size
        if n <= 0:
            return None
        rel = (x - r.left()) / r.width()
        idx = int(rel * n)
        if 0 <= idx < n:
            return idx
        return None

    def mousePressEvent(self, ev):
        if ev.button() == Qt.MouseButton.LeftButton:
            idx = self._indexAtX(ev.position().x())
            if idx is not None:
                self.barClicked.emit(idx)

    def paintEvent(self, _ev):
        p, r = self._begin_paint()

        n = self._counts.size
        if n <= 0:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, "(no hits)")
            return

        cmax = int(self._counts.max()) if self._counts.size else 1
        cmax = max(cmax, 1)
        bar_w = r.width() / n

        for i, c in enumerate(self._counts):
            h = (c / cmax) * r.height()
            x0 = r.left() + i * bar_w
            y0 = r.bottom() - h
            color = QColor(THEME.ACCENT)
            if self._highlight is not None and i == self._highlight:
                color = QColor(THEME.HIGHLIGHT)
            elif c == 0:
                color = QColor(THEME.BUTTON_HOVER)
            p.fillRect(QRectF(x0 + 0.5, y0, max(bar_w - 1.0, 1.0), h), color)

        # y-axis ticks
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 8))
        for frac in (0.0, 0.5, 1.0):
            y = r.bottom() - frac * r.height()
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            p.drawText(
                int(r.left() - 38),
                int(y + 4),
                f"{int(cmax * frac):,}",
            )

        # x-axis ticks
        step = max(1, n // 16)
        for i in range(0, n, step):
            x = r.left() + (i + 0.5) * bar_w
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            p.drawText(int(x - 14), int(r.bottom() + 14), str(i))


class Histogram(_TaggerPlot):
    """1-D histogram painted with QPainter."""

    MARGINS = (65, 20, 10, 40)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._counts: np.ndarray = np.zeros(0, dtype=np.int64)
        self._edges: np.ndarray = np.zeros(0)
        self._xlabel = ""
        self.setMinimumHeight(260)

    def setData(self, counts: np.ndarray, edges: np.ndarray):
        self._counts = np.asarray(counts, dtype=np.int64)
        self._edges = np.asarray(edges, dtype=np.float64)
        self.update()

    def setXLabel(self, label: str):
        self._xlabel = label
        self.update()

    def paintEvent(self, _ev):
        p, r = self._begin_paint()

        n = self._counts.size
        if n <= 0 or self._counts.sum() == 0:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, "(no hits for this channel)")
            return

        cmax = int(self._counts.max())
        cmax = max(cmax, 1)
        bar_w = r.width() / n

        p.setPen(Qt.PenStyle.NoPen)
        for i, c in enumerate(self._counts):
            if c == 0:
                continue
            h = (c / cmax) * r.height()
            x0 = r.left() + i * bar_w
            y0 = r.bottom() - h
            p.fillRect(
                QRectF(x0, y0, max(bar_w, 1.0), h),
                QColor(THEME.ACCENT),
            )

        # y ticks
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 8))
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = r.bottom() - frac * r.height()
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            p.drawText(
                int(r.left() - 46),
                int(y + 4),
                f"{int(cmax * frac):,}",
            )

        # x ticks
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            x = r.left() + frac * r.width()
            val_idx = int(round(frac * (self._edges.size - 1)))
            val = self._edges[val_idx] if self._edges.size > 0 else 0
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            p.drawText(int(x - 30), int(r.bottom() + 14), f"{val:.0f}")

        if self._xlabel:
            p.drawText(
                int(r.center().x() - 60),
                int(r.bottom() + 30),
                self._xlabel,
            )


class Heatmap2D(_TaggerPlot):
    """2-D histogram rendered via a scaled QImage (numpy-built RGB buffer)."""

    MARGINS = (75, 25, 15, 50)

    def __init__(self, parent=None):
        super().__init__(parent)
        self._counts: np.ndarray = np.zeros((0, 0), dtype=np.int64)
        self._xedges: np.ndarray = np.zeros(0)
        self._yedges: np.ndarray = np.zeros(0)
        self._xlabel = "X"
        self._ylabel = "Y"
        self._image: Optional[QImage] = None
        self._rgb_buffer: Optional[np.ndarray] = None  # keep alive for QImage
        self.setMinimumHeight(300)

    def setData(self, counts: np.ndarray, xedges: np.ndarray, yedges: np.ndarray):
        self._counts = np.asarray(counts, dtype=np.int64)
        self._xedges = np.asarray(xedges, dtype=np.float64)
        self._yedges = np.asarray(yedges, dtype=np.float64)
        self._rebuild_image()
        self.update()

    def setLabels(self, xlabel: str, ylabel: str):
        self._xlabel = xlabel
        self._ylabel = ylabel
        self.update()

    def _rebuild_image(self):
        n = self._counts
        if n.size == 0 or n.sum() == 0:
            self._image = None
            self._rgb_buffer = None
            return
        cmax = max(1.0, float(n.max()))
        # np.histogram2d returns shape (nxbins, nybins) with H[i,j] = count in
        # (x bin i, y bin j).  We need an RGB buffer of shape (height, width, 3)
        # where height → y-axis bins (flipped because screen y is downward) and
        # width → x-axis bins.
        nxbins, nybins = n.shape
        t = n.astype(np.float64) / cmax
        # Viridis approximation (same polynomial as geo.js PALETTES.viridis).
        r = np.clip(-0.87 + 4.26*t - 4.85*t*t + 2.5*t*t*t, 0.0, 1.0)
        g = np.clip(-0.03 + 0.77*t + 1.32*t*t - 1.87*t*t*t, 0.0, 1.0)
        b = np.clip( 0.33 + 1.74*t - 4.26*t*t + 3.17*t*t*t, 0.0, 1.0)
        # Zero bins → theme's "no data" surface so they don't pick up the
        # low-t colormap colour. Resolved at paint time so it follows --theme.
        zero_rgb = QColor(THEME.NO_DATA)
        zero = n == 0
        r[zero] = zero_rgb.redF()
        g[zero] = zero_rgb.greenF()
        b[zero] = zero_rgb.blueF()

        # Transpose to (nybins, nxbins) then flip vertically.
        r8 = (r.T[::-1] * 255).astype(np.uint8)
        g8 = (g.T[::-1] * 255).astype(np.uint8)
        b8 = (b.T[::-1] * 255).astype(np.uint8)
        rgb = np.empty((nybins, nxbins, 3), dtype=np.uint8)
        rgb[..., 0] = r8
        rgb[..., 1] = g8
        rgb[..., 2] = b8
        rgb = np.ascontiguousarray(rgb)
        self._rgb_buffer = rgb
        self._image = QImage(
            rgb.data, nxbins, nybins, 3 * nxbins, QImage.Format.Format_RGB888
        )

    def paintEvent(self, _ev):
        p, r = self._begin_paint(frame=False)

        if self._image is None:
            p.drawRect(r)
            p.setPen(QColor(THEME.TEXT_DIM))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter,
                       "(no matched events — set both channels A and B)")
            return

        p.drawImage(r, self._image)
        p.setPen(QPen(QColor(THEME.TEXT)))
        p.drawRect(r)

        # tick labels
        p.setFont(QFont("Monospace", 8))
        p.setPen(QColor(THEME.TEXT_DIM))
        xe, ye = self._xedges, self._yedges
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            x = r.left() + frac * r.width()
            xv = xe[int(round(frac * (xe.size - 1)))] if xe.size > 0 else 0
            p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
            p.drawText(int(x - 26), int(r.bottom() + 14), f"{xv:.0f}")
        for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
            y = r.bottom() - frac * r.height()
            yv = ye[int(round(frac * (ye.size - 1)))] if ye.size > 0 else 0
            p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
            p.drawText(int(r.left() - 48), int(y + 4), f"{yv:.0f}")

        # axis labels
        if self._xlabel:
            p.drawText(int(r.center().x() - 40), int(r.bottom() + 30), self._xlabel)
        if self._ylabel:
            p.drawText(int(r.left() - 50), int(r.top() - 8), self._ylabel)


# Offline file loader (runs in a QThread so the UI stays responsive)


class LoadWorker(QObject):
    """Drives ``load_hits_from_evio`` on a worker thread, emitting
    ``progressed`` every N events; ``request_cancel()`` stops the loop at
    its next progress check."""

    finished   = pyqtSignal(object)         # numpy ndarray
    failed     = pyqtSignal(str)
    progressed = pyqtSignal(int, int)       # (events_seen, hits_seen)

    def __init__(self, path: str, max_events: int,
                 daq_config: str, roc_filter: int):
        super().__init__()
        self._path = path
        self._max = max_events
        self._daq = daq_config
        self._roc = roc_filter
        self._cancel = False

    def request_cancel(self):
        """Flag the loop to stop at its next progress-check point."""
        self._cancel = True

    # --- called from the worker thread ---------------------------------

    def _progress_cb(self, events: int, hits: int) -> bool:
        self.progressed.emit(int(events), int(hits))
        return not self._cancel

    def run(self):
        try:
            hits = load_hits_from_evio(
                self._path,
                max_events=self._max,
                daq_config=self._daq,
                roc_filter=self._roc,
                progress=self._progress_cb,
                progress_every=5000,
            )
        except Exception as exc:                  # noqa: BLE001
            self.failed.emit(f"{type(exc).__name__}: {exc}")
            return
        self.finished.emit(hits)


# Live WebSocket stream (prad2_server)


class LiveStream(QObject):
    """QWebSocket-based live subscriber for prad2_server TDC broadcasts.

    All Qt signals emit on the main GUI thread (QWebSocket runs inside the
    Qt event loop), so callers do not need any additional locking around
    the numpy arrays passed through ``hitsReceived``.
    """

    hitsReceived = pyqtSignal(np.ndarray)   # RECORD_DTYPE rows (one batch)
    stateChanged = pyqtSignal(str)          # free-form label for the status bar
    statsUpdate  = pyqtSignal(dict)         # {rate_hz, dropped, flags}
    subscribed   = pyqtSignal()             # server acknowledged tagger_subscribe
    failed       = pyqtSignal(str)          # raw socket error string

    def __init__(self, parent=None):
        super().__init__(parent)
        self._ws = QWebSocket()
        self._ws.connected.connect(self._on_open)
        self._ws.disconnected.connect(self._on_close)
        self._ws.binaryMessageReceived.connect(self._on_binary)
        self._ws.textMessageReceived.connect(self._on_text)
        try:
            self._ws.errorOccurred.connect(self._on_error)
        except AttributeError:
            # Qt 6.0 spelled this differently; don't crash on older binaries.
            pass

        self._paused = False
        self._reset_stats()

        self._stats_timer = QTimer(self)
        self._stats_timer.setInterval(500)
        self._stats_timer.timeout.connect(self._emit_stats)

    def _reset_stats(self):
        self._total_hits = 0
        self._last_dropped = 0
        self._last_flags = 0
        self._stats_t = _time.monotonic()
        self._stats_hits = 0

    def open(self, url: str):
        self.stateChanged.emit(f"connecting to {url} …")
        self._reset_stats()
        self._ws.open(QUrl(url))
        self._stats_timer.start()

    def is_active(self) -> bool:
        """True between open() and close() or the socket disconnecting."""
        return self._stats_timer.isActive()

    def close(self):
        self._stats_timer.stop()
        if self._ws.state() != self._ws.state().UnconnectedState:
            try:
                self._ws.sendTextMessage(
                    _json.dumps({"type": "tagger_unsubscribe"})
                )
            except Exception:
                pass
        self._ws.close()

    def set_paused(self, paused: bool):
        self._paused = bool(paused)

    def _on_open(self):
        self.stateChanged.emit("connected, subscribing…")
        self._ws.sendTextMessage(_json.dumps({"type": "tagger_subscribe"}))

    def _on_close(self):
        self._stats_timer.stop()
        self.stateChanged.emit("disconnected")

    def _on_text(self, msg: str):
        try:
            d = _json.loads(msg)
        except Exception:
            return
        t = d.get("type")
        if t == "tagger_subscribed":
            n = d.get("subscribers", "?")
            self.stateChanged.emit(f"subscribed ({n} client(s))")
            self.subscribed.emit()

    def _on_binary(self, data):
        if self._paused:
            return
        try:
            hdr, raw = parse_stream_frame(bytes(data))
        except Exception as exc:
            self.stateChanged.emit(f"frame error: {exc}")
            return
        hits = raw_to_record(raw)
        self._total_hits += hits.size
        self._last_dropped = hdr["dropped"]
        self._last_flags = hdr["flags"]
        if hits.size:
            self.hitsReceived.emit(hits)

    def _on_error(self, _err):
        self.stateChanged.emit(f"error: {self._ws.errorString()}")
        self.failed.emit(self._ws.errorString())

    def _emit_stats(self):
        now = _time.monotonic()
        dt = max(now - self._stats_t, 1e-6)
        rate = (self._total_hits - self._stats_hits) / dt
        self._stats_t = now
        self._stats_hits = self._total_hits
        self.statsUpdate.emit({
            "rate_hz": rate,
            "dropped": self._last_dropped,
            "flags":   self._last_flags,
        })


# Event-wise correlation helpers


def _first_hits_for(hits: np.ndarray, slot: int, channel: int,
                    edge_sel: Optional[int]) -> np.ndarray:
    """Return a structured array with one row per event for (slot, channel).

    When a channel fires multiple times within one event (multi-hit TDC),
    we keep the earliest hit (smallest TDC value) — the usual convention for
    timing correlations.
    """
    mask = (hits["slot"] == slot) & (hits["channel"] == channel)
    if edge_sel is not None:
        mask = mask & (hits["edge"] == edge_sel)
    sub = hits[mask]
    if sub.size == 0:
        return sub

    # Sort by (event_num, tdc) → the first occurrence of each event_num is the
    # earliest hit.  np.unique keeps the first index.
    order = np.lexsort((sub["tdc"], sub["event_num"]))
    sub = sub[order]
    _, first_idx = np.unique(sub["event_num"], return_index=True)
    return sub[first_idx]


def _match_pair(hits: np.ndarray,
                a: Tuple[int, int], b: Tuple[int, int],
                edge_sel: Optional[int]) -> Tuple[np.ndarray, np.ndarray]:
    """Inner-join on event_num. Returns (tdc_a, tdc_b) aligned by event."""
    ha = _first_hits_for(hits, a[0], a[1], edge_sel)
    hb = _first_hits_for(hits, b[0], b[1], edge_sel)
    if ha.size == 0 or hb.size == 0:
        return np.empty(0, dtype=np.int64), np.empty(0, dtype=np.int64)
    _, ia, ib = np.intersect1d(
        ha["event_num"], hb["event_num"], return_indices=True, assume_unique=True
    )
    return ha["tdc"][ia].astype(np.int64), hb["tdc"][ib].astype(np.int64)


def _int_range(v: np.ndarray) -> Tuple[int, int]:
    """Integer (min, max) of ``v``, with max raised to min + 1 when equal."""
    lo, hi = int(v.min()), int(v.max())
    return lo, max(hi, lo + 1)


# Resizable error dialog


def show_error_dialog(parent, title: str, heading: str, details: str,
                      width: int = 720, height: int = 420) -> None:
    """Modal error dialog with a read-only, monospace, scrollable text area.

    Unlike ``QMessageBox.critical``, this dialog is resizable by the user so
    long tracebacks / multi-line messages are actually readable.
    """
    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    dlg.resize(width, height)
    dlg.setSizeGripEnabled(True)
    apply_theme_palette(dlg)
    dlg.setStyleSheet(_app_stylesheet())

    lay = QVBoxLayout(dlg)

    if heading:
        lbl = QLabel(heading)
        lbl.setWordWrap(True)
        lbl.setTextInteractionFlags(Qt.TextInteractionFlag.TextSelectableByMouse)
        lay.addWidget(lbl)

    text = QPlainTextEdit()
    text.setReadOnly(True)
    text.setLineWrapMode(QPlainTextEdit.LineWrapMode.NoWrap)
    # Monospace so file paths / tracebacks align.
    mono = QFont("Monospace")
    mono.setStyleHint(QFont.StyleHint.TypeWriter)
    mono.setPointSize(10)
    text.setFont(mono)
    text.setPlainText(details)
    lay.addWidget(text, 1)

    buttons = QDialogButtonBox(QDialogButtonBox.StandardButton.Ok, parent=dlg)
    buttons.accepted.connect(dlg.accept)
    lay.addWidget(buttons)

    dlg.exec()


# Application-wide Qt stylesheet


def _app_stylesheet() -> str:
    """Apple-inspired stylesheet applied to the top-level window.

    Qt stylesheets cascade to child widgets that don't set their own, so
    this is the single place where buttons / inputs / tree / tab chrome
    pick up the active :class:`THEME`. Reads values live, so callers
    should rebuild this string after :func:`set_theme`.
    """
    return (
        # --- Push buttons (8px radius, Apple blue focus) ------------------
        f"QPushButton {{"
        f"  background:{THEME.BUTTON};"
        f"  color:{THEME.TEXT};"
        f"  border:1px solid {THEME.BORDER};"
        f"  border-radius:8px;"
        f"  padding:5px 14px;"
        f"}}"
        f"QPushButton:hover    {{ background:{THEME.BUTTON_HOVER}; }}"
        f"QPushButton:pressed  {{ background:{THEME.BUTTON_HOVER}; }}"
        f"QPushButton:checked  {{"
        f"  background:{THEME.ACCENT_STRONG};"
        f"  color:#ffffff;"
        f"  border:1px solid {THEME.ACCENT_BORDER};"
        f"}}"
        f"QPushButton:disabled {{ color:{THEME.TEXT_MUTED}; }}"
        f"QPushButton:focus    {{ outline:none; border:1px solid {THEME.ACCENT_BORDER}; }}"

        # --- Text inputs / spin / combo (6px radius) ----------------------
        f"QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {{"
        f"  background:{THEME.PANEL};"
        f"  color:{THEME.TEXT};"
        f"  border:1px solid {THEME.BORDER};"
        f"  border-radius:6px;"
        f"  padding:2px 6px;"
        f"  selection-background-color:{THEME.ACCENT};"
        f"  selection-color:#ffffff;"
        f"}}"
        f"QComboBox QAbstractItemView {{"
        f"  background:{THEME.PANEL};"
        f"  color:{THEME.TEXT};"
        f"  border:1px solid {THEME.BORDER};"
        f"  selection-background-color:{THEME.ACCENT};"
        f"}}"

        # --- Check boxes -------------------------------------------------
        f"QCheckBox {{ color:{THEME.TEXT}; spacing:6px; }}"
        f"QCheckBox::indicator {{"
        f"  width:14px; height:14px; border-radius:3px;"
        f"  border:1px solid {THEME.BORDER}; background:{THEME.PANEL};"
        f"}}"
        f"QCheckBox::indicator:hover   {{ border:1px solid {THEME.ACCENT}; }}"
        f"QCheckBox::indicator:checked {{ background:{THEME.ACCENT}; border:1px solid {THEME.ACCENT}; }}"

        # --- Labels (default text colour — frame-shaped labels get a panel)
        f"QLabel {{ color:{THEME.TEXT}; }}"
        f"QLabel[frameShape=\"4\"], QLabel[frameShape=\"1\"] {{"
        f"  background:{THEME.PANEL}; border:1px solid {THEME.BORDER};"
        f"  border-radius:6px; padding:1px 6px;"
        f"}}"

        # --- Tree / Tabs / Splitter / Menu / Status / Tooltip -----------
        f"QTreeWidget, QTreeView {{"
        f"  background:{THEME.PANEL};"
        f"  color:{THEME.TEXT};"
        f"  border:1px solid {THEME.BORDER};"
        f"  alternate-background-color:{THEME.ALT_BASE};"
        f"  selection-background-color:{THEME.ACCENT};"
        f"  selection-color:#ffffff;"
        f"}}"
        f"QHeaderView::section {{"
        f"  background:{THEME.PANEL};"
        f"  color:{THEME.TEXT_DIM};"
        f"  border:0; border-bottom:1px solid {THEME.BORDER};"
        f"  padding:4px 6px;"
        f"}}"
        f"QTabWidget::pane {{ border:1px solid {THEME.BORDER}; border-radius:6px; }}"
        f"QTabBar::tab {{"
        f"  background:transparent;"
        f"  color:{THEME.TEXT_DIM};"
        f"  padding:5px 12px; margin-right:2px;"
        f"  border-top-left-radius:6px; border-top-right-radius:6px;"
        f"}}"
        f"QTabBar::tab:selected {{ color:{THEME.TEXT}; border-bottom:2px solid {THEME.ACCENT}; }}"
        f"QTabBar::tab:hover    {{ color:{THEME.TEXT}; }}"
        f"QSplitter::handle {{ background:{THEME.BORDER}; }}"
        f"QMenuBar {{ background:{THEME.BG}; color:{THEME.TEXT}; }}"
        f"QMenuBar::item:selected {{ background:{THEME.BUTTON_HOVER}; }}"
        f"QMenu {{ background:{THEME.PANEL}; color:{THEME.TEXT}; border:1px solid {THEME.BORDER}; }}"
        f"QMenu::item:selected {{ background:{THEME.ACCENT}; color:#ffffff; }}"
        f"QStatusBar {{ background:{THEME.BG}; color:{THEME.TEXT_DIM}; }}"
        f"QToolTip {{"
        f"  background:{THEME.TOOLTIP};"
        f"  color:{THEME.TEXT};"
        f"  border:1px solid {THEME.BORDER};"
        f"  border-radius:6px; padding:4px 6px;"
        f"}}"

        # --- Plain text view (used by the error dialog) ------------------
        f"QPlainTextEdit, QTextEdit {{"
        f"  background:{THEME.PANEL};"
        f"  color:{THEME.TEXT};"
        f"  border:1px solid {THEME.BORDER};"
        f"  border-radius:6px;"
        f"}}"
    )


# Main window


class TdcViewer(QMainWindow):
    DEFAULT_BINS = 200

    def __init__(
        self,
        hits: Optional[np.ndarray] = None,
        path: str = "",
        *,
        max_events: int = 0,
        daq_config: str = "",
        roc_filter: int = -1,
    ):
        super().__init__()
        self.setWindowTitle("Tagger TDC Viewer")
        self.resize(1280, 800)
        apply_theme_palette(self)
        self.setStyleSheet(_app_stylesheet())

        self._hits: np.ndarray = (
            hits if hits is not None else np.zeros(0, dtype=RECORD_DTYPE)
        )
        self._path = path
        self._current: Optional[Tuple[int, int]] = None
        # Channel A / B for event-wise correlations (Δt, A vs B).
        self._channel_a: Optional[Tuple[int, int]] = None
        self._channel_b: Optional[Tuple[int, int]] = None
        self._ch_names: Dict[Tuple[int, int], str] = load_channel_map()
        self._load_max_events = max_events
        self._load_daq_config = daq_config
        self._load_roc_filter = roc_filter

        # --- live-stream state ---
        self._stream = LiveStream(self)
        self._stream.hitsReceived.connect(self._on_live_hits)
        self._stream.stateChanged.connect(self._on_live_state)
        self._stream.statsUpdate.connect(self._on_live_stats)
        self._reset_live_state()
        # Rolling memory cap; drop the oldest half when exceeded.
        self._max_live_hits = 10_000_000
        self._last_live_url = "ws://localhost:5051"

        self._live_timer = QTimer(self)
        self._live_timer.setInterval(333)   # ~3 Hz repaint
        self._live_timer.timeout.connect(self._live_tick)

        self._build_ui()
        self._make_menu()

        if self._hits.size:
            self._rebuild_index()

    # --- UI layout -------------------------------------------------------

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        main_layout = QVBoxLayout(central)

        # Top row: file label + edge filter + bins
        top = QHBoxLayout()
        self.file_label = QLabel("(no file)")
        top.addWidget(self.file_label, 1)

        top.addWidget(QLabel("Edge:"))
        self.edge_combo = QComboBox()
        self.edge_combo.addItems(["both", "leading (0)", "trailing (1)"])
        self.edge_combo.currentIndexChanged.connect(self._refresh)
        top.addWidget(self.edge_combo)

        top.addWidget(QLabel("Bins:"))
        self.bins_spin = QSpinBox()
        self.bins_spin.setRange(10, 2000)
        self.bins_spin.setValue(self.DEFAULT_BINS)
        self.bins_spin.setSingleStep(10)
        self.bins_spin.valueChanged.connect(self._refresh)
        top.addWidget(self.bins_spin)

        # Live-stream controls — only useful when connected, but cheap to show.
        top.addSpacing(16)
        self.btn_pause = QPushButton("Pause")
        self.btn_pause.setCheckable(True)
        self.btn_pause.setEnabled(False)
        self.btn_pause.toggled.connect(self._toggle_pause)
        top.addWidget(self.btn_pause)
        self.btn_clear_live = QPushButton("Clear buffer")
        self.btn_clear_live.setToolTip("Drop all accumulated hits and start empty")
        self.btn_clear_live.clicked.connect(self._clear_live_buffer)
        top.addWidget(self.btn_clear_live)

        main_layout.addLayout(top)

        # Main splitter: left (tree) | right (plots)
        splitter = QSplitter(Qt.Orientation.Horizontal)

        self.tree = QTreeWidget()
        self.tree.setHeaderLabels(["Slot / Channel", "Hits", "Name"])
        self.tree.setColumnWidth(0, 130)
        self.tree.setColumnWidth(1, 70)
        self.tree.setColumnWidth(2, 80)
        self.tree.itemSelectionChanged.connect(self._on_tree_select)
        splitter.addWidget(self.tree)

        right = QWidget()
        rlay = QVBoxLayout(right)
        rlay.setContentsMargins(4, 4, 4, 4)

        # Pair selector row (applies to Δt and A-vs-B tabs).
        pair_row = QHBoxLayout()
        pair_row.setSpacing(6)

        def pair_selector(which: str) -> QLabel:
            pair_row.addWidget(QLabel(f"{which}:"))
            lbl = QLabel("—")
            lbl.setFrameShape(QFrame.Shape.StyledPanel)
            lbl.setMinimumWidth(110)
            pair_row.addWidget(lbl)
            btn = QPushButton(f"Set {which} ←")
            btn.setToolTip(
                f"Use the currently-selected tree channel as channel {which}")
            btn.clicked.connect(lambda _checked=False: self._set_from_tree(which))
            pair_row.addWidget(btn)
            return lbl

        self.lbl_a = pair_selector("A")
        pair_row.addSpacing(12)
        self.lbl_b = pair_selector("B")

        pair_row.addSpacing(6)
        btn_swap = QPushButton("Swap")
        btn_swap.clicked.connect(self._swap_ab)
        pair_row.addWidget(btn_swap)
        btn_clear = QPushButton("Clear")
        btn_clear.clicked.connect(self._clear_ab)
        pair_row.addWidget(btn_clear)
        pair_row.addStretch(1)
        rlay.addLayout(pair_row)

        self.channel_bar = BarChart()
        self.channel_bar.setTitle("Hits per channel (selected slot) — click a bar")
        self.channel_bar.barClicked.connect(self._on_bar_clicked)
        rlay.addWidget(self.channel_bar)

        # Tabbed plot area: single channel / Δt / A vs B.
        self.plot_tabs = QTabWidget()

        self.tdc_hist = Histogram()
        self.tdc_hist.setTitle("TDC value histogram — select a channel")
        self.tdc_hist.setXLabel("TDC value (LSB = 25 ps after rol2 shift)")
        self.plot_tabs.addTab(self.tdc_hist, "Single channel")

        self.diff_hist = Histogram()
        self.diff_hist.setTitle("Δt = A − B — set channel A and B")
        self.diff_hist.setXLabel("tdc(A) − tdc(B)")
        self.plot_tabs.addTab(self.diff_hist, "Δt = A − B")

        self.scatter_map = Heatmap2D()
        self.scatter_map.setTitle("tdc(A) vs tdc(B) — set channel A and B")
        self.scatter_map.setLabels("tdc(A)", "tdc(B)")
        self.plot_tabs.addTab(self.scatter_map, "A vs B (2D)")

        self.plot_tabs.currentChanged.connect(lambda _i: self._refresh())
        rlay.addWidget(self.plot_tabs, 1)

        splitter.addWidget(right)
        splitter.setStretchFactor(0, 0)
        splitter.setStretchFactor(1, 1)
        main_layout.addWidget(splitter, 1)

        self.setStatusBar(QStatusBar())

    def _make_menu(self):
        m = self.menuBar().addMenu("&File")
        a_open = QAction("&Open file…", self)
        a_open.setShortcut("Ctrl+O")
        a_open.triggered.connect(self._open_dialog)
        m.addAction(a_open)
        m.addSeparator()
        a_live = QAction("&Connect to prad2_server…", self)
        a_live.setShortcut("Ctrl+L")
        a_live.triggered.connect(self._connect_live_dialog)
        m.addAction(a_live)
        a_disc = QAction("&Disconnect", self)
        a_disc.triggered.connect(self._disconnect_live)
        m.addAction(a_disc)
        m.addSeparator()
        a_quit = QAction("&Quit", self)
        a_quit.setShortcut("Ctrl+Q")
        a_quit.triggered.connect(self.close)
        m.addAction(a_quit)

    # --- loading ---------------------------------------------------------

    def _open_dialog(self):
        filt = "EVIO files (*.evio *.evio.*);;All files (*)"
        path, _ = QFileDialog.getOpenFileName(self, "Open evio file", "", filt)
        if path:
            self.load(path)

    def load(self, path: str):
        # Loading a static file is mutually exclusive with the live stream.
        if self._live_timer.isActive():
            self._disconnect_live()

        # A second load while one is in flight cancels the first; its result
        # is ignored when it finishes.
        self._cancel_load()

        self.statusBar().showMessage(f"Loading {path}…")

        # Grey out the main GUI so the user can't interact with stale data
        # (tree clicks, A/B selection, etc.) while the loader runs.
        # Top-level dialogs stay responsive because they're separate windows.
        self.centralWidget().setEnabled(False)
        self.menuBar().setEnabled(False)

        # ``max`` starts at the --max-events cap, or a 1 M placeholder when
        # unlimited.
        pmax = int(self._load_max_events) if self._load_max_events > 0 else 1_000_000
        dlg = QProgressDialog(
            f"Decoding {os.path.basename(path)} …", "Cancel",
            0, pmax, self,
        )
        dlg.setWindowTitle("Loading")
        dlg.setWindowModality(Qt.WindowModality.WindowModal)
        dlg.setMinimumDuration(300)
        dlg.setAutoClose(True)
        dlg.setValue(0)

        worker = LoadWorker(
            path, self._load_max_events,
            self._load_daq_config, self._load_roc_filter,
        )

        # When the user passed an explicit cap, the dialog max is the real
        # limit — don't grow it. In unlimited mode (cap==0, dialog seeded
        # with a 1 M placeholder), grow proactively at 80 % so the bar
        # doesn't visibly snap from 100 % back to ~50 % on the next tick.
        unlimited = (self._load_max_events == 0)
        def _on_progress(events: int, hits: int):
            if unlimited:
                m = dlg.maximum()
                if m > 0 and events * 5 >= m * 4:
                    dlg.setMaximum(max(events * 2, m * 2))
            dlg.setValue(events)
            dlg.setLabelText(
                f"Decoding {os.path.basename(path)}\n"
                f"Events: {events:,}   Hits: {hits:,}"
            )

        worker.progressed.connect(_on_progress)

        self._load_dialog = dlg
        self._load_worker = worker
        self._load_token = path            # cancel-check: must match on finish
        # The dialog's Cancel button stops the worker at its next progress tick.
        self._load_thread = start_worker_thread(
            self, worker,
            lambda hits: self._on_load_finished(path, hits),
            lambda msg: self._on_load_failed(path, msg),
            dialog=dlg, on_thread_finished=dlg.close,
        )

    # --- load callbacks / cancellation -----------------------------------

    def _cancel_load(self):
        """Invalidate any load in flight. The worker may still be mid-loop;
        its eventual completion signal will see a stale token and no-op."""
        self._load_token = None
        if getattr(self, "_load_worker", None) is not None:
            try: self._load_worker.request_cancel()
            except Exception: pass
        if getattr(self, "_load_dialog", None) is not None:
            try: self._load_dialog.close()
            except Exception: pass
            self._load_dialog = None
        self._restore_ui()

    def _restore_ui(self):
        """Re-enable the main window after a load finishes / fails / cancels."""
        try:
            self.centralWidget().setEnabled(True)
            self.menuBar().setEnabled(True)
        except Exception:
            pass

    def _on_load_finished(self, path: str, hits: np.ndarray):
        if self._load_token != path:
            self._restore_ui()
            return                          # superseded by a later load()
        self._hits = hits
        self._path = path
        self._load_token = None
        self._rebuild_index()
        self.statusBar().showMessage(f"Loaded {hits.size:,} hits from {path}")
        self._restore_ui()

    def _on_load_failed(self, path: str, msg: str):
        if self._load_token != path:
            self._restore_ui()
            return
        self._load_token = None
        show_error_dialog(
            self,
            title="Load failed",
            heading=f"Could not load:  {path}",
            details=msg,
        )
        self.statusBar().showMessage("")
        self._restore_ui()

    # --- indexing --------------------------------------------------------

    def _rebuild_index(self, *, reset_selection: bool = True):
        """Rebuild the slot/channel tree from ``self._hits``.

        ``reset_selection=True`` (default) clears _current / A / B — used when
        loading a brand-new file.  Live-stream ticks pass False so the user's
        tree selection and A/B pair survive across rebuilds.
        """
        hits = self._hits
        # In live mode _update_live_status() rewrites this label.
        if self._path:
            self.file_label.setText(
                f"{os.path.basename(self._path)} — {hits.size:,} hits"
            )
        else:
            self.file_label.setText(f"(in-memory) — {hits.size:,} hits")

        saved_current = self._current
        saved_a = self._channel_a
        saved_b = self._channel_b

        self.tree.clear()

        if reset_selection:
            self._current = None
            self._channel_a = None
            self._channel_b = None
            self._update_ab_labels()

        if hits.size == 0:
            self._refresh()
            return

        slots = np.unique(hits["slot"])
        for slot in slots:
            smask = hits["slot"] == slot
            sub = hits[smask]
            slot_item = QTreeWidgetItem([f"slot {int(slot)}", f"{sub.size:,}", ""])
            slot_item.setData(0, Qt.ItemDataRole.UserRole, ("slot", int(slot)))
            self.tree.addTopLevelItem(slot_item)

            chs, counts = np.unique(sub["channel"], return_counts=True)
            for ch, c in zip(chs, counts):
                name = self._name_for(int(slot), int(ch))
                ch_item = QTreeWidgetItem([
                    f"  ch {int(ch):3d}", f"{int(c):,}", name,
                ])
                ch_item.setData(
                    0, Qt.ItemDataRole.UserRole, ("channel", int(slot), int(ch))
                )
                slot_item.addChild(ch_item)
            slot_item.setExpanded(False)

        # Re-select prior item when preserving, otherwise auto-pick busiest slot.
        if not reset_selection and saved_current is not None:
            self._current = saved_current
            self._channel_a = saved_a
            self._channel_b = saved_b
            self._update_ab_labels()
            self._select_current_in_tree()
        else:
            best_slot = int(
                max(slots, key=lambda s: int(np.count_nonzero(hits["slot"] == s)))
            )
            slot_it, _ = self._find_tree_items(best_slot)
            if slot_it is not None:
                self.tree.setCurrentItem(slot_it)

        self._refresh()

    def _find_tree_items(self, slot: int, ch: Optional[int] = None):
        """(slot item, channel item) for slot/ch, either None when missing.
        The slot item is expanded when a channel is looked up."""
        role = Qt.ItemDataRole.UserRole
        for i in range(self.tree.topLevelItemCount()):
            it = self.tree.topLevelItem(i)
            data = it.data(0, role)
            if not data or data[0] != "slot" or data[1] != slot:
                continue
            if ch is None:
                return it, None
            it.setExpanded(True)
            for j in range(it.childCount()):
                cdata = it.child(j).data(0, role)
                if cdata and cdata[0] == "channel" and cdata[2] == ch:
                    return it, it.child(j)
            return it, None
        return None, None

    def _select_current_in_tree(self):
        """Best-effort: find a tree item matching self._current and select it."""
        if self._current is None:
            return
        slot_it, ch_it = self._find_tree_items(*self._current)
        # Slot exists but channel gone — settle for the slot.
        target = ch_it if ch_it is not None else slot_it
        if target is not None:
            self.tree.setCurrentItem(target)

    # --- live streaming -------------------------------------------------

    def _connect_live_dialog(self):
        url, ok = QInputDialog.getText(
            self, "Connect to prad2_server",
            "WebSocket URL:", text=self._last_live_url,
        )
        if not ok or not url.strip():
            return
        self._open_live(url.strip())

    def _connect_live_dialog_auto(self):
        """Called by --live CLI flag after the event loop starts."""
        self._open_live(self._last_live_url)

    def _reset_live_state(self):
        # Batches accumulated between GUI ticks — flushed by _live_timer.
        self._live_batches: list = []
        self._live_rate_hz = 0.0
        self._live_dropped = 0
        self._live_flags = 0

    def _open_live(self, url: str):
        self._last_live_url = url
        # Any in-memory file data is replaced by the live stream.
        self._path = ""
        self._hits = np.zeros(0, dtype=RECORD_DTYPE)
        self._reset_live_state()
        self._rebuild_index(reset_selection=True)
        self.btn_pause.setEnabled(True)
        self.btn_pause.setChecked(False)
        self._stream.set_paused(False)
        self._stream.open(self._last_live_url)
        self._update_live_status()
        self._live_timer.start()

    def _disconnect_live(self):
        self._live_timer.stop()
        self._stream.close()
        self.btn_pause.setEnabled(False)
        self.btn_pause.setChecked(False)

    def _toggle_pause(self, paused: bool):
        self._stream.set_paused(paused)
        self.btn_pause.setText("Resume" if paused else "Pause")

    def _clear_live_buffer(self):
        # Drop everything — applies in both live and file modes.
        self._live_batches.clear()
        self._hits = np.zeros(0, dtype=RECORD_DTYPE)
        self._rebuild_index(reset_selection=False)
        self._update_live_status()

    def _on_live_hits(self, batch: np.ndarray):
        # Called from the Qt event loop on every binary frame. Just queue —
        # the timer flushes to self._hits at a steady 3 Hz.
        self._live_batches.append(batch)

    def _on_live_state(self, state: str):
        self.statusBar().showMessage(f"live: {state}")

    def _on_live_stats(self, s: dict):
        self._live_rate_hz = float(s.get("rate_hz", 0.0))
        self._live_dropped = int(s.get("dropped", 0))
        self._live_flags   = int(s.get("flags", 0))

    def _live_tick(self):
        if not self._live_batches:
            # Nothing new — just refresh stats in the file label / status bar.
            self._update_live_status()
            return
        batches = self._live_batches
        self._live_batches = []
        new_hits = (batches[0] if len(batches) == 1
                    else np.concatenate(batches))
        if self._hits.size == 0:
            self._hits = new_hits
        else:
            self._hits = np.concatenate([self._hits, new_hits])

        if self._hits.size > self._max_live_hits:
            self._hits = self._hits[self._hits.size // 2:].copy()

        self._rebuild_index(reset_selection=False)
        self._update_live_status()

    def _update_live_status(self):
        if not self._stream.is_active():
            return
        extra = ""
        if self._live_flags & 1:
            extra = f"  DROPPED {self._live_dropped} frames"
        self.statusBar().showMessage(
            f"live @ {self._last_live_url}  "
            f"hits={self._hits.size:,}  "
            f"rate={self._live_rate_hz:,.0f}/s{extra}"
        )
        self.file_label.setText(
            f"(live) {self._last_live_url} — {self._hits.size:,} hits  "
            f"({self._live_rate_hz:,.0f}/s)"
        )

    # --- interaction -----------------------------------------------------

    def _current_edge_mask(self) -> Optional[int]:
        idx = self.edge_combo.currentIndex()
        return None if idx == 0 else (idx - 1)

    def _on_tree_select(self):
        items = self.tree.selectedItems()
        if not items:
            return
        data = items[0].data(0, Qt.ItemDataRole.UserRole)
        if not data:
            return
        if data[0] == "slot":
            self._current = (data[1], None)
            self._refresh()
        elif data[0] == "channel":
            self._current = (data[1], data[2])
            self._refresh()

    def _on_bar_clicked(self, idx: int):
        if self._current is None:
            return
        slot, _ = self._current
        self._current = (slot, int(idx))
        _, ch_it = self._find_tree_items(slot, idx)
        if ch_it is not None:
            self.tree.setCurrentItem(ch_it)
            return
        # No such channel: drop the stale tree selection so that clicking
        # the previously selected channel's bar re-selects it.
        self.tree.clearSelection()
        self._refresh()

    # --- channel naming -------------------------------------------------

    def _name_for(self, slot: int, channel: int) -> str:
        """Return the human-readable name for (slot, channel), or '' if
        the map doesn't cover it."""
        return self._ch_names.get((int(slot), int(channel)), "")

    # --- A / B channel pairing ------------------------------------------

    def _fmt_pair(self, p: Optional[Tuple[int, int]]) -> str:
        if p is None:
            return "—"
        name = self._name_for(p[0], p[1])
        return f"slot {p[0]}, ch {p[1]}" + (f"  [{name}]" if name else "")

    def _update_ab_labels(self):
        self.lbl_a.setText(self._fmt_pair(self._channel_a))
        self.lbl_b.setText(self._fmt_pair(self._channel_b))

    def _tree_channel_selection(self) -> Optional[Tuple[int, int]]:
        items = self.tree.selectedItems()
        if not items:
            return None
        data = items[0].data(0, Qt.ItemDataRole.UserRole)
        if not data or data[0] != "channel":
            return None
        return (data[1], data[2])

    def _set_from_tree(self, which: str):
        """Make the selected tree channel channel ``which`` ("A" or "B")."""
        sel = self._tree_channel_selection()
        if sel is None:
            QMessageBox.information(
                self, "Pick a channel",
                "Select a CHANNEL node in the tree (expand a slot first), "
                f"then click Set {which}."
            )
            return
        if which == "A":
            self._channel_a = sel
        else:
            self._channel_b = sel
        self._update_ab_labels()
        self._refresh()

    def _swap_ab(self):
        self._channel_a, self._channel_b = self._channel_b, self._channel_a
        self._update_ab_labels()
        self._refresh()

    def _clear_ab(self):
        self._channel_a = None
        self._channel_b = None
        self._update_ab_labels()
        self._refresh()

    # --- top-level refresh ----------------------------------------------

    def _refresh(self):
        self._refresh_bar()
        # Only repaint the currently visible plot tab — the others will
        # refresh when the user clicks into them (QTabWidget currentChanged
        # is wired to _refresh, so they're up-to-date on demand).
        tab = self.plot_tabs.currentIndex() if hasattr(self, "plot_tabs") else 0
        if tab == 0:
            self._refresh_histogram()
        elif tab == 1:
            self._refresh_diff()
        elif tab == 2:
            self._refresh_scatter()

    def _refresh_bar(self):
        hits = self._hits
        if hits.size == 0 or self._current is None:
            self.channel_bar.setData(np.zeros(0))
            self.channel_bar.setTitle("Hits per channel")
            return

        slot, ch = self._current
        mask = hits["slot"] == slot
        edge_sel = self._current_edge_mask()
        if edge_sel is not None:
            mask = mask & (hits["edge"] == edge_sel)
        sub = hits[mask]

        # Auto-size the channel axis: look at *any* hit in this slot (before
        # the edge cut) so that switching "leading/trailing/both" doesn't
        # reshape the chart under the user.
        slot_hits = hits[hits["slot"] == slot]
        if slot_hits.size > 0:
            max_ch = int(slot_hits["channel"].max())
        else:
            max_ch = 0
        nbars = round_up_channels(max_ch)

        channels = sub["channel"].astype(np.int32)
        counts = np.bincount(channels, minlength=nbars)
        if counts.size > nbars:
            counts = counts[:nbars]

        self.channel_bar.setData(counts)
        self.channel_bar.setTitle(
            f"Hits per channel — slot {slot} ({nbars}-channel axis)"
        )
        self.channel_bar.setHighlight(ch if ch is not None else None)

    def _refresh_histogram(self):
        hits = self._hits
        if hits.size == 0 or self._current is None or self._current[1] is None:
            self.tdc_hist.setData(np.zeros(0), np.zeros(0))
            self.tdc_hist.setTitle("TDC value histogram — select a channel")
            self.statusBar().showMessage("")
            return

        slot, ch = self._current
        mask = (hits["slot"] == slot) & (hits["channel"] == ch)
        edge_sel = self._current_edge_mask()
        if edge_sel is not None:
            mask = mask & (hits["edge"] == edge_sel)
        sub = hits[mask]

        if sub.size == 0:
            self.tdc_hist.setData(np.zeros(0), np.zeros(0))
            self.tdc_hist.setTitle(
                f"TDC histogram — slot {slot}, ch {ch} (no hits)"
            )
            self.statusBar().showMessage("")
            return

        tdc_vals = sub["tdc"].astype(np.int64)
        tmin, tmax = _int_range(tdc_vals)
        nbins = self.bins_spin.value()
        counts, edges = np.histogram(tdc_vals, bins=nbins, range=(tmin, tmax + 1))
        self.tdc_hist.setData(counts, edges)

        edge_name = (
            "both edges"
            if edge_sel is None
            else ("leading" if edge_sel == 0 else "trailing")
        )
        name = self._name_for(slot, ch)
        name_tag = f" [{name}]" if name else ""
        title = (
            f"TDC histogram — slot {slot}, ch {ch}{name_tag}, {edge_name} "
            f"— {sub.size:,} hits, mean={tdc_vals.mean():.1f}, rms={tdc_vals.std():.1f}"
        )
        self.tdc_hist.setTitle(title)

        self.statusBar().showMessage(
            f"slot={slot} ch={ch}{name_tag}  n={sub.size:,}  "
            f"min={tmin}  max={tmax}  mean={tdc_vals.mean():.2f}"
        )

    def _matched_ab(self, widget, empty: tuple, prefix: str):
        """(a, b, tdc_a, tdc_b) for the event-matched A/B pair, or None
        after clearing ``widget`` with ``empty`` and titling it
        "``prefix`` — <reason>"."""
        a, b = self._channel_a, self._channel_b
        if self._hits.size == 0 or a is None or b is None:
            reason = "set both channel A and channel B"
        elif a == b:
            reason = "A and B must differ"
        else:
            t_a, t_b = _match_pair(self._hits, a, b, self._current_edge_mask())
            if t_a.size:
                return a, b, t_a, t_b
            reason = (f"0 matched events "
                      f"(A={self._fmt_pair(a)}, B={self._fmt_pair(b)})")
        widget.setData(*empty)
        widget.setTitle(f"{prefix} — {reason}")
        return None

    def _refresh_diff(self):
        """Event-wise Δt = tdc(A) − tdc(B) histogram."""
        m = self._matched_ab(self.diff_hist, (np.zeros(0), np.zeros(0)),
                             "Δt = A − B")
        if m is None:
            return
        a, b, t_a, t_b = m
        n = t_a.size

        dt = t_a - t_b
        dmin, dmax = _int_range(dt)
        nbins = self.bins_spin.value()
        counts, edges = np.histogram(dt, bins=nbins, range=(dmin, dmax + 1))
        self.diff_hist.setData(counts, edges)
        self.diff_hist.setTitle(
            f"Δt = A − B  (A={self._fmt_pair(a)},  B={self._fmt_pair(b)}) "
            f"— {n:,} events, mean={dt.mean():.1f}, rms={dt.std():.1f}"
        )
        self.statusBar().showMessage(
            f"Δt: {n:,} matched events  "
            f"min={dmin}  max={dmax}  mean={dt.mean():.2f}  rms={dt.std():.2f}"
        )

    def _refresh_scatter(self):
        """Event-wise 2-D histogram of tdc(A) vs tdc(B)."""
        m = self._matched_ab(
            self.scatter_map, (np.zeros((0, 0)), np.zeros(0), np.zeros(0)),
            "tdc(A) vs tdc(B)")
        if m is None:
            return
        a, b, t_a, t_b = m
        n = t_a.size

        # 2-D gets coarser binning than the 1-D plots — each cell needs a
        # reasonable event count or the map looks like shot noise.
        nbins = max(10, self.bins_spin.value() // 4)
        amin, amax = _int_range(t_a)
        bmin, bmax = _int_range(t_b)
        counts, xedges, yedges = np.histogram2d(
            t_a, t_b, bins=nbins, range=[[amin, amax + 1], [bmin, bmax + 1]]
        )
        self.scatter_map.setData(counts, xedges, yedges)
        self.scatter_map.setLabels(
            f"tdc(A) — {self._fmt_pair(a)}",
            f"tdc(B) — {self._fmt_pair(b)}",
        )
        rho = float(np.corrcoef(t_a, t_b)[0, 1]) if n > 1 else 0.0
        self.scatter_map.setTitle(
            f"tdc(A) vs tdc(B)  — {n:,} events, corr ρ={rho:+.3f}"
        )
        self.statusBar().showMessage(
            f"2D: {n:,} events  "
            f"A:[{amin},{amax}]  B:[{bmin},{bmax}]  ρ={rho:+.3f}"
        )

    # --- cleanup --------------------------------------------------------

    def closeEvent(self, ev):
        # Make sure the WebSocket is torn down cleanly so the server sees our
        # unsubscribe and the subs counter reaches zero.
        try:
            if self._live_timer.isActive():
                self._disconnect_live()
        except Exception:
            pass
        # Cancel any in-flight file load so its eventual completion signal
        # is a no-op instead of touching a destroyed viewer.
        self._cancel_load()
        super().closeEvent(ev)


# Entry point


def _parse_roc(value: str) -> int:
    s = value.strip().lower()
    if not s:
        return -1
    return int(s, 0)  # handles 0x..., 0..., decimal


def _cli_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="tagger_viewer.py",
        description="Interactive viewer for V1190 TDC hits (bank 0xE107).",
    )
    p.add_argument(
        "path",
        nargs="?",
        help="Offline mode: evio file (*.evio / *.evio.*) decoded via prad2py.",
    )
    p.add_argument(
        "-n", "--max-events",
        type=int,
        default=1_000_000,
        help="When reading an .evio file, stop after N physics events "
             "(default 1,000,000; pass 0 for unlimited).",
    )
    p.add_argument(
        "-D", "--daq-config",
        default="",
        help="Override daq_config.json path (prad2py uses the installed default otherwise).",
    )
    p.add_argument(
        "--roc",
        default="",
        help="Only keep hits from this parent ROC tag, e.g. --roc 0x8E (default: all).",
    )
    p.add_argument(
        "--live",
        default="",
        metavar="URL",
        help="Auto-connect to a prad2_server WebSocket on startup, "
             "e.g. --live ws://clondaq6:5051",
    )
    p.add_argument(
        "--no-smoke-test",
        action="store_true",
        help="Skip the pre-flight subscribe/ack handshake against --live URL "
             "(the round-trip defaults to ON and aborts startup on failure).",
    )
    p.add_argument(
        "--theme", choices=available_themes(), default="dark",
        help="Colour theme (default: dark)",
    )
    return p


def smoke_test_live(url: str, timeout_ms: int = 5000) -> Optional[str]:
    """Open a WebSocket, send tagger_subscribe, wait for the tagger_subscribed
    acknowledgement.  Returns None on success, a human-readable error string
    otherwise.  Uses a local QEventLoop so it's safe to call after
    QApplication() is constructed but before the main window is shown.
    """
    from PyQt6.QtCore import QEventLoop

    loop = QEventLoop()
    stream = LiveStream()
    state = {"error": "timeout (no tagger_subscribed ack)"}

    def finish(err: Optional[str]):
        state["error"] = err
        if loop.isRunning():
            loop.quit()

    stream.subscribed.connect(lambda: finish(None))
    stream.failed.connect(finish)
    QTimer.singleShot(timeout_ms, lambda: finish(state["error"]))
    stream.open(url)
    loop.exec()

    # Unsubscribe and close so the server doesn't see a dangling
    # subscriber hanging around until we reconnect for real.
    stream.close()
    return state["error"]


def main(argv):
    args = _cli_parser().parse_args(argv[1:])
    roc = _parse_roc(args.roc)

    set_theme(args.theme)

    app = QApplication(argv[:1])

    # Pre-flight check for --live: fast subscribe/ack round-trip.  Bail out
    # before building the GUI if the server isn't reachable or our protocol
    # doesn't match — saves the user from a blank window with no feedback.
    if args.live and not args.no_smoke_test:
        err = smoke_test_live(args.live)
        if err:
            sys.stderr.write(f"tagger_viewer: cannot connect to {args.live}: {err}\n")
            show_error_dialog(
                None,
                title="Live stream unreachable",
                heading=f"Cannot connect to:  {args.live}",
                details=err,
            )
            return 1

    win = TdcViewer(
        hits=None,
        path="",
        max_events=args.max_events,
        daq_config=args.daq_config,
        roc_filter=roc,
    )
    # Kick off file / live connection after the event loop is actually
    # running — that way both use the async load path with the progress
    # dialog rather than freezing the GUI on startup.
    if args.path:
        QTimer.singleShot(0, lambda: win.load(args.path))
    if args.live:
        win._last_live_url = args.live
        QTimer.singleShot(0, win._connect_live_dialog_auto)
    win.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main(sys.argv))
