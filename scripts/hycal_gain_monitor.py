#!/usr/bin/env python3
"""
HyCal Gain Monitor (PyQt6)
==========================
Visualises LMS-based gain factors across runs for all HyCal modules.
Reads text-based ``prad_{:06d}_LMS.dat`` files produced by the offline
gain analysis, displays a colour-coded HyCal geo map, LMS reference
channel stability charts, and a table of irregular (module, run) entries.

Usage
-----
    python scripts/hycal_gain_monitor.py [-dir FOLDER] [--theme THEME]
"""

from __future__ import annotations

import glob
import html
import math
import os
import re
import shutil
import sys
from collections import OrderedDict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

try:
    import uproot
    _UPROOT_OK = True
except ImportError:
    _UPROOT_OK = False

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QComboBox, QLineEdit, QSpinBox,
    QFileDialog, QSplitter, QSizePolicy, QTableWidget,
    QTableWidgetItem, QHeaderView, QAbstractItemView, QMenu,
    QDialog, QFormLayout, QTextEdit, QMessageBox, QCheckBox,
)
from PyQt6.QtCore import (
    Qt, QRectF, QPointF, pyqtSignal, QTimer, QProcess, QProcessEnvironment,
)
from PyQt6.QtGui import (
    QPainter, QColor, QPen, QFont, QPalette,
)

from hycal_geoview import (
    Module, load_modules, HyCalMapWidget, ZoomHistWidget, cmap_qcolor,
    nice_ticks, apply_theme_palette, set_theme, available_themes, THEME,
    themed, OVERLAY_BUTTON_QSS,
)
from hycal_calib import LMSRecord, ModuleRecord, read_lms_dat
from evio_io import (
    EVIO_BYTES_PER_FILE_EST, LOCAL_DATA_BASE, REMOTE_DATA_BASE, REMOTE_HOST,
    check_disk_space, fmt_bytes, free_bytes, local_evio_in_range, scp_bash,
)


# ---- Paths & constants ----

SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = SCRIPT_DIR / ".." / "database"
MODULES_JSON = DB_DIR / "hycal_map.json"

LMS_NAMES = ["LMS1", "LMS2", "LMS3"]              # data keys (file/ROOT/geometry)
LMS_DISPLAY = ["Ref1", "Ref2", "Ref3"]            # user-facing labels for the same modules
LMS_REF_DEFAULT = 1          # index into LMS_NAMES -> "LMS2" (shown as "Ref2")
FILE_PATTERN = re.compile(r"prad_(\d{6})_LMS\.dat$")


def lms_display_name(name: Optional[str]) -> str:
    """Map an LMS reference data key (LMS1/2/3) to its UI label (Ref1/2/3).

    Non-LMS module names pass through unchanged so callers can use this
    blindly when formatting any module name for display."""
    if not name:
        return name or ""
    try:
        return LMS_DISPLAY[LMS_NAMES.index(name)]
    except ValueError:
        return name

# Upper-right charts (fixed set, in display order).  Each tuple is
# (kind, label).  `kind` drives data extraction in _update_line_charts.
#   ratio → ref PMT's LMS/Alpha (selected-ref)
#   ref_lms → ref PMT's LMS peak±σ (selected-ref)
#   ref_alpha → ref PMT's Alpha peak±σ (selected-ref)
#   mod_lms → selected HyCal module's LMS peak±σ
#   mod_gain → selected HyCal module's gain factor wrt selected ref
CHART_PLOTS: List[Tuple[str, str]] = [
    ("ratio",     "Ref LMS/Alpha"),
    ("ref_lms",   "Ref LMS"),
    ("ref_alpha", "Ref Alpha"),
    ("mod_lms",   "Module LMS"),
    ("mod_gain",  "Module Gain"),
]

# Default palette for non-drift modes.
_DEFAULT_PALETTE = "blue-orange"

# Separate palette used only in Run-to-Run Drift mode (not cycled by the user)
DRIFT_PALETTE = [
    (0.00, (0, 210, 230)),   # cyan  — large negative drift
    (0.50, (80, 80, 80)),    # grey  — no drift (always maps to 0 in drift mode)
    (1.00, (249, 115, 22)),  # orange — large positive drift
]

# Stylesheets in legacy dark-theme hex: pass them through themed() when the
# widget is built, since main() selects the theme after import.
_COMBO_QSS = (
    "QComboBox{background:#161b22;color:#c9d1d9;"
    "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
    "QComboBox::drop-down{border:none;width:18px;}"
    "QComboBox::down-arrow{border-left:4px solid transparent;"
    "border-right:4px solid transparent;border-top:5px solid #8b949e;"
    "margin-right:4px;}"
    "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
    "border:1px solid #30363d;selection-background-color:#1f6feb;}")
_EDIT_QSS = (
    "QLineEdit{background:#161b22;color:#c9d1d9;"
    "border:1px solid #30363d;border-radius:3px;padding:2px 4px;}")
_BTN_QSS = (
    "QPushButton{background:#21262d;color:#c9d1d9;"
    "border:1px solid #30363d;padding:4px 8px;"
    "font:bold 11px Consolas;border-radius:3px;}"
    "QPushButton:hover{background:#30363d;}")
_TOGGLE_QSS = _BTN_QSS + (
    "QPushButton:checked{background:#1f6feb;color:white;"
    "border-color:#388bfd;}")


def _slabel(text: str) -> QLabel:
    lbl = QLabel(text)
    lbl.setFont(QFont("Consolas", 10))
    lbl.setStyleSheet(themed("color:#c9d1d9;"))
    return lbl


# ---- Data structures ----

@dataclass
class RunData:
    run_number: int
    lms: Dict[str, LMSRecord] = field(default_factory=dict)
    modules: Dict[str, ModuleRecord] = field(default_factory=dict)


@dataclass
class IrregularEntry:
    name: str
    mod_type: str
    run_number: int
    gain: float
    mean_gain: float
    std_dev: float
    deviation_sigma: float


@dataclass
class DriftEntry:
    name: str
    mod_type: str
    run_number: int
    prev_run_number: int
    gain_current: float
    gain_prev: float
    rel_change: float     # (gain_current - gain_prev) / gain_prev


@dataclass
class SummaryEntry:
    name: str
    mod_type: str
    drift_count: int      # number of consecutive run pairs with |Δ| > threshold
    max_rel_change: float # largest |Δ| seen (absolute value)
    max_run: int          # run where max drift occurred
    max_prev_run: int     # previous run for that pair


# ---- File parsing ----

def parse_dat_file(filepath: str) -> Optional[RunData]:
    """Parse a single prad_NNNNNN_LMS.dat file."""
    m = FILE_PATTERN.search(os.path.basename(filepath))
    if not m:
        return None
    parsed = read_lms_dat(filepath)
    if parsed is None or not (parsed[0] or parsed[1]):
        return None
    return RunData(run_number=int(m.group(1)), lms=parsed[0], modules=parsed[1])


def load_all_runs(folder: str) -> List[RunData]:
    """Scan folder for prad_*_LMS.dat files, parse all, sort by run number."""
    pattern = os.path.join(folder, "prad_*_LMS.dat")
    files = sorted(glob.glob(pattern))
    runs: List[RunData] = []
    for f in files:
        rd = parse_dat_file(f)
        if rd is not None:
            runs.append(rd)
    # files are already sorted alphabetically; zero-padded run numbers preserve numerical order
    return runs


# ---- Outlier detection ----

def gain_stats(runs: List[RunData],
               ref_idx: int) -> Dict[str, Tuple[float, float, int]]:
    """Per-module (mean, population std, number of runs) of the gain factor
    against reference ``ref_idx`` over ``runs``."""
    gains: Dict[str, List[float]] = {}
    for rd in runs:
        for mname, mrec in rd.modules.items():
            gains.setdefault(mname, []).append(mrec.gain_factors[ref_idx])
    stats: Dict[str, Tuple[float, float, int]] = {}
    for mname, glist in gains.items():
        mean = sum(glist) / len(glist)
        variance = sum((g - mean) ** 2 for g in glist) / len(glist)
        stats[mname] = (mean, math.sqrt(variance) if variance > 0 else 0.0,
                        len(glist))
    return stats


def compute_irregular_entries(
    runs: List[RunData],
    ref_idx: int,
    mod_by_name: Dict[str, Module],
    sigma_threshold: float = 3.0,
    min_runs: int = 5,
    stats: Optional[Dict[str, Tuple[float, float, int]]] = None,
) -> List[IrregularEntry]:
    """Find (module, run) pairs with outlier gain factors.  ``stats`` is
    gain_stats(runs, ref_idx), computed here when not given."""
    if stats is None:
        stats = gain_stats(runs, ref_idx)
    entries: List[IrregularEntry] = []
    for rd in runs:
        for mname, mrec in rd.modules.items():
            mean, std, n = stats[mname]
            if n < min_runs or std == 0:
                continue
            gain = mrec.gain_factors[ref_idx]
            dev = abs(gain - mean) / std
            if dev > sigma_threshold:
                mod = mod_by_name.get(mname)
                entries.append(IrregularEntry(
                    name=mname,
                    mod_type=mod.mod_type if mod else "?",
                    run_number=rd.run_number,
                    gain=gain,
                    mean_gain=mean,
                    std_dev=std,
                    deviation_sigma=dev,
                ))

    entries.sort(key=lambda e: (e.name, e.run_number))
    return entries


def _sym_rel_change(g_curr: float, g_prev: float) -> float:
    """|g_curr - g_prev| relative to the smaller of the two gains (inf if
    either is 0): the drift measure the thresholds apply to."""
    denom = min(abs(g_curr), abs(g_prev))
    return math.inf if denom == 0 else abs(g_curr - g_prev) / denom


def _drift_threshold(name: str, thresh_g: float, thresh_w: float) -> float:
    return thresh_g if name.startswith("G") else thresh_w


def compute_drift_entries(
    rd_curr: "RunData",
    rd_prev: "RunData",
    ref_idx: int,
    mod_by_name: Dict[str, "Module"],
    thresh_g: float = 0.10,
    thresh_w: float = 0.05,
) -> List[DriftEntry]:
    """Find modules where gain changed by more than threshold relative to previous run."""
    entries: List[DriftEntry] = []
    for mname, mrec in rd_curr.modules.items():
        prev_mrec = rd_prev.modules.get(mname)
        if prev_mrec is None:
            continue
        g_curr = mrec.gain_factors[ref_idx]
        g_prev = prev_mrec.gain_factors[ref_idx]
        rel_display = math.inf if g_prev == 0 else (g_curr - g_prev) / g_prev
        rel_sym = _sym_rel_change(g_curr, g_prev)
        if (math.isinf(rel_sym)
                or rel_sym > _drift_threshold(mname, thresh_g, thresh_w)):
            mod = mod_by_name.get(mname)
            entries.append(DriftEntry(
                name=mname,
                mod_type=mod.mod_type if mod else "?",
                run_number=rd_curr.run_number,
                prev_run_number=rd_prev.run_number,
                gain_current=g_curr,
                gain_prev=g_prev,
                rel_change=rel_display,
            ))

    def sort_key(e):
        rel = _sym_rel_change(e.gain_current, e.gain_prev)
        return (0 if e.name.startswith("W") else 1,
                math.isinf(e.rel_change),
                0 if math.isinf(rel) else -rel)

    entries.sort(key=sort_key)
    return entries


# ---- HyCal Gain Map Widget ----

# Map legend swatches of each view: (palette position, label).
_LEGEND_ITEMS = {
    "drift": [
        (0.0, "gain decreases"),
        (0.5, "stable"),
        (1.0, "gain increases"),
    ],
    "summary": [
        (0.0, "low drift count"),
        (1.0, "high drift count"),
    ],
    "gain": [
        (0.0, "low gain"),
        (1.0, "high gain"),
    ],
    "deviation": [
        (0.0, "below mean"),
        (0.5, "near mean"),
        (1.0, "above mean"),
    ],
}

# _LEGEND_ITEMS key of each GainMonitorWindow view mode (View combo index).
_VIEW_MODE_LEGEND = ("gain", "deviation", "drift", "summary")


class HyCalGainMapWidget(HyCalMapWidget):
    """Gain-monitor specialisation of the shared HyCal map widget.

    Adds a custom palette override (used by Run-to-Run Drift mode) and a
    legend overlay above the colour bar that explains the active view
    mode.  The LMS cells carry their Ref1/2/3 labels, and a module click
    toggles the selection highlight.
    """

    CB_MAX_WIDTH = 300

    def __init__(self, parent=None):
        super().__init__(parent, shrink=0.90, margin_top=8,
                         enable_zoom_pan=True, include_lms=True,
                         label_types={"LMS"}, toggle_select=True)
        self._palette_override = None
        self._legend_mode: Optional[str] = None

    # -- public API additions --

    def set_gain_data(self, values: Dict[str, float],
                      vmin: float, vmax: float):
        self._values = values
        self._vmin = vmin
        self._vmax = vmax
        self.update()

    def set_palette(self, idx_or_name):
        self._palette_override = None
        super().set_palette(idx_or_name)

    def set_palette_override(self, stops):
        """Use a custom stops list instead of the indexed palette. Pass None to clear."""
        self._palette_override = stops
        self.update()

    def set_legend_mode(self, mode: Optional[str]):
        if mode != self._legend_mode:
            self._legend_mode = mode
            self.update()

    # -- base hooks --

    def palette_stops(self):
        if self._palette_override is not None:
            return self._palette_override
        return super().palette_stops()

    def _fmt_value(self, v: float) -> str:
        return f"{v:.4f}"

    def _tooltip_text(self, name: str) -> str:
        label = lms_display_name(name)
        v = self._values.get(name)
        if v is None:
            return label
        return f"{label}: {v:.5f}"

    def _label_text(self, name: str) -> str:
        return lms_display_name(name)

    def _paint_empty(self, p, w, h):
        if not self._values:
            p.setPen(QColor(THEME.TEXT_MUTED))
            p.setFont(QFont("Consolas", 12))
            p.drawText(QRectF(0, 0, w, h),
                       Qt.AlignmentFlag.AlignCenter, "No data loaded")

    def _paint_after_colorbar(self, p, w, h):
        items = _LEGEND_ITEMS.get(self._legend_mode)
        if not items or self._cb_rect is None:
            return
        cb_y = self._cb_rect.y()
        p.setFont(QFont("Consolas", 9))
        fm = p.fontMetrics()
        swatch = 12
        gap = 5
        item_w = swatch + gap + max(fm.horizontalAdvance(lbl) for _, lbl in items)
        spacing = 18
        total_w = len(items) * item_w + (len(items) - 1) * spacing
        lh = max(swatch, fm.height())
        pad = 4
        lx = (w - total_w) // 2 - pad
        ly = cb_y - lh - 2 * pad - 4
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor(10, 14, 20, 200))
        p.drawRoundedRect(QRectF(lx, ly, total_w + 2 * pad, lh + 2 * pad), 4, 4)
        x = lx + pad
        stops = self.palette_stops()
        for t, label in items:
            sy = ly + pad + (lh - swatch) // 2
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(cmap_qcolor(t, stops))
            p.drawRect(QRectF(x, sy, swatch, swatch))
            p.setPen(QColor(THEME.TEXT))
            p.drawText(QRectF(x + swatch + gap, ly + pad,
                              fm.horizontalAdvance(label), lh),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       label)
            x += item_w + spacing


def _collect_series(runs: List[RunData], records: list, point):
    """(1-based positions, run numbers, values, errors) of a chart series:
    ``point(rec)`` gives (value, error) for the record of each run, or None
    to skip the run; an error of None adds no error bar."""
    indices: List[int] = []
    actual_runs: List[int] = []
    vals: List[float] = []
    errs: List[float] = []
    for idx, (rd, rec) in enumerate(zip(runs, records)):
        pt = point(rec)
        if pt is None:
            continue
        indices.append(idx + 1)
        actual_runs.append(rd.run_number)
        vals.append(pt[0])
        if pt[1] is not None:
            errs.append(pt[1])
    return indices, actual_runs, vals, errs


def _chart_y_range(values: List[float], errors: List[float]) -> Tuple[float, float]:
    """Return y-axis (lo, hi): mean±20%, expanded if any data point falls outside."""
    finite = [(j, v) for j, v in enumerate(values) if math.isfinite(v)]
    if not finite:
        return 0.9, 1.1
    finite_vals = [v for _, v in finite]
    mean = sum(finite_vals) / len(finite_vals)
    y_lo = mean * 0.8
    y_hi = mean * 1.2
    for j, v in finite:
        err = errors[j] if j < len(errors) else 0.0
        y_lo = min(y_lo, v - err)
        y_hi = max(y_hi, v + err)
    return y_lo, y_hi


# ---- LMS Line Chart Widget ----

class LMSLineChartWidget(QWidget):
    """Line chart with error bars vs run number."""

    PAD_L, PAD_R, PAD_T, PAD_B = 60, 16, 24, 32

    runClicked = pyqtSignal(int)              # emits actual run number on left-click
    pointDeleteRequested = pyqtSignal(int)    # emits actual run number after right-click → "Delete run N"
    pointBackupRequested = pyqtSignal(int)    # emits actual run number after right-click → "Move run N to backup"

    def __init__(self, parent=None):
        super().__init__(parent)
        self._run_numbers: List[int] = []
        self._actual_run_numbers: List[int] = []
        self._ratios: List[float] = []
        self._errors: List[float] = []
        self._title: str = ""
        self._hover_idx: int = -1
        self._current_run_number: int = -1
        self._y_range: Optional[Tuple[float, float]] = None
        self._series_color: Optional[QColor] = None
        self.setMinimumHeight(80)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMouseTracking(True)

    def set_y_range(self, lo: float, hi: float):
        self._y_range = (lo, hi)
        self.update()

    def set_series_color(self, color: Optional[QColor]):
        """Override the default ACCENT-blue series colour. Pass None to clear."""
        self._series_color = QColor(color) if color is not None else None
        self.update()

    def set_current_run(self, run_number: int):
        if run_number != self._current_run_number:
            self._current_run_number = run_number
            self.update()

    def set_data(self, run_numbers: List[int], ratios: List[float],
                 errors: List[float], title: str,
                 actual_run_numbers: List[int] = None):
        self._run_numbers = run_numbers
        self._actual_run_numbers = actual_run_numbers if actual_run_numbers is not None else run_numbers
        self._ratios = ratios
        self._errors = errors
        self._title = title
        self._hover_idx = -1
        self.update()

    def _screen_xs(self, w: int) -> List[float]:
        """Return screen x-coordinates for all data points."""
        runs = self._run_numbers
        if not runs:
            return []
        px = self.PAD_L
        pw = w - self.PAD_L - self.PAD_R
        x_min, x_max = runs[0], runs[-1]
        if x_min == x_max:
            x_min -= 1; x_max += 1
        return [px + (r - x_min) / (x_max - x_min) * pw for r in runs]

    def _nearest_index(self, mx: float) -> int:
        """Index of the point nearest to widget x ``mx``; -1 if none lies
        within 20 px."""
        xs = self._screen_xs(self.width())
        if not xs:
            return -1
        i = min(range(len(xs)), key=lambda k: abs(mx - xs[k]))
        return i if abs(mx - xs[i]) < 20 else -1

    def mouseMoveEvent(self, event):
        new_idx = self._nearest_index(event.position().x())
        if new_idx != self._hover_idx:
            self._hover_idx = new_idx
            self.update()

    def mousePressEvent(self, event):
        btn = event.button()
        if btn not in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            return
        i = self._nearest_index(event.position().x())
        if i < 0 or i >= len(self._actual_run_numbers):
            return
        actual_rn = self._actual_run_numbers[i]
        if btn == Qt.MouseButton.LeftButton:
            self.runClicked.emit(actual_rn)
        else:
            menu = QMenu(self)
            backup_act = menu.addAction(f"Move run {actual_rn} to backup")
            backup_act.triggered.connect(
                lambda _checked=False, rn=actual_rn: self.pointBackupRequested.emit(rn))
            del_act = menu.addAction(f"Delete run {actual_rn}")
            del_act.triggered.connect(
                lambda _checked=False, rn=actual_rn: self.pointDeleteRequested.emit(rn))
            menu.exec(event.globalPosition().toPoint())

    def leaveEvent(self, event):
        if self._hover_idx != -1:
            self._hover_idx = -1
            self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor(THEME.CANVAS))

        if self._series_color is not None:
            series_color = QColor(self._series_color)
        else:
            series_color = QColor(THEME.ACCENT)
        p.setPen(series_color)
        p.setFont(QFont("Consolas", 10, QFont.Weight.Bold))
        p.drawText(QRectF(self.PAD_L, 2, w - self.PAD_L - self.PAD_R, 20),
                   Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                   self._title)

        runs = self._run_numbers
        ratios = self._ratios
        errors = self._errors
        if not runs or not ratios:
            p.setPen(QColor(THEME.TEXT_MUTED))
            p.setFont(QFont("Consolas", 10))
            p.drawText(QRectF(0, 0, w, h),
                       Qt.AlignmentFlag.AlignCenter, "No data")
            p.end()
            return

        px = self.PAD_L
        py = self.PAD_T
        pw = w - self.PAD_L - self.PAD_R
        ph = h - self.PAD_T - self.PAD_B
        if pw < 20 or ph < 20:
            p.end()
            return

        x_min, x_max = runs[0], runs[-1]
        if x_min == x_max:
            x_min -= 1
            x_max += 1

        # set_y_range() always follows a non-empty set_data().
        y_lo, y_hi = self._y_range
        if not math.isfinite(y_lo) or not math.isfinite(y_hi):
            p.setPen(QColor(THEME.TEXT_MUTED))
            p.setFont(QFont("Consolas", 10))
            p.drawText(QRectF(0, 0, w, h), Qt.AlignmentFlag.AlignCenter,
                       "No valid fit data")
            p.end()
            return

        if y_hi == y_lo:
            y_lo -= 0.05
            y_hi += 0.05

        def to_sx(v):
            return px + (v - x_min) / (x_max - x_min) * pw

        def to_sy(v):
            return py + ph - (v - y_lo) / (y_hi - y_lo) * ph

        p.setPen(QPen(QColor(THEME.BUTTON), 1, Qt.PenStyle.DotLine))
        y_ticks = nice_ticks(y_lo, y_hi, 5)
        for yt in y_ticks:
            sy = to_sy(yt)
            p.drawLine(QPointF(px, sy), QPointF(px + pw, sy))

        p.setPen(QPen(QColor(THEME.BORDER), 1))
        p.drawLine(QPointF(px, py), QPointF(px, py + ph))
        p.drawLine(QPointF(px, py + ph), QPointF(px + pw, py + ph))

        # y-axis labels
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Consolas", 8))
        for yt in y_ticks:
            sy = to_sy(yt)
            p.drawText(QRectF(0, sy - 8, self.PAD_L - 4, 16),
                       Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                       f"{yt:.3f}")

        # x-axis labels (skip some if too dense)
        max_labels = max(pw // 60, 2)
        step = max(len(runs) // max_labels, 1)
        for i in range(0, len(runs), step):
            sx = to_sx(runs[i])
            p.drawText(QRectF(sx - 30, py + ph + 2, 60, 18),
                       Qt.AlignmentFlag.AlignCenter, str(runs[i]))

        p.setPen(QPen(QColor(THEME.TEXT_DIM), 1))
        cap = 3
        for i in range(len(runs)):
            if i >= len(ratios):
                break
            sx = to_sx(runs[i])
            r = ratios[i]
            err = errors[i] if i < len(errors) else 0
            sy_top = to_sy(r + err)
            sy_bot = to_sy(r - err)
            p.drawLine(QPointF(sx, sy_top), QPointF(sx, sy_bot))
            p.drawLine(QPointF(sx - cap, sy_top), QPointF(sx + cap, sy_top))
            p.drawLine(QPointF(sx - cap, sy_bot), QPointF(sx + cap, sy_bot))

        p.setPen(QPen(series_color, 1.5))
        for i in range(len(runs) - 1):
            if i + 1 >= len(ratios):
                break
            p.drawLine(QPointF(to_sx(runs[i]), to_sy(ratios[i])),
                       QPointF(to_sx(runs[i + 1]), to_sy(ratios[i + 1])))

        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(series_color)
        for i in range(len(runs)):
            if i >= len(ratios):
                break
            radius = 5 if i == self._hover_idx else 3
            p.drawEllipse(QPointF(to_sx(runs[i]), to_sy(ratios[i])), radius, radius)

        # current run red circle
        actual = self._actual_run_numbers
        if self._current_run_number >= 0 and actual:
            for i, rn in enumerate(actual):
                if rn == self._current_run_number and i < len(runs) and i < len(ratios):
                    cx = to_sx(runs[i])
                    cy = to_sy(ratios[i])
                    p.setPen(QPen(QColor(THEME.DANGER), 2))
                    p.setBrush(Qt.BrushStyle.NoBrush)
                    p.drawEllipse(QPointF(cx, cy), 8, 8)
                    break

        # hover tooltip
        hi = self._hover_idx
        if 0 <= hi < len(runs) and hi < len(ratios):
            sx = to_sx(runs[hi])
            sy = to_sy(ratios[hi])
            actual_rn = self._actual_run_numbers[hi] if hi < len(self._actual_run_numbers) else runs[hi]
            tip = f"run {actual_rn}\n{ratios[hi]:.4f}"
            p.setFont(QFont("Consolas", 9))
            fm = p.fontMetrics()
            lines = tip.split("\n")
            tw = max(fm.horizontalAdvance(ln) for ln in lines)
            th = fm.height() * len(lines) + 6
            tx = sx + 10
            ty = sy - th // 2
            if tx + tw + 8 > px + pw:
                tx = sx - tw - 14
            if ty < py:
                ty = py
            if ty + th > py + ph:
                ty = py + ph - th
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(20, 30, 45, 220))
            p.drawRoundedRect(QRectF(tx - 4, ty - 2, tw + 8, th), 4, 4)
            p.setPen(QColor(THEME.TEXT_STRONG))
            for j, ln in enumerate(lines):
                p.drawText(QRectF(tx, ty + j * fm.height(), tw, fm.height()),
                           Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                           ln)

        p.end()


# ---- ROOT histogram display widget ----

_HIST_CACHE_MAX = 20


class RootHistWidget(ZoomHistWidget):
    """Displays a TH1 histogram read from a fitted LMS ROOT file.

    Left-drag to zoom into an x range; right-click → Unzoom resets to the
    histogram edges.
    """

    PAD_B = 40

    _X_DEFAULT_LO = 0.0
    _X_DEFAULT_HI = 1000.0

    # Overlay histogram colours (used for the LMS reference cells, where
    # we render both the LMS distribution and the corresponding alpha
    # peak in the same plot).
    OVERLAY_BAR_COLOR  = "#3a86ff"   # blue — alpha histogram bars
    OVERLAY_FIT_COLOR  = "#ffeb3b"   # yellow — alpha gaussian fit

    # Distinct, well-separated bar colours for stack mode.  Cycles when more
    # than len() runs are stacked.  Alpha is applied at draw time.
    STACK_COLORS = [
        "#ef5350",   # red
        "#42a5f5",   # blue
        "#66bb6a",   # green
        "#ffa726",   # orange
        "#ab47bc",   # purple
        "#26c6da",   # cyan
        "#ffca28",   # amber
        "#ec407a",   # pink
        "#7e57c2",   # deep purple
        "#26a69a",   # teal
    ]
    STACK_ALPHA = 0.45

    def __init__(self, parent=None):
        super().__init__(parent)
        self._gauss: Optional[Tuple[float, float, float, float, float]] = None  # amp,mean,sigma,xmin,xmax
        # Optional overlay histogram (e.g., alpha peak for LMS modules).
        self._ovl_values: List[float] = []
        self._ovl_edges: List[float] = []
        self._ovl_gauss: Optional[Tuple[float, float, float, float, float]] = None
        # Stack mode: keep multiple histograms visible at once.
        self._stack_mode: bool = False
        # Which histogram to feed the stack from when both LMS and Alpha
        # are available (reference PMTs only).  Ignored when no overlay.
        self._stack_src: str = "LMS"
        self._stack_entries: List[Tuple[List[float], List[float], str]] = []
        # Title of the latest set_histogram() run, without the stack prefix.
        self._latest_title: str = ""
        self._x_lo = self._X_DEFAULT_LO
        self._x_hi = self._X_DEFAULT_HI
        self.setMinimumHeight(80)

        # Top-right "Stack" toggle.
        self._stack_btn = QPushButton("Stack: off", self)
        self._stack_btn.setFixedSize(86, 22)
        _f = QFont("Consolas", 9); _f.setBold(True)
        self._stack_btn.setFont(_f)
        self._stack_btn.setToolTip(
            "Stack mode:\n"
            "  off — clicking a run replaces the histogram (default)\n"
            "  on  — clicking different runs accumulates histograms\n"
            "        with distinct colours; fits are hidden")
        self._stack_btn.setStyleSheet(themed(OVERLAY_BUTTON_QSS))
        self._stack_btn.clicked.connect(self._toggle_stack)

        # Source selector for stack mode — only meaningful for the three
        # reference PMTs, which carry both an LMS and an Alpha histogram.
        # For other modules the choice silently falls back to LMS.
        self._src_btn = QPushButton("Src: LMS", self)
        self._src_btn.setFixedSize(86, 22)
        self._src_btn.setFont(_f)
        self._src_btn.setToolTip(
            "Stack source (reference PMTs only):\n"
            "  LMS   — stack the LMS histogram (default)\n"
            "  Alpha — stack the alpha-peak histogram\n"
            "Toggling clears the current stack to avoid mixing\n"
            "LMS and Alpha distributions in one plot.")
        self._src_btn.setStyleSheet(themed(OVERLAY_BUTTON_QSS))
        self._src_btn.clicked.connect(self._toggle_stack_src)

    def set_histogram(self, values, edges, title: str = "",
                      gauss: Optional[Tuple[float, float, float, float, float]] = None,
                      overlay_values=None, overlay_edges=None,
                      overlay_gauss: Optional[Tuple[float, float, float, float, float]] = None):
        # The single-mode payload is kept in stack mode too, so toggling
        # the source can re-seed the stack from the currently shown run.
        self._values = list(values)
        self._edges = list(edges)
        self._ovl_values = list(overlay_values) if overlay_values is not None else []
        self._ovl_edges = list(overlay_edges) if overlay_edges is not None else []
        self._ovl_gauss = overlay_gauss
        self._latest_title = title
        if self._stack_mode:
            # Fits and overlays are hidden in this mode so multiple
            # distributions stay readable.
            entry = self._stack_entry(title)
            if entry:
                self._stack_entries.append(entry)
                if len(self._stack_entries) == 1:
                    self._x_lo = self._X_DEFAULT_LO
                    self._x_hi = self._X_DEFAULT_HI
            self._title = self._stack_title(title)
            self._gauss = None
        else:
            self._title = title
            self._gauss = gauss
            self._x_lo = self._X_DEFAULT_LO
            self._x_hi = self._X_DEFAULT_HI
        self._drag_start = self._drag_cur = None
        self.update()

    def clear(self):
        self._values = []
        self._edges = []
        self._title = ""
        self._latest_title = ""
        self._gauss = None
        self._ovl_values = []
        self._ovl_edges = []
        self._ovl_gauss = None
        self._stack_entries = []
        self._x_lo = self._X_DEFAULT_LO
        self._x_hi = self._X_DEFAULT_HI
        self._drag_start = self._drag_cur = None
        self.update()

    def _stack_entry(self, title: str):
        """Stack entry (values, edges, label) of the current histogram in
        the chosen source; non-reference modules have no alpha overlay and
        fall back to LMS.  None when there is nothing to stack."""
        if (self._stack_src == "Alpha" and self._ovl_values
                and len(self._ovl_edges) >= 2):
            return (list(self._ovl_values), list(self._ovl_edges),
                    f"{title} [α]")
        if self._values and self._edges:
            return (list(self._values), list(self._edges), title)
        return None

    def _stack_title(self, latest: Optional[str] = None) -> str:
        n = len(self._stack_entries)
        title = f"Stack [{self._stack_src}]: {n} run(s)"
        return f"{title} — latest: {latest}" if n and latest is not None else title

    def _reseed_stack(self):
        """Restart the stack from the histogram currently shown."""
        entry = self._stack_entry(self._latest_title)
        self._stack_entries = [entry] if entry else []
        self._title = self._stack_title(self._latest_title)
        self._gauss = None

    def _toggle_stack(self):
        self._stack_mode = not self._stack_mode
        self._stack_btn.setText("Stack: on" if self._stack_mode else "Stack: off")
        if self._stack_mode:
            # Seed the stack with the histogram currently showing so the
            # user doesn't have to re-click the first run.
            self._reseed_stack()
        else:
            # Leaving stack mode — drop the stack; the next set_histogram
            # call will repopulate the single-hist view.
            self._stack_entries = []
        self.update()

    def _toggle_stack_src(self):
        # Switching sources restarts the stack so the two distributions are
        # never mixed in a single plot.
        self._stack_src = "Alpha" if self._stack_src == "LMS" else "LMS"
        self._src_btn.setText(f"Src: {self._stack_src}")
        if self._stack_mode:
            self._reseed_stack()
        self.update()

    def resizeEvent(self, event):
        super().resizeEvent(event)
        # Anchor the Stack button to the top-right corner of the canvas,
        # and the Src selector immediately to its left.
        stack_x = self.width() - self._stack_btn.width() - 6
        self._stack_btn.move(stack_x, 4)
        self._src_btn.move(stack_x - self._src_btn.width() - 6, 4)

    def _extend_menu(self, menu):
        if self._stack_mode and self._stack_entries:
            menu.addAction("Clear stack").triggered.connect(self._clear_stack)

    def _clear_stack(self):
        self._stack_entries = []
        self._title = self._stack_title()
        self.update()

    def _unzoom(self):
        if self._stack_mode and self._stack_entries and self._stack_entries[0][1]:
            edgs = self._stack_entries[0][1]
            self._x_lo, self._x_hi = edgs[0], edgs[-1]
            self.update()
        elif self._edges:
            super()._unzoom()
        else:
            self._x_lo, self._x_hi = self._X_DEFAULT_LO, self._X_DEFAULT_HI
            self.update()

    def _draw_stack_legend(self, p, px, py, pw):
        if not self._stack_entries:
            return
        p.setFont(QFont("Consolas", 8))
        fm = p.fontMetrics()
        swatch = 10
        gap    = 4
        line_h = max(swatch, fm.height()) + 2
        # Truncate long titles so the legend never overruns the plot.
        max_label_w = pw // 2
        rows = []
        for idx, (_v, _e, title) in enumerate(self._stack_entries):
            label = title or f"#{idx+1}"
            while fm.horizontalAdvance(label) > max_label_w and len(label) > 4:
                label = label[:-2] + "…"
            rows.append((idx, label, fm.horizontalAdvance(label)))
        if not rows:
            return
        legend_w = swatch + gap + max(w for _, _, w in rows) + 8
        legend_h = len(rows) * line_h + 6
        # Anchor below the Stack button (which is at y=4, height ~22).
        lx = px + pw - legend_w - 6
        ly = py + 30
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor(10, 14, 20, 200))
        p.drawRoundedRect(QRectF(lx, ly, legend_w, legend_h), 4, 4)
        n_colors = len(self.STACK_COLORS)
        for row_idx, (entry_idx, label, _lw) in enumerate(rows):
            sy = ly + 3 + row_idx * line_h
            col = QColor(self.STACK_COLORS[entry_idx % n_colors])
            col.setAlphaF(self.STACK_ALPHA)
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(col)
            p.drawRect(QRectF(lx + 4, sy + (line_h - swatch) // 2 - 1,
                              swatch, swatch))
            p.setPen(QColor(THEME.TEXT))
            p.drawText(QRectF(lx + 4 + swatch + gap, sy,
                              legend_w - swatch - gap - 8, line_h),
                       Qt.AlignmentFlag.AlignVCenter | Qt.AlignmentFlag.AlignLeft,
                       label)

    def paintEvent(self, event):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(0, 0, self.width(), self.height(), QColor(THEME.CANVAS))
        self._paint_title(p)

        # In stack mode treat any non-empty stack as "have data".
        have_data = (self._values and self._edges and len(self._edges) >= 2)
        have_stack = bool(self._stack_entries)
        if not have_data and not have_stack:
            self._paint_placeholder(p, "No histogram")
            p.end()
            return

        px, py, pw, ph = self._plot_rect()
        if pw < 20 or ph < 20:
            p.end()
            return

        x_lo, x_hi, _ = self._x_view()
        to_sx = self._x_map()[0]

        def _vis_max(values, edges):
            return max(self._visible_values(values, edges), default=0.0)

        # y scale from visible bins only — must fit primary, overlay, and
        # every stacked entry.
        if self._stack_mode and have_stack:
            vis_max = max((_vis_max(v, e) for v, e, _ in self._stack_entries),
                          default=0.0)
        else:
            vis_max = max(_vis_max(self._values, self._edges),
                          _vis_max(self._ovl_values, self._ovl_edges))
        y_hi = vis_max * 1.1 if vis_max > 0 else 1.0

        def to_sy(v):
            return py + ph * (1.0 - v / y_hi)

        self._paint_grid(p)

        if self._stack_mode and have_stack:
            # Each stacked run gets its own colour from STACK_COLORS, with
            # a semi-transparent alpha so overlapping bars remain visible.
            n_colors = len(self.STACK_COLORS)
            for idx, (vals, edgs, _stack_title) in enumerate(self._stack_entries):
                col = QColor(self.STACK_COLORS[idx % n_colors])
                col.setAlphaF(self.STACK_ALPHA)
                self._paint_bars(p, vals, edgs, col, to_sy)
            self._draw_stack_legend(p, px, py, pw)
        else:
            # Overlay histogram (alpha peak for LMS modules) — drawn first
            # so the primary LMS bars sit visually on top.
            if self._ovl_values and self._ovl_edges and len(self._ovl_edges) >= 2:
                ovl_color = QColor(self.OVERLAY_BAR_COLOR)
                ovl_color.setAlphaF(0.55)
                self._paint_bars(p, self._ovl_values, self._ovl_edges,
                                 ovl_color, to_sy)
            # Primary histogram (LMS or generic module) — full opacity.
            self._paint_bars(p, self._values, self._edges,
                             QColor(THEME.HIGHLIGHT), to_sy)

            def _draw_gauss(g, color):
                if g is None:
                    return
                amp, mean, sigma, g_xmin, g_xmax = g
                if sigma <= 0:
                    return
                draw_lo = max(g_xmin, x_lo)
                draw_hi = min(g_xmax, x_hi)
                if draw_hi <= draw_lo:
                    return
                n_pts = max(int((draw_hi - draw_lo) / (x_hi - x_lo) * pw), 2)
                pts = []
                for k in range(n_pts + 1):
                    gx = draw_lo + k / n_pts * (draw_hi - draw_lo)
                    gy = amp * math.exp(-0.5 * ((gx - mean) / sigma) ** 2)
                    pts.append(QPointF(to_sx(gx), to_sy(gy)))
                p.setPen(QPen(color, 2))
                for k in range(len(pts) - 1):
                    p.drawLine(pts[k], pts[k + 1])

            # LMS gaussian fit (cyan) and alpha gaussian fit (yellow).
            _draw_gauss(self._gauss,     QColor("#00bcd4"))
            _draw_gauss(self._ovl_gauss, QColor(self.OVERLAY_FIT_COLOR))

        self._paint_drag(p)
        self._paint_axes(p, self._lin_y_labels(y_hi))
        p.end()


# ---- Irregular channels table ----

class IrregularTableWidget(QWidget):
    """Table of outlier entries: irregular gains (deviation), run-to-run
    drift, or the drift summary."""

    runClicked = pyqtSignal(int)    # emits run number when a row is clicked
    moduleClicked = pyqtSignal(str) # emits module name when a row is clicked

    # mode -> (column headers, count label, row(entry) -> (name colour,
    # values of the columns after Module))
    _MODES = {
        "irregular": (
            ["Module", "Run", "Gain", "Mean", "Std Dev", "Dev (σ)"],
            "irregular gains: {} entries",
            lambda e: (THEME.WARN,
                       [e.run_number, round(e.gain, 5), round(e.mean_gain, 5),
                        round(e.std_dev, 5), round(e.deviation_sigma, 5)])),
        "drift": (
            ["Module", "Curr Run", "Prev Run", "Gain (curr)", "Gain (prev)",
             "Δ (%)"],
            "drifted channels: {} entries",
            lambda e: (THEME.DANGER if e.rel_change < 0 else THEME.SUCCESS,
                       [e.run_number, e.prev_run_number,
                        round(e.gain_current, 5), round(e.gain_prev, 5),
                        round(e.rel_change * 100, 3)])),
        "summary": (
            ["Module", "Drift Counts", "Max |Δ%|", "Worst Run", "Prev Run"],
            "problematic channels: {}",
            lambda e: (THEME.DANGER,
                       [e.drift_count, round(e.max_rel_change * 100, 3),
                        e.max_run, e.max_prev_run])),
    }

    def __init__(self, parent=None):
        super().__init__(parent)
        self._entries: list = []
        self._mode: str = "irregular"
        self._build_ui()

    def _build_ui(self):
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(4)

        # filter bar
        fbar = QHBoxLayout()
        fbar.setSpacing(6)

        fbar.addWidget(_slabel("Search:"))

        self._search = QLineEdit()
        self._search.setPlaceholderText("module name...")
        self._search.setFixedWidth(120)
        self._search.setFont(QFont("Consolas", 10))
        self._search.setStyleSheet(themed(_EDIT_QSS))
        self._search.textChanged.connect(self._apply_filter)
        fbar.addWidget(self._search)

        fbar.addWidget(_slabel("Type:"))

        self._type_filter = QComboBox()
        self._type_filter.addItems(["All", "PbWO4", "PbGlass"])
        self._type_filter.setFixedWidth(100)
        self._type_filter.setFont(QFont("Consolas", 10))
        self._type_filter.setStyleSheet(themed(_COMBO_QSS))
        self._type_filter.currentIndexChanged.connect(
            lambda _: self._apply_filter())
        fbar.addWidget(self._type_filter)

        fbar.addStretch()

        count_lbl = QLabel("")
        count_lbl.setFont(QFont("Consolas", 10))
        count_lbl.setStyleSheet(themed("color:#8b949e;"))
        self._count_lbl = count_lbl
        fbar.addWidget(count_lbl)

        layout.addLayout(fbar)

        self._table = QTableWidget()
        headers = self._MODES[self._mode][0]
        self._table.setColumnCount(len(headers))
        self._table.setHorizontalHeaderLabels(headers)
        self._table.setFont(QFont("Consolas", 10))
        self._table.setStyleSheet(themed(
            "QTableWidget{background:#0d1117;color:#c9d1d9;"
            "gridline-color:#21262d;border:1px solid #30363d;}"
            "QTableWidget::item{padding:2px 6px;}"
            "QHeaderView::section{background:#161b22;color:#58a6ff;"
            "border:1px solid #30363d;font:bold 10pt Consolas;padding:4px;}"))
        self._table.setAlternatingRowColors(True)
        pal = self._table.palette()
        pal.setColor(QPalette.ColorRole.AlternateBase, QColor(THEME.ALT_BASE))
        self._table.setPalette(pal)
        self._table.setEditTriggers(
            QAbstractItemView.EditTrigger.NoEditTriggers)
        self._table.setSelectionBehavior(
            QAbstractItemView.SelectionBehavior.SelectRows)
        self._table.setSortingEnabled(False)
        self._table.cellClicked.connect(self._on_cell_clicked)
        self._table.horizontalHeader().setStretchLastSection(True)
        self._table.verticalHeader().setVisible(False)
        self._table.verticalHeader().setDefaultSectionSize(24)

        hdr = self._table.horizontalHeader()
        hdr.setSectionResizeMode(QHeaderView.ResizeMode.Interactive)
        self._table.setColumnWidth(0, 90)
        self._table.setColumnWidth(1, 75)

        layout.addWidget(self._table)

    def set_data(self, entries: List[IrregularEntry]):
        self._set_mode_data("irregular", entries)

    def set_drift_data(self, entries: List[DriftEntry]):
        self._set_mode_data("drift", entries)

    def set_summary_data(self, entries: List[SummaryEntry]):
        self._set_mode_data("summary", entries)

    def _set_mode_data(self, mode: str, entries: list):
        self._entries = entries
        if self._mode != mode:
            self._mode = mode
            headers = self._MODES[mode][0]
            self._table.setColumnCount(len(headers))
            self._table.setHorizontalHeaderLabels(headers)
        self._apply_filter()

    def _apply_filter(self):
        search = self._search.text().strip().upper()
        type_sel = self._type_filter.currentText()

        filtered = []
        for e in self._entries:
            if search and search not in e.name.upper():
                continue
            if type_sel == "PbWO4" and e.mod_type != "PbWO4":
                continue
            if type_sel == "PbGlass" and e.mod_type != "PbGlass":
                continue
            filtered.append(e)

        _headers, count_fmt, row_of = self._MODES[self._mode]
        self._count_lbl.setText(count_fmt.format(len(filtered)))
        self._table.setRowCount(len(filtered))
        for row, e in enumerate(filtered):
            colour, values = row_of(e)
            item = QTableWidgetItem(e.name)
            item.setForeground(QColor(colour))
            self._table.setItem(row, 0, item)
            for col, val in enumerate(values, 1):
                item_v = QTableWidgetItem()
                item_v.setData(Qt.ItemDataRole.DisplayRole, val)
                self._table.setItem(row, col, item_v)

    def select_module(self, name: str):
        """Highlight and scroll to the first row matching name, or clear selection."""
        for row in range(self._table.rowCount()):
            item = self._table.item(row, 0)
            if item is not None and item.text() == name:
                self._table.selectRow(row)
                self._table.scrollToItem(item)
                return
        self._table.clearSelection()

    def _on_cell_clicked(self, row: int, _col: int):
        name_item = self._table.item(row, 0)
        if name_item is not None:
            self.moduleClicked.emit(name_item.text())
        # Run / Curr Run column; in summary mode it holds Drift Counts.
        item = self._table.item(row, 1)
        if item is not None:
            run_num = item.data(Qt.ItemDataRole.DisplayRole)
            if isinstance(run_num, int):
                self.runClicked.emit(run_num)


# ---- Process dialogs (Analyze Data / Get Data / Do It All) ----

_SCRIPT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             "shell", "run_gain_monitor.sh")
# OUTPUTDIR default of run_gain_monitor.sh.
_DEFAULT_OUTPUT_DIR = "/home/clasrun/prad2_daq/gain_monitoring/gain_monitor_output"

_DIALOG_QSS = (
    "QDialog{background:#0d1117;color:#c9d1d9;}"
    "QLabel{color:#c9d1d9;font-family:Consolas;font-size:10pt;}"
    "QLineEdit{background:#161b22;color:#c9d1d9;"
    "border:1px solid #30363d;border-radius:3px;padding:2px 6px;"
    "font-family:Consolas;font-size:10pt;}"
    "QTextEdit{background:#0a0e14;color:#c9d1d9;"
    "border:1px solid #30363d;font-family:Consolas;font-size:9pt;}"
    "QPushButton{background:#21262d;color:#c9d1d9;"
    "border:1px solid #30363d;padding:4px 12px;"
    "font:bold 10pt Consolas;border-radius:3px;}"
    "QPushButton:hover{background:#30363d;}"
    "QPushButton:disabled{color:#555;}")
_PRIMARY_BTN_QSS = (
    "QPushButton{background:#1f6feb;color:white;border:1px solid #388bfd;"
    "padding:4px 16px;font:bold 10pt Consolas;border-radius:3px;}"
    "QPushButton:hover{background:#388bfd;}"
    "QPushButton:disabled{background:#21262d;color:#555;border-color:#30363d;}")


class _ProcessDialog(QDialog):
    """Popup with an input form, Run/Stop/Close buttons and a console that
    streams the output of a bash process.  Subclasses add the form rows in
    _build_form() and start the process with _start() from _on_run()."""

    def __init__(self, parent, title: str, size: Tuple[int, int],
                 run_label: str):
        super().__init__(parent)
        self.setWindowTitle(title)
        self.resize(*size)
        self.setStyleSheet(themed(_DIALOG_QSS))

        self._process = QProcess(self)
        self._process.readyReadStandardOutput.connect(self._on_stdout)
        self._process.readyReadStandardError.connect(self._on_stderr)
        self._process.finished.connect(self._on_finished)
        self._process.errorOccurred.connect(self._on_error)

        root = QVBoxLayout(self)
        root.setSpacing(8)

        # inputs
        form = QFormLayout()
        form.setSpacing(6)
        self._build_form(form)
        root.addLayout(form)

        # buttons
        btn_row = QHBoxLayout()
        self._run_btn = QPushButton(run_label)
        self._run_btn.setStyleSheet(themed(_PRIMARY_BTN_QSS))
        self._run_btn.clicked.connect(self._on_run)
        self._stop_btn = QPushButton("Stop")
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._on_stop)
        close_btn = QPushButton("Close")
        close_btn.clicked.connect(self.close)
        btn_row.addWidget(self._run_btn)
        btn_row.addWidget(self._stop_btn)
        btn_row.addStretch()
        btn_row.addWidget(close_btn)
        root.addLayout(btn_row)

        self._console = QTextEdit()
        self._console.setReadOnly(True)
        self._console.document().setMaximumBlockCount(5000)
        root.addWidget(self._console, stretch=1)

    def _build_form(self, form: QFormLayout):
        """Hook: add the input rows."""

    def _on_run(self):
        """Hook: check the inputs and _start() the process."""

    # -- form rows --

    def _dir_row(self, caption: str, default: str = "",
                 placeholder: str = "") -> Tuple[QHBoxLayout, QLineEdit]:
        """Directory field with a Browse… button that opens a directory
        dialog titled ``caption``."""
        row = QHBoxLayout()
        edit = QLineEdit(default)
        if placeholder:
            edit.setPlaceholderText(placeholder)
        btn = QPushButton("Browse…")
        btn.setFixedWidth(80)

        def browse():
            d = QFileDialog.getExistingDirectory(self, caption)
            if d:
                edit.setText(d)

        btn.clicked.connect(browse)
        row.addWidget(edit)
        row.addWidget(btn)
        return row, edit

    def _file_range_row(self) -> QHBoxLayout:
        row = QHBoxLayout()
        self._start_edit = QLineEdit("0")
        self._start_edit.setFixedWidth(80)
        self._end_edit = QLineEdit("99")
        self._end_edit.setFixedWidth(80)
        row.addWidget(self._start_edit)
        row.addWidget(QLabel(" to "))
        row.addWidget(self._end_edit)
        row.addStretch()
        return row

    def _file_range(self) -> Optional[Tuple[int, int]]:
        """(first, last) evio file number of the range row (empty fields
        mean 0 and 9999); None after reporting invalid input."""
        start_text = self._start_edit.text().strip() or "0"
        end_text   = self._end_edit.text().strip()   or "9999"
        if not start_text.isdigit() or not end_text.isdigit():
            self._append("<span style='color:#f85149'>File number range must be integers.</span>")
            return None
        f_start = int(start_text)
        f_end   = int(end_text)
        if f_end < f_start:
            self._append("<span style='color:#f85149'>End file number must be ≥ start.</span>")
            return None
        return f_start, f_end

    def _host_row(self) -> QHBoxLayout:
        row = QHBoxLayout()
        self._host_edit    = QLineEdit(REMOTE_HOST)
        self._rembase_edit = QLineEdit(REMOTE_DATA_BASE)
        row.addWidget(self._host_edit)
        row.addWidget(QLabel("  base:"))
        row.addWidget(self._rembase_edit)
        return row

    # -- evio copy checks --

    def _confirm_existing(self, run_dir: str, files: List[str],
                          tail: str) -> bool:
        """List the requested ``files`` already in ``run_dir``; False if
        the user cancels."""
        if not files:
            return True
        box = QMessageBox(self)
        box.setWindowTitle("Files Already Present")
        box.setText(
            f"{len(files)} file(s) in the requested range already exist "
            f"in {run_dir}.\n{tail}")
        box.setDetailedText("\n".join(files))
        box.setStandardButtons(
            QMessageBox.StandardButton.Ok | QMessageBox.StandardButton.Cancel)
        box.setDefaultButton(QMessageBox.StandardButton.Ok)
        return box.exec() != QMessageBox.StandardButton.Cancel

    def _warn_insufficient(self, needed: int, free: int, hint: str = ""):
        box = QMessageBox(self)
        box.setWindowTitle("Insufficient Disk Space")
        box.setIcon(QMessageBox.Icon.Critical)
        box.setText(
            "Not enough disk space to copy the requested files.\n\n"
            f"  Required : {fmt_bytes(needed)}{hint}\n"
            f"  Available: {fmt_bytes(free)}\n\n"
            "Free up space and try again.")
        box.setStandardButtons(QMessageBox.StandardButton.Ok)
        box.exec()

    # -- process --

    def _start(self, args: List[str]):
        self._run_btn.setEnabled(False)
        self._stop_btn.setEnabled(True)
        self._process.start("bash", args)

    def _on_stop(self):
        self._process.kill()

    def _on_stdout(self):
        data = self._process.readAllStandardOutput().data().decode(errors="replace")
        self._append(html.escape(data, quote=False).replace("\n", "<br>"))

    def _on_stderr(self):
        data = self._process.readAllStandardError().data().decode(errors="replace")
        data = html.escape(data, quote=False)
        self._append(f"<span style='color:#f85149'>{data.replace(chr(10), '<br>')}</span>")

    def _on_finished(self, exit_code, exit_status):
        self._run_btn.setEnabled(True)
        self._stop_btn.setEnabled(False)
        color = "#3fb950" if exit_code == 0 else "#f85149"
        self._append(f"<span style='color:{color}'>[Process finished with exit code {exit_code}]</span>")

    def _on_error(self, error):
        # finished is not emitted when bash cannot be started at all.
        if error == QProcess.ProcessError.FailedToStart:
            self._run_btn.setEnabled(True)
            self._stop_btn.setEnabled(False)
            self._append("<span style='color:#f85149'>[Failed to start: "
                         f"{html.escape(self._process.errorString())}]</span>")

    def _append(self, markup: str):
        self._console.moveCursor(self._console.textCursor().MoveOperation.End)
        self._console.insertHtml(themed(markup))
        self._console.moveCursor(self._console.textCursor().MoveOperation.End)

    def closeEvent(self, event):
        if self._process.state() != QProcess.ProcessState.NotRunning:
            self._process.kill()
            self._process.waitForFinished(2000)
        super().closeEvent(event)


class AnalyzeDialog(_ProcessDialog):
    """Popup that runs run_gain_monitor.sh and streams its output."""

    def __init__(self, parent=None):
        super().__init__(parent, "Analyze Data", (700, 500), "Run")

    def _build_form(self, form):
        self._run_edit = QLineEdit()
        self._run_edit.setPlaceholderText("e.g. 023735")
        self._cpu_edit = QLineEdit("25")
        # Empty directories leave the defaults to run_gain_monitor.sh.
        in_row, self._indir_edit = self._dir_row(
            "Select Input Directory", placeholder=LOCAL_DATA_BASE)
        out_row, self._outdir_edit = self._dir_row(
            "Select Output Directory", placeholder=_DEFAULT_OUTPUT_DIR)

        form.addRow("Run number:", self._run_edit)
        form.addRow("Number of CPUs:", self._cpu_edit)
        form.addRow("Input directory:", in_row)
        form.addRow("Output directory:", out_row)

    def _on_run(self):
        run_num = self._run_edit.text().strip()
        n_cpu = self._cpu_edit.text().strip()
        if not run_num:
            self._append("<span style='color:#f85149'>Please enter a run number.</span>")
            return
        if not n_cpu.isdigit() or int(n_cpu) < 1:
            self._append("<span style='color:#f85149'>Number of CPUs must be a positive integer.</span>")
            return
        if not os.path.exists(_SCRIPT_PATH):
            self._append(f"<span style='color:#f85149'>Script not found: {_SCRIPT_PATH}</span>")
            return

        env = QProcessEnvironment.systemEnvironment()
        indir = self._indir_edit.text().strip()
        outdir = self._outdir_edit.text().strip()
        if indir:
            env.insert("INPUTDIR", indir)
        if outdir:
            env.insert("OUTPUTDIR", outdir)
        self._process.setProcessEnvironment(env)
        self._process.setWorkingDirectory(os.path.dirname(_SCRIPT_PATH))

        self._console.clear()
        extra = ""
        if indir:
            extra += f" INPUTDIR={indir}"
        if outdir:
            extra += f" OUTPUTDIR={outdir}"
        self._append(f"<span style='color:#8b949e'>${extra} {_SCRIPT_PATH} {run_num} {n_cpu}</span><br>")
        self._start([_SCRIPT_PATH, run_num, n_cpu])


class GetDataDialog(_ProcessDialog):
    """Popup that scps evio files from the DAQ machine for a given run."""

    def __init__(self, parent=None):
        super().__init__(parent, "Get Data", (700, 520), "Get Data")

    def _build_form(self, form):
        self._run_edit = QLineEdit()
        self._run_edit.setPlaceholderText("e.g. 023739")
        form.addRow("Run number:", self._run_edit)
        form.addRow("File number range:", self._file_range_row())
        in_row, self._localbase_edit = self._dir_row(
            "Select Local Data Directory", LOCAL_DATA_BASE)
        form.addRow("Local data directory:", in_row)
        form.addRow("Remote host:", self._host_row())

    def _on_run(self):
        run_num = self._run_edit.text().strip()
        if not run_num:
            self._append("<span style='color:#f85149'>Please enter a run number.</span>")
            return
        file_range = self._file_range()
        if file_range is None:
            return
        f_start, f_end = file_range

        local_base  = self._localbase_edit.text().strip() or LOCAL_DATA_BASE
        remote_host = self._host_edit.text().strip() or REMOTE_HOST
        remote_base = self._rembase_edit.text().strip() or REMOTE_DATA_BASE

        local_run_dir  = f"{local_base}/prad_{run_num}"
        remote_run_dir = f"{remote_base}/prad_{run_num}"

        existing = local_evio_in_range(local_run_dir, f_start, f_end,
                                       f"prad_{run_num}.evio.*")
        if not self._confirm_existing(
                local_run_dir, existing,
                "They will be skipped; only missing files will be copied."):
            return

        # -- disk space check --
        try:
            needed, free = check_disk_space(
                remote_host, remote_run_dir, local_base, f_start, f_end,
                local_run_dir)
        except Exception as exc:
            box = QMessageBox(self)
            box.setWindowTitle("Disk Space Check Failed")
            box.setIcon(QMessageBox.Icon.Warning)
            box.setText(f"Could not check remote file sizes:\n{exc}\n\nProceed anyway?")
            box.setStandardButtons(
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No)
            box.setDefaultButton(QMessageBox.StandardButton.No)
            if box.exec() != QMessageBox.StandardButton.Yes:
                return
        else:
            if needed > free:
                self._warn_insufficient(needed, free)
                return

        self._console.clear()
        self._append(
            f"<span style='color:#8b949e'>Run {run_num}, files {f_start}–{f_end}"
            f" → {local_run_dir}</span><br>")
        if existing:
            self._append(
                f"<span style='color:#d29922'>{len(existing)} file(s) skipped "
                f"(already present).</span><br>")
        self._start(["-c", scp_bash(remote_host, remote_run_dir,
                                    local_run_dir, f_start, f_end)])


class DoItAllDialog(_ProcessDialog):
    """Runs scp then gain monitor analysis in a single sequential workflow."""

    def __init__(self, parent=None):
        super().__init__(parent, "Do It All", (750, 600), "Run")

    def _build_form(self, form):
        self._run_edit = QLineEdit()
        self._run_edit.setPlaceholderText("e.g. 023739")
        form.addRow("Run number(s):", self._run_edit)

        self._batch_check = QCheckBox(
            "Batch mode (multiple runs; raw evio dirs are deleted after each "
            "run is analysed)")
        self._batch_check.setStyleSheet(themed("color:#c9d1d9;"))
        self._batch_check.toggled.connect(self._on_batch_toggled)
        form.addRow("", self._batch_check)

        form.addRow("File number range:", self._file_range_row())

        self._cpu_edit = QLineEdit("25")
        form.addRow("Number of CPUs:", self._cpu_edit)

        in_row,  self._localbase_edit = self._dir_row(
            "Select Local Data Directory", LOCAL_DATA_BASE)
        out_row, self._outdir_edit    = self._dir_row(
            "Select Output Directory", _DEFAULT_OUTPUT_DIR)
        form.addRow("Local data directory:", in_row)
        form.addRow("Output directory:", out_row)

        form.addRow("Remote host:", self._host_row())

    def _on_batch_toggled(self, on: bool):
        if on:
            self._run_edit.setPlaceholderText(
                "e.g. 023739 023740 023741   (space- or comma-separated)")
        else:
            self._run_edit.setPlaceholderText("e.g. 023739")
    def _on_run(self):
        raw_runs   = self._run_edit.text().strip()
        n_cpu      = self._cpu_edit.text().strip()
        batch_mode = self._batch_check.isChecked()

        if not raw_runs:
            self._append("<span style='color:#f85149'>Please enter a run number.</span>")
            return

        # Parse run numbers (single in normal mode, multi in batch mode).
        if batch_mode:
            run_numbers = [r for r in re.split(r'[,\s]+', raw_runs) if r]
        else:
            run_numbers = [raw_runs]
        if not run_numbers:
            self._append("<span style='color:#f85149'>No run numbers parsed.</span>")
            return
        for r in run_numbers:
            if not r.isdigit():
                self._append(f"<span style='color:#f85149'>"
                             f"Invalid run number: {r}</span>")
                return

        if not n_cpu.isdigit() or int(n_cpu) < 1:
            self._append("<span style='color:#f85149'>Number of CPUs must be a positive integer.</span>")
            return
        file_range = self._file_range()
        if file_range is None:
            return
        f_start, f_end = file_range
        if not os.path.exists(_SCRIPT_PATH):
            self._append(f"<span style='color:#f85149'>Script not found: {_SCRIPT_PATH}</span>")
            return

        local_base  = self._localbase_edit.text().strip() or LOCAL_DATA_BASE
        remote_host = self._host_edit.text().strip()      or REMOTE_HOST
        remote_base = self._rembase_edit.text().strip()   or REMOTE_DATA_BASE
        outdir      = self._outdir_edit.text().strip()    or _DEFAULT_OUTPUT_DIR

        script_dir = os.path.dirname(_SCRIPT_PATH)

        # Existing-files check only makes sense for single-run mode; in batch
        # the bash loop already skips per-file (see "Already exists" branch).
        if not batch_mode:
            run_num = run_numbers[0]
            local_run_dir = f"{local_base}/prad_{run_num}"
            existing = local_evio_in_range(local_run_dir, f_start, f_end,
                                           f"prad_{run_num}.evio.*")
            if not self._confirm_existing(local_run_dir, existing,
                                          "They will be skipped."):
                return

        # Disk-space check: in batch mode only check the first run up front,
        # since raw files for each completed run are removed before the next
        # one starts.  The bash loop also re-checks per-run before scp using
        # a conservative ~2GB-per-file estimate (see runtime block below).
        first_run = run_numbers[0]
        first_remote_run_dir = f"{remote_base}/prad_{first_run}"
        first_local_run_dir = f"{local_base}/prad_{first_run}"
        try:
            needed, free = check_disk_space(
                remote_host, first_remote_run_dir, local_base, f_start, f_end,
                first_local_run_dir)
        except Exception as exc:
            # SSH itself failed — fall back to the ~2GB-per-file estimate so
            # the user isn't asked to "proceed without any check".
            needed = (f_end - f_start + 1) * EVIO_BYTES_PER_FILE_EST
            free = free_bytes(local_base)
            self._append(
                f"<span style='color:#d29922'>Could not reach remote host "
                f"({exc}); estimating {fmt_bytes(needed)} needed using "
                f"~2 GB per file.</span><br>")
        if needed > free:
            self._warn_insufficient(
                needed, free,
                "  (one run at a time; raw dir is wiped after each)"
                if batch_mode else "")
            return

        # Bash: a single loop covers both single-run and batch.  In batch
        # mode each iteration ends with `rm -rf` of the local run dir.
        runs_string = " ".join(run_numbers)
        cleanup_block = (
            "    echo \"\"\n"
            "    echo \"=== Run $RUN: cleanup — removing $LRD ===\"\n"
            "    rm -rf \"$LRD\"\n"
        ) if batch_mode else ""

        bash_cmd = (
            f"set -e\n"
            f"RUNS=\"{runs_string}\"\n"
            f"LOCAL_BASE={local_base}\n"
            f"REMOTE_BASE={remote_base}\n"
            f"REMOTE_HOST={remote_host}\n"
            f"OUTDIR={outdir}\n"
            f"SCRIPT_DIR={script_dir}\n"
            f"SCRIPT={_SCRIPT_PATH}\n"
            f"F_START={f_start}\n"
            f"F_END={f_end}\n"
            f"N_CPU={n_cpu}\n"
            f"BYTES_PER_EVIO={EVIO_BYTES_PER_FILE_EST}\n"
            f"\n"
            f"for RUN in $RUNS; do\n"
            f"    LRD=\"$LOCAL_BASE/prad_$RUN\"\n"
            f"    RRD=\"$REMOTE_BASE/prad_$RUN\"\n"
            f"\n"
            f"    echo \"=== Run $RUN: Step 1 — Copying evio files ===\"\n"
            f"    echo \"Local directory: $LRD\"\n"
            f"    mkdir -p \"$LRD\"\n"
            f"    ALL_FILES=$(ssh \"$REMOTE_HOST\" \"ls $RRD/\" 2>/dev/null | sort)\n"
            f"\n"
            f"    # Disk-space check: count files in the requested range that\n"
            f"    # aren't already present locally, then assume ~2 GB each.\n"
            f"    TO_COPY=0\n"
            f"    while IFS= read -r f; do\n"
            f"        NUM=$(echo \"$f\" | grep -oP '\\.evio\\.\\K[0-9]+')\n"
            f"        [ -z \"$NUM\" ] && continue\n"
            f"        N=$((10#$NUM))\n"
            f"        if [ \"$N\" -lt \"$F_START\" ] || [ \"$N\" -gt \"$F_END\" ]; then continue; fi\n"
            f"        if [ -f \"$LRD/$f\" ]; then continue; fi\n"
            f"        TO_COPY=$((TO_COPY+1))\n"
            f"    done <<< \"$ALL_FILES\"\n"
            f"    NEEDED_BYTES=$((TO_COPY * BYTES_PER_EVIO))\n"
            f"    DF_PATH=\"$LRD\"\n"
            f"    while [ -n \"$DF_PATH\" ] && [ ! -e \"$DF_PATH\" ]; do\n"
            f"        DF_PATH=$(dirname \"$DF_PATH\")\n"
            f"    done\n"
            f"    FREE_BYTES=$(df -B1 --output=avail \"$DF_PATH\" 2>/dev/null | tail -n1 | tr -d ' ')\n"
            f"    NEED_GB=$((NEEDED_BYTES / 1073741824))\n"
            f"    FREE_GB=$(( ${{FREE_BYTES:-0}} / 1073741824))\n"
            f"    echo \"Disk-space check: $TO_COPY file(s) × ~2 GB ≈ ${{NEED_GB}} GB needed, ${{FREE_GB}} GB free.\"\n"
            f"    if [ -n \"$FREE_BYTES\" ] && [ \"$FREE_BYTES\" -lt \"$NEEDED_BYTES\" ]; then\n"
            f"        echo \"ERROR: insufficient disk space for run $RUN — aborting.\" >&2\n"
            f"        exit 1\n"
            f"    fi\n"
            f"\n"
            f"    COPIED=0; ALREADY=0\n"
            f"    while IFS= read -r f; do\n"
            f"        NUM=$(echo \"$f\" | grep -oP '\\.evio\\.\\K[0-9]+')\n"
            f"        [ -z \"$NUM\" ] && continue\n"
            f"        N=$((10#$NUM))\n"
            f"        if [ \"$N\" -lt \"$F_START\" ] || [ \"$N\" -gt \"$F_END\" ]; then continue; fi\n"
            f"        if [ -f \"$LRD/$f\" ]; then\n"
            f"            echo \"  Already exists: $f (skipping)\"\n"
            f"            ALREADY=$((ALREADY+1))\n"
            f"        else\n"
            f"            echo \"  Copying $f\"\n"
            f"            scp \"$REMOTE_HOST:$RRD/$f\" \"$LRD/\"\n"
            f"            COPIED=$((COPIED+1))\n"
            f"        fi\n"
            f"    done <<< \"$ALL_FILES\"\n"
            f"    echo \"Run $RUN download done. Copied $COPIED file(s), $ALREADY already present.\"\n"
            f"\n"
            f"    echo \"\"\n"
            f"    echo \"=== Run $RUN: Step 2 — Analysis ===\"\n"
            f"    cd \"$SCRIPT_DIR\"\n"
            f"    INPUTDIR=\"$LOCAL_BASE\" OUTPUTDIR=\"$OUTDIR\" bash \"$SCRIPT\" \"$RUN\" \"$N_CPU\"\n"
            + cleanup_block +
            "done\n"
            "\n"
            "echo \"\"\n"
            "echo \"=== All runs done ===\"\n"
        )

        self._console.clear()
        if batch_mode and len(run_numbers) > 1:
            summary = f"Batch: {len(run_numbers)} runs ({runs_string})"
        else:
            summary = f"Run {run_numbers[0]}"
        cleanup_note = "  | cleanup raw evio after each" if batch_mode else ""
        self._append(
            f"<span style='color:#8b949e'>{summary} | files {f_start}–{f_end}"
            f" | {n_cpu} CPUs{cleanup_note}</span><br>")
        self._start(["-c", bash_cmd])


# ---- Main window ----

class GainMonitorWindow(QMainWindow):

    def __init__(self):
        super().__init__()
        self._all_modules: List[Module] = []
        self._mod_by_name: Dict[str, Module] = {}
        self._runs: List[RunData] = []
        self._current_run_idx: int = 0
        self._current_ref_idx: int = LMS_REF_DEFAULT
        self._auto_range = True
        self._manual_vmin = 0.9
        self._manual_vmax = 1.1
        self._log_scale = False
        self._selected_module: Optional[str] = None
        # Last clicked HyCal (non-LMS) module — drives the orange mod plots.
        # Held independently of `_selected_module` so picking a reference PMT
        # doesn't blank the mod plots.
        self._selected_hycal_module: Optional[str] = None
        self._start_run_idx: int = 0
        self._end_run_idx: int = 0
        self._thresh_g: float = 0.10
        self._thresh_w: float = 0.05
        self._view_mode: int = 2   # 0 = Gain Factor, 1 = Deviation (σ), 2 = Run-to-Run Drift, 3 = Summary
        self._pairwise_diffs: List = []
        self._pairwise_ref_idx: int = -1
        self._current_folder: str = ""
        # gain_stats() of the active runs, valid for (start, end, ref) == key
        self._gain_stats_cache: Dict[str, Tuple[float, float, int]] = {}
        self._gain_stats_key: Optional[Tuple] = None
        self._file_snapshot: Dict[int, float] = {}   # run_number -> mtime
        # (run_number, module) -> set_histogram() payload (values, edges,
        # gauss, overlay values, overlay edges, overlay gauss)
        self._hist_cache: OrderedDict = OrderedDict()
        self._auto_refresh_timer = QTimer(self)
        self._auto_refresh_timer.timeout.connect(self._auto_refresh_check)

        self._load_geometry()
        self._build_ui()
        self._apply_map_mode()

    def _load_geometry(self):
        self._all_modules = load_modules(MODULES_JSON)
        self._mod_by_name = {m.name: m for m in self._all_modules}

    @property
    def _active_runs(self) -> List[RunData]:
        return self._runs[self._start_run_idx:self._end_run_idx + 1]

    # ---- UI ----

    def _build_ui(self):
        self.setWindowTitle("HyCal Gain Monitor")
        self.resize(1800, 1000)
        apply_theme_palette(self)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # -- top bar --
        top = QHBoxLayout()
        lbl = QLabel("HYCAL GAIN MONITOR")
        lbl.setFont(QFont("Consolas", 14, QFont.Weight.Bold))
        lbl.setStyleSheet(themed("color:#58a6ff;"))
        top.addWidget(lbl)
        top.addStretch()

        self._process_btn = self._make_btn(
            "Process Folder...", THEME.ACCENT, self._on_process_folder)
        top.addWidget(self._process_btn)

        self._analyze_btn = self._make_btn(
            "Analyze Data", THEME.SUCCESS, self._on_analyze_data)
        top.addWidget(self._analyze_btn)

        self._getdata_btn = self._make_btn(
            "Get Data", THEME.WARN, self._on_get_data)
        top.addWidget(self._getdata_btn)

        self._doitall_btn = self._make_btn(
            "Do It All", "#bc8cff", self._on_do_it_all)
        top.addWidget(self._doitall_btn)

        self._refresh_btn = self._make_btn(
            "Refresh", THEME.SUCCESS, self._on_refresh)
        self._refresh_btn.setEnabled(False)
        top.addWidget(self._refresh_btn)

        self._auto_refresh_btn = QPushButton("Auto-Refresh: ON")
        self._auto_refresh_btn.setCheckable(True)
        self._auto_refresh_btn.setChecked(True)
        self._auto_refresh_btn.setEnabled(False)
        self._auto_refresh_btn.setFont(QFont("Consolas", 10))
        self._auto_refresh_btn.setFixedHeight(28)
        # Themed via f-string so the :checked state's green tint comes from
        # THEME.SUCCESS / THEME.BUTTON_HOVER and tracks dark/light correctly.
        self._auto_refresh_btn.setStyleSheet(
            f"QPushButton{{background:{THEME.PANEL};color:{THEME.TEXT_DIM};"
            f"border:1px solid {THEME.BORDER};border-radius:3px;padding:0 8px;}}"
            f"QPushButton:checked{{background:{THEME.BUTTON_HOVER};"
            f"color:{THEME.SUCCESS};border-color:{THEME.SUCCESS};}}"
            f"QPushButton:hover{{background:{THEME.BUTTON};}}")
        self._auto_refresh_btn.toggled.connect(self._on_auto_refresh_toggled)
        top.addWidget(self._auto_refresh_btn)

        top.addWidget(_slabel("every"))
        self._auto_refresh_interval = QSpinBox()
        self._auto_refresh_interval.setRange(5, 3600)
        self._auto_refresh_interval.setValue(10)
        self._auto_refresh_interval.setSuffix(" s")
        self._auto_refresh_interval.setFixedWidth(72)
        self._auto_refresh_interval.setFont(QFont("Consolas", 10))
        self._auto_refresh_interval.setStyleSheet(themed(
            "QSpinBox{background:#161b22;color:#c9d1d9;border:1px solid #30363d;"
            "border-radius:3px;padding:2px 4px;}"
            "QSpinBox::up-button,QSpinBox::down-button{width:16px;}"))
        self._auto_refresh_interval.valueChanged.connect(self._on_refresh_interval_changed)
        top.addWidget(self._auto_refresh_interval)

        self._summary_btn = QPushButton("Summary table")
        self._summary_btn.setCheckable(True)
        self._summary_btn.setChecked(False)
        self._summary_btn.setFont(QFont("Consolas", 10))
        self._summary_btn.setFixedHeight(28)
        self._summary_btn.setStyleSheet(themed(
            "QPushButton{background:#21262d;color:#c9d1d9;"
            "border:1px solid #30363d;padding:4px 12px;border-radius:4px;}"
            "QPushButton:checked{background:#1f6feb;color:#ffffff;"
            "border:1px solid #1f6feb;}"
            "QPushButton:hover{background:#28282a;}"))
        self._summary_btn.toggled.connect(self._on_summary_toggle)
        top.addWidget(self._summary_btn)

        self._status_lbl = QLabel("No data loaded")
        self._status_lbl.setFont(QFont("Consolas", 11))
        self._status_lbl.setStyleSheet(themed("color:#8b949e;"))
        top.addWidget(self._status_lbl)
        root.addLayout(top)

        # -- body splitter (horizontal: map | charts+histogram) --
        body = QSplitter(Qt.Orientation.Horizontal)
        self._body = body

        # The report table lives in a hidden top-level window that the
        # "Summary table" toolbar toggle shows.
        self._irregular_table = IrregularTableWidget()
        self._irregular_table.runClicked.connect(self._on_jump_to_run)
        self._irregular_table.moduleClicked.connect(self._on_module_clicked)

        self._summary_window = QDialog(self,
            Qt.WindowType.Window
            | Qt.WindowType.WindowMinMaxButtonsHint
            | Qt.WindowType.WindowCloseButtonHint)
        self._summary_window.setWindowTitle("Summary Table")
        self._summary_window.resize(420, 600)
        _sw_layout = QVBoxLayout(self._summary_window)
        _sw_layout.setContentsMargins(6, 6, 6, 6)
        _sw_layout.addWidget(self._irregular_table)
        self._summary_window.installEventFilter(self)

        # ---- left panel: HyCal map ----
        left = QWidget()
        left_layout = QVBoxLayout(left)
        left_layout.setContentsMargins(0, 0, 0, 0)
        left_layout.setSpacing(4)

        # controls
        ctrl = QHBoxLayout()
        ctrl.setSpacing(6)

        ctrl.addWidget(_slabel("Ref:"))
        self._ref_combo = QComboBox()
        self._ref_combo.addItems(LMS_DISPLAY)
        self._ref_combo.setCurrentIndex(LMS_REF_DEFAULT)
        self._ref_combo.setFixedWidth(90)
        self._ref_combo.setFont(QFont("Consolas", 10))
        self._ref_combo.setStyleSheet(themed(_COMBO_QSS))
        self._ref_combo.currentIndexChanged.connect(self._on_ref_changed)
        ctrl.addWidget(self._ref_combo)

        ctrl.addSpacing(10)
        ctrl.addWidget(_slabel("View:"))
        self._view_combo = QComboBox()
        self._view_combo.addItems(["Gain Factor", "Deviation (σ)", "Run-to-Run Drift", "Summary"])
        self._view_combo.setFixedWidth(150)
        self._view_combo.setFont(QFont("Consolas", 10))
        self._view_combo.setStyleSheet(themed(_COMBO_QSS))
        self._view_combo.currentIndexChanged.connect(self._on_view_mode_changed)
        ctrl.addWidget(self._view_combo)

        self._view_combo.blockSignals(True)
        self._view_combo.setCurrentIndex(2)
        self._view_combo.blockSignals(False)

        self._thresh_lbl = _slabel("G thresh:")
        ctrl.addSpacing(6)
        ctrl.addWidget(self._thresh_lbl)

        self._thresh_g_input = QLineEdit("10.0")
        self._thresh_g_input.setFixedWidth(46)
        self._thresh_g_input.setFont(QFont("Consolas", 10))
        self._thresh_g_input.setStyleSheet(themed(_EDIT_QSS))
        self._thresh_g_input.editingFinished.connect(self._on_drift_threshold_changed)
        ctrl.addWidget(self._thresh_g_input)
        self._thresh_g_pct = _slabel("%")
        ctrl.addWidget(self._thresh_g_pct)

        self._thresh_w_lbl = _slabel("  W thresh:")
        ctrl.addWidget(self._thresh_w_lbl)

        self._thresh_w_input = QLineEdit("5.0")
        self._thresh_w_input.setFixedWidth(46)
        self._thresh_w_input.setFont(QFont("Consolas", 10))
        self._thresh_w_input.setStyleSheet(themed(_EDIT_QSS))
        self._thresh_w_input.editingFinished.connect(self._on_drift_threshold_changed)
        ctrl.addWidget(self._thresh_w_input)
        self._thresh_w_pct = _slabel("%")
        ctrl.addWidget(self._thresh_w_pct)

        ctrl.addSpacing(10)
        ctrl.addWidget(_slabel("Start:"))
        self._start_combo = QComboBox()
        self._start_combo.setMinimumWidth(100)
        self._start_combo.setFont(QFont("Consolas", 10))
        self._start_combo.setStyleSheet(themed(_COMBO_QSS))
        self._start_combo.currentIndexChanged.connect(self._on_start_run_changed)
        ctrl.addWidget(self._start_combo)

        ctrl.addSpacing(6)
        ctrl.addWidget(_slabel("End:"))
        self._end_combo = QComboBox()
        self._end_combo.setMinimumWidth(100)
        self._end_combo.setFont(QFont("Consolas", 10))
        self._end_combo.setStyleSheet(themed(_COMBO_QSS))
        self._end_combo.currentIndexChanged.connect(self._on_end_run_changed)
        ctrl.addWidget(self._end_combo)

        ctrl.addSpacing(10)
        ctrl.addWidget(_slabel("Run:"))

        self._run_combo = QComboBox()
        self._run_combo.setMinimumWidth(100)
        self._run_combo.setFont(QFont("Consolas", 10))
        self._run_combo.setStyleSheet(themed(_COMBO_QSS))
        self._run_combo.currentIndexChanged.connect(self._on_run_changed)
        ctrl.addWidget(self._run_combo)

        self._prev_btn = self._make_btn("<", THEME.TEXT, self._on_prev_run)
        self._prev_btn.setFixedWidth(30)
        ctrl.addWidget(self._prev_btn)

        self._next_btn = self._make_btn(">", THEME.TEXT, self._on_next_run)
        self._next_btn.setFixedWidth(30)
        ctrl.addWidget(self._next_btn)

        ctrl.addStretch()
        left_layout.addLayout(ctrl)

        # range controls
        rng = QHBoxLayout()
        rng.setSpacing(6)

        rng.addWidget(_slabel("Min:"))
        self._range_min = QLineEdit("0.9")
        self._range_min.setFixedWidth(70)
        self._range_min.setFont(QFont("Consolas", 10))
        self._range_min.setStyleSheet(themed(_EDIT_QSS))
        self._range_min.returnPressed.connect(self._on_apply_range)
        rng.addWidget(self._range_min)

        rng.addWidget(_slabel("Max:"))
        self._range_max = QLineEdit("1.1")
        self._range_max.setFixedWidth(70)
        self._range_max.setFont(QFont("Consolas", 10))
        self._range_max.setStyleSheet(themed(_EDIT_QSS))
        self._range_max.returnPressed.connect(self._on_apply_range)
        rng.addWidget(self._range_max)

        self._apply_btn = QPushButton("Apply")
        self._apply_btn.setFixedWidth(55)
        self._apply_btn.clicked.connect(self._on_apply_range)
        rng.addWidget(self._apply_btn)

        self._log_btn = QPushButton("Log")
        self._log_btn.setFixedWidth(45)
        self._log_btn.setCheckable(True)
        self._log_btn.clicked.connect(self._on_log_toggled)
        rng.addWidget(self._log_btn)

        self._auto_btn = QPushButton("Auto")
        self._auto_btn.setFixedWidth(55)
        self._auto_btn.setCheckable(True)
        self._auto_btn.setChecked(True)
        self._auto_btn.clicked.connect(self._on_auto_range)
        rng.addWidget(self._auto_btn)

        self._apply_btn.setStyleSheet(themed(_BTN_QSS))
        self._log_btn.setStyleSheet(themed(_TOGGLE_QSS))
        self._auto_btn.setStyleSheet(themed(_TOGGLE_QSS))

        rng.addStretch()
        left_layout.addLayout(rng)

        self._map = HyCalGainMapWidget()
        self._map.set_modules(self._all_modules)
        self._map.set_palette(_DEFAULT_PALETTE)
        self._map.moduleHovered.connect(self._on_module_hovered)
        self._map.moduleClicked.connect(self._on_module_clicked)
        self._map.paletteClicked.connect(self._on_cycle_palette)
        self._map.rangeEdited.connect(self._on_map_range_edited)
        left_layout.addWidget(self._map, stretch=1)

        self._info = QLabel("Hover over a module for details")
        self._info.setFont(QFont("Consolas", 10))
        self._info.setStyleSheet(themed(
            "QLabel{background:#161b22;color:#c9d1d9;padding:4px 8px;"
            "border:1px solid #30363d;border-radius:4px;}"))
        self._info.setFixedHeight(26)
        left_layout.addWidget(self._info)

        body.addWidget(left)

        # ---- right panel (vertical splitter: charts top, histogram bottom) ----
        right = QSplitter(Qt.Orientation.Vertical)

        charts = QWidget()
        charts_layout = QVBoxLayout(charts)
        charts_layout.setContentsMargins(0, 0, 0, 0)
        charts_layout.setSpacing(2)

        self._charts: List[LMSLineChartWidget] = []
        for kind, label in CHART_PLOTS:
            chart = LMSLineChartWidget()
            chart.set_data([], [], [], label)
            # The two module-driven plots (LMS peak and gain factor) are
            # tinted orange so they are visually distinct from the three
            # blue reference plots above.
            if kind in ("mod_lms", "mod_gain"):
                chart.set_series_color(QColor(THEME.HIGHLIGHT))
            chart.runClicked.connect(
                lambda rn, k=kind: self._on_chart_run_clicked(k, rn))
            chart.pointDeleteRequested.connect(
                lambda rn: self._remove_run(rn, backup=False))
            chart.pointBackupRequested.connect(
                lambda rn: self._remove_run(rn, backup=True))
            self._charts.append(chart)
            charts_layout.addWidget(chart)

        right.addWidget(charts)
        self._hist_widget = RootHistWidget()
        right.addWidget(self._hist_widget)
        # 5 stacked plots up top still need most of the vertical room, but
        # the histogram below needs enough height to be readable.
        right.setStretchFactor(0, 7)
        right.setStretchFactor(1, 3)
        self._right_splitter = right

        body.addWidget(right)
        body.setStretchFactor(0, 2)  # HyCal map
        body.setStretchFactor(1, 5)  # charts + histogram
        QTimer.singleShot(0, lambda: self._body.setSizes([500, 1100]))
        QTimer.singleShot(0, lambda: self._right_splitter.setSizes([700, 300]))

        root.addWidget(body, stretch=1)

    # ---- helpers ----

    def _make_btn(self, text, fg, slot):
        btn = QPushButton(text)
        btn.setStyleSheet(themed(
            f"QPushButton{{background:#21262d;color:{fg};"
            f"border:1px solid #30363d;padding:4px 12px;"
            f"font:bold 11px Consolas;border-radius:3px;}}"
            f"QPushButton:hover{{background:#30363d;}}"
            f"QPushButton:disabled{{color:#555;}}"))
        btn.clicked.connect(slot)
        return btn

    # ---- keyboard navigation ----

    def keyPressEvent(self, event):
        if event.key() == Qt.Key.Key_Left:
            self._on_prev_run()
        elif event.key() == Qt.Key.Key_Right:
            self._on_next_run()
        else:
            super().keyPressEvent(event)

    # ---- slots ----

    def _on_refresh(self):
        if self._current_folder:
            self._load_folder(self._current_folder)

    def _on_summary_toggle(self, checked: bool):
        if checked:
            self._summary_window.show()
            self._summary_window.raise_()
            self._summary_window.activateWindow()
        else:
            self._summary_window.hide()

    def eventFilter(self, obj, event):
        # Sync the toolbar toggle when the user closes the floating
        # summary window via its window-frame X button.
        if (obj is getattr(self, "_summary_window", None)
                and event.type() == event.Type.Close):
            if self._summary_btn.isChecked():
                self._summary_btn.blockSignals(True)
                self._summary_btn.setChecked(False)
                self._summary_btn.blockSignals(False)
        return super().eventFilter(obj, event)

    def _on_auto_refresh_toggled(self, checked: bool):
        self._auto_refresh_btn.setText("Auto-Refresh: ON" if checked else "Auto-Refresh: OFF")
        if checked:
            interval_ms = self._auto_refresh_interval.value() * 1000
            self._auto_refresh_timer.start(interval_ms)
        else:
            self._auto_refresh_timer.stop()

    def _on_refresh_interval_changed(self, value: int):
        if self._auto_refresh_timer.isActive():
            self._auto_refresh_timer.start(value * 1000)

    @staticmethod
    def _take_file_snapshot(folder: str) -> Dict[int, float]:
        """Return {run_number: mtime} for all LMS dat files in folder."""
        snapshot = {}
        for path in Path(folder).glob("prad_*_LMS.dat"):
            m = FILE_PATTERN.match(path.name)
            if m:
                snapshot[int(m.group(1))] = path.stat().st_mtime
        return snapshot

    def _auto_refresh_check(self):
        if not self._current_folder:
            return
        new_snapshot = self._take_file_snapshot(self._current_folder)
        if new_snapshot != self._file_snapshot:
            self._smart_refresh(new_snapshot)

    def _smart_refresh(self, new_snapshot: Dict[int, float]):
        """Reload data while preserving the current run selection and range."""
        current_run_number = (self._runs[self._current_run_idx].run_number
                              if self._runs else None)
        start_run_number = (self._runs[self._start_run_idx].run_number
                            if self._runs else None)
        end_run_number = (self._runs[self._end_run_idx].run_number
                          if self._runs else None)
        old_last_run_number = (self._runs[-1].run_number if self._runs else None)

        # Only parse files that are new or whose mtime changed — avoids re-reading
        # the entire folder on every auto-refresh tick when the dataset is large.
        changed_runs = {rnum for rnum, mtime in new_snapshot.items()
                        if self._file_snapshot.get(rnum) != mtime}
        if not changed_runs:
            return
        runs_by_number = {rd.run_number: rd for rd in self._runs}
        for rnum in changed_runs:
            path = os.path.join(self._current_folder,
                                f"prad_{rnum:06d}_LMS.dat")
            rd = parse_dat_file(path)
            if rd is not None:
                runs_by_number[rnum] = rd
        # Remove runs whose files have disappeared from the snapshot
        for rnum in list(runs_by_number):
            if rnum not in new_snapshot:
                del runs_by_number[rnum]
        new_runs = sorted(runs_by_number.values(), key=lambda r: r.run_number)
        if not new_runs:
            return
        self._runs = new_runs
        self._file_snapshot = new_snapshot
        self._gain_stats_key = None
        self._hist_cache.clear()

        run_numbers = [rd.run_number for rd in self._runs]
        last = len(self._runs) - 1

        if start_run_number in run_numbers:
            self._start_run_idx = run_numbers.index(start_run_number)
        else:
            self._start_run_idx = 0

        # restore end index — extend to newest run if it was already at the end before refresh
        if end_run_number == old_last_run_number or end_run_number not in run_numbers:
            self._end_run_idx = last
        else:
            self._end_run_idx = run_numbers.index(end_run_number)
        self._refill_range_combos()

        new_runs_added = run_numbers[-1] != old_last_run_number
        if new_runs_added:
            self._current_run_idx = last
        elif current_run_number in run_numbers:
            self._current_run_idx = run_numbers.index(current_run_number)
        else:
            self._current_run_idx = self._end_run_idx

        self._populate_run_combo()
        self._recompute_pairwise_diffs()
        self._status_lbl.setText(f"{len(self._runs)} runs loaded  [auto-refreshed]")
        self._status_lbl.setStyleSheet(themed("color:#3fb950;"))
        self._update_all_views()

    def _on_process_folder(self):
        folder = QFileDialog.getExistingDirectory(
            self, "Select Data Folder")
        if not folder:
            return
        self._load_folder(folder)

    def _on_analyze_data(self):
        dlg = AnalyzeDialog(self)
        dlg.exec()

    def _on_get_data(self):
        dlg = GetDataDialog(self)
        dlg.exec()

    def _on_do_it_all(self):
        dlg = DoItAllDialog(self)
        dlg.exec()

    def _load_folder(self, folder: str):
        self._status_lbl.setText("Loading...")
        self._status_lbl.setStyleSheet(themed("color:#d29922;"))
        QApplication.processEvents()

        self._runs = load_all_runs(folder)
        self._current_folder = folder
        if not self._runs:
            self._status_lbl.setText("No data files found")
            self._status_lbl.setStyleSheet(themed("color:#f85149;"))
            return

        last = len(self._runs) - 1
        self._start_run_idx = 0
        self._end_run_idx = last
        self._current_run_idx = last
        self._refill_range_combos()
        self._selected_module = None
        self._selected_hycal_module = None
        self._map.set_selected(None)
        self._gain_stats_key = None

        self._populate_run_combo()
        self._recompute_pairwise_diffs()
        self._status_lbl.setText(f"{len(self._runs)} runs loaded")
        self._status_lbl.setStyleSheet(themed("color:#3fb950;"))
        self._refresh_btn.setEnabled(True)
        self._auto_refresh_btn.setEnabled(True)
        self._file_snapshot = self._take_file_snapshot(folder)
        if self._auto_refresh_btn.isChecked():
            self._auto_refresh_timer.start(self._auto_refresh_interval.value() * 1000)

        self._update_all_views()

    @staticmethod
    def _set_combo(combo: QComboBox, idx: int,
                   labels: Optional[List[str]] = None):
        """Select ``idx`` (after refilling with ``labels``) without
        emitting signals."""
        combo.blockSignals(True)
        if labels is not None:
            combo.clear()
            combo.addItems(labels)
        combo.setCurrentIndex(idx)
        combo.blockSignals(False)

    def _refill_range_combos(self):
        labels = [str(rd.run_number) for rd in self._runs]
        self._set_combo(self._start_combo, self._start_run_idx, labels)
        self._set_combo(self._end_combo, self._end_run_idx, labels)

    def _populate_run_combo(self):
        labels = [str(rd.run_number) for rd in self._active_runs]
        combo_idx = max(0, self._current_run_idx - self._start_run_idx)
        self._set_combo(self._run_combo, min(combo_idx, len(labels) - 1), labels)

    def _set_run_range(self, start: int, end: int):
        self._start_run_idx, self._end_run_idx = start, end
        self._set_combo(self._start_combo, start)
        self._set_combo(self._end_combo, end)
        self._current_run_idx = min(max(self._current_run_idx, start), end)
        self._populate_run_combo()
        self._update_all_views()

    def _on_start_run_changed(self, index: int):
        if 0 <= index < len(self._runs):
            self._set_run_range(index, max(index, self._end_run_idx))

    def _on_end_run_changed(self, index: int):
        if 0 <= index < len(self._runs):
            self._set_run_range(min(self._start_run_idx, index), index)

    def _on_run_changed(self, index: int):
        if index < 0 or index >= len(self._active_runs):
            return
        if self._view_mode == 3:
            return
        self._current_run_idx = self._start_run_idx + index
        self._update_geo_view()
        if self._view_mode == 2:
            self._update_irregular_table()
        curr_run_number = self._runs[self._current_run_idx].run_number
        for chart in self._charts:
            chart.set_current_run(curr_run_number)
        self._load_lms_hist(self._selected_module)

    def _on_ref_changed(self, index: int):
        self._current_ref_idx = index
        self._update_all_views()

    def _on_prev_run(self):
        combo_idx = self._current_run_idx - self._start_run_idx
        if combo_idx > 0:
            self._run_combo.setCurrentIndex(combo_idx - 1)

    def _on_next_run(self):
        combo_idx = self._current_run_idx - self._start_run_idx
        if combo_idx < len(self._active_runs) - 1:
            self._run_combo.setCurrentIndex(combo_idx + 1)

    def _on_apply_range(self):
        try:
            vmin = float(self._range_min.text())
            vmax = float(self._range_max.text())
        except ValueError:
            return
        if vmin >= vmax:
            return
        self._auto_range = False
        self._auto_btn.setChecked(False)
        self._manual_vmin = vmin
        self._manual_vmax = vmax
        self._update_geo_view()

    def _on_map_range_edited(self, vmin: float, vmax: float):
        # An inline colour-bar edit sets the manual range, like Min/Max.
        self._range_min.setText(f"{vmin:g}")
        self._range_max.setText(f"{vmax:g}")
        self._on_apply_range()

    def _on_log_toggled(self):
        self._log_scale = self._log_btn.isChecked()
        self._map.set_log_scale(self._log_scale)

    def _on_auto_range(self):
        self._auto_range = self._auto_btn.isChecked()
        self._update_geo_view()

    def _on_jump_to_run(self, run_number: int):
        for i, rd in enumerate(self._active_runs):
            if rd.run_number == run_number:
                self._run_combo.setCurrentIndex(i)
                return

    def _on_chart_run_clicked(self, kind: str, run_number: int):
        # Clicking a point on one of the three reference-PMT plots should
        # switch the histogram (and selected module) to the current ref PMT
        # — otherwise the histogram would keep showing whatever HyCal
        # module was last picked on the map.
        if kind in ("ratio", "ref_lms", "ref_alpha"):
            ref_name = LMS_NAMES[self._current_ref_idx]
            if self._selected_module != ref_name:
                self._on_module_clicked(ref_name)
        if self._view_mode != 3:
            self._on_jump_to_run(run_number)
        else:
            self._load_lms_hist(self._selected_module, run_number)

    def _run_file_paths(self, run_number: int) -> List[str]:
        folder = self._current_folder
        return [
            os.path.join(folder, f"prad_{run_number:06d}_LMS.dat"),
            os.path.join(folder, f"prad_{run_number:06d}_LMS.root"),
            os.path.join(folder, f"prad_{run_number:06d}_LMS_fitted.root"),
        ]

    def _drop_run_from_state(self, run_idx: int, run_number: int):
        """Remove a run from in-memory state and refresh combos / views.
        Called by _remove_run once the on-disk action is done."""
        del self._runs[run_idx]
        self._file_snapshot.pop(run_number, None)
        for key in [k for k in self._hist_cache if k[0] == run_number]:
            self._hist_cache.pop(key, None)

        n = len(self._runs)

        def _shift(idx: int) -> int:
            return idx - 1 if idx > run_idx else idx

        self._start_run_idx = max(0, min(_shift(self._start_run_idx), n - 1))
        self._end_run_idx = max(0, min(_shift(self._end_run_idx), n - 1))
        self._current_run_idx = max(0, min(_shift(self._current_run_idx), n - 1))
        if self._start_run_idx > self._end_run_idx:
            self._start_run_idx, self._end_run_idx = (
                self._end_run_idx, self._start_run_idx)
        self._current_run_idx = max(self._start_run_idx,
                                    min(self._current_run_idx,
                                        self._end_run_idx))

        self._refill_range_combos()
        self._populate_run_combo()

        self._recompute_pairwise_diffs()
        self._gain_stats_key = None
        self._status_lbl.setText(f"{len(self._runs)} runs loaded")
        self._update_all_views()

    def _remove_run(self, run_number: int, backup: bool):
        """Delete a run's .dat/.root files on disk, or with ``backup`` move
        them into <folder>/backup/, and drop the run from memory."""
        run_idx = next((i for i, rd in enumerate(self._runs)
                        if rd.run_number == run_number), -1)
        if run_idx < 0:
            return
        if len(self._runs) <= 1:
            if backup:
                QMessageBox.warning(self, "Cannot backup",
                                    "Refusing to remove the last remaining run "
                                    "from the view.")
            else:
                QMessageBox.warning(self, "Cannot delete",
                                    "Refusing to delete the last remaining run.")
            return

        backup_dir = os.path.join(self._current_folder, "backup")
        paths = self._run_file_paths(run_number)
        existing = [p for p in paths if os.path.exists(p)]
        details = ("\n  ".join(existing) if existing
                   else "(no on-disk files found — only in-memory record)")

        box = QMessageBox(self)
        if backup:
            box.setIcon(QMessageBox.Icon.Question)
            box.setWindowTitle("Confirm move to backup")
            box.setText(f"Move run {run_number} to backup?")
            box.setInformativeText(
                "The following files will be moved to:\n  "
                + backup_dir + "\n\n  "
                + details
                + "\n\nThe run will also be removed from the current view.")
        else:
            box.setIcon(QMessageBox.Icon.Warning)
            box.setWindowTitle("Confirm run deletion")
            box.setText(f"Are you sure you want to delete run {run_number}?")
            box.setInformativeText(
                "The following files will be permanently removed from disk:\n  "
                + details
                + "\n\nThis cannot be undone.")
        box.setStandardButtons(
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.Cancel)
        box.setDefaultButton(QMessageBox.StandardButton.Cancel)
        yes_btn = box.button(QMessageBox.StandardButton.Yes)
        if yes_btn is not None:
            yes_btn.setText("Move" if backup else "Delete")
        if box.exec() != QMessageBox.StandardButton.Yes:
            return

        if backup:
            try:
                os.makedirs(backup_dir, exist_ok=True)
            except OSError as e:
                QMessageBox.critical(self, "Backup failed",
                                     f"Could not create backup directory:\n"
                                     f"{backup_dir}\n\n{e}")
                return

        errors: List[str] = []
        for p in paths:
            if not os.path.exists(p):
                continue
            dest = os.path.join(backup_dir, os.path.basename(p))
            try:
                if backup:
                    shutil.move(p, dest)
                else:
                    os.remove(p)
            except (OSError, shutil.Error) as e:
                errors.append(f"{p} → {dest}: {e}" if backup else f"{p}: {e}")
        if errors:
            if backup:
                title, what = "Backup partially failed", "moved"
            else:
                title, what = "Delete partially failed", "deleted"
            QMessageBox.warning(self, title,
                                f"Some files could not be {what}:\n\n"
                                + "\n".join(errors))

        self._drop_run_from_state(run_idx, run_number)

    def _on_cycle_palette(self):
        if self._view_mode != 2:  # drift mode uses a fixed palette
            self._map.cycle_palette()

    def _apply_map_mode(self):
        """Map palette and legend of the current view mode."""
        self._map.set_palette_override(
            DRIFT_PALETTE if self._view_mode == 2 else None)
        self._map.set_legend_mode(_VIEW_MODE_LEGEND[self._view_mode])

    def _on_view_mode_changed(self, index: int):
        prev_mode = self._view_mode
        self._view_mode = index
        uses_threshold = (index in (2, 3))
        for w in (self._thresh_lbl, self._thresh_g_input, self._thresh_g_pct,
                  self._thresh_w_lbl, self._thresh_w_input, self._thresh_w_pct):
            w.setVisible(uses_threshold)
        self._run_combo.setEnabled(index != 3)
        self._apply_map_mode()
        if index == 3:
            self._update_summary_views()
            self._hist_widget.clear()
        else:
            self._update_geo_view()
            if not (prev_mode in (0, 1) and index in (0, 1)):
                self._update_irregular_table()
            self._load_lms_hist(self._selected_module)

    def _on_drift_threshold_changed(self):
        for inp, attr, default in (
            (self._thresh_g_input, "_thresh_g", 0.10),
            (self._thresh_w_input, "_thresh_w", 0.05),
        ):
            inp.blockSignals(True)
            try:
                val = float(inp.text())
                if val <= 0:
                    raise ValueError
                setattr(self, attr, val / 100.0)
            except ValueError:
                inp.setText(f"{getattr(self, attr) * 100:.1f}")
            finally:
                inp.blockSignals(False)
        if self._view_mode == 3:
            self._update_summary_views()
        else:
            self._update_geo_view()
            self._update_irregular_table()

    def _on_module_hovered(self, name: str):
        mod = self._mod_by_name.get(name)
        if not mod:
            return
        disp = lms_display_name(name)
        active = self._active_runs
        if not active:
            self._info.setText(f"{disp} ({mod.mod_type})")
            return
        rd = self._runs[self._current_run_idx]
        mrec = rd.modules.get(name)
        ref = LMS_DISPLAY[self._current_ref_idx]
        if mrec:
            gain = mrec.gain_factors[self._current_ref_idx]
            base = (f"{disp} ({mod.mod_type})  gain[{ref}]: {gain:.5f}"
                    f"  lms_peak: {mrec.lms_peak:.2f}"
                    f"  lms_sigma: {mrec.lms_sigma:.2f}")
            if self._view_mode == 1:
                st = self._active_gain_stats().get(name)
                if st and st[2] > 1:
                    mean, std, _n = st
                    dev = (gain - mean) / std if std > 0 else 0.0
                    base += f"  dev: {dev:+.2f}σ"
            self._info.setText(base)
        else:
            self._info.setText(f"{disp} ({mod.mod_type})  no data")

    def _on_module_clicked(self, name: str):
        new_sel = name if name else None
        self._selected_module = new_sel
        if new_sel and new_sel not in LMS_NAMES:
            self._selected_hycal_module = new_sel
        # Clicking one of the three LMS cells doubles as a reference-PMT
        # picker: feed it through the Ref dropdown so the rest of the UI
        # (ref plots, geo view, irregular table) updates exactly as it
        # would from a manual dropdown change.  setCurrentIndex is a no-op
        # when the index is already current.
        if new_sel in LMS_NAMES:
            new_ref_idx = LMS_NAMES.index(new_sel)
            if new_ref_idx != self._current_ref_idx:
                self._ref_combo.setCurrentIndex(new_ref_idx)
        self._map.set_selected(self._selected_module)
        if self._selected_module:
            self._irregular_table.select_module(self._selected_module)
        else:
            self._irregular_table._table.clearSelection()
        self._update_line_charts()
        self._load_lms_hist(self._selected_module)

    def _load_lms_hist(self, module_name: Optional[str],
                       run_number: Optional[int] = None):
        # In summary mode a run_number must be supplied explicitly (from chart click).
        if self._view_mode == 3 and run_number is None:
            self._hist_widget.clear()
            return
        if not module_name or not self._current_folder or not self._runs:
            self._hist_widget.clear()
            return
        if not _UPROOT_OK:
            self._hist_widget.clear()
            return
        if run_number is None:
            run_number = self._runs[self._current_run_idx].run_number
        key = (run_number, module_name)
        if key in self._hist_cache:
            self._hist_cache.move_to_end(key)
            values, edges, gauss, ovl_values, ovl_edges, ovl_gauss = self._hist_cache[key]
        else:
            root_path = os.path.join(self._current_folder,
                                     f"prad_{run_number:06d}_LMS_fitted.root")
            if not os.path.exists(root_path):
                self._hist_widget.clear()
                return

            def _read_hist(rf, hkey):
                """Return (values, edges, gauss) for hkey or (None, None, None)."""
                if hkey not in rf:
                    return None, None, None
                h = rf[hkey]
                vals = h.values()
                edg  = h.axis().edges()
                g    = None
                try:
                    fns = h.member("fFunctions")
                    for fn in fns:
                        if fn.classname == "TF1":
                            params = fn.member("fFormula").member("fClingParameters")
                            if len(params) >= 3:
                                g = (float(params[0]), float(params[1]),
                                     float(params[2]),
                                     float(fn.member("fXmin")),
                                     float(fn.member("fXmax")))
                            break
                except Exception:
                    pass
                return vals, edg, g

            try:
                with uproot.open(root_path) as rf:
                    values, edges, gauss = _read_hist(rf, f"{module_name}_LMS")
                    if values is None:
                        self._hist_widget.clear()
                        return
                    # For the three reference PMTs, also load the alpha-peak
                    # histogram + fit so they overlay on the same plot.
                    if module_name in LMS_NAMES:
                        ovl_values, ovl_edges, ovl_gauss = _read_hist(
                            rf, f"{module_name}_Alpha")
                    else:
                        ovl_values, ovl_edges, ovl_gauss = None, None, None
            except Exception:
                self._hist_widget.clear()
                return
            self._hist_cache[key] = (values, edges, gauss,
                                     ovl_values, ovl_edges, ovl_gauss)
            if len(self._hist_cache) > _HIST_CACHE_MAX:
                self._hist_cache.popitem(last=False)
        self._hist_widget.set_histogram(values, edges,
                                        f"{lms_display_name(module_name)}  run {run_number}",
                                        gauss=gauss,
                                        overlay_values=ovl_values,
                                        overlay_edges=ovl_edges,
                                        overlay_gauss=ovl_gauss)

    # ---- update views ----

    def _recompute_pairwise_diffs(self):
        """Pre-compute all consecutive-run relative gain changes across ALL runs.
        Stored as flat list of (name, mod_type, rel, pair_idx, curr_run,
        prev_run) where pair_idx is the index of the current run in
        self._runs. Called on load, and by _update_summary_views after a ref
        index change, so threshold/range updates just filter this list."""
        ref_idx = self._current_ref_idx
        runs = self._runs
        mod_by_name = self._mod_by_name
        diffs = []
        for i in range(1, len(runs)):
            rd_curr = runs[i]
            rd_prev = runs[i - 1]
            prev_gains = {mname: mrec.gain_factors[ref_idx]
                          for mname, mrec in rd_prev.modules.items()}
            curr_run = rd_curr.run_number
            prev_run = rd_prev.run_number
            for mname, mrec in rd_curr.modules.items():
                g_prev = prev_gains.get(mname)
                if g_prev is None:
                    continue
                rel = _sym_rel_change(mrec.gain_factors[ref_idx], g_prev)
                mod = mod_by_name.get(mname)
                diffs.append((mname, mod.mod_type if mod else "?", rel, i,
                               curr_run, prev_run))
        self._pairwise_diffs = diffs
        self._pairwise_ref_idx = ref_idx

    def _active_gain_stats(self) -> Dict[str, Tuple[float, float, int]]:
        """gain_stats() of the active runs against the current ref."""
        key = (self._start_run_idx, self._end_run_idx, self._current_ref_idx)
        if self._gain_stats_key != key:
            self._gain_stats_cache = gain_stats(self._active_runs,
                                                self._current_ref_idx)
            self._gain_stats_key = key
        return self._gain_stats_cache

    def _prev_active_run(self, active: List[RunData]) -> Optional[RunData]:
        """The active run just before the current one; None if the current
        run is the first active run or outside the active range."""
        curr = self._runs[self._current_run_idx].run_number
        pos = next((i for i, r in enumerate(active)
                    if r.run_number == curr), None)
        return active[pos - 1] if pos else None

    def _update_all_views(self):
        if self._view_mode == 3:
            self._update_summary_views()
        else:
            self._update_geo_view()
            self._update_irregular_table()
        self._update_line_charts()

    def _update_summary_views(self):
        """Filter pre-computed pairwise diffs by active range + threshold."""
        if not self._runs:
            return
        if self._pairwise_ref_idx != self._current_ref_idx:
            self._recompute_pairwise_diffs()
        start_idx = self._start_run_idx
        end_idx = self._end_run_idx
        thresh_g = self._thresh_g
        thresh_w = self._thresh_w
        stats: Dict[str, List] = {}
        for name, mod_type, rel, pair_idx, curr_run, prev_run in self._pairwise_diffs:
            if pair_idx <= start_idx or pair_idx > end_idx:
                continue
            if rel <= _drift_threshold(name, thresh_g, thresh_w):
                continue
            s = stats.get(name)
            if s is None:
                stats[name] = [1, rel, curr_run, prev_run, mod_type]
            else:
                s[0] += 1
                if rel > s[1]:
                    s[1] = rel
                    s[2] = curr_run
                    s[3] = prev_run

        entries = [SummaryEntry(name=n, mod_type=s[4], drift_count=s[0],
                                max_rel_change=s[1], max_run=s[2], max_prev_run=s[3])
                   for n, s in stats.items()]
        entries.sort(key=lambda e: (
        0 if e.name.startswith("W") else 1,
        0 if e.name.startswith("W") else int(math.isinf(e.max_rel_change)),
        -e.drift_count,
        0 if math.isinf(e.max_rel_change) else -e.max_rel_change,
    ))

        values = {e.name: float(e.drift_count) for e in entries}
        if self._auto_range:
            vmax = max(values.values()) if values else 1.0
            self._range_min.setText("0.0000")
            self._range_max.setText(f"{vmax:.4f}")
            self._map.set_gain_data(values, 0.0, vmax)
        else:
            self._map.set_gain_data(values, self._manual_vmin, self._manual_vmax)
        self._irregular_table.set_summary_data(entries)

    def _update_geo_view(self):
        active = self._active_runs
        if not active:
            return
        rd = self._runs[self._current_run_idx]
        ref_idx = self._current_ref_idx

        if self._view_mode == 1:
            # Deviation (σ): signed (gain − mean) / std across active runs
            stats = self._active_gain_stats()
            values: Dict[str, float] = {}
            for mname, mrec in rd.modules.items():
                mean, std, _n = stats.get(mname, (0.0, 0.0, 0))
                if std > 0:
                    values[mname] = (mrec.gain_factors[ref_idx] - mean) / std
                else:
                    values[mname] = 0.0

            if self._auto_range:
                if values:
                    abs_max = max(abs(v) for v in values.values())
                    abs_max = abs_max if abs_max > 0 else 1.0
                    vmin, vmax = -abs_max, abs_max
                else:
                    vmin, vmax = -3.0, 3.0
                self._range_min.setText(f"{vmin:.4f}")
                self._range_max.setText(f"{vmax:.4f}")
            else:
                vmin = self._manual_vmin
                vmax = self._manual_vmax
        elif self._view_mode == 2:
            # Run-to-Run Drift: (gain_current - gain_prev) / gain_prev
            rd_prev = self._prev_active_run(active)
            values: Dict[str, float] = {}
            if rd_prev is not None:
                for mname, mrec in rd.modules.items():
                    prev_mrec = rd_prev.modules.get(mname)
                    if prev_mrec and prev_mrec.gain_factors[ref_idx] != 0:
                        g_prev = prev_mrec.gain_factors[ref_idx]
                        values[mname] = (mrec.gain_factors[ref_idx] - g_prev) / g_prev

            threshold = max(self._thresh_g, self._thresh_w)
            if self._auto_range:
                vmin, vmax = -threshold, threshold
                self._range_min.setText(f"{vmin:.4f}")
                self._range_max.setText(f"{vmax:.4f}")
            else:
                vmin = self._manual_vmin
                vmax = self._manual_vmax

        else:
            # Gain Factor
            values: Dict[str, float] = {}
            for mname, mrec in rd.modules.items():
                values[mname] = mrec.gain_factors[ref_idx]

            if self._auto_range:
                if values:
                    sorted_v = sorted(v for v in values.values() if math.isfinite(v))
                    n = len(sorted_v)
                    lo_idx = max(0, int(n * 0.02))
                    hi_idx = min(n - 1, int(n * 0.98))
                    vmin = sorted_v[lo_idx] if sorted_v else 0.9
                    vmax = sorted_v[hi_idx] if sorted_v else 1.1
                    if vmin == vmax:
                        vmin -= 0.01
                        vmax += 0.01
                else:
                    vmin, vmax = 0.9, 1.1
                self._range_min.setText(f"{vmin:.4f}")
                self._range_max.setText(f"{vmax:.4f}")
            else:
                vmin = self._manual_vmin
                vmax = self._manual_vmax

        self._map.set_gain_data(values, vmin, vmax)

    def _update_line_charts(self):
        active = self._active_runs
        if not active:
            return

        first_run = active[0].run_number
        ref_idx = self._current_ref_idx
        ref_name = LMS_NAMES[ref_idx]
        curr_run_number = self._runs[self._current_run_idx].run_number

        mname = self._selected_hycal_module
        has_module = bool(mname)

        for i, (kind, _label) in enumerate(CHART_PLOTS):
            indices, actual_runs, vals, errs, title = self._series_for_kind(
                kind, ref_idx, ref_name, mname, has_module, active, first_run)
            chart = self._charts[i]
            chart.set_data(indices, vals, errs, title, actual_runs)
            if vals:
                chart.set_y_range(*_chart_y_range(vals, errs))
            chart.set_current_run(curr_run_number)

    def _series_for_kind(
        self, kind: str, ref_idx: int, ref_name: str,
        mname: Optional[str], has_module: bool,
        active: List[RunData], first_run: int,
    ) -> Tuple[List[int], List[int], List[float], List[float], str]:
        """Build (indices, actual_runs, vals, errs, title) for one of the
        five fixed plots."""
        ref_disp = lms_display_name(ref_name)
        if kind in ("mod_lms", "mod_gain") and not has_module:
            what = ("Module LMS" if kind == "mod_lms"
                    else f"Module Gain[{ref_disp}]")
            return [], [], [], [], (
                f"{what} — select a HyCal module  (run1={first_run})")
        mname_disp = lms_display_name(mname) if mname else mname

        def lms_fit(rec):
            if rec is None or rec.lms_peak == 0:
                return None
            return rec.lms_peak, rec.lms_sigma

        def alpha_fit(rec):
            if rec is None or rec.alpha_peak == 0:
                return None
            return rec.alpha_peak, rec.alpha_sigma

        def lms_alpha_ratio(rec):
            if rec is None or rec.alpha_peak == 0 or rec.lms_peak == 0:
                return None
            ratio = rec.lms_peak / rec.alpha_peak
            rel_lms = rec.lms_sigma / rec.lms_peak
            rel_alpha = rec.alpha_sigma / rec.alpha_peak
            return ratio, ratio * math.sqrt(rel_lms ** 2 + rel_alpha ** 2)

        def gain(mrec):
            return None if mrec is None else (mrec.gain_factors[ref_idx], None)

        ref_recs = [rd.lms.get(ref_name) for rd in active]
        mod_recs = [rd.modules.get(mname) for rd in active]
        # kind -> (per-run records, point, title)
        series = {
            "ratio":     (ref_recs, lms_alpha_ratio, f"{ref_disp}  LMS/Alpha"),
            "ref_lms":   (ref_recs, lms_fit,   f"{ref_disp}  LMS peak±σ"),
            "ref_alpha": (ref_recs, alpha_fit, f"{ref_disp}  Alpha peak±σ"),
            "mod_lms":   (mod_recs, lms_fit,   f"{mname_disp}  LMS peak±σ"),
            "mod_gain":  (mod_recs, gain,      f"{mname_disp}  gain[{ref_disp}]"),
        }
        records, point, title = series[kind]
        return (*_collect_series(active, records, point),
                f"{title}  (run1={first_run})")

    def _update_irregular_table(self):
        active = self._active_runs
        if not active:
            return
        if self._view_mode == 2:
            rd_prev = self._prev_active_run(active)
            if rd_prev is not None:
                entries = compute_drift_entries(
                    self._runs[self._current_run_idx], rd_prev,
                    self._current_ref_idx, self._mod_by_name,
                    self._thresh_g, self._thresh_w)
            else:
                entries = []
            self._irregular_table.set_drift_data(entries)
        else:
            entries = compute_irregular_entries(
                active, self._current_ref_idx, self._mod_by_name,
                stats=self._active_gain_stats())
            self._irregular_table.set_data(entries)


# ---- Entry point ----

def main():
    import argparse
    parser = argparse.ArgumentParser()
    parser.add_argument("-dir", dest="folder", default=None,
                        help="Data folder to load on startup")
    parser.add_argument("--theme", choices=available_themes(), default="dark",
                        help="Colour theme (default: dark)")
    args, qt_args = parser.parse_known_args()

    set_theme(args.theme)

    app = QApplication([sys.argv[0]] + qt_args)
    win = GainMonitorWindow()
    if args.folder:
        win._process_btn.setEnabled(False)
        win._load_folder(args.folder)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
