#!/usr/bin/env python3
"""
HyCal Map Builder (PyQt6)
=========================
Simple HyCal geometry viewer that colour-maps user data loaded from
JSON or plain-text files.

Data formats
------------
* JSON  : {"<module_name>": {"<field>": <value>, ...}, ...}

          Values may also be a list of history entries; the last entry
          of each list is used (so gain_equalization_results.json-style
          per-module history files work directly). Nested dicts are
          flattened with dot notation, e.g. fit.slope / edge.percentage.
          Non-numeric fields (strings like timestamps) are ignored.

* Text  : whitespace / comma / tab delimited rows

            <module_name> <val1> <val2> ...

          Lines starting with '#' are ignored. If the first non-comment
          row has a non-numeric second column it is treated as a header
          naming the columns; otherwise columns get default names
          (col1, col2, ...).

Usage
-----
    python scripts/hycal_map_builder.py                        # empty map
    python scripts/hycal_map_builder.py mydata.json            # auto-load
    python scripts/hycal_map_builder.py mydata.txt --field rms
"""

from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QLineEdit, QComboBox, QSlider, QSizePolicy,
    QFileDialog, QMessageBox, QToolTip,
)
from PyQt6.QtCore import Qt, QRectF, QPointF, QSize, pyqtSignal
from PyQt6.QtGui import (
    QPainter, QColor, QPen, QBrush, QFont, QLinearGradient, QPalette,
)


# ===========================================================================
#  Paths
# ===========================================================================

SCRIPT_DIR = Path(__file__).resolve().parent

def _find_modules_json() -> Path:
    candidates = [
        SCRIPT_DIR / ".." / "database" / "hycal_modules.json",
        Path.cwd() / "database" / "hycal_modules.json",
        Path.cwd() / "hycal_modules.json",
    ]
    for c in candidates:
        if c.is_file():
            return c.resolve()
    return (SCRIPT_DIR / ".." / "database" / "hycal_modules.json").resolve()

MODULES_JSON = _find_modules_json()


# ===========================================================================
#  Module database
# ===========================================================================

class Module:
    __slots__ = ("name", "mod_type", "x", "y", "sx", "sy")
    def __init__(self, name, mod_type, x, y, sx, sy):
        self.name = name
        self.mod_type = mod_type
        self.x = x
        self.y = y
        self.sx = sx
        self.sy = sy


def load_modules(path: Path) -> List[Module]:
    with open(path) as f:
        data = json.load(f)
    return [Module(e["n"], e["t"], e["x"], e["y"], e["sx"], e["sy"])
            for e in data]


# ===========================================================================
#  Data loading
# ===========================================================================

def load_data_file(path: Path) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    """
    Returns (data, fields) where
      data[field][module_name] -> float
      fields is the ordered list of field names
    """
    text = path.read_text()
    try:
        obj = json.loads(text)
        if isinstance(obj, dict):
            return _data_from_json_dict(obj)
    except json.JSONDecodeError:
        pass
    return _data_from_text(text)


def _flatten_dict(d: Dict, prefix: str = "") -> Dict[str, object]:
    out: Dict[str, object] = {}
    for k, v in d.items():
        key = f"{prefix}{k}"
        if isinstance(v, dict):
            out.update(_flatten_dict(v, key + "."))
        else:
            out[key] = v
    return out


def _data_from_json_dict(obj: Dict) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    # Normalise each module's value to a flat dict of scalars.
    #   {name: {key: val, ...}}              -> as-is
    #   {name: [{entry}, {entry}, ...]}      -> last entry (history)
    # Nested dicts are flattened with dot-joined keys.
    per_module: Dict[str, Dict[str, object]] = {}
    for name, entry in obj.items():
        if isinstance(entry, list):
            if not entry:
                continue
            entry = entry[-1]
        if not isinstance(entry, dict):
            continue
        flat = _flatten_dict(entry)
        if flat:
            per_module[str(name)] = flat

    fields: List[str] = []
    seen = set()
    for flat in per_module.values():
        for k in flat:
            if k not in seen:
                seen.add(k)
                fields.append(k)

    data: Dict[str, Dict[str, float]] = {f: {} for f in fields}
    for name, flat in per_module.items():
        for k, v in flat.items():
            try:
                data[k][name] = float(v)
            except (TypeError, ValueError):
                pass

    # Drop fields that ended up with no numeric values (e.g. timestamps).
    data = {k: v for k, v in data.items() if v}
    fields = [f for f in fields if f in data]
    return data, fields


_SPLIT = re.compile(r"[,\s\t]+")

def _is_number(s: str) -> bool:
    try:
        float(s)
        return True
    except ValueError:
        return False

def _data_from_text(text: str) -> Tuple[Dict[str, Dict[str, float]], List[str]]:
    rows: List[List[str]] = []
    for line in text.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        parts = [p for p in _SPLIT.split(line) if p]
        if parts:
            rows.append(parts)
    if not rows:
        return {}, []

    first = rows[0]
    ncols = max(len(r) for r in rows) - 1
    if ncols <= 0:
        return {}, []

    header = not all(_is_number(c) for c in first[1:])
    if header:
        fields = [f or f"col{i+1}" for i, f in enumerate(first[1:])]
        data_rows = rows[1:]
    else:
        fields = [f"col{i+1}" for i in range(ncols)]
        data_rows = rows

    # pad/trim field list to ncols
    while len(fields) < ncols:
        fields.append(f"col{len(fields)+1}")
    fields = fields[:ncols]

    data: Dict[str, Dict[str, float]] = {f: {} for f in fields}
    for row in data_rows:
        if len(row) < 2:
            continue
        name = row[0]
        for i, field in enumerate(fields):
            idx = i + 1
            if idx >= len(row):
                break
            try:
                data[field][name] = float(row[idx])
            except ValueError:
                pass
    return data, fields


# ===========================================================================
#  Palettes
# ===========================================================================

PALETTES = {
    "viridis": [
        (0.00, ( 68,   1,  84)), (0.25, ( 59,  82, 139)),
        (0.50, ( 33, 145, 140)), (0.75, ( 94, 201,  98)),
        (1.00, (253, 231,  37)),
    ],
    "inferno": [
        (0.00, (  0,   0,   4)), (0.25, (120,  28, 109)),
        (0.50, (229,  89,  52)), (0.75, (253, 198,  39)),
        (1.00, (252, 255, 164)),
    ],
    "coolwarm": [
        (0.00, ( 59,  76, 192)), (0.25, (141, 176, 254)),
        (0.50, (221, 221, 221)), (0.75, (245, 148, 114)),
        (1.00, (180,   4,  38)),
    ],
    "hot": [
        (0.00, ( 11,   0,   0)), (0.33, (230,   0,   0)),
        (0.66, (255, 210,   0)), (1.00, (255, 255, 255)),
    ],
    "rainbow": [
        (0.00, ( 30,  58,  95)), (0.25, ( 59, 130, 246)),
        (0.50, ( 45, 212, 160)), (0.75, (234, 179,   8)),
        (1.00, (245, 101, 101)),
    ],
    "greyscale": [
        (0.00, ( 20,  20,  20)), (1.00, (240, 240, 240)),
    ],
}
PALETTE_NAMES = list(PALETTES.keys())


def _lerp(a: int, b: int, t: float) -> int:
    return int(a + (b - a) * t)

def _cmap_qcolor(t: float, stops) -> QColor:
    t = max(0.0, min(1.0, t))
    for i in range(len(stops) - 1):
        t0, c0 = stops[i]
        t1, c1 = stops[i + 1]
        if t <= t1:
            s = (t - t0) / (t1 - t0) if t1 > t0 else 0.0
            return QColor(_lerp(c0[0], c1[0], s),
                          _lerp(c0[1], c1[1], s),
                          _lerp(c0[2], c1[2], s))
    _, c = stops[-1]
    return QColor(*c)


# ===========================================================================
#  HyCal map widget  (zoom / pan / hover / colourbar-click to cycle palette)
# ===========================================================================

class HyCalMapWidget(QWidget):
    moduleHovered = pyqtSignal(str)
    paletteClicked = pyqtSignal()

    _SHRINK = 0.92
    _CLICK_THRESHOLD = 4

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumSize(500, 500)

        self._modules: List[Module] = []
        self._pbglass_names: set = set()
        self._pbglass_alpha: float = 1.0
        self._values: Dict[str, float] = {}
        self._vmin = 0.0
        self._vmax = 1.0
        self._log_scale = False
        self._palette_idx = 0
        self._field_label = ""
        self._hovered: Optional[str] = None
        self._rects: Dict[str, QRectF] = {}
        self._rect_names_rev: List[str] = []
        self._geo_bounds = (0.0, 1.0, 0.0, 1.0)
        self._cb_rect: Optional[QRectF] = None
        self._layout_dirty = True

        # zoom / pan
        self._zoom = 1.0
        self._pan_x = 0.0
        self._pan_y = 0.0
        self._drag_last: Optional[QPointF] = None
        self._drag_origin: Optional[QPointF] = None
        self._dragging = False

        # overlay reset button
        self._reset_btn = QPushButton("Reset", self)
        self._reset_btn.setFixedSize(52, 24)
        _rf = QFont("Consolas", 9); _rf.setBold(True)
        self._reset_btn.setFont(_rf)
        self._reset_btn.setStyleSheet(
            "QPushButton{background:rgba(22,27,34,200);color:#8b949e;"
            "border:1px solid #30363d;border-radius:3px;}"
            "QPushButton:hover{background:rgba(33,38,45,220);color:#c9d1d9;}")
        self._reset_btn.clicked.connect(self.reset_view)

    # -- public API --

    def set_modules(self, modules: List[Module]):
        self._modules = [m for m in modules if m.mod_type != "LMS"]
        self._pbglass_names = {m.name for m in self._modules if m.mod_type == "PbGlass"}
        if self._modules:
            self._geo_bounds = (
                min(m.x - m.sx / 2 for m in self._modules),
                max(m.x + m.sx / 2 for m in self._modules),
                min(m.y - m.sy / 2 for m in self._modules),
                max(m.y + m.sy / 2 for m in self._modules),
            )
        self._layout_dirty = True
        self.update()

    def set_values(self, values: Dict[str, float], label: str = ""):
        self._values = values
        self._field_label = label
        self.update()

    def set_range(self, vmin: float, vmax: float):
        self._vmin = vmin
        self._vmax = vmax
        self.update()

    def set_log_scale(self, on: bool):
        self._log_scale = on
        self.update()

    def set_palette(self, idx: int):
        self._palette_idx = idx % len(PALETTES)
        self.update()

    def palette_idx(self) -> int:
        return self._palette_idx

    def set_pbglass_alpha(self, a: float):
        self._pbglass_alpha = max(0.0, min(1.0, a))
        self.update()

    def reset_view(self):
        self._zoom = 1.0
        self._pan_x = 0.0
        self._pan_y = 0.0
        self._layout_dirty = True
        self.update()

    # -- layout --

    def _recompute_layout(self):
        self._rects.clear()
        if not self._modules:
            self._rect_names_rev = []
            return
        w, h = self.width(), self.height()
        margin, top, bot = 12, 8, 50
        pw, ph = w - 2 * margin, h - top - bot
        x0, x1, y0, y1 = self._geo_bounds
        base_scale = min(pw / (x1 - x0), ph / (y1 - y0))
        sc = base_scale * self._zoom
        dw, dh = (x1 - x0) * sc, (y1 - y0) * sc
        ox = margin + (pw - dw) / 2 + self._pan_x
        oy = top + (ph - dh) / 2 + self._pan_y
        shrink = self._SHRINK
        for m in self._modules:
            mw, mh = m.sx * sc * shrink, m.sy * sc * shrink
            cx = ox + (m.x - x0) * sc
            cy = oy + (y1 - m.y) * sc
            self._rects[m.name] = QRectF(cx - mw / 2, cy - mh / 2, mw, mh)
        self._rect_names_rev = list(self._rects)[::-1]
        self._layout_dirty = False

    def resizeEvent(self, event):
        self._layout_dirty = True
        self._reset_btn.move(self.width() - self._reset_btn.width() - 6, 6)
        super().resizeEvent(event)

    # -- painting --

    def paintEvent(self, event):
        if self._layout_dirty:
            self._recompute_layout()
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor("#0a0e14"))

        if not self._rects:
            p.setPen(QColor("#555"))
            p.setFont(QFont("Consolas", 12))
            p.drawText(QRectF(0, 0, w, h),
                       Qt.AlignmentFlag.AlignCenter,
                       "No modules loaded")
            p.end()
            return

        stops = list(PALETTES.values())[self._palette_idx]
        vmin, vmax = self._vmin, self._vmax
        no_data = QColor("#1a1a2e")

        if self._log_scale:
            log_lo = math.log10(max(vmin, 1e-9))
            log_hi = math.log10(max(vmax, vmin * 10, 1e-8))

        glass_alpha = self._pbglass_alpha
        frame_col_base = QColor(160, 165, 175)

        for name, rect in self._rects.items():
            is_glass = name in self._pbglass_names
            a = glass_alpha if is_glass else 1.0

            v = self._values.get(name)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                fill = QColor(no_data)
            else:
                if self._log_scale:
                    lv = math.log10(max(v, 1e-9))
                    t = (lv - log_lo) / (log_hi - log_lo) if log_hi > log_lo else 0.5
                else:
                    t = ((v - vmin) / (vmax - vmin)) if vmax > vmin else 0.5
                fill = _cmap_qcolor(t, stops)

            if a < 1.0:
                fill = QColor(fill)
                fill.setAlphaF(a)
            if a > 0.0:
                p.fillRect(rect, fill)

            if is_glass and a > 0.0:
                frame = QColor(frame_col_base)
                frame.setAlphaF(a * 0.8)
                p.setPen(QPen(frame, 1.0))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRect(rect)

        # hover highlight
        if self._hovered and self._hovered in self._rects:
            p.setPen(QPen(QColor("#58a6ff"), 2.0))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._hovered])

        # colour bar
        cb_w = min(400, w - 80)
        cb_h = 14
        cb_x = (w - cb_w) / 2
        cb_y = h - 40
        self._cb_rect = QRectF(cb_x, cb_y, cb_w, cb_h)
        grad = QLinearGradient(cb_x, 0, cb_x + cb_w, 0)
        for t, (r, g, b) in stops:
            grad.setColorAt(t, QColor(r, g, b))
        p.fillRect(self._cb_rect, QBrush(grad))
        p.setPen(QPen(QColor("#58a6ff"), 1.0))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(self._cb_rect)

        # labels
        p.setPen(QColor("#8b949e"))
        p.setFont(QFont("Consolas", 9))
        p.drawText(QRectF(cb_x, cb_y + cb_h + 2, 120, 14),
                   Qt.AlignmentFlag.AlignLeft, _fmt(vmin))
        p.drawText(QRectF(cb_x + cb_w - 120, cb_y + cb_h + 2, 120, 14),
                   Qt.AlignmentFlag.AlignRight, _fmt(vmax))
        mid = PALETTE_NAMES[self._palette_idx]
        if self._log_scale:
            mid += "  [log]"
        if self._field_label:
            mid = self._field_label + "  \u2014  " + mid
        p.drawText(QRectF(cb_x, cb_y + cb_h + 2, cb_w, 14),
                   Qt.AlignmentFlag.AlignCenter, mid)
        p.end()

    # -- mouse --

    def mousePressEvent(self, e):
        if e.button() == Qt.MouseButton.MiddleButton:
            self.reset_view()
            return
        if e.button() in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            self._drag_last = e.position()
            self._drag_origin = e.position()
            self._dragging = False

    def mouseReleaseEvent(self, e):
        if e.button() in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            if self._dragging:
                self.setCursor(Qt.CursorShape.ArrowCursor)
            elif e.button() == Qt.MouseButton.LeftButton:
                if self._cb_rect and self._cb_rect.contains(e.position()):
                    self.paletteClicked.emit()
            self._drag_last = None
            self._drag_origin = None
            self._dragging = False

    def mouseMoveEvent(self, e):
        if self._drag_last is not None:
            pos = e.position()
            if not self._dragging:
                dx = pos.x() - self._drag_origin.x()
                dy = pos.y() - self._drag_origin.y()
                if dx * dx + dy * dy > self._CLICK_THRESHOLD ** 2:
                    self._dragging = True
                    self.setCursor(Qt.CursorShape.ClosedHandCursor)
            if self._dragging:
                self._pan_x += pos.x() - self._drag_last.x()
                self._pan_y += pos.y() - self._drag_last.y()
                self._drag_last = pos
                self._layout_dirty = True
                self.update()
            return

        pos = e.position()
        if self._cb_rect and self._cb_rect.contains(pos):
            self.setCursor(Qt.CursorShape.PointingHandCursor)
        else:
            self.setCursor(Qt.CursorShape.ArrowCursor)

        found = None
        for name in self._rect_names_rev:
            if self._rects[name].contains(pos):
                found = name
                break
        if found != self._hovered:
            self._hovered = found
            self.update()
            if found:
                v = self._values.get(found)
                tip = f"{found}: {_fmt(v)}" if v is not None else found
                QToolTip.showText(e.globalPosition().toPoint(), tip, self)
                self.moduleHovered.emit(found)
            else:
                QToolTip.hideText()

    def wheelEvent(self, e):
        factor = 1.15 if e.angleDelta().y() > 0 else 1.0 / 1.15
        new_zoom = max(0.5, min(self._zoom * factor, 20.0))
        if new_zoom == self._zoom:
            return
        pos = e.position()
        ratio = new_zoom / self._zoom
        self._pan_x = pos.x() + (self._pan_x - pos.x()) * ratio
        self._pan_y = pos.y() + (self._pan_y - pos.y()) * ratio
        self._zoom = new_zoom
        self._layout_dirty = True
        self.update()

    def sizeHint(self):
        return QSize(760, 760)


def _fmt(v) -> str:
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "\u2014"
    if v == 0:
        return "0"
    return f"{v:.6g}"


# ===========================================================================
#  Main window
# ===========================================================================

class MapBuilderWindow(QMainWindow):

    def __init__(self, modules: List[Module],
                 data_file: Optional[Path] = None,
                 initial_field: Optional[str] = None):
        super().__init__()
        self._modules = modules
        self._data: Dict[str, Dict[str, float]] = {}
        self._fields: List[str] = []
        self._current_field: Optional[str] = None
        self._data_path: Optional[Path] = None
        self._auto_scale_on = True

        self._build_ui()
        self._map.set_modules(modules)

        if data_file is not None:
            self._load_file(data_file, preferred_field=initial_field)

    # -- ui --

    def _build_ui(self):
        self.setWindowTitle("HyCal Map Builder")
        self.resize(900, 960)
        self._apply_dark_palette()

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # -- top bar: title + file --
        top = QHBoxLayout()
        title = QLabel("HYCAL MAP BUILDER")
        title.setFont(QFont("Monospace", 14, QFont.Weight.Bold))
        title.setStyleSheet("color:#58a6ff;")
        top.addWidget(title)

        top.addStretch()

        self._file_lbl = QLabel("(no file loaded)")
        self._file_lbl.setFont(QFont("Monospace", 10))
        self._file_lbl.setStyleSheet("color:#8b949e;")
        top.addWidget(self._file_lbl)

        top.addWidget(self._make_btn("Open File...", "#58a6ff", self._open_file))
        root.addLayout(top)

        # -- map --
        self._map = HyCalMapWidget()
        self._map.paletteClicked.connect(self._cycle_palette)
        self._map.moduleHovered.connect(self._on_hover)
        root.addWidget(self._map, stretch=1)

        # -- controls: field + palette + auto + range + log --
        ctrl = QHBoxLayout()
        ctrl.addWidget(self._styled_label("Field:"))
        self._field_box = QComboBox()
        self._field_box.setMinimumWidth(160)
        self._field_box.setFont(QFont("Monospace", 11))
        self._field_box.setStyleSheet(
            "QComboBox{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}"
            "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;}")
        self._field_box.currentTextChanged.connect(self._on_field_changed)
        ctrl.addWidget(self._field_box)

        ctrl.addSpacing(12)
        ctrl.addWidget(self._make_btn("Palette \u25B6", "#c9d1d9",
                                      self._cycle_palette))

        ctrl.addSpacing(12)
        ctrl.addWidget(self._styled_label("Range:"))
        self._min_edit = self._styled_edit("0")
        self._max_edit = self._styled_edit("1")
        ctrl.addWidget(self._min_edit)
        ctrl.addWidget(self._styled_label("-"))
        ctrl.addWidget(self._max_edit)
        ctrl.addWidget(self._make_btn("Apply", "#c9d1d9", self._apply_range))

        self._auto_btn = self._make_btn("Auto 2-98%", "#d29922",
                                        self._toggle_auto_range)
        ctrl.addWidget(self._auto_btn)
        self._update_auto_btn()

        self._log_btn = self._make_btn("Log: OFF", "#8b949e", self._toggle_log)
        ctrl.addWidget(self._log_btn)

        ctrl.addSpacing(12)
        ctrl.addWidget(self._styled_label("PbGlass \u03B1:"))
        self._alpha_slider = QSlider(Qt.Orientation.Horizontal)
        self._alpha_slider.setRange(0, 100)
        self._alpha_slider.setValue(100)
        self._alpha_slider.setFixedWidth(120)
        self._alpha_slider.setStyleSheet(
            "QSlider::groove:horizontal{background:#30363d;height:4px;border-radius:2px;}"
            "QSlider::sub-page:horizontal{background:#58a6ff;height:4px;border-radius:2px;}"
            "QSlider::handle:horizontal{background:#58a6ff;width:12px;"
            "margin:-5px 0;border-radius:6px;}"
            "QSlider::handle:horizontal:hover{background:#79b8ff;}")
        self._alpha_slider.valueChanged.connect(self._on_alpha_changed)
        ctrl.addWidget(self._alpha_slider)
        self._alpha_lbl = QLabel("100%")
        self._alpha_lbl.setFont(QFont("Monospace", 10))
        self._alpha_lbl.setStyleSheet("color:#8b949e;")
        self._alpha_lbl.setFixedWidth(40)
        ctrl.addWidget(self._alpha_lbl)

        ctrl.addStretch()
        root.addLayout(ctrl)

        # -- info / stats --
        info_row = QHBoxLayout()
        info_row.setSpacing(6)
        self._info = QLabel("Hover over a module")
        self._info.setFont(QFont("Monospace", 11))
        self._info.setStyleSheet(
            "QLabel{background:#161b22;color:#c9d1d9;padding:4px 8px;"
            "border:1px solid #30363d;border-radius:4px;}")
        self._info.setFixedHeight(28)
        info_row.addWidget(self._info, stretch=1)

        self._stats_lbl = QLabel("")
        self._stats_lbl.setFont(QFont("Monospace", 11))
        self._stats_lbl.setStyleSheet(
            "QLabel{background:#161b22;color:#8b949e;padding:4px 8px;"
            "border:1px solid #30363d;border-radius:4px;}")
        self._stats_lbl.setFixedHeight(28)
        info_row.addWidget(self._stats_lbl)
        root.addLayout(info_row)

    def _make_btn(self, text: str, fg: str, slot) -> QPushButton:
        btn = QPushButton(text)
        f = QFont("Monospace", 11); f.setBold(True)
        btn.setFont(f)
        btn.setStyleSheet(
            f"QPushButton{{background:#21262d;color:{fg};"
            f"border:1px solid #30363d;padding:5px 14px;"
            f"border-radius:4px;}}"
            f"QPushButton:hover{{background:#30363d;}}")
        btn.clicked.connect(slot)
        return btn

    def _styled_label(self, text: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setFont(QFont("Monospace", 11))
        lbl.setStyleSheet("color:#c9d1d9;")
        return lbl

    def _styled_edit(self, text: str) -> QLineEdit:
        e = QLineEdit(text)
        e.setFixedWidth(90)
        e.setFont(QFont("Monospace", 11))
        e.setStyleSheet(
            "QLineEdit{background:#161b22;color:#c9d1d9;"
            "border:1px solid #30363d;border-radius:3px;padding:2px 6px;}")
        e.returnPressed.connect(self._apply_range)
        return e

    def _apply_dark_palette(self):
        pal = self.palette()
        for role, colour in [
            (QPalette.ColorRole.Window,     "#0d1117"),
            (QPalette.ColorRole.WindowText, "#c9d1d9"),
            (QPalette.ColorRole.Base,       "#161b22"),
            (QPalette.ColorRole.Text,       "#c9d1d9"),
            (QPalette.ColorRole.Button,     "#21262d"),
            (QPalette.ColorRole.ButtonText, "#c9d1d9"),
            (QPalette.ColorRole.Highlight,  "#58a6ff"),
        ]:
            pal.setColor(role, QColor(colour))
        self.setPalette(pal)

    # -- actions --

    def _open_file(self):
        start_dir = str(self._data_path.parent) if self._data_path else str(Path.cwd())
        path_str, _ = QFileDialog.getOpenFileName(
            self, "Open data file", start_dir,
            "Data files (*.json *.txt *.dat *.csv *.tsv);;All files (*)")
        if path_str:
            self._load_file(Path(path_str))

    def _load_file(self, path: Path, preferred_field: Optional[str] = None):
        try:
            data, fields = load_data_file(path)
        except Exception as ex:
            QMessageBox.warning(self, "Load failed", f"{path}\n\n{ex}")
            return
        if not fields:
            QMessageBox.warning(self, "No data", f"No usable fields in {path}")
            return

        self._data = data
        self._fields = fields
        self._data_path = path
        self._file_lbl.setText(path.name)
        self._file_lbl.setToolTip(str(path))

        self._field_box.blockSignals(True)
        self._field_box.clear()
        self._field_box.addItems(fields)
        if preferred_field and preferred_field in fields:
            self._field_box.setCurrentText(preferred_field)
        self._field_box.blockSignals(False)
        self._current_field = self._field_box.currentText()
        self._refresh_values()

    def _on_field_changed(self, field: str):
        if not field:
            return
        self._current_field = field
        self._refresh_values()

    def _refresh_values(self):
        if not self._current_field:
            return
        values = self._data.get(self._current_field, {})
        self._map.set_values(values, label=self._current_field)
        if self._auto_scale_on:
            self._do_auto_range()
        else:
            self._map.set_range(self._parse_float(self._min_edit.text(), 0.0),
                                self._parse_float(self._max_edit.text(), 1.0))
        self._update_stats()

    def _update_stats(self):
        if not self._current_field:
            self._stats_lbl.setText("")
            return
        vals = [v for v in self._data.get(self._current_field, {}).values()
                if v is not None and not (isinstance(v, float) and math.isnan(v))]
        if not vals:
            self._stats_lbl.setText("no data")
            return
        arr = np.asarray(vals, dtype=float)
        self._stats_lbl.setText(
            f"N={arr.size}  mean={_fmt(float(arr.mean()))}"
            f"  rms={_fmt(float(arr.std()))}")

    @staticmethod
    def _parse_float(s: str, default: float) -> float:
        try:
            return float(s)
        except ValueError:
            return default

    def _apply_range(self):
        try:
            vmin = float(self._min_edit.text())
            vmax = float(self._max_edit.text())
        except ValueError:
            return
        if vmin >= vmax:
            return
        self._map.set_range(vmin, vmax)
        self._auto_scale_on = False
        self._update_auto_btn()

    def _toggle_auto_range(self):
        self._auto_scale_on = not self._auto_scale_on
        self._update_auto_btn()
        if self._auto_scale_on:
            self._do_auto_range()

    def _do_auto_range(self):
        if not self._current_field:
            return
        vals = list(self._data.get(self._current_field, {}).values())
        arr = np.asarray([v for v in vals
                          if v is not None and not (isinstance(v, float) and math.isnan(v))],
                         dtype=float)
        if arr.size == 0:
            vmin, vmax = 0.0, 1.0
        else:
            vmin = float(np.percentile(arr, 2))
            vmax = float(np.percentile(arr, 98))
            if vmin == vmax:
                pad = abs(vmin) * 0.05 if vmin != 0 else 1.0
                vmin -= pad
                vmax += pad
        self._map.set_range(vmin, vmax)
        self._min_edit.setText(_fmt(vmin))
        self._max_edit.setText(_fmt(vmax))

    def _update_auto_btn(self):
        if self._auto_scale_on:
            self._auto_btn.setStyleSheet(
                "QPushButton{background:#d29922;color:#0d1117;"
                "border:1px solid #d29922;padding:5px 14px;"
                "border-radius:4px;}"
                "QPushButton:hover{background:#e0a82b;}")
        else:
            self._auto_btn.setStyleSheet(
                "QPushButton{background:#21262d;color:#d29922;"
                "border:1px solid #30363d;padding:5px 14px;"
                "border-radius:4px;}"
                "QPushButton:hover{background:#30363d;}")

    def _cycle_palette(self):
        self._map.set_palette(self._map.palette_idx() + 1)

    def _on_alpha_changed(self, v: int):
        self._map.set_pbglass_alpha(v / 100.0)
        self._alpha_lbl.setText(f"{v}%")

    def _toggle_log(self):
        on = not self._map._log_scale
        self._map.set_log_scale(on)
        if on:
            self._log_btn.setText("Log: ON")
            self._log_btn.setStyleSheet(
                self._log_btn.styleSheet().replace("#8b949e", "#58a6ff"))
        else:
            self._log_btn.setText("Log: OFF")
            self._log_btn.setStyleSheet(
                self._log_btn.styleSheet().replace("#58a6ff", "#8b949e"))

    def _on_hover(self, name: str):
        parts = [name]
        for m in self._modules:
            if m.name == name:
                parts.append(f"({m.mod_type})")
                break
        if self._current_field:
            v = self._data.get(self._current_field, {}).get(name)
            if v is not None:
                parts.append(f"{self._current_field} = {_fmt(v)}")
        self._info.setText("    ".join(parts))


# ===========================================================================
#  Main
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(description="HyCal geo-view map builder")
    ap.add_argument("data_file", nargs="?", type=Path,
                    help="JSON or plain-text data file to load at startup")
    ap.add_argument("--field", type=str, default=None,
                    help="Initial field to display (default: first field)")
    ap.add_argument("--modules", type=Path, default=MODULES_JSON,
                    help=f"Path to hycal_modules.json (default: {MODULES_JSON})")
    args = ap.parse_args()

    if not args.modules.is_file():
        print(f"ERROR: hycal modules file not found: {args.modules}",
              file=sys.stderr)
        sys.exit(1)

    modules = load_modules(args.modules)
    print(f"Loaded {len(modules)} HyCal modules from {args.modules}")

    app = QApplication(sys.argv)
    win = MapBuilderWindow(modules,
                           data_file=args.data_file,
                           initial_field=args.field)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
