#!/usr/bin/env python3
"""
FADC Gain Config Editor (PyQt6)
================================
Interactive HyCal geo-view editor for the FADC250 gain config
(``adchycal_gain.cnf``).  Always opens a GUI: load a calibration JSON,
edit per-channel gains by clicking or dragging on the HyCal map, and
save the resulting trigger config.

Workflow inside the GUI
-----------------------
* **Load Calibration…** – read a JSON list of ``{name, factor, ...}``
  entries; gains for matched modules become the displayed colormap.
* **Set / Set All** – left-click / drag paints the value in the line
  edit (or click *Set All* to bulk-apply).  Use ``0`` to mask channels.
* **Safe Cap** – clamp every channel's gain to a configurable
  ``[Min, Max]`` (defaults ``[0.0, 0.15]``); gains stay non-negative.
* **Undo / Reset** – revert the last action / discard all manual edits.
* **Load .cnf… / Save .cnf…** – open / save an ``adchycal_gain.cnf`` file.

Optional CLI shortcuts (still always open the GUI):
    python fadc_gain_config.py
    python fadc_gain_config.py -c database/calibration/adc_to_mev_factors_cosmic.json
    python fadc_gain_config.py --pbwo4-gain 0.15 --pbglass-gain 0.12
    python fadc_gain_config.py -o /path/to/adchycal_gain.cnf
    python fadc_gain_config.py -i existing.cnf
    python fadc_gain_config.py -d /path/to/database --theme light
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Dict, List, Optional, Set, Tuple

from PyQt6.QtCore import Qt, QRectF, pyqtSignal
from PyQt6.QtGui import QColor, QPen, QFont, QDoubleValidator
from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QLineEdit, QTextEdit, QSplitter, QFileDialog,
    QDoubleSpinBox, QGroupBox, QFormLayout, QToolTip, QInputDialog,
    QMessageBox,
)

from daq_common import btn_style, iter_fav3, load_module_info, render_fav3
from prad2_env import find_database_file
from hycal_geoview import (
    HyCalMapWidget as _HyCalMapBase,
    AUX_TYPES, CHANNELS_PER_SLOT, CRATE_NAMES, Module,
    ColorRangeControl,
    THEME, apply_theme_palette, set_theme, available_themes, themed,
    cmap_qcolor,
)


DEFAULT_UNMAPPED_GAIN = 0.0    # nonexistent channel — disable
DEFAULT_LMS_GAIN  = 1.0
DEFAULT_VETO_GAIN = 1.0


# ---- Database auto-discovery ----

def find_database_dir(explicit: Optional[str] = None) -> Path:
    if explicit:
        p = Path(explicit).resolve()
        if not p.is_dir():
            sys.exit(f"error: --database path does not exist: {p}")
        return p

    found = find_database_file("hycal_map.json", use_env=False)
    if found is None:
        sys.exit("error: could not locate database directory "
                 "(looked for hycal_map.json)")
    return found.parent


# ---- Gain source ----

def load_calibration(path: Path) -> Dict[str, float]:
    """Return {module_name: gain_factor} from a calibration JSON file."""
    with open(path) as f:
        data = json.load(f)
    out: Dict[str, float] = {}
    for entry in data:
        name = entry.get("name")
        factor = entry.get("factor")
        if name is None or factor is None:
            continue
        out[name] = float(factor)
    return out


def resolve_gain(name: str,
                 mod_type: Optional[str],
                 cal: Dict[str, float],
                 pbwo4_gain: float,
                 pbglass_gain: float) -> float:
    if name in cal:
        return cal[name]
    if mod_type == "PbWO4":
        return pbwo4_gain
    if mod_type == "PbGlass":
        return pbglass_gain
    if mod_type == "LMS":
        return DEFAULT_LMS_GAIN
    # V1-V4 vetos and anything else
    return DEFAULT_VETO_GAIN


def safe_cap_gains(gains: Dict[str, float],
                   min_cap: float = 0.0,
                   max_cap: float = 0.15) -> Dict[str, float]:
    """Clamp each gain in ``gains`` to ``[max(0, min_cap), max_cap]``.

    The lower bound is floored at 0 — gain factors are physical and must
    never be negative.  Raises ``ValueError`` if ``max_cap`` is below the
    (non-negative) lower bound.
    """
    lo = max(0.0, float(min_cap))
    hi = float(max_cap)
    if hi < lo:
        raise ValueError(f"max_cap ({hi}) must be >= min_cap ({lo})")
    return {k: lo if v < lo else (hi if v > hi else v)
            for k, v in gains.items()}


# ---- Config text rendering / parsing ----

def format_gain(g: float) -> str:
    return f"{g:.6f}"


def render_cnf(daq: List[Tuple[str, int, int, int]],
               gains_by_name: Dict[str, float],
               header_comments: Optional[List[str]] = None) -> str:
    """Build the ``.cnf`` text from a per-module gain dict.

    Channels for which no module is present in the DAQ map at the given
    (crate, slot, channel) get the fallback :data:`DEFAULT_UNMAPPED_GAIN`.
    """
    slots: Dict[Tuple[int, int], Dict[int, Tuple[str, float]]] = {}
    for name, crate, slot, ch in daq:
        gain = gains_by_name.get(name, DEFAULT_UNMAPPED_GAIN)
        slots.setdefault((crate, slot), {})[ch] = (name, gain)

    blocks: Dict[Tuple[int, int], Tuple[str, List[str]]] = {}
    for (crate, slot), ch_map in slots.items():
        gains: List[str] = []
        names: List[str] = []
        for ch in range(CHANNELS_PER_SLOT):
            name, g = ch_map.get(ch, (f"ch{ch}:unmapped",
                                      DEFAULT_UNMAPPED_GAIN))
            gains.append(format_gain(g))
            names.append(name)
        blocks[(crate, slot)] = (f"# slot {slot}: {', '.join(names)}", gains)

    lines: List[str] = ["# adchycal_gain.cnf",
                        "# Generated by fadc_gain_config.py"]
    if header_comments:
        lines.extend(header_comments)
    lines.append("")
    lines += render_fav3(blocks, "FAV3_ALLCH_GAIN")
    return "\n".join(lines)


def parse_cnf_text(text: str,
                   daq: List[Tuple[str, int, int, int]]) -> Dict[str, float]:
    """Parse ``.cnf`` text, return ``{module_name: gain}`` for mapped channels."""
    daq_lookup = {(crate, slot, ch): name for name, crate, slot, ch in daq}
    gains: Dict[str, float] = {}
    for crate, slot, vals in iter_fav3(text, "FAV3_ALLCH_GAIN"):
        for ch, val in enumerate(vals[:CHANNELS_PER_SLOT]):
            name = daq_lookup.get((crate, slot, ch))
            if name:
                try:
                    gains[name] = float(val)
                except ValueError:
                    pass
    return gains


# ---- HyCal geo-view widget ----
# Two interaction modes (selected via the editor's right-panel buttons):
#   * Edit (default): click on a module emits ``moduleEditRequested`` with
#     the current gain.  The editor opens a popup dialog to set a new value.
#   * Set:  drag-paint the value carried by ``_paint_value``.
# At drag end ``paintCommitted`` fires with the batch of
# ``(name, prior_override_or_None)`` tuples so the editor can record the
# action on its undo stack.

PAINT_MODE_EDIT = "edit"
PAINT_MODE_SET = "set"


class _HyCalGainMap(_HyCalMapBase):
    moduleEditRequested = pyqtSignal(str, float)   # name, current gain
    paintCommitted = pyqtSignal(list)              # [(name, prior_value_or_None), ...]

    def __init__(self, modules: List[Module], parent=None):
        super().__init__(parent, shrink=0.92, margin_top=10,
                         margin_bottom=40, include_lms=True,
                         label_types=AUX_TYPES, show_colorbar=True,
                         min_size=(500, 500))
        self._mod_map: Dict[str, Module] = {m.name: m for m in modules}
        self._gains: Dict[str, float] = {}
        self._overrides: Dict[str, float] = {}

        self._paint_mode: str = PAINT_MODE_EDIT
        self._paint_value: float = 0.0
        self._paint_dragging = False
        self._drag_visited: Set[str] = set()
        self._drag_batch: List[Tuple[str, Optional[float]]] = []

        self.set_modules(modules)
        self.set_range(0.0, 1.0)

    # ---- public API ----

    def set_paint_mode(self, mode: str) -> None:
        if mode not in (PAINT_MODE_EDIT, PAINT_MODE_SET):
            return
        self._paint_mode = mode
        self._paint_dragging = False
        self._drag_visited.clear()
        self._drag_batch = []
        self.setCursor(
            Qt.CursorShape.CrossCursor if mode != PAINT_MODE_EDIT
            else Qt.CursorShape.ArrowCursor)

    def set_paint_value(self, v: float) -> None:
        self._paint_value = float(v)

    def set_gains(self, gains: Dict[str, float],
                  overrides: Optional[Dict[str, float]] = None) -> None:
        self._gains = dict(gains)
        self._overrides = dict(overrides) if overrides is not None else {}
        # Share the dict with the base widget's _values so the colormap
        # picks up our edits without separate set_values calls.
        self.set_values(self._gains)
        self.update()

    @property
    def gains(self) -> Dict[str, float]:
        return self._gains

    @property
    def overrides(self) -> Dict[str, float]:
        return self._overrides

    # ---- painting ----

    def _paint_modules(self, p):
        stops = self.palette_stops()
        no_data = self.NO_DATA_COLOR
        null_color = QColor(THEME.DANGER)
        for name, rect in self._rects.items():
            m = self._mod_map.get(name)
            if m is None or m.crate < 0:
                p.fillRect(rect, no_data)
                continue
            v = self._gains.get(name, 0.0)
            if v == 0.0:
                p.fillRect(rect, null_color)
            else:
                p.fillRect(rect, cmap_qcolor(self.value_to_t(v), stops))

    def _paint_overlays(self, p, w, h):
        # Highlight border around modules whose gain was set via the GUI.
        sel_pen = QPen(QColor(THEME.SELECT_BORDER), 1.5)
        sel_pen.setCosmetic(True)
        p.setPen(sel_pen)
        p.setBrush(Qt.BrushStyle.NoBrush)
        for name in self._overrides:
            rect = self._rects.get(name)
            if rect is not None:
                p.drawRect(rect)
        super()._paint_overlays(p, w, h)

    def _paint_after_colorbar(self, p, w, h):
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Monospace", 9))
        n_masked = sum(1 for v in self._gains.values() if v == 0.0)
        info = f"Edits: {len(self._overrides)}    Masked: {n_masked}"
        if self._paint_mode == PAINT_MODE_SET:
            info += f"    [SET={self._paint_value:.6g}]"
        p.drawText(QRectF(8, h - 18, w - 16, 16),
                   Qt.AlignmentFlag.AlignLeft, info)

    # ---- hit / mouse ----

    def _hit(self, pos) -> Optional[str]:
        for name in self._rect_names_rev:
            if self._rects[name].contains(pos):
                m = self._mod_map.get(name)
                if m and m.crate >= 0:
                    return name
        return None

    def _apply_paint(self, name: str) -> None:
        """Apply the current paint value to ``name``.  Records the prior
        override (or ``None``) onto the active drag batch so the editor
        can undo the action."""
        m = self._mod_map.get(name)
        if not m or m.crate < 0:
            return
        if self._paint_mode != PAINT_MODE_SET:
            return
        v = self._paint_value
        # No-op if the cell already holds this exact override.
        if (name in self._overrides
                and self._overrides[name] == v
                and self._gains.get(name) == v):
            return
        prior = self._overrides.get(name)
        self._drag_batch.append((name, prior))
        self._gains[name] = v
        self._overrides[name] = v

    def _tooltip_text(self, name: str) -> str:
        m = self._mod_map.get(name)
        v = self._gains.get(name, 0.0)
        tip = f"{name}: {v:.6g}"
        if v == 0.0:
            tip += "  [masked]"
        if name in self._overrides:
            tip += "  (edit)"
        if m and m.crate >= 0:
            tip += f"\ncrate={m.crate} slot={m.slot} ch={m.channel}"
        if self._paint_mode == PAINT_MODE_SET:
            tip += f"\n(click/drag to set {self._paint_value:.6g})"
        else:
            tip += "\n(click to edit gain)"
        return tip

    def mousePressEvent(self, event):
        if event.button() != Qt.MouseButton.LeftButton:
            return
        pos = event.position()
        if self._check_inline_range_edit_click(pos):
            return
        if self._cb_rect and self._cb_rect.contains(pos):
            self.cycle_palette()
            return
        found = self._hit(pos)
        if not found:
            return
        if self._paint_mode == PAINT_MODE_EDIT:
            self.moduleEditRequested.emit(found, self._gains.get(found, 0.0))
        else:
            self._paint_dragging = True
            self._drag_visited = {found}
            self._drag_batch = []
            self._apply_paint(found)
            self.update()

    def mouseMoveEvent(self, event):
        pos = event.position()
        found = self._hit(pos)
        if found != self._hovered:
            self._hovered = found
            self.update()
            if found:
                QToolTip.showText(event.globalPosition().toPoint(),
                                  self._tooltip_text(found), self)
                self.moduleHovered.emit(found)
            else:
                QToolTip.hideText()
        if (self._paint_mode != PAINT_MODE_EDIT and self._paint_dragging
                and found and found not in self._drag_visited):
            self._drag_visited.add(found)
            self._apply_paint(found)
            self.update()

    def mouseReleaseEvent(self, event):
        if event.button() != Qt.MouseButton.LeftButton:
            return
        if self._paint_dragging and self._drag_batch:
            self.paintCommitted.emit(list(self._drag_batch))
        self._paint_dragging = False
        self._drag_visited.clear()
        self._drag_batch = []


class _GainEditor(QMainWindow):
    def __init__(self,
                 modules: List[Module],
                 daq: List[Tuple[str, int, int, int]],
                 mod_types: Dict[str, str],
                 db_dir: Path,
                 cal: Optional[Dict[str, float]] = None,
                 cal_path: Optional[Path] = None,
                 pbwo4_gain: float = 1.0,
                 pbglass_gain: float = 1.0,
                 initial_overrides: Optional[Dict[str, float]] = None,
                 output_path: Optional[str] = None):
        super().__init__()
        self.setWindowTitle("FADC Gain Editor")
        self.resize(1600, 900)

        self._daq = daq
        self._mod_types = mod_types
        self._db_dir = db_dir
        self._cal: Dict[str, float] = dict(cal or {})
        self._cal_path: Optional[Path] = cal_path
        self._pbwo4 = pbwo4_gain
        self._pbglass = pbglass_gain
        self._output_path = output_path
        self._mod_map: Dict[str, Module] = {m.name: m for m in modules}

        apply_theme_palette(self)
        # Window-scoped stylesheet so QLabel / QDoubleSpinBox / QLineEdit
        # (which don't reliably pick up the QPalette on Windows native style)
        # render text against the dark surfaces correctly.
        self.setStyleSheet(themed(
            f"QLabel{{color:{THEME.TEXT};background:transparent;}}"
            f"QDoubleSpinBox,QSpinBox,QLineEdit{{background:{THEME.PANEL};"
            f"color:{THEME.TEXT};border:1px solid {THEME.BORDER};"
            f"border-radius:4px;padding:2px 6px;}}"
            f"QGroupBox{{color:{THEME.TEXT};background:transparent;"
            f"border:1px solid {THEME.BORDER};border-radius:6px;"
            f"margin-top:10px;padding-top:10px;}}"
            f"QGroupBox::title{{subcontrol-origin:margin;left:10px;"
            f"padding:0 6px;color:{THEME.TEXT_DIM};}}"
            f"QSplitter::handle{{background:{THEME.BORDER};}}"
        ))

        self._base_gains = self._compute_base_gains()
        overrides = self._diff_overrides(initial_overrides or {})

        # Undo stack: each entry is a batch (list) of (name, prior_value_or_None)
        # tuples.  ``prior_value=None`` means "module had no override before"
        # — undoing reverts it to the base gain.
        self._history: List[List[Tuple[str, Optional[float]]]] = []

        # Colormap range control built later in _build_right_panel; until
        # then _after_change() skips it.
        self._range_ctrl: Optional[ColorRangeControl] = None

        self._map = _HyCalGainMap(modules)
        merged = dict(self._base_gains)
        merged.update(overrides)
        self._map.set_gains(merged, overrides)

        self._build_right_panel()

        self._map.moduleEditRequested.connect(self._on_module_edit_requested)
        self._map.paintCommitted.connect(self._on_paint_committed)
        self._map.moduleHovered.connect(self._on_hover)

        self._status = QLabel("Click or drag modules to apply gain")
        self._status.setStyleSheet(themed(
            f"color:{THEME.TEXT_DIM};font:10pt Monospace;padding:4px;"))

        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.addWidget(self._map)
        splitter.addWidget(self._right)
        splitter.setSizes([900, 700])

        central = QWidget()
        layout = QVBoxLayout(central)
        layout.setContentsMargins(0, 0, 0, 0)
        layout.addWidget(splitter, 1)
        layout.addWidget(self._status)
        self.setCentralWidget(central)

        self._update_cal_label()
        self._refresh_text()

    # ---- helpers ----

    def _compute_base_gains(self) -> Dict[str, float]:
        return {name: resolve_gain(name, self._mod_types.get(name), self._cal,
                                   self._pbwo4, self._pbglass)
                for name, _, _, _ in self._daq}

    def _diff_overrides(self,
                        candidate: Dict[str, float]) -> Dict[str, float]:
        out: Dict[str, float] = {}
        for name, val in candidate.items():
            base = self._base_gains.get(name)
            if base is None or abs(base - val) > 1e-12:
                out[name] = val
        return out

    def _after_change(self) -> None:
        """The map's gains changed: tell the range control (it re-fits if
        the Auto button is pinned) and refresh the .cnf preview."""
        if self._range_ctrl is not None:
            self._range_ctrl.notify_values_changed(self._map.gains)
        self._refresh_text()

    def _commit(self, gains: Dict[str, float],
                overrides: Dict[str, float]) -> None:
        self._map.set_gains(gains, overrides)
        self._after_change()

    def _apply_changes(self, changes: Dict[str, float]) -> int:
        """Set ``{name: gain}`` as GUI edits, recorded as one undo batch.
        Channels already overridden to that gain are skipped.  Returns the
        number of channels changed."""
        gains = dict(self._map.gains)
        overrides = dict(self._map.overrides)
        batch: List[Tuple[str, Optional[float]]] = []
        for name, v in changes.items():
            if gains.get(name) == v and overrides.get(name) == v:
                continue
            batch.append((name, overrides.get(name)))
            gains[name] = v
            overrides[name] = v
        if batch:
            self._history.append(batch)
            self._commit(gains, overrides)
        return len(batch)

    def _rebuild_from_base(self,
                           keep_overrides: bool = True) -> None:
        """Recompute base gains and re-merge with overrides on top."""
        self._base_gains = self._compute_base_gains()
        overrides = self._map.overrides if keep_overrides else {}
        merged = dict(self._base_gains)
        merged.update(overrides)
        self._commit(merged, overrides)

    # ---- right panel ----

    def _build_right_panel(self) -> None:
        self._right = QWidget()
        v = QVBoxLayout(self._right)
        v.setContentsMargins(8, 8, 8, 8)
        v.setSpacing(8)

        # ---- Source group ----
        src_grp = QGroupBox("Source")
        sform = QFormLayout(src_grp)

        cal_row = QHBoxLayout()
        self._cal_label = QLabel("(no calibration loaded)")
        self._cal_label.setStyleSheet(themed(
            f"color:{THEME.TEXT_DIM};font:9pt Monospace;"))
        cal_row.addWidget(self._cal_label, 1)
        btn_load_cal = QPushButton("Load Calibration…")
        btn_load_cal.setStyleSheet(btn_style())
        btn_load_cal.clicked.connect(self._load_calibration_file)
        cal_row.addWidget(btn_load_cal)
        sform.addRow(cal_row)

        self._sb_pbwo4 = QDoubleSpinBox()
        self._sb_pbwo4.setRange(0.0, 100.0)
        self._sb_pbwo4.setDecimals(6)
        self._sb_pbwo4.setValue(self._pbwo4)
        self._sb_pbwo4.valueChanged.connect(self._on_default_changed)
        sform.addRow(QLabel("PbWO4 default:"), self._sb_pbwo4)

        self._sb_pbglass = QDoubleSpinBox()
        self._sb_pbglass.setRange(0.0, 100.0)
        self._sb_pbglass.setDecimals(6)
        self._sb_pbglass.setValue(self._pbglass)
        self._sb_pbglass.valueChanged.connect(self._on_default_changed)
        sform.addRow(QLabel("PbGlass default:"), self._sb_pbglass)

        v.addWidget(src_grp)

        # ---- Color Range group ----
        # auto_fit="minmax_nonzero" ignores zero-valued (masked) channels.
        range_grp = QGroupBox("Color Range")
        rlayout = QHBoxLayout(range_grp)
        rlayout.setContentsMargins(8, 4, 8, 4)
        self._range_ctrl = ColorRangeControl(
            self._map,
            auto_fit="minmax_nonzero",
            orientation="horizontal",
        )
        rlayout.addWidget(self._range_ctrl)
        v.addWidget(range_grp)

        # ---- Edit group ----
        edit_grp = QGroupBox("Edit")
        elayout = QHBoxLayout(edit_grp)
        elayout.setContentsMargins(8, 4, 8, 4)
        elayout.setSpacing(6)

        self._set_value_edit = QLineEdit("0.150000")
        self._set_value_edit.setMaximumWidth(110)
        self._set_value_edit.setValidator(
            QDoubleValidator(0.0, 100.0, 6, self._set_value_edit))
        self._set_value_edit.editingFinished.connect(self._on_set_value_changed)

        self._btn_set = QPushButton("Set")
        self._btn_set.setStyleSheet(
            btn_style(checked_color=THEME.ACCENT_STRONG))
        self._btn_set.setCheckable(True)
        self._btn_set.setToolTip(
            "Toggle set mode — click or drag modules to apply the value "
            "(use 0 to mask)")
        self._btn_set.toggled.connect(self._on_set_toggled)

        self._btn_set_all = QPushButton("Set All")
        self._btn_set_all.setStyleSheet(btn_style())
        self._btn_set_all.setToolTip(
            "Apply the value to every DAQ-mapped channel (use 0 to mask all)")
        self._btn_set_all.clicked.connect(self._on_set_all)

        self._btn_undo = QPushButton("Undo")
        self._btn_undo.setStyleSheet(btn_style())
        self._btn_undo.setToolTip("Revert the most recent edit")
        self._btn_undo.clicked.connect(self._undo)

        self._btn_reset = QPushButton("Reset")
        self._btn_reset.setStyleSheet(btn_style())
        self._btn_reset.setToolTip(
            "Discard all manual edits, revert to loaded base")
        self._btn_reset.clicked.connect(self._reset_overrides)

        elayout.addWidget(self._set_value_edit)
        elayout.addSpacing(6)
        elayout.addWidget(self._btn_set)
        elayout.addWidget(self._btn_set_all)
        elayout.addStretch()
        elayout.addWidget(self._btn_undo)
        elayout.addWidget(self._btn_reset)

        v.addWidget(edit_grp)

        # ---- Safe Cap group ----
        cap_grp = QGroupBox("Safe Cap")
        clayout = QHBoxLayout(cap_grp)
        clayout.setContentsMargins(8, 4, 8, 4)
        clayout.setSpacing(6)

        self._sb_cap_min = QDoubleSpinBox()
        self._sb_cap_min.setRange(0.0, 100.0)
        self._sb_cap_min.setDecimals(6)
        self._sb_cap_min.setValue(0.0)
        self._sb_cap_min.setToolTip(
            "Lower bound (clamped to >= 0; gains are non-negative)")

        self._sb_cap_max = QDoubleSpinBox()
        self._sb_cap_max.setRange(0.0, 100.0)
        self._sb_cap_max.setDecimals(6)
        self._sb_cap_max.setValue(0.15)
        self._sb_cap_max.setToolTip("Upper bound")

        self._btn_cap = QPushButton("Apply Safe Cap")
        self._btn_cap.setStyleSheet(btn_style())
        self._btn_cap.setToolTip(
            "Clamp every DAQ-mapped channel's gain to [Min, Max]")
        self._btn_cap.clicked.connect(self._on_apply_safe_cap)

        clayout.addWidget(QLabel("Min:"))
        clayout.addWidget(self._sb_cap_min)
        clayout.addSpacing(6)
        clayout.addWidget(QLabel("Max:"))
        clayout.addWidget(self._sb_cap_max)
        clayout.addStretch()
        clayout.addWidget(self._btn_cap)

        v.addWidget(cap_grp)

        # ---- File row ----
        btn_load_cnf = QPushButton("Load .cnf…")
        btn_load_cnf.setStyleSheet(btn_style())
        btn_load_cnf.clicked.connect(self._load_cnf_file)

        btn_save = QPushButton("Save .cnf…")
        btn_save.setStyleSheet(btn_style())
        btn_save.clicked.connect(self._save_as)

        row = QHBoxLayout()
        row.addStretch()
        row.addWidget(btn_load_cnf)
        row.addWidget(btn_save)
        v.addLayout(row)

        # ---- Live .cnf preview ----
        self._text = QTextEdit()
        self._text.setReadOnly(True)
        self._text.setStyleSheet(themed(
            f"QTextEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"font:9pt Monospace;border:1px solid {THEME.BORDER};}}"))
        v.addWidget(self._text, 1)

    # ---- handlers ----

    def _update_cal_label(self) -> None:
        if self._cal_path is not None:
            self._cal_label.setText(
                f"{self._cal_path.name}  ({len(self._cal)} entries)")
        elif self._cal:
            self._cal_label.setText(f"(in-memory, {len(self._cal)} entries)")
        else:
            self._cal_label.setText("(no calibration loaded)")

    def _on_set_toggled(self, on: bool) -> None:
        if on:
            v = self._read_set_value()
            if v is None:
                self._status.setText(
                    "Set: enter a numeric gain value first")
                self._btn_set.blockSignals(True)
                self._btn_set.setChecked(False)
                self._btn_set.blockSignals(False)
                return
            self._map.set_paint_value(v)
            self._map.set_paint_mode(PAINT_MODE_SET)
            self._status.setText(
                f"Set mode — click or drag modules to apply gain = {v:.6g}")
        else:
            self._map.set_paint_mode(PAINT_MODE_EDIT)
            self._status.setText(
                "Edit mode — click a module to set its gain")

    def _on_set_value_changed(self) -> None:
        v = self._read_set_value()
        if v is None:
            return
        self._map.set_paint_value(v)
        if self._btn_set.isChecked():
            self._status.setText(f"Set value = {v:.6g}")

    def _read_set_value(self) -> Optional[float]:
        try:
            return float(self._set_value_edit.text())
        except ValueError:
            return None

    def _on_set_all(self) -> None:
        v = self._read_set_value()
        if v is None:
            self._status.setText("Set All: enter a numeric gain value first")
            return
        if QMessageBox.question(
                self, "Set all channels?",
                f"Set every DAQ-mapped channel's gain to {v:.6g}?\n"
                "Use Undo to revert.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
        ) != QMessageBox.StandardButton.Yes:
            return
        n = self._apply_changes(
            {name: v for name, m in self._mod_map.items() if m.crate >= 0})
        self._status.setText(f"Set all to {v:.6g} — {n} channel(s)" if n
                             else "Nothing to change")

    def _on_apply_safe_cap(self) -> None:
        try:
            capped = safe_cap_gains(
                self._map.gains,
                self._sb_cap_min.value(),
                self._sb_cap_max.value())
        except ValueError as exc:
            QMessageBox.warning(self, "Invalid range", str(exc))
            return

        min_cap = max(0.0, self._sb_cap_min.value())
        max_cap = self._sb_cap_max.value()
        if QMessageBox.question(
                self, "Apply safe cap?",
                f"Clamp every DAQ-mapped channel's gain to "
                f"[{min_cap:.6g}, {max_cap:.6g}]?\nUse Undo to revert.",
                QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
                QMessageBox.StandardButton.No
        ) != QMessageBox.StandardButton.Yes:
            return

        gains = self._map.gains
        n = self._apply_changes(
            {name: capped[name] for name, m in self._mod_map.items()
             if m.crate >= 0 and name in capped
             and capped[name] != gains.get(name)})
        if not n:
            self._status.setText("Safe Cap: all channels already in range")
            return
        self._status.setText(
            f"Safe-capped to [{min_cap:.6g}, {max_cap:.6g}] — "
            f"{n} channel(s)")

    def _on_default_changed(self, _) -> None:
        self._pbwo4 = self._sb_pbwo4.value()
        self._pbglass = self._sb_pbglass.value()
        self._rebuild_from_base(keep_overrides=True)

    def _on_module_edit_requested(self, name: str, current: float) -> None:
        """Open a popup dialog to set ``name``'s gain."""
        m = self._mod_map.get(name)
        crate_str = (f"  ({CRATE_NAMES[m.crate]} slot {m.slot} ch {m.channel})"
                     if m and m.crate >= 0 else "")
        new_val, ok = QInputDialog.getDouble(
            self, f"Edit gain — {name}",
            f"{name}{crate_str}\nGain (current: {current:.6g}):",
            current, 0.0, 100.0, 6)
        if not ok:
            return
        self._apply_changes({name: new_val})
        self._status.setText(f"Set {name} = {new_val:.6g}")

    def _on_paint_committed(self,
                            batch: List[Tuple[str, Optional[float]]]) -> None:
        """Drag-paint finished — record the batch onto the undo stack."""
        if not batch:
            return
        self._history.append(list(batch))
        self._after_change()
        self._status.setText(f"Masked {len(batch)} module(s)")

    def _undo(self) -> None:
        if not self._history:
            self._status.setText("Nothing to undo")
            return
        batch = self._history.pop()
        gains = dict(self._map.gains)
        overrides = dict(self._map.overrides)
        for name, prior in batch:
            if prior is None:
                overrides.pop(name, None)
                gains[name] = self._base_gains.get(name, 0.0)
            else:
                overrides[name] = prior
                gains[name] = prior
        self._commit(gains, overrides)
        self._status.setText(f"Undone {len(batch)} edit(s)")

    def _on_hover(self, name: str) -> None:
        m = self._mod_map.get(name)
        if not m or m.crate < 0:
            self._status.setText(f"{name}  (no DAQ mapping)")
            return
        v = self._map.gains.get(name, 0.0)
        flags = []
        if name in self._map.overrides:
            flags.append("edit")
        if v == 0.0:
            flags.append("closed")
        tag = ("  [" + ", ".join(flags) + "]") if flags else ""
        self._status.setText(
            f"{name}  ({CRATE_NAMES[m.crate]} slot {m.slot} ch {m.channel}) "
            f" gain={v:.6g}{tag}")

    def _refresh_text(self) -> None:
        cal_line = (f"# calibration : {self._cal_path}"
                    if self._cal_path
                    else "# calibration : (none)")
        text = render_cnf(
            self._daq, self._map.gains,
            [cal_line,
             f"# edits applied : {len(self._map.overrides)}",
             f"# defaults : PbWO4={self._pbwo4}, PbGlass={self._pbglass}"])
        self._text.setPlainText(text)

    def _reset_overrides(self) -> None:
        self._history.clear()
        self._commit(self._base_gains, {})
        self._status.setText("All edits cleared")

    def _load_calibration_file(self) -> None:
        cal_dir = self._db_dir / "calibration"
        start = str(cal_dir if cal_dir.is_dir() else self._db_dir)
        path, _ = QFileDialog.getOpenFileName(
            self, "Load Calibration JSON", start,
            "JSON (*.json);;All Files (*)")
        if not path:
            return
        try:
            cal = load_calibration(Path(path))
        except Exception as exc:
            self._status.setText(f"Calibration load failed: {exc}")
            return
        self._cal = cal
        self._cal_path = Path(path)
        self._update_cal_label()
        # Drop manual edits — the user just changed the base.
        self._history.clear()
        self._rebuild_from_base(keep_overrides=False)
        self._status.setText(
            f"Loaded calibration {self._cal_path.name}: {len(cal)} entries")

    def _save_as(self) -> None:
        path, _ = QFileDialog.getSaveFileName(
            self, "Save FADC Gain Config",
            self._output_path or "adchycal_gain.cnf",
            "Config Files (*.cnf);;All Files (*)")
        if not path:
            return
        self._output_path = path
        with open(path, "w") as f:
            f.write(self._text.toPlainText())
        self._status.setText(f"Saved to {path}")

    def _load_cnf_file(self) -> None:
        path, _ = QFileDialog.getOpenFileName(
            self, "Load FADC Gain Config", "",
            "Config Files (*.cnf);;All Files (*)")
        if not path:
            return
        with open(path) as f:
            text = f.read()
        loaded = parse_cnf_text(text, self._daq)
        overrides = self._diff_overrides(loaded)
        merged = dict(self._base_gains)
        merged.update(loaded)
        self._history.clear()
        self._commit(merged, overrides)
        self._status.setText(
            f"Loaded {Path(path).name}: {len(loaded)} channels, "
            f"{len(overrides)} differ from base")


# ---- Main ----

def main():
    parser = argparse.ArgumentParser(
        description="Interactive FADC gain config editor (always opens GUI)")
    parser.add_argument("-c", "--calibration",
                        help="Pre-load this calibration JSON in the editor "
                             "(equivalent to clicking Load Calibration…)")
    parser.add_argument("--pbwo4-gain", type=float, default=1.0,
                        help="Initial PbWO4 default gain (default: 1.0)")
    parser.add_argument("--pbglass-gain", type=float, default=1.0,
                        help="Initial PbGlass default gain (default: 1.0)")
    parser.add_argument("-d", "--database",
                        help="Database directory (default: auto-search)")
    parser.add_argument("-o", "--output", default="adchycal_gain.cnf",
                        help="Default save path (default: adchycal_gain.cnf)")
    parser.add_argument("-i", "--input",
                        help="Pre-load this .cnf as starting edits")
    parser.add_argument("--theme", default="dark",
                        choices=available_themes(),
                        help="GUI colour theme (default: dark)")
    args = parser.parse_args()

    db_dir = find_database_dir(args.database)
    print(f"database : {db_dir}")

    modules = load_module_info(db_dir)
    mod_types = {m.name: m.mod_type for m in modules}
    daq = [(m.name, m.crate, m.slot, m.channel)
           for m in modules if m.crate >= 0]
    print(f"modules  : {len(mod_types)}   daq entries: {len(daq)}")

    cal: Dict[str, float] = {}
    cal_path: Optional[Path] = None
    if args.calibration:
        p = Path(args.calibration)
        if not p.is_absolute() and not p.is_file():
            alt = db_dir / "calibration" / p.name
            if alt.is_file():
                p = alt
        if not p.is_file():
            sys.exit(f"error: calibration file not found: {args.calibration}")
        cal = load_calibration(p)
        cal_path = p
        print(f"cal file : {cal_path}  ({len(cal)} entries)")

    initial_overrides: Dict[str, float] = {}
    if args.input:
        input_path = Path(args.input)
        if not input_path.is_file():
            sys.exit(f"error: --input file not found: {input_path}")
        with open(input_path) as f:
            initial_overrides = parse_cnf_text(f.read(), daq)
        print(f"input cnf: {input_path}  ({len(initial_overrides)} channels)")

    set_theme(args.theme)
    app = QApplication.instance() or QApplication(sys.argv)
    win = _GainEditor(modules, daq, mod_types, db_dir,
                      cal=cal, cal_path=cal_path,
                      pbwo4_gain=args.pbwo4_gain,
                      pbglass_gain=args.pbglass_gain,
                      initial_overrides=initial_overrides,
                      output_path=args.output)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
