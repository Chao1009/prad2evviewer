"""
Shared GUI code for the HyCal scan GUIs (``hycal_snake_scan.py`` and
``hycal_gain_equalizer.py``): the ``ScanWindowBase`` main-window shell
and the ``run_scan_gui`` entry point, focus-guarded input widgets,
session log file setup, log line formatting, the position check panel,
encoder-drift monitoring, and the EPICS / scaler bring-up.
"""

from __future__ import annotations

import argparse
import html as html_mod
import math
import os
import sys
from datetime import datetime
from typing import Any, Dict, List, Optional, Sequence, Tuple

from PyQt6.QtWidgets import (
    QApplication, QComboBox, QDoubleSpinBox, QFrame, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox, QPushButton,
    QSpinBox, QSplitter, QTextEdit, QVBoxLayout, QWidget,
)
from PyQt6.QtCore import Qt, QTimer, pyqtSignal
from PyQt6.QtGui import QFont

from scan_utils import (
    C, Module, DARK_QSS, BEAM_CENTER_X, BEAM_CENTER_Y, DEFAULT_DB_PATH,
    PATHS_FILE, load_modules, load_profiles, filter_scan_modules,
    module_to_ptrans, ptrans_to_module,
)
from scan_epics import (
    MotorEPICS, ScalerPVGroup, SimulatedMotorEPICS, SimulatedScalerEPICS,
    epics_move_to,
)
from scan_engine import DEFAULT_VELO_X, DEFAULT_VELO_Y, build_scan_path
from scan_geoview import HyCalScanMapWidget, PALETTES, PALETTE_NAMES


# -- Constants shared by both GUI scripts -----------------------------------

POLL_MS = 200             # main UI poll interval (5 Hz)
SCALER_POLL_MS = 5_000    # default scaler poll interval (5 s)

PROFILE_AUTOGEN = "(autogen)"
PROFILE_NONE = "(none)"

ENCODER_DRIFT_WARN = 0.5   # mm — yellow threshold
ENCODER_DRIFT_ERR  = 1.5   # mm — red threshold


# -- Focus-guarded input widgets --------------------------------------------
#  Wheel events change spin-box / combo-box values only when the widget is
#  focused.  Unfocused wheel events are ignored so they propagate to the
#  enclosing QScrollArea — scrolling the control panel does not flip
#  parameters by accident.

class _NoScrollMixin:
    def __init__(self, *a, **kw):
        super().__init__(*a, **kw)
        self.setFocusPolicy(Qt.FocusPolicy.StrongFocus)

    def wheelEvent(self, event):
        if self.hasFocus():
            super().wheelEvent(event)
        else:
            event.ignore()


class NoScrollSpinBox(_NoScrollMixin, QSpinBox):
    pass


class NoScrollDoubleSpinBox(_NoScrollMixin, QDoubleSpinBox):
    pass


class NoScrollComboBox(_NoScrollMixin, QComboBox):
    pass


# -- Session log file -------------------------------------------------------

def open_session_log(tool_prefix: str, simulation: bool, observer: bool):
    """Create / append today's session log file.

    Returns ``None`` in observer mode (read-only — never writes to disk).
    Otherwise opens ``logs/{SIM_}{tool_prefix}_YYYYMMDD.log`` and writes
    a session-start banner.
    """
    if observer:
        return None
    log_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "logs")
    os.makedirs(log_dir, exist_ok=True)
    prefix = "SIM_" if simulation else ""
    name = datetime.now().strftime(f"{prefix}{tool_prefix}_%Y%m%d.log")
    f = open(os.path.join(log_dir, name), "a")
    f.write(
        "\n" + "=" * 70 + "\n"
        f"=== Session start: {datetime.now().strftime('%Y-%m-%d %H:%M:%S')} ===\n"
        + "=" * 70 + "\n")
    f.flush()
    return f


def format_log_line(msg: str, level: str = "info") -> str:
    """Return the timestamped log line both GUIs use."""
    ts = datetime.now().strftime("%H:%M:%S")
    return f"[{ts}] {level.upper().ljust(5)} {msg}"


def html_log_line(line: str, level: str) -> str:
    """Return the HTML span for inserting a log line into a QTextEdit."""
    colors = {"info": C.TEXT, "warn": C.YELLOW, "error": C.RED}
    c = colors.get(level, C.DIM)
    return (
        f'<span style="color:{c};font-family:Consolas;font-size:13pt;">'
        f'{html_mod.escape(line)}</span>')


def append_log_line(text_edit, line: str, level: str) -> None:
    """Append a colour-formatted log line to a QTextEdit and auto-scroll.

    Uses the "scroll-stick" pattern: if the user is already at the
    bottom of the log (within a few pixels) the view follows the new
    entry; if they have scrolled up to inspect history the view is
    left where it is.
    """
    sb = text_edit.verticalScrollBar()
    at_bottom = sb.value() >= sb.maximum() - 4
    text_edit.append(html_log_line(line, level))
    if at_bottom:
        sb.setValue(sb.maximum())


# -- Position-check panel ---------------------------------------------------

def build_position_check_panel(parent_layout) -> Dict[str, QLabel]:
    """Build the standard "Position Check" group box.

    Adds the group box to ``parent_layout`` and returns a dict with
    keys ``target``, ``actual``, ``diff``, ``drift`` mapping to the
    QLabel widgets so the caller can update them later.
    """
    pe = QGroupBox("Position Check")
    lo = QVBoxLayout(pe)
    lbl_target = QLabel("Target: --"); lo.addWidget(lbl_target)
    lbl_actual = QLabel("Actual: --"); lo.addWidget(lbl_actual)
    lbl_diff = QLabel("Diff:   --")
    lbl_diff.setStyleSheet("font: bold 13pt 'Consolas';")
    lo.addWidget(lbl_diff)
    lbl_drift = QLabel("Drift:   --"); lo.addWidget(lbl_drift)
    parent_layout.addWidget(pe)
    return {"target": lbl_target, "actual": lbl_actual,
            "diff": lbl_diff, "drift": lbl_drift}


def update_position_check(labels: Dict[str, QLabel], ep: Any,
                          target_px: Optional[float],
                          target_py: Optional[float],
                          target_name: str = "",
                          scanning: bool = False,
                          pos_threshold: float = 0.5) -> None:
    """Refresh the target / actual / diff labels from live PVs.

    Adds an ETA suffix to the diff label based on the motor velocities.
    If ``scanning`` is True, the diff is colored red/green against
    ``pos_threshold`` so the operator can see at a glance whether the
    motor has reached the target.  Otherwise it stays dim.
    """
    rx = ep.get("x_rbv", 0.0) or 0.0
    ry = ep.get("y_rbv", 0.0) or 0.0
    labels["actual"].setText(f"Actual: ({rx:.3f}, {ry:.3f})")
    if target_px is None or target_py is None:
        labels["target"].setText("Target: --")
        labels["diff"].setText("Diff:   --")
        labels["diff"].setStyleSheet(f"color: {C.DIM}; font: bold 13pt 'Consolas';")
        return
    err = math.sqrt((rx - target_px) ** 2 + (ry - target_py) ** 2)
    name_html = (f' <b style="color:{C.ACCENT}">{target_name}</b>'
                 if target_name else "")
    labels["target"].setText(
        f"Target: ({target_px:.3f}, {target_py:.3f}){name_html}")
    vx = ep.get("x_velo", DEFAULT_VELO_X) or DEFAULT_VELO_X
    vy = ep.get("y_velo", DEFAULT_VELO_Y) or DEFAULT_VELO_Y
    dx, dy = abs(rx - target_px), abs(ry - target_py)
    eta_sec = max(dx / vx if vx > 0 else 0, dy / vy if vy > 0 else 0)
    if eta_sec >= 60:
        eta_str = f" ({int(eta_sec)//60}m {int(eta_sec)%60}s)"
    elif eta_sec >= 1:
        eta_str = f" ({eta_sec:.0f}s)"
    else:
        eta_str = ""
    labels["diff"].setText(f"Diff:   {err:.3f} mm{eta_str}")
    if scanning:
        fg = C.RED if err > pos_threshold else C.GREEN
    else:
        fg = C.DIM
    labels["diff"].setStyleSheet(f"color: {fg}; font: bold 13pt 'Consolas';")


# -- Encoder drift checker --------------------------------------------------

class EncoderDriftChecker:
    """Tracks motor encoder vs RBV drift for the position-check panel.

    Calibrates the encoder offset on the first call (once both encoders
    and RBVs are available), then reports absolute drift in the supplied
    QLabel, color-coded against ``ENCODER_DRIFT_WARN`` / ``_ERR``.
    """

    def __init__(self) -> None:
        self.offset_x: Optional[float] = None
        self.offset_y: Optional[float] = None

    def update(self, ep: Any, log_fn, drift_label: QLabel) -> None:
        enc_x = ep.get("x_encoder", None)
        enc_y = ep.get("y_encoder", None)
        rbv_x = ep.get("x_rbv", None)
        rbv_y = ep.get("y_rbv", None)
        if enc_x is None or enc_y is None or rbv_x is None or rbv_y is None:
            return
        if self.offset_x is None:
            self.offset_x = enc_x - rbv_x
            self.offset_y = enc_y - rbv_y
            log_fn(f"Encoder calibrated: offset "
                   f"X={self.offset_x:.4f} Y={self.offset_y:.4f}")
            return
        dx = abs((enc_x - self.offset_x) - rbv_x)
        dy = abs((enc_y - self.offset_y) - rbv_y)
        fx = (C.RED if dx > ENCODER_DRIFT_ERR
              else C.YELLOW if dx > ENCODER_DRIFT_WARN
              else C.GREEN)
        fy = (C.RED if dy > ENCODER_DRIFT_ERR
              else C.YELLOW if dy > ENCODER_DRIFT_WARN
              else C.GREEN)
        drift_label.setText(
            f'Drift:   X <span style="color:{fx}">{dx:.4f}</span>  '
            f'Y <span style="color:{fy}">{dy:.4f}</span>')


# -- EPICS bring-up ---------------------------------------------------------

def setup_motor_epics(observer: bool, simulation: bool):
    """Create and connect the appropriate motor EPICS group.

    Echoes connection counts and disconnected PVs to stdout in non-sim
    mode so the operator sees them at startup.
    """
    if observer:
        ep = MotorEPICS(writable=False)
    elif simulation:
        ep = SimulatedMotorEPICS()
    else:
        ep = MotorEPICS(writable=True)
    n_ok, n_total = ep.connect()
    if not simulation:
        print(f"EPICS: {n_ok}/{n_total} PVs connected")
        for pv in ep.disconnected_pvs():
            print(f"  NOT connected: {pv}")
    return ep


def setup_scaler_epics(simulation: bool, all_modules: List[Module]):
    """Create and connect the appropriate scaler EPICS group."""
    if simulation:
        ep = SimulatedScalerEPICS(all_modules)
    else:
        ep = ScalerPVGroup(all_modules)
    s_ok, s_total = ep.connect()
    if not simulation:
        print(f"Scalers: {s_ok}/{s_total} PVs connected")
    return ep


# -- Main window shell ------------------------------------------------------

class ScanWindowBase(QMainWindow):
    """Main-window shell shared by the snake scan and the gain equalizer.

    Layout::

        [TOP BAR: title | mode | ===BEAM=== | state             ]
        [LEFT half                 | RIGHT half (subclass)       ]
        [  HyCal geo view          |                             ]
        [  legend row              |                             ]
        [  scaler controls         |                             ]

    The base owns the common state (EPICS groups, modules, profiles, the
    scan path and start selection, the move target, the session log), the
    top bar, the map pane with its scaler controls, profile / path
    handling, logging, and the beam and encoder displays.

    Subclasses set ``TITLE``, ``LOG_PREFIX``, ``WINDOW_H``, ``LEGEND`` and
    ``NONE_MSG`` and implement

    * ``_buildRightPane() -> QWidget`` -- the right half.  It must create
      ``_profile_combo`` (activated -> ``_onPathProfileChanged``),
      ``_lg_spin`` (valueChanged -> ``_onLgLayersChanged``),
      ``_start_combo`` (activated -> ``_onStartSelected``), ``_count_spin``
      (valueChanged -> ``_drawPathPreview``) and ``_beam_thresh_spin``,
      and place ``_buildPositionCheck(layout)`` and ``_buildLogGroup()``.
    * ``_disableControls()`` -- observer / PV-error lockout.
    * ``_poll()`` -- the 5 Hz refresh.

    Optional hooks: ``_initState``, ``_onPathSet``, ``_onModuleClicked``,
    ``_beamTripped`` and ``_stopEngine``.
    """

    _logSignal = pyqtSignal(str, str)

    TITLE = ""                    # window title; upper-cased in the top bar
    LOG_PREFIX = ""               # session log file prefix
    WINDOW_H = 900
    LEGEND: Sequence[Tuple[str, str]] = ()   # (label, colour) swatches
    NONE_MSG = "Path: none"       # logged when the (none) profile is picked

    def __init__(self, motor_ep, scaler_ep, simulation, all_modules,
                 profiles=None, observer=False):
        super().__init__()
        self.ep = motor_ep
        self.scaler_ep = scaler_ep
        self.simulation = simulation
        self.observer = observer
        self.all_modules = all_modules
        self._profiles = profiles or {}
        self._active_profile = PROFILE_NONE
        self._lg_layers = 0
        self._mod_by_name = {m.name: m for m in all_modules}
        self._log_file = open_session_log(self.LOG_PREFIX, simulation, observer)

        self.scan_modules: List[Module] = []
        self._scan_names: set = set()
        self._scan_name_to_idx: Dict[str, int] = {}
        self._selected_start_idx = 0
        self._selected_mod_name: Optional[str] = None
        self._encoder_checker = EncoderDriftChecker()

        # target position — set once when a move is commanded
        self._target_px: Optional[float] = None
        self._target_py: Optional[float] = None
        self._target_name: str = ""

        self._initState()
        self._logSignal.connect(self._appendLog)
        self._buildUI()

        if self.observer:
            self._disableControls()
        if not self.simulation and not self.observer:
            disc = self.ep.disconnected_pvs()
            if disc:
                self._disableControls()
                QMessageBox.critical(self, "PV Connection Error",
                    "Not connected:\n" + "\n".join(f"  {p}" for p in disc))

        self._timer = QTimer(self)
        self._timer.timeout.connect(self._poll)
        self._timer.start(POLL_MS)

        self._scaler_timer = QTimer(self)
        self._scaler_timer.timeout.connect(self._pollScalers)
        self._scaler_timer.start(SCALER_POLL_MS)
        self._pollScalers()

    # -- subclass hooks -------------------------------------------------------

    def _initState(self):
        """Create subclass state that ``_buildUI`` relies on (e.g. the engine)."""

    def _onPathSet(self, path):
        """Called by ``_setPath`` right after the new path is stored."""

    def _onModuleClicked(self, name):
        """Called after a map click (other than a deselect) updated the selection."""

    def _beamTripped(self) -> bool:
        """True while the engine is waiting for the beam to recover."""
        return False

    def _stopEngine(self):
        """Called on window close, before the session log is closed."""

    # -- layout ---------------------------------------------------------------

    def _buildUI(self):
        if self.observer:       suffix = "  [OBSERVER]"
        elif self.simulation:   suffix = "  [SIMULATION]"
        else:                   suffix = "  [EXPERT OPERATOR]"
        self.setWindowTitle(self.TITLE + suffix)
        self.setStyleSheet(DARK_QSS)
        self.resize(1600, self.WINDOW_H)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(0, 0, 0, 0)
        root.setSpacing(0)
        root.addWidget(self._buildTopBar())

        body_splitter = QSplitter(Qt.Orientation.Horizontal)
        body_splitter.setContentsMargins(6, 4, 6, 6)
        body_splitter.addWidget(self._buildMapPane())
        body_splitter.addWidget(self._buildRightPane())
        body_splitter.setStretchFactor(0, 1)  # left half
        body_splitter.setStretchFactor(1, 1)  # right half
        root.addWidget(body_splitter, stretch=1)
        self._updateCanvasLabel()

    @staticmethod
    def _beam_ss(fg, pt):
        return (f"color: {fg}; font: bold {pt}pt 'Consolas'; "
                f"background: transparent; border: none;")

    def _buildTopBar(self) -> QWidget:
        top = QWidget()
        top.setFixedHeight(48)
        top.setStyleSheet("background: #0d1520;")
        tl = QHBoxLayout(top)
        tl.setContentsMargins(12, 0, 12, 0)

        lbl = QLabel(self.TITLE.upper())
        lbl.setStyleSheet(f"color: {C.GREEN}; font: bold 17pt 'Consolas'; background: transparent;")
        tl.addWidget(lbl)

        if self.observer:       mt, mf = "OBSERVER", C.ORANGE
        elif self.simulation:   mt, mf = "SIMULATION", C.YELLOW
        else:                   mt, mf = "EXPERT", C.GREEN
        lbl_mode = QLabel(mt)
        lbl_mode.setStyleSheet(f"color: {mf}; font: bold 13pt 'Consolas'; background: transparent;")
        tl.addWidget(lbl_mode)
        tl.addSpacing(16)

        # prominent beam current; the children need border:none because
        # the QFrame rule also matches QLabel
        beam_frame = QFrame()
        beam_frame.setStyleSheet(
            "QFrame { background: #161b22; border: 1px solid #30363d; border-radius: 4px; }")
        beam_frame.setFixedHeight(36)
        bf_layout = QHBoxLayout(beam_frame)
        bf_layout.setContentsMargins(10, 0, 10, 0)
        bf_layout.setSpacing(6)
        beam_icon = QLabel("BEAM")
        beam_icon.setStyleSheet(self._beam_ss(C.DIM, 12))
        bf_layout.addWidget(beam_icon)
        self._lbl_beam_val = QLabel("-- nA")
        self._lbl_beam_val.setStyleSheet(self._beam_ss(C.GREEN, 18))
        self._lbl_beam_val.setMinimumWidth(140)
        bf_layout.addWidget(self._lbl_beam_val)
        self._lbl_beam_status = QLabel("")
        self._lbl_beam_status.setStyleSheet(self._beam_ss("transparent", 13))
        bf_layout.addWidget(self._lbl_beam_status)
        tl.addWidget(beam_frame)

        tl.addStretch()

        self._lbl_state = QLabel("IDLE")
        self._lbl_state.setStyleSheet(
            f"color: {C.DIM}; font: bold 15pt 'Consolas'; background: transparent;")
        tl.addWidget(self._lbl_state)
        return top

    def _buildMapPane(self) -> QWidget:
        left = QWidget()
        left_lo = QVBoxLayout(left)
        left_lo.setContentsMargins(0, 0, 0, 0)
        left_lo.setSpacing(2)

        self._canvas_label = QLabel()
        self._canvas_label.setStyleSheet(f"color: {C.ACCENT}; font: bold 13pt 'Consolas';")
        left_lo.addWidget(self._canvas_label)

        self._map = HyCalScanMapWidget(self.all_modules)
        self._map.moduleClicked.connect(self._onCanvasClick)
        left_lo.addWidget(self._map, stretch=1)

        # reset button overlaid at the bottom-right corner (see eventFilter)
        self._btn_reset_view = QPushButton("Reset", self._map)
        self._btn_reset_view.setFixedSize(56, 28)
        self._btn_reset_view.setStyleSheet(
            f"QPushButton{{background:rgba(22,27,34,220);color:{C.DIM};"
            f"border:1px solid #30363d;border-radius:2px;padding:0;"
            f"font:12pt Consolas;}}"
            f"QPushButton:hover{{color:{C.TEXT};border-color:{C.ACCENT};}}")
        self._btn_reset_view.clicked.connect(self._map.resetView)
        self._map.installEventFilter(self)

        leg = QHBoxLayout()
        leg.setSpacing(4); leg.setContentsMargins(0, 0, 0, 0)
        for label, colour in self.LEGEND:
            sw = QLabel(); sw.setFixedSize(10, 10)
            sw.setStyleSheet(f"background: {colour}; border: none;")
            leg.addWidget(sw)
            ll = QLabel(label); ll.setStyleSheet(f"color: {C.DIM}; font: 12pt 'Consolas';")
            leg.addWidget(ll)
        leg.addStretch()
        left_lo.addLayout(leg)

        sc_row = QHBoxLayout()
        sc_row.setSpacing(4); sc_row.setContentsMargins(0, 2, 0, 0)

        self._btn_scaler_toggle = QPushButton("Scalers: ON")
        self._btn_scaler_toggle.setStyleSheet(self._small_btn_ss(C.GREEN))
        self._btn_scaler_toggle.setFixedHeight(28)
        self._btn_scaler_toggle.clicked.connect(self._toggleScaler)
        sc_row.addWidget(self._btn_scaler_toggle)

        self._btn_scaler_auto = QPushButton("Auto")
        self._btn_scaler_auto.setFixedHeight(28)
        self._btn_scaler_auto.clicked.connect(self._toggleScalerAuto)
        self._scaler_auto_on = True
        self._updateScalerAutoBtn()
        sc_row.addWidget(self._btn_scaler_auto)

        self._scaler_min_edit = self._small_edit("0")
        sc_row.addWidget(self._scaler_min_edit)
        sc_row.addWidget(QLabel("-"))
        self._scaler_max_edit = self._small_edit("1000")
        sc_row.addWidget(self._scaler_max_edit)

        btn_apply = QPushButton("Apply"); btn_apply.setFixedHeight(28)
        btn_apply.clicked.connect(self._applyScalerRange)
        sc_row.addWidget(btn_apply)

        self._btn_scaler_log = QPushButton("Log: OFF"); self._btn_scaler_log.setFixedHeight(28)
        self._btn_scaler_log.setStyleSheet(self._small_btn_ss(C.DIM))
        self._btn_scaler_log.clicked.connect(self._toggleScalerLog)
        sc_row.addWidget(self._btn_scaler_log)

        self._btn_palette = QPushButton(); self._btn_palette.setFixedSize(90, 28)
        self._btn_palette.setToolTip("Click to cycle colour palette")
        self._btn_palette.clicked.connect(self._cycleScalerPalette)
        self._updatePaletteBtn()
        sc_row.addWidget(self._btn_palette)

        sc_row.addStretch()
        left_lo.addLayout(sc_row)
        return left

    def _buildPositionCheck(self, parent):
        self._pos_labels = build_position_check_panel(parent)

    def _buildLogGroup(self) -> QGroupBox:
        log_group = QGroupBox("Event Log")
        log_layout = QVBoxLayout(log_group)
        log_layout.setContentsMargins(4, 4, 4, 4)
        self._log_text = QTextEdit()
        self._log_text.setReadOnly(True)
        log_layout.addWidget(self._log_text)
        return log_group

    # -- scaler controls ------------------------------------------------------

    @staticmethod
    def _small_btn_ss(fg):
        return (f"QPushButton{{background:#21262d;color:{fg};"
                f"border:1px solid #30363d;padding:1px 8px;"
                f"font:bold 12pt Consolas;border-radius:2px;}}"
                f"QPushButton:hover{{background:#30363d;}}")

    def _small_edit(self, text):
        e = QLineEdit(text); e.setFixedWidth(50); e.setFixedHeight(28)
        e.setFont(QFont("Consolas", 10))
        e.setStyleSheet("QLineEdit{background:#161b22;color:#c9d1d9;"
                        "border:1px solid #30363d;border-radius:2px;padding:1px 4px;}")
        e.returnPressed.connect(self._applyScalerRange)
        return e

    def _toggleScaler(self):
        on = not self._map._scaler_enabled
        self._map.setScalerEnabled(on)
        self._btn_scaler_toggle.setText("Scalers: ON" if on else "Scalers: OFF")
        self._btn_scaler_toggle.setStyleSheet(self._small_btn_ss(C.GREEN if on else C.RED))

    def _toggleScalerAuto(self):
        self._scaler_auto_on = not self._scaler_auto_on
        self._map.setScalerAutoRange(self._scaler_auto_on)
        self._updateScalerAutoBtn()
        if self._scaler_auto_on:
            vmin, vmax = self._map.scalerRange()
            self._scaler_min_edit.setText(f"{vmin:.0f}")
            self._scaler_max_edit.setText(f"{vmax:.0f}")

    def _updateScalerAutoBtn(self):
        if self._scaler_auto_on:
            self._btn_scaler_auto.setStyleSheet(
                "QPushButton{background:#d29922;color:#0d1117;"
                "border:1px solid #d29922;padding:1px 8px;"
                "font:bold 12pt Consolas;border-radius:2px;}"
                "QPushButton:hover{background:#e0a82b;}")
        else:
            self._btn_scaler_auto.setStyleSheet(self._small_btn_ss(C.YELLOW))

    def _applyScalerRange(self):
        try:
            vmin = float(self._scaler_min_edit.text())
            vmax = float(self._scaler_max_edit.text())
            if vmin < vmax:
                self._map.setScalerRange(vmin, vmax)
                self._scaler_auto_on = False
                self._map.setScalerAutoRange(False)
                self._updateScalerAutoBtn()
        except ValueError:
            pass

    def _toggleScalerLog(self):
        on = not self._map.is_log_scale()
        self._map.setScalerLogScale(on)
        self._btn_scaler_log.setText("Log: ON" if on else "Log: OFF")
        self._btn_scaler_log.setStyleSheet(self._small_btn_ss(C.ACCENT if on else C.DIM))

    def _cycleScalerPalette(self):
        self._map.cyclePalette(); self._updatePaletteBtn()

    def _updatePaletteBtn(self):
        """Paint the palette button with a gradient of the current palette."""
        idx = self._map._palette_idx
        stops = list(PALETTES.values())[idx]
        parts = [f"stop:{t:.2f} rgb({r},{g},{b})" for t, (r, g, b) in stops]
        self._btn_palette.setStyleSheet(
            f"QPushButton{{background:qlineargradient(x1:0,y1:0,x2:1,y2:0,{','.join(parts)});"
            f"border:1px solid #30363d;border-radius:2px;color:#c9d1d9;"
            f"font:bold 11pt Consolas;padding:0 4px;}}"
            f"QPushButton:hover{{border-color:#58a6ff;}}")
        self._btn_palette.setText(PALETTE_NAMES[idx])

    def _pollScalers(self):
        vals = self.scaler_ep.get_all()
        if vals:
            self._map.setScalerValues(vals)
            if self._scaler_auto_on:
                vmin, vmax = self._map.scalerRange()
                self._scaler_min_edit.setText(f"{vmin:.0f}")
                self._scaler_max_edit.setText(f"{vmax:.0f}")

    # -- path management ------------------------------------------------------

    def _onStartSelected(self, _):
        name = self._start_combo.currentText()
        for i, m in enumerate(self.scan_modules):
            if m.name == name:
                self._selected_start_idx = i
                self._drawPathPreview()
                self._updateCanvasLabel()
                break

    def _onPathProfileChanged(self, _):
        name = self._profile_combo.currentText()
        if name == self._active_profile: return
        self._active_profile = name
        if name == PROFILE_AUTOGEN:
            self._lg_spin.setEnabled(True); self._onLgLayersChanged(force=True); return
        self._lg_spin.setEnabled(False)
        if name == PROFILE_NONE:
            self._setPath([])
            self._log(self.NONE_MSG)
            return
        path_mods = [self._mod_by_name[n] for n in self._profiles.get(name, [])
                     if n in self._mod_by_name]
        if not path_mods:
            self._log(f"Profile '{name}' empty", level="error"); return
        self._setPath(path_mods)
        self._log(f"Path profile: {name} ({len(path_mods)} modules)")

    def _onLgLayersChanged(self, value=0, force=False):
        if self._active_profile != PROFILE_AUTOGEN: return
        nl = self._lg_spin.value()
        if nl == self._lg_layers and not force: return
        self._lg_layers = nl
        # generate the snake path once at autogen — from now on the order is fixed
        path, n_dropped = build_scan_path(filter_scan_modules(self.all_modules, nl))
        if n_dropped:
            self._log(f"WARNING: {n_dropped} modules left out of the path", level="warn")
        self._setPath(path)
        np_ = sum(1 for m in path if m.mod_type == "PbWO4")
        ng = sum(1 for m in path if m.mod_type == "PbGlass")
        self._log(f"LG layers: {nl} ({np_} PbWO4 + {ng} PbGlass = {len(path)})")

    def _setPath(self, path):
        """Set the scan path. ``path`` is the final ordered list of modules."""
        self.scan_modules = path
        self._scan_names = {m.name for m in path}
        self._scan_name_to_idx = {m.name: i for i, m in enumerate(path)}
        self._onPathSet(path)
        self._selected_start_idx = 0
        ns = [m.name for m in path]
        self._start_combo.clear(); self._start_combo.addItems(ns)
        self._count_spin.setMaximum(len(ns)); self._count_spin.setValue(0)
        if not path:
            self._map.setPathPreview([]); self._map.setDashPreview([])
            self._map.setHighlight(None); self._selected_mod_name = None
        self._updateCanvasLabel()

    # -- canvas ---------------------------------------------------------------

    def _updateCanvasLabel(self):
        path = self.scan_modules
        n_pwo4 = sum(1 for m in path if m.mod_type == "PbWO4")
        n_lg = sum(1 for m in path if m.mod_type == "PbGlass")
        base = f"Scan Path: {n_pwo4} PbWO4 + {n_lg} LG" if n_lg else f"Scan Path: {n_pwo4} PbWO4"
        if not path:
            base = "Scan Path: none"
        elif 0 <= self._selected_start_idx < len(path):
            base += f"  start: {path[self._selected_start_idx].name}"
        self._canvas_label.setText(f" {base} ")

    def _drawPathPreview(self):
        path = self.scan_modules
        s = self._selected_start_idx
        if s >= len(path): self._map.setPathPreview([]); return
        c = self._count_spin.value()
        e = min(s + c, len(path)) if c > 0 else len(path)
        self._map.setPathPreview([self._map.modCenter(path[i]) for i in range(s, e)])

    def _refreshMap(self, colors: Dict[str, str], ahead=None):
        """Push module colours, the path preview and the beam marker to the map.

        ``ahead`` lists the modules a running scan has still to visit,
        drawn as a dashed line; ``None`` shows the solid preview of the
        planned path instead.
        """
        self._map.setModuleColors(colors)
        if ahead is None:
            self._drawPathPreview(); self._map.setDashPreview([])
        else:
            self._map.setPathPreview([])
            self._map.setDashPreview([self._map.modCenter(m) for m in ahead])
        rx, ry = self.ep.get("x_rbv", BEAM_CENTER_X), self.ep.get("y_rbv", BEAM_CENTER_Y)
        self._map.setMarkerPosition(*ptrans_to_module(rx, ry))
        self._map.update()

    def _onCanvasClick(self, name):
        if self._selected_mod_name == name:
            self._selected_mod_name = None; self._map.setHighlight(None)
            self._updateCanvasLabel(); return
        self._selected_mod_name = name
        if name in self._scan_name_to_idx:
            self._selected_start_idx = self._scan_name_to_idx[name]
            idx = self._start_combo.findText(name)
            if idx >= 0: self._start_combo.setCurrentIndex(idx)
            self._drawPathPreview()
        self._map.setHighlight(name); self._updateCanvasLabel()
        self._onModuleClicked(name)

    def eventFilter(self, obj, event):
        # keep the Reset overlay in the map's bottom-right corner
        if obj is self._map and event.type() == event.Type.Resize:
            btn = self._btn_reset_view
            btn.move(self._map.width() - btn.width() - 2,
                     self._map.height() - btn.height() - 2)
        return super().eventFilter(obj, event)

    # -- motion ---------------------------------------------------------------

    def _setTarget(self, px, py, name=""):
        self._target_px = px
        self._target_py = py
        self._target_name = name

    def _moveToModule(self, mod: Module):
        """Direct move that centres the beam on *mod*."""
        px, py = module_to_ptrans(mod.x, mod.y)
        self._log(f"Direct move to {mod.name}  ptrans({px:.3f}, {py:.3f})")
        if epics_move_to(self.ep, px, py):
            self._setTarget(px, py, mod.name)
        else:
            self._log(f"BLOCKED: ptrans({px:.3f}, {py:.3f}) outside limits", level="error")

    # -- logging --------------------------------------------------------------

    def _log(self, msg, level="info"):
        """Thread-safe: engines call this from their worker threads."""
        line = format_log_line(msg, level)
        if self._log_file and not self._log_file.closed:
            self._log_file.write(line + "\n"); self._log_file.flush()
        self._logSignal.emit(line, level)

    def _appendLog(self, line, level):
        append_log_line(self._log_text, line, level)

    # -- beam / encoder displays ----------------------------------------------

    def _updateBeamDisplay(self):
        bc = self.ep.get("beam_cur", None)
        if bc is None:
            self._lbl_beam_val.setText("-- nA")
            self._lbl_beam_val.setStyleSheet(self._beam_ss(C.DIM, 18))
            self._lbl_beam_status.setText("")
            return
        thresh = self._beam_thresh_spin.value()
        if self._beamTripped():
            fg, status = C.RED, "TRIP"
        elif thresh > 0 and bc < thresh:
            fg, status = C.YELLOW, "LOW"
        else:
            fg, status = C.GREEN, ""
        self._lbl_beam_val.setText(f"{bc:.2f} nA")
        self._lbl_beam_val.setStyleSheet(self._beam_ss(fg, 18))
        self._lbl_beam_status.setText(status)
        if status:
            self._lbl_beam_status.setStyleSheet(self._beam_ss(fg, 14))

    def _checkEncoder(self):
        self._encoder_checker.update(self.ep, self._log, self._pos_labels["drift"])

    def closeEvent(self, e):
        self._timer.stop()
        self._scaler_timer.stop()
        self._stopEngine()
        if self._log_file and not self._log_file.closed:
            self._log_file.close()
        self._log_file = None
        super().closeEvent(e)


# -- Entry point ------------------------------------------------------------

def run_scan_gui(window_cls) -> None:
    """Parse the common command line, bring up EPICS, and run *window_cls*.

    Modes: simulation (default), ``--expert`` (full control) and
    ``--observer`` (real reads, no writes).
    """
    parser = argparse.ArgumentParser(description=window_cls.TITLE)
    parser.add_argument("--expert", action="store_true")
    parser.add_argument("--observer", action="store_true")
    parser.add_argument("--database", default=DEFAULT_DB_PATH)
    parser.add_argument("--paths", default=PATHS_FILE)
    args = parser.parse_args()

    all_modules = load_modules(args.database)
    by_type: Dict[str, int] = {}
    for m in all_modules:
        by_type[m.mod_type] = by_type.get(m.mod_type, 0) + 1
    print(f"Loaded {len(all_modules)} modules from {args.database}")
    for t, n in sorted(by_type.items()):
        print(f"  {t}: {n}")

    profiles = load_profiles(args.paths)
    if profiles:
        print(f"Loaded {len(profiles)} path profiles")

    observer = args.observer
    simulation = not args.expert and not observer

    motor_ep = setup_motor_epics(observer, simulation)
    scaler_ep = setup_scaler_epics(simulation, all_modules)

    app = QApplication(sys.argv)
    win = window_cls(motor_ep, scaler_ep, simulation, all_modules,
                     profiles, observer=observer)
    win.show()
    sys.exit(app.exec())
