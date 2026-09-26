#!/usr/bin/env python3
"""
HyCal Snake Scan -- Module Scanner (PyQt6)
==========================================
PyQt6 GUI that drives the HyCal transporter in a snake pattern so the
beam centres on each scanned module, dwells for a configurable time,
then advances to the next module.

Includes a live FADC scaler overlay so the beam spot is visible as a
hot region on the HyCal map.

Usage
-----
    python hycal_snake_scan.py                          # simulation
    python hycal_snake_scan.py --expert                  # expert operator
    python hycal_snake_scan.py --observer                # read-only monitor

--database (hycal_map.json) and --paths (scan path profiles JSON) override
the default input files.

Coordinate system
-----------------
    ptrans_x, ptrans_y = (-126.75, 10.11)  -->  beam at HyCal centre (0,0)
    ptrans_x = BEAM_CENTER_X + module_x
    ptrans_y = BEAM_CENTER_Y - module_y

Writable PVs (the ONLY PVs this tool writes to):
    ptrans_x.VAL / ptrans_y.VAL    -- absolute set-point
    ptrans_x.SPMG / ptrans_y.SPMG  -- motor mode  Stop(0) Pause(1) Move(2) Go(3)

Requirements
------------
    Python 3.8+, PyQt6
    pyepics  (only for --expert / --observer mode)
"""

from __future__ import annotations

from typing import Dict, Optional

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGridLayout, QGroupBox, QPushButton,
    QLabel, QComboBox, QSpinBox, QDoubleSpinBox, QProgressBar, QMessageBox,
    QSplitter, QFrame, QDialog, QScrollArea, QSlider,
)
from PyQt6.QtCore import Qt

from scan_utils import (
    C, Module, module_to_ptrans, ptrans_in_limits, DARK_QSS,
    BEAM_CENTER_X, BEAM_CENTER_Y,
)
from scan_epics import (
    SPMG, SPMG_LABELS, epics_move_to, epics_stop,
)
from scan_engine import (
    ScanState, ScanEngine, estimate_scan_time,
    DEFAULT_DWELL, DEFAULT_POS_THRESHOLD, DEFAULT_BEAM_THRESHOLD,
    DEFAULT_VELO_X, DEFAULT_VELO_Y, MAX_LG_LAYERS,
)
from scan_gui_common import (
    ScanWindowBase, run_scan_gui, PROFILE_AUTOGEN, PROFILE_NONE,
    update_position_check,
)


class ModuleInfoDialog(QDialog):
    """Pop-up showing module details with a Move To button.

    Call :meth:`setModule` to refresh the content for a different module
    without closing and re-opening the dialog.
    """

    _FIELDS = ("Scaler", "Name", "Type", "Sector", "Row/Col", "Size", "HyCal", "Ptrans", "In limits")

    def __init__(self, parent=None):
        super().__init__(parent)
        self._mod: Optional[Module] = None
        self.setStyleSheet(DARK_QSS)
        self.setFixedWidth(360)

        lo = QVBoxLayout(self)

        grid = QGridLayout()
        grid.setSpacing(4)
        self._value_labels: Dict[str, QLabel] = {}
        for r, label in enumerate(self._FIELDS):
            lk = QLabel(f"{label}:")
            lk.setStyleSheet(f"color: {C.DIM}; font: 13pt 'Consolas';")
            lk.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
            grid.addWidget(lk, r, 0)
            lv = QLabel("--")
            grid.addWidget(lv, r, 1)
            self._value_labels[label] = lv
        lo.addLayout(grid)

        lo.addSpacing(8)

        btn_row = QHBoxLayout()
        self._btn_move = QPushButton("Move To")
        self._btn_move.setProperty("cssClass", "accent")
        self._btn_move.clicked.connect(self._doMove)
        btn_row.addWidget(self._btn_move)
        btn_close = QPushButton("Close")
        btn_close.clicked.connect(self.close)
        btn_row.addWidget(btn_close)
        lo.addLayout(btn_row)

    def setModule(self, mod: Module, scaler_value: Optional[float] = None):
        self._mod = mod
        px, py = module_to_ptrans(mod.x, mod.y)
        in_limits = ptrans_in_limits(px, py)
        vals = {
            "Scaler":    f"{scaler_value:.1f}" if scaler_value is not None else "--",
            "Name":      mod.name,
            "Type":      mod.mod_type,
            "Sector":    mod.sector or "--",
            "Row/Col":   f"{mod.row} / {mod.col}" if mod.row else "--",
            "Size":      f"{mod.sx:.2f} x {mod.sy:.2f} mm",
            "HyCal":     f"({mod.x:.2f}, {mod.y:.2f}) mm",
            "Ptrans":    f"({px:.2f}, {py:.2f}) mm",
            "In limits": "Yes" if in_limits else "No",
        }
        for label, lv in self._value_labels.items():
            lv.setText(vals.get(label, "--"))
            if label == "In limits" and not in_limits:
                lv.setStyleSheet(f"color: {C.RED}; font: bold 13pt 'Consolas';")
            elif label == "Scaler" and scaler_value is not None:
                lv.setStyleSheet(f"color: {C.GREEN}; font: bold 13pt 'Consolas';")
            else:
                lv.setStyleSheet("")
        self._btn_move.setText(f"Move To {mod.name}")
        self._btn_move.setEnabled(in_limits)
        self.setWindowTitle(f"Module {mod.name}")

    def _doMove(self):
        if self._mod: self.parent()._moveToModule(self._mod)


class SnakeScanWindow(ScanWindowBase):
    TITLE = "HyCal Snake Scan"
    LOG_PREFIX = "snake_scan"
    WINDOW_H = 900
    LEGEND = [("Todo", C.MOD_TODO), ("Skipped", C.MOD_SKIPPED),
              ("Moving", C.MOD_CURRENT), ("Dwell", C.MOD_DWELL),
              ("Done", C.MOD_DONE), ("Error", C.MOD_ERROR),
              ("Start", C.MOD_SELECTED), ("PbGlass", C.MOD_GLASS)]
    NONE_MSG = "Path: none (direct control only)"

    def _initState(self):
        self.engine = ScanEngine(self.ep, self.scan_modules, self._log)
        self._mod_dlg: Optional[ModuleInfoDialog] = None
        self._status_labels: Dict[str, QLabel] = {}
        self._last_scan_idx: int = -1  # track scan engine moves

    def _onPathSet(self, path):
        self.engine = ScanEngine(self.ep, path, self._log)

    def _beamTripped(self):
        return self.engine.beam_tripped

    # -- layout --------------------------------------------------------------

    def _buildRightPane(self):
        # controls (top) + event log (bottom)
        right_splitter = QSplitter(Qt.Orientation.Vertical)

        # top area: two columns — left: scan/direct control, right: status
        ctrl_columns = QWidget()
        ctrl_cols_lo = QHBoxLayout(ctrl_columns)
        ctrl_cols_lo.setContentsMargins(0, 0, 0, 0)
        ctrl_cols_lo.setSpacing(4)

        # left column: scan control + direct control
        left_col_scroll = QScrollArea()
        left_col_scroll.setWidgetResizable(True)
        left_col_scroll.setFrameShape(QFrame.Shape.NoFrame)
        left_col_widget = QWidget()
        left_col_lo = QVBoxLayout(left_col_widget)
        left_col_lo.setSpacing(4)
        left_col_lo.setContentsMargins(0, 0, 0, 0)
        self._buildScanControl(left_col_lo)
        self._buildDirectControl(left_col_lo)
        left_col_lo.addStretch()
        left_col_scroll.setWidget(left_col_widget)
        ctrl_cols_lo.addWidget(left_col_scroll, stretch=1)

        # right column: position check + motor status + scalers
        right_col_scroll = QScrollArea()
        right_col_scroll.setWidgetResizable(True)
        right_col_scroll.setFrameShape(QFrame.Shape.NoFrame)
        right_col_widget = QWidget()
        right_col_lo = QVBoxLayout(right_col_widget)
        right_col_lo.setSpacing(4)
        right_col_lo.setContentsMargins(0, 0, 0, 0)
        self._buildPositionCheck(right_col_lo)
        self._buildMotorStatus(right_col_lo)
        self._buildScalerControl(right_col_lo)
        right_col_lo.addStretch()
        right_col_scroll.setWidget(right_col_widget)
        ctrl_cols_lo.addWidget(right_col_scroll, stretch=1)

        right_splitter.addWidget(ctrl_columns)
        right_splitter.addWidget(self._buildLogGroup())

        right_splitter.setStretchFactor(0, 3)  # controls
        right_splitter.setStretchFactor(1, 2)  # log
        return right_splitter

    # -- Scan Control --------------------------------------------------------

    def _buildScanControl(self, parent):
        sc = QGroupBox("Scan Control")
        lo = QVBoxLayout(sc)

        r = QHBoxLayout(); r.addWidget(QLabel("Path:"))
        self._profile_combo = QComboBox()
        self._profile_combo.addItems([PROFILE_NONE, PROFILE_AUTOGEN] + sorted(self._profiles.keys()))
        self._profile_combo.setCurrentText(PROFILE_NONE)
        self._profile_combo.activated.connect(self._onPathProfileChanged)
        r.addWidget(self._profile_combo, stretch=1); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel(f"LG layers (0-{MAX_LG_LAYERS}):"))
        self._lg_spin = QSpinBox(); self._lg_spin.setRange(0, MAX_LG_LAYERS)
        self._lg_spin.valueChanged.connect(self._onLgLayersChanged)
        r.addWidget(self._lg_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Start:"))
        names = [m.name for m in self.engine.path]
        self._start_combo = QComboBox(); self._start_combo.setEditable(True)
        self._start_combo.addItems(names); self._start_combo.setMinimumWidth(80)
        self._start_combo.activated.connect(self._onStartSelected)
        r.addWidget(self._start_combo)
        r.addWidget(QLabel("Count:"))
        self._count_spin = QSpinBox(); self._count_spin.setRange(0, len(names))
        self._count_spin.setSpecialValueText("all")
        self._count_spin.valueChanged.connect(lambda _: self._drawPathPreview())
        r.addWidget(self._count_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Dwell (s):"))
        self._dwell_spin = QDoubleSpinBox(); self._dwell_spin.setRange(1, 9999)
        self._dwell_spin.setValue(DEFAULT_DWELL); self._dwell_spin.setDecimals(0)
        r.addWidget(self._dwell_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Pos. threshold (mm):"))
        self._thresh_spin = QDoubleSpinBox(); self._thresh_spin.setRange(0.01, 10.0)
        self._thresh_spin.setValue(DEFAULT_POS_THRESHOLD); self._thresh_spin.setSingleStep(0.1)
        self._thresh_spin.setDecimals(2)
        r.addWidget(self._thresh_spin); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("Beam threshold (nA):"))
        self._beam_thresh_spin = QDoubleSpinBox(); self._beam_thresh_spin.setRange(0.0, 1000.0)
        self._beam_thresh_spin.setValue(DEFAULT_BEAM_THRESHOLD)
        self._beam_thresh_spin.setSingleStep(0.1); self._beam_thresh_spin.setDecimals(2)
        self._beam_thresh_spin.setSpecialValueText("off")
        r.addWidget(self._beam_thresh_spin); lo.addLayout(r)

        bf = QHBoxLayout()
        self._btn_start = QPushButton("Start Scan"); self._btn_start.setProperty("cssClass", "green")
        self._btn_start.clicked.connect(self._cmdStart); bf.addWidget(self._btn_start)
        self._btn_pause = QPushButton("Pause"); self._btn_pause.setProperty("cssClass", "warn")
        self._btn_pause.clicked.connect(self._cmdPause); bf.addWidget(self._btn_pause)
        self._btn_stop = QPushButton("Stop"); self._btn_stop.setProperty("cssClass", "danger")
        self._btn_stop.clicked.connect(self._cmdStop); bf.addWidget(self._btn_stop)
        lo.addLayout(bf)

        bf2 = QHBoxLayout()
        self._btn_skip = QPushButton("Skip Module"); self._btn_skip.clicked.connect(self._cmdSkip)
        bf2.addWidget(self._btn_skip)
        self._btn_ack = QPushButton("Ack Error"); self._btn_ack.setProperty("cssClass", "warn")
        self._btn_ack.clicked.connect(self._cmdAckError); bf2.addWidget(self._btn_ack)
        lo.addLayout(bf2)

        r = QHBoxLayout()
        self._lbl_progress = QLabel("Progress: --/--"); r.addWidget(self._lbl_progress)
        self._progress_bar = QProgressBar(); self._progress_bar.setMaximumWidth(140)
        r.addWidget(self._progress_bar); lo.addLayout(r)

        self._lbl_current = QLabel("Current:  --"); lo.addWidget(self._lbl_current)
        self._lbl_eta = QLabel("ETA:      --")
        self._lbl_eta.setStyleSheet(f"color: {C.DIM};"); lo.addWidget(self._lbl_eta)
        self._lbl_dwell_cd = QLabel("")
        self._lbl_dwell_cd.setStyleSheet(f"color: {C.GREEN};"); lo.addWidget(self._lbl_dwell_cd)

        parent.addWidget(sc)

    def _buildDirectControl(self, parent):
        dc = QGroupBox("Direct Control"); lo = QVBoxLayout(dc)
        self._btn_move = QPushButton("Move to Starting Point")
        self._btn_move.clicked.connect(self._cmdMoveToModule); lo.addWidget(self._btn_move)
        self._btn_reset = QPushButton("Reset to Beam Center")
        self._btn_reset.setProperty("cssClass", "accent")
        self._btn_reset.clicked.connect(self._cmdResetCenter); lo.addWidget(self._btn_reset)
        parent.addWidget(dc)

    def _buildMotorStatus(self, parent):
        ms = QGroupBox("Motor Status"); lo = QVBoxLayout(ms)
        self._motor_state_labels = {}
        for title, axis, fields in [
            ("X Motor", "x", [
                ("Encoder", "x_encoder"), ("RBV", "x_rbv"), ("VAL", "x_val"),
                ("MOVN", "x_movn"), ("SPMG", "x_spmg"), ("VELO", "x_velo"),
                ("ACCL", "x_accl"), ("TDIR", "x_tdir"), ("MSTA", "x_msta"), ("ATHM", "x_athm")]),
            ("Y Motor", "y", [
                ("Encoder", "y_encoder"), ("RBV", "y_rbv"), ("VAL", "y_val"),
                ("MOVN", "y_movn"), ("SPMG", "y_spmg"), ("VELO", "y_velo"),
                ("ACCL", "y_accl"), ("TDIR", "y_tdir"), ("MSTA", "y_msta"), ("ATHM", "y_athm")])]:
            # title row with inline state badge
            tr = QHBoxLayout()
            tl = QLabel(title)
            tl.setStyleSheet(f"color: {C.ACCENT}; font: bold 13pt 'Consolas';")
            tr.addWidget(tl)
            sl = QLabel("Idle")
            sl.setStyleSheet(f"color: {C.DIM}; font: bold 12pt 'Consolas'; "
                             f"background: #21262d; border: 1px solid #30363d; "
                             f"border-radius: 3px; padding: 1px 6px;")
            tr.addWidget(sl)
            tr.addStretch()
            self._motor_state_labels[axis] = sl
            lo.addLayout(tr)
            g = QGridLayout(); g.setSpacing(2)
            half = (len(fields) + 1) // 2
            for i, (label, key) in enumerate(fields):
                c = 0 if i < half else 2; r = i % half
                ln = QLabel(f"{label}:")
                ln.setStyleSheet(f"color: {C.DIM}; font: 12pt 'Consolas';")
                ln.setAlignment(Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter)
                g.addWidget(ln, r, c)
                lv = QLabel("--"); lv.setMinimumWidth(90)
                g.addWidget(lv, r, c + 1)
                self._status_labels[key] = lv
            lo.addLayout(g)
        parent.addWidget(ms)

    def _buildScalerControl(self, parent):
        sc = QGroupBox("Scalers"); lo = QVBoxLayout(sc)

        r = QHBoxLayout()
        btn_refresh = QPushButton("Refresh Now")
        btn_refresh.clicked.connect(self._pollScalers)
        r.addWidget(btn_refresh)
        r.addStretch()
        lo.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Poll:"))
        self._scaler_interval_slider = QSlider(Qt.Orientation.Horizontal)
        self._scaler_interval_slider.setRange(20, 100)  # 2.0s - 10.0s in 0.1s steps
        self._scaler_interval_slider.setValue(50)        # default 5.0s
        self._scaler_interval_slider.setTickPosition(QSlider.TickPosition.TicksBelow)
        self._scaler_interval_slider.setTickInterval(10)
        self._scaler_interval_slider.valueChanged.connect(self._onScalerIntervalChanged)
        r.addWidget(self._scaler_interval_slider)
        self._lbl_scaler_interval = QLabel("5.0 s")
        self._lbl_scaler_interval.setMinimumWidth(40)
        r.addWidget(self._lbl_scaler_interval)
        lo.addLayout(r)

        parent.addWidget(sc)

    def _onScalerIntervalChanged(self, val):
        sec = val / 10.0
        self._lbl_scaler_interval.setText(f"{sec:.1f} s")
        self._scaler_timer.setInterval(int(sec * 1000))

    def _disableControls(self):
        for w in (self._btn_start, self._btn_pause, self._btn_stop,
                  self._btn_skip, self._btn_ack, self._btn_move, self._btn_reset):
            w.setEnabled(False)
        for w in (self._start_combo, self._profile_combo, self._lg_spin,
                  self._count_spin, self._dwell_spin, self._thresh_spin,
                  self._beam_thresh_spin):
            w.setEnabled(False)

    # -- canvas helpers ------------------------------------------------------

    def _updateCanvas(self):
        if self.observer:
            colors = {}
            for m in self.all_modules:
                if m.mod_type == "LMS": continue
                if m.mod_type == "PbGlass": colors[m.name] = C.MOD_GLASS
                elif m.mod_type == "PbWO4": colors[m.name] = C.MOD_PWO4_BG
            for m in self.scan_modules:
                colors[m.name] = C.MOD_TODO
            self._refreshMap(colors); return

        eng = self.engine
        running = eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED, ScanState.ERROR)
        idle = eng.state in (ScanState.IDLE, ScanState.COMPLETED)

        # Build scan-state colours.  When scaler overlay is on, these are
        # drawn as borders (heat map fill stays visible).  When off, they
        # are drawn as fills.
        colors = {}
        for m in self.all_modules:
            if m.name in self._scan_names or m.mod_type == "LMS": continue
            colors[m.name] = C.MOD_EXCLUDED if running else (
                C.MOD_GLASS if m.mod_type == "PbGlass" else C.MOD_PWO4_BG if m.mod_type == "PbWO4" else C.MOD_LMS)

        count = self._count_spin.value()
        si = self._selected_start_idx
        ei = min(si + count, len(eng.path)) if count > 0 else len(eng.path)
        if not idle:
            si = eng.current_idx; ei = getattr(eng, '_end_idx', len(eng.path))
        for i, mod in enumerate(eng.path):
            if i == eng.current_idx and eng.state == ScanState.DWELLING:
                colors[mod.name] = C.MOD_DWELL
            elif i == eng.current_idx and eng.state in (ScanState.MOVING, ScanState.PAUSED):
                colors[mod.name] = C.MOD_CURRENT
            elif i in eng.error_modules: colors[mod.name] = C.MOD_ERROR
            elif i in eng.completed: colors[mod.name] = C.MOD_DONE
            elif idle and i == self._selected_start_idx: colors[mod.name] = C.MOD_SELECTED
            elif i < si or i >= ei: colors[mod.name] = C.MOD_SKIPPED
            else: colors[mod.name] = C.MOD_TODO
        self._refreshMap(colors, None if idle else eng.path[eng.current_idx + 1:ei])

    def _onModuleClicked(self, name):
        idle = self.engine.state in (ScanState.IDLE, ScanState.COMPLETED)
        if idle and not self.observer:
            mod = self._mod_by_name.get(name)
            if mod:
                if self._mod_dlg is None:
                    self._mod_dlg = ModuleInfoDialog(parent=self)
                sv = self._map._values.get(mod.name)
                self._mod_dlg.setModule(mod, sv)
                self._mod_dlg.show()
                self._mod_dlg.raise_()

    # -- commands ------------------------------------------------------------

    def _cmdStart(self):
        self._onStartSelected(0)
        path = self.engine.path; s = self._selected_start_idx; c = self._count_spin.value()
        e = min(s + c, len(path)) if c > 0 else len(path)
        oob = [path[i].name for i in range(s, e) if not ptrans_in_limits(*module_to_ptrans(path[i].x, path[i].y))]
        if oob:
            ns = ", ".join(oob[:5]) + (f" ... ({len(oob)} total)" if len(oob) > 5 else "")
            self._log(f"BLOCKED: {len(oob)} modules outside limits: {ns}", level="error")
            QMessageBox.critical(self, "Out of Bounds", f"{len(oob)} outside limits:\n{ns}"); return
        self.engine.dwell_time = self._dwell_spin.value()
        self.engine.pos_threshold = self._thresh_spin.value()
        self.engine.beam_threshold = self._beam_thresh_spin.value()
        self.engine.start(self._selected_start_idx, count=c)

    def _cmdPause(self):
        eng = self.engine
        if eng.state == ScanState.PAUSED:
            eng.resume_scan(); self._btn_pause.setText("Pause")
        elif eng.state in (ScanState.MOVING, ScanState.DWELLING):
            eng.pause_scan(); self._btn_pause.setText("Resume")

    def _cmdStop(self):
        if self.engine.state != ScanState.IDLE:
            self.engine.stop_scan(); self._btn_pause.setText("Pause")
        else:
            epics_stop(self.ep); self._log("Motors stopped")

    def _cmdSkip(self):      self.engine.skip_module()
    def _cmdAckError(self):  self.engine.acknowledge_error()

    def _cmdMoveToModule(self):
        self._onStartSelected(0)
        if not self.scan_modules: return
        self._moveToModule(self.scan_modules[self._selected_start_idx])

    def _cmdResetCenter(self):
        self._log(f"Resetting to beam centre ptrans({BEAM_CENTER_X}, {BEAM_CENTER_Y})")
        if epics_move_to(self.ep, BEAM_CENTER_X, BEAM_CENTER_Y):
            self._setTarget(BEAM_CENTER_X, BEAM_CENTER_Y, "Beam Center")

    # -- polling (5 Hz) ------------------------------------------------------

    def _poll(self):
        # detect scan engine moving to a new module
        eng = self.engine
        if eng.state in (ScanState.MOVING, ScanState.DWELLING) and eng.current_idx != self._last_scan_idx:
            self._last_scan_idx = eng.current_idx
            mod = eng.current_module
            if mod:
                px, py = module_to_ptrans(mod.x, mod.y)
                self._setTarget(px, py, mod.name)
        elif eng.state == ScanState.IDLE:
            self._last_scan_idx = -1
        self._updateStatus()
        self._updateCanvas()
        self._updateScanInfo()
        self._updateButtons()
        self._updateBeamDisplay()
        self._checkEncoder()

    def _updateStatus(self):
        for key, lbl in self._status_labels.items():
            val = self.ep.get(key, "--")
            if val == "--" or val is None:
                lbl.setText("--"); lbl.setStyleSheet(f"color: {C.DIM};"); continue
            if key.endswith("_msta"):        txt = f"0x{int(val):X}"
            elif key.endswith("_spmg"):      txt = f"{SPMG_LABELS.get(int(val), '?')}({int(val)})"
            elif key.endswith(("_movn", "_athm", "_tdir")): txt = str(int(val))
            elif isinstance(val, float):     txt = f"{val:.3f}"
            else:                            txt = str(val)
            fg = C.TEXT
            if key.endswith("_movn") and int(val) == 1: fg = C.YELLOW
            elif key.endswith("_spmg"):
                sv = int(val)
                fg = C.RED if sv == SPMG.STOP else C.ORANGE if sv == SPMG.PAUSE else C.GREEN
            lbl.setText(txt); lbl.setStyleSheet(f"color: {fg};")

        # update motor state badges
        for axis, sl in self._motor_state_labels.items():
            spmg = self.ep.get(f"{axis}_spmg", None)
            movn = self.ep.get(f"{axis}_movn", None)
            if spmg is not None and int(spmg) == SPMG.STOP:
                sl.setText("Stop"); fg, bg = C.RED, "#3d1214"
            elif spmg is not None and int(spmg) == SPMG.PAUSE:
                sl.setText("Pause"); fg, bg = C.ORANGE, "#3d2a0e"
            elif movn is not None and int(movn) == 1:
                sl.setText("Moving"); fg, bg = C.YELLOW, "#3d3010"
            else:
                sl.setText("Idle"); fg, bg = C.DIM, "#21262d"
            sl.setStyleSheet(f"color: {fg}; font: bold 12pt 'Consolas'; "
                             f"background: {bg}; border: 1px solid #30363d; "
                             f"border-radius: 3px; padding: 1px 6px;")

        running = self.engine.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED, ScanState.ERROR)
        update_position_check(self._pos_labels, self.ep,
                              self._target_px, self._target_py, self._target_name,
                              scanning=running, pos_threshold=self.engine.pos_threshold)

    def _updateScanInfo(self):
        eng = self.engine
        sc = {ScanState.IDLE: C.DIM, ScanState.MOVING: C.YELLOW, ScanState.DWELLING: C.GREEN,
              ScanState.PAUSED: C.ORANGE, ScanState.ERROR: C.RED, ScanState.COMPLETED: C.ACCENT}
        self._lbl_state.setText(eng.state)
        self._lbl_state.setStyleSheet(f"color: {sc.get(eng.state, C.DIM)}; font: bold 15pt 'Consolas'; background: transparent;")
        done = len(eng.completed)
        s = getattr(eng, '_start_idx', 0); e = getattr(eng, '_end_idx', len(eng.path))
        total = e - s
        self._lbl_progress.setText(f"Progress: {done}/{total}")
        self._progress_bar.setMaximum(max(total, 1)); self._progress_bar.setValue(done)
        mod = eng.current_module
        self._lbl_current.setText(f"Current:  {mod.name}" if mod else "Current:  --")
        if eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED):
            eta = eng.eta_seconds; h, rem = divmod(int(eta), 3600); m, s = divmod(rem, 60)
            self._lbl_eta.setText(f"ETA:      {h}h {m:02d}m {s:02d}s")
        elif eng.state == ScanState.IDLE and eng.path:
            vx = self.ep.get("x_velo", DEFAULT_VELO_X) or DEFAULT_VELO_X
            vy = self.ep.get("y_velo", DEFAULT_VELO_Y) or DEFAULT_VELO_Y
            eta = estimate_scan_time(eng.path, self._selected_start_idx,
                                     self._count_spin.value(), self._dwell_spin.value(), vx, vy)
            if eta > 0:
                h, rem = divmod(int(eta), 3600); m, s = divmod(rem, 60)
                self._lbl_eta.setText(f"ETA:      ~{h}h {m:02d}m {s:02d}s")
            else: self._lbl_eta.setText("ETA:      --")
        else: self._lbl_eta.setText("ETA:      --")
        if eng.beam_tripped:
            self._lbl_dwell_cd.setText("Dwell:    BEAM TRIP -- waiting")
            self._lbl_dwell_cd.setStyleSheet(f"color: {C.RED}; font: bold 13pt 'Consolas';")
        elif eng.state == ScanState.DWELLING:
            self._lbl_dwell_cd.setText(f"Dwell:    {eng.dwell_remaining:.1f}s remaining")
            self._lbl_dwell_cd.setStyleSheet(f"color: {C.GREEN};")
        else:
            self._lbl_dwell_cd.setText(""); self._lbl_dwell_cd.setStyleSheet(f"color: {C.GREEN};")

    def _updateButtons(self):
        if self.observer: return
        eng = self.engine
        running = eng.state in (ScanState.MOVING, ScanState.DWELLING, ScanState.PAUSED, ScanState.ERROR)
        has_path = len(eng.path) > 0
        self._btn_start.setEnabled(not running and has_path)
        self._btn_pause.setEnabled(running)
        self._btn_stop.setEnabled(True)
        self._btn_skip.setEnabled(eng.state == ScanState.DWELLING)
        self._btn_ack.setEnabled(eng.state == ScanState.ERROR)
        self._start_combo.setEnabled(not running and has_path)
        self._count_spin.setEnabled(not running and has_path)
        self._profile_combo.setEnabled(not running)
        self._lg_spin.setEnabled(not running and self._active_profile == PROFILE_AUTOGEN)
        # dwell/thresholds are only applied at scan start, so lock them while running
        self._dwell_spin.setEnabled(not running)
        self._thresh_spin.setEnabled(not running)
        self._beam_thresh_spin.setEnabled(not running)
        self._btn_move.setEnabled(not running and has_path)
        self._btn_reset.setEnabled(not running)


def main():
    run_scan_gui(SnakeScanWindow)


if __name__ == "__main__":
    main()
