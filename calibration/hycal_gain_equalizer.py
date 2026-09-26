#!/usr/bin/env python3
"""
HyCal Gain Equalizer (PyQt6)
=============================
Automatic gain equalization for HyCal crystal modules.  Moves the beam
to each module, collects peak height histograms from prad2_server, finds
the right edge of the Bremsstrahlung spectrum, and adjusts HV via
prad2hvd until the edge converges to a target ADC value.

Shares scan_utils, scan_epics, scan_engine, scan_geoview, and
scan_gui_common (the ScanWindowBase window shell) with hycal_snake_scan.

Usage
-----
    python hycal_gain_equalizer.py                     # simulation
    python hycal_gain_equalizer.py --expert             # expert operator
    python hycal_gain_equalizer.py --observer            # read-only monitor

--database (hycal_map.json) and --paths (scan path profiles JSON) override
the default input files.
"""

from __future__ import annotations

import math
from typing import Dict, List, Optional

from PyQt6.QtWidgets import (
    QWidget, QVBoxLayout, QHBoxLayout, QGroupBox, QPushButton, QLabel,
    QProgressBar, QMessageBox, QSplitter, QSizePolicy, QFrame, QLineEdit,
    QScrollArea,
)
from PyQt6.QtCore import Qt, QRectF
from PyQt6.QtGui import QColor, QFont, QPainter, QPen

from scan_utils import C, module_to_ptrans, SCANNABLE_TYPES
from scan_epics import epics_stop
from scan_engine import (
    DEFAULT_POS_THRESHOLD, DEFAULT_BEAM_THRESHOLD, MAX_LG_LAYERS,
)
from gain_scanner import (
    GainScanEngine, GainScanState, ServerClient, HVClient, draw_histogram,
)
from scan_gui_common import (
    ScanWindowBase, run_scan_gui, PROFILE_AUTOGEN, PROFILE_NONE,
    update_position_check,
    NoScrollSpinBox, NoScrollDoubleSpinBox, NoScrollComboBox,
)


class HistogramWidget(QWidget):
    """Lightweight bar chart for peak height histogram display."""

    PAD_L, PAD_R, PAD_T, PAD_B = 50, 12, 28, 24

    def __init__(self, parent=None):
        super().__init__(parent)
        self._bins: List[int] = []
        self._target_bin: Optional[int] = None
        self._edge_bin: Optional[int] = None
        self._bin_min: float = 0.0      # ADC offset of bin 0
        self._bin_step: float = 1.0     # ADC width per bin
        self._title: str = ""
        self._info: str = ""
        self._log_y: bool = True  # log or linear y scale
        self.setMinimumHeight(140)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def setLogY(self, on: bool):
        self._log_y = on; self.update()

    def setBinning(self, bin_min: float, bin_step: float):
        """Set the ADC mapping (bin index → ADC value).

        Used by the x-axis tick labels.  Safe to call once when the
        engine is created — the analyzer's binning is fixed for the
        lifetime of a scan.
        """
        self._bin_min = float(bin_min)
        self._bin_step = float(bin_step)
        self.update()

    def setData(self, bins: List[int], target_bin: Optional[int] = None,
                edge_bin: Optional[int] = None):
        self._bins = bins
        self._target_bin = target_bin
        self._edge_bin = edge_bin
        self.update()

    def setTitle(self, text: str):
        self._title = text; self.update()

    def setInfo(self, text: str):
        self._info = text; self.update()

    def clear(self):
        self._bins = []; self._target_bin = None; self._edge_bin = None
        self._title = ""; self._info = ""
        self.update()

    def paintEvent(self, event):
        p = QPainter(self)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor("#0d1117"))

        L, R, T, B = self.PAD_L, self.PAD_R, self.PAD_T, self.PAD_B
        pw, ph = w - L - R, h - T - B
        # title and info, also when empty
        p.setPen(QColor(C.ACCENT))
        p.setFont(QFont("Consolas", 11, QFont.Weight.Bold))
        p.drawText(QRectF(L, 2, pw, T - 2), Qt.AlignmentFlag.AlignLeft, self._title)
        p.setPen(QColor(C.DIM))
        p.setFont(QFont("Consolas", 10))
        p.drawText(QRectF(L, 2, pw, T - 2), Qt.AlignmentFlag.AlignRight, self._info)
        if pw < 10 or ph < 10 or not self._bins:
            p.end(); return

        bins = self._bins
        n = len(bins)
        # y range labels; drawn before the grid labels, which may overlap them
        p.setPen(QColor(C.DIM))
        p.setFont(QFont("Consolas", 9))
        p.drawText(QRectF(0, T - 2, L - 4, 14),
                   Qt.AlignmentFlag.AlignRight, f"{max(bins) or 1}")
        p.drawText(QRectF(0, T + ph - 7, L - 4, 14),
                   Qt.AlignmentFlag.AlignRight, "1" if self._log_y else "0")
        bar_w = draw_histogram(p, L, T, pw, ph, bins,
                               self._target_bin, self._edge_bin, self._log_y)

        # x-axis tick labels in ADC units (from bin_min / bin_step)
        if self._bin_step > 0:
            adc_min = self._bin_min
            adc_max = self._bin_min + n * self._bin_step
            span = adc_max - adc_min
            # target tick count adapts to plot width so labels don't collide.
            # Budget is ~80 px per label so even when "nice" rounding lands
            # on a smaller step (and we end up with more ticks than the
            # target), adjacent 50-px label boxes still don't overlap.
            target_ticks = max(2, int(pw / 80))
            raw = span / target_ticks
            if raw > 0:
                mag = 10 ** math.floor(math.log10(raw))
                norm = raw / mag
                if   norm < 1.5: nice = 1
                elif norm < 3.0: nice = 2
                elif norm < 7.0: nice = 5
                else:            nice = 10
                tick_step = nice * mag
                # round adc_min UP to a multiple of tick_step (with a tiny
                # epsilon so values that land just barely above an integer
                # multiple aren't pushed to the next tick)
                first_tick = math.ceil(adc_min / tick_step - 1e-9) * tick_step
                p.setFont(QFont("Consolas", 8))
                adc = first_tick
                while adc <= adc_max + 1e-6:
                    bin_idx = (adc - adc_min) / self._bin_step
                    if 0 <= bin_idx <= n:
                        x = L + bin_idx * bar_w
                        p.setPen(QPen(QColor("#30363d"), 1))
                        p.drawLine(int(x), T + ph, int(x), T + ph + 3)
                        p.setPen(QColor(C.DIM))
                        label = f"{int(adc)}" if tick_step >= 1 else f"{adc:g}"
                        # keep edge labels inside the widget bounds:
                        # left ticks → left-align, right ticks → right-align,
                        # interior ticks → centered
                        if x < L + 25:
                            rect = QRectF(x - 2, T + ph + 4, 50, 14)
                            align = Qt.AlignmentFlag.AlignLeft
                        elif x > L + pw - 25:
                            rect = QRectF(x - 48, T + ph + 4, 50, 14)
                            align = Qt.AlignmentFlag.AlignRight
                        else:
                            rect = QRectF(x - 25, T + ph + 4, 50, 14)
                            align = Qt.AlignmentFlag.AlignCenter
                        p.drawText(rect, align, label)
                    adc += tick_step

        p.end()


class GainEqualizerWindow(ScanWindowBase):
    TITLE = "HyCal Gain Equalizer"
    LOG_PREFIX = "gain_eq"
    WINDOW_H = 1000
    LEGEND = [("Converged", C.GREEN), ("Failed", C.RED),
              ("In progress", C.YELLOW), ("Todo", C.MOD_TODO),
              ("Skipped", C.MOD_SKIPPED)]

    def _initState(self):
        self._gain_engine: Optional[GainScanEngine] = None

    def _onPathSet(self, path):
        if not path:
            self._map.setModuleColors({}); self._map.update()

    def _beamTripped(self):
        return bool(self._gain_engine and self._gain_engine.beam_tripped)

    def _stopEngine(self):
        if self._gain_engine:
            self._gain_engine.stop()
            t = getattr(self._gain_engine, '_thread', None)
            if t and t.is_alive():
                t.join(timeout=2.0)

    # -- layout -------------------------------------------------------------

    def _buildRightPane(self):
        # control panel | histogram | event log (vertical splitter)
        right_splitter = QSplitter(Qt.Orientation.Vertical)

        # row 1: control panel (scrollable, vertical only)
        ctrl_scroll = QScrollArea()
        ctrl_scroll.setWidgetResizable(True)
        ctrl_scroll.setFrameShape(QFrame.Shape.NoFrame)
        ctrl_scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        ctrl_w = QWidget(); ctrl_lo = QVBoxLayout(ctrl_w)
        ctrl_lo.setSpacing(4); ctrl_lo.setContentsMargins(0, 0, 0, 0)
        self._buildPathControl(ctrl_lo)
        self._buildGainControl(ctrl_lo)
        self._buildControlPanel(ctrl_lo)
        self._buildPositionCheck(ctrl_lo)
        ctrl_lo.addStretch()
        ctrl_scroll.setWidget(ctrl_w)
        right_splitter.addWidget(ctrl_scroll)

        # row 2: peak height histogram
        hist_group = QGroupBox("Peak Height Histogram")
        hist_lo = QVBoxLayout(hist_group); hist_lo.setContentsMargins(4, 4, 4, 4)
        self._histogram = HistogramWidget()
        hist_lo.addWidget(self._histogram)
        self._hist_group = hist_group
        hist_group.setVisible(False)
        right_splitter.addWidget(hist_group)

        # row 3: event log
        right_splitter.addWidget(self._buildLogGroup())

        right_splitter.setStretchFactor(0, 2)  # controls
        right_splitter.setStretchFactor(1, 3)  # histogram
        right_splitter.setStretchFactor(2, 2)  # log
        return right_splitter

    # -- control panel builders ---------------------------------------------

    def _buildPathControl(self, parent):
        pc = QGroupBox("Scan Path"); lo = QVBoxLayout(pc)
        self._path_group = pc

        # Row 1: Path + LG layers
        r = QHBoxLayout(); r.addWidget(QLabel("Path:"))
        self._profile_combo = NoScrollComboBox()
        self._profile_combo.addItems([PROFILE_NONE, PROFILE_AUTOGEN] + sorted(self._profiles.keys()))
        self._profile_combo.setCurrentText(PROFILE_NONE)
        self._profile_combo.activated.connect(self._onPathProfileChanged)
        r.addWidget(self._profile_combo, stretch=1)
        r.addWidget(QLabel("LG:"))
        self._lg_spin = NoScrollSpinBox(); self._lg_spin.setRange(0, MAX_LG_LAYERS)
        self._lg_spin.setMaximumWidth(60)
        self._lg_spin.valueChanged.connect(self._onLgLayersChanged)
        r.addWidget(self._lg_spin); lo.addLayout(r)

        # Row 2: Start + Count
        r = QHBoxLayout(); r.addWidget(QLabel("Start:"))
        self._start_combo = NoScrollComboBox(); self._start_combo.setEditable(True)
        self._start_combo.setMinimumWidth(70)
        self._start_combo.setMaximumWidth(100)
        self._start_combo.activated.connect(self._onStartSelected)
        r.addWidget(self._start_combo)
        r.addStretch()
        r.addWidget(QLabel("Count:"))
        self._count_spin = NoScrollSpinBox(); self._count_spin.setRange(0, 0)
        self._count_spin.setSpecialValueText("all")
        self._count_spin.setMaximumWidth(110)
        self._count_spin.setMinimumWidth(100)
        self._count_spin.valueChanged.connect(lambda _: self._drawPathPreview())
        r.addWidget(self._count_spin); lo.addLayout(r)

        # Row 3: Thresholds: Pos (mm) + Curr (nA) — Curr right-aligned
        r = QHBoxLayout()
        r.addWidget(QLabel("Thres.  Pos. (mm)"))
        self._thresh_spin = NoScrollDoubleSpinBox(); self._thresh_spin.setRange(0.01, 10.0)
        self._thresh_spin.setValue(DEFAULT_POS_THRESHOLD)
        self._thresh_spin.setSingleStep(0.1); self._thresh_spin.setDecimals(2)
        self._thresh_spin.setMaximumWidth(100)
        r.addWidget(self._thresh_spin)
        r.addStretch()
        r.addWidget(QLabel("Curr. (nA)"))
        self._beam_thresh_spin = NoScrollDoubleSpinBox(); self._beam_thresh_spin.setRange(0.0, 1000.0)
        self._beam_thresh_spin.setValue(DEFAULT_BEAM_THRESHOLD)
        self._beam_thresh_spin.setSingleStep(0.1); self._beam_thresh_spin.setDecimals(2)
        self._beam_thresh_spin.setSpecialValueText("off")
        self._beam_thresh_spin.setMaximumWidth(100)
        r.addWidget(self._beam_thresh_spin)
        lo.addLayout(r)

        parent.addWidget(pc)

    def _buildGainControl(self, parent):
        ge = QGroupBox("Gain Equalization"); lo = QVBoxLayout(ge)
        self._gain_group = ge

        r = QHBoxLayout(); r.addWidget(QLabel("Server:"))
        self._ge_server_edit = QLineEdit("http://clondaq6:5051")
        r.addWidget(self._ge_server_edit); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("HV:"))
        self._ge_hv_edit = QLineEdit("http://clonpc19:8765")
        r.addWidget(self._ge_hv_edit); lo.addLayout(r)

        r = QHBoxLayout(); r.addWidget(QLabel("HV Password:"))
        self._ge_hv_pw = QLineEdit("prad2_admin"); self._ge_hv_pw.setEchoMode(QLineEdit.EchoMode.Password)
        r.addWidget(self._ge_hv_pw); lo.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Target ADC:"))
        self._ge_target = NoScrollSpinBox(); self._ge_target.setRange(500, 4000)
        self._ge_target.setValue(3200); r.addWidget(self._ge_target)
        r.addWidget(QLabel("Min counts:"))
        self._ge_counts = NoScrollSpinBox(); self._ge_counts.setRange(100, 1000000)
        self._ge_counts.setValue(10000); self._ge_counts.setSingleStep(1000)
        r.addWidget(self._ge_counts); lo.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Max iter:"))
        self._ge_maxiter = NoScrollSpinBox(); self._ge_maxiter.setRange(1, 50)
        self._ge_maxiter.setValue(8); r.addWidget(self._ge_maxiter)
        r.addWidget(QLabel("Tolerance:"))
        self._ge_tol = NoScrollSpinBox(); self._ge_tol.setRange(10, 500)
        self._ge_tol.setValue(50); r.addWidget(self._ge_tol); lo.addLayout(r)

        r = QHBoxLayout()
        r.addWidget(QLabel("Edge frac %:"))
        self._ge_edge_frac = NoScrollDoubleSpinBox(); self._ge_edge_frac.setRange(0.1, 20.0)
        self._ge_edge_frac.setValue(5.0); self._ge_edge_frac.setSingleStep(0.5)
        self._ge_edge_frac.setDecimals(1); r.addWidget(self._ge_edge_frac)
        self._ge_log_y = QPushButton("LogY: ON")
        self._ge_log_y.setCheckable(True); self._ge_log_y.setChecked(True)
        self._ge_log_y.clicked.connect(self._toggleLogY)
        self._updateLogYBtnStyle(True)
        r.addWidget(self._ge_log_y); lo.addLayout(r)

        parent.addWidget(ge)

    def _buildControlPanel(self, parent):
        cp = QGroupBox("Control"); lo = QVBoxLayout(cp)

        bf = QHBoxLayout()
        self._btn_start = QPushButton("Start")
        self._btn_start.setProperty("cssClass", "green")
        self._btn_start.clicked.connect(self._cmdStart); bf.addWidget(self._btn_start)
        self._btn_pause = QPushButton("Pause")
        self._btn_pause.setProperty("cssClass", "warn")
        self._btn_pause.clicked.connect(self._cmdPause); bf.addWidget(self._btn_pause)
        self._btn_stop = QPushButton("Stop")
        self._btn_stop.setProperty("cssClass", "danger")
        self._btn_stop.clicked.connect(self._cmdStop); bf.addWidget(self._btn_stop)
        lo.addLayout(bf)

        bf2 = QHBoxLayout()
        self._btn_redo = QPushButton("Redo Current")
        self._btn_redo.setProperty("cssClass", "accent")
        self._btn_redo.clicked.connect(self._cmdRedo); bf2.addWidget(self._btn_redo)
        self._btn_skip = QPushButton("Skip Current")
        self._btn_skip.clicked.connect(self._cmdSkip); bf2.addWidget(self._btn_skip)
        lo.addLayout(bf2)

        self._lbl_progress = QLabel("Progress: --/--"); lo.addWidget(self._lbl_progress)
        self._progress_bar = QProgressBar(); lo.addWidget(self._progress_bar)
        self._lbl_ge_status = QLabel("Idle")
        self._lbl_ge_status.setStyleSheet(f"color:{C.DIM};"); lo.addWidget(self._lbl_ge_status)
        self._lbl_ge_detail = QLabel("")
        self._lbl_ge_detail.setStyleSheet(f"color:{C.DIM};"); lo.addWidget(self._lbl_ge_detail)

        parent.addWidget(cp)

    def _updateLogYBtnStyle(self, on: bool):
        if on:
            self._ge_log_y.setStyleSheet(
                "QPushButton{background:#1f6feb;color:white;"
                "border:1px solid #388bfd;border-radius:3px;padding:5px 12px;"
                "font:bold 13pt 'Consolas';}"
                "QPushButton:hover{background:#388bfd;}")
        else:
            self._ge_log_y.setStyleSheet("")

    def _disableControls(self):
        for w in (self._btn_start, self._btn_pause, self._btn_stop,
                  self._btn_redo, self._btn_skip):
            w.setEnabled(False)

    # -- commands -----------------------------------------------------------

    def _toggleLogY(self):
        on = self._ge_log_y.isChecked()
        self._ge_log_y.setText("LogY: ON" if on else "LogY: OFF")
        self._updateLogYBtnStyle(on)
        self._histogram.setLogY(on)

    def _cmdStart(self):
        if self._gain_engine and self._gain_engine.state not in (
                GainScanState.IDLE, GainScanState.COMPLETED, GainScanState.FAILED):
            return

        # resume from failure or stop — reuse existing engine, retry from current module
        eng = self._gain_engine
        if eng and eng._has_run and eng.state in (GainScanState.FAILED, GainScanState.IDLE):
            self._clearHistogramData()
            eng.start()
            return

        if not self.scan_modules:
            # no path loaded — fall back to single-module mode if the user
            # has clicked a scannable module on the canvas
            if not self._selected_mod_name:
                self._log("Select a path, or click a module on the canvas "
                          "to equalize a single module", level="error")
                return
            mod = self._mod_by_name.get(self._selected_mod_name)
            if mod is None or mod.mod_type not in SCANNABLE_TYPES:
                self._log(f"Module {self._selected_mod_name!r} is not scannable "
                          f"(needs to be PbWO4 or PbGlass)", level="error")
                return
            self._setPath([mod])
            self._log(f"Single-module start: {mod.name}")
        # sync start index from combo selection
        self._onStartSelected(0)

        ro = self.simulation
        server_url = self._ge_server_edit.text().strip()
        hv_url = self._ge_hv_edit.text().strip()
        hv_pw = self._ge_hv_pw.text()

        if not ro and not hv_pw.strip():
            self._log("Expert mode requires HV password", level="error")
            QMessageBox.warning(self, "HV Password Required",
                                "Enter the prad2hvd password before starting in expert mode.")
            return

        # pre-flight: verify server is reachable and build key map
        try:
            test_server = ServerClient(server_url, log_fn=self._log, read_only=ro)
            key_map = test_server.build_key_map()
            mode = "read-only" if ro else "read-write"
            self._log(f"Server OK ({mode}), {len(key_map)} DAQ channels")
        except Exception as e:
            self._log(f"Server error: {e}", level="error"); return

        # pre-flight: verify HV is reachable and the password is accepted
        if not ro:
            try:
                test_hv = HVClient(hv_url, log_fn=self._log, read_only=False)
                test_hv.connect(password=hv_pw)
                test_hv.close()
                self._log("HV pre-flight OK")
            except Exception as e:
                self._log(f"HV error: {e}", level="error")
                QMessageBox.critical(self, "HV Connection Failed", str(e))
                return

        eng = GainScanEngine(
            motor_ep=self.ep,
            server_url=server_url, hv_url=hv_url,
            hv_password=hv_pw, read_only=ro,
            modules=self.scan_modules, log_fn=self._log, key_map=key_map,
            report_prefix="SIM_" if self.simulation else "")
        eng.target_adc = self._ge_target.value()
        eng.min_counts = self._ge_counts.value()
        eng.max_iterations = self._ge_maxiter.value()
        eng.convergence_tol = self._ge_tol.value()
        eng.beam_threshold = self._beam_thresh_spin.value()
        eng.pos_threshold = self._thresh_spin.value()
        eng.analyzer.edge_fraction = self._ge_edge_frac.value() / 100.0
        eng.analyzer.use_log_cumul = self._ge_log_y.isChecked()
        eng.use_log_y = self._ge_log_y.isChecked()
        self._gain_engine = eng
        self._histogram.setBinning(eng.analyzer.bin_min, eng.analyzer.bin_step)
        eng.start(self._selected_start_idx, count=self._count_spin.value())

    def _cmdPause(self):
        eng = self._gain_engine
        if not eng: return
        if eng._paused:
            eng.resume(); self._btn_pause.setText("Pause")
        else:
            eng.pause(); self._btn_pause.setText("Resume")

    def _cmdStop(self):
        eng = self._gain_engine
        if not eng:
            epics_stop(self.ep); self._log("Motors stopped")
            return
        # if currently running → stop the scan (will be in resumable state)
        running = eng.state not in (
            GainScanState.IDLE, GainScanState.COMPLETED, GainScanState.FAILED)
        if running:
            eng.stop()
            self._btn_pause.setText("Pause")
            return
        # already stopped (FAILED / IDLE-with-has_run / COMPLETED) →
        # discard engine, reset to clean state for a fresh start
        self._gain_engine = None
        self._log("Reset to fresh start")

    def _cmdRedo(self):
        if self._gain_engine:
            self._gain_engine.redo_module()
            self._clearHistogramData()

    def _cmdSkip(self):
        if self._gain_engine:
            self._gain_engine.skip_module()
            self._clearHistogramData()

    def _clearHistogramData(self):
        """Drop any visible histogram bars (target just changed).

        Keeps the title/info banners; only the bars are wiped so the
        operator doesn't see prior-module data while the next collection
        spins up.  Resets the live-fetch timer so the next poll fetches
        immediately once new counts start arriving.
        """
        self._histogram.setData([], None, None)
        self._last_hist_fetch = 0

    # -- canvas -------------------------------------------------------------

    def _updateCanvas(self):
        eng = self._gain_engine
        colors: Dict[str, str] = {}

        # non-scan modules: show type colour only when scalers are off
        if not self._map._scaler_enabled:
            for m in self.all_modules:
                if m.name in self._scan_names or m.mod_type == "LMS": continue
                colors[m.name] = (C.MOD_GLASS if m.mod_type == "PbGlass"
                                  else C.MOD_PWO4_BG if m.mod_type == "PbWO4" else C.MOD_LMS)

        # scan path modules
        path = self.scan_modules
        if eng and eng.state not in (GainScanState.IDLE, GainScanState.COMPLETED):
            for i, mod in enumerate(path):
                if i == eng.current_idx and eng.state in (GainScanState.MOVING,):
                    colors[mod.name] = C.YELLOW
                elif i == eng.current_idx:
                    colors[mod.name] = C.ACCENT
                elif i in eng.converged:
                    colors[mod.name] = C.GREEN
                elif i in eng.failed:
                    colors[mod.name] = C.RED
                else:
                    colors[mod.name] = C.MOD_TODO
        else:
            si = self._selected_start_idx
            count = self._count_spin.value()
            ei = min(si + count, len(path)) if count > 0 else len(path)
            if eng:
                for i, mod in enumerate(path):
                    if i in eng.converged: colors[mod.name] = C.GREEN
                    elif i in eng.failed: colors[mod.name] = C.RED
                    elif i < si or i >= ei: colors[mod.name] = C.MOD_SKIPPED
                    elif i == si: colors[mod.name] = C.MOD_SELECTED
                    else: colors[mod.name] = C.MOD_TODO
            else:
                for i, mod in enumerate(path):
                    if i < si or i >= ei: colors[mod.name] = C.MOD_SKIPPED
                    elif i == si: colors[mod.name] = C.MOD_SELECTED
                    else: colors[mod.name] = C.MOD_TODO

        # while running, show a dashed preview of modules still ahead;
        # idle / completed / failed → solid preview of the planned path.
        running = eng and eng.state not in (GainScanState.IDLE,
                                            GainScanState.COMPLETED,
                                            GainScanState.FAILED)
        ahead = (eng.path[eng.current_idx + 1:getattr(eng, '_end_idx', len(eng.path))]
                 if running else None)
        self._refreshMap(colors, ahead)

    # -- polling (5 Hz) -----------------------------------------------------

    def _poll(self):
        self._updateGainStatus()
        self._updatePositionCheck()
        self._updateCanvas()
        self._updateBeamDisplay()
        self._checkEncoder()

    def _updateGainStatus(self):
        eng = self._gain_engine
        if eng is None:
            self._path_group.setVisible(True)
            self._gain_group.setVisible(True)
            self._hist_group.setVisible(False)
            self._btn_start.setText("Start")
            self._btn_start.setEnabled(True)
            self._btn_pause.setEnabled(False)
            self._btn_stop.setText("Stop")
            self._btn_stop.setEnabled(False)
            self._btn_redo.setEnabled(False)
            self._btn_redo.setText("Redo Current")
            self._btn_skip.setEnabled(False)
            self._btn_skip.setText("Skip Current")
            return

        running = eng.state not in (
            GainScanState.IDLE, GainScanState.COMPLETED, GainScanState.FAILED)
        # "resumable": scan was started, then stopped/failed mid-way
        resumable = (not running) and eng._has_run and \
                    eng.state != GainScanState.COMPLETED

        # detect transition to a new module → drop stale histogram bars
        # (covers natural progression and the post-skip transition once
        # the engine processes the skip flag)
        prev_idx = getattr(self, '_hist_last_idx', None)
        if running and eng.current_idx != prev_idx:
            self._clearHistogramData()
        self._hist_last_idx = eng.current_idx if running else None

        # panel visibility
        self._path_group.setVisible(not running and not resumable)
        self._gain_group.setVisible(not running and not resumable)
        has_data = bool(eng.last_bins) or bool(eng.iteration_history)
        self._hist_group.setVisible(running or has_data)

        # Start button: "Resume" while resumable, "Start" otherwise
        self._btn_start.setText("Resume" if resumable else "Start")
        self._btn_start.setEnabled(not running)
        # Stop button: "Reset" while resumable, "Stop" while running
        self._btn_stop.setText("Reset" if resumable else "Stop")
        self._btn_stop.setEnabled(running or resumable)
        self._btn_pause.setEnabled(running)
        self._btn_redo.setEnabled(running)
        self._btn_skip.setEnabled(running)

        # show current module name on the per-module action buttons
        cur = eng.current_module
        cur_name = cur.name if (running and cur is not None) else ""
        self._btn_redo.setText(f"Redo Current ({cur_name})" if cur_name else "Redo Current")
        self._btn_skip.setText(f"Skip Current ({cur_name})" if cur_name else "Skip Current")

        # sync combo to current_idx so operator sees resume point
        if resumable and 0 <= eng.current_idx < len(eng.path):
            name = eng.path[eng.current_idx].name
            if self._start_combo.currentText() != name:
                idx = self._start_combo.findText(name)
                if idx >= 0:
                    self._start_combo.setCurrentIndex(idx)
                self._selected_start_idx = eng.current_idx
                self._updateCanvasLabel()
        for w in (self._ge_server_edit, self._ge_hv_edit, self._ge_hv_pw,
                  self._ge_target, self._ge_counts, self._ge_maxiter, self._ge_tol,
                  self._ge_edge_frac, self._ge_log_y):
            w.setEnabled(not running and not resumable)
        sc = {GainScanState.IDLE: C.DIM, GainScanState.MOVING: C.YELLOW,
              GainScanState.COLLECTING: C.ACCENT, GainScanState.ANALYZING: C.ACCENT,
              GainScanState.ADJUSTING: C.ORANGE, GainScanState.CONVERGED: C.GREEN,
              GainScanState.FAILED: C.RED, GainScanState.COMPLETED: C.GREEN}
        self._lbl_state.setText(eng.state)
        self._lbl_state.setStyleSheet(f"color:{sc.get(eng.state, C.DIM)};font:bold 15pt 'Consolas';background:transparent;")
        self._lbl_ge_status.setText(eng.state)
        self._lbl_ge_status.setStyleSheet(f"color:{sc.get(eng.state, C.DIM)};")

        done = len(eng.converged) + len(eng.failed)
        total = getattr(eng, '_end_idx', len(eng.path)) - getattr(eng, '_start_idx', 0)
        self._lbl_progress.setText(f"Progress: {done}/{total}")
        self._progress_bar.setMaximum(max(total, 1)); self._progress_bar.setValue(done)

        # edge / ΔV / rate fragments shared by the detail line and the histogram info
        collecting = eng.state == GainScanState.COLLECTING
        edge_dv = []
        if eng.last_edge_adc is not None:
            edge_dv.append(f"edge={eng.last_edge_adc:.0f}")
        if eng.last_dv is not None:
            edge_dv.append(f"ΔV={eng.last_dv:+.0f}")
        rate = [f"{eng.collect_rate:.0f} Hz"] if collecting and eng.collect_rate > 0 else []

        mod = eng.current_module
        parts = []
        if mod: parts.append(mod.name)
        if eng.current_iteration > 0:
            parts.append(f"iter {eng.current_iteration}/{eng.max_iterations}")
        parts += edge_dv
        if collecting:
            parts.append(f"counts={eng.module_counts}")
        parts += rate
        parts.append(f"[{len(eng.converged)}ok {len(eng.failed)}fail]")
        self._lbl_ge_detail.setText("  ".join(parts))

        if mod:
            px, py = module_to_ptrans(mod.x, mod.y)
            if self._target_name != mod.name:
                self._setTarget(px, py, mod.name)

        # update histogram display
        mod_name = mod.name if mod else ""
        target_bin = eng.analyzer.adc_to_bin(eng.target_adc)
        title_parts = [mod_name]
        if eng.last_vmon is not None:
            title_parts.append(f"VMon={eng.last_vmon:.1f}")
        if eng.last_vset is not None:
            title_parts.append(f"VSet={eng.last_vset:.1f}")
        self._histogram.setTitle("  ".join(title_parts))
        self._histogram.setInfo("  ".join(edge_dv + rate))

        # fetch live histogram during collection for preview (~every 2s)
        if collecting and mod and eng.server:
            import time as _time
            now = _time.time()
            if now - getattr(self, '_last_hist_fetch', 0) > 2.0:
                self._last_hist_fetch = now
                key = eng.key_map.get(mod.name)
                if key and eng.module_counts > 0:
                    try:
                        hist = eng.server.get_height_histogram(key)
                        live_bins = hist.get("bins", [])
                        if live_bins:
                            self._histogram.setData(live_bins, target_bin, None)
                    except Exception:
                        pass
        elif eng.last_bins:
            self._histogram.setData(eng.last_bins, target_bin, eng.last_edge_bin)

    def _updatePositionCheck(self):
        update_position_check(self._pos_labels, self.ep,
                              self._target_px, self._target_py,
                              self._target_name)


def main():
    run_scan_gui(GainEqualizerWindow)


if __name__ == "__main__":
    main()
