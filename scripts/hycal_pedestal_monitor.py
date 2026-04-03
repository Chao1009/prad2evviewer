#!/usr/bin/env python3
"""
HyCal Pedestal Monitor (PyQt6)
==============================
Measures FADC250 pedestals on all 7 HyCal crates via SSH, displays
colour-coded HyCal maps, and reports channels with irregular sigma.

RMS (sigma) is parsed from the faV3peds stdout, not from saved files.
Saved .cnf files contain only pedestal means.

Usage
-----
    python hycal_pedestal_monitor.py            # view existing data
    python hycal_pedestal_monitor.py --sim       # test with simulated data
"""

from __future__ import annotations

import argparse
import json
import os
import random
import re
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel, QTextEdit, QMessageBox, QSplitter, QSizePolicy,
    QFileDialog,
)
from PyQt6.QtCore import Qt, QRectF, pyqtSignal, QThread
from PyQt6.QtGui import (
    QPainter, QColor, QPen, QBrush, QFont, QLinearGradient, QPalette,
)


# ===========================================================================
#  Paths & constants
# ===========================================================================

SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = SCRIPT_DIR / ".." / "database"
MODULES_JSON = DB_DIR / "hycal_modules.json"
DAQ_MAP_JSON = DB_DIR / "daq_map.json"
PEDESTALS_DIR = SCRIPT_DIR / ".." / "pedestals"
ORIGINAL_PED_DIR = Path("/usr/clas12/release/2.0.0/parms/fadc250/peds")

NUM_CRATES = 7
CRATE_NAMES = [f"adchycal{i}" for i in range(1, NUM_CRATES + 1)]
CHANNELS_PER_SLOT = 16

# Thresholds for flagging irregular channels  (adjust as needed)
THRESH_PED_MIN = 50.0       # acceptable pedestal mean lower bound
THRESH_PED_MAX = 300.0      # acceptable pedestal mean upper bound
THRESH_DEAD_AVG = 1.0       # avg below this AND rms below THRESH_DEAD_RMS -> DEAD
THRESH_DEAD_RMS = 0.1       # rms below this AND avg below THRESH_DEAD_AVG -> DEAD
THRESH_HIGH_RMS = 1.5       # rms above this -> HIGH RMS
THRESH_DRIFT = 3.0          # |current - original| above this -> DRIFT


# ===========================================================================
#  Data structures
# ===========================================================================

@dataclass
class Module:
    name: str
    mod_type: str   # "PbWO4", "PbGlass", "LMS"
    x: float
    y: float
    sx: float
    sy: float


# ===========================================================================
#  Data loading
# ===========================================================================

def load_modules(path: Path) -> List[Module]:
    with open(path) as f:
        data = json.load(f)
    return [Module(e["n"], e["t"], e["x"], e["y"], e["sx"], e["sy"])
            for e in data]


def load_daq_map(path: Path) -> Dict[Tuple[int, int, int], str]:
    """(crate_index, slot, channel) -> module_name."""
    with open(path) as f:
        data = json.load(f)
    return {(d["crate"], d["slot"], d["channel"]): d["name"] for d in data}


# ===========================================================================
#  Pedestal .cnf parser  (means only -- RMS is NOT in saved files)
# ===========================================================================

def parse_pedestal_file(filepath: Path) -> Dict[int, List[float]]:
    """Return  slot_number -> [16 pedestal means]."""
    slots: Dict[int, List[float]] = {}
    cur_slot: Optional[int] = None
    vals: List[float] = []
    reading = False

    def _flush():
        nonlocal reading, vals
        if cur_slot is not None and reading and vals:
            slots[cur_slot] = vals[:CHANNELS_PER_SLOT]
        vals = []
        reading = False

    with open(filepath) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith("FADC250_CRATE"):
                _flush()
            elif line.startswith("FADC250_SLOT"):
                _flush()
                cur_slot = int(line.split()[1])
            elif line.startswith("FADC250_ALLCH_PED"):
                _flush()
                reading = True
                vals = [float(v) for v in
                        line[len("FADC250_ALLCH_PED"):].split()]
                if len(vals) >= CHANNELS_PER_SLOT:
                    _flush()
            elif reading:
                try:
                    vals.extend(float(v) for v in line.split())
                    if len(vals) >= CHANNELS_PER_SLOT:
                        _flush()
                except ValueError:
                    _flush()
    _flush()
    return slots


def read_all_pedestals(
    ped_dir: Path, suffix: str,
    daq_map: Dict[Tuple[int, int, int], str],
) -> Dict[str, float]:
    """Read all 7 crate files.  Returns  module_name -> pedestal_mean."""
    result: Dict[str, float] = {}
    for ci, cname in enumerate(CRATE_NAMES):
        fp = ped_dir / f"{cname}{suffix}"
        if not fp.exists():
            continue
        for slot, peds in parse_pedestal_file(fp).items():
            for ch, val in enumerate(peds):
                mod = daq_map.get((ci, slot, ch))
                if mod is not None:
                    result[mod] = val
    return result


# ===========================================================================
#  Parse faV3peds stdout for per-channel RMS
# ===========================================================================

_PED_RE = re.compile(
    r"faV3MeasureChannelPedestal:\s*slot\s*(\d+),\s*chan\s*(\d+)\s*=>\s*"
    r"avg\s+([\d.]+),\s*rms\s+([\d.]+),\s*min\s+(\d+),\s*max\s+(\d+)"
)


def parse_measurement_stdout(
    text: str, crate_idx: int,
    daq_map: Dict[Tuple[int, int, int], str],
) -> Dict[str, dict]:
    """Parse stdout lines.  Returns  module_name -> {avg, rms, min, max}."""
    result: Dict[str, dict] = {}
    for m in _PED_RE.finditer(text):
        slot, chan = int(m.group(1)), int(m.group(2))
        mod = daq_map.get((crate_idx, slot, chan))
        if mod is not None:
            result[mod] = {
                "avg": float(m.group(3)), "rms": float(m.group(4)),
                "min": int(m.group(5)),   "max": int(m.group(6)),
            }
    return result


# ===========================================================================
#  Irregular channel detection
# ===========================================================================

def find_irregular_channels(
    measured: Dict[str, dict],
    original: Dict[str, float],
    daq_map: Dict[Tuple[int, int, int], str],
) -> List[str]:
    """Return formatted lines describing flagged channels."""
    # Reverse map:  module_name -> (crate_name, slot, ch)
    rev: Dict[str, Tuple[str, int, int]] = {}
    for (ci, slot, ch), name in daq_map.items():
        rev[name] = (CRATE_NAMES[ci], slot, ch)

    issues: List[str] = []

    for mod, d in sorted(measured.items(), key=lambda kv: rev.get(kv[0], ("", 0, 0))):
        cname, slot, ch = rev.get(mod, ("???", 0, 0))
        loc = f"{cname} slot {slot:2d} ch {ch:2d}"
        avg, rms = d["avg"], d["rms"]

        if avg < THRESH_DEAD_AVG and rms < THRESH_DEAD_RMS:
            issues.append(f"  DEAD          {mod:<6s}  {loc}  "
                          f"avg={avg:.2f}  rms={rms:.3f}")
        elif avg < THRESH_PED_MIN or avg > THRESH_PED_MAX:
            issues.append(f"  OUT OF RANGE  {mod:<6s}  {loc}  "
                          f"avg={avg:.2f}  rms={rms:.3f}  "
                          f"(valid: {THRESH_PED_MIN:.0f}-{THRESH_PED_MAX:.0f})")
        elif rms > THRESH_HIGH_RMS:
            issues.append(f"  HIGH RMS      {mod:<6s}  {loc}  "
                          f"avg={avg:.2f}  rms={rms:.3f}")

        # Drift check (skip dead channels)
        if mod in original and not (avg < THRESH_DEAD_AVG and rms < THRESH_DEAD_RMS):
            delta = avg - original[mod]
            if abs(delta) > THRESH_DRIFT:
                issues.append(f"  DRIFT         {mod:<6s}  {loc}  "
                              f"cur={avg:.2f}  orig={original[mod]:.2f}  "
                              f"delta={delta:+.2f}")

    return issues


# ===========================================================================
#  Colour helpers  (no numpy / matplotlib)
# ===========================================================================

def _lerp(a: int, b: int, t: float) -> int:
    return int(a + (b - a) * t)


_VIRIDIS = [
    (0.00, (68,   1,  84)),
    (0.25, (59,  82, 139)),
    (0.50, (33, 145, 140)),
    (0.75, (94, 201,  98)),
    (1.00, (253, 231, 37)),
]

_RDBU = [
    (0.00, ( 33, 102, 172)),
    (0.25, (103, 169, 207)),
    (0.50, (247, 247, 247)),
    (0.75, (239, 138,  98)),
    (1.00, (178,  24,  43)),
]


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


def _percentile(data: List[float], p: float) -> float:
    if not data:
        return 0.0
    s = sorted(data)
    k = (len(s) - 1) * p / 100.0
    lo = int(k)
    hi = min(lo + 1, len(s) - 1)
    return s[lo] + (k - lo) * (s[hi] - s[lo])


# ===========================================================================
#  HyCal map widget
# ===========================================================================

class HyCalMapWidget(QWidget):
    """Draws a colour-coded HyCal module map using QPainter."""

    module_hovered = pyqtSignal(str)

    _SHRINK = 0.92

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumSize(400, 400)

        self._modules: List[Module] = []
        self._values: Dict[str, float] = {}
        self._title = ""
        self._cmap = "viridis"
        self._center_zero = False
        self._vmin = 0.0
        self._vmax = 1.0
        self._hovered: Optional[str] = None
        self._rects: Dict[str, QRectF] = {}
        self._layout_dirty = True

    # -- public API --

    def set_data(self, modules: List[Module], values: Dict[str, float],
                 title: str, cmap: str = "viridis",
                 center_zero: bool = False):
        self._modules = modules
        self._values = values
        self._title = title
        self._cmap = cmap
        self._center_zero = center_zero
        self._layout_dirty = True

        live = [v for v in values.values() if v != 0.0] or list(values.values())
        if live:
            self._vmin = _percentile(live, 2)
            self._vmax = _percentile(live, 98)
        else:
            self._vmin, self._vmax = 0.0, 1.0
        if center_zero:
            mx = max(abs(self._vmin), abs(self._vmax), 1e-9)
            self._vmin, self._vmax = -mx, mx
        self.update()

    # -- layout --

    def _recompute_layout(self):
        self._rects.clear()
        det = [m for m in self._modules if m.mod_type != "LMS"]
        if not det:
            return
        w, h = self.width(), self.height()
        margin, top, bot = 12, 30, 50
        pw, ph = w - 2 * margin, h - top - bot
        x0 = min(m.x - m.sx / 2 for m in det)
        x1 = max(m.x + m.sx / 2 for m in det)
        y0 = min(m.y - m.sy / 2 for m in det)
        y1 = max(m.y + m.sy / 2 for m in det)
        sc = min(pw / (x1 - x0), ph / (y1 - y0))
        dw, dh = (x1 - x0) * sc, (y1 - y0) * sc
        ox = margin + (pw - dw) / 2
        oy = top + (ph - dh) / 2
        shrink = self._SHRINK
        for m in det:
            mw, mh = m.sx * sc * shrink, m.sy * sc * shrink
            cx = ox + (m.x - x0) * sc
            cy = oy + (y1 - m.y) * sc
            self._rects[m.name] = QRectF(cx - mw / 2, cy - mh / 2, mw, mh)
        self._layout_dirty = False

    def resizeEvent(self, event):
        self._layout_dirty = True
        super().resizeEvent(event)

    # -- painting --

    def paintEvent(self, event):
        if self._layout_dirty:
            self._recompute_layout()
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, QColor("#0a0e14"))

        # Title
        p.setPen(QColor("#c9d1d9"))
        p.setFont(QFont("Monospace", 11, QFont.Weight.Bold))
        p.drawText(QRectF(0, 4, w, 24), Qt.AlignmentFlag.AlignCenter,
                   self._title)

        if not self._rects:
            if not self._values:
                p.setPen(QColor("#555555"))
                p.setFont(QFont("Monospace", 12))
                p.drawText(QRectF(0, 0, w, h), Qt.AlignmentFlag.AlignCenter,
                           "No data")
            p.end()
            return

        stops = _VIRIDIS if self._cmap == "viridis" else _RDBU
        vmin, vmax = self._vmin, self._vmax

        # Modules
        for name, rect in self._rects.items():
            v = self._values.get(name)
            if v is not None:
                t = ((v - vmin) / (vmax - vmin)) if vmax > vmin else 0.5
                p.fillRect(rect, _cmap_qcolor(t, stops))
            else:
                p.fillRect(rect, QColor("#1a1a2e"))

        # Hover highlight
        if self._hovered and self._hovered in self._rects:
            p.setPen(QPen(QColor("#58a6ff"), 2.0))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._hovered])

        # Colour bar
        cb_w = min(300, w - 80)
        cb_h, cb_x, cb_y = 14, (w - min(300, w - 80)) / 2, h - 40
        grad = QLinearGradient(cb_x, 0, cb_x + cb_w, 0)
        for t, (r, g, b) in stops:
            grad.setColorAt(t, QColor(r, g, b))
        p.fillRect(QRectF(cb_x, cb_y, cb_w, cb_h), QBrush(grad))
        p.setPen(QPen(QColor("#555"), 0.5))
        p.drawRect(QRectF(cb_x, cb_y, cb_w, cb_h))
        p.setPen(QColor("#8b949e"))
        p.setFont(QFont("Monospace", 9))
        p.drawText(QRectF(cb_x, cb_y + cb_h + 2, 60, 14),
                   Qt.AlignmentFlag.AlignLeft, f"{vmin:.1f}")
        p.drawText(QRectF(cb_x + cb_w - 60, cb_y + cb_h + 2, 60, 14),
                   Qt.AlignmentFlag.AlignRight, f"{vmax:.1f}")
        mid = (vmin + vmax) / 2
        p.drawText(QRectF(cb_x + cb_w / 2 - 30, cb_y + cb_h + 2, 60, 14),
                   Qt.AlignmentFlag.AlignCenter, f"{mid:.1f}")
        p.end()

    # -- mouse tracking --

    def mouseMoveEvent(self, event):
        pos = event.position()
        found = None
        for name, rect in self._rects.items():
            if rect.contains(pos):
                found = name
                break
        if found != self._hovered:
            self._hovered = found
            self.update()
            if found:
                self.module_hovered.emit(found)


# ===========================================================================
#  Measurement thread
# ===========================================================================

class MeasureThread(QThread):
    """Runs faV3peds on each crate via SSH, emitting stdout per crate."""
    progress = pyqtSignal(int, str)        # crate_idx, status_msg
    crate_done = pyqtSignal(int, str)      # crate_idx, combined stdout+stderr
    crate_error = pyqtSignal(int, str)     # crate_idx, error_msg

    def run(self):
        for i, cname in enumerate(CRATE_NAMES):
            self.progress.emit(
                i, f"Measuring {cname} ({i + 1}/{NUM_CRATES})...")
            cmd = (f'ssh {cname} '
                   f'"cd ~/prad2_daq/prad2evviewer/pedestals; '
                   f'faV3peds {cname}_latest.cnf"')
            try:
                r = subprocess.run(cmd, shell=True, capture_output=True,
                                   text=True, timeout=180)
                self.crate_done.emit(i, r.stdout + "\n" + r.stderr)
            except subprocess.TimeoutExpired:
                self.crate_error.emit(i, f"{cname}: TIMEOUT (180 s)")
            except Exception as e:
                self.crate_error.emit(i, f"{cname}: {e}")


# ===========================================================================
#  Main window
# ===========================================================================

class PedestalMonitorWindow(QMainWindow):

    def __init__(self, modules: List[Module],
                 daq_map: Dict[Tuple[int, int, int], str],
                 sim: bool = False):
        super().__init__()
        self._modules = modules
        self._daq_map = daq_map
        self._sim = sim

        self._original: Dict[str, float] = {}
        self._latest: Dict[str, float] = {}
        self._measured: Dict[str, dict] = {}   # from stdout: name->{avg,rms,...}

        self._build_ui()
        self._load_data()

    # ---- UI construction ----

    def _build_ui(self):
        self.setWindowTitle("HyCal Pedestal Monitor")
        self.resize(1400, 820)
        self._apply_dark_palette()

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # -- top bar --
        top = QHBoxLayout()
        lbl = QLabel("HYCAL PEDESTAL MONITOR")
        lbl.setFont(QFont("Monospace", 14, QFont.Weight.Bold))
        lbl.setStyleSheet("color:#58a6ff;")
        top.addWidget(lbl)
        top.addStretch()

        self._measure_btn = self._make_btn(
            "Measure Pedestals", "#3fb950", self._on_measure)
        top.addWidget(self._measure_btn)
        self._reload_btn = self._make_btn(
            "Reload Files", "#c9d1d9", self._load_data)
        top.addWidget(self._reload_btn)
        self._save_btn = self._make_btn(
            "Save Report", "#d29922", self._on_save_report)
        top.addWidget(self._save_btn)
        root.addLayout(top)

        # -- maps + report splitter --
        splitter = QSplitter(Qt.Orientation.Vertical)

        maps = QWidget()
        ml = QHBoxLayout(maps)
        ml.setContentsMargins(0, 0, 0, 0)
        ml.setSpacing(8)
        self._map_left = HyCalMapWidget()
        self._map_right = HyCalMapWidget()
        self._map_left.module_hovered.connect(self._on_hover)
        self._map_right.module_hovered.connect(self._on_hover)
        ml.addWidget(self._map_left)
        ml.addWidget(self._map_right)
        splitter.addWidget(maps)

        self._report = QTextEdit()
        self._report.setReadOnly(True)
        self._report.setFont(QFont("Monospace", 10))
        self._report.setStyleSheet(
            "QTextEdit{background:#161b22;color:#8b949e;"
            "border:1px solid #30363d;border-radius:4px;}")
        self._report.setMaximumHeight(220)
        splitter.addWidget(self._report)
        splitter.setStretchFactor(0, 4)
        splitter.setStretchFactor(1, 1)
        root.addWidget(splitter)

        # -- info bar --
        self._info = QLabel("Hover over a module for details")
        self._info.setFont(QFont("Monospace", 11))
        self._info.setStyleSheet(
            "QLabel{background:#161b22;color:#c9d1d9;padding:4px 8px;"
            "border:1px solid #30363d;border-radius:4px;}")
        self._info.setFixedHeight(28)
        root.addWidget(self._info)

        self.statusBar().setStyleSheet(
            "QStatusBar{color:#8b949e;font:10px Monospace;}")

    def _make_btn(self, text: str, fg: str, slot) -> QPushButton:
        btn = QPushButton(text)
        btn.setStyleSheet(
            f"QPushButton{{background:#21262d;color:{fg};"
            f"border:1px solid #30363d;padding:6px 16px;"
            f"font:bold 12px Monospace;border-radius:4px;}}"
            f"QPushButton:hover{{background:#30363d;}}"
            f"QPushButton:disabled{{color:#555;}}")
        btn.clicked.connect(slot)
        return btn

    def _apply_dark_palette(self):
        pal = self.palette()
        for role, colour in [
            (QPalette.ColorRole.Window,       "#0d1117"),
            (QPalette.ColorRole.WindowText,   "#c9d1d9"),
            (QPalette.ColorRole.Base,         "#161b22"),
            (QPalette.ColorRole.Text,         "#c9d1d9"),
            (QPalette.ColorRole.Button,       "#21262d"),
            (QPalette.ColorRole.ButtonText,   "#c9d1d9"),
            (QPalette.ColorRole.Highlight,    "#58a6ff"),
        ]:
            pal.setColor(role, QColor(colour))
        self.setPalette(pal)

    # ---- data loading ----

    def _load_data(self):
        if self._sim:
            self._load_sim_data()
            return

        # Original pedestals
        if ORIGINAL_PED_DIR.exists():
            self._original = read_all_pedestals(
                ORIGINAL_PED_DIR, "_ped.cnf", self._daq_map)
        else:
            self._original = {}

        # Latest measured pedestals
        if PEDESTALS_DIR.exists():
            self._latest = read_all_pedestals(
                PEDESTALS_DIR, "_latest.cnf", self._daq_map)
        else:
            self._latest = {}

        n_orig = len(self._original)
        n_lat = len(self._latest)
        self.statusBar().showMessage(
            f"Loaded {n_orig} original, {n_lat} latest channels", 5000)
        self._update_maps()
        self._update_report()

    def _load_sim_data(self):
        rng = random.Random(42)
        names = [m.name for m in self._modules if m.mod_type != "LMS"]
        self._original.clear()
        self._latest.clear()
        self._measured.clear()
        for m in self._modules:
            if m.mod_type == "LMS":
                continue
            o = rng.gauss(160, 25)
            l = o + rng.gauss(0, 1.0)
            self._original[m.name] = o
            self._latest[m.name] = l
            self._measured[m.name] = {
                "avg": l, "rms": abs(rng.gauss(0.7, 0.15)),
                "min": int(l) - 3, "max": int(l) + 3,
            }
        # Dead channels (avg ~ 0, rms ~ 0)
        for n in rng.sample(names, 15):
            self._original[n] = 0.0
            self._latest[n] = 0.0
            self._measured[n].update(avg=0.0, rms=0.0)
        # Out-of-range channels
        for n in rng.sample(names, 3):
            val = rng.choice([rng.uniform(10, 40), rng.uniform(320, 500)])
            self._original[n] = val
            self._latest[n] = val + rng.gauss(0, 0.5)
            self._measured[n].update(avg=self._latest[n], rms=abs(rng.gauss(0.7, 0.2)))
        # High-RMS channels
        for n in rng.sample(names, 5):
            if self._measured[n]["avg"] >= THRESH_DEAD_AVG:
                self._measured[n]["rms"] = rng.uniform(1.8, 5.0)
        # Drift channels
        for n in rng.sample(names, 4):
            if self._measured[n]["avg"] >= THRESH_PED_MIN:
                drift = rng.choice([-1, 1]) * rng.uniform(4.0, 12.0)
                self._latest[n] = self._original[n] + drift
                self._measured[n]["avg"] = self._latest[n]
        self._update_maps()
        self._update_report()

    # ---- update views ----

    def _update_maps(self):
        has_latest = bool(self._latest)
        cur = self._latest if has_latest else self._original
        label = "Current" if has_latest else "Original"

        self._map_left.set_data(
            self._modules, cur, f"{label} Pedestal Mean")

        if has_latest and self._original:
            delta = {n: cur[n] - self._original[n]
                     for n in cur if n in self._original}
            self._map_right.set_data(
                self._modules, delta,
                "Mean Difference (Current \u2212 Original)",
                cmap="rdbu", center_zero=True)
        else:
            self._map_right.set_data(
                self._modules, {}, "Mean Difference (no comparison data)")

    def _update_report(self):
        lines: List[str] = []

        def _stats(label: str, peds: Dict[str, float]):
            vals = list(peds.values())
            live = [v for v in vals if v >= THRESH_DEAD_AVG]
            dead = sum(1 for v in vals if v < THRESH_DEAD_AVG)
            if not vals:
                lines.append(f"{label}: no data")
                return
            if live:
                avg = sum(live) / len(live)
                lines.append(
                    f"{label}: {len(vals)} ch, {dead} dead, "
                    f"mean={avg:.1f}  min={min(live):.1f}  max={max(live):.1f}")
            else:
                lines.append(f"{label}: {len(vals)} ch, ALL dead")

        if self._original:
            _stats("Original", self._original)
        if self._latest:
            _stats("Current ", self._latest)

        if self._measured:
            issues = find_irregular_channels(
                self._measured, self._original, self._daq_map)
            lines.append("")
            if issues:
                lines.append(f"IRREGULAR CHANNELS  ({len(issues)} flagged):")
                lines.extend(issues)
            else:
                lines.append("All channels within normal parameters.")

        self._report.setPlainText("\n".join(lines))

    # ---- hover info ----

    def _on_hover(self, name: str):
        parts = [name]
        for m in self._modules:
            if m.name == name:
                parts.append(f"({m.mod_type})")
                break
        if name in self._latest:
            parts.append(f"ped: {self._latest[name]:.2f}")
        elif name in self._original:
            parts.append(f"ped: {self._original[name]:.2f}")
        if name in self._original and name in self._latest:
            delta = self._latest[name] - self._original[name]
            parts.append(f"orig: {self._original[name]:.2f}")
            parts.append(f"delta: {delta:+.2f}")
        if name in self._measured:
            parts.append(f"rms: {self._measured[name]['rms']:.3f}")
        self._info.setText("    ".join(parts))

    # ---- save report ----

    def _on_save_report(self):
        text = self._report.toPlainText().strip()
        if not text:
            QMessageBox.information(self, "Save Report", "Nothing to save.")
            return
        path, _ = QFileDialog.getSaveFileName(
            self, "Save Irregular Channel Report",
            str(PEDESTALS_DIR / "pedestal_report.txt"),
            "Text files (*.txt);;All files (*)")
        if not path:
            return
        with open(path, "w") as f:
            f.write(text + "\n")
        self.statusBar().showMessage(f"Report saved to {path}", 5000)

    # ---- measurement ----

    def _on_measure(self):
        reply = QMessageBox.warning(
            self, "Pedestal Measurement",
            "WARNING: Pedestal measurement will INTERRUPT DAQ running!\n"
            "Only proceed when DAQ is IDLE.\n\n"
            "Proceed with measurement on all 7 crates?",
            QMessageBox.StandardButton.Yes | QMessageBox.StandardButton.No,
            QMessageBox.StandardButton.No)
        if reply != QMessageBox.StandardButton.Yes:
            return

        self._measure_btn.setEnabled(False)
        self._reload_btn.setEnabled(False)
        self._measure_btn.setText("Measuring...")
        self._measured.clear()

        self._thread = MeasureThread()
        self._thread.progress.connect(
            lambda _i, s: self.statusBar().showMessage(s))
        self._thread.crate_done.connect(self._on_crate_done)
        self._thread.crate_error.connect(self._on_crate_error)
        self._thread.finished.connect(self._on_measure_finished)
        self._thread.start()

    def _on_crate_done(self, idx: int, stdout: str):
        parsed = parse_measurement_stdout(stdout, idx, self._daq_map)
        self._measured.update(parsed)
        self._report.append(
            f"  {CRATE_NAMES[idx]}: {len(parsed)} channels measured")

    def _on_crate_error(self, idx: int, msg: str):
        self._report.append(f"  ERROR: {msg}")

    def _on_measure_finished(self):
        self._measure_btn.setEnabled(True)
        self._reload_btn.setEnabled(True)
        self._measure_btn.setText("Measure Pedestals")
        self.statusBar().showMessage("Measurement complete", 5000)

        # Re-read latest files (now on shared mount)
        if PEDESTALS_DIR.exists():
            self._latest = read_all_pedestals(
                PEDESTALS_DIR, "_latest.cnf", self._daq_map)
        self._update_maps()
        self._update_report()


# ===========================================================================
#  Main
# ===========================================================================

def main():
    ap = argparse.ArgumentParser(description="HyCal Pedestal Monitor")
    ap.add_argument("--sim", action="store_true",
                    help="Use simulated data for testing")
    ap.add_argument("--modules-db", type=Path, default=MODULES_JSON)
    ap.add_argument("--daq-map", type=Path, default=DAQ_MAP_JSON)
    args = ap.parse_args()

    modules = load_modules(args.modules_db)
    daq_map = load_daq_map(args.daq_map)

    app = QApplication(sys.argv)
    win = PedestalMonitorWindow(modules, daq_map, sim=args.sim)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
