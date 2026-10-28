#!/usr/bin/env python3
"""
HyCal FADC Scaler Map (PyQt6)
=============================
Polls EPICS scaler channels (B_DET_HYCAL_FADC_<name>) for every HyCal
module and displays a live colour-coded geo map.

Usage
-----
    python scripts/hycal_scaler_map.py              # real EPICS
    python scripts/hycal_scaler_map.py --sim         # simulation (random)

Options: --database FILE (default database/hycal_map.json), --theme THEME.
"""

from __future__ import annotations

import argparse
import random
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout,
    QPushButton, QLabel,
)
from PyQt6.QtCore import QTimer, Qt
from PyQt6.QtGui import QFont, QColor, QPen

from hycal_geoview import (
    Module, load_modules, HyCalMapWidget,
    apply_theme_palette, set_theme, available_themes, THEME,
    ColorRangeControl, make_info_label,
)


# ---- Paths & constants -----------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = SCRIPT_DIR / ".." / "database"
MODULES_JSON = DB_DIR / "hycal_map.json"

SCALER_PV = "B_DET_HYCAL_FADC_{label}:c"
# Module types that have a scaler PV; Veto (V1-V4) has none.
SCALER_TYPES = ("PbWO4", "PbGlass", "LMS")
POLL_INTERVAL_MS = 2_500


# ---- EPICS interfaces ------------------------------------------------------

class RealScalerEPICS:
    """Read scaler PVs via pyepics."""

    def __init__(self, modules: List[Module]):
        import epics as _epics
        self._pvs: Dict[str, object] = {}
        for m in modules:
            pv = _epics.PV(SCALER_PV.format(label=m.name), connection_timeout=2.0)
            self._pvs[m.name] = pv

    def get(self, name: str) -> Optional[float]:
        pv = self._pvs.get(name)
        if pv and pv.connected:
            return pv.get()
        return None

    def connection_count(self) -> Tuple[int, int]:
        n = sum(1 for pv in self._pvs.values() if pv.connected)
        return n, len(self._pvs)


class SimulatedScalerEPICS:
    """Return random values for testing."""

    def __init__(self, modules: List[Module]):
        self._rng = random.Random(0)
        self._names = [m.name for m in modules]

    def get(self, name: str) -> Optional[float]:
        return self._rng.uniform(0, 1000)

    def connection_count(self) -> Tuple[int, int]:
        return len(self._names), len(self._names)


# ---- HyCal map widget (CoG crosshair overlay) ------------------------------

class ScalerMapWidget(HyCalMapWidget):
    """Simple value → colour map with palette cycle and log scale.

    Adds a center-of-gravity crosshair overlay: call ``set_cog(x, y)``
    with the rate-weighted centroid (mm, in the same coords as the
    Module x/y); pass ``None`` to hide.
    """

    def __init__(self, parent=None):
        super().__init__(parent, min_size=(500, 500), include_lms=True)
        self._vmax = 1000.0
        self._cog: Optional[Tuple[float, float]] = None

    def _fmt_value(self, v: float) -> str:
        return f"{v:.0f}"

    def set_cog(self, xy: Optional[Tuple[float, float]]):
        self._cog = xy
        self.update()

    def _paint_overlays(self, p, w, h):
        # Hover highlight first so the crosshair sits on top of it.
        super()._paint_overlays(p, w, h)
        if self._cog is None:
            return
        cx, cy = self.geo_to_canvas(self._cog[0], self._cog[1]).x(), \
                 self.geo_to_canvas(self._cog[0], self._cog[1]).y()
        # Crosshair: black outline + ACCENT inner so it's readable on any
        # palette.
        L = 18
        for color, width in ((QColor("#000000"), 4.0),
                             (QColor(THEME.ACCENT), 2.0)):
            p.setPen(QPen(color, width, Qt.PenStyle.SolidLine,
                          Qt.PenCapStyle.RoundCap))
            p.drawLine(int(cx - L), int(cy), int(cx + L), int(cy))
            p.drawLine(int(cx), int(cy - L), int(cx), int(cy + L))
        p.setPen(Qt.PenStyle.NoPen)
        p.setBrush(QColor(THEME.ACCENT))
        p.drawEllipse(int(cx - 3), int(cy - 3), 6, 6)


# ---- Main window -----------------------------------------------------------

class ScalerMapWindow(QMainWindow):

    def __init__(self, modules: List[Module], epics_source, simulation: bool):
        super().__init__()
        self._modules = modules
        self._ep = epics_source
        self._simulation = simulation
        self._values: Dict[str, float] = {}
        self._polling = True

        self._build_ui()

        self._timer = QTimer(self)
        self._timer.timeout.connect(self._refresh)
        self._timer.start(POLL_INTERVAL_MS)
        self._refresh()

    def _build_ui(self):
        self.setWindowTitle("HyCal Scaler Map" +
                            ("  [SIMULATION]" if self._simulation
                             else "  [REALTIME]"))
        self.resize(800, 860)
        apply_theme_palette(self)

        central = QWidget()
        self.setCentralWidget(central)
        root = QVBoxLayout(central)
        root.setContentsMargins(8, 8, 8, 8)
        root.setSpacing(6)

        # -- top bar --
        top = QHBoxLayout()
        lbl = QLabel("HYCAL SCALER MAP")
        lbl.setFont(QFont("Monospace", 14, QFont.Weight.Bold))
        lbl.setStyleSheet(f"color:{THEME.ACCENT};")
        top.addWidget(lbl)

        mode = "SIMULATION" if self._simulation else "REALTIME"
        mode_clr = THEME.WARN if self._simulation else THEME.SUCCESS
        mode_lbl = QLabel(mode)
        mode_lbl.setFont(QFont("Monospace", 10, QFont.Weight.Bold))
        mode_lbl.setStyleSheet(f"color:{mode_clr};")
        top.addWidget(mode_lbl)
        top.addStretch()

        self._poll_btn = self._make_btn("Polling: ON", THEME.SUCCESS,
                                        self._toggle_polling)
        top.addWidget(self._poll_btn)
        top.addWidget(self._make_btn("Refresh Now", THEME.TEXT,
                                     self._refresh))
        root.addLayout(top)

        # -- map --
        self._map = ScalerMapWidget()
        self._map.set_modules(self._modules)
        self._map.moduleHovered.connect(self._on_hover)
        self._map.paletteClicked.connect(self._map.cycle_palette)
        root.addWidget(self._map, stretch=1)

        # -- range controls --
        # Starts pinned so the colormap tracks live EPICS data until the
        # user opts out (single-click Auto, or edit a field).
        ctrl = QHBoxLayout()
        self._range_ctrl = ColorRangeControl(
            self._map,
            auto_fit="minmax",
            include_log=True,
            start_pinned=True,
        )
        ctrl.addWidget(self._range_ctrl)
        ctrl.addStretch()

        self._conn_lbl = QLabel("EPICS: --")
        self._conn_lbl.setFont(QFont("Monospace", 10))
        self._conn_lbl.setStyleSheet(f"color:{THEME.TEXT_DIM};")
        ctrl.addWidget(self._conn_lbl)
        root.addLayout(ctrl)

        # -- info bar --
        self._info = make_info_label("Hover over a module")
        root.addWidget(self._info)

    def _make_btn(self, text: str, fg: str, slot) -> QPushButton:
        btn = QPushButton(text)
        btn.setStyleSheet(
            f"QPushButton{{background:{THEME.BUTTON};color:{fg};"
            f"border:1px solid {THEME.BORDER};padding:5px 14px;"
            f"font:bold 11px Monospace;border-radius:8px;}}"
            f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}")
        btn.clicked.connect(slot)
        return btn

    def _refresh(self):
        # Collect rates and compute the rate-weighted center of gravity
        # over all calorimeter crystals (PbWO4 + PbGlass).  LMS pulses
        # don't represent beam hits so they're excluded from the CoG.
        w_total_hz = 0.0           # PbWO4-only sum (for status line)
        cog_w  = 0.0               # weight sum (Σ rate)
        cog_xw = 0.0               # Σ rate · x
        cog_yw = 0.0               # Σ rate · y
        for m in self._modules:
            v = self._ep.get(m.name)
            if v is None:
                continue
            self._values[m.name] = float(v)
            if m.mod_type in ("PbWO4", "PbGlass") and v > 0.0:
                cog_w  += v
                cog_xw += v * m.x
                cog_yw += v * m.y
                if m.mod_type == "PbWO4":
                    w_total_hz += v

        self._map.set_values(self._values)

        if cog_w > 0.0:
            cog_x = cog_xw / cog_w
            cog_y = cog_yw / cog_w
            self._map.set_cog((cog_x, cog_y))
        else:
            cog_x = cog_y = float("nan")
            self._map.set_cog(None)

        if self._values:
            # Pin handles re-fit when on; otherwise no-op.
            self._range_ctrl.notify_values_changed(self._values)

        n_ok, n_total = self._ep.connection_count()
        fg = THEME.SUCCESS if n_ok == n_total else (
             THEME.WARN if n_ok > 0 else THEME.TEXT_DIM)
        self._conn_lbl.setText(f"EPICS: {n_ok}/{n_total}")
        self._conn_lbl.setStyleSheet(f"color:{fg};font:10px Monospace;")

        if self._values:
            lo = min(self._values.values()) / 1000.0
            hi = max(self._values.values()) / 1000.0
            w_total_khz = w_total_hz / 1000.0
            ave_khz = w_total_khz / 1152.0
            cog_str = (f"[{cog_x:+.1f}, {cog_y:+.1f}] mm"
                       if cog_w > 0 else "[—, —]")
            self.statusBar().showMessage(
                f"Data: {lo:.0f}kHz .. {hi:.0f}kHz  "
                f"Channels: {len(self._values)}  "
                f"PbWO4 Total: {w_total_khz:.2f}kHz  "
                f"Ave: {ave_khz:.3f}kHz  "
                f"CoG: {cog_str}")

    def _toggle_polling(self):
        self._polling = not self._polling
        if self._polling:
            self._timer.start(POLL_INTERVAL_MS)
            self._poll_btn.setText("Polling: ON")
            self._poll_btn.setStyleSheet(
                self._poll_btn.styleSheet().replace(THEME.DANGER, THEME.SUCCESS))
        else:
            self._timer.stop()
            self._poll_btn.setText("Polling: OFF")
            self._poll_btn.setStyleSheet(
                self._poll_btn.styleSheet().replace(THEME.SUCCESS, THEME.DANGER))

    def _on_hover(self, name: str):
        parts = [name]
        for m in self._modules:
            if m.name == name:
                parts.append(f"({m.mod_type})")
                break
        v = self._values.get(name)
        if v is not None:
            parts.append(f"{v:.1f}")
        self._info.setText("    ".join(parts))


# ---- Main ------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser(description="HyCal FADC Scaler Map")
    ap.add_argument("--sim", action="store_true",
                    help="Simulation mode (random values, no EPICS)")
    ap.add_argument("--database", type=Path, default=MODULES_JSON,
                    help="Path to hycal_map.json")
    ap.add_argument("--theme", choices=available_themes(), default="dark",
                    help="Colour theme (default: dark)")
    args = ap.parse_args()

    set_theme(args.theme)

    modules = load_modules(args.database)
    print(f"Loaded {len(modules)} modules")
    # Poll and draw only modules with a scaler PV; a Veto tile would only
    # ever show no data.
    modules = [m for m in modules if m.mod_type in SCALER_TYPES]

    if args.sim:
        ep = SimulatedScalerEPICS(modules)
    else:
        try:
            ep = RealScalerEPICS(modules)
        except ImportError:
            print("ERROR: pyepics not available. Use --sim or install pyepics.")
            sys.exit(1)

    app = QApplication(sys.argv)
    win = ScalerMapWindow(modules, ep, simulation=args.sim)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
