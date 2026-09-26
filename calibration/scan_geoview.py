"""
HyCal geo-view widget with integrated scaler overlay.

``HyCalScanMapWidget`` extends ``scripts/hycal_geoview.HyCalMapWidget`` with
the scan state layer (module colours, path preview, limit box, motor
crosshair), drawn over the live scaler rates.
"""
from __future__ import annotations

import os
import sys
from typing import Dict

from PyQt6.QtCore import Qt, QRectF, QPointF
from PyQt6.QtGui import QColor, QPen

_SCRIPTS = os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "scripts")
if _SCRIPTS not in sys.path:
    sys.path.append(_SCRIPTS)

from hycal_geoview import HyCalMapWidget, PALETTES, PALETTE_NAMES, cmap_qcolor  # noqa: E402

from scan_utils import (  # noqa: E402
    C, BEAM_CENTER_X, BEAM_CENTER_Y,
    PTRANS_X_MIN, PTRANS_X_MAX, PTRANS_Y_MIN, PTRANS_Y_MAX,
)

__all__ = ["HyCalScanMapWidget", "PALETTES", "PALETTE_NAMES"]


class HyCalScanMapWidget(HyCalMapWidget):
    BG_COLOR = QColor("#0a0e14")
    NO_DATA_COLOR = QColor("#15181d")
    _EXCLUDED_COLOR = QColor(C.MOD_EXCLUDED)
    LOG_FLOOR = 1e-6

    def __init__(self, all_modules, parent=None):
        super().__init__(parent, shrink=0.90, margin=8, margin_top=8,
                         margin_bottom=8, show_colorbar=False,
                         enable_zoom_pan=True, enable_inline_range_edit=False)
        # the GUIs place their own Reset button
        self._reset_btn.deleteLater()
        self._reset_btn = None
        self._colors = {}           # scan state colours
        self._path_line = []
        self._dash_line = []
        self._marker_hx = self._marker_hy = None
        self._highlight = None
        self._lim_hx_min = PTRANS_X_MIN - BEAM_CENTER_X
        self._lim_hx_max = PTRANS_X_MAX - BEAM_CENTER_X
        self._lim_hy_min = BEAM_CENTER_Y - PTRANS_Y_MAX
        self._lim_hy_max = BEAM_CENTER_Y - PTRANS_Y_MIN
        self._scaler_enabled = True
        self._scaler_auto = True
        self._scaler_auto_exclude = 5  # exclude top N hot channels from auto range
        self._vmax = 1000.0
        self.set_palette("viridis")
        self.set_modules(all_modules)

    # -- scan state (setters don't repaint; the GUI calls update()) ------

    def setModuleColors(self, c):  self._colors = c
    def setPathPreview(self, p):   self._path_line = p
    def setDashPreview(self, p):   self._dash_line = p
    def setHighlight(self, n):     self._highlight = n

    def setMarkerPosition(self, hx, hy):
        self._marker_hx = hx; self._marker_hy = hy

    def modCenter(self, m):
        if self._layout_dirty:
            self._recompute_layout()
        return self.geo_to_canvas(m.x, m.y)

    # -- scaler overlay --------------------------------------------------

    def setScalerValues(self, vals: Dict[str, float]):
        self._values = vals
        if self._scaler_auto and vals:
            self._computeAutoRange()
        self.update()

    def setScalerEnabled(self, on: bool):
        self._scaler_enabled = on
        self.update()

    def setScalerRange(self, vmin, vmax): self.set_range(vmin, vmax)
    def setScalerLogScale(self, on):      self.set_log_scale(on)
    def cyclePalette(self):               self.cycle_palette()
    def scalerRange(self):                return self._vmin, self._vmax
    def resetView(self):                  self.reset_view()

    def setScalerAutoRange(self, on: bool):
        self._scaler_auto = on
        if on and self._values:
            self._computeAutoRange()
        self.update()

    def _computeAutoRange(self):
        """Auto-range excluding the top N hottest channels."""
        v = sorted(self._values.values())
        if not v: return
        trimmed = v[:max(len(v) - self._scaler_auto_exclude, 1)]
        self._vmin = trimmed[0]
        self._vmax = trimmed[-1]
        if self._vmin == self._vmax:
            self._vmax = self._vmin + 1.0

    # -- painting --------------------------------------------------------

    def _paint_modules(self, p):
        stops = self.palette_stops()
        show_scaler = self._scaler_enabled and bool(self._values)
        colors = self._colors
        scaler_vals = self._values
        value_to_t = self.value_to_t
        no_data = self.NO_DATA_COLOR
        excluded = self._EXCLUDED_COLOR
        BORDER_W = 2.0
        # cache QColor objects for state colours to avoid re-creating per module
        _qcolor_cache: Dict[str, QColor] = {}
        for name, r in self._rects.items():
            sc_hex = colors.get(name)
            if sc_hex:
                qc = _qcolor_cache.get(sc_hex)
                if qc is None:
                    qc = _qcolor_cache[sc_hex] = QColor(sc_hex)
            if show_scaler:
                sv = scaler_vals.get(name)
                p.fillRect(r, cmap_qcolor(value_to_t(sv), stops) if sv is not None else no_data)
                if sc_hex:
                    p.setPen(QPen(qc, BORDER_W))
                    p.setBrush(Qt.BrushStyle.NoBrush)
                    p.drawRect(r)
            else:
                p.fillRect(r, qc if sc_hex else excluded)

    def _paint_overlays(self, p, w, h):
        # limit box
        p.setPen(QPen(QColor(C.RED), 1, Qt.PenStyle.DashLine))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(QRectF(self.geo_to_canvas(self._lim_hx_min, self._lim_hy_max),
                          self.geo_to_canvas(self._lim_hx_max, self._lim_hy_min)))

        # path lines (white, visible over heat map)
        for pts, style, lw in [(self._path_line, Qt.PenStyle.SolidLine, 1.8),
                                (self._dash_line, Qt.PenStyle.DashLine, 1.0)]:
            if len(pts) >= 2:
                p.setPen(QPen(QColor(C.PATH_LINE), lw, style))
                for i in range(len(pts) - 1):
                    p.drawLine(pts[i], pts[i + 1])

        if self._highlight and self._highlight in self._rects:
            p.setPen(QPen(QColor(C.ACCENT), 2))
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.drawRect(self._rects[self._highlight])

        # beam position marker (red dot + crosshair)
        if self._marker_hx is not None:
            c = self.geo_to_canvas(self._marker_hx, self._marker_hy)
            R, ARM = 5, 10
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(QColor(C.RED))
            p.drawEllipse(c, R, R)
            p.setPen(QPen(QColor(C.RED), 1.5))
            p.drawLine(QPointF(c.x() - ARM, c.y()), QPointF(c.x() + ARM, c.y()))
            p.drawLine(QPointF(c.x(), c.y() - ARM), QPointF(c.x(), c.y() + ARM))

    # -- mouse -----------------------------------------------------------

    def _handle_click(self, pos):
        name = self._hit(pos)
        if name:
            self.moduleClicked.emit(name)

    def _tooltip_text(self, name):
        v = self._values.get(name)
        return f"{name}: {v:.1f}" if v is not None else name
