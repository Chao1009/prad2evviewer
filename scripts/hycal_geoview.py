"""Shared HyCal geo-view widget for PyQt6 scripts.

Provides the shared theme, a common ``Module`` record with the
hycal_map.json / daq_config.json loaders and HyCal grid helpers, atomic
JSON writes, colour palettes and plot helpers, the ``ZoomHistWidget``
histogram base, an extensible ``HyCalMapWidget`` base class with the
``ColorRangeController`` / ``ColorRangeControl`` range helpers, config
tuning-dock helpers and a QThread worker launcher.  The HyCal GUIs in
this directory (and calibration/scan_geoview.py) subclass the widget to
add overlays, custom fills, or different mouse behaviour.

Typical usage:

    class MyMap(HyCalMapWidget):
        def _paint_modules(self, p):
            # optional custom fill; default uses value colormap
            ...

    w = MyMap(enable_zoom_pan=True)
    w.set_modules(load_modules(MODULES_JSON))
    w.set_values({name: value, ...})
    w.set_range(vmin, vmax)
"""
from __future__ import annotations

import json
import math
import os
import sys
import tempfile
from pathlib import Path
from typing import Callable, Dict, Iterable, List, Optional, Tuple, Union

from PyQt6.QtWidgets import (
    QWidget, QPushButton, QSizePolicy, QToolTip,
    QLineEdit, QLabel, QHBoxLayout, QVBoxLayout, QApplication, QMenu,
    QCheckBox, QComboBox, QDockWidget, QDoubleSpinBox, QFormLayout, QSpinBox,
)
from PyQt6.QtCore import (
    Qt, QRectF, QPointF, QSize, QTimer, QObject, QEvent, QThread, pyqtSignal,
)
from PyQt6.QtGui import (
    QPainter, QColor, QPen, QBrush, QFont, QFontMetricsF, QLinearGradient,
    QPalette, QDoubleValidator,
)


# ---- Shared theme ----
#
# The ``THEME`` class body is the "dark" theme; a new field goes into
# ``THEME`` and into every other ``_THEMES`` entry.  Scripts typically wire
# the theme up via a ``--theme`` CLI flag:
#     parser.add_argument("--theme", choices=available_themes(),
#                         default="dark")
#     ...
#     set_theme(args.theme)           # before constructing any window
#     apply_theme_palette(window)


class THEME:
    """Active palette — class attrs are overwritten by :func:`set_theme`.

    Apple-inspired: binary dark / light surfaces with a single blue accent
    reserved for interactive elements.
    """

    # --- surfaces ---
    BG            = "#000000"   # window background (Pure Black)
    BG_SUBTLE     = "#161b22"   # inset plot/panel tile
    CANVAS        = "#000000"   # chart / HyCal map canvas
    PANEL         = "#1d1d1f"   # input surfaces (text edits, tables, combos)
    BUTTON        = "#1d1d1f"   # raised controls (Primary Dark)
    BUTTON_HOVER  = "#28282a"   # button :hover background (Dark Surface 3)
    ALT_BASE      = "#242426"   # alternating table rows (Dark Surface 5)
    TOOLTIP       = "#2a2a2d"   # hover/info tooltip background (Dark Surface 4)

    # --- lines ---
    BORDER        = "#424245"   # subtle border — Apple rarely uses borders
    GRID          = "#1d1d1f"   # chart gridlines (very faint on dark)

    # --- text ---
    TEXT          = "#ffffff"
    TEXT_STRONG   = "#ffffff"
    TEXT_DIM      = "#86868b"   # Apple secondary grey
    TEXT_MUTED    = "#6e6e73"   # tertiary / disabled

    # --- semantic / state ---
    ACCENT        = "#2997ff"   # Bright Blue — links/highlights on dark
    ACCENT_STRONG = "#0071e3"   # Apple Blue — primary CTA
    ACCENT_BORDER = "#0071e3"   # focus ring
    SUCCESS       = "#30d158"   # iOS green (system green dark)
    WARN          = "#ff9f0a"   # iOS orange
    DANGER        = "#ff453a"   # iOS red
    HIGHLIGHT     = "#ff9f0a"   # orange emphasis / drift
    NO_DATA       = "#1d1d1f"   # map-fill when no value

    # --- misc ---
    SELECT_BORDER = "#ffffff"   # selected-module border (white on dark)


_THEMES: Dict[str, Dict[str, str]] = {
    # Snapshot taken at import, before any set_theme() mutates THEME.
    "dark": {k: v for k, v in vars(THEME).items() if k.isupper()},
    "light": {
        "BG":            "#ffffff",
        "BG_SUBTLE":     "#f5f5f7",
        "CANVAS":        "#f5f5f7",
        "PANEL":         "#ffffff",
        "BUTTON":        "#fafafc",
        "BUTTON_HOVER":  "#ededf2",
        "ALT_BASE":      "#f5f5f7",
        "TOOLTIP":       "#ffffff",
        "BORDER":        "#d2d2d7",
        "GRID":          "#e5e5ea",
        "TEXT":          "#1d1d1f",
        "TEXT_STRONG":   "#000000",
        "TEXT_DIM":      "#6e6e73",
        "TEXT_MUTED":    "#86868b",
        "ACCENT":        "#0066cc",
        "ACCENT_STRONG": "#0071e3",
        "ACCENT_BORDER": "#0071e3",
        "SUCCESS":       "#248a3d",
        "WARN":          "#c93400",
        "DANGER":        "#d70015",
        "HIGHLIGHT":     "#bf5700",
        "NO_DATA":       "#e5e5ea",
        "SELECT_BORDER": "#000000",
    },
}


def available_themes() -> List[str]:
    """Theme names acceptable to :func:`set_theme`."""
    return list(_THEMES.keys())


def set_theme(name: str) -> None:
    """Activate one of :data:`_THEMES` by mutating :class:`THEME` in place.

    Call **before** any window is constructed; stylesheets in this project
    are plain f-strings that read ``THEME.*`` once, at widget creation time.
    Switching themes after the UI is built will not re-render existing
    stylesheets.
    """
    try:
        values = _THEMES[name]
    except KeyError as exc:
        raise ValueError(
            f"Unknown theme {name!r}. Available: {available_themes()}"
        ) from exc
    for key, value in values.items():
        setattr(THEME, key, value)


def apply_theme_palette(widget) -> None:
    """Install the active theme's :class:`QPalette` on ``widget``.

    Sets QPalette roles used by top-level windows throughout the
    ``scripts/`` GUIs. Idempotent.
    """
    pal = widget.palette()
    for role, colour in (
        (QPalette.ColorRole.Window,     THEME.BG),
        (QPalette.ColorRole.WindowText, THEME.TEXT),
        (QPalette.ColorRole.Base,       THEME.PANEL),
        (QPalette.ColorRole.Text,       THEME.TEXT),
        (QPalette.ColorRole.Button,     THEME.BUTTON),
        (QPalette.ColorRole.ButtonText, THEME.TEXT),
        (QPalette.ColorRole.Highlight,  THEME.ACCENT),
    ):
        pal.setColor(role, QColor(colour))
    widget.setPalette(pal)


# ---- themed(qss): legacy dark hex codes -> active theme ----

_QSS_MAP: Dict[str, str] = {
    # base surfaces
    "#0a0e14": "CANVAS",
    "#0d1117": "BG",
    "#161b22": "PANEL",
    "#21262d": "BUTTON",
    "#30363d": "BORDER",
    "#131820": "ALT_BASE",
    # text
    "#c9d1d9": "TEXT",
    "#e6edf3": "TEXT_STRONG",
    "#8b949e": "TEXT_DIM",
    "#555555": "TEXT_MUTED",
    "#555":    "TEXT_MUTED",
    # accents / state
    "#58a6ff": "ACCENT",
    "#1f6feb": "ACCENT_STRONG",
    "#388bfd": "ACCENT_BORDER",
    "#3fb950": "SUCCESS",
    "#d29922": "WARN",
    "#f85149": "DANGER",
    "#f97316": "HIGHLIGHT",
    "#ff2222": "DANGER",
    # widget fills
    "#1a1a2e": "NO_DATA",
}


def themed(qss: str) -> str:
    """Rewrite the legacy dark-theme hex codes of ``_QSS_MAP`` in a Qt
    stylesheet to the active :class:`THEME`, so one stylesheet string serves
    every theme without per-call f-strings::

        label.setStyleSheet(themed("QLabel{background:#161b22;color:#c9d1d9;}"))

    Other colour literals pass through unchanged (one-off accents).
    """
    out = qss
    for hex_code, key in _QSS_MAP.items():
        out = out.replace(hex_code, getattr(THEME, key))
    return out


# Translucent button floating over a canvas (mode / stack toggles).  Raw
# legacy-hex QSS: pass it through themed() when the widget is built, not at
# import, so --theme takes effect.
OVERLAY_BUTTON_QSS = (
    "QPushButton{background:rgba(29,29,31,220);color:#c9d1d9;"
    "border:1px solid #30363d;border-radius:4px;}"
    "QPushButton:hover{background:#28282a;color:#e6edf3;}")

# Right-click menu of a plot canvas; raw like OVERLAY_BUTTON_QSS.
MENU_QSS = (
    "QMenu{background:#161b22;color:#c9d1d9;border:1px solid #30363d;}"
    "QMenu::item:selected{background:#1f6feb;}")


def make_info_label(text: str = "", dim: bool = False) -> QLabel:
    """Fixed-height monospace info bar (hover details, stats) on a PANEL tile."""
    lbl = QLabel(text)
    lbl.setFont(QFont("Monospace", 11))
    lbl.setStyleSheet(
        f"QLabel{{background:{THEME.PANEL};"
        f"color:{THEME.TEXT_DIM if dim else THEME.TEXT};padding:4px 8px;"
        f"border:1px solid {THEME.BORDER};border-radius:8px;}}")
    lbl.setFixedHeight(28)
    return lbl


# ---- Modules, map loaders and HyCal grid ----

# FADC250 crates of the HyCal DAQ; ``daq.crate`` in hycal_map.json is an
# index into CRATE_NAMES.
NUM_CRATES = 7
CRATE_NAMES = [f"adchycal{i}" for i in range(1, NUM_CRATES + 1)]
CHANNELS_PER_SLOT = 16


class Module:
    """A HyCal detector module with geometric size and position.

    ``row``/``col`` are the 1-based grid indices from hycal_map.json (0 for
    LMS/Veto and when unknown); ``crate``/``slot``/``channel`` are the DAQ
    address, -1 when the module is not read out.
    """
    __slots__ = ("name", "mod_type", "x", "y", "sx", "sy",
                 "row", "col", "crate", "slot", "channel")

    def __init__(self, name: str, mod_type: str,
                 x: float, y: float, sx: float, sy: float,
                 row: int = 0, col: int = 0,
                 crate: int = -1, slot: int = -1, channel: int = -1):
        self.name = name
        self.mod_type = mod_type
        self.x = x
        self.y = y
        self.sx = sx
        self.sy = sy
        self.row = row
        self.col = col
        self.crate = crate
        self.slot = slot
        self.channel = channel


def _daq_address(rec: dict) -> Tuple[int, int, int]:
    """(crate, slot, channel) of a hycal_map.json record; -1s if unmapped."""
    d = rec.get("daq") or {}
    try:
        return (int(d.get("crate", -1)), int(d.get("slot", -1)),
                int(d.get("channel", -1)))
    except (TypeError, ValueError):
        return (-1, -1, -1)


def load_modules(path: Path) -> List[Module]:
    """Load modules from a hycal_map.json file.

    Per-record schema: {n, t, geo:{sx,sy,x,y,row,col,...}, daq:{...}?, bst:{...}?}.
    Records without a ``geo`` block (none in practice) are skipped.
    """
    with open(path, encoding="utf-8") as f:
        data = json.load(f)
    out: List[Module] = []
    for e in data:
        g = e.get("geo")
        if not g:
            continue
        out.append(Module(e["n"], e["t"], g["x"], g["y"], g["sx"], g["sy"],
                          g.get("row", 0), g.get("col", 0), *_daq_address(e)))
    return out


def load_daq_map(path: Path) -> Dict[Tuple[int, int, int], str]:
    """(crate, slot, channel) -> module name from hycal_map.json.

    Records without a DAQ address (boosters, PRad-1 V1-V4) are skipped.
    """
    with open(path, encoding="utf-8") as f:
        entries = json.load(f)
    out: Dict[Tuple[int, int, int], str] = {}
    for e in entries:
        addr = _daq_address(e)
        if addr[0] >= 0:
            out[addr] = e["n"]
    return out


def load_roc_tag_map(path: Path) -> Dict[int, int]:
    """ROC tag -> crate index for the ``type == "roc"`` entries of
    daq_config.json's ``roc_tags`` (ti_slave, tdc, gem, ... excluded).

    Tags may be hex strings ("0x80") or integers; malformed entries are
    skipped.
    """
    with open(path, encoding="utf-8") as f:
        cfg = json.load(f)
    out: Dict[int, int] = {}
    for r in cfg.get("roc_tags", []):
        if r.get("type") != "roc":
            continue
        try:
            tag = r["tag"]
            tag = int(tag, 16) if isinstance(tag, str) else int(tag)
            out[tag] = int(r["crate"])
        except (KeyError, TypeError, ValueError):
            continue
    return out


def atomic_json_write_many(items: Iterable[Tuple[Path, object]], *,
                           indent: int = 2, **dump_kw) -> None:
    """Write each ``(path, value)`` as indented JSON plus a trailing newline.

    Every file is written to a temp file next to its destination before any
    destination is replaced, and the temp files are removed on failure.
    Extra keywords go to ``json.dump`` (e.g. ``sort_keys=True``).
    """
    temporary: List[Tuple[str, Path]] = []
    try:
        for path, value in items:
            path = Path(path)
            path.parent.mkdir(parents=True, exist_ok=True)
            fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.",
                                             dir=path.parent)
            temporary.append((temp_name, path))
            with os.fdopen(fd, "w", encoding="utf-8") as stream:
                json.dump(value, stream, indent=indent, **dump_kw)
                stream.write("\n")
        for temp_name, path in temporary:
            os.replace(temp_name, path)
    except BaseException:
        for temp_name, _ in temporary:
            try:
                os.unlink(temp_name)
            except OSError:
                pass
        raise


def atomic_json_write(path: Path, value, **kw) -> None:
    """Single-file :func:`atomic_json_write_many`."""
    atomic_json_write_many([(path, value)], **kw)


# LMS / Veto modules sit far outside the calorimeter in hycal_map.json;
# map views that show them re-place them in a row just below HyCal.
AUX_TYPES = frozenset({"LMS", "Veto"})
AUX_ROW_Y = -640.0
AUX_ROW_SIZE = 50.0


def place_aux_row(modules: List[Module],
                  x_by_name: Dict[str, float]) -> List[Module]:
    """Move every module named in ``x_by_name`` to (x, AUX_ROW_Y) with an
    AUX_ROW_SIZE square cell.  Mutates the records in place and returns
    ``modules``."""
    for m in modules:
        x = x_by_name.get(m.name)
        if x is not None:
            m.x, m.y, m.sx, m.sy = x, AUX_ROW_Y, AUX_ROW_SIZE, AUX_ROW_SIZE
    return modules


# PbWO4 (W) grid: 34 x 34 cells, 1-based row/col, with the 2x2 beam hole at
# rows/cols 17-18.
W_GRID = 34
W_HOLE = (17, 18)


def hole_ring(row: int, col: int) -> int:
    """Chebyshev distance of a 1-based W cell from the beam hole: 0 = hole,
    1 = the ring right around it (C++ kInnerBound, under the absorber)."""
    lo, hi = W_HOLE
    return max(lo - row, 0, row - hi, lo - col, col - hi)


def edge_depth(row: int, col: int, n: int = W_GRID) -> int:
    """Ring index of a 1-based cell counted from the outer edge of an
    ``n`` x ``n`` grid: 1 = outermost ring."""
    return min(row, col, n + 1 - row, n + 1 - col)


# ---- Colour palettes ----

# First two palettes (``rainbow`` and ``blue-yellow``) match the
# corresponding palettes in prad2hvmon's web monitor (resources/
# monitor_geo_view.js: ``rainbow``, ``darkblue``) so the desktop and web
# views use the same color language.  Order matters — palette cycle
# starts at index 0.
PALETTES: Dict[str, List[Tuple[float, Tuple[int, int, int]]]] = {
    "rainbow": [
        (0.00, (30,   58,  95)), (0.25, (59,  130, 246)),
        (0.50, (45,  212, 160)), (0.75, (234, 179,   8)),
        (1.00, (245, 101, 101)),
    ],
    "blue-yellow": [
        (0.00, (11,   22,  40)),
        (0.50, (59,  158, 255)),
        (1.00, (234, 179,   8)),
    ],
    "viridis": [
        (0.00, (68,   1,  84)), (0.25, (59,  82, 139)),
        (0.50, (33, 145, 140)), (0.75, (94, 201,  98)),
        (1.00, (253, 231,  37)),
    ],
    "inferno": [
        (0.00, (0,     0,   4)), (0.25, (120,  28, 109)),
        (0.50, (229,  89,  52)), (0.75, (253, 198,  39)),
        (1.00, (252, 255, 164)),
    ],
    "coolwarm": [
        (0.00, (59,   76, 192)), (0.25, (141, 176, 254)),
        (0.50, (221, 221, 221)), (0.75, (245, 148, 114)),
        (1.00, (180,   4,  38)),
    ],
    "hot": [
        (0.00, (11,   0,   0)), (0.33, (230,   0,   0)),
        (0.66, (255, 210,   0)), (1.00, (255, 255, 255)),
    ],
    "blue-orange": [
        (0.00, (10,   42, 110)), (0.25, (30,   90, 180)),
        (0.50, (80,   80,  80)), (0.75, (220, 120,  30)),
        (1.00, (249, 115,  22)),
    ],
    "greyscale": [
        (0.00, (20,   20,  20)), (1.00, (240, 240, 240)),
    ],
}
PALETTE_NAMES: List[str] = list(PALETTES.keys())


def _lerp(a: int, b: int, t: float) -> int:
    return int(a + (b - a) * t)


def cmap_qcolor(t: float, stops) -> QColor:
    """Map ``t`` in [0, 1] to a QColor along the given palette stops."""
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


def cmap_rgb_array(t, stops):
    """Array form of :func:`cmap_qcolor`: ``t`` (numpy array of any shape)
    -> uint8 RGB array of shape ``t.shape + (3,)``, values clamped to the
    end stops."""
    import numpy as np
    xs = np.array([s for s, _ in stops], dtype=np.float64)
    rgb = np.array([c for _, c in stops], dtype=np.float64)
    return np.stack([np.interp(t, xs, rgb[:, k]) for k in range(3)],
                    axis=-1).astype(np.uint8)


# Categorical colours for peaks / clusters / channels, same order as ``PC``
# in resources/viewer.js.
SERIES_COLORS: Tuple[str, ...] = (
    "#00b4d8", "#ff6b6b", "#51cf66", "#ffd43b",
    "#cc5de8", "#ff922b", "#20c997", "#f06595",
)


def series_qcolor(i: int) -> QColor:
    return QColor(SERIES_COLORS[i % len(SERIES_COLORS)])


# ---- Plot helpers ----

def fmt_value(v) -> str:
    """Compact number for labels and tooltips; an em dash for None / NaN."""
    if v is None or (isinstance(v, float) and math.isnan(v)):
        return "—"
    if v == 0:
        return "0"
    return f"{v:.6g}"


def nice_ticks(lo: float, hi: float, max_ticks: int = 6) -> List[float]:
    """Axis ticks at 1/2/2.5/5 x 10^k steps inside [lo, hi] (at most about
    ``max_ticks``); [] for a non-finite range, [lo] for an empty one."""
    if not math.isfinite(lo) or not math.isfinite(hi):
        return []
    if hi <= lo:
        return [lo]
    raw = (hi - lo) / max(max_ticks - 1, 1)
    mag = 10 ** math.floor(math.log10(raw)) if raw > 0 else 1
    step = mag
    for c in (1, 2, 2.5, 5, 10):
        if c * mag >= raw:
            step = c * mag
            break
    ticks: List[float] = []
    v = math.ceil(lo / step) * step
    # The length cap stops a step below the float resolution of lo from
    # looping forever.
    while v <= hi + step * 0.01 and len(ticks) < 1000:
        ticks.append(v)
        v += step
    return ticks


def draw_wave_axes(p: QPainter, r: QRectF, ymin: float, ymax: float,
                   n: int, clk_mhz: float, pad_l: float) -> None:
    """Ticks and labels of an FADC waveform plot drawn in ``r``: five ADC
    ticks on the left (labels start ``pad_l`` px left of ``r``) and the
    time axis in ns for ``n`` samples at ``clk_mhz``."""
    p.setPen(QColor(THEME.TEXT_DIM))
    p.setFont(QFont("Monospace", 8))
    for frac in (0.0, 0.25, 0.5, 0.75, 1.0):
        y = r.bottom() - frac * r.height()
        p.drawLine(int(r.left() - 3), int(y), int(r.left()), int(y))
        val = ymin + frac * (ymax - ymin)
        p.drawText(int(r.left() - pad_l + 2), int(y + 4), f"{val:.0f}")
    tick_every = max(1, n // 8)
    for i in range(0, n, tick_every):
        x = r.left() + i / max(1, n - 1) * r.width()
        p.drawLine(int(x), int(r.bottom()), int(x), int(r.bottom() + 3))
        ns = i * 1e3 / clk_mhz
        p.drawText(int(x - 18), int(r.bottom() + 14), f"{ns:g}")
    p.setFont(QFont("Monospace", 9))
    p.drawText(int(r.left() + r.width() / 2 - 10),
               int(r.bottom() + 26), "ns")


# ---- Zoomable 1-D histogram canvas ----

class ZoomHistWidget(QWidget):
    """Base of the painter-drawn 1-D histograms: left-drag zooms into an x
    range, right-click offers Unzoom.

    Subclasses keep ``_values`` / ``_edges`` (bin contents and edges),
    ``_title`` and the visible range ``_x_lo`` / ``_x_hi`` current and
    paint with the ``_paint_*`` helpers.  ``_log_x`` makes the x axis
    logarithmic while the range is positive.  Hooks: ``_on_view_changed()``
    after a zoom or unzoom, ``_extend_menu(menu)`` for more menu items.
    """

    PAD_L, PAD_R, PAD_T, PAD_B = 55, 16, 24, 36
    N_YTICKS = 5
    BAR_GAP = 1.0        # px left empty right of each bar
    X_LABEL_W = 50       # px, width of an x tick label box
    X_LABEL_FMT = ".0f"  # default x tick labels

    def __init__(self, parent=None):
        super().__init__(parent)
        self._values: List[float] = []
        self._edges: List[float] = []
        self._title = ""
        self._x_lo = 0.0
        self._x_hi = 1.0
        self._log_x = False
        self._drag_start: Optional[float] = None   # data x of the drag
        self._drag_cur: Optional[float] = None
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMouseTracking(True)

    # -- geometry --

    def _plot_rect(self) -> Tuple[int, int, int, int]:
        w, h = self.width(), self.height()
        return (self.PAD_L, self.PAD_T,
                w - self.PAD_L - self.PAD_R, h - self.PAD_T - self.PAD_B)

    def _x_view(self) -> Tuple[float, float, bool]:
        """Visible ``(x_lo, x_hi, log_x)``, with x_hi kept above x_lo."""
        lo = self._x_lo
        hi = self._x_hi if self._x_hi > lo else lo + 1
        return lo, hi, self._log_x and lo > 0

    def _x_map(self) -> Tuple[Callable[[float], float], Callable[[float], float]]:
        """``(to_sx, from_sx)``: data x to widget x and back."""
        px, _py, pw, _ph = self._plot_rect()
        lo, hi, log_x = self._x_view()
        if log_x:
            u_lo, u_hi = math.log10(lo), math.log10(hi)

            def to_sx(v):
                if v <= 0:
                    return px - 1
                return px + (math.log10(v) - u_lo) / (u_hi - u_lo) * pw

            def from_sx(sx):
                if pw <= 0:
                    return lo
                return 10.0 ** (u_lo + (sx - px) / pw * (u_hi - u_lo))
        else:
            def to_sx(v):
                return px + (v - lo) / (hi - lo) * pw

            def from_sx(sx):
                if pw <= 0:
                    return lo
                return lo + (sx - px) / pw * (hi - lo)
        return to_sx, from_sx

    def _visible_values(self, values, edges) -> List[float]:
        """Contents of the bins overlapping the visible x range."""
        lo, hi, _ = self._x_view()
        return [v for i, v in enumerate(values)
                if i + 1 < len(edges) and edges[i + 1] > lo and edges[i] < hi]

    # -- zoom --

    def mousePressEvent(self, event):
        if event.button() == Qt.MouseButton.LeftButton:
            px, py, pw, ph = self._plot_rect()
            mx, my = event.position().x(), event.position().y()
            if px <= mx <= px + pw and py <= my <= py + ph + self.PAD_B:
                self._drag_start = self._x_map()[1](mx)
                self._drag_cur = self._drag_start
                self.update()

    def mouseMoveEvent(self, event):
        if self._drag_start is not None:
            self._drag_cur = self._x_map()[1](event.position().x())
            self.update()

    def mouseReleaseEvent(self, event):
        if (event.button() == Qt.MouseButton.LeftButton
                and self._drag_start is not None):
            to_sx, from_sx = self._x_map()
            sx = event.position().x()
            d_start, d_end = self._drag_start, from_sx(sx)
            self._drag_start = self._drag_cur = None
            # Minimum span in pixels, so that it also holds on a log axis.
            if abs(sx - to_sx(d_start)) > 0.01 * self._plot_rect()[2]:
                self._x_lo = min(d_start, d_end)
                self._x_hi = max(d_start, d_end)
                self._on_view_changed()
            self.update()

    def contextMenuEvent(self, event):
        menu = QMenu(self)
        menu.setStyleSheet(themed(MENU_QSS))
        menu.addAction("Unzoom").triggered.connect(self._unzoom)
        self._extend_menu(menu)
        menu.exec(event.globalPos())

    def _extend_menu(self, menu: QMenu) -> None:
        """Hook: add context-menu items after Unzoom."""

    def _unzoom(self):
        """Show every bin; on a log x axis start at the first positive edge."""
        if self._edges:
            self._x_lo, self._x_hi = self._edges[0], self._edges[-1]
            if self._log_x:
                pos = [e for e in self._edges if e > 0]
                if pos:
                    self._x_lo = pos[0]
        self._on_view_changed()
        self.update()

    def _on_view_changed(self) -> None:
        """Hook: the visible x range changed."""

    # -- painting --

    def _paint_title(self, p: QPainter) -> None:
        if self._title:
            p.setPen(QColor(THEME.ACCENT))
            p.setFont(QFont("Consolas", 10, QFont.Weight.Bold))
            p.drawText(QRectF(self.PAD_L, 2,
                              self.width() - self.PAD_L - self.PAD_R, 20),
                       Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                       self._title)

    def _paint_placeholder(self, p: QPainter, text: str) -> None:
        p.setPen(QColor(THEME.TEXT_MUTED))
        p.setFont(QFont("Consolas", 10))
        p.drawText(QRectF(0, 0, self.width(), self.height()),
                   Qt.AlignmentFlag.AlignCenter, text)

    def _paint_grid(self, p: QPainter) -> None:
        """Dotted lines at the y ticks."""
        px, py, pw, ph = self._plot_rect()
        p.setPen(QPen(QColor(THEME.BUTTON), 1, Qt.PenStyle.DotLine))
        for i in range(self.N_YTICKS + 1):
            sy = py + ph * i / self.N_YTICKS
            p.drawLine(QPointF(px, sy), QPointF(px + pw, sy))

    def _paint_bars(self, p: QPainter, values, edges, color: QColor,
                    to_sy: Callable[[float], float]) -> None:
        """Bars of the visible bins; ``to_sy`` maps a content to widget y."""
        px, py, pw, ph = self._plot_rect()
        lo, hi, _ = self._x_view()
        to_sx = self._x_map()[0]
        p.setPen(Qt.PenStyle.NoPen)
        for i, v in enumerate(values):
            if i + 1 >= len(edges):
                break
            b_lo, b_hi = edges[i], edges[i + 1]
            if b_hi <= lo or b_lo >= hi:
                continue
            sx1 = max(to_sx(b_lo), px)
            sx2 = min(to_sx(b_hi), px + pw)
            bar_top = to_sy(v)
            bar_h = (py + ph) - bar_top
            if bar_h > 0 and sx2 > sx1:
                p.fillRect(QRectF(sx1, bar_top, sx2 - sx1 - self.BAR_GAP, bar_h),
                           color)

    def _paint_drag(self, p: QPainter) -> None:
        """The x band of a drag zoom in progress."""
        if self._drag_start is None or self._drag_cur is None:
            return
        px, py, pw, ph = self._plot_rect()
        to_sx = self._x_map()[0]
        sx1 = max(to_sx(min(self._drag_start, self._drag_cur)), px)
        sx2 = min(to_sx(max(self._drag_start, self._drag_cur)), px + pw)
        if sx2 > sx1:
            p.fillRect(QRectF(sx1, py, sx2 - sx1, ph), QColor(255, 255, 100, 50))
            p.setPen(QPen(QColor(255, 255, 100, 180), 1))
            p.drawRect(QRectF(sx1, py, sx2 - sx1, ph))

    def _lin_y_labels(self, y_hi: float) -> List[str]:
        """Tick labels of a linear 0..y_hi axis, top to bottom."""
        n = self.N_YTICKS
        return [f"{y_hi * (n - i) / n:.0f}" for i in range(n + 1)]

    def _paint_axes(self, p: QPainter, y_labels: List[str],
                    x_ticks: Optional[List[Tuple[float, str]]] = None) -> None:
        """Axis lines, ``y_labels`` top to bottom at the grid lines and
        ``x_ticks`` as (widget x, text); by default nice_ticks() over the
        visible range formatted with X_LABEL_FMT."""
        px, py, pw, ph = self._plot_rect()
        p.setPen(QPen(QColor(THEME.BORDER), 1))
        p.drawLine(QPointF(px, py), QPointF(px, py + ph))
        p.drawLine(QPointF(px, py + ph), QPointF(px + pw, py + ph))
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Consolas", 8))
        for i, text in enumerate(y_labels):
            sy = py + ph * i / self.N_YTICKS
            p.drawText(QRectF(0, sy - 8, self.PAD_L - 4, 16),
                       Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                       text)
        if x_ticks is None:
            lo, hi, _ = self._x_view()
            to_sx = self._x_map()[0]
            x_ticks = [(to_sx(xt), format(xt, self.X_LABEL_FMT))
                       for xt in nice_ticks(lo, hi, max(pw // 60, 2))]
        lw = self.X_LABEL_W
        for sx, text in x_ticks:
            p.drawText(QRectF(sx - lw / 2, py + ph + 2, lw, 16),
                       Qt.AlignmentFlag.AlignCenter, text)


# ---- HyCal map base widget ----

class HyCalMapWidget(QWidget):
    """Extensible HyCal geometry view with value → colour mapping.

    Features
    --------
    * Automatic layout: modules laid out in physical coordinates, axis-correct
      (y flipped so positive y is up).
    * Optional colour bar at the bottom (click to cycle palette).
    * Optional zoom/pan (mouse wheel + drag, middle click to reset, overlay
      Reset button top-right).
    * Optional log-scale value mapping.
    * Hover tooltip and module click signal.
    * Optional in-cell name labels for the module types in ``label_types``.
    * Selected-module outline (``set_selected``); with ``toggle_select=True``
      a module click toggles the selection and a click on empty canvas
      clears it.

    Subclass hooks (override to customise)
    --------------------------------------
    * ``_paint_modules(p)``          — per-module fill loop (default uses
                                        ``set_values`` + current palette).
    * ``_paint_before_modules(p, w, h)`` — drawn after background, before modules.
    * ``_paint_overlays(p, w, h)``   — drawn after modules, before colour bar.
                                        Default paints the name labels, the
                                        selection and the hover highlight.
    * ``_paint_after_colorbar(p, w, h)`` — drawn last (legends etc.).
    * ``_fmt_value(v)``              — vmin/vmax label format.
    * ``_tooltip_text(name)``        — tooltip when hovering a module.
    * ``_label_text(name)``          — in-cell label of a labelled module.
    """

    moduleHovered = pyqtSignal(str)
    moduleClicked = pyqtSignal(str)   # "" means deselect
    paletteClicked = pyqtSignal()
    rangeEdited   = pyqtSignal(float, float)   # user-edited via inline editor

    _CLICK_THRESHOLD = 4

    # Lowest value mapped by value_to_t() in log scale (vmin is clamped to it).
    LOG_FLOOR = 1e-9

    SELECT_PEN_WIDTH = 2.5

    # Colour roles resolve from :class:`THEME` at paint time; a subclass may
    # pin one with a plain class attribute (``BG_COLOR = QColor(...)``).

    @property
    def BG_COLOR(self) -> QColor:
        return QColor(THEME.CANVAS)

    @property
    def NO_DATA_COLOR(self) -> QColor:
        return QColor(THEME.NO_DATA)

    @property
    def HOVER_COLOR(self) -> QColor:
        return QColor(THEME.ACCENT)

    @property
    def SELECT_COLOR(self) -> QColor:
        return QColor(THEME.SELECT_BORDER)

    @property
    def CB_BORDER(self) -> QColor:
        return QColor(THEME.ACCENT)

    @property
    def CB_TEXT(self) -> QColor:
        return QColor(THEME.TEXT_DIM)

    @property
    def EMPTY_TEXT(self) -> QColor:
        return QColor(THEME.TEXT_MUTED)

    def __init__(self, parent=None, *,
                 shrink: float = 0.92,
                 margin: int = 12,
                 margin_top: int = 10,
                 margin_bottom: int = 50,
                 include_lms: bool = False,
                 label_types: Iterable[str] = (),
                 toggle_select: bool = False,
                 show_colorbar: bool = True,
                 enable_zoom_pan: bool = False,
                 enable_inline_range_edit: bool = True,
                 min_size: Tuple[int, int] = (400, 400)):
        super().__init__(parent)
        self._shrink = shrink
        self._margin = margin
        self._margin_top = margin_top
        self._margin_bottom = margin_bottom
        self._include_lms = include_lms
        self._label_types = frozenset(label_types)
        self._toggle_select = toggle_select
        self._show_colorbar = show_colorbar
        self._enable_zoom_pan = enable_zoom_pan
        self._enable_inline_range_edit = enable_inline_range_edit

        self.setMouseTracking(True)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)
        self.setMinimumSize(*min_size)

        self._modules: List[Module] = []
        self._values: Dict[str, float] = {}
        self._vmin = 0.0
        self._vmax = 1.0
        self._log_scale = False
        self._palette_idx = 0
        self._hovered: Optional[str] = None
        self._selected: Optional[str] = None
        self._labelled_names: List[str] = []
        self._rects: Dict[str, QRectF] = {}
        self._rect_names_rev: List[str] = []
        self._geo_bounds: Tuple[float, float, float, float] = (0.0, 1.0, 0.0, 1.0)
        self._cb_rect: Optional[QRectF] = None
        self._layout_dirty = True

        # Inline range editor.  ``_cb_min_hit`` / ``_cb_max_hit`` are
        # screen-space rects (set by each colorbar paint) covering the
        # vmin / vmax label + pencil glyph; ``_inline_editor`` is a child
        # QLineEdit, created on first use, shown over whichever was clicked.
        self._cb_min_hit: Optional[QRectF] = None
        self._cb_max_hit: Optional[QRectF] = None
        self._inline_editor: Optional[QLineEdit] = None
        self._inline_which: Optional[str] = None
        self._inline_cancelled = False

        # zoom / pan state (only used when enable_zoom_pan is True)
        self._zoom = 1.0
        self._pan_x = 0.0
        self._pan_y = 0.0
        self._drag_last: Optional[QPointF] = None
        self._drag_origin: Optional[QPointF] = None
        self._dragging = False

        if enable_zoom_pan:
            self._reset_btn = QPushButton("Reset", self)
            self._reset_btn.setFixedSize(52, 24)
            f = QFont("Consolas", 9)
            f.setBold(True)
            self._reset_btn.setFont(f)
            self._reset_btn.setStyleSheet(
                f"QPushButton{{background:{THEME.BUTTON};color:{THEME.TEXT_DIM};"
                f"border:1px solid {THEME.BORDER};border-radius:8px;}}"
                f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};color:{THEME.TEXT};}}")
            self._reset_btn.clicked.connect(self.reset_view)
        else:
            self._reset_btn = None

    # -- Public API --

    def set_modules(self, modules: List[Module]):
        if self._include_lms:
            self._modules = list(modules)
        else:
            self._modules = [m for m in modules if m.mod_type != "LMS"]
        if self._modules:
            self._geo_bounds = (
                min(m.x - m.sx / 2 for m in self._modules),
                max(m.x + m.sx / 2 for m in self._modules),
                min(m.y - m.sy / 2 for m in self._modules),
                max(m.y + m.sy / 2 for m in self._modules),
            )
        self._labelled_names = [m.name for m in self._modules
                                if m.mod_type in self._label_types]
        self._layout_dirty = True
        self.update()

    def set_values(self, values: Dict[str, float]):
        self._values = values
        self.update()

    def set_selected(self, name: Optional[str]):
        """Outline module ``name``; None clears the selection."""
        if name != self._selected:
            self._selected = name
            self.update()

    def set_range(self, vmin: float, vmax: float):
        self._vmin = vmin
        self._vmax = vmax
        self.update()

    def set_palette(self, idx_or_name):
        """Set palette by index or name."""
        if isinstance(idx_or_name, str):
            idx = PALETTE_NAMES.index(idx_or_name)
        else:
            idx = int(idx_or_name)
        self._palette_idx = idx % len(PALETTES)
        self.update()

    def cycle_palette(self):
        self._palette_idx = (self._palette_idx + 1) % len(PALETTES)
        self.update()

    def set_log_scale(self, on: bool):
        self._log_scale = on
        self.update()

    def is_log_scale(self) -> bool:
        return self._log_scale

    def palette_stops(self):
        return list(PALETTES.values())[self._palette_idx]

    def reset_view(self):
        self._zoom = 1.0
        self._pan_x = 0.0
        self._pan_y = 0.0
        self._layout_dirty = True
        self.update()

    def value_to_t(self, v: float) -> float:
        """Map a raw value to [0, 1] using current scale (linear or log)."""
        vmin, vmax = self._vmin, self._vmax
        if self._log_scale:
            floor = max(vmin, self.LOG_FLOOR)
            ceil = max(vmax, floor * 10)
            v = max(v, floor)
            return (math.log10(v) - math.log10(floor)) / \
                   (math.log10(ceil) - math.log10(floor))
        return (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5

    # -- Layout --

    def _recompute_layout(self):
        self._rects.clear()
        if not self._modules:
            self._rect_names_rev = []
            self._layout_dirty = False
            return

        w, h = self.width(), self.height()
        margin, top, bot = self._margin, self._margin_top, self._margin_bottom
        pw, ph = w - 2 * margin, h - top - bot
        x0, x1, y0, y1 = self._geo_bounds
        base_scale = min(pw / max(x1 - x0, 1e-9), ph / max(y1 - y0, 1e-9))
        sc = base_scale * self._zoom
        dw, dh = (x1 - x0) * sc, (y1 - y0) * sc
        ox = margin + (pw - dw) / 2 + self._pan_x
        oy = top + (ph - dh) / 2 + self._pan_y

        # Record layout geometry (useful for subclass overlays)
        self._geo_x0 = x0
        self._geo_y1 = y1
        self._geo_scale = sc
        self._geo_ox = ox
        self._geo_oy = oy

        shrink = self._shrink
        for m in self._modules:
            mw, mh = m.sx * sc * shrink, m.sy * sc * shrink
            cx = ox + (m.x - x0) * sc
            cy = oy + (y1 - m.y) * sc
            self._rects[m.name] = QRectF(cx - mw / 2, cy - mh / 2, mw, mh)
        self._rect_names_rev = list(self._rects)[::-1]
        self._layout_dirty = False

    def geo_to_canvas(self, gx: float, gy: float) -> QPointF:
        """Convert geometry-space coords to widget canvas coords."""
        return QPointF(self._geo_ox + (gx - self._geo_x0) * self._geo_scale,
                       self._geo_oy + (self._geo_y1 - gy) * self._geo_scale)

    def resizeEvent(self, event):
        self._layout_dirty = True
        if self._reset_btn is not None:
            self._reset_btn.move(self.width() - self._reset_btn.width() - 6, 6)
        super().resizeEvent(event)

    # -- Painting --

    def paintEvent(self, event):
        if self._layout_dirty:
            self._recompute_layout()
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, False)
        w, h = self.width(), self.height()
        p.fillRect(0, 0, w, h, self.BG_COLOR)

        if not self._rects:
            self._paint_empty(p, w, h)
            p.end()
            return

        self._paint_before_modules(p, w, h)
        self._paint_modules(p)
        self._paint_overlays(p, w, h)
        if self._show_colorbar:
            self._paint_colorbar(p, w, h)
        self._paint_after_colorbar(p, w, h)
        p.end()

    # -- hook: empty state (no modules loaded) --
    def _paint_empty(self, p: QPainter, w: int, h: int):
        pass

    # -- hook: before modules (title etc.) --
    def _paint_before_modules(self, p: QPainter, w: int, h: int):
        pass

    # -- hook: per-module fill (default: colormap by value) --
    def _paint_modules(self, p: QPainter):
        stops = self.palette_stops()
        no_data = self.NO_DATA_COLOR
        for name, rect in self._rects.items():
            v = self._values.get(name)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                p.fillRect(rect, no_data)
            else:
                p.fillRect(rect, cmap_qcolor(self.value_to_t(v), stops))

    # -- hook: after modules, before colorbar --
    def _paint_overlays(self, p: QPainter, w: int, h: int):
        if self._labelled_names:
            p.setPen(QColor(THEME.TEXT))
            p.setFont(QFont("Monospace", 7, QFont.Weight.Bold))
            for name in self._labelled_names:
                r = self._rects.get(name)
                if r is not None:
                    p.drawText(r, Qt.AlignmentFlag.AlignCenter,
                               self._label_text(name))
        for name, colour, width in (
                (self._selected, self.SELECT_COLOR, self.SELECT_PEN_WIDTH),
                (self._hovered, self.HOVER_COLOR, 2.0)):
            if name and name in self._rects:
                p.setPen(QPen(colour, width))
                p.setBrush(Qt.BrushStyle.NoBrush)
                p.drawRect(self._rects[name])

    # -- hook: after colorbar (legend, extra labels) --
    def _paint_after_colorbar(self, p: QPainter, w: int, h: int):
        pass

    # -- hook: value format in colorbar min/max labels --
    def _fmt_value(self, v: float) -> str:
        return fmt_value(v)

    # -- hook: in-cell label of a module whose type is in label_types --
    def _label_text(self, name: str) -> str:
        return name

    # -- hook: maximum colour bar width --
    CB_MAX_WIDTH = 400

    def _paint_colorbar(self, p: QPainter, w: int, h: int):
        stops = self.palette_stops()
        cb_w = min(self.CB_MAX_WIDTH, w - 80)
        cb_h = 20
        cb_x = (w - cb_w) / 2
        cb_y = h - 50
        self._cb_rect = QRectF(cb_x, cb_y, cb_w, cb_h)

        grad = QLinearGradient(cb_x, 0, cb_x + cb_w, 0)
        for t, (r, g, b) in stops:
            grad.setColorAt(t, QColor(r, g, b))
        p.fillRect(self._cb_rect, QBrush(grad))
        p.setPen(QPen(self.CB_BORDER, 1.0))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(self._cb_rect)

        p.setPen(self.CB_TEXT)
        font = QFont("Consolas", 11)
        p.setFont(font)
        fm = QFontMetricsF(font)
        tick_h   = 5          # tick mark length (px)
        label_y  = cb_y + cb_h + tick_h + 1
        label_h  = fm.height() + 2
        min_str  = self._fmt_value(self._vmin)
        max_str  = self._fmt_value(self._vmax)

        # ── intermediate ticks ────────────────────────────────────────
        # Candidate positions evenly spaced in linear/log, rounded to int.
        # Duplicates, values equal to endpoints, and ticks whose labels
        # would overlap with a neighbour or with the endpoint labels are
        # all discarded.
        N_CANDS = 9     # generate more candidates than we may show
        vmin, vmax = self._vmin, self._vmax
        if self._log_scale and vmin > 0 and vmax > vmin:
            raw_vals = [10 ** (
                math.log10(vmin) + (i + 1) / (N_CANDS + 1) *
                (math.log10(vmax) - math.log10(vmin))
            ) for i in range(N_CANDS)]
        else:
            span = vmax - vmin if vmax > vmin else 1.0
            raw_vals = [vmin + (i + 1) / (N_CANDS + 1) * span
                        for i in range(N_CANDS)]

        # Round to integer and deduplicate while preserving order
        seen = set()
        int_vals = []
        for rv in raw_vals:
            iv = int(round(rv))
            if iv in seen or iv <= int(round(vmin)) or iv >= int(round(vmax)):
                continue
            seen.add(iv)
            int_vals.append(iv)

        # Convert values → pixel positions, then drop ticks whose label
        # would overlap a neighbour (min pixel gap = label_min_gap).
        def val_to_x(v):
            if self._log_scale and vmin > 0 and vmax > vmin:
                t = (math.log10(max(v, 1e-9)) - math.log10(max(vmin, 1e-9))) / \
                    (math.log10(max(vmax, 1e-9)) - math.log10(max(vmin, 1e-9)))
            else:
                t = (v - vmin) / (vmax - vmin) if vmax > vmin else 0.5
            return cb_x + t * cb_w

        label_min_gap = 8   # minimum pixel gap between label edges
        accepted = []   # list of (iv, tx, lbl_w)
        for iv in int_vals:
            lbl  = str(iv)
            lw   = fm.horizontalAdvance(lbl) + 4
            tx   = val_to_x(iv)
            # check against already-accepted ticks and endpoint labels
            ok = True
            # endpoint label widths (approx)
            end_w = 130
            if tx - lw / 2 < cb_x + end_w + label_min_gap:
                ok = False
            elif tx + lw / 2 > cb_x + cb_w - end_w - label_min_gap:
                ok = False
            else:
                for _, ptx, plw in accepted:
                    if abs(tx - ptx) < (lw + plw) / 2 + label_min_gap:
                        ok = False
                        break
            if ok:
                accepted.append((iv, tx, lw))

        p.setPen(QPen(self.CB_TEXT, 1.0))
        for iv, tx, lbl_w in accepted:
            p.drawLine(QPointF(tx, cb_y + cb_h),
                       QPointF(tx, cb_y + cb_h + tick_h))
            lbl = str(iv)
            p.drawText(QRectF(tx - lbl_w / 2, label_y, lbl_w, label_h),
                       Qt.AlignmentFlag.AlignCenter, lbl)

        # ── min / max labels ──────────────────────────────────────────
        p.setPen(self.CB_TEXT)
        p.drawText(QRectF(cb_x, label_y, 130, label_h),
                   Qt.AlignmentFlag.AlignLeft, min_str)
        p.drawText(QRectF(cb_x + cb_w - 130, label_y, 130, label_h),
                   Qt.AlignmentFlag.AlignRight, max_str)

        if self._enable_inline_range_edit:
            # Faint pencil glyph next to each editable label hints
            # "click to edit".  Hit rects cover both the value text and
            # the pencil so the user can click anywhere on either.
            pencil = "✎"   # LOWER RIGHT PENCIL
            pencil_w = fm.horizontalAdvance(pencil) + 4
            min_w = fm.horizontalAdvance(min_str) + 4
            max_w = fm.horizontalAdvance(max_str) + 4
            p.setPen(QColor(THEME.TEXT_DIM))
            self._cb_min_hit = QRectF(
                cb_x - 2, label_y - 1, min_w + pencil_w + 4, label_h + 2)
            p.drawText(QRectF(cb_x + min_w, label_y, pencil_w, label_h),
                       Qt.AlignmentFlag.AlignLeft, pencil)
            self._cb_max_hit = QRectF(
                cb_x + cb_w - max_w - pencil_w - 2, label_y - 1,
                max_w + pencil_w + 4, label_h + 2)
            p.drawText(QRectF(cb_x + cb_w - max_w - pencil_w, label_y,
                              pencil_w, label_h),
                       Qt.AlignmentFlag.AlignLeft, pencil)
        else:
            self._cb_min_hit = None
            self._cb_max_hit = None

    # -- Mouse / hit-test --

    def _hit(self, pos) -> Optional[str]:
        for name in self._rect_names_rev:
            if self._rects[name].contains(pos):
                return name
        return None

    def _tooltip_text(self, name: str) -> str:
        v = self._values.get(name)
        if v is None:
            return name
        return f"{name}: {self._fmt_value(v)}"

    def mousePressEvent(self, e):
        if self._enable_zoom_pan and e.button() == Qt.MouseButton.MiddleButton:
            self.reset_view()
            return
        if e.button() in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            self._drag_last = e.position()
            self._drag_origin = e.position()
            self._dragging = False

    def mouseReleaseEvent(self, e):
        if e.button() not in (Qt.MouseButton.LeftButton, Qt.MouseButton.RightButton):
            return
        if self._dragging:
            self.setCursor(Qt.CursorShape.ArrowCursor)
        elif e.button() == Qt.MouseButton.LeftButton:
            self._handle_click(e.position())
        self._drag_last = None
        self._drag_origin = None
        self._dragging = False

    def _handle_click(self, pos):
        """Default: vmin/vmax label hit → inline edit; colour-bar body →
        paletteClicked; else → moduleClicked.  With ``toggle_select`` the
        click toggles the selection and moduleClicked carries the new one
        ("" when cleared); an empty-canvas click with nothing selected
        emits nothing."""
        if self._check_inline_range_edit_click(pos):
            return
        if self._cb_rect and self._cb_rect.contains(pos):
            self.paletteClicked.emit()
            return
        name = self._hit(pos)
        if not self._toggle_select:
            self.moduleClicked.emit(name or "")
        elif name is not None or self._selected is not None:
            self.set_selected(None if name == self._selected else name)
            self.moduleClicked.emit(self._selected or "")

    def _check_inline_range_edit_click(self, pos) -> bool:
        """Hit-test ``pos`` against the colorbar vmin/vmax labels and open
        the inline editor if it lands on one.  Returns True if handled.

        Subclasses that override ``mousePressEvent`` (e.g. paint editors
        that intercept clicks on press rather than release) should call
        this at the top of their override to keep the inline-edit feature
        working — otherwise their override swallows clicks on the labels.
        """
        if not self._enable_inline_range_edit:
            return False
        if self._cb_min_hit is not None and self._cb_min_hit.contains(pos):
            self._show_inline_editor("min", self._cb_min_hit)
            return True
        if self._cb_max_hit is not None and self._cb_max_hit.contains(pos):
            self._show_inline_editor("max", self._cb_max_hit)
            return True
        return False

    def mouseMoveEvent(self, e):
        # zoom/pan drag
        if self._enable_zoom_pan and self._drag_last is not None:
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

        # hover
        pos = e.position()
        if self._enable_inline_range_edit and (
                (self._cb_min_hit is not None and self._cb_min_hit.contains(pos))
                or (self._cb_max_hit is not None and self._cb_max_hit.contains(pos))):
            self.setCursor(Qt.CursorShape.IBeamCursor)
        elif self._cb_rect and self._cb_rect.contains(pos):
            self.setCursor(Qt.CursorShape.PointingHandCursor)
        else:
            self.setCursor(Qt.CursorShape.ArrowCursor)

        found = self._hit(pos)
        if found != self._hovered:
            self._hovered = found
            self.update()
            if found:
                QToolTip.showText(e.globalPosition().toPoint(),
                                  self._tooltip_text(found), self)
                self.moduleHovered.emit(found)
            else:
                QToolTip.hideText()

    def wheelEvent(self, e):
        if not self._enable_zoom_pan:
            return
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

    # -- Inline range edit (colorbar vmin / vmax labels) --

    def _show_inline_editor(self, which: str, rect: QRectF):
        """Open a child QLineEdit over the clicked vmin/vmax label."""
        if self._inline_editor is None:
            self._inline_editor = QLineEdit(self)
            self._inline_editor.setValidator(
                QDoubleValidator(-1e12, 1e12, 6, self._inline_editor))
            self._inline_editor.setStyleSheet(
                f"QLineEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
                f"border:1px solid {THEME.ACCENT_STRONG};border-radius:3px;"
                f"padding:1px 4px;font:9pt Consolas;}}")
            self._inline_editor.installEventFilter(self)
            self._inline_editor.editingFinished.connect(self._commit_inline_edit)
        self._inline_which = which
        self._inline_cancelled = False
        cur = self._vmin if which == "min" else self._vmax
        self._inline_editor.setText(self._fmt_value(cur))
        ed_w = max(80, int(rect.width()))
        ed_h = max(20, int(rect.height()) + 4)
        if which == "max":
            ed_x = int(rect.right()) - ed_w
        else:
            ed_x = int(rect.left())
        ed_y = int(rect.top()) - 2
        self._inline_editor.setGeometry(ed_x, ed_y, ed_w, ed_h)
        self._inline_editor.show()
        self._inline_editor.raise_()
        self._inline_editor.selectAll()
        self._inline_editor.setFocus()

    def _commit_inline_edit(self):
        # editingFinished fires on Enter and on focus-loss.  Esc handler
        # sets _inline_cancelled and hides the editor; bail in that case.
        if self._inline_editor is None or self._inline_which is None:
            return
        if self._inline_cancelled:
            self._inline_which = None
            return
        try:
            new_val = float(self._inline_editor.text())
        except ValueError:
            self._inline_editor.hide()
            self._inline_which = None
            return
        if self._inline_which == "min":
            new_min, new_max = new_val, self._vmax
        else:
            new_min, new_max = self._vmin, new_val
        if (math.isfinite(new_min) and math.isfinite(new_max)
                and new_max > new_min):
            self.set_range(new_min, new_max)
            self.rangeEdited.emit(new_min, new_max)
        self._inline_editor.hide()
        self._inline_which = None

    def eventFilter(self, obj, event):
        if obj is self._inline_editor and event.type() == QEvent.Type.KeyPress:
            if event.key() == Qt.Key.Key_Escape:
                self._inline_cancelled = True
                self._inline_editor.hide()
                self._inline_which = None
                return True
        return super().eventFilter(obj, event)

    def sizeHint(self):
        return QSize(680, 680)


# ---- Colormap-range control ----
#
# Auto-button gestures of ColorRangeControl:
#   * Single click               → one-shot fit; pin state unchanged.
#   * Double click               → one-shot fit + enter persistent mode
#                                  (button highlights with ACCENT_STRONG).
#   * Click while pinned         → exit persistent mode.
#   * Double-click while pinned  → exit persistent mode (self-cancelling).
#   * Editing a range field      → exits persistent mode automatically.


class _AutoButton(QPushButton):
    """QPushButton that distinguishes single-click from double-click.

    Qt fires ``clicked`` for both presses of a double-click.  We use a
    counter + ``QApplication.doubleClickInterval()`` timer to disambiguate:
    the first ``clicked`` starts the timer; if a second ``clicked`` arrives
    before the timer expires, it's a double-click; otherwise single.
    """

    oneshotRequested  = pyqtSignal()
    pinToggleRequested = pyqtSignal()

    def __init__(self, text: str = "Auto", parent: Optional[QWidget] = None):
        super().__init__(text, parent)
        self._click_count = 0
        self._timer = QTimer(self)
        self._timer.setSingleShot(True)
        self._timer.setInterval(QApplication.doubleClickInterval())
        self._timer.timeout.connect(self._fire_pending)
        self.clicked.connect(self._on_clicked)

    def _on_clicked(self):
        self._click_count += 1
        if self._click_count == 1:
            self._timer.start()
        else:
            self._timer.stop()
            self._click_count = 0
            self.pinToggleRequested.emit()

    def _fire_pending(self):
        if self._click_count == 1:
            self._click_count = 0
            self.oneshotRequested.emit()
        else:
            self._click_count = 0


class ColorRangeController(QObject):
    """Headless controller for a HyCalMapWidget colormap range.

    Holds auto-fit / pin / log behaviour and exposes signals and methods
    so scripts can wire whatever GUI they want — a single "Auto" toolbar
    button, a menu action, the inline edits on the colorbar, or no GUI at
    all.

    Parameters
    ----------
    target
        The ``HyCalMapWidget`` whose range / log scale it drives.
    auto_fit
        How a fit computes the range from the current data dict:
        ``"minmax"`` (full ``min..max``), ``"minmax_nonzero"``
        (``min..max`` of values != 0, for zero as a sentinel such as masked
        channels), ``"percentile"`` (``auto_fit_percentile``), or a
        callable ``f(values) -> (vmin, vmax)``.
    auto_fit_percentile
        ``(lo, hi)`` for the ``"percentile"`` preset.  Default ``(2, 98)``.

    Signals
    -------
    rangeChanged(vmin, vmax) — any range change (auto-fit, set_range, or
                               relayed widget inline edit).
    autoPinned(on)           — persistent auto mode flipped.
    logToggled(on)           — log/linear flipped.

    Public API
    ----------
    notify_values_changed(values)  — call when the data dict changes
                                     (re-fits if pinned).
    auto_fit(values=None)          — one-shot programmatic fit.
    set_range(vmin, vmax)          — push a range; turns off pin.
    set_pinned(on) / is_pinned()
    set_log(on)
    Properties: vmin, vmax.
    """

    rangeChanged = pyqtSignal(float, float)
    autoPinned   = pyqtSignal(bool)
    logToggled   = pyqtSignal(bool)

    _AUTO_FIT_PRESETS = ("minmax", "minmax_nonzero", "percentile")
    AutoFit = Union[str, Callable[[Dict[str, float]], Tuple[float, float]]]

    def __init__(self,
                 target: HyCalMapWidget,
                 *,
                 auto_fit: AutoFit = "minmax",
                 auto_fit_percentile: Tuple[float, float] = (2.0, 98.0),
                 parent: Optional[QObject] = None):
        super().__init__(parent)

        if not isinstance(target, HyCalMapWidget):
            raise TypeError(
                "ColorRangeController target must be a HyCalMapWidget")
        self._map = target

        if not (callable(auto_fit) or auto_fit in self._AUTO_FIT_PRESETS):
            raise ValueError(
                f"auto_fit must be callable or one of {self._AUTO_FIT_PRESETS}; "
                f"got {auto_fit!r}")
        self._auto_fit = auto_fit
        self._auto_pct = auto_fit_percentile

        self._pinned = False
        self._values: Dict[str, float] = {}

        # Relay widget inline edits — they always override pin (manual edit).
        self._map.rangeEdited.connect(self._on_widget_range_edited)

    def _read_target_range(self) -> Tuple[float, float]:
        return self._map._vmin, self._map._vmax

    # ---- public API ------------------------------------------------------

    def notify_values_changed(self, values: Dict[str, float]):
        """Call when the data dict changes; re-fits if pinned."""
        self._values = values or {}
        if self._pinned:
            self._do_auto_fit_and_apply()

    def auto_fit(self, values: Optional[Dict[str, float]] = None):
        """One-shot fit.  Pin state unchanged.  ``values`` overrides cache."""
        if values is not None:
            self._values = values
        self._do_auto_fit_and_apply()

    def set_range(self, vmin: float, vmax: float):
        """Push a range; turns off pin."""
        if not (math.isfinite(vmin) and math.isfinite(vmax)) or vmax <= vmin:
            return
        if self._pinned:
            self._set_pinned(False)
        self._map.set_range(vmin, vmax)
        self.rangeChanged.emit(vmin, vmax)

    def set_pinned(self, on: bool):
        self._set_pinned(bool(on))

    def is_pinned(self) -> bool:
        return self._pinned

    def set_log(self, on: bool):
        on = bool(on)
        self._map.set_log_scale(on)
        self.logToggled.emit(on)

    @property
    def vmin(self) -> float:
        return self._read_target_range()[0]

    @property
    def vmax(self) -> float:
        return self._read_target_range()[1]

    # ---- internal --------------------------------------------------------

    def _on_widget_range_edited(self, vmin: float, vmax: float):
        """Inline edit on the widget's colorbar — always exits pin."""
        if self._pinned:
            self._set_pinned(False)
        self.rangeChanged.emit(vmin, vmax)

    def _set_pinned(self, on: bool):
        if self._pinned == on:
            return
        self._pinned = on
        self.autoPinned.emit(on)

    def _do_auto_fit_and_apply(self):
        vmin, vmax = self._compute_auto_fit()
        if not math.isfinite(vmin):
            return
        if not math.isfinite(vmax) or vmax <= vmin:
            vmax = vmin + max(abs(vmin) * 0.05, 1e-6)
        self._map.set_range(vmin, vmax)
        self.rangeChanged.emit(vmin, vmax)

    def _compute_auto_fit(self) -> Tuple[float, float]:
        values = self._values
        if callable(self._auto_fit):
            return tuple(self._auto_fit(values))
        if not values:
            return self._read_target_range()
        if self._auto_fit == "percentile":
            try:
                import numpy as np
            except ImportError:     # fall back to min..max below
                np = None
            if np is not None:
                arr = np.asarray(list(values.values()), dtype=float)
                arr = arr[np.isfinite(arr)]
                if arr.size == 0:
                    return 0.0, 1.0
                lo, hi = self._auto_pct
                return (float(np.percentile(arr, lo)),
                        float(np.percentile(arr, hi)))
        skip_zero = self._auto_fit == "minmax_nonzero"
        vals = [v for v in values.values()
                if v is not None and not (isinstance(v, float) and math.isnan(v))
                and not (skip_zero and v == 0.0)]
        if not vals:
            return 0.0, 1.0
        return float(min(vals)), float(max(vals))


def _toggle_btn_qss(active: bool, accent: str, idle_fg: str) -> str:
    """QSS of a toggle button: filled with ``accent`` when active, else a
    plain button with ``idle_fg`` text."""
    if active:
        return (f"QPushButton{{background:{accent};color:{THEME.TEXT};"
                f"border:1px solid {accent};padding:5px 14px;"
                f"font:10pt;border-radius:6px;}}")
    return (f"QPushButton{{background:{THEME.BUTTON};color:{idle_fg};"
            f"border:1px solid {THEME.BORDER};padding:5px 14px;"
            f"font:10pt;border-radius:6px;}}"
            f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};"
            f"color:{THEME.TEXT};}}")


class ColorRangeControl(QWidget):
    """Default min/max + Auto button + Log toggle widget.

    Thin shim over :class:`ColorRangeController` that provides the
    classic two-edit row UI.  Construct it for the simple case; for
    custom UIs (single button, inline-only, …) build a
    ``ColorRangeController`` directly and wire your own widgets.

    See :class:`ColorRangeController` for the target / auto-fit / pin
    semantics; this widget passes those parameters straight through.

    Extra parameters
    ----------------
    include_log    — add an inline Log toggle button.
    orientation    — ``"horizontal"`` or ``"vertical"``.
    start_pinned   — start in persistent auto-fit mode.

    The underlying controller is exposed as ``self.controller`` for
    advanced wiring.
    """

    AutoFit = ColorRangeController.AutoFit

    def __init__(self,
                 target,
                 *,
                 auto_fit: AutoFit = "minmax",
                 auto_fit_percentile: Tuple[float, float] = (2.0, 98.0),
                 include_log: bool = False,
                 orientation: str = "horizontal",
                 start_pinned: bool = False,
                 parent: Optional[QWidget] = None):
        super().__init__(parent)
        self._ctrl = ColorRangeController(
            target,
            auto_fit=auto_fit,
            auto_fit_percentile=auto_fit_percentile,
            parent=self,
        )
        self._build_ui(orientation, include_log)
        self._ctrl.rangeChanged.connect(self._set_edits)
        self._ctrl.autoPinned.connect(
            lambda _on: self._update_auto_btn_style())
        self._ctrl.logToggled.connect(self._sync_log_btn)
        self._set_edits(self._ctrl.vmin, self._ctrl.vmax)
        if start_pinned:
            self._ctrl.set_pinned(True)

    # ---- accessors -------------------------------------------------------

    @property
    def controller(self) -> ColorRangeController:
        return self._ctrl

    def notify_values_changed(self, values: Dict[str, float]):
        self._ctrl.notify_values_changed(values)

    def set_range(self, vmin: float, vmax: float):
        self._ctrl.set_range(vmin, vmax)

    # ---- UI construction --------------------------------------------------

    def _make_range_edit(self) -> QLineEdit:
        ed = QLineEdit()
        ed.setMaximumWidth(90)
        ed.setValidator(QDoubleValidator(-1e12, 1e12, 6, ed))
        ed.editingFinished.connect(self._on_edit)
        ed.setStyleSheet(
            f"QLineEdit{{background:{THEME.PANEL};color:{THEME.TEXT};"
            f"border:1px solid {THEME.BORDER};border-radius:4px;"
            f"padding:2px 6px;}}")
        return ed

    def _build_ui(self, orientation: str, include_log: bool):
        self._min_edit = self._make_range_edit()
        self._max_edit = self._make_range_edit()

        self._auto_btn = _AutoButton("Auto", self)
        self._auto_btn.setToolTip(
            "Click: auto-fit once   ·   Double-click: keep auto-fitting")
        self._auto_btn.oneshotRequested.connect(lambda: self._on_auto(False))
        self._auto_btn.pinToggleRequested.connect(lambda: self._on_auto(True))
        self._update_auto_btn_style()

        self._log_btn: Optional[QPushButton] = None
        if include_log:
            self._log_btn = QPushButton("Log")
            self._log_btn.setCheckable(True)
            self._log_btn.toggled.connect(self._on_log_clicked)
            self._update_log_btn_style()

        self.setStyleSheet(
            f"QLabel{{color:{THEME.TEXT};background:transparent;}}")

        if orientation == "vertical":
            outer = QVBoxLayout(self)
            outer.setContentsMargins(0, 0, 0, 0)
            outer.setSpacing(4)
            r1 = QHBoxLayout(); r1.setSpacing(4)
            r1.addWidget(QLabel("min:"))
            r1.addWidget(self._min_edit); r1.addStretch()
            outer.addLayout(r1)
            r2 = QHBoxLayout(); r2.setSpacing(4)
            r2.addWidget(QLabel("max:"))
            r2.addWidget(self._max_edit); r2.addStretch()
            r3 = QHBoxLayout(); r3.setSpacing(6)
            r3.addWidget(self._auto_btn)
            if self._log_btn is not None:
                r3.addWidget(self._log_btn)
            r3.addStretch()
            outer.addLayout(r2); outer.addLayout(r3)
        else:
            row = QHBoxLayout(self)
            row.setContentsMargins(0, 0, 0, 0)
            row.setSpacing(6)
            row.addWidget(QLabel("Range:"))
            row.addWidget(self._min_edit)
            row.addWidget(QLabel("–"))
            row.addWidget(self._max_edit)
            row.addWidget(self._auto_btn)
            if self._log_btn is not None:
                row.addWidget(self._log_btn)
            row.addStretch()

    def _update_auto_btn_style(self):
        self._auto_btn.setStyleSheet(_toggle_btn_qss(
            self._ctrl.is_pinned(), THEME.ACCENT_STRONG, THEME.TEXT))

    def _update_log_btn_style(self):
        self._log_btn.setStyleSheet(_toggle_btn_qss(
            self._log_btn.isChecked(), THEME.ACCENT, THEME.TEXT_DIM))

    def _sync_log_btn(self, on: bool):
        """Show ``on`` on the Log button without re-emitting its toggle."""
        if self._log_btn is None or self._log_btn.isChecked() == on:
            return
        self._log_btn.blockSignals(True)
        self._log_btn.setChecked(on)
        self._log_btn.blockSignals(False)
        self._update_log_btn_style()

    def _set_edits(self, vmin: float, vmax: float):
        self._min_edit.blockSignals(True)
        self._max_edit.blockSignals(True)
        self._min_edit.setText(f"{vmin:.6g}")
        self._max_edit.setText(f"{vmax:.6g}")
        self._min_edit.blockSignals(False)
        self._max_edit.blockSignals(False)

    # ---- handlers --------------------------------------------------------

    def _on_edit(self):
        try:
            vmin = float(self._min_edit.text())
            vmax = float(self._max_edit.text())
        except ValueError:
            return
        self._ctrl.set_range(vmin, vmax)

    def _on_auto(self, pin: bool):
        """Auto click (``pin`` False) or double click: exit persistent mode
        when pinned, else fit once and, on a double click, pin."""
        if self._ctrl.is_pinned():
            self._ctrl.set_pinned(False)
        else:
            self._ctrl.auto_fit()
            if pin:
                self._ctrl.set_pinned(True)

    def _on_log_clicked(self, on: bool):
        self._update_log_btn_style()
        self._ctrl.set_log(on)


# ---- Config tuning docks ----
#
# Tuning docks edit the fields of a prad2py config object (WaveConfig,
# HyCalClusterConfig, GEM ClusterConfig) through one editor per field,
# kept in a {field name: editor} dict.  Editors are QSpinBox (int field),
# QDoubleSpinBox (float), QCheckBox (bool) or QComboBox (int index).

def setup_tuning_dock(dock: QDockWidget) -> QVBoxLayout:
    """Dock at the left or right, closable / movable / floatable, with a
    fresh content widget whose layout is returned."""
    dock.setAllowedAreas(Qt.DockWidgetArea.LeftDockWidgetArea
                         | Qt.DockWidgetArea.RightDockWidgetArea)
    dock.setFeatures(QDockWidget.DockWidgetFeature.DockWidgetClosable
                     | QDockWidget.DockWidgetFeature.DockWidgetMovable
                     | QDockWidget.DockWidgetFeature.DockWidgetFloatable)
    root = QWidget()
    dock.setWidget(root)
    layout = QVBoxLayout(root)
    layout.setContentsMargins(6, 6, 6, 6)
    return layout


def add_config_rows(form: QFormLayout, cfg, specs,
                    on_change: Callable) -> Dict[str, QWidget]:
    """Add a spin-box row to ``form`` per ``(name, lo, hi, step, tip[,
    label])`` in ``specs`` and return the editors by name.

    ``step`` None gives a QSpinBox, otherwise a 3-decimal QDoubleSpinBox.
    The row label defaults to ``name``; the start value is
    ``getattr(cfg, name)`` when ``cfg`` has the field.  Every value change
    calls ``on_change``.
    """
    editors: Dict[str, QWidget] = {}
    for name, lo, hi, step, tip, *label in specs:
        if step is None:
            ed = QSpinBox()
            ed.setRange(lo, hi)
        else:
            ed = QDoubleSpinBox()
            ed.setRange(lo, hi)
            ed.setSingleStep(step)
            ed.setDecimals(3)
        if cfg is not None and hasattr(cfg, name):
            set_editor_value(ed, getattr(cfg, name))
        ed.setToolTip(tip)
        ed.valueChanged.connect(on_change)
        form.addRow(label[0] if label else name, ed)
        editors[name] = ed
    return editors


def editor_value(ed: QWidget):
    """Value of a config editor, typed for its field."""
    if isinstance(ed, QCheckBox):
        return ed.isChecked()
    if isinstance(ed, QComboBox):
        return ed.currentIndex()
    if isinstance(ed, QSpinBox):
        return int(ed.value())
    return float(ed.value())


def set_editor_value(ed: QWidget, v) -> None:
    """Show ``v`` in a config editor without emitting its change signal."""
    ed.blockSignals(True)
    try:
        if isinstance(ed, QCheckBox):
            ed.setChecked(bool(v))
        elif isinstance(ed, QComboBox):
            ed.setCurrentIndex(int(v))
        elif isinstance(ed, QSpinBox):
            ed.setValue(int(v))
        else:
            ed.setValue(float(v))
    finally:
        ed.blockSignals(False)


def config_to_editors(cfg, editors: Dict[str, QWidget]) -> None:
    """Show the fields of ``cfg`` in their editors (fields it lacks are
    left alone), without change signals."""
    for name, ed in editors.items():
        if hasattr(cfg, name):
            set_editor_value(ed, getattr(cfg, name))


def editors_to_config(editors: Dict[str, QWidget], cfg) -> None:
    """Write the editor values into ``cfg``.  A field the binding rejects
    is reported on stderr and skipped."""
    for name, ed in editors.items():
        try:
            setattr(cfg, name, editor_value(ed))
        except Exception as exc:  # noqa: BLE001
            print(f"[config] {type(cfg).__name__}.{name} setattr failed: "
                  f"{exc}", file=sys.stderr)


# ---- Background workers ----

def start_worker_thread(owner: QObject, worker: QObject,
                        on_finished: Callable, on_failed: Callable, *,
                        dialog=None,
                        on_thread_finished: Optional[Callable] = None
                        ) -> QThread:
    """Run ``worker.run`` in a new QThread owned by ``owner`` and start it.

    ``worker`` is a QObject with a ``run()`` slot, ``finished`` / ``failed``
    signals and a ``request_cancel()`` method.  Connect its progress signal
    before calling.  The thread quits when the worker finishes or fails,
    after which both are deleted and ``on_thread_finished`` runs.
    ``dialog.canceled`` (a QProgressDialog) requests cancellation.  The
    worker is unparented, so the caller must keep a Python reference to it
    until the thread has finished.
    """
    thread = QThread(owner)
    worker.moveToThread(thread)
    thread.started.connect(worker.run)
    worker.finished.connect(on_finished)
    worker.failed.connect(on_failed)
    if dialog is not None:
        dialog.canceled.connect(worker.request_cancel)
    worker.finished.connect(thread.quit)
    worker.failed.connect(thread.quit)
    thread.finished.connect(worker.deleteLater)
    thread.finished.connect(thread.deleteLater)
    if on_thread_finished is not None:
        thread.finished.connect(on_thread_finished)
    thread.start()
    return thread
