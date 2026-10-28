#!/usr/bin/env python3
"""Replay Viewer (PyQt6)
=======================
PyQt6 GUI for the PRad-2 online replay pipeline.

Quick Start
-----------
1. Click the "Do It All" button and fill in the popup:
     • Run number       — 6-digit run number, e.g. 024388
     • File index range — range of EVIO file indices to download (default 0–99)
     • Threads          — number of parallel replay threads (default 25)
     • Filter cut JSON  — slow-control cut config for replay_filter (optional)
     • GEM zero-sup     — GEM zero-suppression threshold (default 5)
   Click "Start Pipeline" to run the full pipeline:
     SCP → Replay Recon → hadd → Replay Filter → Quick Check
   Result plots appear automatically in the right panel.

Step-by-Step
------------
   Each step can also be run individually:
     "1. Get Data"       — download EVIO files from clondaq2
     "2. Run Replay"     — run prad2ana_replay_recon
     "3. Replay Filter"  — run prad2ana_replay_filter (also writes a JSON report)
     "4. Quick Check"    — run prad2ana_quick_check
   Detailed parameters for each step are available via the "⚙ Settings…" button.

Filter Report Tab
-----------------
   After Replay Filter completes, the "Filter Report" tab on the right
   automatically displays time-series plots (cut status / livetime+rate / EPICS
   channels).  Use the controls at the top to switch the x-axis (timestamp /
   event number) and toggle individual EPICS channels.

Other Features
--------------
   • "Auto-delete EVIO files" — delete the local EVIO directory after a
                                successful replay
   • "Check Disk"             — estimate disk usage before downloading
   • "Stop"                   — abort the currently running step
   • File → Open ROOT file…   — manually open an existing quick_check ROOT file

Command-Line Usage at Counting House Computer
------------------
   ssh clonfarm11
   source /home/clasrun/prad2_daq/prad2_env.csh
   python /data/soft/prad2evviewer/scripts/replay_viewer.py [quick_check.root]
"""
from __future__ import annotations

import html
import math
import os
import re
import shutil
import sys
from typing import Dict, List, Optional, Tuple

from prad2_env import fix_qt_lib_path

fix_qt_lib_path()   # before the first PyQt6 import

from PyQt6.QtCore import (
    QPointF, QProcess, QProcessEnvironment, QRectF, QThread, Qt, pyqtSignal,
)
from PyQt6.QtGui import (
    QColor, QFont, QImage, QPainter, QPen, QPixmap,
)
from PyQt6.QtWidgets import (
    QApplication, QCheckBox, QComboBox,
    QDialog, QFileDialog, QFormLayout, QHBoxLayout,
    QLabel, QLineEdit, QMainWindow, QMenu, QPushButton,
    QSizePolicy, QSplitter, QSpinBox,
    QTabWidget, QTextEdit, QToolButton, QVBoxLayout, QWidget,
)

# Optional ROOT file reading
try:
    import numpy as np
    import uproot
    HAS_UPROOT = True
except ImportError:
    HAS_UPROOT = False

# Optional matplotlib (for embedded filter report chart)
try:
    from matplotlib import rc_context
    from matplotlib.backends.backend_qtagg import (
        FigureCanvasQTAgg, NavigationToolbar2QT as _MplNavToolbar,
    )
    from matplotlib.figure import Figure as MplFigure
    HAS_MATPLOTLIB = True
except ImportError:
    HAS_MATPLOTLIB = False

# Shared HyCal infrastructure
from evio_io import (  # noqa: E402
    LOCAL_DATA_BASE, REMOTE_DATA_BASE, REMOTE_HOST,
    check_disk_space, fmt_bytes, local_evio_in_range, scp_bash,
)
from hycal_geoview import (  # noqa: E402
    MENU_QSS, PALETTES, PALETTE_NAMES, THEME, HyCalMapWidget, ZoomHistWidget,
    apply_theme_palette, available_themes, cmap_rgb_array,
    load_modules, nice_ticks, set_theme, themed,
)

# Path constants
SCRIPT_DIR   = os.path.dirname(os.path.abspath(__file__))
DB_DIR       = os.path.join(os.path.dirname(SCRIPT_DIR), "database")
MODULES_JSON = os.path.join(DB_DIR, "hycal_map.json")

# Optional replay_report_viewer (for embedded filter report chart)
try:
    from replay_report_viewer import (  # noqa: E402
        checked_epics   as _rrv_checked_epics,
        fill_epics_menu as _rrv_fill_epics_menu,
        load_report     as _rrv_load_report,
        plot_report     as _rrv_plot_report,
        selected_series as _rrv_selected_series,
        title_for       as _rrv_title_for,
    )
    HAS_REPORT_VIEWER = True
except ImportError:
    HAS_REPORT_VIEWER = False

# Replay tools — match the names used by other scripts
_PRAD2_BIN_DIR    = "/data/soft/prad2evviewer/build/bin"
_REPLAY_RECON_CMD   = os.path.join(_PRAD2_BIN_DIR, "prad2ana_replay_recon")
_REPLAY_FILTER_CMD  = os.path.join(_PRAD2_BIN_DIR, "prad2ana_replay_filter")
_QUICK_CHECK_CMD    = os.path.join(_PRAD2_BIN_DIR, "prad2ana_quick_check")
_RECON_BASE       = "/data/replay_recon"


# ---- Helpers ----

def _run_tag(run: str) -> Optional[str]:
    """'24100' -> 'prad_024100'; None if blank or not an integer."""
    try:
        return f"prad_{int(run):06d}"
    except ValueError:
        return None


def _run_dirs(run: str, host: str, remote_base: str, local_base: str
              ) -> Optional[Tuple[str, str, str, str, str]]:
    """(run tag, host, remote run dir, local base, local run dir) for the
    SCP fields, blank fields taking the defaults; None without a valid run."""
    tag = _run_tag(run)
    if tag is None:
        return None
    local_base = local_base.strip() or LOCAL_DATA_BASE
    return (tag, host.strip() or REMOTE_HOST,
            f"{remote_base.strip() or REMOTE_DATA_BASE}/{tag}",
            local_base, os.path.join(local_base, tag))


def _disk_summary_html(needed: int, free: int, prefix: str = "") -> str:
    ok = free >= needed
    return (f"<span style='color:{'#3fb950' if ok else '#f85149'}'>{prefix}"
            f"need {fmt_bytes(needed)}, free {fmt_bytes(free)}"
            f"{'  ✓' if ok else '  ✗ (insufficient)'}</span>")


def _disk_check_html(dirs, f_start: int, f_end: int) -> str:
    """Disk-space check of files [f_start, f_end] for ``dirs`` from
    _run_dirs(), as label HTML."""
    _tag, host, remote_run_dir, local_base, local_run_dir = dirs
    try:
        needed, free = check_disk_space(host, remote_run_dir, local_base,
                                        f_start, f_end, local_run_dir)
    except RuntimeError as exc:
        return f"<span style='color:#f85149'>SSH error: {exc}</span>"
    except Exception as exc:
        return f"<span style='color:#f85149'>Error: {exc}</span>"
    return _disk_summary_html(needed, free)


def _opt(cmd: List[str], flag: str, value: str, skip: Optional[str] = None):
    """Append ``flag value`` to ``cmd`` unless the value is blank or ``skip``."""
    value = value.strip()
    if value and value != skip:
        cmd += [flag, value]


# ---- Background worker ----

class _RootLoader(QThread):
    """Load quick_check ROOT file off the UI thread."""

    finished = pyqtSignal(dict, str)   # data_dict, error_message

    def __init__(self, path: str, parent=None):
        super().__init__(parent)
        self._path = path

    def run(self):
        if not HAS_UPROOT:
            self.finished.emit({}, "uproot not installed. Run: pip install uproot numpy")
            return
        try:
            data = _load_root(self._path)
            self.finished.emit(data, "")
        except Exception as exc:
            self.finished.emit({}, str(exc))


# (histogram path in the quick_check file, ResultsPanel widget attribute)
_HIST_PLOTS = (
    ("hit_pos",                            "_hit_pos_widget"),
    ("energy_plots/one_cluster_energy",    "_h_1cl"),
    ("energy_plots/two_cluster_energy",    "_h_2cl"),
    ("energy_plots/clusters_energy",       "_h_all"),
    ("energy_plots/total_energy",          "_h_tot"),
    ("energy_plots/h2_energy_theta",       "_ev_theta_widget"),
    ("physics_yields/ep_yield",            "_h_ep"),
    ("physics_yields/ee_yield",            "_h_ee"),
    ("physics_yields/yield_ratio",         "_h_ratio"),
    ("moller_analysis/h_moller_z",         "_h_moller_z"),
    ("moller_analysis/h_moller_phi_diff",  "_h_moller_phi"),
    ("moller_analysis/h_moller_x",         "_h_moller_x"),
    ("moller_analysis/h_moller_y",         "_h_moller_y"),
    ("moller_analysis/h2_moller_pos",      "_moller_2arm"),
)


def _load_root(path: str) -> dict:
    """Read a quick_check ROOT output file into numpy arrays."""
    out: dict = {}
    with uproot.open(path) as f:
        # (values, edges) of a 1-D, (values, x edges, y edges) of a 2-D
        # histogram; one that cannot be read is skipped.
        for key, _attr in _HIST_PLOTS:
            obj = f.get(key)
            if obj is None:
                continue
            try:
                out[key] = (obj.values().tolist(),
                            *(ax.edges().tolist() for ax in obj.axes))
            except Exception:
                pass

        # ---- module_energy: per-module hit counts and mean energies ----
        me = f.get("module_energy")
        if me is not None:
            module_counts: Dict[str, float] = {}
            module_means:  Dict[str, float] = {}
            for name in me.keys(cycle=False):
                h = me[name]
                if not hasattr(h, "values"):
                    continue
                # strip leading "h_" → "W432" matches hycal_map name "W432"
                mod_key = name[2:] if name.startswith("h_") else name
                vals  = h.values()
                total = float(vals.sum())
                module_counts[mod_key] = total
                if total > 0:
                    edges = h.axis().edges()
                    mids  = (edges[:-1] + edges[1:]) / 2.0
                    module_means[mod_key] = float((np.asarray(vals) * mids).sum() / total)
                else:
                    module_means[mod_key] = 0.0
            out["module_counts"] = module_counts
            out["module_means"]  = module_means

    return out


# ---- 1-D histogram widget ----

_LOG_BTN_QSS = (
    "QPushButton{background:#21262d;color:#8b949e;border:1px solid #30363d;"
    "border-radius:3px;font:8pt Consolas;padding:0 3px;}"
    "QPushButton:checked{background:#1f6feb;color:#fff;border-color:#388bfd;}"
    "QPushButton:hover{border-color:#58a6ff;color:#c9d1d9;}")


def _log_toggle(parent: QWidget, text: str, slot) -> QPushButton:
    """Small checkable log-scale button floating over a plot."""
    btn = QPushButton(text, parent)
    btn.setCheckable(True)
    btn.setFixedSize(36, 18)
    btn.setStyleSheet(_LOG_BTN_QSS)
    btn.clicked.connect(slot)
    return btn


class Hist1DWidget(ZoomHistWidget):
    """Lightweight 1-D histogram display with zoom drag and right-click unzoom."""

    BAR_COLOR = "#3fb950"
    BAR_GAP = 0.5
    X_LABEL_W = 56
    X_LABEL_FMT = ".4g"

    def __init__(self, title: str = "", parent=None):
        super().__init__(parent)
        self._default_title = title
        self._title = title
        self._log_y = False
        self.setMinimumSize(200, 120)
        self._btn_log_x = _log_toggle(self, "logX", self._toggle_log_x)
        self._btn_log_y = _log_toggle(self, "logY", self._toggle_log_y)
        self._cache_pm: Optional[QPixmap] = None
        self._plotted = False   # the cached pixmap holds a plot (not a placeholder)
        # Crystal Ball auto-fit
        self.auto_cb_fit: bool = False
        # Asymmetric fit window around peak: (left_width, right_width) in data units.
        # None means use the full histogram range.
        self.cb_fit_range: Optional[Tuple[float, float]] = None
        # Asymmetric fit window in units of estimated sigma: (left_nsigma, right_nsigma).
        # Takes priority over cb_fit_range when set.
        self.cb_fit_range_sigma: Optional[Tuple[float, float]] = None
        self._cb_fit_result: Optional[Tuple[float, float, float, float]] = None  # (mean, mean_err, sigma, sigma_err)
        self._cb_fit_curve: Optional[Tuple[List[float], List[float]]] = None  # (xs, ys)

    def resizeEvent(self, event):
        super().resizeEvent(event)
        w = self.width()
        self._btn_log_y.move(w - 42, 3)
        self._btn_log_x.move(w - 80, 3)
        self._cache_pm = None

    def set_data(self, values: List[float], edges: List[float], title: str = ""):
        self._values = list(values)
        self._edges  = list(edges)
        self._title  = title or self._default_title
        if edges:
            self._x_lo = edges[0]
            self._x_hi = edges[-1]
        self._drag_start = self._drag_cur = None
        self._cb_fit_result = None
        self._cb_fit_curve  = None
        if self.auto_cb_fit:
            self._run_cb_fit()
        self._cache_pm = None
        self.update()

    # -- Crystal Ball fit --

    def _run_cb_fit(self):
        """Fit histogram with a Crystal Ball function using scipy.
        Stores fit result and curve; silently ignored if scipy is unavailable.
        """
        if not self._values or not self._edges or len(self._edges) < 3:
            return
        try:
            from scipy.optimize import curve_fit
            import numpy as _np
        except ImportError:
            return

        edges  = _np.array(self._edges)
        values = _np.array(self._values, dtype=float)
        mids   = (edges[:-1] + edges[1:]) / 2.0
        mask   = values > 0
        if mask.sum() < 5:
            return

        def crystal_ball(x, amp, mu, sigma, alpha, n):
            alpha = abs(alpha)
            n     = abs(n)
            t = (x - mu) / sigma
            result = _np.where(
                t > -alpha,
                amp * _np.exp(-0.5 * t * t),
                amp * (n / alpha) ** n * _np.exp(-0.5 * alpha * alpha)
                    / (n / alpha - alpha - t) ** n
            )
            return result

        # initial guesses
        amp0   = float(values.max())
        mu0    = float(mids[values.argmax()])
        # estimate sigma from RMS within half-max region
        half_max_mask = values > amp0 * 0.5
        sigma0 = float(mids[half_max_mask].std()) if half_max_mask.sum() > 1 else (edges[-1] - edges[0]) * 0.05
        if sigma0 <= 0:
            sigma0 = (edges[-1] - edges[0]) * 0.05

        # Determine fit window
        if self.cb_fit_range_sigma is not None:
            left_ns, right_ns = self.cb_fit_range_sigma
            fit_lo = mu0 - left_ns * sigma0
            fit_hi = mu0 + right_ns * sigma0
            fit_mask = mask & (mids >= fit_lo) & (mids <= fit_hi)
        elif self.cb_fit_range is not None:
            left_w, right_w = self.cb_fit_range
            fit_lo = mu0 - left_w
            fit_hi = mu0 + right_w
            fit_mask = mask & (mids >= fit_lo) & (mids <= fit_hi)
        else:
            fit_lo, fit_hi = float(edges[0]), float(edges[-1])
            fit_mask = mask
        if fit_mask.sum() < 5:
            fit_mask = mask  # fall back to full range
            fit_lo, fit_hi = float(edges[0]), float(edges[-1])

        try:
            popt, pcov = curve_fit(
                crystal_ball, mids[fit_mask], values[fit_mask],
                p0=[amp0, mu0, sigma0, 1.5, 3.0],
                bounds=([0, fit_lo, 0, 0.1, 1.1],
                        [_np.inf, fit_hi, (fit_hi - fit_lo), 10.0, 50.0]),
                maxfev=5000,
            )
            perr = _np.sqrt(_np.diag(pcov))
            mu_fit, mu_err     = float(popt[1]), float(perr[1])
            sigma_fit, sigma_err = abs(float(popt[2])), float(perr[2])
            self._cb_fit_result = (mu_fit, mu_err, sigma_fit, sigma_err)
            # build curve for drawing
            xs = _np.linspace(edges[0], edges[-1], 400)
            ys = crystal_ball(xs, *popt)
            self._cb_fit_curve = (xs.tolist(), ys.tolist())
        except Exception:
            pass

    def _on_view_changed(self):
        self._cache_pm = None

    def _toggle_log_x(self):
        self._log_x = self._btn_log_x.isChecked()
        if self._log_x and self._edges:
            pos_lo = next((e for e in self._edges if e > 0), None)
            if pos_lo is not None:
                self._x_lo = max(self._x_lo, pos_lo)
        self._cache_pm = None
        self.update()

    def _toggle_log_y(self):
        self._log_y = self._btn_log_y.isChecked()
        self._cache_pm = None
        self.update()

    # -- paint (QPixmap cache: bars cached, only drag overlay redrawn on mouse-move) --

    def _rebuild_cache(self):
        w, h = self.width(), self.height()
        self._plotted = False
        if w <= 0 or h <= 0:
            self._cache_pm = None
            return
        pm = QPixmap(w, h)
        pm.fill(QColor(THEME.CANVAS))
        self._cache_pm = pm
        p = QPainter(pm)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)
        self._paint_title(p)

        if not self._values or not self._edges or len(self._edges) < 2:
            self._paint_placeholder(p, "No data")
            p.end()
            return

        px, py, pw, ph = self._plot_rect()
        if pw < 10 or ph < 10:
            p.end()
            return

        x_lo, x_hi, use_log_x = self._x_view()
        to_sx = self._x_map()[0]
        vis_vals = self._visible_values(self._values, self._edges)
        vis_max = max(vis_vals, default=0.0)
        n_yticks = self.N_YTICKS

        if self._log_y:
            vis_pos = [v for v in vis_vals if v > 0]
            y_hi_log = math.log10(vis_max * 1.5) if vis_max > 0 else 0.0
            y_lo_log = (math.log10(min(vis_pos)) if vis_pos else y_hi_log - 4)
            y_lo_log = min(y_lo_log, y_hi_log - 1)
            def to_sy(v):
                if v <= 0:
                    return py + ph + 1
                lv = math.log10(v)
                return py + ph * (1.0 - (lv - y_lo_log) / (y_hi_log - y_lo_log))
            y_labels = [f"{10 ** (y_hi_log - (y_hi_log - y_lo_log) * i / n_yticks):.3g}"
                        for i in range(n_yticks + 1)]
        else:
            y_hi_lin = vis_max * 1.1 if vis_max > 0 else 1.0
            def to_sy(v): return py + ph * (1.0 - v / y_hi_lin)
            y_labels = self._lin_y_labels(y_hi_lin)

        # log x: ticks at nice exponents
        x_ticks = None
        if use_log_x:
            u_lo, u_hi = math.log10(x_lo), math.log10(x_hi)
            x_ticks = [(px + (xt - u_lo) / (u_hi - u_lo) * pw, f"10^{xt:.4g}")
                       for xt in nice_ticks(u_lo, u_hi, max(pw // 60, 2))]

        self._paint_grid(p)
        self._paint_bars(p, self._values, self._edges, QColor(self.BAR_COLOR), to_sy)
        self._paint_axes(p, y_labels, x_ticks)

        # Crystal Ball fit curve + annotation
        if self._cb_fit_curve is not None:
            xs, ys = self._cb_fit_curve
            p.setPen(QPen(QColor("#ff7b00"), 1.5))
            pts = [QPointF(to_sx(xv), to_sy(yv))
                   for xv, yv in zip(xs, ys) if x_lo <= xv <= x_hi]
            for i in range(1, len(pts)):
                p.drawLine(pts[i - 1], pts[i])

        if self._cb_fit_result is not None:
            mu, mu_err, sigma, sigma_err = self._cb_fit_result
            label1 = f"\u03bc = {mu:.2f} \u00b1 {mu_err:.2f}"
            label2 = f"\u03c3 = {sigma:.2f} \u00b1 {sigma_err:.2f}"
            p.setFont(QFont("Consolas", 9))
            p.setPen(QColor("#ff7b00"))
            p.drawText(QRectF(px + 4, py + 4, pw - 8, 16),
                       Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                       label1)
            p.drawText(QRectF(px + 4, py + 20, pw - 8, 16),
                       Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                       label2)

        p.end()
        self._plotted = True

    def paintEvent(self, event):
        if self._cache_pm is None or self._cache_pm.size() != self.size():
            self._rebuild_cache()
        if self._cache_pm is None:
            return
        p = QPainter(self)
        p.drawPixmap(0, 0, self._cache_pm)
        # drag-select overlay — not cached, drawn on top each frame
        if self._plotted:
            self._paint_drag(p)
        p.end()


# ---- 2-D histogram widget (hit_pos, h2_energy_theta, h2_moller_pos) ----

class Hist2DWidget(QWidget):
    """Simple 2-D heatmap (painter-based, no external libs)."""

    PAD_L, PAD_R, PAD_T, PAD_B = 55, 80, 28, 36
    CB_W = 18   # colorbar width

    def __init__(self, title: str = "", parent=None):
        super().__init__(parent)
        self._title = title
        self._values: Optional[List[List[float]]] = None  # [ix][iy]
        self._x_edges: List[float] = []
        self._y_edges: List[float] = []
        self._palette_idx = 0   # index into PALETTE_NAMES
        self._log_z = False
        self.setMinimumSize(200, 200)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self._btn_log_z = _log_toggle(self, "logZ", self._toggle_log_z)
        self._cache_pm: Optional[QPixmap] = None

    def resizeEvent(self, event):
        super().resizeEvent(event)
        w = self.width()
        self._btn_log_z.move(w - self.PAD_R + 10, 4)
        self._cache_pm = None

    def set_data(self, values: List[List[float]],
                 x_edges: List[float], y_edges: List[float]):
        self._values = values
        self._x_edges = x_edges
        self._y_edges = y_edges
        self._cache_pm = None
        self.update()

    def contextMenuEvent(self, event):
        menu = QMenu(self)
        menu.setStyleSheet(themed(MENU_QSS))
        pal_menu = menu.addMenu("Palette")
        for i, name in enumerate(PALETTE_NAMES):
            act = pal_menu.addAction(name)
            act.setCheckable(True)
            act.setChecked(i == self._palette_idx)
            idx = i
            act.triggered.connect(lambda _c=False, ii=idx: self._set_palette(ii))
        menu.exec(event.globalPos())

    def _toggle_log_z(self):
        self._log_z = self._btn_log_z.isChecked()
        self._cache_pm = None
        self.update()

    def _set_palette(self, idx: int):
        self._palette_idx = idx
        self._cache_pm = None
        self.update()

    def _rebuild_cache(self):
        """Render entire 2-D heatmap to a QPixmap (cached until data/palette/size changes)."""
        w, h = self.width(), self.height()
        if w <= 0 or h <= 0:
            self._cache_pm = None
            return
        pm = QPixmap(w, h)
        pm.fill(QColor(THEME.CANVAS))
        p = QPainter(pm)
        p.setRenderHint(QPainter.RenderHint.Antialiasing)

        if self._title:
            p.setPen(QColor(THEME.ACCENT))
            p.setFont(QFont("Consolas", 10, QFont.Weight.Bold))
            p.drawText(QRectF(self.PAD_L, 4, w - self.PAD_L - self.PAD_R, 20),
                       Qt.AlignmentFlag.AlignLeft | Qt.AlignmentFlag.AlignVCenter,
                       self._title)

        if self._values is None:
            p.setPen(QColor(THEME.TEXT_MUTED))
            p.setFont(QFont("Consolas", 10))
            p.drawText(QRectF(0, 0, w, h), Qt.AlignmentFlag.AlignCenter, "No data")
            p.end()
            self._cache_pm = pm
            return

        vals = self._values
        xe = self._x_edges
        ye = self._y_edges
        if len(xe) < 2 or len(ye) < 2:
            p.end()
            self._cache_pm = pm
            return

        px = self.PAD_L
        py = self.PAD_T
        cb_x = w - self.PAD_R + 8
        pw = w - self.PAD_L - self.PAD_R
        ph = h - self.PAD_T - self.PAD_B

        if pw < 10 or ph < 10:
            p.end()
            self._cache_pm = pm
            return

        x_lo, x_hi = xe[0], xe[-1]
        y_lo, y_hi = ye[0], ye[-1]

        def to_sx(v): return px + (v - x_lo) / (x_hi - x_lo) * pw
        def to_sy(v): return py + ph - (v - y_lo) / (y_hi - y_lo) * ph

        flat = [v for row in vals for v in row if math.isfinite(v) and v > 0]
        vmax = max(flat) if flat else 1.0
        use_log_z = self._log_z and vmax > 0
        log_vmax = math.log10(vmax) if use_log_z else vmax
        log_vmin = math.log10(min(flat)) if (use_log_z and flat) else 0.0
        if use_log_z and log_vmin >= log_vmax:
            log_vmin = log_vmax - 1

        palette = PALETTES.get(PALETTE_NAMES[self._palette_idx],
                               PALETTES[PALETTE_NAMES[0]])

        # Heatmap and colour bar as QImages; on failure the plot area stays
        # empty and only the axes are drawn.
        try:
            import numpy as _np
            nx, ny = len(xe) - 1, len(ye) - 1
            vals_np = _np.array(vals, dtype=_np.float64)          # (nx, ny)
            if use_log_z:
                with _np.errstate(divide='ignore', invalid='ignore'):
                    lv = _np.where(vals_np > 0,
                                   _np.log10(_np.maximum(vals_np, 1e-300)),
                                   _np.nan)
                span = log_vmax - log_vmin
                t_bins = ((lv - log_vmin) / span) if span > 0 else _np.zeros_like(lv)
            else:
                t_bins = (vals_np / vmax) if vmax > 0 else _np.zeros_like(vals_np)
            t_bins = _np.clip(t_bins, 0.0, 1.0)
            mask   = (vals_np > 0) & _np.isfinite(t_bins)         # (nx, ny)

            # map each screen pixel to its data bin
            xi = _np.clip(
                (_np.arange(pw, dtype=_np.float64) / pw * nx).astype(_np.int32),
                0, nx - 1)                                         # (pw,)
            yi = _np.clip(
                ((1.0 - _np.arange(ph, dtype=_np.float64) / ph) * ny
                 ).astype(_np.int32),
                0, ny - 1)                                         # (ph,)

            t_img   = t_bins[xi[_np.newaxis, :], yi[:, _np.newaxis]]  # (ph, pw)
            msk_img =  mask[xi[_np.newaxis, :], yi[:, _np.newaxis]]

            rgb = cmap_rgb_array(t_img, palette).astype(_np.uint32)
            argb = (_np.where(msk_img, _np.uint32(0xFF000000), _np.uint32(0))
                    | rgb[..., 0] << 16 | rgb[..., 1] << 8 | rgb[..., 2])
            img = QImage(argb.tobytes(), pw, ph, pw * 4,
                         QImage.Format.Format_ARGB32)
            p.drawImage(px, py, img)

            # colorbar as 1×ph QImage scaled to CB_W wide
            t_cb = 1.0 - _np.arange(ph, dtype=_np.float64) / max(ph - 1, 1)
            rgb_cb = cmap_rgb_array(t_cb, palette).astype(_np.uint32)
            argb_cb = (_np.uint32(0xFF000000) | rgb_cb[:, 0] << 16
                       | rgb_cb[:, 1] << 8 | rgb_cb[:, 2])
            cb_img = QImage(argb_cb.tobytes(), 1, ph, 4,
                            QImage.Format.Format_ARGB32)
            p.drawImage(QRectF(cb_x, py, self.CB_W, ph), cb_img,
                        QRectF(0, 0, 1, ph))
        except Exception:
            pass

        # axes + tick labels + colorbar border (always)
        p.setPen(QPen(QColor(THEME.BORDER), 1))
        p.drawLine(QPointF(px, py), QPointF(px, py + ph))
        p.drawLine(QPointF(px, py + ph), QPointF(px + pw, py + ph))
        p.drawRect(QRectF(cb_x, py, self.CB_W, ph))
        p.setPen(QColor(THEME.TEXT_DIM))
        p.setFont(QFont("Consolas", 8))
        for xt in nice_ticks(x_lo, x_hi, max(pw // 60, 2)):
            sx = to_sx(xt)
            p.drawText(QRectF(sx - 28, py + ph + 2, 56, 16),
                       Qt.AlignmentFlag.AlignCenter, f"{xt:.4g}")
        for yt in nice_ticks(y_lo, y_hi, max(ph // 40, 2)):
            sy = to_sy(yt)
            p.drawText(QRectF(0, sy - 8, self.PAD_L - 4, 16),
                       Qt.AlignmentFlag.AlignRight | Qt.AlignmentFlag.AlignVCenter,
                       f"{yt:.4g}")
        if use_log_z:
            p.drawText(QRectF(cb_x, py + ph + 2, self.CB_W + 40, 14),
                       Qt.AlignmentFlag.AlignLeft, f"10^{log_vmin:.2g}")
            p.drawText(QRectF(cb_x, py - 12, self.CB_W + 40, 14),
                       Qt.AlignmentFlag.AlignLeft, f"10^{log_vmax:.2g}")
        else:
            p.drawText(QRectF(cb_x, py + ph + 2, self.CB_W + 20, 14),
                       Qt.AlignmentFlag.AlignLeft, "0")
            p.drawText(QRectF(cb_x, py - 12, self.CB_W + 40, 14),
                       Qt.AlignmentFlag.AlignLeft, f"{vmax:.3g}")
        p.end()
        self._cache_pm = pm

    def paintEvent(self, event):
        if self._cache_pm is None or self._cache_pm.size() != self.size():
            self._rebuild_cache()
        if self._cache_pm is None:
            return
        p = QPainter(self)
        p.drawPixmap(0, 0, self._cache_pm)
        p.end()


# ---- HyCal replay map widget ----

class HyCalReplayMapWidget(HyCalMapWidget):
    """HyCal map of a per-module quantity from quick_check module_energy
    (hit counts or mean energies)."""

    def __init__(self, parent=None):
        super().__init__(parent, shrink=0.92, margin_top=8,
                         enable_zoom_pan=True, include_lms=False,
                         toggle_select=True)

    def set_module_counts(self, counts: Dict[str, float]):
        if not counts:
            self._values = {}
            self._vmin = 0.0
            self._vmax = 1.0
            self.update()
            return
        # Seed all PbWO4 modules with 0 so the full crystal region is
        # visible (min-palette colour) even when the ROOT file lacks an
        # entry for some modules.
        base: Dict[str, float] = {
            m.name: 0.0 for m in self._modules if m.mod_type == "PbWO4"
        }
        base.update(counts)
        self._values = base
        self._vmin = 0.0
        self._vmax = max(counts.values())
        self.update()

    def _fmt_value(self, v: float) -> str:
        return f"{v:.0f}"

    def _tooltip_text(self, name: str) -> str:
        v = self._values.get(name)
        if v is None:
            return name
        return f"{name}: {v:.0f} hits"

    def _paint_empty(self, p, w, h):
        if not self._values:
            p.setPen(QColor(THEME.TEXT_MUTED))
            p.setFont(QFont("Consolas", 12))
            p.drawText(QRectF(0, 0, w, h),
                       Qt.AlignmentFlag.AlignCenter, "Load a ROOT file to view")


# ---- Control Panel (left side) ----

_BTN_PRIMARY = themed(
    "QPushButton{background:#1f6feb;color:white;border:1px solid #388bfd;"
    "padding:5px 16px;font:bold 11pt Consolas;border-radius:3px;}"
    "QPushButton:hover{background:#388bfd;}"
    "QPushButton:disabled{background:#21262d;color:#555;border-color:#30363d;}")

_BTN_NORMAL = themed(
    "QPushButton{background:#21262d;color:#c9d1d9;border:1px solid #30363d;"
    "padding:5px 14px;font:bold 11pt Consolas;border-radius:3px;}"
    "QPushButton:hover{background:#30363d;}"
    "QPushButton:disabled{color:#555;}")

_BTN_DANGER = themed(
    "QPushButton{background:#3d1f22;color:#f85149;border:1px solid #f85149;"
    "padding:5px 14px;font:bold 11pt Consolas;border-radius:3px;}"
    "QPushButton:hover{background:#5a2329;}"
    "QPushButton:disabled{color:#555;}")

_LINEEDIT_SS = themed(
    "QLineEdit{background:#161b22;color:#c9d1d9;"
    "border:1px solid #30363d;border-radius:3px;padding:3px 6px;"
    "font-family:Consolas;font-size:11pt;}")

_SPINBOX_SS = themed(
    "QSpinBox{background:#161b22;color:#c9d1d9;"
    "border:1px solid #30363d;border-radius:3px;padding:3px 6px;"
    "font-family:Consolas;font-size:11pt;}")

_COMBO_SS = themed(
    "QComboBox{background:#161b22;color:#c9d1d9;"
    "border:1px solid #30363d;border-radius:3px;padding:3px 6px;"
    "font-family:Consolas;font-size:11pt;}"
    "QComboBox::drop-down{border:none;width:18px;}"
    "QComboBox::down-arrow{border-left:4px solid transparent;"
    "border-right:4px solid transparent;border-top:5px solid #8b949e;"
    "margin-right:4px;}"
    "QComboBox QAbstractItemView{background:#161b22;color:#c9d1d9;"
    "border:1px solid #30363d;selection-background-color:#1f6feb;}")

_CHK_SS = themed(
    "QCheckBox{color:#c9d1d9;font-family:Consolas;font-size:11pt;spacing:6px;}"
    "QCheckBox::indicator{width:15px;height:15px;"
    "border:1px solid #30363d;border-radius:2px;background:#161b22;}"
    "QCheckBox::indicator:checked{background:#1f6feb;border-color:#388bfd;}")

_LBL_SS  = themed("QLabel{color:#c9d1d9;font-family:Consolas;font-size:11pt;}")
_LBL_MUT = themed("QLabel{color:#8b949e;font-family:Consolas;font-size:10pt;}")


def _le(w: QLineEdit) -> QLineEdit:
    w.setFont(QFont("Consolas", 10))
    w.setStyleSheet(_LINEEDIT_SS)
    return w


def _sp(w: QSpinBox) -> QSpinBox:
    w.setFont(QFont("Consolas", 10))
    w.setStyleSheet(_SPINBOX_SS)
    return w


def _form_dialog(parent: QWidget, title: str, min_width: int = 480
                 ) -> Tuple[QDialog, QFormLayout]:
    """Settings dialog with a right-aligned form layout."""
    dlg = QDialog(parent)
    dlg.setWindowTitle(title)
    dlg.setMinimumWidth(min_width)
    dlg.setStyleSheet(themed("QDialog{background:#0d1117;}"))
    form = QFormLayout(dlg)
    form.setSpacing(8)
    form.setLabelAlignment(Qt.AlignmentFlag.AlignRight)
    return dlg, form


def _ok_cancel_row(dlg: QDialog, form: QFormLayout, on_accept,
                   ok_text: str = "OK") -> None:
    btns = QHBoxLayout()
    ok_btn = QPushButton(ok_text)
    ok_btn.setStyleSheet(_BTN_PRIMARY)
    ok_btn.clicked.connect(on_accept)
    ca_btn = QPushButton("Cancel")
    ca_btn.setStyleSheet(_BTN_NORMAL)
    ca_btn.clicked.connect(dlg.reject)
    btns.addStretch()
    btns.addWidget(ok_btn)
    btns.addWidget(ca_btn)
    form.addRow(btns)


def _browse_row(dlg: QDialog, text: str, title: str,
                file_filter: Optional[str] = None) -> Tuple[QLineEdit, QWidget]:
    """Line edit plus a Browse… button picking a directory, or a file
    matching ``file_filter`` when one is given."""
    edit = _le(QLineEdit(text))
    btn = QPushButton("Browse…")
    btn.setFixedWidth(90)
    btn.setStyleSheet(_BTN_NORMAL)

    def browse():
        if file_filter is None:
            path = QFileDialog.getExistingDirectory(dlg, title, edit.text())
        else:
            path = QFileDialog.getOpenFileName(dlg, title, edit.text(),
                                               file_filter)[0]
        if path:
            edit.setText(path)

    btn.clicked.connect(browse)
    row = QWidget()
    lay = QHBoxLayout(row)
    lay.setContentsMargins(0, 0, 0, 0)
    lay.setSpacing(4)
    lay.addWidget(edit)
    lay.addWidget(btn)
    return edit, row


class ControlPanel(QWidget):
    """Left panel: SCP, replay, quick_check controls + log output."""

    # Emitted when a quick_check ROOT file has been generated (or opened).
    rootFileReady     = pyqtSignal(str)   # path to ROOT file
    # Emitted when a replay-filter JSON report has been generated.
    filterReportReady = pyqtSignal(str)   # path to *.report.json

    def __init__(self, parent=None):
        super().__init__(parent)
        self._process: Optional[QProcess] = None
        self._pending_steps: List[str] = []   # ["scp", "replay", "hadd", "filter", "qcheck"]
        self._current_step: str = ""          # step currently running
        self._evio_dir: str = ""              # set after SCP completes
        self._recon_dir: str = ""             # set after replay completes
        self._hadd_out: str = ""              # merged ROOT file from hadd
        self._hadd_inputs: List[str] = []     # individual files to delete after hadd
        self._filter_out: str = ""            # filtered ROOT file from replay filter
        self._filter_report: str = ""         # JSON report written by replay filter
        self._qcheck_out: str = ""            # final ROOT file path
        self._build_ui()

    # -- UI construction --

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(6, 6, 6, 6)
        root.setSpacing(6)

        # ---- Do It All / Stop row (always visible at top) ----
        do_all_row = QHBoxLayout()
        self._do_all_btn = QPushButton("Do It All  (SCP → Replay → hadd → Filter → Quick Check)")
        self._do_all_btn.setStyleSheet(_BTN_PRIMARY)
        self._do_all_btn.clicked.connect(self._on_do_all_clicked)
        self._stop_btn = QPushButton("Stop")
        self._stop_btn.setStyleSheet(_BTN_DANGER)
        self._stop_btn.setEnabled(False)
        self._stop_btn.clicked.connect(self._on_stop)
        self._auto_delete_evio_chk = QCheckBox("Auto-delete EVIO files")
        self._auto_delete_evio_chk.setStyleSheet(_CHK_SS)
        do_all_row.addWidget(self._do_all_btn)
        do_all_row.addWidget(self._stop_btn)
        do_all_row.addSpacing(12)
        do_all_row.addWidget(self._auto_delete_evio_chk)
        do_all_row.addStretch()
        root.addLayout(do_all_row)

        # ---- Step buttons (compact rows) ----
        steps_widget = QWidget()
        steps_lay = QVBoxLayout(steps_widget)
        steps_lay.setContentsMargins(4, 4, 4, 4)
        steps_lay.setSpacing(6)

        # 1. SCP row
        scp_row = QHBoxLayout()
        self._scp_cfg_btn = QPushButton("⚙  SCP Settings…")
        self._scp_cfg_btn.setStyleSheet(_BTN_NORMAL)
        self._scp_cfg_btn.clicked.connect(self._open_scp_dialog)
        self._scp_btn = QPushButton("1. Get Data")
        self._scp_btn.setStyleSheet(_BTN_PRIMARY)
        self._scp_btn.clicked.connect(lambda: self._start_pipeline(["scp"]))
        self._check_disk_btn = QPushButton("Check Disk")
        self._check_disk_btn.setStyleSheet(_BTN_NORMAL)
        self._check_disk_btn.clicked.connect(self._on_check_disk)
        self._disk_lbl = QLabel("")
        self._disk_lbl.setStyleSheet(_LBL_MUT)
        scp_row.addWidget(self._scp_cfg_btn)
        scp_row.addWidget(self._scp_btn)
        scp_row.addWidget(self._check_disk_btn)
        scp_row.addWidget(self._disk_lbl)
        scp_row.addStretch()
        steps_lay.addLayout(scp_row)

        # 2. Replay row
        rep_row = QHBoxLayout()
        self._rep_cfg_btn = QPushButton("⚙  Replay Settings…")
        self._rep_cfg_btn.setStyleSheet(_BTN_NORMAL)
        self._rep_cfg_btn.clicked.connect(self._open_replay_dialog)
        self._replay_btn = QPushButton("2. Run Replay")
        self._replay_btn.setStyleSheet(_BTN_PRIMARY)
        self._replay_btn.clicked.connect(lambda: self._start_pipeline(["replay"]))
        rep_row.addWidget(self._rep_cfg_btn)
        rep_row.addWidget(self._replay_btn)
        rep_row.addStretch()
        steps_lay.addLayout(rep_row)

        # 3. Replay Filter row
        flt_row = QHBoxLayout()
        self._flt_cfg_btn = QPushButton("⚙  Filter Settings…")
        self._flt_cfg_btn.setStyleSheet(_BTN_NORMAL)
        self._flt_cfg_btn.clicked.connect(self._open_filter_dialog)
        self._filter_btn = QPushButton("3. Replay Filter")
        self._filter_btn.setStyleSheet(_BTN_PRIMARY)
        self._filter_btn.clicked.connect(lambda: self._start_pipeline(["filter"]))
        flt_row.addWidget(self._flt_cfg_btn)
        flt_row.addWidget(self._filter_btn)
        flt_row.addStretch()
        steps_lay.addLayout(flt_row)

        # 4. Quick Check row
        qc_row = QHBoxLayout()
        self._qc_cfg_btn = QPushButton("⚙  QCheck Settings…")
        self._qc_cfg_btn.setStyleSheet(_BTN_NORMAL)
        self._qc_cfg_btn.clicked.connect(self._open_qcheck_dialog)
        self._qcheck_btn = QPushButton("4. Quick Check")
        self._qcheck_btn.setStyleSheet(_BTN_PRIMARY)
        self._qcheck_btn.clicked.connect(lambda: self._start_pipeline(["qcheck"]))
        self._open_root_btn = QPushButton("Open ROOT…")
        self._open_root_btn.setStyleSheet(_BTN_NORMAL)
        self._open_root_btn.clicked.connect(self._browse_root_file)
        qc_row.addWidget(self._qc_cfg_btn)
        qc_row.addWidget(self._qcheck_btn)
        qc_row.addWidget(self._open_root_btn)
        qc_row.addStretch()
        steps_lay.addLayout(qc_row)

        root.addWidget(steps_widget)

        # ---- Hidden parameter widgets (not shown, values read by pipeline) ----
        _hidden = QWidget(); _hidden.hide()
        _hlay = QVBoxLayout(_hidden)

        # SCP params
        self._run_edit = QLineEdit()
        self._host_edit = QLineEdit(REMOTE_HOST)
        self._remote_base_edit = QLineEdit(REMOTE_DATA_BASE)
        self._local_base_edit = QLineEdit(LOCAL_DATA_BASE)
        self._f_start = QSpinBox(); self._f_start.setRange(0, 9999); self._f_start.setValue(0)
        self._f_end   = QSpinBox(); self._f_end.setRange(0, 9999);   self._f_end.setValue(99)
        for w in (self._run_edit, self._host_edit, self._remote_base_edit,
                  self._local_base_edit, self._f_start, self._f_end):
            _hlay.addWidget(w)

        # Replay params
        self._evio_edit   = QLineEdit()
        self._outdir_edit = QLineEdit()
        self._threads_spin = QSpinBox(); self._threads_spin.setRange(1, 256); self._threads_spin.setValue(50)
        self._max_events_rep = QLineEdit("-1")
        self._max_files_spin = QSpinBox(); self._max_files_spin.setRange(-1, 9999); self._max_files_spin.setValue(-1)
        self._max_files_spin.setSpecialValueText("all")
        self._daq_config_edit = QLineEdit()
        self._hycal_map_edit  = QLineEdit()
        self._gem_ped_edit    = QLineEdit()
        self._zerosup_edit    = QLineEdit("5")
        self._prad1_chk = QCheckBox()
        for w in (self._evio_edit, self._outdir_edit, self._threads_spin,
                  self._max_events_rep, self._max_files_spin, self._daq_config_edit,
                  self._hycal_map_edit, self._gem_ped_edit, self._zerosup_edit,
                  self._prad1_chk):
            _hlay.addWidget(w)

        # Replay Filter params
        self._filter_input_edit  = QLineEdit()
        self._filter_output_edit = QLineEdit()
        self._filter_cut_edit    = QLineEdit()   # -c cut JSON for replay filter
        self._max_events_flt     = QLineEdit("-1")
        for w in (self._filter_input_edit, self._filter_output_edit,
                  self._filter_cut_edit, self._max_events_flt):
            _hlay.addWidget(w)

        # Quick Check params
        self._qc_input_edit  = QLineEdit()
        self._qc_output_edit = QLineEdit()
        self._max_events_qc  = QLineEdit("-1")
        for w in (self._qc_input_edit, self._qc_output_edit, self._max_events_qc):
            _hlay.addWidget(w)

        root.addWidget(_hidden)

        # ---- Log output ----
        log_header = QHBoxLayout()
        log_lbl = QLabel("Log")
        log_lbl.setStyleSheet(themed(
            "QLabel{color:#58a6ff;font:bold 10pt Consolas;}"))
        self._status_lbl = QLabel("Ready")
        self._status_lbl.setStyleSheet(_LBL_MUT)
        clr_btn = QPushButton("Clear")
        clr_btn.setFixedWidth(80)
        clr_btn.setStyleSheet(_BTN_NORMAL)
        clr_btn.clicked.connect(lambda: self._console.clear())
        log_header.addWidget(log_lbl)
        log_header.addWidget(self._status_lbl)
        log_header.addStretch()
        log_header.addWidget(clr_btn)
        root.addLayout(log_header)

        self._console = QTextEdit()
        self._console.setReadOnly(True)
        self._console.setFont(QFont("Monospace", 10))
        self._console.document().setMaximumBlockCount(10000)
        self._console.setStyleSheet(themed(
            "QTextEdit{background:#0a0e14;color:#c9d1d9;"
            "border:1px solid #30363d;font-family:Monospace;font-size:10pt;}"))
        root.addWidget(self._console, stretch=1)

    # -- Settings dialogs --

    def _on_do_all_clicked(self):
        """Show a quick-config popup, then launch the full pipeline."""
        dlg, form = _form_dialog(self, "Do It All — Quick Setup")

        run_e = _le(QLineEdit(self._run_edit.text()))
        run_e.setPlaceholderText("e.g. 024388")

        thr_sp = _sp(QSpinBox())
        thr_sp.setRange(1, 256)
        thr_sp.setValue(25)

        cut_e, cut_row = _browse_row(dlg, self._filter_cut_edit.text(),
                                     "Select cut JSON",
                                     "JSON files (*.json);;All files (*)")
        cut_e.setPlaceholderText("(optional) path/to/cut.json")

        zsup_e = _le(QLineEdit(self._zerosup_edit.text()))
        zsup_e.setPlaceholderText("e.g. 5")

        fs_sp = _sp(QSpinBox())
        fs_sp.setRange(0, 9999)
        fs_sp.setValue(self._f_start.value())
        fe_sp = _sp(QSpinBox())
        fe_sp.setRange(0, 9999)
        fe_sp.setValue(self._f_end.value())
        frange_row = QWidget()
        fr = QHBoxLayout(frange_row)
        fr.setContentsMargins(0, 0, 0, 0)
        fr.setSpacing(4)
        fr.addWidget(fs_sp)
        fr.addWidget(QLabel("—"))
        fr.addWidget(fe_sp)
        fr.addStretch()

        form.addRow("Run number:", run_e)
        form.addRow("File index range:", frange_row)
        form.addRow("Threads (-j):", thr_sp)
        form.addRow("Filter cut JSON (-c):", cut_row)
        form.addRow("GEM zero-sup (-z):", zsup_e)

        def accept():
            rn = run_e.text().strip()
            if not rn:
                run_e.setStyleSheet(_LINEEDIT_SS + "border:1px solid #f85149;")
                return
            self._run_edit.setText(rn)
            self._f_start.setValue(fs_sp.value())
            self._f_end.setValue(fe_sp.value())
            self._threads_spin.setValue(thr_sp.value())
            self._filter_cut_edit.setText(cut_e.text())
            self._zerosup_edit.setText(zsup_e.text())
            dlg.accept()

        _ok_cancel_row(dlg, form, accept, "Start Pipeline")

        if dlg.exec() == QDialog.DialogCode.Accepted:
            self._start_pipeline(["scp", "replay", "hadd", "filter", "qcheck"])

    def _open_scp_dialog(self):
        dlg, form = _form_dialog(self, "SCP Settings")

        run_e   = _le(QLineEdit(self._run_edit.text()))
        run_e.setPlaceholderText("e.g. 024100")
        host_e  = _le(QLineEdit(self._host_edit.text()))
        rbas_e  = _le(QLineEdit(self._remote_base_edit.text()))
        lbas_e, lbas_row = _browse_row(dlg, self._local_base_edit.text(),
                                       "Local Base Directory")

        fs_sp   = _sp(QSpinBox()); fs_sp.setRange(0, 9999); fs_sp.setValue(self._f_start.value())
        fe_sp   = _sp(QSpinBox()); fe_sp.setRange(0, 9999); fe_sp.setValue(self._f_end.value())
        frange_row = QHBoxLayout()
        frange_row.addWidget(fs_sp); frange_row.addWidget(QLabel("—")); frange_row.addWidget(fe_sp); frange_row.addStretch()

        disk_lbl = QLabel("(not checked)")
        disk_lbl.setStyleSheet(_LBL_MUT)

        def check_disk():
            dirs = _run_dirs(run_e.text(), host_e.text(), rbas_e.text(), lbas_e.text())
            if dirs is None:
                disk_lbl.setText("<span style='color:#f85149'>Enter run number</span>")
                return
            disk_lbl.setText("Checking…")
            disk_lbl.setText(_disk_check_html(dirs, fs_sp.value(), fe_sp.value()))

        form.addRow("Run number:", run_e)
        form.addRow("Remote host:", host_e)
        form.addRow("Remote base dir:", rbas_e)
        form.addRow("Local base dir:", lbas_row)
        form.addRow("File index range:", frange_row)
        chk_btn = QPushButton("Check Disk"); chk_btn.setStyleSheet(_BTN_NORMAL); chk_btn.clicked.connect(check_disk)
        chk_row = QHBoxLayout(); chk_row.addWidget(chk_btn); chk_row.addWidget(disk_lbl); chk_row.addStretch()
        form.addRow("", chk_row)

        def accept():
            self._run_edit.setText(run_e.text())
            self._host_edit.setText(host_e.text())
            self._remote_base_edit.setText(rbas_e.text())
            self._local_base_edit.setText(lbas_e.text())
            self._f_start.setValue(fs_sp.value())
            self._f_end.setValue(fe_sp.value())
            dlg.accept()

        _ok_cancel_row(dlg, form, accept)
        dlg.exec()

    def _open_replay_dialog(self):
        dlg, form = _form_dialog(self, "Replay Settings", 520)

        evio_e, evio_row   = _browse_row(dlg, self._evio_edit.text(), "EVIO Directory")
        out_e,  out_row    = _browse_row(dlg, self._outdir_edit.text(), "Replay Output Directory")
        thr_sp = _sp(QSpinBox()); thr_sp.setRange(1, 256); thr_sp.setValue(self._threads_spin.value())
        nev_e  = _le(QLineEdit(self._max_events_rep.text())); nev_e.setToolTip("-1 = no limit")
        nf_sp  = _sp(QSpinBox()); nf_sp.setRange(-1, 9999); nf_sp.setValue(self._max_files_spin.value()); nf_sp.setSpecialValueText("all")
        daq_e  = _le(QLineEdit(self._daq_config_edit.text())); daq_e.setPlaceholderText("(default)")
        hmap_e = _le(QLineEdit(self._hycal_map_edit.text())); hmap_e.setPlaceholderText("(default)")
        gem_e  = _le(QLineEdit(self._gem_ped_edit.text())); gem_e.setPlaceholderText("(none)")
        zsup_e = _le(QLineEdit(self._zerosup_edit.text()))
        p1_chk = QCheckBox("PRad-1 mode (-p)"); p1_chk.setStyleSheet(_CHK_SS); p1_chk.setChecked(self._prad1_chk.isChecked())

        form.addRow("EVIO dir / file:", evio_row)
        form.addRow("Output dir:", out_row)
        form.addRow("Threads (-j):", thr_sp)
        form.addRow("Max events (-n):", nev_e)
        form.addRow("Max files (-f):", nf_sp)
        form.addRow("DAQ config (-c):", daq_e)
        form.addRow("HyCal map (-d):", hmap_e)
        form.addRow("GEM pedestal (-g):", gem_e)
        form.addRow("Zero-sup thresh (-z):", zsup_e)
        form.addRow("", p1_chk)

        def accept():
            self._evio_edit.setText(evio_e.text())
            self._outdir_edit.setText(out_e.text())
            self._threads_spin.setValue(thr_sp.value())
            self._max_events_rep.setText(nev_e.text())
            self._max_files_spin.setValue(nf_sp.value())
            self._daq_config_edit.setText(daq_e.text())
            self._hycal_map_edit.setText(hmap_e.text())
            self._gem_ped_edit.setText(gem_e.text())
            self._zerosup_edit.setText(zsup_e.text())
            self._prad1_chk.setChecked(p1_chk.isChecked())
            dlg.accept()

        _ok_cancel_row(dlg, form, accept)
        dlg.exec()

    def _open_filter_dialog(self):
        dlg, form = _form_dialog(self, "Replay Filter Settings")

        inp_e, inp_row = _browse_row(dlg, self._filter_input_edit.text(),
                                     "Filter Input ROOT File",
                                     "ROOT files (*.root);;All files (*)")
        out_e = _le(QLineEdit(self._filter_output_edit.text())); out_e.setPlaceholderText("(auto: prad_XXXXXX_filter.root)")
        nev_e = _le(QLineEdit(self._max_events_flt.text())); nev_e.setToolTip("-1 = no limit")

        form.addRow("Input ROOT (-i):", inp_row)
        form.addRow("Output ROOT (-o):", out_e)
        form.addRow("Max events (-n):", nev_e)

        def accept():
            self._filter_input_edit.setText(inp_e.text())
            self._filter_output_edit.setText(out_e.text())
            self._max_events_flt.setText(nev_e.text())
            dlg.accept()

        _ok_cancel_row(dlg, form, accept)
        dlg.exec()

    def _open_qcheck_dialog(self):
        dlg, form = _form_dialog(self, "Quick Check Settings")

        inp_e, inp_row = _browse_row(dlg, self._qc_input_edit.text(), "Quick Check Input")
        out_e = _le(QLineEdit(self._qc_output_edit.text())); out_e.setPlaceholderText("(auto)")
        nev_e = _le(QLineEdit(self._max_events_qc.text())); nev_e.setToolTip("-1 = no limit")

        form.addRow("Input (dir / file):", inp_row)
        form.addRow("Output ROOT (-o):", out_e)
        form.addRow("Max events (-n):", nev_e)

        def accept():
            self._qc_input_edit.setText(inp_e.text())
            self._qc_output_edit.setText(out_e.text())
            self._max_events_qc.setText(nev_e.text())
            dlg.accept()

        _ok_cancel_row(dlg, form, accept)
        dlg.exec()

    def _browse_root_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open ROOT File", "", "ROOT files (*.root);;All files (*)")
        if path:
            self._log(f"<span style='color:#8b949e'>Opening {path}</span>")
            self.rootFileReady.emit(path)

    # -- Disk space check --

    def _scp_dirs(self):
        """_run_dirs() of the panel's SCP fields."""
        return _run_dirs(self._run_edit.text(), self._host_edit.text(),
                         self._remote_base_edit.text(),
                         self._local_base_edit.text())

    def _on_check_disk(self):
        dirs = self._scp_dirs()
        if dirs is None:
            self._log("<span style='color:#f85149'>Enter a run number first.</span>")
            return
        self._disk_lbl.setText("Checking…")
        self._disk_lbl.setText(
            _disk_check_html(dirs, self._f_start.value(), self._f_end.value()))

    # -- Pipeline logic --

    def _start_pipeline(self, steps: List[str]):
        if self._process is not None and \
                self._process.state() != QProcess.ProcessState.NotRunning:
            self._log("<span style='color:#f85149'>A process is already running.</span>")
            return
        self._pending_steps = list(steps)
        self._set_running(True)
        self._run_next_step()

    def _run_next_step(self):
        if not self._pending_steps:
            self._set_running(False)
            self._log("<span style='color:#3fb950'>[All steps complete]</span>")
            return
        step = self._pending_steps.pop(0)
        self._current_step = step
        if step == "scp":
            self._run_scp()
        elif step == "replay":
            self._run_replay()
        elif step == "hadd":
            self._run_hadd()
        elif step == "filter":
            self._run_filter()
        elif step == "qcheck":
            self._run_qcheck()
        else:
            self._run_next_step()

    # -- SCP step --

    def _run_scp(self):
        dirs = self._scp_dirs()
        if dirs is None:
            self._log("<span style='color:#f85149'>No run number specified.</span>")
            self._set_running(False)
            return
        run_tag, host, remote_run_dir, local_base, local_run_dir = dirs
        f_start = self._f_start.value()
        f_end = self._f_end.value()

        # -- pre-check: find files in range that already exist locally --
        existing = local_evio_in_range(local_run_dir, f_start, f_end,
                                       f"{run_tag}.evio.*")
        if existing:
            self._log(
                f"<span style='color:#d29922'>{len(existing)} file(s) already present "
                f"in {local_run_dir} — will be skipped.</span>")

        # -- disk space check --
        self._log("<span style='color:#8b949e'>Checking disk space…</span>")
        self._status_lbl.setText("Checking disk space…")
        try:
            needed, free = check_disk_space(host, remote_run_dir, local_base,
                                            f_start, f_end, local_run_dir)
            self._log(_disk_summary_html(needed, free, "Disk: "))
            self._disk_lbl.setText(_disk_summary_html(needed, free))
            if free < needed:
                self._log("<span style='color:#f85149'>Insufficient disk space — stopping.</span>")
                self._pending_steps.clear()
                self._set_running(False)
                return
        except RuntimeError as exc:
            self._log(f"<span style='color:#f0883e'>Disk check SSH error: {exc} — continuing anyway.</span>")
        except Exception as exc:
            self._log(f"<span style='color:#f0883e'>Disk check error: {exc} — continuing anyway.</span>")

        self._evio_dir = local_run_dir

        self._log(
            f"<span style='color:#8b949e'>Run {self._run_edit.text().strip()}, "
            f"files {f_start}–{f_end} → {local_run_dir}</span>")
        self._status_lbl.setText("Getting data…")
        self._launch_process(["bash", "-c", scp_bash(
            host, remote_run_dir, local_run_dir, f_start, f_end)])

    # -- Replay step --

    def _run_replay(self):
        evio_path = self._evio_edit.text().strip() or self._evio_dir
        if not evio_path:
            self._log("<span style='color:#f85149'>No EVIO path specified for replay.</span>")
            self._set_running(False)
            return
        out_dir = self._outdir_edit.text().strip()
        if not out_dir:
            out_dir = os.path.join(
                _RECON_BASE, _run_tag(self._run_edit.text())
                or os.path.basename(evio_path.rstrip("/")) + "_recon")
        os.makedirs(out_dir, exist_ok=True)
        self._recon_dir = out_dir

        cmd = [_REPLAY_RECON_CMD, evio_path]
        cmd += ["-o", out_dir]
        cmd += ["-j", str(self._threads_spin.value())]
        _opt(cmd, "-n", self._max_events_rep.text(), skip="-1")
        n_f = self._max_files_spin.value()
        if n_f > 0:
            cmd += ["-f", str(n_f)]
        _opt(cmd, "-c", self._daq_config_edit.text())
        _opt(cmd, "-d", self._hycal_map_edit.text())
        _opt(cmd, "-g", self._gem_ped_edit.text())
        _opt(cmd, "-z", self._zerosup_edit.text())
        if self._prad1_chk.isChecked():
            cmd.append("-p")
        self._run_cmd(cmd, "Running replay recon…")

    # -- hadd merge step --

    def _run_hadd(self):
        import glob
        recon_dir = self._recon_dir
        if not recon_dir or not os.path.isdir(recon_dir):
            self._log("<span style='color:#f85149'>No recon directory for hadd.</span>")
            self._set_running(False)
            return

        root_files = sorted(glob.glob(os.path.join(recon_dir, "*.root")))
        if not root_files:
            self._log("<span style='color:#8b949e'>No ROOT files in recon dir, skipping hadd.</span>")
            self._run_next_step()
            return

        if len(root_files) == 1:
            self._hadd_out = root_files[0]
            self._log(f"<span style='color:#8b949e'>Single ROOT file, skipping hadd: {self._hadd_out}</span>")
            self._run_next_step()
            return

        out_name = self._out_name(root_files[0], "recon", "merged_recon.root")
        hadd_dir = os.path.dirname(recon_dir.rstrip("/")) or _RECON_BASE
        self._hadd_out = os.path.join(hadd_dir, out_name)
        self._hadd_inputs = list(root_files)

        cmd = ["hadd", "-f", self._hadd_out] + root_files
        self._run_cmd(cmd, "Merging ROOT files (hadd)…")

    # -- Replay filter step --

    def _run_filter(self):
        flt_input = self._filter_input_edit.text().strip() or self._hadd_out
        if not flt_input:
            self._log("<span style='color:#f85149'>No input specified for replay filter (run hadd first).</span>")
            self._set_running(False)
            return

        flt_out = self._filter_output_edit.text().strip()
        if not flt_out:
            flt_out = os.path.join(
                os.path.dirname(flt_input),
                self._out_name(flt_input, "filter", "filter_out.root"))
        self._filter_out = flt_out
        # JSON report path: same dir, replace .root → .report.json
        self._filter_report = re.sub(r'\.root$', '.report.json', flt_out,
                                      flags=re.IGNORECASE)

        cmd = [_REPLAY_FILTER_CMD, flt_input, "-o", flt_out,
               "-j", self._filter_report]
        _opt(cmd, "-c", self._filter_cut_edit.text())
        _opt(cmd, "-n", self._max_events_flt.text(), skip="-1")
        self._run_cmd(cmd, "Running replay filter…")

    # -- Quick check step --

    def _run_qcheck(self):
        qc_input = self._qc_input_edit.text().strip() or self._filter_out or self._hadd_out or self._recon_dir
        if not qc_input:
            self._log("<span style='color:#f85149'>No input specified for quick_check.</span>")
            self._set_running(False)
            return

        qc_out = self._qc_output_edit.text().strip()
        if not qc_out:
            tag = _run_tag(self._run_edit.text())
            out_name = f"{tag}_quick.root" if tag else "quick_check_out.root"
            base_dir = (os.path.dirname(qc_input)
                        if not os.path.isdir(qc_input) else qc_input)
            qc_out = os.path.join(base_dir, out_name)
        self._qcheck_out = qc_out

        cmd = [_QUICK_CHECK_CMD, qc_input]
        cmd += ["-o", qc_out]
        _opt(cmd, "-n", self._max_events_qc.text(), skip="-1")
        self._run_cmd(cmd, "Running quick check…")

    def _out_name(self, src: str, suffix: str, fallback: str) -> str:
        """prad_NNNNNN_<suffix>.root for the run field, else for the run in
        the name of ``src``, else ``fallback``."""
        tag = _run_tag(self._run_edit.text())
        if tag is None:
            m = re.match(r'(prad_\d+).*\.root', os.path.basename(src))
            tag = m.group(1) if m else None
        return f"{tag}_{suffix}.root" if tag else fallback

    # -- QProcess management --

    def _launch_process(self, cmd: List[str]):
        proc = QProcess(self)
        proc.readyReadStandardOutput.connect(self._on_stdout)
        proc.readyReadStandardError.connect(self._on_stderr)
        proc.finished.connect(self._on_finished)
        proc.errorOccurred.connect(self._on_process_error)
        proc.setProcessEnvironment(QProcessEnvironment.systemEnvironment())
        self._process = proc
        proc.start(cmd[0], cmd[1:])

    def _run_cmd(self, cmd: List[str], status: str):
        """Log ``cmd``, show ``status`` and launch it."""
        self._log(f"<span style='color:#8b949e'>$ {' '.join(cmd)}</span>")
        self._status_lbl.setText(status)
        self._launch_process(cmd)

    def _on_stdout(self):
        if self._process is None:
            return
        data = self._process.readAllStandardOutput().data().decode(errors="replace")
        self._log(html.escape(data, quote=False).replace("\n", "<br>"))

    def _on_stderr(self):
        if self._process is None:
            return
        data = self._process.readAllStandardError().data().decode(errors="replace")
        data = html.escape(data, quote=False).replace("\n", "<br>")
        self._log(f"<span style='color:#f0883e'>{data}</span>")

    def _on_process_error(self, error):
        labels = {
            QProcess.ProcessError.FailedToStart: "Failed to start (executable not found or no permission)",
            QProcess.ProcessError.Crashed:       "Process crashed",
            QProcess.ProcessError.Timedout:      "Timed out",
            QProcess.ProcessError.ReadError:     "Read error",
            QProcess.ProcessError.WriteError:    "Write error",
            QProcess.ProcessError.UnknownError:  "Unknown error",
        }
        msg = labels.get(error, f"Process error ({error})")
        self._log(f"<span style='color:#f85149'>[Process error] {msg}</span>")
        self._pending_steps.clear()
        self._process = None
        self._set_running(False)

    def _on_finished(self, exit_code: int, _status):
        self._process = None
        color = "#3fb950" if exit_code == 0 else "#f85149"
        self._log(f"<span style='color:{color}'>[Exit {exit_code}]</span>")

        # After hadd: always delete the individual recon ROOT files
        if exit_code == 0 and self._current_step == "hadd" and self._hadd_inputs:
            self._log("<span style='color:#8b949e'>Deleting individual recon ROOT files…</span>")
            for f in self._hadd_inputs:
                if os.path.isfile(f):
                    try:
                        os.remove(f)
                    except Exception as exc:
                        self._log(f"<span style='color:#f85149'>Delete failed ({f}): {exc}</span>")
            self._hadd_inputs = []
            self._log("<span style='color:#3fb950'>Individual recon files deleted.</span>")

        # After replay: optionally delete the downloaded EVIO directory
        if exit_code == 0 and self._current_step == "replay" \
                and self._auto_delete_evio_chk.isChecked():
            evio_path = self._evio_edit.text().strip() or self._evio_dir
            if evio_path and os.path.isdir(evio_path):
                self._log(f"<span style='color:#f0883e'>Auto-deleting EVIO dir: {evio_path}</span>")
                try:
                    shutil.rmtree(evio_path)
                    self._log("<span style='color:#3fb950'>EVIO files deleted.</span>")
                except Exception as exc:
                    self._log(f"<span style='color:#f85149'>Delete failed: {exc}</span>")

        # Emit filter report when filter step succeeded
        if exit_code == 0 and self._current_step == "filter" \
                and self._filter_report and os.path.isfile(self._filter_report):
            self.filterReportReady.emit(self._filter_report)

        if exit_code == 0:
            self._run_next_step()
        else:
            self._pending_steps.clear()
            self._set_running(False)
            # A failed step still loads any quick_check output that exists.
            if self._qcheck_out and os.path.isfile(self._qcheck_out):
                self.rootFileReady.emit(self._qcheck_out)

        # After any successful step, load the quick_check output once its
        # path is set and the file exists (possibly an earlier run's file
        # while quick_check is still running).
        if exit_code == 0 and self._qcheck_out and os.path.isfile(self._qcheck_out):
            self.rootFileReady.emit(self._qcheck_out)

    def _on_stop(self):
        self._pending_steps.clear()
        if self._process is not None:
            self._process.kill()
        self._set_running(False)

    def _set_running(self, running: bool):
        self._stop_btn.setEnabled(running)
        for btn in (self._scp_btn, self._replay_btn, self._filter_btn, self._qcheck_btn, self._do_all_btn):
            btn.setEnabled(not running)
        if not running:
            self._status_lbl.setText("Ready")

    # -- Console helpers --

    def _log(self, html: str):
        self._console.moveCursor(self._console.textCursor().MoveOperation.End)
        self._console.insertHtml(html)
        self._console.insertHtml("<br>")
        self._console.moveCursor(self._console.textCursor().MoveOperation.End)


# ---- Filter Report Widget ----

def _mpl_theme_rc() -> dict:
    """matplotlib rcParams for the active THEME."""
    return {"axes.facecolor": THEME.CANVAS, "axes.edgecolor": THEME.BORDER,
            "axes.labelcolor": THEME.TEXT, "text.color": THEME.TEXT,
            "xtick.color": THEME.TEXT, "ytick.color": THEME.TEXT,
            "legend.facecolor": THEME.PANEL, "legend.edgecolor": THEME.BORDER,
            "legend.labelcolor": THEME.TEXT}


class FilterReportWidget(QWidget):
    """Embedded replay-filter JSON report chart.

    Shows three rows (cut status / livetime+rate / EPICS) drawn by
    replay_report_viewer.plot_report.  Requires matplotlib and
    replay_report_viewer to be importable; shows a placeholder otherwise.
    """

    def __init__(self, parent=None):
        super().__init__(parent)
        v = QVBoxLayout(self)
        v.setContentsMargins(4, 4, 4, 4)
        v.setSpacing(4)

        if not (HAS_REPORT_VIEWER and HAS_MATPLOTLIB):
            lbl = QLabel(
                "Filter report viewer requires matplotlib.\n"
                "Install it with:  pip install matplotlib")
            lbl.setAlignment(Qt.AlignmentFlag.AlignCenter)
            lbl.setStyleSheet("QLabel{color:#8b949e;font-size:11pt;}")
            v.addWidget(lbl)
            self._available = False
            return

        self._available = True
        self._report = None

        # ── top controls ──────────────────────────────────────────────────
        top = QHBoxLayout()

        top.addWidget(QLabel("x-axis:"))
        self._cmb_x = QComboBox()
        self._cmb_x.addItem("associated_timestamp", "time")
        self._cmb_x.addItem("associated_evn",       "evn")
        self._cmb_x.setStyleSheet(_COMBO_SS)
        self._cmb_x.setFont(QFont("Consolas", 10))
        self._cmb_x.currentIndexChanged.connect(self._replot)
        top.addWidget(self._cmb_x)

        top.addSpacing(12)
        self._btn_epics = QToolButton()
        self._btn_epics.setText("EPICS \u25be")
        self._btn_epics.setPopupMode(QToolButton.ToolButtonPopupMode.InstantPopup)
        self._menu_epics = QMenu(self._btn_epics)
        self._btn_epics.setMenu(self._menu_epics)
        top.addWidget(self._btn_epics)

        top.addStretch()
        self._lbl_title = QLabel("")
        self._lbl_title.setStyleSheet(_LBL_MUT)
        top.addWidget(self._lbl_title)

        v.addLayout(top)

        # ── matplotlib canvas ─────────────────────────────────────────────
        self._fig = MplFigure(constrained_layout=True, facecolor=THEME.CANVAS)
        self._canvas = FigureCanvasQTAgg(self._fig)
        self._canvas.setSizePolicy(
            QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        v.addWidget(self._canvas, 1)

        self._nav = _MplNavToolbar(self._canvas, self)
        v.addWidget(self._nav)

        self._selected_epics: set = set()

    def load_report(self, path: str) -> None:
        """Load a *.report.json written by prad2ana_replay_filter."""
        if not self._available:
            return
        try:
            self._report = _rrv_load_report(path)
        except Exception as exc:
            self._lbl_title.setText(
                f"<span style='color:#f85149'>Load error: {exc}</span>")
            return

        self._selected_epics = _rrv_fill_epics_menu(
            self._menu_epics, self._report, self, self._on_epics_toggled)

        self._lbl_title.setText(_rrv_title_for(self._report))
        self._replot()

    def _on_epics_toggled(self) -> None:
        self._selected_epics = _rrv_checked_epics(self._menu_epics)
        self._replot()

    def _replot(self) -> None:
        self._fig.clear()
        self._fig.patch.set_facecolor(THEME.CANVAS)
        r = self._report
        if r is not None:
            sel = _rrv_selected_series(r, self._selected_epics)
            with rc_context(_mpl_theme_rc()):
                _rrv_plot_report(self._fig, r, self._cmb_x.currentData() or "time",
                                 sel, note_color=THEME.TEXT_DIM)
        self._canvas.draw_idle()


# ---- Results Panel (right side) ----

class ResultsPanel(QWidget):
    """Tabbed display of quick_check ROOT file contents."""

    def __init__(self, parent=None):
        super().__init__(parent)
        self._modules = []
        self._loader: Optional[_RootLoader] = None
        self._build_ui()
        self._try_load_modules()

    def _try_load_modules(self):
        if os.path.isfile(MODULES_JSON):
            try:
                self._modules = load_modules(MODULES_JSON)
                self._map_widget.set_modules(self._modules)
                self._map_widget.set_palette("viridis")
            except Exception:
                pass

    def _build_ui(self):
        root = QVBoxLayout(self)
        root.setContentsMargins(4, 4, 4, 4)
        root.setSpacing(4)

        # header
        hdr = QHBoxLayout()
        self._file_lbl = QLabel("No file loaded")
        self._file_lbl.setStyleSheet(themed(
            "QLabel{color:#8b949e;font-family:Consolas;font-size:9pt;}"))
        self._reload_btn = QPushButton("Reload")
        self._reload_btn.setFixedWidth(90)
        self._reload_btn.setStyleSheet(_BTN_NORMAL)
        self._reload_btn.clicked.connect(self._reload)
        self._reload_btn.setEnabled(False)
        hdr.addWidget(self._file_lbl)
        hdr.addStretch()
        hdr.addWidget(self._reload_btn)
        root.addLayout(hdr)

        # tabs
        self._tabs = QTabWidget()
        self._tabs.setStyleSheet(themed(
            "QTabWidget::pane{border:1px solid #30363d;background:#0d1117;}"
            "QTabBar::tab{background:#161b22;color:#8b949e;"
            "padding:4px 14px;border:1px solid #30363d;"
            "border-bottom:none;font-family:Consolas;font-size:10pt;}"
            "QTabBar::tab:selected{background:#0d1117;color:#c9d1d9;}"
            "QTabBar::tab:hover{background:#21262d;color:#c9d1d9;}"))

        # Tab 0: HyCal map
        self._map_widget = HyCalReplayMapWidget()
        map_container = QWidget()
        map_lay = QVBoxLayout(map_container)
        map_lay.setContentsMargins(4, 4, 4, 4)

        map_ctrl = QHBoxLayout()
        map_mode_lbl = QLabel("Show:")
        map_mode_lbl.setStyleSheet(_LBL_SS)
        self._map_mode_combo = QComboBox()
        self._map_mode_combo.addItems(["Module hits (count)", "Mean energy"])
        self._map_mode_combo.setStyleSheet(_COMBO_SS)
        self._map_mode_combo.setFont(QFont("Consolas", 10))
        self._map_mode_combo.currentIndexChanged.connect(self._refresh_map)
        pal_lbl = QLabel("Palette:")
        pal_lbl.setStyleSheet(_LBL_SS)
        self._pal_combo = QComboBox()
        self._pal_combo.addItems(PALETTE_NAMES)
        self._pal_combo.setCurrentText("viridis")
        self._pal_combo.setStyleSheet(_COMBO_SS)
        self._pal_combo.setFont(QFont("Consolas", 10))
        self._pal_combo.currentIndexChanged.connect(
            lambda i: self._map_widget.set_palette(PALETTE_NAMES[i]))
        map_ctrl.addWidget(map_mode_lbl)
        map_ctrl.addWidget(self._map_mode_combo)
        map_ctrl.addSpacing(12)
        map_ctrl.addWidget(pal_lbl)
        map_ctrl.addWidget(self._pal_combo)
        map_ctrl.addStretch()

        map_lay.addLayout(map_ctrl)
        map_lay.addWidget(self._map_widget, stretch=1)
        self._tabs.addTab(map_container, "HyCal Map")

        # Tab 1: Hit Position 2D
        self._hit_pos_widget = Hist2DWidget("Hit Position")
        self._tabs.addTab(self._hit_pos_widget, "Hit Position")

        # Tab 2: Energy Spectra
        energy_tab = QWidget()
        eg = QVBoxLayout(energy_tab)
        eg.setContentsMargins(4, 4, 4, 4)
        eg.setSpacing(4)
        self._h_1cl   = Hist1DWidget("1-cluster energy")
        self._h_2cl   = Hist1DWidget("2-cluster energy")
        self._h_all   = Hist1DWidget("All clusters energy")
        self._h_tot   = Hist1DWidget("Total energy")
        # Give each spectrum a distinct, readable colour
        self._h_1cl.BAR_COLOR = "#3fb950"   # green
        self._h_2cl.BAR_COLOR = "#58a6ff"   # blue
        self._h_all.BAR_COLOR = "#d29922"   # amber
        self._h_tot.BAR_COLOR = "#f85149"   # red
        for h in (self._h_1cl, self._h_2cl, self._h_all, self._h_tot):
            eg.addWidget(h)
        self._tabs.addTab(energy_tab, "Energy Spectra")

        # Tab 3: Energy vs Theta 2D
        self._ev_theta_widget = Hist2DWidget("Energy vs θ")
        self._tabs.addTab(self._ev_theta_widget, "Energy vs Theta")

        # Tab 4: Moller Analysis
        moller_tab = QWidget()
        mg = QVBoxLayout(moller_tab)
        mg.setContentsMargins(4, 4, 4, 4)
        mg.setSpacing(4)
        mol_top = QHBoxLayout()
        self._h_moller_z    = Hist1DWidget("Moller Z vertex")
        self._h_moller_z.auto_cb_fit = True
        self._h_moller_z.cb_fit_range_sigma = (3.0, 1.5)  # left 3σ, right 1.5σ
        self._h_moller_phi  = Hist1DWidget("Moller Φ diff")
        mol_top.addWidget(self._h_moller_z)
        mol_top.addWidget(self._h_moller_phi)
        mol_bot = QHBoxLayout()
        self._h_moller_x = Hist1DWidget("Moller X center")
        self._h_moller_x.auto_cb_fit = True
        self._h_moller_x.cb_fit_range = (5.0, 2.5)  # left 5 mm, right 2.5 mm from peak
        self._h_moller_y = Hist1DWidget("Moller Y center")
        self._h_moller_y.auto_cb_fit = True
        self._h_moller_y.cb_fit_range = (5.0, 2.5)  # left 5 mm, right 2.5 mm from peak
        mol_bot.addWidget(self._h_moller_x)
        mol_bot.addWidget(self._h_moller_y)
        mg.addLayout(mol_top, stretch=1)
        mg.addLayout(mol_bot, stretch=1)
        # 2-arm Moller position 2D
        self._moller_2arm = Hist2DWidget("2-arm Moller position")
        mg.addWidget(self._moller_2arm, stretch=2)
        self._tabs.addTab(moller_tab, "Moller")

        # Tab 5: Physics Yields
        yields_tab = QWidget()
        yg = QVBoxLayout(yields_tab)
        yg.setContentsMargins(4, 4, 4, 4)
        yg.setSpacing(4)
        self._h_ep    = Hist1DWidget("ep yield")
        self._h_ee    = Hist1DWidget("ee yield")
        self._h_ratio = Hist1DWidget("ep/ee ratio")
        self._h_ep.BAR_COLOR = "#3a86ff"
        self._h_ee.BAR_COLOR = "#ff6b6b"
        self._h_ratio.BAR_COLOR = "#ffd166"
        yg.addWidget(self._h_ep)
        yg.addWidget(self._h_ee)
        yg.addWidget(self._h_ratio)
        self._tabs.addTab(yields_tab, "Physics Yields")

        # Tab 6: Filter Report
        self._filter_report_widget = FilterReportWidget()
        self._tabs.addTab(self._filter_report_widget, "Filter Report")

        root.addWidget(self._tabs, stretch=1)

        # status bar
        self._loading_lbl = QLabel("")
        self._loading_lbl.setStyleSheet(_LBL_MUT)
        root.addWidget(self._loading_lbl)

    # -- Public API --

    def load_file(self, path: str):
        if not HAS_UPROOT:
            self._file_lbl.setText("uproot not installed — cannot read ROOT files")
            self._loading_lbl.setText(
                "Install: pip install uproot numpy")
            return
        if not os.path.isfile(path):
            self._loading_lbl.setText(f"File not found: {path}")
            return
        self._current_path = path
        self._file_lbl.setText(os.path.basename(path))
        self._reload_btn.setEnabled(True)
        self._loading_lbl.setText("Loading…")
        self._loader = _RootLoader(path, self)
        self._loader.finished.connect(self._on_loaded)
        self._loader.start()

    def load_filter_report(self, path: str):
        """Load a replay-filter JSON report and switch to the Filter Report tab."""
        self._filter_report_widget.load_report(path)
        idx = self._tabs.indexOf(self._filter_report_widget)
        if idx >= 0:
            self._tabs.setCurrentIndex(idx)

    def _reload(self):
        if hasattr(self, "_current_path"):
            self.load_file(self._current_path)

    def _on_loaded(self, data: dict, error: str):
        self._loader = None
        if error:
            self._loading_lbl.setText(f"Error: {error}")
            return
        self._data = data
        self._loading_lbl.setText(
            f"Loaded: {len(data)} datasets")
        self._populate(data)

    # -- Populate widgets from loaded data --

    def _populate(self, data: dict):
        self._refresh_map()
        for key, attr in _HIST_PLOTS:
            widget = getattr(self, attr)
            # skip a histogram whose dimension does not fit the widget
            if len(data.get(key, ())) == (3 if isinstance(widget, Hist2DWidget) else 2):
                widget.set_data(*data[key])

    def _refresh_map(self):
        if not hasattr(self, "_data"):
            return
        data = self._data
        mode = self._map_mode_combo.currentIndex()
        if mode == 1 and data.get("module_means"):
            self._map_widget.set_module_counts(data["module_means"])
        else:
            self._map_widget.set_module_counts(data.get("module_counts", {}))


# ---- Main Window ----

_MAIN_QSS = (
    "QMainWindow{background:#0d1117;}"
    "QWidget{background:#0d1117;color:#c9d1d9;}"
    "QSplitter::handle{background:#21262d;}"
    "QLabel{color:#c9d1d9;}")

_MENUBAR_QSS = (
    "QMenuBar{background:#161b22;color:#c9d1d9;"
    "font-family:Consolas;font-size:10pt;}"
    "QMenuBar::item:selected{background:#21262d;}"
    "QMenu{background:#161b22;color:#c9d1d9;border:1px solid #30363d;}"
    "QMenu::item:selected{background:#1f6feb;}")


class MainWindow(QMainWindow):
    def __init__(self, initial_root: str = ""):
        super().__init__()
        self.setWindowTitle("PRad-2 Replay Viewer")
        self.resize(1600, 900)
        self._apply_main_qss()

        # menu bar
        mb = self.menuBar()
        file_menu = mb.addMenu("File")
        open_act = file_menu.addAction("Open ROOT file…")
        open_act.triggered.connect(self._open_root)
        file_menu.addSeparator()
        quit_act = file_menu.addAction("Quit")
        quit_act.triggered.connect(self.close)

        view_menu = mb.addMenu("View")
        theme_menu = view_menu.addMenu("Theme")
        for t in available_themes():
            act = theme_menu.addAction(t.capitalize())
            act.triggered.connect(lambda _c=False, tn=t: self._change_theme(tn))

        # main splitter
        splitter = QSplitter(Qt.Orientation.Horizontal)
        splitter.setHandleWidth(4)

        self._ctrl = ControlPanel()
        self._results = ResultsPanel()

        self._ctrl.rootFileReady.connect(self._results.load_file)
        self._ctrl.filterReportReady.connect(self._results.load_filter_report)

        splitter.addWidget(self._ctrl)
        splitter.addWidget(self._results)
        splitter.setSizes([420, 1180])
        splitter.setCollapsible(0, False)
        splitter.setCollapsible(1, False)

        self.setCentralWidget(splitter)

        if initial_root:
            self._results.load_file(initial_root)

    def _open_root(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open ROOT File", "", "ROOT files (*.root);;All files (*)")
        if path:
            self._results.load_file(path)

    def _apply_main_qss(self):
        self.setStyleSheet(themed(_MAIN_QSS))
        self.menuBar().setStyleSheet(themed(_MENUBAR_QSS))

    def _change_theme(self, name: str):
        set_theme(name)
        apply_theme_palette(QApplication.instance())
        self._apply_main_qss()


# ---- Entry point ----

def main():
    app = QApplication(sys.argv)
    app.setApplicationName("PRad-2 Replay Viewer")

    set_theme("dark")
    apply_theme_palette(app)

    initial = sys.argv[1] if len(sys.argv) > 1 else ""
    win = MainWindow(initial_root=initial)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
