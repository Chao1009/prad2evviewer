#!/usr/bin/env python3
"""
energy_corr_viewer.py — HyCal PWO4 energy-correction viewer
============================================================

ep and ee have **independent** correction coefficients.

Left  : HyCal map coloured by the ep or ee mean correction (toggle).
        Top-right overlay: two stacked 5×5 mini-grids — ep (top) / ee (bottom)
        — for the selected module, coloured by their respective corrections.

Right : Two 5×5 histogram blocks (ep block on top, ee block on bottom) for the
        selected module.  Each cell shows its own Gaussian fit and correction.
        Drag on any subplot to set a custom fit range (auto-refit on release).

Controls
  • Beam-energy buttons: 3.5 / 2.2 / 0.7 GeV
  • Map colour mode: ep mean / ee mean / ep max|dev| / ee max|dev|
  • Click a histogram cell to select it
  • Edit ep and ee corrections independently in the control panel
  • R / 1 resets the ep correction to 1.0 for the selected cell
  • E resets the ee correction to 1.0 for the selected cell
  • Export JSON writes separate ep and ee corrections

Logic follows energy_corr.cpp:
  ep fit range : [E_ep ± 2.5σ]  where σ = 0.035·E / √(E/GeV)
  ee fit range : [E_ee ± 2.5σ_ee]
  correction   = E_expected / fit_mean, capped to [0.92, 1.08]

Usage:
  python scripts/energy_corr_viewer.py [root_file] [--json corr.json] [--theme dark|light]
"""
from __future__ import annotations

import argparse
import json
import math
import re
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

from PyQt6.QtCore import Qt
from PyQt6.QtGui import QFont, QKeySequence, QShortcut
from PyQt6.QtWidgets import (
    QApplication, QComboBox, QFileDialog, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMainWindow, QPushButton,
    QScrollArea, QSizePolicy, QSplitter, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector
import matplotlib.patches as mpatches
import matplotlib.colors as mcolors
import matplotlib.cm as cm
import numpy as np

try:
    import uproot
    HAS_UPROOT = True
except ImportError:
    HAS_UPROOT = False

try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from hycal_geoview import (
    HyCalMapWidget, Module, load_modules,
    apply_theme_palette, set_theme, available_themes, THEME,
    ColorRangeControl,
)

# ── Database ──────────────────────────────────────────────────────────────────
_DB_CANDIDATES = [
    _SCRIPT_DIR.parent / "database" / "hycal_map.json",
    _SCRIPT_DIR.parent / "build" / "database" / "hycal_map.json",
]
MODULES_JSON = next((p for p in _DB_CANDIDATES if p.is_file()), None)

# ── Physics ───────────────────────────────────────────────────────────────────
M_PROTON    = 938.272
M_ELECTRON  = 0.511
Z_HYCAL     = 6270.0
RESOLUTION  = 0.035
GRID        = 5
CORR_CAP_HI = 1.08
CORR_CAP_LO = 0.92
EP_COLOR    = "#ff6655"
EE_COLOR    = "#6688ff"

BEAM_CONFIGS: List[Tuple[str, str, float]] = [
    ("3.5 GeV", "3p5GeV", 3485.41),
    ("2.2 GeV", "2p2GeV", 2239.51),
    ("0.7 GeV", "0p7GeV",  728.9),
]


def _energy_loss(theta_deg: float) -> float:
    theta = math.radians(theta_deg)
    cos_t = math.cos(theta)
    sec   = 1.0 / cos_t if cos_t > 0.01 else 100.0
    return (0.500 * 1.6 + 0.120 * 1.6 + 0.100 * 2.0 + 0.480 * 1.8) * sec


def _gauss(x, amp, mu, sig):
    return amp * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def fit_gaussian(counts: np.ndarray, edges: np.ndarray,
                 mu0: float, half_width: float) -> Tuple[float, float, float, bool]:
    if not HAS_SCIPY or counts is None or len(counts) == 0:
        return 0., 0., 0., False
    centers = 0.5 * (edges[:-1] + edges[1:])
    lo, hi  = mu0 - half_width, mu0 + half_width
    mask    = (centers >= lo) & (centers <= hi) & (counts > 0)
    if mask.sum() < 4:
        return 0., 0., 0., False
    xf, yf = centers[mask], counts[mask].astype(float)
    bw  = float(centers[1] - centers[0]) if len(centers) > 1 else 1.0
    sig0 = max(half_width / 3.0, bw)
    try:
        popt, _ = curve_fit(
            _gauss, xf, yf,
            p0=[float(yf.max()), mu0, sig0],
            bounds=([0., lo, 1e-6], [np.inf, hi, half_width * 2.]),
            maxfev=3000,
        )
        mean, sigma, amp = float(popt[1]), abs(float(popt[2])), float(popt[0])
        if math.isfinite(mean) and math.isfinite(sigma) and sigma > 0.:
            return mean, sigma, amp, True
    except Exception:
        pass
    return 0., 0., 0., False


# ── Data structures ───────────────────────────────────────────────────────────

class CellData:
    """Per-cell histogram with independent ep and ee Gaussian fits/corrections."""
    __slots__ = [
        "hist_counts", "hist_edges", "expected_ep", "expected_ee", "entries",
        "fit_mean_ep", "fit_sigma_ep", "fit_amp_ep", "fit_valid_ep",
        "corr_ep", "corr_raw_ep", "span_lo_ep", "span_hi_ep",
        "fit_mean_ee", "fit_sigma_ee", "fit_amp_ee", "fit_valid_ee",
        "corr_ee", "corr_raw_ee", "span_lo_ee", "span_hi_ee",
    ]

    def __init__(self):
        self.hist_counts: Optional[np.ndarray] = None
        self.hist_edges:  Optional[np.ndarray] = None
        self.expected_ep: float = 0.; self.expected_ee: float = 0.
        self.entries:     int   = 0
        self.fit_mean_ep:  float = 0.; self.fit_sigma_ep: float = 0.
        self.fit_amp_ep:   float = 0.; self.fit_valid_ep: bool  = False
        self.corr_ep:      float = 1.; self.corr_raw_ep:  float = 0.
        self.span_lo_ep:   Optional[float] = None
        self.span_hi_ep:   Optional[float] = None
        self.fit_mean_ee:  float = 0.; self.fit_sigma_ee: float = 0.
        self.fit_amp_ee:   float = 0.; self.fit_valid_ee: bool  = False
        self.corr_ee:      float = 1.; self.corr_raw_ee:  float = 0.
        self.span_lo_ee:   Optional[float] = None
        self.span_hi_ee:   Optional[float] = None

    def corr(self, kind: str) -> float:
        return self.corr_ep if kind == "ep" else self.corr_ee

    def sigma_exp(self, kind: str) -> float:
        E = self.expected_ep if kind == "ep" else self.expected_ee
        if E <= 0.:
            return 1.
        return RESOLUTION * E / math.sqrt(max(E, 1.) / 1000.)

    def refit(self, kind: str) -> bool:
        """Three-stage fit matching energy_corr.cpp / nonlinearity.cpp:
        1. Weighted centroid in [Eexp ± 6σ]
        2. Gaussian fit in [mean ± 3σ]
        3. Final Gaussian fit in [mean ± 2σ]
        If a custom span is set, use a single fit in that range instead.
        """
        Eexp    = self.expected_ep if kind == "ep" else self.expected_ee
        span_lo = self.span_lo_ep  if kind == "ep" else self.span_lo_ee
        span_hi = self.span_hi_ep  if kind == "ep" else self.span_hi_ee

        if self.hist_counts is None or Eexp <= 0.:
            if kind == "ep": self.fit_valid_ep = False
            else:            self.fit_valid_ee = False
            return False

        sig_e = self.sigma_exp(kind)

        if span_lo is not None and span_hi is not None:
            # User-specified range: single Gaussian fit
            mu0    = 0.5 * (span_lo + span_hi)
            half_w = 0.5 * (span_hi - span_lo)
            mean, sigma, amp, ok = fit_gaussian(
                self.hist_counts, self.hist_edges, mu0, half_w)
        else:
            # Auto: 3-stage approach
            centers = 0.5 * (self.hist_edges[:-1] + self.hist_edges[1:])

            # Stage 1: weighted centroid in [Eexp ± 6σ]
            mask6 = ((centers >= Eexp - 6. * sig_e) &
                     (centers <= Eexp + 6. * sig_e) &
                     (self.hist_counts > 0))
            if not mask6.any():
                if kind == "ep": self.fit_valid_ep = False
                else:            self.fit_valid_ee = False
                return False
            wsum = float(self.hist_counts[mask6].sum())
            mean_seed = float(
                (centers[mask6] * self.hist_counts[mask6]).sum() / wsum)

            # Stage 2: Gaussian fit in [mean_seed ± 3σ]
            sig_seed = RESOLUTION * mean_seed / math.sqrt(max(mean_seed, 1.) / 1000.)
            m2, _s2, _a2, ok2 = fit_gaussian(
                self.hist_counts, self.hist_edges, mean_seed, 3. * sig_seed)
            if ok2:
                mean_seed = m2
                sig_seed  = RESOLUTION * mean_seed / math.sqrt(max(mean_seed, 1.) / 1000.)

            # Stage 3: final Gaussian fit in [mean_seed ± 2σ]
            mean, sigma, amp, ok = fit_gaussian(
                self.hist_counts, self.hist_edges, mean_seed, 2. * sig_seed)

        if kind == "ep":
            self.fit_valid_ep = ok
            if ok:
                self.fit_mean_ep = mean; self.fit_sigma_ep = sigma
                self.fit_amp_ep  = amp
                raw = Eexp / mean if mean > 0. else 0.
                self.corr_raw_ep = raw
                self.corr_ep = (max(CORR_CAP_LO, min(CORR_CAP_HI, raw))
                                if raw > 0. else 1.)
        else:
            self.fit_valid_ee = ok
            if ok:
                self.fit_mean_ee = mean; self.fit_sigma_ee = sigma
                self.fit_amp_ee  = amp
                raw = Eexp / mean if mean > 0. else 0.
                self.corr_raw_ee = raw
                self.corr_ee = (max(CORR_CAP_LO, min(CORR_CAP_HI, raw))
                                if raw > 0. else 1.)
        return ok


class ModuleEnergyData:
    def __init__(self):
        self.cells: Dict[Tuple[int, int], CellData] = {
            (r, c): CellData() for r in range(GRID) for c in range(GRID)
        }

    def mean_correction(self, kind: str) -> float:
        attr = "corr_ep" if kind == "ep" else "corr_ee"
        vals = [getattr(cd, attr) for cd in self.cells.values()
                if math.isfinite(getattr(cd, attr)) and getattr(cd, attr) > 0.]
        return float(np.mean(vals)) if vals else 1.0

    def max_deviation(self, kind: str) -> float:
        attr = "corr_ep" if kind == "ep" else "corr_ee"
        devs = [abs(getattr(cd, attr) - 1.) for cd in self.cells.values()
                if math.isfinite(getattr(cd, attr)) and getattr(cd, attr) > 0.]
        return float(max(devs)) if devs else 0.

    def valid_count(self, kind: str) -> int:
        attr = "fit_valid_ep" if kind == "ep" else "fit_valid_ee"
        return sum(1 for cd in self.cells.values() if getattr(cd, attr))


class ModuleData:
    def __init__(self, name: str, x: float, y: float):
        self.name  = name
        self.x     = float(x)
        self.y     = float(y)
        self.theta = math.degrees(math.atan2(math.hypot(x, y), Z_HYCAL))
        self.energies: Dict[str, ModuleEnergyData] = {
            sd: ModuleEnergyData() for _, sd, _ in BEAM_CONFIGS
        }


# ── Data loading ──────────────────────────────────────────────────────────────

def load_root(root_path: Path) -> Dict[str, ModuleData]:
    data: Dict[str, ModuleData] = {}
    if MODULES_JSON is not None:
        try:
            with open(MODULES_JSON) as f:
                geo = json.load(f)
            for e in geo:
                if str(e.get("t", "")).lower() != "pbwo4":
                    continue
                g = e.get("geo", {})
                data[str(e["n"])] = ModuleData(str(e["n"]),
                                               float(g["x"]), float(g["y"]))
        except Exception as ex:
            print(f"Warning: geometry: {ex}", file=sys.stderr)

    if not root_path.is_file():
        return data
    if not HAS_UPROOT:
        print("Warning: uproot not installed.", file=sys.stderr)
        return data

    try:
        with uproot.open(str(root_path)) as f:
            for _, subdir, _ in BEAM_CONFIGS:
                if subdir not in f:
                    continue
                sd = f[subdir]

                # ── ep fits from TTree ──
                if "fit_results" in sd:
                    t = sd["fit_results"]
                    try:
                        mod_names = [str(s) for s in
                                     t["module_name"].array(library="np").tolist()]
                    except Exception:
                        mod_names = [str(s) for s in
                                     t["module_name"].array(library="ak").tolist()]
                    rows_a  = t["row"].array(library="np")
                    cols_a  = t["col"].array(library="np")
                    ep_a    = t["expected_energy_ep"].array(library="np")
                    ee_a    = t["expected_energy_ee"].array(library="np")
                    corr_a  = t["correction"].array(library="np")
                    valid_a = t["fit_valid"].array(library="np")
                    ent_a   = t["entries"].array(library="np")
                    mean_a  = t["fit_mean"].array(library="np")
                    sigma_a = t["fit_sigma"].array(library="np")

                    # ee fit branches — present only in newer ROOT files
                    try:
                        corr_ee_a  = t["correction_ee"].array(library="np")
                        valid_ee_a = t["fit_valid_ee"].array(library="np")
                        mean_ee_a  = t["fit_mean_ee"].array(library="np")
                        sigma_ee_a = t["fit_sigma_ee"].array(library="np")
                        has_ee_branch = True
                    except Exception:
                        has_ee_branch = False

                    for i in range(len(mod_names)):
                        name = mod_names[i]
                        if name not in data:
                            continue
                        r, c = int(rows_a[i]), int(cols_a[i])
                        cell = data[name].energies[subdir].cells[(r, c)]
                        cell.expected_ep = float(ep_a[i])
                        cell.expected_ee = float(ee_a[i])
                        cell.entries     = int(ent_a[i])

                        # ep
                        raw_ep = float(corr_a[i])
                        cell.fit_valid_ep = bool(valid_a[i])
                        cell.corr_raw_ep  = raw_ep
                        if cell.fit_valid_ep and raw_ep > 0. and math.isfinite(raw_ep):
                            cell.fit_mean_ep  = float(mean_a[i])
                            cell.fit_sigma_ep = float(sigma_a[i])
                            cell.corr_ep = max(CORR_CAP_LO, min(CORR_CAP_HI, raw_ep))
                        else:
                            cell.corr_ep = 1.0

                        # ee (from TTree if available)
                        if has_ee_branch:
                            raw_ee = float(corr_ee_a[i])
                            cell.fit_valid_ee = bool(valid_ee_a[i])
                            cell.corr_raw_ee  = raw_ee
                            if cell.fit_valid_ee and raw_ee > 0. and math.isfinite(raw_ee):
                                cell.fit_mean_ee  = float(mean_ee_a[i])
                                cell.fit_sigma_ee = float(sigma_ee_a[i])
                                cell.corr_ee = max(CORR_CAP_LO, min(CORR_CAP_HI, raw_ee))
                            else:
                                cell.corr_ee = 1.0

                # ── per-cell histograms ──
                for key_obj in sd:
                    key_str = str(key_obj).split(";")[0]
                    if not key_str.startswith("h_energy_"):
                        continue
                    m = re.match(r"^h_energy_(.+)_r(\d+)_c(\d+)$", key_str)
                    if not m:
                        continue
                    mod_name = m.group(1)
                    r, c = int(m.group(2)), int(m.group(3))
                    if mod_name not in data:
                        continue
                    try:
                        counts, edges = sd[str(key_obj)].to_numpy()
                        cell = data[mod_name].energies[subdir].cells[(r, c)]
                        cell.hist_counts = counts.astype(float)
                        cell.hist_edges  = edges.astype(float)
                        if counts.max() > 0:
                            cell.fit_amp_ep = float(counts.max())
                    except Exception:
                        pass

                # ── refit all cells on load (ep + ee) using 3-stage method ──
                n_ep_ok = n_ee_ok = 0
                for mod in data.values():
                    for cell in mod.energies[subdir].cells.values():
                        if cell.hist_counts is not None:
                            if cell.expected_ep > 0. and cell.refit("ep"):
                                n_ep_ok += 1
                            if cell.expected_ee > 0. and cell.refit("ee"):
                                n_ee_ok += 1
                print(f"[{subdir}] refit on load: ep={n_ep_ok} ee={n_ee_ok}",
                      file=sys.stderr)

    except Exception as ex:
        print(f"Warning: ROOT read: {ex}", file=sys.stderr)

    return data


def apply_json_corrections(data: Dict[str, ModuleData], json_path: Path) -> None:
    try:
        with open(json_path) as f:
            jdata = json.load(f)
        for mod_name, mod_json in jdata.items():
            if mod_name not in data:
                continue
            for _, subdir, _ in BEAM_CONFIGS:
                if subdir not in mod_json:
                    continue
                ep_corr = mod_json[subdir].get("ep", {}).get("correction")
                ee_corr = mod_json[subdir].get("ee", {}).get("correction")
                for r in range(GRID):
                    for c in range(GRID):
                        cell = data[mod_name].energies[subdir].cells[(r, c)]
                        if ep_corr and r < len(ep_corr) and c < len(ep_corr[r]):
                            cell.corr_ep = float(ep_corr[r][c])
                        if ee_corr and r < len(ee_corr) and c < len(ee_corr[r]):
                            cell.corr_ee = float(ee_corr[r][c])
    except Exception as ex:
        print(f"Warning: JSON read: {ex}", file=sys.stderr)


def _fmt(v: float, dec: int) -> str:
    s = f"{v:.{dec}f}".rstrip("0")
    if s.endswith("."):
        s += "0"
    return s


def save_json(data: Dict[str, ModuleData], json_path: Path) -> None:
    id_order: Dict[str, int] = {}
    if MODULES_JSON is not None:
        try:
            with open(MODULES_JSON) as f:
                geo = json.load(f)
            for idx, e in enumerate(geo):
                id_order[str(e.get("n", ""))] = idx
        except Exception:
            pass

    def skey(name: str) -> int:
        return id_order.get(name, 999999)

    lines: List[str] = ["{\n"]
    mod_list = sorted(data.keys(), key=skey)
    for mi, mod_name in enumerate(mod_list):
        mod = data[mod_name]
        lines.append(f'  "{mod_name}": {{\n')
        for ei, (_, subdir, _) in enumerate(BEAM_CONFIGS):
            ed = mod.energies[subdir]
            lines.append(f'    "{subdir}": {{\n')
            for pi, pkey in enumerate(("ep", "ee")):
                lines.append(f'      "{pkey}": {{\n')
                for field in ("correction", "expected_energy"):
                    mat_rows = []
                    for r in range(GRID):
                        vals = []
                        for c in range(GRID):
                            cell = ed.cells[(r, c)]
                            if field == "correction":
                                v = cell.corr_ep if pkey == "ep" else cell.corr_ee
                                vals.append(_fmt(v, 4))
                            else:
                                v = (cell.expected_ep if pkey == "ep"
                                     else cell.expected_ee)
                                vals.append(_fmt(v, 1))
                        mat_rows.append("          [" + ", ".join(vals) + "]")
                    mat_str = "[\n" + ",\n".join(mat_rows) + "\n        ]"
                    comma = "," if field == "correction" else ""
                    lines.append(f'        "{field}": {mat_str}{comma}\n')
                comma_p = "," if pi == 0 else ""
                lines.append(f'      }}{comma_p}\n')
            comma_e = "," if ei < len(BEAM_CONFIGS) - 1 else ""
            lines.append(f'    }}{comma_e}\n')
        comma_m = "," if mi < len(mod_list) - 1 else ""
        lines.append(f'  }}{comma_m}\n')
    lines.append("}\n")
    with open(json_path, "w") as fout:
        fout.writelines(lines)
    print(f"Saved {len(mod_list)} modules to {json_path}", file=sys.stderr)


# ── Mini 5×5 × 2 overlay canvas ──────────────────────────────────────────────

class MiniGridCanvas(FigureCanvas):
    """Two side-by-side 5×5 correction grids: ep (left) and ee (right)."""
    FIXED_W = 340
    FIXED_H = 200

    def __init__(self, parent=None):
        bg = "#1a1a1a"
        self.fig = Figure(figsize=(3.4, 2.0), facecolor=bg)
        self.ax_ep = self.fig.add_axes([0.03, 0.10, 0.44, 0.82])
        self.ax_ee = self.fig.add_axes([0.53, 0.10, 0.44, 0.82])
        for ax in (self.ax_ep, self.ax_ee):
            ax.set_facecolor(bg)
            ax.axis("off")
            ax.set_xlim(-0.5, GRID - 0.5)
            ax.set_ylim(-0.5, GRID - 0.5)
        self.fig.text(0.25, 0.97, "ep", ha="center", va="top",
                      color=EP_COLOR, fontsize=8, fontweight="bold")
        self.fig.text(0.75, 0.97, "ee", ha="center", va="top",
                      color=EE_COLOR, fontsize=8, fontweight="bold")
        self._title = self.fig.text(0.5, 0.01, "", ha="center", va="bottom",
                                    color="#aaaaaa", fontsize=7)

        super().__init__(self.fig)
        if parent is not None:
            self.setParent(parent)
        self.setFixedSize(self.FIXED_W, self.FIXED_H)
        self.setStyleSheet("background: rgba(20,20,20,200); border: 1px solid #555;")

        self._cmap = cm.RdYlGn
        self._norm = mcolors.Normalize(vmin=CORR_CAP_LO, vmax=CORR_CAP_HI)
        self._rects: Dict[str, Dict[Tuple[int, int], mpatches.Rectangle]] = {
            "ep": {}, "ee": {}
        }
        for kind, ax in (("ep", self.ax_ep), ("ee", self.ax_ee)):
            for r in range(GRID):
                for c in range(GRID):
                    rect = mpatches.Rectangle(
                        (c - 0.48, (GRID - 1 - r) - 0.48), 0.96, 0.96,
                        linewidth=0.5, edgecolor="#666666", facecolor="#333333",
                        transform=ax.transData,
                    )
                    ax.add_patch(rect)
                    self._rects[kind][(r, c)] = rect
            ax.set_xlim(-0.5, GRID - 0.5)
            ax.set_ylim(-0.5, GRID - 0.5)

        self.cellClicked = None  # callable(r, c, kind)
        self.mpl_connect("button_press_event", self._on_click)

    def update_data(self, ed: Optional[ModuleEnergyData],
                    mod_name: str = "",
                    selected: Optional[Tuple[int, int, str]] = None) -> None:
        for kind in ("ep", "ee"):
            for (r, c), rect in self._rects[kind].items():
                if ed is not None:
                    corr = ed.cells[(r, c)].corr(kind)
                    rect.set_facecolor(self._cmap(self._norm(corr)))
                else:
                    rect.set_facecolor("#333333")
                sel = (selected is not None
                       and selected[:2] == (r, c)
                       and selected[2] == kind)
                rect.set_edgecolor("#ffffff" if sel else "#555555")
                rect.set_linewidth(2.0 if sel else 0.5)
        self._title.set_text(mod_name)
        self.fig.canvas.draw_idle()

    def _on_click(self, event):
        if self.cellClicked is None:
            return
        for ax, kind in ((self.ax_ep, "ep"), (self.ax_ee, "ee")):
            if event.inaxes is ax:
                c = int(round(event.xdata))
                r = GRID - 1 - int(round(event.ydata))
                if 0 <= r < GRID and 0 <= c < GRID:
                    self.cellClicked(r, c, kind)
                return


# ── Map container with floating overlay ──────────────────────────────────────

class MapContainer(QWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._map  = HyCalMapWidget(self, enable_zoom_pan=True, min_size=(400, 400))
        self._mini = MiniGridCanvas(self)

    def main_map(self) -> HyCalMapWidget:
        return self._map

    def mini(self) -> MiniGridCanvas:
        return self._mini

    def resizeEvent(self, event):
        super().resizeEvent(event)
        w, h = self.width(), self.height()
        self._map.setGeometry(0, 0, w, h)
        mw, mh = MiniGridCanvas.FIXED_W, MiniGridCanvas.FIXED_H
        self._mini.setGeometry(w - mw - 4, 4, mw, mh)
        self._mini.raise_()


# ── 5×5 × 2 histogram canvas ──────────────────────────────────────────────────

class GridHistCanvas(FigureCanvas):
    """Two 5×5 blocks: ep (top) with ep fit+correction, ee (bottom) with ee fit+correction."""

    def __init__(self):
        fc = getattr(THEME, "CANVAS", "#1e1e1e")
        self.fig = Figure(figsize=(10, 14), facecolor=fc)
        height_ratios = [1] * GRID + [0.22] + [1] * GRID
        gs = self.fig.add_gridspec(
            GRID * 2 + 1, GRID,
            height_ratios=height_ratios,
            hspace=0.58, wspace=0.35,
            left=0.055, right=0.985, top=0.975, bottom=0.025,
        )
        self.axes_ep: List[List] = []
        self.axes_ee: List[List] = []
        for r in range(GRID):
            row_ep, row_ee = [], []
            for c in range(GRID):
                ax_ep = self.fig.add_subplot(gs[r, c])
                ax_ee = self.fig.add_subplot(gs[GRID + 1 + r, c])
                self._style_ax(ax_ep)
                self._style_ax(ax_ee)
                row_ep.append(ax_ep)
                row_ee.append(ax_ee)
            self.axes_ep.append(row_ep)
            self.axes_ee.append(row_ee)

        self.fig.text(0.005, 0.988, "e-p", color=EP_COLOR,
                      fontsize=10, fontweight="bold", va="top")
        self.fig.text(0.005, 0.49,  "e-e", color=EE_COLOR,
                      fontsize=10, fontweight="bold", va="top")

        super().__init__(self.fig)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
        self._spans: List = []
        self._selected: Optional[Tuple[int, int, str]] = None
        self.cellClicked = None   # callable(r, c, kind)
        self.spanChanged = None   # callable(r, c, kind, lo, hi)
        self.mpl_connect("button_press_event", self._on_click)

    def _style_ax(self, ax):
        bg = getattr(THEME, "AXES", "#1e1e1e")
        ax.set_facecolor(bg)
        ax.tick_params(colors="#888888", labelsize=5, length=2, pad=1)
        for sp in ax.spines.values():
            sp.set_color("#444444")
            sp.set_linewidth(0.5)

    def clear_all(self):
        self._spans.clear()
        for block in (self.axes_ep, self.axes_ee):
            for row in block:
                for ax in row:
                    ax.cla(); self._style_ax(ax)
        self.fig.canvas.draw_idle()

    def draw_module(self, mod: ModuleData, subdir: str,
                    selected: Optional[Tuple[int, int, str]] = None):
        self._spans.clear()
        ed = mod.energies[subdir]
        for r in range(GRID):
            for c in range(GRID):
                cell = ed.cells[(r, c)]
                ax_ep = self.axes_ep[r][c]; ax_ee = self.axes_ee[r][c]
                ax_ep.cla(); ax_ee.cla()
                self._style_ax(ax_ep); self._style_ax(ax_ee)
                self._draw_cell(ax_ep, cell, "ep", selected == (r, c, "ep"))
                self._draw_cell(ax_ee, cell, "ee", selected == (r, c, "ee"))
                if c == 0:
                    ax_ep.set_ylabel(f"r{r}", fontsize=4,
                                     color="#777777", labelpad=1)
                    ax_ee.set_ylabel(f"r{r}", fontsize=4,
                                     color="#777777", labelpad=1)
                if r == 0:
                    ax_ep.set_title(f"c{c}", fontsize=4, color="#777777", pad=2)

                sp_ep = SpanSelector(
                    ax_ep,
                    lambda lo, hi, _r=r, _c=c: self._on_span("ep", _r, _c, lo, hi),
                    "horizontal", useblit=True,
                    props={"alpha": 0.2, "facecolor": EP_COLOR},
                    interactive=False, drag_from_anywhere=True,
                )
                sp_ee = SpanSelector(
                    ax_ee,
                    lambda lo, hi, _r=r, _c=c: self._on_span("ee", _r, _c, lo, hi),
                    "horizontal", useblit=True,
                    props={"alpha": 0.2, "facecolor": EE_COLOR},
                    interactive=False, drag_from_anywhere=True,
                )
                self._spans.extend([sp_ep, sp_ee])

        self._selected = selected
        self.fig.canvas.draw_idle()

    def _draw_cell(self, ax, cell: CellData, kind: str, selected: bool):
        color = EP_COLOR if kind == "ep" else EE_COLOR
        for sp in ax.spines.values():
            sp.set_color("#ffffff" if selected else "#444444")
            sp.set_linewidth(1.5 if selected else 0.5)

        if cell.hist_counts is None or len(cell.hist_counts) == 0:
            ax.text(0.5, 0.5, "—", ha="center", va="center",
                    color="#555555", fontsize=6, transform=ax.transAxes)
            return

        centers = 0.5 * (cell.hist_edges[:-1] + cell.hist_edges[1:])
        ax.step(centers, cell.hist_counts, where="mid",
                color="#aaaaaa", linewidth=0.5)
        ax.fill_between(centers, cell.hist_counts, step="mid",
                        color="#444444", alpha=0.5)

        Eexp = cell.expected_ep if kind == "ep" else cell.expected_ee
        if Eexp > 0.:
            ax.axvline(Eexp, color="#44cc44", linewidth=0.7,
                       linestyle="--", alpha=0.8)

        fit_valid = cell.fit_valid_ep if kind == "ep" else cell.fit_valid_ee
        fit_mean  = cell.fit_mean_ep  if kind == "ep" else cell.fit_mean_ee
        fit_sigma = cell.fit_sigma_ep if kind == "ep" else cell.fit_sigma_ee
        fit_amp   = cell.fit_amp_ep   if kind == "ep" else cell.fit_amp_ee
        span_lo   = cell.span_lo_ep   if kind == "ep" else cell.span_lo_ee
        span_hi   = cell.span_hi_ep   if kind == "ep" else cell.span_hi_ee
        corr      = cell.corr_ep      if kind == "ep" else cell.corr_ee

        if fit_valid and fit_mean > 0. and fit_sigma > 0.:
            xc = np.linspace(
                max(centers[0], fit_mean - 3. * fit_sigma),
                min(centers[-1], fit_mean + 3. * fit_sigma), 200)
            ax.plot(xc, _gauss(xc, fit_amp, fit_mean, fit_sigma),
                    color=color, linewidth=0.9)

        if span_lo is not None and span_hi is not None:
            ax.axvspan(span_lo, span_hi, alpha=0.15, color=color)

        ax.text(0.97, 0.96, f"{corr:.4f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=4.5, color=color)
        if cell.entries > 0:
            ax.text(0.03, 0.96, f"{cell.entries}", transform=ax.transAxes,
                    ha="left", va="top", fontsize=4, color="#888888")

        if Eexp > 0.:
            sig_e = cell.sigma_exp(kind)
            ax.set_xlim(Eexp - 5. * sig_e, Eexp + 5. * sig_e)
        ax.set_yticks([])
        # Adjust y-range to the max count within the visible x window
        xlim = ax.get_xlim()
        mask = (centers >= xlim[0]) & (centers <= xlim[1])
        ymax = float(cell.hist_counts[mask].max()) if mask.any() else 1.
        ax.set_ylim(0., ymax * 1.15)

    def redraw_cell(self, r: int, c: int, kind: str, cell: CellData):
        block = self.axes_ep if kind == "ep" else self.axes_ee
        ax = block[r][c]
        ax.cla(); self._style_ax(ax)
        self._draw_cell(ax, cell, kind, self._selected == (r, c, kind))
        self.fig.canvas.draw_idle()

    def set_selected(self, r: int, c: int, kind: str):
        if self._selected is not None:
            pr, pc, pk = self._selected
            ax = (self.axes_ep if pk == "ep" else self.axes_ee)[pr][pc]
            for s in ax.spines.values():
                s.set_color("#444444"); s.set_linewidth(0.5)
        self._selected = (r, c, kind)
        ax = (self.axes_ep if kind == "ep" else self.axes_ee)[r][c]
        for s in ax.spines.values():
            s.set_color("#ffffff"); s.set_linewidth(1.5)
        self.fig.canvas.draw_idle()

    def update_corr_text(self, r: int, c: int, kind: str, cell: CellData):
        block = self.axes_ep if kind == "ep" else self.axes_ee
        ax = block[r][c]
        color = EP_COLOR if kind == "ep" else EE_COLOR
        for t in [x for x in ax.texts
                  if x.get_ha() == "right"
                  and abs(x.get_position()[1] - 0.96) < 0.05]:
            t.remove()
        corr = cell.corr_ep if kind == "ep" else cell.corr_ee
        ax.text(0.97, 0.96, f"{corr:.4f}", transform=ax.transAxes,
                ha="right", va="top", fontsize=4.5, color=color)
        self.fig.canvas.draw_idle()

    def _on_click(self, event):
        if event.inaxes is None or self.cellClicked is None:
            return
        for r in range(GRID):
            for c in range(GRID):
                for kind, block in (("ep", self.axes_ep), ("ee", self.axes_ee)):
                    if event.inaxes is block[r][c]:
                        self.cellClicked(r, c, kind); return

    def _on_span(self, kind: str, r: int, c: int, lo: float, hi: float):
        if abs(hi - lo) < 1.0:
            return
        if self.spanChanged is not None:
            self.spanChanged(r, c, kind, lo, hi)


# ── Main window ───────────────────────────────────────────────────────────────

class EnergyCorrViewerWindow(QMainWindow):

    def __init__(self, root_path: Optional[Path] = None,
                 json_path: Optional[Path] = None):
        super().__init__()
        self.setWindowTitle("HyCal Energy Correction Viewer")
        self.resize(1680, 980)
        self._data:      Dict[str, ModuleData] = {}
        self._cur:       Optional[ModuleData]  = None
        self._subdir:    str = BEAM_CONFIGS[0][1]
        self._json_path: Optional[Path] = json_path
        self._sel:       Optional[Tuple[int, int, str]] = None
        self._w_modules: List[Module] = []
        self._build_ui()
        if root_path is not None:
            self._load_root(root_path)
        if json_path is not None and self._data:
            apply_json_corrections(self._data, json_path)
            self._refresh_map()

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root_lay = QVBoxLayout(central)
        root_lay.setContentsMargins(4, 4, 4, 4)
        root_lay.setSpacing(4)

        top = QHBoxLayout()
        btn_root = QPushButton("Open ROOT…")
        btn_root.clicked.connect(self._open_root)
        top.addWidget(btn_root)
        btn_json = QPushButton("Open JSON…")
        btn_json.clicked.connect(self._open_json)
        top.addWidget(btn_json)
        self._file_lbl = QLabel("No file loaded")
        self._file_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM};")
        top.addWidget(self._file_lbl, stretch=1)
        top.addWidget(QLabel("Beam:"))
        self._beam_btns: List[QPushButton] = []
        for label, subdir, _ in BEAM_CONFIGS:
            btn = QPushButton(label); btn.setCheckable(True)
            btn.clicked.connect(lambda _, s=subdir: self._set_beam(s))
            top.addWidget(btn); self._beam_btns.append(btn)
        self._beam_btns[0].setChecked(True)
        top.addSpacing(10)
        top.addWidget(QLabel("Map color:"))
        self._color_cb = QComboBox()
        self._color_cb.addItems([
            "ep mean corr", "ee mean corr",
            "ep max |dev|", "ee max |dev|",
        ])
        self._color_cb.currentIndexChanged.connect(self._refresh_map)
        top.addWidget(self._color_cb)
        top.addSpacing(10)
        self._export_btn = QPushButton("Export JSON")
        self._export_btn.setEnabled(False)
        self._export_btn.clicked.connect(self._export_json)
        top.addWidget(self._export_btn)
        root_lay.addLayout(top)

        splitter = QSplitter(Qt.Orientation.Horizontal)
        root_lay.addWidget(splitter, stretch=1)

        self._map_container = MapContainer()
        self._map  = self._map_container.main_map()
        self._mini = self._map_container.mini()
        self._map.moduleClicked.connect(self._on_module_click)
        self._mini.cellClicked = self._on_mini_cell_click
        map_wrap = QWidget()
        map_vlay = QVBoxLayout(map_wrap)
        map_vlay.setContentsMargins(0, 0, 0, 0)
        map_vlay.addWidget(self._map_container, stretch=1)
        self._range_ctrl = ColorRangeControl(self._map, orientation="horizontal",
                                             parent=map_wrap)
        map_vlay.addWidget(self._range_ctrl)
        splitter.addWidget(map_wrap)

        right = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(4, 0, 4, 4)
        right_lay.setSpacing(4)
        self._mod_lbl = QLabel("← click a module on the map")
        self._mod_lbl.setFont(QFont("Monospace", 10, QFont.Weight.Bold))
        self._mod_lbl.setStyleSheet(f"color: {THEME.ACCENT};")
        right_lay.addWidget(self._mod_lbl)
        self._hist_canvas = GridHistCanvas()
        self._hist_canvas.cellClicked = self._on_hist_cell_click
        self._hist_canvas.spanChanged = self._on_span_changed
        scroll = QScrollArea()
        scroll.setWidget(self._hist_canvas)
        scroll.setWidgetResizable(False)
        scroll.setHorizontalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        scroll.setVerticalScrollBarPolicy(Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        right_lay.addWidget(scroll, stretch=1)
        right_lay.addWidget(self._build_controls())
        splitter.addWidget(right)
        splitter.setSizes([460, 1110])

        QShortcut(QKeySequence("R"), self).activated.connect(
            lambda: self._reset_sel_to_1("ep"))
        QShortcut(QKeySequence("1"), self).activated.connect(
            lambda: self._reset_sel_to_1("ep"))
        QShortcut(QKeySequence("E"), self).activated.connect(
            lambda: self._reset_sel_to_1("ee"))

    def _build_controls(self) -> QGroupBox:
        box = QGroupBox("Cell corrections")
        lay = QVBoxLayout(box)
        lay.setContentsMargins(8, 6, 8, 6)
        lay.setSpacing(4)
        self._cell_lbl = QLabel("No cell selected")
        self._cell_lbl.setFont(QFont("Monospace", 9))
        lay.addWidget(self._cell_lbl)

        for kind, attr, shortcut in (("ep", "_ep_edit", "R"), ("ee", "_ee_edit", "E")):
            color = EP_COLOR if kind == "ep" else EE_COLOR
            row = QHBoxLayout()
            lbl = QLabel(f"<b style='color:{color}'>{kind}</b> corr:")
            lbl.setTextFormat(Qt.TextFormat.RichText)
            row.addWidget(lbl)
            edit = QLineEdit(); edit.setFixedWidth(78)
            edit.setPlaceholderText("1.0000")
            edit.returnPressed.connect(lambda k=kind: self._apply_edit(k))
            setattr(self, attr, edit)
            row.addWidget(edit)
            btn_apply = QPushButton("Apply"); btn_apply.setFixedWidth(52)
            btn_apply.clicked.connect(lambda _, k=kind: self._apply_edit(k))
            row.addWidget(btn_apply)
            btn_r1 = QPushButton(f"=1  [{shortcut}]"); btn_r1.setFixedWidth(65)
            btn_r1.clicked.connect(lambda _, k=kind: self._reset_sel_to_1(k))
            row.addWidget(btn_r1)
            btn_refit = QPushButton(f"Refit {kind}"); btn_refit.setFixedWidth(65)
            btn_refit.clicked.connect(lambda _, k=kind: self._refit_cell(k))
            row.addWidget(btn_refit)
            row.addStretch()
            lay.addLayout(row)

        bulk = QHBoxLayout()
        for label, fn in (
            ("Refit All (ep+ee)", self._refit_all),
            ("Reset All ep → 1", lambda: self._reset_all("ep")),
            ("Reset All ee → 1", lambda: self._reset_all("ee")),
            ("Clear Span", self._clear_span),
        ):
            b = QPushButton(label); b.clicked.connect(fn); bulk.addWidget(b)
        bulk.addStretch()
        lay.addLayout(bulk)

        self._status_lbl = QLabel("")
        self._status_lbl.setFont(QFont("Monospace", 8))
        self._status_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM};")
        lay.addWidget(self._status_lbl)
        return box

    # ── loading ───────────────────────────────────────────────────────────

    def _open_root(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open energy_corr ROOT file", "", "ROOT files (*.root);;All (*)")
        if path:
            self._load_root(Path(path))

    def _open_json(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open correction JSON", "", "JSON files (*.json);;All (*)")
        if path:
            self._json_path = Path(path)
            if self._data:
                apply_json_corrections(self._data, self._json_path)
                self._refresh_map()
                if self._cur is not None:
                    self._draw_current()
                self._status("Loaded JSON: " + path)

    def _load_root(self, path: Path):
        self._file_lbl.setText(f"Loading {path.name}…")
        QApplication.processEvents()
        self._data = load_root(path)
        self._cur  = None; self._sel = None
        self._file_lbl.setText(str(path))
        self._mod_lbl.setText("← click a module on the map")
        self._export_btn.setEnabled(bool(self._data))
        self._hist_canvas.clear_all()
        self._mini.update_data(None)
        if MODULES_JSON is not None:
            mods = load_modules(MODULES_JSON)
            self._w_modules = [m for m in mods if m.mod_type.lower() == "pbwo4"]
        else:
            self._w_modules = [
                Module(d.name, "PbWO4", d.x, d.y, 20.77, 20.75)
                for d in self._data.values()
            ]
        self._map.set_modules(self._w_modules)
        self._refresh_map()
        self._status(f"Loaded {len(self._data)} modules from {path.name}")

    # ── map ───────────────────────────────────────────────────────────────

    def _set_beam(self, subdir: str):
        self._subdir = subdir
        for i, (_, sd, _) in enumerate(BEAM_CONFIGS):
            self._beam_btns[i].setChecked(sd == subdir)
        self._refresh_map()
        if self._cur is not None:
            self._draw_current()

    def _refresh_map(self):
        if not self._data:
            return
        mode = self._color_cb.currentIndex()
        kind    = "ep" if mode % 2 == 0 else "ee"
        use_dev = mode >= 2
        values: Dict[str, float] = {}
        for name, mod in self._data.items():
            ed = mod.energies[self._subdir]
            values[name] = (ed.max_deviation(kind) if use_dev
                            else ed.mean_correction(kind))
        self._map.set_values(values)
        if use_dev:
            self._map.set_range(0., 0.08)
        else:
            self._map.set_range(CORR_CAP_LO, CORR_CAP_HI)

    def _on_module_click(self, name: str):
        if not name or name not in self._data:
            return
        self._cur = self._data[name]
        self._sel = None
        self._cell_lbl.setText(
            f"Module: {self._cur.name}   θ={self._cur.theta:.2f}°")
        self._ep_edit.clear(); self._ee_edit.clear()
        self._draw_current()

    def _draw_current(self):
        if self._cur is None:
            return
        ed = self._cur.energies[self._subdir]
        label, _, _ = next(b for b in BEAM_CONFIGS if b[1] == self._subdir)
        n_ep = ed.valid_count("ep"); n_ee = ed.valid_count("ee")
        mc_ep = ed.mean_correction("ep"); mc_ee = ed.mean_correction("ee")
        self._mod_lbl.setText(
            f"Module: {self._cur.name}   θ={self._cur.theta:.2f}°   "
            f"Beam: {label}   "
            f"ep fits: {n_ep}/25 mean={mc_ep:.4f}   "
            f"ee fits: {n_ee}/25 mean={mc_ee:.4f}")
        self._hist_canvas.draw_module(self._cur, self._subdir, self._sel)
        self._mini.update_data(ed, self._cur.name, self._sel)

    # ── cell selection ────────────────────────────────────────────────────

    def _on_hist_cell_click(self, r: int, c: int, kind: str):
        self._sel = (r, c, kind)
        if self._cur is None:
            return
        cell = self._cur.energies[self._subdir].cells[(r, c)]
        self._hist_canvas.set_selected(r, c, kind)
        self._mini.update_data(self._cur.energies[self._subdir],
                               self._cur.name, (r, c, kind))
        self._cell_lbl.setText(
            f"Cell ({r},{c})  entries: {cell.entries}  "
            f"ep: E={cell.expected_ep:.1f} fit={cell.fit_mean_ep:.1f}  "
            f"ee: E={cell.expected_ee:.1f} fit={cell.fit_mean_ee:.1f}")
        self._ep_edit.setText(f"{cell.corr_ep:.4f}")
        self._ee_edit.setText(f"{cell.corr_ee:.4f}")

    def _on_mini_cell_click(self, r: int, c: int, kind: str):
        self._on_hist_cell_click(r, c, kind)

    # ── span / refit ──────────────────────────────────────────────────────

    def _on_span_changed(self, r: int, c: int, kind: str, lo: float, hi: float):
        if self._cur is None:
            return
        cell = self._cur.energies[self._subdir].cells[(r, c)]
        if kind == "ep":
            cell.span_lo_ep = lo; cell.span_hi_ep = hi
        else:
            cell.span_lo_ee = lo; cell.span_hi_ee = hi
        if cell.refit(kind):
            self._hist_canvas.redraw_cell(r, c, kind, cell)
            self._mini.update_data(self._cur.energies[self._subdir],
                                   self._cur.name, self._sel)
            self._refresh_map()
            if self._sel and self._sel[:2] == (r, c):
                if kind == "ep": self._ep_edit.setText(f"{cell.corr_ep:.4f}")
                else:            self._ee_edit.setText(f"{cell.corr_ee:.4f}")
            c_val = cell.corr_ep if kind == "ep" else cell.corr_ee
            m_val = cell.fit_mean_ep if kind == "ep" else cell.fit_mean_ee
            self._status(f"Refit ({r},{c}) {kind}: fit={m_val:.1f} corr={c_val:.4f}")
        else:
            self._status(f"Refit ({r},{c}) {kind}: failed")

    def _refit_cell(self, kind: str):
        if self._cur is None or self._sel is None:
            return
        r, c = self._sel[0], self._sel[1]
        cell = self._cur.energies[self._subdir].cells[(r, c)]
        if cell.refit(kind):
            self._hist_canvas.redraw_cell(r, c, kind, cell)
            if kind == "ep": self._ep_edit.setText(f"{cell.corr_ep:.4f}")
            else:            self._ee_edit.setText(f"{cell.corr_ee:.4f}")
            self._mini.update_data(self._cur.energies[self._subdir],
                                   self._cur.name, self._sel)
            self._refresh_map()
            self._status(f"Refit ({r},{c}) {kind} OK")
        else:
            self._status(f"Refit ({r},{c}) {kind}: failed")

    def _refit_all(self):
        if self._cur is None:
            return
        n_ep = n_ee = 0
        for cell in self._cur.energies[self._subdir].cells.values():
            if cell.refit("ep"): n_ep += 1
            if cell.refit("ee"): n_ee += 1
        self._draw_current()
        self._refresh_map()
        self._status(
            f"Refitted ep:{n_ep}/25 ee:{n_ee}/25 "
            f"for {self._cur.name} [{self._subdir}]")

    # ── correction editing ────────────────────────────────────────────────

    def _apply_edit(self, kind: str):
        if self._cur is None or self._sel is None:
            return
        r, c = self._sel[0], self._sel[1]
        edit = self._ep_edit if kind == "ep" else self._ee_edit
        try:
            val = float(edit.text())
        except ValueError:
            self._status("Invalid number"); return
        cell = self._cur.energies[self._subdir].cells[(r, c)]
        if kind == "ep": cell.corr_ep = val
        else:            cell.corr_ee = val
        self._hist_canvas.update_corr_text(r, c, kind, cell)
        self._mini.update_data(self._cur.energies[self._subdir],
                               self._cur.name, self._sel)
        self._refresh_map()
        self._status(f"Set ({r},{c}) {kind} = {val:.4f}")

    def _reset_sel_to_1(self, kind: str):
        if self._cur is None or self._sel is None:
            return
        r, c = self._sel[0], self._sel[1]
        cell = self._cur.energies[self._subdir].cells[(r, c)]
        if kind == "ep": cell.corr_ep = 1.0; self._ep_edit.setText("1.0")
        else:            cell.corr_ee = 1.0; self._ee_edit.setText("1.0")
        self._hist_canvas.update_corr_text(r, c, kind, cell)
        self._mini.update_data(self._cur.energies[self._subdir],
                               self._cur.name, self._sel)
        self._refresh_map()
        self._status(f"Reset ({r},{c}) {kind} → 1.0")

    def _reset_all(self, kind: str):
        if self._cur is None:
            return
        for cell in self._cur.energies[self._subdir].cells.values():
            if kind == "ep": cell.corr_ep = 1.0
            else:            cell.corr_ee = 1.0
        self._draw_current(); self._refresh_map()
        self._status(
            f"Reset all {kind} for {self._cur.name} [{self._subdir}] → 1.0")

    def _clear_span(self):
        if self._cur is None or self._sel is None:
            return
        r, c, kind = self._sel
        cell = self._cur.energies[self._subdir].cells[(r, c)]
        cell.span_lo_ep = None; cell.span_hi_ep = None
        cell.span_lo_ee = None; cell.span_hi_ee = None
        self._hist_canvas.redraw_cell(r, c, "ep", cell)
        self._hist_canvas.redraw_cell(r, c, "ee", cell)
        self._status(f"Cleared spans for ({r},{c})")

    # ── JSON export ───────────────────────────────────────────────────────

    def _export_json(self):
        if not self._data:
            return
        path = self._json_path
        if path is None:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save correction JSON", "corrections.json",
                "JSON files (*.json);;All (*)")
            if not path:
                return
            path = Path(path)
        save_json(self._data, path)
        self._json_path = path
        self._status(f"Exported to {path}")

    def _status(self, msg: str):
        self._status_lbl.setText(msg)


# ── Entry point ───────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(description="HyCal energy-correction viewer")
    ap.add_argument("root_file", nargs="?", type=Path,
                    help="ROOT file from prad2ana_energy_corr")
    ap.add_argument("--json", type=Path, default=None,
                    help="Correction JSON (read + write)")
    ap.add_argument("--theme", choices=available_themes(), default="dark")
    args = ap.parse_args()

    set_theme(args.theme)
    app = QApplication(sys.argv)
    app.setApplicationName("Energy Corr Viewer")
    win = EnergyCorrViewerWindow(root_path=args.root_file, json_path=args.json)
    apply_theme_palette(win)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
