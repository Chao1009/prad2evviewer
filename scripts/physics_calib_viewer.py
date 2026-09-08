#!/usr/bin/env python3
"""Viewer for ``physics_calib.cpp`` output.

Manual fits update ``calib_result_iterN.json`` and the accepted factor is
written back to ``calib_factor_iterN.json``.
"""
from __future__ import annotations

import argparse
import json
import math
import os
import re
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
from PyQt6.QtCore import Qt, QThread, pyqtSignal
from PyQt6.QtGui import QColor, QFont, QPen
from PyQt6.QtWidgets import (
	QApplication, QButtonGroup, QComboBox, QFileDialog, QGroupBox,
	QHBoxLayout, QLabel, QLineEdit, QMainWindow, QMessageBox, QPushButton,
	QRadioButton, QSpinBox, QSplitter, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.colors import LogNorm
from matplotlib.widgets import RectangleSelector, SpanSelector

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

SCRIPT_DIR = Path(__file__).resolve().parent
if str(SCRIPT_DIR) not in sys.path:
	sys.path.insert(0, str(SCRIPT_DIR))

from hycal_geoview import (  # noqa: E402
	HyCalMapWidget, ColorRangeControl, Module, THEME, apply_theme_palette,
	available_themes, load_modules, set_theme,
)

DB_PATH = SCRIPT_DIR.parent / "database" / "hycal_map.json"


@dataclass
class Result:
	module_id: int
	peak: float
	expected_peak: float
	sigma: float
	chi2: float
	ratio: float
	old_factor: float
	new_factor: float
	fit_good: bool
	is_dead: bool = False
	is_dead_neighbor: bool = False


@dataclass
class Iteration:
	run: str
	number: int
	directory: Path
	factors: List[dict] = field(default_factory=list)
	factor_by_id: Dict[int, dict] = field(default_factory=dict)
	results: Dict[int, Result] = field(default_factory=dict)
	result_rows: List[dict] = field(default_factory=list)
	histograms: Dict[str, Tuple[np.ndarray, np.ndarray]] = field(default_factory=dict)
	root_path: Optional[Path] = None

	@property
	def factor_path(self) -> Path:
		return self.directory / f"calib_factor_iter{self.number}.json"

	@property
	def result_path(self) -> Path:
		return self.directory / f"calib_result_iter{self.number}.json"

def module_id_from_name(name: str) -> Optional[int]:
	match = re.fullmatch(r"W(\d+)", name)
	return 1000 + int(match.group(1)) if match else None


def module_name(module_id: int) -> str:
	return f"W{module_id - 1000}"


def _read_json(path: Path, default):
	try:
		with path.open() as stream:
			return json.load(stream)
	except (OSError, json.JSONDecodeError):
		return default


def scan_iterations(base: Path) -> Dict[str, Dict[int, Iteration]]:
	found: Dict[str, Dict[int, Iteration]] = {}
	if not base.is_dir():
		return found
	for run_dir in sorted(p for p in base.iterdir() if p.is_dir()):
		run_items: Dict[int, Iteration] = {}
		for result_path in sorted(run_dir.glob("calib_result_iter*.json")):
			if result_path.name.endswith(".manual.json"):
				continue
			match = re.fullmatch(r"calib_result_iter(\d+)\.json", result_path.name)
			if not match:
				continue
			number = int(match.group(1))
			item = Iteration(run_dir.name, number, run_dir)
			item.root_path = run_dir / f"calib_result_iter{number}.root"
			item.factors = _read_json(item.factor_path, [])
			item.factor_by_id = {
				module_id_from_name(entry.get("name", "")): entry
				for entry in item.factors
				if module_id_from_name(entry.get("name", "")) is not None
			}
			item.result_rows = _read_json(result_path, [])
			for raw in item.result_rows:
				try:
					module_id = int(raw["module_id"])
					item.results[module_id] = Result(
						module_id, float(raw.get("peak", 0)),
						float(raw.get("expected_peak", 0)),
						float(raw.get("sigma", 0)),
						float(raw.get("chi2/ndf", 0)),
						float(raw.get("ratio", 0)),
						float(raw.get("old_factor", 0)),
						float(raw.get("new_factor", 0)),
						bool(raw.get("fit_good", False)),
						bool(raw.get("is_dead", False)),
						bool(raw.get("is_deadNeighbor", False)),
					)
				except (KeyError, TypeError, ValueError):
					continue
			run_items[number] = item
		if run_items:
			found[run_dir.name] = run_items
	return found


def load_histogram(root, key: str):
	try:
		converted = root[key].to_numpy(flow=False)
		if len(converted) == 2:
			values, edges = converted
			return np.asarray(values, dtype=float), np.asarray(edges, dtype=float)
		values, xedges, yedges = converted
		return np.asarray(values, dtype=float), (
			np.asarray(xedges, dtype=float), np.asarray(yedges, dtype=float))
	except (KeyError, OSError, ValueError, TypeError):
		return None


def load_root_data(item: Iteration, hist_mode: str) -> None:
	if not HAS_UPROOT or item.root_path is None or not item.root_path.is_file():
		return
	try:
		with uproot.open(item.root_path) as root:
			prefix = "modules_5by5/h1_E_mod_" if hist_mode == "5by5" else "modules_island/h1_E_mod_"
			suffix = "_merged" if hist_mode == "5by5" else "_island_merged"
			for module_id in range(1001, 2157):
				data = load_histogram(root, f"{prefix}{module_id}{suffix}")
				if data is not None and np.any(data[0] > 0):
					item.histograms[module_name(module_id)] = data
			for key in ("h2_energy_theta_merged", "hit_pos_merged",
						"h_E_1cl_merged", "h_center_energy_fraction",
						"h_center_energy", "h_fit_peak_energy",
						"h_fit_peak_ratio", "h_fit_peak_chi2ndf",
						"h_fit_peak_sigma"):
				data = load_histogram(root, key)
				if data is not None:
					item.histograms[key] = data
	except Exception:
		return


def gaussian(x, amplitude, mean, sigma):
	return amplitude * np.exp(-0.5 * ((x - mean) / sigma) ** 2)


def fit_histogram(counts, edges, xmin=None, xmax=None, expected=0.0):
	centers = 0.5 * (edges[:-1] + edges[1:])
	mask = counts > 0
	if xmin is not None:
		mask &= centers >= xmin
	if xmax is not None:
		mask &= centers <= xmax
	if mask.sum() < 4:
		raise ValueError("fit range contains fewer than four non-empty bins")
	x = centers[mask]
	y = counts[mask].astype(float)
	peak_guess = expected if expected > 0 else float(x[np.argmax(y)])
	if expected > 0:
		nearby = (x >= expected * 0.8) & (x <= expected * 1.2)
		if np.any(nearby):
			peak_guess = float(x[nearby][np.argmax(y[nearby])])
	peak_guess = float(np.clip(peak_guess, x.min(), x.max()))
	half = max(float(edges[1] - edges[0]) * 2, abs(peak_guess) * 0.03)
	if xmin is None and xmax is None:
		boundary = y.max() * 0.4
		peak_index = int(np.argmax(y))
		left = peak_index
		right = peak_index
		while left > 0 and y[left - 1] >= boundary:
			left -= 1
		while right + 1 < len(y) and y[right + 1] >= boundary:
			right += 1
		x = x[left:right + 1]
		y = y[left:right + 1]
		half = max((x[-1] - x[0]) / 2, abs(peak_guess) * 0.015)
	if not HAS_SCIPY:
		weights = y / y.sum()
		mean = float((x * weights).sum())
		sigma = float(np.sqrt(max(((x - mean) ** 2 * weights).sum(), 1e-12)))
		return mean, sigma, 0.0, float(y.max())
	bandwidth = float(edges[1] - edges[0])
	sigma0 = max(half / 2, bandwidth)
	low = [0.0, float(x.min()), bandwidth]
	high = [float(y.max()) * 10.0, float(x.max()), max(half * 4, bandwidth * 2)]
	params, _ = curve_fit(gaussian, x, y,
						  p0=[float(y.max()), peak_guess, sigma0],
						  bounds=(low, high), maxfev=10000)
	residual = y - gaussian(x, *params)
	chi2 = float((residual ** 2 / np.maximum(y, 1.0)).sum()) / max(len(y) - 3, 1)
	amplitude, mean, sigma = map(float, params)
	return mean, sigma, chi2, amplitude


def damped_ratio(expected: float, peak: float) -> float:
	if expected <= 0 or peak <= 0:
		raise ValueError("expected peak and fitted peak must be positive")
	return min(2.0, max(0.5, 1.0 + 0.7 * (expected / peak - 1.0)))


def atomic_json_write(path: Path, value) -> None:
	path.parent.mkdir(parents=True, exist_ok=True)
	fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
	try:
		with os.fdopen(fd, "w") as stream:
			json.dump(value, stream, indent=2)
			stream.write("\n")
		os.replace(temp_name, path)
	except Exception:
		try:
			os.unlink(temp_name)
		except OSError:
			pass
		raise


def atomic_json_write_many(items: List[Tuple[Path, object]]) -> None:
	"""Prepare several JSON files before replacing any destination."""
	temporary: List[Tuple[str, Path]] = []
	try:
		for path, value in items:
			path.parent.mkdir(parents=True, exist_ok=True)
			fd, temp_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
			with os.fdopen(fd, "w") as stream:
				json.dump(value, stream, indent=2)
				stream.write("\n")
			temporary.append((temp_name, path))
		for temp_name, path in temporary:
			os.replace(temp_name, path)
	except Exception:
		for temp_name, _ in temporary:
			try:
				os.unlink(temp_name)
			except OSError:
				pass
		raise


class PhysicsMap(HyCalMapWidget):
	selectionChanged = pyqtSignal(set)

	def __init__(self, parent=None):
		super().__init__(parent, enable_zoom_pan=True, min_size=(500, 500))
		self.selected = set()
		self.multi_select = False
		self.marked = set()
		self.preview = set()

	def set_multi_select(self, enabled):
		self.multi_select = bool(enabled)
		if not enabled:
			self.selected.clear()
			self.selectionChanged.emit(set())
		self.update()

	def clear_selection(self):
		self.selected.clear()
		self.selectionChanged.emit(set())
		self.update()

	def set_preview_modules(self, names):
		self.preview = set(names)
		self.update()

	def _handle_click(self, pos):
		if self._check_inline_range_edit_click(pos):
			return
		name = self._hit(pos)
		if self.multi_select and name:
			if name in self.selected:
				self.selected.remove(name)
			else:
				self.selected.add(name)
			self.selectionChanged.emit(set(self.selected))
			self.update()
		else:
			self.moduleClicked.emit(name or "")

	def _paint_overlays(self, painter, width, height):
		super()._paint_overlays(painter, width, height)
		painter.setBrush(Qt.BrushStyle.NoBrush)
		for name, color, size in ((self.marked, THEME.DANGER, 1.8),
								  (self.selected, THEME.ACCENT, 2.2)):
			painter.setPen(QPen(QColor(color), size))
			for module_name in name:
				rect = self._rects.get(module_name)
				if rect is None:
					continue
				painter.drawEllipse(rect.center(), min(rect.width(), rect.height()) * 0.3,
									min(rect.width(), rect.height()) * 0.3)
		painter.setBrush(Qt.BrushStyle.NoBrush)
		painter.setPen(QPen(QColor(THEME.WARN), 2.2, Qt.PenStyle.DashLine))
		for module_name in self.preview:
			rect = self._rects.get(module_name)
			if rect is not None:
				painter.drawEllipse(rect.center(), min(rect.width(), rect.height()) * 0.42,
									min(rect.width(), rect.height()) * 0.42)

	def _tooltip_text(self, name):
		value = self._values.get(name)
		return f"{name}: {self._fmt_value(value)}" if value is not None else f"{name}: no data"


class Canvas(FigureCanvas):
	def __init__(self, parent=None):
		self.figure = Figure(facecolor=THEME.CANVAS, tight_layout=True)
		self.ax = self.figure.add_subplot(111)
		self.ax.set_facecolor(THEME.CANVAS)
		self._colorbar = None
		self._zoom_selector = None
		self._auto_xlim = None
		self._auto_ylim = None
		super().__init__(self.figure)
		self.setParent(parent)

	def style(self, title, xlabel="", ylabel="", dark=True):
		text = THEME.TEXT if dark else "#202124"
		text_dim = THEME.TEXT_DIM if dark else "#4a4a4a"
		self.figure.set_facecolor(THEME.CANVAS if dark else "#ffffff")
		self.ax.set_facecolor(THEME.CANVAS if dark else "#ffffff")
		self.ax.set_title(title, color=text)
		self.ax.set_xlabel(xlabel, color=text_dim)
		self.ax.set_ylabel(ylabel, color=text_dim)
		self.ax.tick_params(colors=text_dim)
		for spine in self.ax.spines.values():
			spine.set_color(THEME.BORDER if dark else "#777777")

	def clear_colorbar(self):
		if self._colorbar is not None:
			self._colorbar.remove()
			self._colorbar = None

	def enable_drag_zoom(self):
		if self._zoom_selector is not None:
			self._zoom_selector.set_active(False)
		self._zoom_selector = RectangleSelector(
			self.ax, self._drag_zoom_selected, useblit=False,
			button=[1, 3], minspanx=5, minspany=5,
			spancoords="pixels", interactive=False,
			props={"facecolor": "#5b9bd5", "edgecolor": "#1f4e79", "alpha": 0.25})

	def set_auto_limits(self, xlim, ylim=None):
		self._auto_xlim = tuple(xlim) if xlim else None
		self._auto_ylim = tuple(ylim) if ylim else None

	def reset_zoom(self):
		if self._auto_xlim:
			self.ax.set_xlim(*self._auto_xlim)
		if self._auto_ylim:
			self.ax.set_ylim(*self._auto_ylim)
		self.draw_idle()

	def _drag_zoom_selected(self, start, end):
		if start.xdata is None or start.ydata is None or end.xdata is None or end.ydata is None:
			return
		x0, x1 = sorted((start.xdata, end.xdata))
		y0, y1 = sorted((start.ydata, end.ydata))
		if x1 <= x0 or y1 <= y0:
			return
		button = getattr(start, "button", 1)
		if button == 3:
			current_x = self.ax.get_xlim()
			current_y = self.ax.get_ylim()
			cx, cy = (x0 + x1) / 2.0, (y0 + y1) / 2.0
			xspan = current_x[1] - current_x[0]
			yspan = current_y[1] - current_y[0]
			factor_x = max(xspan / (x1 - x0), 1.0)
			factor_y = max(yspan / (y1 - y0), 1.0)
			new_x = xspan * factor_x
			new_y = yspan * factor_y
			self.ax.set_xlim(cx - new_x / 2.0, cx + new_x / 2.0)
			self.ax.set_ylim(cy - new_y / 2.0, cy + new_y / 2.0)
		else:
			self.ax.set_xlim(x0, x1)
			self.ax.set_ylim(y0, y1)
		self.draw_idle()


class RootWorker(QThread):
	ready = pyqtSignal(object)

	def __init__(self, item, hist_mode, parent=None):
		super().__init__(parent)
		self.item = item
		self.hist_mode = hist_mode

	def run(self):
		load_root_data(self.item, self.hist_mode)
		self.ready.emit(self.item)


class Viewer(QMainWindow):
	MODES = ("Has data", "fit_good", "module flag", "chi2/ndf", "sigma", "delta E", "|delta E/expected|", "ratio")

	def __init__(self, initial_dir=None, hist_mode="5by5"):
		super().__init__()
		self.base = None
		self.scan = {}
		self.current: Optional[Iteration] = None
		self.hist_mode = hist_mode
		self.current_module = ""
		self.worker = None
		self.span = None
		self.modules = []
		self.geometry = {}
		self._manual_overlays = {}
		self._rebin = 1
		self._rebinned_modules = set()
		self._module_rebin = {}
		self._build_ui()
		self.setWindowTitle(f"Physics Calibration Viewer ({hist_mode})")
		self.resize(1500, 940)
		self._style()
		if DB_PATH.is_file():
			self.modules = load_modules(DB_PATH)
			self._load_geometry()
			self.map.set_modules(self.modules)
		if initial_dir:
			self.load_directory(Path(initial_dir))

	def _load_geometry(self):
		data = _read_json(DB_PATH, [])
		self.geometry = {e.get("n"): e.get("geo", {}) for e in data
						 if e.get("n") and e.get("geo")}

	def _style(self):
		t = THEME
		self.setStyleSheet(f"""
		QMainWindow, QWidget {{ background:{t.BG}; color:{t.TEXT}; }}
		QComboBox, QLineEdit, QSpinBox {{ background:{t.PANEL}; color:{t.TEXT};
			border:1px solid {t.BORDER}; padding:3px; }}
		QSpinBox::up-button, QSpinBox::down-button {{
			subcontrol-origin: border; width:22px; background:{t.BUTTON_HOVER};
			border-left:1px solid {t.BORDER};
		}}
		QSpinBox::up-button {{ subcontrol-position: top right; border-bottom:1px solid {t.BORDER}; }}
		QSpinBox::down-button {{ subcontrol-position: bottom right; }}
		QSpinBox::up-button:hover, QSpinBox::down-button:hover {{ background:{t.ACCENT_STRONG}; }}
		QSpinBox::up-arrow {{
			width:0px; height:0px; border-left:5px solid transparent;
			border-right:5px solid transparent; border-bottom:6px solid white;
		}}
		QSpinBox::down-arrow {{
			width:0px; height:0px; border-left:5px solid transparent;
			border-right:5px solid transparent; border-top:6px solid white;
		}}
		QPushButton {{ background:{t.BUTTON}; color:{t.TEXT}; border:1px solid {t.BORDER};
			padding:4px 9px; border-radius:3px; }}
		QPushButton:hover {{ background:{t.BUTTON_HOVER}; }}
		QGroupBox {{ color:{t.TEXT_DIM}; border:1px solid {t.BORDER}; margin-top:6px; }}
		QGroupBox::title {{ subcontrol-origin:margin; left:8px; }}
		""")

	def _build_ui(self):
		root = QWidget()
		self.setCentralWidget(root)
		outer = QVBoxLayout(root)
		toolbar = QHBoxLayout()
		open_button = QPushButton("Open Physics_calib")
		open_button.clicked.connect(self._browse)
		toolbar.addWidget(open_button)
		self.dir_label = QLabel("(no directory)")
		toolbar.addWidget(self.dir_label, 1)
		toolbar.addWidget(QLabel("Run:"))
		self.run_box = QComboBox()
		self.run_box.currentTextChanged.connect(self._run_changed)
		toolbar.addWidget(self.run_box)
		toolbar.addWidget(QLabel("Iteration:"))
		self.iter_box = QComboBox()
		self.iter_box.currentTextChanged.connect(self._iter_changed)
		toolbar.addWidget(self.iter_box)
		toolbar.addWidget(QLabel("Histogram:"))
		hist_box = QComboBox()
		hist_box.addItems(("5by5", "island"))
		hist_box.setCurrentText(self.hist_mode)
		hist_box.currentTextChanged.connect(self._hist_mode_changed)
		toolbar.addWidget(hist_box)
		outer.addLayout(toolbar)

		mode_bar = QHBoxLayout()
		self.mode_group = QButtonGroup(self)
		for index, mode in enumerate(self.MODES):
			button = QRadioButton(mode)
			button.clicked.connect(lambda _, value=mode: self._set_mode(value))
			self.mode_group.addButton(button, index)
			mode_bar.addWidget(button)
			if index == 0:
				button.setChecked(True)
		mode_bar.addStretch()
		outer.addLayout(mode_bar)

		self.map = PhysicsMap()
		self.map.moduleClicked.connect(self._module_clicked)
		self.map.moduleHovered.connect(self._module_hovered)
		self.map.selectionChanged.connect(self._selection_changed)
		self.range_control = ColorRangeControl(self.map, auto_fit="minmax", include_log=True)
		outer.addWidget(self.range_control)

		splitter = QSplitter(Qt.Orientation.Horizontal)
		outer.addWidget(splitter, 1)
		left = QGroupBox("HyCal map")
		left_layout = QVBoxLayout(left)
		left_layout.addWidget(self.map)
		batch = QHBoxLayout()
		self.multi_button = QPushButton("Multi-select")
		self.multi_button.setCheckable(True)
		self.multi_button.toggled.connect(self.map.set_multi_select)
		batch.addWidget(self.multi_button)
		self.selection_label = QLabel("0 selected")
		batch.addWidget(self.selection_label)
		batch.addWidget(QLabel("Factor:"))
		self.factor_edit = QLineEdit()
		self.factor_edit.setFixedWidth(82)
		batch.addWidget(self.factor_edit)
		apply_button = QPushButton("Apply selected")
		apply_button.clicked.connect(self._apply_selected)
		batch.addWidget(apply_button)
		restore_button = QPushButton("Restore selected old_factor")
		restore_button.clicked.connect(self._restore_selected)
		batch.addWidget(restore_button)
		batch.addWidget(QLabel("Outer layers:"))
		self.layers = QSpinBox()
		self.layers.setRange(1, 12)
		self.layers.setValue(1)
		batch.addWidget(self.layers)
		self.outer_shape = QComboBox()
		self.outer_shape.addItems(("Square", "Circle"))
		self.outer_shape.setToolTip("Select square row/column rings or circular radial rings")
		batch.addWidget(self.outer_shape)
		outer_restore = QPushButton("Restore outer W")
		outer_restore.clicked.connect(self._restore_outer)
		batch.addWidget(outer_restore)
		self.outer_rebin_spin = QSpinBox()
		self.outer_rebin_spin.setRange(1, 100)
		self.outer_rebin_spin.setValue(1)
		self.outer_rebin_spin.setPrefix("Rebin ")
		self.outer_rebin_spin.setToolTip("Combine this many adjacent bins for outer-layer histograms")
		batch.addWidget(self.outer_rebin_spin)
		preview_outer = QPushButton("Preview outer")
		preview_outer.clicked.connect(self._preview_outer)
		preview_outer.setToolTip("Show the selected outer modules on the HyCal map")
		batch.addWidget(preview_outer)
		rebin_outer = QPushButton("Rebin outer")
		rebin_outer.clicked.connect(self._rebin_outer)
		batch.addWidget(rebin_outer)
		left_layout.addLayout(batch)
		splitter.addWidget(left)

		right = QSplitter(Qt.Orientation.Vertical)
		splitter.addWidget(right)
		detail = QGroupBox("Module detail and manual fit")
		detail_layout = QVBoxLayout(detail)
		self.info = QLabel("Click a module")
		self.info.setWordWrap(True)
		self.info.setMinimumHeight(82)
		self.info.setFont(QFont("Monospace", 12, QFont.Weight.Bold))
		self.info.setStyleSheet(f"color: {THEME.TEXT}; padding: 8px; background: {THEME.PANEL};")
		detail_layout.addWidget(self.info)
		self.canvas = Canvas()
		detail_layout.addWidget(self.canvas, 1)
		hist_rebin_row = QHBoxLayout()
		hist_rebin_row.addWidget(QLabel("Histogram rebin:"))
		self.rebin_spin = QSpinBox()
		self.rebin_spin.setRange(1, 100)
		self.rebin_spin.setValue(1)
		self.rebin_spin.setToolTip("Fit and display the current histogram after combining adjacent bins")
		self.rebin_spin.valueChanged.connect(self._hist_rebin_changed)
		hist_rebin_row.addWidget(self.rebin_spin)
		restore_module = QPushButton("Restore old factor")
		restore_module.clicked.connect(self._restore_current_module)
		restore_module.setToolTip("Write this module's producer old_factor back to the factor JSON")
		hist_rebin_row.addWidget(restore_module)
		hist_rebin_row.addStretch()
		detail_layout.addLayout(hist_rebin_row)
		fit_bar = QHBoxLayout()
		fit_bar.addWidget(QLabel("Range:"))
		self.xmin = QLineEdit()
		self.xmax = QLineEdit()
		self.xmin.setPlaceholderText("xmin")
		self.xmax.setPlaceholderText("xmax")
		self.xmin.setFixedWidth(75)
		self.xmax.setFixedWidth(75)
		fit_bar.addWidget(self.xmin)
		fit_bar.addWidget(self.xmax)
		fit_button = QPushButton("Run fit")
		fit_button.clicked.connect(self._run_fit)
		fit_bar.addWidget(fit_button)
		apply_fit = QPushButton("Apply fit + save")
		apply_fit.clicked.connect(self._apply_fit)
		fit_bar.addWidget(apply_fit)
		detail_layout.addLayout(fit_bar)
		self.fit_status = QLabel("")
		self.fit_status.setWordWrap(True)
		detail_layout.addWidget(self.fit_status)
		right.addWidget(detail)

		stats = QGroupBox("Global diagnostics")
		stats_layout = QVBoxLayout(stats)
		self.stats_choice = "Energy vs theta"
		self.stats_buttons = QButtonGroup(self)
		stats_bar = QHBoxLayout()
		for index, choice in enumerate(("Energy vs theta", "Hit position", "One cluster energy",
										"Peak ratio", "Fit chi2/ndf", "Fit sigma")):
			button = QPushButton(choice)
			button.setCheckable(True)
			button.clicked.connect(lambda _, value=choice: self._set_global_choice(value))
			self.stats_buttons.addButton(button, index)
			stats_bar.addWidget(button)
			if index == 0:
				button.setChecked(True)
		reset_zoom = QPushButton("Reset zoom")
		reset_zoom.clicked.connect(self._reset_global_zoom)
		reset_zoom.setToolTip("Restore the automatic chart limits")
		stats_bar.addWidget(reset_zoom)
		stats_layout.addLayout(stats_bar)
		self.stats_canvas = Canvas()
		stats_layout.addWidget(self.stats_canvas)
		right.addWidget(stats)
		splitter.setSizes([650, 780])
		right.setSizes([600, 300])

	def _set_global_choice(self, choice):
		self.stats_choice = choice
		self._draw_global()

	def _reset_global_zoom(self):
		self.stats_canvas.reset_zoom()

	def _browse(self):
		selected = QFileDialog.getExistingDirectory(self, "Physics_calib directory")
		if selected:
			self.load_directory(Path(selected))

	def load_directory(self, path):
		self.base = path
		self.scan = scan_iterations(path)
		self.dir_label.setText(str(path))
		self.run_box.blockSignals(True)
		self.run_box.clear()
		self.run_box.addItems(sorted(self.scan))
		self.run_box.blockSignals(False)
		if self.run_box.count():
			self._run_changed(self.run_box.currentText())

	def _run_changed(self, run):
		self.iter_box.blockSignals(True)
		self.iter_box.clear()
		if run in self.scan:
			self.iter_box.addItems(str(number) for number in sorted(self.scan[run]))
		self.iter_box.blockSignals(False)
		if self.iter_box.count():
			self._iter_changed(self.iter_box.currentText())

	def _iter_changed(self, text):
		if not text or self.run_box.currentText() not in self.scan:
			return
		self.current = self.scan[self.run_box.currentText()][int(text)]
		self._manual_overlays = {}
		self._rebin = 1
		self._rebinned_modules = set()
		self._module_rebin = {}
		self.rebin_spin.setValue(1)
		self.map.clear_selection()
		self.map.set_preview_modules(set())
		self.map.marked = set()
		self._refresh_map(auto_range=True)
		self._draw_global()
		self._start_root_load()

	def _hist_mode_changed(self, mode):
		self.hist_mode = mode
		if self.current:
			self._manual_overlays = {}
			self._rebin = 1
			self._rebinned_modules = set()
			self._module_rebin = {}
			self.rebin_spin.setValue(1)
			self.current.histograms.clear()
			self._start_root_load()

	def _start_root_load(self):
		if self.current is None:
			return
		if self.worker and self.worker.isRunning():
			self.worker.quit()
		self.worker = RootWorker(self.current, self.hist_mode, self)
		self.worker.ready.connect(self._root_ready)
		self.worker.start()

	def _root_ready(self, item):
		if item is self.current:
			self._draw_global()
			if self.current_module:
				self._show_module(self.current_module)

	def _set_mode(self, mode):
		self._map_mode = mode
		self._refresh_map(auto_range=True)

	def _refresh_map(self, auto_range=False):
		if self.current is None:
			return
		values = {}
		mode = getattr(self, "_map_mode", self.MODES[0])
		for module_id, result in self.current.results.items():
			name = module_name(module_id)
			if mode == "Has data":
				values[name] = 1.0
			elif mode == "fit_good":
				values[name] = 1.0 if result.fit_good else 0.0
			elif mode == "module flag":
				values[name] = 2.0 if result.is_dead else (1.0 if result.is_dead_neighbor else 0.0)
			elif mode == "chi2/ndf":
				values[name] = result.chi2
			elif mode == "sigma":
				values[name] = result.sigma
			elif mode == "delta E":
				values[name] = result.peak - result.expected_peak
			elif mode == "|delta E/expected|":
				values[name] = (
					abs(result.peak - result.expected_peak) / result.expected_peak
					if result.expected_peak > 0 else 0.0
				)
			elif mode == "ratio":
				values[name] = result.ratio
		self.map.set_values(values)
		self.map.set_map_label(mode) if hasattr(self.map, "set_map_label") else None
		if auto_range:
			self.range_control.notify_values_changed(values)
			if mode in ("Has data", "fit_good"):
				self.range_control.set_range(0.0, 1.0)
			elif mode == "module flag":
				self.range_control.set_range(0.0, 2.0)
			else:
				self.range_control.controller.auto_fit(values)

	def _hist_rebin_changed(self, factor):
		self._rebin = max(1, int(factor))
		if self.current_module:
			self._module_rebin[self.current_module] = self._rebin
			self._draw_module(self.current_module)

	def _module_hovered(self, name):
		if self.current and (module_id_from_name(name) in self.current.results):
			result = self.current.results[module_id_from_name(name)]
			self.statusBar().showMessage(
				f"{name}: peak={result.peak:.2f}, expected={result.expected_peak:.2f}, "
				f"ratio={result.ratio:.5f}")

	def _module_clicked(self, name):
		if name:
			self.current_module = name
			self.rebin_spin.blockSignals(True)
			self.rebin_spin.setValue(self._module_rebin.get(name, 1))
			self.rebin_spin.blockSignals(False)
			self._rebin = self.rebin_spin.value()
			self._show_module(name)

	def _show_module(self, name):
		if self.current is None:
			return
		module_id = module_id_from_name(name)
		result = self.current.results.get(module_id)
		if result is None:
			self.info.setText(f"{name}: no result JSON entry")
		else:
			current_factor = self.current.factor_by_id.get(module_id, {}).get("factor", 0.0)
			self.info.setText(
				f"{name} (id {module_id})\n"
				f"peak={result.peak:.2f}  expected={result.expected_peak:.2f}  "
				f"sigma={result.sigma:.2f}  chi2/ndf={result.chi2:.4f}\n"
				f"ratio={result.ratio:.6f}  old_factor={result.old_factor:.7f}  "
				f"new_factor={result.new_factor:.7f}  current factor={float(current_factor):.7f}")
		self._draw_module(name)

	def _draw_module(self, name, overlay=None):
		self.canvas.ax.clear()
		data = self.current.histograms.get(name) if self.current else None
		if data is None:
			self.canvas.style(f"{name} energy histogram", "Energy (MeV)", "Counts")
			self.canvas.draw_idle()
			return
		counts, edges = self._display_histogram(name)
		centers = 0.5 * (edges[:-1] + edges[1:])
		self.canvas.ax.bar(centers, counts, width=np.diff(edges), color=THEME.ACCENT,
						   alpha=0.72, linewidth=0)
		result = self.current.results.get(module_id_from_name(name))
		if overlay is None:
			overlay = self._manual_overlays.get(name)
		fit_peak = overlay[0] if overlay else (result.peak if result else None)
		if result and result.expected_peak > 0:
			self.canvas.ax.axvline(
				result.expected_peak, color=THEME.SUCCESS, linewidth=2.0,
				linestyle=":", label=f"Expected peak = {result.expected_peak:.2f}")
		if result and result.peak > 0 and result.sigma > 0:
			x = np.linspace(result.peak - 4 * result.sigma,
							result.peak + 4 * result.sigma, 500)
			amp = float(counts.max())
			self.canvas.ax.plot(x, amp * np.exp(-0.5 * ((x - result.peak) / result.sigma) ** 2),
								"--", color=THEME.TEXT_DIM, label="producer fit")
		if overlay:
			peak, sigma, amp = overlay
			x = np.linspace(peak - 4 * sigma, peak + 4 * sigma, 500)
			self.canvas.ax.plot(x, amp * np.exp(-0.5 * ((x - peak) / sigma) ** 2),
								color=THEME.WARN, linewidth=2, label="manual fit")
		if fit_peak is not None and fit_peak > 0:
			self.canvas.ax.axvline(
				fit_peak, color=THEME.WARN, linewidth=1.8,
				linestyle="--", label=f"Fit peak = {fit_peak:.2f}")
		if result or overlay:
			self.canvas.ax.legend(
				loc="upper left", facecolor=THEME.PANEL,
				labelcolor=THEME.TEXT, framealpha=0.9)
		x_candidates = []
		nonzero = np.flatnonzero(counts > 0)
		if nonzero.size:
			x_candidates.append((float(edges[nonzero[0]]),
								  float(edges[nonzero[-1] + 1])))
		if result and result.peak > 0 and result.sigma > 0:
			x_candidates.append((result.peak - 4.0 * result.sigma,
								 result.peak + 4.0 * result.sigma))
		if overlay:
			x_candidates.append((overlay[0] - 4.0 * overlay[1],
								 overlay[0] + 4.0 * overlay[1]))
		if not x_candidates:
			x_candidates.append((float(edges[0]), float(edges[-1])))
		xmin = min(pair[0] for pair in x_candidates)
		xmax = max(pair[1] for pair in x_candidates)
		if xmax > xmin:
			pad = max((xmax - xmin) * 0.06, float(edges[1] - edges[0]))
			self.canvas.ax.set_xlim(xmin - pad, xmax + pad)
		self.canvas.style(f"{name} energy histogram ({self.hist_mode})", "Energy (MeV)", "Counts")
		self.canvas.draw_idle()
		if self.span:
			self.span.disconnect_events()
		self.span = SpanSelector(self.canvas.ax, self._span_selected, "horizontal",
								 useblit=True, props={"alpha": 0.25, "facecolor": THEME.WARN})

	def _span_selected(self, xmin, xmax):
		if xmax > xmin:
			self.xmin.setText(f"{xmin:.2f}")
			self.xmax.setText(f"{xmax:.2f}")
			self._run_fit()

	def _run_fit(self):
		if not self.current or not self.current_module:
			return
		module_id = module_id_from_name(self.current_module)
		result = self.current.results.get(module_id)
		data = self._display_histogram(self.current_module)
		if result is None or data is None:
			self.fit_status.setText("No result or histogram for this module")
			return
		try:
			xmin = float(self.xmin.text()) if self.xmin.text().strip() else None
			xmax = float(self.xmax.text()) if self.xmax.text().strip() else None
			peak, sigma, chi2, amplitude = fit_histogram(
				data[0], data[1], xmin, xmax, result.expected_peak)
			ratio = damped_ratio(result.expected_peak, peak)
			new_factor = result.old_factor * ratio
			self._fit_value = (peak, sigma, chi2, amplitude, ratio, new_factor)
			self._manual_overlays[self.current_module] = (peak, sigma, amplitude)
			self.fit_status.setText(
				f"manual peak={peak:.3f}, sigma={sigma:.3f}, chi2/ndf={chi2:.5f}; "
				f"ratio={ratio:.7f}, new_factor={new_factor:.8f}")
			self._draw_module(self.current_module, (peak, sigma, amplitude))
		except (ValueError, RuntimeError) as exc:
			self.fit_status.setText(f"Fit failed: {exc}")

	def _apply_fit(self):
		if not hasattr(self, "_fit_value") or not self.current or not self.current_module:
			return
		module_id = module_id_from_name(self.current_module)
		result = self.current.results.get(module_id)
		if result is None:
			return
		peak, sigma, chi2, _amplitude, ratio, new_factor = self._fit_value
		result_rows = [dict(row) for row in self.current.result_rows]
		result_row = next((row for row in result_rows
						  if int(row.get("module_id", -1)) == module_id), None)
		if result_row is None:
			self.fit_status.setText("Result JSON has no entry for this module")
			return
		result_row.update({
			"peak": peak,
			"sigma": sigma,
			"chi2/ndf": chi2,
			"ratio": ratio,
			"new_factor": new_factor,
			"fit_good": bool(peak > 0 and sigma > 0 and chi2 < 1.8),
		})
		factors = [dict(entry) for entry in self.current.factors]
		updated = False
		for entry in factors:
			if entry.get("name") == module_name(module_id):
				entry["factor"] = new_factor
				updated = True
				break
		if not updated:
			self.fit_status.setText("Factor JSON has no entry for this module")
			return
		try:
			atomic_json_write_many([
				(self.current.factor_path, factors),
				(self.current.result_path, result_rows),
			])
		except OSError as exc:
			self.fit_status.setText(f"Save failed; memory unchanged: {exc}")
			return
		self.current.factors = factors
		self.current.factor_by_id = {module_id_from_name(e.get("name", "")): e
									 for e in factors if module_id_from_name(e.get("name", ""))}
		self.current.result_rows = result_rows
		self.current.results[module_id] = Result(
			module_id, peak, result.expected_peak, sigma, chi2, ratio,
			result.old_factor, new_factor,
			bool(peak > 0 and sigma > 0 and chi2 < 1.8),
			result.is_dead, result.is_dead_neighbor)
		self.map.marked.add(self.current_module)
		self._refresh_map()
		self.fit_status.setText("Applied: result and factor JSON saved")
		self._show_module(self.current_module)

	def _selection_changed(self, selected):
		self.selection_label.setText(f"{len(selected)} selected")

	def _finish_multi_selection(self, names):
		self.map.marked.update(names)
		self.map.clear_selection()
		self.multi_button.setChecked(False)

	def _rebin_histogram(self, counts, edges, factor):
		factor = max(1, int(factor))
		if factor == 1 or len(counts) < factor:
			return counts, edges
		usable = (len(counts) // factor) * factor
		rebinned = counts[:usable].reshape(-1, factor).sum(axis=1)
		new_edges = edges[:usable + 1:factor]
		if len(new_edges) != len(rebinned) + 1:
			new_edges = np.r_[new_edges, edges[usable]]
		return rebinned, new_edges

	def _display_histogram(self, name):
		data = self.current.histograms.get(name) if self.current else None
		if data is None:
			return None
		counts, edges = data
		factor = self._module_rebin.get(name, self._rebin)
		return self._rebin_histogram(counts, edges, factor)

	def _outer_module_names(self):
		w_modules = [(m, self.geometry.get(m.name, {})) for m in self.modules
					 if m.name.startswith("W") and self.geometry.get(m.name)]
		if not w_modules:
			return set()
		layers = self.layers.value()
		if self.outer_shape.currentText() == "Square":
			rows = [int(g["row"]) for _, g in w_modules if "row" in g]
			cols = [int(g["col"]) for _, g in w_modules if "col" in g]
			if not rows or not cols:
				return set()
			max_row, max_col = max(rows), max(cols)
			return {m.name for m, geo in w_modules
					if "row" in geo and "col" in geo and
					min(int(geo["row"]) - 1, max_row - int(geo["row"]),
						int(geo["col"]) - 1, max_col - int(geo["col"])) + 1 <= layers}
		xs = np.asarray([float(m.x) for m, _ in w_modules])
		ys = np.asarray([float(m.y) for m, _ in w_modules])
		cx, cy = float((xs.min() + xs.max()) / 2.0), float((ys.min() + ys.max()) / 2.0)
		pitch = float(np.median([max(m.sx, m.sy) for m, _ in w_modules]))
		distances = np.hypot(xs - cx, ys - cy)
		max_radius = float(distances.max())
		# One radial shell is approximately one crystal pitch wide.  The small
		# tolerance keeps corner modules in the requested outer shell.
		cut = max_radius - layers * pitch * 1.15
		return {m.name for m, distance in zip((m for m, _ in w_modules), distances)
				if distance >= cut}

	def _apply_factor_map(self, names, factor_by_name):
		if not self.current or not names:
			return
		factors = [dict(entry) for entry in self.current.factors]
		changed = set()
		for entry in factors:
			name = entry.get("name")
			if name in names and name in factor_by_name:
				entry["factor"] = factor_by_name[name]
				changed.add(name)
		if not changed:
			return
		try:
			atomic_json_write(self.current.factor_path, factors)
		except OSError as exc:
			self.statusBar().showMessage(f"Save failed: {exc}", 6000)
			return
		self.current.factors = factors
		self.current.factor_by_id = {module_id_from_name(e.get("name", "")): e
									 for e in factors if module_id_from_name(e.get("name", ""))}
		self.map.marked.update(changed)
		self._finish_multi_selection(changed)
		self.statusBar().showMessage(f"Saved {len(changed)} factor(s)", 4000)

	def _apply_selected(self):
		try:
			factor = float(self.factor_edit.text())
			if factor <= 0:
				raise ValueError
		except ValueError:
			self.statusBar().showMessage("Factor must be positive", 4000)
			return
		self._apply_factor_map(self.map.selected, {name: factor for name in self.map.selected})

	def _restore_selected(self):
		if not self.current:
			return
		values = {}
		for name in self.map.selected:
			result = self.current.results.get(module_id_from_name(name))
			if result and result.old_factor > 0:
				values[name] = result.old_factor
		self._apply_factor_map(values.keys(), values)

	def _restore_current_module(self):
		if not self.current or not self.current_module:
			self.statusBar().showMessage("Select a module first", 4000)
			return
		module_id = module_id_from_name(self.current_module)
		result = self.current.results.get(module_id)
		if result is None or result.old_factor <= 0:
			self.statusBar().showMessage(
				f"{self.current_module} has no valid old_factor in result JSON", 4000)
			return
		self._apply_factor_map(
			{self.current_module}, {self.current_module: result.old_factor})
		self._show_module(self.current_module)

	def _restore_outer(self):
		if not self.current:
			return
		names = self._outer_module_names()
		values = {}
		for name in names:
			result = self.current.results.get(module_id_from_name(name))
			if result and result.old_factor > 0:
				values[name] = result.old_factor
		self._apply_factor_map(values.keys(), values)

	def _preview_outer(self):
		if not self.current:
			return
		if self.map.preview:
			self.map.set_preview_modules(set())
			self.statusBar().showMessage("Outer-layer preview cleared", 4000)
			return
		names = self._outer_module_names()
		self.map.set_preview_modules(names)
		self.statusBar().showMessage(
			f"Previewing {len(names)} {self.outer_shape.currentText().lower()} outer-layer modules", 4000)

	def _rebin_outer(self):
		if not self.current:
			return
		self._rebin = self.outer_rebin_spin.value()
		outer_names = self._outer_module_names()
		self._rebinned_modules = {
			name for name in outer_names if name in self.current.histograms
		}
		for name in self._rebinned_modules:
			self._module_rebin[name] = self._rebin
		if self.current_module in outer_names:
			self._draw_module(self.current_module)
		self.statusBar().showMessage(
			f"Rebinned {len(outer_names)} outer module histogram(s) by {self._rebin}", 4000)

	def _draw_global(self):
		if not self.current:
			return
		if self.stats_canvas._zoom_selector is not None:
			self.stats_canvas._zoom_selector.set_active(False)
		self.stats_canvas.clear_colorbar()
		self.stats_canvas.ax.clear()
		choice = self.stats_choice
		keys = {"Energy vs theta": "h2_energy_theta_merged",
				"Hit position": "hit_pos_merged",
				"One cluster energy": "h_E_1cl_merged",
				"Peak ratio": "h_fit_peak_ratio",
				"Fit chi2/ndf": "h_fit_peak_chi2ndf",
				"Fit sigma": "h_fit_peak_sigma"}
		data = self.current.histograms.get(keys[choice])
		if data is not None:
			counts, edges = data
			if counts.ndim == 1:
				centers = 0.5 * (edges[:-1] + edges[1:])
				self.stats_canvas.ax.bar(centers, counts, width=np.diff(edges),
										 color=THEME.ACCENT)
				nonzero = np.flatnonzero(counts > 0)
				if nonzero.size:
					xmin = float(edges[nonzero[0]])
					xmax = float(edges[nonzero[-1] + 1])
					pad = max((xmax - xmin) * 0.06,
							  float(edges[1] - edges[0]))
					self.stats_canvas.ax.set_xlim(xmin - pad, xmax + pad)
			else:
				xedges, yedges = edges
				positive = counts[counts > 0]
				if positive.size:
					image = self.stats_canvas.ax.pcolormesh(
						xedges, yedges, counts.T, shading="auto", cmap="viridis",
						norm=LogNorm(vmin=float(positive.min()),
									 vmax=float(positive.max())))
					self.stats_canvas._colorbar = self.stats_canvas.figure.colorbar(
						image, ax=self.stats_canvas.ax, pad=0.02, label="Counts")
					self.stats_canvas._colorbar.ax.tick_params(colors="#4a4a4a")
					self.stats_canvas._colorbar.ax.yaxis.label.set_color("#202124")
			self.stats_canvas.style(choice, dark=False)
			if counts.ndim == 1:
				self.stats_canvas.set_auto_limits(
					self.stats_canvas.ax.get_xlim(), self.stats_canvas.ax.get_ylim())
			else:
				self.stats_canvas.set_auto_limits(
					(float(edges[0][0]), float(edges[0][-1])),
					(float(edges[1][0]), float(edges[1][-1])))
			self.stats_canvas.enable_drag_zoom()
		self.stats_canvas.draw_idle()


def main():
	parser = argparse.ArgumentParser(description="physics_calib output viewer")
	parser.add_argument("calib_dir", nargs="?", type=Path,
						help="Physics_calib directory")
	parser.add_argument("--hist-mode", choices=("5by5", "island"), default="5by5")
	parser.add_argument("--theme", choices=available_themes(), default="dark")
	args = parser.parse_args()
	set_theme(args.theme)
	app = QApplication(sys.argv)
	win = Viewer(args.calib_dir, args.hist_mode)
	apply_theme_palette(win)
	win.show()
	sys.exit(app.exec())


if __name__ == "__main__":
	main()
