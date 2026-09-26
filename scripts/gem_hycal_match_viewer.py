#!/usr/bin/env python3
"""
gem_hycal_match_viewer.py — interactive PyQt6 GEM↔HyCal matching viewer.

Open an EVIO file event-by-event and inspect the GEM↔HyCal coincidence:
  * front view (X-Y) — HyCal geometry with cluster centroids, GEM hits
    projected through the target onto the HyCal plane (one mark per
    detector, color-coded), and a dashed line from each HC cluster to its
    matched GEM hit;
  * side view (Z-Y) — target / GEM planes / HyCal face with hit markers
    and the HyCal→GEM line for each matched HC cluster;
  * match table — one row per (HC cluster × GEM detector) with residual
    and sigma_total;
  * toolbar — First/Prev/Next/Goto/Last for navigation, plus a
    "Next matched" search controlled by two thresholds (N hits per
    detector, K detectors with ≥N hits).  The search is a foreground scan
    with a cancellable progress dialog.
  * show/hide — checkboxes per detector + HyCal cluster overlay.

Usage:
    python scripts/gem_hycal_match_viewer.py [file.evio.00000] [-r RUN]
        [--db DIR] [--theme THEME]

Configuration is read from $PRAD2_DATABASE_DIR when set, else from the
--db directory (default <repo>/database):
  * daq_config.json            (DAQ + raw decoding)
  * reconstruction_config.json (runinfo pointer + matching constants)
"""

from __future__ import annotations

import argparse
import os
import sys
import time
from pathlib import Path
from typing import List, Optional, Tuple

_SCRIPT_DIR = Path(__file__).resolve().parent
_REPO_DIR = _SCRIPT_DIR.parent

from prad2_env import database_dir, import_prad2py
prad2py, PRAD2PY_ERROR = import_prad2py()
HAVE_PRAD2PY = prad2py is not None

# analysis/pyscripts/_common.py — parametric matching helpers shared with
# the offline TSV/CSV writer (gem_hycal_matching.py).
_ANA_PY = _REPO_DIR / "analysis" / "pyscripts"
if _ANA_PY.is_dir() and str(_ANA_PY) not in sys.path:
    sys.path.insert(0, str(_ANA_PY))

try:
    import _common as C  # type: ignore
    HAVE_COMMON = True
except Exception as _exc:
    C = None  # type: ignore
    HAVE_COMMON = False
    PRAD2PY_ERROR = (PRAD2PY_ERROR + "\n" if PRAD2PY_ERROR else "") + \
                    f"_common import: {type(_exc).__name__}: {_exc}"

from PyQt6.QtCore import Qt, QPointF
from PyQt6.QtGui import QAction, QBrush, QColor, QFont, QKeySequence, QPainter, QPen
from PyQt6.QtWidgets import (
    QApplication, QCheckBox, QDoubleSpinBox, QFileDialog,
    QGroupBox, QHBoxLayout, QLabel, QMainWindow, QMessageBox, QProgressDialog,
    QPushButton, QSizePolicy, QSpinBox, QSplitter, QStatusBar, QTableWidget,
    QTableWidgetItem, QToolBar, QVBoxLayout, QWidget,
)

import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.figure import Figure

# hycal_geoview is a sibling — provides the HyCal front-view widget and the
# theme system reused across event viewers.
from hycal_geoview import (
    HyCalMapWidget, apply_theme_palette, available_themes,
    load_modules as load_geo_modules, set_theme,
)
from evio_io import EvioCursor, iter_physics_records, open_evio

# Per-detector colour palette — same as the web viewer (resources/gem.js).
GEM_COLORS = [
    QColor("#1f77b4"),  # GEM0 — blue
    QColor("#ff7f0e"),  # GEM1 — orange
    QColor("#2ca02c"),  # GEM2 — green
    QColor("#d62728"),  # GEM3 — red
]
GEM_RGB = [c.getRgbF()[:3] for c in GEM_COLORS]   # matplotlib colours
GEM_NAMES = ["GEM0", "GEM1", "GEM2", "GEM3"]


# ---- Data structures -------------------------------------------------------

class HCCluster:
    __slots__ = ("idx", "x", "y", "energy", "lab_x", "lab_y", "lab_z")

    def __init__(self, idx, h, lab):
        self.idx = idx
        self.x = float(h.x)
        self.y = float(h.y)
        self.energy = float(h.energy)
        self.lab_x, self.lab_y, self.lab_z = lab


class GEMHit:
    __slots__ = ("det_id", "lab_x", "lab_y", "lab_z")

    def __init__(self, det_id, lab):
        self.det_id = det_id
        self.lab_x, self.lab_y, self.lab_z = lab


class Match:
    __slots__ = ("hc_idx", "det_id", "gem_idx", "proj_x", "proj_y",
                 "residual", "sigma_total")

    def __init__(self, hc_idx, det_id, gem_idx, proj_x, proj_y, residual, sigma_total):
        self.hc_idx = hc_idx
        self.det_id = det_id
        self.gem_idx = gem_idx
        self.proj_x = float(proj_x)
        self.proj_y = float(proj_y)
        self.residual = float(residual)
        self.sigma_total = float(sigma_total)


class EventResult:
    """Everything we need to render one event."""
    def __init__(self):
        self.event_num = 0
        self.trigger_bits = 0
        self.hc: List[HCCluster] = []
        self.gem: List[List[GEMHit]] = [[], [], [], []]
        self.matches: List[Match] = []   # one entry per (hc, det) best match

    def match_pair(self, m: Match) -> Optional[Tuple[HCCluster, GEMHit]]:
        """The (HC cluster, GEM hit) of a match, or None if out of range."""
        if m.hc_idx >= len(self.hc) or m.gem_idx >= len(self.gem[m.det_id]):
            return None
        return self.hc[m.hc_idx], self.gem[m.det_id][m.gem_idx]


# ---- Match counting --------------------------------------------------------

def matches_per_det(matches: List[Match]) -> List[int]:
    """Number of matched hits on each GEM detector."""
    counts = [0, 0, 0, 0]
    for m in matches:
        counts[m.det_id] += 1
    return counts


def event_passes(matches: List[Match], min_hits_per_det: int, min_dets: int) -> bool:
    """Event qualifies if at least `min_dets` GEM detectors each have at
    least `min_hits_per_det` matched hits (counted across all HC clusters)."""
    n_pass = sum(1 for c in matches_per_det(matches) if c >= min_hits_per_det)
    return n_pass >= min_dets


# ---- Reconstruction pipeline -----------------------------------------------

class Pipeline:
    """Wraps `_common.setup_pipeline` so the viewer reconstructs identically
    to the offline TSV/CSV writer.  Owns the matching constants too."""

    def __init__(self, db_dir: Path, run_num: int, evio_path: Path):
        self.match_nsigma = 3.0
        # The helper finds the database via $PRAD2_DATABASE_DIR; point it at
        # db_dir unless the variable is already set.
        os.environ.setdefault("PRAD2_DATABASE_DIR", str(db_dir))
        self._p = C.setup_pipeline(
            evio_path=str(evio_path),
            run_num=run_num,
        )
        self.daq_cfg      = self._p.cfg
        self.gem_sys      = self._p.gem_sys
        self.geo          = self._p.geo

        # Matching config (parametric sigma).
        self.match_abc, self.gem_pos_res, _ = C.load_matching_config(self._p)

    def reconstruct(self, fadc_evt, ssp_evt) -> EventResult:
        ev = EventResult()
        ev.event_num = int(fadc_evt.info.event_number)
        ev.trigger_bits = int(fadc_evt.info.trigger_bits)

        # HyCal: waveform → energy → cluster (same logic as gem_hycal_matching.py).
        hc_raw = C.reconstruct_hycal(self._p, fadc_evt)
        for k, h in enumerate(hc_raw):
            ev.hc.append(HCCluster(k, h, C.hycal_to_lab(self._p, h)))

        # GEM: pedestal + CM + ZS → 1D + 2D
        C.reconstruct_gem(self._p, ssp_evt)
        for d in range(min(4, self.gem_sys.get_n_detectors())):
            xform = self._p.gem_xforms[d]
            for g in self.gem_sys.get_hits(d):
                ev.gem[d].append(GEMHit(d, xform.to_lab(g.x, g.y)))

        ev.matches = self.match(ev.hc, ev.gem)
        return ev

    def match(self, hc: List[HCCluster], gem: List[List[GEMHit]]) -> List[Match]:
        """Closest GEM hit per (HC cluster × GEM detector) pair inside
        `match_nsigma · σ_total` at the GEM plane (see C.best_gem_matches)."""
        hc_lab = [(h.lab_x, h.lab_y, h.lab_z, h.energy) for h in hc]
        gem_lab = [[(g.lab_x, g.lab_y, g.lab_z) for g in gl] for gl in gem]
        return [Match(k, d, gi, px, py, dr, st)
                for k, d, gi, px, py, dr, st, _ in C.best_gem_matches(
                    hc_lab, gem_lab, self.match_abc, self.gem_pos_res,
                    self.match_nsigma)]


# ---- Event + visibility state shared by the front and side views -----------

class _MatchOverlayState:
    """Mixin for FrontView / SideView; the view supplies _refresh()."""

    def _init_overlay_state(self):
        self._evt: Optional[EventResult] = None
        self._show_hc = True
        self._show_gem = [True, True, True, True]
        self._show_matches = True
        self._z_hc = 6225.0   # default; updated when geometry loads
        self._z_gem = [5400.0] * 4

    def set_event(self, evt: Optional[EventResult]):
        self._evt = evt
        self._refresh()

    def set_zs(self, z_hc: float, z_gem: List[float]):
        if z_hc > 0:
            self._z_hc = z_hc
        for i, z in enumerate(z_gem[:4]):
            if z > 0:
                self._z_gem[i] = z

    def set_show_hc(self, on: bool):
        self._show_hc = on; self._refresh()

    def set_show_gem(self, det_id: int, on: bool):
        if 0 <= det_id < 4:
            self._show_gem[det_id] = on; self._refresh()

    def set_show_matches(self, on: bool):
        self._show_matches = on; self._refresh()


# ---- HyCal front view: GEM-projected hits + match lines --------------------

class FrontView(_MatchOverlayState, HyCalMapWidget):
    def __init__(self, parent=None):
        super().__init__(parent)
        self._init_overlay_state()

    def _refresh(self):
        self.update()

    def _to_hc_plane(self, g: GEMHit) -> Tuple[float, float]:
        """GEM hit projected through the target onto the HyCal plane:
        (lab_x, lab_y) · z_hc / z_gem."""
        z_gem = g.lab_z if g.lab_z > 0 else self._z_gem[g.det_id]
        return C.project_to_z(g.lab_x, g.lab_y, z_gem, self._z_hc)[:2]

    def _paint_overlays(self, p: QPainter, w: int, h: int):
        super()._paint_overlays(p, w, h)
        if not self._evt:
            return
        ev = self._evt

        # HyCal cluster centroids (HyCal-local x,y).
        if self._show_hc:
            p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
            for c in ev.hc:
                pt = self.geo_to_canvas(c.x, c.y)
                p.setPen(QPen(QColor("#ffffff"), 2.0))
                p.setBrush(QBrush(QColor(255, 255, 255, 90)))
                p.drawEllipse(pt, 6.0, 6.0)
                p.setPen(QPen(QColor("#000000"), 1.0))
                font = QFont(p.font()); font.setPointSize(8); p.setFont(font)
                p.drawText(pt + QPointF(8, -2), f"HC{c.idx}: {c.energy:.0f} MeV")

        # GEM hits projected through the target onto HyCal-local x,y.
        # HyCal-local equals lab x,y for an untilted HyCal centred at (0,0,z_hc),
        # which is the standard PRad-II geometry; if hycal_x/y or tilts are non-
        # zero the runinfo loader has already absorbed them into lab coords.
        for d in range(4):
            if not self._show_gem[d]:
                continue
            color = GEM_COLORS[d]
            for g in ev.gem[d]:
                qp = self.geo_to_canvas(*self._to_hc_plane(g))
                p.setPen(QPen(color, 1.4))
                p.setBrush(QBrush(QColor(color.red(), color.green(), color.blue(), 130)))
                p.drawEllipse(qp, 4.0, 4.0)

        # Dashed line from each HC cluster to its best-match GEM hit
        # projected onto the HC plane.
        if self._show_matches and ev.matches:
            for m in ev.matches:
                pair = ev.match_pair(m)
                if not self._show_gem[m.det_id] or pair is None:
                    continue
                c, g = pair
                a = self.geo_to_canvas(c.x, c.y)
                b = self.geo_to_canvas(*self._to_hc_plane(g))
                pen = QPen(GEM_COLORS[m.det_id], 1.6)
                pen.setStyle(Qt.PenStyle.DashLine)
                p.setPen(pen)
                p.drawLine(a, b)


# ---- Side view (Z-Y) — matplotlib canvas -----------------------------------

class SideView(_MatchOverlayState, FigureCanvasQTAgg):
    def __init__(self, parent=None):
        self._fig = Figure(figsize=(6, 4), tight_layout=True)
        super().__init__(self._fig)
        if parent is not None:
            self.setParent(parent)
        self._ax = self._fig.add_subplot(111)
        self._ax.set_xlabel("z (mm)")
        self._ax.set_ylabel("y (mm)")
        self._init_overlay_state()
        self._y_size_gem = [600.0] * 4

    def set_geom(self, z_hc: float, z_gem: List[float],
                 y_size_gem: List[float]):
        self.set_zs(z_hc, z_gem)
        for i, y in enumerate(y_size_gem[:4]):
            if y > 0:
                self._y_size_gem[i] = y

    def _refresh(self):
        self.redraw()

    def redraw(self):
        ax = self._ax
        ax.clear()
        ax.set_xlabel("z (mm)")
        ax.set_ylabel("y (mm)")
        ax.grid(True, color="#888", alpha=0.2, linewidth=0.5)

        # Detector frames (dashed) — GEM and HyCal active areas in y.
        for d in range(4):
            yh = self._y_size_gem[d] / 2
            color = GEM_RGB[d]
            ax.plot([self._z_gem[d], self._z_gem[d]], [-yh, yh],
                    "--", color=color, alpha=0.6, linewidth=1.0)
            ax.text(self._z_gem[d], yh + 20, GEM_NAMES[d],
                    color=color, fontsize=8, ha="center")
        ax.axvline(self._z_hc, color="#cccccc", linestyle="--", linewidth=1.0)
        ax.text(self._z_hc, 0, "HyCal", color="#cccccc", fontsize=8,
                ha="center", va="bottom", rotation=90)
        ax.axvline(0.0, color="#888", linestyle=":", linewidth=0.8)
        ax.text(0, 0, "T", color="#888", fontsize=8, ha="right")

        evt = self._evt
        if evt:
            # HyCal cluster markers (z_hc, lab_y)
            if self._show_hc:
                for c in evt.hc:
                    ax.plot([c.lab_z], [c.lab_y], "s", color="white",
                            markersize=8, markeredgecolor="black",
                            markeredgewidth=0.8)
            # GEM hits per detector
            for d in range(4):
                if not self._show_gem[d]:
                    continue
                color = GEM_RGB[d]
                xs = [g.lab_z for g in evt.gem[d]]
                ys = [g.lab_y for g in evt.gem[d]]
                if xs:
                    ax.plot(xs, ys, "o", color=color, markersize=5,
                            markeredgecolor="black", markeredgewidth=0.4,
                            alpha=0.85)
            # Matched HC↔GEM lines
            if self._show_matches:
                for m in evt.matches:
                    pair = evt.match_pair(m)
                    if not self._show_gem[m.det_id] or pair is None:
                        continue
                    c, g = pair
                    ax.plot([g.lab_z, c.lab_z], [g.lab_y, c.lab_y],
                            "-", color=GEM_RGB[m.det_id], linewidth=1.0, alpha=0.7)

        # Y range — pad around the largest detector size.
        y_max = max(max(self._y_size_gem) / 2, 50.0)
        ax.set_ylim(-y_max * 1.1, y_max * 1.1)
        ax.set_xlim(-200, self._z_hc * 1.05)
        self.draw_idle()


# ---- Main window -----------------------------------------------------------

class GemHycalMatchViewer(QMainWindow):
    def __init__(self, evio_path: Optional[Path] = None,
                 db_dir: Optional[Path] = None,
                 run_num: int = -1):
        super().__init__()
        self.setWindowTitle("GEM↔HyCal Matching Viewer")
        self.resize(1500, 900)

        if not HAVE_PRAD2PY or not HAVE_COMMON:
            QMessageBox.critical(self, "prad2py / _common not available",
                                 PRAD2PY_ERROR or "Cannot find prad2py or "
                                 "analysis/pyscripts/_common.py — build the "
                                 "Python bindings and re-run.")
            sys.exit(1)

        self._db_dir = Path(db_dir).resolve() if db_dir else database_dir()
        self._run_num = run_num
        self._pipeline: Optional[Pipeline] = None
        self._cursor: Optional[EvioCursor] = None
        self._physics_index: List[Tuple[int, int]] = []  # (record_idx, sub_idx)
        self._cur_idx = -1
        self._cur_event: Optional[EventResult] = None

        self._build_ui()

        if evio_path is not None:
            self._open_file(Path(evio_path))

    # ---- UI ---------------------------------------------------------------
    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        outer = QVBoxLayout(central); outer.setContentsMargins(4, 4, 4, 4)

        self._build_toolbar()
        self._build_search_row(outer)
        self._build_visibility_row(outer)
        self._build_views(outer)
        self._build_match_table(outer)

        self._status = QStatusBar()
        self.setStatusBar(self._status)
        self._update_status_bar()

        m = self.menuBar().addMenu("&File")
        act_open = QAction("&Open EVIO…", self)
        act_open.setShortcut(QKeySequence.StandardKey.Open)
        act_open.triggered.connect(self._on_open)
        m.addAction(act_open)
        act_quit = QAction("&Quit", self)
        act_quit.setShortcut("Ctrl+Q")
        act_quit.triggered.connect(self.close)
        m.addAction(act_quit)

        apply_theme_palette(self)

    def _build_toolbar(self):
        tb = QToolBar()
        tb.setMovable(False)
        self.addToolBar(tb)

        def addBtn(text, slot, shortcut=None, tip=None):
            act = QAction(text, self)
            if shortcut:
                act.setShortcut(shortcut)
            if tip:
                act.setToolTip(tip)
            act.triggered.connect(slot)
            tb.addAction(act)
            return act

        addBtn("⏮ First", self._first, "Home", "Jump to first physics event")
        addBtn("◀ Prev",  self._prev,  "Left", "Previous physics event")
        tb.addSeparator()
        tb.addWidget(QLabel(" Event #"))
        self._sb_idx = QSpinBox()
        self._sb_idx.setRange(0, 0)
        self._sb_idx.setKeyboardTracking(False)
        self._sb_idx.editingFinished.connect(
            lambda: self._goto(self._sb_idx.value()))
        tb.addWidget(self._sb_idx)
        addBtn("Goto", lambda: self._goto(self._sb_idx.value()))
        tb.addSeparator()
        addBtn("Next ▶", self._next, "Right", "Next physics event")
        addBtn("Last ⏭", self._last, "End", "Jump to last physics event")

    def _build_search_row(self, outer: QVBoxLayout):
        box = QGroupBox("Next-matched search")
        lay = QHBoxLayout(box); lay.setContentsMargins(8, 4, 8, 4)
        lay.addWidget(QLabel("≥"))
        self._sb_N = QSpinBox(); self._sb_N.setRange(1, 99); self._sb_N.setValue(1)
        lay.addWidget(self._sb_N)
        lay.addWidget(QLabel("matched hits per detector,"))
        lay.addWidget(QLabel("≥"))
        self._sb_K = QSpinBox(); self._sb_K.setRange(1, 4); self._sb_K.setValue(2)
        lay.addWidget(self._sb_K)
        lay.addWidget(QLabel("detectors satisfied"))
        btn = QPushButton("Find next ▶▶"); btn.setShortcut("Shift+Right")
        btn.clicked.connect(self._find_next_matched)
        lay.addWidget(btn)
        lay.addStretch(1)

        # nsigma override
        lay.addWidget(QLabel("  nσ:"))
        self._sb_ns = QDoubleSpinBox()
        self._sb_ns.setRange(0.5, 20.0); self._sb_ns.setSingleStep(0.5)
        self._sb_ns.setValue(3.0)
        self._sb_ns.valueChanged.connect(self._on_nsigma_changed)
        lay.addWidget(self._sb_ns)

        outer.addWidget(box)

    def _build_visibility_row(self, outer: QVBoxLayout):
        row = QHBoxLayout()
        row.setContentsMargins(8, 0, 8, 0)
        self._cb_hc = QCheckBox("HyCal"); self._cb_hc.setChecked(True)
        self._cb_hc.toggled.connect(self._on_show_hc)
        row.addWidget(self._cb_hc)
        self._cb_gem: List[QCheckBox] = []
        for d in range(4):
            cb = QCheckBox(GEM_NAMES[d]); cb.setChecked(True)
            cb.setStyleSheet(
                f"color: {GEM_COLORS[d].name()}; font-weight: bold;")
            cb.toggled.connect(lambda on, dd=d: self._on_show_gem(dd, on))
            self._cb_gem.append(cb)
            row.addWidget(cb)
        self._cb_match = QCheckBox("Matches"); self._cb_match.setChecked(True)
        self._cb_match.toggled.connect(self._on_show_matches)
        row.addWidget(self._cb_match)
        row.addStretch(1)
        outer.addLayout(row)

    def _build_views(self, outer: QVBoxLayout):
        split = QSplitter(Qt.Orientation.Horizontal)
        outer.addWidget(split, stretch=1)

        self._front = FrontView()
        self._front.setSizePolicy(QSizePolicy.Policy.Expanding,
                                  QSizePolicy.Policy.Expanding)
        split.addWidget(self._front)

        self._side = SideView()
        split.addWidget(self._side)
        split.setSizes([700, 700])

    def _build_match_table(self, outer: QVBoxLayout):
        self._tbl = QTableWidget(0, 7)
        self._tbl.setHorizontalHeaderLabels(
            ["HC#", "GEM", "GEM x (mm)", "GEM y (mm)",
             "residual (mm)", "σ_total (mm)", "ratio"])
        self._tbl.horizontalHeader().setStretchLastSection(True)
        self._tbl.setEditTriggers(QTableWidget.EditTrigger.NoEditTriggers)
        self._tbl.setMaximumHeight(160)
        outer.addWidget(self._tbl)

    # ---- File handling ----------------------------------------------------
    def _on_open(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open EVIO file", str(_REPO_DIR),
            "EVIO files (*.evio *.evio.*);;All files (*)")
        if path:
            self._open_file(Path(path))

    def _open_file(self, path: Path):
        if not path.is_file():
            QMessageBox.warning(self, "EVIO not found", str(path))
            return
        try:
            self._pipeline = Pipeline(self._db_dir, self._run_num, path)
        except (Exception, SystemExit) as exc:
            QMessageBox.critical(self, "Pipeline init failed",
                                 f"{type(exc).__name__}: {exc}")
            return

        if self._cursor is not None:
            self._cursor.close()
            self._cursor = None
        self._physics_index = []
        self._cur_idx = -1
        self._cur_event = None

        # Index physics events on a channel of their own; navigation then
        # seeks a cursor kept open on the file.
        try:
            self._cursor = EvioCursor(path, self._pipeline.daq_cfg)
            ch, is_ra = open_evio(path, self._pipeline.daq_cfg)
        except RuntimeError as exc:
            QMessageBox.critical(self, "Cannot open EVIO", str(exc))
            return
        self._physics_index = self._index_physics_events(ch, is_ra)
        ch.close()
        if not self._physics_index:
            QMessageBox.information(self, "No physics events",
                                    "Scanned the file but found no physics records.")
            return

        # Push detector geometry into the views.
        if self._pipeline.geo:
            geo = self._pipeline.geo
            self._front.set_zs(geo.hycal_z, geo.gem_z)
            y_size = []
            for d in range(min(4, self._pipeline.gem_sys.get_n_detectors())):
                dets = self._pipeline.gem_sys.get_detectors()
                y_size.append(float(dets[d].plane_y.size))
            while len(y_size) < 4:
                y_size.append(600.0)
            self._side.set_geom(geo.hycal_z, geo.gem_z, y_size)

        # Push HyCal modules into the front view.
        map_path = self._db_dir / "hycal_map.json"
        if map_path.is_file():
            try:
                modules = load_geo_modules(map_path)
                self._front.set_modules(modules)
            except Exception as exc:
                self.statusBar().showMessage(f"hycal_map.json: {exc}")

        self._sb_idx.setRange(0, len(self._physics_index) - 1)
        self.setWindowTitle(f"GEM↔HyCal Matching Viewer — {path.name}")
        self._goto(0)

    def _index_physics_events(self, ch, is_ra: bool) -> List[Tuple[int, int]]:
        idx: List[Tuple[int, int]] = []
        dlg = QProgressDialog("Indexing physics events…", "Cancel", 0, 0, self)
        dlg.setWindowModality(Qt.WindowModality.ApplicationModal)
        dlg.setMinimumDuration(250)
        dlg.show()
        last_ui = time.monotonic()
        for rec in iter_physics_records(ch, is_ra, dlg.wasCanceled):
            for sub in range(ch.get_n_events()):
                idx.append((rec, sub))
            now = time.monotonic()
            if now - last_ui > 0.1:
                dlg.setLabelText(f"Indexing physics events… {len(idx):,} found")
                QApplication.processEvents()
                last_ui = now
        dlg.close()
        return idx

    # ---- Navigation -------------------------------------------------------
    def _first(self): self._goto(0)
    def _prev(self):  self._goto(self._cur_idx - 1)
    def _next(self):  self._goto(self._cur_idx + 1)
    def _last(self):  self._goto(len(self._physics_index) - 1)

    def _goto(self, idx: int):
        if not self._physics_index:
            return
        idx = max(0, min(idx, len(self._physics_index) - 1))
        if idx == self._cur_idx and self._cur_event is not None:
            return
        self._cur_idx = idx
        self._cur_event = self._decode_at(idx)
        self._show_current()

    def _show_current(self):
        """Show _cur_event in the views and match table, and sync the
        event spinbox and status bar to _cur_idx."""
        ev = self._cur_event
        if ev is not None:
            self._front.set_event(ev)
            self._side.set_event(ev)
            self._populate_match_table(ev)
        self._sb_idx.blockSignals(True)
        self._sb_idx.setValue(self._cur_idx); self._sb_idx.blockSignals(False)
        self._update_status_bar()

    def _decode_at(self, idx: int) -> Optional[EventResult]:
        if idx < 0 or idx >= len(self._physics_index):
            return None
        rec, sub = self._physics_index[idx]
        return self._decode_sub(sub) if self._load_record(rec) else None

    def _load_record(self, rec: int) -> bool:
        """Seek the cursor to record `rec` and scan it."""
        try:
            self._cursor.seek(rec)
        except RuntimeError:
            return False
        return bool(self._cursor.ch.scan())

    def _decode_sub(self, sub: int) -> Optional[EventResult]:
        """Reconstruct sub-event `sub` of the loaded record."""
        decoded = self._cursor.ch.decode_event(sub, with_ssp=True)
        if not decoded.get("ok"):
            return None
        return self._pipeline.reconstruct(decoded["event"], decoded["ssp"])

    # ---- Search -----------------------------------------------------------
    def _find_next_matched(self):
        if not self._physics_index:
            return
        N = self._sb_N.value(); K = self._sb_K.value()
        start = self._cur_idx + 1
        end   = len(self._physics_index)
        dlg = QProgressDialog(
            f"Searching for an event with ≥{K} detectors at ≥{N} hit(s)…",
            "Cancel", start, end, self)
        dlg.setWindowModality(Qt.WindowModality.ApplicationModal)
        dlg.setMinimumDuration(250)
        # Each record is loaded and scanned once for all of its sub-events.
        loaded, ok = -1, False
        last_ui = time.monotonic()
        found = -1
        for i in range(start, end):
            rec, sub = self._physics_index[i]
            if rec != loaded:
                loaded, ok = rec, self._load_record(rec)
            evr = self._decode_sub(sub) if ok else None
            if evr is not None and event_passes(evr.matches, N, K):
                found = i
                self._cur_event = evr
                break
            now = time.monotonic()
            if now - last_ui > 0.1:
                dlg.setValue(i + 1); QApplication.processEvents()
                last_ui = now
                if dlg.wasCanceled():
                    break
        dlg.close()
        if found >= 0:
            self._cur_idx = found
            self._show_current()
        else:
            self.statusBar().showMessage(
                "No matching event found before EOF.", 5000)

    # ---- Visibility -------------------------------------------------------
    def _on_show_hc(self, on: bool):
        self._front.set_show_hc(on); self._side.set_show_hc(on)

    def _on_show_gem(self, det_id: int, on: bool):
        self._front.set_show_gem(det_id, on); self._side.set_show_gem(det_id, on)

    def _on_show_matches(self, on: bool):
        self._front.set_show_matches(on); self._side.set_show_matches(on)

    def _on_nsigma_changed(self, v: float):
        if self._pipeline:
            self._pipeline.match_nsigma = float(v)
        # Re-run matching on the current event without re-decoding.
        if self._cur_event:
            self._cur_event.matches = self._pipeline.match(
                self._cur_event.hc, self._cur_event.gem)
            self._show_current()

    # ---- Match table ------------------------------------------------------
    def _populate_match_table(self, evt: EventResult):
        self._tbl.setRowCount(len(evt.matches))
        for r, m in enumerate(evt.matches):
            _, g = evt.match_pair(m)
            ratio = m.residual / m.sigma_total if m.sigma_total > 0 else float("inf")
            cells = [str(m.hc_idx), GEM_NAMES[m.det_id],
                     f"{g.lab_x:.2f}", f"{g.lab_y:.2f}",
                     f"{m.residual:.2f}", f"{m.sigma_total:.2f}",
                     f"{ratio:.2f}σ"]
            for c, txt in enumerate(cells):
                item = QTableWidgetItem(txt)
                if c == 1:
                    item.setForeground(GEM_COLORS[m.det_id])
                self._tbl.setItem(r, c, item)

    # ---- Status -----------------------------------------------------------
    def _update_status_bar(self):
        if not self._physics_index:
            self._status.showMessage("No file loaded.")
            return
        if self._cur_event is None:
            self._status.showMessage(
                f"event {self._cur_idx + 1} / {len(self._physics_index)}  (decoding…)")
            return
        ev = self._cur_event
        per_det = matches_per_det(ev.matches)
        n_dets = sum(1 for c in per_det if c > 0)   # detectors with ≥1 match
        msg = (f"event {self._cur_idx + 1}/{len(self._physics_index)}  "
               f"(#{ev.event_num})  trig=0x{ev.trigger_bits:08X}  "
               f"HC={len(ev.hc)}  matches={len(ev.matches)} on {n_dets} det "
               f"[{per_det[0]},{per_det[1]},{per_det[2]},{per_det[3]}]")
        self._status.showMessage(msg)


def main(argv: Optional[List[str]] = None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    ap.add_argument("evio", nargs="?", help="EVIO file to open on startup")
    ap.add_argument("--db", default=None,
                    help="database directory (default: $PRAD2_DATABASE_DIR or "
                         "<repo>/database)")
    ap.add_argument("-r", "--run", type=int, default=-1,
                    help="run number for runinfo lookup (default: sniff filename)")
    ap.add_argument("--theme", choices=available_themes(), default="dark",
                    help="colour theme")
    args = ap.parse_args(argv)

    set_theme(args.theme)
    app = QApplication(sys.argv)
    w = GemHycalMatchViewer(
        evio_path=Path(args.evio) if args.evio else None,
        db_dir=Path(args.db) if args.db else None,
        run_num=args.run,
    )
    w.show()
    return app.exec()


if __name__ == "__main__":
    sys.exit(main())
