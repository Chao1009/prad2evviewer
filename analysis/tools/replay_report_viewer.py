#!/usr/bin/env python3
"""
replay_report_viewer.py — viewer for replay_filter JSON reports.

Reads the JSON written by `prad2ana_replay_filter -j …` and shows it as a
stack of long line charts with the rejected regions dimmed so the user can
see at a glance which parts of a run survived the slow-control cuts.

Two modes:
  * GUI (default): three linked rows — cut status, livetime + data rate,
    EPICS values — all sharing the x-axis (associated_timestamp by default,
    switchable to associated_evn).  EPICS channels are picked from a
    drop-down checklist.  Pan/zoom via the matplotlib navigation toolbar.
    Buttons to open another report or save the current view as PNG/PDF.
  * CLI (`--cli`): renders one static figure with a row per channel —
    status row on top, then livetime + rate, then one row per EPICS
    channel.  Total figure height scales with the number of EPICS rows so
    the chart stays readable for any cut JSON.  Default x-axis is
    associated_timestamp; pass `--evn` to use associated_evn instead.

Backend choice: matplotlib's Qt6 backend (no new dependency on top of
PyQt6 + matplotlib) so the same plot code feeds both modes.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Optional

import numpy as np


# ============================================================================
# Report loading
# ============================================================================

@dataclass
class ChannelSeries:
    """One channel's per-checkpoint trace, aligned to ReportData.evns."""
    name: str             # report key, e.g. "livetime", "epics:hallb_IPM2C21A_CUR"
    label: str            # display label (strips "epics:" prefix)
    is_livetime: bool
    is_epics: bool
    values: np.ndarray    # float, NaN for missing
    pass_mask: np.ndarray # bool


@dataclass
class ReportData:
    path:       str
    summary:    dict
    stats:      dict
    cuts:       dict
    keep_intervals_evn: list[tuple[int, int]]
    evns:       np.ndarray            # int64, sorted ascending
    times:      np.ndarray            # float seconds, NaN if missing
    overall_pass: np.ndarray          # bool, AND across channels
    channels:   dict[str, ChannelSeries] = field(default_factory=dict)
    run_number: Optional[int] = None
    live_charge: Optional[dict] = None  # the report's "live_charge" block, or
                                        # None when charge integration is off


def load_report(path: str) -> ReportData:
    with open(path) as f:
        j = json.load(f)
    points = j.get("points") or []
    if not points:
        raise SystemExit(f"{path}: report has no 'points' array")

    by_evn: dict[int, dict[str, tuple[bool, float]]] = {}
    times_by_evn: dict[int, float] = {}
    channel_names: set[str] = set()
    for p in points:
        evn = int(p["associated_evn"])
        ch  = str(p["channel"])
        channel_names.add(ch)
        st  = (p.get("status") == "pass")
        v   = p.get("value")
        v_f = float("nan") if v is None else float(v)
        by_evn.setdefault(evn, {})[ch] = (st, v_f)
        if evn not in times_by_evn:
            t = p.get("associated_timestamp")
            if t is not None:
                times_by_evn[evn] = float(t)

    sorted_evns = sorted(by_evn.keys())
    n_cp = len(sorted_evns)
    evns_arr = np.asarray(sorted_evns, dtype=np.int64)
    times_arr = np.asarray(
        [times_by_evn.get(e, float("nan")) for e in sorted_evns],
        dtype=np.float64,
    )

    # Stable ordering: livetime first, then EPICS alphabetically, then anything
    # else.  Keeps the status row's row labels identical between renders.
    def _key(name: str) -> tuple[int, str]:
        if name == "livetime":         return (0, name)
        if name.startswith("epics:"):  return (1, name)
        return (2, name)
    ordered = sorted(channel_names, key=_key)

    channels: dict[str, ChannelSeries] = {}
    for ch in ordered:
        vals = np.full(n_cp, np.nan)
        passes = np.zeros(n_cp, dtype=bool)
        for i, evn in enumerate(sorted_evns):
            st_v = by_evn[evn].get(ch)
            if st_v is not None:
                passes[i] = st_v[0]
                vals[i]   = st_v[1]
        label = ch.split(":", 1)[1] if ch.startswith("epics:") else ch
        channels[ch] = ChannelSeries(
            name=ch, label=label,
            is_livetime=(ch == "livetime"),
            is_epics=ch.startswith("epics:"),
            values=vals, pass_mask=passes,
        )

    overall = np.ones(n_cp, dtype=bool)
    for s in channels.values():
        overall &= s.pass_mask

    keep_raw = j.get("keep_intervals") or []
    keep_intervals_evn: list[tuple[int, int]] = []
    for p in keep_raw:
        if isinstance(p, (list, tuple)) and len(p) == 2:
            try:
                keep_intervals_evn.append((int(p[0]), int(p[1])))
            except (TypeError, ValueError):
                pass

    run_number = j.get("run_number")
    if run_number is not None:
        try: run_number = int(run_number)
        except (TypeError, ValueError): run_number = None

    lc = j.get("live_charge")
    if not isinstance(lc, dict):
        lc = None

    return ReportData(
        path=path,
        summary=j.get("summary") or {},
        stats=j.get("stats") or {},
        cuts=j.get("cuts") or {},
        keep_intervals_evn=keep_intervals_evn,
        evns=evns_arr, times=times_arr,
        overall_pass=overall,
        channels=channels,
        run_number=run_number,
        live_charge=lc,
    )


# ============================================================================
# Derived quantities + shared plot helpers
# ============================================================================

def get_x(report: ReportData, x_kind: str) -> tuple[np.ndarray, str]:
    """Return (x_array, x_label) for the requested axis kind."""
    if x_kind == "evn":
        return report.evns.astype(float), "associated event number"
    return report.times, "associated timestamp [s]"


def compute_datarate_hz(evns: np.ndarray, times: np.ndarray) -> np.ndarray:
    """Per-checkpoint physics rate in Hz, bridged across checkpoints whose
    associated_timestamp is NaN.

    A NaN time arises when a slow row's `event_number_at_arrival` is not
    present in the events/recon tree (e.g. EPICS arrived before the first
    physics event, or anchored on a physics event that was filtered out
    at replay time).  Naïvely diffing would propagate the NaN into *two*
    adjacent intervals — the one ending at the NaN and the one leaving
    it — so the rate line would gap on both sides of every NaN-time
    checkpoint.

    Bridging instead: rate[i] = (evn[i] − evn[j]) / (t[i] − t[j]) where j
    is the most recent finite-time predecessor.  NaN-time checkpoints
    keep their rate as NaN; the plot helper masks them out so the line
    stays unbroken across the bridge."""
    n = len(evns)
    out = np.full(n, np.nan)
    last = -1
    for i in range(n):
        ti = times[i]
        if not np.isfinite(ti):
            continue
        if last >= 0:
            de = float(evns[i] - evns[last])
            dt = ti - times[last]
            if dt > 0 and de >= 0:
                out[i] = de / dt
        last = i
    return out


def _finite_xy(x: np.ndarray, y: np.ndarray
               ) -> tuple[np.ndarray, np.ndarray]:
    """(x, y) pairs restricted to finite-on-both indices.  matplotlib starts
    a new line segment at any NaN, so we drop NaN rows before plotting to
    keep the bridged rate trace continuous on both axis kinds."""
    m = np.isfinite(x) & np.isfinite(y)
    return x[m], y[m]


def reject_segments(report: ReportData, x: np.ndarray
                    ) -> list[tuple[float, float]]:
    """Half-open intervals on the x-axis where physics events would be
    *rejected*.  An interval (x_i, x_{i+1}) is rejected unless both
    endpoints' overall verdicts pass — same rule replay_filter uses to
    build keep_intervals."""
    out: list[tuple[float, float]] = []
    n = len(x)
    if n < 2:
        return out
    op = report.overall_pass
    for i in range(n - 1):
        if op[i] and op[i + 1]:
            continue
        a, b = x[i], x[i + 1]
        if math.isnan(a) or math.isnan(b):
            continue
        out.append((a, b))
    return out


def shade_rejected(ax, segs: Iterable[tuple[float, float]],
                   color: str = "0.55", alpha: float = 0.22) -> None:
    """Dim the rejected segments on `ax`.  Drawn at zorder 0 so plot data
    sits on top.  Single-colour shading, low alpha — the pass regions are
    the visually prominent ones."""
    for a, b in segs:
        ax.axvspan(a, b, color=color, alpha=alpha, lw=0, zorder=0)


# ----- Per-row drawing primitives -------------------------------------------

_YLABEL_KW = dict(rotation=0, ha="right", va="center", labelpad=8)


def _draw_status_row(ax, report: ReportData, x: np.ndarray) -> None:
    """One step trace per channel, vertically offset by channel index.
    pass = top of band (i+1), fail = bottom (i).  Y-tick labels are the
    channel display labels.  Drops NaN-x samples so the step doesn't
    break across NaN-time checkpoints (matplotlib breaks the line at any
    NaN in either x or y)."""
    n_ch = len(report.channels)
    for i, s in enumerate(report.channels.values()):
        y = s.pass_mask.astype(float) + i
        xx, yy = _finite_xy(x, y)
        ax.step(xx, yy, where="post", lw=1.1, color=f"C{i % 10}")
        ax.axhline(i, color="0.85", lw=0.5, zorder=0)
    ax.set_yticks([i + 0.5 for i in range(n_ch)])
    ax.set_yticklabels([s.label for s in report.channels.values()])
    ax.set_ylim(-0.15, n_ch + 0.15)
    ax.set_ylabel("cut status", **_YLABEL_KW)
    ax.grid(axis="x", which="major", alpha=0.25)


def _draw_livetime_rate(ax_lt, report: ReportData, x: np.ndarray):
    """Livetime (left axis) + data rate (right axis), single x.  Returns
    the right axis so the caller can shade it and link x with siblings."""
    lt = report.channels.get("livetime")
    if lt is not None:
        x_lt, y_lt = _finite_xy(x, lt.values)
        ax_lt.plot(x_lt, y_lt, color="C0", lw=1.0)
        ax_lt.set_ylabel("livetime\n[%]", color="C0", **_YLABEL_KW)
        ax_lt.tick_params(axis="y", labelcolor="C0")
    else:
        ax_lt.set_ylabel("livetime\n[%] (n/a)", **_YLABEL_KW)

    rates = compute_datarate_hz(report.evns, report.times)
    x_rt, y_rt = _finite_xy(x, rates)
    ax_rt = ax_lt.twinx()
    ax_rt.plot(x_rt, y_rt, color="C3", lw=1.0)
    # Right-side label: keep horizontal but anchor on the left side of the
    # text so it doesn't overlap the next column.
    ax_rt.set_ylabel("data rate\n[Hz]", color="C3",
                     rotation=0, ha="left", va="center", labelpad=8)
    ax_rt.tick_params(axis="y", labelcolor="C3")
    ax_lt.grid(axis="x", which="major", alpha=0.25)
    return ax_rt


def _title_for(report: ReportData) -> str:
    pieces: list[str] = []
    if report.run_number is not None:
        pieces.append(f"run {report.run_number}")
    summ = report.summary or {}
    n_in   = summ.get("n_physics_in")
    n_pass = summ.get("n_physics_pass")
    rate   = summ.get("physics_pass_rate")
    if n_in is not None and n_pass is not None:
        pct = f"{rate * 100:.1f}%" if isinstance(rate, (int, float)) else "?"
        pieces.append(f"physics {n_pass:,} / {n_in:,} kept ({pct})")
    n_slow = summ.get("n_slow_events")
    if n_slow is not None:
        pieces.append(f"{n_slow:,} checkpoints")
    lc = report.live_charge
    if isinstance(lc, dict):
        # value_nC is the canonical key; tolerate the older "value" name
        # so reports written before the unit assumption was pinned down
        # still render with a charge label.
        v = lc.get("value_nC")
        if not isinstance(v, (int, float)):
            v = lc.get("value")
        if isinstance(v, (int, float)):
            pieces.append(f"Q_live = {v:.3g} nC")
    return " — ".join(pieces) if pieces else os.path.basename(report.path)


# ============================================================================
# CLI rendering
# ============================================================================

def render_cli(report: ReportData, out_path: str, x_kind: str) -> None:
    """One stacked figure: status / livetime+rate / each EPICS channel."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    epics = [s for s in report.channels.values() if s.is_epics]
    n_rows = 2 + len(epics)

    # Heights tuned to keep status + livetime panels readable while EPICS
    # rows stay scannable when many channels are configured.
    heights = [1.5] + [1.4] + [1.1] * len(epics)
    fig_h = max(4.5, sum(heights) * 0.95)
    fig, axes = plt.subplots(
        n_rows, 1, sharex=True,
        figsize=(15, fig_h),
        gridspec_kw={"height_ratios": heights},
        constrained_layout=True,
    )
    if n_rows == 1:
        axes = np.array([axes])

    x, xlabel = get_x(report, x_kind)
    segs = reject_segments(report, x)

    _draw_status_row(axes[0], report, x)
    shade_rejected(axes[0], segs)

    ax_rt = _draw_livetime_rate(axes[1], report, x)
    shade_rejected(axes[1], segs)
    shade_rejected(ax_rt, segs)

    for k, s in enumerate(epics):
        ax = axes[2 + k]
        xx, yy = _finite_xy(x, s.values)
        ax.plot(xx, yy, lw=1.0, color="C2")
        ax.set_ylabel(s.label, fontsize=9, **_YLABEL_KW)
        ax.grid(axis="x", which="major", alpha=0.25)
        shade_rejected(ax, segs)

    axes[-1].set_xlabel(xlabel)
    # Align horizontal ylabels at a common x-position so every panel's
    # plot area starts at the same column (constrained_layout sizes the
    # left margin once for the widest label).
    fig.align_ylabels(axes)
    fig.suptitle(_title_for(report), fontsize=11)
    fig.savefig(out_path, dpi=140)
    plt.close(fig)
    print(f"replay_report_viewer: wrote {out_path}", file=sys.stderr)


# ============================================================================
# GUI
# ============================================================================

def run_gui(initial_path: Optional[str]) -> int:
    import matplotlib
    matplotlib.use("QtAgg")
    from matplotlib.backends.backend_qtagg import (
        FigureCanvasQTAgg, NavigationToolbar2QT,
    )
    from matplotlib.figure import Figure

    from PyQt6.QtCore import Qt
    from PyQt6.QtGui import QAction
    from PyQt6.QtWidgets import (
        QApplication, QComboBox, QFileDialog, QHBoxLayout, QLabel,
        QMainWindow, QMenu, QPushButton, QSizePolicy, QToolButton,
        QVBoxLayout, QWidget,
    )

    app = QApplication.instance() or QApplication(sys.argv)

    class Viewer(QMainWindow):
        def __init__(self) -> None:
            super().__init__()
            self.setWindowTitle("replay_filter report viewer")
            self.resize(1500, 900)

            self.report: Optional[ReportData] = None
            self.x_kind: str = "time"
            self.selected_epics: set[str] = set()

            central = QWidget()
            self.setCentralWidget(central)
            v = QVBoxLayout(central)

            # ── Top toolbar ────────────────────────────────────────────────
            top = QHBoxLayout()
            v.addLayout(top)

            self.btn_open = QPushButton("Open report…")
            self.btn_open.clicked.connect(self.on_open)
            top.addWidget(self.btn_open)

            top.addSpacing(12)
            top.addWidget(QLabel("x-axis:"))
            self.cmb_x = QComboBox()
            self.cmb_x.addItem("associated_timestamp", "time")
            self.cmb_x.addItem("associated_evn",       "evn")
            self.cmb_x.currentIndexChanged.connect(self._on_xaxis)
            top.addWidget(self.cmb_x)

            top.addSpacing(12)
            self.btn_epics = QToolButton()
            self.btn_epics.setText("EPICS channels ▾")
            self.btn_epics.setPopupMode(
                QToolButton.ToolButtonPopupMode.InstantPopup)
            self.menu_epics = QMenu(self.btn_epics)
            self.btn_epics.setMenu(self.menu_epics)
            top.addWidget(self.btn_epics)

            top.addStretch(1)

            self.btn_save = QPushButton("Save screenshot…")
            self.btn_save.clicked.connect(self.on_save)
            top.addWidget(self.btn_save)

            # ── Plot canvas ────────────────────────────────────────────────
            self.fig = Figure(constrained_layout=True)
            self.canvas = FigureCanvasQTAgg(self.fig)
            self.canvas.setSizePolicy(
                QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)
            v.addWidget(self.canvas, 1)

            self.nav = NavigationToolbar2QT(self.canvas, self)
            v.addWidget(self.nav)

            self.lbl_status = QLabel("No report loaded")
            self.lbl_status.setTextInteractionFlags(
                Qt.TextInteractionFlag.TextSelectableByMouse)
            v.addWidget(self.lbl_status)

            if initial_path:
                self._load(initial_path)

        # ── Slots ─────────────────────────────────────────────────────────
        def on_open(self) -> None:
            path, _ = QFileDialog.getOpenFileName(
                self, "Open replay_filter report", "",
                "JSON reports (*.report.json *.json);;All files (*)")
            if path:
                self._load(path)

        def on_save(self) -> None:
            path, _ = QFileDialog.getSaveFileName(
                self, "Save figure", "",
                "PNG (*.png);;PDF (*.pdf);;SVG (*.svg)")
            if path:
                self.fig.savefig(path, dpi=150)
                self.lbl_status.setText(f"saved {path}")

        def _on_xaxis(self) -> None:
            data = self.cmb_x.currentData()
            self.x_kind = data or "time"
            self._replot()

        def _on_epics_toggled(self) -> None:
            self.selected_epics = {
                a.data() for a in self.menu_epics.actions() if a.isChecked()
            }
            self._replot()

        # ── Loading + drawing ─────────────────────────────────────────────
        def _load(self, path: str) -> None:
            try:
                self.report = load_report(path)
            except Exception as e:
                self.lbl_status.setText(f"error loading {path}: {e}")
                return

            self.menu_epics.clear()
            for s in self.report.channels.values():
                if not s.is_epics:
                    continue
                act = QAction(s.label, self)
                act.setCheckable(True)
                act.setChecked(True)
                act.setData(s.name)
                act.toggled.connect(lambda _checked: self._on_epics_toggled())
                self.menu_epics.addAction(act)
            self.selected_epics = {
                s.name for s in self.report.channels.values() if s.is_epics
            }

            self.lbl_status.setText(
                f"{path} — {_title_for(self.report)}")
            self._replot()

        def _replot(self) -> None:
            self.fig.clear()
            r = self.report
            if r is None:
                self.canvas.draw_idle()
                return

            x, xlabel = get_x(r, self.x_kind)
            segs = reject_segments(r, x)

            axes = self.fig.subplots(
                3, 1, sharex=True,
                gridspec_kw={"height_ratios": [1.4, 1.4, 1.7]},
            )

            _draw_status_row(axes[0], r, x)
            shade_rejected(axes[0], segs)

            ax_rt = _draw_livetime_rate(axes[1], r, x)
            shade_rejected(axes[1], segs)
            shade_rejected(ax_rt, segs)

            ax_ep = axes[2]
            sel = [r.channels[n] for n in self.selected_epics
                   if n in r.channels]
            sel.sort(key=lambda s: s.label)
            if sel:
                for s in sel:
                    xx, yy = _finite_xy(x, s.values)
                    ax_ep.plot(xx, yy, lw=1.0, label=s.label)
                ax_ep.set_ylabel("EPICS", **_YLABEL_KW)
                ax_ep.legend(loc="upper right", fontsize=8, ncol=min(len(sel), 4))
                ax_ep.grid(axis="x", which="major", alpha=0.25)
            else:
                ax_ep.text(
                    0.5, 0.5,
                    "No EPICS channels selected — pick from the menu above.",
                    transform=ax_ep.transAxes, ha="center", va="center",
                    color="0.4")
                ax_ep.set_yticks([])
            shade_rejected(ax_ep, segs)

            axes[-1].set_xlabel(xlabel)
            self.fig.align_ylabels(axes)
            self.fig.suptitle(_title_for(r), fontsize=11)
            self.canvas.draw_idle()

    w = Viewer()
    w.show()
    return app.exec()


# ============================================================================
# Entry point
# ============================================================================

def main(argv: Optional[list[str]] = None) -> int:
    ap = argparse.ArgumentParser(
        description="View a replay_filter JSON cut report.")
    ap.add_argument("report", nargs="?", default=None,
                    help="Path to the *.report.json file.")
    ap.add_argument("--cli", action="store_true",
                    help="Render a static figure and exit (no GUI).")
    ap.add_argument("-o", "--out", default=None,
                    help="(--cli) output figure path "
                         "(default: <report stem>.png).")
    ap.add_argument("--evn", action="store_true",
                    help="Use associated_evn for the x-axis "
                         "(default: associated_timestamp).")
    args = ap.parse_args(argv)

    if args.cli:
        if not args.report:
            ap.error("--cli requires a report path")
        out = args.out or str(Path(args.report).with_suffix(".png"))
        report = load_report(args.report)
        render_cli(report, out, "evn" if args.evn else "time")
        return 0

    return run_gui(args.report)


if __name__ == "__main__":
    sys.exit(main())
