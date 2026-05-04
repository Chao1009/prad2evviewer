"""
GEM rendering — data shaping + QPainter drawing, shared by:

* ``gem_event_viewer.py`` (interactive GUI via ``GemEventCanvas``)
* ``gem_cluster_view.py`` (JSON → PNG batch thin wrapper)
* ``gem_layout.py`` (strip layout PNG thin wrapper)

Two draw entry points:

* ``draw_event_panels`` — N detectors side by side, coloured strips + cluster
  markers + 2D hits + legend + colorbars.  Consumes ``process_zs_hits``
  output + ``build_det_list_from_gemsys`` output.
* ``draw_layout`` — one detector, strip positions from ``build_strip_layout``,
  APV boundary dashed lines, beam hole.

Shared colour LUTs (``CMAP_WINTER_RGB`` / ``CMAP_AUTUMN_RGB``) reproduce
matplotlib's ``cm.winter`` / ``cm.autumn`` two-stop gradients — linear so
an analytic lookup suffices.
"""

from __future__ import annotations

import json
from collections import defaultdict
from typing import Dict, Iterable, List, Optional, Tuple

from PyQt6.QtCore import QPointF, QRectF, Qt
from PyQt6.QtGui import (
    QBrush,
    QColor,
    QFont,
    QFontDatabase,
    QFontMetrics,
    QLinearGradient,
    QPainter,
    QPen,
    QPolygonF,
)

from gem_strip_map import map_apv_strips, map_strip


def _font(size: float, bold: bool = False) -> QFont:
    """Default GUI font at the requested size.  Using QFont('Sans', …) yields
    empty-glyph boxes on Windows and many offscreen Qt installs because
    "Sans" is not a real face — ask the system for whatever its GUI font is."""
    f = QFontDatabase.systemFont(QFontDatabase.SystemFont.GeneralFont)
    f.setPointSizeF(size)
    if bold:
        f.setWeight(QFont.Weight.Bold)
    return f


# =============================================================================
# PRad-II GEM mechanical layout (transcribed from
# mpd_gem_view_ssp/gui/experiment_setup/PRadSetup.cpp).
#
# Each PRad GEM chamber has six horizontal G10 spacers and two vertical
# spacers; positions are given for det_pos == 0 (left chamber).  The right
# chamber is the same physical part rotated 180°, so we mirror around the
# chamber centre.
# =============================================================================

_SPACER_VERT_X_LEFT = (183.8, 366.2)                                    # mm
_SPACER_HORIZ_Y_LEFT = (171.0, 347.0, 523.0, 729.8, 920.8, 1096.8)      # mm

# APV box geometry, all in chamber-mm.  ``MARGIN`` is what we reserve in the
# panel transform to keep boxes inside the panel bounds.
_APV_THICKNESS_MM = 40.0
_APV_GAP_MM = 4.0
_APV_MARGIN_MM = 50.0
_SPLIT_BOX_RATIO = 0.5      # split (half-active) APVs draw at this fraction of thickness

# Per-chamber accent colors — mirror PRadSetup's red (left) / dark-green
# (right) convention.  Saturated so they read on dark and light backgrounds.
_DET_POS_FILL = {
    0: (217,  83,  79),  # left chamber red
    1: ( 61, 139,  61),  # right chamber green
}


# =============================================================================
# Colour LUTs
# =============================================================================

# matplotlib cm.winter / cm.autumn are linear two-stop gradients.  Reproduce
# them as 256-entry RGB tables; callers index with fraction*255.
def _linear_lut(start: Tuple[int, int, int],
                end:   Tuple[int, int, int],
                n: int = 256) -> List[Tuple[int, int, int]]:
    return [(
        int(round(start[0] + (end[0] - start[0]) * i / (n - 1))),
        int(round(start[1] + (end[1] - start[1]) * i / (n - 1))),
        int(round(start[2] + (end[2] - start[2]) * i / (n - 1))),
    ) for i in range(n)]


# cm.winter: (0,0,255) -> (0,255,127) — blue to teal
CMAP_WINTER_RGB: List[Tuple[int, int, int]] = _linear_lut((0, 0, 255), (0, 255, 127))
# cm.autumn: (255,0,0) -> (255,255,0) — red to yellow
CMAP_AUTUMN_RGB: List[Tuple[int, int, int]] = _linear_lut((255, 0, 0), (255, 255, 0))


def _lut_color(lut: List[Tuple[int, int, int]], frac: float) -> QColor:
    if frac < 0.0:
        frac = 0.0
    elif frac > 1.0:
        frac = 1.0
    return QColor(*lut[int(round(frac * (len(lut) - 1)))])


# =============================================================================
# Data shaping — gem_map + layout
# =============================================================================


def load_gem_map(path: str):
    """Parse gem_map.json → (layers, apvs, hole, raw)."""
    with open(path, encoding="utf-8") as f:
        raw = json.load(f)
    layers = {l["id"]: l for l in raw["layers"]}
    apvs = [e for e in raw["apvs"] if "crate" in e]
    hole = raw.get("hole", None)
    return layers, apvs, hole, raw


def build_strip_layout(layers, apvs, hole, raw):
    """Build per-detector strip positions (for the static layout view).

    X strips near the beam hole are shortened (split APVs drop 16 channels);
    Y strips crossing the hole split into top/bottom segments.
    """
    detectors = {}
    strips_per_apv = raw.get("apv_channels", 128)

    for det_id, layer in layers.items():
        x_pitch = layer["x_pitch"]
        y_pitch = layer["y_pitch"]
        y_size = layer["y_apvs"] * strips_per_apv * y_pitch
        detectors[det_id] = {
            "name": layer["name"],
            "x_size": 0,       # computed from actual strips below
            "y_size": y_size,
            "x_pitch": x_pitch,
            "y_pitch": y_pitch,
            "x_strips": [],
            "y_strips": [],
            "x_apv_edges": set(),
            "y_apv_edges": set(),
            "x_apv_boxes": [],     # APV chips on the chamber's top/bottom edges
            "y_apv_boxes": [],     # APV chips on the chamber's left/right edges
            "det_pos": 0,          # 0 = left chamber, 1 = right chamber (mirrored)
            "spacer_x": [],        # x positions of vertical spacer dashed lines
            "spacer_y": [],        # y positions of horizontal spacer dashed lines
        }

    apv_ch = raw.get("apv_channels", 128)
    ro_center = raw.get("readout_center", 32)

    all_x_strips = {det_id: set() for det_id in detectors}
    match_strips = {det_id: [] for det_id in detectors}
    apv_data = []

    for apv in apvs:
        det_id = apv["det"]
        if det_id not in detectors:
            continue
        plane = apv["plane"]
        plane_strips = map_apv_strips(apv, apv_channels=apv_ch, readout_center=ro_center)
        apv_data.append((apv, det_id, plane, plane_strips))
        if plane == "X":
            all_x_strips[det_id].update(plane_strips)
            if apv.get("match", ""):
                match_strips[det_id].extend(plane_strips)

    for det_id, strips in all_x_strips.items():
        if strips:
            det = detectors[det_id]
            det["x_size"] = (max(strips) + 1) * det["x_pitch"]

    if hole:
        hw = hole["width"]
        hh = hole["height"]
    else:
        hw = hh = 0

    ref_det_id = min(detectors.keys())
    ref_det = detectors[ref_det_id]
    if hole and match_strips[ref_det_id]:
        ms = match_strips[ref_det_id]
        hx = (min(ms) + max(ms) + 1) / 2 * ref_det["x_pitch"]
        hy = ref_det["y_size"] / 2
        hole_x0, hole_x1 = hx - hw / 2, hx + hw / 2
        hole_y0, hole_y1 = hy - hh / 2, hy + hh / 2
        hole["x_center"] = hx
        hole["y_center"] = hy
    else:
        hole_x0 = hole_x1 = hole_y0 = hole_y1 = -1

    for apv, det_id, plane, plane_strips in apv_data:
        det = detectors[det_id]
        match = apv.get("match", "")
        is_split = bool(match)
        pos = apv.get("pos", 0)
        orient = apv.get("orient", 0)
        det_pos = apv.get("det_pos", 0)
        det["det_pos"] = det_pos  # consistent across APVs of the same chamber
        if plane == "X":
            pitch = det["x_pitch"]
            strip_positions = sorted(set(plane_strips))
            x_min = min(strip_positions) * pitch
            x_max = (max(strip_positions) + 1) * pitch
            if match == "+Y" and hole:
                y0_edge, y1_edge = hole_y1, det["y_size"]
            elif match == "-Y" and hole:
                y0_edge, y1_edge = 0, hole_y0
            else:
                y0_edge, y1_edge = 0, det["y_size"]
            det["x_apv_edges"].add((x_min, y0_edge, y1_edge))
            det["x_apv_edges"].add((x_max, y0_edge, y1_edge))
            for s in plane_strips:
                det["x_strips"].append((s * pitch, y0_edge, y1_edge))
            # APV chip box outside the chamber frame.  Convention: orient
            # 0 = readout cables on top of chamber, orient 1 = bottom.
            det["x_apv_boxes"].append({
                "pos": pos,
                "side": "top" if orient == 0 else "bottom",
                "x_min": x_min,
                "x_max": x_max,
                "label": str(pos),
                "split": is_split,
                "match": match,
            })
        elif plane == "Y":
            pitch = det["y_pitch"]
            strip_positions = sorted(set(plane_strips))
            y_min = min(strip_positions) * pitch
            y_max = (max(strip_positions) + 1) * pitch
            det["y_apv_edges"].add(y_min)
            det["y_apv_edges"].add(y_max)
            for s in plane_strips:
                strip_y = s * pitch
                if hole and hole_y0 < strip_y < hole_y1:
                    det["y_strips"].append((strip_y, 0, hole_x0))
                    det["y_strips"].append((strip_y, hole_x1, det["x_size"]))
                else:
                    det["y_strips"].append((strip_y, 0, det["x_size"]))
            # Y APVs all sit on the chamber's outer side: left edge for the
            # left chamber, right edge for the right chamber.
            det["y_apv_boxes"].append({
                "pos": pos,
                "side": "left" if det_pos == 0 else "right",
                "y_min": y_min,
                "y_max": y_max,
                "label": str(pos),
                "split": is_split,
                "match": match,
            })

    # Spacer positions — mirrored for the right chamber (det_pos == 1) since
    # the same physical chamber is rotated 180° when installed on that side.
    for det in detectors.values():
        x_size = det["x_size"] if det["x_size"] > 0 else 1.0
        y_size = det["y_size"] if det["y_size"] > 0 else 1.0
        if det.get("det_pos", 0) == 1:
            det["spacer_x"] = [x_size - x for x in _SPACER_VERT_X_LEFT]
            det["spacer_y"] = [y_size - y for y in _SPACER_HORIZ_Y_LEFT]
        else:
            det["spacer_x"] = list(_SPACER_VERT_X_LEFT)
            det["spacer_y"] = list(_SPACER_HORIZ_Y_LEFT)

    return detectors


# =============================================================================
# Data shaping — per-event (zero-suppressed hits)
# =============================================================================


def build_apv_map(gem_map_apvs: Iterable[dict]) -> Dict[Tuple[int, int, int], dict]:
    """(crate, mpd, adc) -> APV entry from gem_map.json."""
    return {(a["crate"], a["mpd"], a["adc"]): a
            for a in gem_map_apvs if "crate" in a}


def process_zs_hits(zs_apvs, apv_map, detectors, hole, raw):
    """Convert zero-suppressed APV channels to drawable strip segments.

    Returns dict det_id -> {"x": [...], "y": [...]}; each entry is a
    ``(strip_pos, line_start, line_end, charge, cross_talk)`` tuple.
    """
    apv_ch = raw.get("apv_channels", 128)
    ro_center = raw.get("readout_center", 32)

    if hole and "x_center" in hole:
        hx, hy = hole["x_center"], hole["y_center"]
        hw, hh = hole["width"], hole["height"]
        hole_x0, hole_x1 = hx - hw / 2, hx + hw / 2
        hole_y0, hole_y1 = hy - hh / 2, hy + hh / 2
    else:
        hole_x0 = hole_x1 = hole_y0 = hole_y1 = -1

    result: Dict[int, Dict[str, list]] = defaultdict(lambda: {"x": [], "y": []})

    for apv_entry in zs_apvs:
        key = (apv_entry["crate"], apv_entry["mpd"], apv_entry["adc"])
        props = apv_map.get(key)
        if props is None:
            continue

        det_id = props["det"]
        plane = props["plane"]
        match = props.get("match", "")
        pos = props["pos"]
        orient = props["orient"]
        pin_rotate = props.get("pin_rotate", 0)
        shared_pos = props.get("shared_pos", -1)
        hybrid_board = props.get("hybrid_board", True)

        if det_id not in detectors:
            continue
        det = detectors[det_id]

        for ch_str, ch_data in apv_entry.get("channels", {}).items():
            ch = int(ch_str)
            _, plane_strip = map_strip(
                ch, pos, orient,
                pin_rotate=pin_rotate, shared_pos=shared_pos,
                hybrid_board=hybrid_board,
                apv_channels=apv_ch, readout_center=ro_center)

            charge = ch_data["charge"]
            cross_talk = ch_data.get("cross_talk", False)

            if plane == "X":
                strip_pos = plane_strip * det["x_pitch"]
                if match == "+Y" and hole_y1 > 0:
                    s0, s1 = hole_y1, det["y_size"]
                elif match == "-Y" and hole_y0 > 0:
                    s0, s1 = 0, hole_y0
                else:
                    s0, s1 = 0, det["y_size"]
                result[det_id]["x"].append((strip_pos, s0, s1, charge, cross_talk))

            elif plane == "Y":
                strip_pos = plane_strip * det["y_pitch"]
                if hole_y0 > 0 and hole_y0 < strip_pos < hole_y1:
                    result[det_id]["y"].append((strip_pos, 0, hole_x0, charge, cross_talk))
                    result[det_id]["y"].append((strip_pos, hole_x1, det["x_size"], charge, cross_talk))
                else:
                    result[det_id]["y"].append((strip_pos, 0, det["x_size"], charge, cross_talk))

    return dict(result)


def charge_range(det_hits: Dict[int, Dict[str, list]]) -> Tuple[float, float]:
    """Span charges across every hit → (vmin, vmax)."""
    vmax = 0.0
    for h in det_hits.values():
        for arr in (h["x"], h["y"]):
            for item in arr:
                if item[3] > vmax:
                    vmax = item[3]
    return 0.0, vmax if vmax > 0 else 1.0


# =============================================================================
# Data shaping — pull from live GemSystem
# =============================================================================


def build_zs_apvs_from_gemsys(gsys) -> List[dict]:
    """Build a ``zs_apvs`` list (same shape gem_dump emits) from a post-
    ProcessEvent GemSystem."""
    zs_thres = gsys.zero_sup_threshold
    xt_thres = gsys.cross_talk_threshold

    out: List[dict] = []
    n_apvs = gsys.get_n_apvs()
    n_ts = 6  # SSP_TIME_SAMPLES
    for idx in range(n_apvs):
        if not gsys.has_apv_zs_hits(idx):
            continue
        cfg = gsys.get_apv_config(idx)
        channels: Dict[str, dict] = {}
        for ch in range(128):
            if not gsys.is_channel_hit(idx, ch):
                continue
            ts = [gsys.get_processed_adc(idx, ch, t) for t in range(n_ts)]
            max_charge = max(ts)
            max_tb = ts.index(max_charge)
            ped = cfg.pedestal(ch)
            xtalk = (max_charge < ped.noise * xt_thres) and \
                    (max_charge > ped.noise * zs_thres)
            channels[str(ch)] = {
                "charge": max_charge,
                "max_timebin": max_tb,
                "cross_talk": bool(xtalk),
                "ts_adc": ts,
            }
        if channels:
            out.append({
                "crate": cfg.crate_id,
                "mpd": cfg.mpd_id,
                "adc": cfg.adc_ch,
                "channels": channels,
            })
    return out


def build_det_list_from_gemsys(gsys) -> List[dict]:
    """Build the per-detector list (x_clusters, y_clusters, hits_2d) that
    draw_event_panels expects, reading from a post-Reconstruct GemSystem."""
    out: List[dict] = []
    dets = gsys.get_detectors()
    for d in range(gsys.get_n_detectors()):
        det = dets[d]
        entry = {
            "id":       d,
            "name":     det.name,
            "x_pitch":  det.plane_x.pitch,
            "y_pitch":  det.plane_y.pitch,
            "x_strips": det.plane_x.n_apvs * 128,
            "y_strips": det.plane_y.n_apvs * 128,
        }
        for p, pre in ((0, "x"), (1, "y")):
            cls = gsys.get_plane_clusters(d, p)
            entry[pre + "_clusters"] = [
                {
                    "position":     cl.position,
                    "peak_charge":  cl.peak_charge,
                    "total_charge": cl.total_charge,
                    "max_timebin":  cl.max_timebin,
                    "cross_talk":   cl.cross_talk,
                    "size":         len(cl.hits),
                    "hit_strips":   [h.strip for h in cl.hits],
                } for cl in cls
            ]
        entry["hits_2d"] = [
            {"x": h.x, "y": h.y,
             "x_charge": h.x_charge, "y_charge": h.y_charge,
             "x_peak":   h.x_peak,   "y_peak":   h.y_peak,
             "x_size":   h.x_size,   "y_size":   h.y_size}
            for h in gsys.get_hits(d)
        ]
        out.append(entry)
    return out


# =============================================================================
# QPainter drawing — common helpers
# =============================================================================


def _panel_transform(world_w: float, world_h: float,
                     panel: QRectF,
                     *, margin_x: float = 0.0, margin_y: float = 0.0,
                     fit_factor: float = 0.92) -> Tuple[float, float, float]:
    """Uniform fit of a (margin-padded) world box into the panel rect.

    With ``margin_x = margin_y = 0`` (default), behaves as before: world
    (0, 0) sits at ``(ox, oy + world_h * scale)`` in panel coords (Y
    flipped so up = +y).

    Non-zero margins reserve that many world-units on every side of the
    chamber, which is where APV boxes / external markers go.  ``ox/oy``
    still point at the chamber's top-left in panel coords."""
    ext_w = world_w + 2 * margin_x
    ext_h = world_h + 2 * margin_y
    if ext_w <= 0 or ext_h <= 0:
        return 1.0, panel.x(), panel.y()
    scale = min(panel.width() / ext_w, panel.height() / ext_h) * fit_factor
    dw, dh = ext_w * scale, ext_h * scale
    ox_box = panel.x() + (panel.width() - dw) / 2
    oy_box = panel.y() + (panel.height() - dh) / 2
    ox = ox_box + margin_x * scale
    oy = oy_box + margin_y * scale
    return scale, ox, oy


# =============================================================================
# QPainter drawing — APV boxes & mechanical spacers
# =============================================================================


def _draw_spacers(p: QPainter, det: dict,
                  scale: float, ox: float, oy: float,
                  x_size: float, y_size: float,
                  *, alpha: int = 80, width: float = 1.0):
    """Mechanical G10 spacer dashed lines (6 horizontal + 2 vertical per chamber)."""
    pen = QPen(QColor(150, 150, 150, alpha), width)
    pen.setStyle(Qt.PenStyle.DashLine)
    p.setPen(pen)
    p.setBrush(Qt.BrushStyle.NoBrush)

    for x in det.get("spacer_x", []):
        if 0.0 < x < x_size:
            xa = ox + x * scale
            p.drawLine(QPointF(xa, oy), QPointF(xa, oy + y_size * scale))
    for y in det.get("spacer_y", []):
        if 0.0 < y < y_size:
            ya = oy + (y_size - y) * scale
            p.drawLine(QPointF(ox, ya), QPointF(ox + x_size * scale, ya))


def _draw_apv_boxes(p: QPainter, det: dict,
                    scale: float, ox: float, oy: float,
                    x_size: float, y_size: float,
                    *,
                    show_labels: bool = True,
                    thickness_mm: float = _APV_THICKNESS_MM,
                    gap_mm: float = _APV_GAP_MM,
                    split_ratio: float = _SPLIT_BOX_RATIO,
                    fill_alpha: int = 200,
                    border_color: Optional[QColor] = None,
                    border_width: float = 0.6,
                    label_color: Optional[QColor] = None,
                    font_size: float = 7.5):
    """Draw APV chip boxes around the chamber frame with index labels.

    X APVs use ``orient`` (0 → top, 1 → bottom); Y APVs use ``det_pos``
    (0 → left, 1 → right).  Split APVs (those with a ``match`` field)
    are drawn at half thickness as visual stubs."""
    if border_color is None:
        border_color = QColor("#222")
    if label_color is None:
        label_color = QColor("white")

    rgb = _DET_POS_FILL.get(det.get("det_pos", 0), _DET_POS_FILL[0])
    fill = QColor(*rgb); fill.setAlpha(fill_alpha)

    border_pen = QPen(border_color, border_width)

    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)

    def _draw_box(rect: QRectF, label: str):
        p.setPen(border_pen)
        p.setBrush(fill)
        p.drawRect(rect)
        if show_labels and rect.width() > 6 and rect.height() > 6:
            p.setPen(label_color)
            p.setBrush(Qt.BrushStyle.NoBrush)
            p.setFont(_font(font_size, bold=True))
            p.drawText(rect, Qt.AlignmentFlag.AlignCenter, label)

    # X APVs — boxes above (top) or below (bottom) the chamber.
    for box in det.get("x_apv_boxes", []):
        thick = thickness_mm * (split_ratio if box["split"] else 1.0)
        if box["side"] == "top":
            y_lo, y_hi = y_size + gap_mm, y_size + gap_mm + thick
        else:
            y_lo, y_hi = -gap_mm - thick, -gap_mm
        rx = ox + box["x_min"] * scale
        ry = oy + (y_size - y_hi) * scale
        rw = (box["x_max"] - box["x_min"]) * scale
        rh = (y_hi - y_lo) * scale
        _draw_box(QRectF(rx, ry, rw, rh), box["label"])

    # Y APVs — boxes to the left or right of the chamber.
    for box in det.get("y_apv_boxes", []):
        thick = thickness_mm * (split_ratio if box["split"] else 1.0)
        if box["side"] == "left":
            x_lo, x_hi = -gap_mm - thick, -gap_mm
        else:
            x_lo, x_hi = x_size + gap_mm, x_size + gap_mm + thick
        rx = ox + x_lo * scale
        ry = oy + (y_size - box["y_max"]) * scale
        rw = (x_hi - x_lo) * scale
        rh = (box["y_max"] - box["y_min"]) * scale
        _draw_box(QRectF(rx, ry, rw, rh), box["label"])

    p.setRenderHint(QPainter.RenderHint.Antialiasing, False)


def _w2p(scale: float, ox: float, oy: float, world_h: float,
         x: float, y: float) -> QPointF:
    """World → panel point (Y flipped)."""
    return QPointF(ox + x * scale, oy + (world_h - y) * scale)


# =============================================================================
# QPainter drawing — event view
# =============================================================================


_LEGEND_ENTRIES = [
    ("patch-cmap",    "strip hits"),
    ("tri-up-blue",   "X cluster"),
    ("tri-right-red", "Y cluster"),
    ("plus-fg",       "2D hit"),
    ("dash-gray",     "Cross-talk"),
]


def _pick_charge_lut(bg: QColor) -> List[Tuple[int, int, int]]:
    """Dark backgrounds → winter (cool); light → autumn (warm).

    X- and Y-strip orientation already separates the two planes visually
    (vertical vs horizontal lines), so a single colormap for both carries
    enough information and frees up a colorbar's worth of canvas space.
    """
    # Rec. 601 luma — darker than ~45% luma counts as a dark theme.
    luma = 0.299 * bg.red() + 0.587 * bg.green() + 0.114 * bg.blue()
    return CMAP_WINTER_RGB if luma < 115 else CMAP_AUTUMN_RGB


def draw_event_panels(painter: QPainter, canvas: QRectF,
                      detectors: Dict[int, dict],
                      det_list: List[dict],
                      det_hits: Dict[int, Dict[str, list]],
                      hole: Optional[dict],
                      *, title: Optional[str] = None,
                      det_filter: int = -1,
                      bg: Optional[QColor] = None,
                      fg: Optional[QColor] = None):
    """Paint one event (all visible detectors) into ``canvas``.

    Layout: title top (~28 px), legend bottom (~32 px), one colorbar right
    (~80 px).  Remaining area divided horizontally into N panels, one per
    detector.  Strip colormap is theme-dependent — winter on dark, autumn
    on light.
    """
    if bg is None: bg = QColor("white")
    if fg is None: fg = QColor("#222")
    lut = _pick_charge_lut(bg)

    painter.fillRect(canvas, bg)

    if det_filter >= 0:
        det_list = [d for d in det_list if d["id"] == det_filter]

    title_h = 32 if title else 8
    legend_h = 36
    cb_w = 90
    margin = 12

    # Title
    if title:
        painter.setPen(fg)
        painter.setFont(_font(12, bold=True))
        painter.drawText(
            QRectF(canvas.x(), canvas.y() + 4, canvas.width(), title_h - 8),
            Qt.AlignmentFlag.AlignCenter, title)

    # Colorbar (right margin) — based on global charge range
    vmin, vmax = charge_range(det_hits)
    cb_area = QRectF(canvas.right() - cb_w - margin,
                     canvas.y() + title_h,
                     cb_w,
                     canvas.height() - title_h - legend_h - margin)
    _paint_colorbar(painter, cb_area, lut, vmin, vmax, "strip charge", fg)

    # Panel strip
    panels_area = QRectF(canvas.x() + margin,
                         canvas.y() + title_h,
                         canvas.width() - cb_w - 3 * margin,
                         canvas.height() - title_h - legend_h - margin)

    n = len(det_list)
    if n == 0:
        painter.setPen(fg)
        painter.setFont(_font(12))
        painter.drawText(panels_area, Qt.AlignmentFlag.AlignCenter,
                         "(no detectors in event)")
    else:
        ref_key = min(detectors.keys()) if detectors else None
        panel_w = panels_area.width() / n
        for i, dd in enumerate(det_list):
            panel = QRectF(panels_area.x() + i * panel_w,
                           panels_area.y(), panel_w, panels_area.height())
            did = dd["id"]
            dg = detectors.get(did, detectors.get(ref_key, {}))
            _draw_event_panel(painter, panel, dg, dd,
                              det_hits.get(did, {"x": [], "y": []}),
                              hole, vmin, vmax, fg, lut)

    # Legend
    legend = QRectF(canvas.x(), canvas.bottom() - legend_h,
                    canvas.width(), legend_h)
    _paint_legend(painter, legend, _LEGEND_ENTRIES, fg, lut)


def _draw_event_panel(p: QPainter, panel: QRectF,
                      geom: dict, det_data: dict, hits: Dict[str, list],
                      hole: Optional[dict],
                      vmin: float, vmax: float, fg: QColor,
                      lut: List[Tuple[int, int, int]]):
    x_size = geom.get("x_size", 1.0)
    y_size = geom.get("y_size", 1.0)
    x_pitch = geom.get("x_pitch", 1.0)
    y_pitch = geom.get("y_pitch", 1.0)

    x_plane_size = det_data.get("x_strips", 0) * det_data.get("x_pitch", x_pitch)
    y_plane_size = det_data.get("y_strips", 0) * det_data.get("y_pitch", y_pitch)
    if x_plane_size == 0: x_plane_size = x_size
    if y_plane_size == 0: y_plane_size = y_size

    # Leave room above for the title; bottom padding is just enough for the
    # X-cluster triangles (drawn ~8 px below the detector frame).
    TITLE_H = 18
    BOTTOM_PAD = 10
    inner = QRectF(panel.x() + 6, panel.y() + TITLE_H,
                   panel.width() - 12,
                   panel.height() - TITLE_H - BOTTOM_PAD)
    scale, ox, oy = _panel_transform(x_size, y_size, inner,
                                     margin_x=_APV_MARGIN_MM,
                                     margin_y=_APV_MARGIN_MM)

    # Detector outline
    p.setPen(QPen(QColor("#888"), 0))  # cosmetic
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawRect(QRectF(ox, oy, x_size * scale, y_size * scale))

    # Mechanical spacers — light dashed lines, drawn before strips so any
    # strip hits render on top.
    _draw_spacers(p, geom, scale, ox, oy, x_size, y_size,
                  alpha=70, width=0.8)

    # APV chip boxes around the chamber (drawn before strips so cluster
    # markers and 2D-hit pluses can land on top without occlusion).
    _draw_apv_boxes(p, geom, scale, ox, oy, x_size, y_size,
                    show_labels=True, fill_alpha=180,
                    border_width=0.5, font_size=6.5)

    # Beam hole — round in the actual detector; draw as ellipse.
    if hole and "x_center" in hole:
        hx, hy = hole["x_center"], hole["y_center"]
        hw, hh = hole["width"], hole["height"]
        rx = ox + (hx - hw / 2) * scale
        ry = oy + (y_size - (hy + hh / 2)) * scale
        rw, rh = hw * scale, hh * scale
        p.setPen(QPen(QColor("#cc3333"), 1.2))
        p.setBrush(QColor(255, 204, 0, 24))
        p.drawEllipse(QRectF(rx, ry, rw, rh))

    # Strip segments (solid first, then cross-talk dashed).  Both planes
    # use the active colormap — orientation alone identifies X vs Y.
    p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
    span = max(vmax - vmin, 1e-9)
    _draw_charge_strips(p, hits.get("x", []), "X", lut,
                        vmin, span, scale, ox, oy, y_size, x_size)
    _draw_charge_strips(p, hits.get("y", []), "Y", lut,
                        vmin, span, scale, ox, oy, y_size, x_size)

    # Clusters + 2D hits
    _draw_clusters_and_2d(p, det_data, x_plane_size, y_plane_size,
                          x_pitch, y_pitch,
                          scale, ox, oy, y_size, x_size, fg)
    p.setRenderHint(QPainter.RenderHint.Antialiasing, False)

    # Per-panel title
    n_xh = len(hits.get("x", []))
    n_yh = len(hits.get("y", []))
    n_xcl = len(det_data.get("x_clusters", []))
    n_ycl = len(det_data.get("y_clusters", []))
    n_2d = len(det_data.get("hits_2d", []))
    name = det_data.get("name", "GEM?")
    p.setPen(fg)
    p.setFont(_font(9, bold=True))
    p.drawText(QRectF(panel.x(), panel.y() + 2, panel.width(), TITLE_H - 2),
               Qt.AlignmentFlag.AlignCenter,
               f"{name}  X:{n_xh}/{n_xcl}cl   Y:{n_yh}/{n_ycl}cl   2D:{n_2d}")


def _draw_charge_strips(p: QPainter, hits, plane: str, lut,
                        vmin: float, span: float,
                        scale: float, ox: float, oy: float,
                        y_size: float, x_size: float):
    """Draw strip hit segments; cross-talk drawn dashed at low alpha."""
    normal_pen_data: List[Tuple[QColor, float, float, float, float]] = []
    xtalk_pen_data:  List[Tuple[QColor, float, float, float, float]] = []

    for (pos, s0, s1, charge, xtalk) in hits:
        frac = (charge - vmin) / span
        col = _lut_color(lut, frac)
        if plane == "X":
            x1 = ox + pos * scale
            y1 = oy + (y_size - s0) * scale
            y2 = oy + (y_size - s1) * scale
            rec = (col, x1, y1, x1, y2)
        else:
            y1 = oy + (y_size - pos) * scale
            x1 = ox + s0 * scale
            x2 = ox + s1 * scale
            rec = (col, x1, y1, x2, y1)
        (xtalk_pen_data if xtalk else normal_pen_data).append(rec)

    pen = QPen()
    pen.setWidthF(1.4)
    for col, xa, ya, xb, yb in normal_pen_data:
        col = QColor(col); col.setAlpha(230)
        pen.setColor(col); pen.setStyle(Qt.PenStyle.SolidLine)
        p.setPen(pen)
        p.drawLine(QPointF(xa, ya), QPointF(xb, yb))

    pen.setWidthF(0.7)
    pen.setStyle(Qt.PenStyle.DashLine)
    for col, xa, ya, xb, yb in xtalk_pen_data:
        col = QColor(col); col.setAlpha(110)
        pen.setColor(col)
        p.setPen(pen)
        p.drawLine(QPointF(xa, ya), QPointF(xb, yb))


def _draw_clusters_and_2d(p: QPainter, det_data: dict,
                          x_plane_size: float, y_plane_size: float,
                          x_pitch: float, y_pitch: float,
                          scale: float, ox: float, oy: float,
                          y_size: float, x_size: float,
                          fg: QColor):
    # X cluster centres — blue ▲ along the bottom edge
    p.setPen(QPen(QColor("#1f6feb"), 1.2))
    p.setBrush(QColor("#1f6feb"))
    tri_s = 5.0
    for cl in det_data.get("x_clusters", []):
        cx = ox + (cl["position"] + x_plane_size / 2 - x_pitch / 2) * scale
        cy = oy + y_size * scale + 4
        pts = QPolygonF([QPointF(cx, cy - tri_s),
                         QPointF(cx - tri_s, cy + tri_s),
                         QPointF(cx + tri_s, cy + tri_s)])
        p.drawPolygon(pts)

    # Y cluster centres — red ▶ along the left edge
    p.setPen(QPen(QColor("#c03a2b"), 1.2))
    p.setBrush(QColor("#c03a2b"))
    for cl in det_data.get("y_clusters", []):
        cy = oy + (y_size - (cl["position"] + y_plane_size / 2 - y_pitch / 2)) * scale
        cx = ox - 4
        pts = QPolygonF([QPointF(cx - tri_s, cy - tri_s),
                         QPointF(cx, cy),
                         QPointF(cx - tri_s, cy + tri_s)])
        p.drawPolygon(pts)

    # 2D hits — "+" marker drawn in the theme foreground so it's visible on
    # both light (dark ink) and dark (light ink) backgrounds.
    pen = QPen(fg, 2.4)
    pen.setCapStyle(Qt.PenCapStyle.RoundCap)
    p.setPen(pen)
    plus_s = 7.0
    for h in det_data.get("hits_2d", []):
        hx = ox + (h["x"] + x_plane_size / 2 - x_pitch / 2) * scale
        hy = oy + (y_size - (h["y"] + y_plane_size / 2 - y_pitch / 2)) * scale
        p.drawLine(QPointF(hx - plus_s, hy), QPointF(hx + plus_s, hy))
        p.drawLine(QPointF(hx, hy - plus_s), QPointF(hx, hy + plus_s))


def _paint_colorbar(p: QPainter, area: QRectF,
                    lut: List[Tuple[int, int, int]],
                    vmin: float, vmax: float,
                    label: str, fg: QColor):
    """One vertical colourbar spanning ``area`` (min at bottom, max at top)."""
    bar_w = 14
    label_w = area.width() - bar_w - 8
    bar = QRectF(area.x() + 6, area.y() + 4, bar_w, area.height() - 8)

    grad = QLinearGradient(0, bar.bottom(), 0, bar.top())
    n_samples = 6
    for k in range(n_samples):
        t = k / (n_samples - 1)
        r, g, b = lut[int(round(t * (len(lut) - 1)))]
        grad.setColorAt(t, QColor(r, g, b))
    p.fillRect(bar, QBrush(grad))
    p.setPen(QPen(fg, 0))
    p.setBrush(Qt.BrushStyle.NoBrush)
    p.drawRect(bar)

    p.setPen(fg)
    p.setFont(_font(8))
    text_x = bar.right() + 4
    p.drawText(QRectF(text_x, bar.top() - 2, label_w, 12),
               Qt.AlignmentFlag.AlignLeft, f"{vmax:.0f}")
    p.drawText(QRectF(text_x, bar.center().y() - 6, label_w, 12),
               Qt.AlignmentFlag.AlignLeft, label)
    p.drawText(QRectF(text_x, bar.bottom() - 10, label_w, 12),
               Qt.AlignmentFlag.AlignLeft, f"{vmin:.0f}")


def _paint_legend(p: QPainter, area: QRectF, entries, fg: QColor,
                  lut: List[Tuple[int, int, int]]):
    """Horizontal legend row inside ``area`` — equal-width cells."""
    p.setPen(fg)
    p.setFont(_font(9))
    fm = QFontMetrics(p.font())
    n = len(entries)
    if n == 0:
        return
    cell_w = area.width() / n
    mark_w = 18
    for i, (kind, label) in enumerate(entries):
        cx = area.x() + i * cell_w + 8
        cy = area.center().y()
        _draw_legend_glyph(p, kind, cx, cy, fg, lut)
        p.setPen(fg)
        p.drawText(QPointF(cx + mark_w + 4, cy + fm.ascent() / 2 - 2), label)


def _draw_legend_glyph(p: QPainter, kind: str, cx: float, cy: float,
                       fg: QColor, lut: List[Tuple[int, int, int]]):
    if kind == "patch-cmap":
        c = QColor(*lut[len(lut) // 2])
        p.setPen(QPen(c, 0)); p.setBrush(c)
        p.drawRect(QRectF(cx, cy - 5, 14, 10))
    elif kind == "tri-up-blue":
        p.setPen(QPen(QColor("#1f6feb"), 0)); p.setBrush(QColor("#1f6feb"))
        s = 5
        p.drawPolygon(QPolygonF([QPointF(cx + 7, cy - s),
                                 QPointF(cx + 7 - s, cy + s),
                                 QPointF(cx + 7 + s, cy + s)]))
    elif kind == "tri-right-red":
        p.setPen(QPen(QColor("#c03a2b"), 0)); p.setBrush(QColor("#c03a2b"))
        s = 5
        p.drawPolygon(QPolygonF([QPointF(cx + 2, cy - s),
                                 QPointF(cx + 2 + s * 1.4, cy),
                                 QPointF(cx + 2, cy + s)]))
    elif kind == "plus-fg":
        pen = QPen(fg, 2.4); pen.setCapStyle(Qt.PenCapStyle.RoundCap)
        p.setPen(pen)
        s = 6
        p.drawLine(QPointF(cx + 7 - s, cy), QPointF(cx + 7 + s, cy))
        p.drawLine(QPointF(cx + 7, cy - s), QPointF(cx + 7, cy + s))
    elif kind == "dash-gray":
        pen = QPen(QColor("#888"), 1.2); pen.setStyle(Qt.PenStyle.DashLine)
        p.setPen(pen)
        p.drawLine(QPointF(cx, cy), QPointF(cx + 14, cy))


# =============================================================================
# QPainter drawing — static layout view
# =============================================================================


def draw_layout(painter: QPainter, canvas: QRectF,
                det: dict, hole: Optional[dict],
                *, show_every: int = 8,
                title: Optional[str] = None,
                bg: Optional[QColor] = None,
                fg: Optional[QColor] = None):
    """Paint a single detector's strip layout into ``canvas``."""
    if bg is None: bg = QColor("white")
    if fg is None: fg = QColor("#222")

    painter.fillRect(canvas, bg)

    title_h = 32 if title else 8
    legend_h = 32
    margin = 16

    if title:
        painter.setPen(fg)
        painter.setFont(_font(13, bold=True))
        painter.drawText(
            QRectF(canvas.x(), canvas.y() + 4, canvas.width(), title_h - 8),
            Qt.AlignmentFlag.AlignCenter, title)

    panel = QRectF(canvas.x() + margin, canvas.y() + title_h,
                   canvas.width() - 2 * margin,
                   canvas.height() - title_h - legend_h - margin)

    x_size = det["x_size"]
    y_size = det["y_size"]
    scale, ox, oy = _panel_transform(x_size, y_size, panel,
                                     margin_x=_APV_MARGIN_MM,
                                     margin_y=_APV_MARGIN_MM)

    painter.setRenderHint(QPainter.RenderHint.Antialiasing, True)

    # Detector frame
    painter.setPen(QPen(QColor("#555"), 1.5))
    painter.setBrush(Qt.BrushStyle.NoBrush)
    painter.drawRect(QRectF(ox, oy, x_size * scale, y_size * scale))

    # Mechanical spacers (dashed light gray, drawn before strips so the
    # strip lines render on top and remain readable).
    _draw_spacers(painter, det, scale, ox, oy, x_size, y_size,
                  alpha=110, width=1.0)

    # Beam hole — round in the actual detector; draw as ellipse.
    if hole and "x_center" in hole:
        hx, hy = hole["x_center"], hole["y_center"]
        hw, hh = hole["width"], hole["height"]
        rx = ox + (hx - hw / 2) * scale
        ry = oy + (y_size - (hy + hh / 2)) * scale
        rw, rh = hw * scale, hh * scale
        painter.setPen(QPen(QColor("#cc3333"), 2.0))
        painter.setBrush(QColor(255, 204, 0, 36))
        painter.drawEllipse(QRectF(rx, ry, rw, rh))

    # X strips (vertical blue lines, decimated by show_every within each extent group)
    x_by_extent: Dict[Tuple[float, float], list] = {}
    for (x, y0, y1) in det["x_strips"]:
        x_by_extent.setdefault((y0, y1), []).append(x)
    pen = QPen(QColor(70, 130, 180, 150), 0.6)  # steelblue
    painter.setPen(pen)
    for key in sorted(x_by_extent):
        ys = sorted(x_by_extent[key])
        y0, y1 = key
        for i, x in enumerate(ys):
            if i % show_every != 0:
                continue
            xa = ox + x * scale
            ya = oy + (y_size - y0) * scale
            yb = oy + (y_size - y1) * scale
            painter.drawLine(QPointF(xa, ya), QPointF(xa, yb))

    # Y strips (horizontal red lines)
    y_by_extent: Dict[Tuple[float, float], list] = {}
    for (y, x0, x1) in det["y_strips"]:
        y_by_extent.setdefault((x0, x1), []).append(y)
    pen = QPen(QColor(205, 92, 92, 150), 0.6)  # indianred
    painter.setPen(pen)
    for key in sorted(y_by_extent):
        ys = sorted(y_by_extent[key])
        x0, x1 = key
        for i, y in enumerate(ys):
            if i % show_every != 0:
                continue
            ya = oy + (y_size - y) * scale
            xa = ox + x0 * scale
            xb = ox + x1 * scale
            painter.drawLine(QPointF(xa, ya), QPointF(xb, ya))

    # APV boundaries — dashed blue vertical, dashed red horizontal
    pen = QPen(QColor(70, 130, 180, 128), 1.0)
    pen.setStyle(Qt.PenStyle.DashLine)
    painter.setPen(pen)
    for (x, y0, y1) in sorted(det["x_apv_edges"]):
        xa = ox + x * scale
        ya = oy + (y_size - y0) * scale
        yb = oy + (y_size - y1) * scale
        painter.drawLine(QPointF(xa, ya), QPointF(xa, yb))

    pen.setColor(QColor(205, 92, 92, 128))
    painter.setPen(pen)
    for y in sorted(det["y_apv_edges"]):
        ya = oy + (y_size - y) * scale
        painter.drawLine(QPointF(ox, ya), QPointF(ox + x_size * scale, ya))

    painter.setRenderHint(QPainter.RenderHint.Antialiasing, False)

    # APV chip boxes around the chamber perimeter (drawn last so they sit
    # on top of strip lines that would otherwise show through the gap).
    _draw_apv_boxes(painter, det, scale, ox, oy, x_size, y_size,
                    show_labels=True, font_size=8.0)

    # Legend
    n_x_apvs = len(det.get("x_apv_boxes", []))
    n_y_apvs = len(det.get("y_apv_boxes", []))
    legend_entries = [
        ("patch-steelblue", f"X strips ({len(det['x_strips'])})"),
        ("patch-indianred", f"Y strip segments ({len(det['y_strips'])})"),
        ("patch-apv",       f"APVs (X:{n_x_apvs}/Y:{n_y_apvs}, split = stub)"),
        ("dash-spacer",     "Spacers"),
    ]
    if hole:
        legend_entries.append(("patch-yellow", "Beam hole"))
    legend = QRectF(canvas.x(), canvas.bottom() - legend_h,
                    canvas.width(), legend_h)
    _paint_layout_legend(painter, legend, legend_entries, fg, det)


def _paint_layout_legend(p: QPainter, area: QRectF, entries, fg: QColor,
                         det: Optional[dict] = None):
    p.setPen(fg)
    p.setFont(_font(9))
    fm = QFontMetrics(p.font())
    n = len(entries)
    if n == 0:
        return
    cell_w = area.width() / n
    apv_rgb = _DET_POS_FILL.get((det or {}).get("det_pos", 0), _DET_POS_FILL[0])
    for i, (kind, label) in enumerate(entries):
        cx = area.x() + i * cell_w + 16
        cy = area.center().y()
        if kind == "patch-steelblue":
            c = QColor(70, 130, 180, 200)
            p.setPen(QPen(c, 0)); p.setBrush(c)
            p.drawRect(QRectF(cx, cy - 5, 14, 10))
        elif kind == "patch-indianred":
            c = QColor(205, 92, 92, 200)
            p.setPen(QPen(c, 0)); p.setBrush(c)
            p.drawRect(QRectF(cx, cy - 5, 14, 10))
        elif kind == "patch-yellow":
            p.setPen(QPen(QColor("#cc3333"), 1.2))
            p.setBrush(QColor(255, 204, 0, 48))
            p.drawEllipse(QRectF(cx, cy - 5, 14, 10))
        elif kind == "patch-apv":
            c = QColor(*apv_rgb, 220)
            p.setPen(QPen(QColor("#222"), 0.8)); p.setBrush(c)
            # full-thickness sample + half-thickness stub side-by-side.
            p.drawRect(QRectF(cx, cy - 5, 9, 10))
            p.drawRect(QRectF(cx + 11, cy - 2, 4, 4))
        elif kind == "dash-spacer":
            pen = QPen(QColor(150, 150, 150, 160), 1.0)
            pen.setStyle(Qt.PenStyle.DashLine)
            p.setPen(pen)
            p.drawLine(QPointF(cx, cy), QPointF(cx + 14, cy))
        p.setPen(fg)
        p.drawText(QPointF(cx + 20, cy + fm.ascent() / 2 - 2), label)
