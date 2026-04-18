"""
GEM cluster rendering — shared between the JSON CLI (gem_cluster_view.py)
and the interactive GUI (gem_event_viewer.py).

The module has two layers:

1. Pure geometry + data shaping — APV-driven conversion of zero-suppressed
   hits (+ optional cluster / 2D-hit lists) into drawable per-detector
   structures.  No matplotlib imports here; callable from any context.

2. Matplotlib rendering — ``plot_detector`` / ``draw_event`` take the
   shaped data and paint into an Axes / Figure.

The rendering deliberately matches the look of the original
``gem_cluster_view.py`` JSON-file front-end: blue→winter X strips,
red→autumn Y strips, cross-talk dashed at low alpha, triangle markers
for cluster centres, ``+`` markers for 2-D hits, yellow beam-hole patch.
"""

from __future__ import annotations

from collections import defaultdict
from typing import Dict, Iterable, List, Tuple

import matplotlib.cm as cm
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
from matplotlib.colors import Normalize

from gem_strip_map import map_strip


# -----------------------------------------------------------------------------
# Data shaping
# -----------------------------------------------------------------------------


def build_apv_map(gem_map_apvs: Iterable[dict]) -> Dict[Tuple[int, int, int], dict]:
    """(crate, mpd, adc) -> APV entry from gem_map.json."""
    return {(a["crate"], a["mpd"], a["adc"]): a
            for a in gem_map_apvs if "crate" in a}


def process_zs_hits(zs_apvs, apv_map, detectors, hole, raw):
    """Convert zero-suppressed APV channels to drawable strip segments.

    Matches the original gem_cluster_view logic: APV-driven geometry via
    ``gem_strip_map.map_strip`` + per-APV ``match`` attribute to decide
    half-strip extents near the beam hole.

    Inputs
    ------
    zs_apvs : list of {crate, mpd, adc, channels: {ch_str: {charge, cross_talk, ...}}}
    apv_map : dict (crate,mpd,adc) -> gem_map APV entry
    detectors : dict det_id -> detector geometry
    hole : dict or None (beam-hole geometry with x_center, y_center, width, height)
    raw : gem_map_json raw dict (for apv_channels, readout_center)

    Returns
    -------
    dict det_id -> {"x": [...], "y": [...]} where each entry is a
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


# -----------------------------------------------------------------------------
# Matplotlib drawing
# -----------------------------------------------------------------------------


def plot_detector(ax, det_geom, det_data, det_hits, hole, norm):
    """Draw one detector panel (strips + cluster markers + 2D hits)."""
    x_size = det_geom["x_size"]
    y_size = det_geom["y_size"]
    x_pitch = det_geom["x_pitch"]
    y_pitch = det_geom["y_pitch"]

    x_plane_size = det_data.get("x_strips", 0) * det_data.get("x_pitch", x_pitch)
    y_plane_size = det_data.get("y_strips", 0) * det_data.get("y_pitch", y_pitch)
    if x_plane_size == 0:
        x_plane_size = x_size
    if y_plane_size == 0:
        y_plane_size = y_size

    ax.add_patch(plt.Rectangle((0, 0), x_size, y_size,
                                fill=False, edgecolor="gray", linewidth=1.5))

    if hole and "x_center" in hole:
        hx, hy = hole["x_center"], hole["y_center"]
        hw, hh = hole["width"], hole["height"]
        ax.add_patch(plt.Rectangle((hx - hw / 2, hy - hh / 2), hw, hh,
                                    fill=True, facecolor="#ffcc0018",
                                    edgecolor="#ffcc00", linewidth=1.5,
                                    linestyle="-", zorder=1))

    x_hits = det_hits.get("x", [])
    y_hits = det_hits.get("y", [])

    if not x_hits and not y_hits:
        ax.set_title(f"{det_data.get('name', 'GEM?')} -- no hits")
        _format_axes(ax, x_size, y_size)
        return

    _draw_strips(ax, x_hits, "X", cm.winter, norm)
    _draw_strips(ax, y_hits, "Y", cm.autumn, norm)

    for cl in det_data.get("x_clusters", []):
        cx = cl["position"] + x_plane_size / 2 - x_pitch / 2
        ax.plot(cx, -y_size * 0.02, "^", color="blue", markersize=6,
                clip_on=False, zorder=6)

    for cl in det_data.get("y_clusters", []):
        cy = cl["position"] + y_plane_size / 2 - y_pitch / 2
        ax.plot(-x_size * 0.02, cy, ">", color="red", markersize=6,
                clip_on=False, zorder=6)

    for h in det_data.get("hits_2d", []):
        hx = h["x"] + x_plane_size / 2 - x_pitch / 2
        hy = h["y"] + y_plane_size / 2 - y_pitch / 2
        ax.plot(hx, hy, "+", color="black", markersize=16,
                markeredgewidth=3, zorder=7)

    n_xh = len(x_hits)
    n_yh = len(y_hits)
    n_xcl = len(det_data.get("x_clusters", []))
    n_ycl = len(det_data.get("y_clusters", []))
    n_2d = len(det_data.get("hits_2d", []))
    ax.set_title(f"{det_data.get('name', 'GEM?')} -- "
                 f"X: {n_xh} hits / {n_xcl} cl   "
                 f"Y: {n_yh} hits / {n_ycl} cl   "
                 f"2D: {n_2d}", fontsize=10)
    _format_axes(ax, x_size, y_size)


def _draw_strips(ax, hits, plane, cmap, norm):
    """Draw strip hit segments as colored lines.

    ``hits`` is a list of (strip_pos, line_start, line_end, charge, cross_talk)
    tuples (produced by ``process_zs_hits``).  Cross-talk strips are drawn
    dashed at reduced alpha.
    """
    normal_lines, normal_colors = [], []
    xtalk_lines, xtalk_colors = [], []

    for (pos, s0, s1, charge, xtalk) in hits:
        if plane == "X":
            line = [(pos, s0), (pos, s1)]
        else:
            line = [(s0, pos), (s1, pos)]
        color = cmap(norm(charge))
        if xtalk:
            xtalk_lines.append(line)
            xtalk_colors.append(color)
        else:
            normal_lines.append(line)
            normal_colors.append(color)

    if normal_lines:
        ax.add_collection(LineCollection(normal_lines, colors=normal_colors,
                                          linewidths=1.2, alpha=0.9, zorder=2))
    if xtalk_lines:
        ax.add_collection(LineCollection(xtalk_lines, colors=xtalk_colors,
                                          linewidths=0.6, linestyles="dashed",
                                          alpha=0.4, zorder=2))


def _format_axes(ax, x_size, y_size):
    ax.set_xlim(-x_size * 0.06, x_size * 1.06)
    ax.set_ylim(-y_size * 0.06, y_size * 1.06)
    ax.set_aspect("equal")
    ax.set_xlabel("X (mm)")
    ax.set_ylabel("Y (mm)")


def add_legend(fig):
    handles = [
        mpatches.Patch(color="teal", alpha=0.8, label="X strip hits"),
        mpatches.Patch(color="orangered", alpha=0.8, label="Y strip hits"),
        plt.Line2D([], [], marker="^", color="blue", linestyle="None",
                   markersize=6, label="X cluster center"),
        plt.Line2D([], [], marker=">", color="red", linestyle="None",
                   markersize=6, label="Y cluster center"),
        plt.Line2D([], [], marker="+", color="black", linestyle="None",
                   markeredgewidth=3, markersize=12, label="2D hit"),
        plt.Line2D([], [], color="gray", linestyle="--", linewidth=0.6,
                   alpha=0.5, label="Cross-talk hit"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=4, fontsize=14,
               framealpha=0.9)


def charge_norm(det_hits: Dict[int, Dict[str, list]]) -> Normalize:
    """Build a Normalize spanning the charge range across every hit."""
    all_q: List[float] = []
    for h in det_hits.values():
        all_q += [x[3] for x in h["x"]] + [x[3] for x in h["y"]]
    return Normalize(vmin=0, vmax=max(all_q) if all_q else 1)


def draw_event(fig, detectors, det_list, det_hits, hole, *, title=None,
               det_filter: int = -1):
    """Render one event into an existing Figure.

    Clears the figure, rebuilds one subplot per detector (filtered by
    ``det_filter`` if non-negative), and draws legend + colorbars.

    Returns the list of created Axes (length == number of visible detectors).
    """
    fig.clear()

    if det_filter >= 0:
        det_list = [d for d in det_list if d["id"] == det_filter]

    n = len(det_list)
    if n == 0:
        ax = fig.add_subplot(1, 1, 1)
        ax.set_title("(no detectors in event)")
        ax.axis("off")
        return [ax]

    axes = fig.subplots(1, n)
    if n == 1:
        axes = [axes]
    else:
        axes = list(axes.flat) if hasattr(axes, "flat") else list(axes)

    norm = charge_norm(det_hits)

    ref = detectors[min(detectors.keys())]
    for i, dd in enumerate(det_list):
        did = dd["id"]
        dg = detectors.get(did, ref)
        plot_detector(axes[i], dg, dd,
                      det_hits.get(did, {"x": [], "y": []}), hole, norm)

    for cmap_obj, label in [(cm.winter, "X charge (ADC)"),
                            (cm.autumn, "Y charge (ADC)")]:
        sm = cm.ScalarMappable(cmap=cmap_obj, norm=norm); sm.set_array([])
        cb = fig.colorbar(sm, ax=axes, shrink=0.4, pad=0.01,
                          aspect=30, location="right")
        cb.set_label(label, fontsize=11); cb.ax.tick_params(labelsize=10)

    if title:
        fig.suptitle(title, fontsize=14)
    add_legend(fig)

    return axes


# -----------------------------------------------------------------------------
# GemSystem-backed event builder — for the live GUI (no JSON roundtrip)
# -----------------------------------------------------------------------------


def build_zs_apvs_from_gemsys(gsys) -> List[dict]:
    """Build a ``zs_apvs`` list (same shape gem_dump emits) from a post-
    ProcessEvent GemSystem.

    Walks every APV, reads the zero-suppression mask + processed ADC
    samples, and emits one entry per APV with any surviving channels.
    The per-channel ``charge`` / ``cross_talk`` flags mirror gem_dump's
    thresholds (max over time samples; xtalk if peak below xt_thres × noise
    but above zs_thres × noise).
    """
    zs_thres = gsys.zero_sup_threshold
    xt_thres = gsys.cross_talk_threshold

    out: List[dict] = []
    n_apvs = gsys.get_n_apvs()
    # SSP_TIME_SAMPLES is fixed at 6 in the C++ layer; expose via the first
    # APV's processed-adc accessor to avoid hard-coding the constant here.
    n_ts = 6
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
    ``plot_detector`` expects, reading from a post-Reconstruct GemSystem.
    """
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
