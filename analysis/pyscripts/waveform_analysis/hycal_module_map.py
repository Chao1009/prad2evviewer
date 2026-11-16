#!/usr/bin/env python3
"""hycal_module_map.py -- 2D map of the HyCal modules: each module a rectangle of its size at its x-y position, named inside.

    python3 waveform_analysis/hycal_module_map.py

Geometry: database/hycal_map.json, the standard current map (the default map of _common.py and template_extract.py);
PRAD2_DATABASE_DIR overrides the database directory. geo x, y are module centres, sx, sy the sizes (mm).
Drawn: PbWO4 modules (W*, light blue) and lead-glass modules (G*, light orange). Veto (V*) and LMS entries are not HyCal
modules and are not drawn. Axes in mm, front view as database/hycal_module.txt (+X right, +Y up). The beam hole is the
empty grid cells inside the PbWO4 block, outlined dashed.
Output: waveform_analysis/figures/hycal_module_map.pdf (vector, square 40 in page) and hycal_module_map.png (preview).
Standalone: matplotlib only.
"""
import json
import os
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

HERE = Path(__file__).resolve().parent
DB_DIR = Path(os.environ.get("PRAD2_DATABASE_DIR") or HERE.parent.parent.parent / "database")
MAP_FILE = DB_DIR / "hycal_map.json"
OUT_DIR = HERE / "figures"
OUT_STEM = "hycal_module_map"

PAGE_IN = 40.0                                      # square page, inches (large so names are legible when zoomed)
PAD_MM = 40.0                                       # margin around the outermost module edge
FONT_PT = 8.0                                       # module name; checked to fit inside every cell
COLOR = {"PbWO4": "#cfe6ff", "PbGlass": "#ffe0b3"}  # light blue, light orange
LABEL = {"PbWO4": "PbWO4 (W)", "PbGlass": "lead glass (G)"}
EDGE = "#444444"

plt.rcParams["pdf.fonttype"] = 42                   # TrueType text in the PDF (selectable, editable)


def beam_hole(pbwo4):
    """(x0, x1, y0, y1) of the empty grid cells inside the PbWO4 block, or None if there is no gap."""
    sx, sy = pbwo4[0]["geo"]["sx"], pbwo4[0]["geo"]["sy"]
    pos = {(round(e["geo"]["x"], 2), round(e["geo"]["y"], 2)) for e in pbwo4}
    xs = sorted({p[0] for p in pos})
    ys = sorted({p[1] for p in pos})
    gaps = [(x, y) for x in xs for y in ys if (x, y) not in pos]
    if not gaps:
        return None
    gx = [g[0] for g in gaps]
    gy = [g[1] for g in gaps]
    return min(gx) - sx / 2, max(gx) + sx / 2, min(gy) - sy / 2, max(gy) + sy / 2


def main():
    entries = json.loads(MAP_FILE.read_text())
    mods = [e for e in entries if e["t"] in COLOR]
    skipped = [e["n"] for e in entries if e["t"] not in COLOR]
    counts = {t: sum(e["t"] == t for e in mods) for t in COLOR}
    pbwo4 = [e for e in mods if e["t"] == "PbWO4"]

    reach = max(max(abs(e["geo"]["x"]) + e["geo"]["sx"] / 2, abs(e["geo"]["y"]) + e["geo"]["sy"] / 2) for e in mods)
    lim = reach + PAD_MM
    tick_max = int(lim // 100) * 100
    ticks = list(range(-tick_max, tick_max + 1, 100))

    fig = plt.figure(figsize=(PAGE_IN, PAGE_IN))
    ax = fig.add_axes([0.07, 0.07, 0.72, 0.86])
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_aspect("equal", adjustable="box")
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    ax.tick_params(labelsize=18)
    ax.set_xlabel("x (mm)", fontsize=22)
    ax.set_ylabel("y (mm)", fontsize=22)

    texts = []
    for e in mods:
        g = e["geo"]
        ax.add_patch(Rectangle((g["x"] - g["sx"] / 2, g["y"] - g["sy"] / 2), g["sx"], g["sy"],
                               facecolor=COLOR[e["t"]], edgecolor=EDGE, linewidth=0.4))
        texts.append((e["n"], ax.text(g["x"], g["y"], e["n"], ha="center", va="center", fontsize=FONT_PT), g))

    hole = beam_hole(pbwo4)
    if hole:
        x0, x1, y0, y1 = hole
        ax.add_patch(Rectangle((x0, y0), x1 - x0, y1 - y0, fill=False, linestyle="--",
                               edgecolor="black", linewidth=1.2))
        ax.text((x0 + x1) / 2, (y0 + y1) / 2, "beam hole", ha="center", va="center",
                fontsize=FONT_PT, style="italic")

    handles = [Patch(facecolor=COLOR[t], edgecolor=EDGE, label=f"{LABEL[t]}   n = {counts[t]}") for t in COLOR]
    handles.append(Patch(fill=False, linestyle="--", edgecolor="black", label="beam hole (no module)"))
    fig.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.82, 0.92), fontsize=20, frameon=False)
    fig.suptitle(f"HyCal modules: {counts['PbWO4']} PbWO4 (W) + {counts['PbGlass']} lead glass (G) "
                 f"= {len(mods)}   [{MAP_FILE.name}]", fontsize=28, y=0.975)

    # label fit: text extent (data units) against its own cell; ratio > 1 means the name overflows the cell
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    inv = ax.transData.inverted()
    worst, bad = 0.0, []
    for name, txt, g in texts:
        bb = txt.get_window_extent(renderer)
        (u0, v0), (u1, v1) = inv.transform([(bb.x0, bb.y0), (bb.x1, bb.y1)])
        ratio = max((u1 - u0) / g["sx"], (v1 - v0) / g["sy"])
        worst = max(worst, ratio)
        if ratio > 1:
            bad.append(name)

    OUT_DIR.mkdir(parents=True, exist_ok=True)
    pdf = OUT_DIR / f"{OUT_STEM}.pdf"
    png = OUT_DIR / f"{OUT_STEM}.png"
    fig.savefig(pdf)
    fig.savefig(png, dpi=100)
    plt.close(fig)

    print(f"source: {MAP_FILE}")
    print(f"drawn: PbWO4 {counts['PbWO4']} + lead glass {counts['PbGlass']} = {len(mods)}; "
          f"not drawn (not HyCal modules): {len(skipped)} ({' '.join(skipped)})")
    if hole:
        print("beam hole (mm): x %.2f .. %.2f, y %.2f .. %.2f" % hole)
    print(f"label fit: worst name/cell ratio {worst:.2f} (<1 fits); not fitting: {len(bad)} {' '.join(bad)}")
    print(f"pdf: {pdf}")
    print(f"png: {png}")


if __name__ == "__main__":
    main()
