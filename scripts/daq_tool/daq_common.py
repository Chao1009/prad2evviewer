"""
Shared pieces of the daq_tool editors (trigger_mask_editor, fadc_gain_config):
the HyCal module list as the editors display it, FAV3 crate-config block
parsing / rendering, and the common button style.

Importing this module puts the parent scripts/ directory on ``sys.path``,
so a tool must import it before ``hycal_geoview`` (which holds the crate
constants and the Module record).
"""

from __future__ import annotations

import sys
from pathlib import Path
from typing import Dict, Iterator, List, Optional, Tuple

_SCRIPTS_DIR = Path(__file__).resolve().parent.parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from hycal_geoview import (  # noqa: E402
    CRATE_NAMES, NUM_CRATES, THEME, Module, load_modules, place_aux_row,
    themed,
)

CRATE_INDEX: Dict[str, int] = {n: i for i, n in enumerate(CRATE_NAMES)}

# Centre-x of the LMS / V cells in the row below HyCal.
AUX_ROW_X: Dict[str, float] = {
    "LMS1": -200.0, "LMS2": -145.0, "LMS3": -90.0,
    "V1":     35.0, "V2":     90.0, "V3":   145.0, "V4":  200.0,
}


def load_module_info(db_dir: Path) -> List[Module]:
    """All modules of ``db_dir``/hycal_map.json with their DAQ address and
    the LMS / V cells moved into the row below HyCal."""
    return place_aux_row(load_modules(Path(db_dir) / "hycal_map.json"),
                         AUX_ROW_X)


def iter_fav3(text: str, directive: str) -> Iterator[Tuple[int, int, List[str]]]:
    """Walk FAV3 crate-config text and yield ``(crate_index, slot, tokens)``
    for every ``directive`` line (e.g. ``"FAV3_ALLCH_GAIN"``) that follows a
    known ``FAV3_CRATE`` and a valid ``FAV3_SLOT``.  ``tokens`` are the raw
    value tokens after the directive; ``#`` starts a comment.
    """
    crate = -1
    slot = -1
    for raw in text.splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        parts = line.split()
        kw = parts[0]
        if kw == "FAV3_CRATE":
            crate = CRATE_INDEX.get(parts[1] if len(parts) > 1 else "", -1)
        elif kw == "FAV3_SLOT" and len(parts) > 1:
            try:
                slot = int(parts[1])
            except ValueError:
                slot = -1
        elif kw == directive and crate >= 0 and slot >= 0:
            yield crate, slot, parts[1:]


def render_fav3(blocks: Dict[Tuple[int, int], Tuple[str, List[str]]],
                directive: str) -> List[str]:
    """Emit FAV3 crate-config lines for ``{(crate_index, slot): (comment_line,
    values)}``: per crate in index order, each slot (sorted) as its comment
    line, ``FAV3_SLOT n`` and ``<directive> v0 v1 ...``, closed by
    ``FAV3_CRATE end`` and a blank line.  Crates without blocks are omitted.
    """
    lines: List[str] = []
    for ci in range(NUM_CRATES):
        slots = sorted(s for (c, s) in blocks if c == ci)
        if not slots:
            continue
        lines.append(f"FAV3_CRATE {CRATE_NAMES[ci]}")
        for slot in slots:
            comment, values = blocks[(ci, slot)]
            lines.append(comment)
            lines.append(f"FAV3_SLOT {slot}")
            lines.append(f"{directive} {' '.join(values)}")
        lines.append("FAV3_CRATE end")
        lines.append("")
    return lines


def btn_style(checked_color: Optional[str] = None) -> str:
    base = themed(
        f"QPushButton{{background:{THEME.BUTTON};color:{THEME.TEXT};"
        f"border:1px solid {THEME.BORDER};padding:6px 14px;"
        f"font:10pt;border-radius:8px;}}"
        f"QPushButton:hover{{background:{THEME.BUTTON_HOVER};}}")
    if checked_color:
        base += themed(
            f"QPushButton:checked{{background:{checked_color};"
            f"color:{THEME.TEXT};border:1px solid {checked_color};}}")
    return base
