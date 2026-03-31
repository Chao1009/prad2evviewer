#!/usr/bin/env python3
"""
HyCal Cabling GUI
==================
Standalone Tkinter tool for mapping HyCal signal cables to FADC250 crate
channels.  Cables are loaded from a CSV (cable_group_numbers.csv) and
assigned to (crate, slot, channel) positions interactively.

Usage
-----
    python hycal_cabling.py                                    # defaults
    python hycal_cabling.py --cables ../../cable_group_numbers.csv
    python hycal_cabling.py --connections cable_connections.json  # resume

Interaction
-----------
    1. Select a cable in the right-hand table (click or type in filter).
    2. Click an unconnected (red) dot on the crate view.
    3. The cable is connected; the dot turns green.
    4. "Auto-Finish" continues the bundle on the same board.
    Right-click a green dot to disconnect.
"""

from __future__ import annotations

import argparse
import csv
import json
import os
import tkinter as tk
from tkinter import ttk, messagebox, filedialog
from dataclasses import dataclass, field
from typing import Dict, List, Optional, Tuple

# ============================================================================
#  Theme colours (standalone -- no external dependency)
# ============================================================================

class C:
    BG       = "#0d1117"
    PANEL    = "#161b22"
    BORDER   = "#30363d"
    TEXT     = "#c9d1d9"
    DIM      = "#8b949e"
    ACCENT   = "#58a6ff"
    GREEN    = "#3fb950"
    YELLOW   = "#d29922"
    RED      = "#f85149"
    ORANGE   = "#db6d28"
    EMPTY    = "#111418"


# ============================================================================
#  Constants
# ============================================================================

CRATE_NAMES = [f"adchycal{i}" for i in range(1, 8)]
GROUP_TO_CRATE = {i: f"adchycal{i + 1}" for i in range(7)}
CHANNELS_PER_BOARD = 16

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
DEFAULT_BOARD_CONFIG = os.path.join(_SCRIPT_DIR, "board_config.json")
DEFAULT_CABLE_CSV = os.path.join(_SCRIPT_DIR, "cable_group_numbers.csv")


# ============================================================================
#  Data model
# ============================================================================

@dataclass
class Cable:
    board: int          # Cable Board Label
    sub: int            # Sub Label
    number: int         # Cable Number
    group: int          # Group Label (0-6 -> crate hint)
    vpcb_board: int = 0       # VPCB.board (not shown in GUI)
    vpcb_connector: int = 0   # VPCB.connector (not shown in GUI)

    @property
    def label(self) -> str:
        return f"{self.board}.{self.sub}.{self.number}"

    @property
    def bundle(self) -> str:
        return f"{self.board}.{self.sub}"

    def sort_key(self) -> Tuple[int, int, int]:
        return (self.board, self.sub, self.number)


def parse_cable_label(s: str) -> Tuple[int, int, int]:
    parts = s.split(".")
    return int(parts[0]), int(parts[1]), int(parts[2])


class CablingModel:
    """Holds cables, board configuration, and the connection map."""

    def __init__(self):
        self.cables: Dict[str, Cable] = {}
        self.board_config: Dict[str, List[int]] = {}
        # cable_label -> (crate, slot, channel)
        self.connections: Dict[str, Tuple[str, int, int]] = {}
        # (crate, slot, channel) -> cable_label
        self._reverse: Dict[Tuple[str, int, int], str] = {}

    # -- loading ------------------------------------------------------------

    def load_cables_csv(self, path: str) -> int:
        """Load cables from the group-numbers CSV.  Returns count loaded."""
        count = 0
        with open(path, newline="") as f:
            reader = csv.reader(f)
            for row in reader:
                if not row or row[0].startswith("#"):
                    continue
                if len(row) < 6:
                    continue
                # skip header rows
                if row[2] == "Cable Board Label" or row[0] == "VPCB.board":
                    continue
                board_s, sub_s, group_s, num_s = row[2], row[3], row[4], row[5]
                if not board_s or not sub_s or not group_s or not num_s:
                    continue
                try:
                    cable = Cable(int(board_s), int(sub_s),
                                  int(num_s), int(group_s),
                                  int(row[0]), int(row[1]))
                except ValueError:
                    continue
                self.cables[cable.label] = cable
                count += 1
        return count

    def load_board_config(self, path: str):
        with open(path) as f:
            self.board_config = json.load(f)

    # -- queries ------------------------------------------------------------

    def get_cable_at(self, crate: str, slot: int, ch: int) -> Optional[str]:
        return self._reverse.get((crate, slot, ch))

    def has_board(self, crate: str, slot: int) -> bool:
        return slot in self.board_config.get(crate, [])

    def is_connected(self, label: str) -> bool:
        return label in self.connections

    def unconnected_cables(self, filter_str: str = "",
                           group: Optional[int] = None) -> List[Cable]:
        result = []
        filt = filter_str.strip().lower()
        for label, cable in self.cables.items():
            if label in self.connections:
                continue
            if group is not None and cable.group != group:
                continue
            if filt and filt not in label.lower():
                continue
            result.append(cable)
        result.sort(key=Cable.sort_key)
        return result

    def crate_stats(self, crate: str) -> Tuple[int, int]:
        """Return (connected, total) channels for a crate."""
        slots = self.board_config.get(crate, [])
        total = len(slots) * CHANNELS_PER_BOARD
        connected = sum(1 for (cr, sl, ch) in self._reverse
                        if cr == crate)
        return connected, total

    # -- mutations ----------------------------------------------------------

    def connect(self, label: str, crate: str, slot: int, ch: int) -> bool:
        if label not in self.cables:
            return False
        if label in self.connections:
            return False
        if (crate, slot, ch) in self._reverse:
            return False
        if not self.has_board(crate, slot):
            return False
        self.connections[label] = (crate, slot, ch)
        self._reverse[(crate, slot, ch)] = label
        return True

    def disconnect(self, label: str) -> bool:
        if label not in self.connections:
            return False
        key = self.connections.pop(label)
        self._reverse.pop(key, None)
        return True

    def auto_finish(self, label: str) -> List[Tuple[str, str, int, int]]:
        """Auto-connect remaining cables in the same bundle on the same board.

        Starting from the cable *label* that was just connected, fills
        subsequent empty channels on the same board with the next
        unconnected cables in the same bundle (ascending cable number).

        Returns list of (cable_label, crate, slot, channel) created.
        """
        if label not in self.connections:
            return []
        cable = self.cables[label]
        crate, slot, start_ch = self.connections[label]

        # Gather unconnected cables in same bundle with number > this one
        bundle = cable.bundle
        candidates = sorted(
            [c for c in self.cables.values()
             if c.bundle == bundle
             and c.number > cable.number
             and c.label not in self.connections],
            key=lambda c: c.number)

        # Gather empty channels on same board, ch > start_ch
        empty_chs = sorted(
            ch for ch in range(start_ch + 1, CHANNELS_PER_BOARD)
            if (crate, slot, ch) not in self._reverse)

        result = []
        for c, ch in zip(candidates, empty_chs):
            self.connect(c.label, crate, slot, ch)
            result.append((c.label, crate, slot, ch))
        return result

    # -- persistence --------------------------------------------------------

    def save_json(self, path: str):
        """Save connection map as JSON, one entry per line.

        Each entry includes the FADC connection and the VPCB origin:
          "board.sub.number": ["crate", slot, channel, vpcb_board, vpcb_connector]
        """
        entries = sorted(self.connections.items(),
                         key=lambda kv: self.cables[kv[0]].sort_key())
        with open(path, "w") as f:
            f.write("{\n")
            for i, (label, (crate, slot, ch)) in enumerate(entries):
                c = self.cables[label]
                comma = "," if i < len(entries) - 1 else ""
                f.write(f'  "{label}": ["{crate}", {slot}, {ch}, '
                        f'{c.vpcb_board}, {c.vpcb_connector}]{comma}\n')
            f.write("}\n")

    def load_json(self, path: str) -> Tuple[int, int]:
        """Load connections from JSON.  Returns (loaded, skipped).

        Accepts both 3-element [crate, slot, ch] and 5-element
        [crate, slot, ch, vpcb_board, vpcb_connector] entries.
        """
        with open(path) as f:
            data = json.load(f)
        loaded = skipped = 0
        for label, val in data.items():
            crate, slot, ch = val[0], int(val[1]), int(val[2])
            if self.connect(label, crate, slot, ch):
                loaded += 1
            else:
                skipped += 1
        return loaded, skipped


# ============================================================================
#  GUI
# ============================================================================

class CablingGUI:

    COL_W  = 31
    ROW_H  = 26
    DOT_R  = 9
    LEFT_M = 28
    TOP_M  = 12
    BOT_M  = 28
    CANVAS_W = 28 + 20 * 31 + 12   # 660
    CANVAS_H = 12 + 16 * 26 + 28   # 456

    def __init__(self, root: tk.Tk, model: CablingModel):
        self.root = root
        self.model = model
        self._current_crate = CRATE_NAMES[0]
        self._selected_cable: Optional[str] = None
        self._last_connection: Optional[str] = None
        self._dots: Dict[Tuple[int, int], int] = {}
        self._build_ui()

    # ------------------------------------------------------------------
    #  UI construction
    # ------------------------------------------------------------------

    def _build_ui(self):
        self.root.title("HyCal Cable Manager")
        self.root.configure(bg=C.BG)

        style = ttk.Style()
        style.theme_use("clam")
        style.configure(".", background=C.BG, foreground=C.TEXT,
                         fieldbackground=C.PANEL, bordercolor=C.BORDER)
        style.configure("TLabel", background=C.BG, foreground=C.TEXT)
        style.configure("TLabelframe", background=C.BG, foreground=C.ACCENT)
        style.configure("TLabelframe.Label", background=C.BG,
                         foreground=C.ACCENT, font=("Consolas", 9, "bold"))
        style.configure("TButton", background=C.PANEL, foreground=C.TEXT,
                         padding=4)
        style.map("TButton",
                  background=[("active", C.BORDER)],
                  foreground=[("disabled", "#484f58")])
        style.configure("Accent.TButton", background="#1f6feb",
                         foreground="white")
        style.map("Accent.TButton", background=[("active", "#388bfd")])
        style.configure("Danger.TButton", background="#da3633",
                         foreground="white")
        style.map("Danger.TButton", background=[("active", "#f85149")])
        style.configure("Green.TButton", background="#238636",
                         foreground="white")
        style.map("Green.TButton", background=[("active", "#3fb950")])
        style.configure("TCombobox", fieldbackground=C.PANEL,
                         background=C.BORDER, foreground=C.TEXT,
                         selectbackground=C.BORDER, selectforeground=C.TEXT)
        style.map("TCombobox",
                  fieldbackground=[("readonly", C.PANEL),
                                   ("disabled", C.BG)],
                  foreground=[("disabled", "#484f58")])
        # Treeview
        style.configure("Treeview", background=C.PANEL, foreground=C.TEXT,
                         fieldbackground=C.PANEL, rowheight=22,
                         font=("Consolas", 9))
        style.configure("Treeview.Heading", background=C.BORDER,
                         foreground=C.TEXT, font=("Consolas", 9, "bold"))
        style.map("Treeview",
                  background=[("selected", C.ACCENT)],
                  foreground=[("selected", "white")])

        # Combobox listbox
        self.root.option_add("*TCombobox*Listbox.background", C.PANEL)
        self.root.option_add("*TCombobox*Listbox.foreground", C.TEXT)
        self.root.option_add("*TCombobox*Listbox.selectBackground", C.ACCENT)
        self.root.option_add("*TCombobox*Listbox.selectForeground", "white")

        # -- top bar ---------------------------------------------------------
        top = tk.Frame(self.root, bg="#0d1520", height=32)
        top.pack(fill="x")
        tk.Label(top, text="  HYCAL CABLE MANAGER  ", bg="#0d1520",
                 fg=C.GREEN,
                 font=("Consolas", 13, "bold")).pack(side="left", padx=8)
        self._lbl_total = tk.Label(top, text="", bg="#0d1520", fg=C.DIM,
                                    font=("Consolas", 9))
        self._lbl_total.pack(side="right", padx=12)

        # -- main area -------------------------------------------------------
        main = tk.Frame(self.root, bg=C.BG)
        main.pack(fill="both", expand=True, padx=6, pady=4)
        main.columnconfigure(1, weight=1)
        main.rowconfigure(0, weight=1)

        self._build_crate_view(main)
        self._build_cable_panel(main)

        # -- log -------------------------------------------------------------
        log_frame = ttk.LabelFrame(self.root, text=" Log ")
        log_frame.pack(fill="x", padx=6, pady=(0, 4))
        self._log_text = tk.Text(log_frame, height=4, wrap="word",
                                  bg=C.PANEL, fg=C.TEXT, bd=0,
                                  font=("Consolas", 8),
                                  insertbackground=C.TEXT, state="disabled")
        self._log_text.pack(fill="x", padx=4, pady=4)
        self._log_text.tag_configure("info", foreground=C.TEXT)
        self._log_text.tag_configure("ok", foreground=C.GREEN)
        self._log_text.tag_configure("warn", foreground=C.YELLOW)
        self._log_text.tag_configure("error", foreground=C.RED)

        # -- initial draw ----------------------------------------------------
        self._draw_crate()
        self._refresh_cable_table()
        self._update_stats()

    # -- crate view ----------------------------------------------------------

    def _build_crate_view(self, parent):
        frame = ttk.LabelFrame(parent, text=" Crate View ")
        frame.grid(row=0, column=0, sticky="ns", padx=(0, 4))

        # Crate selector
        sel = tk.Frame(frame, bg=C.BG)
        sel.pack(fill="x", padx=6, pady=4)
        tk.Label(sel, text="Crate:", bg=C.BG, fg=C.TEXT,
                 font=("Consolas", 9)).pack(side="left")
        self._crate_var = tk.StringVar(value=self._current_crate)
        self._crate_combo = ttk.Combobox(
            sel, textvariable=self._crate_var, values=CRATE_NAMES,
            width=14, state="readonly", font=("Consolas", 9))
        self._crate_combo.pack(side="right")
        self._crate_combo.bind("<<ComboboxSelected>>", self._on_crate_changed)

        # Canvas
        self._canvas = tk.Canvas(frame, width=self.CANVAS_W,
                                  height=self.CANVAS_H,
                                  bg=C.BG, highlightthickness=0)
        self._canvas.pack(padx=6, pady=(0, 2))

        # Tooltip (floats over canvas)
        self._tooltip = tk.Label(self._canvas, bg=C.BORDER, fg=C.TEXT,
                                  font=("Consolas", 8), padx=4, pady=2,
                                  justify="left")

        # Stats
        self._lbl_crate_stats = tk.Label(frame, text="", bg=C.BG, fg=C.DIM,
                                          font=("Consolas", 9))
        self._lbl_crate_stats.pack(padx=6, pady=(0, 4))

        # Instruction
        self._lbl_hint = tk.Label(frame, text="Select a cable, then click a"
                                  " red dot to connect",
                                  bg=C.BG, fg=C.DIM, font=("Consolas", 8))
        self._lbl_hint.pack(padx=6, pady=(0, 4))

    def _draw_crate(self):
        c = self._canvas
        c.delete("all")
        self._dots.clear()
        self._tooltip.place_forget()

        crate = self._current_crate
        slots = self.model.board_config.get(crate, [])

        # Channel labels (left)
        for ch in range(CHANNELS_PER_BOARD):
            row = 15 - ch
            cy = self.TOP_M + row * self.ROW_H + self.ROW_H // 2
            c.create_text(14, cy, text=str(ch), fill=C.DIM,
                          font=("Consolas", 7))

        for slot_idx in range(20):
            slot = slot_idx + 1
            x0 = self.LEFT_M + slot_idx * self.COL_W + 1
            x1 = x0 + self.COL_W - 2
            y0 = self.TOP_M
            y1 = self.TOP_M + 16 * self.ROW_H
            xc = (x0 + x1) // 2

            if slot in slots:
                # Board background
                c.create_rectangle(x0, y0, x1, y1,
                                   fill=C.PANEL, outline=C.BORDER, width=1)
                # Channel dots
                for ch in range(CHANNELS_PER_BOARD):
                    row = 15 - ch
                    cy = self.TOP_M + row * self.ROW_H + self.ROW_H // 2
                    cable = self.model.get_cable_at(crate, slot, ch)
                    fill = C.GREEN if cable else C.RED

                    tag = f"s{slot}c{ch}"
                    dot_id = c.create_oval(
                        xc - self.DOT_R, cy - self.DOT_R,
                        xc + self.DOT_R, cy + self.DOT_R,
                        fill=fill, outline="", tags=(tag, "dot"))
                    self._dots[(slot, ch)] = dot_id

                    c.tag_bind(tag, "<Button-1>",
                               lambda e, _s=slot, _c=ch:
                                   self._on_dot_click(_s, _c))
                    c.tag_bind(tag, "<Button-3>",
                               lambda e, _s=slot, _c=ch:
                                   self._on_dot_right_click(e, _s, _c))
                    c.tag_bind(tag, "<Enter>",
                               lambda e, _s=slot, _c=ch:
                                   self._on_dot_enter(e, _s, _c))
                    c.tag_bind(tag, "<Leave>",
                               lambda e: self._on_dot_leave())
            else:
                # Empty slot
                c.create_rectangle(x0, y0, x1, y1,
                                   fill=C.EMPTY, outline="")

            # Slot label
            c.create_text(xc, y1 + 12, text=str(slot), fill=C.DIM,
                          font=("Consolas", 7))

        self._update_crate_stats()

    def _refresh_dot(self, slot: int, ch: int):
        dot_id = self._dots.get((slot, ch))
        if dot_id is None:
            return
        cable = self.model.get_cable_at(self._current_crate, slot, ch)
        self._canvas.itemconfigure(dot_id, fill=C.GREEN if cable else C.RED)

    # -- cable panel ---------------------------------------------------------

    def _build_cable_panel(self, parent):
        frame = ttk.LabelFrame(parent, text=" Cables ")
        frame.grid(row=0, column=1, sticky="nsew", padx=(0, 0))

        # Filter row
        filt = tk.Frame(frame, bg=C.BG)
        filt.pack(fill="x", padx=6, pady=4)
        tk.Label(filt, text="Filter:", bg=C.BG, fg=C.TEXT,
                 font=("Consolas", 9)).pack(side="left")
        self._filter_var = tk.StringVar()
        self._filter_var.trace_add("write", lambda *_: self._refresh_cable_table())
        tk.Entry(filt, textvariable=self._filter_var, bg=C.PANEL, fg=C.TEXT,
                 insertbackground=C.TEXT, bd=1, relief="solid",
                 font=("Consolas", 9)).pack(side="left", fill="x",
                                             expand=True, padx=(4, 0))

        # Show connected toggle
        self._show_connected_var = tk.BooleanVar(value=False)
        tk.Checkbutton(filt, text="Show connected", variable=self._show_connected_var,
                       bg=C.BG, fg=C.DIM, selectcolor=C.PANEL,
                       activebackground=C.BG, activeforeground=C.TEXT,
                       font=("Consolas", 8),
                       command=self._refresh_cable_table).pack(side="right", padx=4)

        # Filter hint
        tk.Label(frame, text="e.g. board:1,sub:2  or  1.3.",
                 bg=C.BG, fg=C.DIM,
                 font=("Consolas", 7)).pack(padx=6, anchor="w")

        # Treeview
        tree_frame = tk.Frame(frame, bg=C.BG)
        tree_frame.pack(fill="both", expand=True, padx=6)

        cols = ("board", "sub", "group", "number", "status")
        self._tree = ttk.Treeview(tree_frame, columns=cols, show="headings",
                                   selectmode="browse")
        self._tree.heading("board", text="Board", anchor="w",
                           command=lambda: self._sort_tree("board"))
        self._tree.heading("sub", text="Sub", anchor="w",
                           command=lambda: self._sort_tree("sub"))
        self._tree.heading("group", text="Group", anchor="w",
                           command=lambda: self._sort_tree("group"))
        self._tree.heading("number", text="Number", anchor="w",
                           command=lambda: self._sort_tree("number"))
        self._tree.heading("status", text="Status", anchor="w",
                           command=lambda: self._sort_tree("status"))
        self._tree.column("board", width=55, minwidth=45)
        self._tree.column("sub", width=45, minwidth=35)
        self._tree.column("group", width=50, minwidth=40)
        self._tree.column("number", width=60, minwidth=45)
        self._tree.column("status", width=140, minwidth=80)

        vsb = ttk.Scrollbar(tree_frame, orient="vertical",
                             command=self._tree.yview)
        self._tree.configure(yscrollcommand=vsb.set)
        self._tree.pack(side="left", fill="both", expand=True)
        vsb.pack(side="right", fill="y")

        self._tree.bind("<<TreeviewSelect>>", self._on_cable_selected)

        # Selected cable display
        self._lbl_selected = tk.Label(frame, text="Selected: (none)",
                                       bg=C.BG, fg=C.ACCENT,
                                       font=("Consolas", 9, "bold"))
        self._lbl_selected.pack(padx=6, pady=(4, 2), anchor="w")

        # Buttons
        bf = tk.Frame(frame, bg=C.BG)
        bf.pack(fill="x", padx=6, pady=4)

        ttk.Button(bf, text="Auto-Finish Bundle", style="Green.TButton",
                   command=self._cmd_auto_finish
                   ).pack(side="left", fill="x", expand=True, padx=2)
        ttk.Button(bf, text="Clear All", style="Danger.TButton",
                   command=self._cmd_clear_all
                   ).pack(side="left", fill="x", expand=True, padx=2)

        bf2 = tk.Frame(frame, bg=C.BG)
        bf2.pack(fill="x", padx=6, pady=(0, 4))

        ttk.Button(bf2, text="Save Map", style="Accent.TButton",
                   command=self._cmd_save
                   ).pack(side="left", fill="x", expand=True, padx=2)
        ttk.Button(bf2, text="Load Map",
                   command=self._cmd_load
                   ).pack(side="left", fill="x", expand=True, padx=2)

    @staticmethod
    def _parse_filter(filt: str) -> dict:
        """Parse filter string into column constraints.

        Supports two styles:
          - Column filters:  "board:1,sub:2"  or  "group:6,number:18"
          - Plain text:      "1.3."  (matched against the full label)

        Returns dict with keys 'board','sub','group','number','text'.
        """
        result: Dict[str, str] = {}
        filt = filt.strip()
        if not filt:
            return result
        # Check for column:value syntax
        if ":" in filt:
            for part in filt.split(","):
                part = part.strip()
                if ":" in part:
                    key, _, val = part.partition(":")
                    key = key.strip().lower()
                    val = val.strip()
                    if key in ("board", "sub", "group", "number") and val:
                        result[key] = val
        else:
            result["text"] = filt.lower()
        return result

    def _refresh_cable_table(self):
        tree = self._tree
        tree.delete(*tree.get_children())
        filters = self._parse_filter(self._filter_var.get())
        show_connected = self._show_connected_var.get()

        cables = sorted(self.model.cables.values(), key=Cable.sort_key)
        for c in cables:
            connected = self.model.is_connected(c.label)
            if not show_connected and connected:
                continue
            # Apply filters
            if "text" in filters:
                if filters["text"] not in c.label.lower():
                    continue
            if "board" in filters:
                if str(c.board) != filters["board"]:
                    continue
            if "sub" in filters:
                if str(c.sub) != filters["sub"]:
                    continue
            if "group" in filters:
                if str(c.group) != filters["group"]:
                    continue
            if "number" in filters:
                if str(c.number) != filters["number"]:
                    continue

            if connected:
                cr, sl, ch = self.model.connections[c.label]
                status = f"{cr}:{sl}.{ch}"
            else:
                status = ""
            tree.insert("", "end", iid=c.label,
                        values=(c.board, c.sub, c.group, c.number, status))

    def _sort_tree(self, col: str):
        """Sort Treeview by column (toggle ascending/descending)."""
        items = [(self._tree.set(iid, col), iid)
                 for iid in self._tree.get_children()]
        if col in ("board", "sub", "group", "number"):
            items.sort(key=lambda x: int(x[0]) if x[0] else 999)
        else:
            items.sort(key=lambda x: x[0])
        for idx, (_, iid) in enumerate(items):
            self._tree.move(iid, "", idx)

    # ------------------------------------------------------------------
    #  Event handlers
    # ------------------------------------------------------------------

    def _on_crate_changed(self, _event):
        self._current_crate = self._crate_var.get()
        self._draw_crate()

    def _on_cable_selected(self, _event):
        sel = self._tree.selection()
        if sel:
            self._selected_cable = sel[0]
            c = self.model.cables.get(sel[0])
            if c:
                self._lbl_selected.configure(
                    text=f"Selected: {c.label}  (Board {c.board}, "
                         f"Sub {c.sub}, Group {c.group}, #{c.number})")
            else:
                self._lbl_selected.configure(text=f"Selected: {sel[0]}")
        else:
            self._selected_cable = None
            self._lbl_selected.configure(text="Selected: (none)")

    def _on_dot_click(self, slot: int, ch: int):
        crate = self._current_crate
        # If already connected, just show info
        existing = self.model.get_cable_at(crate, slot, ch)
        if existing:
            self._log(f"Channel {crate}:{slot}.{ch} already has cable "
                      f"{existing}", "warn")
            return
        # Need a selected cable
        if not self._selected_cable:
            self._log("No cable selected -- pick one from the table first",
                      "warn")
            return
        label = self._selected_cable
        if self.model.is_connected(label):
            self._log(f"Cable {label} is already connected", "warn")
            return

        if self.model.connect(label, crate, slot, ch):
            self._refresh_dot(slot, ch)
            self._last_connection = label
            self._log(f"Connected {label} -> {crate}:{slot}.{ch}", "ok")
            # Advance selection to next unconnected cable in same bundle
            self._advance_selection(label)
            self._refresh_cable_table()
            self._update_stats()
        else:
            self._log(f"Failed to connect {label} -> {crate}:{slot}.{ch}",
                      "error")

    def _on_dot_right_click(self, event, slot: int, ch: int):
        crate = self._current_crate
        label = self.model.get_cable_at(crate, slot, ch)
        if not label:
            return
        menu = tk.Menu(self.root, tearoff=0, bg=C.PANEL, fg=C.TEXT,
                       activebackground=C.ACCENT, activeforeground="white",
                       font=("Consolas", 9))
        menu.add_command(
            label=f"Disconnect {label}",
            command=lambda: self._disconnect(label, slot, ch))
        menu.tk_popup(event.x_root, event.y_root)

    def _on_dot_enter(self, event, slot: int, ch: int):
        crate = self._current_crate
        cable = self.model.get_cable_at(crate, slot, ch)
        text = f"Slot {slot}, Ch {ch}"
        if cable:
            text += f"\nCable: {cable}"
            c = self.model.cables.get(cable)
            if c:
                text += f"  (group {c.group})"
        self._tooltip.configure(text=text)
        # Position near the cursor but within canvas bounds
        tx = min(event.x + 15, self.CANVAS_W - 120)
        ty = max(event.y - 30, 5)
        self._tooltip.place(x=tx, y=ty)

    def _on_dot_leave(self):
        self._tooltip.place_forget()

    # ------------------------------------------------------------------
    #  Commands
    # ------------------------------------------------------------------

    def _disconnect(self, label: str, slot: int, ch: int):
        if self.model.disconnect(label):
            self._refresh_dot(slot, ch)
            self._log(f"Disconnected {label} from "
                      f"{self._current_crate}:{slot}.{ch}", "warn")
            self._refresh_cable_table()
            self._update_stats()

    def _cmd_auto_finish(self):
        if not self._last_connection:
            self._log("No recent connection to auto-finish from", "warn")
            return
        label = self._last_connection
        if not self.model.is_connected(label):
            self._log(f"Cable {label} is not connected", "warn")
            return
        results = self.model.auto_finish(label)
        if not results:
            self._log("Auto-finish: no more cables/channels to fill", "warn")
            return
        # Update dots if we are viewing the same crate
        crate = self.model.connections[label][0]
        for clbl, cr, sl, ch in results:
            if cr == self._current_crate:
                self._refresh_dot(sl, ch)
        self._last_connection = results[-1][0]
        self._log(f"Auto-finished {len(results)} cables in bundle "
                  f"{self.model.cables[label].bundle} on "
                  f"{crate}:slot {self.model.connections[label][1]}",
                  "ok")
        self._refresh_cable_table()
        self._update_stats()

    def _cmd_clear_all(self):
        n = len(self.model.connections)
        if n == 0:
            return
        if not messagebox.askyesno("Clear All",
                                    f"Disconnect all {n} cables?"):
            return
        labels = list(self.model.connections.keys())
        for label in labels:
            self.model.disconnect(label)
        self._last_connection = None
        self._draw_crate()
        self._refresh_cable_table()
        self._update_stats()
        self._log(f"Cleared all {n} connections", "warn")

    def _cmd_save(self):
        path = filedialog.asksaveasfilename(
            initialdir=_SCRIPT_DIR,
            defaultextension=".json",
            filetypes=[("JSON", "*.json"), ("All", "*.*")],
            initialfile="cable_connections.json")
        if not path:
            return
        self.model.save_json(path)
        self._log(f"Saved {len(self.model.connections)} connections to "
                  f"{os.path.basename(path)}", "ok")

    def _cmd_load(self):
        path = filedialog.askopenfilename(
            initialdir=_SCRIPT_DIR,
            filetypes=[("JSON", "*.json"), ("All", "*.*")])
        if not path:
            return
        loaded, skipped = self.model.load_json(path)
        self._draw_crate()
        self._refresh_cable_table()
        self._update_stats()
        msg = f"Loaded {loaded} connections from {os.path.basename(path)}"
        if skipped:
            msg += f" ({skipped} skipped)"
        self._log(msg, "ok")

    # ------------------------------------------------------------------
    #  Helpers
    # ------------------------------------------------------------------

    def _advance_selection(self, label: str):
        """Select the next unconnected cable in the same bundle."""
        cable = self.model.cables.get(label)
        if not cable:
            return
        bundle = cable.bundle
        candidates = sorted(
            [c for c in self.model.cables.values()
             if c.bundle == bundle
             and c.number > cable.number
             and not self.model.is_connected(c.label)],
            key=lambda c: c.number)
        if candidates:
            nxt = candidates[0].label
            self._selected_cable = nxt
            self._lbl_selected.configure(text=f"Selected: {nxt}")
            # Select in tree if visible
            if self._tree.exists(nxt):
                self._tree.selection_set(nxt)
                self._tree.see(nxt)

    def _update_stats(self):
        total_conn = len(self.model.connections)
        total_cables = len(self.model.cables)
        self._lbl_total.configure(
            text=f"Connected: {total_conn} / {total_cables}")
        conn, tot = self.model.crate_stats(self._current_crate)
        self._lbl_crate_stats.configure(
            text=f"{self._current_crate}:  {conn} / {tot} channels used")
        self._update_crate_stats()

    def _update_crate_stats(self):
        conn, tot = self.model.crate_stats(self._current_crate)
        self._lbl_crate_stats.configure(
            text=f"{self._current_crate}:  {conn} / {tot} channels used")

    def _log(self, msg: str, level: str = "info"):
        self._log_text.configure(state="normal")
        self._log_text.insert("end", msg + "\n", level)
        self._log_text.see("end")
        self._log_text.configure(state="disabled")


# ============================================================================
#  Entry point
# ============================================================================

def main():
    parser = argparse.ArgumentParser(
        description="HyCal Cable Manager -- map cables to FADC250 channels")
    parser.add_argument("--cables", default=DEFAULT_CABLE_CSV,
                        help="Path to cable_group_numbers.csv")
    parser.add_argument("--boards", default=DEFAULT_BOARD_CONFIG,
                        help="Path to board_config.json")
    parser.add_argument("--connections",
                        help="Load existing connection map to resume")
    args = parser.parse_args()

    model = CablingModel()

    # Load board config
    if os.path.exists(args.boards):
        model.load_board_config(args.boards)
        print(f"Board config: {args.boards}")
    else:
        print(f"WARNING: board config not found: {args.boards}")
        print("Using default: slots 3-10, 13-20 for all crates")
        default_slots = [3, 4, 5, 6, 7, 8, 9, 10,
                         13, 14, 15, 16, 17, 18, 19, 20]
        model.board_config = {name: list(default_slots)
                              for name in CRATE_NAMES}

    # Load cables
    if os.path.exists(args.cables):
        n = model.load_cables_csv(args.cables)
        print(f"Loaded {n} cables from {args.cables}")
    else:
        print(f"ERROR: cable CSV not found: {args.cables}")
        return

    # Load existing connections
    if args.connections and os.path.exists(args.connections):
        loaded, skipped = model.load_json(args.connections)
        print(f"Loaded {loaded} existing connections"
              f" ({skipped} skipped)" if skipped else
              f"Loaded {loaded} existing connections")

    root = tk.Tk()
    CablingGUI(root, model)
    root.mainloop()


if __name__ == "__main__":
    main()
