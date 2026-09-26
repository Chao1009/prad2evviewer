#!/usr/bin/env python3
"""
Scintillator–HyCal Coincidence Monitor
=======================================
Connects to a running prad2_server (HTTP REST API, port 5051 by default),
iterates through every event in the loaded file, and accumulates per-module
coincidence statistics between the four upstream scintillators (V1–V4) and
every HyCal module.

Two event-selection modes are provided:

  AND mode — rate(V_i, M) = N(V_i fired AND M is best) / N(M is best)
             Only events where at least one scintillator fired are counted.
             Gives P(V_i fired | M was the highest-ADC module).

  OR  mode — rate(V_i, M) = N(V_i fired AND M is best) / N(V_i fired)
             All HyCal cluster events counted regardless of scintillator.
             Gives P(M is best | V_i fired).

A channel "fired" when it has at least one FADC peak whose height (above
pedestal) exceeds the user-specified threshold.

The bottom half of the window shows individual waveforms: all four
scintillators (V1–V4) overlaid and the HyCal module last clicked on the map,
both fetched from the server for the event number entered in the Event Browser.

Usage
-----
    python scripts/coincidence_monitor.py [--url http://HOST:PORT]
                                          [--theme dark|light]
"""
from __future__ import annotations

import argparse
import json as json_mod
import math
import sys
import threading
import time
import multiprocessing
import os
from concurrent.futures import ThreadPoolExecutor, as_completed
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Dict, List, Optional

import urllib.request
import urllib.error

import numpy as np

from PyQt6.QtWidgets import (
    QApplication, QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QGridLayout,
    QPushButton, QLabel, QDoubleSpinBox, QLineEdit, QProgressBar,
    QSplitter, QSizePolicy, QButtonGroup, QRadioButton, QGroupBox,
    QSpinBox, QFrame, QMessageBox, QCheckBox, QFileDialog, QScrollArea,
    QListWidget, QAbstractItemView,
)
from PyQt6.QtCore import Qt, QThread, QTimer, pyqtSignal, QRectF, QPointF
from PyQt6.QtGui import (
    QColor, QFont, QPen, QPainter, QPolygonF,
)

from hycal_geoview import (
    load_modules, load_daq_map, load_roc_tag_map, hole_ring,
    HyCalMapWidget, cmap_qcolor, series_qcolor, draw_wave_axes,
    apply_theme_palette, set_theme, available_themes, THEME, themed,
)
from prad2_env import import_prad2py

_prad2py, _ = import_prad2py(build_first=False)
_HAVE_PRAD2PY = _prad2py is not None


# Paths & constants

SCRIPT_DIR = Path(__file__).resolve().parent
DB_DIR = SCRIPT_DIR / ".." / "database"
MODULES_JSON  = DB_DIR / "hycal_map.json"
DAQ_CFG_JSON  = DB_DIR / "daq_config.json"

DEFAULT_URL = "http://localhost:5051"

# Scintillator names (channel keys resolved at runtime from daq_map + daq_config)
SCINTILLATORS = ("V1", "V2", "V3", "V4")

# Default thresholds (FADC peak height above pedestal, ADC counts)
DEFAULT_SCINT_THR = 500.0
DEFAULT_HYCAL_THR = 500.0

# Default signal time windows (ns); scintillators and HyCal share the FADC window
DEFAULT_SCINT_TMIN = 160.0
DEFAULT_SCINT_TMAX = 200.0
DEFAULT_HYCAL_TMIN = 160.0
DEFAULT_HYCAL_TMAX = 200.0

# Default max number of local ADC maxima allowed in HyCal (1 = single cluster only)
DEFAULT_MAX_LOCAL_MAXIMA = 1

# Event-selection modes
MODE_AND = "AND"   # rate = N(Vi AND M_best) / N(M_best)  — P(Vi fired | M is best)
MODE_OR  = "OR"    # rate = N(Vi AND M_best) / N(Vi fired) — P(M is best | Vi fired)

# Skip calibration/background events by trigger bit OR by module occupancy.
LMS_TRIGGER_BIT   = 1 << 24   # bit 24 — LMS light source
ALPHA_TRIGGER_BIT = 1 << 25   # bit 25 — alpha source
SKIP_TRIGGER_MASK = LMS_TRIGGER_BIT | ALPHA_TRIGGER_BIT
LMS_MAX_MODULES   = 1000      # more than this many modules above threshold → LMS

# Cluster-size cut: events with fewer than this many HyCal modules above threshold
# are discarded as isolated electronic noise / discharge artefacts.
DEFAULT_MIN_CLUSTER_MODS = 2

# Top-level display modes
DISPLAY_COINC   = "coinc"    # coincidence rate / occupancy statistics
DISPLAY_INSTANT = "instant"  # live event-by-event ADC display

# Map display modes (used within DISPLAY_COINC)
VIEW_COINC   = "coincidence"  # colour = coincidence rate with selected scintillator
VIEW_OCC     = "occupancy"    # colour = number of events module fired above threshold
VIEW_INSTANT = "instant"      # colour = per-module ADC signal for current event

# Veto scintillator motor PV prefixes (suffixes: see VetoMotorController)
VETO_PV_BASES: Dict[str, str] = {
    "V1": "prad:veto1",
    "V2": "prad:veto2",
    "V3": "prad:veto3",
    "V4": "prad:veto4",
}
VETO_POLL_MS = 500   # 2 Hz

N_WORKERS  = 16    # parallel HTTP workers for event scanning
BATCH_SIZE = 64    # events per processing batch (small keeps in-flight JSON bounded)

CLK_MHZ = 250.0    # FADC clock (for x-axis in ns)


# Veto motor position reader

class VetoMotorController:
    """Reads and writes position PVs for the four veto scintillator motors.

    PVs per motor:
        prad:vetoN.VAL  — setpoint (write to command a move)
        prad:vetoN.RBV  — actual read-back position
        prad:vetoN.MOVN — 1 while moving, 0 when at rest

    Degrades gracefully when pyepics is not installed.
    """

    def __init__(self) -> None:
        self._pvs: Dict[str, Any] = {}
        self._epics_ok = False
        try:
            import epics as _epics
            for vname, base in VETO_PV_BASES.items():
                self._pvs[f"{vname}_val"]  = _epics.PV(f"{base}.VAL")
                self._pvs[f"{vname}_rbv"]  = _epics.PV(f"{base}.RBV")
                self._pvs[f"{vname}_movn"] = _epics.PV(f"{base}.MOVN")
            self._epics_ok = True
        except ImportError:
            pass

    def get(self, key: str) -> Optional[float]:
        pv = self._pvs.get(key)
        if pv is None or not pv.connected:
            return None
        v = pv.get()
        return float(v) if v is not None else None

    def put(self, key: str, value: float) -> bool:
        """Write value to a VAL PV. Returns True on success."""
        pv = self._pvs.get(key)
        if pv is None or not pv.connected:
            return False
        pv.put(value)
        return True

    @property
    def available(self) -> bool:
        return self._epics_ok


# HTTP helpers

def _http_get(url: str, timeout: float = 5.0) -> Optional[dict]:
    try:
        with urllib.request.urlopen(url, timeout=timeout) as resp:
            return json_mod.loads(resp.read())
    except Exception:
        return None


def _build_neighbor_map(modules) -> Dict[str, frozenset]:
    """Return name→frozenset-of-neighbor-names for HyCal modules.

    Two modules are neighbors when they share an edge (touching sides, not just
    a corner).  Works for mixed W/PbGl geometries because it checks actual
    center-distance against the sum of half-widths with a 1 mm tolerance.
    """
    phys = [(m.name, m.x, m.y, m.sx, m.sy)
            for m in modules if m.mod_type != "LMS"]
    neighbors: Dict[str, List[str]] = {name: [] for name, *_ in phys}
    tol = 1.0   # mm
    for i, (n1, x1, y1, sx1, sy1) in enumerate(phys):
        for n2, x2, y2, sx2, sy2 in phys[i + 1:]:
            dx = abs(x1 - x2)
            dy = abs(y1 - y2)
            share_v = (abs(dx - (sx1 + sx2) / 2) < tol
                       and dy < min(sy1, sy2) / 2 + tol)
            share_h = (abs(dy - (sy1 + sy2) / 2) < tol
                       and dx < min(sx1, sx2) / 2 + tol)
            if share_v or share_h:
                neighbors[n1].append(n2)
                neighbors[n2].append(n1)
    return {name: frozenset(nbs) for name, nbs in neighbors.items()}


def _count_local_maxima(module_adc: Dict[str, float],
                        neighbor_map: Dict[str, frozenset]) -> int:
    """Count HyCal modules that are strictly higher than all their neighbors.

    Neighbors absent from module_adc are treated as ADC = 0 (below threshold).
    """
    count = 0
    for name, adc in module_adc.items():
        if all(adc > module_adc.get(nb, 0.0) for nb in neighbor_map.get(name, ())):
            count += 1
    return count


def _build_tuple_role_map(scint_keys: Dict[str, str],
                          module_keys: Dict[str, str]) -> Dict[tuple, tuple]:
    """Build (roc_tag, slot, chan) → (kind, name) lookup.

    kind=0 → scintillator, kind=1 → HyCal module.
    Avoids f-string formatting in the per-channel hot loop.
    """
    role: Dict[tuple, tuple] = {}
    for sname, key in scint_keys.items():
        roc, s, c = key.split("_")
        role[(int(roc), int(s), int(c))] = (0, sname)
    for mname, key in module_keys.items():
        roc, s, c = key.split("_")
        role[(int(roc), int(s), int(c))] = (1, mname)
    return role


# Event selection & counting

@dataclass(frozen=True)
class CoincCuts:
    """Event-selection cuts of a coincidence run (picklable, so the parallel
    file workers get it as is)."""
    scint_thr: float
    hycal_thr: float
    mode: str = MODE_AND
    min_mods: int = DEFAULT_MIN_CLUSTER_MODS
    scint_t_min: float = -math.inf
    scint_t_max: float = math.inf
    hycal_t_min: float = -math.inf
    hycal_t_max: float = math.inf
    neighbor_map: Dict[str, frozenset] = field(default_factory=dict)
    max_lm: int = 0   # max HyCal local maxima; 0 = no cut

    @property
    def require_scint(self) -> bool:
        """Events without a fired scintillator are dropped in AND mode, and
        in either mode while the scint time cut is active."""
        return self.mode == MODE_AND or self.scint_t_min > -math.inf


def _select_best(ma: Dict[str, float], cuts: CoincCuts) -> Optional[str]:
    """The highest-ADC module of an event with module ADCs ``ma`` (modules
    above threshold), or None when the event fails the LMS-occupancy,
    cluster-size or local-maxima cut."""
    n_ma = len(ma)
    if n_ma > LMS_MAX_MODULES or n_ma < cuts.min_mods:
        return None
    # Single-module clusters are trivially 1 local max — skip the scan.
    if (cuts.max_lm > 0 and n_ma > 1
            and _count_local_maxima(ma, cuts.neighbor_map) > cuts.max_lm):
        return None
    return max(ma, key=ma.__getitem__)


class Counters:
    """Coincidence counts of one worker.  Callers count every event in
    ``processed`` themselves; add() counts the events that passed the cuts.
    """

    def __init__(self, module_names) -> None:
        self._mod_names = list(module_names)
        self._reset()

    def _reset(self) -> None:
        mods = self._mod_names
        self.module_hits:    Dict[str, int] = {m: 0 for m in mods}
        self.scint_hits:     Dict[str, int] = {s: 0 for s in SCINTILLATORS}
        self.scint_hits_any: Dict[str, int] = {s: 0 for s in SCINTILLATORS}
        self.coincidences: Dict[str, Dict[str, int]] = {
            s: {m: 0 for m in mods} for s in SCINTILLATORS
        }
        self.processed = 0

    def add(self, sf: Dict[str, bool], best: str, require_scint: bool) -> bool:
        """Count an event whose best module is ``best``; ``sf`` maps
        scintillator → fired.

        Individual scintillator fires are counted first, so the stats panel
        always shows per-scintillator rates.  The event itself is dropped
        (returns False) when ``require_scint`` and no scintillator fired.
        Otherwise module_hits[best] += 1 and, for each fired Vi,
        scint_hits[Vi] and coincidences[Vi][best] += 1.
        """
        for sname, fired in sf.items():
            if fired:
                self.scint_hits_any[sname] += 1
        if require_scint and not any(sf.values()):
            return False
        self.module_hits[best] += 1
        for sname, fired in sf.items():
            if fired:
                self.scint_hits[sname] += 1
                self.coincidences[sname][best] += 1
        return True

    def snapshot(self, mode: str) -> dict:
        rates: Dict[str, Dict[str, float]] = {}
        for sname in SCINTILLATORS:
            rates[sname] = {}
            s_denom = self.scint_hits[sname]   # N(V_i fired AND HyCal cluster)
            for mname, m_denom in self.module_hits.items():
                ncoinc = self.coincidences[sname][mname]
                if mode == MODE_AND:
                    rates[sname][mname] = (ncoinc / m_denom
                                           if m_denom > 0 else math.nan)
                else:
                    rates[sname][mname] = (ncoinc / s_denom
                                           if s_denom > 0 else math.nan)
        return {
            "rates": rates,
            "mode": mode,
            "module_hits": dict(self.module_hits),
            "scint_hits": dict(self.scint_hits),
            "scint_hits_any": dict(self.scint_hits_any),
            "processed": self.processed,
        }

    def take_delta(self) -> dict:
        """The counts as plain dicts (for another process); counting
        restarts from zero in new dicts, so the returned ones stay as they
        are while a Queue feeder thread pickles them."""
        delta = {
            "module_hits":    self.module_hits,
            "scint_hits":     self.scint_hits,
            "scint_hits_any": self.scint_hits_any,
            "coincidences":   self.coincidences,
            "processed":      self.processed,
        }
        self._reset()
        return delta

    def merge(self, delta: dict) -> None:
        """Add a take_delta() result."""
        for m, v in delta["module_hits"].items():
            if v: self.module_hits[m] += v
        for s, v in delta["scint_hits"].items():
            if v: self.scint_hits[s] += v
        for s, v in delta["scint_hits_any"].items():
            if v: self.scint_hits_any[s] += v
        for s, ccol in delta["coincidences"].items():
            row = self.coincidences[s]
            for m, v in ccol.items():
                if v: row[m] += v
        self.processed += delta["processed"]


# Parallel local-EVIO worker (subprocess)

# Worker-process global for the update queue.  ``forkserver``/``spawn``
# refuse to pickle a Queue in apply_async args ("Queue objects should
# only be shared between processes through inheritance"), so we install
# it once per worker via Pool(initializer=...).
_WORKER_UPDATE_Q = None


def _init_worker(q):
    global _WORKER_UPDATE_Q
    _WORKER_UPDATE_Q = q


def _process_files(args: dict) -> None:
    """Process a list of EVIO files in a subprocess.

    Periodically pushes incremental count deltas to the worker's update
    queue (see _init_worker) so the parent QThread can update the UI while
    subprocesses are still running.  Returns None — all results flow
    through the queue.
    """
    paths           = args["paths"]
    tuple_to_role   = args["tuple_to_role"]
    cuts: CoincCuts = args["cuts"]
    worker_id       = args["worker_id"]
    n_files_in_chunk = len(paths)
    update_q        = _WORKER_UPDATE_Q
    push_interval   = args.get("push_interval", 0.5)
    scint_thr       = cuts.scint_thr
    hycal_thr       = cuts.hycal_thr

    # Counts since the last push to the parent.
    delta = Counters(args["mod_names"])

    # Per-file progress state — included in every push so the parent UI can
    # render one progress bar per worker.
    cur_file_idx     = -1
    cur_basename     = ""
    cur_records_done = 0
    cur_records_tot  = 0
    finished_chunk   = False

    # Split the channel list into a tiny scint table (visited first so we can
    # short-circuit non-coincidence events) and a HyCal table.  The ROC tag
    # is bucket-indexed up front to save dict.get() overhead per channel.
    scint_locs: List[tuple] = []   # [(roc_tag, slot, chan, name), ...]
    hycal_locs_by_roc: Dict[int, List[tuple]] = {}
    for (rt, sn, cn), (kind, name) in tuple_to_role.items():
        if kind == 0:
            scint_locs.append((rt, sn, cn, name))
        else:
            hycal_locs_by_roc.setdefault(rt, []).append((sn, cn, name))
    # Sort each ROC's HyCal channel list by (slot, chan) so consecutive
    # iterations stay on the same slot — lets us cache slot.channel() lookups.
    for _lst in hycal_locs_by_roc.values():
        _lst.sort()

    # With require_scint an event whose scintillators don't fire is never
    # counted, so skipping HyCal in that case is safe and saves all
    # per-channel work for those events.
    require_scint = cuts.require_scint

    def _flush():
        if update_q is None:
            return
        update_q.put({
            **delta.take_delta(),
            "worker_id":       worker_id,
            "n_files":         n_files_in_chunk,
            "file_idx":        cur_file_idx,
            "file_basename":   cur_basename,
            "records_done":    cur_records_done,
            "records_total":   cur_records_tot,
            "finished":        finished_chunk,
        })

    try:
        import prad2py as _p2       # type: ignore
    except ImportError as e:
        sys.stderr.write(f"[worker {worker_id} pid={os.getpid()}] "
                         f"prad2py import failed: {e}\n")
        sys.stderr.flush()
        finished_chunk = True
        _flush()
        return None

    dec      = _p2.dec
    cfg      = dec.load_daq_config()
    ch       = dec.EvChannel()
    ch.set_config(cfg)
    # WaveAnalyzer not needed in the parallel worker — the per-event hot
    # loop uses a numpy-only peak-in-window heuristic instead.

    # Sample-index bounds for the time-window peak search, on the
    # WaveAnalyzer clock.  N_PED is the number of leading samples used to
    # estimate pedestal — chosen large enough to be statistically stable but
    # small enough never to overlap a physics signal (those arrive after
    # ~100 ns at the earliest).
    N_PED = 30
    clk_mhz = cfg.wave_cfg.clk_mhz
    ns_per_sample = 1e3 / (clk_mhz if clk_mhz > 0 else CLK_MHZ)

    def _sample_window(t_min: float, t_max: float) -> tuple:
        lo = (max(N_PED, int(t_min / ns_per_sample))
              if t_min > -math.inf else N_PED)
        hi = int(t_max / ns_per_sample) + 1 if t_max < math.inf else 1 << 30
        return lo, hi

    s_lo, s_hi = _sample_window(cuts.scint_t_min, cuts.scint_t_max)
    h_lo, h_hi = _sample_window(cuts.hycal_t_min, cuts.hycal_t_max)

    last_push = time.monotonic()

    for fidx, path in enumerate(paths):
        cur_file_idx     = fidx
        cur_basename     = os.path.basename(path)
        cur_records_done = 0
        cur_records_tot  = 0

        st = ch.open_auto(path)
        if st != dec.Status.success:
            sys.stderr.write(f"[worker {worker_id} pid={os.getpid()}] "
                             f"open_auto({path}) -> {st}\n")
            sys.stderr.flush()
            ch.close()
            _flush()   # let parent know this file was skipped
            continue

        if ch.is_random_access():
            cur_records_tot = ch.get_random_access_event_count()

        # Push initial progress so the bar appears immediately.
        _flush()
        last_push = time.monotonic()

        while True:
            if ch.read() != dec.Status.success:
                break
            cur_records_done += 1
            if not ch.scan():
                continue
            if ch.get_event_type() != dec.EventType.Physics:
                continue

            for si in range(ch.get_n_events()):
                ch.select_event(si)
                info = ch.info()
                delta.processed += 1

                if int(info.trigger_bits) & SKIP_TRIGGER_MASK:
                    continue

                fadc_evt = ch.fadc()
                sf: Dict[str, bool]  = {s: False for s in SCINTILLATORS}
                ma: Dict[str, float] = {}

                # Build a roc_tag -> roc-object map once per event so the
                # scint and HyCal passes can both look up ROCs without
                # iterating fadc_evt twice.
                roc_by_tag: Dict[int, object] = {}
                for r in range(fadc_evt.nrocs):
                    rr = fadc_evt.roc(r)
                    roc_by_tag[int(rr.tag)] = rr

                scint_any = False
                for rt, sn, cn, name in scint_locs:
                    rr = roc_by_tag.get(rt)
                    if rr is None:
                        continue
                    chan = rr.slot(sn).channel(cn)
                    if chan.nsamples < N_PED + 1:
                        continue
                    samples = chan.samples
                    if samples.max() - samples.min() < scint_thr:
                        continue
                    ped    = samples[:N_PED].mean()
                    window = samples[s_lo:s_hi]
                    if window.size == 0:
                        continue
                    if window.max() - ped > scint_thr:
                        sf[name]  = True
                        scint_any = True

                if require_scint and not scint_any:
                    continue

                for rt, hy_locs in hycal_locs_by_roc.items():
                    rr = roc_by_tag.get(rt)
                    if rr is None:
                        continue
                    last_slot_num = -1
                    slot = None
                    for sn, cn, name in hy_locs:
                        if sn != last_slot_num:
                            slot = rr.slot(sn)
                            last_slot_num = sn
                        chan = slot.channel(cn)
                        if chan.nsamples < N_PED + 1:
                            continue
                        samples = chan.samples
                        if samples.max() - samples.min() < hycal_thr:
                            continue
                        ped    = samples[:N_PED].mean()
                        window = samples[h_lo:h_hi]
                        if window.size == 0:
                            continue
                        height = float(window.max() - ped)
                        if height > hycal_thr:
                            ma[name] = height

                best = _select_best(ma, cuts)
                if best is not None:
                    delta.add(sf, best, require_scint)

            # Time-based partial push (between records, not per event,
            # to keep time.monotonic() out of the hottest inner loop).
            now = time.monotonic()
            if now - last_push >= push_interval:
                _flush()
                last_push = now

        # Final push for this file so the parent sees the 100% mark.
        _flush()
        ch.close()

    finished_chunk = True
    _flush()
    return None


# Waveform collector

class WaveformCollector:
    """Thread-safe accumulator of per-event waveform records for coincidences.

    For each qualifying event, stores ADC samples from the fired scintillator
    and the HyCal module with the highest ADC.  Samples are read from the
    inline 's' field if present (ET mode); otherwise a separate
    /api/waveform/{ev}/{key} request is made (file mode).

    Call save(path) to write a compressed .npz file when done.
    """

    def __init__(self, server_url: str, max_records: int) -> None:
        self._url  = server_url
        self._max  = max_records
        self._buf: List[dict] = []
        self._lock = threading.Lock()

    @property
    def count(self) -> int:
        with self._lock:
            return len(self._buf)

    @property
    def full(self) -> bool:
        with self._lock:
            return len(self._buf) >= self._max

    def _get_samples(self, ev_num: int, ch_data: dict, key: str) -> List[int]:
        sp = ch_data.get("s", [])
        if sp:
            return list(sp)
        wf = _http_get(f"{self._url}/api/waveform/{ev_num}/{key}")
        return wf.get("s", []) if wf and "error" not in wf else []

    def collect(self, data: dict,
                scint_name: str, scint_key: str,
                hycal_name: str, hycal_key: str) -> bool:
        """Try to add one record.  Returns True if a record was stored."""
        with self._lock:
            if len(self._buf) >= self._max:
                return False
        channels  = data.get("channels", {})
        ev_num    = data.get("event_number", data.get("event", 0))
        scint_ch  = channels.get(scint_key, {})
        hycal_ch  = channels.get(hycal_key, {})
        scint_sp  = self._get_samples(ev_num, scint_ch, scint_key)
        hycal_sp  = self._get_samples(ev_num, hycal_ch, hycal_key)
        scint_int = max((pk.get("i", 0.0) for pk in scint_ch.get("pk", [])), default=0.0)
        hycal_int = max((pk.get("i", 0.0) for pk in hycal_ch.get("pk", [])), default=0.0)
        with self._lock:
            if len(self._buf) >= self._max:
                return False
            self._buf.append({
                "event":         ev_num,
                "scint_name":    scint_name,
                "scint_key":     scint_key,
                "scint_int":     scint_int,
                "scint_pm":      float(scint_ch.get("pm", 0) or 0),
                "scint_samples": scint_sp,
                "hycal_name":    hycal_name,
                "hycal_key":     hycal_key,
                "hycal_int":     hycal_int,
                "hycal_pm":      float(hycal_ch.get("pm", 0) or 0),
                "hycal_samples": hycal_sp,
            })
            return True

    def save(self, path: Path) -> int:
        """Write accumulated records to a compressed .npz file.
        Returns the number of records written."""
        with self._lock:
            buf = list(self._buf)
        if not buf:
            return 0
        path.parent.mkdir(parents=True, exist_ok=True)
        n = len(buf)
        sp_len = max(
            max((len(e["scint_samples"]) for e in buf), default=0),
            max((len(e["hycal_samples"]) for e in buf), default=0),
            1,
        )
        scint_sp = np.zeros((n, sp_len), dtype=np.int16)
        hycal_sp = np.zeros((n, sp_len), dtype=np.int16)
        for i, e in enumerate(buf):
            ss, hs = e["scint_samples"], e["hycal_samples"]
            scint_sp[i, :len(ss)] = ss
            hycal_sp[i, :len(hs)] = hs
        np.savez_compressed(
            str(path),
            event_numbers   = np.array([e["event"]      for e in buf], dtype=np.int64),
            scint_names     = np.array([e["scint_name"]  for e in buf]),
            scint_keys      = np.array([e["scint_key"]   for e in buf]),
            scint_integrals = np.array([e["scint_int"]   for e in buf], dtype=np.float32),
            scint_ped_means = np.array([e["scint_pm"]    for e in buf], dtype=np.float32),
            scint_samples   = scint_sp,
            hycal_names     = np.array([e["hycal_name"]  for e in buf]),
            hycal_keys      = np.array([e["hycal_key"]   for e in buf]),
            hycal_integrals = np.array([e["hycal_int"]   for e in buf], dtype=np.float32),
            hycal_ped_means = np.array([e["hycal_pm"]    for e in buf], dtype=np.float32),
            hycal_samples   = hycal_sp,
        )
        return n


# Coincidence scan worker

def _fetch_event(server_url: str, ev: int) -> Optional[dict]:
    return _http_get(f"{server_url}/api/event/{ev}")


def _channel_fired(ch_data: dict, threshold: float,
                   t_min: float = -math.inf, t_max: float = math.inf) -> bool:
    """True if any peak has height > threshold and time within [t_min, t_max] ns."""
    for pk in ch_data.get("pk", []):
        if pk.get("h", 0.0) > threshold and t_min <= pk.get("t", 0.0) <= t_max:
            return True
    return False


def _collect_wfm(collector: "WaveformCollector", data: dict,
                 scint_fired: Dict[str, bool],
                 scint_keys: Dict[str, str],
                 hycal_name: str, hycal_key: str) -> None:
    """Pick the highest-ADC fired scintillator and hand off to the collector."""
    channels = data.get("channels", {})
    best_sname, best_skey, best_sint = None, None, -1.0
    for sname, skey in scint_keys.items():
        if not scint_fired.get(sname):
            continue
        integral = max(
            (pk.get("i", 0.0) for pk in channels.get(skey, {}).get("pk", [])),
            default=0.0,
        )
        if integral > best_sint:
            best_sint, best_sname, best_skey = integral, sname, skey
    if best_sname is None:
        return
    collector.collect(data,
                      scint_name=best_sname, scint_key=best_skey,
                      hycal_name=hycal_name, hycal_key=hycal_key)


class _CoincWorker(QThread):
    """Common part of the coincidence workers: channel keys, cuts, optional
    waveform collection and the progress / stats / finished reporting."""
    progress         = pyqtSignal(int, int)   # (done, total); total -1 = ET, 0 = unknown
    stats_update     = pyqtSignal(dict)
    finished         = pyqtSignal(str)
    waveforms_saved  = pyqtSignal(int, str)   # (count, file_path)

    def __init__(self, module_keys: Dict[str, str],
                 scint_keys: Dict[str, str],
                 cuts: CoincCuts,
                 wfm_collector: Optional["WaveformCollector"] = None,
                 wfm_save_path: Optional[Path] = None,
                 parent=None):
        super().__init__(parent)
        self._mod_keys   = module_keys
        self._scint_keys = scint_keys
        self._cuts       = cuts
        self._wfm_coll   = wfm_collector
        self._wfm_path   = wfm_save_path
        self._stop_evt   = threading.Event()

    def stop(self):
        self._stop_evt.set()

    def _emit(self, counters: Counters, done: int, total: int) -> None:
        self.progress.emit(done, total)
        self.stats_update.emit(counters.snapshot(self._cuts.mode))

    def _finish(self, counters: Counters, done: int, total: int,
                err: str = "") -> None:
        """Last progress/stats update, the waveform file, then finished()
        with ``err``, or 'stopped' / ''."""
        self._emit(counters, done, total)
        if self._wfm_coll and self._wfm_path and self._wfm_coll.count > 0:
            n = self._wfm_coll.save(self._wfm_path)
            self.waveforms_saved.emit(n, str(self._wfm_path))
        self.finished.emit(err or ("stopped" if self._stop_evt.is_set() else ""))

    def _count_server_event(self, data: dict, counters: Counters) -> None:
        """Select and count one /api/event JSON (server WaveAnalyzer peaks)."""
        # Reject LMS / alpha events by trigger bit (fast path)
        if data.get("trigger_bits", 0) & SKIP_TRIGGER_MASK:
            return
        cuts = self._cuts
        channels = data.get("channels", {})
        scint_fired = {
            sname: (skey in channels
                    and _channel_fired(channels[skey], cuts.scint_thr,
                                       cuts.scint_t_min, cuts.scint_t_max))
            for sname, skey in self._scint_keys.items()
        }

        # Per-module ADC: max peak height within the time cut.
        module_adc: Dict[str, float] = {}
        for mname, mkey in self._mod_keys.items():
            if mkey in channels:
                peaks = channels[mkey].get("pk", [])
                adc = max(
                    (float(pk.get("h", 0.0)) for pk in peaks
                     if cuts.hycal_t_min <= pk.get("t", 0.0) <= cuts.hycal_t_max),
                    default=0.0)
                if adc > cuts.hycal_thr:
                    module_adc[mname] = adc

        best = _select_best(module_adc, cuts)
        if best is None or not counters.add(scint_fired, best,
                                            cuts.require_scint):
            return
        if self._wfm_coll and not self._wfm_coll.full:
            _collect_wfm(self._wfm_coll, data, scint_fired,
                         self._scint_keys, best, self._mod_keys[best])


class ProcessWorker(_CoincWorker):
    """Scans events 1..n_events of the file loaded in the server."""

    def __init__(self, server_url: str, n_events: int,
                 module_keys: Dict[str, str],
                 scint_keys: Dict[str, str],
                 cuts: CoincCuts,
                 wfm_collector: Optional["WaveformCollector"] = None,
                 wfm_save_path: Optional[Path] = None,
                 parent=None):
        super().__init__(module_keys, scint_keys, cuts,
                         wfm_collector, wfm_save_path, parent)
        self._url = server_url
        self._n   = n_events

    def run(self):
        counters  = Counters(self._mod_keys)
        last_emit = time.monotonic()

        with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
            batch_start = 1
            while batch_start <= self._n and not self._stop_evt.is_set():
                batch_end = min(batch_start + BATCH_SIZE - 1, self._n)
                futures = {
                    pool.submit(_fetch_event, self._url, ev): ev
                    for ev in range(batch_start, batch_end + 1)
                }

                for fut in as_completed(futures):
                    if self._stop_evt.is_set():
                        break
                    data = fut.result()
                    counters.processed += 1
                    if data and "error" not in data:
                        self._count_server_event(data, counters)

                batch_start = batch_end + 1
                now = time.monotonic()
                if now - last_emit > 0.5:
                    self._emit(counters, counters.processed, self._n)
                    last_emit = now

        self._finish(counters, counters.processed, self._n)


class ProcessWorkerET(_CoincWorker):
    """Accumulates coincidence statistics from live ET events.

    Polls /api/ring every POLL_INTERVAL, fetches each new sequence number
    from the ring buffer, and processes it exactly once.  Runs until stopped.
    """

    POLL_INTERVAL = 0.05  # seconds between /api/ring polls (20 Hz)

    def __init__(self, server_url: str,
                 module_keys: Dict[str, str],
                 scint_keys: Dict[str, str],
                 cuts: CoincCuts,
                 max_rate_hz: float = 0.0,
                 wfm_collector: Optional["WaveformCollector"] = None,
                 wfm_save_path: Optional[Path] = None,
                 parent=None):
        super().__init__(module_keys, scint_keys, cuts,
                         wfm_collector, wfm_save_path, parent)
        self._url      = server_url
        self._max_rate = max_rate_hz

    def run(self):
        counters  = Counters(self._mod_keys)
        last_seq  = 0    # max sequence number seen
        last_emit = time.monotonic()

        with ThreadPoolExecutor(max_workers=N_WORKERS) as pool:
            while not self._stop_evt.is_set():
                batch_t0  = time.monotonic()
                ring_data = _http_get(f"{self._url}/api/ring", timeout=2.0)
                if ring_data is None:
                    time.sleep(self.POLL_INTERVAL)
                    continue

                # Sequence numbers increase monotonically, so anything above
                # last_seq has not been processed yet.
                new_seqs = sorted(
                    s for s in ring_data.get("ring", []) if s > last_seq)

                if new_seqs:
                    futures = {
                        pool.submit(_fetch_event, self._url, seq): seq
                        for seq in new_seqs
                        if not self._stop_evt.is_set()
                    }
                    for fut in as_completed(futures):
                        if self._stop_evt.is_set():
                            break
                        data = fut.result()
                        counters.processed += 1
                        if data and "error" not in data:
                            self._count_server_event(data, counters)

                    last_seq = max(new_seqs)

                now = time.monotonic()
                if now - last_emit > 0.5:
                    self._emit(counters, counters.processed, -1)
                    last_emit = now

                # Rate limiting: sleep at least POLL_INTERVAL; sleep longer
                # if the user set a max rate and the batch finished too fast.
                batch_elapsed = time.monotonic() - batch_t0
                if self._max_rate > 0 and new_seqs:
                    target = len(new_seqs) / self._max_rate
                    sleep_t = max(self.POLL_INTERVAL, target - batch_elapsed)
                else:
                    sleep_t = self.POLL_INTERVAL
                time.sleep(sleep_t)

        self._finish(counters, counters.processed, -1)


# Local EVIO worker  (no server — uses prad2py directly)

class ProcessWorkerLocal(_CoincWorker):
    """Reads an EVIO file directly via prad2py and accumulates coincidence
    statistics.  Requires prad2py to be installed (_HAVE_PRAD2PY == True).

    Optimised for throughput: a (roc_tag, slot, chan) lookup means only
    channels that belong to a known scintillator or HyCal module are
    analysed; raw waveform samples are never copied to Python lists in the
    stats path.
    """
    coincidence_event = pyqtSignal(dict)       # step-through: per-event waveform data
    worker_progress   = pyqtSignal(dict)       # parallel mode: per-worker file progress
    workers_setup     = pyqtSignal(int)        # parallel mode: number of workers about to run

    def __init__(self, evio_paths: List[str],
                 module_keys: Dict[str, str],
                 scint_keys: Dict[str, str],
                 cuts: CoincCuts,
                 wfm_collector: Optional[WaveformCollector] = None,
                 wfm_save_path: Optional[Path] = None,
                 step_through: bool = False,
                 n_workers: int = 1,
                 parent=None):
        super().__init__(module_keys, scint_keys, cuts,
                         wfm_collector, wfm_save_path, parent)
        self._paths        = list(evio_paths)
        self._step_through = step_through
        self._n_workers    = max(1, int(n_workers))
        self._continue_evt = threading.Event()
        self._tuple_to_role = _build_tuple_role_map(scint_keys, module_keys)

    def stop(self):
        self._stop_evt.set()
        self._continue_evt.set()   # unblock any pending step-through wait

    def step_continue(self):
        self._continue_evt.set()

    @staticmethod
    def _fetch_wfm_channels(fadc_evt, analyzer, needed_keys: frozenset) -> Dict[str, dict]:
        """Build a minimal channel dict (with raw samples) for a small set of keys.
        Called only for coincidence events when waveform saving or
        step-through is enabled."""
        channels: Dict[str, dict] = {}
        for r in range(fadc_evt.nrocs):
            roc     = fadc_evt.roc(r)
            roc_tag = int(roc.tag)
            for s in roc.present_slots():
                slot = roc.slot(s)
                for c in slot.present_channels():
                    key = f"{roc_tag}_{s}_{c}"
                    if key not in needed_keys:
                        continue
                    samples = slot.channel(c).samples
                    if samples.size < 4:
                        continue
                    ped_mean, ped_rms, peaks = analyzer.analyze(samples)
                    channels[key] = {
                        "pm": float(ped_mean),
                        "pr": float(ped_rms),
                        "s":  list(samples),
                        "pk": [{"i": float(p.integral), "h": float(p.height),
                                 "t": float(p.time),    "p": int(p.pos),
                                 "l": int(p.left),      "r": int(p.right),
                                 "o": int(p.overflow)}
                               for p in peaks],
                    }
        return channels

    def run(self):
        if not _HAVE_PRAD2PY:
            self.finished.emit("error: prad2py not available")
            return

        cuts     = self._cuts
        counters = Counters(self._mod_keys)

        # Parallel path: the files are split across worker subprocesses.
        # Skipped when waveform collection or step-through is enabled
        # (those need per-event UI).
        if (len(self._paths) > 1
                and not self._wfm_coll
                and not self._step_through):
            n_files = len(self._paths)
            # Partition files into n_workers near-equal chunks.  With 10 files
            # and 3 workers → sizes [4, 3, 3]; the first ``extra`` chunks each
            # take one extra file.
            n_workers = max(1, min(self._n_workers, n_files))
            base, extra = divmod(n_files, n_workers)
            chunks: List[List[str]] = []
            start = 0
            for i in range(n_workers):
                size = base + (1 if i < extra else 0)
                if size > 0:
                    chunks.append(list(self._paths[start:start + size]))
                start += size

            # "fork" inherits Qt state and is fragile from a Qt thread.
            # "forkserver" forks from a clean intermediate process, which
            # is much more reliable for GUI apps.
            try:
                ctx = multiprocessing.get_context("forkserver")
            except ValueError:
                ctx = multiprocessing.get_context("fork")
            update_q = ctx.Queue()
            pool     = ctx.Pool(processes=len(chunks),
                                initializer=_init_worker,
                                initargs=(update_q,))

            # Stash worker errors so they can be surfaced as a finished("error: …")
            # message — silent crashes are how a "processing ends immediately"
            # bug usually looks.
            worker_errors: List[str] = []
            def _on_err(exc):
                import traceback
                msg = "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
                sys.stderr.write(f"[coinc worker error]\n{msg}\n")
                sys.stderr.flush()
                worker_errors.append(repr(exc))

            file_args = [
                {
                    "paths":         chunk,
                    "tuple_to_role": self._tuple_to_role,
                    "cuts":          cuts,
                    "mod_names":     list(self._mod_keys),
                    "push_interval": 0.5,
                    "worker_id":     i,
                }
                for i, chunk in enumerate(chunks)
            ]
            self.workers_setup.emit(len(chunks))
            async_results = [pool.apply_async(_process_files, (a,),
                                              error_callback=_on_err)
                             for a in file_args]
            pool.close()

            # Coalesce per-worker progress so we emit at most one signal per
            # worker per UI cycle even when the worker pushes more often.
            latest_progress: Dict[int, dict] = {}

            last_emit = time.monotonic()
            try:
                while not self._stop_evt.is_set():
                    drained = False
                    while True:
                        try:
                            msg = update_q.get(timeout=0.1)
                        except Exception:
                            break
                        counters.merge(msg)
                        if "worker_id" in msg:
                            latest_progress[msg["worker_id"]] = msg
                        drained = True

                    n_chunks_total = len(async_results)
                    n_done = sum(1 for ar in async_results if ar.ready())

                    now = time.monotonic()
                    if drained and now - last_emit > 0.5:
                        self._emit(counters, n_done, n_chunks_total)
                        for prog in latest_progress.values():
                            self.worker_progress.emit(prog)
                        latest_progress.clear()
                        last_emit = now

                    if n_done >= n_chunks_total:
                        break
            finally:
                if self._stop_evt.is_set():
                    pool.terminate()

                # Wait for worker processes to fully exit.  Each worker's queue
                # feeder thread flushes pending messages on process exit, so
                # this is when the *final* per-worker _flush() lands in the
                # OS pipe buffer.  Without this wait, the parent would see
                # ar.ready()==True and proceed before the worker's last delta
                # message arrived, silently losing the tail of each chunk.
                pool.join()

                # Drain everything that's now sitting in the queue, including
                # the final flushes from every worker.
                while True:
                    try:
                        msg = update_q.get_nowait()
                    except Exception:
                        break
                    counters.merge(msg)
                    if "worker_id" in msg:
                        latest_progress[msg["worker_id"]] = msg
                for prog in latest_progress.values():
                    self.worker_progress.emit(prog)
                latest_progress.clear()

                # Surface any exceptions raised inside the workers.
                for ar in async_results:
                    if ar.ready():
                        try:
                            ar.get(timeout=0)
                        except Exception as e:
                            if not worker_errors:
                                _on_err(e)

            n_chunks_total = len(async_results)
            err = ""
            if worker_errors and not self._stop_evt.is_set():
                err = f"error: worker failed — {worker_errors[0]}"
            self._finish(counters, n_chunks_total, n_chunks_total, err)
            return

        # Sequential path: single file, or wfm/step-through requested.
        dec      = _prad2py.dec
        cfg      = dec.load_daq_config()
        ch       = dec.EvChannel()
        ch.set_config(cfg)
        analyzer = dec.WaveAnalyzer(cfg.wave_cfg)

        # Pre-scan to get total EVIO record count for the progress bar.
        total_records = 0
        for path in self._paths:
            if ch.open_auto(path) == dec.Status.success:
                if ch.is_random_access():
                    total_records += ch.get_random_access_event_count()
                else:
                    total_records = 0
                    ch.close()
                    break
                ch.close()
            else:
                ch.close()

        scint_thr     = cuts.scint_thr
        hycal_thr     = cuts.hycal_thr
        scint_t_min   = cuts.scint_t_min
        scint_t_max   = cuts.scint_t_max
        hycal_t_min   = cuts.hycal_t_min
        hycal_t_max   = cuts.hycal_t_max
        require_scint = cuts.require_scint
        wfm_coll      = self._wfm_coll
        tuple_to_role = self._tuple_to_role
        scint_key_set = frozenset(self._scint_keys.values())

        last_emit = time.monotonic()
        n_records = 0

        for path in self._paths:
            if self._stop_evt.is_set():
                break

            if ch.open_auto(path) != dec.Status.success:
                continue

            while not self._stop_evt.is_set():
                if ch.read() != dec.Status.success:
                    break
                n_records += 1
                if not ch.scan():
                    continue
                if ch.get_event_type() != dec.EventType.Physics:
                    continue

                for si in range(ch.get_n_events()):
                    ch.select_event(si)
                    info = ch.info()
                    counters.processed += 1

                    if int(info.trigger_bits) & SKIP_TRIGGER_MASK:
                        continue

                    fadc_evt = ch.fadc()
                    sf: Dict[str, bool]  = {}
                    ma: Dict[str, float] = {}

                    for r in range(fadc_evt.nrocs):
                        roc     = fadc_evt.roc(r)
                        roc_tag = int(roc.tag)
                        for s in roc.present_slots():
                            slot = roc.slot(s)
                            for c in slot.present_channels():
                                role = tuple_to_role.get((roc_tag, s, c))
                                if role is None:
                                    continue
                                samples = slot.channel(c).samples
                                if samples.size < 4:
                                    continue
                                # Cheap pre-filter — peak height above
                                # pedestal is bounded by max - min, so if
                                # that is below threshold we can skip the
                                # full waveform analysis entirely.
                                kind, name = role
                                thr_quick = scint_thr if kind == 0 else hycal_thr
                                if samples.max() - samples.min() < thr_quick:
                                    continue
                                _, _, peaks = analyzer.analyze(samples)
                                if kind == 0:
                                    sf[name] = any(
                                        p.height > scint_thr
                                        and scint_t_min <= p.time <= scint_t_max
                                        for p in peaks
                                    )
                                else:
                                    height = max(
                                        (p.height for p in peaks
                                         if hycal_t_min <= p.time <= hycal_t_max),
                                        default=0.0)
                                    if height > hycal_thr:
                                        ma[name] = height

                    best = _select_best(ma, cuts)
                    if best is None:
                        continue
                    counters.add(sf, best, require_scint)

                    # Waveforms are only shown / saved for events in which a
                    # scintillator fired.
                    want_wfm = wfm_coll is not None and not wfm_coll.full
                    if not (any(sf.values())
                            and (want_wfm or self._step_through)):
                        continue
                    best_key = self._mod_keys[best]
                    wfm_ch   = self._fetch_wfm_channels(
                        fadc_evt, analyzer, scint_key_set | {best_key})

                    if want_wfm:
                        data_stub = {
                            "trigger_bits": int(info.trigger_bits),
                            "event_number": int(info.event_number),
                            "event":        int(info.event_number),
                            "channels":     wfm_ch,
                        }
                        _collect_wfm(wfm_coll, data_stub, sf,
                                     self._scint_keys, best, best_key)

                    if self._step_through:
                        self._emit(counters, n_records, total_records)
                        self._continue_evt.clear()
                        self.coincidence_event.emit({
                            "event_number": int(info.event_number),
                            "best":         best,
                            "best_adc":     ma[best],
                            "scint_fired":  dict(sf),
                            "wfm_channels": wfm_ch,
                        })
                        while not self._stop_evt.is_set():
                            if self._continue_evt.wait(timeout=0.05):
                                break
                        if self._stop_evt.is_set():
                            break

                now = time.monotonic()
                if now - last_emit > 0.5:
                    self._emit(counters, n_records, total_records)
                    last_emit = now

            ch.close()

        self._finish(counters, total_records or n_records, total_records)


class InstantDisplayWorker(QThread):
    """Live event-by-event display from the ET ring buffer.

    Polls /api/ring at ~10 Hz, fetches each new latest event, computes
    per-module ADC values (max peak height above pedestal), and emits
    event_ready with the raw channels dict so the main thread can render
    both the map and the waveform panels.
    """
    event_ready = pyqtSignal(dict)   # {adc_vals, channels, seq, max_mod}
    finished    = pyqtSignal(str)

    POLL_INTERVAL = 0.1   # 10 Hz

    def __init__(self, server_url: str,
                 module_keys: Dict[str, str],
                 hycal_thr: float = 0.0,
                 min_mods: int = DEFAULT_MIN_CLUSTER_MODS,
                 parent=None):
        super().__init__(parent)
        self._url       = server_url
        self._mod_keys  = module_keys
        self._hycal_thr = hycal_thr
        self._min_mods  = min_mods
        self._stop_evt  = threading.Event()

    def stop(self):
        self._stop_evt.set()

    def run(self):
        last_seq = 0

        while not self._stop_evt.is_set():
            ring_data = _http_get(f"{self._url}/api/ring", timeout=2.0)
            if ring_data is None:
                time.sleep(self.POLL_INTERVAL)
                continue

            latest = ring_data.get("latest", 0)
            if latest == 0 or latest <= last_seq:
                time.sleep(self.POLL_INTERVAL)
                continue

            last_seq = latest
            data = _http_get(f"{self._url}/api/event/{latest}")
            if not data or "error" in data:
                time.sleep(self.POLL_INTERVAL)
                continue

            # Reject LMS / alpha events by trigger bit.
            if data.get("trigger_bits", 0) & SKIP_TRIGGER_MASK:
                time.sleep(self.POLL_INTERVAL)
                continue

            channels = data.get("channels", {})

            adc_vals: Dict[str, float] = {}
            for mname, mkey in self._mod_keys.items():
                ch = channels.get(mkey, {})
                peaks = ch.get("pk", [])
                val = max((float(pk.get("h", 0.0)) for pk in peaks), default=0.0)
                adc_vals[mname] = val if val > self._hycal_thr else 0.0

            # Skip events with too few modules (noise) or too many (LMS occupancy).
            n_above = sum(1 for v in adc_vals.values() if v > 0.0)
            if n_above < max(self._min_mods, 1) or n_above > LMS_MAX_MODULES:
                time.sleep(self.POLL_INTERVAL)
                continue

            max_mod = max(adc_vals, key=adc_vals.get) if adc_vals else ""

            self.event_ready.emit({
                "adc_vals":  adc_vals,
                "channels":  channels,
                "seq":       latest,
                "max_mod":   max_mod,
            })

            time.sleep(self.POLL_INTERVAL)

        self.finished.emit("" if not self._stop_evt.is_set() else "stopped")


# Waveform fetcher

class WaveformFetcher(QThread):
    """Fetches scintillator and HyCal module waveforms for one event.

    File mode: calls /api/waveform/{ev}/{key} for each channel.
    ET mode:   calls /api/event/latest and extracts channel data from
               the full event JSON (ring events include samples).
    """
    waveform_ready = pyqtSignal(str, dict)   # label ("scint"|"module"), data

    def __init__(self, url: str, event_n: int,
                 scint_keys: Dict[str, str], module_key: str,
                 et_mode: bool = False,
                 parent=None):
        super().__init__(parent)
        self._url         = url
        self._ev          = event_n
        self._scint_keys  = scint_keys   # name → channel key, all 4 scints
        self._module_key  = module_key
        self._et_mode     = et_mode

    def run(self):
        if self._et_mode:
            event_data = _http_get(f"{self._url}/api/event/latest")
            if not event_data or "error" in event_data:
                err = {"error": "no ET event available"}
                self.waveform_ready.emit("scint",  {n: err for n in self._scint_keys})
                self.waveform_ready.emit("module", err)
                return
            channels = event_data.get("channels", {})
            scint_wfms = {}
            for name, key in self._scint_keys.items():
                if not key:
                    scint_wfms[name] = {"error": "no channel"}
                else:
                    ch = channels.get(key)
                    scint_wfms[name] = ch if ch else {"error": f"channel {key} not in event"}
            self.waveform_ready.emit("scint", scint_wfms)
            if not self._module_key:
                self.waveform_ready.emit("module", {"error": "no channel"})
            else:
                ch = channels.get(self._module_key)
                self.waveform_ready.emit(
                    "module",
                    ch if ch else {"error": f"channel {self._module_key} not in event"})
        else:
            scint_wfms = {}
            for name, key in self._scint_keys.items():
                if not key:
                    scint_wfms[name] = {"error": "no channel"}
                else:
                    data = _http_get(f"{self._url}/api/waveform/{self._ev}/{key}")
                    scint_wfms[name] = data or {"error": "no response"}
            self.waveform_ready.emit("scint", scint_wfms)
            if not self._module_key:
                self.waveform_ready.emit("module", {"error": "no channel"})
            else:
                data = _http_get(f"{self._url}/api/waveform/{self._ev}/{self._module_key}")
                self.waveform_ready.emit("module", data or {"error": "no response"})


# Waveform plot widgets

class MultiWavePanel(QWidget):
    """FADC waveform display: multiple channels overlaid with a colour legend.

    Primary use: show all scintillators (V1-V4) on one shared axis so the
    operator can compare timing and amplitude across all channels at once.

    Each channel is the server's waveform JSON (/api/waveform/<n>/<key>):
        {"s": [int, ...], "pm": float, "pr": float,
         "pk": [{"p": int, "h": float, "i": float,
                 "l": int, "r": int, "t": float, "o": int}, ...]}
    """

    PAD_L, PAD_R, PAD_T, PAD_B = 52, 14, 28, 32

    _CHAN_COLORS = ["#4a9eff", "#ff6b6b", "#51cf66", "#ffd43b",
                    "#cc5de8", "#ff922b", "#20c997", "#f06595"]

    # Presentation switches (WavePanel sets the single-channel look).
    PER_PEAK_COLORS   = False   # colour peaks by index instead of by channel
    SHOW_PEDESTAL     = False   # dashed line at channel 0's pedestal
    SHOW_INFO         = False   # ped / rms / peak-count box for channel 0
    SHOW_LEGEND       = True
    TITLE_FIRED       = False   # [FIRED] tag of channel 0 after the title
    REJECT_FILL_ALPHA = 0.3     # fill alpha of peaks failing the cuts

    def __init__(self, label: str = "", parent=None):
        super().__init__(parent)
        self._label       = label
        self._title       = label
        self._channels: List[dict] = []
        self._threshold   = 0.0
        self._t_min       = -math.inf
        self._t_max       = math.inf
        self._placeholder = "No data — select an event and click Fetch"
        self.setMinimumHeight(120)
        self.setSizePolicy(QSizePolicy.Policy.Expanding,
                           QSizePolicy.Policy.Expanding)

    # Public API

    def set_multi_data(self, channels: dict, threshold: float,
                       scint_order: list,
                       title: Optional[str] = None,
                       t_min: float = -math.inf,
                       t_max: float = math.inf) -> None:
        """Overlay multiple waveforms. *channels* maps name → wave_json dict."""
        self._title     = title or self._label
        self._threshold = threshold
        self._t_min     = t_min
        self._t_max     = t_max
        self._channels  = []
        for i, name in enumerate(scint_order):
            wj = channels.get(name, {})
            if "error" in wj or "s" not in wj:
                wj = {}
            peaks = list(wj.get("pk", []))
            self._channels.append({
                "name": name, "color": self._chan_color(i),
                "samples":  list(wj.get("s", [])), "peaks": peaks,
                "ped_mean": float(wj.get("pm", 0)),
                "ped_rms":  float(wj.get("pr", 0)),
                "fired":    any(self._passes(pk) for pk in peaks),
            })
        self.update()

    def clear(self, title: Optional[str] = None,
              placeholder: Optional[str] = None) -> None:
        self._title    = title or self._label
        self._channels = []
        if placeholder is not None:
            self._placeholder = placeholder
        self.update()

    def _chan_color(self, i: int) -> str:
        return self._CHAN_COLORS[i % len(self._CHAN_COLORS)]

    def _peak_color(self, ch: dict, i: int) -> QColor:
        return series_qcolor(i) if self.PER_PEAK_COLORS else QColor(ch["color"])

    def _passes(self, pk: dict) -> bool:
        """Peak above threshold and inside the time window."""
        return (pk.get("h", 0) > self._threshold
                and self._t_min <= pk.get("t", 0.0) <= self._t_max)

    # Painting

    def paintEvent(self, _ev):
        p = QPainter(self)
        p.setRenderHint(QPainter.RenderHint.Antialiasing, True)
        p.fillRect(self.rect(), QColor(THEME.BG))

        r = QRectF(self.PAD_L, self.PAD_T,
                   max(1.0, self.width() - self.PAD_L - self.PAD_R),
                   max(1.0, self.height() - self.PAD_T - self.PAD_B))
        p.setPen(QColor(THEME.BORDER))
        p.setBrush(Qt.BrushStyle.NoBrush)
        p.drawRect(r)

        self._draw_title(p, r)

        drawn = [ch for ch in self._channels if ch["samples"]]
        n = max((len(ch["samples"]) for ch in drawn), default=0)
        if n < 2:
            p.setPen(QColor(THEME.TEXT_DIM))
            p.setFont(QFont("Monospace", 10))
            p.drawText(r, Qt.AlignmentFlag.AlignCenter, self._placeholder)
            if self.SHOW_LEGEND:
                self._draw_legend(p, r)
            return

        ymin, ymax = self._y_range([ch["samples"] for ch in drawn])

        def sx(i: float) -> float:
            return r.left() + i / max(1, n - 1) * r.width()

        def sy(v: float) -> float:
            return r.bottom() - (v - ymin) / max(1e-6, ymax - ymin) * r.height()

        if self.SHOW_PEDESTAL:
            self._draw_pedestal(p, r, sy)
        for ch in drawn:
            self._draw_peak_fills(p, sx, sy, ch)
        for ch in drawn:
            self._draw_waveform(p, sx, sy, ch)
        for ch in drawn:
            self._draw_peak_markers(p, sx, sy, ch)
        draw_wave_axes(p, r, ymin, ymax, n, CLK_MHZ, self.PAD_L)
        if self.SHOW_LEGEND:
            self._draw_legend(p, r)
        if self.SHOW_INFO:
            self._draw_info(p, r)

    def _y_range(self, all_samples):
        flat = [v for s in all_samples for v in s]
        ymin, ymax = float(min(flat)), float(max(flat))
        if ymax - ymin < 5.0:
            ymax = ymin + 5.0
        pad = (ymax - ymin) * 0.06
        return ymin - pad, ymax + pad

    def _draw_title(self, p: QPainter, r: QRectF):
        f = QFont("Monospace", 10)
        f.setBold(True)
        p.setFont(f)
        p.setPen(QColor(THEME.TEXT))
        p.drawText(int(r.left()), int(r.top() - 8), self._title)
        if self.TITLE_FIRED and self._channels and self._channels[0]["samples"]:
            fired = self._channels[0]["fired"]
            tw = p.fontMetrics().horizontalAdvance(self._title)
            fired_txt = "  [FIRED]" if fired else "  [—]"
            p.setPen(QColor(THEME.SUCCESS) if fired else QColor(THEME.TEXT_DIM))
            p.drawText(int(r.left() + tw), int(r.top() - 8), fired_txt)

    def _draw_pedestal(self, p: QPainter, r: QRectF, sy):
        ped = self._channels[0]["ped_mean"]
        if not ped:
            return
        y_ped = sy(ped)
        p.setPen(QPen(QColor(THEME.TEXT_DIM), 1, Qt.PenStyle.DashLine))
        p.drawLine(int(r.left()), int(y_ped), int(r.right()), int(y_ped))

    def _draw_peak_fills(self, p, sx, sy, ch):
        ped = ch["ped_mean"]
        if not ped:
            return
        y_ped    = sy(ped)
        samples  = ch["samples"]
        cn       = len(samples)
        for i, pk in enumerate(ch["peaks"]):
            base = self._peak_color(ch, i)
            # dim peaks that fail height threshold or fall outside the time window
            if not self._passes(pk):
                base.setAlphaF(self.REJECT_FILL_ALPHA)
            fill = QColor(base)
            fill.setAlphaF(fill.alphaF() * 0.25)
            lft  = max(0, int(pk.get("l", pk["p"])))
            rgt  = min(cn - 1, int(pk.get("r", pk["p"])))
            poly = QPolygonF()
            for k in range(lft, rgt + 1):
                poly.append(QPointF(sx(k), sy(samples[k])))
            poly.append(QPointF(sx(rgt), y_ped))
            poly.append(QPointF(sx(lft), y_ped))
            p.setPen(Qt.PenStyle.NoPen)
            p.setBrush(fill)
            p.drawPolygon(poly)

    def _draw_waveform(self, p, sx, sy, ch):
        samples = ch["samples"]
        cn      = len(samples)
        p.setPen(QPen(QColor(ch["color"]), 1.4))
        p.setBrush(Qt.BrushStyle.NoBrush)
        for i in range(cn - 1):
            p.drawLine(int(sx(i)),     int(sy(samples[i])),
                       int(sx(i + 1)), int(sy(samples[i + 1])))

    def _draw_peak_markers(self, p, sx, sy, ch):
        samples = ch["samples"]
        cn      = len(samples)
        for i, pk in enumerate(ch["peaks"]):
            pos = int(pk.get("p", 0))
            if pos < 0 or pos >= cn:
                continue
            col = self._peak_color(ch, i)
            if not self._passes(pk):
                col.setAlphaF(0.4)
            p.setPen(QPen(col, 1.2))
            p.setBrush(col)
            cx, cy = sx(pos), sy(float(samples[pos]))
            diamond = QPolygonF([
                QPointF(cx,     cy - 4),
                QPointF(cx + 4, cy),
                QPointF(cx,     cy + 4),
                QPointF(cx - 4, cy),
            ])
            p.drawPolygon(diamond)
            # height label above the diamond
            p.setFont(QFont("Monospace", 8))
            p.setPen(col)
            p.drawText(int(cx - 16), int(cy - 7), f"{pk.get('h', 0):.0f}")

    def _draw_legend(self, p: QPainter, r: QRectF):
        if not self._channels:
            return
        f = QFont("Monospace", 9)
        p.setFont(f)
        fm = p.fontMetrics()
        lh = fm.height() + 3

        entries = []
        for ch in self._channels:
            peaks_in  = [pk for pk in ch["peaks"] if self._passes(pk)]
            t_str     = f" t={peaks_in[0].get('t', 0):.0f}ns" if peaks_in else ""
            fired_txt = "[FIRED]" if ch["fired"] else "[—]"
            label     = f"{ch['name']} {fired_txt}{t_str}"
            entries.append((label, ch["color"], ch["fired"]))

        max_w   = max(fm.horizontalAdvance(e[0]) for e in entries) + 20
        total_h = lh * len(entries) + 8

        lx = int(r.right() - max_w - 6)
        ly = int(r.top() + 4)

        bg = QColor(THEME.BG)
        bg.setAlphaF(0.82)
        p.fillRect(QRectF(lx - 2, ly, max_w + 4, total_h), bg)

        y = ly + 4
        for label, color, fired in entries:
            col = QColor(color)
            p.fillRect(QRectF(lx, y + 1, 12, lh - 2), col)
            p.setPen(QColor(THEME.SUCCESS) if fired else QColor(THEME.TEXT_DIM))
            p.setFont(f)
            p.drawText(int(lx + 16), int(y + fm.ascent()), label)
            y += lh

    def _draw_info(self, p, r):
        ch = self._channels[0]
        above = sum(1 for pk in ch["peaks"] if self._passes(pk))
        info = (f"ped={ch['ped_mean']:.1f}  rms={ch['ped_rms']:.2f}"
                f"  peaks={len(ch['peaks'])} ({above} above thr)")
        p.setFont(QFont("Monospace", 9))
        fm = p.fontMetrics()
        tw = fm.horizontalAdvance(info)
        th = fm.height()
        box = QRectF(r.right() - tw - 8, r.top() + 4, tw + 6, th + 2)
        bg = QColor(THEME.BG)
        bg.setAlphaF(0.75)
        p.fillRect(box, bg)
        p.setPen(QColor(THEME.TEXT_DIM))
        p.drawText(box, Qt.AlignmentFlag.AlignCenter, info)


class WavePanel(MultiWavePanel):
    """Single-channel waveform display: peaks colour-coded by index, the
    pedestal line, a [FIRED] tag and a ped / rms / peaks info box."""

    PER_PEAK_COLORS   = True
    SHOW_PEDESTAL     = True
    SHOW_INFO         = True
    SHOW_LEGEND       = False
    TITLE_FIRED       = True
    REJECT_FILL_ALPHA = 0.4

    def _chan_color(self, i: int) -> str:
        return THEME.ACCENT

    def set_data(self, wave_json: dict, threshold: float,
                 title: Optional[str] = None,
                 t_min: float = -math.inf,
                 t_max: float = math.inf) -> None:
        """Load waveform from the server JSON response; an error response
        becomes the placeholder text."""
        if "error" in wave_json:
            self.clear(title or self._label,
                       wave_json.get("error", "channel not found"))
            return
        self.set_multi_data({self._label: wave_json}, threshold,
                            [self._label], title, t_min, t_max)


# Coincidence map widget

class CoincidenceMapWidget(HyCalMapWidget):
    """HyCal map coloured by coincidence rate, with informative tooltip."""

    def __init__(self, parent=None):
        super().__init__(parent, enable_zoom_pan=True, show_colorbar=True)
        self._snapshot: dict = {}
        self._active_scint = "V1"
        self._view_mode = VIEW_COINC
        self._instant_max_mod = ""
        self.set_palette("rainbow")

    def set_snapshot(self, snap: dict, scint: str,
                     view_mode: str = VIEW_COINC) -> None:
        self._snapshot    = snap
        self._active_scint = scint
        self._view_mode   = view_mode

        if view_mode == VIEW_OCC:
            hits = snap.get("module_hits", {})
            self.set_values(hits if hits else {})
            max_hits = max(hits.values()) if hits else 1
            self.set_range(0, max_hits)
        else:
            rates = snap.get("rates", {}).get(scint, {})
            valid = {m: v for m, v in rates.items() if not math.isnan(v)}
            self.set_values(valid if valid else {})
            self.set_range(0.0, 1.0)
        self.update()

    def set_instant_event(self, adc_vals: dict, max_mod: str = "") -> None:
        """Display per-module ADC values for a single live event."""
        self._view_mode = VIEW_INSTANT
        self._instant_max_mod = max_mod
        self.set_values(adc_vals)
        max_val = max(adc_vals.values(), default=0.0)
        self.set_range(0.0, max(max_val, 1.0))
        self.update()

    # Dark purple used for modules that have been seen (module_hits > 0) but
    # have zero coincidence rate or zero occupancy count.  Distinct from the
    # NO_DATA_COLOR (dark grey) which means the module was never seen at all.
    _SEEN_ZERO = QColor(55, 0, 75)

    def _paint_modules(self, p: QPainter):
        if self._view_mode not in (VIEW_COINC, VIEW_OCC):
            super()._paint_modules(p)
            return
        stops = self.palette_stops()
        no_data = self.NO_DATA_COLOR
        seen_zero = self._SEEN_ZERO
        for name, rect in self._rects.items():
            v = self._values.get(name)
            if v is None or (isinstance(v, float) and math.isnan(v)):
                p.fillRect(rect, no_data)
            elif v == 0.0:
                p.fillRect(rect, seen_zero)
            else:
                p.fillRect(rect, cmap_qcolor(self.value_to_t(v), stops))

    def _tooltip_text(self, name: str) -> str:
        if self._view_mode == VIEW_INSTANT:
            adc = self._values.get(name, 0.0)
            suffix = "  ★ max" if name == self._instant_max_mod else ""
            return f"{name}\nADC: {adc:.0f}{suffix}"

        snap  = self._snapshot
        scint = self._active_scint
        nhits = snap.get("module_hits", {}).get(name, 0)

        if self._view_mode == VIEW_OCC:
            processed = snap.get("processed", 0)
            frac = nhits / processed if processed > 0 else float("nan")
            if processed == 0:
                return f"{name}\nNo data"
            return (f"{name}\n"
                    f"Module hits: {nhits:,}\n"
                    f"Occupancy: {frac:.4f}")

        rate = snap.get("rates", {}).get(scint, {}).get(name, math.nan)
        if math.isnan(rate):
            return f"{name}\nNo data"
        mode = snap.get("mode", MODE_AND)
        if mode == MODE_AND:
            denom = nhits
            denom_label = "Module hits"
        else:
            denom = snap.get("scint_hits", {}).get(scint, 0)
            denom_label = f"{scint} hits"
        ncoinc = round(rate * denom)
        return (f"{name}\n"
                f"Coincidence rate: {rate:.4f}\n"
                f"Coincidences: {ncoinc:,}\n"
                f"{denom_label}: {denom:,}")


# Main window

class MainWindow(QMainWindow):

    def __init__(self, server_url: str):
        super().__init__()
        self.setWindowTitle("Scintillator–HyCal Coincidence Monitor")
        self.resize(1400, 900)

        self._server_url      = server_url
        self._n_events        = 0
        self._local_evio_paths: List[str] = []   # files queued for local analysis
        self._stats_worker:   Optional[_CoincWorker]         = None
        self._instant_worker: Optional[InstantDisplayWorker] = None
        self._fetcher:  Optional[WaveformFetcher] = None
        self._snapshot: dict = {}
        self._selected_module: str = ""   # module last clicked on the map

        # Channel keys "<roc>_<slot>_<channel>": the event JSON of the C++
        # server uses the actual ROC tag, not the crate index, for <roc>.
        try:
            crate_to_roc = {crate: tag for tag, crate
                            in load_roc_tag_map(DAQ_CFG_JSON).items()}
        except Exception:
            crate_to_roc = {}
        daq = {name: f"{crate_to_roc.get(crate, crate)}_{slot}_{chan}"
               for (crate, slot, chan), name
               in load_daq_map(MODULES_JSON).items()}
        self._modules = load_modules(MODULES_JSON)
        # W module → ring around the beam hole (1 = innermost)
        self._w_layers: Dict[str, int] = {
            m.name: hole_ring(m.row, m.col) for m in self._modules
            if m.name.startswith("W") and m.row
        }
        self._neighbor_map = _build_neighbor_map(self._modules)
        physics_names = {m.name for m in self._modules if m.mod_type != "LMS"}
        self._mod_keys: Dict[str, str] = {
            name: key for name, key in daq.items()
            if name in physics_names and name not in SCINTILLATORS
        }
        # Scintillator channel keys resolved from daq_map (V1-V4 are listed there)
        self._scint_keys: Dict[str, str] = {
            sname: daq[sname] for sname in SCINTILLATORS if sname in daq
        }

        self._veto_ctrl    = VetoMotorController()
        self._et_mode      = False         # True when server is in online/ET mode
        self._map_view     = VIEW_COINC
        self._display_mode = DISPLAY_COINC # top-level mode switch
        self._display_paused = False       # True while instant display is frozen

        self._build_ui()
        self._map.set_modules(self._modules)
        self._map.set_palette("rainbow")
        self._map.set_range(0.0, 1.0)

        self._veto_timer = QTimer(self)
        self._veto_timer.timeout.connect(self._poll_veto_positions)
        self._veto_timer.start(VETO_POLL_MS)

    # UI construction

    def _build_ui(self):
        h_split = QSplitter(Qt.Orientation.Horizontal)

        # ---- left control panel ----------------------------------------
        left = QWidget()
        left.setStyleSheet(themed(f"QWidget{{background:{THEME.PANEL};}}"))
        lv = QVBoxLayout(left)
        lv.setContentsMargins(10, 10, 10, 10)
        lv.setSpacing(8)

        # Server
        srv_box = QGroupBox("Server")
        srv_box.setStyleSheet(self._groupbox_style())
        sv = QVBoxLayout(srv_box)
        sv.setSpacing(4)
        self._url_edit = QLineEdit(self._server_url)
        self._url_edit.setStyleSheet(self._input_style())
        sv.addWidget(self._url_edit)
        self._connect_btn = QPushButton("Connect")
        self._connect_btn.setStyleSheet(self._btn_style())
        self._connect_btn.clicked.connect(self._on_connect)
        sv.addWidget(self._connect_btn)
        self._conn_label = self._dim_label("Not connected")
        sv.addWidget(self._conn_label)
        self._mode_btn = QPushButton("Switch to ET Mode")
        self._mode_btn.setEnabled(False)
        self._mode_btn.setStyleSheet(self._btn_style())
        self._mode_btn.clicked.connect(self._on_mode_switch)
        sv.addWidget(self._mode_btn)
        lv.addWidget(srv_box)

        # Local EVIO files (prad2py direct — no server needed)
        local_box = QGroupBox("Local EVIO Files (no server)")
        local_box.setStyleSheet(self._groupbox_style())
        loc = QVBoxLayout(local_box)
        loc.setSpacing(4)

        # Button row: Add / Remove selected / Clear all
        loc_btn_row = QHBoxLayout()
        loc_btn_row.setSpacing(4)
        self._browse_btn = QPushButton("Add files…")
        self._browse_btn.setStyleSheet(self._btn_style())
        if not _HAVE_PRAD2PY:
            self._browse_btn.setEnabled(False)
            self._browse_btn.setToolTip("prad2py not found — rebuild with -DBUILD_PYTHON=ON")
        self._browse_btn.clicked.connect(self._on_add_evio_files)
        loc_btn_row.addWidget(self._browse_btn)
        self._remove_file_btn = QPushButton("Remove")
        self._remove_file_btn.setStyleSheet(self._btn_style())
        self._remove_file_btn.setEnabled(False)
        self._remove_file_btn.clicked.connect(self._on_remove_evio_file)
        loc_btn_row.addWidget(self._remove_file_btn)
        self._clear_files_btn = QPushButton("Clear")
        self._clear_files_btn.setStyleSheet(self._btn_style())
        self._clear_files_btn.setEnabled(False)
        self._clear_files_btn.clicked.connect(self._on_clear_evio_files)
        loc_btn_row.addWidget(self._clear_files_btn)
        loc.addLayout(loc_btn_row)

        # List showing selected files (basename; full path in tooltip)
        self._file_list = QListWidget()
        self._file_list.setMaximumHeight(80)
        self._file_list.setSelectionMode(
            QAbstractItemView.SelectionMode.SingleSelection)
        self._file_list.setStyleSheet(
            themed(f"QListWidget{{background:{THEME.PANEL};"
                   f"color:{THEME.TEXT};border:1px solid {THEME.BORDER};"
                   f"font-size:10px;}}"
                   f"QListWidget::item:selected{{background:{THEME.ACCENT};}}"))
        self._file_list.itemSelectionChanged.connect(self._on_file_selection_changed)
        loc.addWidget(self._file_list)

        self._local_file_label = QLabel("No files selected")
        self._local_file_label.setStyleSheet(
            f"color:{THEME.TEXT_DIM};font-size:10px;")
        self._local_file_label.setWordWrap(True)
        loc.addWidget(self._local_file_label)
        lv.addWidget(local_box)

        # Display Mode
        disp_box = QGroupBox("Display Mode")
        disp_box.setStyleSheet(self._groupbox_style())
        dpv = QVBoxLayout(disp_box)
        dpv.setSpacing(2)
        self._disp_group = QButtonGroup(self)
        self._rb_disp_coinc   = QRadioButton("Coincidence Stats")
        self._rb_disp_instant = QRadioButton("Instant Event Display")
        for rb in (self._rb_disp_coinc, self._rb_disp_instant):
            rb.setStyleSheet(f"QRadioButton{{color:{THEME.TEXT};font-size:12px;}}")
            self._disp_group.addButton(rb)
            dpv.addWidget(rb)
        self._rb_disp_coinc.setChecked(True)
        self._disp_group.buttonClicked.connect(self._on_display_mode_changed)
        lv.addWidget(disp_box)

        # Scintillator selector
        sci_box = QGroupBox("Scintillator")
        sci_box.setStyleSheet(self._groupbox_style())
        scv = QVBoxLayout(sci_box)
        scv.setSpacing(2)
        self._scint_group = QButtonGroup(self)
        for name in SCINTILLATORS:
            rb = QRadioButton(name)
            rb.setStyleSheet(f"QRadioButton{{color:{THEME.TEXT};font-size:12px;}}")
            self._scint_group.addButton(rb)
            scv.addWidget(rb)
            if name == "V1":
                rb.setChecked(True)
        self._scint_group.buttonClicked.connect(self._on_scint_changed)
        lv.addWidget(sci_box)

        # Map view mode (Coincidence Stats only)
        self._mapview_box = QGroupBox("Map View")
        self._mapview_box.setStyleSheet(self._groupbox_style())
        mvv = QVBoxLayout(self._mapview_box)
        mvv.setSpacing(2)
        self._mapview_group = QButtonGroup(self)
        self._rb_coinc = QRadioButton("Coincidence Rate")
        self._rb_occ   = QRadioButton("Occupancy (hit count)")
        for rb in (self._rb_coinc, self._rb_occ):
            rb.setStyleSheet(
                f"QRadioButton{{color:{THEME.TEXT};font-size:12px;}}")
            self._mapview_group.addButton(rb)
            mvv.addWidget(rb)
        self._rb_coinc.setChecked(True)
        self._mapview_group.buttonClicked.connect(self._on_mapview_changed)
        lv.addWidget(self._mapview_box)

        # Thresholds
        thr_box = QGroupBox("Thresholds (ADC peak height)")
        thr_box.setStyleSheet(self._groupbox_style())
        tv = QVBoxLayout(thr_box)
        tv.setSpacing(4)
        tv.addWidget(self._dim_label("Scintillator:"))
        self._scint_thr_spin = self._spin(0, 100000, DEFAULT_SCINT_THR, 50,
                                          decimals=0)
        tv.addWidget(self._scint_thr_spin)

        tv.addWidget(self._dim_label("Scint time cut (ns):"))
        self._scint_tcut_min = self._spin(0, 10000, DEFAULT_SCINT_TMIN, 4,
                                          decimals=0, suffix=" ns")
        self._scint_tcut_max = self._spin(0, 10000, DEFAULT_SCINT_TMAX, 4,
                                          decimals=0, suffix=" ns")
        tv.addLayout(self._range_row(self._scint_tcut_min,
                                     self._scint_tcut_max))

        tv.addWidget(self._dim_label("HyCal module:"))
        self._hycal_thr_spin = self._spin(0, 100000, DEFAULT_HYCAL_THR, 10,
                                          decimals=0)
        tv.addWidget(self._hycal_thr_spin)

        tv.addWidget(self._dim_label("HyCal time cut (ns):"))
        self._hycal_tcut_min = self._spin(0, 10000, DEFAULT_HYCAL_TMIN, 4,
                                          decimals=0, suffix=" ns")
        self._hycal_tcut_max = self._spin(0, 10000, DEFAULT_HYCAL_TMAX, 4,
                                          decimals=0, suffix=" ns")
        tv.addLayout(self._range_row(self._hycal_tcut_min,
                                     self._hycal_tcut_max))

        tv.addWidget(self._dim_label(
            "Max HyCal local maxima (1 = single cluster):"))
        self._max_lm_spin = self._spin(1, 20, DEFAULT_MAX_LOCAL_MAXIMA)
        tv.addWidget(self._max_lm_spin)

        n_cpu_max = max(1, os.cpu_count() or 1)
        tv.addWidget(self._dim_label(
            f"Parallel CPUs for local files (1–{n_cpu_max}):"))
        self._cpu_spin = self._spin(1, n_cpu_max, min(4, n_cpu_max))
        tv.addWidget(self._cpu_spin)

        tv.addWidget(self._dim_label(
            "Inner W exclusion layers (0 = none, 16 = all W):"))
        self._excl_spin = self._spin(0, 16, 0)
        tv.addWidget(self._excl_spin)
        tv.addWidget(self._dim_label(
            "Min HyCal modules fired (cluster cut, 1 = off):"))
        self._min_mods_spin = self._spin(1, 50, DEFAULT_MIN_CLUSTER_MODS)
        tv.addWidget(self._min_mods_spin)
        lv.addWidget(thr_box)

        # Event selection mode (Coincidence Stats only)
        self._selmode_box = QGroupBox("Event Selection Mode")
        self._selmode_box.setStyleSheet(self._groupbox_style())
        mv = QVBoxLayout(self._selmode_box)
        mv.setSpacing(4)
        self._mode_grp = QButtonGroup(self)
        self._rb_and = QRadioButton("AND — selected veto + HyCal fire")
        self._rb_or  = QRadioButton("OR  — all HyCal events (no veto gate)")
        for rb in (self._rb_and, self._rb_or):
            rb.setStyleSheet(
                f"QRadioButton{{color:{THEME.TEXT};font-size:12px;}}")
            self._mode_grp.addButton(rb)
            mv.addWidget(rb)
        self._rb_and.setChecked(True)
        lv.addWidget(self._selmode_box)

        # Processing
        proc_box = QGroupBox("Processing")
        proc_box.setStyleSheet(self._groupbox_style())
        pv = QVBoxLayout(proc_box)
        pv.setSpacing(6)

        self._start_btn = QPushButton("Start")
        self._start_btn.setEnabled(False)
        self._start_btn.setStyleSheet(self._btn_style(accent=True))
        self._start_btn.clicked.connect(self._on_start_stop)
        pv.addWidget(self._start_btn)

        pause_row = QHBoxLayout()
        pause_row.setSpacing(4)
        self._pause_btn = QPushButton("Pause")
        self._pause_btn.setEnabled(False)
        self._pause_btn.setStyleSheet(self._btn_style())
        self._pause_btn.clicked.connect(self._on_pause)
        pause_row.addWidget(self._pause_btn)
        self._resume_btn = QPushButton("Resume")
        self._resume_btn.setEnabled(False)
        self._resume_btn.setStyleSheet(self._btn_style())
        self._resume_btn.clicked.connect(self._on_resume)
        pause_row.addWidget(self._resume_btn)
        pv.addLayout(pause_row)

        self._screenshot_btn = QPushButton("Save screenshot…")
        self._screenshot_btn.setStyleSheet(self._btn_style())
        self._screenshot_btn.setToolTip(
            "Save a PNG snapshot of the entire window.")
        self._screenshot_btn.clicked.connect(self._on_screenshot)
        pv.addWidget(self._screenshot_btn)

        self._progress = QProgressBar()
        self._progress.setRange(0, 100)
        self._progress.setValue(0)
        self._progress.setStyleSheet(
            f"QProgressBar{{background:{THEME.PANEL};border:1px solid "
            f"{THEME.BORDER};border-radius:4px;height:14px;}}"
            f"QProgressBar::chunk{{background:{THEME.ACCENT_STRONG};"
            f"border-radius:3px;}}")
        pv.addWidget(self._progress)

        # Per-worker progress rows (parallel local-EVIO mode).  Hidden until a
        # parallel job kicks off; populated with one row per worker process.
        self._worker_prog_container = QWidget()
        wpc_layout = QVBoxLayout(self._worker_prog_container)
        wpc_layout.setSpacing(2)
        wpc_layout.setContentsMargins(0, 0, 0, 0)
        self._worker_prog_layout: QVBoxLayout = wpc_layout
        self._worker_prog_rows: List[tuple] = []   # [(label, bar), ...]
        self._worker_prog_container.setVisible(False)
        pv.addWidget(self._worker_prog_container)

        rate_row = QHBoxLayout()
        rate_row.setSpacing(6)
        rate_row.addWidget(self._dim_label("Max ET rate:"))
        self._max_rate_spin = self._spin(0, 10000, 0, 50, suffix=" ev/s")
        self._max_rate_spin.setSpecialValueText("unlimited")
        self._max_rate_spin.setToolTip(
            "Maximum events per second consumed from the ET ring buffer.\n"
            "0 = process as fast as possible (unlimited).\n"
            "Only effective in ET mode.")
        rate_row.addWidget(self._max_rate_spin, 1)
        pv.addLayout(rate_row)

        # Waveform saving
        sep = QFrame()
        sep.setFrameShape(QFrame.Shape.HLine)
        sep.setStyleSheet(f"color:{THEME.BORDER};")
        pv.addWidget(sep)

        self._save_wfm_chk = QCheckBox("Save coincidence waveforms")
        self._save_wfm_chk.setStyleSheet(
            f"QCheckBox{{color:{THEME.TEXT};font-size:12px;}}")
        pv.addWidget(self._save_wfm_chk)

        wfm_row = QHBoxLayout()
        wfm_row.setSpacing(6)
        wfm_row.addWidget(self._dim_label("Max records:"))
        self._wfm_max_spin = self._spin(1, 100000, 200, 100)
        self._wfm_max_spin.setToolTip(
            "Maximum number of coincidence waveform pairs to save.")
        wfm_row.addWidget(self._wfm_max_spin, 1)
        pv.addLayout(wfm_row)

        self._wfm_label = QLabel("")
        self._wfm_label.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:10px;")
        self._wfm_label.setWordWrap(True)
        pv.addWidget(self._wfm_label)

        self._status_label = self._dim_label("Ready")
        self._status_label.setWordWrap(True)
        pv.addWidget(self._status_label)
        lv.addWidget(proc_box)

        # Event browser
        ev_box = QGroupBox("Event Browser")
        ev_box.setStyleSheet(self._groupbox_style())
        ev = QVBoxLayout(ev_box)
        ev.setSpacing(4)

        ev_row = QHBoxLayout()
        ev_row.setSpacing(4)
        ev_row.addWidget(self._dim_label("Event #"))
        self._ev_spin = self._spin(1, 1, 1)
        ev_row.addWidget(self._ev_spin, 1)
        ev.addLayout(ev_row)

        self._fetch_btn = QPushButton("Fetch waveforms")
        self._fetch_btn.setEnabled(False)
        self._fetch_btn.setStyleSheet(self._btn_style())
        self._fetch_btn.clicked.connect(self._on_fetch_waveforms)
        ev.addWidget(self._fetch_btn)

        self._sel_mod_label = self._dim_label("Module: (click map)")
        self._sel_mod_label.setWordWrap(True)
        ev.addWidget(self._sel_mod_label)
        lv.addWidget(ev_box)

        # Statistics
        stats_box = QGroupBox("Statistics")
        stats_box.setStyleSheet(self._groupbox_style())
        stv = QVBoxLayout(stats_box)
        self._stats_label = self._dim_label("—")
        self._stats_label.setWordWrap(True)
        stv.addWidget(self._stats_label)
        lv.addWidget(stats_box)

        # Veto motor positions — two sub-rows per motor:
        #   row A: name  |  "RBV:"  |  readback value  |  moving indicator
        #   row B: ""    |  "Set:"  |  spinbox          |  Move button
        veto_box = QGroupBox("Veto Motor Positions")
        veto_box.setStyleSheet(self._groupbox_style())
        vg = QGridLayout(veto_box)
        vg.setHorizontalSpacing(6)
        vg.setVerticalSpacing(2)
        vg.setContentsMargins(6, 8, 6, 8)
        self._veto_widgets: Dict[str, dict] = {}
        for i, vname in enumerate(VETO_PV_BASES):
            ra = i * 3       # row A: readback
            rb = i * 3 + 1  # row B: setpoint
            rc = i * 3 + 2  # row C: thin separator

            lbl_name = QLabel(vname)
            lbl_name.setStyleSheet(
                f"color:{THEME.ACCENT};font-size:12px;font-weight:bold;")
            vg.addWidget(lbl_name, ra, 0)

            vg.addWidget(self._dim_label("RBV:"), ra, 1)

            lbl_rbv = QLabel("—")
            lbl_rbv.setStyleSheet(
                f"color:{THEME.TEXT};font-size:12px;font-weight:bold;")
            lbl_rbv.setMinimumWidth(70)
            vg.addWidget(lbl_rbv, ra, 2)

            lbl_movn = QLabel("● moving")
            lbl_movn.setStyleSheet(
                f"color:{THEME.TEXT_DIM};font-size:10px;")
            vg.addWidget(lbl_movn, ra, 3)

            vg.addWidget(self._dim_label("Set:"), rb, 1)

            spin = self._spin(-9999, 9999, 0.0, 0.1, decimals=2)
            vg.addWidget(spin, rb, 2)

            btn = QPushButton("Move")
            btn.setStyleSheet(self._btn_style(accent=True))
            btn.clicked.connect(lambda _, v=vname: self._move_veto(v))
            vg.addWidget(btn, rb, 3)

            sep = QFrame()
            sep.setFrameShape(QFrame.Shape.HLine)
            sep.setStyleSheet(f"color:{THEME.BORDER};")
            sep.setFixedHeight(6)
            vg.addWidget(sep, rc, 0, 1, 4)

            self._veto_widgets[vname] = {
                "rbv": lbl_rbv, "movn": lbl_movn,
                "spin": spin, "btn": btn,
            }

        if not self._veto_ctrl.available:
            warn = QLabel("pyepics not available")
            warn.setStyleSheet(f"color:{THEME.WARN};font-size:10px;")
            vg.addWidget(warn, len(VETO_PV_BASES) * 3, 0, 1, 4)
        lv.addWidget(veto_box)

        lv.addStretch(1)

        # Wrap left panel in a scroll area so the window can resize freely
        # vertically even when the panel content is taller than the window.
        left_scroll = QScrollArea()
        left_scroll.setWidget(left)
        left_scroll.setWidgetResizable(True)
        left_scroll.setHorizontalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAlwaysOff)
        left_scroll.setVerticalScrollBarPolicy(
            Qt.ScrollBarPolicy.ScrollBarAsNeeded)
        left_scroll.setMinimumWidth(340)
        left_scroll.setMaximumWidth(500)
        left_scroll.setStyleSheet(
            themed(f"QScrollArea{{background:{THEME.PANEL};"
                   f"border:none;}}"))

        # ---- right side: map on top, waveforms on bottom ---------------
        v_split = QSplitter(Qt.Orientation.Vertical)

        # Thin toolbar above the HyCal map for step-through controls.
        map_container = QWidget()
        map_vl = QVBoxLayout(map_container)
        map_vl.setContentsMargins(0, 0, 0, 0)
        map_vl.setSpacing(2)

        step_bar = QWidget()
        step_bar.setStyleSheet(
            themed(f"QWidget{{background:{THEME.PANEL};"
                   f"border-bottom:1px solid {THEME.BORDER};}}"))
        step_hl = QHBoxLayout(step_bar)
        step_hl.setContentsMargins(8, 4, 8, 4)
        step_hl.setSpacing(8)

        self._step_chk = QCheckBox("Step through coincidence events")
        self._step_chk.setStyleSheet(
            f"QCheckBox{{color:{THEME.TEXT};font-size:12px;}}"
            f"QCheckBox::indicator{{width:14px;height:14px;"
            f"background:#ffffff;border:2px solid #888888;border-radius:3px;}}"
            f"QCheckBox::indicator:checked{{background:#4a9eff;"
            f"border:2px solid #4a9eff;}}"
        )
        self._step_chk.setToolTip(
            "Pause on each coincidence event found in the local EVIO file\n"
            "and show its waveforms in the panels below.\n"
            "Click 'Continue →' to advance to the next coincidence.\n"
            "Only available with local EVIO files.")
        step_hl.addWidget(self._step_chk)

        self._step_continue_btn = QPushButton("Continue →")
        self._step_continue_btn.setEnabled(False)
        self._step_continue_btn.setStyleSheet(self._btn_style(accent=True))
        self._step_continue_btn.setToolTip("Advance to the next coincidence event.")
        self._step_continue_btn.clicked.connect(self._on_step_continue)
        step_hl.addWidget(self._step_continue_btn)
        step_hl.addStretch(1)

        map_vl.addWidget(step_bar)

        self._map = CoincidenceMapWidget()
        self._map.moduleClicked.connect(self._on_module_clicked)
        map_vl.addWidget(self._map)

        v_split.addWidget(map_container)

        wave_row = QWidget()
        wave_row.setMinimumHeight(280)
        wrl = QHBoxLayout(wave_row)
        wrl.setContentsMargins(0, 0, 0, 0)
        wrl.setSpacing(4)
        self._wave_scint  = MultiWavePanel("Scintillators")
        self._wave_module = WavePanel("HyCal Module")
        wrl.addWidget(self._wave_scint)
        wrl.addWidget(self._wave_module)
        v_split.addWidget(wave_row)

        v_split.setStretchFactor(0, 1)
        v_split.setStretchFactor(1, 1)

        h_split.addWidget(left_scroll)
        h_split.addWidget(v_split)
        h_split.setStretchFactor(0, 2)
        h_split.setStretchFactor(1, 3)

        self.setCentralWidget(h_split)

    # Style helpers

    def _groupbox_style(self) -> str:
        return themed(
            f"QGroupBox{{color:{THEME.TEXT_DIM};font-size:11px;font-weight:bold;"
            f"border:1px solid {THEME.BORDER};border-radius:6px;margin-top:6px;"
            f"padding-top:4px;}}"
            f"QGroupBox::title{{subcontrol-origin:margin;left:8px;padding:0 4px;}}")

    def _input_style(self) -> str:
        return themed(
            f"QLineEdit,QDoubleSpinBox,QSpinBox{{background:{THEME.PANEL};"
            f"color:{THEME.TEXT};border:1px solid {THEME.BORDER};"
            f"border-radius:4px;padding:3px 6px;font-size:12px;}}"
            f"QLineEdit:focus,QDoubleSpinBox:focus,QSpinBox:focus{{"
            f"border-color:{THEME.ACCENT_BORDER};}}")

    @staticmethod
    def _dim_label(text: str) -> QLabel:
        lbl = QLabel(text)
        lbl.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        return lbl

    def _spin(self, lo: float, hi: float, value: float, step: float = 1,
              decimals: Optional[int] = None, suffix: str = ""):
        """A QDoubleSpinBox with ``decimals`` when given, else a QSpinBox,
        in the input style."""
        if decimals is None:
            sp = QSpinBox()
        else:
            sp = QDoubleSpinBox()
            sp.setDecimals(decimals)
        sp.setRange(lo, hi)
        sp.setValue(value)
        sp.setSingleStep(step)
        if suffix:
            sp.setSuffix(suffix)
        sp.setStyleSheet(self._input_style())
        return sp

    def _range_row(self, lo: QWidget, hi: QWidget) -> QHBoxLayout:
        """``lo`` "to" ``hi`` in one row."""
        row = QHBoxLayout()
        row.setSpacing(4)
        row.addWidget(lo)
        row.addWidget(self._dim_label("to"))
        row.addWidget(hi)
        return row

    def _btn_style(self, accent: bool = False) -> str:
        bg  = THEME.ACCENT_STRONG if accent else THEME.BUTTON
        hov = THEME.ACCENT        if accent else THEME.BUTTON_HOVER
        fg  = "#ffffff"           if accent else THEME.TEXT
        return themed(
            f"QPushButton{{background:{bg};color:{fg};border:1px solid "
            f"{THEME.BORDER};border-radius:6px;padding:5px 10px;font-size:12px;}}"
            f"QPushButton:hover{{background:{hov};}}"
            f"QPushButton:disabled{{background:{THEME.PANEL};"
            f"color:{THEME.TEXT_MUTED};}}")

    # Slots — connection & scan

    def _on_connect(self):
        url = self._url_edit.text().rstrip("/")
        self._server_url = url
        self._conn_label.setText("Connecting…")
        self._connect_btn.setEnabled(False)
        QApplication.processEvents()

        cfg = _http_get(f"{url}/api/config")
        if cfg is None:
            self._conn_label.setText("Connection failed")
            self._conn_label.setStyleSheet(f"color:{THEME.DANGER};font-size:11px;")
            self._connect_btn.setEnabled(True)
            return

        # Server connection takes priority — clear any local file selection
        # so Start uses the server, not the previously queued EVIO files.
        self._reset_file_queue()

        self._apply_config(cfg)
        self._connect_btn.setEnabled(True)

    def _on_add_evio_files(self):
        """Open a multi-select file dialog and add EVIO files to the queue."""
        paths, _ = QFileDialog.getOpenFileNames(
            self, "Add EVIO Files", "",
            "EVIO files (*.evio *.evio.*);;All files (*)")
        if not paths:
            return

        dec     = _prad2py.dec
        cfg_dec = dec.load_daq_config()
        ch      = dec.EvChannel()
        ch.set_config(cfg_dec)

        for path in paths:
            if path in self._local_evio_paths:
                continue
            if ch.open_auto(path) != dec.Status.success:
                ch.close()
                continue   # skip unreadable files silently
            n_records = (ch.get_random_access_event_count()
                         if ch.is_random_access() else 0)
            ch.close()

            self._local_evio_paths.append(path)
            fname = Path(path).name
            n_str = f"{n_records:,} rec" if n_records else "seq"
            item_text = f"{fname}  ({n_str})"
            from PyQt6.QtWidgets import QListWidgetItem
            item = QListWidgetItem(item_text)
            item.setToolTip(path)
            self._file_list.addItem(item)

        if not self._local_evio_paths:
            return

        n = len(self._local_evio_paths)
        self._local_file_label.setText(f"{n} file(s) queued")
        self._local_file_label.setStyleSheet(
            f"color:{THEME.SUCCESS};font-size:10px;")
        self._clear_files_btn.setEnabled(True)

        # Local file mode takes priority — clear any server state.
        self._n_events = 0
        self._et_mode  = False
        self._conn_label.setText("Using local files (no server)")
        self._conn_label.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")
        self._start_btn.setEnabled(True)
        self._fetch_btn.setEnabled(False)
        self._mode_btn.setEnabled(False)

    def _on_remove_evio_file(self):
        """Remove the currently selected file from the queue."""
        row = self._file_list.currentRow()
        if row < 0 or row >= len(self._local_evio_paths):
            return
        self._file_list.takeItem(row)
        self._local_evio_paths.pop(row)
        self._on_file_selection_changed()
        n = len(self._local_evio_paths)
        if n == 0:
            self._reset_file_queue()
            self._start_btn.setEnabled(
                self._n_events > 0 or self._et_mode)
        else:
            self._local_file_label.setText(f"{n} file(s) queued")

    def _on_clear_evio_files(self):
        """Remove all files from the queue."""
        self._reset_file_queue()
        self._start_btn.setEnabled(self._n_events > 0 or self._et_mode)
        self._conn_label.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")

    def _reset_file_queue(self):
        self._local_evio_paths.clear()
        self._file_list.clear()
        self._local_file_label.setText("No files selected")
        self._local_file_label.setStyleSheet(
            f"color:{THEME.TEXT_DIM};font-size:10px;")
        self._remove_file_btn.setEnabled(False)
        self._clear_files_btn.setEnabled(False)

    def _on_file_selection_changed(self):
        self._remove_file_btn.setEnabled(
            self._file_list.currentRow() >= 0)

    def _apply_config(self, cfg: dict) -> None:
        """Update UI state from a /api/config response dict."""
        n            = cfg.get("event_count", 0)
        mode         = cfg.get("mode", "unknown")
        et_connected = cfg.get("et_connected", False)
        self._n_events = n
        self._et_mode  = (mode == "online")

        if n > 0:
            self._conn_label.setText(f"Connected  ({n:,} events, mode={mode})")
            self._conn_label.setStyleSheet(f"color:{THEME.SUCCESS};font-size:11px;")
            self._start_btn.setEnabled(True)
            self._fetch_btn.setEnabled(True)
            self._ev_spin.setRange(1, n)
        elif self._et_mode and et_connected:
            self._conn_label.setText("Connected — ET online (live accumulation)")
            self._conn_label.setStyleSheet(f"color:{THEME.SUCCESS};font-size:11px;")
            self._start_btn.setEnabled(True)
            self._fetch_btn.setEnabled(True)   # fetches latest ring event
        elif self._et_mode:
            self._conn_label.setText("Connected — ET online (waiting for DAQ)")
            self._conn_label.setStyleSheet(f"color:{THEME.WARN};font-size:11px;")
            self._start_btn.setEnabled(True)
            self._fetch_btn.setEnabled(False)
        else:
            self._conn_label.setText(f"Connected — no file loaded (mode={mode})")
            self._conn_label.setStyleSheet(f"color:{THEME.WARN};font-size:11px;")
            self._start_btn.setEnabled(False)
            self._fetch_btn.setEnabled(False)

        self._mode_btn.setEnabled(True)
        if self._et_mode:
            self._mode_btn.setText("Switch to File Mode")
        else:
            self._mode_btn.setText("Switch to ET Mode")

    def _any_worker_running(self) -> bool:
        return (
            (self._stats_worker   is not None and self._stats_worker.isRunning()) or
            (self._instant_worker is not None and self._instant_worker.isRunning())
        )

    def _stop_workers_and_wait(self, ms: int) -> None:
        if self._any_worker_running():
            self._stop_worker()
            if self._stats_worker:   self._stats_worker.wait(ms)
            if self._instant_worker: self._instant_worker.wait(ms)

    def _on_mode_switch(self):
        """Toggle the server between online (ET) and file mode."""
        self._stop_workers_and_wait(2000)

        endpoint = "/api/mode/file" if self._et_mode else "/api/mode/online"
        self._mode_btn.setEnabled(False)
        self._mode_btn.setText("Switching…")
        QApplication.processEvents()

        result = _http_get(f"{self._server_url}{endpoint}")
        if result is None:
            self._conn_label.setText("Mode switch failed")
            self._conn_label.setStyleSheet(f"color:{THEME.DANGER};font-size:11px;")
            self._mode_btn.setEnabled(True)
            self._mode_btn.setText(
                "Switch to File Mode" if self._et_mode else "Switch to ET Mode")
            return

        cfg = _http_get(f"{self._server_url}/api/config")
        if cfg:
            self._apply_config(cfg)

    def _active_scint(self) -> str:
        btn = self._scint_group.checkedButton()
        return btn.text() if btn else "V1"


    def _on_display_mode_changed(self):
        instant = self._rb_disp_instant.isChecked()
        self._display_mode = DISPLAY_INSTANT if instant else DISPLAY_COINC

        # Pause/Resume only active in Instant mode while the instant worker runs.
        worker_live = self._instant_worker is not None and self._instant_worker.isRunning()
        if instant and worker_live:
            self._set_pause_buttons(not self._display_paused,
                                    self._display_paused)
        else:
            self._display_paused = False
            self._set_pause_buttons(False, False)

        # Map view (Coinc Rate / Occupancy) is only meaningful in stats mode.
        # AND/OR mode stays visible because it governs the background stats worker
        # that accumulates in both display modes.
        self._mapview_box.setVisible(not instant)

        # Switch the map view without stopping workers or resetting stats.
        if instant:
            self._wave_module.clear("Max ADC Module")
            # If no live instant data yet, blank the map until the next event.
            if not (self._instant_worker and self._instant_worker.isRunning()):
                self._map.set_values({})
                self._map.update()
        else:
            self._wave_module.clear("HyCal Module")
            # Restore the accumulated stats immediately.
            if self._snapshot:
                self._map.set_snapshot(self._snapshot, self._active_scint(),
                                       self._map_view)
            else:
                self._map.set_values({})
                self._map.update()

    def _on_scint_changed(self):
        if self._snapshot and self._display_mode == DISPLAY_COINC:
            self._map.set_snapshot(self._snapshot, self._active_scint(), self._map_view)
            self._update_stats_label(self._snapshot)
        # refresh scintillator waveform panel title & re-fetch if event is available
        if self._n_events > 0 or self._et_mode:
            self._wave_scint.clear("Scintillators")
            self._on_fetch_waveforms()

    def _on_mapview_changed(self):
        self._map_view = VIEW_OCC if self._rb_occ.isChecked() else VIEW_COINC
        if self._snapshot:
            self._map.set_snapshot(self._snapshot, self._active_scint(), self._map_view)
            self._update_stats_label(self._snapshot)

    def _on_start_stop(self):
        if self._any_worker_running():
            self._stop_worker()
        else:
            self._start_worker()

    def _on_pause(self):
        self._display_paused = True
        self._set_pause_buttons(False, True)

    def _on_resume(self):
        self._display_paused = False
        self._set_pause_buttons(True, False)

    def _set_pause_buttons(self, pause: bool, resume: bool) -> None:
        self._pause_btn.setEnabled(pause)
        self._resume_btn.setEnabled(resume)

    def _set_run_controls_enabled(self, on: bool) -> None:
        """Cut and processing settings are locked while a run is going."""
        for w in (self._scint_thr_spin, self._scint_tcut_min,
                  self._scint_tcut_max, self._hycal_thr_spin,
                  self._hycal_tcut_min, self._hycal_tcut_max,
                  self._max_lm_spin, self._cpu_spin, self._excl_spin,
                  self._min_mods_spin, self._max_rate_spin,
                  self._rb_and, self._rb_or, self._save_wfm_chk,
                  self._wfm_max_spin, self._step_chk):
            w.setEnabled(on)

    def _on_screenshot(self):
        import datetime
        ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
        default_path = str(Path.cwd() / f"coinc_monitor_{ts}.png")
        path, _ = QFileDialog.getSaveFileName(
            self, "Save screenshot", default_path,
            "PNG image (*.png);;All files (*)")
        if not path:
            return
        if not path.lower().endswith(".png"):
            path += ".png"
        pixmap = self.grab()
        if pixmap.save(path, "PNG"):
            self._status_label.setText(f"Screenshot saved → {path}")
        else:
            QMessageBox.warning(self, "Screenshot failed",
                                f"Could not write PNG to:\n{path}")

    def _start_worker(self):
        self._snapshot = {}
        self._display_paused = False
        self._map.set_values({})
        self._map.update()

        self._connect_btn.setEnabled(False)
        self._mode_btn.setEnabled(False)
        self._browse_btn.setEnabled(False)
        self._start_btn.setText("Stop")
        self._start_btn.setStyleSheet(self._btn_style(accent=False))

        # Apply inner-layer exclusion: drop W modules whose layer ≤ excl_layers.
        excl = self._excl_spin.value()
        if excl > 0:
            active_mod_keys = {
                name: key for name, key in self._mod_keys.items()
                if not name.startswith('W') or self._w_layers.get(name, 0) > excl
            }
        else:
            active_mod_keys = self._mod_keys

        sel_mode = MODE_AND if self._rb_and.isChecked() else MODE_OR
        cuts = CoincCuts(
            scint_thr=self._scint_thr_spin.value(),
            hycal_thr=self._hycal_thr_spin.value(),
            mode=sel_mode,
            min_mods=self._min_mods_spin.value(),
            scint_t_min=self._scint_tcut_min.value(),
            scint_t_max=self._scint_tcut_max.value(),
            hycal_t_min=self._hycal_tcut_min.value(),
            hycal_t_max=self._hycal_tcut_max.value(),
            neighbor_map=self._neighbor_map,
            max_lm=self._max_lm_spin.value(),
        )

        # --- Waveform collector (optional) ---
        wfm_coll = wfm_path = None
        if self._save_wfm_chk.isChecked():
            import datetime
            ts = datetime.datetime.now().strftime("%Y%m%d_%H%M%S")
            wfm_path = Path(f"coinc_wfm_{ts}.npz")
            wfm_coll = WaveformCollector(
                server_url=self._server_url,
                max_records=self._wfm_max_spin.value(),
            )
            self._wfm_label.setText(f"Will save to: {wfm_path}")

        # --- Stats worker (always) ---
        common = dict(module_keys=active_mod_keys, scint_keys=self._scint_keys,
                      cuts=cuts, wfm_collector=wfm_coll,
                      wfm_save_path=wfm_path, parent=self)
        if self._local_evio_paths:
            self._stats_worker = ProcessWorkerLocal(
                evio_paths=self._local_evio_paths,
                step_through=self._step_chk.isChecked(),
                n_workers=self._cpu_spin.value(),
                **common)
            self._progress.setRange(0, 100)
        elif self._et_mode:
            self._stats_worker = ProcessWorkerET(
                server_url=self._server_url,
                max_rate_hz=float(self._max_rate_spin.value()),
                **common)
            self._progress.setRange(0, 0)
        else:
            self._stats_worker = ProcessWorker(
                server_url=self._server_url,
                n_events=self._n_events,
                **common)
            self._progress.setRange(0, 100)

        self._stats_worker.progress.connect(self._on_progress)
        self._stats_worker.stats_update.connect(self._on_stats_update)
        self._stats_worker.finished.connect(self._on_finished)
        self._stats_worker.waveforms_saved.connect(self._on_waveforms_saved)
        if isinstance(self._stats_worker, ProcessWorkerLocal):
            self._stats_worker.coincidence_event.connect(self._on_coincidence_event)
            self._stats_worker.workers_setup.connect(self._setup_worker_prog_rows)
            self._stats_worker.worker_progress.connect(self._on_worker_progress)
        self._stats_worker.start()

        # Instant display worker: ET mode only (not available for local files).
        if self._et_mode and not self._local_evio_paths:
            self._instant_worker = InstantDisplayWorker(
                server_url=self._server_url,
                module_keys=active_mod_keys,
                hycal_thr=cuts.hycal_thr,
                min_mods=cuts.min_mods,
                parent=self,
            )
            self._instant_worker.event_ready.connect(self._on_instant_event)
            self._instant_worker.finished.connect(self._on_finished)
            self._instant_worker.start()
        else:
            self._instant_worker = None

        # Pause/Resume only meaningful in Instant Event Display mode (ET).
        self._set_pause_buttons(
            self._et_mode and not self._local_evio_paths
            and self._display_mode == DISPLAY_INSTANT, False)

        n_local = len(self._local_evio_paths)
        src = (f"{n_local} local file(s)" if self._local_evio_paths
               else ("ET live" if self._et_mode else "file"))
        self._status_label.setText(
            f"{'Accumulating' if self._et_mode else 'Processing'}… "
            f"({src}, mode: {sel_mode})")
        self._set_run_controls_enabled(False)

    def _stop_worker(self):
        if self._stats_worker:   self._stats_worker.stop()
        if self._instant_worker: self._instant_worker.stop()
        self._status_label.setText("Stopping…")
        self._start_btn.setEnabled(False)

    # Per-worker progress (parallel local-EVIO mode)

    def _setup_worker_prog_rows(self, n_workers: int):
        self._clear_worker_prog_rows()
        if n_workers <= 0:
            return
        for i in range(n_workers):
            row = QWidget()
            rh  = QHBoxLayout(row)
            rh.setSpacing(4)
            rh.setContentsMargins(0, 0, 0, 0)

            lbl = QLabel(f"W{i}: waiting…")
            lbl.setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:10px;")
            lbl.setMinimumWidth(180)
            lbl.setMaximumWidth(220)
            lbl.setWordWrap(False)
            lbl.setToolTip("")
            rh.addWidget(lbl, 0)

            bar = QProgressBar()
            bar.setRange(0, 0)   # indeterminate until first update
            bar.setValue(0)
            bar.setTextVisible(True)
            bar.setFormat("%p%")
            bar.setFixedHeight(12)
            bar.setStyleSheet(
                f"QProgressBar{{background:{THEME.PANEL};border:1px solid "
                f"{THEME.BORDER};border-radius:3px;font-size:9px;}}"
                f"QProgressBar::chunk{{background:{THEME.ACCENT};"
                f"border-radius:2px;}}")
            rh.addWidget(bar, 1)

            self._worker_prog_layout.addWidget(row)
            self._worker_prog_rows.append((lbl, bar))
        self._worker_prog_container.setVisible(True)

    def _clear_worker_prog_rows(self):
        while self._worker_prog_layout.count():
            item = self._worker_prog_layout.takeAt(0)
            w = item.widget()
            if w is not None:
                w.deleteLater()
        self._worker_prog_rows = []
        self._worker_prog_container.setVisible(False)

    def _on_worker_progress(self, info: dict):
        wid = int(info.get("worker_id", -1))
        if wid < 0 or wid >= len(self._worker_prog_rows):
            return
        lbl, bar = self._worker_prog_rows[wid]
        n_files   = int(info.get("n_files", 0))
        fidx      = int(info.get("file_idx", -1))
        basename  = str(info.get("file_basename", ""))
        rec_done  = int(info.get("records_done", 0))
        rec_tot   = int(info.get("records_total", 0))
        finished  = bool(info.get("finished", False))

        if finished:
            lbl.setText(f"W{wid}: done")
            bar.setRange(0, 100)
            bar.setValue(100)
            return

        if fidx >= 0 and basename:
            lbl.setText(f"W{wid} [{fidx + 1}/{n_files}] {basename}")
            lbl.setToolTip(basename)
        else:
            lbl.setText(f"W{wid}: starting…")

        if rec_tot > 0:
            bar.setRange(0, rec_tot)
            bar.setValue(min(rec_done, rec_tot))
        else:
            bar.setRange(0, 0)   # indeterminate
            bar.setValue(0)

    def _on_progress(self, processed: int, total: int):
        if self._display_mode == DISPLAY_INSTANT:
            return   # status label owned by _on_instant_event in this mode
        if total < 0:   # ET mode — no total
            self._status_label.setText(f"Accumulated {processed:,} events")
        elif total == 0:   # local EVIO sequential (unknown total)
            self._status_label.setText(f"Processed {processed:,} events…")
        else:
            pct = min(100, int(100 * processed / total))
            self._progress.setValue(pct)
            self._status_label.setText(f"Processed {processed:,} / {total:,}")

    def _on_stats_update(self, snap: dict):
        self._snapshot = snap
        # Only repaint the map when the user is viewing coinc/occ stats.
        if self._display_mode == DISPLAY_COINC:
            self._map.set_snapshot(snap, self._active_scint(), self._map_view)
        self._update_stats_label(snap)

    def _update_stats_label(self, snap: dict):
        scint = self._active_scint()
        mhits = snap.get("module_hits", {})
        processed = snap.get("processed", 0)
        if self._map_view == VIEW_OCC:
            total_hits = sum(mhits.values())
            max_hits   = max(mhits.values()) if mhits else 0
            nonzero    = sum(1 for v in mhits.values() if v > 0)
            mean_hits  = total_hits / nonzero if nonzero else 0.0
            self._stats_label.setText(
                f"Events processed: {processed:,}\n"
                f"Module hits (total): {total_hits:,}\n"
                f"Non-zero modules: {nonzero}\n"
                f"Mean hits: {mean_hits:.1f}\n"
                f"Max hits:  {max_hits:,}"
            )
        else:
            shits     = snap.get("scint_hits",     {}).get(scint, 0)
            scint_hits_any = snap.get("scint_hits_any", {})
            rates = snap.get("rates", {}).get(scint, {})
            valid = [v for v in rates.values() if not math.isnan(v) and v > 0]
            mean_rate = sum(valid) / len(valid) if valid else 0.0
            max_rate  = max(valid)              if valid else 0.0
            # Show all 4 scintillators' individual (pre-AND) hit counts at once
            # so the user can immediately spot which scintillator is not firing.
            indiv_parts = "  ".join(
                f"{s}={scint_hits_any.get(s, 0):,}" for s in SCINTILLATORS
            )
            and_label = (f"{scint} AND HyCal hits: {shits:,}"
                         if self._rb_and.isChecked()
                         else f"{scint} hits: {shits:,}")
            self._stats_label.setText(
                f"Indiv. hits: {indiv_parts}\n"
                f"{and_label}\n"
                f"Module hits (total): {sum(mhits.values()):,}\n"
                f"Non-zero modules: {len(valid)}\n"
                f"Mean rate: {mean_rate:.4f}\n"
                f"Max rate:  {max_rate:.4f}"
            )

    def _on_coincidence_event(self, data: dict):
        """Step-through mode: show waveforms for a single coincidence event and wait."""
        sf         = data["scint_fired"]
        scint_name = self._active_scint()

        # Only pause when the currently-selected scintillator itself fired in-window.
        if not sf.get(scint_name, False):
            if isinstance(self._stats_worker, ProcessWorkerLocal):
                self._stats_worker.step_continue()
            return

        ev_num  = data["event_number"]
        best    = data["best"]
        best_adc = data.get("best_adc", 0.0)
        wfm_ch  = data["wfm_channels"]
        t_min = self._scint_tcut_min.value()
        t_max = self._scint_tcut_max.value()
        scint_wfms = {
            sname: wfm_ch.get(skey, {"error": "no waveform"})
            for sname, skey in self._scint_keys.items()
        }
        self._wave_scint.set_multi_data(
            scint_wfms, self._scint_thr_spin.value(),
            list(SCINTILLATORS),
            f"Scintillators  (event {ev_num})",
            t_min=t_min, t_max=t_max)

        best_key   = self._mod_keys.get(best, "")
        module_wfm = wfm_ch.get(best_key, {"error": "no waveform"})
        self._wave_module.set_data(
            module_wfm, self._hycal_thr_spin.value(),
            f"{best}  ADC={best_adc:.0f}  (event {ev_num})")

        fired_names = [s for s, f in sf.items() if f]
        self._status_label.setText(
            f"[step] Event {ev_num}  |  best: {best}  |  fired: "
            f"{', '.join(fired_names) if fired_names else 'none'}\n"
            f"Click 'Continue →' for next coincidence.")
        self._step_continue_btn.setEnabled(True)

    def _on_step_continue(self):
        self._step_continue_btn.setEnabled(False)
        if isinstance(self._stats_worker, ProcessWorkerLocal):
            self._stats_worker.step_continue()

    def _on_instant_event(self, data: dict):
        """Handle a new live event from InstantDisplayWorker."""
        if self._display_mode != DISPLAY_INSTANT:
            return   # accumulating in background; don't touch the map or panels
        if self._display_paused:
            return   # paused: freeze display while stats keep accumulating

        adc_vals = data.get("adc_vals", {})
        channels = data.get("channels", {})
        seq      = data.get("seq", 0)
        max_mod  = data.get("max_mod", "")

        self._map.set_instant_event(adc_vals, max_mod)

        scint_wfms = {}
        for sname, skey in self._scint_keys.items():
            ch = channels.get(skey, {}) if skey else {}
            scint_wfms[sname] = ch if ch else {"error": "no signal"}
        self._wave_scint.set_multi_data(
            scint_wfms, self._scint_thr_spin.value(),
            list(SCINTILLATORS),
            f"Scintillators  (seq {seq})")

        if max_mod:
            max_key = self._mod_keys.get(max_mod, "")
            max_ch  = channels.get(max_key, {}) if max_key else {}
            max_adc = adc_vals.get(max_mod, 0.0)
            self._wave_module.set_data(
                max_ch if max_ch else {"error": "no waveform"},
                self._hycal_thr_spin.value(),
                f"{max_mod}  ADC={max_adc:.0f}  (seq {seq})")

        self._status_label.setText(
            f"Seq {seq}  |  max: {max_mod}  {adc_vals.get(max_mod, 0):.0f}")

    def _on_finished(self, msg: str):
        # Both workers emit finished; wait until neither is running.
        if self._any_worker_running():
            return
        self._display_paused = False
        self._set_pause_buttons(False, False)
        self._progress.setRange(0, 100)
        self._progress.setValue(100 if msg == "" else self._progress.value())
        self._clear_worker_prog_rows()
        self._start_btn.setText("Start")
        self._start_btn.setStyleSheet(self._btn_style(accent=True))
        self._start_btn.setEnabled(True)
        self._connect_btn.setEnabled(True)
        self._mode_btn.setEnabled(True)
        if _HAVE_PRAD2PY:
            self._browse_btn.setEnabled(True)
            self._remove_file_btn.setEnabled(self._file_list.currentRow() >= 0)
            self._clear_files_btn.setEnabled(len(self._local_evio_paths) > 0)
        self._set_run_controls_enabled(True)
        self._step_continue_btn.setEnabled(False)
        if msg == "stopped":
            self._status_label.setText("Stopped.")
        else:
            n = self._snapshot.get("processed", 0)
            self._status_label.setText(f"Done — {n:,} events processed.")

    def _on_waveforms_saved(self, count: int, path: str):
        self._wfm_label.setText(f"Saved {count} waveform pairs → {path}")

    # Slots — waveform browser

    def _on_module_clicked(self, name: str):
        if not name:
            return
        if name in SCINTILLATORS:
            # Scintillator clicked on map → select its radio button and refresh left panel
            for btn in self._scint_group.buttons():
                if btn.text() == name:
                    btn.setChecked(True)
                    break
            self._on_scint_changed()
        else:
            # HyCal module clicked → update right waveform panel
            self._selected_module = name
            self._sel_mod_label.setText(f"Module: {name}")
            self._on_fetch_waveforms()

    def _on_fetch_waveforms(self):
        if self._n_events == 0 and not self._et_mode:
            return

        if self._fetcher and self._fetcher.isRunning():
            self._fetcher.wait(500)

        ev         = self._ev_spin.value()
        module_key = self._mod_keys.get(self._selected_module, "")

        ev_label = "latest ET event" if self._et_mode else f"event {ev}"
        self._wave_scint.clear("Scintillators", f"Fetching {ev_label}…")
        self._wave_module.clear(
            self._selected_module or "HyCal Module",
            f"Fetching {ev_label}…" if self._selected_module
            else "Click a module on the map")

        self._fetcher = WaveformFetcher(
            self._server_url, ev, self._scint_keys, module_key,
            et_mode=self._et_mode, parent=self)

        self._fetcher.waveform_ready.connect(self._on_waveform_ready)
        self._fetcher.start()

    def _on_waveform_ready(self, label: str, data: dict):
        ev = self._ev_spin.value()
        if label == "scint":
            t_min = self._scint_tcut_min.value()
            t_max = self._scint_tcut_max.value()
            title = f"Scintillators  (event {ev})"
            self._wave_scint.set_multi_data(
                data, self._scint_thr_spin.value(),
                list(SCINTILLATORS), title,
                t_min=t_min, t_max=t_max)
        else:
            mod_name = self._selected_module or "HyCal Module"
            title = f"{mod_name}  (event {ev})"
            self._wave_module.set_data(
                data, self._hycal_thr_spin.value(), title)

    # Veto position polling

    def _poll_veto_positions(self):
        for vname in VETO_PV_BASES:
            rbv  = self._veto_ctrl.get(f"{vname}_rbv")
            movn = self._veto_ctrl.get(f"{vname}_movn")
            w = self._veto_widgets[vname]

            if rbv is not None:
                w["rbv"].setText(f"{rbv:.2f}")
                w["rbv"].setStyleSheet(f"color:{THEME.TEXT};font-size:11px;")
            else:
                w["rbv"].setText("—")
                w["rbv"].setStyleSheet(f"color:{THEME.TEXT_DIM};font-size:11px;")

            moving = movn is not None and int(movn) == 1
            w["movn"].setText("● moving" if moving else "")
            w["movn"].setStyleSheet(
                f"color:{THEME.SUCCESS};font-size:10px;" if moving
                else "color:transparent;font-size:10px;")

    def _move_veto(self, vname: str):
        """Write the spinbox setpoint to the motor VAL PV."""
        target = self._veto_widgets[vname]["spin"].value()
        ok = self._veto_ctrl.put(f"{vname}_val", target)
        if not ok:
            self._status_label.setText(
                f"{vname} move failed — PV not connected")

    # Cleanup

    def closeEvent(self, event):
        self._veto_timer.stop()
        self._stop_workers_and_wait(3000)
        if self._fetcher and self._fetcher.isRunning():
            self._fetcher.wait(1000)
        event.accept()


# Entry point

def main():
    parser = argparse.ArgumentParser(
        description="Scintillator–HyCal coincidence rate monitor")
    parser.add_argument("--url", default=DEFAULT_URL,
                        help=f"prad2_server base URL (default: {DEFAULT_URL})")
    parser.add_argument("--theme", choices=available_themes(), default="dark")
    args = parser.parse_args()

    app = QApplication(sys.argv)
    app.setApplicationName("Coincidence Monitor")

    set_theme(args.theme)

    win = MainWindow(server_url=args.url)
    apply_theme_palette(win)
    win.setStyleSheet(themed(f"QMainWindow{{background:{THEME.BG};}}"))
    win.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
