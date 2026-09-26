"""
EPICS interface layer for HyCal calibration tools.

Provides base PV wrappers, motor-control classes, scaler readers, and
simulated stand-ins so every calibration GUI can run offline without
a live EPICS environment.

Tier 1 -- base utilities:
    PVGroup, ScalerPVGroup

Tier 2 -- application-specific:
    SPMG, PV (constants), MOTOR_PV_MAP,
    MotorEPICS (read-only observer when writable=False),
    SimulatedMotorEPICS, SimulatedScalerEPICS,
    epics_move_to, epics_stop, epics_pause, epics_resume,
    epics_is_moving, epics_read_rbv, epics_wait_move_done
"""

from __future__ import annotations

import math
import random
import threading
import time
from enum import IntEnum
from typing import Any, Callable, Dict, List, Tuple

from scan_utils import (
    Module,
    BEAM_CENTER_X,
    BEAM_CENTER_Y,
    SCANNABLE_TYPES,
    ptrans_in_limits,
)


# -- Tier 1: base utility classes ------------------------------------------

class PVGroup:
    """Thin wrapper around *pyepics* PV objects.

    Parameters
    ----------
    pv_map : list of (key, pv_name) pairs
        Friendly key to EPICS PV name mapping.
    writable : bool
        If False, :meth:`put` is a no-op that always returns False.
    timeout : float
        Connection timeout in seconds passed to each ``epics.PV``.
    """

    def __init__(
        self,
        pv_map: List[Tuple[str, str]],
        writable: bool = False,
        timeout: float = 5.0,
    ) -> None:
        self._pv_map = list(pv_map)
        self._writable = writable
        self._timeout = timeout
        self._pvs: Dict[str, Any] = {}
        self._all_connected = False

    def connect(self) -> Tuple[int, int]:
        """Create PV objects, wait for connections, return (connected, total)."""
        import epics as _epics
        for key, pvname in self._pv_map:
            self._pvs[key] = _epics.PV(pvname, connection_timeout=self._timeout)
        time.sleep(2.0)
        n = sum(1 for p in self._pvs.values() if p.connected)
        self._all_connected = n == len(self._pvs)
        return n, len(self._pvs)

    def get(self, key: str, default: Any = None) -> Any:
        """Read a PV by its friendly *key*."""
        pv = self._pvs.get(key)
        if pv and pv.connected:
            v = pv.get()
            return v if v is not None else default
        return default

    def put(self, key: str, value: Any) -> bool:
        """Write *value* to the PV identified by *key*.

        Returns True on success, False if not writable or not connected.
        """
        if not self._writable or not self._all_connected:
            return False
        pv = self._pvs.get(key)
        if pv and pv.connected:
            pv.put(value)
            return True
        return False

    def disconnected_pvs(self) -> List[str]:
        """Return PV *names* (not keys) that failed to connect."""
        return [
            pvname
            for key, pvname in self._pv_map
            if key in self._pvs and not self._pvs[key].connected
        ]

    def stop(self) -> None:  # noqa: D401
        """No-op in the base class; subclasses may override."""


class ScalerPVGroup(PVGroup):
    """Reads FADC scaler PVs for PbWO4 and PbGlass modules.

    PV pattern: ``B_DET_HYCAL_FADC_{label}:c`` where *label* is the
    module name (e.g. ``W1``, ``G100``).
    """

    _PATTERN = "B_DET_HYCAL_FADC_{}:c"

    def __init__(self, modules: List[Module]) -> None:
        super().__init__([(m.name, self._PATTERN.format(m.name))
                          for m in modules if m.mod_type in SCANNABLE_TYPES])

    def get_all(self) -> Dict[str, float]:
        """Batch-read all connected scaler PVs."""
        result: Dict[str, float] = {}
        for label, _ in self._pv_map:
            v = self.get(label)
            if v is not None:
                result[label] = float(v)
        return result


# -- Tier 2: application-specific constants --------------------------------

class SPMG(IntEnum):
    STOP = 0; PAUSE = 1; MOVE = 2; GO = 3  # noqa: E702


SPMG_LABELS = {0: "Stop", 1: "Pause", 2: "Move", 3: "Go"}


class PV:
    """EPICS PV name constants for the HyCal transporter motors."""
    X_VAL  = "ptrans_x.VAL";   Y_VAL  = "ptrans_y.VAL"          # noqa: E222,E702
    X_SPMG = "ptrans_x.SPMG";  Y_SPMG = "ptrans_y.SPMG"        # noqa: E222,E702
    X_ENCODER = "hallb_ptrans_x_encoder"
    Y_ENCODER = "hallb_ptrans_y1_encoder"
    BEAM_CUR = "hallb_IPM2C21A_CUR"
    X_RBV  = "ptrans_x.RBV";   Y_RBV  = "ptrans_y.RBV"          # noqa: E222,E702
    X_MOVN = "ptrans_x.MOVN";  Y_MOVN = "ptrans_y.MOVN"         # noqa: E222,E702
    X_VELO = "ptrans_x.VELO";  Y_VELO = "ptrans_y.VELO"         # noqa: E222,E702
    X_ACCL = "ptrans_x.ACCL";  Y_ACCL = "ptrans_y.ACCL"         # noqa: E222,E702
    X_TDIR = "ptrans_x.TDIR";  Y_TDIR = "ptrans_y.TDIR"         # noqa: E222,E702
    X_MSTA = "ptrans_x.MSTA";  Y_MSTA = "ptrans_y.MSTA"         # noqa: E222,E702
    X_ATHM = "ptrans_x.ATHM";  Y_ATHM = "ptrans_y.ATHM"         # noqa: E222,E702
    X_PREC = "ptrans_x.PREC";  Y_PREC = "ptrans_y.PREC"         # noqa: E222,E702
    X_BVEL = "ptrans_x.BVEL";  Y_BVEL = "ptrans_y.BVEL"         # noqa: E222,E702
    X_BACC = "ptrans_x.BACC";  Y_BACC = "ptrans_y.BACC"         # noqa: E222,E702
    X_VBAS = "ptrans_x.VBAS";  Y_VBAS = "ptrans_y.VBAS"         # noqa: E222,E702
    X_BDST = "ptrans_x.BDST";  Y_BDST = "ptrans_y.BDST"         # noqa: E222,E702
    X_FRAC = "ptrans_x.FRAC";  Y_FRAC = "ptrans_y.FRAC"         # noqa: E222,E702


MOTOR_PV_MAP: List[Tuple[str, str]] = [
    ("x_val", PV.X_VAL), ("y_val", PV.Y_VAL),
    ("x_spmg", PV.X_SPMG), ("y_spmg", PV.Y_SPMG),
    ("x_encoder", PV.X_ENCODER), ("y_encoder", PV.Y_ENCODER),
    ("beam_cur", PV.BEAM_CUR),
    ("x_rbv", PV.X_RBV), ("y_rbv", PV.Y_RBV),
    ("x_movn", PV.X_MOVN), ("y_movn", PV.Y_MOVN),
    ("x_velo", PV.X_VELO), ("y_velo", PV.Y_VELO),
    ("x_accl", PV.X_ACCL), ("y_accl", PV.Y_ACCL),
    ("x_tdir", PV.X_TDIR), ("y_tdir", PV.Y_TDIR),
    ("x_msta", PV.X_MSTA), ("y_msta", PV.Y_MSTA),
    ("x_athm", PV.X_ATHM), ("y_athm", PV.Y_ATHM),
    ("x_prec", PV.X_PREC), ("y_prec", PV.Y_PREC),
    ("x_bvel", PV.X_BVEL), ("y_bvel", PV.Y_BVEL),
    ("x_bacc", PV.X_BACC), ("y_bacc", PV.Y_BACC),
    ("x_vbas", PV.X_VBAS), ("y_vbas", PV.Y_VBAS),
    ("x_bdst", PV.X_BDST), ("y_bdst", PV.Y_BDST),
    ("x_frac", PV.X_FRAC), ("y_frac", PV.Y_FRAC),
]


# -- Tier 2: motor EPICS classes -------------------------------------------

class MotorEPICS(PVGroup):
    """Live EPICS interface for the HyCal transporter motors.

    Extends :class:`PVGroup` with :data:`MOTOR_PV_MAP` and an emergency
    ``stop()`` that commands SPMG.STOP on both axes.  With
    ``writable=False`` it is a pure observer: ``put()`` and ``stop()``
    never touch the hardware.
    """

    def __init__(self, writable: bool = False, timeout: float = 5.0) -> None:
        super().__init__(MOTOR_PV_MAP, writable=writable, timeout=timeout)

    def stop(self) -> None:
        """Emergency stop -- set SPMG to STOP on both axes.

        Bypasses the all-connected check of :meth:`put` so a partially
        connected motor can still be stopped.
        """
        if not self._writable:
            return
        for key in ("x_spmg", "y_spmg"):
            pv = self._pvs.get(key)
            if pv and pv.connected:
                pv.put(int(SPMG.STOP))


# -- Tier 2: simulated stand-ins -------------------------------------------

class SimulatedMotorEPICS:
    """Offline motor simulator with the same interface as :class:`MotorEPICS`."""

    def __init__(self) -> None:
        self._lock = threading.Lock()
        self._x = BEAM_CENTER_X
        self._y = BEAM_CENTER_Y
        self._tx = self._x
        self._ty = self._y
        self._x_spmg = int(SPMG.GO)
        self._y_spmg = int(SPMG.GO)
        self._x_movn = 0
        self._y_movn = 0
        self._x_speed = 0.5
        self._y_speed = 0.2
        self._moving = False

    def connect(self) -> Tuple[int, int]:
        return (0, 0)

    def disconnected_pvs(self) -> List[str]:
        return []

    def get(self, key: str, default: Any = None) -> Any:
        with self._lock:
            return {
                "x_encoder": round(self._x + random.gauss(0, 0.002), 4),
                "y_encoder": round(self._y + random.gauss(0, 0.002), 4),
                "x_rbv": round(self._x, 3), "y_rbv": round(self._y, 3),
                "x_val": round(self._tx, 3), "y_val": round(self._ty, 3),
                "x_movn": self._x_movn, "y_movn": self._y_movn,
                "x_spmg": self._x_spmg, "y_spmg": self._y_spmg,
                "x_velo": self._x_speed, "y_velo": self._y_speed,
                "x_accl": 0.2, "y_accl": 1.0,
                "x_tdir": 1 if self._tx >= self._x else 0,
                "y_tdir": 1 if self._ty >= self._y else 0,
                "x_msta": 0x10B, "y_msta": 0x10B,
                "x_athm": int(abs(self._x - BEAM_CENTER_X) < 1.0),
                "y_athm": int(abs(self._y - BEAM_CENTER_Y) < 1.0),
                "x_prec": 3, "y_prec": 3,
                "x_bvel": 1.0, "y_bvel": 1.0,
                "x_bacc": 1.0, "y_bacc": 1.0,
                "x_vbas": 0.5, "y_vbas": 0.5,
                "x_bdst": 0.0, "y_bdst": 0.0,
                "x_frac": 1.0, "y_frac": 1.0,
                "beam_cur": 50.0,
            }.get(key, default)

    def put(self, key: str, value: Any) -> bool:
        with self._lock:
            if key == "x_val":
                self._tx = float(value)
            elif key == "y_val":
                self._ty = float(value)
            elif key == "x_spmg":
                self._x_spmg = int(value)
            elif key == "y_spmg":
                self._y_spmg = int(value)
            else:
                return False
            self._evaluate_motion()
        return True

    def stop(self) -> None:
        self.put("x_spmg", int(SPMG.STOP))
        self.put("y_spmg", int(SPMG.STOP))

    def _evaluate_motion(self) -> None:
        halted = (SPMG.STOP, SPMG.PAUSE)
        if self._x_spmg in halted or self._y_spmg in halted:
            self._moving = False; self._x_movn = 0; self._y_movn = 0       # noqa: E702
        elif self._x_spmg == SPMG.GO and self._y_spmg == SPMG.GO:
            if not self._moving:
                self._moving = True
                threading.Thread(target=self._run_move, daemon=True).start()

    def _run_move(self) -> None:
        dt = 0.02
        while True:
            with self._lock:
                if not self._moving:
                    self._x_movn = 0; self._y_movn = 0                     # noqa: E702
                    return
                dx, dy = self._tx - self._x, self._ty - self._y
                if abs(dx) > 0.001:
                    self._x += math.copysign(min(self._x_speed * dt, abs(dx)), dx)
                if abs(dy) > 0.001:
                    self._y += math.copysign(min(self._y_speed * dt, abs(dy)), dy)
                self._x_movn = 1 if abs(self._tx - self._x) > 0.001 else 0
                self._y_movn = 1 if abs(self._ty - self._y) > 0.001 else 0
                if self._x_movn == 0 and self._y_movn == 0:
                    self._x, self._y = self._tx, self._ty
                    self._moving = False
                    return
            time.sleep(dt)


class SimulatedScalerEPICS:
    """Offline scaler simulator matching the :class:`ScalerPVGroup` interface.

    Uses a seeded :class:`random.Random` instance (seed=0) so values are
    reproducible across runs.
    """

    def __init__(self, modules: List[Module]) -> None:
        self._labels: List[str] = [
            m.name for m in modules if m.mod_type in SCANNABLE_TYPES
        ]
        self._rng = random.Random(0)

    def connect(self) -> Tuple[int, int]:
        return self.connection_count()

    def get_all(self) -> Dict[str, float]:
        return {label: self._rng.uniform(0, 1000) for label in self._labels}

    def connection_count(self) -> Tuple[int, int]:
        return len(self._labels), len(self._labels)


# -- Tier 2: helper functions ----------------------------------------------

def epics_move_to(ep: Any, x: float, y: float) -> bool:
    """Command the transporter to move to (*x*, *y*).

    Returns False if the target is outside travel limits.
    """
    if not ptrans_in_limits(x, y):
        return False
    ep.put("x_val", x); ep.put("y_val", y)                                 # noqa: E702
    ep.put("x_spmg", int(SPMG.GO)); ep.put("y_spmg", int(SPMG.GO))        # noqa: E702
    return True


def epics_stop(ep: Any) -> None:
    """Emergency stop on both axes."""
    ep.stop()


def epics_pause(ep: Any) -> None:
    """Pause motion on both axes."""
    ep.put("x_spmg", int(SPMG.PAUSE)); ep.put("y_spmg", int(SPMG.PAUSE))  # noqa: E702


def epics_resume(ep: Any) -> None:
    """Resume motion on both axes."""
    ep.put("x_spmg", int(SPMG.GO)); ep.put("y_spmg", int(SPMG.GO))        # noqa: E702


def epics_is_moving(ep: Any) -> bool:
    """True if either axis is currently in motion."""
    return bool(ep.get("x_movn", 0)) or bool(ep.get("y_movn", 0))


def epics_read_rbv(ep: Any) -> Tuple[float, float]:
    """Return the current (x, y) read-back values."""
    return (ep.get("x_rbv", 0.0), ep.get("y_rbv", 0.0))


def epics_wait_move_done(ep: Any, x: float, y: float, *,
                         pos_threshold: float, timeout: float,
                         aborted: Callable[[], bool],
                         hold_if_paused: Callable[[], None],
                         log: Callable[..., None]) -> bool:
    """Wait for the transporter to stop at (*x*, *y*).

    "On position" requires BOTH MOVN=0 AND RBV within *pos_threshold*
    (mm) of the target.  Checking only MOVN is unsafe because MOVN is
    still 0 in the brief window after issuing a move command before
    the IOC has processed it — a check in that window would declare
    the move "done" before it started.

    Polls every 0.1 s.  *aborted* ends the wait early; *hold_if_paused*
    blocks while the operator has paused the scan.  Returns False on
    abort or after *timeout* seconds (logged via *log*); the caller
    decides which happened.
    """
    t0 = time.time()
    while not aborted():
        hold_if_paused()
        if aborted():
            return False
        if not epics_is_moving(ep):
            rx, ry = epics_read_rbv(ep)
            if math.sqrt((rx - x) ** 2 + (ry - y) ** 2) <= pos_threshold:
                return True
        if time.time() - t0 > timeout:
            log(f"MOVE TIMEOUT after {timeout:.0f}s", level="error")
            return False
        time.sleep(0.1)
    return False
