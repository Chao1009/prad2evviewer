"""
Shared GEM strip mapping — thin wrapper over prad2py.det.map_strip /
map_apv_strips.

The 6-step pipeline that maps an APV25 channel index to a plane-wide strip
number lives in C++ (prad2det/src/GemSystem.cpp: gem::MapStrip) and is
shared between on-line reconstruction and these off-line scripts.

Build requirement: the ``prad2py`` pybind11 module must be built and
importable.  Configure with ``-DBUILD_PYTHON=ON`` and either install it
(``cmake --build build --target install``) or prepend ``build/python/`` to
``PYTHONPATH``; otherwise importing this module raises ImportError.
"""

from __future__ import annotations

import os as _os
import sys as _sys


_SCRIPTS = _os.path.join(_os.path.dirname(_os.path.realpath(__file__)), "..", "scripts")
if _SCRIPTS not in _sys.path:
    _sys.path.append(_SCRIPTS)
from prad2_env import PRAD2PY_HINT, import_prad2py  # noqa: E402


def _resolve_prad2py():
    """prad2py.det from PYTHONPATH, else from this checkout's build tree."""
    mod, err = import_prad2py(build_first=False)
    if mod is None:
        raise ImportError(f"gem_strip_map requires prad2py.det ({err}).\n" + PRAD2PY_HINT)
    return mod.det


_det = _resolve_prad2py()


def map_strip(ch, plane_index, orient, pin_rotate=0, shared_pos=-1,
              hybrid_board=True, apv_channels=128, readout_center=32):
    """Map APV channel to plane-wide strip number.

    Returns ``(local_strip, plane_strip)``; ``local_strip`` is the
    plane-wide strip with the plane offset undone.
    """
    plane = _det.map_strip(ch=ch,
                           plane_index=plane_index,
                           orient=orient,
                           pin_rotate=pin_rotate,
                           shared_pos=shared_pos,
                           hybrid_board=hybrid_board,
                           apv_channels=apv_channels,
                           readout_center=readout_center)

    # Mirrors the plane_shift math in gem::MapStrip (steps 4-6).
    eff_pos = shared_pos if shared_pos >= 0 else plane_index
    plane_shift = (eff_pos - plane_index) * apv_channels - pin_rotate
    local = plane - (plane_shift + plane_index * apv_channels)
    return local, plane


def map_apv_strips(apv, apv_channels=128, readout_center=32):
    """Map all channels of an APV entry (from gem_map.json) to plane strip
    numbers.

    Returns a list of length ``apv_channels`` — the plane-wide strip number
    for each APV channel.  Accepts the same dict shape the JSON loader
    produces (keys: ``pos``, ``orient``, optional ``pin_rotate``,
    ``shared_pos``, ``hybrid_board``).
    """
    return _det.map_apv_strips(
        plane_index=apv["pos"],
        orient=apv["orient"],
        pin_rotate=apv.get("pin_rotate", 0),
        shared_pos=apv.get("shared_pos", -1),
        hybrid_board=apv.get("hybrid_board", True),
        apv_channels=apv_channels,
        readout_center=readout_center,
    )
