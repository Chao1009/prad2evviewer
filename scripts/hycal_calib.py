"""Qt-free HyCal calibration helpers for the scripts/ viewers: elastic
e-p / Moller kinematics with the upstream energy loss (mirrors
analysis/src/PhysicsTools.cpp) and the LMS gain-factor table reader.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Dict, Optional, Tuple

M_PROTON = 938.272        # MeV
M_ELECTRON = 0.51099895   # MeV


def energy_loss(theta_deg: float) -> float:
    """Energy lost upstream of HyCal (MeV) at polar angle ``theta_deg``."""
    cos_t = math.cos(math.radians(theta_deg))
    sec = 1.0 / cos_t if cos_t > 0.01 else 100.0
    eloss = 0.500 * 1.6 * sec    # Al window
    eloss += 0.120 * 1.6 * sec   # GEM window Al foils (2 GEMs)
    eloss += 0.100 * 2.0 * sec   # GEM foils (2 GEMs)
    eloss += 0.480 * 1.8 * sec   # kapton cover
    return eloss


def expected_energy(theta_deg: float, ebeam: float, kind: str = "ep") -> float:
    """Scattered-electron energy (MeV) at HyCal for beam energy ``ebeam``
    (MeV): elastic e-p (``kind="ep"``) or Moller (``"ee"``), minus
    energy_loss() and clamped at 0; 0 for an unknown kind."""
    cos_t = math.cos(math.radians(theta_deg))
    if kind == "ep":
        e = ebeam * M_PROTON / (M_PROTON + ebeam * (1.0 - cos_t))
    elif kind == "ee":
        gamma = ebeam / M_ELECTRON
        num = (gamma + 1.0) + (gamma - 1.0) * cos_t * cos_t
        den = (gamma + 1.0) - (gamma - 1.0) * cos_t * cos_t
        if den <= 0.0:
            return 0.0
        e = M_ELECTRON * num / den
    else:
        return 0.0
    return max(0.0, e - energy_loss(theta_deg))


@dataclass
class LMSRecord:
    alpha_peak: float
    alpha_sigma: float
    alpha_chi2ndf: float
    lms_peak: float
    lms_sigma: float
    lms_chi2ndf: float


@dataclass
class ModuleRecord:
    lms_peak: float
    lms_sigma: float
    lms_chi2ndf: float
    gain_factors: Tuple[float, float, float]


def read_lms_dat(path) -> Optional[Tuple[Dict[str, LMSRecord],
                                         Dict[str, ModuleRecord]]]:
    """Read a prad_<run>_LMS.dat gain-factor table into
    ``(reference PMTs, modules)`` keyed by name; None if it cannot be read.

    Each row is a name and six numbers, separated by blanks or commas.
    LMS<n> rows are the reference PMTs with their LMS and alpha peak fits:
    refGain_produce writes a "Name ..." header line and the LMS fit first,
    gain_fitter writes no header and the alpha fit first.  Every other row
    is a module's LMS fit followed by its gain factors against reference
    PMT 1/2/3.  Rows that do not parse are skipped.  The C++ reader is
    prad2::LoadRefGainFile (analysis/include/gain_factor.h).
    """
    try:
        with open(path) as f:
            rows = [line.replace(",", " ").split() for line in f]
    except OSError:
        return None
    lms_first = bool(rows) and rows[0][:1] == ["Name"]
    lms: Dict[str, LMSRecord] = {}
    modules: Dict[str, ModuleRecord] = {}
    for parts in rows:
        if len(parts) < 7:
            continue
        try:
            v = [float(x) for x in parts[1:7]]
        except ValueError:
            continue
        name = parts[0]
        if name.startswith("LMS"):
            alpha, lmsfit = (v[3:], v[:3]) if lms_first else (v[:3], v[3:])
            lms[name] = LMSRecord(*alpha, *lmsfit)
        else:
            modules[name] = ModuleRecord(v[0], v[1], v[2], (v[3], v[4], v[5]))
    return lms, modules
