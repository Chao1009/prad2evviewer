#!/usr/bin/env python3
"""
nonlinearity_viewer.py — HyCal PWO4 non-linearity analysis viewer
==================================================================

Reads the ROOT output of prad2ana_nonlinearity:
  energy_3p5GeV/h_energy_{name}_3p5   per-module TH1F at 3.5 GeV
  energy_2p2GeV/h_energy_{name}_2p2   per-module TH1F at 2.2 GeV
  energy_0p7GeV/h_energy_{name}_0p7   per-module TH1F at 0.7 GeV

Left:  HyCal map (W modules only), coloured by non-linearity parameter nl.
Right: for the selected module, the e-p and e-e peaks at each beam energy
       (2x3 histogram grid) and the linearity plot (E_rec/E_exp vs E_rec),
       plus controls to re-fit individual peaks and re-run the nl fit.

The peak and nl fits follow analysis/tools/nonlinearity.cpp.

Usage:
  python scripts/nonlinearity_viewer.py [root_file] [--theme dark|light]
"""
from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

from PyQt6.QtCore import Qt, QTimer
from PyQt6.QtGui import QFont, QKeySequence, QShortcut
from PyQt6.QtWidgets import (
    QApplication, QCheckBox, QComboBox, QFileDialog, QGridLayout, QGroupBox,
    QHBoxLayout, QLabel, QLineEdit, QMainWindow, QPushButton,
    QSizePolicy, QSpinBox, QDoubleSpinBox, QSplitter, QVBoxLayout, QWidget,
)
import matplotlib
matplotlib.use("QtAgg")
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.figure import Figure
from matplotlib.widgets import SpanSelector

try:
    import uproot
    HAS_UPROOT = True
except ImportError:
    HAS_UPROOT = False

try:
    from scipy.optimize import curve_fit
    HAS_SCIPY = True
except ImportError:
    HAS_SCIPY = False

from hycal_calib import expected_energy
from hycal_geoview import (
    HyCalMapWidget, Module, load_modules,
    apply_theme_palette, set_theme, available_themes, THEME,
    ColorRangeControl, edge_depth, hole_ring,
)
from prad2_env import database_dir, find_database_file

MODULES_JSON = find_database_file("hycal_map.json", use_env=False)


# Tool constants

Z_HYCAL    = 6270.0    # Legacy ROOT files without run-geometry metadata.
RESOLUTION = 0.03      # 3% / sqrt(E/1000), matching PhysicsTools::findPeakWindow.


def estimate_sigma(E: float) -> float:
    """Expected peak width used to initialize PhysicsTools::fitPeak."""
    return RESOLUTION * math.sqrt(max(E, 0.0) * 1000.0)


# Gaussian fit helpers

def _gauss(x: np.ndarray, amp: float, mu: float, sig: float) -> np.ndarray:
    return amp * np.exp(-0.5 * ((x - mu) / sig) ** 2)


def _gauss_fit_range(
    centers: np.ndarray,
    counts: np.ndarray,
    lo: float,
    hi: float,
    sig0: float,
    sig_max: float,
    mu0: Optional[float] = None,
    maxfev: int = 2000,
) -> Optional[Tuple[float, float, float, float]]:
    """Fit a Gaussian with its mean bounded to [lo, hi] to the non-empty bins
    in that range; ``mu0`` defaults to their weighted mean.  Returns
    (mean, sigma, amp, mean_error), or None with fewer than 4 bins.  curve_fit errors
    propagate."""
    mask = (centers >= lo) & (centers <= hi) & (counts > 0)
    if mask.sum() < 4:
        return None
    xf, yf = centers[mask], counts[mask].astype(float)
    if mu0 is None:
        mu0 = float((xf * yf).sum() / yf.sum())
    popt, pcov = curve_fit(
        _gauss, xf, yf,
        p0=[float(yf.max()), mu0, sig0],
        bounds=([0.0, lo, 1e-6], [np.inf, hi, sig_max]),
        maxfev=maxfev,
    )
    mean_error = math.sqrt(max(float(pcov[1, 1]), 0.0))
    return float(popt[1]), abs(float(popt[2])), float(popt[0]), mean_error


def _find_peak_window(
    counts: np.ndarray,
    edges: np.ndarray,
    Eexp: float,
) -> Optional[Tuple[float, float, float, float, float]]:
    """Mirror PhysicsTools::findPeakWindow for the Gaussian peak fit."""
    n_bins = len(counts)
    if n_bins < 4:
        return None

    first = max(0, int(np.searchsorted(edges, 0.8 * Eexp, side="right") - 1))
    last = min(n_bins - 1, int(np.searchsorted(edges, 1.2 * Eexp, side="right") - 1))
    peak_bin = -1
    peak_height = 0.0
    if first <= last:
        for i in range(first, last + 1):
            left = counts[i - 1] if i > 0 else counts[i]
            right = counts[i + 1] if i < n_bins - 1 else counts[i]
            if counts[i] > 0.0 and counts[i] >= left and counts[i] >= right and counts[i] > peak_height:
                peak_bin, peak_height = i, float(counts[i])

    if peak_bin < 0:
        peak_bin = int(np.argmax(counts))
        peak_height = float(counts[peak_bin])
    if peak_height <= 0.0:
        return None

    centers = 0.5 * (edges[:-1] + edges[1:])
    peak0 = float(centers[peak_bin])
    sigma0 = estimate_sigma(peak0)
    left_bin = right_bin = peak_bin
    threshold = 0.4 * peak_height

    while left_bin > 0 and counts[left_bin] > threshold:
        left_bin -= 1
    while right_bin < n_bins - 1 and counts[right_bin] > threshold:
        right_bin += 1
    while centers[right_bin] - centers[left_bin] < 2.5 * sigma0 and (
        left_bin > 0 or right_bin < n_bins - 1
    ):
        if left_bin > 0:
            left_bin -= 1
        if right_bin < n_bins - 1:
            right_bin += 1
    while centers[right_bin] - centers[left_bin] > 4.5 * sigma0 and (
        left_bin > 0 or right_bin < n_bins - 1
    ):
        left_bin += 1
        right_bin -= 1

    lo, hi = float(centers[left_bin]), float(centers[right_bin])
    if not (hi > lo and math.isfinite(sigma0) and sigma0 > 0.0):
        return None
    return peak_height, peak0, sigma0, lo, hi


def fit_peak_full(
    counts: np.ndarray,
    edges: np.ndarray,
    Eexp: float,
    entries: float,
) -> Tuple[float, float, float, float]:
    """
    Adaptive Gaussian fit matching PhysicsTools::fitPeak(h, Eexp, true).
    Returns (mean_MeV, sigma_MeV, amplitude, mean_error_MeV).
    mean=0 on failure.
    """
    if entries < 100.0 or not HAS_SCIPY:
        return 0.0, 0.0, 0.0, 0.0
    window = _find_peak_window(counts, edges, Eexp)
    if window is None:
        return 0.0, 0.0, 0.0, 0.0

    height, peak0, sigma0, lo, hi = window
    centers = 0.5 * (edges[:-1] + edges[1:])
    mask = (centers >= lo) & (centers <= hi) & (counts > 0.0)
    if np.count_nonzero(mask) < 4:
        return 0.0, 0.0, 0.0, 0.0

    xf = centers[mask]
    yf = counts[mask].astype(float)
    try:
        params, pcov = curve_fit(
            _gauss, xf, yf,
            p0=[height, peak0, sigma0],
            sigma=np.sqrt(yf),
            absolute_sigma=True,
            maxfev=2000,
        )
    except (RuntimeError, ValueError, FloatingPointError):
        return 0.0, 0.0, 0.0, 0.0
    amp, mean, sigma = (float(v) for v in params)
    mean_error = math.sqrt(max(float(pcov[1, 1]), 0.0))
    if not all(math.isfinite(v) for v in (amp, mean, sigma, mean_error)) or abs(sigma) <= 0.0:
        return 0.0, 0.0, 0.0, 0.0
    return mean, abs(sigma), amp, mean_error


# Non-linearity model & fit  (match nonlinearity.cpp)
_NL_FIT_ENERGY_MAX = 3500.0
_NL_FIT_MAX_DEVIATION = 0.2
_NL_MANUAL_COEFFICIENT_LIMIT = 1_000_000.0
_NL_FIXED_VERTEX_ENERGY = 3000.0


def _nl_coefficient_bounds(E_base: float, order: int) -> Tuple[List[float], List[float]]:
    max_abs_t = max(E_base, abs(_NL_FIT_ENERGY_MAX - E_base)) / 1000.0
    deviation = _NL_FIT_MAX_DEVIATION if order == 1 else _NL_FIT_MAX_DEVIATION / 2.0
    nl1_limit = deviation / max_abs_t
    if order == 1:
        return [-nl1_limit], [nl1_limit]
    nl2_limit = deviation / max_abs_t**2
    return [-nl1_limit, -nl2_limit], [nl1_limit, nl2_limit]


def nl_ratio_1st(E_rec: np.ndarray, nl1: float, E_base: float) -> np.ndarray:
    """1st order: E_rec/E_exp = 1 + nl1*(E_rec-E_base)/1000"""
    return 1.0 + nl1 * (E_rec - E_base) / 1000.0


def nl_ratio_2nd(E_rec: np.ndarray, nl1: float, nl2: float, E_base: float) -> np.ndarray:
    """2nd order: E_rec/E_exp = 1 + nl1*t + nl2*t^2, t=(E_rec-E_base)/1000"""
    t = (E_rec - E_base) / 1000.0
    return 1.0 + nl1 * t + nl2 * t * t


def fit_nonlinearity(
    E_exp: np.ndarray,
    E_rec: np.ndarray,
    E_base: float,
    order: int = 1,
    nl1_init: float = 0.01,
    sigma_Erec: Optional[np.ndarray] = None,
    point_weights: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray, float, int]:
    """Fit E_rec/E_exp vs E_rec using propagated peak-mean uncertainties.
    ``point_weights`` multiplies each point's chi-square weight.
    Returns (params, errors, chi2, ndf) with ``order`` params each;
    all zeros when the fit cannot run or fails."""
    failed = (np.zeros(order), np.zeros(order), 0.0, 0)
    if not HAS_SCIPY or len(E_exp) < order + 1 or E_base <= 0.0:
        return failed
    if order == 1:
        model = lambda x, nl1: nl_ratio_1st(x, nl1, E_base)
        bounds = _nl_coefficient_bounds(E_base, order)
        nl1_limit = bounds[1][0]
        p0 = [float(np.clip(nl1_init, -nl1_limit * 0.999, nl1_limit * 0.999))]
    else:
        model = lambda x, nl1, nl2: nl_ratio_2nd(x, nl1, nl2, E_base)
        bounds = _nl_coefficient_bounds(E_base, order)
        nl1_limit = bounds[1][0]
        p0 = [
            float(np.clip(nl1_init, -nl1_limit * 0.999, nl1_limit * 0.999)),
            0.0,
        ]
    safe_E_exp = np.where(E_exp > 0, E_exp, 1.0)
    ratio = E_rec / safe_E_exp
    if sigma_Erec is None:
        sigma_ratio = np.ones_like(ratio)
    else:
        errors = np.asarray(sigma_Erec, dtype=float)
        errors = np.where(np.isfinite(errors) & (errors > 0.0), errors, 1.0)
        sigma_ratio = errors / safe_E_exp
    if point_weights is not None:
        sigma_ratio /= np.sqrt(np.maximum(np.asarray(point_weights, dtype=float), 1.0))
    sigma_ratio = np.maximum(sigma_ratio, 1e-12)
    try:
        popt, pcov = curve_fit(model, E_rec, ratio,
                               p0=p0, bounds=bounds,
                               sigma=sigma_ratio,
                               absolute_sigma=True, maxfev=5000)
        perr = np.sqrt(np.maximum(np.diag(pcov), 0.0))
        chi2 = float(np.sum(((ratio - model(E_rec, *popt)) / sigma_ratio) ** 2))
        return popt, perr, chi2, len(E_rec) - order
    except Exception:
        return failed


# Per-module data

# Beam tag (ROOT dir energy_<tag>GeV, histogram h_energy_<name>_<tag>) and
# beam energy [MeV].
_BEAMS = (("3p5", 3485.41), ("2p2", 2239.51), ("0p7", 728.9))
_BASE_BEAMS = (("3p5", "3.5"), ("2p2", "2.2"), ("0p7", "0.7"))
_KINDS = (("ep", "e-p", "tomato"), ("ee", "e-e", "cornflowerblue"))
_BEAM_TAGS = [tag for tag, _ in _BEAMS]

# Kind-major order (ep_3p5, ep_2p2, ep_0p7, ee_3p5, ...) fills the 2x3 plot grid.
_PEAK_CONFIGS: List[Dict[str, object]] = []
for _kind, _name, _color in _KINDS:
    for _tag, _ebeam in _BEAMS:
        _energy = f"{_tag.replace('p', '.')} GeV"
        _cfg = {
            "key": f"{_kind}_{_tag}", "label": f"{_name}  {_energy}", "ebeam": _ebeam,
            "kind": _kind, "color": _color, "tag": _tag, "energy_label": _energy,
        }
        if _kind == "ee":
            _cfg.update(skip_if_close_to=f"ep_{_tag}",
                        min_separation=170.0 if _tag == "0p7" else 0.0)
        _PEAK_CONFIGS.append(_cfg)
del _kind, _name, _color, _tag, _ebeam, _energy, _cfg
_PEAK_KEYS = [cfg["key"] for cfg in _PEAK_CONFIGS]
_PEAK_BY_KEY = {cfg["key"]: cfg for cfg in _PEAK_CONFIGS}


def peak_cfg(key: str) -> Dict[str, object]:
    return _PEAK_BY_KEY[key]


def hist_for_key(d: "ModuleData", key: str) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    return d.hists[str(peak_cfg(key)["tag"])]


def should_fit_peak(d: "ModuleData", key: str) -> bool:
    cfg = peak_cfg(key)
    other_key = cfg.get("skip_if_close_to")
    min_sep = cfg.get("min_separation")
    if isinstance(other_key, str) and isinstance(min_sep, (float, int)):
        if abs(d.e_exp[other_key] - d.e_exp[key]) <= float(min_sep):
            return False
    return True


class ModuleData:
    """All non-linearity data for one PbWO4 module."""

    def __init__(
        self, name: str, x: float, y: float, row: int = -1, col: int = -1,
        hycal_z_by_beam: Optional[Dict[str, float]] = None,
    ):
        self.name  = name
        self.x     = float(x)
        self.y     = float(y)
        self.row   = int(row)
        self.col   = int(col)
        hycal_z_by_beam = hycal_z_by_beam or {}
        self.theta = math.degrees(math.atan2(
            math.hypot(x, y), hycal_z_by_beam.get("3p5", Z_HYCAL)))

        self.e_exp:     Dict[str, float] = {}
        self.sigma_exp: Dict[str, float] = {}
        for cfg in _PEAK_CONFIGS:
            key = str(cfg["key"])
            z = hycal_z_by_beam.get(str(cfg["tag"]), Z_HYCAL)
            theta = math.degrees(math.atan2(math.hypot(x, y), z))
            E = expected_energy(theta, float(cfg["ebeam"]), str(cfg["kind"]))
            self.e_exp[key]     = E
            self.sigma_exp[key] = estimate_sigma(E) if E > 0.0 else 1.0

        # fitted peak results (0 = not found / invalid)
        self.peaks:      Dict[str, float] = {k: 0.0 for k in _PEAK_KEYS}
        self.peak_sigma: Dict[str, float] = {k: 0.0 for k in _PEAK_KEYS}
        self.peak_mean_error: Dict[str, float] = {k: 0.0 for k in _PEAK_KEYS}
        self.peak_amp:   Dict[str, float] = {k: 0.0 for k in _PEAK_KEYS}
        self.peak_enabled: Dict[str, bool] = {k: True for k in _PEAK_KEYS}
        self.manual_peaks: Dict[str, float] = {k: 0.0 for k in _PEAK_KEYS}
        self.fit_order: str = "2nd"

        # non-linearity fit result (1st order)
        self.nl:      float = 0.0
        self.nl_err:  float = 0.0
        self.chi2_nl: float = 0.0
        self.ndf_nl:  int   = 0
        # non-linearity fit result (2nd order)
        self.nl2_1:      float = 0.0
        self.nl2_1_err:  float = 0.0
        self.nl2_2:      float = 0.0
        self.nl2_2_err:  float = 0.0
        self.chi2_nl2:   float = 0.0
        self.ndf_nl2:    int   = 0
        self.nl_coefficients_manual: bool = False
        self.nl2_vertex_fixed: bool = False
        self.mean_abs_data_fit_deviation_2nd = math.nan
        self.base_ee_corr_residual_2nd = math.nan

        # geometry flags (0-based row/col)
        self.is_absorber = self.row >= 0 and hole_ring(self.row + 1, self.col + 1) <= 1

        # raw histograms loaded from ROOT (counts, edges), keyed by beam tag
        self.hists: Dict[str, Optional[Tuple[np.ndarray, np.ndarray]]] = {
            t: None for t in _BEAM_TAGS}
        self.hist_entries: Dict[str, float] = {t: 0.0 for t in _BEAM_TAGS}

    def measured_peak(self, key: str) -> float:
        """Manual peak center when set, else the fitted peak (0 = none)."""
        if self.manual_peaks.get(key, 0.0) > 0.0:
            return self.manual_peaks[key]
        return self.peaks[key]

    def is_outer_layer(self, layers: int) -> bool:
        return (layers > 0 and self.row >= 0 and self.col >= 0
                and edge_depth(self.row + 1, self.col + 1) <= layers)

    def linearity_points(self, base_beam_tag: str = "2p2") -> List[Tuple[str, float, float]]:
        """(key, E_expected, E_rec) for every enabled peak with a valid measurement."""
        if self.measured_peak(f"ep_{base_beam_tag}") <= 0.0:
            return []
        pts = []
        for key in _PEAK_KEYS:
            p = self.measured_peak(key)
            if self.peak_enabled.get(key, True) and p > 0.0:
                pts.append((key, self.e_exp[key], p))
        return pts

    def auto_fit_peaks(self) -> bool:
        """Adaptive PhysicsTools Gaussian fit of every peak; True if any succeeded."""
        any_ok = False
        for key in _PEAK_KEYS:
            hist = hist_for_key(self, key)
            E = self.e_exp[key]
            if hist is None or E <= 0.0 or not should_fit_peak(self, key):
                continue
            entries = self.hist_entries[str(peak_cfg(key)["tag"])]
            mu, sigma_f, amp, mean_error = fit_peak_full(hist[0], hist[1], E, entries)
            if mu > 0.0:
                self.peaks[key], self.peak_sigma[key], self.peak_amp[key] = mu, sigma_f, amp
                self.peak_mean_error[key] = mean_error
                any_ok = True
        return any_ok

    def refit_nl(self, outer_zero_layers: int = 1, base_beam_tag: str = "2p2") -> None:
        """Re-run the 1st and 2nd order non-linearity fits from current measured peaks."""
        pts = self.linearity_points(base_beam_tag)
        E_exp = np.array([p[1] for p in pts])
        E_rec = np.array([p[2] for p in pts])
        # A manually entered peak without covariance gets a conservative 1 MeV error.
        sigma_Erec = np.array([
            error if math.isfinite(error) and error > 0.0 else 1.0
            for error in (self.peak_mean_error[p[0]] for p in pts)
        ])
        point_weights = np.array([
            2.0 if peak_cfg(p[0])["tag"] == base_beam_tag else 1.0
            for p in pts
        ])
        E_base = self.e_exp[f"ep_{base_beam_tag}"]
        force_zero = self.is_outer_layer(outer_zero_layers) or self.is_absorber
        (nl,), (err,), self.chi2_nl, self.ndf_nl = fit_nonlinearity(
            E_exp, E_rec, E_base, sigma_Erec=sigma_Erec,
            point_weights=point_weights)
        self.nl     = 0.0 if force_zero else float(nl)
        self.nl_err = float(err)
        (nl1, nl2), (nl1_err, nl2_err), self.chi2_nl2, self.ndf_nl2 = fit_nonlinearity(
            E_exp, E_rec, E_base, order=2, nl1_init=self.nl,
            sigma_Erec=sigma_Erec, point_weights=point_weights)
        self.nl2_1     = 0.0 if force_zero else float(nl1)
        self.nl2_1_err = float(nl1_err)
        self.nl2_2     = 0.0 if force_zero else float(nl2)
        self.nl2_2_err = float(nl2_err)
        self.nl_coefficients_manual = False
        self.nl2_vertex_fixed = False

        self.update_correction_residuals(
            E_base, E_exp, E_rec, [p[0] for p in pts], base_beam_tag)

    def fit_nl2_fixed_vertex(self, base_beam_tag: str) -> Tuple[float, float]:
        """Fit base-beam e-p/e-e points with R(E_base)=1 and a fixed vertex."""
        keys = (f"ep_{base_beam_tag}", f"ee_{base_beam_tag}")
        E_base = self.e_exp[keys[0]]
        if E_base <= 0.0:
            raise ValueError("Base-energy e-p expected energy is unavailable.")

        vertex_t = (_NL_FIXED_VERTEX_ENERGY - E_base) / 1000.0
        basis = []
        residuals = []
        weights = []
        for key in keys:
            E_exp = self.e_exp[key]
            E_rec = self.measured_peak(key)
            if E_exp <= 0.0 or E_rec <= 0.0:
                raise ValueError(
                    f"{peak_cfg(key)['label']} peak is unavailable; fit both base-energy peaks first.")
            t = (E_rec - E_base) / 1000.0
            basis.append(t * (t - 2.0 * vertex_t))
            residuals.append(E_rec / E_exp - 1.0)
            error = self.peak_mean_error[key]
            sigma_Erec = error if math.isfinite(error) and error > 0.0 else 1.0
            sigma_ratio = sigma_Erec / E_exp
            weights.append(2.0 / max(sigma_ratio * sigma_ratio, 1e-24))

        basis_arr = np.asarray(basis)
        residuals_arr = np.asarray(residuals)
        weights_arr = np.asarray(weights)
        denominator = float(np.sum(weights_arr * basis_arr * basis_arr))
        if not math.isfinite(denominator) or denominator <= 0.0:
            raise ValueError("Base-energy peaks do not constrain the fixed-vertex fit.")

        nl2 = float(np.sum(weights_arr * basis_arr * residuals_arr) / denominator)
        if not math.isfinite(nl2) or nl2 >= 0.0:
            raise ValueError(
                "The base-energy peaks do not produce a concave-down fit with a maximum.")
        nl1 = -2.0 * nl2 * vertex_t
        if max(abs(nl1), abs(nl2)) > _NL_MANUAL_COEFFICIENT_LIMIT:
            raise ValueError("The fitted coefficients exceed the manual coefficient limit.")
        return nl1, nl2

    def update_correction_residuals(
        self, E_base: float, E_exp: np.ndarray, E_rec: np.ndarray,
        keys: List[str], base_beam_tag: str,
    ) -> None:
        # Mean |corrected ratio - 1| for map colouring
        self.mean_corr_resid_1st: float = 0.0
        self.mean_corr_resid_2nd: float = 0.0
        self.mean_abs_data_fit_deviation_2nd = math.nan
        self.base_ee_corr_residual_2nd = math.nan
        if len(E_exp) > 0 and E_base > 0.0:
            resids1, resids2 = [], []
            fit_deviations_2nd = []
            base_ee_residual = None
            has_2nd_fit = self.ndf_nl2 > 0 or self.nl_coefficients_manual
            for key, Ee, Er in zip(keys, E_exp, E_rec):
                if Ee <= 0.0:
                    continue
                d1 = nl_ratio_1st(Er, self.nl, E_base)
                if d1 != 0.0:
                    Ec1 = Er / d1
                    resids1.append(abs(Ec1 / Ee - 1.0))
                d2 = nl_ratio_2nd(Er, self.nl2_1, self.nl2_2, E_base)
                if d2 != 0.0:
                    Ec2 = Er / d2
                    resids2.append(abs(Ec2 / Ee - 1.0))
                    if has_2nd_fit:
                        fit_deviations_2nd.append(
                            abs(Er / Ee - d2))
                        if key == f"ee_{base_beam_tag}":
                            base_ee_residual = Ec2 / Ee - 1.0
            if resids1:
                self.mean_corr_resid_1st = float(np.mean(resids1))
            if resids2:
                self.mean_corr_resid_2nd = float(np.mean(resids2))
            if fit_deviations_2nd:
                self.mean_abs_data_fit_deviation_2nd = float(
                    np.mean(fit_deviations_2nd))
            if base_ee_residual is not None:
                self.base_ee_corr_residual_2nd = float(base_ee_residual)


# Data loading

def _load_json_geometry(path: Path) -> Dict[str, Tuple[float, float, int, int]]:
    """Read hycal_map.json and return {name: (x, y, row_0based, col_0based)}
    for in-fiducial PbWO4 modules, matching nonlinearity.cpp."""
    result: Dict[str, Tuple[float, float, int, int]] = {}
    try:
        result = {m.name: (float(m.x), float(m.y), m.row - 1, m.col - 1)
                  for m in load_modules(path)
                  if m.mod_type == "PbWO4"
                  and (abs(m.x) > 2.0 * 20.77 or abs(m.y) > 2.0 * 20.75)
                  and abs(m.x) < 16.0 * 20.77
                  and abs(m.y) < 16.0 * 20.75}
    except Exception as e:
        print(f"Warning: cannot read geometry from {path}: {e}", file=sys.stderr)
    return result


def _read_hycal_z(root_dir) -> Optional[float]:
    """Read the per-beam HyCal Z metadata written by nonlinearity.cpp."""
    if root_dir is None:
        return None
    try:
        value = float(root_dir["hycal_z_mm"].member("fTitle"))
    except (KeyError, AttributeError, TypeError, ValueError):
        return None
    return value if math.isfinite(value) and value > 0.0 else None


def _read_base_beam_tag(root_path: Path) -> Optional[str]:
    """Read the selected base beam written by nonlinearity.cpp."""
    if not HAS_UPROOT:
        return None
    try:
        with uproot.open(str(root_path)) as root_file:
            label = str(root_file["base_beam_energy_GeV"].member("fTitle")).strip()
    except (OSError, KeyError, AttributeError, TypeError, ValueError):
        return None
    return next((tag for tag, energy in _BASE_BEAMS if label == energy), None)


def load_all(root_path: Optional[Path]) -> Dict[str, ModuleData]:
    """Load module geometry, histograms from ROOT, and fit all peaks."""
    if MODULES_JSON is None:
        print("Warning: hycal_map.json not found", file=sys.stderr)
        return {}

    hycal_z_by_beam = {}
    root_dirs = {}
    if root_path is not None and root_path.is_file() and HAS_UPROOT:
        try:
            with uproot.open(str(root_path)) as f:
                root_dirs = {t: f.get(f"energy_{t}GeV") for t in _BEAM_TAGS}
                for tag, root_dir in root_dirs.items():
                    z = _read_hycal_z(root_dir)
                    if z is not None:
                        hycal_z_by_beam[tag] = z
            if not hycal_z_by_beam:
                print(
                    "Warning: ROOT file has no hycal_z_mm metadata; "
                    f"using legacy {Z_HYCAL:g} mm geometry.",
                    file=sys.stderr,
                )
        except Exception as e:
            print(f"ROOT read error: {e}", file=sys.stderr)

    geo = _load_json_geometry(MODULES_JSON)
    data: Dict[str, ModuleData] = {}
    for name, (x, y, row, col) in geo.items():
        data[name] = ModuleData(name, x, y, row, col, hycal_z_by_beam)

    if root_path is not None and root_path.is_file() and HAS_UPROOT:
        try:
            with uproot.open(str(root_path)) as f:
                root_dirs = {t: f.get(f"energy_{t}GeV") for t in _BEAM_TAGS}
                for name, d in data.items():
                    for t, root_dir in root_dirs.items():
                        if root_dir is None:
                            continue
                        try:
                            hist = root_dir[f"h_energy_{name}_{t}"]
                            d.hists[t] = hist.to_numpy()
                            d.hist_entries[t] = float(hist.member("fEntries"))
                        except Exception:
                            pass
        except Exception as e:
            print(f"ROOT read error: {e}", file=sys.stderr)
    elif root_path is not None and not HAS_UPROOT:
        print("Warning: uproot not installed — cannot read ROOT file.", file=sys.stderr)

    for d in data.values():
        d.auto_fit_peaks()

    return data


# Matplotlib canvas

class PlotCanvas(FigureCanvas):
    def __init__(self):
        fc = THEME.CANVAS
        self.fig = Figure(facecolor=fc)
        self.ax = self.fig.add_subplot(111)
        self.ax.set_facecolor(fc)
        super().__init__(self.fig)
        self.setSizePolicy(QSizePolicy.Policy.Expanding, QSizePolicy.Policy.Expanding)

    def _style(self, ax, title="", xlabel="", ylabel=""):
        tc = THEME.TEXT
        dc = THEME.TEXT_DIM
        ax.set_title(title, color=tc, fontsize=9, pad=3)
        ax.set_xlabel(xlabel, color=dc, fontsize=8)
        ax.set_ylabel(ylabel, color=dc, fontsize=8)
        ax.tick_params(colors=dc, labelsize=7)
        for sp in ax.spines.values():
            sp.set_color(dc)

    def redraw(self):
        self.fig.canvas.draw_idle()


# Main window

class NLViewerWindow(QMainWindow):

    def __init__(self, root_path: Optional[Path] = None):
        super().__init__()
        self.setWindowTitle("HyCal Non-linearity Viewer")
        self.resize(1480, 940)

        self._data: Dict[str, ModuleData] = {}
        self._cur:  Optional[ModuleData]  = None
        self._w_modules: List[Module] = []   # W-only Module list for the map

        # SpanSelector state: per peak key → (xmin, xmax) or None
        self._span_range: Dict[str, Optional[Tuple[float, float]]] = {
            k: None for k in _PEAK_KEYS}
        self._span_selectors: List = []   # keep references to avoid GC

        self._peak_use_cb: Dict[str, QCheckBox] = {}
        self._peak_center_edit: Dict[str, QLineEdit] = {}
        self._manual_nl_spins: Dict[str, QDoubleSpinBox] = {}
        self._hist_canvases: Dict[str, PlotCanvas] = {}
        self._linearity_canvas: Optional[PlotCanvas] = None
        self._fit_order_cb: Optional[QComboBox] = None
        self._base_energy_cb: Optional[QComboBox] = None
        self._base_beam_tag = "2p2"

        # scan state
        self._scan_timer  = QTimer(self)
        self._scan_timer.setInterval(500)   # ms per module
        self._scan_timer.timeout.connect(self._scan_step)
        self._scan_index  = 0
        self._scan_names: List[str] = []

        # track active map mode to detect changes and preserve user range
        self._current_map_mode: str = ""

        self._build_ui()
        if root_path is not None:
            self._load(root_path)

    # ── UI construction ────────────────────────────────────────────────────

    def _build_ui(self):
        central = QWidget()
        self.setCentralWidget(central)
        root_lay = QVBoxLayout(central)
        root_lay.setContentsMargins(4, 4, 4, 4)
        root_lay.setSpacing(4)

        top = QHBoxLayout()
        btn = QPushButton("Open ROOT file…")
        btn.clicked.connect(self._open_file)
        top.addWidget(btn)
        self._file_lbl = QLabel("No file loaded")
        self._file_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM};")
        top.addWidget(self._file_lbl, stretch=1)
        top.addWidget(QLabel("Base energy:"))
        self._base_energy_cb = QComboBox()
        for tag, energy in _BASE_BEAMS:
            self._base_energy_cb.addItem(f"{energy} GeV", tag)
        self._base_energy_cb.setCurrentIndex(1)
        self._base_energy_cb.setToolTip(
            "Beam energy whose e-p expected energy defines E_base for both nl fits.")
        self._base_energy_cb.currentIndexChanged.connect(self._on_base_energy_changed)
        top.addWidget(self._base_energy_cb)
        top.addWidget(QLabel("Color map:"))
        self._mode_cb = QComboBox()
        self._mode_cb.addItems([
            "Non-linearity (nl)",
            "Mean |corr/exp-1| (1st order)",
            "Mean |corr/exp-1| (2nd order)",
            "Mean |data - 2nd-order fit| (E_rec/E_exp)",
            "Base-energy e-e (corr/exp - 1), 2nd order",
            "Flat",
        ])
        self._mode_cb.currentIndexChanged.connect(self._on_mode_changed)
        top.addWidget(self._mode_cb)

        top.addSpacing(12)
        top.addWidget(QLabel("Zero outer W layers:"))
        self._outer_layer_spin = QSpinBox()
        self._outer_layer_spin.setRange(0, 17)
        self._outer_layer_spin.setValue(1)
        self._outer_layer_spin.setFixedWidth(54)
        self._outer_layer_spin.setToolTip(
            "Set non-linearity coefficients to 0 for this many outer W-module layers.\n"
            "0 keeps all non-absorber W modules fitted; 1 matches the historical behavior.")
        self._outer_layer_spin.valueChanged.connect(self._on_outer_layers_changed)
        top.addWidget(self._outer_layer_spin)

        top.addSpacing(12)
        self._scan_btn = QPushButton("▶ Scan")
        self._scan_btn.setToolTip(
            "Auto-scan through all W modules.\n"
            "Press Enter or click the button again to stop.")
        self._scan_btn.setEnabled(False)
        self._scan_btn.clicked.connect(self._toggle_scan)
        top.addWidget(self._scan_btn)

        top.addWidget(QLabel("Speed:"))
        self._scan_speed_cb = QComboBox()
        self._scan_speed_cb.setFixedWidth(90)
        for label, ms in [("0.2 s", 200), ("0.5 s", 500), ("1 s", 1000),
                          ("2 s", 2000), ("5 s", 5000)]:
            self._scan_speed_cb.addItem(label, ms)
        self._scan_speed_cb.setCurrentIndex(1)   # default 0.5 s
        self._scan_speed_cb.currentIndexChanged.connect(self._on_scan_speed_changed)
        top.addWidget(self._scan_speed_cb)

        self._scan_lbl = QLabel("")
        self._scan_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM}; font-size: 11px;")
        top.addWidget(self._scan_lbl)

        top.addSpacing(12)
        self._export_btn = QPushButton("↩ Export nl → JSON")
        self._export_btn.setToolTip(
            "Read a base calibration JSON and save nl1/nl2 coefficients from the current\n"
            "fit or manual adjustment (2nd order by default) and selected base beam energy.")
        self._export_btn.setEnabled(False)
        self._export_btn.clicked.connect(self._export_nl)
        top.addWidget(self._export_btn)

        root_lay.addLayout(top)

        stop_sc = QShortcut(QKeySequence(Qt.Key.Key_Return), self)
        stop_sc.activated.connect(self._stop_scan)
        stop_sc2 = QShortcut(QKeySequence(Qt.Key.Key_Enter), self)
        stop_sc2.activated.connect(self._stop_scan)
        prev_module_sc = QShortcut(QKeySequence(Qt.Key.Key_PageUp), self)
        prev_module_sc.activated.connect(lambda: self._shortcut_jump_module(-1))
        next_module_sc = QShortcut(QKeySequence(Qt.Key.Key_PageDown), self)
        next_module_sc.activated.connect(lambda: self._shortcut_jump_module(1))

        # horizontal splitter: map | plots+controls
        splitter = QSplitter(Qt.Orientation.Horizontal)
        root_lay.addWidget(splitter, stretch=1)

        # ── left: map ─────────────────────────────────────────────────────
        map_wrap = QWidget()
        map_lay  = QVBoxLayout(map_wrap)
        map_lay.setContentsMargins(0, 0, 0, 0)
        self._map = HyCalMapWidget(map_wrap, enable_zoom_pan=True, min_size=(440, 440))
        self._map.moduleClicked.connect(self._on_click)
        map_lay.addWidget(self._map)
        self._range_ctrl = ColorRangeControl(self._map, orientation="horizontal",
                                             parent=map_wrap)
        map_lay.addWidget(self._range_ctrl)
        splitter.addWidget(map_wrap)

        # ── right: plots + controls ────────────────────────────────────────
        right     = QWidget()
        right_lay = QVBoxLayout(right)
        right_lay.setContentsMargins(4, 0, 4, 4)
        right_lay.setSpacing(4)

        module_nav = QHBoxLayout()
        self._mod_lbl = QLabel("← click a W module on the map")
        self._mod_lbl.setFont(QFont("Monospace", 12, QFont.Weight.Bold))
        self._mod_lbl.setStyleSheet(f"color: {THEME.ACCENT};")
        module_nav.addWidget(self._mod_lbl, stretch=1)
        self._prev_module_btn = QPushButton("← Previous module")
        self._prev_module_btn.setEnabled(False)
        self._prev_module_btn.setToolTip("Select the previous module by W number (PageUp).")
        self._prev_module_btn.clicked.connect(lambda: self._jump_module(-1))
        module_nav.addWidget(self._prev_module_btn)
        self._next_module_btn = QPushButton("Next module →")
        self._next_module_btn.setEnabled(False)
        self._next_module_btn.setToolTip("Select the next module by W number (PageDown).")
        self._next_module_btn.clicked.connect(lambda: self._jump_module(1))
        module_nav.addWidget(self._next_module_btn)
        right_lay.addLayout(module_nav)

        right_lay.addWidget(self._build_refit_panel())
        right_lay.addWidget(self._build_histogram_grid(), stretch=3)
        self._build_manual_nl_controls(right_lay)
        self._linearity_canvas = PlotCanvas()
        right_lay.addWidget(self._linearity_canvas, stretch=2)

        splitter.addWidget(right)
        splitter.setSizes([460, 900])

    def _build_manual_nl_controls(self, parent_layout: QVBoxLayout) -> None:
        group = QGroupBox("Manual non-linearity coefficients")
        layout = QVBoxLayout(group)
        layout.setContentsMargins(8, 3, 8, 6)
        layout.setSpacing(3)
        controls = QHBoxLayout()
        controls.setSpacing(6)
        controls.addWidget(QLabel("1st nl₁:"))
        controls.addWidget(self._new_nl_spin("nl1"))
        controls.addSpacing(8)
        controls.addWidget(QLabel("2nd nl₁:"))
        controls.addWidget(self._new_nl_spin("nl2_1"))
        controls.addWidget(QLabel("nl₂:"))
        controls.addWidget(self._new_nl_spin("nl2_2"))
        controls.addStretch()
        layout.addLayout(controls)

        actions = QHBoxLayout()
        self._fit_fixed_vertex_btn = QPushButton("Fit 2nd (vertex = 3000 MeV)")
        self._fit_fixed_vertex_btn.setEnabled(False)
        self._fit_fixed_vertex_btn.setToolTip(
            "Fit only the selected base-energy e-p and e-e peaks.\n"
            "Constrain the 2nd-order curve to R(E_base)=1 and its maximum to 3000 MeV.")
        self._fit_fixed_vertex_btn.clicked.connect(self._do_fit_fixed_vertex)
        hint = QLabel("Edit or use arrows; curves and corrected points update live.")
        hint.setStyleSheet(f"color: {THEME.TEXT_DIM}; font-size: 10px;")
        actions.addWidget(self._fit_fixed_vertex_btn)
        actions.addWidget(hint)
        actions.addStretch()
        layout.addLayout(actions)
        parent_layout.addWidget(group)

    def _build_histogram_grid(self) -> QWidget:
        grid_widget = QWidget()
        grid = QGridLayout(grid_widget)
        grid.setContentsMargins(0, 0, 0, 0)
        grid.setHorizontalSpacing(8)
        grid.setVerticalSpacing(4)
        for idx, key in enumerate(_PEAK_KEYS):
            cell = QWidget()
            cell_lay = QVBoxLayout(cell)
            cell_lay.setContentsMargins(0, 0, 0, 0)
            cell_lay.setSpacing(0)
            canvas = PlotCanvas()
            canvas.setMinimumHeight(100)
            self._hist_canvases[key] = canvas
            cell_lay.addWidget(canvas, stretch=1)
            cell_lay.addWidget(self._build_peak_control_cell(key))
            grid.addWidget(cell, idx // 3, idx % 3)
        return grid_widget

    def _build_peak_control_cell(self, key: str) -> QWidget:
        cell = QWidget()
        cell_lay = QHBoxLayout(cell)
        cell_lay.setContentsMargins(2, 0, 2, 0)
        cell_lay.setSpacing(4)

        use_cb = QCheckBox("Use")
        use_cb.setChecked(True)
        use_cb.setEnabled(False)
        use_cb.setToolTip("Include this peak point in the linearity fit.")
        use_cb.toggled.connect(
            lambda checked, _key=key: self._on_peak_use_changed(_key, checked))

        center_edit = QLineEdit()
        center_edit.setPlaceholderText("μ (MeV)")
        center_edit.setFixedWidth(78)
        center_edit.setEnabled(False)
        center_edit.setToolTip(
            "Manual peak center in MeV. Leave empty to use the fitted peak.")
        center_edit.editingFinished.connect(
            lambda _key=key: self._on_peak_center_edit(_key))

        cell_lay.addStretch()
        cell_lay.addWidget(use_cb)
        cell_lay.addWidget(center_edit)
        cell_lay.addStretch()

        self._peak_use_cb[key] = use_cb
        self._peak_center_edit[key] = center_edit
        return cell

    def _new_nl_spin(self, coefficient: str) -> QDoubleSpinBox:
        spin = QDoubleSpinBox()
        spin.setDecimals(6)
        spin.setRange(-_NL_MANUAL_COEFFICIENT_LIMIT, _NL_MANUAL_COEFFICIENT_LIMIT)
        spin.setSingleStep(0.001)
        spin.setKeyboardTracking(True)
        spin.setFixedWidth(105)
        spin.setEnabled(False)
        spin.setToolTip(
            "Manual coefficients can exceed the automatic fit bounds. "
            "Large values may move the curve outside the visible plot range.")
        spin.valueChanged.connect(
            lambda value, name=coefficient: self._on_manual_nl_changed(name, value))
        self._manual_nl_spins[coefficient] = spin
        return spin

    def _build_refit_panel(self) -> QGroupBox:
        box = QGroupBox("Peak re-fit  (drag on histogram to auto-fit)")
        lay = QVBoxLayout(box)
        lay.setContentsMargins(8, 10, 8, 6)
        lay.setSpacing(4)

        row1 = QHBoxLayout()

        hint_lbl = QLabel("Drag on any peak histogram to auto-fit that peak.")
        hint_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM}; font-size: 11px;")
        row1.addWidget(hint_lbl)
        row1.addSpacing(12)
        row1.addWidget(QLabel("Range:"))
        self._range_lbl = QLabel("(drag to set → auto-fit)")
        self._range_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM}; font-size: 11px;")
        row1.addWidget(self._range_lbl)

        self._clear_range_btn = QPushButton("✕ Clear")
        self._clear_range_btn.setFixedWidth(60)
        self._clear_range_btn.setEnabled(False)
        self._clear_range_btn.setToolTip(
            "Clear the drag-selected fit range and revert to the adaptive PhysicsTools fit window.")
        self._clear_range_btn.clicked.connect(self._clear_span_range)
        row1.addWidget(self._clear_range_btn)

        row1.addSpacing(16)

        self._refit_peak_btn = QPushButton("Refit All Peaks")
        self._refit_peak_btn.setEnabled(False)
        self._refit_peak_btn.setToolTip(
            "Re-run the adaptive PhysicsTools Gaussian fit on all 6 peaks (clears manual drag ranges).\n"
            "Updates the linearity plot.")
        self._refit_peak_btn.clicked.connect(self._do_refit_all_peaks)
        row1.addWidget(self._refit_peak_btn)

        row1.addSpacing(16)

        self._refit_nl_btn = QPushButton("Refit Linearity")
        self._refit_nl_btn.setEnabled(False)
        self._refit_nl_btn.setToolTip(
            "Re-run the non-linearity fit using current measurement points.\n"
            "Weighted by each peak-fit mean uncertainty; peaks at the selected base beam have 2x weight.\n"
            "1st: E_rec/E_exp = 1 + nl1*(E_rec-E_base)/1000\n"
            "2nd: + nl2*((E_rec-E_base)/1000)^2")
        self._refit_nl_btn.clicked.connect(self._do_refit_nl)
        row1.addWidget(self._refit_nl_btn)

        row1.addSpacing(12)
        row1.addWidget(QLabel("Export fit:"))
        self._fit_order_cb = QComboBox()
        self._fit_order_cb.addItems(["2nd order", "1st order"])
        self._fit_order_cb.setCurrentIndex(0)
        self._fit_order_cb.setEnabled(False)
        self._fit_order_cb.setToolTip(
            "Coefficient order used when exporting JSON for this module.\n"
            "2nd order is the default and matches the C++ calibration output.\n"
            "1st order writes the linear coefficient and nl2 = 0.")
        self._fit_order_cb.currentIndexChanged.connect(self._on_fit_order_changed)
        row1.addWidget(self._fit_order_cb)

        row1.addStretch()
        lay.addLayout(row1)

        self._status_lbl = QLabel("(select a module first)")
        self._status_lbl.setFont(QFont("Monospace", 9))
        self._status_lbl.setWordWrap(True)
        self._status_lbl.setStyleSheet(f"color: {THEME.TEXT_DIM};")
        lay.addWidget(self._status_lbl)

        return box

    # ── file loading ───────────────────────────────────────────────────────

    def _open_file(self):
        path, _ = QFileDialog.getOpenFileName(
            self, "Open nonlinearity ROOT file", "",
            "ROOT files (*.root);;All files (*)")
        if path:
            self._load(Path(path))

    def _after_load(self):
        self._scan_btn.setEnabled(bool(self._data))
        self._export_btn.setEnabled(bool(self._data))
        self._prev_module_btn.setEnabled(bool(self._data))
        self._next_module_btn.setEnabled(bool(self._data))

    def _set_peak_controls_enabled(self, enabled: bool):
        for key in _PEAK_KEYS:
            self._peak_use_cb[key].setEnabled(enabled)
            self._peak_center_edit[key].setEnabled(enabled)
        if self._fit_order_cb is not None:
            self._fit_order_cb.setEnabled(enabled)

    def _sync_peak_controls(self, d: Optional[ModuleData]):
        enabled = d is not None
        self._set_peak_controls_enabled(enabled)
        for key in _PEAK_KEYS:
            cb = self._peak_use_cb[key]
            edit = self._peak_center_edit[key]
            cb.blockSignals(True)
            edit.blockSignals(True)
            if d is None:
                cb.setChecked(True)
                edit.setText("")
            else:
                cb.setChecked(d.peak_enabled.get(key, True))
                manual = d.manual_peaks.get(key, 0.0)
                edit.setText(f"{manual:.2f}" if manual > 0.0 else "")
            cb.blockSignals(False)
            edit.blockSignals(False)
        if self._fit_order_cb is not None:
            self._fit_order_cb.blockSignals(True)
            if d is None:
                self._fit_order_cb.setCurrentIndex(0)
            else:
                self._fit_order_cb.setCurrentIndex(0 if d.fit_order == "2nd" else 1)
            self._fit_order_cb.blockSignals(False)
        self._sync_manual_nl_controls(d)

    def _sync_manual_nl_controls(self, d: Optional[ModuleData]) -> None:
        if d is None:
            for spin in self._manual_nl_spins.values():
                spin.setEnabled(False)
            self._fit_fixed_vertex_btn.setEnabled(False)
            return
        values = {"nl1": d.nl, "nl2_1": d.nl2_1, "nl2_2": d.nl2_2}
        for key, spin in self._manual_nl_spins.items():
            spin.blockSignals(True)
            spin.setRange(-_NL_MANUAL_COEFFICIENT_LIMIT, _NL_MANUAL_COEFFICIENT_LIMIT)
            spin.setValue(float(np.clip(
                values[key], -_NL_MANUAL_COEFFICIENT_LIMIT,
                _NL_MANUAL_COEFFICIENT_LIMIT)))
            spin.setEnabled(True)
            spin.blockSignals(False)
        self._fit_fixed_vertex_btn.setEnabled(True)

    def _on_manual_nl_changed(self, coefficient: str, value: float) -> None:
        d = self._cur
        if d is None:
            return
        if coefficient == "nl1":
            d.nl = float(value)
        elif coefficient == "nl2_1":
            d.nl2_1 = float(value)
            d.nl2_vertex_fixed = False
        else:
            d.nl2_2 = float(value)
            d.nl2_vertex_fixed = False
        d.nl_coefficients_manual = True
        points = d.linearity_points(self._base_beam_tag)
        E_exp = np.array([p[1] for p in points])
        E_rec = np.array([p[2] for p in points])
        E_base = d.e_exp[f"ep_{self._base_beam_tag}"]
        d.update_correction_residuals(
            E_base, E_exp, E_rec, [p[0] for p in points], self._base_beam_tag)
        self._refresh_map()
        self._refresh_current_fit_display(refit=False)
        self._status_lbl.setText(
            f"{d.name}: manual coefficients updated; click Refit Linearity to restore fitted values.")

    def _do_fit_fixed_vertex(self) -> None:
        d = self._cur
        if d is None:
            return
        try:
            nl1, nl2 = d.fit_nl2_fixed_vertex(self._base_beam_tag)
        except ValueError as exc:
            self._status_lbl.setText(f"Fixed-vertex fit failed: {exc}")
            return

        d.nl2_1, d.nl2_2 = nl1, nl2
        d.nl_coefficients_manual = True
        d.nl2_vertex_fixed = True
        points = d.linearity_points(self._base_beam_tag)
        E_base = d.e_exp[f"ep_{self._base_beam_tag}"]
        d.update_correction_residuals(
            E_base,
            np.array([p[1] for p in points]),
            np.array([p[2] for p in points]),
            [p[0] for p in points],
            self._base_beam_tag,
        )
        self._refresh_map()
        self._refresh_current_fit_display(
            status=(
                f"{d.name}: fitted base-energy e-p/e-e only; "
                f"2nd-order maximum fixed at {_NL_FIXED_VERTEX_ENERGY:.0f} MeV "
                f"(nl₁={nl1:.6f}, nl₂={nl2:.6f}); R(E_base)=1."),
            refit=False,
        )

    def _refresh_current_fit_display(self, status: Optional[str] = None, refit: bool = True):
        d = self._cur
        if d is None:
            return
        if refit:
            d.refit_nl(self._outer_zero_layers(), self._base_beam_tag)
            self._refresh_map()
        self._sync_manual_nl_controls(d)
        self._draw_module(d)
        if d.nl_coefficients_manual:
            self._mod_lbl.setText(
                f"Module: {d.name}   θ={d.theta:.3f}°   "
                f"manual 1st nl₁={d.nl:.6f}; "
                f"2nd (nl₁, nl₂)=({d.nl2_1:.6f}, {d.nl2_2:.6f})")
        else:
            self._mod_lbl.setText(
                f"Module: {d.name}   θ={d.theta:.3f}°   "
                f"nl = {d.nl:.4f} ± {d.nl_err:.4f}   "
                f"χ²/ndf = {d.chi2_nl:.2f}/{d.ndf_nl}")
        if status is not None:
            self._status_lbl.setText(status)

    def _on_peak_use_changed(self, key: str, checked: bool):
        d = self._cur
        if d is None:
            return
        d.peak_enabled[key] = bool(checked)
        self._refresh_current_fit_display(
            f"{peak_cfg(key)['label']} {'included in' if checked else 'removed from'} linearity fit.")

    def _on_peak_center_edit(self, key: str):
        d = self._cur
        if d is None:
            return
        text = self._peak_center_edit[key].text().strip()
        if not text:
            d.manual_peaks[key] = 0.0
            self._refresh_current_fit_display(
                f"{peak_cfg(key)['label']} manual center cleared; using fitted peak.")
            return
        try:
            value = float(text)
        except ValueError:
            self._peak_center_edit[key].setText(
                f"{d.manual_peaks[key]:.2f}" if d.manual_peaks.get(key, 0.0) > 0.0 else "")
            self._status_lbl.setText("Manual peak center must be a number in MeV.")
            return
        if value <= 0.0:
            d.manual_peaks[key] = 0.0
            self._peak_center_edit[key].setText("")
            self._refresh_current_fit_display(
                f"{peak_cfg(key)['label']} manual center cleared; using fitted peak.")
            return
        d.manual_peaks[key] = value
        self._peak_center_edit[key].setText(f"{value:.2f}")
        self._refresh_current_fit_display(
            f"{peak_cfg(key)['label']} center manually set to {value:.2f} MeV.")

    def _on_fit_order_changed(self, _idx: int):
        d = self._cur
        if d is None or self._fit_order_cb is None:
            return
        d.fit_order = "2nd" if self._fit_order_cb.currentIndex() == 0 else "1st"
        self._status_lbl.setText(
            f"{d.name} export will use {self._fit_order_cb.currentText()} coefficients.")

    def _on_base_energy_changed(self, _idx: int):
        if self._base_energy_cb is None:
            return
        self._base_beam_tag = str(self._base_energy_cb.currentData())
        if not self._data:
            return
        self._refit_all_nl()
        self._refresh_map()
        self._refresh_current_fit_display(refit=False)
        energy = dict(_BASE_BEAMS)[self._base_beam_tag]
        self._status_lbl.setText(f"Base beam energy set to {energy} GeV.")

    def _load(self, path: Path):
        self._file_lbl.setText(f"Loading {path.name} …")
        QApplication.processEvents()
        root_base = _read_base_beam_tag(path)
        if root_base is not None and self._base_energy_cb is not None:
            self._base_energy_cb.setCurrentIndex(
                next(i for i, (tag, _) in enumerate(_BASE_BEAMS) if tag == root_base))
        else:
            self._base_beam_tag = str(self._base_energy_cb.currentData())
        self._data   = load_all(path)
        self._refit_all_nl()
        self._cur    = None
        self._file_lbl.setText(str(path))
        self._mod_lbl.setText("← click a W module on the map")
        self._refit_peak_btn.setEnabled(False)
        self._refit_nl_btn.setEnabled(False)
        self._sync_peak_controls(None)
        self._status_lbl.setText("")
        self._reset_spans()

        if MODULES_JSON is not None:
            all_mods = load_modules(MODULES_JSON)
            self._w_modules = [m for m in all_mods if m.mod_type == "PbWO4"]
        else:
            # Fallback: build synthetic Module objects from ModuleData
            self._w_modules = [
                Module(d.name, "PbWO4", d.x, d.y, 20.8, 20.8)
                for d in self._data.values()
            ]

        self._map.set_modules(self._w_modules)
        self._refresh_map()
        self._clear_plots()
        self._after_load()

    # ── map ───────────────────────────────────────────────────────────────

    def _on_mode_changed(self):
        """Called only when the user switches the colour-map combo.
        Forces an auto range-reset."""
        self._current_map_mode = ""   # flag: treat as mode change
        self._refresh_map()

    def _outer_zero_layers(self) -> int:
        return int(self._outer_layer_spin.value())

    def _refit_all_nl(self):
        layers = self._outer_zero_layers()
        for d in self._data.values():
            d.refit_nl(layers, self._base_beam_tag)

    def _on_outer_layers_changed(self, _value: int):
        if not self._data:
            return
        self._refit_all_nl()
        self._refresh_map()
        self._refresh_current_fit_display(refit=False)
        self._status_lbl.setText(
            f"Outer {self._outer_zero_layers()} W layer(s) set to nl=0.")

    def _refresh_map(self):
        """Update map values.  Auto-resets the colour range only when the
        mode has actually changed; otherwise preserves the user-set range
        so manual edits survive data updates (peak refit, etc.)."""
        if not self._data:
            return
        mode = self._mode_cb.currentText()
        mode_changed = (mode != self._current_map_mode)
        self._current_map_mode = mode

        if "Non-linearity" in mode:
            values = {name: d.nl for name, d in self._data.items()}
            if mode_changed:
                vals = [v for v in values.values() if v != 0.0]
                if vals:
                    arr  = np.array(vals)
                    mean = float(np.mean(arr))
                    std  = float(np.std(arr))
                    self._range_ctrl.set_range(mean - 3.0 * std,
                                               mean + 3.0 * std)
        elif "data - 2nd-order fit" in mode:
            values = {
                name: d.mean_abs_data_fit_deviation_2nd
                for name, d in self._data.items()
            }
            if mode_changed:
                vals = [v for v in values.values() if math.isfinite(v)]
                vmax = float(np.percentile(np.array(vals), 95)) if vals else 0.01
                self._range_ctrl.set_range(0.0, max(vmax, 1e-6))
        elif "Base-energy e-e" in mode:
            values = {
                name: d.base_ee_corr_residual_2nd
                for name, d in self._data.items()
            }
            if mode_changed:
                vals = [v for v in values.values() if math.isfinite(v)]
                vmax = (float(np.percentile(np.abs(vals), 95))
                        if vals else 0.01)
                vmax = max(vmax, 1e-6)
                self._range_ctrl.set_range(-vmax, vmax)
        elif "order" in mode:
            attr = "mean_corr_resid_1st" if "1st" in mode else "mean_corr_resid_2nd"
            values = {name: getattr(d, attr) for name, d in self._data.items()}
            if mode_changed:
                vals = [v for v in values.values() if v > 0.0]
                vmax = float(np.percentile(np.array(vals), 95)) if vals else 0.01
                self._range_ctrl.set_range(0.0, vmax)
        else:
            values = {name: 0.0 for name in self._data}
            if mode_changed:
                self._range_ctrl.set_range(0.0, 1.0)
        self._map.set_values(values)

    def _sorted_module_names(self) -> List[str]:
        return sorted(self._data.keys(),
                      key=lambda n: int(n[1:]) if n[1:].isdigit() else 0)

    def _current_index(self, names: List[str], default: int) -> int:
        cur = self._cur.name if self._cur else ""
        return names.index(cur) if cur in names else default

    def _show_module_at(self, names: List[str], idx: int):
        """Highlight names[idx] on the map and select it like a click."""
        name = names[idx]
        self._map._hovered = name
        self._map.update()
        self._on_click(name)
        self._scan_lbl.setText(f"{name}  ({idx + 1}/{len(names)})")

    def _jump_module(self, step: int):
        names = self._sorted_module_names()
        if not names:
            return
        self._stop_scan()
        idx = self._current_index(names, -1 if step > 0 else 0)
        self._show_module_at(names, (idx + step) % len(names))

    def _shortcut_jump_module(self, step: int):
        if isinstance(QApplication.focusWidget(), QLineEdit):
            return
        self._jump_module(step)

    def keyPressEvent(self, event):
        focus = QApplication.focusWidget()
        if isinstance(focus, QLineEdit):
            super().keyPressEvent(event)
            return
        if not event.isAutoRepeat():
            key = event.key()
            scan = event.nativeScanCode()
            virt = event.nativeVirtualKey()
            if key == Qt.Key.Key_AltGr or scan in (108, 100) or virt in (108, 100):
                self._jump_module(1)
                event.accept()
                return
            if key == Qt.Key.Key_Alt and (scan in (64, 56) or virt in (64, 56)):
                self._jump_module(-1)
                event.accept()
                return
        super().keyPressEvent(event)

    def _on_click(self, name: str):
        if not name:
            return
        d = self._data.get(name)
        if d is None:
            return
        self._cur = d
        self._reset_spans()
        self._refit_peak_btn.setEnabled(True)
        self._refit_nl_btn.setEnabled(True)
        self._sync_peak_controls(d)
        self._refresh_current_fit_display("", refit=False)

    def _reset_spans(self):
        """Drop all drag-selected fit ranges."""
        self._span_range = {k: None for k in _PEAK_KEYS}
        self._range_lbl.setText("(auto — drag to set)")
        self._clear_range_btn.setEnabled(False)

    def _clear_span_range(self):
        """Clear all drag-selected fit ranges and refit the peaks automatically."""
        self._reset_spans()
        d = self._cur
        if d is None:
            return
        self._refresh_current_fit_display(
            "Fit ranges cleared; peaks refit with adaptive windows.",
            refit=d.auto_fit_peaks())

    # ── drawing ────────────────────────────────────────────────────────────

    def _clear_plots(self):
        canvases = [*self._hist_canvases.values()]
        if self._linearity_canvas is not None:
            canvases.append(self._linearity_canvas)
        for canvas in canvases:
            ax = canvas.ax
            ax.cla()
            ax.set_facecolor(THEME.CANVAS)
            canvas.redraw()

    def _draw_module(self, d: ModuleData):
        for key in _PEAK_KEYS:
            canvas = self._hist_canvases[key]
            ax = canvas.ax
            ax.cla()
            ax.set_facecolor(THEME.CANVAS)
            self._draw_energy_hist_single(ax, d, key)
            cfg = peak_cfg(key)
            canvas._style(ax,
                title=f"{cfg['label']}  ({cfg['energy_label']})",
                xlabel="E (MeV)", ylabel="Counts")
            canvas.redraw()

        if self._linearity_canvas is not None:
            ax = self._linearity_canvas.ax
            ax.cla()
            ax.set_facecolor(THEME.CANVAS)
            self._draw_linearity(ax, d)
            self._linearity_canvas._style(
                ax, title="Non-linearity",
                xlabel="E_rec (MeV)", ylabel="E_rec / E_exp")
            self._linearity_canvas.redraw()
        self._setup_span_selectors(d)

    def _draw_energy_hist_single(self, ax, d: ModuleData, key: str):
        """Draw histogram for one peak key in the expected-energy search range."""
        cfg = peak_cfg(key)
        hist = hist_for_key(d, key)
        if hist is None:
            ax.text(0.5, 0.5, "No data", ha="center", va="center",
                    transform=ax.transAxes, color=THEME.TEXT_DIM, fontsize=8)
            return

        counts, edges = hist
        centers = 0.5 * (edges[:-1] + edges[1:])
        color = str(cfg["color"])
        Eexp  = d.e_exp[key]
        sig   = d.sigma_exp[key]

        # Show a 50% wider display window for e-e peaks while keeping the
        # actual peak-search and fit range aligned with PhysicsTools::fitPeak.
        display_half_width = 0.3 if key.startswith("ee_") else 0.2
        xlo = max(float(edges[0]), (1.0 - display_half_width) * Eexp)
        xhi = min(float(edges[-1]), (1.0 + display_half_width) * Eexp)
        ax.set_xlim(xlo, xhi)

        ax.step(centers, counts, where="mid",
                color=THEME.TEXT, linewidth=0.9, alpha=0.85)
        visible = (centers >= xlo) & (centers <= xhi)
        ymax = float(counts[visible].max()) if np.any(visible) else float(counts.max())
        ax.set_ylim(0.0, max(1.0, ymax * 1.15))

        # Drag-selected span range band
        span = self._span_range.get(key)
        if span is not None:
            ax.axvspan(span[0], span[1], alpha=0.18,
                       color=THEME.ACCENT, linewidth=0, zorder=2)
            ax.axvline(span[0], color=THEME.ACCENT,
                       linewidth=1.2, linestyle="-", alpha=0.8, zorder=3)
            ax.axvline(span[1], color=THEME.ACCENT,
                       linewidth=1.2, linestyle="-", alpha=0.8, zorder=3)

        if Eexp > 0.0:
            # Candidate peak search range [0.8 Eexp, 1.2 Eexp].
            ax.axvspan(0.8 * Eexp, 1.2 * Eexp,
                       alpha=0.07, color=color, linewidth=0)
            # Expected energy marker
            ax.axvline(Eexp, color="orange", linewidth=1.5, linestyle="-", alpha=0.95)
            ax.text(0.97, 0.97, f"exp={Eexp:.0f}\nσ={sig:.0f}",
                    transform=ax.transAxes, fontsize=6, ha="right", va="top",
                    color=THEME.TEXT_DIM)

        # Fitted Gaussian overlay
        mu   = d.peaks[key]
        amp  = d.peak_amp[key]
        sigf = d.peak_sigma[key]
        manual_mu = d.manual_peaks.get(key, 0.0)
        if manual_mu <= 0.0 and mu > 0.0 and amp > 0.0 and sigf > 0.0:
            xg = np.linspace(max(xlo, mu - 4.0 * sigf),
                             min(xhi, mu + 4.0 * sigf), 300)
            yg = _gauss(xg, amp, mu, sigf)
            ax.plot(xg, yg, color=color, linewidth=2.0)
            ax.axvline(mu, color=color, linewidth=1.2, linestyle="--", alpha=0.7)
            ax.annotate(f"{mu:.1f} MeV",
                        xy=(mu, amp), xytext=(0, 6), textcoords="offset points",
                        color=color, fontsize=7, ha="center", va="bottom")
        elif Eexp > 0.0:
            ax.axvline(Eexp, color="orange", linewidth=1.5, linestyle="-", alpha=0.95)
        if manual_mu > 0.0:
            ax.axvline(manual_mu, color=color, linewidth=1.5, linestyle="--", alpha=0.9)
            ax.annotate(f"manual {manual_mu:.1f}",
                        xy=(manual_mu, ax.get_ylim()[1] * 0.82),
                        xytext=(0, 0), textcoords="offset points",
                        color=color, fontsize=7, ha="center", va="bottom")
        if not d.peak_enabled.get(key, True):
            ax.text(0.03, 0.92, "excluded", transform=ax.transAxes,
                    fontsize=7, color=THEME.TEXT_DIM, ha="left", va="top")

    def _draw_linearity(self, ax, d: ModuleData):
        ax.set_xlim(0, 4500)
        E_base = d.e_exp[f"ep_{self._base_beam_tag}"]
        pts = d.linearity_points(self._base_beam_tag)
        if not pts:
            ax.text(0.5, 0.5, "No valid peaks found", ha="center", va="center",
                    transform=ax.transAxes, color=THEME.TEXT_DIM, fontsize=9)
            return

        E_rec_arr = np.array([p[2] for p in pts])
        has_1st_coefficients = d.ndf_nl > 0 or d.nl_coefficients_manual
        has_2nd_coefficients = d.ndf_nl2 > 0 or d.nl_coefficients_manual

        ax.set_ylim(0.85, 1.15)

        xlo, xhi = ax.get_xlim()
        xplot = np.linspace(xlo, xhi, 500)

        ax.axhline(1.0, color="tomato", linestyle="--",
                   linewidth=1.5, alpha=0.75, label="y = 1 (ideal)")

        if E_base > 0.0 and len(pts) >= 2 and has_1st_coefficients:
            y1fit = nl_ratio_1st(xplot, d.nl, E_base)
            ax.plot(xplot, y1fit, color="cornflowerblue", linewidth=2.0,
                    label="1st order fit")

        if (E_base > 0.0 and (len(pts) >= 3 or d.nl_coefficients_manual)
                and has_2nd_coefficients):
            y2fit = nl_ratio_2nd(xplot, d.nl2_1, d.nl2_2, E_base)
            ax.plot(xplot, y2fit, color="violet", linewidth=2.0,
                    linestyle="--", label="2nd order fit")

        # Scatter: measured points
        for key, xe, ye_rec in pts:
            cfg = peak_cfg(key)
            ratio = ye_rec / xe if xe > 0 else 0.0
            ax.scatter(ye_rec, ratio, color=str(cfg["color"]), s=70, zorder=5,
                       edgecolors="white", linewidths=0.5)
            ax.annotate(str(cfg["label"]), xy=(ye_rec, ratio),
                        xytext=(6, 4), textcoords="offset points",
                        color=str(cfg["color"]), fontsize=6.5)

        # Corrected points — 1st order: E_corr = E_rec / (1 + nl*(E_rec-E_base)/1000)
        if E_base > 0.0 and len(pts) >= 2 and has_1st_coefficients:
            for key, xe, ye_rec in pts:
                denom = nl_ratio_1st(ye_rec, d.nl, E_base)
                E_corr = ye_rec / denom if denom != 0.0 else ye_rec
                ratio_corr = E_corr / xe if xe > 0 else 0.0
                ax.scatter(E_corr, ratio_corr, color="limegreen", s=55, zorder=6,
                           marker="s", edgecolors="white", linewidths=0.5)

        # Corrected points — 2nd order
        if (E_base > 0.0 and (len(pts) >= 3 or d.nl_coefficients_manual)
                and has_2nd_coefficients):
            for key, xe, ye_rec in pts:
                denom = nl_ratio_2nd(ye_rec, d.nl2_1, d.nl2_2, E_base)
                E_corr2 = ye_rec / denom if denom != 0.0 else ye_rec
                ratio_corr2 = E_corr2 / xe if xe > 0 else 0.0
                ax.scatter(E_corr2, ratio_corr2, color="orange", s=55, zorder=6,
                           marker="s", edgecolors="white", linewidths=0.5)

        # Dummy legend entries for point types
        ax.scatter([], [], color="tomato",         s=50, marker="o", label="e-p (meas.)")
        ax.scatter([], [], color="cornflowerblue", s=50, marker="o", label="e-e (meas.)")
        ax.scatter([], [], color="limegreen",      s=50, marker="s", label="corrected (1st)")
        ax.scatter([], [], color="orange",         s=50, marker="s", label="corrected (2nd)")

        # Fit result annotation
        if E_base > 0.0:
            chi1_str = (f"χ²/ndf={d.chi2_nl:.2f}/{d.ndf_nl}"
                        if d.ndf_nl > 0 else "")
            chi2_str = (f"χ²/ndf={d.chi2_nl2:.2f}/{d.ndf_nl2}"
                        if d.ndf_nl2 > 0 else "")
            if d.nl_coefficients_manual:
                note = (
                    f"MANUAL COEFFICIENTS\n"
                    f"1st: nl₁={d.nl:.6f}\n"
                    f"2nd: nl₁={d.nl2_1:.6f}  nl₂={d.nl2_2:.6f}\n"
                    f"{'vertex=3000 MeV\n' if d.nl2_vertex_fixed else ''}"
                    f"E_base={E_base:.1f} MeV "
                    f"({dict(_BASE_BEAMS)[self._base_beam_tag]} GeV)"
                )
            else:
                note = (
                    f"1st: nl₁={d.nl:.4f}±{d.nl_err:.4f}  {chi1_str}\n"
                    f"2nd: nl₁={d.nl2_1:.4f}±{d.nl2_1_err:.4f}"
                    f"  nl₂={d.nl2_2:.4f}±{d.nl2_2_err:.4f}\n"
                    f"     {chi2_str}\n"
                    f"E_base={E_base:.1f} MeV ({dict(_BASE_BEAMS)[self._base_beam_tag]} GeV)"
                )
            ax.text(0.03, 0.97, note,
                    transform=ax.transAxes, fontsize=7.0,
                    va="top", ha="left", color="lightcyan",
                    bbox=dict(boxstyle="round,pad=0.35",
                              fc="#1a2a3a", ec=THEME.BORDER, alpha=0.85))

        # Two side-by-side legends in the upper-right corner
        _leg_kw = dict(facecolor=THEME.PANEL, edgecolor=THEME.BORDER,
                       labelcolor=THEME.TEXT, fontsize=13)
        all_h, all_l = ax.get_legend_handles_labels()
        _fit_labels  = {"y = 1 (ideal)", "1st order fit", "2nd order fit"}
        _pt_labels   = {"e-p (meas.)", "e-e (meas.)", "corrected (1st)", "corrected (2nd)"}
        h_fit = [h for h, l in zip(all_h, all_l) if l in _fit_labels]
        l_fit = [l for l in all_l if l in _fit_labels]
        h_pt  = [h for h, l in zip(all_h, all_l) if l in _pt_labels]
        l_pt  = [l for l in all_l if l in _pt_labels]
        # Left legend: fits + reference  (anchored at 77% from left)
        leg_left = ax.legend(h_fit, l_fit,
                             loc="upper right", bbox_to_anchor=(0.77, 1.0),
                             **_leg_kw)
        ax.add_artist(leg_left)
        # Right legend: point types  (anchored at right edge)
        ax.legend(h_pt, l_pt,
                  loc="upper right", bbox_to_anchor=(1.0, 1.0),
                  **_leg_kw)

        y_values = []
        for line in ax.lines:
            x_data = np.asarray(line.get_xdata(), dtype=float)
            y_data = np.asarray(line.get_ydata(), dtype=float)
            visible = np.isfinite(x_data) & np.isfinite(y_data) & (x_data >= xlo) & (x_data <= xhi)
            if np.any(visible):
                y_values.append(y_data[visible])
        for collection in ax.collections:
            offsets = np.asarray(collection.get_offsets(), dtype=float)
            if offsets.ndim == 2 and offsets.shape[1] >= 2:
                visible = (
                    np.isfinite(offsets[:, 0]) & np.isfinite(offsets[:, 1])
                    & (offsets[:, 0] >= xlo) & (offsets[:, 0] <= xhi)
                )
                if np.any(visible):
                    y_values.append(offsets[visible, 1])
        if y_values:
            visible_y = np.concatenate(y_values)
            y_min, y_max = float(np.min(visible_y)), float(np.max(visible_y))
            y_span = max(y_max - y_min, 0.3)
            padding = 0.05 * y_span
            current_min, current_max = 0.85, 1.15
            ax.set_ylim(
                min(current_min, y_min - padding),
                max(current_max, y_max + padding),
            )

    # ── span selector ──────────────────────────────────────────────────────

    def _setup_span_selectors(self, d: ModuleData):
        """Attach SpanSelectors to all 6 peak histogram axes.
        Must be called after every canvas redraw (axes are recreated)."""
        self._span_selectors.clear()
        for key in _PEAK_KEYS:
            hist = hist_for_key(d, key)
            if hist is None:
                continue
            ax = self._hist_canvases[key].ax
            sel = SpanSelector(
                ax,
                lambda xmin, xmax, _key=key: self._on_span_select(_key, xmin, xmax),
                direction="horizontal",
                useblit=True,
                props=dict(alpha=0.25, facecolor=THEME.ACCENT),
                interactive=True,
                drag_from_anywhere=True,
            )
            self._span_selectors.append(sel)

    def _on_span_select(self, key: str, xmin: float, xmax: float):
        if xmax - xmin < 1.0:
            return
        d = self._cur
        if d is None:
            return
        self._span_range[key] = (xmin, xmax)
        self._range_lbl.setText(
            f"{peak_cfg(key)['label']}:  [{xmin:.1f},  {xmax:.1f}] MeV")
        self._clear_range_btn.setEnabled(True)
        self._refresh_current_fit_display(refit=self._fit_span_peak(d, key))

    # ── re-fit actions ─────────────────────────────────────────────────────

    def _fit_span_peak(self, d: ModuleData, key: str) -> bool:
        """Fit peak *key* in its drag-selected range; True on success. Does NOT redraw."""
        cfg = peak_cfg(key)
        hist = hist_for_key(d, key)
        if hist is None:
            self._status_lbl.setText(f"No histogram for {cfg['label']}.")
            return False
        if d.e_exp[key] <= 0.0:
            return False
        if not HAS_SCIPY:
            self._status_lbl.setText("Manual range fit requires scipy.")
            return False
        counts, edges = hist
        xmin, xmax = self._span_range[key]
        centers = 0.5 * (edges[:-1] + edges[1:])
        try:
            r = _gauss_fit_range(centers, counts, xmin, xmax,
                                 max((xmax - xmin) / 4.0, float(edges[1] - edges[0])),
                                 xmax - xmin, maxfev=3000)
        except Exception as exc:
            self._status_lbl.setText(f"Fit failed: {exc}")
            return False
        if r is None:
            self._status_lbl.setText(
                f"Not enough bins in [{xmin:.0f}, {xmax:.0f}] MeV — try wider.")
            return False
        mu, sigma_f, amp, mean_error = r
        range_note = f"  (drag [{xmin:.0f}, {xmax:.0f}] MeV)"
        if mu <= 0.0:
            self._status_lbl.setText(
                f"Fit FAILED for {cfg['label']}{range_note}.")
            return False
        d.peaks[key]      = mu
        d.peak_sigma[key] = sigma_f
        d.peak_amp[key]   = amp
        d.peak_mean_error[key] = mean_error
        self._status_lbl.setText(
            f"{cfg['label']}: μ={mu:.2f} MeV, σ={sigma_f:.2f} MeV{range_note}")
        return True

    def _do_refit_all_peaks(self):
        """Refit all 6 peaks with adaptive windows (clears manual ranges)."""
        d = self._cur
        if d is None:
            return
        self._reset_spans()
        ok = d.auto_fit_peaks()
        self._refresh_current_fit_display(
            "All 6 peaks refit with adaptive windows." if ok else None, refit=ok)

    def _do_refit_nl(self):
        d = self._cur
        if d is None:
            return
        self._refresh_current_fit_display()
        self._status_lbl.setText(
            f"Linearity refit:  nl = {d.nl:.4f} ± {d.nl_err:.4f},  "
            f"χ²/ndf = {d.chi2_nl:.2f}/{d.ndf_nl}")

    # ── export ────────────────────────────────────────────────────────────

    def _export_nl(self):
        """Export nl coefficients to a new calibration JSON.

        Workflow:
          1. Ask user to select a *base* calibration JSON (for factor/base_energy values).
          2. Patch W-module nl1/nl2/base_energy fields from current fits.
          3. Save as a new file (same dir, same name with _nonLinear suffix, or user-chosen path).
        """
        if not self._data:
            return

        # ── Step 1: choose base file ──────────────────────────────────────────
        calib_dir = database_dir(use_env=False) / "calibration"
        calib_dir = str(calib_dir if calib_dir.is_dir() else Path("."))
        src_path, _ = QFileDialog.getOpenFileName(
            self, "Select base calibration JSON", calib_dir,
            "JSON files (*.json);;All files (*)")
        if not src_path:
            return
        src_path = Path(src_path)

        try:
            with open(src_path) as f:
                entries = json.load(f)
        except Exception as exc:
            self._status_lbl.setText(f"Failed to read {src_path.name}: {exc}")
            return

        if not isinstance(entries, list):
            self._status_lbl.setText(
                f"{src_path.name} is not a calibration JSON array.")
            return

        named_entries = [e for e in entries if isinstance(e, dict) and "name" in e]
        if not named_entries:
            self._status_lbl.setText(
                f"{src_path.name} has no entries with a module name.")
            return
        calib_like_entries = [
            e for e in named_entries
            if "factor" in e and "base_energy" in e
        ]
        if not calib_like_entries:
            self._status_lbl.setText(
                f"{src_path.name} does not look like a calibration JSON (missing factor/base_energy).")
            return

        # ── Step 2: patch nl1 / nl2 ──────────────────────────────────────────
        n_patched = 0
        n_skipped = 0
        for entry in entries:
            if not isinstance(entry, dict):
                n_skipped += 1
                continue
            name = entry.get("name", "")
            if not isinstance(name, str) or not name.startswith("W"):
                continue
            d = self._data.get(name)
            if d is None:
                n_skipped += 1
                continue

            has_fit = d.ndf_nl2 > 0 or d.nl_coefficients_manual
            has_anchor = d.measured_peak(f"ep_{self._base_beam_tag}") > 0.0
            force_zero = has_anchor and (
                d.is_outer_layer(self._outer_zero_layers()) or d.is_absorber)
            if not (has_fit or force_zero):
                n_skipped += 1
                continue

            entry["base_energy"] = round(float(d.e_exp[f"ep_{self._base_beam_tag}"]), 6)
            if force_zero and not d.nl_coefficients_manual:
                nl1, nl2 = 0.0, 0.0
            elif d.fit_order == "1st":
                nl1, nl2 = d.nl, 0.0
            else:
                nl1, nl2 = d.nl2_1, d.nl2_2
            entry["nl1"] = round(float(nl1), 6)
            entry["nl2"] = round(float(nl2), 6)
            n_patched += 1

        if n_patched == 0:
            self._status_lbl.setText(
                f"No W modules were patched from {src_path.name}; check the base JSON and fit results.")
            return

        # ── Step 3: choose output path ────────────────────────────────────────
        stem = src_path.stem
        if not stem.endswith("_nonLinear"):
            stem += "_nonLinear"
        default_out = str(src_path.parent / (stem + ".json"))
        out_path, _ = QFileDialog.getSaveFileName(
            self, "Save new calibration JSON", default_out,
            "JSON files (*.json);;All files (*)")
        if not out_path:
            return

        try:
            with open(out_path, "w") as f:
                json.dump(entries, f, indent=2)
            self._status_lbl.setText(
                f"Exported nl1/nl2 for {n_patched} W modules; skipped {n_skipped}  →  {Path(out_path).name}")
        except Exception as exc:
            self._status_lbl.setText(f"Save failed: {exc}")

    # ── scan ──────────────────────────────────────────────────────────────

    def _on_scan_speed_changed(self, _idx: int):
        ms = self._scan_speed_cb.currentData()
        self._scan_timer.setInterval(ms)

    def _toggle_scan(self):
        if self._scan_timer.isActive():
            self._stop_scan()
        else:
            self._start_scan()

    def _start_scan(self):
        all_names = self._sorted_module_names()
        if not all_names:
            return
        self._scan_index = self._current_index(all_names, 0)
        self._scan_names = all_names
        self._scan_btn.setText("■ Stop")
        self._scan_timer.setInterval(self._scan_speed_cb.currentData())
        self._scan_timer.start()
        self._scan_step()   # show first module immediately

    def _stop_scan(self):
        if not self._scan_timer.isActive():
            return
        self._scan_timer.stop()
        self._scan_btn.setText("▶ Scan")
        self._scan_lbl.setText("")

    def _scan_step(self):
        names = self._scan_names
        if not names or self._scan_index >= len(names):
            self._stop_scan()
            return
        self._show_module_at(names, self._scan_index)
        total = len(names)
        self._scan_index += 1
        if self._scan_index >= total:
            self._scan_timer.stop()
            self._scan_btn.setText("▶ Scan")
            self._scan_lbl.setText(f"Done ({total} modules)")


# Entry point

def main():
    ap = argparse.ArgumentParser(description="HyCal PWO4 non-linearity analysis viewer")
    ap.add_argument("root_file", nargs="?", type=Path,
                    help="Path to nonlinearity_results.root (output of prad2ana_nonlinearity)")
    ap.add_argument("--theme", choices=available_themes(), default="dark")
    args = ap.parse_args()

    set_theme(args.theme)

    app = QApplication(sys.argv)
    app.setApplicationName("NL Viewer")
    win = NLViewerWindow(root_path=args.root_file)
    apply_theme_palette(win)
    win.show()
    sys.exit(app.exec())


if __name__ == "__main__":
    main()
