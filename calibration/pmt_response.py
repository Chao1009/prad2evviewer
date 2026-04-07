"""
PMT Gain Response Model
=======================
Models the relationship between a HyCal PMT's high voltage and its
observed signal strength (the right edge of the Bremsstrahlung
spectrum) as a power law:

    edge = A * V^k

which is a straight line in log-log space:

    log(edge) = log(A) + k * log(V)

A :class:`PMTGainModel` instance accumulates ``(Vmon, edge)``
measurements for one module, fits the power law by ordinary least
squares in log-log space, and proposes a ΔV that brings the next
measurement to a target ADC value.  When the fit cannot be trusted —
too few points, slope non-physical, or quality below ``r2_min`` — the
proposal falls back to a static lookup table keyed on the magnitude
of the ADC error.

Pure Python, no third-party dependencies, so the model is trivially
unit-testable in isolation.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import List, Optional, Tuple


# ----------------------------------------------------------------------
#  Power-law fit result (immutable)
# ----------------------------------------------------------------------

@dataclass(frozen=True)
class PMTFitResult:
    """Result of a power-law fit ``edge = A * V^k`` (linear in log-log).

    Attributes
    ----------
    a, k
        Power-law parameters in linear units.
    log_a
        Natural log of ``a`` — the intercept of the log-log line.
    n_points
        Number of measurements used in the fit.
    rss
        Residual sum of squares in log-log space.
    r_squared
        Coefficient of determination in log-log space.  Equals 1.0
        when ``n_points == 2`` (the line passes through both points).
    """
    a: float
    k: float
    log_a: float
    n_points: int
    rss: float
    r_squared: float

    def predict_edge(self, vmon: float) -> float:
        """Edge ADC predicted by the fit at voltage ``vmon`` (V > 0)."""
        return self.a * (vmon ** self.k)

    def predict_voltage(self, edge: float) -> float:
        """Voltage at which the fit predicts the given ``edge`` (edge > 0).

        Raises
        ------
        ValueError
            If the fit slope ``k`` is zero (a horizontal line cannot
            be inverted).  Callers should reject ``k <= 0`` fits via
            :meth:`PMTGainModel.is_good_fit` before calling this.
        """
        if self.k == 0:
            raise ValueError("predict_voltage: fit slope k is 0 (not invertible)")
        return math.exp((math.log(edge) - self.log_a) / self.k)


# ----------------------------------------------------------------------
#  Static lookup table (used when only one point is available)
# ----------------------------------------------------------------------

#: Step rules ``(min_diff_strict, dv_volts)`` checked top-down.  An
#: entry matches when ``|target - current| > min_diff_strict``.  The
#: final ``None`` entry is the catch-all and uses ``diff/20`` instead.
_LOOKUP_RULES: Tuple[Tuple[Optional[float], float], ...] = (
    (1000.0, 50.0),
    (500.0,  30.0),
    (200.0,  20.0),
    (100.0,   5.0),
    (None,    0.0),   # catch-all → diff / 20
)


def lookup_delta_v(target_edge: float, current_edge: float) -> float:
    """Static voltage-step lookup keyed on the absolute ADC error.

    Returns a ΔV (V), positive when the PMT needs more gain
    (``current_edge < target_edge``).  Step magnitudes:

    ::

        |diff| > 1000   →  50 V
        500 <  |diff|   →  30 V
        200 <  |diff|   →  20 V
        100 <  |diff|   →   5 V
            |diff| <=100 →  diff / 20 V
    """
    diff = target_edge - current_edge
    sign = 1.0 if diff >= 0 else -1.0
    ad = abs(diff)
    for threshold, step in _LOOKUP_RULES:
        if threshold is None:
            return round(sign * ad / 20.0, 1)
        if ad > threshold:
            return round(sign * step, 1)
    return 0.0  # unreachable; satisfies the type checker


# ----------------------------------------------------------------------
#  PMT gain model
# ----------------------------------------------------------------------

class PMTGainModel:
    """Accumulate ``(Vmon, edge)`` measurements and propose ΔV.

    One instance models one PMT.  Call :meth:`clear` between modules
    or when restarting a module from scratch (Redo).

    A fit is considered usable only if it is **physical** (the
    log-log slope ``k`` is positive — gain rises with voltage) and of
    **good quality** (``R² >= r2_min``, applied when there are at
    least 3 points; for 2 points the line is exact and R² carries no
    information).  When either check fails,
    :meth:`delta_v_to_target` falls back to the static lookup table
    on the most recent point.
    """

    #: Default minimum R² for a 3+ point fit to be trusted.
    DEFAULT_R2_MIN: float = 0.90

    def __init__(self, r2_min: float = DEFAULT_R2_MIN) -> None:
        self._points: List[Tuple[float, float]] = []
        self._fit: Optional[PMTFitResult] = None
        self.r2_min: float = r2_min

    # =====================================================================
    #  Data point management
    # =====================================================================

    def add_point(self, vmon: float, edge: float) -> bool:
        """Record one ``(vmon, edge)`` measurement.

        Both must be strictly positive and finite — the model lives
        in log space.  Non-positive, ``None``, ``NaN``, and ``±inf``
        inputs are silently dropped to keep callers on the happy
        path.  Returns ``True`` if the point was kept.
        """
        if vmon is None or edge is None:
            return False
        try:
            v, e = float(vmon), float(edge)
        except (TypeError, ValueError):
            return False
        if v <= 0 or e <= 0:
            return False
        if not (math.isfinite(v) and math.isfinite(e)):
            return False
        self._points.append((v, e))
        self._fit = None  # invalidate cached fit
        return True

    def clear(self) -> None:
        """Forget all measurements and the cached fit."""
        self._points.clear()
        self._fit = None

    # =====================================================================
    #  Read access (no computation)
    # =====================================================================

    @property
    def points(self) -> List[Tuple[float, float]]:
        """All recorded ``(vmon, edge)`` measurements (defensive copy)."""
        return list(self._points)

    @property
    def n_points(self) -> int:
        return len(self._points)

    @property
    def latest(self) -> Optional[Tuple[float, float]]:
        """The most recent ``(vmon, edge)`` point, or ``None`` if empty."""
        return self._points[-1] if self._points else None

    @property
    def fit(self) -> Optional[PMTFitResult]:
        """The cached fit, or ``None`` if none has been computed.

        Touching this property does **not** trigger a fit — call
        :meth:`linear_fit` if you need one computed on demand.  The
        cache is invalidated by :meth:`add_point` and :meth:`clear`.
        """
        return self._fit

    # =====================================================================
    #  Fitting
    # =====================================================================

    def linear_fit(self) -> Optional[PMTFitResult]:
        """Fit the power law to the recorded points.

        Performs ordinary least squares on ``(log V, log edge)``.
        Returns the fit, or ``None`` when there are fewer than 2
        points or all voltages are identical (slope undefined).
        Subsequent calls return the cached result until new points
        are added.
        """
        if self._fit is not None:
            return self._fit

        n = len(self._points)
        if n < 2:
            return None

        xs = [math.log(v) for v, _ in self._points]
        ys = [math.log(e) for _, e in self._points]
        mean_x = sum(xs) / n
        mean_y = sum(ys) / n
        sxx = sum((x - mean_x) ** 2 for x in xs)
        if sxx <= 0:
            return None  # all voltages identical → slope undefined
        sxy = sum((x - mean_x) * (y - mean_y) for x, y in zip(xs, ys))

        k = sxy / sxx
        log_a = mean_y - k * mean_x

        rss = sum((y - (log_a + k * x)) ** 2 for x, y in zip(xs, ys))
        tss = sum((y - mean_y) ** 2 for y in ys)
        r2 = 1.0 - rss / tss if tss > 0 else 1.0

        self._fit = PMTFitResult(
            a=math.exp(log_a),
            k=k,
            log_a=log_a,
            n_points=n,
            rss=rss,
            r_squared=r2,
        )
        return self._fit

    # =====================================================================
    #  Fit diagnostics
    # =====================================================================

    def _fit_rejection_reason(self,
                              fit: Optional[PMTFitResult]) -> Optional[str]:
        """Return a short rejection string, or ``None`` if the fit is OK.

        This is the single source of truth for fit acceptance — both
        :meth:`is_good_fit` and :meth:`delta_v_to_target` route through
        it so the rules cannot drift apart.
        """
        if fit is None:
            return "degenerate fit"
        if fit.k <= 0:
            return f"non-physical k={fit.k:.2f}"
        if fit.n_points >= 3 and fit.r_squared < self.r2_min:
            return (f"poor fit R²={fit.r_squared:.3f} "
                    f"< {self.r2_min:.2f}")
        return None

    def is_good_fit(self,
                    fit: Optional[PMTFitResult] = None) -> bool:
        """Return ``True`` if ``fit`` (or the cached fit) is usable.

        See :meth:`_fit_rejection_reason` for the precise rules.
        """
        if fit is None:
            fit = self._fit
        return self._fit_rejection_reason(fit) is None

    # =====================================================================
    #  ΔV proposal
    # =====================================================================

    def delta_v_to_target(self,
                          target_edge: float,
                          current_vmon: float,
                          current_edge: float) -> Tuple[float, str]:
        """Compute the voltage change needed to reach ``target_edge``.

        The current state is passed in explicitly rather than read
        from the latest accumulated point — this keeps the model a
        pure analyzer and lets the caller control the reference
        point (typically the just-measured ``(vmon, edge)`` pair).

        Parameters
        ----------
        target_edge
            The desired edge ADC value.
        current_vmon
            The voltage the PMT is currently at.  Used by the
            fit-based path as the reference for ``dv = v_target -
            current_vmon``.
        current_edge
            The most recently measured edge ADC.  Used by the
            lookup-table path as the reference for the ADC error.

        Returns
        -------
        (dv, mode_tag)
            ``dv`` — voltage step in V, rounded to 0.1 V to match
            HV crate granularity.  Positive when more gain is needed.
            ``mode_tag`` — short string for log output describing
            which method was selected and, if a fallback was used,
            why.

        Selection rule (based on the model's recorded points):

        * ``n_points >= 2`` with a fit accepted by
          :meth:`is_good_fit` → invert the fit at ``target_edge``
          and subtract ``current_vmon``.
        * ``n_points == 1`` → static lookup table on the
          ``(target_edge, current_edge)`` ADC error.
        * ``n_points == 0`` → :class:`RuntimeError`.  Callers must
          record at least one measurement first.

        A multi-point fit that fails :meth:`is_good_fit` (degenerate,
        non-physical, or poor R²) falls back to the lookup table,
        with ``mode_tag`` carrying the reason.
        """
        n = len(self._points)
        if n == 0:
            raise RuntimeError(
                "PMTGainModel.delta_v_to_target called with no data points"
            )

        if n == 1:
            return (lookup_delta_v(target_edge, current_edge),
                    "lookup (1 point)")

        fit = self.linear_fit()
        reason = self._fit_rejection_reason(fit)
        if reason is not None:
            dv = lookup_delta_v(target_edge, current_edge)
            return dv, f"lookup ({reason}, n={n})"

        # accepted fit — fit is non-None and physical
        v_target = fit.predict_voltage(target_edge)
        dv = round(v_target - current_vmon, 1)
        return dv, (f"fit n={n} k={fit.k:.2f} "
                    f"A={fit.a:.3g} R²={fit.r_squared:.3f}")

    # =====================================================================
    #  Dunder
    # =====================================================================

    def __len__(self) -> int:
        return len(self._points)

    def __bool__(self) -> bool:
        return bool(self._points)

    def __repr__(self) -> str:
        if self._fit is not None:
            return (f"PMTGainModel(n={len(self._points)}, "
                    f"k={self._fit.k:.3f}, A={self._fit.a:.3g})")
        return f"PMTGainModel(n={len(self._points)}, no fit)"
