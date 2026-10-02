"""
Effective population size (N_e) estimation from LD decay.

Uses the Hill (1981) model: E[r²] = 1 / (1 + 4·N_e·c)
"""

import numpy as np
from scipy.optimize import curve_fit


def equilibrium_r2(distance_bp, ne, recomb_rate=4e-8):
    """
    Calculate expected r² under drift-recombination equilibrium.

    Parameters
    ----------
    distance_bp : float or array
        Physical distance in base pairs.
    ne : float
        Effective population size.
    recomb_rate : float
        Per-base-pair recombination rate.

    Returns
    -------
    float or array
        Expected r² at the given distance.
    """
    c = distance_bp * recomb_rate
    return 1.0 / (1.0 + 4.0 * ne * c)


def estimate_ne(distances, r2_values, recomb_rate=4e-8, n_samples=None, max_distance=None,
                p0=10000.0, upper=1e7, floor=1e-10):
    """
    Estimate N_e by nonlinear least squares on the Hill (1981) curve E[r2] = 1/(1 + 4 N_e c).

    Defaults reproduce the rice/soybean analyses (scripts 03, 27): r2 - 1/n correction when
    `n_samples` is given, all distance bins, p0 = 1e4, bounds (1, 1e7), corrected r2 floored at 1e-10.

    Parameters
    ----------
    distances : array-like    Bin mid-points (bp).
    r2_values : array-like    Mean r2 per bin (unadjusted).
    recomb_rate : float       Per-bp recombination rate.
    n_samples : int, optional Sample size for the r2 - 1/n correction.
    max_distance : float      Use only bins <= max_distance (bp). Default: all.

    Returns
    -------
    dict {'ne', 'ne_se', 'r2_predicted'}
    """
    d = np.asarray(distances, dtype=float)
    r2 = np.asarray(r2_values, dtype=float)
    mask = (d > 0) & np.isfinite(r2) & (r2 > 0)
    if max_distance is not None:
        mask &= d <= max_distance
    d, r2 = d[mask], r2[mask]
    if n_samples is not None and n_samples > 1:
        r2 = np.maximum(r2 - 1.0 / n_samples, floor)
    if len(d) < 3:
        raise ValueError("Insufficient data points for N_e estimation")

    def model(dist, ne):
        return 1.0 / (1.0 + 4.0 * ne * dist * recomb_rate)

    try:
        popt, pcov = curve_fit(model, d, r2, p0=[p0], bounds=(1, upper))
    except RuntimeError as e:
        raise RuntimeError(f"N_e fitting failed: {e}")
    return {"ne": float(popt[0]), "ne_se": float(np.sqrt(pcov[0, 0])), "r2_predicted": model(d, popt[0])}
