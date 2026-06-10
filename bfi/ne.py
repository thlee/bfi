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


def estimate_ne(distances, r2_values, recomb_rate=4e-8, n_samples=None):
    """
    Estimate N_e by nonlinear least-squares fitting of LD decay curve.

    Parameters
    ----------
    distances : array-like
        Physical distances in base pairs.
    r2_values : array-like
        Observed mean r² at each distance.
    recomb_rate : float
        Per-base-pair recombination rate.
    n_samples : int, optional
        Sample size for correction (r²_corrected = r² - 1/n).

    Returns
    -------
    dict
        {'ne': float, 'ne_se': float, 'r2_predicted': array}
    """
    distances = np.asarray(distances, dtype=float)
    r2_values = np.asarray(r2_values, dtype=float)

    # Sample size correction
    if n_samples is not None and n_samples > 1:
        r2_values = r2_values - 1.0 / n_samples
        r2_values = np.maximum(r2_values, 0.001)

    # Filter valid data
    mask = (distances > 0) & np.isfinite(r2_values) & (r2_values > 0)
    d = distances[mask]
    r2 = r2_values[mask]

    if len(d) < 3:
        raise ValueError("Insufficient data points for N_e estimation")

    def model(dist, ne):
        c = dist * recomb_rate
        return 1.0 / (1.0 + 4.0 * ne * c)

    try:
        popt, pcov = curve_fit(model, d, r2, p0=[1000], bounds=(1, 1e7))
        ne = popt[0]
        ne_se = np.sqrt(pcov[0, 0])
        r2_pred = model(d, ne)
        return {"ne": ne, "ne_se": ne_se, "r2_predicted": r2_pred}
    except RuntimeError as e:
        raise RuntimeError(f"N_e fitting failed: {e}")
