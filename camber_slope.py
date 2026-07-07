"""
camber_slope.py
===============
Purpose  : Compute and plot the slope of the camber line dy/dx at any x.
Inputs   : Camber line function y = f(x) from camber_line.py.
Outputs  : Slope value(s) dy/dx; a plot of dy/dx vs x/c.
Assumptions: Central difference with step dx=1e-5. Endpoints clamped to [1e-7, 1-1e-7].
"""

import numpy as np
import matplotlib.pyplot as plt


def compute_slope(camber_func, x, dx=1e-5):
    """
    Compute dy/dx at chord position(s) x using central finite difference.

    Formula: dy/dx = [ f(x + dx) - f(x - dx) ] / (2 dx)

    x+dx and x-dx are clamped to stay inside [0,1] to avoid
    evaluating the camber function outside the airfoil.

    Parameters
    ----------
    camber_func : callable       y = f(x)
    x           : float/ndarray  Chord position(s) in [0, 1]
    dx          : float          Step size (default 1e-5)

    Returns
    -------
    float or ndarray   dy/dx at the given x position(s)
    """
    x  = np.asarray(x, dtype=float)
    xp = np.clip(x + dx, 1e-7, 1.0 - 1e-7)   # forward point, clamped
    xm = np.clip(x - dx, 1e-7, 1.0 - 1e-7)   # backward point, clamped
    step = xp - xm                             # actual step after clamping

    try:
        return (camber_func(xp) - camber_func(xm)) / step
    except Exception:
        zp = np.array([float(camber_func(xi)) for xi in np.atleast_1d(xp)])
        zm = np.array([float(camber_func(xi)) for xi in np.atleast_1d(xm)])
        result = (zp - zm) / np.atleast_1d(step)
        return float(result[0]) if x.ndim == 0 else result


def plot_slope(camber_func, label='Camber Line Slope', n_points=400):
    """
    Plot dy/dx vs x/c for the given camber line.

    Parameters
    ----------
    camber_func : callable   y = f(x)
    label       : str        Plot title label
    n_points    : int        Number of evaluation points

    Returns
    -------
    x_arr, slope_arr : ndarrays
    """
    x_arr     = np.linspace(0.005, 0.995, n_points)   # interior only
    slope_arr = compute_slope(camber_func, x_arr)

    fig, ax = plt.subplots(figsize=(9, 3))
    ax.plot(x_arr, slope_arr, 'r-', lw=2, label='dy/dx')
    ax.axhline(0, color='gray', lw=0.8, linestyle='--')
    ax.set_xlabel("x/c")
    ax.set_ylabel("dy/dx")
    ax.set_title(f"Camber Line Slope — {label}")
    ax.set_xlim(0, 1)
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()

    return x_arr, slope_arr
