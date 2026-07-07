"""
camber_line.py
==============
Purpose  : Generate and plot the camber line y = f(x) for an airfoil.
Inputs   : NACA 4-digit string  OR  a math expression string in terms of x.
Outputs  : Arrays of x and y values; a plot of the camber line.
Assumptions: Normalised chord — leading edge at x=0, trailing edge at x=1.
"""

import numpy as np
import matplotlib.pyplot as plt


def naca4_camber(naca_str):
    """
    Return camber line function for a NACA 4-digit airfoil.

    NACA MPTT:
      m = M/100  (max camber)
      p = P/10   (position of max camber along chord)

    Camber formula:
      y = (m/p^2)(2px - x^2)            for 0 <= x <= p
      y = (m/(1-p)^2)(1-2p+2px-x^2)    for p < x <= 1
    """
    m = int(naca_str[0]) / 100.0   # maximum camber fraction
    p = int(naca_str[1]) / 10.0    # position of maximum camber

    if m == 0 or p == 0:
        return lambda x: np.zeros_like(np.asarray(x, dtype=float))

    def camber_func(x):
        x  = np.asarray(x, dtype=float)
        sc = (x.ndim == 0)
        x  = np.atleast_1d(x)
        y  = np.where(
            x <= p,
            (m / p**2) * (2*p*x - x**2),
            (m / (1-p)**2) * (1 - 2*p + 2*p*x - x**2)
        )
        return float(y[0]) if sc else y

    return camber_func


def custom_camber(expr_str):
    """
    Return a camber line function from a user math expression string.

    Expression must be in terms of x. Available: sin, cos, pi, sqrt, exp, tan.
    Examples:
      '0.08 * sin(pi * x)'
      '0.06 * x * (1 - x) * (1 - 2*x)'
      '0.07 * (1 - cos(pi * x)) / 2'
    """
    ns = {
        'sin': np.sin, 'cos': np.cos, 'tan': np.tan,
        'sqrt': np.sqrt, 'exp': np.exp, 'pi': np.pi, 'abs': np.abs,
    }

    def camber_func(x):
        env = dict(ns, x=x)
        return eval(expr_str, {"__builtins__": {}}, env)

    return camber_func


def generate_camber_line(camber_func, label='Airfoil', n_points=300):
    """
    Generate and plot camber line points y = f(x).

    Parameters
    ----------
    camber_func : callable   y = f(x)
    label       : str        Plot title label
    n_points    : int        Number of sample points

    Returns
    -------
    x_arr, y_arr : ndarrays of chord positions and camber heights
    """
    x_arr = np.linspace(0.0, 1.0, n_points)

    try:
        y_arr = np.asarray(camber_func(x_arr), dtype=float)
    except Exception:
        y_arr = np.array([float(camber_func(xi)) for xi in x_arr])

    # Plot
    fig, ax = plt.subplots(figsize=(9, 3))
    ax.plot(x_arr, y_arr, 'b-', lw=2, label=label)
    ax.axhline(0, color='gray', lw=0.8, linestyle='--', label='Chord line')
    ax.set_xlabel("x/c")
    ax.set_ylabel("y/c")
    ax.set_title(f"Camber Line — {label}")
    ax.set_xlim(0, 1)
    ax.set_aspect('equal')
    ax.legend()
    ax.grid(True)
    plt.tight_layout()
    plt.show()

    return x_arr, y_arr
