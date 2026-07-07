"""
fourier_coefficients.py
=======================
Purpose  : Compute Fourier coefficients A0, A1, A2 and aerodynamic
           coefficients Cl and Cm using Thin Airfoil Theory.
Inputs   : Camber line function, angle of attack alpha [degrees].
Outputs  : A0, A1, A2, Cl, Cm_LE as a dict.
Assumptions:
  - Incompressible, inviscid, 2D flow.
  - Glauert substitution: x = (c/2)(1 - cos theta), theta in [0, pi].
  - Midpoint quadrature with 1000 points avoids leading-edge singularity.
  - Cm computed about the leading edge.
"""

import numpy as np
from camber_slope import compute_slope


def compute_coefficients(camber_func, alpha_deg, n_quad=1000, n_terms=20):
    """
    Compute A0, A1, A2, Cl, Cm for a given angle of attack.

    Glauert substitution: x = (c/2)(1 - cos theta)
    Midpoint quadrature (avoids singularity at theta = 0):
      theta_k = (k - 0.5) * d_theta,   k = 1 ... n_quad

    Fourier coefficients:
      A0 = alpha - (1/pi) * SUM (dy/dx)_k * d_theta
      An = (2/pi) * SUM (dy/dx)_k * cos(n * theta_k) * d_theta

    Aerodynamic coefficients:
      Cl    = pi * (2*A0 + A1)
      Cm_LE = -(pi/2) * (A0 + A1 - A2/2)

    Parameters
    ----------
    camber_func : callable   y = f(x)
    alpha_deg   : float      Angle of attack [degrees]
    n_quad      : int        Quadrature points (default 1000)
    n_terms     : int        Number of An terms (default 20)

    Returns
    -------
    dict : A0, A1, A2, Cl, Cm_LE, An (full array)
    """
    alpha_rad = np.deg2rad(alpha_deg)
    d_theta   = np.pi / n_quad

    # Midpoint theta values and corresponding chord positions
    theta  = np.linspace(d_theta / 2.0, np.pi - d_theta / 2.0, n_quad)
    x_vals = 0.5 * (1.0 - np.cos(theta))   # chord = 1, so c/2 = 0.5

    # Camber slope at each quadrature point
    slope = compute_slope(camber_func, x_vals)

    # A0
    A0 = alpha_rad - (1.0 / np.pi) * np.sum(slope) * d_theta

    # A1 to AN
    An = np.zeros(n_terms)
    for n in range(1, n_terms + 1):
        An[n-1] = (2.0 / np.pi) * np.sum(slope * np.cos(n * theta)) * d_theta

    A1 = An[0]
    A2 = An[1]

    Cl    = np.pi * (2.0 * A0 + A1)
    Cm_LE = -(np.pi / 2.0) * (A0 + A1 - A2 / 2.0)

    return {'A0': A0, 'A1': A1, 'A2': A2,
            'Cl': Cl, 'Cm_LE': Cm_LE, 'An': An}
