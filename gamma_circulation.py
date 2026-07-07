"""
gamma_circulation.py
====================
Purpose  : Compute the vortex/circulation distribution gamma(theta) and
           the total bound circulation by analytically integrating gamma.
Inputs   : Fourier coefficients A0, An; freestream velocity V_inf.
Outputs  : Arrays of theta, x, gamma(theta); total bound circulation Gamma.
Assumptions:
  - gamma(theta) = 2*V_inf * [A0*(1+cos t)/sin t + SUM An*sin(n*t)]
  - Trailing edge Kutta condition: gamma(pi) = 0.
  - Leading edge singularity: gamma -> infinity as theta -> 0 (integrable).
  - Bound circulation: Gamma = pi * c * V_inf * (A0 + A1/2)
"""

import numpy as np


def compute_gamma(A0, An, V_inf, chord=1.0, n_points=500):
    """
    Compute vortex-sheet strength distribution gamma(theta).

    gamma(theta) = 2*V_inf * [ A0*(1+cos theta)/sin theta
                               + SUM_n An*sin(n*theta) ]

    Parameters
    ----------
    A0      : float      Zeroth Fourier coefficient
    An      : ndarray    A1...AN (index 0 = A1)
    V_inf   : float      Freestream velocity [m/s]
    chord   : float      Chord length (default 1.0)
    n_points: int        Number of evaluation points

    Returns
    -------
    theta : ndarray   Theta values [rad]
    x     : ndarray   Chord positions x/c
    gamma : ndarray   gamma(theta) [m/s]
    """
    # Avoid theta=0 (leading edge singularity)
    theta = np.linspace(np.pi / n_points, np.pi * (1 - 1.0/n_points), n_points)

    # A0 contribution: (1 + cos theta) / sin theta
    gamma = A0 * (1.0 + np.cos(theta)) / np.sin(theta)

    # An contributions: SUM An * sin(n*theta)
    for n, a_n in enumerate(An, start=1):
        gamma += a_n * np.sin(n * theta)

    gamma *= 2.0 * V_inf                              # scale by 2*V_inf
    x      = (chord / 2.0) * (1.0 - np.cos(theta))   # chord positions

    return theta, x, gamma


def compute_bound_circulation(A0, A1, V_inf, chord=1.0):
    """
    Total bound circulation by analytically integrating gamma.

    Gamma = INT_0^c gamma(x) dx = pi * c * V_inf * (A0 + A1/2)

    Consistency via Kutta-Joukowski: L = rho*V_inf*Gamma => Cl = 2*Gamma/(V_inf*c)

    Parameters
    ----------
    A0, A1 : Fourier coefficients
    V_inf  : float   Freestream velocity [m/s]
    chord  : float   Chord length

    Returns
    -------
    float   Total bound circulation Gamma [m^2/s]
    """
    return np.pi * chord * V_inf * (A0 + A1 / 2.0)
