"""
plotting.py
===========
Purpose  : All plotting functions for the Thin Airfoil Theory tool.
           Centralises every matplotlib call in one place.
Inputs   : Camber functions, Fourier coefficients, gamma arrays — from
           fourier_coefficients.py and gamma_circulation.py.
Outputs  : Matplotlib plots (displayed via plt.show()).
Assumptions: Same as fourier_coefficients.py and gamma_circulation.py.
"""

import numpy as np
import matplotlib.pyplot as plt
from fourier_coefficients import compute_coefficients
from gamma_circulation     import compute_gamma, compute_bound_circulation


def plot_coefficients(camber_func, alpha_min=-3.0, alpha_max=12.0,
                      n_alpha=50, label='Airfoil', cfd_data=None):
    """
    Plot A0/A1/A2 vs alpha, Cl vs alpha, and Cm vs alpha.

    Parameters
    ----------
    camber_func : callable   y = f(x)
    alpha_min   : float      Start angle [degrees]
    alpha_max   : float      End angle [degrees]
    n_alpha     : int        Number of angle steps
    label       : str        Legend label
    cfd_data    : dict       Optional CFD overlay — keys: 'alpha', 'Cl', 'Cm'
    """
    alphas = np.linspace(alpha_min, alpha_max, n_alpha)
    A0_arr, A1_arr, A2_arr = [], [], []
    Cl_arr, Cm_arr         = [], []

    for a in alphas:
        c = compute_coefficients(camber_func, a)
        A0_arr.append(c['A0']); A1_arr.append(c['A1']); A2_arr.append(c['A2'])
        Cl_arr.append(c['Cl']); Cm_arr.append(c['Cm_LE'])

    A0_arr = np.array(A0_arr); A1_arr = np.array(A1_arr); A2_arr = np.array(A2_arr)
    Cl_arr = np.array(Cl_arr); Cm_arr = np.array(Cm_arr)

    # ── Plot A0, A1, A2 vs alpha ─────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(alphas, A0_arr, 'b-o', ms=3, lw=1.8, label='A0')
    ax.plot(alphas, A1_arr, 'r-s', ms=3, lw=1.8, label='A1')
    ax.plot(alphas, A2_arr, 'g-^', ms=3, lw=1.8, label='A2')
    ax.axhline(0, color='gray', lw=0.7, linestyle='--')
    ax.set_xlabel("α (deg)"); ax.set_ylabel("Fourier Coefficient")
    ax.set_title(f"Fourier Coefficients A0, A1, A2 vs α — {label}")
    ax.legend(); ax.grid(True)
    plt.tight_layout(); plt.show()

    # ── Plot Cl vs alpha ─────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(alphas, Cl_arr, 'b-', lw=2, label=f'TAT — {label}')
    if cfd_data:
        ax.plot(cfd_data['alpha'], cfd_data['Cl'], 'ro--', ms=5, lw=1.5, label='CFD')
    ax.axhline(0, color='gray', lw=0.6, linestyle='--')
    ax.set_xlabel("α (deg)"); ax.set_ylabel("Cl")
    ax.set_title(f"Cl vs α — {label}"); ax.legend(); ax.grid(True)
    plt.tight_layout(); plt.show()

    # ── Plot Cm vs alpha ─────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(9, 4))
    ax.plot(alphas, Cm_arr, 'r-', lw=2, label=f'TAT Cm_LE — {label}')
    if cfd_data and 'Cm' in cfd_data:
        ax.plot(cfd_data['alpha'], cfd_data['Cm'], 'bo--', ms=5, lw=1.5, label='CFD')
    ax.axhline(0, color='gray', lw=0.6, linestyle='--')
    ax.set_xlabel("α (deg)"); ax.set_ylabel("Cm (about LE)")
    ax.set_title(f"Cm vs α — {label}"); ax.legend(); ax.grid(True)
    plt.tight_layout(); plt.show()

    return alphas, Cl_arr, Cm_arr


def plot_gamma(A0, An, V_inf, label='Airfoil', chord=1.0):
    """
    Plot gamma(theta) vs theta and gamma vs x/c.
    Also computes and prints the total bound circulation.

    Parameters
    ----------
    A0    : float     Zeroth Fourier coefficient
    An    : ndarray   A1...AN
    V_inf : float     Freestream velocity [m/s]
    label : str       Plot title label
    chord : float     Chord length

    Returns
    -------
    Gamma : float   Total bound circulation [m^2/s]
    """
    theta, x, gamma = compute_gamma(A0, An, V_inf, chord)

    fig, axes = plt.subplots(1, 2, figsize=(12, 4))

    # ── gamma vs theta ───────────────────────────────────────────────────────
    axes[0].plot(np.rad2deg(theta), gamma / V_inf, 'm-', lw=2)
    axes[0].set_xlabel("θ (deg)")
    axes[0].set_ylabel("γ(θ) / V∞")
    axes[0].set_title("Circulation Distribution γ(θ) vs θ")
    axes[0].set_xlim(0, 180)
    axes[0].grid(True)

    # ── gamma vs x/c ─────────────────────────────────────────────────────────
    axes[1].plot(x, gamma / V_inf, 'c-', lw=2)
    axes[1].set_xlabel("x/c")
    axes[1].set_ylabel("γ(x) / V∞")
    axes[1].set_title("Circulation Distribution γ(x) vs x/c")
    axes[1].set_xlim(0, 1)
    axes[1].grid(True)

    fig.suptitle(f"Circulation Distribution — {label}")
    plt.tight_layout()
    plt.show()

    Gamma = compute_bound_circulation(A0, An[0], V_inf, chord)
    print(f"  Bound Circulation Γ (distribution integral) = {Gamma:.5f}  m²/s")

    return Gamma
