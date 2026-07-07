"""
main.py
=======
Purpose  : Main script for AE244 Assignment 2.
           Takes user inputs from inputs.py and calls all functions in order.

Usage    : python main.py


"""

import numpy as np

from camber_line import naca4_camber, custom_camber, generate_camber_line
from camber_slope import plot_slope, compute_slope
from fourier_coefficients import compute_coefficients
from gamma_circulation import compute_bound_circulation
from vector_field import plot_vector_field
from plotting import plot_coefficients, plot_gamma
from inputs import camber_func, LABEL,N_POINTS, DESIGNS, V_INF, ALPHA_DEG, ALPHA_MIN, ALPHA_MAX, N_ALPHA, cfd_data

print("=" * 55)
print(f"  AE244 Assignment 2 — Thin Airfoil Theory")
print(f"  Airfoil : {LABEL}")
print(f"  V_inf   : {V_INF} m/s    Alpha : {ALPHA_DEG} deg")
print("=" * 55)

# Camber line
print("\nCamber line plotted")
generate_camber_line(camber_func, label=LABEL, n_points=N_POINTS)

# Camber slope
print("Camber slope plotted")
plot_slope(camber_func, label=LABEL)

# Fourier coefficients, Cl and Cm vs alpha
print(" Fourier coefficients, Cl and Cm vs alpha")
plot_coefficients(camber_func,
                  alpha_min=ALPHA_MIN, alpha_max=ALPHA_MAX, n_alpha=N_ALPHA,
                  label=LABEL)

# Single-point coefficients at ALPHA_DEG
coeff = compute_coefficients(camber_func, ALPHA_DEG)
A0, A1, A2 = coeff['A0'], coeff['A1'], coeff['A2']
An = coeff['An']
print(f"\n  At α = {ALPHA_DEG}°:")
print(f"  A0={A0:.5f}  A1={A1:.5f}  A2={A2:.5f}")
print(f"  Cl={coeff['Cl']:.5f}  Cm_LE={coeff['Cm_LE']:.5f}")

# Circulation distribution + bound circulation
print(f"\nCirculation distribution and bound circulation at α = {ALPHA_DEG}°")
Gamma_dist = plot_gamma(A0, An, V_INF, label=LABEL )

# Velocity field + line integral
print(f"\nVelocity field and line integral at α = {ALPHA_DEG}° ...")
Gamma_line = plot_vector_field(camber_func, A0, An, V_INF, ALPHA_DEG, label=LABEL)

# Compare both circulation values
err = abs(Gamma_line - Gamma_dist) / abs(Gamma_dist) * 100
print(f"\nCirculation comparison at α = {ALPHA_DEG}°")
print(f"  Γ distribution integral  = {Gamma_dist:.5f}  m²/s")
print(f"  Γ velocity line integral = {Gamma_line:.5f}  m²/s")
print(f"  Difference               = {err:.3f} %")



print('\nDone.')