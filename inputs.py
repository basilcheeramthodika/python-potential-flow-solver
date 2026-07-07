"""
inputs.py
=========
Purpose  : Single file where ALL user inputs are defined.
           Edit only this file before running main.py
"""

from camber_line import naca4_camber, custom_camber

# ── Primary airfoil ───────────────────────────────────────────────────────────
# Option 1: NACA 4-digit
NACA_NUMBER = '5412'
camber_func = naca4_camber(NACA_NUMBER)
LABEL       = f'NACA {NACA_NUMBER}'
N_POINTS    = 300

# Option 2: Custom expression — uncomment and comment out Option 1
# EXPRESSION  = '0.08 * sin(pi * x)'
# camber_func = custom_camber(EXPRESSION)
# LABEL       = f'Custom: {EXPRESSION}'

# ── Three custom designs for Section 4 ───────────────────────────────────────
DESIGNS = [
    ('High Camber Sine', custom_camber('0.08 * sin(pi * x)')),
    ('Reflexed S-shape', custom_camber('0.06 * x * (1 - x) * (1 - 2*x)')),
    ('Forward Loaded',   custom_camber('0.07 * (1 - cos(pi * x)) / 2')),
]

# ── Flight conditions ─────────────────────────────────────────────────────────
V_INF     = 20.0   # freestream velocity [m/s]
ALPHA_DEG = 3.0    # angle of attack for single-point plots [degrees]

# ── Alpha sweep range ─────────────────────────────────────────────────────────
ALPHA_MIN = -3.0
ALPHA_MAX = 12.0
N_ALPHA   = 50

# ── Optional CFD overlay ──────────────────────────────────────────────────────
# Replace None with numpy arrays if you have CFD data, e.g.:
# import numpy as np
# cfd_data = {
#     'alpha': np.array([-3, 0, 3, 6, 9, 12]),
#     'Cl':    np.array([0.25, 0.52, 0.78, 1.05, 1.30, 1.50]),
#     'Cm':    np.array([-0.16, -0.24, -0.31, -0.37, -0.43, -0.48])
# }
cfd_data = {'alpha': None, 'Cl': None, 'Cm': None}
