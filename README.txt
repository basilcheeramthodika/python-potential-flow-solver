================================================================================
                    AE244 ASSIGNMENT 2 — THIN AIRFOIL THEORY TOOL
                              README AND USER GUIDE
================================================================================

TEAM
-------
Shaurya Gupta — 24B0051
Harshvardhan Verma — 24B0058
Basil Cheeramthodika — 24B0068

--------------------------------------------------------------------------------
OVERVIEW
--------------------------------------------------------------------------------

This program implements Thin Airfoil Theory (TAT) to analytically compute
aerodynamic properties of any 2-D airfoil defined by its camber line y = f(x).
It validates results against CFD (Ansys Fluent) simulations from Assignment 1
and explores custom airfoil designs.

--------------------------------------------------------------------------------
FILES
--------------------------------------------------------------------------------
  inputs.py               — Enter all the inputs here before running main.py
  camber_line.py          — Function (a): Generate and plot camber line y = f(x)
  camber_slope.py         — Function (b): Compute slope dy/dx at any chord point
  fourier_coefficients.py — Function (c): Compute A0, A1, A2, Cl, Cm via TAT
  gamma_circulation.py    — Compute gamma(theta) distribution and bound circulation
  vector_field.py         — Function (d): Velocity field plot + line integral
  plotting.py             — All matplotlib plots (plot_coefficients, plot_gamma)
  main.py                 — Function (e): Main script + all output


  README.txt              — This file


--------------------------------------------------------------------------------
REQUIREMENTS
--------------------------------------------------------------------------------

  Python   : 3.8 or higher
  NumPy    : pip install numpy
  Matplotlib : pip install matplotlib 
  Scipy    : pip install Scipy

  All libraries are free and open source.
  Install all at once:
      pip install numpy matplotlib scipy

--------------------------------------------------------------------------------
HOW TO RUN
--------------------------------------------------------------------------------

  Step 1 — Open main.py in any text editor.

  Step 2 — Edit the USER INPUTS section at the top of main.py:

      a) Primary airfoil:
         Option 1 (NACA):    Set NACA_NUMBER = '2412'  (or your airfoil)
         Option 2 (Custom):  Uncomment EXPRESSION line and comment out Option 1

      b) Flight conditions — set V_INF (m/s) and ALPHA_DEG (degrees).

      c) Alpha sweep range — set ALPHA_MIN, ALPHA_MAX, N_ALPHA.

      d) CFD data (optional) — set CFD_DATA as a dict with keys
         'alpha', 'Cl', 'Cm' (numpy arrays). Set to None to skip overlay.

         Example:
           CFD_DATA = {
               'alpha': np.array([-3, 0, 3, 6, 9, 12]),
               'Cl':    np.array([0.05, 0.32, 0.60, 0.87, 1.14, 1.40]),
               'Cm':    np.array([-0.04, -0.04, -0.04, -0.04, -0.04, -0.04])
           }

  Step 3 — Run the script from terminal:
      python main.py

  Step 4 — Plots will appear one at a time. Close each plot window to see
           the next one. Numerical results are printed to the terminal.


--------------------------------------------------------------------------------
OUTPUTS GENERATED (in order)
--------------------------------------------------------------------------------

  NACA Airfoil:
    Camber line y vs x
    Camber slope dy/dx vs x
    Fourier coefficients A0, A1, A2 vs alpha
    Cl vs alpha (with optional CFD overlay)
    Cm vs alpha (with optional CFD overlay)
    Gamma(theta) distribution — two subplots (vs theta and vs x/c)
    Bound circulation printed to terminal
    Velocity vector field (4c x 3c domain)
    Circulation from line integral printed to terminal
    Comparison of both circulation values printed to terminal



--------------------------------------------------------------------------------
CAMBER LINE EXPRESSION SYNTAX
--------------------------------------------------------------------------------

  Use standard Python math with these available functions:
    sin, cos, tan, sqrt, exp, pi, abs

  Variable must be named 'x'.
  x ranges from 0 (leading edge) to 1 (trailing edge).

  Examples:
    '0.08 * sin(pi * x)'
    '0.06 * x * (1 - x) * (1 - 2*x)'
    '0.07 * (1 - cos(pi * x)) / 2'
    '0.05 * x**2 * (1 - x)'

  NACA 4-digit camber is handled automatically by naca4_camber('2412').

--------------------------------------------------------------------------------
NOTES
--------------------------------------------------------------------------------

  * Velocity field and line integral computations involve a 250-element
    Biot-Savart sum over a 55x42 grid — expect ~10-20 seconds per plot.

  * The line integral uses a rectangular contour with 0.4c margin on all
    sides. The result should match the distribution integral to within 0.01%.

  * Thin Airfoil Theory assumes small angles, thin airfoils, and attached
    flow. Results will deviate from CFD at high angles of attack (stall).

  * Moment coefficient Cm is computed about the LEADING EDGE in this code


