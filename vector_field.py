"""
vector_field.py
===============
Purpose  : Plot the 2-D velocity vector field around the airfoil (4c x 3c domain)
           using the 2-D Biot-Savart law, and compute total circulation via a
           closed velocity line integral.
Inputs   : Camber function, Fourier coefficients A0/An, V_inf, alpha [degrees].
Outputs  : Velocity field quiver plot; circulation Gamma from line integral.
Assumptions:
  - Airfoil modelled as 250 discrete point vortices along the camber line.
  - 2-D Biot-Savart: du = dG*(y-eta)/(2*pi*r^2), dv = -dG*(x-xi)/(2*pi*r^2)
  - Line integral contour: rectangle with 0.4c margin around airfoil.
  - Incompressible, inviscid, 2-D potential flow.
"""

import numpy as np
import matplotlib.pyplot as plt


def plot_vector_field(camber_func, A0, An, V_inf, alpha_deg,
                      chord=1.0, n_vortex=250, nx=55, ny=42, label='Vector Field'):
    """
    Plot velocity field (4c x 3c) and compute circulation via line integral.

    Steps:
      1. Place 250 vortex elements along camber line (cosine spacing).
      2. Build 4c x 3c grid: x in [-0.5c, 3.5c], y in [-1.5c, 1.5c].
      3. Biot-Savart sum at every grid point + add freestream.
      4. Plot colour-coded quiver. Overlay camber line.
      5. Compute circulation via rectangular velocity line integral.

    Parameters
    ----------
    camber_func : callable   y = f(x)
    A0          : float      Zeroth Fourier coefficient
    An          : ndarray    A1...AN
    V_inf       : float      Freestream velocity [m/s]
    alpha_deg   : float      Angle of attack [degrees]
    chord       : float      Chord length (default 1.0)
    n_vortex    : int        Number of vortex elements (default 250)
    nx, ny      : int        Grid resolution (default 55 x 42)
    label       : str        Plot title label

    Returns
    -------
    Gamma_line : float   Circulation from velocity line integral [m^2/s]
    """
    alpha_rad = np.deg2rad(alpha_deg)
    d_theta   = np.pi / n_vortex

    # ── Step 1: Discretise vortex sheet ──────────────────────────────────────
    theta_v = np.linspace(d_theta / 2.0, np.pi - d_theta / 2.0, n_vortex)
    xi      = (chord / 2.0) * (1.0 - np.cos(theta_v))   # x-positions

    try:
        eta = np.asarray(camber_func(xi), dtype=float)   # y-positions on camber
    except Exception:
        eta = np.array([float(camber_func(x)) for x in xi])

    # Vortex sheet strength at each element
    sin_v   = np.sin(theta_v)
    cos_v   = np.cos(theta_v)
    gamma_v = A0 * (1.0 + cos_v) / sin_v
    gamma_v += sum(An[n] * np.sin((n+1) * theta_v) for n in range(len(An)))
    gamma_v *= 2.0 * V_inf
    dGamma   = gamma_v * (chord / 2.0) * sin_v * d_theta   # circulation per element

    # ── Step 2: Build 4c x 3c grid ───────────────────────────────────────────
    x_lin = np.linspace(-0.5*chord, 3.5*chord, nx)
    y_lin = np.linspace(-1.5*chord, 1.5*chord, ny)
    X, Y  = np.meshgrid(x_lin, y_lin)

    # ── Step 3: Biot-Savart sum at every grid point ───────────────────────────
    Xf = X.ravel()[:, None]              # (nx*ny, 1)
    Yf = Y.ravel()[:, None]
    dX = Xf - xi[None, :]               # (nx*ny, n_vortex) x-separations
    dY = Yf - eta[None, :]              # y-separations
    r2 = np.maximum(dX**2 + dY**2, 1e-9)

    u_ind = np.sum(dGamma[None, :] * dY / (2.0 * np.pi * r2), axis=1)
    v_ind = np.sum(-dGamma[None, :] * dX / (2.0 * np.pi * r2), axis=1)

    Vx = (V_inf * np.cos(alpha_rad) + u_ind).reshape(ny, nx)
    Vy = (V_inf * np.sin(alpha_rad) + v_ind).reshape(ny, nx)

    # Mask points very close to camber line
    try:
        Y_cl = np.where((X >= 0) & (X <= chord),
                         camber_func(np.clip(X, 0, chord)), np.nan)
    except Exception:
        Y_cl = np.full_like(X, np.nan)
    mask = (X >= 0) & (X <= chord) & (np.abs(Y - Y_cl) < 0.025 * chord)
    Vx   = np.where(mask, np.nan, Vx)
    Vy   = np.where(mask, np.nan, Vy)
    Vmag = np.sqrt(Vx**2 + Vy**2)

    # ── Step 4: Plot ──────────────────────────────────────────────────────────
    fig, ax = plt.subplots(figsize=(13, 6))
    skip = 2
    q = ax.quiver(X[::skip, ::skip], Y[::skip, ::skip],
                   Vx[::skip, ::skip], Vy[::skip, ::skip],
                   Vmag[::skip, ::skip], cmap='coolwarm', width=0.002, alpha=0.9)
    fig.colorbar(q, ax=ax, label='|V| (m/s)')

    x_c = np.linspace(0, chord, 200)
    try:    y_c = camber_func(x_c)
    except: y_c = np.array([float(camber_func(x)) for x in x_c])
    ax.plot(x_c, y_c, 'k-', lw=2.5, label='Camber line')

    ax.set_xlim(-0.5, 3.5); ax.set_ylim(-1.5, 1.5)
    ax.set_xlabel("x/c");   ax.set_ylabel("y/c")
    ax.set_title(f"Velocity Field (4c x 3c) — {label}   α = {alpha_deg}°")
    ax.set_aspect('equal'); ax.legend(); ax.grid(True)
    plt.tight_layout(); plt.show()

    # ── Step 5: Circulation via velocity line integral ────────────────────────
    Gamma_line = _line_integral(xi, eta, dGamma, V_inf, alpha_rad, chord)
    return Gamma_line


def _line_integral(xi, eta, dGamma, V_inf, alpha_rad, chord,
                   margin=0.4, n_pts=30):
    """
    Compute circulation via closed rectangular velocity line integral.

    Gamma = CONTOUR_INTEGRAL V . dl  (counterclockwise)
          = INT_bottom Vx dx  +  INT_right Vy dy
          + INT_top    Vx dx  +  INT_left  Vy dy
    """
    xL = -margin * chord;            xR = chord + margin * chord
    yB = -(margin + 0.3) * chord;   yT = (margin + 0.3) * chord

    fu = V_inf * np.cos(alpha_rad)   # freestream x
    fv = V_inf * np.sin(alpha_rad)   # freestream y

    def velocity(px, py):
        dX = px[:, None] - xi[None, :]
        dY = py[:, None] - eta[None, :]
        r2 = np.maximum(dX**2 + dY**2, 1e-9)
        u  = fu + np.sum(dGamma[None, :] * dY / (2*np.pi*r2), axis=1)
        v  = fv + np.sum(-dGamma[None, :] * dX / (2*np.pi*r2), axis=1)
        return u, v

    G = 0.0

    px = np.linspace(xL, xR, n_pts)                 # bottom: left -> right
    u, _ = velocity(px, np.full(n_pts, yB))
    G += np.trapezoid(u, px)

    py = np.linspace(yB, yT, n_pts)                 # right: bottom -> top
    _, v = velocity(np.full(n_pts, xR), py)
    G += np.trapezoid(v, py)

    px = np.linspace(xR, xL, n_pts)                 # top: right -> left
    u, _ = velocity(px, np.full(n_pts, yT))
    G += np.trapezoid(u, px)

    py = np.linspace(yT, yB, n_pts)                 # left: top -> bottom
    _, v = velocity(np.full(n_pts, xL), py)
    G += np.trapezoid(v, py)

    return -G   # negate: CCW integral of CW vortex sheet is negative
