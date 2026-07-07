# Python Potential Flow Solver

A Python implementation of Thin Airfoil Theory for analyzing 2D cambered airfoils.

## Features

- Generate camber lines
- Compute camber slope
- Calculate Fourier coefficients
- Predict lift coefficient (Cl)
- Predict pitching moment coefficient (Cm)
- Compute circulation
- Visualize velocity vector fields
- Compare predictions with ANSYS CFD

## Technologies

- Python
- NumPy
- Matplotlib

## Project Structure

```
├── main.py
├── inputs.py
├── camber_line.py
├── camber_slope.py
├── fourier_coefficients.py
├── gamma_circulation.py
├── plotting.py
└── vector_field.py
```

## Results

The solver implements Thin Airfoil Theory and produces:

- Camber line
- Camber slope
- Fourier coefficients
- Lift coefficient
- Moment coefficient
- Circulation distribution
- Velocity vector field

The numerical predictions were compared against ANSYS Fluent simulations and showed good agreement for small angles of attack.

## Authors

Basil Cheeramthodika  
Indian Institute of Technology Bombay

