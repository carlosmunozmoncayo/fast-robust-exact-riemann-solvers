A Comparative Study of Iterative Riemann Solvers for the Shallow Water and Euler Equations
===========================================================================
This is a reproducibility repository for the paper **A Comparative Study of Iterative Riemann Solvers for the Shallow Water and Euler Equations** by Carlos Muñoz Moncayo, Manuel Quezada de Luna, and David I. Ketcheson.

### Dependencies
  - python (v3.8.9)
  - numpy (v1.18.4)
  - f2py (v1.22.3)
  - A Fortran compiler (GNU Fortran compiler v11.2.0)
  - Clawpack (v5.8.2)

### Subdirectories
The subdirectories **SWEs** and **EEs** contain code to produce the results presented in Sections 2 and 3 respectively. Each one is organized in the following way:

  - **Initial_guesses_comparison**: Code for Table 1 (SWEs) and Table 3 (EEs).
  - **Iterative_method_comparison**: Code for Table 2 (SWEs) and Table 4 (EEs).
  - **Finite_volumes_simulations**: Code for Tables 5, 7 (SWEs) and Tables 6, 8 (EEs).
