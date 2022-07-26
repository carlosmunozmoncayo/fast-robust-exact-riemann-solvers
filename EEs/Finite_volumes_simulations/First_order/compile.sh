#/bin/bash

python3 -m numpy.f2py -c rp1_euler_newton.f90 -m euler_exact_1D
