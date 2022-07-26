#!/bin/bash

python3	 -m numpy.f2py -c rp1_shallow_exact.f90 -m sw_exact_1D
