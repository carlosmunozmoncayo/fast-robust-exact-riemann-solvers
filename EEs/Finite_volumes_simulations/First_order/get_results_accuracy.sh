#!/bin/bash
rm -rf *_output_{exact,HLLE,Roe}

for method in {exact,HLLE,Roe}
do
    for grid in {50,150,450,1350,4050}
    do
        mkdir -p _output_${method}/${grid}
        ipython3 simulation_EE_Pyclaw_accuracy.py riemann_solver=${method} mx=${grid} outdir=./_output_${method}/${grid}
    done
done
