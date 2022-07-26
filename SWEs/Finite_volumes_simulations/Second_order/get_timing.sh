#!/bin/bash
touch times.log
echo -n "" >times.log
for method in {hlle,roe,exact}
do
    echo -n "Execution time $method" >> times.log
    { time python3 simulation_SWEs_Pyclaw_time.py riemann_solver=$method > /dev/null ; } 2>> times.log
done
cat times.log
