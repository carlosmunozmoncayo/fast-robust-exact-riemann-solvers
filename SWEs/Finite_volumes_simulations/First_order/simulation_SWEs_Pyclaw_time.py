
#!/usr/bin/env python
# encoding: utf-8

r"""
Shallow water flow
==================

Solve the one-dimensional shallow water equations:

.. math::
    h_t + (hu)_x & = 0 \\
    (hu)_t + (hu^2 + \frac{1}{2}gh^2)_x & = 0.

Here h is the depth, u is the velocity, and g is the gravitational constant.
The default initial condition used here models a dam break.
"""

from __future__ import absolute_import
import numpy as np
from clawpack import riemann
from clawpack.riemann.shallow_roe_with_efix_1D_constants import depth, momentum, num_eqn
import sw_exact_1D
import gc

def setup(outdir='./_output',solver_type='classic',riemann_solver='exact',mx=4050, order=1):

    from clawpack import pyclaw

    if riemann_solver.lower() == 'exact':
        rs = sw_exact_1D
    elif riemann_solver.lower() == 'hlle':
        rs = riemann.shallow_hlle_1D
    elif riemann_solver.lower() == 'roe':
        rs = riemann.shallow_roe_with_efix_1D
 
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(rs)
        solver.order = order
        solver.limiters = pyclaw.limiters.tvd.vanleer
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)

    solver.num_waves = 2
    solver.num_eqn = 2
    
    solver.cfl_desired = 0.9
    solver.cfl_max = 1.0
    #solver.dt_variable=False
    #solver.dt_initial=0.0001

    xlower = -5.0
    xupper = 5.0
    #mx = 500
    x = pyclaw.Dimension(xlower,xupper,mx,name='x')
    domain = pyclaw.Domain(x)
    state = pyclaw.State(domain,num_eqn)

    # Gravitational constant
    state.problem_data['grav'] = 1.0
    state.problem_data['dry_tolerance'] = 1e-3
    state.problem_data['sea_level'] = 0.0

    xc = state.grid.x.centers

    IC='WC'
    x0=0.

    if IC=='WC':
        #Woodward-Collela blast wave problem adapted for the SWEs
        hl = 30.
        ul = 0.
        hr = 50.
        ur = 0.
        hm=1.
        um=0.
        state.q[depth,:] = hl * (xc <= x0-2) + hr * (xc > x0+2)+ hm*(xc>x0-2)*(xc<=x0+2)
        state.q[momentum,:] = hl*ul * (xc <= x0-2) + hr*ur * (xc > x0+2)+ um*(xc>x0-2)*(xc<=x0+2)
        solver.bc_lower[0] = pyclaw.BC.wall
        solver.bc_upper[0] = pyclaw.BC.wall

    claw = pyclaw.Controller()
    claw.keep_copy = False
    claw.tfinal = 10.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.num_output_times= 1
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot
    claw.output_format = None

    return claw

def setplot(plotdata):
    return plotdata

if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)
    del output
    gc.collect()
