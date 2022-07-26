
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

MIN_h = 1E10
MAX_u = -1E10 
MIN_dt = 1E10

def setup(outdir='./_output',solver_type='classic',riemann_solver='exact',mx=50):

    from clawpack import pyclaw

    if riemann_solver.lower() == 'exact':
        rs = sw_exact_1D
    elif riemann_solver.lower() == 'hlle':
        rs = riemann.shallow_hlle_1D
    elif riemann_solver.lower() == 'roe':
        rs = riemann.shallow_roe_with_efix_1D
 
    if solver_type == 'classic':
        solver = pyclaw.ClawSolver1D(rs)
        solver.order = 1
        solver.limiters = pyclaw.limiters.tvd.vanleer
    elif solver_type == 'sharpclaw':
        solver = pyclaw.SharpClawSolver1D(rs)

    solver.num_waves = 2
    solver.num_eqn = 2
    
    solver.cfl_desired = 0.9
    solver.cfl_max = 1.0
    #solver.max_steps=5000000

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

    state.problem_data['MIN_h'] = MIN_h
    state.problem_data['MIN_dt'] = MIN_dt
    state.problem_data['MAX_u'] = MAX_u

    
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
    claw.keep_copy = True
    claw.tfinal = 10.0
    claw.solution = pyclaw.Solution(state,domain)
    claw.solver = solver
    claw.outdir = outdir
    claw.setplot = setplot
    #claw.output_format = None

    return claw


#--------------------------
def setplot(plotdata):
#--------------------------
    """ 
    Specify what is to be plotted at each frame.
    Input:  plotdata, an instance of visclaw.data.ClawPlotData.
    Output: a modified version of plotdata.
    """ 
    plotdata.clearfigures()  # clear any old figures,axes,items data

    # Figure for depth
    plotfigure = plotdata.new_plotfigure(name='Water height', figno=0)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Water height'
    plotaxes.axescmd = 'subplot(311)'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = depth
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}

    # Figure for momentum[1]
    #plotfigure = plotdata.new_plotfigure(name='Momentum', figno=1)

    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(312)'
    plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Momentum'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = momentum
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}
    
    #Trying to plot velocity
    def velocity(current_data):
        return current_data.q[1, :] / current_data.q[0, :]
        
    # Set up for axes in this figure:
    plotaxes = plotfigure.new_plotaxes()
    plotaxes.axescmd = 'subplot(313)'
    plotaxes.xlimits = [-5.0,5.0]
    plotaxes.title = 'Velocity'

    # Set up for item on these axes:
    plotitem = plotaxes.new_plotitem(plot_type='1d')
    plotitem.plot_var = velocity
    plotitem.plotstyle = '-'
    plotitem.color = 'b'
    plotitem.kwargs = {'linewidth':3}
    
    return plotdata


if __name__=="__main__":
    from clawpack.pyclaw.util import run_app_from_main
    output = run_app_from_main(setup,setplot)

  

