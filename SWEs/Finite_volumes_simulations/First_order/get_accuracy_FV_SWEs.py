from clawpack.pyclaw.solution import Solution
import matplotlib
matplotlib.use('Agg')
import matplotlib.pylab as pl
from matplotlib import rc
#rc('text', usetex=True)
import numpy as np
import os

solver_list=["HLLE","Roe","exact"]
path_list=["./_output_"+i for i in ["HLLE","Roe","exact"]]
spacing_list=[(40,80),(13,26),(4,8),(1,2)]

grid_list=[50*3**i for i in range(4)]

def get_last_h(path):
    sol=Solution(frame=10,read_aux=False,path=path,file_prefix='fort')
    h=np.copy(sol.state.q[0,:])
    x=sol.state.grid.x.centers
    nx=len(x)
    dx=x[1]-x[0]
    return h,x,dx,nx

def path_fine(solver):
    return("./_output_"+solver+"/4050")

def coarse_from_fine(fine,spacing,len_coarse):
   #spacing should be a tuple
   a,b=spacing
   coarse=np.array([fine[a+i*(b+1)] for i in range(len_coarse)])
   return coarse

def L2_grid_function_norm(x,dx):
    return np.sqrt(dx)*np.linalg.norm(x,ord=2)

def L_infty_grid_function_norm(x):
    return np.linalg.norm(x,ord=np.inf)

def main():
    print(f"L2 norm")
    for idx, grid in enumerate(grid_list):
        print(f"####Accuracy {grid} points####")
        for solver in solver_list:
            path="./_output_"+solver+"/"+str(grid)
            h,x,dx,nx = get_last_h(path)
            h_fine,x_fine,dx_fine,nx_fine = get_last_h(path_fine(solver))
            exact_coarse = coarse_from_fine(fine=h_fine,spacing=spacing_list[idx],len_coarse=len(h))
            diff = h-exact_coarse         
            print(f"{solver}: {100*round(L2_grid_function_norm(x=diff,dx=dx)/L2_grid_function_norm(x=exact_coarse,dx=dx),4)}%")
    print(f"L infty norm")
    for idx, grid in enumerate(grid_list):
        print(f"####Accuracy {grid} points####")
        for solver in solver_list:
            path="./_output_"+solver+"/"+str(grid)
            h,x,dx,nx=get_last_h(path)
            h_fine,x_fine,dx_fine,nx_fine=get_last_h(path_fine(solver))
            exact_coarse = coarse_from_fine(fine=h_fine,spacing=spacing_list[idx],len_coarse=len(h))
            diff=h-exact_coarse
            print(f"{solver}: {100*round(L_infty_grid_function_norm(x=diff)/L_infty_grid_function_norm(x=exact_coarse),4)}%")

if __name__=="__main__":
    main()



