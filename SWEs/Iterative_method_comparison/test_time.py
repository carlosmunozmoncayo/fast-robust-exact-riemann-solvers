####################
#IMPORTING MODULES
####################
print("Importing modules...")
import numpy as np
import random
import modules.newton_TS as newton_TS
import modules.mod_newton_TS as mod_newton_TS
import modules.ostrowski_newton_TS as ostrowski_newton_TS
import modules.ostrowski_TS as ostrowski_TS
import modules.quad_TR as quad_TR
import modules.quad_below_TR as quad_below_TR
import modules.linear_TR as linear_TR
print("Done importing modules.")

random.seed(1)

#as in the repo
lenmatrix_strong=2*10**6
lenmatrix_weak=lenmatrix_strong*4
lenmatrix=lenmatrix_strong+lenmatrix_weak

print("Generating dataset of Riemann problems...")
#Riemann problems with strong waves
#matrix_strong=[[10**random.uniform(-4,4),10**random.uniform(-4,4),10**random.uniform(-2,2),-random.uniform(-2,2)] for i in range(lenmatrix_strong)]
matrix_strong=[[10**random.uniform(-4,4),10**random.uniform(-4,4),10**random.uniform(-2,2),-10**random.uniform(-2,2)] for i in range(lenmatrix_strong)]
#Riemann problems with weak waves
matrix_weak=[[random.uniform(0.1,1),random.uniform(0.1,1),0,0] for i in range(lenmatrix_weak)]

matrix=matrix_strong+matrix_weak
print("Dataset generated")


modules=[quad_TR,newton_TS,ostrowski_newton_TS,ostrowski_TS,mod_newton_TS,quad_below_TR,linear_TR]
names=["quad_TR","newton_TS","ostrowski_newton_TS","ostrowski_TS","mod_newton_TS","quad_below_TR","linear_TR"]
n_runs=7
grav=1.
conv_criteria=1
tol=1.0E-12
print("Results tol 1.e-12")
for name, modules in zip(names,modules):
    print("############\nRunning "+name+"...")
    print("Time "+name+": ",min([modules.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(n_runs)]))
    print("############")
tol=1.0E-6
print("Results tol 1.e-6")
for name, modules in zip(names,modules):
    print("############\nRunning "+name+"...")
    print("Time "+name+": ",min([modules.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(n_runs)]))
    print("############")