####################
#IMPORTING MODULES
####################

import numpy as np
import random



import modules.newton_TS as newton_TS
import modules.mod_newton_TS as mod_newton_TS
import modules.ostrowski_newton_TS as ostrowski_newton_TS
import modules.ostrowski_TS as ostrowski_TS
import modules.quad_TR as quad_TR


random.seed(1)

#as in the repo
lenmatrix_strong=2*10**6
lenmatrix_weak=lenmatrix_strong*4
lenmatrix=lenmatrix_strong+lenmatrix_weak


#Riemann problems with strong waves
#matrix_strong=[[10**random.uniform(-4,4),10**random.uniform(-4,4),10**random.uniform(-2,2),-random.uniform(-2,2)] for i in range(lenmatrix_strong)]
matrix_strong=[[10**random.uniform(-4,4),10**random.uniform(-4,4),10**random.uniform(-2,2),-10**random.uniform(-2,2)] for i in range(lenmatrix_strong)]
#Riemann problems with weak waves
matrix_weak=[[random.uniform(0.1,1),random.uniform(0.1,1),0,0] for i in range(lenmatrix_weak)]

matrix=matrix_strong+matrix_weak



grav=1.
conv_criteria=1
tol=1.0E-12
print("Results tol 1.e-12")
print("Time quad_TR: ",min([quad_TR.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time ostrowski_TS: ",min([ostrowski_TS.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time ostrowski_newton_TS: ",min([ostrowski_newton_TS.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_TS: ",min([newton_TS.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time mod_newton_TS: ",min([mod_newton_TS.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))

tol=1.0E-6
print("Results tol 1.e-6")
print("Time quad_TR: ",min([quad_TR.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time ostrowski_TS: ",min([ostrowski_TS.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time ostrowski_newton_TS: ",min([ostrowski_newton_TS.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_TS: ",min([newton_TS.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time mod_newton_TS: ",min([mod_newton_TS.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
