####################
#IMPORTING MODULES
####################

import numpy as np
import random


import modules.newton_AV as newton_AV
import modules.newton_CC as newton_CC
import modules.newton_QA as newton_QA
import modules.newton_TR as newton_TR
import modules.newton_TS as newton_TS
import modules.newton_PV as newton_PV
import modules.newton_HLLE as newton_HLLE



random.seed(1)

#as in the repo
lenmatrix_strong=2*10**6
lenmatrix_weak=lenmatrix_strong*4
lenmatrix=lenmatrix_strong+lenmatrix_weak


#Riemann problems with strong waves
matrix_strong=[[10**random.uniform(-4,4),10**random.uniform(-4,4),10**random.uniform(-2,2),-10**random.uniform(-2,2)] for i in range(lenmatrix_strong)]

#Riemann problems with weak waves
matrix_weak=[[random.uniform(0.1,1),random.uniform(0.1,1),0,0] for i in range(lenmatrix_weak)]

matrix=matrix_strong+matrix_weak




grav=1.
conv_criteria=1

tol=1.0E-12
print("Results tol 1.e-12")
print("Time newton_CC: ",min([newton_CC.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_TS: ",min([newton_TS.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_TR: ",min([newton_TR.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_AV: ",min([newton_AV.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_QA: ",min([newton_QA.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_PV: ",min([newton_PV.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_HLLE: ",min([newton_HLLE.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))


tol=1.0E-6
print("Results tol 1.e-6")
print("Time newton_CC: ",min([newton_CC.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_TS: ",min([newton_TS.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_TR: ",min([newton_TR.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_AV: ",min([newton_AV.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_QA: ",min([newton_QA.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_PV: ",min([newton_PV.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_HLLE: ",min([newton_HLLE.hstar(tol=tol,grav=grav,conv_criteria=conv_criteria,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
