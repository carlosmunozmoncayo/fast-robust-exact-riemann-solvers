####################
#IMPORTING MODULES
####################

import numpy as np
import random


import modules.newton_AV as newton_AV
import modules.newton_CC as newton_CC
import modules.newton_PV as newton_PV
import modules.newton_TR as newton_TR
import modules.newton_TS as newton_TS
import modules.newton_HLLE as newton_HLLE
import modules.adaptive as adaptive


random.seed(1)

#as in the repo
lenmatrix_strong=2*10**6
lenmatrix_weak=lenmatrix_strong*4
lenmatrix=lenmatrix_strong+lenmatrix_weak


####################
#INITIALIZING MATRIX
####################

#Harder problems
matrix=[[10**random.uniform(-2,2),10**random.uniform(-1,1),random.uniform(0.01,0.9),10**random.uniform(-2,2),-10**random.uniform(-1,1),random.uniform(0.01,0.9)] for i in range(lenmatrix_strong)]

#Adding easier problems
matrix=matrix+[[random.uniform(0.1,1),0.0,random.uniform(0.01,0.9),random.uniform(0.1,1),0.0,random.uniform(0.01,0.9)] for i in range(lenmatrix_weak)]

i=0
while False:
	pl=random.uniform(0.1,1)
	pr=random.uniform(0.1*pl,pl)
	rhol=random.uniform(0.1,1)
	rhor=random.uniform(0.1,1)
	cl=np.sqrt(pl/(rhol*1.4))
	cr=np.sqrt(pr/(rhor*1.4))
	if 0.05<(cr/(cl+cr)) and (cr/(cl+cr))<0.95:
		matrix.append([pl,0,rhol,pr,0,rhor])
		i+=1
	if i==lenmatrix:
		break

tol=1.0E-12
gamma=1.4

tol=1.0E-6
print("Results tolerance 1.e-6")
print("Time newton_TS: ",min([newton_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_TR: ",min([newton_TR.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_AV: ",min([newton_AV.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_PV: ",min([newton_PV.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_CC: ",min([newton_CC.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_HLLE: ",min([newton_HLLE.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time adaptive: ",min([adaptive.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))

tol=1.0E-12
print("Results tolerance 1.e-12")
print("Time newton_TS: ",min([newton_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_TR: ",min([newton_TR.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_AV: ",min([newton_AV.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_PV: ",min([newton_PV.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_CC: ",min([newton_CC.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time newton_HLLE: ",min([newton_HLLE.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time adaptive: ",min([adaptive.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))

