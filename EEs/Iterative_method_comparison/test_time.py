####################
#IMPORTING MODULES
####################

import numpy as np
import random



import modules.ostrowski_newton_TS as ostrowski_newton_TS
import modules.ostrowski_TS as ostrowski_TS

import modules.newton_TS as newton_TS
import modules.gottlieb as gottlieb
import modules.van_leer_TS as van_leer_TS
import modules.mod_newton_TS as mod_newton_TS
import modules.quad_TR as quad_TR
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

gamma=1.4

tol=1.0E-6
print("Results tolerance 1.e-6")
print("Time newton_TS: ",min([newton_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time ostrowski_newton_TS: ",min([ostrowski_newton_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time ostrowski_TS: ",min([ostrowski_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time quad_TR: ",min([quad_TR.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time Gottlieb: ",min([gottlieb.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time Mod Newton TS",min([mod_newton_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time Van Leer TS: ",min([van_leer_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))

tol=1.0E-12
print("Results tolerance 1.e-12")
print("Time newton_TS: ",min([newton_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time ostrowski_newton_TS: ",min([ostrowski_newton_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time ostrowski_TS: ",min([ostrowski_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time quad_TR: ",min([quad_TR.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time Gottlieb: ",min([gottlieb.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time Mod Newton TS",min([mod_newton_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
print("Time Van Leer TS: ",min([van_leer_TS.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(7)]))
