####################
#IMPORTING MODULES
####################
print("Importing modules...")
import numpy as np
import random



import modules.ostrowski_newton_TS as ostrowski_newton_TS
import modules.ostrowski_TS as ostrowski_TS

import modules.newton_TS as newton_TS
import modules.gottlieb as gottlieb
import modules.van_leer_TS as van_leer_TS
import modules.mod_newton_TS as mod_newton_TS
import modules.quad_TR as quad_TR
import modules.quad_below_TR as quad_below_TR
import modules.linear_TR as linear_TR
print("Modules imported")
random.seed(1)

#as in the repo
lenmatrix_strong=2*10**6
lenmatrix_weak=lenmatrix_strong*4
lenmatrix=lenmatrix_strong+lenmatrix_weak


####################
#INITIALIZING MATRIX
####################
print("Generating dataset of Riemann problems...")
#Harder problems
matrix=[[10**random.uniform(-2,2),10**random.uniform(-1,1),random.uniform(0.01,0.9),10**random.uniform(-2,2),-10**random.uniform(-1,1),random.uniform(0.01,0.9)] for i in range(lenmatrix_strong)]

#Adding easier problems
matrix=matrix+[[random.uniform(0.1,1),0.0,random.uniform(0.01,0.9),random.uniform(0.1,1),0.0,random.uniform(0.01,0.9)] for i in range(lenmatrix_weak)]
print("Dataset generated")

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
n_runs=7


modules=[quad_TR,newton_TS,ostrowski_newton_TS,ostrowski_TS,
gottlieb,mod_newton_TS,van_leer_TS,quad_below_TR,linear_TR]
names=["quad_TR","newton_TS","ostrowski_newton_TS","ostrowski_TS",
"gottlieb","mod_newton_TS","van_leer_TS","quad_below_TR","linear_TR"]


tol=1.0E-6
print("Results tolerance 1.e-6")
for name,module in zip(names,modules):
	print("############\nRunning "+name+"...")
	print("Time "+name+": ",min([module.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(n_runs)]))
	print("############")

tol=1.0E-12
print("Results tolerance 1.e-12")
for name,module in zip(names,modules):
	print("############\nRunning "+name+"...")
	print("Time "+name+": ",min([module.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix,n_data=lenmatrix) for i in range(n_runs)]))
	print("############")