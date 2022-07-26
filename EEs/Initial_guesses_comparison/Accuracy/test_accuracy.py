####################
#IMPORTING MODULES
####################

import numpy as np
import random



import modules.accuracy as accuracy


random.seed(1)

#as in the repo
lenmatrix_strong=2*10**6
lenmatrix_weak=lenmatrix_strong*4
lenmatrix=lenmatrix_strong+lenmatrix_weak


####################
#INITIALIZING MATRIX
####################

#Harder problems
matrix_strong=[[10**random.uniform(-2,2),10**random.uniform(-1,1),random.uniform(0.01,0.9),10**random.uniform(-2,2),-10**random.uniform(-1,1),random.uniform(0.01,0.9)] for i in range(lenmatrix_strong)]

#Adding easier problems
matrix_weak=[[random.uniform(0.1,1),0.0,random.uniform(0.01,0.9),random.uniform(0.1,1),0.0,random.uniform(0.01,0.9)] for i in range(lenmatrix_weak)]




gamma=1.4
tol=1.e-12

print(f"#####\nTesting accuracy EEs with strong waves\n#####")
accuracy.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix_strong,n_data=lenmatrix_strong)

print(f"#####\nTesting accuracy EEs with weak waves\n#####")
accuracy.pstar(tol=tol,gamma=gamma,conv_criteria=1,rp_data=matrix_weak,n_data=lenmatrix_weak)

