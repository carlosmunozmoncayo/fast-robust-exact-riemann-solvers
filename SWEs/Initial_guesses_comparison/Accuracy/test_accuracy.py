####################
#IMPORTING MODULES
####################

import numpy as np
import random

#importing approximate solvers for accuracy testing
import modules.accuracy as accuracy

random.seed(1)

lenmatrix_strong=2*10**6
lenmatrix_weak=lenmatrix_strong*4
lenmatrix=lenmatrix_strong+lenmatrix_weak


#Riemann problems with strong waves
matrix_strong=[[10**random.uniform(-4,4),10**random.uniform(-4,4),10**random.uniform(-2,2),-10**random.uniform(-2,2)] for i in range(lenmatrix_strong)]

#Riemann problems with weak waves
matrix_weak=[[random.uniform(0.1,1),random.uniform(0.1,1),0,0] for i in range(lenmatrix_weak)]

#For strong waves
print(f"#####\nTesting accuracy SWEs with strong waves\n#####")
accuracy.hstar(tol=1.e-12,grav=1.0,rp_data=matrix_strong,n_data=lenmatrix_strong)

print(f"#####\nTesting accuracy SWEs with weak waves\n#####")
accuracy.hstar(tol=1.e-12,grav=1.0,rp_data=matrix_weak,n_data=lenmatrix_weak)
