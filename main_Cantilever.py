import numpy as np
import os
import time
from misc_functions import GetEigenVal
from abaqus_functions import GenerateRealization,AssignPropertiesRun,ReadSaveOutput
#from matplotlib import pyplot as plt

# Number of Monte Carlo Simulations
NMCS = 2

# Initial parameters

Lw= 52             # moving window size
step = 1           # moving window step as a fraction of window size
Lx= 100            # domain size x-direction
Ly= 100            # domain size y-direction
p = 1              # Matern smoothness parameter, nu = p + 1/2
n = 8              # Number of random field components
element_size = 10  # (stochastic) element size, usually identical to finite element size

# Load random field data #
filename = 'InputData\Params' + str(Lw) + '.txt'
Params = np.loadtxt(filename,delimiter=',')
bx_list = Params[:,2]
by_list = Params[:,3]

filename = 'InputData\Rcoeff_Lw_' + str(Lw) + '.txt'
rcoeff_matrix = np.loadtxt(filename,delimiter=',')

# Start Monte Carlo loop
start_time = time.time()
[eigenvalue_vector,eigenvector_matrix] = GetEigenVal(element_size,Lx,Ly,step,n,p,bx_list,by_list,rcoeff_matrix) # compute eigenvalues for random field generations, runs once
print("Elapsed time: {:.2f} seconds".format(time.time() - start_time))
for ID in range(1,NMCS+1):
    start_time = time.time()
    ABD_list = GenerateRealization(eigenvalue_vector,eigenvector_matrix,Params,n) # Generate a realization of the random field of ABD matrix components
    AssignPropertiesRun(ABD_list,Lx,Ly,element_size) # Creates a model and assigns random ABD properties to elements
    ReadSaveOutput(ID,Lx,Ly,Lw)
    print("Elapsed time: {:.2f} seconds".format(time.time() - start_time))
    print("ITERATION " + str(ID) + " COMPLETED")



