
from math import factorial
import numpy as np
from numpy import linalg as la

def expMattern(p,x,y,bx,by):
    r = np.sqrt( (x/bx)**2 + (y/by)**2)       
    
    p1 = np.exp(-np.sqrt(2*p+1)*r)
    p2 = factorial(p)/float(factorial(2*p)) # float for numerical errors in abaqus
    sump= 0
    for i in range(0,p+1):
        p3 = factorial(p+i)/(factorial(i)*factorial(p-i))
        p4 = (2*r*np.sqrt(2*p+1))**(p-i)
        sump = p3*p4 + sump
    return p1*p2*sump

def nearestPD(A):

    B = (A + A.T) / 2
    _, s, V = la.svd(B)
    H = np.dot(V.T, np.dot(np.diag(s), V))
    A2 = (B + H) / 2
    A3 = (A2 + A2.T) / 2
    if isPD(A3):
        return A3
    spacing = np.spacing(la.norm(A))
    I = np.eye(A.shape[0])
    k = 1
    while not isPD(A3):
        mineig = np.min(np.real(la.eigvals(A3)))
        A3 += I * (-mineig * k**2 + spacing)
        k += 1
    print("Max imag eig:", np.max(np.abs(np.imag(la.eigvals(A3)))))
    print("Min real eig:", np.min(np.real(la.eigvals(A3))))
    print("Symmetry error:", np.max(np.abs(A3 - A3.T)))
    return A3

def isPD(B):
    """Returns true when input is positive-definite, via Cholesky"""
    try:
        _ = la.cholesky(B)
        #_ = la.cholesky(B + np.eye(B.shape[0]) * 1e-12) # can fix numerical issues
        return True
    except la.LinAlgError:
        return False





def GetEigenVal(element_size,Lx,Ly,step,n,p,bx_list,by_list,rcoeff_matrix): # Runs only once

    ########################## Generate stochastic element grid ##############################
    tol = 0.000001
    xc =  np.arange(element_size/2,Lx-element_size/2+tol,element_size/step)  
    yc =  np.arange(element_size/2,Ly-element_size/2+tol,element_size/step)  
    nw1 = len(xc)
    nw2 = len(yc)
    nw = nw1*nw2
    xc_array,yc_array= np.meshgrid(xc,yc)         # stochastic element grid points # CHECK IF MATCH WITH ELEMENT CENTROIDS
    x,y= xc_array.flatten(),yc_array.flatten()
    # calculate all distances in matrix form, statistically anisotropic field #
    dist_x = np.zeros((nw,nw))
    dist_y = np.zeros((nw,nw))
    for i in range(0,nw):
        for j in range(0,nw):
            dist_x[i,j]= abs(x[i]-x[j])
            dist_y[i,j]= abs(y[i]-y[j])
    ############################### Compute correlation matrix based on distance matrix ############################################

    Corr_matrix = np.zeros((nw*n,nw*n))
    for i in range (0,n):
        Corr_matrix[(i)*nw:(i+1)*nw, (i)*nw:(i+1)*nw] = expMattern(p,dist_x,dist_y,bx_list[i],by_list[i])  # autocorrelation 
        j = i 
        #print(Corr_matrix[(i)*nw:(i+1)*nw, (i)*nw:(i+1)*nw])
        while j < n: 
            Corr_matrix[(i)*nw:(i+1)*nw, (j)*nw:(j+1)*nw] = rcoeff_matrix[i,j]*expMattern(p,dist_x,dist_y,(bx_list[i] + bx_list[j])/2,(by_list[i] + by_list[j]/2)) # cross-correlations
            j= j+1
    Corr_matrix = np.triu(Corr_matrix) # Computed matrix is only upper triangular, make it symmetric
    Corr_matrix = (Corr_matrix +np.transpose(Corr_matrix)) - np.eye(len(Corr_matrix))* Corr_matrix
    Corr_matrix = nearestPD(Corr_matrix) # convert to pos def in case it is not 

    ############### Compute eigenvalues of correlation matrix and generate RF #################
    eigenvalue_vector, eigenvector_matrix = la.eig(Corr_matrix)    
    return [eigenvalue_vector,eigenvector_matrix]