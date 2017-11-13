################################################################################
import numpy as np
import scipy as sp
import splipy as spl
import sys
################################################################################

#
def Prepare_Data(T, p):
                
    basis = spl.BSplineBasis(p, T)
    
    #number of weights of knots
    n = basis.num_functions()
    
    print(n)
    
    #Add extra knot to make even number of knots
    if n%2==1:
        t = (basis.knots[p]+basis.knots[p+1])/2
        basis.insert_knot(t)
        n = basis.num_functions()
        
    T = basis.knots
        
    #Convert n to format specified in [1].
    n = int(n/2)
    
    G = np.array(basis.greville())
    I = (T[p:] - T[:-(p)])/(p+1) 

    X = np.empty((n, ))
    W = np.empty((n, ))
    
#    print(G)
#    print(I)

    for i in range(n):
        X[i] = (G[2*i+1] + G[2*i])/2
        W[i] = I[2*i+1] + I[2*i]
        
#    print(X)
#    print(W)
    
    return basis, I, W, X, n
        
    
def Assembly(basis,I,W,X,n):
    
    N = basis.evaluate(X)
    print("N = ",N ,sep = "\n")
    dN = basis.evaluate(X, d=1)
    print("dN = ", dN, sep = "\n")
    #print(X, N, dN, W,sep="\n\n")
    
    F = np.empty((2*n, ))
    
    J = np.empty((2*n, 2*n))
    
    for j in range(2*n):
        for i in range(n): 
            #print(i, j)
            J[i, j] = N[i, j]
            J[i+n, j] = np.dot(W[i], dN[i, j])
        
        F[j] = np.dot(W, N[:, j]) - I[j]
    
    print("F = ", F, sep="\n")
        
    print("J = ", J, sep="\n")
    
    #Run Newton
    S = np.linalg.solve(J, F)
    
    W = W - S[:n]
    X = X - S[n:]
    
    print("X = {}".format(X))
    
    return X
    
    

def Spline_Quadrature():
    T = [0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4]
    p = 3
    
    Errors = []
    TOL = 1e-10
    improvement = float('inf')
    basis, I, W, X, n = Prepare_Data(T, p)
    print("X0 = {}".format(X))
    
    while improvement > TOL:
        print("---")
        X_last = X
        X = Assembly(basis, I, W, X, n)
        Errors.append(X)
        improvement = np.linalg.norm(X-X_last)
    
    print(Errors)


t1 = np.array([0, 0, 0, 1, 2, 3, 4, 4, 4])
Spline_Quadrature()

    
    
