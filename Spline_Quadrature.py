################################################################################
import numpy as np
import scipy as sp
import splipy as spl
import sys
################################################################################
import time


def Prepare_Data(T,p):
    """Prepare_Data, for creating the initial condition and augmenting the knot
    vector. It must also return n and the constant integral vector in F.
    """
    basis = spl.BSplineBasis(p+1, knots=T)
    
    #number of weights of knots
    n = len(basis)-p-1
    
    #augmenting knot vector
    if n%2==1:
        t = (basis.knots[p]+basis.knots[p+1])/2
        basis.insert_knot(t) #(τp+1 + τp+2)/2 between τp+1 and τp+2
        n+=1
        
    #constant integral vector in F.
    I = (basis.knots[p+1:]-basis.knots[:-p-1])/(p+1)
    
    #create initial condition
    G = np.array(basis.greville())
    X = np.array(G[1::2]+G[::2])/2
    W = I[1::2]+I[::2]
    
    return basis, I, W, X, n

def Assembly(basis,I,W,X,n):
    '''updating Fn and ∂Fn every time in the Newton iteration.
    Recall that ∂Fn must be permuted and made sparse.
    '''
    N = basis.evaluate(X).T
    dN = basis.evaluate(X, d=1).T

    F = np.zeros(n)
    F[:] = np.dot(N, W) - I
    
    
    dFde = dN*sp.sparse.diags(W)
    
    J = np.concatenate((N, dFde), axis=1) #dN*diag(sparse(W))
        
    #Uncomment if we want sparse
#    J = sp.sparse.csr_matrix(J)    

    return F, J

basis, I, W, X, n = Prepare_Data(t1, 2)


def Spline_Quadrature(T, p):
    """Newton iteration"""
    
    #Tolerance
    TOL = 1e-11

    #Prepare_data
    basis, I, W, X, n = Prepare_Data(T, p)
    #print('Starting at: \n W: ', W, '\n X: ', X)
    
    #Assembly
    F, J = Assembly(basis, I, W, X, n)
    
    #Initialize improvement to inf and iteration count to 0
    norm = float('inf')
    itcount = 0
    
    while norm>TOL and itcount < 20:

        F, J = Assembly(basis,I,W,X,n)
        #Sparse solver
#        delta = sp.sparse.linalg.spsolve(J, F)
        #Dense solver        
        delta = np.linalg.solve(J, F)
        
        W -= delta[:int(n/2)]
        X -= delta[int(n/2):]

        norm = np.linalg.norm(delta)

        itcount+=1
        
        #print('Iteration ', itcount, '\t X = ',X ,'\t W = ', W, '\t Norm = %0.2E' % norm)
        
        #If goes outside boundaries, we have a singular matrix
        if min(X)<T[0] or max(X)>T[-1]:
            print('SINGULAR MATRIX!')
            break

     #Does not converge within 20 iterations
    if abs(norm) > TOL:
        itcount = -1
        
    return W, X, itcount


t1 = np.array([0, 0, 0, 1, 2, 3, 4, 4, 4])
t2 = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4])
t3 = np.array([0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4])
t4 = np.array([0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4])

t = time.clock()
Spline_Quadrature(t1, 2)
Spline_Quadrature(t2, 2)
Spline_Quadrature(t3, 3)
Spline_Quadrature(t4, 3)
print("Time spent: {}".format(int((time.clock() - t)*1000)), "ms.", sep="")