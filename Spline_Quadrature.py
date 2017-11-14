################################################################################
import numpy as np
import splipy as spl
import scipy as sp
################################################################################
import time


def Prepare_Data(T,p):
    """Prepare_Data, create the initial condition and augment the knot
    vector.
    """
    
    #Create basis
    basis = spl.BSplineBasis(p+1, knots=T)
    
    #As in the paper, 2*n = len(T) - p - 1, for ease of coding written here 
    #just as n. See below if (len(T) - p - 1) is odd. 
    n = len(T)-p-1
    
    #If n is odd, augment knot vector
    if n%2==1:

        t = (T[p]+T[p+1])/2
        basis.insert_knot(t)

        n+=1
        
    #constant integral vector in F.
    I = (basis.knots[p+1:]-basis.knots[:-p-1])/(p+1)
    
    #create initial condition
    G = np.array(basis.greville())
    X = (G[1::2]+G[::2])/2
    W = I[1::2]+I[::2]
    
    return basis, I, W, X, n

def Assembly(basis, I, W, X, n, dense = True):
    """update F and Jacobian every time in the Newton iteration.
    """
    
    N = basis.evaluate(X).T
    dN = basis.evaluate(X, d=1).T

    F = np.zeros(n)
    F[:] = np.dot(N, W) - I
    
    #I have here not utilized the band structure of J
    #Doing this might improve running time.
    
    #Dense solver
    #------------
    if dense:
        dFde = dN*np.diag(W)
        J = np.concatenate((N, dFde), axis=1) #dN*diag(sparse(W))
    
    #Sparse solver
    #------------
    else:
        dFde = dN*sp.sparse.diags(W)
        J = np.concatenate((N, dFde), axis=1) #dN*diag(sparse(W))
        J = sp.sparse.csr_matrix(J)    

    return J, F


def Spline_Quadrature(T, p, dense = True):
    """Newton iteration"""
    
    #Tolerance
    TOL = 1e-11

    #Prepare_data
    basis, I, W, X, n = Prepare_Data(T, p)

    #print('Starting at: \n W: ', W, '\n X: ', X)
        
    #Initialize improvement to inf and iteration count to 0
    norm = float('inf')
    itcount = 0
    
    while norm>TOL and itcount < 20:

        J, F = Assembly(basis,I ,W ,X ,n , dense)

        #Dense solver        
        if dense:
            delta = np.linalg.solve(J, F)
        
        #Sparse solver
        else:
            delta = sp.sparse.linalg.spsolve(J, F)

        
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


def testRunningTime(dense = True):
    
    t1 = np.array([0, 0, 0, 1, 2, 3, 4, 4, 4])
    t2 = np.array([0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4])
    t3 = np.array([0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4])
    t4 = np.array([0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4])
#    t5 = np.array([0, 0, 0, 0, 0, 1, 1, 1, 
    
    t = time.clock()
    for _ in range(100):
        Spline_Quadrature(t1, 2, dense)
        Spline_Quadrature(t2, 2, dense)
        Spline_Quadrature(t3, 3, dense)
        Spline_Quadrature(t4, 3, dense)
    print("Time spent: {}".format(int((time.clock() - t)*1000)), "ms.", sep="")
    
testRunningTime()