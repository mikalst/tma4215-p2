################################################################################
import sys
import Spline_Quadrature as SPQ
from splipy.io import G2
import numpy as np
################################################################################


def Curve_Integral(path):
    with G2(path) as file:
        #Make a splipy curve object from file
        curve = file.read()[0]
        
        #Extract knots and order
        Tau = curve.knots(0, True)
        p = curve.order(0) - 1
        
        #Calculate optimal quadrature points and weights
        W, X, itcount = SPQ.Spline_Quadrature(Tau, p, False)
        
#        #This is equivalent of the one-liner on line 33, but line 33 is much
#        #faster and so we have used the latter 
#        path = 0
#        for i, x in enumerate(X):
#            dd = curve.derivative(x)
#            path += W[i]*np.linalg.norm(dd)
        
        
        #curve.derivative(X) evaluates the derivative in the points X. As this
        # returns a (n*2) array in the format of n*(dx/dt, dy/dt) we use 
        # np.linalg.norm with axis=1 to take the 2-norm s.t we get 
        # ds = sqrt((dx/dt)**2 + (dy/dt)**2). We then dot the weights with the
        # ds-vector. 
        path = np.dot(W, np.linalg.norm(curve.derivative(X), axis = 1))


        return path


def main(*args, **kwargs):
    
    args = args if args else sys.argv[1:]
    
    if not args:
        path = "Curve.g2"
    else:    
        path = args[0]
    
    print("Numpy.dot:", Curve_Integral(path))
    
    
    
if __name__ == "__main__":
    main()