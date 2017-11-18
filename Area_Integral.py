################################################################################
import sys
import numpy as np
import scipy as sp

import Spline_Quadrature as SPQ
from splipy.io import G2
import splipy as spl
################################################################################





def Area_Integral(path):
    with G2(path) as file:
        surface = file.read()[0]
        
        TauX = surface.knots(0, True)
        TauY = surface.knots(1, True)

        px = surface.order(0) - 1
        px = surface.order(1) - 1
        
        Wx, Xx, itCountx = SPQ.Spline_Quadrature(TauX, px, False)
        Wy, Xy, itCounty = SPQ.Spline_Quadrature(TauY, px, False)
                
        area = 0
        
        for i, x in enumerate(Xx):
            for j, y in enumerate(Xy):
                dd = np.linalg.norm(surface.derivative(x, y))
                print("dd = ", dd)
                area += dd*Wx[i]*Wy[j]
                #print(area)
                
        print(area)
        
#        
#        dy = surface.derivative(Xy)
#        
#        #print(dx)
#        print(dy)
#        
#        J = np.cross(surface.derivative(Xx), surface.derivative(Xy))
#        
#        print(J.shape, J)
#        
#        print(W.shape, W)
#        
#        path = np.dot(W, np.cross(surface.derivative(Xx), surface.derivative(Xy)))
#        
#        print(path)
        
        
        
        
#        print(Tau, p, boundaries, sep="\n")
#        

        return

    

def main(*args, **kwargs):
    
    args = args if args else sys.argv[1:]
    
    if not args:
        path = "Area.g2"
    else:    
        path = args[0]
    
    print(Area_Integral(path))
    
    
    
if __name__ == "__main__":
    main()