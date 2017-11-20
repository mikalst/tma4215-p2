################################################################################
import sys
import numpy as np
import Spline_Quadrature as SPQ
from splipy.io import G2
################################################################################


def Area_Integral(path):
    with G2(path) as file:
        surface = file.read()[0]
        
        TauX = surface.knots(0, True)
        TauY = surface.knots(1, True)
        
        px = surface.order(0) - 1
        py = surface.order(1) - 1
        
        Wx, Xx, itCountx = SPQ.Spline_Quadrature(TauX, px, False)
        Wy, Xy, itCounty = SPQ.Spline_Quadrature(TauY, py, False)
                
        area = 0
        for i in range(len(Wx)):
            for j in range(len(Wy)):

                Dx = surface.derivative(Xx[i], Xy[j], d = (1,0))
                Dy = surface.derivative(Xx[i], Xy[j], d = (0,1))
                
                J = np.abs(np.cross(Dx, Dy))
                
                area += Wx[i] * Wy[j] * J

        return area

    

def main(*args, **kwargs):
    
    args = args if args else sys.argv[1:]
    
    if not args:
        path = "Area.g2"
    else:    
        path = args[0]
    
    print(Area_Integral(path))
    
    
    
if __name__ == "__main__":
    main()