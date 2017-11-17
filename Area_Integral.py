################################################################################
import sys
import numpy as np
import scipy as sp

import Spline_Quadrature as SPQ
from splipy.io import G2
import splipy.surface_factory as spf
################################################################################





def Area_Integral(path):
    with G2(path) as file:
        pass

        return

    

def main(*args, **kwargs):
    
    args = args if args else sys.argv[1:]
    
    if not args:
        path = "Curve.g2"
    else:    
        path = args[0]
    
    print(Area_Integral(path))
    
    
    
if __name__ == "__main__":
    main()