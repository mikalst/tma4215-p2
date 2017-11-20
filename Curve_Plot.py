################################################################################
import numpy as np
#import splipy as spl
from splipy.io import G2
import matplotlib.pyplot as plt
#import sys
################################################################################


def Curve_Plot():
    with G2('Curve.g2') as file:
        curve = (file.read()[0])
        
        N = 180
        t = np.linspace(curve.start(), curve.end(), N)        
        
        print("Start: {}   End: {}".format(curve.start(), curve.end()))
        
        V = curve.evaluate(t)
        X = V[:, 0]
        Y = V[:, 1]
        
        plt.figure()
        plt.gcf().subplots_adjust(left=0.18)
        plt.plot(X, Y)
        plt.xlabel("$x$")
        plt.ylabel("$y$")#, rotation='horizontal')
        plt.savefig("Curve_Plot.pdf")
        plt.show()
        

if __name__ == "__main__":    
    Curve_Plot()