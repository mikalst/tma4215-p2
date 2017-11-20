################################################################################
import numpy as np
import splipy as spl
from splipy.io import G2
import splipy.surface_factory as spf
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
################################################################################


def Area_Plot():
    with G2("Area.g2") as file:
        surface = file.read()[0]

        #2D        
        N = 30
        u = np.linspace(surface.start('u'), surface.end('u'), N)
        v = np.linspace(surface.start('v'), surface.end('v'), N)
        X = surface(u, v)
        
        plt.figure()
        plt.plot(X[:,:,0], X[:,:,1], '-', color="red")
#        plt.plot(X[:,:,0], X[:,:,1].T, '-', color="blue")
                
        plt.gcf().subplots_adjust(bottom=0.18)
        plt.xticks(rotation=45)
        plt.axis('equal')
        plt.savefig("Area_Plot.pdf")
        plt.show()
        
        

        

if __name__ == "__main__":
    Area_Plot()
        