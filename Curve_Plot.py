################################################################################
import numpy as np
#import splipy as spl
from splipy.io import *
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
#import sys
################################################################################


def Curve_Plot():
    with G2('Curve.g2') as file:
        curve = file.read()[0]
        
        N = 1000
        
        T = np.linspace(curve.start(), curve.end(), N)        

        V = curve.evaluate(T)
        X = V[:, 0]
        Y = V[:, 1]

        # set up a list of (x,y) points
        points = np.array([X,Y]).transpose().reshape(-1,1,2)
        
        # set up a list of segments
        segs = np.concatenate([points[:-1],points[1:]],axis=1)

        # make the collection of segments
        lc = LineCollection(segs, cmap=plt.get_cmap('jet'))
        lc.set_array(T) # color the segments by our parameter
        
        # plot the collection
        plt.gca().add_collection(lc)
        plt.xlim(X.min(), X.max())
        plt.ylim(Y.min(), Y.max())
        cb = plt.colorbar(lc, label="$t$")
        
        plt.xlabel("$x$")
        plt.ylabel("$y$")

        plt.show()
        
Curve_Plot()