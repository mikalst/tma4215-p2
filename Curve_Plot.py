################################################################################
import numpy as np
#import splipy as spl
from splipy.io import G2
import matplotlib.pyplot as plt
from matplotlib.collections import LineCollection
#import sys
################################################################################


def Curve_Plot():
    with G2('Curve.g2') as file:
        curve = file.read()[0]
        
        N = 150
        t = np.linspace(curve.start(), curve.end(), N)        
        V = curve.evaluate(t)
        X = V[:, 0]
        Y = V[:, 1]

        # set up a list of (x,y) points
        points = np.array([X,Y]).transpose().reshape(-1,1,2)
        
        # set up a list of segments
        segs = np.concatenate([points[:-1],points[1:]],axis=1)

        # make the collection of segments
        lc = LineCollection(segs, cmap=plt.get_cmap('jet'))
        lc.set_array(t) # color the segments by our parameter
        
        # plot the collection
        plt.gca().add_collection(lc)
        plt.xlim(X.min(), X.max())
        plt.ylim(Y.min(), Y.max())
        plt.colorbar(lc, label="$t$")
        
        plt.axis('equal')
        
        plt.xlabel("$x$")
        plt.ylabel("$y$")

        plt.show()
        
Curve_Plot()