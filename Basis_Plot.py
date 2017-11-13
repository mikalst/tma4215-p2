################################################################################
import numpy as np
import splipy as spl
import matplotlib.pyplot as plt
import splipy
################################################################################


def Basis_Plot():
    Tau1 = [0, 0, 0, 1, 2, 3, 4, 4, 4]
    Tau2 = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]
    Tau3 = [0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4]
    Tau4 = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4]
    
    a = splipy.BSplineBasis(2, knots=Tau1)
    a.num_functions()
        
    X = np.linspace(0, 4, 200)
    plt.plot(X, a(X))
    plt.show()
    
    a = splipy.BSplineBasis(2, knots=Tau2)    
    X = np.linspace(0, 4, 200)
    plt.plot(X, a(X))
    plt.show()
    
    a = splipy.BSplineBasis(2, knots=Tau3)    
    X = np.linspace(0, 4, 200)
    plt.plot(X, a(X))
    plt.show()
    
    a = splipy.BSplineBasis(2, knots=Tau4)    
    X = np.linspace(0, 4, 200)
    plt.plot(X, a(X))
    plt.show()
    
    
    
if __name__ == "__main__":
    Basis_Plot()
    

