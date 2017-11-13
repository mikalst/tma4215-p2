################################################################################
import numpy as np
import splipy as spl
import matplotlib.pyplot as plt
################################################################################


def Basis_Plot():
    Tau1 = [0, 0, 0, 1, 2, 3, 4, 4, 4]# p = 9-4-1 = 5
    Tau2 = [0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4]# p = 12-4-1 = 7
    Tau3 = [0, 0, 0, 0, 1, 2, 3, 4, 4, 4, 4]# p = 11 - 1 - 4
    Tau4 = [0, 0, 0, 0, 1, 1, 2, 2, 3, 3, 4, 4, 4, 4]
    
    a = spl.BSplineBasis(3, knots=Tau1)
    a.num_functions()
        
    X = np.linspace(0, 4, 200)
    plt.plot(X, a(X))
    plt.savefig("2a_T1.pdf")
    plt.show()
    
    a = spl.BSplineBasis(3, knots=Tau2)    
    X = np.linspace(0, 4, 200)
    plt.plot(X, a(X))
    plt.savefig("2a_T2.pdf")
    plt.show()
    
    a = spl.BSplineBasis(4, knots=Tau3)    
    X = np.linspace(0, 4, 200)
    plt.plot(X, a(X))
    plt.savefig("2a_T3.pdf")
    plt.show()
    
    a = spl.BSplineBasis(4, knots=Tau4)    
    X = np.linspace(0, 4, 200)
    plt.plot(X, a(X))
    plt.savefig("2a_T4.pdf")
    plt.show()
    
    
if __name__ == "__main__":
    Basis_Plot()
    

