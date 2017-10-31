################################################################################
import xml.etree.ElementTree as et
import numpy as np
from numpy.polynomial.legendre import *
import matplotlib.pyplot as plt
################################################################################


def XML_Extraction(XMLFILE):
    tree = et.parse(XMLFILE)
    root = tree.getroot()
    
    f = lambda x: eval(root[0].text)
    exactIntegral = float(root[1].text)
    
    return f, exactIntegral
    
def Gauss_Legendre_Data(n):
    
    X = [0 for _ in range(n)]
    W = [0 for _ in range(n)]
    
    for i in range(1, n):
        x0 = 1/2*(np.cos((2*i-1)*np.pi/(2*n+1)) + np.cos((2*i*np.pi)/(2*n+1)))
        X[i] = Olver(n, x0)
        l1 = Legendre_1(n, X[i])
        W[i] = 2/((1-X[i]**2)*l1**2)
    return X, W
    

def Olver(n,x0):
    
    TOL = 1e-2
    x = x0
    s = TOL+1

    while (abs(s) >= TOL):   
        l0 = Legendre_0(n, x)
        l1 = Legendre_1(n, x)
        l2 = Legendre_2(n, x)
                 
        s = np.float128(l0/l1 + l2*l0**2/(2*l1**3))
        x = np.float128(x - s)
        
        print(x)
        
    return x
    
#def evaluateLegendreSeries(coeff, x):
#    n = len(coeff)
#    res = [None for _ in range(n)]
#    res[0] = coeff[0]*1
#    res[1] = coeff[1]*x
#    
#    for i in range(1, n-1):
#        res[i+1] = coeff[i+1]*((2*i+1)*x*res[i] - i*res[i-1])/(i+1)
#        
#    return sum(res)
def evaluateLegendreSeries(coeff, x):
    a = Legendre(coeff)
    return a(x)


def Legendre_0(n,x):
    
    c = [1 for _ in range(n)]
            
    return evaluateLegendreSeries(c, x)

def Legendre_1(n,x):
    
    c = [1 for _ in range(n)]
    
    a = np.polynomial.legendre.Legendre(c)
    b = a.deriv()
    
    return b(x)


def Legendre_2(n,x):
    
    c = [1 for _ in range(n)]
    
    a = np.polynomial.legendre.Legendre(c)
    b = a.deriv(2)
    
    return b(x)

def Gauss_Legendre_Quadrature(n,G,f):
    return

def Return_Quadrature(XMLFILE,n):
    return

    
if __name__ == "__main__":
    
    l = Legendre((1, 1, 1))
    
    X = np.linspace(-1, 1, 100)
    Y = [l(x) for x in X]
    plt.plot(X, Y)
    
    
    #print(Gauss_Legendre_Data(5))
        
    #print(np.polynomial.legendre.leggauss(5))
    
    
