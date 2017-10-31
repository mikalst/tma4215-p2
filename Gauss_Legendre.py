################################################################################
import xml.etree.ElementTree as et
import sys
import numpy as np
import matplotlib.pyplot as plt
################################################################################


def XML_Extraction(XMLFILE):
    tree = et.parse(XMLFILE)
    root = tree.getroot()
    
    f = lambda x: eval(root[0].text)
    exactIntegral = float(root[1].text)
    
    return f, exactIntegral
    
def Gauss_Legendre_Data(n):
    return
    

def Olver(n,x0):
    pass

def Legendre_0(n,x):
    
    if n == 1:
        return 1
    
    res = [0 for k in range(n)]
    res[0] = 1
    res[1] = x
    
    if n == 2:
        return sum(res)
    
    for i in range(1, n-1):
        res[i+1] = ( (2*i+1)*x*res[i] - i*res[i-1] )/(i+1)
        
    return sum(res)

def Legendre_1(n,x):
    return

def Legendre_2(n,x,L0,L1):
    return

def Gauss_Legendre_Quadrature(n,G,f):
    return

def Return_Quadrature(XMLFILE,n):
    return

    
if __name__ == "__main__":
    k = range(7, 10)
    x = np.linspace(-1, 1, 100)
    
    for j in k:    
        plt.plot(x, [Legendre_0(j, z) - np.polynomial.legendre.Legendre(([1]*j))(z) for z in x])
    
    
