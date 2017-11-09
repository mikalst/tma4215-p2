################################################################################
import xml.etree.ElementTree as et
import numpy as np
from numpy import sin, cos, exp # Need this for the xml files
import time
################################################################################


def XML_Extraction(XMLFILE):
    tree = et.parse(XMLFILE)
    root = tree.getroot()
    
    f = lambda x: eval(root[0].text)
    exactIntegral = eval(root[1].text)
    
    return f, exactIntegral
    
def Gauss_Legendre_Data(n):
    
    G = np.zeros((n, 2))
    
    for k in range(1, n+1): # 1 in {1, 2, ... n}
        xl = np.cos((2*k-1)*np.pi/(2*n+1))
        xh = np.cos(2*k*np.pi/(2*n+1))
        x0 = 1/2*(xl+ xh)
        x = Olver(n, x0)
        
        #print("{} -> {}".format((xl, xh), x))
        
        l1 = Legendre_1(n, x)[1][-1]
        w = 2/((1-x**2)*l1**2)
        
        G[k-1, 0] = x
        G[k-1, 1] = w

    
    return G
    

def Olver(n,x0):
    
    TOL = 1e-14
    x = x0
    s = TOL+1

    while (abs(s) >= TOL):   
        l0 = Legendre_0(n, x)[1][-1]
        l1 = Legendre_1(n, x)[1][-1]
        l2 = Legendre_2(n, x)[1][-1]

        s = np.float128(l0/l1 - l2*l0**2/(2*l1**3))
        
        x = np.float128(x - s)
        
    return x
    

def Legendre_0(n,x):
    if n == 0:
        return 1, [1]
    elif n == 1:
        return x, [x]
    
    l0 = [None for _ in range(n+1)]
    l0[0] = 1
    l0[1] = x
    
    for i in range(1, n):
        a = ((2*i+1)*x*l0[i] - i*l0[i-1])/(i+1)
        l0[i+1] = a
        
    #print("x = {} \nl0 = {}".format(x, l0))
        
    return sum(l0), l0


def Legendre_1(n,x):
    if n == 0:
        return 0, [0]
    elif n == 1:
        return 1, [1]
    
    t, l0 = Legendre_0(n, x)

    l1 = [0 for _ in range(n+1)]    
    l1[0] = 0
    l1[1] = 1
    
    for i in range(2, n+1):
        l1[i] = (2*i-1)/i*l0[i-1] + (2*i-1)/i*x*l1[i-1] - (i-1)/i * l1[i-2]
        
    #print("x = {} \nl0 = {}\n l1 = {}".format(x, l0, l1))
    #print("l1 = {}".format(l1[-1]))
    
    return sum(l1), l1



def Legendre_2(n,x):
    if n == 0:
        return 0, [0]
    elif n == 1:
        return 0, [0]
    elif n == 2:
        return 3, [3]
    
    t, l1 = Legendre_1(n, x)
    l2 = [0 for _ in range(n+1)]
    
    l2[0] = 0
    l2[1] = 0
    l2[2] = 3
    
    for i in range(3, n+1):
        l2[i] = 2*(2*i-1)/i*l1[i-1] + (2*i-1)/i*x*l2[i-1] - (i-1)/i*l2[i-2]
        
    #print("x = {} \nl1 = {}\n l2 = {}".format(x, l1, l2))
    
    return sum(l2), l2
    

def Gauss_Legendre_Quadrature(n,G,f):
    total = 0
    for (x, w) in G:
        total += w*f(x)
    
    return total

def Return_Quadrature(XMLFILE,n):
    G = Gauss_Legendre_Data(n)
    f, exactIntegral = XML_Extraction(XMLFILE)
    
    numericalIntegral = Gauss_Legendre_Quadrature(n, G, f)
    
    relError = np.abs((numericalIntegral - exactIntegral)/exactIntegral)
    
#    print("numerical = {}, analytic = {}, difference = {}"
#          .format(numericalIntegral, exactIntegral, (numericalIntegral - exactIntegral)))    

    return numericalIntegral, exactIntegral, relError
        
    
if __name__ == "__main__":
    Return_Quadrature("functions/f6.xml", 1)

    