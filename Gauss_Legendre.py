################################################################################
import xml.etree.ElementTree as et
import numpy as np
from numpy import sin, cos, exp # Need this for the xml files
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
        l0, l0n = Legendre_0(n, x)
        l1 = Legendre_1(n, x, l0)
        w = 2/((1-x**2)*l1**2)
        
        G[k-1, 0] = x
        G[k-1, 1] = w

    
    return G
    

def Olver(n,x0):
    
    TOL = 1e-14
    x = x0
    s = TOL+1

    while (abs(s) > TOL):   
        l0, l0n = Legendre_0(n, x)
        l1n = Legendre_1(n, x, l0)
        l2n = Legendre_2(n, x, l0)
        #print(l0n, l1n, l2n)
        s = l0n/l1n - l2n*l0n**2/(2*l1n**3)
        #print(l0s)
        
        x = x - s
        
    return x
    

def Legendre_0(n,x):
    if n == 0:
        return [1], 1
    elif n == 1:
        return [x], x
    
    l0 = np.zeros((n+1, 1), dtype=np.float128)
    l0[0] = 1
    l0[1] = x
    
    for i in range(1, n):
        a = ((2*i+1)*x*l0[i] - i*l0[i-1])/(i+1)
        l0[i+1] = a
        
    #print("x = {} \nl0 = {}".format(x, l0))
    #print(l0[-1])
        
    return l0, l0[-1]


def Legendre_1(n, x, l0):
    if n == 0:
        return 0
    elif n == 1:
        return 1
    
    indices = filter(lambda k: (k+n)%2==1, range(0, n))
    l1n = sum([(2*k+1)*l0[k] for k in indices])
    
    return l1n


def Legendre_2(n, x, l0):
    if n == 0:
        return 0
    elif n == 1:
        return 0
    elif n == 2:
        return 3
        
    indices = filter(lambda k: (k+n)%2==0, range(0, n-1))
    l2n = sum([(k+1/2)*(n*(n+1)-k*(k+1))*l0[k] for k in indices])
    
    return l2n
    

def Gauss_Legendre_Quadrature(n,G,f):
    total = 0
    for (x, w) in G:
        total += w*f(x)
    
    return total

def Return_Quadrature(XMLFILE,n):

    G = Gauss_Legendre_Data(n)
    f, exactIntegral = XML_Extraction(XMLFILE)
    numInt = Gauss_Legendre_Quadrature(n, G, f)
        
    relErr = np.abs((numInt - exactIntegral)/exactIntegral)
    
    return numInt, exactIntegral, relErr

        
    
if __name__ == "__main__":
    print(Return_Quadrature("functions/f6.xml", 10))


    