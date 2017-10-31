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
    
    X = [0 for _ in range(n)]
    W = [0 for _ in range(n)]
    
    for i in range(1, n):
        x0 = 1/2*(np.cos((2*i-1)*np.pi/(2*n+1)) + np.cos((2*i*np.pi)/(2*n+1)))
        X[i] = Olver(n, x0)
        l1 = Legendre_1(n, X[i])
        W[i] = 2/((1-X[i]**2)*l1**2)
    return X, W
    

def Olver(n,x0):
    
    TOL = 1e-10
    x = x0
    s = TOL+1

    while (abs(s) >= TOL):   
        l0, t = Legendre_0(n, x)
        l1, t = Legendre_1(n, x)
#        l2 = Legendre_2(n, x)
                 
#        s = np.float128(l0/l1 + l2*l0**2/(2*l1**3))
        s = np.float128(l0/l1)
        x = np.float128(x - s)
        print(x)    
        
    return x
    

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
        
    return sum(res), res

def Legendre_1(n,x):
    total, res = Legendre_0(n, x)
    
    # Pick out the odd values from {0, 1, ... n-1}
    #indices = filter(lambda k: (k+n)%2==1, range(0, n-1))
    
    l1res = [0 for k in range(n)]
    
    l1res[0] = 0
    l1res[1] = 1
#    
#    for k in indices:
#        l1 += (2*k+1)*res[k]

    for i in range(1, n-2):
        l1res[i+1] = (n+1)/(x**2-1)*(x*res[i]-res[i-1]) 

    return sum(l1res), l1res

#def Legendre_2(n,x):
#    l1tot, l1res = Legendre_1(n, x)
#    
#    l2res = [0 for k in range(n)]
#    
#    l2res[0] = 0
#    l2res[1] = 0
##    
##    for k in indices:
##        l1 += (2*k+1)*res[k]
#
#    for i in range(1, n-2):
#        l1res[i+1] = (n+1)/(x**2-1)*(x*res[i]-res[i-1]) 
#
#    return sum(l1res), l1res
#
#    
#    return l2

def Gauss_Legendre_Quadrature(n,G,f):
    return

def Return_Quadrature(XMLFILE,n):
    return

    
if __name__ == "__main__":
    
    
    print(Gauss_Legendre_Data(5))
        
    #print(np.polynomial.legendre.leggauss(5))
    
    
