import xml.etree.ElementTree as et

import numpy as np
from numpy import sin, cos, exp  # Need this for the xml files


def XML_Extraction(XMLFILE):
    tree = et.parse(XMLFILE)
    root = tree.getroot()

    f = eval('lambda x: ' + root[0].text)
    exactIntegral = eval(root[1].text)

    return f, exactIntegral


def Gauss_Legendre_Data(n: int):

    G = np.zeros((n, 2), dtype=np.float128)
    
    if n%2==0:
        #Polynomial is even.
        #Due to symmetry of the zeroes, we only need to solve on x in [-1, 0].
        for k in range(1, int(n/2)+1): # 1 in {1, 2, ... n/2 + 1}
            xl = np.cos((2*k-1)*np.pi/(2*n+1))
            xh = np.cos(2*k*np.pi/(2*n+1))
            
            x0 = np.float128(1/2*(xl+ xh))
            x = Olver(n, x0)
            
            l0, l0n = Legendre_0(n, x)
            l1 = Legendre_1(n, x, l0)
            w = 2/((1-x**2)*l1**2)
            
            G[k-1, 0] = x
            G[k-1, 1] = w
            
            G[n-k, 0] = -x
            G[n-k, 1] = w
            
    else:
        #Polynomial is odd.
        #Same as for even, but we are guaranteed to have a zero at x = 0.
        l0, l0n = Legendre_0(n, 0)
        l1 = Legendre_1(n, 0, l0)
        w = 2/(l1**2)
        G[int((n-1)/2), 0] = 0
        G[int((n-1)/2), 1] = w
        
        for k in range(1, int((n+1)/2)): # 1 in {1, 2, ... (n+1)/2}
            xl = np.cos((2*k-1)*np.pi/(2*n+1))
            xh = np.cos(2*k*np.pi/(2*n+1))
            
            x0 = np.float128(1/2*(xl+ xh))
            x = Olver(n, x0)
            
            l0, l0n = Legendre_0(n, x)
            l1 = Legendre_1(n, x, l0)
            w = 2/((1-x**2)*l1**2)
            
            G[k-1, 0] = x
            G[k-1, 1] = w
            
            G[n-k, 0] = -x
            G[n-k, 1] = w


    return G


def Olver(n,x0):

    TOL = 1e-14
    x = x0
    s = TOL+1

    while (abs(s) > TOL):
        l0, l0n = Legendre_0(n, x)
        l1n = Legendre_1(n, x, l0)
        l2n = Legendre_2(n, x, l0)

        s = l0n/l1n - l2n*l0n**2/(2*l1n**3)

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
        l0n = ((2*i+1)*x*l0[i] - i*l0[i-1])/(i+1)# Now l0[i] and l0[i-1] are float128
        l0[i+1] = l0n

    return l0, l0n


def Legendre_1(n, x, l0):
    if n == 0:
        return 0
    elif n == 1:
        return 1

    if n % 2 == 0:
        indices = np.arange(1, n, step=2)
    else:
        indices = np.arange(0, n, step=2)

    return np.sum((2 * indices + 1) * l0[indices].T)


def Legendre_2(n, x, l0):
    if n == 0:
        return 0
    elif n == 1:
        return 0
    elif n == 2:
        return 3

    indices = np.arange(n % 2, n-1, step=2)
    return np.sum((indices+1/2)*(n*(n+1)-indices*(indices+1))*l0[indices].T)


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
    print(Return_Quadrature("functions/f6.xml", 100))
