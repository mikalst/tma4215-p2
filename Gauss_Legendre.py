import xml.etree.ElementTree as et

import numpy as np
from numpy import sin, cos, exp  # Need this for the xml files


def XML_Extraction(XMLFILE):
    tree = et.parse(XMLFILE)
    root = tree.getroot()

    f = eval('lambda x: ' + root[0].text)
    exactIntegral = eval(root[1].text)

    return f, exactIntegral


def Gauss_Legendre_Data(n):

    G = np.zeros((n, 2), dtype=np.float128)

    for k in range(1, n+1): # 1 in {1, 2, ... n}

        xl = np.cos((2*k-1)*np.pi/(2*n+1))
        xh = np.cos(2*k*np.pi/(2*n+1))

        x0 = np.float128(1/2*(xl+ xh))
        x = Olver(n, x0)

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

    l1 = np.zeros((n+1, 1), dtype=np.float128)
    indices = filter(lambda k: (k+n)%2==1, range(0, n))

    for k in indices:
        l1[k] = (2*k+1)*l0[k]

    l1n = np.sum(l1)

    return l1n


def Legendre_2(n, x, l0):
    if n == 0:
        return 0
    elif n == 1:
        return 0
    elif n == 2:
        return 3

    l2 = np.zeros((n+1, 1), dtype=np.float128)
    indices = filter(lambda k: (k+n)%2==0, range(0, n-1))

    for k in indices:
        l2[k] = (k+1/2)*(n*(n+1)-k*(k+1))*l0[k]

    l2n = np.sum(l2)

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
    print(Return_Quadrature("functions/f6.xml", 100))
