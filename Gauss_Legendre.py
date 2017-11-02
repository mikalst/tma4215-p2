################################################################################
import xml.etree.ElementTree as et
import numpy as np
import matplotlib.pyplot as plt
################################################################################


def XML_Extraction(XMLFILE):
    tree = et.parse(XMLFILE)
    root = tree.getroot()
    
    f = lambda x: eval(root[0].text)
    exactIntegral = eval(root[1].text)
    
    return f, exactIntegral
    
def Gauss_Legendre_Data(n):
    
    X = []
    W = []
    
    for k in range(1, n+1): # 1 in {1, 2, ... n}
        xl = np.cos((2*k-1)*np.pi/(2*n+1))
        xh = np.cos(2*k*np.pi/(2*n+1))
        x0 = 1/2*(xl+ xh)
        x = Olver(n, x0)
        print("{} -> {}".format((xl, xh), x))
        X.append(x)
        if (x):    
            l1 = Legendre_1(n, x)[1][-1]
            W.append(2/((1-x**2)*l1**2))
    
    G = np.stack((X, W), axis = 1)
    
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
        return 1
    elif n == 1:
        return x
    
    l0 = [None for _ in range(n+1)]
    l0[0] = 1
    l0[1] = x
    
    for i in range(1, n):
        a = ((2*i+1)*x*l0[i] - i*l0[i-1])/(i+1)
        l0[i+1] = a
        
    #print("x = {} \nl0 = {}".format(x, l0))
        
    return sum(l0), l0


def Legendre_1(n,x):
    
    t, l0 = Legendre_0(n, x)

    l1 = [0 for _ in range(n+1)]    
    l1[0] = 0
    l1[1] = 1
    
    for i in range(2, n+1):
        l1[i] = (2*i-1)/i*l0[i-1] + (2*i-1)/i*x*l1[i-1] - (i-1)/i * l1[i-2]
        
    #print("x = {} \nl0 = {}\n l1 = {}".format(x, l0, l1))
    
    return sum(l1), l1



def Legendre_2(n,x):
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
    
    print("numerical = {}, analytic = {}, difference = {}"
          .format(numericalIntegral, exactIntegral, (numericalIntegral - exactIntegral)))
    
    return numericalIntegral


def TestLegendreSeries():

    N = 8
    
    X = np.linspace(-1.0, 1.0, 50)
    
    
    a = np.polynomial.legendre.Legendre([1]*(N+1))
    b = a.deriv()
    c = a.deriv(2)
    
    Y1 = [Legendre_0(N, x)[0] for x in X]
    Y2 = [a(x) for x in X]
    #print(b(-1), Legendre_1(N, -1))
    
    
    dY1 = [Legendre_1(N, x)[0] for x in X]
    dY2 = [b(x) for x in X]
    
    ddY1 = [Legendre_2(N, x)[0] for x in X]
    ddY2 = [c(x) for x in X]
    
    
    plt.style.use("ggplot")
    plt.plot(X, Y1)
    plt.plot(X, Y2)
    plt.legend(["Y Ours","Y Numpy"])
    plt.show()
    
    plt.plot(X, dY1)
    plt.plot(X, dY2)
    plt.legend(["dY Ours","dY Numpy"])
    plt.show()
    
    plt.plot(X, ddY1)
    plt.plot(X, ddY2)
    plt.legend(["ddY Ours","ddY Numpy"])
    plt.show()
        
    
if __name__ == "__main__":
#    
#    N = 5
#    
#    X = np.linspace(-1.0, 1.0, 100)
#    Y = [Legendre_0(N, x)[1][-1] for x in X]
#    dY = [Legendre_1(N, x)[1][-1] for x in X]
#    ddY = [Legendre_2(N, x)[1][-1] for x in X]
#    
#    plt.plot(X, Y)
#    plt.show()
#    plt.plot(X, dY)
#    plt.show()
#    plt.plot(X, ddY)
#    
#    X = np.linspace(-1.0, 1.0, 100)
#    Y = [Legendre_0(20, x)[1][-1] for x in X]
#    plt.plot(X, Y)
#    G = Gauss_Legendre_Data(20)
#    for el in G:
#        print(el[0], el[1])
#        plt.bar(el[0], el[1], width=0.1, color="black")

    Return_Quadrature("functions/f1.xml", 6)

    

    #print(np.polynomial.legendre.leggauss(5))
    