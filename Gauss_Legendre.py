################################################################################
import xml.etree.ElementTree as et
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
    
    for i in range(1, n+1):
        x0 = 1/2*(np.cos((2*i-1)*np.pi/(2*n+1)) + np.cos((2*i*np.pi)/(2*n+1)))
        X[i-1] = Olver(n, x0)
        l1 = Legendre_1(n, X[i-1])
        W[i-1] = 2/((1-X[i-1]**2)*l1**2)
    return X, W
    

def Olver(n,x0):
    
    TOL = 1e-2
    x = x0
    s = TOL+1

    while (abs(s) >= TOL):   
        l0, t = Legendre_0(n, x)
        l1 = Legendre_1(n, x)
        l2 = Legendre_2(n, x)

                 
        s = np.float128(l0/l1 + l2*l0**2/(2*l1**3))
        x = np.float128(x - s)
        
        
    return x
    

def Legendre_0(n,x):
    
    if n == 1:
        return 1
    if n == 2:
        return x
    
    res = [None for _ in range(n)]
    res[0] = 1
    res[1] = x
    
    for i in range(1, n-1):
        res[i+1] = ((2*i+1)*x*res[i] - i*res[i-1])/(i+1)
        
    return res[-1], res


def Legendre_1(n,x):
    t, res = Legendre_0(n, x)
    indices = filter(lambda k: (k+n)%2==1, range(0, n-1))
    
    l1 = [(2*k +1)*res[k] for k in indices]
    
    return sum(l1)


def Legendre_2(n,x):
    t, res = Legendre_0(n, x)
    indices = filter(lambda k: (k+n)%2==0, range(0, n-2))
    
    l2 = [(k+1/2)*(n*(n+1) - k*(k+1))*res[k] for k in indices]
    
    return sum(l2)
    

def Gauss_Legendre_Quadrature(n,G,f):
    return

def Return_Quadrature(XMLFILE,n):
    return

    
if __name__ == "__main__":
#    
#    X = np.linspace(-0.95, 0.95, 100)
#    
#    Y = [Legendre_0(3, x)[0] for x in X]
#    dY = [Legendre_1(3, x) for x in X]
#    ddY = [Legendre_2(3, x) for x in X]
#    
#    plt.plot(X, Y)
#    plt.plot(X, dY)
#    plt.plot(X, ddY)
#    

    
    
    
    print(Gauss_Legendre_Data(3))
        
    #print(np.polynomial.legendre.leggauss(5))
    
    
