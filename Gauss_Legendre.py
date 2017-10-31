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
    
    X = []
    W = []
    
    for k in range(1, n+1):
        xhigh = -np.cos((2*k-1)*np.pi/(2*n+1))
        xlow = -np.cos((2*k)*np.pi/(2*n+1))
        x0 = 1/2*(xlow + xhigh)
        x = Olver(n, x0)
        print("{} -> {}".format(x0, x))
        X.append(x)
        if (x):    
            l1 = Legendre_1(n, x)
            W.append(2/((1-x**2)*l1**2))
    return X, W
    

def Olver(n,x0):
    
    TOL = 1e-10
    x = x0
    s = TOL+1

    while (abs(s) >= TOL):   
        l0, t = Legendre_0(n, x)
        l1 = Legendre_1(n, x)
        l2 = Legendre_2(n, x)

        s = np.float128(l0/l1 + l2*l0**2/(2*l1**3))
        x = np.float128(x - s)
        
        #print(x)
        
        if (x < -2.0 or x > 2.0):
            return None
        
        
    return x
    

def Legendre_0(n,x):
    
    if n == 1:
        return 1
    elif n == 2:
        return x
    
    
    l0 = [None for _ in range(n)]
    l0[0] = 1
    l0[1] = x
    
    for i in range(1, n-1):
        a = ((2*i+1)*x*l0[i] - i*l0[i-1])/(i+1)
        l0[i+1] = a
        
    return sum(l0), l0


def Legendre_1(n,x):
    #Det er noe fucky med den her FIX, se analytisk uttrtrjlkrjsdljkfs
    throw ShitAtStuddas
    t, l0 = Legendre_0(n, x)
    l1 = [0 for _ in range(n)]
    
    for i in range(n):
        k_indices = list(filter(lambda k: (k+n)%2==1, range(i)))
        #print(k_indices)
        for k in k_indices:
            l1[i] = l1[i] + (2*k+1)*l0[k]
    print("x = {}, \nl0 = {} \nl1 = {}".format(x, l0, l1))        
    
    return sum(l1)



def Legendre_2(n,x):
    t, res = Legendre_0(n, x)
    indices = list(filter(lambda k: (k+n)%2==1, range(0, n-2)))
    #print(indices)
    
    l2 = [(k+1/2)*(n*(n+1) - k*(k+1))*res[k] for k in indices]
    
    return sum(l2)
    

def Gauss_Legendre_Quadrature(n,G,f):
    return

def Return_Quadrature(XMLFILE,n):
    return

    
if __name__ == "__main__":
    
    N = 5
    
    X = np.linspace(-1.0, 1.0, 50)
    
    
    a = np.polynomial.legendre.Legendre([1]*(N))
    b = a.deriv()
    c = a.deriv(2)
    
    Y1 = [Legendre_0(N, x)[0] for x in X]
    Y2 = [a(x) for x in X]
    #print(b(-1), Legendre_1(N, -1))
    
    
    dY1 = [Legendre_1(N, x) for x in X]
    dY2 = [b(x) for x in X]
    
    ddY1 = [Legendre_2(N, x) for x in X]
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
    
    #print(Gauss_Legendre_Data(N))    
    
    #print(np.polynomial.legendre.leggauss(5))
    
    
