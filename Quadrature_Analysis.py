################################################################################
import sys
import os
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import Gauss_Legendre as Leg
################################################################################


def Repeated_Quadrature(pathToXml, pathToErr, n1, n2):
    relErrors = {}
    
    for n in range(n1, n2+1):
        nI, eI, relErr = Leg.Return_Quadrature(pathToXml, n)
        relErrors[n] = relErr
        
    with open(pathToErr, 'w') as file:
        for n, error in relErrors.items():
            file.write(str(n)+" "+str(error)+"\n")
    
    print("Saved to "+pathToErr)
    
    return 
    

def Convergence_Graph(pathToErr, pathToPlot):
    
    relErrors = {}
    
    f_index = os.path.split(pathToErr)[1]
    
    
    with open (pathToErr, 'r') as file:
        for line in file.readlines():
            n, err = line.split(" ")
            relErrors[n] = err
    
    X = list(relErrors.keys())
    Y = list(relErrors.values())
    
    plt.style.use("ggplot")
    plt.figure()
    plt.semilogy(X, Y, marker=".")    
    
    plt.title("$(I_{num} - I_{ex})/I_{ex}$ for "+f_index)
    plt.xlabel("$n$")
    plt.ylabel("Relative Error")
    plt.savefig(pathToPlot)
    plt.show()
    

def main():
    
    name, *args = sys.argv
    
    try:
        XMLFILE = str(args[0])
        PROGRAM = int(args[1])
        if PROGRAM == 1:
            n1 = int(args[2])
            n2 = int(args[3])
    except IndexError:
        print("Couldn't understand argument list:", args)
        print("Should be: functionIndex Program n1 n2")
        print("Example 1: Quadrature Analysis f1 1 1 10")
        print("Example 2: Quadrature Analysis f1 2")
        return
    except ValueError:
        print("One or more values are not in the correct format", args)
        print("XMLFILE=file.xml PROGRAM={1,2} N1={1, ...} N2={N1+1, ...}")
        return
    
    pathToXML = os.path.join("functions", XMLFILE)+".xml"
    pathToPlot = os.path.join("plots", XMLFILE)+".pdf"
    pathToErr = os.path.join("errors", XMLFILE)+".txt"

    if PROGRAM == 1:
        assert(n2 > n1)
        Repeated_Quadrature(pathToXML, pathToErr, n1, n2)
        
    elif PROGRAM == 2:
        Convergence_Graph(pathToErr, pathToPlot)
        
    else:
        print("Program index out of range, N_p = "+str(PROGRAM))
        
        
if __name__ == "__main__":
    main()
