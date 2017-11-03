################################################################################
import re
import sys
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import Gauss_Legendre as Leg
################################################################################


def Repeated_Quadrature(METHOD,XMLFILE,n1,n2):
    relErrors = {}
    
    ascName = "errors/" + re.search('[a-z]+[0-9]', XMLFILE).group(0) + ".txt"
    
    for n in range(n1, n2+1):
        nI, eI, relErr = Leg.Return_Quadrature(XMLFILE, n)
        relErrors[n] = relErr
        
    with open(ascName, 'w') as file:
        for n, error in relErrors.items():
            file.write(str(n)+" "+str(error)+"\n")
    
    print("Saved to "+ascName)
    
    return 
    

def Convergence_Graph(ASCFILE):
    
    relErrors = {}
    
    # Determine function number for easy plotting
    index = re.search('[0-9]', ASCFILE).group(0)
    
    with open (ASCFILE, 'r') as file:
        for line in file.readlines():
            n, err = line.split(" ")
            relErrors[n] = err
    
    X = list(relErrors.keys())
    Y = list(relErrors.values())
    
    plt.style.use("ggplot")
    plt.figure()
    plt.semilogy(X, Y, marker=".")    
    
    plt.title("$(I_{num} - I_{ex})/I_{ex}$ for $f_{"+index+"}$")
    plt.xlabel("$n$")
    plt.ylabel("Relative Error")
    plt.savefig("plots/f"+str(index)+".pdf")
    plt.show()
    

def main():
    
    name, *args = sys.argv
    
    try:
        XMLFILE = str(args[0])
        PROGRAM = int(args[1])
        n1 = int(args[2])
        n2 = int(args[3])
    except IndexError:
        print("Couldn't undertand argument list", args)
        print("Should be: XMLFILE PROGRAM N1 N2")
    except TypeError:
        print("One or more values are not in the correct format", args)
        print("XMLFILE=file.xml PROGRAM={1,2} N1={1, ...} N2={N1+1, ...}")
    
    filename = re.search("(?<=functions/)f+[0-9]", XMLFILE).group(0)
    ascfile = "errors/"+filename+".txt"

    if PROGRAM == 1:
        assert(n2 > n1)
        Repeated_Quadrature("method?", XMLFILE, n1, n2)
        
    elif PROGRAM == 2:
        Convergence_Graph(ascfile)
        
    else:
        print("Program index out of range, N_p = "+str(PROGRAM))
        
        
if __name__ == "__main__":
    main()
