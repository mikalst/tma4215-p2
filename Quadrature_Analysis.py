import sys
import os

import matplotlib.pyplot as plt
import Gauss_Legendre as Leg


def Repeated_Quadrature(pathToXml, pathToErr, n1, n2):
    relErrors = {}

    for n in range(n1, n2+1):
        nI, eI, relErr = Leg.Return_Quadrature(pathToXml, n)
        relErrors[n] = relErr

    with open(pathToErr, 'w') as file:
        for n, error in relErrors.items():
            file.write(str(n)+" "+str(error)+"\n")

    print("Saved to", pathToErr)

    return


def Convergence_Graph(pathToErr, pathToPlot):

    relX = []
    relErrors = []

    f_index = os.path.split(pathToErr)[1]


    with open (pathToErr, 'r') as file: 
        mylist = file.read().splitlines()
        for el in mylist:
            n, err = el.split(" ")
            relX.append(int(n))
            relErrors.append(float(err))

    plt.style.use("ggplot")
    plt.figure()
    plt.semilogy(relX, relErrors, marker=".")

    plt.title("$(I_{num} - I_{ex})/I_{ex}$ for "+f_index)
    plt.xlabel("$n$")
    plt.ylabel("Relative Error")
    plt.savefig(pathToPlot)
    plt.show()


def main(*args, **kwargs):
    # Program can accept STDIN arguments or parameters
    args = args if args else sys.argv[1:]

    if not args:
        XMLFILE = "f1"
        PROGRAM = 1
        n1 = 1
        n2 = 10

    else:
        try:
            XMLFILE = str(args[0])
            PROGRAM = int(args[1])
            if PROGRAM == 1 or PROGRAM == 12:
                n1 = int(args[2])
                n2 = int(args[3])
        except IndexError or ValueError:
            print("Argument list is incomplete or assigned incorrect values:", args)
            print("Should be: functionIndex Program n1 n2")
            print("Example 1: Quadrature Analysis f1 1 1 10")
            print("Example 2: Quadrature Analysis f1 2")
            return

    if not(os.path.isfile(XMLFILE)):
        f_index, extension = os.path.splitext(XMLFILE)

        pathToXML = os.path.join("functions", f_index)+".xml"
        pathToPlot = os.path.join("plots", f_index)+".pdf"
        pathToErr = os.path.join("errors", f_index)+".txt"

    else:
        path, fileName = os.path.split(XMLFILE)
        f_index, extension = os.path.splitext(fileName)

        pathToXML = XMLFILE
        pathToPlot = os.path.join("plots", f_index)+".pdf"
        pathToErr = os.path.join("errors", f_index)+".txt"


    if PROGRAM == 1:
        assert(n2 > n1)
        print("Running Leg-Gauss on", pathToXML,
              "with n = {}, {}, ... {}".format(n1, (n1+1), n2))
        Repeated_Quadrature(pathToXML, pathToErr, n1, n2)

    elif PROGRAM == 2:
        print("Plotting from", pathToPlot)
        Convergence_Graph(pathToErr, pathToPlot)

    elif PROGRAM  == 12:
        assert(n2 > n1)
        print("Running Leg-Gauss on", pathToXML,
              "with n = {}, {}, ... {}".format(n1, (n1+1), n2))
        Repeated_Quadrature(pathToXML, pathToErr, n1, n2)
        print("Plotting from", pathToPlot)
        Convergence_Graph(pathToErr, pathToPlot)



if __name__ == "__main__":
    main()
