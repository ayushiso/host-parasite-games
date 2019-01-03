import matplotlib.pyplot as plt
import numpy as np
import sys

# strInt takes a list of strings and converts all elements to int
def strInt(strlist):
    intlist = []
    for item in strlist:
        intlist.append(int(item))
    return intlist

# strFlt takes a list of strings and converts all elements to float
def strFlt(strlist):
    fltlist = []
    for item in strlist:
        fltlist.append(float(item))
    return fltlist

# plotData takes a file and plots it, and saves the plot to file
def plotData(filename):

    file = filename + ".txt"

    f = open(file, 'r')

    # reading data from file
    lines = f.readlines()

    header = lines[0]

    # parsing data from file
    hostpop = strInt(lines[1].split())
    hostkfr = strFlt(lines[2].split())
    parpop = strInt(lines[3].split())
    parkfr = strFlt(lines[4].split())

    # plotting the data
    t = np.arange(len(hostpop))
    data1 = parpop
    data2 = hostpop
    data3 = parkfr
    data4 = hostkfr

    fig, ax1 = plt.subplots()

    color1 = 'r'
    color2 = 'b'
    color3 = 'k'

    plt.title(header)  # add the histogram title

    ax1.set_xlabel('Generation')
    ax1.set_ylabel('Population size', color=color3)
    ax1.plot(t, data1, color=color1, label='parasite population')
    ax1.plot(t, data2, color=color2, label='host population')
    ax1.tick_params(axis='y', labelcolor=color3)

    plt.legend()

    ax2 = ax1.twinx()  # instantiate a second axes that shares the same x-axis

    ax2.set_ylabel('Killer frequency', color=color3)  # already handled the x-label with ax1
    ax2.plot(t, data3, dashes=[6, 2], color=color1, label='parasite killer frequency')
    ax2.plot(t, data4, dashes=[6, 2], color=color2, label='host killer frequency')
    ax2.tick_params(axis='y', labelcolor=color3)

    plt.legend()

    fig.tight_layout()

    plt.savefig(filename + ".png")
    print("Finished drawing!")


if len(sys.argv) < 2:  # requires two arguments, graphing.py and mode, to begin
    print("Please give correct command line arguments!")

else:
    if sys.argv[1] == "compare":  # compare mode
        plotData("matching")
        plotData("replicator")

    elif sys.argv[1] == "sandbox":  # sandbox mode
        plotData("sandbox")

    elif sys.argv[1] == "analyze":  # analyze mode
        if len(sys.argv) != 3:
            print("Please input the number of graphs to be made!")
        else:
            n = int(sys.argv[2])
            for i in range(n):  # generate a graph for every different value of the parameter being analyzed
                filename = "parameter" + str(i+1)
                plotData(filename)

    else:
        print("Wrong input!")
