import sys
import numpy as np
import matplotlib.pyplot as plt

def CV_plot(time, potential, current):


    plt.plot(potential, current, "k")
    plt.title("Cyclic Voltammetry")
    plt.xlabel("Potential [V]")
    plt.ylabel("Current [A]")
    plt.grid()
    #plt.legend()
    plt.tight_layout()
    plt.show()

def CD_plot(time, potential, current):

    plt.plot(time, potential, "k")
    plt.title("Charge/Discharge")
    plt.xlabel("Time [s]")
    plt.ylabel("Potential [V]")
    plt.grid()
    #plt.legend()
    plt.tight_layout()
    plt.show()

def Rate_plot(time, potential, current):

    plt.plot(time, potential, "k")
    plt.title("Rate Testing")
    plt.xlabel("Time [s]")
    plt.ylabel("Potential [V]")
    plt.grid()
    #plt.legend()
    plt.tight_layout()
    plt.show()



try:
    dir_file = sys.argv[1]
    dir_file_split = dir_file.split("/")
    dir = dir_file_split[0]
    filename = dir_file_split[1]
except IndexError:
    print("\nERROR: Please input the directory and filename when you run the code.")
    sys.exit(1)

time = []
potential = []
current = []

with open(dir_file) as infile:
    infile.readline()
    lines = infile.readlines()
    for line in lines:
        words = line.split("\t")
        time.append(words[0])
        potential.append(words[1])
        current.append(words[2])

time = np.array(time)
potential = np.array(potential)
current = np.array(current)

print("Beware! It might take time to plot the whole file.")

measurement = {"CV": CV_plot, "CD": CD_plot, "Rate": Rate_plot}[dir](time, potential, current)
