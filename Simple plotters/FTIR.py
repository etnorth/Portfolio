import numpy as np # For arrays and math
import matplotlib.pyplot as plt # Plotting module
import matplotlib.axes as ax # Customizing tick marks on plot

def read(filename):
    """
    Reads an inputfile and returns its wavenumber and transmittance as a nested list
    File as followed:
    1 info
    2 units
    3 data1 data2
    3 ...
    """
    wavenumber = []
    transmittance = []
    with open(filename) as infile: # Opens the file and reads it
        lines = infile.readlines() # Reads and stores all data

        if "," in lines[0]: # If the separator is a comma
            for line in lines: # Iterates over each line
                line = line.replace(",", ".") # Replaces comma with period
                words = line.split() # Splits the line into words for each whitespace it finds
                wavenumber.append(float(words[0])) # Stores the wavenumber as a number to the list "wavenumber"
                transmittance.append(float(words[1])) # Stores the transmittance as a number to the list "transmittance"
        elif "." in lines[0]: # If the separator is a period
            for line in lines: # Iterates over each line
                words = line.split() # Splits the line into words for each whitespace it finds
                wavenumber.append(float(words[0])) # Stores the wavenumber as a number to the list "wavenumber"
                transmittance.append(float(words[1])) # Stores the transmittance as a number to the list "transmittance"
        else:
            print("Could not determine if separator was \",\" or \".\". Defaulting to \".\".")
            for line in lines: # Iterates over each line
                words = line.split() # Splits the line into words for each whitespace it finds
                wavenumber.append(float(words[0])) # Stores the wavenumber as a number to the list "wavenumber"
                transmittance.append(float(words[1])) # Stores the transmittance as a number to the list "transmittance"
    return wavenumber, transmittance

def plot(x_and_y, title="../FTIR/NO_TITLE.txt"):
    """
    Plots a graph based on a list/array of x- and y-values.
    """

    title = title[8:-4] # Removes path to file and filetype from title

    fig = plt.figure() # Creates a figure to plot in
    ax = fig.gca() # Used for defining tick-marks and direction
    ax.tick_params(bottom=True, top=True, left=True, right=True) # Tick-marks on all sides
    ax.tick_params(axis="x", direction="in") # Point tick-marks inwards (x-axes)
    ax.tick_params(axis="y", direction="in") # Point tick-marks inwards (y-axes)

    # x_and_y = (x,y) where x and y are lists
    x = x_and_y[0] # Splits x_and_y into x
    y = x_and_y[1] # Splits x_and_y into y

    ax.plot(x, y, "k") # Plots a black line
    plt.title(title) # Graph title
    plt.xlabel("Wavenumber [nm$^{-1}$]") # x-axis name
    plt.ylabel("Transmittance [\%]")
    plt.tight_layout() # Creates a tight layout
    plt.savefig("FTIR_images/" + title + ".png") # Saves the figure to the specified location and name
    plt.show() # Shows the figure


# Reads UV-Vis files
ETN4015_Y2Qz3_normal = np.array(read("../FTIR/ETN4015_Y2Qz3_normaldetector_res4_scantime64_reflec.dpt"))
ElectroSteelBackground_normal = np.array(read("../FTIR/ElectroSteelBackground_normaldetector_res4_scantime64_reflec.dpt"))
ETN4015_Y2Qz3_N2 = np.array(read("../FTIR/ETN4015_Y2Qz3_N2detector_res4_scantime64_reflec.dpt"))
ElectroSteelBackground_N2 = np.array(read("../FTIR/ElectroSteelBackground_N2detector_res4_scantime64_reflec.dpt"))

#Subtracts the background from the measurement
ETN4015_Y2Qz3_normal[1] = ETN4015_Y2Qz3_normal[1]-ElectroSteelBackground_normal[1][6119:8001] # Checked indexes manually
ETN4015_Y2Qz3_N2[1] = ETN4015_Y2Qz3_N2[1]-ElectroSteelBackground_N2[1][6119:8001] # Checked indexes manually


# Plots the UV-Vis files
plot(ETN4015_Y2Qz3_normal, "../FTIR/ETN4015_Y2Qz3_normaldetector_res4_scantime64_reflec.dpt")
plot(ETN4015_Y2Qz3_N2, "../FTIR/ETN4015_Y2Qz3_N2detector_res4_scantime64_reflec.dpt")
plot((ElectroSteelBackground_normal[0][6119:8001],ElectroSteelBackground_normal[1][6119:8001]), "../FTIR/ElectroSteelBackground_normaldetector_res4_scantime64_reflec.dpt")
plot((ElectroSteelBackground_N2[0][6119:8001],ElectroSteelBackground_N2[1][6119:8001]), "../FTIR/ElectroSteelBackground_N2detector_res4_scantime64_reflec.dpt")
