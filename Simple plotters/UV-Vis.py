import numpy as np # For arrays and math
import matplotlib.pyplot as plt # Plotting module
import matplotlib.axes as ax # Customizing tick marks on plot

def read(filename):
    """
    Reads an inputfile and returns its wavelength and transmission as a nested list
    File as followed:
    1 info
    2 units
    3 data1 data2
    3 ...
    """
    wavelength = []
    transmission = []
    with open(filename) as infile: # Opens the file and reads it
        infile.readline() # Reads and discards the 1st line (info)
        infile.readline() # Reads and discards the 2nd line (units)
        lines = infile.readlines() # Reads and stores all data

        if "," in lines[0]: # If the separator is a comma
            for line in lines: # Iterates over each line
                line = line.replace(",", ".") # Replaces comma with period
                words = line.split() # Splits the line into words for each whitespace it finds
                wavelength.append(float(words[0])) # Stores the wavelength as a number to the list "wavelength"
                transmission.append(float(words[1])) # Stores the transmission as a number to the list "transmission"
        elif "." in lines[0]: # If the separator is a period
            for line in lines: # Iterates over each line
                words = line.split() # Splits the line into words for each whitespace it finds
                wavelength.append(float(words[0])) # Stores the wavelength as a number to the list "wavelength"
                transmission.append(float(words[1])) # Stores the transmission as a number to the list "transmission"
        else:
            print("Could not determine if separator was \",\" or \".\". Defaulting to \".\".")
            for line in lines: # Iterates over each line
                words = line.split() # Splits the line into words for each whitespace it finds
                wavelength.append(float(words[0])) # Stores the wavelength as a number to the list "wavelength"
                transmission.append(float(words[1])) # Stores the transmission as a number to the list "transmission"
    return wavelength, transmission

def plot(x_and_y, title="../UV-Vis/NO_TITLE.txt"):
    """
    Plots a graph based on a list/array of x- and y-values.
    """

    title = title[10:-4] # Removes path to file and filetype from title

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
    plt.xlabel("Wavelength [nm]") # x-axis name

    # Sets y-axis name based on type of measurement
    if "trans" in title:
        plt.ylabel("Transmission [T\%]")
    elif "reflect" in title:
        plt.ylabel("Reflectance [R\%]")
    else:
        print("Could not determine if measurement was trans or reflect, defaulting axix name to transmission")
        plt.ylabel("Transmission")

    plt.tight_layout() # Creates a tight layout
    plt.savefig("UV_Vis_images/" + title + ".png") # Saves the figure to the specified location and name
    plt.show() # Shows the figure


# Reads UV-Vis files
ETN4013_Yb2Qz3 = read("../UV-Vis/ETN4013_Yb2Qz3_210226_250-850_M_trans_glass.txt")
ETN4029_Y2Qz3 = read("../UV-Vis/ETN4029_Y2Qz3_210428_250-850_M_trans_glass.txt")
ETN4029_Y2Qz3_trans_sphere = read("../UV-Vis/ETN4029_Y2Qz3_210428_250-850_M_trans_glass_sphere.txt")
ETN3029_Y2Qz3_reflect_sphere = read("../UV-Vis/ETN4029_Y2Qz3_210428_250-850_M_reflect_glass_sphere.txt")

# Plots the UV-Vis files
plot(ETN4013_Yb2Qz3, "../UV-Vis/ETN4013_Yb2Qz3_210226_250-850_M_trans_glass.txt")
plot(ETN4029_Y2Qz3, "../UV-Vis/ETN4029_Y2Qz3_210428_250-850_M_trans_glass.txt")
plot(ETN4029_Y2Qz3_trans_sphere, "../UV-Vis/ETN4029_Y2Qz3_210428_250-850_M_trans_glass_sphere.txt")
plot(ETN3029_Y2Qz3_reflect_sphere, "../UV-Vis/ETN4029_Y2Qz3_210428_250-850_M_reflect_glass_sphere.txt")
