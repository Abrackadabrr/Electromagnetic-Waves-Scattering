import matplotlib.pyplot as plt
import numpy as np
import csv

filename = ("/media/evgen/SecondLinuxDisk/4_level/"
            "Electromagnetic-Waves-Scattering/vtk_files/examples/cylinder/sigma_back_2002.csv")

sigma = []
angle = []

with open(filename) as csvfile:
    r = csv.reader(csvfile, delimiter=',')
    next(r)
    for line in r:
        sigma.append(float(line[0]))
        angle.append(float(line[1]))

polar = False

if polar:
    fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
    ax.plot(angle, sigma, color='red')
    ax.set_rmax(max(sigma) * 1.05)
    ax.set_rlabel_position(-90)  # Move radial labels away from plotted line
    ax.grid(True)

    ax.set_title("Effective dispersion area, plane (1m x 1m)", va='bottom')
else:
    fig, ax = plt.subplots()
    ax.plot(angle[::1], 10 * np.log10(sigma[::1]), label=f"sigma")
    ax.grid(True)
    ax.set_xlabel("Angle")
    ax.set_ylabel("Effective dispersion area")
    ax.minorticks_on()
    ax.set_title("Plane (1m x 1m), lambda = 25 sm", va='bottom')
plt.show()
