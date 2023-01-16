#! /usr/bin/python

from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import string
import numpy as np
from spectra_converter import read_spectrum

def read_3d_vectors(filename: string) -> tuple:
    energies = []
    origins = []
    directions = []
    with open(filename) as px_f:
        lines = px_f.readlines();
        for line in lines:
            line = line.replace('(', ' ').replace(')', ' ').replace(',', ' ')
            values = line.split()
            values = np.asfarray(values)
            energies.append(values[0])
            origins.append(values[1:4])
            directions.append(values[4:7])
    return energies, origins, directions

def plot_space(origins: list, directions: list, threshold:int = 200):
    plt.figure()
    ax = plt.axes(projection ='3d')

    each = 1
    if len(origins) > threshold: each = len(origins) // threshold
    
    t = np.linspace(0, 2, 3)
    for i in range(len(directions)):
        if not i % each == 0: continue
        ax.plot3D(directions[i][0]*t, directions[i][1]*t, directions[i][2]*t, color='red')

    plt.figure()
    ax = plt.axes(projection ='3d')
    for i in range(len(directions)):
        if not i % each == 0: continue
        ax.scatter(origins[i][0], origins[i][1], origins[i][2])

if __name__ == '__main__':

    energies, origins, directions = read_3d_vectors("ph_vis_data/photons_finite.txt")
    
    plot_space(origins, directions)

    plt.figure()
    plt.hist(energies, bins=150)

    energies_should_be, photons_should_be = read_spectrum("spectra/SPECTRA_120kVp_17deg_1Al.txt")
    plt.figure()
    plt.plot(energies_should_be, photons_should_be)

    plt.show()
    