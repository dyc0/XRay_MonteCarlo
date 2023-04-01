#! /usr/bin/python

from matplotlib import pyplot as plt
import sys


def parse_spectrum(path: str) -> list:
    with open(path) as in_file:
        str_energies = in_file.read().split(" ")[:-2]
        energies = [float(i) for i in str_energies]
        return energies


def read_spectrum(in_filename: str) -> tuple:
    energies = []
    photons = []

    with open(in_filename) as f:
        lines = f.readlines()

        for line in lines[18:]:
            substrs = line.split("  ")
            energies.append(float(substrs[0]))
            photons.append(float(substrs[1][:-1]))
        
    return energies, photons


if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Zadati ime ulaznog fajla')
        quit(-1)
    
    energies = parse_spectrum("spectrum_" + sys.argv[1] + ".txt")

    plt.figure()
    plt.subplot(1,2,1)
    plt.hist(energies, bins=150)
    plt.title('Generated photons')

    energies_should_be, photons_should_be = read_spectrum("theoretical_" + sys.argv[1] + ".txt")
    plt.subplot(1,2,2)
    plt.plot(energies_should_be, photons_should_be)
    plt.title('Expected ideal')

    plt.show()
