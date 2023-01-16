#! /usr/bin/python

import numpy as np
from matplotlib import pyplot as plt

if __name__ == '__main__':
    result = np.loadtxt("RESULTS.txt")
    flat = np.loadtxt("FLAT_FIELD.txt")

    corrected = result / flat

    plt.figure()
    plt.imshow(result, cmap='gray')
    plt.title('RESULT')

    plt.figure()
    plt.imshow(corrected, cmap='gray')
    plt.title('RESULT, corrected')
    plt.show()