#! /usr/bin/python

import numpy as np
from matplotlib import pyplot as plt
import sys

if __name__ == '__main__':
    if len(sys.argv) < 2:
        print('Potrebno je uneti ulaznu sliku.')
        exit(-1)
    
    
    result = np.loadtxt(sys.argv[1])
    flat = np.loadtxt("FLAT_FIELD.txt")
    ellipsoid = np.loadtxt("ELLIPSOID.txt")

    corrected = result / flat
    corrected2 = result/flat - ellipsoid/flat

    plt.figure()

    plt.subplot(2,2,1)
    plt.imshow(result, cmap='gray')
    plt.title('RESULT')

    plt.subplot(2,2,2)
    plt.imshow(corrected, cmap='gray')
    plt.title('RESULT, corrected')
    
    plt.subplot(2,2,3)
    plt.imshow(corrected2, cmap='gray')
    plt.title("Corrected for ellipsoid shadow")

    plt.show()