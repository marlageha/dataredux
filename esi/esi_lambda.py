#################################################################
#       Finds the wavelength solution for each order.           #
#       Must first put the directory order_lists/,              #
#       which holds the emission line lists for each order,     #
#       into Calibs/. Line lists were prepared tediusly,        #
#       mostly by hand, and match what the roughpeaks.py        #
#       script finds.                                           #
#                                                               #
#       Writes polynomial wavelength solution to                #
#       Calibs/lambda_solutions.p                               #
#                                                               #
#       Kareem El-Badry, 07/15/2014                             #
#################################################################

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pyfits
import scipy.ndimage 
from scipy.optimize import curve_fit
import pdb

def esi_lambda():
    
    #Get lines - detected and expected.
    wavelengths = []
    pixels = []
    
    for order_num in range(10):
        
        #retrieve the line lists, originally from esi but sorted to be good
        lines = open('Calibs/line_lists/order_lists/wavelengths/order_' + str(order_num), 'r')
        data1 = lines.readlines()
        lines.close()
        
        lambdas = []
        usable = []
        for line in data1:
            p = line.split()
            lambdas.append(float(p[0]))
            usable.append(p[1])
        
        info = []
        for line in range(len(usable)):
            info.append([lambdas[line], usable[line]])
            
        good_lines = [info[line][0] for line in range(len(lambdas)) if 'yes' in info[line][1]]
        
        wavelengths.append(good_lines)
        
        #retrieve the pixel locations of those lines, found by esi_roughpeaks.py()
        lines = open('Calibs/line_lists/order_lists/pixels/order_' + str(order_num), 'r')
        data1 = lines.readlines()
        lines.close()
        
        pix = []
        usable = []
        for line in data1:
            p = line.split()
            pix.append(float(p[0])) #p[1] is obsolete
            usable.append(p[2])
        
        info = []
        for line in range(len(usable)):
            info.append([pix[line], usable[line]])
            
        good_pix = [info[line][0] for line in range(len(pix)) if 'yes' in info[line][1]]
        
        pixels.append(good_pix)
        
    #Now find the best wavelength solution
    solutions = []
    
    for ord_num in range(len(pixels)):
        
        pixs = pixels[ord_num]
        angstroms = wavelengths[ord_num]
        
        fit = np.polyfit(pixs, angstroms, 5)
        solution = np.poly1d(fit)
        
        solutions.append(solution)

    pickle.dump(solutions, open('Calibs/lambda_solutions.p', 'wb'))
    
    
