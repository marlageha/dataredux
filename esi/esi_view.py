##############################################################################
#                                                                            #
#   Easy way to view the spectrum or signal to noise of a given object       #
#   which has already been fully reduced.                                    #
#                                                                            #
#   Kareem El-Badry, 01/28/2015                                              #
#                                                                            #
##############################################################################

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
from scipy.ndimage import gaussian_filter1d


def esi_view(date, obj_id):
    
    # Load spectrum
    liness = open(str(date)+'/Final/spectra/'+str(obj_id)+'_spectrum.dat', 'r')
    data1 = liness.readlines()
    liness.close()
    
    spec_wavelengths = []
    spec_intensities = []
    
    for line in data1:
        p = line.split()
        spec_wavelengths.append(float(p[0]))
        spec_intensities.append(float(p[1]))

    # Load signal to noise:
    liness = open(str(date)+'/Final/sig_2_noise/'+str(obj_id)+'_signoise.dat', 'r')
    data1 = liness.readlines()
    liness.close()
    
    signoise_wavelengths = []
    signoise_intensities = []
    
    for line in data1:
        p = line.split()
        signoise_wavelengths.append(float(p[0]))
        signoise_intensities.append(float(p[1]))
        
    signoise_med = np.median(signoise_intensities)
    spec_med = np.median(spec_intensities)
    
    spec_intensities = gaussian_filter1d(spec_intensities, 3) 
    
    
    #Plot Spectrum 
    f = plt.figure()
    ax = f.add_subplot(111)
    ax.plot(spec_wavelengths, spec_intensities, 'k')
    ax.set_title(str(obj_id))
    ax.set_xlabel(r'$\lambda (\AA)$')
    ax.set_ylabel('normalized flux')
    ax.set_ylim([-1*spec_med,2*spec_med])
    
    #Plot Signal to Noise
    g = plt.figure()
    ax = g.add_subplot(111)
    ax.plot(signoise_wavelengths, signoise_intensities, 'k')
    ax.set_title(str(obj_id) + ' signal to noise')
    ax.set_xlabel(r'$\lambda (\AA)$')
    ax.set_ylabel('$S/N$')
    ax.set_ylim([-1*signoise_med, 3*signoise_med])