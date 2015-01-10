#############################################################################
#      Subtracts sky in wavelength space                                    #
#      Reads in skymasks and 2d wavelenght solutions from pickle file       #
#      Assigns wavelenght to each pixel                                     #
#      Fits a 1d sky I(lambda) to the sky pixels as a spline                #
#      Applies spline to all pixels and subtracts sky                       #
#      Writes to Calibs/sky_sub/objid_skysub.fits                           #
#                                                                           #
#      Kareem El-Badry, 10/07/2014                                          #
#############################################################################

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pyfits
from scipy.signal import argrelextrema
from scipy.interpolate import LSQUnivariateSpline

def esi_lambdaskysub():
    
    print "reading masks, wavelengths solutions..."
    solutions2d = pickle.load(open('Calibs/solution2d.p', 'rb'))
    all_order_mask = pickle.load(open('Calibs/all_order_masks.p', 'rb'))
    sky_mask = pickle.load(open('Calibs/sky_mask.p', 'rb'))
    background = pickle.load(open('Calibs/background_mask.p', 'rb'))
    
    im1 = open('Logs/esi_info.dat','r')
    data1 = im1.readlines()
    im1.close()

    filename = []
    dateobs = []
    objname = []
    imgtype = []
    ra = []
    dec = []
    exptime = []
    usable = []

    for line in data1:
        p = line.split()
        filename.append(p[0])
        dateobs.append(p[1])
        objname.append(p[2])
        imgtype.append(p[3])
        ra.append(p[4])
        dec.append(p[5])
        exptime.append(p[6])
        usable.append(p[7])

    #Rewrite in a more convenient array with format array[line][element]    
    alldata = []

    for line in range(len(usable)):
        alldata.append([filename[line],dateobs[line],objname[line],imgtype[line],
                        ra[line],dec[line],exptime[line],usable[line]])

    #Find good files:
    good = []
    for line in range(len(alldata)):
       if "yes" in alldata[line][7]:
           good.append(alldata[line])

    #Find list of lamps 
    names = []
    for line in range(len(alldata)):
        if "Object" in alldata[line][3] and float(alldata[line][6]) > 600:
            names.append(alldata[line][2])

    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place

    for obj_id in objects:
        
        print "sky subtracting " +str(obj_id)
        
        lamp = pyfits.getdata('Calibs/reduced/' + str(obj_id) + '_mean.fits')
        noise = pyfits.getdata('Calibs/variance/'+str(obj_id)+'_noise.fits')
    
        #to hold skysubtracked image
        blanc = np.zeros((4096, 2045))
        blanc_error = np.zeros((4096, 2045))
    
        #loop through orders
        for num in range(10):
        
            print "working on order " + str(num)
            sol = solutions2d[num]
        
            heights = lamp[all_order_mask[num]]
            sky_heights = lamp[sky_mask[num]]
            
            error_heights = noise[all_order_mask[num]]
            skyerror_heights = noise[sky_mask[num]]
            
            if 0 < num < 3:
    
                med = np.median(sky_heights)
                bad_mask = np.abs(sky_heights) > 15
    
                sky_heights[bad_mask] = med
            
            y, x = np.mgrid[:4096, :2045]
            
            #gives the wavelengths of all pixels in order
            lambdas = sol(x[all_order_mask[num]], y[all_order_mask[num]]) 
            sky_lambdas = sol(x[sky_mask[num]], y[sky_mask[num]])
            
            #Make knot sequence 
            blank = np.zeros((4096, 2045))

            blank[all_order_mask[num]] = lamp[all_order_mask[num]]

            ord_tot = []
            for line in range(4096):
                row = np.array(blank[line])
                nonzero = row[np.where(row != 0)[0]]
                row_avg = np.median(nonzero)
                ord_tot.append(row_avg)
            ord_tot = np.array(ord_tot)

            peaks = argrelextrema(ord_tot, np.greater)[0]
            peaks = [peaks[line] for line in range(len(peaks)) if ord_tot[peaks[line]]>30]
            peaks = np.array(peaks)

            #find x coords of order center
            centers = []
            for line in range(len(peaks)):
                row = np.array(blank[peaks[line]])
                centers.append(np.where(row != 0)[0][50])
            centers = np.array(centers)

            peak_lambdas = sol(centers, peaks)

            #make line mask
            masks = []
            for line in range(len(peak_lambdas)):
                mask_below = lambdas > peak_lambdas[line]-4
                mask_above = lambdas < peak_lambdas[line]+4
                mask = mask_below*mask_above
                masks.append(mask)
            tot_mask = np.sum(masks, axis = 0, dtype = bool)
    
            if 9 > num > 6:
                spaceing = 0.18
                deg = 1
            elif num < 7:
                spaceing = 0.13
                deg = 1
            elif num > 8:
                spaceing = 0.23
            
            cont_space = (sky_lambdas[-1] - sky_lambdas[0])/4096
            
            
            #continuum knots
            cont_knots = np.arange(sky_lambdas[0]+.1, sky_lambdas[-1]-.1, cont_space)
            line_knots = np.arange(sky_lambdas[0]+.1, sky_lambdas[-1]-.1, spaceing)
            
            #fit splines
            sp = LSQUnivariateSpline(sky_lambdas, sky_heights, k=3, t = cont_knots)
            bsp = LSQUnivariateSpline(sky_lambdas, sky_heights, k=deg, t = line_knots)
            
            sp_error = LSQUnivariateSpline(sky_lambdas, skyerror_heights, k=3, t = cont_knots)
            bsp_error = LSQUnivariateSpline(sky_lambdas, skyerror_heights, k=deg, t = line_knots)
            
            skei = sp(lambdas)
            skei_error = sp_error(lambdas)
            
            skei[tot_mask] = bsp(lambdas)[tot_mask]
            skei_error[tot_mask] = bsp_error(lambdas)[tot_mask]
            
            corrected = heights - skei
            corrected_error = np.sqrt(error_heights**2 + skei_error**2)

            blanc[all_order_mask[num]] = corrected
            blanc_error[all_order_mask[num]] = corrected_error
            #End order loop
            
        blanc[background] = lamp[background]
        blanc_error[background] = noise[background]
        
        fits = pyfits.PrimaryHDU(blanc)
        fits.writeto('Calibs/sky_sub/'+str(obj_id)+'_skysub.fits', clobber = "True")
        
        fits = pyfits.PrimaryHDU(blanc_error)
        fits.writeto('Calibs/variance/'+str(obj_id)+'_noise.fits', clobber = "True")
        