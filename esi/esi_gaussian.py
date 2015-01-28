##############################################################################
#                                                                            #
#    Gives the center of a line more accurately than roughpeaks(). Fits a    #
#    single gaussian to all peaks detected by the roughpeaks() algorithm     #
#    and returns the center of the gaussian (corresponding to the 'mean      #
#    variable.) Writes these new centers to better_order[obj_id].dat in the  #
#    same directory as the other line lists.                                 #
#                                                                            #
#    One thing worth playing with is fitting a sum of gaussians rather than  #
#    a single gaussian to each line. See fitting.custom_model1d.             #
#                                                                            #
#    Kareem El-Badry, 07/25/2014                                             #
#                                                                            #
##############################################################################

from __future__ import division
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pyfits
import scipy.ndimage 

def esi_gaussian(date):
    
    #Get edge masks from file
    print 'loading masks...'
    orders_mask = pickle.load(open(str(date)+'/Calibs/orders_mask_'+str(date)+'.p', 'rb'))
    background_mask = pickle.load(open(str(date)+'/Calibs/background_mask_'+str(date)+'.p', 'rb')) 
    all_order_mask = pickle.load(open(str(date)+'/Calibs/all_order_masks_'+str(date)+'.p', 'rb'))

    #get names of lamps, first read in log:
    im1 = open(str(date)+'/Logs/esi_info_'+str(date)+'.dat','r')
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
        if "Line" in alldata[line][3]:
            names.append(alldata[line][2])
        
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place

    obj_id = objects[0] #HgNeXeCuAr.  

    lamp = pyfits.getdata(str(date)+'/Calibs/reduced/' + str(obj_id) + '_mean.fits')

    print "finding peaks..."
    gaus_orders = []
    for num in range(10): 

        #load order to otherwise blank array
        blank = np.zeros((4096, 2045))
        blank[all_order_mask[num]] = lamp[all_order_mask[num]]

        #average over x pixels
        ord_tot = []
        for line in range(4096):
            row_avg = np.mean(blank[line])
            ord_tot.append(row_avg)
        ord_tot = np.array(ord_tot)

        #fit a polynomial to baseline with iterative fit
        x = range(4096)
        y = ord_tot

        done = False
        while not done:
            done = True
            fit = np.poly1d(np.polyfit(x, y, 10))
            resid = y - fit(x)
            std = np.std(resid)
            badindices = np.where(resid > 1.5*std)[0]
            if badindices.size > 0:
                done = False
                x = np.delete(x, badindices)
                y = np.delete(y, badindices)
        #subtract baseline
        ord_tot = ord_tot - fit(range(4096))
    

        #smooth slightly, find slope
        smoothed = scipy.ndimage.gaussian_filter1d(ord_tot, 5)

        smooted = scipy.ndimage.gaussian_filter1d(ord_tot, 1.5)
        deriv = scipy.ndimage.sobel(smooted)#mayber deriv of smoothed?

        #find local maxes
        mask = np.r_[True, smoothed[1:] > smoothed[:-1]] & np.r_[smoothed[:-1] > smoothed[1:], True]
        xcenr = np.arange(4096)[mask]
        mcenr = smoothed[mask]

        #makes easier to work with for present
        xcenr = list(xcenr)
        mcenr = list(mcenr)

        #Make list of list of slopes before and after each peak:
        slopes_after = []
        slopes_before = []
        mins = []
        maxes = []

        for line in range(len(xcenr)):
            try: 
                slopes_before.append([deriv[xcord] for xcord in range(xcenr[line]-5, xcenr[line])])
                slopes_after.append([deriv[xcord] for xcord in range(xcenr[line]+1, xcenr[line]+6)])
            except IndexError:
                slopes_after.append(deriv[xcenr]) #if its at the edge; no good anyway
                slopes_before.append(deriv[xcenr])

            mins.append(np.min(slopes_after[line]))
            maxes.append(np.max(slopes_before[line]))

        mins = np.array(mins)
        maxes = np.array(maxes)

        #slopes at least steep enough on either side
        good_max = [line for line in range(len(xcenr)) if mins[line] < -0.001 and maxes[line] > 0.001]
        xcenr = np.array(xcenr)
        mcenr = np.array(mcenr)

        xcenr = xcenr[good_max]
        mcenr = mcenr[good_max]

        #remove peaks at ends of orders:
        mcenr = mcenr[(xcenr > 5) & (xcenr < 4090)] #important that mcenr is first!
        xcenr = xcenr[(xcenr > 5) & (xcenr < 4090)] #otherwise the first line changes the second.
 
        #Fit gaussians to each peak
        gaus_means = []
        fit_g = fitting.LevMarLSQFitter() #nonlinear least-square minimization

        for line in range(len(xcenr)):
    
            amp = mcenr[line]
            mean = xcenr[line]
            std = 5
            g_init = models.Gaussian1D(amplitude=amp, mean = mean, stddev = std)
    
            #Set constraints
            g_init.stddev.bounds = (0,12)
            g_init.mean.bounds = (xcenr[line+0]-0.6, xcenr[line+0]+0.6)
            g_init.amplitude.bounds = (0.99*mcenr[line+0], 1.01*mcenr[line+0])
    
            #We only fit the 7 pixels centered on the center. This gives best fit. 
            g = fit_g(g_init, range(xcenr[line]-3, xcenr[line]+4), smoothed[xcenr[line]-3: xcenr[line]+4])
            #plt.plot(range(4096), g(range(4096)))
    
            gaus_means.append(g.mean.value)
            #now gaus_means is a better version of xcenr. We want to feed it 
            #into order line lists
        gaus_orders.append(gaus_means)
    
    #now gaus_orders has the important info. Next load order_nums and match.
    #Get lines - detected and expected.
    wavelengths = []
    pixels = []

    for order_num in range(10):
        #retrieve the line lists, originally from esi but sorted to be good
        
        lines = open(str(date)+'/Calibs/line_lists/order_lists/wavelengths/order_' 
                    + str(order_num) +'_'+str(date)+'.dat', 'r')
        data1 = lines.readlines()
        lines.close()
    
        
        #FIX FOR OLD FORMAT;
        try:
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
            
        #NEW FORMAT
        except (ValueError, IndexError):
            good_lines = []
            
            for line in data1:
                p = line.split()
            
                if len(p) == 1: # if there's no # sign
                    good_lines.append(float(p[0]))

    
        wavelengths.append(good_lines)
    
        #retrieve the pixel locations of those lines, found by esi_roughpeaks.py()
        lines = open(str(date)+'/Calibs/line_lists/order_lists/pixels/order_' 
                    + str(order_num) +'_'+str(date)+'.dat', 'r')
        data1 = lines.readlines()
        lines.close()
        
        #FIX FOR OLD FORMAT;
        try:
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
            
        except (ValueError, IndexError):
            good_pix = []
            
            for line in data1:
                p = line.split()
            
                if len(p) == 2: # if there's no # sign
                    good_pix.append(float(p[0]))
    
        pixels.append(good_pix)
        
        if len(pixels[order_num]) != len(wavelengths[order_num]):
            print "Oh no. You have a different number of lines pixels than line is order ", str(order_num)

    #Now match gaus_orders to pixels
    all_good_peaks = []
    for order_num in range(10):
    
        gaus_order = gaus_orders[order_num]
        ord_pixels = pixels[order_num]
        good_peaks = []
    
        #loop through detected gaussian peaks
        #will work as long as no lines are within a pixel of eachother. 
        for peak in gaus_order:
            for line in ord_pixels:
                dist = peak - line
                if np.abs(dist) < 1:
                
                    #hold good, revevant peaks
                    good_peaks.append(peak)
                
        #hold all good peaks, all orders            
        all_good_peaks.append(good_peaks)
    
    #write thes to file:
    for num in range(len(all_good_peaks)):
    
        file = open(str(date)+'/Calibs/line_lists/order_lists/pixels/better_order_'
                    +str(num)+'_'+str(date)+'.dat', 'w')

        for line in range(len(all_good_peaks[num])):
            file.write(' '.join([str(all_good_peaks[num][line]), 'yes', '\n']))
        file.close()