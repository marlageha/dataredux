##############################################################################
#                                                                            #
#    Fits a 2-dimensional wavelength solution to each order, which is        #
#    necessary because the emission lines are not quite horizontal -- that   #
#    is, on a given order, wavelength is not constant for a given y coord.   #
#                                                                            #
#    Breaks each order into about 15 colums that are 10 pixels wide and run  #
#    parallel to the edge of the order. Finds the locations of the peaks in  #
#    in each column, first with the algorithm from roughpeaks() and then     #
#    more accurately using the gaussian-fitting technique developed in       #
#    gauss().                                                                #
#                                                                            #
#    For each order, creates x_coords and y_coords, which are arrays of all  #
#    the locations of all the measured peaks. Also creates wavelengths2d,    #
#    an array of all the wavelengths corresponding to the locations of the   #
#    peaks. Then fits a 2-d polynomial to these points, which can tell us    #
#    wavelength of any pixel in the order. Fit can be improved if we find    #
#    more lines.                                                             #
#                                                                            #
#    Writes array of 10 2d wavelength solutions, one for each order, to      #
#    Calibs/solution2d.p                                                     #
#                                                                            #
#    Kareem El-Badry, 07/25/2014                                             #
#                                                                            #
##############################################################################
 
 
from __future__ import division
from astropy.modeling import models, fitting
import matplotlib.pyplot as plt
import numpy as np
import pdb
import pickle
import pyfits
import scipy.ndimage 

def esi_2dsolution(date):
    
    print "loading masks..."
    all_order_mask = pickle.load(open(str(date)+'/Calibs/all_order_masks_'+str(date)+'.p', 'rb'))
    orders_mask = pickle.load(open(str(date)+'/Calibs/orders_mask_'+str(date)+'.p', 'rb'))

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

    #Find list of lines 
    names = []
    for line in range(len(alldata)):
        if "Line" in alldata[line][3]:
            names.append(alldata[line][2])
    lines = np.array(list(set(names)))
    lines.sort() #ascending order, modify in place

    #Find list of objects 
    names = []
    for line in range(len(alldata)):
        if ("Object" in alldata[line][3] and float(alldata[line][6]) > 600):
            names.append(alldata[line][2])
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place
    
    esiorders = pickle.load(open(str(date)+'/Calibs/order_edges_'+str(date)+'.p', 'rb'))
    
    #Load saved order linelists
    #now gaus_orders has the important info. Next load order_nums and match.
    #Get lines - detected and expected.
    wavelengths = []
    pixels = []

    for order_num in range(10):
        #retrieve the line lists, originally from esi but sorted to be good
        
        liness = open(str(date)+'/Calibs/line_lists/order_lists/wavelengths/order_' 
                    + str(order_num) +'_'+str(date)+'.dat', 'r')
        data1 = liness.readlines()
        liness.close()
    
        
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
        liness = open(str(date)+'/Calibs/line_lists/order_lists/pixels/order_' 
                    + str(order_num) +'_'+str(date)+'.dat', 'r')
        data1 = liness.readlines()
        liness.close()
        
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
    
    #find peaks along different columns in each order 
    pickle_obj = []
    
    for num in range(10): 
        
        if num == 9:
            poly_fit_ord = 5
        else:
            poly_fit_ord = 5
        
        print "reducing order " +str(num) +'...'
        
        obj_id = lines[0]
        sky_id = objects[0]
            
        lamp = pyfits.getdata(str(date)+'/Calibs/reduced/' + str(obj_id) + '_mean.fits')
        
        column_width = 20

        order_left = esiorders[num].xl
        order_right = esiorders[num].xr

        width = order_right(range(4096))-order_left(range(4096))
        width = np.abs(width)
        max_wid = np.max(width)

        #make masks at 5 pixel intervals
        lx = 2045
        ly = 4096
        X, Y = np.ogrid[0:ly, 0:lx]

        gaus_orders = []
        for line in range(int(max_wid/column_width)-1): #leave off the bad part
            print "column " + str(line) + '...'
            mask_left = order_left(X) + column_width*line < Y
            mask_left5 = order_left(X) + column_width*(line+1)  > Y
            mask_right = order_right(X) > Y

            mask = mask_left*mask_left5*mask_right

            #Now work down mask (5 pixels wide) and add into a 1d spectrum
            #load order to otherwise blank array
            blank = np.zeros((4096, 2045))
            blank[mask] = lamp[mask]

            #average over x pixels
            ord_tot = []
            for line in range(4096):
                row_avg = np.mean(blank[line])*max_wid #normalized to before
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
                g_init.stddev.bounds = (2,12)
                g_init.mean.bounds = (xcenr[line+0]-0.6, xcenr[line+0]+0.6)
                g_init.amplitude.bounds = (0.99*mcenr[line+0], 1.01*mcenr[line+0])

                #We only fit the 7 pixels centered on the center. This gives best fit. 
                g = fit_g(g_init, range(xcenr[line]-3, xcenr[line]+4), smoothed[xcenr[line]-3: xcenr[line]+4])
                #plt.plot(range(4096), g(range(4096)))

                gaus_means.append(g.mean.value)
                #now gaus_means is a better version of xcenr. We want to feed it 
                #into order line lists
            gaus_orders.append(gaus_means)

        #Now match gaus_orders to pixels for each column
        all_good_peaks = []
        for col_num in range(int(max_wid/column_width)-1):
    
            gaus_order = gaus_orders[col_num]
            ord_pixels = pixels[num]
            good_peaks = []
    
            #loop through gaussian peaks
            for peak in gaus_order:
                for line in ord_pixels:
                    dist = peak - line
                    if np.abs(dist) < 2.5:
                
                        #hold good, revevant peaks
                        good_peaks.append(peak)
                
            #hold all good peaks, all orders            
            all_good_peaks.append(good_peaks)
    
        #make sure all have the right number of lines    
        a = [len(all_good_peaks[col]) for col in range(int(max_wid/column_width)-1)]
        while len(set(a)) > 1:
            print "Looks like we have a bit of a fitting problem... I'll try and fix it."
    
            #Find peaks that appear in all orders
            useless = []
            for col in range(len(a)):
                for line in range(len(all_good_peaks[col])):
                    dists = dict()
            
                    for col1 in range(len(a)):
                        dists[col1] = []
                
                        for line1 in range(len(all_good_peaks[col1])):
                            dists[col1].append(np.abs(all_good_peaks[col1][line1]-all_good_peaks[col][line]))
                    
                        dists[col1] = np.min(dists[col1])
                
                    dist = dists.values()
                    mindist = np.max(dist)
            
                    if mindist > 2.5:
                        useless.append((col, line))
    
            #reformat
            bad_cols = []
            bad_lines = []
            for line in range(len(useless)):
                bad_cols.append(useless[line][0])
                bad_lines.append(useless[line][1])
                bad_lines.sort()
                bad_cols.sort()
                
            #now remove the offensive items, looping backwards    
            for point in range(len(useless)-1,-1,-1):
                all_good_peaks[useless[point][0]].pop(useless[point][1])
            a = [len(all_good_peaks[col]) for col in range(int(max_wid/column_width)-1)]
        print 'number of lines detected in each column (should all be the same)\n', a
            #god this piece was a pain to write. Hallelujah. 

        #Now find the stupid broken lines.
        #define maximum line offset allowed before fixing: TODO make this less arbitrary
        if num == 0 or num == 1:
            threshold = 1.8
        elif num == 2 or num == 3: 
            thershold = 1.5
        elif num == 4 or num == 5:
            threshold = 1.0
        elif num == 6 or num == 7:
            threshold = 0.9
        elif num == 8:
            threshold = 0.7
        else:
            threshold = 1.0
        
        badlines = []
        avgs = []
        for peak in range(len(all_good_peaks[0])):
            cens = [all_good_peaks[col][peak] for col in range(int(max_wid/column_width)-1)]
            maxx = np.max(cens) - np.min(cens)
            if maxx > threshold: #just for this order
                badlines.append(peak)
                avgs.append(np.mean(cens))
   
        for line in range(len(badlines)):
            #find nearest good line
            nearests = []
            for peak in range(0,5): #loop through nearby lines

                if badlines[line]+peak not in badlines and 0 < badlines[line]+peak < len(all_good_peaks[0]):
                    nearests.append(badlines[line]+peak)
    
            if len(nearests) == 0:
                for peak in range(-5,0):
                    if badlines[line]+peak not in badlines and 0 < badlines[line]+peak < len(all_good_peaks[0]):
                        nearests.append(badlines[line]+peak)
        
            if len(nearests) == 0:
                for peak in range(5,10): 
                    if badlines[line]+peak not in badlines and 0 < badlines[line]+peak < len(all_good_peaks[0]):
                        nearests.append(badlines[line]+peak)
        
            if len(nearests) == 0:
                for peak in range(-10, -5):
                    if badlines[line]+peak not in badlines and 0 < badlines[line]+peak < len(all_good_peaks[0]):
                        nearests.append(badlines[line]+peak)

                
            nearest = nearests[0]        

            last_peak = [all_good_peaks[col][nearest] for col in range(int(max_wid/column_width)-1)]
            last_slope = last_peak - np.mean(last_peak)

            for col in range(int(max_wid/column_width)-1):
                all_good_peaks[col][badlines[line]] = avgs[line]+last_slope[col]
    
    
        #cool. This seems to fix the ugly lines
        #now make an array to hold all the peaks in their locations:
        peak_coords = []
        for col in range(len(all_good_peaks)):
            for line in range(len(all_good_peaks[0])):
                peak_coords.append(((col+0.5)*column_width, all_good_peaks[col][line]))
        peak_coords = np.array(peak_coords)
    
        x_coords = peak_coords[:,0]
        y_coords = peak_coords[:,1]

        #we still need to transform these x coords into real x coords
        x_coords = order_left(y_coords) + x_coords
        
        if num == 9: #Because I'm using sky lines for the reddest order
            
            sky = pyfits.getdata(str(date)+'/Calibs/reduced/' + str(sky_id) + '_mean.fits')
        
            column_width = 20

            order_left = esiorders[num].xl
            order_right = esiorders[num].xr

            width = order_right(range(4096))-order_left(range(4096))
            width = np.abs(width)
            max_wid = np.max(width)

            #make masks at 5 pixel intervals
            lx = 2045
            ly = 4096
            X, Y = np.ogrid[0:ly, 0:lx]

            gaus_orders = []
            for line in range(int(max_wid/column_width)-1): #leave off the bad part
                print "column " + str(line) + '...'
                mask_left = order_left(X) + column_width*line < Y
                mask_left5 = order_left(X) + column_width*(line+1)  > Y
                mask_right = order_right(X) > Y

                mask = mask_left*mask_left5*mask_right

                #Now work down mask (5 pixels wide) and add into a 1d spectrum
                #load order to otherwise blank array
                blank = np.zeros((4096, 2045))
                blank[mask] = lamp[mask]

                #average over x pixels
                ord_tot = []
                for line in range(4096):
                    row_avg = np.mean(blank[line])*max_wid #normalized to before
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
                    g_init.stddev.bounds = (2,12)
                    g_init.mean.bounds = (xcenr[line+0]-0.6, xcenr[line+0]+0.6)
                    g_init.amplitude.bounds = (0.99*mcenr[line+0], 1.01*mcenr[line+0])

                    #We only fit the 7 pixels centered on the center. This gives best fit. 
                    g = fit_g(g_init, range(xcenr[line]-3, xcenr[line]+4), smoothed[xcenr[line]-3: xcenr[line]+4])
                    #plt.plot(range(4096), g(range(4096)))

                    gaus_means.append(g.mean.value)
                    #now gaus_means is a better version of xcenr. We want to feed it 
                    #into order line lists
                gaus_orders.append(gaus_means)

            #Now match gaus_orders to pixels for each column
            all_good_peaks_sky = []
            for col_num in range(int(max_wid/column_width)-1):
    
                gaus_order = gaus_orders[col_num]
                ord_pixels = pixels[num]
                good_peaks = []
    
                #loop through gaussian peaks
                for peak in gaus_order:
                    for line in ord_pixels:
                        dist = peak - line
                        if np.abs(dist) < 2.5:
                
                            #hold good, revevant peaks
                            good_peaks.append(peak)
                
                #hold all good peaks, all orders            
                all_good_peaks_sky.append(good_peaks)
    
            #make sure all have the right number of lines    
            a = [len(all_good_peaks_sky[col]) for col in range(int(max_wid/column_width)-1)]
            while len(set(a)) > 1:
                print "Looks like we have a bit of a fitting problem... I'll try and fix it."
    
                #Find peaks that appear in all orders
                useless = []
                for col in range(len(a)):
                    for line in range(len(all_good_peaks_sky[col])):
                        dists = dict()
            
                        for col1 in range(len(a)):
                            dists[col1] = []
                
                            for line1 in range(len(all_good_peaks_sky[col1])):
                                dists[col1].append(np.abs(all_good_peaks_sky[col1][line1]-all_good_peaks_sky[col][line]))
                    
                            dists[col1] = np.min(dists[col1])
                
                        dist = dists.values()
                        mindist = np.max(dist)
            
                        if mindist > 2.5:
                            useless.append((col, line))
    
                #reformat
                bad_cols = []
                bad_lines = []
                for line in range(len(useless)):
                    bad_cols.append(useless[line][0])
                    bad_lines.append(useless[line][1])
                    bad_lines.sort()
                    bad_cols.sort()
                
                #now remove the offensive items, looping backwards    
                for point in range(len(useless)-1,-1,-1):
                    all_good_peaks_sky[useless[point][0]].pop(useless[point][1])
                a = [len(all_good_peaks_sky[col]) for col in range(int(max_wid/column_width)-1)]
            print 'number of lines detected in each column (should all be the same)\n', a
                #god this piece was a pain to write. Hallelujah. 

            #Now find the stupid broken lines.
            #define maximum line offset allowed before fixing: TODO make this less arbitrary
            if num == 0 or num == 1:
                threshold = 1.8
            elif num == 2 or num == 3: 
                thershold = 1.5
            elif num == 4 or num == 5:
                threshold = 1.0
            elif num == 6 or num == 7:
                threshold = 0.9
            elif num == 8:
                threshold = 0.7
            else:
                threshold = 1.0
        
            badlines = []
            avgs = []
            for peak in range(len(all_good_peaks_sky[0])):
                cens = [all_good_peaks_sky[col][peak] for col in range(int(max_wid/column_width)-1)]
                maxx = np.max(cens) - np.min(cens)
                if maxx > threshold: #just for this order
                    badlines.append(peak)
                    avgs.append(np.mean(cens))
   
            for line in range(len(badlines)):
                #find nearest good line
                nearests = []
                for peak in range(0,5): #loop through nearby lines

                    if badlines[line]+peak not in badlines and 0 < badlines[line]+peak < len(all_good_peaks[0]):
                        nearests.append(badlines[line]+peak)
    
                if len(nearests) == 0:
                    for peak in range(-5,0):
                        if badlines[line]+peak not in badlines and 0 < badlines[line]+peak < len(all_good_peaks[0]):
                            nearests.append(badlines[line]+peak)
        
                if len(nearests) == 0:
                    for peak in range(5,10): 
                        if badlines[line]+peak not in badlines and 0 < badlines[line]+peak < len(all_good_peaks[0]):
                            nearests.append(badlines[line]+peak)
        
                if len(nearests) == 0:
                    for peak in range(-10, -5):
                        if badlines[line]+peak not in badlines and 0 < badlines[line]+peak < len(all_good_peaks[0]):
                            nearests.append(badlines[line]+peak)

                
                nearest = nearests[0]        

                last_peak = [all_good_peaks_sky[col][nearest] for col in range(int(max_wid/column_width)-1)]
                last_slope = last_peak - np.mean(last_peak)

                for col in range(int(max_wid/column_width)-1):
                    all_good_peaks_sky[col][badlines[line]] = avgs[line]+last_slope[col]
    
    
            #cool. This seems to fix the ugly lines
            #now make an array to hold all the peaks in their locations:
            peak_coords = []
            for col in range(len(all_good_peaks_sky)):
                for line in range(len(all_good_peaks[0])):
                    peak_coords.append(((col+0.5)*column_width, all_good_peaks_sky[col][line]))
            peak_coords = np.array(peak_coords)
    
            x_coords_sky = peak_coords[:,0]
            y_coords_sky = peak_coords[:,1]

            #we still need to transform these x coords into real x coords
            x_coords_sky = order_left(y_coords_sky) + x_coords_sky
            
            x_coords = np.array(list(x_coords) + list(x_coords_sky))
            y_coords = np.array(list(y_coords) + list(y_coords_sky))
            all_good_peaks = np.array(list(all_good_peaks) + list(all_good_peaks_sky))
            #END IF(ORDER == 9)
            

        #now find the wavelengths corresponding to these
        g_wavelengths = []
        for line in range(len(wavelengths[num])):
            dist = [np.abs(pixels[num][line] - all_good_peaks[0][peak]) for peak in range(len(all_good_peaks[0]))]
            dist = np.min(dist)
    
            if dist < 2.5:
                g_wavelengths.append(wavelengths[num][line])
        
        wavelengths2d = []
        for col in range(len(all_good_peaks)):
            for line in range(len(all_good_peaks[0])):
                wavelengths2d.append(g_wavelengths[line])
        
        p_init = models.Polynomial2D(degree=poly_fit_ord)
        fit_p = fitting.LinearLSQFitter() #possibly should switch to LinearLSQFitter()
        p = fit_p(p_init, x_coords, y_coords, wavelengths2d)
        pickle_obj.append(p)
    
    pickle.dump(pickle_obj, open(str(date)+'/Calibs/solution2d_'+str(date)+'.p', 'wb'))
        
    

        