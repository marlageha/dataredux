##########################################################
#       Finds the rough emission lines in lamp           #
#       Works order by order and lamp by lamp            #
#       Saves y coordinates of each emission line        #
#       and corresponding intensity                      #
#       In pickle files as lamp[order[0]] = ycords       #
#                          lamp[order[1]] = intensity    #
#                          lamp[order[2]] = spectrum     #
#                          lamp[order[3]] = smoothed     #
#                                                        #
#       writes to Calibs/lamp_peaks/obj_id_roughpeaks    #
#                                                        #
#       Kareem El-Badry, 07/14/2014                      #
##########################################################
from __future__ import division
import numpy as np
import os
import pickle
import pyfits
import scipy
from scipy.ndimage.filters import gaussian_filter1d

def esi_roughpeaks(date):
    
    #Get edge masks from file
    print 'loading masks...'
    orders_mask = pickle.load(open(str(date)+'/Calibs/orders_mask_'+str(date)+'.p', 'rb'))
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

    #Read in lamps:
    for obj_id in lines:
        
        print "reducing " + str(obj_id) + '...'
        
        #load file
        lamp = pyfits.getdata(str(date)+'/Calibs/reduced/' + str(obj_id) + '_mean.fits')
        sky = pyfits.getdata(str(date)+'/Calibs/reduced/' + str(objects[0]) + '_mean.fits')
        
        pickle_obj = []
        #loop through orders:
        for num in range(len(all_order_mask)):
                
            print 'working on order ' + str(num)
             
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
                badindices = np.where(resid > 3*std)[0]
                if badindices.size > 0:
                    done = False
                    x = np.delete(x, badindices)
                    y = np.delete(y, badindices)
                    
            #subtract baseline
            ord_tot = ord_tot - fit(range(4096))
            
            #smooth slightly, find slope
            smoothed = gaussian_filter1d(ord_tot, 5)
            
            smooted = gaussian_filter1d(ord_tot, 1.5)
            deriv = scipy.ndimage.sobel(smooted)#mayber deriv of smoothed?
            
            #find local maxes. Could also just use the scipy max finder, but this is faster. 
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
            
            if num == 9: #NEED SKYLINES FOR LAST ORDER
                
                #load order to otherwise blank array
                blank = np.zeros((4096, 2045))
                blank[all_order_mask[num]] = sky[all_order_mask[num]]
            
                #average over x pixels
                ord_tot_sky = []
                for line in range(4096):
                    row_avg = np.mean(blank[line])
                    ord_tot_sky.append(row_avg)
                ord_tot_sky = np.array(ord_tot_sky)
            
                #fit a polynomial to baseline with iterative fit
                x = range(4096)
                y = ord_tot_sky

                done = False
                while not done:
                    done = True
                    fit = np.poly1d(np.polyfit(x, y, 10))
                    resid = y - fit(x)
                    std = np.std(resid)
                    badindices = np.where(resid > 3*std)[0]
                    if badindices.size > 0:
                        done = False
                        x = np.delete(x, badindices)
                        y = np.delete(y, badindices)
                    
                #subtract baseline
                ord_tot_sky = ord_tot_sky - fit(range(4096))
            
                #smooth slightly, find slope
                smoothed_sky = gaussian_filter1d(ord_tot_sky, 5)
            
                smooted_sky = gaussian_filter1d(ord_tot_sky, 1.5)
                deriv_sky = scipy.ndimage.sobel(smooted_sky)#mayber deriv of smoothed?
            
                #find local maxes. Could also just use the scipy max finder, but this is faster. 
                mask = np.r_[True, smoothed_sky[1:] > smoothed_sky[:-1]] & np.r_[smoothed_sky[:-1] > smoothed_sky[1:], True]
                xcenr_sky = np.arange(4096)[mask]
                mcenr_sky = smoothed_sky[mask]
            
                #makes easier to work with for present
                xcenr_sky = list(xcenr_sky)
                mcenr_sky = list(mcenr_sky)
            
                #Make list of list of slopes before and after each peak:
                slopes_after = []
                slopes_before = []
                mins = []
                maxes = []

                for line in range(len(xcenr_sky)):
                    try: 
                        slopes_before.append([deriv_sky[xcord] for xcord in range(xcenr_sky[line]-5, xcenr_sky[line])])
                        slopes_after.append([deriv_sky[xcord] for xcord in range(xcenr_sky[line]+1, xcenr_sky[line]+6)])
                    except IndexError:
                        slopes_after.append(deriv_sky[xcenr_sky]) #if its at the edge; no good anyway
                        slopes_before.append(deriv_sky[xcenr_sky])
        
                    mins.append(np.min(slopes_after[line]))
                    maxes.append(np.max(slopes_before[line]))
    
                mins = np.array(mins)
                maxes = np.array(maxes)
    
                #slopes at least steep enough on either side
                good_max = [line for line in range(len(xcenr_sky)) if mins[line] < -0.001 and maxes[line] > 0.001]
                xcenr_sky = np.array(xcenr_sky)
                mcenr_sky = np.array(mcenr_sky)

                xcenr_sky = xcenr_sky[good_max]
                mcenr_sky = mcenr_sky[good_max]
            
                #remove peaks at ends of orders:
                mcenr_sky = mcenr_sky[(xcenr_sky > 5) & (xcenr_sky < 4090)] #important that mcenr is first!
                xcenr_sky = xcenr_sky[(xcenr_sky > 5) & (xcenr_sky < 4090)] #otherwise the first line changes the second.
            
                xcenr = list(xcenr) + list(xcenr_sky) 
                mcenr = list(mcenr) + list(mcenr_sky) 
                ord_tot = list(ord_tot) + list(ord_tot_sky) 
                smoothed = list(smoothed) + list(smoothed_sky) 
                
                xcenr.sort()

            #for each order
            pickle_order = [xcenr, mcenr, ord_tot, smoothed]
            
            #for entire lamp
            pickle_obj.append(pickle_order)
            
        #Make directory to hold decosmicified files
        if not os.path.exists(str(date)+'/Calibs/lamp_peaks/'):
            os.makedirs(str(date)+'/Calibs/lamp_peaks/')
            
        #write to file
        pickle.dump(pickle_obj, open(str(date) + '/Calibs/lamp_peaks/'+ str(obj_id) + '_roughpeaks_'+str(date)+'.p', 'wb'))
        
        # TODO: for reddest order, find sky lines (from a science exposure) to get a better 
        # wavelength solution
        
        
        
            
            
        