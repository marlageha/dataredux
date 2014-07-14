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

import matplotlib.pyplot as plt
import numpy as np
import pickle
import pyfits
import scipy
from scipy.ndimage.filters import gaussian_filter1d

def esi_roughpeaks():
    
    #Get edge masks from file
    print 'loading masks...'
    orders_mask = pickle.load(open('Calibs/orders_mask.p', 'rb'))
    background_mask = pickle.load(open('Calibs/background_mask.p', 'rb')) 
    all_order_mask = pickle.load(open('Calibs/all_order_masks.p', 'rb'))
    
    #get names of lamps, first read in log:
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

    #Find list of objects 
    names = []
    for line in range(len(alldata)):
        if "Line" in alldata[line][3]:
            names.append(alldata[line][2])
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place

    #Read in lamps:
    for obj_id in objects:
        
        print "reducing " + str(obj_id) + '...'
        
        #load file
        lamp = pyfits.getdata('Calibs/reduced/' + str(obj_id) + '_mean.fits')
        
        pickle_obj = []
        #loop through orders:
        for num in range(len(all_order_mask)):
            
            print 'working on order ' + str(num + 1)
             
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
            smoothed = gaussian_filter1d(ord_tot, 3)
            deriv = scipy.ndimage.sobel(ord_tot)#mayber deriv of smoothed?
            
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
            good_max = [line for line in range(len(xcenr)) if mins[line] < -0.004 and maxes[line] > 0.004]
            xcenr = np.array(xcenr)
            mcenr = np.array(mcenr)

            xcenr = xcenr[good_max]
            mcenr = mcenr[good_max]
            
            #remove peaks at ends of orders:
            mcenr = mcenr[(xcenr > 5) & (xcenr < 4090)] #important that mcenr is first!
            xcenr = xcenr[(xcenr > 5) & (xcenr < 4090)] #otherwise the first line changes the second. 
            
            #for each order
            pickle_order = [xcenr, mcenr, ord_tot, smoothed]
            
            #for entire lamp
            pickle_obj.append(pickle_order)
            
        #write to file
        pickle.dump(pickle_obj, open('Calibs/lamp_peaks/'+ str(obj_id) + '_roughpeaks.p', 'wb'))
        
            
            
        