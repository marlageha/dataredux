##############################################################################
#                                                                            #
#    Subtracts sky, defined according to pixels [3:15] per line one each     #
#    side of each order and subtracted line by line. For each order, a sky   #
#    value is defined for each line as the median of the 2*(15-3) pixels     #
#    on the edges of the order. A sky image is then made with this value in  #
#    each line of each order. That sky image is then smoothed in the y       #
#    direction with a length 5 median filter.                                #
#                                                                            #
#    Mask for sky and orders are read from pickle files which are written    #
#    in esi_traceflat.py                                                     #
#                                                                            #
#    Sky subtracted images are written to Calibs/sky_sub.                    #
#                                                                            #
#    Kareem El-Badry, 07/17/2014                                             #
#                                                                            #
##############################################################################

from __future__ import division
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pyfits
import scipy.ndimage
import pdb

def esi_skysub():

    #Get edge masks from file
    print 'loading masks...'
    orders_mask = pickle.load(open('Calibs/orders_mask.p', 'rb'))
    all_order_mask = pickle.load(open('Calibs/all_order_masks.p', 'rb'))
    sky_mask = pickle.load(open('Calibs/sky_mask.p', 'rb'))
    

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
    for line in range(len(good)):
        if "Object" in good[line][3] and float(good[line][6]) > 600:
            names.append(good[line][2])
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place

    '''
    sky_objs = ['81315', '113209', '119887', '36363', 'Beowulf1',
                '37836', '50778'] #dim objects that are mostly sky
                
    skys = []
    for obj_id in sky_objs:
        
        img = pyfits.getdata('Calibs/reduced/'+str(obj_id)+'_mean.fits')
        skys.append(img)
        
    rough_sky = np.median(skys, axis = 0)
    fits = pyfits.PrimaryHDU(rough_sky)
    fits.writeto('Calibs/skytest.fits', clobber = True)'''
    
    
    rough_sky = pyfits.getdata('Calibs/skytest.fits')
    rough_sky = scipy.ndimage.median_filter(rough_sky, size = (1,5))#smooth a bit
    #loop through objects 
    for obj_id in objects:
        print "sky subtracting " +str(obj_id) + '...'

        #Get object
        image = pyfits.getdata('Calibs/reduced/'+str(obj_id)+'_mean.fits')
        #image = image - rough_sky
        
        #start with a blank array; fill it order by order
        sky = np.zeros((4096, 2045))
        for num in range(len(sky_mask)):

            #load order to otherwise blank array
            blank = np.zeros((4096, 2045))
            blank[sky_mask[num]] = image[sky_mask[num]] 

            #average over x pixels
            ord_tot = []
            for line in range(4096):
                row = np.array(blank[line])
                nonzero = row[np.where(row != 0)[0]]
                row_avg = np.median(nonzero)
                ord_tot.append(row_avg)
            ord_tot = np.array(ord_tot)
    
            #make sky image
            dummy = np.zeros((4096, 2045))
            dummy[all_order_mask[num]] = image[all_order_mask[num]]
    
            for line in range(4096):
                dummy[line] = ord_tot[line]
    
            sky[all_order_mask[num]] = dummy[all_order_mask[num]]
    
        sky = scipy.ndimage.median_filter(sky, size = (3,1))
    
        image[orders_mask] = (image - sky)[orders_mask]
        #image[-orders_mask] = 0
        
        
        #Get rid of pesky dead column nonsense:
        mask1 = np.abs(image) > 2
        mask2 = np.zeros((4096, 2045), dtype = 'bool')
        mask2[2645:4096, 410:450] = True
        mask = mask1*mask2
        #pickle.dump(mask, open('Calibs/bad_pix.p', 'wb'))
        #pdb.set_trace()
        '''
        print "removing bad pixels..."
        x_zero = np.where(image < -5)[0]
        y_zero = np.where(image < -5)[1]

        for line in range(len(x_zero)):
            
            #ok, but slow
            adjacent = [(x_zero[line], y_zero[line]+1), (x_zero[line], y_zero[line]-1), 
                        (x_zero[line]+1, y_zero[line]), (x_zero[line]-1, y_zero[line]),
                        (x_zero[line]+1, y_zero[line]+1), (x_zero[line]-1, y_zero[line]-1),
                        (x_zero[line]+1, y_zero[line]-1), (x_zero[line]-1, y_zero[line]+1),
                        (x_zero[line], y_zero[line]+2), (x_zero[line], y_zero[line]-2), 
                        (x_zero[line]+2, y_zero[line]), (x_zero[line]-2, y_zero[line]),
                        (x_zero[line]+2, y_zero[line]+2), (x_zero[line]-2, y_zero[line]-2),
                        (x_zero[line]+2, y_zero[line]-2), (x_zero[line]-2, y_zero[line]+2)]
                                        
            try:
                adj_vals = [image[x] for x in adjacent if image[x] > -5]
        
            except IndexError:
                 #If it's at the edge and there's no adjacent:
                adj_vals = [np.median(image) for x in range(len(adjacent))]  
        
    
            image[x_zero[line], y_zero[line]] = np.median(adj_vals)
        '''
        image[mask] = np.median(image)
        
        fits = pyfits.PrimaryHDU(image)
        fits.writeto('Calibs/sky_sub/'+str(obj_id)+'_skysub.fits', clobber = True)
    
    
    
