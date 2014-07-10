######################################################
#       Applies, bias, flats, median combines        #
#       Science data and lamps                       #
#       Makes variance imagees for each science      #
#       Saves 1/(sigma)^2 to Calibs/variance         #
#                                                    #
#       Kareem El-Badry, 07/10/2014                  #
######################################################

from __future__ import division
import numpy as np
import pickle
import pyfits
import scipy.ndimage 

def esi_reduce():

    #Get edge masks from file
    orders_mask = pickle.load(open('Calibs/orders_mask.p', 'rb'))
    background_mask = pickle.load(open('Calibs/background_mask.p', 'rb'))

    #Bias
    bias = pyfits.getdata('Calibs/bias.fits')[:, 25:2070]

    #Normalized Flat
    flat = pyfits.getdata('Calibs/norm_flat.fits')

    #READ LOG
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

    #Find list of objects:
    names = []
    for line in range(len(alldata)):
        if "Object" in alldata[line][3] and float(alldata[line][6]) > 600:
            names.append(alldata[line][2])
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place 

    #Find list of lamps/lines
    lines = []
    for line in range(len(alldata)):
        if "Line" in alldata[line][3]:
            lines.append(alldata[line][2])
    lamps = np.array(list(set(lines)))
    lamps.sort()

    #Reduce for each lamp:
    for lamp in lamps:
        print "reducing " + str(lamp) + '...'
    
        #Find files related to lamp:
        alamp = []
        for line in range(len(good)):
            if str(lamp) in good[line][2]:
                alamp.append(good[line])
        print str(len(alamp)), "files for "+ str(lamp)
    
        #Write path to each lamp's files:
        lamp_locs = ['Raw/' + str(alamp[line][0]) for line in range(len(alamp))]
    
        #Read in data
        all_lamp = []
        for line in lamp_locs:
            sci = pyfits.getdata(line)[:, 25:2070]
            sci = (sci - bias)/flat
            all_lamp.append(sci)
        
        #Median Combine
        med = np.median(all_lamp, axis = 0)
        fits = pyfits.PrimaryHDU(med)
        fits.writeto('Calibs/reduced/'+str(lamp)+'_med.fits', clobber = True)

    #Reduce for each object
    for obj_id in objects:
    
        print "reducing " + str(obj_id) + '...'
    
        #Find files related to obj_id:
        aobj_id = []
        for line in range(len(good)):
            if str(obj_id) in good[line][2]:
                aobj_id.append(good[line])
        print str(len(aobj_id)),"files for " + str(obj_id)

        #Write Path to each object's files
        obj_locs = ["Raw/" + str(aobj_id[line][0]) for line in range(len(aobj_id))]

        #Read in Data
        all_sci = []
        for line in obj_locs:
            sci = pyfits.getdata(line)[:, 25:2070]
            sci = (sci - bias)/flat
            all_sci.append(sci)
    
        #Median Combine
        med = np.median(all_sci, axis = 0)
        fits = pyfits.PrimaryHDU(med)
        fits.writeto('Calibs/reduced/'+str(obj_id)+'_med.fits', clobber = True)
    



    #Make variance image for each science image:
    #The break in the amplifiers is at 1022-1023. Everything increases at 1023. 

    for obj_id in objects:
    
        print "making variance for " + str(obj_id) + '...'
    
        #Find files related to obj_id:
        aobj_id = []
        for line in range(len(good)):
            if str(obj_id) in good[line][2]:
                aobj_id.append(good[line])
        print str(len(aobj_id)),"files for " + str(obj_id)

        #Write Path to each object's files
        obj_locs = ["Raw/" + str(aobj_id[line][0]) for line in range(len(aobj_id))]
    
        #CALCULATE EACH NOISE CONTRIBUTION SEPARATELY
    
        #Read in Data
        all_var = []
        for line in range(len(obj_locs)):
            im = pyfits.getdata(obj_locs[line])[:, 25:2070]
        
            rn = np.zeros((4096, 2045))
            #make gain mask: 
            right_mask = np.empty(flat.shape,dtype=bool)
            right_mask[:, :1022] = 0
            right_mask[:, 1023:] = 1
            left_mask = -right_mask

            #MAKE A READ NOISE FRAME
            #random bias image
            rand_bias = pyfits.getdata('Raw/e140306_0012.fits.gz')[:, 25:2070]

            read_noise = rand_bias - bias

            left_rn = np.std(read_noise[2000:2100, 920:1020])
            right_rn = np.std(read_noise[2000:2100, 1030:1130])

            rn[left_mask] = left_rn
            rn[right_mask] = right_rn
            #as expected, the read nose is a bit higher on the right side. 

            #GAIN AND POISSON NOISE
            #gain on right side == 1.29 e-/DN
            poisson = np.zeros((4096, 2045))
            left_bias = np.mean(bias[2000:2100, 920:1020])
            right_bias = np.mean(bias[2000:2100, 1030:1130])

            left_gain = 1.29

            right_gain = left_gain*left_bias/right_bias
            #as expected, gain a little lower on right side
        
            poisson[left_mask] = np.sqrt(im[left_mask]/left_gain)
            poisson[right_mask] = np.sqrt(im[right_mask]/right_gain)
        
            noise = np.sqrt(rn**2+poisson**2)
        
            variance = 1/noise**2
        
            #write to file
            fits = pyfits.PrimaryHDU(variance)
            fits.writeto('Calibs/variance/'+str(obj_id)+'_'+str(line + 1)+'_var.fits', clobber = True)
    