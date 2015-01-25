######################################################
#       Applies, bias, flats, median combines        #
#       Science data and lamps                       #
#       Throws out remaining cosmic rays             #
#       Makes variance imagees for each science      #
#       Saves 1/(sigma)^2 to Calibs/variance         #
#                                                    #
#       Kareem El-Badry, 07/10/2014                  #
######################################################

from __future__ import division
from astropy.stats import sigma_clip
import numpy as np
import os
import pickle
import pyfits
import scipy.ndimage 
import pdb

def esi_reduce(date):

    #Get edge masks from file
    orders_mask = pickle.load(open(str(date)+'/Calibs/orders_mask_'+str(date)+'.p', 'rb'))
    background_mask = pickle.load(open(str(date)+'/Calibs/background_mask_'+str(date)+'.p', 'rb'))

    #Bias
    bias = pyfits.getdata(str(date)+'/Calibs/bias_'+str(date)+'.fits')[:, 25:2070]

    #Normalized Flat
    flat = pyfits.getdata(str(date)+'/Calibs/norm_flat_'+str(date)+'.fits')

    #READ LOG
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

    #Find list of objects 
    names = []
    for line in range(len(alldata)):
        if ("Object" in alldata[line][3] and float(alldata[line][6]) > 600):
            names.append(alldata[line][2])
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place 
    
    #Make directory to hold decosmicified files
    if not os.path.exists(str(date)+'/Calibs/reduced/'):
        os.makedirs(str(date)+'/Calibs/reduced/')

    #Reduce for each object 
    for obj_id in objects:

        print "reducing " + str(obj_id) + '...'

        #Find files related to obj_id:
        aobj_id = []
        for line in range(len(good)):
            if str(obj_id) in good[line][2]:
                aobj_id.append(good[line])
        print str(len(aobj_id)), "files for " + str(obj_id)
        
        if len(aobj_id) == 0:
            continue

        #Write Path to each object's files
        obj_locs = [str(date)+'/Calibs/cosmicless/' + str(obj_id) + '_' + str(line+1) + 
                    'decos.fits' for line in range(len(aobj_id))]

        #Read in Data
        all_sci = []
        for line in obj_locs:
            sci = pyfits.getdata(line)[:, 25:2070]
            all_sci.append(sci)
        
        #Get a mask for remaining cosmic rays; mean combine
        print "averaging with rejection..."
        filtered = sigma_clip(all_sci, axis=0, copy = False, varfunc = np.median, sig = 15)
        filtered = (filtered - bias)/flat
    
        mean = np.mean(filtered, axis = 0)
        mean.data[background_mask] = 0
        
        hdu = pyfits.CompImageHDU(mean.data)
        hdu.writeto(str(date)+'/Calibs/reduced/'+str(obj_id)+'_mean.fits', clobber = True)
        #pdb.set_trace()

    #Reduce lines
    names = []
    for line in range(len(alldata)):
        if ("Line" in alldata[line][3]) or ("*" in alldata[line][2]):
            names.append(alldata[line][2])
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place 


    #Reduce for each line or star (because not decosmfied)
    for obj_id in objects:

        print "reducing " + str(obj_id) + '...'

        #Find files related to obj_id:
        aobj_id = []
        for line in range(len(good)):
            if str(obj_id) in good[line][2]:
                aobj_id.append(good[line])
        print str(len(aobj_id)),"files for " + str(obj_id)
        
        if len(aobj_id) == 0:
            continue

        #Write Path to each object's files
        obj_locs = [str(date) + '/Raw/' + str(aobj_id[line][0]) for line in range(len(aobj_id))]

        #Read in Data
        all_sci = []
        for line in obj_locs:
            sci = pyfits.getdata(line)[:, 25:2070]
            sci = (sci - bias)/flat
            all_sci.append(sci)
        
        #median, not mean, for lines
        mean = np.median(all_sci, axis = 0)
        mean[background_mask] = 0    
        hdu = pyfits.CompImageHDU(mean)
        hdu.writeto(str(date)+'/Calibs/reduced/'+str(obj_id)+'_mean.fits', clobber = True)

    #Make variance image for each science image:
    #The break in the amplifiers is at 1022-1023. Everything increases at 1023. 
    
    names = []
    for line in range(len(alldata)):
        if ("Object" in alldata[line][3] and float(alldata[line][6]) > 600) or ("*" in alldata[line][2]):
            names.append(alldata[line][2])
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place 
    
    for obj_id in objects:

        print "making variance for " + str(obj_id) + '...'

        #Find files related to obj_id:
        aobj_id = []
        for line in range(len(good)):
            if str(obj_id) in good[line][2]:
                aobj_id.append(good[line])
        print str(len(aobj_id)),"files for " + str(obj_id)
        
        if len(aobj_id) == 0:
            continue

        #Write Path to each object's files
        obj_locs = [str(date)+'/Calibs/cosmicless/' + str(obj_id) + '_' + str(line+1) + 
                    'decos.fits' for line in range(len(aobj_id))]

        #CALCULATE EACH NOISE CONTRIBUTION SEPARATELY

        #Read in Data
        all_var = []
        for line in range(len(obj_locs)):
            im = pyfits.getdata(obj_locs[line])[:, 25:2070]
            mask = bias > im #to avoid taking sqrt(negative)
            im = im - bias
    
            rn = np.zeros((4096, 2045))
            #make gain mask: 
            right_mask = np.empty(flat.shape,dtype=bool)
            right_mask[:, :1022] = 0
            right_mask[:, 1023:] = 1
            left_mask = -right_mask

            #MAKE A READ NOISE FRAME
            #random bias image
            
            biases = []
            for line in range(len(good)):
                if "Bias" in good[line][3] and "Bias" in good[line][2]:
                    biases.append(good[line])
            rand_path = str(date)+'/Raw/'+str(biases[0][0]) #pick the first bias. 
            rand_bias = pyfits.getdata(rand_path)[:, 25:2070]
            
            read_noise = rand_bias - bias
            
            #Big enough square to be representative
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
            poisson[mask] = np.median(poisson)
    
            noise = np.sqrt(rn**2+poisson**2)
            noise = noise/flat
    
            #variance = 1/noise**2
            
            #noise = (noise - bias)/flat
            #noise[background_mask] = 0
            all_var.append(noise)
            #write to file
            #fits = pyfits.PrimaryHDU(variance)
            #fits.writeto('Calibs/variance/'+str(obj_id)+'_'+str(line + 1)+'_var.fits', clobber = True)
            
        #Make directory to hold decosmicified files
        if not os.path.exists(str(date)+'/Calibs/variance/'):
            os.makedirs(str(date)+'/Calibs/variance/')
            
        #Combine variance images for each object
        tot_noise = np.sqrt(np.sum([(siggma/len(all_var))**2 for siggma in all_var], axis = 0))
        tot_noise[background_mask] = 0
        fits = pyfits.CompImageHDU(tot_noise)
        fits.writeto(str(date)+'/Calibs/variance/'+str(obj_id)+'_noise.fits', clobber = True)
    
    