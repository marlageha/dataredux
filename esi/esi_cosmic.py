##############################################################################
#                                                                            #
#       Reads in all Raw Science and Lamps                                   #
#       Uses LA Cosmics (van Dokkum)                                         #
#       All credit goes to Pieter van Dokkum for writing the original        #
#       IRAF script and Malte Tewes for translating it to Python             #
#       Before running, save cosmics.py in PYTHONPATH directory              #
#       Makes directory Calibs/cosmicless; this writes to there              #
#       Runs in about 2 minutes per image. Can take a while.                 #
#       Kareem El-Badry, 07/10/2014                                          #
#                                                                            #
##############################################################################


from __future__ import division
import cosmics
import numpy as np
import os
import pickle
import pyfits
import scipy


def esi_cosmic(date):
    

    #Get edge masks from file
    orders_mask = pickle.load(open(str(date)+'/Calibs/orders_mask_'+str(date)+'.p', 'rb'))
    background_mask = pickle.load(open(str(date)+'/Calibs/background_mask_'+str(date)+'.p', 'rb'))

    #Bias
    bias = pyfits.getdata(str(date)+'/Calibs/bias_'+str(date)+'.fits')

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
       if "yes" in alldata[line][ 7]:
           good.append(alldata[line])
       
    #Find list of objects:
    names = []
    for line in range(len(alldata)):
        if ("Object" in alldata[line][3] and float(alldata[line][6]) > 600) or "*" in alldata[line][2]:
            names.append(alldata[line][2])
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place 
    
    #Find list of lines
    names = []
    for line in range(len(alldata)):
        if "Line" in alldata[line][3]:
            names.append(alldata[line][2])
    lines = np.array(list(set(names)))
    lines.sort() #ascending order, modify in place
    
    #Find list of stars:
    names = []
    for line in range(len(alldata)):
        if "*" in alldata[line][2]:
            names.append(alldata[line][2])
    stars = np.array(list(set(names)))
    stars.sort() #ascending order, modify in place 
    
    #Make directory to hold decosmicified files
    if not os.path.exists(str(date)+'/Calibs/cosmicless/'):
        os.makedirs(str(date)+'/Calibs/cosmicless/')
    
    #Remove cosmics from each lamp
    for obj_id in lines:
        print "reducing " + str(obj_id) + '...'
    
        #Find files related to this object
        aobj_id = []
        for line in range(len(good)):
            if str(obj_id) in good[line][2]:
                aobj_id.append(good[line])
        print str(len(aobj_id)), "files for "+ str(obj_id)
    
        #Write path to each objects files:
        obj_locs = [str(date)+'/Raw/' + str(aobj_id[line][0]) for line in range(len(aobj_id))]
    
        #Read in and de-cosmify
        for line in range(len(obj_locs)):
        
            array, header = cosmics.fromfits(obj_locs[line])
            #array = array - bias #backwards for some reason
            c = cosmics.cosmicsimage(array, gain = 1.29, readnoise = 2.2, sigclip = 6, objlim = 5.0, sigfrac = 0.7, satlevel = 1e4)
            c.run(maxiter = 3) # can increase up to 4 to improve precision, but takes longer
            cosmics.tofits(str(date)+'/Calibs/cosmicless/' + str(obj_id)+'_'+str(line+1)+'decos.fits', c.cleanarray, header)
            
            #now zip it (to save space)
            f = pyfits.getdata(str(date)+'/Calibs/cosmicless/' + str(obj_id)+'_'+str(line+1)+'decos.fits')
            hdu = pyfits.CompImageHDU(f)
            hdu.writeto(str(date)+'/Calibs/cosmicless/' + str(obj_id)+'_'+str(line+1)+'decos.fits', clobber=True)

    #Remove cosmics from each object
    for obj_id in objects:
        print "reducing " + str(obj_id) + '...'
    
        #Find files related to this object
        aobj_id = []
        for line in range(len(good)):
            if str(obj_id) in good[line][2]:
                aobj_id.append(good[line])
        print str(len(aobj_id)), "files for "+ str(obj_id)
    
        #Write path to each objects files:
        obj_locs = [str(date)+'/Raw/' + str(aobj_id[line][0]) for line in range(len(aobj_id))]
    
        #Read in and de-cosmify
        for line in range(len(obj_locs)):
        
            array, header = cosmics.fromfits(obj_locs[line])
            #array = array - bias #backwards for some reason
            c = cosmics.cosmicsimage(array, gain = 1.29, readnoise = 2.2, sigclip = 3.8, objlim = 3.0, sigfrac = 0.7, satlevel = -1)
            c.run(maxiter = 3) # can increase up to 4 to improve precision, but takes longer
            cosmics.tofits(str(date)+'/Calibs/cosmicless/' + str(obj_id)+'_'+str(line+1)+'decos.fits', c.cleanarray, header)
            
            #now zip it (to save space)
            f = pyfits.getdata(str(date)+'/Calibs/cosmicless/' + str(obj_id)+'_'+str(line+1)+'decos.fits')
            hdu = pyfits.CompImageHDU(f)
            hdu.writeto(str(date)+'/Calibs/cosmicless/' + str(obj_id)+'_'+str(line+1)+'decos.fits', clobber=True)
            
    #Remove cosmics from each star
    for obj_id in stars:
        print "reducing " + str(obj_id) + '...'
    
        #Find files related to this object
        aobj_id = []
        for line in range(len(good)):
            if str(obj_id) == good[line][2]:
                aobj_id.append(good[line])
        print str(len(aobj_id)), "files for "+ str(obj_id)
    
        #Write path to each objects files:
        obj_locs = [str(date)+'/Raw/' + str(aobj_id[line][0]) for line in range(len(aobj_id))]
    
        #Read in and de-cosmify
        for line in range(len(obj_locs)):
        
            array, header = cosmics.fromfits(obj_locs[line])
            #array = array - bias #backwards for some reason
            c = cosmics.cosmicsimage(array, gain = 1.29, readnoise = 2.2, sigclip = 10, objlim = 7.0, sigfrac = 0.7, satlevel = 15000)
            c.run(maxiter = 3) # can increase up to 4 to improve precision, but takes longer
            cosmics.tofits(str(date)+'/Calibs/cosmicless/' + str(obj_id)+'_'+str(line+1)+'decos.fits', c.cleanarray, header)
            
            #now zip it (to save space)
            f = pyfits.getdata(str(date)+'/Calibs/cosmicless/' + str(obj_id)+'_'+str(line+1)+'decos.fits')
            hdu = pyfits.CompImageHDU(f)
            hdu.writeto(str(date)+'/Calibs/cosmicless/' + str(obj_id)+'_'+str(line+1)+'decos.fits', clobber=True)
                
