#######################################################################
#       Reads in all Raw Science and Lamps                            #
#       Uses LA Cosmics (van Dokkum)                                  #
#       Before running, save cosmics.py in PYTHONPATH directory       #
#       Make directory Calibs/cosmicless; this writes to that         #
#       Runs in about an hour for ~100 images                         #
#       Kareem El-Badry, 07/10/2014                                   #
#######################################################################


from __future__ import division
import cosmics
import matplotlib.pyplot as plt
import numpy as np
import pickle
import pyfits
import scipy

def esi_cosmic():

    #Get edge masks from file
    orders_mask = pickle.load(open('Calibs/orders_mask.p', 'rb'))
    background_mask = pickle.load(open('Calibs/background_mask.p', 'rb'))

    #Bias
    bias = pyfits.getdata('Calibs/bias.fits')

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
       if "yes" in alldata[line][ 7]:
           good.append(alldata[line])
       
    #Find list of objects and lines:
    names = []
    for line in range(len(alldata)):
        if "Object" in alldata[line][3] and float(alldata[line][6]) > 600 or "Line" in alldata[line][3]:
            names.append(alldata[line][2])
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place 

    #Remove cosmics from each objects/lamp
    for obj_id in objects:
        print "reducing " + str(obj_id) + '...'
    
        #Find files related to this object
        aobj_id = []
        for line in range(len(good)):
            if str(obj_id) in good[line][2]:
                aobj_id.append(good[line])
        print str(len(aobj_id)), "files for "+ str(obj_id)
    
        #Write path to each objects files:
        obj_locs = ['Raw/' + str(aobj_id[line][0]) for line in range(len(aobj_id))]
    
        #Read in and de-cosmify
        for line in range(len(obj_locs)):
        
            array, header = cosmics.fromfits(obj_locs[line])
            #array = array - bias #backwards for some reason
            c = cosmics.cosmicsimage(array, gain = 1.29, readnoise = 2.2, sigclip = 3.8, objlim = 3.0, sigfrac = 0.7, satlevel = -1)
            c.run(maxiter = 3) # can increase up to 4 to improve precision, but takes longer
            cosmics.tofits('Calibs/cosmicless/' + str(obj_id)+'_'+str(line+1)+'decos.fits', c.cleanarray, header)