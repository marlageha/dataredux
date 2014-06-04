# Starts to reduce the science images for target 121130. 
# Run in directory above Raw, Calibs, Final

import pyfits
import numpy as np

def whirc_121130():
    
    #Read in the bias file to be subtracted from flats
    bias = pyfits.getdata("Calibs/bias.fits")
    bias = bias[0:2048,0:2048]
    
    #Read in the dark file to be subtracted from flats
    dark = pyfits.getdata("Calibs/dark.fits")
    dark = dark[0:2048,0:2048]
    
    #Read in data from "whirc_info.dat" (the observation log)
    im1 = open('whirc_info.dat','r')
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
        alldata.append([filename[line],dateobs[line],objname[line],imgtype[line],ra[line],dec[line],exptime[line],usable[line]])
    
    #Find junk files and good files:
    junk = []
    good = []
    for line in range(len(alldata)):
       if "no" in alldata[line][7]:
           junk.append(alldata[line])
       if "yes" in alldata[line][7]:
           good.append(alldata[line])
           
    #Find files related to 121130 in good files
    a121130 = []
    for line in range(len(good)):
        if "121130" in good[line][2]:
            a121130.append(good[line])
    print str(len(a121130)),"files for 121130"
    
    #Find files through J_filter
    a121130J = []
    for line in range(len(a121130)):
        if "J" in a121130[line][2]:
            a121130J.append(a121130[line])
    print str(len(a121130J)), "J files for 121130"
    
    #Find files through K_filter
    a121130K = []
    for line in range(len(a121130)):
        if "K" in a121130[line][2]:
            a121130K.append(a121130[line])
    print str(len(a121130K)), "K files for 121130"
    
    #WRITE PATH TO J_FILES
    J121130_locs = []
    for line in range(len(a121130J)):
        J121130_locs.append("Raw/"+str(a121130J[line][0]))
        
    #WRITE PATH TO K_FILES
    K121130_locs = []
    for line in range(len(a121130K)):
        K121130_locs.append("Raw/"+str(a121130K[line][0]))
        
    #FOR EACH J_FILE, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    all121130J = []
    norm = np.median(pyfits.getdata(J121130_locs[0])[500:1500,500:1500]) #median of first flat
    for element in J121130_locs:
        im121130j = pyfits.getdata(element)
        im121130j = im121130j[0:2048,0:2048]
        im121130j = im121130j - bias - dark #because 60 sec exposure
        im121130j = im121130j*(norm/np.median(im121130j[500:1500,500:1500])) #normalize all flats to have same median as first flat 
        all121130J.append(im121130j)
        print np.shape(all121130J), np.median(im121130j)
        
    #FOR EACH K_FILE, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    all121130K = []
    norm = np.median(pyfits.getdata(K121130_locs[0])[500:1500,500:1500]) #median of first flat
    for element in K121130_locs:
        im121130k = pyfits.getdata(element)
        im121130k = im121130k[0:2048,0:2048]
        im121130k = im121130k - bias - dark #because 60 sec exposure
        im121130k = im121130k*(norm/np.median(im121130k[500:1500,500:1500])) #normalize all flats to have same median as first flat 
        all121130K.append(im121130k)
        print np.shape(all121130K), np.median(im121130k)
        
    #MEDIAN COMBINE
    print "median combining..."
    J121130 = np.median(all121130J, axis = 0)
    K121130 = np.median(all121130K, axis = 0)
    
    #WRITE TO FILE
    j_sky = pyfits.PrimaryHDU(J121130)
    k_sky = pyfits.PrimaryHDU(K121130)
    
    j_sky.writeto('Calibs/sky/sky_121130J.fits', clobber = True)
    k_sky.writeto('Calibs/sky/sky_121130K.fits', clobber = True)
    