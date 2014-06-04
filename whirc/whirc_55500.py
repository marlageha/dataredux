# Starts to reduce the science images for target 55500. 
# Run in directory above Raw, Calibs, Final

import pyfits
import numpy as np

def whirc_55500():
    
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
    
    #Find files related to 55500 in good files
    a55500 = []
    for line in range(len(good)):
        if "55500" in good[line][2]:
            a55500.append(good[line])
    print str(len(a55500)),"files for 55500"

    
    #Find files through J_filter
    a55500J = []
    for line in range(len(a55500)):
        if "J" in a55500[line][2]:
            a55500J.append(a55500[line])
    print str(len(a55500J)), "J files for 55500"
    
    #Find files through K_filter
    a55500K = []
    for line in range(len(a55500)):
        if "K" in a55500[line][2]:
            a55500K.append(a55500[line])
    print str(len(a55500K)), "K files for 55500"
    
    #Find J files from Night 1
    a55500JN1 = []
    for line in range(len(a55500J)):
        if "N1" in a55500J[line][0]:
            a55500JN1.append(a55500J[line])
    print str(len(a55500JN1)), "J files in N1"
    
    #Fink K files from Night 1
    a55500KN1 = []
    for line in range(len(a55500K)):
        if "N1" in a55500K[line][0]:
            a55500KN1.append(a55500K[line])
    print str(len(a55500KN1)), "K files in N1"
    
    #Find J files from Night 2
    a55500JN2 = []
    for line in range(len(a55500J)):
        if "N2" in a55500J[line][0]:
            a55500JN2.append(a55500J[line])
    print str(len(a55500JN2)), "J files in N2"
    
    #Find K files from Night 2
    a55500KN2 = []
    for line in range(len(a55500K)):
        if "N2" in a55500K[line][0]:
            a55500KN2.append(a55500K[line])
    print str(len(a55500KN2)), "K files in N2"
    
    #WRITE PATH TO J_FILES_N1
    J55500N1_locs = []
    for line in range(len(a55500JN1)):
        J55500N1_locs.append("Raw/"+str(a55500JN1[line][0]))
        
    #WRITE PATH TO K_FILES_N1
    K55500N1_locs = []
    for line in range(len(a55500KN1)):
        K55500N1_locs.append("Raw/"+str(a55500KN1[line][0]))
        
    #WRITE PATH OF J_FILES_N2
    J55500N2_locs = []
    for line in range(len(a55500JN2)):
        J55500N2_locs.append("Raw/"+str(a55500JN2[line][0]))
    
    #WRITE PATH OF K_FILES_N2
    K55500N2_locs = []
    for line in range(len(a55500KN2)):
        K55500N2_locs.append("Raw/"+str(a55500KN2[line][0]))
        
    #FOR EACH J_FILE_N1, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    all55500JN1 = []
    norm = np.median(pyfits.getdata(J55500N1_locs[0])[500:1500,500:1500]) #median of first flat
    for element in J55500N1_locs:
        im55500j = pyfits.getdata(element)
        im55500j = im55500j[0:2048,0:2048]
        im55500j = im55500j - bias - dark #because 60 sec exposure
        im55500j = im55500j*(norm/np.median(im55500j[500:1500,500:1500])) #normalize all flats to have same median as first flat 
        all55500JN1.append(im55500j)
        print np.shape(all55500JN1), np.median(im55500j)
        
    #FOR EACH K_FILE_N1, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    all55500KN1 = []
    norm = np.median(pyfits.getdata(K55500N1_locs[0])[500:1500,500:1500]) #median of first flat
    for element in K55500N1_locs:
        im55500k = pyfits.getdata(element)
        im55500k = im55500k[0:2048,0:2048]
        im55500k = im55500k - bias - dark #because 60 sec exposure
        im55500k = im55500k*(norm/np.median(im55500k[500:1500,500:1500])) #normalize all flats to have same median as first flat 
        all55500KN1.append(im55500k)
        print np.shape(all55500KN1), np.median(im55500k)
        
    #FOR EACH J_FILE_N2, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    all55500JN2 = []
    norm = np.median(pyfits.getdata(J55500N2_locs[0])[500:1500,500:1500]) #median of first flat
    for element in J55500N2_locs:
        im55500j = pyfits.getdata(element)
        im55500j = im55500j[0:2048,0:2048]
        im55500j = im55500j - bias - dark #because 60 sec exposure
        im55500j = im55500j*(norm/np.median(im55500j[500:1500,500:1500])) #normalize all flats to have same median as first flat 
        all55500JN2.append(im55500j)
        print np.shape(all55500JN2), np.median(im55500j)
        
    #FOR EACH K_FILE_N2, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    all55500KN2 = []
    norm = np.median(pyfits.getdata(K55500N2_locs[0])[500:1500,500:1500]) #median of first flat
    for element in K55500N2_locs:
        im55500k = pyfits.getdata(element)
        im55500k = im55500k[0:2048,0:2048]
        im55500k = im55500k - bias - dark #because 60 sec exposure
        im55500k = im55500k*(norm/np.median(im55500k[500:1500,500:1500])) #normalize all flats to have same median as first flat 
        all55500KN2.append(im55500k)
        print np.shape(all55500KN2), np.median(im55500k)
        
    #MEDIAN COMBINE
    J55500N1 = np.median(all55500JN1, axis = 0)
    K55500N1 = np.median(all55500KN1, axis = 0)
    J55500N2 = np.median(all55500JN2, axis = 0)
    K55500N2 = np.median(all55500KN2, axis = 0)
    
    #WRITE TO FILE
    j_sky1 = pyfits.PrimaryHDU(J55500N1)
    k_sky1 = pyfits.PrimaryHDU(K55500N1)
    j_sky2 = pyfits.PrimaryHDU(J55500N2)
    k_sky2 = pyfits.PrimaryHDU(K55500N2)
    
    j_sky1.writeto('Calibs/sky/sky_55500J_N1.fits', clobber = True)
    k_sky1.writeto('Calibs/sky/sky_55500K_N1.fits', clobber = True)
    j_sky2.writeto('Calibs/sky/sky_55500J_N2.fits', clobber = True)
    k_sky2.writeto('Calibs/sky/sky_55500K_N2.fits', clobber = True)
    
    
    