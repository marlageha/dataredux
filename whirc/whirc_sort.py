# Right now, this makes flats from raw data
# Bias and dark are subracted
# Run in the directory above Raw/, Calibs/, and Final/ 

import pyfits
import numpy as np
import matplotlib.pyplot as plt

def whirc_sort():
    
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
           
    #Find J-Flats
    J_flats = []
    for line in range(len(good)):
        if "J-flat" in good[line][2]:
            J_flats.append(good[line])
    
    #Find K-Flats
    K_flats = []
    for line in range(len(good)):
        if "K-flat" in good[line][2]:
            K_flats.append(good[line])
            
    #Find on & off flats in J_flats
    J_flat_on = []
    J_flat_off = []
    for line in range(len(J_flats)):
       if "off" in J_flats[line][2]:
           J_flat_off.append(J_flats[line])
       else:
           J_flat_on.append(J_flats[line])
    
    #Find on & off flats in K_flats
    K_flat_on = []
    K_flat_off = []
    for line in range(len(K_flats)):
        if "off" in K_flats[line][2]:
            K_flat_off.append(K_flats[line])
        else:
            K_flat_on.append(K_flats[line])
   
    #Gives path to J-flats-on
    J_flat_on_locs = []
    for line in range(len(J_flat_on)):
        J_flat_on_locs.append("Raw/"+str(J_flat_on[line][0]))
    
    #Gives path to J-flats-off
    J_flat_off_locs = []
    for line in range(len(J_flat_off)):
        J_flat_off_locs.append("Raw/"+str(J_flat_off[line][0]))
    
    #Gives path to K-flats-on
    K_flat_on_locs = []
    for line in range(len(K_flat_on)):
        K_flat_on_locs.append("Raw/"+str(K_flat_on[line][0]))
    
    #Gives path to K-flats-off
    K_flat_off_locs = []
    for line in range(len(K_flat_off)):
        K_flat_off_locs.append("Raw/"+str(K_flat_off[line][0]))

    #FOR EACH J-FLAT-ON, READ IMAGE, SUBTRACT BIAS, AND ADD TO ARRAY
    all_J_flat_on = [] 
    norm = np.median(pyfits.getdata(J_flat_on_locs[0])[500:1500,500:1500]) #median of first flat
    for element in J_flat_on_locs:
        im_J_on = pyfits.getdata(element)
        im_J_on = im_J_on[0:2048,0:2048]
        im_J_on = im_J_on - bias - (5.0/60.0)*dark #because flats are only 5 seconds and darks are 60 seconds
        im_J_on = im_J_on*(norm/np.median(im_J_on[500:1500,500:1500])) #normalize all flats to have same median as first flat 
        all_J_flat_on.append(im_J_on)
        print np.shape(all_J_flat_on), np.median(im_J_on)
        
    #Median combine
    JFO = np.median(all_J_flat_on, axis = 0)
    
    #Write to directory
    JFFO = pyfits.PrimaryHDU(JFO)
    JFFO.writeto('Calibs/J_flat_on.fits', clobber = True)
    print "Done with J_flat_on"
        
    #FOR EACH J-FLAT-OFF, READ IMAGE, SUBTRACT BIAS, AND ADD TO ARRAY
    all_J_flat_off = [] 
    norm = np.median(pyfits.getdata(J_flat_off_locs[1])) #median of second flat (first one is weird)
    for element in J_flat_off_locs:
        im_J_off = pyfits.getdata(element)
        im_J_off = im_J_off[0:2048,0:2048]
        im_J_off = im_J_off - bias - (5.0/60.0)*dark #because flats are only 5 seconds and darks are 60 seconds
        im_J_off = im_J_off*(norm/np.median(im_J_off)) #normalize all flats to have same median as first flat 
        all_J_flat_off.append(im_J_off)
        print np.shape(all_J_flat_off), np.median(im_J_off)
    
    #Median combine
    JFOF = np.median(all_J_flat_off, axis = 0)
    
    #Write to directory
    JFFOF = pyfits.PrimaryHDU(JFOF)
    JFFOF.writeto('Calibs/J_flat_off.fits', clobber = True)
    print "Done with J_flat_off"
    
    #FOR EACH K-FLAT-ON, READ IMAGE, SUBTRACT BIAS, AND ADD TO ARRAY
    all_K_flat_on = [] 
    norm = np.median(pyfits.getdata(K_flat_on_locs[1])[500:1500,500:1500]) #median of second flat
    for element in K_flat_on_locs:
        im_K_on = pyfits.getdata(element)
        im_K_on = im_K_on[0:2048,0:2048]
        im_K_on = im_K_on - bias - (5.0/60.0)*dark #because flats are only 5 seconds and darks are 60 seconds
        im_K_on = im_K_on*(norm/np.median(im_K_on[500:1500,500:1500])) #normalize all flats to have same median as first flat 
        all_K_flat_on.append(im_K_on)
        print np.shape(all_K_flat_on), np.median(im_K_on)
        
    #Median combine
    KFO = np.median(all_K_flat_on, axis = 0)
    
    #Write to directory
    KFFO = pyfits.PrimaryHDU(KFO)
    KFFO.writeto('Calibs/K_flat_on.fits', clobber = True)
    print "Done with K_flat_on"
    
    #FOR EACH K-FLAT-OFF, READ IMAGE, SUBTRACT BIAS, AND ADD TO ARRAY
    all_K_flat_off = [] 
    norm = np.median(pyfits.getdata(K_flat_off_locs[1])) #median of second flat
    for element in K_flat_off_locs:
        im_K_off = pyfits.getdata(element)
        im_K_off = im_K_off[0:2048,0:2048]
        im_K_off = im_K_off - bias - (5.0/60.0)*dark #because flats are only 5 seconds and darks are 60 seconds
        im_K_off = im_K_off*(norm/np.median(im_K_off)) #normalize all flats to have same median as first flat 
        all_K_flat_off.append(im_K_off)
        print np.shape(all_K_flat_off), np.median(im_K_off)
    
    #Median combine
    KFOF = np.median(all_K_flat_off, axis = 0)
    
    #Write to directory
    KFFOF = pyfits.PrimaryHDU(KFOF)
    KFFOF.writeto('Calibs/K_flat_off.fits', clobber = True)
    print "Done with K_flat_off"
    
    #SUBTRACTING DARKS FROM LIGHTS
    J_master = JFO - JFOF
    K_master = KFO - KFOF
    
    J_mas = pyfits.PrimaryHDU(J_master)
    J_mas.writeto('Calibs/J_master_flat.fits', clobber = True)
    
    K_mas = pyfits.PrimaryHDU(K_master)
    K_mas.writeto('Calibs/K_master_flat.fits', clobber = True)
        