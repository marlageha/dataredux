# Sorts data from Raw/ folder (in N1/, N2/, N3) into
# several different folders for flats, bias, different objects, etc

# Run in the directory above Raw/, Calibs/, and Final/ 

import pyfits
import numpy as np
import matplotlib.pyplot as plt

def whirc_sort():
    
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

    #FOR EACH J-FLAT-ON, READ IMAGE AND ADD TO ARRAY
    J_flat_on_median = []
    all_J_flat_on = [] 
    for element in J_flat_on_locs:
        im_J_on = pyfits.getdata(element)
        im_J_on = im_J_on[0:2048,0:2048]
        all_J_flat_on.append(im_J_on)
        J_flat_on_median.append(float(np.median(im_J_on)))
        print np.shape(all_J_flat_on), np.median(im_J_on)
        
    #Median combine
    JFO = np.median(all_J_flat_on, axis = 0)
    
    #Write to directory
    JFFO = pyfits.PrimaryHDU(JFO)
    JFFO.writeto('Calibs/J_flat_on.fits', clobber = True)
    
    plt.hist(J_flat_on_median)
        