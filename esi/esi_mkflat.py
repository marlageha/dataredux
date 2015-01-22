#! /usr/bin/python
##################################################
#  ESI_mkflat
#
#  USING ESI SUMMARY TABLE, COADD ALL DOME FLAT FRAMES
#  SUBTRACT BIAS AND NORMALIZE
#  WRITE TO FILE flat.fits
#  Kareem El-Badry 07/25/2014
##################################################

import pyfits
import numpy as np

def esi_mkflat(date):
    

    # READ BIAS FRAME
    bias = pyfits.getdata(str(date)+'/Calibs/bias_'+str(date)+'.fits')
    
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
        alldata.append([filename[line],dateobs[line],objname[line],imgtype[line],ra[line],dec[line],exptime[line],usable[line]])
            
    #Find junk files and good files:
    junk = []
    good = []
    for line in range(len(alldata)):
       if "no" in alldata[line][7]:
           junk.append(alldata[line])
       if "yes" in alldata[line][7]:
           good.append(alldata[line])
           
    #Find dome flats 
    dflat = []
    for line in range(len(good)):
        if "DmFlat" in good[line][3] and "Dome" in good[line][2]:
            dflat.append(good[line])
    print 'Number of dome flats to combine = ' + str(len(dflat))
    
    #Find pinhole flats
    pflat = []
    for line in range(len(good)):
        if "DmFlat" in good[line][3] and "Pinhole" in good[line][2]:
            pflat.append(good[line])
    print 'Number of pinhole flats to combine = ' + str(len(pflat)) + '\n'
            
    #Find path to domeflats
    domepath = []
    for line in range(len(dflat)):
        domepath.append(str(date)+"/Raw/"+str(dflat[line][0]))
    
    #Find path to pinflats
    pinpath = [str(date)+"/Raw/"+str(pflat[line][0]) for line in range(len(pflat))]
    
	# FOR EACH FILENAME, READ IMAGE AND ADD TO ARRAY
    print "flat cube shape | flat slice median   (if medians very different, consider getting rid of those flats)"
    alldflat=[]
    for i in domepath:    
        im = pyfits.getdata(i)
        im_med = (im - bias)/np.median(im - bias) #normalize 
        alldflat.append(im_med)
        print np.shape(alldflat),np.median(im)
    
    allpflat=[]
    for i in pinpath:    
        im = pyfits.getdata(i)
        im_med = (im - bias)/np.median(im - bias) #normalize 
        allpflat.append(im_med)
        print np.shape(allpflat),np.median(im)

    # MEDIAN COMBINE
    if len(alldflat) > 0:
        dome_flat = np.median(alldflat, axis=0)
        
        # WRITE TO DIRECTORY
        dfits = pyfits.PrimaryHDU(dome_flat[:, 25:2070])   # TRIM FLAT FIELD!!
        dfits.writeto(str(date)+'/Calibs/dome_flat_'+str(date)+'.fits', clobber=True)
        
    if len(allpflat) > 0:
        pin_flat = np.median(allpflat, axis=0)

        pfits = pyfits.PrimaryHDU(pin_flat[:, 25:2070])
        pfits.writeto(str(date)+'/Calibs/pinhole_flat_'+str(date)+'.fits', clobber=True)

    


