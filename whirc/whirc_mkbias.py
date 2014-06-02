#! /usr/bin/python
##################################################
#  WHIRC_mkbias
#
#  USING WHIRC SUMMARY TABLE, COADD ALL BIAS FRAMES
#  WRITE TO FILE bias.fits
#  MG 4/14
##################################################

import pyfits
import numpy as np

def whirc_mkbias():

##THIS PART FIXED TO WORK OFF UPDATED LOG; SOME 'IMGTYPE = ZERO' FILES ARE NOT BIAS##
	# READ DATA SUMMARY
	#wfile = pyfits.open('whirc_data.fits')
	#whirc = wfile[1].data

	# FIND ALL BIAS FILES
	#biastrue = (whirc.field('imgtype') == 'zero')
	#biasfile = 'Raw/'+whirc[biastrue].field('filename')
	#print 'Number of Bias to combine = ' + str(len(biasfile))
    
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
           
    #Find bias files
    bias = []
    for line in range(len(good)):
        if "bias" in good[line][2]:
            bias.append(good[line])
    print 'Number of Bias to combine = ' + str(len(bias))
    
    #Find dark files
    darks = []
    for line in range(len(good)):
        if "dark" in good[line][2]:
            darks.append(good[line])
    print 'Number of darks to combine = ' + str(len(darks))
            
    #Find path to bias
    biaspath = []
    for line in range(len(bias)):
        biaspath.append("Raw/"+str(bias[line][0]))
   
    #Find path to darks
    darkpath = []
    for line in range(len(darks)):
        darkpath.append("Raw/"+str(darks[line][0]))

    # FOR EACH BIAS, READ IMAGE AND ADD TO ARRAY
    allbias=[]
    for i in biaspath:    
        im = pyfits.getdata(i)
        im = im[0:2048,0:2048]
        allbias.append(im)
        print np.shape(allbias),np.median(im)

	# MEDIAN COMBINE
    bias = np.median(allbias, axis=0)

	# WRITE TO DIRECTORY
    fits = pyfits.PrimaryHDU(bias)
    fits.writeto('Calibs/bias.fits',clobber=True)
    
    print "writing darks..."
    # FOR EACH DARK, READ IMAGE AND ADD TO ARRAY
    alldarks = []
    for i in darkpath:
        im = pyfits.getdata(i)
        im = im[0:2048,0:2048]
        alldarks.append(im)
        print np.shape(alldarks), np.median(im)
        
    #WRITE TO DIRECTORY


