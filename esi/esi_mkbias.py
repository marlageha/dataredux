##################################################
#  ESI_mkbias
#
#  USING ESI SUMMARY TABLE, COADD ALL BIAS FRAMES
#  WRITE TO FILE bias.fits
#  Kareem El-Badry
#  07/25/2014
##################################################

import pyfits
import numpy as np

def esi_mkbias(date):
    '''    
    ##THIS PART FIXED TO WORK OFF UPDATED LOG; SOME 'IMGTYPE = ZERO' FILES ARE NOT BIAS##
	# READ DATA SUMMARY
	efile = pyfits.open('esi_data.fits')
	esi = efile[1].data

	# FIND ALL BIAS FILES
	biastrue = (esi.field('exptime') == 0)
	biasfile = 'Raw/'+esi[biastrue].field('filename')
	print 'Number of Bias to combine = ' + str(len(biasfile))
    '''

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
           
    #Find bias files
    bias = []
    for line in range(len(good)):
        if "Bias" in good[line][3]:
            bias.append(good[line])
    print 'Number of Bias to combine = ' + str(len(bias))
            
    #Find path to bias
    biaspath = []
    for line in range(len(bias)):
        biaspath.append(str(date)+"/Raw/"+str(bias[line][0]))
        
    
	# FOR EACH FILENAME, READ IMAGE AND ADD TO ARRAY
    allbias=[]
    for i in biaspath:    
	    im = pyfits.getdata(i)
	    allbias.append(im)
	    print np.shape(allbias),np.median(im)

	# MEDIAN COMBINE
    bias = np.median(allbias, axis=0)

    # WRITE TO DIRECTORY
    fits = pyfits.PrimaryHDU(bias)
    fits.writeto(str(date)+'/Calibs/bias_'+str(date)+'.fits',clobber=True)


