#! /usr/bin/python
##################################################
#  ESI_mkbias
#
#  USING ESI SUMMARY TABLE, COADD ALL BIAS FRAMES
#  WRITE TO FILE bias.fits
#  MG 4/14
##################################################

import pyfits
import numpy as np

def esi_mkbias():

	# READ DATA SUMMARY
	efile = pyfits.open('esi_data.fits')
	esi = efile[1].data

	# FIND ALL BIAS FILES
	biastrue = (esi.field('exptime') == 0)
	biasfile = 'Raw/'+esi[biastrue].field('filename')
	print 'Number of Bias to combine = ' + str(len(biasfile))


	# FOR EACH FILENAME, READ IMAGE AND ADD TO ARRAY
	allbias=[]
	for i in biasfile:    
	    im = pyfits.getdata(i)
	    allbias.append(im)
	    print np.shape(allbias),np.median(im)

	# MEDIAN COMBINE
	bias = np.median(allbias, axis=0)

	# WRITE TO DIRECTORY
	fits = pyfits.PrimaryHDU(bias)
	fits.writeto('Calibs/bias.fits',clobber=True)


