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

	# READ DATA SUMMARY
	wfile = pyfits.open('whirc_data.fits')
	whirc = wfile[1].data

	# FIND ALL BIAS FILES
	biastrue = (whirc.field('imgtype') == 'zero')
	biasfile = 'Raw/'+whirc[biastrue].field('filename')
	print 'Number of Bias to combine = ' + str(len(biasfile))


	# FOR EACH FILENAME, READ IMAGE AND ADD TO ARRAY
	allbias=[]
	for i in biasfile:    
	    im = pyfits.getdata(i)
	    im = im[0:2048,0:2048]
	    allbias.append(im)
	    print np.shape(allbias),np.median(im)

	# MEDIAN COMBINE
	bias = np.median(allbias, axis=0)

	# WRITE TO DIRECTORY
	fits = pyfits.PrimaryHDU(bias)
	fits.writeto('Calibs/bias.fits',clobber=True)


