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

# READ DATA SUMMARY
efile = pyfits.open('esi_data.fits')
esi = efile[1].data


# FIND ALL SCIENCE FILES
exp = np.array([esi.field('exptime')])
obj = np.array([esi.field('objname')])
sciname = obj[exp > 600.0]   
unique_obj = np.array([list(set(sciname))])


for targ in np.nditer(unique_obj):
	print targ
	scitrue = (esi.field('objname') == targ)
	scifile = 'Raw/'+esi[scitrue].field('filename')
#	sci=np.array(scifile)

	# FOR EACH FILENAME, READ IMAGE AND ADD TO ARRAY
	allsci=[]
	for i in scifile:
		print i
		im = pyfits.getdata(i)
		allsci.append(im)
		print np.shape(allsci)

	# MEDIAN COMBINE
	sci = np.median(allsci, axis=0)


	# WRITE TO DIRECTORY
	filename=str(targ)+'.fits'
	fits = pyfits.PrimaryHDU(sci)
	fits.writeto(filename,clobber=True)


