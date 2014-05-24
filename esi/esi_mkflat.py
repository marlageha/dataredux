#! /usr/bin/python
##################################################
#  ESI_mkflat
#
#  USING ESI SUMMARY TABLE, COADD ALL DOME FLAT FRAMES
#  SUBTRACT BIAS AND NORMALIZE
#  WRITE TO FILE flat.fits
#  MG 4/14
##################################################

import pyfits
import numpy as np

def esi_mkflat():

    # READ DATA SUMMARY
    efile = pyfits.open('esi_data.fits')
    esi = efile[1].data

    # READ BIAS FRAME
    bfile = pyfits.open('Calibs/bias.fits')
    bias = bfile[0].data

    # FIND ALL DOME FLAT FILES
    flattrue = (esi.field('exptime') == 150)  # CAREFUL!!
    flatfile = 'Raw/'+esi[flattrue].field('filename')
    print len(flatfile)


    # FOR EACH FILENAME, READ IMAGE AND ADD TO ARRAY
    allflat=[]
    flatfile=flatfile[0:7]
    for i in flatfile:    
        im = pyfits.getdata(i)
        im = (im - bias)/np.median(im-bias)
        allflat.append(im)
        print i, np.shape(allflat),np.median(im)

    # MEDIAN COMBINE
    flat = np.median(allflat, axis=0)

    # WRITE TO DIRECTORY
    fits = pyfits.PrimaryHDU(flat[:,25:2070])   # TRIM FLAT FIELD!!
    fits.writeto('Calibs/flat.fits',clobber=True)


