#! /usr/bin/python
##################################################
#  whirc_mkflat
#
#  USING whirc SUMMARY TABLE, COADD ALL DOME FLAT FRAMES
#  SUBTRACT BIAS AND NORMALIZE
#  WRITE TO FILE flat.fits
#
#  MG 4/14
##################################################

import pyfits
import numpy as np

def whirc_prepflat():

    # READ DATA SUMMARY

    efile = pyfits.open('whirc_data.fits')
    whirc = efile[1].data


    # READ BIAS FRAME
    bfile = pyfits.open('Calibs/bias.fits')
    bias = bfile[0].data

    # FIND ALL DOME FLAT FILES
    flattrue = (whirc.field('exptime') == 150)  # CAREFUL!!
    flatfile = 'Raw/'+whirc[flattrue].field('filename')
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


