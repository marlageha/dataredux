#Run in same director as whirc_data.fits
# one directory above Raw/, Calibs/, Final/

import pyfits 
import numpy as np
import os

def whirc_info():
    
    hdulist = pyfits.open('whirc_data.fits')
    table = hdulist[1].data #because table is first extension
    
    #Define Columns
    filename = []
    dateobs = []
    objname = []
    imgtype = []
    ra = []
    dec = []
    exptime = []
    
    #Fill Columns from FITS table 
    for line in range(len(table)):
        filename.append(table[line][0])
        dateobs.append(table[line][1])
        objname.append(table[line][2])
        imgtype.append(table[line][3])
        ra.append(table[line][4])
        dec.append(table[line][5])
        exptime.append(table[line][6])
        
    print "Columns: filename|dateobs|objname|imgtype|ra|dec|exptime"
    
    file = open('whirc_info.dat','w')
    for line in range(len(filename)):
        file.write(" ".join([str(filename[line]), str(dateobs[line]), str(objname[line]), 
        str(imgtype[line]),str(ra[line]),str(dec[line]),str(exptime[line]),"yes",'\n']))
    file.close()