#! /usr/bin/python
##################################################
#  ESI_setup
#
#  USING ESI SUMMARY TABLE, COADD ALL BIAS FRAMES
#  WRITE TO FILE bias.fits
#  MG 4/14
##################################################

import pyfits
from pyfits import Column
import numpy as np
import os

# DEFINE HEADER INFO TO SAVE
filename = []
dateobs=[]
objname=[]
ra=[]
dec=[]
exptime=[]
e2=np.zeros(1)

# READ ALL FILES IN /RAW
for file in os.listdir(os.getcwd()+'/Raw'):

    # IF FILENAME ENDS WITH .gz, READ AND PARSE HEADERS
    if file.endswith('.gz'):
        hdr = pyfits.getheader('Raw/'+file)
        filename.append(file)
        print hdr['RA'],hdr['TARGNAME'],hdr['EXPOSURE'],hdr['OBJECT']
        dateobs.append(hdr['DATE-OBS'])
        objname.append(hdr['OBJECT'])
        ra.append(hdr['RA'])
        dec.append(hdr['DEC'])
        exptime.append(hdr['EXPOSURE'])
        np.append(e2,hdr['EXPOSURE'])

# PARSE OBJECT NAMES
exp = np.array([exptime])
obj = np.array([objname])
sciname = obj[exp > 600.0]    
unique_obj = np.array([list(set(sciname))])


# ENUMERATE EACH SET OF SCIECE OBJECTS
#objnum = np.zeros(len(obj))
#num=0
#for target in np.nditer(unique_obj):
#    num  += 1
#    print target,num
#    q = np.where(obj == target)
#    objname[q] = num
    


# CREATE DATA TABLE
c1 = Column(name='filename', format='20A', array=filename)
c2 = Column(name='dateobs', format='A10', array=dateobs)
c3 = Column(name='objname', format='A10',array=objname)
c4 = Column(name='ra', format='A14',array=ra)
c5 = Column(name='dec', format='A14',array=dec)
c6 = Column(name='exptime',format='D',array=exptime)





# WRITE TO DIRECTORY
esi = pyfits.new_table([c1, c2, c3, c4, c5,c6])
esi.writeto('esi_data.fits',clobber=True)
