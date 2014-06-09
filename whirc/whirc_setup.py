#! /usr/bin/python
##################################################
#  whirc_setup
#
#  READ FITS HEADERS, CREATE SUMMARY FITS TABLE
#  >> mkdir Raw/
#  >> mkdir Calibs/ Final/
#  >
#  > import mgwhirc
#  > mgwhirc.whirc_setup() 
#
#
#  MG 4/14
##################################################

import pyfits
from pyfits import Column
import numpy as np
import os

def whirc_setup():

    # DEFINE HEADER INFO TO SAVE
    filename = []
    dateobs=[]
    objname=[]
    imgtype=[]
    ra=[]
    dec=[]
    exptime=[]
    rotangle = []
    ra_offset = []
    dec_offset = []
    e2=np.zeros(1)

    # READ ALL FILES IN /RAW
    for dr in os.listdir(os.getcwd()+'/Raw/'):
        dir = dr +'/'
        for file in os.listdir(os.getcwd()+'/Raw/'+dir):

            print dir+file

            # IF FILENAME ENDS WITH .gz, READ AND PARSE HEADERS
            if file.endswith('.gz'):
                hdr = pyfits.getheader('Raw/'+dir+file)
                filename.append(dir+file)
                dateobs.append(hdr['DATE-OBS'])
                objname.append(hdr['OBJECT'])
                imgtype.append(hdr['IMGTYPE'])
                ra.append(hdr['RA'])
                dec.append(hdr['DEC'])
                exptime.append(hdr['EXPTIME'])
                rotangle.append(hdr['ROTANGLE'])
                ra_offset.append(hdr['RAOFFST'])
                dec_offset.append(hdr['DECOFFST'])
                np.append(e2,hdr['EXPTIME'])
                print hdr['OBJECT'],hdr['IMGTYPE'],hdr['EXPTIME'],hdr['ROTANGLE'], hdr['RAOFFST'],hdr['DECOFFST']


    # PARSE OBJECT NAMES
    exp = np.array([exptime])
    obj = np.array([objname])
    sciname = obj[exp == 60.0]    
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
    c1 = Column(name='filename', format='30A', array=filename)
    c2 = Column(name='dateobs', format='A10', array=dateobs)
    c3 = Column(name='objname', format='A10',array=objname)
    c4 = Column(name='imgtype', format='A10',array=imgtype)
    c5 = Column(name='ra', format='A14',array=ra)
    c6 = Column(name='dec', format='A14',array=dec)
    c7 = Column(name='exptime',format='D',array=exptime)
    c8 = Column(name='rotangle',format='A14',array=rotangle)
    c9 = Column(name='raoffset',format='A14',array=ra_offset)
    c10= Column(name='decoffset',format='A14',array=dec_offset)



    # WRITE TO DIRECTORY
    whirc = pyfits.new_table([c1, c2, c3, c4, c5,c6,c7,c8,c9,c10])
    whirc.writeto('whirc_data.fits',clobber=True)
        
