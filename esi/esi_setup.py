#! /usr/bin/python
##################################################
#  ESI_setup
#
#  READ FITS HEADERS, CREATE SUMMARY FITS TABLE
#  >> mkdir Raw/
#  >> mkdir Calibs/ Final/
#  >
#  > import mgesi
#  > mgesi.esi_setup() 
#
#
#  MG 4/14
##################################################

import pyfits
from pyfits import Column
import numpy as np
import os

def esi_setup():

    # DEFINE HEADER INFO TO SAVE
    filename = []
    dateobs=[]
    objname=[]
    imgtype = []
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
            dateobs.append(hdr['DATE-OBS'])
            objname.append(hdr['OBJECT'])
            imgtype.append(hdr['OBSTYPE'])
            ra.append(hdr['RA'])
            dec.append(hdr['DEC'])
            exptime.append(hdr['EXPOSURE'])
            np.append(e2,hdr['EXPOSURE'])
            print hdr['OBJECT'], hdr['OBSTYPE'], hdr['EXPOSURE']

    # PARSE OBJECT NAMES
    exp = np.array([exptime])
    obj = np.array([objname])
    sciname = obj[exp > 600.0]    
    unique_obj = np.array([list(set(sciname))])

    print unique_obj


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
    c4 = Column(name='imgtype', format='A10',array=imgtype)
    c5 = Column(name='ra', format='A14',array=ra)
    c6 = Column(name='dec', format='A14',array=dec)
    c7 = Column(name='exptime',format='D',array=exptime)


    # WRITE TO DIRECTORY
    esi = pyfits.new_table([c1, c2, c3, c4, c5, c6, c7])
    esi.writeto('Logs/esi_data.fits',clobber=True)
    
    #NOW WRITE TO MORE EASILY EDITABLE DATA FILE
    hdulist = pyfits.open('Logs/esi_data.fits')
    table = hdulist[1].data
    
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
        
    print "Columns: filename|dateobs|objname|imgtype|ra|dec|exptime|yes/no"
    
    file = open('Logs/esi_info.dat','w')
    for line in range(len(filename)):
        file.write(" ".join([str(filename[line]), str(dateobs[line]), str(objname[line]), 
        str(imgtype[line]),str(ra[line]),str(dec[line]),str(exptime[line]),"yes",'\n']))
        
    file.close()
