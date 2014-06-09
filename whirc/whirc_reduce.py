# Applies flats, darks, bias, and sky subtraction
# To all exposures of obj_id

def whirc_reduce(obj_id):
    
    import pyfits
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.coordinates as coord
    import astropy.units as u
    
    #READ IN DATA AND OBSERVING LOG
    #Read in the bias file to be subtracted from flats
    bias = pyfits.getdata("Calibs/bias.fits")
    bias = bias[0:2048,0:2048]
    
    #Read in the dark file to be subtracted from flats
    dark = pyfits.getdata("Calibs/dark.fits")
    dark = dark[0:2048,0:2048]
    
    #Read in flats
    J_flat = pyfits.getdata('Calibs/New_J_Flat.fits')
    K_flat = pyfits.getdata('Calibs/New_K_Flat.fits')
    
    #Read in data from "whirc_info.dat" (the observation log)
    im1 = open('whirc_info.dat','r')
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
    rotangle = []
    raoffset = []
    decoffset = []
    
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
        rotangle.append(p[8])
        raoffset.append(p[9])
        decoffset.append(p[10])
        
    
    #Rewrite in a more convenient array with format array[line][element]    
    alldata = []
    
    for line in range(len(usable)):
        alldata.append([filename[line],dateobs[line],objname[line],imgtype[line],ra[line],dec[line],exptime[line],usable[line],rotangle[line],raoffset[line],decoffset[line]])
    
    #Find junk files and good files:
    junk = []
    good = []
    for line in range(len(alldata)):
       if "no" in alldata[line][7]:
           junk.append(alldata[line])
       if "yes" in alldata[line][7]:
           good.append(alldata[line])
    
    #Find files related to obj_id in good files
    aobj_id = []
    for line in range(len(good)):
        if str(obj_id) in good[line][2]:
            aobj_id.append(good[line])
    print str(len(aobj_id)),"files for "+str(obj_id)
    
    #Find files through J_filter
    aobj_idJ = []
    for line in range(len(aobj_id)):
        if "J" in aobj_id[line][2]:
            aobj_idJ.append(aobj_id[line])
    print str(len(aobj_idJ)), "J files for "+str(obj_id)
    
    #Find files through K_filter
    aobj_idK = []
    for line in range(len(aobj_id)):
        if "K" in aobj_id[line][2]:
            aobj_idK.append(aobj_id[line])
    print str(len(aobj_idK)), "K files for "+str(obj_id)
    
    
    #Find J files from Night 1
    aobj_idJN1 = []
    for line in range(len(aobj_idJ)):
        if "N1" in aobj_idJ[line][0]:
            aobj_idJN1.append(aobj_idJ[line])
    print str(len(aobj_idJN1)), "J files in N1"
    
    #Fink K files from Night 1
    aobj_idKN1 = []
    for line in range(len(aobj_idK)):
        if "N1" in aobj_idK[line][0]:
            aobj_idKN1.append(aobj_idK[line])
    print str(len(aobj_idKN1)), "K files in N1"
    
    #Find J files from Night 2
    aobj_idJN2 = []
    for line in range(len(aobj_idJ)):
        if "N2" in aobj_idJ[line][0]:
            aobj_idJN2.append(aobj_idJ[line])
    print str(len(aobj_idJN2)), "J files in N2"
    
    #Find K files from Night 2
    aobj_idKN2 = []
    for line in range(len(aobj_idK)):
        if "N2" in aobj_idK[line][0]:
            aobj_idKN2.append(aobj_idK[line])
    print str(len(aobj_idKN2)), "K files in N2"
    
    #Find J files from Night 3
    aobj_idJN3 = []
    for line in range(len(aobj_idJ)):
        if "N3" in aobj_idJ[line][0]:
            aobj_idJN3.append(aobj_idJ[line])
    print str(len(aobj_idJN3)), "J files in N3"
    
    #Find K files from Night 3
    aobj_idKN3 = []
    for line in range(len(aobj_idK)):
        if "N3" in aobj_idK[line][0]:
            aobj_idKN3.append(aobj_idK[line])
    print str(len(aobj_idKN3)), "K files in N3"
    
    #WRITE PATH TO J_FILES_N1
    Jobj_idN1_locs = []
    for line in range(len(aobj_idJN1)):
        Jobj_idN1_locs.append("Raw/"+str(aobj_idJN1[line][0]))
        
    #WRITE PATH TO K_FILES_N1
    Kobj_idN1_locs = []
    for line in range(len(aobj_idKN1)):
        Kobj_idN1_locs.append("Raw/"+str(aobj_idKN1[line][0]))
        
    #WRITE PATH OF J_FILES_N2
    Jobj_idN2_locs = []
    for line in range(len(aobj_idJN2)):
        Jobj_idN2_locs.append("Raw/"+str(aobj_idJN2[line][0]))
    
    #WRITE PATH OF K_FILES_N2
    Kobj_idN2_locs = []
    for line in range(len(aobj_idKN2)):
        Kobj_idN2_locs.append("Raw/"+str(aobj_idKN2[line][0]))
        
    #WRITE PATH OF J_FILES_N3
    Jobj_idN3_locs = []
    for line in range(len(aobj_idJN3)):
        Jobj_idN3_locs.append("Raw/"+str(aobj_idJN3[line][0]))
    
    #WRITE PATH OF K_FILES_N3
    Kobj_idN3_locs = []
    for line in range(len(aobj_idKN3)):
        Kobj_idN3_locs.append("Raw/"+str(aobj_idKN3[line][0]))
        
    #READ IN SKY FILES
    #READ IN SCIENCE EXPOSURES
    
    #J_filter, Night 1
    print "reading in JN1..."
    if len(aobj_idJN1) > 0:
        skyJN1 = pyfits.getdata('Calibs/sky/sky_'+str(obj_id)+'J_N1.fits')
        
        for line in range(len(aobj_idJN1)):
            sci = pyfits.getdata(Jobj_idN1_locs[line])[0:2048,0:2048]
            reduced_sci = (sci - bias - dark - skyJN1)/J_flat
            file = pyfits.PrimaryHDU(reduced_sci)
            file.writeto('Calibs/reduced/'+str(obj_id)+'_'+'JN1_'+str(line + 1)+'.fits', clobber = True)
            
    #K_filter, Night 1
    print "reading in KN1..."    
    if len(aobj_idKN1) > 0:
        skyKN1 = pyfits.getdata('Calibs/sky/sky_'+str(obj_id)+'K_N1.fits')
        
        for line in range(len(aobj_idKN1)):
            sci = pyfits.getdata(Kobj_idN1_locs[line])[0:2048,0:2048]
            reduced_sci = (sci - bias - dark - skyKN1)/K_flat
            file = pyfits.PrimaryHDU(reduced_sci)
            file.writeto('Calibs/reduced/'+str(obj_id)+'_'+'KN1_'+str(line + 1)+'.fits', clobber = True)
        
    #J_filter, Night 2
    print "reading in JN2..."    
    if len(aobj_idJN2) > 0:
        skyJN2 = pyfits.getdata('Calibs/sky/sky_'+str(obj_id)+'J_N2.fits')
        
        for line in range(len(aobj_idJN2)):
            sci = pyfits.getdata(Jobj_idN2_locs[line])[0:2048,0:2048]
            reduced_sci = (sci - bias - dark - skyJN2)/J_flat
            file = pyfits.PrimaryHDU(reduced_sci)
            file.writeto('Calibs/reduced/'+str(obj_id)+'_'+'JN2_'+str(line + 1)+'.fits', clobber = True)

    #K_filter, Night 2
    print "reading in KN2..."
    if len(aobj_idKN2) > 0:
        skyKN2 = pyfits.getdata('Calibs/sky/sky_'+str(obj_id)+'K_N2.fits')
        
        for line in range(len(aobj_idKN2)):
            sci = pyfits.getdata(Kobj_idN2_locs[line])[0:2048,0:2048]
            reduced_sci = (sci - bias - dark - skyKN2)/K_flat
            file = pyfits.PrimaryHDU(reduced_sci)
            file.writeto('Calibs/reduced/'+str(obj_id)+'_'+'KN2_'+str(line + 1)+'.fits', clobber = True)
        
    #J_filter, Night 3
    print "reading in JN3..."  
    if len(aobj_idJN3) > 0:
        skyJN3 = pyfits.getdata('Calibs/sky/sky_'+str(obj_id)+'J_N3.fits')
        
        for line in range(len(aobj_idJN3)):
            sci = pyfits.getdata(Jobj_idN3_locs[line])[0:2048,0:2048]
            reduced_sci = (sci - bias - dark - skyJN3)/J_flat
            file = pyfits.PrimaryHDU(reduced_sci)
            file.writeto('Calibs/reduced/'+str(obj_id)+'_'+'JN3_'+str(line + 1)+'.fits', clobber = True)
        
    if len(aobj_idKN3) > 0:
        skyKN3 = pyfits.getdata('Calibs/sky/sky_'+str(obj_id)+'K_N3.fits')
        
        for line in range(len(aobj_idKN3)):
            sci = pyfits.getdata(Kobj_idN3_locs[line])[0:2048,0:2048]
            reduced_sci = (sci - bias - dark - skyKN3)/K_flat
            file = pyfits.PrimaryHDU(reduced_sci)
            file.writeto('Calibs/reduced/'+str(obj_id)+'_KN3_'+str(line + 1)+'.fits', clobber = True)