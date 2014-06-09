# ALLIGNS DITHERED IMAGES OF OBJ_ID SO WE CAN COADD THEM
# RUN IN DIRECTORY ABOVE RAW, CALIBS, FINAL

def whirc_align(obj_id):

    import pyfits
    import numpy as np
    import matplotlib.pyplot as plt
    import astropy.coordinates as coord
    import astropy.units as u
    from scipy.ndimage import interpolation as interp 
    pi = 3.1415926
    
    #READ IN DATA AND OBSERVING LOG
    #Read in the bias file to be subtracted from flats
    bias = pyfits.getdata("Calibs/bias.fits")
    bias = bias[0:2048,0:2048]
    
    #Read in the dark file to be subtracted from flats
    dark = pyfits.getdata("Calibs/dark.fits")
    dark = dark[0:2048,0:2048]
    
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
    print str(len(aobj_id)),"files for " + str(obj_id)

    raobj_id = []
    for line in range(len(aobj_id)):
        raobj_id.append(aobj_id[line][4])
        
    decobj_id = []
    for line in range(len(aobj_id)):
        decobj_id.append(aobj_id[line][5])
        
    rotobj_id = []
    for line in range(len(aobj_id)):
        rotobj_id.append(aobj_id[line][8])
    
    roffobj_id = []
    for line in range(len(aobj_id)):
        roffobj_id.append(aobj_id[line][9])
        
    doffobj_id = []
    for line in range(len(aobj_id)):
        doffobj_id.append(aobj_id[line][10])
    
    # Converts declination to a scalar (degrees)    
    decobj_iddegs = []
    for line in range(len(decobj_id)):
        dec = coord.Angle(decobj_id[line], unit=u.degree)
        decobj_iddegs.append(float(dec.degree))
    decobj_id = np.array(decobj_iddegs)
    #print decobj_id
    
    # Conversts RA to a scalar (hrs)
    raobj_iddegs = []
    for line in range(len(raobj_id)):
        ra = coord.Angle(raobj_id[line], unit=u.hour)
        raobj_iddegs.append(float(ra.degree))
    raobj_id = np.array(raobj_iddegs)
    #print raobj_id
    
    #Convert raoffset to a scalar(degrees)
    raoff = []
    for line in range(len(roffobj_id)):
        roff = coord.Angle(roffobj_id[line],unit=u.hour)
        raoff.append(float(roff.degree))
        if np.abs(raoff[line]) > 1:
            raoff[line] = -(360 - raoff[line])
    raoffobj_id = np.array(raoff)
    #print raoffobj_id
    
    #Convert decoffset to a scalar (degrees)
    decoff =[]
    for line in range(len(doffobj_id)):
        doff = coord.Angle(doffobj_id[line],unit=u.degree)
        decoff.append(float(doff.degree))
    decoffobj_id = np.array(decoff)
    #print decoffobj_id[0:10]
    
    # Converts rotangle to a scalar (deg)
    rot_degs = []
    for line in range(len(rotobj_id)):
        rot = coord.Angle(rotobj_id[line], unit=u.degree)
        rot_degs.append(float(rot.degree))
    rotobj_id = np.array(rot_degs)
    rot_ang = -1*rotobj_id*pi/180.0
    #print rot_ang
        
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
        Jobj_idN1_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'JN1_'+str(line + 1)+'.fits')
        
    #WRITE PATH TO K_FILES_N1
    Kobj_idN1_locs = []
    for line in range(len(aobj_idKN1)):
        Kobj_idN1_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'KN1_'+str(line + 1)+'.fits')
        
    #WRITE PATH OF J_FILES_N2
    Jobj_idN2_locs = []
    for line in range(len(aobj_idJN2)):
        Jobj_idN2_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'JN2_'+str(line + 1)+'.fits')
    
    #WRITE PATH OF K_FILES_N2
    Kobj_idN2_locs = []
    for line in range(len(aobj_idKN2)):
        Kobj_idN2_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'KN2_'+str(line + 1)+'.fits')
        
    #WRITE PATH OF J_FILES_N3
    Jobj_idN3_locs = []
    for line in range(len(aobj_idJN3)):
        Jobj_idN3_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'JN3_'+str(line + 1)+'.fits')
    
    #WRITE PATH OF K_FILES_N3
    Kobj_idN3_locs = []
    for line in range(len(aobj_idKN3)):
        Kobj_idN3_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'KN3_'+str(line + 1)+'.fits')
        
    for line in range(len(raoffobj_id)):
        if raoffobj_id[line] < -0.001:
            raoffobj_id[line] = float(295)
        elif raoffobj_id[line] > 0.001:
            raoffobj_id[line] = float(-295)
        elif np.abs(raoffobj_id[line]) < 0.001:
            raoffobj_id[line] = float(0.0)
        
    for line in range(len(decoffobj_id)):
        if decoffobj_id[line] < -0.001:
            decoffobj_id[line] = float(308)
        elif decoffobj_id[line] > 0.001:
            decoffobj_id[line] = float(-308)
        elif np.abs(decoffobj_id[line]) < 0.001:
            decoffobj_id[line] = float(0.0)
    
    xy_shift = []
    for line in range(len(raoffobj_id)):
        shift = [raoffobj_id[line],decoffobj_id[line]]
        xy_shift.append(shift)
        
    xy_shift_JN1 = []
    xy_shift_KN1 = []
    xy_shift_JN2 = []
    for line in range(len(xy_shift)):
        if 0 <= line < len(aobj_idJN1):
            xy_shift_JN1.append(xy_shift[line])
        elif len(aobj_idJN1) <= line <  len(aobj_idJN1)+len(aobj_idKN1):
            xy_shift_KN1.append(xy_shift[line])
        elif len(aobj_idJN1)+len(aobj_idKN1) <= line < len(aobj_idJN1)+len(aobj_idKN1) + len(aobj_idJN2):
            xy_shift_JN2.append(xy_shift[line])
        elif len(aobj_idJN1)+len(aobj_idKN1) + len(aobj_idJN2) <= line < len(aobj_idJN1)+len(aobj_idKN1) + len(aobj_idJN2) + len(aobj_idKN2):
            xy_shift_KN2.append(xy_shift[line])
        
    print xy_shift_JN1
    print xy_shift_KN1
    print xy_shift_JN2
    print xy_shift_KN2
        
    #Calculate necessary shift
    #ra_shift = raobj_id[4] - raobj_id
    #dec_shift = decobj_id[4] - decobj_id
    #print ra_shift*36511
    #print dec_shift*35892
    

    #WHIRC PIXEL SCALE IS 0.10030 X 0.0986 arcsec per pixel
    #translates to 35892 x 36511 pixel per degree
    
    #x_pix_scale = 35892
    #y_pix_scale = 36511
    
    #ra_off_pix = -1.0*raoffobj_id*36511
    #dec_off_pix = -1.0*decoffobj_id*35892
    
    #print ra_off_pix[0:6]
    #print dec_off_pix[0:6]
    
    #rot_ang = pi/2
    
    #x_offset = ra_off_pix*np.cos(rot_ang) - dec_off_pix*np.sin(rot_ang)
    #y_offset = ra_off_pix*np.sin(rot_ang) + dec_off_pix*np.cos(rot_ang)
    
    #READ IN SKY FILES
    #READ IN SCIENCE EXPOSURES
    
    #J_filter, Night 1
        
    print "reading in JN1..."
    if len(aobj_idJN1) > 0:
        shifteds = []
        for line in range(len(aobj_idJN1)):
            sci = pyfits.getdata(Jobj_idN1_locs[line])
            shifted = interp.shift(sci, xy_shift_JN1[line] ,order = 0)
            shifteds.append(shifted)
            file = pyfits.PrimaryHDU(shifted)
            file.writeto('Calibs/shifted/'+str(obj_id)+'_JN1_'+str(line+1)+'.fits',clobber = True)
    combined = sum(shifteds)
    median = np.median(shifteds, axis = 0)
    
    
    file = pyfits.PrimaryHDU(combined)
    mfile = pyfits.PrimaryHDU(median)
    file.writeto('Calibs/shifted/'+str(obj_id)+'_JN1_combined.fits',clobber=True)
    mfile.writeto('Calibs/shifted/'+str(obj_id)+'_JN1_median.fits',clobber=True)
    
    print "reading in KN1..."
    if len(aobj_idKN1) > 0:
        shifteds = []
        for line in range(len(aobj_idKN1)):
            sci = pyfits.getdata(Kobj_idN1_locs[line])
            shifted = interp.shift(sci, xy_shift_KN1[line] ,order = 0)
            shifteds.append(shifted)
            file = pyfits.PrimaryHDU(shifted)
            file.writeto('Calibs/shifted/'+str(obj_id)+'_KN1_'+str(line+1)+'.fits',clobber = True)
    combined = sum(shifteds)
    median = np.median(shifteds, axis = 0)
    
    
    file = pyfits.PrimaryHDU(combined)
    mfile = pyfits.PrimaryHDU(median)
    file.writeto('Calibs/shifted/'+str(obj_id)+'_KN1_combined.fits',clobber=True)
    mfile.writeto('Calibs/shifted/'+str(obj_id)+'_KN1_median.fits',clobber=True)
    