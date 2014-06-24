# Makes Sky Images 
# Takes an argument of object name;
# e.g > whirc_mksky(55000) will make sky images 
# for obj_id. 

def whirc_mksky(obj_id):
    
    import pyfits
    import numpy as np
    
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
    
    #Rewrite in a more convenient array with format array[line][element]    
    alldata = []
    
    for line in range(len(usable)):
        alldata.append([filename[line],dateobs[line],objname[line],imgtype[line],ra[line],dec[line],exptime[line],usable[line]])
    
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
        
    #FOR EACH J_FILE_N1, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    print "reading in J_N1..."
    allobj_idJN1 = []
    if len(Jobj_idN1_locs) > 0:
        #norm = np.median(pyfits.getdata(Jobj_idN1_locs[0])[500:1500,500:1500]) #median of first flat
        
        for element in Jobj_idN1_locs:
            imobj_idj = pyfits.getdata(element)
            imobj_idj = imobj_idj[0:2048,0:2048]
            imobj_idj = imobj_idj - bias - dark #because 60 sec exposure
            #imobj_idj = imobj_idj*(norm/np.median(imobj_idj[500:1500,500:1500])) #normalize all flats to have same median as first flat 
            allobj_idJN1.append(imobj_idj)
            print np.shape(allobj_idJN1), np.median(imobj_idj)
    
        
        #MEDIAN COMBINE
        #WRITE TO FILE
    
        print "median combining J_N1..."
        Jobj_idN1 = np.median(allobj_idJN1, axis = 0)
        j_sky1 = pyfits.PrimaryHDU(Jobj_idN1)
        j_sky1.writeto('Calibs/sky/sky_'+str(obj_id)+'J_N1.fits', clobber = True)
    
        
    #FOR EACH K_FILE_N1, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    print "reading in K_N1..."
    allobj_idKN1 = []
    if len(Kobj_idN1_locs) > 0:
        #norm = np.median(pyfits.getdata(Kobj_idN1_locs[0])[500:1500,500:1500]) #median of first flat
        
        for element in Kobj_idN1_locs:
            imobj_idk = pyfits.getdata(element)
            imobj_idk = imobj_idk[0:2048,0:2048]
            imobj_idk = imobj_idk - bias - dark #because 60 sec exposure
            #imobj_idk = imobj_idk*(norm/np.median(imobj_idk[500:1500,500:1500])) #normalize all flats to have same median as first flat 
            allobj_idKN1.append(imobj_idk)
            print np.shape(allobj_idKN1), np.median(imobj_idk)
            
        print "median combining K_N1..."
        Kobj_idN1 = np.median(allobj_idKN1, axis = 0)
        k_sky1 = pyfits.PrimaryHDU(Kobj_idN1)
        k_sky1.writeto('Calibs/sky/sky_'+str(obj_id)+'K_N1.fits', clobber = True)
            

        
    #FOR EACH J_FILE_N2, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    print "reading in J_N2..."
    allobj_idJN2 = []
    if len(Jobj_idN2_locs) > 0:
        #norm = np.median(pyfits.getdata(Jobj_idN2_locs[0])[500:1500,500:1500]) #median of first flat
        
        for element in Jobj_idN2_locs:
            imobj_idj = pyfits.getdata(element)
            imobj_idj = imobj_idj[0:2048,0:2048]
            imobj_idj = imobj_idj - bias - dark #because 60 sec exposure
            #imobj_idj = imobj_idj*(norm/np.median(imobj_idj[500:1500,500:1500])) #normalize all flats to have same median as first flat 
            allobj_idJN2.append(imobj_idj)
            print np.shape(allobj_idJN2), np.median(imobj_idj)
            
        print "median combining J_N2..."
        Jobj_idN2 = np.median(allobj_idJN2, axis = 0)
        j_sky2 = pyfits.PrimaryHDU(Jobj_idN2)
        j_sky2.writeto('Calibs/sky/sky_'+str(obj_id)+'J_N2.fits', clobber = True)
        
    #FOR EACH K_FILE_N2, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    print "reading in K_N2..."
    allobj_idKN2 = []
    if len(Kobj_idN2_locs) > 0:
        #norm = np.median(pyfits.getdata(Kobj_idN2_locs[0])[500:1500,500:1500]) #median of first flat
        
        for element in Kobj_idN2_locs:
            imobj_idk = pyfits.getdata(element)
            imobj_idk = imobj_idk[0:2048,0:2048]
            imobj_idk = imobj_idk - bias - dark #because 60 sec exposure
            #imobj_idk = imobj_idk*(norm/np.median(imobj_idk[500:1500,500:1500])) #normalize all flats to have same median as first flat 
            allobj_idKN2.append(imobj_idk)
            print np.shape(allobj_idKN2), np.median(imobj_idk)
            
        print "median combining K_N2..."
        Kobj_idN2 = np.median(allobj_idKN2, axis = 0)
        k_sky2 = pyfits.PrimaryHDU(Kobj_idN2)
        k_sky2.writeto('Calibs/sky/sky_'+str(obj_id)+'K_N2.fits', clobber = True)
        
    #FOR EACH J_FILE_N3, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    print "reading in J_N3..."
    allobj_idJN3 = []
    if len(Jobj_idN3_locs) > 0:
        #norm = np.median(pyfits.getdata(Jobj_idN3_locs[0])[500:1500,500:1500]) #median of first flat
        
        for element in Jobj_idN3_locs:
            imobj_idj = pyfits.getdata(element)
            imobj_idj = imobj_idj[0:2048,0:2048]
            imobj_idj = imobj_idj - bias - dark #because 60 sec exposure
            #imobj_idj = imobj_idj*(norm/np.median(imobj_idj[500:1500,500:1500])) #normalize all flats to have same median as first flat 
            allobj_idJN3.append(imobj_idj)
            print np.shape(allobj_idJN3), np.median(imobj_idj)
        
        print "median combining J_N3..."
        Jobj_idN3 = np.median(allobj_idJN3, axis = 0)
        j_sky3 = pyfits.PrimaryHDU(Jobj_idN3)
        j_sky3.writeto('Calibs/sky/sky_'+str(obj_id)+'J_N3.fits', clobber = True)
        
        
    #FOR EACH K_FILE_N3, READ IMAGE, SUBTRACT BIAS AND DARKS, AND ADD TO ARRAY
    print "reading in K_N3..."
    allobj_idKN3 = []
    if len(Kobj_idN3_locs) > 0:
        #norm = np.median(pyfits.getdata(Kobj_idN3_locs[0])[500:1500,500:1500]) #median of first flat


        for element in Kobj_idN3_locs:
            imobj_idk = pyfits.getdata(element)
            imobj_idk = imobj_idk[0:2048,0:2048]
            imobj_idk = imobj_idk - bias - dark #because 60 sec exposure
            #imobj_idk = imobj_idk*(norm/np.median(imobj_idk[500:1500,500:1500])) #normalize all flats to have same median as first flat 
            allobj_idKN3.append(imobj_idk)
            print np.shape(allobj_idKN3), np.median(imobj_idk)
        
        print "median combining K_N3..."
        Kobj_idN3 = np.median(allobj_idKN3, axis = 0)
        k_sky3 = pyfits.PrimaryHDU(Kobj_idN3)
        k_sky3.writeto('Calibs/sky/sky_'+str(obj_id)+'K_N3.fits', clobber = True)
        
    
    
    
