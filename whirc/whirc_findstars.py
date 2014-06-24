# Rejects Cosmic Rays
# Must first move cosmics.py to dataredux/whirc
# Then finds coordinates of stars
# Uses the cosmics.py module, courtesy of P. G. van Dokkum, 2001, PASP, 113, 1420 
# availabe for download at  http://www.astro.yale.edu/dokkum/lacosmic/

def whirc_findstars(obj_id): 

    from astropy.convolution import convolve, Tophat2DKernel, AiryDisk2DKernel, MexicanHat2DKernel, Box2DKernel
    import astropy.coordinates as coord
    import astropy.units as u
    import scipy.ndimage as snd 
    from scipy.ndimage import interpolation as interp 
    import pyfits
    import numpy as np
    import pdb
    import pickle
    
    #READ IN DATA AND OBSERVING LOG
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
    
    # Change RA and Dec stuff to degrees
    for line in range(len(aobj_id)):
        raoff = aobj_id[line][9]
        decoff = aobj_id[line][10]
        roff = coord.Angle(raoff,unit=u.hour)
        doff = coord.Angle(decoff, unit=u.degree)
        raoff = (roff.degree)
        decoff = doff.degree
        if np.abs(raoff) > 1:
            raoff = -(360 - raoff)
        aobj_id[line][9] = raoff
        aobj_id[line][10] = decoff        
        #print aobj_id[line][9], aobj_id[line][10]
    
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
        Jobj_idN1_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'JN1_cos'+str(line + 1)+'.fits')
        
    #WRITE PATH TO K_FILES_N1
    Kobj_idN1_locs = []
    for line in range(len(aobj_idKN1)):
        Kobj_idN1_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'KN1_cos'+str(line + 1)+'.fits')
        
    #WRITE PATH OF J_FILES_N2
    Jobj_idN2_locs = []
    for line in range(len(aobj_idJN2)):
        Jobj_idN2_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'JN2_cos'+str(line + 1)+'.fits')
    
    #WRITE PATH OF K_FILES_N2
    Kobj_idN2_locs = []
    for line in range(len(aobj_idKN2)):
        Kobj_idN2_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'KN2_cos'+str(line + 1)+'.fits')
        
    #WRITE PATH OF J_FILES_N3
    Jobj_idN3_locs = []
    for line in range(len(aobj_idJN3)):
        Jobj_idN3_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'JN3_cos'+str(line + 1)+'.fits')
    
    #WRITE PATH OF K_FILES_N3
    Kobj_idN3_locs = []
    for line in range(len(aobj_idKN3)):
        Kobj_idN3_locs.append('Calibs/reduced/'+str(obj_id)+'_'+'KN3_cos'+str(line + 1)+'.fits')
        
    #TELL FINDSTARS HOW MANY STARS TO LOOK FOR:
    if obj_id == 55500:
        num_star = 13
    elif obj_id == 77610:
        num_star = 6
    elif obj_id == 54655:
        num_star = 10
    elif obj_id == 119887:
        num_star = 10
    elif obj_id == 67565:
        num_star = 11
    elif obj_id == 36363:
        num_star = 10
    elif obj_id == 120659:
        num_star = 11
    elif obj_id == 37836: #But this is useless anyway...
        num_star = 7
    elif obj_id == 51306:
        num_star = 8 
    elif obj_id == 122277:
        num_star = 8
    elif obj_id == 20700:
        num_star = 16
    elif obj_id == 57476:
        num_star = 6
    elif obj_id == 38329:
        num_star = 7
    elif obj_id == 121130:
        num_star = 13
    else:
        num_star = 12
    
    
        
    #J_filter, Night 1
    if len(aobj_idJN1) > 0: 
        print "processing J_N1..."
    
        tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
        star_coords = []
        centrals = []
        
        #### MORE TESTING ####
        
        # FIND REFERENCE IMAGE
        centrals = []
        for line in range(len(aobj_idJN1)):
            if aobj_idJN1[line][9] == 0 and aobj_idJN1[line][10] == 0:
                centrals.append(line)
        
        if len(centrals) == 0: 
            centrals.append(4) #in case no object has zero offset, align all objects to fifth object
            
        print 'centrals: ', centrals
        
        #### /MORE TESTING ######
        
        for line in range(len(aobj_idJN1)):
            
            if line in centrals:
                img = pyfits.getdata(Jobj_idN1_locs[line])
        
                #Smoothing
                smooth = convolve(img, tophat)
                threshold = np.mean(smooth)+20 # bit arbitrary, but it works
                print line, threshold
                labels, num = snd.label(smooth > threshold, np.ones((3,3)))
                #find centroids of peaks
                centers = snd.center_of_mass(smooth, labels, range(1, num+1))
                #print centers
        
                cens = []
                badcens = []
                for line in range(len(centers)):
                    if centers[line][0] > 2047 or centers[line][1] > 2047:
                        badcens.append(line)
            
                if len(badcens) > 0:
                    for line in range(len(badcens)-1,-1,-1): #looping backwards
                        centers.pop(badcens[line])
            
                for line in range(len(centers)):
                    cens.append(np.array([[centers[line][0], centers[line][1]],smooth[centers[line][0],centers[line][1]]]))
            
                badpix = []
                for line in range(len(cens)):
                    if 256 < cens[line][0][0] < 276 and 326 < cens[line][0][1] < 346:
                        badpix.append(line)
                    
                if len(badpix) > 0:
                    del cens[badpix[0]]
                        
                cens = np.array(cens)
                #sort in descending order according to 3rd column
                cens = cens[cens[:,1].argsort()[::-1]]
            
                #only include central stars (that will be in all exposures)
                good_cens = []
                for line in range(len(cens)):
                
                        #remove edges
                    #if 50 < cens[line][0][0] < 2000 and 50 < cens[line][0][1] < 2000:
                    good_cens.append(cens[line])
                
                cens = np.array(good_cens)[0:num_star] #remove fainter stars
                    
                cens = cens[:,0]
                print cens
                cens = cens.tolist()

            
        
                star_coords.append(cens)
                print np.shape(star_coords)
                
            #if line =/= center
            else: 
                img = pyfits.getdata(Jobj_idN1_locs[line])
        
                #Smoothing
                smooth = convolve(img, tophat)
                threshold = np.mean(smooth)+20 # bit arbitrary, but it works
                print line, threshold
                labels, num = snd.label(smooth > threshold, np.ones((3,3)))
                #find centroids of peaks
                centers = snd.center_of_mass(smooth, labels, range(1, num+1))
                #print centers
        
                cens = []
                badcens = []
                for line in range(len(centers)):
                    if centers[line][0] > 2047 or centers[line][1] > 2047:
                        badcens.append(line)
            
                if len(badcens) > 0:
                    for line in range(len(badcens)-1,-1,-1): #looping backwards
                        centers.pop(badcens[line])
            
                for line in range(len(centers)):
                    cens.append(np.array([[centers[line][0], centers[line][1]],smooth[centers[line][0],centers[line][1]]]))
            
                badpix = []
                for line in range(len(cens)):
                    if 256 < cens[line][0][0] < 276 and 326 < cens[line][0][1] < 346:
                        badpix.append(line)
                    
                if len(badpix) > 0:
                    del cens[badpix[0]]
                        
                cens = np.array(cens)
                #sort in descending order according to 3rd column
                cens = cens[cens[:,1].argsort()[::-1]]
            
                #only include central stars (that will be in all exposures)
                good_cens = []
                for line in range(len(cens)):
                
                        #remove edges
                    #if 320 < cens[line][0][0] < 1740 and 320 < cens[line][0][1] < 1740:
                    good_cens.append(cens[line])
                
                cens = np.array(good_cens)[0:num_star] #remove fainter stars
                    
                cens = cens[:,0]
                print cens
                cens = cens.tolist()

            
        
                star_coords.append(cens)
                print np.shape(star_coords)
                
            
            #print star_coords
            #print star_coords
            #if aobj_idJN1[line][9] == 0 and aobj_idJN1[line][10] == 0:
            #    print line
            #    centrals.append(line)
            
        pickle.dump(star_coords, open("Calibs/starfields/star_coords"+str(obj_id)+'JN1.p',"wb"))
        #pdb.set_trace()
    
            
        #centrals.append(4) #in case no object has zero offset, align all objects to fifth object
        #center = centrals[0]
        #print star_coords
            
        #xy_offsets_JN1 = []
        #for line in range(len(star_coords)):
        #    offset = [star_coords[center][0][0] - star_coords[line][0][0], star_coords[center][0][1] - star_coords[line][0][1]]
        #    xy_offsets_JN1.append(offset)
        #print xy_offsets_JN1
    
    
    #The next bit is currently not useful but might be later. 
    '''
    print "reading in JN1..."
    if len(aobj_idJN1) > 0:
        shifteds = []
        for line in range(len(aobj_idJN1)):
            sci = pyfits.getdata(Jobj_idN1_locs[line])
            shifted = interp.shift(sci, xy_offsets_JN1[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = sum(shifteds)
        median = np.median(shifteds, axis = 0)
    
    
        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_JN1_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_JN1_median.fits',clobber=True)
        
    
    
    #K_filter, Night 1
    print "processing K_N1..."
    
    tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
    star_coords = []
    centrals = []
    for line in range(len(aobj_idKN1)):
        img = pyfits.getdata(Kobj_idN1_locs[line])
        
        #Smoothing
        smooth = convolve(img, tophat)
        threshold = 0.99*np.max(smooth) # bit arbitrary, but it works
        print threshold 
        labels, num = snd.label(smooth > threshold, np.ones((3,3)))
        #find centroids of peaks
        centers = snd.center_of_mass(smooth, labels, range(1, num+1))
        star_coords.append(centers)
        if aobj_idKN1[line][9] == 0 and aobj_idKN1[line][10] == 0:
            print line, centers
            centrals.append(line)
            
    centrals.append(4) #in case no object has zero offset, align all objects to fifth object
    center = centrals[0]
    print star_coords
            
    xy_offsets_KN1 = []
    for line in range(len(star_coords)):
        offset = [star_coords[center][0][0] - star_coords[line][0][0], star_coords[center][0][1] - star_coords[line][0][1]]
        xy_offsets_KN1.append(offset)
    print xy_offsets_KN1
    
    print "reading in KN1..."
    if len(aobj_idKN1) > 0:
        shifteds = []
        for line in range(len(aobj_idKN1)):
            sci = pyfits.getdata(Kobj_idN1_locs[line])
            shifted = interp.shift(sci, xy_offsets_KN1[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = sum(shifteds)
        median = np.median(shifteds, axis = 0)
    
    
        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_KN1_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_KN1_median.fits',clobber=True)
    
    
    #J_filter, Night 2
    print "processing J_N2..."
    
    tophat = Tophat2DKernel(4) #Radius can be changed to increase sharpess 
    
    star_coords = []
    centrals = []
    for line in range(len(aobj_idJN2)):
        img = pyfits.getdata(Jobj_idN2_locs[line])
        
        #Smoothing
        smooth = convolve(img, tophat)
        threshold = 0.99*np.max(smooth) # bit arbitrary, but it works
        print threshold 
        labels, num = snd.label(smooth > threshold, np.ones((3,3)))
        #find centroids of peaks
        centers = snd.center_of_mass(smooth, labels, range(1, num+1))
        star_coords.append(centers)
        if aobj_idJN2[line][9] == 0 and aobj_idJN2[line][10] == 0:
            print line
            centrals.append(line)
            
    centrals.append(4) #in case no object has zero offset, align all objects to fifth object
    center = centrals[0]
    print star_coords
            
    xy_offsets_JN2 = []
    for line in range(len(star_coords)):
        offset = [star_coords[center][0][0] - star_coords[line][0][0], star_coords[center][0][1] - star_coords[line][0][1]]
        xy_offsets_JN2.append(offset)
    print xy_offsets_JN2
    
    print "reading in JN2..."
    if len(aobj_idJN2) > 0:
        shifteds = []
        for line in range(len(aobj_idJN2)):
            sci = pyfits.getdata(Jobj_idN2_locs[line])
            shifted = interp.shift(sci, xy_offsets_JN2[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = sum(shifteds)
        median = np.median(shifteds, axis = 0)
    
    
        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_JN2_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_JN2_median.fits',clobber=True)
        
        
    #K_filter, Night 2
    print "processing K_N2..."
    
    tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
    star_coords = []
    centrals = []
    for line in range(len(aobj_idKN2)):
        img = pyfits.getdata(Kobj_idN2_locs[line])
        
        #Smoothing
        smooth = convolve(img, tophat)
        threshold = 0.99*np.max(smooth) # bit arbitrary, but it works
        print threshold 
        labels, num = snd.label(smooth > threshold, np.ones((3,3)))
        #find centroids of peaks
        centers = snd.center_of_mass(smooth, labels, range(1, num+1))
        star_coords.append(centers)
        if aobj_idKN2[line][9] == 0 and aobj_idKN2[line][10] == 0:
            print line
            centrals.append(line)
            
    centrals.append(4) #in case no object has zero offset, align all objects to fifth object
    center = centrals[0]
    print star_coords
            
    xy_offsets_KN2 = []
    for line in range(len(star_coords)):
        offset = [star_coords[center][0][0] - star_coords[line][0][0], star_coords[center][0][1] - star_coords[line][0][1]]
        xy_offsets_KN2.append(offset)
    print xy_offsets_KN2
    
    print "reading in KN2..."
    if len(aobj_idKN2) > 0:
        shifteds = []
        for line in range(len(aobj_idKN2)):
            sci = pyfits.getdata(Kobj_idN2_locs[line])
            shifted = interp.shift(sci, xy_offsets_KN2[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = sum(shifteds)
        median = np.median(shifteds, axis = 0)
    
    
        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_KN2_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_KN2_median.fits',clobber=True)
        
    #J_filter, Night 3
    print "processing J_N3..."
    
    tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
    star_coords = []
    centrals = []
    for line in range(len(aobj_idJN3)):
        img = pyfits.getdata(Jobj_idN3_locs[line])
        
        #Smoothing
        smooth = convolve(img, tophat)
        threshold = 0.99*np.max(smooth) # bit arbitrary, but it works
        print threshold 
        labels, num = snd.label(smooth > threshold, np.ones((3,3)))
        #find centroids of peaks
        centers = snd.center_of_mass(smooth, labels, range(1, num+1))
        print centers
        star_coords.append(centers)
        if aobj_idJN3[line][9] == 0 and aobj_idJN3[line][10] == 0:
            print line
            centrals.append(line)
            
    centrals.append(5) #in case no object has zero offset, align all objects to fifth object
    center = centrals[0]
    print star_coords
    
            
    xy_offsets_JN3 = []
    for line in range(len(star_coords)):
        offset = [star_coords[center][0][0] - star_coords[line][0][0], star_coords[center][0][1] - star_coords[line][0][1]]
        xy_offsets_JN3.append(offset)
    print xy_offsets_JN3
    
    print "reading in JN3..."
    if len(aobj_idJN3) > 0:
        shifteds = []
        for line in range(len(aobj_idJN3)):
            sci = pyfits.getdata(Jobj_idN3_locs[line])
            shifted = interp.shift(sci, xy_offsets_JN3[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = sum(shifteds)
        median = np.median(shifteds, axis = 0)
    
    
        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_JN3_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_JN3_median.fits',clobber=True)
        
    #K_filter, Night 3
    print "processing K_N3..."
    
    tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
    star_coords = []
    centrals = []
    for line in range(len(aobj_idKN3)):
        img = pyfits.getdata(Kobj_idN3_locs[line])
        
        #Smoothing
        smooth = convolve(img, tophat)
        threshold = 0.99*np.max(smooth) # bit arbitrary, but it works
        print threshold 
        labels, num = snd.label(smooth > threshold, np.ones((3,3)))
        #find centroids of peaks
        centers = snd.center_of_mass(smooth, labels, range(1, num+1))
        star_coords.append(centers)
        if aobj_idKN3[line][9] == 0 and aobj_idKN3[line][10] == 0:
            print line
            centrals.append(line)
    
    centrals.append(4) #in case no object has zero offset, align all objects to fifth object
    center = centrals[0]
    print star_coords
            
    xy_offsets_KN3 = []
    for line in range(len(star_coords)):
        offset = [star_coords[center][0][0] - star_coords[line][0][0], star_coords[center][0][1] - star_coords[line][0][1]]
        xy_offsets_KN3.append(offset)
    print xy_offsets_KN3
    
    print "reading in KN3..."
    if len(aobj_idKN3) > 0:
        shifteds = []
        for line in range(len(aobj_idKN3)):
            sci = pyfits.getdata(Kobj_idN3_locs[line])
            shifted = interp.shift(sci, xy_offsets_KN3[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = sum(shifteds)
        median = np.median(shifteds, axis = 0)
    
    
        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_KN3_combined.fits', clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_KN3_median.fits', clobber=True)
        '''
        
    #Here I continue the useful code:
    #K_filter, Night 1
    if len(aobj_idKN1) > 0:
        print "processing K_N1..."
    
        tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
        star_coords = []
        centrals = []
        for line in range(len(aobj_idKN1)):
            img = pyfits.getdata(Kobj_idN1_locs[line])
        
            #Smoothing
            smooth = convolve(img, tophat)
            threshold = np.mean(smooth)+20 # bit arbitrary, but it works
            print line, threshold
            labels, num = snd.label(smooth > threshold, np.ones((3,3)))
            #find centroids of peaks
            centers = snd.center_of_mass(smooth, labels, range(1, num+1))
            #print centers
        
            cens = []
            badcens = []
            for line in range(len(centers)):
                if centers[line][0] > 2047 or centers[line][1] > 2047:
                    badcens.append(line)
            
            if len(badcens) > 0:
                for line in range(len(badcens)-1,-1,-1): #looping backwards
                    centers.pop(badcens[line])
                    
                    
            for line in range(len(centers)):
                cens.append(np.array([[centers[line][0], centers[line][1]],smooth[centers[line][0],centers[line][1]]]))
            
            cens = np.array(cens)
            #sort in descending order according to 3rd column
            cens = cens[cens[:,1].argsort()[::-1]] 
            
            #only include central stars (that will be in all exposures)
            good_cens = []
            for line in range(len(cens)):
                if 100 < cens[line][0][0] < 1950 and 100 < cens[line][0][1] < 1950:
                    good_cens.append(cens[line])
            cens = np.array(good_cens)[0:num_star] #remove fainter stars
            cens = cens[:,0]
            print "stars", cens
            cens = cens.tolist()
            

            star_coords.append(cens)
            print np.shape(star_coords)
            #print star_coords
            #if aobj_idJN1[line][9] == 0 and aobj_idJN1[line][10] == 0:
            #    print line
            #    centrals.append(line)
            
        pickle.dump(star_coords, open("Calibs/starfields/star_coords"+str(obj_id)+'KN1.p',"wb"))
        
    #J_filter, Night 2
    if len(aobj_idJN2) > 0: 
        print "processing J_N2..."
    
        tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
        star_coords = []
        centrals = []
        for line in range(len(aobj_idJN2)):
            img = pyfits.getdata(Jobj_idN2_locs[line])
        
            #Smoothing
            smooth = convolve(img, tophat)
            threshold = np.mean(smooth)+20 # bit arbitrary, but it works
            print line, threshold
            labels, num = snd.label(smooth > threshold, np.ones((3,3)))
            #find centroids of peaks
            centers = snd.center_of_mass(smooth, labels, range(1, num+1))
            #print centers
        
            cens = []
            badcens = []
            for line in range(len(centers)):
                if centers[line][0] > 2047 or centers[line][1] > 2047:
                    badcens.append(line)
            
            if len(badcens) > 0:
                if len(badcens) > 0:
                    for line in range(len(badcens)-1,-1,-1): #looping backwards
                        centers.pop(badcens[line])
                
            for line in range(len(centers)):
                cens.append(np.array([[centers[line][0], centers[line][1]],smooth[centers[line][0],centers[line][1]]]))
            
            cens = np.array(cens)
            #sort in descending order according to 3rd column
            cens = cens[cens[:,1].argsort()[::-1]] 
            #only include central stars (that will be in all exposures)
            good_cens = []
            for line in range(len(cens)):
                if 100 < cens[line][0][0] < 1950 and 100 < cens[line][0][1] < 1950:
                    good_cens.append(cens[line])
                    
            cens = np.array(good_cens)[0:num_star] #remove fainter stars
            cens = cens[:,0]
            print 'stars: ',cens
            cens = cens.tolist()

            
        
            star_coords.append(cens)
            print np.shape(star_coords)
            #print star_coords
            #print star_coords
            #if aobj_idJN1[line][9] == 0 and aobj_idJN1[line][10] == 0:
            #    print line
            #    centrals.append(line)
            
        pickle.dump(star_coords, open("Calibs/starfields/star_coords"+str(obj_id)+'JN2.p',"wb"))
        
    #K_filter, Night 2
    if len(aobj_idKN2) > 0:
        print "processing K_N2..."
    
        tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
        star_coords = []
        centrals = []
        for line in range(len(aobj_idKN2)):
            img = pyfits.getdata(Kobj_idN2_locs[line])
        
            #Smoothing
            smooth = convolve(img, tophat)
            threshold = np.mean(smooth)+20 # bit arbitrary, but it works
            print line, threshold
            labels, num = snd.label(smooth > threshold, np.ones((3,3)))
            #find centroids of peaks
            centers = snd.center_of_mass(smooth, labels, range(1, num+1))
            #print centers
        
            cens = []
            badcens = []
            for line in range(len(centers)):
                if centers[line][0] > 2047 or centers[line][1] > 2047:
                    badcens.append(line)
            
            if len(badcens) > 0:
                if len(badcens) > 0:
                    for line in range(len(badcens)-1,-1,-1): #looping backwards
                        centers.pop(badcens[line])
                
            for line in range(len(centers)):
                cens.append(np.array([[centers[line][0], centers[line][1]],smooth[centers[line][0],centers[line][1]]]))
            
            cens = np.array(cens)
            #sort in descending order according to 3rd column
            cens = cens[cens[:,1].argsort()[::-1]] 
            #only include central stars (that will be in all exposures)
            good_cens = []
            for line in range(len(cens)):
                if 100 < cens[line][0][0] < 1950 and 100 < cens[line][0][1] < 1950:
                    good_cens.append(cens[line])
                    
            cens = np.array(good_cens)[0:num_star] #remove fainter stars
            cens = cens[:,0]
            print 'stars: ', cens
            cens = cens.tolist()

            star_coords.append(cens)
            print np.shape(star_coords)
            #print star_coords
            #print star_coords
            #if aobj_idJN1[line][9] == 0 and aobj_idJN1[line][10] == 0:
            #    print line
            #    centrals.append(line)
            
        pickle.dump(star_coords, open("Calibs/starfields/star_coords"+str(obj_id)+'KN2.p',"wb"))
        
    #J_filter, Night 3
    if len(aobj_idJN3) > 0: 
        print "processing J_N3..."
    
        tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
        star_coords = []
        centrals = []
        for line in range(len(aobj_idJN3)):
            img = pyfits.getdata(Jobj_idN3_locs[line])
        
            #Smoothing
            smooth = convolve(img, tophat)
            threshold = np.mean(smooth)+20 # bit arbitrary, but it works
            print line, threshold
            labels, num = snd.label(smooth > threshold, np.ones((3,3)))
            #find centroids of peaks
            centers = snd.center_of_mass(smooth, labels, range(1, num+1))
            #print centers
        
            cens = []
            badcens = []
            for line in range(len(centers)):
                if centers[line][0] > 2047 or centers[line][1] > 2047:
                    badcens.append(line)
            
            if len(badcens) > 0:
                if len(badcens) > 0:
                    for line in range(len(badcens)-1,-1,-1): #looping backwards
                        centers.pop(badcens[line])
            
            for line in range(len(centers)):
                cens.append(np.array([[centers[line][0], centers[line][1]],smooth[centers[line][0],centers[line][1]]]))
            
            cens = np.array(cens)
            #sort in descending order according to 3rd column
            cens = cens[cens[:,1].argsort()[::-1]] 
            #only include central stars (that will be in all exposures)
            good_cens = []
            for line in range(len(cens)):
                if 100 < cens[line][0][0] < 1950 and 100 < cens[line][0][1] < 1950:
                    good_cens.append(cens[line])
                    
            cens = np.array(good_cens)[0:num_star] #remove fainter stars
            cens = cens[:,0]
            print 'stars: ', cens
            cens = cens.tolist()

            
        
            star_coords.append(cens)
            print np.shape(star_coords)
            #print star_coords
            #print star_coords
            #if aobj_idJN1[line][9] == 0 and aobj_idJN1[line][10] == 0:
            #    print line
            #    centrals.append(line)
            
        pickle.dump(star_coords, open("Calibs/starfields/star_coords"+str(obj_id)+'JN3.p',"wb"))
        
    #K_filter, Night 3
    if len(aobj_idKN3) > 0: 
        print "processing K_N3..."
    
        tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
        star_coords = []
        centrals = []
        for line in range(len(aobj_idKN3)):
            img = pyfits.getdata(Kobj_idN3_locs[line])
        
            #Smoothing
            smooth = convolve(img, tophat)
            threshold = np.mean(smooth)+20 # bit arbitrary, but it works
            print line, threshold
            labels, num = snd.label(smooth > threshold, np.ones((3,3)))
            #find centroids of peaks
            centers = snd.center_of_mass(smooth, labels, range(1, num+1))
            #print centers
        
            cens = []
            badcens = []
            for line in range(len(centers)):
                if centers[line][0] > 2047 or centers[line][0] < 0 or centers[line][1] <0 or centers[line][1] > 2047:
                    badcens.append(line)
            
            if len(badcens) > 0:
               if len(badcens) > 0:
                   for line in range(len(badcens)-1,-1,-1): #looping backwards
                       centers.pop(badcens[line])
               
            for line in range(len(centers)):
                cens.append(np.array([[centers[line][0], centers[line][1]],smooth[centers[line][0],centers[line][1]]]))
            
            cens = np.array(cens)
            #sort in descending order according to 3rd column
            cens = cens[cens[:,1].argsort()[::-1]] 
            #only include central stars (that will be in all exposures)
            good_cens = []
            for line in range(len(cens)):
                if 100 < cens[line][0][0] < 1950 and 100 < cens[line][0][1] < 1950:
                    good_cens.append(cens[line])
                    
            cens = np.array(good_cens)[0:num_star] #remove fainter stars
            cens = cens[:,0]
            print 'stars: ', cens
            cens = cens.tolist()

            
        
            star_coords.append(cens)
            print np.shape(star_coords)
            #print star_coords
            #print star_coords
            #if aobj_idJN1[line][9] == 0 and aobj_idJN1[line][10] == 0:
            #    print line
            #    centrals.append(line)
            
        pickle.dump(star_coords, open("Calibs/starfields/star_coords"+str(obj_id)+'KN3.p',"wb"))