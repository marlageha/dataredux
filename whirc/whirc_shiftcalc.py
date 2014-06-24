# READS IN STARFIELD FROM PICKLE FILE 
# CALCULATES THE NECESSARY SHIFT TO ALLIGN STARS FOR TO A REFERRENCE IMAGE

def whirc_shiftcalc(obj_id):
    
    import pyfits
    import numpy as np
    import pdb
    import pickle
    import astropy.coordinates as coord
    import astropy.units as u
    import scipy.ndimage as snd 
    from scipy.ndimage import interpolation as interp 
    import pyfits
    import matplotlib.pyplot as plt
    
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
        
    
    #READ IN STARFIELDS FROM PICKLE FILE
    print "reading in JN1"
    
    if len(aobj_idJN1) > 0:
        
        
        xy_shifts_JN1 = []
                
        # FIND REFERENCE IMAGE
        centrals = []
        for line in range(len(aobj_idJN1)):
            if aobj_idJN1[line][9] == 0 and aobj_idJN1[line][10] == 0:
                centrals.append(line)
        
        if len(centrals) == 0: 
            centrals.append(4) #in case no object has zero offset, align all objects to fifth object
        center = centrals[0]
        
        # READ IN AND INITIALIZE DISTANCE MINIMIZATION
        star_coords = pickle.load(open( "Calibs/starfields/star_coords"+str(obj_id)+'JN1.p', 'rb'))
        
        ref_im = star_coords[center]
        
        def distance(x0,y0,x1,y1):
            return np.sqrt((x0-x1)**2+(y0-y1)**2)
            
        #LOOP THROUGH EACH EXPOSURE IN JN1
        for shf_num in range(len(aobj_idJN1)):
            
            shift_im = np.array(star_coords[shf_num]) 
            
            #----------testing possible xrange to improve time------------#
            xranges = []
            yranges = []
            for sline in range(len(shift_im)):
                for rline in range(len(ref_im)):
                    xranges.append(float(shift_im[sline][0]-ref_im[rline][0]))
                    yranges.append(float(shift_im[sline][1]-ref_im[rline][1]))
                                
            #pdb.set_trace()
            #x_shift = np.array(xranges)
            #y_shift = np.array(yranges)
                    
                
            
            #---------------end testing ---------------------------------------#
            
            #FIRST ITERATION, TO GET ROUGH APPROXIMATION
            x_shift = np.arange(-1000, 1000, 20)
            y_shift = np.arange(-1000, 1000, 20)
            datai = np.zeros((len(x_shift),len(y_shift)))
            
            #Loop through different x and y shifts
            for x in range(len(x_shift)):
               for y in range(len(y_shift)):
                   tot_points = []
                   
                   #Apply shift
                   shift_im = np.array(star_coords[shf_num]) + [x_shift[x],y_shift[y]]
                   dist = dict() #creates a separate variable for each star in shift_im
            
                   #loop through stars in image to shift
                   for shift_line in range(len(shift_im)):
                       dist[shift_line] = [] 
                       
                       #loop through stars in reference image
                       for ref_line in range(len(ref_im)):
                           dist[shift_line].append(distance(ref_im[ref_line][0], ref_im[ref_line][1], shift_im[shift_line][0], shift_im[shift_line][1]))
                           
                       #find distance to closest star in ref_im for each star in shift_im
                       dist[shift_line] = np.min(dist[shift_line])
                       
                       if dist[shift_line] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift[mins[0]]),np.median(y_shift[mins[1]])]
            #print mina
       
            #pdb.set_trace()
            
            #-----------------------------Now Repeat--------------------------------#
            
            #SECOND ITERATION, TO GET BETTER APPROXIMATION       
            x_shift2 = np.arange(float(mina[0])-150, float(mina[0])+150, 10)
            y_shift2 = np.arange(float(mina[1])-150, float(mina[1])+150, 10)
            
            datai = np.zeros((len(x_shift2),len(y_shift2)))

            for x in range(len(x_shift2)):
               for y in range(len(y_shift2)):
                   
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift2[x], y_shift2[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       
                       if dist[shift_star] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift2[mins[0]]),np.median(y_shift2[mins[1]])]
            #print [float(mina[0]),float(mina[1])]

            #THIRD ITERATION, TO GET BETTER APPROXIMATION
            x_shift3 = np.arange(float(mina[0])-20, float(mina[0])+20, 1)
            y_shift3 = np.arange(float(mina[1])-20, float(mina[1])+20, 1)
            datai = np.zeros((len(x_shift3),len(y_shift3)))

            for x in range(len(x_shift3)):
               for y in range(len(y_shift3)):
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift3[x],y_shift3[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       if dist[shift_star] < 4:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift3[mins[0]]),np.median(y_shift3[mins[1]])]
            print mina
            
            #Append optimal x,y shift for each exposure 
            xy_shifts_JN1.append([float(mina[0]),float(mina[1])])
            
        print xy_shifts_JN1
        
        #APPLY SHIFTS TO DATA
        shifteds = []
        for line in range(len(aobj_idJN1)):
            sci = pyfits.getdata(Jobj_idN1_locs[line])
            shifted = interp.shift(sci, xy_shifts_JN1[line] ,order = 0)
            #shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = np.mean(shifteds, axis = 0) #trying instead of sum
        #for i in shifteds:
            
        
        ##### MAKING PLOTS #########
        '''
        boxs = [shifteds[a][500:1500] for a in range(len(Jobj_idN1_locs))]
        meds = [np.median(a) for a in boxs]
        means = [np.mean(a) for a in boxs]
        fig = plt.subplot(1,1,1)
        xvar = range(len(Jobj_idN1_locs))
        plt.plot(xvar, meds, 'o', color = 'blue', label = 'median')
        plt.plot(xvar, means, 'o', color = 'red', label = 'mean')
        fig.legend(loc = 2)
        fig.set_xlabel('exposure number')
        fig.set_ylabel('counts')
        
        
        pdb.set_trace()
        '''
        
        median = np.median(shifteds, axis = 0)


        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_JN1_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_JN1_median.fits',clobber=True)
        
        #--------------------------------------------------------------------------#
        
    print "reading in KN1"
    
    if len(aobj_idKN1) > 0:
        
        
        xy_shifts_KN1 = []
                
        # FIND REFERENCE IMAGE
        centrals = []
        for line in range(len(aobj_idKN1)):
            if aobj_idKN1[line][9] == 0 and aobj_idKN1[line][10] == 0:
                centrals.append(line)
        
        if len(centrals) == 0: 
            centrals.append(4) #in case no object has zero offset, align all objects to fifth object
        center = centrals[0]
        
        # READ IN AND INITIALIZE DISTANCE MINIMIZATION
        star_coords = pickle.load(open( "Calibs/starfields/star_coords"+str(obj_id)+'KN1.p', 'rb'))
        
        ref_im = star_coords[center]
        
        def distance(x0,y0,x1,y1):
            return np.sqrt((x0-x1)**2+(y0-y1)**2)
            
        #LOOP THROUGH EACH EXPOSURE IN JN1
        for shf_num in range(len(aobj_idKN1)):
            
            shift_im = np.array(star_coords[shf_num]) 
            
            #----------testing possible xrange to improve time------------#
            xranges = []
            yranges = []
            for sline in range(len(shift_im)):
                for rline in range(len(ref_im)):
                    xranges.append(float(shift_im[sline][0]-ref_im[rline][0]))
                    yranges.append(float(shift_im[sline][1]-ref_im[rline][1]))
                                
            #pdb.set_trace()
            #x_shift = np.array(xranges)
            #y_shift = np.array(yranges)
                    
                
            
            #---------------end testing ---------------------------------------#
            
            #FIRST ITERATION, TO GET ROUGH APPROXIMATION
            x_shift = np.arange(-1000, 1000, 20)
            y_shift = np.arange(-1000, 1000, 20)
            datai = np.zeros((len(x_shift),len(y_shift)))
            
            #Loop through different x and y shifts
            for x in range(len(x_shift)):
               for y in range(len(y_shift)):
                   tot_points = []
                   
                   #Apply shift
                   shift_im = np.array(star_coords[shf_num]) + [x_shift[x],y_shift[y]]
                   dist = dict() #creates a separate variable for each star in shift_im
            
                   #loop through stars in image to shift
                   for shift_line in range(len(shift_im)):
                       dist[shift_line] = [] 
                       
                       #loop through stars in reference image
                       for ref_line in range(len(ref_im)):
                           dist[shift_line].append(distance(ref_im[ref_line][0], ref_im[ref_line][1], shift_im[shift_line][0], shift_im[shift_line][1]))
                           
                       #find distance to closest star in ref_im for each star in shift_im
                       dist[shift_line] = np.min(dist[shift_line])
                       
                       if dist[shift_line] < 10:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift[mins[0]]),np.median(y_shift[mins[1]])]
            #print mina
       
            #pdb.set_trace()
            
            #-----------------------------Now Repeat--------------------------------#
            
            #SECOND ITERATION, TO GET BETTER APPROXIMATION       
            x_shift2 = np.arange(float(mina[0])-150, float(mina[0])+150, 10)
            y_shift2 = np.arange(float(mina[1])-150, float(mina[1])+150, 10)
            
            datai = np.zeros((len(x_shift2),len(y_shift2)))

            for x in range(len(x_shift2)):
               for y in range(len(y_shift2)):
                   
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift2[x], y_shift2[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       
                       if dist[shift_star] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift2[mins[0]]),np.median(y_shift2[mins[1]])]
            #print [float(mina[0]),float(mina[1])]

            #THIRD ITERATION, TO GET BETTER APPROXIMATION
            x_shift3 = np.arange(float(mina[0])-20, float(mina[0])+20, 1)
            y_shift3 = np.arange(float(mina[1])-20, float(mina[1])+20, 1)
            datai = np.zeros((len(x_shift3),len(y_shift3)))

            for x in range(len(x_shift3)):
               for y in range(len(y_shift3)):
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift3[x],y_shift3[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       if dist[shift_star] < 4:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift3[mins[0]]),np.median(y_shift3[mins[1]])]
            print mina
            
            #Append optimal x,y shift for each exposure 
            xy_shifts_KN1.append([float(mina[0]),float(mina[1])])
            
        print xy_shifts_KN1
        
        #APPLY SHIFTS TO DATA
        shifteds = []
        for line in range(len(aobj_idKN1)):
            sci = pyfits.getdata(Kobj_idN1_locs[line])
            shifted = interp.shift(sci, xy_shifts_KN1[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = sum(shifteds)
        median = np.median(shifteds, axis = 0)


        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_KN1_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_KN1_median.fits',clobber=True)

    
    #---------------------------------------------------------------------------------#
            
    print "reading in JN2"
    
    if len(aobj_idJN2) > 0:
        
        xy_shifts_JN2 = []
                
        # FIND REFERENCE IMAGE
        centrals = []
        for line in range(len(aobj_idJN2)):
            if aobj_idJN2[line][9] == 0 and aobj_idJN2[line][10] == 0:
                centrals.append(line)
        
        if len(centrals) == 0: 
            centrals.append(4) #in case no object has zero offset, align all objects to fifth object
        center = centrals[0]
        
        # READ IN AND INITIALIZE DISTANCE MINIMIZATION
        star_coords = pickle.load(open( "Calibs/starfields/star_coords"+str(obj_id)+'JN2.p', 'rb'))
        
        ref_im = star_coords[center]
        
        def distance(x0,y0,x1,y1):
            return np.sqrt((x0-x1)**2+(y0-y1)**2)
            
        #LOOP THROUGH EACH EXPOSURE IN JN1
        for shf_num in range(len(aobj_idJN2)):
            
            shift_im = np.array(star_coords[shf_num]) 
            
            #----------testing possible xrange to improve time------------#
            xranges = []
            yranges = []
            for sline in range(len(shift_im)):
                for rline in range(len(ref_im)):
                    xranges.append(float(shift_im[sline][0]-ref_im[rline][0]))
                    yranges.append(float(shift_im[sline][1]-ref_im[rline][1]))
                                
            #pdb.set_trace()
            #x_shift = np.array(xranges)
            #y_shift = np.array(yranges)
                    
                
            
            #---------------end testing ---------------------------------------#
            
            #FIRST ITERATION, TO GET ROUGH APPROXIMATION
            x_shift = np.arange(-1000, 1000, 20)
            y_shift = np.arange(-1000, 1000, 20)
            datai = np.zeros((len(x_shift),len(y_shift)))
            
            #Loop through different x and y shifts
            for x in range(len(x_shift)):
               for y in range(len(y_shift)):
                   tot_points = []
                   
                   #Apply shift
                   shift_im = np.array(star_coords[shf_num]) + [x_shift[x],y_shift[y]]
                   dist = dict() #creates a separate variable for each star in shift_im
            
                   #loop through stars in image to shift
                   for shift_line in range(len(shift_im)):
                       dist[shift_line] = [] 
                       
                       #loop through stars in reference image
                       for ref_line in range(len(ref_im)):
                           dist[shift_line].append(distance(ref_im[ref_line][0], ref_im[ref_line][1], shift_im[shift_line][0], shift_im[shift_line][1]))
                           
                       #find distance to closest star in ref_im for each star in shift_im
                       dist[shift_line] = np.min(dist[shift_line])
                       
                       if dist[shift_line] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift[mins[0]]),np.median(y_shift[mins[1]])]
            #print mina
       
            #pdb.set_trace()
            
            #-----------------------------Now Repeat--------------------------------#
            
            #SECOND ITERATION, TO GET BETTER APPROXIMATION       
            x_shift2 = np.arange(float(mina[0])-150, float(mina[0])+150, 10)
            y_shift2 = np.arange(float(mina[1])-150, float(mina[1])+150, 10)
            
            datai = np.zeros((len(x_shift2),len(y_shift2)))

            for x in range(len(x_shift2)):
               for y in range(len(y_shift2)):
                   
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift2[x], y_shift2[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       
                       if dist[shift_star] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift2[mins[0]]),np.median(y_shift2[mins[1]])]
            #print [float(mina[0]),float(mina[1])]

            #THIRD ITERATION, TO GET BETTER APPROXIMATION
            x_shift3 = np.arange(float(mina[0])-20, float(mina[0])+20, 1)
            y_shift3 = np.arange(float(mina[1])-20, float(mina[1])+20, 1)
            datai = np.zeros((len(x_shift3),len(y_shift3)))

            for x in range(len(x_shift3)):
               for y in range(len(y_shift3)):
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift3[x],y_shift3[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       if dist[shift_star] < 4:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift3[mins[0]]),np.median(y_shift3[mins[1]])]
            print mina
            
            #Append optimal x,y shift for each exposure 
            xy_shifts_JN2.append([float(mina[0]),float(mina[1])])
            
        print xy_shifts_JN2
        
        #APPLY SHIFTS TO DATA
        shifteds = []
        for line in range(len(aobj_idJN2)):
            sci = pyfits.getdata(Jobj_idN2_locs[line])
            shifted = interp.shift(sci, xy_shifts_JN2[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = np.mean(shifteds, axis =0) #trying instead of sum
        median = np.median(shifteds, axis = 0)


        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_JN2_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_JN2_median.fits',clobber=True)
        
        #--------------------------------------------------------------------------#
        
    print "reading in KN2"
    
    if len(aobj_idKN2) > 0:
        
        xy_shifts_KN2 = []
                
        # FIND REFERENCE IMAGE
        centrals = []
        for line in range(len(aobj_idKN2)):
            if aobj_idKN2[line][9] == 0 and aobj_idKN2[line][10] == 0:
                centrals.append(line)
        
        if len(centrals) == 0: 
            centrals.append(4) #in case no object has zero offset, align all objects to fifth object
        center = centrals[0]
        
        # READ IN AND INITIALIZE DISTANCE MINIMIZATION
        star_coords = pickle.load(open( "Calibs/starfields/star_coords"+str(obj_id)+'KN2.p', 'rb'))
        
        ref_im = star_coords[center]
        
        def distance(x0,y0,x1,y1):
            return np.sqrt((x0-x1)**2+(y0-y1)**2)
            
        #LOOP THROUGH EACH EXPOSURE IN JN1
        for shf_num in range(len(aobj_idKN2)):
            
            shift_im = np.array(star_coords[shf_num]) 
            
            #----------testing possible xrange to improve time------------#
            xranges = []
            yranges = []
            for sline in range(len(shift_im)):
                for rline in range(len(ref_im)):
                    xranges.append(float(shift_im[sline][0]-ref_im[rline][0]))
                    yranges.append(float(shift_im[sline][1]-ref_im[rline][1]))
                                
            #pdb.set_trace()
            #x_shift = np.array(xranges)
            #y_shift = np.array(yranges)
                    
                
            
            #---------------end testing ---------------------------------------#
            
            #FIRST ITERATION, TO GET ROUGH APPROXIMATION
            x_shift = np.arange(-1000, 1000, 20)
            y_shift = np.arange(-1000, 1000, 20)
            datai = np.zeros((len(x_shift),len(y_shift)))
            
            #Loop through different x and y shifts
            for x in range(len(x_shift)):
               for y in range(len(y_shift)):
                   tot_points = []
                   
                   #Apply shift
                   shift_im = np.array(star_coords[shf_num]) + [x_shift[x],y_shift[y]]
                   dist = dict() #creates a separate variable for each star in shift_im
            
                   #loop through stars in image to shift
                   for shift_line in range(len(shift_im)):
                       dist[shift_line] = [] 
                       
                       #loop through stars in reference image
                       for ref_line in range(len(ref_im)):
                           dist[shift_line].append(distance(ref_im[ref_line][0], ref_im[ref_line][1], shift_im[shift_line][0], shift_im[shift_line][1]))
                           
                       #find distance to closest star in ref_im for each star in shift_im
                       dist[shift_line] = np.min(dist[shift_line])
                       
                       if dist[shift_line] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift[mins[0]]),np.median(y_shift[mins[1]])]
            #print mina
       
            #pdb.set_trace()
            
            #-----------------------------Now Repeat--------------------------------#
            
            #SECOND ITERATION, TO GET BETTER APPROXIMATION       
            x_shift2 = np.arange(float(mina[0])-150, float(mina[0])+150, 10)
            y_shift2 = np.arange(float(mina[1])-150, float(mina[1])+150, 10)
            
            datai = np.zeros((len(x_shift2),len(y_shift2)))

            for x in range(len(x_shift2)):
               for y in range(len(y_shift2)):
                   
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift2[x], y_shift2[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       
                       if dist[shift_star] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift2[mins[0]]),np.median(y_shift2[mins[1]])]
            #print [float(mina[0]),float(mina[1])]

            #THIRD ITERATION, TO GET BETTER APPROXIMATION
            x_shift3 = np.arange(float(mina[0])-20, float(mina[0])+20, 1)
            y_shift3 = np.arange(float(mina[1])-20, float(mina[1])+20, 1)
            datai = np.zeros((len(x_shift3),len(y_shift3)))

            for x in range(len(x_shift3)):
               for y in range(len(y_shift3)):
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift3[x],y_shift3[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       if dist[shift_star] < 4:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift3[mins[0]]),np.median(y_shift3[mins[1]])]
            print mina
            
            #Append optimal x,y shift for each exposure 
            xy_shifts_KN2.append([float(mina[0]),float(mina[1])])
            
        print xy_shifts_KN2
        
        #APPLY SHIFTS TO DATA
        shifteds = []
        for line in range(len(aobj_idKN2)):
            sci = pyfits.getdata(Kobj_idN2_locs[line])
            shifted = interp.shift(sci, xy_shifts_KN2[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = sum(shifteds)
        median = np.median(shifteds, axis = 0)


        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_KN2_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_KN2_median.fits',clobber=True)
        
    #--------------------------------------------------------------------------------#
            
    print "reading in JN3"
    
    if len(aobj_idJN3) > 0:
        
        xy_shifts_JN3 = []
                
        # FIND REFERENCE IMAGE
        centrals = []
        for line in range(len(aobj_idJN3)):
            if aobj_idJN3[line][9] == 0 and aobj_idJN3[line][10] == 0:
                centrals.append(line)
        
        if len(centrals) == 0: 
            centrals.append(4) #in case no object has zero offset, align all objects to fifth object
        center = centrals[0]
        
        # READ IN AND INITIALIZE DISTANCE MINIMIZATION
        star_coords = pickle.load(open( "Calibs/starfields/star_coords"+str(obj_id)+'JN3.p', 'rb'))
        
        ref_im = star_coords[center]
        
        def distance(x0,y0,x1,y1):
            return np.sqrt((x0-x1)**2+(y0-y1)**2)
            
        #LOOP THROUGH EACH EXPOSURE IN JN1
        for shf_num in range(len(aobj_idJN3)):
            
            shift_im = np.array(star_coords[shf_num]) 
            
            #----------testing possible xrange to improve time------------#
            xranges = []
            yranges = []
            for sline in range(len(shift_im)):
                for rline in range(len(ref_im)):
                    xranges.append(float(shift_im[sline][0]-ref_im[rline][0]))
                    yranges.append(float(shift_im[sline][1]-ref_im[rline][1]))
                                
            #pdb.set_trace()
            #x_shift = np.array(xranges)
            #y_shift = np.array(yranges)
                    
                
            
            #---------------end testing ---------------------------------------#
            
            #FIRST ITERATION, TO GET ROUGH APPROXIMATION
            x_shift = np.arange(-1000, 1000, 20)
            y_shift = np.arange(-1000, 1000, 20)
            datai = np.zeros((len(x_shift),len(y_shift)))
            
            #Loop through different x and y shifts
            for x in range(len(x_shift)):
               for y in range(len(y_shift)):
                   tot_points = []
                   
                   #Apply shift
                   shift_im = np.array(star_coords[shf_num]) + [x_shift[x],y_shift[y]]
                   dist = dict() #creates a separate variable for each star in shift_im
            
                   #loop through stars in image to shift
                   for shift_line in range(len(shift_im)):
                       dist[shift_line] = [] 
                       
                       #loop through stars in reference image
                       for ref_line in range(len(ref_im)):
                           dist[shift_line].append(distance(ref_im[ref_line][0], ref_im[ref_line][1], shift_im[shift_line][0], shift_im[shift_line][1]))
                           
                       #find distance to closest star in ref_im for each star in shift_im
                       dist[shift_line] = np.min(dist[shift_line])
                       
                       if dist[shift_line] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift[mins[0]]),np.median(y_shift[mins[1]])]
            #print mina
       
            #pdb.set_trace()
            
            #-----------------------------Now Repeat--------------------------------#
            
            #SECOND ITERATION, TO GET BETTER APPROXIMATION       
            x_shift2 = np.arange(float(mina[0])-150, float(mina[0])+150, 10)
            y_shift2 = np.arange(float(mina[1])-150, float(mina[1])+150, 10)
            
            datai = np.zeros((len(x_shift2),len(y_shift2)))

            for x in range(len(x_shift2)):
               for y in range(len(y_shift2)):
                   
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift2[x], y_shift2[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       
                       if dist[shift_star] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift2[mins[0]]),np.median(y_shift2[mins[1]])]
            #print [float(mina[0]),float(mina[1])]

            #THIRD ITERATION, TO GET BETTER APPROXIMATION
            x_shift3 = np.arange(float(mina[0])-20, float(mina[0])+20, 1)
            y_shift3 = np.arange(float(mina[1])-20, float(mina[1])+20, 1)
            datai = np.zeros((len(x_shift3),len(y_shift3)))

            for x in range(len(x_shift3)):
               for y in range(len(y_shift3)):
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift3[x],y_shift3[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       if dist[shift_star] < 4:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift3[mins[0]]),np.median(y_shift3[mins[1]])]
            print mina
            
            #Append optimal x,y shift for each exposure 
            xy_shifts_JN3.append([float(mina[0]),float(mina[1])])
            
        print xy_shifts_JN3
        
        #APPLY SHIFTS TO DATA
        shifteds = []
        for line in range(len(aobj_idJN3)):
            sci = pyfits.getdata(Jobj_idN3_locs[line])
            shifted = interp.shift(sci, xy_shifts_JN3[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = np.mean(shifteds, axis =0) #trying instead of sum
        median = np.median(shifteds, axis = 0)


        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_JN3_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_JN3_median.fits',clobber=True)
        
        #--------------------------------------------------------------------------#
        
    print "reading in KN3"
    
    if len(aobj_idKN3) > 0:
        
        xy_shifts_KN3 = []
                
        # FIND REFERENCE IMAGE
        centrals = []
        for line in range(len(aobj_idKN3)):
            if aobj_idKN3[line][9] == 0 and aobj_idKN3[line][10] == 0:
                centrals.append(line)
        
        if len(centrals) == 0: 
            centrals.append(4) #in case no object has zero offset, align all objects to fifth object
        center = centrals[0]
        
        # READ IN AND INITIALIZE DISTANCE MINIMIZATION
        star_coords = pickle.load(open( "Calibs/starfields/star_coords"+str(obj_id)+'KN3.p', 'rb'))
        
        ref_im = star_coords[center]
        
        def distance(x0,y0,x1,y1):
            return np.sqrt((x0-x1)**2+(y0-y1)**2)
            
        #LOOP THROUGH EACH EXPOSURE IN JN1
        for shf_num in range(len(aobj_idKN3)):
            
            shift_im = np.array(star_coords[shf_num]) 
            
            #----------testing possible xrange to improve time------------#
            xranges = []
            yranges = []
            for sline in range(len(shift_im)):
                for rline in range(len(ref_im)):
                    xranges.append(float(shift_im[sline][0]-ref_im[rline][0]))
                    yranges.append(float(shift_im[sline][1]-ref_im[rline][1]))
                                
            #pdb.set_trace()
            #x_shift = np.array(xranges)
            #y_shift = np.array(yranges)
                    
                
            
            #---------------end testing ---------------------------------------#
            
            #FIRST ITERATION, TO GET ROUGH APPROXIMATION
            x_shift = np.arange(-1000, 1000, 20)
            y_shift = np.arange(-1000, 1000, 20)
            datai = np.zeros((len(x_shift),len(y_shift)))
            
            #Loop through different x and y shifts
            for x in range(len(x_shift)):
               for y in range(len(y_shift)):
                   tot_points = []
                   
                   #Apply shift
                   shift_im = np.array(star_coords[shf_num]) + [x_shift[x],y_shift[y]]
                   dist = dict() #creates a separate variable for each star in shift_im
            
                   #loop through stars in image to shift
                   for shift_line in range(len(shift_im)):
                       dist[shift_line] = [] 
                       
                       #loop through stars in reference image
                       for ref_line in range(len(ref_im)):
                           dist[shift_line].append(distance(ref_im[ref_line][0], ref_im[ref_line][1], shift_im[shift_line][0], shift_im[shift_line][1]))
                           
                       #find distance to closest star in ref_im for each star in shift_im
                       dist[shift_line] = np.min(dist[shift_line])
                       
                       if dist[shift_line] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift[mins[0]]),np.median(y_shift[mins[1]])]
            #print mina
       
            #pdb.set_trace()
            
            #-----------------------------Now Repeat--------------------------------#
            
            #SECOND ITERATION, TO GET BETTER APPROXIMATION       
            x_shift2 = np.arange(float(mina[0])-150, float(mina[0])+150, 10)
            y_shift2 = np.arange(float(mina[1])-150, float(mina[1])+150, 10)
            
            datai = np.zeros((len(x_shift2),len(y_shift2)))

            for x in range(len(x_shift2)):
               for y in range(len(y_shift2)):
                   
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift2[x], y_shift2[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       
                       if dist[shift_star] < 20:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift2[mins[0]]),np.median(y_shift2[mins[1]])]
            #print [float(mina[0]),float(mina[1])]

            #THIRD ITERATION, TO GET BETTER APPROXIMATION
            x_shift3 = np.arange(float(mina[0])-20, float(mina[0])+20, 1)
            y_shift3 = np.arange(float(mina[1])-20, float(mina[1])+20, 1)
            datai = np.zeros((len(x_shift3),len(y_shift3)))

            for x in range(len(x_shift3)):
               for y in range(len(y_shift3)):
                   tot_points = []
                   shift_im = np.array(star_coords[shf_num]) + [x_shift3[x],y_shift3[y]]
                   dist = dict()
                   for shift_star in range(len(shift_im)):
                       dist[shift_star] = [] 
                       for ref_star in range(len(ref_im)):
                           dist[shift_star].append(distance(ref_im[ref_star][0], ref_im[ref_star][1], shift_im[shift_star][0],shift_im[shift_star][1]))
                       dist[shift_star] = np.min(dist[shift_star])
                       if dist[shift_star] < 4:
                           tot_points.append(1.0)
        
                   totpoints = sum(tot_points)      
                   #tot_dist =  sum(dist.values())
                   datai[x, y] = totpoints
                   #shift_im = array(star_coords[0])
            mins = np.where(datai==datai.max())
            
            mina = [np.median(x_shift3[mins[0]]),np.median(y_shift3[mins[1]])]
            print mina
            
            #Append optimal x,y shift for each exposure 
            xy_shifts_KN3.append([float(mina[0]),float(mina[1])])
            
        print xy_shifts_KN3
        
        #APPLY SHIFTS TO DATA
        shifteds = []
        for line in range(len(aobj_idKN3)):
            sci = pyfits.getdata(Kobj_idN3_locs[line])
            shifted = interp.shift(sci, xy_shifts_KN3[line] ,order = 0)
            shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        combined = sum(shifteds)
        median = np.median(shifteds, axis = 0)


        file = pyfits.PrimaryHDU(combined)
        mfile = pyfits.PrimaryHDU(median)
        file.writeto('Calibs/shifted/'+str(obj_id)+'_KN3_combined.fits',clobber=True)
        mfile.writeto('Calibs/shifted/'+str(obj_id)+'_KN3_median.fits',clobber=True)

        #--------------------------------------------------------------------------#
        
        
        
        
        
            
            
            
            
            

        
        
        
        
        