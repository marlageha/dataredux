#Finds the already combined final images from each night and combines them. 

def whirc_combine(obj_id):
    
    import os
    import pyfits
    import numpy as np
    import scipy.ndimage as snd 
    from scipy.ndimage import interpolation as interp
    from astropy.convolution import convolve, Tophat2DKernel, AiryDisk2DKernel, MexicanHat2DKernel, Box2DKernel
    
    '''
    
    try: 
        JN1 = pyfits.getdata('Calibs/shifted/'+str(obj_id)+'_JN1_median.fits')
        
    except IOError:
        print "no JN1..."
        
    try: 
        JN2 = pyfits.getdata('Calibs/shifted/'+str(obj_id)+'_JN2_median.fits')
        
    except IOError:
        print "no JN2..."
        
    try: 
        JN3 = pyfits.getdata('Calibs/shifted/'+str(obj_id)+'_JN3_median.fits')
        
    except IOError:
        print "no JN3..."
        
        
    try: 
        KN1 = pyfits.getdata('Calibs/shifted/'+str(obj_id)+'_KN1_median.fits')
        
    except IOError:
        print "no KN1..."
        
    try: 
        KN2 = pyfits.getdata('Calibs/shifted/'+str(obj_id)+'_KN2_median.fits')
        
    except IOError:
        print "no KN2..."
        
    try: 
        KN3 = pyfits.getdata('Calibs/shifted/'+str(obj_id)+'_KN3_median.fits')
        
    except IOError:
        print "no KN3..."
    '''
    
    if obj_id == 55500:
        num_star = 13
    elif obj_id == 35979:
        num_star = 10
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
    
    files = []
    for file in os.listdir(os.getcwd()+'/Calibs/shifted/'):
        if str(obj_id) in file and "median" in file:
            files.append(file)
    print files
            
    Jfiles = []    
    for line in range(len(files)):
        if "J" in files[line]:
            Jfiles.append('Calibs/shifted/'+files[line])
    print Jfiles 
    
    Kfiles = []    
    for line in range(len(files)):
        if "K" in files[line]:
            Kfiles.append('Calibs/shifted/'+files[line])
    
    if len(Jfiles) > 1:
        print "processing Jfiles..."
    
        tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
        star_coords = []
        centrals = []
        for line in range(len(Jfiles)):
            img = pyfits.getdata(Jfiles[line])
        
            #Smoothing
            smooth = convolve(img, tophat)
            threshold = np.mean(smooth)+3 # bit arbitrary, but it works
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
            
    print "calculationg shifts"
    
    if len(Jfiles) > 1:
        
        xy_shifts_J = []
                
        # FIND REFERENCE IMAGE
        center = 0 #referenc is first image
        
        # READ IN AND INITIALIZE DISTANCE MINIMIZATION
        
        ref_im = star_coords[center]
        
        def distance(x0,y0,x1,y1):
            return np.sqrt((x0-x1)**2+(y0-y1)**2)
            
        #LOOP THROUGH EACH EXPOSURE IN JN1
        for shf_num in range(len(Jfiles)):
            
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
            xy_shifts_J.append([float(mina[0]),float(mina[1])])
            
        print xy_shifts_J
        
        #APPLY SHIFTS TO DATA
        shifteds = []
        for line in range(len(Jfiles)):
            sci = pyfits.getdata(Jfiles[line])
            shifted = interp.shift(sci, xy_shifts_J[line] ,order = 0)
            #shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        median = np.median(shifteds, axis = 0)

        mfile = pyfits.PrimaryHDU(median)
        mfile.writeto('Final/'+str(obj_id)+'_J_median.fits',clobber=True)
        
    if len(Kfiles) > 1:
        print "processing Kfiles..."
    
        tophat = Tophat2DKernel(5) #Radius can be changed to increase sharpess 
    
        star_coords = []
        centrals = []
        for line in range(len(Kfiles)):
            img = pyfits.getdata(Kfiles[line])
        
            #Smoothing
            smooth = convolve(img, tophat)
            threshold = np.mean(smooth)+3 # bit arbitrary, but it works
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
            
    print "calculationg shifts"
    
    if len(Kfiles) > 1:
        
        xy_shifts_K = []
                
        # FIND REFERENCE IMAGE
        center = 0 #referenc is first image
        
        # READ IN AND INITIALIZE DISTANCE MINIMIZATION
        
        ref_im = star_coords[center]
        
        def distance(x0,y0,x1,y1):
            return np.sqrt((x0-x1)**2+(y0-y1)**2)
            
        #LOOP THROUGH EACH EXPOSURE IN JN1
        for shf_num in range(len(Kfiles)):
            
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
            xy_shifts_K.append([float(mina[0]),float(mina[1])])
            
        print xy_shifts_K
        
        #APPLY SHIFTS TO DATA
        shifteds = []
        for line in range(len(Kfiles)):
            sci = pyfits.getdata(Kfiles[line])
            shifted = interp.shift(sci, xy_shifts_K[line] ,order = 0)
            #shifted = shifted - np.median(shifted)
            shifteds.append(shifted)
        median = np.median(shifteds, axis = 0)

        mfile = pyfits.PrimaryHDU(median)
        mfile.writeto('Final/'+str(obj_id)+'_K_median.fits',clobber=True)
    
    