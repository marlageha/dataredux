#Displays galaxies and fits nicely

def whirc_draw(obj_id):
    
    import pyfits
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.ndimage as snd 
    
    
    #FIND RELEVANT IMAGEBLOCK (FROM GALFIT)
    
    #J File
    try: 
        img = pyfits.open('Calibs/galfits/'+str(obj_id)+'_J_block.fits')
        
        galaxy = img[1].data
        fit = img[2].data
        resid = img[3].data
        
        width = np.shape(galaxy)[0]
        height = np.shape(galaxy)[1]

        lima = 0.98*np.max(galaxy[0.25*width : 0.75*width, 0.25*height : 0.75*height])
        limb = 0.98*np.max(resid[0.25*width : 0.75*width, 0.25*height : 0.75*height])

        plt.figure()
    

        plt.subplot(131)
        plt.imshow(galaxy, cmap=plt.cm.gray, vmin=0, vmax=lima)
        plt.axis('off')

        plt.subplot(132)
        plt.imshow(fit, cmap=plt.cm.gray)
        plt.title(str(obj_id)+ ' J Filter')
        plt.axis('off')

        plt.subplot(133)
        plt.imshow(resid, cmap=plt.cm.gray, vmin=0, vmax=limb)
        plt.axis('off')
    

        plt.subplots_adjust(wspace=0, hspace=0., top=0.99, bottom=0.01, left=0.05, right=0.99)
        plt.savefig('Calibs/pics/'+str(obj_id)+'_J.png', clobber = True)
        '''
        #INTENSITY PLOT
        def Sersic(r):
            return brightness*np.exp(-bn*((r/r_e)**(1/ind)-1))
            
        def Sersic2(r):
            return brightness2*np.exp(-bn2*((r/r_e2)**(1/ind2)-1))
    
        
        #model = img[2].data
        hdr = img[2].header

        #FIRST COMPONENT
        profile = hdr['COMP_1']

        brightness = hdr['1_MAG']
        p = brightness.split()
        brightness = float(p[0])

        #effective radius
        r_e = float(hdr['1_RE'].split()[0])

        #sersic index
        ind = float(hdr['1_N'].split()[0])

        bn = 2*ind - 0.327 #from Ciotti 1991, good approximation w/o solving Gamma functions 

        #axis ratio (b/a)
        ar = float(hdr['1_AR'].split()[0])
         
        brlist = []   
        try:
            profile2 = hdr['COMP_2']
            brightness2 = float(hdr['2_MAG'].split()[0])
            r_e2 = float(hdr['2_RE'].split()[0])
            ind2 = float(hdr['2_N'].split()[0])
            bn2 = 2*ind2 - 0.327 #from Ciotti 1991, good approximation w/o solving Gamma functions 
            ar2 = float(hdr['2_AR'].split()[0])
            
            brlist.append(brightness2)
        except KeyError:
            print 'just one component...'
  
        
        if len(brlist) > 0:
            if brightness < brightness2:
                xr = np.arange(0, 2*r_e, 0.1)
                
                f = plt.figure()
                
                ax = f.add_subplot(111)
        
                plt.semilogy(xr, Sersic(xr),'b', label = '1st component')
                plt.title(str(obj_id)+' J filter')
                plt.xlabel('Radius (pix)')
                plt.ylabel('Intensity')
                plt.text(0.7, 0.9,'{$n$ = '+str(ind)+', $R_{e} =$ '+str(r_e)+'}', ha='center', va='center', transform=ax.transAxes, color = 'blue')
                plt.semilogy(xr, Sersic2(xr),'r', label = '2nd component')
                plt.semilogy(xr, Sersic(xr) + Sersic2(xr),'k', label = 'total')
                plt.text(0.7, 0.8,'{$n$ = '+str(ind2)+', $R_{e} =$ '+str(r_e2)+'}', ha='center', va='center', transform=ax.transAxes, color = 'red')
                ax.legend(loc = 3)
          
            else:
                xr = np.arange(0, 2*r_e2, 0.1)
                f = plt.figure()
                
                ax = f.add_subplot(111)
        
                plt.semilogy(xr, Sersic2(xr),'b', label = '1st component')
                plt.title(str(obj_id)+' J filter')
                plt.xlabel('Radius (pix)')
                plt.ylabel('Intensity')
                plt.text(0.7, 0.9,'{$n$ = '+str(ind2)+', $R_{e} =$ '+str(r_e2)+'}', ha='center', va='center', transform=ax.transAxes, color = 'blue')
                plt.semilogy(xr, Sersic(xr),'r', label = '2nd component')
                plt.semilogy(xr, Sersic(xr) + Sersic2(xr),'k', label = 'total')
                plt.text(0.7, 0.85,'{$n$ = '+str(ind)+', $R_{e} =$ '+str(r_e)+'}', ha='center', va='center', transform=ax.transAxes, color = 'red')
                ax.legend(loc = 3)
 
        else: 
            xr = np.arange(0,2*r_e,0.1)
            f = plt.figure()
            ax = f.add_subplot(111)
    
            plt.semilogy(xr, Sersic(xr),'k', label = '1st component')
            plt.title(str(obj_id)+' J filter')
            plt.xlabel('Radius (pix)')
            plt.ylabel('Intensity')
            plt.text(0.7, 0.9,'{$n$ = '+str(ind)+', $R_{e} =$ '+str(r_e)+'}', ha='center', va='center', transform=ax.transAxes, color = 'blue')
            plt.semilogy(xr, Sersic2(xr),'r', label = '2nd component')
            plt.semilogy(xr, Sersic(xr) + Sersic2(xr),'k', label = 'total')
            plt.text(0.7, 0.8,'{$n$ = '+str(ind2)+', $R_{e} =$ '+str(r_e2)+'}', ha='center', va='center', transform=ax.transAxes, color = 'red')
            ax.legend(loc = 3)
        
        
     
        
        plt.savefig('Calibs/pics/'+str(obj_id)+'_J_plot.png', clobber = True) '''
        
        
    except IOError:
        print "no J File..."
        
    #KN1
    try: 
        img = pyfits.open('Calibs/galfits/'+str(obj_id)+'_K_block.fits')
        
        galaxy = img[1].data
        fit = img[2].data
        resid = img[3].data
        
        width = np.shape(galaxy)[0]
        height = np.shape(galaxy)[1]

        lima = 0.98*np.max(galaxy[0.25*width : 0.75*width, 0.25*height : 0.75*height])
        limb = 0.98*np.max(resid[0.25*width : 0.75*width, 0.25*height : 0.75*height])

        plt.figure()
    

        plt.subplot(131)
        plt.imshow(galaxy, cmap=plt.cm.gray, vmin=0, vmax=lima)
        plt.axis('off')

        plt.subplot(132)
        plt.imshow(fit, cmap=plt.cm.gray)
        plt.title(str(obj_id) +' K Filter')
        plt.axis('off')

        plt.subplot(133)
        plt.imshow(resid, cmap=plt.cm.gray, vmin=0, vmax=limb)
        plt.axis('off')
    

        plt.subplots_adjust(wspace=0, hspace=0., top=0.99, bottom=0.01, left=0.05, right=0.99)
        plt.savefig('Calibs/pics/'+str(obj_id)+'_K.png', clobber = True)
        '''
        #INTENSITY PLOT
        def Sersic(r):
            return brightness*np.exp(-bn*((r/r_e)**(1/ind)-1))
    
        #model = img[2].data
        hdr = img[2].header

        #FIRST COMPONENT
        profile = hdr['COMP_1']

        brightness = hdr['1_MAG']
        p = brightness.split()
        brightness = float(p[0])

        #effective radius
        r_e = float(hdr['1_RE'].split()[0])

        #sersic index
        ind = float(hdr['1_N'].split()[0])

        bn = 2*ind - 0.327 #from Ciotti 1991, good approximation w/o solving Gamma functions 

        #axis ratio (b/a)
        ar = float(hdr['1_AR'].split()[0])
        
        xr = np.arange(0,2*r_e,0.1)
        
        f = plt.figure()
        ax = f.add_subplot(111)
        
        plt.semilogy(xr, Sersic(xr),'k')
        plt.title(str(obj_id)+' K Filter')
        plt.xlabel('Radius (pix)')
        plt.ylabel('Intensity')
        plt.text(0.7, 0.9,'{$n$ = '+str(ind)+', $R_{e} =$ '+str(r_e)+'}', ha='center', va='center', transform=ax.transAxes)
        
        try:
            profile2 = hdr['COMP_2']
            brightness2 = float(hdr['2_MAG'].split()[0])
            r_e2 = float(hdr['2_RE'].split()[0])
            ind2 = float(hdr['2_N'].split()[0])
            bn2 = 2*ind2 - 0.327 #from Ciotti 1991, good approximation w/o solving Gamma functions 
            ar2 = float(hdr['2_AR'].split()[0])
            
            plt.semilogy(xr, Sersic2(xr), 'b')
            
        except KeyError:
            print 'just one component...'
        
        plt.savefig('Calibs/pics/'+str(obj_id)+'_K_plot.png', clobber = True)
        '''
        
    except IOError:
        print "no K File..."
        
        
    #KN1
    try: 
        img = pyfits.open('Calibs/galfits/'+str(obj_id)+'_r_block.fits')
        
        galaxy = img[1].data
        fit = img[2].data
        resid = img[3].data
        
        width = np.shape(galaxy)[0]
        height = np.shape(galaxy)[1]

        lima = 0.98*np.max(galaxy[0.25*width : 0.75*width, 0.25*height : 0.75*height])
        limb = 0.98*np.max(resid[0.25*width : 0.75*width, 0.25*height : 0.75*height])

        plt.figure()

        plt.subplot(131)
        plt.imshow(galaxy, cmap=plt.cm.gray, vmin=0, vmax=lima)
        plt.axis('off')

        plt.subplot(132)
        plt.imshow(fit, cmap=plt.cm.gray)
        plt.title(str(obj_id) +' SSDS r-band')
        plt.axis('off')

        plt.subplot(133)
        plt.imshow(resid, cmap=plt.cm.gray, vmin=0, vmax=limb)
        plt.axis('off')
    

        plt.subplots_adjust(wspace=0, hspace=0., top=0.99, bottom=0.01, left=0.05, right=0.99)
        plt.savefig('Calibs/pics/'+str(obj_id)+'_r.png', clobber = True)
    
    except IOError:
        print 'no SDSS file'
    
    
        
    
    
        
    
    
    