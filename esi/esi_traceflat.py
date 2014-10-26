##################################################
#      Gets orders from flat                     #
#      Creates normalized flat                   #
#      Writes orders masks to file               #
#                                                #
#      Kareem El-Badry, 07/07/2014               #
##################################################

#so we can use a function from another program in this module.
from esi import esi_trace_cen
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
import pdb
import pickle
import pyfits

    
class orders:
    def __init__(self, xsize, ysize):  
       self.xl = np.zeros(xsize) #left edge, x coords
       self.xr = np.zeros(xsize) #right edge, x coords
       self.yl = np.zeros(ysize)
       self.yr = np.zeros(ysize)

def esi_traceflat():
    
    # READ FLAT FRAME
    ffile = pyfits.open('Calibs/dome_flat.fits')
    flat = ffile[0].data

    # SMOOTH FLAT SLIGHTLY
    sflat = ndimage.gaussian_filter(flat, 3)

    lx = len(sflat[0,:])
    ly = len(sflat[:,0])

    # FIND CENTRAL PEAKS
    print "finding central gradient peaks..."
    xcen, ycen = esi_trace_cen()

    # CREATE ARRAY FOR ESI ORDERS -- REDDEST ORDER IS FIRST
    xtrace=np.zeros([len(xcen),4096])
    ytrace=np.zeros([len(ycen),4096])

    xtrace[:,2000] = xcen
    ytrace[:,2000] = ycen


    #  TRACE IN FOUR PARTS, STARTING FROM CENTRAL ROW
    #  (1)  WORK DOWN FOR FULL ARRAY
    #  (2)  WORK UP FOR FULL ARRAY
    '''
    #  (3)  FINISH WORKING UP, REMOVE BLUE-MOST ORDER
    #  (4)  FINISH WORKING DOWN, REMOVE BLUE-MOST ORDER
    '''

    # (1) WORK DOWN FOR FULL ARRAY
    print "extending peaks..."
    for row in xrange(1999,0,-1):
    	# SUM +/- 5 PIXELS TO DETERMINE PEAK FOR CURRENT ROW
        tflat = sflat[row-5:row+5,:]


        # FIND PEAKS FOR CURRENT ROW
        dx = ndimage.sobel(tflat, 0)  # horizontal derivative
        dy = ndimage.sobel(tflat, 1)  # vertical derivative
        mag = np.hypot(dx, dy) 		  # magnitude
        mag1d = np.sum(mag,axis=0)
        x = np.arange(0,len(mag1d))

        for n in xrange(0,20,1):        
            xt = np.int(np.round(xtrace[n,row+1]))
            xarr = np.arange(xt-5,xt+5)
            local_peak = mag1d[xarr]
            mx = local_peak.argmax() #highest point in +/- 5 pixels

            xtrace[n,row] = mx-5+xt
            ytrace[n,row] = row
            #pdb.set_trace()


    # (2)  WORK UP FOR FULL ARRAY
    for row in range(2001,4096):

    	# SUM +/- 5 PIXELS TO DETERMINE PEAK FOR CURRENT ROW
        tflat = sflat[row-5:row+5,:]


    	# FIND PEAKS FOR CURRENT ROW
        dx = ndimage.sobel(tflat, 0)  # horizontal derivative
        dy = ndimage.sobel(tflat, 1)  # vertical derivative
        mag = np.hypot(dx, dy) 		  # magnitude
        mag1d = np.sum(mag,axis=0)
        x = np.arange(0,len(mag1d))

        for n in xrange(0,20,1):
            xt = np.int(np.round(xtrace[n,row-1]))
            xarr = np.arange(xt-5,xt+5)
            local_peak = mag1d[xarr]
    
            #just trying alternative:
            mx = local_peak.argmax() #highest point instead of local max finding.
    		#mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1

            xtrace[n,row] = mx-5+xt
            ytrace[n,row] = row

    '''

    # (3)  FINISH WORKING UP, REMOVE BLUE-MOST ORDER
    for row in range(3550,4005):

    	# SUM +/- 5 PIXELS TO DETERMINE PEAK FOR CURRENT ROW
        tflat = sflat[row-5:row+5,:]



        # FIND PEAKS FOR CURRENT ROW
        dx = ndimage.sobel(tflat, 0)  # horizontal derivative
        dy = ndimage.sobel(tflat, 1)  # vertical derivative
        mag = np.hypot(dx, dy) 		  # magnitude
        mag1d = np.sum(mag,axis=0)
        x = np.arange(0,len(mag1d))

        for n in xrange(2,16,1):
            xt = np.int(np.round(xtrace[n,row-1]))#
            xarr = np.arange(xt-5,xt+5)
            local_peak = mag1d[xarr]#
            #just trying alternative:
            mx = local_peak.argmax() #highest point instead of local max finding.
    
            #mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1
            #print row, n ,mx
            xtrace[n,row] = mx-5+xt
            ytrace[n,row] = row

        if (row < 3700):
            for n in xrange(0,2,1):
                xt = np.int(np.round(xtrace[n,row-1]))#
                xarr = np.arange(xt-5,xt+5)
                local_peak = mag1d[xarr]#
                #just trying alternative:
                mx = local_peak.argmax() #highest point instead of local max finding.
        
                #mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1
                #print row, n ,mx

                xtrace[n,row] = mx-5+xt
                ytrace[n,row] = row


    # (4) FINISH WORKING DOWN FOR FULL ARRAY
    for row in xrange(920,50,-1):
    	# SUM +/- 5 PIXELS TO DETERMINE PEAK FOR CURRENT ROW
        tflat = sflat[row-5:row+5,:]


    	# FIND PEAKS FOR CURRENT ROW
        dx = ndimage.sobel(tflat, 0)  # horizontal derivative
        dy = ndimage.sobel(tflat, 1)  # vertical derivative
        mag = np.hypot(dx, dy) 		  # magnitude
        mag1d = np.sum(mag,axis=0)
        x = np.arange(0,len(mag1d))


        for n in xrange(0,14,1):
            xt = np.int(np.round(xtrace[n,row+1]))
            xarr = np.arange(xt-5,xt+5)
            local_peak = mag1d[xarr]
    
            #just trying alternative:
            mx = local_peak.argmax() #highest point instead of local max finding.
    
    		#mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1
    #		plt.plot(xarr,local_peak,'r+')#
    #		plt.show()

            xtrace[n,row] = mx-5+xt
            ytrace[n,row] = row
        if (row > 550):
            for n in xrange(14,16,1):
                xt = np.int(np.round(xtrace[n,row+1]))
                xarr = np.arange(xt-5,xt+5)
                local_peak = mag1d[xarr]
        
                #just trying alternative:
                mx = local_peak.argmax() #highest point instead of local max finding.
    			#mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1

                xtrace[n,row] = mx-5+xt
                ytrace[n,row] = row


    #####################################################
    print "plotting and polyfitting..."
    f = plt.figure()
    ax = f.add_subplot(111)

    #sflat = pyfits.getdata('Raw/e140306_0046.fits.gz')
    #sflat = sflat[:,25:2070]
    #plt.imshow(sflat,cmap='gray', vmax = np.median(sflat)+0.1*np.std(sflat))

    plt.imshow(sflat, cmap = 'gray')
    for n in xrange(0,20,1): #loop through all 20 edges
        #plt.plot(xtrace[n],ytrace[n], '.', color = 'blue', markersize = 1) #actual points
        #returns array [x0, x1, x2, x3]
        p = np.polyfit(ytrace[n,1000:3400],xtrace[n,1000:3400],3) 
        #print p
        xp = np.poly1d(p)
        y=np.arange(0,4096,1)	
        #print xp(y)
        plt.plot(xp(y),y, 'b')
    #plt.plot(xcen, ycen, 'o', color = 'green')
    plt.title('orders')
    plt.xlabel('x pixels')
    plt.ylabel('y pixels')
    ax.set_xlim([0, 2044])
    ax.set_ylim([0, 4096])
    plt.savefig('Calibs/pics/orders2.png')
    '''

    #SAVE ORDERS ACCORDING TO POLYFIT
    edge_fits = [np.polyfit(ytrace[n,1000:3400], xtrace[n,1000:3400],3) for n in range(20)]

    polyx = []
    polyy = []
    y=np.arange(0,4096,1)	

    for line in range(20)[::-1]: #go from left to right
        polyx.append(np.poly1d(edge_fits[line]))# add a (y) if you want an array
        polyy.append(np.arange(0,4096,1))

    esiorders = []
    for line in range(len(edge_fits)/2):
        esiorders.append(orders(2045, 4096))
        esiorders[line].xl = polyx[2*line]
        esiorders[line].xr = polyx[2*line + 1]
        esiorders[line].yl = polyy[2*line]
        esiorders[line].yr = polyy[2*line + 1]
    '''
    for order in range(len(esiorders)): 
        num = order

        order_left = esiorders[num].xl
        order_right = esiorders[num].xr
        width = order_right - order_left

        lx = len(sflat[0,:])
        ly = len(sflat[:,0])
        X, Y = np.ogrid[0:ly, 0:lx]

        mask_left = order_left(X) < Y 
        mask_right = order_right(X) > Y

        #only accept what passes both masks:
        mask = mask_left*mask_right
        z = np.zeros((ly, lx))
        z[mask] = flat[mask]

        ord_tot = []
        for line in range(ly):
            line_sum = np.mean(z[line])
            ord_tot.append(line_sum)
    
        #p = np.polyfit(y, ord_tot, 5)
        #xp = poly1d(p)

        smoothed = gaussian_filter1d(ord_tot, 100)

        plt.plot(y, ord_tot, '.', markersize = 1, color = 'black')
        #plt.plot(y, xp(y), color = 'magenta')
        plt.plot(y, smoothed, color = 'magenta')
    plt.xlabel('y - pixel')
    plt.ylabel('intensity')
    plt.title('Gaussian Smoothing with $\sigma = 100$')
    plt.savefig('Calibs/pics/flat_int_ordspline.png')
    '''

    #MAKE MASK FOR ALL ORDERS
    #master_mask == all the orders; background_mask == background (doh!)
    print 'writing masks...'
    #edges 
    pickle.dump(esiorders, open('Calibs/order_edges.p', 'wb'))
    
    
    master_mask = np.zeros((ly, lx), dtype = bool)
    all_order_masks = []
    sky_mask = []
    for order in range(len(esiorders)):
    
        num = order

        order_left = esiorders[num].xl
        order_right = esiorders[num].xr
    
        lx = len(sflat[0,:])
        ly = len(sflat[:,0])
        X, Y = np.ogrid[0:ly, 0:lx]

        mask_left = order_left(X) < Y
        mask_left3 = order_left(X) +13 < Y #because there's less light at the edges
        mask_left15 = order_left(X)+25 > Y #15 pixels to the right of order edge 
        mask_right = order_right(X) > Y
        mask_right15 = order_right(X)-25 < Y
        mask_right3 = order_right(X)-3 > Y

        #only accept what passes both masks:
        mask = mask_left*mask_right #for one mask
        
        smask_left = mask_left3*mask_left15
        smask_right = mask_right3*mask_right15
        
        smask = smask_left + smask_right
        
        all_order_masks.append(mask)
        sky_mask.append(smask)
                
        #add all orders
        master_mask = master_mask + mask 
    
    background_mask = -master_mask
    
    #make rough bad pixel map
    some_flat = pyfits.getdata('Raw/e140306_0006.fits.gz')[:, 25:2070]
    
    too_hot = some_flat > 2500
    mask1 = np.zeros((4096, 2045), dtype = 'bool')
    mask1[2600:4096, 410:450] = True
    bad_pix = too_hot*mask1
    
    #write masks to file
    pickle.dump(master_mask, open('Calibs/orders_mask.p', 'wb'))
    pickle.dump(background_mask, open('Calibs/background_mask.p', 'wb'))
    pickle.dump(all_order_masks, open('Calibs/all_order_masks.p', 'wb'))
    pickle.dump(sky_mask, open('Calibs/sky_mask.p', 'wb'))
    pickle.dump(bad_pix, open('Calibs/bad_pix.p', 'wb'))

    # DIVIDE EACH ORDER BY SMOOTHED PROFILE

    #f = plt.figure()
    #ax = f.add_subplot(111)
    #ax.set_ylim([-1, 2])
    print "normalizing flat..."
    z = np.zeros((ly, lx)) #to hold divided orders
    order_polys = []
    for order in range(len(esiorders)):

        num = order #looping through orders

        order_left = esiorders[num].xl
        order_right = esiorders[num].xr
        width = order_right - order_left

        #making a maks
        lx = len(sflat[0,:])
        ly = len(sflat[:,0])
        X, Y = np.ogrid[0:ly, 0:lx]

        mask_left = order_left(X) < Y 
        mask_right = order_right(X) > Y

        #only accept what passes both masks:
        mask = mask_left*mask_right

        #dummy order to hold each order while it is divided
        blank = np.zeros((ly, lx))
        blank[mask] = flat[mask]

        #making the profile for each order
        ord_tot = []
        for line in range(ly):
            line_mean = np.mean(blank[line])
            ord_tot.append(line_mean)

        #fit 10th order polynomial
        p = np.polyfit(y, ord_tot, 10) 
        xp = np.poly1d(p)
        order_polys.append(xp)

        #dividing each order line by line
        for line in range(ly):
           blank[line] = blank[line]/xp(line)

        z[mask] = blank[mask]
    pickle.dump(order_polys, open('Calibs/order_polys.p', 'wb'))
    z[background_mask] = 1.0 
    
    #Get rid of zeros in flat
    x_zero = np.where(z == 0)[0]
    y_zero = np.where(z == 0)[1]

    for line in range(len(x_zero)):
        adjacent = [(x_zero[line], y_zero[line]+1), (x_zero[line], y_zero[line]-1), 
                    (x_zero[line]+1, y_zero[line]), (x_zero[line]-1, y_zero[line])]
                
        try:
            adj_vals = [z[x] for x in adjacent]
        
        except IndexError:
             #If it's at the edge and there's no adjacent:
            adj_vals = [12.0 for x in range(4)] #Like whatevz. 
        
    
        z[x_zero[line], y_zero[line]] = np.mean(adj_vals)
    
    '''
    f = plt.figure()
    plt.imshow(z[1400:2800, :], cmap = 'gray')
    plt.savefig('Calibs/pics/normflat.tiff')
    '''

    #Divide old flat by new flat, write to file
    new_flat = np.zeros((ly, lx))
    for order in range(len(esiorders)):

        num = order

        order_left = esiorders[num].xl
        order_right = esiorders[num].xr

        lx = len(sflat[0,:])
        ly = len(sflat[:,0])
        X, Y = np.ogrid[0:ly, 0:lx]

        mask_left = order_left(X) < Y 
        mask_right = order_right(X) > Y

        #only accept what passes both masks:
        mask = mask_left*mask_right

        #dummy = np.zeros((ly, lx))
        new_flat[mask] = flat[mask]/z[mask]
    
    
    '''
    test = np.zeros((ly, lx))
    test[master_mask] = 1
    test[background_mask] = 0.5
    q = plt.figure()
    plt.imshow(test, cmap = 'gray')
    # turns out you don't have to do this, you can just do 
    # plt.imshow(master_mask)


    g = plt.figure()
    plt.imshow(new_flat[1400:2800, :], cmap = 'gray')
    plt.savefig('Calibs/pics/dividedflat.tiff')
    '''

    divided = pyfits.PrimaryHDU(new_flat)
    divider = pyfits.PrimaryHDU(z)

    divided.writeto('Calibs/divided_flat.fits', clobber = True)
    divider.writeto('Calibs/norm_flat.fits', clobber = True)

