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

def esi_traceflat(date):
    
    #READ LOG
    im1 = open(str(date)+'/Logs/esi_info_'+str(date)+'.dat','r')
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
    
    # READ FLAT FRAME
    ffile = pyfits.open(str(date)+'/Calibs/dome_flat_'+str(date)+'.fits')
    flat = ffile[0].data

    # SMOOTH FLAT SLIGHTLY
    sflat = ndimage.gaussian_filter(flat, 3)

    lx = len(sflat[0,:])
    ly = len(sflat[:,0])

    # FIND CENTRAL PEAKS
    print "finding central gradient peaks..."
    xcen, ycen = esi_trace_cen(date)

    # CREATE ARRAY FOR ESI ORDERS -- REDDEST ORDER IS FIRST
    xtrace=np.zeros([len(xcen),4096])
    ytrace=np.zeros([len(ycen),4096])

    xtrace[:,2000] = xcen
    ytrace[:,2000] = ycen


    #  TRACE IN TWO PARTS, STARTING FROM CENTRAL ROW
    #  (1)  WORK DOWN FOR FULL ARRAY
    #  (2)  WORK UP FOR FULL ARRAY


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


    #MAKE MASK FOR ALL ORDERS (and other masks needed for reduction)
    print 'writing masks...'
    
    #Edges, to make other masks easily later. Take less space
    pickle.dump(esiorders, open(str(date)+'/Calibs/order_edges_'+str(date)+'.p', 'wb'))
    
    # Will hold all orders in one mask. Opposite of background.
    master_mask = np.zeros((ly, lx), dtype = bool)
    
    # Will hold all orders in separate masks in an array. Lets you choose one mask
    all_order_masks = []
    
    # For sky subraction
    sky_mask = [] 
    
    # Rejects pixels near edge of order. Can change exact definition. 
    cen_mask = []
    
    for num in range(len(esiorders)):
    
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
        mask_right3 = order_right(X)-13 > Y
        cen_left = order_left(X) + 30 < Y
        cen_right = order_right(X) - 30 > Y

        cen_mas = cen_left*cen_right #inner ~60 pixels or so
        
        #only accept what passes both masks:
        mask = mask_left*mask_right #for one order mask
        
        #sky mask
        smask_left = mask_left3*mask_left15
        smask_right = mask_right3*mask_right15
        
        #add together right and left sides of sky
        smask = smask_left + smask_right
        
        #build order by order
        all_order_masks.append(mask)
        sky_mask.append(smask)
        cen_mask.append(cen_mas)
                
        #add all orders
        master_mask = master_mask + mask 
    
    background_mask = -master_mask
    
    #make a bad pixel map. First find a random domeflat
    dflat = []
    for line in range(len(good)):
        if "DmFlat" in good[line][3] and "Dome" in good[line][2]:
            dflat.append(good[line])
    some_path = str(date)+'/Raw/'+str(dflat[0][0]) #pick the first flat. 
    some_flat = pyfits.getdata(some_path)[:, 25:2070]
    
    too_hot = some_flat > 2500
    mask1 = np.zeros((4096, 2045), dtype = 'bool')
    mask1[2600:4096, 410:450] = True
    bad_pix = too_hot*mask1
    
    #write masks to file
    pickle.dump(master_mask, open(str(date)+'/Calibs/orders_mask_'+str(date)+'.p', 'wb'))
    pickle.dump(background_mask, open(str(date)+'/Calibs/background_mask_'+str(date)+'.p', 'wb'))
    pickle.dump(all_order_masks, open(str(date)+'/Calibs/all_order_masks_'+str(date)+'.p', 'wb'))
    pickle.dump(sky_mask, open(str(date)+'/Calibs/sky_mask_'+str(date)+'.p', 'wb'))
    pickle.dump(bad_pix, open(str(date)+'/Calibs/bad_pix_'+str(date)+'.p', 'wb'))
    pickle.dump(cen_mask, open(str(date)+'/Calibs/cen_mask_'+str(date)+'.p', 'wb'))
    

    print "normalizing flat..."
    
    # DIVIDE EACH ORDER BY SMOOTHED PROFILE
    z = np.zeros((ly, lx)) #to hold divided orders
    order_polys = []
    
    for num in range(len(esiorders)):

        order_left = esiorders[num].xl
        order_right = esiorders[num].xr
        width = order_right - order_left

        #making a mask
        lx = len(sflat[0,:])
        ly = len(sflat[:,0])
        X, Y = np.ogrid[0:ly, 0:lx]

        mask_left = order_left(X) < Y 
        mask_right = order_right(X) > Y

        #only accept what passes both masks:
        mask = mask_left*mask_right

        # dummy order to hold each order while it is divided
        # can probably divide it straight away but this is more straight forward
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

        #dividing each row of pixels by corresponding polynomial
        for line in range(ly):
           blank[line] = blank[line]/xp(line)

        z[mask] = blank[mask]
        
    pickle.dump(order_polys, open(str(date)+'/Calibs/order_polys_'+str(date)+'.p', 'wb'))
    z[background_mask] = 1.0 
    
    #Get rid of zeros in flat
    x_zero = np.where(z == 0)[0]
    y_zero = np.where(z == 0)[1]

    # Valiant, but probably unnecessary effort to set to value of nearby pixels
    for line in range(len(x_zero)):
        adjacent = [(x_zero[line], y_zero[line]+1), (x_zero[line], y_zero[line]-1), 
                    (x_zero[line]+1, y_zero[line]), (x_zero[line]-1, y_zero[line])]
                
        try:
            adj_vals = [z[x] for x in adjacent]
        
        except IndexError:
             #If it's at the edge and there's no adjacent:
            adj_vals = [12.0 for x in range(4)] #setting to 12. These edges don't matter. 
        
    
        z[x_zero[line], y_zero[line]] = np.mean(adj_vals)
    

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
    
    
    divided = pyfits.PrimaryHDU(new_flat)
    divider = pyfits.PrimaryHDU(z)

    divided.writeto(str(date)+'/Calibs/divided_flat_'+str(date)+'.fits', clobber = True)
    divider.writeto(str(date)+'/Calibs/norm_flat_'+str(date)+'.fits', clobber = True)

