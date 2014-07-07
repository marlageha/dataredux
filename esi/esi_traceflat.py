#! /usr/bin/python
##################################################
#  ESI_traceflat
#
#  USING FLAT FIELD, FIND X/Y FOR ORDER EDGES  
#
#  MG 4/14
##################################################

#so we can use a function from another program in this module.
from esi import esi_trace_cen
import pyfits
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
import pdb

    
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
    '''

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
   
    #SAVE ORDERS ACCORDING TO POLYFIT
    edge_fits = [np.polyfit(ytrace[n,1000:3400], xtrace[n,1000:3400],3) for n in range(20)]

    polyx = []
    polyy = []

    for line in range(20)[::-1]: #go from left to right
        polyx.append(np.poly1d(edge_fits[line])(y))
        polyy.append(np.arange(0,4096,1))

    esiorders = []
    for line in range(len(edge_fits)/2):
        esiorders.append(orders(2044, 4096))
        esiorders[line].xl = polyx[2*line]
        esiorders[line].xr = polyx[2*line + 1]
        esiorders[line].yl = polyy[2*line]
        esiorders[line].yr = polyy[2*line + 1]