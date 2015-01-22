##################################################
#  ESI_trace_cen
#
#  USING FLAT FIELD, FIND X/Y FOR ORDER EDGES  
#
#  Kareem El-Badry, 07/25/2014
##################################################


import pyfits
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage
import pdb

def esi_trace_cen(date):

    # READ FLAT FRAME
    ffile = pyfits.open(str(date)+'/Calibs/dome_flat_'+str(date)+'.fits')
    flat = ffile[0].data

    # SMOOTH FLAT SLIGHTLY
    sflat = ndimage.gaussian_filter(flat, 3)


    # FIND EDGES @ CENTER OF CHIP
    tflat = sflat[1995:2005,:]

    dx = ndimage.sobel(tflat, 0)  # horizontal derivative
    dy = ndimage.sobel(tflat, 1)  # vertical derivative
    mag = np.hypot(dx, dy)  # magnitude of derivative
    mag1d = np.sum(mag, axis=0)
    
    '''
    plt.plot(range(len(mag1d)), mag1d)
    plt.title('One Dimensional Gradient')
    plt.xlabel('X pixel')
    plt.ylabel('gradient')
    plt.savefig('esi_1d_mag.png')    '''
    
    x = np.arange(0,len(mag1d)) #[0,1,2,3,4...2048]
    y = np.ones(len(mag1d))*2000. #[2000,2000,2000....2000]

    #This line returns an array which gives the (element + 1) of mag1d where:
    # the difference of the sign of the difference between adjacent element
    # of mag1d (always either 0, -2, or 2, becuase sign is +/- 1) is -2. 
    # In other words, the (element+1) where mag1d goes from increasing to 
    # decreasasing. This corresponds to a local max. 
    mx = (np.diff(np.sign(np.diff(mag1d))) < 0).nonzero()[0] + 1 # local max

    mg = mag1d[mx] # the heigh at the maxima 


    # THIS GETS RED PEAKS
    pred  = np.max(mag1d) * 0.1 #threshold all red peaks are above 
    xcenr = mx[mag1d[mx] > pred] #finds locations of red peaks
    ycenr = y[mag1d[mx]  > pred] # all 2000s..
    mcenr = mg[mag1d[mx] > pred] # heigh of gradient at red peaks

    # THIS GETS BLUE PEAKS
    pblue = 2
    xmax  = np.min(xcenr)

    # & gives intersection of two sets
    xcenb = mx[(mag1d[mx] > pblue)&(mx < xmax)] 
    ycenb = y[(mag1d[mx]  > pblue) & (mx < xmax)]
    mcenb = mg[(mag1d[mx]  > pblue) & (mx < xmax)]


    # PUT THEM TOGETHER, SHOULD BE 20 PEAKS
    xcen = np.concatenate([xcenr,xcenb])
    ycen = np.concatenate([ycenr,ycenb])
    mcen = np.concatenate([mcenr,mcenb])

    # PLOT EDGES
    '''
    plt.close()
    plt.plot(x,mag1d, label = 'gradient')
    plt.plot(x[mx], mag1d[mx], "go", label="erroneous maxima")
    plt.plot(xcen,mcen, "ro", label = 'edges of orders')
    plt.xlabel('x pixel')
    plt.ylabel('gradient')
    plt.legend(loc = 0)
    plt.savefig('Calibs/pics/cen_tes_max.png')
    '''
    xcen = xcen[np.argsort(xcen)][::-1]    # reverse sort


    return xcen, ycen





