#! /usr/bin/python
##################################################
#  ESI_trace_cen
#
#  USING FLAT FIELD, FIND X/Y FOR ORDER EDGES  
#
#  MG 4/14
##################################################



class orders(object):    
    def __init__(self, x1, x2, y1,y2):  
        self.x1 = x1
        self.ycen = ycen

import pyfits
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

plt.close()


# READ FLAT FRAME
ffile = pyfits.open('flat.fits')
flat = ffile[0].data

# SMOOTH FLAT SLIGHTLY
sflat = ndimage.gaussian_filter(flat, 3)


# FIND EDGES @ CENTER OF CHIP
tflat = sflat[1995:2005,:]

dx = ndimage.sobel(tflat, 0)  # horizontal derivative
dy = ndimage.sobel(tflat, 1)  # vertical derivative
mag = np.hypot(dx, dy)  # magnitude
mag1d = np.sum(mag,axis=0)
x = np.arange(0,len(mag1d))
y = np.ones(len(mag1d))*2000.


mx = (np.diff(np.sign(np.diff(mag1d))) < 0).nonzero()[0] + 1 # local max
mg = mag1d[mx]


# THIS GETS RED PEAKS
pred  = np.max(mag1d) * 0.1
xcenr = mx[mag1d[mx] > pred]
ycenr = y[mag1d[mx]  > pred]
mcenr = mg[mag1d[mx] > pred]

# THIS GETS BLUE PEAKS
pblue = 2
xmax  = np.min(xcenr)
xcenb = mx[(mag1d[mx] > pblue)&(mx < xmax)]
ycenb = y[(mag1d[mx]  > pblue) & (mx < xmax)]
mcenb = mg[(mag1d[mx]  > pblue) & (mx < xmax)]


# PUT THEM TOGETHER, SHOULD BE 20 PEAKS
xcen=np.concatenate([xcenr,xcenb])
ycen=np.concatenate([ycenr,ycenb])
mcen=np.concatenate([mcenr,mcenb])
#	print len(xcen)
print xcen
print ycen


# PLOT EDGES
#plt.xlim(1550,1600)
plt.plot(x,mag1d)
plt.plot(x[mx], mag1d[mx], "o", label="max")
plt.plot(xcen,mcen, "ro")
plt.plot(x, np.ones(len(x))*pred)

plt.show()



esiorders = [orders() for n in xrange(10)]


for i in xrange(10):
   esiorders[i].x1=xcen[2*i]
   esiorders[i].x1=xcen[2*i+1]



