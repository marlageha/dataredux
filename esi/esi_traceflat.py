#! /usr/bin/python
##################################################
#  ESI_traceflat
#
#  USING FLAT FIELD, FIND X/Y FOR ORDER EDGES  
#
#  MG 4/14
##################################################


class orders:
    def __init__(self, xsize, ysize):  
       self.xr = np.zeros(xsize)
       self.xb = np.zeros(xsize)
       self.yr = np.zeros(ysize)
       self.yb = np.zeros(ysize)


import mgesi
import pyfits
import matplotlib.pyplot as plt
import numpy as np
from scipy import ndimage

def esi_traceflat():

	# READ FLAT FRAME
	ffile = pyfits.open('Calibs/flat.fits')
	flat = ffile[0].data

	# SMOOTH FLAT SLIGHTLY
	sflat = ndimage.gaussian_filter(flat, 3)


	# FIND CENTRAL PEAKS
	xcen, ycen = mgesi.esi_trace_cen()
	print xcen
	print ycen

	# CREATE ARRAT FOR ESI ORDERS -- REDDEST ORDER IS FIRST
	xtrace=np.zeros([len(xcen),4096])
	ytrace=np.zeros([len(ycen),4096])

	xtrace[:,2000] = xcen
	ytrace[:,2000] = ycen


	#  TRACE IN FOUR PARTS, STARTING FROM CENTRAL ROW
	#  (1)  WORK DOWN FOR FULL ARRAY
	#  (2)  WORK UP FOR FULL ARRAY
	#  (3)  FINISH WORKING UP, REMOVE BLUE-MOST ORDER
	#  (4)  FINISH WORKING DOWN, REMOVE BLUE-MOST ORDER


	# (1) WORK DOWN FOR FULL ARRAY
	for row in xrange(1999,920,-1):
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
			mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1
	#		plt.plot(xarr,local_peak,'r+')#
	#		plt.show()

			xtrace[n,row] = mx-5+xt
			ytrace[n,row] = row


	# (2)  WORK UP FOR FULL ARRAY
	for row in range(2001,3550):

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
			mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1

			xtrace[n,row] = mx-5+xt
			ytrace[n,row] = row
		


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
			mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1
			print row, n ,mx
			xtrace[n,row] = mx-5+xt
			ytrace[n,row] = row

		if (row < 3700):
			for n in xrange(0,2,1):
				xt = np.int(np.round(xtrace[n,row-1]))#
				xarr = np.arange(xt-5,xt+5)
				local_peak = mag1d[xarr]#
				mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1
				print row, n ,mx

				xtrace[n,row] = mx-5+xt
				ytrace[n,row] = row


	# (4) FINISH WORKing DOWN FOR FULL ARRAY
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
			mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1
	#		plt.plot(xarr,local_peak,'r+')#
	#		plt.show()

			xtrace[n,row] = mx-5+xt
			ytrace[n,row] = row
		if (row > 550):
			for n in xrange(14,16,1):
				xt = np.int(np.round(xtrace[n,row+1]))
				xarr = np.arange(xt-5,xt+5)
				local_peak = mag1d[xarr]
				mx = (np.diff(np.sign(np.diff(local_peak))) < 0).nonzero()[0] + 1

				xtrace[n,row] = mx-5+xt
				ytrace[n,row] = row


#####################################################
 	plt.close()
	plt.imshow(sflat,cmap='gray')
	for n in xrange(0,19,1):
#		plt.plot(xtrace[n],ytrace[n])
		p = np.polyfit(ytrace[n,:],xtrace[n,:],3)
		print p
		xp = np.poly1d(p)
		y=np.arange(0,4096,1)	
		print xp(y)
		plt.plot(xp,y)
	#return xtrace, yt#race


#	plt.plot(xtrace,ytrace,'r.')
	plt.show()


#	esiorders = [orders(2000,2000) for n in xrange(10)]


#	for i in xrange(10):
#	   esiorders[i].xr=xcen[2*i]
#	   esiorders[i].xb=xcen[2*i+1]


#	print esiorders[0].xr[100]
	# PARSE INTO CCD ECHELL ORDERS
	#for n in xrange(0,19,2):
	#	order[n/2].x1 = xtrace[n,:]
	#	order[n/2].x2 = xtrace[n+1,:]
	#	order[n/2].y1 = ytrace[n,:]
	#	order[n/2].y2 = ytrace[n+1,:]


#	esiorders = [orders() for n in xrange(10)]#
#	for i in xrange(10):
#	   esiorders[i].x1=xcen[2*i]
#	   esiorders[i].x1=xcen[2*i+1]



