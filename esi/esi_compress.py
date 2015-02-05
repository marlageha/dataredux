##############################################################################
#                                                                            #
#   Combine orders into a single spectrum. In regions of overlap, weight     #
#   each order by its polynomial fit. Propagate error through entire         #
#   reduction process. Smooth out badly subtracted sky, setting the error    #
#   to 1000 in region where "manual" correction is necessary. Write spectra  #
#   to Final/spectra/ and signal-to-noise to Final/sig_2_noise/. Identify    #
#   some lines, if so desired.                                               #
#                                                                            #
#   Kareem El-Badry, 01/19/2015                                              #
##############################################################################

from __future__ import division
import numpy as np
import os
import pickle
import pyfits
import scipy.ndimage
import scipy.signal
from scipy.ndimage.filters import gaussian_filter1d
from scipy.ndimage.filters import sobel
from scipy.signal import argrelextrema
import warnings

def esi_compress(date):
    
    #load masks and wavelength solutions
    cen_mask = pickle.load(open(str(date)+'/Calibs/cen_mask_'+str(date)+'.p','rb')) 
    #solutions = pickle.load(open(str(date)+'/Calibs/lambda_solutions_'+str(date)+'.p', 'rb'))
    #pix_to_ang = solutions[0]
    solutions2d = pickle.load(open(str(date)+'/Calibs/solution2d_'+str(date)+'.p', 'rb'))
    esiorders = pickle.load(open(str(date)+'/Calibs/order_edges_'+str(date)+'.p', 'rb'))
    
    
    #get names of lamps, first read in log:
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
        alldata.append([filename[line],dateobs[line],objname[line],imgtype[line],
                        ra[line],dec[line],exptime[line],usable[line]])

    #Find good files:
    good = []
    for line in range(len(alldata)):
       if "yes" in alldata[line][7]:
           good.append(alldata[line])

    #Find list of objects 
    names = []
    for line in range(len(good)):
        if ("Object" in good[line][3] and float(good[line][6]) > 600) or ("*" in good[line][2]):
            names.append(good[line][2])
    objects = np.array(list(set(names)))
    objects.sort() #ascending order, modify in place
    
    #Make directory to hold spectrum files
    if not os.path.exists(str(date)+'/Final/spectra/'):
        os.makedirs(str(date)+'/Final/spectra/')
        
    if not os.path.exists(str(date)+'/Final/sig_2_noise/'):
        os.makedirs(str(date)+'/Final/sig_2_noise/')

    #TODO: Autodmate this. 
    full_slit = ['119887', '36363', '37836', '54655', '81315']

    for obj_id in objects: 
        print 'reducing ' +str(obj_id)+ '...' 

        #load object
        if obj_id in full_slit: 
            
            #should maybe use un_sky_subtracted. Then /reduced/_mean.fits
            img = pyfits.getdata(str(date)+'/Calibs/sky_sub/'+str(obj_id)+'_skysub.fits')
        else: 
            img = pyfits.getdata(str(date)+'/Calibs/sky_sub/'+str(obj_id)+'_skysub.fits')

        noise = pyfits.getdata(str(date)+'/Calibs/variance/'+str(obj_id)+'_noise.fits')
    
        #Get rid of bogus values
        med = np.median(img)
        noise_med = np.median(noise)
        
        noise_nan_mask = np.isnan(noise)
        nan_mask = np.isnan(img)
        

        
        img[nan_mask] = med 
        noise[noise_nan_mask] = noise_med 
        
        #Fit polynomials to each order:
        warnings.simplefilter('ignore', np.RankWarning)

        new_points = np.linspace(3800, 11000, 2*40960)

        spec = 0*new_points    
        error_spec = 0*new_points
    
        sum_weights = 0*new_points

        flat = pyfits.getdata(str(date)+'/Calibs/dome_flat_'+str(date)+'.fits')
        order_polys = []

        for num in range(len(cen_mask)):
            print "compressing order " + str(num) + '...'
    
            sol = solutions2d[num] # load 2d wavelength solution
            y, x = np.mgrid[:4096, :2045] # holds x and y coordinates of all pixels
            heights = img[cen_mask[num]] # the actual values of the pixels
            
            #smooth = gaussian_filter1d(heights, 20)
            #peaks = argrelextrema(smooth, np.greater)[0]
            #peak_lambdas = np.array([lambdas[i] for i in peaks])
            
            #get rid of nonsense 
            #bad_peaks = np.diff(peak_lambdas) < 0.5*np.median(np.diff(peak_lambdas))
            
            
            #   FIND WAVELENGTHS OF CENTRAL PIXELS IN EACH ROW
            bin_widths = (esiorders[num].xr(range(4096)) - esiorders[num].xl(range(4096)))/2
            midpoints = esiorders[num].xl(range(4096))+bin_widths
            bin_cens = sol(midpoints, range(4096))
            
            
            errors = noise[cen_mask[num]] # corresponding noise
        
            lambdas = sol(x[cen_mask[num]], y[cen_mask[num]]) # wavelenghts of pixels
            spacing = (lambdas.max() - lambdas.min())/4096
            
            bin_edges = []
            for i in range(len(bin_cens)):
                bin_edges.append(bin_cens[i]-0.5*spacing)
            bin_edges.append(bin_cens[i]+0.5*spacing) # last edge
    
            #breaks = np.linspace(lambdas.min(), lambdas.max(), 4095) 
            #bin_cens = [breaks[line] + 0.5*spacing for line in range(len(breaks)-1)] 
            
            #bin_cens = peak_lambdas
            #bin_cens = np.array(bin_cens) # give wavenlength of center of each pixel row

            binned = np.histogram(lambdas, bin_edges, weights = heights)[0]
            means = binned/np.histogram(lambdas, bin_edges)[0] # mean 
        
            #add errors in quadrature, as usual
            binned_errors = np.sqrt(np.histogram(lambdas, bin_edges, weights = errors**2)[0])
            error_means = binned_errors/np.histogram(lambdas, bin_edges)[0]
        
            #flat_heights = flat[cen_mask[num]]
            #binned_flat =  np.histogram(lambdas, breaks, weights = flat_heights)[0]
            #flat_means = binned_flat/np.histogram(lambdas, breaks)[0]
    
            #Fit polynomial - alterantive, can fit polynomial to flats
            x = bin_cens
            y = means #or flat_means
            
            poly_med = np.median(y)
            
            # gets rid of dead pixels and other nonsense.
            x = x[np.abs(y) < 10*poly_med] 
            y = y[np.abs(y) < 10*poly_med]
    
            
            done = False
            while not done:
                done = True
                fit = np.poly1d(np.polyfit(x, y, 10))
                resid = y - fit(x)
                std = np.std(resid)
                badindices = np.where(np.abs(resid) > 3*std)[0]
                
                if badindices.size > 0 and len(x) - len(badindices) > 3500:
                    done = False
                    x = np.delete(x, badindices)
                    y = np.delete(y, badindices)
                    
            # fit = np.poly1d(np.polyfit(x, y, 10))
            order_polys.append(np.poly1d(fit))    
    
            #interpolate each order to full spectrum
            #Nearly equivalent to adding across order
            interp_order = np.interp(new_points, bin_cens, means)
            interp_order_error = np.sqrt(np.interp(new_points, bin_cens, error_means**2))
        
            #interpolate weights
            interp_weights = np.interp(new_points, bin_cens, order_polys[num](bin_cens))
            interp_weights[(new_points < lambdas.min()) + (new_points > lambdas.max())] = 0 #useless data
            
            interp_weights[interp_weights < 0.1*interp_weights.max()] = 0 
        
            sum_weights = sum_weights + interp_weights**2 #because we've multiplied by the weight twice
    
            interped = interp_weights * interp_order 
            interped[np.isnan(interped)] = 0
            spec = spec + interped # add to spec for each order (each iteration of for loop)
        
            interped_error = interp_weights * interp_order_error
            error_spec = np.sqrt(error_spec**2 + interped_error**2) #always combine in quadrature
            #error_spec[np.isnan(error_spec)] = 0
            
            
        spec = spec/sum_weights
        error_spec = error_spec/sum_weights

        new_lambdas = np.linspace(3800,10000, 30000) #throwing out edges on both sides
        new_spec = np.interp(new_lambdas, new_points, spec)
        new_error = np.sqrt(np.interp(new_lambdas,new_points, error_spec**2))
    
        new_spec[np.isnan(new_spec)] = 1
    

        if obj_id not in full_slit:
            #Now get rid of parts that that are too steep 
            
            #Throw out the parts that aren't gaussian
            filtered = scipy.ndimage.gaussian_filter1d(new_spec, 3)
            resid = new_spec - filtered
            new_spec[np.abs(resid > 0.2)] = 1

            #deriv = sobel(new_spec)


            #new_spec[np.abs(deriv3) > 1] = 1

            
            new_error[new_spec == 1] = 1000
    
        if obj_id in full_slit:
        
            filtered = scipy.ndimage.gaussian_filter1d(new_spec, 3)
            resid = new_spec - filtered
            new_spec[np.abs(resid > 0.3)] = 1

            #deriv = sobel(new_spec)
            #deriv2 = sobel(deriv)
            #deriv3 = sobel(deriv2)
            new_error[new_spec == 1] = 1000
            
            #new_spec[np.abs(deriv3) > 1] = 1
            #new_spec[np.abs(deriv2) > 0.8] = 1
            #new_spec[np.abs(deriv) > 0.7] = 1
            
        
        spec_file = open(str(date)+'/Final/spectra/'+str(obj_id)+'_spectrum.dat', 'w')
        sn_file = open(str(date)+'/Final/sig_2_noise/'+str(obj_id)+'_signoise.dat', 'w')
        for line in range(len(new_lambdas)):
            spec_file.write(" ".join([str(new_lambdas[line]), str(new_spec[line]), '\n']))
            sn_file.write(" ".join([str(new_lambdas[line]), str((new_spec/new_error)[line]), '\n']))
        spec_file.close()
        sn_file.close()
        
        '''
        #Interesting lines
        lines = ["Ca II K", "Ca II H", "$H\\beta$", "Mg Ib", "Mg Ib", "Mg Ib", 
                 "$H\\alpha$", "Ca II", "Ca II", "Ca II", "[S II]", "[S II]",
                  "[N II]", "[N II]", "O I", "Na D2", "Na D1","Ca + Fe"]
        line_lambdas = [3933.663, 3968.468, 4861.32, 5167.32, 5172.68, 5183.60, 
                        6562.80, 8498.0, 8542.1, 8662.2, 6716.42, 6730.78, 6583.6, 
                        6548.1, 8446, 5889.950, 5895.924, 6495]
                
        line_lambdas = np.array(line_lambdas)
                
        
        if obj_id == '35979':
            z = 0.019651
        elif obj_id == '55500':
            z = 0.0207593
        elif obj_id == '113209':
            z = 0.0118022
        elif obj_id == '119887':
            z = 0.0223354
        elif obj_id == '120659':
            z = 0.027492
        elif obj_id == '122277':
            z = 0.0244582
        elif obj_id == '20700':
            z = 0.0241752
        elif obj_id == '3478':
            z = 0.0244122
        elif obj_id ==  '36363':
            z = 0.00470846
        elif obj_id == '37836':
            z = 0.00950647
        elif obj_id == '38329':
            z = 0.0103071
        elif obj_id == '38465':
            z = 0.0258144	
        elif obj_id == '50778':
            z = 0.00413376
        elif obj_id == '51306':
            z = 0.025976
        elif obj_id == '54655':
            z = 0.00232896
        elif obj_id == '67565':
            z = 0.0231906
        elif obj_id == '81315': 
            z = 0.00161523
        else: 
            z = 0.02
            print "redshift unknown... guessing z = 0.02"
        '''