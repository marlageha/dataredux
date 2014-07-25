1) Make directory ESI. Inside ESI, make directories Raw/, Calibs/, Final/, Logs/.

   Inside Calibs/, make directory pics/.  Move data to Raw/. Extract data from headers
   in Python,
   >> import esi
   >> esi.esi_setup()

2) Cross-check logs from headers. I've already done this and saved the results in 
   esi_info.dat. Move esi_info.dat to Logs/ (from GitHub) AFTER running esi_setup.
   (Otherwise running esi_setup() will overwrite it with something less useful). 
   This tells which files are bad.

3) Find bias files from esi_info and median combine. Write bias.fits to Calibs/.
   >> esi.esi_mkbias()

4) Find flats from esi_info and median combine. Write dome_flat.fits and
   pinhole_flat.fits to Calibs/.
   >> esi.esi_mkflats()

5) From flat, find orders and fit a polynomial to them. (This also requires the esi 
   module esi_trace_cen() to be in the pythonpath directory). Normalize each order 
   by dividing it by its polynomial. Write Calibs/norm_flat.fits. Create masks for
   the orders and the background; write orders_mask.p and background_mask.p to 
   Calibs/.  
   >> esi.esi_traceflat()

6) Make a directory Calibs/cosmicless. From http://obswww.unige.ch/~tewes/cosmics_dot_py/,
   (or just from the git_hub repository), download cosmics.py and save it to the 
   directory of your PYTHONPATH variable in .bash_profile. Remove cosmic rays from each 
   of the raw science images and the lamps. Write each file to Calibs/cosmicless.
   >> esi.esi_cosmic()

7) In Calibs, make directories reduced/ and variance/. Read in the spectra from
   each object and average with rejection. Write each one to Calibs/reduced/objid_med.fits.
   Do the same for lamps. Make a variance image for each science image whose pixels
   are 1/(sigma)^2. Write to Calibs/variance/. 
   >> esi.esi_reduce() 

8) In Calibs, make directory Calibs/lamp_peaks. Read in spectra from each object, 
   sum over x coordiates in each order, and find the peaks (to +/- 1 pixel) in 
   each order. Write peak files to Calibs/lamp_peaks/.
   >> esi.esi_roughpeaks()

9) In Calibs, make directory Calibs/line_lists. Copy the directory order_lists/ from 
   GitHub into line_lists. Fit a wavelength solution to each order using line lists.
   Write wavelength solutions and good lines to lambda_solutions.p
   >> esi.esi_lambda()

10) In Calibs, make directory Calibs/sky_sub. Using masks made in esi_traceflat(), 
    measure sky values at edges of orders and subtract. Only very rough at the moment.
    Write (roughly) sky-subtracted stuff to Calibs/sky_sub.
    >> esi.esi_skysub()

11) Fit a gaussian to each peak in the lamps and record its center, which is more
    accurate than the estimate from esi_roughpeaks(). Compare to rough line lists
    and only keep the ones that we've already found and know the actual wavelengths 
    of. Write improved line lists to Calibs/line_lists/order_lists/pixels/better_
    order[obj_id].dat. 
    >> esi.esi_gauss()

12) Fit a 2-dimensional wavelength solution to each order. Break each order into
    ~15 columns that are 10 pixels wide. Find all the peaks in each column, so we
    can see how the peak moves in y-coordinate as we move across the order. This give us
    a whole bunch of points in (x,y,wavelength) space, to which we fit a 2d polynomial 
    surface. Write the polynomial solutions to Calibs/solution2d.p
    >> esi.esi_solution2d()

