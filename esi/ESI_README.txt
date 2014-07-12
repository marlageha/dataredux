1) Make directory ESI. Inside ESI, make directories Raw/, Calibs/, Final/, Logs/.
   Inside Logs/, make directory pics/.  Move data to Raw/. Extract data from headers
   in Python,
   >> import esi
   >> esi.esi_setup()

2) Cross-check logs from headers. I've already done this and saved the results in 
   esi_info.dat. Move esi_info.dat to Logs/ (from GitHub) AFTER running esi_setup. 
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
