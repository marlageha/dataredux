---------------------------------------
Set environment variable:

bash$ emacs .bash_profile
add:

	PYTHONPATH="$PYTHONPATH:<full path to /dataredux/ on local directory>" 
	export PYTHONPATH

0. Create directories /Raw /Calibs and /Final
   Move rawdata into /Raw. Run all scripts in directory above Raw, Calibs, Final.

1. In python:
   > import whirc

2. Run script to organize data.
   a) Read all headers and create fits table 'whirc_data.fits'
   >  whirc.whirc_setup() 
   b) Move whirc_info.dat, the digital observing log, to directory 
      above /Raw, /Calibs, /Final.
   c) whirc_info.py is no longer needed, it was used to create a draft of the observing
      log which was then edited and cross-checked with the data. 
  
3. Co-add all bias frames, bias.fits is written to /Calibs
   > whirc.whirc_mkbias()
   
4. Co-add all dark frames, dark.fits is written to /Calibs
   > whirc.whirc_mkdark()
   
5. Co-add all flats (subtracting appropriate dark and bias), subtract off_flats from
   on_flats, normalize. J_master_flat.fits and K_master_flat.fits are written to /Calibs
   > whirc.whirc_sort()

6. From http://www.noao.edu/kpno/manuals/whirc/datared.html, download 'bpix.whirc.fits'.
   (the bad pixel map) Move to directory above Raw, Calibs, Final. 
   
7. Smooth background of master flats and subtract to get rid of pupil ghost. 
   "New_J_Flat.fits," "New_K_Flat.fits" are written to /Calibs/. (These are the best 
   flats.)
   > whirc.whirc_rmpupil()

8. Inside Calibs, make a directory "sky." E.g. there should be a directory Calibs/sky.
   Make a sky image for each flat, each night, each object by median combining science 
   images.
   > whirc.whirc_mksky(obj_id) 
   i.e. if you want to make a sky image for 55500, > whirc.whirc_mksky(55500)

9. Inside Calibs, make a directory “reduced.” 
   From http://www.astro.yale.edu/dokkum/lacosmic/, download and install the “cosmics.py”  
   module. Move cosmics.py to dataredux/whirc
   For a particular obj_id, apply bias, dark, sky subtraction, and flat. Make copy of 
   images w/o cosmic rays (for alignment). Write reducedimages to Calibs/reduced.
   > whirc.whirc_reduce(obj_id)

10. Inside Calibs, make a directory “starfields” For a particular obj_id, find brightest 
    stars in all exposures and write starfield to a pickle file.   
    > whirc.whirc_align(obj_id)

11. Inside Calibs, make a directory “shifted.” Use starfields from last step to align 
    images for each object, night, and filter and median combine. Write to Calibs/shifted.
    > whirc.whirc_shiftcalc(obj_id)

12. Find master science images and coadd them (including shifts) for objects with more 
    than one night’s worth of data. Write finals J and K images to /Final/.
    > whirc.whirc_combine(obj_id)

13. Get GALFIT from http://users.obs.carnegiescience.edu/peng/work/galfit/galfit.html. 
    Install and set path to directory where GalFit is. Following the template in
    EXAMPLE.INPUT, make an input file for each galaxy. In Calibs, make directories 
    Calibs/galfits and Calibs/inputs; put input files in /inputs. Guess fittable 
    parameters and run GALFIT. Fit Sersic profiles and sky. 

    galfit obj_id_inputj 

14. Extract plots and images from GalFit output. In Calibs, create directory Calibs/pics. 
    > whirc.whirc_draw(obj_id)

15. Plot some stuff. To view in context with the NSA catalog, download the NSA catalog to 
    Downloads/ from http://www.nsatlas.org/data. Then run
    > whirc.whirc_hist()

16. Plot this dataset in context of Geha 2012. Download NSA catalog for M < 1e10 M_sun
    from https://www.dropbox.com/s/b1kqrucrrew8num/geha2012.fits.gz?n=229870315.Then run
    > whirc.whirc_plt()
 

