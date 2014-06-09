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
   "New_J_Flat.fits," "New_K_Flat.fits" are written to /Calibs/. (These are the best flats.)
   > whirc.whirc_rmpupil()

8. Inside Calibs, make a directory "sky." E.g. there should be a directory Calibs/sky.
   Make a sky image for each flat, each night, each object by median combining science 
   images.
   > whirc.whirc_mksky(obj_id) 
   i.e. if you want to make a sky image for 55500, > whirc.whirc_mksky(55500)

9. Inside Calibs, make a directory “reduced.” For a particular obj_id, apply bias, dark, 
   sky subtraction, and flat. Write reduced images to Calibs/reduced.
   > whirc.whirc_reduce(obj_id)

10. Inside Calibs, make a directory “shifted.” For a particular obj_id, align images and 
    median combine. Write aligned images and median combined image to Calibs/shifted.
    > whirc.whirc_align(obj_id)
 

