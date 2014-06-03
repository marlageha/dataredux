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
 

