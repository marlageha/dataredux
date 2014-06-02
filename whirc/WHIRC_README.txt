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
   b) I then used > whirc.whirc_info() to create an easy-to-edit digital log of the  
   observations. I then manually corrected this log; the whirc_info.dat file is the
   most up-to-date version. Don't run whirc_info(); it will overwrite the manually 
   corrected log with one that hasn't been corrected yet. 
   c) Move whirc_info.dat to directory above /Raw, /Calibs, /Final.
  
3. Co-add all bias frames, bias.fits is written to /Calibs
   > whirc.whirc_mkbias()
   
4. Co-add all dark frames, dark.fits is written to /Calibs
   > whirc.whirc_mkdark()
   
5. Co-add all flats, subtract off_flats from on_flats, normalize. J_master_flat.fits
   and K_master_flat.fits are written to /Calibs
   > whirc.whirc_sort()

