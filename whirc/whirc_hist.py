#Run in directory above Raw/, Calibs/, Final/
#NSA catalog can be downloaded at http://www.nsatlas.org/data. 
#Save to Downloads/

def whirc_hist():
    
    import pyfits
    import numpy as np
    import matplotlib.pyplot as plt
    import scipy.ndimage as snd 
    import os

    # FIND FILES, MAKE LISTS OF THEIR PATHS
    files = []
    for file in os.listdir(os.getcwd()+'/Calibs/galfits/'):
        if "block" in file:
            files.append(file)
            
    Jfiles = []    
    for line in range(len(files)):
        if "J" in files[line]:
            Jfiles.append('Calibs/galfits/'+files[line])

    Kfiles = []    
    for line in range(len(files)):
        if "K" in files[line]:
            Kfiles.append('Calibs/galfits/'+files[line])

    Rfiles = []    
    for line in range(len(files)):
        if "r" in files[line]:
            Rfiles.append('Calibs/galfits/'+files[line])
            
    #GET DATA FROM MODEL HEADERS (obj_id_J_block.fits[2])
    Jmag = []
    Jre = []
    Jind = []
    Jar = []

    for line in range(len(Jfiles)):
        img = pyfits.open(Jfiles[line])
        hdr = img[2].header
    
        profile = hdr['COMP_1']
        Jmag.append(float(hdr['1_MAG'].split()[0]))
        Jre.append(float(hdr['1_RE'].split()[0]))
        Jind.append(float(hdr['1_N'].split()[0]))
        Jar.append(float(hdr['1_AR'].split()[0]))
    
    Kmag = []
    Kre = []
    Kind = []
    Kar = []

    for line in range(len(Kfiles)):
        img = pyfits.open(Kfiles[line])
        hdr = img[2].header
    
        profile = hdr['COMP_1']
        Kmag.append(float(hdr['1_MAG'].split()[0]))
        Kre.append(float(hdr['1_RE'].split()[0]))
        Kind.append(float(hdr['1_N'].split()[0]))
        Kar.append(float(hdr['1_AR'].split()[0]))
    
    Rmag = []
    Rre = []
    Rind = []
    Rar = []

    for line in range(len(Rfiles)):
        img = pyfits.open(Rfiles[line])
        hdr = img[2].header
    
        profile = hdr['COMP_1']
        Rmag.append(float(hdr['1_MAG'].split()[0]))
        Rre.append(float(hdr['1_RE'].split()[0]))
        Rind.append(float(hdr['1_N'].split()[0]))
        Rar.append(float(hdr['1_AR'].split()[0]))

    Jrre = np.array(Jre)
    Krre = np.array(Kre)
    Rrre = np.array(Rre)
    
    #Pixel scale for WHIRC and SDSS 
    Jre_arcsec = Jrre*0.1
    Kre_arcsec = Krre*0.1
    Rre_arcsec = Rrre*0.4 
    
   
    '''    
    #PLOT SOME STUFF. USE AS MUCH AS YOU WANT
    
    f = plt.figure()
    ax = f.add_subplot(121)
        
    plt.hist(Jind)
    plt.title('Sersic Index Fits: J Filter')
    plt.xlabel('Sersic Index')
    plt.ylabel('frequency')

    ax = f.add_subplot(122)
        
    plt.hist(Kind)
    plt.title('Sersic Index Fits: K Filter')
    plt.xlabel('Sersic Index')
    plt.ylabel('frequency')
    plt.savefig('Calibs/pics/index_hist.png')
    
    g = plt.figure()
    
    ax = g.add_subplot(311)
    plt.hist(Jmag)
    plt.title('Integrated Magnitude Fits: J Filter')
    plt.xlabel('Magnitude')
    plt.ylabel('frequency')

    ax = g.add_subplot(312)
    plt.hist(Kmag)
    plt.title('Integrated Magnitude Fits: K Filter')
    plt.xlabel('Magnitude')
    plt.ylabel('frequency')
    plt.savefig('Calibs/pics/mag_hist.png')
    
    ax = g.add_subplot(313)
    plt.hist(Rmag)
    plt.title('Integrated Magnitude Fits: SDSS r-band')
    plt.xlabel('Magnitude')
    plt.ylabel('frequency')
    plt.savefig('Calibs/pics/mag_hist.png')
    plt.tight_layout()
    
    h = plt.figure()
    ax = h.add_subplot(121)

    plt.hist(Jre)
    plt.title('Effective Radius: J Filter')
    plt.xlabel('$R_{e}$ (pixels)')
    plt.ylabel('frequency')

    ax = h.add_subplot(122)
        
    plt.hist(Kre)
    plt.title('Effective Radius: K Filter')
    plt.xlabel('$R_{e}$ (pixels)')
    plt.ylabel('frequency')
    plt.savefig('Calibs/pics/re_hist.png')
    '''
    
    
    #NOW GETTING INFO FROM THE NSA CATALOG
    
    NSA_cat = pyfits.open('nsa_v0_1_2.fits')

    table = NSA_cat[1].data

    NSAID = table['NSAID']
    sersic_n = table['sersic_n']
    r50 = table['SERSIC_TH50']
    mass_light = table['MTOL']
    mass_light = [mass_light[i][6] for i in range(len(mass_light))]
    Z = table['Z']
    mass = table['MASS'] #from K-correct
    nsa_ba = table['SERSIC_BA']

    obj_ids = [119887, 120659, 121130, 122277, 20700, 35979, 36363, 38329,
               51306, 54655, 55500, 57467, 67565, 77610]
               
    locs = []
    for line in range(len(obj_ids)):
        locs.append(int(np.where(NSAID == obj_ids[line])[0]))

    sersic_ns = [sersic_n[line] for line in locs]
    r50s = [r50[line] for line in locs]
    m_to_l = [mass_light[line] for line in locs]
    redshift =  [Z[line] for line in locs]
    stellar_mass = [mass[line] for line in locs]
    nsa_rar = [nsa_ba[line] for line in locs]

    #for line in range(len(nsa_rar)):
        #print nsa_rar[line]
        
    
    #PLOT SOME STUFF (JUST FROM GALITS, NOT FROM NSA CATALOG RIGHT NOW)
    f = plt.figure()
    ax = f.add_subplot(111)

    plt.plot(Jre_arcsec, Jind, 'o', color = 'black', label = 'Infrared (WHIRC)',
            markersize = 8)
    plt.plot(Kre_arcsec, Kind, 'o', color = 'black', markersize = 8)
    #plt.plot(r50s, sersic_ns, 'o', color = 'green', label = 'NSA catalog')
    plt.plot(Rre_arcsec, Rind, 'o', color = 'red', label = 'Optical (SDSS)',
             markersize = 8)
    plt.xlabel('Effective Radius (arcsec)', fontsize = 14)
    plt.ylabel('Sersic Profile Steepness', fontsize = 14)
    plt.title('SDSS / WHIRC comparison', fontsize = 20)
    ax.set_xlim([0,5.33])
    ax.legend()
    savefig('Whirc_vs_sdss.png', bbox_inches = 'tight')
    
    '''
    g = plt.figure()
    ax = g.add_subplot(111)
    plt.plot(Rmag, Rar, 'o', color = 'red', label = 'SDSS Optical Fit')
    plt.plot(Jmag, Jar, 'o', color = 'gray', label = 'WHIRC Infrared Fit')
    plt.plot(Kmag, Kar, 'o', color = 'gray', label = 'WHIRC Infrared Fit')
    plt.xlabel('Magnitude')
    plt.ylabel('Axis ration (b/a)')
    plt.title('SDSS / WHRIC comparison')
    ax.legend(loc = 4, fontsize = 11)
    plt.savefig('Calibs/pics/mag_fit.png')'''
    
    
    