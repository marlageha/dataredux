#Plots some stuff from this data sample to put stuff in context!

from __future__ import division
import pyfits
import matplotlib.pyplot as plt
import numpy as np
import pdb

def whirc_plot():
    
    #Here I get the data from the Geha 2012 paper, which has M < 1e10 M_sun
    hdulist = pyfits.open('Calibs/geha2012.fits.gz')
    ntable = hdulist[1].data

    NSAID = ntable['NSAID']
    ndist = ntable['dhost']
    nmass = ntable['mass']
    nha = ntable['haew']
    nd4000 = ntable['d4000']

    #FIND ISOLATED STUFF AND NOT ISOLATED STUFF
    #mass
    iso_mass = [nmass[line] for line in range(len(nmass)) if ndist[line] > 1.5]
    niso_mass = [nmass[line] for line in range(len(nmass)) if ndist[line] < 1.5]

    #D4000
    iso_d4000 = [nd4000[line] for line in range(len(nmass)) if ndist[line] > 1.5]
    niso_d4000 = [nd4000[line] for line in range(len(nmass)) if ndist[line] < 1.5]

    #Ha
    iso_ha = [nha[line] for line in range(len(nmass)) if ndist[line] > 1.5]
    niso_ha = [nha[line] for line in range(len(nmass)) if ndist[line] < 1.5]
    
        
    #Define Quenching Criterion
    a, b = 0.6, 0.1
    def dbreak(m):
        return a + b*m
 

    #FIND QUENCHED STUFF AND NOT QUENCHED STUFF
    #Dn4000
    quenched_d4000 = [nd4000[line] for line in range(len(nd4000)) if nd4000[line] > dbreak(nmass[line]) and nha[line] < 2]
    unquen_d4000 = [nd4000[line] for line in range(len(nd4000)) if nd4000[line] < dbreak(nmass[line]) and nha[line] > 2]

    #Meet one criterion but not the other
    K_A_d1 = [nd4000[line] for line in range(len(nd4000)) if nha[line] < 2 and nd4000[line] < dbreak(nmass[line])]
    K_A_d2 = [nd4000[line] for line in range(len(nd4000)) if nha[line] > 2 and nd4000[line] > dbreak(nmass[line])]

    #H_alpha
    quenched_ha = [nha[line] for line in range(len(nha)) if nha[line] < 2 and nd4000[line] > dbreak(nmass[line])]
    unquen_ha = [nha[line] for line in range(len(nha)) if nha[line] > 2 and nd4000[line] < dbreak(nmass[line])] 

    #Meet one criterion but not the other
    K_A_ha1 = [nha[line] for line in range(len(nha)) if nha[line] < 2 and nd4000[line] < dbreak(nmass[line])]
    K_A_ha2 = [nha[line] for line in range(len(nha)) if nha[line] > 2 and nd4000[line] > dbreak(nmass[line])]

    #Mass
    quenched_mass = [nmass[line] for line in range(len(nmass)) if nd4000[line] > dbreak(nmass[line]) and nha[line] < 2]
    unquen_mass = [nmass[line] for line in range(len(nmass)) if nd4000[line] < dbreak(nmass[line]) and nha[line] > 2]

    #Meet one criterion but not the other
    K_A_m1 = [nmass[line] for line in range(len(nmass)) if nha[line] < 2 and nd4000[line] < dbreak(nmass[line])]
    K_A_m2 = [nmass[line] for line in range(len(nmass)) if nha[line] > 2 and nd4000[line] > dbreak(nmass[line])]


    #Find the objects for which I have WHRIC or ESI data 
    obj_ids = [113209, 119887, 120659, 122277, 20700, 3478, 35979, 36363,
               37836, 38329, 38465, 46677, 50778, 51306, 54655, 55500,  
               67565, 81315]
           
    locs = []
    for line in range(len(obj_ids)):
        locs.append(int(np.where(NSAID == obj_ids[line])[0]))

    sample_dist = [ndist[line] for line in locs]
    sample_mass = [nmass[line] for line in locs]
    sample_ha = [nha[line] for line in locs]
    sample_d4000 = [nd4000[line] for line in locs]
    sample_names = [NSAID[line] for line in locs]  


    #Find quenched objects in Geha 2012
    dwarf_quenched_ha = [nha[a] for a in range(len(nha)) if nha[a] < 2.0]
    dwarf_notquen_ha = [nha[a] for a in range(len(nha)) if nha[a] > 2.0]

    dwarf_quenched_m = [nmass[a] for a in range(len(nmass)) if nha[a] < 2.0]
    dwarf_notquen_m = [nmass[a] for a in range(len(nmass)) if nha[a] > 2.0]


    #Here I get the data from the full NSA catalog
    NSA_cat = pyfits.open('../Downloads/nsa_v0_1_2.fits')
    table = NSA_cat[1].data

    NSAID = table['NSAID']
    mass = table['MASS'] 
    d4000 = table['d4000']
    ha = table['haew']

    #Only plotting higher masses from NSA catalog
    logmass = [np.log10(mass[a]) for a in range(len(mass)) if 1e10 < mass[a]]
    d4000 = [d4000[a] for a in range(len(mass)) if 1e10 < mass[a]]
    ha = [ha[a] for a in range(len(mass)) if 1e10 < mass[a]]


    #Find quenched objects in NSA catalog
    big_quenched_ha = [ha[a] for a in range(len(ha)) if ha[a] < 2.0]
    big_notquen_ha = [ha[a] for a in range(len(ha)) if ha[a] > 2.0]

    big_quenched_m = [logmass[a] for a in range(len(logmass)) if ha[a] < 2.0]
    big_notquen_m = [logmass[a] for a in range(len(logmass)) if ha[a] > 2.0]

    #Define Quenching Criterion
    a, b = 0.6, 0.1
    def dbreak(m):
        return a + b*m
    
    
    #Recreate Geha 2012 Figure 2
    f = plt.figure()
    ax = f.add_subplot(111)
    xrang = np.arange(5, 12, 0.01)

    plt.plot(logmass, d4000, '.', color = '0.75', markersize = 1, label = 'NSA catalog')
    plt.plot(nmass, nd4000, '.', color = 'black', markersize = 1, label = 'Geha 2012')
    plt.plot(xrang, dbreak(xrang), 'r', label = '$D_n4000$ break')
    plt.xlabel('$\log_{10}[M_{stellar}$ $(M_{\odot})]$')
    plt.ylabel('$D_n4000$')
    plt.title('$D_n4000$ Quenching Criterion')
    ax.set_ylim([0.8, 2.2])
    ax.set_xlim([11.2, 7.6])
    ax.legend(loc = 1, frameon = False, markerscale = 8, fontsize = 11)
    plt.savefig('Calibs/pics/geha_2012_fig2')

    #Plot Mass vs HaEW
    g = plt.figure()
    ax = g.add_subplot(111)

    plt.plot(dwarf_quenched_m, dwarf_quenched_ha, '.', color = 'red', markersize = 1, label = 'Quenched ($H\\alpha < 2$)')
    plt.plot(dwarf_notquen_m, dwarf_notquen_ha, '.', color = 'black', markersize = 1, label = 'Star Forming')
    plt.plot(big_quenched_m, big_quenched_ha, '.', color = 'red', markersize = 1)
    plt.plot(big_notquen_m, big_notquen_ha, '.', color = '0.75', markersize = 1)
    plt.xlabel('$\log_{10}[M_{stellar}$ $(M_{\odot})]$')
    plt.ylabel('$H\\alpha$ EW')
    plt.title('$H\\alpha$ EW Quenching Criterion')
    ax.set_ylim([0, 100])
    ax.set_xlim([11.2, 7.6])
    ax.legend(loc = 2, frameon = False, markerscale = 8, fontsize = 11)
    plt.savefig('Calibs/pics/geha_2012_mass_vs_HaEW')


    #Plot Ha vs D4000, showing quenched, K+A, and not quenched
    h = plt.figure()
    ax = h.add_subplot(111)

    plt.plot(quenched_ha, quenched_d4000, '.', color = 'magenta', markersize = 1, label = 'Quenched')
    plt.plot(K_A_ha1, K_A_d1, '.', color = 'cyan', markersize = 2, label = 'K+A')
    plt.plot(unquen_ha, unquen_d4000, '.', color = 'black', markersize = 1, label = 'Not Quenched')
    plt.plot(K_A_ha2, K_A_d2, '.', color = 'cyan', markersize = 2)
    plt.plot(sample_ha, sample_d4000, 'o', color = 'red', markersize = 5, label = 'This Sample' )
    plt.xlabel('$H\\alpha$ EW')
    plt.ylabel('$D_n4000$')
    plt.title('Sample in Context')
    ax.set_ylim([0.8, 2.2])
    ax.set_xlim([-1.5, 5])
    ax.legend(loc = 1, frameon = True, markerscale = 1, fontsize = 11)
    plt.savefig('Calibs/pics/geha_2012_haEW_vs_d40001')


    #Plot mass vs D4000, showing quenched, K+A, and not quenched
    j = plt.figure()
    ax = j.add_subplot(111)

    plt.plot(quenched_mass, quenched_d4000, '.', color = 'magenta', markersize = 1, label = 'Quenched')
    plt.plot(K_A_m1, K_A_d1, '.', color = 'cyan', markersize = 2, label = 'K+A')
    plt.plot(K_A_m2, K_A_d2, '.', color = 'cyan', markersize = 2)
    plt.plot(unquen_mass, unquen_d4000, '.', color = 'black', markersize = 1, label = 'Not Quenched')
    plt.plot(sample_mass, sample_d4000, 'o', color = 'red', markersize = 5, label = 'This Sample' )
    plt.xlabel('$\log_{10}[M_{stellar}$ $(M_{\odot})]$')
    plt.ylabel('$D_n4000$')
    plt.title('Sample in Context')
    ax.set_ylim([0.8, 2.2])
    ax.set_xlim([10, 6.95])
    ax.legend(loc = 1, frameon = True, markerscale = 2, fontsize = 11)
    plt.savefig('Calibs/pics/geha_2012_mass_vs_d4000')


    #FIND QUENCHED STUFF AND NOT QUENCHED STUFF (IN ISO vs NON_ISO)
    #Dn4000
    quenched_d4000_iso = [iso_d4000[line] for line in range(len(iso_d4000)) if iso_d4000[line] > dbreak(iso_mass[line]) and iso_ha[line] < 2]
    unquen_d4000_iso = [iso_d4000[line] for line in range(len(iso_d4000)) if iso_d4000[line] < dbreak(iso_mass[line]) and iso_ha[line] > 2]

    quenched_d4000_niso = [niso_d4000[line] for line in range(len(niso_d4000)) if niso_d4000[line] > dbreak(niso_mass[line]) and niso_ha[line] < 2]
    unquen_d4000_niso = [niso_d4000[line] for line in range(len(niso_d4000)) if niso_d4000[line] < dbreak(niso_mass[line]) and niso_ha[line] > 2]


    #Meet one criterion but not the other
    K_A_d1_iso = [iso_d4000[line] for line in range(len(iso_d4000)) if iso_ha[line] < 2 and iso_d4000[line] < dbreak(iso_mass[line])]
    K_A_d2_iso = [iso_d4000[line] for line in range(len(iso_d4000)) if iso_ha[line] > 2 and iso_d4000[line] > dbreak(iso_mass[line])]

    K_A_d1_niso = [niso_d4000[line] for line in range(len(niso_d4000)) if niso_ha[line] < 2 and niso_d4000[line] < dbreak(niso_mass[line])]
    K_A_d2_niso = [niso_d4000[line] for line in range(len(niso_d4000)) if niso_ha[line] > 2 and niso_d4000[line] > dbreak(niso_mass[line])]


    #H_alpha
    quenched_ha_iso = [iso_ha[line] for line in range(len(iso_d4000)) if iso_d4000[line] > dbreak(iso_mass[line]) and iso_ha[line] < 2]
    unquen_ha_iso = [iso_ha[line] for line in range(len(iso_d4000)) if iso_d4000[line] < dbreak(iso_mass[line]) and iso_ha[line] > 2]

    quenched_ha_niso = [niso_ha[line] for line in range(len(niso_d4000)) if niso_d4000[line] > dbreak(niso_mass[line]) and niso_ha[line] < 2]
    unquen_ha_niso = [niso_ha[line] for line in range(len(niso_d4000)) if niso_d4000[line] < dbreak(niso_mass[line]) and niso_ha[line] > 2]


    #Meet one criterion but not the other
    K_A_ha1_iso = [iso_ha[line] for line in range(len(iso_d4000)) if iso_ha[line] < 2 and iso_d4000[line] < dbreak(iso_mass[line])]
    K_A_ha2_iso = [iso_ha[line] for line in range(len(iso_d4000)) if iso_ha[line] > 2 and iso_d4000[line] > dbreak(iso_mass[line])]

    K_A_ha1_niso = [niso_ha[line] for line in range(len(niso_d4000)) if niso_ha[line] < 2 and niso_d4000[line] < dbreak(niso_mass[line])]
    K_A_ha2_niso = [niso_ha[line] for line in range(len(niso_d4000)) if niso_ha[line] > 2 and niso_d4000[line] > dbreak(niso_mass[line])]


    #Mass
    quenched_mass_iso = [iso_mass[line] for line in range(len(iso_d4000)) if iso_d4000[line] > dbreak(iso_mass[line]) and iso_ha[line] < 2]
    unquen_mass_iso = [iso_mass[line] for line in range(len(iso_d4000)) if iso_d4000[line] < dbreak(iso_mass[line]) and iso_ha[line] > 2]

    quenched_mass_niso = [niso_mass[line] for line in range(len(niso_d4000)) if niso_d4000[line] > dbreak(niso_mass[line]) and niso_ha[line] < 2]
    unquen_mass_niso = [niso_mass[line] for line in range(len(niso_d4000)) if niso_d4000[line] < dbreak(niso_mass[line]) and niso_ha[line] > 2]


    #Meet one criterion but not the other
    K_A_m1_iso = [iso_mass[line] for line in range(len(iso_d4000)) if iso_ha[line] < 2 and iso_d4000[line] < dbreak(iso_mass[line])]
    K_A_m2_iso = [iso_mass[line] for line in range(len(iso_d4000)) if iso_ha[line] > 2 and iso_d4000[line] > dbreak(iso_mass[line])]

    K_A_m1_niso = [niso_mass[line] for line in range(len(niso_d4000)) if niso_ha[line] < 2 and niso_d4000[line] < dbreak(niso_mass[line])]
    K_A_m2_niso = [niso_mass[line] for line in range(len(niso_d4000)) if niso_ha[line] > 2 and niso_d4000[line] > dbreak(niso_mass[line])]

    sample_d4000_iso = [sample_d4000[line] for line in range(len(sample_d4000)) if sample_dist[line] > 1.5]
    sample_d4000_niso = [sample_d4000[line] for line in range(len(sample_d4000)) if sample_dist[line] < 1.5]

    sample_mass_iso = [sample_mass[line] for line in range(len(sample_mass)) if sample_dist[line] > 1.5]
    sample_mass_niso = [sample_mass[line] for line in range(len(sample_mass)) if sample_dist[line] < 1.5]

    #Plot mass vs D4000, separated as isolated or not
    k = plt.figure()
    ax = k.add_subplot(121)

    plt.title('Isolated ($d_{host} > 1.5$ Mpc)')
    plt.xlabel('$\log_{10}[M_{stellar}$ $(M_{\odot})]$')
    plt.ylabel('$D_n4000$')
    ax.set_ylim([0.4, 2.2])
    ax.set_xlim([10, 6.95])
    plt.plot(quenched_mass_iso, quenched_d4000_iso, '.', color = 'magenta', markersize = 1, label = 'Quenched')
    plt.plot(K_A_m1_iso, K_A_d1_iso, '.', color = 'cyan', markersize = 1, label = 'K+A')
    plt.plot(K_A_m2_iso, K_A_d2_iso, '.', color = 'cyan', markersize = 1)
    plt.plot(unquen_mass_iso, unquen_d4000_iso, '.', color = 'black', markersize = 1, label = 'Not Quenched')
    plt.plot(sample_mass_iso, sample_d4000_iso, 'o', color = 'red', markersize = 5, label = 'This Sample')
    ax.legend(loc = 1, markerscale = 1, fontsize = 11)

    bx = k.add_subplot(122)
    plt.title('Social ($d_{host} < 1.5$ Mpc)')
    plt.xlabel('$\log_{10}[M_{stellar}(M_{\odot})]$')
    plt.ylabel('$D_n4000$')
    bx.set_ylim([0.4, 2.2])
    bx.set_xlim([10, 6.95])
    plt.plot(quenched_mass_niso, quenched_d4000_niso, '.', color = 'magenta', markersize = 1, label = 'Quenched')
    plt.plot(K_A_m1_niso, K_A_d1_niso, '.', color = 'cyan', markersize = 1, label = 'K+A')
    plt.plot(K_A_m2_niso, K_A_d2_niso, '.', color = 'cyan', markersize = 1)
    plt.plot(unquen_mass_niso, unquen_d4000_niso, '.', color = 'black', markersize = 1, label = 'Not Quenched')
    plt.plot(sample_mass_niso, sample_d4000_niso, 'o', color = 'red', markersize = 5, label = 'This Sample')

    plt.savefig('Calibs/pics/Geha_2012_iso_vs_niso')
    
    #BIN ISOLATED STUFF ACCORDING TO MASS
    mass_bins = np.arange(7, 10.5, 0.5)
    plotbins = np.arange(7.5, 10.5, 0.5)
    binned_isomass = np.histogram(iso_mass, mass_bins)[0]
    
    binned_iso_quenched = np.histogram(quenched_mass_iso, mass_bins)[0]
    iso_quen_frac = [binned_iso_quenched[line]/binned_isomass[line] for line in range(len(binned_isomass))]
    
    binned_iso_unquen = np.histogram(unquen_mass_iso, mass_bins)[0]
    iso_unquen_frac = [binned_iso_unquen[line]/binned_isomass[line] for line in range(len(binned_isomass))]
    
    
    binned_iso_KA1 = np.histogram(K_A_m1_iso, mass_bins)[0]
    iso_KA1_frac = [binned_iso_KA1[line]/binned_isomass[line] for line in range(len(binned_isomass))]
    
    
    binned_iso_KA2 = np.histogram(K_A_m2_iso, mass_bins)[0]
    iso_KA2_frac = [binned_iso_KA2[line]/binned_isomass[line] for line in range(len(binned_isomass))]
    
    print binned_isomass

    l = plt.figure()
    ax = l.add_subplot(211)
    
    plt.plot(plotbins, iso_quen_frac, 'o-', color = 'magenta', label = '$H\\alpha<2$, $D_n4000 >$ break')
    plt.plot(plotbins, iso_KA1_frac, 'o-', color = 'cyan', label = '$H\\alpha<2$, $D_n4000 <$ break')
    plt.plot(plotbins, iso_KA2_frac, 'o-', color = 'yellow', label = '$H\\alpha>2$, $D_n4000 >$ break')
    plt.title('Quenched / K+A fraction')
    plt.xlabel('Stellar Mass Bin: $\log_{10}[M_{stellar}$ $(M_{\odot})]$')
    plt.ylabel('Fraction of total $M_{stellar}$')
    ax.set_xlim([10.1, 7.4])
    plt.legend(loc = 1, fontsize = 11)
    
    bx = l.add_subplot(212)
    plt.plot(plotbins, iso_unquen_frac, 'o-', color = 'black', label = '$H\\alpha>2$, $D_n4000$ < break')
    plt.title('Star Forming Fraction')
    plt.xlabel('Stellar Mass Bin: $\log_{10}[M_{stellar}$ $(M_{\odot})]$')
    plt.ylabel('Fraction of total $M_{stellar}$')
    bx.set_xlim([10.1, 7.4])
    plt.legend(loc = 1, fontsize = 11)
    plt.tight_layout()
    plt.savefig('Calibs/pics/Geha_2012_quenched_fraction_binned.png')
    
    small_ha_high_d4000 = [iso_ha[line] for line in range(len(iso_ha)) if iso_mass[line] < 10.0 and iso_d4000[line] > dbreak(iso_mass[line])]
    high_d4000_small_ha = [iso_d4000[line] for line in range(len(iso_d4000)) if iso_mass[line] < 10.0 and iso_ha[line] < 2.0]
    
    k = plt.figure()
    ax = k.add_subplot(211)
    
    plt.hist(small_ha_high_d4000, range=(-5,20), bins = 20)
    plt.title('$M_{stellar}<10$ and $D_n4000 >$ break')
    plt.xlabel('H$\\alpha$EW')
    plt.ylabel('frequency')
    
    bx = k.add_subplot(212)
    plt.hist(high_d4000_small_ha, range = (0.7, 2.5), bins = 20)
    plt.title('$M_{stellar}<10$ and H$\\alpha<2$')
    plt.xlabel('$D_n4000$')
    plt.ylabel('frequency')
    plt.tight_layout()
    plt.savefig('Calibs/pics/Geha_2012_hist_of_high_low_dn_Ha')
    