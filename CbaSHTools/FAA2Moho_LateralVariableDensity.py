#!/usr/bin/env python

"""
Invert for crustal thickness map
"""
from __future__ import absolute_import, division, print_function

import os
import sys

sys.path.append(r'/Users/dingmin/tutorial_env/lib/python3.7/site-packages') 

import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
sys.path.append(os.path.join(os.path.dirname(__file__), "../../.."))

import pyshtools
from pyshtools import shio
from pyshtools import expand
from pyshtools import gravmag
from pyshtools import constant

pyshtools.utils.figstyle()


# ==== MAIN FUNCTION ====

def main():
    InvMoho()

# ==== TEST FUNCTIONS ====


def InvMoho():
    """
    Invert for the Moho relief from FAA map
    """
    delta_max = 5.0 # Default = 5 [m] for convergence criteria
    nmax = 10 #7
    degmax = 600 # [310]
    lmax = 1200 # [930]2*degmax ~ nmax*degmax
    #rho_c = 2550.0 # crustal density for use in the first iteration
    rho_m = 3220.0
    filter_type = 1
    half = 80 #80;  max=110
    #    Tc_mean = 44e3 # [m] assumed mean crustal thickness
    Tc_mean = 35e3 # [m] assumed mean crustal thickness
    
    # load gravity file
    #gravfile = '../../ExampleDataFiles/gmm3_120_sha.tab'
    #gravfile = 'JGGRAIL_1200C12A_SHA.TAB'
    gravfile = 'gggrx_1200a_sha.tab'
    pot, lmaxp, header = shio.shread(gravfile, lmax=degmax, header=True)
    gm = float(header[1]) * 1.e9
    mass = gm / constant.G.value
    r_grav = float(header[0]) * 1.e3
    print(r_grav, gm, mass, lmaxp)
    
    # load topography file
    #topofile = '../../ExampleDataFiles/MarsTopo719.shape'
    topofile = 'lro_ltm05_2050_sha.tab'
    hlm, lmaxt = shio.shread(topofile, lmax=degmax)
    r0 = hlm[0, 0, 0] # degree-0 radius (reference)
    d = r0 - Tc_mean
    print(r0, lmaxt)
    
    # load crusatal denisty variation grid (from the current folder)
    rho_grid = np.loadtxt('rhoc.in')
    
    # first expand the density grid to spherical harmonics
    rho = expand.SHExpandDH(rho_grid, sampling=2, lmax_calc=degmax)
#    # consider varible crustal density with resolution better than 5 degree (l = 36)
#    lrhomax = 36
#    rho[:,lrhomax+1:,:]=0 
#    # back to grid file
#    rho_grid = expand.MakeGridDH(rho, lmax=lmax, sampling=2,lmax_calc=lrhomax)
    
#    fig_map = plt.figure()
#    im = plt.imshow(rho_grid,cmap='jet')
#    fig_map.colorbar(im, orientation='horizontal')
#    fig_map.savefig('rhoc_tmp.png')
    
    # change the reference of the gravity file from r_grav to r0
    for l in range(2, lmaxp + 1):
        pot[:, l, :l + 1] = pot[:, l, :l + 1] * (r_grav / r0)**l

    # expand the toopgraphy
    topo_grid = expand.MakeGridDH(hlm, lmax=lmax, sampling=2,lmax_calc=degmax)
    print("Maximum radius (km) = ", topo_grid.max() / 1.e3)
    print("Minimum radius (km) = ", topo_grid.min() / 1.e3)


    # calculate the Bouguer correction
    #bc, r0 = gravmag.CilmPlusDH(topo_grid, nmax, mass, rho_c, lmax=degmax)
    bc, r0 = gravmag.CilmPlusRhoHDH(topo_grid, nmax, mass, rho_grid, lmax=degmax)
    
    # calculate the gravitational contribution from the lateral variations in density 
    # of the shell between the two radii r0 and d, referenced to r0
    # calculate the gravity terms (from Wieczoerk 2003 SM Equation 19)
    shellc = np.zeros([2, degmax + 1, degmax + 1], dtype=float)
    for l in range(1, degmax + 1):
        shellc[:, l, :l + 1] = rho[:, l, :l + 1]* 4 * np.pi * r0**3 / mass / \
                                    (2 * l + 1) / (l + 3) * (1 - (d / r0)**(l + 3))
    
#    shellc_grid,_, _, _, _ = gravmag.MakeGravGridDH(shellc, gm, r0, lmax=degmax)
#    shellc_grid = -shellc_grid*1e5
#
#    shellc1 = np.load('Clm_rho.npy')
#    shellc = np.zeros([2, degmax + 1, degmax + 1], dtype=float)
#    shellc[:,0:311,0:311] = shellc1
#
#    shellc1_grid,_, _, _, _ = gravmag.MakeGravGridDH(shellc, gm, r0, lmax=degmax)
#    shellc1_grid = -shellc1_grid*1e5
#
#
#    fig, (ax1,ax2) = plt.subplots(2,1)
#    pos = ax1.imshow(shellc_grid,cmap='jet') 
#    cbar = fig.colorbar(pos, ax=ax1,orientation='vertical')
#
#    pos = ax2.imshow(shellc1_grid,cmap='jet')
#    cbar = fig.colorbar(pos,ax=ax2,orientation='vertical')
#
#    fig.savefig('DensityShell.png')
#
#    
#
#    print('Maximum misfit = ', abs(shellc_grid-shellc1_grid).max())


##    """
##    Forward modeling of gravity attraction due to mare loading 
##    Pasted from ForwardMare_VariableDrhom.py
##    
##    """
#    # load mare bottom topography from the current folder
#    marethick_grid = np.loadtxt('marethick.in')
#    
#    # load crusatal for shallow part
#    rho_shallow_grid = np.loadtxt('rhoc_shallow.in')
#    
#
#    # calculate the gravity attaction due to surface with variable density
#    gt, r1 = gravmag.CilmPlusRhoHDH(topo_grid, nmax, mass, 3150-rho_shallow_grid, lmax=degmax)
#
#    # calculate the gravity attaction due to mare bottom topography with variable density
#    gb, r2 = gravmag.CilmPlusRhoHDH(topo_grid-marethick_grid, nmax, mass, 3150-rho_shallow_grid, lmax=degmax)
#
#    # calculate the overall contribution  
#    grav = gt - gb
#    
#    grav_grid,_, _, _, _ = gravmag.MakeGravGridDH(grav, gm, r0, lmax=degmax)
#    grav_grid = -grav_grid*1e5 # convert unit from m/s^2 to mGal; convert upward positive to downward positive
#
##    fig_map = plt.figure()
##    im = plt.imshow(grav_grid,cmap='jet')
##    plt.clim(0, 100)
##    fig_map.colorbar(im, orientation='horizontal')
##    fig_map.savefig('grav_centralmare.png')
##    np.savetxt('grav_centralmare.out',grav_grid)

## save BA coefficients before correcting for mare 
#    ba_tmp = pot - bc 
#    for l in range(2, lmaxp + 1):
#        ba_tmp[:, l, :l + 1] = ba_tmp[:, l, :l + 1] * (r0 / r_grav)**l
#    
#    Cba_grid = expand.MakeGridDH(ba_tmp, lmax=degmax, sampling=2,lmax_calc=degmax)
#    np.savetxt('BAcoef.out',Cba_grid)
#
## save BA coefficients after correcting for mare 
#    ba_tmp = pot - bc - 0.3*grav
#    for l in range(2, lmaxp + 1):
#        ba_tmp[:, l, :l + 1] = ba_tmp[:, l, :l + 1] * (r0 / r_grav)**l
#    
#    Cba_grid = expand.MakeGridDH(ba_tmp, lmax=degmax, sampling=2,lmax_calc=degmax)
#    np.savetxt('BAcoef_MareCorrected.out',Cba_grid)

#    # calculate the BA after subtracting bc & bs 
#    ba = pot - bc - shellc - 0.3*grav
#
#
#    """
#    End of mare loading modeling
#    """
    ba = pot - bc - shellc



    ba_grid,_, _, _, _ = gravmag.MakeGravGridDH(ba, gm, r0, lmax=degmax)
    ba_grid = -ba_grid*1e5
#    np.savetxt('BA.out',ba_grid)


#    """
#    Plot the three gravity components
#    """
#    
#    bc_grid,_, _, _, _ = gravmag.MakeGravGridDH(bc, gm, r0, lmax=degmax)
#    bc_grid = -bc_grid*1e5
#
#    shellc_grid,_, _, _, _ = gravmag.MakeGravGridDH(shellc, gm, r0, lmax=degmax)
#    shellc_grid = -shellc_grid*1e5
#
#    ba_grid,_, _, _, _ = gravmag.MakeGravGridDH(ba, gm, r0, lmax=degmax)
#    ba_grid = -ba_grid*1e5
#
#    pot_grid,_, _, _, _ = gravmag.MakeGravGridDH(pot, gm, r0, lmax=degmax)
#    pot_grid = -pot_grid*1e5
#
#    fig, (ax1,ax2,ax3,ax4) = plt.subplots(4,1)
#    pos = ax1.imshow(bc_grid,cmap='jet',vmin=-500,vmax=500)
#    cbar = fig.colorbar(pos, ax=ax1,orientation='vertical')
#
#    pos = ax2.imshow(shellc_grid,cmap='jet',vmin=-500,vmax=500)
#    cbar = fig.colorbar(pos,ax=ax2,orientation='vertical')
#
#    pos = ax3.imshow(ba_grid,cmap='jet',vmin=-500,vmax=500)
#    cbar = fig.colorbar(pos,ax=ax3,orientation='vertical')
#
#    pos = ax4.imshow(pot_grid,cmap='jet',vmin=-500,vmax=500)
#    cbar = fig.colorbar(pos,ax=ax4,orientation='vertical')
#
#    fig.savefig('GravComponets.png')
#
#
#    """
#    end plot
#    """

    moho_c = np.zeros([2, degmax + 1, degmax + 1], dtype=float)
    
    # initial h*drho, degree 1 to degmax (h is not moho, moho = d+h)
    for l in range(1, degmax + 1):
        if filter_type == 0:
            moho_c[:, l, :l + 1] = ba[:, l, :l + 1] * mass * (2 * l + 1) * \
                                   ((r0 / d)**l) \
                                   / (4.0 * np.pi * d**2)
        elif filter_type == 1:
            moho_c[:, l, :l + 1] = gravmag.DownContFilterMA(l, half, r0, d) * \
                                   ba[:, l, :l + 1] * mass * (2 * l + 1) * \
                                   ((r0 / d)**l) / \
                                   (4.0 * np.pi * d**2)
        else:
            moho_c[:, l, :l + 1] = gravmag.DownContFilterMC(l, half, r0, d) * \
                                   ba[:, l, :l + 1] * mass * (2 * l + 1) *\
                                   ((r0 / d)**l) / \
                                   (4.0 * np.pi * d**2)

    # expand to grid of h, divide by density map to derive initial moho relief  
    moho_grid3 = expand.MakeGridDH(moho_c, lmax=lmax, sampling=2,lmax_calc=degmax)
    moho_grid3 = np.divide(moho_grid3,rho_m-rho_grid)
    
    # add degree-0 term to h, which gives moho relief (radius)
    moho_c[0, 0, 0] = d
    moho_grid3 = moho_grid3 + d

    print('Maximum Crustal thickness (km) = ',(topo_grid - moho_grid3).max() / 1.e3)
    print('Minimum Crustal thickness (km) = ',(topo_grid - moho_grid3).min() / 1.e3)
    print('Average Crustal thickness (km) = ',(topo_grid - moho_grid3).mean() / 1.e3)


#    moho_c = gravmag.BAtoHilmDH(ba, moho_grid3, nmax, mass, r0,
#                                (rho_m - rho_c), lmax=lmax,
#                                filter_type=filter_type, filter_deg=half,
#                                lmax_calc=degmax)

    moho_c = gravmag.BAtoHilmRhoHDH(ba, moho_grid3, rho_m-rho_grid, nmax, mass, r0, 
                                filter_type=filter_type, filter_deg=half, 
                                lmax=lmax,lmax_calc=degmax)

    moho_grid2 = expand.MakeGridDH(moho_c, lmax=lmax, sampling=2,lmax_calc=degmax)
    print('Delta (km) = ', abs(moho_grid3 - moho_grid2).max() / 1.e3)

    temp_grid = topo_grid - moho_grid2
    print('Maximum Crustal thickness (km) = ', temp_grid.max() / 1.e3)
    print('Minimum Crustal thickness (km) = ', temp_grid.min() / 1.e3)
    print('Average Crustal thickness (km) = ', temp_grid.mean() / 1.e3)

    iter = 0
    delta = 1.0e9

    while delta > delta_max:
        iter += 1
        print('Iteration ', iter)

        moho_grid = (moho_grid2 + moho_grid3) / 2.0
        print("Delta (km) = ", abs(moho_grid - moho_grid2).max() / 1.e3)

        temp_grid = topo_grid - moho_grid
        print('Maximum Crustal thickness (km) = ', temp_grid.max() / 1.e3)
        print('Minimum Crustal thickness (km) = ', temp_grid.min() / 1.e3)
        print('Average Crustal thickness (km) = ', temp_grid.mean() / 1.e3)
        
        moho_grid3 = moho_grid2
        moho_grid2 = moho_grid

        iter += 1
        print('Iteration ', iter)

#        moho_c = gravmag.BAtoHilmDH(ba, moho_grid2, nmax, mass, r0,
#                                    rho_m - rho_c, lmax=lmax,
#                                    filter_type=filter_type, filter_deg=half,
#                                    lmax_calc=degmax)
        moho_c = gravmag.BAtoHilmRhoHDH(ba, moho_grid2, rho_m-rho_grid, nmax, mass, r0, 
                                    filter_type=filter_type, filter_deg=half, 
                                    lmax=lmax,lmax_calc=degmax)                                    

        moho_grid = expand.MakeGridDH(moho_c, lmax=lmax, sampling=2,
                                      lmax_calc=degmax)
                                      
        delta = abs(moho_grid - moho_grid2).max()
        print('Delta (km) = ', delta / 1.e3)

        temp_grid = topo_grid - moho_grid
        print('Maximum Crustal thickness (km) = ', temp_grid.max() / 1.e3)
        print('Minimum Crustal thickness (km) = ', temp_grid.min() / 1.e3)
        print('Average Crustal thickness (km) = ', temp_grid.mean() / 1.e3)

        moho_grid3 = moho_grid2
        moho_grid2 = moho_grid

        if temp_grid.max() > 100.e3:
            print('Not converging')
            exit(1)
    
    temp_grid = temp_grid*1e-3
#    fig_map = plt.figure()
#    im = plt.imshow(temp_grid,cmap='jet')
#    fig_map.colorbar(im, orientation='horizontal')
#    fig_map.savefig('InvCrustalThickness_test.png')
    np.savetxt('CrustalThickness_Lambda80.out',temp_grid)

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
