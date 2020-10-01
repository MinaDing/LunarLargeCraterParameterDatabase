#!/usr/bin/env python
"""
This script forward models the contribution of mare basalts (only in the central r<0.2*R region) to the gravity anomaly 
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
    ForwardGrav()

# ==== TEST FUNCTIONS ====


def ForwardGrav():
    """
    forward modeling of gravity attraction due to post-impact mare deposits
    """
    nmax = 7
    degmax = 600
    lmax = 2*degmax # 2*degmax ~ nmax*degmax
    gm = 4.9035e12
    mass = gm / constant.G.value # [kg]
    # load topography file
#   topofile = '../../ExampleDataFiles/MarsTopo719.shape'
    topofile = '/Users/dingmin/Data/MoonShape/lro_ltm05_2050_sha.tab'
    hlm, lmaxt = shio.shread(topofile,lmax=degmax)
    r0 = hlm[0, 0, 0]
    #    hlm[:,0:16,:] = 0 # set deg smaller than 16 to be 0
    
    # expand the toopgraphy    
    topo_grid = expand.MakeGridDH(hlm, lmax=lmax, sampling=2,lmax_calc=degmax)
    
    # load \deltadenisty variation grid from the current folder
    drho_grid = np.loadtxt('drho_mare.in')

#    # transfer the density grid to spherical harmonics
#    rho = expand.SHExpandDH(drho_grid, sampling=2, lmax_calc=degmax)
    
    # load mare bottom topography from the current folder
    marethick_grid = np.loadtxt('marethick.in')
    marethick_grid = 0.5*marethick_grid
    
#    # transfer the mare bottom topography to spherical harmonics
#    hlm_bot = expand.SHExpandDH(marebot, sampling=2, lmax_calc=degmax)

    # calculate the gravity attaction due to surface with variable density
    gt, r0 = gravmag.CilmPlusRhoHDH(topo_grid, nmax, mass, drho_grid, lmax=degmax)

    # calculate the gravity attaction due to mare bottom topography with variable density
    gb, r0 = gravmag.CilmPlusRhoHDH(topo_grid-marethick_grid, nmax, mass, drho_grid, lmax=degmax)
    
    # calculate the overall contribution  
    grav = gt - gb

    # expand the grav
    Cgrav_centralmare = expand.MakeGridDH(grav, lmax=degmax, sampling=2,lmax_calc=degmax)
    np.savetxt('Cgrav_centralmare_50percent.out',Cgrav_centralmare)

#    grav_grid,theta, phi, total, pot = gravmag.MakeGravGridDH(grav, gm, r0, lmax=degmax)
#    grav_grid = -grav_grid*1e5 # convert unit from m/s^2 to mGal; convert upward positive to downward positive
#    
#    fig_map = plt.figure()
#    im = plt.imshow(grav_grid,cmap='jet')
#    plt.clim(0, 100)
#    fig_map.colorbar(im, orientation='horizontal')
#    fig_map.savefig('grav_centralmare.png')
#    np.savetxt('grav_centralmare.out',grav_grid)

# ==== EXECUTE SCRIPT ====
if __name__ == "__main__":
    main()
