from astropy.coordinates import SkyCoord
import astropy.units as u
from marxs import source, optics
from marxs.simulator import Sequence
from marxs import source, simulator
from marxs.optics import MultiAperture, CircleAperture,FlatStack
import numpy as np
import matplotlib.pyplot as plt
from marxs.optics.marx import MarxMirror
import os
import marxs

marxm = marxs.optics.marx.MarxMirror(parfile='test_source.par') #parameter file with the path to .rdb file 

from marxs.optics import FlatDetector
from marxs.design import uncertainties
#det = optics.FlatDetector(position=[6710, 0,0 ], zoom=[1, 150, 150], pixsize=0.0015)
det = optics.FlatDetector(position=[-310, 0, 0 ], zoom=[1, 250, 250], pixsize=0.00025)


off =np.array([0.0,0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7])
ra = 30
dec = 30
pos = simulator.KeepCol('pos')
instrum = simulator.Sequence(elements=[marxm,det], postprocess_steps=[pos])
ncols=4
nrows=2
fig, ax = plt.subplots(nrows, ncols, figsize=(28,10))

for i in range(len(off)):
    star = marxs.source.PointSource(coords= SkyCoord(ra+off[i], dec+off[i], unit="deg"), flux= 1e4 / u.s / u.cm**2, energy= 1.5* u.keV )
    pointing = marxs.source.FixedPointing(coords=SkyCoord(ra, dec, unit='deg'))
    photons = star.generate_photons(1e3 * u.s)
    photons = pointing(photons)
    photons = instrum(photons)
    #ind = (photons['unreflected'] == False)
    #photons = photons[ind]
    posdat = pos.format_positions()
    #posdat = pos.format_positions()[ind,:,:]
    print(photons)

    x = photons['det_x']
    y = photons['det_y']

    if (i<4):
        nrow = 0
        ncol = i
    else:
        nrow = 1 
        ncol = i-4
    print("Plotting for off-axis {} (in deg)".format(off[i]))

    xlim = x.min(), x.max()
    ylim = y.min(), y.max()

    hb = ax[nrow, ncol].hexbin(x, y, gridsize=1000, cmap='gnuplot2')
    ax[nrow, ncol].set(xlim=xlim, ylim=ylim)
    ax[nrow, ncol].set_title("Off-axis {}".format(off[i]))
    cb = fig.colorbar(hb, ax=ax[nrow, ncol], label='counts')

plt.savefig("xmm_offaxis_misalgn_all_ph.png")