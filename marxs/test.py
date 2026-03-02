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
marxm = marxs.optics.marx.MarxMirror(parfile="test_source.par")
print("Marx Mirror is ready.")
from marxs.optics import FlatDetector
det = optics.FlatDetector(position=[6710, 0,0 ], zoom=[1, 600, 600], pixsize=1)
optics.FlatDetector.display['color'] = 'green'
off =np.array([0.0,1.0,2.0,3.0,4.0,5.0, 6.0, 7.0])
ra = 30
dec = 30
pos = simulator.KeepCol('pos')
instrum = simulator.Sequence(elements=[marxm,det], postprocess_steps=[pos])
ncols=4
nrows=2
fig, ax = plt.subplots(nrows, ncols, figsize=(12,5), layout="constrained")

for i in range(len(off)):
    star = marxs.source.DiskSource(coords= SkyCoord(ra+off[i], dec+off[i], unit="deg"), a_outer = 0.5*u.degree, flux = 1e2 / u.s / u.cm**2, energy=3* u.keV )
    pointing = marxs.source.FixedPointing(coords=SkyCoord(ra, dec, unit='deg'))
    photons = star.generate_photons(1e2 * u.s)
    photons = pointing(photons)
    photons = instrum(photons)
    ind = (photons['unreflected'] == False)
    photons = photons[ind]
    #posdat = pos.format_positions()
    posdat = pos.format_positions()[ind,:,:]
    #print(photons)

    x = photons['detpix_x']
    y = photons['detpix_y']

    if (i<4):
        row = 0
        col = i
    else:
        row = 1 
        col = i-4
    print("Plotting for off-axis {} (in deg)".format(off[i]))
    ax[row, col].plot(x, y, ',')
plt.savefig("off_axis2.png")    
plt.show()
