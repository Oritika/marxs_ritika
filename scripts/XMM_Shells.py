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
det = optics.FlatDetector(position=[-310, 0, 0 ], zoom=[1, 150, 150], pixsize=0.0015)
off = 0.2

pos = simulator.KeepCol('pos')
star = marxs.source.PointSource(coords= SkyCoord(30, 30, unit="deg"), 
                               flux = 1e5 / u.s / u.cm**2, energy=1* u.keV )
pointing = marxs.source.FixedPointing(coords=SkyCoord(30+off, 30+off, unit='deg'))
photons = star.generate_photons(1e3* u.s)
instrum = simulator.Sequence(elements=[marxm, det], postprocess_steps=[pos])
photons = pointing(photons)
photons = instrum(photons)

ind = (photons['unreflected'] == False)
photons = photons[ind]
posdat = pos.format_positions()[ind,:,:]
#posdat = pos.format_positions()

x= photons['det_x']
y= photons['det_y']

xlim = x.min(), x.max()
ylim = y.min(), y.max()
fig, ax = plt.subplots()
hb = ax.hexbin(x, y, gridsize=3000, bins= 'log', cmap='gnuplot2')
ax.set(xlim=xlim, ylim=ylim)
ax.set_title("Off-axis {}".format(off))
cb = fig.colorbar(hb, ax=ax, label='counts')
plt.savefig("xmmm_test05.png")
