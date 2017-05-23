import numpy as np
import astropy.io.fits as fits
import apex_helper as ah
import input_parameters
from input_parameters import pathto, name, loc_fits

line = fits.open(loc_fits+"apec_v3.0.2_line.fits")
continuum = fits.open(loc_fits+"apec_v3.0.2_coco.fits")

for ii in range(len(name)):
    print name[ii]
    sim = ah.open_files(name=name[ii])
    eps = ah.lines(sim, line)
    eps_c = ah.continua(sim, continuum)
    Llambda = sim.Ox*eps[0]+sim.Fe*eps[1]+sim.H*eps[2]+sim.He*eps[3]
    Clambda = sim.Ox*eps_c[0]+sim.Fe*eps_c[1]+sim.H*eps_c[2]+sim.He*eps_c[3]
    lambdaTOT = Llambda+Clambda # line + continuum

    Lx = (lambdaTOT*sim.ne*sim.ni*sim.vol).sum()
    print lambdaTOT
    print "X--ray luminosity in the 0.5-2keV band  is %e erg s^-1 \n \n" %(Lx)

















