import galdens
import numpy as np

from astropy import units as u

regfile = 'input/MOO_1142+1527.reg'
regdata = galdens.readregion(regfile)

sigma = 1*u.arcmin

hdu = galdens.getmap(regdata[:,0],regdata[:,1],sigma=sigma,cdelt=0.5*u.arcsec)