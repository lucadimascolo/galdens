import matplotlib; matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

from astropy import constants as const
from astropy import units as u

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian2DKernel

import numpy as np

import aplpy

# Extract RA/Dec info from ds9 region file
# ------------------------------------------------------
def readregion(name):
  f = open(name,'r'); f.readline()
  p = []
  for line in f:
    l = line.replace('(',',').replace(')',',').split(',')
    p.append(np.array([float(val) for val in l[1:3]]))
  f.close()
  return np.asarray(p)

# Round up naxis
# ------------------------------------------------------
def roundnaxis(radius,cdelt):
  naxis = 2.00*radius/cdelt
  scale = 10**int(np.log10(naxis))
  return int(scale*np.round((naxis/scale)+0.5))

# Build a galaxy density map
# ------------------------------------------------------
def getmap(xpos,ypos,sigma=1.00*u.arcmin,cdelt=1.00*u.arcsec):
  xmin, xmax = xpos.min(), xpos.max()
  ymin, ymax = ypos.min(), ypos.max()

  xmid = np.median(xpos); xrad = np.maximum(np.abs(xmax-xmid),np.abs(xmin-xmid))
  ymid = np.median(ypos); yrad = np.maximum(np.abs(ymax-ymid),np.abs(ymin-ymid))

  cdelt = cdelt.to(u.deg).value
  sigma = sigma.to(u.deg).value

  cradi = np.maximum(xrad,yrad)
  naxis = roundnaxis(cradi,cdelt)

  chist = fits.PrimaryHDU(data=np.zeros((naxis,naxis)))
  chist.header['CRVAL1'], chist.header['CRPIX1'], chist.header['CDELT1'] = xmid, 1.00+0.50*naxis, -cdelt
  chist.header['CRVAL2'], chist.header['CRPIX2'], chist.header['CDELT2'] = ymid, 1.00+0.50*naxis,  cdelt

  chist.header['CTYPE1'], chist.header['CUNIT1'] = 'RA---SIN', 'deg'
  chist.header['CTYPE2'], chist.header['CUNIT2'] = 'DEC--SIN', 'deg'

  chist.header['RADESYS'] = 'ICRS'

  wcs = WCS(chist.header)
  coords = SkyCoord(ra=xpos*u.deg,dec=ypos*u.deg,frame='icrs')
  idx = wcs.world_to_array_index(coords)

  for i in range(len(idx[0])):
    chist.data[idx[0][i],idx[1][i]] += 1

  chist.data = convolve_fft(chist.data,Gaussian2DKernel(sigma/cdelt,mode='linear_interp'))

  return chist
