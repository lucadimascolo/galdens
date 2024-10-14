from astropy import constants as const
from astropy import units as u

from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.convolution import convolve_fft
from astropy.convolution import Gaussian2DKernel

import numpy as np
import scipy.stats

from . import plotting

# WCS 2D mesh grid
# ----------------------------------------------------------------------
def generateWCSMesh(header,ext=False):
  headerWCS = header.copy()

  if ext:
    headerWCS['CRPIX1'] = headerWCS['CRPIX1']+int(headerWCS['NAXIS1']/2)
    headerWCS['CRPIX2'] = headerWCS['CRPIX2']+int(headerWCS['NAXIS2']/2)
    headerWCS['NAXIS1'] = headerWCS['NAXIS1']+2*int(headerWCS['NAXIS1']/2)
    headerWCS['NAXIS2'] = headerWCS['NAXIS2']+2*int(headerWCS['NAXIS2']/2)

  cdelt1 = headerWCS[list(filter(lambda x: x in headerWCS,['CDELT1','CD1_1']))[0]]
  cdelt2 = headerWCS[list(filter(lambda x: x in headerWCS,['CDELT2','CD2_2']))[0]]

  gridWCS = WCS(headerWCS)
  gridmx, gridmy = np.meshgrid(np.arange(headerWCS['NAXIS1']),np.arange(headerWCS['NAXIS2']))
  gridwx, gridwy = gridWCS.all_pix2world(gridmx,gridmy,0)
  
  if (np.abs(gridwx.max()-gridwx.min()-3.6e2)<np.abs(2.0*cdelt1)): 
    gridix = np.where(gridwx>headerWCS['CRVAL1']+cdelt1*(headerWCS['NAXIS1']-headerWCS['CRPIX1']+1)+3.6e2)
    gridwx[gridix] = gridwx[gridix]-3.6e2
  return np.array([gridwx,gridwy])


# Extract RA/Dec info from ds9 region file
# ----------------------------------------------------------------------
def readregion(name):
  f = open(name,'r'); f.readline()
  p = []
  for line in f:
    l = line.replace('(',',').replace(')',',').split(',')
    p.append(np.array([float(val) for val in l[1:3]]))
  f.close()
  return np.asarray(p)


# Round up naxis
# ----------------------------------------------------------------------
def roundnaxis(radius,cdelt):
  naxis = 2.00*radius/cdelt
  scale = 10**int(np.log10(naxis))
  return int(scale*np.round((naxis/scale)+0.5))


# Build a galaxy density map
# ----------------------------------------------------------------------
def getmap(xpos,ypos,cdelt=1.00*u.arcsec,method='convolve',**kwargs):
  xmin, xmax = xpos.min(), xpos.max()
  ymin, ymax = ypos.min(), ypos.max()

  xmid = np.median(xpos); xrad = np.maximum(np.abs(xmax-xmid),np.abs(xmin-xmid))
  ymid = np.median(ypos); yrad = np.maximum(np.abs(ymax-ymid),np.abs(ymin-ymid))

  cdelt = cdelt.to(u.deg).value

  cradi = np.maximum(xrad,yrad)*kwargs.get('rfac',1.00)
  naxis = roundnaxis(cradi,cdelt)

  chist = fits.PrimaryHDU(data=np.zeros((naxis,naxis)))
  chist.header['CRVAL1'], chist.header['CRPIX1'], chist.header['CDELT1'] = xmid, 1.00+0.50*naxis, -cdelt
  chist.header['CRVAL2'], chist.header['CRPIX2'], chist.header['CDELT2'] = ymid, 1.00+0.50*naxis,  cdelt

  chist.header['CTYPE1'], chist.header['CUNIT1'] = 'RA---SIN', 'deg'
  chist.header['CTYPE2'], chist.header['CUNIT2'] = 'DEC--SIN', 'deg'

  chist.header['RADESYS'] = 'ICRS'

  wcs = WCS(chist.header)
  
  if method=='histogram':
    sigma = kwargs.get('sigma',None)

    if sigma is None:
      raise ValueError('If using method="histogram", you must provide a sigma value.')
    else:
      sigma = sigma.to(u.deg).value

    coords = SkyCoord(ra=xpos*u.deg,dec=ypos*u.deg,frame='icrs')
    idx = wcs.world_to_array_index(coords)

    for i in range(len(idx[0])):
      chist.data[idx[0][i],idx[1][i]] += 1

    chist.data = convolve_fft(chist.data,Gaussian2DKernel(sigma/cdelt,mode='linear_interp'))
  elif method=='kde':
    xgrid, ygrid = generateWCSMesh(chist.header,ext=False)

    kernel = kwargs.get('kernel','gaussian')  
    bw = kwargs.get('bw_method','scott')

    if kernel=='gaussian':
      kde = scipy.stats.gaussian_kde(np.vstack([ypos,xpos]),bw_method=bw)
      out = kde(np.vstack([ygrid.ravel(),xgrid.ravel()]))
      out = out.reshape(xgrid.shape)
      chist.data = out.copy()
    else:
      ValueError('For the moment, the KDE method only supports a Gaussian kernel.')

  return chist
