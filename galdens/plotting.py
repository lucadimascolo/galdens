import matplotlib; matplotlib.use('TkAgg')
matplotlib.rc('font',**{'family':'serif','sans-serif':['Times'],'size':11})
matplotlib.rc('text',usetex=True)

from matplotlib import gridspec
import matplotlib.cm as cm
import matplotlib.colors as colors
import matplotlib.pyplot as plt
plt.rc('text.latex', preamble=r'\usepackage{newtxtext}')
plt.rc('text.latex', preamble=r'\usepackage{newtxmath}')
plt.rc('font',family='serif',size=11)
plt.rc('text',usetex=True)

from palettable.cartocolors.sequential import Mint_7 as pal4

cmap = pal4.mpl_colormap

import numpy as np
import scipy.special
import scipy.ndimage

from astropy.io import fits
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=70.00,Om0=0.30)

import aplpy
#######################################################################################################################
perc = {'1sigma': [0.50-0.50*scipy.special.erf(1.00/np.sqrt(2.00)),0.50+0.50*scipy.special.erf(1.00/np.sqrt(2.00))],
        '2sigma': [0.50-0.50*scipy.special.erf(2.00/np.sqrt(2.00)),0.50+0.50*scipy.special.erf(2.00/np.sqrt(2.00))],
        '3sigma': [0.50-0.50*scipy.special.erf(3.00/np.sqrt(2.00)),0.50+0.50*scipy.special.erf(3.00/np.sqrt(2.00))]}

#######################################################################################################################
dpi = 72.27*390.00/504.00

factorx = 0.40 # 0.60
factory = 0.40 # 0.50
figsize = (factorx*504.00/dpi,factory*504.00/dpi)
#######################################################################################################################
spifact = 1.00E+00

def plot_galdens(hdu,coord_array,spimin=None,spimax=None,spidlt=1.00):
    fig = plt.figure(figsize=figsize)

    spidata, spihead = hdu.data, hdu.header
    spidata *= spifact/np.abs(spihead['CDELT1']*60.00)/np.abs(spihead['CDELT2']*60.00)

    spiimg = aplpy.FITSFigure(data=fits.PrimaryHDU(data=spidata,header=spihead),figure=fig)
    spiimg.recenter(x=spihead['CRVAL1'],y=spihead['CRVAL2'],radius=2.40/60.00)
    
    if spimin is None: spimin = np.round(spidata.min()-0.50)
    if spimax is None: spimax = np.round(spidata.max()+0.50)

    spilev = np.linspace(spimin,spimax,1+int((spimax-spimin)/spidlt)) 

    spiimg.show_contour(data=fits.PrimaryHDU(data=spidata,header=spihead),cmap=cmap,filled=True,levels=spilev,extend='both')
    spicon = spiimg._layers['contour_set_1']
    spicon.cmap.set_under('white')

    spiimg.ax.tick_params(axis='both',which='both',direction='in')

    spibar = plt.colorbar(spicon)
    spibar.set_label(r'Galaxy density [$\mathrm{arcmin^{-2}}$]',labelpad=10)
    #spibar.ax.set_yticklabels(['{0:2d}'.format(int(np.round(float(tick.get_text()[1:-1])))) for tick in spibar.ax.yaxis.get_ticklabels()])

    spiimg.show_markers(coord_array[:,0],coord_array[:,1],coords_frame='world',marker='o',s=5,alpha=0.30,edgecolor='black',facecolor='none',linewidths=0.50)

    plt.show()
    #plt.savefig('figures/{0}_density_nopts_l15.pdf'.format(moo['name']),format='pdf',dpi=300)
    #plt.close()