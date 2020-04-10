import pyfits as pf
from scipy import interpolate
import pandas as pd
import numpy as np
from scipy.signal import find_peaks
from sklearn import linear_model
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from QSmooth import open_calibrate_fits, mask_SDSS, smooth

# Adjust the filename and path here:
filename = 'spec-5488-56013-0860.fits'
path = 'example_fits/'

###################################

[loglam, flux, err] = open_calibrate_fits(filename=filename,path=path)
mask = mask_SDSS(filename=filename,path=path)
[loglam_smooth, flux_smooth] = smooth(loglam,flux,err,mask=mask)

# Plot the raw data and the QSmooth flux fit
plt.rc('text',usetex=True)
font = {'family':'sans-serif','sans-serif':['Helvetica'],'size':8}
plt.rc('font',**font)
fig, axs = plt.subplots(2,1,sharey=True,figsize=(6.97,3.31))
axs[0].plot(10**loglam,flux,c=(0.25,0.25,0.25),label=r'${\rm \ SDSS\ J151727.68+133358.60\ raw\ data}$')
axs[1].plot(10**loglam,flux,c=(0.25,0.25,0.25))
axs[0].plot(10**loglam,1/err,c='r',linewidth=1,label=r'${\rm \ flux\ errors}$')
axs[1].plot(10**loglam,1/err,c='r',linewidth=1)
axs[0].plot(10**loglam_smooth,flux_smooth,color='c',label=r'${\rm QSmooth\ fit}$')
axs[1].plot(10**loglam_smooth,flux_smooth,color='c')
axs[0].set_ylabel(r'${\rm flux\ [erg\ s}$'+r'$^{-1} {\rm cm} ^{-2}\mathrm{\AA}$'+r'${\rm ]}$')
axs[1].set_ylabel(r'${\rm flux\ [erg\ s}$'+r'$^{-1} {\rm cm} ^{-2}\mathrm{\AA}$'+r'${\rm ]}$')
axs[1].set_xlabel(r'${\rm rest\ wavelength\ [}$' r'$\mathrm{\AA}$' r'${\rm ]}$')
axs[0].set_xlim([1150, 2900])
axs[1].set_xlim([1150, 1500])
axs[0].set_yticks([0, 100])
axs[1].set_yticks([0, 100])
axs[0].set_ylim(bottom=-5)
axs[1].set_ylim(bottom=-5)
axs[0].legend(frameon=False)
plt.savefig('example_plots/SDSSJ151727.68+133358.60_example.png',bbox_inches='tight',dpi=400)
