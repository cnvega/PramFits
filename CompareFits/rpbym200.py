#! /usr/bin/env python

import h5py
import matplotlib as mpl

mpl.use('PS')

import matplotlib.pyplot as plt
import numpy as np


def set_style(style='book', Hratio=1.0, Wfrac=1.0):
   if style == 'talk':
      size, fsize = 6, 16
   if style == 'book':
      size, fsize = 5.39, 12
   if style == 'mnras':
      size, fsize = 3.32, 8
   if style == 'mnras-fw':
      size, fsize = 6.97, 8

   mpl.rcParams['figure.figsize'] = (size*Wfrac, (size*Wfrac)*Hratio)
   mpl.rcParams['font.family'] = 'serif'
   mpl.rcParams['font.serif'] = ['Times', 'Liberation Serif', 'Times New Roman']
   mpl.rcParams['font.size'] = fsize
   mpl.rcParams['legend.fontsize'] = 'medium'
   mpl.rcParams['legend.frameon'] = False
   mpl.rcParams['text.usetex'] = True
   mpl.rcParams['axes.linewidth'] = 1.0
   try:
      mpl.rcParams['xtick.minor.visible'] = True
      mpl.rcParams['ytick.minor.visible'] = True
   except: pass
   try:
      mpl.rcParams['xtick.top'] = True
      mpl.rcParams['ytick.right'] = True
   except: pass
   try:
      mpl.rcParams['xtick.direction'] = 'in'
      mpl.rcParams['ytick.direction'] = 'in'
   except: pass

set_style(style='mnras', Hratio=1., Wfrac=1.)

f = h5py.File("rampressure_092.h5", "r")

rampress = f['rp'][:]
rpfit_tecce = f['rp_fit'][:]
rpfit_vega  = f['rp_fit_new'][:]
M200 = f['M200_host'][:]
R200 = f['R200_host'][:]
rrel = f['rrel'][:]

x = rrel/R200

massbin = [16,15,14,13,12]

fig, axs = plt.subplots(ncols=2, nrows=2, sharex='col', sharey='row')
for i1 in range(2): 
   for i2 in range(2):
      i = 2*i1+i2
      mmax = massbin[i]
      mmin = massbin[i+1]
      flt = (10**mmin <= M200)&(M200 < 10**mmax)
      axs[i1][i2].plot(x, np.log10(rampress), ".", ms=2, color="#D0D0FF")
      axs[i1][i2].plot(x[flt], np.log10(rampress[flt]), ".",color="purple", ms=1)
      axs[i1][i2].set_xlim((0,1.0))
      axs[i1][i2].set_ylim((-16,-8))

      if 0 != i:    
         label = r"$"+str(mmin)+"\leq \log M_{200} < "+str(mmax)+"$"
      else:
         label = r"$"+str(mmin)+"\leq \log M_{200} $"

      axs[i1][i2].text(0.1, 0.05, label, transform=axs[i1][i2].transAxes)

      if 1 == i1:
          axs[i1][i2].set_xlabel(r'$r/R_{200}$')
      if 0 == i2:
          axs[i1][i2].set_ylabel(r'$\log ({\rm P}_{\rm ram} [h^2 {\rm dyn\;cm}^{-2}])$')

fig.tight_layout()
fig.subplots_adjust(hspace=0.1, wspace=0.2)
fig.savefig("figs/rpbymass.eps")

