#! /usr/bin/env python

import h5py
import matplotlib as mpl

mpl.use('PS')

import matplotlib.pyplot as plt
import numpy as np

addprofiles = False

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

# I put this in a class only to encapsulate.
class RPFits():
   def __init__(self):
      # The units:
      self.L_cm = 3.085678e21
      self.V_cms = 1.0e5
      self.M_gr  = 1.989e43
      self.factor = self.L_cm**3/self.V_cms**2/self.M_gr
   
   def Pram_beta(self, r, p0, rs, beta):
       return p0*(1.+(r/rs)**2.)**(-3./2.*beta)
   
   # T11 fit:
   def tomas(self, m200, r200, aexp, r):
       ap = -0.8 + 1.2*(np.log10(m200) - 12.)
       bp = 1.2 - 0.4*(np.log10(m200) - 12.)
       ar = 0.59 - 0.14*(np.log10(m200) - 12.)
       br = -0.44 + 0.12*(np.log10(m200) - 12.)
       ab, bb = 0.92, -0.4
       beta = ab + bb*(aexp - 0.25)
       rs = r200*(ar + br*(aexp - 0.25))
       p0  = (1e-12)*10**(ap+bp*(aexp - 0.25))
       return self.Pram_beta(r, p0, rs, beta)
   
   # My new fit: 
   def cvm(self, m200, r200, aexp, r):
       rrel  = r/r200
       p0 = (1e-12)*10**(7.01*aexp**(-0.122) - 10 + 0.83)
       rs = -3.4*aexp**(-0.42) + 10.2
       bm  = 3.36e-3*aexp**(1.33) + 0.512
       bn  = -5.6 
       beta = bm*np.log10(m200) + bn
   
       return p0*(rrel/rs)**(-3./2.*beta)


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

profx = np.arange(0.1,2.2,0.05)
fits = RPFits()

fig, axs = plt.subplots(ncols=2, nrows=2, sharex='col', sharey='row')
for i1 in range(2): 
   for i2 in range(2):
      i = 2*i1+i2
      mmax = massbin[i]
      mmin = massbin[i+1]
      flt = (10**mmin <= M200)&(M200 < 10**mmax)
      axs[i1][i2].plot(x, np.log10(rampress), ".", ms=2, color="#E0E0FF")
      axs[i1][i2].plot(x[flt], np.log10(rampress[flt]), ".",color="purple", ms=1)
      
      if addprofiles:
         #median = np.median(np.unique(M200[flt]))
         mlow = np.min(np.unique(M200[flt]))
         mhigh = np.max(np.unique(M200[flt]))
         
         # The mass bins are too large for using just the medians...
         #axs[i1][i2].plot(profx, np.log10( fits.cvm(median, 1.0, 1.0, profx)  ), "-r")
         #axs[i1][i2].plot(profx, np.log10( fits.tomas(median, 1.0, 1.0, profx)  ), "-b")
         
         #axs[i1][i2].plot(profx, np.log10( fits.cvm(mlow, 1.0, 1.0, profx)  ), ":r")
         #axs[i1][i2].plot(profx, np.log10( fits.cvm(mhigh, 1.0, 1.0, profx)  ), ":r")
         axs[i1][i2].fill_between(profx, np.log10(fits.tomas(mhigh, 1.0, 1.0, profx)), 
               np.log10(fits.tomas(mlow, 1.0, 1.0, profx)), 
               color='red', alpha=0.3, zorder=9, label="T11 profile")
               #color='orangered', alpha=0.35, zorder=9)
         
         axs[i1][i2].fill_between(profx, np.log10(fits.cvm(mhigh, 1.0, 1.0, profx)), 
               np.log10(fits.cvm(mlow, 1.0, 1.0, profx)), 
               color='cyan', alpha=0.35, zorder=10, label="New profile")
               #color='turquoise', alpha=0.4, zorder=10)

      axs[i1][i2].set_xlim((0,1.8))
      axs[i1][i2].set_ylim((-16,-8))

      if 0 != i:    
         label = r"$"+str(mmin)+"\leq \log M_{200} < "+str(mmax)+"$"
      else:
         label = r"$"+str(mmin)+"\leq \log M_{200} $"

      axs[i1][i2].text(0.05, 0.07, label, transform=axs[i1][i2].transAxes)

      if 1 == i1:
          axs[i1][i2].set_xlabel(r'$r/R_{200}$')
          axs[i1][i2].set_xticks([0,0.5,1,1.5])
          if 1 == i2 and addprofiles:
             axs[i1][i2].legend(frameon=False, loc='upper right',
                  fontsize='small', handlelength=2, labelspacing=0.3)
      if 0 == i2:
          axs[i1][i2].set_ylabel(r'$\log ({\rm P}_{\rm ram} [h^2 {\rm dyn\;cm}^{-2}])$')

fig.tight_layout()
fig.subplots_adjust(hspace=0.1, wspace=0.1)

if not addprofiles:
   fig.savefig("figs/rpbymass.eps")
else:
   fig.savefig("figs/rpbymass_prof.pdf")
