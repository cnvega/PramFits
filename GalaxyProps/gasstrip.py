#!/usr/bin/env python 

import matplotlib as mpl
mpl.use('Agg')

import read_data as rd
import numpy as np
import matplotlib.pyplot as plt

def set_style(style='book', Hratio=1.0, Wfrac=1.0):
   if style == 'talk':
      size, fsize = 6, 16
   if style == 'book':
      size, fsize = 5.39, 12
   if style == 'mnras':
      size, fsize = 3.3, 8
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

##### Here it starts the script:

set_style('mnras', Hratio=1.0, Wfrac=1.0)

models = ['SAGBeta1.3-NoTS-RPsim', 
          'SAGBeta1.3-NoTS',
          'SAGBeta1.3-NoTS-oldRP']

modelnames = ['SPH', 'This work', 'T11 profile']
formats = [':', '-', '--']
colors = ['navy', 'turquoise', 'orangered']

figs = './figs_paper/'

snap = 91

Nsnaps = rd.redshifts.shape[0]


SAG_sim, Groups, ParentGroup, Contam = rd.read_data('SAGBeta1.3-NoTS-cal-RPsim', snap)
SAG_fit,      _,           _,      _ = rd.read_data('SAGBeta1.3-NoTS-cal', snap)
SAG_badfit,   _,           _,      _ = rd.read_data('SAGBeta1.3-NoTS-cal-oldRP', snap)

SAGdata = [SAG_sim, SAG_fit, SAG_badfit]

## 'SAG_xxx': para extraer datasets de los archivos SAG (Ngal x [1|snaps])
## 'Groups': diccionario con las propiedades de main host halos (Ngroups)
##           NOTE: arrays en unidades de SUBFIND!
## 'ParentGroup': Indice de los main hosts del SAG, en Groups (Ngal x 1)
## 'Contam': En 'True' las galaxias en grupos contaminados 

units = SAG_fit.readUnits()

GalTypes = SAG_sim.readDataset('Galaxy_Type')

# Hacemos arrays de tamano Ngalx1 con las M200 y R200 de los main hosts:

parentM200 = np.zeros(GalTypes.shape)
parentR200 = np.zeros(GalTypes.shape)
for i, par in enumerate(ParentGroup):
    parentM200[i] = Groups['GroupM200'][par]
    parentR200[i] = Groups['GroupR200'][par]

# Unidades de SUBFIND: 1e10 Msun/h, 1.0 kpc/h, 1.0 km/s
parentM200 *= 1e10    # in Msun/h

# Hacemos filtros para todas las galaxias tipo 1,2 que viven en cumulos 
# no contaminados, para diferentes rangos de masa log(M200)=(13, 14, 15)

fltmax,_ = np.where((1e15<=parentM200)&
                    (0<GalTypes)&(False==Contam))

fltmid,_ = np.where((1e14<=parentM200)&(1e15>parentM200)&
                    (0<GalTypes)&(False==Contam))

fltmin,_ = np.where((1e13<=parentM200)&(1e14>parentM200)&
                    (0<GalTypes)&(False==Contam))

#filtros = [fltmax, fltmid, fltmin]
### me quedo solo con 2 pq dan muy similares todos...

filtros = [fltmax, fltmin]
fltlabels = [r'$\log (M_{200}h^{-1}[\mathrm{M}_\odot]) > 15$', r'$13 \leq \log (M_{200}h^{-1}[\mathrm{M}_\odot]) < 14$']

# Number of selected systems:

Ngroups = [float(len(np.unique(ParentGroup[flt]))) for flt in filtros]

print(Ngroups)

## Vamos a calcular la media de masa de estos dos samples:
meds = []
flttmp,_ = np.where((1e15<=parentM200)&
                    (0==GalTypes)&(False==Contam))
meds.append(parentM200[flttmp].mean())
print(len(flttmp))
flttmp,_ = np.where((1e13<=parentM200)&(1e14>parentM200)&
                    (0==GalTypes)&(False==Contam))
meds.append(parentM200[flttmp].mean())
print(len(flttmp))
del flttmp
print("{:e} {:e}".format(meds[0],meds[1]))

## Array de tiempos (ojo que combina todos los cumulos, asi que nos
## quedamos solo con el primero). Le cambio la shape para que no genere
## matrices al operar con arrays 1-D (N,):

dt = SAG_sim.readDataset('Histories/DeltaT_List')[:Nsnaps]
dt = dt.reshape((dt.shape[0], ))*units.time.yr  # yr/h

print("Selected galaxies: ", len(fltmax), len(fltmid), len(fltmin))

fig, axs = plt.subplots(nrows=2, sharex=True)
figc, axsc = plt.subplots(nrows=2, sharex=True)

print('redshift.shape: ', rd.redshifts.shape)

for sag, tag, fmt, col in zip(SAGdata, modelnames, formats, colors):

   stripHG_all = sag.readDataset('Histories/HotGasRPejected')
   stripCG_all = sag.readDataset('Histories/ColdGasRPejected')
   stripGas_all = stripHG_all + stripCG_all
   del stripHG_all, stripCG_all
  
   rpsz0 = []

   for ax, axc, flt, flttag, ng in zip(axs, axsc, filtros, fltlabels, Ngroups):
      strip_sum = np.zeros(Nsnaps)
      for i in range(Nsnaps):
         strip_sum[i] = stripGas_all[flt, i].sum()

      ax.plot(rd.redshifts[25:], np.log10(strip_sum[25:]/dt[25:]/ng), ls=fmt, c=col, 
               lw=1.5, label=tag)
       
      ax.set_ylabel(r'$\log (\mathrm{d}M_\mathrm{RPS}/\mathrm{d}t/\mathrm{N_{gr}}[\mathrm{M_\odot\;yr^{-1}}])$')
      ax.set_xlim([0,3])

      ax.text(0.5, 0.1, flttag, transform= ax.transAxes, fontsize='small')

      strip_tot = np.zeros(Nsnaps)
      strip_tot[0] = strip_sum[0]
      for i in range(1, Nsnaps):
         strip_tot[i] = strip_sum[i] + strip_tot[i-1]

      rpsz0.append(strip_tot[-1])
  
      axc.plot(rd.redshifts[25:], np.log10(strip_tot[25:]/ng), ls=fmt, c=col, 
            lw=1.5, label=tag)
   
      axc.set_ylabel(r'$\log (M_{\rm RPS} (> z) h^{-1}/\mathrm{N_{gr}} [{\rm M}_\odot])$')
      axc.set_xlim([0,3])

      axc.text(0.08, 0.1, flttag, transform= axc.transAxes, fontsize='small')

   print("{:e}, {:e} = {:.2f}, {:.2f}".format(*rpsz0,
      100.*rpsz0[0]/meds[0]/Ngroups[0], 100.*rpsz0[1]/meds[1]/Ngroups[1]))

axs[1].set_xlabel('$z$')
axs[0].legend(frameon=False, loc='upper right',
      fontsize='small', handlelength=2, labelspacing=0.3)
fig.tight_layout()

fig.savefig(figs+'GasStripped.eps')


axsc[1].set_xlabel('$z$')
axsc[0].legend(frameon=False, fontsize='small', handlelength=2,
               labelspacing=0.3)
figc.tight_layout()
figc.savefig(figs+'GasStripped_acc.eps')

