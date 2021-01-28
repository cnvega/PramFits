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


models = ['SAGBeta1.3-NoTS-RPsim', 
          'SAGBeta1.3-NoTS',
          'SAGBeta1.3-NoTS-oldRP']

modelnames = ['SPH', 'This work', 'T11 profile']
formats = [':', '-', '--']
colors = ['navy', 'turquoise', 'orangered']

rd.Path_Clusters = "data/subfind" 
rd.Path_SAG = "data/galaxies"

figs = './figs_paper/'

snap = 91
#snap = 92

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

fmax_c,_ = np.where((1e15<=parentM200)&
                  (0==GalTypes)&(False==Contam))

fmid_c,_ = np.where((1e14<=parentM200)&(1e15>parentM200)&
                  (0==GalTypes)&(False==Contam))

fmin_c,_ = np.where((1e13<=parentM200)&(1e14>parentM200)&
                  (0==GalTypes)&(False==Contam))

fmax_s1,_ = np.where((1e15<=parentM200)&
                  (1==GalTypes)&(False==Contam))

fmid_s1,_ = np.where((1e14<=parentM200)&(1e15>parentM200)&
                  (1==GalTypes)&(False==Contam))

fmin_s1,_ = np.where((1e13<=parentM200)&(1e14>parentM200)&
                    (1==GalTypes)&(False==Contam))

#filtros = [fltmax, fltmid, fltmin]
### me quedo solo con 2 pq dan muy similares todos...

filtros_c = [fmax_c, fmid_c, fmin_c]
filtros_s1 = [fmax_s1, fmid_s1, fmin_s1]
fltlabels = [r'$\log (M_{200}) > 15$', r'$13 \leq \log (M_{200}) < 14$']

# Number of selected systems:

#Ngroups = [float(len(np.unique(ParentGroup[flt]))) for flt in filtros]

#print(Ngroups)

## Array de tiempos (ojo que combina todos los cumulos, asi que nos
## quedamos solo con el primero). Le cambio la shape para que no genere
## matrices al operar con arrays 1-D (N,):

dt = SAG_sim.readDataset('Histories/DeltaT_List')[:Nsnaps]
dt = dt.reshape((dt.shape[0], ))*units.time.yr  # yr/h

print("Selected galaxies: ", len(fmax_s1), len(fmid_s1), len(fmin_s1))

set_style('mnras', Hratio=1.3)
figh, axsh = plt.subplots(nrows=3, sharex=True)
figl, axsl = plt.subplots(nrows=3, sharex=True)

for sag, tag, axh, axl in zip(SAGdata, modelnames, axsh, axsl):

   disc = sag.readDataset('M_star_disk')*units.mass.Msun 
   bulge = sag.readDataset('M_star_bulge')*units.mass.Msun
   mstar = disc + bulge
   del disc, bulge

   coldgas_d = sag.readDataset('M_gas_disk')*units.mass.Msun
   coldgas_b = sag.readDataset('M_gas_bulge')*units.mass.Msun
   mgas = coldgas_d + coldgas_b
   del coldgas_d, coldgas_b
   
   mhot = sag.readDataset('M_hot')*units.mass.Msun

   galtype = sag.readDataset('Galaxy_Type')
   m200 = sag.readDataset('Halo/M200c')*units.mass.Msun

   fbar = (mstar + mgas + mhot)/m200

   del mgas, mhot, m200

   axh.scatter(np.log10(mstar[fmax_s1]), fbar[fmax_s1], s=3, color=(0.5,0.5,1.0))
   axh.plot(np.log10(mstar[fmax_c]), fbar[fmax_c], marker="^", ms=6, lw=0,
         mec="#FF0000", mfc="None")
   #print(fmax_c)
   #print(mstar[fmax_c])
   axh.plot([8,13], [0.1333]*2, '--k', lw=1)
   axh.set_ylabel(r"$f_\mathrm{bar} = (M_\star + M_\mathrm{gas})/M_\mathrm{tot}$")
   axh.set_ylim([0.,0.4])
   axh.set_xlim([8,13])

   
   axl.scatter(np.log10(mstar[fmin_s1]), fbar[fmin_s1], s=3, color=(0.5,0.5,1.0))
   axl.plot(np.log10(mstar[fmin_c]), fbar[fmin_c], marker="^", ms=6, lw=0,
         mec="#FF0000", mfc="None")
   axl.plot([8,13], [0.1333]*2, '--k', lw=1)
   axl.set_ylabel(r"$f_\mathrm{bar} = (M_\star + M_\mathrm{gas})/M_\mathrm{tot}$")
   axl.set_ylim([0.,0.4])
   axl.set_xlim([8,13])
   
   del fbar, mstar
      
axsh[-1].set_xlabel(r"$\log (M_\star)$")
axsl[-1].set_xlabel(r"$\log (M_\star)$")

figh.tight_layout()
figh.savefig(figs+'fbar_max.eps')

figl.tight_layout()
figl.savefig(figs+'fbar_min.eps')


## Just in case the referee ask for it, let's do a plot with the 
#  there mass samples of the SPH RP model.

filtros_c = [fmax_c, fmid_c, fmin_c]
filtros_s1 = [fmax_s1, fmid_s1, fmin_s1]
fltlabels = [r'$\log (M_{200}) > 15$',
             r'$14 \leq \log (M_{200}) < 15$',
             r'$13 \leq \log (M_{200}) < 14$']
fig, axs = plt.subplots(nrows=3, sharex=True)

sag = SAG_sim

disc = sag.readDataset('M_star_disk')*units.mass.Msun 
bulge = sag.readDataset('M_star_bulge')*units.mass.Msun
mstar = disc + bulge
del disc, bulge
coldgas_d = sag.readDataset('M_gas_disk')*units.mass.Msun
coldgas_b = sag.readDataset('M_gas_bulge')*units.mass.Msun
mgas = coldgas_d + coldgas_b
del coldgas_d, coldgas_b
mhot = sag.readDataset('M_hot')*units.mass.Msun
galtype = sag.readDataset('Galaxy_Type')
m200 = sag.readDataset('Halo/M200c')*units.mass.Msun

fbar = (mstar + mgas + mhot)/m200

del mgas, mhot, m200

for ax, flttag, fltc, flts in zip(axs, fltlabels, filtros_c, filtros_s1):

   ax.scatter(np.log10(mstar[flts]), fbar[flts], s=3, color=(0.5,0.5,1.0))
   ax.plot(np.log10(mstar[fltc]), fbar[fltc], marker="^", ms=4, lw=0,
         mec="#FF0000", mfc="None")
   ax.plot([8,13], [0.1333]*2, '--k', lw=1)
   ax.set_ylabel(r"$(M_\star + M_\mathrm{gas})/M_\mathrm{sub}$")
   ax.set_ylim([0.,0.38])
   ax.set_xlim([8,13])
   
   ax.text(0.05, 0.88, flttag, transform=ax.transAxes, fontsize='small')

del fbar, mstar
      
axs[-1].set_xlabel(r"$\log (M_\star \; h^{-1} [\mathrm{M}_\odot])$")

fig.tight_layout()
fig.subplots_adjust(hspace=0.05)
fig.savefig(figs+'fbar_SHP.eps')


