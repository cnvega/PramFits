#!/usr/bin/env python 

import matplotlib as mpl

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


masscut = 1e8

Nsnaps = rd.redshifts.shape[0]

figs = './figs_paper/'


def select_galaxies_1model():
   
   SAG_sim, Groups, ParentGroup, Contam = rd.read_data('SAGBeta1.3-NoTS-cal-RPsim', 91)
   GalTypes = SAG_sim.readDataset('Galaxy_Type')
   
   # Hacemos arrays de tamano Ngalx1 con las M200 y R200 de los main hosts:
   parentM200 = np.zeros(GalTypes.shape)
   parentR200 = np.zeros(GalTypes.shape)
   for i, par in enumerate(ParentGroup):
       parentM200[i] = Groups['GroupM200'][par]
       parentR200[i] = Groups['GroupR200'][par]
   # Unidades de SUBFIND: 1e10 Msun/h, 1.0 kpc/h, 1.0 km/s
   parentM200 *= 1e10    # in Msun/h
   
   # Seleccionamos las 3 galaxias centrales que viven en cumulos 
   centrals,_ = np.where((1e15<=parentM200)&
                       (0==GalTypes)&(False==Contam))
   GalIDs_sim = SAG_sim.readDataset('GalaxyStaticID')

   print("Selected galaxies: ", len(centrals))
   print("galaxies: ")
   print(GalIDs_sim[centrals])
   return GalIDs_sim[centrals]

def plot_satmassfunc():
   
   snaps = [91, 36, 26]
   redshifts = ['0.0', '1.5', '3.0']
  
   ncols = len(snaps)
   
   set_style('mnras-fw', Hratio=0.25)
   fig, axs = plt.subplots(nrows=1, ncols=ncols, sharey='row')

   cen_g15 = select_galaxies_1model()

   for snap, z, ax_z in zip(snaps, redshifts, range(len(snaps))):
      SAG_sim, Groups, ParentGroup, Contam = rd.read_data("SAGBeta1.3-NoTS-cal-RPsim", snap)
      units = SAG_sim.readUnits()
      
      SAGdata = SAG_sim   # to avoid modify to much code...

      galids = SAG_sim.readDataset('GalaxyStaticID')
      tipos = SAG_sim.readDataset('Galaxy_Type')
      discmass = SAG_sim.readDataset('M_star_disk')
      bulgemass = SAG_sim.readDataset('M_star_bulge')
      mstar_all = (discmass + bulgemass)*units.mass.Msun/units.h
      del discmass, bulgemass
      
      ax_g15 = axs[ax_z]
      
      # the main host halos:
      cenmask = np.in1d(galids, np.array(cen_g15))
      hosts = ParentGroup[cenmask]
      # Selecting the satellites:
      sats = np.in1d(ParentGroup, hosts)
      sats = sats.reshape((len(sats),1))
      # Discard the centrals:
      sats[tipos == 0] = False
      print("Selected "+str(len(sats))+" satellites at z = "+z)
      logmstar_g15 = np.log10(mstar_all[sats])
      del cenmask, hosts, sats
      
      phi_i, bins = np.histogram(logmstar_g15, bins=20, range=[6,13])
      phi = phi_i.astype(float)/(bins[1]-bins[0])/3.
      x = (bins[1:]+bins[:-1])/2.0
      ax_g15.plot(x[phi>0], np.log10(phi[phi>0]), '-', color='navy', lw=2, zorder=10)
      ax_g15.plot(np.log10([masscut]*2), [-1, 4], ':r', lw=1.5, zorder=1)
      del phi, phi_i, bins, x

      ax_g15.set_xlim((6,12.5))
      ax_g15.set_ylim((-0.1, 4)) 
      
      ax_g15.set_xlabel(r'$\log (M_\star[{\rm M}_\odot])$')
      ax_g15.text(0.05, 0.1, r'$z = '+z+'$', transform=ax_g15.transAxes)

   axs[0].set_ylabel(r'$\log ({\rm d N_{gal}/d}\log M_\star /{\rm N_{gr}})$')
   axs[0].text(0.48, 0.85, r"$\log (M_{200}/{\rm M}_\odot) > 15$", 
               transform=axs[0].transAxes)
   
   fig.tight_layout()
   fig.subplots_adjust(wspace=0.1)
   fig.savefig(figs+'SatMassFunc_m15.eps')


if __name__ == '__main__':
   
   plot_satmassfunc()



