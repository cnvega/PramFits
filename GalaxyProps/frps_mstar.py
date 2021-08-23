#!/usr/bin/env python 

import matplotlib as mpl
mpl.use('PS')

import read_data as rd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import colors as mplcolors

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

models = ['SAGBeta1.3-NoTS-cal-RPsim', 
          'SAGBeta1.3-NoTS-cal',
          'SAGBeta1.3-NoTS-cal-oldRP']

modelnames = ['SPH', 'New profile', 'T11 profile']
formats    = [':', '-', '--']
colors     = ['navy', 'turquoise', 'orangered']


masscut = 5e7
#masscut = 1e8

ScatterPlot = False

Nsnaps = rd.redshifts.shape[0]

figs = './figs_paper/'

##### Here it starts the script:

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

###############################################################################

def plot_frps_mstar(centrals):
  
   set_style('mnras', Hratio=0.618)
   fig, ax = plt.subplots() 

   snap = 91
   ztag = '0.0'
   
   # these are for the loop across models:
   colors     = ['navy', 'turquoise', 'orangered']

   SAG_sim, Groups, ParentGroup, Contam = rd.read_data('SAGBeta1.3-NoTS-cal-RPsim', snap)
   SAG_fit,      _,           _,      _ = rd.read_data('SAGBeta1.3-NoTS-cal', snap)
   SAG_badfit,   _,           _,      _ = rd.read_data('SAGBeta1.3-NoTS-cal-oldRP', snap)
   units = SAG_fit.readUnits()
   
   SAGdata = [SAG_sim, SAG_fit, SAG_badfit]
   
   # finding the central galaxies and their main host halos:
   galids = SAG_sim.readDataset('GalaxyStaticID')
   cenmask = np.in1d(galids, np.array(centrals))
   del galids
   hosts = ParentGroup[cenmask]
   
   # Selecting the satellites:
   tipos = SAG_sim.readDataset('Galaxy_Type')
   sats = np.in1d(ParentGroup, hosts)
   sats = sats.reshape((len(sats),1))
   # Discard the centrals:
   sats[tipos == 0] = False
   del tipos

   # Finally we apply a mass cut according to the CSMF behavior. These
   # selected galaxies will be the same in the three models:
   discmass = SAG_sim.readDataset('M_star_disk')
   bulgemass = SAG_sim.readDataset('M_star_bulge')
   mstar = (discmass + bulgemass)*units.mass.Msun/units.h
   sats[mstar < masscut] = False
   del discmass, bulgemass, mstar
   
   # We add a new filter to discard the galaxies with Mstar=0 in the 
   # other two models:
   discmass = SAG_fit.readDataset('M_star_disk')
   bulgemass = SAG_fit.readDataset('M_star_bulge')
   mstar = (discmass + bulgemass)*units.mass.Msun/units.h
   sats[mstar == 0] = False
   del discmass, bulgemass, mstar
   
   discmass = SAG_badfit.readDataset('M_star_disk')
   bulgemass = SAG_badfit.readDataset('M_star_bulge')
   mstar = (discmass + bulgemass)*units.mass.Msun/units.h
   sats[mstar == 0] = False
   del discmass, bulgemass, mstar

   ymax = []
   
   delta = 0.3
   mstarbins = np.arange(7.5, 12, delta)

   # Now we get the data and plot the different models: 
   for sag, tag, fmt, col in zip(SAGdata, modelnames, formats, colors):
      discmass = sag.readDataset('M_star_disk')
      bulgemass = sag.readDataset('M_star_bulge')
      mstar_all = (discmass + bulgemass)*units.mass.Msun/units.h
      del discmass, bulgemass
     
      # Select only the sample:
      mstar = mstar_all[sats]
      del mstar_all

      # Load the mass stripped by ram pressure:
      stripHG_all = sag.readDataset('Histories/HotGasRPejected')
      stripCG_all = sag.readDataset('Histories/ColdGasRPejected')
      sats_arr = sats.reshape((len(sats), ))
      stripGas = stripHG_all[sats_arr,:] + stripCG_all[sats_arr,:]
      del stripHG_all, stripCG_all

      # Get the fraction of total stripped mass w/r to the stellar mass
      strip_sum = np.zeros(mstar.shape)
      for gal in range(len(mstar)):
         strip_sum[gal] = stripGas[gal,:].sum()
      del stripGas 
      
      logmstar = np.log10(mstar)
      logfrps = np.log10(strip_sum/mstar)

      # A few statistics:
      print("log(Ms)>=8: ", logmstar[logmstar>=8].shape)
      print("log(Ms)>10.4: ", logmstar[logmstar>10.4].shape)
   
      frpsmed = np.zeros(len(mstarbins)-1)
      frps10 = np.zeros(len(mstarbins)-1)
      frps90 = np.zeros(len(mstarbins)-1)
      
      x = mstarbins[:-1]+delta/2.
      
      for i in range(len(mstarbins)-1):
         msmin, msmax = mstarbins[i], mstarbins[i+1]
     
         #print(msmin, msmax)
         flt = (msmin<=logmstar)&(logmstar<msmax)
         
         frpsmed[i] = np.median(logfrps[flt])
         frps10[i] = np.quantile(logfrps[flt], 0.1)
         frps90[i] = np.quantile(logfrps[flt], 0.9)

       
      # a scatter plot:
      if 'SPH'==tag:
         if ScatterPlot:
            ax.scatter(logmstar, logfrps, s=1, color='#E0E0E0')
         else:
            ax.fill_between(x, frps90, frps10, color='#f0f0ff')
      
      ax.plot(x, frpsmed, ls=fmt, c=col, lw=1.5, label=tag)
      ax.plot(x, frps10, ls=fmt, c=col, lw=0.5)
      ax.plot(x, frps90, ls=fmt, c=col, lw=0.5)

   ax.text(0.04, 0.08, r"$z="+ztag+"$", transform=ax.transAxes)
   
   # some clearing to save memory:
   del sats
   del SAG_sim, SAG_fit, SAG_badfit, Groups, ParentGroup, Contam
   
   ax.set_xlim((8, 11))
   #ax.set_xlim((7.5, 11))
   ax.set_ylim((-1.1, 1.6))
   #ax.set_ylim((-3, 1.6))
  
   ax.set_ylabel(r'$\log(f_{\rm RPS})$')
   ax.set_xlabel(r'$\log(M_\star [{\rm M}_\odot])$')
   
   ax.legend(frameon=False, loc="upper right",
             numpoints=1, fontsize='small', handlelength=2.0, labelspacing=0.3)
   
   fig.tight_layout()
   if ScatterPlot:
      fig.savefig(figs+"frps_mstar_z0_scatter.eps")
   else:
      fig.savefig(figs+"frps_mstar_z0.eps")



###############################################################################
###############################################################################
###############################################################################


if __name__ == '__main__':
    
   centrals = select_galaxies_1model()
   
   plot_frps_mstar(centrals)
   

