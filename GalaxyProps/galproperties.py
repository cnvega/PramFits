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

PDFhists = True

models = ['SAGBeta1.3-NoTS-cal-RPsim', 
          'SAGBeta1.3-NoTS-cal',
          'SAGBeta1.3-NoTS-cal-oldRP']

modelnames = ['SPH', 'New profile', 'T11 profile']
formats    = [':', '-', '--']
colors     = ['navy', 'turquoise', 'orangered']

masscut = 1e8

Nsnaps = rd.redshifts.shape[0]

figs = './figs_paper/'

s2z = {91:'0.0', 36:'1.5', 26:'3.0'}
z2s = {'0.0':91, '1.5':36, '3.0':26}
   

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

def plot_frps_comb(centrals):
  
   set_style('mnras', Hratio=1.4)
   fig, axs = plt.subplots(nrows=3, sharex=True) 

   #snaps = [92, 36, 26]
   snaps = [91, 36, 26]
   redshifts = ['0.0', '1.5', '3.0']
   
   # these are for the loop across models:
   colors     = ['navy', 'turquoise', 'orangered']

   fout = open("frps_quantiles.txt", "w")
   fout.write("#z qnt  SPH  NewFit T11fit\n")

   for snap, ztag, ax in zip(snaps, redshifts, axs):

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
      dw = 0.0

      q50arr, q75arr = [], []

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
         logfrps = np.log10(strip_sum/mstar)
         
         # Now we calculate the distribution of this in the galaxies:
         xbins = np.arange(-3, 3, 0.35)
         phi_i, bins = np.histogram(logfrps, bins=xbins)
         x = (bins[1:]+bins[:-1])/2.0
         w = bins[1]-bins[0]
        
         # averaged number density per cluster:
         phi = phi_i.astype(float)/w/float(len(centrals))
            
         if 'SPH'==tag:   # in this case we add a shaded region
            ax.bar(x, phi, width=w-dw, ec=None, fc='#f0f0ff', fill=True)
         
         ax.step(x, phi, where='mid', ls=fmt, c=col, lw=1.5, label=tag)
         
         q50arr.append(np.quantile(logfrps, 0.50))
         q75arr.append(np.quantile(logfrps, 0.75))

         dw += 0.2*w
         ax.set_ylabel(r'${\rm dN_{gal} / d} \log (f_{\rm RPS}) / {\rm N_{gr}}$')
         ymax.append(phi.max())
      
      fout.write("{:s} 50 {:.3f} {:.3f} {:.3f}\n".format(ztag, *q50arr))
      fout.write("{:s} 75 {:.3f} {:.3f} {:.3f}\n".format(ztag, *q75arr))

      # Now we add the marks of the quantiles 50 and 75: 
      for q50, q75, fmt, col in zip(q50arr, q75arr, formats, colors):
         yq = np.array([1.08, 1.15, 1.2])*max(ymax)
         xq50 = np.array([q50]*3)
         xq75 = np.array([q75]*3)
         
         ax.plot(xq50, yq, ls='-', c=col, lw=1.5)
         #ax.plot(xq75, yq, ls='--', c=col, lw=0.9)

      ax.set_xlim((-2.5, 2.5))
      if phi.max()>0:
         ax.set_ylim((0.0, max(ymax)*1.2))
      
      ax.text(0.04, 0.88, r"$z="+ztag+"$", transform=ax.transAxes)
      
      # some clearing to save memory:
      del sats
      del SAG_sim, SAG_fit, SAG_badfit, Groups, ParentGroup, Contam
   
   axs[2].set_xlabel(r'$\log(f_{\rm RPS})$')
   axs[2].legend(frameon=False, loc="upper right",
             numpoints=1, fontsize='small', handlelength=2.0, labelspacing=0.3)
   
   fout.close()
   fig.tight_layout()
   fig.subplots_adjust(hspace=0, wspace=0)
   fig.savefig(figs+"histfrps_comb.eps")


###############################################################################
###############################################################################
###############################################################################

def plot_mstarhist_qnt(centrals, cuts, sfx, snap=92, axptr=None):
 
   if not axptr:
      set_style('mnras', Hratio=0.75)
      fig, ax = plt.subplots()
   else:
      ax = axptr

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

   print(hosts)
   
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

   # Now we get the data and plot the different models: 
   for sag, tag, fmt, col, frpscut in zip(SAGdata, modelnames, formats, colors, cuts):
            
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
      logfrps = np.log10(strip_sum/mstar)

      # Select the galaxies according to the fRPS cut:
      smp = logfrps >= frpscut
      
      xbins = np.arange(7.4, 10, 0.15)
      #phi_i, bins = np.histogram(np.log10(mstar[smp]), bins=14, range=[7.5,10])
      phi_i, bins = np.histogram(np.log10(mstar[smp]), bins=xbins)
      x = (bins[1:]+bins[:-1])/2.0
      w = bins[1]-bins[0]
      phi = phi_i.astype(float)/w/float(len(centrals))
      
      if 'SPH'==tag:   # in this case we add a shaded region
         ax.bar(x, phi, width=w, ec=None, fc='#f0f0ff', fill=True)
      
      ax.step(x, phi, where='mid', ls=fmt, c=col, lw=1.5, label=tag)
      
      #phi = phi_i.astype(float)/float(len(centrals))
      #ax.bar(x, phi, width=w-dw, ec=None, fc=col, fill=True, alpha=0.15)
      #ax.bar(x, phi, width=w-dw, ec=col, fill=False, ls=fmt, lw=1, alpha=0.5)
      #ax.bar([1], [0], width=w-dw, ec=col, fc=col, lw=1, ls=fmt, fill=True, alpha=0.4, label=tag)
      #ax.set_ylabel(r'$N_{\rm gal} / N_{\rm gr}$')
      
      ax.set_ylabel(r'${\rm dN / d} \log (M_\star) / {\rm N_{gr}}$')
      ymax.append(phi.max())

   ax.set_xlim((7.5, 9.8))
   if phi.max()>0:
      ax.set_ylim((0.0, max(ymax)*1.2))
  
   ax.text(0.05, 0.88, r"$z="+s2z[snap]+"$", transform=ax.transAxes)
   
   if not axptr:
      ax.set_xlabel(r'$\log (M_\star[{\rm M}_\odot])$')
      ax.legend(frameon=False, #loc="upper left",
             numpoints=1, fontsize='small', handlelength=2, labelspacing=0.3)
   
      fig.tight_layout()
   
      fname = figs+'histmstar_s'+str(snap)+'frps_'+sfx+'.eps'
      fig.savefig(fname)

###############################################################################
###############################################################################
###############################################################################

def plot_ssfrhist_qnt(centrals, cuts, sfx, snap=92, axptr=None):
 
   if not axptr:
      set_style('mnras', Hratio=0.75)
      fig, ax = plt.subplots()
   else:
      ax = axptr
   
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

   print(hosts)
   
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

   # Now we get the data and plot the different models: 
   for sag, tag, fmt, col, frpscut in zip(SAGdata, modelnames, formats, colors, cuts):
            
      discmass = sag.readDataset('M_star_disk')
      bulgemass = sag.readDataset('M_star_bulge')
      sfr = sag.readDataset('SFR')
      sfr /= units.h
      mstar_all = (discmass + bulgemass)*units.mass.Msun/units.h
      del discmass, bulgemass
     
      # Select only the sample:
      mstar = mstar_all[sats]
      ssfr = sfr[sats]/mstar_all[sats]
      del mstar_all, sfr

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
      logfrps = np.log10(strip_sum/mstar)

      # Here we apply a manual correction to include all the quenched galaxies
      # (with ssfr=0):
      ssfr[ssfr<1e-13] = 1e-13
      logssfr = np.log10(ssfr)
      del ssfr

      # Select the galaxies according to the fRPS cut:
      smp = logfrps >= frpscut

      xbins = np.arange(-13.875, -7, 0.25)
      phi_i, bins = np.histogram(logssfr[smp], bins=xbins)
      x = (bins[1:]+bins[:-1])/2.0
      w = bins[1]-bins[0]
      phi = phi_i.astype(float)/w/float(len(centrals))
      
      if 'SPH'==tag:   # in this case we add a shaded region
         ax.bar(x, phi, width=w, ec=None, fc='#f0f0ff', fill=True)
      
      ax.step(x, phi, where='mid', ls=fmt, c=col, lw=1.5, label=tag)
     
      ax.set_ylabel(r'${\rm dN / d} \log ({\rm sSFR) / N_{gr}}$')
      ymax.append(phi.max())
   
   ax.set_xlim((-13.5, -8))
   #ax.set_xticklabels([r'$\leq -13\;\;\;$', r'$-12$', r'$-10$', r'$-8$'])
   if phi.max() > 0:
      ax.set_ylim((0.0, max(ymax)*1.2))
   
   ax.text(0.05, 0.88, r"$z="+s2z[snap]+"$", transform=ax.transAxes)
  
   if not axptr:
      ax.set_xlabel(r'$\log ({\rm sSFR [yr}^{-1}])$')
      ax.legend(frameon=False, #loc="upper left",
             numpoints=1, fontsize='small', handlelength=2, labelspacing=0.3)
   
      fig.tight_layout()

      fname = figs+'histssfr_s'+str(snap)+'frps_'+sfx+'.eps'
      fig.savefig(fname)


###############################################################################
###############################################################################
###############################################################################

def plot_colorhist_qnt(centrals, cuts, sfx, snap=92, axptr=None):
 
   if not axptr:
      set_style('mnras', Hratio=0.75)
      fig, ax = plt.subplots()
   else:
      ax = axptr
   
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

   print(hosts)
   
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

   # Now we get the data and plot the different models: 
   for sag, tag, fmt, col, frpscut in zip(SAGdata, modelnames, formats, colors, cuts):
            
      discmass = sag.readDataset('M_star_disk')
      bulgemass = sag.readDataset('M_star_bulge')
      mstar_all = (discmass + bulgemass)*units.mass.Msun/units.h
      del discmass, bulgemass
      mag_g = sag.readDataset('Magnitudes/Mag_gS_dust1')
      mag_r = sag.readDataset('Magnitudes/Mag_rS_dust1')
      g_r_all = mag_g - mag_r
      del mag_g, mag_r
    
      # Select only the sample:
      mstar = mstar_all[sats]
      g_r = g_r_all[sats]
      del mstar_all, g_r_all

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
      logfrps = np.log10(strip_sum/mstar)

      # Select the galaxies according to the fRPS cut:
      smp = logfrps >= frpscut

      xbins = np.arange(0, 1.05, 0.05)
      phi_i, bins = np.histogram(g_r[smp], bins=xbins)
      x = (bins[1:]+bins[:-1])/2.0
      w = bins[1]-bins[0]
      phi = phi_i.astype(float)/w/float(len(centrals))

      if 'SPH'==tag:   # in this case we add a shaded region
         ax.bar(x, phi, width=w, ec=None, fc='#f0f0ff', fill=True)
      
      ax.step(x, phi, where='mid', ls=fmt, c=col, lw=1.5, label=tag)

      ax.set_ylabel(r'${\rm dN / d} M / {\rm N_{gr}}$')
      ymax.append(phi.max())
   
   ax.set_xlim((0.0,0.9))
   if phi.max()>0:
      ax.set_ylim((0.0, max(ymax)*1.2))
      
   ax.text(0.05, 0.88, r"$z="+s2z[snap]+"$", transform=ax.transAxes)
  
   if not axptr: 
      ax.set_xlabel(r'$(g - r)$')
      fig.tight_layout()

      ax.legend(frameon=False, #loc="upper left",
                numpoints=1, fontsize='small', handlelength=2, labelspacing=0.3)
      
      fig.tight_layout()
         
      fname = figs+'histcolor_s'+str(snap)+'frps_'+sfx+'.eps'
      
      fig.savefig(fname)

###############################################################################
###############################################################################
###############################################################################


if __name__ == '__main__':
    
   centrals = select_galaxies_1model()
   
   ### This creates a file with the quantile values
   plot_frps_comb(centrals)
   #exit()
   
   fin = open("frps_quantiles.txt", "r")
   fin.readline()
   
   ### Individual plots: 
   #for line in fin.readlines():
   #   ztag = line.split()[0]
   #   qnt = 'q'+line.split()[1]
   #   cuts = [float(i) for i in line.split()[2:]]

   #   #plot_mstarhist_qnt(centrals, cuts, qnt, snap=z2s[ztag])
   #   #plot_ssfrhist_qnt(centrals, cuts, qnt, snap=z2s[ztag])
   #   #plot_colorhist_qnt(centrals, cuts, qnt, snap=z2s[ztag])

   #   qnt = 'q'+line.split()[1]+'SPH'
   #   cuts = [float(line.split()[2])]*3
   #   
   #   plot_mstarhist_qnt(centrals, cuts, qnt, snap=z2s[ztag])
   #   plot_ssfrhist_qnt(centrals, cuts, qnt, snap=z2s[ztag])
   #   plot_colorhist_qnt(centrals, cuts, qnt, snap=z2s[ztag])

   ## Combining all the plots (50 qnt)
   set_style('mnras-fw', Hratio=0.618)
   fig, axs = plt.subplots(nrows=3, ncols=3, sharex='col')
   pltr = {'0.0':0, '1.5':1, '3.0':2}

   fin.seek(0)
   fin.readline()
   for line in fin.readlines():
      ztag = line.split()[0]
      if line.split()[1] == '50':
         qnt = 'q'+line.split()[1]+'SPH'
         cuts = [float(line.split()[2])]*3

         plot_mstarhist_qnt(centrals, cuts, qnt, snap=z2s[ztag], axptr=axs[pltr[ztag],0])
         plot_ssfrhist_qnt(centrals, cuts, qnt, snap=z2s[ztag], axptr=axs[pltr[ztag],1])
         plot_colorhist_qnt(centrals, cuts, qnt, snap=z2s[ztag], axptr=axs[pltr[ztag],2])
     
   axs[pltr['3.0'],0].set_xlabel(r'$\log (M_\star[{\rm M}_\odot])$')
   axs[pltr['3.0'],1].set_xlabel(r'$\log ({\rm sSFR [yr}^{-1}])$')
   axs[pltr['3.0'],2].set_xlabel(r'$(g - r)$')

   axs[pltr['0.0'],0].legend(frameon=False, #loc="upper right",
             numpoints=1, fontsize='small', handlelength=2, labelspacing=0.3)
   axs[pltr['1.5'],1].legend(frameon=False, #loc="upper left",
             numpoints=1, fontsize='small', handlelength=2, labelspacing=0.3)
   axs[pltr['3.0'],2].legend(frameon=False, #loc="upper left",
             numpoints=1, fontsize='small', handlelength=2, labelspacing=0.3)
   
   axs[pltr['3.0'],1].set_xticks([-13, -12, -11, -10, -9, -8])
   axs[pltr['3.0'],1].text(-0.02, -0.085, r"$\leq$", 
                           transform=axs[pltr['3.0'],1].transAxes)
   
   fig.tight_layout()
   fig.subplots_adjust(hspace=0.05)

   fname = figs+'hist_galpropZ_frps_q50SPH.eps'
   fig.savefig(fname)

   ### Combining all the plots (75 qnt)
   set_style('mnras-fw', Hratio=0.618)
   fig, axs = plt.subplots(nrows=3, ncols=3, sharex='col')
   pltr = {'0.0':0, '1.5':1, '3.0':2}

   fin.seek(0)
   fin.readline()
   for line in fin.readlines():
      ztag = line.split()[0]
      if line.split()[1] == '75':
         qnt = 'q'+line.split()[1]+'SPH'
         cuts = [float(line.split()[2])]*3

         plot_mstarhist_qnt(centrals, cuts, qnt, snap=z2s[ztag], axptr=axs[pltr[ztag],0])
         plot_ssfrhist_qnt(centrals, cuts, qnt, snap=z2s[ztag], axptr=axs[pltr[ztag],1])
         plot_colorhist_qnt(centrals, cuts, qnt, snap=z2s[ztag], axptr=axs[pltr[ztag],2])
     
   axs[pltr['3.0'],0].set_xlabel(r'$\log (M_\star[{\rm M}_\odot])$')
   axs[pltr['3.0'],1].set_xlabel(r'$\log ({\rm sSFR [yr}^{-1}])$')
   axs[pltr['3.0'],2].set_xlabel(r'$(g - r)$')

   axs[pltr['0.0'],0].legend(frameon=False, #loc="upper right",
             numpoints=1, fontsize='small', handlelength=2, labelspacing=0.3)
   axs[pltr['1.5'],1].legend(frameon=False, #loc="upper left",
             numpoints=1, fontsize='small', handlelength=2, labelspacing=0.3)
   axs[pltr['3.0'],2].legend(frameon=False, #loc="upper left",
             numpoints=1, fontsize='small', handlelength=2, labelspacing=0.3)
   
   axs[pltr['3.0'],1].set_xticks([-13, -12, -11, -10, -9, -8])
   axs[pltr['3.0'],1].text(-0.02, -0.085, r"$\leq$", 
                           transform=axs[pltr['3.0'],1].transAxes)
   
   fig.tight_layout()
   fig.subplots_adjust(hspace=0.05)

   fname = figs+'hist_galpropZ_frps_q75SPH.eps'
   fig.savefig(fname)

