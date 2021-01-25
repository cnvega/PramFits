#!/usr/bin/env python 

import matplotlib as mpl

#mpl.use('Agg')
mpl.use('PS')

import matplotlib.pylab as pl

import numpy as np
import h5py

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

snaps = [92, 36, 26]
scales = [1.0, 0.402041, 0.248]
redshifts = ['0.0','1.5','3.0']

massbins = [11,12,13,14,15,16]

lshape = [':',(0,(4,2,1,2,1,2)),'-.','--','-']
lscols = ['c','orange','r','m','b']

rdelta = 0.1
rbins = np.arange(0.0, 1.2, rdelta)
print(rbins)

set_style(style='mnras', Hratio=1.3, Wfrac=1)

fig, axs = pl.subplots(nrows=3, sharex=True)

l10 = np.log(10)

for snap, ztag, ax in zip(snaps, redshifts, axs):

   f = h5py.File("rampressure_"+str(snap).zfill(3)+".h5", "r")
   
   rpf = f['rp_fit'][:] 
   rampress = f['rp'][:]
   M200 = f['M200_host'][:]
   R200 = f['R200_host'][:] 
   rrel = f['rrel'][:] 
   f.close()

   rpdiff = rpf/rampress
   rdist  = rrel/R200

   for i in range(len(massbins)-1):
      mmin = massbins[i]
      mmax = massbins[i+1]
      mflt = np.where((mmin<np.log10(M200))&(np.log10(M200)<mmax))[0]

      rpdiff_mb = rpdiff[mflt]
      rdist_mb  = rdist[mflt]

      rp = np.zeros(len(rbins)-1)
      rpd = np.zeros(len(rbins)-1)
      r = np.zeros(len(rbins)-1)
      
      for b in range(len(rbins)-1):
         rmin, rmax = rbins[b], rbins[b+1]
         rflt = np.where((rdist_mb >= rmin)&(rdist_mb < rmax))[0]
         rp[b] = rpdiff_mb[rflt].mean()
         rpd[b] = rpdiff_mb[rflt].std()
         r[b]  = (rmin+rmax)/2.
      
      if (snap==92 and i==(len(massbins)-2)) or \
         (snap==36 and i==(len(massbins)-3)) or \
         (snap==26 and i==(len(massbins)-4)):
         x1, x2 = r[1]-rdelta/2, r[1]+rdelta/2
         y1,y2=np.log10(rp[1])-rpd[1]/rp[1]/l10,np.log10(rp[1])+rpd[1]/rp[1]/l10
         ax.plot([x1,x2], [y1,y1],c=lscols[i],lw=0.5)
         ax.plot([x1,x2], [y2,y2],c=lscols[i],lw=0.5)
         ax.plot([x1,x1], [y1,y2],c=lscols[i],lw=0.5)
         ax.plot([x2,x2], [y1,y2],c=lscols[i],lw=0.5)
         ax.plot(r[1], np.log10(rp[1]),'o',c=lscols[i], ms=1.5)
      if i==0:
         x1, x2 = r[-3]-rdelta/2, r[-3]+rdelta/2
         y1,y2=np.log10(rp[-3])-rpd[-3]/rp[-3]/l10,np.log10(rp[-3])+rpd[-3]/rp[-3]/l10
         ax.plot([x1,x2], [y1,y1],'c',lw=0.5)
         ax.plot([x1,x2], [y2,y2],'c',lw=0.5)
         ax.plot([x1,x1], [y1,y2],'c',lw=0.5)
         ax.plot([x2,x2], [y1,y2],'c',lw=0.5)
         ax.plot(r[-3], np.log10(rp[-3]), 'oc', ms=1.5)
     
      if i!=(len(massbins)-2):
         label = r"$"+str(mmin)+"\leq \log M_{200} < "+str(mmax)+"$"
      else:
         label = r"$"+str(mmin)+"\leq \log M_{200}$"
      ax.plot(r, np.log10(rp), c=lscols[i], linestyle=lshape[i], lw=1, label=label)
      ax.plot([0,1],[0,0],'k',lw=0.3)
      ax.text(0.85, 0.87, r"$z = "+ztag+"$", transform=ax.transAxes)

   ax.set_xlim((0.,1.))
   ax.set_ylabel(r'$\log ({\rm P}_{\rm ram}^{\rm \;\; fit} / {\rm P}_{\rm ram}^{\rm \;\; sim})$')
   if snap == 26:
      ax.set_xlabel(r'$ r/R_{200}$')
      leg = ax.legend(fontsize='small', frameon=False, handlelength=2.5,
                      labelspacing=0.1, loc='upper left')

axs[0].set_ylim((-2.5,2.5))
axs[1].set_ylim((-2.5,2.5))
axs[2].set_ylim((-2.5,2.5))
fig.tight_layout()
fig.savefig("figs/rpcompcgs.eps")

