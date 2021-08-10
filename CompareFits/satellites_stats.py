#!/usr/bin/env python3

import matplotlib as mpl
mpl.use('PS')
import matplotlib.pyplot as plt

import h5py 
import numpy as np

f = h5py.File("rampressure_092.h5", "r")

m200 = f['M200_host'][:]
rrel = f['rrel'][:]/f['R200_host'][:]

f.close()

edges = np.arange(11, 16.5)
hist, _ = np.histogram(np.log10(np.unique(m200)), bins=edges)

print(edges)
print(hist)

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

set_style('mnras', Hratio=0.75)

fig, ax = plt.subplots()
figa, axa = plt.subplots()

lshape = [':',(0,(4,2,1,2,1,2)),'-.','--','-']
lscols = ['c','orange','r','m','b']

massbins = [11,12,13,14,15,16]
deltax = 0.2
rbins = np.arange(0.0, 4.0, deltax)
x = (rbins[1:]+rbins[:-1])/2.0

for i in range(len(massbins)-1):
   mmin = massbins[i]
   mmax = massbins[i+1]
   mflt = np.where((mmin<np.log10(m200))&(np.log10(m200)<mmax))[0]
   
   hist, _ = np.histogram(rrel[mflt], bins=rbins)
   pdf = hist.astype(float)/hist.sum()

   #print(pdf.sum())

   cumsum = np.cumsum(pdf)
      
   if i!=(len(massbins)-2):
      label = r"$"+str(mmin)+"\leq \log M_{200} < "+str(mmax)+"$"
   else:
      label = r"$"+str(mmin)+"\leq \log M_{200}$"
   
   ax.plot(x, pdf, c=lscols[i], ls=lshape[i], lw=1, label=label)
   axa.plot(x, cumsum, c=lscols[i], ls=lshape[i], lw=1, label=label)

ax.set_xlim((0,3.0))
ax.set_ylim((0,0.2))
ax.set_yticks(np.arange(0,0.21, 0.05))
ax.set_xlabel(r"$r/R_{200}$")
ax.set_ylabel(r"PDF($N_\mathrm{sat} | r/R_{200}$)")
ax.legend(fontsize='small', frameon=False, handlelength=2.5,
               labelspacing=0.1, loc='upper right')

fig.tight_layout()
fig.savefig("figs/satellites.eps")


axa.set_xlim((0,3.0))
axa.set_ylim((0,1))
axa.set_xlabel(r"$r/R_{200}$")
axa.set_ylabel(r"$f_\mathrm{sat}(<r/R_{200})$")
axa.legend(fontsize='small', frameon=False, handlelength=2.5,
               labelspacing=0.1, loc='lower right')

figa.tight_layout()
figa.savefig("figs/satellites_acc.eps")

