
# coding: utf-8

# Vemos si podemos cargar los icmpropeties con el SAGreader:

import SAGreader as SAG
import numpy as np
import os 
import h5py

import groups

import matplotlib.pylab as pl

from scipy.optimize import curve_fit
import sys

# Leemos el archivo con los parámetros de escala por snapshot:

snapList, scales, redshifts = np.loadtxt("scales.txt", unpack=True)

snapList = snapList.astype(np.int)

for snap in snapList[10:]:

   print("### Snapshot: ", snap)
   sys.stdout.flush()
   scale   = scales[snapList==snap]
   
   fname   = "icmpropertiesJ25_NFW_itf_"
   
   figs = 'figs/snap'+str(snap).zfill(3)+'_'
   
   datapath = "/data/DolagClusters/"
   cllist = os.listdir(datapath)
   datafiles = []
   propfiles = []
   contfiles = []
   clnames = []
   for i in cllist:
       if 'ovisc' in i:
           datafiles.append(datapath+i+"/Postprocessing/ICMProperties/"+fname+str(snap).zfill(3)+".hdf5")
           propfiles.append(datapath+i+"/Postprocessing/Properties/properties_"+str(snap).zfill(3))
           contfiles.append(datapath+i+"/Postprocessing/Contamination/contamination_"+str(snap).zfill(3))
           clnames.append(i.replace("_ovisc",""))
   
   
   # Cargamos los datos del ICM con la runtina del SAG
   
   fdata = SAG.SAGdata(None, None)
   
   for fi in datafiles:
          fdata.addFile(fi)
   
   tipos_all = fdata.readDataset('Type')
   
   
   # Cargamos también las propiedades de los grupos incluyendo la contaminación 
   # para quedarnos con todos los clusters y sólo descartar los contaminados.
   
   
   for i, fnm in enumerate(propfiles):
       if 0 == i:
           grp = groups.read_properties(fnm)
           cont = groups.read_contamination(contfiles[i])
           grp['GroupContCount'] = cont['GroupContCount']
           grp['Ngroups_cluster'] = [grp['Ngroups']]
           grp['cluster_name'] = [clnames[i]]
       else:
           tmp = groups.read_properties(fnm)
           for k in grp.keys():
               if k not in ['Ngroups','cluster_name','Ngroups_cluster', 'GroupContCount']:
                   grp[k] = np.concatenate([grp[k], tmp[k][1:]])
           tmp = groups.read_contamination(contfiles[i])
           grp['GroupContCount'] = np.concatenate([grp['GroupContCount'],tmp['GroupContCount'][1:]])
           grp['Ngroups'] += tmp['Ngroups']
           grp['Ngroups_cluster'].append(tmp['Ngroups'])
           grp['cluster_name'].append(clnames[i])
   
   
   # Para leer y asignar bien las GroupM200/R200 tenemos que hacer un ParentGroup contextual, 
   # ad-hoc al nuevo sistema de ordenamiento que incluye todos los clusters y los trata como uno solo.
   
   ngal_cluster = []
   for fd in fdata.dataList:
       ngal_cluster.append(len(fd['ParentGroup']))
   ngal_cluster
   
   parentCorrect = np.zeros((len(tipos_all),1), dtype=np.int32)
   
   totsumgal = 0
   for i, ng in enumerate(grp['Ngroups_cluster'][:-1]):
       totsumgal += ngal_cluster[i]
       parentCorrect[totsumgal:] += ng
   
   parentgroup_all = fdata.readDataset('ParentGroup')
   parentgroup_all += parentCorrect
   
   
   # Hacemos un filtro que incluya sólo las galaxias tipo 1 de halos no contaminados.
   
   contgal = np.zeros((len(tipos_all),1), dtype=np.bool)
   
   for i in range(len(tipos_all)):
       if 0 < grp['GroupContCount'][parentgroup_all[i]]:
           contgal[i] = True
           
   flt = np.where((tipos_all==1)&(contgal==False))[0]
   
   parentgroup = parentgroup_all[flt]
   
   # Cargamos los datasets usando ese filtro:
   
   medianRhoICM = fdata.readDataset('medianRhoICM', idxfilter=flt)
   relvel2ICM = fdata.readDataset('relvel2ICM', idxfilter=flt)
   pos = fdata.readDataset('Pos', idxfilter=flt) # comoving
   pos *= scale  # physical
   
   rampress = np.log10(medianRhoICM*relvel2ICM)   # in physical units
   
   # Y asignamos las M200/R200
   
   parentM200 = np.zeros(len(rampress))
   parentR200 = np.zeros(len(rampress))
   for i, par in enumerate(parentgroup):
       parentM200[i] = grp['GroupM200'][par]
       parentR200[i] = grp['GroupR200'][par]
   
   
   # Necesitamos además las posiciones de sus respectivas centrales y la distancia a las mismas.
   
   # La posición de las centrales de cada tipo 1
   parentPos = np.zeros(pos.shape)
   
   # Creamos un array con las posiciones de todas las centrales en orden
   pos_central = fdata.readDataset('Pos', idxfilter=np.where(tipos_all==0)[0])
   pos_central *= scale
   
   # y otro con la ID del grupo de cada una de ellas
   parent_central = parentgroup_all[tipos_all == 0]
   
   for i, par in enumerate(parentgroup):
       parentPos[i] = pos_central[np.where(parent_central==par)[0]]    
   
   rrel = np.zeros(len(parentgroup))
   for i, p in enumerate(pos):
       rrel[i] = ((p[0]-parentPos[i][0])**2 + (p[1]-parentPos[i][1])**2 + (p[2]-parentPos[i][2])**2 )**0.5
   
   # and here we have rrel and parentPos in physical units.
   
   
   # Vemos la distribución de masas de la muestra de galaxias tipo 1
   
   
   counts, edges, patch = pl.hist(np.log10(parentM200*1e10), bins=10, range=[9.5,15.5], log=True)
   pl.ylabel(r'$N_{\rm gal}$')
   pl.xlabel(r'$M_{200} [h^{-1} M_{\odot}]$')
   pl.savefig(figs+'masses.eps')
   pl.close()
   
   # Vemos el comportamiento de los puntos a diferentes bins de masa:
   
   pl.figure()
   maxbin = int(max(np.log10(parentM200*1e10)))+1
   mb = []
   for i in range(4):
       mb.append((maxbin-i-1, maxbin-i))
   plts  = [(0,0), (0,1), (1,0), (1,1)]
   
   fig, axs = pl.subplots(2,2, sharex='col', sharey='row')
   for i in range(4):
       massflt_bin = np.where((mb[i][0]<=np.log10(parentM200*1e10))
                          &(np.log10(parentM200*1e10)<mb[i][1])&
                          ((rrel/parentR200) <= 2.5) )[0]
       ax = axs[plts[i]]
       ax.plot(rrel/parentR200, rampress, 'y.')
       ax.plot(rrel[massflt_bin]/parentR200[massflt_bin], rampress[massflt_bin], 'b.')
       ax.set_xlim((0,2.1))
       if (plts[i][1] == 0):
           ax.set_ylabel(r"$\log_{10}(P_{\rm ram})$")
       if (plts[i][0] == 1):
           ax.set_xlabel(r"$r/R_{200}$")
       s = r"$"+str(mb[i][0])+" < \log(M_{200}) < "+str(mb[i][1])+"$"
       ax.text(0.43, 0.85, s, transform=ax.transAxes)
   #pl.tight_layout()
   fig.savefig(figs+"RP.eps")
   pl.close()
   
   # Fiteamos todo junto:
   
   
   ### el fit original me causa un overflow. Le escalo rrel a rrel/R200 para corregirlo:
   dflt = np.where(rrel/parentR200 <= 1.)[0]
   xx = np.array([rrel[dflt]/parentR200[dflt], parentM200[dflt]])
   yy = rampress[dflt].T[0]
   
   def Pram(rrel_mass, p0, rs, b_m):
       r    = rrel_mass[0]
       b_n = -5.6
       beta = b_m*np.log10(rrel_mass[1]*1e10)+b_n
       
       #return 10**p0*(1+(r/rs)**2)**(-3./2.*beta)
       return p0-3./2.*beta*np.log10(r/rs)
   
   #ppot, pcov = curve_fit(Pram, xx, yy, bounds=([-10,0,-100,-100],[10,10,100,100]))
   ppot, pcov = curve_fit(Pram, xx, yy, bounds=([-10,0,0],[10,10,10]))

   print("Params: ", ppot)
   print("Cov. matrix: ", pcov)
   sys.stdout.flush()
   
   pl.figure()
   
   maxbin = int(max(np.log10(parentM200*1e10)))+1
   mb = []
   
   for i in range(4):
       mb.append((maxbin-i-1, maxbin-i))
   
   plts  = [(0,0), (0,1), (1,0), (1,1)]
   fig, axs = pl.subplots(2,2, sharex='col')
   for i in range(4):
       massflt_bin = np.where((mb[i][0]<=np.log10(parentM200*1e10))
                          &(np.log10(parentM200*1e10)<mb[i][1])&
                          ((rrel/parentR200) <= 2.5) )[0]
       x_bin      = np.array([rrel[massflt_bin]/parentR200[massflt_bin], parentM200[massflt_bin]])
       pram_curve = Pram(x_bin, *ppot )
       
       ax = axs[plts[i]]
       ax.plot(rrel/parentR200, rampress, 'y.')
       ax.plot(rrel[massflt_bin]/parentR200[massflt_bin], rampress[massflt_bin], 'b.')
       ax.plot(x_bin[0], pram_curve, 'r.')
       ax.set_xlim((0,2.5))
       if (plts[i][1] == 0):
           ax.set_ylabel(r"$\log_{10}(P_{\rm ram})$")
       if (plts[i][0] == 1):
           ax.set_xlabel(r"$r/R_{200}$")
       s = r"$"+str(mb[i][0])+" < \log(M_{200}) < "+str(mb[i][1])+"$"
       ax.text(0.43, 0.85, s, transform=ax.transAxes)
   fig.savefig(figs+"RPfit.eps")
   pl.close()
   
   # Guardamos este fit en el hdf5:
   
   fout = h5py.File("RPfits.h5")
   
   grname = str(snap)
   
   if grname in fout.keys(): del fout[grname]
   
   fout.create_group(grname)
   
   fout[grname+'/params'] = np.array(ppot)
   fout[grname+'/cov_mtx'] = np.array(pcov)
   fout[grname].attrs['Params'] = "[p_0, r_s, b_m]"
   fout[grname].attrs['numpoints'] = np.array(len(xx[0]))
   fout[grname].attrs['scale'] = scale
   fout[grname].attrs['redshift'] = redshifts[snapList==snap]
   
   fout.close()
   

