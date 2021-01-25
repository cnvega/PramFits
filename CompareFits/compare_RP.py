#!/usr/bin/env python 

import SAGreader as SAG
import numpy as np
import os 
import h5py

import groups

def Pram_beta(r, p0, rs, beta):
    return p0*(1.+(r/rs)**2.)**(-3./2.*beta)

# The units:
L_cm = 3.085678e21
V_cms = 1.0e5
M_gr  = 1.989e43
factor = L_cm*L_cm/V_cms/V_cms*L_cm/M_gr

# El fit de Tomas:
def rpfit_tomas(m200, r200, aexp, r):
    ap = -0.8 + 1.2*(np.log10(m200) - 12.)
    bp = 1.2 - 0.4*(np.log10(m200) - 12.)
    ar = 0.59 - 0.14*(np.log10(m200) - 12.)
    br = -0.44 + 0.12*(np.log10(m200) - 12.)
    ab, bb = 0.92, -0.4
    
    beta = ab + bb*(aexp - 0.25)
    rs = r200*(ar + br*(aexp - 0.25))
    p0  = (1e-12)*10**(ap+bp*(aexp - 0.25))
    
    return Pram_beta(r, p0, rs, beta)


def rpfit_cvm(m200, r200, aexp, r):
    rrel  = r/r200

    p0 = (1e-12)*10**(7.01*aexp**(-0.122) - 10 + 0.83)
    rs = -3.4*aexp**(-0.42) + 10.2
    bm  = 3.36e-3*aexp**(1.33) + 0.512
    bn  = -5.6 # fixed
    beta = bm*np.log10(m200) + bn

    return p0*(rrel/rs)**(-3./2.*beta)

snaps = [92, 36, 26]
scales = [1.0, 0.402041, 0.248]
redshifts = ['0.0','1.5','3.0']


fname   = "icmpropertiesJ25_NFW_itf_"

datapath = "/data1/DolagClusters/"
cllist = os.listdir(datapath)

for snap, scale, redshift in zip(snaps, scales, redshifts):

   print("Snapshot "+str(snap))
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
   
   
   # Cargamos también las propiedades de los grupos incluyendo la contaminación para quedarnos con todos los clusters y sólo descartar los contaminados.
   
   
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
   
   
   grp['Ngroups'], grp['Ngroups_cluster'], grp['cluster_name']
   
   # Para leer y asignar bien las GroupM200/R200 tenemos que hacer un ParentGroup 
   # contextual, ad-hoc al nuevo sistema de ordenamiento que incluye todos los clusters
   # y los trata como uno solo.
   
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
   
   medianRhoICM = fdata.readDataset('medianRhoICM', idxfilter=flt)
   relvel2ICM = fdata.readDataset('relvel2ICM', idxfilter=flt)

   unit = fdata.dataList[0]['medianRhoICM'].attrs['Unit']
   medianRhoICM *= unit   ## h^2 gr cm^-3
   
   unit = fdata.dataList[0]['relvel2ICM'].attrs['Unit']
   relvel2ICM *= unit    ## (cm/s)^2

   pos = fdata.readDataset('Pos', idxfilter=flt) # comoving
   pos *= scale  # physical
   
   rampress = medianRhoICM*relvel2ICM   # in physical units
   
   
   # Y asignamos las M200/R200
   
   parentM200 = np.zeros(rampress.shape)
   parentR200 = np.zeros(rampress.shape)
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
   
   
   rrel = np.zeros(parentgroup.shape)
   for i, p in enumerate(pos):
       rrel[i] = ((p[0]-parentPos[i][0])**2 + (p[1]-parentPos[i][1])**2 + 
                  (p[2]-parentPos[i][2])**2 )**0.5
   
   # and here we have rrel and parentPos in physical units.
   
   
   # Vemos como da el fit original:
   rpf     = rpfit_tomas(parentM200*1e10, parentR200, scale, rrel)
   rpf_new = rpfit_cvm(parentM200*1e10, parentR200, scale, rrel)
   
   f = h5py.File("rampressure_"+str(snap).zfill(3)+".h5", "w")
   
   f['rp_fit'] = rpf
   f['rp_fit_new'] = rpf_new
   f['rp'] = rampress
   f['M200_host'] = 1e10*parentM200
   f['R200_host'] = parentR200
   f['rrel']  = rrel
   f.attrs['redshift'] = redshift
   f.attrs['scale'] = scale
   
   f.close()

   del fdata


