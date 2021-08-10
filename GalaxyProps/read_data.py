
import SAGreader as SAG
import numpy as np
import os 
import h5py

import groups

import matplotlib.pylab as pl

import sys

snapList, scales, redshifts = np.loadtxt("scales.txt", unpack=True)
snapList = snapList.astype(np.int)

Path_Clusters = "data/subfind" 
Path_SAG = "data/galaxies"

def _read_catalogues(snap):
   
   print("### Reading catalogues, snapshot:", snap)
   sys.stdout.flush()
   scale   = scales[snapList==snap]
   
   datapath = Path_Clusters+'/'
   cllist = sorted(os.listdir(datapath))
   datafiles = []
   propfiles = []
   contfiles = []
   clnames = []
   for i in cllist:
       if 'ovisc' in i:
           propfiles.append(datapath+i+"/Postprocessing/Properties/properties_"+str(snap).zfill(3))
           contfiles.append(datapath+i+"/Postprocessing/Contamination/contamination_"+str(snap).zfill(3))
           clnames.append(i.replace("_ovisc",""))
   
   
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
   
   return grp   


def _read_sag_model(modelname, snap):

   print("### Reading SAG outputs, snapshot: ", snap)
   sys.stdout.flush()
   scale   = scales[snapList==snap]
    
   fname   = "gal_itf_"+str(snap).zfill(3)+"_"+modelname+".hdf5"
   
   
   datapath = Path_SAG+'/'
   cllist = sorted(os.listdir(datapath))
   datafiles = []
   for i in cllist:
       if 'ovisc' in i:
           datafiles.append(datapath+i+"/"+modelname+"/"+fname)
   
   fdata = SAG.SAGdata(None, None)
   
   for fi in datafiles:
          fdata.addFile(fi)
   
   return fdata

def read_data(modelname, snap):
   
   grp =  _read_catalogues(snap)

   sag =  _read_sag_model(modelname, snap)

   ngal_cluster = []
   for fd in sag.dataList:
       ngal_cluster.append(len(fd['MainHaloID']))
   ngal_cluster
   
   parentCorrect = np.zeros((sum(ngal_cluster),1), dtype=np.int32)
   
   totsumgal = 0
   for i, ng in enumerate(grp['Ngroups_cluster'][:-1]):
       totsumgal += ngal_cluster[i]
       parentCorrect[totsumgal:] += ng
   
   parentgroup = sag.readDataset('MainHaloID')
   parentgroup += parentCorrect
  
   # The contamination flag:
   contgal = np.zeros(parentgroup.shape, dtype=np.bool)
   for i in range(len(parentgroup)):
       if 0 < grp['GroupContCount'][parentgroup[i]]:
           contgal[i] = True

   return sag, grp, parentgroup, contgal

