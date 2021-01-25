#! /usr/bin/env python 

import numpy as np
import struct 

def read_properties(fname):
   
   f = open(fname, "rb")
   p = {}
   
   Ngroups = np.fromfile(f, dtype=np.int32, count=1)[0]
   p['Ngroups'] = Ngroups

   p['GroupM200'] = np.zeros(Ngroups+1, dtype=np.float32)
   p['GroupM200'][1:] = np.fromfile(f, dtype=np.float32, count=Ngroups)

   p['GroupR200'] = np.zeros(Ngroups+1, dtype=np.float32)
   p['GroupR200'][1:] = np.fromfile(f, dtype=np.float32, count=Ngroups)

   p['GroupMvir'] = np.zeros(Ngroups+1, dtype=np.float32)
   p['GroupMvir'][1:] = np.fromfile(f, dtype=np.float32, count=Ngroups)

   p['GroupRvir'] = np.zeros(Ngroups+1, dtype=np.float32)
   p['GroupRvir'][1:] = np.fromfile(f, dtype=np.float32, count=Ngroups)

   p['GroupCenter'] = np.zeros((Ngroups+1,3), dtype=np.float32)
   for i in range(Ngroups):
      p['GroupCenter'][i+1] = np.fromfile(f, dtype=np.float32, count=3)

   p['GroupMostBound'] = np.zeros(Ngroups+1, dtype=np.int32)
   p['GroupMostBound'][1:] = np.fromfile(f, dtype=np.int32, count=Ngroups)

   p['GroupPosCM'] = np.zeros((Ngroups+1,3), dtype=np.float32)
   for i in range(Ngroups):
      p['GroupPosCM'][i+1] = np.fromfile(f, dtype=np.float32, count=3)

   p['GroupVelCM'] = np.zeros((Ngroups+1,3), dtype=np.float32)
   for i in range(Ngroups):
      p['GroupVelCM'][i+1] = np.fromfile(f, dtype=np.float32, count=3)
  
   p['GroupEkin'] = np.zeros(Ngroups+1, dtype=np.float32)
   p['GroupEkin'][1:] = np.fromfile(f, dtype=np.float32, count=Ngroups)

   p['GroupEpot'] = np.zeros(Ngroups+1, dtype=np.float32)
   p['GroupEpot'][1:] = np.fromfile(f, dtype=np.float32, count=Ngroups)

   p['GroupLambda'] = np.zeros(Ngroups+1, dtype=np.float32)
   p['GroupLambda'][1:] = np.fromfile(f, dtype=np.float32, count=Ngroups)

   p['GroupL'] = np.zeros((Ngroups+1,3), dtype=np.float32)
   for i in range(Ngroups):
      p['GroupL'][i+1] = np.fromfile(f, dtype=np.float32, count=3)

   p['GroupVmax'] = np.zeros(Ngroups+1, dtype=np.float32)
   p['GroupVmax'][1:] = np.fromfile(f, dtype=np.float32, count=Ngroups)

   f.close()

   return p


def read_contamination(fname):
   
   f = open(fname, "rb")
   p = {}
   
   Ngroups = np.fromfile(f, dtype=np.int32, count=1)[0]
   p['Ngroups'] = Ngroups

   p['GroupContCount'] = np.zeros(Ngroups+1, dtype=np.int32)
   p['GroupContCount'][1:] = np.fromfile(f, dtype=np.int32, count=Ngroups)

   p['GroupMinDist'] = np.zeros(Ngroups+1, dtype=np.float32)
   p['GroupMinDist'][1:] = np.fromfile(f, dtype=np.float32, count=Ngroups)

   f.close()

   return p



