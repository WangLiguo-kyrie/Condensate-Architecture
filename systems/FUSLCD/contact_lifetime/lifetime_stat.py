#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
for the calculation of life time of interchain contact
"""

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import os,glob
import itertools
import pickle

    
life_files=sorted(glob.glob('contacts_life/*.p'))
life_collect=[]
life_collect_frag=[]
for index, life_file in enumerate(life_files):
    contact_dic = pickle.load( open( life_file, "rb" ))
    print(life_file)
    print('Process: index {}/{}'.format(index,len(life_files)))
    for keys, values in contact_dic.items():
        #print(keys)
        #print(values)
        life_collect_frag=life_collect_frag+values
    if ((index+2)%10)==0:
        print(life_collect_frag)
        life_collect=life_collect+life_collect_frag
        life_collect_frag=[]
        print(len(life_collect))
life_collect=np.array(life_collect)
pickle.dump(life_collect, open( "./contacts_life_collect.p", "wb" ))
print(life_collect.shape)
#life_collect = pickle.load( open( "./contacts_life_collect.p", "rb" ))

