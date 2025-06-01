#! /usr/bin/python

import sys
import os
import numpy as np
import MDAnalysis as mda
from MDAnalysis.analysis import contacts
from MDAnalysis import *
import math
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager

trajectory_dir='/grain/liguo/MYC/IDP-conf-change/condensate/7-M3IDP-dih7-SC15/HSPB2/more_copy/'

u = mda.Universe(trajectory_dir+'production-protein.tpr',trajectory_dir+'production2-Verlet-pbc-protein.xtc')
condensate=u.select_atoms('all')
#condensate=u.select_atoms('resid 61-146')   #only use folded domain to define COG to avoid COG fluctuation due to IDR motion
ag = u.select_atoms('name BB')
print(len(ag.atoms))
fragments = ag.atoms.fragments

#discard initial 2us 
start = 4000
end = -1
#every 200ns
step = 400
monomer_length=182

dist_result=[]
for index, ts in enumerate(u.trajectory[start:end:step]):  
    dist_collect=[]
    dist_collect.append(ts.time)
    for i, frag in enumerate(fragments):
        foldD_start=i*monomer_length+61
        foldD_end=i*monomer_length+146
        foldD=frag.select_atoms('resid {}-{}'.format(foldD_start,foldD_end))
        print(len(foldD.atoms))
        
        #distance based on folded domain COG to condensate COG to determine if belongs to interface region
        dist = contacts.distance_array(condensate.center_of_geometry(), foldD.center_of_geometry())
        dist_collect.append(dist[0][0])
    dist_result.append(dist_collect)
dist_result=np.array(dist_result)
print(dist_result.shape)
np.savetxt('droplet_foldD_distr.xvg',dist_result,delimiter='	')

interface_list=[]
core_list=[]
for k in range(100):
    radial_dist=dist_result[:,k+1]
    if radial_dist.mean()>85:
        interface_list.append(k)
    elif radial_dist.mean()<75:
        core_list.append(k)
np.savetxt('interface_list.xvg',interface_list,delimiter='	')
np.savetxt('core_list.xvg',core_list,delimiter='	')

fig, ax = plt.subplots(figsize=(8,6))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
tick_prop = font_manager.FontProperties(fname=font_path, size=16)
legend_prop = font_manager.FontProperties(fname=font_path, size=15) 
cm = plt.cm.get_cmap('tab20')
for i in range(100):
    color_index=i%20
    ax.plot(dist_result[:,0]/1000,dist_result[:,i+1]/10,alpha=1,color=cm.colors[color_index])
ax.hlines(y=7.5,xmin=2000,xmax=30000,color='gray')
ax.hlines(y=8.5,xmin=2000,xmax=30000,color='gray')
plt.xlim(2000,30000)
ax.set_xlabel('Time (ns)',fontproperties=font_prop)
ax.set_ylabel('Radial Distance (nm)',fontproperties=font_prop)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)
#plt.legend(prop=legend_prop)
plt.savefig('Droplet_foldD_distr.png',dpi=600,bbox_inches='tight')
plt.show()

    
    



