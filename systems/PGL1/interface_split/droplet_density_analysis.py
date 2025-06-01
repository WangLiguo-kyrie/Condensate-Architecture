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

trajectory_dir='/grain/liguo/MYC/IDP-conf-change/condensate/7-M3IDP-dih7-SC15/PGL1/more_copy/'

u = mda.Universe(trajectory_dir+'production-protein.tpr',trajectory_dir+'production2-Verlet-pbc-protein.xtc')
condensate=u.select_atoms('all')

n_res=730
copy=33

#discard initial 2us 
start = 4000
end = -1
#every 1ns
step = 2

hist_collect=np.zeros(79)
density=np.zeros(79)
for index, ts in enumerate(u.trajectory[start:end:step]): 
    dist = contacts.distance_array(condensate.center_of_geometry(), condensate.positions)
    #bindwidth=0.3nm
    hist,edge=np.histogram(dist,bins=np.arange(0,240,3))
    hist_collect=hist_collect+hist
frames=len(u.trajectory[start:end:step])
hist_ave=hist_collect/frames

bins=[int((edge[i]+edge[i+1])/2) for i in range(79)]
print(bins)
bins=np.array(bins)
density=[hist_ave[i]/(4*np.pi*bins[i]**2) for i in range(79)]
np.savetxt('droplet_density.xvg',density,delimiter='	')


fig, ax = plt.subplots(figsize=(8,6))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
tick_prop = font_manager.FontProperties(fname=font_path, size=16)
legend_prop = font_manager.FontProperties(fname=font_path, size=15)  
ax.plot(bins[8:]/10,density[8:],label='Protein',alpha=1,color='green')
ax.set_ylabel('Radial Bead Density',fontproperties=font_prop)
ax.set_xlabel('Radial Bead Distance (nm)',fontproperties=font_prop)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)
#plt.legend(prop=legend_prop)
plt.savefig('Droplet_density.png',dpi=600,bbox_inches='tight')
plt.show()

    
    



