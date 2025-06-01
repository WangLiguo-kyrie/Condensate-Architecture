#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
for the calculation of life time of interchain contact
"""

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
from numpy import mean
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
import seaborn as sns
import os,glob
import itertools
import pickle

#follow this seq
contact_mode=['IDR-IDR','IDR-Folded','Folded-Folded']
condensate_mean=[]
dimer_mean=[]
condensate_median=[]
dimer_median=[]

life_collect_d2d = pickle.load( open( "./contacts_life_collect_d2d.p", "rb" ))
condensate_mean.append(life_collect_d2d.mean())
condensate_median.append(np.median(life_collect_d2d))
life_collect_d2f = pickle.load( open( "./contacts_life_collect_d2f.p", "rb" ))
condensate_mean.append(life_collect_d2f.mean())
condensate_median.append(np.median(life_collect_d2f))
life_collect_f2f = pickle.load( open( "./contacts_life_collect_f2f.p", "rb" ))
condensate_mean.append(life_collect_f2f.mean())
condensate_median.append(np.median(life_collect_f2f))

HSPB2_d2d=pickle.load( open( "/grain/liguo/MYC/IDP-conf-change/Project2_Scaffold/monomer_conformation/Singlechain_sim/HSPB2/HSPB2_dimer/contacts_life_collect_d2d.p", "rb" ))
dimer_mean.append(HSPB2_d2d.mean())
dimer_median.append(np.median(HSPB2_d2d))
HSPB2_d2f=pickle.load( open( "/grain/liguo/MYC/IDP-conf-change/Project2_Scaffold/monomer_conformation/Singlechain_sim/HSPB2/HSPB2_dimer/contacts_life_collect_d2f.p", "rb" ))
dimer_mean.append(HSPB2_d2f.mean())
dimer_median.append(np.median(HSPB2_d2f))
HSPB2_f2f=pickle.load( open( "/grain/liguo/MYC/IDP-conf-change/Project2_Scaffold/monomer_conformation/Singlechain_sim/HSPB2/HSPB2_dimer/contacts_life_collect_f2f.p", "rb" ))
dimer_mean.append(HSPB2_f2f.mean())
dimer_median.append(np.median(HSPB2_f2f))

print(condensate_mean)
print(dimer_mean)
print(condensate_median)
print(dimer_median)

mean_incre=np.array(condensate_mean)/np.array(dimer_mean)
median_incre=np.array(condensate_median)/np.array(dimer_median)

font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=18)
tick_prop = font_manager.FontProperties(fname=font_path, size=15)
legend_prop = font_manager.FontProperties(fname=font_path, size=13)  
fig, ax = plt.subplots(figsize=(6,4.7))
ax2 = ax.twinx()
width = 0.3  # the width of the bars
pos=np.arange(len(condensate_mean))

bar1=ax.bar(pos, [round(x,1) for x in dimer_mean], width, label='HSPB2 Dimer',color='#a5a7ab',alpha=0.8,log=True)
ax.bar_label(bar1, padding=3,fontproperties=tick_prop)
bar2=ax.bar(pos + width, [round(x,1) for x in condensate_mean], width, label='HSPB2 Condensate',color='#3d5a9c',alpha=0.7,log=True)
ax.bar_label(bar2, padding=3,fontproperties=tick_prop)
ax.legend(prop=legend_prop)
ax2.plot(pos+width/2,mean_incre,linewidth=1.5,color='#b03d56',marker='D')
   
ax.set_ylabel('Mean contact life time (ns)',fontproperties=font_prop)
ax2.set_ylabel('Increment',fontproperties=font_prop)
ax.set_xticks(pos+width/2,labels=contact_mode,fontproperties=tick_prop)
ax.set_yticks([1,10,100,1000],labels=[1,10,100,1000],fontproperties=tick_prop)
ax2.set_yticks([2,4,6,8,10],labels=[2,4,6,8,10],fontproperties=tick_prop)   
ax.set_ylim(1,6000)
plt.savefig('lifetime_mode_mean.png',dpi=600,bbox_inches='tight')
plt.show()


fig, ax = plt.subplots(figsize=(6,4.7))
ax2 = ax.twinx()
width = 0.3  # the width of the bars
pos=np.arange(len(condensate_median))
bar1=ax.bar(pos, [round(x,1) for x in dimer_median], width, label='HSPB2 Dimer',color='#a5a7ab',alpha=0.8,log=True)
ax.bar_label(bar1, padding=3,fontproperties=tick_prop)
bar2=ax.bar(pos + width, [round(x,1) for x in condensate_median], width, label='HSPB2 Condensate',color='#3d5a9c',alpha=0.7,log=True)
ax.bar_label(bar2, padding=3,fontproperties=tick_prop)
ax.legend(prop=legend_prop)
ax2.plot(pos+width/2,median_incre,linewidth=1.5,color='#b03d56',marker='D')
   
ax.set_ylabel('Median contact life time (ns)',fontproperties=font_prop)
ax2.set_ylabel('Increment',fontproperties=font_prop)
ax.set_xticks(pos+width/2,labels=contact_mode,fontproperties=tick_prop)
ax.set_yticks([1,10,100],labels=[1,10,100],fontproperties=tick_prop)
ax2.set_yticks([2,4,6,8],labels=[2,4,6,8],fontproperties=tick_prop)   
ax.set_ylim(1,500)
plt.savefig('lifetime_mode_median.png',dpi=600,bbox_inches='tight')
plt.show()
