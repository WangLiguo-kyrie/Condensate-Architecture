#! /usr/bin/python
# -*- coding: utf-8 -*-
# calculate gyrate distribution in different trajectory interval to analyze convergence and conformation ensemble
import sys
import os
import numpy as np
import argparse
import pandas as pd
import seaborn as sns
#sns.set(color_codes=True)
import matplotlib

import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from statistics import mean
from scipy.stats import pearsonr,spearmanr

#for i in range(36):    
#    inter = np.loadtxt('../interchain_contact/interchain_contact_monomer{}.xvg'.format(i))
#    intra = np.loadtxt('../intrachain_contact/intrachain_contact_monomer{}.xvg'.format(i))
#    inter_contact=inter[:,1]
#    intra_contact=intra[:,1]
#    inter_intra=inter
#    inter_intra[:,1]=inter_contact+intra_contact
#    np.savetxt('interintra_contact_monomer{}.xvg'.format(i),inter_intra,delimiter='	') 


fig, ax = plt.subplots(figsize=(10,6))
font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=20)
tick_prop = font_manager.FontProperties(fname=font_path, size=16)
legend_prop = font_manager.FontProperties(fname=font_path, size=16) 
cm = plt.cm.get_cmap('tab20')
contact_collect=[]
for i in range(36):    
    data = np.loadtxt('interintra_contact_monomer{}.xvg'.format(i))
    contact = data[:,1]
    resid = data[:,0]      
    if i==0:contact_collect=contact
    else:contact_collect=np.vstack((contact_collect,contact))
        
    color_index=i%20
    ax.plot(resid,contact, color=cm.colors[color_index],alpha=0.5) 
print(contact_collect.shape)  
mean_value=contact_collect.mean(axis=0)
std_value=contact_collect.std(axis=0)
ax.plot(resid,mean_value, color='#3d5a9c', linewidth=3,label='Condensate')
ax.errorbar(resid,mean_value,yerr=std_value, fmt='-',elinewidth=1)
print('contact number per residue: {}'.format(mean(mean_value)))

singlechain= np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project2_Scaffold/monomer_conformation/Singlechain_sim/FUSLCD/intrachain_contact_1D.xvg')
contact_single = singlechain[:,1]
resid_single = singlechain[:,0] 
ax.plot(resid_single,contact_single, color='black', linewidth=3,label='Dilute') 

corr_pearson, _ = pearsonr(mean_value,contact_single)
corr_spearman, _ = spearmanr(mean_value,contact_single)
print(corr_pearson)
print(corr_spearman)
ax.text( 0.15,13,'Corr: {:.2f}'.format(corr_pearson),horizontalalignment='left',fontproperties=font_prop)

ax.spines[['right', 'top']].set_visible(False)
ax.set_ylabel('Contact Number',fontproperties=font_prop)
ax.set_xlabel('FUSLCD Residue',fontproperties=font_prop)
#plt.xlim(0.2,0.45)
plt.xticks(fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
plt.legend(prop=legend_prop)
plt.savefig('interchain_contact_1D.png',dpi=600,bbox_inches='tight')
plt.show()  





