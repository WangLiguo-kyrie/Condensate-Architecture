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

# Nter 1:64 , ACD: 65-147, Cter:148-182
folded_domain=list(range(65,148))
print(folded_domain)
IDR_domain=list(range(1,65))+list(range(148,183))
print(IDR_domain)

#life_files=glob.glob('contacts_life/*.p')
#life_collect_f2f=[]
#life_collect_d2d=[]
#life_collect_d2f=[]
#life_collect_f2f_frag=[]
#life_collect_d2d_frag=[]
#life_collect_d2f_frag=[]
#for i, life_file in enumerate(life_files):
#    index=i+1
#    contact_dic = pickle.load( open( life_file, "rb" ))
#    print(life_file)
#    print('Process: index {}/{}'.format(index,len(life_files)))
#    for keys, values in contact_dic.items():
#        print(keys)
#        positions=keys.split('-')
#        #print(positions)
#        #print(positions[0])
        
#        if ((int(positions[0]) in folded_domain) &(int(positions[1]) in folded_domain)):
#            life_collect_f2f_frag=life_collect_f2f_frag+values
#        elif ((int(positions[0]) in IDR_domain) &(int(positions[1]) in IDR_domain)):
#            life_collect_d2d_frag=life_collect_d2d_frag+values
#        else: life_collect_d2f_frag=life_collect_d2f_frag+values
#    if (index%50)==0:
#        print(life_collect_f2f_frag)
#        print(life_collect_d2d_frag)
#        print(life_collect_d2f_frag)
#        life_collect_f2f=life_collect_f2f+life_collect_f2f_frag
#        life_collect_d2d=life_collect_d2d+life_collect_d2d_frag
#        life_collect_d2f=life_collect_d2f+life_collect_d2f_frag
#        life_collect_f2f_frag=[]
#        life_collect_d2d_frag=[]
#        life_collect_d2f_frag=[]
#        print(len(life_collect_f2f))
#        print(len(life_collect_d2d))
#        print(len(life_collect_d2f))
               
#print(life_collect_f2f)
#print(life_collect_d2d)
#print(life_collect_d2f)
#life_collect_f2f=np.array(life_collect_f2f)
#pickle.dump(life_collect_f2f, open( "./contacts_life_collect_f2f.p", "wb" ))
#life_collect_d2d=np.array(life_collect_d2d)
#pickle.dump(life_collect_d2d, open( "./contacts_life_collect_d2d.p", "wb" ))
#life_collect_d2f=np.array(life_collect_d2f)
#pickle.dump(life_collect_d2f, open( "./contacts_life_collect_d2f.p", "wb" ))

life_collect_f2f = pickle.load( open( "./contacts_life_collect_f2f.p", "rb" ))
life_collect_d2d = pickle.load( open( "./contacts_life_collect_d2d.p", "rb" ))
life_collect_d2f = pickle.load( open( "./contacts_life_collect_d2f.p", "rb" ))
#life_collect = pickle.load( open( "./contacts_life_collect.p", "rb" ))

font_path = '/grain/liguo/biocondensates/PMF-test/pulldim-YYY/Arial.ttf'
font_prop = font_manager.FontProperties(fname=font_path, size=18)
tick_prop = font_manager.FontProperties(fname=font_path, size=15)
legend_prop = font_manager.FontProperties(fname=font_path, size=13)  


data_folded_dimer=pickle.load( open( "/grain/liguo/MYC/IDP-conf-change/Project2_Scaffold/monomer_conformation/Singlechain_sim/HSPB2/foldD_dimer/contacts_life_stepsize400ps/contacts_life_collect.p", "rb" ))
data_HSPB2_dimer=pickle.load( open( "/grain/liguo/MYC/IDP-conf-change/Project2_Scaffold/monomer_conformation/Singlechain_sim/HSPB2/HSPB2_dimer/contacts_life_collect_f2f.p", "rb" ))
fig, ax = plt.subplots(figsize=(6,4.7))

#sns.violinplot([life_collect_f2f,data_HSPB2_dimer,data_folded_dimer],palette="Accent",inner_kws=dict(box_width=10, whis_width=1.5),log_scale=True)
sns.violinplot([life_collect_f2f,data_HSPB2_dimer],palette="Accent",inner_kws=dict(box_width=10, whis_width=1.5),log_scale=True)

ax.spines[['right', 'top']].set_visible(False)       
ax.set_ylabel('Contact life time (ns)',fontproperties=font_prop)
ax.set_xlabel('HSPB2 Folded-Folded contact',fontproperties=font_prop)
#ax.set_yscale('log')
#plt.xticks(ticks=[0,1,2],labels=['Condensate','Dimer','Folded ACD Dimer'],fontproperties=tick_prop)
plt.xticks(ticks=[0,1],labels=['Condensate','Dimer'],fontproperties=tick_prop)
plt.yticks(fontproperties=tick_prop)    
#plt.legend(prop=legend_prop)
plt.savefig('contact_life_comp_f2f_violin.png',dpi=600,bbox_inches='tight')
plt.show()
