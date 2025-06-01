#! /usr/bin/python

import sys
import os
import numpy as np

import math
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
from scipy.stats import pearsonr,spearmanr
from scipy.stats import entropy
from scipy.spatial import distance

os.chdir('2D')
n_res=97
inter_matrix_collect=np.zeros((n_res, n_res))
for i in range(48):
    data=np.loadtxt('2D_interchain_contact_monomer{}.xvg'.format(i)) 
    inter_matrix_collect=inter_matrix_collect+data
inter_matrix_condensate=inter_matrix_collect/48
flat_condensate=inter_matrix_condensate.reshape(-1,1)

data1=np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project2_Scaffold/monomer_conformation/Singlechain_sim/TIA1/TIA1_dimer/replica1/2D_interchain_contact.xvg')
data2=np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project2_Scaffold/monomer_conformation/Singlechain_sim/TIA1/TIA1_dimer/replica2/2D_interchain_contact.xvg')
data3=np.loadtxt('/grain/liguo/MYC/IDP-conf-change/Project2_Scaffold/monomer_conformation/Singlechain_sim/TIA1/TIA1_dimer/replica3/2D_interchain_contact.xvg')
data_dilute=(data1+data2+data3)/3
data_dilute_ave=(data_dilute+data_dilute.T)/2
flat_dilute=data_dilute_ave.reshape(-1,1)

KL=entropy(flat_condensate,flat_dilute)
print(KL)
JS=distance.jensenshannon(flat_condensate,flat_dilute)
print(JS**2)
with open('2D_interchain_KL.txt','w') as g:
    print('KL divergence: {}'.format(KL[0]),file=g)
    print('JS divergence: {}'.format(JS**2),file=g) 
