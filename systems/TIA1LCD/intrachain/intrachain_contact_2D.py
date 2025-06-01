#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import os

    
def contactsMatrix_within_cutoff(group_a, group_b,box_dim, radius=10):
    '''
    BB distance cutoff 1nm instead of 0.6nm distance for whole residue 
    '''
    # calculate distances between group_a and group_b
    dist = contacts.distance_array(group_a.positions, group_b.positions,box=box_dim)
    # determine which distances <= radius
    contacts_matrix = contacts.contact_matrix(dist, radius)*1
    return contacts_matrix   

#trajectory pbc and protein extraction
trajectory_dir='/grain/liguo/MYC/IDP-conf-change/condensate/7-M3IDP-dih7-SC15/TIA1-lowconc1/'

u = mda.Universe(trajectory_dir+'production-protein.tpr',trajectory_dir+'production-pbc-protein.xtc')
ag = u.select_atoms('name BB')
print(len(ag.atoms))
fragments = ag.atoms.fragments

n_res=97

#discard initial 8us 
start = 16000
end = -1
#every 1ns
step = 2

dilutes=np.loadtxt(trajectory_dir+'phase_Exchange/Phase_Exchange.txt')
Dilute_monomer_list=np.unique(dilutes[:,1])
print(Dilute_monomer_list)

#loop over each fragment
for i, frag in enumerate(fragments): 
    print('***** Protein Monomer {} *****'.format(i))       
    intra_matrix_collect=np.zeros((n_res, n_res))
    timeseries = []
    BB=frag.select_atoms('name BB')
    print(BB.atoms)    
      
    if i in Dilute_monomer_list:
        dilute_interval=np.loadtxt(trajectory_dir+'phase_Exchange/Dilute_monoer{}.txt'.format(i))
        for index, ts in enumerate(u.trajectory[start:end:step]): 
            if ts.time not in dilute_interval[:,0]:
                intra_matrix=contactsMatrix_within_cutoff(BB, BB,ts.dimensions)
                intra_matrix_collect=intra_matrix_collect+intra_matrix
                timeseries.append(ts.time)
                print(' Time {}ps'.format(ts.time))

    else:   
        for index, ts in enumerate(u.trajectory[start:end:step]):    
            intra_matrix=contactsMatrix_within_cutoff(BB, BB,ts.dimensions)
            intra_matrix_collect=intra_matrix_collect+intra_matrix
            timeseries.append(ts.time)
            print(' Time {}ps'.format(ts.time))

    frames=len(timeseries)  
    print(frames)
    contacts_2D=intra_matrix_collect/frames  
    contacts_strength=contacts_2D.sum()   
    print( contacts_strength)
    print(contacts_2D.shape)
    np.savetxt('2D_intrachain_contact_monomer{}.xvg'.format(i),contacts_2D,delimiter='	')
  

    


