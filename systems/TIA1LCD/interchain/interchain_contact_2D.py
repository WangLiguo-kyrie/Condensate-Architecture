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
copy=48

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
    frag_list=list(range(0,copy))
    # select other protein copy by fragments
    frag_list.remove(i)
    print(frag_list)
    
    print('***** Protein Monomer {} *****'.format(i))       
    inter_matrix_collect=np.zeros((n_res, n_res))
    timeseries = []
    
    monomer_A=frag.select_atoms('name BB')
    print(monomer_A.atoms) 
    #iterately select other protein copys 
    for j in frag_list:
        monomer_B=fragments[j].select_atoms('name BB')
        print(monomer_B.atoms) 
        inter_matrix_monomer=np.zeros((n_res, n_res))

        #directly construct the contact matrix between monomers, and then add all matrix together
        if i in Dilute_monomer_list:
            dilute_interval=np.loadtxt(trajectory_dir+'phase_Exchange/Dilute_monoer{}.txt'.format(i))
            for index, ts in enumerate(u.trajectory[start:end:step]): 
                if ts.time not in dilute_interval[:,0]:
                    inter_matrix=contactsMatrix_within_cutoff(monomer_A, monomer_B,ts.dimensions)
                    inter_matrix_monomer=inter_matrix_monomer+inter_matrix
                    timeseries.append(ts.time)
                    print(' Time {}ps'.format(ts.time))
                        
        else:
            for index, ts in enumerate(u.trajectory[start:end:step]):    
                inter_matrix=contactsMatrix_within_cutoff(monomer_A, monomer_B,ts.dimensions)
                inter_matrix_monomer=inter_matrix_monomer+inter_matrix
                timeseries.append(ts.time)
                print(' Time {}ps'.format(ts.time))
            
        inter_matrix_collect=inter_matrix_collect+inter_matrix_monomer
    frames=len(timeseries)  
    print(frames)
    contacts_2D=inter_matrix_collect*(copy-1)/frames  
    contacts_strength=contacts_2D.sum()   
    print( contacts_strength)
    print(contacts_2D.shape)
    np.savetxt('2D/2D_interchain_contact_monomer{}.xvg'.format(i),contacts_2D,delimiter='	')

    


