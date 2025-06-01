#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import os

    
def contacts_within_cutoff(group_a, group_b,box_dim, radius=10):
    '''
    BB distance cutoff 1nm instead of 0.6nm distance for whole residue 
    '''
    # calculate distances between group_a and group_b
    dist = contacts.distance_array(group_a.positions, group_b.positions,box=box_dim)
    # determine which distances <= radius
    n_contacts = contacts.contact_matrix(dist, radius).sum()
    return n_contacts  

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
    results_contact=[]
    for k in range(1,n_res+1):
        res_A=k	
        # neighboring i-3,i-2,i-1,i, i+1,i+2,i+3 residues
        neighbor=[x for x in range(res_A-3,res_A+4)]
        neighbor_filter=[x for x in neighbor if (98> x >0)]
        
        resA=frag.select_atoms('resid {} and name BB'.format(res_A+i*n_res))
        neighbor_list=''
        for j in neighbor_filter:
            neighbor_list=neighbor_list+' '+str(j++i*n_res)
        print(neighbor_list)
        resB=frag.select_atoms('not resid {} and name BB '.format(neighbor_list))
        print(resA.atoms)
        print(resB.atoms)
        timeseries = []
      
        if i in Dilute_monomer_list:
            dilute_interval=np.loadtxt(trajectory_dir+'phase_Exchange/Dilute_monoer{}.txt'.format(i))
            for index, ts in enumerate(u.trajectory[start:end:step]): 
                if ts.time not in dilute_interval[:,0]:
                    ca=contacts_within_cutoff(resA, resB,ts.dimensions)
                    timeseries.append([ts.time, ca])
                    print('Res {}, Time {}ps, Contacts {}'.format(res_A,ts.time,ca))

        else:   
            for index, ts in enumerate(u.trajectory[start:end:step]):    
                ca=contacts_within_cutoff(resA, resB,ts.dimensions)
                timeseries.append([ts.time, ca])
                print('Res {}, Time {}ps, Contacts {}'.format(res_A,ts.time,ca)) 
        timeseries=np.array(timeseries)

        ca_contact=timeseries[:,1].sum()  
        frames=timeseries.shape[0] 
        print(frames)
        contacts_strength=ca_contact/frames
        results_contact.append([res_A,contacts_strength])    
    results_contact=np.array(results_contact)    
    np.savetxt('intrachain_contact_monomer{}.xvg'.format(i),results_contact,delimiter='	')
  

    


