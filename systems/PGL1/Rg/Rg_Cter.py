#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import os

    
def contacts_within_cutoff(group_a, group_b, radius=6):
    '''
    determine if one protein monomer is outside condensate
    '''
    # calculate distances between group_a and group_b
    dist = contacts.distance_array(group_a.positions, group_b.positions)
    # determine which distances <= radius
    n_contacts = contacts.contact_matrix(dist, radius).sum()
    return n_contacts    

#trajectory pbc and protein extraction
trajectory_dir='/grain/liguo/MYC/IDP-conf-change/condensate/7-M3IDP-dih7-SC15/PGL1/more_copy/'

u = mda.Universe(trajectory_dir+'production-protein.tpr',trajectory_dir+'production2-Verlet-pbc-protein.xtc')
ag = u.select_atoms('name BB')
print(len(ag.atoms))
fragments = ag.atoms.fragments

monomer_length=730

#discard initial 2us 
start = 4000
end = -1
#every 1ns
step = 2

results = {}

#dilutes=np.loadtxt(trajectory_dir+'phase_Exchange/Phase_Exchange.txt')
#Dilute_monomer_list=np.unique(dilutes[:,1])
#print(Dilute_monomer_list)
Dilute_monomer_list=[]

#loop over each fragment
for i, frag in enumerate(fragments): 
    print('***** Protein Monomer {} *****'.format(i))       
    dilute_fragment_result = []
    condensate_fragment_result = []
    #res 450-730
    Cter_start=i*monomer_length+450
    Cter_end=i*monomer_length+730
    print(Cter_start)
    print(Cter_end)
    Cter=frag.select_atoms('resid {}-{}'.format(Cter_start,Cter_end))
    print(len(Cter.atoms))
      
    if i in Dilute_monomer_list:
        dilute_interval=np.loadtxt(trajectory_dir+'phase_Exchange/Dilute_monoer{}.txt'.format(i))
        for index, ts in enumerate(u.trajectory[start:end:step]): 
            if ts.time in dilute_interval[:,0]:
                print('Dilute Protein Monomer {}, Time {}, Rg {}'.format(i,ts.time,Cter.radius_of_gyration()))
                dilute_fragment_result.append([i,ts.time,ts.frame,Cter.radius_of_gyration()])
            else:
                print('Condensate Protein Monomer {}, Time {}, Rg {}'.format(i,ts.time,Cter.radius_of_gyration()))
                condensate_fragment_result.append([i,ts.time,ts.frame,Cter.radius_of_gyration()])
    else:   
        for index, ts in enumerate(u.trajectory[start:end:step]):    
            print('Condensate Protein Monomer {}, Time {}, Rg {}'.format(i,ts.time,Cter.radius_of_gyration()))
            condensate_fragment_result.append([i,ts.time,ts.frame,Cter.radius_of_gyration()])                 

    np.savetxt('Cter_Condensate_phase_monomer{}.xvg'.format(i),condensate_fragment_result,delimiter='	')
    if len(dilute_fragment_result)>0:
        np.savetxt('Cter_Dilute_phase_monomer{}.xvg'.format(i),dilute_fragment_result,delimiter='	')     

     


