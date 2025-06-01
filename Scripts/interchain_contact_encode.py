#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
encode contact state of each residue pair
and calculate the contact time of each association event
"""

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import os
import itertools
import pickle

    
def contactsMatrix_within_cutoff(group_a, group_b,box_dim, radius1=8,radius2=12):
    '''
    dual BB distance cutoff 0.8/1.2nm for a robust creteria of association/dissociation process 
    '''
    # calculate distances between group_a and group_b, matrix shape(a.len,b.len)
    dist = contacts.distance_array(group_a.positions, group_b.positions,box=box_dim)
    # value2: <0.8nm;  value1: 0.8< <1.2nm;  value0: >1.2nm 
    contacts_matrix1 = contacts.contact_matrix(dist, radius1)*1
    contacts_matrix2 = contacts.contact_matrix(dist, radius2)*1
    contacts_matrix=contacts_matrix1+contacts_matrix2
    contacts_matrix=contacts_matrix.reshape(contacts_matrix.shape[0],contacts_matrix.shape[1],1)
    return contacts_matrix  
    
def contact_life(data, stepsize=1):
    
    '''
    split contact seq by value 0 to generate the subtraj in contact, no value=0 inside subtraj
    and then calculate the time length from first value 2 (i.e. become bound state) to final frame in subtraj (where value become 0, meaning totally dissociation)
    '''
    contacts_traj=[]
    contacts_life=[]
    time_unit=0.5   #0.5ns each frame
    #group consecutive elements of the same value 0 to generate subtraj withou value 0 (unbound) inside
    for key, group in itertools.groupby(data,lambda x: x==0):
        if key==False:
            contacts_traj.append(list(group))
    
    for subtraj in contacts_traj:
        #calculate from first value 2
        if 2 in subtraj:
            idxs= [i for i, v in enumerate(subtraj) if v == 2]
            subtraj_clean=subtraj[idxs[0]:]
            print(subtraj_clean)
            subtraj_time=len(subtraj_clean)*time_unit*stepsize
            contacts_life.append(subtraj_time)
    
    return contacts_life
    
           
#trajectory pbc and protein extraction
trajectory_dir='/grain/liguo/MYC/IDP-conf-change/condensate/7-M3IDP-dih7-SC15/FUSLCD-lowconc2/'

u = mda.Universe(trajectory_dir+'production-protein.tpr',trajectory_dir+'production-pbc-protein.xtc')
ag = u.select_atoms('name BB')
print(len(ag.atoms))
fragments = ag.atoms.fragments

n_res=163
copy=36

#discard initial 4us 
start = 8000
end = -1
#every 0.5ns
step = 1

dilutes=np.loadtxt(trajectory_dir+'phase_Exchange/Phase_Exchange.txt')
Dilute_monomer_list=np.unique(dilutes[:,1])
print(Dilute_monomer_list)

#os.mkdir('contacts_life')
#loop over each fragment
for i, frag in enumerate(fragments):
    #simply remove the dilute-condensate exchange monomers here
    if ((2<i<8) & (i not in Dilute_monomer_list)):
        #remove reverse duplicates
        frag_list=list(range(i+1,copy))
        # select other protein copy by fragments
        print(frag_list)
    
        print('***** Protein Monomer {} *****'.format(i))         
        monomer_A=frag.select_atoms('name BB ')
        print(monomer_A.atoms) 
        #iterately select other protein copys 
        for j in frag_list:
            monomer_B=fragments[j].select_atoms('name BB')
            print(monomer_B.atoms) 

            #directly construct the contact matrix between monomers, and then collect each frame matrix 
            inter_matrix_monomer=np.zeros((n_res, n_res,1))
            #save every 250ns-block to accelerate 
            inter_matrix_block=np.zeros((n_res, n_res,1))            
            for index, ts in enumerate(u.trajectory[start:end:step]):    
                inter_matrix=contactsMatrix_within_cutoff(monomer_A, monomer_B,ts.dimensions)
                inter_matrix_block=np.concatenate((inter_matrix_block,inter_matrix),axis=2)
                if (ts.time % 250000)==0:
                    print(' Time {}ns'.format(ts.time/1000))
                    inter_matrix_monomer=np.concatenate((inter_matrix_monomer,inter_matrix_block[:,:,1:]),axis=2)
                    print(inter_matrix_monomer.shape)
                    inter_matrix_block=np.zeros((n_res, n_res,1))
            print(inter_matrix_monomer.shape)
            
            dic={}
            #filter the existing contacts during trajectory
            for res_i in range(n_res):
                for res_j in range(n_res):
                    # the contact state evolution of each pair
                    contact_seq=inter_matrix_monomer[res_i,res_j,:]
                    #print(contact_seq)
                
                    #only save the contacts existing with value 2
                    if 2 in contact_seq:
                        lifetime= contact_life(contact_seq)
                        contact_name='{}-{}'.format(res_i+1,res_j+1)
                        dic_contact={contact_name:lifetime}
                        print(dic_contact)
                        dic.update(dic_contact)

                    else:
                        print('{}_{}-{}_{} contact does not exist'.format(i,res_i+1,j,res_j+1))
            print(dic)
            pickle.dump(dic, open( "contacts_life/monomer{}_monomer{}.p".format(i,j), "wb" ))
                


    


