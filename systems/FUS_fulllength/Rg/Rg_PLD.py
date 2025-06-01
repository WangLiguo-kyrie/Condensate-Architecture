#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import MDAnalysis as mda
from MDAnalysis.analysis import contacts
import numpy as np
import matplotlib.pyplot as plt
import os


#trajectory pbc and protein extraction
trajectory_dir='/grain/liguo/MYC/IDP-conf-change/condensate/7-M3IDP-dih7-SC15/Folded_with_IDR/FUS/cubic2/'
#os.system('echo "1\n 1\nq" |gmx trjconv -f {}production.xtc -s {}production.tpr -pbc mol -center -o {}production-pbc-protein.xtc'.format(trajectory_dir,trajectory_dir,trajectory_dir))
#os.system('echo "1\n 1\nq" |gmx trjconv -f {}production.gro -s {}production.tpr -pbc mol -center -o {}production-pbc-protein.gro'.format(trajectory_dir,trajectory_dir,trajectory_dir))
#os.system('gmx grompp -f /grain/liguo/MYC/IDP-conf-change/condensate/martini_nvt.mdp -c {}production-pbc-protein.gro -p {}topol-protein.top -o {}production-protein.tpr'.format(trajectory_dir,trajectory_dir,trajectory_dir))

u = mda.Universe(trajectory_dir+'production-protein.tpr',trajectory_dir+'production-pbc-protein.xtc')
ag = u.select_atoms('name BB')
print(len(ag.atoms))
fragments = ag.atoms.fragments

monomer_length=526


#discard initial 4us 
start = 20000
end = -1
#every 1ns
step = 5

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
    PLD_start=i*monomer_length+1
    PLD_end=i*monomer_length+163
    print(PLD_start)
    print(PLD_end)
    PLD=frag.select_atoms('resid {}-{}'.format(PLD_start,PLD_end))
    print(len(PLD.atoms))
      
    if i in Dilute_monomer_list:
        dilute_interval=np.loadtxt(trajectory_dir+'phase_Exchange/Dilute_monoer{}.txt'.format(i))
        for index, ts in enumerate(u.trajectory[start:end:step]): 
            if ts.time in dilute_interval[:,0]:
                print('Dilute Protein Monomer {}, Time {}, Rg {}'.format(i,ts.time,PLD.radius_of_gyration()))
                dilute_fragment_result.append([i,ts.time,ts.frame,PLD.radius_of_gyration()])
            else:
                print('Condensate Protein Monomer {}, Time {}, Rg {}'.format(i,ts.time,PLD.radius_of_gyration()))
                condensate_fragment_result.append([i,ts.time,ts.frame,PLD.radius_of_gyration()])
    else:   
        for index, ts in enumerate(u.trajectory[start:end:step]):    
            print('Condensate Protein Monomer {}, Time {}, Rg {}'.format(i,ts.time,PLD.radius_of_gyration()))
            condensate_fragment_result.append([i,ts.time,ts.frame,PLD.radius_of_gyration()])                 

    np.savetxt('PLD_Condensate_phase_monomer{}.xvg'.format(i),condensate_fragment_result,delimiter='	')
    if len(dilute_fragment_result)>0:
        np.savetxt('PLD_Dilute_phase_monomer{}.xvg'.format(i),dilute_fragment_result,delimiter='	')     

     


