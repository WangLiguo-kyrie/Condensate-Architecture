#! /usr/bin/python

import sys
import os
import numpy as np

import pandas as pd

trajectory_dir='/grain/liguo/MYC/IDP-conf-change/condensate/7-M3IDP-dih7-SC15/HSPB2/more_copy/'
copy=100
res=182
Natoms=438
folded_s=155 #resid 65
folded_e=361 #resid 147

os.system('echo "del 2-14\ndel 0\nq" | gmx make_ndx -f {}/production2-Verlet.gro -o MSD_folded.ndx '.format(trajectory_dir))
for i in range(copy):
    folded_index_s=i*Natoms+folded_s
    folded_index_e=i*Natoms+folded_e
    os.system('echo "a {}-{}\nname {} folded{}\nq" | gmx make_ndx -f {}/production2-Verlet.gro -o MSD_folded.ndx -n MSD_folded.ndx'.format(folded_index_s,folded_index_e,str(i+1),str(i+1),trajectory_dir))

#directly use gromacs trajectory without trjconv processing.
#Compare the msd result with/without whole+nojump processing, got the same results, so gmx msd should consider the pbc by default.
# in new version 2024.2 gmx, -pbc -rmpbc could be set explicitly with "yes" bydefault, the msd result is same to 2021.7 gmx msd (without -rmcomm). But this -rmcomm option is not implemeted in 2024.4
#add rmcomm to deal with center of mass motion of condnesate  
for i in range(copy):
#for i in [98,99]:
    os.system('echo "{}\n0" |gmx msd -f {}/production2-Verlet.xtc -s {}/production2-Verlet.tpr -n MSD_folded.ndx -o MSD_folded{}.xvg -beginfit 0 -endfit 20000  -b 29500000 -rmcomm'.format(str(i+1),trajectory_dir,trajectory_dir,str(i+1)))
    



