#! /usr/bin/python

import sys
import os
import numpy as np

import pandas as pd

trajectory_dir='/grain/liguo/MYC/IDP-conf-change/condensate/7-M3IDP-dih7-SC15/PGL1/more_copy/'
copy=33
res=730
Natoms=1629
folded_s=501 #resid 223
folded_e=1032 #resid 449

os.system('echo "del 2-14\ndel 0\nq" | gmx make_ndx -f {}/production2-Verlet.gro -o MSD_CDD.ndx '.format(trajectory_dir))
for i in range(copy):
    folded_index_s=i*Natoms+folded_s
    folded_index_e=i*Natoms+folded_e
    os.system('echo "a {}-{}\nname {} CDD{}\nq" | gmx make_ndx -f {}/production2-Verlet.gro -o MSD_CDD.ndx -n MSD_CDD.ndx'.format(folded_index_s,folded_index_e,str(i+1),str(i+1),trajectory_dir))


#directly use gromacs trajectory without trjconv processing.
#Compare the msd result with/without whole+nojump processing, got the same results, so gmx msd should consider the pbc by default.
# in new version 2024.2 gmx, -pbc -rmpbc could be set explicitly with "yes" bydefault, the msd result is same to 2021.7 gmx msd (without -rmcomm). But this -rmcomm option is not implemeted in 2024.4
#add rmcomm to deal with center of mass motion of condnesate  
for i in range(copy):
    os.system('echo "{}\n0" |gmx msd -f {}/production2-Verlet.xtc -s {}/production2-Verlet.tpr -n MSD_CDD.ndx -o msd_CDD{}.xvg -beginfit 0 -endfit 20000  -b 29500000 -rmcomm'.format(str(i+1),trajectory_dir,trajectory_dir,str(i+1)))
    



