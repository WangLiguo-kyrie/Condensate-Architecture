#! /usr/bin/python

import sys
import os
import numpy as np

import pandas as pd

trajectory_dir='/grain/liguo/MYC/IDP-conf-change/condensate/7-M3IDP-dih7-SC15/PGL1/more_copy/'
copy=33
res=730
Natoms=1629
folded_s=1 #resid 1
folded_e=488 #resid 217

os.system('echo "del 2-14\ndel 0\nq" | gmx make_ndx -f {}/production2-Verlet.gro -o MSD_NtDD.ndx '.format(trajectory_dir))
for i in range(copy):
    folded_index_s=i*Natoms+folded_s
    folded_index_e=i*Natoms+folded_e
    os.system('echo "a {}-{}\nname {} NtDD{}\nq" | gmx make_ndx -f {}/production2-Verlet.gro -o MSD_NtDD.ndx -n MSD_NtDD.ndx'.format(folded_index_s,folded_index_e,str(i+1),str(i+1),trajectory_dir))


#Periodic boundary conditions were dealt by default in the MSD calculations by gmx msd versions before 2022 (for example Gromacs versions 2021.7 used here) with raw trajectory input.
#Compare the msd result of traj input with/without whole+nojump PBC processing, got the same results, so gmx msd should consider the pbc by default.
#new version 2023 gmx, -pbc -rmpbc could be set explicitly with "yes", -rmcomm option was removed. We could use gmx trjconv -fit translation to remove the center of mass motion first. 
for i in range(copy):
    os.system('echo "{}\n0" |gmx msd -f {}/production2-Verlet.xtc -s {}/production2-Verlet.tpr -n MSD_NtDD.ndx -o msd_NtDD{}.xvg -beginfit 0 -endfit 20000  -b 29500000 -rmcomm'.format(str(i+1),trajectory_dir,trajectory_dir,str(i+1)))
    



