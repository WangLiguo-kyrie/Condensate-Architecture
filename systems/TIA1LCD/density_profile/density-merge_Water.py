#! /usr/bin/python
# -*- coding: utf-8 -*-
# @Date    : 2018-01-30 14:05:51
# @Author  : WDD 
# @Link    : https://github.com/dongdawn
# @Version : v1
import sys
import os
import numpy as np

import pandas as pd


intervals=[0,4000,8000,12000,16000,20000,24000,28000,32000,36000]
data=[]
for i in range(9):
    former=intervals[i]
    latter=intervals[i+1]

    filename='density-{}_{}ns.xvg'.format(former,latter)
    print(filename)
    os.system("sed -i 's/^@/#/g' %s " %filename)
    data_interval = np.loadtxt(filename)
    ax_z = data_interval[:,0]
    ligands_density = np.reshape(data_interval[:,2],(-1,1))
    print(ligands_density.shape)
    if i==0:
        data=data_interval[:,(0,2)]
    else: data=np.concatenate((data,ligands_density),axis=1)
np.savetxt('density_Water.xvg',data,delimiter='	')


