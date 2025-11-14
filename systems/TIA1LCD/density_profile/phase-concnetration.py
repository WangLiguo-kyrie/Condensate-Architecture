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


intervals=[8000,12000,16000,20000,24000,28000,32000,36000]
factor=10903/12492
data=[]
with open('phase.txt','w') as g:
    print('Begin\tEnd\tDiluteconcentration\tDenseconcentration',file=g)
    for i in range(7):
        former=intervals[i]
        latter=intervals[i+1]
        print('density interval: {} to {}'.format(former,latter))

        filename='density-{}_{}ns.xvg'.format(former,latter)
        print(filename)
        os.system("sed -i 's/^@/#/g' %s " %filename)
        data_interval = np.loadtxt(filename)
        print(data_interval.shape)
        #Dilute phase
        dilute1=data_interval[(data_interval[:,0]>20)]
        dilute2=data_interval[(data_interval[:,0]<-20)]
        dilute_conc=(np.mean(dilute1[:,1])+np.mean(dilute2[:,1]))*factor/2
        #Dense phase
        dense=data_interval[((data_interval[:,0]>-7)&(data_interval[:,0]<7))]*factor
        print(dense)
        dense_conc=(np.mean(dense[:,1]))
        print('{}\t{}\t{}\t{}'.format(former,latter,dilute_conc,dense_conc),file=g)



