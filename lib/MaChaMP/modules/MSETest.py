# -*- coding: utf-8 -*-
"""
Created on Sat Apr 15 11:31:53 2017

@author: Sir Thomas
"""

import numpy as np


def MSE_test_loc(seq,loc,lookup='none',offset=2,p='no',Reg=10):
    T_1 = []
    T_2 = []
    
    intervals = [0] + loc + [len(seq)]
    
    for left,right in zip(intervals[:-1],intervals[1:]):
        T_1.append(np.var(seq[left:right-1])*(right-left))
        T_2.append(np.mean(seq[left:right-1]))
    
    theta = sum([np.abs(i-j) for i,j in zip(T_2[:-1],T_2[1:])])
    
    if  theta > Reg:
        return 10000
    else:
        return sum(T_1)/len(seq)
