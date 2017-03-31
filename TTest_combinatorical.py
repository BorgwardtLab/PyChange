# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:59:19 2016

@author: Sir Thomas
"""


import numpy as np
import math
import bottleneck as bn
from collections import namedtuple
import scipy.stats as stats
import copy


def TTest_combinatorical(seq,no,offset):
    if (no+1)*3 > len(seq):
        return [0.]
    else:
        length = len(seq)
        loc = [(x+1)*3 for x in range(no)]
        best = [loc,_test_loc(seq,loc,offset)]
        while next(loc,length,no,offset):
            if loc[-1] > length:
                break
            p = _test_loc(seq,loc,offset)
            if best[1] > p:
                best = [copy.copy(loc),p]
        return best[0]


def next(loc,length,no,offset):
    i = -1
    while loc[i] >= length+offset*i:
        if abs(i) >= no:
            break
        else:
            i = i -1
    loc[i]=loc[i]+1
    for j in range(abs(i)-1):
        loc[i+j+1] = loc[i+j]+offset
    return 1

def _test_loc_p(seq,loc,offset):
    p_ttest = []
    p_out = []
    
    i = 0
    while i < len(loc): #intervals defined by loc. Take them and compare them seq[k],seq[l],seq[m]
        if i == 0:
            k = 0
        else:
            k = loc[i-1]
        l = loc[i]
        if i == len(loc)-1:
            m = len(seq)
        else:
            m = loc[i+1]
        
        if l-k < offset or m-l < offset:
            return 5.
        p_ttest.append(stats.ttest_ind(seq[k:l], seq[l:m], axis=0, equal_var=True)[1])
        
        i = i + 1
        
    return p_ttest

def _test_loc(seq,loc,offset):
    p_ttest = []
    p_out = []
    
    i = 0
    while i < len(loc): #intervals defined by loc. Take them and compare them seq[k],seq[l],seq[m]
        if i == 0:
            k = 0
        else:
            k = loc[i-1]
        l = loc[i]
        if i == len(loc)-1:
            m = len(seq)
        else:
            m = loc[i+1]
        
        if l-k < offset or m-l < offset:
            return 5.
        """ If taking only seond data point in the middle
        if k == 0:
            p_ttest.append(stats.ttest_ind(seq[k:l], seq[l+1:m:2], axis=0, equal_var=True)[1]) #take only every second data point
        elif m == len(seq):
            p_ttest.append(stats.ttest_ind(seq[k:l:2], seq[l:m], axis=0, equal_var=True)[1]) #take only every second data point
        else:
            p_ttest.append(stats.ttest_ind(seq[k:l:2], seq[l+1:m:2], axis=0, equal_var=True)[1]) #take only every second data point
        """
        p_ttest.append(stats.ttest_ind(seq[k:l], seq[l:m], axis=0, equal_var=True)[1])
        
        i = i + 1
        
    return stats.chisqprob(-2.*np.sum([np.log(x) for x in p_ttest]),2*len(loc)) 


