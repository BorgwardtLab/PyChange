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
import pandas as pd

def mc_no(seq,no):
    if no == 1:
        return [TTest(seq)[0]]
    else:
        cp = []
        p_vals = []
        gen_id = []
        dividereturn(seq,cp,p_vals,0,no+1,gen_id)
        backtracking(cp,p_vals,gen_id,no) #remove until best "first" changes are found.
        return sorted(cp)

def backtracking(cp,p_vals,gen_id,no):
    maxgen = np.max(gen_id)
    while len(cp) > no:   
        lowest_p = [i[0] for i in sorted(enumerate(p_vals), key=lambda x:-x[1])]
        gen_max_id = [i for i,x in enumerate(gen_id) if x == maxgen]
        lowest_p_gen = [i for i in lowest_p if i in gen_max_id ]
        while len(cp) > no and not len(lowest_p_gen) == 0:
            cp.pop(lowest_p_gen[0])
            p_vals.pop(lowest_p_gen[0])
            gen_id.pop(lowest_p_gen[0])
            lowest_p_gen.pop(0)
        maxgen = maxgen - 1

def dividereturn(seq,cp,p_vals,offset,maximum,gen_id,gen=0):
    loc,p_v = TTest(seq)
    if (len(cp)-1 <= (maximum-gen) and (loc!=0)):# or (p_v < 0.005):
        cp.append(loc+offset)
        p_vals.append(p_v)
        gen_id.append(gen)
        dividereturn(seq[:loc],cp,p_vals,offset,maximum,gen_id,gen+1)
        dividereturn(seq[loc:],cp,p_vals,offset+loc,maximum,gen_id,gen+1)




def mc_p(seq,p=0.05):
    cp = []
    divide_return(seq,p,cp,0)
    return sorted(cp)

def divide_return(seq,p,cp,offset):
    loc = T_Test(seq,p)
    if loc !=  0:
        cp.append(loc+offset)
        if loc > 7:
            divide_return(seq[:loc],p,cp,offset)
        if len(seq)-loc > 7:
            divide_return(seq[loc:],p,cp,offset+loc)

    
#Multiple change detection using T-Test with or withou6t fixed p
def TTest(seq):
    offset = 2
    if len(seq) < offset*2+1:
        return 0,1
    else:    
        stop = 0
        A = []
        i = offset
        length = len(seq)
        while i < length-offset:
            A.append(stats.ttest_ind(seq[0:i], seq[i:length], axis=0, equal_var=True)[1])
            i = i+1
        A = pd.Series(A).abs()
        stop = A.idxmin()+offset
            
        return stop,A.min()


def T_Test(seq,p=0.01): #min length 7
    stop = 0
    A = []
    offset = 3
    i = offset
    length = len(seq)
    while i < length-offset:
        A.append(stats.ttest_ind(seq[0:i], seq[i:length], axis=0, equal_var=True)[1])
        i = i+1
    A = pd.Series(A).abs()
    if A.min() < p:#/length: Bonferroni correction for multiple hypothesis testing most probable cp
        stop = A.idxmin()+offset
        
    return stop