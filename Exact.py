# -*- coding: utf-8 -*-
"""
Created on Mon Apr 18 11:59:19 2016

@author: Sir Thomas
"""


import numpy as np
import scipy as scipy
import copy
import itertools


import KSTest as KS
import TTest as T 
import MWTest as MW
import ADTest as AD 
import LLTest as LL
import P_value as P

def flatten(l):
  out = []
  for item in l:
    if isinstance(item, (list, tuple)):
      out.extend(flatten(item))
    else:
      out.append(item)
  return out

def Exact(seq,test_loc='TT',no=-1,offset=2,thres=0.05, constraints=[]):
    #Initialization:
    if test_loc == 'KS':
        test_loc = KS.KS_test_loc
    elif test_loc == 'AD':
        test_loc = AD.AD_test_loc
    elif test_loc == 'MW':
        test_loc = MW.MW_test_loc
    elif test_loc == 'LL':
        test_loc = LL.LL_test_loc
    else:
        test_loc = T.T_test_loc


    if constraints == []:
        constraints = range(len(seq))
    else:
        if type(constraints[0]) == int:
            constraints = [constraints]
        constraints = [range(i[0],i[1]) for i in constraints]
        constraints = flatten(constraints)


    if test_loc == T.T_test_loc:
        if no >= 2:
            Table = T.lookup_table(seq)  #create lookup table (takes a while)
        else:
            Table= 'none'
    else:
        Table = 'none'


    if no == -1:
        #Look for the optimal number of changes
        #Look for one change:
        loc1 = Combinatorical(seq,1,offset,test_loc,Table,constraints)
        p1 = test_loc(seq,loc1)
        #print "Best first change" , p1, loc1


        #Look for two changes
        loc2 = Combinatorical(seq,2,offset,test_loc,Table,constraints)
        p2 = test_loc(seq,loc1)
        #print "Best second change", p2, loc2

        #if there is only one change
        if p2 > p1:
            if p1 < thres:
                return loc1,P.Permutation(p1,loc1,seq,test_loc)
            else:
                return [],1.

        #Check for three changes
        middle,loc3 = best_middle(seq,offset,test_loc,Table,constraints)
        p3 = test_loc(seq,loc3)
        #print "Best Third change",p3, loc3

        #if there are only two changes
        if p3 > p2:
            if p2 < thres:
                return loc2,P.Permutation(p2,loc1,seq,test_loc)
            else:
                return [],1.
        #Look for more than three changes
        indices = prune_min(seq,middle,offset)
        locn,pn = Find_Best_Combination(seq,indices,no,offset,test_loc, Table)

        #Is is three changes the optimum?
        if pn > p3:
            if p3 < thres:
                return loc3,P.Permutation(p3,loc1,seq,test_loc)
            else:
                return [],1.
        else:
            if pn < thres:
                return locn,P.Permutation(pn,loc1,seq,test_loc)
            else:
                return [],1.
    else:
        if no > 3:
            middle,loc3 = best_middle(seq,offset,test_loc,Table,constraints)
            indices = prune_min(seq,middle,offset)
            locn,pn = Find_Best_Combination(seq,indices,no,offset,test_loc, Table)
            return locn, P.Permutation(pn,locn,seq,test_loc)
        else:
            loc1 = Combinatorical(seq,no,offset,test_loc,Table,constraints)
            p1 = P.Permutation(test_loc(seq,loc1),loc1,seq,test_loc)
            return loc1,p1



#From the chandidate indices, pick a combination that gives the best score
def Find_Best_Combination(seq,indices,no=-1,offset=2,test_loc=T.T_test_loc,Table='none'):
    p = 100.
    best_loc = []
    if no == -1:
        for no_cp in range(4,len(indices)+1):
            #Iterate over number of changepoints and check every possible configuration of local minima
            for comb in itertools.combinations(indices,no_cp):
                loc = list(comb)
                spacing = [loc[j]-loc[i] >= offset for i,j in zip(range(len(loc)-2),range(1,len(loc)-1))]
                if all(spacing) and len(spacing)>0:
                    if test_loc(seq,loc,Table,offset)< p:
                        p = test_loc(seq,loc,Table,offset)
                        best_loc = loc
            if len(best_loc) == no_cp-1:
                return best_loc, test_loc(seq,best_loc)
    else:
        for comb in itertools.combinations(indices,no):
            loc = list(comb)
            spacing = [loc[j]-loc[i] >= offset for i,j in zip(range(len(loc)-1),range(1,len(loc)))]
            if all(spacing) and len(spacing)>0:
                if test_loc(seq,loc,Table,offset)< p:
                    p = test_loc(seq,loc,Table,offset)
                    best_loc = loc
    return best_loc,test_loc(seq,best_loc)

#Look for local minimum in dictionary on checked changepoints
def prune_min(seq,middle,offset):
    p_dir = 1
    p_local = 1.
    save_i0 = 0
    save_i1 = 0
    save_i2 = 0
    indices = []
    min_start = len(seq)-offset
    max_stop = offset

    for i in range(len(seq)):
        if str(i) in middle.keys():
            if middle[str(i)][2] < p_local:
                p_dir = -1
                p_local = middle[str(i)][2]
                save_i1 = i
                save_i0 = middle[str(i)][0]
                save_i2 = middle[str(i)][1]
            elif middle[str(i)][2] > p_local:
                if p_dir == -1:
                    indices.append(save_i0)
                    indices.append(save_i1)
                    indices.append(save_i2)
                    #print save_i, middle[str(save_i)]
                p_dir = 1
                p_local = middle[str(i)][2]
    if min_start not in indices:
        indices.append(min_start)
    if max_stop not in indices:
        indices.append(max_stop)
    return list(np.unique(indices))


#This checks all configurations for a given no
#Theoretically it suffices to scan for no=3 and imput good combinations from the results of all tested subintervals
def Combinatorical(seq,no,offset,test_loc=T.T_test_loc,Table='none',constraints=[]):
    if constraints == []:
        constraints = seq
    if (no+1)*offset > len(seq):
        raise NameError('Minimum spacing not possible')
    else:
        length = len(seq)
        loc = [(x+1)*offset for x in range(no)] #initial spacing (all to the right)
        best = [loc,test_loc(seq,loc,Table,offset)] #keep track of best configuration so far
        while next(loc,length,no,offset): 
            if loc[-1] > length:
                break
            p = test_loc(seq,loc,Table,offset)
            if best[1] > p:
                best = [copy.copy(loc),p]
        return best[0]

def best_middle(seq,offset=2,test_loc=T.T_test_loc,Table='none',constraints=[]):
    if constraints == []:
        constraints = seq
    length = len(seq)
    middle = {}
    no = 3
    loc = [(x+1)*offset for x in range(no)] #initial spacing (all to the right)
    best = [loc,test_loc(seq,loc,Table,offset)] #keep track of best configuration so far
    while next(loc,length,no,offset): 
        if loc[-1] > length:
            break
        if all( [inst in constraints for inst in loc]): 
            p = test_loc(seq,loc,Table,offset)
            if best[1] > p:
                best = [copy.copy(loc),p]
            if str(loc[1]) in middle.keys():
                if middle[str(loc[1])][2] > p:
                    middle[str(loc[1])] = [loc[0],loc[2],p]
            else:
                middle[str(loc[1])] = [loc[0],loc[2],p]
    return middle, best[0]



#move changes along the sequence
def next(loc,length,no,offset):
    i = -1
    #Find the index to change and not to violate mimimum spacing condition
    while loc[i] >= length+offset*i:
        #If all number of changes have been exhausted, break the loop in TTest_combinatorical
        if abs(i) >= no:
            return 0
            break
        else:
            i = i -1
    loc[i]=loc[i]+1
    for j in range(abs(i)-1):
        loc[i+j+1] = loc[i+j]+offset
    return 1
