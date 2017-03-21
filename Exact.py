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
import scipy as scipy
import copy
import itertools


#Cose your method to correct for multiple hypothesis testing
def Correction(p,no,seq,offset):
    p_corr = p
    #p_corr = p*np.math.factorial(no)*len(seq) #theoretical lower bound on number of tests
    #p_corr = p*scipy.misc.comb(len(seq)-offset*no,no) #actual tests performed
    return p_corr

def Exact(seq,no=-1,offset=2,thres=0.05):
    #Initialization:
    Table = lookup_table(seq)  #create lookup table (takes a while)


    #Look for one change:
    loc1 = Combinatorical(seq,1,offset,Table)
    p1 = Correction(test_loc(seq,loc1,Table, offset),1,seq,offset)
    #print "Best first change" , p1, loc1
    if no == 1:
        return loc1,p1


    #Look for two changes
    loc2 = Combinatorical(seq,2,offset,Table)
    p2 = Correction(test_loc(seq,loc2,Table, offset),2,seq,offset)
    #print "Best second change", p2, loc2
    if no == 2:
        return loc2,p2

    #if there is only one change
    if p2 > p1:
        if p1 < thres:
            return loc1,p1
        else:
            return [],1.

    #Check for three changes
    middle,loc3 = best_middle(seq,Table,offset)
    p3 = Correction(test_loc(seq,loc3,Table, offset),3,seq,offset)
    #print "Best Third change",p3, loc3
    if no == 3:
        return loc3,p3

    #if there are only two changes
    if p3 > p2:
        if p2 < thres:
            return loc2,p2
        else:
            return [],1.

    #Look for more than three changes
    indices = prune_min(seq,middle)
    locn,pn = Find_Best_Combination(seq,indices, Table,no)
    if no == len(locn):
        return locn,pn

    #Is is three changes the optimum?
    if pn > p3:
        if p3 < thres:
            return loc3,p3
        else:
            return [],1.
    else:
        if pn < thres:
            return locn,pn
        else:
            return [],1.

#From the chandidate indices, pick a combination that gives the best score
def Find_Best_Combination(seq,indices, Table,no=-1,offset=2):
    p = 1.
    best_loc = []
    if no != -1:
        for no_cp in range(4,len(indices)+1):
            #Iterate over number of changepoints and check every possible configuration of local minima
            for comb in itertools.combinations(indices,no_cp):
                loc = list(comb)
                spacing = [loc[j]-loc[i] >= offset for i,j in zip(range(len(loc)-2),range(1,len(loc)-1))]
                if all(spacing) and len(spacing)>0:
                    if Correction(test_loc(seq,loc,Table,offset),len(loc),seq,offset )< p:
                        p = Correction(test_loc(seq,loc,Table,offset),len(loc),seq,offset )
                        best_loc = loc
            #print no_cp-1 , best_loc, p
            if len(best_loc) == no_cp-1:
                return best_loc, p
    else:
        loc = list( itertools.combinations(indices,no))
        spacing = [loc[j]-loc[i] >= offset for i,j in zip(range(len(loc)-2),range(1,len(loc)-1))]
        if all(spacing) and len(spacing)>0:
            if Correction(test_loc(seq,loc,Table,offset),len(loc),seq,offset )< p:
                p = Correction(test_loc(seq,loc,Table,offset),len(loc),seq,offset )
                best_loc = loc
    return best_loc,p

#Look for local minimum in dictionary on checked changepoints
def prune_min(seq,middle):
    p_dir = 1
    p_local = 1.
    save_i = 0
    indices = []

    for i in range(len(seq)):
        for j in range(len(seq)):
            if str(i) in middle.keys() and str(j) in middle.keys():
                if j in middle[str(i)]:
                    if middle[str(i)][2] < p_local:
                        p_dir = -1
                        p_local = middle[str(i)][2]
                        save_i = i
                    elif middle[str(i)][2] > p_local:
                        if p_dir == -1:
                            indices.append(save_i)
                            #print save_i, middle[str(save_i)]
                        p_dir = 1
                        p_local = middle[str(i)][2]
    return indices



def Full_combinatorical(seq,offset=2,thres=0.05):
    #Initialization:
    no = 1 #start at 1 changepoint
    p_prev = 1. #significance of previous no
    loc_prev = [] #there are no previously found changepoints
    Table = lookup_table(seq)  #create lookup table (takes a while)

    #print "Lets go"

    #Look for one change:
    loc = Combinatorical(seq,no,offset,Table)
    p = Correction(test_loc(seq,loc,Table, offset),no,seq,offset)
    #print p, loc

    #Iterate over number of changes until an additional change does not benefit the score
    while p < p_prev:
        p_prev = p
        no = no + 1
        loc_prev = loc
        loc = Combinatorical(seq,no,offset,Table)
        p = Correction(test_loc(seq,loc,Table, offset),no,seq,offset)
        #print p, loc

    if p_prev < thres:
        return loc_prev #report findings only if threshold condition is satisfied
    else:
        return []


#This checks all configurations for a given no
#Theoretically it suffices to scan for no=3 and imput good combinations from the results of all tested subintervals
def Combinatorical(seq,no,offset,Table):
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

def best_middle(seq,Table,offset=2):
    length = len(seq)
    middle = {}
    no = 3
    loc = [(x+1)*offset for x in range(no)] #initial spacing (all to the right)
    best = [loc,test_loc(seq,loc,Table,offset)] #keep track of best configuration so far
    while next(loc,length,no,offset): 
        if loc[-1] > length:
            break
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

def online_variance(data):
    n = 0
    mean = 0.0
    M2 = 0.0
    v = []
    m = []
     
    for x in data:
        n += 1
        delta = x - mean
        mean += delta/n
        M2 += delta*(x - mean)
        if n <= 1:
            v.append('nan')
            m.append(mean)
        if n > 1:
            v.append(M2/(n-1))
            m.append(mean)
    return m,v

#Fill in lookup table of mean and variances by computing subintervals starting at i in an online fashion.
def lookup_table(seq):
    Table= {}
    i = 0
    while i < len(seq):
        m,v = online_variance(seq[i:])
        j = 0
        while j < len(m):
            Table[str((i,i+j+1))] = [m[j],v[j]]
            j = j + 1
        i = i + 1
    return Table


#Compute optimization objective for loc
def test_loc(seq,loc,lookup,offset):
    p_ttest = []
    p_out = []
    
    i = 0
    while i < len(loc): #intervals defined by loc. Take them and compare them seq[index1],seq[index2],seq[index3]
        if i == 0:
            index1 = 0
        else:
            index1 = loc[i-1]
        index2 = loc[i]
        if i == len(loc)-1:
            index3 = len(seq)
        else:
            index3 = loc[i+1]

        if index2-index1 < offset or index3-index2 < offset:
            raise NameError('The minimum spacing conditions are not fulfilled for the testd loc '+str(loc))

        m1, v1 = lookup[str((index1,index2))]
        m2, v2 =  lookup[str((index2,index3))]
        
        p_ttest.append(ttest_ind_from_stats(m1, math.sqrt(v1), index2-index1, m2, math.sqrt(v2), index3-index2)[1])
        
        i = i + 1
        
    return stats.chisqprob(-2.*np.sum([np.log(x) for x in p_ttest]),2*len(loc)) 






Ttest_indResult = namedtuple('Ttest_indResult', ('statistic', 'pvalue'))
def ttest_ind_from_stats(mean1, std1, nobs1, mean2, std2, nobs2):
    """
    T-test for means of two independent samples from descriptive statistics.
    This is a two-sided test for the null hypothesis that 2 independent samples
    have identical average (expected) values.
    Parameters
    ----------
    mean1 : array_like
        The mean(s) of sample 1.
    std1 : array_like
        The standard deviation(s) of sample 1.
    nobs1 : array_like
        The number(s) of observations of sample 1.
    mean2 : array_like
        The mean(s) of sample 2
    std2 : array_like
        The standard deviations(s) of sample 2.
    nobs2 : array_like
        The number(s) of observations of sample 2.
    equal_var : bool, optional
        If True (default), perform a standard independent 2 sample test
        that assumes equal population variances [1]_.
        If False, perform Welch's t-test, which does not assume equal
        population variance [2]_.
    Returns
    -------
    statistic : float or array
        The calculated t-statistics
    pvalue : float or array
        The two-tailed p-value.
    See also
    --------
    scipy.stats.ttest_ind
    Notes
    -----
    .. versionadded:: 0.16.0
    References
    ----------
    .. [1] http://en.wikipedia.org/wiki/T-test#Independent_two-sample_t-test
    .. [2] http://en.wikipedia.org/wiki/Welch%27s_t_test
    """
    df, denom = _equal_var_ttest_denom(std1**2, nobs1, std2**2, nobs2)

    res = _ttest_ind_from_stats(mean1, mean2, denom, df)
    return Ttest_indResult(*res) 


def _equal_var_ttest_denom(v1, n1, v2, n2):
    df = n1 + n2 - 2.0
    svar = ((n1 - 1) * v1 + (n2 - 1) * v2) / df
    denom = math.sqrt(svar * (1.0 / n1 + 1.0 / n2))
    return df, denom

def _ttest_ind_from_stats(mean1, mean2, denom, df):

    d = mean1 - mean2
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.divide(d, denom)
    t, prob = _ttest_finish(df, t,'two-sided')

    return (t, prob)


def _ttest_finish(df, t,alternative):
    """Common code between all 3 t-test functions."""
    dist_gen = stats.distributions.t
    if alternative == 'two-sided':
        prob = dist_gen.sf(np.abs(t), df) * 2 # use np.abs to get upper alternative
    elif alternative == 'greater':
        prob = dist_gen.sf(t, df)
    elif alternative == 'less':
        prob = dist_gen.cdf(t, df)
    else:
        raise ValueError("Unknown alternative %r" % alternative)

    t_isnan = np.isnan(t)
    if np.any(t_isnan) and __scipy_prior0101:
        # older scipy's would return 0 for nan values of the argument
        # which is incorrect
        if np.isscalar(prob):
            prob = np.nan
        else:
            prob[t_isnan] = np.nan

    if t.ndim == 0:
        t = np.asscalar(t)

    return t, prob




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