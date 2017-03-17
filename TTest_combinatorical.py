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
    if (no+1)*offset > len(seq):
        return [0.]
    else:
        length = len(seq)
        loc = [(x+1)*offset for x in range(no)]
        best = [loc,_test_loc(seq,loc,offset)]
        while next(loc,length,no):
            if loc[-1] > length:
                break
            p = _test_loc(seq,loc,offset)
            if best[1] > p:
                best = [copy.copy(loc),p]
        return best[0]

def Full_combinatorical(seq,offset=2):
    no = 1
    p_prev = 1.
    loc_prev = []
    loc = TTest_combinatorical(seq,no,offset)
    p = _test_loc(seq,loc,offset)
    while p < p_prev:
        print loc,p
        p_prev = p
        no = no + 1
        loc_prev = loc
        loc = TTest_combinatorical(seq,no,offset)
        p = _test_loc(seq,loc,offset)
    return loc_prev

def Optimal_TTest(seq,thres,offset):
    length = len(seq)
    cp = []
    p = 1.
    ptemp = 0.
    no = 1

    while no < 5:
        cptemp = TTest_combinatorical(seq,no,offset)
        if any(a > thres for a in _test_loc_p(seq,cptemp,offset)):
            return cp
        else:
            ptemp = _test_loc(seq,cptemp,offset)
            if ptemp < p:
                cp = cptemp
                p = ptemp
            else:
                return cp
        no = no + 1

def next(loc,length,no):
    i = -1
    while loc[i] >= length+3*i:
        if abs(i) >= no:
            break
        else:
            i = i -1
    loc[i]=loc[i]+1
    for j in range(abs(i)-1):
        loc[i+j+1] = loc[i+j]+3
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
        
    return p_ttest

def _test_loc(seq,loc,offset):
    p_ttest = []
    p_out = []

    if len(loc) == 0:
        return 1.
    
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


def ttest_ind(a, b, axis=0, equal_var=True):
    """
    Calculates the T-test for the means of TWO INDEPENDENT samples of scores.
    This is a two-sided test for the null hypothesis that 2 independent samples
    have identical average (expected) values. This test assumes that the
    populations have identical variances.
    Parameters
    ----------
    a, b : array_like
        The arrays must have the same shape, except in the dimension
        corresponding to `axis` (the first, by default).
    axis : int, optional
        Axis can equal None (ravel array first), or an integer (the axis
        over which to operate on a and b).
    equal_var : bool, optional
        If True (default), perform a standard independent 2 sample test
        that assumes equal population variances [1]_.
        If False, perform Welch's t-test, which does not assume equal
        population variance [2]_.
        .. versionadded:: 0.11.0
    Returns
    -------
    t : float or array
        The calculated t-statistic.
    prob : float or array
        The two-tailed p-value.
    """
    
    a, b, axis = _chk2_asarray(a, b, axis)
    if a.size == 0 or b.size == 0:
        return (np.nan, np.nan)

    v1 = np.var(a, axis, ddof=1)
    v2 = np.var(b, axis, ddof=1)
    n1 = a.shape[axis]
    n2 = b.shape[axis]

    if (equal_var):
        df = n1 + n2 - 2
        svar = ((n1 - 1) * v1 + (n2 - 1) * v2) / float(df)
        denom = np.sqrt(svar * (1.0 / n1 + 1.0 / n2))
    else:
        vn1 = v1 / n1
        vn2 = v2 / n2
        df = ((vn1 + vn2)**2) / ((vn1**2) / (n1 - 1) + (vn2**2) / (n2 - 1))

        # If df is undefined, variances are zero (assumes n1 > 0 & n2 > 0).
        # Hence it doesn't matter what df is as long as it's not NaN.
        df = np.where(np.isnan(df), 1, df)
        denom = np.sqrt(vn1 + vn2)

    d = np.mean(a, axis) - np.mean(b, axis)
    t = np.divide(d, denom)
    t, prob = _ttest_finish(df, t)

    return t, prob

if __name__ == "__main__":
    np.random.seed(2)
    length = 20
    noise = 3
    step = 1
    test_sequence = np.concatenate((np.random.randn(length)*noise, np.random.randn(length)*noise+step,np.random.randn(length)*noise))
    loc = Full_combinatorical(test_sequence)
    print loc
