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
from scipy import stats
from scipy import special

def TTest_ANOVA(seq,no,offset):
    """
    INPUT: 
    seq: list of doubles, time series with KNOWN number of change-points
    no: integer, number of change-points
    offset: minimum distance between changes. If offset < 2, an error arises because variance can be computed
    OUTPUT:
    best[0] is an indexing array for seq for the most likely locations for change in mean, according to ANOVA
    """
    if (no+1)*3 > len(seq):
        return [0.]
    else:
        length = len(seq)
        loc = [(x+1)*offset for x in range(no)]
        best = [loc,_test_loc(seq,loc,offset)]
        while next(loc,length,no,offset):
            if loc[-1] > length:
                break
            p = _test_loc(seq,loc,offset)
            if best[1] > p:
                best = [copy.copy(loc),p]
        return best[0]

def Optimal_ANOVA(seq,thres,offset):
    length = len(seq)
    cp = []
    p = 1.
    ptemp = 0.
    no = 1

    while no < 5:
        cptemp = TTest_ANOVA(seq,no,offset)
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

    return p_ttest

F_onewayResult = namedtuple('F_onewayResult', ('statistic', 'pvalue'))

def _sum_of_squares(a, axis=0):
    """
    Squares each element of the input array, and returns the sum(s) of that.
    Parameters
    ----------
    a : array_like
        Input array.
    axis : int or None, optional
        Axis along which to calculate. Default is 0. If None, compute over
        the whole array `a`.
    Returns
    -------
    sum_of_squares : ndarray
        The sum along the given axis for (a**2).
    See also
    --------
    _square_of_sums : The square(s) of the sum(s) (the opposite of
    `_sum_of_squares`).
    """
    a, axis = _chk_asarray(a, axis)
    return np.sum(a*a, axis)


def _chk_asarray(a, axis):
    if axis is None:
        a = np.ravel(a)
        outaxis = 0
    else:
        a = np.asarray(a)
        outaxis = axis

    if a.ndim == 0:
        a = np.atleast_1d(a)

    return a, outaxis

def _square_of_sums(a, axis=0):
    """
    Sums elements of the input array, and returns the square(s) of that sum.
    Parameters
    ----------
    a : array_like
        Input array.
    axis : int or None, optional
        Axis along which to calculate. Default is 0. If None, compute over
        the whole array `a`.
    Returns
    -------
    square_of_sums : float or ndarray
        The square of the sum over `axis`.
    See also
    --------
    _sum_of_squares : The sum of squares (the opposite of `square_of_sums`).
    """
    a, axis = _chk_asarray(a, axis)
    s = np.sum(a, axis)
    if not np.isscalar(s):
        return s.astype(float) * s
    else:
        return float(s) * s

#From scipy stats with the minor difference that the argument of ANOVA contains data insted of *args
def oneway_ANOVA(data):
    args = [np.asarray(arg, dtype=float) for arg in data]

    # ANOVA on N groups, each in its own array
    num_groups = len(args)
    alldata = np.concatenate(args)
    bign = len(alldata)

    # Determine the mean of the data, and subtract that from all inputs to a
    # variance (via sum_of_sq / sq_of_sum) calculation.  Variance is invariance
    # to a shift in location, and centering all data around zero vastly
    # improves numerical stability.
    offset = alldata.mean()
    alldata -= offset

    sstot = _sum_of_squares(alldata) - (_square_of_sums(alldata) / float(bign))
    ssbn = 0
    for a in args:
        ssbn += _square_of_sums(a - offset) / float(len(a))

    # Naming: variables ending in bn/b are for "between treatments", wn/w are
    # for "within treatments"
    ssbn -= (_square_of_sums(alldata) / float(bign))
    sswn = sstot - ssbn
    dfbn = num_groups - 1
    dfwn = bign - num_groups
    msb = ssbn / float(dfbn)
    msw = sswn / float(dfwn)
    f = msb / msw

    prob = special.fdtrc(dfbn, dfwn, f)   # equivalent to stats.f.sf

    return F_onewayResult(f, prob)


def _test_loc(seq,loc,offset):
    p_ttest = []
    p_out = []

    data =[list(seq[:loc[0]])]
    
    i = 0
    while i < len(loc): #intervals defined by loc. Take them and compare them seq[k],seq[l],seq[m]
        if i == 0:
            k = 0
        else:
            k = loc[i-1]
        l = loc[i]
        if i == len(loc)-1:
            m = len(seq) - 1
        else:
            m = loc[i+1]
        
        if l-k < offset or m-l < offset:
            return 5.
        data.append(list(seq[l:m]))
        #p_ttest.append(stats.ttest_ind(seq[k:l], seq[l:m], axis=0, equal_var=True)[1])
        
        i = i + 1
    return oneway_ANOVA(data)[1]  #this is giving an error, because the input of f_oneway is suppoed to be of the form (data1,data2,...)


