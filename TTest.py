import numpy as np 
from scipy import stats


#Easiest parametric test
#Otehrwise http://www.jmlr.org/papers/volume13/gretton12a/gretton12a.pdf

""" [m,n] = [0,n] - [0,m]
# eq 15 of http://www.jstor.org/stable/2683386
def two_samples(Tm,Tn,Sm,Sn,m,n):
	if n > m:
		Tx = Tn - Tm
		Sx = Sn - Sm - 1.*m/(n*(n-m)) * ((1.*(n-m)/m)*Tm - Tx)**2
	else:
		Tx = Tm - Tn
		Sx = Sm - Sn - 1. * n / (m*(m-n)) * ((1.*(m-n)/n)*Tn - Tx)**2
	return Tx,Sx,np.abs(n-m)

def convert(T,S,n):
	return 1.*T/n, 1.*S/(n-1)

"""


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

#Equation 1.5 from Chan, Golub, LeVeque
def updating(m,m1,v,v1,n,n1):
    T = m*n
    T1 = m1*n1
    S = v*(n-1)
    S1 = v1*(n1-1)
    T2 = T - T1
    S2 = S - S1 - 1.* ((1.*T1*n/n1 - T)**2)  *n1/(n*(n-n1)) 
    m2 = T2/(n-n1)
    v2 = S2/(n-n1)
    return m2,v2

#Fill in lookup table of mean and variances by computing subintervals starting at i in an online fashion.
#This is the fastest way of doing it and it is O(n**2)
#def lookup_table(seq):
#    Table= {}
#    i = 0
#    while i < len(seq):
#        m,v = online_variance(seq[i:])
#        j = 0
#        while j < len(m):
#            Table[str((i,i+j+1))] = [m[j],v[j]]
#            j = j + 1
#        i = i + 1
#    return Table
#Lookup table of intervals from [,j] for j in [0,n]
def lookup_table(seq):
    Table = {}
    i = 0
    m,v = online_variance(seq)
    j = 0
    while j < len(seq):
        Table[str(j)] = [m[j],v[j]]
        j = j + 1
    return Table

def get_m_v_l(seq,loc,Table='none'):
    m = []
    v = []
    lengths=[]
    if Table != 'none':
        j = 0
        while j < len(loc):
            if j > 1:
                m_full, v_full = Table[str(loc[j])]
                m_first, v_first = Table[str(loc[j-1])]
                m_second, v_second = updating(m_full, m_first, v_full,v_first,loc[j],loc[j-1] )
                m.append(m_second)
                v.append(v_second)
                lengths.append(loc[j]-loc[j-1])
            else:
                m_full, v_full = Table[str(loc[j])]
                m.append(m_full)
                v.append(v_full)
                lengths.append(loc[j])
            j=j+1
    else:
        m.append(np.mean(seq[:loc[0]]))
        v.append(np.var(seq[:loc[0]]))
        lengths.append(loc[0])
        j = 1
        while j < len(loc):
            m.append(np.mean(seq[loc[j-1]:loc[j]]))
            v.append(np.var(seq[loc[j-1]:loc[j]]))
            lengths.append(loc[j]-loc[j-1])  
            j = j + 1
    return m,v,lengths



def T_test_loc(seq,loc,lookup='none',offset=2,p='no'):
    p_ttest = []

    m,s,lengths = get_m_v_l(seq,loc,lookup)

        
    for m1,s1,l1,m2,s2,l2 in zip(m[:-1],s[:-1],lengths[:-1],m[1:],s[1:],lengths[1:]):
        p_ttest.append(ttest_ind_from_stats(m1, s1, l1, m2, s2, l2))

    if p=='no':
        return stats.chisqprob(-2.*np.sum([np.log(x) for x in p_ttest]),2*len(loc))
    else:
        return p_ttest



def ttest_ind_from_stats(mean1, s1, nobs1, mean2, s2, nobs2):
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
    df, denom = _unequal_var_ttest_denom(s1, nobs1, s2, nobs2)

    res = _ttest_ind_from_stats(mean1, mean2, denom, df)
    return res[1]


def _unequal_var_ttest_denom(v1, n1, v2, n2):
    vn1 = v1 / n1
    vn2 = v2 / n2
    with np.errstate(divide='ignore', invalid='ignore'):
        df = (vn1 + vn2)**2 / (vn1**2 / (n1 - 1) + vn2**2 / (n2 - 1))

    # If df is undefined, variances are zero (assumes n1 > 0 & n2 > 0).
    # Hence it doesn't matter what df is as long as it's not NaN.
    df = np.where(np.isnan(df), 1, df)
    denom = np.sqrt(vn1 + vn2)
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
    if np.any(t_isnan):# and __scipy_prior0101:
        # older scipy's would return 0 for nan values of the argument
        # which is incorrect
        if np.isscalar(prob):
            prob = np.nan
        else:
            prob[t_isnan] = np.nan

    if t.ndim == 0:
        t = np.asscalar(t)

    return t, prob
