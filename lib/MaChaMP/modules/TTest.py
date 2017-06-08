import numpy as np
from scipy import stats
from os.path import join

import logging
loggerT = logging.getLogger(__name__)
formatterT = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
file_handlerT = logging.FileHandler(join('meta', 'TTest_log.log'))
file_handlerT.setFormatter(formatterT)
loggerT.addHandler(file_handlerT)
loggerT.setLevel(logging.DEBUG)


def online_variance(data):
    """
    Equation 1.5 from Chan, Golub, LeVeque

    INPUT:
    data: time series

    OUTPUT:
    m: list of means for all subintervals in data starting at 0
    v: list of variances for all subintervals in data starting at 0
    """

    # Initialization
    n = 0.
    mean = 0.0
    M2 = 0.0
    v = []
    m = []

    for x in data:
        n += 1
        delta = x - mean
        mean += delta / n
        M2 += delta * (x - mean)
        if n <= 1:
            v.append('nan')
            m.append(mean)
        if n > 1:
            v.append(M2 / (n))  # biased
            m.append(mean)
    return m, v


def updating(m, m1, v, v1, n, n1):
    """
    Equation 2.1a and 2.1b from Chan, Golub, LeVeque
    [m1,v1] + [m2,v2] = [m,v]

    INPUT:
    m: mean of full interval
    m1: mean of one part of the full interval
    v: variance of the full interval
    v1: variance of one part of the full interval
    n: size of the full interval
    n1: size of one part of the inveral

    OUTPUT:
    Mean and variance of the other part of the interval
    """
    T = m * n
    T1 = m1 * n1
    S = v * (n)             # biased
    S1 = v1 * (n1)      # biased
    T2 = T - T1
    S2 = S - S1 - 1. * ((1. * T1 * n / n1 - T)**2) * n1 / (n * (n - n1))
    m2 = T2 / (n - n1)  # biased
    v2 = S2 / (n - n1)  # biased

    return m2, v2


def lookup_table(seq):
    """
    A dictionary containing mean and variances from all subintervals of seq satrting at 0

    INPUT:
    seq: Time series

    OUTPUT:
    Table: dictionary of [mean,variance]
    """
    Table = {}
    m, v = online_variance(seq)
    j = 1
    while j <= len(seq):
        Table[str(j)] = [m[j - 1], v[j - 1]]
        j = j + 1
    return Table


def get_m_v_l(seq, loc, Table='none'):
    """
    For a division into subintervals, get all mean, varaicnes and lengths of those

    INPUT:
    seq: time series
    loc: changepoints as division into subintervals
    Table: lookup table of all subintervals starting at 0

    OUTPUT:
    m,v,l: lists of means variances and lengths

    >>> [i==j for i,j in  zip(get_m_v_l([1,1,2,1,5,5,6,5],[4], lookup_table([1,1,2,1,5,5,6,5]))[1],  get_m_v_l([1,1,2,1,5,5,6,5],[4])[1])]
    [True, True]
    """

    # Inizialize return values
    m = []                          # List of means
    v = []                           # List of variances
    lengths = []                 # List of intervals lengths
    loc = loc + [len(seq)]  # The last interval ends at the end of the sequence

    # If a lookup Table has  been precomputed
    if Table != 'none':

        # For the first interval: Just look up the values
        j = 0
        m_full, v_full = Table[str(loc[j])]
        m.append(m_full)
        v.append(v_full)
        lengths.append(loc[j])

        # Iterate over the remaining len(loc) intervals
        j = 1
        while j < len(loc):
            m_full, v_full = Table[str(loc[j])]
            m_first, v_first = Table[str(loc[j - 1])]
            m_second, v_second = updating(m_full, m_first, v_full, v_first, loc[j], loc[j - 1])
            m.append(m_second)
            v.append(v_second)
            lengths.append(loc[j] - loc[j - 1])
            j = j + 1
    # if there has not been a lookup table precomputed, do means and variances now by hand
    else:
        loggerT.warning('Runtime warning: No Table precomputed for ' + str(loc))
        m.append(np.mean(seq[:loc[0]]))
        v.append(np.var(seq[:loc[0]]))
        lengths.append(loc[0])
        j = 1
        while j < len(loc):
            m.append(np.mean(seq[loc[j - 1]:loc[j]]))
            v.append(np.var(seq[loc[j - 1]:loc[j]]))
            lengths.append(loc[j] - loc[j - 1])
            j = j + 1

    loggerT.debug('For changes ' + str(loc) + ', mean=' + str(m) + ', variances=' + str(v) + ', lengths=' + str(lengths))
    return m, v, lengths  # gives biased (numpy) mean and variance


def T_test_loc(seq, loc, lookup='none', offset=2, p='no', Reg=False):
    """
    Evaluation of Welch-Fisher optimization objective without Bonferroni-correction

    INPUT:
    seq: time series
    loc: Division into intervals at incices stored in loc
    lookup: table computed by lookup_table (see above)
    offset: minimal interval size
    p: if true, not not apply Fishers method and return list of single p-values
    Reg: dummy variable required for implementing other optimization objectives

    OUTPUT:
    evaluated Welch-Fisher objective
    """

    # Initialize list of p-values for single t-Test and obtain means and variances
    p_ttest = []
    m, s, lengths = get_m_v_l(seq, loc, lookup)

    # For every consecutive pair of intervals, evaluate Welchs t-Test
    for m1, s1, l1, m2, s2, l2 in zip(m[:-1], s[:-1], lengths[:-1], m[1:], s[1:], lengths[1:]):
        p_ttest.append(ttest_ind_from_stats(m1, s1, l1, m2, s2, l2))

    if p == 'no':
        # If demanded to evaluate completely and there are more than one tests, apply Fishers methods
        if len(p_ttest) > 1:
            return stats.chisqprob(-2. * np.sum([np.log(x) for x in p_ttest]), 2 * (len(loc)))  # without Bonferroni correction
        else:
            return p_ttest[0]
    else:
        return p_ttest


# The following code is copied from scipy. sp.stats version of ttest does not allow the input of tests statistics
# It is therefore slow, because it does not make use of our speedup trick.

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
    df, denom = _unequal_var_ttest_denom(s1, nobs1, s2, nobs2)  # Denominator and d.o.f. for ttest

    res = _ttest_ind_from_stats(mean1, mean2, denom, df)
    return res[1]


def _unequal_var_ttest_denom(v1, n1, v2, n2):
    vn1 = v1 / n1
    vn2 = v2 / n2
    with np.errstate(all='ignore'):
        den = (vn1**2 / (n1 - 1) + vn2**2 / (n2 - 1))
        if den != 0:
            df = (vn1 + vn2)**2 / den
        if den == 0:
            df = 1

    # If df is undefined, variances are zero (assumes n1 > 0 & n2 > 0).
    # Hence it doesn't matter what df is as long as it's not NaN.
    df = np.where(np.isnan(df), 1, df)
    denom = np.sqrt(vn1 + vn2)
    return df, denom


def _ttest_ind_from_stats(mean1, mean2, denom, df):

    d = mean1 - mean2
    with np.errstate(divide='ignore', invalid='ignore'):
        t = np.divide(d, denom)
    t, prob = _ttest_finish(df, t)

    return (t, prob)


def _ttest_finish(df, t):
    """Common code between all 3 t-test functions."""
    dist_gen = stats.distributions.t
    prob = dist_gen.sf(np.abs(t), df) * 2  # two-sided, use np.abs to get upper alternative

    t_isnan = np.isnan(t)
    if np.any(t_isnan):  # and __scipy_prior0101:
        # older scipy's would return 0 for nan values of the argument
        # which is incorrect
        if np.isscalar(prob):
            prob = np.nan
        else:
            prob[t_isnan] = np.nan

    if t.ndim == 0:
        t = np.asscalar(t)

    return t, prob
