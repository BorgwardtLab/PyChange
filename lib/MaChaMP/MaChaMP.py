"""
Author: Thomas Gumbsch
Date: 20.5.2017
Email: thomas.gumbsch@bsse.ethz.ch
"""

import numpy as np
import scipy.misc as misc

import P_value as P
import Exact as Ex

import EstimateScale.Scale as S
from EstimateScale.DM import dissimilarity_manifold as DM


import modules.KSTest as KS
import modules.TTest as T
import modules.MWTest as MW
import modules.ADTest as AD
import modules.LLTest as LL
import modules.ETest as E
import modules.MSETest as M

from os.path import join


import warnings
warnings.filterwarnings("ignore")

import logging
logger = logging.getLogger(__name__)
formatter = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
file_handler = logging.FileHandler(join('meta', 'MaChaMP_log.log'))
file_handler.setFormatter(formatter)
logger.addHandler(file_handler)
logger.setLevel(logging.INFO)


def threshold_min(stat_seq, window, offset, Fisher=True, constraints=[]):
    """
    Find the lgn maxima in a sequence.

    INPUT:
    stat_seq: the sequence
    window: the window size that generated the sequence of test statistics
    offset: for the computation of a single test statistic
    Fisher: The threshold is set to 99th percentile if Fishers method fails
    constraints: Restricted region of sequence all changepoints are constrained to be confined in

    OUTPUT:
    A tuple of
    - a list containing the sorted lg2 indices and
    - the largest (lowest p-value) in that sequence

    >>> threshold_min([5,5,5,5,5,5,5,5,10,5,5,5,5,5],0,0)
    ([8], 10)
    >>> threshold_min([3,3,3,3,3,10,3,3,3,3,3,11,3,3,3],3,2)
    ([10, 16], 11)
    >>> threshold_min(list(np.random.randn(100))+[1000],0,0,False)
    ([100], 1000)
    """

    # If the p-values with Fishers methods are too low to be displayed, we increase the threshold. This is a hack for large changes
    if Fisher:
        thres = np.percentile(stat_seq, 50)
    else:
        thres = np.percentile(stat_seq, 99)
    # If no contraints are passed, all indices are allowed
    if len(constraints) <= 1:
        constraints = range(len(stat_seq) + 2 * window + 2 * offset)

    # Variable definition:
    min_index = []              # The list of indices from local minima (in sequence space)
    min_values = []            # The list of values corresponding to that list (in dissimilarity manifold)
    i = 1                             # counter variable
    candidate = [thres, 0]  # Running local maximum with [value,index]
    # Print to file: Setup for finding local minima
    logger.debug('Local minima finding: Window=' + str(window) + '; offset=' + str(offset) + '; threhold=' + str(thres) + ';')

    # One pass through the sequence
    while i < len(stat_seq):

        # Update running maximum
        if stat_seq[i] > candidate[0]:
            candidate = [stat_seq[i], i]

        # If threshold is surpassed, store and reset current local maximimum.
        if stat_seq[i] <= thres or i == len(stat_seq) - 1:
            if candidate[0] > thres:
                # Note all local maxima
                logger.debug('The list of maxima gets appended by ' + str(candidate) + ' at index ' + str(i))
                min_index.append(candidate[1] + window + offset)
                min_values.append(candidate[0])
            candidate = [thres, 0]

        # Increase counter
        i = i + 1

    # Postprocessing of list of local maxima: Rewrtite this to be iof order O(n)
    number_of_cases = int(np.log2(len(stat_seq) + 2 * offset + 2 * window) + 1.5)  # The number of maxima allowed to be reported
    index_of_min_index = sorted(range(len(min_values)), key=lambda k: -min_values[k])  # Sort values of all local maxima by their significance.
    return_min_index = [min_index[item] for item in index_of_min_index]  # Sorted list of indices of local amxima in time series space
    return_min_index = [item for item in return_min_index if item in constraints]  # Only take valid candidates from contraints
    # Prepare return list: Only top lg2n local maxima
    if number_of_cases <= len(return_min_index):
        candidates = return_min_index[:number_of_cases]
    else:
        candidates = return_min_index

    # Log the performance of the method and return the list of indices
    logger.info('From the ' + str(number_of_cases) + ' possible extrama ' + str(len(candidates)) + ' were found.')
    logger.info('These candidates are ' + str(candidates) + ' with maximum ' + str(max(stat_seq)))
    return sorted(candidates), max(stat_seq)


# Alternative extraction method:
# sort local maxima by their topological isolation and report them in this order
def topographic_isolation(stat_seq):
    pass


def MaChaMP(seq, returntype='changes', test_loc=T.T_test_loc):
    """
    Find the most significant configuration of changepoints in linear time

    INPUT:
    seq: Time series, a list

    returntype: what MaChaMP returns. Possible choices:
    'changes': The list of changepoints, default
    'both': Tulpe of list fo changepoints and overall combined p-value
    'debug': Tuple of list fo changepoints, the scale, the significance of changes and the threshold for reporting changepoints

    test_loc: test statistic. Possible choices:
    T.T_test_loc: Default for fast computation is Welchs t-Test
    MW.MW_test_loc: Mann-Whitney nonparametric test
    KS.KS_test_loc: Kolmogorov-Smirnov
    AD.AD_test_loc: Anderson darling distributional test (beta)
    LL.LL_test_loc: Log-likelihood penalty function, no p-values (beta)


    OUTPUT:
    List of changepoints

    >>> np.abs(sum(MaChaMP(list(np.random.randn(100))+list(np.random.randn(50)+3.)+list(np.random.randn(100))))-250) < 5
    True
    """

    # Create lookup Table for fast computation of test statsitcis
    Table = T.lookup_table(seq)
    # find scale at which dissimilarity manifold is of lowest dimension
    window, cheat = S.Estimate_Scale(seq, test_loc, Table)

    # Did the numerics of computing p-value break in the extraction of the dissimilarity manifold?
    Fisher = True
    if cheat != 1.:
        # Do not perform Fishers method, rather work in the test-statistic space.
        Fisher = False

    # Extract configuration of changepoints taking interdependence into account
    # For using log-likelyhood, permutation and regularization parameetrs have to be set to true
    Locs, p__, p_tilde = P_tilde(seq, test_loc, window, 2, Fisher, Table, [], False, False)

    # Return changes. Also possible returning scale, significance of changes and trhehold for reporting
    if returntype == 'changes':
        return Locs  # , window, p__, p_tilde
    elif returntype == 'both':
        return Locs, p__
    elif returntype == 'debug':
        return Locs, window, p__, p_tilde
    else:
        return Locs


def P_tilde(seq, test_loc=T.T_test_loc, window=20, offset=2, Fisher=True, Table='none', constraints=[], permutation=False, Reg=False):
    """
    Solve the universal optimization objective, given a scale

    INPUT:
    seq: time series, a list

    test_loc: test statistic. Possible choices:
    T.T_test_loc: Default for fast computation is Welchs t-Test
    MW.MW_test_loc: Mann-Whitney nonparametric test
    KS.KS_test_loc: Kolmogorov-Smirnov
    AD.AD_test_loc: Anderson darling distributional test (beta)
    LL.LL_test_loc: Log-likelihood penalty function, no p-values (beta)

    window: The scale at which the dissimilarity manifold has lowest dimension
    offset: Minimum number of data points in an interval requried to compute a test statistic
    Fisher: Did the numeric of computing p-values work for EstimateScale?
    Table: Lookup table, mean and variance of all subintervals starting at index 0

    constraints: Additional restriction of solution to subintervals
    []: No additional constraints, default
    [a,b]: All reported changes lie between a and b
    [[a1,a2],[b1,b2]]: All reported changes are between a1 and a2 or between b1 and b2

    permutation: Computing p-values for the objective
    False: Use Bonferroni correction to p-values, deault
    True: Compute p-values of confirguation with a permutation test, use if non test-based objectives aer used

    Reg: Input of regularization parameter to non test-based objective, if used

    OUTPUT:
    loc: The most significant configuration of changes
    list_p: List of significance of all single changes
    thres: p-value trheshold of reporting p-values

    >>> np.abs(sum(P_tilde(list(np.random.randn(100))+list(np.random.randn(50)+3.)+list(np.random.randn(100)))[0])-250) < 5
    True
    """

    # Initialization: Constraints, lookup table
    # Set additional constraints to location of changepoints, if they are specified
    if constraints == []:
        constraints = range(len(seq))
    else:
        if type(constraints[0]) == int:
            constraints = [constraints]
        constraints = [range(i[0], i[1]) for i in constraints]
        constraints = Ex.flatten(constraints)
    # Compute lookup table if that has not been done before
    if test_loc == T.T_test_loc and Table == 'none':
        Table = T.lookup_table(seq)
    elif test_loc != T.T_test_loc:
        Table = 'none'

    # Extract Candidates with sliding window:
    # Compute negative dissimilarity surface
    stat_seq, placeholder = DM(seq, window, test_loc, Table, offset)
    stat_seq = [-np.log(s) for s in stat_seq]
    # Extract local minima and threshold from dissimilarity surface
    constraints, thres = threshold_min(stat_seq, window, offset, Fisher, constraints)
    thres = np.exp(thres)  # convert to a p-value

    # Loop variables for iterationg over m:
    if Fisher:
        list_p = [thres]        # The ubber bound threshold for reporting is stored as the significance of m=0
    else:                           # i.e. there is at least one change
        list_p = [10**300]
    list_loc = [[]]                # The reported most significant configruations from m=0 to m=len(constriants)
    m = 1                          # Loop over number fo changes

    while m <= len(constraints):  # Iterate over all m
        loc_m, p_tilde = Ex.Combinatorical(seq, m, offset, test_loc, Table, constraints, Reg, Fisher)  # get best solution at m
        number_of_tests = Combinations(len(seq), m, offset)  # Number of possible configurations of m changes in seq with offset
        p_tilde *= number_of_tests  # apply Bonferroni correction
        list_p.append(p_tilde)
        list_loc.append(loc_m)
        m += 1
    logger.info('The solutions for all m are ' + str(list_loc) + ' with p-values ' + str(list_p))

    # Post-processing of all solutions for m changes
    best_m = np.argmin(list_p)  # The best number of changes is the one where the corrected p-value is lowest
    if any(np.isinf(list_p[1:])):  # If the configuration of changes gives numericalle too large values, just report all candidates
        best_m = -1
    return list_loc[best_m], list_p[best_m], thres


def Combinations(n, m, offset):
    """
    Possible configurations of m change with minimal padding offset in a sequence of n=length

    INPUT:
    n: number of data points in time series
    m: number of changes
    offset: Minimum number of data points in an interval requried to compute a test statistic

    OUTPUT:
    The number of possible combinations placing m changepoints in a sequence of length n with offset as the minimial interval size

    >>> Combinations(9,2,3)
    1.0
    >>> Combinations(18,5,3)
    1.0
    """
    logger.info('Correction for ' + str(m) + ' changes in a sequence with length ' + str(n) + ' with omega = ' + str(offset) + ' yields ' + str(misc.comb(n - (m + 1) * offset + m, m)))
    return misc.comb(n - (m + 1) * offset + m, m)
