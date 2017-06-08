# -*- coding: utf-8 -*-

import numpy as np
import copy
import itertools
from matplotlib import pyplot as plt
from os.path import join

import modules.TTest as T

import logging
loggerE = logging.getLogger(__name__)
formatterE = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
file_handlerE = logging.FileHandler(join('meta', 'Exact_log.log'))
file_handlerE.setFormatter(formatterE)
loggerE.addHandler(file_handlerE)
# loggerE.setLevel(logging.INFO)


def flatten(l):
    """
    flatten a list

    INPUT:
    l: list containing lists or tuples

    OUTPUT:
    flattened l

    >>> flatten([1,[2,3],(4,5,6), [7,[8,[9,10]]]])
    [1, 2, 3, 4, 5, 6, 7, 8, 9, 10]
    """
    out = []
    for item in l:
        if isinstance(item, (list, tuple)):
            out.extend(flatten(item))  # Recursive call in the case of item consisting of more than one dimension
        else:
            out.append(item)
    return out


def valid(n, loc, offset):
    """
    Indicator function wherther a configuration respecting the minimal spacing condition

    INPUT:
    n: sequence length
    loc: configuration of changepoints
    offset: miniumum interval size

    OUTPUT:
    True if the configuration is valid
    Remark: The changepoints are defined to be the beginning of a new interval

    >>> valid(9,[3,6],3)
    True
    >>> valid(8,[3,6],3)
    False
    """
    return all([i - j >= offset for i, j in zip(loc + [n], [0] + loc)])


# This checks all configurations for a given no
# Theoretically it suffices to scan for no=3 and imput good combinations from the results of all tested subintervals
def Combinatorical(seq, no, offset=3, test_loc=T.T_test_loc, Table='none', constraints=[], Reg=False, Fisher=True):
    """
    Check all configurations of no changepoints of constraints in seq and reports the most significant accoring to test_loc

    INPUT:
    seq: Time series
    no: m, number of changes
    offset: minimum interval size
    Table: lookup table for test_loc
    constraints: Candidates, if they have been extracted
    Reg: Regularization parameter if test_loc needs one
    Fisher: Indicator whether p-values are numerically displayable

    OUTPUT:
    best: configuration,test_loc(configuration)  best configuration of changepoints and its evaluated optimization objective
    """

    # If there are no contraints passed, explore all possible options
    if constraints == []:
        loggerE.warning('Runtime warning: All possible configruations checked')
        constraints = range(len(seq))

    # Check if conditions for running are satisfied
    if (no + 1) * offset > len(seq):
        loggerE.critical('Minimum spacing not possible, len(can)= ' + str(len(constraints)))
    elif no > len(constraints):
        loggerE.critical('You look for more changes that there are candidates.')
    else:

        best = [[], 1.]  # Initialization of running best config
        for loc in itertools.combinations(constraints, no):
            loc = list(loc)  # convert tuple to list

            # Only consider valid configurations (i.e. where the minimum spacing by offset if respected)
            if valid(len(seq), loc, offset) is True:
                # Update if a configuration gives a better score
                if Fisher:  # For p-values, compare the non-Bonferroni corrected
                    p = test_loc(seq, loc, Table, offset, 'no', Reg)
                    if p < best[1]:
                        best = [copy.copy(loc), p]
                else:  # For non-displayable combination of p-values, compare test statistics
                    p = 2. * np.sum([np.log(x) for x in test_loc(seq, loc, Table, offset, 'yes', Reg)])
                    if p > best[1]:
                        best = [copy.copy(loc), p]
            else:
                loggerE.info('Non-valid configuration ' + str(loc))

        return best[0], best[1]
