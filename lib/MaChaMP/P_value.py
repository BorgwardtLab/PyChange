import numpy as np
import modules.TTest as T
import matplotlib.pyplot as plt
from os.path import join

import logging
loggerP = logging.getLogger(__name__)
formatterP = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
file_handlerP = logging.FileHandler(join('meta', 'Permutation_log.log'))
file_handlerP.setFormatter(formatterP)
loggerP.addHandler(file_handlerP)
# loggerP.setLevel(logging.INFO)


# Permutation bast estimation of p-value (#P'<p)/R-1
# Why is this giving me too low p-vlaues?
def Permutation(p, loc, seq, offset, test_loc=T.T_test_loc, Reg=False):
    """
    Estimate p-value of a configurtaion of changepoint with Westfall-Young

    INPUT:
    p: Threhold given as evalutation of optimization objective
    loc: configruation of changes
    seq: time series
    offset: Minimum inetrval size
    test_loc: Compute optimization objective for a single configuration without bonferroni
    Reg: regularization parameter for test_loc, if needed

    OUTPUT:
    Permutation corrected p
    """

    # Inizialization:
    no = 1000 * (len(seq)**len(loc))      # The number of tests
    # if no > 100000:
    #     no = 100000
    counts = 1.                                    # How many tests performed better
    maximum = p                                # Worst test result, to normalize if no test result is better than p

    # Run no permutations of sequence and test loc
    for i in range(no):
        r = np.random.RandomState(i)
        tseq = r.permutation(seq)
        test = test_loc(tseq, loc, 'none', offset, 'no', Reg)
        if test <= p:
            counts += 1.
        if test > maximum:
            maximum = test

    # Post-processing
    # If there is at least one configuration better than loc on seq
    if counts >= 1.:
        loggerP.info('Permutation test worked with counts = ' + str(counts) + ' from no_tests= ' + str(no))
        p_corr = 1. * counts / no
    # elif Reg == False:
    #     print loc, "Bonferroni correction"
    #     p_corr = p*(len(seq)**len(loc)) #If there are no counts, take bonferroni correction
    #     print p, p_corr
    # Otherwise, normalize to "no signal" value and report that
    else:
        loggerP.warning('Permutation test failed with p = ' + str(p) + ' and maximum = ' + str(maximum))
        p_corr = p / maximum

    return p_corr
