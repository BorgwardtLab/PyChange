import numpy as np
from os.path import join

import logging
loggerD = logging.getLogger(__name__)
formatterD = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
file_handlerD = logging.FileHandler(join('meta', 'DM_log.log'))
file_handlerD.setFormatter(formatterD)
loggerD.addHandler(file_handlerD)
# loggerD.setLevel(logging.INFO)


def dissimilarity_manifold(seq, window, test_loc, lookup, offset=2):
    """
    Transform sequence into dissimilarity manfold (DM).

    INPUT:
    seq: consists of data points
    window: window size / scale for DM
    offset: for the computation of a single test statistic
    test_loc: Function that returns $p$-value based on similarity after division into segements
    lookup: mean and variance of all subintervals starting at index 0

    OUTPUT:
    Dissimilarity manifold of all i in [window+offset,n-window-offset]
    Did the computation fo p-values fail at one point?
    """
    stat_seq = []
    # start at earliest value and end at last possible value
    for l in range(window + offset, len(seq) - window - offset):
        # Compute T(c_2) for c_1 = l-window-offset, c_2 = l, c_3 = l+window+offset
        stat_seq.append(test_loc(seq, [l - window, l, l + window], lookup, offset, 'yes')[1] * (len(seq) - (2 * offset + 1)))  # Bonferroni correction
    # Log an indication of: How much does it make sense?
    # logger.info('The transformation yields ' + str(np.argmin(stat_seq) + window + offset) + ' as the location for the most likely single change.')

    # Exclude all values bigger than 1E300 (also inf). If there are any, Fishers method will fail.
    length_before_cut = len(stat_seq)
    stat_seq = [i for i in stat_seq if i < 1E300]
    length_after_cut = len(stat_seq)
    if length_after_cut == 0:
        length_after_cut = 1

    loggerD.info('The fraction of DM surviving after excluding very large values is ' + str(1. * length_before_cut / (length_after_cut)))
    loggerD.info('The transformation yields ' + str(np.argmin(stat_seq) + window + offset) + ' as the location for the most likely single change.')

    return stat_seq, (1. * length_before_cut / length_after_cut == 1.)
