
import numpy as np
import Hausdorff as h
from DM import dissimilarity_manifold as DM
from os.path import join

import logging
loggerS = logging.getLogger(__name__)
formatterS = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
file_handlerS = logging.FileHandler(join('meta', 'Scale_log.log'))
file_handlerS.setFormatter(formatterS)
loggerS.addHandler(file_handlerS)
# loggerS.setLevel(logging.INFO)


def Estimate_Scale(seq, test_loc, Table):
    """
    Get the scale at whcih dissimilarity surface if of lowest dimension

    INPUT:
    seq: time series
    test_loc: Compute optimization obejctive  from seq and a configuration of cp (i.e. Welch-Fisher)
    Table: For Welch-Fisher, precomputed lookup table for test_loc

    OUTPUT:
    scale at which Box counting dimension of DM is lowest
    True/False if the computation of p-values failed at one point
    """
    window_exp = int(np.log2(len(seq)) - 1.5)
    dim = [502, 501, 500]
    cheat = [1., 1., 1.]
    delta = 1
    window = [int(2**(window_exp)), int(2**(window_exp))]

    while window_exp > 2:  # not(dim[-2] < dim[-1]) and window_exp > 2:

        stat_seq, P_value_cheat = DM(seq, 2**window_exp, test_loc, Table, 2)
        local = h.Box_counting(stat_seq)
        dim.append(local)
        cheat.append(P_value_cheat)
        window_exp = window_exp - delta
        window.append(int(2**(window_exp)))

    window.append(int(2**(window_exp - delta)))
    i = np.argmin(dim)
    loggerS.info(dim, i)

    if window[i] < 10:
        loggerS.warning('warning: Window size small.')
    elif window[i] > 1000:
        loggersS.warning('warning: Window size large.')
    return window[i], all([i for i in cheat])
