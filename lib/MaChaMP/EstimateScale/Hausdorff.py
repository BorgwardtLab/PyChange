

import numpy as np
from os.path import join


import logging
loggerH = logging.getLogger(__name__)
formatterH = logging.Formatter('%(asctime)s:%(levelname)s:%(message)s')
file_handlerH = logging.FileHandler(join('meta', 'Hausdorff_log.log'))
file_handlerH.setFormatter(formatterH)
loggerH.addHandler(file_handlerH)
# loggerH.setLevel(logging.INFO)


def Boxlen(k, n, u, big_k):
    """
    Compute lengths of Box

    INPUT:
    number of data points n
    largest power of two that fits the time series big_K,
    variation u, level

    OUTPUT:
    xlen and ylen of box
    """
    level = 2**(k - big_k)
    return level * n, level * u


def get_u(seq):
    """
    Compute the  approximate variation of the sequence
    Note: We use here only the log2n data points which will be looked at later.

    INPUT: time series
    OUTPUT: max-min, min
    """

    # Exact log2 n data points from teh sequence
    r = np.random.RandomState(len(seq))
    ind = r.choice(len(seq), int(np.log2(len(seq)) + 1.5))
    logseq = [seq[i] for i in ind]

    # Get mean and std
    mean = np.mean(logseq)
    std = np.std(logseq)

    # minimum cannot be lower than 0, has to be at least the lowest value extracted and is otherwise 1.5 times the standar devaiation below mean
    m = min(min(logseq), max(mean - 1.5 * std, 0.))
    # The spread is at least what has been observed in the log2n sample data and is other wise 3 times the standard deviation
    u = max(logseq) - min(logseq)# max(np.abs(max(logseq) - min(logseq)), 3 * std)

    loggerH.info('The spread and variance ' + str(u) + ' ' + str(m) + ' have been chosen from the data ' + str(logseq))

    return u, m


def get_big_k(n):
    """
    Largest box-eponent possible

    Input: length of time series
    Output: largest power of two that fits into n

    >>> get_big_k(10)
    3
    """
    return int(np.floor(np.log2(n)))


def Box_index(i, j, x, y, m):
    """
    Find the index of the box from looking at one point of the time series

    INPUT:
    i: timepoint currently under observation
    j: seq[i] value at that point
    x: x-size of box
    y: y-size of box
    m: minimum of time series (lower bound on j)

    OUTPUT:
    x,y-index of Box corresponding to i,j

    >>> Box_index(4.973,99.923,2,2,95)
    (2, 2)
    """
    ep = 1E-10  # take care of rounding errors in y-dimension
    return int(np.ceil(1. * i / x) - 1), int(np.ceil((j - m) / y - ep) - 1)


def Ne(seq, k, n, big_k, u, m):
    """
    Get approximation of number of occupied boxes

    INPUT:
    seq: time series
    k: exponent of box size
    n: length of data
    big_k: Largest box-exponent possible
    u: spread of data
    m: minimum of seq

    OUTPUT:
    Number of occupied boxes N_epsilon
    """

    # For zero width, all Boxes are filled
    ep = 1E-10  # Take care of rounding errors in y-dimension
    if u == 0:
        x, y = Boxlen(k, n, u, big_k)
        return int(np.ceil(1. * n / x))
    else:
        x, y = Boxlen(k, n, u, big_k)
        # Inizialize a matrix of Box indices to zero occupation
        Boxes = np.zeros((int(np.ceil(1. * n / x)), int(np.ceil(1. * u / y + ep) + 1)))
        # Speedup: only look at log2n points in the sequence mark corresponding box as filled
        r = np.random.RandomState(len(seq))
        ind = r.choice(len(seq), int(np.log2(len(seq)) + 1.5))
        logseq = [seq[i] for i in ind]
        offset = min(logseq)
        logseq = [i - offset for i in logseq]
        for i, j in zip(ind, logseq):
            a, b = Box_index(i, j, x, y, m)
            loggerH.debug(a, b, x, y, m, u)
            if Boxes[a, b] == 0:
                Boxes[a, b] = 1

        loggerH.debug(Boxes)

        return 2**(np.sum(Boxes))


def D_Box(s, Nepsilon):
    """
    Compute slope of linear regression logN vs loge

    INPUT:
    s: list of logepsilon
    Nepsilon: list of number of occupied boxes correponidng to epsilon boxsixe

    OUTPUT:
    negative slope of linear regression between log(e) vs log(Ne)

    >>> D_Box([1,2,3],[2,4,8])
    -1.0
    """
    ssxm, ssxym, ssyxm, ssym = np.cov(s, np.log2(Nepsilon)).flat  # flat iterates over array
    return - ssxym / ssxm  # -temp_nom / temp_denom  # Linar regression slope alternatively ssxym/ssxm


def k_to_logepsilon(k, big_k, n):
    """
    Transform index of box-counting exponent to scaled log of box-size

    INPUT:
    k: index of box-counting exponent
    big_k: Largest box-exponent possible
    n: sequence length

    OUTPUT:
    scaled log of box-size
    """
    return np.log2(n * 2**(k - big_k))


def Box_counting(seq):
    """
    Monte carlo approximation to Box-counting algorithm

    INPUT:
    seq: time series

    OUTPUT:
    Box-counting dimension
    """

    # Initialization
    exclude_exponents = 2  # Recommended by Liebovic and Toth, exclude two largest and smallest Box sizes
    n = len(seq)  # Length of sequence
    u, m = get_u(seq)  # Get spread and minimum (determine it faster on the fly from log2n evaluated elements)
    big_k = get_big_k(n)  # Largest box size exponent
    listk = [big_k - exclude_exponents]  # List of evaluated box-size exponents
    logepsilon = [listk[-1] - big_k]  # In epsilon notation
    Nepsilon = [Ne(seq, listk[-1], n, big_k, u, m)]  # Get number of occupied Boxes

    # reduce Box size and count number of occupied boxes
    while len(Nepsilon) < big_k - 2 * exclude_exponents:  # Nepsilon[-1] < int(np.log2(n)): # formley n/4 for non-speedup
        listk.append(listk[-1] - 1)
        Nepsilon.append(Ne(seq, listk[-1], n, big_k, u, m))

    loggerH.info(Nepsilon)

    logepsilon = [k_to_logepsilon(i, big_k, n) for i in listk]

    return D_Box(logepsilon, Nepsilon)
