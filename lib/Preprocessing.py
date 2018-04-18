"""
Author: Thomas Gumbsch
Date: 16.3.2017
Lab: Machine learning and computational biology
Insitute: D-BSSE at ETHZ
"""
import numpy as np
import pandas as pd


def preprocessing(seq, transform):
    """
    Preprocessing of time series
    """
    if transform == 'diff':
        seq = diff(seq)
    elif transform == 'logdiff':
        seq = logdiff(seq)
    elif transform == 'percdiff':
        seq = percdiff(seq)
    elif transform == 'logpercdiff':
        seq = logpercdiff(seq)
    elif transform == 'std':
        seq = standardize(seq)
    return seq


def standardize(seq):
    mean_zero = [i - np.mean(seq) for i in seq]
    var_one = [i / np.std(seq) for i in mean_zero]
    return var_one


def diff(seq):
    """
    Take differences from sequence
    """
    return [i - j for i, j in zip(seq[1:], seq[:-1])]


def logpercdiff(seq):
    """
    Take log percentile differences from sequence
    """
    return [[np.log((np.abs(i - j) * (i + j) / 2.) + 1.) for i, j in zip(seq[1:], seq[:-1])]]


def logdiff(seq):
    """
    Take log differences from sequence
    """
    return [[np.log(np.abs((i - j)) + 1.) for i, j in zip(seq[1:], seq[:-1])]]


def precdiff(seq):
    """
    Take percentile diffrerence from sequence
    """
    return [[(i - j) * (i + j) / 2. for i, j in zip(seq[1:], seq[:-1])]]
