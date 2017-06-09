"""
Author: Thomas Gumbsch
Date: 16.3.2017
Lab: Machine learning and computational biology
Insitute: D-BSSE at ETHZ
"""
import numpy as np
import pandas as pd
import argparse


from lib.MaChaMP.MaChaMP import MaChaMP
from lib.E_Divise.E_Divise import E_Divise
from lib.PELT.PELT import PELT
from lib.WBS.WBS import WBS
from lib.SMUCE.SMUCE import SMUCE
from lib.BCP.BCP import BCP
from lib.CPM.CPM import CPM as Lepage


def PyChange(seq, transform='none', method='MaChaMP'):
    """
    Changepoint detection of input sequence
    """
    tseq = preprocessing(seq, transform)
    cp = solve(tseq, method)
    return cp


def solve(seq, method):
    """
    Apply method to sequence
    """
    cp = []
    if method == 'MaChaMP':
        cp = MaChaMP(seq)
    elif method == 'PELT':
        cp = PELT(seq)
    elif method == 'WBS':
        cp = WBS(seq)
    elif method == 'SMUCE':
        cp = SMUCE(seq)
    elif method == 'E-Divise':
        cp = E_Divise(seq)
    elif method == 'BCP':
        cp = BCP(seq)
    elif method == 'Lepage':
        cp = Lepage(seq)
    else:
        print "Not a known module"
        cp = [0]
    return cp


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
    return seq


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


def init_random_csv():
    """
    Sample csv data set
    """
    d = pd.DataFrame({'A': np.concatenate((np.cumsum(np.random.randn(50)), np.cumsum(np.random.randn(50) + 3), np.cumsum(np.random.randn(100)))), 'B': np.concatenate((['C1'] * 100, ['C2'] * 100)), 'T': range(200)})
    d.to_csv('random.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description='Computes changes in time sereis with python with different methods')
    parser.add_argument('--filename', help='Filename with .csv extension as from QtFy.')
    parser.add_argument('--cell', help='Cellid column name, separate time series have diffrent cell ids.')
    parser.add_argument('--values', help='Expression level column name')
    parser.add_argument('--time', help='timepoint column name')
    parser.add_argument('--method', choices=['MaChaMP', 'PELT', 'WBS', 'SMUCE', 'E-Divise', 'BCP', 'Lepage'], default='MaChaMP', help='Changepoint detection method')
    parser.add_argument('--preprocessing', choices=['none', 'diff', 'logdiff', 'percdiff', 'logpercdiff'], default='none', help='transformation of time series')
    args = parser.parse_args()
    name = str(args.filename)
    cell = str(args.cell)
    TF = str(args.values)
    method = str(args.method)
    time = str(args.time)
    transform = str(args.preprocessing)

    # init_random_csv()

    data = pd.read_csv(name)
    cells = data[cell].unique()
    changes = pd.DataFrame({'CellID': [], 'CP': []})

    for c in cells:
        seq = data[data[cell] == c][TF].dropna().values.tolist()
        cp = PyChange(seq, transform, method)
        print len(seq), cp
        changes = pd.concat((changes, pd.DataFrame({'CellID': [c] * len(cp), 'CP': cp, 'Timepoint': [data[(data[cell] == c) & (data[TF] == seq[change])][time].values.tolist()[0] for change in cp]})), ignore_index=True)

    changes.to_csv('Changes' + name)

    # Iterate over all unique Cell ids, get time series and apply Test
    # save in changes
