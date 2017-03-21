"""
Author: Thomas Gumbsch
Date: 16.3.2017
Lab: Machine learning and computational biology
Insitute: D-BSSE at ETHZ
"""
import numpy as np
import pandas as pd
import TTest_combinatorical as TC
import Exact as EX
import TTest_divide as TD
import TTest_ANOVA as TA
import argparse

def solve(seq,method='single_diff'):
	if method == 'single':
		cp, p = TD.TTest(seq)
		cp = [cp]
		p = p * len(seq)
	elif method == 'single_diff':
		cp, p = TD.TTest(diff(seq))
		cp = [cp]
		p = p * len(seq)
	elif method == 'DP':
		cp = TD.mc_p(seq,0.05/len(seq))
		p = TC._test_loc(seq,cp,2)
	elif method == 'DP_diff':
		cp = TD.mc_p(diff(seq),0.05/len(seq))
		p = TC._test_loc(seq,cp,2)
	elif method == 'exact':
		cp,p = EX.Exact(seq)
	elif method == 'exact_diff':
		cp,p = EX.Exact(diff(seq))
	elif 'exact' in method:
		no = int(method[-1])
		if 'diff' in method:
			seq = diff(seq)
		cp,p = EX.Exact(seq,no)
	elif 'ANOVA' in method:
		if 'diff' in method:
			seq = diff(seq)
		cp = TA.Optimal_ANOVA(seq,0.05,2)
		p = TC._test_loc(seq,cp,2)
	else:
		print "Not a known module"
		cp = [0]
		p = 1.

	if p > 1:
		p = 1.
	return cp, p

def diff(seq):
	return [i-j for i,j in zip(seq[1:],seq[:-1])]


def init():
	d = pd.DataFrame({'A': np.concatenate((np.cumsum(np.random.randn(50)),np.cumsum(np.random.randn(50)+3), np.cumsum(np.random.randn(100)))), 'B': np.concatenate((['C1']*100, ['C2']*100)), 'T': range(200)})
	d.to_csv('random.csv')

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description = 'Computes changes in time sereis with python with different methods')
    parser.add_argument('Filename', help='Filename with .csv extension as from QtFy.')
    parser.add_argument('Cell', help='CellId column name')
    parser.add_argument('TF', help='Expression level column name')
    parser.add_argument('Time', help='timepoint column name')
    parser.add_argument('Type', help='Standard: single_diff. Type of method: single, DP, ANOVA and exact, possible with _diff for differences. For exact also possible exact_diff_3 for three changes on differences')
    args = parser.parse_args()
    name = str(args.Filename)
    cell = str(args.Cell)
    TF = str(args.TF)
    method = str(args.Type)
    time = str(args.Time)

    #init()

    data = pd.read_csv(name)
    cells = data[cell].unique()
    changes = pd.DataFrame({'CellID': [], 'CP': [], 'P': []})

    for c in cells:
    	seq = data[data[cell] == c][TF].dropna().values.tolist()
    	cp, p = solve(seq,method)
    	changes = pd.concat((changes,pd.DataFrame({'CellID': [c]*len(cp), 'CP': cp, 'P': [p]*len(cp), 'Timepoint': [data[(data[cell]==c) & (data[TF] ==seq[change])][time].values.tolist()[0] for change in cp]})),ignore_index=True)

    changes.to_csv('Changes'+name)

    #Iterate over all unique Cell ids, get time series and apply Test
    # save in changes





