import numpy as np 
from scipy import stats

#Other optimization objective: Log-likelihood

def negative_ll(seq):
    s = np.std(seq)
    return (np.log(2*np.pi) + np.log(s**2))*len(seq)/2.


def LL_test_loc(seq,loc,lookup='none',offset=2,p='no',Reg=False):
    p_test = []
    p_out = []
    
    index = 0

    seq_div = []
    
    while index <= len(loc): #intervals defined by loc. Take them and compare them seq[k],seq[l],seq[m]
        if index == 0:
            first = 0
        else:
            first = loc[index-1]
        if index == len(loc):
            second = len(seq)
        else:
            second = loc[index]
        
        if second-first < offset:
            print len(seq)
            raise NameError('The minimum spacing conditions are not fulfilled for the testd loc '+str(loc))
        
        seq_div.append(seq[first:second])
            
        index = index + 1


    nll = 0.    
    for s in seq_div:
        if len(np.unique(s)) < 2:
            s = list(s)
            s.append(1.)
        p_test.append(negative_ll(s))

    if p=='no':
        return sum(p_test)
    else:
        return p_test

