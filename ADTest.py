import numpy as np 
from scipy import stats

#Anderson-Darling for shape of Distributions, works not very well....


def AD_test_loc(seq,loc,lookup='none',offset=2,p='no'):
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
        
    for s1,s2 in zip(seq_div[:-1],seq_div[1:]):
        if len(np.unique(s1)) < 2 or len(np.unique(s2)) < 2:
            s1 = list(s1)
            s2 = list(s2)
            s1.append(1.)
            s2.append(1.)
        p_test.append(stats.anderson_ksamp([s1,s2])[2])

    if p=='no':
        return stats.chisqprob(-2.*np.sum([np.log(x) for x in p_test]),2*len(loc))
    else:
        return p_test

