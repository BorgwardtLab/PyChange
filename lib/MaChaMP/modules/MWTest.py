import numpy as np 
from scipy import stats

#Easiest non-parametric test
#Better: http://download.springer.com/static/pdf/571/art%253A10.1007%252Fs00357-005-0012-9.pdf?originUrl=http%3A%2F%2Flink.springer.com%2Farticle%2F10.1007%2Fs00357-005-0012-9&token2=exp=1491238778~acl=%2Fstatic%2Fpdf%2F571%2Fart%25253A10.1007%25252Fs00357-005-0012-9.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Farticle%252F10.1007%252Fs00357-005-0012-9*~hmac=fb934cfba263a901715e224cf3ba3787078b8b2c4b5c96a4fbdb36d8a839d217


def MW_test_loc(seq,loc,lookup='none',offset=2,p='no',Reg=False):
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
        p_test.append(stats.mannwhitneyu(s1,s2)[1])

    if p=='no':
        return stats.chisqprob(-2.*np.sum([np.log(x) for x in p_test]),2*len(loc))
    else:
        return p_test

