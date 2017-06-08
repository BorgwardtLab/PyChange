from scipy import misc
import numpy as np

#E-Distance from http://download.springer.com/static/pdf/571/art%253A10.1007%252Fs00357-005-0012-9.pdf?originUrl=http%3A%2F%2Flink.springer.com%2Farticle%2F10.1007%2Fs00357-005-0012-9&token2=exp=1492672053~acl=%2Fstatic%2Fpdf%2F571%2Fart%25253A10.1007%25252Fs00357-005-0012-9.pdf%3ForiginUrl%3Dhttp%253A%252F%252Flink.springer.com%252Farticle%252F10.1007%252Fs00357-005-0012-9*~hmac=9a25132020641e9783ca873fc5c60f38f57f2f987b38c742a367686856e087d8
#used in Matteson and James (2014) equation (5) and (6)
def E_test_loc(seq,loc,lookup='none',offset=2,p='no',alpha=1):
    E = []
    
    intervals = [0] + loc + [len(seq)]
    
    for left,cp,right in zip(intervals[:-2],intervals[1:-1],intervals[2:]):
        n = cp-left
        m = right - cp #In the paper, they let right= [min,len(seq)]

        E1 = 0
        for i in range(n):
            for j in range(m):
                E1 += np.abs(seq[left+i]-seq[j+cp])**alpha
        E1 = E1*2./(n*m)

        E2 = 0
        for k in range(n):
            for i in range(k):
                E2 += np.abs(seq[left+i]-seq[left+k])**alpha
        E2 = E2/misc.comb(n, 2)

        E3 = 0
        for k in range(m):
            for j in range(k):
                E3 += np.abs(seq[cp+j]-seq[cp+k])**alpha
        E3 = E3/misc.comb(m, 2)

        E.append(    (-1)*m*n*(E1 - E2 - E3)/(m+n) )

    return sum(E)
