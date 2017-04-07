import numpy as np
import TTest as T
import scipy.stats as stats
#import matplotlib.pyplot as plt 


#Permutation bast estimation of p-value (#P'<p)/R-1
#Why is this giving me too low p-vlaues?
def Permutation(p,loc,seq,offset,test_loc=T.T_test_loc):
    no = len(seq)*int(np.log(1./p)+1.)
    if no > 10000:
        no = 10000
    #print np.log(1./p), no
    counts = 1.
    #p_vals = []
    for i in range(no):
        r = np.random.RandomState(i)
        tseq = r.permutation(seq)
        #p_vals.append(test_loc(tseq,loc) )
        if test_loc(tseq,loc) <= p:
            counts += 1.
    #plt.hist(p_vals)
    #plt.show(block=True)
    #plt.clf()
    if counts > 1.:
        #print "Permutation correction"
        p_corr = 1. * counts / no
    else:
        #print "Bonferroni correction"
        p_corr = p*(len(seq)**len(loc)) #If there are no counts, take bonferroni correction
    return p_corr

