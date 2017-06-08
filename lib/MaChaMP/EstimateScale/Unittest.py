import numpy as np 
from CP_Module import P_tilde
import modules.KSTest as KS
import modules.TTest as T 
import modules.MWTest as MW
import modules.ADTest as AD 
import modules.LLTest as LL
import P_value as P
import Exact as Ex
import modules.MSETest as M

import unittest
import scipy.stats as st
#"""

class TestComputation(unittest.TestCase):
    def test_table(self):
        r = np.random.RandomState(0)
        a = r.randn(10000)
        Table = T.lookup_table(a)
        for j in range(3,len(a)):  
            m,v =Table[str(j)]
            self.assertAlmostEqual(np.var(a[:j])*j/(j-1),v,places=3) #numpy is biased
            self.assertAlmostEqual(np.mean(a[:j]),m,places=3)
    def test_locs(self):
        r = np.random.RandomState(0)
        a = r.randn(100000)
        Table = T.lookup_table(a)
        for i in range(100):
            loc = r.choice(range(3,len(a)-3),1)[0]
            m,v,l = T.get_m_v_l(a,[loc,len(a)],Table)
            self.assertAlmostEqual(np.var(a[loc:]),v[1],places=3)
            self.assertAlmostEqual(np.mean(a[loc:]),m[1],places=3)




class TestTTest(unittest.TestCase):			
    def test_ttest(self):
         for i in range(10):
             #Check if some sequences are correct
             r = np.random.RandomState(i)
             a = list(r.randn(100))
             b = list(r.randn(10*i+5))
             self.assertAlmostEqual(T.T_test_loc(a+b,[100]),st.ttest_ind(a,b,equal_var=False)[1], places=1)
    def test_fisher(self):
        for i in range(10):
            r = np.random.RandomState(i)
            a = list(r.randn(100))
            b = list(r.randn(10*i+5))
            c = list(r.randn(100))
            first = st.ttest_ind(a,b,equal_var=False)[1]
            second =  st.ttest_ind(b,c,equal_var=False)[1]
            self.assertAlmostEqual(T.T_test_loc(a+b+c,[100,100+10*i+5]),st.chisqprob(-2.*np.sum([np.log(x) for x in [first,second]]),4) , places=1)

class TestExact(unittest.TestCase):
    def test_one(self):
        for i in range(10):
            r = np.random.RandomState(i)
            a = list(r.randn(100))
            b = list(r.randn(100)+10.)
            self.assertEqual(Ex.Combinatorical(a+b,1)[0],[100])
    def test_two(self):
        for i in range(10):
            r = np.random.RandomState(i)
            a = list(r.randn(50))
            b = list(r.randn(50)+20.)
            self.assertEqual(Ex.Combinatorical(a+b+a,2)[0],[50,100])
    #def test_hand_computed(self):
    #    a = [0.]*10+[0.1]+[0.]*10
    #    b = [1.]*10+[1.1]+[1.]*10
    #    self.assertEqual(Ex.Combinatorical(a+b+a,2)[0],[21,42])
    #    self.assertEqual(Ex.Combinatorical(a+b+a+b,3)[0],[21,42,63])

class TestP_tilde(unittest.TestCase):
    def test_one(self):
        for i in range(1):
            r = np.random.RandomState(i)
            a = list(r.randn(100))
            b = list(r.randn(100)+10.)
            self.assertEqual(P_tilde(a+b,permutation=True)[0],[100])
    #def test_hand_computed(self):
    #    a = [0.]*10+[0.1]+[0.]*2
    #    b = [1.]*10+[1.1]+[1.]*2
    #    self.assertEqual(P_tilde(a+b,permutation=True)[0],[13])
        
  
        



if __name__ == '__main__':
	unittest.main()