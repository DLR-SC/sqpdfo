# -*- coding: utf-8 -*-
"""
Created on Fri Nov 07 11:52:38 2014

@author: jaco_da
"""
import sys
sys.path.append("../")

import unittest
from blls import *
#import helper

class Test_blls(unittest.TestCase):
    """
      Reminder :
      This class is a test for blls which is a solver for bound-counstrained linear least-squares problems.
    """ 
    def setUp(self):
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;

    def test_blls_1(self):
        """
              test with  some random data. Results have been compared with matlab results : it is OK
        """
        self.A = matlabarray([[6,     0,     4],
                [1,     5,     7],
                 [0,     2,     2],
                 [3,     7,     0]])

        self.b = matlabarray([1, 5, 4, 3]).T

        self.lb = matlabarray([-1, -1, -1]).T
        self.ub = -self.lb
        
        s, resn, opt, exitc = blls_(self.A,self.b, self.lb, self.ub)

        correctS =matlabarray([ -0.171259842519685, 0.516235951275322, 0.466316710411199]).T
        self.assertTrue(compare_matlabarray(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(resn,2.152108290072786, places=10)
        self.assertAlmostEqual(opt,3.631770138003687e-14, places=10)
        self.assertEqual(exitc, 0)
         

if __name__ == '__main__':
    unittest.main()