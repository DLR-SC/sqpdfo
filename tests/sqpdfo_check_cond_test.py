# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 11:53:25 2014

@author: jaco_da
"""

import unittest
from sqpdfo.sqpdfo_check_cond import *
from sqpdfo.runtime import compare_array, cond_
from sqpdfo import helper
from numpy import array


class Test_sqpdfo_check_cond(unittest.TestCase):
    """
      Reminder :
      This class is a test for sqpdfo_check_cond which check the the condition of the matrix does not exceed a certain threshold, and if it does,
      the matrix is modified so  that it does not anymore. The returned matrix is also symmetric when either original matrix is symmetric or when it has been modified.
    """ 
    def setUp(self):
        self.A = array([[ 1,  0, -1], [ 0,  1,  0], [ 1,  0,  1]])
        self.B = array([[1e12,0],[0,1e-8]])
        self.dummyOptions = helper.dummyOptions() 
        self.dummyOptions.verbose=0;
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;
        #print self.A

    def test_sqpdfo_check_cond_1(self):
        """
              test with a matrix chose condition number is sqrt(2). Results have been compared with matlab results : it is OK
        """
        res = sqpdfo_check_cond_(self.A, 1.41421356237, self.dummyOptions)
        self.assertTrue(res[1])
        self.assertFalse(compare_array(self.A, res[0], self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(res[0].T, res[0], self.abs_tol, self.rel_tol))
        res = sqpdfo_check_cond_(self.A, 1.41421356238, self.dummyOptions)
        self.assertFalse(res[1])
        self.assertTrue(compare_array(self.A, res[0], self.abs_tol, self.rel_tol))
        #print "np cond", np.linalg.cond(self.A)
        #print "sqpdfo  check cond", res

    def test_sqpdfo_check_cond_2(self):
        """
              test with a matrix chose condition number is 1e20 and which has very small singular values. Results have been compared with matlab results : it is OK
        """
        res = sqpdfo_check_cond_(self.B, 1e19, self.dummyOptions)
        self.assertTrue(res[1])
        self.assertFalse(compare_array(self.B, res[0], self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(res[0].T, res[0], self.abs_tol, self.rel_tol))
        res = sqpdfo_check_cond_(self.B, 1e21, self.dummyOptions)
        self.assertFalse(res[1])
        self.assertTrue(compare_array(self.B, res[0], self.abs_tol, self.rel_tol))


if __name__ == '__main__':
    unittest.main()