# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 11:26:49 2014
check if function returns convex matrix for convex matrix (1) and non convex matrix (2) and if matrix is non symmetric (3)
@author: jaco_da
"""

import sys
sys.path.append("../")
#sys.path.append("tests/")

import unittest
from sqpdfo_check_convex import *
import helper
import numpy as np
from runtime import compare_array, eig_
from random import random
from numpy import array

class Test_sqpdfo_check_convex(unittest.TestCase):
    """
      Reminder :
      This class is a test for sqpdfo_check_convex which check if the matrix is convex and if not modify it so it becomes convex.
      The returned matrix is also symmetric when either original matrix is symmetric or when it has been modified.
    """ 
    def setUp(self):
    #    self.dummyOptions = helper.dummyOptions()
        self.abs_tol=1e-15;
        self.rel_tol=1e-15;

    def test_sqpdfo_check_convex1(self):
        """
        With a convex matrix
        """
        A = array([[ 2,  -1, 0], [ -1,  2,  -1], [ 0,  -1,  2]])

        res = sqpdfo_check_convex_(A)
        #print "A", A
        #print "sqpdfo  check cond", str(res)

        self.assertTrue(compare_array(A, res, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(res, A, self.abs_tol, self.rel_tol))
#
    def test_sqpdfo_check_convex2(self): 
        """
        With a non convex matrix
        """
        B =  array([[ 1,  2], [ 2,  1]])

        res = sqpdfo_check_convex_(B)
        #print "B", B
        #print "sqpdfo  check cond", str(res)
        self.assertFalse(compare_array(B, res, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(res, res.T, self.abs_tol, self.rel_tol))
        
        matlabres=array([[3.000000000000000,  0.000000001000000]])
        self.assertTrue(compare_array(matlabres, eig_(res), 1e-09, 1e-09))
        
    def test_sqpdfo_check_convex3(self):
        """
        With a non symmetric matrix non convex
        """
        C = array([[1,2,3],[4,5,6],[7,8,9]])
        
        res = sqpdfo_check_convex_(C, helper.dummyOptions())

        self.assertFalse(compare_array(C, res, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(res, res.T, self.abs_tol, self.rel_tol))

        matlabres=array([[  16.116843969881806, 0.000000000925233,  0.000000001000001]])        
        self.assertTrue(compare_array(matlabres, eig_(res), 1e-8, 1e-8))
        
        

if __name__ == '__main__':
#    suite = unittest.TestLoader().loadTestsFromTestCase(Test_sqpdfo_check_convex)
#    unittest.TextTestRunner(verbosity=2).run(suite)
    unittest.main()
