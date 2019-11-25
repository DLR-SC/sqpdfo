# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 14:55:53 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from sqpdfo_func import *
from sqpdfo_global_variables import set_prob
from numpy import array, double


class Test_sqpdfo_func(unittest.TestCase):
    """
          Reminder :
          This class is a test for sqpdfo_func which computes f, ci and ce according to a given problem
    """
    def setUp(self):
        self.x =  array([[ -1.717000000000000],   [1.595700000000000],   [1.827200000000000],  [-0.763600000000000],  [-0.763600000000000]])
        self.abs_tol=1e-11
        self.rel_tol=1e-11
        pass

    
    def test_sqpdfo_func(self):
        """
          Test with problem 5, results compared with matlab
        """
        set_prob(5)
        
        msg,f,ci,ce = sqpdfo_func_(self.x)
        
        correctmsg = 0
        correctf =0.053985698962084
        correctce = array([  -0.000822750000001, 0.000238240000000, 0.001195859492999])
        
        self.assertEqual(correctmsg, msg)
        self.assertAlmostEqual(correctf, double(f), 7)
        self.assertTrue(isempty_(ci))
        self.assertTrue(compare_array(correctce, ce, self.abs_tol, self.rel_tol))

if __name__ == '__main__':
    unittest.main()