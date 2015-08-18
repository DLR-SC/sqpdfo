# -*- coding: utf-8 -*-
"""
Created on Tue Nov 11 11:22:07 2014

@author: jaco_da
"""
import sys
sys.path.append("../")

import unittest
from sqplab_bfgs import *
import numpy as np
import helper
from numpy import array

class dummyInfo():
    """
          Reminder :
          This class is a test for sqplab_bfgs which computes the BFGS update of the matrix M, which is
           supposed to be positive definite approximation of some Hessian. If
           y'*s is not sufficiently positive and options.algo_descent is set to
           values.powell, Powell's correction is applied to y
    """
    def __init__(self):
        self.g = array([[  -0.009783923878659,  1.000000000000000, 0.0]]).T
        self.ai = array([])
        self.ae = array([[1.000000000000000,   0.999999999999999,   1.000000000000000],
   [1.000000000000000,   1.999999999999999,   3.000000000000000]])
        self.hl = array([])
        self.niter = 1
        self.flag = 0
        self.nsimul = array([0, 5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.ci = array([])
        self.ce = array([[0.312467441560872, -0.056489622187450]]).T
        self.f = 0.099563123926544
        self.glag = array([[-0.333333333013890,  0.666666666736111, -0.333333333513889]]).T
        self.glagn = 0.666666666736111
        self.feasn = 3
        self.compl = 0

class Test_sqplab_bfgs(unittest.TestCase):

    def setUp(self):
        self.M = array(np.eye(3))
        self.y = array([[-0.009783923878659,  -0.0, 0.0]]).T
        self.s = array([[-0.503054026564966,  -1.0, -0.184478531874161]]).T
        self.first = 1
        
        self.info = dummyInfo()
        self.options = helper.dummyOptions()
        #set hessian approximation to bfgs
        self.options.hess_approx = 130
        self.values = helper.dummyValues()
        

    def test_sqplab_bfgs_null_step(self):
        """ 
        Test with a null step
        """
        self.s = array(np.zeros(3))
        _, _, info, values = sqplab_bfgs_(self.M,self.y,self.s,self.first,self.info,self.options,self.values)
        
        self.assertEqual(info.flag, values.fail_strange)
        
    def test_sqplab_bfgs_negative_definite(self):
        """
        Test with a negative definite step
        """
        #self.M = array([[1,2,3],[2,1,2],[3,2,1]])
        self.M = array(-np.eye(3))
        #print "M", self.M
        #self.s = array([1.0,1.0])
        _, _, info, values = sqplab_bfgs_(self.M,self.y,self.s,self.first,self.info,self.options,self.values)
                
        self.assertEqual(info.flag, values.fail_strange)
    
    def test_sqplab_bfgs_test(self):
        """
     Test to compare the results with matlab, using data from problem 3 options.hess_approx = 'bfgs'. All data (from info and values attributes)
     have been compared to matlab (though tests are not written) : results are OK
        """
        M, pc, info, values = sqplab_bfgs_(self.M,self.y,self.s,self.first,self.info,self.options,self.values)
        #print "M:\n", M
        correctM = array([[0.205243860649550,   0.003553460424288,   0.000655537162145],
   [0.003553460424288,   0.195307727402361,  -0.000901167227317],
   [0.000655537162145,  -0.000901167227317,   0.200026425015353]])

        self.assertAlmostEqual(0.803070936030031, pc, places=10)
        self.assertTrue(compare_array(M, correctM, 1e-12, 1e-12))
if __name__ == '__main__':
    unittest.main()