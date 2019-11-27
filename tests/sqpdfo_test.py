# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 17:23:18 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from sqpdfo_evalfgh import sqpdfo_evalfgh_
from runtime import compare_array, isempty_
from sqpdfo import sqpdfo_
from sqpdfo_global_variables import set_prob, set_threshold, set_fileoutput,set_simul_not_initialized, set_check_condition
import numpy as np
from numpy import array


class myOptions():
    def __init__(self):
        self.hess_approx = 'model'
        self.bfgs_restart = 0
        self.algo_descent = 'powell'
        self.tol_grad = 1.000000000000000e-06
        self.tol_feas = 1.000000000000000e-06
        self.tol_bnds = 1.000000000000000e-06
        self.dxmin = 1.000000000000000e-16
        self.miter = 500
        self.msimul = 500
        self.verbose = 0
        self.whichmodel = 'subbasis'
        self.final_degree = 'quadratic'


class Test_sqpdfo(unittest.TestCase):
    """
      Reminder :
      This class is a test for sqpdfo 
    """
    def setUp(self):

        set_prob(3)
        set_threshold(1e-08)
        set_fileoutput(1)
        set_simul_not_initialized(1)
        set_check_condition(0)
        self.options = myOptions()
        self.abs_tol=1e-10
        self.rel_tol=1e-10

        pass

    def test_sqpdfo(self):
        """
         Test which compare python and matlab results
         test runs with whichmodel = 0 in sqpdfo.py
        """
        x,lm,info = sqpdfo_(self.options)
        correctx = array([ -0.5, 0.0, 0.5])
        correctlm = array([[0,-0.000005713064576,0,2.,-1.]])

        self.assertTrue(compare_array(correctx, x, 1e-9, 1e-9))
        self.assertTrue(compare_array(correctlm, lm, 1e-5, 1e-5))
        self.assertEqual(info.niter,4)
        self.assertEqual(info.flag,0)
        self.assertAlmostEqual(info.f, 0.5, places=9)
        self.assertTrue(compare_array(info.g, array([-1.,0.,1.]),1e-5, 1e-5))
        self.assertTrue(compare_array(info.ae, array([[1.,1.,1.],[1.,2.,3.]]), 1e-9, 1e-9))
        self.assertTrue(isempty_(info.hl))
        self.assertAlmostEqual(2.220446049250313e-16,info.feasn,places=9)
        self.assertTrue(compare_array(info.ce, array([2.220446049250313e-16,-1.110223024625157e-16]),self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.glag, array([0.,0.,0.]), 1e-9, 1e-9))
        self.assertAlmostEqual(0.,info.glagn,places=8)
        

if __name__ == '__main__':
    unittest.main()