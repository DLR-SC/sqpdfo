# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 10:45:22 2014
INPUT VALUES:

info = 

        g: []
       ai: []
       ae: []
       hl: []
    niter: 0
                
                options = 

           algo_method: 'quasi-Newton'
    algo_globalization: 'trust regions'
           hess_approx: 'model'
          bfgs_restart: 0
          algo_descent: 'Powell'
                   tol: [1.000000000000000e-05 1.000000000000000e-05 1.000000000000000e-05]
                 dxmin: 1.000000000000000e-06
                 miter: 500
                msimul: 500
               verbose: 2
OUTPUT VALUES.

info = 

        g: []
       ai: []
       ae: []
       hl: []
    niter: 0
     flag: 0


options = 

           algo_method: 101
    algo_globalization: 112
           hess_approx: 131
          bfgs_restart: 0
          algo_descent: 120
                   tol: [1.000000000000000e-05 1.000000000000000e-05 1.000000000000000e-05]
                 dxmin: 1.000000000000000e-06
                 miter: 500
                msimul: 500
               verbose: 2
                  fout: 1
                   inf: Inf
                   df1: 0


values = 

                       success: 0
              fail_on_argument: 1
               fail_on_problem: 2
                 fail_on_simul: 3
                 stop_on_simul: 4
              stop_on_max_iter: 5
             stop_on_max_simul: 6
                 stop_on_dxmin: 7
          fail_on_non_decrease: 8
            fail_on_ascent_dir: 9
           fail_on_max_ls_iter: 10
              fail_on_ill_cond: 11
    stop_on_small_trust_region: 15
             fail_on_null_step: 20
         fail_on_infeasible_QP: 21
          fail_on_unbounded_QP: 22
                  fail_strange: 99
                    nsimultype: 16
                max_null_steps: 1
                        newton: 100
                  quasi_newton: 101
            cheap_quasi_newton: 102
                 unit_stepsize: 110
                    linesearch: 111
                 trust_regions: 112
                        powell: 120
                         wolfe: 121
                          bfgs: 130
                         model: 131

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from sqplab_options import *
from runtime import *
#import numpy as np
import helper

class dummyInfo():
    def __init__(self):

        self.g= matlabarray([])
        self.ai= matlabarray([])
        self. ae= matlabarray([])
        self.hl= matlabarray([])
        self.niter= 0

class Test_sqplab_options(unittest.TestCase):
    """
      Reminder :
      This class is a test for sqplab_options which sets the options of the optimization solver 'ecdfo'
    """
    def setUp(self):
        self.options = helper.dummyUnionStruct()
        self.options.algo_method = 'newton'
        self.options.algo_descent = 'powell'
        self.options.hess_approx = 'bfgs'
                
        self.info = dummyInfo()
        #self.values = helper.dummyValues()


    def test_sqplab_options1(self):
        """
        First little test
        """
        info,options,values = sqplab_options_(self.info,self.options)
        #print "algo_method, values", self.options.algo_method, values.newton
        self.assertEqual(options.algo_method, values.newton)
    
    def test_sqplab_options2(self):
        """ 
        Second little test
        """
        #print "---------------TEST 2-----------------"
        self.options.algo_globalization = 'trust regions'
        info,options,values = sqplab_options_(self.info,self.options)
        
        self.assertEqual(options.algo_descent, values.powell)
        self.assertEqual(options.hess_approx, values.bfgs)
        

if __name__ == '__main__':
    unittest.main()