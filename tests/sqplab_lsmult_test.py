# -*- coding: utf-8 -*-
"""
Created on Fri Nov 07 15:30:36 2014
lsmult
K>> x

x =

   0.500000000000000
   1.000000000000000
   0.500000000000000

K>> lb

lb =

     []

K>> ub

ub =

     []

K>> info

info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 0
      flag: 0
    nsimul: [0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 1.500000000000000

K>> options

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

K>> info.g

ans =

     0
     1
     0

K>> info.ae

ans =

     1     1     1
     1     2     3

K>> info.ce

ans =

     2
     3

K>> values

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
                         dline: '--------------------------------------------------------------------------------------'
                         eline: '======================================================================================'
                         sline: '**************************************************************************************'
                                                                                                    
----------OUTPUT VALUES----------------
lm

lm =

                   0
                   0
                   0
  -0.333333332763891
  -0.000000000249999

K>> info

info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 0
      flag: 0
    nsimul: [0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 1.500000000000000

K>> infor.g
Undefined variable "infor" or class "infor.g".
 
Did you mean:
K>> info.g

ans =

     0
     1
     0

K>> info.ae

ans =

     1     1     1
     1     2     3

K>> info.ce

ans =

     2
     3                                                                                                
@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from sqplab_lsmult import *
import numpy as np
import helper

class dummyInfo():
    def __init__(self):
        #g: [3x1 double]
        self.g = matlabarray([0, 1, 0]).T
        self.ai = matlabarray([])
        #ae : [2x3 double]
        self.ae = matlabarray([[1,     1,     1],    [ 1 ,    2,     3]])
        self.hl = matlabarray([])
        self.niter = 0
        self.flag = 0
        self.nsimul = matlabarray([0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.ci = matlabarray([])
        #ce: [2x1 double]
        self.ce = matlabarray([2, 3]).T
        self.f = 1.5
                                                                                                
class Test_sqplab_lsmult(unittest.TestCase):
    """
      Reminder :
      This class is a test for sqplab_lsmult which  computes an exact least-squares multiplier 'lm'. It solves in lm the
       quadratic optimization problem:
           min || g+A'*lm ||^2
           subject to possible bounds on lm,
    """ 
    def setUp(self):
        self.x = matlabarray([0.5, 1.0, 0.5])
        self.lb = matlabarray([])
        self.ub = matlabarray([])
        self.values = helper.dummyValues()
        self.options = helper.dummyOptions()
        self.info = dummyInfo()
        self.abs_tol=1e-8;
        self.rel_tol=1e-8;

    def test_sqplab_lsmult(self):
        """
            Test comparing python results with the matlab results
        """
        lm,info = sqplab_lsmult_(x=self.x,lb=self.lb,ub=self.ub,info=self.info,options=self.options,values=self.values)
    
        correctlm = matlabarray([0, 0, 0, -0.333333332763891, -0.000000000249999]).T
        self.assertTrue(compare_matlabarray(correctlm, lm, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(info.g, self.info.g, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(info.ae, self.info.ae, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(info.ce, self.info.ce, self.abs_tol, self.rel_tol))

        
        #print "lm", lm
        #print "info", info

if __name__ == '__main__':
    unittest.main()