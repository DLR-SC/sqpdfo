# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 14:21:56 2014

INPUT VALUES
nb,mi,me,info,options,values

nb =

     2


mi =

     0


me =

     2


info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 4
      flag: 0
    nsimul: [0 11 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 0.500000000000000
      glag: [3x1 double]
     glagn: 6.661338147750939e-16
     feasn: 2.220446049250313e-16
     compl: 0


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
                         dline: '--------------------------------------------------------------------------------------'
                         eline: '======================================================================================'
                         sline: '**************************************************************************************'

K>> info.g

ans =

  -1.000001597463464
   0.000000267687146
   0.999990706694728

K>> info.ae

ans =

   1.000000000001365   1.000000000001130   0.999999999995923
   0.999999999985469   1.999999999999835   2.999999999990382

K>> info.ce

ans =

   1.0e-15 *

   0.222044604925031
  -0.111022302462516

K>> info.glag

ans =

   1.0e-15 *

   0.666133814775094
   0.661787688722732
  -0.555111512312578
@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from ecdfo_finish import *
#import numpy as np
import helper

class dummyInfo():
    def __init__(self):
        values = helper.dummyValues()
        self.flag = values.stop_on_small_trust_region
        self.f = 0.500000000000000
        self.glagn = 6.661338147750939e-16
        self.feasn = 2.220446049250313e-16
        self.compl = 0
        self.niter=  4
        self.nsimul = matlabarray([0, 11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])

class Test_ecdfo_finish(unittest.TestCase):
    """
      Reminder :
      This class is a test for ecdfo_finish which prints output status
    """ 
    def setUp(self):
        self.options = helper.dummyOptions()
        self.values = helper.dummyValues()

        pass

    def test_ecdfo_finish(self):
        ecdfo_finish_(nb=2,mi=0,me=2,info=dummyInfo(),options=self.options,values=self.values)

if __name__ == '__main__':
    unittest.main()
    