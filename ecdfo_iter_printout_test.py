# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 15:24:01 2014

INPUT VALUES
info,old_delta,norms,pc,itype,values,nb,mi,options,constrained_pbl,merit

info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 1
      flag: 0
    nsimul: [0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 1.500000000000000
      glag: [3x1 double]
     glagn: 0.666666666736111
     feasn: 3
     compl: 0


old_delta =

     1


norms =

     0


pc =

     0


itype =

 


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


nb =

     2


mi =

     0


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


constrained_pbl =

     2


merit =

   5.105551275463990

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

K>> info.glag

ans =

  -0.333333333013890
   0.666666666736111
  -0.333333333513889
		
OUTPUT VALUES
none
@author: jaco_da
"""

import unittest
from ecdfo_iter_printout import *
#import numpy as np
import helper

class dummyInfo():
	def __init__(self):
		self.g = matlabarray([0,1,0]).T
		self.ai = matlabarray([])
		self.ae = ([[1,     1,     1],[ 1,     2,     3]])
		self.hl = matlabarray([])
		self.niter = 1
		self.flag = 0
		self.nsimul = matlabarray( [0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
		self.ci = matlabarray([])
		self.ce = matlabarray([2, 3]).T
		self.f = 1.500000000000000
		self.glag = matlabarray([  -0.333333333013890,  0.666666666736111,  -0.333333333513889]).T
		self.glagn = 0.666666666736111
		self.feasn = 3
		self.compl = 0
		
class Test_ecdfo_iter_printout(unittest.TestCase):

	def setUp(self):
		self.info = dummyInfo()
		self.old_delta =  1
		self.norms =    0
		self.pc =    0
		self.itype = None
		self.options = helper.dummyOptions()
		self.nb = 2
		self.mi =  0
		self.values = helper.dummyValues()
		self.constrained_pbl =  2
		self.merit =   5.105551275463990

	def test_ecdfo_iter_printout(self):
		ecdfo_iter_printout_(self.info,self.old_delta,self.norms,self.pc,
		self.itype,self.values,self.nb, self.mi,self.options,self.constrained_pbl,self.merit)

if __name__ == '__main__':
	unittest.main()