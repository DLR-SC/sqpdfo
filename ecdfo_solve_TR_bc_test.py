# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 10:22:26 2014
INPUT VALUES

 simul,x,lb,ub,delta,mi,me,M,prec_r,prec_t,info,options,values,...
    radius_has_been_rejected,lm,ceY,ciY,gx

simul = 

    @evalfgh


x =

   0.500000000000000
   1.000000000000000
   0.500000000000000


lb =

  -0.500000000000000
                   0
                -Inf


ub =

   Inf
   Inf
   Inf


delta =

     1


mi =

     0


me =

     2


M =

     1     0     0
     0     1     0
     0     0     1


prec_r =

     1.000000000000000e-06


prec_t =

     1.000000000000000e-06


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


radius_has_been_rejected =

     0


lm =

                   0
                   0
                   0
  -0.333333332763891
  -0.000000000249999


ceY =

     2     1     1     1
     3     2     1     0


ciY =

     []


gx =

     0
     1
     0

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

K>> xnew,deltaTR,rpred,active_r,active_t,lm_computed,lm,info

xnew =

  -0.003054026564966
                   0
   0.315521468125839


deltaTR =

  -0.134502394164804


rpred =

   3.288018632338156


active_r =

     0


active_t =

     0


lm_computed =

     0


lm =

                   0
  -0.999999999499998
                   0
  -0.000000000375002
   0.000000000083334


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


@author: jaco_da
"""

import unittest
from ecdfo_solve_TR_bc import *
from evalfgh import evalfgh_
#import numpy as np
import helper

class dummyInfo():
	def __init__(self):
		self.g =  matlabarray([0, 1, 0]).T
		self.ai = matlabarray([])
		self.ae = matlabarray([[1,     1,     1], [1,     2,     3]])
		self.hl = matlabarray([])
		self.niter = 1
		self.flag = 0
		self.nsimul = matlabarray( [0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
		self.ci=  matlabarray([])
		self.ce = matlabarray([2,3]).T
		self.f = 1.500000000000000
		self.glag = matlabarray([-0.333333333013890, 0.666666666736111, -0.333333333513889]).T
		self.glagn = 0.666666666736111
		self.feasn = 3
		self.compl = 0
		
class Test_ecdfo_solve_TR_bc(unittest.TestCase):

	def setUp(self):
		self.simul = evalfgh_
		self.x = matlabarray([0.500000000000000, 1.000000000000000, 0.500000000000000]).T
		self.lb = matlabarray([ -0.500000000000000, 0, -np.Inf]).T
		self.ub = matlabarray([np.Inf, np.Inf, np.Inf]).T
		self.delta =  1
		self.mi =  0
		self.me =  2
		self.M = matlabarray([[1,     0,     0], [0,     1,     0], [0,     0,     1]])
		self.prec_r = 1.000000000000000e-06
		self.prec_t = 1.000000000000000e-06
		
		self.info = dummyInfo()
		self.options = helper.dummyOptions()
		self.values = helper.dummyValues()
	
		self.radius_has_been_rejected = 0
		self.lm = matlabarray([0, 0, 0, -0.333333332763891, -0.000000000249999]).T
		self.ceY = matlabarray([[ 2,     1,     1,     1],[ 3,     2,     1,     0]])
		self.ciY = matlabarray([])
		self.gx =matlabarray([ 0, 1, 0]).T

	@unittest.expectedFailure
	def test_ecdfo_solve_TR_bc(self):
		set_Threshold(1.000000000000000e-08)
		xnew,delta,rpred,active_r,active_t,lm_computed,lm,info = ecdfo_solve_TR_bc_(self.simul,self.x,self.lb,self.ub,self.delta,
		self.mi,self.me,self.M,self.prec_r,self.prec_t,self.info,
		self.options,self.values,self.radius_has_been_rejected,self.lm,
		self.ceY,self.ciY,self.gx)
		
		
		correctxnew =matlabarray([  -0.003054026564966, 0, 0.315521468125839]).T
		correctdelta = -0.134502394164804
		correctrpred =  3.288018632338156
		correctactive_r =  0
		correctactive_t =  0
		correctlm_computed = 0
		correctlm = matlabarray([0, -0.999999999499998, 0, -0.000000000375002, 0.000000000083334]).T
		
		print "THRU"
		print xnew,delta,rpred,active_r,active_t,lm_computed,lm,info
		print "Warning: Correctness of Code is questionable"
		self.assertEqual(correctxnew, xnew)
		self.assertEqual(correctdelta, delta)
		self.assertEqual(correctrpred, rpred)
		self.assertEqual(correctactive_r, active_r)
		self.assertEqual(correctactive_t, active_t)
		self.assertEqual(correctlm_computed, lm_computed)
		self.assertEqual(correctlm, lm)

if __name__ == '__main__':
	unittest.main()