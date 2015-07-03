# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 17:23:18 2014
-----INPUT VALUES-----
func = 

    @evalfgh


x0 =

                   0
                   0
   0.500000000000000


lm0 =

     []


lb =

  -0.500000000000000                   0                -Inf


ub =

   Inf   Inf   Inf


options = 

           algo_method: 'quasi-Newton'
    algo_globalization: 'trust regions'
           hess_approx: 'model'
          bfgs_restart: 0
          algo_descent: 'powell'
                   tol: [1.000000000000000e-05 1.000000000000000e-05 1.000000000000000e-05]
                 dxmin: 1.000000000000000e-06
                 miter: 500
                msimul: 500
               verbose: 2
															
															
-----OUTPUT VALUES-----

K>> x,lm,info

x =

  -0.500000000000000
                   0
   0.500000000000000


lm =

                   0
  -0.000005713064576
                   0
   1.999997749517402
  -0.999996152071198


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

@author: jaco_da
"""

import unittest
from evalfgh import evalfgh_
from runtime import matlabarray, char
from ecdfo import ecdfo_
import numpy as np
import helper

class myOptions():
	def __init__(self):
		self.algo_method = char('quasi-newton')
		self.algo_globalization = char('trust regions')
		self.hess_approx = char('model')
		self.bfgs_restart = 0
		self.algo_descent = char('powell')
		self.tol = matlabarray( [1.000000000000000e-05, 1.000000000000000e-05, 1.000000000000000e-05])
		self.dxmin = 1.000000000000000e-06
		self.miter = 500
		self.msimul = 500
		self.verbose = 2

class Test_ecdfo(unittest.TestCase):

	def setUp(self):

		self.func = evalfgh_
		self.x0 = matlabarray([0, 0, 0.500000000000000]).T
		self.lm0 = matlabarray( [])
		self.lb = matlabarray([ -0.500000000000000, 0, -np.Inf])
		self.ub = matlabarray([np.Inf,   np.Inf,   np.Inf])
		
		#options = 

           #algo_method: 'quasi-Newton'
    #algo_globalization: 'trust regions'
     #      hess_approx: 'model'
     #     bfgs_restart: 0
      #    algo_descent: 'powell'
       #            tol: [1.000000000000000e-05 1.000000000000000e-05 1.000000000000000e-05]
        #         dxmin: 1.000000000000000e-06
         #        miter: 500
          #      msimul: 500
           #    verbose: 2		
		
		self.options = myOptions()#helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_ecdfo(self):
		x,lm,info = ecdfo_(self.func,self.x0,self.lm0,self.lb,self.ub,self.options)
		
		print "x\n", x
		print "lm\n", lm
		
		

if __name__ == '__main__':
	unittest.main()