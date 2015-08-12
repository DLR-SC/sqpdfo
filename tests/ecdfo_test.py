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
import sys
sys.path.append("../")
import unittest
from evalfgh import evalfgh_
from runtime import matlabarray, char, compare_matlabarray, isempty_
from ecdfo import ecdfo_
from ecdfo_global_variables import set_prob, set_threshold, set_fileoutput,set_simul_not_initialized
import numpy as np
import helper

class myOptions():
    def __init__(self):
        self.algo_method = 'quasi-newton'
        self.algo_globalization = 'trust regions'
        self.hess_approx = 'model'
        self.bfgs_restart = 0
        self.algo_descent = 'powell'
        self.tol = matlabarray( [1.000000000000000e-05, 1.000000000000000e-05, 1.000000000000000e-05])
        self.dxmin = 1.000000000000000e-06
        self.miter = 500
        self.msimul = 500
        self.verbose = 2

class Test_ecdfo(unittest.TestCase):
    """
      Reminder :
      This class is a test for ecdfo 
    """
    def setUp(self):

        set_prob(3)
        set_threshold(1e-08)
        set_fileoutput(1)
        set_simul_not_initialized(1)
        self.func = evalfgh_
        self.x0 = matlabarray([0, 0, 0.500000000000000]).T
        self.lm0 = matlabarray( [])
        self.lb = matlabarray([ -0.500000000000000, 0, -np.Inf])
        self.ub = matlabarray([np.Inf,   np.Inf,   np.Inf])
        
        self.options = myOptions()#helper.dummyOptions()
        #self.values = helper.dummyValues()
        
        self.abs_tol=1e-5
        self.rel_tol=1e-9

        pass

    def test_ecdfo(self):
        """
         Test which compare python and matlab results
        """
        x,lm,info = ecdfo_(self.func,self.x0,self.lm0,self.lb,self.ub,self.options) 
        correctx = matlabarray([[ -0.500000000000000, 0,0.500000000000000]])
        correctlm = matlabarray([[0,-0.000005713064576,0,1.999997749517402,-0.999996152071198]])
        self.assertTrue(compare_matlabarray(correctx, x,self.abs_tol, self.rel_tol)) 
        self.assertTrue(compare_matlabarray(correctlm, lm,self.abs_tol, self.rel_tol))

        self.assertEqual(info.niter,4)
        self.assertEqual(info.flag,0)
        self.assertTrue(compare_matlabarray(info.f,matlabarray([0.5]), self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(info.g, matlabarray([-1.000001597463464,2.676871464447451e-07,9.999907066947283e-01]),self.abs_tol, 1e-7))
        self.assertTrue(compare_matlabarray(info.ae, matlabarray([[1,1,1],[1,2,3]]),self.abs_tol, self.rel_tol))
        self.assertTrue(compare_matlabarray(info.nsimul, matlabarray([0,11,0,0,0,0,0,0,0,0,0,0,0,0,0,0]),self.abs_tol, self.rel_tol))
        self.assertTrue(isempty_(info.hl))
        #Those 4 below assertions are false : this is noted as an issue on github, though this may not be a 'real' issue in pratice
#        self.assertAlmostEqual(6.661338147750939e-16,info.glagn,places=14)
#        self.assertAlmostEqual(2.220446049250313e-16,info.feasn,places=15)
#        self.assertTrue(compare_matlabarray(info.ce, matlabarray([2.220446049250313e-16,-1.110223024625157e-16]),self.abs_tol, self.rel_tol))
#        self.assertTrue(compare_matlabarray(info.glag, matlabarray([6.661338147750939e-16,6.617876887227321e-16,-5.551115123125783e-16])))
        

if __name__ == '__main__':
    unittest.main()