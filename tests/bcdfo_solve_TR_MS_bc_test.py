# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 13:25:23 2014

@author: jaco_da
"""

import unittest
from bcdfo_solve_TR_MS_bc import *
from runtime import matlabarray, compare_matlabarray
#import numpy as np
#import helper

class Test_bcdfo_solve_TR_MS_bc(unittest.TestCase):

    """
      Reminder :
      This class is a test for solve_TR_MS which as its name indicates, solves the Trust Region minimization based on the More-Sorensen algorithm subject to bound  constraints
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
        self.abs_tol=1e-13;
        self.rel_tol=1e-13;
        pass

    def test_bcdfo_solve_TR_MS_bc(self):
        """
            Tests initially written in the matlab code.  We compare the results from python with the results from matlab and verify that it is correct
        """
        gx = matlabarray([[2,3]])
        H = matlabarray([[4, 6], [6, 5]])
        lb = matlabarray([[-10, -10]])
        ub = matlabarray([[10,10]])
        Delta = 1.0
        eps_D = 0.001
        stratLam = 1
        
        s, lamb, norms, value, gplus, nfact, neigd, msg= bcdfo_solve_TR_MS_bc_( gx, H, lb, ub, Delta, eps_D, stratLam )
        correctS = matlabarray( [0.515282741049029, -0.8575287994226390])
        correctgplus=matlabarray( [0.9159581676602855, 4.804052449180984])
        self.assertTrue(compare_matlabarray(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(lamb, 2.103780596518304, places=15)
        self.assertAlmostEqual(norms, 1.000435877536504, places=15)
        self.assertAlmostEqual(value, -1.823817946895660, places=15)
        self.assertTrue(compare_matlabarray(correctgplus, gplus, self.abs_tol, self.rel_tol))
        self.assertEqual(nfact,8)
        self.assertEqual(neigd,0)
        self.assertEqual(str(msg), '(partly) interior solution')
        
        lb = matlabarray([[-0.1, -0.1]])
        s, lamb, norms, value, gplus, nfact, neigd, msg= bcdfo_solve_TR_MS_bc_( gx, H, lb, ub, Delta, eps_D, stratLam )
        correctS = matlabarray( [-0.1, -0.1])
        correctgplus=matlabarray( [0.9159581676602855, 4.804052449180984])
        self.assertTrue(compare_matlabarray(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(lamb, 21.92900284196444, places=15)
        self.assertAlmostEqual(norms, 0.1414213562373095, places=15)
        self.assertAlmostEqual(value, -0.3950, places=15)
        self.assertTrue(compare_matlabarray(correctgplus, gplus, self.abs_tol, self.rel_tol))
        self.assertEqual(nfact,16)
        self.assertEqual(neigd,0)
        self.assertEqual(str(msg), 'boundary solution')
        
        lb = matlabarray([[-10, -10]])
        ub = matlabarray([[0, 0]])
        s, lamb, norms, value, gplus, nfact, neigd, msg= bcdfo_solve_TR_MS_bc_( gx, H, lb, ub, Delta, eps_D, stratLam )
        correctS = matlabarray( [0, -0.6])
        correctgplus=matlabarray( [0.9159581676602855, 4.804052449180984])
        self.assertTrue(compare_matlabarray(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(lamb, 0, places=15)
        self.assertAlmostEqual(norms, 0.6, places=15)
        self.assertAlmostEqual(value, -0.9, places=15)
        self.assertTrue(compare_matlabarray(correctgplus, gplus, self.abs_tol, self.rel_tol))
        self.assertEqual(nfact,9)
        self.assertEqual(neigd,0)
        self.assertEqual(str(msg), '(partly) interior solution')
              

if __name__ == '__main__':
    unittest.main()