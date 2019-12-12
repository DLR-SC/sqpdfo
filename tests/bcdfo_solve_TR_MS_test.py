# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 16:15:18 2014

@author: jaco_da
"""

import unittest
from sqpdfo.bcdfo_solve_TR_MS import bcdfo_solve_TR_MS_
from sqpdfo.runtime import compare_array
from numpy import array, double


class Test_bcdfo_solve_TR_MS(unittest.TestCase):
    """
      Reminder :
      This class is a test for solve_TR_MS which as its name indicates, solves the Trust Region minimization based on the More-Sorensen algorithm
    """ 
    def setUp(self):
         self.abs_tol=1e-13;
         self.rel_tol=1e-13;
         pass

    def test_bcdfo_solve_TR_MS(self):
        """
            Tests initially written in the matlab code.  We compare the results from python with the results from matlab and verify that it is correct
        """
        s, lamb, norms, value, gplus, nfact, neigd, msg, hardcase = bcdfo_solve_TR_MS_( array([[ 2] ,[ 3] ]), array([[ 4, 6], [6, 5 ]]), 1.0, 0.001 )
        correctS = array( [0.515282741049029, -0.8575287994226390])
        correctgplus=array( [-1.084041832339715, 1.804052449180985])
        self.assertTrue(compare_array(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(lamb, 2.103780596518304, places=15)
        self.assertAlmostEqual(norms, 1.000435877536504, places=14)
        self.assertAlmostEqual(double(value), -1.823817946895660, places=14)
        self.assertTrue(compare_array(correctgplus, gplus, self.abs_tol, self.rel_tol))
        self.assertEqual(nfact,8)
        self.assertEqual(neigd,0)
        self.assertEqual(msg, 'boundary solution')
        self.assertEqual(hardcase,0)
        
        s, lamb, norms, value, gplus, nfact, neigd, msg, hardcase =bcdfo_solve_TR_MS_( array([ [2] , [0] ]), array([[ 4, 0], [0, -15 ]]), 1.0, 0.001 )
        correctS = array( [-0.1052631578947368,0.9944444014574307])
        correctgplus=array( [1.578947368421053, -14.91666602186146])
        self.assertTrue(compare_array(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(lamb, 15, places=15)
        self.assertAlmostEqual(norms, 1, places=15)
        self.assertAlmostEqual(double(value), -7.605263157894736, places=14)
        self.assertTrue(compare_array(correctgplus, gplus, self.abs_tol, self.rel_tol))
        self.assertEqual(nfact,45)
        self.assertEqual(neigd,1)
        self.assertEqual(msg, 'boundary solution ( 45 factorizations, 1 eigen decomposition, lambda = 15.0 )')
        self.assertEqual(hardcase,1)

if __name__ == '__main__':
    unittest.main()