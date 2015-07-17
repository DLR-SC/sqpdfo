# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 16:15:18 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from bcdfo_solve_TR_MS import bcdfo_solve_TR_MS_
#from bcdfo_solve_TR_MS_Anke import bcdfo_solve_TR_MS_

from runtime import matlabarray, compare_matlabarray
#import numpy as np
#import helper

class Test_bcdfo_solve_TR_MS(unittest.TestCase):
    """
      Reminder :
      This class is a test for solve_TR_MS which as its name indicates, solves the Trust Region minimization based on the More-Sorensen algorithm
    """ 
    def setUp(self):
        #self.options = helper.dummyOptions()
        #self.values = helper.dummyValues()
         self.abs_tol=1e-13;
         self.rel_tol=1e-13;
         pass

    def test_bcdfo_solve_TR_MS(self):
        """
            Tests initially written in the matlab code.  We compare the results from python with the results from matlab and verify that it is correct
        """
        s, lamb, norms, value, gplus, nfact, neigd, msg, hardcase = bcdfo_solve_TR_MS_( matlabarray([[ 2] ,[ 3] ]), matlabarray([[ 4, 6], [6, 5 ]]), 1.0, 0.001 )
        correctS = matlabarray( [0.515282741049029, -0.8575287994226390])
        correctgplus=matlabarray( [-1.084041832339715, 1.804052449180985])
        self.assertTrue(compare_matlabarray(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(lamb, 2.103780596518304, places=15)
        self.assertAlmostEqual(norms, 1.000435877536504, places=15)
        self.assertAlmostEqual(value, -1.823817946895660, places=15)
        self.assertTrue(compare_matlabarray(correctgplus, gplus, self.abs_tol, self.rel_tol))
        self.assertEqual(nfact,8)
        self.assertEqual(neigd,0)
        self.assertEqual(str(msg), 'boundary solution')
        self.assertEqual(hardcase,0)
        
        s, lamb, norms, value, gplus, nfact, neigd, msg, hardcase =bcdfo_solve_TR_MS_( matlabarray([ [2] , [0] ]), matlabarray([[ 4, 0], [0, -15 ]]), 1.0, 0.001 )
        correctS = matlabarray( [-0.1052631578947368,0.9944444014574307])
        correctgplus=matlabarray( [1.578947368421053, -14.91666602186146])
        self.assertTrue(compare_matlabarray(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(lamb, 15, places=15)
        self.assertAlmostEqual(norms, 1, places=15)
        self.assertAlmostEqual(value, -7.605263157894736, places=15)
        self.assertTrue(compare_matlabarray(correctgplus, gplus, self.abs_tol, self.rel_tol))
        self.assertEqual(nfact,45)
        self.assertEqual(neigd,1)
        self.assertEqual(str(msg), "['boundary solution ( 45 factorizations, 1 eigen decomposition, lambda = 15.0 )']")
#        self.assertEqual(str(msg), "[[ 'boundary solution ( 45 factorizations, 1 eigen decomposition, lambda = 15.0 )']]") # for Anke's translation
        self.assertEqual(hardcase,1)

if __name__ == '__main__':
    unittest.main()