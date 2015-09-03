# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 13:25:23 2014

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from bcdfo_solve_TR_MS_bc import *
from runtime import compare_array
from numpy import array
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
        self.abs_tol=1e-14;
        self.rel_tol=1e-14;
        pass

    def test_bcdfo_solve_TR_MS_bc_1(self):
        """
            Tests initially written in the matlab code.  We compare the results from python with the results from matlab and verify that it is correct
        """
        gx = array([[2.],[3.]])
        H = array([[4., 6.], [6., 5.]])
        lb = array([[-10.], [-10.]])
        ub = array([[10.],[10.]])
        Delta = 1.0
        eps_D = 0.001
        stratLam = 1
        
        s, lamb, norms, value, gplus, nfact, neigd, msg= bcdfo_solve_TR_MS_bc_( gx, H, lb, ub, Delta, eps_D, stratLam )
        correctS = array( [0.515282741049029, -0.8575287994226390])
        correctgplus=array( [0.9159581676602855, 4.804052449180984])
        self.assertTrue(compare_array(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(lamb, 2.103780596518304, places=15)
        self.assertAlmostEqual(norms, 1.000435877536504, places=14)
        self.assertAlmostEqual(value, -1.823817946895660, places=14)
        self.assertTrue(compare_array(correctgplus, gplus, self.abs_tol, self.rel_tol))
        self.assertEqual(nfact,8)
        self.assertEqual(neigd,0)
        self.assertEqual(str(msg), '(partly) interior solution')
        
        lb = array([[-0.1], [-0.1]])
        s, lamb, norms, value, gplus, nfact, neigd, msg= bcdfo_solve_TR_MS_bc_( gx, H, lb, ub, Delta, eps_D, stratLam )
        correctS = array( [-0.1, -0.1])
        correctgplus=array( [0.9159581676602855, 4.804052449180984])
        self.assertTrue(compare_array(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(lamb, 21.92900284196444, places=15)
        self.assertAlmostEqual(norms, 0.1414213562373095, places=15)
        self.assertAlmostEqual(value, -0.3950, places=15)
        self.assertTrue(compare_array(correctgplus, gplus, self.abs_tol, self.rel_tol))
        self.assertEqual(nfact,16)
        self.assertEqual(neigd,0)
        self.assertEqual(str(msg), 'boundary solution')
        
        lb = array([[-10], [-10]])
        ub = array([[0], [0]])
        s, lamb, norms, value, gplus, nfact, neigd, msg= bcdfo_solve_TR_MS_bc_( gx, H, lb, ub, Delta, eps_D, stratLam )
        correctS = array( [0, -0.6])
        correctgplus=array( [0.9159581676602855, 4.804052449180984])
        self.assertTrue(compare_array(correctS, s, self.abs_tol, self.rel_tol))
        self.assertAlmostEqual(lamb, 0, places=15)
        self.assertAlmostEqual(norms, 0.6, places=15)
        self.assertAlmostEqual(value, -0.9, places=15)
        self.assertTrue(compare_array(correctgplus, gplus, self.abs_tol, self.rel_tol))
        self.assertEqual(nfact,9)
        self.assertEqual(neigd,0)
        self.assertEqual(str(msg), '(partly) interior solution')
#        
    def test_bcdfo_solve_TR_MS_bc_2(self):
        """
            Test to debug bcdfo_solve_TR_MS_bc because it does not yield the same result as matlab for theses values. This error has been found while testing ecdfo_solve_TR_MS_bc
        """
        g1=array([[ 5.], [ 8.], [11.]])
        H1=array([[ 2., 3, 4.], [ 3., 5., 7.], [ 4., 7., 10.]])
        lb_r=array([[ -1.], [ -1.], [-inf]])
        ub_r=array([[ inf], [ inf], [ inf]])   
        delta_r=1.0
        prec_r=1e-6
        stratLam=1
        s, lamb, norms, value, gplus, nfact, neigd, msg=bcdfo_solve_TR_MS_bc_(g1,H1,lb_r,ub_r,delta_r,prec_r,stratLam,nargout=8)
        correctS = array([-7.217986001584447e-01,-5.625108528130431e-01,-4.032231054676398e-01])
        correctgplus=array( [5.255977819373422e+00, 8.199488197185971,1.114299857499852e+01])
        
        #The relatively low accuracy of the results, compared to matlab results, comes from a different result in the cholesky
#        factorization during the call to bcdfo_solve_TR_MS. The last number of the matrix has its eigth digit different in python and matlab,
#        therefore we loose so much in accuracy. But other than that, I believe that the test is correct.
        self.assertTrue(compare_array(correctS, s, 1e-6, 1e-6))
        self.assertAlmostEqual(lamb, 3.546388415234277e-01, places=5)
        self.assertAlmostEqual(norms, 1.000000275753019, places=6)
        self.assertAlmostEqual(value, -6.449586510274759, places=6)
        self.assertTrue(compare_array(correctgplus, gplus, 1e-6, 1e-6))
        self.assertEqual(nfact,7)
        self.assertEqual(neigd,0)
        self.assertEqual(str(msg), '(partly) interior solution')

if __name__ == '__main__':
    unittest.main()