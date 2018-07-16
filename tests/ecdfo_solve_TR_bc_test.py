# -*- coding: utf-8 -*-
"""
Created on Wed Nov 19 10:22:26 2014

@author: jaco_da
"""

import sys
sys.path.append("../")
import unittest
from ecdfo_solve_TR_bc import *
from evalfgh import evalfgh_
import helper
from ecdfo_global_variables import set_threshold
from numpy import array


class dummyInfo():
    def __init__(self):
        self.g =  array([[0, 1, 0]]).T
        self.ai = array([])
        self.ae = array([[1,     1,     1], [1,     2,     3]])
        self.hl = array([])
        self.niter = 1
        self.flag = 0
        self.nsimul = array( [0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.ci=  array([])
        self.ce = array([[2,3]]).T
        self.f = 1.500000000000000
        self.glag = array([[-0.333333333013890, 0.666666666736111, -0.333333333513889]]).T
        self.glagn = 0.666666666736111
        self.feasn = 3
        self.compl = 0


class Test_ecdfo_solve_TR_bc(unittest.TestCase):
    """
      Reminder :
      This class is a test for ecdfo_solve_TR_bc which compute the new SQP-trust-region step inside delta and subject to simple bounds
      accordingly.
    """ 
    def setUp(self):

        self.simul = evalfgh_
        self.x = array([[0.500000000000000, 1.000000000000000, 0.500000000000000]]).T
        self.lb = array([[ -0.500000000000000, 0, -np.Inf]]).T
        self.ub = array([[np.Inf, np.Inf, np.Inf]]).T
        self.delta =  1
        self.mi =  0
        self.me =  2
        self.M = array([[1,     0,     0], [0,     1,     0], [0,     0,     1]])
        self.prec_r = 1.000000000000000e-06
        self.prec_t = 1.000000000000000e-06
        
        self.info = dummyInfo()
        self.options = helper.dummyOptions()
        self.values = helper.dummyValues()
    
        self.radius_has_been_rejected = 0
        self.lm = array([[0, 0, 0, -0.333333332763891, -0.000000000249999]]).T
        self.ceY = array([[ 2,     1,     1,     1],[ 3,     2,     1,     0]])
        self.ciY = array([])
        self.gx =array([[ 0, 1, 0]]).T
        
        self.abs_tol=1e-9
        self.rel_tol=1e-9

    def test_ecdfo_solve_TR_bc(self):
        """
            This is a test with some values and compared with matlab results.
        """
        set_check_condition(0)                
        set_threshold(1.000000000000000e-08)
        xnew,delta,rpred,active_r,active_t,lm_computed,lm,info = ecdfo_solve_TR_bc_(self.simul,self.x,self.lb,self.ub,self.delta,
        self.mi,self.me,self.M,self.prec_r,self.prec_t,self.info,
        self.options,self.values,self.radius_has_been_rejected,self.lm,
        self.ceY,self.ciY,self.gx)
        
        
        correctxnew =array([[  -3.054026564966283e-03, 0, 3.155214681258386e-01]]).T
        correctdelta = 0.43274889405056605
        correctrpred =  3.288018632338156e+00
        correctactive_r =  0
        correctactive_t =  0
        correctlm_computed = 0
        correctlm = array([[0, -1, 0, 0, 0]]).T
        
        self.assertTrue(compare_array(correctxnew, xnew, 1e-6, 1e-3))
        self.assertAlmostEqual(correctdelta, delta, places=6)
        self.assertAlmostEqual(correctrpred, rpred, places=6)
        self.assertEqual(correctactive_r, active_r)
        self.assertEqual(correctactive_t, active_t)
        self.assertEqual(correctlm_computed, lm_computed)
        self.assertTrue(compare_array(correctlm,lm, self.abs_tol, self.rel_tol))

if __name__ == '__main__':
#<<<<<<< HEAD
    unittest.main()
#
#=======
#    unittest.main()
#>>>>>>> a959f16419e1ef0052e0d04109e5c0b7e711e3f0
