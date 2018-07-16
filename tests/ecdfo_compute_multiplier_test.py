# -*- coding: utf-8 -*-
"""
@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from ecdfo_compute_multiplier import *
from ecdfo_global_variables import set_check_condition
import helper
from numpy import array


class dummyInfo():
    def __init__(self):
        #g: [3x1 double]
        self.g = array([[0, 1, 0]]).T
        self.ai = array([])
        #ae : [2x3 double]
        self.ae = array([[1,     1,     1],    [ 1 ,    2,     3]])
        self.hl = array([])
        self.niter = 0
        self.flag = 0
        self.nsimul = array([0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
        self.ci = array([])
        #ce: [2x1 double]
        self.ce = array([[2, 3]]).T
        self.f = 1.5

class Test_ecdfo_compute_multiplier(unittest.TestCase):
    """
      Reminder :
      This class is a test for ecdfo_compute_multiplier which computes exact least-squares multipliers 'lm'.
      It solves in lm the quadratic optimization problem:
           min || g+A'*lm ||^2
      subject to possible bounds on lm.
    """
    def setUp(self):
        set_check_condition(0)
        self.x = array([[0.5, 1.0, 0.5]]).T
        self.lb = array([])
        self.ub = array([])
        self.values = helper.dummyValues()
        self.options = helper.dummyOptions()
        self.info = dummyInfo()
        self.abs_tol=1e-8;
        self.rel_tol=1e-8;

    def test_ecdfo_compute_multiplier(self):
        """
            Test comparing python results with the matlab results
        """
        lm,info = ecdfo_compute_multiplier_(self.x,self.lb,self.ub,self.info,self.options,self.values)

        correctlm = array([0, 0, 0, -0.333333332763891, -0.000000000249999]).T

        self.assertTrue(compare_array(correctlm, lm, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.g, self.info.g, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.ae, self.info.ae, self.abs_tol, self.rel_tol))
        self.assertTrue(compare_array(info.ce, self.info.ce, self.abs_tol, self.rel_tol))


if __name__ == '__main__':
    unittest.main()
