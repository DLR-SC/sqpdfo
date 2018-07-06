# -*- coding: utf-8 -*-

import sys
sys.path.append("../")
import unittest
from ecdfo_options import *
from numpy import array
import helper

class dummyInfo():
    def __init__(self):
        self.g= array([])
        self.ai= array([])
        self. ae= array([])
        self.hl= array([])
        self.niter= 0


class Test_ecdfo_options(unittest.TestCase):
    """
      Reminder :
      This class is a test for ecdfo_options which sets the options of the optimization solver 'ecdfo'
    """

    def setUp(self):
        self.options = helper.dummyUnionStruct()
        self.options.dxmin = 1e-6
        self.options.algo_descent = 'powell'
        self.options.hess_approx = 'bfgs'
        self.options.fout = -1

        self.info = dummyInfo()

    def test_sqplab_options1(self):
        """
        First little test
        """
        info, options, values = ecdfo_options_(self.info, self.options)
        self.assertEqual(options.dxmin, 1e-6)
        self.assertEqual(options.fout, 1)

    def test_sqplab_options2(self):
        """
        Second little test
        """
        info, options, values = ecdfo_options_(self.info, self.options)
        self.assertEqual(options.algo_descent, values.powell)
        self.assertEqual(options.hess_approx, values.bfgs)


if __name__ == '__main__':
    unittest.main()