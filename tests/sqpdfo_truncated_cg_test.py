# -*- coding: utf-8 -*-

import unittest
from sqpdfo.sqpdfo_truncated_cg import *
from numpy import array


class Test_sqpdfo_truncated_cg(unittest.TestCase):
    """
      Reminder :
      This class is a test for sqpdfo_truncated_cg which solves A x = b for x, by Steihaug's conjugate gradient (CG) method,
       where A is a symmetric (possibly indefinite) matrix.
    """

    def setUp(self):
        # self.options = helper.dummyOptions()
        # self.values = helper.dummyValues()

        self.A = array([[1]])
        self.b = array([[0.816496580927725]])
        self.delta = 1
        self.max_iter = 20
        self.tol = 1.000000000000000e-06
        self.plevel = 0
        self.fout = 1

    def test_sqpdfo_truncated_cg(self):
        """
               Test comparing matlab results with python results
        """
        u, info_t = sqpdfo_truncated_cg_(self.A, self.b, self.delta, self.max_iter, self.tol)

        self.assertEqual(u, 0.816496580927725)

        self.assertEqual(info_t.flag, 0)
        self.assertEqual(info_t.iter, 2)
        self.assertEqual(info_t.prec, 0)
        self.assertEqual(info_t.curv, 1)


if __name__ == '__main__':
    unittest.main()