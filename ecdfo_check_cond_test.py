# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 11:53:25 2014

@author: jaco_da
"""

import unittest
from ecdfo_check_cond import *
import numpy as np

class Test_ecdfo_check_cond(unittest.TestCase):

	def setUp(self):
		self.A = np.array([[ 1,  0, -1], [ 0,  1,  0], [ 1,  0,  1]])
		print self.A

	def test_ecdfo_check_cond1(self):
		res = ecdfo_check_cond_(self.A, 1.41421356237)
		self.assertTrue(res[1])
		print "np cond", np.linalg.cond(A)
		print "ecdfo  check cond", res

	def test_ecdfo_check_cond2(self):
		res = ecdfo_check_cond_(self.A, 1.41421356238)
		self.assertFalse(res[1])
		print "ecdfo  check cond", res


if __name__ == '__main__':
	unittest.main()