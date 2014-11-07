# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 11:53:25 2014

@author: jaco_da
"""

import unittest
from ecdfo_check_cond import *
import numpy as np
import helper

class Test_ecdfo_check_cond(unittest.TestCase):

	def setUp(self):
		self.A = np.array([[ 1,  0, -1], [ 0,  1,  0], [ 1,  0,  1]])
		self.dummyOptions = helper.dummyOptions()
		print self.A

	def test_ecdfo_check_cond1(self):
		res = ecdfo_check_cond_(self.A, 1.41421356237, self.dummyOptions)
		self.assertTrue(res[1])
		print "np cond", np.linalg.cond(self.A)
		print "ecdfo  check cond", res

	def test_ecdfo_check_cond2(self):
		res = ecdfo_check_cond_(self.A, 1.41421356238, self.dummyOptions)
		self.assertFalse(res[1])
		#self.assertFalse(True)
		print "ecdfo  check cond", res


if __name__ == '__main__':
	unittest.main()