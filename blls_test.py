# -*- coding: utf-8 -*-
"""
Created on Fri Nov 07 11:52:38 2014

@author: jaco_da
"""

import unittest
from blls import *
import numpy as np
#import helper

class Test_blls(unittest.TestCase):

	def setUp(self):
		self.A = matlabarray([[6,     0,     4],
				[1,     5,     7],
			     [0,     2,     2],
			     [3,     7,     0]])

		self.b = matlabarray([1, 5, 4, 3]).T

#>> blls(A,b, lb, ub)

		self.ans =matlabarray([-0.1713, 0.5162, 0.4663]).T

		self.lb = matlabarray([-1, -1, -1]).T
		self.ub = -self.lb

	def test_blls(self):
		res = blls_(self.A,self.b, self.lb, self.ub)
		
		#print "res", res[0]
		#print "ans", self.ans
		np.testing.assert_almost_equal(res[0], self.ans, 4)
		

if __name__ == '__main__':
	unittest.main()