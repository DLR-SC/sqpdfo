# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 16:15:18 2014

@author: jaco_da
"""

import unittest
from bcdfo_solve_TR_MS import bcdfo_solve_TR_MS_
from runtime import matlabarray
#import numpy as np
#import helper

class Test_bcdfo_solve_TR_MS(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_bcdfo_solve_TR_MS(self):
		#  TEST:

		ans = bcdfo_solve_TR_MS_( matlabarray([ 2 , 3 ]), matlabarray([[ 4, 6], [6, 5 ]]), 1.0, 0.001 )
		print ans
		correctans = matlabarray( [0.5153, -0.8575])
		self.assertTrue((abs(ans[0] - correctans) < 1e-4).all())
#  should give
#    0.5153
#   -0.8575
#		bcdfo_solve_TR_MS( [ 2 ; 0 ], [ 4 0; 0 -15 ], 1.0, 0.001 )
#  should give
#   -0.1053
#    0.9944

if __name__ == '__main__':
	unittest.main()