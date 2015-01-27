# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 15:04:23 2014

@author: jaco_da
"""

import unittest
from bcdfo_evalP import bcdfo_evalP_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray
import numpy as np
#import helper

class Test_bcdfo_evalP(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_bcdfo_evalP(self):
#		TEST:
		Y = matlabarray([[ 1, 2, 1, 3, 3, 1],[1, 2, 2, 1, 1.01, 3 ]]) 
		QZ, RZ, xbase, scale  = bcdfo_build_QR_of_Y_( Y, 0, 1 , 1,1, 1e15)
		#print "scale", scale
		#print "scale[1]", scale[1][1]
		model = ( QZ * np.linalg.solve( RZ.T , matlabarray([[1], [2], [3], [4], [5], [6] ]) ) ).T
		ans = bcdfo_evalP_( model, matlabarray([[1],[3]]), xbase, scale, 1 )
		#print "ans", ans
		self.assertAlmostEqual(float(ans), 6.0)
#  should give 
#     6.0

if __name__ == '__main__':
	unittest.main()