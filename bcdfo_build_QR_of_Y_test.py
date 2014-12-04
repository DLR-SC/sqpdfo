

# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 11:58:39 2014

@author: jaco_da
"""

import unittest
from bcdfo_build_QR_of_Y import *
import numpy as np
from runtime import *
import helper

class Test_bcdfo_build_QR_of_Y(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_bcdfo_build_QR_of_Y(self):
		Y = matlabarray([[ 1, 2, 1, 3, 3, 1],  [1, 2, 2, 1, 1.01, 3 ]])
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
		print "QZ", QZ
		print "RZ", RZ
		model = ( QZ * np.linalg.solve( RZ.T , matlabarray([1, 2, 3, 4, 5, 6 ]).T ) ).T
		print "model", model
		res = np.dot(helper.convert(model) , bcdfo_evalZ( (np.array([[1,3]]).T-helper.convert(xbase))*helper.convert(scale[2]),6))
		print "res", res
#  should give 6.0
		self.assertAlmostEqual(float(res), 6)

if __name__ == '__main__':
	unittest.main()

