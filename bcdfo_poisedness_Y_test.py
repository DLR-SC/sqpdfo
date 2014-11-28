# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 13:54:53 2014

@author: jaco_da
"""
import unittest
from bcdfo_poisedness_Y import bcdfo_poisedness_Y_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y
import numpy as np
from runtime import matlabarray
#import helper

class Test_bcdfo_poisedness_Y(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	@unittest.expectedFailure
	def test_bcdfo_poisedness_Y(self):
		print "Self-test of the function with shifting:"
		Y = np.array([ [ 0, 1, 0, 2, 1, 0 ], [ 0, 0, 1, 0, 0.01, 2] ])
		print Y
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y( Y, 0, 1, 1, 1, 1e15 )
		QZ = matlabarray(QZ)
		RZ = matlabarray(RZ)
		xbase = matlabarray(xbase)
		scale = matlabarray(scale)
		
		print "Warning: Not enough input arguments"
		lambd ,Y_radius  = bcdfo_poisedness_Y_( QZ, RZ, Y, 0.001, xbase, 1, 0, scale, 1 )
		print lambd, Y_radius
#  should give:
#  lambd =
#  204.8586
#  Y_radius =
#     2

if __name__ == '__main__':
	unittest.main()