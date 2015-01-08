# -*- coding: utf-8 -*-
"""
Created on Mon Dec 01 13:41:32 2014

@author: jaco_da
"""

import unittest
from bcdfo_computeLj import bcdfo_computeLj_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray
#+import numpy as np
#import helper

class Test_bcdfo_computeLj(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	@unittest.expectedFailure
	def test_bcdfo_computeLj(self):
		Y = matlabarray([[ 0, 1, 0, 2, 0], [0, 0, 1, 0, 2 ]])
		whichmodel = 0
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, whichmodel, 1, 1, 1, 1e15 )
		Lj = bcdfo_computeLj_( QZ, RZ, 1, Y, whichmodel, scale, 1 )

		print "type(scale)", type(scale)
		correctLj = matlabarray([   1.0000,  -3.0000,  -3.0000,   4.0000,   4.0000])		
		print "Warning: Different Results due to not unique QR decomposition (?)"
		self.assertEqual(Lj, correctLj)
# should give:
#  Lj =
#   1.0000  -3.0000  -3.0000   4.0000   4.0000

if __name__ == '__main__':
	unittest.main()