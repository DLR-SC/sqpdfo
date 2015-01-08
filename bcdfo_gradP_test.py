# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 11:59:20 2014

@author: jaco_da
"""

import unittest
from bcdfo_gradP import bcdfo_gradP_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray
import numpy as np
#import helper

class Test_bcdfo_gradP(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_bcdfo_gradP(self):
		Y = matlabarray([[ 1, 2, 1, 3, 3, 1],[1, 2, 2, 1, 2, 3 ]])
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 0, 1,1, 1e15 )
		model = ( QZ * np.linalg.solve( RZ.T , matlabarray([1, 2, 3, 4, 5, 6 ]).T ) ).T
		#ans = bcdfo_gradP_( model, matlabarray([0,0]).T, xbase, scale, 0 )
		ans = bcdfo_gradP_( model, matlabarray([[0],[0]]), xbase, scale, 0 )
		print "ans", ans
		correctans = matlabarray([ -6.0000, 1.0000]).T
		
		self.assertTrue((abs(ans - correctans) < 1e-8).all())

if __name__ == '__main__':
	unittest.main()