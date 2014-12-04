# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 12:24:11 2014

@author: jaco_da
"""

import unittest
from runtime import matlabarray
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from bcdfo_hessP import bcdfo_hessP_
import numpy as np
#import helper

class Test_bcdfo_hessP(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_bcdfo_hessP(self):
		Y = matlabarray([[ 1, 2, 1, 3, 3, 1],[ 1, 2, 2, 1, 2, 3 ]])
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 0, 1, 1, 1e15 );
		model = ( QZ * np.linalg.solve( RZ.T , matlabarray([1, 2, 3, 4, 5, 6 ]).T ) ).T
		ans = bcdfo_hessP_( model, matlabarray([[0],[0]]) ,xbase, scale, 0 )
		
		correctans = matlabarray([[    4.0000,   -0.5000],
						[-0.5000,    1.0000]])
		
		self.assertTrue((abs(ans - correctans) < 1e-8).all())
		print "ans", ans
  
#  ans =
#
#    4.0000   -0.5000
#   -0.5000    1.0000

if __name__ == '__main__':
	unittest.main()