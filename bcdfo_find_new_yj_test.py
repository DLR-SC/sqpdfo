# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 11:17:28 2014

@author: jaco_da
"""

import unittest
from bcdfo_find_new_yj import bcdfo_find_new_yj_, bcdfo_find_new_yj
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_, bcdfo_build_QR_of_Y
from runtime import matlabarray
import numpy as np
#import helper

class Test_bcdfo_find_new_yj(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def est_bcdfo_find_new_yj(self):
		Y = matlabarray([[ 3, 1, 0, 2, 1, 0],[0, 0, 1, 0, 0.01, 2 ]])
		whichmodel = 0
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y , whichmodel, 0, 1,1, 1e15 )
		ynew, improvement = bcdfo_find_new_yj_( QZ, RZ, Y, 5, 1.0, 0.001, xbase, 1, whichmodel, scale, 0 )
		
		print "ynew", ynew
		print "improvement", improvement
	
	@unittest.expectedFailure	
	def test_bcdfo_find_new_yjOriginal(self):
		Y = np.array([[ 3, 1, 0, 2, 1, 0],[0, 0, 1, 0, 0.01, 2 ]])
		whichmodel = 0
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y( Y , whichmodel, 0, 1,1, 1e15 )
		ynew, improvement = bcdfo_find_new_yj( QZ, RZ, Y, 5, 1.0, 0.001, xbase, 1, whichmodel, scale, 1 )
		

		print "QZ\n", QZ
		print "RZ\n", RZ		
		
		correctynew = matlabarray([[3],[0]])
		print "ynew", ynew
		self.assertEqual(ynew, correctynew)
		print "improvement", improvement
		correctimprovement = 0
		self.assertEqual(improvement, correctimprovement)
		
#  should give
#  ynew =
#
#     3
#     0
#
#
#improvement =
#
#     0


if __name__ == '__main__':
	unittest.main()