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

	def test_bcdfo_find_new_yj(self):
		#changed example to reproduce results from MATLAB, due to optimization
		Y = matlabarray([[ 3.0, 1.0, 1.0, 2.0, 1.0, 0.0],[0.0, 0.0, 1.0, 0.0, 0.01, 2.0 ]])
		whichmodel = 0
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y , whichmodel, 1, 1,1, 1e15 )
		ynew, improvement = bcdfo_find_new_yj_( QZ, RZ, Y, 5, 1.0, 0.001, xbase, 1, whichmodel, scale, 1 )


		#QZ = matlabarray(
    #[[1.0000,         0,         0,         0,         0,         0],
    #[     0,   -0.9636,   -0.0783,   -0.2555,         0,   -0.0000],
    # [    0,         0,   -0.7308,    0.2240,    0.5993,    0.2378],
    #  [   0,    0.2673,   -0.2823,   -0.9213,         0,   -0.0000],
    #   [  0,         0,   -0.1013,    0.0311,   -0.4806,    0.8705],
    #    [ 0,         0,    0.6081,   -0.1863,    0.6402 ,   0.4309]])

		#RZ = matlabarray(

   #[ [1.0000,    1.0000,    1.0000,    1.0000,    1.0000,    1.0000],
   #  [    0,    0.5756,    0.8943,    0.2775,    0.5756,    0.8943],
   #   [   0,         0  , -0.3795,    0.0109,   -0.0030,   -0.7342],
   #    [  0,         0 ,        0,    0.0354,    0.0009,    0.1087],
   #     [ 0,         0,         0,         0,    0.0007,   -0.0370],
   #     [ 0,        0,         0 ,        0,         0 ,   0.0670]])

		correctynew = matlabarray([3.2866, 0.9581]).T
		correctimprovement = 217.2211
		#print "ynew", ynew
		#print "improvement", improvement
		
		self.assertAlmostEqual(correctimprovement, improvement, 4)
		#print "abs", abs(ynew - correctynew)
		self.assertTrue((abs(correctynew - ynew) < 1e-3).all())

	
	#@unittest.expectedFailure	
	def est_bcdfo_find_new_yjOriginal(self):
		Y = np.array([[ 3.0, 1.0, 0.0, 2.0, 1.0, 0.0],[0.0, 0.0, 1.0, 0.0, 0.01, 2.0 ]])
		whichmodel = 0
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y( Y , whichmodel, 1, 1,1, 1e15 )
		ynew, improvement = bcdfo_find_new_yj( QZ, RZ, Y, 5-1, 1.0, 0.001, xbase, 1, whichmodel, scale, 1 )
		

		#print "QZ\n", QZ
		#print "RZ\n", RZ		
		
		#correctynew = matlabarray([[3],[0]])
		#print "ynew", ynew
		#self.assertEqual(ynew, correctynew)
		#print "improvement", improvement
		#correctimprovement = 0
		#self.assertEqual(improvement, correctimprovement)
		
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