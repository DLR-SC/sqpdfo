# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 15:16:08 2014

@author: jaco_da
"""

import unittest
from bcdfo_repair_Y import bcdfo_repair_Y_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray
#import numpy as np
#import helper

class Test_bcdfo_repair_Y(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_bcdfo_repair_Y(self):
#		TEST:
		Y = matlabarray([[ 0, 1, 0, 2, 1, 0],[0, 0, 1, 0, 0.01, 2 ]])
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, 0, 1, 1, 1, 1e15 )
		QZplus, RZplus, Yplus, replaced, maximprove, Y_radius, xbase, scale = bcdfo_repair_Y_( QZ, RZ, Y, 0.7, 10, 1.0e-10, 1.1, 0.001, xbase, 1, 0, 0, matlabarray([-10,-10]), matlabarray([10,10]), matlabarray([1,2]), 1, scale, 1, 1, 1e15 )
		#print QZplus, RZplus, Yplus, replaced, maximprove, Y_radius, xbase, scale
		
		correctQZplus = matlabarray([

   [-1.0000,         0,         0,        0 ,        0 ,        0],
   [      0,    0.7017,    0.6633,    0.1499,   -0.2060,   -0.0529],
   [      0,    0.6810,   -0.7144,    0.0589,    0.0925,   -0.1178],
   [      0,   -0.0881,    0.0262,    0.9003,    0.4120,    0.1058],
   [      0,   -0.0830,    0.1119,   -0.0352,    0.2945,   -0.9448],
   [      0,   -0.1710,   -0.1910,    0.4028,   -0.8322,   -0.2820]])


		correctRZplus = matlabarray([

  [-1.0000,  -1.0000,  -1.0000,  -1.0000,  -1.0000,  -1.0000],
  [      0,  -0.3580,  -0.0762,   0.6576,  -0.0867,   0.6395],
  [      0,        0,   0.3491,   0.6764,  -0.3138,  -0.6584],
  [      0,        0,        0,   0.6001,  -0.0159,   0.0413],
  [      0,        0,        0,        0,   0.1465,   0.2397],
  [      0,        0,        0,        0,        0,  -0.5902]])

		correctYplus = matlabarray([

        [0,  -0.5023,   0.3561,   2.0000,  -0.6030,        0],
        [0,  -0.4875,  -0.6026,        0,   0.3555,   2.0000]])

		correctreplaced = matlabarray([1,    2,    4]) # changed to python indices

		correctmaximprove =  1.0002

		correctY_radius =  2
  
		correctxbase = matlabarray([
    [0],
    [0]])

		correctscale = matlabarray([
   [1.0000],
   [0.5000],
   [0.5000],
   [0.2500],
   [0.2500],
   [0.2500]])
			
		#self.assertTrue((abs(correctQZplus - QZplus) < 1e-3).all())
		#self.assertTrue((abs(correctRZplus - RZplus) < 1e-3).all())
		#print "Warning not unique QR (?)"
		#print "Yplus\n", Yplus
		#print "correctYplus\n", correctYplus			
		self.assertTrue((abs(correctYplus - Yplus) < 1e-4).all())
		
		self.assertEqual(correctreplaced, replaced)
		self.assertAlmostEqual(correctmaximprove, maximprove, 3)
		
#-------these are correct---------:		
		self.assertTrue((abs(xbase - correctxbase) < 1e-4).all())
		self.assertEqual(correctY_radius, Y_radius)
		self.assertTrue((abs(correctscale - scale) < 1e-4).all())
		

if __name__ == '__main__':
	unittest.main()