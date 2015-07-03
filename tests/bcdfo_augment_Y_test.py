# -*- coding: utf-8 -*-
"""
Created on Mon Dec 01 11:12:48 2014

@author: jaco_da
"""

import unittest
from bcdfo_augment_Y import bcdfo_augment_Y_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray
#import numpy as np
#import helper

class Test_bcdfo_augment_Y(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_bcdfo_augment_Y(self):
		Y = matlabarray([[ 0, 1, 0, 2, 0], [0, 0, 1, 0, 2 ]])
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, 0, 0, 1, 1, 1e15 )
		p1, QZ, RZ, Y, xbase, scale = bcdfo_augment_Y_( matlabarray([1, 0.01]).T, Y, 0, 0, 1, 1, 1e15 )
		#print "p1, QZ, RZ, Y, xbase, scale", p1, QZ, RZ, Y, xbase, scale
		
		correctY = matlabarray([
         [0,    1.0000,         0,    2.0000,         0,    1.0000],
         [0,         0,    1.0000,         0,    2.0000,    0.0100]])

		correctxbase = matlabarray([ 0, 0]).T
		correctscale = matlabarray([1, 1, 1, 1, 1, 1]).T
		
		self.assertEqual(p1, 6)
		self.assertEqual(Y, correctY)
		self.assertEqual(xbase, correctxbase)
		self.assertEqual(scale, correctscale)
		
#  gives
#  QZ =
#
#    1.0000         0         0         0         0         0
#         0   -0.8944         0   -0.4472         0         0
#         0         0   -0.8944         0   -0.4472         0
#         0   -0.4472         0    0.8944         0         0
#         0         0   -0.4472         0    0.8944         0
#         0         0         0         0         0    1.0000
#
#  RZ =
#
#    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000
#         0   -1.1180         0   -2.6833         0   -1.1180
#         0         0   -1.1180         0   -2.6833   -0.0090
#         0         0         0    0.8944         0         0
#         0         0         0         0    0.8944   -0.0044
#         0         0         0         0         0    0.0100


if __name__ == '__main__':
	unittest.main()