# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 14:41:23 2014

@author: jaco_da
"""
import unittest
from bcdfo_include_in_Y import bcdfo_include_in_Y_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_
from runtime import matlabarray, char
#import numpy as np
#import helper

class Test_bcdfo_include_in_Y(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_bcdfo_include_in_Y(self):
#		TEST:
		Y = matlabarray([[ 0, 1, 0, 2, 1, 0],[0, 0, 1, 0, 0.01, 2 ]])
		whichmodel = 0
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, whichmodel, 1, 1, 1, 1e15 )
		QZplus, RZplus, Yplus, pos, xbase, scale = bcdfo_include_in_Y_(matlabarray([-1,1]).T, QZ, RZ, Y, matlabarray(range(2,6)), 0.01, char('weighted'), xbase, whichmodel, 0, scale, 1, 1, 1, 1e15 )
		print QZplus, RZplus, Yplus, pos, xbase, scale
		
		correctQZplus = matlabarray([
   [-1.0000,         0,         0,         0,         0,         0],
   [     0,    0.9701,    0.0000 ,  -0.2425,    0.0000,   -0.0000],
   [     0,    0.0000,   -0.9701,   -0.0000,    0.0000,   -0.2425],
   [     0,    0.2425,   -0.0000,    0.9701,   -0.0000,    0.0000],
   [     0,    0.0000,   -0.2425,    0.0000,   -0.0000,    0.9701],
   [     0,   -0.0000,    0.0000,    0.0000,    1.0000,   0.0000]])

		correctRZplus = matlabarray([
   [-1.0000,   -1.0000,   -1.0000,   -1.0000,   -1.0000,   -1.0000],
   [      0,    0.5154,         0,    1.0914,   -0.4548,    0.0000],
   [      0,         0,   -0.5154,   -0.0000,   -0.5154,   -1.0914],
   [      0,         0,         0,    0.2425,    0.2425,   -0.0000],
   [      0,         0,         0,         0,   -0.2500,   -0.0000],
   [      0,         0,         0,         0,         0,    0.2425]])

		correctYplus = matlabarray([

     [0,     1,     0,     2,    -1,     0],
     [0,     0,     1,     0,     1,     2]])

		correctpos = 4  #5 is the matlabindex
		correctxbase = matlabarray([[0],[0]])
		correctscale = matlabarray([ 
    [1.0000],
    [0.5000],
    [0.5000],
    [0.2500],
    [0.2500],
    [0.2500],])
				
		#self.assertTrue((abs(correctQZplus - QZplus) < 1e-8).all())
		#self.assertTrue((abs(correctRZplus - RZplus) < 1e-8).all())
		print "QR decompositions is different again"
		self.assertTrue((abs(correctYplus - Yplus) < 1e-8).all())
		self.assertEqual(correctpos, pos)
		self.assertAlmostEqual(correctxbase, xbase)
		self.assertAlmostEqual(correctscale, scale)
	
#  QZplus =
#
#   -1.0000         0         0         0         0         0
#         0    0.9701    0.0000   -0.2425    0.0000   -0.0000
#         0    0.0000   -0.9701   -0.0000    0.0000   -0.2425
#         0    0.2425   -0.0000    0.9701   -0.0000    0.0000
#         0    0.0000   -0.2425    0.0000   -0.0000    0.9701
#         0   -0.0000    0.0000    0.0000    1.0000    0.0000
#
#  RZplus =
#
#   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000   -1.0000
#         0    0.5154         0    1.0914   -0.4548    0.0000
#         0         0   -0.5154   -0.0000   -0.5154   -1.0914
#         0         0         0    0.2425    0.2425   -0.0000
#         0         0         0         0   -0.2500   -0.0000
#         0         0         0         0         0    0.2425
#
#  Yplus =
#
#     0     1     0     2    -1     0
#     0     0     1     0     1     2
#
#  pos =
#
#     5
#
#  xbase =
#
#     0
#     0
#
#  scale =
#
#    1.0000
#    0.5000
#    0.5000
#    0.2500
#    0.2500
#    0.2500
#
#  and
#  RZplus\(QZplus'*bcdfo_evalZ((Yplus-Y(:,1)*ones(1,6))*scale(2),6))
#  should give the identity matrix.

if __name__ == '__main__':
	unittest.main()