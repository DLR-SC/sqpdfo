# -*- coding: utf-8 -*-
"""
Created on Mon Dec 01 13:41:32 2014

@author: jaco_da
"""

import unittest
from bcdfo_computeLj import *#bcdfo_computeLj_
from bcdfo_build_QR_of_Y import *#bcdfo_build_QR_of_Y_
from runtime import matlabarray
import numpy as np
#import helper

class Test_bcdfo_computeLj(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	#@unittest.expectedFailure
	def test_bcdfo_computeLj(self):
		Y = matlabarray([[ 0, 1, 0, 2, 0], [0, 0, 1, 0, 2 ]])#matlabarray([[14.0, 32.0, 50.0], [32.0, 77.0, 122.0]])#, [50.0, 122.0, 194.0]])
		whichmodel = 0
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(  Y, whichmodel, 1, 1, 1, 1e15 )
		Lj = bcdfo_computeLj_( QZ, RZ, 1, Y, whichmodel, scale, 1 )
		print "Lj", Lj

		#print "type(scale)", type(scale)
		#print "Lj", Lj
		correctLj = matlabarray([1.0000,  -3.0000,  -3.0000,   4.0000,   4.0000])#matlabarray([0.833333333333333,  -4.455569580721621,  -0.146502319909943])
		#print "Warning: Different Results due to not unique QR decomposition (?)"
		self.assertTrue((abs(Lj - correctLj) < 1e-8).all())
		
	def est_bcdfo_computeLj_Original(self):
		Y = np.array([[ 0.0, 1.0, 0.0, 2.0, 0.0], [0.0, 0.0, 1.0, 0.0, 2.0 ]])
		whichmodel = 0
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y(  Y, whichmodel, 1, 1, 1, 1e15 )
		#print "QZ\n", QZ
		#print "RZ\n", RZ 
		Lj = bcdfo_computeLj( QZ, RZ, 1, Y, whichmodel, scale, 1 )
#  should give:
		correctLj = np.array([1.0000,  -3.0000,  -3.0000,   4.0000,   4.0000])
		print Lj #- correctLj
			
			
			
# should give:
#  Lj =
#   1.0000  -3.0000  -3.0000   4.0000   4.0000

if __name__ == '__main__':
	unittest.main()