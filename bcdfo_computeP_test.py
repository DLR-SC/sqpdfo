# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 13:39:12 2014

@author: jaco_da
"""


import unittest
from bcdfo_computeP import *
from bcdfo_build_QR_of_Y import *
from runtime import *
#import numpy as np
#import helper

class Test_bcdfo_computeP(unittest.TestCase):

	def setUp(self):
		pass
	
	def test_bcdfo_computeP(self):
		Y = matlabarray([[ 0, 1, 0, 2, 0],[0, 0, 1, 0, 2]])
		fY =matlabarray([[ 1, 2, 3, 4, 5 ]])
		whichmodel = 0
		
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_(Y, whichmodel, 1, 1, 1, 1e15 )
		#QZ = matlabarray(QZ)
		#RZ = matlabarray(RZ)
		
		#xbase = matlabarray(xbase)
		#scale = matlabarray(scale)
		P = bcdfo_computeP_( QZ, RZ, Y, fY, whichmodel, matlabarray([0, 0, 0, 0, 0]), 0,
			1, 1, matlabarray([0, 0]), scale, 1, 1, matlabarray([1, 2]), matlabarray([]), matlabarray([[1, 0],[0, 1]]), 1e-5, None )

		#print "P\n", P
		#  should give:
#  P =
#   1.0000   1.0000   4.0000   4.0000       0
		
	def est_bcdfo_computePOriginal(self):
		#print "Testing original func"
		Y = np.array([[ 0, 1, 0, 2, 0],[0, 0, 1, 0, 2]])
		fY = np.array([ 1, 2, 3, 4, 5 ])
		whichmodel = 0
		
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y(Y, whichmodel, 1, 1, 1, 1e15 )
		#print "QZ\n", QZ
		#print "RZ\n", RZ
	
		P = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, np.array([0, 0, 0, 0, 0]), 0,
		1, 1, np.array([0, 0]), scale, 1, 1, np.array([1, 2]), np.array([]), np.array([[1, 0],[0, 1]]), 1e-5, None )

		#print "P\n", P
	

if __name__ == '__main__':
	unittest.main()