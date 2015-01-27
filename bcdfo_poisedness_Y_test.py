# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 13:54:53 2014

@author: jaco_da
"""
import unittest
from bcdfo_poisedness_Y import *#bcdfo_poisedness_Y_
from bcdfo_build_QR_of_Y import *#bcdfo_build_QR_of_Y
import numpy as np
from runtime import matlabarray
#import helper

class Test_bcdfo_poisedness_Y(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass
	
	def test_bcdfo_poisedness_Y(self):
		
		Y = matlabarray([ [ 0, 1, 0, 2, 1, 0 ], [ 0, 0, 1, 0, 0.01, 2] ])
		#print Y
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1, 1, 1e15 )
		#QZ = matlabarray(QZ)
		#RZ = matlabarray(RZ)
		#xbase = matlabarray(xbase)
		#scale = matlabarray(scale)
		


		#0, matlabarray([-10,-10]), matlabarray([10,10]), matlabarray([1,2]), 1

		lSolver = 1
		
		hardcons = 0 
		xl = matlabarray([-10,-10])
		xu = matlabarray([10,10])
		indfree = matlabarray([1,2])
		stratLam = 1
		
		#print "Warning: Not enough input arguments"
		lambd ,Y_radius  = bcdfo_poisedness_Y_( QZ, RZ, Y, 0.001, xbase, lSolver, 1, hardcons, xl, xu, indfree, stratLam,  scale, 1 )
		#print "lambd = ", lambd
		#print "Y_radius = ", Y_radius
		
		correctlambda = 402.6381
		correctY_radius = 2
		self.assertAlmostEqual(Y_radius, correctY_radius)
		self.assertAlmostEqual(lambd, correctlambda, 4)
		

	#@unittest.expectedFailure
	def est_bcdfo_poisedness_YOriginal(self):
		#print "Self-test of the function with shifting:"
		Y = np.array([ [ 0, 1, 0, 2, 1, 0 ], [ 0, 0, 1, 0, 0.01, 2] ])
		#print Y
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y( Y, 0, 1, 1, 1, 1e15 )
		#QZ = matlabarray(QZ)
		#RZ = matlabarray(RZ)
		#xbase = matlabarray(xbase)
		#scale = matlabarray(scale)
		


		#0, matlabarray([-10,-10]), matlabarray([10,10]), matlabarray([1,2]), 1

		lSolver = 1
		
		hardcons = 0 
		xl = np.array([-10,-10])
		xu = np.array([10,10])
		indfree = np.array([1,2])
		stratLam = 1
		
		#print "Warning: Not enough input arguments"
		lambd ,Y_radius  = bcdfo_poisedness_Y( QZ, RZ, Y, 0.001, xbase, lSolver, 1, hardcons, xl, xu, indfree, stratLam,  scale, 1 )
		#print "lambd = ", lambd
		#print "Y_radius = ", Y_radius

		
		#print lambd, Y_radius
		#self.asserTrue(False)
#  should give:
#  lambd =
#  204.8586
#  Y_radius =
#     2

if __name__ == '__main__':
	unittest.main()