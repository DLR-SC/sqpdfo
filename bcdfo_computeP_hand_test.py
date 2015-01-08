# -*- coding: utf-8 -*-
"""
Created on Tue Nov 25 13:39:12 2014

@author: jaco_da
"""


import unittest
from bcdfo_computeP_hand import *
from bcdfo_build_QR_of_Y import *
from runtime import *
#import numpy as np
#import helper

class Test_bcdfo_computeP(unittest.TestCase):

	def setUp(self):
		pass
	
	def est_bcdfo_computeP(self):
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

		print "P\n", P
		#  should give:
#  P =
#   1.0000   1.0000   4.0000   4.0000       0
		
	def test_bcdfo_computePOriginal(self):
		print "Testing original func"
		Y = np.array([[ 0, 1, 0, 2, 0],[0, 0, 1, 0, 2]])
		fY = np.array([ 1, 2, 3, 4, 5 ])
		whichmodel = 0
		
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y(Y, whichmodel, 1, 1, 1, 1e15 )
		print "QZ\n", QZ
		print "RZ\n", RZ
	
		P = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, np.array([0, 0, 0, 0, 0]), 0,
		1, 1, np.array([0, 0]), scale, 1, 1)#, np.array([1, 2]), np.array([]), np.array([[1, 0],[0, 1]]), 1e-5, None )

		print "P\n", P
		
	def test_bcdfo_computePOriginal2(self):
		QZ = np.array([
     [1,     0,     0,     0],
     [0,     1,     0,     0],
     [0,     0,     1,     0],
     [0,     0,     0,     1]])


		RZ = np.array([

     [1,     1,     1,     1],
     [0,    -1,     0,     0],
     [0,     0,    -1,     0],
     [0,     0,     0,    -1]])

		Y = np.array([

   [0.500000000000000,  -0.500000000000000,   0.500000000000000,   0.500000000000000],
   [1.000000000000000,   1.000000000000000,                   0,   1.000000000000000],
   [0.500000000000000,   0.500000000000000,   0.500000000000000,  -0.500000000000000]])

		fY = np.array([

   [1.500000000000000,   1.500000000000000,   0.500000000000000,   1.500000000000000],
   [2.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000],
   [3.000000000000000,   2.000000000000000,   1.000000000000000,                   0]])

		whichmodel = 0

		P_old = np.array([0,     0,     0,     0,     0,     0,     0,     0,     0,     0])
		ind_Y = np.array( [1,     2,     3,     4])
		i_xold = 0
		i_xplus = 0
		g = np.NaN
		scale = np.array( [  [1],  [1],  [1],  [1]])
		shift_Y = 1
		Delta0 = 1
		
		P = bcdfo_computeP( QZ, RZ, Y, fY, whichmodel, P_old, ind_Y, i_xold, i_xplus, g, scale, shift_Y, Delta0 )
		print P
		
		correctP =  np.array([[1.500000000000000,                   0 ,  1.000000000000000,                   0],
						[2.000000000000000,   1.000000000000000,   1.000000000000000,   1.000000000000000],
						[3.000000000000000,   1.000000000000000,   2.000000000000000,   3.000000000000000]])
						
		self.assertTrue((abs(P - correctP) < 1e-8).all())						
	

if __name__ == '__main__':
	unittest.main()