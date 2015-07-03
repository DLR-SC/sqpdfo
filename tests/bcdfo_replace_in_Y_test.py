# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 15:23:24 2014

@author: jaco_da
"""

import unittest
from bcdfo_replace_in_Y import *#bcdfo_replace_in_Y_
from bcdfo_build_QR_of_Y import *#bcdfo_build_QR_of_Y_
from runtime import matlabarray
import numpy as np
#import helper

class Test_bcdfo_replace_in_Y(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	#@unittest.expectedFailure
	def test_bcdfo_replace_in_Y(self):
		Y = matlabarray([[ 1.0, 1.0, 0.0, 0.0, 0.0, -1.0], [1.0, 0.0, -1.0, 1.0, 0.0, 0.0, ]])
		whichmodel = 0
		QZ,RZ,xbase,scale = bcdfo_build_QR_of_Y( Y , 0, 1, 1, 1, 1e15 )
		#print "QZ\n", QZ
		#print "RZ\n", RZ

		#RZ = np.array([

    #[1.0000,    1.0000,    1.0000,    1.0000,   1.0000,    1.0000],
    # [    0.0,   -0.4583,   -0.9602,         0,  -0.4583,   -0.4583],
    #  [   0.0,         0,    0.6388,    0.3287,    0.4540,    0.9392],
    #  [   0,         0,         0,   -0.3193,   -0.1903,   -0.4111],
    #  [   0,         0,         0,         0,    0.0877,    0.1754],
    #  [   0,         0,         0,         0,         0,    0.1952		]])
		
		ynew = 0.25*Y[:,3].T
		#print "ynew", ynew
		QZplus, RZplus, Yplus, xbase, scale = bcdfo_replace_in_Y( QZ, RZ, ynew, Y, 3, xbase, whichmodel, scale, 1, 1, 1, 1e15 )
		
		#print "QZplus:\n", QZplus, "\nRZplus:\n",  RZplus, "\nYplus:\n", Yplus
		
		correctYplus = matlabarray([
    [1.0000,    1.0000,         0,         0,         0,   -1.0000],
    [1.0000,         0,   -0.2500,    1.0000,         0,         0]])
				
		self.assertTrue((abs(Yplus - correctYplus) < 1e-8).all())
		
				
#		correctQZplus = matlabarray([
#   [-0.4714,   -0.4714,    0.7205,    0.1527,    0.1149,    0.0000],
#   [-0.4714,   -0.4714,   -0.5764,   -0.1221,   -0.0919,    0.4472],
#   [-0.4714,    0.4714,   -0.1891,    0.6333,    0.3446,    0.0000],
#   [-0.2357,   -0.2357,   -0.2882,   -0.0611,   -0.0459,   -0.8944],
#   [-0.2357,    0.2357,    0.1081,    0.1813,   -0.9189,   -0.0000],
#   [-0.4714,    0.4714,    0.1351,   -0.7240,    0.1149,   -0.0000]])

#		correctRZplus = matlabarray([
#   [-2.1213,   -1.0607,   -0.3609,   -1.0607 ,  -0.4714,   -0.1179],
#   [      0,   -1.0607,   -0.5819,    0.1179,   -0.4714,   -0.1179],
#   [      0,         0,  0.7711 ,   0.5854 ,   0.7205 ,   1.1527],
#   [      0,         0,        0,    0.8766 ,   0.1527,    0.2442],
#   [      0,         0,         0,         0,    0.1149,    0.1838],
#   [      0,        0,         0 ,        0,         0 ,  -0.8944]])

#		correctYplus = matlabarray([
#    [1.0000,    1.0000,         0,         0,         0,   -1.0000],
#    [1.0000,         0,   -0.2500,    1.0000,         0,         0]])
		
		
	def est_bcdfo_replace_in_YOriginal(self):
		Y = np.array([[ 1.0, 1.0, 0.0, 0.0, 0.0, -1.0], [1.0, 0.0, -1.0, 1.0, 0.0, 0.0, ]])
		whichmodel = 0
		QZ,RZ,xbase,scale = bcdfo_build_QR_of_Y( Y , 0, 1, 1, 1, 1e15 )
		#print "QZ\n", QZ
		#print "RZ\n", RZ

		#RZ = np.array([

    #[1.0000,    1.0000,    1.0000,    1.0000,   1.0000,    1.0000],
    # [    0.0,   -0.4583,   -0.9602,         0,  -0.4583,   -0.4583],
    #  [   0.0,         0,    0.6388,    0.3287,    0.4540,    0.9392],
    #  [   0,         0,         0,   -0.3193,   -0.1903,   -0.4111],
    #  [   0,         0,         0,         0,    0.0877,    0.1754],
    #  [   0,         0,         0,         0,         0,    0.1952		]])
		
		ynew = np.array([0.25*Y[:,3-1]]).T
		#print "ynew", ynew
		QZplus, RZplus, Yplus, xbase, scale = bcdfo_replace_in_Y( QZ, RZ, ynew, Y, 3-1, xbase, whichmodel, scale, 1, 1, 1, 1e15 )
		
		#print "QZplus:\n", QZplus, "\nRZplus:\n",  RZplus, "\nYplus:\n", Yplus
		
		correctYplus = np.array([
    [1.0000,    1.0000,         0,         0,         0,   -1.0000],
    [1.0000,         0,   -0.2500,    1.0000,         0,         0]])
				
		self.assertTrue((abs(Yplus - correctYplus) < 1e-8).all())

if __name__ == '__main__':
	unittest.main()