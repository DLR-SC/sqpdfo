# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 14:54:04 2014

@author: jaco_da
"""

import unittest
from bcdfo_evalL import bcdfo_evalL_
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_, bcdfo_build_QR_of_Y
from runtime import matlabarray
import numpy as np
#import helper

class Test_bcdfo_evalL(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_bcdfo_evalL(self):
		Y = matlabarray([[ 0, 1, 0, 2, 1, 0],[0, 0, 1, 0, 0.01, 2 ]]) 
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y_( Y, 0, 1, 1,1,1e15 )
		values = bcdfo_evalL_( QZ, RZ, Y, matlabarray(range(1,6)), matlabarray([[-1],[1]]), xbase, 0, scale, 1 )
		print "values", values
		correctvalues = matlabarray([ [0],   [97.0000],    [2.9900],    [1.0000], [-100.0000],   [-0.4950] ])
		self.assertTrue((abs(values - correctvalues) < 1e-8).all())
		
	def est_bcdfo_evalLOriginal(self):
		Y = np.array([[ 0, 1, 0, 2, 1, 0],[0, 0, 1, 0, 0.01, 2 ]]) 
		QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y( Y, 0, 1, 1,1,1e15 )
		values = bcdfo_evalL_( QZ, RZ, Y, np.array(range(1,6)), np.array([[-1], [1]]), xbase, 0, scale, 1 )
		print "values", values
		

#Ankes Test
#Y = np.array([[ 0, 1, 0, 2, 1, 0],[ 0, 0, 1, 0, 0.01, 2 ]])
#n,p1=np.shape(Y)
#I=np.eye(p1)
#lc     = len(np.array([1,2,3,4,5]))
#q=( ( n + 1 ) * ( n + 2 ) ) / 2
#values = np.zeros((p1,1))
#from bcdfo_build_QR_of_Y import *
#[QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y, 0, 1, 1, 1, 1e15 )
#from bcdfo_evalL import *
#x = np.array([[-1.],[1.]])
#values = bcdfo_evalL( QZ, RZ, Y, np.array([1,2,3,4,5]), x, xbase, 0, scale, 1 )	

##  should give 
##     values =
##         0   97.0000    2.9900    1.0000 -100.0000   -0.4950 

if __name__ == '__main__':
	unittest.main()