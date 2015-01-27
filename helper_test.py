# -*- coding: utf-8 -*-
"""
Created on Wed Dec 03 11:34:49 2014

@author: jaco_da
"""

import unittest
#from helper import *
from runtime import matlabarray
#import numpy as np
import helper

class Test_convertArray(unittest.TestCase):

	def setUp(self):
		self.x = matlabarray([1,2,3,4]).T
		self.y = matlabarray([5,6,7,8])
		self.z = matlabarray([1.1478])
		
		self.A = matlabarray([[1,2,3],[4,5,6],[7,8,9]])

	def est_convertArray_reshape(self):
		x = helper.convertArray(self.x)
		y = helper.convertArray(self.y)
		z = helper.convertArray(self.z)
		
		A = helper.convertArray(self.A)
		
		#print "Conversion done:"
		#print "x\n", x
		#print "y\n", y
		#print "z\n", z
		
		#print "A\n", A
		
		x.reshape((4,1))
		y.reshape((2,2))
		A.reshape((1,4))
		
		#print "Reshape done:"
		#print "x\n", x
		#print "y\n", y
		#print "z\n", z
		
		#print "A\n", A


	def test_convertArray_resize(self):
		x = helper.convertArray(self.x)
		y = helper.convertArray(self.y)
		z = helper.convertArray(self.z)
		
		A = helper.convertArray(self.A)
		
		#print "Conversion done:"
		#print "x\n", x
		#print "y\n", y
		#print "z\n", z
		
		#print "A\n", A
		
		x.resize((4,2))
		y.resize((2,4))
		z.resize((2,2))
		A.resize((4,4))
		
		#print "Resize done:"
		#print "x\n", x
		#print "y\n", y
		#print "z\n", z
		
		#print "A\n", A


if __name__ == '__main__':
	unittest.main()