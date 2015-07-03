# -*- coding: utf-8 -*-
"""
Created on Tue Jan 27 13:31:30 2015

@author: jaco_da
"""

import unittest
from runtime import *
import numpy as np
#import helper

class Test_runtime(unittest.TestCase):

	def setUp(self):
		pass

	def test_runtime_boolean_indexing1(self):
		B2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
		indfree2 = matlabarray([[0,1,2], [3,4,5],[6, 7,8]])
		#ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
		ilb2 = matlabarray([[False, False, False], [False, False, False], [False, False, False]])


		print "B2\n", B2
		print "indfree2\n", indfree2
		print "ilb2\n", ilb2


		ret = indfree2[ilb2]
		print "indfree2[ilb2]\n", ret
		
		B2[ilb2] = indfree2[ilb2]
		
		print "B2\n", B2
		
		
	def test_runtime_boolean_indexing2(self):
		B2 = matlabarray([[10, 11, 12],[13, 14, 15], [16, 17, 18]])
		indfree2 = matlabarray([[0,1,2], [3,4,5],[6, 7,8]])
		#ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
		ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])


		print "B2\n", B2
		print "indfree2\n", indfree2
		print "ilb2\n", ilb2


		ret = indfree2[ilb2]
		print "indfree2[ilb2]\n", ret
		
		B2[ilb2] = indfree2[ilb2]
		
		print "B2\n", B2
		
	def test_runtime_boolean_indexing3(self):
		B2 = matlabarray([[10, 11, 12, 13, 14, 15, 16, 17, 18]]).T
		indfree2 = matlabarray([[0,1,2, 3,4,5, 6, 7,8]]).T
		#ilb2 = matlabarray([[True, False, True], [False, False, False], [True, False, True]])
		ilb2 = matlabarray([[False, True, False, True, True, True, False, True, False]]).T


		print "B2\n", B2
		print "indfree2\n", indfree2
		print "ilb2\n", ilb2


		ret = indfree2[ilb2]
		print "indfree2[ilb2]\n", ret
		
		B2[ilb2] = indfree2[ilb2]
		
		print "B2\n", B2

	def test_runtime_boolean_indexing4(self):
		lbounds = matlabarray([[-np.inf, -np.inf, -np.inf]]).T
		ilb = matlabarray([[False, True, False]]).T
		lb = matlabarray([[-0.5, 0.0, -np.inf]]).T
		indfree = matlabarray([[1,2,3]]).T

		print "lbounds\n", lbounds
		print "ilb\n", ilb
		print "lb\n", lb
		print "indfree\n", indfree

		lbounds[ilb]=lb[indfree[ilb]]

		print "lboundsafter\n", lbounds

if __name__ == '__main__':
	unittest.main()