# -*- coding: utf-8 -*-
"""
Created on Thu Nov 06 11:26:49 2014
check if function returns convex matrix for convex matrix (1) and non convex matrix (2) and if matrix is non symmetric (3)
@author: jaco_da
"""
import unittest
from ecdfo_check_convex import *
import helper
import numpy as np

class Test_ecdfo_check_convex(unittest.TestCase):

	#def setUp(self):
	#    self.dummyOptions = helper.dummyOptions()

	def test_ecdfo_check_convex1(self):
		A = np.array([[ 2,  -1, 0], [ -1,  2,  -1], [ 0,  -1,  2]])

		res = ecdfo_check_convex_(A)
		#print "A", A
		#print "ecdfo  check cond", str(res)

		self.assertTrue(np.all(np.linalg.eigvals(res) > 0))

	def test_ecdfo_check_convex2(self):
		B = np.array([[ 1,  2], [ 2,  1]])

		res = ecdfo_check_convex_(B)
		#print "B", B
		#print "ecdfo  check cond", str(res)

		self.assertTrue(np.all(np.linalg.eigvals(res) > 0))
		
	def test_ecdfo_check_convex3(self):
		C = np.array([[1,2,3],[5,4,6],[7,8,9]])
		
		res = ecdfo_check_convex_(C, helper.dummyOptions())
		
		#print "nonsymmetric res", res
		
		self.assertTrue(np.all(np.linalg.eigvals(res) > 0))
		self.assertTrue((res.T == res).all())

if __name__ == '__main__':
#	suite = unittest.TestLoader().loadTestsFromTestCase(Test_ecdfo_check_convex)
#	unittest.TextTestRunner(verbosity=2).run(suite)
	unittest.main()
