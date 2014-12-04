# -*- coding: utf-8 -*-
"""
Created on Tue Dec 02 15:13:04 2014

@author: jaco_da
"""

import unittest
from bcdfo_projgrad import bcdfo_projgrad_
from runtime import matlabarray
#import numpy as np
#import helper

class Test_bcdfo_projgrad(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_bcdfo_projgrad(self):
		print "No test specified..."
#		  n     : dimension
#  x     : current iterate
#  g     : gradient at x
#  bl    : lower bounds
#  bu    : upper bounds
		n = 2
		x = matlabarray([0.5,0.5]).T
		g = matlabarray([-1,-1])#.T
		bl = matlabarray([0,0]).T
		bu = matlabarray([2,2]).T
		
		gnorm, gn = bcdfo_projgrad_(n,x,g,bl,bu)
		
		print "gnorm", gnorm
		print "gn", gn
		
		correctgn = matlabarray([-1,-1]).T
		correctgnorm = 1
		
		self.assertTrue((abs(correctgn - gn) < 1e-8).all())
		self.assertEqual(gnorm, correctgnorm)

if __name__ == '__main__':
	unittest.main()