# -*- coding: utf-8 -*-
"""
Created on Wed Nov 26 13:35:53 2014

@author: jaco_da
"""

import unittest
from bcdfo_evalZ import *
#import numpy as np
#import helper

class Test_bcdfo_evalZ(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	@unittest.expectedFailure
	def test_bcdfo_evalZ(self):
		res = bcdfo_evalZ_(matlabarray([[ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]]), 6)
		print "res", res
		correctres = matlabarray([
    [1.0000,    1.0000,    1.0000,    1.0000,    1.0000,    1.0000],
     [    0,    1.0000,         0,    2.0000,    1.0000,         0],
      [   0,         0,    1.0000,         0,    0.0100,    2.0000],
       [  0,    0.5000,         0,    2.0000,    0.5000,         0],
        [ 0,         0,    0.5000,         0,    0.0001,    2.0000],
         [0,         0,         0,         0,    0.0100,         0]])
			
		print "abs:\n\n", abs(res - correctres)
		print "5e-5  -> 1e-1 ?"
		self.assertTrue((abs(res - correctres) < 1e-8).all())


#  should give
# ans =
#
#    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000
#         0    1.0000         0    2.0000    1.0000         0
#         0         0    1.0000         0    0.0100    2.0000
#         0    0.5000         0    2.0000    0.5000         0
#         0         0    0.5000         0    0.0001    2.0000
#         0         0         0         0    0.0100         0
#

if __name__ == '__main__':
	unittest.main()