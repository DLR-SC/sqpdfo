# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 16:28:06 2014
INPUT VALUES
key,xy,lm

key =

     2


xy =

    -2
     2
     2
     1
     1

lm = none (Not enough input arguments.)
OUTPUT VALUES
 outdic,fvalue,info.ci,info.ce

outdic =

     0


fvalue =

     3.354626279025119e-04


ans =

     []


ans =

     4
    -1
     1

@author: jaco_da
"""

import unittest
from evalfgh import *
from ecdfo_func import *
#import numpy as np
#import helper

class Test_evalfgh(unittest.TestCase):

	def setUp(self):
		self.key = 2
		self.xy = matlabarray([-2, 2, 2, 1, 1]).T
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()
		
		set_simul_not_initialized(0)
		set_prob(5)

	def test_evalfgh(self):
		msg,out2,out3,out4 = evalfgh_(self.key,self.xy,lm=None)

		#print "msg,out2,out3,out4", msg,out2,out3,out4
		
		correctmsg = 0
		correctout2 = 3.354626279025119e-04
		correctout3 =  matlabarray([])
		correctout4 =    matlabarray([4, -1, 1]).T
		
		self.assertEqual(correctmsg, msg)
		self.assertAlmostEqual(correctout2, out2, 8)
		self.assertEqual(correctout3, out3)
		
		self.assertTrue((abs(correctout4 - out4) < 1e-7).all())
		
		#self.assertEqual(correctout4, out4)


if __name__ == '__main__':
	unittest.main()
	