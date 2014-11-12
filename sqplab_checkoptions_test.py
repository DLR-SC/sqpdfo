# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 09:57:44 2014
INPUT VALUES:

K>> info

info = 

         g: [3x1 double]
        ai: []
        ae: [2x3 double]
        hl: []
     niter: 0
      flag: 0
    nsimul: [0 4 0 0 0 0 0 0 0 0 0 0 0 0 0 0]
        ci: []
        ce: [2x1 double]
         f: 1.500000000000000
      glag: [3x1 double]

K>> info.g

ans =

     0
     1
     0

K>> info.ae

ans =

     1     1     1
     1     2     3

K>> info.ce

ans =

     2
     3

K>> info.glag

ans =

  -0.333333333013890
   0.666666666736111
  -0.333333333513889

---------------------------

nb =

     2


mi =

     0


me =

     2


ms =

     0
@author: jaco_da
"""

import unittest
from sqplab_checkoptions import *
#import numpy as np
import helper

class dummyInfo():
	def __init__(self):

		self.g = matlabarray([0, 1,   0]).T
		self.ai = matlabarray([])
		self.ae =      matlabarray([[1,     1,     1],[1,     2,     3]])
		self.hl = matlabarray([])
		self.niter = 0
		self.flag = 0
		self.nsimul = matlabarray([0, 4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0])
		self.ci = matlabarray( [])
		self.ce =  matlabarray([2, 3]).T
		self.f = 1.500000000000000
		self.glag = matlabarray([-0.333333333013890, 0.666666666736111, -0.333333333513889]).T

class Test_sqplab_checkoptions(unittest.TestCase):

	def setUp(self):
		self.options = helper.dummyOptions()
		self.values = helper.dummyValues()
		self.info = dummyInfo()
		self.nb = 2
		self.mi = 0
		self.me = 2
		self.ms = 0

	def test_sqplab_checkoptions_(self):
		#print "running test"
		self.values.unit_stepsize = self.options.algo_globalization
		self.options.algo_method = self.values.quasi_newton
		
		info,options = sqplab_checkoptions_(self.nb,self.mi,self.me,self.ms,self.info,self.options,self.values)
		print "options", options.algo_descent
		print "value", self.values.powell
		self.assertTrue(options.algo_descent==self.values.powell)
		
if __name__ == '__main__':
	#print "hello"
	unittest.main()

