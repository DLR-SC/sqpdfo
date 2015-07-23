# -*- coding: utf-8 -*-
"""
Created on Wed Nov 12 14:00:19 2014
INPUT VALUES

A,b,delta,max_iter,tol,plevel,fout

A =

     1


b =

   0.816496580927725


delta =

     1


max_iter =

    20


tol =

     1.000000000000000e-06


plevel =

     0


fout =

     1
					
OUTPUT VALUES.

 u,info_t

u =

   0.816496580927725


info_t = 

    flag: 0
    iter: 2
    prec: 0
    curv: 1
					

@author: jaco_da
"""
import sys
sys.path.append("../")
import unittest
from sqplab_tcg import *
#import numpy as np
#import helper

class Test_sqplab_tcg(unittest.TestCase):

	def setUp(self):
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()
	
		self.A = matlabarray([ 1])
		self.b = matlabarray([ 0.816496580927725])
		self.delta =1
		self.max_iter = 20
		self.tol =1.000000000000000e-06
		self.plevel = 0
		self.fout = 1

	def test_sqplab_tcg(self):
		u,info_t = sqplab_tcg_(self.A,self.b,self.delta,self.max_iter,self.tol,self.plevel,self.fout)
		
		#print "u", u
		self.assertEqual(u, 0.816496580927725)

		self.assertEqual(info_t.flag, 0)
		self.assertEqual(info_t.iter, 2)
		self.assertEqual(info_t.prec, 0)
		self.assertEqual(info_t.curv, 1)

if __name__ == '__main__':
	unittest.main()