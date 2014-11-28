# -*- coding: utf-8 -*-
"""
Created on Thu Nov 13 11:35:42 2014
TEST:
%  X = [ 0 1 0 ; 0 0 1 ];
%  fX = [ 1 100 101 ];
%  y = [ 2; 4 ];
%
%  To put y on 4-th position in the set X call:
%  [X, fX, neval, xstatus, sstatus, dstatus] = bcdfo_augmX_evalf(@banana, ...
%     y, 4, X, fX, 0, [0;0], [], [1 2 3], 1e25, 3, 100, 1, ...
%     0, [1 1 1], 1, [1 1 1], [ 0 0 0], 0, [ 1, 1])
%  where
%     function fx = banana( x )
%     fx  = 100 * ( x(2) - x(1)^2 ) ^2 + (1-x(1))^2;
%
%  which gives
%  X =
%     0     1     0     2
%     0     0     1     4
%  fX =
%     1     100     101   1
%  neval =
%     4
%  xstatus =
%     1     1     1     1
%  sstatus =
%     1     1     1     1
%  dstatus =
%     0     0     0     0
@author: jaco_da
"""

import unittest
from ecdfo_augmX_evalf import *
#import numpy as np
#import helper

class dummyInfo():
	def __init__(self):
		self.nsimul = matlabarray([0,0,0])
		self.ci = matlabarray([])
		self.ce = matlabarray([])

class Test_ecdfo_augmX_evalf(unittest.TestCase):

	def setUp(self):
		self.X = matlabarray([[ 0, 1, 0],[0, 0, 1] ])
		self.fX = matlabarray([ 1, 100, 101 ])
		self.y = matlabarray([ 2, 4 ])
		
		self.info = dummyInfo()
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()
	def banana(self, x, y):
		#return outdic,fvalue,info.ci,info.c 
		print "x", x
		print "y", y
		return None,100 * ( y[2] - y[1]^2 ) ^2 + (1-y[1])^2,matlabarray([]),matlabarray([])

#	@unittest.expectedFailure
	def test_ecdfo_augmX_evalf(self):
		#(f=None,y=None,m=None,X=None,fX=None,ciX=None,ceX=None,nfix=None,xfix=None,indfix=None,
		#indfree=None,fxmax=None,neval=None,xstatus=None,xstatus_val=None,sstatus=None,dstatus=None,
		#scaleX=None,scalefacX=None,info=None,options=None,values=None)
		#returns X, fX, neval, xstatus, sstatus, dstatus
	
		X,fX,ciX,ceX,neval,xstatus,sstatus,dstatus,info,outdic = ecdfo_augmX_evalf_(
		f=self.banana,y=self.y, m=4, X=self.X,
		fX=self.fX,ciX=0, ceX=matlabarray([[0],[0]]), nfix=matlabarray([]), xfix=matlabarray([1, 2, 3]), 
		indifx=1e25, indfree=3, fxmax=100, neval=1, xstatus_val=0, xstatus=matlabarray([1, 1, 1]),
		sstatus=matlabarray([1, 1, 1]), dstatus=matlabarray([ 0, 0, 0]), scaleX=0, 
		scalefacX=matlabarray([ 1, 1]), info=self.info)
		
		#sstatus=1, dstatus=matlabarray([1, 1, 1]), scaleX=matlabarray([ 0, 0, 0]), scalefacX=0, 
		#info=matlabarray([ 1, 1]))
		correctX = matlabarray([[0,     1,     0,     2],[0,     0,     1,     4]])  
		correctfX = matlabarray([1,     100,     101,   1])
		correctneval = 4
		correctxstatus =matlabarray([ 1,     1,     1,     1])
		correctsstatus =matlabarray([ 1,     1,     1,     1])
		correctdstatus = matlabarray([0,     0,     0,     0])
		
		#print "type correctX", type(correctX)
		#print "type X", type(X)
		self.assertEqual(X, correctX)
		self.assertEqual(sstatus, correctsstatus)
		self.assertEqual(dstatus, correctdstatus)
		
		print "Warning: Specified Test does not match function"
		#return unittest.skip("Thespecified Test does not match function...")
		self.assertEqual(xstatus, correctxstatus)
		self.assertEqual(fX, correctfX)
		self.assertEqual(neval, correctneval)

if __name__ == '__main__':
	unittest.main()