# -*- coding: utf-8 -*-
"""
Created on Fri Nov 14 18:53:38 2014
%  TEST:
%
%  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; ind_Y = [1 2 3 4 5 6];
%  fY = [2 5 1 3 2 6];
%  [QZ,RZ,xbase,scale] = bcdfo_build_QR_of_Y( Y, 0, 0, 1, 1, 1e15 );
%  Z = QZ*RZ;
%  [ QZ, RZ, Y, xbase, scale ] = bcdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
%     fY, xbase, 0, scale, 0, 1, 1, 1e15 )
%  [ QZ, RZ, Y, xbase, scale ] = bcdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
%     fY, xbase, 0, scale, 0, 1, 1, 1e15 )
%  norm( Z - QZ*RZ)
%  should give something very small.
%
%  The same holds for the scaled version:
%  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; ind_Y = [1 2 3 4 5 6];
%  fY = [2 5 1 3 2 6];
%  [QZ,RZ,xbase,scale] = bcdfo_build_QR_of_Y( Y, 0, 1, 1, 1, 1e15 );
%  Z = QZ*RZ;
%  [ QZ, RZ, Y, xbase, scale ] = bcdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
%     fY, xbase, 0, scale, 1, 1, 1, 1e15 )
%  [ QZ, RZ, Y, xbase, scale ] = bcdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
%     fY, xbase, 0, scale, 1, 1, 1, 1e15 )
%  norm( Z - QZ*RZ)
@author: jaco_da
"""

import unittest
from ecdfo_swap_in_Y import *
from runtime import *

from bcdfo_build_QR_of_Y import *
#import numpy as np
#import helper

class Test_ecdfo_swap_in_Y(unittest.TestCase):

	def setUp(self):
		self.Y = matlabarray([[ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]])
		self.ind_Y = matlabarray([1, 2, 3, 4, 5, 6])
		self.fY = matlabarray([2, 5, 1, 3, 2, 6])
		#self.options = helper.dummyOptions()
		#self.values = helper.dummyValues()

		pass

	def test_ecdfo_swap_in_Y(self):

		#self.Y = np.array([[ 0, 1, 0, 2, 1, 0], [0, 0, 1, 0, 0.01, 2 ]])
		QZ,RZ,xbase,scale = bcdfo_build_QR_of_Y_( self.Y, 0, 0, 1, 1, 1e15 )
		Z = QZ*RZ
		
		print "QZ", QZ
		print "type(QZ)", type(QZ)
		print "RZ", RZ
		print "type(RZ)", type(RZ)
		print "Z", Z
		print "xbase", xbase
		print "type(xbase)", type(xbase)		
		#self.fY, matlabarray(xbase), 0, matlabarray(scale), 0, 1, 1, 1e15 )

		print "scale", scale
		print "tyoe(scale)", type(scale)
		
		#QZ, RZ, Y, xbase, scale = ecdfo_swap_in_Y_( 1, 3, matlabarray(QZ), matlabarray(RZ), self.Y, self.ind_Y,
		ceY = self.Y
		ciY = self.Y
		QZ,RZ,Y,ind_Y,fY,ciY,ceY,xbase,scale = ecdfo_swap_in_Y_( 1, 3, QZ, RZ, self.Y, self.ind_Y,
		self.fY, ciY, ceY, xbase, 0, scale, 0, 1, 1, 1e15 )
		ret = norm_( Z - QZ*RZ)		
		self.assertTrue(ret < 2.0)
		
		#ecdfo_swap_in_Y_(i=None,j=None,QZ=None,RZ=None,Y=None,ind_Y=None,fY=None,
		#ciY=None,ceY=None,xbase=None,whichmodel=None,scale=None,shift_Y=None,Delta=None,
		#normgx=None,kappa_ill=None)

if __name__ == '__main__':
	unittest.main()