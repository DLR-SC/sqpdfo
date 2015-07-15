# -*- coding: utf-8 -*-
"""
Created on Thu Nov 20 14:02:48 2014

@author: jaco_da
"""

# -*- coding: utf-8 -*-
"""
Created on Mon Nov 17 14:08:04 2014

@author: jaco_da
"""
import sys
import numpy as np
from runtime import *
#from bcdfo_build_QR_of_Y import *
from bcdfo_solve_TR_MS_bc import *




#print bcdfo_solve_TR_MS_bc( np.array([2, 3]).T, [4 6; 6 5], [-10; -10], [10; 10], 1.0, 0.001, 1 )
#  should give
#    0.5153
#   -0.8575
#
#print bcdfo_solve_TR_MS_bc( [2; 3], [4 6; 6 5], [-0.1; -0.1], [10; 10], 1.0, 0.001, 1 )
#  should give
#   -0.1
#   -0.1
#
#bcdfo_solve_TR_MS_bc( gx, H, lb, ub, Delta, eps_D, stratLam ):
#gx = matlabarray(np.array([2, 3])).T

#n = np.array([5,6,7,8])

#a = matlabarray([1,2,3,4])
#a = matlabarray([1,2,3,4]).T
#a = matlabarray([[1,2,3,4],[5,6,7,8]])

#print "a 1 1", a[1,2]

#dim = len(a.shape)
#print "a.shape", a.shape
#print "dim", dim

#print "b.shape", b.shape
#--------------------------------------------
#if a.shape[0] == 1:
	#b = np.zeros(a.shape)[0]
	#for i in range(len(b)):
	#	b[i] = a[i+1]
#else:
#	b = np.zeros(a.shape)#[0]
#	#print "b", b
#	#print "b.shape", b.shape
#	for i in range(b.shape[0]):
#		#print "i", i
#		for j in range(b.shape[1]):
#			#print "j", j
#			b[i][j] = a[i+1,j+1]#if a.shape[]:
#----------------------------------------------			
	#pass
#else:
#	raise NotImplemented


#print "a", a
#print "type(a)", type(a)

#print "b", b
#print "type(b)", type(b)

#print "a.shape", a.shape
#print "size_(a)", size_(a)#.size
#b = np.zeros(size_(a))
#print "b", b
#print "type(b)", type(b)

#aravel = a.ravel()
#print "aravel", aravel
#print "type(aravel)", type(aravel)


#b = a.__class__
#b = super(type(a), a)

#print "a", a
#print "type(a)", type(a)
#print "a.original", a.original
#print "type(a.original)", type(a.original)
#a.original[0][0] = 42
#print "a", a
#print "a.original", a.original

#print "n",n
#print "type(n)", n
#print "n.__class__.__dict__", n.__class__.__dict__, "\n"
#print "a", a
#print "type(a)", type(a)
#print "b", b
##print "b dict", b.__dict__, "\n"
#print "type(b)", type(b)

#a.__class__ = np.ndarray
#print "changed a", a
#print "changed type(a)", type(a)
#print "b[0]", b[0]
#print "b[2]", b[2]
#print "b[4]", b[4]

#sys.exit(0)

gx = matlabarray([2, 3]).T
print "gx", gx
#H = matlabarray(np.array([[4, 6],[6, 5]]))
H = matlabarray([[4, 6],[6, 5]])
print "H", H
#lb = np.array([-10, -10]).T
lb = matlabarray([-10, -10]).T
print "lb", lb
#ub = np.array([0, 0]).T
ub = matlabarray([0, 0]).T
print "ub", ub
#bcdfo_solve_TR_MS_bc( [2; 3], [4 6; 6 5], [-10; -10], [0; 0], 1.0, 0.001, 1 )
print bcdfo_solve_TR_MS_bc( gx, H, lb, ub, matlabarray([1.0]), matlabarray([0.001]), matlabarray([1]) )
#  should give
#    0
#   -0.6





#Y = np.array([ [1, 2, 1, 3, 3, 1], [1, 2, 2, 1, 1.01, 3 ]])
#QZ, RZ, xbase, scale = bcdfo_build_QR_of_Y( Y, 0, 1, 1, 1, 1e15 )
#print "QZ", QZ
#print "RZ", RZ
#print "xbase", xbase		
#print "scale", scale
#model = ( QZ * ( np.linalg.solve(RZ.T, np.array([1, 2, 3, 4, 5, 6 ]).T) ) ).T
#print "model", model
##print "scale[2]", scale[2]
#print "bcdfo_evalZ( (np.array([1,3]).T-xbase)*scale[2],6)", bcdfo_evalZ( (np.array([1,3]).T-xbase)*scale[2],6)
#print np.dot(model, bcdfo_evalZ( (np.array([1,3]).T-xbase)*scale[2],6))

#  should give 6.0
#bcdfo_build_QR_of_Y( Y, whichmodel, shift_Y, Delta, normgx, kappa_ill )
###############################################################################
#
#  Computes the QR factorization of the (possibly shifted) matrix containing
#  the polynomial expansion of the interpolation points. If shifting is
#  required, (i.e. if shift_Y is true) the matrix being factorized has columns
#  containing (Y(:,j)-Y(:,1))/ scale, where scale is the max( norm(Y(:,j)-Y(:,1)).
#
#  INPUT:
#
#  Y           : a matrix whose columns contain the current interpolation points
#  whichmodel  : kind of model to build
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#  Delta       : trust-region radius
#  normgx      : infinity norm of the projected gradient
#  kappa_ill   : threshold to declare a system matrix as ill-conditioned
#
#  OUTPUT:
#
#  QZ, RZ      : the QR factors of the (possibly shifted) matrix containing
#                the polynomial expansion of the interpolation points,
#  xbase       : the base point,
#  scale       : the model diagonal scaling.
#
#  PROGRAMMING: A. Troeltzsch, Ph. Toint, S. Gratton, 2009-2011.
#               (This version 14 I 2011)
#
#  DEPENDENCIES: bcdfo_evalZ
#
#  TEST:
#  Y = [ 1 2 1 3 3 1 ; 1 2 2 1 1.01 3 ];
#  [QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y, 0, 1, 1, 1, 1e15 );
#  model = ( QZ * ( RZ' \ [1 2 3 4 5 6 ]' ) )';
#  model * bcdfo_evalZ( ([1;3]-xbase)*scale(2),6)
#  should give 6.0
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
################################################################################







#arr = np.array([[1,2,3,4], [5,6,7,8], [9,10,11,12]])
#print "arr\n", arr

#print np.delete(arr, 1,1)

#detele last row
#arr1 = np.delete(arr, 1, 0)

#arr1 = np.delete(arr, 0, 1)
#print "arr1\n", arr1

#arr1 = np.delete(arr.T, 1, 0).T
#print "arr1\n", arr1

#print "arr1\n", arr1

#arr2 = np.delete(arr, np.s_[::2], 1)
#print "arr2\n", arr2
#b = np.s_[::2]
#print "nps", b
#a = slice(None, None, 2)
#print "a == b", a == b
#print "[::2]\n", arr[::2]
#print "arr[slice]", arr[b]

#print "arr 34\n", arr[2][:]

#arr3 = np.delete(arr, [1,3,5], None)
#print "arr3\n", arr3


