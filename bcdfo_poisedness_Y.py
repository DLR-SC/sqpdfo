#!/usr/local/bin/python
from numpy import *
from bcdfo_find_new_yj import *
import helper
from runtime import matlabarray

def bcdfo_poisedness_Y_( QZ, RZ, Y, eps_L, xbase, lSolver, whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y, nargout=None ):
	QZ = helper.convert(QZ)
	RZ = helper.convert(RZ)
	Y = helper.convert(Y)
	eps_L = helper.convert(eps_L)
	xbase = helper.convert(xbase)
	lSolver = helper.convert(lSolver)
	whichmodel = helper.convert(whichmodel)
	hardcons = helper.convert(hardcons)
	xl = helper.convert(xl)
	xu = helper.convert(xu)
	
	#print "indfree before = ", indfree	
	#print "len indfree", len(indfree)
	#print "ones", ones(len(indfree))
	indfree = helper.convert(indfree - 1)
	#print "indfree after = ", indfree
	
	stratLam = helper.convert(stratLam)
	scale = helper.convert(scale)
	shift_Y = helper.convert(shift_Y)
	
	
	lambd, Y_radius = bcdfo_poisedness_Y( QZ, RZ, Y, eps_L, xbase, lSolver, whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y )
	
	lambd = matlabarray(lambd)
	Y_radius = matlabarray(Y_radius)
	
	return lambd, Y_radius
	
	

def bcdfo_poisedness_Y( QZ, RZ, Y, eps_L, xbase, lSolver, whichmodel, hardcons, xl, xu, indfree, stratLam, scale, shift_Y ):
   """
   Computes the poisedness of the interpolation set Y in a ball of radius
   Delta centered at Y(:,1), assuming that a QR factorization of the
   interpolation system is known (and given by QZ, RZ).
   Poisedness is defined here as the maximum aboslute value of
   the Lagrange polynomials taken over the ball and for all polynomials.

   INPUT:

   QZ          : the Q matrix of the QR decomposition of Z(Y)
   RZ          : the R matrix of the QR decomposition of Z(Y)
   Y           : a matrix whose columns contain the current interpolation points
   eps_L       : the relative accuracy on the trust-region constraint for
                maximization of the Lagrange polynomials
   xbase       : the current base point
   lSolver     : linear solver used for the minimization of the model
   whichmodel  : kind of model/Lagrange polynomial to compute
   scale       : the current interpolation set scaling
   shift_Y     : 0 if no shift in interpolation points, 1 otherwise

   OUTPUT:

   lambda      : the poisedness of the interpolation set Y
   Y_radius    : the poisedness radius for Y, that is the largest distance
                from Y(:,i) to Y(:,1) (the base point)

   PROGRAMMING: Ph. Toint and S. Gratton, April 2009. (This version 22 VI 2009)

   DEPENDENCIES: bcdfo_find_new_yj, bcdfo_find_new_yj_bc

   CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
   """

   lambd = 0
   [n, p1] = shape( Y )
   
   #  Compute the radius of the poisedness ball.

   Y_radius = 0;
   lst = range(1,p1)
   for j in lst:
      Y_radius = max( Y_radius, linalg.norm( Y[:,j] - Y[:,0] ) )

   #  Loop on all the possible replacements to find the best one.

   for j in lst:
      if ( hardcons == 1 ):
   
         [ y, improvement ] = bcdfo_find_new_yj_bc( QZ, RZ, Y, j, Y_radius, 
          eps_L, xbase, lSolver, whichmodel, xl, xu, indfree, 
          stratLam, scale, shift_Y );
               
      else:

         [ y, improvement ] = bcdfo_find_new_yj( QZ, RZ, Y, j, Y_radius, eps_L, 
         xbase, lSolver, whichmodel, scale, shift_Y );
        

      #  Remember the current polynomial value, index and replacement point if
      #  this is the best so far.
   
      lambd = max( improvement, lambd );

   return lambd, Y_radius

#----------------------------------

if __name__ == "__main__":
   from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y
   import numpy
   #print "Self-test of the function without shifting:"
   Y = numpy.array([ [ 0, 1, 0, 2, 1, 0 ], [ 0, 0, 1, 0, 0.01, 2] ])
   #print Y
   [QZ, RZ, xbase, scale] = bcdfo_build_QR_of_Y.bcdfo_build_QR_of_Y( Y, 0, 0 )
   [ lambd ,Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, 0.001, xbase, 1, 0, scale, 0 )
   #print lambd, Y_radius
#  should give:
#  lambd =
#  204.8586
#  Y_radius =
#     2

   #print "Self-test of the function with shifting:"
   Y = numpy.array([ [ 0, 1, 0, 2, 1, 0 ], [ 0, 0, 1, 0, 0.01, 2] ])
   #print Y
   [QZ, RZ, xbase, scale] = bcdfo_build_QR_of_Y( Y, 0, 1 )
   [ lambd ,Y_radius ] = bcdfo_poisedness_Y( QZ, RZ, Y, 0.001, xbase, 1, 0, scale, 1 )
   #print lambd, Y_radius
#  should give:
#  lambd =
#  204.8586
#  Y_radius =
#     2