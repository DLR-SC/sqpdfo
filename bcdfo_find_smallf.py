#!/usr/local/bin/python
from numpy import *
from math import *
from bcdfo_swap_in_Y import *

def bcdfo_find_smallf(c, QZ, RZ, Y, fY, ind_Y, i_xbest, cur_degree,
      indfree, x, xl, xu, fx, dstatus, whichmodel, scale, shift_Y, Delta, 
      normgx, kappa_ill):
#############################################################################
#
#  Subroutine finds the smallest value in fY, which are the associated function
#  values of the set Y. There are only points inside the bounds considered
#  and those which are and no dummy points!
#  Exchanges the best point in the current interpolation set Y if a smaller
#  function value was found.
#
#  INPUT:
#
#  c           : contains a bunch of constants
#  QZ, RZ      : the QR factors of the (possibly shifted) matrix containing
#                the polynomial expansion of the interpolation points
#  Y           : interpolation set
#  fY          : function values associated to the current interpolation points
#  ind_Y       : indices of the points in Y out of X, the set of all points
#  i_xbest     : index of the best point
#  cur_degree  : number of interpolation points
#  indfree     : number of free variables
#  x           : best point
#  xl, xu      : lower/upper bounds on the variables
#  fx          : best function value
#  dstatus     : status vector of dummy points in X
#  whichmodel  : kind of model to build
#  scale       : model diagonal scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#  Delta       : trust-region radius
#  normgx      : infinity norm of the projected gradient
#  kappa_ill   : threshold to declare a system matrix as ill-conditioned
#
#  OUTPUT:
#
#  (possibly) updated INPUT values
#
#  PROGRAMMING: A. Troeltzsch, August 2010.
#
#  DEPENDENCIES: bcdfo_swap_in_Y
#
#  TEST:
#  Y = [ 0 1 0 2 0 ; 0 0 1 0 2 ]; x = [ 0; 0];
#  [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y(  Y, 0, 0, 1, 1, 1e15 );
#  fY = [ 1.5 2 3 4 1 ]; fx = 1.5;
#  ind_Y = [ 10 11 8 7 12 ]; i_xbest = 10;
#  cur_degree = 5;
#  indfree = [ 1 2 ];
#  xl = [ -1e20; -1e20 ]; xu = [ 1e20; 1e20 ];
#  dstatus = [ 0; 0; 0; 0; 0 ];
#  c.dummy = 1;
#  [x, fx, QZ, RZ, Y, fY, ind_Y, i_xbest, scale] = ...
#  bcdfo_find_smallf(c, QZ, RZ, Y, fY, ind_Y, i_xbest, cur_degree, ...
#      indfree, x, xl, xu, fx, dstatus, 0, scale, 1, 1, 1, 1e15 )
#  which gives:
#  x =
#     0
#     2
#  fx =
#     1
#  QZ =
#   1.00000000000000                  0                  0                  0                  0
#                  0  -0.42519520276219  -0.82880426938540  -0.10330005312161  -0.34874291623146
#                  0   0.85039040552437  -0.49972022124708   0.11662909223407   0.11624763874382
#                  0  -0.07516460280028  -0.14651327978970  -0.73515420756377   0.65759594922143
#                  0  -0.30065841120113  -0.20468767029443   0.65975377601877   0.65759594922143
#  RZ =
#   1.00000000000000   1.00000000000000   1.00000000000000   1.00000000000000   1.00000000000000
#                  0  -0.83150841847813  -0.31944956190120  -0.99593098710375  -0.67648142520255
#                  0                  0   0.16388479917652  -0.32049779953996   0.30218363956625
#                  0                  0                  0  -0.17436349794782   0.08246922200235
#                  0                  0                  0                  0   0.08219949365268
#  Y =
#     0     1     0     2     0
#     2     0     1     0     0
#  fY =
#   1.00000000000000   2.00000000000000   3.00000000000000   4.00000000000000   1.50000000000000
#  ind_Y =
#    12    11     8     7    10
#  i_xbest =
#    12
#  scale =
#   1.00000000000000
#   0.35355339059327
#   0.35355339059327
#   0.12500000000000
#   0.12500000000000
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
#############################################################################

   dummy_set = nonzero( dstatus == c.dummy() )[0]
   ind_insideBounds = array([],dtype=int)
   
   #  Select indices of all points inside the bounds and which are not a dummy point

   for i in range( 0, cur_degree ):
      if ( len(nonzero(Y[:,i] < xl[0,indfree])[0]) == 0 and len(nonzero(Y[:,i] > xu[0,indfree])[0] ) == 0 and (len(nonzero( dummy_set == ind_Y[i] )[0]) == 0)):
         ind_insideBounds = append( ind_insideBounds, i )
      else:
         ind_insideBounds = append( ind_insideBounds, 0 )
   
   #  Find the smallest function value among them

   fmin  = min( fY[ind_insideBounds] )

   #  Exchange best point only if the new smallest fmin is strictly smaller than fx
   #  (to avoid a change of the center of the interpolation at equal function values)

   if ( fmin < fx):
      imin = argmin(fY)
      [ QZ, RZ, Y, ind_Y, fY, x, scale ] = bcdfo_swap_in_Y( 0, imin, QZ, RZ, Y,
           ind_Y, fY, x, whichmodel, scale, shift_Y, Delta, normgx, kappa_ill );

      fx       = fY[0]
      i_xbest  = ind_Y[0]
      if ( not shift_Y ):
         x = Y[:,0]

   return x, fx, QZ, RZ, Y, fY, ind_Y, i_xbest, scale
   
