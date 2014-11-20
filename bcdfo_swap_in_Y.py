#!/usr/local/bin/python
from numpy import *
from bcdfo_build_QR_of_Y import *

def bcdfo_swap_in_Y( i, j, QZ, RZ, Y, ind_Y, fY, xbase, whichmodel, scale, 
      shift_Y, Delta, normgx, kappa_ill ):

###############################################################################
#
#  Swaps the position of interpolation points i and j in Y, and updates the
#  factorization of Z(Y) accordingly.
#
#  INPUT:
#
#  i           : the position in which Y(:,j) should be inserted
#  j           : the position in which Y(:,i) should be inserted
#  QZ, RZ      : the QR decomposition of Z(Y)
#  Y           : a matrix whose columns contain the current interpolation points
#  ind_Y		   : set of indices which points are in Y out of set X
#  fY				: corresponding function values
#  xbase       : the current base point
#  whichmodel  : kind of model/Lagrange polynomial to compute
#  scale       : the current interpolation set scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#  Delta       : trust-region radius
#  normgx      : infinity norm of the projected gradient
#  kappa_ill   : threshold to declare a system matrix as ill-conditioned
#
#  OUTPUT:
#
#  QZ          : the updated Q matrix of the QR decomposition of Z(Y)
#  RZ          : the updated R matrix of the QR decomposition of Z(Y)
#  Y           : the updated set of interpolation points (column-wise)
#  xbase       : the base point after the exchange
#  scale       : the interpolation set scaling after the exchange
#
#  PROGRAMMING: Ph. Toint and A. Troeltzsch, January 2009.
#               ( This version 13 IX 2010 )
#
#  DEPENDENCIES: bcdfo_build_QR_of_Y
#
#  TEST:
#
#  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; ind_Y = [1 2 3 4 5 6];
#  fY = [2 5 1 3 2 6];
#  [QZ,RZ,xbase,scale] = bcdfo_build_QR_of_Y( Y, 0, 0, 1, 1, 1e15 );
#  Z = QZ*RZ;
#  [ QZ, RZ, Y, xbase, scale ] = bcdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
#     fY, xbase, 0, scale, 0, 1, 1, 1e15 )
#  [ QZ, RZ, Y, xbase, scale ] = bcdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
#     fY, xbase, 0, scale, 0, 1, 1, 1e15 )
#  norm( Z - QZ*RZ)
#  should give something very small.
#
#  The same holds for the scaled version:
#  Y = [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ]; ind_Y = [1 2 3 4 5 6];
#  fY = [2 5 1 3 2 6];
#  [QZ,RZ,xbase,scale] = bcdfo_build_QR_of_Y( Y, 0, 1, 1, 1, 1e15 );
#  Z = QZ*RZ;
#  [ QZ, RZ, Y, xbase, scale ] = bcdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
#     fY, xbase, 0, scale, 1, 1, 1, 1e15 )
#  [ QZ, RZ, Y, xbase, scale ] = bcdfo_swap_in_Y( 1, 3, QZ, RZ, Y, ind_Y, ...
#     fY, xbase, 0, scale, 1, 1, 1, 1e15 )
#  norm( Z - QZ*RZ)
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
################################################################################

   #  Ensure that ii is smaller than jj.

   if ( i > j ):
      ii = j;
      jj = i;
   elif ( i < j ):
      ii = i;
      jj = j;
   else:
      return QZ, RZ, Y, ind_Y, fY, xbase, scale

   #  Permute the columns of Y, the indices in ind_Y and the values in fY.

   y         = copy(Y[:,ii])
   Y[:,ii]   = copy(Y[:,jj])
   Y[:,jj]   = y;

   ind       = copy(ind_Y[ii])
   ind_Y[ii] = copy(ind_Y[jj])
   ind_Y[jj] = ind

   f         = copy(fY[ii])
   fY[ii]    = copy(fY[jj])
   fY[jj]    = f;

   #  a new factorization must be computed

   [ QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y, whichmodel, shift_Y, 
      Delta, normgx, kappa_ill )

   return QZ, RZ, Y, ind_Y, fY, xbase, scale
