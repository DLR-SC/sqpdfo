#!/usr/local/bin/python
from numpy import *

def bcdfo_augmX_evalf(f, y, m, X, fX, nfix, xfix, indfix, indfree, fxmax, 
      neval, maxeval, verbose, show_errg, xstatus, xstatus_val, sstatus,
      dstatus, scaleX, scalefacX, compute_g, gmat ):
#################################################################################
#
#  function bcdfo_augmX_evalf
#
#  A new point y is added to the set of all points X. The objective function is
#  evaluated at y and the function value is added to fX (or the model value
#  of the dummy point y is replaced by its real function value in fX).
#
#  INPUTS:
#
#  f           : function handle of the function to evaluate
#  y           : the point which is added to X
#  m           : index where to place y in X
#  X           : set of all points
#  fX          : function (or model) values associated to X
#  nfix        : number of fixed variables
#  xfix        : contains the values of the fixed variables
#  indfix      : fixed indices of x
#  indfree     : free indices of x
#  fxmax       : maximal allowed function value
#                (to avoid infinity in the calculation)
#  neval       : actual number of function evaluations
#  maxeval     : maximal allowed number of function evaluations
#  verbose     : level of printing (from calling routine)
#  show_errg   : gradient error estimate is computed or not
#  xstatus     : status vector if points contained in current interpolation set
#  xstatus_val : value, if new point y will be contained in interpolation set or not
#  sstatus     : status vector if points in lie in current subspace
#  dstatus     : status vector if points are dummy points
#  scaleX      : scaling of variables is applied
#  scalefacX   : scaling factors
#  compute_g   : compute function AND gradient at current x
#  gmat        : a matrix containing all computed gradients
#
#  OUTPUTS:
#
#  X           : updated set of all points
#  fX          : updated set of function values associated with X
#  neval       : updated number of function evaluations
#  xstatus     : updated status vector if points contained in current interpolation set
#  sstatus     : updated status vector if points in lie in current subspace
#  dstatus     : updated status vector if points are dummy points
#  gmat        : a matrix containing all computed gradients
#  msg         : error message if any
#
#  PROGRAMMING: A. Troeltzsch, August 2010.
#              (This version: 3 XI 2011)
#
#  DEPENDENCIES: -
#
#  TEST:
#  To append point y on set X:
#  X = [ 1; 0 ];
#  fX = [ 100 ];
#  y = [ 2; 4 ];
#  [X, fX, neval, xstatus, sstatus, dstatus, gmat, msg] = bcdfo_augmX_evalf( @banana, ...
#     y, 2, X, fX, 0, [0;0], [], [1 2], 1e25, 1, 100, 1, ...
#     0, [1], 1, [1], [0], 0, [1,1], 0, [] )
#  where
#     function fx = banana( x )
#     fx  = 100 * ( x(2) - x(1)^2 ) ^2 + (1-x(1))^2;
#  which gives
#  X =
#     1     2
#     0     4
#  fX =
#    100    1
#  neval =
#     2
#  xstatus =
#     1     1
#  sstatus =
#     1     1
#  dstatus =
#     0     0
#  gmat =
#     []
#  msg =
#     empty message
#
#  To compute and store also the gradient in gmat at a point y:
#  y = [ 0; 1 ];
#  [X, fX, neval, xstatus, sstatus, dstatus, gmat, msg] = bcdfo_augmX_evalf( @banana_g, ...
#     y, 3, X, fX, 0, [0;0], [], [1 2], 1e25, neval, 100, 1, ...
#     0, xstatus, 1, sstatus, dstatus, 0, [1,1], 1, [400 2; -200 0] )
#  where
#     function [fx,gx] = banana_g( x )
#     fx  = 100 * ( x(2) - x(1)^2 )^2 + (1-x(1))^2;
#     if (nargout > 1) # gradient required
#        gx = [-400*(x(2)-x(1)^2)*x(1)-2*(1-x(1)); 200*(x(2)-x(1)^2)];
#     end
#  which gives
#  X =
#     1     2     0
#     0     4     1
#  fX =
#   100     1   101
#  neval =
#     4
#  xstatus =
#     1     1     1
#  sstatus =
#     1     1     1
#  dstatus =
#     0     0     0
#  gmat =
#   400     2    -2
#  -200     0   200
#  msg =
#     empty message
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
###############################################################################
   
   msg     = 'unexpected error'
   full_n  = len( xfix )
   I       = eye( full_n )
   
   #  Check dimension of y
   
   if ( ndim(y) > 1 ):
      if ( y.shape[1] == 1 and y.shape[0] > 1 ):
         y = y[:,0] 
      else:
         y = y[0] 
         
   #  Set status of the new point
   
   xstatus = append( xstatus, xstatus_val )
   sstatus = append( sstatus, 1 )         # point is contained in current subspace
   dstatus = append( dstatus, 0 )         # point is not a dummy point

   #  Augment X with full-dimensional y, scale if user-defined and evaluate f at y

   if ( nfix > 0 ):

      yfull  = dot( I[:,indfix], xfix[indfix] ) + dot( I[:,indfree], y )

      X = append( X, yfull, axis=1 )

      if ( scaleX ):
         yfull = yfull / scalefacX;

      if ( compute_g == 1 ):
         [fvalue, g] = f( yfull )
         gmat = append( gmat, g, axis=1 )
      else:
         fvalue = f( yfull )

   else:

      if len( X[0] ) == 0:
         X = y.reshape(full_n,1)
      else:
         X = c_[ X, y ]

      if ( scaleX ):
         y = y / scalefacX;

      if ( compute_g == 1 ):
         [fvalue, g] = f( y )
         gmat = append( gmat, g, axis=1 )
      else:
         fvalue = f( y )

   #  Augment fX with new function value

   fX = append( fX, min( fxmax, real( fvalue ) ) )

   #  Update evaluation counter

   if ( compute_g == 1 ):
      neval = neval + 2
   else:
      neval = neval + 1

   #  Check for maximum number of evaluations

   if ( neval >= maxeval ):
      if ( verbose >= 1 ):
         if ( show_errg ):
            disp( ' !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! MAX. EVAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!')
         else:
            disp( ' !!!!!!!!!!!!!!!!!!!!!!!!!!! MAX. EVAL !!!!!!!!!!!!!!!!!!!!!!!!!')

      msg  = 'Error: Maximum number of '+str( maxeval )+' function evaluations reached.'
      return X, fX, neval, xstatus, sstatus, dstatus, gmat, msg

   return X, fX, neval, xstatus, sstatus, dstatus, gmat, msg
   
