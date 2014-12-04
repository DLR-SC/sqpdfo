#! /usr/bin/env python
from numpy import *
import helper

@helper.convertingDecorator
def bcdfo_gradP_( P, x, xbase, scale, shift_Y ):
	return bcdfo_gradP( P, x, xbase, scale, shift_Y )

def bcdfo_gradP( P, x, xbase, scale, shift_Y ):

###############################################################################

#  Computes the gradient of the polynomial P at x, where P is represented by
#  the row vector containing its coefficients for the successive monomials.
#  More specifically, these values are:
#  P(1)        : constant coefficient,
#  P(2:n+1)    : coefficients for the linear terms in x(1)... x(n),
#  P(n+2:2n+2) : coefficients for the squared terms in x(1)^2 ... x(n)^2
#  P(2n+3,3n+2): coefficients for the quadratic terms of the first subdiagonal:
#                in x(1)*x(2) ... x(n-1)*x(n)
#  (3n+3,4n+1): coefficients for the quadratic terms of the second subdiagonal:
#                in x(1)*x(3) ... x(n-2)*x(n)
#  etc.
#
#  INPUT:
#
#  P       : a row vector contains the coefficients of the polynomial
#  x       : the point at which the gradient must be evaluated
#  xbase   : the current base point
#  scale   : the current interpolation set scaling
#  shift_Y : 0 if no shift in interpolation points, 1 otherwise
#
#  OUTPUT:
#
#  g : the gradient of P at x
#
#  PROGRAMMING: Ph. Toint, S. Gratton and A. Troeltzsch, April 2009.
#               (This version 30 IV 2009)
#
#  DEPENDENCIES: -
#
#  TEST:
#
#  Y = [ 1 2 1 3 3 1 ; 1 2 2 1 2 3 ];
#  [QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y, 0, 0 );
#  model = ( QZ * ( RZ' \ [1 2 3 4 5 6 ]' ) )';
#  bcdfo_gradP( model, [0;0], xbase, scale, 0 )
#
#  ans =
#
#   -6.0000
#    1.0000
#
#  the same result is obtained by
#
#  Y = [ 1 2 1 3 3 1 ; 1 2 2 1 2 3 ];
#  [QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y, 0, 1 );
#  model = ( QZ * ( RZ' \ [1 2 3 4 5 6 ]' ) )';
#  bcdfo_gradP( model, [0;0], xbase, scale, 1 )

#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.

###############################################################################

   n    = len( x )
   p1   = len( P[0] )
   P_cp = copy(P)
   x_cp = copy(x.transpose())
   print "x\n", x			
   
   if ( shift_Y ):
      P_cp = P_cp * scale.transpose()[0]
      x_cp = x_cp - xbase.transpose();
   
   #  First-order terms in the polynomial

   ng         = min( n, p1-1 )
   g          = zeros(( n,1 ),float).transpose()
   g[0][0:ng] = P_cp[0][1:ng+1]

   #  Second-order terms

   nquad = p1 - n - 1
   
   if( nquad > 0 ):

       # diagonal

       ndiag        = min( nquad, n )
       print "g\n", g
       print "P_cp\n", P_cp
       print "x_cp\n", x_cp							
       g[0][0:ndiag] = g[0][0:ndiag] + P_cp[0][n+1:n+1+ndiag].transpose()*x_cp[0][0:ndiag]
       nquad        = nquad - ndiag

       # subdiagonals

       if ( nquad > 0 ):
          k = 2 * n + 1
          for i in range(0,n-1):
             nsd = min( n-i-1, nquad )
             
             if ( nsd > 0 ):
                g[0][i+1:i+nsd+1] = g[0][i+1:i+nsd+1] + P_cp[0][k:k+nsd] * x_cp[0][0:nsd]
                g[0][0:nsd]       = g[0][0:nsd]   + P_cp[0][k:k+nsd] * x_cp[0][i+1:i+nsd+1]
                k     = k + nsd
                nquad = nquad - nsd

             if ( nquad == 0 ):
                break

   return g
   
