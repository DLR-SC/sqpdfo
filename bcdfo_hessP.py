#! /usr/bin/env python
from numpy import *
import helper

@helper.convertingDecorator
def  bcdfo_hessP_( P, x, xbase, scale, shift_Y ):
	return bcdfo_hessP( P, x, xbase, scale, shift_Y )

def  bcdfo_hessP( P, x, xbase, scale, shift_Y ):
###############################################################################

#  Computes the Hessian of the polynomial P at x, where P is represented by
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

#  INPUT:

#  P : a row vector containing the coefficient of the polynomial
#  x : the vector at which the Hessian of P is to be evaluated (not really
#      needed except for finding the space dimension, since we consider at most
#      quadratic polynomials).
#  xbase   : the current base point
#  scale   : the current interpolation set scaling
#  shift_Y : 0 if no shift in interpolation points, 1 otherwise

#  OUTPUT:

#  H : the Hessian matrix of P at x

#  PROGRAMMING: Ph. Toint,S. Gratton and A. Troeltzsch, April 2009.
#               (This version 22 VI 2009)

#  DEPENDENCIES: -

#  TEST:

#  Y = [ 1 2 1 3 3 1 ; 1 2 2 1 2 3 ];
#  [QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y, 0, 0 );
#  model = ( QZ * ( RZ' \ [1 2 3 4 5 6 ]' ) )';
#  bcdfo_hessP( model, [0;0],xbase, scale, 0 )
#  
#  ans =
#
#    4.0000   -0.5000
#   -0.5000    1.0000
#
#  the same result is obtained by
#
#  [QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y, 0, 1 );
#  model = ( QZ * ( RZ' \ [1 2 3 4 5 6 ]' ) )';
#  bcdfo_hessP( model, [0;0],xbase, scale, 1 )

#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.

###############################################################################

    n     = len( x )
    p1    = len( P[0] )
    nquad = p1 - n - 1

    P2 = copy(P)

    if ( shift_Y ):
        P2 = P2 * scale.transpose()[0]

    if( nquad > 0 ):

        # diagonal

        ndiag = min( nquad, n )
        H1 = P2[0][n+1:n+1+ndiag]
        H1 = append(H1,zeros(( 1,n-ndiag ))[0], axis=0)

        H     = diag( H1 )
        nquad = nquad - ndiag

        # subdiagonals

        if ( nquad > 0 ):
            k = 2*n+1;
            for i in range(1,n):
                nsd = min( n-i, nquad);
                if ( nsd > 0 ):
                    for j in range(0,nsd):
                        H[ i+j, j ] = P2[0][ k+j ] 
                        H[ j, i+j ] = P2[0][ k+j ] 

                    k = k + nsd;
                    nquad = nquad - nsd;

                if ( nquad == 0 ):
                    break;

    else:
        H = zeros(( n, n ))

    return H
