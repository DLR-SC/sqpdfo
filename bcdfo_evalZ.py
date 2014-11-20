#! /usr/bin/env python
from matrixalgebra import *
from numpy import *

def bcdfo_evalZ( X, q ):
###############################################################################
#
#  Compute the matrix Z(X), where X is a matrix whose columns contains points
#  of the underlying space.  The vector Z(x) is, for a vector x, given by the
#  sequence of the values of the (at most quadratic) monomials taken at x.  
#  More specifically, these values are:
#  Z(x)(1)        : 1,
#  Z(x)(2:n+1)    : the linear terms x(1)... x(n),
#  Z(x)(n+2:2n+2) : the diagonal terms of the quadratic: x(1)^2 ... x(n)^2
#  Z(x)(2n+3,3n+2): the first subdiagonal of the quadratic: 
#                    x(1)*x(2) ... x(n-1)*x(n)
#  Z(x)(3n+3,4n+1): the second subdiagonal of the quadratic: 
#                    x(1)*x(3) ... x(n-2)*x(n)
#  etc.
#
#  INPUTS:
#
#  X          : the matrix whose columns contains the points at which the monomials
#               should be evaluated.
#  q          : the number of monomials considered (q <= (n+1)*(n+2)/2)
#
#  OUTPUT:
#
#  Z  : the matrix Z(X), of size q x m.
#
#  PROGRAMMING: Ph. Toint, January 2009.
#               (This version 12 IV 2010)
#
#  DEPENDENCIES: -
#
#  TEST:
#  bcdfo_evalZ( [ 0 1 0 2 1 0 ; 0 0 1 0 0.01 2 ], 6 )
#  should give
# ans =
#
#    1.0000    1.0000    1.0000    1.0000    1.0000    1.0000
#         0    1.0000         0    2.0000    1.0000         0
#         0         0    1.0000         0    0.0100    2.0000
#         0    0.5000         0    2.0000    0.5000         0
#         0         0    0.5000         0    0.0001    2.0000
#         0         0         0         0    0.0100         0
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
###############################################################################

   [ n, m ] = shape( X );             # [ dimension of the space, number of points in X ]
   nlin     = min( n+1, q );         # number of constant and linear terms
   nquad    = max( 0, q-nlin );      # number of quadratic terms
   nlin     = nlin - 1;              # number of linear terms
   Z        = zeros(( q, m ), float );
   
   if ( q == 1 ):
      Z = ones(( 1, m ), float);                       # constant terms
   elif ( q <= n+1 ):
      Z = ones(( 1, m ));
      Z = append(Z, X[0:nlin,0:m], axis=0);    # constant and linear
   else:
      ndiag = min( n, nquad );
      Z = ones(( 1, m ), float)
      Z = append( Z, X[0:n,0:m], axis=0 ); # same + diagonal
      Z = append( Z, 0.5*X[0:ndiag,0:m]**2, axis=0 );
      nquad = nquad - ndiag;
      if ( nquad > 0 ):
         for k in range(0,n-1):                        # the (i+1)-th subdiagonal
            nsd = min( n-k-1, nquad );

            if ( nsd > 0 ):            
               Z = append( Z, X[k+1:k+1+nsd,0:m]*X[0:nsd,0:m], axis=0 )
               nquad = nquad - nsd;
         
            if ( nquad == 0 ):
               break;

   return Z


