# -*- coding: utf-8 -*-

from sqpdfo.runtime import *
    
def bcdfo_gradP_(P=None,x=None,xbase=None,scale=None,shift_Y=None,*args,**kwargs):
    """
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

    """
#    varargin = cellarray(args)
#    nargin = 5-[P,x,xbase,scale,shift_Y].count(None)+len(args)
#
    n=length_(x)
    p1=length_(P)
    if (shift_Y):
        P=P * scale.T
        x=x - xbase
    P=P.reshape(-1)
    x=x.reshape(-1)
#  First-order terms in the polynomial
    ng=min_(n,p1 - 1)
    g=zeros_(n,1).reshape(-1)
    g[0:ng]=P[1:ng + 1]
#  Second-order terms
    nquad=p1 - n - 1
    if (nquad > 0):
#     diagonal
        ndiag=min_(nquad,n)
        g[0:ndiag]=g[0:ndiag] + P[n + 1:n + 1 + ndiag] *(x[0:ndiag])
        nquad=nquad - ndiag
#         subdiagonals
        if (nquad > 0):
            k=2 * n + 1
            for i in range(0,n - 1):
                nsd=min_(n - i -1,nquad)
                if (nsd > 0):
                    g[i + 1:i + nsd+1]=g[i + 1:i + nsd+1] + P[k :k + nsd] *(x[0:nsd])
                    g[0:nsd]=g[0:nsd] + P[k :k + nsd] * (x[i + 1:i + nsd+1])
                    k=k + nsd
                    nquad=nquad - nsd
                if (nquad == 0):
                    break
    return g.reshape(-1,1)