# -*- coding: utf-8 -*-
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from copy import copy    
from numpy import diag, zeros 
    
def bcdfo_hessP_(P=None,x=None,xbase=None,scale=None,shift_Y=None,*args,**kwargs):
    
    """
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

    """
#    varargin = cellarray(args)
#    nargin = 5-[P,x,xbase,scale,shift_Y].count(None)+len(args)

    n=length_(x)
    p1=length_(P)
    nquad=p1 - n - 1
    if (shift_Y):
        P=P*scale.T
        x=x - xbase
    P=P.reshape(-1)
#     diagonal
    if (nquad > 0):
        ndiag=min_(nquad,n)
        H=diag(concatenate_([P[n + 1:n + 1 + ndiag],zeros(n - ndiag)],axis=1))
        nquad=nquad - ndiag
#        subdiagonals
        if (nquad > 0):
            k=2 * n + 1
            for i in range(0,n - 1):
                nsd=min_(n - i - 1,nquad)
                if (nsd > 0):
                    for j in range(0,nsd):
                        H[i + j+1,j]=P[k + j]
                        H[j,i + j+1]=P[k + j]
                    k=k + nsd
                    nquad=nquad - nsd
                if (nquad == 0):
                    break
    else:
        H=zeros_(n,n)
    return H