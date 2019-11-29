# -*- coding: utf-8 -*-

from sqpdfo.runtime import *
from sqpdfo.bcdfo_evalZ import bcdfo_evalZ_

def bcdfo_evalP_(P=None,x=None,xbase=None,scale=None,shift_Y=None,*args,**kwargs):
    """
#
#  Computes the value of the polynomial P at x, where P is represented by
#  the row vector containing its coefficients for the successive monomials.
#  More specifically, these values are:
#  P(1)        : constant coefficient,
#  P(2:n+1)    : coefficients for the linear terms in x(1)... x(n),
#  P(n+2:2n+2) : coefficients for the squared terms in x(1)^2 ... x(n)^2
#  P(2n+3,3n+2): coefficients for the quadratic terms of the first subdiagonal: 
#                in x(1)*x(2) ... x(n-1)*x(n)
#  (3n+3,4n+1) : coefficients for the quadratic terms of the second subdiagonal: 
#                in x(1)*x(3) ... x(n-2)*x(n)
#  etc.
#
#  INPUT:
#
#  model       : the current model
#  x           : the point at which the model must be evaluated
#  xbase       : the current base point
#  whichmodel  : kind of model to build
#  scale       : the current model scaling
#  shift_Y     : 0 if no shift in interpolation points, 1 otherwise
#
#  OUTPUT:
#
#  value       : the value of the model at x.
#
#  PROGRAMMING: Ph. Toint, S. Gratton, April 2009. (This version 22 VI 2009)
#
#  DEPENDENCIES: bcdfo_evalZ
#
#  TEST:
#  Y = [ 1 2 1 3 3 1 ; 1 2 2 1 1.01 3 ]; 
#  [QZ, RZ, xbase, scale ] = bcdfo_build_QR_of_Y( Y, 0, 1 );
#  model = ( QZ * ( RZ' \ [1 2 3 4 5 6 ]' ) )';
#  bcdfo_evalP( model, [1;3], xbase, scale, 1 )
#  should give 
#     6.0
#
#  CONDITIONS OF USE: Use at your own risk! No guarantee of any kind given.
#
    """
#    varargin = cellarray(args)
#    nargin = 5-[P,x,xbase,scale,shift_Y].count(None)+len(args)

    if (shift_Y):
        value=P.dot(  bcdfo_evalZ_((x - xbase) * scale[1],size_(P,2)))
    else:
        value=P.dot(  bcdfo_evalZ_(x,size_(P,2)))
    return value