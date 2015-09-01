# -*- coding: utf-8 -*-
#! /usr/bin/env python
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *


def bcdfo_evalZ_(X=None,q=None,*args,**kwargs):
    """
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
    """

#    varargin = cellarray(args)
#    nargin = 2-[X,q].count(None)+len(args)

    n,m=size_(X,nargout=2) # [ dimension of the space, number of points in X ]
    nlin=min_(n + 1,q) #  number of constant and linear terms
    nquad=max_(0,q - nlin) # number of quadratic terms
    nlin=nlin - 1 # number of linear terms
    Z=zeros_(q,m) 
    if (q == 1):
        Z=ones_(1,m) # constant terms
    elif (q <= n + 1):
        Z=concatenate_((ones_(1,m),X[0:nlin,0:m]),axis=0) #   constant and linear
    else:
        ndiag=min_(n,nquad)
        Z=concatenate_((ones_(1,m),X[0:n,0:m],0.5 * X[0:ndiag,0:m] ** 2), axis=0) # same + diagonal
        nquad=nquad - ndiag
        if (nquad > 0):
            for k in range(1,n):                          # the (i+1)-th subdiagonal
                nsd=min_(n - k,nquad)
                if (nsd > 0):
                    Z=concatenate_((Z,X[k:k + nsd,0:m]*(X[0:nsd,0:m])),axis=0)
                    nquad=nquad - nsd
                if (nquad == 0):
                    break
    return Z


