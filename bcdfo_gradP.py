# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 13:40:58 2015

@author: lien_ol
"""
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
def bcdfo_gradP_(P=None,x=None,xbase=None,scale=None,shift_Y=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 5-[P,x,xbase,scale,shift_Y].count(None)+len(args)
#
    n=length_(x)
    p1=length_(P)
    if (shift_Y):
        P=P.dot(scale.T)
        x=x - xbase
    ng=min_(n,p1 - 1)
    g=zeros_(n,1)
    g[1:ng]=P[2:ng + 1]
    nquad=p1 - n - 1
    if (nquad > 0):
        ndiag=min_(nquad,n)
        g[1:ndiag]=g[1:ndiag] + P[n + 2:n + 1 + ndiag].T.dot(x[1:ndiag])
        nquad=nquad - ndiag
        if (nquad > 0):
            k=2 * n + 1
            for i in arange_(1,n - 1).reshape(-1):
                nsd=min_(n - i,nquad)
                if (nsd > 0):
                    g[i + 1:i + nsd]=g[i + 1:i + nsd] + P[k + 1:k + nsd].T.dot(x[1:nsd])
                    g[1:nsd]=g[1:nsd] + P[k + 1:k + nsd].T.dot(x[i + 1:i + nsd])
                    k=k + nsd
                    nquad=nquad - nsd
                if (nquad == 0):
                    break
    return g