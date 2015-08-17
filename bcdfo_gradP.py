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
from copy import copy    
    
def bcdfo_gradP_(P_=None,x_=None,xbase=None,scale=None,shift_Y=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 5-[P,x,xbase,scale,shift_Y].count(None)+len(args)
#
    P=copy(P_)
    x=copy(x_)

    n=length_(x)
    p1=length_(P)
    if (shift_Y):
        P=P * scale.T
        x=x - xbase
    P=P.reshape(-1)
    x=x.reshape(-1)
    ng=min_(n,p1 - 1)
    g=zeros_(n,1).reshape(-1)
    g[0:ng]=P[1:ng + 1]
    nquad=p1 - n - 1
    if (nquad > 0):
        ndiag=min_(nquad,n)
        g[0:ndiag]=g[0:ndiag] + P[n + 1:n + 1 + ndiag] *(x[0:ndiag])
        nquad=nquad - ndiag
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