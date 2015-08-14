# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 13:40:57 2015

@author: lien_ol
"""
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
def bcdfo_hessP_(P_=None,x_=None,xbase=None,scale=None,shift_Y=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 5-[P,x,xbase,scale,shift_Y].count(None)+len(args)

    #Copy only if necessary
    if (shift_Y):
        P=copy_(P_)
        x=copy_(x_)
    else:
        P=P_
        x=x_

    n=length_(x)
    p1=length_(P)
    nquad=p1 - n - 1
    if (shift_Y):
        P=P.dot(scale.T)
        x=x - xbase
    if (nquad > 0):
        ndiag=min_(nquad,n)
        H=diag_(concatenate_([P[n + 2:n + 1 + ndiag],zeros_(1,n - ndiag)],axis=1))
        nquad=nquad - ndiag
        if (nquad > 0):
            k=2 * n + 1
            for i in arange_(1,n - 1).reshape(-1):
                nsd=min_(n - i,nquad)
                if (nsd > 0):
                    for j in arange_(1,nsd).reshape(-1):
                        H[i + j,j]=P[k + j]
                        H[j,i + j]=P[k + j]
                    k=k + nsd
                    nquad=nquad - nsd
                if (nquad == 0):
                    break
    else:
        H=zeros_(n,n)
    return H