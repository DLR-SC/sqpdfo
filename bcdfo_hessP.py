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
from copy import copy    
from numpy import diag, zeros 
    
def bcdfo_hessP_(P_=None,x_=None,xbase=None,scale=None,shift_Y=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 5-[P,x,xbase,scale,shift_Y].count(None)+len(args)

    #Copy only if necessary
    if (shift_Y):
        x=copy(x_)
    else:
        x=x_

    P=copy(P_)
    n=length_(x)
    p1=length_(P)
    nquad=p1 - n - 1
    if (shift_Y):
        P=P*scale.T
        x=x - xbase
    P=P.reshape(-1)
    if (nquad > 0):
        ndiag=min_(nquad,n)
        H=diag(concatenate_([P[n + 1:n + 1 + ndiag],zeros(n - ndiag)],axis=1))
        nquad=nquad - ndiag
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