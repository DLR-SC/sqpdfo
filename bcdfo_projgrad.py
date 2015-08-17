# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:25:42 2015

@author: lien_ol
"""

from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from numpy import array
    
def bcdfo_projgrad_(n=None,x=None,g=None,bl=None,bu=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 5-[n,x,g,bl,bu].count(None)+len(args)
    gn=zeros_(n,1)
    gnorm=0.0
    for i in range(0,n):
        gi=g[i,0]
        if gi < 0.0:
            gn[i,0]=- min_(abs(bu[i,0] - x[i,0]),- gi)
        else:
            gn[i,0]=min_(abs(bl[i,0] - x[i,0]),gi)
        gnorm=max_(gnorm,abs(gn[i,0]))
    return gnorm, gn