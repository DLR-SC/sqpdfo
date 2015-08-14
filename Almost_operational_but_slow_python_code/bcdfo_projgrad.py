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
    
def bcdfo_projgrad_(n=None,x=None,g=None,bl=None,bu=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 5-[n,x,g,bl,bu].count(None)+len(args)
    gn=matlabarray()
    gnorm=0.0
    for i in arange_(1,n).reshape(-1):
        gi=g[i]
        if gi < 0.0:
            gn[i]=- min_(abs_(bu[i] - x[i]),- gi)
        else:
            gn[i]=min_(abs_(bl[i] - x[i]),gi)
        gnorm=max_(gnorm,abs_(gn[i]))
    return gnorm,gn