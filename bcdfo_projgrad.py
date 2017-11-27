# -*- coding: utf-8 -*-

from runtime import *
    
def bcdfo_projgrad_(n=None,x=None,g=None,bl=None,bu=None,*args,**kwargs):
    """
#  This function computes the projected gradient and its
#  infinity norm.
#
#  INPUTS:
#
#  n     : dimension
#  x     : current iterate
#  g     : gradient at x
#  bl    : lower bounds
#  bu    : upper bounds
#
#  OUTPUTS:
#
#  gnorm : infinity norm of the projected gradient
#  gn    : projected gradient vector
#
    """
    
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
        # infinity-norm
        gnorm=max_(gnorm,abs(gn[i,0]))
    return gnorm, gn