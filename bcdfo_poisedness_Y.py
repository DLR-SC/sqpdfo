# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 16:25:40 2015

@author: lien_ol
"""

from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
    
from bcdfo_find_new_yj import bcdfo_find_new_yj_
    
def bcdfo_poisedness_Y_(QZ=None,RZ=None,Y=None,eps_L=None,xbase=None,lSolver=None,whichmodel=None,hardcons=None,xl=None,xu=None,indfree=None,stratLam=None,scale=None,shift_Y=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 14-[QZ,RZ,Y,eps_L,xbase,lSolver,whichmodel,hardcons,xl,xu,indfree,stratLam,scale,shift_Y].count(None)+len(args)

    _lambda=0
    n,p1=size_(Y,nargout=2)
    Y_radius=0
    for j in range(1,p1): #it is indexed range(2,p1) in matlab
        Y_radius=max_(Y_radius,norm_(Y[:,j] - Y[:,0]))
    for j in range(1,p1): #it is indexed range(2,p1) in matlab
        if (hardcons == 1):
            y,improvement, msgTR=bcdfo_find_new_yj_bc_(QZ,RZ,Y,j,Y_radius,eps_L,xbase,lSolver,whichmodel,xl,xu,indfree,stratLam,scale,shift_Y,nargout=2)
        else:
            y,improvement, msgTR=bcdfo_find_new_yj_(QZ,RZ,Y,j,Y_radius,eps_L,xbase,lSolver,whichmodel,scale,shift_Y,nargout=2)
        _lambda=max_(improvement,_lambda)
    return _lambda,Y_radius