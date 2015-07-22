# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 15:32:32 2015

@author: lien_ol
"""
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_


def bcdfo_replace_in_Y_(QZ_=None,RZ_=None,ynew=None,Y_=None,j=None,xbase_=None,whichmodel=None,scale_=None,shift_Y=None,Delta=None,normgx=None,kappa_ill=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 12-[QZ,RZ,ynew,Y,j,xbase,whichmodel,scale,shift_Y,Delta,normgx,kappa_ill].count(None)+len(args)

    Y=copy_(Y_)
    QZ=copy_(QZ_)
    RZ=copy_(RZ_)
    xbase=copy_(xbase_)
    scale=copy_(scale_)

    Y[:,j]=ynew
    QZ,RZ,xbase,scale=bcdfo_build_QR_of_Y_(copy_(Y),whichmodel,shift_Y,Delta,normgx,kappa_ill,nargout=4)
    return QZ,RZ,Y,xbase,scale