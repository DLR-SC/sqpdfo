# -*- coding: utf-8 -*-
"""
Created on Tue Jul 21 13:23:32 2015

@author: lien_ol
"""

from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
    
from bcdfo_build_QR_of_Y import bcdfo_build_QR_of_Y_

def bcdfo_augment_Y_(Ynew=None,Y_=None,whichmodel=None,shift_Y=None,Delta=None,normgx=None,kappa_ill=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 7-[Ynew,Y,whichmodel,shift_Y,Delta,normgx,kappa_ill].count(None)+len(args)

    Y=copy_(Y_)

    n,p1=size_(Y,nargout=2)
    if ((p1 >= ((n + 1) * (n + 2)) / 2) and (whichmodel != 3)):
        disp_(char(' === augment_Y: warning!!! The interpolation is already fully quadratic!'))
        disp_(char('     Ignoring augmentation...'))
        QZ,RZ,xbase,scale=bcdfo_build_QR_of_Y_(Y,whichmodel,shift_Y,Delta,normgx,kappa_ill,nargout=4)
    else:
        Y=concatenate_([Y, Ynew], axis=1)
#        Y=matlabarray([Y,Ynew])
        p1=p1 + size_(Ynew,2)
        QZ,RZ,xbase,scale=bcdfo_build_QR_of_Y_(Y,whichmodel,shift_Y,Delta,normgx,kappa_ill,nargout=4)
    return p1,QZ,RZ,Y,xbase,scale