# -*- coding: utf-8 -*-
"""
Created on Mon Jul 20 13:04:32 2015

@author: lien_ol
"""
from __future__ import division
try:
    from runtime import *
except ImportError:
    from smop.runtime import *
    
from bcdfo_evalP import bcdfo_evalP_
import numpy
    
def bcdfo_evalL_(QZ=None,RZ=None,Y=None,choice_set=None,x=None,xbase=None,whichmodel=None,scale=None,shift_Y=None,*args,**kwargs):
#    varargin = cellarray(args)
#    nargin = 9-[QZ,RZ,Y,choice_set,x,xbase,whichmodel,scale,shift_Y].count(None)+len(args)

    n,p1=size_(Y,nargout=2)
    I=eye_(p1)
    lc=length_(choice_set)
    q=((n + 1) * (n + 2)) / 2
    values=zeros_(p1,1)
    if (whichmodel == 3 and p1 < q):
        whichmodel=2
    if (whichmodel == 0):
        values[choice_set]=bcdfo_evalP_(I[choice_set,:] * (numpy.linalg.solve(RZ,QZ.T)),x,xbase,scale,shift_Y)
    else:
        if (whichmodel == 1):
            if (p1 == n + 1 or p1 == q):
                values[choice_set]=bcdfo_evalP_(I[choice_set,:] * (QZ / RZ.T),x,xbase,scale,shift_Y)
            else:
                if (shift_Y):
                    xbase=Y[:,1]
                    scaleY=0
                    for i in arange_(1,p1).reshape(-1):
                        Y[:,i]=Y[:,i] - xbase
                        scaleY=max_(scaleY,norm_(Y[:,i]))
                    Y=Y / scaleY
                M=bcdfo_evalZ_(Y,q).T
                MQ=M[:,n + 2:q]
                if (shift_Y):
                    phi=bcdfo_evalZ_((x - xbase) * scale[2],q)
                else:
                    phi=bcdfo_evalZ_(x,q)
                values[choice_set]=[I[choice_set,:],zeros_(lc,n + 1)] * (QZ * (numpy.linalg.solve(RZ.T,[[MQ * phi[n + 2:q]],[phi[1:n + 1]]])))
        else:
            if (whichmodel == 2):
                if (p1 < q):
                    values[choice_set]=bcdfo_evalP_(I[choice_set,:] * (pinv_(RZ) * QZ.T),x,xbase,scale,shift_Y)
                else:
                    values[choice_set]=bcdfo_evalP_(I[choice_set,:] * (numpy.linalg.solve(RZ,QZ.T)),x,xbase,scale,shift_Y)
            else:
                if (whichmodel == 3):
                    values[choice_set]=bcdfo_evalP_(I[choice_set,:] * (QZ * pinv_(RZ.T)),x,xbase,scale,shift_Y)
    return values